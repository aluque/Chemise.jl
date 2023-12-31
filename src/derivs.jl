function _rate_expr(S, F, L, idx, pref)
    if isnothing(pref)
        expr = :(*(evalk(rs.reactions[$idx], args...)))
    else
        v = Symbol("prefetch_" * pref)
        expr = :(*(evalk(rs.reactions[$idx], args...; prefetch=$v)))
    end
    
    for (spec, coeff) in L
        i = speciesindex(S, spec)
        if !isnothing(i)            
            push!(expr.args, coeff == 1 ? :(n[$i]) : :(n[$i]^$coeff))
            continue
        end
        
        i = speciesindex(F, spec)
        if !isnothing(i)
            push!(expr.args, coeff == 1 ?
                :(evalfixed(rs.fixedval[$i], x)) :
                :(evalfixed(rs.fixedval[$i], x)^$coeff))
        end
    end
    return expr
end

"""
Obtain a list of expressions for the derivatives of each of the species in `S` given reactions in 
`R`.
"""
function _derivs_expr_list(S, R)    
    exprs = Any[]
    for (i, spec) in enumerate(S)
        expri = :(+())
        for (j, r) in enumerate(R.parameters)
            c = 0
            for (spec1, coeff) in lhs(r)
                spec1 == spec || continue
                c -= coeff
            end
            for (spec1, coeff) in rhs(r)
                spec1 == spec || continue
                c += coeff
            end
            if c != 0
                push!(expri.args, :($c * _rates[$j]))
            end            
        end
        push!(exprs, length(expri.args) > 1 ? expri : 0)
    end

    return exprs
end

"""
Build an expression that constructs an StaticArray with the derivatives of each species.
"""
function _derivs_expr_sa(S, R)
    expr_sa = :(SA[])
    exprs = _derivs_expr_list(S, R)
    for item in exprs
        push!(expr_sa.args, item)
    end
    
    return expr_sa
end

"""
Build an expression that sets the derivatives of each species into a vector variable dn
"""
function _derivs_expr_set(S, R, dn_var)
    expr_set = quote end
    exprs = _derivs_expr_list(S, R)
    for (i, item) in enumerate(exprs)
        push!(expr_set.args, :(($dn_var)[$i] = $item))
    end
    
    return expr_set
end


function _prefetching(R, G)
    expr = quote end
    prefetched = Set{String}()
    for (i, r) in enumerate(R.parameters)
        isactive(r, G) || continue
        
        v = varname(r)
        if !isnothing(varname(r)) && !(v in prefetched)
            push!(prefetched, v)
            v = Symbol("prefetch_" * v)
            push!(expr.args, :($v = prefetch(rs.reactions[$i], args...)))
        end        
    end
    return expr
end

"""
Compute derivatives for the reaction set `rs` with given densities `n` and external variables
`args...` (e.g. electric field). `groups` is a possible set of reaction groups that we want to
(exclusively) activate; it must be provided as `Val((:group1, :group2...))` or Val(group).
`x` is a location; it can be ignored or used to set inhomogeneous fixed values.
"""
@generated function derivs(n::AbstractVector, rs::ReactionSet{S, R, F}, groups::Val{G}, 
                           args...; x=nothing) where {S, R, F, G}
    rates = :(())
    for (idx, r) in enumerate(R.parameters)
        if isactive(r, G)
            pref = varname(r)
            push!(rates.args, _rate_expr(S, F, lhs(r), idx, pref))
        else
            push!(rates.args, 0)
        end
    end
    pref = _prefetching(R, G)
    dexpr = _derivs_expr_sa(S, R)
    expr = quote
        $pref
        _rates = $rates
        _dn = $dexpr
        return _dn
    end
    return expr
end

derivs(n::AbstractVector, rs::ReactionSet{S, R, F}, args...; x=nothing) where {S, R, F} = derivs(n, rs, Val{()}(), args...; x)


"""
Store in `dn` the derivatives for the reaction set `rs` with given densities `n` and external variables
`args...` (e.g. electric field). `groups` is a possible set of reaction groups that we want to
(exclusively) activate; it must be provided as `Val((:group1, :group2...))` or Val(group).
`x` is a location; it can be ignored or used to set inhomogeneous fixed values.
"""
@generated function derivs!(dn, n::AbstractVector, rs::ReactionSet{S, R, F}, groups::Val{G},
                            args...; x=nothing) where {S, R, F, G}
    rates = :(())
    for (idx, r) in enumerate(R.parameters)
        if isactive(r, G)
            pref = varname(r)
            push!(rates.args, _rate_expr(S, F, lhs(r), idx, pref))
        else
            push!(rates.args, 0)
        end
    end
    pref = _prefetching(R, G)
    dexpr = _derivs_expr_set(S, R, :dn)

    expr = quote
        $pref
        _rates = $rates
        _dn = $dexpr
        return _dn
    end
    return expr
end

derivs!(dn, n::AbstractVector, rs::ReactionSet{S, R, F}, args...; x=nothing) where {S, R, F} = derivs!(dn, n, rs, Val{()}(), args...; x)
