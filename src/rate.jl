abstract type AbstractRate; end

evalk(r::Reaction, args...; prefetch=nothing) = _evalk(prefetch, r.k, args...)
_evalk(::Nothing, r, args...) = evalk(r, args...)
_evalk(prefetch, r, args...) = evalk(r, args...; prefetch)
       
evalk(f::Real, args...) = f
evalk(f::Function, args...) = f(args...)

varname(::Type{Reaction{S, K, G}}) where {S, K, G} = varname(K)

varname(::Type{T}) where T <: Real = nothing
varname(::Type{T}) where T <: Function = nothing
varname(::Type{T}) where T = nothing


"""
Wrap a rate and a bibliography reference.
"""
struct Biblio{K, S <: AbstractString} <: AbstractRate
    k::K
    ref::Vector{S}
end

Biblio(k, s::AbstractString) = Biblio(k, [s])

evalk(f::Biblio, args...) = evalk(f.k, args...)


"""
An expression index. It can store anything that can be transformed into latex.
"""
const EXPR_INDEX = Any[]

"""
Wrap a rate with a julia expression.
"""
struct RateExpr{K} <: AbstractRate
    k::K

    # A julia Expr is mutable and not concrete so if we store it within the struct we turn it
    # into non-concrete, with a performance penaly. What we do is keeping a global list of expressions
    # and store here only an index into it.
    idx::Int

    function RateExpr(k, expr)
        idx = lastindex(push!(EXPR_INDEX, expr))
        return new{typeof(k)}(k, idx)
    end    
end


"""
Build a `RateExpr` rate from a julia expression. The expression is evaluated immediately to set
the rate coefficient but it is also saved as a julia Expr to be eventually converted into LaTeX.
"""
macro kexpr(x)
    :(RateExpr($(esc(x)), $(Expr(:quote, x))))
end

evalk(k::RateExpr) = evalk(k.k)


"""
Traverse a full `ReactionSet` definition and replace rate coefficient expressions by the
corresponding `RateExpr`. It is equivalent to passing rate coefficient arithmetic
expressions through a `@kexpr` macro.
"""
macro kexprs(s)
    out = postwalk(s) do x
        if @capture(x, s_String_string => k_ .. ref_)
            if k isa Number
                return x
            elseif Meta.isexpr(k, :call) && k.args[1] == :RateLookup
                return x
            else
                return :($s => Chemise.RateExpr($k, $(Expr(:quote, k))) .. $ref)
            end
        elseif @capture(x, s_String_string => k_)
            if k isa Number
                return x
            elseif Meta.isexpr(k, :call) && k.args[1] == :RateLookup
                return x
            else
                return :($s => Chemise.RateExpr($k, $(Expr(:quote, k))))
            end
        else
            return x
        end
    end
    return esc(out)
end

"""
Convenience notation for rates with a reference. `1.0 .. "ref"` returns `Biblio(1.0, "ref")`.
This avoids cumbersome repeats of the `Biblio(...)` constructor when defining a ReactionSet.
"""
..(k, ref::AbstractString) = Biblio(k, ref)
..(k, ref::Vector{<:AbstractString}) = Biblio(k, ref)


"""
Convenience macro to set reference(s) for a number of reaction rates.  Use it inside a list
of reactions such as

```
ReactionSet([ "e + O2 + O2 -> O2- + O2" => 0.1,
              @withref("Einstein1905",
                     "e + O3 -> O- + O2" => 7.4e-18,
                     "e + O3 -> O2- + O" => 1.24e-18),
              "O- + O -> e + O2" => 5e-16])
```
"""
macro withref(ref, exprs...)
    exprs1 = :([])
    for expr in exprs
        @capture(expr, a_ => b_) || error("all arguments to @withref must be pairs a => b")
        push!(exprs1.args, :($(esc(a)) => Biblio($(esc(b)), $(esc(ref)))))
    end
    res = :(($exprs1)...)
    
    return :(($exprs1)...)
end


"""
A rate taken from a lookup table.

The parameter `H` is a hash of the lookup table; this allows to for compile-time prefetch of the lookup.
"""
struct RateLookup{L, I, H} <: AbstractRate
    lookup::L
    index::I

    function RateLookup(lookup::L, index=nothing) where L
        H = hash(L)
        new{L, typeof(index), H}(lookup, index)
    end
end

function RateLookup(lookup, fname::String)
    index = addcol(lookup, fname)
    return RateLookup(lookup, index)
end

function RateLookup(lookup, index::Symbol)
    i = findfirst(==(index), lookup.colnames)
    isnothing(i) && throw(ArgumentError("Column $(index) not found in the lookup table"))

    RateLookup(lookup, i)
end


@inline evalk(f::RateLookup, args...; prefetch=nothing) = @inline f.lookup(args..., f.index; prefetch)
@inline prefetch(f::RateLookup, args...) = prefetch(f.lookup, args...)
@inline prefetch(r::Reaction, args...) = prefetch(r.k, args...)

varname(::Type{RateLookup{L, H}}) where {L, H} = string(H)


"""
A wrapping struct to allow cunstructions such as f * Rate(...) or Rate(...) / f.
Note that for this to work Rate must be <: AbstractRate (or, trivially, a number).
It will not work for functions because I do not want to type-pirate *(::Real, ::Function).
In that case use ScaledRate(f, ::Function).
"""
struct ScaledRate{T, R}
    factor::T
    rate::R
end


@inline evalk(f::ScaledRate, args...; kw...) = @inline f.factor * evalk(f.rate, args...; kw...)
@inline prefetch(f::ScaledRate, args...) = prefetch(f.rate, args...)
varname(::Type{ScaledRate{T, R}}) where {T, R} = varname(R)

Base.:*(factor::Real, k::AbstractRate) = ScaledRate(factor, k)
Base.:/(k::AbstractRate, denom::Real) = ScaledRate(1/denom, k)

