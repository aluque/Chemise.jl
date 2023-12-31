"""
A struct to contain a reaction set.

The species are encoded in the type `S` whereas reactions are a tuple encoded in `R`. `F` contains
a possible set of 'fixed' densities with values in the `fixed` field.
"""
struct ReactionSet{S, R, F, FT}
    reactions::R
    fixedval::FT
    
    function ReactionSet(species::Tuple{Vararg{Symbol}},
                         reactions::Tuple{Vararg{<:Reaction}},
                         fixed::Tuple{Vararg{Symbol}},
                         fixedval)
        new{species, typeof(reactions), fixed, typeof(fixedval)}(reactions, fixedval)
    end
end

function ReactionSet(reactions::Pair...; fix=())
    fixed = first.(fix)
    fixedval = last.(fix)

    _reactions = vcat(map(((g, r),) -> map(((sig, k),) -> Reaction(sig, k, g), r), reactions)...)
    species = tuple(list_of_species(_reactions, fixed)...)

    return ReactionSet(species, tuple(_reactions...), Tuple(fixed), Tuple(fixedval))
end

ReactionSet(reactions::Vector; kw...) = ReactionSet(nothing => reactions; kw...)


nspecies(rs::ReactionSet{S}) where {S} = length(S)

@generated function species_charge(::ReactionSet{S}) where S
    guess_charge(y) = y == :e ? -1 : count(==('+'), string(y)) - count(==('-'), string(y))
    expr = :(SA[])
    for spec in S
        push!(expr.args, guess_charge(spec))
    end
    return expr
end

function list_of_species(reactions, fixed)
    a = Set{Symbol}()
    for r in reactions
        for (s, _) in lhs(r)
            if isnothing(speciesindex(fixed, s))
                push!(a, s)
            end
        end

        for (s, _) in rhs(r)
            if isnothing(speciesindex(fixed, s))
                push!(a, s)
            end
        end
    end

    # Electrons are always first
    if :e in a
        delete!(a, :e)
        return [:e; sort!(collect(a))]
    else
        return sort!(collect(a))
    end
end


species(rs::ReactionSet{S}) where S = S
idx(rs::ReactionSet{S}, s::Symbol) where S = speciesindex(S, s)
idx(rs::ReactionSet{S}, s::String) where S = speciesindex(S, Symbol(s))

fixed_species(rs::ReactionSet{S, R, F}) where {S, R, F} = F
speciesindex(S, s::Symbol) = findfirst(==(s), S)

evalfixed(v::Real, x) = v
evalfixed(f::Function, x) = f(x)

function init(rs::ReactionSet{S}, dict, default=zero(valtype(dict))) where S
    SVector(map(k -> get(dict, Symbol(k), default), S))
end

function init!(u, rs::ReactionSet{S}, dict, default=zero(valtype(dict))) where S
    for (i, spec) in enumerate(S)
        u[i] = get(dict, Symbol(spec), default)
    end
end

getdens(u, rs::ReactionSet, spec) = u[idx(rs, spec)]
setdens!(u, rs::ReactionSet, spec, value) = u[idx(rs, spec)] = value
