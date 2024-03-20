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
struct Biblio{K, S <: AbstractString}
    k::K
    ref::S
end

evalk(f::Biblio, args...) = evalk(f.k, args...)


"""
Convenience notation for rates with a reference. `1.0 .. "ref"` returns `Biblio(1.0, "ref")`.
This avoids cumbersome repeats of the `Biblio(...)` constructor when defining a ReactionSet.
"""
..(k, ref::AbstractString) = Biblio(k, ref)


macro withref(ref, exprs...)
    exprs1 = :([])
    for expr in exprs
        @capture(expr, a_ => b_) || error("all arguments to @withref must be pairs a => b")
        push!(exprs1.args, :($a => Biblio($b, $ref)))
    end
    res = :(($exprs1)...)
    @show res
    
    return :(($exprs1)...)
end


"""
A rate taken from a lookup table.

The parameter `H` is a hash of the lookup table; this allows to for compile-time prefetch of the lookup.
"""
struct RateLookup{L, H}
    lookup::L
    index::Int

    function RateLookup(lookup::L, index::Int) where L
        H = hash(L)
        new{L, H}(lookup, index)
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
