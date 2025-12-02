function Adapt.adapt_structure(to, rs::ReactionSet{S, R, F, FT}) where {S, R, F, FT}
    reactions = map(r -> Adapt.adapt(to, r), rs.reactions)
    fixedval = map(f -> Adapt.adapt(to, f), rs.fixedval)

    return ReactionSet(S, reactions, F, fixedval)
end

function Adapt.adapt_structure(to, r::Reaction{S, K, G}) where {S, K, G}
    k = Adapt.adapt(to, r.k)
    return Reaction(r.signature, k, G)
end

# Biblio rates just ignore bibliography when passing to CUDA.
function Adapt.adapt_structure(to, k::Biblio)
    return Adapt.adapt(to, k.k)
end

function Adapt.adapt_structure(to, k::RateLookup)
    r = RateLookup(Adapt.adapt(to, k.lookup), k.index; k.narg)
    return r
end

function Adapt.adapt_structure(to, l::LookupTable)
    fx = Adapt.adapt_structure(to, l.fx)
    if l.gy isa Vector{<:Real}
        gy = Adapt.adapt(to, l.gy)
    else
        gy = Tuple(map(itm -> Adapt.adapt(to, itm), l.gy))
    end
    f = l.f
    g = l.g
    ginv = l.ginv
    
    # When passing to a kernel the column names are useless and cause problems so we strip them.
    colnames = nothing
    extrapol = l.extrapol
    
    return LookupTable(fx, gy; f, g, ginv, colnames, extrapol)    
end
