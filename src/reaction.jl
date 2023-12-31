"""
A reaction with a rate.

The type G may contain a symbol that identifies different 'groups' of reactions. Sometimes you may want
to activate only certain groups, for example when computing derivatives pre/post photo-ioniation.
"""
struct Reaction{S <: Signature, K, G}
    signature::S
    k::K
    function Reaction(signature::S, k::K, G::Union{Symbol, Bool, Nothing}) where {S <: Signature, K}
        new{S, K, G}(signature, k)
    end
    function Reaction(signature::S, k::K) where {S <: Signature, K}
        new{S, K, nothing}(signature, k)
    end
end

Reaction(s::String, k, g=nothing) = Reaction(parse(Signature, s), k, g)

isactive(::Type{Reaction{S, H, G}}, groups::Union{Tuple, Int}) where {S, H, G} = isempty(groups) || G in groups
isactive(::Type{Reaction{S, H, G}}, groups::Symbol) where {S, H, G} = groups == G

signature(::Reaction{S}) where S = S
signature(::Type{Reaction{S, K, G}}) where {S, K, G} = S
lhs(::Reaction{S}) where {S}= lhs(S)
lhs(::Type{Reaction{S, K, G}}) where {S, K, G} = lhs(S)
rhs(::Reaction{S}) where {S} = rhs(S)
rhs(::Type{Reaction{S, K, G}}) where {S, K, G} = rhs(S)
