"""
A struct to contain the reactants and products, with their multiplicities.

To have performant schemes, we encode the reactants/products and their stochiometric coefficient
in the type, in tuples of 2-element tuples.
"""
struct Signature{L, R}
    function Signature(L::Tuple{Vararg{Tuple{Symbol, Int}}},
                       R::Tuple{Vararg{Tuple{Symbol, Int}}})
        return new{L, R}()
    end    
end

lhs(::Signature{L, R}) where {L, R} = L
lhs(::Type{Signature{L, R}}) where {L, R} = L
rhs(::Signature{L, R}) where {L, R} = R
rhs(::Type{Signature{L, R}}) where {L, R} = R

function Base.show(io::IO, r::Signature{L, R}) where {L, R}
    maybemult(m, s) = m > 1 ? "$m * $s" : "$s"
    
    lhs_str = join(map(((s, m),) -> maybemult(m, s), L), " + ")
    rhs_str = join(map(((s, m),) -> maybemult(m, s), R), " + ")
    print(io, "react\"", lhs_str, " -> ", rhs_str, "\"")
end

