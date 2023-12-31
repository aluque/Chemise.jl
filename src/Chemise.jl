module Chemise

export ReactionSet, Biblio, @react_str, species, idx, getdens, setdens!, derivs, derivs!,
    nspecies

using StaticArrays

include("signature.jl")
include("reaction.jl")
include("rate.jl")
include("reactionset.jl")
include("derivs.jl")
include("parse.jl")

end
