module Chemise

export ReactionSet, Biblio, @react_str, species, idx, getdens, setdens!, derivs, derivs!,
    nspecies, loadtable, RateLookup, .., @withref

using StaticArrays
using DocStringExtensions
using MacroTools

@template DEFAULT =
    """
    $(SIGNATURES)
    $(DOCSTRING)
    """

include("signature.jl")
include("reaction.jl")
include("rate.jl")
include("reactionset.jl")
include("derivs.jl")
include("parse.jl")
include("lookup.jl")

end
