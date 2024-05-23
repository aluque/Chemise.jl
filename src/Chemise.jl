module Chemise

export ReactionSet, Biblio, @react_str, species, idx, getdens, setdens!, derivs, derivs!,
    nspecies, loadtable, RateLookup, .., @withref, @kexpr

using Printf
using StaticArrays
using DocStringExtensions
using MacroTools
using LaTeXStrings
using Latexify

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
include("latexify.jl")

end
