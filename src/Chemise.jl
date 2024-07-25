module Chemise

export ReactionSet, Biblio, @react_str, species, idx, getdens, setdens!, derivs, derivs!,
    nspecies, loadtable, RateLookup, .., @withref, @kexpr, @kexprs, writelatex, AbstractRate

using Printf
using StaticArrays
using DocStringExtensions
using MacroTools
using MacroTools: @capture, postwalk
using LaTeXStrings
using Latexify
using REPL

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
