# Code for producing latex

@latexrecipe function f(x::RateLookup)
    return L"$k_\text{table}$"
end

@latexrecipe function f(x::RateExpr)
    return EXPR_INDEX[x.idx]
end

@latexrecipe function f(x::Biblio)
    # @info "In latexify Biblio" x
    # p = latexify(x.k)
    # @show p
    fmt --> FancyNumberFormatter()
    return x.k
end

@latexrecipe function f(r::Signature{L, R}) where {L, R}
    maybemult(m, s) = m > 1 ? "$m * $s" : "$s"
    
    lhs_str = join(map(((s, m),) -> maybemult(m, s), L), " + ")
    rhs_str = join(map(((s, m),) -> maybemult(m, s), R), " + ")

    return LaTeXString("\\cee{$lhs_str -> $rhs_str}")
end

"""
Write a LaTeX table to the `io` stream containing a summary of the `rs` reaction set.
"""
function writelatex(io::IO, rs::ReactionSet; env="longtable", preamble="", preheader="")
    rows = String[]
    print(io,
          """
          $(preamble)
          \\begin{$env}{p{2em}lllp{2.5cm}}
            $(preheader)
            \\rule{0pt}{2em} & \\textbf{Reaction} & & \\textbf{Rate (\\si{m^{3(n-1)}s^{-1}})} & \\textbf{Reference} \\\\ \\hline
        """)

    for (i, r) in enumerate(rs.reactions)
        s = replace(unicode2latex(latexify(r.signature, colsplit=false)), "->" => "& ->")
        k = unicode2latex(latexify(r.k, fmt=numformatter))
        c = refs(r.k)

        print(io, "\\textbf{$i} & $s & $k & $c \\\\\n")
    end

    print(io, """
        \\hline      
        \\end{$env}
        """)
    
end

writelatex(rs::ReactionSet) = writelatex(stdout, rs)

# This is to keep compatibility with Latexify and allow tables in jupyter 
@latexrecipe function f(rs::ReactionSet)
    lst = Any[]
    for (i, r) in enumerate(rs.reactions)
        s = unicode2latex(latexify(r.signature))
        k = unicode2latex(latexify(r.k, fmt=numformatter))
        c = LaTeXString(refs(r.k))

        push!(lst, [i, s, k, c]) 
    end
    arr = permutedims(hcat(lst...))

    env --> :mdtable
    adjustment --> [:c, :l, :l, :l]
    
    return arr
end


refs(x) = L"\text{Unknown}"
function refs(b::Biblio)
    inner = join(b.ref, ",")
    return "\\cite{$inner}"
end

const float_regex = r"(?'mantissa'(?'before_dp'(?'sign'-?)(?'before_dp_nosign'\d+))(\.(?'after_dp'\d+))?)(?'e_or_E'e)(?'raw_exp'(?'sign_exp'-?)\+?0*(?'mag_exp'\d+))"i

function numformatter(x)
    s = @sprintf("%.4g", x)
    m = match(float_regex, s)
    if isnothing(m)
        return s
    end

    if m[:mantissa] == "1"
        return "\\num{e$(m[:raw_exp])}"
    else
        return "\\num{$(m[:mantissa])e$(m[:raw_exp])}"
    end
end

"""
Replace unicode chars in `str` with their LaTeX equivalent.
"""
function unicode2latex(str)
    # As a table for latex replacements we use the julia REPL table.  Note however that this has
    # to be reversed into unicode => latex pairs.
    ls = collect(REPL.REPLCompletions.latex_symbols)
    return foldr((itm, s) -> replace(s, reverse(itm)), ls; init=str)
end
