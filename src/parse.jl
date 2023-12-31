function Base.parse(::Type{Signature}, signature)
    arrow = r"\s*[-=]+>\s*"
    plus = r"\s+\+\s+"

    (lhs, rhs) = (filter(s -> length(s) > 0, split(x, plus)) for x in split(signature, arrow))
    lhs, rhs = [[_multiplicity(item) for item in x] for x in (lhs, rhs)]

    return Signature(Tuple(lhs), Tuple(rhs))
end

function _multiplicity(code)
    rx = r"((\d+)\s*\*\s*)?([\w\(\).^+*-]+)"
    m = match(rx, code)
    (n1, s) = (m.captures[2], m.captures[3])
    
    n = isnothing(n1) ? 1 : parse(Int, n1)

    return (Symbol(s), n)
end

macro react_str(str)
    parse(Signature, str)
end
