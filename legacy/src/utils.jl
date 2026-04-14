"""
    CantorPair(x, y) -> Int

Generate a unique integer from two non-negative integers using the Cantor pairing function.

# Arguments
- `x::Int`: The first non-negative integer.
- `y::Int`: The second non-negative integer.

# Returns
- `Int`: A unique integer generated from the two input integers using the Cantor pairing function. Returns 0 if `y` is 0.
"""
function CantorPair(x, y)
    y == 0 && return 0
    return trunc(Int, (x^2 + 3 * x + 2 * x * y + y + y^2) / 2)
end