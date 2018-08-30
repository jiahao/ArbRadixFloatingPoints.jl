# ArbRadixFloatingPoints.jl
Floating point numbers with arbitrary radixes (may be negative or nonreal)

## Some usage examples

```jl
julia> using ArbRadixFloatingPoints

julia> radix = 2 - im; precision = 16; num = 5 + 3im;

julia> y = convert(ArbRadixFloat{radix, precision}, num) #Convert a number into an ArbRadixFloat
(1̅.232222222222150)₂ ₋ ₁ᵢₘ × (2 - 1im)²

julia> float(y) #Conver an ArbRadixFloat back to a standard floating point number
4.999999999999999 + 3.0im

julia> convert(ArbRadixFloat{6im, 5}, 1000) #Pretty printing uses base36 with overbars for negative digits
(0.s̅08̅0)₀ ₊ ₆ᵢₘ × (0 + 6im)³

julia> z = convert(ArbRadixFloat{y, precision}, 1000 - 39im) #A stranger example with an ArbRadixFloat radix
([-8, 49, 25, 11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])₍₁̅.₂₃₂₂₂₂₂₂₂₂₂₂₁₅₀₎₂ ₋ ₁ᵢₘ × ₍₂ ₋ ₁ᵢₘ₎² × ((1̅.232222222222150)₂ ₋ ₁ᵢₘ × (2 - 1im)²)³

julia> float(z) #Still works!
999.9999999999999 - 38.99999999999979im

```

```jl
julia> using ArbRadixFloatingPoints, Quaternions

julia> q1 = Quaternion(1, -2, 0.9, 0); q2 = Quaternion(10, -3.1, 50, -1);

julia> Base.float(q::Quaternion{T}) where T<:AbstractFloat = q #Define missing Quaternion method

julia> Q = convert(ArbRadixFloat{q1,16}, q2) #Printing is messy
(1̅.02̅1̅1̅2̅1̅22̅21̅2̅1111)Qᵤₐₜₑᵣₙᵢₒₙ{Fₗₒₐₜ₆₄}₍₁.₀, ₋₂.₀, ₀.₉, ₀.₀, fₐₗₛₑ₎ × (Quaternion{Float64}(1.0, -2.0, 0.9, 0.0, false))⁴

julia> abs(float(Q) - q2) #pretty lossy but doable
44.33522390835951
```

## Implementation note

Quite a lot of papers have been written on non-integer numerations, complex base systems, and the like.

Strangely, though, the non-integer floating point numbers do not appear to be written down anywhere, even though the generalization from ordinary integer-radixed floating point numbers is obvious to define.

Converting an arbitrary-radixed floating point number to an ordinary binary float is quite straightforward.
The converse, however, is quite tricky to get right!
The current algorithm is implemented in the (unexported) function `factorize_leastsquares`.
If we relax the problem of finding a digit expansion to a least squares problem, then the coefficients of the expansion can be found by solving a (massively rank-deficient) linear system of equations.
This package implements a hacky solution to obtain the most significant digit by rounding off the first component of the solution,
and then deflating the problem.
I don't think there is a guarantee of finding the minimum norm solution after rounding off.
In practice, this seems to work quite well.

My first attempt as this conversion is `factorize_greedy`, which is quite fragile and generates obviously worse solutions that the current method.

## Limitations

- By default, the digits are allowed to take any value that fits in machine `Int`s, _including negative values_. Restrictions to nonnegative digits, binary digits, or similar have not been implemented, so many numeration systems, such as Knuth's quater-imaginary system or the balanced ternary system, are _not_ computed by this package.
- No arithmetic on these numbers have been implemented. There's so much to do!
- The conversion from ordinary binary floating point to `ArbRadixFloat` is not guaranteed to be the best possible norm preserving conversion.
- I have no idea how to think about rounding modes for `ArbRadixFloat`s
