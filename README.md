# ArbRadixFloatingPoints.jl
Floating point numbers with arbitrary radixes (may be negative or nonreal)

```jl
julia> using ArbRadixFloatingPoints

julia> radix = 2 - im; precision = 16; num = 5 + 3im;

julia> y = convert(ArbRadixFloat{radix, precision}, num) #Convert a number into an ArbRadixFloat
(1̅.232222222222150)₂ ₋ ₁ᵢₘ × (2 - 1im)²

julia> float(y) #Conver an ArbRadixFloat back to a standard floating point number
4.999999999999999 + 3.0im

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
