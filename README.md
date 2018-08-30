# ArbRadixFloatingPoints.jl
Floating point numbers with arbitrary radixes (may be negative or nonreal)

```jl
julia> using ArbRadixFloatingPoints

julia> radix = 2 - im; precision = 16; num = 5 + 3im;

julia> y = convert(ArbRadixFloat{radix, precision}, num) #Convert a number into an ArbRadixFloat
ArbRadixFloat{2 - 1im,16,Int64}([-1, 2, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 5, 0], 2)

julia> float(y) #Conver an ArbRadixFloat back to a standard floating point number
4.999999999999999 + 3.0im

julia> z = convert(ArbRadixFloat{y, precision}, 1000 - 39im) #A stranger example with an ArbRadixFloat radix
ArbRadixFloat{ArbRadixFloat{2 - 1im,16,Int64}([-1, 2, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 5, 0], 2),16,Int64}([-8, 49, 25, 11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 3)

julia> float(z) #Still works!
999.9999999999999 - 38.99999999999979im
```

```jl
julia> using ArbRadixFloatingPoints, Quaternions

julia> q1 = Quaternion(1, -2, 0.9, 0); q2 = Quaternion(10, -3.1, 50, -1);

julia> Base.float(q::Quaternion{T}) where T<:AbstractFloat = q #Define missing Quaternion method

julia> Q = convert(ArbRadixFloat{q1,16}, q2)
ArbRadixFloat{Quaternion{Float64}(1.0, -2.0, 0.9, 0.0, false),16,Int64}([-1, 0, -2, -1, -1, -2, -1, 2, -2, 2, -1, -2, 1, 1, 1, 1], 4)

julia> abs(float(Q) - q2) #pretty lossy but doable
44.33522390835951
```
