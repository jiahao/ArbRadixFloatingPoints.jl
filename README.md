# ArbRadixFloatingPoints.jl
Floating point numbers with arbitrary radixes (may be negative or nonreal)

```
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

