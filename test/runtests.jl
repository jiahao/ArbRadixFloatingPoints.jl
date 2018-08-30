using ArbRadixFloatingPoints
using Test

let x = 2im, z = convert(ArbRadixFloat{x,2}, x)
    @test z.exponent == 1
    @test z.significand == [1, 0]
    @test float(z) == x
end

let x = 2im, z = convert(ArbRadixFloat{-0.5im,2}, x)
    @test z.exponent == 0
    @test z.significand == [0, 1]
    @test float(z) == x
end

