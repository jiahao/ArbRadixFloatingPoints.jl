module ArbRadixFloatingPoints

using LinearAlgebra
using StaticArrays

export ArbRadixFloat

import Base: +, -, *, /, \, ^, abs, conj, convert, float, imag, precision, real
import LinearAlgebra: norm

"""
Floating point with arbitrary radix

# parameters
radix : "base" of floating point number
precision: number of digits
Tdigit : type of digit. If not specified in convert(), defaults to Int

significand and exponent are stored as fields
"""
struct ArbRadixFloat{radix,precision,Tdigit} <: Number
    significand :: SVector{precision,Tdigit}
    exponent :: Int
end

include("convert.jl")
include("show.jl")

# Basic floating point interface

precision(num::ArbRadixFloat{radix,precision,Tdigit}) where {radix,precision,Tdigit} = precision

conj(num::ArbRadixFloat{radix,precision,Tdigit}) where {radix,precision,Tdigit} =
    ArbRadixFloat{conj(radix),precision,Tdigit}(num.significand, num.exponent)

-(num::ArbRadixFloat{radix,precision,Tdigit}) where {radix,precision,Tdigit} =
    ArbRadixFloat{conj(radix),precision,Tdigit}(-num.significand, num.exponent)

function float(num::ArbRadixFloat{radix,precision,Tdigit}) where radix where precision where Tdigit
    anstype = promote_type(typeof(float(radix)), Tdigit)
    convert(anstype, ArbRadixFloat{float(radix),precision,Tdigit}(num.significand, num.exponent))
end

# These operations are just placeholders in lieu of actually doing the computations natively

for func in (:abs, :imag, :real, :norm)
    @eval begin
        function ($func)(num::ArbRadixFloat{radix,precision,Tdigit}) where {radix,precision,Tdigit}
            $(func)(convert(promote_type(typeof(radix), Tdigit), num))
        end
    end
end

for func in (:+, :-, :*, :/, :\, :^)
    @eval begin
        function ($op)(x::ArbRadixFloat{radix,precision,Tdigit}, y::ArbRadixFloat{radix,precision,Tdigit}) where {radix,precision,Tdigit}
            Twork = promote_type(typeof(radix), Tdigit)
            convert(ArbRadixFloat{radix,precision,Tdigit}, ($op)(convert(Twork, x), convert(Twork, y)))
        end
    end
end


end # module

