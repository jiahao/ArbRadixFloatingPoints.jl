module ArbRadixFloatingPoint

using StaticArrays

export ArbRadixFloat

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

# Naive conversion and addition
function Base.convert(::Type{T}, num::ArbRadixFloat{radix,precision,Tdigit}) where T <: AbstractFloat where radix where precision where Tdigit
    y = zero(T)
    monomial = radix^num.exponent
    for digit in num.significand
        y += digit * monomial
        monomial /= radix
    end
    y
end

Base.convert(::Type{ArbRadixFloat{radix,precision}}, num :: Number) where radix where precision =
    Base.convert(ArbRadixFloat{radix,precision,Int}, num)

#TODO Bad things happen for bases of magnitude <= 1
function Base.convert(::Type{ArbRadixFloat{radix,precision,Tdigit}},
                      num :: Number) where radix where precision where Tdigit

    log_abs_radix = log(abs(radix))
    if log_abs_radix == 0
        #Do not attempt to normalize significand if radix has magnitude 1
        exponent = 0
    elseif  #log_abs_radix != 0
        exponent = floor(Int, log(abs(num))/log_abs_radix)
        if log_abs_radix < 0 #"Inverted radix with magnitude < 1
            exponent -= precision
        end
    end
    significand = zeros(Tdigit, precision)

    monomial = radix ^ exponent
    for digit_id in 1:precision
        quotient = real(conj(num)*monomial)/real(conj(monomial)*monomial)
        digit = round(Tdigit, quotient)
        num -= digit*monomial
        monomial /= radix
        significand[digit_id] = digit
        #println(digit_id, '\t', digit, '\t', quotient, '\t', abs(num))
        if num == 0 #Exactly representable
            break
        end
    end
    ArbRadixFloat{radix,precision,Tdigit}(significand,exponent)
end


Base.precision(num::ArbRadixFloat{radix,precision,Tdigit}) where radix where precision where Tdigit = precision

end # module

