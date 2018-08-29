module ArbRadixFloatingPoint

using LinearAlgebra
using StaticArrays

export ArbRadixFloat

import Base: -, *, conj, convert, abs, precision, real, float

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
function convert(::Type{T}, num::ArbRadixFloat{radix,precision,Tdigit}) where T <: Number where radix where precision where Tdigit
    y = zero(T)
    monomial = radix^num.exponent
    for digit in num.significand
        y += digit * monomial
        monomial /= radix
    end
    y
end

convert(::Type{ArbRadixFloat{radix,precision}}, num :: Number, verbose::Bool=false) where radix where precision =
    convert(ArbRadixFloat{radix,precision,Int}, num)

convert(::Type{ArbRadixFloat{radix,precision,Tdigit}}, num :: ArbRadixFloat{radix,precision,Tdigit}, verbose::Bool=false) where radix where precision where Tdigit = num


#TODO Bad things happen for bases of magnitude <= 1
function convert(::Type{ArbRadixFloat{radix,precision,Tdigit}},
                      num :: Number) where radix where precision where Tdigit

    exponent = find_exponent(num, radix, precision)
    significand = factorize_leastsquares(num, radix, precision, Tdigit, exponent)
    ArbRadixFloat{radix,precision,Tdigit}(significand,exponent)
end

function find_exponent(num, radix, precision)
    log_abs_radix = log(abs(radix))
    if log_abs_radix == 0
        #Do not attempt to normalize significand if radix has magnitude 1
        exponent = 0
    else #log_abs_radix != 0
        exponent = floor(Int, log(abs(num))/log_abs_radix)
        if log_abs_radix < 0 #Inverted radix with magnitude < 1
            exponent += precision+1
        end
    end
    exponent
end

function factorize_leastsquares(num, radix, precision, Tdigit, exponent;
                   verbose::Bool=false)

    significand = zeros(Tdigit, precision)

    #Super expensive greedy exhaustion following a full least squares solve
    #This is completely overkill because
    # D is low rank (rank 1 for real radix, 2 for complex)
    # we only need the first entry of signif
    for k=1:precision
        D = zeros(typeof(real(radix)), precision+1-k, precision+1-k)
        c = zeros(typeof(real(radix)), precision+1-k)
        r = float(radix)
        rconj = conj(r)
        for i=1:precision+1-k, j=1:precision+1-k
            z = rconj^-(i+k-1)*r^-(j+k-1)
            D[i,j] = real(z+conj(z))
        end
        for i=1:precision+1-k
            z = rconj^-(i+k-1)*num*r^-(exponent+1)
            c[i] = real(z+conj(z))
        end
        signif = pinv(D)*c
        significand[k] = d = round(Tdigit, signif[1])
        num -= d * r^(exponent+1-k)
    end
    significand
end

function factorize_greedy(num, radix, precision, Tdigit, exponent;
                   verbose::Bool=false)

    significand = zeros(Tdigit, precision)
    monomial = radix ^ exponent
    if verbose
        println("Starting norm: ", abs(num))
        println("Starting exponent: ", exponent)
        println("Power\t| Digit\t| Quotient\t| Residual norm\t| Change")
        println("-----\t| -----\t| --------\t| -------------\t| ------")
    end
    for digit_id in 1:precision
        quotient = real(conj(num)*monomial)/real(conj(monomial)*monomial)
        digit = round(Tdigit, quotient)
        num -= digit*monomial
        monomial /= radix
        significand[digit_id] = digit
        if verbose
            println(exponent+1-digit_id, "\t|", digit, "\t|", quotient, "\t|", abs(num), "\t|", abs(digit*monomial))
        end
        if num == 0 #Exactly representable
            if verbose
                println("Exact representation found")
            end
            break
        end
    end
    significand
end

precision(num::ArbRadixFloat{radix,precision,Tdigit}) where radix where precision where Tdigit = precision

conj(num::ArbRadixFloat{radix,precision,Tdigit}) where radix where precision where Tdigit =
    ArbRadixFloat{conj(radix),precision,Tdigit}(num.significand, num.exponent)

function -(num::ArbRadixFloat{radix,precision,Tdigit}) where radix where precision where Tdigit
    ArbRadixFloat{conj(radix),precision,Tdigit}(-num.significand, num.exponent)
end

# These operations are just placeholders in lieu of actually doing the computations natively

function abs(num::ArbRadixFloat{radix,precision,Tdigit}) where radix where precision where Tdigit
    abs(convert(promote_type(typeof(radix), Tdigit), num))
end

function real(num::ArbRadixFloat{radix,precision,Tdigit}) where radix where precision where Tdigit
    real(convert(promote_type(typeof(radix), Tdigit), num))
end

function float(num::ArbRadixFloat{radix,precision,Tdigit}) where radix where precision where Tdigit
    anstype = promote_type(typeof(float(radix)), Tdigit)
    convert(anstype, ArbRadixFloat{float(radix),precision,Tdigit}(num.significand, num.exponent))
end

function *(x::ArbRadixFloat{radix,precision,Tdigit}, y::ArbRadixFloat{radix,precision,Tdigit}) where radix where precision where Tdigit
    Twork = promote_type(typeof(radix), Tdigit)
    convert(ArbRadixFloat{radix,precision,Tdigit}, convert(Twork, x) * convert(Twork, y))
end

end # module

