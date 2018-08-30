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


function convert(::Type{ArbRadixFloat{radix,precision,Tdigit}},
                      num :: Number) where radix where precision where Tdigit
    if abs(radix) >= 1
        exponent = find_exponent(num, radix)
        significand = factorize_leastsquares(num, radix, precision, Tdigit, exponent)
        ArbRadixFloat{radix,precision,Tdigit}(significand, exponent)
    else #Compute using inverse radix and then convert back
        iradix = inv(radix)
        exponent = find_exponent(num, iradix)
        significand = factorize_leastsquares(num, iradix, precision, Tdigit, exponent)
        ArbRadixFloat{radix,precision,Tdigit}(reverse!(significand), precision-1-exponent)
    end
end

function find_exponent(num, radix)
    log_abs_radix = log(abs(radix))
    if log_abs_radix == 0
        #Do not attempt to normalize significand if radix has magnitude 1
        exponent = 0
    else #log_abs_radix != 0
        exponent = floor(Int, log(abs(num))/log_abs_radix)
    end
    exponent
end

function factorize_leastsquares(num, radix, precision, Tdigit, exponent;
                   verbose::Bool=false)

    significand = zeros(Tdigit, precision)

    # Super expensive greedy exhaustion following a full least squares solve
    # This is completely overkill because
    # D is low rank and symmetric (rank 1 for real radix, 2 for complex)
    # we only need the first entry of signif
    # and also "wrong" because the least squares solution is "undesirable"
    # because if you look at unit radixes like im, you will see the solution
    # norm spread out over many elements but really you just want the sparsest one
    # for unit radixes there is usually no unique numeration
    # and this method will eventually converge on one particular numeration
    # with the mass focused toward the last few 'digits'
    # but in some sense it's not the intuitive one; want the mass focused toward
    # the front
    if verbose
        println("Starting norm: ", abs(num))
        println("Starting exponent: ", exponent)
        println("Power\t| Digit\t| Solution \t| Residual norm\t| Change")
        println("-----\t| -----\t| -------- \t| -------------\t| ------")
    end

    #Element type of working intermediate problem
    eltype = typeof(float(real(radix)))

    for digit_id=1:precision
        nk = precision+1-digit_id
        D = zeros(eltype, nk, nk)
        c = zeros(eltype, nk)
        r = float(radix)
        rconj = conj(r)
        for i=1:nk, j=1:nk
            z = rconj^-(i+digit_id-1)*r^-(j+digit_id-1)
            D[i,j] = real(z+conj(z))
        end
        for i=1:nk
            z = rconj^-(i+digit_id-1)*num*r^-(exponent+1)
            c[i] = real(z+conj(z))
        end
        signif = pinv(D)*c
        significand[digit_id] = digit = round(Tdigit, signif[1])
        powr = exponent+1-digit_id
        Δ = digit * r^powr
        num -= Δ
        if verbose
            println(powr, "\t|", digit, "\t|", signif, "\t|", abs(num), "\t|", abs(Δ))
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

