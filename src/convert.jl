function convert(::Type{T}, num::ArbRadixFloat{radix,precision,Tdigit}) where {T<:Number,radix,precision,Tdigit}
    y = zero(T)
    monomial = radix^num.exponent
    for digit in num.significand
        y += digit * monomial
        monomial /= radix
    end
    y
end

convert(::Type{ArbRadixFloat{radix,precision}}, num :: Number, verbose::Bool=false) where {radix,precision} =
    convert(ArbRadixFloat{radix,precision,Int}, num)

convert(::Type{ArbRadixFloat{radix,precision,Tdigit}}, num :: ArbRadixFloat{radix,precision,Tdigit},
        verbose::Bool=false) where {radix,precision,Tdigit} = num


function convert(::Type{ArbRadixFloat{radix,precision,Tdigit}},
                      num :: Number) where {radix,precision,Tdigit}
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

##############################################################################
# Unexported private methods for conversion
##############################################################################

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
        ir = inv(r)
        irconj = conj(ir)
        for i=1:nk, j=1:nk
            z = irconj^(i+digit_id-1)*ir^(j+digit_id-1)
            D[i,j] = real(z+conj(z))
        end
        for i=1:nk
            z = irconj^(i+digit_id-1)*num*ir^(exponent+1)
            c[i] = real(z+conj(z))
        end
        signif = pinv(D)*c
        significand[digit_id] = digit = round(Tdigit, signif[1])
        powr = exponent+1-digit_id
        if powr >= 0
            Δ = digit * r^powr
        else #Some number types don't have x^-n defined, manually rewrite as inv(x)^n
            Δ = digit * ir^-powr
        end
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

