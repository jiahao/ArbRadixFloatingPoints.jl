import Base.show

function show(io::IO, z::ArbRadixFloat{radix,precision,Tdigit}) where {radix,precision,Tdigit}
    print(io, '(')
    print(io, prettystring_significand(z.significand))
    print(io, ')')
    print(io, subscriptstring(radix))
    print(io, " × (")
    print(io, radix)
    print(io, ')')
    print(io, superscriptstring(z.exponent))
end

"This character gets used to print negative digits"
const overline = '\u305'

function prettystring_significand(significand)
    s = IOBuffer()

    print_point = true
    #Test range
    maxdigit = maximum(abs.(extrema(significand)))
    if maxdigit >= 36 #Do not attempt to pretty print
        print(s, significand)
    else# maxdigit < 36, use base36
        for digit in significand
            absdigit = abs(digit)
            if absdigit > 10
                print(s, 'a'+absdigit-10)
            else
                print(s, absdigit)
            end
            if digit < 0
                print(s, overline)
            end
            if print_point
                print(s, '.')
                print_point = false
            end
        end
    end
    String(take!(s))
end

const superscriptize = Dict(
    '1' => '¹',
    '2' => '²',
    '3' => '³',
    '4' => '⁴',
    '5' => '⁵',
    '6' => '⁶',
    '7' => '⁷',
    '8' => '⁸',
    '9' => '⁹',
    '0' => '⁰',

    '-' => '\u207b',
    '!' => 'ᵎ',

    'A' => 'ᴬ',
    'B' => 'ᴮ',
    'D' => 'ᴰ',
    'E' => 'ᴱ',
    'G' => 'ᴳ',
    'H' => 'ᴴ',
    'I' => 'ᴵ',
    'J' => 'ᴶ',
    'K' => 'ᴷ',
    'L' => 'ᴸ',
    'M' => 'ᴹ',
    'N' => 'ᴺ',
    'O' => 'ᴼ',
    'P' => 'ᴾ',
    'R' => 'ᴿ',
    'T' => 'ᵀ',
    'U' => 'ᵁ',
    'V' => 'ⱽ',
    'W' => 'ᵂ',

    'a' => 'ᵃ',
    'b' => 'ᵇ',
    'c' => 'ᶜ',
    'd' => 'ᵈ',
    'e' => 'ᵉ',
    'f' => 'ᶠ',
    'g' => 'ᵍ',
    'h' => 'ʰ',
    'i' => 'ⁱ',
    'j' => 'ʲ',
    'k' => 'ᵏ',
    'l' => 'ˡ',
    'm' => 'ᵐ',
    'n' => 'ⁿ',
    'o' => 'ᵒ',
    'p' => 'ᵖ',
    'r' => 'ʳ',
    's' => 'ˢ',
    't' => 'ᵗ',
    'u' => 'ᵘ',
    'v' => 'ᵛ',
    'w' => 'ʷ',
    'x' => 'ˣ',
    'y' => 'ʸ',
    'z' => 'ᶻ',
)
function superscriptstring(num)
    s = IOBuffer()
    for d in string(num)
        print(s, get(superscriptize, d, d))
    end
    String(take!(s))
end

#Most letters don't have a subscript counterpart...
const subscriptize = Dict(
    '1' => '₁',
    '2' => '₂',
    '3' => '₃',
    '4' => '₄',
    '5' => '₅',
    '6' => '₆',
    '7' => '₇',
    '8' => '₈',
    '9' => '₉',
    '0' => '₀',

    '+' => '₊',
    '-' => '₋',
    '=' => '₌',
    '(' => '₍',
    ')' => '₎',

    'a' => 'ₐ',
    'e' => 'ₑ',
    'h' => 'ₕ',
    'i' => 'ᵢ',
    'j' => 'ⱼ',
    'k' => 'ₖ',
    'l' => 'ₗ',
    'm' => 'ₘ',
    'n' => 'ₙ',
    'o' => 'ₒ',
    'p' => 'ₚ',
    'r' => 'ᵣ',
    's' => 'ₛ',
    't' => 'ₜ',
    'u' => 'ᵤ',
    'v' => 'ᵥ',
    'x' => 'ₓ',

    'β' => 'ᵦ',
    'γ' => 'ᵧ',
    'ρ' => 'ᵨ',
    'ϕ' => 'ᵩ',
    'φ' => 'ᵩ',
    'χ' => 'ᵪ',
)
function subscriptstring(num)
    s = IOBuffer()
    for d in string(num)
        print(s, get(subscriptize, d, d))
    end
    String(take!(s))
end

