a = 1
λ = 1
β = 1
μ1 = 1
ϵ_cut = 3
ħ = 1
Ω = 10
ω = 1

g1 = 2*ħ * Ω*a


function Lerch(z, j)
    Φ = 0
    for i in 1:1000
        Φ += z^(i-1)/(j+z)
    end
    return Φ
end

function gamma(μ)
    gam = 0
    for j in 1:1000
        Φ = Lerch.(exp.( β * μ - 2*β*ϵ_cut), j)
        gam += exp.(β * μ*(j+1) - 2 * β * ϵ_cut *j ) * Φ ^ 2
    end
    γ = 8*a^2 / λ^2 * gam
    return γ
end

function Me(μ)
    M = 16* pi *a^2/(exp(β*( μ - ϵ_cut ) ) -1)
    return M
end