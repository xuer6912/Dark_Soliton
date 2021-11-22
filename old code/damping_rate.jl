using Plots,LaTeXStrings
##
a = 5.8e-9 #m
T = 20e-9 #k
Tc = 60e-9 #k
m = 87 *1.6605402e-27#kg
kb = 1.380649e-23 #jk-1
μ1 = kb*T #j
ħ = 1.054571817e-34 #
Ω = 2*pi*100 #hz
ω = 2 * pi* 6 #hz
λ = sqrt(2*pi*ħ^2/(m*kb*T)) #m
ϵ_cut = 3*μ1
ϵ_cut2(μ) = 3*μ
g1 = 2*ħ * Ω*a
β = 1/(kb*T)
σ = sqrt(ħ/m/Ω)
##

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
        Φ = Lerch.(exp.( β * μ - 2*β*ϵ_cut2(μ)), j)
        gam += exp.(β * μ*(j+1) - 2 * β * ϵ_cut2(μ) *j ) * Φ ^ 2
    end
    γ = 8*a^2 / λ^2 * gam
    return γ
end
function Γ_γ(μ)
    Γ_γ =gamma(μ)/3 *μ/ħ
    return Γ_γ
end  

function Γ_M(μ)
    M1 = (16* pi *a^2/(exp(β*(ϵ_cut2(μ)-μ))-1)) / sqrt(8*pi*σ^2)
    Γ_M = M1 * (2/15) * μ.^2 /ħ/g1
    return Γ_M
end

##
μr = LinRange(0.0001, 2.8,100)*μ1
p1 = plot(μr/(ħ*ω), Γ_γ.(μr), xaxis=:log, yaxis=:log,label=L"\Gamma_\gamma")
plot!(μr/(ħ*ω) , Γ_M.(μr), xaxis=:log, yaxis=:log,label=L"\Gamma_M")
#
xlabel!(L"\mu/\hbar\omega")
ylabel!(L"\Gamma")
##
savefig(L"\Gamma")
##
