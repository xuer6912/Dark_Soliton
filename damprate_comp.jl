using Plots, LaTeXStrings
#physical constant
amu = 1.66053e-27
a0 = 5.29e-11 # Bohr radius
ħ = 1.055e-34
kB = 1.381e-23
##
#Rubidium 87
m = 87 * amu
a = 95 * a0
g = 4 * π * ħ^2 * a / m
##
#trap
wz = 2 * pi * 3
wr = 2 * pi * 300
az = sqrt(ħ / m / wz)
ar = sqrt(ħ / m / wr)
g1 = 2 * ħ * wr * a
##
#temperature
T = 50e-9
Tc = 60e-9
β(T) = 1 / (kB * T)
ϵc(μ) = 3 * μ
μ1 = ħ * wz * 10 # entering Thomas-Fermi regime
λdB(T) = sqrt(2π * ħ^2 / (m * kB * T))
λ = λdB(T)
function Lerch(z, j)
    Φ = 0
    for i in 1:1000
        Φ += z^(i - 1) / (j + z)
    end
    return Φ
end

begin
    function gamma(μ)
        gam = 0
        for j in 1:1000
            Φ = Lerch.(exp.(β(T) * μ - 2 * β(T) * ϵc(μ)), j)
            gam += exp.(β(T) * μ * (j + 1) - 2 * β(T) * ϵc(μ) * j) * Φ^2
        end
        γ = 8 * a^2 / λ^2 * gam
        return γ
    end
    function Γ_γ(μ)
        Γ_γ = gamma(μ) / 3 * μ / ħ
        return Γ_γ
    end

    function Γ_M(μ)
        M1 = (16 * pi * a^2 / (exp(β(T) * (ϵc(μ) - μ)) - 1)) / sqrt(8 * pi * ar^2)
        Γ_M = M1 * (2 / 15) * μ .^ 2 / ħ / g1
        return Γ_M
    end
end
########
T = 30e-9
ϵc(μ) = 2 * μ
μr = LinRange(1, 10, 100) * μ1
G1 = Γ_γ.(μr) * sqrt(2) / wz
m1 = Γ_M.(μr) * sqrt(2) / wz
ϵc(μ) = 2.5 * μ
#μr = LinRange(.5, 10,100)*μ1
G2 = Γ_γ.(μr) * sqrt(2) / wz
m2 = Γ_M.(μr) * sqrt(2) / wz
ϵc(μ) = 3 * μ
G3 = Γ_γ.(μr) * sqrt(2) / wz
m3 = Γ_M.(μr) * sqrt(2) / wz
plot(μr / (ħ * wz), G1, xaxis = :log, yaxis = :log, label = L"\gamma:2\mu", line = (:dash, :blue))
plot!(μr / (ħ * wz), G2, label = L"\gamma:2.5\mu", line = (:dash, :green))
plot!(μr / (ħ * wz), G3, label = L"\gamma:3\mu", line = (:dash, :red))
plot!(μr / (ħ * wz), m1, label = L"M_1:2\mu", line = (:blue))
plot!(μr / (ħ * wz), m2, label = L"M_1:2.5\mu", line = (:green))
plot!(μr / (ħ * wz), m3, label = L"M_1:3\mu", line = (:red), scale = :log10, minorticks = true, grid = false)
##
vline!([0.5wr / wz], label = L"50\hbar\omega_x")
xlabel!(L"\mu/\hbar\omega_x")
ylabel!(L"\Gamma/\omega_s")
savefig("dampratecut.pdf")
plot(μr / (ħ * wz), m1 ./ G1, xaxis = :log, yaxis = :log, label = L"2\mu", size = (600, 300), line = (:blue), legend = :bottomright)
plot!(μr / (ħ * wz), m2 ./ G2, xaxis = :log, yaxis = :log, label = L"2.5\mu", size = (600, 300), line = (:green))
plot!(μr / (ħ * wz), m3 ./ G3, xaxis = :log, yaxis = :log, label = L"3\mu", size = (600, 300), line = (:red), scale = :log10, minorticks = true, grid = false)
xlabel!(L"\mu/\hbar\omega_x")
ylabel!(L"\Gamma_{M_1}/ \Gamma_\gamma")
savefig("dampratecut_comp.pdf")

#################

ϵc(μ) = 3 * μ

T = 12e-9
GT1 = Γ_γ.(μr) * sqrt(2) / wz
MT1 = Γ_M.(μr) * sqrt(2) / wz
μr = LinRange(1, 10, 100) * μ1



T = 30e-9
GT2 = Γ_γ.(μr) * sqrt(2) / wz
MT2 = Γ_M.(μr) * sqrt(2) / wz



T = 48e-9
GT3 = Γ_γ.(μr) * sqrt(2) / wz
MT3 = Γ_M.(μr) * sqrt(2) / wz

p1 = plot(μr / (ħ * wz), GT1, xaxis = :log, yaxis = :log, label = L"\gamma:12nK", legend = :topright, line = (:dash, :blue))
plot!(μr / (ħ * wz), GT2, xaxis = :log, yaxis = :log, label = L"\gamma:30nK", line = (:dash, :green))
plot!(μr / (ħ * wz), GT3, xaxis = :log, yaxis = :log, label = L"\gamma:48nK", line = (:dash, :red))
plot!(μr / (ħ * wz), MT1, xaxis = :log, yaxis = :log, label = L"M_1:12nK", size = (600, 300), line = (:blue))
plot!(μr / (ħ * wz), MT2, xaxis = :log, yaxis = :log, label = L"M_1:30nK", size = (600, 300), line = (:green))
plot!(μr / (ħ * wz), MT3, xaxis = :log, yaxis = :log, label = L"M_1:48nK", size = (600, 300), line = (:red), scale = :log10, minorticks = true, grid = false)

vline!([0.5wr / wz], label = L"50\hbar\omega_x")
xlabel!(L"\mu/\hbar\omega_x")
ylabel!(L"\Gamma/\omega_s")


savefig("damprateT.pdf")

plot(μr / (ħ * wz), MT1 ./ GT1, xaxis = :log, yaxis = :log, label = L"12nK", size = (600, 300), line = (:blue), legend = :bottomright)
plot!(μr / (ħ * wz), MT2 ./ GT2, xaxis = :log, yaxis = :log, label = L"30nK", size = (600, 300), line = (:green))
plot!(μr / (ħ * wz), MT3 ./ GT3, xaxis = :log, yaxis = :log, label = L"48nK", size = (600, 300), line = (:red), scale = :log10, minorticks = true, grid = false)
xlabel!(L"\mu/\hbar\omega_x")
ylabel!(L"\Gamma_{M_1}/ \Gamma_\gamma")
savefig("damprateT_comp.pdf")