using Plots, FourierGPE, LaTeXStrings, VortexDistributions, FFTW
gr(fmt = "pdf", legend = true, display = true, titlefontsize = 12, size = (500, 200), grid = false, colorbar = false);
## convenient plotting method
plot(rand(10, 10))


##system size
ωs = 1 / sqrt(2)
L = (40.0,)
N = (512,)
sim = Sim(L, N)
@unpack_Sim sim;
μ = 25.0
R = sqrt(2 * μ)
include("function.jl")
## Declaring the potential
import FourierGPE.V
V(x, t) = 0.5 * x^2

## Initial condition
ψ0(x, μ, g) = sqrt(μ / g) * sqrt(max(1.0 - V(x, 0.0) / μ, 0.0) + im * 0.0)
x = X[1];
ψi = ψ0.(x, μ, g)
ϕi = kspace(ψi, sim) #sim uses Fourier transforms that are norm-preserving
@pack_Sim! sim

## harmonic
M = 0.00
sol = runsim(sim);
ϕg = sol[end]
ψg = xspace(ϕg, sim)

K2 = k2(K)
dx = diff(x)[1]
dt = diff(t)[1]
##
c = sqrt(μ)
ξ = 1 / c
vi = 0.1 * c
xs = 0
f = sqrt(1 - (vi / c)^2)

Nt = 800
tf = 16 * pi / sqrt(2);
t = LinRange(ti, tf, Nt);
ψs = @. ψg * (f * tanh(f * (x - xs) / ξ) + im * vi / c)
ϕi = kspace(ψs, sim)
#ψf = xspace(sol[end],sim)


γ = 0;
M = 0.000
simSoliton = Sim(sim, γ = γ, tf = tf, t = t, ϕi = ϕi)
@time sols = runsim(simSoliton);
##

ϕf = sols[50]
ψf = xspace(ϕf, simSoliton)



psi_0 = XField(ψf, X, K, K2)
plot(x, velocity2(psi_0))

#plot(xnt)
Ms = [0.0000]
ΓMs = μ * μ / g * 2 / 15 .* Ms
##
ts = t[1:300]
xnt, Ekit, Ept, Ns = solitondynamics(sols, sim, ts)
#xnt2,Ekit2,Ept2,Ns2 = solitondynamics(sols2,sim,ts)
#xnt3,Ekit3,Ept3,Ns3 = solitondynamics(sols3,sim,ts)


p1 = plot(ts, xnt, legend = :bottomright, label = L"x_n")
plot!(ts, xat(ΓMs[1], ts), label = L"x_a")
plot!(ts, vi * sqrt(2) * exp.(ΓMs[1] * ts), label = L"e^{\Gamma t}")
xlabel!(L"t/t_0")
ylabel!(L"x/x_0")
ylims!(-R, R)
#hline!([R/2],label=L"R/2")
#title!(L"M = 0.0001")
plot(ts, xnt, legend = :topleft, label = "M = 0.0001")
plot!(ts, xnt2, label = "M = 0.00015")
plot!(ts, xnt3, label = "M = 0.0002")
xlabel!(L"t/t_0")
ylabel!(L"x/x_0")
ylims!(-10, 10)
#savefig("position_num.pdf")
plot(ts, xat(ΓMs[1], ts), legend = :topleft, label = "M = 0.0001")
#plot!(ts,xat(ΓMs[2],ts),label="M = 0.00015")
#plot!(ts,xat(ΓMs[3],ts),label="M = 0.0002")
xlabel!(L"t/t_0")
ylabel!(L"x/x_0")
#savefig("position_ana.pdf")
plot(ts, vat(ΓMs[1] / c, ts), legend = :bottomright, label = "M = 0.0001")
xlabel!(L"t/t_0")
ylabel!(L"v/c")

savefig("velocity_ana.pdf")



#########energy numeric








#savefig("Etot.pdf")
plot(p1, p2, p3, layout = (3, 1), size = (700, 800))
savefig("Energy_num.pdf")

######energy analytic
ts = t[1:300]
xt1 = xat(ΓMs[1], ts)

vt1 = vat(ΓMs[1], ts)
Ept, Ept_e = trapdynamics(sols, sim, ts)

plot(ts, Ept / μ, legend = :bottomright, label = "numerical")
#plot!(ts, Ep.(xt1, vt1) / μ, label = "analytic")
plot!(ts, (cos.(ts * (1 + 1 / sqrt(2))) .- 1) / 2*0.65, label = L"\omega+\omega_s")
#plot!(ts, -Ept_e,label="expectation")
xlabel!(L"t/t_0")
ylabel!(L"E_{p}/\mu")
savefig("guess.pdf")



######compare E




p1 = plot(ts, Eki.(xt1, vt1) / μ, legend = :bottomleft, label = "analytic")
plot!(ts, Ekit / μ, label = "numeric")
p2 = plot(ts, Ep.(xt1, vt1) / μ, legend = :bottomleft, label = "analytic")
plot!(ts, Ept / μ, label = "numeric")
p3 = plot(ts, E.(xt1, vt1) / μ, legend = :bottomleft, label = "analytic")
plot!(ts, Ept / μ + Ekit / μ, label = "numeric")

plot(p1, p2, p3, layout = (3, 1), size = (800, 800))
xlabel!(L"t/t_0")
ylabel!(L"E/\mu")
title!(L"M=0.0001")
savefig("Energy_comp.pdf")

v_scan = LinRange(0, 1, 100) * c
p1 = plot(v_scan / c, E.(0, v_scan) / μ, label = L"x_i = 0")
xlabel!(L"v/c")
ylabel!(L"E/\mu")


x_scan = LinRange(0, R, 100)
p2 = plot(x_scan / R, E.(x_scan, 0) / μ, label = L"v_i = 0")
xlabel!(L"x/R")
ylabel!(L"E/\mu")

plot(p1, p2, layout = (1, 2), size = (600, 250))
savefig("solitonEnergy.pdf")





#########
p0 = plot(ts, xat(ΓMs[1], ts), label = L"x_a", legend = :bottomright)
plot(ts, xat(ΓMs[1], ts) .^ 2, label = L"x_a^2")
plot!(ts, xat(ΓMs[1], ts) .^ 4, label = L"v_a")
plot!(ts, Ep.(xt1, vt1) / μ, label = "analytic")
plot!(ts, Ept / μ, label = "numeric")








ET, dJt = solitondynamics2(sols, sim, ts)
ET2, dJt2 = solitondynamics2(sols2, sim, ts)
ET3, dJt3 = solitondynamics2(sols3, sim, ts)
plot1 = plot(ts, ET / μ, legend = :bottomleft, label = L"M=0.0001")
plot!(ts, ET2 / μ, label = L"M=0.00015")
plot!(ts, ET3 / μ, label = L"M=0.0002")
xlabel!(L"t/t_0")
ylabel!(L"E_{tot}/\mu")


plot2 = plot(ts, xt1 .^ 2, legend = :topleft, label = L"M=0.0001")
plot!(ts, xt2 .^ 2, label = L"M=0.00015")
plot!(ts, xt3 .^ 2, label = L"M=0.0002")
xlabel!(L"t/t_0")
ylabel!(L"x^2/x_0^2")


plot3 = plot(ts, dJt, legend = :topleft, label = L"M=0.0001")
plot!(ts, dJt2, label = L"M=0.00015")
plot!(ts, dJt3, label = L"M=0.0002")
xlabel!(L"t/t_0")
ylabel!(L"\int dx|\nabla_x J|^2")

plot(plot1, plot2, plot3, layout = (3, 1), size = (600, 600))
savefig("total_energy_damp.pdf")






ψ1 = xspace(sols[1], simSoliton)
psi = XField(ψ1, X, K, K2)
Ekit = Energy(psi)[1]
Ept = Energy(psi)[2]
p1 = plot(x, Ept / μ)
xlims!(-10, 10)
tn = t[1]
title!("t = $tn ")
xlabel!(L"x/a_x");
ylabel!(L"E_p/\mu");

ψ2 = xspace(sols[51], simSoliton)
psi = XField(ψ2, X, K, K2)
Ekit = Energy(psi)[1]
Ept = Energy(psi)[2]
p2 = plot(x, Ept / μ)
xlims!(-10, 10)
tn = t[51]
title!("t = $tn ")
xlabel!(L"x/a_x");
ylabel!(L"E_p/\mu");

ψ3 = xspace(sols[101], simSoliton)
psi = XField(ψ3, X, K, K2)
Ekit = Energy(psi)[1]
Ept = Energy(psi)[2]
p3 = plot(x, Ept / μ)
xlims!(-10, 10)
tn = t[101]
title!("t = $tn ")
xlabel!(L"x/a_x");
ylabel!(L"E_p/\mu");

ψ4 = xspace(sols[152], simSoliton)
psi = XField(ψ4, X, K, K2)
Ekit = Energy(psi)[1]
Ept = Energy(psi)[2]
p4 = plot(x, Ept / μ)
xlims!(-10, 10)
tn = t[152]
title!("t = $tn ")
xlabel!(L"x/a_x");
ylabel!(L"E_p/\mu");

ψ5 = xspace(sols[199], simSoliton)
psi = XField(ψ5, X, K, K2)
Ekit = Energy(psi)[1]
Ept = Energy(psi)[2]
p5 = plot(x, Ept / μ)
xlims!(-10, 10)
tn = t[200]
title!("t = $tn ")
xlabel!(L"x/a_x");
ylabel!(L"E_p/\mu");

plot(p2, p3, p4, p5, layout = (4, 1), size = (1000, 1000), legend = false)
savefig("potential_dist.pdf")

anim = @animate for i in 1:length(t)-4 #make it periodic by ending early
    #ψi = ψ0.(x,μ,g)
    ψ = xspace(sols[i], simSoliton)
    psi = XField(ψ, X, K, K2)
    Ept = Energy(psi)[2]
    plot(x, Ept)
    xlims!(-10, 10)
    ylims!(-1250, 500)
    tn = t[i]
    title!("$tn")
    xlabel!(L"x/a_x")
end
filename = "pt.gif"
gif(anim, filename, fps = 30)