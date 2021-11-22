using Plots, FourierGPE, LaTeXStrings
gr(fmt="png",legend=true,titlefontsize=12,size=(500,200),grid=false,colorbar=false);
## convenient plotting method


##system size

L = (40.0,)
N = (1028,)
sim = Sim(L,N)
@unpack_Sim sim;
μ = 25.0

## Declaring the potential
import FourierGPE.V
V(x,t) = 0.5*x^2

## Initial condition
##Let's define a useful Thomas-Fermi wavefunction
ψ0(x,μ,g) = sqrt(μ/g)*sqrt(max(1.0-V(x,0.0)/μ,0.0)+im*0.0)
x = X[1];
ψi = ψ0.(x,μ,g)
ϕi = kspace(ψi,sim) #sim uses Fourier transforms that are norm-preserving
@pack_Sim! sim
## harmonic
sol = runsim(sim);
ϕg = sol[end]
ψg = xspace(ϕg,sim)
p=plot(x,g*abs2.(ψg),fill=(0,0.2),size=(500,200),label=L"gn(x)")
plot!(x,one.(x)*μ,label=L"\mu")
plot!(x,V.(x,0.0),label=L"V(x)",legend=:topright)
xlims!(-10,10); ylims!(0,1.3*μ)
title!(L"\textrm{local}\; \mu(x)")
xlabel!(L"x/a_x"); ylabel!(L"\mu(x)/\hbar\omega_x")
plot(p)


## Soliton
ψf = xspace(sol[end],sim)
c = sqrt(μ)
ξ = 1/c
v = 0#.01*c
xi = -0.5
f = sqrt(1-(v/c)^2)
ψs = @. ψf*(f*tanh(f*(x-xi)/ξ)+im*v/c)
#showpsi(x,ψs)
xlims!(-10,10)
γ = 0.0
tf = 8*pi/sqrt(2); t = LinRange(ti,tf,Nt)
dt = 0.01π/μ
ϕi = kspace(ψs,sim)
simSoliton = Sim(sim;γ=γ,tf=tf,t=t,ϕi=ϕi)
@time sols = runsim(simSoliton);

ϕf = sols[end-4]
ψf = xspace(ϕf,simSoliton)
#showpsi(x,ψf)

##
anim = @animate for i in 1:length(t)-4 #make it periodic by ending early
    ψ = xspace(sols[i],simSoliton)
    y = g*abs2.(ψ)
    plot(x,y,fill=(0,0.2),size=(500,200))
    xlims!(-10,10); ylims!(0,1.3*μ)
    title!(L"\textrm{local}\; \mu(x)")
    xlabel!(L"x/a_x"); ylabel!(L"\mu(x)/\hbar\omega_x")
end
filename = "darksoliton.gif"
gif(anim,filename,fps=30);
##
sum(x .* abs2.(ψs) *diff(x)[1]) 
Ms = -4.0*μ/g/c
-sum(x .* abs2.(ψs) *diff(x)[1]) /(f*Ms/2.0)



##

Ms2 = -4.0*abs2.(ψs)/c
 #-4.0*μ/g/c
##
minimum(Ms2)
##
xs=soliton_X(sols,sim,x,t)
x0=maximum(xs)

plot(t,xs, label = "numerical")
plot!(t,xi*cos.(t/sqrt(2)),label= "analytical")
xlabel!(L"t/t_0")
ylabel!(L"x/\xi") 
title!("soliton position")
##
savefig("solitonposition")
##


##

vn = diff(xi*cos.(t/sqrt(2))/diff(t)[1])
v2 = -xi/sqrt(2)*sin.(t/sqrt(2))
plot(t[1:end-1],vn, label = "density minimum")
plot!(t,va, label = "analytic")
##
xsa,Ms2, xe = soliton_X2(sols,sim,x,t,xs)
plot(t,xs, label = "numerical")
plot!(t,xsa,label= "analytical")
plot!(t,xe/100,label= "expectation")
##
plot(t,Ms2)
##

##Ps= momentum_exp(ψf)
##
Pa = ff.(va).*va/2*Ms/N
Pn = soliton_P(sols, sim, vn, t)/N
plot(t,Pn,label="numerical")
plot!(t,Pa, label="analytic")
xlabel!(L"t")
ylabel!(L"P")
title!("Momentum")
##
savefig("Momentum")


## Energy

KE(ψ) = 0.5*sum(abs2.( diff(ψ))) / diff(x)[1]
PE(ψ) = sum(V.(x,0) .* (abs2.(ψ) .- abs2.(ψi)) * diff(x)[1])
IE(ψ) = 0.5*g*sum( ( abs2.(ψ)-abs2.(ψi) ).^2 ) * diff(x)[1]

##


function Soliton_E(sol, sim, x, t)
    Ek = zero.(t)
    Ep = zero.(t)
    Ei = zero.(t)
    #xsa = zero.(t)
    #dx = diff(x)[1]
    for i= 1 : length(t)
        ψ = xspace(sol[i],sim)
        Ek[i] = KE(ψ)
        Ep[i] = PE(ψ)
        Ei[i] = IE(ψ)
        #xsa[i] =  @. -8 * sum(x .* abs2.(ψ) ./ Ms2 ) /f * dx
        #xsa[i] =  @. 2 * sum(x .* abs2.(ψ) ./ Ms2 ) /f * dx
    end
    return Ek,Ep,Ei
end

##
Ek1,Ep1,Ei1=Soliton_E(sols, sim, x, t)


##
using Statistics
plot(t,Ek1 )
##
plot(t,Ep1)
##
plot(t,Ei1)


















#################


##








#plot(angle.(ψs).-pi/2) 
#ϕ = angle.(ψs).-pi/2
#findmax(ϕ)
#n1=findmin(abs.(ϕ[Int(length(ϕ)/2-200) : Int(length(ϕ)/2+200)]))[2]

##round(length(ϕ)/2-200)|>Int

#ψs[Int(length(ϕ)/2-200)-1+n1]
##

###########
