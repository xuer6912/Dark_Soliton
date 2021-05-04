using Plots, FourierGPE, LaTeXStrings,VortexDistributions, FFTW
gr(fmt="png",legend=true,titlefontsize=12,size=(500,200),grid=false,colorbar=false);
## convenient plotting method
function showpsi(x,ψ)
    p1 = plot(x,abs2.(ψ))
    xlabel!(L"x/a_x");ylabel!(L"|\psi|^2")
    p2 = plot(x,angle.(ψ))
    xlabel!(L"x/a_x");ylabel!(L"\textrm{phase}(\psi)")
    p = plot(p1,p2,layout=(2,1),size=(600,400))
    return p
end
function velocity2(psi::XField{1})
	@unpack psiX,K = psi; kx = K[1]; ψ = psiX
	rho = abs2.(ψ)
    ψx = gradient(psi::XField{1})
	vx = @. imag(conj(ψ)*ψx)/rho; @. vx[isnan(vx)] = zero(vx[1])
	return vx
end
S(ψ) =  @. real( -im/2*(log(ψ)-log(conj(ψ))))
DS(ψ) = sum(diff(unwrap(S(ψ))))
##system size

L = (40.0,)
N = (2048,)
sim = Sim(L,N)
@unpack_Sim sim;
μ = 25.0

## Declaring the potential
import FourierGPE.V
V(x,t) = 0.5*x^2

## Initial condition
ψ0(x,μ,g) = sqrt(μ/g)*sqrt(max(1.0-V(x,0.0)/μ,0.0)+im*0.0)
x = X[1];
ψi = ψ0.(x,μ,g)
ϕi = kspace(ψi,sim) #sim uses Fourier transforms that are norm-preserving
@pack_Sim! sim

## harmonic
sol = runsim(sim);
ϕg = sol[end]
ψg = xspace(ϕg,sim)
#Soliton imprinting
ψf = xspace(sol[end],sim)
c = sqrt(μ)
ξ = 1/c
v = 0#.01*c
xs = -0.5
f = sqrt(1-(v/c)^2)
ψs = @. ψf*(f*tanh(f*(x-xs)/ξ)+im*v/c)
xlims!(-10,10)
γ = 0.0
tf = 8*pi/sqrt(2); t = LinRange(ti,tf,Nt)
#dt = 0.01π/μ
ϕi = kspace(ψs,sim)
simSoliton = Sim(sim;γ=γ,tf=tf,t=t,ϕi=ϕi)
@time sols = runsim(simSoliton);
##
ϕf = sols[152]
ψf = xspace(ϕf,simSoliton)
dx=diff(x)[1]
K2=k2(K)
#psi=XField(ψf,X,K,K2)
##numerical method for position
xt=zero(t)#(nearest grid point)
for i in 1:length(t) #make it periodic by ending early
    #ψi = ψ0.(x,μ,g)
    ψ = xspace(sols[i],simSoliton)
    psi=XField(ψ,X,K,K2)
    v=velocity2(psi)
    v_mask = v[g*abs2.(ψi).>0.1*μ]
    x_mask = x[g*abs2.(ψi).>0.1*μ]
    #plot(x,y)
    #plot!(x[g*abs2.(ψi).>0.1*μ],abs2.(dϕ[g*abs2.(ψi[1:end-1]).>0.1*μ]))
    xt[i] = x_mask[findmax(abs.(v_mask))[2]]
end
dt=t[2]-t[1]
plot(t[2:end],xt[2:end])
##
ΔS=zero(t[1:end])
for i in 1:length(t) #make it periodic by ending early
    #ψi = ψ0.(x,μ,g)
    ψd = xspace(sols[i],simSoliton)
    ΔS[i] = DS(ψd[g*abs2.(ψi).>0.1*μ])
    if ΔS[i] < 0
        ΔS[i] += 0
    else
        ΔS[i] += -2*pi
    end
end
##
savefig("velocity comparison")
##
p1=plot(t,ΔS,label=false,size=(600,400),title="phase change")
p2=plot(t,cos.(ΔS/2),label=false,size=(600,400),title="analytic velocity")
p3=plot(t[2:end-1],diff((xt[2:end]))/dt,label=false,size=(600,400), title="numerical velocity")
p4=plot(t[2:end],xt[2:end],size=(600,400))
plot(p1,p2,p3,p4,layout=(4,1))
#title!("velocity vs time")


