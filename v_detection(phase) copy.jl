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
#p=plot(x,g*abs2.(ψg),fill=(0,0.2),size=(500,200),label=L"gn(x)")
#plot!(x,one.(x)*μ,label=L"\mu")
#plot!(x,V.(x,0.0),label=L"V(x)",legend=:topright)
#xlims!(-10,10); ylims!(0,1.3*μ)
#title!(L"\textrm{local}\; \mu(x)")
#xlabel!(L"x/a_x"); ylabel!(L"\mu(x)/\hbar\omega_x")
#plot(p)

## Soliton
ψf = xspace(sol[end],sim)
c = sqrt(μ)
ξ = 1/c
v = 0#.01*c
xs = -0.5
f = sqrt(1-(v/c)^2)
ψs = @. ψf*(f*tanh(f*(x-xs)/ξ)+im*v/c)
#showpsi(x,ψs)
xlims!(-10,10)
γ = 0.0
Nt=800
tf = 8*pi/sqrt(2); t = LinRange(ti,tf,Nt)
dt=diff(t)[1]
#dt = 0.01π/μ
ϕi = kspace(ψs,sim)
simSoliton = Sim(sim;γ=γ,tf=tf,t=t,ϕi=ϕi)
##
@time sols = runsim(simSoliton);
##
ϕf = sols[152]
ψf = xspace(ϕf,simSoliton)
#ϕs=angle.(ψf)
#ϕs2=unwrap(ϕs)
dx=diff(x)[1]

anim = @animate for i in 1:length(t)-4 #make it periodic by ending early
    #ψi = ψ0.(x,μ,g)
    ψ = xspace(sols[i],simSoliton)
    y = g*abs2.(ψ)
    ϕ=angle.(ψ)
    ϕ2=unwrap(ϕ)
    dx=diff(x)[1]
    dϕ=diff(ϕ2)/dx
    plot(x,y)
    plot!(x[g*abs2.(ψi).>0.1*μ],abs2.(dϕ[g*abs2.(ψi[1:end-1]).>0.1*μ]))
    xlims!(-10,10); ylims!(0,1.3*μ)
    title!(L"\textrm{local}\; \mu(x)")
    xlabel!(L"x/a_x"); ylabel!(L"\mu(x)/\hbar\omega_x")
end
filename = "darksoliton_x.gif"
gif(anim,filename,fps=30);

##
K

K2=k2(K)
##
psi=XField(ψf,X,K,K2)
##
v=velocity2(psi)
##
vt=zeros(t[1:end-4])
anim = @animate for i in 1:length(t)-4 #make it periodic by ending early
    #ψi = ψ0.(x,μ,g)
    ψ = xspace(sols[i],simSoliton)
    psi=XField(ψ,X,K,K2)
    v=velocity2(psi)
    y = g*abs2.(ψ)
    ϕ=angle.(ψ)
    ϕ2=unwrap(ϕ)
    dx=diff(x)[1]
    dϕ=diff(ϕ2)/dx
    #plot(x,y)
    #plot!(x[g*abs2.(ψi).>0.1*μ],abs2.(dϕ[g*abs2.(ψi[1:end-1]).>0.1*μ]))
    plot(x[g*abs2.(ψi).>0.1*μ],v[g*abs2.(ψi).>0.1*μ])

    xlims!(-10,10); ylims!(-100,100)
    title!("velocity")
    xlabel!(L"x/a_x"); ylabel!(L"v/c")
end
filename = "darksoliton_v.gif"


##
vt=zero(t[1:end])
for i in 1:length(t) #make it periodic by ending early
    #ψi = ψ0.(x,μ,g)
    ψ = xspace(sols[i],simSoliton)
    psi=XField(ψ,X,K,K2)
    v=velocity2(psi)
    y = g*abs2.(ψ)
    ϕ=angle.(ψ)
    ϕ2=unwrap(ϕ)
    dx=diff(x)[1]
    dϕ=diff(ϕ2)/dx
    v_mask = v[g*abs2.(ψi).>0.1*μ]
    #plot(x,y)
    #plot!(x[g*abs2.(ψi).>0.1*μ],abs2.(dϕ[g*abs2.(ψi[1:end-1]).>0.1*μ]))
    vt[i] = v_mask[findmax(abs.(v_mask))[2]]
end
##
plot(t,vt)
##
maximum(v[g*abs2.(ψi).>0.1*μ])
##
v_mask = v[g*abs2.(ψi).>0.1*μ]
v_mask[findmax(abs.(v_mask))[2]]
##
function velocity2(psi::XField{1})
	@unpack psiX,K = psi; kx = K[1]; ψ = psiX
	rho = abs2.(ψ)
    ψx = gradient(psi::XField{1})
	vx = @. imag(conj(ψ)*ψx)/rho; @. vx[isnan(vx)] = zero(vx[1])
	return vx
end

##
plot(x,v)



##
xt=zero(t[1:end])
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
S(ψ) =  @. real( -im/2*(log(ψ)-log(conj(ψ))))
DS(ψ) = sum(diff(unwrap(S(ψ))))
##

ΔS=zero(t[1:end])
for i in 1:length(t) #make it periodic by ending early
    #ψi = ψ0.(x,μ,g)
    ψd = xspace(sols[i],simSoliton)
    ΔS[i] = DS(ψd[g*abs2.(ψi).>0.1*μ])
end

##
p1=plot(t,ΔS,label=false,size=(600,400),title="phase change")
p2=plot(t,cos.(ΔS/2),label=false,size=(600,400),title="analytic velocity")
p3=plot(t[2:end-1],diff((xt[2:end]))/dt,label=false,size=(600,400), title="numerical velocity")

plot(p1,p2,p3,layout=(3,1))
#title!("velocity vs time")
##


##
##
plot(t[3:end],diff(xt[2:end])/dt)
##

S2(ψ) = angle.(ψ)
DS2(ψ) = sum(diff(unwrap(S2(ψ))))
##

ΔS2=zero(t[1:end-4])
for i in 1:length(t)-4 #make it periodic by ending early
    #ψi = ψ0.(x,μ,g)
    ψd = xspace(sols[i],simSoliton)
    ΔS2[i] = DS2(ψd[g*abs2.(ψi).>0.1*μ])
   
end

##
plot(cos.(0.5*ΔS2))
##
p1=plot(x,S(ψf))
##
p2=plot(x,S2(ψf))
plot(p1,p2)
##
p3=plot(x,unwrap(S(ψf)))
##
p4=plot(x,unwrap(S2(ψf)))
plot(p3,p4)
##
p1=plot(imag(-im/2*log.(ψf ./ conj.(ψf))))
minimum(imag(-im/2*log.(ψf ./ conj.(ψf))))
##
##
p5=plot(x[g*abs2.(ψi).>0.01*μ],unwrap(S(ψf[g*abs2.(ψi).>0.01*μ])))
p6=plot(x[g*abs2.(ψi).>0.01*μ],unwrap(S2(ψf[g*abs2.(ψi).>0.01*μ])))
##
plot(p5,p6)

##
plot(real(-im/2*(log.(ψf )-log.( conj.(ψf)))))
##

S3(ψ) =  
##
plot(S3.(ψf))
##
plot(imag.(log.(ψf)))
plot!(imag.(log.(conj.(ψf))))
##
plot(imag.(log.(ψf))-imag.(log.(conj.(ψf))))

##



using JLD2
##
@save "phase.jld2" ψf
##
ψf=nothing
##
@load "phase.jld2" ψf
##
plot(abs2.(ψf))
