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

S(ψ) =  @. real( -im/2*(log(ψ)-log(conj(ψ))))# phase of wave

DS(ψ) = sum(diff(unwrap(S(ψ)))) #phase change
##
#new method for velocity
S(ψ) =  @. real( -im/2*(log(ψ)-log(conj(ψ))))# phase of wave

DS2(ψ) = sum(diff(unwrap(S(ψ))) .*abs2.(ψi[1:end-1]))/(μ/g)
##system size

L = (80.0,)
N = (2048,)
sim = Sim(L,N)
@unpack_Sim sim;
μ = 25

## Declaring the potential
import FourierGPE.V
V(x,t) = 0.5*x^2

## Initial condition
ψ0(x,μ,g) = sqrt(μ/g)*sqrt(max(1.0-V(x,0.0)/μ,0.0)+im*0.0)
x = X[1];
ψi = ψ0.(x,μ,g)
ϕi = kspace(ψi,sim) #sim uses Fourier transforms that are norm-preserving
@pack_Sim! sim

# harmonic
sol = runsim(sim);
ϕg = sol[end]
ψg = xspace(ϕg,sim)
##
#Soliton imprinting &damping
ψf = xspace(sol[end],sim)
c = sqrt(μ)
ξ = 1/c
v = .01*c
xs = 0
f = sqrt(1-(v/c)^2)
ψs = @. ψf*(f*tanh(f*(x-xs)/ξ)+im*v/c)
γ = 0.1
tf = pi/sqrt(2); t = LinRange(ti,tf,Nt)
#dt = 0.01π/μ
ϕi = kspace(ψs,sim)
simSoliton = Sim(sim;γ=γ,tf=tf,t=t,ϕi=ϕi)
@time sols = runsim(simSoliton);



##Soliton evolution
ϕs = sols[end]
ψs = xspace(ϕs,simSoliton)
plot(abs2.(ψs))
##
γ = 0.
Nt=400
tf = 16*pi/sqrt(2); t = LinRange(ti,tf,Nt)
#dt = 0.01π/μ
ϕi = kspace(ψs,sim)
simSoliton = Sim(sim;γ=γ,tf=tf,t=t,ϕi=ϕi)
@time sols2 = runsim(simSoliton);
##
ϕf = sols2[152]
ψf = xspace(ϕf,simSoliton)
dx=diff(x)[1]
dt=diff(t)[1]
K2=k2(K)
#psi=XField(ψf,X,K,K2)
##numerical method for position
xt=zero(t)#(nearest grid point)
for i in 1:length(t) #make it periodic by ending early
    #ψi = ψ0.(x,μ,g)
    ψ = xspace(sols2[i],simSoliton)
    psi=XField(ψ,X,K,K2)
    v=velocity2(psi)
    v_mask = v[g*abs2.(ψi).>0.1*μ]
    x_mask = x[g*abs2.(ψi).>0.1*μ]
    #plot(x,y)
    #plot!(x[g*abs2.(ψi).>0.1*μ],abs2.(dϕ[g*abs2.(ψi[1:end-1]).>0.1*μ]))
    xt[i] = x_mask[findmax(abs.(v_mask))[2]]
end

##
psi=XField(ψf,X,K,K2)
v=velocity2(psi)
plot(v[g*abs2.(ψi).>0.1*μ])
##
sum(v[g*abs2.(ψi).>0.1*μ])*dx/c
##
plot(t[2:end],xt[2:end])
##
##
p1 = plot(unwrap(angle.(ψf[g*abs2.(ψi).>0.1*μ])),label=false,title="numerical phase")
p2 = plot(unwrap(S(ψf[g*abs2.(ψi).>0.1*μ])),label=false,title="analytical phase")
plot(p1,p2, layout=(2,1))

##
p1 = plot(abs2.(ψf[g*abs2.(ψi).>0.1*μ]),label=false,title="density")
p2 = plot(unwrap(S(ψf[g*abs2.(ψi).>0.1*μ])),label=false,title="analytical phase")
ϕ =    unwrap(S(ψf[g*abs2.(ψi).>0.1*μ]))
p3 = plot(diff(ϕ)/dx,label=false,title="phase gradient")
plot(p1,p2,p3,layout=(3,1),size=(600,500))

##
ΔS = DS(ψf[g*abs2.(ψi).>0.1*μ])
dϕ = diff(ϕ)/(dx*ΔS)
xm=x[g*abs2.(ψi).>0.1*μ]
plot(xm[1:end-1],dϕ)
##
xa = xm[1:end-1]'*dϕ *dx
##
xat = zero(t)#(nearest grid point)
xnt = zero(t)
for i in 1:length(t) #make it periodic by ending early
    ψ = xspace(sols2[i],simSoliton)
    psi=XField(ψ,X,K,K2)
    v=velocity2(psi)
    mask = g*abs2.(ψi).>0.1*μ
    v_mask = v[mask]
    x_mask = x[mask]
    xnt[i] = x_mask[findmax(abs.(v_mask))[2]]
    ΔS = DS(ψ[mask])
    ϕ =    unwrap(S(ψ[mask]))
    dϕ = diff(ϕ)/(dx*ΔS)
    xm = x[mask]
    xat[i] = xm[1:end-1]'*dϕ *dx
  
end

##
xmax=findmax(xnt)[1]
plot(t[3:end], xat[3:end], label="analytic")
#plot!(t[2:end], xnt[2:end],label ="numerical")
#plot!(t,xmax*sin.(t/sqrt(2)))
##
savefig("position comp")

##
N = sum(abs2.(ψi)*dx)
##
ΔSt2=zero(t)
for i in 1:length(t) #make it periodic by ending early
    #ψi = ψ0.(x,μ,g)
    ψd = xspace(sols2[i],simSoliton)
    ΔSt2[i] = DS2(ψd)
    if ΔSt2[i] < 0
        ΔSt2[i] += 0
    else
        ΔSt2[i] += -2*pi
    end

end
##
ΔSt=zero(t[1:end])
for i in 1:length(t) #make it periodic by ending early
    #ψi = ψ0.(x,μ,g)
    ψd = xspace(sols2[i],simSoliton)
    ΔSt[i] = DS(ψd[g*abs2.(ψi).>0.001*μ])
    if ΔSt[i] < 0
        ΔSt[i] += 0
    else
        ΔSt[i] += -2*pi
    end
end
##

p1=plot(t,ΔSt,label=false,size=(600,400),title="phase change")
p2=plot(t,cos.(ΔSt/2),label=false,size=(600,400),title="analytic velocity")
p3=plot(t,cos.( ΔSt2/2))
p4=plot(t[2:end-1],diff((xt[2:end]))/dt,label=false,size=(600,400), title="numerical velocity")

plot(p1,p2,p3,p4,layout=(4,1))
##

#title!("velocity vs time")
savefig("velocity comparison")


##
a=findmin(vat[1:100])[2] ;b=findmin(vat[101:200])[2]
period = (t[100+b-1]-t[a-1])/(2*pi)
##

anim = @animate for i in 1:length(t)-4 #make it periodic by ending early
    #ψi = ψ0.(x,μ,g)
    ψ = xspace(sols2[i],simSoliton)

    mask = g*abs2.(ψi).>0
    a = (1:length(mask))[mask]
    mask2 = a[1]-25:a[end]+25
    ϕ1 = unwrap(angle.(ψ[mask2]))
    ϕ2 = unwrap(S(ψ[mask2]))
    p1=plot(x[mask2],ϕ1)
    #plot!(x[mask2],ϕ2)
    plot!(x[g*abs2.(ψi).>0.1*μ],abs2.(dϕ[g*abs2.(ψi[1:end-1]).>0.1*μ]))
    xlims!(-10,10); ylims!(-2*pi,2*pi)
    #title!(L"\textrm{local}\; \mu(x)")
    xlabel!(L"x/a_x"); ylabel!(L"phase")
    y = g*abs2.(ψ)
    p2 = plot(x,y)
    xlims!(-10,10); ylims!(0,1.3*μ)
    title!(L"\textrm{local}\; \mu(x)")
    xlabel!(L"x/a_x"); ylabel!(L"\mu(x)/\hbar\omega_x")
    plot(p1,p2,layout=(2,1),size=(400,400))
end
filename = "phase.gif"
gif(anim,filename,fps=30)
##
p1=plot((@.real(f*tanh(f*(x-xs)/ξ)+im*v/c)))
plot!(abs2.(ψi)/(μ/g))
p2=plot((@.imag(f*tanh(f*(x-xs)/ξ)+im*v/c)))
plot(p1,p2,layout=(2,1),size=(400,400))
##
mask = g*abs2.(ψi).>0
a = (1:length(mask))[mask]


##
anim = @animate for i in 1:length(t)-4 #make it periodic by ending early
    #ψi = ψ0.(x,μ,g)
    ψ = xspace(sols2[i],simSoliton)
    
    mask = g*abs2.(ψi).>0.1*μ
    a = findmin(abs2.(ψ[mask]))[2]
    ϕ = unwrap(angle.(ψ[mask]))
    xmask = x[mask]
    p1=plot(xmask[a-100:a+100],ϕ[a-100:a+100])
    xlims!(-10,10); ylims!(-2*pi,2*pi)

    xlabel!(L"x/a_x"); ylabel!(L"phase")
 
end
filename = "phase2.gif"
gif(anim,filename,fps=30)
##
mask = g*abs2.(ψi).>0.1*μ
a = findmin(abs2.(ψf[mask]))[2]


ΔSt3=zero(t[1:end])
for i in 1:length(t) #make it periodic by ending early
    #ψi = ψ0.(x,μ,g)
    ψ = xspace(sols2[i],simSoliton)
    mask = g*abs2.(ψi).>0.1*μ
    ψm = ψ[mask]
    a = findmin(abs2.(ψm))[2]
    #ϕ = unwrap(angle.(ψ[mask]))
    xmask = x[mask]

    ΔSt3[i] = DS(ψm[a-50:a+50])
    if ΔSt3[i] < 0
        ΔSt3[i] += 0
    else
        ΔSt3[i] += -2*pi
    end
end
##

p1=plot(t,cos.(ΔSt/2),label=false,size=(600,400),title=L"v_a(1)")
p2=plot(t,cos.(ΔSt3/2),label=false,size=(600,400),title="symmetric sampling")
p3=plot(t,cos.( ΔSt2/2),label=false,title=L"(n0)")
p4=plot(t[2:end-1],diff((xt[2:end]))/dt,label=false,size=(600,400), title="numerical velocity")

plot(p1,p2,p3,p4,layout=(4,1))
##
savefig("velocity comp_3")
##