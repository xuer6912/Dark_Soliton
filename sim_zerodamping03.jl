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
## Soliton
ψf = xspace(sol[end],sim)
c = sqrt(μ)
ξ = 1/c
v = .3*c
xs = 0
f = sqrt(1-(v/c)^2)
ψs = @. ψf*(f*tanh(f*(x-xs)/ξ)+im*v/c)
function nlin!(dϕ,ϕ,sim::Sim{1},t)
    @unpack g,X,K,V0 = sim; x = X[1]; kx = K[1]
    dϕ .= ϕ
    xspace!(dϕ,sim)
    #Ve =  -M*diffcurrent(dϕ,kx)
    @. dϕ *= V0 + V(x,t) + g*abs2(dϕ) #+ Ve
    kspace!(dϕ,sim)
    return nothing
end

γ = 0.00

Nt=800
tf = 16*pi/sqrt(2); t = LinRange(ti,tf,Nt)
dt=diff(t)[1]
ϕi = kspace(ψs,sim)
simSoliton = Sim(sim;γ=γ,tf=tf,t=t,ϕi=ϕi)
@time sols = runsim(simSoliton);

ϕf = sols[152]
ψf = xspace(ϕf,simSoliton)
##
dx= diff(x)[1]
dt= diff(t)[1]
K2=k2(K)
xat = zero(t)#(nearest grid point)
xnt = zero(t)
for i in 1:length(t) #make it periodic by ending early
    ψ = xspace(sols[i],simSoliton)
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
#plot(t[3:end], xat[3:end], label="analytic",xlims=(0,25),ylims=(-5,5))
plot(t[2:end], xnt[2:end],label ="numerical",xlims=(0,25),ylims=(-5,5))
xi=xat[2]
plot!(t[3:end],xi*exp.(μ/3*γ*t[3:end]))
plot!(t[3:end],-xi*exp.(μ/3*γ*t[3:end]),legend=:false)
##




savefig("xs_numberdamp.png")
##
savefig("gamma=$gamma.png")
##
anim = @animate for i in 1:length(t)-4 #make it periodic by ending early
    #ψi = ψ0.(x,μ,g)
    ψ = xspace(sols[i],simSoliton)
    y = g*abs2.(ψ)
    plot(x,y)
    xlims!(-10,10); ylims!(0,1.3*μ)
    title!(L"\textrm{local}\; \mu(x)")
    xlabel!(L"x/a_x"); ylabel!(L"\mu(x)/\hbar\omega_x")
end


filename = "decay.gif"
gif(anim,filename,fps=30)

##


 #make it periodic by ending early
function position_index()
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


Es_a(ψ) = 4/3 * abs2.(ψ) * c - xs^2 *abs2.(ψ) * ξ




##Energy distribution


function Energy(psi::XField{1})
	@unpack psiX,K = psi; kx = K[1]; ψ = psiX
	rho = abs2.(ψ)
    ψx = gradient(psi::XField{1})
	Ek = 0.5 * abs2.(ψx) 
    Ep = 0.5 * (x.^2) .* (rho -ψi.^ 2)
    Ei = 0.5 * g * (rho -ψi.^2)
	return Ek, Ep, Ei
end
##
##
psi=XField(ψf,X,K,K2)

Ek , Ep, Ei = Energy(psi)

##
#plot(Real.(Ek))
plot(Real.(Ep))
plot!(Real.(Ei))

##
Ekt =  zeros(length(t),N[1])
Ept = zeros(length(t),N[1])
Eit =  zeros(length(t),N[1])

for i in 1:length(t) #make it periodic by ending early
    #ψi = ψ0.(x,μ,g)
    ψ = xspace(sols[i],simSoliton)
    psi=XField(ψf,X,K,K2)
    Ekt[i,:] , Ept[i,:], Eit[i,:] = Energy(psi)
    #xlims!(-10,10); ylims!(-1000,1000)
    #title!(L"\textrm{local}\; \mu(x)")
    #xlabel!(L"x/a_x"); ylabel!(L"\mu(x)/\hbar\omega_x")
end

heatmap(t,x,Ekt')
ylims!(-2.5,2.5)


##
ϕ1 = sols[1]
ψ1 = xspace(ϕf,simSoliton)
p1 = plot(x,g*abs2.(ψf), legend = false)
xlims!(-10,10); ylims!(0,1.3*μ)
tn = t[1]
title!("t = $tn ")
xlabel!(L"x/a_x"); ylabel!(L"\mu(x)/\hbar\omega_x")

ψ2 = xspace(sols[51],simSoliton)
p2 = plot(x,g*abs2.(ψ2))
xlims!(-10,10); ylims!(0,1.3*μ)
tn = t[51]
title!("t = $tn ")
xlabel!(L"x/a_x"); ylabel!(L"\mu(x)/\hbar\omega_x")

ψ3 = xspace(sols[101],simSoliton)
p3 = plot(x,g*abs2.(ψ3))
xlims!(-10,10); ylims!(0,1.3*μ)
tn = t[101]
title!("t = $tn ")
xlabel!(L"x/a_x"); ylabel!(L"\mu(x)/\hbar\omega_x")

ψ4 = xspace(sols[152],simSoliton)
p4 = plot(x,g*abs2.(ψ4))
xlims!(-10,10); ylims!(0,1.3*μ)
tn = t[152]
title!("t = $tn ")
xlabel!(L"x/a_x"); ylabel!(L"\mu(x)/\hbar\omega_x")

ψ5 = xspace(sols[199],simSoliton)
p5 = plot(x,g*abs2.(ψ5))
xlims!(-10,10); ylims!(0,1.3*μ)
tn = t[200]
title!("t = $tn ")
xlabel!(L"x/a_x"); ylabel!(L"\mu(x)/\hbar\omega_x")


plot(p1,p2,p3,p4,p5,layout=(5,1),size=(1000,1500), legend = false)
savefig("evolution03")
anim = @animate for i in 1:length(t)-4 #make it periodic by ending early
    #ψi = ψ0.(x,μ,g)
    ψ = xspace(sols[i],simSoliton)
    y = g*abs2.(ψ)
    plot(x,y)
    xlims!(-10,10); ylims!(0,1.3*μ)
    tn = t[i]
    title!("$tn")
    xlabel!(L"x/a_x"); ylabel!(L"\mu(x)/\hbar\omega_x")
end
filename = "z03.gif"
gif(anim,filename,fps=30)



ts = t[1:100]
n = zeros(length(t),N[1])
phase = zeros(length(t),N[1])
##
for i in 1:length(t) #make it periodic by ending early
    #ψi = ψ0.(x,μ,g)
    ψ = xspace(sols[i],simSoliton)
    n[i,:] = abs2.(ψ) 
    phase[i,:] = angle.(ψ)
    #xlims!(-10,10); ylims!(-1000,1000)
    #title!(L"\textrm{local}\; \mu(x)")
    #xlabel!(L"x/a_x"); ylabel!(L"\mu(x)/\hbar\omega_x")
end

##
plot1 = heatmap(t,x,n')
ylims!(-10,10)
xlabel!(L"t/t_0")
ylabel!(L"x/x_0")
title!(L"v=0.3c")

##
#plot!(t[2:end], xnt[2:end],label ="numerical")
#xi=xat[2]
#plot!(t[3:end],xi*exp.(μ/3*γ*t[3:end]))
#plot!(t[3:end],-xi*exp.(μ/3*γ*t[3:end]),legend=:false)
##
plot2 = heatmap(t,x,phase')
ylims!(-10,10)
xlabel!(L"t/t_0")
ylabel!(L"x/x_0")
title!(L"v=0.3c")

plot(plot1,plot2,layout = (2,1),size=(600,800))

savefig("nodamping03")