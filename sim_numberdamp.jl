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
v = .1*c
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

γ = 0.01

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
plot(t[3:end], xat[3:end], label="analytic",xlims=(0,25),ylims=(-5,5))
plot(t[2:end], xnt[2:end],label ="numerical",xlims=(0,25),ylims=(-5,5))
xi=xnt[53]
plot!(t[53:end],xi*exp.(μ/3*γ*t[53:end]))
plot!(t[53:end],-xi*exp.(μ/3*γ*t[53:end]),legend=:false)
##


ts = t[1:100]
n = zeros(length(t),N[1])
##
for i in 1:length(t) #make it periodic by ending early
    #ψi = ψ0.(x,μ,g)
    ψ = xspace(sols[i],simSoliton)
    n[i,:] = abs2.(ψ) 
    #xlims!(-10,10); ylims!(-1000,1000)
    #title!(L"\textrm{local}\; \mu(x)")
    #xlabel!(L"x/a_x"); ylabel!(L"\mu(x)/\hbar\omega_x")
end

##
heatmap(t,x,n')
#ylims!(-6,6)
xlabel!(L"t/t_0")
ylabel!(L"x/x_0")
##
plot!(t[2:end], xnt[2:end],label ="numerical")
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
Es_n(ψ) = sum(abs2.(diff(ψ)))/dx + g/2*sum(abs2.(ψ) - abs2.(ψi))*dx

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

Ekt =  zeros(length(t))
Ept = zeros(length(t))
Eit =  zeros(length(t))

for i in 1:length(t) #make it periodic by ending early
    #ψi = ψ0.(x,μ,g)
    ψ = xspace(sols[i],simSoliton)
    psi=XField(ψ,X,K,K2)
    Ekt[i] =sum( Energy2(psi)[1])*dx
    Ept[i] = sum( Energy2(psi)[2])*dx
    Eit[i] = sum( Energy2(psi)[3])*dx
    #xlims!(-10,10); ylims!(-1000,1000)
    #title!(L"\textrm{local}\; \mu(x)")
    #xlabel!(L"x/a_x"); ylabel!(L"\mu(x)/\hbar\omega_x")
end

plot(t,Ekt,label=L"E_k")
plot(t,Ept,label=L"E_p")
plot!(t,Eit,label=L"E_i")
plot!(t,Ekt+Ept+Eit,label=L"E_{tot}")
savefig("Energy_t_nd.pdf")


ϕf = sols[584]
ψf = xspace(ϕf,simSoliton)
plot(abs2.(ψf))

psi=XField(ψf,X,K,K2)

Ek , Ep, Ei = Energy(psi)

##
#plot(Real.(Ek))
plot(x,Ep,label=L"E_p")
plot!(x,Ek,label=L"E_k")
plot!(x,Ei,label=L"E_i")
xlabel!(L"x/x_0")
ylabel!(L"E")
savefig("Energysclice")