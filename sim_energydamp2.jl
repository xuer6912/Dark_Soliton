using Plots, FourierGPE, LaTeXStrings,VortexDistributions, FFTW
gr(fmt="png",legend=true,titlefontsize=12,size=(500,200),grid=false,colorbar=false);
## convenient plotting method
function showpsi(x,ψ)
    p1 = plot(x,g*abs2.(ψ))
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
N = (512,)
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
M = 0.00
sol = runsim(sim);
##
import FourierGPE.nlin!
M = 0.0002
function nlin!(dϕ,ϕ,sim::Sim{1},t)
    @unpack g,X,K,V0 = sim; x = X[1]; kx = K[1]
    dϕ .= ϕ
    xspace!(dϕ,sim)
    Ve =  -M*diffcurrent(dϕ,kx)
    @. dϕ *= V0 + V(x,t) + g*abs2(dϕ) + Ve
    kspace!(dϕ,sim)
    return nothing
end

##
ϕg = sol[end]
ψg = xspace(ϕg,sim)
## Soliton
K2=k2(K)
#psi=XField(ψf,X,K,K2)
#jx = diffcurrent(psi)
ψf = xspace(sol[end],sim)
c = sqrt(μ)
ξ = 1/c
v = .01*c
xs = 0
f = sqrt(1-(v/c)^2)
ψs = @. ψf*(f*tanh(f*(x-xs)/ξ)+im*v/c)
γ = 0; 
#gamma = γ
Nt=800
tf = 16*pi/sqrt(2); t = LinRange(ti,tf,Nt)
dt=diff(t)[1]
ϕi = kspace(ψs,sim)
simSoliton = Sim(sim;γ=γ,tf=tf,t=t,ϕi=ϕi)
@time sols = runsim(simSoliton);
##
ϕf = sols[1]
ψf = xspace(ϕf,simSoliton)
showpsi(x,ψf)
##
kx= K[1]
#plot(x,(diffcurrent(ψf,kx)))
#xlims!(-9,9)
##
#plot(x,J(ψf,kx))

##
dx= diff(x)[1]
dt= diff(t)[1]

xat = zero(t)#(nearest grid point)
xnt = zero(t)
K2=k2(K)
for i in 1:length(t)#make it periodic by ending early
    ψ = xspace(sols[i],sim)
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
ylims!(-5,5)
xlabel!(L"t/t_0")
ylabel!(L"x/x_0")
##
savefig("xs_ed")
##
plot!(t[44:end], xat[44:end], label="analytic",xlims=(0,25),ylims=(-2,2))
#plot(t,xnt)
xi=xat[44]
plot!(t[44:end],xi*exp.(μ*μ/g*2/15*M*t[44:end]),legend=false)
#plot!(t[3:end],xi*exp.(μ*2/15*M*t[3:end]*125),legend=false)
plot!(t[44:end],-xi*exp.(μ*μ/g*2/15*M*t[44:end]),legend=false)

savefig("xs_energydamp")
ts= t
Ekt =  zeros(length(ts))
Ept = zeros(length(ts))
Eit =  zeros(length(ts))
N =  zeros(length(ts))
for i in 1:length(ts) #make it periodic by ending early
    #ψi = ψ0.(x,μ,g)
    ψ = xspace(sols[i],simSoliton)
    psi=XField(ψ,X,K,K2)
    #mask = g*abs2.(ψi).>0.1*μ
    Ekt[i] =sum( Energy(psi)[1])*dx
    Ept[i] = sum( Energy(psi)[2])*dx
    Eit[i] = sum( Energy(psi)[3])*dx
    
    #xlims!(-10,10); ylims!(-1000,1000)
    #title!(L"\textrm{local}\; \mu(x)")
    #xlabel!(L"x/a_x"); ylabel!(L"\mu(x)/\hbar\omega_x")
end

plot(ts,Ekt+Eit,label=L"E_{ki}")
plot(ts,Ept,label=L"E_p")
plot(ts,Ekt+Ept+Eit,label=L"E_{tot}")

savefig("Energy_t_nd.pdf")

for i in 1:length(ts) #make it periodic by ending early
    #ψi = ψ0.(x,μ,g)
    ψ = xspace(sols[i],simSoliton)
    N[i] = sum(abs2.(ψ)-abs2.(ψg))*dx
    #xlims!(-10,10); ylims!(-1000,1000)
    #title!(L"\textrm{local}\; \mu(x)")
    #xlabel!(L"x/a_x"); ylabel!(L"\mu(x)/\hbar\omega_x")
end

plot(N)

R=sqrt(2*μ)
E(x,v) = 4/3 *μ/g*(1-x.^2/R^2)*c*(1-v.^2/c^2).^(3/2)-2*v*μ/g*(1-x.^2/R^2)*ξ*c*(1-v.^2/c^2).^(1/2)
Eki(x,v) = 4/3 *μ/g*(1-x.^2/R^2)*c*(1-v.^2/c^2).^(3/2)
Ep(x,v) = -x.^2*μ/g*(1-x.^2/R^2)*ξ*c*(1-v.^2/c^2).^(1/2)
ΓM = μ*μ/g*2/15*M
ωs = 1/sqrt(2)
xt = 0.084*exp.(ΓM*t) .* sin.(sqrt(ωs^2-ΓM^2)*t)
vt = 0.084*ΓM*exp.(ΓM*t) .* sin.(sqrt(ωs^2-ΓM^2)*t) + 0.084*sqrt(ωs^2-ΓM^2)*exp.(ΓM*t) .* cos.(sqrt(ωs^2-ΓM^2)*t)
plot(t,xt)
plot!(t,xnt)


plot(t,vt)
ind = 400
plot(t[1:ind],E.(xt[1:ind],vt[1:ind]))
plot(ts[1:ind],Ekt[1:ind]+Eit[1:ind],label=L"E_{ki}")
plot!(t[1:ind],Eki.(xt[1:ind],vt[1:ind]))
plot(ts[1:ind],Ept[1:ind],label=L"E_p")
plot!(t[1:ind],Ep.(xt[1:ind],vt[1:ind]))


C(x) = sqrt(μ*(1-x.^2/R^2))

E(x,v) = 4/3 *μ/g*(1-x.^2/R^2)*C.(x) *(1-v.^2 ./C.(x)^2).^(3/2)-2*v*μ/g*(1-x.^2/R^2)*ξ*C.(x)*(1-v.^2 ./C.(x)^2).^(1/2)
Eki(x,v) = 4/3 *μ/g*(1-x.^2/R^2)*C.(x)*(1-v.^2 ./C.(x)^2).^(3/2)
Ep(x,v) = -x.^2*μ/g*(1-x.^2/R^2)*ξ .*C.(x)*(1-v.^2 ./C.(x)^2).^(1/2)