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

function J(ψ,kx)
    ϕ= fft(ψ)
	ψx = ifft(im*kx.*ϕ)
	j = @. imag(conj(ψ)*ψx)
    return real.(j)
end

function diffcurrent(ψ,kx)
    ϕ= fft(ψ)
	ψx = ifft(im*kx.*ϕ)
	j = @. imag(conj(ψ)*ψx)
    jx = ifft(im*kx.* fft(j)) # current direvative wrt x
	return real.(jx)
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
μ =25.0

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
M = 0.0000
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
#ϕg = sol[end]
#ψg = xspace(ϕg,sim)
## Soliton
K2=k2(K)
#psi=XField(ψf,X,K,K2)
#jx = diffcurrent(psi)
ψf = xspace(sol[end],sim)
c = sqrt(μ)
ξ = 1/c
v = .5*c
xs = 0
f = sqrt(1-(v/c)^2)
ψs = @. ψf*(f*tanh(f*(x-xs)/ξ)+im*v/c)
γ = 0; 
#gamma = γ
Nt=800
tf = 8*pi/sqrt(2); t = LinRange(ti,tf,Nt)
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
plot(x,(diffcurrent(ψf,kx)))
xlims!(-9,9)
##
plot(x,J(ψf,kx))

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
plot(t[1:end], xat[1:end], label="analytic",xlims=(0,25),ylims=(-5,5))

##
#plot!(t[2:end], xnt[2:end],label ="numerical",xlims=(0,25),ylims=(-5,5))
xi=xat[3]
#plot!(t[3:end],xi*exp.(μ/3*γ*t[3:end]))
#plot!(t[3:end],-xi*exp.(μ/3*γ*t[3:end]),legend=:false)
##
##
anim = @animate for i in 1:length(t) #make it periodic by ending early
    #ψi = ψ0.(x,μ,g)
    ψ = xspace(sols[i],simSoliton)
    showpsi(x,ψ)
    #plot(x,y)
    #xlims!(-10,10); ylims!(0,1.3*μ)
    #title!(L"\textrm{local}\; \mu(x)")
    #xlabel!(L"x/a_x"); ylabel!(L"\mu(x)/\hbar\omega_x")
end
filename = "decay.gif"
gif(anim,filename,fps=30)

##
diffcurrent(ψg,K[1])

##
anim = @animate for i in 1:length(t) #make it periodic by ending early
    #ψi = ψ0.(x,μ,g)
    ψ = xspace(sols[i],simSoliton)
    plot(x,J(ψ,kx))
    xlims!(-10,10); ylims!(-100,100)
    #title!(L"\textrm{local}\; \mu(x)")
    #xlabel!(L"x/a_x"); ylabel!(L"\mu(x)/\hbar\omega_x")
end
filename = "current.gif"
gif(anim,filename,fps=30)

##
anim = @animate for i in 1:length(t) #make it periodic by ending early
    #ψi = ψ0.(x,μ,g)
    ψ = xspace(sols[i],simSoliton)
    plot(x,diffcurrent(ψ,kx))
    xlims!(-10,10); ylims!(-1000,1000)
    #title!(L"\textrm{local}\; \mu(x)")
    #xlabel!(L"x/a_x"); ylabel!(L"\mu(x)/\hbar\omega_x")
end
filename = "diffcurrent.gif"
gif(anim,filename,fps=30)










