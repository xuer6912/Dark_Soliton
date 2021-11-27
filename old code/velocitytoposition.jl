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
γ = 0; 
xs = 0
Nt=800
tf = 8*pi/sqrt(2); t = LinRange(ti,tf,Nt)
dt=diff(t)[1]
kx= K[1]
dx= diff(x)[1]


#gamma = γ


##
function vscan(vv)
    f = sqrt(1-(vv/c)^2)
    ψs = @. ψf*(f*tanh(f*(x-xs)/ξ)+im*vv/c)
    ϕi = kspace(ψs,sim)
    simSoliton = Sim(sim;γ=γ,tf=tf,t=t,ϕi=ϕi)
    @time sols = runsim(simSoliton);
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
    x1=findmax(xnt)[1]
    x2=findmax(xat)[1]
    plot(t, xat, label="analytic",xlims=(0,25),ylims=(-5,5))
    plot!(t,xnt)
    return x1,x2
end
##

##
#plot!(t[3:end],xi*exp.(μ*2/15*M*t[3:end]*60))
#plot!(t[3:end],-xi*exp.(μ*2/15*M*t[3:end]),legend=:false)
vrange = LinRange(0.01:0.01:0.5) .* c




#
Xa = zero(vrange)
Xn = zero(vrange)
##
for i in 1: length(vrange)
    Xa[i], Xn[i]=vscan(vrange[i])
end
##
plot(Xa)
plot(vrange[1:45]/c, Xn[1:45],legend = false)
plot!(vrange[1:45]/c,vrange[1:45]*sqrt(2))
xlabel!(L"v/c")
ylabel!(L"x/x_0")
savefig("xv1")

plot(vrange[1:40]/c, Xn[1:40],legend = false)
plot!(vrange[1:45]/c,vrange[1:45]*sqrt(2))
xlabel!(L"v/c")
ylabel!(L"x/x_0")
savefig("xv2")