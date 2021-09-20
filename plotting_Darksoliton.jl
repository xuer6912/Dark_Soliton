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
K2=k2(K)
ψf = xspace(sol[end],sim)
c = sqrt(μ)
ξ = 1/c
v = .05*c
xs = 0
f = sqrt(1-(v/c)^2)
##
fs(vs) = sqrt(1-(vs/c)^2)
ψds(vs) = @. ψf*(fs(vs)*tanh(fs(vs)*(x-xs)/ξ)+im*vs/c)
##

ψ1 = ψds.(0.1*c)
ψ2 = ψds.(0.5*c)
ψ3 = ψds.(0.9*c)
p1 = plot(x,g*abs2.(ψ1),label=L"0.1c")
plot!(x,g*abs2.(ψ2),label=L"0.5c")
plot!(x,g*abs2.(ψ3),label=L"0.9c")
xlabel!(L"x/a_x");ylabel!(L"|\psi|^2")
xlims!(-10,10)
##
p2 = plot(x,angle.(ψ1),label=L"0.1c")
plot!(x,angle.(ψ2),label=L"0.5c")
plot!(x,angle.(ψ3),label=L"0.9c")
xlabel!(L"x/a_x");ylabel!(L"\textrm{phase}(\psi)")
xlims!(-10,10)
ylims!(-0.05,3.5)
p = plot(p1,p2,layout=(2,1),size=(600,400))



##########
 