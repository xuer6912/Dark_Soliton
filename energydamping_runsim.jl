using FourierGPE, VortexDistributions, FFTW
##
#ϕ = fft(ψ)
function diffcurrent2(ϕ,kx)
    ψ = ifft(ϕ)
	ψx = ifft(im*kx.*ϕ)
	j = @. imag(conj(ψ)*ψx)
    jx = im*kx.* fft(j)# current direvative wrt x
	return jx
end
##
import FourierGPE.V
V(x,ψ,kx,t) = 0.5*x^2 
function initsim!(sim;flags=FFTW.MEASURE)
    @unpack L,N = sim
    X,K,dX,dK,DX,DK,T = makearraystransforms(L,N)
    espec = 0.5*k2(L,N)
    @pack! sim = T,X,K,espec
    return nothing
end
##
import FourierGPE.nlim!,FourierGPE.Lgp!
function nlin!(dϕ,ϕ,sim::Sim{1},t)
    @unpack g,X,K,V0 = sim; x = X[1]; kx = K[1]
    dϕ .= ϕ
    xspace!(dϕ,sim)
    Ve = - M*diffcurrent(dϕ,kx)
    @. dϕ *= V0 + V(x,t) + g*abs2(dϕ) + Ve
    kspace!(dϕ,sim)
    return nothing
end
##

function Lgp!(dϕ,ϕ,sim,t)
    @unpack γ,μ,espec = sim
    nlin!(dϕ,ϕ,sim,t)
    @. dϕ = -im*(1.0 - im*γ)*(dϕ + (espec - μ)*ϕ)
    return nothing
end
##
function diffcurrent(ψ,kx)
    ϕ= fft(ψ)
	ψx = ifft(im*kx.*ϕ)
	j = @. imag(conj(ψ)*ψx)
    jx = ifft(im*kx.* fft(j)) # current direvative wrt x
	return real.(jx)
end