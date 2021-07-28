using FourierGPE, VortexDistributions, FFTW
##
function diffcurrent(ψ,kx)
    ϕ = fft(ψ)
	ψx = ifft(im*kx.*ϕ)
	j = @. imag(conj(ψ)*ψx)
    jx = ifft(im*kx.* fft(j)) # current direvative wrt x
	return jx
end
##
import FourierGPE.V
V(x,ψ,kx,t) = 0.5*x^2 - M*diffcurrent(ψ,kx) 
function initsim!(sim;flags=FFTW.MEASURE)
    @unpack L,N = sim
    X,K,dX,dK,DX,DK,T = makearraystransforms(L,N)
    espec = 0.5*k2(L,N)
    @pack! sim = T,X,K,espec
    return nothing
end
##
function nlin2!(dϕ,ϕ,sim::Sim{1},t)
    @unpack g,X,K,V0 = sim; x = X[1]; kx = K[1]
    dϕ .= ϕ
    xspace!(dϕ,sim)
    @. dϕ *= V0 + V(x,ψ,kx,t) + g*abs2(dϕ)
    kspace!(dϕ,sim)
    return nothing
end
##
function Lgp2!(dϕ,ϕ,sim,t)
    @unpack γ,μ,espec = sim
    nlin2!(dϕ,ϕ,sim,t)
    @. dϕ = -im*(1.0 - im*γ)*(dϕ + (espec - μ)*ϕ)
    return nothing
end
##
function runsim2(sim,ϕ=sim.ϕi;info=true,tplot=false,nfiles=false)
    @unpack nfiles,path,filename = sim

    function savefunction(ψ...)
        isdir(path) || mkpath(path)
        i = findfirst(x->x== ψ[2],sim.t)
        padi = lpad(string(i),ndigits(length(sim.t)),"0")
        info && println("⭆ Save $i at t = $(trunc(ψ[2];digits=3))")
        # tofile = path*"/"*filename*padi*".jld2"
        tofile = joinpath(path,filename*padi*".jld2")
        save(tofile,"ψ",ψ[1],"t",ψ[2])
    end

    savecb = FunctionCallingCallback(savefunction;
                     funcat = sim.t, # times to save at
                     func_everystep=false,
                     func_start = true,
                     tdir=1)

    prob = ODEProblem(Lgp2!,ϕ,(sim.ti,sim.tf),sim)
    info && @info "⭆ 𝒅𝜳 Evolving in kspace"
    info && @info "⭆ Damping γ = $(sim.γ)"
    (info && nfiles) && @info "⭆ Saving to "*path
    nfiles ?
    (sol = solve(prob,alg=sim.alg,saveat=sim.t[end],reltol=sim.reltol,callback=savecb,dense=false,maxiters=1e10,progress=true)) :
    (sol = solve(prob,alg=sim.alg,saveat=sim.t,reltol=sim.reltol,dense=false,maxiters=1e10,progress=true))
    info && @info "⭆ Finished."
return sol
end