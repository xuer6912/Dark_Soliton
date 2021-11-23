function showpsi(x,ψ)
    p1 = plot(x,g*abs2.(ψ))
    xlabel!(L"x/a_x");ylabel!(L"|\psi|^2")
    p2 = plot(x,angle.(ψ))
    xlabel!(L"x/a_x");ylabel!(L"\textrm{phase}(\psi)")
    p = plot(p1,p2,layout=(2,1),size=(600,400))
    return p
end




ff(v)=sqrt(1-v^2/c^2)
momentum_exp(ψ) = im/2.0*(-ψ[1:end-1]'*diff(ψ) + diff(ψ)'*ψ[1:end-1])


function solitonposition(ψ)
    L = length(ψ)
    n1 = (L/2-150)|>Int
    n2 = (L/2+150)|>Int
    ψ2 = ψ[n1:n2]
    n3 = findmin(abs2.(ψ2))[2] + n1 -1
    return n3
end
function soliton_X(sol, sim, x, t)
    f = sqrt(1-(v/c)^2)
    xs = zero.(t)
    #xsa = zero.(t)
    #dx = diff(x)[1]
    for i= 1 : length(t)
        ψ = xspace(sol[i],sim)
        Ms2 = -4.0*abs2.(ψ)/c
        ns = solitonposition(ψ)
        xs[i] = x[ns]
        #xsa[i] =  @. -8 * sum(x .* abs2.(ψ) ./ Ms2 ) /f * dx
        #xsa[i] =  @. 2 * sum(x .* abs2.(ψ) ./ Ms2 ) /f * dx
    end
    return xs
end
function soliton_X2(sol, sim, x, t, xs)
    ff = @. sqrt(1.0 - (v2 ./ c).^2)
    xsa = zero.(t)
    #Ms = - 4.0*μ/g/c
    Ms2 = zero.(t)
    dx = diff(x)[1]
    xe = zero.(t)
    for i = 1 : length(t)
        Ms2[i] = - 4.0*abs2.(ψ0(xs[i],μ,g))/c
        f = ff[i]
        ψ = xspace(sol[i],sim)
        xsa[i] =   2 * sum(x .* abs2.(ψ)  * dx) /(f*Ms2[i]) 
        xe[i] = sum(x .* abs2.(ψ)  * dx)
    end
    return xsa, Ms2, xe
end

function soliton_P(sol, sim, v, t)
    f = @. sqrt(1-(v/c)^2)
    Ps = zero.(t)
    #xsa = zero.(t)
    #dx = diff(x)[1]
    for i= 1 : length(t)
        ψ = xspace(sol[i],sim)
        Ps[i] = momentum_exp(ψ)
        #xsa[i] =  @. -8 * sum(x .* abs2.(ψ) ./ Ms2 ) /f * dx
        #xsa[i] =  @. 2 * sum(x .* abs2.(ψ) ./ Ms2 ) /f * dx
    end
    return Ps
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

#
function Energy(psi::XField{1})
	@unpack psiX,K = psi; kx = K[1]; ψ = psiX
	rho = abs2.(ψ)
    rhoi = abs2.(ψg)
    ψx = gradient(psi::XField{1})
	Ek = 0.5 * abs2.(ψx) 
    Ep = 0.5 * (x.^2) .* (rho - rhoi)
    Ei = 0.5 * g * (rho - rhoi).^2
	return Ek, Ep, Ei
end
##
C(x) = sqrt(μ*(1-x.^2/R^2))

E(x,v) = 4/3 *μ/g*(1-x.^2/R^2)*C.(x) *(1-v.^2 ./C.(x)^2).^(3/2)-2*v*μ/g*(1-x.^2/R^2)*ξ*C.(x)*(1-v.^2 ./C.(x)^2).^(1/2)
Eki(x,v) = 4/3 *μ/g*(1-x.^2/R^2)*C.(x)*(1-v.^2 ./C.(x)^2).^(3/2)
Ep(x,v) = -2*v*μ/g*(1-x.^2/R^2)*ξ .*C.(x)*(1-v.^2 ./C.(x)^2).^(1/2)


