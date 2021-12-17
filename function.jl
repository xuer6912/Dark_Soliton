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

function diffcurrent(ψ,kx)
    ϕ= fft(ψ)
	ψx = ifft(im*kx.*ϕ)
	j = @. imag(conj(ψ)*ψx)
    
    jx = ifft(im*kx.* fft(j)) # current direvative wrt x
	return real.(jx)
end

function current(ψ,kx)
    ϕ= fft(ψ)
	ψx = ifft(im*kx.*ϕ)
	j = @. imag(conj(ψ)*ψx)
     # current direvative wrt x
	return j
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
    psig=XField(ψg,X,K,K2)
    ψx = gradient(psi::XField{1})
    ψgx = gradient(psig::XField{1})
	Ek = 0.5 * abs2.(ψx) - 0.5 * abs2.(ψgx) 
    Ep = 0.5 * (x.^2) .* (rho - rhoi)
    Ei = 0.5 * g * (rho - rhoi).^2
    Eki = Ek+Ei
	return Eki, Ep
end

S(ψ) =  @. real( -im/2*(log(ψ)-log(conj(ψ))))# phase of wave

DS(ψ) = sum(diff(unwrap(S(ψ)))) #phase change
##
##

#E(x,v) = 4/3 *μ/g*(1-x.^2/R^2)*c*(1-v.^2/c^2).^(3/2)-2*v*μ/g*(1-x.^2/R^2)*ξ*c*(1-v.^2/c^2).^(1/2)
#Eki(x,v) = 4/3 *μ/g*(1-x.^2/R^2)*c*(1-v.^2/c^2).^(3/2)
#Ep(x,v) = -x.^2*μ/g*(1-x.^2/R^2)*ξ*c*(1-v.^2/c^2).^(1/2)
xat(ΓM,t) = vi*sqrt(2)*exp.(ΓM*t) .* sin.(sqrt(ωs^2-ΓM^2)*t)
vat(ΓM,t) = vi*sqrt(2)*ΓM*exp.(ΓM*t) .* sin.(sqrt(ωs^2-ΓM^2)*t) + vi*sqrt(2)*sqrt(ωs^2-ΓM^2)*exp.(ΓM*t) .* cos.(sqrt(ωs^2-ΓM^2)*t)
C(x) = sqrt(μ*(1-x.^2/R^2))
#E(x,v) = 4/3 *μ/g*(1-x.^2/R^2)*C.(x) *(1-v.^2 ./C.(x)^2).^(3/2)-x.^2*μ/g*(1-x.^2/R^2)*ξ*C.(x)*(1-v.^2 ./C.(x)^2).^(1/2)
Eki(x,v) = 4/3 *μ/g*(1-x.^2/R^2)*C.(x)*(1-v.^2 ./C.(x)^2).^(3/2)
Ep(x,v) = -x.^2*μ/g*(1-x.^2/R^2)*ξ *(1-v.^2 ./C.(x)^2).^(1/2)
E(x,v) = Eki(x,v)+ Ep(x,v)


function solitondynamics(sols,sim,t)
ts=t
#xft = zero(t)#(nearest grid point)
xnt = zero(t)
xt2 = zero(t)
Ekit =  zeros(length(ts))
Ept = zeros(length(ts))
#Eit =  zeros(length(ts))
Ns =  zeros(length(ts))
for i in 1:length(t)
    ψ = xspace(sols[i],sim)
    psi=XField(ψ,X,K,K2)
    mask = g*abs2.(ψi).>0.1*μ
    v=velocity2(psi)
    vm = v[mask]
    xm = x[mask]
    xnt[i] = xm[findmax(abs.(vm))[2]]
    if i == 1
        xt2[i] = xnt[1]
    else
        j = findall(a->a==xt2[i-1],x)[1] 
        mask2 = zero.(ψ)
        mask2[j-2:j+2] .= 1
        xt2[i] = xm[findmax(abs.(vm .* mask2[mask]))[2]]
    end

    #xnt2[i] = xm[findmin(abs2.(ψm))[2]]
    Ekit[i] =sum( Energy(psi)[1])*dx
    Ept[i] = sum( Energy(psi)[2])*dx
    #Eit[i] = sum( Energy(psi)[3])*dx
    Ns[i] = sum(abs2.(ψ)-abs2.(ψg))*dx
end
    return xt2,Ekit,Ept,Ns
end

function Energy2(psi::XField{1})
	@unpack psiX,K = psi; kx = K[1]; ψ = psiX
	rho = abs2.(ψ)
    #rhoi = abs2.(ψg)
    #psig=XField(ψg,X,K,K2)
    ψx = gradient(psi::XField{1})
    #ψgx = gradient(psig::XField{1})
	Ek = 0.5 * abs2.(ψx) #- 0.5 * abs2.(ψgx) 
    Ep = 0.5 * (x.^2) .* (rho )
    Ei = 0.5 * g * rho .^2
	return Ek, Ep, Ei
end


function solitondynamics2(sols,sim,t)
    ts=t
    Ekt =  zeros(length(ts))
    Ept = zeros(length(ts))
    Eit =  zeros(length(ts))
    dJt = zeros(length(ts))
    #Eit =  zeros(length(ts))
    #Ns =  zeros(length(ts))
    for i in 1:length(t)
        ψ = xspace(sols[i],sim)
        psi=XField(ψ,X,K,K2)
        Ekt[i] =sum( Energy2(psi)[1])*dx
        Ept[i] = sum( Energy2(psi)[2])*dx
        Eit[i] = sum( Energy2(psi)[3])*dx
        dJt[i] = sum(diffcurrent(ψ,kx).^2)*dx
    end
    ET = Ekt+Ept+Eit
    return ET,dJt
end