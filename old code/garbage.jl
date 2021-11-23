#
using LsqFit
##
@. model(x,p) = p[1]*sin(p[2]*x-p[3])
p0 = [0.5, 0.5]
fit = curve_fit(model, t, xs', p0)
##
plot(t,xs)
ylims!(-0.75,-0.6)
xlims!(6.5,7.2)
#savefig("xs")
##
plot(t,xs)
ylims!(-0.75,-0.6)
xlims!(15.1,15.9)
##

gr(fmt="png",legend=true,titlefontsize=12,size=(500,200),grid=false,colorbar=false)

########
function Energy2(psi::XField{1})
	@unpack psiX,K = psi; kx = K[1]; ψ = psiX
	rho = abs2.(ψ)
    rhoi = abs2.(ψg)
    ψx = gradient(psi::XField{1})
	Ek = 0.5 * abs2.(ψx) 
    Ep = 0.5 * (x.^2) .* (rho )
    Ei = 0.5 * g * (rho- rhoi).^2
	return Ek, Ep, Ei
end
##
##
psi=XField(ψf,X,K,K2)

Ek , Ep, Ei = Energy(psi)






########
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