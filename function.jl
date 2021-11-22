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












##
anim = @animate for i in 1:length(t)-4 #make it periodic by ending early
    ψ = xspace(sols[i],simSoliton)
    y = g*abs2.(ψ)
    plot(x,y,fill=(0,0.2),size=(500,200))
    xlims!(-10,10); ylims!(0,1.3*μ)
    title!(L"\textrm{local}\; \mu(x)")
    xlabel!(L"x/a_x"); ylabel!(L"\mu(x)/\hbar\omega_x")
end
filename = "darksoliton.gif"
gif(anim,filename,fps=30);
##

anim = @animate for i in 1:length(t)-4 #make it periodic by ending early
    ψ = xspace(sols[i],simSoliton)

    ϕ = angle.(ψ)
    dϕ= diff(ϕ)
    plot(x[1:end-1],dϕ,fill=(0,0.2),size=(500,200))
    xlims!(-10,10)
    title!(L"\textrm{local}\; \mu(x)")
    #xlabel!(L"x/a_x"); ylabel!(L"\mu(x)/\hbar\omega_x")
end
filename = "darksolitonphase.gif"
gif(anim,filename,fps=30);







##
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
##
psi=XField(ψf,X,K,K2)

Ek , Ep, Ei = Energy(psi)

##
#plot(Real.(Ek))


##
Ekt1 =  zeros(length(t),N[1])
Ept1 = zeros(length(t),N[1])
Eit1 =  zeros(length(t),N[1])

for i in 1:length(t) #make it periodic by ending early
    #ψi = ψ0.(x,μ,g)
    ψ = xspace(sols[i],simSoliton)
    psi=XField(ψ,X,K,K2)
    Ekt1[i,:] , Ept1[i,:], Eit1[i,:] = Energy(psi)
    #xlims!(-10,10); ylims!(-1000,1000)
    #title!(L"\textrm{local}\; \mu(x)")
    #xlabel!(L"x/a_x"); ylabel!(L"\mu(x)/\hbar\omega_x")
end

p1 = heatmap(t,x,Ekt1',xlabel=L"t/t_0",ylabel=L"x/x_0",title=L"E_k" )
#xlabel!(L"t/t_0")
#ylabel!(L"x/x_0")
#title!(L"E_k")
p2 = heatmap(t,x,Ept1',xlabel=L"t/t_0",ylabel=L"x/x_0",title=L"E_p" )
p3 = heatmap(t,x,Eit1',xlabel=L"t/t_0",ylabel=L"x/x_0",title=L"E_i" )
p4 = heatmap(t,x,Ekt1'+Ept1'+Eit1',xlabel=L"t/t_0",ylabel=L"x/x_0",title=L"E_{tot}" )
plot(p1,p2,p3,p4,layout=(4,1),size=(600,800))

savefig("Energydistrobution_n.pdf")


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

