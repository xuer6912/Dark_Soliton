using Plots, FourierGPE, LaTeXStrings,VortexDistributions
gr(fmt="png",legend=true,titlefontsize=12,size=(500,200),grid=false,colorbar=false);
## convenient plotting method
function showpsi(x,ψ)
    p1 = plot(x,abs2.(ψ))
    xlabel!(L"x/a_x");ylabel!(L"|\psi|^2")
    p2 = plot(x,angle.(ψ))
    xlabel!(L"x/a_x");ylabel!(L"\textrm{phase}(\psi)")
    p = plot(p1,p2,layout=(2,1),size=(600,400))
    return p
end

##system size

L = (40.0,)
N = (2058,)
sim = Sim(L,N)
@unpack_Sim sim;
μ = 25.0

## Declaring the potential
import FourierGPE.V
V(x,t) = 0.5*x^2

## Initial condition
##Let's define a useful Thomas-Fermi wavefunction
ψ0(x,μ,g) = sqrt(μ/g)*sqrt(max(1.0-V(x,0.0)/μ,0.0)+im*0.0)
x = X[1];
ψi = ψ0.(x,μ,g)
ϕi = kspace(ψi,sim) #sim uses Fourier transforms that are norm-preserving
@pack_Sim! sim

## harmonic
sol = runsim(sim);
ϕg = sol[end]
ψg = xspace(ϕg,sim)
p=plot(x,g*abs2.(ψg),fill=(0,0.2),size=(500,200),label=L"gn(x)")
plot!(x,one.(x)*μ,label=L"\mu")
plot!(x,V.(x,0.0),label=L"V(x)",legend=:topright)
xlims!(-10,10); ylims!(0,1.3*μ)
title!(L"\textrm{local}\; \mu(x)")
xlabel!(L"x/a_x"); ylabel!(L"\mu(x)/\hbar\omega_x")
plot(p)



## Soliton
ψf = xspace(sol[end],sim)
c = sqrt(μ)
ξ = 1/c
v = 0#.01*c
xs = -0.5
f = sqrt(1-(v/c)^2)
ψs = @. ψf*(f*tanh(f*(x-xs)/ξ)+im*v/c)
#showpsi(x,ψs)
xlims!(-10,10)
γ = 0.0
tf = 8*pi/sqrt(2); t = LinRange(ti,tf,Nt)
dt = 0.01π/μ
ϕi = kspace(ψs,sim)
simSoliton = Sim(sim;γ=γ,tf=tf,t=t,ϕi=ϕi)
@time sols = runsim(simSoliton);

ϕf = sols[end-4]
ψf = xspace(ϕf,simSoliton)
#showpsi(x,ψf)

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
sum(x .* abs2.(ψs) *diff(x)[1]) 
#
Ms = -4.0*μ/g/c
##

-sum(x .* abs2.(ψs) *diff(x)[1]) /(f*Ms/2.0)



##
ϕs=angle.(ψs)
ϕs2=unwrap(ϕs)
dx=diff(x)[1]
dϕs=diff(ϕs)/dx
dϕs2=diff(ϕs2)/dx
#sum(x[250:750].* abs2.(dϕs[250:750]) *dx)
c=sum( abs2.(dϕs2[520:1550]) *dx)
sum(x[520:1550].* abs2.(dϕs2[520:1550]) *dx)/c
##
plot(x[520:1550],abs2.(dϕs2[520:1550])/c)

##
c=sum( abs2.(dϕs2[250:750]) *dx)



##
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
        #Ms2 = -4.0*abs2.(ψ)/c
        ns = solitonposition(ψ)
        xs[i] = x[ns]
        #xsa[i] =  @. -8 * sum(x .* abs2.(ψ) ./ Ms2 ) /f * dx
        #xsa[i] =  @. 2 * sum(x .* abs2.(ψ) ./ Ms2 ) /f * dx
    end
    return xs
end
##
#Ms2 = -4.0*abs2.(ψs)/c
 -4.0*μ/g/c
##
minimum(Ms2)
##
xs=soliton_X(sols,sim,x,t)
x0=maximum(xs)
##
plot(t,xs, label = "numerical")
plot!(t,-x0*cos.(t/sqrt(2)),label= "analytical")
xlabel!(L"t/t_0")
ylabel!(L"x/\xi")
title!("soliton position")
##
savefig("solitonposition")
##
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
##
xsa,Ms2, xe = soliton_X2(sols,sim,x,t,xs)
plot(t,xs, label = "numerical")
#plot!(t,xsa,label= "analytical")
#plot!(t,xe/100,label= "expectation")
##
plot(t,Ms2)

##
ff = @. sqrt(1.0 - (v2 ./ c).^2)
plot(t,ff)
##
v1 = diff(-x0*cos.(t/sqrt(2))/diff(t)[1])
v2 = x0/sqrt(2)*sin.(t/sqrt(2))
plot(t[1:end-1],v1)
plot!(t,v2)
##



#################

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





#plot(angle.(ψs).-pi/2) 
#ϕ = angle.(ψs).-pi/2
#findmax(ϕ)
#n1=findmin(abs.(ϕ[Int(length(ϕ)/2-200) : Int(length(ϕ)/2+200)]))[2]

##round(length(ϕ)/2-200)|>Int

#ψs[Int(length(ϕ)/2-200)-1+n1]
##

###########
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