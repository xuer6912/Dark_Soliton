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

#E(x,v) = 4/3 *μ/g*(1-x.^2/R^2)*c*(1-v.^2/c^2).^(3/2)-2*v*μ/g*(1-x.^2/R^2)*ξ*c*(1-v.^2/c^2).^(1/2)
#Eki(x,v) = 4/3 *μ/g*(1-x.^2/R^2)*c*(1-v.^2/c^2).^(3/2)
#Ep(x,v) = -2*v*μ/g*(1-x.^2/R^2)*ξ*c*(1-v.^2/c^2).^(1/2)

C(x) = sqrt(μ*(1-x.^2/R^2))
E(x,v) = 4/3 *μ/g*(1-x.^2/R^2)*C.(x) *(1-v.^2 ./C.(x)^2).^(3/2)-2*v*μ/g*(1-x.^2/R^2)*ξ*C.(x)*(1-v.^2 ./C.(x)^2).^(1/2)
Eki(x,v) = 4/3 *μ/g*(1-x.^2/R^2)*C.(x)*(1-v.^2 ./C.(x)^2).^(3/2)
Ep(x,v) = -2*v*μ/g*(1-x.^2/R^2)*ξ .*C.(x)*(1-v.^2 ./C.(x)^2).^(1/2)