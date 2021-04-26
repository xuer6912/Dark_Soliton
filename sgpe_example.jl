using Plots, Interact, LaTeXStrings
gr(grid=false,size=(400,200),lw=1,c=:red)
linspace(a,b,n) = LinRange(a,b,n) |> collect

βa(x)=x*(1-tan(x)/x)
ϵ = 0.01
function boundStates(x,ϵ)
    betaa=βa.(x)
    x1=minimum(x)
    x2=maximum(x)
    dx=x[2]-x[1]
    n1 = 0
    n2 = ceil(0.5*(2*x2/π-1))
    n = collect(n1:n2)
    xbs = (1 .+2*n)π/2
    for j in 1:length(xbs)
        i = findall(abs.(x .-xbs[j]).<ϵ)
        betaa[i] .= NaN
    end
    return betaa,xbs 
end

βb = linspace(0,10,500)
bs,xbs = boundStates(βb,ϵ)

p1 = plot(βb,bs,c=:blue,lw=2,legend=false,left_margin = -2mm,bottom_margin=-2mm)
ylims!(-20,20)
xlims!(0,10)
xlabel!(L"βb")
ylabel!(L"\beta a");

ybs = linspace(-20,20,length(βb))
plot!(xbs[1]*one.(βb),ybs,ls=:dash)
plot!(xbs[2]*one.(βb),ybs,ls=:dash)
plot!(xbs[3]*one.(βb),ybs,ls=:dash)

annotate!([(1.8, 0, text(L"$\frac{\pi}{2}$",8,:black,:center))])
annotate!([(5.0, 0, text(L"$\frac{3\pi}{2}$",8,:black,:center))])
annotate!([(8.15, 0, text(L"$\frac{5\pi}{2}$",8,:black,:center))])


rm = 8
r = linspace(0,rm,500);

function kcotd(k,b,β)
    q=sqrt(k^2+β^2)
    return (q*cot(q*b)+k*tan(k*b))/(1-(q/k)*cot(q*b)*tan(k*b))
end

function u1(r,β,b,k)
    q=sqrt(k^2+β^2)
    return sin(q*r)
end

function u2(r,β,b,k,δ₀)
    q=sqrt(k^2+β^2)
    return sin(q*b)/sin(k*b+δ₀)*sin.(k*r+δ₀)
end

d0(k,b,β)=acot(k^(-1)*kcotd(k,b,β))
δ0(k,b,β)=acot(-β/βa(β*b)/k)

b=1.
β=4.8
k=.1
d0(k,b,β),δ0.(k,b,β)

@manipulate for k=0.05:0.02:3.0, β=.5:.1:5.5*pi/2 
    #draw potential
    p2 = plot(size=(500,200),legend=false)
    i = findall(r.<1)
    j = findall(r.>1)
    rb = r[i]
    rc = r[j]
    l1 = -β^2*one.(rb)
    l2 = linspace(-β^2,0,100)
    l3 = zero.(rc)
    l4 = zero.(rb)
    l5 = -β^2*one.(rc)
    plot!(rb,l1)
    plot!(one.(l2),l2)
    plot!(rc,l3)
    amp = .7*β*sin.(k*b .+d0(k,b,β))/sin.(sqrt.(k^2 .+β^2)*b)
    #amp = 1
    #draw wavefunction
    y1=k.^2 .+amp^2*abs2.(u1.(rb,β,b,k)) 
    y2=k.^2 .+amp^2*abs2.(u2.(rc,β,b,k,d0.(k,b,β))) 
    y3=k.^2 .+amp*u1.(rb,β,b,k)
    y4=k.^2 .+amp*u2.(rc,β,b,k,d0.(k,b,β))
    plot!(rb,y1,c=:grey)
    plot!(rc,y2,c=:grey,alpha=1)
    plot!(rb,y3,c=:blue)
    plot!(rc,y4,c=:blue)

    ymax2=maximum(y2)
    ymax4=maximum(y4)
    yminU=-β^2
    ymin = yminU
    ymax=maximum([ymax2,ymax4])
    ylims!(1.1*ymin,1.5*ymax2)
    xlims!(0,rm)
    xlabel!(L"r/b")
    ylabel!(L"Energy");

  #draw bound state energies
    nmax = ceil(0.5*(2*β/π-1))
    #show(nmax)
    n = 0:nmax
    eb = (2*n .+1)π/2
    bline=one.(rb)
    for j=1:length(eb)
        plot!(rb,l1 .+eb[j].^2,c=:black,ls=:dash)
    end
   #plot bound states
    j = findall(eb.^2 .<β^2);
    ebn=eb[j];
    for i=1:length(ebn)
        z1 = l1 .+ebn[i].^2 .+amp^2*abs2.(u1.(rb,im*β,b,k))
        z2 = l5 .+ebn[i].^2 .+amp^2*abs2.(u2.(rc,im*β,b,k,d0.(k,b,im*β)))
        #plot(rb,z1)
    end
    plot(p2)
end