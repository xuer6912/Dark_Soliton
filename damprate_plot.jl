using Plots,LaTeXStrings
#physical constant
amu = 1.66053e-27
a0 = 5.29e-11 # Bohr radius
ħ = 1.055e-34
kB = 1.381e-23
##
#Rubidium 87
m = 87*amu 
a = 95*a0 
g = 4*π*ħ^2*a/m
##
#trap
wz = 2*pi*3
wr = 2*pi*300
az = sqrt(ħ/m/wz)
ar = sqrt(ħ/m/wr)
g1 = 2*ħ*wr*a 
##
#temperature
T = 50e-9
Tc = 60e-9
β(T) = 1/(kB*T)
ϵc(μ) = 3*μ
μ1 = ħ*wz*10 # entering Thomas-Fermi regime
λdB(T)=sqrt(2π*ħ^2/(m*kB*T))
λ=λdB(T)
function Lerch(z, j)
    Φ = 0
    for i in 1:1000
        Φ += z^(i-1)/(j+z)
    end
    return Φ
end

begin
	function gamma(μ)
	    gam = 0
	    for j in 1:1000
	        Φ = Lerch.(exp.( β(T) * μ - 2*β(T)*ϵc(μ)), j)
	        gam += exp.(β(T) * μ*(j+1) - 2 * β(T) * ϵc(μ) *j ) * Φ ^ 2
	    end
	    γ = 8*a^2 / λ^2 * gam
	    return γ
	end
	function Γ_γ(μ)
	    Γ_γ =gamma(μ)/3 *μ/ħ
	    return Γ_γ
	end  
	
	function Γ_M(μ)
	    M1 = (16* pi *a^2/(exp(β(T)*(ϵc(μ)-μ))-1)) / sqrt(8*pi*ar^2)
	    Γ_M = M1 * (2/15) * μ.^2 /ħ/g1
	    return Γ_M
	end
end
##

μr = LinRange(.5, 10,100)*μ1
p1 = plot(μr/(ħ*wz), Γ_γ.(μr)*sqrt(2)/wz, xaxis=:log, yaxis=:log,label=L"\Gamma_\gamma",legend=:topleft)
plot!(μr/(ħ*wz) , Γ_M.(μr)*sqrt(2)/wz, xaxis=:log, yaxis=:log,label=L"\Gamma_M")

vline!([0.5wr/wz])
xlabel!(L"\mu/\hbar\omega_z")
ylabel!(L"\Gamma\sqrt{2}/\omega_z")
##


	##
	# savefig(L"\Gamma")
##
N = 1.21e5
##
g2d(az) = g/sqrt(2*pi*az^2)
μtf(wr,az,a,N) = ħ*wr*sqrt(sqrt(8/pi)*a/az*N)  
Rtf(μ,w) = sqrt(2*μ/m/w^2)
##
μ  = μtf(wr,az,a,N)
n0 = μ/g 
c = sqrt(μ/m)
ξ = ħ/m/c
R = Rtf(μ,wr)
μ/(ħ*wr),μ/(ħ*wz),R/ar,ξ/ar
ξ
c
n0 






T = 30e-9
ϵc(μ) = 2*μ
μr = LinRange(.5, 10,100)*μ1
p1 = plot(μr/(ħ*wz), Γ_γ.(μr)*sqrt(2)/wz, xaxis=:log, yaxis=:log,label=L"\Gamma_\gamma",legend=:topleft)
plot!(μr/(ħ*wz) , Γ_M.(μr)*sqrt(2)/wz, xaxis=:log, yaxis=:log,:dash,label=L"\Gamma_M",size = (600,300),linesty)

vline!([0.5wr/wz],label =L"50\hbar\omega_z")
xlabel!(L"\mu/\hbar\omega_z")
ylabel!(L"\Gamma\sqrt{2}/\omega_z")


##
ϵc(μ) = 2.5*μ
μr = LinRange(.5, 10,100)*μ1
p2 = plot(μr/(ħ*wz), Γ_γ.(μr)*sqrt(2)/wz, xaxis=:log, yaxis=:log,label=L"\Gamma_\gamma",legend=:topleft)
plot!(μr/(ħ*wz) , Γ_M.(μr)*sqrt(2)/wz, xaxis=:log, yaxis=:log,label=L"\Gamma_M",size = (600,300))

vline!([0.5wr/wz],label =L"50\hbar\omega_z")
xlabel!(L"\mu/\hbar\omega_z")
ylabel!(L"\Gamma\sqrt{2}/\omega_z")

ϵc(μ) = 3*μ
μr = LinRange(.5, 10,100)*μ1
p3 = plot(μr/(ħ*wz), Γ_γ.(μr)*sqrt(2)/wz, xaxis=:log, yaxis=:log,label=L"\Gamma_\gamma",legend=:topleft)
plot!(μr/(ħ*wz) , Γ_M.(μr)*sqrt(2)/wz, xaxis=:log, yaxis=:log,label=L"\Gamma_M",size = (600,300))

vline!([0.5wr/wz],label =L"50\hbar\omega_z")
xlabel!(L"\mu/\hbar\omega_z")
ylabel!(L"\Gamma\sqrt{2}/\omega_z")

plot(p1,p2,p3,layout=(3,1),size=(600,800))
savefig("dampratecut")



#################

ϵc(μ) = 3*μ

T = 12e-9
μr = LinRange(.5, 10,100)*μ1
p1 = plot(μr/(ħ*wz), Γ_γ.(μr)*sqrt(2)/wz, xaxis=:log, yaxis=:log,label=L"\Gamma_\gamma",legend=:topleft)
plot!(μr/(ħ*wz) , Γ_M.(μr)*sqrt(2)/wz, xaxis=:log, yaxis=:log,label=L"\Gamma_M",size = (600,300))

vline!([0.5wr/wz],label =L"50\hbar\omega_z")
xlabel!(L"\mu/\hbar\omega_z")
ylabel!(L"\Gamma\sqrt{2}/\omega_z")

T = 30e-9
μr = LinRange(.5, 10,100)*μ1
p2 = plot(μr/(ħ*wz), Γ_γ.(μr)*sqrt(2)/wz, xaxis=:log, yaxis=:log,label=L"\Gamma_\gamma",legend=:topleft)
plot!(μr/(ħ*wz) , Γ_M.(μr)*sqrt(2)/wz, xaxis=:log, yaxis=:log,label=L"\Gamma_M",size = (600,300))

vline!([0.5wr/wz],label =L"50\hbar\omega_z")
xlabel!(L"\mu/\hbar\omega_z")
ylabel!(L"\Gamma\sqrt{2}/\omega_z")

T = 48e-9
μr = LinRange(.5, 10,100)*μ1
p3 = plot(μr/(ħ*wz), Γ_γ.(μr)*sqrt(2)/wz, xaxis=:log, yaxis=:log,label=L"\Gamma_\gamma",legend=:topleft)
plot!(μr/(ħ*wz) , Γ_M.(μr)*sqrt(2)/wz, xaxis=:log, yaxis=:log,label=L"\Gamma_M",size = (600,300))

vline!([0.5wr/wz],label =L"50\hbar\omega_z")
xlabel!(L"\mu/\hbar\omega_z")
ylabel!(L"\Gamma\sqrt{2}/\omega_z")

plot(p1,p2,p3,layout=(3,1),size=(600,800))
savefig("damprateT")