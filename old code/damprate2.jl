


T = 30e-9
ϵc(μ) = 2*μ
μr = LinRange(.5, 10,100)*μ1
plot(μr/(ħ*wz), Γ_γ.(μr)*sqrt(2)/wz, xaxis=:log, yaxis=:log,label=L"\gamma:2\mu",legend=:topright)
plot!(μr/(ħ*wz) , Γ_M.(μr)*sqrt(2)/wz, xaxis=:log, yaxis=:log,label=L"M:2\mu",size = (600,300))




##
ϵc(μ) = 2.5*μ
#μr = LinRange(.5, 10,100)*μ1
plot!(μr/(ħ*wz), Γ_γ.(μr)*sqrt(2)/wz, xaxis=:log, yaxis=:log,label=L"\gamma:2.5\mu")
plot!(μr/(ħ*wz) , Γ_M.(μr)*sqrt(2)/wz, xaxis=:log, yaxis=:log,label=L"M:2.5\mu",size = (600,300))

ϵc(μ) = 3*μ
plot!(μr/(ħ*wz), Γ_γ.(μr)*sqrt(2)/wz, xaxis=:log, yaxis=:log,label=L"\gamma:3\mu")
plot!(μr/(ħ*wz) , Γ_M.(μr)*sqrt(2)/wz, xaxis=:log, yaxis=:log,label=L"M:3\mu",size = (600,300))

vline!([0.5wr/wz],label =L"50\hbar\omega_z")
xlabel!(L"\mu/\hbar\omega_z")
ylabel!(L"\Gamma/\omega_s")
savefig("dampratecut.pdf")



#################

ϵc(μ) = 3*μ

T = 12e-9
μr = LinRange(.5, 10,100)*μ1
p1 = plot(μr/(ħ*wz), Γ_γ.(μr)*sqrt(2)/wz, xaxis=:log, yaxis=:log,label=L"\gamma:12nK",legend=:topright)
plot!(μr/(ħ*wz) , Γ_M.(μr)*sqrt(2)/wz, xaxis=:log, yaxis=:log,label=L"M:12nK",size = (600,300))

T = 30e-9
plot!(μr/(ħ*wz), Γ_γ.(μr)*sqrt(2)/wz, xaxis=:log, yaxis=:log,label=L"\gamma:30nK")
plot!(μr/(ħ*wz) , Γ_M.(μr)*sqrt(2)/wz, xaxis=:log, yaxis=:log,label=L"M:30nK",size = (600,300))

T = 48e-9
plot!(μr/(ħ*wz), Γ_γ.(μr)*sqrt(2)/wz, xaxis=:log, yaxis=:log,label=L"\gamma:48nK")
plot!(μr/(ħ*wz) , Γ_M.(μr)*sqrt(2)/wz, xaxis=:log, yaxis=:log,label=L"M:48nK",size = (600,300))

vline!([0.5wr/wz],label =L"50\hbar\omega_z")
xlabel!(L"\mu/\hbar\omega_z")
ylabel!(L"\Gamma/\omega_s")


savefig("damprateT.pdf")