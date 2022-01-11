function trapdynamics2(sols, sim, t)
    ts = t
    Ept = zero(t)
    #Ept_e = zero(t)
    Ept_s = zero(t)
    Ept_b = zero(t)
    xt = zero(t)
    xnt = zero(t)
    #Eit =  zeros(length(ts))
    #Ns =  zeros(length(ts))
    for i in 1:length(t)
        ψ = xspace(sols[i], sim)
        psi = XField(ψ, X, K, K2)
        mask = g * abs2.(ψi) .> 0.1 * μ
        v = velocity2(psi)
        vm = v[mask]
        xm = x[mask]
        xnt[i] = xm[findmax(abs.(vm))[2]]
        if i == 1
            xt[i] = xnt[1]
        else
            j = findall(a -> a == xt[i-1], x)[1]
            mask2 = zero.(ψ)
            mask2[j-2:j+2] .= 1
            xt[i] = xm[findmax(abs.(vm .* mask2[mask]))[2]]
        end
        rhoi = abs2.(ψg)
        rho = abs2.(ψ)
        #psi = XField(ψ, X, K, K2)
        Ep = 0.5 * (x .^ 2) .* (rho - rhoi)
        Ept[i] = sum(Ep) * dx
       #  X_e = sum(x .* (rho - rhoi)) .* dx
        #norm = sum(rho) .* dx
        #Ept_e[i] = 0.5 * X_e .^ 2 / norm
        mask_s= abs.(x .-xt[i] ) .< 1
        mask_b= abs.(x .-xt[i] ) .>1
        Ept_s[i] = sum(Ep.*mask_s) * dx
        Ept_b[i] = sum(Ep.*mask_b) * dx
    end
    return Ept, Ept_s,Ept_b,xt
end


function PE(psi::XField{1})
    @unpack psiX, K = psi
    #kx = K[1]
    ψ = psiX
    rho = abs2.(ψ)
    rhoi = abs2.(ψg)
    Ep = 0.5 * (x .^ 2) .* (rho - rhoi)
    return Ep
end


function Kone_mode(sols, sim, t)
    ts = t
    Ept = zero.(ts)
    Ept_k = zero.(ts)
    #Eit =  zeros(length(ts))
    #Ns =  zeros(length(ts))
    for i in 1:length(t)
        ψ = xspace(sols[i], sim)
        rho = abs2.(ψ)
        rhoi = abs2.(ψg)
        X_e = (sum(x .* (rho)) .* dx) / (sum(rho) .* dx)
        ψk = ψ0.(x .- X_e, μ, g)
        rhok = abs2.(ψk)
        #psi = XField(ψ, X, K, K2)
        Ep = 0.5 * (x .^ 2) .* (rho - rhoi)
        Ept[i] = sum(Ep) * dx
        Ep = 0.5 * (x .^ 2) .* (rho - rhok)
        Ept_k[i] = sum(Ep) * dx

        #Ept_k[i] = 0.5 * X_e .^ 2 / norm
    end
    return Ept, Ept_k
end
function position(sols, sim, t)
    ts = t
    xnt = zero(t)
    xt2 = zero(t)

    X_e = zero(t)
    Ns = zeros(length(ts))
    for i in 1:length(t)
        ψ = xspace(sols[i], sim)
        psi = XField(ψ, X, K, K2)
        rho = abs2.(ψ)
        mask = g * abs2.(ψi) .> 0.1 * μ
        v = velocity2(psi)
        vm = v[mask]
        xm = x[mask]
        xnt[i] = xm[findmax(abs.(vm))[2]]
        if i == 1
            xt2[i] = xnt[1]
        else
            j = findall(a -> a == xt2[i-1], x)[1]
            mask2 = zero.(ψ)
            mask2[j-2:j+2] .= 1
            xt2[i] = xm[findmax(abs.(vm .* mask2[mask]))[2]]
        end
        X_e[i] = (sum(x .* (rho)) .* dx) / (sum(rho) .* dx)
    end
    return xt2, X_e
end

function trapdynamics(sols, sim, t)
    ts = t
    Ept = zeros(length(ts))
    Ept_e = zeros(length(ts))
    #Eit =  zeros(length(ts))
    #Ns =  zeros(length(ts))
    for i in 1:length(t)
        ψ = xspace(sols[i], sim)
        rhoi = abs2.(ψg)
        rho = abs2.(ψ)
        #psi = XField(ψ, X, K, K2)
        Ep = 0.5 * (x .^ 2) .* (rho - rhoi)
        Ept[i] = sum(Ep) * dx
        X_e = sum(x .* (rho - rhoi)) .* dx
        norm = sum(rho) .* dx
        Ept_e[i] = 0.5 * X_e .^ 2 / norm
    end
    return Ept, Ept_e
end


function Energy2(psi::XField{1})
    @unpack psiX, K = psi
    kx = K[1]
    ψ = psiX
    rho = abs2.(ψ)
    #rhoi = abs2.(ψg)
    #psig=XField(ψg,X,K,K2)
    ψx = gradient(psi::XField{1})
    #ψgx = gradient(psig::XField{1})
    Ek = 0.5 * abs2.(ψx) #- 0.5 * abs2.(ψgx) 
    Ep = 0.5 * (x .^ 2) .* (rho)
    Ei = 0.5 * g * rho .^ 2
    return Ek, Ep, Ei
end

function solitondynamics(sols, sim, t)
    ts = t
    #xft = zero(t)#(nearest grid point)
    xnt = zero(t)
    xt2 = zero(t)
    Ekit = zeros(length(ts))
    Ept = zeros(length(ts))
    #Eit =  zeros(length(ts))
    Ns = zeros(length(ts))
    for i in 1:length(t)
        ψ = xspace(sols[i], sim)
        psi = XField(ψ, X, K, K2)
        mask = g * abs2.(ψi) .> 0.1 * μ
        v = velocity2(psi)
        vm = v[mask]
        xm = x[mask]
        xnt[i] = xm[findmax(abs.(vm))[2]]
        if i == 1
            xt2[i] = xnt[1]
        else
            j = findall(a -> a == xt2[i-1], x)[1]
            mask2 = zero.(ψ)
            mask2[j-2:j+2] .= 1
            xt2[i] = xm[findmax(abs.(vm .* mask2[mask]))[2]]
        end

        #xnt2[i] = xm[findmin(abs2.(ψm))[2]]
        Ekit[i] = sum(Energy(psi)[1]) * dx
        Ept[i] = sum(Energy(psi)[2]) * dx
        #Eit[i] = sum( Energy(psi)[3])*dx
        Ns[i] = sum(abs2.(ψ) - abs2.(ψg)) * dx
    end
    return xt2, Ekit, Ept, Ns
end

