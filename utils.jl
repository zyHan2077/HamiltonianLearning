module hamLearning_utils
    using Statistics
    export hamLearningProc, pauliparametersReconstruction, reconstructedParameters

    ## ======================begin of phase 1==============================
    ⊕ = (x,y)->mod.(x+y,2)

    function generateJn()
        m = zeros(Int, 2*n, 2*n)
        for i = 1:2*n-1
            m[i+1, i] = mod(i, 2)
            m[i, i+1] = mod(i, 2)
        end
        return m
    end

    function int2bin(num, len)
        return map(x->parse(Int,x),collect(reverse(string(num,base=2,pad=len))))
    end
    
    function bin2int(s)
        return parse(Int, reverse(join(s)); base=2)
    end
    
    function normalBinDot(x, y)
        return (-1)^(mod((x'*(y)), 2))
    end
    
    function innerBinDot(x, y)
        return (-1)^(mod((x'*(Jn*y)), 2))
    end
    
    # j from 0 to 2^b-1
    function bPointSubBin(b, M, j, d, fOracle)
        j = int2bin(j, b)
        ans = 0.0
        Mp = Jn * M
        for l in 0:2^b-1
            binl = int2bin(l, b)
            k = (Mp * binl) ⊕ d
            ans += normalBinDot(j, binl) * fOracle(bin2int(k))
        end
        return ans / (2^b)
    end
    
    function linearDecoder(U)
        baseSign = signbit(U[1])
        return Jn*Int.(signbit.(U[2:2*n+1]))
    end
    
    function binDetector(b, M, U, ds, fOracle, params)
        if mean(U.^2) < params["Zerosensitivity"]
            return ("zero-ton", nothing, nothing)
        end
        m = linearDecoder(U)
        ls = [(innerBinDot(m, ds[d]))*U[d] for d in 1:length(ds)]
        pm = mean(ls)
        v = (U - [innerBinDot(m, d) for d in ds] * pm)
        if mean(v.^2) < params["single"]
            return ("single-ton", bin2int(m), pm)
        end
        return ("multi-ton", nothing, nothing)
    end

    function doPeel(b, fOracle, params, cycles, rounds, verbose=false)
        total_sum = 0.0
        ps_reconstructed = Dict()
        M=[]
        U=[]
        d=[]
        for i=1:cycles
            ds=vcat(
                [[0 for _ = 1:2*n]],
                [map(x->parse(Int,x),collect(reverse(string(i,base=2,pad=2*n)))) for i in [2^b for b=0:2*n-1]],
                [map(x->parse(Int,x),collect(reverse(string(rand(0:4^n-1),base=2,pad=2*n)))) for i=1:3*n ])
            push!(M,rand(0:1, 2*n, b))
            push!(d,ds)
            UC=[]
            for j=1:2^b
                f = x -> bPointSubBin(b, M[i], j-1, x, fOracle)
                Uc = f.(ds)
                push!(UC,Uc)
            end
            push!(U,UC)
        end
        flag=false
        for k =1:rounds
            if flag==true
                if verbose
                    println("End at round ", k-1)
                end
                return ps_reconstructed
            end
            flag=true
            for i = 1:cycles
                for j = 1:2^b
                    tonType, m, pm = binDetector(b, M[i], U[i][j], d[i], fOracle, params)
                    if tonType == "single-ton"
                        if !haskey(ps_reconstructed, m)
                            if verbose
                                println("m=", m, " p_m=", pm, " at cycle ", i, " at round ",k)
                            end
                            ps_reconstructed[m] = pm
                            flag=false
                            for t=1:cycles
                                if t!=i
                                    mbin = int2bin(m, 2*n)
                                    j_occur = bin2int(mod.(M[t]'*mbin, 2))
                                    g = x -> innerBinDot(x,mbin)
                                    U[t][j_occur+1] = U[t][j_occur+1] - pm * g.(d[t])
                                end
                            end
                        end
                    end
                end
            end
        end    
        println("reconstruction failed")
        return ps_reconstructed
    end

    ## =============================end of phase 1 ===================================================

    ## ======================== utils for phase 2 ======================================================


    function int2quad(a)
        return map(x->parse(Int,x),collect(reverse(string(a,base=4,pad=n))))
    end
    
    function determineAllCoefficients(αs, β, γ)
        α0 = bin2int(int2bin(β, 2*n) ⊕ int2bin(γ, 2*n))
        betaLs = int2quad(β)
        alpha0Ls = int2quad(α0)
        len = length(αs)
        row = zeros(len)
        for i in 1:len
            α = αs[i]
            alphaLs = int2quad(α)
            validBits = [(alphaLs[k] == alpha0Ls[k]) || (alphaLs[k] == betaLs[k]) for k in 1:n]
            if !(false in validBits)
                if innerBinDot(int2bin(α, 2*n), int2bin(β, 2*n)) == -1
                    row[i] = -basisProjection(α, β)
                end
            end
        end
        return row
    end

    function basisProjection(a, b)
        ans = 1.0 + 0.0im
        for i = 1:n
            r1 = (a % 4)
            r2 = (b % 4)
            a = a ÷ 4
            b = b ÷ 4
            if r1 !=0 && r2 != 0
                if (r1, r2) in [(1, 3), (3, 2), (2, 1)]
                    ans *= 1im
                else
                    ans *= -1im
                end
            end
        end
        return imag(ans) * 2.0
    end

    function pauliparametersReconstruction(αs, cycles, oracle)
        s = Dict()
        As = []
        b = Vector{Float64}()
        for α in αs
            for i in 1:cycles
                γ = parse(Int64, join(rand(1:3, n)), base=4)
                β = bin2int(int2bin(α, 2*n) ⊕ int2bin(γ, 2*n))
                #  assure that each α can be included
                if innerBinDot(int2bin(α, 2*n), int2bin(β, 2*n)) == -1
                    push!(As, determineAllCoefficients(αs, β, γ))
                    push!(b, oracle(β, γ))
                end
            end
        end
        A = As[1]'
        for j in 2:length(As)
            A = vcat(A, As[j]')
        end
        try
            vals = A \ b
            result = Dict()
            for k in 1:length(αs)
                result[αs[k]+1] = vals[k]
            end
            return result
        catch e
            println("bad condition, retry ...", e)
            return "failed"
        end
    end

    function reconstructedParameters(ps, alphas)
        result = Dict()
        for (key, val) in alphas
            result[key] = sign(val) * sqrt(ps[key-1])
        end
        return result
    end

    ## =============== initialization and phase 1 =======================
    struct hamLearningProc
        # params
        dopeel
        # pauliReconstruction
        function hamLearningProc(n_)
            global n = n_
            global Jn = generateJn()
            dopeel = (b, fOracle, params, cycles, rounds, verbose) -> doPeel(b, fOracle, params, cycles, rounds, verbose)
            new(dopeel)
        end
    end

end