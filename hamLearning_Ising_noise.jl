import Distributions.Normal
using Distributed
include("./utils.jl")
include("./noisy_oracle.jl")
include("./hamsGen.jl")
using LinearAlgebra, Statistics, .hamLearning_utils, .noisyOracle
using FileIO, JSON


function printboth(s, fs)
    maps = (k) -> (haskey(s, k) ? s[k] : nothing)
    mapfs = (k) -> (haskey(fs, k) ? fs[k] : nothing)
    ks = union(keys(s), keys(fs))
    for k in ks
        println(k, "=>", "(", maps(k), ", ", mapfs(k), ")")
    end
end

function parameterDis(s, s_re)
    maps = (k) -> (haskey(s, k) ? s[k] : 0.0)
    mapsre = (k) -> (haskey(s_re, k) ? s_re[k] : 0.0)
    ans = 0.0
    for key in union(keys(s), keys(s_re))
        ans+=abs(maps(key) - mapsre(key))
    end
    idea=0.0
    for (key,val) in s
       idea+= abs(val)
    end
    return (ans)/(idea)
end

# (noise, single, phase2Lines, Zerosensitivity)

modifiedParams = [
    (1e-4, 1e-3, 100, 2e-2),
    (1e-3, 4e-3, 100, 4e-2),
    (3e-3, 0.5, 100, 10e-2),
    (6e-3, 0.5, 100, 0.2),
    (1e-2, 0.5, 100, 0.6),
    (0.015, 0.5, 100, 0.6)
]

n = 4
verbose = parse(Bool, ARGS[1])


params = Dict(
    "Cutoff"=>0.995,
    "Zerosensitivity"=>2e-2,
    "noise"=>0.001,
    "largeEnoughOffset"=>100,
    "single"=>1e-3,
    "phase2Lines"=>200
);

A = params["largeEnoughOffset"]

totalHams = 50
eachRounds = 10

for j in 2:2
    savingItems = []
    ave=0.0
    params["noise"] = modifiedParams[j][1]
    params["single"] = modifiedParams[j][2]
    params["phase2Lines"] = modifiedParams[j][3]
    params["Zerosensitivity"] = modifiedParams[j][4]
    failCount = 0
    for hamsCount in 1:totalHams
        s = TFIsingHamiltonian1d(n)
        push!(savingItems, s);
        ham0, oracle_f = construct_oracle_f(s, n, params["noise"])
        # real oracle being actually used
        function fOracle(k)
            if haskey(F,k)
                return F[k]
            end
            global totalCall = totalCall + 1
            F[k] = oracle_f(k+1) + A
            return F[k]
        end
        oracle_s = construct_oracle_s(ham0, n, params["noise"])

        function sOracle(β, γ)
            global totalCall2 = totalCall2 + 1
            return oracle_s(β, γ)
        end

        mainProc = hamLearningProc(n)
        for i=1:eachRounds
            global F = Dict()
            # t = [j*0.01 for j=1:timeSteps]
            # global X = [ones(timeSteps) t.^2];
            # global Ut = exp.(-1im * fill(ham0, timeSteps) .* t);
            # global UtDagger = conj.(transpose.(Ut));
            # totalCall = 0
            estimateCalls = 4^n
            global totalCall = 0
            global totalCall2 = 0
            # progbar = Progress(estimateCalls, desc="sample")
            try
                global ps_re = mainProc.dopeel(4, fOracle, params, 3, 10, verbose)
                # ps_re = doPeel(4, fOracle2, params, 3,10)
                if verbose
                    println("=======================")
                end
            
                nonZeroAlphas = filter!(x->x≠0, [a for a in keys(ps_re)])

                global s_re = pauliparametersReconstruction(nonZeroAlphas, params["phase2Lines"], sOracle)
                if verbose
                    println("total Calls to oracle at round ", i + (totalHams - 1)*eachRounds, " with qubit ",n ,": " ,totalCall+totalCall2)
                end
                global fs = reconstructedParameters(ps_re, s_re)
            catch
                println("error in do peel")
                failCount += 1
                i = i-1
                continue
            end
            if verbose
                printboth(s, fs)
                println("distance at round ", i ,": ",parameterDis(s,fs))
                println("=======================")
            end
            ave += parameterDis(s,fs)
            push!(savingItems, (fs, totalCall, totalCall2))
        end
    end
    ave /= (totalHams * eachRounds)
    println("average error: ", ave)
    println("fail rate ", failCount/(totalHams * eachRounds))
    
    open("./data/RandomIsingRounds=500_noise="*string(params["noise"])*".json", "w") do f
        JSON.print(f, savingItems, 4)
    end
end

# println(savingItems)