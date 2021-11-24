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
#         println(maps(k))
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

s = Dict{Any, Any}(5 => 1.0, 41 => 1.0, 65 => 1.0, 11 => 1.0, 2 => 1.0, 161 => 1.0, 17 => 1.0)

# Dict{Any, Any}(5 => 1.0, 41 => 1.0, 65 => 1.0, 11 => 1.0, 2 => 1.0, 161 => 1.0, 641 => 1.0, 17 => 1.0, 257 => 1.0)

modifiedParams = [(1e-4, 4e-3, 100, 2e-2), (1e-3, 4e-3, 100, 2e-2), (3e-3, 0.1, 100, 0.2), (6e-3, 0.2, 100, 0.2), (1e-2, 0.3, 100, 0.4), (0.015, 0.5, 100, 0.6)]

# n = parse(Int64, ARGS[1])
n = 4
verbose = parse(Bool, ARGS[2])
totalRounds = parse(Int64, ARGS[3])
filename = ARGS[4]


params = Dict(
    "Cutoff"=>0.995,
    "Zerosensitivity"=>2e-2,
    "noise"=>0.001,
    "largeEnoughOffset"=>100,
    "single"=>1e-3,
    "phase2Lines"=>200
);

A = params["largeEnoughOffset"]


for j in 6:6
    savingItems = []
    ave=0.0
    params["noise"] = modifiedParams[j][1]
    params["single"] = modifiedParams[j][2]
    params["phase2Lines"] = modifiedParams[j][3]
    params["Zerosensitivity"] = modifiedParams[j][4]

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
    failCount = 0
    for i=1:totalRounds
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
            println("=======================")
        
            nonZeroAlphas = filter!(x->x≠0, [a for a in keys(ps_re)])

            global s_re = pauliparametersReconstruction(nonZeroAlphas, params["phase2Lines"], sOracle)
            println("total Calls to oracle at round ",i," with qubit ",n ,": " ,totalCall+totalCall2)
            global fs = reconstructedParameters(ps_re, s_re)
        catch
            println("error in do peel")
            failCount += 1
            i = i-1
            continue
        end
        println(printboth(s, fs))
        println("distance at round ",i,": ",parameterDis(s,fs))
        println("=======================")
        ave += parameterDis(s,fs)
        push!(savingItems, (fs, totalCall, totalCall2))
    end
    ave/=totalRounds
    println(ave)
    println("fail rate ", failCount/totalRounds)
    
    open("./data/IsingRounds=500_noise="*string(params["noise"])*".json", "w") do f
        JSON.print(f, savingItems, 4)
    end
end
# println(savingItems)