# generate the result for reconstruction
# of random transverse-field Ising Hamiltonian
#
# for each qubit number n in 1:7
# generate [totalHams] TFIsing Hamiltonians (see hamsGen.jl)
# each is reconstructed [eachRounds] times
# output as json
#
# usage: julia hamLearning_RandomIsing.jl [verbose (true|false)]
#
# 


import Distributions.Normal
using Distributed
include("./noisy_oracle.jl")
include("./utils.jl")
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

# modifiedParams:
# (singleSensitivity, linesInPhase2, zeroSensitivity, b, n, m, deltat, timeSteps)
modifiedParams = [
    (5e-3, 100, 3e-2, 4, 4, 2, 0.01, 10),
    (50e-3, 100, 6e-2, 4, 4, 4, 0.003, 24),
    (80e-3, 100, 10e-2, 4, 4, 6, 0.0008, 100),]

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

for j in 3:3

    savingItems = [] # will be the json output

    global n = modifiedParams[j][5]
    m = modifiedParams[j][6]
    params["single"] = modifiedParams[j][1]
    params["phase2Lines"] = modifiedParams[j][2]
    params["Zerosensitivity"] = modifiedParams[j][3]
    b = modifiedParams[j][4]
    failCount = 0
    ave=0.0
    for hamsCount in 1:totalHams
        s = TFIsingHamiltonian1d(n)
        push!(savingItems, s);
        ham0, oracle_f = construct_oracle_f(s, n, params["noise"], m=m, deltat=modifiedParams[j][7], timeSteps=modifiedParams[j][8])

        # real oracle being actually used
        # total calls is counted
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
            global totalCall = 0
            global totalCall2 = 0
            # try
                global ps_re = mainProc.dopeel(b, fOracle, params, 3, 10, verbose)
                if verbose
                    println("=======================")
                end
            
                nonZeroAlphas = filter!(x->x≠0, [a for a in keys(ps_re)])

                global s_re = pauliparametersReconstruction(nonZeroAlphas, params["phase2Lines"], sOracle)
            # catch
            #     println("error in do peel")
            #     failCount += 1
            #     i = i-1
            #     continue
            # end
            if verbose
                println("total Calls to oracle at round ",i," with qubit ",n ,": " ,totalCall+totalCall2)
            end
            fs = reconstructedParameters(ps_re, s_re)
            if verbose
                printboth(s, fs)
                println("distance at round ",(hamsCount-1)*eachRounds + i,": ",parameterDis(s,fs))
                println("=======================")
            end
            ave += parameterDis(s,fs)
            push!(savingItems, (fs, totalCall, totalCall2))
        end
    end
    ave/=(eachRounds * totalHams)
    println("average error: ", ave)
    println("fail rate ", failCount/(eachRounds * totalHams))
    
    open("./data/RandomIsingRounds=500_m="*string(m)*".json", "w") do f
        JSON.print(f, savingItems, 4)
    end
end