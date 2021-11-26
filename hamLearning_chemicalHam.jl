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


s = Dict(
    # 1  =>  -2.7104083700726536 , # const in the Hamiltonian has no effect
    3073  =>  -1.745198975210092 ,
    769  =>  -1.2726305026333158 ,
    3841  =>  0.5647364031520397 ,
    193  =>  -0.12777587071339788 ,
    3265  =>  0.484399127629246 ,
    961  =>  0.41584938456426535 ,
    49  =>  -1.745198975210092 ,
    3121  =>  0.7823637778985213 ,
    817  =>  0.6659936243343323 ,
    241  =>  0.5389855881574944 ,
    13  =>  -1.272630502633316 ,
    3085  =>  0.6659936243343323 ,
    781  =>  0.7400280578673019 ,
    205  =>  0.5684250585068191 ,
    61  =>  0.5647364031520397 ,
    4  =>  -0.12777587071339735 ,
    3076  =>  0.5389855881574944 ,
    772  =>  0.5684250585068191 ,
    196  =>  1.0024488587320461 ,
    52  =>  0.484399127629246 ,
    16  =>  0.41584938456426535 ,
    321  =>  0.3806607628552998 ,
    3393  =>  -0.029793806500249334 ,
    641  =>  0.3806607628552998 ,
    3713  =>  -0.029793806500249334 ,
    369  =>  -0.09039934735000948 ,
    689  =>  -0.09039934735000948 ,
    333  =>  -0.036634285428473336 ,
    653  =>  -0.036634285428473336 ,
    324  =>  0.22383338673048747 ,
    644  =>  0.22383338673048747 ,
    1301  =>  0.10125722118229263 ,
    2581  =>  0.10125722118229263 ,
    1321  =>  0.10125722118229263 ,
    2601  =>  0.10125722118229263 ,
    1877  =>  0.06060554084976015 ,
    2965  =>  0.06060554084976015 ,
    1897  =>  0.06060554084976015 ,
    2985  =>  0.06060554084976015 ,
    1310  =>  0.06060554084976015 ,
    2590  =>  0.06060554084976015 ,
    1327  =>  0.06060554084976015 ,
    2607  =>  0.06060554084976015 ,
    1886  =>  0.05458646052824838 ,
    2974  =>  0.05458646052824838 ,
    1903  =>  0.05458646052824838 ,
    2991  =>  0.05458646052824838 ,
    6  =>  0.38066076285529965 ,
    3078  =>  -0.09039934735000948 ,
    774  =>  -0.036634285428473336 ,
    198  =>  0.22383338673048747 ,
    54  =>  -0.029793806500249334 ,
    11  =>  0.38066076285529965 ,
    3083  =>  -0.09039934735000948 ,
    779  =>  -0.036634285428473336 ,
    203  =>  0.22383338673048747 ,
    59  =>  -0.029793806500249334 ,
    326  =>  0.15257567394255375 ,
    646  =>  0.15257567394255375 ,
    331  =>  0.15257567394255375 ,
    163  =>  0.15257567394255375 ,
)

n = 6
verbose = parse(Bool, ARGS[1])


params = Dict(
    "Cutoff"=>0.995,
    "Zerosensitivity"=>4e-3,
    "noise"=>0.001,
    "largeEnoughOffset"=>100,
    "single"=>4e-2,
    "phase2Lines"=>200
);

A = params["largeEnoughOffset"]

# totalHams = 50
eachRounds = 200

for b in 4:8
    savingItems = []
    ave=0.0
    failCount = 0
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
            global ps_re = mainProc.dopeel(b, fOracle, params, 3, 10, verbose)
            # ps_re = doPeel(4, fOracle2, params, 3,10)
            if verbose
                println("=======================")
            end
        
            nonZeroAlphas = filter!(x->x≠0, [a for a in keys(ps_re)])

            global s_re = pauliparametersReconstruction(nonZeroAlphas, params["phase2Lines"], sOracle)
            if verbose
                println("total Calls to oracle at round ", i, " with qubit ",n ,": " ,totalCall+totalCall2)
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
    ave /= (eachRounds)
    println("average error: ", ave)
    println("fail rate ", failCount/(eachRounds))
    
    open("./data/LiH4Rounds=500_b="*string(b)*".json", "w") do f
        JSON.print(f, savingItems, 4)
    end
end

# println(savingItems)