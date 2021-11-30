module noisyOracle
    using LinearAlgebra
    using Statistics
    import Distributions.Normal
    export construct_oracle_f, construct_oracle_s, ⊗

    ⊗ = kron
    single_paulis = [
        [1.0 0; 0 1.0],
        [0 1.0; 1.0 0],
        [1 0; 0 -1],
        [0 -1.0im; 1.0im 0]
    ]; # I, X, Z, Y

    function construct_pauli(x, n)
        x = x - 1
        m = 1
        for i = 1:n
            r = (x % 4) + 1
            x = x ÷ 4
            m = single_paulis[r] ⊗ m
        end
        return m
    end

    # oracle for phase 1
    # s: Dict, terms in the Hamiltonian
    # n: qubit number
    # σ: standard deviation of noise
    function construct_oracle_f(s, n, σ; m=2, deltat=0.01, timeSteps=8)

        ham = zeros(2^n, 2^n)
        for (x, val) in s
            ham += val*construct_pauli(x, n)
        end

        t = [j*deltat for j=1:timeSteps]
        ms = vcat([0, (2:m)...])
        X = hcat([t.^i for i in ms]...);
        Ut = exp.(-1im * fill(ham, timeSteps) .* t);
        UtDagger = conj.(transpose.(Ut));
        d = Normal(0.0, σ);

        function oracle_f(k)
            P = fill(construct_pauli(k, n), timeSteps)
            ft = real(tr.(P .* Ut .* P .* UtDagger)) / (2.0^n) .+ rand(d, timeSteps);
            return (X \ ft)[2]
        end

        return ham, oracle_f
    end

    # oracle for phase 2
    # ham: reuse the matrix returned by oracle_f
    function construct_oracle_s(ham, n, σ)
        eigenStates = [
                        [0.5 0.5; 0.5 0.5],
                        [1.0 0.0; 0.0 0.0],
                        [0.5 -0.5im; 0.5im 0.5]
        ]
    
        function construct_states(x)
            m = 1
            for i = 1:n
                r = (x % 4)
                x = x ÷ 4
                m = eigenStates[r] ⊗ m
            end
            return m
        end
    
        timeSteps = 8
        t = [j*0.01 for j=1:timeSteps]
        X = [ones(timeSteps) t];
        Ut = exp.(-1im * fill(ham, timeSteps) .* t);
        UtDagger = conj.(transpose.(Ut));
        d = Normal(0.0, σ);
    
        function oracle_s(β, γ)
            ρ = construct_states(γ)
            M = construct_pauli(β+1, n)
            ρs = fill(ρ, timeSteps)
            Ms = fill(M, timeSteps)
            ft = real(tr.(Ms .* Ut .* ρs .* UtDagger)) .+ rand(d, timeSteps)
            return (X \ ft)[2]
        end
    
        return oracle_s
    end
end