import Distributions.Normal

function TFIsingHamiltonian1d(n)
    s = Dict()
    d = Normal(0.0, 1.0)
    # each term with absolute value
    # at least [cutoff]
    cutoff = 0.5
    for i = 0:n-2
        j = 2 * (4^i) + 2 * (4^(i+1)) + 1 # "II...IZZII..II"
        a = rand(d)
        s[j] = a + sign(a)*cutoff
        a=rand(d)
        s[4^i + 1] = a + sign(a)*cutoff
    end
    a = rand(d)
    s[4^(n-1) + 1] =a + sign(a)*cutoff
    return s
end