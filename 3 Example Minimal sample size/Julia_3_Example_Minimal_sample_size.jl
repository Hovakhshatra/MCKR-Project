# The Monte Carlo integrations mentioned in the subsection "Minimal sample size" of the paper.

import Distributions: Uniform # We need the "Uniform" function from the "Distribution" package.

function samplor(a, b) # To generate a uniform random number between a and b.
    return rand(Uniform(a, b))
end

function chi(a, b, x) # The indicator function on the interval [a,b].
    if x < a
        return 0
    elseif x > b
        return 0
    else
        return 1
    end
end

@fastmath function example1(B,NN)
    # Simple Monte Carlo integration when U(B) is chosen as the probability distribution for x and the sample size is NN.
    # B is a list of two numbers.
    # NN is an integer.
    # The output is a list of three numbers; I-hat, e-hat and the computation time in seconds.
    S=0
    Ians=0
    st1=time_ns()
    COEFF=B[2]-B[1]
    x=samplor(B[1],B[2])
    Ians+=chi(0,2,x)
    Ians=COEFF*Ians
    @inbounds for i=2:NN
        x=samplor(B[1],B[2])
        delto=COEFF*chi(0,2,x)-Ians
        Ians+=delto/i
        S+=delto^2*(i-1)/i
    end
    Estandardo=S/NN
    Estandardo=Estandardo/(NN-1)
    Estandardo=sqrt(Estandardo)
    st2=time_ns()
    return Ians,Estandardo,(st2-st1)/10^9
end

# Now asking Julia to compute the Monte Karlo integration for NN=10^n for n=1,...,8, first with B=[0,1000], then for B=[0,2].
print("\nexample 1 with B=[0,10000]\n")
for n=1:8
    print("N=10^$n\n")
    X1=example1([0.0,10000.0],10^n)
    print("With simple Monte Carlo: Ians = $(X1[1]), Standard error = $(X1[2]), Computation time = $(X1[3]).\n")
end
print("\nexample 1 with B=[0,2]\n")
for n=1:8
    print("N=10^$n\n")
    X1=example1([0.0,2.0],10^n)
    print("With simple Monte Carlo: Ians = $(X1[1]), Standard error = $(X1[2]), Computation time = $(X1[3]).\n")
end

# When we compiled it, the output was as in below.
#=
example 1 with B=[0,10000]
N=10^1
With simple Monte Carlo: Ians = 0.0, Standard error = 0.0, Computation time = 2.901e-6.
N=10^2
With simple Monte Carlo: Ians = 0.0, Standard error = 0.0, Computation time = 4.499e-6.
N=10^3
With simple Monte Carlo: Ians = 0.0, Standard error = 0.0, Computation time = 2.61e-5.
N=10^4
With simple Monte Carlo: Ians = 0.0, Standard error = 0.0, Computation time = 0.000255299.
N=10^5
With simple Monte Carlo: Ians = 1.5000000000000009, Standard error = 0.3872712225172398, Computation time = 0.0023557.
N=10^6
With simple Monte Carlo: Ians = 1.9900000000000573, Standard error = 0.1410533934227002, Computation time = 0.0285829.
N=10^7
With simple Monte Carlo: Ians = 1.9440000000001039, Standard error = 0.04408653173886777, Computation time = 0.2324725.
N=10^8
With simple Monte Carlo: Ians = 1.9952000000000127, Standard error = 0.014123745741299128, Computation time = 2.2977433.

example 1 with B=[0,2]
N=10^1
With simple Monte Carlo: Ians = 2.0, Standard error = 0.0, Computation time = 1.4e-6.
N=10^2
With simple Monte Carlo: Ians = 2.0, Standard error = 0.0, Computation time = 3.0e-6.
N=10^3
With simple Monte Carlo: Ians = 2.0, Standard error = 0.0, Computation time = 2.67e-5.
N=10^4
With simple Monte Carlo: Ians = 2.0, Standard error = 0.0, Computation time = 0.0002344.
N=10^5
With simple Monte Carlo: Ians = 2.0, Standard error = 0.0, Computation time = 0.002338101.
N=10^6
With simple Monte Carlo: Ians = 2.0, Standard error = 0.0, Computation time = 0.0231234.
N=10^7
With simple Monte Carlo: Ians = 2.0, Standard error = 0.0, Computation time = 0.2326413.
N=10^8
With simple Monte Carlo: Ians = 2.0, Standard error = 0.0, Computation time = 2.310959501.
=#

# The end of the file.
