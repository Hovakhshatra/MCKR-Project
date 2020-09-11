# Computations for the double HK network in the paper. In this file we keep all 14 parameters free.
# For the case that all parameters other than k12 and k14 are fixed, see the Julia file Julia_12_Example_Double_HK_2.jl

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

@fastmath function g(t,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k13,k14)
    return -(t^5*k1*k4*k5*k6*k7*k10*k11 + t^5*k1*k4*k5*k6*k8*k10*k11 + t^5*k2*k4*k5*k6*k7*k10*k11 + t^5*k2*k4*k5*k6*k8*k10*k11 - t^4*k1*k4*k5*k6*k7*k10*k11*k13 - t^4*k1*k4*k5*k6*k8*k10*k11*k13 + t^4*k1*k4*k5*k7*k8*k10*k11*k14 - t^4*k2*k4*k5*k6*k7*k10*k11*k13 - t^4*k2*k4*k5*k6*k8*k10*k11*k13 + t^4*k2*k4*k5*k7*k8*k10*k11*k14 + t^4*k1*k2*k5*k6*k7*k10*k11 + t^4*k1*k2*k5*k6*k8*k10*k11 + t^4*k1*k3*k5*k6*k7*k10*k11 + t^4*k1*k3*k5*k6*k8*k10*k11 + t^4*k1*k4*k5*k6*k7*k8*k11 + t^4*k1*k4*k5*k6*k7*k9*k11 + t^4*k2*k4*k5*k6*k7*k8*k11 + t^4*k2*k4*k5*k6*k7*k9*k11 - t^3*k1*k2*k5*k6*k7*k10*k11*k13 - t^3*k1*k2*k5*k6*k8*k10*k11*k13 + t^3*k1*k2*k5*k7*k8*k10*k11*k14 - t^3*k1*k3*k5*k6*k7*k10*k11*k13 - t^3*k1*k3*k5*k6*k8*k10*k11*k13 + t^3*k1*k3*k5*k7*k8*k10*k11*k14 - t^3*k1*k4*k5*k6*k7*k8*k11*k13 - t^3*k1*k4*k5*k6*k7*k9*k11*k13 + t^3*k1*k4*k5*k7*k8*k9*k11*k14 - t^3*k2*k4*k5*k6*k7*k8*k11*k13 - t^3*k2*k4*k5*k6*k7*k9*k11*k13 + t^3*k2*k4*k5*k7*k8*k9*k11*k14 + t^3*k1*k2*k3*k6*k7*k10*k11 + t^3*k1*k2*k3*k6*k8*k10*k11 + t^3*k1*k2*k5*k6*k7*k8*k11 + t^3*k1*k2*k5*k6*k7*k9*k11 + t^3*k1*k3*k5*k6*k7*k8*k11 + t^3*k1*k3*k5*k6*k7*k9*k11 + t^3*k1*k4*k5*k6*k7*k8*k9 + t^3*k2*k4*k5*k6*k7*k8*k9 - t^2*k1*k2*k3*k6*k7*k10*k11*k13 - t^2*k1*k2*k3*k6*k8*k10*k11*k13 + t^2*k1*k2*k3*k7*k8*k10*k11*k14 - t^2*k1*k2*k5*k6*k7*k8*k11*k13 - t^2*k1*k2*k5*k6*k7*k9*k11*k13 + t^2*k1*k2*k5*k7*k8*k9*k11*k14 - t^2*k1*k3*k5*k6*k7*k8*k11*k13 - t^2*k1*k3*k5*k6*k7*k9*k11*k13 + t^2*k1*k3*k5*k7*k8*k9*k11*k14 - t^2*k1*k4*k5*k6*k7*k8*k9*k13 - t^2*k2*k4*k5*k6*k7*k8*k9*k13 + t^2*k1*k2*k3*k6*k7*k8*k11 + t^2*k1*k2*k3*k6*k7*k9*k11 + t^2*k1*k2*k5*k6*k7*k8*k9 + t^2*k1*k3*k5*k6*k7*k8*k9 - t*k1*k2*k3*k6*k7*k8*k11*k13 - t*k1*k2*k3*k6*k7*k9*k11*k13 + t*k1*k2*k3*k7*k8*k9*k11*k14 - t*k1*k2*k5*k6*k7*k8*k9*k13 - t*k1*k3*k5*k6*k7*k8*k9*k13 + t*k1*k2*k3*k6*k7*k8*k9 - k1*k2*k3*k6*k7*k8*k9*k13)/(t*k1*k2*k5*(t^3*k4*k7*k10*k11 + t^3*k4*k8*k10*k11 + t^2*k3*k7*k10*k11 + t^2*k3*k8*k10*k11 + t^2*k4*k7*k8*k11 + t^2*k4*k7*k9*k11 + t*k3*k7*k8*k11 + t*k3*k7*k9*k11 + t*k4*k7*k8*k9 + k3*k7*k8*k9))
end
@fastmath function J(t,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14)
    return abs(5*(-k1*k4*k5*k6*k7*k10*k11 - k1*k4*k5*k6*k8*k10*k11 - k2*k4*k5*k6*k7*k10*k11 - k2*k4*k5*k6*k8*k10*k11)*t^4 + 4*(-k1*k2*k4*k5*k7*k10*k11*k12 - k1*k2*k4*k5*k8*k10*k11*k12 + k1*k4*k5*k6*k7*k10*k11*k13 + k1*k4*k5*k6*k8*k10*k11*k13 - k1*k4*k5*k7*k8*k10*k11*k14 + k2*k4*k5*k6*k7*k10*k11*k13 + k2*k4*k5*k6*k8*k10*k11*k13 - k2*k4*k5*k7*k8*k10*k11*k14 - k1*k2*k5*k6*k7*k10*k11 - k1*k2*k5*k6*k8*k10*k11 - k1*k3*k5*k6*k7*k10*k11 - k1*k3*k5*k6*k8*k10*k11 - k1*k4*k5*k6*k7*k8*k11 - k1*k4*k5*k6*k7*k9*k11 - k2*k4*k5*k6*k7*k8*k11 - k2*k4*k5*k6*k7*k9*k11)*t^3 + 3*(-k1*k2*k3*k5*k7*k10*k11*k12 - k1*k2*k3*k5*k8*k10*k11*k12 - k1*k2*k4*k5*k7*k8*k11*k12 - k1*k2*k4*k5*k7*k9*k11*k12 + k1*k2*k5*k6*k7*k10*k11*k13 + k1*k2*k5*k6*k8*k10*k11*k13 - k1*k2*k5*k7*k8*k10*k11*k14 + k1*k3*k5*k6*k7*k10*k11*k13 + k1*k3*k5*k6*k8*k10*k11*k13 - k1*k3*k5*k7*k8*k10*k11*k14 + k1*k4*k5*k6*k7*k8*k11*k13 + k1*k4*k5*k6*k7*k9*k11*k13 - k1*k4*k5*k7*k8*k9*k11*k14 + k2*k4*k5*k6*k7*k8*k11*k13 + k2*k4*k5*k6*k7*k9*k11*k13 - k2*k4*k5*k7*k8*k9*k11*k14 - k1*k2*k3*k6*k7*k10*k11 - k1*k2*k3*k6*k8*k10*k11 - k1*k2*k5*k6*k7*k8*k11 - k1*k2*k5*k6*k7*k9*k11 - k1*k3*k5*k6*k7*k8*k11 - k1*k3*k5*k6*k7*k9*k11 - k1*k4*k5*k6*k7*k8*k9 - k2*k4*k5*k6*k7*k8*k9)*t^2 + 2*(-k1*k2*k3*k5*k7*k8*k11*k12 - k1*k2*k3*k5*k7*k9*k11*k12 + k1*k2*k3*k6*k7*k10*k11*k13 + k1*k2*k3*k6*k8*k10*k11*k13 - k1*k2*k3*k7*k8*k10*k11*k14 - k1*k2*k4*k5*k7*k8*k9*k12 + k1*k2*k5*k6*k7*k8*k11*k13 + k1*k2*k5*k6*k7*k9*k11*k13 - k1*k2*k5*k7*k8*k9*k11*k14 + k1*k3*k5*k6*k7*k8*k11*k13 + k1*k3*k5*k6*k7*k9*k11*k13 - k1*k3*k5*k7*k8*k9*k11*k14 + k1*k4*k5*k6*k7*k8*k9*k13 + k2*k4*k5*k6*k7*k8*k9*k13 - k1*k2*k3*k6*k7*k8*k11 - k1*k2*k3*k6*k7*k9*k11 - k1*k2*k5*k6*k7*k8*k9 - k1*k3*k5*k6*k7*k8*k9)*t - k1*k2*k3*k5*k7*k8*k9*k12 + k1*k2*k3*k6*k7*k8*k11*k13 + k1*k2*k3*k6*k7*k9*k11*k13 - k1*k2*k3*k7*k8*k9*k11*k14 + k1*k2*k5*k6*k7*k8*k9*k13 + k1*k3*k5*k6*k7*k8*k9*k13 - k1*k2*k3*k6*k7*k8*k9)
end
@fastmath function IntegrandWithoutCOEFF(t,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,kk12,k13,k14)
    return J(t,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,g(t,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k13,k14),k13,k14)*chi(kk12[1],kk12[2],g(t,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k13,k14))/(t*k1*k2*k5*(t^3*k4*k7*k10*k11 + t^3*k4*k8*k10*k11 + t^2*k3*k7*k10*k11 + t^2*k3*k8*k10*k11 + t^2*k4*k7*k8*k11 + t^2*k4*k7*k9*k11 + t*k3*k7*k8*k11 + t*k3*k7*k9*k11 + t*k4*k7*k8*k9 + k3*k7*k8*k9))
end
@fastmath function sumo_with_SE_t(B,NN)
    kk1=B[1]
    kk2=B[2]
    kk3=B[3]
    kk4=B[4]
    kk5=B[5]
    kk6=B[6]
    kk7=B[7]
    kk8=B[8]
    kk9=B[9]
    kk10=B[10]
    kk11=B[11]
    kk12=B[12]
    kk13=B[13]
    kk14=B[14]
    S=0
    Ians=0
    st1=time_ns()
    COEFF=kk12[2]/(kk12[2]-kk12[1])
    t=samplor(0,kk12[2])
    k1=samplor(kk1[1],kk1[2])
    k2=samplor(kk2[1],kk2[2])
    k3=samplor(kk3[1],kk3[2])
    k4=samplor(kk4[1],kk4[2])
    k5=samplor(kk5[1],kk5[2])
    k6=samplor(kk6[1],kk6[2])
    k7=samplor(kk7[1],kk7[2])
    k8=samplor(kk8[1],kk8[2])
    k9=samplor(kk9[1],kk9[2])
    k10=samplor(kk10[1],kk10[2])
    k11=samplor(kk11[1],kk11[2])
    k13=samplor(kk13[1],kk13[2])
    k14=samplor(kk14[1],kk14[2])
    Ians+=IntegrandWithoutCOEFF(t,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,kk12,k13,k14)
    Ians=COEFF*Ians
    @inbounds for i=2:NN
        t=samplor(0,kk12[2])
        k1=samplor(kk1[1],kk1[2])
        k2=samplor(kk2[1],kk2[2])
        k3=samplor(kk3[1],kk3[2])
        k4=samplor(kk4[1],kk4[2])
        k5=samplor(kk5[1],kk5[2])
        k6=samplor(kk6[1],kk6[2])
        k7=samplor(kk7[1],kk7[2])
        k8=samplor(kk8[1],kk8[2])
        k9=samplor(kk9[1],kk9[2])
        k10=samplor(kk10[1],kk10[2])
        k11=samplor(kk11[1],kk11[2])
        k13=samplor(kk13[1],kk13[2])
        k14=samplor(kk14[1],kk14[2])
        delto=COEFF*IntegrandWithoutCOEFF(t,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,kk12,k13,k14)-Ians
        Ians+=delto/i
        S+=delto^2*(i-1)/i
    end
    Estandardo=S/NN
    Estandardo=Estandardo/(NN-1)
    Estandardo=sqrt(Estandardo)
    st2=time_ns()
    return Ians,Estandardo,(st2-st1)/10^9
end
@fastmath function sumo_antithetic_with_SE_t(B,NN)
    kk1=B[1]
    kk2=B[2]
    kk3=B[3]
    kk4=B[4]
    kk5=B[5]
    kk6=B[6]
    kk7=B[7]
    kk8=B[8]
    kk9=B[9]
    kk10=B[10]
    kk11=B[11]
    kk12=B[12]
    kk13=B[13]
    kk14=B[14]
    S=0
    Ians=0
    st1=time_ns()
    centero=[(0+kk12[2])/2,(kk1[1]+kk1[2])/2,(kk2[1]+kk2[2])/2,(kk3[1]+kk3[2])/2,(kk4[1]+kk4[2])/2,(kk5[1]+kk5[2])/2,(kk6[1]+kk6[2])/2,(kk7[1]+kk7[2])/2,(kk8[1]+kk8[2])/2,(kk9[1]+kk9[2])/2,(kk10[1]+kk10[2])/2,(kk11[1]+kk11[2])/2,(kk13[1]+kk13[2])/2,(kk14[1]+kk14[2])/2]
    COEFF=kk12[2]/(kk12[2]-kk12[1])
    t=samplor(0,kk12[2])
    k1=samplor(kk1[1],kk1[2])
    k2=samplor(kk2[1],kk2[2])
    k3=samplor(kk3[1],kk3[2])
    k4=samplor(kk4[1],kk4[2])
    k5=samplor(kk5[1],kk5[2])
    k6=samplor(kk6[1],kk6[2])
    k7=samplor(kk7[1],kk7[2])
    k8=samplor(kk8[1],kk8[2])
    k9=samplor(kk9[1],kk9[2])
    k10=samplor(kk10[1],kk10[2])
    k11=samplor(kk11[1],kk11[2])
    k13=samplor(kk13[1],kk13[2])
    k14=samplor(kk14[1],kk14[2])
    Ians+=IntegrandWithoutCOEFF(t,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,kk12,k13,k14)
    Ians=COEFF*Ians
    t=2*centero[1]-t
    k1=2*centero[2]-k1
    k2=2*centero[3]-k2
    k3=2*centero[4]-k3
    k4=2*centero[5]-k4
    k5=2*centero[6]-k5
    k6=2*centero[7]-k6
    k7=2*centero[8]-k7
    k8=2*centero[9]-k8
    k9=2*centero[10]-k9
    k10=2*centero[11]-k10
    k11=2*centero[12]-k11
    k13=2*centero[13]-k13
    k14=2*centero[14]-k14
    delto=COEFF*IntegrandWithoutCOEFF(t,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,kk12,k13,k14)-Ians
    Ians+=delto/2
    S+=delto^2/2
    @inbounds for i=2:floor(NN/2)
        t=samplor(0,kk12[2])
        k1=samplor(kk1[1],kk1[2])
        k2=samplor(kk2[1],kk2[2])
        k3=samplor(kk3[1],kk3[2])
        k4=samplor(kk4[1],kk4[2])
        k5=samplor(kk5[1],kk5[2])
        k6=samplor(kk6[1],kk6[2])
        k7=samplor(kk7[1],kk7[2])
        k8=samplor(kk8[1],kk8[2])
        k9=samplor(kk9[1],kk9[2])
        k10=samplor(kk10[1],kk10[2])
        k11=samplor(kk11[1],kk11[2])
        k13=samplor(kk13[1],kk13[2])
        k14=samplor(kk14[1],kk14[2])
        delto=COEFF*IntegrandWithoutCOEFF(t,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,kk12,k13,k14)-Ians
        Ians+=delto/(2*i+1)
        S+=(delto^2)*2*i/(2*i+1)
        t=2*centero[1]-t
        k1=2*centero[2]-k1
        k2=2*centero[3]-k2
        k3=2*centero[4]-k3
        k4=2*centero[5]-k4
        k5=2*centero[6]-k5
        k6=2*centero[7]-k6
        k7=2*centero[8]-k7
        k8=2*centero[9]-k8
        k9=2*centero[10]-k9
        k10=2*centero[11]-k10
        k11=2*centero[12]-k11
        k13=2*centero[13]-k13
        k14=2*centero[14]-k14
        delto=COEFF*IntegrandWithoutCOEFF(t,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,kk12,k13,k14)-Ians
        Ians+=delto/(2*i+2)
        S+=(delto^2)*(2*i+1)/(2*i+2)
    end
    Estandardo=S/NN
    Estandardo=Estandardo/(NN-1)
    Estandardo=sqrt(Estandardo)
    st2=time_ns()
    return Ians,Estandardo,(st2-st1)/10^9
end

# Simple and antithetic Monte Carlo integrations on the box [0,100]^{14}. This is not mentioned on the paper.
Binput=[[0.0,100.0],[0.0,100.0],[0.0,100.0],[0.0,100.0],[0.0,100.0],[0.0,100.0],[0.0,100.0],[0.0,100.0],[0.0,100.0],[0.0,100.0],[0.0,100.0],[0.0,100.0],[0.0,100.0],[0.0,100.0]]
print("\nSimple Monte Carlo for the double HK network with B = $Binput\n")

for n=1:8
    print("N=10^$n\n")
    X1=sumo_with_SE_t(Binput,10^n)
    print("Ians = $(X1[1]), Standard error = $(X1[2]), Computation time = $(X1[3]).\n")
end

# output
#=
Simple Monte Carlo for the double HK network with B = [[0.0, 100.0], [0.0, 100.0], [0.0, 100.0], [0.0, 100.0], [0.0, 100.0], [0.0, 100.0], [0.0, 100.0], [0.0, 100.0], [0.0, 100.0], [0.0, 100.0], [0.0, 100.0], [0.0, 100.0], [0.0, 100.0], [0.0, 100.0]]
N=10^1
Ians = 0.2289013418102575, Standard error = 0.16287482921382365, Computation time = 1.5e-5.
N=10^2
Ians = 0.6075167005656528, Standard error = 0.24217857913655613, Computation time = 7.34e-5.
N=10^3
Ians = 0.605468532341809, Standard error = 0.064734937290308, Computation time = 0.0006557.
N=10^4
Ians = 1.0748171601255523, Standard error = 0.25001178091506737, Computation time = 0.0061756.
N=10^5
Ians = 0.85937462386121, Standard error = 0.035338565250521335, Computation time = 0.0584935.
N=10^6
Ians = 0.9762076316948526, Standard error = 0.023705597136096566, Computation time = 0.575894.
N=10^7
Ians = 0.9866778717971968, Standard error = 0.011355576692058274, Computation time = 5.7569765.
N=10^8
Ians = 0.9965856174263001, Standard error = 0.004466006384699092, Computation time = 58.3942104.
=#

# Simple and antithetic Monte Carlo integrations on the box B in the paper.
Binput=[[0.0,2.0],[0.0,100.0],[0.0,50.0],[0.0,10.0],[0.0,100.0],[0.0,50.0],[0.0,5.0],[50.0,100.0],[0.0,100.0],[0.0,100.0],[0.0,100.0],[0.0,100.0],[0.0,50.0],[0.0,100.0]]
print("\nSimple Monte Carlo for the double HK network with B = $Binput\n")
for n=1:8
    print("N=10^$n\n")
    X1=sumo_with_SE_t(Binput,10^n)
    print("Ians = $(X1[1]), Standard error = $(X1[2]), Computation time = $(X1[3]).\n")
end
print("\nAntithetic Monte Carlo for the double HK network with B = $Binput\n")
for n=1:8
    print("N=10^$n\n")
    X1=sumo_antithetic_with_SE_t(Binput,10^n)
    print("Ians = $(X1[1]), Standard error = $(X1[2]), Computation time = $(X1[3]).\n")
end

# output
#=
Simple Monte Carlo for the double HK network with B = [[0.0, 2.0], [0.0, 100.0], [0.0, 50.0], [0.0, 10.0], [0.0, 100.0], [0.0, 50.0], [0.0, 5.0], [50.0, 100.0], [0.0, 100.0], [0.0, 100.0], [0.0, 100.0], [0.0, 100.0], [0.0, 50.0], [0.0, 100.0]]
N=10^1
Ians = 3.9581028460088534, Standard error = 3.782118272384846, Computation time = 5.6e-5.
N=10^2
Ians = 0.589605909697324, Standard error = 0.2902622502720907, Computation time = 7.32e-5.
N=10^3
Ians = 0.7140678283535945, Standard error = 0.1773641198543826, Computation time = 0.000598701.
N=10^4
Ians = 0.8229949422878612, Standard error = 0.08243057250361314, Computation time = 0.005898699.
N=10^5
Ians = 1.358263316984987, Standard error = 0.2011848997694101, Computation time = 0.059403.
N=10^6
Ians = 1.158746585378628, Standard error = 0.06737975085076364, Computation time = 0.580036399.
N=10^7
Ians = 1.2268710260627829, Standard error = 0.045539648082807326, Computation time = 6.107687399.
N=10^8
Ians = 1.219883154558217, Standard error = 0.011381225608800206, Computation time = 57.6215047.

Antithetic Monte Carlo for the double HK network with B = [[0.0, 2.0], [0.0, 100.0], [0.0, 50.0], [0.0, 10.0], [0.0, 100.0], [0.0, 50.0], [0.0, 5.0], [50.0, 100.0], [0.0, 100.0], [0.0, 100.0], [0.0, 100.0], [0.0, 100.0], [0.0, 50.0], [0.0, 100.0]]
N=10^1
Ians = 0.0, Standard error = 0.0, Computation time = 2.17e-5.
N=10^2
Ians = 0.6882582878600081, Standard error = 0.4345402382164831, Computation time = 0.000108301.
N=10^3
Ians = 1.0382623531913004, Standard error = 0.43824402150440184, Computation time = 0.000541.
N=10^4
Ians = 0.782517478594255, Standard error = 0.07761314880794201, Computation time = 0.0055485.
N=10^5
Ians = 1.2553104005941382, Standard error = 0.11983555981866702, Computation time = 0.0521221.
N=10^6
Ians = 1.2347775094598001, Standard error = 0.08245748286069149, Computation time = 0.474456999.
N=10^7
Ians = 1.2281609706932848, Standard error = 0.03754188585646872, Computation time = 4.781883001.
N=10^8
Ians = 1.2130190093760627, Standard error = 0.008647718311013693, Computation time = 48.446972701.
=#

# Bisect search
function SameBox(B)
    return(deepcopy(B))
end

@fastmath function fun_doubleHK(B) # The same as the sumo_antithetic_with_SE_t, but the output doesn't containt e-hat and the computation time.
    NN=10^8
    kk1=B[1]
    kk2=B[2]
    kk3=B[3]
    kk4=B[4]
    kk5=B[5]
    kk6=B[6]
    kk7=B[7]
    kk8=B[8]
    kk9=B[9]
    kk10=B[10]
    kk11=B[11]
    kk12=B[12]
    kk13=B[13]
    kk14=B[14]
    Ians=0
    centero=[(0+kk12[2])/2,(kk1[1]+kk1[2])/2,(kk2[1]+kk2[2])/2,(kk3[1]+kk3[2])/2,(kk4[1]+kk4[2])/2,(kk5[1]+kk5[2])/2,(kk6[1]+kk6[2])/2,(kk7[1]+kk7[2])/2,(kk8[1]+kk8[2])/2,(kk9[1]+kk9[2])/2,(kk10[1]+kk10[2])/2,(kk11[1]+kk11[2])/2,(kk13[1]+kk13[2])/2,(kk14[1]+kk14[2])/2]
    COEFF=kk12[2]/(kk12[2]-kk12[1])
    t=samplor(0,kk12[2])
    k1=samplor(kk1[1],kk1[2])
    k2=samplor(kk2[1],kk2[2])
    k3=samplor(kk3[1],kk3[2])
    k4=samplor(kk4[1],kk4[2])
    k5=samplor(kk5[1],kk5[2])
    k6=samplor(kk6[1],kk6[2])
    k7=samplor(kk7[1],kk7[2])
    k8=samplor(kk8[1],kk8[2])
    k9=samplor(kk9[1],kk9[2])
    k10=samplor(kk10[1],kk10[2])
    k11=samplor(kk11[1],kk11[2])
    k13=samplor(kk13[1],kk13[2])
    k14=samplor(kk14[1],kk14[2])
    Ians+=IntegrandWithoutCOEFF(t,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,kk12,k13,k14)
    Ians=COEFF*Ians
    t=2*centero[1]-t
    k1=2*centero[2]-k1
    k2=2*centero[3]-k2
    k3=2*centero[4]-k3
    k4=2*centero[5]-k4
    k5=2*centero[6]-k5
    k6=2*centero[7]-k6
    k7=2*centero[8]-k7
    k8=2*centero[9]-k8
    k9=2*centero[10]-k9
    k10=2*centero[11]-k10
    k11=2*centero[12]-k11
    k13=2*centero[13]-k13
    k14=2*centero[14]-k14
    delto=COEFF*IntegrandWithoutCOEFF(t,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,kk12,k13,k14)-Ians
    Ians+=delto/2
    @inbounds for i=2:floor(NN/2)
        t=samplor(0,kk12[2])
        k1=samplor(kk1[1],kk1[2])
        k2=samplor(kk2[1],kk2[2])
        k3=samplor(kk3[1],kk3[2])
        k4=samplor(kk4[1],kk4[2])
        k5=samplor(kk5[1],kk5[2])
        k6=samplor(kk6[1],kk6[2])
        k7=samplor(kk7[1],kk7[2])
        k8=samplor(kk8[1],kk8[2])
        k9=samplor(kk9[1],kk9[2])
        k10=samplor(kk10[1],kk10[2])
        k11=samplor(kk11[1],kk11[2])
        k13=samplor(kk13[1],kk13[2])
        k14=samplor(kk14[1],kk14[2])
        delto=COEFF*IntegrandWithoutCOEFF(t,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,kk12,k13,k14)-Ians
        Ians+=delto/(2*i+1)
        t=2*centero[1]-t
        k1=2*centero[2]-k1
        k2=2*centero[3]-k2
        k3=2*centero[4]-k3
        k4=2*centero[5]-k4
        k5=2*centero[6]-k5
        k6=2*centero[7]-k6
        k7=2*centero[8]-k7
        k8=2*centero[9]-k8
        k9=2*centero[10]-k9
        k10=2*centero[11]-k10
        k11=2*centero[12]-k11
        k13=2*centero[13]-k13
        k14=2*centero[14]-k14
        delto=COEFF*IntegrandWithoutCOEFF(t,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,kk12,k13,k14)-Ians
        Ians+=delto/(2*i+2)
    end
    return Ians
end

# Now the bisect search computation mentioned in the paper.

@fastmath function NewBisectSearch_t(fun, B, minimum_size, sh_terminate)
    st1=time_ns()
    Bchosen=SameBox(B)
    anschosen=fun(Bchosen)
    if anschosen > 2.95
        return(Bchosen,anschosen,1,0,(time_ns()-st1)/(10^9))
    end
    d = length(B)
    shlist=[]
	for i=1:d
		push!(shlist,i)
	end
	for i=1:d-1
		push!(shlist,i)
	end
    sh2=0
    sh1=0
    @inbounds while sh2 != sh_terminate
        sh1+=1
        if sh1 >= d+1
            sh1=1
        end
        sizeend=1
		for i=sh1:sh1+d-1
			if Bchosen[shlist[i]][2]-Bchosen[shlist[i]][1] >= 2*minimum_size[shlist[i]]
                #print("\n$(B[shlist[i]][2])-$(B[shlist[i]][1])=$(B[shlist[i]][2]-B[shlist[i]][1]) >= $(2*minimum_size[shlist[i]]) is it true?$(B[shlist[i]][2]-B[shlist[i]][1] >= 2*minimum_size[shlist[i]])")
				sizeend=0
                sh2+=1
                B1 = SameBox(Bchosen)
                B1[shlist[i]][2] = (Bchosen[shlist[i]][1] + Bchosen[shlist[i]][2]) / 2
                ans1=fun(B1)
                B2 = SameBox(Bchosen)
                B2[shlist[i]][1] = (Bchosen[shlist[i]][1] + Bchosen[shlist[i]][2]) / 2
                ans2=fun(B2)
                if ans1 >= ans2
                    Bchosen=SameBox(B1)
                    anschosen=ans1
                else
                    Bchosen=SameBox(B2)
                    anschosen=ans2
                end
                if anschosen > 2.95
                    return(Bchosen,anschosen,1+2*sh2,sh2,(time_ns()-st1)/(10^9))
                end
                sh1=shlist[i]
				break
			end
		end
		if sizeend==1
			return(Bchosen,anschosen,1+2*sh2,sh2,(time_ns()-st1)/(10^9))
		end
    end
    st2=time_ns()
    return(Bchosen,anschosen,1+2*sh2,sh2,(st2-st1)/(10^9))
end

bisectsteps=50
Binput=[[0.0,2.0],[0.0,100.0],[0.0,50.0],[0.0,10.0],[0.0,100.0],[0.0,50.0],[0.0,5.0],[50.0,100.0],[0.0,100.0],[0.0,100.0],[0.0,100.0],[0.0,100.0],[0.0,50.0],[0.0,100.0]]
minimumsizeinput=[0.1,10,5,1,10,5,0.5,5,10,10,10,25,25,50]
print("\nThe bisect search for the double HK 1v 14p with the initial box\nB = $Binput\nand with restriction to at most $bisectsteps steps of bisection and minimum size of the sub-boxes' edges' lengths $minimumsizeinput.\n")
X=NewBisectSearch_t(fun_doubleHK,Binput,minimumsizeinput,bisectsteps)
print("\nThe output of the search is the following sub-box\n$(X[1])\nwith the integral equal to $(X[2]) and computed at $(X[5]) seconds.\nNumber of the integrals computed is $(X[3]) in $(X[4]) bisect steps.")

# output
#=
The bisect search for the double HK 1v 14p with the initial box
B = [[0.0, 2.0], [0.0, 100.0], [0.0, 50.0], [0.0, 10.0], [0.0, 100.0], [0.0, 50.0], [0.0, 5.0], [50.0, 100.0], [0.0, 100.0], [0.0, 100.0], [0.0, 100.0], [0.0, 100.0], [0.0, 50.0], [0.0, 100.0]]
and with restriction to at most 50 steps of bisection and minimum size of the sub-boxes' edges' lengths [0.1, 10.0, 5.0, 1.0, 10.0, 5.0, 0.5, 5.0, 10.0, 10.0, 10.0, 25.0, 25.0, 50.0].

The output of the search is the following sub-box
[[0.25, 0.5], [75.0, 87.5], [43.75, 50.0], [3.75, 5.0], [87.5, 100.0], [18.75, 25.0], [1.25, 2.5], [62.5, 75.0], [75.0, 100.0], [25.0, 50.0], [75.0, 100.0], [75.0, 100.0], [25.0, 50.0], [50.0, 100.0]]
with the integral equal to 2.9581494511520465 and computed at 3082.6467648 seconds.
Number of the integrals computed is 65 in 32 bisect steps.
=#

# Box B.
print(sumo_antithetic_with_SE_t([[0.0,2.0],[0.0,100.0],[0.0,50.0],[0.0,10.0],[0.0,100.0],[0.0,50.0],[0.0,5.0],[50.0,100.0],[0.0,100.0],[0.0,100.0],[0.0,100.0],[0.0,100.0],[0.0,50.0],[0.0,100.0]],10^8))
# (1.2146747624743393, 0.010183440402819955, 46.9366298))

# Box C.
print(sumo_antithetic_with_SE_t([[0.25, 0.5], [75.0, 87.5], [43.75, 50.0], [3.75, 5.0], [87.5, 100.0], [18.75, 25.0], [1.25, 2.5], [62.5, 75.0], [75.0, 100.0], [25.0, 50.0], [75.0, 100.0], [75.0, 100.0], [25.0, 50.0], [50.0, 100.0]],10^8))
# (2.9548415358024123, 0.008737881094094789, 47.252372201)

# The end of the file.
