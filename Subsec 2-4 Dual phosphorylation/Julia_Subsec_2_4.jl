# The Monte Carlo integrations of the example in subsection 2.4 of the paper.

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

@fastmath function g1(x,e,f,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12)
    return k12*k10*(k2 + k3)*f^2*x/(k3*k1*(k11 + k12)) + x + k6*k4*(k8 + k9)*e^2*x/(k9*k7*(k5 + k6)) + k12*k10*e*f^2*x/(k3*(k11 + k12)) + k4*e^2*f*x/(k5 + k6) + k4*e^2*f*x/(k9*(k5 + k6)) + k10*e*f^2*x/(k11 + k12)
end
@fastmath function g2(x,e,f,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12)
    return e^2*f + k12*k10*e*f^2*x/(k3*(k11 + k12)) + k4*e^2*f*x/(k5 + k6)
end
@fastmath function g3(x,e,f,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12)
    return e*f^2 + k4*e^2*f*x/(k9*(k5 + k6)) + k10*e*f^2*x/(k11 + k12)
end
@fastmath function J(x,e,f,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12)
    return abs(2*(((3*f*(k1*k4*(k7*(k5 + k6)*(x*k10 + k11 + k12)*f + (((-x*k4 + 3*k5 + 3*k6)*k12 + (5*x*k10 + 3*k11)*k6 + x*(x*k4 + 5*k5)*k10 - k11*(x*k4 - 3*k5))*k6)/3)*(k11 + k12)*e^2 + k7*k10*f^2*k1*(k5 + k6)*(x*k4 + k5 + k6)*(k11 + k12)*e + k7*(k12*((5*x*k4 + 3*k5 + 3*k6)*k12 + (-x*k10 + 3*k11)*k6 + x*(x*k4 - k5)*k10 + 5*(x*k4 + (3*k5)/5)*k11)*k10*f^2 + 3*k1*(x*k4 + k5 + k6)*(k11 + k12)*(x*k10 + k11 + k12))*(k5 + k6)/3)*k9^2)/2 + e*(x*k1*k4*k6*(k11 + k12)*e^2 + (3*f*(k7*(k5 + k6)*(k11 + k12)*f + (((-x*k4 + 3*k5 + 3*k6)*k12 + (5*x*k10 + 3*k11)*k6 + x*(x*k4 + 5*k5)*k10 - k11*(x*k4 - 3*k5))*k8*k6)/3)*k1*e)/2 - 2*k7*k10*k12*f^2*x*(k5 + k6))*k4*(k11 + k12)*k9 + e^3*k8*x*k1*k4^2*k6*(k11 + k12)^2)*e*k3^2 - (3*k12*(-f*(k5 + k6)*((k7*(k5 + k6)*f - (4*x*k4*k6)/3)*k1*(k11 + k12)*e^2 + k7*k2*((5*x*k4 + 3*k5 + 3*k6)*k12 + (-x*k10 + 3*k11)*k6 + x*(x*k4 - k5)*k10 + 5*(x*k4 + (3*k5)/5)*k11)*f*e/3 + (2*k7*k10*k12*f^2*x*(k5 + k6))/3)*k9^2 + x*e*k4*((k7*(k5 + k6)*f + x*k4*k6/3)*k1*(k11 + k12)*e^2 + f*(k5 + k6)*(k11 + k12)*(k7*f*k1 + 4/3*k8*k1*k6 + 4/3*k7*k2)*e + x*k7*(k5 + k6)*(k10*k12*f^2 + 3*k1*(k11 + k12))/3)*k9 + e^3*x^2*k1*k4^2*k6*k8*(k11 + k12)/3)*f*k10*k3)/2 - x*k7*k2*k9*k12^2*(-2*f*(k5 + k6)*k9 + k4*e*x)*f^3*k10^2*(k5 + k6)/2)*e*f/((k5 + k6)^2*k7*k9^2*(k11 + k12)^2*k1*k3^2))
end
@fastmath function IntegrandWithoutCOEFF(x,e,f,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,TT1,TT2,TT3)
    return J(x,e,f,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12)*chi(TT1[1],TT1[2],g1(x,e,f,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12))*chi(TT2[1],TT2[2],g2(x,e,f,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12))*chi(TT3[1],TT3[2],g3(x,e,f,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12))
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
    TT1=B[13]
    TT2=B[14]
    TT3=B[15]
    S=0
    Ians=0
    st1=time_ns()
    COEFF=TT1[2]*TT2[2]*TT3[2]/((TT1[2]-TT1[1])*(TT2[2]-TT2[1])*(TT3[2]-TT3[1]))
    x=samplor(0,TT1[2])
    e=samplor(0,TT2[2])
    f=samplor(0,TT3[2])
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
    k12=samplor(kk12[1],kk12[2])
    Ians+=IntegrandWithoutCOEFF(x,e,f,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,TT1,TT2,TT3)
    Ians=COEFF*Ians
    @inbounds for i=2:NN
        x=samplor(0,TT1[2])
        e=samplor(0,TT2[2])
        f=samplor(0,TT3[2])
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
        k12=samplor(kk12[1],kk12[2])
        delto=COEFF*IntegrandWithoutCOEFF(x,e,f,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,TT1,TT2,TT3)-Ians
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
    TT1=B[13]
    TT2=B[14]
    TT3=B[15]
    S=0
    Ians=0
    st1=time_ns()
    centero=[(0+TT1[2])/2,(0+TT2[2])/2,(0+TT3[2])/2,(kk1[1]+kk1[2])/2,(kk2[1]+kk2[2])/2,(kk3[1]+kk3[2])/2,(kk4[1]+kk4[2])/2,(kk5[1]+kk5[2])/2,(kk6[1]+kk6[2])/2,(kk7[1]+kk7[2])/2,(kk8[1]+kk8[2])/2,(kk9[1]+kk9[2])/2,(kk10[1]+kk10[2])/2,(kk11[1]+kk11[2])/2,(kk12[1]+kk12[2])/2]
    COEFF=TT1[2]*TT2[2]*TT3[2]/((TT1[2]-TT1[1])*(TT2[2]-TT2[1])*(TT3[2]-TT3[1]))
    x=samplor(0,TT1[2])
    e=samplor(0,TT2[2])
    f=samplor(0,TT3[2])
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
    k12=samplor(kk12[1],kk12[2])
    Ians+=IntegrandWithoutCOEFF(x,e,f,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,TT1,TT2,TT3)
    Ians=COEFF*Ians
    x=2*centero[1]-x
    e=2*centero[2]-e
    f=2*centero[3]-f
    k1=2*centero[4]-k1
    k2=2*centero[5]-k2
    k3=2*centero[6]-k3
    k4=2*centero[7]-k4
    k5=2*centero[8]-k5
    k6=2*centero[9]-k6
    k7=2*centero[10]-k7
    k8=2*centero[11]-k8
    k9=2*centero[12]-k9
    k10=2*centero[13]-k10
    k11=2*centero[14]-k11
    k12=2*centero[15]-k12
    delto=COEFF*IntegrandWithoutCOEFF(x,e,f,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,TT1,TT2,TT3)-Ians
    Ians+=delto/2
    S+=delto^2/2
    @inbounds for i=2:floor(NN/2)
        x=samplor(0,TT1[2])
        e=samplor(0,TT2[2])
        f=samplor(0,TT3[2])
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
        k12=samplor(kk12[1],kk12[2])
        delto=COEFF*IntegrandWithoutCOEFF(x,e,f,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,TT1,TT2,TT3)-Ians
        Ians+=delto/(2*i+1)
        S+=(delto^2)*2*i/(2*i+1)
        x=2*centero[1]-x
        e=2*centero[2]-e
        f=2*centero[3]-f
        k1=2*centero[4]-k1
        k2=2*centero[5]-k2
        k3=2*centero[6]-k3
        k4=2*centero[7]-k4
        k5=2*centero[8]-k5
        k6=2*centero[9]-k6
        k7=2*centero[10]-k7
        k8=2*centero[11]-k8
        k9=2*centero[12]-k9
        k10=2*centero[13]-k10
        k11=2*centero[14]-k11
        k12=2*centero[15]-k12
        delto=COEFF*IntegrandWithoutCOEFF(x,e,f,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,TT1,TT2,TT3)-Ians
        Ians+=delto/(2*i+2)
        S+=(delto^2)*(2*i+1)/(2*i+2)
    end
    Estandardo=S/NN
    Estandardo=Estandardo/(NN-1)
    Estandardo=sqrt(Estandardo)
    st2=time_ns()
    return Ians,Estandardo,(st2-st1)/10^9
end

# Simple and Antithetic Monte Carlo Integration for the first box which we expect to be equal to 1.
Binput=[[500.0,1000.0],[25.0,50.0],[25.0,50.0],[5.0,10.0],[5.0,10.0],[5.0,10.0],[5.0,10.0],[1.0,2.0],[1.0,2.0],[5.0,10.0],[50.0,100.0],[50.0,100.0],[1.0,2.0],[2.0,4.0],[2.0,4.0]]
print("\nSimple Monte Carlo for the 2-site network with B = $Binput\n")

for n=1:8
    print("N=10^$n\n")
    X1=sumo_with_SE_t(Binput,10^n)
    print("Ians = $(X1[1]), Standard error = $(X1[2]), Computation time = $(X1[3]).\n")
end

print("\nAntithetic Monte Carlo for the 2-site network with B = $Binput\n")

for n=1:8
    print("N=10^$n\n")
    X1=sumo_antithetic_with_SE_t(Binput,10^n)
    print("Ians = $(X1[1]), Standard error = $(X1[2]), Computation time = $(X1[3]).\n")
end

# output
#=
Simple Monte Carlo for the 2-site network with B = [[500.0, 1000.0], [25.0, 50.0], [25.0, 50.0], [5.0, 10.0], [5.0, 10.0], [5.0, 10.0], [5.0, 10.0], [1.0, 2.0], [1.0, 2.0], [5.0, 10.0], [50.0, 100.0], [50.0, 100.0], [1.0, 2.0], [2.0, 4.0], [2.0, 4.0]]
N=10^1
Ians = 0.0, Standard error = 0.0, Computation time = 7.15e-5.
N=10^2
Ians = 0.0, Standard error = 0.0, Computation time = 7.59e-5.
N=10^3
Ians = 1.9246691575156003, Standard error = 0.980449858527484, Computation time = 0.0005663.
N=10^4
Ians = 1.2688530103402198, Standard error = 0.2629169441190915, Computation time = 0.0056389.
N=10^5
Ians = 0.9933094949012656, Standard error = 0.07240617222050184, Computation time = 0.065468.
N=10^6
Ians = 0.9664963985076391, Standard error = 0.022245033654315842, Computation time = 0.4414969.
N=10^7
Ians = 0.991698386800196, Standard error = 0.007162336897001736, Computation time = 4.3432891.
N=10^8
Ians = 0.9991422988787995, Standard error = 0.002269119957907759, Computation time = 46.527284599.

Antithetic Monte Carlo for the 2-site network with B = [[500.0, 1000.0], [25.0, 50.0], [25.0, 50.0], [5.0, 10.0], [5.0, 10.0], [5.0, 10.0], [5.0, 10.0], [1.0, 2.0], [1.0, 2.0], [5.0, 10.0], [50.0, 100.0], [50.0, 100.0], [1.0, 2.0], [2.0, 4.0], [2.0, 4.0]]
N=10^1
Ians = 0.0, Standard error = 0.0, Computation time = 2.99e-5.
N=10^2
Ians = 0.0, Standard error = 0.0, Computation time = 8.24e-5.
N=10^3
Ians = 1.6409068367455981, Standard error = 0.8400241364979224, Computation time = 0.0007407.
N=10^4
Ians = 0.7707592847006028, Standard error = 0.20385294086516184, Computation time = 0.0041595.
N=10^5
Ians = 1.0030164214522057, Standard error = 0.0711536372779283, Computation time = 0.0364178.
N=10^6
Ians = 1.0205200826261092, Standard error = 0.023102462559498965, Computation time = 0.3366299.
N=10^7
Ians = 1.0063110119512777, Standard error = 0.007193884708942684, Computation time = 3.281446799.
N=10^8
Ians = 1.0027470413694393, Standard error = 0.0022735907562196666, Computation time = 28.6470403.
=#

# Simple and Antithetic Monte Carlo Integration for the second box which we expect to be more than 1. 
Binput=[[0.5,1.5],[509.5,510.5],[1.5,2.5],[1.5,2.5],[0.5,1.5],[0.5,1.5],[1.5,2.5],[0.5,1.5],[0.5,1.5],[1.5,2.5],[0.5,1.5],[0.5,1.5],[110.0,120.0],[20.0,30.0],[15.0,25.0]]
print("\nSimple Monte Carlo for the 2-site network with B = $Binput\n")

for n=1:11
    print("N=10^$n\n")
    X1=sumo_antithetic_with_SE_t(Binput,10^n)
    print("Ians = $(X1[1]), Standard error = $(X1[2]), Computation time = $(X1[3]).\n")
end

# Output
#=
Simple Monte Carlo for the 2-site network with B = [[0.5, 1.5], [509.5, 510.5], [1.5, 2.5], [1.5, 2.5], [0.5, 1.5], [0.5, 1.5], [1.5, 2.5], [0.5, 1.5], [0.5, 1.5], [1.5, 2.5], [0.5, 1.5], [0.5, 1.5], [110.0, 120.0], [20.0, 30.0], [15.0, 25.0]]
N=10^1
Ians = 0.0, Standard error = 0.0, Computation time = 0.020249399.
N=10^2
Ians = 0.0, Standard error = 0.0, Computation time = 3.9501e-5.
N=10^3
Ians = 0.0, Standard error = 0.0, Computation time = 0.0002897.
N=10^4
Ians = 0.0, Standard error = 0.0, Computation time = 0.002881401.
N=10^5
Ians = 0.0, Standard error = 0.0, Computation time = 0.028806.
N=10^6
Ians = 0.0, Standard error = 0.0, Computation time = 0.282131099.
N=10^7
Ians = 0.3868309256925489, Standard error = 0.21207959450202113, Computation time = 2.7677205.
N=10^8
Ians = 1.7082305969777736, Standard error = 0.5570161065045919, Computation time = 28.5033501.
N=10^9
Ians = 1.2802405840617102, Standard error = 0.12926969998729035, Computation time = 292.426348399.
N=10^10
Ians = 1.4366314383717866, Standard error = 0.0493954801911099, Computation time = 2945.7095192.
N=10^11
Ians = 1.4516855424490653, Standard error = 0.015918731804291437, Computation time = 27563.9494109.
=#

# The end of the file.
