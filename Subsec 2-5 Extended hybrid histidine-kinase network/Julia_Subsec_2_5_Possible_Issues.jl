# Computations for the double HK network in the paper. In this file only parameters k12 and k14 are free. The rest are fixed.

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

import Distributions: Uniform
function samplor(a, b)
    return rand(Uniform(a, b))
end
function chi(a, b, x)
    if x < a
        return 0
    elseif x > b
        return 0
    else
        return 1
    end
end

# Here all parameters other than k12 and k14 are fixed to the vales given in equation (35) of the paper.

@fastmath function h2p12(t)
    return -468.4598789 * t^4 - 46854.40101 * t^3 - 1087.040656 * t^2 - 24572.83200 * t
end
@fastmath function g2p12(t, k14)
    return (-93770052422884376700 * t^5 - 187539354688350000 * t^4 * k14 + 1463950974630357138909 * t^4 - 204244474710835000000 * t^3 * k14 - 360657036215560929100 * t^3 - 133247984200000000000 * t^2 * k14 + 22142366768922618000000 * t^2 - 2860512000000000000000 * t * k14 + 206772287920000000000 * t + 11214585600000000000000) / (-468.4598789 * t^4 - 46854.40101 * t^3 - 1087.040656 * t^2 - 24572.83200 * t)
end
@fastmath function detJf2p(t, k12, k14)
    return -23442.51311 * t^4 + 4 * (-468.4598789 * k12 - 9.376967734 * k14 + 73197.54873) * t^3 + 3 * (-46854.40101 * k12 - 1021222374 * k14 - 18032.85181) * t^2 + 2 * (-1087.040656 * k12 - 6662.399210 * k14 + 1.107118338 * 10^6) * t - 24572.83200 * k12 - 143025.6000 * k14 + 10338.61440;
end
@fastmath function J2p(t, k14)
    return abs(detJf2p(t,g2p12(t, k14),k14))
end
@fastmath function IntegrandWithoutCOEFF2p(t, ak12, bk12, k14)
    return J2p(t, k14) * chi(ak12,bk12,g2p12(t,k14))/h2p12(t)
end
@fastmath function sumo_with_SE_t_2p(kk12,kk14,NN)
	S=0
	Ians=0
	st1=time_ns()
	COEFF=kk14[2]/(kk12[2]-kk12[1])
	t=samplor(0,kk14[2])
	k14=samplor(kk14[1],kk14[2])
	Ians+=IntegrandWithoutCOEFF2p(t,kk12[1],kk12[2],k14)
	Ians=COEFF*Ians
	@inbounds for i = 2:NN
		t=samplor(0,kk14[2])
		k14=samplor(kk14[1],kk14[2])
		delto=COEFF*IntegrandWithoutCOEFF2p(t,kk12[1],kk12[2],k14)-Ians
		Ians+=delto/i
		S+=delto^2*(i-1)/i
	end
	Estandardo=S/NN
	Estandardo=Estandardo/(NN-1)
	Estandardo=sqrt(Estandardo)
	st2=time_ns()
	return Ians,Estandardo,(st2-st1)/10^9
end

@fastmath function sumo_antithetic_with_SE_t_2p(kk12,kk14,NN)
	S=0
	Ians=0
	st1=time_ns()
	centero=[(0+kk14[2])/2,(kk14[1]+kk14[2])/2]
	COEFF=kk14[2]/(kk12[2]-kk12[1])
	t=samplor(0,kk14[2])
	k14=samplor(kk14[1],kk14[2])
	Ians+=IntegrandWithoutCOEFF2p(t,kk12[1],kk12[2],k14)
	Ians=COEFF*Ians
	t=2*centero[1]-t
	k14=2*centero[2]-k14
	delto=COEFF*IntegrandWithoutCOEFF2p(t,kk12[1],kk12[2],k14)-Ians
	Ians+=delto/2
	S+=delto^2/2
	@inbounds for i = 1:floor(NN/2)
		t=samplor(0,kk14[2])
		k14=samplor(kk14[1],kk14[2])
		delto=COEFF*IntegrandWithoutCOEFF2p(t,kk12[1],kk12[2],k14)-Ians
		Ians+=delto/(2*i+1)
		S+=(delto^2)*2*i/(2*i+1)
		t=2*centero[1]-t
		k14=2*centero[2]-k14
		delto=COEFF*IntegrandWithoutCOEFF2p(t,kk12[1],kk12[2],k14)-Ians
		Ians+=delto/(2*i+2)
		S+=(delto^2)*(2*i+1)/(2*i+2)
	end
	Estandardo=S/NN
	Estandardo=Estandardo/(NN-1)
	Estandardo=sqrt(Estandardo)
	st2=time_ns()
	return Ians,Estandardo,(st2-st1)/10^9
end

print("DoubleHK, 1v2p, x[5], k[12], k[14], B=[6.3,6.4]x[7.8,7.9]\n")

sumo_with_SE_t_2p([6.3,6.4],[7.8,7.9],100)
sumo_antithetic_with_SE_t_2p([6.3,6.4],[7.8,7.9],100)

for n = 1:8
	X1=sumo_with_SE_t_2p([6.3,6.4],[7.8,7.9],10^n)
	X2=sumo_antithetic_with_SE_t_2p([6.3,6.4],[7.8,7.9],10^n)
	print("N=",10^n,"\n")
	print("With simple Monte Carlo: \n",X1[1]," = Ians\n",X1[2]," = Standard","Error\n",X1[3]," = Time of computation\n")
	print("With antithetic Monte Carlo: \n",X2[1]," = Ians\n",X2[2]," =","Standard Error\n",X2[3]," = Time of computation\n")
	if (X2[3]*X2[2])!=0
		print("The efficiency of the antithetic to the simple method=",(X1[3]*X1[2])/(X2[3]*X2[2]),"\n")
	end
end

#= Output
DoubleHK, 1v2p, x[5], k[12], k[14], B=[6.3,6.4]x[7.8,7.9]
N=10
With simple Monte Carlo:
0.0 = Ians
0.0 = StandardError
2.2e-6 = Time of computation
With antithetic Monte Carlo:
0.0 = Ians
0.0 =Standard Error
4.9e-6 = Time of computation
N=100
With simple Monte Carlo:
0.0 = Ians
0.0 = StandardError
9.9e-6 = Time of computation
With antithetic Monte Carlo:
0.0 = Ians
0.0 =Standard Error
9.601e-6 = Time of computation
N=1000
With simple Monte Carlo:
0.0 = Ians
0.0 = StandardError
6.8299e-5 = Time of computation
With antithetic Monte Carlo:
0.0 = Ians
0.0 =Standard Error
6.87e-5 = Time of computation
N=10000
With simple Monte Carlo:
0.0 = Ians
0.0 = StandardError
0.000661001 = Time of computation
With antithetic Monte Carlo:
0.0 = Ians
0.0 =Standard Error
0.000660599 = Time of computation
N=100000
With simple Monte Carlo:
0.0 = Ians
0.0 = StandardError
0.007019101 = Time of computation
With antithetic Monte Carlo:
0.0 = Ians
0.0 =Standard Error
0.0066788 = Time of computation
N=1000000
With simple Monte Carlo:
0.0 = Ians
0.0 = StandardError
0.0680479 = Time of computation
With antithetic Monte Carlo:
0.0 = Ians
0.0 =Standard Error
0.0653157 = Time of computation
N=10000000
With simple Monte Carlo:
0.0 = Ians
0.0 = StandardError
0.647055201 = Time of computation
With antithetic Monte Carlo:
0.0 = Ians
0.0 =Standard Error
0.653433701 = Time of computation
N=100000000
With simple Monte Carlo:
0.0 = Ians
0.0 = StandardError
6.377053401 = Time of computation
With antithetic Monte Carlo:
0.0 = Ians
0.0 =Standard Error
6.3629198 = Time of computation
=#

# Here all parameters other than k12 and k14 are fixed and equal to 1.

@fastmath function h2p12(t)
    return 2*t^4 + 4*t^3 + 3*t^2 + t
end
@fastmath function g2p12(t, k14)
    return (4*t^5 - (-2*k14 - 4)*t^4 + 4*k14*t^3 - (-3*k14 + 4)*t^2 - (-k14 + 3)*t - 1) / (-2*t^4 - 4*t^3 - 3*t^2 - t)
end
@fastmath function detJf2p(t, k12, k14)
    return -20*t^4 + 4*(-2*k12 - 2*k14 - 4)*t^3 + 3*(-4*k12 - 4*k14)*t^2 + 2*(-3*k12 - 3*k14 + 4)*t - k12 - k14 + 3;
end
@fastmath function J2p(t, k14)
    return abs(detJf2p(t,g2p12(t, k14),k14))
end
@fastmath function IntegrandWithoutCOEFF2p(t, ak12, bk12, k14)
    return J2p(t, k14) * chi(ak12,bk12,g2p12(t,k14))/h2p12(t)
end
@fastmath function sumo_with_SE_t_2p(kk12,kk14,NN)
	S=0
	Ians=0
	st1=time_ns()
	COEFF=kk14[2]/(kk12[2]-kk12[1])
	t=samplor(0,kk14[2])
	k14=samplor(kk14[1],kk14[2])
	Ians+=IntegrandWithoutCOEFF2p(t,kk12[1],kk12[2],k14)
	Ians=COEFF*Ians
	@inbounds for i = 2:NN
		t=samplor(0,kk14[2])
		k14=samplor(kk14[1],kk14[2])
		delto=COEFF*IntegrandWithoutCOEFF2p(t,kk12[1],kk12[2],k14)-Ians
		Ians+=delto/i
		S+=delto^2*(i-1)/i
	end
	Estandardo=S/NN
	Estandardo=Estandardo/(NN-1)
	Estandardo=sqrt(Estandardo)
	st2=time_ns()
	return Ians,Estandardo,(st2-st1)/10^9
end

@fastmath function sumo_antithetic_with_SE_t_2p(kk12,kk14,NN)
	S=0
	Ians=0
	st1=time_ns()
	centero=[(0+kk14[2])/2,(kk14[1]+kk14[2])/2]
	COEFF=kk14[2]/(kk12[2]-kk12[1])
	t=samplor(0,kk14[2])
	k14=samplor(kk14[1],kk14[2])
	Ians+=IntegrandWithoutCOEFF2p(t,kk12[1],kk12[2],k14)
	Ians=COEFF*Ians
	t=2*centero[1]-t
	k14=2*centero[2]-k14
	delto=COEFF*IntegrandWithoutCOEFF2p(t,kk12[1],kk12[2],k14)-Ians
	Ians+=delto/2
	S+=delto^2/2
	@inbounds for i = 1:floor(NN/2)
		t=samplor(0,kk14[2])
		k14=samplor(kk14[1],kk14[2])
		delto=COEFF*IntegrandWithoutCOEFF2p(t,kk12[1],kk12[2],k14)-Ians
		Ians+=delto/(2*i+1)
		S+=(delto^2)*2*i/(2*i+1)
		t=2*centero[1]-t
		k14=2*centero[2]-k14
		delto=COEFF*IntegrandWithoutCOEFF2p(t,kk12[1],kk12[2],k14)-Ians
		Ians+=delto/(2*i+2)
		S+=(delto^2)*(2*i+1)/(2*i+2)
	end
	Estandardo=S/NN
	Estandardo=Estandardo/(NN-1)
	Estandardo=sqrt(Estandardo)
	st2=time_ns()
	return Ians,Estandardo,(st2-st1)/10^9
end

print("DoubleHK, 1v2p, x[5], k[12], k[14], B=[0,1]x[0,1]\n")

sumo_with_SE_t_2p([0.0,1,0],[0.0,1,0],100)
sumo_antithetic_with_SE_t_2p([0.0,1,0],[0.0,1,0],100)

for n = 1:8
	X1=sumo_with_SE_t_2p([0.0,1,0],[0.0,1,0],10^n)
	X2=sumo_antithetic_with_SE_t_2p([0.0,1,0],[0.0,1,0],10^n)
	print("N=",10^n,"\n")
	print("With simple Monte Carlo: \n",X1[1]," = Ians\n",X1[2]," = Standard","Error\n",X1[3]," = Time of computation\n")
	print("With antithetic Monte Carlo: \n",X2[1]," = Ians\n",X2[2]," =","Standard Error\n",X2[3]," = Time of computation\n")
	if (X2[3]*X2[2])!=0
		print("The efficiency of the antithetic to the simple method=",(X1[3]*X1[2])/(X2[3]*X2[2]),"\n")
	end
end

#= Output
DoubleHK, 1v2p, x[5], k[12], k[14], B=[0,1]x[0,1]
N=10
With simple Monte Carlo:
0.8188702465909725 = Ians
0.5474142657534327 = StandardError
1.501e-6 = Time of computation
With antithetic Monte Carlo:
1.166624532399677 = Ians
0.6195340720195684 =Standard Error
4.4e-6 = Time of computation
The efficiency of the antithetic to the simple method=0.3014247649254351
N=100
With simple Monte Carlo:
0.7408903294637547 = Ians
0.15361416590913227 = StandardError
6.301e-6 = Time of computation
With antithetic Monte Carlo:
0.9787503269968922 = Ians
0.17557786512690038 =Standard Error
8.8e-6 = Time of computation
The efficiency of the antithetic to the simple method=0.6264527361833737
N=1000
With simple Monte Carlo:
1.0221290022392455 = Ians
0.05171897357387131 = StandardError
6.12e-5 = Time of computation
With antithetic Monte Carlo:
0.9708488553857052 = Ians
0.051436705956723915 =Standard Error
5.8699e-5 = Time of computation
The efficiency of the antithetic to the simple method=1.0483286826570084
N=10000
With simple Monte Carlo:
0.9947666521502923 = Ians
0.016540890594517466 = StandardError
0.0006026 = Time of computation
With antithetic Monte Carlo:
0.9907279935736022 = Ians
0.01652289738822547 =Standard Error
0.000577 = Time of computation
The efficiency of the antithetic to the simple method=1.0455047193233158
N=100000
With simple Monte Carlo:
1.005963590562607 = Ians
0.005222391872369607 = StandardError
0.005846 = Time of computation
With antithetic Monte Carlo:
1.0018512391694216 = Ians
0.005209283655814571 =Standard Error
0.0056701 = Time of computation
The efficiency of the antithetic to the simple method=1.0336167611923879
N=1000000
With simple Monte Carlo:
0.9971364163316753 = Ians
0.001645765018550593 = StandardError
0.057759 = Time of computation
With antithetic Monte Carlo:
1.0014993784222017 = Ians
0.001648920082402065 =Standard Error
0.054200001 = Time of computation
The efficiency of the antithetic to the simple method=1.06362513209224
N=10000000
With simple Monte Carlo:
1.0007033816028228 = Ians
0.0005210941716519036 = StandardError
0.5604329 = Time of computation
With antithetic Monte Carlo:
1.0005540215317057 = Ians
0.0005211243389428338 =Standard Error
0.5284066 = Time of computation
The efficiency of the antithetic to the simple method=1.0605478000871031
N=100000000
With simple Monte Carlo:
1.0001509617674702 = Ians
0.00016476975123632455 = StandardError
5.6886389 = Time of computation
With antithetic Monte Carlo:
1.0002096804765566 = Ians
0.00016476384012355244 =Standard Error
5.3652026 = Time of computation
The efficiency of the antithetic to the simple method=1.0603221185293665
=#

# The end of the file.
