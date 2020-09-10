# The Monte Carlo integrations used in Table 1 of the paper.

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

@fastmath function alok1(x1,x2,k2,k3,k4)
	return (2*k2*x1^3+k3*x1*x2^2-2*k4*x2^3)/(x2*x1^2)
end
@fastmath function aloT(x1,x2,k2,k3,k4)
	return x1+x2
end
@fastmath function detJf(x1,x2,k1,k2,k3,k4)
	return -(6*k2+k1)*x1^2-(6*k4+k3)*x2^2+2*(k1+k3)*x1*x2
end
@fastmath function J(x1,x2,k2,k3,k4)
	return abs(detJf(x1,x2,alok1(x1,x2,k2,k3,k4),k2,k3,k4))
end

# First the computations for part (a) Uniform distributions

@fastmath function IntegrandWithoutCOEFF(x1,x2,k2,k3,k4,ak1,bk1,aT,bT)
	return J(x1,x2,k2,k3,k4)*chi(ak1,bk1,alok1(x1,x2,k2,k3,k4))*chi(aT,bT,aloT(x1,x2,k2,k3,k4))/(x2*x1^2)
end
@fastmath function sumo_with_SE_t(kk1,kk2,kk3,kk4,TT,NN) # Each kk is an ordered pair; a_i, b_i.
	S=0
	Ians=0
	t1=time_ns()
	COEFF=(TT[2]^2)/((TT[2]-TT[1])*(kk1[2]-kk1[1]))
	x1=samplor(0,TT[2])
	x2=samplor(0,TT[2])
	k2=samplor(kk2[1],kk2[2])
	k3=samplor(kk3[1],kk3[2])
	k4=samplor(kk4[1],kk4[2])
	Ians+=IntegrandWithoutCOEFF(x1,x2,k2,k3,k4,kk1[1],kk1[2],TT[1],TT[2])
	Ians=COEFF*Ians
	@inbounds for i = 2:NN
		x1=samplor(0,TT[2])
		x2=samplor(0,TT[2])
		k2=samplor(kk2[1],kk2[2])
		k3=samplor(kk3[1],kk3[2])
		k4=samplor(kk4[1],kk4[2])
		delto=COEFF*IntegrandWithoutCOEFF(x1,x2,k2,k3,k4,kk1[1],kk1[2],TT[1],TT[2])-Ians
		Ians+=delto/(i+1)
		S+=(delto^2)*i/(i+1)
	end
	Estandardo=S/NN
	Estandardo=Estandardo/(NN-1)
	Estandardo=sqrt(Estandardo)
	t2=time_ns()
	return Ians,Estandardo,(t2-t1)/10^9
end
@fastmath function sumo_antithetic_with_SE_t(kk1,kk2,kk3,kk4,TT,NN)
	S=0
	Ians=0
	t1=time_ns()
	centero=[(0+TT[2])/2,(0+TT[2])/2,(kk2[1]+kk2[2])/2,(kk3[1]+kk3[2])/2,(kk4[1]+kk4[2])/2]
	COEFF=(TT[2]^2)/((TT[2]-TT[1])*(kk1[2]-kk1[1]))
	x1=samplor(0,TT[2])
	x2=samplor(0,TT[2])
	k2=samplor(kk2[1],kk2[2])
	k3=samplor(kk3[1],kk3[2])
	k4=samplor(kk4[1],kk4[2])
	Ians+=IntegrandWithoutCOEFF(x1,x2,k2,k3,k4,kk1[1],kk1[2],TT[1],TT[2])
	Ians=COEFF*Ians
	x1=2*centero[1]-x1
	x2=2*centero[2]-x2
	k2=2*centero[3]-k2
	k3=2*centero[4]-k3
	k4=2*centero[5]-k4
	delto=COEFF*IntegrandWithoutCOEFF(x1,x2,k2,k3,k4,kk1[1],kk1[2],TT[1],TT[2])-Ians
	Ians+=delto/2
	S+=(delto^2)*1/2
	@inbounds for i = 1:floor(NN/2)
		x1=samplor(0,TT[2])
		x2=samplor(0,TT[2])
		k2=samplor(kk2[1],kk2[2])
		k3=samplor(kk3[1],kk3[2])
		k4=samplor(kk4[1],kk4[2])
		delto=COEFF*IntegrandWithoutCOEFF(x1,x2,k2,k3,k4,kk1[1],kk1[2],TT[1],TT[2])-Ians
		Ians+=delto/(2*i+1)
		S+=(delto^2)*2*i/(2*i+1)
		x1=2*centero[1]-x1
		x2=2*centero[2]-x2
		k2=2*centero[3]-k2
		k3=2*centero[4]-k3
		k4=2*centero[5]-k4
		delto=COEFF*IntegrandWithoutCOEFF(x1,x2,k2,k3,k4,kk1[1],kk1[2],TT[1],TT[2])-Ians
		Ians+=delto/(2*i+2)
		S+=(delto^2)*(2*i+1)/(2*i+2)
	end
	Estandardo=S/NN
	Estandardo=Estandardo/(NN-1)
	Estandardo=sqrt(Estandardo)
	t2=time_ns()
	return Ians,Estandardo,(t2-t1)/10^9
end


for n = 1:9
	X1=sumo_with_SE_t([0,100],[0,2],[0,200],[0,100],[0,2],10^n)
	X2=sumo_antithetic_with_SE_t([0,100],[0,2],[0,200],[0,100],[0,2],10^n)
	print("N=",10^n,"\n")
	print("With simple Monte Carlo: \n",X1[1]," = Ians\n",X1[2]," = Standard","Error\n",X1[3]," = Time of computation\n")
	print("With antithetic Monte Carlo: \n",X2[1]," = Ians\n",X2[2]," =","Standard Error\n",X2[3]," = Time of computation\n")
	if (X2[3]*X2[2])!=0
		print("The efficiency of the antithetic to the simple method=",(X1[3]*X1[2])/(X2[3]*X2[2]),"\n")
	end
end

# The output
#=
N=10
With simple Monte Carlo:
0.527070881255435 = Ians
0.5826987836820811 = StandardError
2.2899e-5 = Time of computation
With antithetic Monte Carlo:
0.1994769270422469 = Ians
0.1716493769623022 =Standard Error
8.001e-6 = Time of computation
The efficiency of the antithetic to the simple method=9.715700684204101
N=100
With simple Monte Carlo:
2.542496024686625 = Ians
1.8836177917318668 = StandardError
2.0801e-5 = Time of computation
With antithetic Monte Carlo:
1.099524942885649 = Ians
0.41486771322711696 =Standard Error
1.48e-5 = Time of computation
The efficiency of the antithetic to the simple method=6.381248319125964
N=1000
With simple Monte Carlo:
2.470114858573769 = Ians
0.744916211380292 = StandardError
0.000165101 = Time of computation
With antithetic Monte Carlo:
0.9423223445789625 = Ians
0.14734894066697932 =Standard Error
9.4999e-5 = Time of computation
The efficiency of the antithetic to the simple method=8.785997550820658
N=10000
With simple Monte Carlo:
1.4678351234404488 = Ians
0.17541589805824356 = StandardError
0.001570899 = Time of computation
With antithetic Monte Carlo:
1.662102537905749 = Ians
0.3327997097674981 =Standard Error
0.0009501 = Time of computation
The efficiency of the antithetic to the simple method=0.871495115767444
N=100000
With simple Monte Carlo:
1.9898912973529528 = Ians
0.5953702420779486 = StandardError
0.015469 = Time of computation
With antithetic Monte Carlo:
1.392013810905736 = Ians
0.05525623997994873 =Standard Error
0.008919999 = Time of computation
The efficiency of the antithetic to the simple method=18.685435676064067
N=1000000
With simple Monte Carlo:
1.4321610174913961 = Ians
0.030766932466041595 = StandardError
0.1506252 = Time of computation
With antithetic Monte Carlo:
1.4490360503534478 = Ians
0.03357545713593764 =Standard Error
0.090292401 = Time of computation
The efficiency of the antithetic to the simple method=1.5286522461172367
N=10000000
With simple Monte Carlo:
1.4216616899614227 = Ians
0.007325828328597585 = StandardError
1.5362033 = Time of computation
With antithetic Monte Carlo:
1.419472226036245 = Ians
0.008804104156489781 =Standard Error
0.939077999 = Time of computation
The efficiency of the antithetic to the simple method=1.3611895092993156
N=100000000
With simple Monte Carlo:
1.413392972223649 = Ians
0.0029837457963567556 = StandardError
15.494148199 = Time of computation
With antithetic Monte Carlo:
1.4151109371603663 = Ians
0.0027409518159804093 =Standard Error
10.040993001 = Time of computation
The efficiency of the antithetic to the simple method=1.6797763421641456
N=1000000000
With simple Monte Carlo:
1.418665401538655 = Ians
0.0014011376894688375 = StandardError
155.544266099 = Time of computation
With antithetic Monte Carlo:
1.4182581455624133 = Ians
0.0012250901172471875 =Standard Error
92.3336835 = Time of computation
The efficiency of the antithetic to the simple method=1.926666844237346
=#

# Now the computations for part (b) Truncated Normal Distributions

# We need to define two functions TNpdf and samplorTN in below.
using Distributions
function samplor(a, b)
   	return rand(Uniform(a, b))
end
function TNpdf(a, b, m, v, x)
   	return(Distributions.pdf(Distributions.Truncated(Normal(m, v), a, b), x))
end
function samplorTN(a, b, m, v)
   	return(rand(Distributions.Truncated(Normal(m, v), a, b)))
end

@fastmath function IntegrandWithoutCOEFF(x1, x2, k2, k3, k4, ak1, bk1, mk1, vk1, aT, bT, mT, vT)
   	return J(x1, x2, k2, k3, k4) * TNpdf(ak1, bk1, mk1, vk1, alok1(x1, x2, k2, k3, k4)) * TNpdf(aT, bT, mT, vT, aloT(x1, x2, k2, k3, k4)) / (x2 * x1^2)
end
@fastmath function sumo_with_SE_t(kk1, kk2, kk3, kk4, TT, NN) # Each kk here is a 4-tuple; a_i, b_i, mu_i, sigma_i^2.
   	S = 0
   	Ians = 0
   	t1 = time_ns()
   	COEFF = (TT[2]^2)
   	x1 = samplor(0, TT[2])
   	x2 = samplor(0, TT[2])
   	k2 = samplorTN(kk2[1], kk2[2], kk2[3], kk2[4])
   	k3 = samplorTN(kk3[1], kk3[2], kk3[3], kk3[4])
   	k4 = samplorTN(kk4[1], kk4[2], kk4[3], kk4[4])
   	Ians += IntegrandWithoutCOEFF(x1, x2, k2, k3, k4, kk1[1], kk1[2], kk1[3], kk1[4], TT[1], TT[2], TT[3], TT[4])
   	Ians = COEFF * Ians
   	@inbounds for i = 2:NN
      		x1 = samplor(0, TT[2])
      		x2 = samplor(0, TT[2])
      		k2 = samplorTN(kk2[1], kk2[2], kk2[3], kk2[4])
      		k3 = samplorTN(kk3[1], kk3[2], kk3[3], kk3[4])
      		k4 = samplorTN(kk4[1], kk4[2], kk4[3], kk4[4])
      		delto = COEFF * IntegrandWithoutCOEFF(x1, x2, k2, k3, k4, kk1[1], kk1[2], kk1[3], kk1[4], TT[1], TT[2], TT[3], TT[4]) - Ians
      		Ians += delto / (i + 1)
      		S += (delto^2) * i / (i + 1)
   	end
   	Estandardo = S / NN
   	Estandardo = Estandardo / (NN - 1)
   	Estandardo = sqrt(Estandardo)
   	t2 = time_ns()
   	return Ians, Estandardo, (t2 - t1) / 10^9
end
@fastmath function sumo_antithetic_with_SE_t(kk1, kk2, kk3, kk4, TT, NN)
   	S = 0
   	Ians = 0
   	t1 = time_ns()
   	if (kk2[1] + kk2[2]) / 2 != kk2[3] || (kk3[1] + kk3[2]) / 2 != kk3[3] || (kk4[1] + kk4[2]) / 2 != kk4[3]
      		return error # Because for antithetic Monte-Carlo we need symmetric distribtuion (with respect to its support, not only the function alone).
   	end
   	centero = [(0 + TT[2]) / 2,(0 + TT[2]) / 2,kk2[3],kk3[3],kk4[3]]
   	COEFF = (TT[2]^2)
   	x1 = samplor(0, TT[2])
   	x2 = samplor(0, TT[2])
   	k2 = samplorTN(kk2[1], kk2[2], kk2[3], kk2[4])
   	k3 = samplorTN(kk3[1], kk3[2], kk3[3], kk3[4])
   	k4 = samplorTN(kk4[1], kk4[2], kk4[3], kk4[4])
   	Ians += IntegrandWithoutCOEFF(x1, x2, k2, k3, k4, kk1[1], kk1[2], kk1[3], kk1[4], TT[1], TT[2], TT[3], TT[4])
   	Ians = COEFF * Ians
   	x1 = 2 * centero[1] - x1
   	x2 = 2 * centero[2] - x2
   	k2 = 2 * centero[3] - k2
   	k3 = 2 * centero[4] - k3
   	k4 = 2 * centero[5] - k4
   	delto = COEFF * IntegrandWithoutCOEFF(x1, x2, k2, k3, k4, kk1[1], kk1[2], kk1[3], kk1[4], TT[1], TT[2], TT[3], TT[4]) - Ians
   	Ians += delto / 2
   	S += (delto^2) * 1 / 2
   	@inbounds for i = 1:floor(NN / 2)
      		x1 = samplor(0, TT[2])
      		x2 = samplor(0, TT[2])
      		k2 = samplorTN(kk2[1], kk2[2], kk2[3], kk2[4])
      		k3 = samplorTN(kk3[1], kk3[2], kk3[3], kk3[4])
      		k4 = samplorTN(kk4[1], kk4[2], kk4[3], kk4[4])
      		delto = COEFF * IntegrandWithoutCOEFF(x1, x2, k2, k3, k4, kk1[1], kk1[2], kk1[3], kk1[4], TT[1], TT[2], TT[3], TT[4]) - Ians
      		Ians += delto / (2 * i + 1)
      		S += (delto^2) * 2 * i / (2 * i + 1)
      		x1 = 2 * centero[1] - x1
      		x2 = 2 * centero[2] - x2
      		k2 = 2 * centero[3] - k2
      		k3 = 2 * centero[4] - k3
      		k4 = 2 * centero[5] - k4
      		delto = COEFF * IntegrandWithoutCOEFF(x1, x2, k2, k3, k4, kk1[1], kk1[2], kk1[3], kk1[4], TT[1], TT[2], TT[3], TT[4]) - Ians
      		Ians += delto / (2 * i + 2)
      		S += (delto^2) * (2 * i + 1) / (2 * i + 2)
   	end
   	Estandardo = S / NN
   	Estandardo = Estandardo / (NN - 1)
   	Estandardo = sqrt(Estandardo)
   	t2 = time_ns()
   	return Ians, Estandardo, (t2 - t1) / 10^9
end
#

for n = 1:9
   	X1 = sumo_with_SE_t([0,100,50,1], [0,2,1,0.01], [0,200,100,1], [0,100,50,1], [0,2,1,0.01], 10^n)
   	X2 = sumo_antithetic_with_SE_t([0,100,50,1], [0,2,1,0.01], [0,200,100,1], [0,100,50,1], [0,2,1,0.01], 10^n)
   	print("N=", 10^n, "\n")
   	print("With simple Monte Carlo: \n", X1[1], " = Ians\n", X1[2], " = Standard", "Error\n", X1[3], " = Time of computation\n")
   	print("With antithetic Monte Carlo: \n", X2[1], " = Ians\n", X2[2], " =", "Standard Error\n", X2[3], " = Time of computation\n")
   	if (X2[3] * X2[2]) != 0
      		print("The efficiency of the antithetic to the simple method=", (X1[3] * X1[2]) / (X2[3] * X2[2]), "\n")
   	end
end

# The output
#=
N=10
With simple Monte Carlo:
2.1961157682347245e-136 = Ians
2.4278973330890967e-136 = StandardError
1.59e-5 = Time of computation
With antithetic Monte Carlo:
0.0 = Ians
0.0 =Standard Error
2.3201e-5 = Time of computation
N=100
With simple Monte Carlo:
1.9467579201407053e-80 = Ians
1.9663238179910386e-80 = StandardError
6.02e-5 = Time of computation
With antithetic Monte Carlo:
2.8325060689863976e-105 = Ians
2.8894422879249166e-105 =Standard Error
5.73e-5 = Time of computation
The efficiency of the antithetic to the simple method=7.14961876178925e24
N=1000
With simple Monte Carlo:
1.041302972936955e-23 = Ians
1.042344797082421e-23 = StandardError
0.000542801 = Time of computation
With antithetic Monte Carlo:
4.107059637836882e-9 = Ians
4.115277868281302e-9 =Standard Error
0.000391101 = Time of computation
The efficiency of the antithetic to the simple method=3.515312601768633e-15
N=10000
With simple Monte Carlo:
1.1023945490046478 = Ians
1.1010910189271206 = StandardError
0.0055866 = Time of computation
With antithetic Monte Carlo:
0.0021414127502171042 = Ians
0.0020972653391067323 =Standard Error
0.004718499 = Time of computation
The efficiency of the antithetic to the simple method=621.603642299938
N=100000
With simple Monte Carlo:
0.127164054554445 = Ians
0.06905929664847914 = StandardError
0.0546727 = Time of computation
With antithetic Monte Carlo:
2.0267111175175963 = Ians
0.8170046601520861 =Standard Error
0.0380528 = Time of computation
The efficiency of the antithetic to the simple method=0.12144552609639689
N=1000000
With simple Monte Carlo:
1.0190028920329663 = Ians
0.17792830375521115 = StandardError
0.520236099 = Time of computation
With antithetic Monte Carlo:
1.0210253983591764 = Ians
0.17132686876815298 =Standard Error
0.3727226 = Time of computation
The efficiency of the antithetic to the simple method=1.4495537318881115
N=10000000
With simple Monte Carlo:
0.9630442062629708 = Ians
0.055539584867976326 = StandardError
5.2181101 = Time of computation
With antithetic Monte Carlo:
0.9645204050548902 = Ians
0.05688799256503448 =Standard Error
3.681650801 = Time of computation
The efficiency of the antithetic to the simple method=1.3837341674453802
N=100000000
With simple Monte Carlo:
1.019522590247999 = Ians
0.018544074771678576 = StandardError
52.1472254 = Time of computation
With antithetic Monte Carlo:
1.033604754848135 = Ians
0.01861847521801407 =Standard Error
37.9521037 = Time of computation
The efficiency of the antithetic to the simple method=1.3685365808123355
N=1000000000
With simple Monte Carlo:
1.0107017602003563 = Ians
0.005810387622929225 = StandardError
537.260416701 = Time of computation
With antithetic Monte Carlo:
0.9886805877030328 = Ians
0.005727196093538311 =Standard Error
371.116398899 = Time of computation
The efficiency of the antithetic to the simple method=1.4687157497766636
=#

# The end of the file.
