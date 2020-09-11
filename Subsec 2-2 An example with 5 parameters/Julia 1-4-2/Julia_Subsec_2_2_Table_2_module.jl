# The simple and antithetic Monte Carlo integrations for the Network 2v5p with unform distibution on the parameters.
# This julia file is written so that the simple and antithetic Monte Carlo functions from here being called in the other julia file for the parallelization.
# This file should be in the same directory as the Julia file Julia_Subsec_2_2_Table_2.jl

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
@fastmath function sumo_with_S(kk1,kk2,kk3,kk4,TT,NN) # Each kk is an ordered pair; a_i, b_i.
	S=0
	Ians=0
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
	return Ians,S
end
@fastmath function sumo_antithetic_with_S(kk1,kk2,kk3,kk4,TT,NN)
	S=0
	Ians=0
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
	return Ians,S
end

# The end of the module file.
