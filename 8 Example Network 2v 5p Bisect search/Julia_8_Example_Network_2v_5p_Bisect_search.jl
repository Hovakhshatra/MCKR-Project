# From Table 1 (a), we conclude that the suitables sample size for 2 digits of significant is 10^7 and the antithetic method takes 5/6 time of the Simple method.
# So we choose the Antithetic Monte Carlo with sample size 10^7 in all of the following integrations.

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
@fastmath function IntegrandWithoutCOEFF(x1,x2,k2,k3,k4,ak1,bk1,aT,bT)
	return J(x1,x2,k2,k3,k4)*chi(ak1,bk1,alok1(x1,x2,k2,k3,k4))*chi(aT,bT,aloT(x1,x2,k2,k3,k4))/(x2*x1^2)
end

@fastmath function Antithetic(B)
    # Antithetic Monte Carlo integration for the HK 1v 2p network with B as the parameter region and sample size equal to 10^7.
    # B is a list of 2 lists. Each sublist consists of two numbers.
    # The output is a numbers; I-hat.
    NN = 10^7
	kk1=B[1]
	kk2=B[2]
	kk3=B[3]
	kk4=B[4]
	TT=B[5]
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
	@inbounds for i = 1:floor(NN/2)
		x1=samplor(0,TT[2])
		x2=samplor(0,TT[2])
		k2=samplor(kk2[1],kk2[2])
		k3=samplor(kk3[1],kk3[2])
		k4=samplor(kk4[1],kk4[2])
		delto=COEFF*IntegrandWithoutCOEFF(x1,x2,k2,k3,k4,kk1[1],kk1[2],TT[1],TT[2])-Ians
		Ians+=delto/(2*i+1)
		x1=2*centero[1]-x1
		x2=2*centero[2]-x2
		k2=2*centero[3]-k2
		k3=2*centero[4]-k3
		k4=2*centero[5]-k4
		delto=COEFF*IntegrandWithoutCOEFF(x1,x2,k2,k3,k4,kk1[1],kk1[2],TT[1],TT[2])-Ians
		Ians+=delto/(2*i+2)
	end
   	return Ians
end

function SameBox(B) # To copy a box without dependency between the copy and the original.
    return(deepcopy(B))
end

@fastmath function BisectSearch(fun, B, minimum_size,sh_terminate)
    # The bisect search for finding a subbox with value of "fun" equal to 3.
    # This function should be used when M_min=1 and M_max=3.
    # fun is the function which associate a number to each subbox, B is the initial box, minimum_size is the minimum size condition for edges of the subboxes, sh_terminate is the upper bound for number of bisect steps.
    # fun is a function, B is a list of list of numbers, minimum_size is a list of numbers, sh_terminater is an integer.
    # The output is a 3-tuple or 4-tuple depending on finding a subbox inside the multistationary region.
    # A string message saying it found multistationary subbox, found no multistationary subbox, needs more bisect steps.
	# Number of times function "fun" is evaluated at a subbox.
    # The amount of time the computation took.
	# In case of finding a subbox inside the multistationary region, a list consisting of the bisect step, the subbox and the value of "fun" on this subbox.
    # The termination condition is determined by number of bisect steps, all edges being smaller than minimum sizes, or the current chosen subbox having value of fun greater than 2.95 or less than 1.05.
    st1=time_ns()
    d = length(B)
	shlist=[]
	for i=1:d
		push!(shlist,i)
	end
	for i=1:d-1
		push!(shlist,i)
	end
	sh1=0
    sh2=0
    ans0 = fun(B)
    integral_number=1
    if ans0 < 1.05
        return "No multistationary subbox",integral_number,(time_ns()-st1)/(10^9)
    elseif ans0 > 2.95
        return "Found Multistationary subbox",integral_number,(time_ns()-st1)/(10^9),[sh2,B,ans0]
    end
    @inbounds while sh2 != sh_terminate
        if sh1 == d
            sh1 = 1
        else
            sh1 += 1
        end
		sizeend=1
		for i=sh1:sh1+d-1
			if B[shlist[i]][2]-B[shlist[i]][1] > minimum_size[shlist[i]]
				sizeend=0
		    	sh2 += 1
				B1 = SameBox(B)
		        B1[shlist[i]][2] = (B[shlist[i]][1] + B[shlist[i]][2]) / 2
			    ans1=fun(B1)
			    B2 = SameBox(B)
		        B2[shlist[i]][1] = (B[shlist[i]][1] + B[shlist[i]][2]) / 2
			    ans2=fun(B2)
				integral_number+=2
				if ans1 > ans2
					B=B1
					ans0=ans1
				else
					B=B2
					ans0=ans2
				end
				if ans0 < 1.05
			        return "No multistationary subbox",integral_number,(time_ns()-st1)/(10^9)
			    elseif ans0 > 2.95
			        return "Found Multistationary subbox",integral_number,(time_ns()-st1)/(10^9),[sh2,B,ans0]
			    end
				break
			end
		end
		if sizeend==1
			return "Termination condition on size of box occured",integral_number,(time_ns()-st1)/(10^9)
		end
    end
    st2=time_ns()
    return "Termination condition on the number of bisect steps occured",integral_number,(time_ns()-st1)/(10^9)
end

Binput=[[0.0,100.0],[0.0,2.0],[0.0,200.0],[0.0,100.0],[0.0,2.0]]
bisectlimit=10
minimumsizeinput=[10.0,0.1,10.0,10.0,0.1]
print("\nThe bisect search for Network 2v 5p with initial box\n$Binput\nand the constrain on the minimum size of each interval\n$minimumsizeinput\nand a limit on the maximum number of bisect steps $bisectlimit.\n")
X1=BisectSearch(Antithetic,Binput,minimumsizeinput,bisectlimit)
print("\nIt took $(X1[3]) seconds and computed $(X1[2]) integrals to conclude that\n$(X1[1])\n")
if X1[1]=="Found Multistationary subbox"
	print("The following subbox is almost inside the multistationary region.\n$(X1[4][2])\nWith I-hat equal to $(X1[4][3])\n")
end

# The output
#=
The bisect search for Network 2v 5p with initial box
[[0.0, 100.0], [0.0, 2.0], [0.0, 200.0], [0.0, 100.0], [0.0, 2.0]]
and the constrain on the minimum size of each interval
[10.0, 0.1, 10.0, 10.0, 0.1]
and a limit on the maximum number of bisect steps 10.

It took 14.8843844 seconds and computed 19 integrals to conclude that
Found Multistationary subbox
The following subbox is almost inside the multistationary region.
[[50.0, 75.0], [1.0, 1.5], [150.0, 200.0], [0.0, 25.0], [1.0, 2.0]]
With I-hat equal to 2.9726942915759573
=#

# The end of the file.
