# The computations for the HK 1v 8p of the paper. This example is used at the end of subsection 2.1.

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

@fastmath function h1(t, k1, k2, k3, k4, k5, k6)
   	return(k1 * k2 * k4 * k5 * t^2 + k1 * k2 * k3 * k5 * t)
end
@fastmath function beta1(t, k1, k2, k3, k4, k5, k6)
   	return((k4 * k5 * (k1 + k2) * t^2 + k1 * k5 * (k2 + k3) * t + k1 * k2 * k3) * k6)
end
@fastmath function g1(t, T2, k1, k2, k3, k4, k5, k6)
   	return beta1(t, k1, k2, k3, k4, k5, k6) * (T2 - t) / h1(t, k1, k2, k3, k4, k5, k6)
end
@fastmath function detJf(t, T1, T2, k1, k2, k3, k4, k5, k6)
   	return(T1 * k1 * k2 * k3 * k5 - T2 * k1 * k2 * k5 * k6 - T2 * k1 * k3 * k5 * k6 + k1 * k2 * k3 * k6 + 3 * t^2 * (k1 * k4 * k5 * k6
	       + k2 * k4 * k5 * k6) + 2 * t * (T1 * k1 * k2 * k4 * k5 - T2 * k1 * k4 * k5 * k6 - T2 * k2 * k4 * k5 * k6 + k1 * k2 * k5 * k6 + k1 * k3 * k5 * k6))
end
@fastmath function J(t, T2, k1, k2, k3, k4, k5, k6)
   	return abs(detJf(t, g1(t, T2, k1, k2, k3, k4, k5, k6), T2, k1, k2, k3, k4, k5, k6))
end
@fastmath function IntegrandWithoutCOEFF(t, aT1, bT1, T2, k1, k2, k3, k4, k5, k6)
   	return J(t, T2, k1, k2, k3, k4, k5, k6) * chi(aT1, bT1, g1(t, T2, k1, k2, k3, k4, k5, k6)) / h1(t, k1, k2, k3, k4, k5, k6)
end

@fastmath function sumo_with_SE_t(kk1, kk2, kk3, kk4, kk5, kk6, TT1, TT2, NN)
   	S = 0
   	Ians = 0
   	st1 = time_ns()
   	COEFF = TT2[2] / (TT1[2] - TT1[1])
   	t = samplor(0, TT2[2])
   	k1 = samplor(kk1[1], kk1[2])
   	k2 = samplor(kk2[1], kk2[2])
   	k3 = samplor(kk3[1], kk3[2])
   	k4 = samplor(kk4[1], kk4[2])
   	k5 = samplor(kk5[1], kk5[2])
   	k6 = samplor(kk6[1], kk6[2])
   	T2 = samplor(TT2[1], TT2[2])
   	Ians += IntegrandWithoutCOEFF(t, TT1[1], TT1[2], T2, k1, k2, k3, k4, k5, k6)
   	Ians = COEFF * Ians
   	@inbounds for i = 2:NN
      		t = samplor(0, TT2[2])
      		k1 = samplor(kk1[1], kk1[2])
      		k2 = samplor(kk2[1], kk2[2])
      		k3 = samplor(kk3[1], kk3[2])
      		k4 = samplor(kk4[1], kk4[2])
      		k5 = samplor(kk5[1], kk5[2])
      		k6 = samplor(kk6[1], kk6[2])
      		T2 = samplor(TT2[1], TT2[2])
      		delto = COEFF * IntegrandWithoutCOEFF(t, TT1[1], TT1[2], T2, k1, k2, k3, k4, k5, k6) - Ians
      		Ians += delto / i
      		S += delto^2 * (i - 1) / i
   	end
   	Estandardo = S / NN
   	Estandardo = Estandardo / (NN - 1)
   	Estandardo = sqrt(Estandardo)
   	st2 = time_ns()
   	return Ians, Estandardo, (st2 - st1) / 10^9
end

@fastmath function sumo_antithetic_with_SE_t(kk1, kk2, kk3, kk4, kk5, kk6, TT1, TT2, NN)
   	S = 0
   	Ians = 0
   	st1 = time_ns()
   	centero = [(0 + TT2[2]) / 2,(kk1[1] + kk1[2]) / 2,(kk2[1] + kk2[2]) / 2,(kk3[1] + kk3[2]) / 2,(kk4[1] + kk4[2]) / 2,(kk5[1] + kk5[2]) / 2,(kk6[1] + kk6[2]) / 2
	        ,(TT2[1] + TT2[2]) / 2]
   	COEFF = TT2[2] / (TT1[2] - TT1[1])
   	t = samplor(0, TT2[2])
   	k1 = samplor(kk1[1], kk1[2])
   	k2 = samplor(kk2[1], kk2[2])
   	k3 = samplor(kk3[1], kk3[2])
   	k4 = samplor(kk4[1], kk4[2])
   	k5 = samplor(kk5[1], kk5[2])
   	k6 = samplor(kk6[1], kk6[2])
   	T2 = samplor(TT2[1], TT2[2])
   	Ians += IntegrandWithoutCOEFF(t, TT1[1], TT1[2], T2, k1, k2, k3, k4, k5, k6)
   	Ians = COEFF * Ians
   	t = 2 * centero[1] - t
   	k1 = 2 * centero[2] - k1
   	k2 = 2 * centero[3] - k2
   	k3 = 2 * centero[4] - k3
   	k4 = 2 * centero[5] - k4
   	k5 = 2 * centero[6] - k5
   	k6 = 2 * centero[7] - k6
   	T2 = 2 * centero[8] - T2
   	delto = COEFF * IntegrandWithoutCOEFF(t, TT1[1], TT1[2], T2, k1, k2, k3, k4, k5, k6) - Ians
   	Ians += delto / 2
   	S += delto^2 / 2
   	@inbounds for i = 1:floor(NN / 2)
      		t = samplor(0, TT2[2])
      		k1 = samplor(kk1[1], kk1[2])
      		k2 = samplor(kk2[1], kk2[2])
      		k3 = samplor(kk3[1], kk3[2])
      		k4 = samplor(kk4[1], kk4[2])
      		k5 = samplor(kk5[1], kk5[2])
      		k6 = samplor(kk6[1], kk6[2])
      		T2 = samplor(TT2[1], TT2[2])
      		delto = COEFF * IntegrandWithoutCOEFF(t, TT1[1], TT1[2], T2, k1, k2, k3, k4, k5, k6) - Ians
      		Ians += delto / (2 * i + 1)
      		S += (delto^2) * 2 * i / (2 * i + 1)
      		t = 2 * centero[1] - t
      		k1 = 2 * centero[2] - k1
      		k2 = 2 * centero[3] - k2
      		k3 = 2 * centero[4] - k3
      		k4 = 2 * centero[5] - k4
      		k5 = 2 * centero[6] - k5
      		k6 = 2 * centero[7] - k6
      		T2 = 2 * centero[8] - T2
      		delto = COEFF * IntegrandWithoutCOEFF(t, TT1[1], TT1[2], T2, k1, k2, k3, k4, k5, k6) - Ians
      		Ians += delto / (2 * i + 2)
      		S += (delto^2) * (2 * i + 1) / (2 * i + 2)
   	end
   	Estandardo = S / NN
   	Estandardo = Estandardo / (NN - 1)
   	Estandardo = sqrt(Estandardo)
   	st2 = time_ns()
   	return Ians, Estandardo, (st2 - st1) / 10^9
end

# To find the suitable minimum sample size for 2 digits of significant.
#
for n = 1:7
   	X1 = sumo_with_SE_t([0,1], [0,200], [0,100], [0,100], [0,200], [0,10], [0,5], [0,5], 10^n)
   	X2 = sumo_antithetic_with_SE_t([0,1], [0,200], [0,100], [0,100], [0,200], [0,10], [0,5], [0,5], 10^n)
   	print("N=", 10^n, "\n")
   	print("With simple Monte Carlo: \n", X1[1], " = Ians\n", X1[2], " = Standard", "Error\n", X1[3], " = Time of computation\n")
   	print("With antithetic Monte Carlo: \n", X2[1], " = Ians\n", X2[2], " =", "Standard Error\n", X2[3], " = Time of computation\n")
   	if (X2[3] * X2[2]) != 0
      		print("The efficiency of the antithetic to the simple method=", (X1[3] * X1[2]) / (X2[3] * X2[2]), "\n")
   	end
end

# The output.
#=
N=10
With simple Monte Carlo:
0.1954108972008636 = Ians
0.19541089720086358 = StandardError
9.2e-6 = Time of computation
With antithetic Monte Carlo:
0.5650376359322211 = Ians
0.38523630068540066 =Standard Error
9.199e-6 = Time of computation
The efficiency of the antithetic to the simple method=0.5073045802720323
N=100
With simple Monte Carlo:
0.9498514174653341 = Ians
0.5836710188820787 = StandardError
2.8e-5 = Time of computation
With antithetic Monte Carlo:
1.6333827936868879 = Ians
1.1430041673698372 =Standard Error
2.1299e-5 = Time of computation
The efficiency of the antithetic to the simple method=0.671303843502122
N=1000
With simple Monte Carlo:
1.745736119322389 = Ians
0.5152339299322087 = StandardError
0.000249701 = Time of computation
With antithetic Monte Carlo:
0.7615845767474352 = Ians
0.08258096139488895 =Standard Error
0.0001426 = Time of computation
The efficiency of the antithetic to the simple method=10.925096371582214
N=10000
With simple Monte Carlo:
1.532046751566549 = Ians
0.29618568164617926 = StandardError
0.0024938 = Time of computation
With antithetic Monte Carlo:
1.1849084220735275 = Ians
0.18559954858011854 =Standard Error
0.0014095 = Time of computation
The efficiency of the antithetic to the simple method=2.8234732961922817
N=100000
With simple Monte Carlo:
1.1033651460190443 = Ians
0.03191546624303278 = StandardError
0.024954299 = Time of computation
With antithetic Monte Carlo:
1.1540643642467294 = Ians
0.05147506710727833 =Standard Error
0.013797801 = Time of computation
The efficiency of the antithetic to the simple method=1.1213463383891138
N=1000000
With simple Monte Carlo:
1.1916698979398175 = Ians
0.017609697817659253 = StandardError
0.2417435 = Time of computation
With antithetic Monte Carlo:
1.189237875669009 = Ians
0.014341282407222914 =Standard Error
0.1384544 = Time of computation
The efficiency of the antithetic to the simple method=2.143936703096952
N=10000000
With simple Monte Carlo:
1.204680866055838 = Ians
0.0122599168120903 = StandardError
2.405281401 = Time of computation
With antithetic Monte Carlo:
1.1994527938616892 = Ians
0.006783265430812669 =Standard Error
1.3649523 = Time of computation
The efficiency of the antithetic to the simple method=3.184909749761219
=#

# The conclusion is that 10^7 is a good sample size.

@fastmath function Antithetic(B)
   	NN = 10^7
   	kk1 = B[1]
   	kk2 = B[2]
   	kk3 = B[3]
   	kk4 = B[4]
   	kk5 = B[5]
   	kk6 = B[6]
   	TT1 = B[7]
   	TT2 = B[8]
   	Ians = 0
   	centero = [(0 + TT2[2]) / 2,(kk1[1] + kk1[2]) / 2,(kk2[1] + kk2[2]) / 2,(kk3[1] + kk3[2]) / 2,(kk4[1] + kk4[2]) / 2,(kk5[1] + kk5[2]) / 2,(kk6[1] + kk6[2]) / 2
	        ,(TT2[1] + TT2[2]) / 2]
   	COEFF = TT2[2] / (TT1[2] - TT1[1])
   	t = samplor(0, TT2[2])
   	k1 = samplor(kk1[1], kk1[2])
   	k2 = samplor(kk2[1], kk2[2])
   	k3 = samplor(kk3[1], kk3[2])
   	k4 = samplor(kk4[1], kk4[2])
   	k5 = samplor(kk5[1], kk5[2])
   	k6 = samplor(kk6[1], kk6[2])
   	T2 = samplor(TT2[1], TT2[2])
   	Ians += IntegrandWithoutCOEFF(t, TT1[1], TT1[2], T2, k1, k2, k3, k4, k5, k6)
   	Ians = COEFF * Ians
   	t = 2 * centero[1] - t
   	k1 = 2 * centero[2] - k1
   	k2 = 2 * centero[3] - k2
   	k3 = 2 * centero[4] - k3
   	k4 = 2 * centero[5] - k4
   	k5 = 2 * centero[6] - k5
   	k6 = 2 * centero[7] - k6
   	T2 = 2 * centero[8] - T2
   	delto = COEFF * IntegrandWithoutCOEFF(t, TT1[1], TT1[2], T2, k1, k2, k3, k4, k5, k6) - Ians
   	Ians += delto / 2
   	@inbounds for i = 1:floor(NN / 2)
      		t = samplor(0, TT2[2])
      		k1 = samplor(kk1[1], kk1[2])
      		k2 = samplor(kk2[1], kk2[2])
      		k3 = samplor(kk3[1], kk3[2])
      		k4 = samplor(kk4[1], kk4[2])
      		k5 = samplor(kk5[1], kk5[2])
      		k6 = samplor(kk6[1], kk6[2])
      		T2 = samplor(TT2[1], TT2[2])
      		delto = COEFF * IntegrandWithoutCOEFF(t, TT1[1], TT1[2], T2, k1, k2, k3, k4, k5, k6) - Ians
      		Ians += delto / (2 * i + 1)
      		t = 2 * centero[1] - t
      		k1 = 2 * centero[2] - k1
      		k2 = 2 * centero[3] - k2
      		k3 = 2 * centero[4] - k3
      		k4 = 2 * centero[5] - k4
      		k5 = 2 * centero[6] - k5
      		k6 = 2 * centero[7] - k6
      		T2 = 2 * centero[8] - T2
      		delto = COEFF * IntegrandWithoutCOEFF(t, TT1[1], TT1[2], T2, k1, k2, k3, k4, k5, k6) - Ians
      		Ians += delto / (2 * i + 2)
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

# Now doing a bisect research.
Binput=[[0.0,1.0],[0.0,200.0],[0.0,100.0],[0.0,100.0],[0.0,200.0],[0.0,10.0],[0.0,5.0],[0.0,5.0]]
bisectlimit=30
minimumsizeinput=[0.1,10.0,10.0,10.0,10.0,1.0,0.5,0.5]
print("\nThe bisect search for HK 1v 8p with initial box\n$Binput\nand the constrain on the minimum size of each interval\n$minimumsizeinput\nand a limit on the maximum number of bisect steps $bisectlimit.\n")
X1=BisectSearch(Antithetic,Binput,minimumsizeinput,bisectlimit)
print("\nIt took $(X1[3]) seconds and computed $(X1[2]) integrals to conclude that\n$(X1[1])\n")
if X1[1]=="Found Multistationary subbox"
	X2=sumo_antithetic_with_SE_t(X1[4][2][1],X1[4][2][2],X1[4][2][3],X1[4][2][4],X1[4][2][5],X1[4][2][6],X1[4][2][7],X1[4][2][8],10^7)
	print("The following subbox is almost inside the multistationary region.\n$(X1[4][2])\nWith I-hat equal to $(X1[4][3]). Recomputing I-hat; $(X2[1]) with standard error $(X2[2]).\n")
end

# The output
#=
The bisect search for HK 1v 8p with initial box
[[0.0, 1.0], [0.0, 200.0], [0.0, 100.0], [0.0, 100.0], [0.0, 200.0], [0.0, 10.0], [0.0, 5.0], [0.0, 5.0]]
and the constrain on the minimum size of each interval
[0.1, 10.0, 10.0, 10.0, 10.0, 1.0, 0.5, 0.5]
and a limit on the maximum number of bisect steps 30.

It took 58.2698015 seconds and computed 45 integrals to conclude that
Found Multistationary subbox
The following subbox is almost inside the multistationary region.
[[0.125, 0.25], [125.0, 150.0], [75.0, 87.5], [12.5, 25.0], [175.0, 200.0], [2.5, 3.75], [3.75, 5.0], [3.75, 5.0]]
With I-hat equal to 3.0028351242975724. Recomputing I-hat; 2.9722792561449305 with standard error 0.008981473280893281.
=#

# The output for another try. Remember there is a randomness involved and the accuaracy used here is 2 digits of significant. Hence different boxes might get chosen in different tries.
#=
The bisect search for HK 1v 8p with initial box
[[0.0, 1.0], [0.0, 200.0], [0.0, 100.0], [0.0, 100.0], [0.0, 200.0], [0.0, 10.0], [0.0, 5.0], [0.0, 5.0]]
and the constrain on the minimum size of each interval
[0.1, 10.0, 10.0, 10.0, 10.0, 1.0, 0.5, 0.5]
and a limit on the maximum number of bisect steps 30.

It took 69.786018299 seconds and computed 51 integrals to conclude that
Found Multistationary subbox
The following subbox is almost inside the multistationary region.
[[0.0625, 0.125], [100.0, 125.0], [87.5, 100.0], [25.0, 37.5], [175.0, 200.0], [1.25, 2.5], [3.125, 3.75], [2.5, 3.125]]
With I-hat equal to 2.9754039962531924. Recomputing I-hat; 2.9742987047875595 with standard error 0.013676431991910536.
=#

# Now a bisect description.
@fastmath function BisectDescription_t(fun, B, minimum_size,sh_terminate)
    st1=time_ns()
    d = length(B)
	shlist=[]
	for i=1:d
		push!(shlist,i)
	end
	for i=1:d-1
		push!(shlist,i)
	end
    PartitionList=[]
    UndecidedList=[]
    ans0 = fun(B)
    integral_number=1
    if ans0 < 1.05
        push!(PartitionList, [SameBox(B),ans0])
    elseif ans0 > 2.95
        push!(PartitionList, [SameBox(B),ans0])
    else
        sh1 = 1
		sizeend=1
		for i=sh1:sh1+d-1
			if B[shlist[i]][2]-B[shlist[i]][1] > minimum_size[shlist[i]]
				sizeend=0
				B1 = SameBox(B)
		        B1[shlist[i]][2] = (B[shlist[i]][1] + B[shlist[i]][2]) / 2
		        push!(UndecidedList, [SameBox(B1),shlist[i],1])
		        B2 = SameBox(B)
		        B2[shlist[i]][1] = (B[shlist[i]][1] + B[shlist[i]][2]) / 2
		        push!(UndecidedList, [SameBox(B2),shlist[i],1])
				break
			end
		end
		if sizeend==1
			push!(PartitionList, [SameBox(B),ans0])
		end
	end
    @inbounds while length(UndecidedList) != 0
        Bsho = UndecidedList[1]
        deleteat!(UndecidedList, 1)
        B = Bsho[1]
        sh1 = Bsho[2]
        sh2 = Bsho[3]
        ans0 = fun(B)
        integral_number+=1
        if ans0 < 1.05
            push!(PartitionList, [SameBox(B),ans0])
        elseif ans0 > 2.95
            push!(PartitionList, [SameBox(B),ans0])
        elseif sh2 == sh_terminate
            push!(PartitionList, [SameBox(B),ans0])
        else
            sh2 += 1
            if sh1 == d
                sh1 = 1
            else
                sh1 += 1
            end
			sizeend=1
			for i=sh1:sh1+d-1
				if B[shlist[i]][2]-B[shlist[i]][1] > minimum_size[shlist[i]]
					sizeend=0
					B1 = SameBox(B)
			        B1[shlist[i]][2] = (B[shlist[i]][1] + B[shlist[i]][2]) / 2
			        push!(UndecidedList, [SameBox(B1),shlist[i],sh2])
			        B2 = SameBox(B)
			        B2[shlist[i]][1] = (B[shlist[i]][1] + B[shlist[i]][2]) / 2
			        push!(UndecidedList, [SameBox(B2),shlist[i],sh2])
					break
				end
			end
			if sizeend==1
				push!(PartitionList, [SameBox(B),ans0])
			end
        end
    end
    st2=time_ns()
    return PartitionList,integral_number,(st2-st1)/(10^9),UndecidedList
end

# On the example.
filename="Output_Julia_Subsec_2_1_with_8_parameters_Problem_I"
Binput=[[0.125,0.375],[100.0,125.0],[75.0,100.0],[12.5,37.5],[150.0,200.0],[1.25,3.75],[0.0,5.0],[0.0,5]]
bisectlimit=12
minimumsizeinput=[0.125,12.5,12.5,12.5,25.0,1.25,1.25,1.25]
print("\nThe bisect description for HK 1v 8p with initial box\n$Binput\nand the constrain on the minimum size of each interval\n$minimumsizeinput\nand a limit on the maximum number of bisect steps $bisectlimit.\n")
X1=BisectDescription_t(Antithetic,Binput,minimumsizeinput,bisectlimit)
count1=0
count3=0
countbetween=0
for i=1:length(X1[1])
	if X1[1][i][2] < 1.05
		global count1+=1
	elseif X1[1][i][2] > 1.95
		global count3+=1
	else
		global countbetween+=1
	end
end
print("\nIt took $(X1[3]) seconds and computed $(X1[2]) integrals. Number of sub-boxes in the output is $(length(X1[1])) from which\n$count1 sub-boxes have I-hat less than 1.05,\n$count3 sub-boxes have I-hat greater than 2.95,\n$countbetween sub-boxes have I-hat between 1.05 and 2.95.\n\nAll sub-boxes and their I-hat are saved in the output text file.\n")
File=open(joinpath(pwd(),"$filename.txt"),"w")
write(File,"Output_BisectDescription_HK1v8p_1\n\n\nThe bisect description for HK 1v 8p with initial box\n$Binput\nand the constrain on the minimum size of each interval\n$minimumsizeinput\nand a limit on the maximum number of bisect steps $bisectlimit.\n\nIt took $(X1[3]) seconds and computed $(X1[2]) integrals. Number of sub-boxes in the output is $(length(X1[1])) from which\n$count1 sub-boxes have I-hat less than 1.05,\n$count3 sub-boxes have I-hat greater than 2.95,\n$countbetween sub-boxes have I-hat between 1.05 and 2.95.\n\nAll sub-boxes and their I-hat are listed below.\n\n")
for i=1:length(X1[1])
	write(File,"$(X1[1][i][1])    Ians = $(X1[1][i][2])\n")
end
write(File,"\nThe end of the report")
close(File)

# The output
#=
The bisect description for HK 1v 8p with initial box
[[0.125, 0.375], [100.0, 125.0], [75.0, 100.0], [12.5, 37.5], [150.0, 200.0], [1.25, 3.75], [0.0, 5.0], [0.0, 5.0]]
and the constrain on the minimum size of each interval
[0.125, 12.5, 12.5, 12.5, 25.0, 1.25, 1.25, 1.25]
and a limit on the maximum number of bisect steps 12.

It took 2163.8358728 seconds and computed 1787 integrals. Number of sub-boxes in the output is 894 from which
204 sub-boxes have I-hat less than 1.05,
334 sub-boxes have I-hat greater than 2.95,
356 sub-boxes have I-hat between 1.05 and 2.95.

All sub-boxes and their I-hat are saved in the output text file.
=#

# The end of the file.
