# The Monte Carlo integrations used in Figure 2 of the paper.

# This file can be understood better if read after the file "Julia_4_Example_HK_1v_2p_Figure_1.jl".

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

@fastmath function h(t)
   	return 366450 * t^2 + 537142.41 * t
end
@fastmath function beta(t)
   	return 2518322.5 * t^2 + 63502.1205 * t + 26857.1205
end
@fastmath function g1(t, T2) # g_{T_2}(t)
   	return beta(t) * (T2 - t) / h(t)
end
@fastmath function detJf(t, T1, T2)
   	return 537142.41 * T1 - 63502.1205 * T2 + 7554967.5 * t^2 + 2 * t * (366450 * T1 - 2518322.5 * T2 + 63502.1205) + 26857.1205
end
@fastmath function J(t, T2)
   	return abs(detJf(t, g1(t, T2), T2))
end
@fastmath function IntegrandWithoutCOEFF(t, aT1, bT1, T2)
   	return J(t, T2) * chi(aT1, bT1, g1(t, T2)) / h(t)
end

@fastmath function Antithetic(B)
    # Antithetic Monte Carlo integration for the HK 1v 2p network with B as the parameter region and sample size equal to 10^7.
    # B is a list of 2 lists. Each sublist consists of two numbers.
    # The output is a numbers; I-hat.
    NN = 10^7
    TT1 = B[1]
    TT2 = B[2]
    Ians = 0
    centero = [(0 + TT2[2]) / 2,(TT2[1] + TT2[2]) / 2]
    COEFF = TT2[2] / (TT1[2] - TT1[1])
    t = samplor(0, TT2[2])
    T2 = samplor(TT2[1], TT2[2])
    Ians += IntegrandWithoutCOEFF(t, TT1[1], TT1[2], T2)
    Ians = COEFF * Ians
    t = 2 * centero[1] - t
   	T2 = 2 * centero[2] - T2
   	delto = COEFF * IntegrandWithoutCOEFF(t, TT1[1], TT1[2], T2) - Ians
    Ians += delto / 2
    @inbounds for i = 1:floor(NN / 2)
        t = samplor(0, TT2[2])
   		T2 = samplor(TT2[1], TT2[2])
   		delto = COEFF * IntegrandWithoutCOEFF(t, TT1[1], TT1[2], T2) - Ians
        Ians += delto / (2 * i + 1)
        t = 2 * centero[1] - t
       	T2 = 2 * centero[2] - T2
        delto = COEFF * IntegrandWithoutCOEFF(t, TT1[1], TT1[2], T2) - Ians
        Ians += delto / (2 * i + 2)
    end
   	return Ians
end

function SameBox(B) # To copy a box without dependency between the copy and the original.
    return(deepcopy(B))
end

@fastmath function BisectSearch_SeeEverySteps(fun, B, minimum_size,sh_terminate)
    # The bisect search for finding a subbox with value of "fun" equal to 3.
    # This function should be used when M_min=1 and M_max=3.
    # fun is the function which associate a number to each subbox, B is the initial box, minimum_size is the minimum size condition for edges of the subboxes, sh_terminate is the upper bound for number of bisect steps.
    # fun is a function, B is a list of list of numbers, minimum_size is a list of numbers, sh_terminater is an integer.
    # The output is a 4-tuple or 5-tuple depending on finding a subbox inside the multistationary region.
    # A string message saying it found multistationary subbox, found no multistationary subbox, needs more bisect steps.
	# A list of lists. Each list contains the number of the bisect step and the two/one subbox in this step and the value of "fun" on them.
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
	SearchList=[[sh2,SameBox(B),ans0]]
    integral_number=1
    if ans0 < 1.05
        return "No multistationary subbox",SearchList,integral_number,(time_ns()-st1)/(10^9)
    elseif ans0 > 2.95
        return "Found Multistationary subbox",SearchList,integral_number,(time_ns()-st1)/(10^9),[sh2,B,ans0]
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
				push!(SearchList,[sh2,SameBox(B1),ans1])
			    B2 = SameBox(B)
		        B2[shlist[i]][1] = (B[shlist[i]][1] + B[shlist[i]][2]) / 2
			    ans2=fun(B2)
				push!(SearchList,[sh2,SameBox(B2),ans2])
				integral_number+=2
				if ans1 > ans2
					B=B1
					ans0=ans1
				else
					B=B2
					ans0=ans2
				end
				if ans0 < 1.05
			        return "No multistationary subbox",SearchList,integral_number,(time_ns()-st1)/(10^9)
			    elseif ans0 > 2.95
			        return "Found Multistationary subbox",SearchList,integral_number,(time_ns()-st1)/(10^9),[sh2,B,ans0]
			    end
				break
			end
		end
		if sizeend==1
			return "Termination condition on size of box occured",SearchList,integral_number,(time_ns()-st1)/(10^9)
		end
    end
    st2=time_ns()
    return "Termination condition on the number of bisect steps occured",SearchList,integral_number,(time_ns()-st1)/(10^9)
end

# The information of table in Figure 2 (a) which is also the data needed to plot Figure 2 (b) by Maple in the file "Maple_Subsec_2_1_Figure_2.mw" is from the following computation.

# We ask Julia to write the report on a txt file as well.
filename="Output_Julia_Subsec_2_1_Figure_2"
Binput=[[1.0,3.0],[2.0,4.0]]
bisectlimit=10
minimumsizeinput=[0.5,0.5]
print("\nThe bisect search for HK 1v 2p with initial box\n$Binput\nand the constrain on the minimum size of each interval\n$minimumsizeinput\nand a limit on the maximum number of bisect steps $bisectlimit.\n")
X1=BisectSearch_SeeEverySteps(Antithetic,Binput,minimumsizeinput,bisectlimit)
print("\nIt took $(X1[4]) seconds and computed $(X1[3]) integrals to conclude that\n$(X1[1])\nAll sub-boxes in the search steps and their I-hat are saved in the output text file.\n")
if X1[1]=="Found Multistationary subbox"
	print("The following subbox is almost inside the multistationary region.\n$(X1[5][2])\nWith I-hat equal to $(X1[5][3])\n")
end
File=open(joinpath(pwd(),"$filename.txt"),"w")
write(File,"\nThe bisect search for HK 1v 2p with initial box\n$Binput\nand the constrain on the minimum size of each interval\n$minimumsizeinput\nand a limit on the maximum number of bisect steps $bisectlimit.\n\nIt took $(X1[4]) seconds and computed $(X1[3]) integrals to conclude that\n$(X1[1])\nAll sub-boxes in the search steps and their I-hat are saved in the output text file.\n")
if X1[1]=="Found Multistationary subbox"
	write(File,"The following subbox is almost inside the multistationary region.\n$(X1[5][2])\nWith I-hat equal to $(X1[5][3])\n\nNow the whole detailes of the search.\n")
end
for i=1:length(X1[2])
	write(File,"Bisect step no. $(X1[2][i][1]), B = $(X1[2][i][2]), Ians = $(X1[2][i][3])\n")
end
write(File,"\nThe end of the report")
close(File)

# The output
#=
The bisect search for HK 1v 2p with initial box
[[1.0, 3.0], [2.0, 4.0]]
and the constrain on the minimum size of each interval
[0.5, 0.5]
and a limit on the maximum number of bisect steps 10.

It took 4.0847225 seconds and computed 9 integrals to conclude that
Found Multistationary subbox
All sub-boxes in the search steps and their I-hat are saved in the output text file.
The following subbox is almost inside the multistationary region.
[[2.5, 3.0], [2.0, 2.5]]
With I-hat equal to 2.992759053764498
=#

# The end of the file.
