# The Monte Carlo integrations used in Figure 1 of the paper.

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

# First we check which of Simple or Antithetic Monte Carlo is better for this example and how big the sample size is suitable.
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

@fastmath function Simple_SE_t(B,NN)
    # Simple Monte Carlo integration for the HK 1v 2p network with B as the parameter region and NN as the sample size.
    # B is a list of 2 lists. Each sublist consists of two numbers.
    # NN is an integer.
    # The output is a list of three numbers; I-hat, e-hat and the computation time.
    TT1 = B[1]
    TT2 = B[2]
    Ians = 0
    S = 0
    st1 = time_ns()
    COEFF = TT2[2] / (TT1[2] - TT1[1])
    t = samplor(0, TT2[2])
    T2 = samplor(TT2[1], TT2[2])
    Ians += IntegrandWithoutCOEFF(t, TT1[1], TT1[2], T2)
    Ians = COEFF * Ians
    @inbounds for i = 2:NN
        t = samplor(0, TT2[2])
   		T2 = samplor(TT2[1], TT2[2])
   		delto = COEFF * IntegrandWithoutCOEFF(t, TT1[1], TT1[2], T2) - Ians
        Ians += delto / i
      	S += delto^2 * (i - 1) / i
    end
    Estandardo = S / NN
   	Estandardo = Estandardo / (NN - 1)
   	Estandardo = sqrt(Estandardo)
   	st2 = time_ns()
   	return Ians, Estandardo, (st2 - st1) / 10^9
end

@fastmath function Antithetic_SE_t(B,NN)
    # Antithetic Monte Carlo integration for the HK 1v 2p network with B as the parameter region and NN as the sample size.
    # B is a list of 2 lists. Each sublist consists of two numbers.
    # NN is an integer.
    # The output is a list of three numbers; I-hat, e-hat and the computation time.
    TT1 = B[1]
    TT2 = B[2]
    Ians = 0
    S = 0
    st1 = time_ns()
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
   	S += delto^2 / 2
    @inbounds for i = 1:floor(NN / 2)
        t = samplor(0, TT2[2])
   		T2 = samplor(TT2[1], TT2[2])
   		delto = COEFF * IntegrandWithoutCOEFF(t, TT1[1], TT1[2], T2) - Ians
        Ians += delto / (2 * i + 1)
      	S += (delto^2) * 2 * i / (2 * i + 1)
        t = 2 * centero[1] - t
       	T2 = 2 * centero[2] - T2
        delto = COEFF * IntegrandWithoutCOEFF(t, TT1[1], TT1[2], T2) - Ians
        Ians += delto / (2 * i + 2)
      	S += (delto^2) * (2 * i + 1) / (2 * i + 2)
    end
    Estandardo = S / NN
   	Estandardo = Estandardo / (NN - 1)
   	Estandardo = sqrt(Estandardo)
   	st2 = time_ns()
   	return Ians, Estandardo, (st2 - st1) / 10^9
end

# Trying the two Monte Carlo methods with NN=10^n for n=1,...,8 on the parameter box B=[0,5]^2.
print("\nB=[0,5]^2.\n")
print("\nWith Simple Monte Carlo:\n")
for n=1:8
    print("N=10^$n\n")
    X1=Simple_SE_t([[0.0,5.0],[0.0,5.0]],10^n)
    print("Ians = $(X1[1]), Standard error = $(X1[2]), Computation time = $(X1[3]).\n")
end
print("\nWith Antithetic Monte Carlo:\n")
for n=1:8
    print("N=10^$n\n")
    X1=Antithetic_SE_t([[0.0,5.0],[0.0,5.0]],10^n)
    print("Ians = $(X1[1]), Standard error = $(X1[2]), Computation time = $(X1[3]).\n")
end

# The output.
#=
B=[0,5]^2.

With Simple Monte Carlo:
N=10^1
N=10^1
N=10^1
Ians = 0.7446344938603167, Standard error = 0.49903360357101906, Computation time = 4.9e-6.
N=10^2
Ians = 2.047503493580139, Standard error = 0.6780052573359778, Computation time = 9.1e-6.
N=10^3
Ians = 1.1154310780138168, Standard error = 0.09269148888609961, Computation time = 7.2801e-5.
N=10^4
Ians = 1.2979001112153017, Standard error = 0.06546965873186893, Computation time = 0.0006951.
N=10^5
Ians = 1.37435399597473, Standard error = 0.04548089909669083, Computation time = 0.0069045.
N=10^6
Ians = 1.3557356166107013, Standard error = 0.017076681973511542, Computation time = 0.0687607.
N=10^7
Ians = 1.3645991532555035, Standard error = 0.003689604999780314, Computation time = 0.6719161.
N=10^8
Ians = 1.3647211786802147, Standard error = 0.0014695839509436884, Computation time = 6.724865001.

With Antithetic Monte Carlo:
N=10^1
Ians = 0.666902196109972, Standard error = 0.5467176767261079, Computation time = 0.0106682.
N=10^2
Ians = 1.8023098856550501, Standard error = 0.702157982677682, Computation time = 1.18e-5.
N=10^3
Ians = 1.262883091683313, Standard error = 0.14204791712487155, Computation time = 5.37e-5.
N=10^4
Ians = 1.2449037010837478, Standard error = 0.05901648335374041, Computation time = 0.000528101.
N=10^5
Ians = 1.4223770942806266, Standard error = 0.08176279034781933, Computation time = 0.004989401.
N=10^6
Ians = 1.3770340942017982, Standard error = 0.014566186652320299, Computation time = 0.0498795.
N=10^7
Ians = 1.3611053086169516, Standard error = 0.004285848275309889, Computation time = 0.484691399.
N=10^8
Ians = 1.3686797629788985, Standard error = 0.0017684942861203343, Computation time = 5.1144983.

=#

# The conclusion is that the suitables sample size for 2 digits of significant is 10^7 in both methods.
# The Antithetic method takes 2/3 time of the Simple method.
# So we choose the Antithetic Monte Carlo with sample size 10^7 in all of the following integrations.

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

# The 100 integrations needed for Figure 1 (b).
print("\nThe average number of solutions for the subboxes [(i-1)*h,i*h]x[(j-1)*h,j*h] where h=(5-0)/10 and i,j=1,...,10.\n")
GridDescription=zeros(Float64,10,10)
hh=(5-0)/10 # h is a predefined constant in Julia.
for i1=1:10
    m1=0+(i1-1)*hh
    for i2=1:10
        m2=0+(i2-1)*hh
        global GridDescription[i1,i2]=Antithetic([[m1,m1+hh],[m2,m2+hh]])
    end
end
print("\n",GridDescription,"\n")

# The output
#=
The average number of solutions for the subboxes [(i-1)*h,i*h]x[(j-1)*h,j*h] where h=(5-0)/10 and i,j=1,...,10.

[1.002607779361504 1.000207477696043 0.9979342685358096 1.0007488614204831 0.9996300303692774 1.0007933479067517 1.0004584573991566 1.0032150144667193 1.000274884351453 0.9970005105522692; 1.003838450251899 1.1078909702070254 1.0484848160399125 1.000585469292574 1.0013261150439494 0.9996303614383619 1.0027069687529782 0.9988411939833595 1.0006953315930314 0.9971718189789102; 0.9975455867359654 1.003835313930864 1.7896085779582291 1.0402070734269182 1.000072266940248 1.002323007095388 1.001807448068555 1.003472206995995 1.0013747948112988 0.9995087763366592; 1.004540225101247 0.9996982072699336 1.3614509725038022 2.33927949198663 1.0337987169653202 1.0017511627823696 0.9988507606145585 0.9989415757666538 1.0048558037487059 0.9987388814545003; 0.9955898287525814 0.9933639312722956 0.9968909964799657 2.396863801491822 2.3203165175750735 1.0272689415421605 0.9989741469658004 0.9977291506042594 0.9989627260146906 1.001087516157902; 1.0028691410448238 1.0093767412002561 1.0006158367813047 1.518945060363419 2.9994369368019966 2.2944950925445893 1.0223067855429029 0.9982475313774372 1.0004876621315402 0.9984519216580742; 0.9871699817775198 1.0154104234848025 1.0089629049302764 1.0013320806796882 2.6968227349079315 2.9940203702043067 2.2833407985220564 1.021101338595854 0.9996647897265939 0.9973678174112808; 1.0048332146030425 0.9876389185915351 1.008339391158617 0.996523984176243 1.9069905213175813 3.0065410068199827 2.999904334503541 2.254161089755452 1.018471221436173 0.9961324068193056; 0.9819871814606389 0.9941744106628547 1.0032187966375896 1.0080859486417777 1.194665357806602 2.9618541376005516 3.007977456197367 2.9997510876010414 2.229277491537151 1.0136516979706223; 0.9705402535935708 0.9972597879861738 0.9975501545963363 1.0012869153443105 1.0082579399027594 2.44419904817351 3.0048594336797962 3.0019792713421376 3.0003309352175136 2.201463087510853]
=#

# The above result is used to plot Figure 1 (b) by Maple in the file "Maple_4_Example_HK_1v_2p_Figure_1.mw".

function SameBox(B) # To copy a box without dependency between the copy and the original.
    return(deepcopy(B))
end

@fastmath function BisectDescription_t(fun, B, minimum_size,sh_terminate)
    # Using the bisect strategy to make an approximation of the parameter region similar to the grid description.
    # This function should be used when M_min=1 and M_max=3.
    # fun is the function which associate a number to each subbox, B is the initial box, minimum_size is the minimum size condition for edges of the subboxes, sh_terminate is the upper bound for number of bisect steps.
    # fun is a function, B is a list of list of numbers, minimum_size is a list of numbers, sh_terminater is an integer.
    # The output is a 4-tuple.
    # A list of subboxes that are not going to be bisected anymore together with their value of "fun".
    # Number of times function "fun" is evaluated at a subbox.
    # The amount of time the computation took.
    # A list of subboxes that might be unchecked by any reason.
    # The termination condition is determined by number of bisect steps, all edges being smaller than minimum sizes, or the subbozes having value of fun greater than 2.95 or less than 1.05.
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

# Here we ask Julia to show only the summary of the report on terminal and list the details of the subboxes in a txt file.
filename="Output_4_Example_HK_1v_2p_Figure_1_c"
Binput=[[0.0,5.0],[0.0,5.0]]
bisectlimit=10
minimumsizeinput=[0.5,0.5]
print("\nThe bisect description for HK 1v 2p with initial box\n$Binput\nand the constrain on the minimum size of each interval\n$minimumsizeinput\nand a limit on the maximum number of bisect steps $bisectlimit.\n")
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
write(File,"Output_BisectDescription_HK1v2p\n\n\nThe bisect description for HK 1v 2p with initial box\n$Binput\nand the constrain on the minimum size of each interval\n$minimumsizeinput\nand a limit on the maximum number of bisect steps $bisectlimit.\n\nIt took $(X1[3]) seconds and computed $(X1[2]) integrals. Number of sub-boxes in the output is $(length(X1[1])) from which\n$count1 sub-boxes have I-hat less than 1.05,\n$count3 sub-boxes have I-hat greater than 2.95,\n$countbetween sub-boxes have I-hat between 1.05 and 2.95.\n\nAll sub-boxes and their I-hat are listed below.\n\n")
for i=1:length(X1[1])
	write(File,"$(X1[1][i][1])    Ians = $(X1[1][i][2])\n")
end
write(File,"\nThe end of the report")
close(File)

# The output
#=
The bisect description for HK 1v 2p with initial box
[[0.0, 5.0], [0.0, 5.0]]
and the constrain on the minimum size of each interval
[0.5, 0.5]
and a limit on the maximum number of bisect steps 10.

It took 58.554622099 seconds and computed 123 integrals. Number of sub-boxes in the output is 62 from which
20 sub-boxes have I-hat less than 1.05,
33 sub-boxes have I-hat greater than 2.95,
9 sub-boxes have I-hat between 1.05 and 2.95.

All sub-boxes and their I-hat are saved in the output text file.
=#

# The result is used to plot Figure 1 (c) by Maple in the file "Maple_4_Example_HK_1v_2p_Figure_1.mw".
# BisectDescriptionList in the Maple file is the following.
print(X1[1])
# The output
#=
[[[[0.0, 2.5], [2.5, 5.0]], 0.9992956190455282], Any[[[0.0, 1.25], [0.0, 2.5]], 1.04547777365172], Any[[[1.25, 2.5], [0.0, 1.25]], 1.0125608727776307], Any[[[2.5, 3.75], [0.0, 1.25]], 1.0032982559521606], Any[[[3.75, 5.0], [0.0, 1.25]], 0.9922549022446409], Any[[[2.5, 3.75], [3.75, 5.0]], 1.0013539707927226], Any[[[4.375, 5.0], [1.25, 2.5]], 0.9981669991569515], Any[[[3.75, 4.375], [2.5, 3.75]], 2.9998536448640816], Any[[[1.25, 1.875], [1.875, 2.5]], 1.021291127357108], Any[[[3.125, 3.75], [1.25, 1.875]], 0.9958082468378373], Any[[[3.75, 4.375], [1.25, 1.875]], 1.0080786635174117], Any[[[2.5, 3.125], [3.125, 3.75]], 1.0138581488337346], Any[[[3.125, 3.75], [2.5, 3.125]], 2.997052600896121], Any[[[4.375, 5.0], [3.125, 3.75]], 2.992705169150627], Any[[[3.75, 4.375], [4.375, 5.0]], 1.011193672628307], Any[[[4.375, 5.0], [3.75, 4.375]], 2.9940643385664827], Any[[[2.8125, 3.125], [1.25, 1.875]], 0.9982874539940778], Any[[[2.5, 2.8125], [1.875, 2.5]], 3.004570709397327], Any[[[1.25, 1.5625], [1.25, 1.5625]], 2.4564676995861436], Any[[[1.25, 1.5625], [1.5625, 1.875]], 1.1001353068532693], Any[[[1.5625, 1.875], [1.25, 1.5625]], 2.0570180653101304], Any[[[1.5625, 1.875], [1.5625, 1.875]], 2.5226268788708675], Any[[[1.875, 2.1875], [1.25, 1.5625]], 1.1682269225606967], Any[[[1.875, 2.1875], [1.5625, 1.875]], 2.9229335543257378], Any[[[2.1875, 2.5], [1.25, 1.5625]], 1.001929341189325], Any[[[2.1875, 2.5], [1.5625, 1.875]], 2.160914005756467], Any[[[1.875, 2.1875], [1.875, 2.1875]], 2.5007586651643345], Any[[[1.875, 2.1875], [2.1875, 2.5]], 1.0812067806833159], Any[[[2.1875, 2.5], [1.875, 2.1875]], 3.004475251191103], Any[[[2.1875, 2.5], [2.1875, 2.5]], 2.478591416985572], Any[[[2.5, 2.8125], [1.25, 1.5625]], 0.9964943775822951], Any[[[2.5, 2.8125], [1.5625, 1.875]], 1.289620445842181], Any[[[2.8125, 3.125], [1.875, 2.1875]], 2.4382252847141297], Any[[[2.8125, 3.125], [2.1875, 2.5]], 3.002260631057597], Any[[[3.125, 3.4375], [1.875, 2.1875]], 1.633474409740284], Any[[[3.125, 3.4375], [2.1875, 2.5]], 2.9991306478720925], Any[[[3.4375, 3.75], [1.875, 2.1875]], 1.02655351719798], Any[[[3.4375, 3.75], [2.1875, 2.5]], 2.80134235403927], Any[[[3.75, 4.0625], [1.875, 2.1875]], 0.9904359872325619], Any[[[3.75, 4.0625], [2.1875, 2.5]], 2.057075192189435], Any[[[4.0625, 4.375], [1.875, 2.1875]], 1.0063039421350919], Any[[[4.0625, 4.375], [2.1875, 2.5]], 1.3220238388036287], Any[[[2.5, 2.8125], [2.5, 2.8125]], 2.4584154521450765], Any[[[2.5, 2.8125], [2.8125, 3.125]], 1.065783482121664], Any[[[2.8125, 3.125], [2.5, 2.8125]], 2.9954992737028516], Any[[[2.8125, 3.125], [2.8125, 3.125]], 2.440568040496486], Any[[[3.125, 3.4375], [3.125, 3.4375]], 2.4173825340123423], Any[[[3.125, 3.4375], [3.4375, 3.75]], 1.0584075129774577], Any[[[3.4375, 3.75], [3.125, 3.4375]], 2.9976751718522423], Any[[[3.4375, 3.75], [3.4375, 3.75]], 2.3995027141595764], Any[[[4.375, 4.6875], [2.5, 2.8125]], 2.58435679630783], Any[[[4.375, 4.6875], [2.8125, 3.125]], 2.9923128863938064], Any[[[4.6875, 5.0], [2.5, 2.8125]], 1.8825341573572987], Any[[[4.6875, 5.0], [2.8125, 3.125]], 3.0083527486246173], Any[[[3.75, 4.0625], [3.75, 4.0625]], 2.3799956888430946], Any[[[3.75, 4.0625], [4.0625, 4.375]], 1.0420884084401296], Any[[[4.0625, 4.375], [3.75, 4.0625]], 2.995225659098818], Any[[[4.0625, 4.375], [4.0625, 4.375]], 2.365443116170008], Any[[[4.375, 4.6875], [4.375, 4.6875]], 2.338353729405777], Any[[[4.375, 4.6875], [4.6875, 5.0]], 1.0338710227891945], Any[[[4.6875, 5.0], [4.375, 4.6875]], 3.0095272242501125], Any[[[4.6875, 5.0], [4.6875, 5.0]], 2.322891445029989]]
=#

# The end of the file.
