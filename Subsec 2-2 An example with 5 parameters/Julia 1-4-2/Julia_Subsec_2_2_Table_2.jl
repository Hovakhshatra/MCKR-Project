# In this file the Monte Carlo integrations of Table 2 will be done with and without a paralellization using 2 workers.

# Loading the package for parallel computation.
using Distributed
addprocs(2) # I added two extra workers (using two more cpus than the one already used as the main worker of Julia).
@everywhere include(joinpath(pwd(),"Julia_Subsec_2_2_Table_2_module.jl"))
print("\nConsider the Monte Carlo integration for the Network 2v5p example on the box\n[0,100],[0,2],[0,200],[0,100],[0,2]\nwith uniform distributions and N=10^9.\n")
print("\nWithout parallelization.\n")
function test_simple_nonparal(NNinput)
    st=time_ns()
    ans1=sumo_with_S([0,100],[0,2],[0,200],[0,100],[0,2],NNinput)
    standard_error=ans1[2]/NNinput
    standard_error=standard_error/(NNinput-1)
    standard_error=sqrt(standard_error)
    st=time_ns()-st
    return(ans1[1],standard_error,st/10^9)
end
function test_antithetic_nonparal(NNinput)
    st=time_ns()
    ans1=sumo_antithetic_with_S([0,100],[0,2],[0,200],[0,100],[0,2],NNinput)
    standard_error=ans1[2]/NNinput
    standard_error=standard_error/(NNinput-1)
    standard_error=sqrt(standard_error)
    st=time_ns()-st
    return(ans1[1],standard_error,st/10^9)
end
# warm-up round for Julia to digest the function test_t1.
test_simple_nonparal(1000)
test_antithetic_nonparal(1000)
# Now serious round.
out1_1=test_simple_nonparal(10^9)
print("Simple Monte Carlo:\nI-hat = $(out1_1[1]), e-hat = $(out1_1[2]), time = $(out1_1[3])\n")
out1_2=test_antithetic_nonparal(10^9)
print("Antithetic Monte Carlo:\nI-hat = $(out1_2[1]), e-hat = $(out1_2[2]), time = $(out1_2[3])\n")
print("\nWith parallelization and using only 2 workers (workers number 2 and 3).\n")
function test_simple_paral(NNinput)
    st=time_ns()
    w1ans1=remotecall(sumo_with_S,2,[0,100],[0,2],[0,200],[0,100],[0,2],Int(floor(NNinput/2)))
    w2ans1=remotecall(sumo_with_S,3,[0,100],[0,2],[0,200],[0,100],[0,2],Int(floor(NNinput/2)))
    ans2w1=fetch(w1ans1)
    ans2w2=fetch(w2ans1)
    ans2_1=(ans2w1[1]+ans2w2[1])/2
    ans2_2=ans2w1[2]+ans2w2[2]
    standard_error=ans2_2/NNinput
    standard_error=standard_error/(NNinput-1)
    standard_error=sqrt(standard_error)
    st=time_ns()-st
    return(ans2_1,standard_error,st/10^9)
end
function test_antithetic_paral(NNinput)
    st=time_ns()
    w1ans1=remotecall(sumo_antithetic_with_S,2,[0,100],[0,2],[0,200],[0,100],[0,2],Int(floor(NNinput/2)))
    w2ans1=remotecall(sumo_antithetic_with_S,3,[0,100],[0,2],[0,200],[0,100],[0,2],Int(floor(NNinput/2)))
    ans2w1=fetch(w1ans1)
    ans2w2=fetch(w2ans1)
    ans2_1=(ans2w1[1]+ans2w2[1])/2
    ans2_2=ans2w1[2]+ans2w2[2]
    standard_error=ans2_2/NNinput
    standard_error=standard_error/(NNinput-1)
    standard_error=sqrt(standard_error)
    st=time_ns()-st
    return(ans2_1,standard_error,st/10^9)
end
# warm-up round for Julia to digest the function test_t2.
test_simple_paral(1000)
test_antithetic_paral(1000)
# Now serious round.
out2_1=test_simple_paral(10^9)
print("Simple Monte Carlo:\nI-hat = $(out2_1[1]), e-hat = $(out2_1[2]), time = $(out2_1[3])\n")
out2_2=test_antithetic_paral(10^9)
print("Antithetic Monte Carlo:\nI-hat = $(out2_2[1]), e-hat = $(out2_2[2]), time = $(out2_2[3])\n")

#= Output
Consider the Monte Carlo integration for the Network 2v5p example on the box
[0,100],[0,2],[0,200],[0,100],[0,2]
with uniform distributions and N=10^9.

Without parallelization.
Simple Monte Carlo:
I-hat = 1.4221338619269406, e-hat = 0.0028305697695442784, time = 156.2836928
Antithetic Monte Carlo:
I-hat = 1.4190118168649908, e-hat = 0.0012688965814635043, time = 88.815760299

With parallelization and using only 2 workers (workers number 2 and 3).
Simple Monte Carlo:
I-hat = 1.4187853771514467, e-hat = 0.0015029248979047275, time = 76.165509701
Antithetic Monte Carlo:
I-hat = 1.4212984605157672, e-hat = 0.001935801329076283, time = 47.355208499
=#
rmprocs(2,waitfor=0)
rmprocs(3,waitfor=0)

# The end of the file.
