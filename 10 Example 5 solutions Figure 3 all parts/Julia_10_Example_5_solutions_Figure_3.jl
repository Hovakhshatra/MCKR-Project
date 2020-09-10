# The Monte Carlo integrations used in Figure 3 (c) of the paper.

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

@fastmath function g2(t, k1)
    return 25 / 4 - 100 * t^5 + 100 * t^4 * k1 + 450 * t^4 - 450 * t^3 * k1 - 525 * t^3 + 575 * t^2 * k1 - 75 / 2 * t^2 - 375 / 2 * k1 * t + 575 / 2 * t;
end

@fastmath function J(t, k1)
    return abs(5 * t^4 + 4 * (-k1 - 9 / 2) * t^3 + 3 * ((9 * k1) / 2 + 21 / 4) * t^2 + 2 * (-(23 * k1) / 4 + 3 / 8) * t + (15 * k1) / 8 - 23 / 8)
end

@fastmath function SummandWithoutCOEFF(t, k1, ak2, bk2)
    return 100 * J(t, k1) * chi(ak2, bk2, g2(t, k1))
end

# The simple Monte-Carlo summation with standard error and the computation time.
@fastmath function sumo_with_SE_t(kk1, kk2, NN)
    S = 0
    Ians = 0
    st1 = time_ns()
    COEFF = 1 / (kk2[2] - kk2[1])
    t = samplor(0, 1)
    k1 = samplor(kk1[1], kk1[2])
    Ians += SummandWithoutCOEFF(t, k1, kk2[1], kk2[2]) + SummandWithoutCOEFF(1 / t, k1, kk2[1], kk2[2]) / t^2
    Ians = COEFF * Ians
    @inbounds for i = 2:NN
        t = samplor(0, 1)
        k1 = samplor(kk1[1], kk1[2])
        delto = COEFF * ( SummandWithoutCOEFF(t, k1, kk2[1], kk2[2]) + SummandWithoutCOEFF(1 / t, k1, kk2[1], kk2[2]) ) / t^2 - Ians
        Ians += delto / i
        S += delto^2 * (i - 1) / i
    end
    Estandardo = S / NN
    Estandardo = Estandardo / (NN - 1)
    Estandardo = sqrt(Estandardo)
    st2 = time_ns()
    return Ians, Estandardo, (st2 - st1) / 10^9
end

# The Antithetic Monte-Carlo summation with standard error and the computation time.
@fastmath function sumo_antithetic_with_SE_t(kk1, kk2, NN)
    S = 0
    Ians = 0
    st1 = time_ns()
    centero = [1 / 2,(kk1[1] + kk1[2]) / 2]
    COEFF = 1 / (kk2[2] - kk2[1])
    t = samplor(0, 1)
    k1 = samplor(kk1[1], kk1[2])
    Ians += SummandWithoutCOEFF(t, k1, kk2[1], kk2[2]) + SummandWithoutCOEFF(1 / t, k1, kk2[1], kk2[2]) / t^2
    Ians = COEFF * Ians
    t = 2 * centero[1] - t
    k1 = 2 * centero[2] - k1
    delto = SummandWithoutCOEFF(t, k1, kk2[1], kk2[2]) + SummandWithoutCOEFF(1 / t, k1, kk2[1], kk2[2]) / t^2 - Ians
    Ians += delto / 2
    S += delto^2 / 2
    @inbounds for i = 1:floor(NN / 2)
        t = samplor(0, 1)
        k1 = samplor(kk1[1], kk1[2])
        delto = SummandWithoutCOEFF(t, k1, kk2[1], kk2[2]) + SummandWithoutCOEFF(1 / t, k1, kk2[1], kk2[2]) / t^2 - Ians
        Ians += delto / (2 * i + 1)
        S += 2 * i * delto^2 / (2 * i + 1)
        t = 2 * centero[1] - t
        k1 = 2 * centero[2] - k1
        delto = SummandWithoutCOEFF(t, k1, kk2[1], kk2[2]) + SummandWithoutCOEFF(1 / t, k1, kk2[1], kk2[2]) / t^2 - Ians
        Ians += delto / (2 * i + 2)
        S += (2 * i + 1) * delto^2 / (2 * i + 2)
    end
    Estandardo = S / NN
    Estandardo = Estandardo / (NN - 1)
    Estandardo = sqrt(Estandardo)
    st2 = time_ns()
    return Ians, Estandardo, (st2 - st1) / 10^9
end

# We use the antithetic Monte Carlo integration with sample size 10^9 for the Kac-Rice integrals of the grid description in Figure 3 (a).
# We also check the standard error of each.
sumo_with_SE_t([4,5], [1,2], 100)
stt = time_ns()
for i = 1:10
    a1 = (i - 1) / 2
    b1 = i / 2
    for j = 1:10
        a2 = j - 1
        b2 = j
        print("Subbox = [", a1, ",", b1, "]x[", a2, ",", b2, "]\n")
        X = sumo_with_SE_t([a1,b1], [a2,b2], 10^9)
        print("With the simple Monte-Carlo: \n", X[1], " = Ians\n", X[2], " = Standard Error\n", X[3], " = Time of computation\n\n")
    end
end
sttlength = time_ns() - stt
print("Total time = ", sttlength / 10^9)

# The output
#=
Subbox = [0.0,0.5]x[0,1]
With the simple Monte-Carlo:
0.9961026912300243 = Ians
0.003120025943018131 = Standard Error
65.087235701 = Time of computation

Subbox = [0.0,0.5]x[1,2]
With the simple Monte-Carlo:
1.0031035112802575 = Ians
0.003125205273631686 = Standard Error
65.420406 = Time of computation

Subbox = [0.0,0.5]x[2,3]
With the simple Monte-Carlo:
1.0015855210514977 = Ians
0.00311784707225434 = Standard Error
65.997850399 = Time of computation

Subbox = [0.0,0.5]x[3,4]
With the simple Monte-Carlo:
0.9976586935139492 = Ians
0.0031062780881145648 = Standard Error
66.0568566 = Time of computation

Subbox = [0.0,0.5]x[4,5]
With the simple Monte-Carlo:
0.9937969352320325 = Ians
0.0030957176491064227 = Standard Error
66.5295158 = Time of computation

Subbox = [0.0,0.5]x[5,6]
With the simple Monte-Carlo:
0.9999805068475218 = Ians
0.0031001010270952276 = Standard Error
65.261243 = Time of computation

Subbox = [0.0,0.5]x[6,7]
With the simple Monte-Carlo:
1.7547124930926001 = Ians
0.003129406550448934 = Standard Error
64.815809 = Time of computation

Subbox = [0.0,0.5]x[7,8]
With the simple Monte-Carlo:
1.997869938139312 = Ians
0.003124686443171286 = Standard Error
65.5823033 = Time of computation

Subbox = [0.0,0.5]x[8,9]
With the simple Monte-Carlo:
1.9987413932953089 = Ians
0.0031216801964509874 = Standard Error
65.755461699 = Time of computation

Subbox = [0.0,0.5]x[9,10]
With the simple Monte-Carlo:
1.9941082157052217 = Ians
0.003108788740043067 = Standard Error
65.530633299 = Time of computation

Subbox = [0.5,1.0]x[0,1]
With the simple Monte-Carlo:
0.9974359872478978 = Ians
0.0028876034725404853 = Standard Error
65.3507844 = Time of computation

Subbox = [0.5,1.0]x[1,2]
With the simple Monte-Carlo:
0.9945789417177416 = Ians
0.002877242636992255 = Standard Error
65.4030904 = Time of computation

Subbox = [0.5,1.0]x[2,3]
With the simple Monte-Carlo:
0.995960597741682 = Ians
0.002873692394426725 = Standard Error
65.2974057 = Time of computation

Subbox = [0.5,1.0]x[3,4]
With the simple Monte-Carlo:
0.9994364609116491 = Ians
0.0028727709933251554 = Standard Error
64.7354112 = Time of computation

Subbox = [0.5,1.0]x[4,5]
With the simple Monte-Carlo:
1.0017903808330813 = Ians
0.002870067470392339 = Standard Error
64.4709642 = Time of computation

Subbox = [0.5,1.0]x[5,6]
With the simple Monte-Carlo:
0.9966809477211206 = Ians
0.002856862144457869 = Standard Error
65.4586645 = Time of computation

Subbox = [0.5,1.0]x[6,7]
With the simple Monte-Carlo:
1.7495628411193607 = Ians
0.0028739355992350897 = Standard Error
65.3429244 = Time of computation

Subbox = [0.5,1.0]x[7,8]
With the simple Monte-Carlo:
2.000908546486986 = Ians
0.002878381597029554 = Standard Error
65.4902307 = Time of computation

Subbox = [0.5,1.0]x[8,9]
With the simple Monte-Carlo:
2.004424447170691 = Ians
0.002877300220148343 = Standard Error
65.0825319 = Time of computation

Subbox = [0.5,1.0]x[9,10]
With the simple Monte-Carlo:
1.998955104558242 = Ians
0.0028638554660063323 = Standard Error
65.0895606 = Time of computation

Subbox = [1.0,1.5]x[0,1]
With the simple Monte-Carlo:
0.9987536538019167 = Ians
0.002655682781884911 = Standard Error
74.6859563 = Time of computation

Subbox = [1.0,1.5]x[1,2]
With the simple Monte-Carlo:
0.9977266252382966 = Ians
0.0026476256173974763 = Standard Error
66.241061899 = Time of computation

Subbox = [1.0,1.5]x[2,3]
With the simple Monte-Carlo:
0.9960434492716684 = Ians
0.002638618869770148 = Standard Error
65.7294411 = Time of computation

Subbox = [1.0,1.5]x[3,4]
With the simple Monte-Carlo:
1.0010465079190625 = Ians
0.0026383163876978503 = Standard Error
65.641717299 = Time of computation

Subbox = [1.0,1.5]x[4,5]
With the simple Monte-Carlo:
0.9990839683050481 = Ians
0.002628731297879795 = Standard Error
65.8823216 = Time of computation

Subbox = [1.0,1.5]x[5,6]
With the simple Monte-Carlo:
1.0018878755956966 = Ians
0.0026259703224399567 = Standard Error
66.425663601 = Time of computation

Subbox = [1.0,1.5]x[6,7]
With the simple Monte-Carlo:
1.7478729232938397 = Ians
0.002622213063347329 = Standard Error
65.915383499 = Time of computation

Subbox = [1.0,1.5]x[7,8]
With the simple Monte-Carlo:
2.002330890132737 = Ians
0.0026271295912966655 = Standard Error
65.931032 = Time of computation

Subbox = [1.0,1.5]x[8,9]
With the simple Monte-Carlo:
2.0017312684714077 = Ians
0.002622732603094959 = Standard Error
65.730363001 = Time of computation

Subbox = [1.0,1.5]x[9,10]
With the simple Monte-Carlo:
2.0000099564740745 = Ians
0.00261534050476698 = Standard Error
65.438974 = Time of computation

Subbox = [1.5,2.0]x[0,1]
With the simple Monte-Carlo:
1.0038073982295654 = Ians
0.0024674685806088947 = Standard Error
65.4985269 = Time of computation

Subbox = [1.5,2.0]x[1,2]
With the simple Monte-Carlo:
1.0003267951116446 = Ians
0.0024550444467411514 = Standard Error
64.809176 = Time of computation

Subbox = [1.5,2.0]x[2,3]
With the simple Monte-Carlo:
1.00005340036099 = Ians
0.0024468079583856517 = Standard Error
64.716249799 = Time of computation

Subbox = [1.5,2.0]x[3,4]
With the simple Monte-Carlo:
1.0007396352676956 = Ians
0.0024397768936242708 = Standard Error
65.5875119 = Time of computation

Subbox = [1.5,2.0]x[4,5]
With the simple Monte-Carlo:
1.0900934671635567 = Ians
0.0024235713951124478 = Standard Error
65.8575441 = Time of computation

Subbox = [1.5,2.0]x[5,6]
With the simple Monte-Carlo:
1.7373934276300158 = Ians
0.0024302458385546113 = Standard Error
65.313687199 = Time of computation

Subbox = [1.5,2.0]x[6,7]
With the simple Monte-Carlo:
2.1085687746332873 = Ians
0.002423132995241592 = Standard Error
64.984989201 = Time of computation

Subbox = [1.5,2.0]x[7,8]
With the simple Monte-Carlo:
1.9986333661681486 = Ians
0.002419555728656342 = Standard Error
64.996817201 = Time of computation

Subbox = [1.5,2.0]x[8,9]
With the simple Monte-Carlo:
2.0028172511470808 = Ians
0.0024194418133826044 = Standard Error
65.2000719 = Time of computation

Subbox = [1.5,2.0]x[9,10]
With the simple Monte-Carlo:
1.9989175837715518 = Ians
0.0024099903093966236 = Standard Error
64.601035399 = Time of computation

Subbox = [2.0,2.5]x[0,1]
With the simple Monte-Carlo:
1.239574242226472 = Ians
0.0024692229434371497 = Standard Error
65.204969001 = Time of computation

Subbox = [2.0,2.5]x[1,2]
With the simple Monte-Carlo:
1.6270019387328334 = Ians
0.0024607214223924235 = Standard Error
65.689231199 = Time of computation

Subbox = [2.0,2.5]x[2,3]
With the simple Monte-Carlo:
2.0447827681661788 = Ians
0.0024577039786416823 = Standard Error
65.8255261 = Time of computation

Subbox = [2.0,2.5]x[3,4]
With the simple Monte-Carlo:
2.5092356154462467 = Ians
0.002460499582981208 = Standard Error
65.3203147 = Time of computation

Subbox = [2.0,2.5]x[4,5]
With the simple Monte-Carlo:
2.9465845066472904 = Ians
0.0024640979752183225 = Standard Error
64.3456638 = Time of computation

Subbox = [2.0,2.5]x[5,6]
With the simple Monte-Carlo:
3.000993540479677 = Ians
0.002467467099773682 = Standard Error
64.398752 = Time of computation

Subbox = [2.0,2.5]x[6,7]
With the simple Monte-Carlo:
2.248864482016979 = Ians
0.002441436764298013 = Standard Error
64.7627866 = Time of computation

Subbox = [2.0,2.5]x[7,8]
With the simple Monte-Carlo:
1.9972639452758256 = Ians
0.0024254223300966713 = Standard Error
67.7808043 = Time of computation

Subbox = [2.0,2.5]x[8,9]
With the simple Monte-Carlo:
2.002432606363487 = Ians
0.002424421105805044 = Standard Error
67.283961201 = Time of computation

Subbox = [2.0,2.5]x[9,10]
With the simple Monte-Carlo:
1.9978934607051493 = Ians
0.0024121544988388484 = Standard Error
66.0237631 = Time of computation

Subbox = [2.5,3.0]x[0,1]
With the simple Monte-Carlo:
3.0019673980006085 = Ians
0.003226953884589115 = Standard Error
65.949598399 = Time of computation

Subbox = [2.5,3.0]x[1,2]
With the simple Monte-Carlo:
2.996438450595444 = Ians
0.003216785723160081 = Standard Error
66.4024089 = Time of computation

Subbox = [2.5,3.0]x[2,3]
With the simple Monte-Carlo:
3.001191506394412 = Ians
0.003218289432924406 = Standard Error
67.337732401 = Time of computation

Subbox = [2.5,3.0]x[3,4]
With the simple Monte-Carlo:
2.9969573279670296 = Ians
0.0032077550035628974 = Standard Error
67.5842184 = Time of computation

Subbox = [2.5,3.0]x[4,5]
With the simple Monte-Carlo:
2.996608997991329 = Ians
0.003203496199314992 = Standard Error
65.993113401 = Time of computation

Subbox = [2.5,3.0]x[5,6]
With the simple Monte-Carlo:
3.003987800589343 = Ians
0.003211076436885126 = Standard Error
65.8419024 = Time of computation

Subbox = [2.5,3.0]x[6,7]
With the simple Monte-Carlo:
2.2761838350870627 = Ians
0.0031650438664656885 = Standard Error
65.600894899 = Time of computation

Subbox = [2.5,3.0]x[7,8]
With the simple Monte-Carlo:
2.08609150596282 = Ians
0.0031634618047066483 = Standard Error
66.640968399 = Time of computation

Subbox = [2.5,3.0]x[8,9]
With the simple Monte-Carlo:
2.130539628053404 = Ians
0.0031486812071631143 = Standard Error
65.8277665 = Time of computation

Subbox = [2.5,3.0]x[9,10]
With the simple Monte-Carlo:
2.1851738095204043 = Ians
0.0031450364720642856 = Standard Error
65.810147801 = Time of computation

Subbox = [3.0,3.5]x[0,1]
With the simple Monte-Carlo:
4.723428625623352 = Ians
0.005376883175085313 = Standard Error
66.0542899 = Time of computation

Subbox = [3.0,3.5]x[1,2]
With the simple Monte-Carlo:
4.763707753466185 = Ians
0.0053479770551129995 = Standard Error
65.783898599 = Time of computation

Subbox = [3.0,3.5]x[2,3]
With the simple Monte-Carlo:
4.8151022058479365 = Ians
0.005347063415078264 = Standard Error
66.137046401 = Time of computation

Subbox = [3.0,3.5]x[3,4]
With the simple Monte-Carlo:
4.882713037271855 = Ians
0.005386106374722426 = Standard Error
65.5372897 = Time of computation

Subbox = [3.0,3.5]x[4,5]
With the simple Monte-Carlo:
4.922299872303559 = Ians
0.005367278605222867 = Standard Error
65.9027798 = Time of computation

Subbox = [3.0,3.5]x[5,6]
With the simple Monte-Carlo:
4.980541073034036 = Ians
0.005380764617931791 = Standard Error
65.9322817 = Time of computation

Subbox = [3.0,3.5]x[6,7]
With the simple Monte-Carlo:
4.250212877450532 = Ians
0.005350629032863543 = Standard Error
67.4341954 = Time of computation

Subbox = [3.0,3.5]x[7,8]
With the simple Monte-Carlo:
4.003900020210291 = Ians
0.005361453009239169 = Standard Error
66.2864618 = Time of computation

Subbox = [3.0,3.5]x[8,9]
With the simple Monte-Carlo:
3.995220256213367 = Ians
0.0053443324662447896 = Standard Error
74.550557801 = Time of computation

Subbox = [3.0,3.5]x[9,10]
With the simple Monte-Carlo:
3.998933517710548 = Ians
0.005342433872680428 = Standard Error
65.910435599 = Time of computation

Subbox = [3.5,4.0]x[0,1]
With the simple Monte-Carlo:
5.014984956272623 = Ians
0.008949640255831082 = Standard Error
65.823620801 = Time of computation

Subbox = [3.5,4.0]x[1,2]
With the simple Monte-Carlo:
5.004013598783203 = Ians
0.008909921967613402 = Standard Error
65.994103701 = Time of computation

Subbox = [3.5,4.0]x[2,3]
With the simple Monte-Carlo:
5.01099801242691 = Ians
0.008944177216737799 = Standard Error
66.4510516 = Time of computation

Subbox = [3.5,4.0]x[3,4]
With the simple Monte-Carlo:
5.00714301301951 = Ians
0.00892204017163232 = Standard Error
66.5403592 = Time of computation

Subbox = [3.5,4.0]x[4,5]
With the simple Monte-Carlo:
5.015017393992495 = Ians
0.008937625177874222 = Standard Error
66.767908001 = Time of computation

Subbox = [3.5,4.0]x[5,6]
With the simple Monte-Carlo:
5.004736227921822 = Ians
0.008920054806601883 = Standard Error
65.772662999 = Time of computation

Subbox = [3.5,4.0]x[6,7]
With the simple Monte-Carlo:
4.247028517357584 = Ians
0.008855628083337955 = Standard Error
66.157324899 = Time of computation

Subbox = [3.5,4.0]x[7,8]
With the simple Monte-Carlo:
3.998483280052057 = Ians
0.008863381861031328 = Standard Error
66.3256372 = Time of computation

Subbox = [3.5,4.0]x[8,9]
With the simple Monte-Carlo:
3.9957455223574048 = Ians
0.00881780355857697 = Standard Error
66.066956499 = Time of computation

Subbox = [3.5,4.0]x[9,10]
With the simple Monte-Carlo:
4.004636766580441 = Ians
0.008870497325043534 = Standard Error
65.930974601 = Time of computation

Subbox = [4.0,4.5]x[0,1]
With the simple Monte-Carlo:
4.995773911946785 = Ians
0.013885976576503421 = Standard Error
67.443459801 = Time of computation

Subbox = [4.0,4.5]x[1,2]
With the simple Monte-Carlo:
4.976681590701036 = Ians
0.013739284108570452 = Standard Error
66.572933099 = Time of computation

Subbox = [4.0,4.5]x[2,3]
With the simple Monte-Carlo:
4.998816397137433 = Ians
0.01389112332500087 = Standard Error
67.795565601 = Time of computation

Subbox = [4.0,4.5]x[3,4]
With the simple Monte-Carlo:
5.002711169218696 = Ians
0.013908699233508072 = Standard Error
68.118679901 = Time of computation

Subbox = [4.0,4.5]x[4,5]
With the simple Monte-Carlo:
5.01301536401011 = Ians
0.013956030386087916 = Standard Error
67.2934957 = Time of computation

Subbox = [4.0,4.5]x[5,6]
With the simple Monte-Carlo:
5.00393760535643 = Ians
0.013924891736340246 = Standard Error
65.124175 = Time of computation

Subbox = [4.0,4.5]x[6,7]
With the simple Monte-Carlo:
4.258047749344354 = Ians
0.01392454588889124 = Standard Error
68.7682854 = Time of computation

Subbox = [4.0,4.5]x[7,8]
With the simple Monte-Carlo:
4.004441225902264 = Ians
0.013932363553753619 = Standard Error
65.172571899 = Time of computation

Subbox = [4.0,4.5]x[8,9]
With the simple Monte-Carlo:
3.998580898423362 = Ians
0.013867612114168868 = Standard Error
65.7454109 = Time of computation

Subbox = [4.0,4.5]x[9,10]
With the simple Monte-Carlo:
4.01135081814241 = Ians
0.013933144287791173 = Standard Error
65.979236 = Time of computation

Subbox = [4.5,5.0]x[0,1]
With the simple Monte-Carlo:
4.9527860837408655 = Ians
0.02007357874548212 = Standard Error
66.995969101 = Time of computation

Subbox = [4.5,5.0]x[1,2]
With the simple Monte-Carlo:
4.983973801825871 = Ians
0.02052386009212918 = Standard Error
67.8620472 = Time of computation

Subbox = [4.5,5.0]x[2,3]
With the simple Monte-Carlo:
4.9929990826932205 = Ians
0.020514688127955365 = Standard Error
67.063210199 = Time of computation

Subbox = [4.5,5.0]x[3,4]
With the simple Monte-Carlo:
5.02075513566439 = Ians
0.020693624541805653 = Standard Error
67.275764999 = Time of computation

Subbox = [4.5,5.0]x[4,5]
With the simple Monte-Carlo:
5.0055673091322115 = Ians
0.020634113541372606 = Standard Error
66.5915571 = Time of computation

Subbox = [4.5,5.0]x[5,6]
With the simple Monte-Carlo:
4.98933721183469 = Ians
0.02043052919753083 = Standard Error
66.9963904 = Time of computation

Subbox = [4.5,5.0]x[6,7]
With the simple Monte-Carlo:
4.222807716761765 = Ians
0.020233861885235887 = Standard Error
66.2309054 = Time of computation

Subbox = [4.5,5.0]x[7,8]
With the simple Monte-Carlo:
3.9865149866445537 = Ians
0.020447400020759433 = Standard Error
65.8592979 = Time of computation

Subbox = [4.5,5.0]x[8,9]
With the simple Monte-Carlo:
4.016666513850456 = Ians
0.020691223446424583 = Standard Error
66.496677799 = Time of computation

Subbox = [4.5,5.0]x[9,10]
With the simple Monte-Carlo:
4.024092088790899 = Ians
0.020741007834271147 = Standard Error
66.853410799 = Time of computation

Total time = 6614.788276399 =#

# The end of the file.
