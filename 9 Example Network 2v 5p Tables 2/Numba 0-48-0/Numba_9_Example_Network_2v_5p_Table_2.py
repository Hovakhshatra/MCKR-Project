# The Monte Carlo integration with Numba. Table 2 of the paper.
#
from numba import jit, njit, prange
import math
import random
import time
@njit(fastmath=True)
def samplor(a,b):
    return(random.uniform(a,b))
@njit(fastmath=True)
def chi(a,b,x):
    if x<a:
        return(0)
    elif x>b:
        return(0)
    else:
        return(1)
@njit(fastmath=True)
def alok1(x1,x2,k2,k3,k4):
    return((2*k2*(x1**3)+k3*x1*(x2**2)-2*k4*(x2**3))/((x1**2)*x2))
@njit(fastmath=True)
def aloT(x1,x2,k2,k3,k4):
    return(x1+x2)
@njit(fastmath=True)
def detJf(x1,x2,k1,k2,k3,k4):
    return(-(6*k2+k1)*(x1**2)-(6*k4+k3)*(x2**2)+2*(k1+k3)*x1*x2)
@njit(fastmath=True)
def J(x1,x2,k2,k3,k4):
    return(math.fabs(detJf(x1,x2,alok1(x1,x2,k2,k3,k4),k2,k3,k4)))
@njit(fastmath=True)
def IntegrandWithoutCOEFF(x1,x2,k2,k3,k4,ak1,bk1,aT,bT):
    return(J(x1,x2,k2,k3,k4)*chi(ak1,bk1,alok1(x1,x2,k2,k3,k4))*\
           chi(aT,bT,aloT(x1,x2,k2,k3,k4))/((x1**2)*x2))
@njit(fastmath=True)
def sumo(kk1,kk2,kk3,kk4,TT,NN):
    sumor=0
    COEFF=((TT[1]**2)/((TT[1]-TT[0])*(kk1[1]-kk1[0])))
    for i in range(NN):
        x1=samplor(0,TT[1])
        x2=samplor(0,TT[1])
        k2=samplor(kk2[0],kk2[1])
        k3=samplor(kk3[0],kk3[1])
        k4=samplor(kk4[0],kk4[1])
        sumor+=IntegrandWithoutCOEFF(x1,x2,k2,k3,k4,kk1[0],kk1[1],TT[0],TT[1])
    return(COEFF*sumor/NN)
@njit(fastmath=True)
def sumo_with_SE(kk1,kk2,kk3,kk4,TT,NN):
    S=0
    sumor=0
    COEFF=((TT[1]**2)/((TT[1]-TT[0])*(kk1[1]-kk1[0])))
    x1=samplor(0,TT[1])
    x2=samplor(0,TT[1])
    k2=samplor(kk2[0],kk2[1])
    k3=samplor(kk3[0],kk3[1])
    k4=samplor(kk4[0],kk4[1])
    sumor+=IntegrandWithoutCOEFF(x1,x2,k2,k3,k4,kk1[0],kk1[1],TT[0],TT[1])
    Ians=COEFF*sumor
    for i in range(1,NN):
        x1=samplor(0,TT[1])
        x2=samplor(0,TT[1])
        k2=samplor(kk2[0],kk2[1])
        k3=samplor(kk3[0],kk3[1])
        k4=samplor(kk4[0],kk4[1])
        delto=COEFF*\
               IntegrandWithoutCOEFF(x1,x2,k2,k3,k4,kk1[0],kk1[1],TT[0],TT[1])-\
               Ians
        Ians+=delto/(i+1)
        S+=(delto**2)*i/(i+1)
    Estandardo=math.sqrt(S/(NN*(NN-1)))
    return(Ians,Estandardo)
@njit(fastmath=True)
def sumo_antithetic(kk1,kk2,kk3,kk4,TT,NN):
    sumor=0
    centero=[(0+TT[1])/2,(0+TT[1])/2,(kk2[0]+kk2[1])/2,(kk3[0]+kk3[1])/2,\
             (kk4[0]+kk4[1])/2]
    COEFF=((TT[1]**2)/((TT[1]-TT[0])*(kk1[1]-kk1[0])))
    for i in range(math.floor(NN/2)):
        x1=samplor(0,TT[1])
        x2=samplor(0,TT[1])
        k2=samplor(kk2[0],kk2[1])
        k3=samplor(kk3[0],kk3[1])
        k4=samplor(kk4[0],kk4[1])
        sumor+=IntegrandWithoutCOEFF(x1,x2,k2,k3,k4,kk1[0],kk1[1],TT[0],TT[1])
        x1=2*centero[0]-x1
        x2=2*centero[1]-x2
        k2=2*centero[2]-k2
        k3=2*centero[3]-k3
        k4=2*centero[4]-k4
        sumor+=IntegrandWithoutCOEFF(x1,x2,k2,k3,k4,kk1[0],kk1[1],TT[0],TT[1])
    return(COEFF*sumor/NN)
@njit(fastmath=True)
def sumo_antithetic_with_SE(kk1,kk2,kk3,kk4,TT,NN):
    S=0
    sumor=0
    centero=[(0+TT[1])/2,(0+TT[1])/2,(kk2[0]+kk2[1])/2,(kk3[0]+kk3[1])/2,\
             (kk4[0]+kk4[1])/2]
    COEFF=((TT[1]**2)/((TT[1]-TT[0])*(kk1[1]-kk1[0])))
    x1=samplor(0,TT[1])
    x2=samplor(0,TT[1])
    k2=samplor(kk2[0],kk2[1])
    k3=samplor(kk3[0],kk3[1])
    k4=samplor(kk4[0],kk4[1])
    sumor+=IntegrandWithoutCOEFF(x1,x2,k2,k3,k4,kk1[0],kk1[1],TT[0],TT[1])
    Ians=COEFF*sumor
    x1=2*centero[0]-x1
    x2=2*centero[1]-x2
    k2=2*centero[2]-k2
    k3=2*centero[3]-k3
    k4=2*centero[4]-k4
    delto=COEFF*IntegrandWithoutCOEFF(x1,x2,k2,k3,k4,kk1[0],kk1[1],TT[0],TT[1])\
           -Ians
    Ians+=delto/2
    S+=(delto**2)*1/2
    for i in range(1,math.floor(NN/2)):
        x1=samplor(0,TT[1])
        x2=samplor(0,TT[1])
        k2=samplor(kk2[0],kk2[1])
        k3=samplor(kk3[0],kk3[1])
        k4=samplor(kk4[0],kk4[1])
        delto=COEFF*\
               IntegrandWithoutCOEFF(x1,x2,k2,k3,k4,kk1[0],kk1[1],TT[0],TT[1])\
               -Ians
        Ians+=delto/(2*i+1)
        S+=(delto**2)*2*i/(2*i+1)
        x1=2*centero[0]-x1
        x2=2*centero[1]-x2
        k2=2*centero[2]-k2
        k3=2*centero[3]-k3
        k4=2*centero[4]-k4
        delto=COEFF*\
               IntegrandWithoutCOEFF(x1,x2,k2,k3,k4,kk1[0],kk1[1],TT[0],TT[1])-\
               Ians
        Ians+=delto/(2*i+2)
        S+=(delto**2)*(2*i+1)/(2*i+2)
    Estandardo=math.sqrt(S/(NN*(NN-1)))
    return(Ians,Estandardo)
# # # # #
st=time.time()
ans=sumo((0,100),(0,2),(0,200),(0,100),(0,2),1000)
st=time.time()-st
print('Warm-up: just compiling for the first time. N=',1000,', sumo='\
      ,ans,',and st=',st)
st=time.time()
ans=sumo_with_SE((0,100),(0,2),(0,200),(0,100),(0,2),1000)
st=time.time()-st
print('Warm-up: just compiling for the first time. N=',1000,', sumo_with_SE='\
      ,ans,',and st=',st)
st=time.time()
ans=sumo_antithetic((0,100),(0,2),(0,200),(0,100),(0,2),1000)
st=time.time()-st
print('Warm-up: just compiling for the first time. N=',1000,', sumo_antithetic='\
      ,ans,',and st=',st)
st=time.time()
ans=sumo_antithetic_with_SE((0,100),(0,2),(0,200),(0,100),(0,2),1000)
st=time.time()-st
print('Warm-up: just compiling for the first time. N=',1000,', sumo_antithetic_with_SE='\
      ,ans,',and st=',st)
print('\n Now the main computation.\n')
for n in range(1,10):
    st1=time.time()
    X1=sumo_with_SE((0,100),(0,2),(0,200),(0,100),(0,2),10**n)
    st1=time.time()-st1
    st2=time.time()
    X2=sumo_antithetic_with_SE((0,100),(0,2),(0,200),(0,100),(0,2),10**n)
    st2=time.time()-st2
    print('N=',10**n,'\n')
    print('With simple Monte Carlo: \n',X1[0],' = Ians\n',X1[1],' = Standard',
          'Error\n',st1,' = Time of computation\n')
    print('With antithetic Monte Carlo: \n',X2[0],' = Ians\n',X2[1],' =',
          'Standard Error\n',st2,' = Time of computation\n')
    if (st2*X2[1])!=0:
        print('The efficiency of the antithetic to the simple method=',
              (st1*X1[1])/(st2*X2[1]),'\n')
