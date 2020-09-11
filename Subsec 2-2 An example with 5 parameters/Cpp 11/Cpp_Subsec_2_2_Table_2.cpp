/* 
The Monte Carlo integration with C++. Table 2 of the paper.
*/
# include <iostream>
# include <cmath>
# include <chrono>
# include <random>
using namespace std;
using namespace std::chrono;
std::default_random_engine & global_urng()
{
	static std::default_random_engine u{};
	return u;
}
void randomize()
{
	static std::random_device rd{};
	global_urng().seed( rd() );
}
double samplor(double a, double b)
{
  static std::uniform_real_distribution<> d{};
  using parm_t=decltype(d)::param_type;
  return d( global_urng(), parm_t{a,b} );
}
int chi(double a, double b, double x)
{
	if (x<a)
	  return 0;
	else if (x>b)
	  return 0;
	else
	  return 1;
}
double alok1(double x1, double x2, double k2, double k3, double k4)
{
	return (2*k2*pow(x1,3)+k3*x1*pow(x2,2)-2*k4*pow(x2,3))/(pow(x1,2)*x2);
}
double aloT(double x1, double x2, double k2, double k3, double k4)
{
	return x1+x2;
}
double detJf(double x1, double x2, double k1, double k2, double k3, double k4)
{
	return -(6*k2+k1)*pow(x1,2)-(6*k4+k3)*pow(x2,2)+2*(k1+k3)*x1*x2;
}
double J(double x1, double x2, double k2, double k3, double k4)
{
	return abs(detJf(x1,x2,alok1(x1,x2,k2,k3,k4),k2,k3,k4));
}
double IntegrandWithoutCOEFF(double x1, double x2, double k2, double k3, double k4, double ak1, double bk1, double aT, double bT)
{
	return J(x1,x2,k2,k3,k4)*chi(ak1,bk1,alok1(x1,x2,k2,k3,k4))*chi(aT,bT,aloT(x1,x2,k2,k3,k4))/(pow(x1,2)*x2);
}
double *sumo_t (double kk1[], double kk2[], double kk3[], double kk4[], double TT[], int NN)
{
	double sumor=0;
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	double COEFF=pow(TT[1],2)/((TT[1]-TT[0])*(kk1[1]-kk1[0]));
	for (int i=0;i<NN;i+=1)
	{
		double x1=samplor(0,TT[1]);
		double x2=samplor(0,TT[1]);
		double k2=samplor(kk2[0],kk2[1]);
		double k3=samplor(kk3[0],kk3[1]);
		double k4=samplor(kk4[0],kk4[1]);
		sumor+=IntegrandWithoutCOEFF(x1,x2,k2,k3,k4,kk1[0],kk1[1],TT[0],TT[1]);
	}
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	auto st = duration_cast<microseconds>( t2 - t1 ).count();
	double result[2]={};
	result[0]=COEFF*sumor/NN;
	result[1]=st;
	double *resulto=result;
	return resulto;
}
double *sumo_with_SE_t(double kk1[], double kk2[], double kk3[], double kk4[], double TT[], int NN)
{
	double Ians=0,x1,x2,k2,k3,k4,delto,S=0;
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	double COEFF=pow(TT[1],2)/((TT[1]-TT[0])*(kk1[1]-kk1[0]));
	x1=samplor(0,TT[1]);
	x2=samplor(0,TT[1]);
	k2=samplor(kk2[0],kk2[1]);
	k3=samplor(kk3[0],kk3[1]);
	k4=samplor(kk4[0],kk4[1]);
	Ians+=IntegrandWithoutCOEFF(x1,x2,k2,k3,k4,kk1[0],kk1[1],TT[0],TT[1]);
	Ians*=COEFF;
	for (int i=1;i<NN;i+=1)
	{
		x1=samplor(0,TT[1]);
		x2=samplor(0,TT[1]);
		k2=samplor(kk2[0],kk2[1]);
		k3=samplor(kk3[0],kk3[1]);
		k4=samplor(kk4[0],kk4[1]);
		delto=COEFF*IntegrandWithoutCOEFF(x1,x2,k2,k3,k4,kk1[0],kk1[1],TT[0],TT[1])-Ians;
		Ians+=delto/(i+1);
		S+=pow(delto,2)*i/(i+1);
	}
	double Estandardo=S/NN;
	Estandardo/=NN-1;
	Estandardo=sqrt(Estandardo);
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	auto st = duration_cast<microseconds>( t2 - t1 ).count();
	double result[3]={};
	result[0]=Ians;
	result[1]=Estandardo;
	result[2]=st;
	double *resulto=result;
	return resulto;
}
double *sumo_antithetic_with_SE_t(double kk1[], double kk2[], double kk3[], double kk4[], double TT[], int NN)
{
	double Ians=0,x1,x2,k2,k3,k4,delto,S=0;
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	double centero[5]={(0+TT[1])/2,(0+TT[1])/2,(kk2[0]+kk2[1])/2,(kk3[0]+kk3[1])/2,(kk4[0]+kk4[1])/2};
	double COEFF=pow(TT[1],2)/((TT[1]-TT[0])*(kk1[1]-kk1[0]));
	x1=samplor(0,TT[1]);
	x2=samplor(0,TT[1]);
	k2=samplor(kk2[0],kk2[1]);
	k3=samplor(kk3[0],kk3[1]);
	k4=samplor(kk4[0],kk4[1]);
	Ians+=IntegrandWithoutCOEFF(x1,x2,k2,k3,k4,kk1[0],kk1[1],TT[0],TT[1]);
	Ians*=COEFF;
	x1=2*centero[0]-x1;
	x2=2*centero[1]-x2;
	k2=2*centero[2]-k2;
	k3=2*centero[3]-k3;
	k4=2*centero[4]-k4;
	delto=COEFF*IntegrandWithoutCOEFF(x1,x2,k2,k3,k4,kk1[0],kk1[1],TT[0],TT[1])-Ians;
	Ians+=delto/2;
	S+=pow(delto,2)/2;
	for (int i=1;i<floor(NN/2);i+=1)
	{
		x1=samplor(0,TT[1]);
		x2=samplor(0,TT[1]);
		k2=samplor(kk2[0],kk2[1]);
		k3=samplor(kk3[0],kk3[1]);
		k4=samplor(kk4[0],kk4[1]);
		delto=COEFF*IntegrandWithoutCOEFF(x1,x2,k2,k3,k4,kk1[0],kk1[1],TT[0],TT[1])-Ians;
		Ians+=delto/(2*i+1);
		S+=pow(delto,2)*2*i/(2*i+1);
		x1=2*centero[0]-x1;
		x2=2*centero[1]-x2;
		k2=2*centero[2]-k2;
		k3=2*centero[3]-k3;
		k4=2*centero[4]-k4;
		delto=COEFF*IntegrandWithoutCOEFF(x1,x2,k2,k3,k4,kk1[0],kk1[1],TT[0],TT[1])-Ians;
		Ians+=delto/(2*i+2);
		S+=pow(delto,2)*(2*i+1)/(2*i+2);
	}
	double Estandardo=S/NN;
	Estandardo/=NN-1;
	Estandardo=sqrt(Estandardo);
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	auto st = duration_cast<microseconds>( t2 - t1 ).count();
	double result[3]={};
	result[0]=Ians;
	result[1]=Estandardo;
	result[2]=st;
	double *resulto=result;
	return resulto;
}

int main()
{
    double kk1[2]={0,100};
	double kk2[2]={0,2};
	double kk3[2]={0,200};
	double kk4[2]={0,100};
	double TT[2]={0,2};
	for (int i=1;i<9;i+=1)
	{
		auto y8=sumo_with_SE_t(kk1,kk2,kk3,kk4,TT,pow(10,i));
	    double y8o[3]={y8[0],y8[1],y8[2]};
		auto y9=sumo_antithetic_with_SE_t(kk1,kk2,kk3,kk4,TT,pow(10,i));
	    double y9o[3]={y9[0],y9[1],y9[2]};
	    cout << "\nN=" << pow(10,i) << ".\n" \
		<< "With simple Monte Carlo: \n" \
		<< " " << y8o[0] << " = Ians,\n" \
		<< " " << y8o[1] << " = Standard Error,\n" \
		<< " " << y8o[2]/pow(10,6) << " seconds = Time of the computation.\n" \
		<< "With antithetic Monte Carlo: \n" \
		<< " " << y9o[0] << " = Ians,\n" \
		<< " " << y9o[1] << " = Standard Error,\n" \
		<< " " << y9o[2]/pow(10,6) << " seconds = Time of the computation.\n";
		if (y9o[2]*y9o[1]!=0)
		{
			cout << "The efficiency of the antithetic to the simple method= " << (y8o[2]*y8o[1])/(y9o[2]*y9o[1]) << ".\n";
		}
	}
	
	getchar();
	return 0;
}
