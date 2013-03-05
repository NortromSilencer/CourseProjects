#include <fftw3.h>
#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include <string>

using namespace std;
const int T  = 50000;
const double k = 0.6;
double a[T];
double Sxx(double omega);
void noiseGenerator();
double randomGenerator(double rmin, double rmax);
bool saveArray( const double* pdata, size_t length, const string& file_path );

int main()
{
	double sum = 0;
	double omega = 0.3;
	for (int i=0;i<1000;i++)
	{
		noiseGenerator();
		sum+=Sxx(omega);
		cout<<i<<endl;
	}
	cout<<sum/1000.0<<endl;
	cout<<1/12.0/(1+k*k-2*k*cos(omega))<<endl;
	cout<<"Press any key to continue."<<endl;
	getchar();
	return 0;
}
void noiseGenerator()
{

	a[0] = 0;
	for (int i=1;i<T;i++)
	{
		double Rt =randomGenerator(-0.5,0.5);
		a[i] = k*a[i-1]+Rt;
	}
	return;
}
double randomGenerator(double rmin, double rmax)
{
	random_device rd;
	if (rd.entropy() == 0)
	{
		cout<<"real random number device does not exist on this system."<<endl;
		return 0;
	}
	double width = rmax-rmin;
	if (width <= 0)
	{
		cout<<"Error: rmin should not >= rmax"<<endl;
		return 0;
	}
	double r  = rd()/(double)rd.max()*width+rmin; //rmin to rmax uniform distribution;  
	return r;
}
double Sxx(double omega)
{
	double result;

	double Re = 0;
	double Im = 0;
	for(int i=0;i<T;i++)
	{
		Re += a[i]*(cos(i*omega));
		Im += a[i]*(sin(i*omega));
	}
	result = 1.0/(double)T*(Re*Re+Im*Im);
	return result;
}
