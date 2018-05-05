#include <cstdio>
#include <cmath>
#include <fftw3.h>

int main(){
	//fftw_plan fftw_plan_r2r_1d(int n, double *in, double *out,
    //                     fftw_r2r_kind kind, unsigned flags);

	// Specify simulation parameters
	double theta_i = 20.0;
	double theta_v = -20.0;
	double sigma = 100.0;
	double deltaX = 2e-4;
	double deltaTheta = 0.02;
	//Calculate other parameters
	double deltaT = deltaTheta / sigma;
	double maxT = 2 * fabs(theta_v - theta_i) / sigma;
	double maxX = 6*sqrt(maxT);
	int n_space = (int)( maxX / deltaX ); // number of spacesteps
	int m_time = (int)( maxT / deltaT ); // number of timesteps
	// Calculate Thomas coefficients
	double lambda = deltaT / (deltaX*deltaX);
	double alpha = -lambda;
	double beta = 2.0*lambda + 1.0;
	double gamma = -lambda;
    // Set the dynamic potential 
	double Theta = theta_i;


	//Use fft solver to solve the concentration profile

    // Size of the matrix
    const int n = n_space;
    //Discretized source term and allocate fft memory
	double* const f = fftw_alloc_real(n);
    //Discretized solution term and allocate fft memory
    double* const v = fftw_alloc_real(n);
    //Store eignvalue array
    double* eigen = new double[n];
    //Make FFT plans
    fftw_plan plan_a = fftw_plan_r2r_1d(n,f,v,FFTW_RODFT00,FFTW_MEASURE);

    //initialize the source term
    for (int i=0; i<n; ++i) f[i] = 1;
    f[0] = 1.0 / (1.0 + exp(-Theta));
	f[n-1] = 1;

	


    fftw_destroy_plan(plan_a);
    delete [] eigen;
    fftw_free(v);
    fftw_free(f);
}