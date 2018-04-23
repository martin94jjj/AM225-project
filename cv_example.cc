#include <fstream> // Header for file output
#include <vector> // Header for vectors
#include <cmath> // Header for sqrt(), fabs()
#include "omp.h"
int main()
	{
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
	int n = (int)( maxX / deltaX ); // number of spacesteps
	int m = (int)( maxT / deltaT ); // number of timesteps
	// Calculate Thomas coefficients
	double lambda = deltaT / (deltaX*deltaX);
	double alpha = -lambda;
	double beta = 2.0*lambda + 1.0;
	double gamma = -lambda;

	printf("%g %g %g \n",alpha,beta,gamma);
	// Create containers
	std::vector<double> g_mod(n, 0);
	std::vector<double> C(n, 1.0); // concentration profile
	// Modify gamma coefficients
	g_mod[0] = 0; // boundary condition
	for(int i=1; i<n-1; i++) {
		g_mod[i] = gamma / (beta - g_mod[i-1] * alpha);
	}
	// Open file to output CV
	std::ofstream CV("CV_Output.txt");

	// BEGIN SIMULATION
	double Theta = theta_i;
	double t0 = omp_get_wtime();
	for(int k=0; k<m; k++)
	{
		if(k<m/2) Theta -= deltaTheta; 
		else Theta += deltaTheta; 
		// Forward sweep - create modified deltas
		C[0] = 1.0 / (1.0 + exp(-Theta));
		for(int i=1; i<n-1; i++) {
		C[i] = (C[i]-C[i-1]*alpha)/(beta- g_mod[i-1]*alpha);
		}
		// Back Substitution
		C[n-1] = 1.0;
		for(int i=n-2; i>=0; i--) 
			C[i] = C[i] - g_mod[i] * C[i+1];
		//Output current
		double flux = -(-C[2] + 4*C[1] -3*C[0]) / (2*deltaX);
		//printf("%g %g %g %g\n", Theta, C[2],C[1],C[0]);
		CV << Theta << "\t" << flux << "\n";
	}
	double t1 = omp_get_wtime();
	printf("Time elapsed is %g \n", t1-t0);

	// END SIMULATION
}