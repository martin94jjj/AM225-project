#include <fstream> // Header for file output
#include <vector> // Header for vectors
#include <cmath> // Header for sqrt(), fabs()
#include "omp.h"
int main() {
		// Specify simulation parameters
	for(int scaling = 0 ; scaling<8;scaling++){
		double theta_i = 20.0;
		double theta_v = -20.0;
		double sigma = 100.0;
		double deltaTheta = 0.02;
		double trans_coeff = 4;


		// Calculate other parameters
		double deltaT = deltaTheta / sigma;
		double maxT = 2 * fabs(theta_v - theta_i) / sigma;
		double maxX = 6 * sqrt(maxT);
		int m = (int)( maxT / deltaT ); // number of timesteps
		double k_not = 1*exp(-scaling);//aqds is around 500
		printf("#index = %d, k_not = %g\n",scaling,k_not);

		// Arrays storing peak values
		double* theta_array = new double[m];
		double* flux_array = new double[m];

		//unequally spaced grid
		double omega = 1.001;
		double h = 0.00016;//optimized point
		std::vector<double> X;
		X.push_back(0.0);
		while(X.back() < maxX) // where maxX = 6*sqrt(maxT)
		{
			X.push_back(X.back() + h);
			h *= omega;
		}
		int n = X.size(); //number of spacesteps

		// Calculate Thomas coefficients
		std::vector<double> alpha(n, 0);
		std::vector<double> beta(n, 0);
		std::vector<double> gamma(n, 0);
		std::vector<double> f_theta(m,0);
		for(int i=1; i<n-1; i++)
		{
			double delX_m = X[i] - X[i-1];
			double delX_p = X[i+1] - X[i];
			alpha[i] = -(2* deltaT) / ( delX_m * (delX_m + delX_p) );
			gamma[i] = -(2* deltaT) / ( delX_p * (delX_m + delX_p) );
			beta[i] = 1 - alpha[i] - gamma[i];
		}


		// Create containers
		std::vector<double> g_mod(n-1, 0);
		std::vector<double> C(n, 1.0); // concentration profile

		// Modify gamma coefficients
		// g_mod[0] = 0; // boundary condition


		// BEGIN SIMULATION
		//printf("#K^0=%g\n",k_not);
		double Theta = theta_i;
		for(int k=0; k<m; k++) {
			if ( k < m / 2 ) { 
				Theta -= deltaTheta; 
			} else { 
				Theta += deltaTheta; 
			}
			f_theta[k] = exp(-trans_coeff*Theta);

			// Forward sweep - create modified deltas
			
			alpha[0] = 0;
			beta[0] = 1+ (X[1] - X[0])*f_theta[k]*k_not*(1+exp(Theta));
			gamma[0] = -1;
			C[0] = (X[1] - X[0])*f_theta[k]*k_not*exp(Theta);//delta
			g_mod[0] = gamma[0]/beta[0];

			C[0] = C[0]/beta[0];//delta prime

			//Subce g_mod[0] is no longer constant, needs to upgrate other g_mod accordingly
			for(int i=1; i<n-1; i++) {
				g_mod[i] = gamma[i] / (beta[i] - g_mod[i-1] * alpha[i]);
			}
			//C[0] = 1.0 / (1.0 + exp(-Theta));//Nernst

			for(int i=1; i<n-1; i++) {
				C[i] = ( C[i] - C[i-1]*alpha[i] ) / ( beta[i] - g_mod[i-1] * alpha[i] );
			}
			// Back Substitution
			// C[n-1] = 1.0;
			for(int i=n-2; i>=0; i--) {
				C[i] = C[i] - g_mod[i] * C[i+1];
			}
			
			// Output current
			double flux = -(-C[2] + 4*C[1] -3*C[0]) / (2*(X[1]-X[0]));

			//print results
			printf("%g %g \n",Theta*-0.04-0.15,(flux*-1)/9.23*5);

			theta_array[k] = Theta;
			flux_array[k] = flux;
		}
		printf("\n\n");
	}

	// END SIMULATION
}