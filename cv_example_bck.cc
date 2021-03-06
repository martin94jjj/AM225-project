#include <fstream> // Header for file output
#include <vector> // Header for vectors
#include <cmath> // Header for sqrt(), fabs()

int main() {

	// Specify simulation parameters
	double theta_i = 20.0;
	double theta_v = -20.0;
	double sigma = 100.0;
	double deltaX = 2e-2;
	double deltaTheta = 0.02;

	// Calculate other parameters
	double deltaT = deltaTheta / sigma;
	double maxT = 2 * fabs(theta_v - theta_i) / sigma;
	double maxX = 6 * sqrt(maxT);
	int n = (int)( maxX / deltaX ); // number of spacesteps
	int m = (int)( maxT / deltaT ); // number of timesteps
	// printf("maxX=%g, deltaX=%g, n=%d, (n-1)*deltaX=%g\n", maxX, deltaX, n, (n-1)*deltaX);
	// printf("maxT=%g, deltaT=%g, m=%d, (m-1)*deltaT=%g\n", maxT, deltaT, m, (m-1)*deltaT);

	// Calculate Thomas coefficients
	double lambda = deltaT / (deltaX*deltaX);
	double alpha = -lambda;
	double beta = 2.0*lambda + 1.0;
	double gamma = -lambda;

	// printf("lambda=%g\n", lambda);

	// Create containers
	std::vector<double> g_mod(n-1, 0);
	std::vector<double> C(n, 1.0); // concentration profile

	// Modify gamma coefficients
	// g_mod[0] = 0; // boundary condition
	for(int i=1; i<n-1; i++) {
		g_mod[i] = gamma / (beta - g_mod[i-1] * alpha);
	}

	// Open file to output CV
	std::ofstream CV("CV_Output1.txt");

	// BEGIN SIMULATION
	double Theta = theta_i;
	for(int k=0; k<m; k++) {
		if ( k < m / 2 ) { 
			Theta -= deltaTheta; 
		} else { 
			Theta += deltaTheta; 
		}
	
		// Forward sweep - create modified deltas
		C[0] = 1.0 / (1.0 + exp(-Theta));
		for(int i=1; i<n-1; i++) {
			C[i] = ( C[i] - C[i-1]*alpha ) / ( beta - g_mod[i-1] * alpha );
		}
		
		// Back Substitution
		// C[n-1] = 1.0;
		for(int i=n-2; i>=1; i--) {
			C[i] = C[i] - g_mod[i] * C[i+1];
		}
		
		// Output current
		double flux = -(-C[2] + 4*C[1] -3*C[0]) / (2*deltaX);
		CV << Theta << "\t" << flux << "\t" << C[0] << "\n";

		// Output concentration profile at maxT/4, maxT/2, maxT*3/4, maxT
		// if ((k==m/4)|(k==m/2)|(k==m*3/4)|(k==m-1)) {
		// 	for(int i=0; i<n; i++) {
		// 		printf("%g %g\n", i*deltaX, C[i]);
		// 	}
		// 	printf("\n\n");
		// }
		
	}
// END SIMULATION
}