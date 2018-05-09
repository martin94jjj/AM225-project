#include <fstream> // Header for file output
#include <vector> // Header for vectors
#include <cmath> // Header for sqrt(), fabs()

int main() {

	// Specify simulation parameters
	double theta_i = 20.0;
	double theta_v = -20.0;
	double sigma = 100.0;
	//double deltaX = 2e-5;
	double deltaTheta = 0.02;

	// Calculate other parameters
	double deltaT = deltaTheta / sigma;
	double maxT = 2 * fabs(theta_v - theta_i) / sigma;
	double maxX = 6 * sqrt(maxT);
	int m = (int)( maxT / deltaT ); // number of timesteps
	//printf("maxX=%g, deltaX=%g, n=%d, (n-1)*deltaX=%g\n", maxX, deltaX, n, (n-1)*deltaX);

	// Arrays storing peak values
	double* theta_array = new double[m];
	double* flux_array = new double[m];

	//unequally spaced grid
	double omega = 1.001;
	double h = 1e-4;
	std::vector<double> X;
	X.push_back(0.0);
	while(X.back() < maxX) // where maxX = 6*sqrt(maxT)
	{
		X.push_back(X.back() + h);
		h *= omega;
	}
	int n = X.size(); //number of spacesteps
	//printf("maxX=%g, deltaX=%g, n=%d, (n-1)*deltaX=%g\n", maxX, deltaX, n, (n-1)*deltaX);

	// Calculate Thomas coefficients
	std::vector<double> alpha(n, 0);
	std::vector<double> beta(n, 0);
	std::vector<double> gamma(n, 0);
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
	for(int i=1; i<n-1; i++) {
		g_mod[i] = gamma[i] / (beta[i] - g_mod[i-1] * alpha[i]);
	}

	// Open file to output CV
	std::ofstream CV("CV_Output_changing_X.txt");

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
			C[i] = ( C[i] - C[i-1]*alpha[i] ) / ( beta[i] - g_mod[i-1] * alpha[i] );
		}
		
		// Back Substitution
		// C[n-1] = 1.0;
		for(int i=n-2; i>=1; i--) {
			C[i] = C[i] - g_mod[i] * C[i+1];
		}
		
		// Output current
		double flux = -(-C[2] + 4*C[1] -3*C[0]) / (2*(X[1]-X[0]));
		CV << Theta << "\t" << flux << "\n";
		theta_array[k] = Theta;
		flux_array[k] = flux;
	}
		double peak_pos = 0;
		double peak_neg = 0;
		double theta_pos = 0;
		double theta_neg = 0;
		double peak_separation = 0;
		for(int i = 0; i < m; i++){
			if(flux_array[i]>peak_pos){
				peak_pos = flux_array[i];
				theta_pos = theta_array[i];
			}
			if(flux_array[i]<peak_neg){
				peak_neg = flux_array[i];
				theta_neg = theta_array[i];
			}
		}

		peak_separation = theta_pos-theta_neg;
		printf("perfecly reversible peak separation %g\n",peak_separation);
// END SIMULATION
}