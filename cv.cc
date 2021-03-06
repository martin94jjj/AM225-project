#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>
#include <cstring>
#include <complex>   // std::complex
#include "omp.h"

typedef std::complex<double> dcomplex;

extern "C" {
    int zgbsv_(const int *n, const int *kl, const int *ku, const int *nrhs,
               dcomplex *ab, const int *ldab, int *ipiv, dcomplex *b,
               const int *ldb, int *info);
}

int main(){
	// Specify simulation parameters
    for(int power = 1; power < 12;power++){
    	double theta_i = 20.0;
    	double theta_v = -20.0;
    	double sigma = pow(2,power);
    	double deltaX = 2e-4;
    	double deltaTheta = 0.02;
    	//Calculate other parameters
    	double deltaT = deltaTheta / sigma;
    	double maxT = 2 * fabs(theta_v - theta_i) / sigma;
    	double maxX = 6*sqrt(maxT);
    	int n_space = (int)( maxX / deltaX ); // number of spacesteps
    	int m_time = (int)( maxT / deltaT ); // number of timesteps
    	// Calculate coefficients for banded solver
    	double lambda = deltaT / (deltaX*deltaX);
    	double alpha = -lambda;
    	double beta = 2.0*lambda + 1.0;
    	double gamma = -lambda;

    	//Use banded solver to solve the concentration profile

        // Size of the matrix
        const int n = n_space;
        // Number of lower and upper diagonals
        const int kl = 1, ku =1;
        // Number of source terms
        const int nrhs = 1;
        // Number of rows in banded storage
        const int ldab = 2*kl+ku+1;
        // Matrix in banded storage
        dcomplex* ab = new dcomplex[ldab*n];
        // Reuse ab
        dcomplex* ab_temp = new dcomplex[ldab*n];
        // Source term
        dcomplex* b = new dcomplex[n];
        // Pivot storage
        int* ipiv   = new int[n];
        // Error flag
        int info;
        // Set the dynamic potential 
    	double Theta = theta_i;

        // Matrix in row-major banded storage
        // Note that it is easier to construct the matrix on a per-diagonal basis
        // if we store things in row-major order
        std::vector<dcomplex> diags(ldab*n, 0.);

        // The first kl rows of banded storage are for internal use by LAPACK
        std::vector<dcomplex>::iterator start = diags.begin()+n*kl;

        // Construct the diagonals
        // Refer to LAPACK's banded storage format for details
        std::fill_n(start+n*kl,         n, dcomplex(beta,0));
        std::fill_n(start+n*(kl-1)+1, n-1, dcomplex(alpha,0));
        std::fill_n(start+n*(kl+1),   n-1, dcomplex(gamma,0));

        //Set the boundary conditions
        diags[n*2] = 1;
        diags[n+1] = 0;

        diags[n*2+n-1] = 1;
        diags[(ldab-1)*n+n-2] = 0;

        //Check the matrix construction
        /*
        for (int i=0; i<ldab; i++){
        	for(int j=0;j<n;j++){
        		printf("%g ",diags[n*i+j].real());
        	}
        	printf("\n");
        }*/

        // Transpose to column-major order for LAPACK
        for (int i=0; i<ldab; ++i)
            for (int j=0; j<n; ++j)
                ab[ldab*j+i] = diags[n*i+j];


        memcpy(ab_temp,ab,ldab*n*sizeof(dcomplex));
        for (int i=0; i<n; ++i) b[i] = dcomplex(1,0);
        //Boundary condition
        b[0] = dcomplex(1.0 / (1.0 + exp(-Theta)),0);
    	b[n-1] = dcomplex(1.0,0);
    //////////////////////////////////////////////////////////////////
        //Perform simulation and count elapsed time
        double t0 = omp_get_wtime();
        printf("#sigma is=%g \n \n \n",sigma);
        for (int i=0 ; i<m_time;i++){
    	    // Print the matrix
    	    /*
    	    for (int j=0; j<n; j++) {
    	        printf(j==0?"A=[":"  [");
    	        for (int i=0; i<n; i++) {
    	            int k=i-j;
    	            dcomplex x = k<-kl||k>ku?0:ab[kl+ku+i*ldab-k];
    	            printf("%2.g ", x.real());
    	        }
    	        puts("]");
    	    }

    	    // Print the source term
    	    printf("\nb=[%2.g ]\n", b->real());
    	    for (int j=1; j<n; j++) printf("  [%2.g ]\n", b[j].real());
    	    */

    		if(i<m_time/2) Theta -= deltaTheta; 
    		else Theta += deltaTheta; 
    	    b[0] = dcomplex(1.0 / (1.0 + exp(-Theta)),0);
    		b[n-1] = dcomplex(1.0,0);

    	    // Call LAPACK
    	    zgbsv_(&n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &n, &info);
    		/*
    	    // Process errors
    	    if (info < 0) {
    	        printf("\nInvalid argument at position %d\n", -info);
    	    } else if (info > 0) {
    	        printf("\nMatrix is singular at position %d\n", info);
    	    } else {
    	        // Print the solution
    	        printf("\nx=[%2.g]\n", b->real());
    	        for (int j=1; j<n; j++) printf("  [%2.g]\n", b[j].real());
    	    }*/
    		double flux = -(-b[2].real()+4*b[1].real()-3*b[0].real())/(2*deltaX);

            //print resutls-> save for gnuplot
    		printf("%g %g\n", Theta, flux);
    		//for(int i =0;i<n;i++) printf("%g \n",b[i].real());
    	    //Boundary condition

    		//reset ab (because zgbsv changes the composition of ab)
        	memcpy(ab,ab_temp,ldab*n*sizeof(dcomplex));

    	}
    	double t1 = omp_get_wtime();
    	//printf("Time elapsed is %g\n",t1-t0);
        delete [] ab;
        delete [] ab_temp;
        delete [] b;
        delete [] ipiv;

    }
}