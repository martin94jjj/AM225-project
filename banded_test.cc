#include <cstdio>
#include <vector>
#include <algorithm> // std::fill_n
#include <complex>   // std::complex

typedef std::complex<double> dcomplex;

extern "C" {
    int zgbsv_(const int *n, const int *kl, const int *ku, const int *nrhs,
               dcomplex *ab, const int *ldab, int *ipiv, dcomplex *b,
               const int *ldb, int *info);
}

int main()
{
    // Size of the matrix
    const int n = 6;
    // Number of lower and upper diagonals
    const int kl = 1, ku =1;
    // Number of source terms
    const int nrhs = 1;
    // Number of rows in banded storage
    const int ldab = 2*kl+ku+1;
    // Matrix in banded storage
    dcomplex* ab = new dcomplex[ldab*n];
    // Source term
    dcomplex* b = new dcomplex[n];
    // Pivot storage
    int* ipiv   = new int[n];
    // Error flag
    int info;

    // Matrix in row-major banded storage
    // Note that it is easier to construct the matrix on a per-diagonal basis
    // if we store things in row-major order
    std::vector<dcomplex> diags(ldab*n, 0.);

    // The first kl rows of banded storage are for internal use by LAPACK
    std::vector<dcomplex>::iterator start = diags.begin()+n*kl;

    // Construct the diagonals
    // Refer to LAPACK's banded storage format for details
    std::fill_n(start+n*kl,         n, dcomplex(-2,-1));
    std::fill_n(start+n*(kl-1)+1, n-1, dcomplex(0,1));
    std::fill_n(start+n*(kl+1),   n-1, dcomplex(0,1));

    // Transpose to column-major order for LAPACK
    for (int i=0; i<ldab; ++i)
        for (int j=0; j<n; ++j)
            ab[ldab*j+i] = diags[n*i+j];

    // Fill the source term
    for (int i=0; i<n; ++i) b[i] = dcomplex(1,1);

    // Print the matrix
    for (int j=0; j<n; j++) {
        printf(j==0?"A=[":"  [");
        for (int i=0; i<n; i++) {
            int k=i-j;
            dcomplex x = k<-kl||k>ku?0:ab[kl+ku+i*ldab-k];
            printf("%2.g%+2.gi ", x.real(), x.imag());
        }
        puts("]");
    }

    // Print the source term
    printf("\nb=[%2.g%+2.gi ]\n", b->real(), b->imag());
    for (int j=1; j<n; j++) printf("  [%2.g%+2.gi ]\n", b[j].real(), b[j].imag());

    // Call LAPACK
    zgbsv_(&n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &n, &info);

    // Process errors
    if (info < 0) {
        printf("\nInvalid argument at position %d\n", -info);
    } else if (info > 0) {
        printf("\nMatrix is singular at position %d\n", info);
    } else {
        // Print the solution
        printf("\nx=[%2.g%+2.gi ]\n", b->real(), b->imag());
        for (int j=1; j<n; j++) printf("  [%2.g%+2.gi ]\n", b[j].real(), b[j].imag());
    }

    delete [] ab;
    delete [] b;
    delete [] ipiv;
}