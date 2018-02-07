#ifndef blas_api_h
#define blas_api_h 1


extern "C" {

// sgemm - BLAS routine - single precision !!!
//---------------------------------------------
// used to calculate general matrix-matrix multiplications of the form
// C = alpha * A * B + beta * C
static inline int sgemm(char transa, char transb, int m, int n, int k, float alpha, float* a,  int lda, float* b, int ldb, float beta, float* c, int ldc) {
  extern int sgemm_(char* transa, char* transb, int* m, int* n, int* k, float* alpha, float* a, int* lda, float* b, int* ldb, float* beta, float* c, int* ldc);
  return sgemm_(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}


// dgemm - BLAS routine
//---------------------
// used to calculate general matrix-matrix multiplications of the form
// C = alpha * A * B + beta * C
static inline int dgemm(char transa, char transb, int m, int n, int k, double alpha, double* a, int lda, double* b, int ldb, double beta, double* c, int ldc) {
  extern int dgemm_(char* transa, char* transb, int* m, int* n, int* k, double* alpha, double* a, int* lda, double* b, int* ldb, double* beta, double* c, int* ldc);
  return dgemm_(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}

}


void multiplyMatrix(bool transA, bool transB, int rowsA, int colsB, int k, double* mA, double* mB, double* mC)
{

// rowsA und colsB sind die Werte NACH einem möglichen Transponieren!!!
    
	dgemm(transA? 'T' : 'N', transB? 'T' : 'N', rowsA, colsB, k, 1., mA, transA? k : rowsA, mB, transB? colsB : k, 0., mC, rowsA);
}   

void multiplyMatrix(bool transA, bool transB, int rowsA, int colsB, int k, float* mA, float* mB, float* mC)
{

// rowsA und colsB sind die Werte NACH einem möglichen Transponieren!!!
    
	sgemm(transA? 'T' : 'N', transB? 'T' : 'N', rowsA, colsB, k, 1., mA, transA? k : rowsA, mB, transB? colsB : k, 0., mC, rowsA);
}   


#endif
