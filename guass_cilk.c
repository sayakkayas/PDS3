/*
 * Original author:  Sandhya Dwarkadas, 2002 and before.
 * Modified by Grant Farmer, 2003 and Kai Shen, 2010.
 * Minor cleanup by Michael Scott, 2017.
Modified by Sayak Chakraborti ,March 2017 for parallelism using cilk
Compile-gcc guass_cilk.c -fcilkplus -lcilkrts -O3 -o cilk -w
Run-./cilk matrix_2000.dat
 */
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/time.h>
#include <math.h>
#include <assert.h>
//#include <omp.h>
#include <cilk/cilk.h>
/* #define DEBUG */

#define SWAP(a, b)      {double tmp = a; a = b; b = tmp;}
#define SWAPPTR(a,b)    {void* tmp = a; a = b; b = tmp;}
//#define CHUNK 100

/*
 * The code in this file takes as input a matrix of the following format:
 *    line 1: rows cols rows*cols
 *      (we require that rows == cols)
 *    last line: 0 0 0
 *    intermediate lines: i j v
 *      indicating that matrix[i-1][j-1] == v
 *
 * We then choose a right-hand side R such that the vector [1 2 3 ... rows]
 * will be a solution to the equation
 *      matrix * X = R
 *
 * Finally, we solve the equation and verify that we got the expected
 * solution.
 */

double **matrix, *X, *R;

double *X__;        // pre-set solution

// Initialize the matrix.
//
int initMatrix(const char *fname) 
{
  FILE *file;
  int l1, l2, l3;
  double d;
  int nsize;
  int i, j;
  double *tmp;
  char buffer[1024];

  if ((file = fopen(fname, "r")) == NULL)        //open matrix file
  {
    fprintf(stderr, "Matrix file open error\n");        
    exit(-1);
  }

  // Parse the first line to get the matrix size:
  fgets(buffer, 1024, file);
  assert(sscanf(buffer, "%d %d %d", &l1, &l2, &l3) == 3);
  nsize = l1;
  assert(l2 == l1 && l3 == l1 * l2);
  #ifdef DEBUG
  fprintf(stdout, "matrix size is %d\n", nsize);
  #endif

  // Initialize the space and set all elements to zero:
  matrix = (double**) malloc(nsize * sizeof(double*));
  assert(matrix != NULL);
  tmp = (double*) malloc(nsize * nsize * sizeof(double));
  assert(tmp != NULL);
  for (i = 0; i < nsize; i++) 
  {
    matrix[i] = tmp;
    tmp += nsize;
  }
  
  //initialize the matrix with zeros
  
  for (i = 0; i < nsize; i++) 
  {
    for (j = 0; j < nsize; j++) 
    {
      matrix[i][j] = 0.0;
    }
  }

  // Parse the rest of the input file to fill the matrix:
  for (;;) 
  {
    fgets(buffer, 1024, file);
    assert(sscanf(buffer, "%d %d %lf", &l1, &l2, &d) == 3);
    if (l1 == 0) 
    {
      assert(l2 == 0 && d == 0);
      break;
    }
    assert(0 < l1 && l1 <= nsize);
    assert(0 < l2 && l1 <= nsize);

    matrix[l1-1][l2-1] = d;

#ifdef DEBUG
    fprintf(stdout, "row %d column %d of matrix is %e\n",
        l1-1, l2-1, matrix[l1-1][l2-1]);
#endif

  }

  fclose(file);
  return nsize;
}

// Initialize the right-hand-side following the pre-set solution.
//
void initRHS(int nsize)
 {
  int i, j;

  X__ = (double*) malloc(nsize * sizeof(double));
  assert(X__ != NULL);
  
  for (i = 0; i < nsize; i++) 
  {
    X__[i] = i+1;
  }

  R = (double*) malloc(nsize * sizeof(double));
  assert(R != NULL);
  for (i = 0; i < nsize; i++) 
  {
    R[i] = 0.0;
    for (j = 0; j < nsize; j++) 
    {
      R[i] += matrix[i][j] * X__[j];  //as stated R=matrix*X
    }
  }
}

// Initialize the results.
//
void initResult(int nsize) 
{
  int i;

  X = (double*) malloc(nsize * sizeof(double));
  assert(X != NULL);
  for (i = 0; i < nsize; i++) 
  {
    X[i] = 0.0;
  }
}

// Get the pivot - make sure the next row is the one whose initial value
// has the largest absolute value.
//
void getPivot(int nsize, int currow)
 {
  int i, pivotrow;

  pivotrow = currow;
  for (i = currow+1; i < nsize; i++) 
  {
    if (fabs(matrix[i][currow]) > fabs(matrix[pivotrow][currow])) 
    {
      pivotrow = i;
    }
  }

  if (fabs(matrix[pivotrow][currow]) == 0.0) 
  {
    fprintf(stderr, "The matrix is singular\n");
    exit(-1);
  }

  if (pivotrow != currow) 
  {
    #ifdef DEBUG
    fprintf(stdout, "pivot row at step %5d is %5d\n", currow, pivotrow);
    #endif
    SWAPPTR(matrix[pivotrow], matrix[currow]);
    SWAP(R[pivotrow], R[currow]);
  }
}

// For all the rows, get the pivot and eliminate all rows and columns
// for that particular pivot row.
//
void computeGauss(int nsize,int num) 
{
 
   int i;
  //
  int num_threads,tid;
  double pivotval;
 
   __cilkrts_set_param("nworkers","4");        //set number of threads

   printf("Number of threads %d\n",__cilkrts_get_nworkers());  //print number of threads

    for (i = 0; i < nsize; i++) 
    {
      getPivot(nsize, i);
    
      // Scale the main row:
      pivotval = matrix[i][i];   //got pivot value in the main diagonal
   
     //if pivot value not equal to 1 then divide the entire row by it to make it 1
        if (pivotval != 1.0) 
          {
              matrix[i][i] = 1.0;
              int j;
              for (j = i + 1; j < nsize; j++) 
              {
                  matrix[i][j] /= pivotval;
              }
            
              R[i] /= pivotval;
          }

        // Factorize the rest of the matrix:
 
        
                                      
        #pragma grainsize=100           //set grain size
         int j,k;//double piv;           //private variable
        cilk_for (j = i + 1; j < nsize; j++) 
        {
            pivotval = matrix[j][i];

            matrix[j][i] = 0.0;
      
            for (k = i + 1; k < nsize; k++) 
            {
                matrix[j][k] -= pivotval * matrix[i][k];
            }
              R[j] -= pivotval * R[i];
        }

        //__clkrts_end_cilk();
      }
  
     
}

// Solve the equation.
//
void solveGauss(int nsize) 
{
  int i, j;

  X[nsize-1] = R[nsize-1];
  for (i = nsize - 2; i >= 0; i --) 
  {
    X[i] = R[i];
    for (j = nsize - 1; j > i; j--) 
    {
      X[i] -= matrix[i][j] * X[j];
    }
  }

  #ifdef DEBUG
  fprintf(stdout, "X = [");
  for (i = 0; i < nsize; i++) 
  {
    fprintf(stdout, "%.6f ", X[i]);
   }
  fprintf(stdout, "];\n");
#endif
}

int main(int argc, char *argv[])
 {
   int i,iter;
   struct timeval start, finish;
   int nsize = 0;
   double error;

  if (argc != 2)
   {
    fprintf(stderr, "usage: %s <matrixfile>\n", argv[0]);
    exit(-1);
   }

 // for(iter=4;iter<=104;iter=iter+4)
  {	nsize = initMatrix(argv[1]);
  	initRHS(nsize);
  	initResult(nsize);

 	 gettimeofday(&start, 0);
  		computeGauss(nsize,iter);                //parallelize this
  	 gettimeofday(&finish, 0);

  	solveGauss(nsize);

  fprintf(stdout, "Time:  %f seconds\n", (finish.tv_sec - start.tv_sec)
      + (finish.tv_usec - start.tv_usec) * 0.000001);

  error = 0.0;
  for (i = 0; i < nsize; i++) 
  {
    double error__ = (X__[i]==0.0) ? 1.0 : fabs((X[i]-X__[i]) / X__[i]);
    if (error < error__) 
    {
      error = error__;
    }
  }
  fprintf(stdout, "Error: %e\n", error);
 }
  return 0;
}
