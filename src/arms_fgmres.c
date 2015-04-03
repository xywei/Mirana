#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>
// ITSOL2:
#include "globheads.h"
#include "defs.h"
#include "protos.h"
#include "ios.h"

#define TOL_DD 0.7   /* diagonal dominance tolerance for */
                     /* independent sets                 */

/** Solve a (square) system using ARMS preconditioned fgmres by calling ITSOL2.
 *  \param n integer, number of rows.
 *  \param val real NNZ-by-1 array, coefficients.
 *  \param col_ind integer NNZ-by-1 array, column indices.
 *  \param row_ptr integer (N+1)-by-1 array, with the last element row_ptr[N] = NNZ+1.
 *  \param rhs real N-by-1 array, the right hand side.
 *  \param solu real N-by-1 array, place for the solution.
 *  \param ierr 0 if no error occurs.
 *  \note When calling this function from Fortran, one must notice:
 *  - case sensitivity (since C is case sensitive but Fortran is not);
 *  - possible additional trailing underscore ("_") in C function names (the underscore is added by some Fortran compilers);
 *  - possible additional Fortran libraries that need to be linked to build the executable;
 *  - fortran subroutines are always "call by reference" while C can be either "call by reference" or "call by value". 
 */
void arms_fgmres_(int * n, double val[], int col_ind[], int row_ptr[], \
		 double rhs[], double solu[], int * ierr)
{
  int N = *n;
  int i,j;
  int NNZ = row_ptr[N]-1;
  printf("*---------------------------------------------------------\n");
  printf("* Solving linear system using ARMS preconditioned FGMRES,\n");
  printf("* Matrix information: N = %d, NNZ = %d.\n",N,NNZ);
  
  ierr[0] = 0;
}


