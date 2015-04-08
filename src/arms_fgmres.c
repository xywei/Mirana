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

#define TOL 1.0e-8    /* tolerance for stopping iteration */

#define TOL_DD 0.7   /* diagonal dominance tolerance for */
                     /* independent sets                 */

#define MAXITS 200   /* Maximum number of outer iterations */

#define MAX_NUM_LEV 10  /* maximum number of levels for arms */

#define IM 60    /* dimension of Krylov subspace in (outer) FGMRES */

#define PERM_TYPE 0           /* Indset (0) / PQ (1)    permutation   */
                              /* note that these refer to completely  */
                              /* different methods for reordering A   */
                              /* 0 = standard ARMS independent sets   */
                              /* 1 = arms with ddPQ ordering          */
                              /* 2 = acoarsening-based ordering [new] */

#define BLOCK_SIZE 30 /* smallest size allowed for last schur comp. */

#define DIAG_SCALE 1 /* diagonal scaling  0:no 1:yes */

#define LFIL0 50    /* initial fill-in parameter */

/** Solve a (square) system using ARMS preconditioned fgmres by calling ITSOL2. Adapted from "mainARMS.c" in ITSOL2/TESTS.
 *  \param N integer, number of rows.
 *  \param val real NNZ-by-1 array, coefficients.
 *  \param col_ind integer NNZ-by-1 array, column indices.
 *  \param row_ptr integer (N+1)-by-1 array, with the last element row_ptr[N] = NNZ+1.
 *  \param rhs real N-by-1 array, the right hand side.
 *  \param sol real N-by-1 array, place for the solution.
 *  \param init real N-by-1 array, initial guess.
 *  \param ie 0 if no error occurs.
 *  \note When calling this function from Fortran, one must notice:
 *  - case sensitivity (since C is case sensitive but Fortran is not);
 *  - possible additional trailing underscore ("_") in C function names (the underscore is added by some Fortran compilers);
 *  - possible additional Fortran libraries that need to be linked to build the executable;
 *  - fortran subroutines are always "call by reference" while C can be either "call by reference" or "call by value". 
 */
void arms_fgmres_(int * N, double val[], int col_ind[], int row_ptr[], \
		  double rhs[], double sol[], double init[], int * ie)
{
  csptr mat = NULL;    /* matrix structure */
  arms ArmsSt = NULL;    /* arms preconditioner structure */
  SMatptr MAT = NULL;    /* Matrix structure for matvecs    */
  SPreptr PRE = NULL;    /* general precond structure       */
  double * tmp;
  int ierr, job=1, i, k, l;
  int n = *N;
  int nnz = row_ptr[n] - 1;
  double fillfact, terr;
  double tolind = TOL_DD;

  int ipar[18];
  double droptol[7], dropcoef[7];
  int lfil_arr[7];
  int its;

  FILE *flog = stdout;

  fprintf(flog, "Entering linear solver, initializing (n=%d,nnz=%d)..\n",n,nnz);
  
  /*-------------------- setup data structure for mat (csptr) struct */
  tmp = (double *) malloc(n*sizeof(double));
  mat = (csptr) malloc( sizeof(SparMat) ); 
  mat->n = n;
  mat->nzcount = (int *)malloc( n*sizeof(int));
  mat->ja = (int **) Malloc( n*sizeof(int *), "setupCS" );
  mat->ma = (double **) Malloc( n*sizeof(double *), "setupCS" );
  
  /* Construct C-Style matrix from [val,col_ind,row-ptr] */
  for (k=0; k<n; k++)
    mat->nzcount[k] = row_ptr[k+1] - row_ptr[k];
  for (k=0; k<n; k++)
    {
     l = mat->nzcount[k];
     i = row_ptr[k] - 1;
     if (l > 0)
       {
	 mat->ja[k] = &(col_ind[i]);
	 mat->ma[k] = &(val[i]);
       }
     else
       {
	 mat->ja[k] = NULL;
	 mat->ma[k] = NULL;
       }
    }
  
  /*-------------------- Preconditioning using ARMS */
  fprintf(flog, "Done.\nSetting up ARMS preconditionor..\n");
  for (i=0; i<17; i++)
    ipar[i] = 0;
  
  ipar[0] = MAX_NUM_LEV;
  ipar[1] = PERM_TYPE;
  ipar[2] = BLOCK_SIZE;
  ipar[3] = 0; /* whether or not to print statistics */

  /*-------------------- interlevel methods */  
  ipar[10] = 0;       /* Always do permutations - currently not used  */    
  ipar[11] = 0;       /* ILUT or ILUTP - currently only ILUT is implemented */
  ipar[12] = DIAG_SCALE;  /* diagonal row scaling before PILUT */
  ipar[13] = DIAG_SCALE;  /* diagonal column scaling before PILUT */

  /*-------------------- last level methods */
  ipar[14] = 1;       /* Always do permutations at last level */
  ipar[15] = 1;       /* ILUTP for last level(0 = ILUT at last level) */
  ipar[16] = DIAG_SCALE;  /* diagonal row scaling  0:no 1:yes */
  ipar[17] = DIAG_SCALE;  /* diagonal column scaling  0:no 1:yes */

  /*--------- dropcoef (droptol[k] = tol0*dropcoef[k]) ----- */
  dropcoef[0] = 1.6;     /* dropcoef for L of B block */
  dropcoef[1] = 1.6;     /* dropcoef for U of B block */
  dropcoef[2] = 1.6;     /* dropcoef for L\ inv F */
  dropcoef[3] = 1.6;     /* dropcoef for U\ inv E */
  dropcoef[4] = 0.004;   /* dropcoef for forming schur comple. */
  dropcoef[5] = 0.004;   /* dropcoef for last level L */
  dropcoef[6] = 0.004;   /* dropcoef for last level U */

  for (i=0; i<7; i++)
    {
      lfil_arr[i] = LFIL0 * ((int) nnz/n);
      droptol[i] = TOL * dropcoef[i];
    }

  fprintf(flog, "begin arms\n");
  
  ArmsSt = (arms) malloc(sizeof(armsMat));
  setup_arms(ArmsSt);

  ierr = arms2(mat, ipar, droptol, lfil_arr, tolind, ArmsSt, flog);

  fillfact = (double)nnz_arms(ArmsSt, flog)/(double)(nnz + 1);
  fprintf(flog, "ARMS ends, fill factor (mem used) = %f\n", fillfact);

  /*---------------- get rough idea of cond number - exit if too big */
  if(condestArms(ArmsSt, sol, flog) != 0) {
    fprintf(flog, "Not attempting iterative solution: Cond number too big!\n");
    ie[0] = 3;
    return;
  }

  /*-------------------- initial guess */
  for (i=0; i < n; i++)
    sol[i] = init[i];
    //sol[i] = 0.0;
    
  /*-------------------- set up the structs before calling fgmr */
  MAT = (SMatptr) malloc(sizeof(SMat));
  PRE = (SPreptr) malloc(sizeof(SPre));
  MAT->n = n;
  MAT->CS = mat;
  MAT->matvec = matvecCSR;
  PRE->ARMS = ArmsSt;
  PRE->precon = preconARMS;

  /*-------------------- Solving using FGEMRES */
  its = MAXITS;
  ierr = fgmr(MAT, PRE, rhs, sol, TOL, IM, &its, NULL);

  if(its < MAXITS)
    fprintf(flog, "FGMR OK: converged in %d steps...\n", its);
  else
    {
      fprintf(flog, "FGMR not converged in %d steps...\n", MAXITS);
      ie[0] = 4;
      return;
    }

  /*-------------------- calculate residual norm from generated rhs */
  matvec(mat, sol, tmp);
  terr = 0.0;
  for(i = 0; i < n; i++)
    terr += (rhs[i] - tmp[i]) * (rhs[i] - tmp[i]);
  terr = sqrt(terr);
  fprintf(flog, "Residual norm: %lf.\n", terr);
  
  ie[0] = 0;
  free(tmp);
  free(ArmsSt);
  free(MAT);
  free(PRE);
  return;
}
