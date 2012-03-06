/******************************************************************************/
/* For all support and information please contact:                            */
/*                                                                            */
/*   Sigurd Enghoff                                                           */
/*   The Salk Institute, CNL                                                  */
/*   enghoff@salk.edu                                                         */
/*                                                                            */
/* Additional ICA software:                                                   */
/*   http://www.cnl.salk.edu/~enghoff/                                        */
/*                                                                            */
/******************************************************************************/

#ifndef __mpiica_h
#define __mpiica_h

#include "mpi.h"
#include "ica.h"

/* Command to be excuted to spawn new ICA processes */
#define SPAWN_ICA "ica"

/* Nice value used by slave processes */
#define SLAVE_NICE 20

/* Number of integer elements in the assign buffer */
#define ASSIGN_NINT 13

/* Integer elements in the assign buffer */
#define ASSIGN_ID         0
#define ASSIGN_CHANS      1
#define ASSIGN_FRAMES     2
#define ASSIGN_EPOCHS     3
#define ASSIGN_BIAS       4
#define ASSIGN_SIGNS      5
#define ASSIGN_EXTENDED   6
#define ASSIGN_EXTBLOCKS  7
#define ASSIGN_PDFSIZE    8 
#define ASSIGN_NSUB       9
#define ASSIGN_VERBOSE    10 
#define ASSIGN_BLOCK      11
#define ASSIGN_MAXSTEPS   12


/* Number of double elements in the MPI assign buffer */
#define ASSIGN_NDBL       5

/* Double elements in the MPI assign buffer */
#define ASSIGN_LRATE      0
#define ASSIGN_ANNEALSTEP 1
#define ASSIGN_ANNEALDEG  2
#define ASSIGN_NOCHANGE   3
#define ASSIGN_MOMENTUM   4

/* Structure processed by the assign send and receive routines */
typedef struct {
	integer     id;

	doublereal *data;
	doublereal *weights;
	integer     chans;
	integer     frames;
	integer     epochs;
	integer     bias;      /* boolean value */
	integer     signs;     /* boolean value */

	integer     extended;
	integer     extblocks;
	integer     pdfsize;
	integer     nsub;
	integer     verbose;
	integer     block;
	integer     maxsteps;
	doublereal  lrate;
	doublereal  annealstep;
	doublereal  annealdeg;
	doublereal  nochange;
	doublereal  momentum;
} tassign;


/* Number of integer elements in the MPI result buffer */
#define RESULT_NINT  4

/* Integer elements in the MPI result buffer */
#define RESULT_ID    0
#define RESULT_CHANS 1
#define RESULT_BIAS  2
#define RESULT_SIGNS 3


/* Number of double elements in the result buffer */
#define RESULT_NDBL  1

/* Double elements in the result buffer */
#define RESULT_LRATE 0

/* Structure processed by the result send and receive routines */
typedef struct {
	integer     id;

	doublereal *weights;
	integer     chans;
	doublereal *bias;
	integer    *signs;

	doublereal  lrate;
} tresult;

extern int slave();
extern void mpiica(doublereal*, doublereal*, doublereal*, doublereal*, integer, integer, integer, integer, int*, doublereal*, integer*, char**);

#endif
