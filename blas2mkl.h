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

#ifndef BLAS2MKL_INCLUDE
#define BLAS2MKL_INCLUDE

/* BLAS Level 1 */
#define drotg_  DROTG
#define drotmg_ DROTMG
#define drot_   DROT
#define drotm_  DROTM
#define dswap_  DSWAP
#define dscal_  DSCAL
#define dcopy_  DCOPY
#define daxpy_  DAXPY
#define ddot_   DDOT
#define ddotu_  DDOTU
#define ddotc_  DDOTC
#define dnrm2_  DNRM2
#define dasum_  DASUM
#define idamax_ IDAMAX

/* BLAS Level 2 */
#define dgemv_  DGEMV
#define dgbmv_  DGBMV
#define dhemv_  DHEMV
#define dhbmv_  DHBMV
#define dhpmv_  DHPMV
#define dsymv_  DSYMV
#define dsbmv_  DSBMV
#define dspmv_  DSPMV
#define dtrmv_  DTRMV
#define dtbmv_  DTBMV
#define dtpmv_  DTPMV
#define dtrsv_  DTRMV
#define dtbsv_  DTBMV
#define dtpsv_  DTPMV

#define dger_   DGER
#define dgeru_  DGERU
#define dgerc_  DGERC
#define dher_   DHER
#define dhrp_   DHRP
#define dher2_  DHER2
#define dhrp2_  DHRP2
#define dsyr_   DSYR
#define dspr_   DSPR
#define dsyr2_  DSYR2
#define dspr2_  DSPR2

/* BLAS Level 3 */
#define dgemm_  DGEMM
#define dsymm_  DSYMM
#define dhemm_  DHEMM
#define dsyrk_  DSYRK
#define dherk_  DHERK
#define dsyr2k_ DSYR2K
#define dher2k_ DHER2K
#define dtrmm_  DTRMM
#define dtrsm_  DTRSM

#endif
