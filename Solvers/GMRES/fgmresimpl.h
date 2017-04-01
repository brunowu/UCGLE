#if !defined(__GMRES)
#define __GMRES

#include "petsc/private/kspimpl.h"        /*I "petscksp.h" I*/

#define KSPGMRESHEADER                                                  \
  /* Hessenberg matrix and orthogonalization information. */            \
  PetscScalar *hh_origin;   /* holds hessenburg matrix that has been multiplied by plane rotations (upper tri) */ \
  PetscScalar *hes_origin;  /* holds the original (unmodified) hessenberg matrix which may be used to estimate the Singular Values of the matrix */ \
  PetscScalar *cc_origin;   /* holds cosines for rotation matrices */   \
  PetscScalar *ss_origin;   /* holds sines for rotation matrices */     \
  PetscScalar *rs_origin;   /* holds the right-hand-side of the Hessenberg system */ \
                                                                        \
  PetscScalar *orthogwork; /* holds dot products computed in orthogonalization */ \
                                                                        \
  /* Work space for computing eigenvalues/singular values */            \
  PetscReal   *Dsvd;                                                    \
  PetscScalar *Rsvd;                                                    \
                                                                        \
                                                                        \
  PetscReal haptol;                                       /* tolerance for happy ending */ \
  PetscInt  max_k;                                        /* number of vectors in Krylov space, restart size */ \
  PetscInt  nextra_vecs;                                  /* number of extra vecs needed, e.g. for a pipeline */ \
                                                                        \
  PetscErrorCode (*orthog)(KSP,PetscInt);                    \
  KSPGMRESCGSRefinementType cgstype;                                    \
                                                                        \
  Vec      *vecs;                                        /* the work vectors */ \
  PetscInt q_preallocate;    /* 0=don't preallocate space for work vectors */ \
  PetscInt delta_allocate;    /* number of vectors to preallocaate in each block if not preallocated */ \
  PetscInt vv_allocated;      /* number of allocated gmres direction vectors */ \
  PetscInt vecs_allocated;                              /*   total number of vecs available */ \
  /* Since we may call the user "obtain_work_vectors" several times, we have to keep track of the pointers that it has returned */ \
  Vec      **user_work;                                              \
  PetscInt *mwork_alloc;       /* Number of work vectors allocated as part of  a work-vector chunck */ \
  PetscInt nwork_alloc;        /* Number of work vector chunks allocated */ \
                                                                        \
  /* Information for building solution */                               \
  PetscInt    it;              /* Current iteration: inside restart */  \
  PetscScalar *nrs;            /* temp that holds the coefficients of the Krylov vectors that form the minimum residual solution */ \
  Vec         sol_temp;        /* used to hold temporary solution */

typedef struct {
  KSPGMRESHEADER
} KSP_GMRES;

PETSC_INTERN PetscErrorCode KSPView_GMRES(KSP,PetscViewer);
PETSC_INTERN PetscErrorCode KSPSetUp_GMRES(KSP);
PETSC_INTERN PetscErrorCode KSPSetFromOptions_GMRES(PetscOptions *PetscOptionsObject,KSP);
PETSC_INTERN PetscErrorCode KSPComputeExtremeSingularValues_GMRES(KSP,PetscReal*,PetscReal*);
PETSC_INTERN PetscErrorCode KSPComputeEigenvalues_GMRES(KSP,PetscInt,PetscReal*,PetscReal*,PetscInt*);
PETSC_INTERN PetscErrorCode KSPReset_GMRES(KSP);
PETSC_INTERN PetscErrorCode KSPDestroy_GMRES(KSP);
PETSC_INTERN PetscErrorCode KSPGMRESGetNewVectors(KSP,PetscInt);

typedef PetscErrorCode (*FCN)(KSP,PetscInt); /* force argument to next function to not be extern C*/

PETSC_INTERN PetscErrorCode KSPGMRESSetHapTol_GMRES(KSP,PetscReal);
PETSC_INTERN PetscErrorCode KSPGMRESSetPreAllocateVectors_GMRES(KSP);
PETSC_INTERN PetscErrorCode KSPGMRESSetRestart_GMRES(KSP,PetscInt);
PETSC_INTERN PetscErrorCode KSPGMRESGetRestart_GMRES(KSP,PetscInt*);
PETSC_INTERN PetscErrorCode KSPGMRESSetOrthogonalization_GMRES(KSP,FCN);
PETSC_INTERN PetscErrorCode KSPGMRESGetOrthogonalization_GMRES(KSP,FCN*);
PETSC_INTERN PetscErrorCode KSPGMRESSetCGSRefinementType_GMRES(KSP,KSPGMRESCGSRefinementType);
PETSC_INTERN PetscErrorCode KSPGMRESGetCGSRefinementType_GMRES(KSP,KSPGMRESCGSRefinementType*);

#endif


#if !defined(__FGMRES)
#define __FGMRES

#include <petsc/private/kspimpl.h>
#define KSPGMRES_NO_MACROS

typedef struct {
  KSPGMRESHEADER

  /* new storage for fgmres */
  Vec *prevecs;                  /* holds the preconditioned basis vectors for fgmres.
                                    We will allocate these at the same time as vecs
                                    above (and in the same "chunks". */
  Vec **prevecs_user_work;       /* same purpose as user_work above, but this one is
                                    for our preconditioned vectors */

  /* we need a function for interacting with the pcfamily */

  PetscErrorCode (*modifypc)(KSP,PetscInt,PetscInt,PetscReal,void*);    /* function to modify the preconditioner*/
  PetscErrorCode (*modifydestroy)(void*);
  void *modifyctx;
} KSP_FGMRES;

#define HH(a,b)  (fgmres->hh_origin + (b)*(fgmres->max_k+2)+(a))
/* HH will be size (max_k+2)*(max_k+1)  -  think of HH as
   being stored columnwise for access purposes. */
#define HES(a,b) (fgmres->hes_origin + (b)*(fgmres->max_k+1)+(a))
/* HES will be size (max_k + 1) * (max_k + 1) -
   again, think of HES as being stored columnwise */
#define CC(a)    (fgmres->cc_origin + (a)) /* CC will be length (max_k+1) - cosines */
#define SS(a)    (fgmres->ss_origin + (a)) /* SS will be length (max_k+1) - sines */
#define RS(a)    (fgmres->rs_origin + (a)) /* RS will be length (max_k+2) - rt side */

// /* vector names */
#define VEC_OFFSET     2
#define VEC_TEMP       fgmres->vecs[0]               /* work space */
#define VEC_TEMP_MATOP fgmres->vecs[1]               /* work space */
#define VEC_VV(i)      fgmres->vecs[VEC_OFFSET+i]    /* use to access
                                                        othog basis vectors */
#define PREVEC(i)      fgmres->prevecs[i]            /* use to access
                                                        preconditioned basis */

#endif
