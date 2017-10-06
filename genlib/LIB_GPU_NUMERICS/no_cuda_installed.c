#include <stdio.h>
#include <stdlib.h>

 ///////////
 // TOOLS //
 ///////////

extern "C" void  cufftPlan1d_C2C_(int *plan, int *nx,  int *batch){};
extern "C" void  cufftPlan1d_Z2Z_(int *plan, int *nx,  int *batch){};
extern "C" void  cufftPlan1d_Z2D_(int *plan, int *nx,  int *batch){};
extern "C" void  cufftPlan1d_D2Z_(int *plan, int *nx,  int *batch){};
extern "C" void  cufftPlan2d_C2C_(int *plan, int *nx, int *ny, int *batch) {};
extern "C" void  cufftPlan2d_Z2Z_(int *plan, int *nx, int *ny, int *batch)       {};
extern "C" void  cufftPlan2d_Z2D_(int *plan, int *nx, int *ny, int *batch)       {};
extern "C" void  cufftPlan2d_D2Z_(int *plan, int *nx, int *ny, int *batch)       {};
extern "C" void  cufftPlan3d_C2C(int *plan, int *nx, int *ny, int *nz, int *batch) {};
extern "C" void  cufftPlan3d_Z2Z(int *plan, int *nx, int *ny, int *nz, int *batch)        {};
extern "C" void  cufftPlan3d_Z2D(int *plan, int *nx, int *ny, int *nz, int *batch)        {};
extern "C" void  cufftPlan3d_D2Z(int *plan, int *nx, int *ny, int *nz, int *batch)        {};
extern "C" void  cufftDestroy_(int *plan) {};
extern "C" void  cufftExecC2C_(int *plan, float *idata, float *odata, int *direction, int *nx, int *batch){};
extern "C" void  cufftExecZ2Z_(int *plan, double *idata, double *odata, int *direction, int *nx, int *batch){};
extern "C" void  cufftExecZ2D_(int *plan, double *idata, double *odata, int *direction, int *nx, int *batch){};
extern "C" void  cufftExecD2Z_(int *plan, double *idata, double *odata, int *direction, int *nx, int *batch){};

extern "C" void cuinitf_(int *dev) { }
extern "C" void cudevicegetcountf_(int *dev) { }
extern "C" void cusetdevicef_(int *dev) { }
extern "C" void cugetdevicef_(int *dev) { }

extern "C" void compare_cuda_(){} 
extern "C" void cuda_complex_invert_(int *pBlockSize, double** A_, double** invA_, int *pndim){}
extern "C" void cuda_c_invert_(int *pBlockSize, double** A_, double** A_i, double** invA_, double** invA_i, int *pndim){}
extern "C" void matmul_gpu_fortran_(int* pBLOCKSIZE, double** A, double** B, double** C,  int* hA, int* wA, int* wB){}
extern "C" void matmul_gpu_cxx(int* pBLOCKSIZE, double** A, double** B, double** C,  int* hA, int* wA, int* wB){}
extern "C" void matmul_gpu_fortran_c_(int* pBLOCKSIZE, double** A, double** B, double** C,  int* hA, int* wA, int* wB){}
extern "C" void matmul_gpu_cxx_c(int* pBLOCKSIZE, double** A, double** B, double** C,  int* hA, int* wA, int* wB){}
extern "C" void Mul(int size, const double* A, const double* B, double* C){}
extern "C" void cuda_invert__(int *pBlockSize, double** A_, double** invA_, int *pndim){}
extern "C" void cufftplan1d_(int* plan, int* nx, int* cucu, int* batch) { }
extern "C" void cufftexecc2c_(int* plan, float** data, float** data2, int* cucu, int* nx, int* batch) { }
extern "C" void cufftdestroy_(int* plan) { }

extern "C" void  magma_inv_routine_double_all_gpu_single_(int* nn,float* AAA,float* BBB) {}

extern "C" void cublas_matrix_multiply_(double *A, double *B, double *C,int *m, int *n, int *q){};
extern "C" void magma_matrix_multiply_(double *A, double *B, double *C,int *m, int *n, int *q){};
extern "C" void cublas_complex_matrix_multiply_(double *A, double *B, double *C,int *m, int *n, int *q){};
extern "C" void magma_complex_matrix_multiply_single_(double *A, double *B, double *C,int *m, int *n, int *q){};
extern "C" void magma_complex_matrix_multiply_(double *A, double *B, double *C,int *m, int *n, int *q){};
extern "C" void invert_my_matrix_(int *nb, int *n,double *A){};


extern void invertspd__(double * A, double *invAA,int lda, int n){};
extern void invertge__(double * A, double *invAA, int n){};
extern "C" void invert_ge_(int *nn, double *AA, double *invAA){};
extern "C" void invert_spd_(int *nn, double *AA, double *invAA){};
extern "C" void invert_magma_ge_(int *nn, double *AA){};
extern "C" void invert_magma_double_ge_(int *nn, double *AA){};
extern "C" void  magma_inv_routine_double_all_gpu_(int* nn,double* AAA,double* BBB) {};
extern "C" void freemem_(){};
extern "C" void  magma_inv_routine_comp_all_gpu_(int* nn,double* AAA,double* BBB){};
extern "C" void  magma_inv_routine_comp_all_gpu_single_(int* nn,double* AAA,double* BBB){};


 //////////////////  
 // LANCZOS REAL //
 //////////////////

 
extern "C" void hmult_sz_real_cuda_rout_(int *pblocksize, int *offdiasize, int *roff,
int* ntot, double* QUART, double* diagsz, double* vec_in, double* vec_out, int* noffsz, int* rankoffsz, double* offdiagsz, int *rank ){};

extern "C" void lanczos_test_gpu_(int *blocksize, int *ntot){};

extern "C" void lanczos_real_dynamic_cuda_(int *pblocksize, int *Niter_lanczos_,int *offdiasize, int *roff, int *ntot,
double *QUART, double *diagsz, int *noffsz, int *rankoffsz, double *offdiagsz, double *diag, double *subdiag , double *vecinit, int *rank){};

extern "C" void lanczos_real_cuda_(int *pblocksize,  int *Niter_lanczos_,int *offdiasize, int *roff, int *ntot,
double *QUART, double *diagsz, int *noffsz, int *rankoffsz, double *offdiagsz, double *diag, double *subdiag, int *rank ){};

extern "C" void lanczos_real_get_gs_cuda_(int *pblocksize, int *Niter_lanczos_,int *offdiasize, int *roff,
int *ntot, double *QUART, double *diagsz, int *noffsz, int *rankoffsz, double *offdiagsz, double *vecp, double *GS, int *rank){}

extern "C" void lanczos_real_fly_cuda_(int *dimen_in, int *pblocksize, int* norbs_in, int *Niter_lanczos_, int* ntot, double* quart,
double* diag, double* subdiag, double* Eb_in, double* Ec_in, double* Vbc_in, int* sector_states_in, int* sector_ranks_in
,int* bathnorbs_in,int* impnorbs_in,int* imporbs_in, int* bathorbs_in, int* maskEb_in, int* maskEc_in, int* maskVbc_in, int *rank){}

extern "C" void lanczos_real_fly_dynamic_cuda_(int *dimen_in, int *pblocksize, int* norbs_in, int *Niter_lanczos_, int* ntot, double* quart,
double* diag, double* subdiag, double* Eb_in, double* Ec_in, double* Vbc_in, int* sector_states_in, int* sector_ranks_in
,int* bathnorbs_in,int* impnorbs_in,int* imporbs_in, int* bathorbs_in, int* maskEb_in, int* maskEc_in, int* maskVbc_in, double* vecinit, int *rank){}

extern "C" void lanczos_real_fly_gs_cuda_(int *dimen_in, int *pblocksize, int* norbs_in, int *Niter_lanczos_, int* ntot, double* quart,
double* Eb_in, double* Ec_in, double* Vbc_in, int* sector_states_in, int* sector_ranks_in ,int* bathnorbs_in,int* impnorbs_in,
int* imporbs_in, int* bathorbs_in, int* maskEb_in, int* maskEc_in, int* maskVbc_in,double *vecp,double *GS, int *rank){}

extern "C" void lanczos_real_updo_cuda_(int *norbs, int *pblocksize, int *Niter_lanczos_,int *offdiasizeup, int* offdiasizedn,
int *roffup, int *roffdn, int *ntot,  int *sizup, int *sizdn,
double *QUART, double *diagup, double *diagdn, int *noffup, int *noffdn, int *rankoffup, int *rankoffdn,
double *offdiagup, double *offdiagdn, double *diag, double *subdiag, int *UMASK, int *statesup, int *statesdn , int *iorbup, int *iorbdn, int *rank){};

extern "C" void lanczos_real_updo_dynamic_cuda_(int *norbs, int *pblocksize, int *Niter_lanczos_,int *offdiasizeup, int* offdiasizedn,
int *roffup, int *roffdn, int *ntot,  int *sizup, int *sizdn,
double *QUART, double *diagup, double *diagdn, int *noffup, int *noffdn, int *rankoffup, int *rankoffdn,
double *offdiagup, double *offdiagdn, double *diag, double *subdiag, int *UMASK, int *statesup, int *statesdn ,
int *iorbup, int *iorbdn, double *vecinit, int *rank){};

extern "C" void lanczos_real_updo_gs_cuda_(int *norbs, int *pblocksize, int *Niter_lanczos_,int *offdiasizeup, int* offdiasizedn,
   int *roffup, int *roffdn, int *ntot,  int *sizup, int *sizdn,
   double *QUART, double *diagup, double *diagdn, int *noffup, int *noffdn, int *rankoffup, int *rankoffdn,
   double *offdiagup, double *offdiagdn, double *diag, double *subdiag, int *UMASK, int *statesup, int *statesdn ,
   int *iorbup, int *iorbdn, double *vecp, double *GS, int *rank){};




 /////////////////////
 // LANCZOS COMPLEX //
 /////////////////////


extern "C" void hmult_sz_complex_cuda_rout_(int *pblocksize, int *offdiasize, int *roff, int* ntot,
double* QUART, double* diagsz, double* vec_in, double* vec_out, int* noffsz, int* rankoffsz, double* offdiagsz ,int *rank){};

extern "C" void lanczos_dynamic_cuda_complex_(int *pblocksize, int *Niter_lanczos_,int *offdiasize, int *roff,
int *ntot, double *QUART, double *diagsz, int *noffsz, int *rankoffsz, double *offdiagsz, double *diag,
double *subdiag , double *vecinit, int *rank){};

extern "C" void lanczos_cuda_complex_(int *pblocksize,  int *Niter_lanczos_,int *offdiasize, int *roff, int *ntot,
double *QUART, double *diagsz, int *noffsz, int *rankoffsz, double *offdiagsz, double *diag, double *subdiag, int *rank ){};

extern "C" void lanczos_get_gs_cuda_complex_(int *pblocksize, int *Niter_lanczos_,int *offdiasize, int *roff,
int *ntot, double *QUART, double *diagsz, int *noffsz, int *rankoffsz, double *offdiagsz, double *vecp, double *GS, int *rank){};

extern "C" void lanczos_complex_fly_gs_cuda_(int *dimen_in, int *pblocksize, int* norbs_in, int *Niter_lanczos_, int* ntot, double* quart,
  double* Eb_in, double* Ec_in, double* Vbc_in, int* sector_states_in, int* sector_ranks_in
 ,int* bathnorbs_in,int* impnorbs_in,int* imporbs_in, int* bathorbs_in, int* maskEb_in, int* maskEc_in, int* maskVbc_in,double *vecp,double *GS, int *rank){};

extern "C" void lanczos_complex_fly_dynamic_cuda_(int *dimen_in, int *pblocksize, int* norbs_in, int *Niter_lanczos_, int* ntot, 
double* quart, double* diag, double* subdiag, double* Eb_in, double* Ec_in, double* Vbc_in, 
int* sector_states_in, int* sector_ranks_in ,int* bathnorbs_in,int* impnorbs_in,int* imporbs_in, int* bathorbs_in, int* maskEb_in, 
int* maskEc_in, int* maskVbc_in, double* vecinit, int *rank){};

extern "C" void lanczos_complex_fly_cuda_(int *dimen_in, int *pblocksize, int* norbs_in, int *Niter_lanczos_, int* ntot, double* quart,
double* diag, double* subdiag, double* Eb_in, double* Ec_in, double* Vbc_in, int* sector_states_in, 
int* sector_ranks_in ,int* bathnorbs_in,int* impnorbs_in,int* imporbs_in, int* bathorbs_in, int* maskEb_in, int* maskEc_in, 
int* maskVbc_in, int *rank) {};

extern "C" void lanczos_complex_updo_gs_cuda_(int *norbs, int *pblocksize,
   int *Niter_lanczos_,int *offdiasizeup, int* offdiasizedn, int *roffup, int *roffdn, int *ntot,  int *sizup, int *sizdn,
   double *QUART, double *diagup, double *diagdn, int *noffup, int *noffdn, int *rankoffup, int *rankoffdn,
   double *offdiagup, double *offdiagdn,  int *UMASK, int *statesup, int *statesdn,
   int *iorbup, int *iorbdn, double *vecp, double *GS, int *rank){}

extern "C" void lanczos_complex_updo_dynamic_cuda_(int *norbs, int *pblocksize,
   int *Niter_lanczos_,int *offdiasizeup, int* offdiasizedn, int *roffup, int *roffdn, int *ntot,  int *sizup, int *sizdn,
   double *QUART, double *diagup, double *diagdn, int *noffup, int *noffdn, int *rankoffup, int *rankoffdn,
   double *offdiagup, double *offdiagdn, double *diag, double *subdiag, int *UMASK, int *statesup, int *statesdn ,
   int *iorbup, int *iorbdn, double *vecinit, int *rank){}

extern "C" void lanczos_complex_updo_cuda_(int *norbs, int *pblocksize, int *Niter_lanczos_,
   int *offdiasizeup, int* offdiasizedn, int *roffup, int *roffdn, int *ntot,  int *sizup, int *sizdn,
   double *QUART, double *diagup, double *diagdn, int *noffup, int *noffdn, int *rankoffup, int *rankoffdn,
   double *offdiagup, double *offdiagdn, double *diag, double *subdiag, int *UMASK, int *statesup, int *statesdn , int *iorbup, int *iorbdn, int *rank){}


 /////////////////////
 // SUM OF INVERSE  //
 /////////////////////

extern "C" void sum_of_inverse_frequ_(int* nnn_, int* nfrequ_,  double* Eb_, double* totsum_ , double *frequ_){}
extern "C" void sum_of_inverse_frequ_collect_(int* nnn_, int* nfrequ_,  double* Eb_, double* totsum_ , double *frequ_){}
extern "C" void sum_of_inverse_frequ_array_(int* nnn_, int* nfrequ_, double* totsum_){}
extern "C" void sum_of_inverse_frequ_complex_(int* nnn_, int* nfrequ_,  double* Eb_, double* totsum_ , double *frequ_, int* firstlast){}
extern "C" void sum_of_inverse_frequ_complex_collect_(int* nnn_, int* nfrequ_,  double* Eb_, double* collect_ , double *frequ_, int* firstlast){}
extern "C" void sum_of_inverse_frequ_complex_array_(int* nnn_, int* nfrequ_, double* collect_ , int* firstlast){}


// MAGMA 

 extern "C" void magma_dsygvd_( int *NNN, double *h_A, double *h_S, double *w){}
 extern "C" void magma_ssygvd_( int *NNN, float *h_A, float *h_S, float *w){}
 extern "C" void magma_dsyevd_( int *NNN, double *h_A, double *w){}
 extern "C" void magma_ssyevd_( int *NNN, float *h_A, float *w){}
 extern "C" void magma_zhegvd_type_(int *ibtype,  int *NNN, double *h_A, double *h_S, double *w, double *vec){}
 extern "C" void magma_chegvd_type_(int *ibtype,  int *NNN, double *h_A, double *h_S, float *w, double *vec){}
 extern "C" void magma_dsygvd_type_( int *ibtype, int *NNN, double *h_A, double *h_S, double *w, double *vec){}
 extern "C" void magma_ssygvd_type_( int *ibtype, int *NNN, float *h_A, float *h_S, float *w, float *vec){}



