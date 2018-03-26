   subroutine lanc_simp()

      use globalvar_ed_solver, only: cutoff_min_lanczos_vec, USE_TRANSPOSE_TRICK_MPI

      implicit none
 
      !---------------------------!
      ! 1 LANCZOS RECURSION STEP  !
      !---------------------------!

      if(verbose) write(*,*) 'iter normvin normvout : ', iter,normv_in,normv_out

      !----------------------------------------------------------------------!
      if(iter==1)then
         if(abs(normv_in)<cutoff_min_lanczos_vec) then 
            write(*,*) 'lanczos vector is : ', vec_in
            write(*,*) 'norm of vector is : ', normv_in
            write(*,*) 'vec in            : ', vec_in
            write(*,*) 'iter              : ', iter
            write(*,*) 'use mpi trick     : ', USE_TRANSPOSE_TRICK_MPI
            stop 'Lanczos ED solver starting vector with 0 norm'
         endif
         aa=1.d0/normv_in
         vec_in  = vec_in*aa
      else
         if(abs(normv_out)<cutoff_min_lanczos_vec)then
            write(*,*) 'Danger, 0 norm in Lanczos, iter=',iter
            write(*,*) 'normv_in, normv_out : ', normv_in, normv_out
            write(*,*) 'norm v in calculated = ', SQRT(ABS(MPI_DOT_PRODUCT(vec_in,vec_in,split=USE_TRANSPOSE_TRICK_MPI)))
            stop 'critical'
         else  
            aa=1.d0/normv_out
            vec_in  = vec_in*aa
         endif
      endif
      !----------------------------------------------------------------------!

      CALL Hmult__(size(vec_in),vec_out,vec_in)

      if(verbose) write(*,*) 'max vecout : ', maxval(abs(vec_out))

      if(iter>1) vec_out         =  vec_out - tri%subdiag(iter)*tmp
      tmp             =  vec_in
      tri%diag(iter)  =  DBLE(MPI_DOT_PRODUCT(vec_out,vec_in,split=USE_TRANSPOSE_TRICK_MPI))
      vec_out         =  vec_out - tri%diag(iter) * tmp
      normv_out       =  SQRT(ABS(MPI_DOT_PRODUCT(vec_out,vec_out,split=USE_TRANSPOSE_TRICK_MPI)))
      IF(iter<tri%N) tri%subdiag(iter+1)  = normv_out
      vec_in          =  vec_out

   end subroutine

   subroutine lanc_opt()

      use globalvar_ed_solver, only: cutoff_min_lanczos_vec, USE_TRANSPOSE_TRICK_MPI

      implicit none

!$OMP PARALLEL PRIVATE(i) 
      if(iter==1)then
         if(abs(normv_in)<cutoff_min_lanczos_vec)then
            write(*,*) 'cutoff min Lanczos : ', cutoff_min_lanczos_vec
            write(*,*) 'normv_in           : ', normv_in
            write(*,*) 'vec in             : ', vec_in
            write(*,*) 'iter               : ', iter
            write(*,*) 'use mpi trick      : ', USE_TRANSPOSE_TRICK_MPI
            stop 'starting lanczos with null vector, critical'
         endif
         aa=1.d0/normv_in
!$OMP DO
         do i=1,size(vec_in)
            vec_in(i)  = vec_in(i)*aa
         enddo
!$OMP END DO
      else
         if(abs(normv_out)<cutoff_min_lanczos_vec)then
            write(*,*) 'vec in  : ', vec_in
            write(*,*) 'vec out : ', vec_out
            write(*,*) 'norm    : ', normv_out
            write(*,*) 'iter    : ', iter
            write(*,*) 'dimen   : ', size(vec_in)
            stop '0 norm vector in lanczos, critical'
         endif
         aa=1.d0/normv_out
!$OMP DO
         do i=1,size(vec_in)
            vec_in(i)  = vec_in(i)*aa
         enddo
!$OMP END DO
      endif
!$OMP END PARALLEL

      CALL Hmult__(size(vec_in),vec_out,vec_in)

      ddiag=0.d0
#ifndef OPENMP_MPI_SAFE
!$OMP PARALLEL PRIVATE(i) REDUCTION(+:ddiag)
#endif
      if(iter>1) then
#ifndef OPENMP_MPI_SAFE
!$OMP DO 
#endif
         do i=1,size(vec_out)
            vec_out(i)     =  vec_out(i) - tri%subdiag(iter)*tmp(i)
#ifdef _complex
            ddiag =  ddiag + dble(conjg(vec_out(i))*vec_in(i))
#else
            ddiag =  ddiag + vec_out(i)*vec_in(i)
#endif
         enddo
#ifndef OPENMP_MPI_SAFE
!$OMP END DO
#endif
      else
#ifndef OPENMP_MPI_SAFE
!$OMP DO
#endif
         do i=1,size(vec_out)
#ifdef _complex
            ddiag =  ddiag + dble(conjg(vec_out(i))*vec_in(i))
#else
            ddiag =  ddiag + vec_out(i)*vec_in(i)
#endif
         enddo
#ifndef OPENMP_MPI_SAFE
!$OMP END DO
#endif
      endif
#ifndef OPENMP_MPI_SAFE
!$OMP END PARALLEL 
#endif
      tri%diag(iter)=ddiag

      normv_out = 0.d0
#ifndef OPENMP_MPI_SAFE
!$OMP PARALLEL PRIVATE(i) REDUCTION(+:normv_out)
!$OMP DO
#endif
      do i=1,size(vec_out) 
         vec_out(i)      =  vec_out(i) - tri%diag(iter) * vec_in(i)
#ifdef _complex
         normv_out       =  normv_out + dble(conjg(vec_out(i))*vec_out(i))
#else
         normv_out       =  normv_out + vec_out(i)*vec_out(i)
#endif
      enddo
#ifndef OPENMP_MPI_SAFE
!$OMP END DO
!$OMP END PARALLEL
#endif
      normv_out=sqrt(normv_out)

      IF(iter<tri%N)     tri%subdiag(iter+1)  = normv_out

!$OMP PARALLEL PRIVATE(i)
!$OMP DO
      do i=1,size(vec_in)
         tmp(i)          =  vec_in(i)
         vec_in(i)       =  vec_out(i)
      enddo
!$OMP END DO
!$OMP END PARALLEL

   end subroutine

