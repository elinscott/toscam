 module lockmod

  public          :: maxt

  INTEGER         :: NTHREADS, TID, MAXT, PROCS, CHUNK
  LOGICAL         :: INPAR, DYNAMIC, NESTED

  type omp_lock_t
    sequence   !!! Based on the standards POINTER(IL,LOCK)
    integer, pointer :: omp_lock_t_ptr
  end type

 end module
