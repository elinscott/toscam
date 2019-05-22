program groupdimers

   use StringManip, only: StrInt2, toString

   implicit none

   integer             :: i, j, k, l
   character(300)      :: tmp, list(1000)
   character(3000)     :: ttt
   integer, allocatable :: dimer_table(:, :)
   logical, allocatable :: dimer_crossed(:)
   logical             :: check

   INQUIRE (file='mask_dimer', EXIST=check)
   if (.not. check) then
      write (*, *) 'doing single site DMFT'
      stop
   endif

   i = 0
   open (unit=313, file='mask_dimer')
   do
      read (313, *, end=54)
      i = i + 1
   enddo
54 continue
   allocate (dimer_table(i, 2))
   dimer_table = 0
   rewind (313)
   do j = 1, i
      write (*, *) 'j, max      : ', j, i
      write (*, *) 'shape dimer : ', shape(dimer_table)
      read (313, *) dimer_table(j, 1), dimer_table(j, 2)
      write (*, *) 'DIMER TABLE = ', dimer_table(j, 1:2)
   enddo
   close (313)

   allocate (dimer_crossed(i))

   write (*, *) 'green spin up'
  call system(" echo > collect_rank_files ; for i in  `ls green_output*1 2> /dev/null` ; do echo $i >> collect_rank_files; done ; ")
   call goforit(1)

   write (*, *) 'green spin down'
  call system(" echo > collect_rank_files ; for i in  `ls green_output*2 2> /dev/null` ; do echo $i >> collect_rank_files; done ; ")
   call goforit(1)

   write (*, *) 'sigma spin up'
  call system(" echo > collect_rank_files ; for i in  `ls sigma_output*1 2> /dev/null` ; do echo $i >> collect_rank_files; done ; ")
   call goforit(2)

   write (*, *) 'sigma spin down'
  call system(" echo > collect_rank_files ; for i in  `ls sigma_output*2 2> /dev/null` ; do echo $i >> collect_rank_files; done ; ")
   call goforit(2)

contains

!----------------------------------------------------------!
   subroutine goforit(kkk)
      implicit none
      integer :: iii, kkk, kk, maxj

!kkk=1 green
!kkk=2 sigma

      write (*, *) 'opening file collect_rank_files'
      i = 0
      maxj = -1000
      open (unit=1010, file='collect_rank_files')
      read (1010, *, end=66)
      do
         read (1010, *, end=66) tmp
         i = i + 1
         call get_atom_number(j, tmp)
         list(j) = tmp
         if (j > maxj) maxj = j
         write (*, *) 'file is :', i, list(j)
      enddo
66    continue
      close (1010)

      if (i == 0) then
         write (*, *) 'there are simply no files'
         return
      endif

      dimer_crossed = .false.

      write (*, *) 'THERE ARE [X] SITES (not all correlated) : ', maxj

      do k = 1, maxj

         do kk = 1, size(dimer_table, 1)
            if (k == dimer_table(kk, 1)) then
               goto 52
            endif
         enddo
         goto 53
52       continue

         if (.not. dimer_crossed(kk)) then

            dimer_crossed(kk) = .true.

            if (dimer_table(kk, 2) > 0) then
               do iii = 1, size(dimer_table, 1)
                  if (dimer_table(iii, 1) == dimer_table(kk, 2)) then
                     dimer_crossed(iii) = .true.
                     goto 94
                  endif
               enddo
               write (*, *) 'ERROR in onetep_group_dimer, could not find conjugate site of the dimer'
               stop
94             continue
            endif

            if (dimer_table(kk, 2) > 0) then
               call system(" echo >> i_am_dimer"//trim(adjustl(tostring(dimer_table(kk, 1)))))
   ttt=TRIM(ADJUSTL(list(dimer_table(kk,1))))//" "//TRIM(ADJUSTL(list( dimer_table(kk,2) )))//" "//TRIM(ADJUSTL(list(dimer_table(kk,1))))//"_dimer "//TRIM(ADJUSTL(list( dimer_table(kk,2) )))//"_dimer "
            else
               ttt = TRIM(ADJUSTL(list(dimer_table(kk, 1))))
            endif

            write (*, *) 'files are : ', TRIM(ADJUSTL(ttt))

            if (kkk == 1) then
               call system(" dmft_project.out  "//TRIM(ADJUSTL(ttt)))
            else
               call system(" dmft_project_sigma.out "//TRIM(ADJUSTL(ttt)))
            endif
         endif

53       continue
      enddo

      call system(" rm collect_rank_files 2> /dev/null ")
   end subroutine
!----------------------------------------------------------!

   subroutine get_atom_number(newatom, filename)
      implicit none
      integer       :: newatom, j, i
      character*(*) :: filename
      do i = LEN(filename), 1, -1
         if (filename(i:i) == "t") then
            exit
         endif
      enddo
      i = i + 1
      do j = i, 2000
         if (filename(j:j) == "_") then
            exit
         endif
      enddo
      j = j - 1
      write (*, *) 'ATOM NUMBER : ', filename(i:j)
      newatom = StrInt2(TRIM(ADJUSTL(filename(i:j))))
   end subroutine
!----------------------------------------------------------!
end program
