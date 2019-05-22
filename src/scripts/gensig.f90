program gensigg
   integer :: i, j, k, icol, ifrequ
   real(8), allocatable :: gensig(:)
   real(8) :: a1, a2, beta
   integer :: imatsu

   write (*, *) 'real 0 or matsu 1 ?'
   read (*, *) imatsu
   write (*, *) 'how many columns'
   read (*, *) icol
   write (*, *) 'how many frequ'
   read (*, *) ifrequ
   allocate (gensig(ifrequ))
   write (*, *) 'start,end'
   read (*, *) a1, a2

   if (imatsu == 1) then
      write (*, *) 'beta'
      read (*, *) beta
      gensig = (/(dacos(-1.d0)/beta*(2*i - 1), i=1, ifrequ)/)
   else
      gensig = (/(a1 + (a2 - a1)*(i - 1)/(ifrequ - 1), i=1, ifrequ)/)
   endif

   do i = 1, ifrequ
      write (10, '(100f10.4)') gensig(i), (0.D0, j=1, icol)
   enddo

end program
