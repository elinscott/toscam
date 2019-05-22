program units
   real(8) :: dd

   dd = 1.d0

   dd = (1.88972)**2*2.d0*dacos(-1.d0)*(0.0002434134)*(1973.26)**2
   dd = dd*(0.036749)**2*(1.d0/1.d-8)*(1973.26)**2/((510998.9098)**2)
   write (*, *) 'DD is : ', dd

end program

