
       subroutine drdfun(n,d2, par)
       double precision d2(n), par(2), dtemp
       integer n
       if( int(par(2)).eq.0) then

         do 5 k =1,n
         d2(k)= par(1)*(d2(k))**( par(1)-1)
   5     continue
        else 
         do 6 k=1,n
          dtemp= d2(k)
          if( dtemp.GE.1e-35)  then
           d2(k)=  (par(1)*log(dtemp) +1)*(dtemp)**( par(1)-1)
          else
           d2(k)=0.0
          endif
   6   continue
       endif
        return
        end
