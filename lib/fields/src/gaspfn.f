
       subroutine  gaspfn(n,d2, par)
       real*8 d2(n), par(1)
       integer n

         do 5 k =1,n
         d2(k)=  exp(-1*d2(k))
   5     continue

        return
        end
