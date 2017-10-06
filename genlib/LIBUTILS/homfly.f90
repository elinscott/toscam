module HOMFLYPOL

 use geometry
 use linalg
 use sorting
 use random
 use StringManip

real(8),dimension(:,:),allocatable,private :: rback

contains
 
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************

function sizeanalysis2(rr,STRING,step,keep_segment,change,icoord1,icoord2,howmuch)
implicit none
integer :: sizeanalysis2,change,howmuch,icoord1,icoord2,segmentsinv
real(8) :: AA1(3),AA2(3),BB1(3),BB2(3)
integer::n,m_compt,i,j_compt,k,l,ii,ij,ik,in,ierr,long,longueur,kkk,kkkk
integer::step,start,isize,step_size,keep_segment,compt_size,sample
integer::counter,stp,failure,sn,snn,maxcross,nn,idm,irm,kdm,jdm,lp,lrev
integer::j,m,kst,jst,ia,ib,junk,lst,lmax,iminus,npts,ind,ist,kk
integer,parameter::sts=4
integer::tag(1001),r_tag(1001),rotag(1001),inlet(10001),onlet(10001),num(10001)
integer::ith(10001,4)
real(8)::tvalues(2001),svalues(2001)
real(8)::tintix(10001),t(999),inter(10001),ss(6),sols(2),coords(3)
real(8)::points(999,3),setpts(999,3),pts(999,3),d_pts(999,3),p_end(10001,6)
real(8)::reverse(10001,2),p_app(10001,6),orso(2001,2),sutions(20001,2)
real(8)::th(10001,6),r_pts(999,3),r_delta(999,3),rotp(999,3),rotd(999,3)
real(8)::p(999,6),solutions(2001,2),coordinates(2001,3)
real(8)::ordsols(2001,2,2001)
real(8)::writhe
real(8),intent(in) ::rr(:,:)
real(8),dimension(:,:),allocatable::r_size,r_temp,mat,vect,norm,rot,matrice1,matrice2,matrice3
real(8)::teta,phi,func,acn_value,check_norm
real(8)::diff,x1,y1,z1,zj,t1,tj,sig1,sgg1,qq1,qq2,pp1,pp2,subangle,totalcross,xj,yj
integer :: NNN
logical :: okey
character(1)::ilet(0:10001),olet(0:10001),cdwr(0:10001,0:9,4),lett(10001,24)    
character(3)::ch1,ch2,ch4,ch5,ch6 
character(4)::cd(0:10001,0:9),minus
logical::vrai,boolean,logic
character*(*),intent(inout) :: STRING

NNN=size(rr(:,1))
if(change>0.and.icoord1>0.and.icoord2>0)then
 if(icoord1-1>=1) AA1=rr(icoord1-1,:)
 if(icoord1-1<1)  AA1=rr(NNN,:)
                  AA2=rr(icoord1,:)

 if(icoord2-1>=1) BB1=rr(icoord2-1,:)
 if(icoord2-1<1)  BB1=rr(NNN,:)
                  BB2=rr(icoord2,:)
endif
howmuch=0
segmentsinv=0

do i=1,10000
write(ilet(i),'(a1)') ' '
write(olet(i),'(a1)') ' '
enddo
write(ch1,'(a3)') '   '
write(ch2,'(a3)') '   '
write(ch4,'(a3)') '   '
write(ch5,'(a3)') '   '
write(ch6,'(a3)') '   '
write(minus,'(a4)') '    '

!**********************************************************************
!**********************************************************************
!**********************************************************************
sizeanalysis2=200
maxcross=999
call initString(STRING)
!**********************************************************************
!**********************************************************************
!**********************************************************************
!STDIN STDOUT:
if(allocated(rback))deallocate(rback)
allocate(rback(10000,3))

!**********************************************************************         
!**********************************************************************
!**********************************************************************
 
          if(step>=keep_segment) then
               step_size=step
               if(allocated(r_size))deallocate(r_size)
               allocate(r_size(step_size,3))
               r_size=rr
               k=1
               compt_size=0
               logic=.true.
               do while (logic)
                    if (k.eq.1) then
                         call delta2(r_size,step_size,isize)
                         k=k+1
                         if(allocated(r_size))deallocate(r_size)
                    endif
     
                    if (isize.eq.step_size) then
                         compt_size=compt_size+1
                         if (compt_size.le.isize) then
                              step_size=isize
                              if(allocated(r_size)) deallocate(r_size)
                              if(allocated(r_temp)) deallocate(r_temp)
                              allocate(r_size(step_size,3))
                              allocate(r_temp(step_size,3))
                              do i=1,step_size
                                r_size(i,:)=rback(i,:) 
                              enddo
                              r_temp=cshift(r_size,-1,1)
                              call delta2(r_temp,step_size,isize)
                              if(allocated(r_size))deallocate(r_size)
                              if(Allocated(r_temp))deallocate(r_temp)
                         else
                              logic=.false.
                         endif
                    else
                         compt_size=0
                         step_size=isize
                         if(allocated(r_size))deallocate(r_size)
                         allocate(r_size(step_size,3)) 
                         do i=1,step_size
                              r_size(i,:)=rback(i,:)
                         enddo
                         call delta2(r_size,step_size,isize)
                         if(allocated(r_size))deallocate(r_size)
                    endif
               enddo
          
               step_size=isize
               if(allocated(r_size))deallocate(r_size)
               allocate(r_size(step_size,3))
               do i=1,step_size
                    r_size(i,:)=rback(i,:)
               enddo
          else
               step_size=step
               if(allocated(r_size))deallocate(r_size)
               allocate(r_size(step_size,3))
               do i=1,step_size
                    r_size(i,:)=rr(i,:)
               enddo
          endif
          
          counter=0
          lp=1
          lrev=1
          
!**********************************************************************
!**********************************************************************
!**********************************************************************
!RESET VARIABLE TO 0:

          !INIT ROUTINE ET SWAP (inutile.....)
          sutions=0.d0
          reverse=0.d0
          p_app=0.d0
          p_end=0.d0
          inter=0.d0
          num=0
          tintix=0.d0
          inlet=0
          onlet=0
          ilet=' '
          olet=' '
          ith=0
          th=0.d0
          cd=' '
          cdwr(:,:,4)=' '
          lett=' '
          r_tag=0
          tag=0
          rotag=0
          setpts=0.d0
          pts=0.d0
          d_pts=0.d0
          r_pts=0.d0
          r_delta=0.d0
          rotp=0.d0
          rotd=0.d0
!536       format(a80)

          do i=1,step_size
           pts(i,1)=r_size(i,1)
           pts(i,2)=r_size(i,2)
           pts(i,3)=r_size(i,3)
          enddo

          nn=step_size
          counter=counter+1

          do j=1,(nn-1)
           d_pts(j,:)=pts(j+1,:)-pts(j,:)
          enddo
          d_pts(nn,:)=pts(1,:)-pts(nn,:)

          tag= (/ (j, j=1,nn) /)
          
          do j=1,(nn-1)
               r_pts(j,:)=pts(j+1,:)
               r_delta(j,:)=d_pts(j+1,:)
          enddo
          r_pts(nn,:)=pts(1,:)
          r_delta(nn,:)=d_pts(1,:)

          do i=1,(nn-1)
               r_tag(i)=tag(i+1)
          enddo
          r_tag(nn)=tag(1)

          ind=0
          stp=nn
          ordsols=0.d0
          snn =0
          sn=0

704      continue !!!! REVIENT DEPUIS LE DEBUT APRES SWAP

          if ( ind.ge.stp) goto 703
          ind=ind+1
          do j=1,nn
            rotp(j,:)=r_pts(j,:)
            rotd(j,:)=r_delta(j,:)
          enddo
          rotag=r_tag
          npts=0
          p=0.d0
          solutions=0.d0
          orso=0.d0
          coordinates=0.d0
          tvalues=0.d0
          svalues=0.d0
          t=0.d0
          t1=0.d0
          tj=0.d0
          x1=0.d0
          y1=0.d0
          z1=0.d0
          zj=0.d0
          qq1=0.d0
          pp1=0.d0
          pp2=0.d0
          ss=0.d0
          sols=0.d0
          coords=0.d0
!     end of initialize: p,solutions,coordinates

!**********************************************************************
!**********************************************************************
!**********************************************************************



!****************END INTERSECTION UP DN***********************

          !calculate intersection of vector-line k with 1.

          do k=3,nn-1

          diff=-r_delta(1,1)*r_delta(k,2)+r_delta(k,1)*r_delta(1,2)
   
               if (diff==0.d0) then
                  write(*,*) 'Singular Knot',r_tag(1),r_tag(k)
                  write(*,*) (r_pts(1,j),j=1,3)
                  write(*,*) (r_pts(k,j),j=1,3)
                  goto 40
                endif
                     
                t(r_tag(1))=-((r_pts(k,1)-r_pts(1,1))*r_delta(k,2) &
                & -(r_pts(k,2)-r_pts(1,2))*r_delta(k,1) )/diff
                t(r_tag(k))=((r_pts(k,2)-r_pts(1,2))*r_delta(1,1) &
                & -(r_pts(k,1)-r_pts(1,1))*r_delta(1,2) )/diff

                t1=t(r_tag(1))
                tj=t(r_tag(k))

                !intersection points in plane
                !TRUE intersection point can be up or dn

             x1=r_pts(1,1)+r_delta(1,1)*t1
             y1=r_pts(1,2)+r_delta(1,2)*t1
             z1=r_pts(1,3)+r_delta(1,3)*t1

             xj=r_pts(k,1)+r_delta(k,1)*tj
             yj=r_pts(k,2)+r_delta(k,2)*tj
             zj=r_pts(k,3)+r_delta(k,3)*tj

             sig1=z1-zj 

             if(change>0)then
              !1: AA1,BB1
              !2: AA1,BB2
              !3: AA2,BB1
              !4: AA2,BB2

              if(segmentsinv==0)then
               if((r_pts(1,:).eqt.AA1).and.(r_pts(k,:).eqt.BB1)) segmentsinv=1
               if((r_pts(1,:).eqt.AA1).and.(r_pts(k,:).eqt.BB2)) segmentsinv=2
               if((r_pts(1,:).eqt.AA2).and.(r_pts(k,:).eqt.BB1)) segmentsinv=3
               if((r_pts(1,:).eqt.AA2).and.(r_pts(k,:).eqt.BB2)) segmentsinv=4

               if((r_pts(k,:).eqt.AA1).and.(r_pts(1,:).eqt.BB1)) segmentsinv=1
               if((r_pts(k,:).eqt.AA1).and.(r_pts(1,:).eqt.BB2)) segmentsinv=2
               if((r_pts(k,:).eqt.AA2).and.(r_pts(1,:).eqt.BB1)) segmentsinv=3
               if((r_pts(k,:).eqt.AA2).and.(r_pts(1,:).eqt.BB2)) segmentsinv=4

               if(segmentsinv>0)then
                sig1=-sig1
                howmuch=howmuch+1
               endif

              else

                okey=.false.
                if((r_pts(1,:).eqt.AA1).and.(r_pts(k,:).eqt.BB1).and.segmentsinv==1) okey=.true.
                if((r_pts(1,:).eqt.AA1).and.(r_pts(k,:).eqt.BB2).and.segmentsinv==2) okey=.true.
                if((r_pts(1,:).eqt.AA2).and.(r_pts(k,:).eqt.BB1).and.segmentsinv==3) okey=.true.
                if((r_pts(1,:).eqt.AA2).and.(r_pts(k,:).eqt.BB2).and.segmentsinv==4) okey=.true.

                if((r_pts(k,:).eqt.AA1).and.(r_pts(1,:).eqt.BB1).and.segmentsinv==1) okey=.true.
                if((r_pts(k,:).eqt.AA1).and.(r_pts(1,:).eqt.BB2).and.segmentsinv==2) okey=.true.
                if((r_pts(k,:).eqt.AA2).and.(r_pts(1,:).eqt.BB1).and.segmentsinv==3) okey=.true.
                if((r_pts(k,:).eqt.AA2).and.(r_pts(1,:).eqt.BB2).and.segmentsinv==4) okey=.true.

                if(okey)then
                 sig1=-sig1      
                 howmuch=howmuch+1
                endif

              endif 
             endif

               !CHANGE AVEC UNE INTERSECTION
               !sig1 dessus ou dessous
                if(sig1<0.) sn=-1  
                if(sig1==0.)sn=0
                if(sig1>0.) sn=1 

                qq1=r_delta(1,1)
                qq2=r_delta(1,2)
                pp1=r_delta(k,1)
                pp2=r_delta(k,2)

                !chiralite intersection gauche ou droite
                sgg1=qq1*pp2-qq2*pp1
                if(sgg1<0 ) snn=-1
                if(sgg1==0) snn= 0
                if(sgg1>0 ) snn= 1

               !se coupe ou a cote

               if(abs(t1)<1.d-8.or.abs(t1-1.d0)<1.d-8.or.abs(tj)<1.d-8.or.abs(tj-1.d0)<1.d-8)then
                if(rank==0) write(77,*) '***********************************************'
                if(rank==0) write(77,*) 'problem for intersection' 
                if(rank==0) writE(77,*) 'points too close to decide in homfly.f90'
                if(rank==0) write(77,*) 'r_pts1 : ' , r_pts(1,1:2)
                if(rank==0) write(77,*) 'r_pts1 : ' , r_pts(2,1:2)
                if(rank==0) write(77,*) 'r_pts2 : ' , r_pts(k,1:2)
                if(rank==0) write(77,*) 'r_pts2 : ' , r_pts(k,1)+r_delta(k,1),r_pts(k,2)+r_delta(k,2) 
                if(rank==0) write(77,*) 't1,t2  : ' , t1,tj
                if(rank==0) write(77,*) '***********************************************'
                return
               endif

               if (0.d0<t1.and.1.0>t1) then
               if (0.d0<tj.and.1.0>tj) then

                    if(sn==0) goto 40
                    failure=0                   
                    npts=npts+1
                    ss(1)=0.d0
                    ss(2)=snn
                    ss(3)=sn
                    ss(4)=x1
                    ss(5)=y1
                    ss(6)=t1
                    p(npts,1:6)=ss(1:6)
                    sols(1)=t1
                    sols(2)=tj
                    solutions(npts,1:2)=sols(1:2)
                    orso(1:npts,1)=solutions(1:npts,1)
                    coords(1)=x1
                    coords(2)=y1
                    coords(3)=z1
                    coordinates(npts,1:3)=coords(1:3)
                else
                  junk=0 !coupe a l exterieur
                endif
                else
                  junk=0 !coupe a l exterieur
                endif
                  
40              continue
                sn=1
                failure=2
                if (t1==1.0) goto 42 
                failure=0
                if (t1==0.d0) goto 42 
                failure=0
                if (tj==1.0) goto 42 
                failure=0
                if (tj==0.d0) goto 42 
                failure=0
                if (sn==0) goto 42 
                failure=0
            enddo
42         continue




!****************END INTERSECTION UP DN***********************


!**********************************************************************
!**********************************************************************
!**********************************************************************

          if (npts.ge.1998) then
           write(*,*) 'attention'
           write(*,*) 'npts=',npts
          endif


          !afficher p : 3updn 4,5:x,y
          !afficher vertex : r_pts(i,:) .........

          do j=1,npts
            do i=1,6
               p_app(lp,i)=p(j,i)
            enddo
            lp=lp+1
          enddo
     
             if(npts.gt.0) then
                lp=lp-1
             endif
             if (npts.gt.1) then      
                 do j=1,npts      
                  inter(j)=p_app(lp-npts+j,6)
                 enddo
               call ssort(npts,inter)
                   do j=1,npts 
                      do  k=1,npts 
                       if(inter(j).eq.p_app(lp-npts+k,6)) then
                           do m=1,6
                            p_end(lp-npts+j,m)=p_app(lp-npts+k,m)
                           enddo
                        endif
                     enddo
                   enddo
             elseif (npts.eq.1) then
                  do i=1,6
                    p_end(lp,i)=p_app(lp,i)
                  enddo
             endif
             continue
             if (npts.gt.0) then
               lp=lp+1
             endif
             do i=1,npts
                     tvalues(i)=p(i,6)
             enddo
             do i=1,npts
                     svalues(i)=tvalues(i)
             enddo
             if (npts.gt.1) then
              call ssort(npts,svalues)
              call ssort(npts,orso(1:npts,1))
             endif
                do i=1,npts
                   do j=1,npts
                    if (orso(i,1)==solutions(j,1))  then
                        orso(i,2)=solutions(j,2)
                    endif
                   enddo
                enddo

                do k=1,npts
                   ordsols(ind,1,k)=orso(k,1)
                   ordsols(ind,2,k)=orso(k,2)
                enddo

               if (ind.lt.stp) then
                   do j=1,(nn-1)
                    r_pts(j,1)=rotp(j+1,1)
                    r_pts(j,2)=rotp(j+1,2)
                    r_pts(j,3)=rotp(j+1,3)
                    r_delta(j,1)=rotd(j+1,1)
                    r_delta(j,2)=rotd(j+1,2)
                    r_delta(j,3)=rotd(j+1,3)
                   enddo
                    r_pts(nn,1)=rotp(1,1)
                    r_pts(nn,2)=rotp(1,2)
                    r_pts(nn,3)=rotp(1,3)
                    r_delta(nn,1)=rotd(1,1)
                    r_delta(nn,2)=rotd(1,2)
                    r_delta(nn,3)=rotd(1,3)
                     do i=1,(nn-1)
                      r_tag(i)=rotag(i+1)
                     enddo
                     r_tag(nn)=rotag(1)
                     goto 704 
                     !!!! REVIENT AU DEBUT APRES SWAP
                endif
703          continue


!*********************************************************
!*********************************************************
!*********************************************************
!*********************************************************
!*********************************************************


          if (failure.ne.0) goto 999
          junk=0

!     means not in Gen Posn


          !SI intersection SUR point gauche segment 
          l=1
          do i=1,step_size 
               do k=1,step_size-1 
                    if( ordsols(i,1,k)==0.d0) then
                         goto 71
                    else 
                         sutions(l,1)=ordsols(i,1,k)
                         sutions(l,2)=ordsols(i,2,k)
                         l=l+1
                    endif
               enddo
71          continue
          enddo


!         lmax greater that 1998, increase memory lmax
          lmax=l-1
          if (lmax.ge.2*maxcross) then
                   goto 999
           endif
           failure=1



 1902          format(' #      newknot  ')
               if (lmax.eq.0) then
               call putinString(STRING,' #      newknot ')
               call LineBreakString(STRING)

               !write(7,'(a10)') '1+1b1a1d1c '
               call putinString(STRING,'1+1b1a1d1c ')
               call LineBreakString(STRING)

               j_compt=j_compt+1
               !write(7,1795)
               write(STRING,1795)
               goto 999
               endif
                     
               failure=0
                     
              if (lmax.le.4) then
               call putinString(STRING,' #      newknot  ')
               call LineBreakString(STRING)
               !write(7,'(a10)') '1+1b1a1d1c'
               call putinString(STRING,'1+1b1a1d1c')
               call LineBreakString(STRING)
               j_compt=j_compt+1
               !write(7,1795)
               call putinString(STRING,' ')
               goto 999
             endif

!     means has 0 crossings

                !sutions contient t1,t2 parametre intersections croisement
                !noeud dans sens contraire, reverse:
                do k=1,lmax
                     reverse(k,1)=sutions(k,2)
                     reverse(k,2)=sutions(k,1)
                enddo
                
                !REVERSE??????
                lrev=1
                do k=1,lmax
                   do j=1,lmax
                           if(sutions(j,1)==reverse(k,1)) then
                                  tintix(lrev)=j 
                                  lrev= lrev+1 !reverse orientation
                           endif
                   enddo
                enddo

                do i=1,lmax
                     th(i,1)=p_end(i,2) !+-1
                     th(i,2)=p_end(i,3) !+-1!haut bas
                     th(i,3)=i          ! i
                     th(i,4)=tintix(i)  !index
                enddo


                do  i=1,lmax
                     ith(i,1)=nint(th(i,1))
                     ith(i,2)=nint(th(i,2))
                     ith(i,3)=nint(th(i,3))
                     ith(i,4)=nint(th(i,4))
                enddo

!     calculate the standard (Millett/Ewing convention) matrix 
!     of the link, starting with the Thistlethwaite Representation in ith()

                do k=1,lmax
                     if(mod(ith(k,3),2)==0) then
                        num(k)=ith(k,3)/2
                     else
                        num(k)=ith(k,4)/2
                     endif
                enddo 

                do k=1,lmax
                     if (ith(k,2).eq.1) then
                        ilet(k)='c'
                        inlet(k)=5
                        olet(k)='a'
                        onlet(k)=3
                     else
                        if (ith(k,1).eq.1) then
                                ilet(k)='d'
                                inlet(k)=6
                                olet(k)='b'
                                onlet(k)=4
                        else
                                ilet(k)='b'
                                inlet(k)=4
                                olet(k)='d'
                                onlet(k)=6
                        endif
                     endif
               enddo

!*******************************************
!*******************************************

!     start of the while5 loop

              write(ch2,'(i3)' ) num(1)
              write(ch4,'(i3)' ) num(lmax)

               do i=1,lmax
               if (ith(i,2).eq.1) then
                        j=num(i)
                    write(ch1,'(i3)' ) j
                    if (i.lt.lmax) then
                             write(ch5,'(i3)') num(i+1)
                           else
                             write(ch5,'(i3)') num(1)
                           endif

                           if (i.gt.1) then
                             write(ch6,'(i3)') num(i-1)
                           else
                             write(ch6,'(i3)') num(lmax)
                           endif 

                        cd(j,1)=ch1 
                        iminus=ith(i,1)*ith(i,2)
                        write(minus,'(i3)') iminus
                        cd(j,2)=minus 

                          if (cd(j,2)(1:1).eq.'-') then
                                  goto 1071
                             elseif (cd(j,2)(2:2).eq.'-') then
                                  goto 1071
                             elseif (cd(j,2)(3:3).eq.'-') then
                                  goto 1071
                             elseif (cd(j,2)(4:4).eq.'-') then
                                  goto 1071
                          else
                                  cd(j,2)='+'
                                  goto 1072
                          endif

1071                          cd(j,2)='-'

1072                    continue

                        if (i.eq.lmax) then
                                cd(j,3)=ch2//ilet(1)
                        else 
                                cd(j,3)=ch5//ilet(i+1)
                    endif
               
                        if (i.eq.1) then
                                cd(j,5)=ch4//olet(lmax)
                        else
                                cd(j,5)=ch6//olet(i-1)
                        endif
                   
                        if (i.eq.lmax) then
                                cd(num(1),inlet(1))=ch1//'a'
                        else
                                cd(num(i+1),inlet(i+1))=ch1//'a'
                        endif
                        if (i.eq.1) then
                                cd(num(lmax),onlet(lmax))=ch1//'c'
                        else
                                cd(num(i-1),onlet(i-1))=ch1//'c'
                        endif
                      else 
                        if ((ith(i,1)*ith(i,2)).eq.1) then
                             j=num(i)

                             write(ch1,'(i3)') j

                                if (i.lt.lmax) then
                                  write(ch5,'(i3)') num(i+1)
                                else
                                  write(ch5,'(i3)') num(1)
                                endif

                                if (i.gt.1) then
                                  write(ch6,'(i3)') num(i-1)
                                else
                                  write(ch6,'(i3)') num(lmax)
                                endif

                                if (i.eq.lmax) then
                                       cd(j,6)=ch2//ilet(1)
                                else
                                       cd(j,6)=ch5//ilet(i+1)
                                endif
                                     
                                if (i.eq.1) then
                                       cd(j,4)=ch4//olet(lmax)
                                else
                                       cd(j,4)=ch6//olet(i-1)
                                endif
                                     
                                if (i.eq.lmax) then
                                       cd(num(1),inlet(1))=ch1//'d'
                                else
                                       cd(num(i+1),inlet(i+1))=ch1//'d'
                                endif
                                
                                if (i.eq.1) then
                                       cd(num(lmax),onlet(lmax))=ch1//'b'
                                else
                                       cd(num(i-1),onlet(i-1))=ch1//'b'
                                endif
                                
                        else 
                                j=num(i)

                             write(ch1,'(i3)' ) j

                                if (i.lt.lmax) then
                                  write(ch5,'(i3)') num(i+1)
                                else
                                  write(ch5,'(i3)') num(1)
                                endif

                                if (i.gt.1) then
                                  write(ch6,'(i3)') num(i-1)
                                else
                                  write(ch6,'(i3)') num(lmax)
                                endif

                                if (i.eq.lmax) then
                                       cd(j,4)=ch2//ilet(1)
                                else
                                       cd(j,4)=ch5//ilet(i+1)
                                endif
                                     
                                if (i.eq.1) then
                                       cd(j,6)=ch4//olet(lmax)
                                else
                                       cd(j,6)=ch6//olet(i-1)
                                endif
                                     
                                if (i.eq.lmax) then
                                       cd(num(1),inlet(1))=ch1//'b'
                                else
                                       cd(num(i+1),inlet(i+1))=ch1//'b'
                                endif
                                     
                                if (i.eq.1) then
                                       cd(num(lmax),onlet(lmax))=ch1//'d'
                                else
                                       cd(num(i-1),onlet(i-1))=ch1//'d'
                                endif
                        endif 
                      endif
               enddo
  
!     end of while#5 loop
!*******************************************
!*******************************************     

 
            if (failure.ne.0) then
                  goto 999
                else
                  junk=0
                endif

!     the next section prints out the matrix and calculates
!     the writhe start of while #6 loop

                 do ist=1,lmax/2
                      do jst=1,6
                           cdwr(ist,jst,1)=cd(ist,jst)(1:1)
                           cdwr(ist,jst,2)=cd(ist,jst)(2:2) 
                           cdwr(ist,jst,3)=cd(ist,jst)(3:3)
                           cdwr(ist,jst,4)=cd(ist,jst)(4:4)
                      enddo
                 enddo

                totalcross=totalcross+lmax/2
                call putinString(STRING,' #      newknot  ')
                call LineBreakString(STRING)
                do ist=1,lmax/2
                   ia=1
                      do jst=1,6
                           do  kst=1,4 
                                if (cdwr(ist,jst,kst).ne.' ') then
                                    lett(ist,ia)=cdwr(ist,jst,kst)
                                    ia=ia+1
                                endif
                           enddo
                      enddo
               ib=ia-1
               !write(7,1793) (lett(ist,lst),lst=1,ib)

               call putinString(STRING, (/ (lett(ist,lst),lst=1,ib) /)  )
               call LineBreakString(STRING)
1793           format(24a1)
               enddo
     
                !write(7,1795)
                call putinString(STRING,' ')
1795            format(a1)

                 if (failure.eq.0) then
                    failure=5 
                 else
                    junk=0
                 endif

!     next section is for when the knot has no crossings
!     or is not in general position

            if (failure.lt.5) then
                if (failure.eq.1) then
                 write(*,*) 'NO CROSSINGS'
                else
                 junk=0
                endif

               if (failure.eq.2) then
                      write(*,*) 'Not in General Position'
               else
                     junk=0
               endif

               if (failure.eq.3) then
                    write(*,*) 'The same point was chosen twice. Execution stopped'
               else
                    junk=0
               endif
           endif


999  continue
1011 continue
1001 continue

     if(allocated(r_size))deallocate(r_size)
     sizeanalysis2=1
     call finalizeString(STRING)
     call LineBreakString(STRING)
     return
end function

!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************

subroutine delta2(r,long,size)
implicit none
integer::i,j,long,icount
integer,intent(out)::size
real(8),dimension(3,3)::matrix_s,matrix_t,matrix_u,matrix_0
real(8),dimension(long,3)::r
real(8)::u,s,t
logical::boolean,check

icount=1;size=0;i=1;j=0

if (long>=5) then
     boolean =.true.
     do while (boolean)
          if (i>=long-2) then
               if (i>=long-4) then
                    j=i+2
                    do while (j>=long-2)
                         j=j+1
                         matrix_s=0.d0
                         matrix_t=0.d0
                         matrix_u=0.d0
                         matrix_0=0.d0
     
                         matrix_t(1,:)=(/r(j,1)-r(i,1),r(i+1,1)-r(i,1),r(j,1)-r(j+1,1)/)
                         matrix_t(2,:)=(/r(j,2)-r(i,2),r(i+1,2)-r(i,2),r(j,2)-r(j+1,2)/)
                         matrix_t(3,:)=(/r(j,3)-r(i,3),r(i+1,3)-r(i,3),r(j,3)-r(j+1,3)/)

                         matrix_u(1,:)=(/r(i+2,1)-r(i,1),r(j,1)-r(i,1),r(j,1)-r(j+1,1)/)
                         matrix_u(2,:)=(/r(i+2,2)-r(i,2),r(j,2)-r(i,2),r(j,2)-r(j+1,2)/)
                         matrix_u(3,:)=(/r(i+2,3)-r(i,3),r(j,3)-r(i,3),r(j,3)-r(j+1,3)/)
          
                         matrix_s(1,:)=(/r(i+2,1)-r(i,1),r(i+1,1)-r(i,1),r(j,1)-r(i,1)/)
                         matrix_s(2,:)=(/r(i+2,2)-r(i,2),r(i+1,2)-r(i,2),r(j,2)-r(i,2)/)
                         matrix_s(3,:)=(/r(i+2,3)-r(i,3),r(i+1,3)-r(i,3),r(j,3)-r(i,3)/)
          
                         matrix_0(1,:)=(/r(i+2,1)-r(i,1),r(i+1,1)-r(i,1),r(j,1)-r(j+1,1)/)
                         matrix_0(2,:)=(/r(i+2,2)-r(i,2),r(i+1,2)-r(i,2),r(j,2)-r(j+1,2)/)
                         matrix_0(3,:)=(/r(i+2,3)-r(i,3),r(i+1,3)-r(i,3),r(j,3)-r(j+1,3)/)
          
                         if (determinant(matrix_0).ne.0.d0) then
                              t=determinant(matrix_t)/determinant(matrix_0)
                              u=determinant(matrix_u)/determinant(matrix_0)
                              s=determinant(matrix_s)/determinant(matrix_0)
          
                              if (s.ge.0.d0.and.s.le.1.d0&
                & .and.u.ge.0.d0.and.u.le.1.d0.and.t.ge.0.d0.and.t.le.(1.d0-u)) boolean=.false.     
                         endif
                    enddo
               endif
                
               if (i.ge.2) then
                    if (i.le.long-3) then
                         matrix_s=0.d0
                         matrix_t=0.d0
                         matrix_u=0.d0
                         matrix_0=0.d0
                         
                         matrix_t(1,:)=(/r(long,1)-r(i,1),r(i+1,1)-r(i,1),r(long,1)-r(1,1)/)
                         matrix_t(2,:)=(/r(long,2)-r(i,2),r(i+1,2)-r(i,2),r(long,2)-r(1,2)/)
                         matrix_t(3,:)=(/r(long,3)-r(i,3),r(i+1,3)-r(i,3),r(long,3)-r(1,3)/)

                         matrix_u(1,:)=(/r(i+2,1)-r(i,1),r(long,1)-r(i,1),r(long,1)-r(1,1)/)
                         matrix_u(2,:)=(/r(i+2,2)-r(i,2),r(long,2)-r(i,2),r(long,2)-r(1,2)/)
                         matrix_u(3,:)=(/r(i+2,3)-r(i,3),r(long,3)-r(i,3),r(long,3)-r(1,3)/)
          
                         matrix_s(1,:)=(/r(i+2,1)-r(i,1),r(i+1,1)-r(i,1),r(long,1)-r(i,1)/)
                         matrix_s(2,:)=(/r(i+2,2)-r(i,2),r(i+1,2)-r(i,2),r(long,2)-r(i,2)/)
                         matrix_s(3,:)=(/r(i+2,3)-r(i,3),r(i+1,3)-r(i,3),r(long,3)-r(i,3)/)
          
                         matrix_0(1,:)=(/r(i+2,1)-r(i,1),r(i+1,1)-r(i,1),r(long,1)-r(1,1)/)
                         matrix_0(2,:)=(/r(i+2,2)-r(i,2),r(i+1,2)-r(i,2),r(long,2)-r(1,2)/)
                         matrix_0(3,:)=(/r(i+2,3)-r(i,3),r(i+1,3)-r(i,3),r(long,3)-r(1,3)/)
          
                if(determinant(matrix_0).ne.0.d0) then
                  t=determinant(matrix_t)/determinant(matrix_0)
                  u=determinant(matrix_u)/determinant(matrix_0)
                  s=determinant(matrix_s)/determinant(matrix_0)
          
                if(s.ge.0.d0.and.s.le.1.d0.and.&
            & u.ge.0.d0.and.u.le.1.d0.and.t.ge.0.d0.and.t.le.(1.d0-u)) &
            & boolean=.false.
                 endif
                endif
                    
                    if (i.gt.2) then
                         do j=1,i-2
                              matrix_s=0.d0
                              matrix_t=0.d0
                              matrix_u=0.d0
                              matrix_0=0.d0
                                   
                              matrix_t(1,:)=(/r(j,1)-r(i,1),r(i+1,1)-r(i,1),r(j,1)-r(j+1,1)/)
                              matrix_t(2,:)=(/r(j,2)-r(i,2),r(i+1,2)-r(i,2),r(j,2)-r(j+1,2)/)
                              matrix_t(3,:)=(/r(j,3)-r(i,3),r(i+1,3)-r(i,3),r(j,3)-r(j+1,3)/)

                              matrix_u(1,:)=(/r(i+2,1)-r(i,1),r(j,1)-r(i,1),r(j,1)-r(j+1,1)/)
                              matrix_u(2,:)=(/r(i+2,2)-r(i,2),r(j,2)-r(i,2),r(j,2)-r(j+1,2)/)
                              matrix_u(3,:)=(/r(i+2,3)-r(i,3),r(j,3)-r(i,3),r(j,3)-r(j+1,3)/)
          
                              matrix_s(1,:)=(/r(i+2,1)-r(i,1),r(i+1,1)-r(i,1),r(j,1)-r(i,1)/)
                              matrix_s(2,:)=(/r(i+2,2)-r(i,2),r(i+1,2)-r(i,2),r(j,2)-r(i,2)/)
                              matrix_s(3,:)=(/r(i+2,3)-r(i,3),r(i+1,3)-r(i,3),r(j,3)-r(i,3)/)
          
                              matrix_0(1,:)=(/r(i+2,1)-r(i,1),r(i+1,1)-r(i,1),r(j,1)-r(j+1,1)/)
                              matrix_0(2,:)=(/r(i+2,2)-r(i,2),r(i+1,2)-r(i,2),r(j,2)-r(j+1,2)/)
                              matrix_0(3,:)=(/r(i+2,3)-r(i,3),r(i+1,3)-r(i,3),r(j,3)-r(j+1,3)/)
          
                      if (determinant(matrix_0).ne.0.d0) then
                       t=determinant(matrix_t)/determinant(matrix_0)
                       u=determinant(matrix_u)/determinant(matrix_0)
                       s=determinant(matrix_s)/determinant(matrix_0)
                       if (s.ge.0.d0.and.s.le.1.d0.and.u.ge.0.d0.and.u&
                       & .le.1.d0.and.t.ge.0.d0.and.t.le.(1.d0-u)) &
                       & boolean=.false.
                      endif

                   enddo
                         
                  endif


               endif
          

               if (boolean.eqv..false.) then
                    if (i.lt.long-2) then
                         rback(icount,:)=r(i,:)
                         icount=icount+1
                         i=i+1
                         size=size+1
                         boolean=.true.
                    else
                         rback(icount,:)=r(long-2,:)
                         icount=icount+1
                         rback(icount,:)=r(long-1,:)
                         icount=icount+1
                         rback(icount,:)=r(long,:)
                         icount=icount+1
                         size=size+3
                    endif
               else
                    
                    if (i.lt.long-2) then
                         boolean=.false.
                         rback(icount,:)=r(i,:)
                         icount=icount+1
                         size=size+1
                         do j=i+2,long
                              rback(icount,:)=r(j,:)
                              icount=icount+1
                              size=size+1
                         enddo
                    else
                         boolean=.false.
                         rback(icount,:)=r(long-2,:)
                         icount=icount+1
                         rback(icount,:)=r(long,:)
                         icount=icount+1
                         size=size+2
                    endif
               endif
          endif
     enddo
else
     do i=1,long
          rback(icount,:)=r(i,:)
          icount=icount+1
          size=size+1
     enddo
endif

return
end subroutine



!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************

end module

