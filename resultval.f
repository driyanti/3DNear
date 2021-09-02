      SUBROUTINE resultval(omegacnt,gv)

      IMPLICIT NONE
      
      INCLUDE 'bounds.h'            
      INCLUDE 'num.h'           
      INCLUDE 'param.h'

      INTEGER ii,jj,cnt,kcnt,omegacnt,dcnt
      
      REAL*8  fr(xdim),fi(xdim)
      
      COMPLEX greens_p(xdim,xdim,dim,dim),wlin(ns)
      
      COMPLEX gv(max_rcv,max_rcv,maxksamp,dim,dim)
                                                     
c      CALL kxkyvalue     
         
      CALL input
      
c      CALL gridvalue
            
c      CALL linwave(wlin)

      DO dcnt=1,ndat      
            
      CALL matgreen(omegacnt,sx(dcnt),sy(dcnt),greens_p,
     .              rcvpos(3),sz(dcnt))
            
       do kcnt=1,dim              
        do cnt=1,dim
       
        do ii=1,xdim
         do jj=1,xdim 
           fr(jj) = real(greens_p(ii,jj,cnt,kcnt))
           fi(jj) = imag(greens_p(ii,jj,cnt,kcnt))               
         end do
        
         call fftpot(2,mlog,xdim,fr,fi)
         
         do jj=1,max_rcv
           greens_p(ii,jj,cnt,kcnt) = cmplx(fr(jj),fi(jj))         
         end do
        
        end do 
       
        do jj=1,max_rcv
	 
         do ii=1,xdim
           fr(ii)= real(greens_p(ii,jj,cnt,kcnt))
	   fi(ii)= imag(greens_p(ii,jj,cnt,kcnt))
         end do
    	   
         call fftpot(2,mlog,xdim,fr,fi)
	 
	 do ii=1,max_rcv
           
           greens_p(ii,jj,cnt,kcnt) = 1./float(xdim)*
     &                                cmplx(fr(ii),fi(ii))
     
	   gv(ii,jj,dcnt,cnt,kcnt) = greens_p(ii,jj,cnt,kcnt) 
            
         end do	  	         
        end do
       
         end do  !end cnt
        end do !end kcnt              
       end do !end ndat
        
      end subroutine      
cccccccccccccccccccccccccccccccccccccccccccccccc
