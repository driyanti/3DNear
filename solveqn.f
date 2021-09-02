      SUBROUTINE solveeqnval(omegacnt,utot)
        
      IMPLICIT NONE

      INCLUDE 'bounds.h'               
      INCLUDE 'num.h'               
      INCLUDE 'param.h'	
	
      REAL*8     rt,rcond
      
      COMPLEX    newh(maxksamp),newg(maxksamp,maxksamp),utot(maxksamp)
                        
      COMPLEX*16 b(maxksamp),bb(maxksamp,maxksamp),solv(maxksamp)
     &           ,zz(maxksamp)
        
      INTEGER    ntot,ii,jj,ipvt(maxksamp),omegacnt,lda

      DATA lda/maxksamp/
                         
      CALL input 
        
      ntot=dim*ndat

      CALL matBH(omegacnt,newh)  
		        
      CALL matBG(omegacnt,newg)
               
      DO ii=1,ntot	     	   
       DO jj=1,ntot	      
        b(jj) = newh(jj)		          
	bb(ii,jj) = newg(ii,jj)
c        print*,'hg',ii,jj,newh(jj),newg(ii,jj)              
       END DO	   
      END DO
                
      CALL zgeco(bb,lda,ntot,ipvt,rcond,zz)
	
      WRITE(*,*)'rcond',rcond
	
      rt = 1.0 + rcond	
	
      IF (rt .EQ. 1.0) GO TO 90
	
      CALL zgesl(bb,lda,ntot,ipvt,b,0)

      DO ii=1,ntot       
       solv(ii) = b(ii) 
       utot(ii) = solv(ii)          
c       print*,'solve',ii,solv(ii)       
      END DO
	      	
90    WRITE(*,99)

99    FORMAT(40H matrix is singular to working precision)
             
      END SUBROUTINE      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
