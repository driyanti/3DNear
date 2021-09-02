       SUBROUTINE matBG(omegacnt,newg)
             
       IMPLICIT NONE

       INCLUDE 'bounds.h'
       INCLUDE 'num.h'
       INCLUDE 'param.h'
        
       INTEGER ii,jj,cnt,kcnt,omegacnt,ntot
                
       COMPLEX greens_p(maxksamp,maxksamp,dim,dim)
     .         ,newg(maxksamp,maxksamp)
       
       ntot = dim*ndat
        
       DO cnt=1,dim
        DO kcnt=1,dim

        CALL mkgmat(omegacnt,cnt,kcnt,newg)        
        
        DO ii=1,ndat
         DO jj=1,ndat
         
         greens_p(ii,jj,cnt,kcnt)=newg(ii,jj)
c          print*,'newg',ii,jj,kcnt,cnt,newg(ii,jj)            

         END DO
        END DO
        
        END DO
       END DO
              
       DO cnt=1,ndat
        DO kcnt=1,ndat 
         DO ii=(cnt-1)*dim+1,cnt*dim
          DO jj=(kcnt-1)*dim+1,kcnt*dim 
         
           newg(ii,jj) = greens_p(cnt,kcnt,ii-(cnt-1)
     .                   *dim,jj-(kcnt-1)*dim)

          END DO
         END DO
        END DO
       END DO 

       END SUBROUTINE            
          
