        subroutine matBH(omegacnt,newhmat)

        implicit none

        INCLUDE 'bounds.h'
        INCLUDE 'num.h'
        INCLUDE 'param.h'
        
        integer nclm,omegacnt,cnt

        complex hmat(maxksamp,dim),newhmat(maxksamp)
                          
        do cnt=1,ndat
          
         CALL matH(omegacnt,cnt,hmat)

         do nclm=(cnt-1)*dim+1,cnt*dim
          newhmat(nclm) = hmat(cnt,nclm-(cnt-1)*dim)
         end do         
        
        end do
          
        end subroutine
        
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc            
         
