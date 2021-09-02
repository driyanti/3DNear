      SUBROUTINE mkgfree(omegacnt,rx,ry,rz,DDx,rsum,
     .                   green_f,rcvlevel,srclevel)
                  
      IMPLICIT NONE
      
      INCLUDE 'bounds.h'
      INCLUDE 'num.h'
      INCLUDE 'param.h'
      
      REAL*8  rcvlevel,srclevel,omega,rsum,rx,ry,rz,DDx
                
      COMPLEX green_kf(dim,dim),green_f(dim,dim),matDS(dim*dim),
     .        matDP(dim*dim),kS,kP,factorS,mkcoshS,mksinhS,matSS,mu
     .        ,factorP,mkcoshP,mksinhP,matPP 
      
      INTEGER omegacnt
           
      CALL input
      	       
      omega = 2.0*pi*float(omegacnt)
	                                          	                  
      mu = beta(1)**2*rho(1)
                  
cccc----starting compute the free green's tensor------
c-cccccc--DDx=min(dx,dy,dz)-------
      
      kS = omega/beta(1)
      
      kP = omega/alpha(1)
                  
      factorS = 1d0/3d0*pi*(kS*DDx)**2*rsum
      
      mkcoshS = 0.5d0*(exp(ci*kS*DDx)+exp(-ci*kS*DDx))
      
      mksinhS = 0.5d0*(exp(ci*kS*DDx)-exp(-ci*kS*DDx))
      
      matSS = exp(ci*kS*rsum)*(mksinhS/(0.5d0*ci*kS)-mkcoshS)/factorS
      
      factorP = 1d0/3d0*pi*(kP*DDx)**2*rsum
      
      mkcoshP = 0.5d0*(exp(ci*kP*DDx)+exp(-ci*kP*DDx))
      
      mksinhP = 0.5d0*(exp(ci*kP*DDx)-exp(-ci*kP*DDx))
      
      matPP = exp(ci*kP*rsum)*(mksinhP/(0.5d0*ci*kP)-mkcoshP)/factorP 
      
      matDS(1)=exp(ci*kS*rsum)/(-ci*kS*rsum**3)*(3d0*rx**2/rsum**2-1d0 
     .         -ci*kS*rsum*(3d0*rx**2/rsum**2-1d0) + (-ci*kS*rx)**2)
     
      matDS(2)=exp(ci*kS*rsum)/(-ci*kS*rsum**3)*(3d0*rx*ry/rsum**2-1d0 
     .         -ci*kS*rsum*(3d0*rx*ry/rsum**2-1d0) + (-ci*kS)**2*rx*ry)

      matDS(3)=exp(ci*kS*rsum)/(-ci*kS*rsum**3)*(3d0*rx*rz/rsum**2-1d0 
     .         -ci*kS*rsum*(3d0*rx*rz/rsum**2-1d0) + (-ci*kS)**2*rx*rz)

      matDS(4) = matDS(2)
        
      matDS(5)=exp(ci*kS*rsum)/(-ci*kS*rsum**3)*(3d0*ry**2/rsum**2-1d0 
     .        - ci*kS*rsum*(3d0*ry**2/rsum**2-1d0) + (-ci*kS*ry)**2)
     
      matDS(6)=exp(ci*kS*rsum)/(-ci*kS*rsum**3)*(3d0*ry*rz/rsum**2-1d0 
     .       - ci*kS*rsum*(3d0*ry*rz/rsum**2-1d0) + (-ci*kS)**2*ry*rz)
  
      matDS(7) = matDS(3)
      
      matDS(8) = matDS(6)
      
      matDS(9)=exp(ci*kS*rsum)/(-ci*kS*rsum**3)*(3d0*rz**2/rsum**2-1d0 
     .         -ci*kS*rsum*(3d0*rz**2/rsum**2-1d0) + (-ci*kS*rz)**2)
       
      matDP(1)=exp(ci*kP*rsum)/(-ci*kP*rsum**3)*(3d0*rx**2/rsum**2-1d0 
     .         -ci*kP*rsum*(3d0*rx**2/rsum**2-1d0) + (-ci*kP*rx)**2)
     
      matDP(2)= exp(ci*kP*rsum)/(-ci*kP*rsum**3)*(3d0*rx*ry/rsum**2-1d0 
     .         -ci*kP*rsum*(3d0*rx*ry/rsum**2-1d0) + (-ci*kP)**2*rx*ry)

      matDP(3)= exp(ci*kP*rsum)/(-ci*kP*rsum**3)*(3d0*rx*rz/rsum**2-1d0 
     .        -ci*kP*rsum*(3d0*rx*rz/rsum**2-1d0) + (-ci*kP)**2*rx*rz)

      matDP(4) = matDP(2)
        
      matDP(5)=exp(ci*kP*rsum)/(-ci*kP*rsum**3)*(3d0*ry**2/rsum**2-1d0 
     .         -ci*kP*rsum*(3d0*ry**2/rsum**2-1d0) + (-ci*kP*ry)**2)
     
      matDP(6)=exp(ci*kP*rsum)/(-ci*kP*rsum**3)*(3d0*ry*rz/rsum**2-1d0 
     .       - ci*kP*rsum*(3d0*ry*rz/rsum**2-1d0) + (-ci*kP)**2*ry*rz)
  
      matDP(7) = matDP(3)
      
      matDP(8) = matDP(6)
      
      matDP(9)=exp(ci*kP*rsum)/(-ci*kP*rsum**3)*(3d0*rz**2/rsum**2-1d0 
     .         -ci*kP*rsum*(3d0*rz**2/rsum**2-1d0) + (-ci*kP*rz)**2)
      
      green_f(1,1) = matSS/(mu*4*pi) - (matDS(1) - matDP(1))
     .               /(4.*pi*rho(1)*omega)
      
      green_f(1,2) = -(matDS(2) - matDP(2))/(4.*pi*rho(1)*omega)
      
      green_f(1,3) = -(matDS(3) - matDP(3))/(4.*pi*rho(1)*omega)
      
      green_f(2,1) = green_f(1,2)
      
      green_f(2,2) = matSS/(4*pi*mu) - (matDS(5) - matDP(5))
     .               /(4.*pi*rho(1)*omega)
      
      green_f(2,3) = -(matDS(6) - matDP(6))/(4.*pi*rho(1)*omega)
      
      green_f(3,1) = green_f(1,3)
      
      green_f(3,2) = green_f(2,3)
      
      green_f(3,3) = matSS/(4.*pi*mu) - (matDS(9) - matDP(9))/
     .               (4.*pi*rho(1)*omega)
                                              
      END SUBROUTINE
