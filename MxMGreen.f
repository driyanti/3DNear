      SUBROUTINE matgreenz(omegacnt,r,xi,greens_p,greens_h,rcvlevel,
     .                     srclevel,delta_k)
      
c
c author     :  Auke Ditzel
c modified by: Christina Dwi Riyanti 
c last changes:   January 2003
c
c subroutine to calculate Greens function in slowness domain(radial domain)
c 
c September 2000: implementation force direction (mksource)
c 
      IMPLICIT NONE

      INCLUDE 'bounds.h'
      INCLUDE 'num.h'
      INCLUDE 'param.h'
      
      INTEGER  kcnt,layer,wavecntx,ii,omegacnt,azm
      
      REAL*8   omega,pabs,psum,rcvlevel,srclevel   
      
      REAL*8   J0,J1,J2,length,delta_k,k_x(ndim),fnu,r,p,xi   
      
      COMPLEX*16  matD(4,4),matDh(2,2),Su(3,3),Sd(3,3),
     .		  v_src_min(6),v_src_plus(6),v_rcv(6),b_rcv(6,3),hlpvar
          
      COMPLEX     greens_p(ndim,dim,dim),greens_h(ndim,dim,dim),arg
     .            ,cy(10),gr,gxi,green_h(dim,dim)
      
      INTEGER ierr,nz
            
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c v_src_min 	  	wavefield above src level
c
c v_src_plus	  	wavefield below src level
c
c matDcombi 	  	combination of matrices D
c     	    	  	(P-SV) and (SH)
c 
c SU,SD     	  	scattering matrices (Kennett)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c do NOT calculate for omega==0.0
c
c      nk_totpos = 0
      
      omega = 2d0*pi*float(omegacnt)
       
      IF (omega.eq.0.d0) GOTO 750
      
      length  = pi
      delta_k = length/float(ndim)
      
c length_k_1=length_k_2 denotes the length interval in the wavenumber domain c(k_x,k_y)
           
      do ii=1,ndim
	k_x(ii) = float(ii)*delta_k
      end do
                                    
      DO wavecntx = 1,ndim  
          
          k_x(wavecntx) = float(wavecntx-1)*delta_k                  
    	 
          psum = k_x(wavecntx)/omega
                              
          p=psum**2
         	        
	  pabs = sqrt(p)

c
c initialize q_alpha and q_beta with appropriate signs
c sign of q_alpha and q_beta should correspond to the fourier transform
c see De Hoop (waves in elastic media, dissipative medium and our fourier trnsf)
c
	DO layer = 1,nlay+1
          hlpvar        = sqrt(-psum**2 + 1.d0/alpha(layer)**2)
	  qalpha(layer) = cmplx(abs(dble(hlpvar)),dble(-ci*hlpvar)) 
          hlpvar        = sqrt(-psum**2 + 1.d0/beta(layer)**2)
	  qbeta(layer)  = cmplx(abs(dble(hlpvar)),dble(-ci*hlpvar)) 
          
        END DO	  ! end layer loop

        CALL mkgfreepz(omegacnt,r,xi,psum,green_h,rcvlevel,srclevel)
         
         greens_h(wavecntx,1,1) = green_h(1,1)
         
         greens_h(wavecntx,1,2) = green_h(1,2)
         
         greens_h(wavecntx,1,3) = green_h(1,3)
         
         greens_h(wavecntx,2,1) = green_h(2,1)
         
         greens_h(wavecntx,3,1) = green_h(3,1)
         
         greens_h(wavecntx,2,2) = green_h(2,2)
         
         greens_h(wavecntx,2,3) = green_h(2,3)
         
         greens_h(wavecntx,3,3) = green_h(3,3)
         
         greens_h(wavecntx,3,2) = green_h(3,2)  

c construct matDcombi matrix to be used later
c
 	CALL mkD(pabs,matD,rcvlayer)    
                  	
        CALL mkDh(matDh,rcvlayer)
        
c
c matrices SU and SD are constructed
c matrices RD, RU, TD, TU are constructed as well
c for each layer
c
                  
	CALL matrix_proloog(pabs,omega,Su,Sd,srclevel)
         
	DO kcnt = 1,dim	! for each force direction
        
	  k_x(wavecntx) = float(wavecntx-1)*delta_k     
        
          psum = k_x(wavecntx)/omega
             
          p = psum**2
                     
c at source level the field in terms of up- and downgoing waves 
c is constructed.
c output: v_src_min and v_src_plus represent field above and below 
c the source level
c kcnt : direction force
c                                       
         DO azm=1,3
                
	 CALL mksource(psum,omega,Su,Sd,v_src_min,v_src_plus,kcnt
     .                 ,azm)
                        
         CALL Gcalcfield(psum,omega,v_src_min,v_src_plus,v_rcv,
     .         rcvlevel,srclevel)
     
c
c DISPLACEMENTS
c

     	  b_rcv(1,azm) = matD(1,3)*(v_rcv(4) - v_rcv(1)) +
     .                   matD(1,4)*(v_rcv(2) + v_rcv(5))

	  b_rcv(2,azm) = matD(2,1)*(v_rcv(4) + v_rcv(1)) +
     .                   matD(2,4)*(v_rcv(5) - v_rcv(2))

	  b_rcv(3,azm) = matDh(1,1)*(v_rcv(3) + v_rcv(6))
          
         END DO 
         
c
c STRESSES
c
	
!	  b_rcv(4) = matD(3,1)*(v_rcv(4) + v_rcv(1)) +
!     .	  	     matD(3,4)*(v_rcv(5) - v_rcv(2))

!	  b_rcv(5) = matD(4,3)*(v_rcv(4) - v_rcv(1)) +
!     .	  	     matD(4,4)*(v_rcv(5) + v_rcv(2))

!	  b_rcv(6) = matDh(2,2)*(-v_rcv(3) + v_rcv(6))

c
c Transformation (U,V,W) to (u,v,w) (slowness domain) 
c
         
         IF (r .eq. 0. .or. xi .eq. 0.) THEN
          
          r = sqrt(dx*dy/(2.*pi))
          
          arg = omega*psum*r
          
          fnu = 0.0
          
          CALL CBESH (arg,fnu,1,1,3,cy,nz,ierr)  
                          
          J1 = real(cy(2)) 
              
          IF (psum .eq. 0.d0) THEN
                           
          greens_p(wavecntx,kcnt,1) = J1*(b_rcv(2,2)-b_rcv(2,3)
     .                + ci*(b_rcv(3,2)+ b_rcv(3,3)))*r*pi/(2.*pi)


c          greens_p(wavecntx,kcnt,2) =  J1*(ci*(b_rcv(2,2)-b_rcv(2,3))
          greens_p(wavecntx,kcnt,2) =  J1*(ci*(b_rcv(2,2)+b_rcv(2,3))
     .                - (b_rcv(3,2)- b_rcv(3,3)))*r*pi/(2.*pi)
     
          greens_p(wavecntx,kcnt,3) =  r*2.*pi*J1*b_rcv(1,1)/(2.*pi)
                                        
	  ELSE !NOT ((p2.eq.0.d0).and.(p1.eq.0.d0))

          greens_p(wavecntx,kcnt,1) = J1*(b_rcv(2,2)-b_rcv(2,3)
     .     + ci*(b_rcv(3,2)+ b_rcv(3,3)))*r*pi/(2.*pi)
          
          greens_p(wavecntx,kcnt,2) = J1*(ci*(b_rcv(2,2)-b_rcv(2,3))
     .     - (b_rcv(3,2)- b_rcv(3,3)))*r*pi/(2.*pi)
                    
          greens_p(wavecntx,kcnt,3) = r*2.*pi*J1*b_rcv(1,1)/(2.*pi)
              
          END IF              
         
         ELSE
     
          arg = omega*psum*r
          
          fnu = 0.0
          
          CALL CBESH (arg,fnu,1,1,3,cy,nz,ierr)  
          
          J0 = real(cy(1))
                
          J1 = real(cy(2)) 
          
          J2 = real(cy(3))

          IF (psum .eq. 0.d0) THEN
                           
          greens_p(wavecntx,kcnt,1) = (psum*(-b_rcv(2,1))*J1+ 
     .           exp(ci*xi)*(b_rcv(2,2)*((-J2*psum)+1./(omega*r)*J1)  
     .           + b_rcv(3,2)*ci/(omega*r)*J1)+ 
     .           exp(-ci*xi)*(b_rcv(2,3)*(-J0*psum+1./(omega*r)*J1) + 
     .           b_rcv(3,3)*(ci)/(omega*r)*J1))*omega
          
          gr = greens_p(wavecntx,kcnt,1)     
     
          greens_p(wavecntx,kcnt,2) =  (psum*b_rcv(3,1)*J1 + 
     .           exp(ci*xi)*(-b_rcv(3,2)*(-J2*psum+1./(omega*r)*J1) + 
     .           b_rcv(2,2)*ci/(omega*r)*J1) +
     .           exp(-ci*xi)*(-b_rcv(3,3)*(-J0*psum+1./(omega*r)*J1)+
     .           b_rcv(2,3)*(ci)/(omega*r)*J1))*omega
     
          gxi = greens_p(wavecntx,kcnt,2)
     
          greens_p(wavecntx,kcnt,3) =  psum*(b_rcv(1,1)*J0
     .           + b_rcv(1,2)*J1*exp(ci*xi)-b_rcv(1,3)*J1*exp(-ci*xi))
     .           *omega/(2.*pi)
     
          greens_p(wavecntx,kcnt,1)=(gr*cos(xi) - gxi*sin(xi))/(2.*pi)
     
          greens_p(wavecntx,kcnt,2)=(gr*sin(xi) + gxi*cos(xi))/(2.*pi)
                                   
	  ELSE !NOT ((p2.eq.0.d0).and.(p1.eq.0.d0))
     
     	  greens_p(wavecntx,kcnt,1) = (psum*(-b_rcv(2,1))*J1+
     .           exp(ci*xi)*(b_rcv(2,2)*(-J2*psum+1./(omega*r)*J1) 
     .           + b_rcv(3,2)*ci/(omega*r)*J1)+
     .           exp(-ci*xi)*(b_rcv(2,3)*(-J0*psum+1./(omega*r)*J1) + 
     .           b_rcv(3,3)*ci/(omega*r)*J1))*omega

          gr = greens_p(wavecntx,kcnt,1)     
     
          greens_p(wavecntx,kcnt,2) =  (psum*b_rcv(3,1)*J1 + 
     .           exp(ci*xi)*(-b_rcv(3,2)*(-J2*psum+1./(omega*r)*J1) + 
     .           b_rcv(2,2)*ci/(omega*r)*J1) +
     .           exp(-ci*xi)*(-b_rcv(3,3)*(-J0*psum+1./(omega*r)*J1)+
     .           b_rcv(2,3)*(ci)/(omega*r)*J1))*omega

          gxi = greens_p(wavecntx,kcnt,2)     
     
          greens_p(wavecntx,kcnt,3) =  psum*(b_rcv(1,1)*J0
     .           + b_rcv(1,2)*J1*exp(ci*xi)-b_rcv(1,3)*J1*exp(-ci*xi))
     .           *omega/(2.*pi)
     
          greens_p(wavecntx,kcnt,1)=(gr*cos(xi) - gxi*sin(xi))/(2.*pi)
     
          greens_p(wavecntx,kcnt,2)=(gr*sin(xi) + gxi*cos(xi))/(2.*pi)
        
          END IF
         END IF
        END DO  !end kcnt loop
	          	           
       END DO   !end p=sqrt(p1^2+p2^2) loop      
	      
750   END SUBROUTINE
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
