
      integer ndim,mlog,ns,slog,ypos,xpos,nf
      
      double precision h,r1,r2,c01,c11,t,xs,ys,zs,z,zps,delta,cdm
      
      double precision zps1,drx

      parameter (ndim=256,mlog=8,ns=256,slog=8,r1=2500d0,drx=1d0)
      
      parameter (r2=2500d0,c01=400d0,c11=600d0,t=6.28d0,h=5d0)
     
      parameter (xs=256d0,ys=256d0,zs=100d0,z=0d0,zps=6.25d0)
      
      parameter (ypos=256,xpos=256,cdm=400d0,delta=0.1d0,zps1=6.25d0)

c========================================================================
cNote:
c     mlog : The 2 log of ndim
c     ndim : Number of points of the FFT in the wavenumber domain
c     slog : The 2 log of ns 
c     ns   : Number of points of the FFT in the frequency domain
c     tt    : Length interval of the wavenumber domain
c     w    : Length interval of the frequency domain
c     r1   : The density of the background of the first layer
c     r2   : The density of the background of the second layer
c     c01  : The velocity of the backgroud of the first layer
c     c11  : The velocity of the backgroud of the second layer
c     cdm  : The velocity of the impedance function 
c     sy   : The position of the source in y-direction
c     sx   : The position of the source in x-direction
c     h    : The depth of the layer
c     zs   : The depth of the source
c     z    : The depth of the receiver
c     zps  : The depth of the scatterer
c     ypos : Position of the scatterer in y-direction
c     delta: Pertubation of the velocity
c     dz   : Distance between h and zps

c==========================================================================
c========================================================================== 
