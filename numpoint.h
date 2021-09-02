      integer ndim,mlog,ns,slog,dim,maxfreq
      
      double precision Lx,Ly,csr,cpr,zrcv,zsrc,xsrc,ysrc,rho,delta
      
      double precision dlx,dly
            
      parameter (ndim=256,mlog=8,ns=256,slog=8,Lx=6.28d0,rho=1600d0)
      
      parameter (xsrc=0d0,ysrc=0d0,delta=0.1d0,dlx=Lx/ndim,dim=3)
      
      parameter (zsrc=0d0,maxfreq=47)
      
      parameter (Ly=6.28d0,csr=400d0,cpr=600d0,zrcv=5d0,dly=Ly/ndim)

