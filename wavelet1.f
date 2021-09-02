        program  wavelet1
       
       include 'number.h'
       
       integer i,j,k
       
       double precision wl(ns,ns),wr(ns),wi(ns),df(ns)
       
       complex Ws(ns),ptaper(ns)
       
       common/wfreq/Ws
       common/pfreq/ptaper
       
       open(unit=14,file='wave.dat',form='formatted') 

       call wsvalue
       call wsvalue1 
       
       do i=1,ns
        write(14,*)real(Ws(i)),imag(Ws(i))
       end do	
       
       do i=1,ns
	wr(i)=real(Ws(i))
        wi(i)=imag(WS(i))           
       end do
       
       call fftpot(1,slog,ns,wr,wi)
       
       do i=1,ns
        Ws(i)=cmplx(wr(i),wi(i))
       end do
       
       do i=1,ns
        write(14,*)real(Ws(i)),imag(Ws(i))
       end do	

      print*,Ws(2)
c340     format(2(e12.6,2x))
    
       close(14)
       close(10)
       close(9)
       stop
       end

ccccccccccccccccccccccccccccccccccccccccccccc

        subroutine wsvalue
        
	include 'number.h'

	double precision tt,dx,x,pi,w

	complex Ws(ns)

	integer i,nf1,no,n1,n2,n3

	parameter(x=10.,w=ns)

        common/wfreq/Ws
	
	common/nsample/no,n1,n2,n3
	 	
	pi=4d0*datan(1d0)

	nf1=x*ns/w
	
	print*, nf1,ns,w,x
	
	dx=w/ns
	
	no = (5)*ns/w
	
	n1 = (16)*ns/w+no
	
	
	n2=(32)*ns/w+no

	
	n3=(47)*ns/w+no
	
	do i=1,no
	 Ws(i)=cmplx(0.0,0.0)
	end do

	do i = no+1,n1
	  Ws(i)=cmplx(sin((float(i-no))*dx/x)*
     1	  sin((float(i-no))*dx/x),0.0)
        end do

	do i = n1+1,n2-1
	Ws(i)=cmplx(1.0,0.0)   
        end do

	do i = n2,n3
         Ws(i)=cmplx(cos((float(i-no))*dx/x)*
     1	   cos((float(i-no))*dx/x),0.0)
        end do

	do i = n3+1,ns
           Ws(i)=cmplx(0.0,0.0)
        end do

	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine wsvalue1
        
	include 'number.h'

	double precision tt,dx,x,pi,w

	complex ptaper(ns)

	integer counter,nk_total, i,nf1,no,n1,n2,n3

	parameter(x=10.,w=ns)

        common/pfreq/ptaper
	
	common/nsample/no,n1,n2,n3
	 	
	pi=4d0*datan(1d0)

	nf1=x*ns/w
	
	print*, nf1,ns,w,x
	
	dx=w/ns
	
	no=(5+10)*ns/w
	
	n1=(16)*ns/w+no
	
	n2=(32)*ns/w+no
	
	n3=(47)*ns/w+no
	   
	   nk_total = 100
	   
          factor = 1.0
	        
        DO counter=1,nk_total
         ptaper(counter) = factor
        END DO

      DO counter=1,int(nk_total/16.d0)
        ptaper(counter) = 0.d0
      END DO      

      DO counter=int(15.d0/16.d0*nk_total)+1,nk_total
        ptaper(counter) = 0.d0
      END DO

      DO counter=int(6.d0/8.d0*nk_total)+1,int(15.d0/16.d0*nk_total)
         ptaper(counter) =factor* 
     1   (cos(
     2    dble(counter-6.d0/8.d0*nk_total-1)/
     3    dble(nk_total*15.d0/16.d0-6.d0/8.d0*nk_total-1)*
     4    pi/2.d0))**2
      END DO
      
      DO counter=int(nk_total/16.d0)+1,int(nk_total/8.d0)
        ptaper(counter) =factor* 
     1   (sin(dble(counter-1-nk_total/16.d0)/
     3    dble(nk_total/8.d0-1-nk_total/16.d0)*pi/2.d0))**2
      END DO
      
      END SUBROUTINE


cccccccccccccccccccccccccccccccccccccccccccccccc
    	subroutine fftpot (kdft,mlog,ndim,rre,rim)
c
cccc    implicit double precision (a-h,o-z)
        double precision rre(ndim),rim(ndim)
        double precision cosh,freq,pi,rimh,rreh,sinh,wi,wih,wr,wrh
     0  integer index,k,kdft,mh,mhlk,ml,mlk,mlog,
     1          nd,ndh,ndim,ndimm1,nh,nhdim
c
c
        pi=4d0*datan(1d0)
	nhdim=ndim/2
	ndimm1=ndim-1
	index=1
	
	do 10 nd=1,ndimm1
	  if (nd.ge.index) goto 4
	  rreh=rre(index)
	  rimh=rim(index)
	  rre(index)=rre(nd)
	  rim(index)=rim(nd)
	  rre(nd)=rreh
	  rim(nd)=rimh
4	  k=nhdim
6	  if (k.ge.index) goto 8
	  index=index-k
	  k=k/2
	  goto 6
8	  index=index+k
10	continue
	
	do 40 ml=1,mlog
	  mlk=2**ml
	  mhlk=mlk/2
	  wr=1.
	  wi=0.
	  freq=pi/mhlk
	  cosh=dcos(freq)
	  sinh=dsin(freq)
	  do 30 mh=1,mhlk
	     do 20 nd=mh,ndim,mlk
		ndh=nd+mhlk
		rreh=wr*rre(ndh)-wi*rim(ndh)
		rimh=wi*rre(ndh)+wr*rim(ndh)
		rre(ndh)=rre(nd)-rreh
		rim(ndh)=rim(nd)-rimh
		rre(nd)=rre(nd)+rreh
		rim(nd)=rim(nd)+rimh
20	  continue
	  wrh=wr
	  wih=wi
	  wr=wrh*cosh-wih*sinh
	  wi=wrh*sinh+wih*cosh
30	  continue
40	  continue

	  do 60 nd=1,ndim
	  rre(nd)=rre(nd)/dsqrt(dfloat(ndim))
	  rim(nd)=rim(nd)/dsqrt(dfloat(ndim))
60	  continue
	  if (kdft.eq.2) goto 64
          do 50 nh=2,nhdim
	  index=ndim-nh+2
	  rreh=rre(index)
	  rimh=rim(index)
	  rre(index)=rre(nh)
          rim(index)=rim(nh)
          rre(nh)=rreh
          rim(nh)=rimh
50        continue
64        return
          end
	  
ccccccccccccccccccccccccccccccccccccccccccc
