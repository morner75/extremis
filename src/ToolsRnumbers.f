
c=======================================================================                  
c=======================================================================                  
c     SUBROUTINES AND FUNCTIONS FOR RN GENERATION
c
c      Alejandro Jara
c      Department of Statistics
c      Facultad de Matematicas
c      Pontificia Universidad Catolica de Chile
c      Casilla 306, Correo 22 
c      Santiago
c      Chile
c      Voice: +56-2-3544506  URL  : http://www.mat.puc.cl/~ajara
c      Fax  : +56-2-3547729  Email: atjara@uc.cl
c
c=======================================================================                  
c=======================================================================                  


c=======================================================================                  
      real function runif()
c=======================================================================                  
c     This function generates a uniform random variable
c     A.J.V., 2006
      real ranf
      runif=ranf()
      return
      end

c=======================================================================                        
      double precision function rexpo(lambda)
c=======================================================================                  
c     This function generates a exp(lambda) random variable  
c     A.J.V., 2005
      implicit none
      double precision lambda
      real runif
      rexpo=-log(1.d0-dble(runif()))/lambda
      return
      end     

c=======================================================================                        
      subroutine rdisc(imin,imax,evali)
c=======================================================================                  
c     This subroutine generates a discrite uniform random variable in 
c     {imin,...,imax}
c     A.J.V., 2006
      implicit none 
      integer imin,imax,evali,ignuin

      evali=ignuin(imin,imax)
      if(evali.lt.imin)evali=imin
      if(evali.gt.imax)evali=imax
      return
      end    

c=======================================================================                  
c      double precision function rgamma(alpha,beta)
c=======================================================================                  
c     This function generates a random gamma value.
c     The parametrization is such that E(X)=alpha/beta
c     A.J.V., 2006 
c      implicit none 
c      double precision beta,alpha
c      real a,r,gengam
c      a=beta
c      r=alpha
c      rgamma = gengam(a,r)
c      return
c      end      


c=======================================================================                  
      double precision function rgamma(alpha,beta)
c=======================================================================                  
c     This function generates a random gamma value.
c     It call gamdv which is a modified version of Mike West's gamma
c     generator.
c     gamdv(alpha)  ~ gamma(alpha,1)
c     rangam(alpha,beta) = gamdv(alpha)/beta ~ gamma(alpha,beta)
c     A.J.V., 2005 
      implicit none 
      double precision beta,alpha,gamdv
      rgamma = gamdv(alpha)
      rgamma = rgamma/beta
      return
      end      

c=======================================================================                  
      double precision function gamdv(a)
c=======================================================================            
c     This subrountine is from Mike West's code. I slightly modified it
c     to use runif as the uniform random number generator.
c     Generates a random gamma variable with shape a>0 and scale=1
c     ix : random seed
c     requires uniform random generator, runif
c     A.J.V., 2005
      implicit none
      double precision a,aa,ea,u0,u1,u2,c1,c2,c3,c4,c5,w
      real runif
      
      aa=a
      if(aa.eq.1.d0)then
         gamdv = -log(dble(runif()))
         return
      endif
      if(aa.lt.1.d0)then
         ea=2.7182818d0/(aa+2.7182818d0)
11       u0=dble(runif())
         u1=dble(runif())
         if(u0.le.ea)then
            gamdv=(u0/ea)**(1.d0/aa)
            if(u1.gt.exp(-gamdv))then
               go to 11
              else
               return
            end if
           else
            gamdv=-log((1.d0-u0)/(ea*aa))
            if(u1.gt.gamdv**(aa-1.d0))then
               go to 11
              else
               return
              end if
            end if
           else
            c1=aa-1.d0
            c2=(aa-1.d0/(6.d0*aa))/c1
            c3=2.d0/c1
            c4=c3+2.d0
            c5=1.d0/sqrt(aa)
12          u1=dble(runif())
            u2=dble(runif())
            u1=u2+c5*(1.d0-1.86d0*u1)
            if(abs(u1-0.5d0).ge.0.5d0)go to 12
            w=c2*u2/u1
            if((c3*log(u1)-log(w)+w).ge.1.d0) then
               go to 12
              else
               gamdv=c1*w
            end if
        end if
        return
        end

c=======================================================================                        
      double precision function rbeta(a0,b0)
c=======================================================================                  
c     This function generates a beta random variable  
c     A.J.V., 2006
      implicit none
      double precision a0,b0
      real genbet,a,b
      a=a0
      b=b0
      rbeta=genbet(a,b)
      return
      end         

c======================================================================      
      double precision function rtbeta(alpha,beta,a,b,ainf,binf)
c=======================================================================            
c     generate truncated(a,b) Beta(alpha,beta)
c     a,b  = end points of interval (ainf = binf = .false.)   
c     ainf = true, if left endpoint is 0; otherwise = false
c     binf = true, if right endpoint is 1; otherwise = false      
c     A.J.V., 2006
      implicit none
      double precision alpha,beta,a,b
      double precision rbeta,invcdfbetas,cdfbetas
      double precision uni,tmp,tmp1,tmp2
      logical ainf,binf
      real runif

      uni=dble(runif())
      rtbeta=0.d0
      if(ainf.and.binf) go to 110
      if(ainf.or.binf) go to 100
      
      if(a.gt.b)then
        call rexit("error in limits rtbeta")
        rtbeta=a
        return
      end if  

      tmp1=cdfbetas(a,alpha,beta,1,0)
      tmp2=cdfbetas(b,alpha,beta,1,0)
      tmp=tmp1+uni*(tmp2-tmp1)
     
      rtbeta=invcdfbetas(tmp,alpha,beta,1,0)
      go to 120

100   if(ainf)then

         tmp2=cdfbetas(b,alpha,beta,1,0)
         tmp=uni*tmp2
         rtbeta=invcdfbetas(tmp,alpha,beta,1,0)
         go to 120
      end if
      if(binf)then
         tmp1=cdfbetas(a,alpha,beta,1,0)      
         tmp=uni+(1.d0-uni)*tmp1
         rtbeta=invcdfbetas(tmp,alpha,beta,1,0)
         go to 120
      end if

110   rtbeta=rbeta(alpha,beta)

120   continue

      return
      end      
      
c======================================================================      
      double precision function rtbeta2(alpha,beta,a,b)
c=======================================================================            
c     generate truncated(a,b) Beta(alpha,beta) using a AR
c     a,b  = end points of interval
c     A.J.V., 2006
      implicit none
      integer exit
      double precision alpha,beta,a,b
      double precision rbeta

      exit=0
      do while(exit.eq.0)
         rtbeta2=rbeta(alpha,beta)
         if(rtbeta2.gt.a.and.rtbeta2.le.b)exit=1
      end do

      return
      end      


c=======================================================================                        
      subroutine dirichlet(alpha,kreal,k,x)
c=======================================================================                  
c     This subroutine generates a dirichlet random vector x  
c     A.J.V., 2005
      implicit none
      integer i,k,kreal
      double precision alpha(kreal),x(kreal),tmp,rgamma,a0

      tmp=0.d0

      do i=1,k
         x(i)=0.d0 
         a0=alpha(i)
         if(a0.gt.0.d0)then
            x(i)=rgamma(a0,1.d0) 
            tmp=tmp+x(i)
         end if   
      end do
      
      do i=1,k
         x(i)=x(i)/tmp
      end do   
      return
      end         

c=======================================================================
      subroutine simdisc(prob,n,m,val)
c=======================================================================
c     generates a sample from a discrete distribution
c     n= real dimension
c     m= used dimension      
c     A.J.V., 2006
      implicit none 
      integer n,m,val,i1,ok
      double precision prob(n),temp1,u,total
      real runif
      
      total=0.d0
      do i1=1,m
         total=total+prob(i1)
      end do

      if(total.eq.0.d0)then
        call rdisc(1,m,val) 
        return
      end if  

c++++ Generating the rn
      temp1=0.d0
      
      u=dble(runif())
      i1=1
      ok=1
      do while(ok.eq.1.and.i1.le.m)
         temp1=temp1+(prob(i1)/total)
         if(u.lt.temp1)then
            val=i1
            ok=0 
         end if
         i1=i1+1
      end do
     
      return
      end           

c=======================================================================
      subroutine simdiscint(prob,n,imin,imax,val)
c=======================================================================
c     generates a sample from a discrete distribution
c     n= real dimension
c     m= used dimension      
c     A.J.V., 2006
      implicit none 
      integer n,val,i1,ok,imin,imax
      double precision prob(n),temp1,u,total
      real runif
      
      total=0.d0
      do i1=imin,imax
         total=total+prob(i1)
      end do

      if(total.eq.0.d0)then
        call rdisc(imin,imax,val) 
        return
      end if  

c++++ Generating the rn
      temp1=0.d0
      
      u=dble(runif())
      i1=imin
      ok=1
      do while(ok.eq.1.and.i1.le.imax)
         temp1=temp1+(prob(i1)/total)
         if(u.lt.temp1)then
            val=i1
            ok=0 
         end if
         i1=i1+1
      end do
      
      return
      end      

      
c=======================================================================            	        
      double precision function rtslogistic(ind,eta)
c=======================================================================            	  
c     This function gerenates a truncated logistic(alpha=0,beta=1) 
c     random variable. The truncation region is (-Inf,eta] if ind=1 and
c     (eta,+Inf) if ind=0.
c     A.J.V., 2005
      implicit none 
      integer ind
      double precision uni,eta,cdfslogistic,invcdfslogistic
      real runif
      uni=runif()
      rtslogistic=0.d0
      if(ind.eq.1)then
         rtslogistic=invcdfslogistic(uni*cdfslogistic(eta))
      end if
      if(ind.eq.0)then
         rtslogistic=invcdfslogistic(uni+(1-uni)*cdfslogistic(eta))
      end if
      return
      end
      
c=======================================================================            	        
      double precision function rtslogistic2(ainf,binf,a,b)
c=======================================================================            	  
c     This function gerenates a truncated logistic(alpha=0,beta=1) 
c     random variable. The truncation region is (a,b) 
c     a,b  = end points of interval (ainf = binf = .false.)   
c     ainf = true, if left endpoint is infinite; otherwise = false
c     binf = true, if right endpoint is infinite; otherwise = false      
c     A.J.V., 2005
      implicit none 
      double precision uni,cdfslogistic,invcdfslogistic,a,b
      real runif
      logical ainf,binf
      
      uni=dble(runif())
      rtslogistic2=0.d0
      if(ainf.and.binf) go to 110
      if(ainf.or.binf) go to 100
      
      if(a.gt.b)then
        call rexit("error in limits rtslogistic2")
        rtslogistic2=a
        return
      end if  

      rtslogistic2=invcdfslogistic(cdfslogistic(a)+
     &             uni*(cdfslogistic(b)-cdfslogistic(a)))
      go to 120

100   if(ainf)then
         rtslogistic2=invcdfslogistic(uni*cdfslogistic(b))
         go to 120
      end if
      if(binf)then
         rtslogistic2=invcdfslogistic(uni+(1.d0-uni)*cdfslogistic(a))
         go to 120
      end if

110   rtslogistic2=invcdfslogistic(uni)

120   continue
      return
      end      
      

c=======================================================================            
      double precision function rnorm(mu,sd)
c=======================================================================            
c     This function generates a N(mu,sd^2) random values.
c     A.J.V., 2006
      implicit none
      double precision mu,sd
      real gennor,av0,sd0
      
      av0=mu
      sd0=sd
      rnorm = gennor(av0,sd0)
      return
      end

c======================================================================      
      double precision function rtnorm(mu,sd,a,b,ainf,binf)
c=======================================================================            
c     generate truncated(a,b) Normal(mu,sd**2) using the Geweke's 
c     algorithm.
c     mu is the mean of TN distribution
c     sd is standard deviation of TN distribution
c     a,b  = end points of interval (ainf = binf = .false.)   
c     ainf = true, if left endpoint is infinite; otherwise = false
c     binf = true, if right endpoint is infinite; otherwise = false      
c     A.J.V., 2006
      implicit none
      double precision mu,sd,a,b,a1,b1,rtsnorm
      logical ainf,binf
      a1=(a-mu)/sd
      b1=(b-mu)/sd
      rtnorm=mu+sd*rtsnorm(a1,b1,ainf,binf)
      return
      end      
      
c======================================================================            
      double precision function rtsnorm(a,b,la,lb)
c======================================================================            
c     generates a N(0,1) random variable
c     subject to the constraint that it be in an interval
c     (a,b), where the endpoints may be finite or infinite.
c     a, b    endpoints of interval; a < b if la = lb = .false.
c     la      .true. if left endpoint is - infinity; in this
c             case A is ignored.
c     lb      .true. if right endpoint is + infinity; in this
c             case B is ignored.
c     A.J.V., 2006
      implicit double precision (a-h,o-z)
      logical la,lb,lflip
      real runif
      double precision dexpone,rnorm
      data eps,t1,t2,t3,t4/2.0d0,.375d0,2.18d0,.725d0,.45d0/

      if(la.and.lb)go to 160
      lflip=.false.
      if(la.or.lb)go to 100
      if(b.le.a)go to 170
c ******* Finite interval
      c1=a
      c2=b
      if((c1*c2).gt.0.0d0)go to 30
c ++++ (A,B) includes 0
      if((c1.gt.-t1).and.(c2.lt.t1))go to 20
c -- F(A) or F(B) small: full normal with rejection
   10 x=rnorm(0.d0,1.d0)
      if(x.lt.c1)go to 10
      if(x.gt.c2)go to 10
      GO TO 150
c -- F(A) and F(B) large: uniform importance sampling
   20 cdel=c2-c1
   25 x=c1+cdel*dble(runif())
      if(dble(runif()).gt.dexpone(x))go to 25
      go to 150
c ++++ (A,B) excludes 0
c -- Transform to both positive
   30 if(c1.gt.0.0d0)go to 40
      c=c1
      c1=-c2
      c2=-c
      lflip=.true.
   40 f1=dexpone(c1)
      f2=dexpone(c2)
      if(f2.lt.eps)go to 60
      if((f1/f2).gt.t2)go to 60
c  -- F(A)/F(B) not large: uniform importance sampling
      cdel=c2-c1
   55 x=c1+cdel*runif()
      if(dble(runif()).gt.(dexpone(x)/f1))go to 55
      go to 140
   60 if(c1.gt.t3)go to 80
c -- P(X>A) and F(A)/F(B) large: half-normal with rejection
   70 x=abs(rnorm(0.d0,1.d0))
      if(x.lt.c1)go to 70
      if(x.gt.c2)go to 70
      go to 140
c -- P(X>A) small, F(A)/F(B) large: exponential importance
c    sampling with rejection
   80 c=c2-c1
   90 z=-log(runif())/c1
      if(z.gt.c)go to 90
      if(dble(runif()).gt.dexpone(z))GO TO 90
      x=c1+z
      go to 140
c ****** Half-line interval
  100 c1=a
c -- Transform to bound from below if A = -infinity
      if(lb)go to 110
      c1=-b
      lflip=.true.
  110 if(c1.gt.t4)go to 130
c -- A not large: full normal with rejection
  120 x=rnorm(0.d0,1.d0)
      if(x.lt.c1)go to 120
      go to 140
c -- A small: exponential importance sampling
  130 z=-log(runif())/c1
      if(dble(runif()).gt.dexpone(z))go to 130
      x=c1+z
  140 if(lflip)x=-x
  150 rtsnorm=X
      return
c ****** Whole interval
  160 rtsnorm=rnorm(0.d0,1.d0)
      return
c  ***** Singleton
  170 rtsnorm=A
      return
      end

c=======================================================================                        
      double precision function dexpone(x)
c=======================================================================                       
c     evaluate a exponential function
c     A.J.V., 2006
      implicit none
      double precision x,expin
      expin=-.5d0*x**2
      if (expin .le. -50.0d0) then
        dexpone=0.0d0
      else
        dexpone=dexp(expin)
      end if
      return
      end
      
    
c=======================================================================      
      subroutine normalvec(n,work)
c=======================================================================
c     generates a vector of normal variables
c     A.J.V., 2006
      implicit none
      integer i,n
      double precision work(n),rnorm
      
      do i=1,n
         work(i)=rnorm(0.d0,1.d0)
      end do
      return
      end


c=======================================================================                        
      subroutine rbinom(n,p,evali)
c=======================================================================                  
c     This subroutine generates a Binomial(n,p) random variable 
c     A.J.V., 2006
      implicit none 
      integer n,evali,ignbin
      double precision p
      real pp
      
      pp=p

      evali=ignbin(n,pp) 

      return
      end        

c======================================================================      
      double precision function rtlnorm(mu,sd,a,b,ainf,binf)
c=======================================================================            
c     generate truncated(a,b) LogNormal(mu,sd**2) using the Geweke's 
c     algorithm.
c     mu is the mean of log(variable)
c     sd is standard deviation of log(variable)
c     a,b  = end points of interval (ainf = binf = .false.)   
c     ainf = true, if left endpoint is infinite; otherwise = false
c     binf = true, if right endpoint is infinite; otherwise = false      
c     A.J.V., 2006
      implicit none
      double precision mu,sd,a,b,a1,b1,rtsnorm,x
      logical ainf,binf
      
      if(ainf)then
         a1=0.d0
        else
         a1=(log(a)-mu)/sd  
      end if  

      if(binf)then
         b1=0.d0
        else
         b1=(log(b)-mu)/sd  
      end if  

      x=mu+sd*rtsnorm(a1,b1,ainf,binf)
      rtlnorm=exp(x)
      
      return
      end   
      
c=======================================================================                  
      double precision function rchisq(nu)
c=======================================================================                  
c     This function generates a random chi2 value.
c     A.J.V., 2006
      implicit none 
      double precision nu
      double precision rgamma
      rchisq=rgamma(0.5d0*nu,0.5d0)
      return
      end            
      

c=======================================================================                  
      subroutine samalph(alpha,aa0,ab0,ndis,k)
c=======================================================================            
c     This routine samples another alpha value using the technique
c     in Escobar and West (TR 533). The prior distribution for alpha
c     is a gamma(aa0, ab0).
c
c     ndis: the number of clusters.
c     k   : Sample size
c     A.J.V., 2006
      integer ndis,k
      double precision alpha, aa0, ab0, xalp,s,e,rgamma
      real runif

      xalp=rgamma(1.d0+alpha,1.d0)
      xalp=xalp/(xalp+rgamma(dble(k),1.d0))
      s=ab0-dlog(xalp)
      e=aa0+dble(ndis)-1.d0
      if(dble(runif()).lt.e/(e+dble(k)*s))then
         e=e+1.d0
      end if
      alpha=rgamma(e,s)
      return
      end


c=======================================================================      
      integer function rpois(mu)
c=======================================================================      
c     This function gerenates a Poisson(mu) random variable. 
c     A.J.V., 2006
      integer ignpoi
      double precision  mu
      real mean

      mean=mu
      rpois=ignpoi(mean)
      return
      end
c======================================================================      
      subroutine rperm(nl,n,p)
c=======================================================================            
c     generate a random permutation of the first n integers
c     A.J.V., 2007
      implicit none
      integer n,nl,p(nl),i,j,k,ipj,itemp,m
      double precision u(100)
      real runif

      do i=1,n
         p(i)=i
      end do
      
      do i=1,n,100
         m=min(n-i+1,100)
         do j=1,100
            u(j)=dble(runif())
         end do
         
         do j=1,m
            ipj=i+j-1
            k=int(u(j)*(n-ipj+1))+ipj
            itemp=p(ipj)
            p(ipj)=p(k)
            p(k)=itemp
         end do
      end do
      return
      end

c=======================================================================                        
      integer function rpoiss(mu)
c=======================================================================                  
c     This function generates a poisson random variable  
c     A.J.V., 2006
      implicit none
      double precision mu
      integer ignpoi
      real mur
      
      mur=mu
      rpoiss=ignpoi(mur)
      return
      end         



c======================================================================      
      integer function rpoiss2(mu)
c======================================================================            
c     generate a Poisson(mu) random number
c     A.J.V., 2007
c======================================================================
      implicit none
      integer i
      double precision mu,uni,p,f
      real runif

      i=0
      p=dexp(-mu)
      f=p
      uni=dble(runif())

      rpoiss2=i
      do while(uni.ge.f)
         p=mu*p/dble(i+1)
         f=f+p
         i=i+1
         rpoiss2=i
      end do

      return
      end      


c======================================================================      
      subroutine rtbetas(alpha,beta,a,b,ainf,binf,val)
c=======================================================================            
c     generate truncated(a,b) Beta(alpha,beta)
c     a,b  = end points of interval (ainf = binf = .false.)   
c     ainf = true, if left endpoint is 0; otherwise = false
c     binf = true, if right endpoint is 1; otherwise = false      
c     A.J.V., 2006
      implicit none
      double precision alpha,beta,a,b
      double precision bound      
      double precision rbeta
      double precision val
      double precision uni,tmp,tmp1,tmp2,tmp3,tmp4
      logical ainf,binf
      real runif
      integer status

      uni=dble(runif())
      val=0.d0
      if(ainf.and.binf) go to 110
      if(ainf.or.binf) go to 100
      
      if(a.gt.b)then
        call rexit("error in limits rtbetas")
        val=a
        return
      end if  


      call cdfbet(1,tmp1,tmp3,a,1.d0-a,alpha,beta,status,bound)
      call cdfbet(1,tmp2,tmp3,b,1.d0-b,alpha,beta,status,bound)
      tmp=tmp1+uni*(tmp2-tmp1)
      
      tmp2=1.d0-tmp
      call cdfbet(2,tmp,tmp2,val,tmp4,alpha,beta,status,bound)

      go to 120

100   if(ainf)then

         call cdfbet(1,tmp2,tmp3,b,1.d0-b,alpha,beta,status,bound)
         tmp=uni*tmp2
         tmp2=1.d0-tmp
         call cdfbet(2,tmp,tmp2,val,tmp4,alpha,beta,status,bound)
         go to 120
      end if
      if(binf)then

         call cdfbet(1,tmp1,tmp3,a,1.d0-a,alpha,beta,status,bound)
         tmp=uni+(1.d0-uni)*tmp1
      
         tmp2=1.d0-tmp
         call cdfbet(2,tmp,tmp2,val,tmp4,alpha,beta,status,bound)
         go to 120
      end if

110   val=rbeta(alpha,beta)

120   continue

      return
      end   







