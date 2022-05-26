
c=======================================================================                  
c=======================================================================                  
c     SUBROUTINES AND FUNCTIONS TO EVALUATE THE DENSITY AND CDF OF RVs
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
      subroutine dgamma(y,mu,v,eval)        
c=======================================================================
c     return the log of a gamma distribution
c     A.J.V., 2006
      double precision y,mu,v,eval
      double precision dgamlog
      eval=(v-1)*log(y)-y*(v/mu)-v*log(mu/v)-dgamlog(v)
      return
      end


c=======================================================================
      subroutine dgamma2(y,alpha,beta,eval)        
c=======================================================================
c     return the log of a gamma distribution
c     A.J.V., 2006
      double precision y,alpha,beta,eval
      double precision dgamlog
      eval=alpha*log(beta)+(alpha-1.d0)*log(y) - beta*y - dgamlog(alpha)
      return
      end


c=======================================================================
      subroutine dgammai(y,alpha,beta,eval)        
c=======================================================================
c     return the log of a inverted gamma distribution
c     A.J.V., 2006
      double precision y,alpha,beta,eval
      double precision dgamlog
      eval=alpha*log(beta)-(alpha+1.d0)*log(y) - beta/y - dgamlog(alpha)
      return
      end


c=======================================================================            	  
      double precision function cdfslogistic(x)
c=======================================================================            
c     This function evaluate the cdf of a standard logistic distribution.
c     A.J.V., 2005
      implicit none 
      double precision x
      double precision cdflogis
      cdfslogistic=cdflogis(x,0.d0,1.d0,1,0)
      return
      end
     
c=======================================================================            	  
      double precision function invcdfslogistic(p)
c=======================================================================            
c     This function evaluate the cdf of a standard logistic distribution.
c     A.J.V., 2005
      implicit none 
      double precision p
      double precision invcdflogis
      invcdfslogistic=invcdflogis(p,0.d0,1.d0,1,0)
      return
      end



c=======================================================================            	  
      double precision function cdfbeta(x,alpha,beta)
c=======================================================================            
c     This function evaluate the cdf of a Beta(alpha,beta) distribution.
c     A.J.V., 2005
      implicit none 
      double precision x,alpha,beta
      double precision cdfbetas
      cdfbeta=cdfbetas(x,alpha,beta,1,0)
      return
      end

     
c=======================================================================            	  
      double precision function invcdfbeta(p,alpha,beta)
c=======================================================================            
c     This function evaluate the cdf of a Beta(alpha,beta) distribution.
c     A.J.V., 2005
      implicit none 
      double precision p,alpha,beta
      double precision invcdfbetas
      invcdfbeta=invcdfbetas(p,alpha,beta,1,0)
      return
      end


