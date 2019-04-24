c	**********************************************	  
	double precision function sgn(x)
	double precision x

	if(x.lt.0.D0) then
	  sgn=-1.D0
	else
	  sgn=1.D0
	end if

	return
	end

c	**********************************************	  
	double precision function theta(x)
	double precision x

	if(x.lt.0.D0) then
	  theta=0.D0
	else if(x.eq.0.d0) then
	  theta=0.5d0
	else
	  theta=1.d0
	end if

	return
	end

c	**********************************************	  



	double precision function fermic(x,beta)
	double precision x,beta
c	fermi function which avoids underflow/overflow.  If beta=0,
c	it is assumed that T=0 is meant!

	if(beta.eq.0.0D0) then
	  if(x.lt.0.0D0) then	  
	    fermic=1.0D0
	  else if(x.eq.0.0D0) then
	    fermic=0.5D0
	  else
	    fermic=0.0D0
	  end if
	else
	  if(x.lt.0.0D0) then	  
	    fermic=1.0D0/(dexp(beta*x)+1.0D0)
	  else
	    fermic=dexp(-beta*x)/(dexp(-beta*x)+1.0D0)
	  end if
	end if
	return
	end

c	**********************************************	  
	double precision function dfermic(x,beta)
	double precision x,beta
c
c	This function returns the minus derivative of the fermi function
c
c                               1
c	fermic(x,beta)= ----------------
c                       exp(-beta x) + 1
c
c       d fermic(x,beta)     -beta exp(-beta x)
c	----------------  = --------------------- = -dfermic(x,beta)
c              d x           (exp(-beta x) + 1)^2
c
c
	if(x.lt.0.0D0) then	  
	  dfermic=beta*dexp(beta*x)/(dexp(beta*x)+1.0D0)**2
        else
	  dfermic=beta*dexp(-beta*x)/(dexp(-beta*x)+1.0D0)**2
	end if

	return
	end

c	**********************************************	  
c
	INTEGER FUNCTION LOCATE(N,xx,ll,ul,x)
	integer N,ll,ul
	real*8 xx(1:N),x
c
c	given an array xx(1:n),and given a value x, returns a value j 
c 	such that x is between xx(j) and xx(j+1). xx(1:n) must be
c 	monotonic.

	integer jl,jm,ju


	if(x.lt.xx(ll)) then
	  locate=ll
c	  pause'input out of range in locate,left'
	  return
	end if

	if(x.gt.xx(ul)) then
	  locate=ul-1
c	  pause'input out of range in locate,right'
	  return
	end if

	jl=ll
	ju=ul
        do while((ju-jl) .gt. 1) 
	  jm=(ju+jl)/2
	  if (x .ge. xx(jm)) then
	    jl=jm
	  else
	    ju=jm
	 endif
	end do 
	locate=jl
c
	return
	end

c	****************************************************************
	double complex function gauss(z)
c	This block calculates -i*sqrt(pi)*w(z)
	double complex z,ii
	double precision pi2,x,y,t
	logical flag
	parameter(pi2=1.7724539D0)
	ii=dcmplx(0.D0,1.D0)

	t=1.d0
	if(dimag(z).lt.0.0D0) then
	  call wofz(dreal(z),-dimag(z),x,y,flag)
	  gauss=-pi2*ii*dcmplx(x,-y)
	else
	  call wofz(dreal(z),dimag(z),x,y,flag)
	  gauss=-pi2*ii*dcmplx(x,y)
	end if

	if(flag) write(6,*) 'error in cerfjd'
	return
	end


c	**************************************************
        complex*16 function sem(z)
c       This block calculates Hilbert transform for a semicircular 
c       DOS with t^*=t
c       In general sem=8*(z-sqrt(z**2-D**2/4))/D**2
        double complex z,z1,z2,ii,sem1,sem2,sz2
        double precision dband,pi,t,t2

        include 'Glob_cons'

	t=1.d0
!        t2=2.d0*t
	dband=t
        z1=dcmplx(dband,0.D0)
					             
        z2=z**2-z1**2
				                  
        sz2=zsqrt(z2)
												               
        sem1=2.D0/(z+sz2)
												                    
        sem2=2.D0/(z-sz2)
        
   
        if(dimag(sem1).le.0.D0) then
             sem=sem1
        else if(dimag(sem2).le.0.D0) then
             sem=sem2
        else 
             write(6,*) 'no causal root found'
             stop 
        end if
    
        return
        end  


c	*************************************************************
	
	  double complex function flat(z)
	  double complex z,ii
	  double precision rews,imws,rez,imz,pi
	  double precision dband,epsvh,zr,zi,acc,dnorm,rho0

	  common/dos/dband,epsvh,zr,zi,acc,dnorm,nnorm,idos

	
c	  Flat Band from -D/2 to D/2 where 1.D0
	  
	  include 'Glob_cons'

	  rez=dreal(z)
	  imz=dimag(z)
	  
	  if(dabs(rez).eq.dband/2.D0) rez=rez+1.D-6
	  if(dabs(imz).le.1.D-6) then
	    rews=dlog(dabs((rez+dband/2.D0)/(rez-dband/2.D0)))
	    if(dabs(rez).ge.dband/2.D0) then
	      imws=0.D0
	    else
	    imws=-pi
	    end if
	  else
	    rews=0.5D0*dlog(((rez+dband/2.D0)**2+imz**2)/
     &	         	   ((rez-dband/2.D0)**2+imz**2))

	    imws=-(datan((dband/2.D0-rez)/imz)+
     &		  datan((dband/2.D0+rez)/imz))
	  end if

 355	  flat=dcmplx(rews,imws)

 	  return
	  end

c	*************************************************************

C      ALGORITHM 680, COLLECTED ALGORITHMS FROM ACM.
C      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
C      VOL. 16, NO. 1, PP. 47.
      SUBROUTINE WOFZ (XI, YI, U, V, FLAG)
C
C  GIVEN A COMPLEX NUMBER Z = (XI,YI), THIS SUBROUTINE COMPUTES
C  THE VALUE OF THE FADDEEVA-FUNCTION W(Z) = EXP(-Z**2)*ERFC(-I*Z),
C  WHERE ERFC IS THE COMPLEX COMPLEMENTARY ERROR-FUNCTION AND I
C  MEANS SQRT(-1).
C  THE ACCURACY OF THE ALGORITHM FOR Z IN THE 1ST AND 2ND QUADRANT
C  IS 14 SIGNIFICANT DIGITS; IN THE 3RD AND 4TH IT IS 13 SIGNIFICANT
C  DIGITS OUTSIDE A CIRCULAR REGION WITH RADIUS 0.126 AROUND A ZERO
C  OF THE FUNCTION.
C  ALL REAL VARIABLES IN THE PROGRAM ARE DOUBLE PRECISION.
C
C
C  THE CODE CONTAINS A FEW COMPILER-DEPENDENT PARAMETERS :
C     RMAXREAL = THE MAXIMUM VALUE OF RMAXREAL EQUALS THE ROOT OF
C                RMAX = THE LARGEST NUMBER WHICH CAN STILL BE
C                IMPLEMENTED ON THE COMPUTER IN DOUBLE PRECISION
C                FLOATING-POINT ARITHMETIC
C     RMAXEXP  = LN(RMAX) - LN(2)
C     RMAXGONI = THE LARGEST POSSIBLE ARGUMENT OF A DOUBLE PRECISION
C                GONIOMETRIC FUNCTION (DCOS, DSIN, ...)
C  THE REASON WHY THESE PARAMETERS ARE NEEDED AS THEY ARE DEFINED WILL
C  BE EXPLAINED IN THE CODE BY MEANS OF COMMENTS
C
C
C  PARAMETER LIST
C     XI     = REAL      PART OF Z
C     YI     = IMAGINARY PART OF Z
C     U      = REAL      PART OF W(Z)
C     V      = IMAGINARY PART OF W(Z)
C     FLAG   = AN ERROR FLAG INDICATING WHETHER OVERFLOW WILL
C              OCCUR OR NOT; TYPE LOGICAL;
C              THE VALUES OF THIS VARIABLE HAVE THE FOLLOWING
C              MEANING :
C              FLAG=.FALSE. : NO ERROR CONDITION
C              FLAG=.TRUE.  : OVERFLOW WILL OCCUR, THE ROUTINE
C                             BECOMES INACTIVE
C  XI, YI      ARE THE INPUT-PARAMETERS
C  U, V, FLAG  ARE THE OUTPUT-PARAMETERS
C
C  FURTHERMORE THE PARAMETER FACTOR EQUALS 2/SQRT(PI)
C
C  THE ROUTINE IS NOT UNDERFLOW-PROTECTED BUT ANY VARIABLE CAN BE
C  PUT TO 0 UPON UNDERFLOW;
C
C  REFERENCE - GPM POPPE, CMJ WIJERS; MORE EFFICIENT COMPUTATION OF
C  THE COMPLEX ERROR-FUNCTION, ACM TRANS. MATH. SOFTWARE.
C
*
*
*
*
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
*
      LOGICAL A, B, FLAG
      PARAMETER (FACTOR   = 1.12837916709551257388D0,
     *           RMAXREAL = 0.5D+154,
     *           RMAXEXP  = 708.503061461606D0,
     *           RMAXGONI = 3.53711887601422D+15)
*
      FLAG = .FALSE.
*
      XABS = DABS(XI)
      YABS = DABS(YI)
      X    = XABS/6.3
      Y    = YABS/4.4
*
C
C     THE FOLLOWING IF-STATEMENT PROTECTS
C     QRHO = (X**2 + Y**2) AGAINST OVERFLOW
C
      IF ((XABS.GT.RMAXREAL).OR.(YABS.GT.RMAXREAL)) GOTO 100
*
      QRHO = X**2 + Y**2
*
      XABSQ = XABS**2
      XQUAD = XABSQ - YABS**2
      YQUAD = 2*XABS*YABS
*
      A     = QRHO.LT.0.085264D0
*
      IF (A) THEN
C
C  IF (QRHO.LT.0.085264D0) THEN THE FADDEEVA-FUNCTION IS EVALUATED
C  USING A POWER-SERIES (ABRAMOWITZ/STEGUN, EQUATION (7.1.5), P.297)
C  N IS THE MINIMUM NUMBER OF TERMS NEEDED TO OBTAIN THE REQUIRED
C  ACCURACY
C
        QRHO  = (1-0.85*Y)*DSQRT(QRHO)
        N     = IDNINT(6 + 72*QRHO)
        J     = 2*N+1
        XSUM  = 1.0/J
        YSUM  = 0.0D0
        DO 10 I=N, 1, -1
          J    = J - 2
          XAUX = (XSUM*XQUAD - YSUM*YQUAD)/I
          YSUM = (XSUM*YQUAD + YSUM*XQUAD)/I
          XSUM = XAUX + 1.0/J
 10     CONTINUE
        U1   = -FACTOR*(XSUM*YABS + YSUM*XABS) + 1.0
        V1   =  FACTOR*(XSUM*XABS - YSUM*YABS)
        DAUX =  DEXP(-XQUAD)
        U2   =  DAUX*DCOS(YQUAD)
        V2   = -DAUX*DSIN(YQUAD)
*
        U    = U1*U2 - V1*V2
        V    = U1*V2 + V1*U2
*
      ELSE
C
C  IF (QRHO.GT.1.O) THEN W(Z) IS EVALUATED USING THE LAPLACE
C  CONTINUED FRACTION
C  NU IS THE MINIMUM NUMBER OF TERMS NEEDED TO OBTAIN THE REQUIRED
C  ACCURACY
C
C  IF ((QRHO.GT.0.085264D0).AND.(QRHO.LT.1.0)) THEN W(Z) IS EVALUATED
C  BY A TRUNCATED TAYLOR EXPANSION, WHERE THE LAPLACE CONTINUED FRACTION
C  IS USED TO CALCULATE THE DERIVATIVES OF W(Z)
C  KAPN IS THE MINIMUM NUMBER OF TERMS IN THE TAYLOR EXPANSION NEEDED
C  TO OBTAIN THE REQUIRED ACCURACY
C  NU IS THE MINIMUM NUMBER OF TERMS OF THE CONTINUED FRACTION NEEDED
C  TO CALCULATE THE DERIVATIVES WITH THE REQUIRED ACCURACY
C
*
        IF (QRHO.GT.1.0) THEN
          H    = 0.0D0
          KAPN = 0
          QRHO = DSQRT(QRHO)
          NU   = IDINT(3 + (1442/(26*QRHO+77)))
        ELSE
          QRHO = (1-Y)*DSQRT(1-QRHO)
          H    = 1.88*QRHO
          H2   = 2*H
          KAPN = IDNINT(7  + 34*QRHO)
          NU   = IDNINT(16 + 26*QRHO)
        ENDIF
*
        B = (H.GT.0.0)
*
        IF (B) QLAMBDA = H2**KAPN
*
        RX = 0.0
        RY = 0.0
        SX = 0.0
        SY = 0.0
*
        DO 11 N=NU, 0, -1
          NP1 = N + 1
          TX  = YABS + H + NP1*RX
          TY  = XABS - NP1*RY
          C   = 0.5/(TX**2 + TY**2)
          RX  = C*TX
          RY  = C*TY
          IF ((B).AND.(N.LE.KAPN)) THEN
            TX = QLAMBDA + SX
            SX = RX*TX - RY*SY
            SY = RY*TX + RX*SY
            QLAMBDA = QLAMBDA/H2
          ENDIF
 11     CONTINUE
*
        IF (H.EQ.0.0) THEN
          U = FACTOR*RX
          V = FACTOR*RY
        ELSE
          U = FACTOR*SX
          V = FACTOR*SY
        END IF
*
        IF (YABS.EQ.0.0) U = DEXP(-XABS**2)
*
      END IF
*
*
C
C  EVALUATION OF W(Z) IN THE OTHER QUADRANTS
C
*
      IF (YI.LT.0.0) THEN
*
        IF (A) THEN
          U2    = 2*U2
          V2    = 2*V2
        ELSE
          XQUAD =  -XQUAD
*
C
C         THE FOLLOWING IF-STATEMENT PROTECTS 2*EXP(-Z**2)
C         AGAINST OVERFLOW
C
          IF ((YQUAD.GT.RMAXGONI).OR.
     *        (XQUAD.GT.RMAXEXP)) GOTO 100
*
          W1 =  2*DEXP(XQUAD)
          U2  =  W1*DCOS(YQUAD)
          V2  = -W1*DSIN(YQUAD)
        END IF
*
        U = U2 - U
        V = V2 - V
        IF (XI.GT.0.0) V = -V
      ELSE
        IF (XI.LT.0.0) V = -V
      END IF
*
      RETURN
*
  100 FLAG = .TRUE.
      RETURN
*
      END

c	************************************************************

      SUBROUTINE LINT(XA,YA,Y2A,N,KLO,X,Y) 
      INTEGER N,KLO,KHI,K
      DOUBLE PRECISION XA(N),YA(N),Y2A(N),X,Y,H,B
      
      IF(KLO.EQ.0) THEN		! MEANS UNSET IN MAIN PROGRAM
      KLO=1 
      END IF
      KHI=N 
1     IF (KHI-KLO.GT.1) THEN 
        K=(KHI+KLO)/2 
        IF(XA(K).GT.X)THEN 
          KHI=K 
        ELSE 
          KLO=K 
        ENDIF 
      GOTO 1 
      ENDIF 
      H=XA(KHI)-XA(KLO) 
      IF (H.EQ.0.D0) then
         write(6,*) 'Bad XA input.' 
         stop
      END IF
      B=(X-XA(KLO))
      Y=YA(KLO)+Y2A(KLO)*B
      RETURN 
      END 

c	**********************************************	  
      SUBROUTINE LINE(X,Y,N,Y2) 
      PARAMETER (NMAX=5000) 
      INTEGER N,I
      DOUBLE PRECISION X(N),Y(N),Y2(N)

      DO I=1,N-1
        Y2(I)=(Y(I+1)-Y(I))/(X(I+1)-X(I))
      END DO
      Y2(N)=Y2(N-1)
        
      RETURN 
      END 

c	**********************************************	  

!  Inversion of a complex matrix a(m,m) 
       subroutine cmatinv(a,ai,m) 
       implicit none
       complex*16:: a(m,m),ai(m,m),zero 
       integer m,i,j,k 
       zero=dcmplx(0.d0,0.d0)
       ai=a
       do k=1,m 
        do j=1,m
        if(j.ne.k) then
            if (ai(k,k).eq.zero) ai(k,k)=dcmplx(1.d-9,1.D-10)
            ai(k,j)=ai(k,j)/ai(k,k) 
        endif
       end do 
       ai(k,k)=1.0d0/ai(k,k) 
       do i = 1,m 
        if (i.ne.k) then 
            do j=1,m 
                if (j.ne.k) ai(i,j)=ai(i,j)-ai(k,j)*ai(i,k)
             end do
         end if
        end do 
        do i=1,m 
        if (i.ne.k) ai(i,k)=-ai(i,k)*ai(k,k) 
       enddo
       end do 
       end subroutine
      SUBROUTINE broydn(x,n,check,funcv)
C     Given an initial guess x(1:n) for a root in n dimensions, find
C     root. The vector of functions to be zeroed, called fvec(1:n)
C     in the routine below, is returned by a user-supplied subroutine
C     that must be called funcv and have the declaration subroutine
C     funcv(n,x,fvec).
      INTEGER n,nn,NP,MAXITS
      DOUBLE PRECISION x(n),fvec,EPS,TOLF,TOLMIN,TOLX,STPMX
      LOGICAL check
      PARAMETER (NP=40,MAXITS=50,EPS=1.d-7,TOLF=1.d-7,
     *TOLMIN=1.d-8,TOLX=EPS,STPMX=100.D0)
      COMMON /newtv/ fvec(NP),nn
CU    USES fdjac,fmin,lnsrch,qrdcmp,qrupdt,rsolv
      INTEGER i,its,j,k
      DOUBLE PRECISION den,f,fold,stpmax,sum,temp,test,c(NP),
     *d(NP),fvcold(NP),g(NP),
     *p(NP),qt(NP,NP),r(NP,NP),s(NP),t(NP),w(NP),xold(NP),
     *fmin
      LOGICAL restrt,sing,skip
      EXTERNAL fmin
      EXTERNAL funcv
      nn=n
      f=fmin(x,funcv)
      test=0.D0
      do 11 i=1,n
        if(dabs(fvec(i)).gt.test)test=dabs(fvec(i))
11    continue
      if(test.lt..01D0*TOLF)return
      sum=0.D0
      do 12 i=1,n
        sum=sum+x(i)**2
12    continue
      stpmax=STPMX*max(dsqrt(sum),dfloat(n))
      restrt=.true.
      do 44 its=1,MAXITS
c
c	write no of iterations
c
c	write(6,*)'its=',its
        if(restrt)then
          call fdjac(n,x,fvec,NP,r,funcv)
          call qrdcmp(r,n,NP,c,d,sing)
          if(sing) write(6,*) 'singular Jacobian in broydn'
          do 14 i=1,n
            do 13 j=1,n
              qt(i,j)=0.D0
13          continue
            qt(i,i)=1.D0
14        continue
          do 18 k=1,n-1
            if(c(k).ne.0.D0)then
              do 17 j=1,n
                sum=0.D0
                do 15 i=k,n
                  sum=sum+r(i,k)*qt(i,j)
15              continue
                sum=sum/c(k)
                do 16 i=k,n
                  qt(i,j)=qt(i,j)-sum*r(i,k)
16              continue
17            continue
            endif
18        continue
          do 21 i=1,n
            r(i,i)=d(i)
            do 19 j=1,i-1
              r(i,j)=0.D0
19          continue
21        continue
        else
          do 22 i=1,n
            s(i)=x(i)-xold(i)
22        continue
          do 24 i=1,n
            sum=0.D0
            do 23 j=i,n
              sum=sum+r(i,j)*s(j)
23          continue
            t(i)=sum
24        continue
          skip=.true.
          do 26 i=1,n
            sum=0.D0
            do 25 j=1,n
              sum=sum+qt(j,i)*t(j)
25          continue
            w(i)=fvec(i)-fvcold(i)-sum
        if(dabs(w(i)).ge.EPS*(dabs(fvec(i))+dabs(fvcold(i))))then
              skip=.false.
            else
              w(i)=0.D0
            endif
26        continue
          if(.not.skip)then
            do 28 i=1,n
              sum=0.D0
              do 27 j=1,n
                sum=sum+qt(i,j)*w(j)
27            continue
              t(i)=sum
28          continue
            den=0.D0
            do 29 i=1,n
              den=den+s(i)**2
29          continue
            do 31 i=1,n
              s(i)=s(i)/den
31          continue
            call qrupdt(r,qt,n,NP,t,s)
            do 32 i=1,n
              if(r(i,i).eq.0.D0) write(6,*) 'r singular in broydn'
              d(i)=r(i,i)
32          continue
          endif
        endif
        do 34 i=1,n
          sum=0.D0
          do 33 j=1,n
            sum=sum+qt(i,j)*fvec(j)
33        continue
          g(i)=sum
34      continue
        do 36 i=n,1,-1
          sum=0.D0
          do 35 j=1,i
            sum=sum+r(j,i)*g(j)
35        continue
          g(i)=sum
36      continue
        do 37 i=1,n
          xold(i)=x(i)
          fvcold(i)=fvec(i)
37      continue
        fold=f
        do 39 i=1,n
          sum=0.D0
          do 38 j=1,n
            sum=sum+qt(i,j)*fvec(j)
38        continue
          p(i)=-sum
39      continue
        call rsolv(r,n,NP,d,p)
        call lnsrch(n,xold,fold,g,p,x,f,stpmax,check,fmin,funcv)
        test=0.D0
        do 41 i=1,n
          if(dabs(fvec(i)).gt.test)test=dabs(fvec(i))
41      continue
        if(test.lt.TOLF)then
          check=.false.
          return
        endif
        if(check)then
          if(restrt)then
            return
          else
            test=0.D0
            den=max(f,.5D0*n)
            do 42 i=1,n
              temp=dabs(g(i))*max(dabs(x(i)),1.D0)/den
              if(temp.gt.test)test=temp
42          continue
            if(test.lt.TOLMIN)then
              return
            else
              restrt=.true.
            endif
          endif
        else
          restrt=.false.
          test=0.D0
          do 43 i=1,n
            temp=(dabs(x(i)-xold(i)))/max(dabs(x(i)),1.D0)
            if(temp.gt.test)test=temp
43        continue
          if(test.lt.TOLX)return
        endif
44    continue
	write(6,*) 'MAXITS exceeded in broydn'
c      pause 'MAXITS exceeded in broydn'
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 3ki-i]3|-9.
c
c
c
c
      SUBROUTINE fdjac(n,x,fvec,np,df,funcv)
      INTEGER n,np,NMAX
      DOUBLE PRECISION df(np,np),fvec(n),x(n),EPS
      PARAMETER (NMAX=40,EPS=1.d-4)
CU    USES funcv
      INTEGER i,j
      DOUBLE PRECISION h,temp,f(NMAX)
      EXTERNAL funcv
      do 12 j=1,n
        temp=x(j)
        h=EPS*dabs(temp)
        if(h.eq.0.D0)h=EPS
        x(j)=temp+h
        h=x(j)-temp
        call funcv(n,x,f)
        x(j)=temp
        do 11 i=1,n
          df(i,j)=(f(i)-fvec(i))/h
11      continue
12    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 3ki-i]3|-9.
c
c
      FUNCTION fmin(x,funcv)
      INTEGER n,NP
      DOUBLE PRECISION fmin,x(*),fvec
      PARAMETER (NP=40)
      COMMON /newtv/ fvec(NP),n
      SAVE /newtv/
CU    USES funcv
      INTEGER i
      DOUBLE PRECISION sum
      EXTERNAL funcv
      call funcv(n,x,fvec)
      sum=0.D0
      do 11 i=1,n
        sum=sum+fvec(i)**2
11    continue
      fmin=0.5D0*sum
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 3ki-i]3|-9.
c
c

      SUBROUTINE lnsrch(n,xold,fold,g,p,x,f,stpmax,check,func,func1)
      INTEGER n
      LOGICAL check
      DOUBLE PRECISION f,fold,stpmax,g(n),p(n),x(n),xold(n),func,ALF,
     *TOLX
      PARAMETER (ALF=1.d-4,TOLX=1.d-7)
      EXTERNAL func
      EXTERNAL func1
CU    USES func
      INTEGER i
      DOUBLE PRECISION a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,
     *rhs2,slope,sum,temp,
     *test,tmplam
      check=.false.
      sum=0.D0
      do 11 i=1,n
        sum=sum+p(i)*p(i)
11    continue
      sum=dsqrt(sum)
      if(sum.gt.stpmax)then
        do 12 i=1,n
          p(i)=p(i)*stpmax/sum
12      continue
      endif
      slope=0.D0
      do 13 i=1,n
        slope=slope+g(i)*p(i)
13    continue
      test=0.D0
      do 14 i=1,n
        temp=dabs(p(i))/max(dabs(xold(i)),1.D0)
        if(temp.gt.test)test=temp
14    continue
      alamin=TOLX/test
      alam=1.D0
1     continue
        do 15 i=1,n
          x(i)=xold(i)+alam*p(i)
15      continue
        f=func(x,func1)
        if(alam.lt.alamin)then
          do 16 i=1,n
            x(i)=xold(i)
16        continue
          check=.true.
          return
        else if(f.le.fold+ALF*alam*slope)then
          return
        else
          if(alam.eq.1.D0)then
            tmplam=-slope/(2.D0*(f-fold-slope))
          else
            rhs1=f-fold-alam*slope
            rhs2=f2-fold2-alam2*slope
            a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
            b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/(alam-alam2)
            if(a.eq.0.D0)then
              tmplam=-slope/(2.D0*b)
            else
              disc=b*b-3.D0*a*slope
              tmplam=(-b+dsqrt(disc))/(3.D0*a)
            endif
            if(tmplam.gt..5D0*alam)tmplam=.5D0*alam
          endif
        endif
        alam2=alam
        f2=f
        fold2=fold
        alam=max(tmplam,.1D0*alam)
      goto 1
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 3ki-i]3|-9.
c
c
      SUBROUTINE qrdcmp(a,n,np,c,d,sing)
      INTEGER n,np
      DOUBLE PRECISION a(np,np),c(n),d(n)
      LOGICAL sing
      INTEGER i,j,k
      DOUBLE PRECISION scale,sigma,sum,tau
      sing=.false.
      scale=0.D0
      do 17 k=1,n-1
        do 11 i=k,n
          scale=max(scale,dabs(a(i,k)))
11      continue
        if(scale.eq.0.D0)then
          sing=.true.
          c(k)=0.D0
          d(k)=0.D0
        else
          do 12 i=k,n
            a(i,k)=a(i,k)/scale
12        continue
          sum=0.D0
          do 13 i=k,n
            sum=sum+a(i,k)**2
13        continue
          sigma=dsign(dsqrt(sum),a(k,k))
          a(k,k)=a(k,k)+sigma
          c(k)=sigma*a(k,k)
          d(k)=-scale*sigma
          do 16 j=k+1,n
            sum=0.D0
            do 14 i=k,n
              sum=sum+a(i,k)*a(i,j)
14          continue
            tau=sum/c(k)
            do 15 i=k,n
              a(i,j)=a(i,j)-tau*a(i,k)
15          continue
16        continue
        endif
17    continue
      d(n)=a(n,n)
      if(d(n).eq.0.D0)sing=.true.
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 3ki-i]3|-9.
c
c
      SUBROUTINE qrupdt(r,qt,n,np,u,v)
      INTEGER n,np
      DOUBLE PRECISION r(np,np),qt(np,np),u(np),v(np)
CU    USES rotate
      INTEGER i,j,k
      do 11 k=n,1,-1
        if(u(k).ne.0.D0)goto 1
11    continue
      k=1
1     do 12 i=k-1,1,-1
        call rotate(r,qt,n,np,i,u(i),-u(i+1))
        if(u(i).eq.0.D0)then
          u(i)=dabs(u(i+1))
        else if(dabs(u(i)).gt.dabs(u(i+1)))then
          u(i)=dabs(u(i))*dsqrt(1.D0+(u(i+1)/u(i))**2)
        else
          u(i)=dabs(u(i+1))*dsqrt(1.D0+(u(i)/u(i+1))**2)
        endif
12    continue
      do 13 j=1,n
        r(1,j)=r(1,j)+u(1)*v(j)
13    continue
      do 14 i=1,k-1
        call rotate(r,qt,n,np,i,r(i,i),-r(i+1,i))
14    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 3ki-i]3|-9.
c
c
      SUBROUTINE rsolv(a,n,np,d,b)
      INTEGER n,np
      DOUBLE PRECISION a(np,np),b(n),d(n)
      INTEGER i,j
      DOUBLE PRECISION sum
      b(n)=b(n)/d(n)
      do 12 i=n-1,1,-1
        sum=0.D0
        do 11 j=i+1,n
          sum=sum+a(i,j)*b(j)
11      continue
        b(i)=(b(i)-sum)/d(i)
12    continue
      return
      end

c      1986-92 Numerical Recipes Software 3ki-i]3|-9.
c
c
      SUBROUTINE rotate(r,qt,n,np,i,a,b)
      INTEGER n,np,i
      DOUBLE PRECISION a,b,r(np,np),qt(np,np)
      INTEGER j
      DOUBLE PRECISION c,fact,s,w,y
      if(a.eq.0.D0)then
        c=0.D0
        s=dsign(1.D0,b)
      else if(dabs(a).gt.dabs(b))then
        fact=b/a
        c=dsign(1.D0/dsqrt(1.D0+fact**2),a)
        s=fact*c
      else
        fact=a/b
        s=dsign(1.D0/dsqrt(1.D0+fact**2),b)
        c=fact*s
      endif
      do 11 j=i,n
        y=r(i,j)
        w=r(i+1,j)
        r(i,j)=c*y-s*w
        r(i+1,j)=s*y+c*w
11    continue
      do 12 j=1,n
        y=qt(i,j)
        w=qt(i+1,j)
        qt(i,j)=c*y-s*w
        qt(i+1,j)=s*y+c*w
12    continue
      return
      END

c	************************************************************
