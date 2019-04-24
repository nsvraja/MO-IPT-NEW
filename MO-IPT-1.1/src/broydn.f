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
          if(sing) pause 'singular Jacobian in broydn'
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
              if(r(i,i).eq.0.D0) pause 'r singular in broydn'
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
