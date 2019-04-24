!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         real*8 function  broyden(tmp_muc)
         use Global
         implicit none
         include'mkl_rci.fi'
         integer M,L,K
         PARAMETER(L = 1)
         PARAMETER(M = 1)
         real*8 x(L),fjac(m,L),eps,f(m),D(L),U_(L)
         real*8 C,Y(L,L),Z(L,L),xold(L),tmp_muc
         integer ipiv(N),info,lwork
         PARAMETER(lwork=L*L)
         integer work(lwork)
         external fcn
         eps=1.0d-8
         x(1)=tmp_muc
         call fcn(m,l,x,f)

         if(DJACOBI (fcn,l,m,fjac,x,eps)
     .  /=TR_SUCCESS)then
        write(*,*) '|ERROR IN DJACOBI'
        endif

        call dgetrf(M,L,fjac,L,ipiv,info)

        if(info==0) then
        call dgetri(M,fjac,L,ipiv,work,lwork,info)
        else
        write(*,*) "problem with matrix inversion",info
        endif

        IF (info.NE.0) THEN
        stop 'Matrix inversion failed!'
        ELSE
        PRINT '(" Inverse Successful ")'
        ENDIF 


        D=-matmul(fjac,f) 

        x=x+D

!        write(*,"(2a,3F15.4)") "X values", "0", x


        k=0

        do
        call fcn(m,l,x,f)
        U_=matmul(fjac,f)
        C=dot_product(D,U_+D)

        Y=0.d0;Z=0.d0

        Y(1,:)=D
        Z(1,:)=U_

       Y=transpose(Y)
       fjac=fjac-(1/C)*matmul(matmul(Z,Y),fjac)
       k=k+1
       D=-matmul(fjac,f)
       xold=x
       x=x+D
       if(dot_product(x-xold,x-xold)<10.d-7) exit
       enddo

        broyden=x(1)

        

         return
         end 

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         subroutine fcn(size_f,size_x,var,ff)
         use Global
         implicit none
         integer :: size_f,size_x,io
         real*8 var(size_x),ff(size_f)
         real*8 nf_tot_tmp
         include 'Glob_cons'
         mu_c=var(1)
         call findgf()
         call occupancies()
          nf_tot_tmp=0.d0
          do io=1,Norbs
          nf_tot_tmp=nf_tot_tmp+nf(io)
         end do

         ff(1)=ntot-nf_tot_tmp
         write(*,*)'funcv1 mu_c=',mu_c,' fvec',ff(1)              

        return
        end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
