
subroutine convolve_fft(n,w,dw,rho1,rho2,rho3)
integer :: n
real*8 :: w(-n:n),rho1(-n:n), rho2(-n:n), rho3(-n:n),dw
integer :: nft,m,n2
complex*16, allocatable :: z1(:),z2(:),z3(:)
complex*16,parameter :: zero=(0.d0,0.d0)

n2=2*(2*n+1)
nft=2
m=1
do while(nft < n2)
 nft=nft*2
 m=m+1
end do
nfti=nft-1

!write(6,*) 'nft=',nft,' m=',m
allocate(z1(0:nft))
allocate(z2(0:nft))
allocate(z3(0:nft))

do i=0,nft
 z1(i)=zero
 z2(i)=zero
 z3(i)=zero
end do

do i=-n,n
  z1(i+n+1)=dcmplx(rho1(i),0.d0)
  z2(i+n+1)=dcmplx(rho2(i),0.d0)
end do

call dffteu(nft,z1,m,0)
call dffteu(nft,z2,m,0)

do i=1,nfti
  z3(i)=dw*z1(nft-i)*z2(i)
end do
z3(0)=dw*z1(0)*z2(0)

call dffteu(nft,z3,m,-1)

do i=-n,-1
  rho3(i) = dreal(z3(nfti+i+1))
  write(11,*) w(i),rho3(i)
end do
do i=0,n
  rho3(i)=dreal(z3(i))
!  write(11,*) w(i),rho3(i)
end do

!close(11)

deallocate(z1)
deallocate(z2)
deallocate(z3)
return
end subroutine convolve_fft



!-------------------------------------------------------------c
!                                                             c
!  Subroutine dffteu( n, z, m, itype )                     c
!                                                             c
!  This routine is a slight modification of a complex split   c
!  radix FFT routine presented by C.S. Burrus.  The original  c
!  program header is shown below.                             c
!                                                             c
!  Arguments:                                                 c
!     z - complex*16 array containing transform sequence
!     x - real*8 array containing real*8 parts of transform   c
!              sequence (in/out)                              c
!     y - real*8 array containing imag parts of transform     c
!              sequence (in/out)                              c
!     n - integer length of transform (in)                    c
!     m - integer such that n = 2**m  (in)                    c
!     itype - integer job specifier (in)                      c
!              itype .ne. -1 --> foward transform             c
!              itype .eq. -1 --> backward transform           c
!                                                             c
!  The forward transform computes                             c
!     X(k) = sum_{j=0}^{N-1} x(j)*exp(-2ijk*pi/N)             c
!                                                             c
!  The backward transform computes                            c
!     x(j) = (1/N) * sum_{k=0}^{N-1} X(k)*exp(2ijk*pi/N)      c
!                                                             c
!                                                             c
!  Requires standard FORTRAN functions - DSIN, DCOS             c
!                                                             c
!  Steve Kifowit, 9 July 1997                                 c
!                                                             c
!-------------------------------------------------------------C
!  A Duhamel-Hollman Split-Radix DIF FFT                      C
!  Reference:  Electronics Letters, January 5, 1984           C
!  Complex input and output in data arrays X and Y            C
!  Length is N = 2**M                                         C
!                                                             C
!  C.S. Burrus          Rice University         Dec 1984      C
!-------------------------------------------------------------C
!
      SUBROUTINE DFFTEU( N, Z, M, ITYPE )
      INTEGER  ::  N, M, ITYPE
      real*8 :: E, A, A3, CC1, SS1, CC3, SS3
      real*8 ::  R1, R2, S1, S2, S3, XT,T1,T2
      INTRINSIC  DSIN, DCOS
      COMPLEX*16 :: Z(N)
      REAL*8 :: X(N), Y(N)
      INTEGER::  I, J, K, N1, N2, N4, IS, ID, I0, I1, I2, I3
      REAL*8, PARAMETER ::   TWOPI = 6.2831853071795864769d0 

!      CALL CPU_TIME(T1)
      !WRITE(*,*)"Entering DFFTEU"
!
      IF ( N .EQ. 1 ) RETURN
      !IF ( NMAX .LT. N ) RETURN
!
      DO I = 1, N
	    X(I) = DREAL(Z(I))
	    Y(I) = DIMAG(Z(I))
      END DO
!
      IF ( ITYPE .EQ. -1 ) THEN
	 DO 1, I = 1, N
	    Y(I) = - Y(I)
 1       CONTINUE
      ENDIF
!
      N2 = 2 * N
      DO 10, K = 1, M-1
	 N2 = N2 / 2
	 N4 = N2 / 4
	 E = TWOPI / N2
	 A = 0.0
	 DO 20, J = 1, N4
	    A3 = 3 * A
	    CC1 = DCOS( A )
	    SS1 = DSIN( A )
	    CC3 = DCOS( A3 )
	    SS3 = DSIN( A3 )
	    A = J * E
	    IS = J
	    ID = 2 * N2
 40         DO 30, I0 = IS, N-1, ID
	       I1 = I0 + N4
	       I2 = I1 + N4
	       I3 = I2 + N4
	       R1 = X(I0) - X(I2)
	       X(I0) = X(I0) + X(I2)
	       R2 = X(I1) - X(I3)
	       X(I1) = X(I1) + X(I3)
	       S1 = Y(I0) - Y(I2)
	       Y(I0) = Y(I0) + Y(I2)
	       S2 = Y(I1) - Y(I3)
	       Y(I1) = Y(I1) + Y(I3)
	       S3 = R1 - S2
	       R1 = R1 + S2
	       S2 = R2 - S1
	       R2 = R2 + S1
	       X(I2) = R1 * CC1 - S2 * SS1
	       Y(I2) = - S2 * CC1 - R1 * SS1
	       X(I3) = S3 * CC3 + R2 * SS3
	       Y(I3) = R2 * CC3 - S3 * SS3
 30         CONTINUE
	    IS = 2 * ID - N2 + J
	    ID = 4 * ID
	    IF ( IS .LT. N ) GOTO 40
 20      CONTINUE
 10   CONTINUE
!
!--------LAST STAGE, LENGTH-2 BUTTERFLY ----------------------C
!
      IS = 1
      ID = 4
 50   DO 60, I0 = IS, N, ID
	 I1 = I0 + 1
	 R1 = X(I0)
	 X(I0) = R1 + X(I1)
	 X(I1) = R1 - X(I1)
	 R1 = Y(I0)
	 Y(I0) = R1 + Y(I1)
	 Y(I1) = R1 - Y(I1)
 60   CONTINUE
      IS = 2 * ID - 1
      ID = 4 * ID
      IF ( IS .LT. N ) GOTO 50
!
!-------BIT REVERSE COUNTER-----------------------------------C
!
 100  J = 1
      N1 = N - 1
      DO 104, I = 1, N1
	 IF ( I .GE. J ) GOTO 101
	 XT = X(J)
	 X(J) = X(I)
	 X(I) = XT
	 XT = Y(J)
	 Y(J) = Y(I)
	 Y(I) = XT
 101     K = N / 2
 102     IF ( K .GE. J ) GOTO 103
	 J = J - K
	 K = K / 2
	 GOTO 102
 103     J = J + K
 104  CONTINUE
!
      IF ( ITYPE .EQ. -1 ) THEN
	 DO 2, I = 1, N
	    X(I) = X(I) / N
	    Y(I) = - Y(I) / N
 2       CONTINUE
      ENDIF
!
      DO I = 1, N
	Z(I)=DCMPLX(X(I),Y(I))
      END DO

      

!      CALL CPU_TIME(T2)
!      WRITE(*,*)"DFFTEU took",T2-T1,"seconds"

      RETURN
!
! ... End of subroutine DFFTEU ...
!
      END
