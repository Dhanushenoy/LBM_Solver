module initial

contains

!Initialisation for D2Q9 scheme
subroutine D2Q9

      use definition
      use lecture

      implicit none

      call lect(1)

      allocate(f(0:8,0:n,0:m))
      allocate(feq(0:8,0:n,0:m))
      allocate(rho(0:n,0:m))
      allocate(u(0:n,0:m),v(0:n,0:m))
      allocate(x(0:n))
      allocate(y(0:m))


x(0)=0
y(0)=0
do i=1,n
  x(i)=x(i-1)+dx
end do

do j=1,m
  y(j)=y(j-1)+dy
end do

!for dx=dy
csk=dx*dx/(dt*dt)
omega=1.0/((3.*alpha/(csk*dt))+0.5)

Re=u0*m/alpha
print*,"Re=",Re
!mstep=20000
print*,"mstep=",mstep
print*,'omega=',omega

!weighting factor assaignment
  w(0)=4./9.
  do k=1,4
  w(k)=1./9.
  end do
  do k=5,8
  w(k)=1./36.
  end do

time=0

!Velocity components___________________________
  cx(0)=0.
  cx(1)=1.
  cx(2)=0.
  cx(3)=-1.
  cx(4)=0.
  cx(5)=1.
  cx(6)=-1.
  cx(7)=-1.
  cx(8)=1.
  cy(0)=0.
  cy(1)=0.
  cy(2)=1.
  cy(3)=0.
  cy(4)=-1.
  cy(5)=1.
  cy(6)=1.
  cy(7)=-1.
  cy(8)=-1.

!Rho init
do i=0,n
 do j=0,m
     rho(i,j)=rhou
 end do
 end do



 do i=0,n
 do j=0,m
 u(i,j)=0.0
 v(i,j)=0.0
 end do
 end do


 do j=1,m-1
 u(0,j)=u0
 v(0,j)=0.
 end do


  do j=0,m
      do i=0,n
          do k=0,8

          f(k,i,j)=w(k)*rho(i,j)

          end do
      end do
  end do

print*,"end of D2Q9"
end subroutine D2Q9

end module initial
