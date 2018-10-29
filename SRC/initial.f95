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
      allocate(fmom(0:8,0:n,0:m))
      allocate(rho(0:n,0:m))
      allocate(u(0:n,0:m),v(0:n,0:m))
      allocate(x(0:n))
      allocate(y(0:m))


x(0)=0
y(0)=0
do i=1,n
  x(i)=x(i-1)+dx
!x(i)=i*dx/n
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

w(:)=(/4./9.,1./9.,1./9.,1./9.,1./9.,1./36.,1./36.,1./36.,1./36./)
time=0

!Velocity components___________________________
cx(:)=(/0.,1.,0.,-1.,0.,1.,-1.,-1.,1./)
cy(:)=(/0.,0.,1.,0.,-1.,1.,1.,-1.,-1./)

!MRT Schemes-----------------------------------
tm(0,:)=(/1.,1.,1.,1.,1.,1.,1.,1.,1./)
tm(1,:)=(/-4.,-1.,-1.,-1.,-1.,2.,2.,2.,2./)
tm(2,:)=(/4.,-2.,-2.,-2.,-2.,1.,1.,1.,1./)
tm(3,:)=(/0.,1.,0.,-1.,0.,1.,-1.,-1.,1./)
tm(4,:)=(/0.,-2.,0.,2.,0.,1.,-1.,-1.,1./)
tm(5,:)=(/0.,0.,1.,0.,-1.,1.,1.,-1.,-1./)
tm(6,:)=(/0.,0.,-2.,0.,2.,1.,1.,-1.,-1./)
tm(7,:)=(/0.,1.,-1.,1.,-1.,0.,0.,0.,0./)
tm(8,:)=(/0.,0.,0.,0.,0.,1.,-1.,1.,-1./)

!MRT inverse matrix----------------------------
a1=1./36.
tminv(0,:)=(/4.*a1,-4.*a1,4.*a1,0.,0.,0.,0.,0.,0./)
tminv(1,:)=(/4.*a1,-a1,-2.*a1,6.*a1,-6.*a1,0.,0.,9.*a1,0./)
tminv(2,:)=(/4.*a1,-a1,-2.*a1,0.,0.,6.*a1,-6.*a1,-9.*a1,0./)
tminv(3,:)=(/4.*a1,-a1,-2.*a1,-6.*a1,6.*a1,0.,0.,9.*a1,0./)
tminv(4,:)=(/4.*a1,-a1,-2.*a1,0.,0.,-6.*a1,6.*a1,-9.*a1,0./)
tminv(5,:)=(/4.*a1,2.*a1,a1,6.*a1,3.*a1,6.*a1,3.*a1,0.,9.*a1/)
tminv(6,:)=(/4.*a1,2.*a1,a1,-6.*a1,-3.*a1,6.*a1,3.*a1,0.,-9.*a1/)
tminv(7,:)=(/4.*a1,2.*a1,a1,-6.*a1,-3.*a1,-6.*a1,-3.*a1,0.,9.*a1/)
tminv(8,:)=(/4.*a1,2.*a1,a1,6.*a1,3.*a1,-6.*a1,-3.*a1,0.,-9.*a1/)

do i = 0,8
  do j=0,8
    sumcc=0.0
    do l=0,8
      sumcc=sumcc+tminv(i,l)*tm(l,j)
    end do
    ev(i,j)=sumcc
  end do
end do

tau=1./omega
sm(:)=(/1.0,1.4,1.4,1.0,1.2,1.0,1.2,tau,tau/)

do i=0,8
  do j=0,8
    stmiv(i,j)=tminv(i,j)*sm(j)
  end do
end do







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
