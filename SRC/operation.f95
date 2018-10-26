module operation

contains
subroutine collision
  use definition
  use lecture
  use initial

  do i=0,n
   do j=0,m

     s1=u(i,j)*u(i,j)+v(i,j)*v(i,j)

     do k=0,8

        s2=cx(k)*u(i,j)+cy(k)*v(i,j)
        feq(k,i,j)=rho(i,j)*w(k)*(1.+3.*s2+4.5*s2*s2-1.5*s1)
        f(k,i,j)=feq(k,i,j)*omega+(1.-omega)*f(k,i,j)

      end do
   end do
 end do
end subroutine collision

subroutine streaming
  use definition
  use lecture
  use initial
  implicit none
  !streaming
      do i=n,1,-1
          do j=0,m
      f(1,i,j)=f(1,i-1,j)
          end do
      end do

          do i=0,n
          do j=m,1,-1
      f(2,i,j)=f(2,i,j-1)
          end do
      end do

          do i=0,n-1
          do j=0,m
      f(3,i,j)=f(3,i+1,j)
          end do
      end do

        do i=0,n
          do j=0,m-1
      f(4,i,j)=f(4,i,j+1)
          end do
      end do

      do i=n,1,-1
          do j=m,1,-1
      f(5,i,j)=f(5,i-1,j-1)
          end do
      end do

      do i=0,n-1
          do j=m,1,-1
      f(6,i,j)=f(6,i+1,j-1)
          end do
      end do

      do i=0,n-1
          do j=0,m-1
      f(7,i,j)=f(7,i+1,j+1)
          end do
      end do

      do i=n,1,-1
          do j=0,m-1
      f(8,i,j)=f(8,i-1,j+1)
          end do
      end do
  end subroutine streaming

  subroutine boundary

    use definition
    use lecture
    use initial

implicit none
real::rhow

!-----------west inlet------------------------
    do j=0,m
!For inlet velocity (pipe flow)
!  rhow=(f(0,0,j)+f(2,0,j)+f(4,0,j)+2.*(f(3,0,j)+f(6,0,j)+f(7,0,j)))/(1.-u0)
!  f(1,0,j)=f(3,0,j)+2.*rhow*u0/3.
!  f(5,0,j)=f(7,0,j)+rhow*u0/6.-(0.5*(f(2,0,j)-f(4,0,j)))
!  f(8,0,j)=f(6,0,j)+rhow*u0/6.+(0.5*(f(2,0,j)-f(4,0,j)))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!bounce back on west
  f(1,0,j)=f(3,0,j)
  f(5,0,j)=f(7,0,j)
  f(8,0,j)=f(6,n,j)
  end do



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!----------------South-------------------------!
  !bounce on south phase
  do i=0,n
  f(2,i,0)=f(4,i,0)
  f(6,i,0)=f(8,i,0)
  f(5,i,0)=f(7,i,0)
  end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----------------north-------------------------
  do i=1,n-1
  rhow=f(0,i,m)+f(1,i,m)+f(3,i,m)+2.*(f(2,i,m)+f(6,i,m)+f(5,i,m))
  f(4,i,m)=f(2,i,m)
  f(8,i,m)=f(6,i,m)+rhow*u0/6.0
  f(7,i,m)=f(5,i,m)-rhow*u0/6.0
  end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------EAST--------------------------
  do j=0,m

    !open condition east
    !ux=-1+(f(0,n,j)+f(2,n,j)+f(4,n,j)+2.*(f(1,n,j)+f(5,n,j)+f(8,n,j)))/rhou
    !f(3,n,j)=f(1,n,j)-2.*rhou*ux/3.
    !f(7,n,j)=f(5,n,j)+0.5*(f(2,n,j)-f(4,n,j))-(rhou*ux)/6.
    !f(6,n,j)=f(8,n,j)-0.5*(f(2,n,j)-f(4,n,j))-(rhou*ux)/6.
    !Interpolation scheme
!  f(1,n,j)=2.*f(1,n-1,j)-f(1,n-2,j)
!  f(5,n,j)=2.*f(5,n-1,j)-f(5,n-2,j)
!  f(8,n,j)=2.*f(8,n-1,j)-f(8,n-2,j)
!Bounceback
  f(3,n,j)=f(1,n,j)
  f(7,n,j)=f(5,n,j)
  f(6,n,j)=f(8,n,j)
  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--------------obstacle----------------------
! Mention the coordinates of obstacle in 2D
! For now only simply geometries
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine boundary

subroutine Rhouv

  use definition
  use lecture
  use initial

  implicit none
    real::usum,vsum,tot
  !rhouv
    do i=0,n
    do j=0,m
    tot=0.0
    do k=0,8
    tot=tot+f(k,i,j)
    end do
    rho(i,j)=tot
    end do
    end do


    do i=0,n
    do j=0,m
    usum=0.0
    vsum=0.0
    do k=0,8
    usum=usum+f(k,i,j)*cx(k)
    vsum=vsum+f(k,i,j)*cy(k)
    end do
    u(i,j)=usum/rho(i,j)
    v(i,j)=vsum/rho(i,j)

    end do
    end do

!    do j=0,m
!    v(n,j)=0.0
!    end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--------obstacle force velocity---------
! Set the velocity inside obj as zero
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    do j=15,25
!    do i=100,110
!    u(i,j)=0.0
!    v(i,j)=0.0
!    end do
!    end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end subroutine Rhouv


end module operation
