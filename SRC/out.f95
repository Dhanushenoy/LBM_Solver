module out

contains
subroutine output(out)
use definition
use lecture
use initial
use operation

implicit none
integer::out
real::gradx,grady,gradm
real::mag
real ::strf(0:n,0:m)
real::rhoav,rhom
real::R(0:n-1,0:m-1)
real::vort
real::Phi,u3
real::magV(0:n,0:m)


!open(unit=5,file='streamf')
!strf(0,0)=0
!do i=0,n
!  rhoav=0.5*(rho(i-1,0)+rho(i,0))
!  if (i.ne.0) then
!  strf(i,0)=strf(i-1,0)-rhoav*0.5*(v(i-1,0)+v(i,0))
!endif
!  do j=1,m
!  rhom=0.5*(rho(i,j)+rho(i,j-1))
!  strf(i,j)=strf(i,j-1)+rhom*0.5*(u(i,j-1)+u(i,j))
!  end do
!  end do

!open(unit=10,file='gradient')

!do i=0,n-1
!do j=0,m-1
!    u3=sqrt(u(i,j)*u(i,j)+v(i,j)*v(i,j))
!    gradx=u(i,j)*(u(i+1,j)-u(i,j))
!    grady=v(i,j)*(v(i,j+1)-v(i,j))
!    R(i,j)=u3*u3*u3/(u(i,j)*grady-v(i,j)*gradx)
!    vort=(u(i+1,j)-u(i,j))-(v(i,j+1)-v(i,j))

!    Phi=2*u3*vort/R(i,j)

!    if (Phi<0)then
!        Phi=-1
!    elseif (Phi==0)then
!    Phi=0
!    else
!        Phi=1
!    end if

!write(10,*)i,j,Phi
!end do
!write(10,*)
!end do


!open(unit=2,file='all')
!open(unit=4,file='floats')
!write(2,*)"VARIABLES =X, Y, U, V, S"
!write(2,*)"ZONE","I=",n+1,"J=",m+1,"","","F=BLOCK"
!do j=0,m
!write(2,*)(i,i=0,n)
!end do
!do j=0,m
!write(2,*)(j,i=0,n)
!end do
!do j=0,m
!write(2,*)(u(i,j),i=0,n)
!end do
!do j=0,m
!write(2,*)(v(i,j),i=0,n)
!end do
!do j=0,m
!write(2,*)(strf(i,j),i=0,n)
!end do
!do j=0,m
!write(3,*)j/float(m),u(n/2,j)/u0,u(n/4,j)/u0,u(3*n/4,j)/u0
!end do
!do i=0,n
!write(4,*) i/float(n),v(i,m/2)/u0
!end do

open(unit=10,file='output')
write(10,*)'x,y,u,v,mag'
do i=0,n
do j=0,m
magV(i,j)=(u(i,j)*u(i,j)+v(i,j)*v(i,j))**0.5
write(10,*)i,j,u(i,j),v(i,j),magV(i,j)
end do
write(10,*)
end do

open(unit=15,file='plot25')
open(unit=16,file='plot100')
open(unit=17,file='plot250')

do i=0,n
  if (i==25) then
    do j=0,m
    write(15,*)j,u(i,j)/u0
    end do

  elseif (i==100) then
      do j=0,m
      write(16,*)j,u(i,j)/u0
      end do

    elseif (i==250) then
        do j=0,m
        write(17,*)j,u(i,j)/u0
        end do
  end if

end do
end subroutine output



subroutine VTKout(vtk)
  use definition
  use lecture
  use initial
  use operation

  implicit none
  integer::vtk
  character(len=10) :: String1
  real::magV

open(unit=vtk,file='out.vtk')!,form='formatted',access='stream',status='replace')

write(vtk,'(A)')'# vtk DataFile Version 3.0'
write(vtk,'(A)')'output.vtk'
write(vtk,'(A)')'ASCII'
write(vtk,'(A)')'DATASET STRUCTURED_GRID'
write(vtk,'(A,I5,I5,I5)')'DIMENSIONS',size(x),size(y),1
write(vtk,'(A,I10,A)')'POINTS',(size(x)*size(y)),' float'
  do i=0,n
    do j=0,m
write(vtk,'(F10.2,F10.2,F10.2)') x(i),y(j),1.0
    end do
  end do

write(vtk,'(A,I10)')'POINT_DATA',(size(x)*size(y))
write(vtk,'(A)')'SCALARS mag_velocity float'
write(vtk,'(A)')'LOOKUP_TABLE table_u'
do i=0,n
  do j=0,m
    magV=(u(i,j)*u(i,j)+v(i,j)*v(i,j))**0.5
write(vtk,'(F16.10)')u(i,j)
  end do
end do

!write(vtk,'(A)')'SCALARS v_velocity float'
!write(vtk,'(A)')'LOOKUP_TABLE my_table 1.0'
!do i=0,n
!  do j=0,m
!write(vtk,'(F16.10)')v(i,j)
!  end do
!end do


end subroutine VTKout
end module out
