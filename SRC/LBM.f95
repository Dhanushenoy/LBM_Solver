program LBM






use definition
use lecture
use initial
use operation
use out

implicit none

!******For cpu time*********!
real::start,finish          !
call cpu_time(start)        !
!***************************!

call lect(1)

!**********************************|
! D2Q9 is used for 2d flows        |
! 3D schemes will be used later    |
!**********************************|
call D2Q9

!Calculation loop ---------------
do while(time<=mstep)

!***********************************|
! The collision step contains       |
! both MRT not sigle relaxation     |
! Just use collsion or collisionMRT |
!***********************************|
call collision

call streaming
!Boundary conditions----------------
call boundary

!Collecting all the terms-----------
call Rhouv


time=time+dt
if(mod(time,100)==0)then
print*,"time=",time
end if
end do

call VTKout(100)
call output(10)

!***********CPU time ****************************!
call cpu_time(finish)                            !
print '("time=",f16.10," seconds")',finish-start !
!************************************************!
end program LBM
