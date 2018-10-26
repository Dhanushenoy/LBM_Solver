program LBM

use definition
use lecture
use initial
use operation
use out

implicit none

call lect(1)
call D2Q9
!call VTKout(3)

do while(time<=mstep)

call collision
call streaming
call boundary
call Rhouv


time=time+dt
if(mod(time,100)==0)then
print*,"time=",time
end if
end do

call VTKout(100)
!call output(10)
end program LBM
