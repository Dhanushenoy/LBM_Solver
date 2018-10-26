module lecture

contains
  subroutine lect(nfile)

    use definition

    implicit none
    integer::nfile

    open(nfile,file='input',form='formatted')
    read(nfile,*)n,m
    print*,n,m
    read(nfile,*)dx,dy,dt
    print*,dx,dy,dt
    read(nfile,*)Alpha
    print*,Alpha
    read(nfile,*)u0
    print*,u0
    read(nfile,*)rhou
    print*,rhou
    read(nfile,*)mstep
    print*,mstep

    close(nfile)

  end subroutine lect




end module lecture
