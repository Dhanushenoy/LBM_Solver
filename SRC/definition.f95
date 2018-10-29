module definition

  implicit none
  integer::n,m,time,mstep
  real,allocatable::f(:,:,:),feq(:,:,:)
  real,allocatable::rho(:,:)
  real,dimension(0:8)::w,cx,cy
  real,allocatable::u(:,:),v(:,:)
  real,allocatable::x(:),y(:)
  integer::i,j,k
  real::Re    !Lattice Reynold's
  real::Alpha,rhou   !Viscousity
  real::csk,omega
  real::dx,dy,dt
  real::u0
  real::s1,s2
  
!MRT Schemes
  integer::l
  real::a1,sumcc,tau
  real::suma,sumb
  real,dimension(0:8,0:8)::ev,tm,tminv,stmiv
  real,dimension(0:8)::sm
  real,allocatable::fmom(:,:,:)


end module definition
