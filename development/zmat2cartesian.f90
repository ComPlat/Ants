!
! This code is adapted from the xtb source code.
! Original source: https://github.com/grimme-lab/xtb/src/intmodes.f90
! License: LGPL-3 | GPL-3
!


pure subroutine gmetry(natoms, geo, coord, na,nb,nc,fail)
  use iso_c_binding
  implicit none

  integer, intent(in)    :: natoms
  real(c_double),intent(inout) :: coord(3,natoms)
  real(c_double),intent(in)    :: geo(3,natoms)
  integer, intent(in)    :: na(natoms), nb(natoms), nc(natoms)
  logical, intent(out)   :: fail
  real(c_double) :: ccos
  real(c_double) :: xa,ya,za
  real(c_double) :: xb,yb,zb
  real(c_double) :: xd,yd,zd
  real(c_double) :: rbc,xpa,xpb,xyb
  real(c_double) :: ypa,xqa,zqa,yza
  real(c_double) :: ypd,zpd,xpd,zqd,xqd,yqd,xrd
  real(c_double) :: costh,sinth,sinph,cosph
  real(c_double) :: coskh,sinkh,sina,cosa,sind,cosd
  integer  :: i,k
  integer  :: ma,mb,mc

  fail=.false.

  coord(1,1)=0.0d00
  coord(2,1)=0.0d00
  coord(3,1)=0.0d00
  coord(1,2)=geo(1,2)
  coord(2,2)=0.0d00
  coord(3,2)=0.0d00
  if(natoms.eq.2) goto 110
  ccos=cos(geo(2,3))
  if(na(3).eq.1)then
    coord(1,3)=coord(1,1)+geo(1,3)*ccos
  else
    coord(1,3)=coord(1,2)-geo(1,3)*ccos
  endif
  coord(2,3)=geo(1,3)*sin(geo(2,3))
  coord(3,3)=0.0d00
  do i=4,natoms
    cosa=cos(geo(2,i))
    mb=nb(i)
    mc=na(i)
    xb=coord(1,mb)-coord(1,mc)
    yb=coord(2,mb)-coord(2,mc)
    zb=coord(3,mb)-coord(3,mc)
    rbc=1.0d00/dsqrt(xb*xb+yb*yb+zb*zb)
    if (abs(cosa).lt.0.9999999999d0) go to 40
    !
    !     atoms mc, mb, and (i) are collinear
    !
    rbc=geo(1,i)*rbc*cosa
    coord(1,i)=coord(1,mc)+xb*rbc
    coord(2,i)=coord(2,mc)+yb*rbc
    coord(3,i)=coord(3,mc)+zb*rbc
    cycle
    !
    !     the atoms are not collinear
    !
    40    ma=nc(i)
    xa=coord(1,ma)-coord(1,mc)
    ya=coord(2,ma)-coord(2,mc)
    za=coord(3,ma)-coord(3,mc)
    !
    !     rotate about the z-axis to make yb=0, and xb positive.  if xyb is
    !     too small, first rotate the y-axis by 90 degrees.
    !
    xyb=dsqrt(xb*xb+yb*yb)
    k=-1
    if (xyb.gt.0.1d00) go to 50
    xpa=za
    za=-xa
    xa=xpa
    xpb=zb
    zb=-xb
    xb=xpb
    xyb=dsqrt(xb*xb+yb*yb)
    k=1
    !
    !     rotate about the y-axis to make zb vanish
    !
    50    costh=xb/xyb
    sinth=yb/xyb
    xpa=xa*costh+ya*sinth
    ypa=ya*costh-xa*sinth
    sinph=zb*rbc
    cosph=dsqrt(abs(1.d00-sinph*sinph))
    xqa=xpa*cosph+za*sinph
    zqa=za*cosph-xpa*sinph
    !
    !     rotate about the x-axis to make za=0, and ya positive.
    !
    yza=dsqrt(ypa*ypa+zqa*zqa)
    if(yza.lt.1.d-4 )then
      if(yza.lt.1.d-4)goto 70
      !           write(*,'(/9x,'' atoms'',i3,'','',i3,'', and'',i3,
      !    1'' are within'',f7.4,'' angstroms of a straight line'')')
      !    2mc,mb,ma,yza
      fail=.true.
      return
    endif
    coskh=ypa/yza
    sinkh=zqa/yza
    goto 80
    70    continue
    !
    !   angle too small to be important
    !
    coskh=1.d0
    sinkh=0.d0
    80    continue
    !
    !     coordinates :-   a=(xqa,yza,0),   b=(rbc,0,0),  c=(0,0,0)
    !     none are negative.
    !     the coordinates of i are evaluated in the new frame.
    !
    sina=sin(geo(2,i))
    sind=-sin(geo(3,i))
    cosd=cos(geo(3,i))
    xd=geo(1,i)*cosa
    yd=geo(1,i)*sina*cosd
    zd=geo(1,i)*sina*sind
    !
    !     transform the coordinates back to the original system.
    !
    ypd=yd*coskh-zd*sinkh
    zpd=zd*coskh+yd*sinkh
    xpd=xd*cosph-zpd*sinph
    zqd=zpd*cosph+xd*sinph
    xqd=xpd*costh-ypd*sinth
    yqd=ypd*costh+xpd*sinth
    if (k.lt.1) go to 90
    xrd=-zqd
    zqd=xqd
    xqd=xrd
    90    coord(1,i)=xqd+coord(1,mc)
    coord(2,i)=yqd+coord(2,mc)
    coord(3,i)=zqd+coord(3,mc)
  end do
  110 continue

end subroutine gmetry


subroutine zmat2cart(n,n2,at,geo,xyz,na,nb,nc,fail)
  use iso_c_binding
  implicit none
  integer, intent(in) :: n,n2
  integer, intent(in) :: na(n2),nb(n2),nc(n2)
  integer, intent(in) :: at(n2)
  real(c_double),intent(inout) :: xyz(3,n)
  real(c_double),intent(in) :: geo(3,n2)
  logical, intent(out) :: fail

  integer k,m
  real(c_double),allocatable :: tmp(:,:)
  allocate(tmp(3,n2))

  call gmetry(n2,geo,tmp,na,nb,nc,fail)

  m=0
  do k=1,n2
    if(at(k).lt.100) then
      m=m+1
      xyz(1:3,m)=tmp(1:3,k)
    endif
  enddo

end subroutine

subroutine get_element_symbol(Z, symbol)
   implicit none
   integer, intent(in) :: Z
   character(len=2), intent(out) :: symbol
   character(len=2), dimension(1:6) :: table = ['H ','C ','N ','O ','F ','S ']
   if (Z >= 1 .and. Z <= 6) then
      symbol = table(Z)
   else
      symbol = 'X '
   end if
end subroutine

subroutine write_xyz(filename, nat, at, xyz)
   use iso_c_binding
   implicit none
   character(len=*), intent(in) :: filename
   integer, intent(in) :: nat
   integer, intent(in) :: at(nat)
   real(c_double), intent(in) :: xyz(3, nat)
   integer :: i
   character(len=2) :: symbol
   integer :: unit

   open(newunit=unit, file=filename, status='replace', action='write')
   write(unit, *) nat
   write(unit, *) 'Generated from zmat2cart'

   do i = 1, nat
      call get_element_symbol(at(i), symbol)
      write(unit, '(A2,3F12.6)') trim(symbol), xyz(1,i), xyz(2,i), xyz(3,i)
   end do

   close(unit)
end subroutine

program test

use iso_c_binding
implicit none
integer :: n ! atomic numbers of atoms with atomic number < 100
integer :: n2 ! total number of atoms in the molecule
integer, allocatable :: at(:) ! atomic numbers for all n2
real(c_double), allocatable :: geo(:,:) ! row 1 = bond lengths; row 2 = bond angles; row 3 = dihedral angles
integer, allocatable :: na(:),nb(:),nc(:) ! atom indices defining the internal coordinates for each atom
! na(i) atom to which atom i is bound (= bond length)
! nb(i) containing the angle between atom i and nb(i)
! nc(i) = atom defining dihedral with nb(i) (defines dihedral angle)
logical :: fail = .false.
real(c_double), allocatable :: xyz(:,:)
integer :: i
character(len=10) :: istr

n = 18
n2 = 18
allocate(at(n2), na(n2), nb(n2), nc(n2))
allocate(geo(3, n2))
allocate(xyz(3, n2))

! Atomic numbers: C = 6, H = 1
at = [6, 6, 6, 1, 1, 6, 1, 1, 6, 1, 1, 6, 1, 1, 1, 1, 1, 1]

! Internal coordinate references (Fortran is 1-based):
na = [0, 1, 2, 2, 2, 3, 3, 3, 6, 6, 6, 9, 9, 9,12,12, 1, 1]
nb = [0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 6, 6, 6, 9, 9, 2, 2]
nc = [0, 0, 0, 3, 3, 1, 1, 1, 2, 2, 2, 3, 3, 3, 6, 6, 3, 3]

! Geometry input (bond length, angle [deg], dihedral [deg]):
geo(:,1) = [0.0d0, 0.0d0, 0.0d0]  ! First atom: position is defined arbitrarily
geo(:,2) = [1.5240d0, 0.0d0, 0.0d0]  ! Second atom: just bond length
geo(:,3) = [1.5532d0, 109.83d0, 0.0d0] ! Third atom: bond + angle
geo(:,4) = [1.1109d0, 113.14d0, 237.32d0]
geo(:,5) = [1.1130d0, 106.46d0, 119.29d0]
geo(:,6) = [1.5532d0, 113.52d0, 30.94d0]
geo(:,7) = [1.1115d0, 111.38d0, 270.84d0]
geo(:,8) = [1.1117d0, 106.43d0, 153.81d0]
geo(:,9) = [1.5240d0, 109.74d0, 32.14d0]
geo(:,10)= [1.1130d0, 110.20d0, 275.07d0]
geo(:,11)= [1.1108d0, 109.49d0, 156.86d0]
geo(:,12)= [1.5532d0, 109.83d0, 295.79d0]
geo(:,13)= [1.1109d0, 113.14d0, 173.11d0]
geo(:,14)= [1.1130d0, 106.46d0, 55.08d0]
geo(:,15)= [1.1115d0, 111.38d0, 270.84d0]
geo(:,16)= [1.1117d0, 106.43d0, 153.81d0]
geo(:,17)= [1.1130d0, 106.58d0, 55.11d0]
geo(:,18)= [1.1108d0, 113.12d0, 173.20d0]

! conversion from degree to radian
do i = 4, n2
  geo(2,i) = geo(2,i) * (acos(-1.0d0) / 180.0d0)
  geo(3,i) = geo(3,i) * (acos(-1.0d0) / 180.0d0)
end do

call zmat2cart(n,n2,at,geo,xyz,na,nb,nc,fail)

print *, 'fail =', fail
do i = 1, n
   write(istr, '(I0)') i  ! convert integer to string
   print '(A,3F12.6)', 'Atom '//trim(istr)//': ', xyz(1,i), xyz(2,i), xyz(3,i)
end do

call write_xyz("test.xyz", n2, at, xyz)

deallocate(at, na, nb, nc)
deallocate(geo, xyz)

end program test