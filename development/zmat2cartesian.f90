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




program test_zmat2cart
use iso_c_binding
use xtb_intmodes
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
real(c_double), parameter :: pi = acos(-1.0d0)

n = 18
n2 = 18
allocate(at(n2), na(n2), nb(n2), nc(n2))
allocate(geo(3, n2))
allocate(xyz(3, n))

! Atomic numbers: C = 6, H = 1
at = [6, 6, 6, 1, 1, 6, 1, 1, 6, 1, 1, 6, 1, 1, 1, 1, 1, 1]

! Internal coordinate references (Fortran is 1-based):
! the zeros are ignored here
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
  geo(2,i) = geo(2,i) * (pi / 180.0d0)
  geo(3,i) = geo(3,i) * (pi / 180.0d0)
  ! works also
  ! geo(2,i) = geo(2,i) * (acos(-1.0d0) / 180.0d0)
  ! geo(3,i) = geo(3,i) * (acos(-1.0d0) / 180.0d0)
end do

call zmat2cart(n,n2,at,geo,xyz,na,nb,nc,fail)

print *, 'fail =', fail
do i = 1, n
   write(istr, '(I0)') i  ! convert integer to string
   print '(A,3F12.6)', 'Atom '//trim(istr)//': ', xyz(1,i), xyz(2,i), xyz(3,i)
end do

call write_xyz("test.xyz", n, at, xyz)

deallocate(at, na, nb, nc)
deallocate(geo, xyz)

end program test_zmat2cart

