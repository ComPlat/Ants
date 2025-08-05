
program xyz2zmat
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: max_atoms = 200
  character(len=2) :: symbols(max_atoms)
  integer :: nat, i, j, na(max_atoms), nb(max_atoms), nc(max_atoms)
  real(dp) :: xyz(3, max_atoms), geo(3, max_atoms)
  logical :: bonded(max_atoms, max_atoms)
  character(len=200) :: line
  character(len=100) :: infile, outline

  real(8), dimension(118) :: covalent_radii

  call init_radii()
  call read_xyz("input.xyz", nat, symbols, xyz)
  call get_bonded_atoms(nat, symbols, xyz, bonded, covalent_radii)
  call determine_references(nat, bonded, na, nb, nc)
  call compute_geometry(nat, xyz, na, nb, nc, geo)
  call write_zmatrix("output.zmat", nat, symbols, na, nb, nc, geo)

contains

  subroutine init_radii()
    covalent_radii = 0.0d0
    covalent_radii(1)  = 0.31d0
    covalent_radii(6)  = 0.76d0
    covalent_radii(7)  = 0.71d0
    covalent_radii(8)  = 0.66d0
    covalent_radii(9)  = 0.57d0
    covalent_radii(16) = 1.05d0
    covalent_radii(17) = 1.02d0
    covalent_radii(26) = 1.32d0
    covalent_radii(28) = 1.24d0
  end subroutine

    subroutine read_xyz(file, nat, symbols, xyz)
      character(len=*), intent(in) :: file
      integer, intent(out) :: nat
      character(len=2), intent(out) :: symbols(max_atoms)
      real(dp), intent(out) :: xyz(3, max_atoms)
      integer :: i, ios
      character(len=200) :: line
      open(unit=10, file=file, status='old')
      read(10,*) nat
      read(10,*)
      do i = 1, nat
        read(10,'(A)', iostat=ios) line
        if (ios /= 0) exit
        read(line,*) symbols(i), xyz(1,i), xyz(2,i), xyz(3,i)
      end do
      close(10)
    end subroutine


  subroutine get_bonded_atoms(nat, symbols, xyz, bonded, radii)
    integer, intent(in) :: nat
    character(len=2), intent(in) :: symbols(nat)
    real(dp), intent(in) :: xyz(3, nat)
    logical, intent(out) :: bonded(nat, nat)
    real(dp), intent(in) :: radii(:)
    integer :: i, j, Zi, Zj
    real(dp) :: dist, r_cov, tol
    tol = 0.45d0
    bonded = .false.
    do i = 1, nat
      Zi = get_Z(symbols(i))
      do j = i+1, nat
        Zj = get_Z(symbols(j))
        dist = norm2(xyz(:,i) - xyz(:,j))
        r_cov = radii(Zi) + radii(Zj) + tol
        if (dist <= r_cov) then
          bonded(i,j) = .true.
          bonded(j,i) = .true.
        end if
      end do
    end do
  end subroutine

  subroutine determine_references(nat, bonded, na, nb, nc)
    integer, intent(in) :: nat
    logical, intent(in) :: bonded(nat,nat)
    integer, intent(out) :: na(nat), nb(nat), nc(nat)
    integer :: i, j
    na = 0; nb = 0; nc = 0
    do i = 2, nat
      do j = 1, i-1
        if (bonded(i,j)) then
          na(i) = j
          exit
        end if
      end do
      if (na(i) == 0) na(i) = 1
      if (i >= 3) then
        do j = 1, i-1
          if (j /= na(i) .and. bonded(na(i),j)) then
            nb(i) = j
            exit
          end if
        end do
        if (nb(i) == 0) nb(i) = merge(1, 2, na(i)==1)
      end if
      if (i >= 4) then
        do j = 1, i-1
          if (j /= na(i) .and. j /= nb(i) .and. bonded(nb(i),j)) then
            nc(i) = j
            exit
          end if
        end do
        if (nc(i) == 0) nc(i) = 1
      end if
    end do
  end subroutine

  subroutine compute_geometry(nat, xyz, na, nb, nc, geo)
    integer, intent(in) :: nat, na(nat), nb(nat), nc(nat)
    real(dp), intent(in) :: xyz(3,nat)
    real(dp), intent(out) :: geo(3,nat)
    integer :: i
    geo = 0.0d0
    do i = 1, nat
      if (na(i) > 0) geo(1,i) = norm2(xyz(:,i) - xyz(:,na(i)))
      if (nb(i) > 0) geo(2,i) = angle(xyz(:,i), xyz(:,na(i)), xyz(:,nb(i)))
      if (nc(i) > 0) geo(3,i) = dihedral(xyz(:,nc(i)), xyz(:,nb(i)), xyz(:,na(i)), xyz(:,i))
    end do
  end subroutine

  subroutine write_zmatrix(file, nat, symbols, na, nb, nc, geo)
    character(len=*), intent(in) :: file
    integer, intent(in) :: nat, na(nat), nb(nat), nc(nat)
    character(len=2), intent(in) :: symbols(nat)
    real(dp), intent(in) :: geo(3,nat)
    integer :: i
    open(unit=11, file=file, status='replace')
    do i = 1, nat
      select case(i)
      case(1)
        write(11,'(A2)') symbols(i)
      case(2)
        write(11,'(A2,1X,I3,1X,F12.6)') symbols(i), na(i), geo(1,i)
      case(3)
        write(11,'(A2,1X,I3,1X,F12.6,1X,I3,1X,F12.6)') symbols(i), na(i), geo(1,i), nb(i), geo(2,i)*180.0d0/acos(-1.0d0)
      case default
        write(11,'(A2,1X,I3,1X,F12.6,1X,I3,1X,F12.6,1X,I3,1X,F12.6)') &
              symbols(i), na(i), geo(1,i), nb(i), geo(2,i)*180.0d0/acos(-1.0d0), nc(i), geo(3,i)*180.0d0/acos(-1.0d0)
      end select
    end do
    close(11)
  end subroutine

  integer function get_Z(sym)
    character(len=2), intent(in) :: sym
    select case (trim(sym))
    case ("H"); get_Z = 1
    case ("C"); get_Z = 6
    case ("N"); get_Z = 7
    case ("O"); get_Z = 8
    case ("F"); get_Z = 9
    case ("S"); get_Z = 16
    case ("Cl"); get_Z = 17
    case ("Fe"); get_Z = 26
    case ("Ni"); get_Z = 28
    case default; get_Z = 1
    end select
  end function

  real(dp) function norm2(v)
    real(dp), intent(in) :: v(3)
    norm2 = sqrt(sum(v**2))
  end function

  real(dp) function angle(a, b, c)
    real(dp), intent(in) :: a(3), b(3), c(3)
    real(dp) :: ba(3), bc(3), cosa
    ba = a - b
    bc = c - b
    cosa = dot_product(ba, bc)/(norm2(ba)*norm2(bc))
    cosa = max(min(cosa,1.0d0), -1.0d0)
    angle = acos(cosa)
  end function

  real(dp) function dihedral(p0, p1, p2, p3)
    real(dp), intent(in) :: p0(3), p1(3), p2(3), p3(3)
    real(dp) :: b0(3), b1(3), b2(3), v(3), w(3), x, y
    b0 = p0 - p1
    b1 = p2 - p1
    b2 = p3 - p2
    b1 = b1 / norm2(b1)
    v = b0 - dot_product(b0,b1)*b1
    w = b2 - dot_product(b2,b1)*b1
    x = dot_product(v, w)
    y = dot_product(cross_product(b1, v), w)
    dihedral = atan2(y, x)
  end function

  function cross_product(a, b) result(c)
    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: c(3)
    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
  end function

end program xyz2zmat
