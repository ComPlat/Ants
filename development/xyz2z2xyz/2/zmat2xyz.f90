program zmat2xyz
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer :: i, nat, ios
  character(len=200) :: line, infile, outfile
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)

  character(len=2), allocatable :: symbols(:)
  integer, allocatable :: na(:), nb(:), nc(:)
  real(dp), allocatable :: bond(:), angle(:), dihed(:)
  real(dp), allocatable :: geo(:,:), coord(:,:)
  logical :: fail

  ! === Eingabedateiname holen ===
  call get_command_argument(1, infile)
  if (infile == "") then
    print *, "Usage: ./zmat2xyz <input.zmat>"
    stop
  end if

  ! === Anzahl Zeilen bestimmen ===
  open(10, file=trim(infile), status='old')
  nat = 0
  do
    read(10, '(A)', iostat=ios)
    if (ios /= 0) exit
    nat = nat + 1
  end do
  rewind(10)

  ! === Speicher allokieren ===
  allocate(symbols(nat), na(nat), nb(nat), nc(nat))
  allocate(bond(nat), angle(nat), dihed(nat))
  allocate(geo(3, nat), coord(3, nat))

  ! === Z-Matrix einlesen ===
  do i = 1, nat
    read(10, '(A)') line
    if (i == 1) then
      read(line, *) symbols(i)
      na(i)=0; nb(i)=0; nc(i)=0
      bond(i)=0.0_dp; angle(i)=0.0_dp; dihed(i)=0.0_dp
    else if (i == 2) then
      read(line, *) symbols(i), na(i), bond(i)
      nb(i)=0; nc(i)=0; angle(i)=0.0_dp; dihed(i)=0.0_dp
    else if (i == 3) then
      read(line, *) symbols(i), na(i), bond(i), nb(i), angle(i)
      nc(i)=0; dihed(i)=0.0_dp
    else
      read(line, *) symbols(i), na(i), bond(i), nb(i), angle(i), nc(i), dihed(i)
    end if
  end do
  close(10)

  ! === Umwandlung in Radiant und 0-basierte Indizes
  do i = 1, nat
    if (na(i) > 0) na(i) = na(i) - 1
    if (nb(i) > 0) nb(i) = nb(i) - 1
    if (nc(i) > 0) nc(i) = nc(i) - 1
    angle(i) = angle(i) * pi / 180.0_dp
    dihed(i) = dihed(i) * pi / 180.0_dp
  end do

  geo(1,:) = bond(:)
  geo(2,:) = angle(:)
  geo(3,:) = dihed(:)

  ! === Koordinaten berechnen
  call gmetry(nat, geo, na, nb, nc, coord, fail)

  if (fail) then
    print *, "Fehler bei Rücktransformation"
    stop
  end if

  ! === Ausgabe schreiben
  outfile = infile(1:len_trim(infile)-5) // "_reconstructed.xyz"
  open(20, file=trim(outfile), status='replace')
  write(20,*) nat
  write(20,'(A)') "Reconstructed from Z-Matrix using xtb-style gmetry"
  do i = 1, nat
    write(20,'(A2,3F12.6)') symbols(i), coord(1,i), coord(2,i), coord(3,i)
  end do
  close(20)
  print *, "Wrote reconstructed XYZ to ", trim(outfile)

contains

  subroutine gmetry(natoms, geo, na, nb, nc, coord, fail)
    integer, intent(in) :: natoms
    real(dp), intent(in) :: geo(3, natoms)
    integer, intent(in) :: na(natoms), nb(natoms), nc(natoms)
    real(dp), intent(out) :: coord(3, natoms)
    logical, intent(out) :: fail

    integer :: i, ma, mb, mc, k
    real(dp) :: cosa, sina, sind, cosd, xd, yd, zd
    real(dp) :: xb, yb, zb, xa, ya, za
    real(dp) :: xpa, ypa, xqa, zqa, xpb
    real(dp) :: ypd, zpd, xpd, zqd, xqd, yqd
    real(dp) :: rbc, rbc2, xyb, yza
    real(dp) :: costh, sinth, cosph, sinph, coskh, sinkh

    fail = .false.
    coord(:,1) = 0.0_dp
    coord(1,2) = geo(1,2)
    coord(2,2) = 0.0_dp
    coord(3,2) = 0.0_dp

    if (natoms == 2) return

    cosa = cos(geo(2,3))
    if (na(3) == 0) then
      coord(1,3) = coord(1,1) + geo(1,3) * cosa
    else
      coord(1,3) = coord(1,2) - geo(1,3) * cosa
    end if
    coord(2,3) = geo(1,3) * sin(geo(2,3))
    coord(3,3) = 0.0_dp

    do i = 4, natoms
      cosa = cos(geo(2,i))
      mb = nb(i)
      mc = na(i)

      xb = coord(1,mb+1) - coord(1,mc+1)
      yb = coord(2,mb+1) - coord(2,mc+1)
      zb = coord(3,mb+1) - coord(3,mc+1)

      rbc = 1.0_dp / sqrt(xb**2 + yb**2 + zb**2)

      if (abs(cosa) >= 0.9999999999_dp) then
        rbc2 = geo(1,i) * rbc * cosa
        coord(1,i) = coord(1,mc+1) + xb * rbc2
        coord(2,i) = coord(2,mc+1) + yb * rbc2
        coord(3,i) = coord(3,mc+1) + zb * rbc2
        cycle
      end if

      ma = nc(i)

      xa = coord(1,ma+1) - coord(1,mc+1)
      ya = coord(2,ma+1) - coord(2,mc+1)
      za = coord(3,ma+1) - coord(3,mc+1)

      xyb = sqrt(xb**2 + yb**2)
      k = -1
      if (xyb <= 0.1_dp) then
        xpa = za
        za  = -xa
        xa  = xpa
        xpb = zb
        zb  = -xb
        xb  = xpb
        xyb = sqrt(xb**2 + yb**2)
        k = 1
      end if

      costh = xb / xyb
      sinth = yb / xyb

      xpa = xa * costh + ya * sinth
      ypa = ya * costh - xa * sinth

      sinph = zb * rbc
      cosph = sqrt(abs(1.0_dp - sinph**2))

      xqa = xpa * cosph + za * sinph
      zqa = za * cosph - xpa * sinph

      yza = sqrt(ypa**2 + zqa**2)
      if (yza < 1.0e-4_dp) then
        fail = .true.
        return
      end if

      coskh = ypa / yza
      sinkh = zqa / yza

      sina = sin(geo(2,i))
      sind = -sin(geo(3,i))
      cosd = cos(geo(3,i))

      xd = geo(1,i) * cosa
      yd = geo(1,i) * sina * cosd
      zd = geo(1,i) * sina * sind

      ypd = yd * coskh - zd * sinkh
      zpd = zd * coskh + yd * sinkh
      xpd = xd * cosph - zpd * sinph
      zqd = zpd * cosph + xd * sinph
      xqd = xpd * costh - ypd * sinth
      yqd = ypd * costh + xpd * sinth

      if (k < 1) then
        coord(1,i) = xqd + coord(1,mc+1)
        coord(2,i) = yqd + coord(2,mc+1)
        coord(3,i) = zqd + coord(3,mc+1)
      else
        coord(1,i) = -zqd + coord(1,mc+1)
        coord(2,i) = yqd + coord(2,mc+1)
        coord(3,i) = xqd + coord(3,mc+1)
      end if
    end do
  end subroutine gmetry

end program zmat2xyz
