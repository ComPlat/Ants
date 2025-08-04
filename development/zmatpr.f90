subroutine zmatpr(nat, at, geo, na, nb, nc, molnum)
   use mod_symbols, only: toSymbol
   implicit none
   integer, intent(in) :: nat
   integer, intent(in) :: at(nat)
   integer, intent(in) :: na(nat), nb(nat), nc(nat)
   real(8), intent(in) :: geo(3,nat)
   integer, intent(in) :: molnum
   character(len=20) :: filename
   integer :: i, l, m, n, ich
   real(8) :: bl, ang, dihed
   real(8), parameter :: pi = 3.14159265358979323846d0

   do i = 1, nat
      l = 1; m = 1; n = 1
      if (i == 1) then
         l = 0; m = 0; n = 0
      else if (i == 2) then
         m = 0; n = 0
      else if (i == 3) then
         n = 0
      end if
      bl = geo(1,i)
      ang = geo(2,i)*180.0d0/pi
      dihed = geo(3,i)*180.0d0/pi
      if (dihed > 180.0d0) dihed = dihed - 360.0d0
      write(*,'(i4,2x,a2,f12.6,2x,f10.4,2x,f10.4,i6,2i5)') i, toSymbol(at(i)), bl, ang, dihed, na(i), nb(i), nc(i)
   end do

   write(filename,'("zmatrix",i0,".zmat")') molnum
   open(unit=20, file=filename, status='replace')

   write(20,'(a2)') toSymbol(at(1))
   write(20,'(a2,1x,i0,1x,f8.3)') toSymbol(at(2)), na(2), geo(1,2)
   write(20,'(a2,1x,i0,1x,f8.3,1x,i0,1x,f8.3)') toSymbol(at(3)), na(3), geo(1,3), nb(3), geo(2,3)*180.0d0/pi

   do i = 4, nat
      bl = geo(1,i)
      ang = geo(2,i)*180.0d0/pi
      dihed = geo(3,i)*180.0d0/pi
      if (dihed > 180.0d0) dihed = dihed - 360.0d0
      write(20,'(a2,1x,i0,1x,f8.3,1x,i0,1x,f8.3,1x,i0,1x,f8.3)') toSymbol(at(i)), na(i), bl, nb(i), ang, nc(i), dihed
   end do
   write(20,*)
   close(20)
end subroutine zmatpr

