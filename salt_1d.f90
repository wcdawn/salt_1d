program salt_1d
use saltprops
use func_tools
implicit none

real(8) :: tlo,thi
real(8) :: T,dt
integer :: i,nStep
real(8) :: temp
real(8),parameter :: tol = 1.0d-10

tlo = 470.0d0 + 273.15d0 ! [K]
thi = 650.0d0 + 273.15d0 ! [K]
nStep = 1000
dt = (thi - tlo) / real((nStep - 1),8)
T = tlo
do i = 1,nStep
  temp = calc_T(calc_h(T))
  ! write(*,'(3(e12.6,x),l2)') T, temp, (T - temp), (T .approxeq. temp)
  if (T .approxne. temp) then
    write(*,'(a)') 'T /= calc_T(calc_h(T))'
    write(*,'(a,e12.6)') 'T = ',T
    write(*,'(a,e12.6)') 'calc_T(calc_h(T)) = ', temp
    stop
  endif
  T = T + dt
enddo
write(*,'(a)') 'success'

! open(unit = 11, file = 'h_of_T.out',status = 'replace',action = 'write',iostat = ios)
! if (ios /= 0) then
!   write(*,'(a)') 'error opening h_of_T.out'
!   write(*,'(a,i6)') 'ios = ',ios
! endif


endprogram salt_1d
