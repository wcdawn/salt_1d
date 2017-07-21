program salt_1d
use saltprops
implicit none

real(8) :: tlo,thi
real(8) :: T,dt
integer :: i,nStep

write(*,*) calc_h(500.0d0)
write(*,*) calc_h(600.0d0)
write(*,*) calc_h(700.0d0)
write(*,*) calc_h(800.0d0)
write(*,*) calc_h(900.0d0)
tlo = 470.0d0 + 273.15d0 ! [K]
thi = 650.0d0 + 273.15d0 ! [K]
nStep = 100
dt = (thi - tlo) / real((nStep - 1),8)
T = tlo
do i = 1,nStep
  write(10,'(2e13.6)')  T, calc_h(T)
  T = T + dt
enddo

endprogram salt_1d