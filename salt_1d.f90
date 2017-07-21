program salt_1d
use saltprops
implicit none

integer :: ios
real(8) :: tlo,thi
real(8) :: T,dt
integer :: i,nStep
real(8) :: alpha,beta,gamma,delta
real(8) :: h
! tlo = 470.0d0 + 273.15d0 ! [K]
! thi = 650.0d0 + 273.15d0 ! [K]
open(unit = 11, file = 'h_of_T.out',status = 'replace',action = 'write',iostat = ios)
if (ios /= 0) then
  write(*,'(a)') 'error opening h_of_T.out'
  write(*,'(a,i6)') 'ios = ',ios
endif
tlo = 1.0d0
thi = 1.0d6
nStep = 1000
dt = (thi - tlo) / real((nStep - 1),8)
T = tlo
do i = 1,nStep
  write(11,'(2e13.6)')  T, calc_h(T)
  T = T + dt
enddo

h = 269.0d0
write(*,*) calc_T(h)
! should return 810
alpha = 0.0d0
beta  = 2.0d0
gamma = 4.0d0
delta = -4.0d0
write(*,*) cubic_solve(alpha,beta,gamma,delta,umax=2.0d3,umin=0.0d0)

endprogram salt_1d