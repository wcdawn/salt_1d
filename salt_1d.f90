program salt_1d
use saltprops
implicit none

! real(8) :: tlo,thi
! real(8) :: T,dt
! integer :: i,nStep
real(8) :: alpha,beta,gamma,delta
real(8) :: h
! tlo = 470.0d0 + 273.15d0 ! [K]
! thi = 650.0d0 + 273.15d0 ! [K]
! nStep = 100
! dt = (thi - tlo) / real((nStep - 1),8)
! T = tlo
! do i = 1,nStep
  ! write(10,'(2e13.6)')  T, calc_h(T)
  ! T = T + dt
! enddo

h = 200.0d0
write(*,*) calc_T(h)
alpha = 0.0d0
beta  = 2.0d0
gamma = 4.0d0
delta = -4.0d0
write(*,*) cubic_solve(alpha,beta,gamma,delta,RealPos=.true.)

endprogram salt_1d