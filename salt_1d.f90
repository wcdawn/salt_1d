program salt_1d
use saltprops
implicit none

! real(8) :: tlo,thi
! real(8) :: T,dt
! integer :: i,nStep
real(8) :: alpha,beta,gamma,delta

! tlo = 470.0d0 + 273.15d0 ! [K]
! thi = 650.0d0 + 273.15d0 ! [K]
! nStep = 100
! dt = (thi - tlo) / real((nStep - 1),8)
! T = tlo
! do i = 1,nStep
  ! write(10,'(2e13.6)')  T, calc_h(T)
  ! T = T + dt
! enddo

write(*,*) calc_T(200.0d0)
alpha = 1.0d0
beta  = -3.0d0
gamma = -144.0d0
delta = 432.0d0
write(*,*) cubic_solve(alpha,beta,gamma,delta)

endprogram salt_1d