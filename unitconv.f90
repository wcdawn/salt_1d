module unitconv
implicit none

!-------------------------------------------------------------------------------
! conversion factors and procedures to convert between SI and U.S. units
! e.g.
! m  --> ft : multiply with unit%m_ft
! ft --> m  : multiply with unit%ft_m
!-------------------------------------------------------------------------------
private
! public :: unit

type :: unit_type
  sequence
  real(8) :: ft_in=12.0d0
  ! real(8) :: ft_in = 12.0d0
  ! real(8) :: in_ft = 1.0d0 / ft_in
  ! real(8) :: in_cm = 2.54d0
  ! real(8) :: in_m  = in_cm * 1.0d-2
  ! real(8) :: ft_m  = ft_in * ft_in_m
  ! real(8) :: cm_in = 1.0d0 / in_cm
  ! real(8) :: m_in  = 1.0d0 / in_m
  ! real(8) :: m_ft  = 1.0d0 / ft_m
endtype unit_type

type(unit_type) :: unit

endmodule unitconv
