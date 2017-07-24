module func_tools
implicit none

!-------------------------------------------------------------------------------
! module for performing basic functions
! includes custom operators
!-------------------------------------------------------------------------------

private
public :: operator(.approxeq.)
public :: operator(.approxne.)
public :: epsd
character(*),parameter :: modName = 'func_tools'
real(8),parameter :: epsd = 1.0d-10

! define interfaces for custom operators
interface operator(.approxeq.)
  module procedure approxeq_double
endinterface
interface operator(.approxne.)
  module procedure approxne_double
endinterface

contains

!-------------------------------------------------------------------------------
! approximately equal comparison based on public, global, epsilon
!
! arguments --------------------------------------------------------------------
! a     real(8) input for comparison
! b     real(8) input for comparison
! bool  boolean value for return
! local variables --------------------------------------------------------------
!-------------------------------------------------------------------------------
  function approxeq_double(a,b) result(bool)
    character(*),parameter :: myName = 'approxeq_double'
    real(8),intent(in) :: a,b
    logical :: bool
    bool = (abs(a - b) < epsd)
  endfunction approxeq_double

!-------------------------------------------------------------------------------
! approximately not equal comparison based on public, global, epsilon
!
! arguments --------------------------------------------------------------------
! a     real(8) input for comparison
! b     real(8) input for comparison
! bool  boolean value for return
! local variables --------------------------------------------------------------
!-------------------------------------------------------------------------------
  function approxne_double(a,b) result(bool)
    real(8),intent(in) :: a,b
    logical :: bool
    bool = (abs(a - b) > epsd)
  endfunction approxne_double

endmodule func_tools
