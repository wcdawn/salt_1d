module exception_handler
use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
implicit none
character(500),dimension(10) :: msg

private
public :: raise_fatal
public :: raise_warning
public :: msg

contains

  subroutine raise_fatal(modName,myName,msg)
    character(*),intent(in) :: modName,myName
    character(*),dimension(:),intent(in) :: msg
    integer :: i
    
    write(error_unit,*)
    write(error_unit,'(a)') '********************************************************************************'
    write(error_unit,'(4a)') 'FATAL -- ',modName,' -- ',myName
    do i = 1,size(msg)
      if (trim(msg(i)) == '') then
        exit
      endif
      write(error_unit,'(a)') trim(msg(i))
    enddo
    write(error_unit,'(a)') '********************************************************************************'
    write(error_unit,*)
    stop
  endsubroutine raise_fatal
  
  subroutine raise_warning(modName,myName,msg)
    character(*),intent(in) :: modName,myName
    character(*),dimension(:),intent(in) :: msg
    integer :: i
    
    write(output_unit,*)
    write(output_unit,'(4a)') 'warn -- ',modName,' -- ',myName
    do i = 1,size(msg)
      if (trim(msg(i)) == '') then
        exit
      endif
      write(output_unit,'(a)') trim(msg(i))
    enddo
    write(output_unit,*)
    endsubroutine raise_warning
    
endmodule exception_handler