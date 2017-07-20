module saltprops
use exception_handler
implicit none

! private
character(*),parameter :: modName = 'saltprops'

contains 
! TO-DO: fix the units 
!        work on multiplying by number of moles
!        add abstraction to allow for different functions later
!        add reverse function calc_T
!        mixing
  real(8) function h_nacl(T)
    character(*),parameter :: myName = 'h_nacl'
    real(8) :: T
    real(8) :: a,b,c,dfH0
    
    if ((T >= 298.0d0) .and. (T <= 1500.0d0)) then
      a = 77.7638d0
      b = -0.0075312d0
      c = 0.0d0
      dfH0 = -394.956d0
    elseif ((T > 1500.0d0) .and. (T <= 2000.0d0)) then
      a = 66.944d0
      b = 0.0d0
      c = 0.0d0
      dfH0 = -390.090d0
    else
      msg = ''
      write(msg(1),'(a)') 'Temperature out of bounds. Limit 298 < T < 2000'
      write(msg(2),'(a,e12.6)') 'T = ', T
      call raise_fatal(modName,myName,msg)
    endif
    h_nacl = h_abc(T,a,b,c,dfH0)
  endfunction h_nacl
  
  real(8) function h_ucl3(T)
    real(8) :: T
    real(8) :: a,b,c,dfH0
    
    a = 150.0d0
    b = 0.0d0
    c = 0.0d0
    dfH0 = -846.433
    h_ucl3 = h_abc(T,a,b,c,dfH0)
  endfunction h_ucl3
  
  real(8) function h_pucl3(T)
    real(8) :: T
    real(8) :: a,b,c,dfH0
    
    a = 144.0d0
    b = 0.0d0
    c = 0.0d0
    dfH0 = -931.116d0
    h_pucl3 = h_abc(T,a,b,c,dfH0)
  endfunction h_pucl3

  real(8) function h_abc(T,a,b,c,dfH0)
    real(8) :: T,a,b,c,dfH0
    real(8) :: temp_dep, temp_indep
    
    temp_dep   = a * T + 0.5d0 * b * T ** 2.0d0 + (1.0d0 / 3.0d0) * c * T ** 3.0d0
    temp_indep = a * 298.0d0 + 0.5d0 * b * 298.0d0 ** 2.0d0 + (1.0d0 / 3.0d0) * c * 298.0d0 ** 3.0d0
    h_abc = dfH0 + temp_dep - temp_indep
  endfunction h_abc
  
  real(8) function hmix_ucl3_nacl(x_ucl3)
    real(8) :: x_ucl3
    real(8),dimension(23),parameter :: x_arr = [0.0048d0,0.0505d0,0.0978d0, &
      0.100d0,0.124d0,0.150d0,0.174d0,0.201d0,0.247d0,0.296d0,0.323d0,0.374d0, &
      0.390d0,0.413d0,0.500d0,0.554d0,0.626d0,0.652d0,0.749d0,0.802d0,0.864d0, &
      0.903d0,0.949d0]
    real(8),dimension(23),parameter :: hmix_arr = (-1.0d0) * [0.172d0,1.67d0, &
      2.79d0,3.04d0,3.76d0,3.95d0,4.62d0,5.14d0,5.78d0,5.99d0,6.60d0,7.45d0, &
      7.28d0,7.08d0,7.38d0,6.99d0,7.10d0,6.58d0,4.57d0,3.89d0,1.75d0,1.64d0, &
      0.841d0]
    integer :: ilo
    real(8) :: weight
    
    call return_weight_binary(x_ucl3,x_arr,ilo,weight)
    hmix_ucl3_nacl = linear_interpolate(hmix_arr,ilo,weight)
  endfunction hmix_ucl3_nacl
  
  real(8) function calc_h(T,u_nacl,u_ucl3,u_pucl3)
    character(*),parameter :: myName = 'calc_h'
    real(8) :: T
    real(8) :: x_nacl,x_ucl3,x_pucl3
    real(8),optional :: u_nacl,u_ucl3,u_pucl3
    
    x_nacl = (10.d0 / 19.0d0)
    x_ucl3 = (8.0d0 / 19.0d0)
    x_pucl3 = 1.0d0 - x_nacl - x_ucl3
    if (present(u_nacl)) then
      if ((.not. present(u_pucl3)) .or. (.not. present(u_nacl))) then
        msg = ''
        write(msg(1),'(a)') 'if one user specified mole fraction is present'
        write(msg(2),'(a)') 'then, all must be present'
        call raise_fatal(modName,myName,msg)
      elseif ((u_nacl + u_ucl3 + u_pucl3) /= 1.0d0) then
        msg = ''
        write(msg(1),'(a)') 'sum of user specified mole fraction /= 1.0d0'
        write(msg(2),'(a,e12.6)') 'sum = ',(u_ucl3 + u_pucl3 + u_nacl)
        call raise_fatal(modName,myName,msg)
      endif
      x_nacl = u_nacl
      x_ucl3 = u_ucl3
      x_pucl3 = u_pucl3
    endif
    calc_h = x_nacl * h_nacl(T) + x_ucl3 * h_ucl3(T) + x_pucl3 * h_pucl3(T) + &
      hmix_ucl3_nacl(x_ucl3)
  endfunction calc_h
  
  subroutine return_weight_binary(x,arr,ilo,weight)
    character(*),parameter :: myName = 'return_weight'
    real(8),intent(in) :: x
    real(8),dimension(:),intent(in) :: arr
    integer,intent(out) :: ilo
    real(8),intent(out) :: weight
    
    integer :: lower,upper,i
    logical :: found = .false.
    
    lower = 1
    upper = size(arr)
    do while (.not. found)
      i = (upper + lower) / 2
      if ((x > arr(i)) .and. (x < arr(upper))) then
        upper = upper
        lower = i
      elseif (x < arr(i) .and. (x > arr(lower))) then
        upper = i
        lower = lower
      else
        msg = ''
        write(msg(1),'(a)') 'trouble searching'
        call raise_fatal(modName,myName,msg)
      endif
      if ((x > arr(i)) .and. (x < arr(i + 1))) then
        found = .true.
        ilo = i
        weight = (x - arr(i)) / (arr(i + 1) - arr(i))
        if ((weight > 1.0d0) .or. (weight < 0.0d0)) then
          msg = ''
          write(msg(1),'(a)') 'bad weight'
          write(msg(2),'(a,e12.6)') 'weight = ',weight
          call raise_fatal(modName,myName,msg)
        endif
      endif
    enddo
  endsubroutine return_weight_binary
  
  real(8) function linear_interpolate(arr,ilo,weight)
    real(8),dimension(:) :: arr
    integer :: ilo
    real(8) :: weight
    
    linear_interpolate = arr(ilo) + weight * (arr(ilo + 1) - arr(ilo))
  endfunction linear_interpolate
  
endmodule saltprops