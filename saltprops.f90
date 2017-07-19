module saltprops
implicit none
character(*),parameter :: modName = 'saltprops'
integer :: iout = 0

contains 
! TO-DO: fix the units here! Also, work on multiplying by number of moles
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
      write(0,'(4a)') 'FATAL -- ', modName, ' -- ', myName
      write(0,'(a)') 'Temperature out of bounds. Limit 298 < T < 2000'
      write(0,'(a,e12.6)') 'T = ', T
      stop
    endif
    h_nacl = calc_h(T,a,b,c,dfH0)
  endfunction h_nacl
  
  real(8) function h_ucl3(T)
    real(8) :: T
    real(8) :: a,b,c,dfH0
    
    a = 150.0d0
    b = 0.0d0
    c = 0.0d0
    dfH0 = -846.433
    h_ucl3 = calc_h(T,a,b,c,dfH0)
  endfunction h_ucl3
  
  real(8) function h_pucl3(T)
    real(8) :: T
    real(8) :: a,b,c,dfH0
    
    a = 144.0d0
    b = 0.0d0
    c = 0.0d0
    dfH0 = -931.116d0
    h_pucl3 = calc_h(T,a,b,c,dfH0)
  endfunction h_pucl3

  real(8) function calc_h(T,a,b,c,dfH0)
    real(8) :: T,a,b,c,dfH0
    real(8) :: temp_dep, temp_indep
    temp_dep   = a * T + 0.5d0 * b * T ** 2.0d0 + (1.0d0 / 3.0d0) * c * T ** 3.0d0
    temp_indep = a * 298.0d0 + 0.5d0 * 298.0d0 ** 2.0d0 + (1.0d0 / 3.0d0) * c * 298.0d0 ** 3.0d0
    calc_h = dfH0 + temp_dep - temp_indep
  endfunction
endmodule saltprops