module saltprops
use exception_handler
implicit none

!-------------------------------------------------------------------------------
! module for calculating salt properties
! currently designed with 1*PuCl3-8*UCl3-10*NaCl in mind and calculating only 
! enthalpy and temperature
! the functions were maintained general where practical and should be expanded 
! in the future
!------------------------------------------------------------------------------

private
public :: calc_h
character(*),parameter :: modName = 'saltprops'
contains 
! TO-DO:
!        add abstraction to allow for different functions later
!        add reverse function calc_T

!-------------------------------------------------------------------------------
! calculate molar enthalpy of NaCl
! coefficients from "Thermodynamic evaluation of the NaCl-MgCl2-Ucl3-PuCl3 system"
!   O. Benes, R.J.M. Konings, 2008
!
! arguments --------------------------------------------------------------------
! T       temperature             [K]
! h_nacl  molar enthalpy of NaCl  [kJ/mol]
! local variables --------------------------------------------------------------
! a     coefficient for calculation of cp = a + b*T + c*T**2  [kJ/K/mol]
! b     coefficient for calculation of cp = a + b*T + c*T**2  [kJ/K/mol]
! c     coefficient for calculation of cp = a + b*T + c*T**2  [kJ/K/mol]
!-------------------------------------------------------------------------------
  real(8) function h_nacl(T)
    character(*),parameter :: myName = 'h_nacl'
    real(8) :: T
    real(8) :: a,b,c
    
    if ((T >= 298.0d0) .and. (T <= 1500.0d0)) then
      a = 77.7638d-3
      b = -0.0075312d-3
      c = 0.0d-3
    elseif ((T > 1500.0d0) .and. (T <= 2000.0d0)) then
      a = 66.944d-3
      b = 0.0d-3
      c = 0.0d-3
    else
      msg = ''
      write(msg(1),'(a)') 'Temperature out of bounds. Limit 298 < T < 2000'
      write(msg(2),'(a,e12.6)') 'T = ', T
      call raise_fatal(modName,myName,msg)
    endif
    h_nacl = h_abc(T,a,b,c)
  endfunction h_nacl
  
!-------------------------------------------------------------------------------
! calculate molar enthalpy of UCl3
! coefficients from "Thermodynamic evaluation of the NaCl-MgCl2-Ucl3-PuCl3 system"
!   O. Benes, R.J.M. Konings, 2008
! these coefficients are noted as "Estimated"
! 
! arguments --------------------------------------------------------------------
! T       temperature             [K]
! h_ucl3  molar enthalpy of UCl3  [kJ/mol]
! local variables --------------------------------------------------------------
! a     coefficient for calculation of cp = a + b*T + c*T**2  [kJ/K/mol]
! b     coefficient for calculation of cp = a + b*T + c*T**2  [kJ/K/mol]
! c     coefficient for calculation of cp = a + b*T + c*T**2  [kJ/K/mol]
!-------------------------------------------------------------------------------
  real(8) function h_ucl3(T)
    real(8) :: T
    real(8) :: a,b,c
    
    a = 150.0d-3
    b = 0.0d-3
    c = 0.0d-3
    h_ucl3 = h_abc(T,a,b,c)
  endfunction h_ucl3
  
!-------------------------------------------------------------------------------
! calculate molar enthalpy of PuCl3
! coefficients from "Thermodynamic evaluation of the NaCl-MgCl2-Ucl3-PuCl3 system"
!   O. Benes, R.J.M. Konings, 2008
!
! arguments --------------------------------------------------------------------
! T        temperature              [K]
! h_pucl3  molar enthalpy of PuCl3  [kJ/mol]
! local variables --------------------------------------------------------------
! a     coefficient for calculation of cp = a + b*T + c*T**2  [kJ/K/mol]
! b     coefficient for calculation of cp = a + b*T + c*T**2  [kJ/K/mol]
! c     coefficient for calculation of cp = a + b*T + c*T**2  [kJ/K/mol]
!-------------------------------------------------------------------------------
  real(8) function h_pucl3(T)
    real(8) :: T
    real(8) :: a,b,c
    
    a = 144.0d-3
    b = 0.0d-3
    c = 0.0d-3
    h_pucl3 = h_abc(T,a,b,c)
  endfunction h_pucl3

!-------------------------------------------------------------------------------
! calculate molar enthalpy given coefficients a,b, and c.
! absolute enthalpies referenced to T0=273.15K
! equation of the form h = integral(cp(T),from=273.15K,to=T)
!
! arguments --------------------------------------------------------------------
! T      temperature in Kelvin
! a      coefficient for calculation of cp = a + b*T + c*T**2  [kJ/K/mol]
! b      coefficient for calculation of cp = a + b*T + c*T**2  [kJ/K/mol]
! c      coefficient for calculation of cp = a + b*T + c*T**2  [kJ/K/mol]
! dfH0   standard enthalpy of formation                        [kJ/mol]
! h_abc  molar enthalpy at temperature T                       [kJ/mol]
! local variables --------------------------------------------------------------
! cp_int  integral of cp
! dT      difference between temperature and reference temperature  [K]
!-------------------------------------------------------------------------------
  real(8) function h_abc(T,a,b,c)
    real(8) :: T,a,b,c
    real(8) :: cp_int,dT
    
    dT = T - 273.15d0
    cp_int   = a * dT + 0.5d0 * b * dT ** 2.0d0 + (1.0d0 / 3.0d0) * c * dT ** 3.0d0
    h_abc = cp_int
  endfunction h_abc
  
!-------------------------------------------------------------------------------
! calculate the molar enthalpy of mixing in a molten UCl3-NaCl system
! data from "Enthalpies of mixing in molten UCle-NaCl system"
!   Haruaki Matsurra et al., 2002
! in using this function for a PuCl3-UCl3-NaCl system, PuCl3 mole fraction
! should be included in x_ucl3
!
! arguments --------------------------------------------------------------------
! x_ucl3          mole fraction of UCl3 in the UCl3-NaCl system
! hmix_ucl3_nacl  molar enthalpy of mixing of UCl3-NaCl system   [kJ/mol]
! local variables --------------------------------------------------------------
! x_arr     array of x_ucl3 from experimental data
! hmix_arr  array of molar enthalpy of mixing from experimental data [kJ/mol]
! ilo       index of lower value for interpolation
! weight    weight fraction for interpolation
!-------------------------------------------------------------------------------
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
  
!-------------------------------------------------------------------------------
! calculate molar enthalpy of PuCl3-UCl3-NaCl system
! default mole fractions are assumed fror 1*PuCl3-8*UCl3-10*NaCl as reported in 
! "Reator with very low fission product inventory" M. Taube, W. Heer, 1980
!
! arguments --------------------------------------------------------------------
! T       temperature of system     [K]
! calc_h  specific enthalpy of system  [kJ/kg]
!   u_nacl   user specified mole fraction of NaCl
!   u_ucl3   user specified mole fraction of UCl3
!   u_pucl3  user specified mole fraction of PuCl3
! local variables --------------------------------------------------------------
! x_nacl   mole fraction of NaCl (default=10/19)
! x_ucl3   mole fraction of UCl3 (default=8/19)
! x_pucl3  mole fraction of PuCl3 (default=1/19)
! m_nacl   molar mass of NaCl (58.4428)           [gm/mol]
! m_ucl3   molar mass of UCl3 (344.3879)          [gm/mol]
! m_pucl3  molar mass of PuCl3 (344.4086)         [gm/mol]
! m_mol    molar mass of PuCl3-UCl3-NaCl system   [kg/mol]
!-------------------------------------------------------------------------------
  real(8) function calc_h(T)!,u_nacl,u_ucl3,u_pucl3)
    character(*),parameter :: myName = 'calc_h'
    real(8) :: T
    ! real(8),optional :: u_nacl,u_ucl3,u_pucl3
    real(8) :: x_nacl,x_ucl3,x_pucl3
    real(8) :: m_nacl,m_ucl3,m_pucl3,m_mol
    
    write(20,*) myName,T
    ! calculate molar masses based on elemental consituents
    m_nacl = 22.989769280d0 + 35.4530d0
    m_ucl3 = 238.028910d0 + 3.0d0 * 35.4530d0
    m_pucl3 = 238.0495599d0 + 3.0d0 * 35.4530d0
    ! set default mole fractions
    x_nacl = (10.d0 / 19.0d0)
    x_ucl3 = (8.0d0 / 19.0d0)
    x_pucl3 = 1.0d0 - x_nacl - x_ucl3
    ! if (present(u_nacl)) then
      ! if ((.not. present(u_pucl3)) .or. (.not. present(u_nacl))) then
        ! msg = ''
        ! write(msg(1),'(a)') 'if one user specified mole fraction is' // &
          ! ' present then, all must be present'
        ! call raise_fatal(modName,myName,msg)
      ! elseif ((u_nacl + u_ucl3 + u_pucl3) /= 1.0d0) then
        ! msg = ''
        ! write(msg(1),'(a)') 'sum of user specified mole fraction /= 1.0d0'
        ! write(msg(2),'(a,e12.6)') 'sum = ',(u_ucl3 + u_pucl3 + u_nacl)
        ! call raise_fatal(modName,myName,msg)
      ! endif
      ! overwrite default mole fractions
      ! x_nacl = u_nacl
      ! x_ucl3 = u_ucl3
      ! x_pucl3 = u_pucl3
    ! endif
    ! calculate molar mass of system
    m_mol = (m_nacl * x_nacl + m_ucl3 * x_ucl3 + m_pucl3 * x_pucl3) * 1.0d-3
    ! calculate molar enthalpy of system
    calc_h = x_nacl * h_nacl(T) + x_ucl3 * h_ucl3(T) + x_pucl3 * h_pucl3(T) + &
      hmix_ucl3_nacl(x_ucl3 + x_pucl3)
    ! convert molar enthalpy to specific enthalpy
    calc_h = calc_h * (1.0d0 / m_mol)
  endfunction calc_h
  
!-------------------------------------------------------------------------------
! return the lower index and weight fraction for interpolation
! this method employs a binary search for an ordered list of general length
! TO-DO: consider moving this to another module
!
! arguments --------------------------------------------------------------------
! x       value for which interpolation will be performed
! arr     array of values in-between which x will be interpolated
! ilo     output lower index for interpolation
! weight  output weight fraction for interpolation
! local variables --------------------------------------------------------------
! lower  lower bound for current iteration of binary search
! upper  upper bound for current iteration of binary search
! i      index for current iteration of binary search
! found  logical if interpolation data has been successfully found
!-------------------------------------------------------------------------------
  subroutine return_weight_binary(x,arr,ilo,weight)
    character(*),parameter :: myName = 'return_weight'
    real(8),intent(in) :: x
    real(8),dimension(:),intent(in) :: arr
    integer,intent(out) :: ilo
    real(8),intent(out) :: weight
    integer :: lower,upper,i
    logical :: found = .false.
    
    if ((x > arr(size(arr))) .or. (x < arr(1))) then
      msg = ''
      write(msg(1),'(a)') 'x out of bounds for interpolation'
      write(msg(2),'(a,e12.6)') 'x = ', x
      write(msg(3),'(2(a,e12.6))') 'arr(1) = ',arr(1), &
        'arr(size(arr)) = ',arr(size(arr))
      call raise_fatal(modName,myName,msg)
    endif    
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
  
!-------------------------------------------------------------------------------
! linearly interpolate on an array of data
! TO-DO: consider moving this to another module
!
! arguments --------------------------------------------------------------------
! arr                 array of data for linear interpolation
! ilo                 lower index of array for interpolation
! weight              weight fraction for interpolation
! linear_interpolate  linear interpolation of array between ilo and ilo+1
! local variables --------------------------------------------------------------
!-------------------------------------------------------------------------------
  real(8) function linear_interpolate(arr,ilo,weight)
    real(8),dimension(:) :: arr
    integer :: ilo
    real(8) :: weight
    
    linear_interpolate = arr(ilo) + weight * (arr(ilo + 1) - arr(ilo))
  endfunction linear_interpolate
  
endmodule saltprops