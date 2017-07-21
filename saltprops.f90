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
public :: calc_T
! TO-DO: make this private
public :: cubic_solve
character(*),parameter :: modName = 'saltprops'
character(*),parameter :: f1 = '(a)'
character(*),parameter :: f2 = '(a,e12.6)'
! 101 format(a)       ! plain-text descriptor
! 102 format(a,e12.6) ! plain-text followed by real
contains 
! TO-DO:
!        add abstraction to allow for different functions later
!        add reverse function calc_T
!        add .APPROXEQ. see .../Futility/src/IntrType.f90
!        fix cubic_solve add option for positive real root
!        clarify when h is molar and specific (hm vs. h)

!-------------------------------------------------------------------------------
! return coefficients for calculating cp for NaCl
! coefficients from "Thermodynamic evaluation of the NaCl-MgCl2-Ucl3-PuCl3 system"
!   O. Benes, R.J.M. Konings, 2008
!
! arguments --------------------------------------------------------------------
! T  temperature                                     [K]
! a  coefficient for calculation of cp=a+b*T+c*T**2  [kJ/K/mol]
! b  coefficient for calculation of cp=a+b*T+c*T**2  [kJ/K/mol]
! c  coefficient for calculation of cp=a+b*T+c*T**2  [kJ/K/mol]
! local variables --------------------------------------------------------------
!-------------------------------------------------------------------------------
  subroutine nacl_abc(T,a,b,c)
    character(*),parameter :: myName = 'h_nacl_abc'
    real(8),intent(in) :: T
    real(8),intent(out) :: a,b,c
    
    ! TO-DO: repair
    ! if ((T >= 298.0d0) .and. (T <= 1500.0d0)) then
    if ((T <= 1500.0d0)) then
      a = 77.7638d-3
      b = -0.0075312d-3
      c = 0.0d-3
    ! elseif ((T > 1500.0d0) .and. (T <= 2000.0d0)) then
    elseif ((T > 1500.0d0)) then
      msg = ''
      write(msg(1),f1) 'this set of NaCl properties is not properly ' // &
        'incorporated into the codebase'
      write(msg(2),f1) 'integrating this function will yield a moderatly ' // &
        'incorrect result'
      call raise_warning(modName,myName,msg)
      a = 66.944d-3
      b = 0.0d-3
      c = 0.0d-3
    else
      msg = ''
      write(msg(1),f1) 'Temperature out of bounds. Limit 298 < T < 2000'
      write(msg(2),f2) 'T = ', T
      ! call raise_fatal(modName,myName,msg)
      call raise_warning(modName,myName,msg)
    endif
  endsubroutine nacl_abc

!-------------------------------------------------------------------------------
! calculate molar enthalpy of NaCl
!
! arguments --------------------------------------------------------------------
! T       temperature             [K]
! h_nacl  molar enthalpy of NaCl  [kJ/mol]
! local variables --------------------------------------------------------------
! a     coefficient for calculation of cp=a+b*T+c*T**2  [kJ/K/mol]
! b     coefficient for calculation of cp=a+b*T+c*T**2  [kJ/K/mol]
! c     coefficient for calculation of cp=a+b*T+c*T**2  [kJ/K/mol]
!-------------------------------------------------------------------------------
  real(8) function h_nacl(T)
    real(8) :: T
    real(8) :: a,b,c
    
    call nacl_abc(T,a,b,c)
    h_nacl = h_abc(T,a,b,c)
  endfunction h_nacl

!-------------------------------------------------------------------------------
! return coefficients for calculating cp for UCl3
! coefficients from "Thermodynamic evaluation of the NaCl-MgCl2-Ucl3-PuCl3 system"
!   O. Benes, R.J.M. Konings, 2008
! these coefficients are noted as "Estimated"
!
! arguments --------------------------------------------------------------------
! a  coefficient for calculation of cp=a+b*T+c*T**2  [kJ/K/mol]
! b  coefficient for calculation of cp=a+b*T+c*T**2  [kJ/K/mol]
! c  coefficient for calculation of cp=a+b*T+c*T**2  [kJ/K/mol]
! local variables --------------------------------------------------------------
!-------------------------------------------------------------------------------
  subroutine ucl3_abc(a,b,c)
    real(8),intent(out) :: a,b,c
    
    a = 150.0d-3
    b = 0.0d-3
    c = 0.0d-3
  endsubroutine ucl3_abc
  
!-------------------------------------------------------------------------------
! calculate molar enthalpy of UCl3
! 
! arguments --------------------------------------------------------------------
! T       temperature             [K]
! h_ucl3  molar enthalpy of UCl3  [kJ/mol]
! local variables --------------------------------------------------------------
! a     coefficient for calculation of cp=a+b*T+c*T**2  [kJ/K/mol]
! b     coefficient for calculation of cp=a+b*T+c*T**2  [kJ/K/mol]
! c     coefficient for calculation of cp=a+b*T+c*T**2  [kJ/K/mol]
!-------------------------------------------------------------------------------
  real(8) function h_ucl3(T)
    real(8) :: T
    real(8) :: a,b,c
    
    call ucl3_abc(a,b,c)
    h_ucl3 = h_abc(T,a,b,c)
  endfunction h_ucl3
  
!-------------------------------------------------------------------------------
! return coefficients for calculating cp for PuCl3
! coefficients from "Thermodynamic evaluation of the NaCl-MgCl2-Ucl3-PuCl3 system"
!   O. Benes, R.J.M. Konings, 2008
!
! arguments --------------------------------------------------------------------
! a  coefficient for calculation of cp=a+b*T+c*T**2  [kJ/K/mol]
! b  coefficient for calculation of cp=a+b*T+c*T**2  [kJ/K/mol]
! c  coefficient for calculation of cp=a+b*T+c*T**2  [kJ/K/mol]
! local variables --------------------------------------------------------------
!-------------------------------------------------------------------------------
  subroutine pucl3_abc(a,b,c)
    real(8),intent(out) :: a,b,c
    
    a = 144.0d-3
    b = 0.0d-3
    c = 0.0d-3
  endsubroutine pucl3_abc
  
!-------------------------------------------------------------------------------
! calculate molar enthalpy of PuCl3
!
! arguments --------------------------------------------------------------------
! T        temperature              [K]
! h_pucl3  molar enthalpy of PuCl3  [kJ/mol]
! local variables --------------------------------------------------------------
! a     coefficient for calculation of cp=a+b*T+c*T**2  [kJ/K/mol]
! b     coefficient for calculation of cp=a+b*T+c*T**2  [kJ/K/mol]
! c     coefficient for calculation of cp=a+b*T+c*T**2  [kJ/K/mol]
!-------------------------------------------------------------------------------
  real(8) function h_pucl3(T)
    real(8) :: T
    real(8) :: a,b,c
    
    call pucl3_abc(a,b,c)
    h_pucl3 = h_abc(T,a,b,c)
  endfunction h_pucl3

!-------------------------------------------------------------------------------
! calculate molar enthalpy given coefficients a, b, and c.
! absolute enthalpies referenced to T0=273.15K
! equation of the form h = integral(cp(T),from=273.15K,to=T)
!
! arguments --------------------------------------------------------------------
! T      temperature in Kelvin
! a      coefficient for calculation of cp=a+b*T+c*T**2  [kJ/K/mol]
! b      coefficient for calculation of cp=a+b*T+c*T**2  [kJ/K/mol]
! c      coefficient for calculation of cp=a+b*T+c*T**2  [kJ/K/mol]
! h_abc  molar enthalpy at temperature T                 [kJ/mol]
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
! should be included in x_ucl3 i.e. x_ucl3 = x_UCl3 + x_PuCl3
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
! calculate speicific enthalpy of PuCl3-UCl3-NaCl system
!
! arguments --------------------------------------------------------------------
! T       temperature of system     [K]
! calc_h  specific enthalpy of system  [kJ/kg]
!   u_nacl   user specified mole fraction of NaCl
!   u_ucl3   user specified mole fraction of UCl3
!   u_pucl3  user specified mole fraction of PuCl3
! local variables --------------------------------------------------------------
! x_nacl   mole fraction of NaCl
! x_ucl3   mole fraction of UCl3
! x_pucl3  mole fraction of PuCl3
!-------------------------------------------------------------------------------
  real(8) function calc_h(T,u_nacl,u_ucl3,u_pucl3)
    character(*),parameter :: myName = 'calc_h'
    real(8) :: T
    real(8),optional :: u_nacl,u_ucl3,u_pucl3
    real(8) :: x_nacl,x_ucl3,x_pucl3
    
    ! set default mole fractions
    call x_default(x_nacl,x_ucl3,x_pucl3)
    if (present(u_nacl)) then
      if ((.not. present(u_pucl3)) .or. (.not. present(u_nacl))) then
        msg = ''
        write(msg(1),f1) 'if one user specified mole fraction is' // &
          ' present then, all must be present'
        call raise_fatal(modName,myName,msg)
      elseif ((u_nacl + u_ucl3 + u_pucl3) /= 1.0d0) then
        msg = ''
        write(msg(1),f1) 'sum of user specified mole fraction /= 1.0d0'
        write(msg(2),f2) 'sum = ',(u_ucl3 + u_pucl3 + u_nacl)
        call raise_fatal(modName,myName,msg)
      endif
      ! overwrite default mole fractions
      x_nacl = u_nacl
      x_ucl3 = u_ucl3
      x_pucl3 = u_pucl3
    endif
    ! calculate molar enthalpy of system
    calc_h = x_nacl * h_nacl(T) + x_ucl3 * h_ucl3(T) + x_pucl3 * h_pucl3(T) + &
      hmix_ucl3_nacl(x_ucl3 + x_pucl3)
    ! convert molar enthalpy to specific enthalpy
    calc_h = calc_h * (1.0d0 / molar_mass(x_nacl,x_ucl3,x_pucl3))
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
    logical :: found
    
    if ((x > arr(size(arr))) .or. (x < arr(1))) then
      msg = ''
      write(msg(1),f1) 'x out of bounds for interpolation'
      write(msg(2),f2) 'x = ', x
      write(msg(3),'(2(a,e12.6))') 'arr(1) = ',arr(1), &
        ' arr(size(arr)) = ',arr(size(arr))
      call raise_fatal(modName,myName,msg)
    endif    
    lower = 1
    upper = size(arr)
    found = .false.
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
        write(msg(1),f1) 'trouble searching'
        call raise_fatal(modName,myName,msg)
      endif
      if ((x > arr(i)) .and. (x < arr(i + 1))) then
        found = .true.
        ilo = i
        weight = (x - arr(i)) / (arr(i + 1) - arr(i))
      endif
    enddo
    if ((weight > 1.0d0) .or. (weight < 0.0d0)) then
      msg = ''
      write(msg(1),f1) 'bad weight'
      write(msg(2),f2) 'weight = ',weight
      call raise_fatal(modName,myName,msg)
    elseif ((ilo < 1) .or. (ilo > (size(arr) - 1))) then
      msg = ''
      write(msg(1),f1) 'bad ilo'
      write(msg(2),'(a,i10)') 'ilo = ',ilo
    endif
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
  
!-------------------------------------------------------------------------------
! solve a general cubic function
! alpha*x**3+beta*x**2+gamma*x+delta=0
! currently only supports functions with one multiple-root
! TO-DO: consider moving this to another module
! 
! arguments --------------------------------------------------------------------
! alpha  function coefficient
! beta      "         "
! gamma     "         "
! delta     "         "
! local variables --------------------------------------------------------------
! disc   discriminant for true cubic and true quadratic
! disc0  discriminant0 for true cubic
!-------------------------------------------------------------------------------
  real(8) function cubic_solve(alpha,beta,gamma,delta,umin,umax)
    character(*),parameter :: myName = 'cubic_solve'
    real(8) :: alpha,beta,gamma,delta
    real(8),optional :: umin,umax
    real(8) :: disc,disc0
    real(8) :: r1,r2,r3
    
    if (alpha /= 0.0d0) then
      ! true cubic
      disc = 18.0d0 * alpha * beta * gamma * delta - &
        4.0d0 * beta ** 3.0d0 * delta + beta ** 2.0d0 * gamma ** 2.0d0 - &
        4.0d0 * alpha * gamma ** 3.0d0 - &
        27.0d0 * alpha ** 2.0d0 * delta ** 2.0d0
      disc0 = beta ** 2.0d0 - 3.0d0 * alpha * gamma
      if ((disc == 0.0d0) .and. (disc0 == 0.0d0)) then
        ! equation has a single, triple-root
        cubic_solve = (-1.0d0) * (beta / (3.0d0 * alpha))
      ! elseif (isRealPos) then
        ! TO-DO: add multiple root sols for cubic
      else
        msg = ''
        write(msg(1),f1) 'equation is a true cubic but has multiple sols'
        write(msg(2),f1) 'it is required that disc=disc0=0.0 for unique sol'
        write(msg(3),f2) 'disc  = ',disc
        write(msg(4),f2) 'disc0 = ',disc0
        write(msg(5),f1) 'i can try to force a real & positive solution by' // &
          ' passing RealPos=.true.'
        call raise_fatal(modName,myName,msg)
      endif
    elseif (beta /= 0.0d0) then
      ! true quadratic
      disc = gamma ** 2.0d0 - 4.0d0 * beta * delta
      if (disc == 0.0d0) then
        ! equation has a single, double root
        cubic_solve = (-1.0d0) * (gamma / (2.0d0 * beta))
      elseif (disc < 0.0d0) then
        write(msg(1),f1) 'equation is a true quadratic'
        write(msg(2),f1) 'this quadratic has two imaginary roots' // &
          'because disc < 0.0'
        write(msg(3),f2) 'disc = ',disc
        call raise_fatal(modName,myName,msg)
      else
        r1 = (((-1.0d0) * gamma) + sqrt(disc)) / (2.0d0 * beta)
        r2 = (((-1.0d0) * gamma) - sqrt(disc)) / (2.0d0 * beta)
        if (present(umin)) then
          if (min(r1,r2) < umin) then
            cubic_solve = max(r1,r2)
          elseif (present(umax)) then
            if (max(r1,r2) > umax) then
              cubic_solve = min(r1,r2)
            else
              msg = ''
              write(msg(1),f1) 'equation is true quadratic and user specified' // &
                ' umin and umax'
              write(msg(2),f1) 'however these options didnt eliminate any sols'
              write(msg(3),f2) 'umin = ',umin
              write(msg(4),f2) 'umax = ',umax
              write(msg(5),f2) 'r1 = ',r1
              write(msg(6),f2) 'r2 = ',r2
              call raise_fatal(modName,myName,msg)
            endif
          else
            msg = ''
            write(msg(1),f1) 'equation is true quadratic and user specified umin'
            write(msg(2),f1) 'however this option didnt eliminate any sols'
            write(msg(3),f2) 'umin = ',umin
            write(msg(4),f2) 'r1 = ',r1
            write(msg(5),f2) 'r2 = ',r2
            call raise_fatal(modName,myName,msg)
          endif
        elseif (present(umax)) then
          if (max(r1,r2) > umax) then
            cubic_solve = min(r1,r2)
          else
            msg = ''
            write(msg(1),f1) 'equation is true quadratic and user specified umax'
            write(msg(2),f1) 'however this option didnt eliminate any sols'
            write(msg(3),f2) 'umax = ',umax
            write(msg(4),f2) 'r1 = ',r1
            write(msg(5),f2) 'r2 = ',r2
            call raise_fatal(modName,myName,msg)
          endif
        else
          msg = ''
          write(msg(1),f1) 'equation is true quadratic but has multiple sols'
          write(msg(2),f1) 'it is required that disc=0.0 for unique sol'
          write(msg(3),f2) 'disc = ',disc
          write(msg(4),f1) 'i can try to force a solution by' // &
            ' speicifying umin or umax'
          call raise_fatal(modName,myName,msg)
        endif
      endif
    elseif (gamma /= 0.0d0) then
      ! true linear
      cubic_solve = (-1.0d0) * (delta / gamma)
    else
      msg = ''
      write(msg(1),f1) 'arguments imply a constant function with no root!'
      write(msg(2),f1) 'alpha*x**3+beta*x**2+gamma*x+delta=0'
      write(msg(3),f2) 'alpha = ',alpha
      write(msg(4),f2) 'beta  = ',beta
      write(msg(5),f2) 'gamma = ',gamma
      write(msg(6),f2) 'delta = ',delta
      call raise_fatal(modName,myName,msg)
    endif
  endfunction cubic_solve

!-------------------------------------------------------------------------------
! calculate temperature of PuCl3-UCl3-NaCl system
!
! arguments --------------------------------------------------------------------
! h_in    speicifc enthalpy input  [kJ/kg]
! calc_T  temperature              [K]
!   u_nacl   user specified mole fraction of NaCl
!   u_ucl3   user specified mole fraction of UCl3
!   u_pucl3  user specified mole fraction of PuCl3
! local variables --------------------------------------------------------------
! hm  molar enthalpy [kJ/mol]
! h   specific enthalpy [kJ/kg]
! suma  sum of a coefficinents from cp calculation cp=a+b*T+c*T**2
! sumb  sum of b coefficinents from cp calculation cp=a+b*T+c*T**2
! sumc  sum of c coefficinents from cp calculation cp=a+b*T+c*T**2
! x_nacl   mole fraction of NaCl
! x_ucl3   mole fraction of UCl3
! x_pucl3  mole fraction of PuCl3
!-------------------------------------------------------------------------------
  real(8) function calc_T(hin,u_nacl,u_ucl3,u_pucl3)
    character(*),parameter :: myName = 'calc_T'
    real(8) :: hin
    real(8),optional :: u_nacl,u_ucl3,u_pucl3
    real(8) :: hm,h
    real(8) :: suma,sumb,sumc
    real(8) :: nacl_a,nacl_b,nacl_c    ! used for sums
    real(8) :: ucl3_a,ucl3_b,ucl3_c    ! used for sums
    real(8) :: pucl3_a,pucl3_b,pucl3_c ! used for sums
    real(8) :: x_nacl,x_ucl3,x_pucl3
    
    h = hin
    call x_default(x_nacl,x_ucl3,x_pucl3)
    if (present(u_nacl)) then
      if ((.not. present(u_pucl3)) .or. (.not. present(u_nacl))) then
        msg = ''
        write(msg(1),f1) 'if one user specified mole fraction is' // &
          ' present then, all must be present'
        call raise_fatal(modName,myName,msg)
      elseif ((u_nacl + u_ucl3 + u_pucl3) /= 1.0d0) then
        msg = ''
        write(msg(1),f1) 'sum of user specified mole fraction /= 1.0d0'
        write(msg(2),f2) 'sum = ',(u_ucl3 + u_pucl3 + u_nacl)
        call raise_fatal(modName,myName,msg)
      endif
      ! overwrite default mole fractions
      x_nacl = u_nacl
      x_ucl3 = u_ucl3
      x_pucl3 = u_pucl3
    endif
    ! TO-DO: fix temperature in nacl_abc call
    call nacl_abc(1000.0d0,nacl_a,nacl_b,nacl_c)
    call ucl3_abc(ucl3_a,ucl3_b,ucl3_c)
    call pucl3_abc(pucl3_a,pucl3_b,pucl3_c)
    suma = x_nacl * nacl_a + x_ucl3 * ucl3_a + x_pucl3 * pucl3_a
    sumb = x_nacl * nacl_b + x_ucl3 * ucl3_b + x_pucl3 * pucl3_b
    sumc = x_nacl * nacl_c + x_ucl3 * ucl3_c + x_pucl3 * pucl3_c
    hm = h * molar_mass(x_nacl,x_ucl3,x_pucl3)
    ! 0 = sumc*(1/3)*T**3+sumb*(1/2)*T**2+suma*T+(hmix-h)
    sumc = sumc * (1.0d0 / 3.0d0)
    sumb = sumb * 0.5d0
    calc_T = cubic_solve(sumc,sumb,suma,(hmix_ucl3_nacl(x_ucl3 + x_pucl3) - hm), &
      umin=298.0d0,umax=2.0d3)
    msg = ''
    write(msg(1),f1) 'calc_T not yet supported'
    call raise_warning(modName,myName,msg)
  endfunction calc_T

!-------------------------------------------------------------------------------
! return default mole fractions for soft reactor
! based on 1*PuCl3-8*UCl3-10*NaCl from 
! "Reator with very low fission product inventory" M. Taube, W. Heer, 1980
!
! arguments --------------------------------------------------------------------
! x_nacl   mole fraction NaCl
! x_ucl3   mole fraction UCl3
! x_pucl3  mole fraction PuCl3}
! local variables --------------------------------------------------------------
!-------------------------------------------------------------------------------
  subroutine x_default(x_nacl,x_ucl3,x_pucl3)
    real(8),intent(out) :: x_nacl,x_ucl3,x_pucl3
    
    x_nacl  = (10.0d0 / 19.0d0)
    x_ucl3  = (8.0d0  / 19.0d0)
    x_pucl3 = (1.0d0  / 19.0d0)
  endsubroutine x_default

!-------------------------------------------------------------------------------
! calculate molar mass of a PuCl3-UCl3-NaCl system given mole fractions
! output value in [kg/mol]
!
! arguments --------------------------------------------------------------------
! x_nacl   mole fraction NaCl
! x_ucl3   mole fraction UCl3
! x_pucl3  mole fraction PuCl3
! local variables --------------------------------------------------------------
! m_nacl   molar mass NaCl   [gm/mol]
! m_ucl3   molar mass UCl3   [gm/mol]
! m_pucl3  molar mass PuCl3  [gm/mol]
!-------------------------------------------------------------------------------
  real(8) function molar_mass(x_nacl,x_ucl3,x_pucl3)
    character(*),parameter :: myName = 'molar_mass'
    real(8) :: x_nacl,x_ucl3,x_pucl3
    real(8) :: m_nacl,m_ucl3,m_pucl3
    
    ! check to make sure mole fractions sum to 1.0
    if ((x_nacl + x_ucl3 + x_pucl3) /= 1.0d0) then
      msg = ''
      write(msg(1),f1) 'sum of mole fractions /= 1.0d0'
      write(msg(2),f2) 'sum = ',x_nacl + x_ucl3 + x_pucl3
      call raise_fatal(modName,myName,msg)
    endif
    ! calculate molar masses based on elemental consituents
    m_nacl = 22.989769280d0 + 35.4530d0
    m_ucl3 = 238.028910d0 + 3.0d0 * 35.4530d0
    m_pucl3 = 238.0495599d0 + 3.0d0 * 35.4530d0
    ! calculate molar mass and convert from gm to kg
    molar_mass = (x_nacl * m_nacl + x_ucl3 * m_ucl3 + x_pucl3 * m_pucl3) * (1.0d-3)
  endfunction molar_mass
endmodule saltprops