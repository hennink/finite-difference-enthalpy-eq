! Solver for the non-linear ODE
!       d/dt (rho * h)  =  - lambda * h + Q ,
! 
!       h  =  cp * T - h0
!          =  cp * (T - T0)
! where 
!       lambda,         is a constant,
!       h0, T0, cp      are constants,
!       t               is the independent variable (the time),
!       T               is the temperature (the unknown),
!       h = h(t)        is the specific enthalpy,
!       rho = rho(T)    is the density, which is an arbitrary function of T,
!       Q = Q(t)        is the source term.
! Uses a finite difference method with various techniques to deal with the non-linearity in the ODE. 
! 
! This code accompanies the paper
! 
!     "A Pressure-based Solver for Low-Mach Number Flow using a Discontinuous Galerkin Method", A. Hennink, M. Tiberga, D. Lathouwers
! 
! that was submitted to the Journal of Computational Fluids on 2019-11.
! 
module kinds_mod
    use, intrinsic :: iso_fortran_env
    implicit none

    integer, parameter :: WP = REAL128
    
    real(wp), parameter :: UNINITIALIZED = -huge(1.0_wp)
end module


module props_mod
    use kinds_mod
    implicit none

    real(wp) :: cp = UNINITIALIZED
    real(wp) :: T0 = UNINITIALIZED
    real(wp) :: rho0 = UNINITIALIZED, rho1 = UNINITIALIZED

    procedure(T2rho__affine), pointer :: T2rho => null()
    procedure(T2deriv_rho_T__affine), pointer :: T2deriv_rho_T => null()

contains
    ! h <--> T relationships:
    pure real(wp) function T2h(T) result(h)
        real(wp), intent(in) :: T

        h = cp * (T - T0)
    end function

    pure real(wp) function h2T(h) result(T)
        real(wp), intent(in) :: h

        T = T0 + h / cp
    end function
    
    pure real(wp) function T2deriv_h_T(T) result(deriv_h_T)
        real(wp), intent(in) :: T

        deriv_h_T = cp
    end function
    

    ! 'mixing' fluid prop for rho <--> T:
    pure real(wp) function T2rho__mixing(T) result(rho)
        real(wp), intent(in) :: T

        rho = 1 / (T / rho0  +  (1-T) / rho1)
    end function

    pure real(wp) function T2deriv_rho_T__mixing(T) result(deriv_rho_T)
        real(wp), intent(in) :: T

        deriv_rho_T = - rho0 * rho1 * (rho1 - rho0) / (rho0*(1-T) + rho1*T)**2
    end function


    ! 'affine' fluid prop for rho <--> T:
    pure real(wp) function T2rho__affine(T) result(rho)
        real(wp), intent(in) :: T

        rho = rho0 * T + rho1 * (1-T)
    end function

    pure real(wp) function T2deriv_rho_T__affine(T) result(deriv_rho_T)
        real(wp), intent(in) :: T

        deriv_rho_T = rho0 - rho1
    end function


    ! 'constant' fluid prop for rho:
    pure real(wp) function T2rho__constant(T) result(rho)
        real(wp), intent(in) :: T

        rho = rho0
    end function

    pure real(wp) function T2deriv_rho_T__constant(T) result(deriv_rho_T)
        real(wp), intent(in) :: T

        deriv_rho_T = 0.0_wp
    end function
end module



module mansol_mod
    ! Manufactured solution.
    ! Contains exact functions for T, h, and the corresponding external source (Q).
    use kinds_mod
    use props_mod
    implicit none

    private
    public :: lambda
    public :: ex_T, ex_deriv_T, ex_h, ex_src

    real(wp) :: lambda = UNINITIALIZED
    real(wp), parameter :: PI = 3.14159265358979323846264338327950288419716939937510582097494459_wp
    
    real(wp), parameter :: AMPLITUDE = 0.1_wp
    real(wp), parameter ::  TMIN = 0.5_wp - AMPLITUDE, &
                            TMAX = 0.5_wp + AMPLITUDE   ! adjust manually

contains
    pure real(wp) function ex_T(time)
        ! Don't get too close to 0 or 1, or the extrapolation will give negative densities:
        real(wp), intent(in) :: time

        ex_T = 0.5_wp + AMPLITUDE * sin(2*PI*time)
    end function

    pure real(wp) function ex_deriv_T(time)
        real(wp), intent(in) :: time

        ex_deriv_T = AMPLITUDE * 2 * PI * cos(2*PI*time)
    end function

    pure real(wp) function ex_h(time)
        real(wp), intent(in) :: time

        ex_h = T2h(ex_T(time))
    end function

    pure real(wp) function ex_src(time)
        real(wp), intent(in) :: time

        real(wp) :: T, deriv_T, rho, h, deriv_volh_T

        T = ex_T(time)
        deriv_T = ex_deriv_T(time)

        rho = T2rho(T)
        h = T2h(T)
        deriv_volh_T = rho * T2deriv_h_T(T) + h * T2deriv_rho_T(T)
        
        ex_src = lambda * h + deriv_volh_T * deriv_T
    end function
end module



module time_stepping_mod
    use kinds_mod
    implicit none

    type error_t
        real(wp) :: relative_error
        logical :: pos_deriv_Hh, nonpos_deriv_Hh, &
                   nonpos_implicit_weight, &
                   h_too_low, h_too_high
    end type

contains
    function rel_T_error_after_time_steps(nsteps,derivhr_strategy,order_extrapolation,order_BDF) result(error)
        ! Uses a finite difference scheme to solve the ODE on the domain
        !       0 <= t <= 1 .
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf, ieee_negative_inf
        use props_mod, only: h2T, T2rho, T2deriv_rho_T, T2deriv_h_T
        use mansol_mod, only: lambda, ex_src, ex_h, ex_T
        
        integer, intent(in) :: nsteps
        character(*), intent(in) :: derivhr_strategy
        integer, intent(in) :: order_extrapolation, order_BDF
        type(error_t) :: error
        
        real(wp) :: dt
        real(wp), allocatable :: BDF_weights(:) ! without (1/dt)
        real(wp), allocatable :: extrapolation_weights(:)
        
        real(wp), allocatable :: num_h(:)
        
        integer :: i, k
        real(wp) :: explicit_time_terms, old_h, old_T, old_rho
        real(wp) :: predict_h, predict_rho, predict_T, &
                    predict_deriv_volh_h, predict_beta, predict_deriv_rho_h, predict_deriv_rho_T, predict_cp
        real(wp) :: src_new_time
        real(wp) :: impl_part, expl_part
        real(wp) :: final_ex_T

        ! Initializations:
        ! (Defining plenty of previous values eliminates the startup error.)
        dt = 1.0_wp / nsteps
        allocate(num_h(-10:nsteps))
        do k = lbound(num_h,1), 0
            num_h(k) = ex_h(time=k*dt)
        enddo
        select case(order_BDF)
            case(1);        BDF_weights = [1, -1]
            case(2);        BDF_weights = [3, -4, 1] / 2.0_wp
            case(3);        BDF_weights = [11, -18, 9, -2] / 6.0_wp
            case default;   error stop 'unsupported BDF order'
        end select
        select case(order_extrapolation)
            case(1);        extrapolation_weights = [1]
            case(2);        extrapolation_weights = [2, -1]
            case(3);        extrapolation_weights = [3, -3, 1]
            case(4);        extrapolation_weights = [4, -6, 4, -1]
            case default;   error stop 'unsupported order of extrapolation'
        end select
        error%h_too_low = .false.
        error%h_too_high = .false.
        error%nonpos_implicit_weight = .false.
        error%nonpos_deriv_Hh = .false.
        error%pos_deriv_Hh = .false.

        ! Time stepping:
        do i = 1, nsteps
            ! Terms in the finite difference scheme that depend on the previous time steps:
            explicit_time_terms = 0.0_wp
            do k = 1, size(BDF_weights) - 1
                old_h = num_h(i-k)
                old_T = h2T(old_h)
                old_rho = T2rho(old_T)
                explicit_time_terms = explicit_time_terms + BDF_weights(k+1) * old_rho * old_h
            enddo
            
            ! Find the reference point (the predictor):
            predict_h = 0.0_wp
            do k = 1, size(extrapolation_weights)
                predict_h = predict_h + extrapolation_weights(k) * num_h(i-k)
            enddo
            predict_T = h2T(predict_h)
            predict_rho = T2rho(predict_T)
            predict_cp = T2deriv_h_T(predict_T)
            predict_deriv_rho_T = T2deriv_rho_T(predict_T)
            predict_beta = - predict_deriv_rho_T / predict_rho
            predict_deriv_volh_h = predict_rho + predict_h * predict_deriv_rho_T / predict_cp
            predict_deriv_rho_h = predict_deriv_rho_T / predict_cp
            
            ! The actual 'solve':
            src_new_time = ex_src(time=i*dt)
            select case(derivhr_strategy)
                case('deriv_rh')
                    impl_part = BDF_weights(1) * predict_deriv_volh_h + dt * lambda
                    expl_part = BDF_weights(1) * predict_h**2 * predict_deriv_rho_h  -  explicit_time_terms  +  dt * src_new_time
                case('rho_predictor')
                    impl_part = BDF_weights(1) * predict_rho + dt * lambda
                    expl_part = - explicit_time_terms + dt * src_new_time
                case('mansol_rho')
                    impl_part = BDF_weights(1) * t2rho(ex_T(time=i*dt)) + dt * lambda
                    expl_part = - explicit_time_terms + dt * src_new_time
                case default
                    error stop 'unknown derivhr_strategy'
            end select
            
            ! Check whether certain conditions were met:
            if (predict_h <= -predict_cp / predict_beta) error%h_too_low = .true.
            if (predict_h >= +predict_cp / predict_beta) error%h_too_high = .true.
            if (predict_deriv_volh_h <= 0) error%nonpos_deriv_Hh = .true.
            if (predict_deriv_volh_h > 0) error%pos_deriv_Hh = .true.
            if (impl_part <= 0) error%nonpos_implicit_weight = .true.
            num_h(i) = expl_part / impl_part
        enddo
        
        final_ex_T = ex_T(time=nsteps*dt)
        error%relative_error = (h2T(num_h(nsteps)) - final_ex_T) / final_ex_T
         ! Avoid large exponents, which 64-bit plotting software (e.g., Python) cannot read:
        if (error%relative_error < -100) then
            error%relative_error = ieee_value(error%relative_error,ieee_negative_inf)
        elseif (error%relative_error > 100) then
            error%relative_error = ieee_value(error%relative_error,ieee_positive_inf)
        endif
    end function rel_T_error_after_time_steps
end module


program fd
    ! Prints the error in the time-stepping scheme a series of T0 values.
    ! 
    ! Usage:
    !     ./fd fluid_props=affine lambda=0.1 cp=1.0 rho0=0.5 rho1=2.0    T0=1.0    derivhr_strategy=deriv_rh order_BDF=2 order_EX=3 nsteps=2048
    use kinds_mod
    use props_mod
    use mansol_mod, only: lambda
    use time_stepping_mod, only: error_t, rel_T_error_after_time_steps
    implicit none
    
    character(:), allocatable :: props, derivhr_strategy
    integer :: order_extrapolation, order_BDF
    integer :: nsteps
    type(error_t) :: error
    
    ! Fluid properties:
    props = cl_value("fluid_props")
    lambda = cl_value_as_REAL("lambda")
    cp = cl_value_as_REAL("cp")
    rho0 = cl_value_as_REAL("rho0")
    rho1 = cl_value_as_REAL("rho1")

    T0 = cl_value_as_REAL("T0")

    ! Solution strategy:
    derivhr_strategy = cl_value("derivhr_strategy")
    order_BDF = cl_value_as_int("order_BDF")
    order_extrapolation = cl_value_as_int("order_EX")

    nsteps = cl_value_as_int("nsteps")

    select case(props)
    case('constant')
        T2rho => T2rho__constant
        T2deriv_rho_T => T2deriv_rho_T__constant
    case('affine')
        T2rho => T2rho__affine
        T2deriv_rho_T => T2deriv_rho_T__affine
    case('mixing')
        T2rho => T2rho__mixing
        T2deriv_rho_T => T2deriv_rho_T__mixing
    case default
        error stop 'unknown fluid property'
    end select

    error = rel_T_error_after_time_steps(           &
        nsteps = nsteps,                            &
        derivhr_strategy = derivhr_strategy,        &
        order_extrapolation = order_extrapolation,  &
        order_BDF = order_BDF                       &
    )
    write(*,"(es20.2e4, 7i5)")                              &
        error%relative_error,                               &
        merge(1, 0, [                                       &
            error%pos_deriv_Hh .and. error%nonpos_deriv_Hh, &
            error%nonpos_deriv_Hh,                          &
            error%nonpos_implicit_weight,                   &
            error%h_too_low, error%h_too_high               &
        ])

contains
    real(wp) function cl_value_as_real(option) result(val)
        ! Return the value that was provided for 'option' on the command line.
        ! 
        ! This fails to error with inputs like "23,2" or "2.3;4", but I don't know how to read an 
        ! unknown real format with a stricter format descriptor than '*'.
        character(*), intent(in) :: option

        character(:), allocatable :: as_char

        as_char = cl_value(option)
        read(as_char,*) val
    end function

    integer function cl_value_as_int(option) result(val)
        ! Return the value that was provided for 'option' on the command line.
        character(*), intent(in) :: option

        character(:), allocatable :: as_char
        character(10) :: fmt

        as_char = cl_value(option)
        write(fmt,"( '(i', i0, ')' )") len(as_char)
        read(as_char,fmt) val
    end function

    function cl_value(option)
        ! Return the value that was provided for 'option' on the command line.
        character(*), intent(in) :: option
        character(:), allocatable :: cl_value

        character(:), allocatable, save :: all_args(:,:)
        integer :: i

        if (.not. allocated(all_args)) all_args = cl_args()

        i = findloc(all_args(1,:),option,dim=1)
        call assert(i > 0, "could not find this option in the command line arguments; option=" // option)
        cl_value = trim(all_args(2,i))
    end function

    function cl_args() result(names_vals)
        ! Return an array of (name,value) character pairs of command line arguments.
        ! The command line arguments are expected to be of the form
        !     option=foo  other_option=bar  n=3  flag=T  ...
        character(:), allocatable :: names_vals(:,:)

        character(:), allocatable :: raw_args(:), arg
        integer :: i, n, splitloc
        integer :: status
        integer :: maxlen
        integer, allocatable :: lengths(:)
        
        n = command_argument_count()
        
        ! Get a list of arguments in raw (i.e., unparsed) form:
        allocate(lengths(n))
        do i = 1, n
            call get_command_argument(i,length=lengths(i),status=status)
            call assert(status == 0, "bad status in get_command_argument")
        enddo
        maxlen = maxval(lengths)
        allocate(character(maxlen) :: raw_args(n))
        do i = 1, n
            call get_command_argument(i,value=raw_args(i))
            call assert(status == 0, "bad status in get_command_argument")
        enddo

        ! Split the raw arguments into names and values:
        allocate(character(maxlen) :: names_vals(2,n)) ! over-dimensionalized, but does not matter.
        do i = 1, n
            arg = raw_args(i)
            splitloc = index(arg,'=')
            call assert(splitloc > 0, "could not find '=' in command line argument; arg=" // arg)
            names_vals(1,i) = arg(:splitloc-1)
            names_vals(2,i) = arg(splitloc+1:)
            deallocate(arg)
        enddo
    end function

    subroutine assert(statement,error_msg)
        use, intrinsic :: iso_fortran_env, only: ERROR_UNIT
        logical, intent(in) :: statement
        character(*), intent(in) :: error_msg

        if (.not. statement) then
            write(ERROR_UNIT,"(a)") error_msg
            error stop
        endif
    end subroutine
end program
