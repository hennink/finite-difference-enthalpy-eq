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
! This code accomanies the paper
! 
!     "A Pressure-based Solver for Low-Mach Number Flow using a Discontinuous Galerkin Method", A. Hennink, M. Tiberga, D. Lathouwers
! 
! that was submitted to the Journal of Computational Fluids on 2019-11.
! 
module kinds_mod
    use, intrinsic :: iso_fortran_env
    implicit none

    integer, parameter :: WP = REAL128
    
    real(WP), parameter :: UNINITIALIZED = -huge(1.0_WP)
end module


module props_mod
    use kinds_mod
    implicit none

    real(WP) :: cp = UNINITIALIZED
    real(WP) :: T0 = UNINITIALIZED
    real(WP) :: rho0 = UNINITIALIZED, rho1 = UNINITIALIZED

    procedure(T2rho__affine), pointer :: T2rho => null()
    procedure(T2deriv_rho_T__affine), pointer :: T2deriv_rho_T => null()

contains
    ! h <--> T relationships:
    pure real(WP) function T2h(T) result(h)
        real(WP), intent(in) :: T

        h = cp * (T - T0)
    end function

    pure real(WP) function h2T(h) result(T)
        real(WP), intent(in) :: h

        T = T0 + h / cp
    end function
    
    pure real(WP) function T2deriv_h_T(T) result(deriv_h_T)
        real(WP), intent(in) :: T

        deriv_h_T = cp
    end function
    

    ! 'mixing' fluid prop for rho <--> T:
    pure real(WP) function T2rho__mixing(T) result(rho)
        real(WP), intent(in) :: T

        rho = 1 / (T / rho0  +  (1-T) / rho1)
    end function

    pure real(WP) function T2deriv_rho_T__mixing(T) result(deriv_rho_T)
        real(WP), intent(in) :: T

        deriv_rho_T = - rho0 * rho1 * (rho1 - rho0) / (rho0*(1-T) + rho1*T)**2
    end function


    ! 'affine' fluid prop for rho <--> T:
    pure real(WP) function T2rho__affine(T) result(rho)
        real(WP), intent(in) :: T

        rho = rho0 * T + rho1 * (1-T)
    end function

    pure real(WP) function T2deriv_rho_T__affine(T) result(deriv_rho_T)
        real(WP), intent(in) :: T

        deriv_rho_T = rho0 - rho1
    end function


    ! 'constant' fluid prop for rho:
    pure real(WP) function T2rho__constant(T) result(rho)
        real(WP), intent(in) :: T

        rho = rho0
    end function

    pure real(WP) function T2deriv_rho_T__constant(T) result(deriv_rho_T)
        real(WP), intent(in) :: T

        deriv_rho_T = 0.0_WP
    end function
end module



module mansol_mod
    ! Manufactured solution.
    ! Contains exact functions for T, h, and the corresponding external source (Q).
    use kinds_mod
    use props_mod
    implicit none

    real(WP) :: lambda = UNINITIALIZED
    real(WP), parameter :: PI = 3.1415926535897932384626433_WP
    
    real(WP), parameter :: AMPLITUDE = 0.1_WP
    real(WP), parameter ::  TMIN = 0.5_WP - AMPLITUDE, &
                            TMAX = 0.5_WP + AMPLITUDE   ! adjust manually

contains
    pure real(WP) function ex_T(time)
        ! Don't get too close to 0 or 1, or the extrapolation will give negative densities:
        real(WP), intent(in) :: time

        ex_T = 0.5_WP + AMPLITUDE * sin(2*PI*time)
    end function

    pure real(WP) function ex_deriv_T(time)
        real(WP), intent(in) :: time

        ex_deriv_T = AMPLITUDE * 2 * PI * cos(2*PI*time)
    end function

    pure real(WP) function ex_h(time)
        real(WP), intent(in) :: time

        ex_h = T2h(ex_T(time))
    end function

    pure real(WP) function ex_src(time)
        real(WP), intent(in) :: time

        real(WP) :: T, deriv_T, rho, h, deriv_volh_T

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

contains
    function rel_T_error_after_time_steps(nsteps,derivhr_strategy,order_extrapolation,order_BDF) result(error_msg)
        ! Uses a finite difference scheme to solve the ODE on the domain
        !       0 <= t <= 1 .
        ! Returns a character, which contains
        !       * the error in the temperature (T) at t=1;
        !       * a series of 0/1, which indicate whether certain conditions were met during the calculation.
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf
        use props_mod, only: h2T, T2rho, T2deriv_rho_T, T2deriv_h_T
        use mansol_mod, only: lambda, ex_src, ex_h, ex_T
        
        integer, intent(in) :: nsteps
        character(*), intent(in) :: derivhr_strategy
        integer, intent(in) :: order_extrapolation, order_BDF
        character(1000) :: error_msg
        
        real(WP) :: dt
        real(WP), allocatable :: BDF_weights(:) ! without (1/dt)
        real(WP), allocatable :: extrapolation_weights(:)
        
        real(WP), allocatable :: num_h(:)
        
        integer :: i, k
        real(WP) :: explicit_time_terms, old_h, old_T, old_rho
        real(WP) :: predict_h, predict_rho, predict_T, &
                    predict_deriv_volh_h, predict_beta, predict_deriv_rho_h, predict_deriv_rho_T, predict_cp
        real(WP) :: src_new_time
        real(WP) :: impl_part, expl_part
        real(WP) :: final_ex_T, error
        logical :: h_too_low, h_too_high, nonpos_impl_weight, nonpos_deriv_Hh, pos_deriv_Hh

        ! Initializations:
        ! (Defining plenty of previous values eliminates the startup error.)
        dt = 1.0_wp / nsteps
        allocate(num_h(-10:nsteps))
        do k = lbound(num_h,1), 0
            num_h(k) = ex_h(time=k*dt)
        enddo
        select case(order_BDF)
            case(1);        BDF_weights = [1, -1]
            case(2);        BDF_weights = [3, -4, 1] / 2.0_WP
            case(3);        BDF_weights = [11, -18, 9, -2] / 6.0_WP
            case default;   error stop 'unsupported BDF order'
        end select
        select case(order_extrapolation)
            case(1);        extrapolation_weights = [1]
            case(2);        extrapolation_weights = [2, -1]
            case(3);        extrapolation_weights = [3, -3, 1]
            case(4);        extrapolation_weights = [4, -6, 4, -1]
            case default;   error stop 'unsupported order of extrapolation'
        end select
        h_too_low = .false.
        h_too_high = .false.
        nonpos_impl_weight = .false.
        nonpos_deriv_Hh = .false.
        pos_deriv_Hh = .false.

        ! Time stepping:
        do i = 1, nsteps
            ! Terms in the finite difference scheme that depend on the previous time steps:
            explicit_time_terms = 0.0_WP
            do k = 1, size(BDF_weights) - 1
                old_h = num_h(i-k)
                old_T = h2T(old_h)
                old_rho = T2rho(old_T)
                explicit_time_terms = explicit_time_terms + BDF_weights(k+1) * old_rho * old_h
            enddo
            
            ! Find the reference point (the predictor):
            predict_h = 0.0_WP
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
            if (-predict_h >= predict_cp / predict_beta) h_too_low = .true.
            if (+predict_h >= predict_cp / predict_beta) h_too_high = .true.
            if (predict_deriv_volh_h <= 0) nonpos_deriv_Hh = .true.
            if (predict_deriv_volh_h > 0) pos_deriv_Hh = .true.
            if (impl_part <= 0) nonpos_impl_weight = .true.
            num_h(i) = expl_part / impl_part
        enddo
        
        ! Write error and conditions that were met during the calculation:
        final_ex_T = ex_T(time=nsteps*dt)
        error = (h2T(num_h(nsteps)) - final_ex_T) / final_ex_T
        if (abs(error) > 100) error = error * ieee_value(error,ieee_positive_inf) ! Avoid large exponents, which 64-bit plotting software (e.g., Python) cannot read.
        write(error_msg,'(es20.2e4,7i5)') error,                &
                merge(1,0,nonpos_deriv_Hh .and. pos_deriv_Hh),  &
                merge(1,0,nonpos_deriv_Hh),                     &
                merge(1,0,nonpos_impl_weight),                  &
                merge(1,0,h_too_low),                           &
                merge(1,0,h_too_high)
    end function rel_T_error_after_time_steps
end module


program main
    ! Prints the error in the time-stepping scheme a series of T0 values.
    use kinds_mod
    use props_mod
    use mansol_mod, only: TMIN, TMAX, AMPLITUDE, lambda, ex_src, ex_h, ex_T
    use time_stepping_mod, only: rel_T_error_after_time_steps
    implicit none
    
    character(1000) :: error
    integer :: i
    real(WP), allocatable :: T0_vals(:)
    character(:), allocatable :: props
    real(WP) :: drho
    
    character(*), parameter :: derivhr_strategy = 'deriv_rh' ! 'rho_predictor'
    integer, parameter :: order_extrapolation = 3, &
                          order_BDF = 2

    ! Initialize the props:
    allocate(props,source='affine')
    lambda = 0.1_WP
    cp = 1.0_WP
    ! T0 = ...
    rho0 = 0.5_WP
    rho1 = 2.0_WP

    print *, '#  ======   INPUT:  ========'
    print *, '# EOS (rho <--> T): ', props
    print *, '# rho0, rho1 = ', rho0, rho1
    print *, '# cp = ', cp
    print *, '# lambda = ', lambda
    print *, '# '
    print *, '# exact_T = 0.5 + AMPLITUDE * sin(2*PI*time)'
    print *, '# AMPLITUDE = ', AMPLITUDE
    print *, '# TMIN = ', TMIN
    print *, '# TMAX = ', TMAX
    print *, '# '
    print *, '# derivhr_strategy    = ', derivhr_strategy
    print *, '# order_extrapolation = ', order_extrapolation
    print *, '# order_BDF           = ', order_BDF
    print *, '#  ========================='
    print *, '#'

    drho = rho1 - rho0
    select case(props)
    case('constant')
        print *, '# dH/dh = 0 for no T0'
        T2rho => T2rho__constant
        T2deriv_rho_T => T2deriv_rho_T__constant
    case('affine')
        print *, '# dH/dh = 0         for T0 = ', 2*TMAX - rho1 / (rho1 - rho0)
        print *, '# -h = cp/beta      for T0 = ', (drho + rho0)/drho
        print *, '# +h = cp/beta      for T0 = ', 2*TMIN - 1 - rho0/drho, '...', 2*TMAX - 1 - rho0/drho
        T2rho => T2rho__affine
        T2deriv_rho_T => T2deriv_rho_T__affine
    case('mixing')
        print *, '# dH/dh = 0         for T0 = ', -rho0 / drho
        print *, '# -h = cp/beta      for T0 = ', 2*TMIN + rho0/drho, '...', 2*TMAX + rho0/drho
        print *, '# +h = cp/beta      for T0 = ', -rho0 / drho
        T2rho => T2rho__mixing
        T2deriv_rho_T => T2deriv_rho_T__mixing
    case default
        error stop 'unknown fluid property'
    end select
    print *, '#'

    ! Loop over values of T0:
    print "(a,6(', ',a))",              &
            ' #    T0',                 &
            '  rel_error_T',            & 
            '  (non-uniq (rho*h))',     &
            '  (d/dh (rho*h) <= 0)',    &
            '  (impl_part <= 0)',       &
            '  (-h >= cp/beta)',        &
            '  (+h >= cp/beta)'
    
    ! Loop over T0 values, which is a global variable.
    T0_vals = linspace(-10.0_WP,10.0_WP,step=0.01_WP)
    do i = 1, size(T0_vals)
        T0 = T0_vals(i)
        error = rel_T_error_after_time_steps(           &
            nsteps = 2**11,                             &
            derivhr_strategy = derivhr_strategy,        &
            order_extrapolation = order_extrapolation,  &
            order_BDF = order_BDF                       &
        )
        print '(f20.2,a)', T0, trim(error)
    enddo

contains
    function linspace(x1,x2,step) result(x)
        real(WP), intent(in) :: x1, x2
        real(WP), intent(in) :: step
        real(WP), allocatable :: x(:)

        integer :: i, n

        n = 1 + floor((x2 - x1) / step)
        allocate(x(n))
        do i = 1, n
            x(i) = x1 + step * (i-1)
        enddo
    end function
end program

