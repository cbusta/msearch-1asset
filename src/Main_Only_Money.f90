program main_only_money
    ! --------------------------------------------------------------------------------------
    ! Solving the modified Lagos and Wright model when prefrences are not quasi-linear
    ! Solving the stationary distribution on a finer grid than for the bargaining problem
    ! Bivariate splines for Q and D
    ! Money balances expressed as a ratio wrt to aggregate supply of nominal money
    ! Solving for the price of money using partial adjustment relative to the excess demand 
    ! function for money
    ! Using a utility function in CM with constant Frisch elasticity of labor supply
    ! --------------------------------------------------------------------------------------
    ! This is the benchmark model version
    ! --------------------------------------------------------------------------------------
    ! This program nests the solution as a function of chi (inverse Frisch elasticity)
    ! --------------------------------------------------------------------------------------
    ! Christian Bustamante
    ! June 14, 2018, 04:49

    use lib_kind
    use lib_basic
    use lib_interp
    use lib_rwhdf5
    use get_params
    use par_pass
    use lib_solver
    use tools_fnutil
    use tools_value_fns
    use tools_distributions
    use tools_extra
    use bargaining
    use solve_steady
    implicit none

    integer           :: i,j,run_steadies
    integer,parameter :: Nchi = 6
    integer,parameter :: Nmu  = 4
    real(dp)          :: chivec_num(Nchi),muvec_num(Nmu)
    character(3)      :: chivec_str(Nchi),muvec_str(Nmu)

    ! Values for chi
    ! Need to start above 0.25 (it gets unstable below this value)
    chivec_num = [0.0, 0.25, 0.50, 1.0, 2.0, 3.0]
    chivec_str = ['000','025','050','100','200','300']

    ! Values for mu
    muvec_num = [1.01, 1.02, 1.10, 1.20]**ipyear-1.0
    muvec_str = ['010','020','100','200']

    ! Looping over values for chi
    ! Results are saved in h5 and log files
    run_steadies = 1
    if (run_steadies == 1) then
        do i = 1,Nchi
            do j = 1,Nmu
                call solve_4chimu(chivec_num(i),chivec_str(i),muvec_num(j),muvec_str(j))
            end do
        end do
    end if

end program main_only_money