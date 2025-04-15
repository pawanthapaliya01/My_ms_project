        module parameters
        implicit none

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Chain sequence and parameters for simulation range - edit as needed

        ! sequence L33B9: CCGCCAGCGGCCTTTACTAAAGGCCGCTGCGCC
        integer, parameter :: nmax = 33       ! number of base pairs
        integer, parameter, dimension(nmax) :: dna = [1,1,1,1,1,0,1,1,1,1,1,1,0,0,0,0,1, &
             & 0,0,0,0,1,1,1,1,1,1,0,1,1,1,1,1]  ! 0=AT, 1=GC

        ! Energies are expressed in eV and lengths in angstrom

        integer, parameter   :: dbl = kind(1.d0)       ! dbl makes real numbers double precision

        integer, parameter :: mctherm = 1000000
        integer, parameter :: mctime  = 100
        integer, parameter :: mcrun   = 10000
        integer, parameter :: mol     = 100

        real(dbl), parameter :: tcmin =  30_dbl    ! minimum temperature in Celcius
        real(dbl), parameter :: tcmax = 100_dbl    ! maximum temperature in Celsius
        integer, parameter :: numt = 35            ! number of temperatures in temperature loop 

        real(dbl), parameter :: fcmin =  0_dbl     ! minimum force in pN
        real(dbl), parameter :: fcmax = 40_dbl     ! maximum force in pN
        integer, parameter :: numf = 20            ! number of forces in force loop


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Fixed model parameters

        real(dbl), parameter :: pi = 3.141592653589793_dbl
        real(dbl), parameter :: twopi = 6.283185307179586_dbl

        real(dbl), parameter :: dat = 0.058_dbl        ! strength of Morse potential for AT in eV  
        real(dbl), parameter :: dgc = 1.5_dbl * dat    ! strength of Morse potential for GC
        real(dbl), parameter :: aat = 8.4_dbl          ! 1/range of Morse potential for AT in 1/angstrom
        real(dbl), parameter :: agc = 13.8_dbl         ! 1/range of Morse potential for GC

        real(dbl), parameter :: evj = 1.602177e-19_dbl        ! one electron volt in joules
        real(dbl), parameter :: angstrom = 1e-10_dbl          ! one angstrom in meters
        real(dbl), parameter :: kb = 1.380649e-23_dbl / evj   ! Boltzmann constant in eV/K
        real(dbl), parameter :: pn = 1.0e-22_dbl / evj        ! one pN in eV/angstrom

        real(dbl), parameter :: h0 = 2.9_dbl                  ! distance between base pair planes at T=0 in angstrom
        real(dbl), parameter :: r0 = 10.0_dbl                 ! equilibrium radial distance of bases in angstrom
        real(dbl), parameter :: theta0 = twopi / 10.4_dbl     ! equilibrium twist angle 0.604 for B-DNA in rad 
        real(dbl), parameter :: thetamin = -0.1_dbl           ! minimum twist angle
        real(dbl), parameter :: thetamax =  0.7_dbl           ! maximum twist angle  
        
        real(dbl), parameter :: ropen = 10.25_dbl    ! threshold for radial distance to consider a base pair open in angstrom

        real(dbl), parameter :: ss = 0.1_dbl         ! stacking interaction in ev/angstrom^2 as in Ares
        real(dbl), parameter :: bstack = 0.7_dbl     ! 1/range of stacking interaction in 1/angstrom as in Ares
        
        real(dbl), parameter :: bb = 0.128_dbl    ! elastic constant for backbone stretching in ev / angstrom^2
        real(dbl), parameter :: l0 = 6.853_dbl    ! equilibrium backbone length in angstrom

        real(dbl), parameter :: cc  = 0.028_dbl   ! elastic energy constant in ev/angstrom^2 as in Cocco
        real(dbl), parameter :: bel = 0.10_dbl    ! 1/range of elastic energy in 1/angstrom  

        real(dbl), parameter :: tt = 5.64_dbl     ! energy constant for twist in eV/rad^2
        
        end module parameters
