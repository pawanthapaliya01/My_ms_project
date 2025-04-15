! compile:  
! > gfortran -c parameters.f90
! > gfortran -O3 parameters.o pb.f90 -o pb.out
! execute:
! ./pb.out

        program pb
        use parameters    ! module with global parameters
        implicit none

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Define variables

        real(dbl) :: rcur(0:nmax+1)       ! current radial distances, periodic BCs: 0 = nmax, nmax+1 = 1
        real(dbl) :: lcur(1:nmax+1)       ! current back bone lengths, periodc BCs: nmax+1 = 1
        real(dbl) :: thetacur(1:nmax+1)   ! current twist angles, periodic BCs: nmax+1 = 1

        real(dbl) :: dd(nmax), aa(nmax)   ! Morse potential parameters along chain

        real(dbl) :: temp    ! thermal energy k_B T in ev
        real(dbl) :: forc    ! stretching force in ev/angstrom
        
        integer :: n0             ! selected base pair
        integer :: i, j, k, n     ! indices

        integer :: sim            ! loop index for simulations (molecules)        
        integer :: val            ! loop index for temperatures or forces
        real(dbl) :: dt, df       ! temperature and force steps

        real(dbl) :: rmed(0:nmax+1)       ! <r(n)> = mean radial distances
        real(dbl) :: thetamed(1:nmax+1)   ! <theta(n)> = mean twist angles
        real(dbl) :: lmed                 ! mean fraction of open base pairs
        real(dbl) :: f                    ! fraction of base pairs for which <r(n)> > ropen
        real(dbl) :: ymed                 ! mean displacement per base pair
        real(dbl) :: qmed                 ! mean twist angle per base pair
        real(dbl) :: hmed                 ! mean extension per base pair
        real(dbl) :: nbubmed              ! mean number of bubbles
        real(dbl) :: bubsizemed           ! mean bubble sizes averaged
        
        real(dbl) :: ll(0:200)   ! fraction of open base pairs averaged over all simulations
        real(dbl) :: ff(0:200)   ! fraction of base pairs for which <r(n)> > ropen averaged over all simulations
        real(dbl) :: yy(0:200)   ! mean displacement per base pair averaged over all simulations
        real(dbl) :: qq(0:200)   ! mean twist angle per base pair averaged over all simulations
        real(dbl) :: hh(0:200)   ! mean extension per base pair averaged over all simulations
        real(dbl) :: nn(0:200)   ! mean number of bubbles averaged over all simulations
        real(dbl) :: xx(0:200)   ! mean bubble sizes averaged over all simulations

        integer :: nopen, nbub   ! number of open base pairs and bubbles
        integer :: flag
        
        do n = 1, nmax            ! define Morse potential parameters along chain
           if (dna(n).eq.0) then        ! AT base pair
              dd(n) = dat
              aa(n) = aat
           else if (dna(n).eq.1) then   ! GC base pair
              dd(n) = dgc
              aa(n) = agc
           else
              write(*,*) 'base pair not defined'  ! should never happen
              stop
           endif
        enddo

                
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Set up simulation

        ll(:) = 0.0_dbl
        ff(:) = 0.0_dbl
        yy(:) = 0.0_dbl
        qq(:) = 0.0_dbl
        hh(:) = 0.0_dbl
        nn(:) = 0.0_dbl
        xx(:) = 0.0_dbl 
        
        dt = (tcmax - tcmin) / numt     ! temperature steps
        df = (fcmax - fcmin) / numf     ! force steps  
        
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Monte Carlo Simulation and measurement of output variables 

        do sim = 1, mol       ! loop over simulations (molecules) 

           do val = 0, numt   ! loop over temperatures                  ! comment out as needed
              temp = kb * (273.16_dbl + tcmin) + kb*dt*val              ! comment out as needed
              forc = 0.0_dbl                                            ! comment out as needed

           !do val = 0, numf   ! loop over forces                        ! comment out as needed
           !   forc = pn * fcmin + pn*df*val                             ! comment out as needed
           !   temp = kb * (273.16_dbl + 30.0_dbl)  ! 30 deg Celsius     ! comment out as needed
              
              rcur(:) = r0           ! start from equilibrium configuration at each temperature/force
              lcur(:) = l0
              thetacur(:) = theta0

              rmed(:)     = 0.0_dbl
              thetamed(:) = 0.0_dbl 
              lmed        = 0.0_dbl 
              hmed        = 0.0_dbl
              nbubmed     = 0.0_dbl
              bubsizemed  = 0.0_dbl
              call mcsweeps(rcur,lcur,thetacur,temp,forc,mctherm)    ! do mctherm mcsweeps for initial thermalization

              do i = 1, mcrun                                        ! mcrun = number of saved cofigurations
                 call mcsweeps(rcur,lcur,thetacur,temp,forc,mctime)  ! do mctime mcsweeps between saving configuration
                 rmed(:) = rmed(:) + rcur(:)                         ! sum r(n) over saved configurations
                 thetamed(:) = thetamed(:) + thetacur(:)             ! sum theta(n) over saved configurations
                 nopen = 0
                 nbub = 0
                 do n = 1, nmax                                         
                    hmed = hmed + hn(n,rcur,lcur,thetacur,flag)      ! sum base pair extensions over base pairs and saved configurations
                    if (rcur(n).gt.ropen) then         
                       nopen = nopen + 1                      ! nopen = number of open base pairs
                       if (rcur(n+1).lt.ropen) then
                          nbub = nbub + 1                     ! nbub = number of open bubbles
                       endif    
                    endif
                 enddo
                 lmed = lmed + dble(nopen)                    ! sum number of open base pairs over saved configurations                    
                 nbubmed = nbubmed + dble(nbub)               ! sum number of open bubbles over saved configurations                    
                 if (nbub.gt.0) then
                    bubsizemed = bubsizemed + dble(nopen)/dble(nbub) ! sum bubble sizes over saved configurations
                 endif
              enddo  
              rmed(:) = rmed(:) / mcrun          ! rmed(n) = <r(n)> for current simulation 
              thetamed(:) = thetamed(:) / mcrun  ! thetamed(n) = <theta(n)> for current simulation
              lmed = lmed / (mcrun*nmax)         ! lmed = mean fraction of open base pairs for current simulation
              hmed = hmed / (mcrun*nmax)         ! hmed = mean extension per base pair for current simulation
              nbubmed = nbubmed / mcrun          ! nbubmed = mean number of open bubbles for current simulation
              bubsizemed = bubsizemed / mcrun    ! bubsizemed = mean bubble sizes for current simulation
              
              f = 0.0_dbl
              ymed = 0.0_dbl
              qmed = 0.0_dbl
              do n = 1, nmax           
                 if (rmed(n).gt.ropen) then
                    f = f + 1.0_dbl
                 endif   
                 ymed = ymed + rmed(n) - r0
                 qmed = qmed + thetamed(n)
              enddo
              f = f / nmax           ! f = fraction of base pairs for which <r(n)> > ropen for current simulation
              ymed = ymed / nmax     ! ymed = mean displacement per base pair for current sim, temp
              qmed = qmed / nmax     ! qmed = average theta per base pair for current sim, temp
              
              ll(val) = (dble(sim)-1.0_dbl) / dble(sim) * ll(val) + lmed / dble(sim)        ! calculates progressively avg(lmed) (average over simulations) 
              ff(val) = (dble(sim)-1.0_dbl) / dble(sim) * ff(val) + f    / dble(sim)        ! calculates progressively avg(f)
              yy(val) = (dble(sim)-1.0_dbl) / dble(sim) * yy(val) + ymed / dble(sim)        ! calculates progressively avg(ymed)
              qq(val) = (dble(sim)-1.0_dbl) / dble(sim) * qq(val) + qmed / dble(sim)        ! calculates progressively avg(qmed)
              hh(val) = (dble(sim)-1.0_dbl) / dble(sim) * hh(val) + hmed / dble(sim)        ! calculates progressively avg(hmed)
              nn(val) = (dble(sim)-1.0_dbl) / dble(sim) * nn(val) + nbubmed / dble(sim)     ! calculates progressively avg(nbubmed)
              xx(val) = (dble(sim)-1.0_dbl) / dble(sim) * xx(val) + bubsizemed / dble(sim)  ! calculates progressively avg(bubsizemed)
              
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              ! Screen Output

              write(*,'(A17,I4,F12.3)') 'sim, T:', sim, tcmin+dt*val            ! comment out as needed 
              !write(*,'(A17,I4,F12.3)') 'sim, force:', sim, fcmin+df*val        ! comment out as needed
              write(*,'(A17,2F8.3)') 'l, f:', ll(val), ff(val)
              write(*,'(A17,3F8.3)') 'y, theta, h:', yy(val), qq(val), hh(val)
              write(*,'(A17,2F8.3)') 'nbub, bubsize:', nn(val), xx(val)
              !write(*,*) 'r(n)'                     ! output of r(n)
              !write(*,'(F12.6)') rmed(0)
              !write(*,'(10F12.6)') rmed(1:10)
              !write(*,'(10F12.6)') rmed(11:20)
              !write(*,'(10F12.6)') rmed(21:30)
              !write(*,'(10F12.6)') rmed(31:40)
              !write(*,'(10F12.6)') rmed(41:50)
              !write(*,'(10F12.6)') rmed(51:60)
              !write(*,'(F12.6)') rmed(61)
              !write(*,*) 'theta(n)'                 ! output of theta(n)
              !write(*,'(10F12.6)') thetamed(1:10)
              !write(*,'(10F12.6)') thetamed(11:20)
              !write(*,'(10F12.6)') thetamed(21:30)
              !write(*,'(10F12.6)') thetamed(31:40)
              !write(*,'(10F12.6)') thetamed(41:50)
              !write(*,'(10F12.6)') thetamed(51:60)
              !write(*,'(F12.6)') thetamed(61)
              !write(*,*)
              write(*,*)

           enddo  ! end of temperature/force loop
           write(*,*) ' '

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Output to Files

           open(14,file='l.dat',position='append')
               write(14,*) ll(:)
               write(14,*)
           close(14)
           open(15,file='f.dat',position='append')
               write(15,*) ff(:)
               write(15,*)
           close(15)
           open(16,file='y.dat',position='append')
               write(16,*) yy(:)
               write(16,*)
           close(16)
           open(17,file='theta.dat',position='append')
               write(17,*) qq(:)
               write(17,*)
           close(17)
           open(18,file='h.dat',position='append')
               write(18,*) hh(:)
               write(18,*)
           close(18)
           open(19,file='nbub.dat',position='append')
               write(19,*) nn(:)
               write(19,*)
           close(19)
           open(20,file='bubsize.dat',position='append')
               write(20,*) xx(:)
               write(20,*)
           close(20)    

        enddo ! end of loop over simulations (molecules) 


        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutines

        contains        
        
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Monte Carlo sweeps

        subroutine mcsweeps(rcur,lcur,thetacur,temp,forc,nsteps)   ! perform nsteps Monte Carlo sweeps (one sweep = nmax trial moves)
        implicit none

        real(dbl), intent(inout) :: rcur(0:nmax+1)
        real(dbl), intent(inout) :: lcur(1:nmax+1)
        real(dbl), intent(inout) :: thetacur(1:nmax+1)

        real(dbl), intent(in) :: temp    ! thermal energy k_B T in eV
        real(dbl), intent(in) :: forc    ! stretching force in eV / angstrom
        integer, intent(in) :: nsteps    

        real(dbl) :: rtrl(0:nmax+1)
        real(dbl) :: ltrl(1:nmax+1)
        real(dbl) :: thetatrl(1:nmax+1)

        integer :: err, flag  
        integer :: i, j
        real(dbl) :: ran              ! random number in [0,1)
        integer :: n0                 ! selected base pair       
        real(dbl) :: ucur, utrl       ! energies of current and trial 
        real(dbl) :: sdopen
        real(dbl) :: sdbackbone
        real(dbl) :: sdtwist
        
        integer :: rnum, lnum, tnum
        integer :: racc, lacc, tacc
        integer :: mtype

        sdopen = sqrt(0.5_dbl*temp/ss)   ! standard deviation of the Gaussian for new base pair separation
        sdbackbone = sqrt(temp/bb)       ! standard deviation of the Gaussian for backbone stretching
        sdtwist = sqrt(temp/tt)          ! standard deviation of the Gaussian for twisting

        do i = 1, nsteps
           do j = 1, nmax

              rtrl = rcur
              ltrl = lcur
              thetatrl = thetacur

              ! pick random base pair n0 in (1,...,nmax)
              call random_number(ran)    ! generates random number ran in [0,1)
              n0 = floor(ran*nmax) + 1    
              
              ! make a trial move for base pair n0 
              call random_number(ran)
              if ((ran.lt.0.4_dbl)) then
                 call openbasepair(n0,rtrl,sdopen,err) 
                 if (err.eq.1) cycle     ! trial move not allowed: reject trial move and continue with current
              else if (ran.lt.0.8_dbl) then
                 call twisting(n0,thetatrl,sdtwist,err)
                 if (err.eq.1) cycle     ! trial move not allowed: reject trial move and continue with current
              else
                 call backbonestretching(n0,ltrl,sdbackbone)
              endif

              ! calculate energy of current/trial configurations that depends on base pair n0

              ucur = morse(n0,rcur,dd,aa)   ! Morse potential for n0
              utrl = morse(n0,rtrl,dd,aa)

              ! periodc BCs: if n0=1 then bond is bw 1,0 corresponding to 1,nmax
                 ucur = ucur + bond(n0,rcur,lcur,thetacur,forc,flag)
                 utrl = utrl + bond(n0,rtrl,ltrl,thetatrl,forc,flag)
                 if (flag.eq.1) cycle     ! trial move not allowed: reject trial move and continue with current

              ! periodc BCs: if n0=nmax+1 then bond is bw nmax+1,nmax corresponding to 1,nmax 
                 ucur = ucur + bond(n0+1,rcur,lcur,thetacur,forc,flag)
                 utrl = utrl + bond(n0+1,rtrl,ltrl,thetatrl,forc,flag)
                 if (flag.eq.1) cycle
                 
              ! accept or reject trial configuration using Metropolis criterion
              call random_number(ran) 
              if (exp((ucur-utrl)/temp).lt.ran) then
                 cycle     ! reject trial move, continue with previous current
              endif 

              ! accept trial move as new current
              rcur = rtrl               
              rcur(nmax+1) = rcur(1)  ! periodic BCs
              rcur(0) = rcur(nmax)    ! periodic BCs

              lcur = ltrl
              lcur(nmax+1) = lcur(1)  ! periodic BCs    

              thetacur = thetatrl
              thetacur(nmax+1) = thetacur(1)  ! periodic BCs

           enddo
        enddo
        
        end subroutine mcsweeps


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Energy

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Morse potential for base pair n0

        function morse(n0,r,dd,aa)   
        implicit none
        integer, intent(in) :: n0
        real(dbl), intent(in) :: r(0:nmax+1)
        real(dbl), intent(in) :: dd(nmax), aa(nmax)
        real(dbl) :: morse
        real(dbl) :: mp
        
        mp = exp(-aa(n0)*(r(n0)-r0)) - 1.0_dbl
        morse = dd(n0) * mp**2
        end function morse

        
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Total energy for bond between n0, n0-1        

        function bond(n0,r,l,theta,forc,flag)
        implicit none
        integer, intent(in) :: n0
        real(dbl), intent(in) :: r(0:nmax+1), l(1:nmax+1), theta(1:nmax+1)
        real(dbl), intent(in) :: forc
        integer, intent(out)  :: flag
        real(dbl) :: bond
          
        bond = stacking(n0,r) + backbone(n0,l) + &
             & elastic(n0,r,l,theta,flag) + exforce(n0,r,l,theta,forc,flag)

        ! bond = stacking(n0,r)

        end function bond


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Stacking interaction between n0, n0-1

        function stacking(n0,r)   
        implicit none
        integer, intent(in) :: n0
        real(dbl), intent(in) :: r(0:nmax+1)
        real(dbl) :: stacking
        real(dbl) :: dr, d2r

        dr  = r(n0) - r(n0-1)
        d2r = r(n0) + r(n0-1) - 2.0_dbl * r0
        stacking = ss/2.0_dbl * (1.0_dbl + 2.0_dbl * exp(-bstack*d2r)) * dr**2    ! as in Ares
        end function stacking

        
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Backbone elastic energy between n0, n0-1

        function backbone(n0,l)   
        implicit none
        integer, intent(in) :: n0
        real(dbl), intent(in) :: l(1:nmax+1)
        real(dbl) :: backbone
        real(dbl) :: dl
        
        dl = l(n0) - l0
        backbone = bb/2.0_dbl * dl**2
        end function backbone


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Elastic coupling of twist and base pair opening between n0, n0-1

        function elastic(n0,r,l,theta,flag)   
        implicit none
        integer, intent(in) :: n0
        real(dbl), intent(in) :: r(0:nmax+1), l(1:nmax+1), theta(1:nmax+1)
        integer, intent(out)  :: flag
        real(dbl) :: elastic
        real(dbl) :: dh
        real(dbl) :: d2r

        dh = hn(n0,r,l,theta,flag) - h0 
        d2r = r(n0) + r(n0-1) - 2.0_dbl * r0
        elastic = cc/2.0_dbl * dh**2 * exp(-bel*d2r)
        end function elastic        

        
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Energy due to external stretching force between n0, n0-1
 
        function exforce(n0,r,l,theta,forc,flag)
        implicit none
        integer, intent(in) :: n0
        real(dbl), intent(in) :: r(0:nmax+1), l(1:nmax+1), theta(1:nmax+1)
        real(dbl), intent(in) :: forc
        integer, intent(out)  :: flag
        real(dbl) :: exforce

        exforce = - forc * hn(n0,r,l,theta,flag)
        end function exforce 


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Distance between base pair planes nn, nn-1

        function hn(n0,r,l,theta,flag)
        implicit none
        integer, intent(in)   :: n0
        real(dbl), intent(in) :: r(0:nmax+1), l(1:nmax+1), theta(1:nmax+1)
        integer, intent(out)  :: flag
        real(dbl) :: hn 
        real(dbl) :: r1, r2, ln, thetan
        real(dbl) :: arg
        
        r1 = r(n0)
        r2 = r(n0-1)
        ln = l(n0)
        thetan = theta(n0)
        arg = ln**2 - r1**2 - r2**2 + 2.0_dbl*r1*r2*cos(thetan)
        if (arg.lt.0) then
           flag=1
           hn = 0.0_dbl 
        else 
           flag=0
           hn = sqrt(arg)
        endif
        end function hn        


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Trial moves        

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Trial move for base pair opening

        subroutine openbasepair(n0,r,sdopen,err)

        implicit none
        integer, intent(in) :: n0
        real(dbl), intent(inout) :: r(0:nmax+1)
        real(dbl), intent(in) :: sdopen   ! standard deviation of the Gaussian proposing new base pair separation
        integer, intent(out) :: err 
        real(dbl) :: xi
        integer :: n
        
        xi = ran_gaussian()              ! random number drawn from Gaussian distribution with mean 0 and variance 1
        r(n0) = r(n0) + sdopen * xi      ! displace base pair n0
    
        err = 0
        r(nmax+1) = r(1)
        r(0) = r(nmax) 
        if (minval(r).ge.ropen) then     ! dsDNA condition: at least one base pair must be closed
           err = 1                     
        endif
                        
        end subroutine openbasepair

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Trial move for backbone stretching
        
        subroutine backbonestretching(n0,l,sdbackbone)

        implicit none
        integer, intent(in) :: n0
        real(dbl), intent(inout) :: l(1:nmax+1)
        real(dbl), intent(in) :: sdbackbone
        real(dbl) :: xi

        xi = ran_gaussian()                ! random number drawn from Gaussian distribution with mean 0 and variance 1
        l(n0) = l(n0) + sdbackbone * xi    ! modify backbone length
    
        end subroutine backbonestretching

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Trial move for twisting

        subroutine twisting(n0,theta,sdtwist,err)

        implicit none
        integer, intent(in) :: n0
        real(dbl), intent(inout) :: theta(1:nmax+1)
        real(dbl), intent(in) :: sdtwist
        integer, intent(out) :: err 
        real(dbl) :: xi

        xi = ran_gaussian()                     ! random number drawn from Gaussian distribution with mean 0 and variance 1
        theta(n0) = theta(n0) + sdtwist * xi    ! modify twist angle

        err = 0
        if( (theta(n0).lt.thetamin) .or. (theta(n0).gt.thetamax) ) then  ! theta outside of [thetamin, thetamax]
           err = 1 
        endif    
        end subroutine twisting
        
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Generate random number drawn from Gaussian distribution with mean 0 and variance 1 using the Box-Muller method

        function ran_gaussian()
        implicit none

        real(dbl) :: ran_gaussian
        real(dbl) :: ran
        real(dbl) :: r1
        real(dbl) :: r2

        call random_number(ran)
        r1 = ran
        call random_number(ran)
        r2 = ran
        ran_gaussian = sqrt(-2.0_dbl * log(r1)) * cos(twopi * r2)
        end function ran_gaussian
        
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        

        end program pb


