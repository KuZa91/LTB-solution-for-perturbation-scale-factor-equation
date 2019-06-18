!*****************************************************************************
!                                                                            *
!  PROGRAM: LTBMatterCostant                                                 *
!                                                                            *
!  PURPOSE:  This program will simulate the evolution of semi-isotropic      *
!  secondary perturbations, in a universe filled with matter and cosmological*
!  costant. The equation that describe the evolution of the perturbation,    *
!  has been calcutated starting from a Lemaitre-Tolman-Bondi model, with the *
!  approximation of spherical simmetry. The software will evaluate the       *
!  evolution of the time dependance b(t) of the secondary perturbations ,    *
!  by integrating the equation with a explicit Eulero method.                *
!  As the equation will have some critic point around t=0, the time step     *
!  of the integraion will be variable on each different step, so that the    *
!  precision of the simulation will be the same on each step.                *
!                                                                            *
!*****************************************************************************

! Main Program
    program LTBMatterCostant 

    implicit none

    ! Variables !

  
    integer (kind=4) :: i
    
    real (kind=8) :: prec,t,dt,tend,tapp,TSTART,TSTOP
    real (kind=8) :: a,eta,lambda,var
    ! the secondary scale factor b and its first derivative u will be defined as vector, so that during simulation, second slot refer to step n + 1 and first slot to step n
    real (kind=8), dimension (2) :: b,u
    
    character (len=1) :: c
    
    logical (kind=1) opt,ok
    

    ! Body of LTBMatterCostant !

    ! Defining initial values
     
    opt= .TRUE.

    do while(opt)

        t=0.
        lambda = -1.
        tend = -1.
        dt = 0.001
        prec = -1


    	write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'

      do while(lambda <= 0. .or. lambda > 1.)
           
            write (*,*) 'Please insert the value of Lambda (between 0 and 1) :'
    	      read (*,*) lambda

      end do    

        
      write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'
      
      do while(tend <= 0.)
           
            write (*,*) 'Please insert the ending time for the simulation (greater than 0) :'
            read (*,*) tend

      end do 


      write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'
      
      do while(prec <= 0. .or. prec >= 0.5)
           
            write (*,*) 'Please insert the precision of the simulation (between 0 and 0.5) :'
            read (*,*) prec

      end do  

     
      write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'

      write (*,*) 'The value of a(0) and b(0) will be set by default to 0, the value of u(0) instead will be set to 1'

      a = 0.
      b(1) = 0.
      u(1) = 1.
      eta = 0.

    

    ! Opening support files for the memorization of the results

    	open (unit = 1, file = 'RisultatiLTBMatterCostant.txt')
    	write (1,*) '################################################################'
   	  write (1,*) '#             Risultati numerici della simulazione             #'
    	write (1,*) '################################################################'
    	write (1,*) ''
    	write (1,*) '#        t           a          b            eta               #'
    	write (1,*) ''
101   format (4f20.15)

      write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'

      write (*,*) 'The file of the results has been open'

    ! Writing on file the value of the variables at t = 0

      write (1,101) t,a,b(1),eta

      write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'

    ! Beginning of the simulation

      write (*,*) 'The numerical simulation has started !'
      
    ! Inizialization of system tyme for the determination of total time elapsed
    	
      call CPU_TIME(TSTART)
    	write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'
    	
      ! First cycle, it let times go from t = 0 to tend
      do while (t.lt.tend)
        	
          ! the index i will break the cycle if the simulation doesn't converge 
          i = 0 
          ok = .TRUE.
          dt = 0.001
 
          ! Second cycle, used to set the value of dt that respect the desired precision
		      do while (ok)
              
              i = i +1
              if(i == 5) then
                  write (*,*) 'There was an convergence error!'
                  goto 999
              endif

              tapp = t + dt
              if(tapp>tend) tapp=tend

              a = (DEXP(tapp*(lambda/3.)**(1./2.))*(1. - DEXP(-tapp*(3.*lambda)**(1./2.)))**(2./3.))/&
              (DEXP((lambda/3.)**(1./2.))*(1. - DEXP(-(3.*lambda)**(1./2.)))**(2./3.))

              b(2) = u(1)*dt + b(1)

              u(2) = u(1) - dt*(u(1)*((lambda/3.)**(1./2.) + (2.*DEXP(-tapp*(3.*lambda)**(1./2.))*((lambda)/3.)**(1./2.))/&
              (1. - DEXP(-tapp*(3.*lambda)**(1./2.)))) - b(2)*((2.*lambda)/(3.) + ((2./3.)*DEXP(-tapp*(3.*lambda)**(1./2.))*lambda)/&
              ((1.) - DEXP(-tapp*(3.*lambda)**(1./2.))) + ((2./3.)*DEXP(-tapp*2.*(3.*lambda)**(1./2.))*lambda)/&
              (((1.) - DEXP(-tapp*(3.*lambda)**(1./2.)))**2.)) + (1.)/(2.*a))

              var = (u(2) - u(1))/(u(1))

              if(var >= prec) then

                  dt = ABS((prec * u(1))/(2*(u(1)*((lambda/3.)**(1./2.) + (2.*DEXP(-tapp*(3.*lambda)**(1./2.))*&
                  ((lambda)/3.)**(1./2.))/((1.) - DEXP(-tapp*(3.*lambda)**(1./2.)))) - b(2)*((2.*lambda)/(3.) &
                  +  ((2./3.)*DEXP(-tapp*(3.*lambda)**(1./2.))*lambda)/((1.) - DEXP(-tapp*(3.*lambda)**(1./2.))) +&
                  ((2./3.)*DEXP(-tapp*2.*(3.*lambda)**(1./2.))*lambda)/(((1.) - DEXP(-tapp*(3.*lambda)**(1./2.)))**2.)) &
                  + (1.)/(2.*a))))
              else

                  ok = .FALSE.

              endif

          end do

          t = t + dt

          eta = (b(2)**2)/(a**2)

          write (1,101) t,a,b(2),eta

          u(1) = u(2)
          b(1) = b(2)
      end do

    ! Closing files

      close (unit = 1)        

    ! creation of the graphics
      CALL SYSTEM ('gnuplot -persist LTBMatterCostanta.plt') 
      CALL SYSTEM ('gnuplot -persist LTBMatterCostantb.plt')
      CALL SYSTEM ('gnuplot -persist LTBMatterCostanteta.plt')     
    
    ! Tempo di calcolo finale della CPU
    	call CPU_TIME (TSTOP)
    
    ! Execution Completed!
        write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'
        write (*,*) 'Simulation completed in ',TSTOP-TSTART,' seconds'

999   	write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'
        write (*,*) 'Do you wish to relaunch the program (y/n)?'
    	  read (*,*) c
    	  if(c/='y'.and.c/='n'.and.c/='Y'.and.c/='N') then
		          write (*,*) 'Wrong Immission!'
        	    goto 999
   	    endif
        
        write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'
        if(c=='n'.or.c=='N') then
		        opt = .FALSE.
   	    endif
    end do  

!

    ! End Program !
    stop
    end program LTBMatterCostant
