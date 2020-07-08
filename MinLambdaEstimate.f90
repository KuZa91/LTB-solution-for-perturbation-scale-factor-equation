!******************************************************************************
!                                                                             *
!  PROGRAM: MinLambdaEstimate                                                 *
!                                                                             *
!  PURPOSE:  This program will simulate the evolution of semi-isotropic       *
!  secondary perturbations, in a universe filled with matter and cosmological *
!  costant. The equation that describe the evolution of the perturbation,     *
!  has been calcutated starting from a Lemaitre-Tolman-Bondi model, with the  *
!  approximation of spherical simmetry. Knowing that the solution of the      *
!  perturbation scale factor b(t) will vanish for small Lambda, this          *
!  software will use a "Divide et Impera" method to estimate the minimum      *
!  value of Lambda in order to have non-vanishing inhomogeneous perturbations *
!  across all the time range evolution of the simulated model.                *
!                                                                             *
!******************************************************************************

! Main Program
    program MinLambdaEstimate 

    implicit none

    ! Variables !

  
    integer (kind=4) :: i
    
    real (kind=8) :: t,dt,tend,TSTART,TSTOP
    real (kind=8) :: a,app,lambda,lambda_min,lambda_max,num_err
    ! the secondary scale factor b and its first derivative u will be defined as vector, so that during simulation, second slot refer to step n + 1 and first slot to step n
    real (kind=8), dimension (2) :: b,u
    
    character (len=1) :: c
    
    logical (kind=1) opt, OK
    

    ! Body of LTBMatterCostant !

    ! Defining initial values
     
    opt= .TRUE.

    do while(opt)

      dt = -1
      tend = -1
      i = 0
      t=0.
      lambda_min = 0.
      lambda_max = 0.3
      lambda = (lambda_max + lambda_min)/2.
      num_err = (lambda_max - lambda_min)/2.
        
      write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'
      
      do while(tend <= 0.)
           
            write (*,*) 'Please insert the ending time for the simulation (greater than 0) :'
            read (*,*) tend

      end do

      write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'

      do while(dt <= 0. .or. dt >=tend)
           
            write (*,*) 'Please insert the precision of the simulation (greater than 0) :'
            read (*,*) dt

      end do 

     
      write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'

      write (*,*) 'The starting value of a(0) and b(0) will be set by default to 0 !' 
      write (*,*) 'The value of the first derivative of b in 0, will be set in a way to have the scale factor b,& 
                   following the primary scale factor a, at the start of the simulation'
    
    ! Opening support files for the memorization of the results

    	open (unit = 1, file = 'LambdaMinConvergence.txt')
    	write (1,*) '################################################################'
   	  write (1,*) '#             Numerical Results of the Simulation              #'
    	write (1,*) '################################################################'
    	write (1,*) ''
    	write (1,*) '#        step        Lambda_Min        Numerical Error         #'
    	write (1,*) ''
101   format (I6,A,f10.5,A,f10.5)

      write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'

      write (*,*) 'The file of the results has been open'

    ! Writing on file the iteration at step = 0

      write (1,101) i,',',lambda,',',num_err
      write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'

    ! Beginning of the simulation

      write (*,*) 'The numerical simulation has started !'
      
    ! Inizialization of system tyme for the determination of total time elapsed
    	
      call CPU_TIME(TSTART)
    	write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'
      
    ! Estimate the value untill the percentual variation is bigger than the precision

      do while (num_err >= dt)

          ! Incrementing the index to count the number of iteration steps

          i = i + 1

          if(mod(i,10) == 0) write (*,*) 'We are on iteration number ',i
          
          ! Resetting the values of the parameters for the new simulation
          
          OK = .TRUE.
          t=0.
          a = 0.
          b(1) = 0.

          ! Estimating the value on second step of a

          b(2) = (DEXP(dt*(lambda)**(1./2.))*(1. - DEXP(-3*dt*(lambda)**(1./2.)))**(2./3.))/&
                 (DEXP((lambda)**(1./2.))*(1. - DEXP(-3*(lambda)**(1./2.)))**(2./3.))

          ! Estimating the value of the first derivative of b in order to have the ratio b/a = 1 at the start of the simulation      

          u(1) = (b(2) - b(1))/dt
    	
          ! Cycle for the time span from t = 0 to tend

          do while (t.lt.tend)
        	
              t = t + dt

              if(t>tend) t=tend

              a = (DEXP(t*(lambda)**(1./2.))*(1. - DEXP(-3.*t*(lambda)**(1./2.)))**(2./3.))/&
                  (DEXP((lambda)**(1./2.))*(1. - DEXP(-3.*(lambda)**(1./2.)))**(2./3.))
          
              
              b(2) = u(1)*dt + b(1)

              u(2) = u(1) - dt*(u(1)*((lambda)**(1./2.) + (2.*DEXP(-3.*t*(lambda)**(1./2.))&
                     *(lambda)**(1./2.))/(1. - DEXP(-3.*t*(lambda)**(1./2.)))) - 2.*b(2)&
                     *(lambda + (DEXP(-3.*t*(lambda)**(1./2.))*lambda)/(1. - &
                     DEXP(-3.*t*(lambda)**(1./2.))) + (DEXP(-t*6.*(lambda)**(1./2.))*lambda)/&
                     ((1. - DEXP(-3.*t*(lambda)**(1./2.)))**2.)) + (1.)/(2.*a))
          
              ! If b < 0 we stop the cycle and select the right part of the interval for the new split

              if(b(2) < 0) then
		          OK = .FALSE.
                  lambda_min = lambda
                  EXIT
              endif

              u(1) = u(2)
              b(1) = b(2)

   	      enddo

          if(OK) then
            lambda_max = lambda
          endif

          ! Estimating the new values for the estimated step
          
          app = lambda
          lambda =  (lambda_max + lambda_min)/2.
          num_err = (lambda_max - lambda_min)/2.

          ! Writing on file the iteration results at the actual step

          write (1,101) i,',',lambda,',',num_err


      end do

    ! Closing files

      close (unit = 1)        

    ! creation of the graphics
      CALL SYSTEM ('gnuplot -persist ConvergeToLambdaMin.plt')  
    
    ! Tempo di calcolo finale della CPU
    	call CPU_TIME (TSTOP)
    
    ! Execution Completed!
        write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'
        write (*,*) 'Simulation completed in ',TSTOP-TSTART,' seconds'

        write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'
        write (*,*) 'The estimated value of the minimum Omega Lambda is ',lambda,' +- ',num_err

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
    end program MinLambdaEstimate
