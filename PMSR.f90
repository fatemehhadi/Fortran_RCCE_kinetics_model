
program PaSR

   use RCCE_Mod
   use Matrix_Mod
   use mod_hectars
   use LINPACK_Mod
   use LINPACK_Interface
   use cpreactor_test
   implicit none
!------------------------------------------------------------------------------------------
 !  explicit interface to the spread function
    interface
        function spread(minf,maxf)
           integer, intent(in) :: minf, maxf
        end function
    end interface
!------------------------------------------------------------------------------------------  
 ! Define  variables 
   integer, parameter :: Np=100  ! number of particles
   integer, parameter :: nSpecs=26
   integer, parameter :: nCnstrnts=12    
   integer :: nDt                 ! number of time steps
   integer :: iphiFirst  ! First index to phi array
   integer :: iphiLast   ! Last index to phi array
   integer :: imixFrac   ! index to mixture fraction
   integer :: nphi       ! Lenght of composition vector
   integer :: nvar       ! Number of particle variables
   integer :: nIn        ! Number of particle variables  
   integer :: mOut   
   integer :: solverInt   
   integer :: restartITime,restartTime
   integer :: i,j,k,kk,iTime, iTimeCount
   integer, allocatable, dimension(:)  :: pOut   
   integer, allocatable, dimension(:,:)  :: pOutArray      
   real(8) :: deltaT
   real(8) :: aMtrx(nCnstrnts,nSpecs)
   real(8) :: a11aMtrx(nCnstrnts,nSpecs)
   real(8) :: RHSVec(Np,nCnstrnts+3)
   real(8) :: massFrcPrtl(Np,nSpecs)
   real(8) :: massFrc2Prtl(Np,nSpecs)   
   real(8) :: gammaPrtl(Np,nCnstrnts)
   real(8) :: press(Np)
   real(8) :: rho(Np)    
   real(8) :: Temp(Np)
   real(8) :: Cnstrnts(Np,nCnstrnts)
   real(8) :: CnstrntsAvg(nCnstrnts)   
   real(8) :: hAvg, xiAvg
   real(8) :: t, tout
   real(8) :: t_dassl(Np), tout_dassl(Np)
   real(8) :: MWs(nSpecs)
   real(8) :: massFrac(nSpecs)
   real(8) :: totalMass
   real(8) :: hF, hO2
   real(8) :: Tmix
   real(8) :: Tres 
   real(8) :: t1(Np),t2(Np)  
   real(8) :: massSum
   real(8), dimension(:,:), allocatable :: phiIn
   real(8), dimension(:,:), allocatable :: phi
   real(8), allocatable, dimension(:,:) :: average
   real(8), allocatable :: averageTime(:)
   real(8), allocatable, dimension(:,:)  :: mInn
   type(hectars) :: htr         !hectars variable
   character(20) :: solver
   character(80) :: fmt1, fmt2, fmt3, fmt4, fmt5, fmt6, fmt7, fmt8, fmt9, fmt10, fmt11 
   logical :: restart
   logical :: IncomingPrtl
      
 ! Mohammad RCCE initialization
   real(8) :: atol
   real(8) :: rtol
   real(8) :: anini
   integer :: nfuel
   integer :: nf
   integer :: noxy
   integer :: anox
   integer :: nn2
   real(8) :: an2
   real(8) :: tempdif
   real(8) :: timelimit
   integer :: xnnout
   integer :: xnout
   real(8) :: toutmax
   real(8) :: mmcount
   integer :: itrlim
   real(8) :: ddtt
   real(8) :: tmin 
   real(8) :: tmax
   integer :: tstep
   integer, parameter :: LRW = 500
   integer, parameter :: LIW = 100
   real(8) :: RWORK(Np,LRW)
   integer :: IWORK(Np,LIW)   
!  nc
!  temp
!  pres   
!  dt
!  dt2   
!  gamma
!  a
!  aij

   call random_seed
   call init_MohammadRCCE
   call formats(nCnstrnts,nSpecs,Np,fmt1,fmt2,fmt3,fmt4,fmt5,&
        fmt6,fmt7,fmt8,fmt9,fmt10,fmt11)   
   call header_main

   restart = .false.
   restartITime = 1  ! Assign this number to one if restart = .false.
   restartTime = 0.0d0 ! Assign this number to zero if restart = .false.
      
   solverInt = 5
   select case (solverInt)
      case (1)
         solver = "Euler"
      case (2)
         solver = "rungkutta4"
      case (3)
         solver = "dvodeSolver"
      case (4)
         solver = "dasslSolver" 
      case (5)
         solver = "dlsodaSolver"   
      case (6)
         solver = "MohammadDASSL"              
   end select
   write(*,'(/A,1x,A)') "Your solver is", solver
      
 ! Assigning pressure. P = 1.0e6 ergs/cm3 = 1.013250d5 N/m2 = 1 atm
   !! press = 1.013250d6
   !! 3.
   press(1:Np) = 1.013250d5 ![N/m2]
   
   call init_CH4_JUMK_CNF2009(Np,nvar,nphi,nIn,iphiFirst,iphiLast,iMixFrac,hF,hO2)

   deltaT = 1.0d-5 ! Time step [s]
   nDt = 300        
   t = restartTime
   tout = t + deltaT
   Tmix = 1.d-2
   Tres = 1.d-4 
   mOut=int(Np*deltaT/Tres)
   t_dassl(1:Np)=t
   tout_dassl(1:Np)=tout
   
   allocate(average(nvar,nDt))
   allocate(averageTime(nDt))
   allocate(pOut(mOut))
   allocate(mInn(nvar, mOut))  
   allocate(pOutArray(mOut, nDt)) 
   
   open(unit=100, file = 'CPUTime.txt')
   open(unit=4,   file = 'pOutArray.txt')
   open(unit=14,  file = 'Case5.txt')
   open(unit=600, file = 'massSumMainFileError.txt')
   open(unit=500, file = 'massSumMainFile.txt')
   
   averageTime(:) = 0
   do k=1,nDt
      read(4,*) pOutArray(:,k)
   enddo 
   close(4)   
   do j=1,nSpecs
      read(14,'(1x,d25.18,1x)') massFrac(j)
   enddo 
   read(14,'(1x,d25.18,1x)') Temp(Np)
   close(14)

 ! gets molecular weight of species [kg/kmol]
   !!call get_MWs(cpr, MWs)
   !! 4.
   call getMolecularWeights(htr, MWs)
 
 ! Set mass fractions
   do k=1,Np
      !!call set_YT(cpr, phi(:,k), massFrac, Temperature+2.d0*k)
      !! 5.
      do j=1,nSpecs
         phi(j,k) = massFrac(j)
      enddo
      Temp(k) = Temp(Np) + 2.d0*k
      call setState_TPY(htr,Temp(k),press(k),phi(1:nphi-1,k))
      phi(nphi,k)=enthalpy_mass(htr)
      !! write(*,'(/d25.17)') enthalpy_mass(htr)
   enddo
   
   if (restart) then
      call read_phi(restart,restartITime,restartTime,Np,nSpecs,nphi,phi)
   endif
    
 ! RCCE initial condition calculation
   !! 6.
   call RCCEInit(.false.,iTime,Np,nCnstrnts,nSpecs,nphi,htr,phi,Temp,press,aMtrx,&
        massFrcPrtl,massFrc2Prtl,gammaPrtl)
  
   !! 7.
   do k=1,Np
   call Calc_rhoTemp(k,Np,nphi,htr,phi(:,k),press(k),rho(k),Temp(k))
   enddo
   
 ! z (mixture fraction) = (sYF - YO + Y0O)/(sY0F + Y0O) 
 ! s = 2MO/MF
   do k=1,Np
      phi(nvar,k) = (2.* MWs(2)/MWs(7) * phi(7,k) - phi(2,k) + 1.D0)/ &
                    (2.* MWs(2)/MWs(7) * 1.D0 + 1.D0)
   enddo 
      
 ! RCCE calculation
   do iTime=restartITime,nDt
       
      !!call header_loop(iTime)

      !! 8.
      call save_phi(restart,restartITime,iTime,t,deltaT,&
      Np,nphi,htr,phi,nvar,Temp,nDt,average)

      !! 9.
      call events(Np, mOut, nvar, htr, phi, hF, hO2, press, pOut, &
      mInn, iTime, nDt, pOutArray, Temp, nCnstrnts, gammaPrtl)
      
      do k=1,Np
      call Calc_rhoTemp(k,Np,nphi,htr,phi(:,k),press(k),rho(k),Temp(k))
      enddo
      
      !! 12. 
      call RHS(iTime, Np, nCnstrnts, nSpecs, aMtrx, a11aMtrx, &
           htr, phi, nphi, massFrcPrtl, massFrc2Prtl, Temp, rho, &
           press, Cnstrnts, CnstrntsAvg, hAvg, xiAvg, RHSVec,Tmix)

      do k=1,Np
     
      select case (solver)

         case ("Euler")
            !! 13.
            call Euler(iTime, Np, nCnstrnts, nSpecs, aMtrx, htr, &
                 phi, nphi,massFrc2Prtl, CnstrntsAvg, hAvg, xiAvg, &
                 gammaPrtl, rho, Temp, press(k), deltaT, RHSVec)
         
         case ("rungkutta4")
            !! 14.
            call rungkutta4(iTime, Np, nCnstrnts, nSpecs, aMtrx, htr, &
                 phi, nphi, massFrc2Prtl, CnstrntsAvg, hAvg, xiAvg, &
                 gammaPrtl, rho, Temp, press(k), deltaT, RHSVec)
                 
         case ("dvodeSolver")
            call cpu_time ( t1(k) )
            !! 15.
            call dvodeSolver(k,iTime, Np, nCnstrnts, nSpecs, aMtrx, htr, &
                 phi(:,k), nphi, massFrc2Prtl(k,:), CnstrntsAvg, hAvg, xiAvg, &
                 gammaPrtl(k,:), rho(k), Temp(k), press(k), deltaT, RHSVec(k,:),Tmix)             
            
            call cpu_time ( t2(k) )
            write (100,*) iTime,k,(t2(k)-t1(k)) 
                       
         case ("dasslSolver")
            call cpu_time ( t1(k) )
            !! 16.
            call dasslSolver(k,iTime, Np, nCnstrnts, nSpecs, aMtrx, htr, &
                 phi(:,k), nphi, massFrc2Prtl(k,:), CnstrntsAvg, hAvg, xiAvg, &
                 gammaPrtl(k,:), rho(k), Temp(k), press(k), deltaT, RHSVec(k,:),Tmix)            
            
            call cpu_time ( t2(k) )
            write (100,*) 'Elapsed CPU time = ', t2 - t1
                 
         case ("dlsodaSolver")
            call cpu_time ( t1(k) )
            !! 17.
            call dlsodaSolver(k,iTime, Np, nCnstrnts, nSpecs, aMtrx, htr, &
                 phi(:,k), nphi, massFrc2Prtl(k,:), CnstrntsAvg, hAvg, xiAvg, &
                 gammaPrtl(k,:), rho(k), Temp(k), press(k), deltaT, RHSVec(k,:),Tmix) 
            
            call cpu_time ( t2(k) )
            write (100,*) iTime,k,(t2(k)-t1(k))
                
         case("MohammadDASSL")
            call cpu_time ( t1(k) )
            !RWORK(k,:) = 0.d0 ! Remove in the continuation mode
            !IWORK(k,:) = 0 ! Remove in the continuation mode
            call Rate_Controlled_Constrained_Equilibrium(mOut,pOutArray(:,iTime),k,t_dassl(k),tout_dassl(k),atol,rtol, &
	             Temp(k),1.0d0,tempdif,timelimit, &
                 deltaT,deltaT,ddtt,tmin,tmax, &
                 gammaPrtl(k,:),a11aMtrx,aMtrx,phi(1:nphi-1,k),LRW,LIW,RWORK(k,:),IWORK(k,:),iTime)
 
            call cpu_time ( t2(k) )
            write (100,*) iTime,k,(t2(k)-t1(k))
              
            !!call set_YT(cpr, phi(:,k), massFrc2Prtl(k,:), Temp(k)) 
            !! 18. 
            call setState_TPY(htr,Temp(k),press(k),phi(1:nphi-1,k))
            massFrc2Prtl(k,:)=phi(1:nphi-1,k)
      end select
            
      !! 19.
      if (mod(iTime,10) == 0) then
         write (*,'(A, I10, A, I10, A, d25.18, A, d17.10)') "iTime: ", &
         iTime, " prtlNo = ", k, "  Temp: ", Temp(k), "  deltaT: ", deltaT 
      endif
      
      call check_massSum     

      enddo !! do k=1,Np

      averageTime(iTime)=sum(t2(:)-t1(:))/Np
      write ( 100, * ) 'Average elapsed CPU time in iTime ', iTime, averageTime(iTime)
      write ( 100, * ) 'Average elapsed CPU time in total ', sum(averageTime(:))/iTime
      
      if ((mod(iTime,1000) == 0) .OR. (iTime == nDt)) then
         call writeRestart_phi(iTime,t,Np,nSpecs,nphi,phi)
      endif
               
      if (iTime == nDt) then
         !! 20.
         call save_phi(restart,restartITime,iTime,t,deltaT,Np,nphi,htr,phi,nvar,Temp,nDt,average)
      endif
           
   enddo
   call cpu_time (t2(Np))
   write (*,*) 'Elapsed CPU time = ', t2(Np) 

   
  !===================================================================== 
Contains
  subroutine check_massSum
      massSum = 0
      do j=1,nSpecs
         massSum = massSum + phi(j,k)  
      end do
    
      write(500,'(/i5,i5,d25.18)') iTime, k, massSum

      if (abs(massSum-1) > 0.01) then
         IncomingPrtl = .false.
         do kk=1,mOut
            if (k == pOutArray(kk,iTime)) then
               IncomingPrtl = .true.  
            endif
         enddo
         if (IncomingPrtl == .true.) then
            write(600,'(/i5,i5,d25.18,1x,A)') iTime, k, massSum, "  Incoming particle"
         else
            write(600,'(/i5,i5,d25.18)') iTime, k, massSum
         endif
      endif  
  end subroutine check_massSum
  !===================================================================== 
  subroutine init_MohammadRCCE

 ! Mohammad RCCE initialization
   nfuel = 3
   nf = 2
   noxy = 2
   anox = 1
   nn2 = 1
   xnnout = 100
   xnout = 100
   itrlim = 10
   tstep = 1000
   atol = 1.d-9
   rtol = 1.d-9
   anini = 1.d-30
   an2 = 1.d-25
   tempdif = 0.5d0 
   timelimit = 1.8627d0
   toutmax = 1.d+10
   mmcount = 0.d0
   ddtt = 1.1d0
   tmin = 1.d-3 
   tmax = 1.d-2
!  nc
!  temp
!  pres
!  dt
!  dt2   
!  gamma
!  a
!  aij

  end subroutine init_MohammadRCCE
  !=====================================================================
  !!  Sub: Init_CH4_JUMK_CNF2009
  !!       Initialization for methane oxidation
  !!       
  !!  Arguments: 
  !!             Np: integer(in)
  !!             phiIn(nIn,nphi+1): real(8) (out)
  !!             nvar,nphi,iphiFirst,iphiLast,iMixFrac: integer(out)
  !!             phi(Np,nphi+1): real(8)(out)
  !=====================================================================
  subroutine init_CH4_JUMK_CNF2009(Np,nvar,nphi,nIn,iphiFirst,iphiLast,iMixFrac,hF,hO2)
    implicit none
    
    integer,intent(in)  :: Np                                            
    integer,intent(out) :: nvar 
    integer,intent(out) :: iphiFirst 
    integer,intent(out) :: iphiLast 
    integer,intent(out) :: iMixFrac 
    integer,intent(out) :: nphi,nIn
    real(8),intent(out)  :: hF, hO2

    integer :: i
    integer :: iInF,iInO
    integer ::   iH2, iO2, iCO2, iH2O2, iHO2, iH2O, iCH4, iCH3, &
                 iCH3OOH, iH2CO, iHO, iCO, iH, iCH3O, iHCO, iCH, &
                 iCH2, iC, iO, iCH2OH, iCH3OH, iHOCO, iHOCHO, iOCHO, & 
                 iCH3OO, iCH2OOH, iHOOCHO, iHOOCO, iOOCHO 
                            
    real(8) :: Tfuel, To2
    real(8) :: Y0F, Y0O, Y0N2
    
    !! call init( cpr, "CH4_JUMK_CNF2009.inp" ) 
    !! 2.
    !... init hectars ...
    !! call Init_phase(htr,"gri_29sp.cti","gri29_mix")
    call init_phase(htr,"CH4_JUMK_CNF2009_Modified.cti","gas")
       
    !! nphi = get_nphi( cpr )
    !! 3.
    !! nSpecs=nSpecies(htr)  ! number of species
    nphi=nSpecs+1         ! size of array phi (num. species+enthalpy)
    iMixFrac = nphi+1       ! mixture fraction index
    nvar= nphi+1            ! number of variables in phi array
    iphiFirst = 1
    iphiLast = nphi
    nIn =2     !number of inlet streams=2 (Fuel and Oxidizer)
    
    allocate( phiIn(nIn,nvar))
    allocate( phi(nvar,Np))
    
   !... get indices
    !! 4.
    iH2     = speciesIndex( htr, "H2" )
    iO2     = speciesIndex( htr, "O2" )
    iCO2    = speciesIndex( htr, "CO2" )
    iH2O2   = speciesIndex( htr, "H2O2" )
    iHO2    = speciesIndex( htr, "HO2" )
    iH2O    = speciesIndex( htr, "H2O" )
    iCH4    = speciesIndex( htr, "CH4" )
    iCH3    = speciesIndex( htr, "CH3" )
    iCH3OOH = speciesIndex( htr, "CH3OOH" )
    iH2CO   = speciesIndex( htr, "H2CO" )
    iHO     = speciesIndex( htr, "HO" )
    iCO     = speciesIndex( htr, "CO" )
    iH      = speciesIndex( htr, "H" )
    iCH3O   = speciesIndex( htr, "CH3O" )
    iHCO    = speciesIndex( htr, "HCO" )
    iCH     = speciesIndex( htr, "CH" )
    iCH2    = speciesIndex( htr, "CH2" )
    iC      = speciesIndex( htr, "C" )
    iO      = speciesIndex( htr, "O" )
    iCH2OH  = speciesIndex( htr, "CH2OH" )
    iCH3OH  = speciesIndex( htr, "CH3OH" )
    !!iHOCO   = speciesIndex( htr, "HOCO" ) 
    !!iHOCHO  = speciesIndex( htr, "HOCHO" )
    iOCHO   = speciesIndex( htr, "OCHO" )
    iCH3OO  = speciesIndex( htr, "CH3OO" )
    !!iCH2OOH = speciesIndex( htr, "CH2OOH" )
    iHOOCHO = speciesIndex( htr, "HOOCHO" )
    iHOOCO  = speciesIndex( htr, "HOOCO" )
    iOOCHO  = speciesIndex( htr, "OOCHO" )

    !... mole fractions on the fuel side
    phiIn = 0
    iInF=1
    phiIn(iInF,iCH4) = 1.d0    
    phiIn(iInF,iMixFrac) = 1.d0
    Tfuel = 6.d2 !Temperature [K]
    !! call set_YT( cpr, phiIn(iInF,:), phiIn(iInF,:), Tfuel )
    !! 4.
    call setState_TPY(htr,Tfuel,press(1),phiIn(iInF,1:nphi-1))
    !! hF=get_h( cpr, phiIn(iInF,:) ) !enthalpy [ergs/gm] cgs units 
    !! 5.
    
    hF=enthalpy_mass(htr) ![J/kg]
    phiIn(iInF,nphi)=hF ![J/kg]
    Y0F = phiIn(iInF,iCH4)
    
    !... mole fractions on the oxidizer side
    iInO=2
    phiIn(iInO,iO2 ) = 1.0d0
    phiIn(iInO,iMixFrac) = 0.d0
    To2 = 6.d2 !Temperature [K]
    !! call set_YT( cpr, phiIn(iInO,:), phiIn(iInO,:), To2 )
    !! 6.
    call setState_TPY(htr,To2,press(1),phiIn(iInO,1:nphi-1))
    !! ho2=get_h( cpr, phiIn(iInO,:) ) !enthalpy [ergs/gm] cgs units
    !! 7.
    !!phiIn(iInO,nphi)=Rmix(htr)*temperature(htr) ![J/kg]
    hO2=enthalpy_mass(htr) ![J/kg]
    phiIn(iInO,nphi)=hO2
    Y0O = phiIn(iInO,iO2 )
    
    !... mixture fraction ...
    call random_number(phi(iMixFrac,:)) !uniform random number for mixture fraction
    
    !... init based on fast chemistry ...
    call fast(Np,nvar,iMixFrac,nIn,iInF,iInO,Y0F,Y0O,Y0N2,Tfuel,To2)
    
  ! Check the initialization of Burke Schumann spectrum
    open(unit=1, file='phiInit.txt')
    do i=1,Np
       !! write(1,'(1x,7(e15.9,1x))') phi(7,i),phi(2,i),phi(3,i),phi(6,i),phi(8,i),&
       !! get_T(cpr,phi(:,i)),phi(nvar,i)
       !! 8.      
       write(1,'(1x,7(e15.9,1x))') phi(7,i),phi(2,i),phi(3,i),phi(6,i),phi(8,i),&
       temperature(htr),phi(nvar,i)
    enddo
       
  end subroutine init_CH4_JUMK_CNF2009
  
!========================================================================================
!!   Sub: fast
!========================================================================================

subroutine fast (Np,nvar,iMixFrac,nIn,iInF,iInO,Y0F,Y0O,Y0N2,Tfuel,To2)

   implicit none

   integer, intent(in) :: Np, nvar, iMixFrac, nIn, iInF, iInO
   real(8), intent(in) :: Y0F, Y0O, Y0N2
   real(8), intent(in) :: Tfuel,To2
   real(8) :: MixFracSt, C_co2, hF,hO, MWs(nvar), T0,T1, C_H2O
   real(8) :: massfrac(nvar), p, TInit
   integer :: n, nphi, err
   integer :: iCH4, iO2, iCO2, iH2O, iN2
   real(8) :: h(Np)

   nphi = nvar-1
   TInit = 1.2d3

 ! get indices
   !! 1.
   iCH4  = speciesIndex( htr, "CH4" )
   iO2   = speciesIndex( htr, "O2" )
   iCO2  = speciesIndex( htr, "CO2" )
   iH2O  = speciesIndex( htr, "H2O" )
   iN2   = speciesIndex( htr, "N2" )
   
 ! Molecular weight
   !!call get_MWs(cpr, MW)
   !! 2.
   call getMolecularWeights(htr, MWs)
   
 ! Stocheometric mixture fraction
   !mixFracSt = 1./(1.+2.*MW(iO2)/MW(iCH4))
   mixFracSt = 1./(1.+2.*MWs(iO2)/MWs(iCH4)*Y0F/Y0O)
   C_co2=mixFracSt*MWs(iCO2)/MWs(iCH4)
   C_H2O=mixFracSt*2*MWs(iH2O)/MWs(iCH4)
   
 ! Enthalpy at F and O sides
   !! hF=get_h( cpr, phiIn(iInF,:) ) !enthalpy [ergs/gm] cgs units 
   !! 3.
   call setState_TPY(htr,Tfuel,press(1),phiIn(iInF,1:nphi-1))
   phiIn(iInF,nphi)=enthalpy_mass(htr)   
   hF=phiIn(iInF,nphi)
   call setState_TPY(htr,To2,press(1),phiIn(iInO,1:nphi-1))
   phiIn(iInO,nphi)=enthalpy_mass(htr)   
   hO=phiIn(iInO,nphi)
   
   do n=1,Np
     massfrac=0.
     
   ! total enthalpy
     h(n)=(hF-hO)*phi(iMixFrac,n)+hO
     
   ! Fuel
     if (phi(iMixFrac,n) <= mixFracSt) then
        massfrac(iCH4)=0.
     else
        massfrac(iCH4)=Y0F/(1-mixFracSt)*(phi(iMixFrac,n)-mixFracSt)
     endif
     
   ! O2
     if (phi(iMixFrac,n) <= mixFracSt) then
        massfrac(iO2)=Y0O/mixFracSt*(mixFracSt-phi(iMixFrac,n))
     else
        massfrac(iO2)=0.
     endif
     
   ! CO2
     if (phi(iMixFrac,n) <= mixFracSt) then
        massfrac(iCO2)=C_co2/mixFracSt*phi(iMixFrac,n)
     else
        massfrac(iCO2)=C_co2/(1-mixFracSt)*(mixFracSt-phi(iMixFrac,n))+C_co2
     endif
     
   ! H2O
     if (phi(iMixFrac,n) <= mixFracSt) then
        massfrac(iH2O)=C_H2O/mixFracSt*phi(iMixFrac,n)
     else
        massfrac(iH2O)=C_H2O/(1-mixFracSt)*(mixFracSt-phi(iMixFrac,n))+C_H2O
     endif
     
   ! N2
     massfrac(iN2)=1.-massfrac(iO2)-massfrac(iCO2)-massfrac(iCH4)-massfrac(iH2O)
     !!write (*,'(/A, d25.17)') "N2 mass fraction in fast is ", massfrac(iN2)
     
   ! Mass Fractions
     !! call set_YT( cpr, phi(:,n), massfrac, TInit)
   ! Temperature
     !! call set_h( cpr, phi(:,n), h(n))
     !! 4.
     phi(1:nphi-1,n)=massfrac(:)
     call setMassFractions(htr,phi(1:nphi-1,n))
     call setState_HP(htr,h(n),press(1))     
     phi(nphi,n)=enthalpy_mass(htr)
     
  enddo

end subroutine fast

end program PaSR

!----------------------------------------------------------------------------------------
subroutine RCCEInit(printStatus,iTime,Np,nCnstrnts,nSpecs,nphi,htr,phi,Temp,press,aMtrx,&
           massFrcPrtl,massFrc2Prtl,gammaPrtl)
   
   !! use mod_cpreactor
   !! 1.
   use mod_hectars   
   use RCCE_Mod
   implicit none
   
   logical, intent(in) :: printStatus   
   integer, intent(in) :: iTime,Np,nCnstrnts,nSpecs,nphi
   !! 2.
   type(hectars), intent(inout) :: htr
   real(8), intent(inout) :: phi(nphi+1,Np)
   real(8), intent(inout) :: Temp(Np)
   real(8), intent(in) :: press(Np)
   real(8), intent(out) :: aMtrx(nCnstrnts,nSpecs)
   real(8), intent(out) :: massFrcPrtl(Np,nSpecs)
   real(8), intent(out) :: massFrc2Prtl(Np,nSpecs)    
   real(8), intent(out) :: gammaPrtl(Np,nCnstrnts)  

   integer :: errorPrtcle(Np,1)
   character(8) :: errorSpc(Np,10)
   real(8) :: gamma(nCnstrnts)
   real(8) :: g(nSpecs)  
   integer :: i,j,k 
   character(80) :: fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,fmt7,fmt8,fmt9,fmt10,fmt11

   
   call formats(nCnstrnts,nSpecs,Np,fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,fmt7,fmt8,fmt9,fmt10,fmt11)
         
 ! Calculation of RCCE initial condition for each particle  
   do k=1,Np
   
      if (printStatus) then
      print '(/ "------------------------------------------------------------------------------------------------------------")'
      write (*,'(/A, i5)') "Printing matrices in RCCE Initial Condition formulation for Particle number: ", k
      endif
       
    ! Assigning a very small value to species which have zero mass 
    ! fractions which do not exist in Burke Schumann initial condition
      do j=1,nSpecs
         if (phi(j,k) < 1.0d-16) then
            phi(j,k) = 1.0d-16
         endif
         massFrcPrtl(k,j) = phi(j,k)
      enddo
       
    ! Printing mass fractions from Burke Schumann initial condition
      if (printStatus) then
      write (*,'(/A)') "Initialization of phi "
      write(*, fmt4) "The massFrc matrix is :", "Matrix", &
      (phi(j,k),getSpeciesName(htr,j),j = 1,nSpecs)
      !! (phi(j,k),get_species_name(cpr,j),j = 1, nSpecs) 
      !! 4.
      endif
    
    ! Setting state
      !!call setMassFractions(htr,phi(1:nphi-1,k))
      !!call setState_HP(htr,phi(nphi,k),press)
      call setState_TPY(htr,Temp(k),press(k),phi(1:nphi-1,k))
            
    ! Printing temperature to be checked! 
      !! 5.
      if (printStatus) then
      write(*,'(/A, d25.8)') "Runiv(htr)", Runiv(htr)
      write (*,'(/A, e17.10)') "Temperature in RCCEInit before call to g is ", &
      temperature(htr)
      write (*,'(/A, e17.10)') "Runiv in RCCEInit before call to g is ", &
      Runiv(htr)  ![J/kmole k]    
      endif
       
    ! Calculationg standard state Gibbs free energies of species from Chemckin. 
    ! Stndard state Gibbs free energies are calculated in cgs units. 
      !! call get_g(cpr,phi(:,k),g(:))  ![ergs/mole]
      !! 6.
      call getPureSpeciesGibbsPerMole(htr,g) ![J/kmole]
            
      if (printStatus) then
      write(*, fmt4) "Gibbs Free Energy vector:", "Matrix", (g(j), &
      !! get_species_name(cpr, j), j = 1, nSpecs)
      !! 7.
      getSpeciesName(htr,j), j = 1, nSpecs)
      endif
      
    ! Nondimensionalizing of Stndard state Gibbs free energies of species 
    ! Universal gas constant (R_UNIVERSAL in cpreactor.f file) = 83144720._PREC [ergs/mole.K]   
      do j=1,nSpecs
         !! g(j) = g(j)/(get_Runiv(cpr)*get_T(cpr,phi(:,k)))
         !! 8.
         g(j) = g(j)/(Runiv(htr)*temperature(htr))      
      enddo
            
      if (printStatus) then
      write(*, fmt4) "Gibbs Free Energy vector:", "Matrix", (g(j), &
      !! get_species_name(cpr, j), j = 1, nSpecs)
      !! 7.
      getSpeciesName(htr,j), j = 1, nSpecs)
      endif
             
    ! Call RCCE initial condition subroutine      
      call RCCEInit_Sub(printStatus,iTime, Np, nSpecs, nphi, htr, &
           phi(:,k), k, nCnstrnts, massFrcPrtl(k,:) , gamma, &
           g, aMtrx, massFrc2Prtl(k,:), errorPrtcle, errorSpc, press(k)) !mtrxInput=aT

      call setState_TPY(htr,Temp(k),press(k),phi(1:nphi-1,k))
      phi(nphi,k) = enthalpy_mass(htr)           
     
    ! Saving gamma matrix of each particle into a new variable  
      do i=1,nCnstrnts
         gammaPrtl(k,i) = gamma(i)
      enddo
         
    ! Printion initial RCCE state mass fraction and gamma matrix 
      !! 10. 
      if (printStatus) then
      write(*, fmt4) "RCCE State:", "Matrix", (massFrc2Prtl(k,j), &
      getSpeciesName(htr,j), j = 1, nSpecs)
      write(*, fmt1) "The gamma matrix is:", "Matrix", &
      (gammaPrtl(k,i), i = 1, nCnstrnts)
      write(*, fmt4) "Normalized Free Energy vector vector:", "Matrix", &
      (g(j), getSpeciesName(htr,j), j = 1, nSpecs)
      endif
      
   enddo

end subroutine RCCEInit
   
!---------------------------------------------------------------------------------------------------------------------------
subroutine Calc_rhoTemp(prtlNo,Np,nphi,htr,phi,press,rho,Temp)
   
   !! 1. 
   use mod_hectars
   implicit none
    
   integer, intent(in) :: prtlNo,Np,nphi
   type(hectars), intent(inout) :: htr
   real(8), intent(inout) :: phi(nphi+1)
   real(8), intent(in) :: press
   real(8), intent(out) :: rho  
   real(8), intent(out) :: Temp
   
   integer :: j
   real(8) :: MWs(nphi-1)
   real(8) :: MW

   
 ! Saving temperature of each particle and universal gas constant into new variables 
   !! 2.
   call getMolecularWeights(htr,MWs) ! SI units, kg/mole kg/kmole?
   MW=0
   do j=1,(nphi-1)
      MW=phi(j)/MWs(j)+MW
   end do   
   MW=1.0/MW      
   call setState_TPY(htr,Temp,press,phi(1:nphi-1))
   rho = press / (Runiv(htr)*Temp/MW)
   !!write(*,'(/A, d17.10, A, i10)') "Temperature in Calc_rhoTemp is: ", Temp, "  prtlNo = ", prtlNo
   !!write(*,'(/A, d17.10)') "Density in Calc_rhoTemp is: ", rho
   
   !! Chemkin   
 ! Calculation of the density of each particle  
 ! get_RT returns multiplication of temperature and gas constant of each particle  
 ! R: [ergs/g.K], rho: [g/cm3]

end subroutine Calc_rhoTemp 
!---------------------------------------------------------------------------------------------------------------------------
subroutine Calc_a11aMtrx(Np,nCnstrnts,nSpecs,aMtrx,a11aMtrx)
   
   use Matrix_Mod
   implicit none
   
   integer, intent(in) :: Np,nCnstrnts,nSpecs
   real(8), intent(in) :: aMtrx(nCnstrnts,nSpecs)
   real(8), intent(out) :: a11aMtrx(nCnstrnts,nSpecs)
   
   integer :: i,ii,j
   integer :: rowsC, colsC, ErrCode
   real(8), dimension(nCnstrnts,nCnstrnts) :: mtrxInv, mtrxIn
   real(8) :: a11aMtrxT(nSpecs,nCnstrnts)
   integer :: INDX(nCnstrnts)
   character(80) :: fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,fmt7,fmt8,fmt9,fmt10,fmt11
   
   call formats(nCnstrnts,nSpecs,Np,fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,fmt7,fmt8,fmt9,fmt10,fmt11)
      
   do i=1,nCnstrnts
      do ii=1,nCnstrnts
	     mtrxIn(i,ii) = aMtrx(i,ii)
	  end do
   end do
   
   call MIGS (mtrxIn,nCnstrnts,mtrxInv,INDX)
	
   call MATRIXPRODUCT(mtrxInv, nCnstrnts, nCnstrnts, aMtrx, nCnstrnts, nSpecs, a11aMtrx, rowsC, colsC, ErrCode)
	
 ! For the sake of print only
   do j=1,nSpecs
	  do i=1,nCnstrnts
	     a11aMtrxT(j,i)= a11aMtrx(i,j)
	  end do
   end do

   !!write(*, fmt5) "a11aMtrxT is:", "Matrix", ((a11aMtrxT(i, j), j = 1, size(a11aMtrxT, 2)), i = 1, size(a11aMtrxT, 1))

end subroutine Calc_a11aMtrx
    
!---------------------------------------------------------------------------------------------------------------------------
subroutine RHS(iTime, Np, nCnstrnts, nSpecs, aMtrx, a11aMtrx, &
           htr, phi, nphi, massFrcPrtl, massFrc2Prtl, Temp, rho, &
           press, Cnstrnts, CnstrntsAvg, hAvg, xiAvg, RHSVec,TauMix)
   
   use mod_hectars
   implicit none

   integer, intent(in) :: iTime, Np, nCnstrnts, nSpecs
   real(8), intent(in) :: aMtrx(nCnstrnts,nSpecs)
   real(8), intent(out) :: a11aMtrx(nCnstrnts, nSpecs)
   type(hectars), intent(inout) :: htr
   real(8), intent(in) :: phi(nphi+1,Np)
   integer, intent(in) :: nphi
   real(8), intent(in) :: massFrcPrtl(Np,nSpecs)     
   real(8), intent(in) :: massFrc2Prtl(Np,nSpecs) 
   real(8), intent(in) :: Temp(Np), rho(Np), press(Np)
   real(8), intent(out) :: Cnstrnts(Np,nCnstrnts)
   real(8), intent(out) :: CnstrntsAvg(nCnstrnts) 
   real(8), intent(out) :: hAvg, xiAvg 
   real(8), intent(out) :: RHSVec(Np,nCnstrnts+3) 
   real(8), intent(in) :: TauMix 

   real(8) :: Ch
   real(8) :: MWs(nSpecs)     
   real(8) :: h(Np)
   real(8) :: wdotPrtl(nSpecs,Np)   
   integer :: i,j,k
   character(80) :: fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,fmt7,fmt8,fmt9,fmt10,fmt11

      
   call formats(nCnstrnts,nSpecs,Np,fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,fmt7,fmt8,fmt9,fmt10,fmt11)

   Ch = 1.d0/TauMix
      
   call Calc_a11aMtrx(Np,nCnstrnts,nSpecs,aMtrx,a11aMtrx)
   
 ! Set mass fractions 
 ! gets molecular weight of species [kg/kmol]
   !!call get_MWs(cpr, MWs)
   !! 4.
   call getMolecularWeights(htr, MWs)
   
   do k=1,Np
      do i=1, nCnstrnts
         Cnstrnts(k,i) = 0
         do j=1, nSpecs
            Cnstrnts(k,i) = a11aMtrx(i,j)*phi(j,k)/MWs(j) + Cnstrnts(k,i)
         enddo
      enddo
   enddo
  
   do i=1,nCnstrnts
      CnstrntsAvg(i) = 0
      do k=1,Np
         CnstrntsAvg(i) = Cnstrnts(k,i) + CnstrntsAvg(i)
      enddo
      CnstrntsAvg(i) = CnstrntsAvg(i)/Np
   enddo
       
   hAvg = 0
   do k=1,Np
      call setState_TPY(htr,Temp(k),press(k),phi(1:nphi-1,k))
      h(k) = enthalpy_mass(htr)
      hAvg = h(k) + hAvg
   enddo
   hAvg = hAvg/Np
         
   do k=1,Np
      do i=1,nCnstrnts
         RHSVec(k,i) = 0.d0 !- Ch*(Cnstrnts(k,i)-CnstrntsAvg(i))
      enddo
      RHSVec(k,nCnstrnts+1) = 0.d0 !- Ch*( h(k) - hAvg )
      RHSVec(k,nCnstrnts+2) = 0.d0
   enddo
   
   !!write(*,'(/A,i10)') "iTime in RHS Subroutine:  ", iTime

   call RHSChem(Np, nCnstrnts, nSpecs, htr, phi, nphi, Temp, press, wdotPrtl)

   do k=1,Np
      do i=1,nCnstrnts
         do j=1,nSpecs
            RHSVec(k,i) = RHSVec(k,i) + a11aMtrx(i,j)*wdotPrtl(j,k)/rho(k)
         enddo
      enddo
      !!write(*, '(/A,d17.10)') "rho(k) = ", rho(k)
      !!write(*, fmt2) "The RHSVec matrix is:", "Matrix", (RHSVec(k,i), i = 1, nCnstrnts+2)
   enddo
   
   do k=1,Np
      RHSVec(k,nCnstrnts+1) = RHSVec(k,nCnstrnts+1)
   enddo

   xiAvg = 0
   do k=1,Np
      xiAvg = phi(nphi+1,k) + xiAvg
   enddo
   xiAvg = xiAvg/Np
      
   do k=1,Np
      RHSVec(k,nCnstrnts+3) = 0.d0!- Ch*( phi(nphi+1,k) - xiAvg )
   enddo
     
end subroutine RHS
    
!---------------------------------------------------------------------------------------------------------------------------
subroutine RHSChem(Np, nCnstrnts, nSpecs, htr, phi, nphi, Temp, press, wdotPrtl)

   use mod_hectars
   implicit none
      
   integer, intent(in) :: Np, nCnstrnts, nSpecs
   type(hectars), intent(inout) :: htr
   real(8), intent(in) :: phi(nphi+1,Np)
   integer, intent(in) :: nphi
   real(8), intent(in) :: Temp(Np), press(Np)
   real(8), intent(out) :: wdotPrtl(nSpecs,Np)
   
   integer :: k,j
   real(8) :: wdot(nphi)
   character(80) :: fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,fmt7,fmt8,fmt9,fmt10,fmt11

      
   call formats(nCnstrnts,nSpecs,Np,fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,fmt7,fmt8,fmt9,fmt10,fmt11)
   
 ! wdot chemical production rates of the species 
 ! cgs units, moles/(cm**3*sec)
   
   do k=1,Np
      !! write(*, fmt4) "RCCE State in RHSChem:", "Matrix", (phi(j,1), get_species_name(cpr, j), j = 1, nSpecs)   
      !! call get_wdot( cpr, phi(:,k), wdot )
      !! 1.
      call setState_TPY(htr,Temp(k),press(k),phi(1:nphi-1,k))
      call get_wdot_hectars( htr, phi(:,k), press(k), wdot ) ![moles/m3*sec]
      do j=1,nSpecs
         wdotPrtl(j,k) = wdot(j)
      enddo
      !! write(*, '(/A,i4)') "PrtlNo = ", k
      !! write(*, fmt4) "wdot:", "Matrix", (wdot(j), get_species_name(cpr, j), j = 1, nSpecs)
   enddo
   
end subroutine RHSChem

!---------------------------------------------------------------------------------------------------------------------------
subroutine thermo_variables(PrtlNo,Np,nSpecs,nphi,htr,phi,Temp,press,Runi,MWs,MW,cp,hsp,u,g,P0)

   use mod_hectars
   implicit none
    
   integer, intent(in) :: PrtlNo,Np,nSpecs,nphi
   type(hectars), intent(inout) :: htr
   real(8), intent(in) :: phi(nphi+1)
   real(8), intent(in) :: Temp, press
   real(8), intent(out) :: Runi
   real(8), intent(out) :: MWs(nSpecs)
   real(8), intent(out) :: MW
   real(8), intent(out) :: cp
   real(8), intent(out) :: hsp(nSpecs)
   real(8), intent(out) :: u(nSpecs)    
   real(8), intent(out) :: g(nSpecs) 
   real(8), intent(out) :: P0
   
   integer :: j
   
   call setState_TPY(htr,Temp,press,phi(1:nphi-1))
   
 ! Universal gas constant (R_UNIVERSAL in cpreactor.f file) = 83144720._PREC [ergs/mole.K]  
   !! 1.
   Runi = Runiv(htr) ! 8.3144720_PREC [J/mole.K]
    
 ! Calculationg molecular weights from Chemckin.
 ! mean molecular weights are calculated in cgs units. 
 ! MWs: [gm/mole]         
   !! call get_MWs(cpr,MWs)          ! cgs units, gm/mole  
   !! 2.
   call getMolecularWeights(htr,MWs) ! SI units, kg/mole kg/kmole?
   MW=0
   !! do i=1,nCnstrnts
   do j=1,nSpecs
      MW=phi(j)/MWs(j)+MW
   end do   
   MW=1.0/MW
        
 ! Calculationg mean specific heat at constant pressure from Chemckin.
 ! mean specific heat at constant pressure is calculated in cgs units. 
 ! cp: [ergs/gm*K]    
   !! cp = get_cp(cpr,phi(:))     ! cgs units, ergs/gm*K
   !! 3.
   cp = cp_mass(htr) ![J/kg*K]
   
 ! Calculationg enthalpies of species from Chemckin. 
 ! enthalpies of species are calculated in cgs units. 
 ! hsp: [ergs/gm]        
   !! call get_hsp(cpr,phi(:),hsp(:)) ! cgs units, ergs/gm
   !! 4.
   call getEnthalpies(htr,hsp) ![J/kg]
     
 ! Calculationg internal energies of species from Chemckin. 
 ! internal energies of species are calculated in cgs units. 
 ! u: [ergs/mole]        
   !! call get_u(cpr,phi(:),u(:)) ! cgs units, ergs/mole
   !! 5.
   call getIntEnergyPerMole(htr,u) ![J/kmole]
 
 ! Calculationg standard state Gibbs free energies of species from Chemckin. 
 ! Stndard state Gibbs free energies are calculated in cgs units. 
 ! g: [ergs/mole] 
   !call get_g(cpr,phi(:),g(:))
  !! 6.
  call getPureSpeciesGibbsPerMole(htr,g) ![J/kmole] 
     
! cgs units, dyn·cm−2 or (P0 = 101325 Pa, 1dyn·cm−2 = 0.1 Pa) 
  !! P0 = 1.013250e6  
  !! 7.
  P0 = 1.013250d5 ![N/m2]
         
end subroutine thermo_variables 
           
!---------------------------------------------------------------------------------------------------------------------------
subroutine read_phi(restart,restartITime,restartTime,Np,nSpecs,nphi,phi)
   
   implicit none
   
   logical, intent(in) :: restart 
   integer, intent(in) :: restartITime,Np,nSpecs,nphi
   real(8), intent(in) :: restartTime
   real(8), intent(out) :: phi(nphi+1,Np)

   integer :: k,j   
   character(11) :: str1
   character(4) :: str2
   character(12) :: str3
   character(27) :: str4
   
   if ( restart == .true. ) then
      write(str3,*) restartITime
      str1 = 'Restart_phi'
      str2 = '.txt'
      str4 = str3//str1//str2
      str4 = adjustl(str4)
      str4 = trim(str4)
   else
      str4 ='Init_phi.txt'
    endif
   
   open(unit=5, file = str4)
   do k=1,Np
      write(5,'(1x,d25.18,1x)') restartTime
      do j=1,nSpecs
         read(5,'(1x,d25.18,1x)') phi(j,k)
      enddo
      read(5,'(1x,d25.18,1x)') phi(nphi,k)
      read(5,'(1x,d25.18,1x)') phi(nphi+1,k)
   enddo 
   close(5)
    
end subroutine read_phi
       
!---------------------------------------------------------------------------------------------------------------------------
subroutine save_phi(restart,restartITime,iTime,t,deltaT,Np,nphi,htr,phi,nvar,Temp,nDt,average)

   use mod_hectars
   implicit none
   logical, intent(in) :: restart
   integer, intent(in) :: restartITime,iTime,Np,nphi
   real(8), intent(in) :: t,deltaT
   type(hectars), intent(inout) :: htr
   real(8), intent(in) :: phi(nphi+1,Np)
   integer, intent(in) :: nvar,nDt
   real(8), intent(in) :: Temp(Np)
   real(8), intent(inout) :: average(nvar,nDt)  
    
   integer :: k,j 
   integer ::   iCH4, iO2, iCO2, iH2O
   real(8) :: TempAvg
   
   character(3) :: str1
   character(4) :: str2
   character(12) :: str3
   character(19) :: str4
   
   if (restart == .true. ) then
       write(str3,*) restartITime
      str1 = 'phi'
      str2 = '.txt'
      str4 = str3//str1//str2
      str4 = adjustl(str4)
      str4 = trim(str4)
   else
      str4 = 'phi.txt'
   endif
    
   open(unit=2, file = str4)
   open(unit=3, file = 'TempAvg.txt')
   
   iCH4 = speciesIndex(htr,"CH4")
   iO2  = speciesIndex(htr,"O2")
   iCO2 = speciesIndex(htr,"CO2")
   iH2O = speciesIndex(htr,"H2O")

 ! New Methane
   do k=1,Np
      write(2,'(1x,7(d25.18,1x))') phi(iCH4,k),&
      phi(iO2,k),phi(iCO2,k),phi(iH2O,k),Temp(k),phi(nphi+1,k),iTime*deltaT
   enddo
   
   do j=1,nvar
      average(j,iTime)=0
   enddo
      
   do k=1,Np
      average(iCH4,iTime) = average(iCH4,iTime) + phi(iCH4,k)
      average(iO2,iTime) = average(iO2,iTime) + phi(iO2,k)
      average(iCO2,iTime) = average(iCO2,iTime) + phi(iCO2,k)
      average(iH2O,iTime) = average(iH2O,iTime) + phi(iH2O,k)
      average(nphi,iTime) = average(nphi,iTime) + Temp(k)
   enddo
      
   average(iCH4,iTime) = average(iCH4,iTime)/Np
   average(iO2,iTime) = average(iO2,iTime)/Np
   average(iCO2,iTime) = average(iCO2,iTime)/Np
   average(iH2O,iTime) = average(iH2O,iTime)/Np
   average(nphi,iTime) = average(nphi,iTime)/Np
            
   write(3,'(1x,6(d17.10,2x))') average(iCH4,iTime),average(iO2,iTime),&
   average(iCO2,iTime),average(iH2O,iTime),average(nphi,iTime),iTime*deltaT
          
end subroutine save_phi
        
!---------------------------------------------------------------------------------------------------------------------------
subroutine writeRestart_phi(iTime,t,Np,nSpecs,nphi,phi)
   
   implicit none
      
   integer, intent(in) :: iTime,Np,nSpecs,nphi
   real(8), intent(in) :: t
   real(8), intent(in) :: phi(nphi+1,Np)

   integer :: k,j
   character(11) :: str1
   character(4) :: str2
   character(12) :: str3
   character(27) :: str4
 
   write(str3,*) iTime
     
   str1 = 'Restart_phi'
   str2 = '.txt'
   !str4 = str1//str2
   str4 = str3//str1//str2
   str4 = adjustl(str4)
   str4 = trim(str4)
   
   open(unit=6, file = str4)
   do k=1,Np
      write(6,'(1x,d25.18,1x)') t
      do j=1,nSpecs
         write(6,'(1x,d25.18,1x)') phi(j,k)
      enddo
      write(6,'(1x,d25.18,1x)') phi(nphi,k)
      write(6,'(1x,d25.18,1x)') phi(nphi+1,k)
   enddo 
   close(6)
    
end subroutine writeRestart_phi
       
!---------------------------------------------------------------------------------------------------------------------------       
subroutine formats(nCnstrnts,nSpecs,Np,fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,fmt7,fmt8,fmt9,fmt10,fmt11)
integer, intent(in) :: nCnstrnts, nSpecs, Np
character(80), intent(out) :: fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,fmt7,fmt8,fmt9,fmt10,fmt11
    
    write (fmt1, '("(/ A / A, "" = ["", ", i0, "( ""["",", i0, "(d17.10, 1x), ""]"" / 10x), ""]"" )")') nCnstrnts, 1
    write (fmt2, '("(/ A / A, "" = ["", ", i0, "( ""["",", i0, "(1x,d17.10, 1x), ""]"" / 10x), ""]"" )")') nCnstrnts+2, 1
    write (fmt3, '("(/ A / A, "" = ["", ", i0, "( ""["",", i0, "(d17.10, 1x), ""]"" / 10x), ""]"" )")') nSpecs, 1
    write (fmt4, '("(/ A / A, "" = ["", ", i0, "( ""["",", i0, "(1x,d17.10, 1x), ""]"", 1x, A / 10x), ""]"" )")') nSpecs, 1
    write (fmt5, '("(/ A / A, "" = ["", ", i0, "( ""["",", i0, "(F10.2, 1x), ""]"" / 10x), ""]"" )")') nSpecs, nCnstrnts
    write (fmt6, '("(/ A / A, "" = ["", ", i0, "( ""["",", i0, "(d17.10, 1x), ""]"" / 10x), ""]"" )")') nCnstrnts, nCnstrnts
    write (fmt7, '("(/ A / A, "" = ["", ", i0, "( ""["",", i0, "(1x,d17.10, 1x), ""]"" / 10x), ""]"" )")') nCnstrnts+2, nCnstrnts+2
    write (fmt8, '("(/ A / A, "" = ["", ", i0, "( ""["",", i0, "(d17.10, 1x), ""]"" / 10x), ""]"" )")') Np, 1
    write (fmt9, '("(/ A / A, "" = ["", ", i0, "( ""["",", i0, "(1x,i2, 1x), ""]"", 1x, i5, 2x, 10A / 10x), ""]"" )")') Np, 1
    write (fmt10, '("(/ A / A, "" = ["", ", i0, "( ""["",", i0, "(F10.2, 1x), ""]"" / 10x), ""]"" )")') nSpecs, nCnstrnts
    write (fmt11, '("(/ A / A, "" = ["", ", i0, "( ""["",", i0, "(d17.10, 1x), ""]"" / 10x), ""]"" )")') nSpecs+2, 1
         
end subroutine formats

!---------------------------------------------------------------------------------------------------------------------------
subroutine header_loop(iTime)

   implicit none
   integer, intent(in) :: iTime
       
   print '(/ "------------------------------------------------------------------------------------------------------------")'
   write (*,'(/A, I10)') "iTime = ", iTime
   
end subroutine header_loop
       
!---------------------------------------------------------------------------------------------------------------------------
subroutine header_main

   implicit none
       
   print '(/ "------------------------------------------------------------------------------------------------------------")'
   write (*,'(/A)') "PaSR reactor"
   write (*,'(/A)') "RCCE Formulation"
   write (*,'(/A)') "GRI Mechanism"
   
end subroutine header_main
       
!---------------------------------------------------------------------------------------------------------------------------
subroutine mix (Np,nvar,deltaT,Tmix,nphi,htr,phi)
   use mod_hectars
   implicit none
   integer, intent(in) :: Np, nvar,nphi
   type(hectars), intent(inout) :: htr
   real(8), intent(in) :: Tmix, deltaT
   real(8), intent(inout), dimension(nvar,Np) :: phi

   real(8), dimension(nvar,Np) :: particlesNew
   real(8), dimension(nvar) :: averageMix
   integer :: j, k

   do j=1,nvar
      averageMix(j)=sum(phi(j,:))/Np
   enddo 
      
!   do k=1,Np
!      do j=1,nvar
!       ! (phi(n+1)-phi(n))/deltaT = -1/Tmix( phi(n+1) - phibar(n) )
!         particlesNew(j,k) = (deltaT/Tmix*averageMix(j)+phi(j,k))/(1+deltaT/(Tmix))
!       ! (phi(n+1)-phi(n))/deltaT = -1/Tmix( phi(n) - phibar(n) )
!         !!particlesNew(j,k) = (1-deltaT/Tmix)*phi(j,k)+deltaT/Tmix*averageMix(j)   
!      enddo
!   enddo

   do k=1,Np
      particlesNew(nvar,k) = (deltaT/Tmix*averageMix(nvar)+phi(nvar,k))/(1+deltaT/(Tmix))
   enddo
     
   do k=1,Np
      !!do j=1,nvar
         phi(nvar,k)=particlesNew(nvar,k)
      !!enddo
   enddo

end subroutine mix
!----------------------------------------------------------------------------------------
subroutine events(Np, mOut, nvar, htr, phi, hF, hO2, press, pOut, mInn, iTime, nDt, pOutArray, Temp, nCnstrnts, gammaPrtl)
   use mod_hectars
   implicit none
   integer, intent(in) :: Np, mOut, nvar, iTime, nDt, nCnstrnts
   type(hectars), intent(inout) :: htr
   real(8), intent(inout) :: phi(nvar, Np)
   real(8), intent(in) ::  hF, hO2
   real(8), intent(in) :: press(Np)
   integer, intent(out)  :: pOut(mOut)
   real(8), intent(out) :: mInn(nvar, mOut)
   integer, intent(in)  :: pOutArray(mOut, nDt)
   real(8), intent(inout) :: Temp(Np)
   real(8), intent(inout) :: gammaPrtl(Np,nCnstrnts)
   integer :: i
      
 ! Returns mOut random numbers in particle pairs
   !!call spreadn(1, Np, mOut, pOut)
   do i=1,mOut
      pOut(i) = pOutArray(i,iTime)
   enddo
   call createmIn(Np, mOut, pOut, nvar, htr, phi, hF, hO2, press, mInn, Temp, nCnstrnts, gammaPrtl)

end subroutine events
!----------------------------------------------------------------------------------------
subroutine createmIn(Np, mOut, pOut, nvar, htr, phi, hF, hO2, press, mInn, Temp, nCnstrnts, gammaPrtl)
   use mod_hectars
   implicit none
   integer, intent(in) :: Np, mOut, nvar, nCnstrnts
   type(hectars), intent(inout) :: htr
   real(8), intent(inout) :: phi(nvar, Np)
   integer, intent(in) :: pOut(mOut)
   real(8), intent(in) :: hF, hO2
   real(8), intent(in) :: press(Np)
   real(8), intent(out) :: mInn(nvar, mOut)
   real(8), intent(inout) :: Temp(Np)
   real(8), intent(inout) :: gammaPrtl(Np,nCnstrnts)
   
   integer :: i,j
   real(8) :: Y0F, Y0O
   character(80) :: fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,fmt7,fmt8,fmt9,fmt10,fmt11
   
 ! nCnstrnts = 12, nSpecs = 26
   call formats(12,26,Np,fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,fmt7,fmt8,fmt9,fmt10,fmt11)
     
 ! Creating a candidate pile consisting of 2mInn phi
   do j=1,nvar
      mInn(j,1:mOut)=0.0
   enddo

! RCCE State:
  mInn(1,1:mOut) = 0.1147234754D-06 ! H2                  
  mInn(2,1:mOut) = 0.7991967665D+00 ! O2                  
  mInn(3,1:mOut) = 0.5036263943D-07 ! CO2                 
  mInn(4,1:mOut) = 0.2378620683D-04 ! H2O2                
  mInn(5,1:mOut) = 0.3108129227D-04 ! HO2                 
  mInn(6,1:mOut) = 0.1755176416D-03 ! H2O                 
  mInn(7,1:mOut) = 0.2002517589D+00 ! CH4                 
  mInn(8,1:mOut) = 0.2439233281D-04 ! CH3                 
  mInn(9,1:mOut) = 0.4622658437D-06 ! CH3OOH              
  mInn(10,1:mOut) = 0.2846021510D-03 ! CH2O                
  mInn(11,1:mOut) = 0.4785201590D-07 ! OH                  
  mInn(12,1:mOut) = 0.9383427708D-05 ! CO                  
  mInn(13,1:mOut) = 0.1042345613D-09 ! H                   
  mInn(14,1:mOut) = 0.2496968455D-12 ! CH3O                
  mInn(15,1:mOut) = 0.7304181823D-14 ! HCO                 
  mInn(16,1:mOut) = 0.1000000000D-15 ! CH                  
  mInn(17,1:mOut) = 0.1000000000D-15 ! CH2                 
  mInn(18,1:mOut) = 0.1000000000D-15 ! C                   
  mInn(19,1:mOut) = 0.7089785664D-10 ! O                   
  mInn(20,1:mOut) = 0.3366944785D-10 ! CH2OH               
  mInn(21,1:mOut) = 0.7396517332D-10 ! CH3OH               
  mInn(22,1:mOut) = 0.1000000000D-15 ! OCHO                
  mInn(23,1:mOut) = 0.7497031780D-06 ! CH3OO               
  mInn(24,1:mOut) = 0.1000000000D-15 ! HOOCHO              
  mInn(25,1:mOut) = 0.1000000000D-15 ! HOOCO               
  mInn(26,1:mOut) = 0.1000000000D-15 ! OOCHO
                 
  write(*,'(10i5)') pOut

  do i=1,mOut
     do j=1,nvar
        phi(j,pOut(i))=mInn(j,i)
     enddo
     call setState_TPY(htr,1.5d3,press(pOut(i)),mInn(1:nvar-2,i))
     phi(nvar-1,pOut(i))=enthalpy_mass(htr) !Rmix(htr)*temperature(htr)
     Temp(pOut(i)) = 1.5d3

     gammaPrtl(pOut(i),1) =  0.3200008441D+02  
     gammaPrtl(pOut(i),2) =  0.2818968918D+02  
     gammaPrtl(pOut(i),3) =  0.7905076066D+02  
     gammaPrtl(pOut(i),4) =  0.5507319888D+02  
     gammaPrtl(pOut(i),5) =  0.4113541066D+02  
     gammaPrtl(pOut(i),6) =  0.5392709724D+02  
     gammaPrtl(pOut(i),7) =  0.3453318531D+02  
     gammaPrtl(pOut(i),8) =  0.2623973913D+02  
     gammaPrtl(pOut(i),9) =  0.7525507918D+02  
     gammaPrtl(pOut(i),10) =  0.4788728658D+02  
     gammaPrtl(pOut(i),11) =  0.3827221811D+02  
     gammaPrtl(pOut(i),12) =  0.4725194127D+02  
    
  enddo
   
end subroutine createmIn
!----------------------------------------------------------------------------------------
module numz
   integer, parameter:: b8 = selected_real_kind(14)
   integer gene_size,num_genes
   integer,allocatable :: a_gene(:),many_genes(:,:)
end module
!----------------------------------------------------------------------------------------
function ran1()  !returns random number between 0 - 1 
   use numz
   implicit none
   real ran1,x
   call random_seed()
   call random_number(x) ! built in fortran 90 random number function
   ran1=x
end function ran1
!----------------------------------------------------------------------------------------
function spread(minf,maxf)  !returns random number between min - max
   use numz
   implicit none
   integer, intent(in) :: minf, maxf
   real:: ran1
   integer:: spread
   spread=int((maxf+1-minf)*ran1()+minf)
end function spread
!----------------------------------------------------------------------------------------
subroutine spreadn(minfn,maxfn,n,b) !returns n random numbers between min - max
   use numz
   implicit none
   interface
        function spread(minf,maxf)
           integer, intent(in) :: minf, maxf
        end function
   end interface

   integer, intent(in) :: minfn, maxfn, n
   integer, intent(out), dimension(n) :: b
   
   integer :: i,j
   integer :: check, check2
   
   check=1
   check2=0
   i=0
   !print*,minfn,maxfn
   do while ((i<(n+1)) .AND. (check2==0))
      if (check==1) then
         i=i+1
      end if
      b(i)=spread(minf=minfn,maxf=maxfn)
      !print*,i,b(i)
      check=1
      if (i>1) then
         do j=1,i-1
            if (b(i)==b(j)) then
               check=0
            end if
         enddo
       end if
       if ((check==1) .AND. (i==n)) then
          check2=1
       end if
    enddo 
       
end subroutine spreadn
!----------------------------------------------------------------------------------------