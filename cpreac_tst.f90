module cpreactor_test

   use Matrix_Mod
   use mod_hectars

   implicit none
   
   private

      integer, parameter :: Np=100
      integer, parameter :: nSpecs=26 
      integer, parameter :: nCnstrnts=12 
      integer, parameter :: nphi=nSpecs+1
      integer, parameter :: neq = nCnstrnts+3
      integer :: iTime
      real(8) :: t
      real(8) :: tout
      integer :: PrtlNo
      real(8) :: aMtrx(nCnstrnts,nSpecs)
      type(hectars) :: htr
      real(8) :: phi(nphi+1)
      real(8) :: massFrc2Prtl(nSpecs) 
      real(8) :: CnstrntsAvg(nCnstrnts) 
      real(8) :: hAvg, xiAvg
      real(8) :: gammaPrtl(nCnstrnts)
      real(8) :: rho
      real(8) :: Temp
      real(8) :: press
      real(8) :: RHSVec(neq)
      real(8) :: TauMix
      integer :: NoCalls
    
   public :: dvodeSolver
   public :: dasslSolver
   public :: dlsodaSolver

contains
!---------------------------------------------------------------------------------------------------------------------------
subroutine dlsodaSolver(PrtlNoIn, iTimeIn, NpIn, nCnstrntsIn, nSpecsIn, aMtrxIn, htrIn, &
           phiIn, nphiIn, massFrc2PrtlIn, &
           CnstrntsAvgIn, hAvgIn, xiAvgIn, gammaPrtlIn, rhoIn, TempIn, pressIn, deltaT, RHSVecIn, TauMixIn)
      
  ! Decleration of Variables
    implicit none
    integer, intent(in) :: PrtlNoIn, iTimeIn, NpIn, nCnstrntsIn, nSpecsIn   
    real(8), intent(in) :: aMtrxIn(nCnstrntsIn,nSpecsIn)
    type(hectars), intent(inout) :: htrIn
    real(8), intent(inout) :: phiIn(nphiIn+1)
    integer, intent(in) :: nphiIn
    real(8), intent(inout) :: massFrc2PrtlIn(nSpecsIn) 
    real(8), intent(in) :: CnstrntsAvgIn(nCnstrnts) 
    real(8), intent(in) :: hAvgIn, xiAvgIn    
    real(8), intent(inout) :: gammaPrtlIn(nCnstrntsIn)
    real(8), intent(inout) :: rhoIn
    real(8), intent(inout) :: TempIn
    real(8), intent(in) :: pressIn
    real(8), intent(inout) :: deltaT
    real(8), intent(inout) :: RHSVecIn(neq)
    real(8), intent(in) :: TauMixIn
    
    PrtlNo = PrtlNoIn
    iTime = iTimeIn
    aMtrx = aMtrxIn
    htr = htrIn
    phi = phiIn
    massFrc2Prtl = massFrc2PrtlIn
    CnstrntsAvg = CnstrntsAvgIn 
    hAvg = hAvgIn
    xiAvg = xiAvgIn
    rho = rhoIn
    Temp = TempIn
    press = pressIn
    gammaPrtl = gammaPrtlIn
    RHSVec = RHSVecIn
    TauMix = TauMixIn
    
    call integ_dlsoda(nCnstrnts, gammaPrtl, rho, Temp, deltaT)

    htrIn = htr
    phiIn = phi
    massFrc2PrtlIn = massFrc2Prtl
    rhoIn = rho
    TempIn =Temp 
    gammaPrtlIn = gammaPrtl
   
    
end subroutine dlsodaSolver

!---------------------------------------------------------------------------------------------------------------------------           
! In this part DLSODA is presented.
subroutine integ_dlsoda(nCnstrnts, gammaPrtl, rho, Temp, deltaT)

  ! Decleration of Variables
    implicit none
    integer, intent(in) :: nCnstrnts
    real(8), intent(inout) :: gammaPrtl(nCnstrnts)
    real(8), intent(inout) :: rho
    real(8), intent(inout) :: Temp
    real(8), intent(inout) :: deltaT
 
    real(8) :: y(neq)
    
    integer :: i
    !!integer :: neq 
    integer :: lrw, liw, jac
    integer :: iopt, mf
    integer :: iwork(20+(neq)+10), itask, istate
    real(8) :: rwork(22+(neq)*(neq+9)) ! 22 + NEQ * MAX(16, NEQ + 9)
    real(8) :: atol(neq), rtol
    integer :: itol

    real(8) :: min_dt
    integer :: num_steps, num_fun_eval, num_Jac_eval 

    istate = 1
    itask = 1
    iopt = 0
    lrw = 22+neq*(neq+9)
    liw = 20+neq+10
    mf = 2   
  ! DLSODA: 1/ supply Jacobian matrix as a subroutine, 
  ! 2/the full Jacobian, 5/banded jacobian; for LSODE, MF=22
    
    itol = 2
    iwork(6)=5000
    do i=1, neq
       atol(i) = 1.0d-9
    enddo
    atol(nCnstrnts+1) = 1.0d-9
    rtol = 1.0d-9

    NoCalls = 0
    t = 0.d0
    tout = deltaT
      
    !!call header_devode(t,iTime,PrtlNo,Temp,get_h(cpr, phi),DeltaT)
       
    call GammaRhoTemp2y(Np,nCnstrnts,PrtlNo,y,gammaPrtl,rho,Temp)
    y(nCnstrnts+3) = phi(nphi+1)
            
  ! Run stiff solver DLSODA       
    call DLSODA(fex_cpreactor, neq, y, t, tout, itol, rtol, atol, &
         itask, istate, iopt, rwork, lrw, iwork, liw, jac, mf)
         
    tout = tout + deltaT
            
    call y2GammaRhoTemp(Np,nCnstrnts,PrtlNo,y,gammaPrtl,rho,Temp) 
    phi(nphi+1) = y(nCnstrnts+3)  
    
    call update_variables(iTime, NoCalls, Np, nCnstrnts, nSpecs, PrtlNo, aMtrx, htr, &
         phi, nphi, massFrc2Prtl, CnstrntsAvg, hAvg, xiAvg, &
         gammaPrtl, rho, Temp, press, RHSVec)

  ! Get slover data
    min_dt=rwork(11)
    num_steps=iwork(11)
    num_fun_eval=iwork(12)
    num_Jac_eval=iwork(13)

    open(unit=300, file = 'getSolverData.txt')    
    write ( 300, * ) min_dt,num_steps,num_fun_eval,num_Jac_eval    
        
    open(unit=200, file = 'NumberSubIterations.txt')    
    write ( 200, * ) iTime, PrtlNo, NoCalls
           
end subroutine integ_dlsoda
                 
!---------------------------------------------------------------------------------------------------------------------------
subroutine dasslsolver(PrtlNoIn, iTimeIn, NpIn, nCnstrntsIn, nSpecsIn, aMtrxIn, htrIn, &
           phiIn, nphiIn, massFrc2PrtlIn, &
           CnstrntsAvgIn, hAvgIn, xiAvgIn, gammaPrtlIn, rhoIn, TempIn, pressIn, deltaT, RHSVecIn, TauMixIn)
      
  ! Decleration of Variables
    implicit none
    integer, intent(in) :: PrtlNoIn, iTimeIn, NpIn, nCnstrntsIn, nSpecsIn
    real(8), intent(in) :: aMtrxIn(nCnstrntsIn,nSpecsIn)
    type(hectars), intent(inout) :: htrIn
    real(8), intent(inout) :: phiIn(nphiIn+1)
    integer, intent(in) :: nphiIn
    real(8), intent(inout) :: massFrc2PrtlIn(nSpecsIn) 
    real(8), intent(in) :: CnstrntsAvgIn(nCnstrnts) 
    real(8), intent(in) :: hAvgIn, xiAvgIn    
    real(8), intent(inout) :: gammaPrtlIn(nCnstrntsIn)
    real(8), intent(inout) :: rhoIn
    real(8), intent(inout) :: TempIn
    real(8), intent(in) :: pressIn
    real(8), intent(inout) :: deltaT
    real(8), intent(inout) :: RHSVecIn(neq)
    real(8), intent(in) :: TauMixIn
    
    integer :: i
    character(80) :: fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,fmt7,fmt8,fmt9,fmt10,fmt11
    
    call formats(nCnstrnts,nSpecs,Np,fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,fmt7,fmt8,fmt9,fmt10,fmt11)
    
    PrtlNo = PrtlNoIn    
    iTime = iTimeIn
    aMtrx = aMtrxIn
    htr = htrIn
    phi = phiIn
    massFrc2Prtl = massFrc2PrtlIn
    CnstrntsAvg = CnstrntsAvgIn 
    hAvg = hAvgIn
    xiAvg = xiAvgIn
    rho = rhoIn
    Temp = TempIn
    press = pressIn
    gammaPrtl = gammaPrtlIn
    RHSVec = RHSVecIn
    TauMix = TauMixIn
    
    call integ_dassl(deltaT)
    
    htrIn = htr
    phiIn = phi
    massFrc2PrtlIn = massFrc2Prtl
    rhoIn = rho
    TempIn =Temp 
    gammaPrtlIn = gammaPrtl
   
end subroutine dasslsolver

!---------------------------------------------------------------------------------------------------------------------------
subroutine integ_dassl(deltaT)
    
    implicit none    

  ! Decleration of Variables
    real(8), intent(inout) :: deltaT
    
  ! Internal Variables
    !!integer :: neq
    integer, parameter :: maxord = 5
    integer :: info(15)   
    integer :: idid  
    integer :: lrw, liw  
    integer :: ipar 
    integer :: iwork(20 + (neq) + 4) 
    real(8) :: rpar 
    real(8) :: atol, rtol               
    real(8) :: y(neq)
    real(8) :: ydot(neq)
    real(8) :: rwork(40 + (maxord + 4)*(neq) + (neq)**2 + 10)    
    !!real(8) :: atol(neq), rtol(neq) !If info(2) = 1
    
    integer :: i,j
    real(8) :: t1,t2,t3
    character(80) :: fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,fmt7,fmt8,fmt9,fmt10,fmt11
    integer :: num_steps, num_fun_eval, num_Jac_eval
    
    call formats(nCnstrnts,nSpecs,Np,fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,fmt7,fmt8,fmt9,fmt10,fmt11)

  ! -------------BE CAUCIOUS WITH THIS PART-------------
  ! BEGIN 
    rtol = 1.d-9
    atol = 1.d-9
  ! If info(2) = 1, then you need to asign a vector error control
  ! do i=1,neq
  ! rtol(i) = 1.d-7
  ! atol(i) = 1.d-5
  ! enddo
  ! atol(ncnstrnts+1) = 1.d-7
  ! END       
  ! -------------BE CAUCIOUS WITH THIS PART-------------
    
    do i=1,15
       info(i) = 0
    enddo
    !!info(11) = 1 ! if you want to provide yprime at initial call
    !!info(2) = 1 ! if you want to provide a vector of error control
    !!info(5) = 1 ! if you want to provide jac.
	
	!!if (iTime ==1 ) ! Add in the continuation mode 
	   !!info(1) =0
	!!else
	   !!info(1)=1
	!!endif 
	       
    !!maxord = 5
    liw = 20 + neq + 4 ! 39
    lrw = 40 + (maxord + 4)*neq + neq**2 + 10 ! 410
    
    NoCalls = 0
    t = 0.d0      ! Remove in the continuation mode
    tout = deltaT ! Remove in the continuation mode
            
    !!call header_devode(t,iTime,PrtlNo,Temp,get_h(cpr, phi),DeltaT)
       
     call GammaRhoTemp2y(Np,nCnstrnts,PrtlNo,y,gammaPrtl,rho,Temp)  
     y(nCnstrnts+3) = phi(nphi+1) 
       
     call cpu_time ( t1 )
            
!2    call ydot_external(neq, PrtlNo, Np, nCnstrnts, nSpecs, aMtrx, htr, &
!          phi, nphi, massFrc2Prtl, CnstrntsAvg, hAvg, xiAvg, &
!          rho, Temp, press, gammaPrtl, RHSVec, ydot)
       
     call cpu_time ( t2 )
     !!write ( *, * ) 'Elapsed CPU time after ydot = ', t2 - t1
     
1     call ddassl(res,neq,t,y,ydot,tout,info,rtol,atol,idid,rwork,lrw,iwork,liw,rpar,ipar,jac)

     call cpu_time ( t3 )
     !!write ( *, * ) 'Elapsed CPU time after dasl = ', t3 - t2
        
     if ( (idid == 1) .OR. (idid == 2) .OR. (idid == 3) ) then !.AND. ( PrtlNo == Np)
       
        tout = tout + deltaT

     else if ( (idid == -1) .OR. (idid == -2) .OR. (idid == -3) ) then
       
      ! IDID = -1
      ! The code has taken about 500 steps.
      ! If you want to continue, set INFO(1) = 1 and
      ! call the code again. An additional 500 steps
      ! will be allowed. 
      ! IDID = -2
      ! The error tolerances RTOL, ATOL have been
      ! increased. If you are sure you want to continue
      ! with relaxed error tolerances, set INFO(1)=1 and
      ! call the code again.  
      ! IDID = -3, A solution component is zero and you set the
      ! corresponding component of ATOL to zero. 
      ! I never set atol to zero   
        
        info(1) = 1
        write(*,'(/A,i5)') 'idid = ', idid
        write(*,'(/A,d25.18/)') 'deltaT = ', deltaT
          
        if (idid == -3) then
           write(*, fmt2) "ydot vector is:", "Matrix", (ydot(i), i = 1, (nCnstrnts+3))          
           stop
        endif
        go to 1
       
     else
        
        write(*,'(/A,i5)') 'idid = ', idid
      ! IDID = -6
      ! Repeated error test failures occurred. A singularity in the
      ! solution may be present. you should restart
      ! the integration. (Provide initial values of Y and
      ! YPRIME which are consistent)  
      ! IDID = -7
      ! Repeated convergence test failures occurred. An inaccurate
      ! or ill-conditioned JACOBIAN may be the problem. You
      ! should restart the integration.     
      ! IDID = -8
      ! The matrix of partial derivatives is singular.
      ! Some of your equations may be redundant.
      ! DDASSL cannot solve the problem as stated.
      ! IDID = -9
      ! DDASSL had multiple convergence test
      ! failures, preceeded by multiple error
      ! test failures. You should restart
      ! the integration.  
      ! IDID =-10
      ! IRES=-1. You should restart the integration.
      ! IDID =-11
      ! IRES=-2. Control is being
      ! returned to the calling program. 
      ! IDID= -33
      ! You cannot continue
        
      ! You can restart on idid = -6 or -7 or -9 or -10
      ! You cannot restart on idid = -8 or -11 or -33
                    
!        write(*,'(/A,i5)') 'idid = ', idid
!        go to 2
          
     endif
       
     call y2GammaRhoTemp(Np,nCnstrnts,PrtlNo,y,gammaPrtl,rho,Temp) 
     phi(nphi+1) = y(nCnstrnts+3) 
     
     call update_variables(iTime, NoCalls, Np, nCnstrnts, nSpecs, PrtlNo, aMtrx, htr, &
          phi, nphi, massFrc2Prtl, CnstrntsAvg, hAvg, xiAvg, &
          gammaPrtl, rho, Temp, press, RHSVec)
    
   ! RWORK(3)--Which contains the step size H to be
   ! attempted on the next step.
   
     !!if (mod(iTime,10) == 0) then
        !!write (*,'(/A, I10, A, I10, A, d25.18, A, d17.10/)') "iTime: ", iTime, " prtlNo = ", &
        !!PrtlNo, "  Temp: ", get_T(cpr,phi), "  deltaT: ", deltaT 
     !!endif 

	 num_steps=iwork(11)
	 num_fun_eval=iwork(12)
	 num_Jac_eval=iwork(13)
	 open(unit=300, file = 'getSolverData.txt')    
	 write ( 300, * )num_steps,num_fun_eval,num_Jac_eval,NoCalls
	             
end subroutine integ_dassl

!---------------------------------------------------------------------------------------------------------------------------
subroutine jac
   implicit none
end subroutine jac

!---------------------------------------------------------------------------------------------------------------------------
subroutine res(t,y,ydot,delta,ires,rpar,ipar)
   implicit none
   real(8) :: t, y(neq), ydot(neq), tout, rtol, atol, & !rtol(neq), atol(neq), &
   rwork(15), rpar, delta(neq) 
   integer :: ipar, ires !neq

   call res_external( neq, t, y, ydot, delta, ires, PrtlNo, iTime, Np, nCnstrnts, &
        nSpecs, aMtrx, htr, phi, nphi,  &
        massFrc2Prtl, CnstrntsAvg, hAvg, xiAvg, gammaPrtl, &
        rho, Temp, press, RHSVec, NoCalls, TauMix)
       
end subroutine res

!---------------------------------------------------------------------------------------------------------------------------
subroutine dvodeSolver(PrtlNoIn, iTimeIn, NpIn, nCnstrntsIn, nSpecsIn, aMtrxIn, &
           htrIn, phiIn, nphiIn, massFrc2PrtlIn, &
           CnstrntsAvgIn, hAvgIn, xiAvgIn, gammaPrtlIn, rhoIn, TempIn, pressIn, deltaT, RHSVecIn, TauMixIn)
      
  ! Decleration of Variables
    implicit none
    integer, intent(in) :: PrtlNoIn, iTimeIn, NpIn, nCnstrntsIn, nSpecsIn
    real(8), intent(in) :: aMtrxIn(nCnstrntsIn,nSpecsIn)
    type(hectars), intent(inout) :: htrIn
    real(8), intent(inout) :: phiIn(nphiIn+1)
    integer, intent(in) :: nphiIn
    real(8), intent(inout) :: massFrc2PrtlIn(nSpecsIn) 
    real(8), intent(in) :: CnstrntsAvgIn(nCnstrnts) 
    real(8), intent(in) :: hAvgIn, xiAvgIn    
    real(8), intent(inout) :: gammaPrtlIn(nCnstrntsIn)
    real(8), intent(inout) :: rhoIn
    real(8), intent(inout) :: TempIn
    real(8), intent(in) :: pressIn
    real(8), intent(inout) :: deltaT
    real(8), intent(inout) :: RHSVecIn(neq)
    real(8), intent(in) :: TauMixIn
    
    integer :: k,j
    character(80) :: fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,fmt7,fmt8,fmt9,fmt10,fmt11


    call formats(nCnstrnts,nSpecs,Np,fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,fmt7,fmt8,fmt9,fmt10,fmt11)
        
    PrtlNo = PrtlNoIn
    iTime = iTimeIn
    aMtrx = aMtrxIn
    htr = htrIn
    phi = phiIn
    massFrc2Prtl = massFrc2PrtlIn
    CnstrntsAvg = CnstrntsAvgIn 
    hAvg = hAvgIn
    xiAvg = xiAvgIn
    rho = rhoIn
    Temp = TempIn
    press = pressIn
    gammaPrtl = gammaPrtlIn
    RHSVec = RHSVecIn
    TauMix = TauMixIn
    
    call integ_vode(gammaPrtl, rho, Temp, deltaT)

    htrIn = htr
    phiIn = phi
    massFrc2PrtlIn = massFrc2Prtl
    rhoIn = rho
    TempIn =Temp 
    gammaPrtlIn = gammaPrtl
   
!    print*,"HEY in dvodeSolver"
!    print*,PrtlNo
!    write(*, fmt4) "Mass fractions are :", "Matrix", (phi(j,PrtlNo), &
!                   getSpeciesName(htr, j), j = 1, nSpecs) 
        
end subroutine dvodeSolver

!--------------------------------------------------------------------------------------------------------------------------- 
! In this part DVODE is presented.
subroutine integ_vode(gammaPrtl, rho, Temp, deltaT)

  ! Decleration of Variables
    implicit none
    !!integer, parameter :: neq = nCnstrnts+3
    real(8), intent(inout) :: gammaPrtl(nCnstrnts)
    real(8), intent(inout) :: rho
    real(8), intent(inout) :: Temp
    real(8), intent(inout) :: deltaT
 
    real(8) :: y(neq)
    
    integer :: i 
    integer :: i1,i2,i3,i4,i5,i6,i7, nr, ni
    integer :: max_steps,nrvode,nivode,state
    integer :: ivow(30 + (neq))    
    real(8) :: rvow(32 + 9*(neq) + 2*(neq)*(neq))
    real(8) :: atol(neq), rtol
    
    real(8) :: min_dt
    integer :: num_steps, num_fun_eval, num_Jac_eval 
    real(8) :: t1,t2
   
    max_steps=5000
    do i=1, neq
       atol(i) = 1.0d-9
    enddo
    atol(nCnstrnts+1) = 1.0d-9
    rtol = 1.0d-9
    
    rvow = 0.0
    ivow = 0
    ivow(6) = max_steps
    state = 1
    i1 = 2;  i2 = 1; i3 = 1; i4=0; i5=22; i6=0; i7=0
    nr = size(rvow)
    ni = size(ivow)

    NoCalls = 0
    t = 0.d0
    tout = deltaT         
         
    !!call header_devode(t,iTime,PrtlNo,Temp,get_h(cpr, phi),DeltaT)
       
    call GammaRhoTemp2y(Np,nCnstrnts,PrtlNo,y,gammaPrtl,rho,Temp)
    y(nCnstrnts+3) = phi(nphi+1)

    call cpu_time ( t1 )
    if ((iTime == 1) .OR. (iTime == 50) .OR. (iTime == 100) .OR. (iTime == 150)) then
    if (PrtlNo == 1) then
       write ( *, * ) 'Elapsed CPU time in cpreac_tst before call to DVODE = ', t1
    endif
    endif
     
  ! Run stiff solver DVODE
    call DVODE (fex_cpreactor, neq, y, t, tout, &
         i1, rtol, atol, &
         i2, state, i3, rvow, nr, ivow, &
         ni, i4, i5, i6, i7 )
         
    call cpu_time ( t2 )
    if ((iTime == 1) .OR. (iTime == 50) .OR. (iTime == 100) .OR. (iTime == 150)) then
    if (PrtlNo == 1) then
       write ( *, * ) 'Elapsed CPU time for DVODE = ', t2-t1
    endif
    endif
           
  ! Get slover data
    min_dt=rvow(11)
    num_steps=ivow(11)
    num_fun_eval=ivow(12)
    num_Jac_eval=ivow(13)

    open(unit=300, file = 'getSolverData.txt')    
    write ( 300, * ) min_dt,num_steps,num_fun_eval,num_Jac_eval    
        
    open(unit=200, file = 'NumberSubIterations.txt')    
    write ( 200, * ) iTime, PrtlNo, NoCalls
                
    call y2GammaRhoTemp(Np,nCnstrnts,PrtlNo,y,gammaPrtl,rho,Temp) 
    phi(nphi+1) = y(nCnstrnts+3)
    
    call update_variables(iTime, NoCalls, Np, nCnstrnts, nSpecs, PrtlNo, aMtrx, htr, &
         phi, nphi, massFrc2Prtl, CnstrntsAvg, hAvg, xiAvg, &
         gammaPrtl, rho, Temp, press, RHSVec)
                   
end subroutine integ_vode

!---------------------------------------------------------------------------------------------------------------------------
subroutine fex_cpreactor( neq, t, y, ydot, rw, iw )
  implicit none

  integer :: neq
  real(8) :: t, y(neq), ydot(neq), rw(*), iw(*)
  real(8) :: gamma(neq)

! This routine which is a part of the cpreactor calls an external subroutine to get the wdot
! The header has to be consistent with how DVODE needs it
! phi_means is global in this scope and it passed to rcce
  call fex_external( neq, t, y, ydot, PrtlNo, iTime, Np, nCnstrnts, nSpecs, aMtrx, htr, &
       phi, nphi, massFrc2Prtl, CnstrntsAvg, hAvg, xiAvg, &
       gammaPrtl, rho, Temp, press, RHSVec, NoCalls, TauMix)
   
end subroutine fex_cpreactor

!----------------------------------------------------------------------------------------------------------------------
subroutine update_variables(iTime, NoCalls, Np, nCnstrnts, nSpecs, PrtlNo, aMtrx, &
           htr, phi, nphi, massFrc2Prtl, CnstrntsAvg, hAvg, xiAvg, &
           gammaPrtl, rho, Temp, press, RHSVec)
   
    integer, intent(in) :: iTime, NoCalls, Np, nCnstrnts, nSpecs, PrtlNo
    real(8), intent(in) :: aMtrx(nCnstrnts,nSpecs)
    type(hectars), intent(inout) :: htr
    real(8), intent(inout) :: phi(nphi+1)
    integer, intent(in) :: nphi
    real(8), intent(inout) :: massFrc2Prtl(nSpecs) 
    real(8), intent(in) :: CnstrntsAvg(nCnstrnts) 
    real(8), intent(in) :: hAvg, xiAvg 
    real(8), intent(in) :: gammaPrtl(nCnstrnts)    
    real(8), intent(in) :: rho
    real(8), intent(in) :: Temp
    real(8), intent(in) :: press
    real(8), intent(in) :: RHSVec(neq) 
    
!    real(8) :: tolRT
!    real(8) :: RTOld
!    real(8) :: RTNew
    real(8) :: tolT
    real(8) :: TOld
    real(8) :: TNew
    integer :: maxIter
    integer :: j
    character(80) :: fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,fmt7,fmt8,fmt9,fmt10,fmt11


    call formats(nCnstrnts,nSpecs,Np,fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,fmt7,fmt8,fmt9,fmt10,fmt11)
      
  ! Print variables befor update
!    call print_variables(iTime, NoCalls, Np, nCnstrnts, nSpecs, PrtlNo, aMtrx, htr, &
!         phi, nphi, massFrc2Prtl, CnstrntsAvg, hAvg, xiAvg, &
!         gammaPrtl, rho, Temp, RHSVec)
    
  ! Update Variables
    tolT = 1.0d-5
    maxIter = 0
    
    !! do while ( (tolRT>1.0d-6) .AND. (maxIter<5))
    !! 1.
    do while ( (tolT>1.0d-6) .AND. (maxIter<1))
       
       !! call Calc_RT(Np,PrtlNo,htr,phi,nphi,Temp,RTOld)
       !! 2.
       call setState_TPY(htr,Temp,press,phi(1:nphi-1))
       TOld = enthalpy_mass(htr)

!       write(*, fmt4) "Mass fractions are :", "Matrix", (phi(j), &
!                   getSpeciesName(htr, j), j = 1, nSpecs) 
                     
     ! Update mass fractions
       call update_massFrcs(iTime, NoCalls, Np, nCnstrnts, nSpecs, PrtlNo, aMtrx, htr, &
            phi, nphi, massFrc2Prtl, CnstrntsAvg, hAvg, xiAvg, &
            gammaPrtl, rho, Temp, press, RHSVec)
                   
       !!write(*, fmt4) "The massFrc matrix in update_variables is :", &
       !!"Matrix", (massFrc2Prtl(j), get_species_name(cpr, j), j = 1, nSpecs) 
       !!write(*,'(/A, e17.10 )') "Temp in update_variables is :", Temp
       !!write(*,'(/A,e17.10)') "get_T(cpr,phi) in update_variables is :", &
       !!get_T(cpr,phi)          
       
       !! call Calc_RT(Np,PrtlNo,htr,phi,nphi,Temp,RTNew)
       !! 3.
       call setState_TPY(htr,Temp,press,phi(1:nphi-1))
       TNew = enthalpy_mass(htr)

     ! Update RT 
       !! phi(nphi) = RTNew
       !! 4.     
       phi(nphi) = enthalpy_mass(htr)
    
       !! call tol_RT(RTOld, RTNew, tolRT)
       !! 5.
       call tol_RT(TOld, TNew, tolT)
       maxIter = maxIter + 1 
       
    enddo
    
  ! Print variables after update
!    call print_variables(iTime, NoCalls, Np, nCnstrnts, nSpecs, PrtlNo, aMtrx, htr, &
!         phi, nphi, massFrc2Prtl, CnstrntsAvg, hAvg, xiAvg, &
!         gammaPrtl, rho, Temp, press, RHSVec)
    
end subroutine update_variables

!---------------------------------------------------------------------------------------------------------------------------
subroutine Calc_RT(Np,PrtlNo,htr,phi,nphi,Temp,RT)

   implicit none
   
   integer, intent(in) :: Np,PrtlNo
   type(hectars), intent(inout) :: htr
   real(8), intent(in) :: phi(nphi+1)
   integer, intent(in) :: nphi 
   real(8), intent(in) :: Temp
   real(8), intent(out) :: RT 
   
   real(8) :: MWs(nphi+1)
   real(8) :: Runi
   real(8) :: mlrMassAvg  
   integer :: j 

    
   !! Runi = get_Runi(cpr)
   !! mlrMassAvg = get_MW(cpr,phi)
   !! 1.
   call getMolecularWeights(htr,MWs) 
         
   mlrMassAvg=0
   do j=1,nphi-1
      mlrMassAvg=phi(j)/MWs(j)+mlrMassAvg
   end do
   mlrMassAvg=1.0/mlrMassAvg
   Runi = Runiv(htr)
   RT = Temp*Runi/mlrMassAvg
    
end subroutine Calc_RT

!---------------------------------------------------------------------------------------------------------------------------
subroutine tol_RT(RTOld, RTNew, tolRT)

   implicit none
   
   real(8), intent(in) :: RTOld
   real(8), intent(in) :: RTNew
   real(8), intent(out) :: tolRT
          
   tolRT = abs ( RTNew - RTOld ) / RTOld
   
   !!write(*,'( /A, d17.10)') "tol_RT = ", tolRT
     
end subroutine tol_RT

!---------------------------------------------------------------------------------------------------------------------------
subroutine update_massFrcs(iTime, NoCalls, Np, nCnstrnts, nSpecs, PrtlNo, aMtrx, htr, &
           phi, nphi, massFrc2Prtl, CnstrntsAvg, hAvg, xiAvg, &
           gammaPrtl, rho, Temp, press, RHSVec)
                   
    integer, intent(in) :: iTime, NoCalls, Np, nCnstrnts, nSpecs, PrtlNo
    real(8), intent(in) :: aMtrx(nCnstrnts,nSpecs)
    type(hectars), intent(inout) :: htr
    real(8), intent(inout) :: phi(nphi+1)
    integer, intent(in) :: nphi
    real(8), intent(inout) :: massFrc2Prtl(nSpecs) 
    real(8), intent(in) :: CnstrntsAvg(nCnstrnts) 
    real(8), intent(in) :: hAvg, xiAvg 
    real(8), intent(in) :: gammaPrtl(nCnstrnts)    
    real(8), intent(in) :: rho
    real(8), intent(in) :: Temp
    real(8), intent(in) :: press
    real(8), intent(in) :: RHSVec(neq)
           
    real(8) :: Runi
    real(8) :: MWs(nSpecs)
    real(8) :: MW
    real(8) :: cp
    real(8) :: hsp(nSpecs)
    real(8) :: u(nSpecs)    
    real(8) :: g(nSpecs) 
    real(8) :: P0
    real(8):: a11aMtrxTGlobal(nSpecs,nCnstrnts)
    real(8):: a11aMtrxTGama(nSpecs)
    integer :: j,i
    character(80) :: fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,fmt7,fmt8,fmt9,fmt10,fmt11


    call formats(nCnstrnts,nSpecs,Np,fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,fmt7,fmt8,fmt9,fmt10,fmt11)
            
  ! Update mass fractions
    call thermo_variables(PrtlNo,Np,nSpecs,nphi,htr,phi,Temp,press,Runi,MWs,MW,cp,hsp,u,g,P0)

!    write(*, fmt4) "Mass fractions are :", "Matrix", (phi(j), &
!                   getSpeciesName(htr, j), j = 1, nSpecs) 
                   
    open(unit=6, file="a11aMtrxTGlobal.txt")	
	do j=1,nSpecs
	   do i=1,nCnstrnts
	      read(6,'(1x,e17.10,1x)') a11aMtrxTGlobal(j,i)
	   end do
	end do	
	close(6)
	    
    a11aMtrxTGama = 0
    a11aMtrxTGama = MATMUL(a11aMtrxTGlobal,gammaPrtl)
       
    do j=1,nSpecs

       massFrc2Prtl(j) = MWs(j)/rho*P0/(Runi*Temp)* &
                               exp(-g(j)/(Runi*Temp))*exp(-a11aMtrxTGama(j))
!       massFrc2Prtl(j)=P0/Press*MWs(j)/MW*exp(-a11aMtrxTGama(j)-g(j)/(Runi*Temp))
                              
       if (massFrc2Prtl(j)<1.0d-16) then
          massFrc2Prtl(j) = 1.0d-16
       endif
       if (massFrc2Prtl(j)>1.0) then
          massFrc2Prtl(j) = 1.0d0
       endif       
       phi(j) = massFrc2Prtl(j)
       
    enddo
    
!    call print_speciesVariables(iTime, 101, Np, nCnstrnts, nSpecs, PrtlNo, 1001, &
!         aMtrx, htr, phi, nphi, massFrc2Prtl, &
!         CnstrntsAvg, hAvg, xiAvg, gammaPrtl, rho, Temp, MWs, Runi, g, a11aMtrxTGama)       
       
end subroutine update_massFrcs 

!---------------------------------------------------------------------------------------------------------------------------
subroutine print_variables(iTime, NoCalls, Np, nCnstrnts, nSpecs, PrtlNo, aMtrx, htr, &
           phi, nphi, massFrc2Prtl, CnstrntsAvg, hAvg, xiAvg, &
           gammaPrtl, rho, Temp, press, RHSVec)
       
    integer, intent(in) :: iTime, NoCalls, Np, nCnstrnts, nSpecs, PrtlNo
    real(8), intent(in) :: aMtrx(nCnstrnts,nSpecs)
    type(hectars), intent(inout) :: htr
    real(8), intent(in) :: phi(nphi+1)
    integer, intent(in) :: nphi
    real(8), intent(in) :: massFrc2Prtl(nSpecs) 
    real(8), intent(in) :: CnstrntsAvg(nCnstrnts) 
    real(8), intent(in) :: hAvg, xiAvg 
    real(8), intent(in) :: gammaPrtl(nCnstrnts)    
    real(8), intent(in) :: rho
    real(8), intent(in) :: Temp
    real(8), intent(in) :: press
    real(8), intent(in) :: RHSVec(neq)
    
    integer :: i
    character(80) :: fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,fmt7,fmt8,fmt9,fmt10,fmt11
    
    call formats(nCnstrnts,nSpecs,Np,fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,fmt7,fmt8,fmt9,fmt10,fmt11)
    
    call setState_TPY(htr,Temp,press,phi(1:nphi-1))
    
    write(*,'(/A, e17.10 )') "Temp  = ", Temp
    write(*,'(/A,e17.10)') "temperature(htr)  =", temperature(htr)
    write(*,'(/A,e17.10,A,e17.10)') "phi(nphi,PrtlNo) =", phi(nphi), &
          "  temperature(htr)*Rmix(htr) =", temperature(htr)*Rmix(htr)
    write(*,'(/A, e17.10 )') "rho(PrtlNo,1)  = ", rho
    write(*, fmt1) "The gammaPrtl is:", "Matrix", (gammaPrtl(i), i = 1, nCnstrnts)
    write(*, fmt4) "Mass fractions are :", "Matrix", (massFrc2Prtl(i), &
                   getSpeciesName(htr, i), i = 1, nSpecs)
    write(*,'(/A, e17.10 )') "Runi  = ", Runiv(htr)
    !!write(*,'(/A,e17.10)') "get_MW(cpr,phi) =", get_MW(cpr,phi)
    write(*, fmt2) "The RHSVec matrix is:", "Matrix", (RHSVec(i), i = 1, neq)
    
end subroutine print_variables

!---------------------------------------------------------------------------------------------------------------------------
subroutine y2GammaRhoTemp(Np,nCnstrnts,PrtlNo,yk,gammaPrtl,rho,Temp)

    implicit none

    integer, intent(in) :: Np, nCnstrnts, PrtlNo
    real(8), intent(in) :: yk(neq)
    real(8), intent(out) :: gammaPrtl(nCnstrnts)    
    real(8), intent(out) :: rho
    real(8), intent(out) :: Temp
    
    integer :: i
    
    do i=1,nCnstrnts
       gammaPrtl(i) = yk(i)
    enddo
    rho = yk(nCnstrnts+1)
    Temp = yk(nCnstrnts+2)  
      
end subroutine y2GammaRhoTemp

!---------------------------------------------------------------------------------------------------------------------------
subroutine GammaRhoTemp2y(Np,nCnstrnts,PrtlNo,yk,gammaPrtl,rho,Temp)

    implicit none

    integer, intent(in) :: Np, nCnstrnts, PrtlNo
    real(8), intent(out) :: yk(neq)
    real(8), intent(in) :: gammaPrtl(nCnstrnts)    
    real(8), intent(in) :: rho
    real(8), intent(in) :: Temp
    
    integer :: i
    
    do i=1,nCnstrnts
       yk(i) = gammaPrtl(i)  
    enddo
    yk(nCnstrnts+1) = rho 
    yk(nCnstrnts+2) = Temp   
      
end subroutine GammaRhoTemp2y

!---------------------------------------------------------------------------------------------------------------------------
subroutine header_devode(t,iTime,PrtlNo,Temp,h,DeltaT)

  implicit none
  real(8), intent(in) :: t
  integer, intent(in) :: iTime,PrtlNo
  real(8), intent(in) :: Temp,h,DeltaT

  if (mod(iTime,10) == 0) then
     write(*,'(/A, d17.10)') "t: ", t
     write (*,'(/A, I10, A, d17.10, A, d17.10, A, d17.10)') "iTime: ", iTime, "  Temp: ", &
            Temp, "  Enthalpy: ", h ,"  DeltaT: ", DeltaT!!, "  particle number ", PrtlNo 
  endif

end subroutine header_devode

!---------------------------------------------------------------------------------------------------------------------------  
end module cpreactor_test
