!--------------------------------------------------------------------------------------------------------------------------
! Author: Fatemeh Hadi
! Name:   RCCEInit_Mod
! This module calculates constraint potentials using particle mass fractions to be used in
! RCCE code. 
!--------------------------------------------------------------------------------------------------------------------------
module RCCE_Mod
   use Matrix_Mod
   use mod_hectars
   use LINPACK_Mod
   use LINPACK_Interface
   implicit none
   private
    
   integer, parameter :: nSpeciesGlobal=26 
   integer, parameter :: nCnstrntsGlobal=12
   real(8):: a11aMtrxTGlobal(nSpeciesGlobal,nCnstrntsGlobal)
   real(8) :: time
            
   public :: RCCEInit_Sub
   public :: res
   public :: ydot_rcce
   public :: fex
   public :: case_func
   public :: Euler   
   public :: update_variables
   public :: print_variables
   public :: thermo_variables

contains
!---------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------RCCE Initial Condition Calculation-----------------------------------------
!--------------------------------------------------------------------------------------------------------------------------- 
subroutine RCCEInit_Sub(printStatus,iTime, Np, nSpecs, nphi, htr, phi, prtlNo, nCnstrnts, massFrc , gama, g, aMtrx, &
           massFrc2, errorPrtcle, errorSpc, press)
    logical, intent(in) :: printStatus
    type(hectars), intent(inout) :: htr
    integer, intent(in) :: iTime, Np, nSpecs, nphi 
    real(8), intent(inout) :: phi(nphi+1)
    real(8), intent(in) :: g(nSpecs)
    integer, intent(in) :: prtlNo, nCnstrnts
	real(8), intent(in) :: massFrc(nSpecs)
	real(8), intent(out) ::gama(nCnstrnts)
	real(8), intent(out) :: aMtrx(nCnstrnts,nSpecs)
	real(8), intent(out) :: massFrc2(nSpecs)
    integer, intent(inout):: errorPrtcle(Np,1)
    character(8), intent(inout):: errorSpc(Np,10)
    real(8), intent(in) :: press
    	
    integer :: i,j,k,ErrCode
    integer :: nRows1, nCols1, nRows2, nCols2, nRows3, nCols3, nRows4, nCols4
    integer :: rowsC, colsC
    integer :: INDX(nCnstrnts) 
    real(8) :: mtrxInv(nCnstrnts,nCnstrnts)
    real(8) :: mtrxMlt(nCnstrnts,nCnstrnts)
    real(8) :: mtrxIn(nCnstrnts,nCnstrnts)
	real(8) :: mu11(nCnstrnts)
    real(8) :: mlrMassAvg,P0
    real(8) ::  a11aMtrx(nCnstrnts,nSpecs) 
    real(8) :: a11aMtrxT(nSpecs,nCnstrnts)
    real(8) :: a11aMtrxTGama(nSpecs,1)
    real(8) :: mlFrc2(nSpecs,1)
    real(8) :: Pj(nCnstrnts), LnPj(nCnstrnts)
    real(8) :: MWs(nSpecs)
    real(8) :: massSum
    
    character(80) :: fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,fmt7,fmt8,fmt9,fmt10
    
    
    call formats(nCnstrnts,nSpecs,Np,fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,fmt7,fmt8,fmt9,fmt10)
    
    if (printStatus) then
    write (*,'(/A)') "RCCE Initial Condition "
    write (*,'(/A)') "Reading matrix dimensions ... "
    endif

    open(unit=11, file="mtrxIn2.txt")
    open(unit=6, file="a11aMtrxTGlobal.txt")
    open(unit=20, file="errors.txt")
    open(unit=30, file="massSumRCCEInit.txt")
    open(unit=40, file="massSumRCCEInitError.txt")
        
    do i=1,nCnstrnts
       do j=1,nSpecs
          read (11,*) aMtrx(i,j)
       end do
    end do
    close(11)
    
    do i=1,nCnstrnts
       do j=1,nCnstrnts
          mtrxIn(i,j)=aMtrx(i,j)
       end do
    end do
    
    call getMolecularWeights(htr,MWs)
    mlrMassAvg=0
    do j=1,nSpecs
       mlrMassAvg=phi(j)/MWs(j)+mlrMassAvg
    end do
    mlrMassAvg=1.d0/mlrMassAvg

    P0 = 1.013250d5![N/m2]    

    do i=1,nCnstrnts
       gama(i) = - (g(i)+log(press*phi(i)*mlrMassAvg/(P0*MWs(i))))
    end do
       
	call MIGS (mtrxIn,nCnstrnts,mtrxInv,INDX)
	call MATRIXPRODUCT(mtrxInv, nCnstrnts, nCnstrnts, aMtrx, nCnstrnts, &
	     nSpecs, a11aMtrx, rowsC, colsC, ErrCode)

	do j=1, nSpecs
	   do i=1, nCnstrnts
	      a11aMtrxT(j,i)= a11aMtrx(i,j)
	   end do
	end do

    call MATRIXPRODUCT(a11aMtrxT, nSpecs, nCnstrnts, gama, nCnstrnts, &
         1, a11aMtrxTGama, rowsC, colsC, ErrCode)
    
    do j=1,nSpecs
       phi(j)=P0/press*MWs(j)/mlrMassAvg*exp(-a11aMtrxTGama(j,1)-g(j))
       if (phi(j)<1.0d-16) then
          phi(j) = 1.0d-16
       endif   
    end do

	do j=1, nSpecs
	   do i=1, nCnstrnts
	      a11aMtrxTGlobal(j,i) = a11aMtrx(i,j)
	      write(6,'(1x,e17.10,1x)') a11aMtrxTGlobal(j,i)
	   end do
	end do	
	close(6)
	                                  
    do i=1, nCnstrnts
       errorPrtcle(prtlNo,1)=0
       if (phi(i)>1.0d-16) then
       if (massFrc(i)>1.0d-16) then
          if (abs(massFrc(i)-phi(i))/massFrc(i)*100>0.1) then
              errorPrtcle(prtlNo,1)=1
          endif
       endif
       endif
    enddo

    do k=1,10  
       errorSpc(prtlNo,k) = " "  
    enddo
    k=1
    do i=(nCnstrnts+1),nSpecs
       if (abs(phi(i))>1) then
           errorSpc(prtlNo,k)=getSpeciesName(htr,i)
           k=k+1 
           if (errorPrtcle(prtlNo,1)==1) then
              errorPrtcle(prtlNo,1)=3
           else
              errorPrtcle(prtlNo,1)=2
           endif
       endif
    enddo
    
    if (PrtlNo==Np) then
       write(*,fmt9) "Particles with wrong RCCE state :", "Matrix", ((errorPrtcle(i,j), j=1, &
                      1), i, (errorSpc(i,k), k=1, 10), i=1, Np)
       write(20,'(/A, i10)') "iTime = ", iTime
       write(20, fmt9) "Particles with wrong RCCE state :", "Matrix", ((errorPrtcle(i,j), j=1, &
                      1), i, (errorSpc(i,k), k=1, 10), i=1, Np)               
    endif

    massSum = 0
    do i=1, nSpecs
       massSum = massSum + phi(i)  
    end do
    write(30,'(/i5,d25.18)') prtlNo, massSum
 
    if (abs(massSum-1) > 0.05) then
       write(*,'(/A,i5,A,d25.18)') "PrtlNo = ", prtlNo, "  massSum = ", massSum
       write(40,'(/i5,d25.18)') iTime, PrtlNo, massSum
    endif
             
    if (printStatus) then            
       write(*,fmt1) "The gamma matrix is:", "Matrix", (gama(i), i=1, nCnstrnts) 
       write(*,fmt6) "Input matrix is:", "Matrix", ((mtrxIn(i,j), j=1, nCnstrnts), &
                     i=1, nCnstrnts)
       write(*,fmt6) "Input matrix is:", "Matrix", ((aMtrx(i,j), j=1, nCnstrnts), &
                     i=1, nCnstrnts)              
       write(*,fmt6) "Inverse matrix is:", "Matrix", ((mtrxInv(i, j), j=1, nCnstrnts), &
                     i=1, nCnstrnts)
       write(*,fmt10) "Multiplication of inverse and aMtrx matrix is:", "Matrix",((a11aMtrx(i,j), &
                     j=1, nSpecs), i=1, nCnstrnts)
       write(*,fmt3) "exponential of minus of Multiplication of a11aMtrxT and  gama is:", "Matrix", &
                     (a11aMtrxTGama(j,1), i=1, nSpecs)
       write(*,fmt4) "RCCE State:", "Matrix", (phi(j), getSpeciesName(htr,j), j=1, nSpecs)
       write(*,fmt4) "Molar mass:", "Matrix", (MWs(j), getSpeciesName(htr,j), j=1, nSpecs)
       write(*,'(/A, d25.8)') "Averaged molar mass", mlrMassAvg
    endif

end subroutine RCCEInit_Sub
!---------------------------------------------------------------------------------------------------------------------------
subroutine res( neq, t, y, ydot, delta, ires, PrtlNo, iTime, Np, nCnstrnts, &
                nSpecs, aMtrx, htr, phi, nphi, massFrc2Prtl, &
                CnstrntsAvg, hAvg, xiAvg, gammaPrtl, rho, Temp, press, RHSVec, NoCalls, TauMix)
                           
    use cpreactor_test
    implicit none
  
    integer, intent(in) :: neq
    real(8), intent(inout) :: t
    real(8), intent(in) :: y(neq)
    real(8), intent(inout) :: ydot(neq)
    real(8), intent(out) :: delta(neq)
    integer, intent(out) :: ires
    integer, intent(in) :: PrtlNo
    integer, intent(in) :: iTime, Np, nCnstrnts, nSpecs
    real(8), intent(in) :: aMtrx(nCnstrnts,nSpecs)
    type(hectars), intent(inout) :: htr
    real(8), intent(inout) :: phi(nphi+1)
    integer, intent(in) :: nphi
    real(8), intent(inout) :: massFrc2Prtl(nSpecs) 
    real(8), intent(in) :: CnstrntsAvg(nCnstrnts) 
    real(8), intent(in) :: hAvg, xiAvg 
    real(8), intent(inout) :: gammaPrtl(nCnstrnts)    
    real(8), intent(inout) :: rho
    real(8), intent(inout) :: Temp
    real(8), intent(in) :: press
    real(8), intent(inout) :: RHSVec(neq) 
    integer, intent(inout) :: NoCalls
    real(8), intent(in) :: TauMix
    
    real(8) :: CfcntMtrx(neq,neq)
    integer :: i,ii,j
    character(80) :: fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,fmt7,fmt8,fmt9,fmt10
    
    NoCalls = NoCalls + 1
    
    call formats(nCnstrnts,nSpecs,Np,fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,fmt7,fmt8,fmt9,fmt10)
        
    !!call header(iTime, PrtlNo, NoCalls)

    call y2GammaRhoTemp(Np,nCnstrnts,PrtlNo,y,gammaPrtl,rho,Temp)
    phi(nphi+1) = y(nCnstrnts+3)
    
  ! Update Variables
    if ((NoCalls == 1) .AND. (iTime == 1)) then
       !!write(*,'(/A, i10)') "First iteration and sub iteration, prtlNo = ", prtlNo
    else
       call update_variables(iTime, NoCalls, Np, nCnstrnts, nSpecs, PrtlNo, aMtrx, htr, &
            phi, nphi, massFrc2Prtl, CnstrntsAvg, hAvg, xiAvg, &
            gammaPrtl, rho, Temp, press, RHSVec)
    endif
   
    call RHS(neq, PrtlNo, iTime, Np, nCnstrnts, nSpecs, aMtrx, &
         htr, phi, nphi, massFrc2Prtl, rho, Temp, &
         press, CnstrntsAvg, hAvg, xiAvg, RHSVec,TauMix)
    
    call RCCECoefficient4_Sub(neq, PrtlNo, Np, nCnstrnts, nSpecs, aMtrx, htr, &
         phi, nphi, massFrc2Prtl, CnstrntsAvg, hAvg, xiAvg, &
         gammaPrtl, rho, Temp, press, CfcntMtrx, RHSVec)          
    
    ires = 0        
    do i=1,nCnstrnts+3
       delta(i) = 0
       do ii=1,nCnstrnts+3
          delta(i) = delta(i) + CfcntMtrx(i,ii)*ydot(ii) 
       enddo
       delta(i) = delta(i) - RHSVec(i)
    enddo
    
    !!write(*,'(/A, d25.18)') "t in res is: ", t
    !!write(*, fmt7) "The CfcntMtrx in res is:", "Matrix", ((CfcntMtrx(i, ii), ii = 1, nCnstrnts+3), i = 1, nCnstrnts+3)     
    !!write(*, fmt2) "The ydot matrix in res is:", "Matrix", (ydot(i), i = 1, nCnstrnts+3)
    !!write(*, fmt2) "The delta matrix in res is:", "Matrix", (delta(i), i = 1, nCnstrnts+3)
    
  ! OR
!    delta = MATMUL(CfcntMtrx,ydot) - RHSVec
    
end subroutine res

!---------------------------------------------------------------------------------------------------------------------------
subroutine ydot_rcce(neq, PrtlNo, Np, nCnstrnts, nSpecs, aMtrx, htr, &
           phi, nphi, massFrc2Prtl, CnstrntsAvg, hAvg, xiAvg, &
           rho, Temp, press, gammaPrtl, RHSVec, ydot)

  use cpreactor_test
  implicit none
   
    integer, intent(in) :: neq, PrtlNo
    integer, intent(in) :: Np, nCnstrnts, nSpecs
    real(8), intent(in) :: aMtrx(nCnstrnts,nSpecs)
    type(hectars), intent(inout) :: htr
    real(8), intent(in) :: phi(nphi+1)
    integer, intent(in) :: nphi
    real(8), intent(in) :: massFrc2Prtl(nSpecs) 
    real(8), intent(in) :: CnstrntsAvg(nCnstrnts) 
    real(8), intent(in) :: hAvg, xiAvg 
    real(8), intent(in) :: rho
    real(8), intent(in) :: Temp
    real(8), intent(in) :: press
    real(8), intent(in) :: gammaPrtl(nCnstrnts)    
    real(8), intent(in) :: RHSVec(nCnstrnts+3) 
    real(8), intent(out) :: ydot(nCnstrnts+3)
    
    real(8) :: CfcntMtrx(nCnstrnts+3,nCnstrnts+3)
    real(8) :: Answer(nCnstrnts+3)
    integer :: info      
    integer :: i
    character(80) :: fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,fmt7,fmt8,fmt9,fmt10
    
    call formats(nCnstrnts,nSpecs,Np,fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,fmt7,fmt8,fmt9,fmt10)
    
    !!write(*, fmt2) "The RHSVec matrix is:", "Matrix", (RHSVec(i), i = 1, nCnstrnts+2)
                 
  ! In this part the matrices which are required for RCCE Mixing formulation are calculated.
    call RCCECoefficient4_Sub(neq, PrtlNo, Np, nCnstrnts, nSpecs, aMtrx, htr, &
         phi, nphi, massFrc2Prtl, CnstrntsAvg, hAvg, xiAvg, &
         gammaPrtl, rho, Temp, press, CfcntMtrx, RHSVec)           

    call LINPACK_Sub(PrtlNo, Np, nCnstrnts, CfcntMtrx, RHSVec, Answer, info) 
    
  ! The RHS vector for the ODE
    do i=1,(nCnstrnts+3)
       ydot(i) = Answer(i)
    enddo
    
    !!write(*, fmt2) "The ydot matrix in ydot_rcce is:", "Matrix", (ydot(i), i = 1, nCnstrnts+2)
        
end subroutine ydot_rcce

!---------------------------------------------------------------------------------------------------------------------------
!This routine computes the wdot. It also uses the cpreactor, so cpreactor routines are
!available here, if needed.
subroutine fex(neq, t, y, ydot, PrtlNo, iTime, Np, nCnstrnts, nSpecs, aMtrx, htr, phi, nphi, massFrc2Prtl, &
           CnstrntsAvg, hAvg, xiAvg, gammaPrtl,  rho, Temp, press, RHSVec, NoCalls,TauMix)
           
  use cpreactor_test
  implicit none

    integer, intent(in) :: neq
    real(8), intent(in) :: t
    real(8), intent(in) :: y(neq)
    real(8), intent(inout) :: ydot(neq)
    integer, intent(in) :: PrtlNo
    integer, intent(in) :: iTime, Np, nCnstrnts, nSpecs
    real(8), intent(in) :: aMtrx(nCnstrnts,nSpecs)
    type(hectars), intent(inout) :: htr
    real(8), intent(inout) :: phi(nphi+1)
    integer, intent(in) :: nphi
    real(8), intent(inout) :: massFrc2Prtl(nSpecs) 
    real(8), intent(in) :: CnstrntsAvg(nCnstrnts) 
    real(8), intent(in) :: hAvg, xiAvg 
    real(8), intent(inout) :: gammaPrtl(nCnstrnts)    
    real(8), intent(inout) :: rho
    real(8), intent(inout) :: Temp
    real(8), intent(in) :: press
    real(8), intent(inout) :: RHSVec(neq) 
    integer, intent(inout) :: NoCalls 
    real(8), intent(in) :: TauMix 
          
    real(8) :: CfcntMtrx(neq,neq)
    real(8) :: Answer(neq)
    integer :: info   
    integer :: i,j   
    character(80) :: fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,fmt7,fmt8,fmt9,fmt10
    real(8) :: t1,t2
    logical :: ifPrint
    
    ifPrint = .false.
    if ((iTime == 1) .OR. (iTime == 50) .OR. (iTime == 100) .OR. (iTime == 150)) then
    if (PrtlNo == 1) then
       ifPrint = .true.
    endif
    endif
    
    NoCalls = NoCalls + 1
    if (NoCalls == 1) then
       if (ifPrint) then
          write(*,'(/A, i10, A, d25.17)') "NoCalls = ", NoCalls, "  time = ", t
       endif
    else
       if (ifPrint) then
          write(*,'(/A, i10, A, d25.17)') "NoCalls = ", NoCalls, "  sub iteration = ", t-time
       endif
    end if
    time =t
    
    call formats(nCnstrnts,nSpecs,Np,fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,fmt7,fmt8,fmt9,fmt10)
        
    !!call header(iTime, PrtlNo, NoCalls)
    
    call y2GammaRhoTemp(Np,nCnstrnts,PrtlNo,y,gammaPrtl(:),rho,Temp)
    phi(nphi+1) = y(nCnstrnts+3)
    
    call cpu_time ( t1 )
        
  ! Update Variables
    if ((NoCalls == 1) .AND. (iTime == 1)) then
       !!write(*,'(/A, i10)') "First iteration and sub iteration, prtlNo = ", prtlNo
    else
    
       call update_variables(iTime, NoCalls, Np, nCnstrnts, nSpecs, PrtlNo, aMtrx, htr, &
            phi, nphi, massFrc2Prtl, CnstrntsAvg, hAvg, xiAvg, &
            gammaPrtl, rho, Temp, press, RHSVec)
            
       call cpu_time ( t2 )
       if (ifPrint) then
          write ( *, * ) 'Elapsed CPU time for update_variables = ', t2-t1
       endif
       call cpu_time ( t1 )
    endif

    call RHS(neq, PrtlNo, iTime, Np, nCnstrnts, nSpecs, aMtrx, &
         htr, phi, nphi, massFrc2Prtl, rho, Temp, &
         press, CnstrntsAvg, hAvg, xiAvg, RHSVec,TauMix)
         
    call cpu_time ( t2 )
    if (ifPrint) then
       write ( *, * ) 'Elapsed CPU time for RHS = ', t2-t1
    endif
    call cpu_time ( t1 )
                                         
  ! In this subroutine the matrices which are required for RCCE Mixing formulation are calculated
    call RCCECoefficient4_Sub(neq, PrtlNo, Np, nCnstrnts, nSpecs, aMtrx, htr, &
         phi, nphi, massFrc2Prtl, CnstrntsAvg, hAvg, xiAvg, &
         gammaPrtl, rho, Temp, press, CfcntMtrx, RHSVec)  
                      
    call cpu_time ( t2 )
    if (ifPrint) then
       write ( *, * ) 'Elapsed CPU time for RCCECoefficient= ', t2-t1
    endif
    call cpu_time ( t1 )
    
  ! In this subroutine ydot dector is calculated 
    call LINPACK_Sub(PrtlNo, Np, nCnstrnts, CfcntMtrx, RHSVec, Answer, info) 
    
    call cpu_time ( t2 )
    if (ifPrint) then
       write ( *, * ) 'Elapsed CPU time for LINPACK= ', t2-t1
    endif
       
    do i=1,neq
       ydot(i) = Answer(i)
    enddo

    !write(*, fmt2) "ydot vector in fex is:", "Matrix", (ydot(i), i = 1, (nCnstrnts+2))
  
end subroutine fex

!---------------------------------------------------------------------------------------------------------------------------
subroutine case_func(neq, iTime, NoCalls, Np, nCnstrnts, nSpecs, PrtlNo, aMtrx, &
           htr, phi, nphi, massFrc2Prtl, CnstrntsAvg, hAvg, xiAvg, &
           gammaPrtl, rho, Temp, press, deltaT, RHSVec, yk, ki)

    implicit none

    integer, intent(in) :: neq, iTime, NoCalls, Np, nCnstrnts, nSpecs, PrtlNo
    real(8), intent(in) :: aMtrx(nCnstrnts,nSpecs)
    type(hectars), intent(inout) :: htr
    real(8), intent(inout) :: phi(nphi+1,Np)
    integer, intent(in) :: nphi
    real(8), intent(inout) :: massFrc2Prtl(Np,nSpecs) 
    real(8), intent(in) :: CnstrntsAvg(nCnstrnts) 
    real(8), intent(in) :: hAvg, xiAvg 
    real(8), intent(inout) :: gammaPrtl(Np,nCnstrnts)    
    real(8), intent(inout) :: rho(Np)
    real(8), intent(inout) :: Temp(Np)
    real(8), intent(in) :: press
    real(8), intent(in) :: deltaT
    real(8), intent(in) :: RHSVec(Np,nCnstrnts+3) 
    real(8), intent(in) :: yk(nCnstrnts+3)
    real(8), intent(out) :: ki(nCnstrnts+3)
      
    real(8) :: CfcntMtrx(nCnstrnts+3,nCnstrnts+3)
    real(8) :: Answer(nCnstrnts+3)
    integer :: info      
    integer :: i 
    character(80) :: fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,fmt7,fmt8,fmt9,fmt10
    
    call formats(nCnstrnts,nSpecs,Np,fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,fmt7,fmt8,fmt9,fmt10)
        
    call header(iTime, PrtlNo, NoCalls)
    
    call y2GammaRhoTemp(Np,nCnstrnts,PrtlNo,yk,gammaPrtl(PrtlNo,:),rho(PrtlNo),Temp(PrtlNo))
        
  ! Update Variables
    if ((NoCalls == 1) .AND. (iTime == 1)) then
       !!write(*,'(/A, i10)') "First iteration and sub iteration, prtlNo = ", prtlNo
    else
       call update_variables(iTime, NoCalls, Np, nCnstrnts, nSpecs, PrtlNo, aMtrx, htr, &
            phi(:,PrtlNo), nphi, massFrc2Prtl(PrtlNo,:), CnstrntsAvg, hAvg, xiAvg, &
            gammaPrtl(PrtlNo,:), rho(PrtlNo), Temp(PrtlNo), press, RHSVec(PrtlNo,:))
    endif
       
  ! In this subroutine the matrices which are required for RCCE Mixing formulation are calculated
    call RCCECoefficient4_Sub(neq, PrtlNo, Np, nCnstrnts, nSpecs, aMtrx, htr, &
         phi(:,PrtlNo), nphi, massFrc2Prtl, CnstrntsAvg, hAvg, xiAvg, &
         gammaPrtl(PrtlNo,:), rho(PrtlNo), Temp(PrtlNo), press, CfcntMtrx, RHSVec(PrtlNo,:))            
  
  ! In this subroutine ydot dector is calculated  
    call LINPACK_Sub(PrtlNo, Np, nCnstrnts, CfcntMtrx, RHSVec(PrtlNo,:), Answer, info) 

    do i=1,nCnstrnts+3
       ki(i) = Answer(i)
    enddo
    
    write(*, fmt2) "The ydot is:", "Matrix", (Answer(i), i = 1, nCnstrnts+3)
        
end subroutine case_func

!---------------------------------------------------------------------------------------------------------------------------
subroutine Euler(iTime, Np, nCnstrnts, nSpecs, aMtrx, htr, &
           phi, nphi, massFrc2Prtl, CnstrntsAvg, hAvg, xiAvg, &
           gammaPrtl, rho, Temp, press, deltaT, RHSVec)

   implicit none

   integer, intent(in) :: iTime, Np, nCnstrnts, nSpecs
   real(8), intent(in) :: aMtrx(nCnstrnts,nSpecs)
   type(hectars), intent(inout) :: htr
   real(8), intent(inout) :: phi(nphi+1,Np)
   integer, intent(in) :: nphi
   real(8), intent(inout) :: massFrc2Prtl(Np,nSpecs) 
   real(8), intent(in) :: CnstrntsAvg(nCnstrnts) 
   real(8), intent(in) :: hAvg, xiAvg 
   real(8), intent(inout) :: gammaPrtl(Np,nCnstrnts)    
   real(8), intent(inout) :: rho(Np)
   real(8), intent(inout) :: Temp(Np)
   real(8), intent(in) :: press
   real(8), intent(in) :: deltaT
   real(8), intent(inout) :: RHSVec(Np,nCnstrnts+3) 
  
   real(8) :: CfcntMtrx(nCnstrnts+3,nCnstrnts+3)
   real(8) :: Answer(nCnstrnts+3)
   integer :: info      
   integer :: PrtlNo
   integer :: NoCalls
   character(80) :: fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,fmt7,fmt8,fmt9,fmt10
    
   call formats(nCnstrnts,nSpecs,Np,fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,fmt7,fmt8,fmt9,fmt10)
    
   NoCalls = 1
    
   do PrtlNo=1,Np
    
      call header(iTime, PrtlNo, NoCalls)
       
    ! Print variables befor update
      !!call print_variables(iTime, NoCalls, Np, nCnstrnts, nSpecs, PrtlNo, aMtrx, htr, &
      !!     phi(:,PrtlNo), nphi, massFrc2Prtl(PrtlNo,:), CnstrntsAvg, hAvg, xiAvg, &
      !!     gammaPrtl(PrtlNo,:), rho(PrtlNo), Temp(PrtlNo), RHSVec(PrtlNo,:))    
                  
    ! In this part the matrices which are required for RCCE Mixing formulation are calculated.
      call RCCECoefficient3_Sub(PrtlNo, Np, nCnstrnts, nSpecs, aMtrx, htr, &
           phi, nphi, massFrc2Prtl, CnstrntsAvg, hAvg, xiAvg, &
           gammaPrtl, rho, Temp, press, deltaT, CfcntMtrx, RHSVec)
  
      call LINPACK_Sub(PrtlNo, Np, nCnstrnts, CfcntMtrx, RHSVec(PrtlNo,:), Answer, info) 
    
    ! Update Variables
      call y2GammaRhoTemp(Np,nCnstrnts,PrtlNo,Answer,gammaPrtl(PrtlNo,:),rho(PrtlNo),Temp(PrtlNo))
       
      call update_variables(iTime, NoCalls, Np, nCnstrnts, nSpecs, PrtlNo, aMtrx, htr, &
           phi(:,PrtlNo), nphi, massFrc2Prtl(PrtlNo,:), CnstrntsAvg, hAvg, xiAvg, &
           gammaPrtl(PrtlNo,:), rho(PrtlNo), Temp(PrtlNo), press, RHSVec(PrtlNo,:))
   enddo
   
end subroutine Euler

!---------------------------------------------------------------------------------------------------------------------------
subroutine RCCECoefficient4_Sub(neq, PrtlNo, Np, nCnstrnts, nSpecs, aMtrx, htr, phi, nphi, massFrc2, CnstrntsAvg, hAvg, xiAvg, &
           gammaPrtl, rho, Temp, press, CfcntMtrx, RHSVec)

  ! Variables of RCCEMix_Sub
    integer, intent(in) :: neq, PrtlNo
    integer, intent(in) :: Np, nCnstrnts, nSpecs
    real(8), intent(in) :: aMtrx(nCnstrnts,nSpecs)
    type(hectars), intent(inout) :: htr
    real(8), intent(in) :: phi(nphi+1)
    integer, intent(in) :: nphi
    real(8), intent(in) :: massFrc2(nSpecs)
    real(8), intent(in) :: CnstrntsAvg(nCnstrnts) 
    real(8), intent(in) :: hAvg, xiAvg
    real(8), intent(in) :: gammaPrtl(nCnstrnts)     
    real(8), intent(in) :: rho
    real(8), intent(in) :: Temp 
    real(8), intent(in) :: press
    real(8), intent(out) :: CfcntMtrx(neq,neq)
    real(8), intent(in) :: RHSVec(neq)
    
    real(8) :: Cnstrnts(nCnstrnts)
    real(8) :: CTK(nCnstrnts)
    real(8) :: CKM(nCnstrnts,nCnstrnts)
    real(8) :: ET
    real(8) :: BM(nCnstrnts)  
    real(8) :: FK(nCnstrnts)
    real(8) :: FT
    real(8) :: Runi
    real(8) :: MWs(nSpecs)
    real(8) :: MW
    real(8) :: cp
    real(8) :: hsp(nSpecs)
    real(8) :: u(nSpecs)
    real(8) :: g(nSpecs) 
    real(8) :: P0
    integer :: i, ii, j
    character(80) :: fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,fmt7,fmt8,fmt9,fmt10
    
    call formats(nCnstrnts,nSpecs,Np,fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,fmt7,fmt8,fmt9,fmt10)     
    
  ! In this part the matrices which are required for RCCE Mixing formulation are calculated.
     
    call thermo_variables(PrtlNo,Np,nSpecs,nphi,htr,phi,Temp,press,Runi,MWs,MW,cp,hsp,u,g,P0)
     
    do i=1, nCnstrnts
       Cnstrnts(i) = 0
       do j=1, nSpecs
          Cnstrnts(i) = a11aMtrxTGlobal(j,i)*phi(j)/MWs(j) + Cnstrnts(i)
       enddo
    enddo
       
    do i=1, nCnstrnts
       CTK(i) = 0
       do j=1, nSpecs
          CTK(i) = a11aMtrxTGlobal(j,i)*phi(j)*u(j)/(Runi*Temp)/MWs(j) + CTK(i)
       enddo
    enddo       
    
    do i=1, nCnstrnts
    do ii=1, nCnstrnts
       CKM(i,ii) = 0
       do j=1, nSpecs
          CKM(i,ii) = a11aMtrxTGlobal(j,i)*phi(j)*a11aMtrxTGlobal(j,ii)/MWs(j) + CKM(i,ii)
       enddo
    enddo
    enddo

    ET = 0
    do j=1, nSpecs
       ET = hsp(j)*phi(j)*u(j)/(Runi*Temp) + ET
    enddo
    ET = ET + cp*Temp   

    do i=1, nCnstrnts
       BM(i) = 0
       do j=1, nSpecs
          BM(i) = a11aMtrxTGlobal(j,i)*hsp(j)*phi(j) + BM(i)
       enddo
    enddo

    FT = 0
    do j=1, nSpecs
       FT = phi(j)*u(j)/(MWs(j)*Runi*Temp) + FT
    enddo
       
    do i=1, nCnstrnts
       FK(i) = 0
       do j=1, nSpecs
          FK(i) = a11aMtrxTGlobal(j,i)*phi(j)/MWs(j) + FK(i)
       enddo
    enddo
       
    do i=1, nCnstrnts
    do ii=1, nCnstrnts
       CfcntMtrx(i,ii) = - CKM(i,ii)
    enddo
    enddo
       
    do i=1, nCnstrnts
       CfcntMtrx(i,nCnstrnts+1) = - Cnstrnts(i)/rho
    enddo
       
    do i=1, nCnstrnts
       CfcntMtrx(i,nCnstrnts+2) = CTK(i)/Temp
    enddo
       
    do ii=1, nCnstrnts
       CfcntMtrx(nCnstrnts+1,ii) = -BM(ii)
    enddo 

    do ii=1, nCnstrnts
       CfcntMtrx(nCnstrnts+2,ii) = - MW*FK(ii)
    enddo 
    
    call setState_TPY(htr,Temp,press,phi(1:nphi-1)) 
    
    CfcntMtrx(nCnstrnts+1,nCnstrnts+1) = - enthalpy_mass(htr)/rho
    CfcntMtrx(nCnstrnts+1,nCnstrnts+2) = ET/Temp
    CfcntMtrx(nCnstrnts+2,nCnstrnts+1) = 0
    CfcntMtrx(nCnstrnts+2,nCnstrnts+2) = (MW*FT+1)/Temp
    
    do i=1,nCnstrnts+2
       CfcntMtrx(i,nCnstrnts+3) = 0.d0
       CfcntMtrx(nCnstrnts+3,i) = 0.d0
    enddo
    CfcntMtrx(nCnstrnts+3,nCnstrnts+3) = 1.d0
    
    !!if (PrtlNo == 100) then
    !!call print_Coefficients(PrtlNo,Np,nCnstrnts,nSpecs,Cnstrnts,CTK,CKM,BM,FK,FT,ET,CfcntMtrx,RHSVec)
    !!stop
    !!end if
    
    !!if (PrtlNo == 100) then
       !!write(*, fmt4) "Mass fractions", "Matrix", (phi(j), getSpeciesName(htr, j), j = 1, nSpecs) 
       !!stop
    !!endif
 
end subroutine RCCECoefficient4_Sub
!---------------------------------------------------------------------------------------------------------------------------
subroutine RCCECoefficient3_Sub(PrtlNo, Np, nCnstrnts, nSpecs, aMtrx, htr, phi, nphi, massFrc2, CnstrntsAvg, hAvg, xiAvg, &
           gammaPrtl, rho, Temp, press, deltaT, CfcntMtrx, RHSVec)
    

  ! Variables of RCCEMix_Sub
    integer, intent(in) :: PrtlNo
    integer, intent(in) :: Np, nCnstrnts, nSpecs
    real(8), intent(in) :: aMtrx(nCnstrnts,nSpecs)
    type(hectars), intent(inout) :: htr
    real(8), intent(inout) :: phi(nphi+1,Np)
    integer, intent(in) :: nphi
    real(8), intent(inout) :: massFrc2(Np,nSpecs) 
    real(8), intent(in) :: CnstrntsAvg(nCnstrnts) 
    real(8), intent(in) :: hAvg, xiAvg
    real(8), intent(in) :: gammaPrtl(Np,nCnstrnts)     
    real(8), intent(in) :: rho(Np)
    real(8), intent(in) :: Temp(Np) 
    real(8), intent(in) :: press
    real(8), intent(in) :: deltaT
    real(8), intent(out) :: CfcntMtrx(nCnstrnts+3,nCnstrnts+3)
    real(8), intent(inout) :: RHSVec(Np,nCnstrnts+3)
    
    real(8) :: Cnstrnts(nCnstrnts)
    real(8) :: CTK(nCnstrnts)
    real(8) :: CKM(nCnstrnts,nCnstrnts)
    real(8) :: ET
    real(8) :: BM(nCnstrnts)  
    real(8) :: FK(nCnstrnts)
    real(8) :: FT
    real(8) :: Runi
    real(8) :: MWs(nSpecs)
    real(8) :: MW
    real(8) :: cp
    real(8) :: hsp(nSpecs)
    real(8) :: u(nSpecs)
    real(8) :: g(nSpecs) 
    real(8) :: P0
    integer :: i, ii, j
    character(80) :: fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,fmt7,fmt8,fmt9,fmt10
    
    call formats(nCnstrnts,nSpecs,Np,fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,fmt7,fmt8,fmt9,fmt10)     
    
  ! In this part the matrices which are required for RCCE Mixing formulation are calculated.
     
    call thermo_variables(PrtlNo,Np,nSpecs,nphi,htr,phi(:,PrtlNo),Temp(PrtlNo),press,Runi,MWs,MW,cp,hsp,u,g,P0)
        
  ! In this part the matrices which are required for RCCE Mixing formulation are calculated.
  
    do i=1, nCnstrnts
       Cnstrnts(i) = 0
       do j=1, nSpecs
          Cnstrnts(i) = a11aMtrxTGlobal(j,i)*massFrc2(PrtlNo,j)/MWs(j) + Cnstrnts(i)
       enddo
    enddo
     
    do i=1, nCnstrnts
       CTK(i) = 0
       do j=1, nSpecs
          CTK(i) = a11aMtrxTGlobal(j,i)*massFrc2(PrtlNo,j)*u(j)/(Runi*Temp(PrtlNo))/MWs(j) + CTK(i)
       enddo
    enddo       

    do i=1, nCnstrnts
    do ii=1, nCnstrnts
       CKM(i,ii) = 0
       do j=1, nSpecs
          CKM(i,ii) = a11aMtrxTGlobal(j,i)*massFrc2(PrtlNo,j)*a11aMtrxTGlobal(j,ii)/MWs(j) + CKM(i,ii)
       enddo
    enddo
    enddo

    ET = 0
    do j=1, nSpecs
       ET = hsp(j)*massFrc2(PrtlNo,j)*u(j)/(Runi*Temp(PrtlNo)) + ET
    enddo
    ET = ET + cp*Temp(PrtlNo)   

    do i=1, nCnstrnts
       BM(i) = 0
       do j=1, nSpecs
          BM(i) = a11aMtrxTGlobal(j,i)*hsp(j)*massFrc2(PrtlNo,j) + BM(i)
       enddo
    enddo

    FT = 0
    do j=1, nSpecs
       FT = massFrc2(PrtlNo,j)*u(j)/(MWs(j)*Runi*Temp(PrtlNo)) + FT
    enddo
       
    do i=1, nCnstrnts
       FK(i) = 0
       do j=1, nSpecs
          FK(i) = a11aMtrxTGlobal(j,i)*massFrc2(PrtlNo,j)/MWs(j) + FK(i)
       enddo
    enddo
       
    do i=1, nCnstrnts
    do ii=1, nCnstrnts
       CfcntMtrx(i,ii) = - CKM(i,ii)
    enddo
    enddo
       
    do i=1, nCnstrnts
       CfcntMtrx(i,nCnstrnts+1) = - Cnstrnts(i)/rho(PrtlNo)
    enddo
       
    do i=1, nCnstrnts
       CfcntMtrx(i,nCnstrnts+2) = CTK(i)/Temp(PrtlNo)
    enddo
       
    do ii=1, nCnstrnts
       CfcntMtrx(nCnstrnts+1,ii) = -BM(ii)*1.0d-20
    enddo 
     
    do ii=1, nCnstrnts
       CfcntMtrx(nCnstrnts+2,ii) = - MW*FK(ii)
    enddo 
     
    call setState_TPY(htr,Temp(PrtlNo),press,phi(1:nphi-1,PrtlNo))
      
    CfcntMtrx(nCnstrnts+1,nCnstrnts+1) = - enthalpy_mass(htr)/rho(PrtlNo)*1.0d-20
    CfcntMtrx(nCnstrnts+1,nCnstrnts+2) = ET/Temp(PrtlNo)*1.0d-20
    CfcntMtrx(nCnstrnts+2,nCnstrnts+1) = 0
    CfcntMtrx(nCnstrnts+2,nCnstrnts+2) = (MW*FT+1)/Temp(PrtlNo)

  ! RHSVec
    do i=1,nCnstrnts
       RHSVec(PrtlNo,i) = RHSVec(PrtlNo,i)*deltaT
       do ii=1,nCnstrnts
          RHSVec(PrtlNo,i) = RHSVec(PrtlNo,i) - CKM(i,ii)*gammaPrtl(PrtlNo,ii)
       enddo
       RHSVec(PrtlNo,i) = RHSVec(PrtlNo,i) - Cnstrnts(i) + CTK(i)
    enddo
    
    
    RHSVec(PrtlNo,nCnstrnts+1) = RHSVec(PrtlNo,nCnstrnts+1)/1.0d-20*deltaT
    do ii=1,nCnstrnts
       RHSVec(PrtlNo,nCnstrnts+1) = RHSVec(PrtlNo,nCnstrnts+1) - BM(ii)*gammaPrtl(PrtlNo,ii)
    enddo    
    RHSVec(PrtlNo,nCnstrnts+1) = RHSVec(PrtlNo,nCnstrnts+1) - enthalpy_mass(htr) + ET
    RHSVec(PrtlNo,nCnstrnts+1) = RHSVec(PrtlNo,nCnstrnts+1)*1.0d-20
    
    RHSVec(PrtlNo,nCnstrnts+2) = RHSVec(PrtlNo,nCnstrnts+2)*deltaT
    do ii=1,nCnstrnts
       RHSVec(PrtlNo,nCnstrnts+2) = RHSVec(PrtlNo,nCnstrnts+2) - MW*FK(ii)*gammaPrtl(PrtlNo,ii)
    enddo    
    RHSVec(PrtlNo,nCnstrnts+2) = RHSVec(PrtlNo,nCnstrnts+2) + (MW*FT+1)
   
    !!call print_Coefficients(PrtlNo,Np,nCnstrnts,nSpecs,Cnstrnts,CTK,CKM,BM,FK,FT,ET,CfcntMtrx,RHSVec)
     
end subroutine RCCECoefficient3_Sub

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
   character(80) :: fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,fmt7,fmt8,fmt9,fmt10
   
   call formats(nCnstrnts,nSpecs,Np,fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,fmt7,fmt8,fmt9,fmt10)
      
   do i=1,nCnstrnts
      do ii=1,nCnstrnts
	     mtrxIn(i,ii) = aMtrx(i,ii)
	  end do
   end do
   
   call MIGS(mtrxIn,nCnstrnts,mtrxInv,INDX)
	
   call MATRIXPRODUCT(mtrxInv, nCnstrnts, nCnstrnts, aMtrx, nCnstrnts, nSpecs, a11aMtrx, rowsC, colsC, ErrCode)
	
 ! For the sake of print only
   do j=1,nSpecs
	  do i=1,nCnstrnts
	     a11aMtrxT(j,i)= a11aMtrx(i,j)
	  end do
   end do

   !!write(*, fmt5) "a11aMtrxT is:", "Matrix", ((a11aMtrxT(i, j), j = 1, size(a11aMtrxT, 2)), i = 1, size(a11aMtrxT, 1))

end subroutine Calc_a11aMtrx
    
!--------------------------------------------------------------------------------------------------------------------
subroutine RHS(neq, PrtlNo, iTime, Np, nCnstrnts, nSpecs, aMtrx, &
           htr, phi, nphi, massFrc2Prtl, rho, Temp, press, CnstrntsAvg, hAvg, xiAvg, RHSVec, TauMix)
   
   use mod_hectars
   implicit none

   integer, intent(in) :: neq, PrtlNo, iTime, Np, nCnstrnts, nSpecs
   real(8), intent(in) :: aMtrx(nCnstrnts,nSpecs)
   type(hectars), intent(inout) :: htr
   real(8), intent(in) :: phi(nphi+1)
   integer, intent(in) :: nphi
   real(8), intent(in) :: massFrc2Prtl(nSpecs) 
   real(8), intent(in) :: rho
   real(8), intent(in) :: Temp
   real(8), intent(in) :: press
   real(8), intent(in) :: CnstrntsAvg(nCnstrnts) 
   real(8), intent(in) :: hAvg, xiAvg    
   real(8), intent(out) :: RHSVec(neq) 
   real(8), intent(in) :: TauMix
   
   real(8) :: Ch
   real(8) :: a11aMtrx(nCnstrnts, nSpecs)
   real(8) :: Cnstrnts(nCnstrnts)
   real(8) :: MWs(nSpecs)     
   real(8) :: h
   real(8) :: wdotPrtl(nSpecs)   
   integer :: i,j
   character(80) :: fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,fmt7,fmt8,fmt9,fmt10

   real(8) :: t1,t2
   logical :: ifPrint
   
   ifPrint = .false.
   if ((iTime == 1) .OR. (iTime == 50) .OR. (iTime == 100) .OR. (iTime == 150)) then
   if (PrtlNo == 1) then
      ifPrint = .true.
   endif
   endif
          
   call formats(nCnstrnts,nSpecs,Np,fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,fmt7,fmt8,fmt9,fmt10)

   !!call Calc_a11aMtrx(Np,nCnstrnts,nSpecs,aMtrx,a11aMtrx)
   
   do j=1,nSpecs
      do i=1,nCnstrnts
         a11aMtrx(i,j) = a11aMtrxTGlobal(j,i)
      end do
   end do
    
   !! call get_MWs(cpr,MWs)
   !! 1.
   call getMolecularWeights(htr,MWs) ! SI units, kg/mole kg/kmole?

   Ch = 1.d0/TauMix
       
      do i=1, nCnstrnts
         Cnstrnts(i) = 0
         do j=1, nSpecs
            Cnstrnts(i) = a11aMtrx(i,j)*phi(j)/MWs(j) + Cnstrnts(i)
         enddo
      enddo
      
      call setState_TPY(htr,Temp,press,phi(1:nphi-1))
      
      !! h = get_h(cpr,phi)
      !! 2.
      h = enthalpy_mass(htr)
    
      do i=1,nCnstrnts
         RHSVec(i) = 0.d0 !- Ch*(Cnstrnts(i)-CnstrntsAvg(i))
      enddo
      RHSVec(nCnstrnts+1) = 0.d0 !- Ch*( h - hAvg)
      RHSVec(nCnstrnts+2) = 0.d0

   call cpu_time ( t1 )
    
   call RHSChem(PrtlNo, iTime, Np, nCnstrnts, nSpecs, htr, phi, nphi, press, wdotPrtl)

   call cpu_time ( t2 )
   if (ifPrint) then
      write ( *, * ) 'Elapsed CPU time for RHSChem = ', t2-t1
   endif

   do i=1,nCnstrnts
      do j=1,nSpecs
         RHSVec(i) = RHSVec(i) + a11aMtrx(i,j)*wdotPrtl(j)/rho
      enddo
   enddo
   
   RHSVec(nCnstrnts+3) = 0.d0!- Ch*( phi(nphi+1) - xiAvg )
      
end subroutine RHS
    
!----------------------------------------------------------------------------------------------------------------------
subroutine RHSChem(PrtlNo, iTime, Np, nCnstrnts, nSpecs, htr, phi, nphi, press, wdotPrtl)

   use mod_hectars
   !! use mod_chemrate
   implicit none
      
   integer, intent(in) :: PrtlNo, iTime, Np, nCnstrnts, nSpecs
   type(hectars), intent(inout) :: htr
   real(8), intent(in) :: phi(nphi+1)
   integer, intent(in) :: nphi
   real(8), intent(in) :: press
   real(8), intent(out) :: wdotPrtl(nSpecs)
   
   integer :: j
   real(8) :: wdot(nphi)
   character(80) :: fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,fmt7,fmt8,fmt9,fmt10

      
   call formats(nCnstrnts,nSpecs,Np,fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,fmt7,fmt8,fmt9,fmt10)
   
   !! Note: the state has been set in RHS :-)
   
 ! wdot chemical production rates of the species 
 ! cgs units, moles/(cm**3*sec)

   !! write(*, fmt4) "RCCE State in RHSChem:", "Matrix", (phi(j,1), get_species_name(cpr, j), j = 1, nSpecs)   
   !! call get_wdot( cpr, phi(:,k), wdot )
   !! 1.
   call get_wdot_hectars( htr, phi, press, wdot ) ![moles/m3*sec]
   do j=1,nSpecs
      wdotPrtl(j) = wdot(j)
   enddo
   !!write(*, '(/A,i4)') "PrtlNo = ", k
   !!write(*, fmt4) "wdot:", "Matrix", (wdot(j), get_species_name(cpr, j), j = 1, nSpecs)
   
end subroutine RHSChem

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
    real(8), intent(in) :: RHSVec(nCnstrnts+3) 
    
!    real(8) :: tolRT
!    real(8) :: RTOld
!    real(8) :: RTNew
    real(8) :: tolT
    real(8) :: TOld
    real(8) :: TNew
    integer :: maxIter
    character(80) :: fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,fmt7,fmt8,fmt9,fmt10

      
    call formats(nCnstrnts,nSpecs,Np,fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,fmt7,fmt8,fmt9,fmt10)
       
  ! Print variables befor update
!    call print_variables(iTime, NoCalls, Np, nCnstrnts, nSpecs, PrtlNo, aMtrx, htr, &
!         phi, nphi, massFrc2Prtl, CnstrntsAvg, hAvg, xiAvg, &
!         gammaPrtl, rho, Temp, press, RHSVec)
    
  ! Update Variables
    tolT = 1.0d-5
    maxIter = 0
    
    !! do while ( (tolRT>1.0d-6) .AND. (maxIter<5))
    !! 1.
    do while ( (tolT>1.0d-6) .AND. (maxIter<1))
       
       !! call Calc_RT(Np,PrtlNo,htr,phi,nphi,Temp,RTOld)
       !! 2.
       !!call setState_TPY(htr,Temp,press,phi(1:nphi-1))
       !!TOld = enthalpy_mass(htr)
           
     ! Update mass fractions
       call update_massFrcs(iTime, NoCalls, Np, nCnstrnts, nSpecs, PrtlNo, aMtrx, htr, &
            phi, nphi, massFrc2Prtl, CnstrntsAvg, hAvg, xiAvg, &
            gammaPrtl, rho, Temp, press, RHSVec)

       !!write(*, fmt4) "The massFrc matrix in update_variables is :", &
       !!"Matrix", (phi(j), get_species_name(cpr, j), j = 1, nSpecs) 
       !!write(*,'(/A, e17.10 )') "Temp in update_variables is :", Temp
       !!write(*,'(/A,e17.10)') "get_T(cpr,phi) in update_variables is :", &
       !!get_T(cpr,phi)          
       
       !! call Calc_RT(Np,PrtlNo,htr,phi,nphi,Temp,RTNew)
       !! 3.
       call setState_TPY(htr,Temp,press,phi(1:nphi-1))
       !!TNew = enthalpy_mass(htr)
       
     ! Update RT 
       !! phi(nphi) = RTNew
       !! 4.     
       phi(nphi) = enthalpy_mass(htr)
    
       !! call tol_RT(RTOld, RTNew, tolRT)
       !! 5.
       !!call tol_RT(TOld, TNew, tolT)
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
     
end subroutine

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
    real(8), intent(in) :: RHSVec(nCnstrnts+3)
           
    real(8) :: Runi
    real(8) :: MWs(nSpecs)
    real(8) :: MW
    real(8) :: cp
    real(8) :: hsp(nSpecs)
    real(8) :: u(nSpecs)    
    real(8) :: g(nSpecs) 
    real(8) :: P0
    real(8) :: a11aMtrxTGama(nSpecs)
    real(8) :: massSum
    integer :: j
        
  ! Update mass fractions
    
    call thermo_variables(PrtlNo,Np,nSpecs,nphi,htr,phi,Temp,press,Runi,MWs,MW,cp,hsp,u,g,P0)
    
    a11aMtrxTGama = 0
    a11aMtrxTGama = MATMUL(a11aMtrxTGlobal,gammaPrtl)
       
    do j=1,nSpecs
!       phi(j) = MWs(j)/rho*P0/(Runi*Temp)* &
!                               exp(-g(j)/(Runi*Temp)-a11aMtrxTGama(j))
       phi(j)=P0/Press*MWs(j)/MW*exp(-a11aMtrxTGama(j)-g(j)/(Runi*Temp))
                               
       if (phi(j)<1.0d-16) then
          phi(j) = 1.0d-16
       endif
       if (phi(j)>1.0) then
          phi(j) = 1.0d0
       endif       
    enddo 
      
    !!write(*,'(/A,i5,A,i5)') "NoCalls = ", NoCalls, "   PrtlNo = ", PrtlNo
    
!    call print_speciesVariables(iTime, NoCalls, 2, Np, nCnstrnts, nSpecs, PrtlNo, 16, &
!         aMtrx, htr, phi, nphi, massFrc2Prtl, &
!         CnstrntsAvg, hAvg, xiAvg, gammaPrtl, rho, Temp, press, MWs, Runi, g, a11aMtrxTGama)       
  
end subroutine update_massFrcs 
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
subroutine y2GammaRhoTemp(Np,nCnstrnts,PrtlNo,yk,gammaPrtl,rho,Temp)

    implicit none

    integer, intent(in) :: Np, nCnstrnts, PrtlNo
    real(8), intent(in) :: yk(nCnstrnts+3)
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
    real(8), intent(out) :: yk(nCnstrnts+3)
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
    real(8), intent(in) :: RHSVec(nCnstrnts+3)
    
    integer :: i
    character(80) :: fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,fmt7,fmt8,fmt9,fmt10
    
    call formats(nCnstrnts,nSpecs,Np,fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,fmt7,fmt8,fmt9,fmt10)
    
    call setState_TPY(htr,Temp,press,phi(1:nphi-1))
    
    write(*,'(/A, i5 )') "PrtlNo = ", PrtlNo
    write(*,'(/A, i5 )') "NoCalls = ", NoCalls
    write(*,'(/A, e17.10 )') "Temp  = ", Temp
    write(*,'(/A,e17.10)') "temperature(htr)  =", temperature(htr)
    write(*,'(/A,e17.10,A,e17.10)') "phi(nphi =", phi(nphi), &
          "  enthalpy_mass(htr) =", enthalpy_mass(htr)
    write(*,'(/A, e17.10 )') "rho  = ", rho
    write(*, fmt1) "The gammaPrtl is:", "Matrix", (gammaPrtl(i), i = 1, nCnstrnts)
    write(*, fmt4) "Mass fractions are :", "Matrix", (phi(i), &
                   getSpeciesName(htr, i), i = 1, nSpecs)
    write(*,'(/A, e17.10 )') "Runi  = ", Runiv(htr)
    !!write(*,'(/A,e17.10)') "get_MW(cpr,phi) =", get_MW(cpr,phi)
    write(*, fmt2) "The RHSVec matrix is:", "Matrix", (RHSVec(i), i = 1, nCnstrnts+2)
    
end subroutine print_variables

!----------------------------------------------------------------------------------------------------------------------
subroutine print_speciesVariables(iTime, NoCalls, iTimeP, Np, nCnstrnts, nSpecs, PrtlNo, PrtlNoP, aMtrx, htr, &
           phi, nphi, massFrc2Prtl, &
           CnstrntsAvg, hAvg, xiAvg, gammaPrtl, rho, Temp, press, MWs, Runi, g, a11aMtrxTGama)
       
    integer, intent(in) :: iTime, NoCalls, iTimeP, Np, nCnstrnts, nSpecs, PrtlNo, PrtlNoP
    real(8), intent(in) :: aMtrx(nCnstrnts,nSpecs)
    type(hectars), intent(inout) :: htr ! Modify later to intent(in)
    real(8), intent(in) :: phi(nphi+1)
    integer, intent(in) :: nphi
    real(8), intent(in) :: massFrc2Prtl(nSpecs) 
    real(8), intent(in) :: CnstrntsAvg(nCnstrnts) 
    real(8), intent(in) :: hAvg, xiAvg 
    real(8), intent(in) :: gammaPrtl(nCnstrnts)    
    real(8), intent(in) :: rho
    real(8), intent(in) :: Temp
    real(8), intent(in) :: press
    real(8), intent(in) :: MWs(nSpecs), Runi
    real(8), intent(in) :: g(nSpecs) 
    real(8), intent(in) :: a11aMtrxTGama(nSpecs)
    
    integer :: i,j,jj
    character(80) :: fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,fmt7,fmt8,fmt9,fmt10
    
    call formats(nCnstrnts,nSpecs,Np,fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,fmt7,fmt8,fmt9,fmt10)
    
    call setState_TPY(htr,Temp,press,phi(1:nphi-1))
    
    if ((NoCalls == iTimeP) .AND. (PrtlNo == PrtlNoP)) then
    do j=1, nSpecs 
    
       print '(/ "----------------------------------------------------------------------------------")'
       print '(/ "----------------------------------------------------------------------------------")'
       write(*,'(/A, i4 )') "Iteration number = ", iTime
       write (*,'(/A, I4)') "Particle number: ", PrtlNo
       write (*,'(/A, I4)') "Number of calls: ", NoCalls
       write (*,'(/A, I4)') "Species Number: ", j
       write(*,'(/A,A)') "getSpeciesName(htr, j) =", getSpeciesName(htr, j)
       write(*,'(/A,e17.10)') "MWs(j) =", MWs(j)
       write(*,'(/A,e17.10)') "rho =", rho
       write(*,'(/A,e17.10)') "Runi =", Runi
       write(*,'(/A,e17.10)') "Temp =", Temp
       write(*,'(/A,e17.10)') "Temperature from hectars =", temperature(htr)
       write(*,'(/A,e17.10)') "g(j) =", g(j)
       write(*,'(/A,e17.10)') "-g(j)/(Runi*Temp) =", -g(j)/(Runi*Temp)
       write(*,'(/A,e17.10)') "exp(-g(j)/(Runi*Temp)) =", exp(-g(j)/(Runi*Temp))
       write(*,'(/A,e17.10)') "-a11aMtrxTGama(j) =", -a11aMtrxTGama(j)
       write(*,'(/A,e17.10)') "exp(-a11aMtrxTGama(j)) =", exp(-a11aMtrxTGama(j))
       write(*,'(/A,e17.10)') "phi(j) =", phi(j)
       write(*, fmt1) "gammaPrtl is:", "Matrix", (gammaPrtl(i), i = 1, nCnstrnts )
       write(*, fmt5) "a11aMtrxTGlobal is:", "Matrix", ((a11aMtrxTGlobal(jj,i), i = 1, nCnstrnts), jj = 1, nSpecs) 
       write(*, fmt4) "a11aMtrxTGama is:", "Matrix", (a11aMtrxTGama(jj), getSpeciesName(htr, jj), jj = 1, nSpecs) 
       write(*,'(/A,i4,A,i4,A,i4/)') "Number of rows", size(a11aMtrxTGlobal,1), "  Number of columns", &
          size(a11aMtrxTGlobal,2), "  Number of rows",size(gammaPrtl)
    
    enddo
    endif
    
end subroutine print_speciesVariables

!---------------------------------------------------------------------------------------------------------------------------
subroutine print_Coefficients(PrtlNo,Np,nCnstrnts,nSpecs,Cnstrnts,CTK,CKM,BM,FK,FT,ET,CfcntMtrx,RHSVec)
    
    integer, intent(in) :: PrtlNo,Np,nCnstrnts,nSpecs
    real(8), intent(in) :: Cnstrnts(nCnstrnts)
    real(8), intent(in) :: CTK(nCnstrnts)
    real(8), intent(in) :: CKM(nCnstrnts,nCnstrnts)
    real(8), intent(in) :: BM(nCnstrnts)  
    real(8), intent(in) :: FK(nCnstrnts)
    real(8), intent(in) :: FT
    real(8), intent(in) :: ET    
    real(8), intent(in) :: CfcntMtrx(nCnstrnts+3,nCnstrnts+3)
    real(8), intent(in) :: RHSVec(nCnstrnts+3)
    
    integer :: i,j
    character(80) :: fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,fmt7,fmt8,fmt9,fmt10
    
    call formats(nCnstrnts,nSpecs,Np,fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,fmt7,fmt8,fmt9,fmt10)
        
  ! In this part the matrices which are required for RCCE Mixing formulation are printed.
  ! Particle # PrtlNo matices  
    
  ! Printing
    write(*, fmt1) "The Cnstrnts matrix is:", "Matrix", (Cnstrnts(i), i = 1, nCnstrnts)
    write(*, fmt1) "The CTK matrix is:", "Matrix", (CTK(i), i = 1, nCnstrnts)
    write(*, fmt6) "The CKM matrix is:", "Matrix", ((CKM(i, j), j = 1, nCnstrnts), i = 1, nCnstrnts)
    write(*, fmt1) "The BM matrix is:", "Matrix", (BM(i), i = 1, nCnstrnts)
    write(*, fmt1) "The FK matrix is:", "Matrix", (FK(i), i = 1, nCnstrnts)
    write(*, fmt7) "The CfcntMtrx is:", "Matrix", ((CfcntMtrx(i, j), j = 1, nCnstrnts+2), i = 1, nCnstrnts+2)    
    write(*, fmt2) "The RHSVec matrix is:", "Matrix", (RHSVec(i), i = 1, nCnstrnts+2)
    
end subroutine print_Coefficients

!---------------------------------------------------------------------------------------------------------------------------
subroutine formats(nCnstrnts,nSpecs,Np,fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,fmt7,fmt8,fmt9,fmt10)
integer, intent(in) :: nCnstrnts, nSpecs, Np
character(80), intent(out) :: fmt1,fmt2,fmt3,fmt4,fmt5,fmt6,fmt7,fmt8,fmt9,fmt10
    
    write (fmt1, '("(/ A / A, "" = ["", ", i0, "( ""["",", i0, "(1x,d16.9, 1x), ""]"" / 10x), ""]"" )")') nCnstrnts, 1
    write (fmt2, '("(/ A / A, "" = ["", ", i0, "( ""["",", i0, "(1x,1x,d16.9, 1x), ""]"" / 10x), ""]"" )")') nCnstrnts+3, 1
    write (fmt3, '("(/ A / A, "" = ["", ", i0, "( ""["",", i0, "(1x,d16.9, 1x), ""]"" / 10x), ""]"" )")') nSpecs, 1
    write (fmt4, '("(/ A / A, "" = ["", ", i0, "( ""["",", i0, "(1x,d16.9, 1x), ""]"", 1x, A / 10x), ""]"" )")') nSpecs, 1
    write (fmt5, '("(/ A / A, "" = ["", ", i0, "( ""["",", i0, "(1x,F10.2, 1x), ""]"" / 10x), ""]"" )")') nSpecs, nCnstrnts
    write (fmt6, '("(/ A / A, "" = ["", ", i0, "( ""["",", i0, "(1x,d16.9, 1x), ""]"" / 10x), ""]"" )")') nCnstrnts, nCnstrnts
    write (fmt7, '("(/ A / A, "" = ["", ", i0, "( ""["",", i0, "(1x,d10.2, 1x), ""]"" / 10x), ""]"" )")') nCnstrnts+3, nCnstrnts+3
    write (fmt8, '("(/ A / A, "" = ["", ", i0, "( ""["",", i0, "(1x,d16.9, 1x), ""]"" / 10x), ""]"" )")') Np, 1
    write (fmt9, '("(/ A / A, "" = ["", ", i0, "( ""["",", i0, "(1x,i2, 1x), ""]"", 1x, i5, 2x, 10A / 10x), ""]"" )")') Np, 1
    write (fmt10, '("(/ A / A, "" = ["", ", i0, "( ""["",", i0, "(1x,F10.2, 1x), ""]"" / 10x), ""]"" )")') nSpecs, nCnstrnts
         
end subroutine formats

!---------------------------------------------------------------------------------------------------------------------------
subroutine header(iTime, PrtlNo, NoCalls)
integer, intent(in) :: iTime, PrtlNo, NoCalls
        
    print '(/ "------------------------------------------------------------------------------------------------------------")'
    write(*,'(/A, I10 )') "Iteration number = ", iTime
    write (*,'(/A, I4)') "Particle number = ", PrtlNo
    write (*,'(/A, I10)') "Number of calls in iteration number = ", NoCalls
    
end subroutine header

!---------------------------------------------------------------------------------------------------------------------------
end module RCCE_Mod
