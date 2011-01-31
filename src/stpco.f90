module stpcomod

!$$$ module documentation block
!           .      .    .                                       .
! module:   stpcomod    module for stpco
!  prgmmr:
!
! abstract: module for stpco
!
! program history log:
!   2010-06-16  todling - based on stpco
!
! subroutines included:
!   sub stpco
!   sub stpcolay_
!   sub stpcolev_
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

implicit none

PRIVATE
PUBLIC stpco

contains

!subroutine stpco(cohead,co3lhead,rval,sval,out,sges,nstep)
subroutine stpco(co3lhead,rval,sval,out,sges,nstep)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    stpco       call components to calculate contrib. to
!                            penalty and stepsize for carbon monoxide
!   prgmmr: sienkiewicz     org: GMAO                 date: 2009-01-22
!
! abstract: The routine calls individual components that calculate 
!           contribution to the penalty and step size from layer 
!           and level carbon monoxine measurements
!
! program history log:
!   2010-06-14  todling - based on Sienkiewicz carbon monoxine code
!
!   input argument list:
!     cohead
!     co3lhead
!     rval - search direction for carbon monoxine
!     sval - input carbon monoxine correction field
!     sges - step size estimates (nstep)
!     nstep- number of stepsize estimates (==0 means use outer iteration value)
!
!   output argument list:
!     out(1:nstep) - contribution of carbon monoxine data to penalty sges(1:nstep)
!
! attributes:
!   language: f90
!   machine:
!
!$$$  
  use kinds, only: r_kind,r_quad,i_kind
! use obsmod, only: co_ob_type,co3l_ob_type
  use obsmod, only: co3l_ob_type
  use gridmod, only: latlon1n
  use constants, only: zero_quad
  use gsi_bundlemod, only: gsi_bundle
  implicit none

! Declare passed variables

! type(  co_ob_type),pointer          ,intent(in   ) :: cohead
  type(co3l_ob_type),pointer          ,intent(in   ) :: co3lhead
  integer(i_kind)                     ,intent(in   ) :: nstep
  type(gsi_bundle)                    ,intent(in   ) :: sval
  type(gsi_bundle)                    ,intent(in   ) :: rval
  real(r_kind),dimension(max(1,nstep)),intent(in   ) :: sges
  real(r_quad),dimension(max(1,nstep)),intent(inout) :: out

  out=zero_quad

! call stpcolay_(  cohead,rval,sval,out,sges,nstep)
! call stpcolev_(co3lhead,rval,sval,out,sges,nstep)

  return

end subroutine stpco

end module stpcomod