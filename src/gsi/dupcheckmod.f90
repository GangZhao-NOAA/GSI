module dupcheckmod
!$$$  module documentation block
!                .      .    .                                       .
! module:    dupcheckmod
!   prgmmr: pondeca/morris   org: ncep/emc           date: 2026-09-16
!
! abstract: module containing duplicate observation check code for
!           conventional observations
!
! program history log:
!   2026-09-16  pondeca/morris - original code, extracted from setup routines
!
! subroutines included:
!   def dupcheck - perform duplicate observation check
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$

  implicit none
  private
  public :: dupcheck

contains

subroutine dupcheck(nobs, nele, data, muse, dup, &
                    ier, itime, ilate, ilone, id, &
                    l_closeobs, min_offset, &
                    ipres, skip_pres_check, &
                    ikxx,  skip_otyp_check, &
                    l_print)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    dupcheck  perform duplicate observation check
!   prgmmr: pondeca/morris  org: ncep/emc           date: 2026-09-16
!
! abstract: check for duplicate observations at same location using
!           true lat/lon values (ilate/ilone) and station id matching.
!           duplicate stations can have lat/lon specs differing by as
!           much as epsdup (~0.005 deg).  a station can also appear as
!           a TAC station and BUFR station with lat/lon specs differing
!           by as much as epsdup_2 (~0.1 deg); the station id match
!           option handles this situation.  with l_closeobs, only the
!           observation closest to the analysis time is retained;
!           otherwise a duplicate weighting scheme is applied.
!           an optional pressure-level matching check can be applied
!           per observation via the skip_pres_check argument.
!
! program history log:
!   2026-09-16  pondeca/morris - original code extracted from setup routines
!
!   input argument list:
!     nobs     - number of observations
!     nele     - number of data elements per observation
!     data     - observation data array (nele,nobs)
!     muse     - logical array indicating if obs is used (intent inout)
!     dup      - duplicate weight array (intent inout)
!     ier      - index of observation error in data array
!     itime    - index of observation time in data array
!     ilate    - index of observation latitude (degrees) in data array
!     ilone    - index of observation longitude (degrees) in data array
!     id       - index of station id in data array
!     l_closeobs - if .true., retain only obs closest to analysis time
!     min_offset - offset of analysis time from reference time (minutes)
!     ipres    - (optional) index of pressure in data array
!     skip_pres_check - (optional) logical array of length nobs;
!                       if element k is .true., skip pressure check
!                       for obs k; defaults to .true. for all obs
!                       (i.e., no pressure check applied by default)
!     ikxx     - (optional) index of observation k-type (100~999) in data array
!     skip_otyp_check - (optional) logical array of length nobs;
!                       if element k is .true., skip obs-type check
!                       for obs k; defaults to .true. for all obs
!                       (i.e., no obs-type check applied by default)
!     l_print  - (optional) logical, in .true., print out information for debugging 
!                       default is .false. locally
!
!   output argument list:
!     muse     - updated logical array indicating if obs is used
!     dup      - updated duplicate weight array
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$

  use kinds, only: r_kind, i_kind, r_double
  use constants, only: one, r1000
  use qcmod, only: dfact, dfact1, epsdup, epsdup_2
! epsdup and epsdup_2 are initialized to zero in qcmod and may be set via
! the OBSQC namelist group (e.g., epsdup~0.005 deg, epsdup_2~0.1 deg)

  integer(i_kind), intent(in)    :: nobs, nele
  integer(i_kind), intent(in)    :: ier, itime, ilate, ilone, id
  integer(i_kind), intent(in)    :: min_offset
  real(r_kind),    intent(in)    :: data(nele,nobs)
  logical,         intent(inout) :: muse(nobs)
  real(r_kind),    intent(inout) :: dup(nobs)
  logical,         intent(in)    :: l_closeobs
  integer(i_kind), optional, intent(in) :: ipres
  logical,         optional, intent(in) :: skip_pres_check(:)
  integer(i_kind), optional, intent(in) :: ikxx
  logical,         optional, intent(in) :: skip_otyp_check(:)
  logical,         optional, intent(in) :: l_print

  integer(i_kind)            :: i, k, l, nlen, nlen2, ipres_loc, ikxx_loc
  real(r_kind)               :: hr_offset, tfact
  real(r_double)             :: rstn
  character(len=8)           :: cstn, cstn_tmp
  character(len=8), dimension(nobs) :: stnid_main
  character(len=1), parameter :: cblank = ' '
  logical :: duplogic, duplogic_1, duplogic_2, apply_pres_k
  logical :: l_print_loc

  equivalence(rstn, cstn)

  if (present(ipres)) then
     ipres_loc = ipres
  else
     ipres_loc = 0
  end if
  if (present(ikxx)) then
     ikxx_loc = ikxx
  else
     ikxx_loc = 0
  end if
  if (present(l_print)) then
     l_print_loc = l_print
  else
     l_print_loc = .false.
  end if

  hr_offset = real(min_offset/60.0_r_kind,r_kind)
  dup = one

! Get the main station ID (the string prededing the first blank in original station ID string)  
! Purpose: for mesonet station ids of the form "STNIxxxa", the trailing character "a" 
!          in position 8 (preceded by blanks) is intentionally excluded from the comparison.
  do k=1,nobs
     rstn = data(id,k)
!    cstn = transfer(rstn, cstn)     ! equivalence is obsolete in modern fortran
     cstn_tmp = adjustl(cstn)        ! removing the heading space(s)
!    if first blank was detected, set the following element(s) be blank 
     do i=1,7
        if (cstn_tmp(i:i) == cblank) then
           cstn_tmp(i+1:i+1) = cblank
        end if
     end do
     stnid_main(k) = cstn_tmp
     if (l_print_loc) then
        write(6,'(1x,A9,A8,A4,A8,4x,I8.8)') 'ID Check:', adjustl(cstn), '===>', stnid_main(k), k
     end if
  end do

  kloop: do k=1,nobs

!    Determine if pressure check applies for this observation
!    (requires both ipres and skip_pres_check to be present)
    if (present(skip_pres_check) .and. present(ipres)) then
       apply_pres_k = .not. skip_pres_check(k)
    else
       apply_pres_k = .false.
    end if

    lloop: do l=k+1,nobs

       duplogic_1=abs(data(ilate,k)-data(ilate,l))<=epsdup .and.  &  !duplicate stations can have lat/lon specs
                  abs(data(ilone,k)-data(ilone,l))<=epsdup           !differing by as much as epsdup (~0.005 deg)

       duplogic_2=abs(data(ilate,k)-data(ilate,l))<=epsdup_2 .and.  & !station can appear as TAC station and BUFR station
                  abs(data(ilone,k)-data(ilone,l))<=epsdup_2 .and.  & !with lat/lon specs differing by as much as epsdup_2 (~0.1 deg)
                  stnid_main(k)==stnid_main(l)                        !this logic addresses this situation, but only when the station ids
                                                                      !are the same. when they are different, the duplicate obs will slip in

       if (apply_pres_k) then
          duplogic=(duplogic_1.or.duplogic_2).and.&
          data(ipres_loc,k) == data(ipres_loc,l) .and. &
          data(ier,k) < r1000 .and. data(ier,l) < r1000 .and. &
          muse(k) .and. muse(l)
       else
          duplogic=(duplogic_1.or.duplogic_2).and.&
          data(ier,k) < r1000 .and. data(ier,l) < r1000 .and. &
          muse(k) .and. muse(l)
       end if

       if (duplogic) then
          if(l_closeobs) then
             if(abs(data(itime,k)-hr_offset)<abs(data(itime,l)-hr_offset)) then
                 muse(l)=.false.
             else
                 muse(k)=.false.
                 exit lloop
             endif
          else
             tfact=min(one,abs(data(itime,k)-data(itime,l))/dfact1)
             dup(k)=dup(k)+one-tfact*tfact*(one-dfact)
             dup(l)=dup(l)+one-tfact*tfact*(one-dfact)
          endif
       end if
    end do lloop
  end do kloop

end subroutine dupcheck

end module dupcheckmod
