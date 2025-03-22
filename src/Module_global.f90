module module_global

  implicit none
  integer, parameter  :: ibio=1, iocn=2, ifire=3, ifossil=4
  integer, public :: bnx, bny          ! background flux domains
  real, public, allocatable, DIMENSION(:) :: blat, blon  !background flux
  integer, public :: bnvar
  integer, public :: opt_scheme   !1, bio only, 2, bio and ocn, 3, bio, ocn, and fire
  integer, public :: netflux      !1, optimize netflux on land, 0, no
  real, public :: bunc(3)  !1, bio, 2, ocn, 3, fire
  integer, public :: ensemble, nrt, issrf, issat, isgosat, isoco2, istansat
  character(len=20), public :: casename
  LOGICAL, public :: DEBUG
  integer, public :: nlon, nlat, nlev   ! nlon, nlat: transport model; nx, ny: carbon flux
  real, public, allocatable, DIMENSION(:) :: xlat, xlon   !tranport model
  character(len=100), public :: mozartdir, obsdir_srf, obsdir_gosat, obsdir_oco2, obsdir_tansat, opt_dir, priordir
  integer, public :: resb, resm  ! resb, background, resm, resoluation of transport model
  integer, public :: year, month, day, icdate, rundays, da_step
  integer, public :: erropt
  real, public    :: alpha(8), tune, ermax(8), ermin(8), L_ratio(2), L_scale(2)

  LOGICAL                                             ::    alive

  integer         , PARAMETER                         :: logun = 1111
  REAL            ,PARAMETER                          :: MISSING = -999.
  INTEGER         , PARAMETER                         :: maxobs = 10000


  public :: dmin2

  character(len=256) :: logfile
  data logfile/'ensrf.log'/

  !*******************************Netcdf
  ! NETCDF
  !*******************************Netcdf
  INTEGER                                         ::       nxID,nyID,nspID,enID,area_pointID,NumOfObsID
  INTEGER                                         ::       ncid,retval,VARID, dimid
  INTEGER                                         ::       DIMS2D(2),DIMS3D(3),DIMS4D(4)


contains

  integer function dmin2(a,b,n)
        implicit none
        integer n
        real b(n),a
        integer i
        real m,m1

        m=(b(1)-a)*(b(1)-a)
        dmin2=1
        do i=2, n
           m1=(b(i)-a)*(b(i)-a)
           if(m1<m)then
                 m=m1
                dmin2=i
           endif
        enddo

  end function

 subroutine nextday(sday)
  implicit none
  character(len=8) sday
  integer year, month, day
  integer mon1(12), mon2(12), mon(12)

  data mon1/31,28,31,30,31,30,31,31,30,31,30,31/
  data mon2/31,29,31,30,31,30,31,31,30,31,30,31/

  read(sday(1:4),*) year
  read(sday(5:6),*) month
  read(sday(7:8),*) day
  if(year==2000.or.year==2004.or.year==2008.or.year==2012.or.year==2016.or. &
     year==2020.or.year==2024.or.year==2028.or.year==2032.or.year==2036.or. &
     year==2040.or.year==2044.or.year==2048)then
     mon=mon2
  else
     mon=mon1
  endif
  day=day+1
  if(day>mon(month))then
      day=1
      month=month+1
      if(month>12)then
         month=1
         year=year+1
      endif
  endif

  write(sday,'(i4,i2.2,i2.2)') year, month, day

end subroutine


end module
