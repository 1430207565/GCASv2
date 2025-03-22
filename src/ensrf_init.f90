subroutine ensrf_init()
 use module_global
! ************************************************************

  NAMELIST    /time_control/ year, month, day, rundays, da_step
  NAMELIST    /assimilation/ casename, ensemble, nrt, opt_scheme, netflux, opt_dir, L_scale, L_ratio, debug
  NAMELIST    /transport_model/ resm, mozartdir
  NAMELIST    /background_flux/ resb, priordir, bunc
  NAMELIST    /observations/ issrf, obsdir_srf, issat, isgosat, obsdir_gosat, isoco2, obsdir_oco2, istansat, obsdir_tansat

  OPEN(unit=logun,file=logfile,access='append', form='formatted')
  OPEN(unit=1,file="ensrf.nml",status='old',form='formatted',iostat=ios)
  IF ( ios /= 0 ) STOP "ERROR opening ensrf.nml"
  READ(1,nml=time_control)
  READ(1,nml=assimilation)
  READ(1,nml=transport_model)
  READ(1,nml=background_flux)
  READ(1,nml=observations)
! 
  if(issrf==0.and.issat==0)then
     print*,'***No Observation is set! please check the namelist!***'
     stop
  endif

  if(issrf==0.and.issat==1)then
     if(isgosat+isoco2+istansat==0)then
        print*,'***XCO2 is selected, but no product is set! please check the namelist!***'
        stop
     endif
  endif

  if(opt_scheme==1)then
     bnvar=1
  elseif(opt_scheme==2)then
     bnvar=2
  elseif(opt_scheme==3)then
     bnvar=3
  endif
  icdate=year*10000+month*100+day
  !print*, icdate, year, month, day
  
  if(resm==3)then
    nlon=128   !transport model resolusion
    nlat=64
    nlev=28
  elseif(resm==2)then   !use GEOS-5 data
    nlon=144
    nlat=96
    nlev=56
  else
    nlon=288
    nlat=192
    nlev=28
  endif
  if(resb>6) then
    print*, "the resoluton of background flux is error, it should be <=6!, please check!"
    stop
  endif
  bny=180/resb    !resb<=6
  bnx=360/resb

  allocate (xlat(nlat))
  allocate (xlon(nlon))
  allocate (blat(bny))
  allocate (blon(bnx))
  do i=1, nlat
    xlat(i)=(i-1)*180.0/(nlat-1)-90
  enddo
  do i=1, nlon
    xlon(i)=(i-1)*360.0/nlon
  enddo
  do i=1, bny
    blat(i)=(i-1)*resb+resb/2.0-90
  enddo
  do i=1, bnx
    blon(i)=(i-1)*resb+resb/2.0
  enddo
  if(nrt==1)then
    bnvar=1
  endif


end subroutine

