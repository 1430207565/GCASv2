program  post_mozart
 use module_global
 implicit none
 integer, parameter :: nx=360, ny=180,nt=8
 real flux(nx,ny,nt),bio(nx,ny,nt), ocn(nx,ny,nt), fossil(nx,ny,nt), fire(nx,ny,nt)
 real,allocatable, dimension(:,:,:) :: bio_m, ocn_m
 real,allocatable, dimension(:,:,:) :: fossil_m, fire_m
 real,allocatable, dimension(:,:,:) :: flux_m
 integer, allocatable, dimension(:) :: datetime
 real lat(ny), lon(nx)
 integer sd, ed, ir,i, j, k, t, d,m,r, y, im, id
 character(len=100) inputfile, outfile
 character(len=60) dir
 integer narg, tlen, ncycle
 CHARACTER(len=32) ::arg
 character(len=8) sday

narg=IARGC()
if(narg<2)then
 print*,'the input parameter error!'
 stop
endif
if(narg >0)then
call getarg(1,arg)
sday=arg
call getarg(2,arg)
read(arg,*)tlen
endif

call ensrf_init()

ncycle=tlen/rundays
tlen=rundays*ncycle

print*, sday, ncycle, tlen

dir=trim(opt_dir)//'/'//trim(casename)
! calculate total month and year

if(nrt==0)then
allocate(bio_m(nx,ny,ncycle))
allocate(ocn_m(nx,ny,ncycle))
allocate(fossil_m(nx,ny,ncycle))
allocate(fire_m(nx,ny,ncycle))
endif

allocate(flux_m(nx,ny,ncycle))
allocate(datetime(ncycle))

do i=1, ny
   lat(i)=-89.5+(i-1)
enddo
do i=1, nx
   lon(i)=0.5+(i-1)
enddo


!print*,'read_map end'

r=1
if(nrt==0)then
bio_m=0
ocn_m=0
fossil_m=0
fire_m=0

endif
flux_m=0

do im=1,ncycle

    read(sday,*) datetime(im)
    print*, datetime(im)

    do id=1, rundays

      
      inputfile=trim(dir)//'/emissions.opt.'//sday//'.nc'
      ! write(*,'(a,i4,2i2.2,a3)') 'emissions.opt.',y,m,d,'.nc'
      if(nrt==1)then
          call read_nc_nrt(inputfile,flux,nx,ny,nt)

           flux=flux/6.02e23*12*1.0e4*3600
     
      else
           call read_nc(inputfile,bio,ocn,fossil,fire,nx,ny,nt)
           ! molecules/cm2/s ==> gC/m2/hr
           bio=bio/6.02e23*12*1.0e4*3600 
           ocn=ocn/6.02e23*12*1.0e4*3600
           fossil=fossil/6.02e23*12*1.0e4*3600
           fire=fire/6.02e23*12*1.0e4*3600
           flux=bio+ocn+fossil+fire
      endif
      do t=1, nt
        if(nrt==0)then
           bio_m(:,:,im)=bio_m(:,:,im)+bio(:,:,t)*3
           ocn_m(:,:,im)=ocn_m(:,:,im)+ocn(:,:,t)*3
           fossil_m(:,:,im)=fossil_m(:,:,im)+fossil(:,:,t)*3
           fire_m(:,:,im)=fire_m(:,:,im)+fire(:,:,t)*3
        endif
        flux_m(:,:,im)=flux_m(:,:,im)+flux(:,:,t)*3
      enddo

      call nextday(sday)

    enddo

enddo

if(nrt==0)then
  bio_m=bio_m/rundays/24      !gC/m2/hr
  ocn_m=ocn_m/rundays/24
  fire_m=fire_m/rundays/24
  fossil_m=fossil_m/rundays/24
endif
flux_m=flux_m/rundays/24

if(nrt==0)then
call output_flux_grid(dir,casename,flux_m,bio_m,ocn_m,  &
                      fossil_m,fire_m,  &
                     nx,ny,ncycle,lat, lon, datetime,  opt_scheme)
endif

end program

subroutine handle_err1(status)
                integer status
                if (status .NE. NF_NOERR) then
                        write(*,*), NF_STRERROR(status), status
                        stop 'program stoped, and you must find the error in your source code.'
                endif
end

subroutine output_flux_grid(dir,casename,flux_m,bio_m,ocn_m, fossil_m,fire_m, nx,ny,nm,lat, lon, datetime, opt_scheme)
 include 'netcdf.inc'
 integer nx, ny, nm
 real lat(ny), lon(nx)
 real bio_m(nx,ny,nm), ocn_m(nx,ny,nm)   !gC/m2/month
 real fossil_m(nx,ny,nm), fire_m(nx,ny,nm)   !gC/m2/month
 real flux_m(nx,ny,nm) 
 integer datetime(nm)
 integer biomid, ocnmid, netmid, fossilmid,firemid, fid

 integer dims3d(3), start(3), count(3)
 integer dims2d(2), start2(2), count2(2)
 integer monthid, latid, lonid
 integer timevarid
 character(len=15) casename, var_fname,var_infile
 character(len=60) dir
 character(len=100) fname
 integer rcode, status, opt_scheme

 fname=trim(dir)//'/'//'posterior.fluxes.'//trim(casename)//'_weekly.nc'
 !------------------/----------------------------------------------------------
   ! create the file, clobbering any existing file, add global attributes
   !----------------------------------------------------------------------------
   rcode = nf_create(trim(fname),NF_CLOBBER,fid)
   if (rcode /= NF_NOERR) then
      write(6,*) "ERROR: creating file ",trim(fName)
      write(6,*) "ERROR: ",nf_strerror(rcode)
      stop
   end if

  !... begin to define your dimentions
        status=nf_def_dim(fid, 'weeks', nm, monthID)  ! the time is unlimited in this program
        if (status .NE. nf_noerr) call handle_err1(status)
  !       print*,'monthid',monthid,nm

        status=nf_def_dim(fid, 'lat', ny, LATID)
        if (status .NE. nf_noerr) call handle_err1(status)
  !      print*,'latid',latid,ny

        status=nf_def_dim(fid, 'lon', nx, LONID)
        if (status .NE. nf_noerr) call handle_err1(status)
  !      print*,'lonid',lonid,nx

  !     print*,'dimentions defining end'

        dims3d(3)=monthid
        dims3d(2)=latid
        dims3d(1)=lonid

     print*,'define vars'

     status=nf_def_var(fid,'netflux', nf_float, 3, dims3d, netmid)
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_put_att_text(fid, netmid, 'long_name', 24, 'optmized weekly net flux')
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_put_att_text(fid, netmid, 'units', 8, 'gC/m2/hr')
     if (status .NE. nf_noerr) call handle_err1(status)

     status=nf_def_var(fid,'fossil_imp', nf_float, 3, dims3d, fossilmid)
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_put_att_text(fid, fossilmid, 'long_name', 33, 'prescribed weekly fossil emission')
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_put_att_text(fid, fossilmid, 'units', 8, 'gC/m2/hr')
     if (status .NE. nf_noerr) call handle_err1(status)

     status=nf_def_var(fid,'fire_imp', nf_float, 3, dims3d, firemid)
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_put_att_text(fid, firemid, 'long_name', 31, 'prescribed weekly fire emission')
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_put_att_text(fid, firemid, 'units', 8, 'gC/m2/hr')
     if (status .NE. nf_noerr) call handle_err1(status)



     status=nf_def_var(fid,'bio_opt', nf_float, 3, dims3d, biomid)
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_put_att_text(fid, biomid, 'long_name', 24, 'optmized weekly bio flux')
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_put_att_text(fid, biomid, 'units', 8, 'gC/m2/hr')
     if (status .NE. nf_noerr) call handle_err1(status)

     status=nf_def_var(fid,'ocn_opt', nf_float, 3, dims3d, ocnmid)
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_put_att_text(fid, ocnmid, 'long_name', 24, 'optmized weekly ocn flux')
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_put_att_text(fid, ocnmid, 'units', 8, 'gC/m2/hr')
     if (status .NE. nf_noerr) call handle_err1(status)

     status=nf_def_var(fid, 'lat', nf_float, 1, latid, latvarid)
     if (status .NE. nf_noerr) call handle_err1(status)

     status=nf_def_var(fid, 'lon', nf_float, 1, lonid, lonvarid)
     if (status .NE. nf_noerr) call handle_err1(status)
     
     status=nf_def_var(fid, 'datetime', nf_int, 1, monthid, timevarid)
     if (status .NE. nf_noerr) call handle_err1(status)


      status=nf_enddef(fid)

   count(3)=nm
   count(2)=ny
   count(1)=nx
   start=1

   print*,'start writing'

   status=nf_put_vara_real(fid, netmid, start, count, flux_m)

   status=nf_put_vara_real(fid, biomid, start, count, bio_m)
   status=nf_put_vara_real(fid, ocnmid, start, count, ocn_m)
   status=nf_put_vara_real(fid, fossilmid, start, count, fossil_m)
   status=nf_put_vara_real(fid, firemid, start, count, fire_m)

   status=nf_put_vara_real(fid, latvarid, 1, ny, lat)
   status=nf_put_vara_real(fid, lonvarid, 1, nx, lon)
   status=nf_put_vara_int(fid, timevarid, 1, nm, datetime)
   status=nf_close(fid)

end subroutine
 
