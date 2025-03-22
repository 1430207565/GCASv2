program mztout2initconc
 implicit none
 integer nx, ny, nz, en_size
 real, allocatable, dimension(:) :: lat,lon, hyam, hybm
 real, allocatable, dimension(:,:) :: ps
 real, allocatable, dimension(:,:,:) :: co2
 real p0
 integer date, datesec
 character(len=100) file_in, file_out, arg
 integer i, j, k, narg

 narg=IARGC()
  if(narg<3)then
    print*,'the input parameter error!'
   stop
  endif
  if(narg >0)then
    call getarg(1,arg)
    read(arg,*) file_in
    call getarg(2,arg)
    read(arg,*) file_out
    call getarg(3,arg)
    read(arg,*) en_size
  endif

  call readinfo(file_in, nx, ny, nz)

  print*, 'input file:', trim(file_in)
  print*, 'nx,ny,nz:', nx,ny,nz

  allocate(co2(nx,ny,nz))
  allocate(ps(nx,ny))
  allocate(hyam(nz))
  allocate(hybm(nz))
  allocate(lon(nx))
  allocate(lat(ny))
  
  call read_mzt(file_in, co2, ps, hyam, hybm, P0, date, datesec, lat, lon, nx, ny, nz)
 
  print*, 'date, datesec:', date, datesec
  print*,'output file:', trim(file_out)

  call write_mzt_ic(file_out,co2,ps, hyam, hybm, P0, date, datesec, lat, lon, nx,ny,nz, en_size)

end program


subroutine readinfo(fname, nx,ny,nz)
   include 'netcdf.inc'
   integer nx,ny,nz
   integer ncid
   character(len=100) fname

   nx=80
   ny=60

   print*,trim(fname)

   status=nf_open(trim(fname),nf_nowrite,ncid)
   if ( status/=nf_noerr ) then
      write (*,*) nf_strerror(status)
      write (*,*) '!!Open NetCDF file Error!!'
       stop
   end if

   status=nf_inq_dimid(ncid, 'lon', dimid)
      if ( status/=nf_noerr ) write (*,*) nf_strerror(status)

          status=nf_inq_dimlen(ncid, dimid, nx)
      if ( status/=nf_noerr ) write (*,*) nf_strerror(status)

      status=nf_inq_dimid(ncid, 'lat', dimid)
      if ( status/=nf_noerr ) write (*,*) nf_strerror(status)

      status=nf_inq_dimlen(ncid, dimid, ny)
      if ( status/=nf_noerr ) write (*,*) nf_strerror(status)

      status=nf_inq_dimid(ncid, 'lev', dimid)
      if ( status/=nf_noerr ) write (*,*) nf_strerror(status)

      status=nf_inq_dimlen(ncid, dimid, nz)
      if ( status/=nf_noerr ) write (*,*) nf_strerror(status)

   status = nf_close(ncid)

 endsubroutine
!write_mzt_ic(file_out,co2,ps, hyam, hybm, P0, date, datesec, lat, lon,
!nx,ny,nz, en_size)
subroutine write_mzt_ic(fname, co2, ps, hyam, hybm, P0, date, datesec, lat, lon, nx, ny, nz, en_size)
  include 'netcdf.inc'
  integer nx,ny,nz, en_size
  real ps(nx,ny)
  real co2(nx,ny,nz)
  real hyam(nz), hybm(nz)
  character(len=100) fname
  integer fid, levid, latid, lonid
  integer dims(3), start(3), count(3)
  integer vco2id, vENid(en_size), vPSid, vP0id, vdateid, vdatesecid
  integer vlatid, vlonid, vlevid, vhyamid, vhybmid
  integer date, datesec, i
  real P0
  real lat(ny), lon(nx), lev(nz)
  character(len=10) vname

  do i=1, nz
    lev(i)=i
  enddo

  rcode = nf_create(trim(fname),NF_CLOBBER,fid)
  if (rcode /= NF_NOERR) then
      write(6,*) "ERROR: creating file ",trim(fName)
      write(6,*) "ERROR: ",nf_strerror(rcode)
      stop
  end if
 
  !... begin to define your dimentions
  status=nf_def_dim(fid, 'lev', nz, levid)  
  status=nf_def_dim(fid, 'lat', ny, latID)
  status=nf_def_dim(fid, 'lon', nx, lonID)

  dims(3)=levid
  dims(2)=latid
  dims(1)=lonid

  status=nf_def_var(fid, 'CO2', nf_float, 3, dims, vco2id)
  status=nf_def_var(fid, 'PS', nf_float, 2, dims(1:2), vPSid)
  status=nf_def_var(fid, 'lat', nf_float, 1, latid, vlatid)
  status=nf_def_var(fid, 'lon', nf_float, 1, lonid, vlonid)
  status=nf_def_var(fid, 'lev', nf_float, 1, levid, vlevid)
  status=nf_def_var(fid, 'hyam', nf_float, 1, levid, vhyamid)
  status=nf_def_var(fid, 'hybm', nf_float, 1, levid, vhybmid)
  status=nf_def_var(fid, 'date', nf_int, 0, 0, vdateid)
  status=nf_def_var(fid, 'datesec', nf_int, 0, 0, vdatesecid)
  status=nf_def_var(fid, 'P0', nf_float, 0, 0, vP0id)
  
  do i=1, en_size
    write(vname,'(a,i3.3)') 'EN', i
    status=nf_def_var(fid, vname, nf_float, 3, dims, vENid(i))
  enddo
   
  status=nf_enddef(fid)
  
  count(3)=nz
  count(2)=ny
  count(1)=nx
  start=1
  
  status=nf_put_vara_real(fid, vPSid, start(1:2), count(1:2), ps)
  status=nf_put_vara_real(fid, vco2id, start, count, co2)
  do i=1, en_size
     status=nf_put_vara_real(fid, vENid(i), start, count, co2)
  enddo

  datesec=0

  status=nf_put_vara_real(fid, vlatid, 1, ny, lat)
  status=nf_put_vara_real(fid, vlonid, 1, nx, lon)
  status=nf_put_vara_real(fid, vlevid, 1, nz, lev)
  status=nf_put_vara_real(fid, vhyamid, 1, nz,hyam)
  status=nf_put_vara_real(fid, vhybmid, 1, nz,hybm)
  status=nf_put_vara_real(fid, vP0id, 1, 1, P0)
  status=nf_put_vara_int(fid, vdateid, 1, 1, date)
  status=nf_put_vara_int(fid, vdatesecid, 1, 1, datesec)

  status=nf_close(fid)

end subroutine

!read_mzt(file_in, co2, ps, hyam, hybm, P0, date, datesec, lat, lon, nx, ny, nz)
subroutine read_mzt(file_in, co2,ps, hyam, hybm, P0, date, datesec, lat, lon, nx, ny, nz)
 include 'netcdf.inc'
 integer nx, ny, nz
 real co2(nx,ny,nz), ps(nx, ny)
 real hyam(nz), hybm(nz), lat(ny), lon(nx)
 real P0
 integer date, datesec
 character(len=100) file_in
 integer i, j
 integer start(4), count(4)

 status=nf_open(trim(file_in),nf_nowrite,ncid)
   if ( status/=nf_noerr ) then
      write (*,*) nf_strerror(status)
      write (*,*) '!!Open NetCDF file Error!!'
       stop
    end if

 status= nf_inq_varid(ncid, 'EN001_VMR_avrg',varid)
 start=1
 count(1)=nx
 count(2)=ny
 count(3)=nz
 count(4)=1
 status=nf_get_vara_real(ncid,  varid,start,count,co2)

 status= nf_inq_varid(ncid, 'PS',varid)
 count(3)=1
 status=nf_get_vara_real(ncid,  varid,start(1:3),count(1:3),ps)

 status= nf_inq_varid(ncid, 'lat',varid)
 count(1)=ny
 status=nf_get_vara_real(ncid,  varid,start(1),count(1),lat)
 status= nf_inq_varid(ncid, 'lon',varid)
 count(1)=nx
 status=nf_get_vara_real(ncid,  varid,start(1),count(1),lon)

 status= nf_inq_varid(ncid, 'hyam',varid)
 count(1)=nz
 status=nf_get_vara_real(ncid,  varid,start(1),count(1),hyam)

 status= nf_inq_varid(ncid, 'hybm',varid)
 count(1)=nz
 status=nf_get_vara_real(ncid,  varid,start(1),count(1),hybm)


 status= nf_inq_varid(ncid, 'date',varid)
 count(1)=1
 status=nf_get_vara_int(ncid,  varid,start(1),count(1),date)

 status= nf_inq_varid(ncid, 'datesec',varid)
 count(1)=1
 status=nf_get_vara_int(ncid,  varid,start(1),count(1),datesec)

 status= nf_inq_varid(ncid, 'P0',varid)
 status=nf_get_var_real(ncid,  varid,p0)

 status = nf_close(ncid)

end subroutine


