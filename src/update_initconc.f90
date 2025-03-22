program update_conc
  use module_global
  implicit none
  real*8, allocatable, dimension(:,:,:) :: co2
  integer  duration, en_size, ios
  character(len=8) sday
  character(len=25) varname
  character(len=100) filename, filename1
  integer i, j, ite, its, it, imonth, iday, sdate

  call ensrf_init()

   allocate (co2(nlon,nlat,nlev))
!  sdate=icdate
!  read(sday(5:6),*) imonth
  write(sday,'(i8)') icdate
  do iday = 1,rundays
  call nextday(sday)
!  read(sday,*) sdate
  enddo

!print*, sday

  filename=trim(mozartdir)//'/run_forward/rest/'//trim(casename)//'.mz4.r.'//sday(1:4)//'-'//sday(5:6)//'-'//sday(7:8)//'-00000.nc'
  filename1=trim(mozartdir)//'/run_ensemble/rest/'//trim(casename)//'.mz4.r.'//sday(1:4)//'-'//sday(5:6)//'-'//sday(7:8)//'-00000.nc'
  
  varname='EN001'
  call readnc3d(filename,varname,co2,nlon, nlat, nlev)

  call writenc3d(filename1,co2,nlon,nlat, nlev, ensemble)

end program

!                      time, nlon,nlat,duration*8
subroutine writenc3d(filename,co2,nx,ny,nz,en_size)
 include 'netcdf.inc'
 integer nx, ny, nz, en_size
 real*8 co2(nx,ny,nz)
 integer varid
 integer stepid, latid, lonid
 character(len=100) filename
 character(len=25) vname
 integer rcode, status
 integer i, j, ncid, s
  integer*4   :: start(10)
  integer*4   :: count(10)
  integer     :: dimids(10)! allow up to 10 dimensions
  integer     :: dimid, xtype, ndim
  character(len=31) :: dummy

  status=nf_open(trim(filename),nf_write,ncid)
   if ( status/=nf_noerr ) then
      write (*,*) nf_strerror(status)
      write (*,*) '!!Open NetCDF file Error!!'
       stop
    end if
 do s=1, en_size 
  write(vname,'(a2,i3.3)') 'EN', s
    status= nf_inq_varid(ncid, trim(vname),varid)
    status=nf_inq_var(ncid,  varid,dummy,xtype,ndim,dimids,natts)
           if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      do j=1,ndim
      status=nf_inq_dim(ncid,dimids(j),dummy,len)
           if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      start(j)=1 ; count(j)=len
      end do
    status=nf_put_vara_double(ncid,  varid,start,count,co2)
 enddo

 status = nf_close(ncid)


end subroutine

subroutine readnc3d(filename,varname,a,nx, ny, nt)
  include 'netcdf.inc'
  integer nx, ny, nt
  real*8 a(nx,ny,nt)
  character(len=25) varname
  character(len=100) filename
 
  status=nf_open(trim(filename),nf_nowrite,ncid)
   if ( status/=nf_noerr ) then
      write (*,*) nf_strerror(status)
      write (*,*) '!!Open NetCDF file Error!!'
       stop
    end if
  call get_3d_var(ncid,varname,a,nx,ny,nt)
  status = nf_close(ncid)

end subroutine

subroutine get_3d_var(ncid,vname,a,nx,ny,nt)
  include 'netcdf.inc'
  integer nx,ny,ncid,nt
  real*8 a(nx,ny,nt)
  character(len=25) vname
  integer*4   :: start(10)
  integer*4   :: count(10)
  integer     :: dimids(10)! allow up to 10 dimensions
  integer     :: dimid, xtype, ndim
  character(len=31) :: dummy

   !   Retrieve data for Variable 'XLAT'
      status= nf_inq_varid(ncid, trim(vname),varid)
      status=nf_inq_var(ncid,  varid,dummy,xtype,ndim,dimids,natts)
           if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      do j=1,ndim
      status=nf_inq_dim(ncid,dimids(j),dummy,len)
           if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      start(j)=1 ; count(j)=len
      end do
      status=nf_get_vara_double(ncid,  varid,start,count,a)

  end subroutine
