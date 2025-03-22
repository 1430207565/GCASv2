subroutine check_prior(yy,mm,dd,sday,priordir)
   implicit none
   integer yy, mm, dd
   character(len=100) priordir, FILE_NAME1,FILE_NAME2,FILE_NAME3,FILE_NAME4
   character(len=8) sday
   logical alive1, alive2, alive3, alive4

   write(sday,'(i4,i2.2,i2.2)') yy, mm, dd
   FILE_NAME1=trim(priordir)//'/bio/emissions.bio.'//sday//'.nc'
   FILE_NAME2=trim(priordir)//'/ocn/emissions.ocn.'//sday//'.nc'
   FILE_NAME3=trim(priordir)//'/fire/emissions.fire.'//sday//'.nc'
   FILE_NAME4=trim(priordir)//'/fossil/emissions.fossil.'//sday//'.nc'
   INQUIRE(FILE=trim(FILE_NAME1),EXIST=ALIVE1)
   INQUIRE(FILE=trim(FILE_NAME2),EXIST=ALIVE2)
   INQUIRE(FILE=trim(FILE_NAME3),EXIST=ALIVE3)
   INQUIRE(FILE=trim(FILE_NAME4),EXIST=ALIVE4)
   do while((.not.alive1).or.(.not.alive2).or.(.not.alive3).or.(.not.alive4))
     yy=yy-1
     if(yy<2000)then
         stop
     endif
     write(sday,'(i4,i2.2,i2.2)') yy, mm, dd
     FILE_NAME1=trim(priordir)//'/bio/emissions.bio.'//sday//'.nc'
     FILE_NAME2=trim(priordir)//'/ocn/emissions.ocn.'//sday//'.nc'
     FILE_NAME3=trim(priordir)//'/fire/emissions.fire.'//sday//'.nc'
     FILE_NAME4=trim(priordir)//'/fossil/emissions.fossil.'//sday//'.nc'
     INQUIRE(FILE=trim(FILE_NAME1),EXIST=ALIVE1)
     INQUIRE(FILE=trim(FILE_NAME2),EXIST=ALIVE2)
     INQUIRE(FILE=trim(FILE_NAME3),EXIST=ALIVE3)
     INQUIRE(FILE=trim(FILE_NAME4),EXIST=ALIVE4)
   enddo
end subroutine

!call readavebio(mbio(:,:,1:8),nx,ny,8)
subroutine readavebio(bio,nx,ny,nt)
      implicit none
      include 'netcdf.inc'
      integer nx,ny,nt, ncid, varid, retval
      real bio(nx,ny,nt)
      character(len=100) FILE_NAME
      character(len=50) code
      code='subroutine readavebio'

      FILE_NAME='../input/emissions.bio.multiyearmean.nc'
  !    print*, File_name

      retval = nf_open(FILE_NAME, NF_NOWRITE, ncid)
      if (retval .ne. nf_noerr) call handle_err(retval,code)

      retval = nf_inq_varid    ( ncid, "emissions", varid )
      retval = nf_get_var_real ( ncid, varid, bio )
      retval = NF_CLOSE(NCID)

end subroutine 

subroutine readprior(priordir,bio,ocn,fossil,fire,nx,ny,nt,sday,netf)
      implicit none
      include 'netcdf.inc'
      integer nx,ny,nt, ncid, varid
      character(len=8) sday, sday2
      integer yy, mm, dd, netf
      real bio(nx,ny,nt),ocn(nx,ny,nt), fossil(nx,ny,nt), fire(nx,ny,nt)
      character(len=100) priordir, FILE_NAME
      integer retval
      character(len=50) code
      code='subroutine readprior'
     
      read(sday(1:4),*) yy
      read(sday(5:6),*) mm
      read(sday(7:8),*) dd
      
      call check_prior(yy,mm,dd,sday2,priordir)
      if(sday2/=sday)then
         print*, 'Current prior fluxes not existed, use '//sday2//' data as priori'
      endif
      FILE_NAME=trim(priordir)//'/bio/emissions.bio.'//sday2//'.nc'
  !    print*, File_name
     
      retval = nf_open(FILE_NAME, NF_NOWRITE, ncid)
      if (retval .ne. nf_noerr) call handle_err(retval,code)

      retval = nf_inq_varid    ( ncid, "emissions", varid )
      retval = nf_get_var_real ( ncid, varid, bio )
      retval = NF_CLOSE(NCID)
      IF (retval .NE. NF_NOERR) CALL HANDLE_ERR(retval,code)


      FILE_NAME=trim(priordir)//'/ocn/emissions.ocn.'//sday2//'.nc'
      retval = nf_open(FILE_NAME, NF_NOWRITE, ncid)
      if (retval .ne. nf_noerr) call handle_err(retval,code)

      retval = nf_inq_varid    ( ncid, "emissions", varid )
      retval = nf_get_var_real ( ncid, varid, ocn )
      retval = NF_CLOSE(NCID)
      IF (retval .NE. NF_NOERR) CALL HANDLE_ERR(retval,code)


      FILE_NAME=trim(priordir)//'/fossil/emissions.fossil.'//sday2//'.nc'
      retval = nf_open(FILE_NAME, NF_NOWRITE, ncid)
      if (retval .ne. nf_noerr) call handle_err(retval,code)

      retval = nf_inq_varid    ( ncid, "emissions", varid )
      retval = nf_get_var_real ( ncid, varid,fossil )
      retval = NF_CLOSE(NCID)
      IF (retval .NE. NF_NOERR) CALL HANDLE_ERR(retval, code)

      FILE_NAME=trim(priordir)//'/fire/emissions.fire.'//sday2//'.nc'
      retval = nf_open(FILE_NAME, NF_NOWRITE, ncid)
      if (retval .ne. nf_noerr) call handle_err(retval, code)

      retval = nf_inq_varid    ( ncid, "emissions", varid )
      retval = nf_get_var_real ( ncid, varid, fire )
      retval = NF_CLOSE(NCID)
      IF (retval .NE. NF_NOERR) CALL HANDLE_ERR(retval)
     ! bio=bio+fire
     ! fire=0
      if(netf==1)then
           bio=bio+fire+fossil
           fossil=0
           fire=0
!      else
!           bio=bio+fire
!           fire=0
      endif
end subroutine

subroutine read_nc(fname,bio,ocn,fossil,fire,nx,ny,nt)

include 'netcdf.inc'
  integer nx, ny,nt
  real bio(nx,ny,nt), ocn(nx,ny,nt)   !, zm(nx,ny,nz,nt)   !, zi(nx,ny,nz,nt)
  real fossil(nx,ny,nt), fire(nx,ny,nt)
  integer i, j, ncid
  character(len=100) fname
  character(len=25) vname,var_infile


   status=nf_open(trim(fname),nf_nowrite,ncid)
   if ( status/=nf_noerr ) then
      write (*,*) nf_strerror(status)
      write (*,*) '!!Open NetCDF file Error!!'
       stop
    end if
  call get_3d_var(ncid,'bio',bio,nx,ny,nt)
  call get_3d_var(ncid,'ocn',ocn,nx,ny,nt)
  call get_3d_var(ncid,'fossil',fossil,nx,ny,nt)
  call get_3d_var(ncid,'fire',fire,nx,ny,nt)
 status = nf_close(ncid)

end subroutine

subroutine read_nc_nrt(fname,bio,nx,ny,nt)

include 'netcdf.inc'
  integer nx, ny,nt
  real bio(nx,ny,nt)
  integer i, j, ncid
  character(len=100) fname
  character(len=25) vname,var_infile


   status=nf_open(trim(fname),nf_nowrite,ncid)
   if ( status/=nf_noerr ) then
      write (*,*) nf_strerror(status)
      write (*,*) '!!Open NetCDF file Error!!'
       stop
    end if
  call get_3d_var(ncid,'netflux',bio,nx,ny,nt)
 status = nf_close(ncid)

end subroutine

!createlambda(filename, blat,blon, bnx, bny, bnvar)
subroutine createlambda(filename,lat, lon, nx,ny,nv)
  implicit none
  include 'netcdf.inc'
  integer nx, ny, nv
  real lat(ny), lon(nx)
  integer rcode
  character(len=100) filename
  integer fid, status, latid, lonid, varid, stepid, latvarid, lonvarid, datevarid, lambdaid
  integer nobsid, Lscaleid, xco2id, nxco2id, dfluxid, co2id, nco2id
  integer dims(4), dims3(3)

  rcode = nf_create(trim(filename),NF_CLOBBER,fid)
  if (rcode /= NF_NOERR) then
       write(6,*) "ERROR: creating file ",trim(filename)
       write(6,*) "ERROR: ",nf_strerror(rcode)
       stop
  end if
  !... begin to define your dimentions
  status=nf_def_dim(fid, 'time', nf_unlimited, STEPID)  ! the time is unlimited in this program
  status=nf_def_dim(fid, 'lat', ny, LATID)
  status=nf_def_dim(fid, 'lon', nx, LONID)
  status=nf_def_dim(fid, 'var', nv, VARID)

  dims(4)=stepid
  dims(3)=varid
  dims(2)=latid
  dims(1)=lonid

  dims3(3)=stepid
  dims3(2)=latid
  dims3(1)=lonid

  status=nf_def_var(fid, 'lat', nf_float, 1, latid, latvarid)
  status=nf_def_var(fid, 'lon', nf_float, 1, lonid, lonvarid)
  status=nf_def_var(fid, 'date', nf_int, 1, stepid, datevarid)
  status=nf_def_var(fid, 'lambda', nf_float, 4, dims, lambdaid)
  status=nf_def_var(fid, 'Nobs', nf_float, 4, dims, nobsid)
  status=nf_def_var(fid, 'Lscale', nf_float, 4, dims, Lscaleid)
  status=nf_def_var(fid, 'xbias', nf_float, 3, dims3, xco2id)
  status=nf_def_var(fid, 'xnum', nf_float, 3, dims3, nxco2id)
  status=nf_def_var(fid, 'xbias_org', nf_float, 3, dims3, xco2id)
  status=nf_def_var(fid, 'xbias_opt', nf_float, 3, dims3, nxco2id)
  status=nf_def_var(fid, 'sbias_org', nf_float, 3, dims3, co2id)
  status=nf_def_var(fid, 'sbias_opt', nf_float, 3, dims3, nco2id)
  status=nf_def_var(fid, 'dflux', nf_float, 4, dims, dfluxid)
  status=nf_enddef(fid)
  status=nf_put_vara_real(fid, latvarid, 1, ny, lat)
  status=nf_put_vara_real(fid, lonvarid, 1, nx, lon)
  status=nf_close(fid)
end subroutine

!call writelambda(filename,lambda,bnx,bny,bnvar,sday, da_step)
subroutine writelambda(filename,lambda,dflux,nx,ny,nv,sday,step)
 implicit none
 include 'netcdf.inc'
 integer nx, ny, nv, step, datetime
 real lambda(nx,ny,nv), dflux(nx,ny,nv)
 character(len=100) filename
 character(len=50) code
 character(len=8) sday
 integer start(4), count(4), ncid, varid, retval
 
 read(sday, *) datetime

 code='subroutine read lambda'
 

 retval = nf_open(trim(filename), NF_WRITE, ncid)
  if (retval .ne. nf_noerr) call handle_err(retval,code)
  count(4)=1
  count(3)=nv
  count(2)=ny
  count(1)=nx
  start(4)=step
  start(1:3)=1

  retval = nf_inq_varid    ( ncid, "lambda", varid )
  retval = nf_put_vara_real ( ncid, varid, start, count, lambda )

  retval = nf_inq_varid    ( ncid, "dflux", varid )
  retval = nf_put_vara_real ( ncid, varid, start, count, dflux )

  retval = nf_inq_varid    ( ncid, "date", varid )
  retval = nf_put_vara_int ( ncid, varid, step, 1, datetime )
 
  retval = NF_CLOSE(NCID) 
 
 
end subroutine

!call writelambda(filename,lambda,bnx,bny,bnvar,sday, da_step)
subroutine writenobs(filename,nobs,lscale,nx,ny,nv,step)
 implicit none
 include 'netcdf.inc'
 integer nx, ny, nv, step
 integer nobs(nx,ny,nv)
 real lscale(nx,ny,nv)
 character(len=100) filename
 character(len=50) code
 integer start(4), count(4), ncid, varid, retval
 integer i, j, k

 code='subroutine write nobs'

 retval = nf_open(trim(filename), NF_WRITE, ncid)
  if (retval .ne. nf_noerr) call handle_err(retval,code)
  count(4)=1
  count(3)=nv
  count(2)=ny
  count(1)=nx
  start(4)=step
  start(1:3)=1

  retval = nf_inq_varid    ( ncid, "Nobs", varid )
  retval = nf_put_vara_real ( ncid, varid, start, count, real(nobs) )

  retval = nf_inq_varid    ( ncid, "Lscale", varid )
  retval = nf_put_vara_real ( ncid, varid, start, count, lscale )

  retval = NF_CLOSE(NCID)


end subroutine


!call writelambda(filename,lambda,bnx,bny,bnvar,sday, da_step)
subroutine writebxco2(filename,xco2,nxco2,nx,ny,step)
 implicit none
 include 'netcdf.inc'
 integer nx, ny, step
 real xco2(nx,ny)
 real nxco2(nx,ny)
 character(len=100) filename
 character(len=50) code
 integer start(3), count(3), ncid, varid, retval
 integer i, j, k

 code='subroutine write nobs'

 retval = nf_open(trim(filename), NF_WRITE, ncid)
  if (retval .ne. nf_noerr) call handle_err(retval,code)
  count(3)=1
  count(2)=ny
  count(1)=nx
  start(3)=step
  start(1:2)=1

  retval = nf_inq_varid    ( ncid, "xbias", varid )
  retval = nf_put_vara_real ( ncid, varid, start, count, xco2 )

  retval = nf_inq_varid    ( ncid, "xnum", varid )
  retval = nf_put_vara_real ( ncid, varid, start, count, nxco2 )

  retval = NF_CLOSE(NCID)


end subroutine


subroutine read_ocean_map(b, nx, ny)
 include 'netcdf.inc'
 integer nx, ny, nr
 integer b(nx,ny)
 integer i, j, ncid, varid
 integer start(2), count(2)
 character(len=30) fname
 character(len=25) vname

 fname='../input/region_gcas.nc'
! print*, fname
 status=nf_open(trim(fname),nf_nowrite,ncid)
   if ( status/=nf_noerr ) then
      write (*,*) nf_strerror(status)
      write (*,*) '!!Open NetCDF file Error!!'
       stop
    end if
 vname='region'
 status= nf_inq_varid(ncid, trim(vname),varid)
 start=1
 count(1)=nx
 count(2)=ny
 status=nf_get_vara_int(ncid,  varid,start,count,b) 

 do i=1, nx
   do j=1, ny
     if(b(i,j)>=52) b(i,j)=0
   enddo
 enddo
 
 status = nf_close(ncid)

end subroutine

subroutine get_3d_var(ncid,vname,a,nx,ny,nt)
  include 'netcdf.inc'
  integer nx,ny,ncid,nt
  real a(nx,ny,nt)
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
      status=nf_get_vara_real(ncid,  varid,start,count,a)

  end subroutine

subroutine get_2d_var(ncid,vname,a,nx,ny)
  include 'netcdf.inc'
  integer nx,ny,ncid
  real a(nx,ny)
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
      status=nf_get_vara_real(ncid,  varid,start,count,a)

  end subroutine


!readlambda(filename, lambda_old, bnx, bny, bnvar, da_step-1)
subroutine readlambda(filename, lambda, nx, ny, nv, step)
  implicit none
  include 'netcdf.inc'
  integer nx, ny, nv, step
  real lambda(nx,ny,nv)
  character(len=100) filename
  character(len=50) code
  integer ncid, retval, varid
  integer start(4), count(4)
  code='subroutine read lambda'

!  print*, nx, ny, nv, step
!  print*, trim(filename)

  retval = nf_open(trim(filename), NF_NOWRITE, ncid)
  if (retval .ne. nf_noerr) call handle_err(retval,code)

  count(4)=1
  count(3)=nv
  count(2)=ny
  count(1)=nx
  start(4)=step
  start(1:3)=1

  retval = nf_inq_varid    ( ncid, "lambda", varid )
!  print*, varid
  retval = nf_get_vara_real ( ncid, varid, start, count, lambda )
!  print*, lambda
!  print*,'------------------------------------------'
  retval = NF_CLOSE(NCID)


end subroutine


subroutine writeopt(filename, bio, ocn, fossil, fire, nx,ny,nt, sday)
      implicit none
      include 'netcdf.inc'
      integer nx, ny, nt
      character(len=8) sday
      integer date(nt), datesec(nt), time(nt)
      real fossil(nx,ny,nt), fire(nx,ny,nt), bio(nx,ny,nt), ocn(nx,ny,nt),lat(ny), lon(nx)
      integer ocnid, bioid, fossilid, fireid, latvarid, lonvarid, fid
      integer datevarid, datesecvarid, timevarid
      integer dims(3), start(3), count(3)
      integer stepid, latid, lonid
      character(len=100) filename
      integer rcode, status
      integer i, j, it
 
      do it=1, nt
           datesec(it)=(it-1)*10800+5400
           read(sday,*) date(it)
           time(it)=it
 !          print*, it, datesec(it), date(it), time(it) 
      enddo
      do i=1, nx
            lon(i)=i-0.5
      enddo
      do i=1, ny
            lat(i)=i-90.5
      enddo
     
        rcode = nf_create(trim(filename),NF_CLOBBER,fid)
        if (rcode /= NF_NOERR) then
           write(6,*) "ERROR: creating file ",trim(filename)
           write(6,*) "ERROR: ",nf_strerror(rcode)
           stop
        end if
       !... begin to define your dimentions
             status=nf_def_dim(fid, 'time', nf_unlimited, STEPID)  ! the time is unlimited in this program
     
             status=nf_def_dim(fid, 'lat', ny, LATID)
     
             status=nf_def_dim(fid, 'lon', nx, LONID)
     
             dims(3)=stepid
             dims(2)=latid
             dims(1)=lonid
     
          status=nf_def_var(fid, 'lat', nf_float, 1, latid, latvarid)
          status=nf_def_var(fid, 'lon', nf_float, 1, lonid, lonvarid)
          status=nf_def_var(fid, 'date', nf_int, 1, stepid, datevarid)
          status=nf_def_var(fid, 'datesec', nf_int, 1, stepid, datesecvarid)
          status=nf_def_var(fid, 'time', nf_int, 1, stepid, timevarid)
          status=nf_def_var(fid, 'bio', nf_float, 3, dims, bioid)
          status=nf_def_var(fid, 'ocn', nf_float, 3, dims, ocnid)
          status=nf_def_var(fid, 'fossil', nf_float, 3, dims, fossilid)
          status=nf_def_var(fid, 'fire', nf_float, 3, dims, fireid)
          status=nf_enddef(fid)
     
        count(3)=nt
        count(2)=ny
        count(1)=nx
        start=1
     
        status=nf_put_vara_real(fid, bioid, start, count, bio)
        status=nf_put_vara_real(fid, ocnid, start, count, ocn)
        status=nf_put_vara_real(fid, fossilid, start, count, fossil)
        status=nf_put_vara_real(fid, fireid, start, count, fire)
        status=nf_put_vara_real(fid, latvarid, 1, ny, lat)
        status=nf_put_vara_real(fid, lonvarid, 1, nx, lon)
        status=nf_put_vara_int(fid, datevarid, 1, nt, date)
        status=nf_put_vara_int(fid, datesecvarid, 1, nt, datesec)
        status=nf_put_vara_int(fid, timevarid, 1, nt, time)
        status=nf_close(fid)
     
     
     end subroutine

 subroutine writeopt_tmp(fname,bgopt,bgorg,dflux, lat,lon,nx,ny,nvar)
  implicit none
  include 'netcdf.inc'
  integer nx, ny, nvar, rundays, rcode, status
  real bgopt(nx,ny,nvar), bgorg(nx,ny,nvar), dflux(nx,ny,nvar)
  integer  fid, nvarid, latid, lonid, ensid
  character(len=100) fname
  character(len=30) var
  integer latvarid,lonvarid, bgensid, bgmeanid, bgorgid
  real lat(ny), lon(nx)
  integer start(3), count(3), dims(3)
  integer i, it

  rcode = nf_create(trim(fname),NF_CLOBBER,fid)
  if (rcode /= NF_NOERR) then
     write(6,*) "ERROR: creating file ",trim(fName)
     write(6,*) "ERROR: ",nf_strerror(rcode)
     stop
  end if
  !... begin to define your dimentions
       status=nf_def_dim(fid, 'var', nvar, nvarid)  ! the time is unlimited in this program
       if (status .NE. nf_noerr) call handle_err(status)
  !      print*,'stepid',stepid,nf_unlimited

       status=nf_def_dim(fid, 'lat', ny, LATID)
       if (status .NE. nf_noerr) call handle_err(status)
  !     print*,'latid',latid,ny

       status=nf_def_dim(fid, 'lon', nx, LONID)
       if (status .NE. nf_noerr) call handle_err(status)


  !    print*,'dimentions defining end'
      dims(3)=nvarid
      dims(2)=latid
      dims(1)=lonid

    var='bgopt'
    status=nf_def_var(fid, trim(var), nf_float, 3, dims(1:3), bgmeanid)
    if (status .NE. nf_noerr) call handle_err(status)
    status=nf_put_att_text(fid, bgmeanid, 'long_name', 19, 'optimze carbon flux')
    if (status .NE. nf_noerr) call handle_err(status)
    status=nf_put_att_text(fid, bgmeanid, 'units', 13, 'gC/m2/rundays')
    if (status .NE. nf_noerr) call handle_err(status)

    var='bgorg'
    status=nf_def_var(fid, trim(var), nf_float, 3, dims(1:3), bgorgid)
    if (status .NE. nf_noerr) call handle_err(status)
    status=nf_put_att_text(fid, bgorgid, 'long_name', 17, 'orgin carbon flux')
    if (status .NE. nf_noerr) call handle_err(status)
    status=nf_put_att_text(fid, bgorgid, 'units', 13, 'gC/m2/rundays')
    if (status .NE. nf_noerr) call handle_err(status)

    var='dflux'
    status=nf_def_var(fid, trim(var), nf_float, 3, dims(1:3), bgensid)
    if (status .NE. nf_noerr) call handle_err(status)
    status=nf_put_att_text(fid, bgensid, 'long_name', 17, 'delta carbon flux')
    if (status .NE. nf_noerr) call handle_err(status)
    status=nf_put_att_text(fid, bgensid, 'units', 13, 'gC/m2/rundays')
    if (status .NE. nf_noerr) call handle_err(status)

    status=nf_def_var(fid, 'lat', nf_float, 1, latid, latvarid)
    if (status .NE. nf_noerr) call handle_err(status)
    status=nf_put_att_text(fid, latvarid, 'long_name', 8, 'Latitude')
    if (status .NE. nf_noerr) call handle_err(status)

    status=nf_def_var(fid, 'lon', nf_float, 1, lonid, lonvarid)
    if (status .NE. nf_noerr) call handle_err(status)
    status=nf_put_att_text(fid, lonvarid, 'long_name', 9, 'Longitude')
    if (status .NE. nf_noerr) call handle_err(status)

    status=nf_enddef(fid)
    count(3)=nvar
    count(2)=ny
    count(1)=nx
    start=1
    status=nf_put_vara_real(fid, bgmeanid, start(1:3), count(1:3), bgopt)
    status=nf_put_vara_real(fid, bgorgid, start(1:3), count(1:3), bgorg)
    status=nf_put_vara_real(fid, bgensid, start(1:3), count(1:3),dflux)
    status=nf_put_vara_real(fid, latvarid, 1, ny, lat)
    status=nf_put_vara_real(fid, lonvarid, 1, nx, lon)

    status=nf_close(fid)    

end subroutine
     !writeopt_nrt(filename,flux,nx,ny,nt,sday)
     subroutine writeopt_nrt(filename, flux, nx,ny,nt, sday)
      include 'netcdf.inc'
      integer nx, ny, nt
      character(len=8) sday
      integer date(nt), datesec(nt), time(nt)
      real flux(nx,ny,nt),lat(ny), lon(nx)
      integer ocnid, bioid, fossilid, fireid, latvarid, lonvarid, fid
      integer datevarid, datesecvarid, timevarid
      integer dims(3), start(3), count(3)
      integer stepid, latid, lonid
      character(len=100) filename
      integer rcode, status
!      print*, nt, sday
      do it=1, nt
           datesec(it)=(it-1)*10800+5400
           read(sday,*) date(it)
           time(it)=it
 !          print*, datesec(it), date(it), time(it)
      enddo
      do i=1, nx
            lon(i)=i-0.5
      enddo
      do i=1, ny
            lat(i)=i-90.5
      enddo
     
        rcode = nf_create(trim(filename),NF_CLOBBER,fid)
        if (rcode /= NF_NOERR) then
           write(6,*) "ERROR: creating file ",trim(filename)
           write(6,*) "ERROR: ",nf_strerror(rcode)
           stop
        end if
       !... begin to define your dimentions
             status=nf_def_dim(fid, 'time', nf_unlimited, STEPID)  ! the time is unlimited in this program
     
             status=nf_def_dim(fid, 'lat', ny, LATID)
     
             status=nf_def_dim(fid, 'lon', nx, LONID)
     
             dims(3)=stepid
             dims(2)=latid
             dims(1)=lonid
     
          status=nf_def_var(fid, 'lat', nf_float, 1, latid, latvarid)
          status=nf_def_var(fid, 'lon', nf_float, 1, lonid, lonvarid)
          status=nf_def_var(fid, 'date', nf_int, 1, stepid, datevarid)
          status=nf_def_var(fid, 'datesec', nf_int, 1, stepid, datesecvarid)
          status=nf_def_var(fid, 'time', nf_int, 1, stepid, timevarid)
          status=nf_def_var(fid, 'netflux', nf_float, 3, dims, bioid)
          
          status=nf_enddef(fid)
     
        count(3)=nt
        count(2)=ny
        count(1)=nx
        start=1
     
        status=nf_put_vara_real(fid, bioid, start, count, flux)
        status=nf_put_vara_real(fid, latvarid, 1, ny, lat)
        status=nf_put_vara_real(fid, lonvarid, 1, nx, lon)
        status=nf_put_vara_int(fid, datevarid, 1, nt, date)
        status=nf_put_vara_int(fid, datesecvarid, 1, nt, datesec)
        status=nf_put_vara_int(fid, timevarid, 1, nt, time)
        status=nf_close(fid)
     
     
     end subroutine

     subroutine readopt(filename, bio, ocn, fossil, fire, nx,ny,nt)
      include 'netcdf.inc'
      integer nx, ny, nt, nvar
      integer date(nt), datesec(nt), time(nt)
      real fossil(nx,ny,nt), fire(nx,ny,nt), bio(nx,ny,nt), ocn(nx,ny,nt)
      integer ncid, varid
      character(len=100) filename
      integer retval
      character(len=50) code
      code='subroutine readopt'

      retval = nf_open(trim(FILENAME), NF_NOWRITE, ncid)
      if (retval .ne. nf_noerr) call handle_err(retval,code)

      retval = nf_inq_varid    ( ncid, "bio", varid )
      retval = nf_get_var_real ( ncid, varid, bio )
      retval = nf_inq_varid    ( ncid, "ocn", varid )
      retval = nf_get_var_real ( ncid, varid, ocn )
      retval = nf_inq_varid    ( ncid, "fossil", varid )
      retval = nf_get_var_real ( ncid, varid, fossil )
      retval = nf_inq_varid    ( ncid, "fire", varid )
      retval = nf_get_var_real ( ncid, varid, fire )
      retval = NF_CLOSE(NCID)
     
     
     end subroutine

     subroutine readopt_nrt(filename, flux, nx,ny,nt)
      implicit none
      include 'netcdf.inc'
      integer nx, ny, nt
      real flux(nx,ny,nt)
      integer varid, ncid
      character(len=100) filename
      integer retval
      character(len=50) code
      code='subroutine readopt_nrt'

      retval = nf_open(trim(FILENAME), NF_NOWRITE, ncid)
      if (retval .ne. nf_noerr) call handle_err(retval,code)

      retval = nf_inq_varid    ( ncid, "netflux", varid )
      retval = nf_get_var_real ( ncid, varid, flux )
      retval = NF_CLOSE(NCID)
     
     
     end subroutine

subroutine read_bgmean(filename,xorg,xens,nx,ny,nv,ensize)
  implicit none
  include 'netcdf.inc'
  character(len=100) filename
  integer nx,ny,nv, retval, varid, ncid, ensize, varid1
  real xorg(nx,ny,nv), xens(nx,ny,nv,ensize)
  character(len=50) code
     
     code='subroutine read_bgmean'
     !****************************************************
     !           open ncfiles
     !****************************************************
 !     print*, filename
   
     retval = nf_open(trim(filename), NF_NOWRITE, ncid)
     if (retval .ne. nf_noerr) call handle_err(retval,code)
!     print*, 'open ed' 
     retval = nf_inq_varid    ( ncid, "bgmean_org", varid )
     if (retval .ne. nf_noerr) call handle_err(retval,code)
 !    print*, varid

     retval = nf_get_var_real ( ncid, varid, xorg )
     if (retval .ne. nf_noerr) call handle_err(retval,code)

     retval = nf_inq_varid    ( ncid, "bgens", varid1 )
     if (retval .ne. nf_noerr) call handle_err(retval,code)
 !    print*, varid

     retval = nf_get_var_real ( ncid, varid1, xens )
     if (retval .ne. nf_noerr) call handle_err(retval,code)

   
     !close ncfile
     retval = nf_close(ncid)
end subroutine

subroutine read_uncertainty(filename,xens,nx,ny,nv,ensize)
     implicit none
     include 'netcdf.inc'
     character(len=100) filename
     integer nx,ny,nv, retval, varid, ncid, ensize
     real xens(nx,ny,nv,ensize)
     character(len=50) code

        code='subroutine read_uncertainty'
        !****************************************************
        !           open ncfiles
        !****************************************************
    !     print*, filename
      
        retval = nf_open(trim(filename), NF_NOWRITE, ncid)
        if (retval .ne. nf_noerr) call handle_err(retval,code)
   
        retval = nf_inq_varid    ( ncid, "bgens", varid )
        if (retval .ne. nf_noerr) call handle_err(retval,code)
    !    print*, varid
   
        retval = nf_get_var_real ( ncid, varid, xens )
        if (retval .ne. nf_noerr) call handle_err(retval,code)
   
      
        !close ncfile
        retval = nf_close(ncid)
   end subroutine

!call writeflux(FILE_NAME,year,month,day,rundays,co2,lat,lon,nx,ny)
subroutine writeflux_formzt(fname,year,month,day,rundays,co2,lat,lon,nx,ny)
      use module_global, only: nextday
      implicit none
      include 'netcdf.inc'
      integer year, month, day, rundays
      integer nx, ny, nt
      integer  fid, stepid, latid, lonid
      character(len=100) fname
      character(len=30) var
      integer varid, latvarid, lonvarid, dateid, datesecid, timeid
      real lat(ny), lon(nx), co2(nx,ny,rundays*8+1)
      integer date(rundays*8+1), time(rundays*8+1), datesec(rundays*8+1)
      character(len=8) sday
      integer start(3), count(3), dims(3)
      integer i, it, rcode, status
      integer sec(8)
      data sec/0, 10800, 21600, 32400, 43200, 54000, 64800, 75600/
      character(len=50) code

!     print*, co2(1:10,30,1)

      code='subroutine writeflux'
      nt=rundays*8+1
      write(sday,'(i4,i2.2,i2.2)') year, month, day
      do it=1, nt
          read(sday,*) date(it)
          time(it)=it
          if(mod(it,8)==0)then
          i=8
          call nextday(sday)
          else
          i=it/8
          i=i*8
          i=it-i
          endif
          datesec(it)=sec(i)
      !   print*, date(it), datesec(it), time(it)
      enddo
      rcode = nf_create(trim(fname),NF_CLOBBER,fid)
      if (rcode /= NF_NOERR) then
         write(6,*) "ERROR: creating file ",trim(fName)
         write(6,*) "ERROR: ",nf_strerror(rcode)
         stop
      end if
      !... begin to define your dimentions
           status=nf_def_dim(fid, 'time', nf_unlimited, STEPID)  ! the time is unlimited in this program
           if (status .NE. nf_noerr) call handle_err(status, code)
      !      print*,'stepid',stepid,nf_unlimited
    
           status=nf_def_dim(fid, 'lat', ny, LATID)
           if (status .NE. nf_noerr) call handle_err(status, code)
      !     print*,'latid',latid,ny
    
           status=nf_def_dim(fid, 'lon', nx, LONID)
           if (status .NE. nf_noerr) call handle_err(status,code)
      !     print*,'lonid',lonid,nx
    
      !    print*,'dimentions defining end'
          dims(3)=stepid
          dims(2)=latid
          dims(1)=lonid
        var='CO2'
        status=nf_def_var(fid, trim(var), nf_float, 3, dims, varid)
        if (status .NE. nf_noerr) call handle_err(status,code)
        status=nf_put_att_text(fid, varid, 'long_name', 17, 'total carbon flux')
        if (status .NE. nf_noerr) call handle_err(status,code)
        status=nf_put_att_text(fid, varid, 'units', 15, 'molecules/cm2/s')
        if (status .NE. nf_noerr) call handle_err(status,code)
    
        status=nf_def_var(fid, 'lat', nf_float, 1, latid, latvarid)
        if (status .NE. nf_noerr) call handle_err(status,code)
        status=nf_put_att_text(fid, latvarid, 'long_name', 8, 'Latitude')
        if (status .NE. nf_noerr) call handle_err(status,code)
    
        status=nf_def_var(fid, 'lon', nf_float, 1, lonid, lonvarid)
        if (status .NE. nf_noerr) call handle_err(status,code)
        status=nf_put_att_text(fid, lonvarid, 'long_name', 9, 'Longitude')
        if (status .NE. nf_noerr) call handle_err(status,code)
    
        status=nf_def_var(fid, 'date', nf_int, 1, stepid, dateid)
        if (status .NE. nf_noerr) call handle_err(status,code)
        status=nf_def_var(fid, 'datesec', nf_int, 1, stepid, datesecid)
        if (status .NE. nf_noerr) call handle_err(status,code)
        status=nf_def_var(fid, 'time', nf_int, 1, stepid, timeid)
        if (status .NE. nf_noerr) call handle_err(status,code)
        status=nf_enddef(fid)
    
        status=nf_put_vara_real(fid, latvarid, 1, ny, lat)
        status=nf_put_vara_real(fid, lonvarid, 1, nx, lon)
        count(3)=nt
        count(2)=ny
        count(1)=nx
        start(3)=1
        start(2)=1
        start(1)=1
        status=nf_put_vara_real(fid, varid, start, count, co2)
        status=nf_put_vara_int(fid, timeid, 1, nt, time(1:nt))
        status=nf_put_vara_int(fid, dateid, 1, nt, date(1:nt))
        status=nf_put_vara_int(fid, datesecid, 1, nt, datesec(1:nt))
        status=nf_close(fid)
    
    end subroutine
!writeflux_forassim
    subroutine writeflux_forassim(fname,bgens,bgmean,bgmean_org,bgstdev,lat,lon,nx,ny,ens,nvar,rundays)
      implicit none
      include 'netcdf.inc'
      integer nx, ny, ens, nvar, rundays, rcode, status
      real bgens(nx,ny,nvar,ens), bgmean(nx,ny,nvar), bgmean_org(nx,ny,nvar), bgstdev(nx,ny,nvar)
      integer  fid, nvarid, latid, lonid, ensid
      character(len=100) fname
      character(len=30) var
      integer latvarid,lonvarid, bgstdid, bgensid, bgmeanid, bgorgid
      real lat(ny), lon(nx)
      integer start(4), count(4), dims(4)
      integer i, it
      character(len=50) code
    
      code='subroutine writeflux_forassim'
    
      rcode = nf_create(trim(fname),NF_CLOBBER,fid)
      if (rcode /= NF_NOERR) then
         write(6,*) "ERROR: creating file ",trim(fName)
         write(6,*) "ERROR: ",nf_strerror(rcode)
         stop
      end if
      !... begin to define your dimentions
           status=nf_def_dim(fid, 'var', nvar, nvarid)  ! the time is unlimited in this program
           if (status .NE. nf_noerr) call handle_err(status,code)
      !      print*,'stepid',stepid,nf_unlimited
    
           status=nf_def_dim(fid, 'lat', ny, LATID)
           if (status .NE. nf_noerr) call handle_err(status,code)
      !     print*,'latid',latid,ny
    
           status=nf_def_dim(fid, 'lon', nx, LONID)
           if (status .NE. nf_noerr) call handle_err(status,code)
    
           status=nf_def_dim(fid, 'en_size', ens, ensid)
           if (status .NE. nf_noerr) call handle_err(status,code)
    
      !     print*,'lonid',lonid,nx
    
      !    print*,'dimentions defining end'
          dims(4)=ensid
          dims(3)=nvarid
          dims(2)=latid
          dims(1)=lonid
    
        var='bgens'
        status=nf_def_var(fid, trim(var), nf_float, 4, dims, bgensid)
        if (status .NE. nf_noerr) call handle_err(status,code)
        status=nf_put_att_text(fid, bgensid, 'long_name', 17, 'total carbon flux')
        if (status .NE. nf_noerr) call handle_err(status)
        status=nf_put_att_text(fid, bgensid, 'units', 15, 'gC/grid/rundays')
        if (status .NE. nf_noerr) call handle_err(status,code)
        
        var='bgstdev'
        status=nf_def_var(fid, trim(var), nf_float, 3, dims(1:3), bgstdid)
        if (status .NE. nf_noerr) call handle_err(status,code)
        status=nf_put_att_text(fid, bgensid, 'long_name', 17, 'total carbon flux')
        if (status .NE. nf_noerr) call handle_err(status)
        status=nf_put_att_text(fid, bgensid, 'units', 15, 'gC/grid/rundays')
        if (status .NE. nf_noerr) call handle_err(status,code)        

        var='bgmean'
        status=nf_def_var(fid, trim(var), nf_float, 3, dims(1:3), bgmeanid)
        if (status .NE. nf_noerr) call handle_err(status,code)
        status=nf_put_att_text(fid, bgmeanid, 'long_name', 17, 'total carbon flux')
        if (status .NE. nf_noerr) call handle_err(status,code)
        status=nf_put_att_text(fid, bgmeanid, 'units', 15, 'gC/grid/rundays')
        if (status .NE. nf_noerr) call handle_err(status,code)
    
        var='bgmean_org'
        status=nf_def_var(fid, trim(var), nf_float, 3, dims(1:3), bgorgid)
        if (status .NE. nf_noerr) call handle_err(status,code)
        status=nf_put_att_text(fid, bgorgid, 'long_name', 17, 'total carbon flux')
        if (status .NE. nf_noerr) call handle_err(status,code)
        status=nf_put_att_text(fid, bgorgid, 'units', 15, 'gC/grid/rundays')
        if (status .NE. nf_noerr) call handle_err(status,code)
    
        status=nf_def_var(fid, 'lat', nf_float, 1, latid, latvarid)
        if (status .NE. nf_noerr) call handle_err(status,code)
        status=nf_put_att_text(fid, latvarid, 'long_name', 8, 'Latitude')
        if (status .NE. nf_noerr) call handle_err(status,code)
    
        status=nf_def_var(fid, 'lon', nf_float, 1, lonid, lonvarid)
        if (status .NE. nf_noerr) call handle_err(status,code)
        status=nf_put_att_text(fid, lonvarid, 'long_name', 9, 'Longitude')
        if (status .NE. nf_noerr) call handle_err(status,code)
    
        status=nf_enddef(fid)
        count(4)=ens
        count(3)=nvar
        count(2)=ny
        count(1)=nx
        start=1
        status=nf_put_vara_real(fid, bgensid, start, count, bgens)
        status=nf_put_vara_real(fid, bgmeanid, start(1:3), count(1:3), bgmean)
        status=nf_put_vara_real(fid, bgorgid, start(1:3), count(1:3), bgmean_org)
        status=nf_put_vara_real(fid, bgstdid, start(1:3), count(1:3),bgstdev)
        status=nf_put_vara_real(fid, latvarid, 1, ny, lat)
        status=nf_put_vara_real(fid, lonvarid, 1, nx, lon)
    
        status=nf_close(fid)    
    
    end subroutine

    subroutine write_uncertainty(fname,bgens,lat,lon,nx,ny,nvar,ens)
     implicit none
     include 'netcdf.inc'
     integer nx, ny, ens, nvar, rcode, status
     real bgens(nx,ny,nvar,ens)
     real lat(ny), lon(nx)
     integer  fid, nvarid, latid, lonid, ensid
     character(len=100) fname
     character(len=30) var
     integer latvarid,lonvarid,  bgensid
     integer start(4), count(4), dims(4)
     integer i, it
     character(len=50) code
   
     code='subroutine write_uncertainty'
   
     rcode = nf_create(trim(fname),NF_CLOBBER,fid)
     if (rcode /= NF_NOERR) then
        write(6,*) "ERROR: creating file ",trim(fName)
        write(6,*) "ERROR: ",nf_strerror(rcode)
        stop
     end if
     !... begin to define your dimentions
          status=nf_def_dim(fid, 'var', nvar, nvarid)  ! the time is unlimited in this program
          if (status .NE. nf_noerr) call handle_err(status,code)
     !      print*,'stepid',stepid,nf_unlimited
   
          status=nf_def_dim(fid, 'lat', ny, LATID)
          if (status .NE. nf_noerr) call handle_err(status,code)
     !     print*,'latid',latid,ny
   
          status=nf_def_dim(fid, 'lon', nx, LONID)
          if (status .NE. nf_noerr) call handle_err(status,code)
   
          status=nf_def_dim(fid, 'en_size', ens, ensid)
          if (status .NE. nf_noerr) call handle_err(status,code)
   
     !     print*,'lonid',lonid,nx
   
     !    print*,'dimentions defining end'
         dims(4)=ensid
         dims(3)=nvarid
         dims(2)=latid
         dims(1)=lonid
   
       var='bgens'
       status=nf_def_var(fid, trim(var), nf_float, 4, dims, bgensid)
       if (status .NE. nf_noerr) call handle_err(status,code)
       status=nf_put_att_text(fid, bgensid, 'long_name', 17, 'total carbon flux')
       if (status .NE. nf_noerr) call handle_err(status)
       status=nf_put_att_text(fid, bgensid, 'units', 15, 'gC/grid/rundays')
       if (status .NE. nf_noerr) call handle_err(status,code)
   
       status=nf_def_var(fid, 'lat', nf_float, 1, latid, latvarid)
       if (status .NE. nf_noerr) call handle_err(status,code)
       status=nf_put_att_text(fid, latvarid, 'long_name', 8, 'Latitude')
       if (status .NE. nf_noerr) call handle_err(status,code)
   
       status=nf_def_var(fid, 'lon', nf_float, 1, lonid, lonvarid)
       if (status .NE. nf_noerr) call handle_err(status,code)
       status=nf_put_att_text(fid, lonvarid, 'long_name', 9, 'Longitude')
       if (status .NE. nf_noerr) call handle_err(status,code)
   
       status=nf_enddef(fid)
       count(4)=ens
       count(3)=nvar
       count(2)=ny
       count(1)=nx
       start=1
       status=nf_put_vara_real(fid, bgensid, start, count, bgens)
       status=nf_put_vara_real(fid, latvarid, 1, ny, lat)
       status=nf_put_vara_real(fid, lonvarid, 1, nx, lon)
   
       status=nf_close(fid)    
   
   end subroutine
!
subroutine handle_err(errcode, code)
      implicit none
      include 'netcdf.inc'
      integer errcode
      character(len=50) code
          
      print *, trim(code),' : ', nf_strerror(errcode)
      stop 2
end subroutine
    
