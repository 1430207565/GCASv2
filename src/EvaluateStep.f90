   program evaluation
   use module_global
   implicit none
   include 'netcdf.inc'
   integer, parameter :: nlevg=35,nt=8,nr=46, nobs_max=2000
   integer, parameter :: nx=360, ny=180  !for global_map
   integer nyr, b(nx,ny,nr)
   real lat(ny), lon(nx)
   real bio(nx,ny,nt), ocn(nx,ny,nt), fossil(nx,ny,nt), fire(nx,ny,nt), flux(nx,ny,nt)
   real obio(nx,ny,nt), oocn(nx,ny,nt), ofossil(nx,ny,nt), ofire(nx,ny,nt),oflux(nx,ny,nt)
   real bio_m(nr), obio_m(nr), ocn_m(nr), oocn_m(nr), fire_m(nr), ofire_m(nr), flux_m(nr), oflux_m(nr)
   real, allocatable, dimension(:,:,:,:,:)  :: co2
   real, allocatable, dimension(:,:,:,:) :: Zm
   real, allocatable, dimension(:,:,:) :: psm, bxco2, bco2
   real, allocatable, dimension(:,:)  :: rxco2, rco2
   
   REAL yobs, yobs_xco2
   real HX(2), Hx_xco2(2)
   real ybias(2), xbias(2)
   integer nsrf, nxco2
   real Ax(nlevg,2)
   integer NumOfObs
   integer iday, sdate
   integer find
   real secs(nt), fac
   character(len=8) sday
   character(len=25) varname
   character(len=30) regname(nr)
   character(len=100) filesim, fileobs, filename, file_gosat
   integer nobs, nsite, is, glev
   real slat, slon, stime 
   real xco2(nobs_max),pxco2(nobs_max), xco2unc(nobs_max), ps(nobs_max)
   real p(nlevg,nobs_max), aj(nlevg,nobs_max),hj(nlevg,nobs_max),Yaj(nlevg,nobs_max)
   real latg(nobs_max),long(nobs_max),timeg(nobs_max)
   character(len=30) site_tmp
   integer  date_tmp, ix_tmp, jx_tmp, tx, kx_tmp
   real lat_tmp,sec_tmp,lon_tmp,hgt_tmp, obs_tmp, unc_tmp, Hx_tmp(2), Robs_tmp
   integer s, i, j, t, k
   integer sitevarid, latvarid, lonvarid, datevarid, secvarid, htvarid, obsvarid, uncvarid, simvarid
   integer tnobs
   integer status, rcode
   integer count(2), dims2(2), start(2)
   integer fid, nobsid, en_size_id, iis, ii, jj
   integer, external :: dmin
!   real hyai(nlev+1), hybi(nlev+1)
!   real p1(nlev+1), pp(nlev)
   real, allocatable, dimension(:) :: hyai, hybi, p1, pp
   character(len=30) sitename(300)
   real conc(300,3,100000), time(300,100000), yo, yopt, yorg
   integer no(nobs_max), date(300,100000)
   real lats(300,100000), lons(300,100000)
   real ylats, ylons
   real p0
   integer ic, ec
   integer narg, tlen, ncycle
   integer obs_flag, sample_flag
   CHARACTER(len=32) ::arg
   character(len=15) class1, class2
   integer syear, eyear, smon, emon,sm, em, iy, im, imm, nofmon
   real, external :: ave,meanbias, rmse, corr, dxyp
  data regname/'North_America_Boreal','North_America_Temperate','South_America_Tropical','South_America_Temperate', &
              'North_Africa','South_Africa',    &
              'Asia_Boreal','Asia_Temperate','Southeast_Asia','Austrialia','Europe', &
              'North_Pacific','West_Pacific','East_Pacific','South_Pacific',  &
              'Arctic_Ocean','North_Atlantic','Tropical_Atlantic','South_Atlantic','Southern_Ocean', &
              'North_Indian_Ocean','South_Indian_Ocean',  &
              'China','North_America', 'South_America','Asia','Africa','Atlantic','Pacific',  &
              'Indian_Ocean',  &
              'Tropical_land','Tropical_ocean','Northern_Hemisphere_land','Northern_Hemisphere_ocean', &
              'Southern_Hemisphere_land','Southern_Hemisphere_Ocean','Global_land', &
              'Global_Ocean', 'Globe','USA','EU','India','Russia',  &
              'N._land_gt30N','T._land_30Nto30S','S._land_lt30S'/

  
  call ensrf_init()
  
  write(sday,'(i8)') icdate
  tlen=rundays
  
  do i=1, ny
     lat(i)=i-90.5
  enddo
  do i=1, nx
     lon(i)=i-0.5
  enddo  
 
  call read_map(b,nx,ny,nr) 
 
   p0=1.0e5
   
   allocate(co2(nlon,nlat,nlev,nt,2))
   allocate(Zm(nlon,nlat,nlev,nt))
   allocate(psm(nlon,nlat,nt))
   allocate(bxco2(bnx,bny,4))
   allocate(bco2(bnx,bny,4))
   allocate(rxco2(nr,4))
   allocate(rco2(nr,4))
   allocate(hyai(nlev+1))
   allocate(hybi(nlev+1))
   allocate(p1(nlev+1))
   allocate(pp(nlev))
   do i=1, nt
      secs(i)=(i-1)*10800+5400
   enddo

   rxco2=0
   bxco2=0
   rco2=0
   bco2=0
   nsrf=0
   nxco2=0
   yobs=0
   hx=0
   ybias=0
   yobs_xco2=0
   Hx_xco2=0
   xbias=0
   bio_m=0
   ocn_m=0
   fire_m=0
   flux_m=0
   obio_m=0
   oocn_m=0
   ofire_m=0
   oflux_m=0
   is=0
 !  if(issat==1)then
 !  endif
   do iday=1, tlen
   !     print*, sday

      call readprior(priordir,obio,oocn,ofossil,ofire,nx,ny,nt,sday,netflux)
      oflux=obio+oocn+ofossil+ofire
      oflux=oflux/6.02e23*12*1.0e4*3600
      obio=obio/6.02e23*12*1.0e4*3600
      oocn=oocn/6.02e23*12*1.0e4*3600
      ofossil=ofossil/6.02e23*12*1.0e4*3600
      ofire=ofire/6.02e23*12*1.0e4*3600

      filename=trim(opt_dir)//'/'//trim(casename)//'/emissions.opt.'//sday//'.nc'
    !  write(*,'(a,i4,2i2.2,a3)') 'emissions.opt.',y,m,d,'.nc'
      if(nrt==1)then
      call read_nc_nrt(filename,flux,nx,ny,nt)
      flux=flux/6.02e23*12*1.0e4*3600
      else
      call read_nc(filename,bio,ocn,fossil,fire,nx,ny,nt)
      ! molecules/cm2/s ==> gC/m2/hr
      bio=bio/6.02e23*12*1.0e4*3600
      ocn=ocn/6.02e23*12*1.0e4*3600
      fossil=fossil/6.02e23*12*1.0e4*3600
      fire=fire/6.02e23*12*1.0e4*3600
      flux=bio+ocn+fossil+fire
      endif
      do i=1, nx
         do j=1, ny
            
            do k=1, nr
              fac=b(i,j,k)*dxyp(lat(j),1.0)*1.0e6*3
           !   if(b(i,j,k)==1)then
                do t=1, nt
                 if(nrt==0)then
                    bio_m(k)=bio_m(k)+bio(i,j,t)*fac
                    ocn_m(k)=ocn_m(k)+ocn(i,j,t)*fac
                    fire_m(k)=fire_m(k)+fire(i,j,t)*fac

                    obio_m(k)=obio_m(k)+obio(i,j,t)*fac
                    oocn_m(k)=oocn_m(k)+oocn(i,j,t)*fac
                    ofire_m(k)=ofire_m(k)+ofire(i,j,t)*fac
                endif
                flux_m(k)=flux_m(k)+flux(i,j,t)*fac
                oflux_m(k)=oflux_m(k)+oflux(i,j,t)*fac
              enddo
          !  endif
          enddo
      enddo
    enddo

      filesim=trim(mozartdir)//'/run_forward/hist/'//trim(casename)//'.mz4.h0.'//sday(1:4)//'-'//sday(5:6)//'-'//sday(7:8)//'-10800.nc'
!    print*, trim(filesim)  
    call readmzt(filesim,co2,Zm,psm,hyai,hybi,nlon,nlat,nlev,nt,2)
   !   print*, sday
! read surface CO2 data
      fileobs=trim(obsdir_srf)//'/EvaluObs/srf_'//sday
      inquire(file=trim(fileobs),exist=alive)
      if(alive)then
         open(1,file=trim(fileobs))
         do while (.true.)
            read(1,*,end=100) site_tmp, date_tmp, sec_tmp, lat_tmp, lon_tmp,hgt_tmp, obs_tmp   !, obs_flag, sample_flag
            if(lon_tmp<0)then
               lon_tmp=lon_tmp+360
            endif
            obs_tmp=obs_tmp*1.00079-0.142    !scaled to WMO X2019
            ix_tmp=dmin2(lon_tmp,xlon,nlon)
            jx_tmp=dmin2(lat_tmp,xlat,nlat)
            tx=dmin2(sec_tmp,secs,nt)
            kx_tmp=dmin2(hgt_tmp,Zm(ix_tmp,jx_tmp,:,tx),nlev)
               Hx_tmp(1:2)=co2(ix_tmp,jx_tmp,kx_tmp,tx,1:2)
               Hx_tmp(1:2)=Hx_tmp(1:2)*29/44*1.0e6
!I want change to use 10
             if(abs(Hx_tmp(1)-obs_tmp)>15) cycle

            ii=dmin2(lon_tmp,blon,bnx)
            jj=dmin2(lat_tmp,blat,bny)
               bco2(ii,jj,1)=bco2(ii,jj,1)+obs_tmp      !obs
               bco2(ii,jj,2)=bco2(ii,jj,2)+hx_tmp(1)    !opt
               bco2(ii,jj,3)=bco2(ii,jj,3)+hx_tmp(2)    !org
               bco2(ii,jj,4)=bco2(ii,jj,4)+1

             ii=dmin2(lon_tmp,lon,nx)
             jj=dmin2(lat_tmp,lat,ny)
               do k=1, nr
               if(b(ii,jj,k)==1)then
                   rco2(k,1)=rco2(k,1)+obs_tmp      !obs
                   rco2(k,2)=rco2(k,2)+hx_tmp(1)    !opt
                   rco2(k,3)=rco2(k,3)+hx_tmp(2)    !org
                   rco2(k,4)=rco2(k,4)+1
               endif
             enddo
            yobs=yobs+obs_tmp
            Hx(1:2)=Hx(1:2)+Hx_tmp(1:2)
            ybias(1)=ybias(1)+Hx_tmp(1)-obs_tmp
            ybias(2)=ybias(2)+Hx_tmp(2)-obs_tmp
            nsrf=nsrf+1
         enddo
100 continue
         close(1)
      endif
      !if(issat==1)then   ! assimilation satellite data
      file_gosat=trim(obsdir_oco2)//'/oco2_v11_'//sday//'.nc'
      inquire(file=trim(file_gosat),exist=alive)
      if(alive)then
            call readsat(file_gosat,latg,long,timeg,xco2,aj,pxco2,hj,Yaj,p,ps,xco2unc,nlevg,nobs_max,glev,nobs)
           ! print*,sday, ' obs number in satellite=', nobs
            do i=1, nobs
               lat_tmp=latg(i)
               lon_tmp=long(i)
               if(lon_tmp<0)then
                  lon_tmp=lon_tmp+360     !lon: 0-359
               endif
               hgt_tmp=0
               obs_tmp=xco2(i)*1.00079-0.142   !scaled to WMO X2019
               unc_tmp=xco2unc(i)*1.9
               if(Robs_tmp<1.5) Robs_tmp=1.5
               ix_tmp=dmin2(long(i),xlon,nlon)
               jx_tmp=dmin2(latg(i),xlat,nlat)
               kx_tmp=-1
               tx=dmin2(timeg(i),secs,nt)
               p1=hyai*p0+hybi*psm(ix_tmp,jx_tmp,tx)
               pp=(p1(1:nlev)+p1(2:nlev+1))/2
               pp=pp/100
               call interplot(Ax(1:glev,1:2),p(1:glev,i),pp,    &
               co2(ix_tmp,jx_tmp,1:nlev,tx,1:2)*29/44*1.0e6, glev, nlev, 2) 
               Hx_tmp(1:2)=pxco2(i)
               do k=1, glev
                  Hx_tmp(1:2)=Hx_tmp(1:2)+hj(k,i)*aj(k,i)*(Ax(k,1:2)-Yaj(k,i))
               enddo
               if(abs(Hx_tmp(1)-obs_tmp)>10) cycle
               yobs_xco2=yobs_xco2+obs_tmp
               Hx_xco2(1:2)=Hx_xco2(1:2)+Hx_tmp(1:2)
               xbias(1)=xbias(1)+hx_tmp(1)-obs_tmp
               xbias(2)=xbias(2)+hx_tmp(2)-obs_tmp
               nxco2=nxco2+1
               ii=dmin2(long(i),lon,nx)
               jj=dmin2(latg(i),lat,ny)
               do k=1, nr
               if(b(ii,jj,k)==1)then
                   rxco2(k,1)=rxco2(k,1)+obs_tmp      !obs
                   rxco2(k,2)=rxco2(k,2)+hx_tmp(1)    !opt
                   rxco2(k,3)=rxco2(k,3)+hx_tmp(2)    !org
                   rxco2(k,4)=rxco2(k,4)+1
               endif
             enddo
               ii=dmin2(long(i),blon,bnx)
               jj=dmin2(latg(i),blat,bny)
               bxco2(ii,jj,1)=bxco2(ii,jj,1)+obs_tmp      !obs
               bxco2(ii,jj,2)=bxco2(ii,jj,2)+hx_tmp(1)    !opt
               bxco2(ii,jj,3)=bxco2(ii,jj,3)+hx_tmp(2)    !org
               bxco2(ii,jj,4)=bxco2(ii,jj,4)+1

           enddo
       endif
    !  endif
      call nextday(sday)
   enddo
      nsite=is
    if(nsrf>0)then
      yobs=yobs/nsrf
      Hx=hx/nsrf
      ybias=ybias/nsrf
    else
      yobs=0
      Hx=0
      ybias=0
    endif
    if(nxco2>0)then
      yobs_xco2=yobs_xco2/nxco2
      Hx_xco2=hx_xco2/nxco2
      xbias=xbias/nxco2
    else
      yobs_xco2=0
      Hx_xco2=0
      xbias=0
    endif
   
   bio_m=bio_m/1.0e12
   obio_m=obio_m/1.0e12
   ocn_m=ocn_m/1.0e12
   oocn_m=oocn_m/1.0e12
   fire_m=fire_m/1.0e12
   ofire_m=ofire_m/1.0e12
   flux_m=flux_m/1.0e12
   oflux_m=oflux_m/1.0e12

   ! print*, ybias/nsrf
     do k=1, nr
        if( rxco2(k,4)>0)then
             rxco2(k,1)=rxco2(k,1)/rxco2(k,4)
             rxco2(k,2)=rxco2(k,2)/rxco2(k,4)
             rxco2(k,3)=rxco2(k,3)/rxco2(k,4)
        else
             rxco2(k,1)=0
             rxco2(k,2)=0
             rxco2(k,3)=0
        endif
        if( rco2(k,4)>0)then
             rco2(k,1)=rco2(k,1)/rco2(k,4)
             rco2(k,2)=rco2(k,2)/rco2(k,4)
             rco2(k,3)=rco2(k,3)/rco2(k,4)
        else
             rco2(k,1)=0
             rco2(k,2)=0
             rco2(k,3)=0
        endif
     enddo
   
     do ii=1, bnx
          do jj=1, bny
            if( bxco2(ii,jj,4)>0)then
             bxco2(ii,jj,1)=bxco2(ii,jj,1)/bxco2(ii,jj,4)
             bxco2(ii,jj,2)=bxco2(ii,jj,2)/bxco2(ii,jj,4)
             bxco2(ii,jj,3)=bxco2(ii,jj,3)/bxco2(ii,jj,4)
            else
              bxco2(ii,jj,1:3)=0
            endif
            if( bco2(ii,jj,4)>0)then
             bco2(ii,jj,1)=bco2(ii,jj,1)/bco2(ii,jj,4)
             bco2(ii,jj,2)=bco2(ii,jj,2)/bco2(ii,jj,4)
             bco2(ii,jj,3)=bco2(ii,jj,3)/bco2(ii,jj,4)
            else
             bco2(ii,jj,1:3)=0
            endif
          enddo
     enddo

   if(da_step==1)then
     open(1,file=trim(opt_dir)//'/'//trim(casename)//'/Global_Evaluations_forSteps.txt')
     write(1,*)'Datetime  SurfBias_opt  SurfBias_org   Xco2Bias_opt  Xco2Bias_org'
   else
     open(1,file=trim(opt_dir)//'/'//trim(casename)//'/Global_Evaluations_forSteps.txt',position='append')
   endif

    print*,'---------Global Observation Eualuation Results---------'
    print*, 'Datetime  SurfBias_opt  SurfBias_org   Xco2Bias_opt  Xco2Bias_org '
    write(*,'(i9,2x,f12.2,2x,f12.2,2x,f12.2,2x,f12.2)') icdate, ybias(1:2), xbias(1:2)
    write(1,'(i9,2x,f12.2,2x,f12.2,2x,f12.2,2x,f12.2)') icdate,  ybias(1:2), xbias(1:2)
    close(1)
   ! if(issat==1)then

if(nrt==0)then

    print*,'---------Regional Observation Eualuation Results---------'      
    write(*,'(a)') '                  Regname   surf_opt surf_org xco2_opt xco2_org biof_opt biof_org ocnf_opt ocnf_org fire_opt fire_org flux_opt flux_org'
     do k=1, nr
       write(*,'(a25,2x,4f9.2,8f9.3)') trim(regname(k)), rco2(k,2)-rco2(k,1),rco2(k,3)-rco2(k,1), rxco2(k,2)-rxco2(k,1),  &
      rxco2(k,3)-rxco2(k,1),bio_m(k), obio_m(k), ocn_m(k), oocn_m(k), fire_m(k), ofire_m(k), flux_m(k), oflux_m(k)
       if(da_step==1)then
         open(1,file=trim(opt_dir)//'/'//trim(casename)//'/'//trim(regname(k))//'_Evaluations_forSteps.txt')
         write(1,'(a)') ' Datetime  surf_opt surf_org xco2_opt xco2_org biof_opt biof_org ocnf_opt ocnf_org fire_opt fire_org flux_opt flux_org'
       else
         open(1,file=trim(opt_dir)//'/'//trim(casename)//'/'//trim(regname(k))//'_Evaluations_forSteps.txt',position='append')
       endif
       write(1,'(i9,2x,4f9.2,8f9.3)') icdate, rco2(k,2)-rco2(k,1),rco2(k,3)-rco2(k,1), rxco2(k,2)-rxco2(k,1),rxco2(k,3)-rxco2(k,1),  &
               bio_m(k), obio_m(k), ocn_m(k), oocn_m(k), fire_m(k), ofire_m(k), flux_m(k), oflux_m(k)
       close(1)
     enddo

else

print*,'---------Regional Observation Eualuation Results---------'
    write(*,'(a)') '                  Regname   surf_opt surf_org xco2_opt xco2_org flux_opt flux_org'
     do k=1, nr
       write(*,'(a25,2x,4f9.2,2f9.3)') trim(regname(k)), rco2(k,2)-rco2(k,1),rco2(k,3)-rco2(k,1), rxco2(k,2)-rxco2(k,1),  &
      rxco2(k,3)-rxco2(k,1), flux_m(k), oflux_m(k)
       if(da_step==1)then
         open(1,file=trim(opt_dir)//'/'//trim(casename)//'/'//trim(regname(k))//'_Evaluations_forSteps.txt')
         write(1,'(a)') ' Datetime  surf_opt surf_org xco2_opt xco2_org flux_opt flux_org'
       else
         open(1,file=trim(opt_dir)//'/'//trim(casename)//'/'//trim(regname(k))//'_Evaluations_forSteps.txt',position='append')
       endif
       write(1,'(i9,2x,4f9.2,2f9.3)') icdate, rco2(k,2)-rco2(k,1),rco2(k,3)-rco2(k,1), rxco2(k,2)-rxco2(k,1),rxco2(k,3)-rxco2(k,1),  &
               flux_m(k), oflux_m(k)
       close(1)
     enddo

endif


     filename =trim(opt_dir)//'/'//trim(casename)//'/emissions.opt.lambda.nc' 
     call writeEvaluxco2(filename,bxco2(:,:,2)-bxco2(:,:,1),bxco2(:,:,3)-bxco2(:,:,1),bco2(:,:,2)-bco2(:,:,1),bco2(:,:,3)-bco2(:,:,1),bnx,bny,da_step) 

end program


!call writelambda(filename,lambda,bnx,bny,bnvar,sday, da_step)
subroutine writeEvaluxco2(filename,xbias_opt,xbias_org,sbias_opt,sbias_org,nx,ny,step)
 implicit none
 include 'netcdf.inc'
 integer nx, ny, step
 real xbias_opt(nx,ny), xbias_org(nx,ny)
 real sbias_opt(nx,ny), sbias_org(nx,ny)
 character(len=100) filename
 character(len=50) code
 integer start(3), count(3), ncid, varid, retval
 integer i, j, k

 retval = nf_open(trim(filename), NF_WRITE, ncid)
  count(3)=1
  count(2)=ny
  count(1)=nx
  start(3)=step
  start(1:2)=1

  retval = nf_inq_varid    ( ncid, "xbias_org", varid )
  retval = nf_put_vara_real ( ncid, varid, start, count, xbias_org )

  retval = nf_inq_varid    ( ncid, "xbias_opt", varid )
  retval = nf_put_vara_real ( ncid, varid, start, count, xbias_opt )

   retval = nf_inq_varid    ( ncid, "sbias_org", varid )
  retval = nf_put_vara_real ( ncid, varid, start, count, sbias_org )

  retval = nf_inq_varid    ( ncid, "sbias_opt", varid )
  retval = nf_put_vara_real ( ncid, varid, start, count, sbias_opt )

  retval = NF_CLOSE(NCID)


end subroutine


Function corr(s, o, n)
 implicit none
 integer     :: n
 real        :: s(n),o(n), as, ao
 real        :: corr
 integer     :: i
 real        :: rivl1,rivl2,rivl3
 real        :: critv(4), Sr, t_value
 data critv/1.676,1.664,1.660,1.645/
 !size, 50, 80, 100, 150 
 corr=0

 as=sum(s(1:n))/n
 ao=sum(o(1:n))/n

! cal coef 
 rivl1=0
 rivl2=0
 rivl3=0
 do i=1,n
       rivl1=rivl1+(s(i)-as)*(o(i)-ao)
       rivl2=rivl2+(s(i)-as)*(s(i)-as)
       rivl3=rivl3+(o(i)-ao)*(o(i)-ao)
 enddo

 rivl2=sqrt(rivl2)*sqrt(rivl3)
 if(rivl2>0)then
   corr=rivl1/rivl2
 endif
 return

end function

Function rmse(s,o, n)
  implicit none
  integer n
  real s(n), o(n), rmse
  integer i, m

  rmse=0
  m=0
  do i=1, n
     if(s(i)>0)then
     rmse=rmse+(s(i)-o(i))*(s(i)-o(i))
     m=m+1
    endif
  enddo
  if(m>0)then 
  rmse=sqrt(rmse/m)
  else
  rmse=0
  endif
end function

Function meanbias(s,o,n)

  implicit none
  integer n
  real s(n), o(n), meanbias
  integer i,m
  m=0
  meanbias=0
  do i=1, n
    if(s(i)>0)then
        meanbias=meanbias+s(i)-o(i)
        m=m+1
    endif
  enddo
  if(m>0)then
  meanbias=meanbias/m
  else
  meanbias=0
  endif
end Function

Function ave(s,n)

  implicit none
  integer n
  real s(n), ave
  integer i,m
  m=0
  ave=0
  do i=1, n
    if(s(i)>0)then
        ave=ave+s(i)
        m=m+1
    endif
  enddo
  if(m>0)then
  ave=ave/m
  else
  ave=0
  endif
end Function

Function dxyp( lat, re)
    implicit none
    real, parameter :: ae=6.371e3, pi=3.14159265358979
    real, parameter :: gtor=pi/180
    real re, lat, l
    real                  ::  dxx,dyy
    real dxyp

    dxx = 1.0*gtor*re
    dyy = 1.0*gtor*re
    l = (lat-re/2)*gtor

    dxyp = dxx * (sin(l+dyy)-sin(l))*ae**2

    return

End Function

subroutine read_map(b, nx, ny,nr)

 include 'netcdf.inc'

 integer nx, ny, nr
 real a(nx,ny), cn(nx,ny), lat(ny)
 integer b(nx,ny,nr)   !, zm(nx,ny,nz,nt)   !, zi(nx,ny,nz,nt)
 integer i, j, ncid
 character(len=30) fname
 character(len=25) vname

 b=0

 do j=1, ny
    lat(j)=j-90.5
 enddo

 fname='../input/regions.nc'
! print*, fname, nx, ny, nr
 status=nf_open(trim(fname),nf_nowrite,ncid)
   if ( status/=nf_noerr ) then
      write (*,*) nf_strerror(status)
      write (*,*) '!!Open NetCDF file Error!!'
       stop
    end if
 vname='region'
  call get_2d_var(ncid,vname,a,nx,ny)
 print*, 'read region ok'

 status = nf_close(ncid)

 fname='../input/region_top5.nc'
! print*, fname, nx, ny, nr
 status=nf_open(trim(fname),nf_nowrite,ncid)
   if ( status/=nf_noerr ) then
      write (*,*) nf_strerror(status)
      write (*,*) '!!Open NetCDF file Error!!'
       stop
    end if

 vname='region'
  call get_2d_var(ncid,vname,cn,nx,ny)
!print*,'read china ok'

 status = nf_close(ncid)

do i=1, nx
     do j=1, ny
       if(cn(i,j)==1)then
        b(i,j,23)=1
       elseif(cn(i,j)==2)then
        b(i,j,40)=1
       elseif(cn(i,j)==3)then
        b(i,j,41)=1
       elseif(cn(i,j)==4)then
        b(i,j,42)=1
       elseif(cn(i,j)==5)then
        b(i,j,43)=1
       endif
       if(a(i,j)>0.and.a(i,j)<12)then
          if(lat(j)>=30)then
             b(i,j,44)=1
          elseif(lat(j)>-30.and.lat(j)<30)then
             b(i,j,45)=1
          elseif(lat(j)<30)then
             b(i,j,46)=1
          endif
       endif
     enddo
 enddo



 do i=1, nx
  do j=1, ny
     k=int(a(i,j))
 !    print*, i,j,k
     if(k>0.and.k<23)then
        b(i,j,k)=1
     endif
  enddo
 enddo


 b(:,:,24)=b(:,:,1)+b(:,:,2)
 b(:,:,25)=b(:,:,3)+b(:,:,4)
 b(:,:,26)=b(:,:,7)+b(:,:,8)+b(:,:,9)
 b(:,:,27)=b(:,:,5)+b(:,:,6)
 b(:,:,28)=b(:,:,17)+b(:,:,18)+b(:,:,19)
 b(:,:,29)=b(:,:,12)+b(:,:,13)+b(:,:,14)+b(:,:,15)
 b(:,:,30)=b(:,:,21)+b(:,:,22)
 b(:,:,31)=b(:,:,5)+b(:,:,6)+b(:,:,9)+b(:,:,3)
 b(:,:,32)=b(:,:,21)+b(:,:,13)+b(:,:,14)+b(:,:,18)
 b(:,:,33)=b(:,:,1)+b(:,:,2)+b(:,:,7)+b(:,:,8)+b(:,:,11)
 b(:,:,34)=b(:,:,12)+b(:,:,17)+b(:,:,16)
 b(:,:,35)=b(:,:,4)+b(:,:,10)
 b(:,:,36)=b(:,:,22)+b(:,:,15)+b(:,:,19)+b(:,:,20)
 b(:,:,37)=b(:,:,31)+b(:,:,33)+b(:,:,35)
 b(:,:,38)=b(:,:,32)+b(:,:,34)+b(:,:,36)
 b(:,:,39)=b(:,:,37)+b(:,:,38)

 !print*, 'end'

end subroutine


subroutine interplot(Ax,p1,p2,x2,n,m,en)
  implicit none
  integer n,m,en
  real Ax(n,en), p1(n), p2(m), x2(m,en)
  integer i,j
  real f
  Ax=0

  do i=1,n
    if(p1(i)<p2(1)) then
        Ax(i,:)=x2(1,:)
    elseif(p1(i)>p2(m)) then
        Ax(i,:)=x2(m,:)
    else
      do j=2,m

        if(p1(i)>p2(j-1) .and. p1(i)<=p2(j)) then
           f=(p2(j)-p1(i))/(p1(i)-p2(j-1))
          Ax(i,:)=(x2(j,:)+f*x2(j-1,:))/(1+f)
        endif
      enddo
    endif
enddo

end subroutine

subroutine readsat(filename,lat,lon,time,xco2,aj,pxco2,hj,Yaj,p,ps,xco2unc,nlevin,nobsin,nlev,nobs)

include 'netcdf.inc'
  integer nlevin, nobsin, nlev,nobs
  real lat(nobsin), lon(nobsin), time(nobsin)
  real xco2(nobsin),pxco2(nobsin), xco2unc(nobsin), ps(nobsin)
  real p(nlevin,nobsin), aj(nlevin,nobsin),hj(nlevin,nobsin),Yaj(nlevin,nobsin)
  character(len=25) varname
  character(len=100) filename

  status=nf_open(trim(filename),nf_nowrite,ncid)
   if ( status/=nf_noerr ) then
      write (*,*) nf_strerror(status)
      write (*,*) '!!Open NetCDF file Error!!'
       stop
    end if

  status=nf_inq_dimid(ncid, 'Level', dimid)
  if ( status/=nf_noerr ) write (*,*) nf_strerror(status)

  status=nf_inq_dimlen(ncid, dimid, nlev)
  if ( status/=nf_noerr ) write (*,*) nf_strerror(status)

  status=nf_inq_dimid(ncid, 'nSamples', dimid)
  if ( status/=nf_noerr ) write (*,*) nf_strerror(status)

  status=nf_inq_dimlen(ncid, dimid, nobs)
  if ( status/=nf_noerr ) write (*,*) nf_strerror(status)

  varname='latitude'
  call get_1d_var(ncid,varname,lat(1:nobs),nobs)
  varname='longitude'
  call get_1d_var(ncid,varname,lon(1:nobs),nobs)
  where(lon<0)lon=lon+360
  varname='xco2'
  call get_1d_var(ncid,varname,xco2(1:nobs),nobs)
  varname='time'
  call get_1d_var(ncid,varname,time(1:nobs),nobs)
  time=time*3600

  varname='averagingKernel'
  call get_2d_var(ncid,varname,aj(1:nlev,1:nobs),nlev,nobs)
  varname='axco2'
  call get_1d_var(ncid,varname,pxco2(1:nobs),nobs)
  varname='weightingfunction'
  call get_2d_var(ncid,varname,hj(1:nlev,1:nobs),nlev,nobs)
  varname='constraintvector'
  call get_2d_var(ncid,varname,Yaj(1:nlev,1:nobs),nlev,nobs)
  varname='pressure'
  call get_2d_var(ncid,varname,p(1:nlev,1:nobs),nlev,nobs)
  varname='surfacepressure'
  call get_1d_var(ncid,varname,ps(1:nobs),nobs)
  varname='xco2unc'
  call get_1d_var(ncid,varname,xco2unc(1:nobs),nobs)
  status = nf_close(ncid)


end subroutine



integer function dmin(a, b, n)
 implicit none
 integer n
 real b(n), a
 real d, d1
 integer i
 d=abs(b(2)-b(1))/2

 do i=1, n
    d1=abs(a-b(i))
    if(d1<d)then
        dmin=i
        return
    endif
 enddo
 
end function

subroutine readmzt(filename,co2,Zm,psm,hyai,hybi,nx,ny,nz,nt,en_size)
  implicit none
  include 'netcdf.inc'
  integer nx, ny, nz, nt, en_size,s
  real co2(nx,ny,nz,nt,en_size)
  real Zm(nx,ny,nz,nt), psm(nx,ny,nt)
  real hyai(nz+1), hybi(nz+1)
  character(len=100) filename
  character(len=25) varname
  integer status, ncid

  status=nf_open(trim(filename),nf_nowrite,ncid)
   if ( status/=nf_noerr ) then
      write (*,*) nf_strerror(status)
      write (*,*) '!!Open NetCDF file Error!!'
       stop
    end if
  varname='PS'
  call get_3d_var(ncid,varname,psm,nx,ny,nt)
  varname='Zm'
  call get_4d_var(ncid,varname,Zm,nx,ny,nz,nt)
  varname='hyai'
  call get_1d_var(ncid,varname,hyai,nz+1)
  varname='hybi'
  call get_1d_var(ncid,varname,hybi,nz+1)
  do s=1, en_size
     write(varname,'(a2,i3.3,a9)')'EN',s,'_VMR_avrg'
!     print*, varname
     call get_4d_var(ncid,varname,co2(:,:,:,:,s),nx,ny,nz,nt)
!     print*, co2(1,1,1,1,s)
  enddo 
  status = nf_close(ncid)
 
end subroutine

subroutine readnc4d(filename,varname,a,nx,ny,nz,nt)
 include 'netcdf.inc'
  integer nx,ny,nz,nt
  real a(nx,ny,nz,nt)
  character(len=25) varname
  character(len=100) filename

  status=nf_open(trim(filename),nf_nowrite,ncid)
   if ( status/=nf_noerr ) then
      write (*,*) nf_strerror(status)
      write (*,*) '!!Open NetCDF file Error!!'
       stop
    end if
  call get_4d_var(ncid,varname,a,nx,ny,nz,nt)
  status = nf_close(ncid)
 

end subroutine

subroutine get_4d_var(ncid,vname,a,nx,ny,nz,nt)
  include 'netcdf.inc'
  integer nx,ny,ncid,nt,nz
  real a(nx,ny,nz,nt)
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


subroutine get_1d_var(ncid,vname,a,n)
  include 'netcdf.inc'
  integer n,ncid
  real a(n)
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


