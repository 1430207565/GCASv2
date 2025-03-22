   program evaluation
   use module_global
   implicit none
   include 'netcdf.inc'
   integer, parameter :: nlevg=35,nt=8,nr=46, nobs_max=550
   integer, parameter :: nx=360, ny=180  !for global_map
   integer nyr, b(nx,ny,nr)
   real lat(ny), lon(nx)
   real, allocatable, dimension(:,:,:,:,:)  :: co2
   real, allocatable, dimension(:,:,:,:) :: Zm, rxco2, lxco2
   real, allocatable, dimension(:,:,:) :: psm, rxco2s, lxco2s
   REAL yobs, yobs_xco2
   real HX(2), Hx_xco2(2)
   real ybias(2), xbias(2)
   integer nsrf, nxco2
   real Ax(nlevg,2)
   integer NumOfObs
   integer iday, sdate
   integer find
   real secs(nt)
   character(len=8) sday
   character(len=25) varname
   character(len=30) regname(nr)
   character(len=100) filesim, fileobs, file_gosat
   integer nsite, glev
   real slat, slon, stime 
   real xco2(nobs_max),pxco2(nobs_max), xco2unc(nobs_max), ps(nobs_max)
   real p(nlevg,nobs_max), aj(nlevg,nobs_max),hj(nlevg,nobs_max),Yaj(nlevg,nobs_max)
   real latg(nobs_max),long(nobs_max),timeg(nobs_max)
   character(len=50) site_tmp
   integer  date_tmp, ix_tmp, jx_tmp, tx, kx_tmp
   real lat_tmp,sec_tmp,lon_tmp,hgt_tmp, obs_tmp, unc_tmp, Hx_tmp(2), Robs_tmp
   integer i, j, k
   integer sitevarid, latvarid, lonvarid, datevarid, secvarid, htvarid, obsvarid, uncvarid, simvarid
   integer tnobs
   integer status, rcode
   integer count(2), dims2(2), start(2)
   integer fid, nobsid, en_size_id, ii, jj
   integer iis, s, is, nobs
   integer, external :: dmin
!   real hyai(nlev+1), hybi(nlev+1)
!   real p1(nlev+1), pp(nlev)
   real, allocatable, dimension(:) :: hyai, hybi, p1, pp
   character(len=50) sitename(nobs_max)
   real conc(nobs_max,3,65000), time(nobs_max,65000), yo, yopt, yorg    !3125812
   integer no(nobs_max)
   integer date(nobs_max,65000), assim(nobs_max,65000), obsflag(nobs_max,65000)
   real lats(nobs_max,65000), lons(nobs_max,65000)
   real hgts(nobs_max,65000), uncs(nobs_max,65000)
   real ylats, ylons
   real p0
   integer ic, ec
   integer narg, tlen, ncycle
   integer obs_flag, is_assim, sample_flag
   CHARACTER(len=32) ::arg
   character(len=15) class1, class2
   integer syear, eyear, smon, emon,sm, em, iy, im, imm, nofmon
   real, external :: ave,meanbias, rmse, corr
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

   narg=IARGC()
  if(narg<2)then
    print*,'the input parameter error!'
   stop
  endif
  if(narg >0)then
    call getarg(1,arg)
    read(arg,*)icdate
    call getarg(2,arg)
    read(arg,*)tlen
  endif
 ! print*, icdate, tlen
  write(sday,'(i8)') icdate
  sdate=icdate
  
  call ensrf_init()

  ncycle=tlen/rundays
  tlen=rundays*ncycle
  tlen=tlen-1
  
  syear=sdate/10000
  smon=(sdate-syear*10000)/100
  do iday=1, tlen
    call nextday(sday)
  enddo 
  read(sday,'(i8)') icdate
  eyear=icdate/10000
  emon=(icdate-eyear*10000)/100
  nyr=eyear-syear+1
  if(eyear==syear)then
    nofmon=emon-smon+1
  elseif(eyear==syear+1)then
    nofmon=12-smon+1+emon
  else
    nofmon=12-smon+1+emon+(eyear-syear-1)*12
  endif
  write(sday,'(i8)') sdate
  print*, sday, syear, smon, eyear, emon, nofmon
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
   allocate(rxco2(nr,12,nyr,4))
   allocate(rxco2s(nr,nofmon,4))
   allocate(lxco2(ny,12,nyr,4))
   allocate(lxco2s(ny,nofmon,4))
   allocate(hyai(nlev+1))
   allocate(hybi(nlev+1))
   allocate(p1(nlev+1))
   allocate(pp(nlev))
   do i=1, nt
      secs(i)=(i-1)*10800+5400
   enddo

   rxco2=0
   lxco2=0
   nsrf=0
   nxco2=0
   yobs=0
   hx=0
   ybias=0
   yobs_xco2=0
   Hx_xco2=0
   xbias=0
   is=0
   open(100,file=trim(opt_dir)//'/'//trim(casename)//'/surface_site_evaluation.txt')
   open(110,file=trim(opt_dir)//'/'//trim(casename)//'/aircraft_site_evaluation.txt')
 !  if(issat==1)then
   open(200,file=trim(opt_dir)//'/'//trim(casename)//'/xco2_evaluation.txt')
 !  endif
   do iday=1, tlen
   !     print*, sday
      filesim='/data/home/TommyLv/GCASv3.2_fire/mozart4/run_forward/hist/'//trim(casename)//'.mz4.h0.'//sday(1:4)//'-'//sday(5:6)//'-'//sday(7:8)//'-10800.nc'
    print*, trim(filesim)
    inquire(file=trim(filesim),exist=alive)
    if(.not.alive) cycle  
    call readmzt(filesim,co2,Zm,psm,hyai,hybi,nlon,nlat,nlev,nt,2)
      print*, sday
! read surface CO2 data
!      fileobs=trim(obsdir_srf)//'/Evalu_All/srf_'//sday
!      fileobs='/share/home/zqjiang/GCASv3.2/input/obs/Amazon/srf_'//sday
!      fileobs='/share/home/zqjiang/datasets/obs/obsoco2mip/srf_'//sday
!      inquire(file=trim(fileobs),exist=alive)
!      if(alive)then
!         open(1,file=trim(fileobs))
!         do while (.true.)
!            read(1,*,end=100) site_tmp, date_tmp, sec_tmp, lat_tmp, lon_tmp,hgt_tmp, obs_tmp, unc_tmp, is_assim, obs_flag
      fileobs=trim(obsdir_srf)//'/EvaluObs/srf_'//sday
      inquire(file=trim(fileobs),exist=alive)
      if(alive)then
         open(1,file=trim(fileobs))
         do while (.true.)
            read(1,*,end=100) site_tmp, date_tmp, sec_tmp, lat_tmp,lon_tmp,hgt_tmp, obs_tmp
            if(lon_tmp<0)then
               lon_tmp=lon_tmp+360
            endif
            ix_tmp=dmin2(lon_tmp,xlon,nlon)
            jx_tmp=dmin2(lat_tmp,xlat,nlat)
            tx=dmin2(sec_tmp,secs,nt)
            kx_tmp=dmin2(hgt_tmp,Zm(ix_tmp,jx_tmp,:,tx),nlev)
            Hx_tmp(1:2)=co2(ix_tmp,jx_tmp,kx_tmp,tx,1:2)
            Hx_tmp(1:2)=Hx_tmp(1:2)*29/44*1.0e6
      !      write(*,'(a50,2f12.0,3f12.2)') site_tmp, date_tmp, sec_tmp, obs_tmp,Hx_tmp(1:2)
            if(trim(site_tmp)=='co2_con_aircraft-insitu_42_allvalid') cycle
            if(trim(site_tmp)=='act_aircraft-insitu_428_allvalid-c130') cycle
            if(is==0)then
                 is=is+1
                 sitename(is)=site_tmp
                 lats(is,1)=lat_tmp
                 lons(is,1)=lon_tmp
                 hgts(is,1)=hgt_tmp
                 uncs(is,1)=unc_tmp
                 assim(is,1)=is_assim
                 obsflag(is,1)=obs_flag
                 conc(is,1,1)=obs_tmp
                 conc(is,2,1)=Hx_tmp(1)
                 conc(is,3,1)=Hx_tmp(2)
                 date(is,1)=date_tmp
                 time(is,1)=sec_tmp
                 no(is)=1               
            else
                 find=0
                  do s=1, is
                 !    print*, trim(sitename(s))
                     if(trim(site_tmp)==trim(sitename(s)))then
                        no(s)=no(s)+1
                        iis=no(s)
                 !    if(trim(site_tmp)=='co2_alt_surface-insitu_6_allvalid')then
                 !       write(*,'(a30,i10,2f12.1,i10)') trim(site_tmp), date_tmp, sec_tmp, obs_tmp, iis
                 !    endif
                 !       if(iis>300000)then
                 !          print*, sitename(s), no(s)
                 !          stop
                 !       endif
                        find=1
                        lats(s,iis)=lat_tmp
                        lons(s,iis)=lon_tmp
                        hgts(s,iis)=hgt_tmp
                        uncs(s,iis)=unc_tmp
                        assim(s,iis)=is_assim
                        obsflag(s,iis)=obs_flag
                        conc(s,1,iis)=obs_tmp
                        conc(s,2,iis)=Hx_tmp(1)
                        conc(s,3,iis)=Hx_tmp(2)
                        date(s,iis)=date_tmp
                        time(s,iis)=sec_tmp
                 !        print*, s, iis, sitename(s), no(s)
                     !   pause
                        exit
                     endif
                 enddo
                 if(find==0)then
                  !  do s=1, is
                  !     print*, s, sitename(s), no(s)
                  !  enddo
                  !  print*, 'new ', site_tmp
                   ! pause
                    is=is+1
             !      print*, is
                    sitename(is)=site_tmp
                    lats(is,1)=lat_tmp
                    lons(is,1)=lon_tmp
                    hgts(is,1)=hgt_tmp
                    uncs(is,1)=unc_tmp
                    assim(is,1)=is_assim
                    obsflag(is,1)=obs_flag
                    conc(is,1,1)=obs_tmp
                    conc(is,2,1)=Hx_tmp(1)
                    conc(is,3,1)=Hx_tmp(2)
                    date(is,1)=date_tmp
                    time(is,1)=sec_tmp
                    no(is)=1
                 endif
            endif
        !    endif
         !   yobs=yobs+obs_tmp
         !   Hx(1:2)=Hx(1:2)+Hx_tmp(1:2)
            if(abs(Hx_tmp(1)-obs_tmp)<10)then
               yobs=yobs+obs_tmp
               Hx(1:2)=Hx(1:2)+Hx_tmp(1:2)
               ybias(1)=ybias(1)+Hx_tmp(1)-obs_tmp
               ybias(2)=ybias(2)+Hx_tmp(2)-obs_tmp
               nsrf=nsrf+1
            endif
         enddo
100 continue
         close(1)
      endif
      !if(issat==1)then   ! assimilation satellite data
      file_gosat=trim(obsdir_tansat)//'/TanSat_v2.'//sday//'.nc'
!      file_gosat=trim(obsdir_gosat)//'/acos_v90_'//sday//'.nc'
      print*, trim(file_gosat)
      inquire(file=trim(file_gosat),exist=alive)
      if(alive)then
            call readsat(file_gosat,latg,long,timeg,xco2,aj,pxco2,hj,Yaj,p,ps,xco2unc,nlevg,nobs_max,glev,nobs)
            print*,sday, ' obs number in satellite=', nobs
            do i=1, nobs
               lat_tmp=latg(i)
               lon_tmp=long(i)
               if(lon_tmp<0)then
                  lon_tmp=lon_tmp+360     !lon: 0-359
               endif
               hgt_tmp=0
               obs_tmp=xco2(i)
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
               if(abs(Hx_tmp(1)-obs_tmp)>6) cycle
               yobs_xco2=yobs_xco2+obs_tmp
               Hx_xco2(1:2)=Hx_xco2(1:2)+Hx_tmp(1:2)
               xbias(1)=xbias(1)+hx_tmp(1)-obs_tmp
               xbias(2)=xbias(2)+hx_tmp(2)-obs_tmp
               nxco2=nxco2+1
               ii=dmin2(long(i),lon,nx)
               jj=dmin2(latg(i),lat,ny)
               read(sday(1:4),*) iy
               iy=iy-syear+1
               read(sday(5:6),*) im

               lxco2(jj,im,iy,1)=lxco2(jj,im,iy,1)+obs_tmp      !obs
               lxco2(jj,im,iy,2)=lxco2(jj,im,iy,2)+hx_tmp(1)    !opt
               lxco2(jj,im,iy,3)=lxco2(jj,im,iy,3)+hx_tmp(2)    !org
               lxco2(jj,im,iy,4)=lxco2(jj,im,iy,4)+1               

               do k=1, nr
               if(b(ii,jj,k)==1)then
                   rxco2(k,im,iy,1)=rxco2(k,im,iy,1)+obs_tmp      !obs
                   rxco2(k,im,iy,2)=rxco2(k,im,iy,2)+hx_tmp(1)    !opt
                   rxco2(k,im,iy,3)=rxco2(k,im,iy,3)+hx_tmp(2)    !org
                   rxco2(k,im,iy,4)=rxco2(k,im,iy,4)+1
               endif
             enddo
               write(200,'(2f10.2,a10,4f10.2)') latg(i), long(i), sday, timeg(i), obs_tmp, Hx_tmp(1:2) 

           enddo
       endif
    !  endif
      call nextday(sday)
   enddo
      nsite=is
      yobs=yobs/nsrf
      Hx=hx/nsrf
      ybias=ybias/nsrf
   ! print*, ybias/nsrf
     imm=0
    do iy=syear, eyear
        if(syear==eyear)then
          sm=smon
          em=emon
        else
           if(iy==syear)then
              sm=smon
              em=12
           elseif(iy==eyear)then
              sm=1
              em=emon
           else
              sm=1
              em=12
           endif
        endif
       do im=sm, em
           imm=imm+1
           do k=1, nr
               if( rxco2(k,im,iy-syear+1,4)>0)then
                  rxco2(k,im,iy-syear+1,1)=rxco2(k,im,iy-syear+1,1)/rxco2(k,im,iy-syear+1,4)
                  rxco2(k,im,iy-syear+1,2)=rxco2(k,im,iy-syear+1,2)/rxco2(k,im,iy-syear+1,4)
                  rxco2(k,im,iy-syear+1,3)=rxco2(k,im,iy-syear+1,3)/rxco2(k,im,iy-syear+1,4)
               endif
           enddo
           do k=1, ny
               if(lxco2(k,im,iy-syear+1,4)>0)then
                  lxco2(k,im,iy-syear+1,1)=lxco2(k,im,iy-syear+1,1)/lxco2(k,im,iy-syear+1,4)
                  lxco2(k,im,iy-syear+1,2)=lxco2(k,im,iy-syear+1,2)/lxco2(k,im,iy-syear+1,4)
                  lxco2(k,im,iy-syear+1,3)=lxco2(k,im,iy-syear+1,3)/lxco2(k,im,iy-syear+1,4)
               endif
           enddo

           rxco2s(:,imm,:)=rxco2(:,im,iy-syear+1,:)
           lxco2s(:,imm,:)=lxco2(:,im,iy-syear+1,:)
           
       enddo
    enddo


   open(200,file=trim(opt_dir)//'/'//trim(casename)//'/latitude_xco2_evaluation.dat',form='unformatted',access='direct',recl=ny)
   k=0
   do im=1, nofmon
       k=k+1
       write(200,rec=k) lxco2s(:,im,2)-lxco2s(:,im,1)
       k=k+1
       write(200,rec=k) lxco2s(:,im,3)-lxco2s(:,im,1)
       k=k+1
       write(200,rec=k) lxco2s(:,im,4)
   enddo
   close(200)

    print*,'---------Surface Observation Eualuation Results---------'
    print*, 'Obs    Sim_opt   Sim_org  Bias_opt  Bias_org'
    write(*,'(5f8.2)') yobs,  hx(1:2), ybias(1:2)
     open(200,file=trim(opt_dir)//'/'//trim(casename)//'/Evaluation_results_stat.txt')
 
    do is=1, nsite
       open(1,file=trim(opt_dir)//'/'//trim(casename)//'/'//trim(sitename(is))//'.txt')
       do iis=1, no(is)
         if(abs(conc(is,1,iis)-conc(is,2,iis))<10)then
            write(1,'(2f10.2,f10.1,i10,f10.2,7f10.2,2i4)') lats(is,iis), lons(is,iis),hgts(is,iis), date(is,iis), time(is,iis), &
                   conc(is,1,iis), conc(is,2,iis), conc(is,3,iis), -999.0, -999.0,-999.0, uncs(is,iis), assim(is,iis), obsflag(is,iis)
         else
             write(1,'(2f10.2,f10.1,i10,f10.2,7f10.2,2i4)') lats(is,iis), lons(is,iis),hgts(is,iis), date(is,iis), time(is,iis), &
                   -999.0, -999.0, -999.0, conc(is,1,iis), conc(is,2,iis), conc(is,3,iis), uncs(is,iis), assim(is,iis), obsflag(is,iis)
         endif
       enddo
       close(1)
       yo=sum(conc(is,1,1:no(is)))/no(is)
       yopt=sum(conc(is,2,1:no(is)))/no(is)
       yorg=sum(conc(is,3,1:no(is)))/no(is)
       ylats=sum(lats(is,1:no(is)))/no(is)
       ylons=sum(lons(is,1:no(is)))/no(is)
       write(*,'(a35,2x,2f10.2,9f8.2)') sitename(is), ylats, ylons, ave(conc(is,1,1:no(is)),no(is)),  &
       ave(conc(is,2,1:no(is)),no(is)), ave(conc(is,3,1:no(is)),no(is)), &
       meanbias(conc(is,2,1:no(is)),conc(is,1,1:no(is)),no(is)), meanbias(conc(is,3,1:no(is)),conc(is,1,1:no(is)),no(is)),  &
       rmse(conc(is,2,1:no(is)),conc(is,1,1:no(is)),no(is)), rmse(conc(is,3,1:no(is)),conc(is,1,1:no(is)),no(is)), &
       corr(conc(is,2,1:no(is)),conc(is,1,1:no(is)),no(is)), corr(conc(is,3,1:no(is)),conc(is,1,1:no(is)),no(is))    
       
       write(200,'(a35,2x,2f10.2,9f8.2)') sitename(is), ylats, ylons, ave(conc(is,1,1:no(is)),no(is)),  & 
       ave(conc(is,2,1:no(is)),no(is)), ave(conc(is,3,1:no(is)),no(is)), &
       meanbias(conc(is,2,1:no(is)),conc(is,1,1:no(is)),no(is)),meanbias(conc(is,3,1:no(is)),conc(is,1,1:no(is)),no(is)),  &
       rmse(conc(is,2,1:no(is)),conc(is,1,1:no(is)),no(is)),rmse(conc(is,3,1:no(is)),conc(is,1,1:no(is)),no(is)), &
       corr(conc(is,2,1:no(is)),conc(is,1,1:no(is)),no(is)),corr(conc(is,3,1:no(is)),conc(is,1,1:no(is)),no(is))


    enddo
    write(200,'(a,2x,5f8.2)')'surf mean', yobs,  hx(1:2), ybias(1:2)
 

   ! if(issat==1)then      
    if(nxco2>0)then 
     yobs_xco2=yobs_xco2/nxco2
      Hx_xco2=hx_xco2/nxco2
      xbias=xbias/nxco2
      print*,'-------Satellite Observation Eualuation Results--------'
    print*, 'Obs    Sim_opt   Sim_org  Bias_opt  Bias_org'
    write(*,'(5f8.2)') yobs_xco2,  hx_xco2(1:2), xbias(1:2)
    write(200,'(a,2x,5f8.2)') 'XCO2 Mean', yobs_xco2,  hx_xco2(1:2), xbias(1:2)

   do k=1, nr
      write(*,'(a30,7f8.2)') trim(regname(k)), ave(rxco2s(k,:,1), nofmon),ave(rxco2s(k,:,2), nofmon),ave(rxco2s(k,:,3), nofmon),  &
                       meanbias(rxco2s(k,:,2),rxco2s(k,:,1), nofmon), meanbias(rxco2s(k,:,3),rxco2s(k,:,1), nofmon), &
                      rmse(rxco2s(k,:,2),rxco2s(k,:,1), nofmon), rmse(rxco2s(k,:,3),rxco2s(k,:,1), nofmon)
     write(200,'(a30,7f8.2)') trim(regname(k)), ave(rxco2s(k,:,1), nofmon),ave(rxco2s(k,:,2), nofmon),ave(rxco2s(k,:,2), nofmon),  &
                       meanbias(rxco2s(k,:,2),rxco2s(k,:,1), nofmon), meanbias(rxco2s(k,:,3),rxco2s(k,:,1), nofmon), &
                      rmse(rxco2s(k,:,2),rxco2s(k,:,1), nofmon), rmse(rxco2s(k,:,3),rxco2s(k,:,1), nofmon)

     open(1,file=trim(opt_dir)//'/'//trim(casename)//'/'//'GOSAT_xco2_'//trim(regname(k))//'.txt')
     do iy=syear, eyear
        if(syear==eyear)then
          sm=smon
          em=emon
        else
           if(iy==syear)then
              sm=smon
              em=12
           elseif(iy==eyear)then
              sm=1
              em=emon
           else
              sm=1
              em=12
           endif
        endif
       do im=sm, em
          write(1,'(i6,i4,3f10.2,3f10.0)') iy, im, rxco2(k,im,iy-syear+1,1:4) 
       enddo
    enddo

    close(1) 
  enddo

   endif
   ! endif
    close(200)

end program

Function corr(s, o, n)
 implicit none
 integer     :: n, m
 real        :: s(n),o(n), as, ao
 real        :: corr
 integer     :: i
 real        :: rivl1,rivl2,rivl3
 real        :: critv(4), Sr, t_value
 data critv/1.676,1.664,1.660,1.645/
 !size, 50, 80, 100, 150 
 corr=0

 as=0
 ao=0
 m=0
 do i=1, n
   if(abs(s(i)-o(i))<10)then
     as=as+s(i)
     ao=ao+o(i)
     m=m+1
   endif
 enddo
 as=as/m
 ao=ao/m

! cal coef 
 rivl1=0
 rivl2=0
 rivl3=0
 do i=1,n
    if(abs(s(i)-o(i))<10)then
       rivl1=rivl1+(s(i)-as)*(o(i)-ao)
       rivl2=rivl2+(s(i)-as)*(s(i)-as)
       rivl3=rivl3+(o(i)-ao)*(o(i)-ao)
    endif
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
     if(abs(s(i)-o(i))<10)then
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
    if(abs(s(i)-o(i))<10)then
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

subroutine read_map(b, nx, ny,nr)

 include 'netcdf.inc'

 integer nx, ny, nr
 real a(nx,ny), cn(nx,ny), lat(ny)
 integer b(nx,ny,nr)   !, zm(nx,ny,nz,nt)   !, zi(nx,ny,nz,nt)
 integer i, j, ncid
 character(len=30) fname
 character(len=25) vname

 do j=1, ny
    lat(j)=j-90.5
 enddo

 b=0

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
  do s=1, 2  !en_size
     write(varname,'(a2,i3.3,a9)')'EN',s,'_VMR_avrg'
!     print*, varname
     call get_4d_var(ncid,varname,co2(:,:,:,:,s),nx,ny,nz,nt)
!     print*, co2(1,1,1,1,s)
  enddo 
! co2(:,:,:,:,2)=co2(:,:,:,:,1)

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
