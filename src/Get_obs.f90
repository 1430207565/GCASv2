subroutine get_obs(yobs,Robs,Hx,ylat,ylon,kx,yhgt,ytype,NumOfObs)
   use module_global
   implicit none
   include 'netcdf.inc'
   integer, parameter :: nlevg=35,nt=8, nobs_max=5000
   integer, parameter :: nx=360, ny=180
   real, parameter :: local_scale=30.0, large_scale=200.0
   real, allocatable, dimension(:,:,:,:,:)  :: co2
   real, allocatable, dimension(:,:,:,:) :: Zm
   real, allocatable, dimension(:,:,:) :: psm
   real, allocatable, dimension(:,:)  :: bxco2, bnxco2
   REAL,dimension(maxobs)           ::  yobs,Robs, ylat, ylon, yhgt
   integer, dimension(maxobs)       :: ytype
   real HX(maxobs,ensemble)
   integer no(maxobs), ix(maxobs), jx(maxobs), kx(maxobs)
   real  err1
   real Ax(nlevg,ensemble)
   integer NumOfObs
   integer iday, sdate
   integer find, ic, ec
   real secs(nt)
   character(len=8) sday
   character(len=25) varname
   character(len=100) filesim, fileobs, fileout, file_gosat, file_oco2, file_tansat, filename
   integer nobs,nobs_tmp, is, glev_tmp, glev(nobs_max)
   real L_obs, L_model(nlat)
   real unc_obs, unc_mis, abias_srf, abias_sat
   real slat, slon, stime 

   real xco2(nobs_max),pxco2(nobs_max), xco2unc(nobs_max), ps(nobs_max)
   real p(nlevg,nobs_max), aj(nlevg,nobs_max),hj(nlevg,nobs_max),Yaj(nlevg,nobs_max)
   real latg(nobs_max),long(nobs_max),timeg(nobs_max)

   real xco2_tmp(nobs_max),pxco2_tmp(nobs_max), xco2unc_tmp(nobs_max), ps_tmp(nobs_max)
   real p_tmp(nlevg,nobs_max),aj_tmp(nlevg,nobs_max),hj_tmp(nlevg,nobs_max),Yaj_tmp(nlevg,nobs_max)
   real latg_tmp(nobs_max),long_tmp(nobs_max),timeg_tmp(nobs_max)

   integer glmap(nx, ny)
   character(len=30) site_tmp
   integer  date_tmp, ix_tmp, jx_tmp, tx, kx_tmp
   real sim_tmp,lat_tmp,sec_tmp,lon_tmp,hgt_tmp, obs_tmp, unc_tmp, Hx_tmp(ensemble)
   integer s, i, j, k,ii,jj, js
   integer sitevarid, latvarid, lonvarid, datevarid, secvarid, htvarid, obsvarid, uncvarid, simvarid
   integer tnobs, nobs_srf, nobs_sat
   integer status, rcode
   integer count(2), dims2(2), start(2)
   integer fid, nobsid, en_size_id
   integer obs_flag, sample_flag 
   !obs_flag,1: large spatial scale representation; 0: local/regional influence
   !sample_flag: 1=4-hour avg; 2=1-hour avg; 3=90-min avg; 4=instantaneous
   integer, external :: dmin
   real hyai(nlev+1), hybi(nlev+1)
   real p1(nlev+1), pp(nlev), p0
   character(len=5) xclass(2)
   character(len=15) class1, class2
   real latk(ny), lonk(nx)
   real, external :: dxyp,emin, emax, stdev
   data xclass/'surf','xco2'/
   p0=1.0e5
  
   call read_ocean_map(glmap, nx, ny)

   do i=1, ny
     latk(i)=i-90.5
   enddo
   do i=1, nx
     lonk(i)=i-0.5
   enddo

 
   allocate(co2(nlon,nlat,nlev,nt,ensemble))
   allocate(Zm(nlon,nlat,nlev,nt))
   allocate(psm(nlon,nlat,nt))
   allocate(bxco2(bnx,bny))
   allocate(bnxco2(bnx,bny))
   bxco2=0
   bnxco2=0

   do i=1, nt
      secs(i)=(i-1)*10800+5400
   enddo

   write(sday,'(i8)') icdate

   sdate=icdate

   do i=1, nlat
      L_model(i)=sqrt(dxyp(xlat(i),real(resm)))
      if(L_model(i)<100) L_model(i)=100
   enddo
   is=0
   no=0
   nobs_srf=0
   nobs_sat=0
   do iday=1, rundays
      filesim=trim(mozartdir)//'/run_ensemble/hist/'//trim(casename)//'.mz4.h0.'//sday(1:4)//'-'//sday(5:6)//'-'//sday(7:8)//'-10800.nc'
      call readmzt(filesim,co2,Zm,psm,hyai,hybi,nlon,nlat,nlev,nt,ensemble)

! read surface CO2 data
    if(issrf==1)then  
      fileobs=trim(obsdir_srf)//'/srf_'//sday
      inquire(file=trim(fileobs),exist=alive)
      if(alive)then
         open(1,file=trim(fileobs))
         do while (.true.)
            read(1,*,end=100) site_tmp, date_tmp, sec_tmp, lat_tmp, lon_tmp,hgt_tmp, obs_tmp, unc_tmp    ! obs_flag, sample_flag 
            if(lon_tmp<0)then
               lon_tmp=lon_tmp+360
            endif
            ix_tmp=dmin2(lon_tmp,xlon,nlon)
            jx_tmp=dmin2(lat_tmp,xlat,nlat)
            tx=dmin2(sec_tmp,secs,nt)
           
!            ic=index(site_tmp,'_')+1
!            ec=index(site_tmp,'-')-1
!            if(ec<ic)then
!               ec=index(site_tmp(ic:),'-')-1
!               ec=ec+ic-1
!            endif
!
!            class1=site_tmp(ic:ec)  !tow or surface or aircraft
!            ic=ec+index(site_tmp(ec:),'-')
!            ec=index(site_tmp(ic:),'_')
!            if(ec<=0)then
!                 ec=len(site_tmp)
!            else
!                 ec=ec+ic-2
!            endif
!            class2=site_tmp(ic:ec)  !flask or instu

!            if(class2(1:4)=='insi') then
!               if(mod(int(sec_tmp),1800)/=0) cycle
!            endif
!
!            if(trim(class2)=='flask')then
!               unc_obs=0.1
!            elseif(trim(class2)=='insitu'.or.trim(class2)=='insit')then
!               unc_obs=0.5
!            elseif(trim(class2)=='pfp')then
!               unc_obs=0.1   !pfp ??
!            endif
!            if(obs_flag==1)then
!               L_obs=large_scale
!            else
!               L_obs=local_scale
!            endif
!            unc_mis=unc_obs*sqrt(L_model(jx_tmp)/L_obs)*1.5
            
!            unc_tmp=sqrt(unc_obs*unc_obs+unc_mis*unc_mis)

!            if(unc_tmp<0.5) unc_tmp=0.5
      
            kx_tmp=dmin2(hgt_tmp,Zm(ix_tmp,jx_tmp,:,tx),nlev)
!            if(sample_flag==1)then
!               if(tx==1)then
!                  Hx_tmp(1:ensemble)=(co2(ix_tmp,jx_tmp,kx_tmp,tx,1:ensemble)+co2(ix_tmp,jx_tmp,kx_tmp,tx+1,1:ensemble))/2.0
!               elseif(tx==8)then
!                  Hx_tmp(1:ensemble)=(co2(ix_tmp,jx_tmp,kx_tmp,tx,1:ensemble)+co2(ix_tmp,jx_tmp,kx_tmp,tx-1,1:ensemble))/2.0
!               else
!                  Hx_tmp(1:ensemble)=(co2(ix_tmp,jx_tmp,kx_tmp,tx,1:ensemble)+co2(ix_tmp,jx_tmp,kx_tmp,tx+1,1:ensemble)  &
!                                      +co2(ix_tmp,jx_tmp,kx_tmp,tx-1,1:ensemble))/3.0
!               endif
!            else
               Hx_tmp(1:ensemble)=co2(ix_tmp,jx_tmp,kx_tmp,tx,1:ensemble)
!            endif
               Hx_tmp(1:ensemble)=Hx_tmp(1:ensemble)*29/44*1.0e6
              sim_tmp=sum(Hx_tmp(1:ensemble))/ensemble
            if(abs(sim_tmp-obs_tmp)>15) cycle
            !sample_flag: 1=4-hour avg; 2=1-hour avg; 3=90-min avg;
            !4=instantaneous
            
!            write(*,'(a4,a10,i10,2f8.2,8f8.2)') site_tmp(1:3),site_tmp(ic:ec),date_tmp, sec_tmp, lat_tmp,L_model(jx_tmp), L_obs, unc_obs, unc_mis, unc_tmp,&
!                       obs_tmp,sum(Hx_tmp(1:ensemble))/ensemble, obs_tmp-sum(Hx_tmp(1:ensemble))/ensemble  
            if(is==0)then
               is=is+1
               nobs_srf=nobs_srf+1
               err1=1/unc_tmp/unc_tmp
               Robs(is)=err1
           !    Robs(is)=unc_tmp
               ix(is)=ix_tmp
               jx(is)=jx_tmp
               kx(is)=kx_tmp
               yobs(is)=obs_tmp*err1
               Hx(is,1:ensemble)=Hx_tmp(1:ensemble)*err1
               ylat(is)=lat_tmp
               ylon(is)=lon_tmp
               yhgt(is)=hgt_tmp
               ytype(is)=1
               no(is)=1
            else
               find=0
               do s=1, is
                  if(ix_tmp==ix(s).and.jx_tmp==jx(s).and.kx_tmp==kx(s))then
                     find=1
                     no(s)=no(s)+1
                     err1=1/unc_tmp/unc_tmp
                     yobs(s)=yobs(s)+obs_tmp*err1
                     Robs(s)=Robs(s)+err1
                     Hx(s,1:ensemble)=Hx(s,1:ensemble)+Hx_tmp(1:ensemble)*err1
                     yhgt(s)=yhgt(s)+hgt_tmp
                     exit
             !        ylat(s)=ylat(s)+lat_tmp
             !        if(abs(lon_tmp-ylon(s))>350)then
             !             if(ylon(s)>360)then
             !              lon_tmp=lon_tmp+360
             !             else
             !              lon_tmp=lon_tmp-360
             !             endif 
             !        endif
             !        ylon(s)=ylon(s)+lon_tmp
                  endif
               enddo
               if(find==0)then
                  is=is+1
                  nobs_srf=nobs_srf+1
                  err1=1/unc_tmp/unc_tmp
                  Robs(is)=err1
                  ix(is)=ix_tmp
                  jx(is)=jx_tmp
                  kx(is)=kx_tmp
                  yobs(is)=obs_tmp*err1
                  Hx(is,1:ensemble)=Hx_tmp(1:ensemble)*err1
                  ylat(is)=lat_tmp
                  ylon(is)=lon_tmp
                  yhgt(is)=hgt_tmp
                  ytype(is)=1
                  no(is)=1
               endif
            endif
         enddo
100 continue
         close(1)
      endif
    endif
      if(issat==1)then   ! assimilation satellite data
         file_gosat=trim(obsdir_gosat)//'/acos_v90_'//sday//'.nc'
         file_oco2=trim(obsdir_oco2)//'/oco2_v11_'//sday//'.nc'
         file_tansat=trim(obsdir_tansat)//'/TanSat_v2.'//sday//'.nc'
!real xco2(nobs_max),pxco2(nobs_max), xco2unc(nobs_max), ps(nobs_max)
!real p(nlevg,nobs_max),aj(nlevg,nobs_max),hj(nlevg,nobs_max),Yaj(nlevg,nobs_max)
!real latg(nobs_max),long(nobs_max),timeg(nobs_max)
    !     print*, trim(file_gosat)
         nobs=0
         inquire(file=trim(file_gosat),exist=alive)
         if(alive .and. isgosat==1)then
            call readsat(file_gosat,latg_tmp,long_tmp,timeg_tmp,xco2_tmp,aj_tmp,pxco2_tmp,hj_tmp,Yaj_tmp,p_tmp,ps_tmp,xco2unc_tmp,nlevg,nobs_max,glev_tmp,nobs_tmp)
            latg(nobs+1:nobs+nobs_tmp)=latg_tmp(1:nobs_tmp)
            long(nobs+1:nobs+nobs_tmp)=long_tmp(1:nobs_tmp)
            timeg(nobs+1:nobs+nobs_tmp)=timeg_tmp(1:nobs_tmp)
            xco2(nobs+1:nobs+nobs_tmp)=xco2_tmp(1:nobs_tmp)*1.00079-0.142
            pxco2(nobs+1:nobs+nobs_tmp)=pxco2_tmp(1:nobs_tmp)
            xco2unc(nobs+1:nobs+nobs_tmp)=xco2unc_tmp(1:nobs_tmp)
            ps(nobs+1:nobs+nobs_tmp)=ps_tmp(1:nobs_tmp)
            p(1:glev_tmp,nobs+1:nobs+nobs_tmp)=p_tmp(1:glev_tmp,1:nobs_tmp)
            aj(1:glev_tmp,nobs+1:nobs+nobs_tmp)=aj_tmp(1:glev_tmp,1:nobs_tmp)
            hj(1:glev_tmp,nobs+1:nobs+nobs_tmp)=hj_tmp(1:glev_tmp,1:nobs_tmp)
            Yaj(1:glev_tmp,nobs+1:nobs+nobs_tmp)=Yaj_tmp(1:glev_tmp,1:nobs_tmp)
            glev(nobs+1:nobs+nobs_tmp)=glev_tmp
            nobs=nobs+nobs_tmp
         endif
        
         inquire(file=trim(file_oco2),exist=alive)
         if(alive .and. isoco2==1)then
             call readsat(file_oco2,latg_tmp,long_tmp,timeg_tmp,xco2_tmp,aj_tmp,pxco2_tmp,hj_tmp,Yaj_tmp,p_tmp,ps_tmp,xco2unc_tmp,nlevg,nobs_max,glev_tmp,nobs_tmp)

            do js=1, nobs_tmp
                ix_tmp=dmin2(long_tmp(js),lonk,nx)
                jx_tmp=dmin2(latg_tmp(js),latk,ny)
                if(glmap(ix_tmp,jx_tmp)==0)then
                    xco2unc_tmp(js)=xco2unc_tmp(js)+0.2  !in ocean,err=err+2.0
                endif
            enddo
            latg(nobs+1:nobs+nobs_tmp)=latg_tmp(1:nobs_tmp)
            long(nobs+1:nobs+nobs_tmp)=long_tmp(1:nobs_tmp)
            timeg(nobs+1:nobs+nobs_tmp)=timeg_tmp(1:nobs_tmp)
            xco2(nobs+1:nobs+nobs_tmp)=xco2_tmp(1:nobs_tmp)*1.00079-0.142 !scaled to WMO 2019
            pxco2(nobs+1:nobs+nobs_tmp)=pxco2_tmp(1:nobs_tmp)
            xco2unc(nobs+1:nobs+nobs_tmp)=xco2unc_tmp(1:nobs_tmp)*3.5
            ps(nobs+1:nobs+nobs_tmp)=ps_tmp(1:nobs_tmp)
            p(1:glev_tmp,nobs+1:nobs+nobs_tmp)=p_tmp(1:glev_tmp,1:nobs_tmp)
            aj(1:glev_tmp,nobs+1:nobs+nobs_tmp)=aj_tmp(1:glev_tmp,1:nobs_tmp)
            hj(1:glev_tmp,nobs+1:nobs+nobs_tmp)=hj_tmp(1:glev_tmp,1:nobs_tmp)
            Yaj(1:glev_tmp,nobs+1:nobs+nobs_tmp)=Yaj_tmp(1:glev_tmp,1:nobs_tmp)
            glev(nobs+1:nobs+nobs_tmp)=glev_tmp
            nobs=nobs+nobs_tmp
         endif
  
         inquire(file=trim(file_tansat),exist=alive)
         if(alive .and. istansat==1)then
            call readsat(file_tansat,latg_tmp,long_tmp,timeg_tmp,xco2_tmp,aj_tmp,pxco2_tmp,hj_tmp,Yaj_tmp,p_tmp,ps_tmp,xco2unc_tmp,nlevg,nobs_max,glev_tmp,nobs_tmp)
            latg(nobs+1:nobs+nobs_tmp)=latg_tmp(1:nobs_tmp)
            long(nobs+1:nobs+nobs_tmp)=long_tmp(1:nobs_tmp)
            timeg(nobs+1:nobs+nobs_tmp)=timeg_tmp(1:nobs_tmp)
            xco2(nobs+1:nobs+nobs_tmp)=xco2_tmp(1:nobs_tmp)
            pxco2(nobs+1:nobs+nobs_tmp)=pxco2_tmp(1:nobs_tmp)
            xco2unc(nobs+1:nobs+nobs_tmp)=xco2unc_tmp(1:nobs_tmp)
            ps(nobs+1:nobs+nobs_tmp)=ps_tmp(1:nobs_tmp)
            p(1:glev_tmp,nobs+1:nobs+nobs_tmp)=p_tmp(1:glev_tmp,1:nobs_tmp)
            aj(1:glev_tmp,nobs+1:nobs+nobs_tmp)=aj_tmp(1:glev_tmp,1:nobs_tmp)
            hj(1:glev_tmp,nobs+1:nobs+nobs_tmp)=hj_tmp(1:glev_tmp,1:nobs_tmp)
            Yaj(1:glev_tmp,nobs+1:nobs+nobs_tmp)=Yaj_tmp(1:glev_tmp,1:nobs_tmp)
            glev(nobs+1:nobs+nobs_tmp)=glev_tmp
            nobs=nobs+nobs_tmp

         endif
!       print*, 'nobs=', nobs
              
         if(nobs>0)then
!            call readsat(fileobs,latg,long,timeg,xco2,aj,pxco2,hj,Yaj,p,ps,xco2unc,nlevg,nobs_max,glev,nobs)
           ! print*,sday, ' obs number in satellite=', nobs
            do i=1, nobs
               lat_tmp=latg(i)
               lon_tmp=long(i)
               if(lon_tmp<0)then
                  lon_tmp=lon_tmp+360     !lon: 0-359
               endif
               hgt_tmp=-999
               obs_tmp=xco2(i)
               unc_tmp=xco2unc(i)    !*1.9
            !   if(unc_tmp<1.0) unc_tmp=1.0
               
               ix_tmp=dmin2(long(i),xlon,nlon)
               jx_tmp=dmin2(latg(i),xlat,nlat)
               kx_tmp=-1
               tx=dmin2(timeg(i),secs,nt)
               p1=hyai*p0+hybi*psm(ix_tmp,jx_tmp,tx)
               pp=(p1(1:nlev)+p1(2:nlev+1))/2
               pp=pp/100
               call interplot(Ax(1:glev(i),1:ensemble),p(1:glev(i),i),pp,    &
               co2(ix_tmp,jx_tmp,1:nlev,tx,1:ensemble)*29/44*1.0e6, glev(i), nlev, ensemble) 
               !if(i==1.and.is==0)then
               !   print*, pp
               !   write(*,*) p(1:glev,i)
               !   do k=1, nlev
               !     write(*,'(i4,50f8.2)') k, co2(ix_tmp,jx_tmp,k,tx,1:ensemble)*29/44*1.0e6
               !   enddo
               !   do k=1, glev
               !     write(*,'(i4,50f8.2)')k, Ax(k,1:ensemble)
               !   enddo
               !endif 
               Hx_tmp(1:ensemble)=pxco2(i)
               do k=1, glev(i)
                !  if(i==1.and.is==0)then
                !    print*, k, Hx_tmp(1), Ax(k,1)-Yaj(k,i),hj(k,i),aj(k,i)
                !  endif
                  Hx_tmp(1:ensemble)=Hx_tmp(1:ensemble)+hj(k,i)*aj(k,i)*(Ax(k,1:ensemble)-Yaj(k,i))
               enddo

               sim_tmp=sum(Hx_tmp(1:ensemble))/ensemble
          
              if(abs(sim_tmp-obs_tmp)>15) cycle
               
               ii=dmin2(long(i),blon,bnx)
               jj=dmin2(latg(i),blat,bny)
               bxco2(ii,jj)=bxco2(ii,jj)+(sim_tmp-obs_tmp)
               bnxco2(ii,jj)=bnxco2(ii,jj)+1
 
               if(is==0)then
                  is=is+1
                  nobs_sat=nobs_sat+1
                  err1=1/unc_tmp/unc_tmp
                  Robs(is)=err1
                  !Robs(is)=unc_tmp
                  ix(is)=ix_tmp
                  jx(is)=jx_tmp
                  kx(is)=kx_tmp
                  yobs(is)=obs_tmp*err1
                  !yobs(is)=obs_tmp
                  Hx(is,1:ensemble)=Hx_tmp(1:ensemble)*err1
                  !Hx(is,1:ensemble)=Hx_tmp(1:ensemble)
                  ylat(is)=lat_tmp
                  ylon(is)=lon_tmp
                  yhgt(is)=hgt_tmp
                  ytype(is)=2
                  no(is)=1
               else
                  find=0
                  do s=1, is
                     if(ix_tmp==ix(s).and.jx_tmp==jx(s).and.kx_tmp==kx(s))then
                        find=1
                        no(s)=no(s)+1
                        err1=1/unc_tmp/unc_tmp
                        yobs(s)=yobs(s)+obs_tmp*err1
                        !yobs(s)=yobs(s)+obs_tmp
                        Robs(s)=Robs(s)+err1
                        !Robs(s)=Robs(s)+unc_tmp
                        Hx(s,1:ensemble)=Hx(s,1:ensemble)+Hx_tmp(1:ensemble)*err1
                        !Hx(s,1:ensemble)=Hx(s,1:ensemble)+Hx_tmp(1:ensemble)
                        yhgt(s)=yhgt(s)+hgt_tmp
                        exit
              !          ylat(s)=ylat(s)+lat_tmp
              !          if(abs(lon_tmp-ylon(s))>350)then
              !             if(ylon(s)>360)then
              !              lon_tmp=lon_tmp+360
              !             else
              !              lon_tmp=lon_tmp-360
              !             endif 
              !          endif
              !          ylon(s)=ylon(s)+lon_tmp
                     endif
                  enddo
                  if(find==0)then
                     is=is+1
                     nobs_sat=nobs_sat+1
                     err1=1/unc_tmp/unc_tmp
                     Robs(is)=err1
                     !Robs(is)=unc_tmp
                     ix(is)=ix_tmp
                     jx(is)=jx_tmp
                     kx(is)=kx_tmp
                     yobs(is)=obs_tmp*err1
                     Hx(is,1:ensemble)=Hx_tmp(1:ensemble)*err1
                    ! yobs(is)=obs_tmp
                    ! Hx(is,1:ensemble)=Hx_tmp(1:ensemble)

                     ylat(is)=lat_tmp
                     ylon(is)=lon_tmp
                     yhgt(is)=hgt_tmp
                     ytype(is)=2
                     no(is)=1
                  endif
               endif

         enddo
         endif
      endif
      call nextday(sday)
   enddo
   NumOfObs=is

   do is=1, NumOfObs
      yobs(is)=yobs(is)/Robs(is)
      Hx(is,1:ensemble)=Hx(is,1:ensemble)/Robs(is)
      Robs(is)=sqrt(1/Robs(is))
     ! yobs(is)=yobs(is)/no(is)
     ! Hx(is,1:ensemble)=Hx(is,1:ensemble)/no(is)
     ! Robs(is)=Robs(is)/no(is)
      Robs(is)=Robs(is)*1.9
      if(Robs(is)<1.0) Robs(is)=1.0
      yhgt(is)=yhgt(is)/no(is)
   !   ylat(is)=ylat(is)/no(is)
   !   ylon(is)=ylon(is)/no(is)
   enddo

   do ii=1, bnx
       do jj=1, bny
           if(bnxco2(ii,jj)>0)then
              bxco2(ii,jj)=bxco2(ii,jj)/bnxco2(ii,jj)
           endif
       enddo
   enddo

  filename =trim(opt_dir)//'/'//trim(casename)//'/emissions.opt.lambda.nc'
  call writebxco2(filename,bxco2,bnxco2,bnx,bny,da_step)

  write(*,'(a,i6,a,i6,a,i6)') 'total obs number is ', NumofObs, ' surface site is ', nobs_srf,' sat obs is', nobs_sat
!  NumofObs=1

  write(sday,'(i8)')icdate
  fileout = trim(opt_dir)//'/'//trim(casename)//'/Obs_sim_'//sday//'.txt'
  open(123,file=trim(fileout))
  abias_srf=0
  abias_sat=0 
  do is=1, NumofObs
    write(123,'(2i4,2x, a5,i4,10f8.2)') is, nlev-kx(is)+1, xclass(ytype(is)),no(is), ylat(is), ylon(is), yhgt(is), yobs(is), Robs(is),&
                           sum(Hx(is,1:ensemble))/ensemble,sum(Hx(is,1:ensemble))/ensemble-yobs(is),  &
                           emin(Hx(is,1:ensemble),ensemble),emax(Hx(is,1:ensemble),ensemble), &
                           emax(Hx(is,1:ensemble),ensemble)-emin(Hx(is,1:ensemble),ensemble)
    if(ytype(is)==1)then
         abias_srf=abias_srf+sum(Hx(is,1:ensemble))/ensemble-yobs(is)
    elseif(ytype(is)==2)then
       !  print*,sum(Hx(is,1:ensemble))/ensemble, yobs(is)
         abias_sat=abias_sat+sum(Hx(is,1:ensemble))/ensemble-yobs(is)
    endif 
  enddo 
  close(123)

  if(nobs_srf>0)then
  abias_srf=abias_srf/nobs_srf
  print*, 'The mean bias of surface obserations is ', abias_srf
  endif
  if(nobs_sat>0)then
  abias_sat=abias_sat/nobs_sat 
  print*, 'The mean bias of satellite obserations is ', abias_sat
  endif

end subroutine

function emin(a,n)
  real emin
  integer n
  real a(n)
  integer i
  emin=a(1)
  do i=2, n
     if(a(i)<emin)emin=a(i)
  enddo
  return
end function

function emax(a,n)
  real emax
  integer n
  real a(n)
  integer i
  emax=a(1)
  do i=2, n
     if(a(i)>emax)emax=a(i)
  enddo
  return
end function

function stdev(a, n)
 implicit none
 integer n, i
 real a(n)
 real stdev
 real aa
 
 aa=0
 do i=1, n
    aa=aa+a(i)
 enddo
 aa=aa/n

 stdev=0
 do i=1, n
    stdev=stdev+(a(i)-aa)*(a(i)-aa)
 enddo
 stdev=sqrt(stdev/n)
 return
end function

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


