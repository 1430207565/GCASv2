PROGRAM MAIN
  use module_global
  use ran_mod
  implicit none
  integer, parameter  :: nx=360, ny=180
  real lat(ny), lon(nx)
  integer    ios
  INTEGER    i,j,k,numofen,domain,is,s, ix(nx), jx(ny)

  REAL, ALLOCATABLE, DIMENSION(:,:,:)      ::      bio, mbio, ocn, fossil, fire, co2, co2_1
  REAL, ALLOCATABLE, DIMENSION(:,:,:)      ::      bio_n, ocn_n, fossil_n, fire_n
!mbio save multiyearmean daily (8 hours) bio flux
  REAL, ALLOCATABLE, DIMENSION(:,:,:)      ::      bgmean, bgmean_org, bgstdev
  REAL, ALLOCATABLE, DIMENSION(:,:,:,:)    ::      bgens, bgensA
  REAL, ALLOCATABLE, DIMENSION(:,:,:,:)    ::      lambda
  character(len=100)    FILE_NAME, filename
  INTEGER(4)                                       ::      FLAG
  integer iv, it, nt, its, ite
  real flux, fac
  character(len=8) sday
  real, external :: dxyp, stdev      

  !******************************
! ************************************************************
  
  call ensrf_init() 


  if(nrt==1)then
     write(*,'(a,i3,a,f4.1,a)'),'Start ', ensemble, ' ensembles DA for net flux with std of ',bunc(1),' for each grid'
  else 
  if(opt_scheme==1)then
     write(*,'(a,i3,a,f4.1,a)'),'Start ', ensemble, ' ensembles DA for Bio only with std of ',bunc(1),' for each grid'
  elseif(opt_scheme==2)then
    if(netflux == 0)then
     write(*,'(a,i3,a,f4.1,a,f4.1,a)'),'Start', ensemble, ' ensembles DA for Bio and Ocn with std of ',bunc(1),' and ', bunc(2),' for each grid'
    else
      write(*,'(a,i3,a,f4.1,a,f4.1,a)'),'Start', ensemble, ' ensembles DA for land net flux and Ocn with std of ',bunc(1),' and ', bunc(2),' for each grid'
    endif
  elseif(opt_scheme==3)then
     write(*,'(a,i3,a,f4.1,a,f4.1,a,f4.1,a)'),'Start', ensemble, ' ensembles DA for Bio, Ocn and Fire with std of ',  &
             bunc(1),' ,', bunc(2),' and ',bunc(3),' for each grid'
  endif
  endif
  nt=8*rundays+1

  do i=1, ny
    lat(i)=i-90.5
  enddo
  do i=1, nx
    lon(i)=i-0.5
  enddo
  
 !****************************************************
 !           generate lambda
 !****************************************************
  ALLOCATE(lambda(bnvar,bnx,bny,ensemble))

   do j=1, bny
    do i=1, bnx
      do iv=1, bnvar
        do is=1, ensemble
        ! if(blat(j)>=23.5.and.iv==1.and.month>=5.and.month<=8.and.nrt/=1)then
        !   lambda(iv,i,j,is)=normal(0.0, 1.0)*3.0
        !   
        ! else
           lambda(iv,i,j,is)=normal(0.0, 1.0)*bunc(iv)
        ! endif
        enddo
      enddo
    enddo
  enddo
!    if(blat(j)>23.5)then   !Northern 
!        if(opt_scheme==1.or.opt_scheme==2)then
!          if(month>=11.and.month<=12)then
!            bunc(1)=10
!          elseif(month>=1.and.month<=2)then
!            bunc(1)=10
!          endif
!        endif
!    elseif(blat(j)<-23.5)then  !Southern 
!        if(opt_scheme==1.or.opt_scheme==2)then
!          if(month>=5.and.month<=8)then
!            bunc(1)=10
!          endif
!        endif
!    endif

 !****************************************************
 !           open flux
 !****************************************************
  ALLOCATE(bio(nx,ny,nt))
  allocate(mbio(nx,ny,nt))
  ALLOCATE(ocn(nx,ny,nt))
  ALLOCATE(fossil(nx,ny,nt))
  ALLOCATE(fire(nx,ny,nt))
  allocate(co2(nx,ny,nt))
  allocate(co2_1(nlon,nlat,nt))

  ALLOCATE(bio_n(nx,ny,8))
  ALLOCATE(ocn_n(nx,ny,8))
  ALLOCATE(fossil_n(nx,ny,8))
  ALLOCATE(fire_n(nx,ny,8))

  write(sday,'(i4,i2.2,i2.2)') year, month, day
  do it=1, rundays
  !    OPEN EMISSION
     its=(it-1)*8+1
     ite=it*8
  !    print*, sday
     call readprior(priordir,bio(:,:,its:ite),ocn(:,:,its:ite),   &
                    fossil(:,:,its:ite),fire(:,:,its:ite),nx,ny,8,sday, netflux)   
     
     if(it==1)then
          call readavebio(mbio(:,:,its:ite),nx,ny,8)
     else
          mbio(:,:,its:ite)=mbio(:,:,1:8)
     endif
     
     call nextday(sday)
  enddo
  print*, 'read prior flux ok!'
 ! print*, sday
  call readprior(priordir,bio_n,ocn_n,fossil_n,fire_n,nx,ny,8,sday, netflux)

  mbio(:,:,nt)=mbio(:,:,1) 
  bio(:,:,nt)=bio_n(:,:,1)
  ocn(:,:,nt)=ocn_n(:,:,1)
  fossil(:,:,nt)=fossil_n(:,:,1)
  fire(:,:,nt)=fire_n(:,:,1)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ajusted by lgy
  ! FOR OSSE, scale the flux.
  ! ! OSSE_1
  ! ocn = ocn * 0.5 
  ! bio = bio * 0.6
  ! fire = fire * 0.7
  ! ! OSSE_2
  ! ocn = ocn * 1.5 
  ! bio = bio * 1.4
  ! fire = fire * 1.3
  ! ! OSSE_3
  ! ocn = ocn * 0.5 
  ! bio = bio * 0.6
  ! fire = fire * 1.3
  ! ! OSSE_4
  ! ocn = ocn * 1.5
  ! bio = bio * 1.4
  ! fire = fire * 0.7
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  mbio=bio


 !*******************idea case************************
!  bio=0
!  ocn=0
!  fossil=1.0e15
!  fire=0
!******************************************************
  allocate(bgens(bnx,bny,bnvar,ensemble))
  allocate(bgensA(nx,ny,bnvar,ensemble))
  allocate(bgmean(bnx,bny,bnvar))
  allocate(bgstdev(bnx,bny,bnvar))
  allocate(bgmean_org(bnx,bny,bnvar))
  bgens=0
  bgensA=0
  bgmean=0
  bgmean_org=0
  do i=1, nx
     ix(i)=dmin2(lon(i),blon,bnx)
  enddo
  do j=1, ny 
     jx(j)=dmin2(lat(j),blat,bny)
  enddo
 ! print*, ix
 ! print*, jx
  do it=1, 8*rundays

    do i=1, nx
      do j=1, ny
        if(nrt==1)then
          flux=bio(i,j,it)+ocn(i,j,it)+fossil(i,j,it)+fire(i,j,it)   
        !! molecules/cm2/s ==> gC/grid/3hrs
!**************idea case***********************
!          flux=1.0e15    ! for L_scale test
!**********************************************
          flux=flux/6.02e23*12*1.0e4*3600*3*dxyp(lat(j),1.0)*1.0e6
          bgmean_org(ix(i),jx(j),ibio)=bgmean_org(ix(i),jx(j),ibio)+flux
        else
          do iv=1, bnvar
              if(iv==1)then
                flux=bio(i,j,it)
              elseif(iv==2)then
                flux=ocn(i,j,it)
              elseif(iv==3)then
                flux=fire(i,j,it)
              endif
        !! molecules/cm2/s ==> gC/grid/3hrs
        !      print*, i, j, iv, flux
              flux=flux/6.02e23*12*1.0e4*3600*3*dxyp(lat(j),1.0)*1.0e6
        !      print*, i, j, iv, flux
              bgmean_org(ix(i),jx(j),iv)=bgmean_org(ix(i),jx(j),iv)+flux
          enddo
        endif
      ENDDO
    enddo

  enddo
 ! print*, bgmean_org(:,:,1)

  DO numofen=1,ensemble
    DO j=1,ny
      DO i=1,nx
        if(nrt==1)then
          co2(i,j,:)=(bio(i,j,:)+mbio(i,j,:))*0.5+ocn(i,j,:)+fossil(i,j,:)+fire(i,j,:)
          co2(i,j,:)=co2(i,j,:)+co2(i,j,:)*lambda(ibio,ix(i),jx(j),numofen)
        else
          if(opt_scheme==1)then
             co2(i,j,:)=bio(i,j,:)+(bio(i,j,:)+mbio(i,j,:))*0.5*lambda(ibio,ix(i),jx(j),numofen) + &
                  ocn(i,j,:)+ fossil(i,j,:)+fire(i,j,:)
          elseif(opt_scheme==2)then
             co2(i,j,:)=bio(i,j,:)+(bio(i,j,:)+mbio(i,j,:))*0.5*lambda(ibio,ix(i),jx(j),numofen) + &
                  ocn(i,j,:)+ocn(i,j,:)*lambda(iocn,ix(i),jx(j),numofen) + &
                fossil(i,j,:)+fire(i,j,:)
            ! if(i==120.and.j==120)then
            ! write(*,'(2i4,5e12.4)') i, j, bio(i,j,3), ocn(i,j,3), fossil(i,j,3), fire(i,j,3), co2(i,j,3)
            ! endif
          elseif(opt_scheme==3)then
             co2(i,j,:)=bio(i,j,:)+(bio(i,j,:)+mbio(i,j,:))*0.5*lambda(ibio,ix(i),jx(j),numofen) +  &
                  ocn(i,j,:)+ocn(i,j,:)*lambda(iocn,ix(i),jx(j),numofen) + &
                fire(i,j,:)+fire(i,j,:)*lambda(ifire,ix(i),jx(j),numofen)+  &
                fossil(i,j,:)
          endif
        endif
      ENDDO
    ENDDO
    do it=1, 8*rundays
      do i=1, nx
        do j=1, ny
          if(nrt==1)then
            flux=(bio(i,j,it)+mbio(i,j,it))*0.5+ocn(i,j,it)+fossil(i,j,it)+fire(i,j,it)
            flux=flux+flux*lambda(ibio,ix(i),jx(j),numofen)   
          !! molecules/cm2/s ==> gC/grid/3hrs
            flux=flux/6.02e23*12*1.0e4*3600*3*dxyp(lat(j),1.0)*1.0e6
            bgens(ix(i),jx(j),ibio,numofen)=bgens(ix(i),jx(j),ibio,numofen)+flux
            bgensA(i,j,ibio,numofen)=bgensA(i,j,ibio,numofen)+flux
          else
            do iv=1, bnvar
               if(iv==1)then
               flux=bio(i,j,it)+(bio(i,j,it)+mbio(i,j,it))*0.5*lambda(iv,ix(i),jx(j),numofen)   
               elseif(iv==2)then
               flux=ocn(i,j,it)+ocn(i,j,it)*lambda(iv,ix(i),jx(j),numofen)
               elseif(iv==3)then
               flux=fire(i,j,it)+fire(i,j,it)*lambda(iv,ix(i),jx(j),numofen)
               endif
            !! molecules/cm2/s ==> gC/grid/3hrs
               flux=flux/6.02e23*12*1.0e4*3600*3*dxyp(lat(j),1.0)*1.0e6
               bgens(ix(i),jx(j),iv,numofen)=bgens(ix(i),jx(j),iv,numofen)+flux
               bgensA(i,j,iv,numofen)=bgensA(i,j,iv,numofen)+flux
            enddo
          endif
        ENDDO
      enddo
    enddo
    write(FILE_NAME,'(a,i3.3,a)') trim(mozartdir)//'/data/emis/emissions.EN',numofen,'.surface.T42LR.nc'

    call spatial_interp(lat, lon, xlat, xlon, nx, ny, nlon, nlat, co2, co2_1, nt)
    call writeflux_formzt(FILE_NAME,year,month,day,rundays,co2_1,xlat,xlon,nlon,nlat)
    !    print*, trim(FILE_NAME) 
 end do
  do i=1, bnx
    do j=1, bny
       do iv=1, bnvar
           bgmean(i,j,iv)=sum(bgens(i,j,iv,1:ensemble))/ensemble
           bgstdev(i,j,iv)=stdev(bgens(i,j,iv,1:ensemble),ensemble)
       enddo
    enddo
  enddo

  !convert gC/grid/rundays => molecules/cm2/s
  !flux=flux/6.02e23*12*1.0e4*3600*3*dxyp(lat(j),1.0)*1.0e6
  !do i=1, bnx
  !  do j=1, bny
  !         fac=6.02e23/12/dxyp(blat(j),resb*1.0)/1.0e6/1.0e4/rundays/24/3600
   !        print*, blat(j), dxyp(blat(j),resb*1.0), fac
  !         bgmean(i,j,:)=bgmean(i,j,:)*fac
  !         bgens(i,j,:,:)=bgens(i,j,:,:)*fac
  !         bgmean_org(i,j,:)=bgmean_org(i,j,:)*fac
  !  enddo
  !enddo

  write(sday,'(i4,i2.2,i2.2)') year, month, day
  FILE_NAME=trim(mozartdir)//'/data/emis/emissions.intermediate.'//sday//'.nc'
  call writeflux_forassim(FILE_NAME,bgens,bgmean,bgmean_org,bgstdev,blat,blon,bnx,bny,ensemble,bnvar,rundays)

  FILE_NAME=trim(opt_dir)//'/'//trim(casename)//'/uncertainties.prior.'//sday//'.nc'
  call write_uncertainty(FILE_NAME,bgensA,lat,lon,nx,ny,bnvar,ensemble)

 ! print*, trim(FILE_NAME)
end program

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



!*****************************
