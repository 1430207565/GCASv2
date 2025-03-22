subroutine ensrf2emission(XMEAN,XMEAN_old,XNEW)
  use module_global
  implicit none
  integer, parameter :: nx=360, ny=180, nt=8
  real bio(nx,ny,nt), ocn(nx,ny,nt),fossil(nx,ny,nt),fire(nx,ny,nt), flux(nx,ny,nt)
  real xbio(nx,ny), xocn(nx,ny), xfossil(nx,ny), xfire(nx,ny)
  real opt_bio
  real r(nx,ny,bnvar), lat(ny), lon(nx), r_x(nx,ny,bnvar,ensemble)
  character(len=8) sday
  real, allocatable, dimension(:,:,:,:) :: co2
  real, allocatable, dimension(:,:,:) :: co2_1
  integer ii,jj, ix, jx, it, ntt, its, ite
  INTEGER                                            ::      i,j,sp,k,iofen

  REAL                        ,DIMENSION(bnx,bny,bnvar)  ::  XOPT, XORG, lambda,lambda_old
  REAL                        ,DIMENSION(bnx,bny,bnvar)  ::  XMEAN,XMEAN_old,dxmean
  REAL                                               ::      XNEW(bnx,bny,bnvar,ensemble), XENS(bnx,bny,bnvar,ensemble)
  real                                               :: XENSA(nx,ny,bnvar,ensemble), XENSB(nx,ny,bnvar,ensemble)
  real                                               ::      ratio_x(bnx,bny,bnvar,ensemble)
  CHARACTER(len=4)                                   ::      DIRN(ensemble),string_tmp
  character(len=100)                                  ::      filename
  INTEGER(4)                                         ::      FLAG
  REAL                                               ::      tmp
  real, external  :: globaltotal, dxyp
  real fac, ro, rn, f, maxlambda, minlambda
  character(len=10) varname

  do i=1, ny
    lat(i)=i-90.5
  enddo
  do i=1, nx
    lon(i)=i-0.5
  enddo

  ntt=8*rundays+1

  allocate(co2(nx,ny,ntt,2))
  allocate(co2_1(nlon,nlat,ntt))

  write(sday,'(i8)')icdate
  filename=trim(mozartdir)//'/data/emis/emissions.intermediate.'//sday//'.nc'
   
  call read_bgmean(filename,xorg,xens,bnx,bny,bnvar,ensemble)

  filename=trim(opt_dir)//'/'//trim(casename)//'/uncertainties.prior.'//sday//'.nc'

  call read_uncertainty(filename,xensA,nx,ny,bnvar,ensemble)

  !calculate ratio between posterior and prior ensemble
  do sp=1, bnvar
    do j=1, bny
      do i=1, bnx
        do iofen=1, ensemble
          if(xens(i,j,sp,iofen)/=0)then
              ratio_x(i,j,sp,iofen)=xnew(i,j,sp,iofen)/xens(i,j,sp,iofen)
          else
              ratio_x(i,j,sp,iofen)=0
          endif
        enddo
      enddo
    enddo
  enddo

  !calculate ratio between posterior and prior flux
  XMEAN_old = XMEAN - XMEAN_old
  do sp=1, bnvar
    DO j=1,bny
      DO i=1,bnx
        xopt(i,j,sp)=  xorg(i,j,sp)+XMEAN_old(i,j,sp)         !负值变成0 还有随机的平均不是零 会反演不到的地方 排放会变!
        if(xorg(i,j,sp)==0)then
           lambda(i,j,sp)=1
        else
           lambda(i,j,sp)=xopt(i,j,sp)/xorg(i,j,sp)
        endif
        if(sp==1)then
            maxlambda=10
            minlambda=0.1
        elseif(sp==2)then
            maxlambda=10
            minlambda=0.1
        elseif(sp==3)then
            maxlambda=5
            minlambda=0.2
        endif
          if(lambda(i,j,sp)>0)then
            if(lambda(i,j,sp)>maxlambda)then
              lambda(i,j,sp)=maxlambda
            elseif(lambda(i,j,sp)<minlambda)then
              lambda(i,j,sp)=minlambda
            endif
          elseif(lambda(i,j,sp)<0)then
!              lambda(i,j,sp)=minlambda
           if(sp/=3)then
            if(lambda(i,j,sp)<-1*maxlambda)then
              lambda(i,j,sp)=-1*maxlambda
            elseif(lambda(i,j,sp)>-1*minlambda)then
               lambda(i,j,sp)=-1*minlambda
            endif
           else
             lambda(i,j,sp)=minlambda
           endif
          endif
      enddo
    ENDDO
  ENDDO

!smth9
 
!  call linearsmth9(lambda,bnx,bny,bnvar)

!---------------------lambda--------------------------------------------!
! keep the lambda coherent, not change too much between t-1 and t windows
!
! filename =trim(opt_dir)//'/'//trim(casename)//'/emissions.opt.lambda.nc'
! if(da_step==1)then
!   call createlambda(filename, blat,blon, bnx, bny, bnvar)
! endif
 if(da_step>1)then
  
  filename =trim(opt_dir)//'/'//trim(casename)//'/emissions.opt.lambda.nc'
  call readlambda(filename, lambda_old, bnx, bny, bnvar, da_step-1)
 
  do sp=1, bnvar
     do j=1, bny
         do i=1, bnx
           ro=lambda_old(i,j,sp)
           rn=lambda(i,j,sp)
          ! rn=(rn+ro)/2
      !     lambda(i,j,sp)=(lambda(i,j,sp)+lambda_old(i,j,sp))/2
!use a continous asumption
!  |(lambda(t)-lambda(t-1))/lambda(t-1)|<1.5
          if(ro/=0)then
           f=(rn-ro)/ro
           if(f>1.5)then
             lambda(i,j,sp)=ro*2.5
           elseif(f<-1.5)then
             lambda(i,j,sp)=-0.5*ro
           endif
           if(sp==1)then
            maxlambda=10
            minlambda=0.1
           elseif(sp==2)then
            maxlambda=10
            minlambda=0.1
           elseif(sp==3)then
            maxlambda=5
            minlambda=0.2
           endif
           if(lambda(i,j,sp)>0)then
              if(lambda(i,j,sp)>maxlambda)then
                 lambda(i,j,sp)=maxlambda
              elseif(lambda(i,j,sp)<minlambda)then
                 lambda(i,j,sp)=minlambda
              endif
           elseif(lambda(i,j,sp)<0)then
            !     lambda(i,j,sp)=minlambda
             if(sp/=3)then 
              if(lambda(i,j,sp)<-1*maxlambda)then
                 lambda(i,j,sp)=-1*maxlambda
              elseif(lambda(i,j,sp)>-1*minlambda)then
                 lambda(i,j,sp)=-1*minlambda
              endif
             else
                lambda(i,j,sp)=minlambda
             endif
           endif
         endif
         enddo
     enddo
  enddo
 endif


 do sp=1, bnvar
   do j=1, bny
     do i=1, bnx       
       xopt(i,j,sp) = lambda(i,j,sp)*xorg(i,j,sp)
     enddo
   enddo
 enddo

  print*, 'Species     ORG      OPT'
  do sp=1, bnvar
     if(nrt==1)then
        varname='Netflux'
     else
        if(sp==1)then
          varname='Bio'
        elseif(sp==2)then
          varname='Ocn'
        elseif(sp==3)then
          varname='Fire'
        endif
     endif

    write(*,'(a7,2f8.2,a)') trim(varname),globaltotal(xorg(:,:,sp),bnx,bny)/1.0e12, globaltotal(xopt(:,:,sp),bnx,bny)/1.0e12, ' TgC'
  enddo

  filename = trim(opt_dir)//'/'//trim(casename)//'/emissions.opt.compare.'//sday//'.nc'
  call writeopt_tmp(filename,xopt,xorg,xopt-xorg,blat,blon,bnx,bny,bnvar)

  filename =trim(opt_dir)//'/'//trim(casename)//'/emissions.opt.lambda.nc'
  call writelambda(filename,lambda,xopt-xorg,bnx,bny,bnvar,sday, da_step)

  do i=1, nx
    do j=1, ny
      ix=dmin2(lon(i),blon,bnx)
      jx=dmin2(lat(j),blat,bny)
      r(i,j,:)=lambda(ix,jx,:)
      r_x(i,j,:,:)=ratio_x(ix,jx,:,:)
    enddo
  enddo
!  xbio=0
!  xocn=0
!  xfossil=0
!  xfire=0
  write(sday,'(i4,i2.2,i2.2)') year, month, day
  do it=1, rundays+1
  !    OPEN EMISSION
  !    print*, sday
      call readprior(priordir,bio,ocn,fossil,fire,nx,ny,nt,sday,netflux)
  !update 2023.9.8
      if(it==rundays+1)then   !update 2023.9.8
         co2(:,:,ntt,2)=bio(:,:,1)+ocn(:,:,1)+fossil(:,:,1)+fire(:,:,1)  !update 2023.9.8
      else  !update 2023.9.8
         its=(it-1)*nt+1
         ite=it*nt
         co2(:,:,its:ite,2)=bio(:,:,1:nt)+ocn(:,:,1:nt)+fossil(:,:,1:nt)+fire(:,:,1:nt)
      endif    !!update 2023.9.8
      if(nrt==1)then
         flux=bio+ocn+fossil+fire
      endif


!*******************New Scheme 2021.6.7 jiangf***********************
      
     if(nrt==1)then
         call fluxadjust(flux,r(:,:,1),nx,ny,nt)
     else
        if(opt_scheme==1)then
          call fluxadjust(bio,r(:,:,ibio),nx,ny,nt)
        elseif(opt_scheme==2)then
          call fluxadjust(bio,r(:,:,ibio),nx,ny,nt)
          call fluxadjust(ocn,r(:,:,iocn),nx,ny,nt)
        elseif(opt_scheme==3)then
          call fluxadjust(bio,r(:,:,ibio),nx,ny,nt)
          call fluxadjust(ocn,r(:,:,iocn),nx,ny,nt)
          call fluxadjust(fire,r(:,:,ifire),nx,ny,nt)
        endif
     endif
!******************Old Scheme**************************************** 
!      do i=1, nt
!        if(nrt==1)then
!          flux(:,:,i)=flux(:,:,i)*r(:,:,1)
!        else
!          if(opt_scheme==1)then
!            bio(:,:,i)=bio(:,:,i)*r(:,:,ibio)
!          elseif(opt_scheme==2)then
!            bio(:,:,i)=bio(:,:,i)*r(:,:,ibio)
!            ocn(:,:,i)=ocn(:,:,i)*r(:,:,iocn)
!          elseif(opt_scheme==3)then
!            bio(:,:,i)=bio(:,:,i)*r(:,:,ibio)
!            ocn(:,:,i)=ocn(:,:,i)*r(:,:,iocn)
!            fossil(:,:,i)=fossil(:,:,i)*r(:,:,ifossil)
!          endif
!        endif
!      enddo
!*********************************************************************

!     print*,'ddddddd'
   if(it==rundays+1)then  !update 2023.9.8
      if(nrt==1)then      !update 2023.9.8
        co2(:,:,ntt,1)=flux(:,:,1)    !update 2023.9.8
     else                              !update 2023.9.8
        co2(:,:,ntt,1)=bio(:,:,1)+ocn(:,:,1)+fossil(:,:,1)+fire(:,:,1) !1, optimized, 2, orgined  !update 2023.9.8
     endif     !update 2023.9.8
   else    !update 2023.9.8
     if(nrt==1)then
        co2(:,:,its:ite,1)=flux(:,:,1:nt)    
     else
        co2(:,:,its:ite,1)=bio(:,:,1:nt)+ocn(:,:,1:nt)+fossil(:,:,1:nt)+fire(:,:,1:nt)  !1, optimized, 2, orgined
     endif
      filename = trim(opt_dir)//'/'//trim(casename)//'/emissions.opt.'//sday//'.nc'

!     print*, filename
     if(nrt==1)then
      call writeopt_nrt(filename,flux,nx,ny,nt,sday)
      else
      call writeopt(filename,bio,ocn,fossil,fire,nx,ny,nt,sday)
     endif
  endif  !update 2023.9.8
     call nextday(sday)
  enddo

!  co2(:,:,ntt,1)=co2(:,:,ntt,2)   !update 2023.9.12
  
  !write flux for runing mozart4 again
  call spatial_interp(lat, lon, xlat, xlon, nx, ny, nlon, nlat, co2(:,:,:,1), co2_1, ntt)
  !write optimized flux
  filename=trim(mozartdir)//'/data/emis/emissions.EN001.surface.T42LR.nc'
  call writeflux_formzt(filename,year,month,day,rundays,co2_1,xlat,xlon,nlon,nlat)

  call spatial_interp(lat, lon, xlat, xlon, nx, ny, nlon, nlat, co2(:,:,:,2),co2_1, ntt)
 !write prior flux
  filename=trim(mozartdir)//'/data/emis/emissions.EN002.surface.T42LR.nc'
  call writeflux_formzt(filename,year,month,day,rundays,co2_1,xlat,xlon,nlon,nlat)

  do i=1, nx
    do j=1, ny
        xensB(i,j,:,:)=xensA(i,j,:,:)*r_x(i,j,:,:)
    enddo
  enddo


  write(sday,'(i8)')icdate
  filename=trim(opt_dir)//'/'//trim(casename)//'/uncertainties.posterior.'//sday//'.nc'
  call write_uncertainty(filename,xensB,lat,lon,nx,ny,bnvar,ensemble)

  deallocate( co2 )

end subroutine

subroutine fluxadjust(f,l,nx,ny,nt)
  implicit none
  integer nx, ny, nt
  real f(nx,ny,nt), l(nx,ny)
  real t(nx,ny), abst(nx,ny)
  real x, delta
  integer i, j, it
  t=0
  abst=0
  do i=1, nx
     do j=1, ny
        do it=1, nt
            t(i,j)=t(i,j)+f(i,j,it)
            abst(i,j)=abst(i,j)+abs(f(i,j,it))
        enddo
     enddo
  enddo

  do i=1, nx
    do j=1, ny
      delta=t(i,j)-l(i,j)*t(i,j)
      if(delta/=0)then
          x=abst(i,j)/delta
          if(x/=0)then
            do it=1, nt
               f(i,j,it)=f(i,j,it)-abs(f(i,j,it))/x
            enddo
          endif
      endif
    enddo
  enddo

end subroutine


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

function globaltotal(a,nx,ny)
  implicit none
  integer nx,ny
  real a(nx,ny)
  real globaltotal
  integer i, j
  
  globaltotal=0
  do i=1, nx
    do j=1, ny
      globaltotal=globaltotal+a(i,j)
    enddo
  enddo

end function

subroutine linearsmth9(a,nx,ny,nv)
 implicit none
 integer nx, ny, nv
 real a(nx,ny,nv), b(nx,ny,nv)
 integer i, j, i1, i2, j1, j2
 b=a
 do i=1, nx
   do j=2, ny-1
     i1=i+1
     j1=j+1
     i2=i-1
     j2=j-1
     if(i1>nx) i1=1
     if(i2<1) i2=nx
     a(i,j,:)=(b(i,j,:)+b(i1,j,:)+b(i2,j,:)+b(i,j1,:)+b(i,j2,:)+b(i1,j1,:)+b(i1,j2,:)+b(i2,j1,:)+b(i2,j2,:))/9.0
   enddo
 enddo


end subroutine


