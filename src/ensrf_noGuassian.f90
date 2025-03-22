subroutine ensrf_noGuassian(yobs,Robs,HX,ylat,ylon,kx,yhgt,ytype,NumOfObs,DATA1,XMEAN,XMEAN_old,XNEW)

  USE module_global
  IMPLICIT NONE
  include 'netcdf.inc'
  REAL, PARAMETER                              :: eps  = 1.e-6
  INTEGER                                     ::      iofobs,NumOfObs
  integer, parameter :: nx=360, ny=180
  REAL                                               ::      BZTMP
  REAL                                               ::      tmpr1
  REAL                                               ::      dividend
  REAL                                               ::      localization
  REAL                                               ::      inflation = 1.1
  REAL                                               ::      inflation2 = 0.5
  REAL                                               ::      distance
  REAL                                               ::      kalman_gain, kalman_gain_per, kalman_adjust
  REAL,external                       ::      fcorrelation
  REAL                                ::      Robs2
  REAL                                ::      test_correlation
  REAL, DIMENSION(ensemble)           ::      testx1,testx2
  LOGICAL                             ::      FLAG
  INTEGER                             ::      i,j,k,nn, selected
  INTEGER                             ::      we,sn,ap, is
  REAL                                ::      tmp_r1, corrx
  character(len=100)                                  ::      filename
  REAL,dimension(maxobs)              ::      yobs,Robs,ylat,ylon,yhgt
  integer, dimension(maxobs)          ::      ytype, kx
  real                                ::      HX(maxobs,ensemble)
  REAL                                ::      HXMEAN(NumOfObs)
  REAL                                ::      HXPER(NumOfObs,ensemble),HZPER(NumOfObs,ensemble),HZTMPT(ensemble,1)
  REAL,DIMENSION(bnx,bny,bnvar,ensemble)    ::      DATA1,XPER,XNEW,ZPER
  REAL,DIMENSION(bnx,bny,bnvar)             ::     XMEAN,XMEAN_old, sscal
  integer, dimension(bnx,bny,bnvar)         :: snum
  integer, dimension(bnx,bny,bnvar,50)      :: sid
  real, dimension(bnx,bny,bnvar,50)          :: sdis 
 character(len=8) sday
  INTEGER                            ::       KSTART(4),KCOUNT(4)
  real bias
  real corr_t_test, LscaleMax
  integer, external :: izero
  real, external :: corr
  !******************************


!  call read_transcom_map(glmap, nx, ny)
  !1-11, Land, 12-22, Ocean

!  if(ensemble==50)then
!     corr_t_test=0.2352
!  elseif(ensemble==80)then
!     corr_t_test=0.1852
!  elseif(ensemble==100)then
!     corr_t_test=0.1654
!  elseif(ensemble==150)then
!     corr_t_test=0.1348
!  else
!     corr_t_test=0.15
!  endif
 

  write(sday,'(i8)')icdate

  dividend=sqrt(real(ensemble-1))
  !         CALL en_to_per(HX,HXMEAN,HZ,NumOfObs)
  !HX
  HXMEAN=0.
  do iofobs=1,NumOfOBS
    do nn=1,ensemble
        HXMEAN(iofobs) = HXMEAN(iofobs)+HX(iofobs,nn)
    enddo
    HXMEAN(iofobs) = HXMEAN(iofobs)/REAL(ensemble)
    do nn=1,ensemble
        HXPER(iofobs,nn)=HX(iofobs,nn)-HXMEAN(iofobs)
    enddo
  enddo

  !*********X
  XMEAN=0.
  do k=1,bnvar
    do j = 1,bny
      do i = 1,bnx
        do nn=1,ensemble
            XMEAN(i,j,k) = XMEAN(i,j,k)+DATA1(i,j,k,nn)
        enddo

        XMEAN(i,j,k) = XMEAN(i,j,k)/REAL(ensemble)
        do nn=1,ensemble
            XPER(i,j,k,nn)=DATA1(i,j,k,nn)-XMEAN(i,j,k)
            !                 ZPER(i,j,k,nn)=XPER(i,j,k,nn)/dividend
        enddo

      enddo
    enddo
  enddo
!calc_correl(xlat,xlon,ylat,ylon,ss,snum,Lsc,x,nx,ny,nv,ens,y,nobs,resx)
  call calc_correl(blat,blon,ylat,ylon,sid,sdis,snum,sscal,xper,bnx,bny,bnvar,ensemble,hxper,NumOfOBS,resb)

  XMEAN_old=XMEAN
  XNEW=XPER

  DO ap = 1,bnvar
   
   if(ap==1.or.ap==3)then
          LscaleMax=3000
   else
          LscaleMax=3000
   endif

    DO sn = 1,bny
      LOOP_M:   DO we = 1,bnx
        
        if(snum(we,sn,ap)>0)then
        
        LOOP_OBS:   DO is = 1,snum(we,sn,ap)

            iofobs=sid(we,sn,ap,is)
            distance=sdis(we,sn,ap,is)

          tmp_r1 = distance / LscaleMax     !sscal(we,sn,ap)     !L_scale(ytype(iofobs))

          IF(tmp_r1 .gt. 1) cycle LOOP_OBS
          localization = fcorrelation(tmp_r1)

          !*****************************************************************
          !         HZTMPT(1:ensemble,1)=HZPER(iofobs,1:ensemble)
          !         HZTMPT(1:ensemble,1)=HXPER(iofobs,1:ensemble)/dividend

          BZTMP=0.
          !!***************************
          !         DO i = 1,ensemble
          !            BZTMP=BZTMP + HZTMPT(i,1)*HZTMPT(i,1)
          !         ENDDO
          !********
          DO i = 1,ensemble
            BZTMP=BZTMP + HXPER(iofobs,i)*HXPER(iofobs,i)    !HPBHT
          ENDDO
          BZTMP=BZTMP/real(ensemble-1)
          !Pb background error covariance

        !  fac(we,sn,ap,iofobs,1)=BZTMP    !HPbH in the formula

          Robs2=Robs(iofobs)*Robs(iofobs)
          !********************************
          !

          tmpr1=sqrt(Robs2/(BZTMP+Robs2))
          kalman_adjust=1/(1+tmpr1)



          kalman_gain=0.
          DO i = 1,ensemble
            !               kalman_gain = kalman_gain + ZPER(we,sn,ap,i)*HZTMPT(i,1)
            !               kalman_gain = kalman_gain + XPER(we,sn,ap,i)*HZTMPT(i,1)/dividend
            kalman_gain = kalman_gain + XPER(we,sn,ap,i)*HXPER(iofobs,i)/real(ensemble-1)
          ENDDO

         ! fac(we,sn,ap,iofobs,2)=kalman_gain

          !without localization
          !            kalman_gain=kalman_gain/(BZTMP+R)
          !with localization
          if(ytype(iofobs)==2)then
             if(ap==2) then  !  .or.ap==3)then
!                localization=localization*0.75
!             else
                localization=localization*1.5
             endif
          endif

          kalman_gain=kalman_gain*localization/(BZTMP+Robs2)
          IF(kalman_gain .lt. 0) cycle LOOP_OBS
          kalman_gain_per=kalman_gain*kalman_adjust
          !*******************
          !fac(we,sn,ap,iofobs,3)=kalman_gain   !K in the formula

          !IF(DEBUG) THEN
         !   WRITE(logun,*) "***********************************"
           !  if(blat(sn)==27.and.blon(we)==273.and.ap==1)then
           !     WRITE(logun,'(2f8.2,i4,3e10.2)') ylat(iofobs), ylon(iofobs), kx(iofobs),fac(we,sn,ap,iofobs,1:3)
         !   WRITE(logun,*) Robs2, BZTMP, kalman_gain
           !  endif
         ! ENDIF


          XMEAN(we,sn,ap)=XMEAN(we,sn,ap)+kalman_gain*(yobs(iofobs)-HXMEAN(iofobs))
          DO i = 1,ensemble
            XPER(we,sn,ap,i)=XPER(we,sn,ap,i)-kalman_gain_per*(HXPER(iofobs,i))
            ! XPER(we,sn,ap,i)=XPER(we,sn,ap,i)-inflation*kalman_gain_per*(HXPER(iofobs,i))
            XNEW(we,sn,ap,i)=XMEAN(we,sn,ap)+XPER(we,sn,ap,i)
          ENDDO

        ENDDO LOOP_OBS
       endif
      ENDDO LOOP_M
    ENDDO
  ENDDO

!  do i=1, ensemble
!      XNEW(:,:,:,i)=XMEAN_old(:,:,:)+XPER(:,:,:,i)
!  enddo
 
  
  filename =trim(opt_dir)//'/'//trim(casename)//'/emissions.opt.lambda.nc'
  call writenobs(filename,snum,sscal,bnx,bny,bnvar,da_step)

END subroutine


subroutine EarthSurfacedistance(lat1,lon1,lat2,lon2,dis)
  implicit NONE
  real lat1, lon1, lat2, lon2
  real dis
  real*8, parameter :: DEF_R =6378.137   !R of earth, km
  real*8, parameter :: DEF_PI180= 0.01745329252
  real*8, parameter :: DEF_2PI= 6.28318530712, DEF_PI = 3.14159265359
  real*8 fl, f, g, l, sg, sl, sf, s, c, w, r, d, h1, h2

  if(lat1==lat2.and.lon1==lon2)then
   dis=0
  else

   f = (lat1 + lat2)/2*DEF_PI180
   g = (lat1 - lat2)/2*DEF_PI180
   l = (lon1 - lon2)/2*DEF_PI180

   sg = sin(g)
   sl = sin(l)
   sf = sin(f)

   fl = 1/298.257

   sg = sg*sg
   sl = sl*sl
   sf = sf*sf

   s = sg*(1-sl) + (1-sf)*sl
   c = (1-sg)*(1-sl) + sf*sl

   w = atan(sqrt(s/c))
   r = sqrt(s*c)/w
   d = 2*w*DEF_R
   h1 = (3*r -1)/2/c
   h2 = (3*r +1)/2/s

   dis= d*(1 + fl*(h1*sf*(1-sg) - h2*(1-sf)*sg))
  endif

end subroutine


!*********************************************8

function izero(a,n)
  implicit none
  integer n, izero
  real a(n)
  integer i

  izero=0
  do i=1, n
     if(a(i)/=0)then
        izero=1
        return
     endif
  enddo
  return

end function

Function corr(s, o, n)
 implicit none
 integer     :: n
 real        :: s(n),o(n)
 real        :: corr
 integer     :: i
 real        :: rivl1,rivl2,rivl3
 real        :: critv(4), Sr, t_value
 data critv/1.676,1.664,1.660,1.645/
 !size, 50, 80, 100, 150 
 corr=0

! cal coef 
 rivl1=0
 rivl2=0
 rivl3=0
 do i=1,n
       rivl1=rivl1+s(i)*o(i)
       rivl2=rivl2+s(i)*s(i)
       rivl3=rivl3+o(i)*o(i)
 enddo

 rivl2=sqrt(rivl2)*sqrt(rivl3)
 if(rivl2>0)then
   corr=rivl1/rivl2
 endif

 !t-test
! if(corr>0)then
!    Sr=sqrt((1-corr*corr)/(n-2))
!    t_value=corr/Sr
!    if(n==50)then
!       if(t_value<1.676)then
!         corr=0
!       endif
!    elseif(n==80)then
!       if(t_value<1.664)then
!         corr=0
!       endif
!    elseif(n==100)then
!       if(t_value<1.66)then
!         corr=0
!       endif
!    elseif(n>100)then
!       if(t_value<1.645)then
!         corr=0
!       endif
!    endif
! endif 

 return
End Function


subroutine handle_err_ensrf(errcode)
  implicit none
  include 'netcdf.inc'
  integer errcode

  print *, 'Ensrf For Emission Error: ', nf_strerror(errcode)
  stop 2
end
!***************************

subroutine en_to_per(X,XMEAN,Z,M)


  USE module_global  ,  only   :          MISSING , ensemble

  IMPLICIT NONE
  INTEGER                         ::         M
  REAL          ,DIMENSION(M,ensemble)   ::         X,Z
  REAL          ,DIMENSION(M)     ::         XMEAN
  REAL                            ::         Xave,dividend
  INTEGER                         ::         i,j
  LOGICAL                         ::         FLAG
  !  ####################################
  dividend=sqrt(real(ensemble-1))
  !  ###################################
  do i=1,M

    Xave=0.0
    !
    FLAG = .false.
    do j=1,ensemble
      IF(X(i,j)==MISSING)THEN
        Xave = MISSING
        FLAG = .true.
      ELSE
        Xave=Xave+X(i,j)
      ENDIF
    enddo

    IF(FLAG)THEN
      Z(i,:) = MISSING
      XMEAN(i) = MISSING
    ELSE

      Xave=Xave/real(ensemble)
      XMEAN(i) = Xave
      !
      do j=1,ensemble
        Z(i,j)=(X(i,j)-Xave)/dividend !
      enddo
    ENDIF

  enddo

  return
end subroutine en_to_per


subroutine en_to_per_std(X,XMEAN,Z,M)


  USE module_global  ,  only   :          MISSING , ensemble

  IMPLICIT NONE
  INTEGER                         ::         M
  REAL          ,DIMENSION(M,ensemble)   ::         X,Z
  REAL          ,DIMENSION(M)     ::         XMEAN
  REAL                            ::         Xave,dividend,std
  INTEGER                         ::         i,j
  LOGICAL                         ::         FLAG

  !  ####################################
  dividend=sqrt(real(ensemble-1))
  !  ###################################
  do i=1,M

    Xave=0.0
    std =0.0
    !
    FLAG = .false.
    do j=1,ensemble
      IF(X(i,j)==MISSING)THEN
        Xave = MISSING
        FLAG = .true.
      ELSE
        Xave=Xave+X(i,j)
      ENDIF
    enddo

    IF(FLAG)THEN
      Z(i,:) = MISSING
      XMEAN(i) = MISSING
    ELSE


      Xave=Xave/real(ensemble)
      XMEAN(i) = Xave
      !
      do j=1,ensemble
        !        Z(i,j)=(X(i,j)-Xave)/dividend !
        std = std + (X(i,j)-Xave)*(X(i,j)-Xave)
      enddo

      do j = 1,ensemble
        Z(i,j)=(X(i,j)-Xave)/dividend !
      enddo

    ENDIF

  enddo

  return
end subroutine en_to_per_std



Function fcorrelation(z)
  USE module_global , only : MISSING
  IMPLICIT NONE
  !    Construction of Correlation in 2d or 3d
  !    Gregory Gaspari and Stephen E. Cohn
  REAL       ,PARAMETER       :: EPS  = 1.E5
  REAL                        :: z
  REAL                        :: fcorrelation
  REAL                        :: fff

  IF(z < 0.0)THEN
    fff = 0.
  ELSEIF(z < 1.0)THEN
    fff = -z**5/4. + z**4/2. + 5./8.*z**3 -5./3.*z*z + 1.
  ELSEIF(z < 2.0)THEN
    fff = 1./12. * z**5 - 1./2. *z**4 + 5./8.*z**3 + 5./3. *z**2 -5.*z +4. -2./3./z
  ELSE
    fff = 0.0
  ENDIF

!  if(fff>0)then
  fcorrelation = fff
!  else
!  fcorrelation = 0
!  endif
  !    IF(fff == 0.0) THEN
  !       CorOfConstrunction = MISSING
  !    ELSE
  !       CorOfConstrunction = 1./fff
  !    ENDIF
  !
  !    IF(CorOfConstrunction <= 1.E-4 )THEN
  !       CorOfConstrunction = MISSING
  !    ELSEIF(CorOfConstrunction > EPS)THEN
  !       CorOfConstrunction = MISSING
  !    ENDIF

END FUNCTION fcorrelation
