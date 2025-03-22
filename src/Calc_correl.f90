!calc_correl(blat,blon,ylat,ylon,yhgt,xper,bnx,bny,bnvar,ensemble,hxper,NumOfOBS)
subroutine calc_correl(xlat,xlon,ylat,ylon,ss,ssh,snum,Lsc,x,nx,ny,nv,ens,y,nobs,resx)
 USE module_global, only: opt_dir, nlev, casename
 use qsort_c_module
 implicit none
 real undef 
 integer resx, nx, ny, nv, ens, nobs
 real x(nx,ny,nv,ens), y(nobs,ens), c(nobs), h(nobs)
 real res
 integer id(nobs), id1(nobs)
 integer ss(nx,ny,nv,50), snum(nx,ny,nv)  !selected siteid
 real    ssh(nx,ny,nv,50), Lsc(nx,ny,nv)
 real xlat(ny), xlon(nx), ylat(nobs),ylon(nobs),yhgt(nobs)
 real lc(14,nobs), lcm(14,nobs)
 integer nc(14,nobs)
 real, external :: correl
 integer i, j, iv, is,k,kk
 character(len=5) xclass(2)
 real l,corr_t_test, L_scale, LscaleMax
 integer kold
 real, external :: dxyp
 integer, external :: izero
 data xclass/'surf','xco2'/

 if(ens==50)then
     corr_t_test=0.2352
  elseif(ens==80)then
     corr_t_test=0.1852
  elseif(ens==100)then
     corr_t_test=0.1654
  elseif(ens==150)then
     corr_t_test=0.1348
  else
     corr_t_test=0.15
  endif

 do i=1, nobs
   id(i)=i
 enddo

 undef=-999
 c=0
 do i=1, nx
      do j=1, ny
         res=sqrt(dxyp(xlat(j),real(resx)))
         if(res<100) res=100
      !   print*, res, xlat(j), resx
         do iv=1, nv
           if(iv==1.or.iv==3)then
                LscaleMax=3000
           else
                LscaleMax=3000
           endif
          ! print*, i,j, iv
           if(izero(x(i,j,iv,:),ens)==1)then
             do is=1, nobs
               c(is)=correl(x(i,j,iv,:),y(is,:),ens,undef)
               call EarthSurfacedistance(ylat(is),ylon(is),xlat(j),xlon(i),h(is))

             enddo
           !  print*, i, j, iv, 'start sort'
             id1=id
             call QsortC(c, id1)

            Lsc(i,j,iv)=0
            kk=0
            do is=1,nobs
            !     print*,is,c(is),id1(is), h(id1(is))
                 if(h(id1(is))<500.and.c(is)>0) then   ! corr_t_test)then
                    kk=kk+1
                    if(kk>20) then
                       exit
                    endif
                    ss(i,j,iv,kk)=id1(is)
                    ssh(i,j,iv,kk)=h(id1(is))
           !          print*, kk, c(is), id1(is), h(id1(is))
                 endif
            enddo
           if(kk>20)then
               k=kk
              Lsc(i,j,iv)=500
           else
            l=750
            L_scale=l
            kold=kk
            do while(l<=LscaleMax)
             k=kk
             
             do is=1,nobs
            !     print*,is,c(is),id1(is), h(id1(is))                 
                 if( h(id1(is))>=250 .and. h(id1(is))<=l .and. c(is)>=corr_t_test)then
                     k=k+1
                     if(k>20) then
                       exit
                     endif
                     ss(i,j,iv,k)=id1(is)
                     ssh(i,j,iv,k)=h(id1(is))
                     
                    ! write(*,'(f5.0,i3,f8.2,i4,f8.2,f6.0)'),l, k, c(is), id1(is), h(id1(is))
                 endif
             enddo
           !  print*, k, kold
             if(k/=kold) then
                L_scale=l
                kold=k
            !    print*, 'Not equal',kold, L_scale
             endif
             if(k>20) then
               exit
             endif
             !print*, L_scale
             l=l+250
            enddo
          endif
            if(k>0) then
            Lsc(i,j,iv)=L_scale
            else
            Lsc(i,j,iv)=0
            endif
           ! print*,i,j,iv, k, Lsc(i,j,iv) 
             snum(i,j,iv)=k
           !  pause
           endif
         enddo
     enddo
 enddo

end subroutine 



Function AVE(a, n, undef)
  implicit none
  integer i, n, m
  real a(n)
  real undef
  real ave

  ave=0
  m=0
  do i=1, n
    if(a(i)/=undef)then
          ave=ave+a(i)
          m=m+1
        endif
  enddo
  ave=ave/m

  return
End Function

Function CORREL(s, o, n, undef)
 implicit none
 integer     :: n
 real        :: s(n),o(n)
 real        :: correl
 integer     :: i
 real        :: rivl1,rivl2,rivl3,as,ao, undef
 integer     :: m
 real, external :: ave
 correl=0

 as=ave(s, n, undef)
 ao=ave(o, n, undef)
 
! cal coef 
 rivl1=0
 rivl2=0
 rivl3=0
 do i=1,n
        if(s(i)/=undef .and. o(i)/=undef)then
           rivl1=rivl1+(s(i)-as)*(o(i)-ao)
           rivl2=rivl2+(s(i)-as)*(s(i)-as)
           rivl3=rivl3+(o(i)-ao)*(o(i)-ao)
        end if
 enddo
 
 rivl2=sqrt(rivl2)*sqrt(rivl3)
 if(rivl2==0)then
   correl=0
 else
   correl=rivl1/rivl2
 endif
 return
End Function

