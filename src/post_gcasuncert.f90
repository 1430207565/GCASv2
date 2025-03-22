program  post_mozart_uncert
 use module_global
 implicit none
 integer, parameter :: nx=360, ny=180,nt=8
 integer, parameter :: nr=58
 integer nm, na
 real,allocatable, dimension(:,:,:,:) :: xold, xnew, xbio, xocn, xfoss, xfire, xflux
 real,allocatable, dimension(:,:,:) :: xold_r, xnew_r, stdxold, stdxnew
 real,allocatable, dimension(:,:) :: fossil_m, fossil_a, fire_m, fire_a
 real,allocatable, dimension(:,:) :: flux_m, flux_a, bio_m, bio_a, ocn_m, ocn_a
 real,allocatable, dimension(:,:) :: obio_m, obio_a, oocn_m, oocn_a
 real,allocatable, dimension(:,:) :: ofossil_m, ofossil_a, ofire_m, ofire_a
 real,allocatable, dimension(:,:) :: oflux_m, oflux_a
 real, allocatable, dimension(:,:)  :: stdxold_r, stdxnew_r
 integer, allocatable, dimension(:) :: nmonth, nyear, myear
 real lat(ny), lon(nx)
 integer b(nx,ny,nr)
 real area(nr), area1(nx,ny)
 real ur(nx,ny)
 integer sd, ed, ir,i, j, k, t, d,m,r, y, sm, em, im, iy, iv, iw
 character(len=100) filename1, filename2 ,outfile, regname(nr)
 integer mon(12), mon1(12), mon2(12)
 character(len=60) dir
 integer start_year, start_month, start_day, end_year, end_month, end_day
 real, external :: dxyp, stdev
 integer narg, tlen, ncycle
 CHARACTER(len=32) ::arg
 character(len=8) sday
 logical alive1

 data mon1/31,28,31,30,30,30,31,31,30,31,30,30/
 data mon2/31,29,31,30,31,30,31,31,30,31,30,30/

 data regname/'North_America_Boreal','North_America_Temperate','South_America_Tropical','South_America_Temperate', &
              'North_Africa','South_Africa',    &
              'Asia_Boreal','Asia_Temperate','Southeast_Asia','Austrialia','Europe',    &
              'North_Pacific','West_Pacific','East_Pacific','South_Pacific',  &
              'Arctic_Ocean','North_Atlantic','Tropical_Atlantic','South_Atlantic','Southern_Ocean', &
              'North_Indian_Ocean','South_Indian_Ocean',  &
              'China','North_America', 'South_America','Asia','Africa','Atlantic','Pacific',  &
              'Indian_Ocean',  &
              'Tropical_land','Tropical_ocean','Northern_Hemisphere_land','Northern_Hemisphere_ocean', &
              'Southern_Hemisphere_land','Southern_Hemisphere_Ocean','Global_land', &
              'Global_Ocean', 'Globe','USA','EU','India','Russia',  &
              'N._land_gt30N','T._land_30Nto30S','S._land_lt30S','Africa_india','Asia_pacific',  &
              'R2_N_America','R2_S_America','R2_Russia','R2_Europe','R2_WestAsia','R2_Africa','R2_EastAsia','R2_SouthAsia','R2_trop_Asia','R2_Aus'/


narg=IARGC()
if(narg<2)then
 print*,'the input parameter error!'
 stop
endif
if(narg >0)then
call getarg(1,arg)
read(arg(1:4),*)start_year
read(arg(5:6),*)start_month
read(arg(7:8),*)start_day
call getarg(2,arg)
read(arg,*)tlen
endif

call ensrf_init()

allocate(xold(nx,ny,bnvar,ensemble))
allocate(xnew(nx,ny,bnvar,ensemble))
allocate(stdxold(nx,ny,bnvar))
allocate(stdxnew(nx,ny,bnvar))

ncycle=tlen/rundays
tlen=rundays*ncycle
tlen=tlen-1
call the_end(start_year,start_month,start_day, end_year,end_month,end_day,tlen)

print*, start_year, start_month, start_day
print*, end_year, end_month, end_day

dir=trim(opt_dir)//'/'//trim(casename)
! calculate total month and year
na=end_year-start_year+1
nm=0
do i=start_year, end_year
  if(i==start_year)then
     sm=start_month
  else
     sm=1
  endif
  if(i==end_year)then
     em=end_month
  else
     em=12
  endif
  do j=sm, em
     nm=nm+1
  enddo
enddo

if(nrt==0)then
allocate(bio_m(nr,nm))
allocate(bio_a(nr,na))
allocate(ocn_m(nr,nm))
allocate(ocn_a(nr,na))
allocate(fossil_m(nr,nm))
allocate(fossil_a(nr,na))
allocate(fire_m(nr,nm))
allocate(fire_a(nr,na))

allocate(obio_m(nr,nm))
allocate(obio_a(nr,na))
allocate(oocn_m(nr,nm))
allocate(oocn_a(nr,na))
allocate(ofossil_m(nr,nm))
allocate(ofossil_a(nr,na))
allocate(ofire_m(nr,nm))
allocate(ofire_a(nr,na))

allocate(xbio(nx,ny,na,2))
allocate(xocn(nx,ny,na,2))
allocate(xfoss(nx,ny,na,2))
allocate(xfire(nx,ny,na,2))

endif

allocate(xold_r(nr,bnvar,ensemble))
allocate(xnew_r(nr,bnvar,ensemble))
allocate(stdxold_r(nr,bnvar))
allocate(stdxnew_r(nr,bnvar))

allocate(flux_m(nr,nm))
allocate(flux_a(nr,na))

allocate(oflux_m(nr,nm))
allocate(oflux_a(nr,na))

allocate(xflux(nx,ny,na,2))

allocate(nmonth(nm))
allocate(nyear(na))
allocate(myear(nm))

do i=1, ny
   lat(i)=-89.5+(i-1)
enddo
do i=1, nx
   lon(i)=0.5+(i-1)
enddo

!print*,'latlon end'

call read_map(b,nx,ny,nr)

!print*,'read_map end'

r=1
if(nrt==0)then
bio_m=0
bio_a=0
ocn_m=0
ocn_a=0
fossil_m=0
fossil_a=0
fire_m=0
fire_a=0

xbio=0
xocn=0
xfoss=0
xfire=0

obio_m=0
obio_a=0
oocn_m=0
oocn_a=0
ofossil_m=0
ofossil_a=0
ofire_m=0
ofire_a=0

endif
flux_m=0
flux_a=0

oflux_m=0
oflux_a=0

xflux=0

im=1
k=1
iy=1
iw=0

open(101,file=trim(dir)//'/'//'UR_for_eachstep.dat',form='unformatted',access='direct',recl=nx*ny)

do y=start_year,end_year
  if(start_year==end_year)then
    sm=start_month
    em=end_month
  else
   if(y==start_year)then
      sm=start_month
      em=12
   elseif(y==end_year)then
      sm=1
      em=end_month
   else
     sm=1
     em=12
   endif
  endif
   if(y==2012.or.y==2016.or.y==2000.or.y==2004.or.y==2008.or.y==2020)then
     mon=mon2
   else
     mon=mon1
   endif
    nyear(iy)=y
   do m=sm, em
     ! print*, im
      nmonth(im)=m
      myear(im)=y
    if(start_year==end_year.and.start_month==end_month)then
       sd=start_day
       ed=end_day
    else
      if(y==start_year.and.m==start_month)then
         sd=start_day
         ed=mon(m)
      elseif(y==end_year.and.m==end_month)then
         sd=1
         ed=end_day
      else
         sd=1
         ed=mon(m)
      endif
    endif
      do d=sd, ed
      write(sday,'(i4,2i2.2)') y, m, d

      
      write(filename1,'(a,i4,2i2.2,a3)')trim(dir)//'/uncertainties.prior.',y,m,d,'.nc'
      write(filename2,'(a,i4,2i2.2,a3)')trim(dir)//'/uncertainties.posterior.',y,m,d,'.nc'
!      print*, trim(filename1)
      inquire(file=filename1, exist=alive1) 
      if(alive1)then
        ! print*, 'uncertainties.prior.',y,m,d,'.nc'
         iw=iw+1
      !   write(*,'(i5,a,i4.4,i2.2,i2.2,a)') iw, 'uncertainties.prior.',y,m,d,'.nc'
         xold_r=0
         xnew_r=0
         call read_uncertainty(filename1,xold,nx,ny,bnvar,ensemble)
         call read_uncertainty(filename2,xnew,nx,ny,bnvar,ensemble)    !gC/grid/rundays
         xnew=xnew/1.0e15
         xold=xold/1.0e15                 !gC/grid/rundays--> PgC/grid/rundays
         do iv=1, bnvar
           do i=1, nx
             do j=1, ny
                  stdxold(i,j,iv)=stdev(xold(i,j,iv,:),ensemble)
                  stdxnew(i,j,iv)=stdev(xnew(i,j,iv,:),ensemble)
             enddo
           enddo
         enddo

         do i=1, nx
             do j=1, ny
              if(nrt==0)then
                if((stdxold(i,j,1)+stdxold(i,j,2))>0)then
                   ur(i,j)=(1-(stdxnew(i,j,1)+stdxnew(i,j,2))/(stdxold(i,j,1)+stdxold(i,j,2)))*100
                else
                   ur(i,j)=-999
                endif
              else
                if(stdxold(i,j,1)>0)then
                   ur(i,j)=(1-stdxnew(i,j,1))/stdxold(i,j,1)*100
                else
                   ur(i,j)=-999
                endif
              endif
             enddo
         enddo

         write(101,rec=iw) ur
 
 
         if(nrt==0)then

         xbio(:,:,iy,1)=xbio(:,:,iy,1)+stdxold(:,:,1)
         xbio(:,:,iy,2)=xbio(:,:,iy,2)+stdxnew(:,:,1)

         xocn(:,:,iy,1)=xocn(:,:,iy,1)+stdxold(:,:,2)
         xocn(:,:,iy,2)=xocn(:,:,iy,2)+stdxnew(:,:,2)

         if(bnvar>=3)then
            xfoss(:,:,iy,1)=xfoss(:,:,iy,1)+stdxold(:,:,3)
            xfoss(:,:,iy,2)=xfoss(:,:,iy,2)+stdxnew(:,:,3)
         endif
         if(bnvar>=4)then
            xfire(:,:,iy,1)=xfire(:,:,iy,1)+stdxold(:,:,4)
            xfire(:,:,iy,2)=xfire(:,:,iy,2)+stdxnew(:,:,4)
         endif
         
         else
            xflux(:,:,iy,1)=xflux(:,:,iy,1)+stdxold(:,:,1)
            xflux(:,:,iy,2)=xflux(:,:,iy,2)+stdxnew(:,:,1)
         endif
 
         do ir=1, nr
            do i=1, nx
               do j=1, ny
                  xold_r(ir,:,:)=xold_r(ir,:,:)+xold(i,j,:,:)*b(i,j,ir)
                  xnew_r(ir,:,:)=xnew_r(ir,:,:)+xnew(i,j,:,:)*b(i,j,ir)
               enddo
            enddo

            do iv=1, bnvar
               stdxold_r(ir,iv)=stdev(xold_r(ir,iv,:),ensemble)
               stdxnew_r(ir,iv)=stdev(xnew_r(ir,iv,:),ensemble)
            enddo
         enddo

        do ir=1, nr
          do iv=1, bnvar
            flux_m(ir,im)=flux_m(ir,im)+stdxnew_r(ir,iv)*stdxnew_r(ir,iv)
            oflux_m(ir,im)=oflux_m(ir,im)+stdxold_r(ir,iv)*stdxold_r(ir,iv)
            flux_a(ir,iy)=flux_a(ir,iy)+stdxnew_r(ir,iv)*stdxnew_r(ir,iv)
            oflux_a(ir,iy)=oflux_a(ir,iy)+stdxold_r(ir,iv)*stdxold_r(ir,iv)
          enddo
        enddo
 

         if(nrt==0)then
         do ir=1, nr
            bio_m(ir,im)=bio_m(ir,im)+stdxnew_r(ir,1)*stdxnew_r(ir,1)
            obio_m(ir,im)=obio_m(ir,im)+stdxold_r(ir,1)*stdxold_r(ir,1)
            bio_a(ir,iy)=bio_a(ir,iy)+stdxnew_r(ir,1)*stdxnew_r(ir,1)
            obio_a(ir,iy)=obio_a(ir,iy)+stdxold_r(ir,1)*stdxold_r(ir,1)
       !     print*, iy, im, bio_a(nr,iy,1:5)
            ocn_m(ir,im)=ocn_m(ir,im)+stdxnew_r(ir,2)*stdxnew_r(ir,2)
            oocn_m(ir,im)=oocn_m(ir,im)+stdxold_r(ir,2)*stdxold_r(ir,2)
            ocn_a(ir,iy)=ocn_a(ir,iy)+stdxnew_r(ir,2)*stdxnew_r(ir,2)
            oocn_a(ir,iy)=oocn_a(ir,iy)+stdxold_r(ir,2)*stdxold_r(ir,2)

            if(bnvar>=3)then
               fire_m(ir,im)=fire_m(ir,im)+stdxnew_r(ir,3)*stdxnew_r(ir,3)
               ofire_m(ir,im)=ofire_m(ir,im)+stdxold_r(ir,3)*stdxold_r(ir,3)
               fire_a(ir,iy)=fire_a(ir,iy)+stdxnew_r(ir,3)*stdxnew_r(ir,3)
               ofire_a(ir,iy)=ofire_a(ir,iy)+stdxold_r(ir,3)*stdxold_r(ir,3)
            endif
            if(bnvar>=4)then
               fossil_m(ir,im)=fossil_m(ir,im)+stdxnew_r(ir,4)*stdxnew_r(ir,4)
               ofossil_m(ir,im)=ofossil_m(ir,im)+stdxold_r(ir,4)*stdxold_r(ir,4)
               fossil_a(ir,iy)=fossil_a(ir,iy)+stdxnew_r(ir,4)*stdxnew_r(ir,4)
               ofossil_a(ir,iy)=ofossil_a(ir,iy)+stdxold_r(ir,4)*stdxold_r(ir,4)
            endif
          enddo
         endif

      endif
     
   

    enddo
    im=im+1
  enddo
  iy=iy+1
enddo

close(101)


if (nrt==0)then
 open(1,file=trim(dir)//'/'//'uncertanity.bio.dat',form='unformatted',access='direct',recl=nx*ny)
 do iy=1, na
     write(1,rec=(iy-1)*2+1) xbio(:,:,iy,1)
     write(1,rec=iy*2) xbio(:,:,iy,2)
 enddo
 close(1)
 
 open(1,file=trim(dir)//'/'//'uncertanity.ocn.dat',form='unformatted',access='direct',recl=nx*ny)
 do iy=1, na
     write(1,rec=(iy-1)*2+1) xocn(:,:,iy,1)
     write(1,rec=iy*2) xocn(:,:,iy,2)
 enddo
 close(1)

else
  open(1,file=trim(dir)//'/'//'uncertanity.netflux.dat',form='unformatted',access='direct',recl=nx*ny)
 do iy=1, na
     write(1,rec=(iy-1)*2+1) xflux(:,:,iy,1)
     write(1,rec=iy*2) xflux(:,:,iy,2)
 enddo
 close(1)
endif

do ir=1, nr
   do im=1, nm
     if(nrt==0)then
      obio_m(ir,im)=sqrt(obio_m(ir,im))
      bio_m(ir,im)=sqrt(bio_m(ir,im))
      oocn_m(ir,im)=sqrt(oocn_m(ir,im))
      ocn_m(ir,im)=sqrt(ocn_m(ir,im))
      ofossil_m(ir,im)=sqrt(ofossil_m(ir,im))
      fossil_m(ir,im)=sqrt(fossil_m(ir,im))
      ofire_m(ir,im)=sqrt(ofire_m(ir,im))
      fire_m(ir,im)=sqrt(fire_m(ir,im))
     endif
      oflux_m(ir,im)=sqrt(oflux_m(ir,im))
      flux_m(ir,im)=sqrt(flux_m(ir,im))
   enddo
   do im=1, na
     if(nrt==0)then
      obio_a(ir,im)=sqrt(obio_a(ir,im))
      bio_a(ir,im)=sqrt(bio_a(ir,im))
      oocn_a(ir,im)=sqrt(oocn_a(ir,im))
      ocn_a(ir,im)=sqrt(ocn_a(ir,im))
      ofossil_a(ir,im)=sqrt(ofossil_a(ir,im))
      fossil_a(ir,im)=sqrt(fossil_a(ir,im))
      ofire_a(ir,im)=sqrt(ofire_a(ir,im))
      fire_a(ir,im)=sqrt(fire_a(ir,im))
     endif
      oflux_a(ir,im)=sqrt(oflux_a(ir,im))
      flux_a(ir,im)=sqrt(flux_a(ir,im))

   enddo
enddo

!data regname/'North America Boreal','North America Temperate','South America
!Tropical','South America Temperate', &
!              'North Africa','South Africa',    &
!              'Asia Boreal','Asia Temperate','Southeast
!Asia','Austrialia','Europe',    &
!              'North Pacific','West Pacific','East Pacific','South Pacific',  &
!              'Arctic Ocean','North Atlantic','Tropical Atlantic','South
!Atlantic',' Southern Ocean', &
!              'North Indian Ocean','South Indian Ocean','ice
!sheet','Globe','Global land', &
!              'Global Ocean','Tropical land','Tropical ocean','Northern
!Hemisphere land','Northern Hemisphere ocean', &
!              'Southern Hemisphere land','Southern Hemisphere Ocean', 'China'/
if(nrt==1)then

write(*,'(a)') ' Year  Net_org Net_opt (PgC/a)'
do ir=1, nr
  open(1,file=trim(dir)//'/'//trim(regname(ir))//'_annualuncertainties.txt')
  write(1,'(a)') ' Year  Net_org Net_opt (PgC/a)'
  if(trim(regname(ir))=='Globe') then
         write(*,*) '--------------------Globe--------------------'
      endif
      if(trim(regname(ir))=='Global_land') then
         write(*,*) '---------------Global land--------------------'
      endif
      if(  trim(regname(ir))=='Global_Ocean')then
         write(*,*) '--------------Global Ocean--------------------'
      endif
      if(   trim(regname(ir))=='Europe')then
         write(*,*) '-------------------Europe--------------------'
      endif
      if(   trim(regname(ir))=='North_America_Temperate')then
         write(*,*) '-------------North America Temperate---------'
      endif
      if(   trim(regname(ir))=='China')then
         write(*,*) '--------------------China--------------------'
      endif
   do i=1, na
       write(1,'(i5,2f8.2)')nyear(i), oflux_a(ir,i), flux_a(ir,i)
       if(trim(regname(ir))=='Globe'.or.    &
            trim(regname(ir))=='Global_land'.or.  &
            trim(regname(ir))=='Global_Ocean'.or. &
            trim(regname(ir))=='Europe'.or.    &
            trim(regname(ir))=='North_America_Temperate' .or.   &
            trim(regname(ir))=='China')then
           write(*,'(i5,2f8.2)')nyear(i), oflux_a(ir,i), flux_a(ir,i)
       endif
   enddo
   close(1)
   
enddo
do ir=1, nr
  open(1,file=trim(dir)//'/'//trim(regname(ir))//'_monthlyuncertainties.txt')
  write(1,'(a)') '   ID Year  Mon Net_org Net_opt (PgC/a)'
  do i=1, nm
       write(1,'(3i5,2f8.2)')i,nyear(i),nmonth(i),oflux_m(ir,i), flux_m(ir,i)
   enddo
   close(1)
enddo


elseif(nrt==0)then


  if(opt_scheme==1)then
   write(*,'(a)') ' Year  Bio_org Bio_opt     Ocn  Fossil    Fire Net_org Net_opt (PgC/a)'
   elseif(opt_scheme==2)then
   write(*,'(a)') ' Year  Bio_org Bio_opt Ocn_org Ocn_opt  Fossil    Fire Net_org Net_opt (PgC/a)'
   elseif(opt_scheme==3)then
   write(*,'(a)') ' Year  Bio_org Bio_opt Ocn_org Ocn_opt Fir_org Fir_opt  Fossil Net_org Net_opt (PgC/a)'
   elseif(opt_scheme==4)then
   write(*,'(a)') ' Year  Bio_org Bio_opt Ocn_org Ocn_opt Fsl_org Fsl_opt Fir_org Fir_opt Net_org Net_opt (PgC/a)'
   endif


   do ir=1, nr
      open(1,file=trim(dir)//'/'//trim(regname(ir))//'_annualuncertainties.txt')
      if(opt_scheme==1)then
   write(1,'(a)') ' Year  Bio_org Bio_opt     Ocn  Fossil    Fire Net_org Net_opt (PgC/a)'
   elseif(opt_scheme==2)then
   write(1,'(a)') ' Year  Bio_org Bio_opt Ocn_org Ocn_opt  Fossil    Fire Net_org Net_opt (PgC/a)'
   elseif(opt_scheme==3)then
   write(1,'(a)') ' Year  Bio_org Bio_opt Ocn_org Ocn_opt Fir_org Fir_opt   Fossil Net_org Net_opt (PgC/a)'
   elseif(opt_scheme==4)then
   write(1,'(a)') ' Year  Bio_org Bio_opt Ocn_org Ocn_opt Fsl_org Fsl_opt Fir_org Fir_opt Net_org Net_opt (PgC/a)'
   endif

      if(trim(regname(ir))=='Globe') then
         write(*,*) '--------------------Globe--------------------'
      endif
      if(trim(regname(ir))=='Global_land') then
         write(*,*) '---------------Global land--------------------'
      endif
      if(  trim(regname(ir))=='Global_Ocean')then
         write(*,*) '--------------Global Ocean--------------------'
      endif
      if(   trim(regname(ir))=='Europe')then
         write(*,*) '-------------------Europe--------------------'
      endif
      if(   trim(regname(ir))=='North_America_Temperate')then
         write(*,*) '-------------North America Temperate---------'
      endif
      if(   trim(regname(ir))=='China')then
         write(*,*) '--------------------China--------------------'
      endif
      do i=1, na
     if(opt_scheme==1)then
        write(1,'(i5,7f8.2)')nyear(i),obio_a(ir,i),bio_a(ir,i),0,0,0,oflux_a(ir,i), flux_a(ir,i)
     elseif(opt_scheme==2)then
        write(1,'(i5,8f8.4)')nyear(i),obio_a(ir,i),bio_a(ir,i),oocn_a(ir,i),ocn_a(ir,i),0,0,oflux_a(ir,i), flux_a(ir,i)
     elseif(opt_scheme==3)then
          write(1,'(i5,9f8.2)')nyear(i),obio_a(ir,i),bio_a(ir,i),oocn_a(ir,i),ocn_a(ir,i),ofire_a(ir,i),  &
                             fire_a(ir,i),0,oflux_a(ir,i), flux_a(ir,i)
     elseif(opt_scheme==4)then
          write(1,'(i5,10f8.2)')nyear(i),obio_a(ir,i),bio_a(ir,i),oocn_a(ir,i),ocn_a(ir,i),ofossil_a(ir,i),  &
                                fossil_a(ir,i),ofire_a(ir,i),fire_a(ir,i),oflux_a(ir,i), flux_a(ir,i)
      endif

         if(trim(regname(ir))=='Globe'.or.    &
            trim(regname(ir))=='Global_land'.or.  &
            trim(regname(ir))=='Global_Ocean'.or. &
            trim(regname(ir))=='Europe'.or.    &
            trim(regname(ir))=='North_America_Temperate' .or.   &
            trim(regname(ir))=='China')then
 
     if(opt_scheme==1)then
        write(*,'(i5,7f8.2)')nyear(i),obio_a(ir,i),bio_a(ir,i),0,0,0,oflux_a(ir,i), flux_a(ir,i)
     elseif(opt_scheme==2)then
        write(*,'(i5,8f8.2)')nyear(i),obio_a(ir,i),bio_a(ir,i),oocn_a(ir,i),ocn_a(ir,i),0,0,oflux_a(ir,i),flux_a(ir,i)
     elseif(opt_scheme==3)then
        write(*,'(i5,9f8.2)')nyear(i),obio_a(ir,i),bio_a(ir,i),oocn_a(ir,i),ocn_a(ir,i),ofire_a(ir,i), &
                             fire_a(ir,i),0,oflux_a(ir,i), flux_a(ir,i)
     elseif(opt_scheme==4)then
        write(*,'(i5,10f8.2)')nyear(i),obio_a(ir,i),bio_a(ir,i),oocn_a(ir,i),ocn_a(ir,i),ofossil_a(ir,i), &
                                fossil_a(ir,i),ofire_a(ir,i),fire_a(ir,i),oflux_a(ir,i), flux_a(ir,i)
      endif



        endif
     enddo
     close(1)

   enddo

do ir=1, nr
   open(1,file=trim(dir)//'/'//trim(regname(ir))//'_monthlyuncertainties.txt')
   if(nrt==1)then
   write(1,'(a)') '   ID Year  Mon Net_org Net_opt (TgC/m)'
   else
   if(opt_scheme==1)then
   write(1,'(a)') '   ID Year  Mon Bio_org Bio_opt     Ocn  Fossil    Fire Net_org Net_opt (TgC/m)'
   elseif(opt_scheme==2)then
   write(1,'(a)') '   ID Year  Mon Bio_org Bio_opt Ocn_org Ocn_opt  Fossil    Fire Net_org Net_opt (TgC/m)'
   elseif(opt_scheme==3)then
   write(1,'(a)') '   ID Year  Mon Bio_org Bio_opt Ocn_org Ocn_opt Fir_org Fir_opt  Fossil Net_org Net_opt (TgC/m)'
   elseif(opt_scheme==4)then
   write(1,'(a)') '   ID Year  Mon Bio_org Bio_opt Ocn_org Ocn_opt Fsl_org Fsl_opt Fir_org Fir_opt Net_org Net_opt (TgC/m)'
   endif

   endif

   do i=1, nm
     if(nrt==1)then
         write(1,'(3i5,7f8.2)')i,myear(i),nmonth(i),oflux_m(ir,i), flux_m(ir,i)
     else

     if(opt_scheme==1)then
        write(1,'(3i5,7f8.2)')i,myear(i),nmonth(i), obio_m(ir,i),bio_m(ir,i),0,0,0,oflux_m(ir,i), flux_m(ir,i)
     elseif(opt_scheme==2)then
        write(1,'(3i5,8f8.4)')i,myear(i),nmonth(i),obio_m(ir,i),bio_m(ir,i),oocn_m(ir,i)*1000,ocn_m(ir,i)*1000,0,0,oflux_m(ir,i), flux_m(ir,i) 

     elseif(opt_scheme==3)then
          write(1,'(3i5,9f8.2)')i,myear(i),nmonth(i),obio_m(ir,i),bio_m(ir,i),oocn_m(ir,i),ocn_m(ir,i),ofire_m(ir,i),fire_m(ir,i),  &
                               0,oflux_m(ir,i), flux_m(ir,i)

     elseif(opt_scheme==4)then
          write(1,'(3i5,10f8.2)')i,myear(i),nmonth(i),obio_m(ir,i),bio_m(ir,i),oocn_m(ir,i),ocn_m(ir,i),ofossil_m(ir,i),fossil_m(ir,i),  &
                          ofire_m(ir,i),fire_m(ir,i),oflux_m(ir,i), flux_m(ir,i)
      endif
     endif
   enddo
   close(1)
enddo
endif
end program

subroutine handle_err1(status)
                integer status
                if (status .NE. NF_NOERR) then
                        write(*,*), NF_STRERROR(status), status
                        stop 'program stoped, and you must find the error in your source code.'
                endif
end

subroutine the_end(sy,sm,sd, ey, em, ed, len)
  implicit none
  integer sy, sm, sd, ey, em, ed, len
  integer y, m, d, i
  integer mon(12), mon1(12), mon2(12)
  data mon1/31,28,31,30,31,30,31,31,30,31,30,31/
  data mon2/31,29,31,30,31,30,31,31,30,31,30,31/ 
  
  y=sy
  m=sm
  d=sd
 ! print*, len 
   if(y==2012.or.y==2016.or.y==2000.or.y==2004.or.y==2008.or.y==2020)then
     mon=mon2
   else
     mon=mon1
   endif 
  
  do i=1, len
  !  print*,i, y, m, d
     d=d+1
     if(d>mon(m))then
        d=1
        m=m+1
        if(m>12)then
           m=1
           y=y+1
           if(y==2012.or.y==2016.or.y==2000.or.y==2004.or.y==2008.or.y==2020)then
              mon=mon2
           else
              mon=mon1
           endif
        endif
     endif
  enddo
  ey=y
  em=m
  ed=d
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

subroutine read_map(b, nx, ny,nr)

 include 'netcdf.inc'

 integer nx, ny, nr
 real a(nx,ny), cn(nx,ny), lat(ny), lon(nx)
 real R2(nx,ny)
 integer b(nx,ny,nr)   !, zm(nx,ny,nz,nt)   !, zi(nx,ny,nz,nt)
 integer i, j, ncid
 character(len=30) fname
 character(len=25) vname

! '1_North_America_Boreal','2_North_America_Temperate','3_South_America_Tropical','4_South_America_Temperate', &
! '5_North_Africa','6_South_Africa',    &
! '7_Asia_Boreal','8_Asia_Temperate','9_Southeast_Asia','10_Austrialia','11_Europe', &
! '12_North_Pacific','13_West_Pacific','14_East_Pacific','15_South_Pacific', &
! '16_Arctic_Ocean','17_North_Atlantic','18_Tropical_Atlantic','19_South_Atlantic','20_Southern_Ocean', &
! '21_North_Indian_Ocean','22_South_Indian_Ocean',  &
! '23_China','24_North_America', '25_South_America','26_Asia','27_Africa','28_Atlantic','29_Pacific',  &
!  '30_Indian_Ocean'
! '31_Tropical_land','32_Tropical_ocean','33_Northern_Hemisphere_land','34_Northern_Hemisphere_ocean', &
! '35_Southern_Hemisphere_land','36_Southern_Hemisphere_Ocean','37_Global_land', &
! '38_Global_Ocean', '39_Globe'/

 b=0

do i=1, ny
   lat(i)=-89.5+(i-1)
enddo
do i=1, nx
   lon(i)=0.5+(i-1)
   if(lon(i)>180)lon(i)=lon(i)-360
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
! print*, 'read region ok'

 status = nf_close(ncid)

fname='../input/region_top5_wOcean.nc'
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

fname='../input/region_RECCAP2.nc'
! print*, fname, nx, ny, nr
 status=nf_open(trim(fname),nf_nowrite,ncid)
   if ( status/=nf_noerr ) then
      write (*,*) nf_strerror(status)
      write (*,*) '!!Open NetCDF file Error!!'
       stop
    end if

 vname='region'
  call get_2d_var(ncid,vname,R2,nx,ny)
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
!print*, 'aa'
 
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

 !Africa_india, 47,and Asia_pacific 
 do i=1, nx
    do j=1, ny
      if(lat(j)>=-35.and.lat(j)<=30.and.lon(i)>=-23.and.lon(i)<=90)then
        b(i,j,47)=1
       endif
       if(lat(j)>=-43.and.lat(j)<=40.and.lon(i)>90.and.lon(i)<=160)then
        b(i,j,48)=1
       endif
   enddo
 enddo

 !

 do i=1, nx
     do j=1, ny
       if(R2(i,j)==1)then
        b(i,j,49)=1
       elseif(R2(i,j)==2)then
        b(i,j,50)=1
       elseif(R2(i,j)==3)then
        b(i,j,51)=1
       elseif(R2(i,j)==4)then
        b(i,j,52)=1
       elseif(R2(i,j)==5)then
        b(i,j,53)=1
       elseif(R2(i,j)==6)then
        b(i,j,54)=1
       elseif(R2(i,j)==7)then
        b(i,j,55)=1
       elseif(R2(i,j)==8)then
        b(i,j,56)=1
       elseif(R2(i,j)==9)then
        b(i,j,57)=1
       elseif(R2(i,j)==10)then
        b(i,j,58)=1
       endif
    enddo
 enddo 


end subroutine

function stdev(a, n)
   implicit none
   integer n, i
   real a(n)
   real stdev
   real aa
!   print*, a 
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
!   print*, 'stdev:', stdev
!   pause
   return
  end function



