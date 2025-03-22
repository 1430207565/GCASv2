program  post_mozart
 use module_global
 implicit none
 integer, parameter :: nx=360, ny=180,nt=8
 integer, parameter :: nr=46
 integer nm, na
 real flux(nx,ny,nt),bio(nx,ny,nt), ocn(nx,ny,nt), fossil(nx,ny,nt), fire(nx,ny,nt)
 real oflux(nx,ny,nt),obio(nx,ny,nt), oocn(nx,ny,nt),ofossil(nx,ny,nt),ofire(nx,ny,nt)
 real,allocatable, dimension(:,:,:) :: bio_m, bio_a, ocn_m, ocn_a
 real,allocatable, dimension(:,:,:) :: fossil_m, fossil_a, fire_m, fire_a
 real,allocatable, dimension(:,:,:) :: flux_m, flux_a
 real,allocatable, dimension(:,:) :: regf_m, regf_a 
 real,allocatable, dimension(:,:) :: regbio_m, regbio_a
 real,allocatable, dimension(:,:) :: regocn_m, regocn_a
 real,allocatable, dimension(:,:) :: regfossil_m, regfossil_a
 real,allocatable, dimension(:,:) :: regfire_m, regfire_a
 real,allocatable, dimension(:,:) :: lonocn_a, lonbio_a

 real,allocatable, dimension(:,:) :: olonocn_a, olonbio_a
 real,allocatable, dimension(:,:,:) :: obio_m, obio_a, oocn_m, oocn_a
 real,allocatable, dimension(:,:,:) :: ofossil_m, ofossil_a, ofire_m, ofire_a
 real,allocatable, dimension(:,:,:) :: oflux_m, oflux_a
 real,allocatable, dimension(:,:) :: oregf_m, oregf_a
 real,allocatable, dimension(:,:) :: oregbio_m, oregbio_a
 real,allocatable, dimension(:,:) :: oregocn_m, oregocn_a
 real,allocatable, dimension(:,:) :: oregfossil_m, oregfossil_a
 real,allocatable, dimension(:,:) :: oregfire_m, oregfire_a

 integer, allocatable, dimension(:) :: nmonth, nyear, myear
 real lat(ny), lon(nx)
 real fac
 integer b(nx,ny,nr)
 real area(nr), area1(nx,ny)
 integer sd, ed, ir,i, j, k, t, d,m,r, y, sm, em, im, iy
 character(len=100) inputfile, outfile, regname(nr)
 integer mon(12), mon1(12), mon2(12)
 character(len=60) dir
 integer start_year, start_month, start_day, end_year, end_month, end_day
 real, external :: dxyp
 integer narg, tlen, ncycle
 CHARACTER(len=32) ::arg
 character(len=8) sday
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
              'Global_Ocean','Globe','USA','EU','India','Russia',  &
              'N._land_gt30N','T._land_30Nto30S','S._land_lt30S'/

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
allocate(bio_m(nx,ny,nm))
allocate(bio_a(nx,ny,na))
allocate(ocn_m(nx,ny,nm))
allocate(ocn_a(nx,ny,na))
allocate(fossil_m(nx,ny,nm))
allocate(fossil_a(nx,ny,na))
allocate(fire_m(nx,ny,nm))
allocate(fire_a(nx,ny,na))
allocate(regbio_m(nr,nm))
allocate(regbio_a(nr,na))
allocate(lonbio_a(ny,na))
allocate(regocn_m(nr,nm))
allocate(regocn_a(nr,na))
allocate(lonocn_a(ny,na))
allocate(regfossil_m(nr,nm))
allocate(regfossil_a(nr,na))
allocate(regfire_m(nr,nm))
allocate(regfire_a(nr,na))

allocate(obio_m(nx,ny,nm))
allocate(obio_a(nx,ny,na))
allocate(oocn_m(nx,ny,nm))
allocate(oocn_a(nx,ny,na))
allocate(ofossil_m(nx,ny,nm))
allocate(ofossil_a(nx,ny,na))
allocate(ofire_m(nx,ny,nm))
allocate(ofire_a(nx,ny,na))
allocate(oregbio_m(nr,nm))
allocate(oregbio_a(nr,na))
allocate(oregocn_m(nr,nm))
allocate(oregocn_a(nr,na))
allocate(olonbio_a(ny,na))
allocate(olonocn_a(ny,na))
allocate(oregfossil_m(nr,nm))
allocate(oregfossil_a(nr,na))
allocate(oregfire_m(nr,nm))
allocate(oregfire_a(nr,na))
endif

allocate(flux_m(nx,ny,nm))
allocate(flux_a(nx,ny,na))
allocate(regf_m(nr,nm))
allocate(regf_a(nr,na))

allocate(oflux_m(nx,ny,nm))
allocate(oflux_a(nx,ny,na))
allocate(oregf_m(nr,nm))
allocate(oregf_a(nr,na))

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
im=1
k=1
iy=1
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
    !  print*, sday
      call readprior(priordir,obio,oocn,ofossil,ofire,nx,ny,nt,sday,netflux)
      oflux=obio+oocn+ofossil+ofire
      oflux=oflux/6.02e23*12*1.0e4*3600
      obio=obio/6.02e23*12*1.0e4*3600 
      oocn=oocn/6.02e23*12*1.0e4*3600
      ofossil=ofossil/6.02e23*12*1.0e4*3600
      ofire=ofire/6.02e23*12*1.0e4*3600

      write(inputfile,'(a,i4,2i2.2,a3)')trim(dir)//'/emissions.opt.',y,m,d,'.nc'
      ! write(*,'(a,i4,2i2.2,a3)') 'emissions.opt.',y,m,d,'.nc'
      if(nrt==1)then
      call read_nc_nrt(inputfile,flux,nx,ny,nt)

      flux=flux/6.02e23*12*1.0e4*3600
     
      else
      call read_nc(inputfile,bio,ocn,fossil,fire,nx,ny,nt)
      ! molecules/cm2/s ==> gC/m2/hr
      bio=bio/6.02e23*12*1.0e4*3600 
      ocn=ocn/6.02e23*12*1.0e4*3600
      fossil=fossil/6.02e23*12*1.0e4*3600
      fire=fire/6.02e23*12*1.0e4*3600
      flux=bio+ocn+fossil+fire
      endif
     ! std=0
      do t=1, nt
        if(nrt==0)then
        bio_m(:,:,im)=bio_m(:,:,im)+bio(:,:,t)*3
        ocn_m(:,:,im)=ocn_m(:,:,im)+ocn(:,:,t)*3
        fossil_m(:,:,im)=fossil_m(:,:,im)+fossil(:,:,t)*3
        fire_m(:,:,im)=fire_m(:,:,im)+fire(:,:,t)*3
        bio_a(:,:,iy)=bio_a(:,:,iy)+bio(:,:,t)*3
        ocn_a(:,:,iy)=ocn_a(:,:,iy)+ocn(:,:,t)*3
        fossil_a(:,:,iy)=fossil_a(:,:,iy)+fossil(:,:,t)*3
        fire_a(:,:,iy)=fire_a(:,:,iy)+fire(:,:,t)*3

        obio_m(:,:,im)=obio_m(:,:,im)+obio(:,:,t)*3
        oocn_m(:,:,im)=oocn_m(:,:,im)+oocn(:,:,t)*3
        ofossil_m(:,:,im)=ofossil_m(:,:,im)+ofossil(:,:,t)*3
        ofire_m(:,:,im)=ofire_m(:,:,im)+ofire(:,:,t)*3
        obio_a(:,:,iy)=obio_a(:,:,iy)+obio(:,:,t)*3
        oocn_a(:,:,iy)=oocn_a(:,:,iy)+oocn(:,:,t)*3
        ofossil_a(:,:,iy)=ofossil_a(:,:,iy)+ofossil(:,:,t)*3
        ofire_a(:,:,iy)=ofire_a(:,:,iy)+ofire(:,:,t)*3 
        endif
        flux_m(:,:,im)=flux_m(:,:,im)+flux(:,:,t)*3
        flux_a(:,:,iy)=flux_a(:,:,iy)+flux(:,:,t)*3
        
        oflux_m(:,:,im)=oflux_m(:,:,im)+oflux(:,:,t)*3
        oflux_a(:,:,iy)=oflux_a(:,:,iy)+oflux(:,:,t)*3
      enddo
!      global_std(:,:,im)=global_std(:,:,im)+std
    enddo
    im=im+1
  enddo
  iy=iy+1
enddo

regf_m=0
regf_a=0

oregf_m=0
oregf_a=0

if(nrt==0)then
regbio_m=0
regbio_a=0
regocn_m=0
regocn_a=0
regfossil_m=0
regfossil_a=0
regfire_m=0
regfire_a=0
lonbio_a=0
lonocn_a=0

oregbio_m=0
oregbio_a=0
oregocn_m=0
oregocn_a=0
olonbio_a=0
olonocn_a=0

oregfossil_m=0
oregfossil_a=0
oregfire_m=0
oregfire_a=0
endif

do i=1, nx
   do j=1, ny
      fac=dxyp(lat(j),1.0)*1.0e6
      lonbio_a(j,:)=lonbio_a(j,:)+bio_a(i,j,:)*fac
      lonocn_a(j,:)=lonocn_a(j,:)+ocn_a(i,j,:)*fac
      olonbio_a(j,:)=olonbio_a(j,:)+obio_a(i,j,:)*fac
      olonocn_a(j,:)=olonocn_a(j,:)+oocn_a(i,j,:)*fac
   enddo
enddo

do k=1, nr
   do i=1, nx
      do j=1, ny
        fac=b(i,j,k)*dxyp(lat(j),1.0)*1.0e6
        regf_m(k,:)=regf_m(k,:)+flux_m(i,j,:)*fac
        regf_a(k,:)=regf_a(k,:)+flux_a(i,j,:)*fac

        oregf_m(k,:)=oregf_m(k,:)+oflux_m(i,j,:)*fac
        oregf_a(k,:)=oregf_a(k,:)+oflux_a(i,j,:)*fac

        if(nrt==0)then
        regbio_m(k,:)=regbio_m(k,:)+bio_m(i,j,:)*fac
        regbio_a(k,:)=regbio_a(k,:)+bio_a(i,j,:)*fac
     
        regocn_m(k,:)=regocn_m(k,:)+ocn_m(i,j,:)*fac
        regocn_a(k,:)=regocn_a(k,:)+ocn_a(i,j,:)*fac
        
        regfossil_m(k,:)=regfossil_m(k,:)+fossil_m(i,j,:)*fac
        regfossil_a(k,:)=regfossil_a(k,:)+fossil_a(i,j,:)*fac

        regfire_m(k,:)=regfire_m(k,:)+fire_m(i,j,:)*fac
        regfire_a(k,:)=regfire_a(k,:)+fire_a(i,j,:)*fac

        oregbio_m(k,:)=oregbio_m(k,:)+obio_m(i,j,:)*fac
        oregbio_a(k,:)=oregbio_a(k,:)+obio_a(i,j,:)*fac

        oregocn_m(k,:)=oregocn_m(k,:)+oocn_m(i,j,:)*fac
        oregocn_a(k,:)=oregocn_a(k,:)+oocn_a(i,j,:)*fac

        oregfossil_m(k,:)=oregfossil_m(k,:)+ofossil_m(i,j,:)*fac
        oregfossil_a(k,:)=oregfossil_a(k,:)+ofossil_a(i,j,:)*fac

        oregfire_m(k,:)=oregfire_m(k,:)+ofire_m(i,j,:)*fac
        oregfire_a(k,:)=oregfire_a(k,:)+ofire_a(i,j,:)*fac

        endif
      enddo
   enddo
enddo

regf_m=regf_m/1.0e12    !TgC/region/month
regf_a=regf_a/1.0e15    !PgC/region/year

oregf_m=oregf_m/1.0e12    !TgC/region/month
oregf_a=oregf_a/1.0e15    !PgC/region/year

if(nrt==0)then
regbio_m=regbio_m/1.0e12    !TgC/region/month
regbio_a=regbio_a/1.0e15    !PgC/region/year
regocn_m=regocn_m/1.0e12    !TgC/region/month
regocn_a=regocn_a/1.0e15    !PgC/region/year
regfossil_m=regfossil_m/1.0e12    !TgC/region/month
regfossil_a=regfossil_a/1.0e15    !PgC/region/year
regfire_m=regfire_m/1.0e12    !TgC/region/month
regfire_a=regfire_a/1.0e15    !PgC/region/year

oregbio_m=oregbio_m/1.0e12    !TgC/region/month
oregbio_a=oregbio_a/1.0e15    !PgC/region/year
oregocn_m=oregocn_m/1.0e12    !TgC/region/month
oregocn_a=oregocn_a/1.0e15    !PgC/region/year
oregfossil_m=oregfossil_m/1.0e12    !TgC/region/month
oregfossil_a=oregfossil_a/1.0e15    !PgC/region/year
oregfire_m=oregfire_m/1.0e12    !TgC/region/month
oregfire_a=oregfire_a/1.0e15    !PgC/region/year
olonbio_a=olonbio_a/1.0e15
olonocn_a=olonocn_a/1.0e15
lonbio_a=lonbio_a/1.0e15
lonocn_a=lonocn_a/1.0e15
endif

if(nrt==0)then
call output_flux_grid(dir,casename,flux_m,flux_a,bio_m,bio_a,ocn_m,ocn_a,  &
                      fossil_m,fossil_a,fire_m,fire_a,  &
                      oflux_m,oflux_a,obio_m,obio_a,oocn_m,oocn_a,  &
                      ofossil_m,ofossil_a,ofire_m,ofire_a,  &
                     nx,ny,nm,na,lat, lon, opt_scheme)
else
call output_flux_grid_nrt(dir,casename,flux_m,flux_a,oflux_m,oflux_a, &
                     nx,ny,nm,na,lat, lon)
endif

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
  open(1,file=trim(dir)//'/'//trim(regname(ir))//'_annualflux.txt')
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
       write(1,'(i5,2f8.2)')nyear(i),oregf_a(ir,i),regf_a(ir,i)
       if(trim(regname(ir))=='Globe'.or.    &
            trim(regname(ir))=='Global_land'.or.  &
            trim(regname(ir))=='Global_Ocean'.or. &
            trim(regname(ir))=='Europe'.or.    &
            trim(regname(ir))=='North_America_Temperate' .or.   &
            trim(regname(ir))=='China')then
           write(*,'(i5,2f8.2)')nyear(i),oregf_a(ir,i),regf_a(ir,i)
       endif
   enddo
   close(1)
   
enddo
do ir=1, nr
  open(1,file=trim(dir)//'/'//trim(regname(ir))//'_monthlyflux.txt')
  write(1,'(a)') '   ID Year  Mon Net_org Net_opt (PgC/a)'
  do i=1, nm
       write(1,'(3i5,2f8.2)')i,nyear(i),nmonth(i),oregf_m(ir,i),regf_m(ir,i)
   enddo
   close(1)
enddo


elseif(nrt==0)then

   open(1,file=trim(dir)//'/Longitude_Mean_bio_annualflux.txt')
   do j=1, ny
     write(1,'(15f12.4)') lat(j), lonbio_a(j,1:na)
   enddo
   close(1)
  
   open(1,file=trim(dir)//'/Longitude_Mean_ocn_annualflux.txt')
   do j=1, ny
     write(1,'(15f12.4)') lat(j), lonocn_a(j,1:na)
   enddo
   close(1)

   open(1,file=trim(dir)//'/Longitude_Mean_priorbio_annualflux.txt')
   do j=1, ny
     write(1,'(15f12.4)') lat(j), olonbio_a(j,1:na)
   enddo
   close(1)
  
   open(1,file=trim(dir)//'/Longitude_Mean_priorocn_annualflux.txt')
   do j=1, ny
     write(1,'(15f12.4)') lat(j), olonocn_a(j,1:na)
   enddo
   close(1)

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
      open(1,file=trim(dir)//'/'//trim(regname(ir))//'_annualflux.txt')
      if(opt_scheme==1)then
   write(1,'(a)') ' Year  Bio_org Bio_opt     Ocn  Fossil    Fire Net_org Net_opt (PgC/a)'
   elseif(opt_scheme==2)then
   write(1,'(a)') ' Year  Bio_org Bio_opt Ocn_org Ocn_opt  Fossil    Fire Net_org Net_opt (PgC/a)'
   elseif(opt_scheme==3)then
   write(1,'(a)') ' Year  Bio_org Bio_opt Ocn_org Ocn_opt Fir_org Fir_opt  Fossil Net_org Net_opt (PgC/a)'
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
        write(1,'(i5,7f8.2)')nyear(i),oregbio_a(ir,i), regbio_a(ir,i), &
               regocn_a(ir,i),regfossil_a(ir,i),regfire_a(ir,i),oregf_a(ir,i), regf_a(ir,i)
     elseif(opt_scheme==2)then
        write(1,'(i5,8f8.2)')nyear(i),oregbio_a(ir,i), regbio_a(ir,i), &
               oregocn_a(ir,i),regocn_a(ir,i),regfossil_a(ir,i),regfire_a(ir,i),oregf_a(ir,i), regf_a(ir,i)

     elseif(opt_scheme==3)then
          write(1,'(i5,9f8.2)')nyear(i),oregbio_a(ir,i), regbio_a(ir,i), &
               oregocn_a(ir,i),regocn_a(ir,i),oregfire_a(ir,i), regfire_a(ir,i),regfossil_a(ir,i),oregf_a(ir,i), regf_a(ir,i)

     elseif(opt_scheme==4)then
          write(1,'(i5,10f8.2)')nyear(i),oregbio_a(ir,i), regbio_a(ir,i), &
               oregocn_a(ir,i),regocn_a(ir,i),oregfossil_a(ir,i),regfossil_a(ir,i), &
               oregfire_a(ir,i),regfire_a(ir,i),oregf_a(ir,i), regf_a(ir,i)
      endif

         if(trim(regname(ir))=='Globe'.or.    &
            trim(regname(ir))=='Global_land'.or.  &
            trim(regname(ir))=='Global_Ocean'.or. &
            trim(regname(ir))=='Europe'.or.    &
            trim(regname(ir))=='North_America_Temperate' .or.   &
            trim(regname(ir))=='China')then
 
       if(opt_scheme==1)then
        write(*,'(i5,7f8.2)')nyear(i),oregbio_a(ir,i), regbio_a(ir,i), &
               regocn_a(ir,i),regfossil_a(ir,i),regfire_a(ir,i),oregf_a(ir,i), regf_a(ir,i)
       elseif(opt_scheme==2)then
        write(*,'(i5,8f8.2)')nyear(i),oregbio_a(ir,i), regbio_a(ir,i), &
               oregocn_a(ir,i),regocn_a(ir,i),regfossil_a(ir,i),regfire_a(ir,i),oregf_a(ir,i), regf_a(ir,i)
            
       elseif(opt_scheme==3)then
          write(*,'(i5,9f8.2)')nyear(i),oregbio_a(ir,i), regbio_a(ir,i), &
             oregocn_a(ir,i),regocn_a(ir,i),oregfire_a(ir,i),   &
             regfire_a(ir,i),regfossil_a(ir,i),oregf_a(ir,i), regf_a(ir,i)

       elseif(opt_scheme==4)then
          write(*,'(i5,10f8.2)')nyear(i),oregbio_a(ir,i), regbio_a(ir,i), &
             oregocn_a(ir,i),regocn_a(ir,i),oregfossil_a(ir,i),regfossil_a(ir,i), &
             oregfire_a(ir,i),regfire_a(ir,i),oregf_a(ir,i), regf_a(ir,i)
      endif




        endif
     enddo
     close(1)

   enddo

do ir=1, nr
   open(1,file=trim(dir)//'/'//trim(regname(ir))//'_monthlyflux.txt')
   if(opt_scheme==1)then
   write(1,'(a)') '   ID Year  Mon Bio_org Bio_opt     Ocn  Fossil    Fire Net_org Net_opt (TgC/m)'
   elseif(opt_scheme==2)then
   write(1,'(a)') '   ID Year  Mon Bio_org Bio_opt Ocn_org Ocn_opt  Fossil    Fire Net_org Net_opt (TgC/m)'
   elseif(opt_scheme==3)then
   write(1,'(a)') '   ID Year  Mon Bio_org Bio_opt Ocn_org Ocn_opt Fir_org Fir_opt  Fossil Net_org Net_opt (TgC/m)'
   elseif(opt_scheme==4)then
   write(1,'(a)') '   ID Year  Mon Bio_org Bio_opt Ocn_org Ocn_opt Fsl_org Fsl_opt Fir_org Fir_opt Net_org Net_opt (TgC/m)'
   endif

   do i=1, nm
     if(opt_scheme==1)then
        write(1,'(3i5,7f8.2)')i,myear(i),nmonth(i),oregbio_m(ir,i), regbio_m(ir,i), &
               regocn_m(ir,i),regfossil_m(ir,i),regfire_m(ir,i),oregf_m(ir,i), regf_m(ir,i)
     elseif(opt_scheme==2)then
        write(1,'(3i5,8f8.2)')i,myear(i),nmonth(i),oregbio_m(ir,i), regbio_m(ir,i), &
               oregocn_m(ir,i),regocn_m(ir,i),regfossil_m(ir,i),regfire_m(ir,i),oregf_m(ir,i), regf_m(ir,i)

     elseif(opt_scheme==3)then
          write(1,'(3i5,9f8.2)')i,myear(i),nmonth(i),oregbio_m(ir,i), regbio_m(ir,i), &
               oregocn_m(ir,i),regocn_m(ir,i),oregfire_m(ir,i), regfire_m(ir,i),regfossil_m(ir,i),oregf_m(ir,i), regf_m(ir,i)

     elseif(opt_scheme==4)then
          write(1,'(3i5,10f8.2)')i,myear(i),nmonth(i),oregbio_m(ir,i), regbio_m(ir,i), &
               oregocn_m(ir,i),regocn_m(ir,i),oregfossil_m(ir,i),regfossil_m(ir,i),  &
               oregfire_m(ir,i),regfire_m(ir,i),oregf_m(ir,i), regf_m(ir,i)
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

subroutine output_flux_grid(dir,casename,flux_m,flux_a,bio_m,bio_a,ocn_m,ocn_a,  &
           fossil_m,fossil_a,fire_m,fire_a,    &
           oflux_m,oflux_a,obio_m,obio_a,oocn_m,oocn_a,  &
           ofossil_m,ofossil_a,ofire_m,ofire_a,  &
           nx,ny,nm,na,lat, lon, opt_scheme)
 include 'netcdf.inc'
 integer nx, ny, nm, na
 real lat(ny), lon(nx)
 real bio_m(nx,ny,nm), ocn_m(nx,ny,nm)   !gC/m2/month
 real bio_a(nx,ny,na), ocn_a(nx,ny,na)    !gC/m2/year
 real fossil_m(nx,ny,nm), fire_m(nx,ny,nm)   !gC/m2/month
 real fossil_a(nx,ny,na), fire_a(nx,ny,na)    !gC/m2/year
 real flux_a(nx,ny,na), flux_m(nx,ny,nm) 

 real obio_m(nx,ny,nm), oocn_m(nx,ny,nm)   !gC/m2/month
 real obio_a(nx,ny,na), oocn_a(nx,ny,na)    !gC/m2/year
 real ofossil_m(nx,ny,nm), ofire_m(nx,ny,nm)   !gC/m2/month
 real ofossil_a(nx,ny,na), ofire_a(nx,ny,na)    !gC/m2/year
 real oflux_a(nx,ny,na), oflux_m(nx,ny,nm)

 integer biomid, bioaid, ocnmid, ocnaid, netmid,netaid,    &
         fossilmid,fossilaid,firemid, fireaid, fid, regmid, regaid

 integer obiomid, obioaid, oocnmid, oocnaid, onetmid,onetaid,    &
         ofossilmid,ofossilaid,ofiremid, ofireaid

 integer dims3d(3), start(3), count(3)
 integer dims2d(2), start2(2), count2(2)
 integer monthid, yearid, regionid, latid, lonid
 integer monthvarid, yearvarid, regionvarid
 character(len=15) casename, var_fname,var_infile
 character(len=60) dir
 character(len=100) fname
 integer rcode, status, opt_scheme

 fname=trim(dir)//'/'//'posterior.fluxes.'//trim(casename)//'.nc'
 !------------------/----------------------------------------------------------
   ! create the file, clobbering any existing file, add global attributes
   !----------------------------------------------------------------------------
   rcode = nf_create(trim(fname),NF_CLOBBER,fid)
   if (rcode /= NF_NOERR) then
      write(6,*) "ERROR: creating file ",trim(fName)
      write(6,*) "ERROR: ",nf_strerror(rcode)
      stop
   end if
  !... begin to define your dimentions
        status=nf_def_dim(fid, 'nmonth', nm, monthID)  ! the time is unlimited in this program
        if (status .NE. nf_noerr) call handle_err1(status)
  !       print*,'monthid',monthid,nm
        status=nf_def_dim(fid, 'nyear', na, yearID)  ! the time is unlimited in this program
        if (status .NE. nf_noerr) call handle_err1(status)
  !      print*,'yearid',yearid,na

        status=nf_def_dim(fid, 'lat', ny, LATID)
        if (status .NE. nf_noerr) call handle_err1(status)
  !      print*,'latid',latid,ny

        status=nf_def_dim(fid, 'lon', nx, LONID)
        if (status .NE. nf_noerr) call handle_err1(status)
  !      print*,'lonid',lonid,nx

  !     print*,'dimentions defining end'

        dims3d(3)=monthid
        dims3d(2)=latid
        dims3d(1)=lonid

     status=nf_def_var(fid,'netflux_monthly_org', nf_float, 3, dims3d, onetmid)
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_put_att_text(fid, onetmid, 'long_name', 23, 'orginal net month flux')
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_put_att_text(fid, onetmid, 'units', 11, 'gC/m2/month')
     if (status .NE. nf_noerr) call handle_err1(status)     

     status=nf_def_var(fid,'netflux_monthly_opt', nf_float, 3, dims3d, netmid)
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_put_att_text(fid, netmid, 'long_name', 23, 'optmized net month flux')
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_put_att_text(fid, netmid, 'units', 11, 'gC/m2/month')
     if (status .NE. nf_noerr) call handle_err1(status)


     status=nf_def_var(fid,'bio_monthly_org', nf_float, 3, dims3d, obiomid)
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_put_att_text(fid, obiomid, 'long_name', 23, 'orginal bio month flux')
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_put_att_text(fid, obiomid, 'units', 11, 'gC/m2/month')
     if (status .NE. nf_noerr) call handle_err1(status)

     status=nf_def_var(fid,'ocn_monthly_org', nf_float, 3, dims3d, oocnmid)
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_put_att_text(fid, oocnmid, 'long_name', 23, 'orginal ocn month flux')
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_put_att_text(fid, oocnmid, 'units', 11, 'gC/m2/month')
     if (status .NE. nf_noerr) call handle_err1(status)

    status=nf_def_var(fid,'fossil_monthly_org', nf_float, 3, dims3d, ofossilmid)
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_put_att_text(fid, ofossilmid, 'long_name', 23, 'orginal fos month flux')
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_put_att_text(fid, ofossilmid, 'units', 11, 'gC/m2/month')
     if (status .NE. nf_noerr) call handle_err1(status)

     status=nf_def_var(fid,'fire_monthly_org', nf_float, 3, dims3d, ofiremid)
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_put_att_text(fid, ofiremid, 'long_name', 23, 'orginal fir month flux')
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_put_att_text(fid, ofiremid, 'units', 11, 'gC/m2/month')
     if (status .NE. nf_noerr) call handle_err1(status)



     status=nf_def_var(fid,'bio_monthly_opt', nf_float, 3, dims3d, biomid)
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_put_att_text(fid, biomid, 'long_name', 23, 'optmized bio month flux')
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_put_att_text(fid, biomid, 'units', 11, 'gC/m2/month')
     if (status .NE. nf_noerr) call handle_err1(status)
     if(opt_scheme>1)then
     status=nf_def_var(fid,'ocn_monthly_opt', nf_float, 3, dims3d, ocnmid)
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_put_att_text(fid, ocnmid, 'long_name', 23, 'optmized ocn month flux')
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_put_att_text(fid, ocnmid, 'units', 11, 'gC/m2/month')
     if (status .NE. nf_noerr) call handle_err1(status)
     endif
     if(opt_scheme>2)then
    status=nf_def_var(fid,'fire_monthly_opt', nf_float, 3, dims3d, firemid)
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_put_att_text(fid, firemid, 'long_name', 23, 'optmized fir month flux')
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_put_att_text(fid, firemid, 'units', 11, 'gC/m2/month')
     if (status .NE. nf_noerr) call handle_err1(status)
     endif
     if(opt_scheme>3)then
     status=nf_def_var(fid,'fossil_monthly_opt', nf_float, 3, dims3d, fossilmid)
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_put_att_text(fid, fossilmid, 'long_name', 23, 'optmized fos month flux')
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_put_att_text(fid, fossilmid, 'units', 11, 'gC/m2/month')
     if (status .NE. nf_noerr) call handle_err1(status)
     endif
     dims3d(3)=yearid
     dims3d(2)=latid
     dims3d(1)=lonid

     status=nf_def_var(fid, 'netflux_annual_opt', nf_float, 3, dims3d, netaid)
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_put_att_text(fid, netaid, 'long_name', 24, 'optmized net annual flux')
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_put_att_text(fid, netaid, 'units', 10, 'gC/m2/year')
     if (status .NE. nf_noerr) call handle_err1(status)

     status=nf_def_var(fid, 'netflux_annual_org', nf_float, 3, dims3d, onetaid)
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_put_att_text(fid, onetaid, 'long_name', 24, 'orginal net annual flux')
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_put_att_text(fid, onetaid, 'units', 10, 'gC/m2/year')
     if (status .NE. nf_noerr) call handle_err1(status)


!
     status=nf_def_var(fid, 'bio_annual_org', nf_float, 3, dims3d, obioaid)
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_put_att_text(fid, obioaid, 'long_name', 24, 'orginal bio annual flux')
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_put_att_text(fid, obioaid, 'units', 10, 'gC/m2/year')
     if (status .NE. nf_noerr) call handle_err1(status)

    status=nf_def_var(fid, 'ocn_annual_org', nf_float, 3, dims3d, oocnaid)
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_put_att_text(fid, oocnaid, 'long_name', 24, 'orginal ocn annual flux')
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_put_att_text(fid, oocnaid, 'units', 10, 'gC/m2/year')
     if (status .NE. nf_noerr) call handle_err1(status)

     status=nf_def_var(fid, 'fossil_annual_org', nf_float, 3, dims3d, ofossilaid)
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_put_att_text(fid, ofossilaid, 'long_name', 24, 'orginal fos annual flux')
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_put_att_text(fid, ofossilaid, 'units', 10, 'gC/m2/year')
     if (status .NE. nf_noerr) call handle_err1(status)

     status=nf_def_var(fid, 'fire_annual_org', nf_float, 3, dims3d,ofireaid)
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_put_att_text(fid, ofireaid, 'long_name', 24, 'orginal fire annual flux')
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_put_att_text(fid, ofireaid, 'units', 10, 'gC/m2/year')
     if (status .NE. nf_noerr) call handle_err1(status)



     status=nf_def_var(fid, 'bio_annual_opt', nf_float, 3, dims3d, bioaid)
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_put_att_text(fid, bioaid, 'long_name', 24, 'optmized bio annual flux')
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_put_att_text(fid, bioaid, 'units', 10, 'gC/m2/year')
     if (status .NE. nf_noerr) call handle_err1(status)
    if(opt_scheme>1)then
    status=nf_def_var(fid, 'ocn_annual_opt', nf_float, 3, dims3d, ocnaid)
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_put_att_text(fid, ocnaid, 'long_name', 24, 'optmized ocn annual flux')
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_put_att_text(fid, ocnaid, 'units', 10, 'gC/m2/year')
     if (status .NE. nf_noerr) call handle_err1(status)
    endif
    if(opt_scheme>2)then
     status=nf_def_var(fid, 'fire_annual_opt', nf_float, 3, dims3d, fireaid)
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_put_att_text(fid, fireaid, 'long_name', 24, 'optmized fir annual flux')
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_put_att_text(fid, fireaid, 'units', 10, 'gC/m2/year')
     if (status .NE. nf_noerr) call handle_err1(status)
     endif
    if(opt_scheme>3)then
     status=nf_def_var(fid, 'fossil_annual_opt', nf_float, 3, dims3d, fossilaid)
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_put_att_text(fid, fossilaid, 'long_name', 24, 'optmized foss annual flux')
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_put_att_text(fid, fossilaid, 'units', 10, 'gC/m2/year')
     if (status .NE. nf_noerr) call handle_err1(status)
    endif


     status=nf_def_var(fid, 'lat', nf_float, 1, latid, latvarid)
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_def_var(fid, 'lon', nf_float, 1, lonid, lonvarid)
     if (status .NE. nf_noerr) call handle_err1(status)

      status=nf_enddef(fid)

   count(3)=nm
   count(2)=ny
   count(1)=nx
   start=1

   status=nf_put_vara_real(fid, onetmid, start, count, oflux_m)
   status=nf_put_vara_real(fid, obiomid, start, count, obio_m)
   status=nf_put_vara_real(fid, oocnmid, start, count, oocn_m)
   status=nf_put_vara_real(fid, ofossilmid, start, count, ofossil_m)
   status=nf_put_vara_real(fid, ofiremid, start, count, ofire_m)


   status=nf_put_vara_real(fid, netmid, start, count, flux_m)
   status=nf_put_vara_real(fid, biomid, start, count, bio_m)
   if(opt_scheme>1)then
   status=nf_put_vara_real(fid, ocnmid, start, count, ocn_m)
   endif
   if(opt_scheme>2)then
   status=nf_put_vara_real(fid, firemid, start, count, fire_m)
   endif
   if(opt_scheme>3)then
   status=nf_put_vara_real(fid, fossilmid, start, count, fossil_m)
   endif
   count(3)=na
   count(2)=ny
   count(1)=nx

   status=nf_put_vara_real(fid, onetaid, start, count, oflux_a)
   status=nf_put_vara_real(fid, obioaid, start, count, obio_a)
   status=nf_put_vara_real(fid, oocnaid, start, count, oocn_a)
   status=nf_put_vara_real(fid, ofossilaid, start, count, ofossil_a)
   status=nf_put_vara_real(fid, ofireaid, start, count, ofire_a)

   status=nf_put_vara_real(fid, netaid, start, count, flux_a)
   status=nf_put_vara_real(fid, bioaid, start, count, bio_a)
   if(opt_scheme>1)then
   status=nf_put_vara_real(fid, ocnaid, start, count, ocn_a)
   endif
   if(opt_scheme>2)then
   status=nf_put_vara_real(fid, fireaid, start, count, fire_a)
   endif
   if(opt_scheme>3)then
   status=nf_put_vara_real(fid, fossilaid, start, count, fossil_a)
   endif
   status=nf_put_vara_real(fid, latvarid, 1, ny, lat)
   status=nf_put_vara_real(fid, lonvarid, 1, nx, lon)
   status=nf_close(fid)

end subroutine
 

subroutine output_flux_grid_nrt(dir,casename,bio_m,bio_a, obio_m, obio_a,nx,ny,nm,na,lat, lon)
 include 'netcdf.inc'
 integer nx, ny, nm, na
 real lat(ny), lon(nx)
 real bio_m(nx,ny,nm), obio_m(nx,ny,nm)   !gC/m2/month
 real bio_a(nx,ny,na), obio_a(nx,ny,na)    !gC/m2/year
 integer biomid, bioaid, obiomid, obioaid, fid, regmid, regaid
 integer dims3d(3), start(3), count(3)
 integer dims2d(2), start2(2), count2(2)
 integer monthid, yearid, regionid, latid, lonid
 integer monthvarid, yearvarid, regionvarid
 character(len=15) casename, var_fname,var_infile
 character(len=60) dir
 character(len=100) fname
 integer rcode, status
 
 fname=trim(dir)//'/'//'posterior.fluxes.'//trim(casename)//'.nc'
 !------------------/----------------------------------------------------------
   ! create the file, clobbering any existing file, add global attributes
   !----------------------------------------------------------------------------
   rcode = nf_create(trim(fname),NF_CLOBBER,fid)
   if (rcode /= NF_NOERR) then
      write(6,*) "ERROR: creating file ",trim(fName)
      write(6,*) "ERROR: ",nf_strerror(rcode)
      stop
   end if

  !... begin to define your dimentions
        status=nf_def_dim(fid, 'nmonth', nm, monthID)  ! the time is unlimited in this program
        if (status .NE. nf_noerr) call handle_err1(status)
  !      print*,'monthid',monthid,nm
        status=nf_def_dim(fid, 'nyear', na, yearID)  ! the time is unlimited in this program
        if (status .NE. nf_noerr) call handle_err1(status)
  !      print*,'yearid',yearid,na

        status=nf_def_dim(fid, 'nlat', ny, LATID)
        if (status .NE. nf_noerr) call handle_err1(status)
  !      print*,'latid',latid,ny

        status=nf_def_dim(fid, 'nlon', nx, LONID)
        if (status .NE. nf_noerr) call handle_err1(status)
  !      print*,'lonid',lonid,nx

  !     print*,'dimentions defining end' 

        dims3d(3)=monthid
        dims3d(2)=latid
        dims3d(1)=lonid

     status=nf_def_var(fid,'netflux_monthly_opt', nf_float, 3, dims3d, biomid)
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_put_att_text(fid, biomid, 'long_name', 23, 'optmized net month flux')
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_put_att_text(fid, biomid, 'units', 11, 'gC/m2/month')
     if (status .NE. nf_noerr) call handle_err1(status)

     status=nf_def_var(fid,'netflux_monthly_org', nf_float, 3, dims3d, obiomid)
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_put_att_text(fid, obiomid, 'long_name', 23, 'orginal net month flux')
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_put_att_text(fid, obiomid, 'units', 11, 'gC/m2/month')
     if (status .NE. nf_noerr) call handle_err1(status)

     dims3d(3)=yearid
     dims3d(2)=latid
     dims3d(1)=lonid

     status=nf_def_var(fid, 'netflux_annual_opt', nf_float, 3, dims3d, bioaid)
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_put_att_text(fid, bioaid, 'long_name', 24, 'optmized net annual flux')
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_put_att_text(fid, bioaid, 'units', 10, 'gC/m2/year')
     if (status .NE. nf_noerr) call handle_err1(status)

     status=nf_def_var(fid, 'netflux_annual_org', nf_float, 3, dims3d, obioaid)
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_put_att_text(fid, obioaid, 'long_name', 24, 'orginal net annual flux')
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_put_att_text(fid, obioaid, 'units', 10, 'gC/m2/year')
     if (status .NE. nf_noerr) call handle_err1(status)

     status=nf_def_var(fid, 'lat', nf_float, 1, latid, latvarid)
     if (status .NE. nf_noerr) call handle_err1(status)
     status=nf_def_var(fid, 'lon', nf_float, 1, lonid, lonvarid)
     if (status .NE. nf_noerr) call handle_err1(status)

      status=nf_enddef(fid)

   count(3)=nm
   count(2)=ny
   count(1)=nx
   start=1

   status=nf_put_vara_real(fid, biomid, start, count, bio_m)
   status=nf_put_vara_real(fid, obiomid, start, count, obio_m)

   count(3)=na
   count(2)=ny
   count(1)=nx
   status=nf_put_vara_real(fid, bioaid, start, count, bio_a)
   status=nf_put_vara_real(fid, obioaid, start, count, obio_a)

   status=nf_put_vara_real(fid, latvarid, 1, ny, lat)
   status=nf_put_vara_real(fid, lonvarid, 1, nx, lon)
   status=nf_close(fid)

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
 real a(nx,ny), cn(nx,ny)
 integer b(nx,ny,nr)   !, zm(nx,ny,nz,nt)   !, zi(nx,ny,nz,nt)
 integer i, j, ncid
 real lat(ny)
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
! '38_Global_Ocean', '39_Globe','40USA','41EU','42India','43Russia',  &
!'44N._land_gt30N','45T._land_30Nto30S','46S._land_lt30S'/

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

end subroutine




