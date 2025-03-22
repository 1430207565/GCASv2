program ensrf

  ! ************************************************************
  ! *                                                          *
  ! *      Program for carbon flux inversion using EnSRF DA       *
  ! *                                                          *
  ! ************************************************************
  ! History:
  ! 17/12/2017  Program was firstly created by Shuzhuang Feng
  ! 10/05/2019  Multiple nested inversions updated by Shuzhuang Feng
  ! 26/09/2019  revised to do carbon inversion by Fei Jiang

  use module_global
  implicit none
  include 'netcdf.inc'
  integer i,j,ios,NumOfObs
  real, allocatable, dimension(:,:,:,:) :: DATA,XNEW
  real, allocatable, dimension(:,:,:) :: XMEAN,XMEAN_old
  REAL,dimension(maxobs)           ::      yobs,Robs,ylat, ylon, yhgt
  integer, dimension(maxobs)       ::  ytype, kx   !ytype,1, surface co2; 2, xco2
  real, allocatable, dimension(:,:) :: HX
  character(len=100) filename
  
  print*,'Ensrf initialize ...'   
  call ensrf_init()
  ! ************************************************************
  
  write(*,"(a,i3,a)") "-do carbon assimilation "//trim(casename)//" with ",ensemble,"ensembles "
  write(*,"(a,i8,a, i3)") "-Starting date and time: ", icdate, "with run days:", rundays

  filename =trim(opt_dir)//'/'//trim(casename)//'/emissions.opt.lambda.nc'

  if(da_step==1)then
   call createlambda(filename, blat,blon, bnx, bny, bnvar)
  endif

  print *, "Get ensrf background..."
  allocate (DATA(bnx,bny,bnvar,ensemble))
  call get_background(DATA)

  print*, "Get ensrf obs..."
  allocate (HX(maxobs,ensemble))
  call get_obs(yobs,Robs,HX,ylat,ylon,yhgt,kx,ytype,NumOfObs)

  print*, "Do Ensrf..."
  allocate(XMEAN(bnx,bny,bnvar),XMEAN_old(bnx,bny,bnvar),XNEW(bnx,bny,bnvar,ensemble))
  call ensrf_noGuassian(yobs,Robs,HX,ylat,ylon,yhgt,kx,ytype,NumOfObs,DATA,XMEAN,XMEAN_old,XNEW)

  print*, "Output opt flux..."
  call ensrf2emission(XMEAN,XMEAN_old,XNEW)

  DEALLOCATE(DATA, HX, XMEAN, XMEAN_old, XNEW)

  print*, "-Sucessful! "
  CLOSE(logun)

end program
