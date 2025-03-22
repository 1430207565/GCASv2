 subroutine spatial_interp(lat, lon, lat1, lon1, nx, ny, nx1, ny1, c2, c1, m)
  implicit none
  integer nx, ny, nx1, ny1, m
  real lat(ny), lat1(ny1)
  real lon(nx), lon1(nx1)
  real a(nx,ny),a2(nx,ny), a1(nx1,ny1), a_tmp(nx,ny1)
  real c(nx,ny,m),c2(nx,ny,m), c1(nx1,ny1,m), c_tmp(nx,ny1,m)
  integer ix(nx), jx(ny)
  real  dix(nx,3), djx(ny,3)
  real facx(nx,2), facy(ny,2)
  integer xid(nx,2), yid(ny,2)
  integer, external :: dmin2
  real, external :: dxyp
  integer i, j
  real dx, dy, f
  real tb(m), ta(m)

  c=c2

  do j=1, ny
    a(:,j)=dxyp(lat(j),1.0)
  enddo

  do j=1, ny
     c(:,j,:)=c(:,j,:)*dxyp(lat(j),1.0)
 !    print*, dxyp(lat(j),1.0), lat(j)
  enddo

  tb=0
  do i=1, nx

    do j=1, ny
       tb=tb+c(i,j,:)
    enddo
  enddo

!  print*,'before', tb(1:3)

  dx=abs(lon1(2)-lon1(1))
  dy=abs(lat1(3)-lat1(2))

  do i=1, nx
     ix(i)=dmin2(lon(i),lon1,nx1)
  enddo
  do j=1, ny
     jx(j)=dmin2(lat(j),lat1,ny1)
  enddo
  dix=0  
  do i=1, nx
     dix(i,1)=abs(lon(i)-lon1(ix(i)))
    if(ix(i)>1)then
      dix(i,2)=abs(lon(i)-lon1(ix(i)-1))
    else
      dix(i,2)=abs(lon(i)-lon1(nx1)+360)
    endif
    if(ix(i)<nx1)then
      dix(i,3)=abs(lon(i)-lon1(ix(i)+1))
    else
      dix(i,3)=abs(lon(i)-lon1(1)-360)
    endif
   
  enddo

  djx=0
  do j=1, ny
    djx(j,1)=abs(lat(j)-lat1(jx(j)))
    if(jx(j)>1)then
    djx(j,2)=abs(lat(j)-lat1(jx(j)-1))
    endif
    if(jx(j)<ny1)then
    djx(j,3)=abs(lat(j)-lat1(jx(j)+1))
    endif
  enddo

 do j=1, ny
   if(djx(j,1)+0.5<dy/2) then
      facy(j,1)=1
      facy(j,2)=0
      yid(j,1)=jx(j)
      yid(j,2)=0
   elseif(abs(djx(j,1)+djx(j,2)-dy)<0.001)then
      facy(j,1)=djx(j,2)/dy
      facy(j,2)=djx(j,1)/dy
      yid(j,1)=jx(j)
      yid(j,2)=jx(j)-1
   elseif(abs(djx(j,1)+djx(j,3)-dy)<0.001)then
      facy(j,1)=djx(j,3)/dy
      facy(j,2)=djx(j,1)/dy
      yid(j,1)=jx(j)
      yid(j,2)=jx(j)+1
   endif
 enddo

 do i=1, nx
   if(dix(i,1)+0.5<dx/2) then
      facx(i,1)=1
      facx(i,2)=0
      xid(i,1)=ix(i)
      xid(i,2)=0
   elseif(dix(i,3)+0.5<dx/2)then
      facx(i,1)=1
      facx(i,2)=0
      xid(i,1)=ix(i)+1
      xid(i,2)=0
      if(ix(i)+1>nx1) xid(i,1)=1
   elseif(abs(dix(i,1)+dix(i,2)-dx)<0.001.and.dix(i,2)+0.5>dx/2.and.dix(i,1)+0.5>dx/2)then
      facx(i,1)=dix(i,2)/dx
      facx(i,2)=dix(i,1)/dx
      xid(i,1)=ix(i)
      xid(i,2)=ix(i)-1
      if(ix(i)-1<1) xid(i,2)=nx1
   elseif(abs(dix(i,1)+dix(i,3)-dx)<0.001.and.dix(i,3)+0.5>dx/2.and.dix(i,1)+0.5>dx/2)then
      facx(i,1)=dix(i,3)/dx
      facx(i,2)=dix(i,1)/dx
      xid(i,1)=ix(i)
      xid(i,2)=ix(i)+1
      if(ix(i)+1>nx1) xid(i,2)=1
   endif
 enddo

 c_tmp=0
 a_tmp=0
 do j=1, ny
   c_tmp(:,yid(j,1),:)=c_tmp(:,yid(j,1),:)+facy(j,1)*c(:,j,:)
   a_tmp(:,yid(j,1))=a_tmp(:,yid(j,1))+facy(j,1)*a(:,j)
   if(yid(j,2)/=0)then
     c_tmp(:,yid(j,2),:)=c_tmp(:,yid(j,2),:)+facy(j,2)*c(:,j,:)
     a_tmp(:,yid(j,2))=a_tmp(:,yid(j,2))+facy(j,2)*a(:,j)
   endif
 enddo
 c1=0
 a1=0
 do i=1, nx
   c1(xid(i,1),:,:)=c1(xid(i,1),:,:)+facx(i,1)*c_tmp(i,:,:)
   a1(xid(i,1),:)=a1(xid(i,1),:)+facx(i,1)*a_tmp(i,:)
   if(xid(i,2)/=0)then
     c1(xid(i,2),:,:)=c1(xid(i,2),:,:)+facx(i,2)*c_tmp(i,:,:)
     a1(xid(i,2),:)=a1(xid(i,2),:)+facx(i,2)*a_tmp(i,:)
   endif
 enddo

  ta=0
  do i=1, nx1
    do j=1, ny1
       ta=ta+c1(i,j,:)
    enddo
  enddo
  
  do i=1, m
   f=tb(i)/ta(i) 
   c1(:,:,i)=c1(:,:,i)*f
  enddo

!  print*,'after1,', ta(1:3)
!  ta=0
!  do i=1, nx1
!    do j=1, ny1
!       ta=ta+c1(i,j,:)
!    enddo
!  enddo

! print*,'after2,', ta(1:3)

 do i=1, nx1
   do j=1, ny1
      c1(i,j,:)=c1(i,j,:)/a1(i,j)
   enddo
 enddo

! dy=(dy+dx)/2
!
! do j=1, ny1
!    
!    if(abs(lat1(j)-90)<0.001)then
!     !  print*, dxyp(lat1(j)-dy/2.0,dy/2)
!       c1(:,j,:)=c1(:,j,:)/dxyp(lat1(j)-dy/2.0,dy/2)
!    elseif(abs(lat1(j)+90)<0.001)then
!     !  print*, dxyp(lat1(j)+dy/2.0,dy/2)
!       c1(:,j,:)=c1(:,j,:)/dxyp(lat1(j)+dy/2.0,dy/2)
!    else
!     !  print*, dxyp(lat1(j),dy)
!       c1(:,j,:)=c1(:,j,:)/dxyp(lat1(j),dy)
!    endif 
! enddo

 ! print*, c1(1:10,30,1)
 
! pause

end subroutine

integer function dmin2(a,b,n)
        implicit none
        integer n
        real b(n),a
        integer i
        real m,m1

        m=(b(1)-a)*(b(1)-a)
        dmin2=1
        do i=2, n
           m1=(b(i)-a)*(b(i)-a)
           if(m1<m)then
                 m=m1
                dmin2=i
           endif
        enddo

  end function



