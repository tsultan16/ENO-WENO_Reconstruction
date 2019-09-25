!Piecewise-polynomial Reconstruction from primitive variables implementation
!(using Newton form for the interpolation polynomials)

program prim_reconstruction
implicit none


integer,parameter::nx=4 !number of cells
integer,parameter::np=1000 !number of interpolation function samples
integer,parameter::k=2 !stencil size (set k>=1)
integer,parameter::r=0 !fixed stencil left-shift(no of cells to left of central cell)
integer::is
real::xmin,xmax,dx,x

real::cr(nx,k,2) !co-efficient array, cr(:,:,1)=value of co-efficient,cr(:,:,2)=value of r

integer::i,j
!real::vbar(1:nx)
real::v2(1:nx),Vdiff(nx,nx+k)

xmin=0.
xmax=1.
dx=(xmax-xmin)/real(nx)

open(unit=10,file='output.txt')

!Compute divided differences
do j=1,k
  do i=1,nx
   Vdiff(i,i+j)=Vdiv(i,i+j)  !Vdiff(i,i+j):=V[x_i-1/2,...,x_i+j-1/2]
  end do
end do

!print*,'i  vbar	 V[x_i-1/2,x_i+1/2]   V[x_i-1/2,..,x_i-1/2+1] V[x_i-1/2,..,x_i-1/2+2]'
!do i=1,nx
!  print*,i,vbar(i),Vdiff(i,i+1),Vdiff(i,i+2),Vdiff(i,i+3)
!end do

!Compute Newton form polynomial coefficients
do i=1,nx
  do j=1,k
    cr(i,j,1)=Vdiff(i-r,i-r+j) !ENO interpolation coefficients
    cr(i,j,2)=r
  end do 
end do

!print*,'i   cr(1)   cr(2)   cr(3)'
!do i=1,nx
! print*,i,cr(i,1),cr(i,2),cr(i,3)
!end do

!Compute Values of interpolated function and store in file
do i=1,np
  x=xmin+i*(xmax-xmin)/np
  write(10,*) x,p(x),v(x) 
end do


print*,'Done.'

close(unit=10)

contains


recursive function Vdiv(i,k) result(vx)
integer,intent(in)::i,k
real::vx

if(k>i+1)then
  vx=(Vdiv(i+1,k)-Vdiv(i,k-1))!/dx   !recursive function call
else if(k==i+1)then
  vx=vbar(i)
end if

end

!Original function
function v(x) result(vx)
real,intent(in)::x
real::vx

!step function
if(x<0.35)then
  vx=0.
else if(x>=0.3 .and. x<0.65)then
  vx=1.
else if(x>=0.65)then
  vx=0.
end if

go to 100
!sinusoid
vx=sin(20.*x)
100 continue

end function v 

!Reconstructed function
function p(x) result(px)

real,intent(in)::x
real::px,temp1,temp2
integer::i,j,m,l

i=floor(x/dx)+1

px=0.
do j=1,k
  temp2=0.
  do m=0,j-1
   temp1=1.
   do l=0,j-1
     if(l .ne. m)then 
       temp1=temp1*(x-(xmin+(i-cr(i,j,2)+l)*dx))
     end if
   end do
   temp2=temp2+temp1
  end do
  px=px+cr(i,j,1)*temp2
end do


end function p


function vbar(m) result(vbarx)

integer,intent(in)::m
real::vbarx

vbarx=0.5*(v(xmin+(m-1)*dx)+v(xmin+m*dx))

end function vbar

end program prim_reconstruction
