!Piecewise-polynomial Reconstruction from primitive variables implementation
!(using Newton form for the interpolation polynomials)

program prim_reconstruction
implicit none

integer,parameter::nx=25 !number of cells
integer,parameter::np=1000 !number of interpolation function samples
integer,parameter::k=4 !stencil size (set k>=1)
integer,parameter::r=1 !fixed stencil left-shift(no of cells to left of central cell)
integer::is
real::xmin,xmax,dx,x

real::cr(nx,k,2) !co-efficient array, cr(:,:,1)=value of co-efficient,cr(:,:,2)=value of r

integer::i,j
real::Vdiff(1-r:nx-r,1-r:nx-r+k)

xmin=-1.
xmax=1.
dx=(xmax-xmin)/real(nx)

open(unit=10,file='output.txt')

!Compute divided differences
do j=1,k
  do i=1-r,nx-r
   Vdiff(i,i+j)=Vdiv(i,i+j)  !Vdiff(i,i+j):=V[x_i-1/2,...,x_i+j-1/2]
  end do
end do

!go to 98
print*,'i  vbar	 V[x_i-1/2,x_i+1/2]   V[x_i-1/2,..,x_i-1/2+1] V[x_i-1/2,..,x_i-1/2+2]'
do i=1-r,nx-r
  print*,i,vbar(i),Vdiff(i,i+1),Vdiff(i,i+2),Vdiff(i,i+3)
end do
!98 continue


!Compute Newton form polynomial coefficients
do i=1,nx
  do j=1,k
    cr(i,j,1)=Vdiff(i-r,i-r+j)! interpolation coefficients
    cr(i,j,2)=r
  end do 
end do

!print*,'i   cr(1)   cr(2)   cr(3)'
!do i=1,nx
! print*,i,cr(i,1),cr(i,2),cr(i,3)
!end do

!Compute Values of interpolated function and store in file
do i=1,np
  x=xmin+(i-1)*(xmax-xmin)/np
  write(10,*) x,p(x),v(x) 
end do


print*,'Done.'

close(unit=10)

contains

recursive function Vdiv(i,l) result(vx)
integer,intent(in)::i,l
real::vx

if(l>i+1)then
  vx=(Vdiv(i+1,l)-Vdiv(i,l-1))/((l-i)*dx)   !recursive function call
else if(l==i+1)then
  vx=vbar(i)
end if

end

!Original function
function v(x) result(vx)
real,intent(in)::x
real::vx

!go to 99
!step function
if(x<0.)then
  vx=1.
else if(x>=0.)then 
 vx=0.
end if
!99 continue

go to 100
!step function
if(x<-1./3.)then
  vx=0.
else if(x>=-1./3. .and. x<=1./3.)then
  vx=1.
else if(x>1./3.)then
  vx=0.
end if
100 continue

go to 101
!sinusoid
vx=sin(10.*x)
101 continue

end function v 

!Reconstructed function
function p(x) result(px)

real,intent(in)::x
real::px,temp1,temp2
integer::i,j,m,l

i=floor((x-xmin)/dx)+1

px=0.
do j=1,k
  temp2=0.
  do m=0,j-1
   temp1=1.
   do l=0,j-1
     if(l .ne. m)then 
       !temp1=temp1*(x-(xmin+(i-cr(i,j,2)+l)*dx))
       temp1=temp1*(x-(xmin+(i-cr(i,j,2)+l-1)*dx))
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
integer::l

!vbarx=0.5*(v(xmin+(m-1)*dx)+v(xmin+m*dx))
vbarx=0.
do l=1,20
  vbarx=vbarx+v(xmin+(m-1)*dx+dx*l/20)*(1./20)
end do


end function vbar

end program prim_reconstruction
