!Piecewise-polynomial Reconstruction from primitive variables implementation
!(using Newton form for the interpolation polynomials)

program prim_reconstruction
implicit none


!integer,parameter::nx=500 !number of cells
integer,parameter::np=1000 !number of interpolation function samples
integer,parameter::k=6 !stencil size (set k>=1)
integer,parameter::r=0 !fixed stencil left-shift(no of cells to left of central cell)
integer::nx,is,counter
real::xmin,xmax,dx,x

real::L1err,L1ord

real,allocatable::cr(:,:,:) !co-efficient array, cr(:,:,1)=value of co-efficient,cr(:,:,2)=value of r

integer::i,j
!real::vbar(1:nx)
real,allocatable::v2(:),Vdiff(:,:)

xmin=0.
xmax=1.
!dx=(xmax-xmin)/real(nx)
dx=0.25

open(unit=10,file='output.txt')

do while(counter<=6)

nx=(xmax-xmin)/dx
print*,'nx,dx=',nx,dx

allocate(cr(nx,k,2),v2(1:nx),Vdiff(1-k:nx+k,1-k:nx+k))

!Compute divided differences
do j=1,k
  do i=1-k,nx
   Vdiff(i,i+j)=Vdiv(i,i+j)  !Vdiff(i,i+j):=V[x_i-1/2,...,x_i+j-1/2]
  end do
end do

go to 98
print*,'i  vbar	 V[x_i-1/2,x_i+1/2]   V[x_i-1/2,..,x_i-1/2+1] V[x_i-1/2,..,x_i-1/2+2]'
do i=1-k,nx
  print*,i,vbar(i),Vdiff(i,i+1),Vdiff(i,i+2),Vdiff(i,i+3)
end do
98 continue


!Compute Newton form polynomial coefficients
do i=1,nx
  do j=1,k
    cr(i,j,1)=Vdiff(i-r,i-r+j)! interpolation coefficients
    !cr(i,j,1)=Vdiff(i-r-1,i-r+j-1)
    cr(i,j,2)=r
  end do 
end do

!print*,'i   cr(1)   cr(2)   cr(3)'
!do i=1,nx
! print*,i,cr(i,1),cr(i,2),cr(i,3)
!end do

!Compute Values of interpolated function and store in file
L1err=0.
do i=1,np
  x=xmin+i*(xmax-xmin)/np 
  L1err=L1err+abs(p(x)-v(x))
end do
L1err=L1err/nx

write(10,*) log(dx),log(L1err),5.*log(dx)+20.

!halve the cell interval
dx=dx/2.
counter=counter+1

deallocate(cr,v2,Vdiff)

print*,'counter=',counter

end do 

print*,'Done.'

close(unit=10)

contains


recursive function Vdiv(i,l) result(vx)
integer,intent(in)::i,l
real::vx

if(l>i+1)then
  vx=(Vdiv(i+1,l)-Vdiv(i,l-1))   !recursive function call
else if(l==i+1)then
  vx=vbar(i)
end if

end

!Original function
function v(x) result(vx)
real,intent(in)::x
real::vx


!step function
if(x<0.5)then
  vx=1.
else if(x>=0.5)then 
 vx=0.
end if


go to 99
!step function
if(x<0.35)then
  vx=0.
else if(x>=0.3 .and. x<0.65)then
  vx=1.
else if(x>=0.65)then
  vx=0.
end if
99 continue

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

vbarx=0.5*(v(xmin+(m-1)*dx)+v(xmin+m*dx))

end function vbar

end program prim_reconstruction
