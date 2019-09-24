!Reconstruction from primitive variables implementation
!(following algorithm outlines in Chi Wang Shu's notes)

program prim_reconstruction
implicit none


integer,parameter::nx=10 !number of cells
integer,parameter::k=3 !stencil size
!integer,parameter::r=2 !left-shift, i.e. no of cells to left of central cell
integer::r,is
real::xmin,xmax,dx

real::cr(nx,k) !co-efficient array

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

print*,'i  vbar	 V[x_i-1/2,x_i+1/2]   V[x_i-1/2,..,x_i-1/2+1] V[x_i-1/2,..,x_i-1/2+2]'
do i=1,nx
  print*,i,vbar(i),Vdiff(i,i+1),Vdiff(i,i+2),Vdiff(i,i+3)
end do

!***********************ENO Stencil Selection***************************
!This step involves choosing the "smoothest" 
!stencil from a selection of stencils. Degree of 
!smoothness is measured by the magnitude of the divided differences.
!For a k-point stencils, there will be (k-1) possible stencils to 
!choose from. Among these,the one that has the smallest kth order 
!divided difference will be chosen by the ENO condition. In most
!cases, it will turn out that the chosen stencil is also the one that
!does not contain a discontinuous point (i.e. the ENO scheme is 
!designed to pick the stencil that deliberately avoids discontinuities
!to acvheive a higher degree of accuracy). 
 
do i=1,nx
  !Start with two point stencil
  j=1
  is=i
  cr(i,1)=Vdiff(is,is+j) !lowest-order stencil
  do j=2,k
    !apply ENO condition to the two (j+1)-point 
    !stencils formed by adding a point on either 
    !side of the preivious lower oder stencil
    if(abs(Vdiff(is-1,is+j-1))<abs(Vdiff(is,is+j)))then
      is=is-1
    end if
    cr(i,j)=Vdiff(is,is+j) !ENO interpolation coefficients
  end do 
end do

print*,'i   cr(1)   cr(2)   cr(3)'
do i=1,nx
 print*,i,cr(i,1),cr(i,2),cr(i,3)
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

!sinusoid
!vx=sin(20.*x)


end function v 


function vbar(m) result(vbarx)

integer,intent(in)::m
real::vbarx

vbarx=0.5*(v(xmin+(m-1)*dx)+v(xmin+m*dx))

end function vbar

end program prim_reconstruction
