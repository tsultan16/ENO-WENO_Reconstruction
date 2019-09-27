!Piecewise-polynomial ENO Reconstruction from primitive variables implementation
!(following algorithm outlines in Chi Wang Shu's notes)

program prim_reconstruction_eno
implicit none


integer,parameter::nx=20 !number of cells
integer,parameter::np=1000 !number of interpolation function samples
integer,parameter::k=4 !stencil size (set k>=1)
integer::r,is
real::xmin,xmax,dx,x

!real::cr(nx,k,2) !co-efficient array, cr(:,:,1)=value of co-efficient,cr(:,:,2)=value of r
integer::r(nx)  !eno stencil left-most cell index

integer::i,j
!real::vbar(1:nx)
real::Vdiff(1-k:nx+k,1-k:nx+k)

xmin=-1.
xmax=1.
dx=(xmax-xmin)/real(nx)

open(unit=10,file='output.txt')

!Compute divided differences
do j=1,k
  do i=1-k,nx
   Vdiff(i,i+j)=Vdiv(i,i+j)  !Vdiff(i,i+j):=V[x_i-1/2,...,x_i+j-1/2]
  end do
end do

go to 97
print*,'i  vbar	 V[x_i-1/2,x_i+1/2]   V[x_i-1/2,..,x_i-1/2+1] V[x_i-1/2,..,x_i-1/2+2]'
do i=1,nx
  print*,i,vbar(i),Vdiff(i,i+1),Vdiff(i,i+2),Vdiff(i,i+3)
end do
97 continue

!***********************ENO Stencil Selection***************************
!This step involves choosing the "smoothest" 
!stencil from a selection of stencils. Degree of 
!smoothness is measured by the magnitude of the divided differences.
!For a k-point stencils, there will be (k-1) possible stencils to 
!choose from for a given cell. Among these,the one that has the smallest  
!kth order divided difference will be chosen by the ENO condition. In most
!cases, it will turn out that the chosen stencil is also the one that
!does not contain a discontinuous point (i.e. the ENO scheme is 
!designed to pick the stencil that deliberately avoids discontinuities
!to acheive a higher degree of accuracy). 
!************************************************************************ 
do i=1,nx
  !Start with two point stencil
  j=1
  is=i
  cr(i,1,1)=Vdiff(is,is+j) !lowest-order stencil
  cr(i,1,2)=i-is !r=0 for lowest-order stencil
  do j=2,k
    !apply ENO condition to the two (j+1)-point 
    !stencils formed by adding a point on either 
    !side of the preivious lower oder stencil
    if(abs(Vdiff(is-1,is+j-1))<abs(Vdiff(is,is+j)))then
      is=is-1
    end if
    cr(i,j,1)=Vdiff(is,is+j) !ENO interpolation coefficients
    cr(i,j,2)=i-is
  end do 
end do

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

go to 99
!step function
if(x<0.)then
 vx=1.
else  
 vx=0.
end if
99 continue


!go to 100
!step function
if(x<-1./3.)then
  vx=0.
else if(x>=-1./3. .and. x<=1./3.)then
  vx=1.
else if(x>1./3.)then
  vx=0.
end if
!100 continue

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
integer::l,n

n=50

!vbarx=0.5*(v(xmin+(m-1)*dx)+v(xmin+m*dx))
vbarx=0.
do l=1,n
vbarx=vbarx+v(xmin+(m-1)*dx+dx*(l-1)/n)*(1./n)
end do

end function vbar


end program prim_reconstruction_eno
