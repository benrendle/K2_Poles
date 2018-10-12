!Fortran examples #8
!Continuation of #7
!Matrix manipulation

!SOLVING THE DET. OF A 3x3 MATRIX - MAKE GENERIC
program determinant
implicit none
real, allocatable, dimension(:,:) :: mat,cof
real :: coeff,Ndet,TwoDet,det
integer :: i,j,size
print *,'Please enter the size of the matrix: '
read *,size
allocate(mat(size,size),cof(size,size))
open(1,file='matrix.txt')
read(1,*) mat
mat = transpose(mat)
do i=1,size
write(*,2) (mat(i,j),j=1,size)
2 format(6f10.2)
end do

if (size == 2) then
  det = TwoDet(mat)
  print *,det
else
  det = Ndet(mat,size)
  print *,'Determinant = ',det
end if
deallocate
end program determinant

!#################
function TwoDet(mat)
implicit none
real,dimension(2,2) :: mat
real :: TwoDet
TwoDet = (mat(1,1)*mat(2,2)) - (mat(1,2)*mat(2,1))
end function TwoDet

function Ndet(mat,size)
implicit none
real, dimension(size,size) :: mat,A,B
real, dimension(size) :: f,g
real :: Ndet,c,d
integer :: size,x,y,z
call positive(mat,size,A,f,c)
call negative(mat,size,B,g,d)
Ndet = c-d
end function Ndet

subroutine positive(mat,size,A,f,c)
!Matrix of all elements to be included in first half of multiplication
implicit none
real, dimension(size,size) ::  mat,A
real, dimension(size) :: f
real :: c
integer :: size,x,y,z
do y=1,size
  z=y
  do x=1,size
    A(x,y) = mat(x,z)
    if (z == size) then
      z=0
    end if
    z=z+1
  end do
end do
f = 1.
do y=1,size
  do x=1,size
    f(y) = f(y) * A(x,y)
  end do
end do
c = sum(f)
end subroutine positive

subroutine negative(mat,size,B,g,d)
!Matrix of all elements to be included in first half of multiplication
implicit none
real, dimension(size,size) ::  mat,B
real, dimension(size) :: g
real :: d
integer :: size,x,y,z
!Matrix of all elements to be included in second half of multiplication and
!subtracted from the first half
do y=1,size
  z=size-y+1
  do x=1,size
    B(x,y) = mat(x,z)
    if (z == 1) then
      z=size+1
    end if
    z=z-1
  end do
end do
g = 1.
do y=1,size
  do x=1,size
    g(y) = g(y) * B(x,y)
  end do
end do
d = sum(g)
end subroutine negative
