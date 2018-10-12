!Fortran examples #7
!MATMUL - multiplies matricies
!DOT_PRODUCT - does what it says
!TRANSPOSE - transpose a matrix
!MAXVAL - max of an array or specified matrix dimension
!MINVAL - min of an array or specified matrix dimension
!SUM - sum of an array or specified matrix dimension


!!demonstrates use of matmul array function and dynamic
!!allocation of array
!program matrixmul
!real,  allocatable, dimension(:,:) :: ra1,ra2,ra3
!integer :: size,row,col
!!initialize the arrays
!print*, 'Shows array  manipulation using SQUARE arrays.'
!print*, 'Allocate the  space for the array at run time.'
!print*, 'Enter the  size of your array'
!read *, size
!allocate(ra1(size,size),ra2(size,size),ra3(size,size))
!print*, 'enter matrix  elements for ra1 row by row'
!call  fill_array(size,ra1)
!print*, 'enter matrix  elements for ra2 row by row'
!call  fill_array(size,ra2)
!!echo the arrays
!print *,'ra1'
!call  outputra(size,ra1)
!print *,'ra2'
!call  outputra(size,ra2)
!!demonstrate the use of matmul and transpose intrinsic !functions
!ra3=matmul(ra1,ra2)
!print *,'matmul of  ra1 and ra2'
!call  outputra(size,ra3)
!ra3=transpose(ra1)
!print *,'transpose of  ra1'
!call  outputra(size,ra3)
!
!!SUM columns and then rows
!write(*,1) sum(ra3,dim=1)
!1    format(100f10.2)
!write(*,4) sum(ra3,dim=2)
!4    format(100f10.2)
!
!!min and max values
!write(*,2) minval(ra3,dim=2)
!2    format(100f10.2)
!!min and max values
!write(*,3) maxval(ra3,dim=2)
!3    format(100f10.2)
!
!deallocate(ra1,ra2,ra3)
!end program matrixmul
!!---------------------------------------------------------
!subroutine  outputra(size,ra)
!implicit none
!!will output a real square array nicely
!integer                               :: size,row,col
!real,dimension(size,size)             :: ra
!character                             :: reply*1
!do row =1,size
!   write(*,10) (ra(row,col),col=1,size)
!   10    format(100f10.2)
!   !as we don't know how many numbers are to be output,
!   !specify  !more than we need - the rest are ignored
!end do
!print*,'__________________________________________________'
!print*,'Hit a key and  press enter to continue'
!read *,reply
!end subroutine  outputra
!!---------------------------------------------------------
!subroutine  fill_array(size,ra)
!implicit none
!!fills the array by prompting from keyboard
!integer        :: row,col,size
!real           :: num
!real, dimension(size,size) :: ra
!do row=1,size
!   do col=1,size
!     print *, row,col
!     read *,num
!     ra(row,col)=num
!   end do
!end do
!end subroutine  fill_array


!CONFIRMATION of (AB)T = (B)T * (A)T
!program trans
!implicit none
!real, allocatable, dimension(:,:) :: a,b,c,d,f
!integer :: size
!print*, 'Shows array  manipulation using SQUARE arrays.'
!print*, 'Allocate the  space for the array at run time.'
!print*, 'Enter the  size of your array'
!read *, size
!allocate(a(size,size),b(size,size),c(size,size),d(size,size),f(size,size))
!call fill_array(size,a)
!call fill_array(size,b)
!c = matmul(a,b)
!c = transpose(c)
!d = matmul(transpose(b),transpose(a))
!call diff(c,d,f,size)
!call check(f,size)
!deallocate
!end program trans
!!################
!subroutine  fill_array(size,ra)
!implicit none
!!fills the array by prompting from keyboard
!integer        :: row,col,size
!real           :: num
!real, dimension(size,size) :: ra
!do row=1,size
!   do col=1,size
!     print *, row,col
!     read *,num
!     ra(row,col)=num
!   end do
!end do
!end subroutine  fill_array
!
!subroutine diff(mat1,mat2,mat3,size)
!implicit none
!real, dimension(size,size) :: mat1,mat2,mat3
!integer :: size,i,j
!mat3=0.
!do i=1,size
!  do j=1,size
!    mat3(i,j) = mat1(i,j) - mat2(i,j)
!  end do
!end do
!end subroutine diff
!
!subroutine check(mat,size)
!implicit none
!real, dimension(size,size) :: mat
!integer :: size, i, j
!do i=1,size
!  do j=1,size
!    if (mat(i,j) /= 0) then
!      print *,'(AB)T = (B)T * (A)T does not hold in this case'
!      stop
!    end if
!  end do
!end do
!print *,'Well done! Your matricies satisfy (AB)T = (B)T * (A)T'
!end subroutine check


!Calculate the cofactors of a read in 3x3 matrix
program cofacts
implicit none
real, dimension(3,3) :: mat,cof
real :: coeff
integer :: i,j
open(1,file='matrix.txt')
read(1,*) mat
mat = transpose(mat)
do i=1,3
write(*,2) (mat(i,j),j=1,3)
2 format(6f10.2)
end do
do i=1,3
  do j=1,3
    call cofactor(i,j,mat,cof)
  end do
end do
do i=1,3
  write(*,3) (cof(i,j),j=1,3)
end do
3 format(10f10.2)

end program cofacts
!#################
subroutine cofactor(i,j,mat,cof)
implicit none
real :: mat(3,3),minor(2,2),cof(3,3)
real :: coeff
integer :: elrow,elcol,i,j,x,y
! cof – the cofactor of matrix mat for element i,j
! remove row i, column j and calculate the determinant of the
! resultant matrix m_(i,j)
! cof = (-1)**(i+j)*m_(i,j)
x=1
print *,'minor of M:',i,j
print *,'-----------------------'
do elrow=1,3
  y=1
  do elcol=1,3
    if (elrow /= i) then
      if (elcol /= j) then
        !print *,mat(elrow,elcol),elrow,elcol,x,y
        minor(x,y) = mat(elrow,elcol)
        y=y+1
      end if
    end if
  end do
  if (i == 1) then
    x = elrow-1
  end if
  if (i == 2) then
    x = 1
  end if
  x=x+1
end do
do x=1,2
  write(*,10) (minor(x,y),y=1,2)
end do
10 format(2f10.2)
print *,'-----------------------'
cof(i,j) = coeff(i,j,minor)
end subroutine cofactor

function coeff(i,j,minor)
!Calculate determinant and then cofactor
implicit none
real,dimension(2,2) :: minor
real,dimension(2) :: A,B
real :: coeff,det
integer :: i,j,x,y,z,k,size
size=2
do y=1,1
  z=y
  do x=1,size
    !print *,minor(x,z)
    A(x) = minor(x,z)
    if (z == size) then
      z=0
    end if
    z=z+1
  end do
end do
do y=2,2
  z=y
  do x=1,size
    !print *,minor(x,z)
    B(x) = minor(x,z)
    if (z == size) then
      z=0
    end if
    z=z+1
  end do
end do
coeff = ((-1)**(i+j)) * (A(1)*A(2)) - (B(1)*B(2))
print *,coeff
end function coeff
