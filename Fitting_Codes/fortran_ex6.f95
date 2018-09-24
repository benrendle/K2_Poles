!Fortran example code #6
!Learning how to use a subroutines and functions
!Subroutine is the fortran version of a python function
!Can be defined after they are called in the main code
!Use 'call' to access the routine
!Arguments are the only things that will be returned, therefore they do not
!require values upon input

!Functions do not require the 'call' command to be used
!Functions reassign the value to their name, not the inputs used

!SUBROUTINE FOR USER CONFIRMATION INPUT
!program output1
!implicit none
!real,dimension(3) :: a,b,c
!!initialise arrays
!a = 1.5
!b = 2.5
!c = 3.5
!write(*,1) 'a',a
!call prompt()
!write(*,1) 'b',b
!call prompt()
!write(*,1) 'c',c
!call prompt()
!write(*,1) 'a*b*c',a * b * c
!1        format(a,3f8.3)
!end program output1
!!++++++++++++++++++++++++++++++++++++++++++++++
!subroutine prompt()
!!prompts for a keypress
!implicit none
!character answer*1
!print *, 'type y to continue or any other key to finish'
!read *, answer
!if (answer /= 'y') stop
!end subroutine prompt


!MULTIPLE INPUTS TO THE SUBROUTINE
!program vols
!!Calculates difference in volume of 2 spheres
!implicit none
!real :: rad1,rad2,vol1,vol2
!character :: response
!do
!  print *, 'Please enter the two radii'
!  read *, rad1,rad2
!  call volume(rad1,vol1)
!  call volume(rad2,vol2)
!  write(*,10) 'The difference in volumes is, ',abs(vol1-vol2)
!  10       format(a,2f10.3)
!  print *, 'Any more? - hit Y for yes, otherwise hit any key'
!  read *, response
!  if (response /= 'Y' .and. response /= 'y') stop
!  end do
!end program vols
!!________________________________________________
!subroutine volume(rad,vol)
!implicit none
!real :: rad,vol,pi
!!calculates the volume of a sphere
!pi=4.0*atan(1.0)
!vol=4./3.*pi*rad*rad*rad
!!It's a little quicker in processing to  do r*r*r than r**3!
!end subroutine volume


!DIFFERENCE BETWEEN THE AREA OF TWO TRIANGLES
!program triangles
!implicit none
!real :: b1,h1,b2,h2,vol1,vol2
!character :: response
!do
!  print *,'Please enter the base and height of triangle 1: '
!  read *,b1,h1
!  print *,'Please enter the base and height of triangle 2: '
!  read *,b2,h2
!  call volume(b1,h1,vol1)
!  call volume(b2,h2,vol2)
!  write(*,1) 'Volume 1: ',vol1
!  write(*,1) 'Volume 2: ', vol2
!  write(*,1) 'Difference: ',abs(vol1-vol2)
!  1 format(a15,f10.5)
!  print *, 'Any more? - hit Y for yes, otherwise enter any other key'
!  read *, response
!  if (response /= 'Y' .and. response /= 'y') stop
!end do
!end program triangles
!!#####################
!subroutine volume(b,h,vol)
!implicit none
!real :: b,h,vol
!vol = 0.5*b*h
!end subroutine volume


!DEGREES TO RADIANS CONVERSION ROUTINE
!program func
!!demonstrates use of user defined functions
!implicit none
!integer, parameter ::  ikind=selected_real_kind(p=15) !Increasing precision
!real (kind=ikind)::  deg,rads
!print *, 'Enter an  angle in degrees'
!read  *, deg
!write(*,10) 'sin =  ',sin(rads(deg))
!write(*,10) 'tan =  ',tan(rads(deg))
!write(*,10) 'cos =  ',cos(rads(deg))
!10   format(a,f10.8)
!end program func
!!_____________________________________________
!function rads(degrees)
!implicit none
!integer, parameter ::  ikind=selected_real_kind(p=15)
!    returns radians
!real (kind=ikind) ::  pi,degrees,rads
!pi=4.0_ikind*atan(1.0_ikind) !No inbuilt value for pi therefore calc.
!rads=(degrees*pi/180.0_ikind) !Actual conversion
!end function rads


!program av_func
!!Program to calculate the average values of a list of numbers from a file
!implicit none
!real, allocatable,dimension(n) :: list(:)
!real :: average !MUST DECLARE FUNCTION AS A VARIABLE TYPE
!integer :: n !Number of values in list
!open(1,file='input.txt')
!n=0
!do
!  read(1,*,end=2)
!  n = n+1
!end do
!2 close (1)
!allocate(list(n))
!open(1,file='input.txt')
!read(1,*) list
!print *,n
!print *,'The average of the input list is: ',average(n,list)
!end program av_func
!!################
!function average(n,list)
!implicit none
!real :: average,sum
!integer :: n,i
!real, dimension(n) :: list
!sum = 0.
!do i=1,n
!  sum = sum + list(i)
!end do
!average = sum/n
!deallocate(list)
!end function average


!SUBROUTINE FOR A FINITE DIFFERENCE MATRIX
program finite_diff
implicit none
integer :: n,i !dimensions of the matrix, counter
integer, allocatable,dimension(:,:) :: matrix
print *,'Please input matrix dimensions (this matrix will be square):'
read *,n
allocate(matrix(n,n))
call finite(matrix,n)
do i=1,n
  write(*,1) matrix(i,1:n) !Print one row at a time for simple formatting
  1 format(100i3)
end do
deallocate(matrix)
end program finite_diff
!######################
subroutine finite(mat,n)
implicit none
integer :: n,i
integer, dimension(n,n) :: mat
mat=0
do i=1,n
  mat(i,i) = 2
  if (i .ne. n) then !Check for boundary conditions. .ne. == not equal to
    mat(i+1,i) = -1
    mat(i,i+1) = -1
  end if
end do
end subroutine finite
