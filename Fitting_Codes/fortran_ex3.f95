!Fortran example code #3

!program divide
!implicit none
!integer  :: x
!real  :: y
!x  = 1
!!Required to use a float in multiplicaiton or division
!!otherwise the integer value will be returned
!y = x/3.
!print  *, y
!end program divide

!program check
!!Integer and real arithmetic
!implicit none
!real :: x,y
!integer i
!x=2.0
!i=2
!!float*(int/int)
!!int**int = int
!y=x*((2**i)/3)
!print *,y
!!float*(float/int)
!!float**int = float
!y=x*((2.0**i)/3)
!print *,y
!end program check

!DO LOOP EXAMPLE (for looping in python)
!program doLoop
!implicit none
!!real :: y
!integer :: i,x
!do i=-10,10
!  x=i
!  !Loop unnecessary, but safeguards against zero division errors
!  if (abs(x) < 0.000000000001) then !float point zero definition
!    print *,'Inifinity'
!  else
!    print *,1./x
!  end if
!  y = 1./x
!  print *,y
!end do
!end program doLoop

!NESTED DO (for) LOOPS
!program  xytab
!implicit none
!!constructs a table of z=x/y for values of x from 1 to 2 and
!!y from 1 to 4 in  steps of .5
!real :: x, y, z
!integer :: a, b !loop parameters and increments have to be integers (as with python)
!print *, '           x           y           z'
!do  a = 1,2
!  x = a/1.
!	do b = 1,7
!    y = b/2. + 0.5
!		z = x/y
!		print *, x,y,z
!	end do
!end  do
!end  program xytab

program increment
implicit none
integer :: i
real :: x
x = 1.0
do i=1,10
	x = x*i
	print *, i,x
end do
end program increment
