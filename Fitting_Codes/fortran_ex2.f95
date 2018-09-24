!Fortran example code #2

!program swap
!Swap values of a and b. Require an intermediate step
!implicit none
!real :: a,b,c
!a = 1
!b = 2
!print *,a,b
!c = a
!a = b
!b = c
!print *,a,b
!end program swap

!program calculate
!implicit none
!real :: x,y,z,a1,a2,a3
!print *,'Enter values for x, y and z: '
!read *,x,y,z
!a1 = (x+y)/(x+z)
!print *,'(x+y)/(x+z) = ',a1
!a2 = x*y*z
!print *,'xyz = ',a2
!a3 = x**(y**z)
!print *,'x^y^z = ',a3
!end program calculate

!FUNCTIONS
!sin(), cos(), tan(), atan(), abs(), sqrt(), exp(), log10()

!program trig
!implicit none
!real :: a,pi
!print *,'Enter an  angle between 0 and 90'
!read *, a
!pi = 4.0*atan(1.0)
!print *,'the sine of  ',a,' is ',sin(a*pi/180)
!end program trig

!TESTING CHOICES
!if  (condition is true) then
!	execute this line
!	and this
!	and so on until we get to …
!end if

!program  test
!implicit none
!!use of a simple menu
!real  :: x,y,answer
!integer  :: choice
!set up the menu – the user may enter  1, 2 or 3
!print  *,'Choose an option'
!print  *,'1    Multiply'
!print  *,'2    Divide'
!print  *,'3    Add'
!read  *,choice
!x=3.4
!y=2.9
!if  (choice == 1) then
!	answer=x*y
!	end if
!if (choice == 2) then
!	answer=x/y
!end if
!if (choice == 3) then
!	answer=x+y
!end if
!print *,'result = ',answer
!end program test

!OPERATORS
!==, /= (not equal to), <, <=, >, >=
!.and., .or.
!if (x == 0) then ... is not a satisfactory test for real numbers.
!Compare an absolute value to a very small, pre-defined value.

program location
implicit none
real :: x
print *,'Please enter a number, x:'
read *,x
if (x > 0 .and. x < 1) print *,'x is between 0 and 1 (not inclusive)'
if (x > 1 .and. x < 10)  print *,'x is between 1 and 10 (not inclusive)'
if (x <= 0 .or. x >= 10 .or. x==1)  print *,'x is outside of the ranges 0 < x < 1 and 1 < x < 10'
end program location
