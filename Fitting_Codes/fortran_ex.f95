!Fortran examples code
program fortran_ex
!Define all parameters sequentially.
!Code works in order there require orderly entering of parameters
!Define all parameters to be used at the start of the code

!NUMERICAL TESTING BASICS
!real :: answer,x,y
!real :: z
!real :: a,b,c
!print *,'My first program'
!print  *, 'Enter two numbers'
!read  *, x
!read  *, y
!answer=x+y
!print  *, 'The total is ', answer
!print *, 'enter the values x,y and z'
!read *, x,y,z   !Enter new values everytime this appears
!print *, 'the values you typed are for z,y,x are: ',z,y,x

!print *,'a = b + c; please enter a value for c:'
!read *,c
!b = 5
!a = b + c
!print *,a

!DIFFERENT CHARACTER TYPES
!This example shows the use of integer  and character variables
!implicit none !Tells Fortran to check all parameters are correctly defined
!integer   ::  pounds,pence,total
!character :: name*10 !Number defines the maximum length of the string
!print *,'What is your name?'
!read *,name
!print *, 'Hi ',name,'! Enter number of pounds and  pence'
!read *, pounds,pence
!total =100 * pounds + pence
!print *,'the total money in pence is ',total

!TEST CODE #1
implicit none
real :: a,b,c
integer :: x,y,z
character :: name*10, surname*10
print *,'First name:'
read *,name
print *,'Surname:'
read *,surname
print *,name,surname,', you have been selected to enter some random numbers!'
print *,'Please enter 2 floats and 3 integers:'
read *,b,c,x,y,z
a = b/(x-y) + c**2/(y-z)
print *,'Well done',name,' you have magically generated the number: ',a
end program fortran_ex
