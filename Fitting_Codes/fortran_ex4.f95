!Fortran example code #4

!program readdata
!implicit none
!!reads data from a file called input.txt
!real :: x,y,z
!open(10,file='input.txt')
!        !Assign a device to the file to read (10 in this case)
!        !Must be a positive integer
!        read(10,*) x,y,z
!print *,x,y,z
!end program  readdata

!program evenodd
!implicit none
!real :: x
!integer :: i
!open(10,file='input.txt')
!  do i=1,10
!    read(10,*) x
!    if (mod(x,2.) > 0.) print *, x,' is odd'
!    if (mod(x,2.) == 0) print *, x,' is even'
!  end do
!end program evenodd

!Save out value depending on it being even or odd
!program io
!implicit none
!real :: x
!integer :: nlines,i
!open(1,file='input.txt')
!!loop to determine the length of the file
!nlines = 0
!do
!  read(1,*,end=10)
!  nlines=nlines+1
!end do
!10 close (1)
!print *,nlines
!open(1,file='input.txt')
!open(2,file='odd.txt')
!open(3,file='even.txt')
!!loop to save the values in the relevant locations
!do
!  read(1,*,end=11) x  !Including 'end' allows one to loop to the end of the
!                      !file without having to determine the file lenght
!  if (mod(x,2.) > 0.) write(2,*) x
!  if (mod(x,2.) == 0) write(3,*) x
!end do
!11 close (1)
!end program io

!Allowing for extended precision in calculations
!program extended
!implicit none
!!ikind is a new parameter that have the desired precision. p = 6,15,18
!!parameters can't be changed once declared
!integer, parameter  :: ikind=selected_real_kind(p=15)
!real (kind=ikind)  :: sum,x
!integer :: i
!sum=0.0
!do i=1,100
!  x=i
!  sum = sum +  1.0/(x**6)
!end do
!print *, sum
!end program  extended

!program extendedconstants
!!demonstrates use of extended precision
!implicit none
!integer, parameter ::  ikind=selected_real_kind(p=18)
!real (kind=ikind) :: val,x,y
!val=10/3
!print*,val                  !10/3 calculated as  integer  - wrong!
!x=10.0
!y=3.0
!val=x/y              !x/y assigned to extended  precision - right!
!print*,val
!val=10.0_ikind/3              !extend precision constant -  right!
!print*,val
!val=10.0/3.0                              !real constants - wrong!
!print*,val
!val = .12345678901234567890               !real constants - wrong!
!print *, val
!val = .12345678901234567890_ikind   !ext precision consts - right!
!print *, val
!end program  extendedconstants

!program magnitude
!implicit none
!integer, parameter :: ikind=selected_real_kind(p=15)
!integer i,maxpower,response
!print *,'Effect of magnitude on integer, real or extended precision'
!print *,'Type 1 for integer, 2 for real or 3 for extended precision'
!read  *, response
!print *,'Enter the maximum power'
!read *,maxpower
!do i=1,maxpower
!  if (response==1) then
!     print *,i,2**i
!  else if (response==2) then
!     print *,i,2.0**i
!  else if (response==3) then
!      print *,i,2.0_ikind**i   !Need the _ikind otherwise the 2.0 has 6pt
!                                !precision, limiting the operation
!  else
!      print *,'invalid response'
!      stop
!  end if
!end do
!end program magnitude

!While looping structure - do loop with break criteria
program  whileloop
implicit none
integer, parameter ::  ikind=selected_real_kind(p=15)
real (kind=ikind) ::  sum,previoussum,x,smallnumber,error
integer :: i
sum = 0.0
previoussum = 0.0
smallnumber = 10.0_ikind**(-15.0)
do i=1,1000
   x=i
   sum = sum + 1.0 /(x**6)
   error=abs(sum-previoussum)
   if  (error<smallnumber) then
       print *,'sum  ',sum,' number of loops ',i
       exit
    end if
   previoussum = sum
end do
end program  whileloop
