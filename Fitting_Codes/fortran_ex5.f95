!Fortran example code #5
!Learning to deal with arrays

!Find average of 10 numbers
!program av
!real :: x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,average
!read *, x1,x2,x3,x4,x5,x6,x7,x8,x9,x10
!average = (x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10)/10
!print *, 'the average is ',average
!print *, 'the numbers are:'
!print *, x1
!print *, x2
!print *, x3
!print *, x4
!print *, x5
!print *, x6
!print *, x7
!print *, x8
!print *, x9
!print *, x10
!end program av

!Use array of numbers
!program ar1
!implicit none
!real, parameter :: imax = 5
!real, dimension(imax) :: x
!real :: average, sum
!integer :: i
!print *,'Average of 5 numbers:'
!sum = 0.
!do i=1,imax
!  read *, x(i)
!  sum = sum + x(i)
!end do
!average = sum/imax
!print *, average
!print *,'Inputs were: ',x
!end program ar1

!Read in values from a file and use array entries
!program ar2
!implicit none
!integer, allocatable,dimension(:):: vector
!integer :: elements, i, line
!open(1,file='input.txt')
!line = 0
!do
!  read(1,*,end=2) elements
!  line = line+1
!end do
!2 close (1)
!print *,line
!allocate(vector(line))
!open(1,file='input.txt')
!do i=1,line
!  read(1,*) vector(i)
!end do
!do i=1,line
!  print *,vector(i)
!end do
!deallocate(vector)
!end program ar2

!Two lists with input imax values
!program alloc
!implicit none
!real, allocatable :: a(:), b(:)
!integer :: i,imax
!print *,'Length of arrays:'
!read *,imax
!allocate(a(imax),b(imax))
!do i=1,imax
!  a(i) = i**2.
!  b(i) = i**3
!end do
!print *,'a,     b,      elemental sum'
!do i=1,imax
!  print *,a(i),     b(i),     a(i)+b(i)
!end do
!deallocate(a,b)
!end program alloc

!Operators being applied to all elements of an array at once
!program ramagic
!implicit none
!real ,dimension(10)  :: a,b,c
!open(10,file='input.txt')
!read(10,*) a !Array same length as input file therefore can read straight in
!a = (a-1.)/10.
!b = cos(a)
!c = sin(a)
!print *, 'a= ',a
!print *, 'b= ',b
!print *, 'c= ',c
!end program  ramagic

!Mulit-dimensional arrays
!program twodra
!implicit none
!integer,dimension(4,4)      ::  a
!integer                     ::row,col,count
!count = 0
!!creates an array with 4 cols and 4 rows
!!sets col 1 to 1, col2 to 2 and so on
!do row=1,4
!   do col =1,4
!      if (row == col) then
!        a(row,col)=1
!      else
!        a(row,col)=0
!      end if
!    end do
!end do
!do row=1,4
!   do col =1,4
!      print  *,a(row,col)
!   end do
!end do
!end program twodra

!FORMATTING OUTPUTS
! nIm: right justified; integer; m = # characters; n = # of integers per line
! nFm.d: right justified; float; m = # characters to decimal point;
!        n = # of reals per line; d = # characters for decimals
! nEm.d: same as above but for exponentiated functions
! nAm: character specification; n = # of strings to print; m = # of characters

!program format
!implicit none
!!demonstrates use of the format statement
!integer, parameter ::  ikind=selected_real_kind(p=15)
!real , dimension(4)               :: x
!integer, dimension(4)             :: nums
!integer                             :: i
!real(kind=ikind),dimension(4)      :: computed
!!fill up the arrays with something
!do i = 1,4
!   nums(i)           = i * 10
!   computed(i)       =  cos(0.1*i)
!   x(i)              = computed(i)
!end do
!print *,'nums -  integer'
!write(*,1) nums
!1     format(2i10)
!print *, 'x - real'
!      write(*,2) x
!2     format(f6.2)
!print *, 'computed -  double precision'
!      write(*,3)  computed
!3     format(f20.7)
!end program format

!OUTPUT FORMATTED INTO APPROPRIATE TABLE WHEN SAVED OUT
!program exp_format
!implicit none
!real, dimension(11) :: x,ex
!integer :: i
!do i=1,11
!  x(i) = (i-1)/10.
!end do
!ex = exp(x)
!open(1,file='out.txt')
!write(1,*) 'x   exp   exp'
!do i=1,11
!  write(1,*) x(i),ex(i),ex(i)
!end do
!1 format(2f8.2,e18.5)
!end program exp_format

!4x4 neatly formatted identity matrix
program identity
implicit none
integer :: ra(4,4)
integer :: row,column
!implied do loop:
do row=1,4
  do column=1,4
    if (row == column) then
      ra(row,column) = 1
    else
      ra(row,column) = 0
    end if
  end do
end do
do row=1,4
  write(*,1) (ra(row,column),column=1,4)
end do
1 format(4i2)
end program identity
