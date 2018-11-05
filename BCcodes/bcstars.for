C     This program computes bolometric corrections (hence color indices, read
C     further below) in up to five filters (as chosen in selectbc.data), for a
C     sample of stars at fixed reddening. The sample is read from the file
C     input.sample, where each row identifies a star of known logg, [Fe/H],
C     Teff. Reddening E(B-V) is fixed to whatever value (<=0.72) provided when
C     running get4bctables.
C
C     For each star, BCs are printed on the monitor (unit 6) from left to
C     right, in the same order specified in selectbc.data, one star per row in
C     the same order as in input.sample. 
C     The user can then compute any color index (or combination of those) 
C     by simply doing the difference between the BCs in different bands.
C     NOTICE: since bolometric corrections -BCs- are returned, then a colour 
C     index i-j = BC_j - BC_i
C
C     If getbc_std is called, then standard variation of alpha elements is
C     assumed (1 = ialf in selectbc.data). If models with different alpha
C     content want to be used (at the ialf values available in selectbc.data),
C     then instead of calling getbc_std, one of the following subroutines must
C     be called getbc_p04 (2 = ialf), getbc_p00 (3 = ialf), getbc_p04
C     (4 = ialf). Subroutines are located in bcutil.for
C
C     If arbitrary values of alpha enhacement want to be generated, then 
C     get_std, getbc_m04, getbc_p00, and getbc_p04 must be called as appropriate
C     and interpolation be performed. See e.g. iso2cmd.for for an example of
C     interpolation in alpha. We strongly discourage any extrapolation!
C     
C     bcstars.for can be compiled and linked as follows (e.g., ifort compiler):
C
C     ifort bcstars.for bcutil.for -o bcstars.exe
C
C     where bcutil.for contains the soubroutines get_std, getbc_m04, getbc_p00,
C     and getbc_p04.
C

      real*4 a(3),x(3),bcstd(5)!,bcm04(5),bcp00(5),bcp04(5)
      character*6 fil(5)
      character*8 sys(5)
      integer status
      integer iffe

      iffe=0
      
c 1001 format(F4.2,X,F5.2,X,F6.0)
 1002 format(3x,'BC_1',4x,'BC_2',4x,'BC_3',4x,'BC_4',4x,'BC_5')
 1003 format(5(F7.3,X))

      write(6,1002) 
      
c *** read input file in loop till the end
      open(3,file='input.sample')
      readloop: do
c      read(3,1001,IOSTAT=status) gv,fe,temp
      read(3,*,IOSTAT=status) gv,fe,temp
      if (status /= 0) exit

      teff=alog10(temp)

c *** interpolate the bc_std tables to [Fe/H] 
      call getbc_std(iffe,fe,gv,teff,ebv,bcstd,fil,sys,nbc)
      iffe=-1
c     *** now perform the gv and Teff interpolation
      call getbc_std(1,fe,gv,teff,ebv,bcstd,fil,sys,nbc)

      bc1 = bcstd(1)
      bc2 = bcstd(2)
      bc3 = bcstd(3)
      bc4 = bcstd(4)
      bc5 = bcstd(5)

c *** output is written on the monitor. User can send the output to a file 
c *** by simply modifying the code below
      write(6,1003) bc1,bc2,bc3,bc4,bc5

      end do readloop
      stop
      end

