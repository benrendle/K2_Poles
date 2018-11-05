C *** This program generates "bc_std.data", "bc_p04.data", "bc_p00.data"
C     and "bc_m04.data" for the selected photometric systems and filters
C     identified in "selectbc.data".  These output tables are produced for
C     the E(B-V) value that is specified by the user.
 1000 format(' enter E(B-V) value for which BCs are to be generated')
      write(6,1000)
      read(5,*) ebmv
      ialf=0
      do 10 i=1,4
      ialf=ialf+1
      call getinputdata(ialf,ebmv)
   10 continue
      stop
      end
      subroutine getinputdata(ialf,xred)
C *** This program is used in conjunction with "selectbc.data" to set up the
C     BC tables for the photometric system and filter bandpasses of interest,
C     assuming E(B-V) = 0.0, 0.12, 0.24, 0.36, 0.48, 0.60, and 0.72.  The
C     resultant 7 tables are then used to produce the BC tables for the same
C     filter bandpasses and any selected value of E(B-V) in the range from
C     0.0 to 0.72.
C
C *** Regarding "selectbc.data":
C       1st line: integer is set to a value of 1, 2, 3, or 4 to specify
C                 which variation of [alpha/Fe] with [Fe/H] is to be assumed.
C       2nd line: integer is set to a value of 1, 2, ..., 5 to indicate the
C                 number of filter bandpasses for which BC values will be
C                 defined.  (The maximum value permitted is 5.)
C       next 1, 2, ..., 5 lines: the pairs of integers define the photometric
C                 system and filter for which BC values will be defined: the
C                 menu of possible values are given below the horizontal 
C                 line (e.g., 5 33 selects the sdss photometric system and
C                 the sdss i filter).
C
C *** The output of this program consists of seven data files named
C     inputbc_r00.data, inputbc_r12.data, ..., inputbc_r72.data.  If these
C     files already exist, they will be overwritten.     
C
C *** NOTE that "slash" as defined in the data statement (below) **MUST**
C     be defined to be a forward slash (i.e., "slash/'/'/" if this program
C     is executed on a linux or UNIX computer.  (A backslash is used in the
C     case of a WINDOWS operating system.)
C--------1---------2---------3---------4---------5---------6---------7--
      real*4 bc(31),teff(31),gv(13)
      integer*4 itmm(13)
      character*50 labtab
      character*42 fnme,fzero
      character*26 redden
      character*12 filter(88),filtr1(41),filtr2(47),xfi(5),finme
      character*11 hdr
      character*9 alpha(8),alnme
      character*8 system(27),xps(5),psnme,outtab
      character*5 ftype
      character*3 ebmv(7),rednme
      character*1 slash
      data alpha/'alpha_std','alpha_p04','alpha_p00','alpha_m04',
     1  'al2fe_std','al2fe_p04','al2fe_p00','al2fe_m04'/,ebmv/
     2  'r00','r12','r24','r36','r48','r60','r72'/,system/'2mass   ',
     3  'hst_ab  ','hst_st  ','hst_vega','sdss    ','ubvri12 ',
     4  'ubvri90 ','MIRIab  ','MIRIst  ','MIRIvega','MODAab  ',
     5  'MODAst  ','MODAvega','MODABab ','MODABst ','MODABveg',
     6  'MODBab  ','MODBst  ','MODBvega','skymap  ','tycho   ',
     7  'PanSTARR','GaiaPAB ','GaiaPveg','GaiaRAB ','GaiaRveg',
     8  'GaiaDR2 '/
      Data fzero/'                                          '/,
C *** on linux/UNIX machines, define slash as '/' instead of '\' 
     1  slash/'/'/,outtab/'inputbc_'/,hdr/'header.data'/,ftype/'.data'/
      data filtr1/'2mass_J     ','2mass_H     ','2mass_K     ',
     1  'acs_f435w   ','acs_f475w   ','acs_f555w   ','acs_f606w   ',
     2  'acs_f814w   ','wfc3_f218w  ','wfc3_f225w  ','wfc3_f275w  ',
     3  'wfc3_f336w  ','wfc3_f350lp ','wfc3_f390m  ','wfc3_f390w  ',
     4  'wfc3_f438w  ','wfc3_f475w  ','wfc3_f547m  ','wfc3_f555w  ',
     5  'wfc3_f606w  ','wfc3_f625w  ','wfc3_f775w  ','wfc3_f814w  ',
     6  'wfc3_f850lp ','wfc3ir_f098m','wfc3ir_f110w','wfc3ir_f125w',
     7  'wfc3ir_f140w','wfc3ir_f160w','sdss_u      ','sdss_g      ',
     8  'sdss_r      ','sdss_i      ','sdss_z      ','jc_U        ',
     9  'jc_B        ','jc_V        ','jc_R        ','jc_I        ', 
     *  'jc_UX       ','jc_BX       '/,nfil1/41/,nfil2/47/
      data filtr2/'miri_f560w  ','miri_f770w  ','miri_f1000w ',
     1  'miri_f1130w ','miri_f1280w ','miri_f1500w ','jwst_f070w  ',
     2  'jwst_f090w  ','jwst_f115w  ','jwst_f140m  ','jwst_f150w2 ',
     3  'jwst_f150w  ','jwst_f162m  ','jwst_f182m  ','jwst_f200w  ',
     4  'jwst_f210m  ','jwst_f250m  ','jwst_f277w  ','jwst_f300m  ',
     5  'jwst_f322w2 ','jwst_f335m  ','jwst_f356w  ','jwst_f360m  ',
     6  'jwst_f410m  ','jwst_f430m  ','jwst_f444w  ','jwst_f460m  ',
     7  'jwst_f480m  ','skym_u      ','skym_v      ','skym_g      ',
     8  'skym_r      ','skym_i      ','skym_z      ','Hp_tycho    ',
     9  'B_tycho     ','V_tycho     ','PanS_g      ','PanS_r      ',
     *  'PanS_i      ','PanS_z      ','PanS_y      ','PanS_w      ',
     *  'PanS_o      ','GaiaG_BP    ','GaiaG       ','GaiaG_RP    '/
 1000 format(i3,14x,i2,14x,i2,14x,i2,a26)
 1001 format(13f6.0)
 1002 format(15f4.1)
 1003 format(15i4)
 1004 format(i3,' Teffs',8x,i2,' log gs',7x,i2,' [Fe/H]s',6x,i2,a26)
 1005 format(a50,3x,a8)
 1006 format(10f8.4) 
 1007 format(' files inputbc_r00.data, inputbc_r12.data, etc., have',
     1  1x,'been created for',5(/,10x,a8,': ',a12))
      do 4 i=1,nfil1
      filter(i)=filtr1(i)
    4 continue
      do 5 i=1,nfil2
      filter(i+41)=filtr2(i)
    5 continue
c *** the input file "selectbc.data" is assigned to unit 9 
      fnme(1:13)='selectbc.data'
      call iofile(9,'in',fnme,13)
      read(9,*) idumy
      fnme=fzero
      fnme(1:8)=outtab(1:8)
      fnme(12:16)=ftype(1:5)
c *** the 7 output files (for different reddenings) are assigned to units 10-16 
      nf1=9
      do 10 i=1,7
      nf1=nf1+1
      rednme=ebmv(i)
      fnme(9:11)=rednme(1:3)
      call iofile(nf1,'oo',fnme,16)
   10 continue
c *** the number of filters for which BC tables are requested is read (unit 9)
      read(9,*) numbc
      if(numbc.gt.0.and.numbc.lt.6) go to 11
      stop 'the number of filters must be > 0 and < 6'
c *** the main DO-loop is over the number of individual system/filter names
   11 do 40 ll=1,numbc
c *** integers selecting the photometric system and filter are read (unit 9)
      read(9,*) isys,ifil
      call chkinput(isys,ifil)
c *** the directory for the desired variation of [alpha/Fe] with [Fe/H] is set
      alnme=alpha(ialf)
      if(isys.le.7) go to 12
      alnme=alpha(4+ialf)
   12 psnme=system(isys)
      finme=filter(ifil)
      xps(ll)=psnme
      xfi(ll)=finme
c *** set up the directory portion of the header.data file names
      fnme=fzero
      fnme(1:9)=alnme(1:9)
      fnme(10:10)=slash(1:1)
      fnme(14:14)=slash(1:1)
c *** the header information is written to units 10-16 only when ll=1
      if(ll.ne.1) go to 20
      fnme(15:25)=hdr(1:11)
      nf1=9
      nf2=19
      do 15 j=1,7
      nf1=nf1+1
      nf2=nf2+1      
c *** modify, as appropriate, the reddening part of the header.data file name
      rednme=ebmv(j)
      fnme(11:13)=rednme(1:3)
c *** open the file header.data file for each E(B-V) value in turn 
      call iofile(nf2,'in',fnme,25)
c *** copy over the header.data information into the unit 10-16 data files
      read(nf2,1000) nt,ng,nfe,nbc,redden
      read(nf2,1001) (teff(i),i=1,nt)
      read(nf2,1002) (gv(i),i=1,ng)
      read(nf2,1003) (itmm(i),i=1,ng)
c *** nbc is redefined to be numbc
      nbc=numbc
      write(nf1,1004) nt,ng,nfe,nbc,redden
      write(nf1,1001) (teff(i),i=1,nt)
      write(nf1,1002) (gv(i),i=1,ng)
      write(nf1,1003) (itmm(i),i=1,ng)
   15 continue
c *** select the appropriate directories for the photometric system and filter
   20 call str_trim(8,psnme,nchp)
      fnme(15:14+nchp)=psnme(1:nchp)
      fnme(14+nchp+1:14+nchp+1)=slash
      call str_trim(12,finme,nch)
      fnme(15+nchp+1:15+nchp+nch)=finme(1:nch)
      fnme(16+nchp+nch:15+nchp+nch+5)=ftype(1:5)
      call str_trim(40,fnme,nch)
      nf1=9
      nf2=29
      do 35 j=1,7
      nf1=nf1+1
      nf2=nf2+1
      rednme=ebmv(j)
      fnme(11:13)=rednme(1:3)
c *** files containing the selected BC data are opened
      call iofile(nf2,'in',fnme,nch)
c *** the BC data are copied over to the output files assigned to units 13-16
      do 30 kk=1,nfe
      read(nf2,1005) labtab
      write(nf1,1005) labtab,psnme
      do 25 jj=1,ng
      itm=itmm(jj)
      read(nf2,1006) (bc(i),i=1,itm)
      write(nf1,1006) (bc(i),i=1,itm)
   25 continue
   30 continue
   35 continue
   40 continue
c *** all of the input and output files are kept on disk
      nf1=9
      nf2=19
      nf3=29
      close(unit=9,status='keep')
      do 65 i=1,7
      nf1=nf1+1
      nf2=nf2+1
      nf3=nf3+1
      close(unit=nf1,status='keep')
      close(unit=nf2,status='keep')
      close(unit=nf3,status='keep')
   65 continue
      if(ialf.ne.1) go to 70
      write(6,1007) (xps(i),xfi(i),i=1,numbc)
   70 call interpbcs(ialf,xred,xfi)
      return
      end
      subroutine chkinput(isys,ifil)
 1000 format(' the selected photometric system: integer',i3,/,6x,
     1  'and/or the selected filter: integer',i3,' are not in the menu')
      if(isys.ne.1) go to 10
      if(ifil.ge.1.and.ifil.le.3) go to 75
    5 write(6,1000) isys,ifil
      stop
   10 if(isys.lt.2) go to 5
      if(isys.gt.4) go to 15
      if(ifil.ge.4.and.ifil.le.29) go to 75
      go to 5
   15 if(isys.ne.5) go to 20
      if(ifil.ge.30.and.ifil.le.34) go to 75
      go to 5
   20 if(isys.ne.6) go to 25
      if(ifil.ge.35.and.ifil.le.39) go to 75
      go to 5
   25 if(isys.ne.7) go to 30
      if(ifil.ge.36.and.ifil.le.41) go to 75
      go to 5
   30 if(isys.gt.10) go to 35
      if(ifil.ge.42.and.ifil.le.47) go to 75
      go to 5
   35 if(isys.gt.19) go to 40
      if(ifil.ge.48.and.ifil.le.69) go to 75
      go to 5
   40 if(isys.ne.20) go to 45
      if(ifil.ge.70.and.ifil.le.75) go to 75
      go to 5
   45 if(isys.ne.21) go to 50
      if(ifil.ge.76.and.ifil.le.78) go to 75
      go to 5
   50 if(isys.ne.22) go to 55
      if(ifil.ge.79.and.ifil.le.85) go to 75
      go to 5
   55 if(isys.gt.27) go to 5
      if(ifil.ge.86.and.ifil.le.88) go to 75
      go to 5
   75 return
      end
      subroutine str_trim(nl,str,nch)
      character*(*) str
      nch=nl
      do while (str(nch:nch).eq.' '.and.nch.gt.0)
        nch=nch-1
      end do
      return 
      end
      subroutine iofile(nunit,fstat,fnme,nch)
      integer*4 nunit,nch
      character*42 fnme
      character*2 fstat
 1000 format(' input file assigned to unit',i3,' does not exist',/,a42)
      if(fstat.eq.'in') go to 10
c *** the named file is apparently an output file
      open(unit=nunit,iostat=ierr,file=fnme(1:nch),status='new')
      if(ierr.eq.0) go to 15
c *** if the named file already exists it is deleted and opened as a
c *** new file; i.e., the existing file will be overwritten
      open(unit=nunit,file=fnme(1:nch),status='old')
      close(unit=nunit,status='delete')
      open(unit=nunit,file=fnme(1:nch),status='new')
      go to 15
c *** the named file is an input file 
   10 open(unit=nunit,iostat=ierr,file=fnme(1:nch),status='old')
      if(ierr.eq.0) go to 15
c *** the named input file apparently does not exist: execution terminated
      write(6,1000) nunit,fnme
      stop 'specified input file does not exist'
   15 return
      end
      subroutine interpbcs(ialf,ebmv,xfi)
C *** This program fits Akima splines to the BC tables for E(B-V) = 0.00,
C *** 0.12, 0.24, 0.36, 0.48, 0.60, and 0.72 and and uses these splines to
C *** produce a BC table (with the name "bc_std.data", "bc_p04.data",
C *** "bc_p00.data", or "bc_m04.data" if ialf = 1, 2, 3, or 4, respectively)
C *** for any desired (input) value of E(B-V) within the above range.  Note
C *** that the output file ("bc_***.data") is rewritten each time this
C *** program is executed.
C
      real*8 c(4,7),x(7),y(7),red,xin,yout,deriv
      real*4 b(7,40),ys(40),teff(31),fet(15),aft(15),gv(13)
      character*42 fnme
      character*21 namred
      character*12 xfi(5)
      character*11 fout(4),fiout,namsys
      character*6 name
      integer*4 itmm(13)
      data fout/'bc_std.data','bc_p04.data','bc_p00.data','bc_m04.data'/
 1002 format(15i4) 
 1003 format(' Tables of Bolometric Corrections Derived from the',
     1  1x,'2008 MARCS model atmospheres',/,i3,' Teffs',8x,i2,
     2  ' log gs',7x,i2,' [Fe/H]s',6x,i2,a21,f6.3)
 1004 format(13f6.0)
 1005 format(15f4.1)
 1006 format(' [Fe/H] = ',f5.2,' [alpha/Fe] =',f5.2,' BC value: ',a6,
     1  a11,1x,a12)
 1007 format(10f8.4)
 1008 format(i3,14x,i2,14x,i2,14x,i2,a21)
 1009 format(10x,f5.2,13x,f5.2,11x,a6,a11)
 1010 format(/,1x,a11,' has been created assuming E(B-V) =',f6.3)
c *** read the input reddening (= xx) and define the independent x(i) values
      fnme(1:13)='selectbc.data'
      call iofile(9,'in',fnme,13)
      read(9,*) idumy
      fiout=fout(ialf)
      xin=dble(ebmv)
      red=-0.12d0
      do 5 i=1,7
      red=red+0.12d0
      x(i)=red
    5 continue
c *** assign the files to be interpolated
c     open(unit=10,err=45,file='bctable.data',status='new')
      open(unit=10,iostat=ierr,file=fiout,status='new')
      if(ierr.eq.0) go to 10
      open(unit=10,file=fiout,status='old')
      close(unit=10,status='delete')
      open(unit=10,file=fiout,status='new')
   10 open(unit=30,err=60,file='inputbc_r00.data',status='old')
      open(unit=31,err=60,file='inputbc_r12.data',status='old')
      open(unit=32,err=60,file='inputbc_r24.data',status='old')
      open(unit=33,err=60,file='inputbc_r36.data',status='old')
      open(unit=34,err=60,file='inputbc_r48.data',status='old')
      open(unit=35,err=60,file='inputbc_r60.data',status='old')
      open(unit=36,err=60,file='inputbc_r72.data',status='old')
c *** n is the number of input BC files (normally 7)
      n=7
c *** read the header lines from the 7 input files
      nf1=29
      do 15 i=1,7
      nf1=nf1+1
      read(nf1,1008) nt,ng,nfe,nbc,namred
      read(nf1,1004) (teff(j),j=1,nt)
      read(nf1,1005) (gv(j),j=1,ng)
      read(nf1,1002) (itmm(j),j=1,ng)
   15 continue
c *** write the header lines for the output file
      write(10,1003) nt,ng,nfe,nbc,namred,ebmv
      write(10,1004) (teff(j),j=1,nt)
      write(10,1005) (gv(j),j=1,ng)
      write(10,1002) (itmm(j),j=1,ng)
c *** main do-loops for the spline interpolations
      do 50 ll=1,nbc
      do 45 kk=1,nfe
c *** read the header line for each of the 7 input BC tables
      nf1=29
      do 20 i=1,7
      nf1=nf1+1
      read(nf1,1009) feh,afe,name,namsys
   20 continue
c *** write the corresponding header line for the output file
      write(10,1006) feh,afe,name,namsys,xfi(ll)
c *** read the BC data for all Teff values at a given gravity
      do 40 jj=1,ng
      itm=itmm(jj)
      nf1=29
      do 25 i=1,7
      nf1=nf1+1
      read(nf1,1007) (b(i,j),j=1,itm)
   25 continue
c *** define the array of independent variables, y(i)
      do 35 j=1,itm
      do 30 i=1,7
      y(i)=dble(b(i,j))
   30 continue
c *** determine the spline-fit coefficients
      call akm_cspl(n,x,y,c)
c *** evaluate the splines at the input xx(i) = input reddening values
      call akm_eval(n,x,c,xin,yout,deriv)
      ys(j)=sngl(yout)
   35 continue
c *** write the interpolated BC values in the output file
      write(10,1007) (ys(j),j=1,itm)
   40 continue
   45 continue
   50 continue
      write(6,1010) fiout,ebmv
      close(unit=9,status='keep')
      close(unit=10,status='keep')
      close(unit=30,status='keep')
      close(unit=31,status='keep')
      close(unit=32,status='keep')
      close(unit=33,status='keep')
      close(unit=34,status='keep')
      close(unit=35,status='keep')
      close(unit=36,status='keep')
      return
   60 stop ' problem with one of the input file assignments'
      end
      SUBROUTINE AKM_CSPL (n,x,f,c)
c
c  DESCRIPTION
c  Subroutine akm_cspl is based on the method of interpolation described
c  by Hiroshi Akima (1970 Journal of the Association for Computing Machinery,
c  Vol. 17, pp. 580-602) and also by Carl de Boor (1978 Applied Mathematical
c  Sciences, Vol. 27, "A Practical Guide to Splines").
c
c  This routine computes the Akima spline interpolant with the "not-a-knot
c  end conditions which are useful when there is no information about the
c  derivatives at the endpoints. In this case P(1)=P(2) and P(n-2)=P(n-1);
c  in other words, x(2) and x(n-1) are not knots even though they define
c  the endpoints of a polynomial piece --- the spline must pass through
c  these two points. Thus, it is a requirement that f''' be continuous
c  across x(2) and x(n-1). 
c
c  DIMENSIONS
c  The internal work arrays dd, w, and s are restricted to 2000 elements 
c  by the parameter nmx; consequently the data arrays x and f are also
c  limited to 2000 elements. The spline coefficients are contained in
c  the array c, a 4x2000 2 dimensional array.
c
c  CALL PARAMETERS
c     n      - number of data pairs
c     x(n)   - array, independent variable
c     f(n)   - array, dependent variable
c     c(4,n) - array, spline coefficients
c
      implicit double precision (a-h, o-z)
      parameter (nmx = 3000)
      parameter (zero = 0.0d0, two = 2.0d0, three = 3.0d0)
      dimension x(n),f(n),c(4,n) 
      dimension dd(nmx),w(nmx),s(nmx)
c
      nm1 = n - 1
      nm2 = n - 2
      nm3 = n - 3
c
c  Create two additional points outside the region of the spline by fitting
c  a quadratic to the three adjacent data points. This is a special feature
c  of the Akima spline.
c
      call akm_ends (x(1),x(2),x(3),f(1),f(2),f(3),f0,fm1)
      call akm_ends (x(n),x(nm1),x(nm2),f(n),f(nm1),f(nm2),fnp1,fnp2)
c
c  Compute the divided differences
c
      ddm1 = (f0 - fm1)/(x(2) - x(1))
      dd0 = (f(1) - f0)/(x(3) - x(2))
      do i = 1, nm1
        ip1 = i + 1
        dd(i) = (f(ip1) - f(i))/(x(ip1) - x(i))
      end do
      dd(n) = (f(n) - fnp1)/(x(nm2) - x(nm1))
      ddnp1 = (fnp1 - fnp2)/(x(nm1) - x(n))
c
c  Compute the Akima weights
c
      w0 = abs (dd0 - ddm1)
      w(1) = abs (dd(1) - dd0)
      do i = 2, nm1
        im1 = i - 1
        w(i) = abs(dd(i) - dd(im1))
      end do
      w(n) = abs (dd(n) - dd(nm1))
      wnp1 = abs (ddnp1 - dd(n))
c
c  Compute Akima slopes at interior knots
c
      if (w(2).eq.zero.and.w0.eq.zero) then
        s(1) = 5.0d-1*(dd(1) + dd0)
      else
        s(1) = (w(2)*dd0 + w0*dd(1))/(w0 + w(2))
      end if
      do i = 2, nm1
        im1 = i - 1
        ip1 = i + 1
        if (w(ip1).eq.zero.and.w(im1).eq.zero) then
          s(i) = 5.0d-1*(dd(i) + dd(im1))
        else
          s(i) = (w(ip1)*dd(im1) + w(im1)*dd(i))/(w(im1) + w(ip1))
        end if
      end do
      if (wnp1.eq.zero.and.w(nm1).eq.zero) then
        s(n) = 5.0d-1*(dd(n) + dd(nm1))
      else
        s(n) = (wnp1*dd(nm1) + w(nm1)*dd(n))/(w(nm1) + wnp1)
      end if
c
c  Spline coefficients for all the polynomial pieces
c
      do i = 1, nm1
        ip1 = i + 1
        dx = x(ip1) - x(i)
        c(1,i) = f(i)
        c(2,i) = s(i)
        c(3,i) = (three*dd(i) - two*s(i) - s(ip1))/dx
        c(4,i) = (s(ip1) + s(i) - two*dd(i))/dx/dx
      end do
      return
      end
      SUBROUTINE AKM_ENDS(x1,x2,x3,f1,f2,f3,f00,f01)
c      
c  DESCRIPTION
c  Subroutine akm_ends is based on the method of interpolation described
c  by Hiroshi Akima (1970 Journal of the Association for Computing Machinery,
c  Vol. 17, pp. 580-602).
c
c  This routine is required by the routine AKM_CSPL. It sets up the two 
c  additional external points required by the Akima method of constructing
c  a spline.
c
c  CALL PARAMETERS
c  (x1,f1), (x2,f2) and (x3,f3) are the known data pairs at the ends of the
c  region of interpolation. f00 and f01 are the computed values of the
c  interpolating polynomial external to the region of interpolation.
c
      implicit double precision (a-h,o-z)
      g0 = f1
      df21 = f2 - f1
      df31 = f3 - f1
      dx21 = x2 - x1
      dx31 = x3 - x1
      dx32 = x3 - x2
      den = dx21*dx32*dx31
      g1 = (df21*dx31*dx31 - df31*dx21*dx21)/den
      g2 = (df31*dx21 - df21*dx31)/den
      f00 = g0 - g1*dx32 + g2*dx32*dx32
      f01 = g0 - g1*dx31 + g2*dx31*dx31
      return
      end
      SUBROUTINE AKM_EVAL (n,x,c,xp,fp,dfp)
c     
c  DESCRIPTION
c  This routine is used to evaluate the spline defined by the
c  routine AKM_CSPL and its derivative at an interpolating point.
c
c  CALL PARAMETERS
c     n      - number of data pairs
c     x(n)   - array, independent variable
c     c(4,n) - array, spline coefficients
c     xp     - interpolating point
c     fp     - value of the spline at xp
c     dfp    - first derivative of the spline at xp   
c     
      implicit double precision (a-h, o-z)
      dimension x(n), c(4,n)
      logical more
c
c  Statement function defines cubic spline interpolation
c
      csval(dx,a,b,e,d) = a +dx*(b + dx*(e + dx*d))
c
c  Statement function defines cubic spline derivative
c
      csder(dx,b,e,d) = b + dx*(2.0d0*e + dx*3.0d0*d)
c
c  Check direction of spline
c
      if (x(n).gt.x(1)) then
c
c  Check to see that x(1) < xp < x(n)
c
      if (xp.lt.x(1)) then
        write(*,'('' ***Test point is less than x(1)***'')')
        stop
      else if (xp.gt.x(n)) then
        write(*,'('' ***Test point is greater than x(n)***'')')
        stop
      end if
c
c  Find the polynomial piece containing the test point
c
      i = 2
      more = .true.
      do while (more)
        if (xp.lt.x(i).or.i.eq.n) then
          more = .false.
          nl = i - 1
        else 
          i = i + 1
        end if
      end do
      else
c
c  The spline is backwards
c  Check to see that x(n) < xp < x(1)
c
      if (xp.gt.x(1)) then
        write(*,'('' ***Test point is greater than x(1)***'')')
        stop
      else if (xp.lt.x(n)) then
        write(*,'('' ***Test point is less than x(n)***'')')
        stop
      end if
c
c  Find the polynomial piece containing the test point
c
      i = 2
      more = .true.
      do while (more)
        if (xp.gt.x(i).or.i.eq.n) then
          more = .false.
          nl = i - 1
        else 
          i = i + 1
        end if
      end do
      end if
c
c  Evaluate spline at the test point xp
c
      dx = xp - x(nl)
      fp = csval(dx,c(1,nl),c(2,nl),c(3,nl),c(4,nl))
c
c  Evaluate the derivative at the test point
c
      dfp = csder(dx,c(2,nl),c(3,nl),c(4,nl))
      return
      end        
