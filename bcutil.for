      Subroutine getbc_std(mode,fe,gv,teff,ebv,bc,fil,sys,nbc)
C ----------------------------------------------------------------------
C *** Note that bc and fil are arrays: the program which calls getbc_std
C     must have the specifications "real*4 bc(5)" and "character*6 fil(5)"
C     in it to avoid errors.
C *** This routine interpolates in the tables of bolometric corrections
C     derived from MARCS model atmospheres for up to 5 filter bandpasses
C     of choice.  The standard variation of [alpha/Fe] with [Fe/H] is
C     assumed; i.e., [alpha/Fe] = +0.4 at [Fe/H] <= -1.0, 0.0 at [Fe/H]
C     >= 0.0, and [alpha/Fe] = -0.4[Fe/H] at intermediate iron abundances.
C     The input tables, in a file named `bc_std.data', are assigned to
C     unit 30.
C
C *** If MODE=0, the BC tables are read and interpolated to the desired
C     value of [Fe/H] (= fe in the argument list).  MODE must be set to
C     zero the first time that this subroutine is called.  A message is
C     sent to unit 6 (which is assumed to be attached to the monitor) to
C     specify the [Fe/H] value and the magnitudes (e.g., B, I, F606W) for
C     which bolometric corrections are provided in the input tables.
C
C *** If MODE=1, interpolations are carried out to obtain the bolometric
C     corrections at the input values of log g and log Te (gv and teff in
C     the argument list) and the specified [Fe/H] value.  It is only for
C     this value of MODE that interpolations are performed to obtain BC
C     information.
C
C *** If MODE=-1, the input tables (which were previously read by setting
C     MODE=0) are re-interpolated for a new value of [Fe/H].  If desired,
C     the subroutine may be modified (after the statement `45 continue')
C     so that a message is sent to unit 6 each time an [Fe/H] interpolation
C     is performed.
C
C *** The interpolation code will issue a STOP code if the input [Fe/H]
C     or log Teff values are outside the ranges assumed in "bc_std.data"
C     or the latter file does not exist.  If the input value of log g is
C     outside the tabulated range, the output BC values are all set to a
C     value of "9.9999".
C -------------------------------------------------------------------7--
      real*4 ttt(8),bc(5)
      character*29 alfe(15)
      character*10 namfe
      character*8 sys(5)
      character*6 fil(5)
      common /bcstd/ a(31,13,5,15),b(31,13,5),c(4,5),te(31),t(31),g(13),
     1  feh(15),ttm(8),r(4),af(4),ag(4),at(4),maxt(13),nt,ng,nfe,ndx
      data ttt/4250.,5000.,5500.,6000.,6500.,6750.,7750.,8000./
 1000 format(/,i3,14x,i2,14x,i2,14x,i2,21x,f6.3)
 1001 format(13f6.0)
 1002 format(20f4.1)
 1003 format(10f8.4)
 1004 format(' bc_std BCs for E(B-V) =',f6.3,' interpolated to [Fe/H] ='
     1  ,f6.2,/,16x,'filters: ',5(a6,2x))
 1005 format(' Teff =',f8.2,' is outside the range of the input tables')
 1006 format(a10,f5.2,a29,a6,3x,a8)
 1007 format(' fe =',f5.2,' is outside the range of the input tables')
 1008 format(20i4)
 1009 format(' log g =',f7.4,' is outside the range of the input',1x,
     1  'tables: no interpolations are performed - bc(i) = 9.9999')
      if(mode.gt.0) go to 50
      if(mode.lt.0) go to 25
C *** The input BC tables are read only if MODE=0.
      open(unit=30,iostat=ierr,file='bc_std.data',status='old')
      if(ierr.ne.0) go to 150
      read(30,1000) nt,ng,nfe,ndx,ebv
      read(30,1001) (te(i),i=1,nt)
      read(30,1002) (g(i),i=1,ng)
      read(30,1008) (maxt(i),i=1,ng)
      nbc=ndx
      do 5 i=1,8
      ttm(i)=alog10(ttt(i))
    5 continue
      do 10 i=1,nt
      t(i)=alog10(te(i))
   10 continue
      do 20 k=1,ndx
      do 20 l=1,nfe
      read(30,1006) namfe,feh(l),alfe(l),fil(k),sys(k)
      do 15 j=1,ng
      ntot=maxt(j)
      read(30,1003) (a(i,j,k,l),i=1,ntot)
   15 continue
   20 continue
      close(unit=30,status='keep')
C *** When MODE=0 OR MODE=-1, the input BC table is interpolated to
C *** create a table of BC values for the input value of [Fe/H) (=fe).
   25 if(fe.lt.feh(1).or.fe.gt.feh(nfe)) go to 35
      nfem1=nfe-1
      do 30 m=3,nfem1
      l=m-2
      if(fe.le.feh(m)) go to 40
   30 continue
      go to 40
   35 write(6,1007) fe
      stop ' input [Fe/H] value is outside the range of the tables'
   40 r(1)=feh(l)
      r(2)=feh(l+1)
      r(3)=feh(l+2)
      r(4)=feh(l+3)
      call lgr4(r,af,fe)
      do 45 k=1,ndx
      do 45 j=1,ng
      ntot=maxt(j)
      do 45 i=1,ntot
      b(i,j,k)=af(1)*a(i,j,k,l)+af(2)*a(i,j,k,l+1)+af(3)*a(i,j,k,l+2)+
     1  af(4)*a(i,j,k,l+3)
   45 continue
C *** If desired, the following statement may be commented out so that
C *** the write statement is executed after each [Fe/H] interpolation.
      if(mode.ne.0) go to 165
      write(6,1004) ebv,fe,(fil(k),k=1,ndx)
C *** When MODE=0 or MODE=-1, control returns to the calling program
C *** once the [Fe/H] interpolations have been completed; i.e., no
C *** interpolations are performed within the tables.
      go to 165
C *** MODE=1: interpolations are performed for the input values of
C *** log g (= gv) and log Teff (= teff).
   50 nbc=ndx
      temp=10.**teff
      if(abs(teff-t(nt)).ge.1.e-5) go to 55
      teff=t(nt)
   55 if(abs(teff-t(1)).ge.1.e-5) go to 60 
      teff=t(1)
   60 if(teff.le.t(nt).and.teff.ge.t(1)) go to 70
   65 write(6,1005) temp
      stop ' range of the tables exceeded'
   70 if(gv.lt.-0.5.or.gv.gt.5.5) go to 155
C *** Note that ngmin=1 only if the input value of gv is < 0.0.
      ngmin=2
      if(gv.ge.0.0) go to 75
      ngmin=1
C *** This loop determines the index ngm where ng(ngm) is the lowest
C *** value of log g for which data are available for the input Teff.
   75 do 80 i=ngmin,8
      ngm=i
      if(teff.le.ttm(i)) go to 85
   80 continue
   85 ngl=ngm+2
C *** Note that ngmax=13 only if the input value of gv is > 5.0.
      ngmax=13
      if(gv.gt.5.0) go to 90
      ngmax=12
   90 ngm=ngmax-1
C *** This loop determines the index mg where g(mg) is the lowest 
C *** of the grid values of log g that are used in the interpolations.
      do 95 i=ngl,ngm
      mg=i-2
      if(gv.le.g(i)) go to 105
   95 continue
      mg=ngmax-3
  105 mmx=maxt(mg)
C *** This section increases the value of mg if the maximum value of
C *** Teff at the initial value of mg is less than the input Teff.
      if(teff.le.t(mmx)) go to 110
      mg=mg+1
      mmx=maxt(mg)
C *** This section tests whether the maximum value of Teff at the
C *** revised value of mg is >= the input Teff value.  If this condition
C *** is not satisfied, the execution terminates.
      if(teff.gt.t(mmx)) go to 65
  110 r(1)=g(mg)
      r(2)=g(mg+1)
      r(3)=g(mg+2)
      r(4)=g(mg+3)
C *** The Lagrangian coefficients for the log g interpolation are computed.
      call lgr4(r,ag,gv)
      ntm=nt-1
C *** This loop determines the value of mt where t(mt) is the lowest of
C *** the grid values of log Teff that are used in the interpolations.
      do 115 i=3,ntm
      mt=i-2
      if(teff.le.t(i)) go to 120
  115 continue
      mt=nt-3
C *** This section checks that the input T is <= the maximum value for
C *** which tabulated data are available at the lowest value of mg 
C *** in the interpolations.
  120 limt=maxt(mg)
      if(gv.le.5.0) go to 125
      limt=19
      if(teff.gt.t(limt)) go to 65
  125 if((mt+3).le.limt) go to 130
      mt=limt-3
  130 r(1)=t(mt)
      r(2)=t(mt+1)
      r(3)=t(mt+2)
      r(4)=t(mt+3)
C *** The Lagrangian coefficients for the log Teff interpolations are computed.
      call lgr4(r,at,teff)
      l=mg-1
C *** The desired bolometric corrections are calculated.
      do 140 i=1,4
      l=l+1
      do 135 k=1,ndx
      c(i,k)=at(1)*b(mt,l,k)+at(2)*b(mt+1,l,k)+at(3)*b(mt+2,l,k)+
     1   at(4)*b(mt+3,l,k)
  135 continue
  140 continue
      do 145 i=1,ndx
      bc(i)=ag(1)*c(1,i)+ag(2)*c(2,i)+ag(3)*c(3,i)+ag(4)*c(4,i)
  145 continue
      go to 165
  150 stop ' input file bc_std.data does not exist'
  155 do 160 i=1,ndx
      bc(i)=9.9999
  160 continue
      write(6,1009) gv
  165 return
      end
      Subroutine getbc_p04(mode,fe,gv,teff,ebv,bc,fil,sys,nbc)
C ----------------------------------------------------------------------
C *** Note that bc and fil are arrays: the program which calls getbc_p04
C     must have the specifications "real*4 bc(5)" and "character*6 fil(5)"
C     in it to avoid errors.
C *** This routine interpolates in the tables of bolometric corrections
C     derived from MARCS model atmospheres for up to 5 filter bandpasses
C     of choice.  The input tables, in a file named `bc_p04.data' which is
C     assigned to unit 30, assume that [alpha/Fe] = +0.4 for all [Fe/H]
C     values.
C
C *** If MODE=0, the BC tables are read and interpolated to the desired
C     value of [Fe/H] (= fe in the argument list).  MODE must be set to
C     zero the first time that this subroutine is called.  A message is
C     sent to unit 6 (which is assumed to be attached to the monitor) to
C     specify the [Fe/H] value and the magnitudes (e.g., B, I, F606W) for
C     which bolometric corrections are provided in the input tables.  
C
C *** If MODE=1, interpolations are carried out to obtain the bolometric
C     corrections at the input values of log g and log Te (gv and teff in
C     the argument list) and the specified [Fe/H] value.  It is only for
C     this value of MODE that interpolations are performed to obtain BC
C     information.
C
C *** If MODE=-1, the input tables (which were previously read by setting
C     MODE=0) are re-interpolated for a new value of [Fe/H].  If desired,
C     the subroutine may be modified (after the statement `45 continue')
C     so that a message is sent to unit 6 each time an [Fe/H] interpolation
C     is performed.
C
C *** The interpolation code will issue a STOP code if the input [Fe/H]
C     or log Teff values are outside the ranges assumed in "bc_p04.data"
C     or the latter file does not exist.  If the input value of log g is
C     outside the tabulated range, the output BC values are all set to a
C     value of "9.9999".
C -------------------------------------------------------------------7--
      real*4 ttt(8),bc(5)
      character*29 alfe(15)
      character*10 namfe
      character*8 sys(5)
      character*6 fil(5)
      common /bcp04/ a(31,13,5,15),b(31,13,5),c(4,5),te(31),t(31),g(13),
     1  feh(15),ttm(8),r(4),af(4),ag(4),at(4),maxt(13),nt,ng,nfe,ndx
      data ttt/4250.,5000.,5250.,5500.,5750.,6250.,7750.,8000./
 1000 format(/,i3,14x,i2,14x,i2,14x,i2,21x,f6.3)
 1001 format(13f6.0)
 1002 format(20f4.1)
 1003 format(10f8.4)
 1004 format(' bc_p04 BCs for E(B-V) =',f6.3,' interpolated to [Fe/H] ='
     1  ,f6.2,/,16x,'filters: ',5(a6,2x))
 1005 format(' Teff =',f8.2,' is outside the range of the input tables')
 1006 format(a10,f5.2,a29,a6,3x,a8)
 1007 format(' fe =',f5.2,' is outside the range of the input tables')
 1008 format(20i4)
 1009 format(' log g =',f7.4,' is outside the range of the input',1x,
     1  'tables: no interpolations are performed - bc(i) = 9.9999')
      if(mode.gt.0) go to 50
      if(mode.lt.0) go to 25
C *** The input BC tables are read only if MODE=0.
      open(unit=30,iostat=ierr,file='bc_p04.data',status='old')
      if(ierr.ne.0) go to 150
      read(30,1000) nt,ng,nfe,ndx,ebv
      ng=13
      read(30,1001) (te(i),i=1,nt)
      read(30,1002) (g(i),i=2,ng)
      read(30,1008) (maxt(i),i=2,ng)
      nbc=ndx
      do 5 i=1,8
      ttm(i)=alog10(ttt(i))
    5 continue
      do 10 i=1,nt
      t(i)=alog10(te(i))
   10 continue
      do 20 k=1,ndx
      do 20 l=1,nfe
      read(30,1006) namfe,feh(l),alfe(l),fil(k),sys(k)
      do 15 j=2,ng
      ntot=maxt(j)
      read(30,1003) (a(i,j,k,l),i=1,ntot)
   15 continue
   20 continue
C *** linear extrapolation of tables to log g = -0.5 (Teff from 2600 to 4250 K)
      g(1)=-0.5
      maxt(1)=16
      do 22 i=1,16
      a(i,1,k,l)=2.*(a(i,2,k,l)-a(i,3,k,l))
   22 continue
      close(unit=30,status='keep')
C *** When MODE=0 OR MODE=-1, the input BC table is interpolated to
C *** create a table of BC values for the input value of [Fe/H) (=fe).
   25 if(fe.lt.feh(1).or.fe.gt.feh(nfe)) go to 35
      nfem1=nfe-1
      do 30 m=3,nfem1
      l=m-2
      if(fe.le.feh(m)) go to 40
   30 continue
      go to 40
   35 write(6,1007) fe
      stop ' input [Fe/H] value is outside the range of the tables'
   40 r(1)=feh(l)
      r(2)=feh(l+1)
      r(3)=feh(l+2)
      r(4)=feh(l+3)
      call lgr4(r,af,fe)
      do 45 k=1,ndx
      do 45 j=1,ng
      ntot=maxt(j)
      do 45 i=1,ntot
      b(i,j,k)=af(1)*a(i,j,k,l)+af(2)*a(i,j,k,l+1)+af(3)*a(i,j,k,l+2)+
     1  af(4)*a(i,j,k,l+3)
   45 continue
C *** If desired, the following statement may be commented out so that
C *** the write statement is executed after each [Fe/H] interpolation.
      if(mode.ne.0) go to 165
      write(6,1004) ebv,fe,(fil(k),k=1,ndx)
C *** When MODE=0 or MODE=-1, control returns to the calling program
C *** once the [Fe/H] interpolations have been completed; i.e., no
C *** interpolations are performed within the tables.
      go to 165
C *** MODE=1: interpolations are performed for the input values of
C *** log g (= gv) and log Teff (= teff).
   50 nbc=ndx
      temp=10.**teff
      if(abs(teff-t(nt)).ge.1.e-5) go to 55
      teff=t(nt)
   55 if(abs(teff-t(1)).ge.1.e-5) go to 60 
      teff=t(1)
   60 if(teff.le.t(nt).and.teff.ge.t(1)) go to 70
   65 write(6,1005) temp
      stop ' range of the tables exceeded'
   70 if(gv.lt.-0.5.or.gv.gt.5.5) go to 155
C *** Note that ngmin=1 only if the input value of gv is < 0.0, and
C *** that ngmin=2 only if gv < 0.5 but >= 0.0
      ngmin=3
      if(gv.ge.0.5) go to 75
      ngmin=2
      if(gv.ge.0.0) go to 75
      ngmin=1
C *** This loop determines the index ngm where ng(ngm) is the lowest
C *** value of log g for which data are available for the input Teff.
   75 do 80 i=ngmin,8
      ngm=i
      if(teff.le.ttm(i)) go to 85
   80 continue
   85 ngl=ngm+2
C *** Note that ngmax=13 only if the input value of gv is > 5.0.
      ngmax=13
      if(gv.gt.5.0) go to 90
      ngmax=12
   90 ngm=ngmax-1
C *** This loop determines the index mg where g(mg) is the lowest 
C *** of the grid values of log g that are used in the interpolations.
      do 95 i=ngl,ngm
      mg=i-2
      if(gv.le.g(i)) go to 105
   95 continue
      mg=ngmax-3
  105 mmx=maxt(mg)
C *** This section increases the value of mg if the maximum value of
C *** Teff at the initial value of mg is less than the input Teff.
      if(teff.le.t(mmx)) go to 110
      mg=mg+1
      mmx=maxt(mg)
C *** This section tests whether the maximum value of Teff at the
C *** revised value of mg is >= the input Teff value.  If this condition
C *** is not satisfied, the execution terminates.
      if(teff.gt.t(mmx)) go to 65
  110 r(1)=g(mg)
      r(2)=g(mg+1)
      r(3)=g(mg+2)
      r(4)=g(mg+3)
C *** The Lagrangian coefficients for the log g interpolation are computed.
      call lgr4(r,ag,gv)
      ntm=nt-1
C *** This loop determines the value of mt where t(mt) is the lowest of
C *** the grid values of log Teff that are used in the interpolations.
      do 115 i=3,ntm
      mt=i-2
      if(teff.le.t(i)) go to 120
  115 continue
      mt=nt-3
C *** This section checks that the input T is <= the maximum value for
C *** which tabulated data are available at the lowest value of mg 
C *** in the interpolations.
  120 limt=maxt(mg)
      if(gv.le.5.0) go to 125
      limt=19
      if(teff.gt.t(limt)) go to 65
  125 if((mt+3).le.limt) go to 130
      mt=limt-3
  130 r(1)=t(mt)
      r(2)=t(mt+1)
      r(3)=t(mt+2)
      r(4)=t(mt+3)
C *** The Lagrangian coefficients for the log Teff interpolations are computed.
      call lgr4(r,at,teff)
      l=mg-1
C *** The desired bolometric corrections are calculated.
      do 140 i=1,4
      l=l+1
      do 135 k=1,ndx
      c(i,k)=at(1)*b(mt,l,k)+at(2)*b(mt+1,l,k)+at(3)*b(mt+2,l,k)+
     1   at(4)*b(mt+3,l,k)
  135 continue
  140 continue
      do 145 i=1,ndx
      bc(i)=ag(1)*c(1,i)+ag(2)*c(2,i)+ag(3)*c(3,i)+ag(4)*c(4,i)
  145 continue
      go to 165
  150 stop ' input file bc_p04.data does not exist'
  155 do 160 i=1,ndx
      bc(i)=9.9999
  160 continue
      write(6,1009) gv
  165 return
      end
      Subroutine getbc_p00(mode,fe,gv,teff,ebv,bc,fil,sys,nbc)
C ----------------------------------------------------------------------
C *** Note that bc and fil are arrays: the program which calls getbc_p00
C     must have the specifications "real*4 bc(5)" and "character*6 fil(5)"
C     in it to avoid errors.
C *** This routine interpolates in the tables of bolometric corrections
C     derived from MARCS model atmospheres for up to 5 filter bandpasses
C     of choice.  The input tables, in a file named `bc_p00.data' which is
C     assigned to unit 30, assume that [alpha/Fe] = 0.0 for all [Fe/H]
C     values.
C
C *** If MODE=0, the BC tables are read and interpolated to the desired
C     value of [Fe/H] (= fe in the argument list).  MODE must be set to
C     zero the first time that this subroutine is called.  A message is
C     sent to unit 6 (which is assumed to be attached to the monitor) to
C     specify the [Fe/H] value and the magnitudes (e.g., B, I, F606W) for
C     which bolometric corrections are provided in the input tables.  
C
C *** If MODE=1, interpolations are carried out to obtain the bolometric
C     corrections at the input values of log g and log Te (gv and teff in
C     the argument list) and the specified [Fe/H] value.  It is only for
C     this value of MODE that interpolations are performed to obtain BC
C     information.
C
C *** If MODE=-1, the input tables (which were previously read by setting
C     MODE=0) are re-interpolated for a new value of [Fe/H].  If desired,
C     the subroutine may be modified (after the statement `45 continue')
C     so that a message is sent to unit 6 each time an [Fe/H] interpolation
C     is performed.
C
C *** The interpolation code will issue a STOP code if the input [Fe/H]
C     or log Teff values are outside the ranges assumed in "bc_p00.data"
C     or the latter file does not exist.  If the input value of log g is
C     outside the tabulated range, the output BC values are all set to a
C     value of "9.9999".
C -------------------------------------------------------------------7--
      real*4 ttt(8),bc(5)
      character*29 alfe(15)
      character*10 namfe
      character*8 sys(5)
      character*6 fil(5)
      common /bcp00/ a(31,13,5,15),b(31,13,5),c(4,5),te(31),t(31),g(13),
     1  feh(15),ttm(8),r(4),af(4),ag(4),at(4),maxt(13),nt,ng,nfe,ndx
      data ttt/4250.,5000.,5250.,5500.,5750.,6250.,7750.,8000./
 1000 format(/,i3,14x,i2,14x,i2,14x,i2,21x,f6.3)
 1001 format(13f6.0)
 1002 format(20f4.1)
 1003 format(10f8.4)
 1004 format(' bc_p00 BCs for E(B-V) =',f6.3,' interpolated to [Fe/H] ='
     1  ,f6.2,/,16x,'filters: ',5(a6,2x))
 1005 format(' Teff =',f8.2,' is outside the range of the input tables')
 1006 format(a10,f5.2,a29,a6,3x,a8)
 1007 format(' fe =',f5.2,' is outside the range of the input tables')
 1008 format(20i4)
 1009 format(' log g =',f7.4,' is outside the range of the input',1x,
     1  'tables: no interpolations are performed - bc(i) = 9.9999')
      if(mode.gt.0) go to 50
      if(mode.lt.0) go to 25
C *** The input BC tables are read only if MODE=0.
      open(unit=30,iostat=ierr,file='bc_p00.data',status='old')
      if(ierr.ne.0) go to 150
      read(30,1000) nt,ng,nfe,ndx,ebv
      ng=13
      read(30,1001) (te(i),i=1,nt)
      read(30,1002) (g(i),i=2,ng)
      read(30,1008) (maxt(i),i=2,ng)
      nbc=ndx
      do 5 i=1,8
      ttm(i)=alog10(ttt(i))
    5 continue
      do 10 i=1,nt
      t(i)=alog10(te(i))
   10 continue
      do 20 k=1,ndx
      do 20 l=1,nfe
      read(30,1006) namfe,feh(l),alfe(l),fil(k),sys(k)
      do 15 j=2,ng
      ntot=maxt(j)
      read(30,1003) (a(i,j,k,l),i=1,ntot)
   15 continue
   20 continue
C *** linear extrapolation of tables to log g = -0.5 (Teff from 2600 to 4250 K)
      g(1)=-0.5
      maxt(1)=16
      do 22 i=1,16
      a(i,1,k,l)=2.*(a(i,2,k,l)-a(i,3,k,l))
   22 continue
      close(unit=30,status='keep')
C *** When MODE=0 OR MODE=-1, the input BC table is interpolated to
C *** create a table of BC values for the input value of [Fe/H) (=fe).
   25 if(fe.lt.feh(1).or.fe.gt.feh(nfe)) go to 35
      nfem1=nfe-1
      do 30 m=3,nfem1
      l=m-2
      if(fe.le.feh(m)) go to 40
   30 continue
      go to 40
   35 write(6,1007) fe
      stop ' input [Fe/H] value is outside the range of the tables'
   40 r(1)=feh(l)
      r(2)=feh(l+1)
      r(3)=feh(l+2)
      r(4)=feh(l+3)
      call lgr4(r,af,fe)
      do 45 k=1,ndx
      do 45 j=1,ng
      ntot=maxt(j)
      do 45 i=1,ntot
      b(i,j,k)=af(1)*a(i,j,k,l)+af(2)*a(i,j,k,l+1)+af(3)*a(i,j,k,l+2)+
     1  af(4)*a(i,j,k,l+3)
   45 continue
C *** If desired, the following statement may be commented out so that
C *** the write statement is executed after each [Fe/H] interpolation.
      if(mode.ne.0) go to 165
      write(6,1004) ebv,fe,(fil(k),k=1,ndx)
C *** When MODE=0 or MODE=-1, control returns to the calling program
C *** once the [Fe/H] interpolations have been completed; i.e., no
C *** interpolations are performed within the tables.
      go to 165
C *** MODE=1: interpolations are performed for the input values of
C *** log g (= gv) and log Teff (= teff).
   50 nbc=ndx
      temp=10.**teff
      if(abs(teff-t(nt)).ge.1.e-5) go to 55
      teff=t(nt)
   55 if(abs(teff-t(1)).ge.1.e-5) go to 60 
      teff=t(1)
   60 if(teff.le.t(nt).and.teff.ge.t(1)) go to 70
   65 write(6,1005) temp
      stop ' range of the tables exceeded'
   70 if(gv.lt.-0.5.or.gv.gt.5.5) go to 155
C *** Note that ngmin=1 only if the input value of gv is < 0.0, and
C *** that ngmin=2 only if gv < 0.5 but >= 0.0
      ngmin=3
      if(gv.ge.0.5) go to 75
      ngmin=2
      if(gv.ge.0.0) go to 75
      ngmin=1
C *** This loop determines the index ngm where ng(ngm) is the lowest
C *** value of log g for which data are available for the input Teff.
   75 do 80 i=ngmin,8
      ngm=i
      if(teff.le.ttm(i)) go to 85
   80 continue
   85 ngl=ngm+2
C *** Note that ngmax=13 only if the input value of gv is > 5.0.
      ngmax=13
      if(gv.gt.5.0) go to 90
      ngmax=12
   90 ngm=ngmax-1
C *** This loop determines the index mg where g(mg) is the lowest 
C *** of the grid values of log g that are used in the interpolations.
      do 95 i=ngl,ngm
      mg=i-2
      if(gv.le.g(i)) go to 105
   95 continue
      mg=ngmax-3
  105 mmx=maxt(mg)
C *** This section increases the value of mg if the maximum value of
C *** Teff at the initial value of mg is less than the input Teff.
      if(teff.le.t(mmx)) go to 110
      mg=mg+1
      mmx=maxt(mg)
C *** This section tests whether the maximum value of Teff at the
C *** revised value of mg is >= the input Teff value.  If this condition
C *** is not satisfied, the execution terminates.
      if(teff.gt.t(mmx)) go to 65
  110 r(1)=g(mg)
      r(2)=g(mg+1)
      r(3)=g(mg+2)
      r(4)=g(mg+3)
C *** The Lagrangian coefficients for the log g interpolation are computed.
      call lgr4(r,ag,gv)
      ntm=nt-1
C *** This loop determines the value of mt where t(mt) is the lowest of
C *** the grid values of log Teff that are used in the interpolations.
      do 115 i=3,ntm
      mt=i-2
      if(teff.le.t(i)) go to 120
  115 continue
      mt=nt-3
C *** This section checks that the input T is <= the maximum value for
C *** which tabulated data are available at the lowest value of mg 
C *** in the interpolations.
  120 limt=maxt(mg)
      if(gv.le.5.0) go to 125
      limt=19
      if(teff.gt.t(limt)) go to 65
  125 if((mt+3).le.limt) go to 130
      mt=limt-3
  130 r(1)=t(mt)
      r(2)=t(mt+1)
      r(3)=t(mt+2)
      r(4)=t(mt+3)
C *** The Lagrangian coefficients for the log Teff interpolations are computed.
      call lgr4(r,at,teff)
      l=mg-1
C *** The desired bolometric corrections are calculated.
      do 140 i=1,4
      l=l+1
      do 135 k=1,ndx
      c(i,k)=at(1)*b(mt,l,k)+at(2)*b(mt+1,l,k)+at(3)*b(mt+2,l,k)+
     1   at(4)*b(mt+3,l,k)
  135 continue
  140 continue
      do 145 i=1,ndx
      bc(i)=ag(1)*c(1,i)+ag(2)*c(2,i)+ag(3)*c(3,i)+ag(4)*c(4,i)
  145 continue
      go to 165
  150 stop ' input file bc_p00.data does not exist'
  155 do 160 i=1,ndx
      bc(i)=9.9999
  160 continue
      write(6,1009) gv
  165 return
      end
      Subroutine getbc_m04(mode,fe,gv,teff,ebv,bc,fil,sys,nbc)
C ----------------------------------------------------------------------
C *** Note that bc and fil are arrays: the program which calls getbc_m04
C     must have the specifications "real*4 bc(5)" and "character*6 fil(5)"
C     in it to avoid errors.
C *** This routine interpolates in the tables of bolometric corrections
C     derived from MARCS model atmospheres for up to 5 filter bandpasses
C     of choice.  The input tables, in a file named `bc_m04.data' which is
C     assigned to unit 30, assume that [alpha/Fe] = -0.4 for all [Fe/H]
C     values.
C
C *** If MODE=0, the BC tables are read and interpolated to the desired
C     value of [Fe/H] (= fe in the argument list).  MODE must be set to
C     zero the first time that this subroutine is called.  A message is
C     sent to unit 6 (which is assumed to be attached to the monitor) to
C     specify the [Fe/H] value and the magnitudes (e.g., B, I, F606W) for
C     which bolometric corrections are provided in the input tables.  
C
C *** If MODE=1, interpolations are carried out to obtain the bolometric
C     corrections at the input values of log g and log Te (gv and teff in
C     the argument list) and the specified [Fe/H] value.  It is only for
C     this value of MODE that interpolations are performed to obtain BC
C     information.
C
C *** If MODE=-1, the input tables (which were previously read by setting
C     MODE=0) are re-interpolated for a new value of [Fe/H].  If desired,
C     the subroutine may be modified (after the statement `45 continue')
C     so that a message is sent to unit 6 each time an [Fe/H] interpolation
C     is performed.
C
C *** The interpolation code will issue a STOP code if the input [Fe/H]
C     or log Teff values are outside the ranges assumed in "bc_m04.data"
C     or the latter file does not exist.  If the input value of log g is
C     outside the tabulated range, the output BC values are all set to a
C     value of "9.9999".
C -------------------------------------------------------------------7--
      real*4 ttt(8),bc(5)
      character*29 alfe(15)
      character*10 namfe
      character*8 sys(5)
      character*6 fil(5)
      common /bcm04/ a(31,13,5,15),b(31,13,5),c(4,5),te(31),t(31),g(13),
     1  feh(15),ttm(8),r(4),af(4),ag(4),at(4),maxt(13),nt,ng,nfe,ndx
      data ttt/4250.,5000.,5250.,5500.,5750.,6250.,7750.,8000./
 1000 format(/,i3,14x,i2,14x,i2,14x,i2,21x,f6.3)
 1001 format(13f6.0)
 1002 format(20f4.1)
 1003 format(10f8.4)
 1004 format(' bc_m04 BCs for E(B-V) =',f6.3,' interpolated to [Fe/H] ='
     1  ,f6.2,/,16x,'filters: ',5(a6,2x))
 1005 format(' Teff =',f8.2,' is outside the range of the input tables')
 1006 format(a10,f5.2,a29,a6,3x,a8)
 1007 format(' fe =',f5.2,' is outside the range of the input tables')
 1008 format(20i4)
 1009 format(' log g =',f7.4,' is outside the range of the input',1x,
     1  'tables: no interpolations are performed - bc(i) = 9.9999')
      if(mode.gt.0) go to 50
      if(mode.lt.0) go to 25
C *** The input BC tables are read only if MODE=0.
      open(unit=30,iostat=ierr,file='bc_m04.data',status='old')
      if(ierr.ne.0) go to 150
      read(30,1000) nt,ng,nfe,ndx,ebv
      ng=13
      read(30,1001) (te(i),i=1,nt)
      read(30,1002) (g(i),i=2,ng)
      read(30,1008) (maxt(i),i=2,ng)
      nbc=ndx
      do 5 i=1,8
      ttm(i)=alog10(ttt(i))
    5 continue
      do 10 i=1,nt
      t(i)=alog10(te(i))
   10 continue
      do 20 k=1,ndx
      do 20 l=1,nfe
      read(30,1006) namfe,feh(l),alfe(l),fil(k),sys(k)
      do 15 j=2,ng
      ntot=maxt(j)
      read(30,1003) (a(i,j,k,l),i=1,ntot)
   15 continue
   20 continue
C *** linear extrapolation of tables to log g = -0.5 (Teff from 2600 to 4250 K)
      g(1)=-0.5
      maxt(1)=16
      do 22 i=1,16
      a(i,1,k,l)=2.*(a(i,2,k,l)-a(i,3,k,l))
   22 continue
      close(unit=30,status='keep')
C *** When MODE=0 OR MODE=-1, the input BC table is interpolated to
C *** create a table of BC values for the input value of [Fe/H) (=fe).
   25 if(fe.lt.feh(1).or.fe.gt.feh(nfe)) go to 35
      nfem1=nfe-1
      do 30 m=3,nfem1
      l=m-2
      if(fe.le.feh(m)) go to 40
   30 continue
      go to 40
   35 write(6,1007) fe
      stop ' input [Fe/H] value is outside the range of the tables'
   40 r(1)=feh(l)
      r(2)=feh(l+1)
      r(3)=feh(l+2)
      r(4)=feh(l+3)
      call lgr4(r,af,fe)
      do 45 k=1,ndx
      do 45 j=1,ng
      ntot=maxt(j)
      do 45 i=1,ntot
      b(i,j,k)=af(1)*a(i,j,k,l)+af(2)*a(i,j,k,l+1)+af(3)*a(i,j,k,l+2)+
     1  af(4)*a(i,j,k,l+3)
   45 continue
C *** If desired, the following statement may be commented out so that
C *** the write statement is executed after each [Fe/H] interpolation.
      if(mode.ne.0) go to 165
      write(6,1004) ebv,fe,(fil(k),k=1,ndx)
C *** When MODE=0 or MODE=-1, control returns to the calling program
C *** once the [Fe/H] interpolations have been completed; i.e., no
C *** interpolations are performed within the tables.
      go to 165
C *** MODE=1: interpolations are performed for the input values of
C *** log g (= gv) and log Teff (= teff).
   50 nbc=ndx
      temp=10.**teff
      if(abs(teff-t(nt)).ge.1.e-5) go to 55
      teff=t(nt)
   55 if(abs(teff-t(1)).ge.1.e-5) go to 60 
      teff=t(1)
   60 if(teff.le.t(nt).and.teff.ge.t(1)) go to 70
   65 write(6,1005) temp
      stop ' range of the tables exceeded'
   70 if(gv.lt.-0.5.or.gv.gt.5.5) go to 155
C *** Note that ngmin=1 only if the input value of gv is < 0.0, and
C *** that ngmin=2 only if gv < 0.5 but >= 0.0
      ngmin=3
      if(gv.ge.0.5) go to 75
      ngmin=2
      if(gv.ge.0.0) go to 75
      ngmin=1
C *** This loop determines the index ngm where ng(ngm) is the lowest
C *** value of log g for which data are available for the input Teff.
   75 do 80 i=ngmin,8
      ngm=i
      if(teff.le.ttm(i)) go to 85
   80 continue
   85 ngl=ngm+2
C *** Note that ngmax=13 only if the input value of gv is > 5.0.
      ngmax=13
      if(gv.gt.5.0) go to 90
      ngmax=12
   90 ngm=ngmax-1
C *** This loop determines the index mg where g(mg) is the lowest 
C *** of the grid values of log g that are used in the interpolations.
      do 95 i=ngl,ngm
      mg=i-2
      if(gv.le.g(i)) go to 105
   95 continue
      mg=ngmax-3
  105 mmx=maxt(mg)
C *** This section increases the value of mg if the maximum value of
C *** Teff at the initial value of mg is less than the input Teff.
      if(teff.le.t(mmx)) go to 110
      mg=mg+1
      mmx=maxt(mg)
C *** This section tests whether the maximum value of Teff at the
C *** revised value of mg is >= the input Teff value.  If this condition
C *** is not satisfied, the execution terminates.
      if(teff.gt.t(mmx)) go to 65
  110 r(1)=g(mg)
      r(2)=g(mg+1)
      r(3)=g(mg+2)
      r(4)=g(mg+3)
C *** The Lagrangian coefficients for the log g interpolation are computed.
      call lgr4(r,ag,gv)
      ntm=nt-1
C *** This loop determines the value of mt where t(mt) is the lowest of
C *** the grid values of log Teff that are used in the interpolations.
      do 115 i=3,ntm
      mt=i-2
      if(teff.le.t(i)) go to 120
  115 continue
      mt=nt-3
C *** This section checks that the input T is <= the maximum value for
C *** which tabulated data are available at the lowest value of mg 
C *** in the interpolations.
  120 limt=maxt(mg)
      if(gv.le.5.0) go to 125
      limt=19
      if(teff.gt.t(limt)) go to 65
  125 if((mt+3).le.limt) go to 130
      mt=limt-3
  130 r(1)=t(mt)
      r(2)=t(mt+1)
      r(3)=t(mt+2)
      r(4)=t(mt+3)
C *** The Lagrangian coefficients for the log Teff interpolations are computed.
      call lgr4(r,at,teff)
      l=mg-1
C *** The desired bolometric corrections are calculated.
      do 140 i=1,4
      l=l+1
      do 135 k=1,ndx
      c(i,k)=at(1)*b(mt,l,k)+at(2)*b(mt+1,l,k)+at(3)*b(mt+2,l,k)+
     1   at(4)*b(mt+3,l,k)
  135 continue
  140 continue
      do 145 i=1,ndx
      bc(i)=ag(1)*c(1,i)+ag(2)*c(2,i)+ag(3)*c(3,i)+ag(4)*c(4,i)
  145 continue
      go to 165
  150 stop ' input file bc_m04.data does not exist'
  155 do 160 i=1,ndx
      bc(i)=9.9999
  160 continue
      write(6,1009) gv
  165 return
      end
      subroutine lgran3(x,a,xx)
C *** a 3-point Lagrange interpolation subroutine
C *** x(i) = abscissae values
C *** a(i) = Lagrange coefficients
C *** xx = value of x-coordinate at which an interpolation is to be made
C ***      in the calling program via yy=a(1)*y(1)+a(2)*y(2)+a(3)*y(3)
      real*4 x(3),a(3)
      r1=(x(1)-x(2))*(x(1)-x(3))
      r2=(x(2)-x(1))*(x(2)-x(3))
      r3=(x(3)-x(1))*(x(3)-x(2))
      a(1)=((xx-x(2))*(xx-x(3)))/r1
      a(2)=((xx-x(1))*(xx-x(3)))/r2
      a(3)=((xx-x(1))*(xx-x(2)))/r3
      return
      end
      subroutine lgr4(x,a,xx)
C *** a 4-point Lagrange interpolation subroutine (similar to lgran3)
      real*4 x(4),a(4)
      r1=(x(1)-x(2))*(x(1)-x(3))*(x(1)-x(4))
      r2=(x(2)-x(1))*(x(2)-x(3))*(x(2)-x(4))
      r3=(x(3)-x(1))*(x(3)-x(2))*(x(3)-x(4))
      r4=(x(4)-x(1))*(x(4)-x(2))*(x(4)-x(3))
      a(1)=((xx-x(2))*(xx-x(3))*(xx-x(4)))/r1
      a(2)=((xx-x(1))*(xx-x(3))*(xx-x(4)))/r2
      a(3)=((xx-x(1))*(xx-x(2))*(xx-x(4)))/r3
      a(4)=((xx-x(1))*(xx-x(2))*(xx-x(3)))/r4
      return
      end