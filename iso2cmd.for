C----------------------------- Program iso2cmd-----------------------------
C
C *** This program, which must be linked with "bcutil.for" to produce an
C *** executable module, converts isochrones from the theoretical
C *** (log Teff, log L)- diagram to user-selected observational planes.
C *** To accomplish this, interpolations are made in tables of bolometric
C *** corrections (BCs) derived from the Gustafsson et al (2008) MARCS
C *** model atmospheres, which have been computed for Teffs from 2600 K
C *** to 8000 K.  As these model atmospheres have been computed on the
C *** assumption of 4 different variations of [alpha/Fe] with [Fe/H], it
C *** is possible to interpolate in these results for arbitrary values of
C *** [alpha/Fe] within the range -0.4 <= [alpha/Fe] <= +0.4: see the 2014,
C *** MNRAS, 444, 392 paper by Casagrande & VandenBerg for a complete
C *** description of the many photometric systems and filters that are
C *** considered).  Separate subroutines to interpolate in the associated 
C *** BC tables are provided for these four choices:  
C ***                subroutine           BC table
C                   getbc_std.for        bc_std.data
C                   getbc_m04.for        bc_m04.data
C                   getbc_p00.for        bc_p00.data
C                   getbc_p04.for        bc_p04.data
C *** Note that the interpolation subroutines open and close the respective
C *** data files (which are explicitly "attached" to unit 30) internally.
C *** These data files must be present in the working directory.
C        
C *** The program prompts the user for the name of the input file, which is
C *** assumed to be in the format of a "*.iso" file produced by PBISO.  The
C *** program also provides a default name for the output file, which is
C *** identical to the input name, but with the extension ".cmd" replacing
C *** ".iso".  The user may type in an alternative name for the ouput file
C *** if something other than the default name is preferred.
C *** Example of an input file:        a0zz_p4y25m15.iso
C *** Example of an output file:       a0zz_p4y25m15.cmd
C
C *** Note that the MARCS transformations (by CV14) extend to [Fe/H] = -4.0
C *** (if [alpha/Fe] = +0.4), -2.5 (if [alpha/Fe] = 0.0), and -2.0 (if
C *** [alpha/Fe] = -0.4).  Normally quadratic interpolation is employed
C *** for intermediate [alpha/Fe] values, though in the range -2.5 <
C *** [Fe/H] < -2.0, linear interpolation of the results for [alpha/Fe]
C *** = 0.0 and 0.4 is used to obtain the transformations for intermediate
C *** values of [alpha/Fe].
C 
C -------------------------------------------------------------------7--
C
      real*8 mixx(5),mass
      real*4 xa(11),bc(5),bcm4(5),bcp0(5),bcp4(5),x(3),a(3),mbol
      character infi*60,outfi*60,sol*50,zmx*47,xx*31,yy*31,ele*22,
     1  lb(6)*18,sys(5)*8,fil(5)*6
 1000 format(a18,i3,/,a18,f7.3,/,a18,3f6.2,/,a18,4f6.2,/,a18,4f6.2,/,
     1  a47,/,a31,/,a31,/,a18,f6.3,/,a50)
 1001 format(a22,1p,5d13.6)
 1002 format(//,f6.2,2x,i4)
 1003 format(/,'  Age',3x,'Npts',8x,'E(B-V) =',f6.3,' is assumed in',
     1  1x,'the Bolometric Corrections',/,f6.2,2x,i4,30x,5(1x,a8))
 1004 format(8x,'Mass',7x,'log Te',4x,'log g',3x,'M_bol',3x,5(a6,3x))
 1005 format(5x,2f9.6,f13.10)
 1006 format(i4,f13.10,f9.6,f8.4,f9.5,f8.4,4f9.4)
c1007 format(' enter (*) correction to predicted log Teff values')
c     write(6,1007)
c     read(5,*) delteff
c
c *** enter the name of the input *.iso file
c
      call iofile(10,'in','file.iso',.false.,' ',infi,.true.,
     1  'Input Isochrone File: ')
c
c *** the default name of the output file is the same as the input name
c *** but with the extension ".cmd" ... hit return to accept the default
c *** name, or else type in an alternative name for the output file 
c
      ncdot=index(infi,'.')
      nch=ncdot+3
      outfi(1:nch)=infi(1:ncdot)//'cmd'
      call iofile(11,'out',outfi(1:nch),.true.,'Default Iso-CMD File: ', 
     1  outfi,.false.,'Output  Iso-CMD File: ')
c
c *** copy the header lines from the input file to the output file and
c *** define the values of niso, fe, and al2fe
c *** note that al2fe = [alpha/Fe] is defined in the following statements
c *** (this parameter is very important because it determines whether
c *** interpolations are made in the bc_std.data, bc_p04.data, bc_p00.data,
C *** and/or bc_m04.data files
c
c     write(6,1007)
      read(10,1000) lb(1),niso,lb(2),fe,lb(3),(xa(i),i=1,3),lb(4),
     1  (xa(i),i=4,7),lb(5),(xa(i),i=8,11),zmx,xx,yy,lb(6),amlt,sol
      al2fe=xa(11)
      if(al2fe.ge.-0.4.and.al2fe.le.0.4) go to 5
      stop '[alpha/Fe] is outside the permitted range from -0.4 to +0.4'
    5 write(11,1000) lb(1),niso,lb(2),fe,lb(3),(xa(i),i=1,3),lb(4),
     1  (xa(i),i=4,7),lb(5),(xa(i),i=8,11),zmx,xx,yy,lb(6),amlt,sol
      do 10 i=1,5
      read(10,1001) ele,(mixx(j),j=1,5)
      write(11,1001) ele,(mixx(j),j=1,5)
   10 continue 
c
c *** if [Fe/H] < -2.5, mode = 0 (bc_std.data is used): the execution
c *** continues only if [alpha/Fe] = 0.4
c
      mode=0
      if(fe.ge.-2.5) go to 20
      if(al2fe.ne.0.4) stop '[alpha/Fe] should = 0.4 if [Fe/H] < -2.5'
   15 mode=0
      call getbc_std(0,fe,gv,teff,ebv,bc,fil,sys,nbc)
      go to 45
c
c *** if [alpha/Fe] = -0.4, mode = 1 (bc_m04.data is used)
c
   20 if(al2fe.gt.-0.4) go to 25
      if(fe.lt.-2.0) stop '[Fe/H >= -2.0 required if [alpha/Fe] = -0.4'
      mode=1
      call getbc_m04(0,fe,gv,teff,ebv,bc,fil,sys,nbc)
      go to 45
c
c *** if [alpha/Fe] = 0.4, mode = 3 (bc_p04.data is used)
c
   25 if(al2fe.lt.0.4) go to 30
      mode=3
      call getbc_p04(0,fe,gv,teff,ebv,bc,fil,sys,nbc)
      go to 45
c
c *** if [alpha/Fe] = 0.0, mode = 2 (bc_p00.data is used) - unless
c *** [Fe/H] > 0.0, in which case, mode=0 (bc_std.data, which includes
c *** BCs for log g = -0.5, is used instead)
c
   30 if(al2fe.ne.0.0) go to 35
      if(fe.gt.0.0) go to 15
      mode=2
      if(fe.lt.-2.5) stop '[Fe/H] >= -2.5 required if [alpha/Fe] = 0.0'
      call getbc_p00(0,fe,gv,teff,ebv,bc,fil,sys,nbc)
      go to 45
c
c *** mode = 4 if 3-point interpolations w.r.t. [alpha/Fe] are employed
c *** mode = 5 if linear interpolation w.r.t. [alpha/Fe] is employed
c
   35 mode=4
      if(fe.ge.-2.0) go to 40
      mode=5
      go to 41
   40 call getbc_m04(0,fe,gv,teff,ebv,bc,fil,sys,nbc)
c
c *** calculate the Lagrange coefficients for the [alpha/Fe] interpolation
c
      x(1)=-0.4
      x(2)=0.0
      x(3)=0.4
      call lgran3(x,a,al2fe)
   41 call getbc_p00(0,fe,gv,teff,ebv,bc,fil,sys,nbc)
      call getbc_p04(0,fe,gv,teff,ebv,bc,fil,sys,nbc)
c
c *** interpolate in BC tables for the selected filters along the isochrones
c
   45 do 85 k=1,niso
      read(10,1002) age,npt
      write(11,1003) ebv,age,npt,(sys(i),i=1,nbc)
      write(11,1004) (fil(i),i=1,nbc)
      do 80 i=1,npt
      read(10,1005) xlum,teff,mass
c     teff=teff+delteff
      gv=sngl(dlog10(mass))-10.60917+4.*teff-xlum
      if(mode.gt.0) go to 50
      call getbc_std(1,fe,gv,teff,ebv,bc,fil,sys,nbc)
      go to 75
   50 if(mode.gt.1) go to 55
      call getbc_m04(1,fe,gv,teff,ebv,bc,fil,sys,nbc)
      go to 75
   55 if(mode.gt.2) go to 60
      call getbc_p00(1,fe,gv,teff,ebv,bc,fil,sys,nbc)
      go to 75
   60 if(mode.gt.3) go to 65
      call getbc_p04(1,fe,gv,teff,ebv,bc,fil,sys,nbc)
      go to 75
   65 if(mode.gt.4) go to 71
      call getbc_m04(1,fe,gv,teff,ebv,bcm4,fil,sys,nbc)
      call getbc_p00(1,fe,gv,teff,ebv,bcp0,fil,sys,nbc)
      call getbc_p04(1,fe,gv,teff,ebv,bcp4,fil,sys,nbc)
      do 70 m=1,nbc
      bc(m)=a(1)*bcm4(m)+a(2)*bcp0(m)+a(3)*bcp4(m)
   70 continue
      go to 75
   71 call getbc_p00(1,fe,gv,teff,ebv,bcp0,fil,sys,nbc)
      call getbc_p04(1,fe,gv,teff,ebv,bcp4,fil,sys,nbc)
      do 72 m=1,nbc
      bc(m)=bcp0(m)+2.5*al2fe*(bcp4(m)-bcp0(m))
      if(i.ne.npt) go to 72
   72 continue
c
c *** write the model parametes and bolometric corrections for the
c *** the selected filter bandpasses.
c
   75 mbol=4.75-2.5*xlum
      write(11,1006) i,mass,teff,gv,mbol,(bc(j),j=1,nbc)
   80 continue
   85 continue
      close(unit=10,status='keep')
      close(unit=11,status='keep')
      stop
      end
       SUBROUTINE IOFILE(NUNIT,FSTAT,DNME,DASK,DPROM,FNME,ASK,PROM)
        INTEGER NUNIT
        CHARACTER*(*) FSTAT,FNME,PROM,DPROM,DNME
        CHARACTER*40 TNME*40, ANS*1, FMT*20
        LOGICAL ASK, MORE, DASK, TASK
C
        IF (DASK) THEN
          FMT(1:15) = '(T  ,A,A,A,A,$)'
          LP = LEN(DPROM)
          LD = LEN(DNME)
          WRITE(FMT(3:4),'(i2.2)') 45 - LP - LD
          WRITE(*,FMT(1:15)) DPROM(1:LP),'[',DNME(1:LD),']: '
          READ(*,'(A)') TNME
          IF (TNME.EQ.' ') THEN
            FNME = DNME
          ELSE
            FNME = TNME
          END IF
        END IF
C       
        MORE = .TRUE.
        TASK = ASK
        DO WHILE (MORE)
          IF (TASK) THEN
           FMT(1:9) = '(T  ,A,$)'
           LP = LEN(PROM)
           WRITE(FMT(3:4),'(I2.2)') 49 - LP
           WRITE(*,FMT(1:9)) PROM
            READ(*,'(A)') FNME
          END IF
C         
          NCH=1
          DO WHILE (FNME(NCH:NCH).NE.' ')
            NCH = NCH + 1
          END DO
          NCH = NCH - 1
C
          IF (FSTAT.EQ.'IN'.or.fstat.eq.'in') THEN
            OPEN (NUNIT,iostat=ierr,FILE=FNME,STATUS='OLD')
            if (ierr.eq.0) then
              more=.false.
            else
              WRITE(*,'(t45,a,a,a)') '*** ',FNME(1:NCH),
     1                ' DOES NOT EXIST ***'
              if (.not.ask) task = .true.
            end if
C
          ELSE IF(FSTAT.EQ.'OUT'.or.fstat.eq.'out') THEN
            OPEN (NUNIT,iostat=ierr,FILE=FNME,STATUS='NEW')
            if (ierr.eq.0) then
              more=.false.
            else
              fmt = '(t  ,a,a,a)'
              nt = 40-(nch+23)/2
              write(fmt(3:4),'(i2)') nt
              WRITE(*,fmt(1:11)) 
     1                        '*** ',FNME(1:NCH),' ALREADY EXISTS ***'
              write(*,'(T23,''Do you want to OVERWRITE? '',$)')
              read(*,'(a)') ans
              if (ans.eq.'y'.or.ans.eq.'Y') then
                open(nunit,file=fnme,status='old')
                close(nunit,status='delete')
                open(nunit,file=fnme,status='new')
                more = .false.
              else         
                WRITE(*,'(T34,''New file name: '',$)')
                READ(*,'(A)') FNME
              end if
            end if
          END IF
       end do
       return
       END