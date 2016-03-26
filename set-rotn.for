c     **** prog to set all rotatable coordinates ****
c
c	OH from TYR, SER, THR
c	SH from CYS
c	CONH2 from GLN and ASN
c	ring in HIS
c	water

      implicit real*8 (a-h,o-z)
      include 'readhin.cmn'
      character*80 fname
      character*3 s

c     **** mtyp= max number of residues with rotatable bits ****
c     **** name of each rotatable residue type
c     **** rcon= rotate all connected to this atom
c     **** rvec= around vector between atoms rcon to rvec
c     **** ratom= list of atoms to rotate
c     **** rhb=  list of hydrogen bonding heavy atoms amongst those rotating
c     **** nhvy= number of rotating hb heavy atoms
c     **** npos= number of different position angles 
c     **** pos=  value of possible rotation angles

      parameter (mtyp=19, mratm=8, mrhvy=2, mrpos=72)
      character*3 restyp(mtyp)
      real*8  rpos(mrpos,mtyp)
      integer rcon(mtyp),rvec(mtyp),ratm(mratm,mtyp),rhvy(mrhvy,mtyp),
     $	      nratm(mtyp),nrhvy(mtyp),nrpos(mtyp)

c     **** mrot= max number of rotatable residues
c     **** mhb= max number of connections froma rotatable res to others ****
c     **** lstr,lstm = residue and molecule nbers of rotatable res ****
c     **** lsti= reverse index to give the rotatable residue type ****
c     **** hb= list of rotatable residues near another ****
c     **** nhb= size of list of rotatable residues near another ****
      parameter (mrot=1000, mhb=10)
      integer lstr(mrot),lstm(mrot),hb(mhb,mrot),lsti(mrot),nhb(mrot)

c     **** todo= flag indicating residue yet to be optimized ****
c     **** ntodo= number of connected rotatable residues ****
c     **** ltodo= list of connected rotatable residues ****
c     **** eext= energies of rot res at all poss posns with ext field
c     **** ee= energies for all poss posns of all pairs of rot res 
c     **** real*4 used so as to maximize cache hits ****
c     **** nowpos= index to actual angle that each var is set to ****
      parameter (mtodo=14)
      integer ltodo(mtodo),nowpos(mtodo)
      logical todo(mrot),isin
      real*4 eext(mrpos,mtodo),ee(mrpos,mtodo,mrpos,mtodo)
      real*8 npermt

c     **** externals list for energy ****
      parameter (mextrn=100)
      integer extres(mextrn),extmol(mextrn)

c     **** indices for opt config ****
      integer ival(mtodo),nval(mtodo),ivalmn(mtodo)

c     **** water H orientations ****
      common /h2o/ water(6,72)

c     **** heavy-heavy Hbond cutoff distance ****
      data rrcut/4.0/

c     ********************************************************************

      rrcut2= rrcut**2
      open (4,file='set-rotn.out',status='unknown')

c     ******** read water orientations *******
      open (1,file='makeh2o.out',status='old')
      read (1,*) ((water(i,j),i=1,6),j=1,72)
      close (1)

c     ********* open .hin file *************

10    write (6,*) 'Enter root name of input .hin file'
C      read (5,'(a)',err=500,end=500) fname
      fname= '1jboic'
      if (fname(1:1).eq.'-') goto 500
      nf= 1
      do while (fname(nf:nf).ne.' ')
	  nf= nf + 1
      end do
      fname(nf:nf+3)= '.hin'

      if (.not.readhin (fname,0)) goto 500

c     **** read rotatable residue info ****

      open (1,file='set-rotn.dat',status='unknown')
      read (1,*) ntyp
      if (ntyp.gt.mtyp) stop 'MTYP'
      do i= 1,ntyp
	read (1,*) restyp(i),rvec(i),rcon(i),
     $		   nratm(i),(ratm(j,i),j=1,nratm(i)),
     $		   nrhvy(i),(rhvy(j,i),j=1,nrhvy(i)),
     $		   nrpos(i),(rpos(j,i),j=1,nrpos(i))
	if (nratm(i).gt.mratm) stop 'MRATM'
	if (nrhvy(i).gt.mrhvy) stop 'MRHVY'
	if (nrpos(i).gt.mrpos) stop 'MRPOS'
      end do
      nrpos(1)= 72
      do i= 1,nrpos(1)
	rpos(i,1)= i
      end do
      close (1)

c     **** make a list of all rotatable (parts of) residues ****

      nrot= 0
      do imol= 1,nmol
	do ires= 1,nres(imol)
	  s= resname(ires,imol)
	  do k= 1,ntyp
	    if (restyp(k).eq.s) then
	      nrot= nrot + 1
C	      write (6,'(a,i6,2a,3i6)') restyp(k),k,' ',s,nrot,imol,ires
	      if (nrot.gt.mrot) stop 'MROT'
	      lstr(nrot)= ires
	      lstm(nrot)= imol
	      lsti(nrot)= k
	    end if
	  end do
	end do
      end do
      write (6,*) 'nber of (non distinct) rotatable residues=',nrot

c     **** find h-bonded rotatable residues around each rot res ****

      do irot= 1,nrot
	ires= lstr(irot)
	imol= lstm(irot)
	ityp= lsti(irot)
	i0= iares(ires,imol) - 1
	nhb(irot)= 0
	todo(irot)= .true.

        do jrot= 1,nrot
	 if (irot.ne.jrot) then
	  jres= lstr(jrot)
	  jmol= lstm(jrot)
	  jtyp= lsti(jrot)
	  j0= iares(jres,jmol) - 1

	  do ihv= 1,nrhvy(ityp)
	    i= rhvy(ihv,ityp) + i0
	    do jhv= 1,nrhvy(jtyp)
	      j= rhvy(jhv,jtyp) + j0
	      r2= hbl2 (i,j)
	      if (r2.lt.rrcut2) then
		nhb(irot)= nhb(irot) + 1
		if (nhb(irot).gt.mhb) stop 'MHB'
		hb(nhb(irot),irot)= jrot
		write (4,*) 'hb ',resname(ires,imol),atname(i),
     $			sqrt(r2),' ', resname(jres,jmol),atname(j)
	      end if
	    end do
	  end do

	 end if
	end do

	if (ityp.eq.1) write (4,860) irot,resname(ires,imol),nhb(irot)
860	format (i4,1x,a,' # hbonds=',i4)

      end do

c     ************ scan each rotatable reside, opt all those con to it ********'

      ncalc= 0
      ntodomax= 0

      do irot= 1,nrot
       if (todo(irot)) then

	todo(irot)= .false.
	ntodo= 1
	ltodo(ntodo)= irot

c	**** add connected residues to list ****
	itodo= 1
	do while (itodo.le.ntodo)
	  nhbnew= nhb(ltodo(itodo))
	  do j= 1,nhbnew
c	    *** test to see if this residue already in todo list ****
	    jtest= hb(j,ltodo(itodo))
	    isin= .true.
	    do k= 1,ntodo
	      if (ltodo(k).eq.jtest) isin= .false.
	    end do
	    if (isin) then
c	      **** add to todo list ****
	      ntodo= ntodo + 1
	      if (ntodo.gt.mtodo) stop 'MTODO'
	      ltodo(ntodo)= jtest
	      if (.not.todo(jtest)) stop 'redoing a residue?'
	      todo(jtest)= .false.
	    end if
	  end do
	  itodo= itodo + 1
	end do

c	******** set all res within list ********

	ncalc= ncalc + 1
	nh2o= 0
	nperm= 1
	do i= 1,ntodo
	  ii= ltodo(i)
	  if (lsti(ii).eq.1) then
c	    **** its a water ****
	    nh2o= nh2o + 1
	  else
	    nperm= nperm * nrpos(lsti(ii))
C	    write (6,*) i,ii,lsti(ii),restyp(lsti(ii)),
C     $		nrpos(lsti(ii)),nperm
	  end if
	end do

	npermt= nperm * 72.d0**nh2o
	write (4,820) ncalc,ntodo,nh2o,nperm,npermt
	write(4,821) (ltodo(i),lstr(ltodo(i)),lstm(ltodo(i))
     $	    ,resname(lstr(ltodo(i)),lstm(ltodo(i)))(1:3),i=1,ntodo)
820	format (' calc=',i4,' nber res=',i3,' nber H2O=',i2,' perm=',i9,
     $		' tot=',f15.0)
821	format (5(3i4,1x,a3))

	if (nperm.gt.0) then
C	if (nperm.gt.1000 .or. nh2o.ge.3) then
C	if (ntodo.gt.ntodomax) then
     	  write (6,820) ncalc,ntodo,nh2o,nperm,npermt
	  write(6,821) (ltodo(i),lstr(ltodo(i)),lstm(ltodo(i))
     $	    ,resname(lstr(ltodo(i)),lstm(ltodo(i)))(1:3),i=1,ntodo)
	  ntodomax= ntodo
	  write (4,*) 'new maximum connected number=',ntodo
	  write (6,*) 'new maximum connected number=',ntodo
	  call hdeselect (1,natom)
	  do i= 1,ntodo
	    call hselres (lstr(ltodo(i)),lstm(ltodo(i)),1)
	  end do
	end if

C	if (ncalc.le.158) goto 4444

c	**** form externals list ****

	write (6,*) 'forming externals list'
	next= 0
	do jmol= 1,nmol
	  do jres= 1,nres(jmol)

	    do itodo= 1,ntodo
	      ires= lstr(ltodo(itodo))
	      imol= lstm(ltodo(itodo))
	      if (jmol.eq.imol .and. jres.eq.ires) goto 333
	    end do

	    do itodo= 1,ntodo
	      ires= lstr(ltodo(itodo))
	      imol= lstm(ltodo(itodo))

	      do i= iares(ires,imol),jares(ires,imol)
	        do j= iares(jres,jmol),jares(jres,jmol)
		  r2= hbl2 (i,j)
		  if (r2.lt.16.) then
c		    **** include in list ***
		    next= next + 1
		    if (next.gt.mextrn) stop 'MEXTRN too small'
		    extres(next)= jres
		    extmol(next)= jmol
		    call hselres (jres,jmol,1)
		    goto 333
		  end if
	        end do
	      end do

	    end do
333	    continue
	  end do
	end do

	write (6,*) 'nber of ext residues=',next

c	**** get all pairwise energies at all possible positions ****

	do itodo= 1,ntodo
	  ires= lstr(ltodo(itodo))
	  imol= lstm(ltodo(itodo))
	  ityp= lsti(ltodo(itodo))
	  nval(itodo)= nrpos(ityp)
	  nowpos(itodo)= 1
	  do ipos= 1,nrpos(ityp)
	    call setpos (ires,imol,rcon(ityp),rvec(ityp),rpos(1,ityp),
     $		ipos,nowpos(itodo),ityp,nratm(ityp),ratm(1,ityp))
	    esum= 0.d0
	    do j= 1,next
	      esum= esum + resres (ires,imol,extres(j),extmol(j))
	    end do
	    eext(ipos,itodo)= esum

	    do jtodo= 1,itodo - 1
	      jres= lstr(ltodo(jtodo))
	      jmol= lstm(ltodo(jtodo))
	      jtyp= lsti(ltodo(jtodo))
	      do jpos= 1,nrpos(jtyp)
	        call setpos (jres,jmol,rcon(jtyp),rvec(jtyp),rpos(1,jtyp),
     $		  jpos,nowpos(jtodo),jtyp,nratm(jtyp),ratm(1,jtyp))
		ee(ipos,itodo,jpos,jtodo)= resres (ires,imol,jres,jmol)
		ee(jpos,jtodo,ipos,itodo)= ee(ipos,itodo,jpos,jtodo)
	      end do
	    end do

	  end do
	end do

c	**** find lowest energy str ****

	if (npermt.lt.1.e7) then
	  call opt (ntodo,mtodo,mrpos,eext,ee,ivalmn,ival,nval)
	else
	  call montec (ntodo,mtodo,mrpos,eext,ee,ivalmn,ival,nval)
        end if

        write (6,875) ncalc,(ivalmn(i),i=1,ntodo)
        write (4,875) ncalc,(ivalmn(i),i=1,ntodo)
875	format (i4,' Emin=',20i3)

c	**** set coords to this str ****

	do itodo= 1,ntodo
	  ires= lstr(ltodo(itodo))
	  imol= lstm(ltodo(itodo))
	  ityp= lsti(ltodo(itodo))
	  ipos= ivalmn(itodo)
	  call setpos (ires,imol,rcon(ityp),rvec(ityp),rpos(1,ityp),
     $		ipos,nowpos(itodo),ityp,nratm(ityp),ratm(1,ityp))
	end do

	write (fname,'(4hcalc,i3.3,4h.hin)') ncalc
	if (nh2o.gt.2) call writehin (fname)

4444	continue
C	if (ncalc.eq.85) stop 'ncalc 85'
       end if
      end do

500   continue

      end

c     ****************************************************************

      subroutine opt (ntodo,mtodo,mrpos,eext,ee,ivalmn,ival,nval)

      implicit real*8 (a-h,o-z)
      include 'readhin.cmn'

      real*4 eext(mrpos,mtodo),ee(mrpos,mtodo,mrpos,mtodo)
      integer ival(ntodo),nval(ntodo),ivalmn(ntodo)
      integer ivalout (6)
      data ivalout /1,2,12,1,0,1/

c     **** initial config ****

      etot= 0.d0
      do itodo= 1,ntodo
	ival(itodo)= 1
	ivalmn(itodo)= 1
	etot= etot + eext(1,itodo)
	do jtodo= 1,itodo-1
	  etot= etot + ee(1,itodo,1,jtodo)
	end do
      end do

c     **** loop over variations ****

      emin= 1.e30
100   continue

	etot= 0.d0
        do itodo= 1,ntodo
	  etot= etot + eext(ival(itodo),itodo)
	  do jtodo= 1,itodo-1
	    etot= etot + ee(ival(jtodo),jtodo,ival(itodo),itodo)
	  end do
	end do

	do i= 1,ntodo
	  if (ivalout(i).gt.0 .and. ival(i).ne.ivalout(i)) goto 200
	end do
C	write (6,810) etot,(ival(i),i=1,ntodo)
810	format (' E=   ',f8.4,20i3)
200	continue

	if (etot.lt.emin) then
	  emin= etot
	  do itodo= 1,ntodo
	    ivalmn(itodo)= ival(itodo)
	  end do
	end if

	ival(1)= ival(1) + 1
	i= 1
	do while (i.le.ntodo .and. ival(i).gt.nval(i))
	 ival(i)= 1
	 i= i + 1
	 if (i.le.ntodo) ival(i)= ival(i) + 1
	end do

	if (i.le.ntodo) goto 100


      return
      end

c     ****************************************************************

      function resres (ires,imol,jres,jmol)

c     **** coulomb energy of two residues ****

      implicit real*8 (a-h,o-z)
      include 'readhin.cmn'

      ener= 0.d0
      rmin= 1.e30
      do i= iares(ires,imol),jares(ires,imol)
	do j= iares(jres,jmol),jares(jres,jmol)
	  r= hbl2(i,j)
	  rmin= min (rmin,r)
C	  r= sqrt (hbl2(i,j))
	  ener1= q(j) * q(i) / r
	  ener= ener + ener1
C	  if (i.eq.4381 .and. jres.eq.3175) write (6,'(2i5,2a,5f8.4)')
C     $		i,j,atname(i),atname(j),q(i),q(j),sqrt(r),ener1,ener
	end do
      end do

C      if (ires.eq.297 .and. jres.eq.216)
C     $     write (6,800) resname(ires,imol),resname(jres,jmol),ener,
C     $	   sqrt(rmin)
800   format (2(1x,a),' INTER=',f6.3,2f8.4)

      resres= ener
      return

      end

c     ****************************************************************

      subroutine setpos (ires,imol,rcon,rvec,rpos,ipos,nowpos,
     $		ityp,nratm,ratm)

c     **** rotates the nratm's in list ratm of ires,imol about the atom 
c     **** con through angle rvec in the direction of the vector from
c     **** rcon to rvec
c

      implicit real*8 (a-h,o-z)
      include 'readhin.cmn'
      integer ratm(nratm),rcon,rvec
      real*8 r(3,3),r2(3,3),r3(3,3),xx(3),rpos(*)
      character*80 fname
      common /h2o/ water(6,72)
      data nc/0/

      i0= iares(ires,imol) - 1

      if (ityp.eq.1) then
c	**** set std water config ****
	x(i0+2)= water(1,ipos) + x(i0+1)
	y(i0+2)= water(2,ipos) + y(i0+1)
	z(i0+2)= water(3,ipos) + z(i0+1)
	x(i0+3)= water(4,ipos) + x(i0+1)
	y(i0+3)= water(5,ipos) + y(i0+1)
	z(i0+3)= water(6,ipos) + z(i0+1)

      else
c	**** other residue ****

c	**** vector to rotate about ****

	r(3,1)= x(i0+rcon) - x(i0+rvec)
	r(3,2)= y(i0+rcon) - y(i0+rvec)
	r(3,3)= z(i0+rcon) - z(i0+rvec)

c     **** transformation from vector in r(3,i) to the z axis ****
c     **** determine polar angles theta and phi ****

      st= sqrt ( r(3,1)**2 + r(3,2)**2 + r(3,3)**2)
      r(3,1)= r(3,1) / st
      r(3,2)= r(3,2) / st
      r(3,3)= r(3,3) / st
      ct= r(3,3)
      st= sqrt (1.D0-ct**2)
      if (abs(st).lt.1.e-5) then
	cp= 1.D0
	sp= 0.D0
      else
	cp= r(3,1)/st
	sp= r(3,2)/st
      end if

c     **** rest of rotation matrix of vector to Z ****
      r(1,1)= - sp
      r(1,2)=   cp
      r(1,3)= 0.D0
      r(2,1)= - ct*cp
      r(2,2)= - ct*sp
      r(2,3)=   st
C      write (6,820) 'r',((r(i,j),j=1,3),i=1,3)

c     **** rotation matrix about vector ****
      delang= rpos(ipos) - rpos(nowpos)
      nowpos= ipos
      st= delang/180.d0*3.14159265d0
      ct= cos (st)
      st= sin (st)
      r2(1,1)= ct
      r2(1,2)= -st
      r2(1,3)= 0.D0
      r2(2,1)= st
      r2(2,2)= ct
      r2(2,3)= 0.D0
      r2(3,1)= 0.D0
      r2(3,2)= 0.D0
      r2(3,3)= 1.D0
C      write (6,820) 'r2',((r2(i,j),j=1,3),i=1,3)

c     **** final rotation matrix ****
      call mmult   (r3,r2,r,3,3)
      call mmultta (r2,r,r3,3,3)
C      write (6,820) 'r2',((r2(i,j),j=1,3),i=1,3)
C      write (6,*) 'nratm=',nratm

c     **** rotate connected atoms ****
      do i1= 1,nratm
	i= i0 + ratm(i1)
	xx(1)= x(i) - x(i0+rcon)
	xx(2)= y(i) - y(i0+rcon)
	xx(3)= z(i) - z(i0+rcon)
C	write (6,'(3f10.4)') xx
	call crotn (xx,1,r2)
C	write (6,'(3f10.4)') xx
	x(i)= xx(1) + x(i0+rcon)
	y(i)= xx(2) + y(i0+rcon)
	z(i)= xx(3) + z(i0+rcon)
      end do

820   format (1x,a2 / (1x,3f10.6) )

      end if

C      nc= nc + 1
C      write (fname,'(1hf,i3.3,4h.dat)') nc
C      open (9,file=fname,status='unknown')
C      write (9,'(i3,3f12.6)') (nat(i),x(i),y(i),z(i),
C     $	   i=iares(ires,imol),jares(ires,imol))
C      write (9,*)
C      write (9,*) 'coords for ',resname(ires,imol),
C     $		' new pos=',rpos(ipos)
C      close (9)

      return
      end

c     ****************************************************************

      subroutine montec (ntodo,mtodo,mrpos,eext,ee,ivalmn,ival,nval)

      implicit real*8 (a-h,o-z)
      include 'readhin.cmn'

      real*4 eext(mrpos,mtodo),ee(mrpos,mtodo,mrpos,mtodo)
      real*4 rand
      integer ival(ntodo),nval(ntodo),ivalmn(ntodo)
      logical accept
      data temp1/0.08/, nmove/10000000/

c     **** initial config ****

      etot= 0.d0
      do itodo= 1,ntodo
	ival(itodo)= 1
	ivalmn(itodo)= 1
	etot= etot + eext(1,itodo)
	do jtodo= 1,itodo-1
	  etot= etot + ee(1,itodo,1,jtodo)
	end do
      end do
      emin= etot

      do imove= 1,nmove

c     **** select resdiue and value ****
      itodo= ifix (rand(0)*ntodo) + 1
      if (itodo.le.0 .or. itodo.gt.ntodo) then
	write (6,*) itodo
	stop 'itodo'
      end if
      ipos=  ifix (rand(0)*nval(itodo)) + 1
      if (ipos.le.0 .or. ipos.gt.nval(itodo)) then
	write (6,*) ipos,itodo,nval(itodo)
	stop 'ipos'
      end if
      dele= eext(ipos,itodo) - eext(ival(itodo),itodo)
      do jtodo= 1,ntodo
	if (jtodo.ne.itodo) then
	  dele= dele + ee(ipos,       itodo,ival(jtodo),jtodo) 
     $		     - ee(ival(itodo),itodo,ival(jtodo),jtodo) 
	end if
      end do

      accept= dele.lt.0.d0
      temp= temp1 * float(imove) / nmove
C      write (6,*) 'dele=',dele
      if (.not.accept) accept= exp(-dele/temp).gt.rand(0)
      if (accept) then
C	write (6,*) 'accept', dele
	naccept= naccept + 1
	ival(itodo)= ipos
	etot= etot + dele
      end if

      if (etot.lt.emin) then
	emin= etot
	write (6,900) etot,imove,(ival(i),i=1,ntodo)
900	format (' new min=',f8.4,' mov=',i9,' at',20i3)
	do i= 1,ntodo
	  ivalmn(i)= ival(i)
	end do
      end if

      end do

      write (4,*) 'MC ratio=',float(naccept)/nmove,' E=',emin
      write (6,*) 'MC ratio=',float(naccept)/nmove,' E=',emin

      return
      end
