c crappy program to plot frequency graph for toy star
c
      program freq 
      implicit none
      integer npts,maxstep
      parameter (npts=20)
      parameter (maxstep = 100000)
      integer i,j,nstep
      real fn(0:npts),omega(0:npts)
      real omegasph(0:npts)
      real time(maxstep), ekin(maxstep),period
      character cn*1,cn10*1,filename*20,rootname*20
      
      call getarg(1,rootname)
      if (rootname(1:1).EQ.' ') rootname = 'mtstarn'
      print*,' rootname = ',rootname
      
      do i=0,npts
         fn(i) = i
	 omega(i) = SQRT(0.5*(i+1.)*(i+2.))
	 period = 2.*3.1415926536/omega(i)
         cn = ACHAR(MOD(i,10)+48)
	 if (i.ge.10) then
	    cn10 = ACHAR(i/10 + 48)
	    filename = TRIM(rootname)//cn10//cn//'.ev'
	 else
	    filename = TRIM(rootname)//cn//'.ev'
	 endif
	 print*,' opening ',filename
	 open(unit=1, file=filename,status='old',err=20)
	     do j=1,maxstep
	        read(1,*, END=10) time(j),ekin(j)
	     enddo
	     print*,' reached maxstep = ',maxstep
	     nstep = maxstep-1
	     goto 11
10           print*,' end of file, steps = ',j-1
             nstep = j-1
	     
11           continue
	     call getmax(ekin(1:nstep),time(1:nstep),
     &	                 nstep,period,omegasph(i))
	 close(unit=1)
	 goto 30
20       continue
         print*,' file not found ',filename
	 omegasph(i) = 0.
30       continue	 
      enddo
      
      do i=0,npts
         IF (omega(i).EQ.0.) PRINT*,'eek',i,omega(i)
         print*,' n = ',i,' err = ',abs(omegasph(i)-omega(i))/omega(i)
      enddo
      
      call pgbegin(0,'?',1,1)
      call pgsch(2.0)
      call pgenv(0.0,REAL(npts+1),0.0,SQRT(0.5*(npts+2)*(npts+3)),0,1)
      call pglabel('n','\gw',' ')
      call pgsch(1.1)
      call pgpt(npts,fn,omega,17)
      call pgpt(npts,fn,omegasph,4)
      call pgend
      
      end program freq
c
c--subroutine to get maximums/mins of function
c      
      subroutine getmax(datain,time,max,perin,omega)
      implicit none
      integer i,max,nmax,nmin,maxmax,maxmin,ncount,ncount2
      real datain(max),time(max),timeprev,timeprev2
      real freq,freq2,avper,avper2,period,omega,omega2,pi
      real tol,perin
      
      pi = 3.1415926536
      tol = 0.2*perin
      print*,' tol = ',tol, perin
      
      timeprev = 0.
      timeprev2 = 0.
      avper = 0.
      avper2 = 0.
      nmax = 0
      nmin = 0
      maxmax = 10
      maxmin = 10
      ncount = 0
      ncount2 = 0
      
      do i=2,max-1
         if (((datain(i).GT.datain(i-1)).and.
     &       (datain(i).GT.datain(i+1)).AND.nmax.LT.maxmax)
     &    .and.(time(i)-timeprev).GT.tol) then
     
	    nmax = nmax + 1
	    IF (MOD(nmax,2).EQ.0) THEN	       
	       IF ((time(i)-timeprev).GT.tol) THEN
               period = time(i)-timeprev            
	       print*,' max: time, period = ',time(i),period
	       avper = avper + period
	       ncount = ncount + 1
	       timeprev = time(i)	
	       ENDIF       
	    ENDIF
	 
	 elseif (((datain(i).LT.datain(i-1)).and.
     &       (datain(i).LT.datain(i+1)).AND.nmin.LT.maxmin).and.
     &       (time(i)-timeprev2).GT.tol) then

	    nmin = nmin + 1
	    IF (MOD(nmin,2).EQ.1) THEN
	       IF ((time(i).GT.2.*period).AND.nmin.GT.1) THEN
               period = time(i)-timeprev2            
	       print*,' min: time, period = ',time(i),period
	       avper2 = avper2 + period
	       ncount2 = ncount2 + 1
	       ENDIF
	       timeprev2 = time(i)	       
	    ENDIF
	    
         endif
      enddo
      
      nmin = nmin - 1
      IF ((ncount2.GE.1).AND.(ncount.GE.1)) THEN
         avper = avper/(ncount)
         avper2 = avper2/(ncount2)	 
         IF ((avper.GT.0.).AND.(avper2.GT.0.)) THEN
	    freq = 1./avper
            freq2 = 1./avper2
         ENDIF
         omega = 2.*pi*freq
         omega2 = 2.*pi*freq2      
      ELSE
         omega = 0.
	 omega2 = 0.
      ENDIF
      print*,' ncount = ',ncount,ncount2
      print*,' average max p,f,w = ',avper,freq,omega,nmax
      print*,' average min p,f,w = ',avper2,freq2,omega2,nmin

!      omega = omega2
      
      return
      end subroutine getmax
      
