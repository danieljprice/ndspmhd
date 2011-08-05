        PROGRAM injsse

c-------------------------------------------------
c     Torus setup subroutine.
c     See, eg, Stone, Pringle & Begelman, 1999, MNRAS, 310, 1002
c-------------------------------------------------
      implicit none

      integer jrings,klayers,njk,Ntot_th,npl,i,j,n
      integer nrings,nlayers
      double precision sumA,rp,thetap,zp,vthet,vt
      double precision Mtot_th,msph,r_in,r_out,zzmax,pi,amu
      double precision deltar,deltaz,dd,npoly
      double precision rj,zk,Mjk,fac_njk,AAb
      real sob(3)
      
c     Maximum allowed particles
      parameter(npl=1500)

c     Particle data - normally delcared in include file
      double precision x(npl),y(npl),z(npl),rad(npl)
      double precision vx(npl),vy(npl),vz(npl)
      double precision dens(npl),pmass(npl),spsound(npl),h(npl)

c     Random numbers
      double precision ral(npl),thetaal(npl),zal(npl)
        
c     Mass of star
      amu=1.0d0
c     pi
      pi = 3.14

c     This is 1/(gamma-1)
c     Adiabatic example
        npoly = 3.0d0/2.0d0
c     Disc mass
        Mtot_th = 0.1d0
c     Number of particles - honoured approximately
        Ntot_th = npl
c     SPH particle mass
        msph = Mtot_th/real(Ntot_th)
c     Number of annuli and z-layers
        nrings = 50
        nlayers = 20
c     Maximum and minimum radius computed for
        r_in = 0.1d0
        r_out = 1.5d0
        deltar = (r_out-r_in)/real(nrings)
c     Maximum z computed for (symetrical)
        zzmax = 0.5d0
        deltaz = 2.0d0*zzmax/real(nlayers)
c     Distortion factor 
        dd = 1.1d0

c     Random or quasirandom numbers set in advance here
        call randinit(17)
        call sobseq(-1,sob)

c     Random numbers
c        do i=1,npl
c           ral(i)=rand()
c        enddo
c        call randinit(200)
c        do i=1,npl
c           thetaal(i)=rand()
c        enddo
c        call randinit(3406)
c        do i=1,npl
c           zal(i)=rand()
c        enddo

c     Or quasirandom ones
        do i=1,npl
           call sobseq(3,sob)
           ral(i)=sob(1)
           thetaal(i)=sob(2)
           zal(i)=sob(3)
        enddo


c     Sum total mass
        sumA = 0
        do jrings = 1,nrings
           rj = r_in + (real(jrings) - 1.0d0)*deltar
           do klayers = 1,nlayers
              zk = (real(klayers) - 1.0d0)*deltaz - zzmax
              fac_njk = (rj**2.0d0+zk**2.0d0)**(-0.5d0)-(0.5d0*
     &((rj)**(-2.0d0)))-1.0d0/(2.0d0*dd)
              sumA = sumA + max(rj*fac_njk**npoly,0.0d0)
           enddo
        enddo
        AAb = (2*Pi*deltar*deltaz*sumA/Mtot_th)
        print*,'Using P = ',((AAb)**(1./npoly))/(npoly+1),'*rho^gamma'
       

c     Compute density in cells
        n = 0
        do jrings=1,nrings
           rj = r_in + (real(jrings) - 1.0d0)*deltar
           do klayers = 1,nlayers
              zk = (real(klayers) - 1.0d0)*deltaz - zzmax
              fac_njk = (rj**2+zk**2)**(-0.5)-(0.5d0*
     &((rj)**(-2.0d0)))-1.0d0/(2.0d0*dd)
              Mjk = 2.0d0*Pi*deltar*deltaz*rj*
     &             (fac_njk**npoly)/AAb
              njk = int(Mjk/msph)
              njk = max(njk,0)

c     We have njk particles in this cell
              do i=n+1,n+njk
                 rp=r_in + (real(jrings)+ral(i)-1.5d0)*deltar
                 thetap=2.0d0*pi*thetaal(i)
                 zp=(real(klayers)+zal(i)-1.5d0)*deltaz - zzmax

c     Particle positions
                 x(i)=rp*cos(thetap)
                 y(i)=rp*sin(thetap)
                 z(i)=zp
                 rad(i)=rp

c     Keplerian velocities
                 vx(i)=-sqrt(amu/rp)*sin(thetap)
                 vy(i)=sqrt(amu/rp)*cos(thetap)
                 vz(i)=0.d0

c     Mass, densitiy and sound speed
                 pmass(i)=msph
                 dens(i)=(fac_njk**npoly)/AAb 
                 spsound(i)=sqrt(fac_njk/npoly)

              enddo
              n = n+njk
           enddo
        enddo


        open(66,file='start.out',status='unknown')
        do i=1,n
           write(66,*) x(i),y(i),z(i),rad(i)
     $          ,vx(i),vy(i),vz(i),pmass(i),dens(i),spsound(i)
        enddo
        close(66)
        
        return 
        END

      SUBROUTINE randinit(seed)
C----------------------------------------------------------------------C
C                                                                      C
C  Lagged Fibonacci random number generator RANMAR.                    C
C  Must be initialized with randinit() before use.                     C
C                                                                      C
C  See F. James, Comp. Phys. Comm. 60, 329 (1990), or                  C
C  G. Marsaglia et al., Stat. Prob. Lett. 9, 35 (1990).                C
C                                                                      C
C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C                                                                      C
C This is the initialization routine RMARIN for the random number      C
C     generator RANMAR                                                 C
C                                                                      C
C NOTE: The seed variables can have values between:  0 <= IJ <= 31328  C
C                                                    0 <= KL <= 30081  C
C----------------------------------------------------------------------C
      IMPLICIT NONE
      INTEGER seed
      INTEGER ij,kl, i,j,k,l, ii,jj, m
      DOUBLE PRECISION s,t
      INTEGER Maxseed
      PARAMETER (Maxseed = 900000000)
c      DOUBLE PRECISION u(97), c, cd, cm
      REAL u(97), c, cd, cm
      INTEGER i97, j97, ivec
      COMMON /raset1/ u, c, cd, cm, i97, j97, ivec 
c$OMP THREADPRIVATE (/raset1/)
c      seed = mod(seed,Maxseed)
      ij = seed / 30082
      kl = seed - (30082 * ij)
      i = mod(ij/177, 177) + 2
      j = mod(ij    , 177) + 2
      k = mod(kl/169, 178) + 1
      l = mod(kl,     169)
      DO 2 ii = 1, 97
        s = 0.0
        t = 0.5
        DO 3 jj = 1, 24
          m = mod(mod(i*j, 179)*k, 179)
          i = j
          j = k
          k = m
          l = mod(53*l+1, 169)
          IF (mod(l*m, 64) .ge. 32) then
            s = s + t
          ENDIF
          t = 0.5 * t
    3   CONTINUE
        u(ii) = s
    2 CONTINUE
      c = 362436.0 / 16777216.0
      cd = 7654321.0 / 16777216.0
      cm = 16777213.0 /16777216.0
      i97 = 97
      j97 = 33
      RETURN
      END


      FUNCTION rand()
C----------------------------------------------------------------------C
C                                                                      C
C  Lagged Fibonacci random number generator RANMAR().                  C
C                                                                      C
C----------------------------------------------------------------------C
      IMPLICIT NONE
c      DOUBLE PRECISION u(97), c, cd, cm, uni, rand
      REAL u(97), c, cd, cm, uni, rand
      INTEGER i97, j97, ivec
      COMMON /raset1/ u, c, cd, cm, i97, j97, ivec 
c$OMP THREADPRIVATE (/raset1/)
      uni = u(i97) - u(j97)
      IF( uni .LT. 0.0 ) uni = uni + 1.0
      u(i97) = uni
      i97 = i97 - 1
      IF(i97 .EQ. 0) i97 = 97
      j97 = j97 - 1
      IF(j97 .EQ. 0) j97 = 97
      c = c - cd
      IF( c .LT. 0.0 ) c = c + cm
      uni = uni - c
      IF( uni .LT. 0.0 ) uni = uni + 1.0
      rand = uni
      RETURN
      END

      SUBROUTINE sobseq(n,x)
      INTEGER n,MAXBIT,MAXDIM
      REAL x(*)
      PARAMETER (MAXBIT=30,MAXDIM=6)
      INTEGER i,im,in,ipp,j,k,l,ip(MAXDIM),iu(MAXDIM,MAXBIT),iv(MAXBIT*
     *MAXDIM),ix(MAXDIM),mdeg(MAXDIM)
      REAL fac
      SAVE ip,mdeg,ix,iv,in,fac
      EQUIVALENCE (iv,iu)
      DATA ip /0,1,1,2,1,4/, mdeg /1,2,3,3,4,4/, ix /6*0/
      DATA iv /6*1,3,1,3,3,1,1,5,7,7,3,3,5,15,11,5,15,13,9,156*0/
      if (n.lt.0) then
        do 14 k=1,MAXDIM
          do 11 j=1,mdeg(k)
            iu(k,j)=iu(k,j)*2**(MAXBIT-j)
11        continue
          do 13 j=mdeg(k)+1,MAXBIT
            ipp=ip(k)
            i=iu(k,j-mdeg(k))
            i=ieor(i,i/2**mdeg(k))
            do 12 l=mdeg(k)-1,1,-1
              if(iand(ipp,1).ne.0)i=ieor(i,iu(k,j-l))
              ipp=ipp/2
12          continue
            iu(k,j)=i
13        continue
14      continue
        fac=1./2.**MAXBIT
        in=0
      else
        im=in
        do 15 j=1,MAXBIT
          if(iand(im,1).eq.0)goto 1
          im=im/2
15      continue
        pause 'MAXBIT too small in sobseq'
1       im=(j-1)*MAXDIM
        do 16 k=1,min(n,MAXDIM)
          ix(k)=ieor(ix(k),iv(im+k))
          x(k)=ix(k)*fac
16      continue
        in=in+1
      endif
      return
      END


