c-----------------------------------------------------------------------
      block data
c-----------------------------------------------------------------------
      real i0,i1,i2,j0,m0,m1,m2,k0,k1,k2
      common /coef/b0,b1,b2,b3,b4,c0,c1,c2,d0,a0,a1,a2,a3,a4,a5,
     * f0,f1,f2,f3,g0,g1,g2,i0,i1,i2,j0,m0,m1,m2,e0,e1,e2,e3,e4,
     * h0,h1,h2,h3,k0,k1,k2

c     coefficients
      data b0/+8.24493e-1/,
     * b1/-4.0899 e-3 /,
     * b2/ 7.6438 e-5 /,
     * b3/-8.2467 e-7 /,
     * b4/ 5.3875 e-9 /,
     * c0/-5.72466e-3 /,
     * c1/ 1.0227 e-4 /,
     * c2/-1.6546 e-6 /,
     * d0/ 4.8314 e-4 /,
     * a0/999.842594 /,
     * a1/ 6.793952e-2 /,
     * a2/-9.095290e-3 /,
     * a3/ 1.001685e-4 /,
     * a4/-1.120083e-6 /,
     * a5/ 6.536332e-9 /

      data f0/54.6746 /,
     * f1/-0.603459 /,
     * f2/ 1.09987 e-2 /,
     * f3/-6.1670  e-5 /,
     * g0/ 7.944  e-2 /,
     * g1/ 1.6483 e-2 /,
     * g2/-5.3009 e-4 /,
     * i0/ 2.2838 e-3 /,
     * i1/-1.0981 e-5 /,
     * i2/-1.6078 e-6 /,
     * j0/ 1.91075 e-4 /,
     * m0/-9.9348 e-7 /,
     * m1/ 2.0816 e-8 /,
     * m2/ 9.1697 e-10 /

      data e0/19652.21  /,
     * e1/ 148.4206 /,
     * e2/-2.327105 /,
     * e3/ 1.360477e-2 /,
     * e4/-5.155288e-5 /,
     * h0/ 3.239908 /,
     * h1/ 1.43713 e-3 /,
     * h2/ 1.16092 e-4 /,
     * h3/-5.77905 e-7 /,
     * k0/ 8.50935 e-5 /,
     * k1/-6.12293 e-6 /,
     * k2/ 5.2787  e-8 /


      end


c-----------------------------------------------------------------------
        subroutine shumtaup(shum,press,taup)
c----------------------------------------------------------------------

      real shum,press,taup
      real ee

        ee=press/100.*shum/(0.622+0.378*shum)
        if (ee.ge.6.176) then
          taup=873.23-27.7124*log(ee)
     &        -sqrt((873.23-27.7124*log(ee))**2 -374871.12)
        else
          taup=6141./(24.3-log(ee))
        endif

        end

c-----------------------------------------------------------------------
      subroutine deco1d2d (zf,zc,iwet,indend,fac)
c-----------------------------------------------------------------------
c     creates  2-d arrays from 1-d arrays
c-----------------------------------------------------------------------

      include 'C_model.f'
      parameter (khor1=khor+1)

      dimension iwet(khor1),indend(n)
      dimension zf(m,n),zc(khor)
#ifdef MPI
      common /cord/ I1,I2,J1,J2,IA1,IA2,JA1,JA2,IB1,IB2,JB1,JB2
     1,ID1,ID2,JD1,JD2
      common /lbkhor/ khorl,lzet(khor),lb0(n),le0(n),
     * lb1(n),le1(n),lb2(n),le2(n)
#endif
#ifdef MPI
      do k=JB1,JB2 
      do i=IB1,IB2
#else      
      do k=1,n
      do i=1,m
#endif 
         zf(i,k) = 0.0
      enddo
      enddo

      lwe = 0
#ifdef MPI
       do k=JB1,JB2
       lwa = lb2(k) 
       lwe = le2(k)
#else
      do k=1,n
        lwa = lwe+1
        lwe = indend(k)
#endif 
        do lw=lwa,lwe
          i  = iwet(lw)
          zz = zc(lw)
          zf(i,k) = zz*fac
        enddo
      enddo


      return
      end


c-----------------------------------------------------------------------
      subroutine deco1d2di (iz,izc,iwet,indend,ifac)
c-----------------------------------------------------------------------
c     creates  2-d arrays from 1-d arrays
c----------------------------------------------------------------------
      include 'C_model.f'
      parameter (khor1=khor+1)

      dimension iwet(khor1),indend(n)
      dimension iz(m,n),izc(khor)

      do k=1,n
      do i=1,m
         iz(i,k) = 0
      enddo
      enddo

      lwe = 0
      do k=1,n
        lwa = lwe+1
        lwe = indend(k)
        do lw=lwa,lwe
          i  = iwet(lw)
          izz = izc(lw)
          iz(i,k) = izz*ifac
        enddo
      enddo


      return
      end


c-----------------------------------------------------------------------
      subroutine comp2d1d (zic,zf)
c-----------------------------------------------------------------------
c    compress 2-d arrays (real)
c-----------------------------------------------------------------------

      include 'C_model.f'
      parameter(ilop1=ilo+1)
      parameter(khor1=khor+1)

      dimension zic(khor),zf(m,n)

      common /ind/ iwet(khor1),ldep(khor),lazc(khor),
     *  indend(n),isornr(n),isorsr(n),islab(n)

      lwe = 0
      nwet = 0
      do k=1,n
        lwa = lwe+1
        lwe = indend(k)
        do lw=lwa,lwe
          i = iwet(lw)
          zic(lw) = zf(i,k)
        enddo
      enddo


      return
      end



      subroutine readmet (feldncep,ncep,nincep,ijahr)
      dimension feldncep(ncep)

c      read(nincep,1000) iday,imonth,iyear,ihour,iminute
c       write(*,*) 'sread: Tag, Monat, Jahr, Stunde, Minute:', iday,
c     &               imonth,iyear,ihour,iminute,nincep

       if(ijahr.ge.1998)then
       do ii=1,ncep/6
         read(nincep,1001) (feldncep(n),n=(ii-1)*6+1,ii*6)
       enddo
       ihilf=ncep/6
       rhilf=real(ncep)/6.
       if (real(ihilf).ne.rhilf) then
         read(nincep,1001) (feldncep(n),n=ihilf*6+1,ncep)
       endif
      else
       do ii=1,ncep/6
         read(nincep,1002) (feldncep(n),n=(ii-1)*6+1,ii*6)
       enddo
       ihilf=ncep/6
       rhilf=real(ncep)/6.
       if (real(ihilf).ne.rhilf) then
         read(nincep,1002) (feldncep(n),n=ihilf*6+1,ncep)
       endif

	endif
 1000  format(2X,i2,1X,i2,1X,i2,4X,i2,1X,i2)
 1001   format(6e13.6)
 1002   format(6e13.5)
       return
       end



c-----------------------------------------------------------------------
      function esatw(tt)
c-----------------------------------------------------------------------
c     saturation of vapor pressure at water temperature t in Pa
c     Magnus equation
c-----------------------------------------------------------------------

      ai    = (7.5*tt)/(237.3+tt)
      esatw = 610.7*(10.**ai)


      return
      end


c-----------------------------------------------------------------------
      subroutine estaklim(saca,teca,nndrei)
c-----------------------------------------------------------------------
c     density calculation from T, S,
c     with function sigma (=unesco-formula)
c-----------------------------------------------------------------------
      include 'C_model.f'

      parameter(ilop1=ilo+1)
      parameter(khor1=khor+1)

      dimension saca(nndrei),teca(nndrei)

      common /ind/ iwet(khor1),ldep(khor),lazc(khor),
     * indend(n),isornr(n),isorsr(n),islab(n)
ccc
      common uc(ndrei),vc(ndrei),stc(ndrei),avc(ndrei),z(m,n)
      common zac(khor),wobc(khor),stuvc(khor),fricu(khor),cxc(khor)
      common cyc(khor),pac(khor),txc(khor),tyc(khor)
      common stpc(ndrei),sac(ndrei),tec(ndrei)
      common pres(ilo),wc(ndrei),fricv(khor)
      common /dreh/ sinfu(m),sinfv(m),cosfu(m),cosfv(m),sincx(m),
     * sincy(m),bx(m),by(m),pxu(m),pyv(m),pyu(m),pxv(m)
      common /gitter/ dt,r,g,dl,dlr,dlrh,dln(m),rdln(m),rad,dth,
     * dlvo(m),dlvu(m),gh,rdt,dt2,r4,coru(m),corv(m),dtrdln(m),dtdlr
      common/num/dc(ilo),av(ilo),ad(ilo),dh(ilo),pd(ilo),
     * prd(ilo),pr2d(ilo),r2d(ilo),tkw(ilo),tau(ilo),dd(ilo),
     * qa(ilop1),qbet(ilop1),qn(ilop1),rd(ilop1)
      common /met/ windx,windy,wlam,pnull,rhoq(ilo),refrho(ilo),stress


      nwet = 0
      imal=0
      lwe = 0
      do 160 k=1,n
      lwa = lwe+1
      lwe = indend(k)
      do 160 lw=lwa,lwe
      ldown = lazc(lw)

      dd(1)=pd(1)+zac(lw)
      do kk=2,ldown-1
        dd(kk)=pd(kk)
      enddo
      dd(ldown)=ldep(lw)

      stc(nwet+1) = sigma(saca(nwet+1),teca(nwet+1),pres(1))-refrho(1)
      do j=2,ldown
        stc(nwet+j) = sigma(saca(nwet+j),teca(nwet+j),pres(j))-refrho(j)
        stpc(nwet+j) = sigma(saca(nwet+j),teca(nwet+j),pres(j-1))
     *  -refrho(j-1)
      enddo

c      mixing with refernce densities
      do 172 ii=1,2
        inst=0
        do 171 j=2,ldown
          stcob=stc(nwet+j-1)+rhoq(j-1)
          stcun=stc(nwet+j)+rhoq(j)
          if(stcun.lt.stcob) then

            rhoob=stc(nwet+j-1)+rhoq(j-1)
            rhoun=stc(nwet+j)+rhoq(j)
            dmm=dd(j)*rhoun+dd(j-1)*rhoob

            sm=(saca(nwet+j-1)*dd(j-1)*rhoob+saca(nwet+j)*dd(j)*rhoun)
     *       /dmm
            tm=(teca(nwet+j-1)*dd(j-1)*rhoob+teca(nwet+j)*dd(j)*rhoun)
     *       /dmm

            saca(nwet+j)=sm
            teca(nwet+j)=tm
            saca(nwet+j-1)=sm
            teca(nwet+j-1)=tm

            stc(nwet+j) = sigma(saca(nwet+j),teca(nwet+j),pres(j))
     * 	      -refrho(j)
            stpc(nwet+j) = sigma(saca(nwet+j),teca(nwet+j),pres(j-1))
     *        -refrho(j-1)
            if(j.gt.2) then
              stpc(nwet+j-1) =
     *        sigma(saca(nwet+j-1),teca(nwet+j-1),pres(j-2))-refrho(j-2)
            endif
            stc(nwet+j-1) =
     *       sigma(saca(nwet+j-1),teca(nwet+j-1),pres(j-1))-refrho(j-1)

            inst=inst+1
          endif
  171   continue

        if(inst.le.0) goto 160
        imal=imal+1
  172 continue

  160 nwet = nwet+ldown
c      print*,imal,' unstable'


      return
      end



c-----------------------------------------------------------------------
      subroutine gittop(ltief,dzbod)
c-----------------------------------------------------------------------
      include 'C_model.f'

      parameter(khor1=khor+1)
      real dzbod(m,n)
      integer ltief(m,n)

      common /ind/ iwet(khor1),ldep(khor),lazc(khor),
     *  indend(n),isornr(n),isorsr(n),islab(n)

      do j=1,n
      do i=1,m
        dzbod(i,j) = 0.
        ltief(i,j) =0
      enddo
      enddo

      lwe = 0
      do j=1,n
        lwa = lwe+1
        lwe = indend(j)
        do lw=lwa,lwe
           i = iwet(lw)
           ltief(i,j) = lazc(lw)
           dzbod(i,j) = real(ldep(lw))
        enddo
      enddo


      return
      end








c-----------------------------------------------------------------------
      subroutine comp2d1dn (d2,d1)
c-----------------------------------------------------------------------
c     extend 2-d arrays to 1-d arrays with mwet grid points
c-----------------------------------------------------------------------
      include 'C_model.f'
      parameter(khor1=khor+1)

      dimension d2(m,n),d1(ndrei)
      common /ind/ iwet(khor1),ldep(khor),lazc(khor),
     *  indend(n),isornr(n),isorsr(n),islab(n)

      nwet = 0
      lwe  = 0

      do k=1,n
        lwa = lwe+1
        lwe = indend(k)
        do lw=lwa,lwe
          i = iwet(lw)
          ldown = lazc(lw)
          do j=1,1
            d2(i,k) = d1(nwet+j)
          enddo
          nwet = nwet+ldown
        enddo
      enddo


      return
      end






c-----------------------------------------------------------------------
      function rbackn(ta,ts,cn,eaa,sigem)
c-----------------------------------------------------------------------
c     longwave radiation after Bignami et al. 1995
c-----------------------------------------------------------------------

*-----lufttemperatur in kelvin
      t=273.15+ta
      ta4=t*t*t*t
*-----oberflaechentemperatur in kelvin
      t=273.15+ts
      ts4=t*t*t*t
*-----bewoelkung umrechnen, damit bei 8-achtel cn**2 = 1 ist
      cn1=cn / 8.
*-----emissivitaet wasser,eis
      emiss=0.97
*-----emissivitaet atmosphaere (clouds!)
      emisa=(0.653+0.00535*eaa)*(1.+0.1762*cn1**2)
*     stephan boltzmann konstante
      boltz=5.67e-8
      sigem = emiss*boltz
*-----langwellige nettostrahlung
      rbackn=sigem*ts4 - emisa*boltz*ta4


      return
      end


c-----------------------------------------------------------------------
      function rbacknn(ta,ts,cn)
c-----------------------------------------------------------------------
c     long wave radiation after Idso and Jackson
c-----------------------------------------------------------------------

*-----lufttemperatur in kelvin
      t=273.15+ta
      ta4=t*t*t*t
*-----oberflaechentemperatur in kelvin
      t=273.15+ts
      ts4=t*t*t*t
c      ea=ea/100.
*-----bewoelkung umrechnen, damit bei 8-achtel cn**2 = 1 ist
c      cn=1   !  wrong number , change!
      cn1=cn/8.
*-----emissivitaet wasser,eis
      emiss=0.97
*-----emissivitaet atmosphaere (clouds!)
      emisa=0.26*exp(-.77*(ta+0.15)**2/1000.)
      ffc=(1.-0.75*cn1)
*     stephan boltzmann konstante
      boltz=5.67e-8
      sigem = emiss*boltz
*-----langwellige nettostrahlung
      rbacknn=-emiss*boltz*(ta4*(1. - emisa*ffc)-ts4)

       return
       end




c-----------------------------------------------------------------------
      subroutine settoz(nlar,ngro,large,gross)
c-----------------------------------------------------------------------
c     dummy arrays set to zero
c-----------------------------------------------------------------------

      dimension large(nlar),gross(ngro)
 
      do l = 1,nlar
        large(l) = 0
      enddo
      do l = 1,ngro
        gross(l) = 0.0
      enddo


      return
      end

c-----------------------------------------------------------------------
      subroutine settime(time,ijulu,ihouu,iminu,isecu)
c-----------------------------------------------------------------------
c     time set to rhour,rminu,rsecu
c-----------------------------------------------------------------------

      itime = int(time)
      isecu = mod(itime,60)
      iminu = mod((itime-isecu),3600)/60
      ihouu = mod((itime-(iminu*60)-isecu),86400)/3600
      ijulu = (itime-(ihouu*3600)-(iminu*60)-isecu)/86400+1
      return
      end
cc-----------------------------------------------------------
c      FUNCTION julian_date (yyyy, mm, dd) RESULT (julian)
cc-------------------------------------------------
c      IMPLICIT NONE
c! converts calendar date to Julian date
c! cf Fliegel & Van Flandern, CACM 11(10):657, 1968
c! example: julian_date(1970,1,1)=2440588
c      INTEGER,INTENT(IN) :: yyyy,mm,dd
c       INTEGER :: julian
c      julian = dd-32075+1461*(yyyy+4800+(mm-14)/12)/4 + 
c     & 367*(mm-2-((mm-14)/12)*12)/12- 
c     &3*((yyyy + 4900 + (mm - 14)/12)/100)/4
c      END FUNCTION julian_date
c----------------------------------------------------------------------
      subroutine intpork(weight,nrmess,rlamu,rphiu,ifunc,
     &                    param,dmax,ianzpkt,nera)
c-----------------------------------------------------------------------
      include 'C_model.f'
      include 'C_mpi.f'
      parameter (khor1=m*n,intpmax=9)

      real weight(m,n,ianzpkt)
      integer nrmess(m,n,ianzpkt)
      real weightv(ianzpkt, khor1)
      integer nrmessv(ianzpkt, khor1)
      real lambdamess(nera),phimess(nera),
     &     lambdaziel(khor1),phiziel(khor1)
      real rlamu(nera),rphiu(nera)
      real param,dmax
      integer ifunc

      common /gitter2/ fianf,yambdanf,dphi,dlambda,rearth
     
c        do i=1,m
c        do j=1,n
c       read(nfumet_in,*) (weight(i,j,jj),jj=1,ianzpkt)
c        read(nfumet_in,*) (nrmess(i,j,jj),jj=1,ianzpkt)
c        enddo
c        enddo
c	  GOTO 150

      nziel=0
      pi=3.1415927
      grdrad = pi/180.
      yambog=yambdanf*grdrad
      fibog=fianf*grdrad
       IF(myid .EQ. 0)THEN
       print*,'lambda, phi, dlambda, dphi'
       write(*,*) yambog/grdrad,fibog/grdrad,
     *          dlambda/grdrad,dphi/grdrad
       ENDIF
      do i=1,m
        do j=1,n
          nziel = nziel+1
          lambdaziel(nziel) = yambog+real(j-1)*dlambda
          phiziel(nziel)    = fibog-real(i-1)*dphi
        enddo
      enddo

      do i=1,nera
        lambdamess(i) = rlamu(i)*grdrad
        phimess(i)    = rphiu(i)*grdrad
      enddo

c      do i=1,nera
c      write(*,77) lambdamess(i)/grdrad,phimess(i)/grdrad
c      enddo
c 77   format(2e12.3)

c      do i=1,khor1
c      write(*,77) lambdaziel(i)*180.0/pi,phiziel(i)*180/pi
c      enddo
c      pause


      call ipol2dkg(weightv,nrmessv,lambdaziel,phiziel,lambdamess,
     &      phimess,khor1,nera,rearth,ianzpkt,dmax,ifunc,param)


c     write(*,*) weightv
c     write(*,*) nrmessv

      nziel=0
      do i=1,m
        do j=1,n
          nziel = nziel+1
          do k=1,ianzpkt
            weight(i,j,k) = weightv(k,nziel)
            nrmess(i,j,k) = nrmessv(k,nziel)
          enddo
        enddo
      enddo

c        do i=1,m
c        do j=1,n
c         write(nfumet_in,*) (weight(i,j,jj),jj=1,ianzpkt)
c         write(nfumet_in,*) (nrmess(i,j,jj),jj=1,ianzpkt)
c         enddo
c         enddo
c	     close(nfumet_in)
c		stop

150      continue
c      close(nfumet_in)

	  return
      end

c-----------------------------------------------------------------------
      subroutine ipol2dkg(weight,nrmess,lambdaziel,phiziel,lambdamess,
     &                 phimess,nziel,nmess,r,ianzpkt,dmax,ifunc,param)
c-----------------------------------------------------------------------

c    Interpolation of measured values on pre-defined grid
c        (e.g. model grid)
c    here to be used to interpolate atmospheric boundary conditions
c    to model grid
c    input values and output grid are given spheric coordinates
c
c    Input:
c            nziel        : number of final points
c            lambdaziel(i): (geogr.) longitude of final points i
c            phiziel(i)   : (geograph.) latitude of final points i
c
c            nmess         : number of input grid points
c            lambdamess(j) : (geograph.) longitude of input points j
c            phimess(j)    : (geograph.) latitude of input points j
c
c            r         : Radius der Kugel
c
c            ianzpkt   : Anzahl der Messpunkte, die fuer 1 Messpunkt
c                        maximal zur Interpolation verwendet
c                        werden sollen
c            dmax      : ist der Abstand eines Messpunktes vom
c                        Zielpunkt >dmax, wird der Messpunkt bei
c                        der Interpolation nicht mehr berueck-
c                        sichtigt
c
c            ifunc     : Kennung fuer die Gewichtungsfunktion
c
c                      =1: Gauss-Funktion
c                              param: Standardabweichung
c                      =2: 1/x
c                             param: unbenutzt
c                      =3: 1/x^2
c                             param: unbenutzt
c                      =4: konstant
c                             param: unbenutzt
c                      =5: linear
c                             param: unbenutzt
c
c    Output:
c            weight(j,i) :    Gewichtung des j-naechsten Messpunktes bzgl.
c                             Zielpunkt i
c            nrmess(j,i) :    Index des j-naechsten Messpunktes bzgl.
c                             Zielpunkt i
c
c    Frank Siegismund, 25.09.98
c    H. Bockelmann (DKRZ), 24.07.2018 : cache optimized version implemented
c
      include 'C_mpi.f'
      integer nziel,nmess,ianzpkt
      integer ifunc
      real param
      real weightsum
      real r,pi
      real xn,yn,zn
      real d,dmax,inv_dmax
      integer nrmess(ianzpkt, nziel)
      real weight(ianzpkt, nziel)
      real lambdaziel(nziel),phiziel(nziel)
      real lambdamess(nmess),phimess(nmess)
      real dzielmess(nmess)
      integer nrmerk(nmess)

      real, dimension(nmess) :: coscos_mess, cossin_mess, sin_mess
      real                   :: coscos_ziel, cossin_ziel, sin_ziel

      data pi/3.1415927/

      
      if (nziel.le.0) goto 9990
      if (nmess.le.0) goto 9992

      grdrad = pi/180.
      inv_dmax = 1./dmax

c precalculte phimess, lambdamess stuff
      do j = 1,nmess
        coscos_mess(j) = cos(phimess(j))*cos(lambdamess(j))
        cossin_mess(j) = cos(phimess(j))*sin(lambdamess(j))
        sin_mess(j)    = sin(phimess(j))
      enddo

c loop over all final points
      do 10 i = 1,nziel
        n = 0

        coscos_ziel = cos(phiziel(i))*cos(lambdaziel(i))
        cossin_ziel = cos(phiziel(i))*sin(lambdaziel(i))
        sin_ziel    = sin(phiziel(i))

        do 5 j=1,nmess
          xn = coscos_mess(j) - coscos_ziel
          yn = cossin_mess(j) - cossin_ziel
          zn = sin_mess(j)    - sin_ziel

          d = r*sqrt(xn*xn+yn*yn+zn*zn)
          if (d.le.dmax) then
            n = n+1
            nrmerk(n) = j
            dzielmess(n) = d
          endif
5       continue

        if (n.eq.0) then
          IF(myid .EQ. 0) THEN
            write(*,*) 'Fuer den Zielpunkt phi=',phiziel(i)/grdrad,
     &',lambda=', lambdaziel(i)/grdrad,' keinen Messpunkt gefunden!'
          ENDIF
          do 6 k=1,ianzpkt
            nrmess(k,i) = 9999999.
            weight(k,i) = 0.
6         continue
          goto 10
        endif

        call sort(dzielmess,nrmerk,n,-1)

        weightsum=0
c       limit the number of measurement points from n to ianzpkt
        if (n.gt.ianzpkt) n=ianzpkt

        do 7 k=1,n
           nrmess(k,i) = nrmerk(k)

           select case (ifunc)
c            Gauss-Funktion
             case (1)
               weight(k,i) = exp(-(dzielmess(k)/param)**2)
c            1/x
             case (2)
               if (dzielmess(k).lt.1E-10) then
                 weight(k,i) = 1E10
               else
                 weight(k,i) = 1/dzielmess(k)
               endif
c            1/x^2
             case (3)
               if (dzielmess(k).lt.1E-5) then
                 weight(k,i) = 1E10
               else
                 weight(k,i) = 1/dzielmess(k)**2
               endif
c            konstant
             case (4)
               weight(k,i) = 1
c            linear
             case default
               weight(k,i) = 1-dzielmess(k)/dmax
          end select

          weightsum = weightsum+weight(k,i)
7       continue

        if (weightsum.lt.1E-30) then
          IF(myid .EQ. 0)THEN
            write(*,*) 'Fuer den Zielpunkt phi=',phiziel(i),', lambda=',
     &             lambdaziel(i),' zu kleine Gewichtungen!'
          ENDIF
          do 77 k=1,ianzpkt
            nrmess(k,i)=999999
            weight(k,i)=0
77        continue
          goto 10
        endif

        do 8 k=1,n
          weight(k,i) = weight(k,i)/weightsum
8       continue

        if (n.lt.ianzpkt) then
          do 9 k=n+1,ianzpkt
            nrmess(k,i)=8888888
            weight(k,i)=0
9         continue
        endif
10    continue

      goto 9999

9990  write(*,*) 'Keine Zielpunkte vorgegeben!'
      goto 9999

9992  write(*,*) 'Keine Messpunkte vorgegeben!'
      goto 9999

9999  end


c-----------------------------------------------------------------------
      subroutine sort(zahl,verknpf,n,imode)
c-----------------------------------------------------------------------

c    Sortieren des n-elementigen Vektors zahl. Der mit zahl verkn201pfte
c    Vektor verknpf wird mitsortiert.
c
c    imod >= 0 : zahl wird absteigend sortiert.
c    imod <  0 : zahl wird aufsteigend sortiert.
c
c    Achtung!!: verknpf ist ein Vektor ganzer Zahlen.
c
c    Frank Siegismund, 19.08.98
      include 'C_mpi.f'
      integer n,nn,imode
      real zahl(n),czahl
      integer verknpf(n),cverknpf

      if (n.lt.1) goto 9990
      if (n.eq.1) goto 9995

      if (imode.lt.0) then
         do 1 i=1,n
            zahl(i)=-zahl(i)
1        continue
      endif

      k=0
      nn=n
2     nn=nn/2
      k=k+1
      if (nn.gt.1) goto 2
      if (2**k.lt.n) k=k+1


      do 30 i=1,k
         ig=n/2**i
         if (ig.eq.0) goto 22
         do 20 j=1,ig
            ii=(j-1)*2**i+1
            do 10 l=2**(i-1) +1+(j-1)*2**i,j*2**i
               czahl=zahl(l)
               cverknpf=verknpf(l)
               do 5 ll=ii,j*2**i
                  if (czahl.gt. zahl(ll)) then
                     if (ll.gt.l) goto 10
                     do 3 jj=l-1,ll,-1
                        zahl(jj+1)=zahl(jj)
                        verknpf(jj+1)=verknpf(jj)
3                    continue
                     zahl(ll)=czahl
                     verknpf(ll)=cverknpf
                     ii=ll+1
                     goto 10
                  endif
5              continue
10          continue
20       continue
22       if (n-ig*2**i.gt.2**(i-1)) then
            ii=ig*2**i+1
            do 28 l=2**(i-1)+1+ig*2**i,n
               czahl=zahl(l)
               cverknpf=verknpf(l)
               do 26 ll=ii,n
                  if (czahl.gt.zahl(ll))then
                     if (ll.gt.l) goto 28
                     do 24 jj=l-1,ll,-1
                        zahl(jj+1)=zahl(jj)
                        verknpf(jj+1)=verknpf(jj)
24                   continue
                     zahl(ll)=czahl
                     verknpf(ll)=cverknpf
                     ii=ll+1
                     goto 28
                  endif
26             continue
28          continue
         endif
30    continue

      if (imode.lt.0) then
         do 40 i=1,n
            zahl(i)=-zahl(i)
40       continue
      endif

      goto 9999

9990  continue
      IF(myid .EQ. 0)THEN
      write(*,*) n,' Zahlen sortieren?'
      endif
      goto 9999

9995  continue
c      write(*,*) 'n=1'
      goto 9999

9999  end


c-----------------------------------------------------------------------
      subroutine interpol2(vargitt,varmess,weight,nrmess,ianzpkt,nera)
c-----------------------------------------------------------------------
      include 'C_model.f'
c      include 'C_mpi.f'
      real vargitt(m,n),varmess(nera)
      real weight(m,n,ianzpkt)
      integer nrmess(m,n,ianzpkt)
#ifdef MPI
      common /cord/ I1,I2,J1,J2,IA1,IA2,JA1,JA2,IB1,IB2,JB1,JB2
     1,ID1,ID2,JD1,JD2
#endif 

#ifdef MPI
      do j=JD1,JD2
      do i=ID1,ID2
c      do i=1,m
c      do j=1,n
#else
      do i=1,m
      do j=1,n
#endif
        vargitt(i,j)=0
        do k=1,ianzpkt
        if(nrmess(i,j,k).lt.8888888) then
         vargitt(i,j)=vargitt(i,j)+weight(i,j,k)*varmess(nrmess(i,j,k))
         end if
        enddo
      enddo
      enddo
      end


c-----------------------------------------------------------------------
      Subroutine uv2sr(u,v,speed,dir,lrev,imet)
c-----------------------------------------------------------------------
c
c     input,output:
c
c     u,v       components of current
c     speed     amount     of current
c     dir       direction  of current
c
c     switches:
c
c     lrev      = 0 => transform u,v to speed,dir
c     lrev      = 1 => transform speed,dir to u,v
c
c     imet      = 0 => water-current
c     imet      = 1 => wind-current
c
c     U. Huebner 22.10.98   last modification
c
c-----------------------------------------------------------------------

c      double precision u,v,speed,dir,x,y,rr
c      double precision xdir,xspeed
      real u,v,speed,dir,x,y,rr
      real xdir,xspeed

      pi = 4.*atan(1.)
      EPS= 1.E-10
      FAK= 180./pi

      if(lrev.eq.0)then
        x = u
        y = v

        speed = sqrt(x*x+y*y)

        su = sign(1.,x)
        sv = sign(1.,y)
        if( abs(x) .le. EPS) goto 11
        if( abs(y) .le. EPS) goto 12
        rr = atan(x/y) * fak + 90.*(2.-su*(1.+sv))
        goto 14
 11     if( abs(y) .le. EPS) goto 13
        rr = 90. * (1.-sv)
        goto 14
 12     rr = 90. * (2.-su)
        goto 14
 13     rr = 0
 14     continue

        if(imet.eq.1)then
          rr = rr+180.0
          if(rr.ge.360.0)rr=rr-360.0
          if(x.eq.0.and.y.eq.0)rr=0.0
        endif


        dir = rr

      else

        xdir = dir
        xspeed = speed

        if(imet.eq.1)then
          xdir = xdir-180.0
        endif
        xdir = xdir/FAK

        v = xspeed*cos(xdir)
        u = xspeed*sin(xdir)

        if(abs(u).le.1.0e-6)u=0.0
        if(abs(v).lt.1.0e-6)v=0.0

      endif

      return
      end


c-----------------------------------------------------------------------
      subroutine xstat(para,ianz,pmax,pmin,pmean,psum)
c-----------------------------------------------------------------------

      dimension para(ianz)

      pmax  = -9999999.9
      pmin  =  9999999.9
      psum  = 0.0
      pmean = 0.0

      do i=1,ianz
        pmax = amax1(para(i),pmax)
        pmin = amin1(para(i),pmin)
        psum = psum+para(i)
      enddo
      pmean = psum/ianz

      return
      end

c-----------------------------------------------------------------------
      subroutine xstat2(para,m,n,jjc,pmax,pmin,pmean,psum)
c-----------------------------------------------------------------------

      dimension para(m,n),jjc(m,n)

      pmax  = -9999999.9
      pmin  =  9999999.9
      psum  = 0.0
      pmean = 0.0
      ianz  = 0

      do i=1,m
      do k=1,n
      do j=1,jjc(i,k)
        ianz = ianz+1
        pmax = amax1(para(i,k),pmax)
        pmin = amin1(para(i,k),pmin)
        psum = psum+para(i,k)
      enddo
      enddo
      enddo
      pmean = psum/ianz

      return
      end


      Subroutine Rzero(a,n)
      dimension a(n)
      do 1 i=1,n
1     a(i)=0.
      return
      end


      Subroutine Rzero999(a,n)
      dimension a(n)
      do 1 i=1,n
1     a(i)=-999.
      return
      end

c-----------------------------------------------------------------------
      subroutine estate
c-----------------------------------------------------------------------
c     berechnung der dichte aus salz und temperatur
c     mit function sigma (=unesco-formel)
c-----------------------------------------------------------------------
      include 'C_model.f'
      include 'C_index.f'

      parameter(ilop1=ilo+1)
      parameter(khor1=khor+1)

      common /ind/ iwet(khor1),ldep(khor),lazc(khor),
     *   indend(n),isornr(n),isorsr(n),islab(n)
      common uc(ndrei),vc(ndrei),stc(ndrei),avc(ndrei),z(m,n)
      common zac(khor),wobc(khor),stuvc(khor),fricu(khor),cxc(khor)
      common cyc(khor),pac(khor),txc(khor),tyc(khor)
      common stpc(ndrei),sac(ndrei),tec(ndrei)
      common pres(ilo),wc(ndrei),fricv(khor)
      common /dreh/ sinfu(m),sinfv(m),cosfu(m),cosfv(m),sincx(m),
     *   sincy(m),bx(m),by(m),pxu(m),pyv(m),pyu(m),pxv(m)
      common /gitter/ dt,r,g,dl,dlr,dlrh,dln(m),rdln(m),rad,dth,
     *  dlvo(m),dlvu(m),gh,rdt,dt2,r4,coru(m),corv(m),dtrdln(m),dtdlr
      common/num/dc(ilo),av(ilo),ad(ilo),dh(ilo),pd(ilo),
     *   prd(ilo),pr2d(ilo),r2d(ilo),tkw(ilo),tau(ilo),dd(ilo),
     *   qa(ilop1),qbet(ilop1),qn(ilop1),rd(ilop1)
      common /met/ windx,windy,wlam,pnull,rhoq(ilo),refrho(ilo),stress
#ifdef MPI
      include 'C_mpi.f'
      common /cord/ I1,I2,J1,J2,IA1,IA2,JA1,JA2,IB1,IB2,JB1,JB2
     1,ID1,ID2,JD1,JD2
      common /icecor/ ICEI1,ICEI2,ICEJ1,ICEJ2,ICEIA1,ICEIA2,ICEJA1,
     * ICEJA2,ICEX1,ICEX2,ICEY1,ICEY2,ICEP1,ICEP2,ICEQ1,ICEQ2,ICEX3,
     * ICEX4,ICEY3,ICEY4,ICEX5,ICEX6,ICEY5,ICEY6
      common /lbkhor/ khorl,lzet(khor),lb0(n),le0(n),
     * lb1(n),le1(n),lb2(n),le2(n)
#endif
#ifdef MPI 
      do j=JB1,JB2 
         lwa = lb2(j)
         lwe = le2(j)
#else
         lwa = lb(1)
         lwe = le(n)
#endif
         if (lwa.le.lwe) then
            llb = indwet(lwa)+1
            lle = indwet(lwe)+lazc(lwe)
!CDIR IEXPAND(SIGMA)
!CDIR VECTOR
            do ll = llb,lle
               k = indver(ll)
               stc(ll) = sigma(sac(ll),tec(ll),pres(k))-refrho(k)
               if (k.gt.1) stpc(ll) = sigma(sac(ll),tec(ll),pres(k-1))-
     1          refrho(k-1)
            end do
         end if
#ifdef MPI
      end do
#endif
      return
      end

c-----------------------------------------------------------------------
      function sigma (s,t,p)
c-----------------------------------------------------------------------
c     testrechnung mit autodbl(dblpad) ist erfolgreich gewesen
c     testdaten:  s=0. t=25. p=1000.  rho=1037.90204
c                 s=35.t=5.  p=0.     rho=1027.67547
c                 s=35.t=25. p=0.     rho=1023.34306
c                 s=35.t=25. p=1000.  rho=1062.53817
c                 kai jancke13.01.83
c-----------------------------------------------------------------------
c
c                       i n p u t
c
c **********   p must be given as depth in meters * 0.1  **********
c **********   s must be given as salinity in promille   **********
c **********   t must be given as degree centigrade      **********
c
c-----------------------------------------------------------------------

      real i0,i1,i2,j0,m0,m1,m2,k0,k1,k2,ksbm
      common /coef/b0,b1,b2,b3,b4,c0,c1,c2,d0,a0,a1,a2,a3,a4,a5
     ,,f0,f1,f2,f3,g0,g1,g2,i0,i1,i2,j0,m0,m1,m2,e0,e1,e2,e3,e4
     ,,h0,h1,h2,h3,k0,k1,k2
      s2=s*s
      s3=s2*s
      t2=t*t
      t3=t2*t
      t4=t2*t2

*************** density of reference pure water (snow)
***************  density at standard atmosphere (p=0)

      rhonul=   a1*t+a2*t2+a3*t3+a4*t4+a5*t4*t
     +           +(b0+b1*t+b2*t2+b3*t3+b4*t4)*s
     +           +(c0+c1*t+c2*t2) *sqrt(s3)+d0*s2
***************  secant bulk modulus
      ksbm =  e0+e1*t+e2*t2+e3*t3+e4*t4
     +          + (f0+f1*t+f2*t2+f3*t3)*s + (g0+g1*t+g2*t2)*sqrt(s3)
     +     + (h0+h1*t+h2*t2+h3*t3
     +          + (i0+i1*t+i2*t2)*s       + j0*sqrt(s3)) * p
     +     + (k0+k1*t+k2*t2
     +          + (m0+m1*t+m2*t2)*s) * p * p

      sigma  = 1.e-3*(rhonul+a0*p/ksbm)/(1.-p/ksbm)


      return
      end


c-----------------------------------------------------------------------
      subroutine plant (zf,zc,ic,iwet,indend,fac,lint)
c-----------------------------------------------------------------------
      include 'C_model.f'
      parameter(khor1=khor+1)

      dimension ic(khor),iwet(khor1),indend(n)
      dimension izetf(m,n),zf(m,n),zc(khor)
      ifac = fac+0.5
      if (ifac.le.0) ifac = 1
c     ---------- plint kann dreierlei verschiedenes ---------
c       lint =        -1            0            +1
c     ausdruck =  int(compr),  real(compr),   real(m,n)
c     -------------------------------------------------------
      print*,'dies ist der Faktor',fac
      do 9 k = 1,n
      do 9 i = 1,m
    9 izetf(i,k) = 88888
      idru = 2
      lwe = 0
      do 15 k = 1,n
      lwa = lwe+1
      lwe = indend(k)
      do 15 lw = lwa,lwe
      i = iwet(lw)
      if (lint) 11,12,13
   11 izetf(i,k) = ic(lw)/ifac
      goto 15
   12 zz = zc(lw)*fac
      goto 14
   13 zz = zf(i,k)*fac
   14 izetf(i,k) = nint(zz)
   15 continue
      na = 1
      me = m
      ne = 29
      do 5 idr = 1,idru
      write(6,600) (k,k = na,ne)
      write(7,600) (k,k = na,ne)
      do 4 i = 1,me
      write(6,601) i,(izetf(i,k),k = na,ne) ,i
  4   write(7,601) i,(izetf(i,k),k = na,ne) ,i
      write(6,600) (k,k=na,ne)
      write(7,600) (k,k = na,ne)
      na = ne+1
      ne = na+28
      ne=min(ne,n)
  5   continue

  600 format (/6x,29(i3,1x)/)
  601 format (1x,i3,1x,30i4)


      return
      end



      subroutine IO_daily_physics(nout,mjar,lmon,nday,ivier,ip,nnn) ! nnn write=1 or read=2
      include 'C_model.f'
      include 'C_daily_physics.f'
      common /radis/fqgmit(m,n),fqrmit(khor),fqsmit(khor),fqlmit(khor)

      if(nnn.eq.1)then
          write (nout) mjar,lmon,nday,ivier,ip
          write (nout) zmit
          write (nout) umit
          write (nout) vmit
          write (nout) wcmit
          write (nout) acmit
          write (nout) szmit
            write (nout) tcmit
            write (nout) scmit
          write(nout) frimit
          write(nout) hismit
          write(nout) hisrmit
          write(nout) tismit
          write(nout) uimit
          write(nout) vimit
            write(nout) fqgmit
            write(nout) fqrmit
            write(nout) fqsmit
            write(nout) fqlmit
          write(nout) qois
          write(nout) qiis
	else !nnn=2
          read (nout) mjarcc,lmoncc,ndaycc,iviercc,ipcc
          read (nout) zmit
          read (nout) umit
          read (nout) vmit
          read (nout) wcmit
          read (nout) acmit
          read (nout) szmit
            read (nout) tcmit
            read (nout) scmit
          read(nout) frimit
          read(nout) hismit
          read(nout) hisrmit
          read(nout) tismit
          read(nout) uimit
          read(nout) vimit
            read(nout) fqgmit
            read(nout) fqrmit
            read(nout) fqsmit
            read(nout) fqlmit
          read(nout) qois
          read(nout) qiis
	endif  !     if(nnn.eq.1)then

	return
	end

      subroutine IO_start_physics(nout,mjar,lmon,nday,ivier,ip,nnn)
      include 'C_model.f'
      parameter (nx=m,ny=n,ilop1=ilo+1,khor1=khor+1)
      include 'C_ice.f'
      common /urand/       zday(lrp),zvar(lrp),zalt(lrp)

      if(nnn.eq.1)then
          write (nout) ip,mjar,lmon,nday
          write (nout) zalt,z,zac
          write (nout) txc,tyc,pac
          write (nout) uc
          write (nout) vc
          write (nout) wc
          write (nout) avc
          write (nout) tec
          write (nout) sac
          write (nout) frice,his,hisr,tis
          write (nout) ui,vi
      else !nnn=2
          read (nout) imalcc,mjarcc,monatcc,ltagcc
          read (nout) zalt,z,zac
          read (nout) txc,tyc,pac
          read (nout) uc
          read (nout) vc
          read (nout) wc
          read (nout) avc
          read (nout) tec
          read (nout) sac
          read (nout) frice,his,hisr,tis
          read (nout) ui,vi
      endif
	return
	end



c-----------------------------------------------------------------------
      subroutine feld(bewoel,einstr,einstrice,
     * mjar,time,
     * itest,jjc)

*	call feld(dumm,einstr,einstrice, !===> incoming radiation [W/m**2]
*     *mjar,time,
*     *itest,jjc)
c-----------------------------------------------------------------------
c   subroutine providing solar radition to water and ice
c   feld doing data handling for the radatm and radice calls
c-----------------------------------------------------------------------
      include 'C_model.f'
      include 'C_index.f'

      parameter (nstat=29,khor1=khor+1)

      dimension einstrice(m,n)
      dimension einstr(m,n),bewoel(m,n)
      dimension jjc(m,n)
      dimension para(13)

      integer itage(1),itagen(13)
      character*80 text

      common /wetpos/ indx(nstat),indy(nstat)
      common /ind/ iwet(khor1),ldep(khor),lazc(khor),
     *   indend(n),isornr(n),isorsr(n),islab(n)
      common /radi/fqg(m,n),fqr(khor),fqs(khor),fql(khor)
	include 'C_xy.f'
#ifdef MPI
       include 'C_mpi.f'
       common  /cord/ I1,I2,J1,J2,IA1,IA2,JA1,JA2,IB1,IB2,JB1,JB2
     1,ID1,ID2,JD1,JD2
       common  /lbkhor/ khorl,lzet(khor),lb0(n),le0(n),
     * lb1(n),le1(n),lb2(n),le2(n)
#endif
      ldim= 13

      para(2)=mjar
      para(3)=86400.
      para(4)=1368.
      para(5)=1
      para(6)=1
      para(7)=1
      para(8)=1
      para(9)=1
      para(10)=1
      para(11)=1
      para(12)=1
      para(13)=1

      in = 1
cval      call setita(in,para(2),itage) ! days in current year
      call setita(para(2),itage(1))
#ifdef MPI
      do j=JB1,JB2
         do lw=lb2(j),le2(j)
#else  
      do j=1,n
         do lw=lb(j),le(j)
#endif
            i = iwet(lw)
            itage(1)=366
            para(1) = xt(i)   ! latitude position
            bew = bewoel(i,j)
c            ea = taup(i,j)*esatw(tl(i,j))/100.
!CDIR NEXPAND
            call radatm(1,ldim,time,bew,para,itage,ein,ea)
            einice = 0.
            bew = bewoel(i,j)
!CDIR NEXPAND
            call radice(1,ldim,time,bew,para,itage,einice)
            einstr(i,j) = ein
            fqg(i,j) = ein
            einstrice(i,j) = einice
         enddo
      enddo


      return
      end
c-----------------------------------------------------------------------
      subroutine radatm(in,ldim,time,bewoel,para,itage,einstr,ea)
c-----------------------------------------------------------------------
c      global radiation according to Dobson and Smith
c        journal of the royal meteorological society, 117, pp845-851.
c-----------------------------------------------------------------------
c parameter (in):
c     1) in              = maximale anzahl der jahre
c     2) time            = aktuelle modellzeit (in sec, std oder tagen)
c     3) ldim            = dimension des parameter vektors para(ldim)
c     4) bewoel          = aktueller bewoelkungsgrad (achtel)
c     5) para(ldim)      = parametervektor
c             para(1)  = breite in grad
c             para(2)  = jahr der simulation
c             para(3)  = umrechnung der zeit auf tage (z.b. 86400. fuer
c                        sekunden)
c             para(4)  = solarkonstante (w / m**2)
c             para(5)  = modellauswahl fuer subroutine param
c                      = 0 --> wird nicht beruecksichtigt
c                      = 1 --> dobson & smith (1985)
c                      = 2 --> dobson & smith (1988)
c                      = 3 --> dobson & smith (1989)
c             para(6)  = modellauswahl fuer subroutine tagums
c                      = 0 --> wird nicht beruecksichtigt
c                      = 1 --> dobson & smith (1985 - 1989)
c             para(7)  = modellauswahl fuer subroutine deklin
c             para(8)  = modellauswahl fuer subroutine taglng
c                      = 0 --> wird nicht beruecksichtigt
c                      = 1 --> dobson & smith (1985 - 1989)
c             para(9)  = modellauswahl fuer subroutine stdwnk
c            para(10)  = modellauswahl fuer subroutine sonhoe
c                      = 0 --> wird nicht beruecksichtigt
c                      = 1 --> dobson & smith (1985 - 1989)
c            para(11)  = modellauswahl fuer subroutine meteor
c                      = 0 --> wird nicht beruecksichtigt
c                      = 1 --> dobson & smith (1985 -1989)
c            para(12)  = modellauswahl fuer subroutine water
c                      = 0 --> wird nicht beruecksichtigt
c                      = 1 --> becker (1981)
c                      = 2 --> nach payne fuer breite = 45 grad
c            para(13)  = modellauswahl fuer strahlungsformel
c                      = 1 --> dobson & smith (1985 - 1989)
c
c parameter (out):
c     1) itage(in)       = vektor mit anzahl der tage jedes jahres
c     2) einstr          = simulationsergebnis
c     2) tagums --> setzt den anfangstag entsprechend den modell-
c                   erfordernissen
c     3) astro  --> berechnet die astronomischen groessen der modelle,
c                   ruft seinerseits auf:
c        4) deklin --> berechnet den deklinationswinkel der sonnenbahn
c        5) taglng --> berechnet die tageslaenge (sonnenunter-aufgang)
c        6) stdwnk --> berechnet den jeweiligen stundenwinkel
c        7) sonhoe --> berechnet den sinus des sonnenwinkels
c     8) meteor --> simuliert die meteorologischen einfluesse
c                   --> abhaengigkeit von der bewoelkung
c     9) water  --> simuliert die einfluesse der meeresoberflaeche
c                   --> albedo
c-----------------------------------------------------------------------

      parameter (jp=9)
      dimension itage(in)
      dimension a(jp) , b(jp) , para(ldim)
      common /atm/ rnil,istop
c      data istop/1/
      ijahr = 1
      sinsow = 0.
      xtag = 0.
      akoeff = 0.
      bkoeff = 0.
      pi = atan(1.) * 4.
       cn1=bewoel/8.
cc      if(istop.eq.1) then
cc     write(6,1201)
 1201 format(////1x,'ausdruck der steuerparameter des strahlungsmodells:
     & ',//)
cc      do 23 i=1,ldim
cc        write(6,1200) i,para(i)
cc  23  continue
c c     end if
 1200 format(1x,'para(',i2,') =',f15.3)

cc      if(para(1).lt.0.or.para(1).gt.90) stop 'breite out of range]'
      breite = grdbom(para(1))
      xjahr = para(2)
      itfix = ifix(para(3))

*-----belegung des vektors itage(in) mit anzahl der tage / jahr
ccc      call setita(in,xjahr,itage)
cvalccc      call setita(xjahr,itage(1))
ccc      if(istop.eq.1) print*,'itage (anzahl der tage im jahr)= ',itage(1)

*-----belegung der vektoren a() und b()
      call param(jp,a,b,ifix(para(5)))
c     if(istop.eq.1) write(6,1000) para(1)

************************************************************************
************************************************************************
*                                                                      *
*     berechnung der einstrahlung                                      *
*                                                                      *
************************************************************************
************************************************************************


*-------feststellung des aktuellen tages
        jtag = ifix(time / float(itfix)) + 1
*-->    umsetzen der tage
        call tagums(in,jtag,ijahr,itage,xtag,ifix(para(6)))
*-->    berechnung der tageslaenge und des sinus des sonnenwinkels
        call astro(time,itfix,ldim,breite,xtag,itage(ijahr),para,taglg,
     &             sinsow)
*---------feststelluing des stundenabschnittes eines tages
*-->      simulation der meteorologischen einfluesse
          call meteor(jp,a,b,bewoel,akoeff,bkoeff,ifix(para(11)))
*-->      simulation der einfluesse der meeresoberflaeche
          call water(jtag,sinsow,albkof,ifix(para(12)))
*-->      berechnung der solarkonstanten in abhaengigkeit vom tag
          sk  = para(4)
          tag = real(jtag)
          call solkon(tag,sk,solark)
cc          goto (10,20,30) ifix(para(13))
cc  10      continue

*-->      berechnung der einstrahlung nach dobson & smith (1985 - 1989)
          einstr = solark * sinsow * (akoeff + bkoeff * sinsow) *
c          eaa=ea/100.
c          einstr=solark*(sinsow**2)/
c     &        (1.2*sinsow+(1.+sinsow)*eaa/1000.+0.046)*
c     &        (1-0.6*cn1)*
     &    (1. - albkof)

          tstd = (time / itfix) - (ifix(time/itfix))
          if(akoeff.eq.rnil.or.bkoeff.eq.rnil) einstr = rnil
*-->      ausdruck der ergebnisse
          ipjahr = ifix(xjahr) + 1900
c      if(jtag.eq.1.or.mod(jtag,50).eq.0) then
c         write(6,1100) ipjahr,jtag,tstd,
c    &       xtag,taglg,bewoel,akoeff,bkoeff,albkof,
c    &       bomgrd(asin(sinsow)),einstr,solark
c      end if
c          goto 100
c  20      continue
c          goto 100
c  30      continue

 100  continue
      istop = istop + 1

1000  format(1h1,
     *1x,
     *'ergebnisse der simulation:'/
     *2x,
     *26('-')/
     */
     &2x,
     *'breite = ',
     *f7.2/
     */
     */
     */
     *2x,
     *120('*')/
     *2x,
     *'* jahr  jtag tstd     xtag',
     &'    taglg(std)  bewoel     a       b      alb    sinsow',
     &'    ber. einstr.  solarkonst. *'/
     *2x,
     *120('*'))

1100  format(3x,i5,i5,f7.3,f10.2,f10.2,f9.1,1x,3f8.3,f8.2,f14.4,f15.4)


      return
      end
c-----------------------------------------------------------------------
      subroutine radice(in,ldim,time,bewoel,para,itage,einstr)
c-----------------------------------------------------------------------
c        global radition to sea ice
c        dobson und smith (1985,1988,1989)
c        journal of the royal meteorological society, 117, pp845-851.
c-----------------------------------------------------------------------
c parameter (in):
c     1) in              = maximale anzahl der jahre
c     2) time            = aktuelle modellzeit (in sec, std oder tagen)
c     3) ldim            = dimension des parameter vektors para(ldim)
c     4) bewoel          = aktueller bewoelkungsgrad (achtel)
c     5) para(ldim)      = parametervektor
c             para(1)  = breite in grad
c             para(2)  = jahr der simulation
c             para(3)  = umrechnung der zeit auf tage (z.b. 86400. fuer
c                        sekunden)
c             para(4)  = solarkonstante (w / m**2)
c             para(5)  = modellauswahl fuer subroutine param
c                      = 0 --> wird nicht beruecksichtigt
c                      = 1 --> dobson & smith (1985)
c                      = 2 --> dobson & smith (1988)
c                      = 3 --> dobson & smith (1989)
c             para(6)  = modellauswahl fuer subroutine tagums
c                      = 0 --> wird nicht beruecksichtigt
c                      = 1 --> dobson & smith (1985 - 1989)
c             para(7)  = modellauswahl fuer subroutine deklin
c                      = 0 --> wird nicht beruecksichtigt
c                      = 1 --> dobson & smith (1985 - 1989)
c             para(8)  = modellauswahl fuer subroutine taglng
c                      = 0 --> wird nicht beruecksichtigt
c                      = 1 --> dobson & smith (1985 - 1989)
c             para(9)  = modellauswahl fuer subroutine stdwnk
c                      = 0 --> wird nicht beruecksichtigt
c                      = 1 --> dobson & smith (1985 - 1989)
c            para(10)  = modellauswahl fuer subroutine sonhoe
c                      = 0 --> wird nicht beruecksichtigt
c                      = 1 --> dobson & smith (1985 - 1989)
c            para(11)  = modellauswahl fuer subroutine meteor
c                      = 0 --> wird nicht beruecksichtigt
c                      = 1 --> dobson & smith (1985 -1989)
c            para(12)  = modellauswahl fuer subroutine water
c                      = 0 --> wird nicht beruecksichtigt
c                      = 1 --> becker (1981)
c                      = 2 --> nach payne fuer breite = 45 grad
c            para(13)  = modellauswahl fuer strahlungsformel
c                      = 1 --> dobson & smith (1985 - 1989)
c parameter (out):
c     1) itage(in)       = vektor mit anzahl der tage jedes jahres
c     2) einstr          = simulationsergebnis
c-----------------------------------------------------------------------
c aufgerufene subroutinen:
c     1) param  --> belegung der bewoelkungsparametervektoren a() , b()
c     2) tagums --> setzt den anfangstag entsprechend den modell-
c                   erfordernissen
c     3) astro  --> berechnet die astronomischen groessen der modelle,
c                   ruft seinerseits auf:
c        4) deklin --> berechnet den deklinationswinkel der sonnenbahn
c        5) taglng --> berechnet die tageslaenge (sonnenunter-aufgang)
c        6) stdwnk --> berechnet den jeweiligen stundenwinkel
c        7) sonhoe --> berechnet den sinus des sonnenwinkels
c     8) meteor --> simuliert die meteorologischen einfluesse
c                   --> abhaengigkeit von der bewoelkung
c     9) water  --> simuliert die einfluesse der meeresoberflaeche
c                   --> albedo
c-----------------------------------------------------------------------

      parameter (jp=9)

      dimension itage(in)
      dimension a(jp) , b(jp) , para(ldim)

      common /atm/ rnil,istop

      ijahr = 1
      sinsow = 0.
      xtag = 0.
      akoeff = 0.
      bkoeff = 0.
      pi = atan(1.) * 4.

      if(istop.eq.1) then
c      write(6,1201)
 1201 format(////1x,'ausdruck der steuerparameter des strahlungsmodells:
     & ',//)
      do 23 i=1,ldim
c        write(6,1200) i,para(i)
  23  continue
      end if
 1200 format(1x,'para(',i2,') =',f15.3)

cval      if(para(1).lt.0.or.para(1).gt.90) stop 'breite out of range]'
      breite = grdbom(para(1))
      xjahr = para(2)
      itfix = ifix(para(3))

*-----belegung des vektors itage(in) mit anzahl der tage / jahr
cval      call setita(in,xjahr,itage)
      call setita(xjahr,itage(1))
c      if(istop.eq.1) print*,'itage (anzahl der tage im jahr)= ',itage(1)

*-----belegung der vektoren a() und b()
      call param(jp,a,b,ifix(para(5)))
c     if(istop.eq.1) write(6,1000) para(1)


*     berechnung der einstrahlung


*-------feststellung des aktuellen tages
        jtag = ifix(time / float(itfix)) + 1
*-->    umsetzen der tage
        call tagums(in,jtag,ijahr,itage,xtag,ifix(para(6)))
*-->    berechnung der tageslaenge und des sinus des sonnenwinkels
        call astro(time,itfix,ldim,breite,xtag,itage(ijahr),para,taglg,
     &             sinsow)
*---------feststelluing des stundenabschnittes eines tages
*-->      simulation der meteorologischen einfluesse
          call meteor(jp,a,b,bewoel,akoeff,bkoeff,ifix(para(11)))
*-->      simulation der einfluesse der meeresoberflaeche
          call water(jtag,sinsow,albkof,ifix(para(12)))
*-->      berechnung der solarkonstanten in abhaengigkeit vom tag
          sk  = para(4)
          tag = real(jtag)
          call solkon(tag,sk,solark)
          if (ifix(para(13)).ne.2.and.ifix(para(13)).ne.3) then
*-->      berechnung der einstrahlung nach dobson & smith (1985 - 1989)
          clo = bewoel/8.
          albkoi =  0.65
cc	  *clo**3.0
          einstr= solark * sinsow * (akoeff + bkoeff * sinsow) *
     &             (1. - albkoi)

          tstd = (time / itfix) - (ifix(time/itfix))
          if(akoeff.eq.rnil.or.bkoeff.eq.rnil) einstr = rnil
*-->      ausdruck der ergebnisse
          ipjahr = ifix(xjahr) + 1900
c      if(jtag.eq.1.or.mod(jtag,50).eq.0) then
c         write(6,1100) ipjahr,jtag,tstd,
c    &       xtag,taglg,bewoel,akoeff,bkoeff,albkof,
c    &       bomgrd(asin(sinsow)),einstr,solark
c      end if
          end if
      istop = istop + 1

1000  format(1h1,1x,'ergebnisse der simulation:'/2x,26('-')//
     &,2x,'breite = ',f7.2////2x,120('*')/2x,'* jahr  jtag tstd     xtag
     &    taglg(std)  bewoel     a       b      alb    sinsow'
     &,'    ber. einstr.  solarkonst. *'/2x,120('*'))
1100  format(3x,i5,i5,f7.3,f10.2,f10.2,f9.1,1x,3f8.3,f8.2,f14.4,f15.4)


      return
      end

c-----------------------------------------------------------------------
      subroutine astro(time,itfix,ldim,breite,xtag,jtage,para,taglg,
     * sinsow)
c-----------------------------------------------------------------------
c     berechnet die tageslaenge und den sinus des sonnenwinkels
c-----------------------------------------------------------------------
c parameter(in):
c       1) time         = aktuelle modellzeit
c       3) ldim         = dimension des parametervektors
c       4) breite       = geographische breite in bom
c       5) xtag         = nummer des tages
c       6) jtage        = anzahl der tage des betreffenden jahres
c       7) para(ldim)   = parametervektor (siehe auch subrt. radatm)
c                         liefert hier die nummern fuer die
c                         aufgerufenen subroutinen
c
c parameter(out):
c       1) taglg        = laenge des betreffenden tages
c       2) sinsow       = der berechnete sinus des sonnenwinkels
c-----------------------------------------------------------------------
c aufgerufene subroutines:
c       1) deklin       = berechnung des deklinationswinkels
c       2) taglng       = berechnung der tageslaenge
c       3) stdwnk       = berechnung des stundenwinkels
c       4) sonhoe       = berechnung des sinus des sonnenwinkels
c-----------------------------------------------------------------------

      dimension para(ldim)
      dekli = 0.
      sonauf = 0.
      sonunt = 0.
      stdwin = 0.
      iwrk = jtage
      call deklin(xtag,iwrk,dekli,ifix(para(7)))
      call taglng(breite,dekli,taglg,sonauf,sonunt,ifix(para(8)))
      call stdwnk(time,itfix,stdwin,ifix(para(9)))
      call sonhoe(time,itfix,sonauf,sonunt,stdwin,breite,dekli,sinsow,
     *            ifix(para(10)))


      return
      end

c-----------------------------------------------------------------------
      subroutine taglng(breite,dekli,taglg,sonauf,sonunt,nr)
c-----------------------------------------------------------------------
c     berechnet die tageslaenge , sonnenauf - und untergang
c-----------------------------------------------------------------------
c parameter(in):
c       1) breite        = geographische breite in bom
c       2) dekli         = deklinationswinkel
c       3) nr            = nummer, die das gewuenschte modell anwaehlt
c                        = 0 --> return
c                        = 1 --> dobson & smith (1985 - 1989)
c parameter(out):
c       1) taglng        = tageslaenge in std
c       2) sonauf        = zeit des sonnenaufgangs in std
c       3) sonunt        = zeit des sonnenunterganges in std
c-----------------------------------------------------------------------

      if(nr.eq.1) then

      stdwin = acos(-tan(breite)*tan(dekli))
      stdwin = bomgrd(stdwin)
      sonauf = 12.-(stdwin/15.)
      sonunt = 12.+(stdwin/15.)
      taglg  = sonunt-sonauf

      end if

      return

      end

c-----------------------------------------------------------------------
      subroutine sonhoe(time,itfix,sonauf,sonunt,stdwin,breite,dekli,
     * sinsow,nr)
c-----------------------------------------------------------------------
c     berechnet den sinus des sonnenwinkels
c-----------------------------------------------------------------------
c parameter(in):
c       1) time        = aktuelle modellzeit
c       2) itfix       = umrechnnung der zeit in tage (z.b. 86400 f sec
c       3) sonauf      = zeit des sonnenaufganges in std
c       4) sonunt      = zeit des sonnenunterganges in std
c       5) stdwin      = stundenwinkel des tages in grad
c       6) breite      = geographische breite in bom
c       7) dekli       = deklinationswinkel in bom
c       8) nr          = nummer, die das gewuenschte modell anwaehlt
c                      = 0 --> return
c                      = 1 --> dobson & smith (1985 - 1989)
c parameter(out):
c       1) sinsow      = berechneter sinus des sonnenwinkels
c-----------------------------------------------------------------------

      if(nr . eq . 0)  return
      zeit = ((time / itfix) - (ifix(time/itfix))) * 24.
      if (nr.eq.1) then
         stwin = grdbom(stdwin)
         sinsow = sin(breite)*sin(dekli)+cos(breite)*cos(dekli)*
     1    cos(stwin)
         if( zeit.le.sonauf. or .zeit.ge.sonunt )  sinsow = 0.
      end if

      return


      end
c-----------------------------------------------------------------------
      subroutine deklin(xtag,jtage,dekli,nr)
c-----------------------------------------------------------------------
c     berechnet den deklinationswinkel
c-----------------------------------------------------------------------
c parameter(in):
c       1) xtag          = nummer des betreffenden tages
c       2) jtage         = anzahl der tage des betreffenden jahres
c       3) nr            = nummer, die das gewuenschte modell anwaehlt
c                        = 0 --> return
c                        = 1 --> dobson & smith (1985 - 1989)
c parameter(out):
c       1) dekli         = deklinationswinkel
c-----------------------------------------------------------------------
      if (nr.eq.1) then
         pi = atan(1.) * 4.
         tage = float(jtage)
         dekli = asin(sin(23.5*pi/180.)*sin(xtag*2.*pi/tage))
      end if


      return
      end

c-----------------------------------------------------------------------
      subroutine solkon(tag,sk,solark)
c-----------------------------------------------------------------------
c     berechnet die solarkonstante und deren variation durch die
c     elliptische erdbahn. im perihel groesste sonnennaehe (3. jan.)
c     und im aphel groesste ferne (4. juli).
c     variation vom mittleren abstand : 3.5 %.
c-----------------------------------------------------------------------
c parameter (in):
c        1) tag         = aktueller tag des jahres
c        2) sk          = solarkonstante, (=para(4) aus hp)
c parameter (out):
c        3) solark      = solarkonstante
c-----------------------------------------------------------------------

      dq     = 149.60
      a      = 149.499
      eps    = 0.01675
      t0     = 3.0
      t      = 365.0
      pi     = atan(1.) * 4.
      d      = a*(1.0-eps*cos(2.0*pi/t*(tag-t0)))
      solark = sk*(dq/d)**2


      return
      end

c-----------------------------------------------------------------------
      subroutine stdwnk(time,itfix,stdwin,nr)
c-----------------------------------------------------------------------
c     berechnet die stundenwinkel des sonnenstandes
c-----------------------------------------------------------------------
c parameter(in):
c       1) time         = aktuelle modellzeit
c       2) itfix        = umrechnung der zeit in tage (z.b. 86400 f sec
c       3) nr           = nummer, die das gewuenschte modell anwaehlt
c                       = 0 --> return
c                       = 1 --> dobson & smith (1985 - 1989)
c parameter(out):
c       1) stdwnk       = vektor mit den berechneten stundenwinkeln
c-----------------------------------------------------------------------

      if(nr.eq.1) then
         timsec = ((time/itfix) - (ifix(time/itfix)))
         stdwin = (timsec * 360.) - 180.
      end if


      return
      end

c-----------------------------------------------------------------------
      subroutine setita(ajahr,itage)
cval      subroutine setita(in,ajahr,itage)
c-----------------------------------------------------------------------
c     besetzt itage(in) mit anzahl tagen in 'in' jahren
c     ajahr = anfangsjahr
c-----------------------------------------------------------------------

      ia = ifix(ajahr)-1
      if(ia.gt.1900) ia = ia-1900
      ia = ia+1
      itage = 365
      if((ia/4)*4.eq.ia) itage = 366


      return
      end

c-----------------------------------------------------------------------
      subroutine param(jp,a,b,nr)
c-----------------------------------------------------------------------
c     belegt die vektoren a(), b() mit bewoelkungsparametern
c-----------------------------------------------------------------------
c parameter(in):
c        1) jp  = dimension fuer a() , b()
c        2) nr  = nummer, die das gewuenschte modell anwaehlt
c               = 0 --> return
c               = 1 --> dobson & smith (1985,ows p coeff., okta)
c               = 2 --> dobson & smith (1988,ows p coeff., fog)
c               = 3 --> dobson & smith (1989,si coeff., fog)
c               = 4 --> dobson & smith (1989,sc coeff., fog)
c parameter(out):
c        1) a(jp)       = vektor, der bewoelkungsparameter aufnimmt
c        1) b(jp)       = vektor, der bewoelkungsparameter aufnimmt
c-----------------------------------------------------------------------

      dimension a(jp) , b(jp)

      if(nr.eq.1) then
c-----dobson + smith (1985)
         a(1) = 0.406
         a(2) = 0.533
         a(3) = 0.479
         a(4) = 0.426
         a(5) = 0.385
         a(6) = 0.352
         a(7) = 0.310
         a(8) = 0.235
         a(9) = 0.117
         b(1) = 0.391
         b(2) = 0.308
         b(3) = 0.389
         b(4) = 0.419
         b(5) = 0.473
         b(6) = 0.474
         b(7) = 0.439
         b(8) = 0.388
         b(9) = 0.304
      else if (nr.eq.2) then
c-----dobson + smith (1988)
         a(1) = 0.400
         a(2) = 0.517
         a(3) = 0.474
         a(4) = 0.421
         a(5) = 0.380
         a(6) = 0.350
         a(7) = 0.304
         a(8) = 0.230
         a(9) = 0.106
         b(1) = 0.386
         b(2) = 0.317
         b(3) = 0.381
         b(4) = 0.413
         b(5) = 0.468
         b(6) = 0.457
         b(7) = 0.438
         b(8) = 0.384
         b(9) = 0.285
      else if (nr.eq.3) then
c-----dobson + smith (1989, sable island-coeff., fog)
         a(1) = 0.539
         a(2) = 0.525
         a(3) = 0.497
         a(4) = 0.479
         a(5) = 0.437
         a(6) = 0.382
         a(7) = 0.330
         a(8) = 0.226
         a(9) = 0.112
         b(1) = 0.303
         b(2) = 0.323
         b(3) = 0.356
         b(4) = 0.363
         b(5) = 0.377
         b(6) = 0.410
         b(7) = 0.422
         b(8) = 0.427
         b(9) = 0.343
      else if (nr.eq.4) then
c-----dobson + smith (1989, swift current-coeff., fog)
         a(1) = 0.634
         a(2) = 0.613
         a(3) = 0.589
         a(4) = 0.564
         a(5) = 0.539
         a(6) = 0.500
         a(7) = 0.452
         a(8) = 0.379
         a(9) = 0.296
         b(1) = 0.201
         b(2) = 0.223
         b(3) = 0.234
         b(4) = 0.236
         b(5) = 0.228
         b(6) = 0.231
         b(7) = 0.219
         b(8) = 0.176
         b(9) = 0.048
      end if


      return
      end

c-----------------------------------------------------------------------
      subroutine tagums(in,jtag,ijahr,itage,xtag,nr)
c-----------------------------------------------------------------------
c     setzt den anfangstag fuer die simulation entsprechend den
c     modellerfordernissen
c-----------------------------------------------------------------------
c parameter(in):
c       1) in         = anzahl der simulationsjahre
c       2) jtag       = nummer des betreffenden tages
c       3) ijahr      = nummer des betreffenden jahres
c       4) itage(in)  = vektor mit anzahl der tage jedes jahres
c       5) nr         = nummer, die das gewuenschte modell anwaehlt
c                     = 0 --> return
c                     = 1 --> dobson & smith (1985 - 1989)
c parameter(out):
c       1) xtag       = nummer des tages, die fuer die berechnung
c                       verwendet wird
c-----------------------------------------------------------------------

      dimension itage(in)


      if (nr.eq.1) then
      tagja = float(itage(ijahr))
      korekt = -81
      if(itage(ijahr).eq.365) korekt = -80
      xtag = tagja+korekt+jtag
      if(xtag.gt.tagja) xtag = xtag-tagja
      end if


      return
      end
c-----------------------------------------------------------------------
      subroutine meteor(jp,a,b,bewoel,akoeff,bkoeff,nr)
c-----------------------------------------------------------------------
c     simuliert die meteorologischen einfluesse auf die  globalstrahlung
c-----------------------------------------------------------------------
c parameter(in):                                                       *
c       1) jp          = dimension --> a(jp) , b(jp)                   *
c       2) a(jp)       = vektor mit bedeckungsparametern               *
c       3) b(jp)       = vektor mit bedeckungsparametern               *
c       4) bewoel      = bewoelkungsgrad (achtel)                      *
c       5) nr          = nummer, die das gewuenschte modell anwaehlt   *
c                      = 0 --> return                                  *
c                      = 1 --> dobson & smith (1985 - 1989)            *
c                                                                      *
c parameter(out):                                                      *
c       1) akoeff      = bedeckungsparameter, der benutzt wird         *
c       2) bkoeff      = bedeckungsparameter, der benutzt wird         *
c-----------------------------------------------------------------------

      dimension a(jp) , b(jp)
      common /atm/ rnil,istop

      if (nr.eq.1) then
         if(bewoel.ge.0.and.bewoel.le.8) then
            index = nint(bewoel)
            akoeff = a(index+1)
            bkoeff = b(index+1)
         else
            akoeff = rnil
            bkoeff = rnil
         endif
      endif


      return
      end
c-----------------------------------------------------------------------
      function grdbom(x)
c-----------------------------------------------------------------------
c     wandelt altgrad in bogenmass
c-----------------------------------------------------------------------

      pi = atan(1.) * 4.
      pi180 = pi / 180.
      grdbom = x * pi180


      return
      end
c-----------------------------------------------------------------------
      subroutine water(jtag,sinsow,albkof,nr)
c-----------------------------------------------------------------------
c     simuliert die einfluesse der meeresoberflaeche auf die        *
c     globale einstrahlung                                          *
c-----------------------------------------------------------------------
c parameter(in):                                                       *
c       1) jtag          = nummer des tages                            *
c       2) sinsow        = sinus des sonnenwinkels des aktuellen tages *
c       3) nr            = nummer, die das gewuenschte modell anwaehlt *
c                        = 0 --> return , albkof = rnil                *
c                        = 1 --> nach becker (1981)                    *
c                        = 2 --> nach payne fuer breite = 54 grad      *
c parameter(out):                                                      *
c       1) albkof        = albedo parameter                            *
c-----------------------------------------------------------------------

      parameter(kdim = 12)

      dimension albedo(kdim)
      common /atm/ rnil,istop

      albkof = rnil
      if (nr.eq.1) then
c albedo nach becker (1981)
         ssin=max(sinsow,0.)
         albkof = 0.04 / (0.04 + ssin)
      else if (nr.eq.2) then
c albedo nach payne fuer breite = 60 grad
         albedo(1) = 0.28
         albedo(2) = 0.12
         albedo(3) = 0.09
         albedo(4) = 0.07
         albedo(5) = 0.07
         albedo(6) = 0.07
         albedo(7) = 0.06
         albedo(8) = 0.07
         albedo(9) = 0.07
         albedo(10) = 0.10
         albedo(11) = 0.16
         albedo(12) = 0.44
         albkof = albedo(monat(jtag))
      end if


      return
      end
c-----------------------------------------------------------------------
      function bomgrd(x)
c-----------------------------------------------------------------------
c     wandelt bogenmass in altgrad
c-----------------------------------------------------------------------
      pi = atan(1.) * 4.
      r180pi = 180. / pi
      bomgrd = x * r180pi


      return
      end
c-----------------------------------------------------------------------
      function monat(i)
c-----------------------------------------------------------------------
c     berechnet, in welchem monat der tag i liegt
c-----------------------------------------------------------------------
      dimension mons(12)
      mons(1)=31
      mons(2)=28
      mons(3)=31
      mons(4)=30
      mons(5)=31
      mons(6)=30
      mons(7)=31
      mons(8)=31
      mons(9)=30
      mons(10)=31
      mons(11)=30
      mons(12)=31
      j = i
      monat = 0
!CDIR UNROLL=12
      do m=1,12
        if(j.le.mons(m).and.monat.eq.0) monat = m
        j = j-mons(m)
      enddo
      return
      end
c-------------------------------------------------------------------------
      subroutine schubsp(taubot1)

      include 'C_model.f'
      include 'C_index.f'

      parameter(ilop1=ilo+1)
      parameter(khor1=khor+1)

      common /ind/ iwet(khor1),ldep(khor),lazc(khor),
     *   indend(n),isornr(n),isorsr(n),islab(n)
      common uc(ndrei),vc(ndrei),stc(ndrei),avc(ndrei),z(m,n)
      common /taubot/ taubu(khor), taubv(khor)
      common /met/ windx,windy,wlam,pnull,rhoq(ilo),refrho(ilo),stress
      dimension rho2d(khor),taubot1(khor)
#ifdef MPI
       include 'C_mpi.f'
       common  /cord/ I1,I2,J1,J2,IA1,IA2,JA1,JA2,IB1,IB2,JB1,JB2
     1,ID1,ID2,JD1,JD2
       common  /lbkhor/ khorl,lzet(khor),lb0(n),le0(n),
     * lb1(n),le1(n),lb2(n),le2(n)
#endif

        call estate
c        tumax=maxval(taubu)
c        tvmax=maxval(taubv)
c        print*,'tu,tv ',tumax,tvmax
#ifdef MPI
      do j=J1,J2
         do lw=lb0(j),le0(j)
#else  
      do j=1,n
      do lw=lb(j),le(j)
#endif
            nwet = indwet(lw)
            ldown = lazc(lw)
            rho2d(lw)=(stc(nwet+ldown)+rhoq(ldown))*10.**3.
        end do
        end do

#ifdef MPI
        do I=1,khorl
        lw=lzet(I)
#else
        do lw=1,khor
#endif
        taubu1=taubu(lw)*rho2d(lw)
        taubv1=taubv(lw)*rho2d(lw)
        taubot1(lw)=sqrt(taubu1**2.+taubv1**2.)
        end do
        
        return
        end

      Subroutine Rzero2d(a,m,n)
      dimension a(m,n)
      do i=1,m
        do j=1,n
      a(i,j)=0.
        enddo
      enddo
      return
      end

c-------------------------------------------------------------------
c      subroutine read_riv_bio_input(nfuriv,riv_bio,nbio,nbiox,numeintr)
c-------------------------------------------------------------------
      subroutine read_riv_input(nfuriv,rivinp,numeintr)
      include 'C_mpi.f'
       dimension rivinp(numeintr)
      integer nfuriv 
       real iijahr,iimonth,iiday
      
       inumeintr=numeintr+3  
cKB      read(nfuriv,'(<inumeintr>e16.7)')iijahr,iimonth,iiday,
c      read(nfuriv,'(222e16.7)')iijahr,iimonth,iiday,
      read(nfuriv,'(225e16.7)')iijahr,iimonth,iiday,
     &(rivinp(num),num=1,numeintr)
       if(myid.eq.0)then
       print*,'read riv: ', iijahr,iimonth,iiday,
     &rivinp(1),rivinp(numeintr),numeintr
       endif
       return
       end 

c---------------------------------------------------------------------------------------
#ifdef _FABM_
      subroutine change_units(f,n_choice,inbio)
      include 'C_model.f'
      include 'C_nbio_fabm.f'
      
      parameter (ips=3,ipl=4,ibg=15,izs=5,izl=6,idet=7,inh4=8,
     &idom=9,ino3=10,ipo4=11,isio=12,io2=13,iop=14,ifi1=16,
     &ised1=1,ised2=2,ised3=3,imb1=4)
      dimension f(ndrei,*)
      common /bio_out/bio_out(ndrei,3:nbio+nsed) 
      real REDF(10)
c/////redfield ratios/////
       DATA
     &      REDF(1)/6.625/,  !C_N    PrmREDF(1) mol_C/mol_N=    6.625
     &      REDF(2)/106.0/,  !C_P    PrmREDF(2) mol_C/mol_P=  106.000
     &      REDF(3)/6.625/,  !C_SiO  PrmREDF(3) mol_C/mol_Si=    6.625
     &      REDF(4)/16.0/,  !N_P    PrmREDF(4) mol_N/mol_P=   16.000
     &      REDF(5)/1.0/,   !N_SiO  PrmREDF(5) mol_N/mol_Si=    1.000
     &      REDF(6)/12.010/,   !C_Cmg  PrmREDF(6) mg_C/mmol_C=   12.010
     &      REDF(7)/44.6608009/,   !mm_ml  PrmREDF(7)(mmolO2/m**3)/(mlO2/l)   = 44.6608009
     &      REDF(8)/14.007/,   !  mmolN/mgN   
     &      REDF(9)/30.97/,   ! mmolP/mgP
     &      REDF(10)/28.09/   ! mmolSi/mgSi
       if(n_choice.eq.1)then

*------change units
         do l=1,ndrei
         do ibio=3,7  !Ps, Pl, Zs, Zl, De
         bio_out(l,ibio)=f(l,ibio)        !*REDF(6)  !input in   mgC /m**3 to [mmolC/m**3]
         enddo
         bio_out(l,idom)=bio_out(l,idom)  !dom ==det
c         bio_out(l,idom)=f(l,idom)  !dom ==det
         bio_out(l,inh4)=f(l,inh4)*REDF(1)*REDF(6)  !input in mmolN /m**3 to [mmolC/m**3]
         bio_out(l,ino3)=f(l,ino3)*REDF(1)*REDF(6)  !input in mmolN /m**3 to [mmolC/m**3]
         bio_out(l,ipo4)  =f(l,ipo4)  *REDF(2)*REDF(6)  !input in mmolP /m**3 to [mmolC/m**3]
         bio_out(l,isio)  =f(l,isio)  *REDF(3)*REDF(6)  !input in mmolSi/m**3 to [mmolC/m**3]
         bio_out(l,iop)  =f(l,iop)  *REDF(3)*REDF(6)  !input in mmolSi/m**3 to [mmolC/m**3]

         bio_out(l,io2)  =f(l,io2)  *REDF(7)          !input in milliliters O2/liter to  mmolO2/m**3

        do ibio=15,inbio
         bio_out(l,ibio)  =f(l,ibio)                    !bg
        enddo
        enddo !l=1,nrei

      else
*------change units
         do l=1,ndrei
         do ibio=3,7  !Ps, Pl, Zs, Zl, De
         bio_out(l,ibio)=f(l,ibio)         !*REDF(16)  !output from [mmolC/m**3] to mgC /m**3
         enddo
         bio_out(l,idom)=f(l,idom)
         bio_out(l,inh4)=f(l,inh4)*1./REDF(1)*1./REDF(6)  !output from [mmolC/m**3] to mmolN /m**3
         bio_out(l,ino3)=f(l,ino3)*1./REDF(1)*1./REDF(6)  !output from [mmolC/m**3] to mmolN /m**3
         bio_out(l,ipo4)  =f(l,ipo4)  *1./REDF(2)*1./REDF(6)  !output from [mmolC/m**3] to mmolP /m**3
         bio_out(l,isio)  =f(l,isio)  *1./REDF(3)*1./REDF(6)  !output from [mmolC/m**3] to mmolSi/m**3
         bio_out(l,iop)  =f(l,iop)  *1./REDF(3)*1./REDF(6)  !output from [mmolC/m**3] to mmolSi/m**3

         bio_out(l,io2)  =f(l,io2)  *1./REDF(7)           !output from mmolO2/m**3 to mlO2/l
       do ibio=15,inbio
         bio_out(l,ibio)  =f(l,ibio)                    !bg
       enddo
       enddo !l=1,nrei
       endif
       end 
#endif
