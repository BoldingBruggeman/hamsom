c----------------------------------------------------------------------- 
      subroutine fluxes(ind2d,tl,sphumi,bew,pac,uw,vw,
     &   uc,vc,tec,qw ,id3sur,wspeed,q10,tl10,
     &    fqs,fql,cddu,evap) 

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	
ccc         Calculation of sea surface fluxes depending on stability of the 
ccc         atmospheric boundary layer following Timo Vihma (1995): 
ccc         Atmosphere-surface interactions over polar oceans and 
ccc         heterogenous surfaces, Finnish Marine Research No. 264
ccc                        C. Schrum        March 1997/Feb 2002
ccc                        I. Alekseeva
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      include 'C_mpi.f' 
      include 'C_model.f'
      dimension ind2d(m,n),id3sur(khor) 
      dimension  tl(khor),taup(khor),bew(khor) 
      dimension  uw(khor),vw(khor),pac(khor) 
   
      dimension  qw(khor),wspeed(khor) 
      dimension tec(ndrei),uc(ndrei),vc(ndrei) 
      dimension q10(khor),tl10(khor) 
 
      dimension fqs(khor),fql(khor),cddu(khor) 
  
      real phihq(khor),phim(khor)  
      real  cdd(khor),cdeh(khor),rhol(khor) 
      real  zro(khor),zrtq(khor),evap(khor) 
      real uscale(khor),tscale(khor),qscale(khor) 
      real psim(khor),psihq(khor),diff(khor),tau(khor) 
      real sspeed(khor),hsens(khor),rl(khor) 
      real sphumi(khor),spsur(khor)
      real eps 
      integer iterat(khor) ,ab_tscale(khor)
      logical loxstat,loxstat1
#ifdef MPI
      common /lbkhor/ khorl,lzet(khor),lb0(n),le0(n),
     * lb1(n),le1(n),lb2(n),le2(n)
#endif
      loxstat=.true.
      loxstat1=.false.

c
c---------------------------------------------------------------------------------
c   set level of wind measurements zw, level of humidity measurements zh,
c   level of air-temperature measurements zt
c---------------------------------------------------------------------------------
c
      zw   =   10. 
      zth  =    2. 
      xk   =    0.4 
      cp   = 1005. 
      xlat = 2500.*1000. 
      g    =    9.81 
      bound=    0.1 
      eps  =    1.e-7 
c--------------------------------------------------------------------------------
c   calculation of the density  of air, specific humidity, wind speed
c   and set phim, phih and phiq eq. 1 for neutral stability
c--------------------------------------------------------------------------------
c
#ifdef MPI
       do ii=1,khorl
       i=lzet(ii)
#else
      do i=1,khor
#endif 
        diff(i) = 1000. 
 
c          ee = exp(-6763.6/(taup(i)+273.15)- 
c     *     4.9283*log(taup(i)+273.15)+54.23) 
c          if(taup(i).lt. 0.)then 
c            ee = exp(-6141./(taup(i)+273.15)+24.3) 
c          endif 
 

       if(tec(id3sur(i)).ge.0.)then 
       esur = exp(-6763.6/(tec(id3sur(i))+273.15)- 
     &  4.9283*log(tec(id3sur(i))+273.15)+54.23) 
       else        
       esur = exp(-6141./(tec(id3sur(i))+273.15)+24.3) 
       endif 
     
c        sphumi(i) = 0.622*ee/(10.*pac(i)-0.378*ee) 
 
        spsur(i) = 0.622*esur /(10.*pac(i)-0.378*esur) 
        rhol(i) = 10.*pac(i)/(2.8705*(tl(i)+273.15)* 
     *   (1.+0.608*sphumi(i))) 
 
	wspeed(i) = sqrt(uw(i)**2+vw(i)**2) 
	sspeed(i) = 0. 
	psim(i)   = 1. 
	psihq(i)  = 1. 
	phim(i)   = 0. 
	iterat(i) = 1 
	uscale(i) = 1.
	phihq(i)=0.
      enddo 
c
c--------------------------------------------------------------------------------
c     set roughness lengths zro, zrtq  from V(10m) and first guess for u*
c     t* and q*
c--------------------------------------------------------------------------------
c
#ifdef MPI
       do ii=1,khorl
       i=lzet(ii)
#else
      do i=1,khor 
#endif
        cdbb = (0.8+0.065*wspeed(i))*1.e-3 
        cdee=(0.824+0.041*wspeed(i))*1.e-3
 
        cdd(i)  = cdbb 
        cdeh(i) = cdee 
        zro(i)  = exp(log(10.)-xk/sqrt(cdbb)) 
        zrtq(i) = exp(log(2.)-xk*sqrt(cdbb)/cdee) 
        uscale(i) = cdd(i)*(wspeed(i))**2 
        uscale(i) = max(sqrt(uscale(i)),eps) 
        tscale(i) = (tl(i)-tec(id3sur(i)) ) 
        tscale(i) = cdeh(i)*tscale(i)*wspeed(i) 
        tscale(i) = tscale(i)/(xk*uscale(i)) 
        qscale(i) = (sphumi(i) -spsur(i)) 
        qscale(i) = cdeh(i)*qscale(i)*wspeed(i) 
        qscale(i) = qscale(i)/(xk*uscale(i)) 

        ab_tscale(i)=abs( tscale(i) )
      enddo 

#ifndef MPI
      if(loxstat1)then
	  write(*,*)'1 before iteration'
        call xstat(tscale,khor,zmax,zmin,zmean,zsum)
	  print*,'tscale  min,max,mean,sum ',zmin,zmax,zmean
        call xstat(ab_tscale,khor,zmax,zmin,zmean,zsum)
	  print*,'ab_tscale min,max,mean,sum ',zmin,zmax,zmean

        write(*,*)'sphumi:'
        call xstat(sphumi,khor,zmax,zmin,zmean,zsum)
	  print*,'sphumi min,max,mean,sum ',zmin,zmax,zmean

        write(*,*)'spsur:'
        call xstat(spsur,khor,zmax,zmin,zmean,zsum)
	  print*,'spsur min,max,mean,sum ',zmin,zmax,zmean

        call xstat(evap,khor,zmax,zmin,zmean,zsum)
	  print*,'evap min,max,mean,sum ',zmin,zmax,zmean
      endif
#endif
c------------------------------------------------------------------------------
c     start of the iteration
c-----------------------------------------------------------------------------
c
      iiicount=0
      imax=0 
#ifdef MPI
       do ik=1,khorl
       i=lzet(ik)
#else
      do i=1,khor   
#endif
       imax=max0(iiicount,imax) 
       iiicount=0 
       to=(tl(i)+tec(id3sur(i)))/2.+273.15
       do ii=1,40 
      do iit=1,iterat(i) 
        iiicount=iiicount+1 
c------------------------------------------------------------------------------
c     calculate new 10 m values for air temperature and humidity
c-----------------------------------------------------------------------------
c
        tl10(i) = tl(i)+tscale(i)*log(5.)*psihq(i) 
        q10(i)  = sphumi(i)+qscale(i)*log(5.)*psihq(i) 
c
c-------------------------------------------------------------------------------
c    calculate transfer coefficients for the input-data height
c--------------------------------------------------------------------------------
c       cdeh(i) = xk/( log(2./zrtq(i))+phihq(i) ) 
c        cdd(i)  = xk/(log(2./zro(i))+phim(i))  


           xlog =  log(2./zrtq(i)) +phihq(i) 
            xlog =sign (max(abs(xlog),eps),xlog)
        cdeh(i) = xk/xlog


          xlog =  log(2./zro(i))+phim(i) 
          xlog =sign (max(abs(xlog),eps),xlog)
        cdd(i)  = xk/xlog

        cdd(i)  = cdd(i)*cdd(i) 
        cdeh(i) = cdeh(i)*sqrt(cdd(i)) 
        cdee    = cdeh(i) 
 
c
c------------------------------------------------------------------------------
c    calculate fluxes and 2-m wind speed w2
c-----------------------------------------------------------------------------
c


        w2 = wspeed(i)-uscale(i)*log(5.)*psim(i)/xk 

	  hsens(i) = cdeh(i)*rhol(i)*w2*(-tl(i)+tec(id3sur(i))) 
        hsens(i) = hsens(i)*cp
	  hsens(i)=sign( max(eps,abs(hsens(i))) ,hsens(i) )
	  
        evap(i)  = cdeh(i)*rhol(i)*w2*(-sphumi(i)+spsur(i))
c
c---------------------------------------------------------------------------------
c    calculate z/L  by using the fluxes
c--------------------------------------------------------------------------------
c

          rl(i) = 1.+0.61*to*cp*evap(i)/(hsens(i)) 
          rl(i) = g*xk*hsens(i)*rl(i) 

          u3=max(uscale(i)**3,eps) 
 
          rl(i) = rl(i)/(u3*to*rhol(i)*cp) 


          istab = 1 
c---------------------------------------------------------------------------------
c    calculate the universal functions phim/hq  and psim/hq 
c    for neutral stratification
c--------------------------------------------------------------------------------
        if( abs(rl(i)).le.eps)  then 
          istab    = 0 
          psihq(i) = 1. 
          phihq(i) = 0.  
          psim(i)  = 1. 
          phim(i)  = 0. 
        endif  
 

        do ixx=1,istab 
          istaby = 0 
	    istabn = 0 
          if(rl(i).lt.-eps) istaby = 1 
          if(rl(i).gt.eps) istabn = 1 

          
c
c---------------------------------------------------------------------------------
c    calculate the universal functions phim/hq  and psim/hq 
c    for stable stratification
c--------------------------------------------------------------------------------
c 
 	  
	   do ix=1,istaby 
ccc                 print*,'Im stabilen Bereich'
		zrl=-2.*rl(i)
		 phihq(i)=-0.7*zrl-0.75*(zrl-5./0.35)
     +		 *exp(-0.35*zrl)-0.75*5./0.35
                 zrl= -zro(i)*rl(i)
		 zrl=-0.7*zrl-0.75*(zrl-5./0.35)
     +		 *exp(-0.35*zrl)-0.75*5./0.35

		 zzrl=-10.*rl(i)
		 phim(i)=-0.7*zzrl-0.75*(zzrl-5./0.35)
     +		 *exp(-0.35*zzrl)-0.75*5./0.35
		 phim(i)=(phim(i)-zrl)
		 zrl=-zrtq(i)*rl(i)
		 zrl=-0.7*zrl-0.75*(zrl-5./0.35)
     +		 *exp(-0.35*zrl)-0.75*5./0.35

		 phihq(i)=(phihq(i)-zrl)
          
		 psihq(i)=rl(i)*5.*2.
		 psim(i)=psihq(i)
ccccc         print*,'hier',psim(i),psihq(i),rl(i) 
             enddo 
c
c---------------------------------------------------------------------------------
c    calculate the universal functions phim/hq  and psim/hq 
c    for unstable stratification
c--------------------------------------------------------------------------------
c 
	     do ix=1,istabn		  
c		psim(i) = sqrt(  1.+19.3*2.*max(0.00001,rl(i)) ) 
		psim(i) = sqrt(  1.+19.3*2.*rl(i) ) 
	    psim(i) = 1./sqrt(psim(i)) 
	    zz = sqrt((1.+19.3*zro(i)*rl(i))) 
	    zz = 1./sqrt(zz) 
	    phim(i) = 2.*log(0.5+1./(2.*psim(i))) 
     * 	     +log(0.5+1./(2.*psim(i)**2)) 
     *       -2.*atan(1./psim(i))+3.14/2. 
            phim(i) = phim(i)-2.*log(0.5+1./(2.*zz)) 
     * 	     -log(0.5+1./(2.*zz**2)) 
     *       +2.*atan(1./zz)-3.14/2. 
	    zophihq  = (1.+12.*zrtq(i)*rl(i)) 
	    zophihq  = 1./sqrt(zophihq) 
     
            psihq(i) = (1.+24.*rl(i)) 
            psihq(i) = 1./sqrt(psihq(i)) 
  		phihq(i) = 2.*log(0.5+1./(2.*psihq(i))) 
	    phihq(i) = phihq(i)-2.*log(0.5+1./(2.*zophihq)) 
          enddo 

        enddo !        do ixx=1,istab 


c       
c------------------------------------------------------------------------------
c     calculate new scaling parameters
c-----------------------------------------------------------------------------
c
        tscale(i) = -hsens(i)/(rhol(i)*cp*xk*uscale(i)) 
        qscale(i) = -evap(i)/(rhol(i)*xk*uscale(i)) 
c
c---------------------------------------------------------------------------------
c    continue with estimations of the scaling parameters if ......
c-------------------------------------------------------------------------------
c
        if(abs(diff(i)-10.*rl(i)).lt.bound) iterat(i) = 0 
        diff(i) = 10.*rl(i) 
 
      enddo 
      enddo 
      enddo 

c------------------------------------------------------------------------------
c     End of the iteration
c-----------------------------------------------------------------------------

c-----------------------------------------------------------------------------
c in case of no convergence [iterat(i).gt.0], apply neutral stability parameterization 
c-----------------------------------------------------------------------------
       sum=0 
#ifdef MPI
       do i=1,khorl
       ii=lzet(i)
#else
       do ii=1,khor 
#endif
       sum=sum+iterat(ii) 
       if (iterat(ii).gt.0) then
        cdee=(0.824+0.041*wspeed(ii))*1.e-3
        hsens(ii) = cdee*rhol(ii)*wspeed(ii)*(-tl(ii)+tec(id3sur(ii)))
        hsens(ii) = hsens(ii)*cp
        hsens(ii)=sign( max(eps,abs(hsens(ii))),hsens(ii))
        evap(ii)  = cdee*rhol(ii)*wspeed(ii)*(-sphumi(ii)+spsur(ii))
        tl10(ii) = tl(ii) 
        q10(ii)  = sphumi(ii) 
      end if 
      enddo  
#ifndef MPI    
      if(loxstat1)then
        write(*,*)'3 no convergence '
        write(*,*)'evap:'
        call xstat(evap,khor,zmax,zmin,zmean,zsum)
        print*,'evap min,max,mean ',zmin,zmax,zmean
        write(*,*)'hsens:'
        call xstat(hsens,khor,zmax,zmin,zmean,zsum)
        print*,'hsense min,max,mean ',zmin,zmax,zmean
      endif
#endif
c----------- windstress ---------------------------------
#ifdef MPI
      do ii=1,khorl
       i=lzet(ii)
#else
      do i=1,khor 
#endif 
c================================================================== 
c       Drag coefficient Cd (cdd) for neutral stability in the atmosphere 
c================================================================== 
*1-----------------------------
c       Charnock (linear form) (1955) 
c        cdd(i) = (0.7+0.09*wspeed(i))*1.e-3 
*2----------------------------- 
c       Smith and Banke (1975) 
c        cdd(i) = (0.63+0.066*wspeed(i))*1.e-3 !cd8
c        cdd(i) = (0.63+0.05*wspeed(i))*1.e-3 !cd9
*3----------------------------- 
c       Wu (1980)  
c       cdd(i) = (0.8+0.065*wspeed(i))*1.e-3 
*3----------------------------- 
c       Large and Pond (1981)  
        if(wspeed(i)>=10.)cdd(i) = (0.49+0.065*wspeed(i))*1.e-3 
        if(wspeed(i)<10.)cdd(i) = 1.14*1.e-3 


*4-----------------------------
c      cdd mit Stabilitaetsansatz f�r die Atmosph�re 
c      ansatz nach J. Launianien and T. Vihma (1990) 


c         if(iterat(i).gt.0.and.wspeed(i).gt.22.)then   !Irina
c         if(iterat(i).gt.0.or.wspeed(i).gt.18.)then   !Irina
c        if(iterat(i).gt.0)then                        !Irina
c------in case of no convergence, apply neutral stability parameterization, Wu (1980)
c         cdd(i) = (0.8+0.065*wspeed(i))*1.e-3   
c	 else
c          xlog=log(10./zro(i))-phim(i)
c          xlog =sign(max(abs(xlog),eps),xlog)
c       cdd(i) = (xk/(log(10./zro(i))-phim(i)))**2.0 
c       cdd(i) = (xk/xlog)**2
c       endif
 
c================================================================== 
 
 
        cddu(i) = cdd(i) 


*-------- wind stress-------------------
        uw(i)  = cdd(i)*rhol(i)*(wspeed(i))*(uw(i)) 
        vw(i)  = cdd(i)*rhol(i)*(wspeed(i))*(vw(i)) 
      enddo 


#ifdef MPI
       do i=1,khorl
       ii=lzet(i)
#else
      do ii=1,khor
#endif
        T_lim=(3.*abs(tl(ii))+5.)/( 1.+.1*abs(tl(ii)) )
        if(abs(tl10(ii)).gt.T_lim) then
        print*,'test T_lim'
        cdee=(0.824+0.041*wspeed(ii))*1.e-3
        hsens(ii)= cdee*rhol(ii)*wspeed(ii)*(-tl(ii)+tec(id3sur(ii)))
        hsens(ii)= hsens(ii)*cp
c	  hsens(ii)=sign( max(eps,abs(hsens(ii))),hsens(ii))
        evap(ii) = cdee*rhol(ii)*wspeed(ii)*(-sphumi(ii)+spsur(ii))
        tl10(ii) = tl(ii) 
        q10(ii)  = sphumi(ii) 
        endif
	enddo 

#ifndef MPI
c      if(loxstat1)then
	  write(*,*)'3 After T_lim'
        call xstat(evap,khor,zmax,zmin,zmean,zsum)
	  print*,'evap min,max,mean ',zmin,zmax,zmean
c      endif
#endif 
 
c     turbulent heat fluxes were already calculated in the loop, 
c     longwave radiation in the main  
#ifdef MPI
       do ii=1,khorl
       i=lzet(ii)
#else
      do i=1,khor 
#endif
        qw(i)  = -evap(i)*xlat- hsens(i) 
        fqs(i) = hsens(i)*(-1.0) 
        fql(i) = evap(i)*(-1.0)*xlat 
        q10(i) = 10.*pac(i)*q10(i)/(0.622+0.378*q10(i)) 
      enddo 
c      print*,'maximale Iteration in fluxes:',imax, 
c     +  'Anzahl Pkte. ohne Abbruch',sum  
c      write(7,*) 'max Iteration in fluxes:',imax, 
c     +  'Anzahl Pkte. ohne Abbruch',sum 
      
#ifndef MPI	
      if(loxstat)then
        write(*,*)'--------in the end of fluxes---------'
        call xstat(qw,khor,zmax,zmin,zmean,zsum)
        print*,'qw  min,max,mean ',zmin,zmax,zmean
        call xstat(fqs,khor,zmax,zmin,zmean,zsum)
        print*,'fqs min,max,mean ',zmin,zmax,zmean
        call xstat(fql,khor,zmax,zmin,zmean,zsum)
        print*,'fql min,max,mean ',zmin,zmax,zmean
        call xstat(q10,khor,zmax,zmin,zmean,zsum)
        print*,'q10 min,max,mean ',zmin,zmax,zmean
        write(*,*)'evap:'
        call xstat(evap,khor,zmax,zmin,zmean,zsum)
        print*,'evap min,max,mean ',zmin,zmax,zmean
        write(*,*)'wind:'
        call xstat(uw,khor,zmax,zmin,zmean,zsum)
        print*,'uw min,max,mean ',zmin,zmax,zmean
        call xstat(vw,khor,zmax,zmin,zmean,zsum)
        print*,'vw min,max,mean ',zmin,zmax,zmean
        write(*,*)'-------- End of output, fluxes---------'
        endif
#endif
      return 
      end 
      
c-----------------------------------------------------------------------
      subroutine trmice(iindex,izet,jc,mm,nn,wgesch,bewoel,
     * taup,einstrice,einstr,taair,dtsec,dz,iilo,ltief,qw,kkhor,
     * qoi,qii)
c-----------------------------------------------------------------------
c     thermodynamic ice-model
c-----------------------------------------------------------------------
      include 'C_model.f'
      include 'C_index.f'
      parameter(nx=m,ny=n)
      parameter(khor1=khor+1)  
      
      dimension iindex(mm,nn),izet(mm,nn),jc(mm,nn)
      dimension wgesch(kkhor),bewoel(kkhor),taup(kkhor)
      dimension einstrice(mm,nn),dz(iilo),taair(kkhor)
      dimension einstr(mm,nn),ltief(mm,nn)
      dimension qw(kkhor),qoi(mm,nn),qii(mm,nn)
c      dimension exbio(ndrei)
      real tw(m,n)
      dimension cw(khor)
          
      include 'C_ice.f' 
#ifdef MPI
      common /cord/ I1,I2,J1,J2,IA1,IA2,JA1,JA2,IB1,IB2,JB1,JB2
     1,ID1,ID2,JD1,JD2
      common /icecor/ ICEI1,ICEI2,ICEJ1,ICEJ2,ICEIA1,ICEIA2,ICEJA1,
     * ICEJA2,ICEX1,ICEX2,ICEY1,ICEY2,ICEP1,ICEP2,ICEQ1,ICEQ2,ICEX3,
     * ICEX4,ICEY3,ICEY4,ICEX5,ICEX6,ICEY5,ICEY6
      common /lbkhor/ khorl,lzet(khor),lb0(n),le0(n),
     * lb1(n),le1(n),lb2(n),le2(n)
#endif

c      print*,'i trmice2',maxval(tec),minval(tec),maxval(sac),minval(sac)  

c--cold stronger absprbtion at the surface (cd3)
         xnue1=14.6
         xnue2=0.24
c-- warm less absorbtion at the surface (nc1/cd2)
cwarm       xnue1 = 29.0
cwarm       xnue2 = 0.1392

       aten=(xnue1*xnue2-xnue1)/(xnue2-xnue1)

       dpmax=40.
       avv= 50.1e-4
 
c---------------------------------------------------------------------- 
c start space-loop 
c--------------------------------------------------------------------- 

#ifdef MPI
      do k = J1,J2
      do lw = lb0(k),le0(k)
#else
      do k = 1,n
      do lw = lb(k),le(k)
#endif
      nwet = indwet(lw)
      i = iwet(lw)

c-------------------------------------------------------------------- 
c shortwave radiation depending on the latitude 
c for ice: frsi,    for water: frsw 
c----------------------------------------------- -------------------
c------------------------------------------------------------------ 
c updating temperatures 
c------------------------------------------------------------------ 
      te = tec(nwet+1)  
      sa = sac(nwet+1) 
      tfr = tfreez(sa,0.) 
      xtl = taair(lw)
      tf=tfreez(sa/5.,0.)
      tis(i,k)=amin1(tf,xtl)
      ddzz=dz(1)+zac(lw)
      dzw=ddzz
      dpenetr=ddzz
      if(lazc(izet(i,k)).eq.1) then
        ddzz=zac(izet(i,k))+ldep(izet(i,k))
        dpenetr=dpmax
      endif

c----------------------------------------------------------------- 
c windspeed 
c------------------------------------------------------------------ 

      wspeed=wgesch(lw)
      csens=0.0 
      sigem=0.97*5.67e-8

c-------------------------------------------------------------------- 
c ratio : density ice / density water 
c------------------------------------------------------------------ 

      roi =rois*0.001
      row=(stc(nwet+1)+rhoq(1)) 
      cp1=4.2198-6.6*sa/1000.
      cp2=4.08-2.47*sa/1000.
      atem=te/20.
      atem=min(atem,1.)
      atem=max(atem,0.)
      btem=1-atem
      cp=atem*cp2+btem*cp1
      cp=cp*1000.
      cw(lw)=cp*row*1000.
c---------------------------------------------------------------------------- 
c start ice model (ice or no ice) 
c---------------------------------------------------------------------------- 
      xmelt = 0.0 
      delti = 0.0 
      delhi = 0.0 
      if((frice(i,k).le.epsis).or.
     + (his(i,k).le.epsis)) then 

c no ice, skip the ice calculation 

      xmelt = 0.0 
      delti = 0.0 
      delhi = 0.0 
      else
c---------------------------------------------------------------------- 
c there is ice on water, determine  heat-fluxes, 
c for new ice surface temperature tis(i,k) 
c using the old ice surface temperatures tp 
c--------------------------------------------------------------------- 

      hisa = his(i,k)/frice(i,k)
      tcih = tci/hisa 

c-------------------------------------------------------------------- 
c longwave radiation frl...... 
c-------------------------------------------------------------------- 

      xtaup=taup(lw)
ccc      xtaup=xtaup/100.    ! depends on units 
      ea=esata(xtaup)
cc      xtaup*esati(tf)      ! consider units 
      pp=1013.
      bew=bewoel(lw) 

c--------------------------------------------------------------------
c   determine short wave radiation
c--------------------------------------------------------------------

       frsi=einstrice(i,k)

c------------------------------------------------------------------- 
c determine new heat conduction through the ice fc.... 
c-------------------------------------------------------------------
 
      fc = tcih*(tfr-tis(i,k)) 
      qii(i,k) = fc

c------------------------------------------------------------------- 
c and fw (ice-ocean heat flux) 
c-------------------------------------------------------------------
 
      fw = cw(lw)*avv*(tfr-te)/ddzz 
      qoi(i,k) = fw 

c------------------------------------------------------------------- 
c ocean temperature change due to heat loss to the ice 
c-------------------------------------------------------------------
 
      delti = dt*fw/(ddzz*cw(lw)) 

c------------------------------------------------------------------- 
c ice production or ice melting at the bottom 
c due to (im-)balance between heat flux through the 
c ice and oceanic heat loss (fc and fw) 
c-------------------------------------------------------------------
 
      delhi =  dtroil*(fc+fw) 

c-------------------------------------------------------------------
c ice melting at the surface of the ice (at tf!) 
c using heat fluxes: frl,fs,fl,f,frsi 
c-------------------------------------------------------------------
 
      frl = rback(xtl,tis(i,k),bew,sigem) 
      fs = sens(wspeed,xtl,tis(i,k),ea,
     + esati(tis(i,k)),csens) 
      fl = evapa(wspeed,xtl,tis(i,k),xtaup,
     + esati(tis(i,k)),pp)
      fl=fl/(1.2*2.5)
      fl=fl*0.55*2.835

c------------------------------------------------------------------- 
c melting rate xmelt (only negative --> melting!) 
c-------------------------------------------------------------------
 
       xs =(fs + fl + frl - frsi-fc) * dtroil 
      xmelt = amin1(xs,0.0) 
cc      if(i.eq.17.and.k.eq.158)then 
cc      print*,'fs,frl,fl,frsi,xmelt im eis' 
cc      print*,fs,frl,fl,frsi,xmelt,fw,fc,delhi,delti 
cc      endif 
      end if

c------------------------------------------------------------------------- 
c there is open water, determine atmosph.-ocean heat-flux
c determine temperature change due to heat-flux (loss) 
c to the  atmosphere and consider temperature 
c change due to short wave radiation in the sea surface layer 
c-------------------------------------------------------------------------- 
c      qocean = fs + frl + fl + frsw

      tw(i,k)=qw(lw)+aten/xnue1*
     + (1.-exp(-xnue1*dpenetr))*einstr(i,k)+
     + (1.-aten)/xnue2*
     + (1.-exp(-xnue2*dpenetr))*einstr(i,k)
      deltw = tw(i,k)* dt/(ddzz*cw(lw)) 
 

c----- add biology--------------------------------------------------to test
c	  xnue1=0.08+exbio(iindex(i,k)+1)
c      tw(i,k)=qw(izet(i,k))+1./xnue1*
c     + (1.-exp(-xnue1*ddzz))*einstr(i,k)
c      deltw = tw(i,k)* dt/(ddzz*cw(lw)) 
c----- end of add biology-------------------------------------------to test
  

c--------------------------------------------------------------------------- 
c new water temperature in surface layer 
c due to heat loss to the ice and to the atmosphere 
c--------------------------------------------------------------------------- 

      waneu=te+frice(i,k)*delti+ 
     +        (1.0-frice(i,k))*deltw 

c--------------------------------------------------------------------------- 
c initial ice thickness in open water 
c---------------------------------------------------------------------------
 
      temdel = tfr - waneu 
      eisneu = ddzz*cw(lw)*amax1(0.0,temdel)/roil 

c---------------------------------------------------------------------------
c ice thickness change due to 
c freezing and melting at the bottom (delhi) 
c and melting at the top (xmelt) 
c---------------------------------------------------------------------------
 
      delhis = delhi + xmelt 

c--------------------------------------------------------------------------- 
c average ice growth rate (hisdel) 
c---------------------------------------------------------------------------
 
      hisdel = frice(i,k)*delhis +(1.-frice(i,k))*eisneu 
      hisneu = his(i,k) + hisdel 

c--------------------------------------------------------------------------- 
c change of ice compactness due to freezing or melting 
c --------------------------------------------------------------------------
 
      teiler = amax1(hisneu,0.1) 
      schmel = amin1(delhis,0.0) 
      frier = eisneu / 0.5 
      delfrs = (frice(i,k)/(2.0*teiler)) * schmel 
      delfri = (1.0-frice(i,k)) * frier 
cc      delfra = dt * (delfrs + delfri) 
      delfra = (delfrs + delfri) 

c--------------------------------------------------------------------------- 
c salt increment due to thermodynamic 
c change of ice-thickness
c sice= 20% of sea salinity is remaining in the ice 
c attention: ice salinity is not a state variable   
c---------------------------------------------------------------------------
 
      sice=sa/5.
      delsal=(row*sa*dzw-roi*sice*hisdel)/(row*dzw-roi*hisdel)
      delsal=amax1(1.,delsal)

c---------------------------------------------------------------------------- 
c new variables 
c---------------------------------------------------------------------------
 
      frice(i,k) = frice(i,k) + delfra 
      his(i,k) = hisneu 
      tec(nwet+1) = amax1(tfr,waneu)  
      sac(nwet+1) = delsal
      end do
      end do
c------------------------------------------------------------------------------
c  change lower layer temperature here, short wave radition is penetrating into
c  the water column
c------------------------------------------------------------------------------
c
#ifdef MPI
      do k = J1,J2
      lwa = lb0(k)
      lwe = le0(k)
#else
      do k = 1,n
      lwa = lb(k)
      lwe = le(k)
#endif
      if (lwa.le.lwe) then
         llb = indwet(lwa)+1
         lle = indwet(lwe)+lazc(lwe)
         do ll = llb,lle
            j = ll-indwet(llw(ll))
            lw = llw(ll)
            i = iwet(lw)
            nwet = indwet(lw)
            ltt=min(9,lazc(lw))
      if (j.gt.1.and.j.le.ltt) then
       z0=zac(lw)+dz(j-1)
       dick=(dz(j)-dz(j-1))
       za=zac(lw)+dz(j)
       if(j.eq.ltt) za=max(za,dpmax)
       tec(nwet+j) =tec(nwet+j)+
     +  (1.0-frice(i,k))*
     +  aten*dt/(dick*cw(lw)*xnue1)*
     + (exp(-xnue1*z0)-exp(-xnue1*za))*einstr(i,k)+
     +  (1.0-frice(i,k))*
     +  (1.-aten)*dt/(dick*cw(lw)*xnue2)*
     + (exp(-xnue2*z0)-exp(-xnue2*za))*einstr(i,k)
      end if
c------- add biology------------------
c      z0=0.
c      za=ddzz
c      ltt=min(9,ltief(i,k))
c      do it=2,ltt
c       xnue1= 0.08+exbio(iindex(i,k)+it)
c	    z0=za
c       dick=(dz(it)-dz(it-1))
c       za=za+dick
c       tec(iindex(i,k)+it) =tec(iindex(i,k)+it)+
c     +  (1.0-frice(i,k))*
c     +  dt/(dick*cw(lw)*xnue1)*
c     + (exp(-xnue1*z0)-exp(-xnue1*za))*einstr(i,k)
c  
c       enddo
c----- end of add biology-----------------
            end do
         end if
      end do

#ifdef MPI
      do k = j1,j2
!CDIR VECTOR NODEP
      do lw = lb0(k),le0(k)
#else
!CDIR VECTOR NODEP
      do lw = lb(1),le(n)
      k = jwet(lw)
#endif
      i = iwet(lw)
      if((his(i,k).lt.epsis).or.(frice(i,k).lt.epsis)) then 
      his(i,k) = 0.0 
      frice(i,k) = 0.0 
      endif 
      end do
#ifdef MPI
      end do
#endif
      return 
      end 
c-----------------------------------------------------------------------
      function evapa(u,ta,ts,ea,es,pp)
c-----------------------------------------------------------------------
c    evaporation: currently only used for estimation of evaporation of ice 
c-----------------------------------------------------------------------

      us=u 
      if (u.le.0.3) u=0.3 
      if (u.gt.50.) u=50. 
c-----wind speed classes  
      if (u.le.2.2) goto 20 
      if (u.le.5.0) goto 30 
      if (u.le.8.0) goto 40 
      if (u.le.25.0) goto 50 
c-----w > 25. 
      a=1.68 
      b=-0.016 
      c=0.0 
      p=1. 
      goto 60 
c-----w < 25. 
50    continue 
      a=1.196 
      b=0.008 
      c=-0.0004 
      p=1. 
      goto 60 
c-----w < 8. 
40    continue 
      a=1.18 
      b=0.01 
      c=0.0 
      p=1.0 
      goto 60 
c-----w < 5. 
30    continue 
      a=0.969 
      b=0.0521 
      c=0.0 
      p=1.0 
      goto 60 
c-----w < 2.2 
20    continue 
      a=0.0 
      b=1.23 
      c=0.0 
      p=-0.16 
c-----transport coefficient and stability parameter  
60    continue 
ccc      ea=ha*(esata(ta))
      ce=a+b*(u**p)+c*((u-8.)**2) 
      snul=((ts-ta)+0.0011*(es-ea))/(u**2) 
      s=snul*abs(snul)/(abs(snul)+0.01) 
      if (s.gt.0.0) goto 80 
c-----stable 
      if (s.lt.-3.3) s=-3.3 
      ce=ce*(0.1+0.03*s+0.9*exp(4.8*s))/1000. 
      goto 100 
c-----unstable  
80    continue 
      ce=ce*(1.0+0.63*sqrt(s))/1000. 
c-----latent flux :  instead spec. humidity  --> ro w-dampf] 
100   continue 
      u=us 
      ce=1.6
      ra=ea/(461.*(ta+273.15)) 
      rs=es/(461.*(ts+273.15)) 
      evapa=ce*u*(rs-ra)*(2.5008-0.00272*ta)*1000. 


      return 
      end 



c-----------------------------------------------------------------------
      function sens(u,ta,ts,ea,es,csens)
c-----------------------------------------------------------------------
c   only used to estimate the sensible heat flux to the ice surface    
c-----------------------------------------------------------------------

c-----windspeed tresholds and classes 
      us=u 
      tss=273.15+ts
      taa=273.15+ta
      if (u.le.0.3) u=0.3 
      if (u.gt.50.) u=50. 

      if (u.le.2.2) then
c-----w < 2.2 
         a=0. 
         b=1.185 
         c=0. 
         p=-0.157 
      else if (u.le.5.) then
c-----w < 5 
         a=0.927 
         b=0.0546 
         c=0. 
         p=1. 
      else if (u.le.8.) then
c-----w < 8 
         a=1.15 
         b=0.01 
         c=0. 
         p=1. 
      else if (u.le.25.) then
c-----w < 25 
         a=1.17 
         b=0.0075 
         c=-0.00045 
         p=1. 
      else
c-----w > 25 
         a=1.652 
         b=-0.017 
         c=0.0 
         p=1. 
      end if
c-----transfer coefficients and stability
  
      ch=a+b*(u**p)+c*((u-8.)**2) 
      snul=((ts-ta)+0.0011*(es-ea))/(u**2) 
      s=snul*abs(snul)/(abs(snul)+0.01) 
      if (s.gt.0.) then

c-----unstable  
         ch=ch*(1.0+0.63*sqrt(s))/1000. 
      else
c-----stable
 
         if (s.lt.-3.3) s=-3.3 
         ch=ch*(0.1+0.03*s+0.9*exp(4.8*s))/1000. 
      end if

c-----heat flux: =cl * rol * ch * u *(twater - tair)
 
      u=us 
      csens=1.28*1015.*ch*u 
      sens=csens*(tss-taa) 


      return 
      end 


c-----------------------------------------------------------------------
      function rback(ta,ts,cn,sigem) 
c-----------------------------------------------------------------------
c     longwave radiation after Maykut
c-----------------------------------------------------------------------
 
*-----air temperature in kelvin 
      t=273.15+ta 
      ta4=t*t*t*t 
*-----SST in kelvin 
      t=273.15+ts 
      ts4=t*t*t*t 
*-----clouds  1/8 cn**2 = 1  
      cn1=cn / 8. 
*-----emissivity of water and ice    
      emiss=0.97 
*-----emissivity atmosphere (clouds!) 
      emisa=0.7855*(1.+(0.2232*cn1)**2.75) 
*     stephan boltzmann konst 
      boltz=5.67e-8 
cc      sigem = emiss*boltz 
*-----long wave net radiation 
       rback=emiss*boltz*ts4 - emisa*boltz*ta4 


       return 
       end 

c-----------------------------------------------------------------------
      function esati(tt) 
c-----------------------------------------------------------------------
c     saturation of vapor pressure at water temperature t 
c-----------------------------------------------------------------------
      
      ai    = (9.5*tt)/(265.5+tt) 
      esati = 610.7*(10.**ai) 


      return 
      end 

c-----------------------------------------------------------------------
      real function tfreez (s,p) 
c-----------------------------------------------------------------------
c     freezing point of sea water (millero 1978) 
c-----------------------------------------------------------------------

      parameter (ca=-0.0575,cb=1.710523e-3,cc=-2.154996e-4,cd=-7.53e-3)
      tfreez = s*(ca+cb*sqrt(s)+cc*s)+cd*p 


      return 
      end 

c-----------------------------------------------------------------------
      function esata(tt) 
c-----------------------------------------------------------------------
c     saturation of vapor pressure at atm. temperature t 
c-----------------------------------------------------------------------

      ai    = (7.5*tt)/(237.3+tt) 
      esata = 610.7*(10.**ai) 


      return 
      end 






c----------------------------------------------------------------------- 
      subroutine fluxes_xt3(ind2d,tl,sphumi,bew,pac,uw,vw,
     &   uc,vc,tec,qw ,id3sur,wspeed,q10,tl10,
     &    fqs,fql,cddu,evap) 

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	
ccc         Calculation of sea surface fluxes depending on stability of the 
ccc         atmospheric boundary layer following Timo Vihma (1995): 
ccc         Atmosphere-surface interactions over polar oceans and 
ccc         heterogenous surfaces, Finnish Marine Research No. 264
ccc                        C. Schrum        March 1997/Feb 2002
ccc                        I. Alekseeva
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
      include 'C_model.f'
      dimension ind2d(m,n),id3sur(khor) 
      dimension  tl(khor),taup(khor),bew(khor) 
      dimension  uw(khor),vw(khor),pac(khor) 
   
      dimension  qw(khor),wspeed(khor) 
      dimension tec(ndrei),uc(ndrei),vc(ndrei) 
      dimension q10(khor),tl10(khor) 
 
      dimension fqs(khor),fql(khor),cddu(khor) 
  
      real phihq(khor),phim(khor)  
      real  cdd(khor),cdeh(khor),rhol(khor) 
      real  zro(khor),zrtq(khor),evap(khor) 
      real uscale(khor),tscale(khor),qscale(khor) 
      real psim(khor),psihq(khor),diff(khor),tau(khor) 
      real sspeed(khor),hsens(khor),rl(khor) 
      real sphumi(khor),spsur(khor)
	real eps 
      integer iterat(khor) ,ab_tscale(khor)

      logical loxstat,loxstat1

      loxstat=.true.
      loxstat1=.false.

c
c---------------------------------------------------------------------------------
c   set level of wind measurements zw, level of humidity measurements zh,
c   level of air-temperature measurements zt
c---------------------------------------------------------------------------------
c
      zw   =   10. 
      zth  =    2. 
      xk   =    0.4 
      cp   = 1005. 
      xlat = 2500.*1000. 
      g    =    9.81 
      bound=    0.1 
      eps  =    1.e-7 
c	 
c--------------------------------------------------------------------------------
c   calculation of the density  of air, specific humidity, wind speed
c   and set phim, phih and phiq eq. 1 for neutral stability
c--------------------------------------------------------------------------------
c
      do i=1,khor 
        diff(i) = 1000. 
 
c          ee = exp(-6763.6/(taup(i)+273.15)- 
c     *     4.9283*log(taup(i)+273.15)+54.23) 
c          if(taup(i).lt. 0.)then 
c            ee = exp(-6141./(taup(i)+273.15)+24.3) 
c          endif 
 



	 if(tec(id3sur(i)).ge.0.)then 
	  esur = exp(-6763.6/(tec(id3sur(i))+273.15)- 
     &  4.9283*log(tec(id3sur(i))+273.15)+54.23) 
       else        
	  esur = exp(-6141./(tec(id3sur(i))+273.15)+24.3) 
	 endif 
     
c        sphumi(i) = 0.622*ee/(10.*pac(i)-0.378*ee) 
 
        spsur(i) = 0.622*esur /(10.*pac(i)-0.378*esur) 
 
        rhol(i) = 10.*pac(i)/(2.8705*(tl(i)+273.15)* 
     *   (1.+0.608*sphumi(i))) 
 
	wspeed(i) = sqrt(uw(i)**2+vw(i)**2) 
	sspeed(i) = 0. 
	psim(i)   = 1. 
	psihq(i)  = 1. 
	phim(i)   = 0. 
	iterat(i) = 1 
	uscale(i) = 1.
	phihq(i)=0.
 
      enddo 
c
c--------------------------------------------------------------------------------
c     set roughness lengths zro, zrtq  from V(10m) and first guess for u*
c     t* and q*
c--------------------------------------------------------------------------------
c
      do i=1,khor 
        cdbb = (0.8+0.065*wspeed(i))*1.e-3 
        cdee=(0.824+0.041*wspeed(i))*1.e-3
 
        cdd(i)  = cdbb 
        cdeh(i) = cdee 
        zro(i)  = exp(log(10.)-xk/sqrt(cdbb)) 
        zrtq(i) = exp(log(2.)-xk*sqrt(cdbb)/cdee) 
        uscale(i) = cdd(i)*(wspeed(i))**2 
        uscale(i) = max(sqrt(uscale(i)),eps) 
        tscale(i) = (tl(i)-tec(id3sur(i)) ) 
        tscale(i) = cdeh(i)*tscale(i)*wspeed(i) 
        tscale(i) = tscale(i)/(xk*uscale(i)) 
        qscale(i) = (sphumi(i) -spsur(i)) 
        qscale(i) = cdeh(i)*qscale(i)*wspeed(i) 
        qscale(i) = qscale(i)/(xk*uscale(i)) 

        ab_tscale(i)=abs( tscale(i) )
      enddo 

      if(loxstat1)then
	  write(*,*)'1 before iteration'
        call xstat(tscale,khor,zmax,zmin,zmean,zsum)
	  print*,'tscale  min,max,mean,sum ',zmin,zmax,zmean
        call xstat(ab_tscale,khor,zmax,zmin,zmean,zsum)
	  print*,'ab_tscale min,max,mean,sum ',zmin,zmax,zmean

        write(*,*)'sphumi:'
        call xstat(sphumi,khor,zmax,zmin,zmean,zsum)
	  print*,'sphumi min,max,mean,sum ',zmin,zmax,zmean

        write(*,*)'spsur:'
        call xstat(spsur,khor,zmax,zmin,zmean,zsum)
	  print*,'spsur min,max,mean,sum ',zmin,zmax,zmean

        call xstat(evap,khor,zmax,zmin,zmean,zsum)
	  print*,'evap min,max,mean,sum ',zmin,zmax,zmean
      endif
c------------------------------------------------------------------------------
c     start of the iteration
c-----------------------------------------------------------------------------
c
      iiicount=0
      imax=0 
      do i=1,khor   
       imax=max0(iiicount,imax) 
       iiicount=0 
       to=(tl(i)+tec(id3sur(i)))/2.+273.15
       do ii=1,40 
      do iit=1,iterat(i) 
        iiicount=iiicount+1 
c
c------------------------------------------------------------------------------
c     calculate new 10 m values for air temperature and humidity
c-----------------------------------------------------------------------------
c
        tl10(i) = tl(i)+tscale(i)*log(5.)*psihq(i) 
        q10(i)  = sphumi(i)+qscale(i)*log(5.)*psihq(i) 
c
c-------------------------------------------------------------------------------
c    calculate transfer coefficients for the input-data height
c--------------------------------------------------------------------------------
c       cdeh(i) = xk/( log(2./zrtq(i))+phihq(i) ) 
c        cdd(i)  = xk/(log(2./zro(i))+phim(i))  


           xlog =  log(2./zrtq(i)) +phihq(i) 
            xlog =sign (max(abs(xlog),eps),xlog)
        cdeh(i) = xk/xlog


          xlog =  log(2./zro(i))+phim(i) 
          xlog =sign (max(abs(xlog),eps),xlog)
        cdd(i)  = xk/xlog

        cdd(i)  = cdd(i)*cdd(i) 
        cdeh(i) = cdeh(i)*sqrt(cdd(i)) 
        cdee    = cdeh(i) 
 
c
c------------------------------------------------------------------------------
c    calculate fluxes and 2-m wind speed w2
c-----------------------------------------------------------------------------
c


        w2 = wspeed(i)-uscale(i)*log(5.)*psim(i)/xk 

	  hsens(i) = cdeh(i)*rhol(i)*w2*(-(tl(i)+3.)+tec(id3sur(i))) 
        hsens(i) = hsens(i)*cp
	  hsens(i)=sign( max(eps,abs(hsens(i))) ,hsens(i) )
	  
        evap(i)  = cdeh(i)*rhol(i)*w2*(-sphumi(i)+spsur(i))
c
c---------------------------------------------------------------------------------
c    calculate z/L  by using the fluxes
c--------------------------------------------------------------------------------
c

          rl(i) = 1.+0.61*to*cp*evap(i)/(hsens(i)) 
          rl(i) = g*xk*hsens(i)*rl(i) 

          u3=max(uscale(i)**3,eps) 
 
          rl(i) = rl(i)/(u3*to*rhol(i)*cp) 


          istab = 1 
c---------------------------------------------------------------------------------
c    calculate the universal functions phim/hq  and psim/hq 
c    for neutral stratification
c--------------------------------------------------------------------------------
        if( abs(rl(i)).le.eps)  then 
          istab    = 0 
          psihq(i) = 1. 
          phihq(i) = 0.  
          psim(i)  = 1. 
          phim(i)  = 0. 
        endif  
 

        do ixx=1,istab 
          istaby = 0 
	    istabn = 0 
          if(rl(i).lt.-eps) istaby = 1 
          if(rl(i).gt.eps) istabn = 1 

          
c
c---------------------------------------------------------------------------------
c    calculate the universal functions phim/hq  and psim/hq 
c    for stable stratification
c--------------------------------------------------------------------------------
c 
 	  
	   do ix=1,istaby 
ccc                 print*,'Im stabilen Bereich'
		zrl=-2.*rl(i)
		 phihq(i)=-0.7*zrl-0.75*(zrl-5./0.35)
     +		 *exp(-0.35*zrl)-0.75*5./0.35
                 zrl= -zro(i)*rl(i)
		 zrl=-0.7*zrl-0.75*(zrl-5./0.35)
     +		 *exp(-0.35*zrl)-0.75*5./0.35

		 zzrl=-10.*rl(i)
		 phim(i)=-0.7*zzrl-0.75*(zzrl-5./0.35)
     +		 *exp(-0.35*zzrl)-0.75*5./0.35
		 phim(i)=(phim(i)-zrl)
		 zrl=-zrtq(i)*rl(i)
		 zrl=-0.7*zrl-0.75*(zrl-5./0.35)
     +		 *exp(-0.35*zrl)-0.75*5./0.35

		 phihq(i)=(phihq(i)-zrl)
          
		 psihq(i)=rl(i)*5.*2.
		 psim(i)=psihq(i)
ccccc         print*,'hier',psim(i),psihq(i),rl(i) 
             enddo 
c
c---------------------------------------------------------------------------------
c    calculate the universal functions phim/hq  and psim/hq 
c    for unstable stratification
c--------------------------------------------------------------------------------
c 
	     do ix=1,istabn		  
c		psim(i) = sqrt(  1.+19.3*2.*max(0.00001,rl(i)) ) 
		psim(i) = sqrt(  1.+19.3*2.*rl(i) ) 
	    psim(i) = 1./sqrt(psim(i)) 
	    zz = sqrt((1.+19.3*zro(i)*rl(i))) 
	    zz = 1./sqrt(zz) 
	    phim(i) = 2.*log(0.5+1./(2.*psim(i))) 
     * 	     +log(0.5+1./(2.*psim(i)**2)) 
     *       -2.*atan(1./psim(i))+3.14/2. 
            phim(i) = phim(i)-2.*log(0.5+1./(2.*zz)) 
     * 	     -log(0.5+1./(2.*zz**2)) 
     *       +2.*atan(1./zz)-3.14/2. 
	    zophihq  = (1.+12.*zrtq(i)*rl(i)) 
	    zophihq  = 1./sqrt(zophihq) 
     
            psihq(i) = (1.+24.*rl(i)) 
            psihq(i) = 1./sqrt(psihq(i)) 
  		phihq(i) = 2.*log(0.5+1./(2.*psihq(i))) 
	    phihq(i) = phihq(i)-2.*log(0.5+1./(2.*zophihq)) 
          enddo 

        enddo !        do ixx=1,istab 


c       
c------------------------------------------------------------------------------
c     calculate new scaling parameters
c-----------------------------------------------------------------------------
c
        tscale(i) = -hsens(i)/(rhol(i)*cp*xk*uscale(i)) 
        qscale(i) = -evap(i)/(rhol(i)*xk*uscale(i)) 
c
c---------------------------------------------------------------------------------
c    continue with estimations of the scaling parameters if ......
c-------------------------------------------------------------------------------
c
        if(abs(diff(i)-10.*rl(i)).lt.bound) iterat(i) = 0 
        diff(i) = 10.*rl(i) 
 
      enddo 
      enddo 
      enddo 

c------------------------------------------------------------------------------
c     End of the iteration
c-----------------------------------------------------------------------------

c-----------------------------------------------------------------------------
c in case of no convergence [iterat(i).gt.0], apply neutral stability parameterization 
c-----------------------------------------------------------------------------
       sum=0 
       do ii=1,khor 
       sum=sum+iterat(ii) 
       if (iterat(ii).gt.0) then 
        cdee=(0.824+0.041*wspeed(ii))*1.e-3
*xt3--------------------------------------------------
      hsens(ii)=cdee*rhol(ii)*wspeed(ii)*(-(tl(ii)+3.)+tec(id3sur(ii)))
        hsens(ii) = hsens(ii)*cp
	  hsens(ii)=sign( max(eps,abs(hsens(ii))),hsens(ii))


        evap(ii)  = cdee*rhol(ii)*wspeed(ii)*(-sphumi(ii)+spsur(ii))

        tl10(ii) = tl(ii) 
        q10(ii)  = sphumi(ii) 
      end if 
      enddo  
    
      if(loxstat1)then
	  write(*,*)'3 no convergence '
        write(*,*)'evap:'
        call xstat(evap,khor,zmax,zmin,zmean,zsum)
	  print*,'evap min,max,mean ',zmin,zmax,zmean
 
 
         write(*,*)'hsens:'
        call xstat(hsens,khor,zmax,zmin,zmean,zsum)
	  print*,'hsense min,max,mean ',zmin,zmax,zmean
      endif


c----------- windstress ---------------------------------
      do i=1,khor 
 
c================================================================== 
c       Drag coefficient Cd (cdd) for neutral stability in the atmosphere 
c================================================================== 
*1-----------------------------
c       Charnock (linear form) (1955) 
c        cdd(i) = (0.7+0.09*wspeed(i))*1.e-3 
*2----------------------------- 
c       Smith and Banke (1975) 
c        cdd(i) = (0.63+0.066*wspeed(i))*1.e-3 
*3----------------------------- 
c       Wu (1980)  
c       cdd(i) = (0.8+0.065*wspeed(i))*1.e-3 
 

*4-----------------------------
c      cdd mit Stabilitaetsansatz f�r die Atmosph�re 
c      ansatz nach J. Launianien and T. Vihma (1990) 


c         if(iterat(i).gt.0.and.wspeed(i).gt.22.)then   !Irina
         if(iterat(i).gt.0.or.wspeed(i).gt.18.)then   !Irina
c        if(iterat(i).gt.0)then                        !Irina
c------in case of no convergence, apply neutral stability parameterization, Wu (1980)
	 cdd(i) = (0.8+0.065*wspeed(i))*1.e-3   
	 else
          xlog=log(10./zro(i))-phim(i)
          xlog =sign(max(abs(xlog),eps),xlog)
c       cdd(i) = (xk/(log(10./zro(i))-phim(i)))**2.0 
       cdd(i) = (xk/xlog)**2
       endif
 
c================================================================== 
 
 
        cddu(i) = cdd(i) 


*-------- wind stress-------------------
        uw(i)  = cdd(i)*rhol(i)*(wspeed(i))*(uw(i)) 
        vw(i)  = cdd(i)*rhol(i)*(wspeed(i))*(vw(i)) 
      enddo 



      do ii=1,khor
        T_lim=(3.*abs(tl(ii))+5.)/( 1.+.1*abs(tl(ii)) )
        if(abs(tl10(ii)).gt.T_lim) then
        cdee=(0.824+0.041*wspeed(ii))*1.e-3
*xt3--------------------------------------------------
      hsens(ii)=cdee*rhol(ii)*wspeed(ii)*(-(tl(ii)+3.)+tec(id3sur(ii)))
        hsens(ii)= hsens(ii)*cp
c	  hsens(ii)=sign( max(eps,abs(hsens(ii))),hsens(ii))
        evap(ii) = cdee*rhol(ii)*wspeed(ii)*(-sphumi(ii)+spsur(ii))
        tl10(ii) = tl(ii) 
        q10(ii)  = sphumi(ii) 
        endif
	enddo 


      if(loxstat1)then
	  write(*,*)'3 After T_lim'
        call xstat(evap,khor,zmax,zmin,zmean,zsum)
	  print*,'evap min,max,mean ',zmin,zmax,zmean
      endif
 
 
c     turbulent heat fluxes were already calculated in the loop, 
c     longwave radiation in the main  
      do i=1,khor 
        qw(i)  = -evap(i)*xlat- hsens(i) 
        fqs(i) = hsens(i)*(-1.0) 
        fql(i) = evap(i)*(-1.0)*xlat 
        q10(i) = 10.*pac(i)*q10(i)/(0.622+0.378*q10(i)) 
      enddo 
c      print*,'maximale Iteration in fluxes:',imax, 
c     +  'Anzahl Pkte. ohne Abbruch',sum  
c      write(7,*) 'max Iteration in fluxes:',imax, 
c     +  'Anzahl Pkte. ohne Abbruch',sum 
      
	
	
	
	
	
	
      if(loxstat)then
	  write(*,*)'--------in the end of fluxes---------'
        call xstat(qw,khor,zmax,zmin,zmean,zsum)
	  print*,'qw  min,max,mean ',zmin,zmax,zmean
        call xstat(fqs,khor,zmax,zmin,zmean,zsum)
	  print*,'fqs min,max,mean ',zmin,zmax,zmean
        call xstat(fql,khor,zmax,zmin,zmean,zsum)
	  print*,'fql min,max,mean ',zmin,zmax,zmean
        call xstat(q10,khor,zmax,zmin,zmean,zsum)
	  print*,'q10 min,max,mean ',zmin,zmax,zmean
        write(*,*)'evap:'
	  call xstat(evap,khor,zmax,zmin,zmean,zsum)
	  print*,'evap min,max,mean ',zmin,zmax,zmean
        write(*,*)'wind:'
        call xstat(uw,khor,zmax,zmin,zmean,zsum)
	  print*,'uw min,max,mean ',zmin,zmax,zmean
        call xstat(vw,khor,zmax,zmin,zmean,zsum)
	  print*,'vw min,max,mean ',zmin,zmax,zmean
	  write(*,*)'-------- End of output, fluxes---------'

      endif

	
	
	
	
	
	
	
	
	
	return 
      end




c----------------------------------------------------------------------- 
      subroutine fluxes_xt4(ind2d,tl,taup,bew,pac,uw,vw,
     &   uc,vc,tec,qw ,id3sur,wspeed,q10,tl10,
     &    fqs,fql,cddu,evap) 

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	
ccc         Calculation of sea surface fluxes depending on stability of the 
ccc         atmospheric boundary layer following Timo Vihma (1995): 
ccc         Atmosphere-surface interactions over polar oceans and 
ccc         heterogenous surfaces, Finnish Marine Research No. 264
ccc                        C. Schrum        March 1997/Feb 2002
ccc                        I. Alekseeva
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
      include 'C_model.f'
      dimension ind2d(m,n),id3sur(khor) 
      dimension  tl(khor),taup(khor),bew(khor) 
      dimension  uw(khor),vw(khor),pac(khor) 
   
      dimension  qw(khor),wspeed(khor) 
      dimension tec(ndrei),uc(ndrei),vc(ndrei) 
      dimension q10(khor),tl10(khor) 
 
      dimension fqs(khor),fql(khor),cddu(khor) 
  
      real phihq(khor),phim(khor)  
      real  cdd(khor),cdeh(khor),rhol(khor) 
      real  zro(khor),zrtq(khor),evap(khor) 
      real uscale(khor),tscale(khor),qscale(khor) 
      real psim(khor),psihq(khor),diff(khor),tau(khor) 
      real sspeed(khor),hsens(khor),rl(khor) 
      real sphumi(khor),spsur(khor)
	real eps 
      integer iterat(khor) ,ab_tscale(khor)

      logical loxstat,loxstat1

      loxstat=.true.
      loxstat1=.false.

c
c---------------------------------------------------------------------------------
c   set level of wind measurements zw, level of humidity measurements zh,
c   level of air-temperature measurements zt
c---------------------------------------------------------------------------------
c
      zw   =   10. 
      zth  =    2. 
      xk   =    0.4 
      cp   = 1005. 
      xlat = 2500.*1000. 
      g    =    9.81 
      bound=    0.1 
      eps  =    1.e-7 
c	 
c--------------------------------------------------------------------------------
c   calculation of the density  of air, specific humidity, wind speed
c   and set phim, phih and phiq eq. 1 for neutral stability
c--------------------------------------------------------------------------------
c
      do i=1,khor 
        diff(i) = 1000. 
 
          ee = exp(-6763.6/(taup(i)+273.15)- 
     *     4.9283*log(taup(i)+273.15)+54.23) 
          if(taup(i).lt. 0.)then 
            ee = exp(-6141./(taup(i)+273.15)+24.3) 
          endif 
 
        sphumi(i) = 0.622*ee/(10.*pac(i)-0.378*ee) 



	 if(tec(id3sur(i)).ge.0.)then 
	  esur = exp(-6763.6/(tec(id3sur(i))+273.15)- 
     &  4.9283*log(tec(id3sur(i))+273.15)+54.23) 
       else        
	  esur = exp(-6141./(tec(id3sur(i))+273.15)+24.3) 
	 endif 
     
 
        spsur(i) = 0.622*esur /(10.*pac(i)-0.378*esur) 
 
        rhol(i) = 10.*pac(i)/(2.8705*(tl(i)+273.15)* 
     *   (1.+0.608*sphumi(i))) 
 
	wspeed(i) = sqrt(uw(i)**2+vw(i)**2) 
	sspeed(i) = 0. 
	psim(i)   = 1. 
	psihq(i)  = 1. 
	phim(i)   = 0. 
	iterat(i) = 1 
	uscale(i) = 1.
	phihq(i)=0.
 
      enddo 
c
c--------------------------------------------------------------------------------
c     set roughness lengths zro, zrtq  from V(10m) and first guess for u*
c     t* and q*
c--------------------------------------------------------------------------------
c
      do i=1,khor 
        cdbb = (0.8+0.065*wspeed(i))*1.e-3 
        cdee=(0.824+0.041*wspeed(i))*1.e-3
 
        cdd(i)  = cdbb 
        cdeh(i) = cdee 
        zro(i)  = exp(log(10.)-xk/sqrt(cdbb)) 
        zrtq(i) = exp(log(2.)-xk*sqrt(cdbb)/cdee) 
        uscale(i) = cdd(i)*(wspeed(i))**2 
        uscale(i) = max(sqrt(uscale(i)),eps) 
        tscale(i) = (tl(i)-tec(id3sur(i)) ) 
        tscale(i) = cdeh(i)*tscale(i)*wspeed(i) 
        tscale(i) = tscale(i)/(xk*uscale(i)) 
        qscale(i) = (sphumi(i) -spsur(i)) 
        qscale(i) = cdeh(i)*qscale(i)*wspeed(i) 
        qscale(i) = qscale(i)/(xk*uscale(i)) 

        ab_tscale(i)=abs( tscale(i) )
      enddo 

      if(loxstat1)then
	  write(*,*)'1 before iteration'
        call xstat(tscale,khor,zmax,zmin,zmean,zsum)
	  print*,'tscale  min,max,mean,sum ',zmin,zmax,zmean
        call xstat(ab_tscale,khor,zmax,zmin,zmean,zsum)
	  print*,'ab_tscale min,max,mean,sum ',zmin,zmax,zmean

        write(*,*)'sphumi:'
        call xstat(sphumi,khor,zmax,zmin,zmean,zsum)
	  print*,'sphumi min,max,mean,sum ',zmin,zmax,zmean

        write(*,*)'spsur:'
        call xstat(spsur,khor,zmax,zmin,zmean,zsum)
	  print*,'spsur min,max,mean,sum ',zmin,zmax,zmean

        call xstat(evap,khor,zmax,zmin,zmean,zsum)
	  print*,'evap min,max,mean,sum ',zmin,zmax,zmean
      endif
c------------------------------------------------------------------------------
c     start of the iteration
c-----------------------------------------------------------------------------
c
      iiicount=0
      imax=0 
      do i=1,khor   
       imax=max0(iiicount,imax) 
       iiicount=0 
       to=(tl(i)+tec(id3sur(i)))/2.+273.15
       do ii=1,40 
      do iit=1,iterat(i) 
        iiicount=iiicount+1 
c
c------------------------------------------------------------------------------
c     calculate new 10 m values for air temperature and humidity
c-----------------------------------------------------------------------------
c
        tl10(i) = tl(i)+tscale(i)*log(5.)*psihq(i) 
        q10(i)  = sphumi(i)+qscale(i)*log(5.)*psihq(i) 
c
c-------------------------------------------------------------------------------
c    calculate transfer coefficients for the input-data height
c--------------------------------------------------------------------------------
c       cdeh(i) = xk/( log(2./zrtq(i))+phihq(i) ) 
c        cdd(i)  = xk/(log(2./zro(i))+phim(i))  


           xlog =  log(2./zrtq(i)) +phihq(i) 
            xlog =sign (max(abs(xlog),eps),xlog)
        cdeh(i) = xk/xlog


          xlog =  log(2./zro(i))+phim(i) 
          xlog =sign (max(abs(xlog),eps),xlog)
        cdd(i)  = xk/xlog

        cdd(i)  = cdd(i)*cdd(i) 
        cdeh(i) = cdeh(i)*sqrt(cdd(i)) 
        cdee    = cdeh(i) 
 
c
c------------------------------------------------------------------------------
c    calculate fluxes and 2-m wind speed w2
c-----------------------------------------------------------------------------
c


        w2 = wspeed(i)-uscale(i)*log(5.)*psim(i)/xk 

	  hsens(i) = cdeh(i)*rhol(i)*w2*(-tl(i)+tec(id3sur(i))) 
        hsens(i) = hsens(i)*cp
	  hsens(i)=sign( max(eps,abs(hsens(i))) ,hsens(i) )
	  
        evap(i)  = cdeh(i)*rhol(i)*w2*(-sphumi(i)+spsur(i))
c
c---------------------------------------------------------------------------------
c    calculate z/L  by using the fluxes
c--------------------------------------------------------------------------------
c

          rl(i) = 1.+0.61*to*cp*evap(i)/(hsens(i)) 
          rl(i) = g*xk*hsens(i)*rl(i) 

          u3=max(uscale(i)**3,eps) 
 
          rl(i) = rl(i)/(u3*to*rhol(i)*cp) 


          istab = 1 
c---------------------------------------------------------------------------------
c    calculate the universal functions phim/hq  and psim/hq 
c    for neutral stratification
c--------------------------------------------------------------------------------
        if( abs(rl(i)).le.eps)  then 
          istab    = 0 
          psihq(i) = 1. 
          phihq(i) = 0.  
          psim(i)  = 1. 
          phim(i)  = 0. 
        endif  
 

        do ixx=1,istab 
          istaby = 0 
	    istabn = 0 
          if(rl(i).lt.-eps) istaby = 1 
          if(rl(i).gt.eps) istabn = 1 

          
c
c---------------------------------------------------------------------------------
c    calculate the universal functions phim/hq  and psim/hq 
c    for stable stratification
c--------------------------------------------------------------------------------
c 
 	  
	   do ix=1,istaby 
ccc                 print*,'Im stabilen Bereich'
		zrl=-2.*rl(i)
		 phihq(i)=-0.7*zrl-0.75*(zrl-5./0.35)
     +		 *exp(-0.35*zrl)-0.75*5./0.35
                 zrl= -zro(i)*rl(i)
		 zrl=-0.7*zrl-0.75*(zrl-5./0.35)
     +		 *exp(-0.35*zrl)-0.75*5./0.35

		 zzrl=-10.*rl(i)
		 phim(i)=-0.7*zzrl-0.75*(zzrl-5./0.35)
     +		 *exp(-0.35*zzrl)-0.75*5./0.35
		 phim(i)=(phim(i)-zrl)
		 zrl=-zrtq(i)*rl(i)
		 zrl=-0.7*zrl-0.75*(zrl-5./0.35)
     +		 *exp(-0.35*zrl)-0.75*5./0.35

		 phihq(i)=(phihq(i)-zrl)
          
		 psihq(i)=rl(i)*5.*2.
		 psim(i)=psihq(i)
ccccc         print*,'hier',psim(i),psihq(i),rl(i) 
             enddo 
c
c---------------------------------------------------------------------------------
c    calculate the universal functions phim/hq  and psim/hq 
c    for unstable stratification
c--------------------------------------------------------------------------------
c 
	     do ix=1,istabn		  
c		psim(i) = sqrt(  1.+19.3*2.*max(0.00001,rl(i)) ) 
		psim(i) = sqrt(  1.+19.3*2.*rl(i) ) 
	    psim(i) = 1./sqrt(psim(i)) 
	    zz = sqrt((1.+19.3*zro(i)*rl(i))) 
	    zz = 1./sqrt(zz) 
	    phim(i) = 2.*log(0.5+1./(2.*psim(i))) 
     * 	     +log(0.5+1./(2.*psim(i)**2)) 
     *       -2.*atan(1./psim(i))+3.14/2. 
            phim(i) = phim(i)-2.*log(0.5+1./(2.*zz)) 
     * 	     -log(0.5+1./(2.*zz**2)) 
     *       +2.*atan(1./zz)-3.14/2. 
	    zophihq  = (1.+12.*zrtq(i)*rl(i)) 
	    zophihq  = 1./sqrt(zophihq) 
     
            psihq(i) = (1.+24.*rl(i)) 
            psihq(i) = 1./sqrt(psihq(i)) 
  		phihq(i) = 2.*log(0.5+1./(2.*psihq(i))) 
	    phihq(i) = phihq(i)-2.*log(0.5+1./(2.*zophihq)) 
          enddo 

        enddo !        do ixx=1,istab 


c       
c------------------------------------------------------------------------------
c     calculate new scaling parameters
c-----------------------------------------------------------------------------
c
        tscale(i) = -hsens(i)/(rhol(i)*cp*xk*uscale(i)) 
        qscale(i) = -evap(i)/(rhol(i)*xk*uscale(i)) 
c
c---------------------------------------------------------------------------------
c    continue with estimations of the scaling parameters if ......
c-------------------------------------------------------------------------------
c
        if(abs(diff(i)-10.*rl(i)).lt.bound) iterat(i) = 0 
        diff(i) = 10.*rl(i) 
 
      enddo 
      enddo 
      enddo 

c------------------------------------------------------------------------------
c     End of the iteration
c-----------------------------------------------------------------------------

c-----------------------------------------------------------------------------
c in case of no convergence [iterat(i).gt.0], apply neutral stability parameterization 
c-----------------------------------------------------------------------------
       sum=0 
       do ii=1,khor 
       sum=sum+iterat(ii) 
       if (iterat(ii).gt.0) then 
        cdee=(0.824+0.041*wspeed(ii))*1.e-3
        hsens(ii) = cdee*rhol(ii)*wspeed(ii)*(-tl(ii)+tec(id3sur(ii)))
        hsens(ii) = hsens(ii)*cp
	  hsens(ii)=sign( max(eps,abs(hsens(ii))),hsens(ii))


        evap(ii)  = cdee*rhol(ii)*wspeed(ii)*(-sphumi(ii)+spsur(ii))

        tl10(ii) = tl(ii) 
        q10(ii)  = sphumi(ii) 
      end if 
      enddo  
    
      if(loxstat1)then
	  write(*,*)'3 no convergence '
        write(*,*)'evap:'
        call xstat(evap,khor,zmax,zmin,zmean,zsum)
	  print*,'evap min,max,mean ',zmin,zmax,zmean
 
 
         write(*,*)'hsens:'
        call xstat(hsens,khor,zmax,zmin,zmean,zsum)
	  print*,'hsense min,max,mean ',zmin,zmax,zmean
      endif
      I=1
      IF(I .EQ. 1)STOP

c----------- windstress ---------------------------------
      do i=1,khor 
 
c================================================================== 
c       Drag coefficient Cd (cdd) for neutral stability in the atmosphere 
c================================================================== 
*1-----------------------------
c       Charnock (linear form) (1955) 
c        cdd(i) = (0.7+0.09*wspeed(i))*1.e-3 
*2----------------------------- 
c       Smith and Banke (1975) 
c        cdd(i) = (0.63+0.066*wspeed(i))*1.e-3 
*3----------------------------- 
c       Wu (1980)  
c       cdd(i) = (0.8+0.065*wspeed(i))*1.e-3 
 

*4-----------------------------
c      cdd mit Stabilitaetsansatz f�r die Atmosph�re 
c      ansatz nach J. Launianien and T. Vihma (1990) 


c         if(iterat(i).gt.0.and.wspeed(i).gt.22.)then   !Irina
         if(iterat(i).gt.0.or.wspeed(i).gt.18.)then   !Irina
c        if(iterat(i).gt.0)then                        !Irina
c------in case of no convergence, apply neutral stability parameterization, Wu (1980)
	 cdd(i) = (0.8+0.065*wspeed(i))*1.e-3   
	 else
          xlog=log(10./zro(i))-phim(i)
          xlog =sign(max(abs(xlog),eps),xlog)
c       cdd(i) = (xk/(log(10./zro(i))-phim(i)))**2.0 
       cdd(i) = (xk/xlog)**2
       endif
 
c================================================================== 
 
 
        cddu(i) = cdd(i) 


*-------- wind stress-------------------
        uw(i)  = cdd(i)*rhol(i)*(wspeed(i))*(uw(i)) 
        vw(i)  = cdd(i)*rhol(i)*(wspeed(i))*(vw(i)) 
      enddo 



      do ii=1,khor
        T_lim=(3.*abs(tl(ii))+5.)/( 1.+.1*abs(tl(ii)) )
        if(abs(tl10(ii)).gt.T_lim) then
        cdee=(0.824+0.041*wspeed(ii))*1.e-3
        hsens(ii)= cdee*rhol(ii)*wspeed(ii)*(-tl(ii)+tec(id3sur(ii)))
        hsens(ii)= hsens(ii)*cp
c	  hsens(ii)=sign( max(eps,abs(hsens(ii))),hsens(ii))
        evap(ii) = cdee*rhol(ii)*wspeed(ii)*(-sphumi(ii)+spsur(ii))
        tl10(ii) = tl(ii) 
        q10(ii)  = sphumi(ii) 
        endif
	enddo 


      if(loxstat1)then
	  write(*,*)'3 After T_lim'
        call xstat(evap,khor,zmax,zmin,zmean,zsum)
	  print*,'evap min,max,mean ',zmin,zmax,zmean
      endif
 
 
c     turbulent heat fluxes were already calculated in the loop, 
c     longwave radiation in the main  
      do i=1,khor 
        qw(i)  = -evap(i)*xlat- hsens(i) 
        fqs(i) = hsens(i)*(-1.0) 
        fql(i) = evap(i)*(-1.0)*xlat 
        q10(i) = 10.*pac(i)*q10(i)/(0.622+0.378*q10(i)) 
      enddo 
c      print*,'maximale Iteration in fluxes:',imax, 
c     +  'Anzahl Pkte. ohne Abbruch',sum  
c      write(7,*) 'max Iteration in fluxes:',imax, 
c     +  'Anzahl Pkte. ohne Abbruch',sum 
      
	
	
	return 
      end 
