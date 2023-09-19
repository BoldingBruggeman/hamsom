*trmice_new.F
*��1 druxav
* Change on line 356: du(nwet+isu) = 0.5*(ley(nwet+isu)+ley(lww(0,1,nwet+isu))) !Irina
* Change on line 555: dv(nwet+isv) = 0.5*(ley(nwet+isv)+ley(lww(1,0,nwet+isv))) !Irina
*��2 strom
*��3 konv
*��4-5deco1d3d and comp3d1d
*��6 motmit
* Change on line 2204: d(nwet+is) = 0.5*(ley(nwet+is)+ley(lww(0,1,nwet+is))) !Irina	 
* Change on line 2630: d(nwet+is) = 0.5*(ley(nwet+is)+ley(lww(1,0,nwet+is))) !Irina
*��7 sor
*��8 sorcof
*��9 konti

      subroutine druxav(vtmit,avmax,cdeav,szahl,ukonv,vkonv,nndrei)
c     -------------------------------------------------------------
c     sloving of implicit terms in the equation of motion 
c          -time averaged barotropic pressure gradients in equation of motion 
c          -implicit calculation of vertical exchange (direct elimination, gauss)
c
c     after call of druxav are transports defined for the new time level n+1   
c     Turbulence closure after Kochergin, 1987; pohlmann, 1996 and Schrum, 1997 
c     -------------------------------------------------------------
c
      include 'C_model.f'
      include 'C_index.f'

      parameter(ilop1=ilo+1)
      parameter(khor1=khor+1)
c

      dimension ukonv(nndrei),vkonv(nndrei)
      dimension szahl(ndrei)
      dimension schmiz(ilo)

      common /ind/ iwet(khor1),ldep(khor),lazc(khor),
     cindend(n),isornr(n),isorsr(n),islab(n)
c

      common uc(ndrei),vc(ndrei),stc(ndrei),avc(ndrei),z(m,n) 
      common zac(khor),wobc(khor),stuvc(khor),fricu(khor),cxc(khor) 
      common cyc(khor),pac(khor),txc(khor),tyc(khor) 
      common stpc(ndrei),sac(ndrei),tec(ndrei) 
      common pres(ilo),wc(ndrei),fricv(khor) 

c
      common /dreh/ sinfu(m),sinfv(m),cosfu(m),cosfv(m),sincx(m),
     ssincy(m),bx(m),by(m),pxu(m),pyv(m),pyu(m),pxv(m)
c
      common /gitter/ dt,r,g,dl,dlr,dlrh,dln(m),rdln(m),rad,dth,
     ddlvo(m),dlvu(m),gh,rdt,dt2,r4,coru(m),corv(m),dtrdln(m),dtdlr
c
      common/num/dc(ilo),av(ilo),ad(ilo),dh(ilo),pd(ilo),
     pprd(ilo),pr2d(ilo),r2d(ilo),tkw(ilo),tau(ilo),dd(ilo)
      dimension qa(ndrei),qbet(ndrei),qn(ndrei),rd(ndrei)
c
      common /met/ windx,windy,wlam,pnull,rhoq(ilo),refrho(ilo),stress
c
      dimension su(0:ndrei),rich(ndrei),
     1 sv(0:ndrei),rho(ilo,m,3),bvfr(ndrei),
     2 dudz(ndrei),avv(ilo,m,2),za(m,3),d(m,3),du(ndrei),dv(ndrei),
     3 dsch(ndrei)
      common /auxint/ jc(m,3),laz(m,2),lay(ilo,m,3)
c
c############ please check if this can be removed????
c!!!!!!!!!
c      common/num/dz(ilo),av(ilo),ah(ilo),dh(ilo),pd(ilo),
c     pprd(ilo),pr2d(ilo),r2d(ilo),tkw(ilo),tau(ilo),dd(ilo),
c     qqa(ilop1),qbet(ilop1),qn(ilop1),rd(ilop1)
c
      logical init
      common /vec2/ init,iisu(khor),iisv(khor),ley(0:ndrei),jc2(0:khor),
     1 lw2(-1:1,-1:1,khor),lww(-1:1,-1:1,ndrei),lev(ndrei),lind(ilo+1)
      data init /.true./
      dimension jcc(ilo,m,3),zac0(0:khor),stuvc0(0:khor)
      dimension tieo(khor),tieu(khor),lun(khor),lob(khor),mask(khor)
      dimension kl(ilo),qab(khor),dpu(khor),dpv(khor)
#ifdef MPI
      include 'C_mpi.f'
      include 'C_parallel.f'
#endif 
      mz = m-1
      nz = n-1

c-----------------------------------------------------------------      
c    epsilon minimum layer thickness   --> should be moved to main.f
c-------------------------------------------------------------------

      epsilon=1.

      if (init) then
         init = .false.
c
         do lw = 0,khor
            jc2(lw) = 0
         end do
         do j = 1,n
            do lw = lb(j),le(j)
               ldown = lazc(lw)
               if (ldown.gt.1) jc2(lw) = 1
            end do
         end do
c
         do j = 1,n
            do lw = lb(j),le(j)
               nwet = indwet(lw)
               ldown = lazc(lw)
               do k = 1,ldown-1
                  ley(nwet+k) = pd(k) 
               end do
               ley(nwet+ldown) = ldep(lw)
            end do
         end do
c
         do k = 1,ilo
            kl(k) = 0
         end do
         do j = 1,n
            do lw = lb(j),le(j)
               ldown = lazc(lw)
               do k = 1,ldown
                  kl(k) = kl(k)+1
               end do
            end do
         end do
c
         lind(1) = 0
         do k = 1,ilo
            lind(k+1) = lind(k)+kl(k)
         end do
c
         do k = 1,ilo
            kl(k) = 0
         end do
c
         do j = 1,n
            do lw = lb(j),le(j)
               nwet = indwet(lw)
               ldown = lazc(lw)
               do k = 1,ldown
                  kl(k) = kl(k)+1
                  lev(lind(k)+kl(k)) = nwet+k
               end do
            end do
         end do
c
         do kk = 1,2
            k = kk+1
            do lw = lb(kk),le(kk)
               i = iwet(lw)
               jc(i,k) = lw
            end do
         end do
         do kk = 2,nz
            do l = 2,3
               do i = 1,m
                  jc(i,l-1) = jc(i,l)
               end do
            end do
            do i = 1,m
               jc(i,3) = 0
            end do
            do lw = lb(kk+1),le(kk+1)
               i = iwet(lw)
               jc(i,3) = lw
            end do
            do lw = lb(kk),le(kk)
               i = iwet(lw)
               do jj = -1,1
                  do ii = -1,1
                     lw2(ii,jj,lw) = jc(i+ii,2+jj)
                  end do
               end do
            end do
         end do
c
         do kk = 1,2
            k = kk+1
            do i = 1,m
               do j = 1,ilo
                  jcc(j,i,k) = 0
               end do
            end do
            do lw = lb(kk),le(kk)
               i = iwet(lw)
               ldown = lazc(lw)
               nwet = indwet(lw)
               do j = 1,ldown
                  jcc(j,i,k) = nwet+j
               end do
            end do
         end do
         do kk = 2,nz
            do i = 1,m
               do j = 1,ilo
                  jcc(j,i,1) = jcc(j,i,2)
                  jcc(j,i,2) = jcc(j,i,3)
               end do
            end do
            do i = 1,m
               do j = 1,ilo
                  jcc(j,i,3) = 0
               end do
            end do
            do lw = lb(kk+1),le(kk+1)
               i = iwet(lw)
               ldown = lazc(lw)
               nwet = indwet(lw)
               do j = 1,ldown
                  jcc(j,i,3) = nwet+j
               end do
            end do
            do lw = lb(kk),le(kk)
               i = iwet(lw)
               ldown = lazc(lw)
               nwet = indwet(lw)
               do j = 1,ldown
                  do ii = -1,1
                     do jj = -1,1
                        lww(ii,jj,nwet+j) = jcc(j,i+ii,2+jj)
                     end do
                  end do
               end do
            end do
         end do
c
         do lw = 1,khor
            iisu(lw) = 0
            iisv(lw) = 0
         end do
         do kk = 2,nz
            do lw = lb(kk),le(kk)
               i = iwet(lw)
               if (i.gt.2.and.lw2(0,1,lw).ne.0) then
                  isu = min0(lazc(lw),lazc(lw2(0,1,lw)))
                  iisu(lw) = -isu
                  lole = jc2(lw)*jc2(lw2(0,1,lw))
                  if (lole.ne.0) iisu(lw) = isu 
               end if
               if (i.gt.1.and.kk.gt.2.and.lw2(1,0,lw).ne.0) then
                  isv = min0(lazc(lw),lazc(lw2(1,0,lw)))
                  iisv(lw) = -isv
                  lole = jc2(lw)*jc2(lw2(1,0,lw))
                  if (lole.ne.0) iisv(lw) = isv
               end if
            end do
         end do
      end if
c
      zac0(0) = 0.
      stuvc0(0) = 0.

#ifdef MPI
      do J=JA1,JA2
      do I = IA1,IA2
         lw=izet(I,J)
       if(lw .GT. 0) THEN
       zac0(lw) = zac(lw)
        stuvc0(lw) = stuvc(lw)
       endif
      enddo
      enddo
#else
      do lw = 1,khor
         zac0(lw) = zac(lw)
         stuvc0(lw) = stuvc(lw)
      end do
#endif
c
      su(0) = 0.
      sv(0) = 0.
c
      ley(0) = 0
      do j = 1,n
         do lw = lb(j),le(j)
            nwet = indwet(lw)
            ldown = lazc(lw)
            ley(nwet+ldown) = ldep(lw)
         end do
      end do
c---------------------------
      ape = 0.
      tke = 0.
      avmin = 0.134e-6
      epsi = 1.e-7

c--------------------------------------------------------------------
c   constants for vertical diffusion
c--------------------------------------------------------------------

      cde = cdeav
      rimax = 0.23
      stfak = 1.
      avminp = avmin+0.0001e-4
      tiemax =  100.
      tiemm= 100.
      tiemid = 100.
      gramin = 1.
      zeiko = 0.*3600.+1.e-15
      zkon = dt/zeiko
      zkon = amin1(1.,zkon)
      czkon = 1.-zkon

c------------------------------------------------------------------------
c     column loop preparaption for layer thickness
c------------------------------------------------------------------------
#ifdef MPI
      do 7 kk =ICEJA1,ICEJA2
      lwa = lb1(kk)
      lwe = le1(kk)
#else
      do 7 kk = 2,nz
      lwa = lb(kk)
      lwe = le(kk)
#endif
      if (lwa.le.lwe) then
         llb = indwet(lwa)+1
         lle = indwet(lwe)+lazc(lwe)
         do ll = llb,lle
            j = indver(ll)
            dsch(ll) = pd(j)
         end do
      end if

      do lw = lwa,lwe 
         nwet = indwet(lw)
         ldown = lazc(lw)
         dsch(nwet+ldown) = real(ldep(lw))
         dsch(nwet+1) = dsch(nwet+1)+zac(lw)
      end do
c--------------------------------------------------------------------------
    7 continue

c--------------------------------------------------------------------------
c    loop for u-component
c--------------------------------------------------------------------------

#ifdef MPI
      do kk = ICEJ1,ICEJ2 
      lwa = lb0(kk)
      lwe = le0(kk)
#else
      do kk = 2,nz 
      lwa = lb(kk)
      lwe = le(kk)
#endif
      if (lwa.le.lwe) then
         llb = indwet(lwa)+1
         lle = indwet(lwe)+lazc(lwe)
         do ll = llb,lle
            j = ll-indwet(llw(ll))
            lw = llw(ll)
            if (iisu(lw).ne.0) then
               isu = abs (iisu(lw))
               if (j.le.isu) then
                  du(ll) = pd(j)
               end if
            end if
         end do
      end if
      end do

#ifdef MPI
      do kk =ICEJ1,ICEJ2 
      do lw = lb0(kk),le0(kk)
#else
      do kk = 2,nz
      lwa = lb(kk)
      lwe = le(kk)
      do lw = lwa,lwe
#endif
      if (iisu(lw).ne.0) then
      isu = abs (iisu(lw))
      lole = jc2(lw)*jc2(lw2(0,1,lw))
      dreid = lole
      zweid = 1.-dreid
      nwet = indwet(lw)
      i = iwet(lw)
      du(nwet+isu) = 0.5*(ley(nwet+isu)+ley(lww(0,1,nwet+isu))) !Irina
      du(nwet+1) = du(nwet+1)+0.5*(zac0(lw)+zac0(lw2(0,1,lw)))  !Irina
c      du(nwet+isu) = amin0(ley(nwet+isu),ley(lww(0,1,nwet+isu)))
c      du(nwet+1) = du(nwet+1)+amin1(zac0(lw),zac0(lw2(0,1,lw)))

      gne=4./(min(1,lw2(-1,0,lw))+min(1,lw2(-1,1,lw))
     1 +min(1,lw2(0,0,lw))+min(1,lw2(0,1,lw)))
      gse=4./(min(1,lw2(1,0,lw))+min(1,lw2(1,1,lw))
     1 +min(1,lw2(0,0,lw))+min(1,lw2(0,1,lw)))
      dpu(lw) = dprod(pxu(i),stuvc0(lw2(0,1,lw)))
     1 -dprod(pxu(i),stuvc0(lw))
     2 +dprod(pyu(i),gne*(stuvc0(lw2(-1,1,lw))+stuvc0(lw2(-1,0,lw))
     3 +stuvc0(lw)+stuvc0(lw2(0,1,lw)))
     4 -gse*(stuvc0(lw)+stuvc0(lw2(1,0,lw))+stuvc0(lw2(0,1,lw))
     5 +stuvc0(lw2(1,1,lw))))
      dpu(lw) = dpu(lw)*2.
      rimp = fricu(lw)
      du(nwet+1)=max(epsilon,du(nwet+1))
      uc(nwet+1) = (zweid*rimp+dreid)*(uc(nwet+1)-dpu(lw)*du(nwet+1)
     //(stc(nwet+1)+2.*rhoq(1)+stc(lww(0,1,nwet+1))))
      su(nwet+1) = uc(nwet+1)/du(nwet+1)   ! currents
      end if
         end do
      end do

#ifdef MPI
      do kk = ICEJ1,ICEJ2 
      lwa = lb0(kk)
      lwe = le0(kk)
#else
      do kk = 2,nz 
      lwa = lb(kk)
      lwe = le(kk)
#endif
      if (lwa.le.lwe) then
         llb = indwet(lwa)+1
         lle = indwet(lwe)+lazc(lwe)
         do ll = llb,lle
            j = ll-indwet(llw(ll))
            lw = llw(ll)
      if (iisu(lw).gt.0) then
      isu = iisu(lw) 
      if (j.eq.1) then
         qa(ll) = 0.0
         rd(ll) = 1./du(ll)
      else if (j.le.isu) then
         uc(ll) = uc(ll)-dpu(lw)*du(ll)/(stc(ll)+2.*rhoq(j)
     1    +stc(lww(0,1,ll)))
         rd(ll) = prd(j)
         qa(ll) = vtmit*dt*(avc(ll)+avc(lww(0,1,ll)))/(du(ll-1)+du(ll))
      end if
      end if

         end do
      end if
      end do

#ifdef MPI
      do kk = ICEJ1,ICEJ2
      do lw = lb0(kk),le0(kk)
#else
      do kk = 2,nz
      lwa = lb(kk)
      lwe = le(kk)
      do lw = lwa,lwe
#endif
         if (iisu(lw).gt.0) then
            isu = iisu(lw)
            nwet = indwet(lw)
            rimp = fricu(lw)
            qa(nwet+isu) = rimp*qa(nwet+isu)
            uc(nwet+isu) = uc(nwet+isu)*rimp
            rd(nwet+isu) = 1./du(nwet+isu)
         end if
      end do
      end do
c
c backward sweep for u-component
c
#ifdef MPI
      do kk = ICEJ1,ICEJ2
      do lw = lb0(kk),le0(kk)
#else
      do kk = 2,nz
      do lw = lb(kk),le(kk)
#endif
         if (iisu(lw).gt.0) then
            isu = iisu(lw)
            nwet = indwet(lw)
            qn(nwet+isu) = uc(nwet+isu)
            qbet(nwet+isu) = 1./(1.0+rd(nwet+isu)*qa(nwet+isu))
            qab(lw) = qa(nwet+isu)*qbet(nwet+isu)*rd(nwet+isu)
         end if
      end do
      end do

#ifndef MPI
      do j = ilo,1,-1
!CDIR VECTOR NODEP
         do ih = lind(j)+1,lind(j+1)
            ll = lev(ih)
            lw = llw(ll)
            if (iisu(lw).gt.0) then
               isu = iisu(lw)
               if (j.lt.isu) then
                  qn(ll) = uc(ll)+qn(ll+1)*qab(lw)
                  qbet(ll) = 1./(1.0+rd(ll)*(qa(ll)+qa(ll+1)
     1             *(1.-qab(lw))))
                  qab(lw) = qa(ll)*qbet(ll)*rd(ll)
               end if
            end if
         end do
      end do
#else
        DO J=J1,J2
        DO lw=lb0(j),le0(j)
         nwet=indwet(lw)
         isu = iisu(lw)
         if (isu.gt.0) then
         do l=ilo,1,-1
         if(l.lt.isu) then
         ll=nwet+l
         qn(ll) = uc(ll)+qn(ll+1)*qab(lw)
         qbet(ll) = 1./(1.0+rd(ll)*(qa(ll)+qa(ll+1)
     1   *(1.-qab(lw))))
         qab(lw) = qa(ll)*qbet(ll)*rd(ll)
         endif
         ENDDO
         endif
         ENDDO
         ENDDO
#endif
c
c forward sweep for u-component
c
#ifdef MPI
      do kk = ICEJ1,ICEJ2      
      do lw = lb0(kk),le0(kk)
#else
      do kk = 2,nz      
      do lw = lb(kk),le(kk)
#endif
         if (iisu(lw).ge.2) then
            nwet = indwet(lw)
            uc(nwet+1) = qbet(nwet+1)*qn(nwet+1)
            su(nwet+1) = uc(nwet+1)*rd(nwet+1)
         end if
      end do
      end do
#ifndef MPI
      do j = 2,ilo
!CDIR VECTOR NODEP
         do ih = lind(j)+1,lind(j+1)
            ll = lev(ih)
            lw = llw(ll)
            if (iisu(lw).gt.0) then
               isu = iisu(lw)
               if (j.le.isu) then
                  uc(ll) = (qn(ll)+qa(ll)*su(ll-1))*qbet(ll)
                  su(ll) = uc(ll)*rd(ll)
               end if
            end if
         end do
      end do
#else
        DO J=J1,J2
        DO lw=lb0(j),le0(j)
         nwet=indwet(lw)
         isu = iisu(lw)
         if (isu.gt.0) then
         do l=2,ilo
         if(l.le.isu) then
         ll=nwet+l
         uc(ll) = (qn(ll)+qa(ll)*su(ll-1))*qbet(ll)
         su(ll) = uc(ll)*rd(ll)
         endif
         ENDDO
         endif
         ENDDO
         ENDDO
#endif

#ifdef MPI
      call INTI_COMM_VAR
      call com_horizontal(su(1),ndrei,dtemp,iel,0,1,1,3,1,icdepb0,
     *idepb0,icdept0,
     *idept0,ncsendb0,nsendb0,ncsendt0,nsendt0,nbhlayer0,nthlayer0)
#endif

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     loops for v-component
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c

#ifdef MPI
      do kk = icej1,icej2 
      lwa = lb0(kk)
      lwe = le0(kk)
#else
      do kk = 2,nz 
      lwa = lb(kk)
      lwe = le(kk)
#endif
      if (lwa.le.lwe) then
         llb = indwet(lwa)+1
         lle = indwet(lwe)+lazc(lwe)
         do ll = llb,lle
            j = ll-indwet(llw(ll))
            lw = llw(ll)
            if (iisv(lw).ne.0) then
               isv = abs (iisv(lw))
               if (j.le.isv) then
                  dv(ll) = pd(j)
               end if
            end if
         end do
      end if
      end do

#ifdef MPI
      do kk = ICEJ1,ICEJ2
      do lw = lb0(kk),le0(kk)
#else
      do kk = 2,nz
      do lw = lb(kk),le(kk)
#endif
      if (iisv(lw).ne.0) then
      isv = abs (iisv(lw))
      lole = jc2(lw)*jc2(lw2(1,0,lw))
      dreid = lole
      zweid = 1.-dreid
      nwet = indwet(lw)
      i = iwet(lw)
      dv(nwet+isv) = 0.5*(ley(nwet+isv)+ley(lww(1,0,nwet+isv)))  !Irina
      dv(nwet+1) = dv(nwet+1)+0.5*(zac(lw)+zac(lw2(1,0,lw)))     !Irina
c      dv(nwet+isv) = amin0(ley(nwet+isv),ley(lww(1,0,nwet+isv)))
c      dv(nwet+1) = dv(nwet+1)+amin1(zac(lw),zac(lw2(1,0,lw)))

      gsw=4./(min(1,lw2(1,0,lw))+min(1,lw2(1,-1,lw))+min(1,lw2(0,0,lw))
     1 +min(1,lw2(0,-1,lw)))
      gse=4./(min(1,lw2(1,0,lw))+min(1,lw2(1,1,lw))
     1 +min(1,lw2(0,0,lw))+min(1,lw2(0,1,lw)))
      dpv(lw) = dprod(pyv(i),stuvc0(lw))
     1 -dprod(pyv(i),stuvc0(lw2(1,0,lw)))
     2 +dprod(pxv(i),gse*(stuvc0(lw2(0,1,lw))+stuvc0(lw)
     3 +stuvc0(lw2(1,1,lw))+stuvc0(lw2(1,0,lw)))
     4 -gsw*(stuvc0(lw2(0,-1,lw))+stuvc0(lw)+stuvc0(lw2(1,-1,lw))
     5 +stuvc0(lw2(1,0,lw))))
      dpv(lw) = dpv(lw)*2.
      rimp = fricv(lw)
      dv(nwet+1)=max(epsilon,dv(nwet+1))
      vc(nwet+1) = (rimp*zweid+dreid)*(vc(nwet+1)-dpv(lw)*dv(nwet+1)
     //(stc(nwet+1)+2.*rhoq(1)+stc(lww(1,0,nwet+1))))
      sv(nwet+1) = vc(nwet+1)/dv(nwet+1)    ! current speed
      end if
         end do
      end do


#ifdef MPI
      do kk = ICEJ1,ICEJ2 
      lwa = lb0(kk)
      lwe = le0(kk)
#else
      do kk = 2,nz 
      lwa = lb(kk)
      lwe = le(kk)
#endif
      if (lwa.le.lwe) then
         llb = indwet(lwa)+1
         lle = indwet(lwe)+lazc(lwe)
         do ll = llb,lle
            j = ll-indwet(llw(ll))
            lw = llw(ll)
      if (iisv(lw).gt.0) then
      isv = iisv(lw) 
      if (j.eq.1) then
      rd(ll) = 1./dv(ll)
      qa(ll) = 0.0
      else if (j.le.isv) then
      vc(ll) = vc(ll)-dpv(lw)*dv(ll)/(stc(ll)+2.*rhoq(j)
     1 +stc(lww(1,0,ll)))
      rd(ll) = prd(j)
   14 qa(ll) = vtmit*dt*(avc(ll)+avc(lww(1,0,ll)))
     1 /(dv(ll-1)+dv(ll))
      end if
      end if

         end do
      end if
      end do

#ifdef MPI
      do kk = ICEJ1,ICEJ2
         do lw = lb0(kk),le0(kk)
#else
      do kk = 2,nz
         do lw = lb(kk),le(kk)
#endif
            if (iisv(lw).gt.0) then
               isv = iisv(lw)
               nwet = indwet(lw)
               rimp = fricv(lw)
               qa(nwet+isv) = rimp*qa(nwet+isv)
               vc(nwet+isv) = vc(nwet+isv)*rimp
               rd(nwet+isv) = 1./dv(nwet+isv)
            end if
         end do
      end do
c
c backword sweep for v-component
c
#ifdef MPI
      do kk = ICEJ1,ICEJ2
      do lw = lb0(kk),le0(kk)
#else
      do kk = 2,nz
      do lw = lb(kk),le(kk)
#endif
         if (iisv(lw).gt.0) then
            isv = iisv(lw)
            nwet = indwet(lw)
            qn(nwet+isv) = vc(nwet+isv)
            qbet(nwet+isv) = 1./(1.0+rd(nwet+isv)*qa(nwet+isv))
            qab(lw) = qa(nwet+isv)*qbet(nwet+isv)*rd(nwet+isv)
         end if
      end do
      end do
#ifndef MPI
      do j = ilo,1,-1
!CDIR VECTOR NODEP
         do ih = lind(j)+1,lind(j+1)
            ll = lev(ih)
            lw = llw(ll)
            if (iisv(lw).gt.0) then
               isv = iisv(lw)
               if (j.lt.isv) then
                  qn(ll) = vc(ll)+qn(ll+1)*qab(lw)
                  qbet(ll) = 1./(1.0+rd(ll)*(qa(ll)+qa(ll+1)
     1             *(1.-qab(lw))))
                  qab(lw) = qa(ll)*qbet(ll)*rd(ll)
               end if
            end if
         end do
      end do
#else
        DO J=J1,J2
        DO lw=lb0(j),le0(j)
         nwet=indwet(lw)
         isv = iisv(lw)
         if (isv.gt.0) then
         do l=ilo,1,-1
         if(l.lt.isv) then
         ll=nwet+l
                  qn(ll) = vc(ll)+qn(ll+1)*qab(lw)
                  qbet(ll) = 1./(1.0+rd(ll)*(qa(ll)+qa(ll+1)
     1             *(1.-qab(lw))))
                  qab(lw) = qa(ll)*qbet(ll)*rd(ll)         
         endif
         ENDDO
         endif
         ENDDO
         ENDDO
#endif
c
c forward sweep for v-component
c
#ifdef MPI
      do kk = ICEJ1,ICEJ2      
      do lw = lb0(kk),le0(kk)
#else
      do kk = 2,nz      
      do lw = lb(kk),le(kk)
#endif
         if (iisv(lw).ge.2) then
            nwet = indwet(lw)
            vc(nwet+1) = qn(nwet+1)*qbet(nwet+1)
            sv(nwet+1) = vc(nwet+1)*rd(nwet+1)
         end if
      end do
      end do

#ifndef MPI
      do j = 2,ilo
!CDIR VECTOR NODEP
         do ih = lind(j)+1,lind(j+1)
            ll = lev(ih)
            lw = llw(ll)
            if (iisv(lw).gt.0) then
               isv = iisv(lw)
               if (j.le.isv) then
                  vc(ll) = (qn(ll)+qa(ll)*sv(ll-1))
     1             *qbet(ll)
                  sv(ll) = vc(ll)*rd(ll)
               end if
            end if
         end do
      end do
#else
        DO J=J1,J2
        DO lw=lb0(j),le0(j)
         nwet=indwet(lw)
         isv = iisv(lw)
         if (isv.gt.0) then
         do l=2,ilo
         if(l.le.isv) then
         ll=nwet+l
                  vc(ll) = (qn(ll)+qa(ll)*sv(ll-1))
     1             *qbet(ll)
                  sv(ll) = vc(ll)*rd(ll)
         endif
         ENDDO
         endif
         ENDDO
         ENDDO
#endif

c inter process communication
#ifdef MPI      
      call com_vertical(sv(1),ndrei,dtemp,iel,1,0,1,3,1,lazc)
      call MPI_Waitall(nreq1, ir_r_b, status, ierr)
      INUM3=0
      INUM4=0
      call HOR_UNPACK(su(1),ndrei,dtemp,iel,0,1,nbhlayer0,nthlayer0,
     &idepb0,idept0,icdepb0,icdept0,1,3,iwet)
      INUM1=0
      INUM2=0
      call VER_UNPACK(sv(1),ndrei,dtemp,iel,1,0,1,3,lazc)
#endif
c-------------------------------------------------------------
c
c     update uc,vc and estimate the time-space dependent 
c       turbulent vertical exchange coefficient 
c-------------------------------------------------------------

#ifdef MPI
      do kk = icej1,icej2 
      lwa = lb0(kk)
      lwe = le0(kk)
#else
      do kk = 2,nz 
      lwa = lb(kk)
      lwe = le(kk)
#endif
      if (lwa.le.lwe) then
         llb = indwet(lwa)+1
         lle = indwet(lwe)+lazc(lwe)
         do ll = llb,lle
            j = ll-indwet(llw(ll))

      if (j.eq.1) then
         szahl(ll) = 0.0
      else

      dschm = max(epsilon,0.5*(dsch(ll-1)+dsch(ll)))
      shu = 0.5*(su(lww(0,-1,ll-1))+su(ll-1)
     1 -su(lww(0,-1,ll))-su(ll))
      shv = 0.5*(sv(ll-1)+sv(lww(-1,0,ll-1))
     1 -sv(ll)-sv(lww(-1,0,ll)))
      bvfr(ll) = g*2.*
     * (stpc(ll)+rhoq(j)-stc(ll-1)-rhoq(j-1))/dschm
     //(rhoq(j-1)+stc(ll-1)+rhoq(j)+stc(ll))
      dudz(ll) =  (shu**2+shv**2)/dschm**2
      if (bvfr(ll) .lt. 0.) bvfr(ll) = 15.*bvfr(ll)
      rich(ll) = bvfr(ll)/(dudz(ll)+1.e-15)
      riri = rich(ll)
      riri = amax1(3.e-3,riri)
      riri= amin1(riri,1.e6)
      rfsq = riri**2-0.316*riri+0.0346
      rfsq = amax1(0.,rfsq)
      rrff = 0.725*(riri+0.186-sqrt(rfsq))
      szahl(ll) =rrff/riri

      end if
         end do
      end if
      end do

#ifdef MPI
      do kk = icej1,icej2
         do lw = lb0(kk),le0(kk)
#else
      do kk = 2,nz
         do lw = lb(kk),le(kk)
#endif
            nwet = indwet(lw)
            ldown = lazc(lw)
            tieo(lw) = dsch(nwet+1)
            lob(lw) = ldown
            mask(lw) = 0
         end do
      end do
c
#ifndef MPI
      do j=2,ilo
!CDIR VECTOR NODEP
         do ih = lind(j)+1,lind(j+1)
            ll = lev(ih)
            lw = llw(ll)
      nwet = indwet(lw)
      ldown = lazc(lw)
      if (mask(lw).eq.0.and.j.le.ldown) then
      if (rich(nwet+j).ge.rimax) then
         lob(lw) = j-1
         mask(lw) = 1
      else
         tieo(lw) = tieo(lw)+dsch(nwet+j)
      end if
      end if
         end do
      end do
#else
        DO J=J1,J2
        DO lw=lb0(j),le0(j)
         nwet=indwet(lw)
         ldown=lazc(lw)
         do l=2,ilo
         if(mask(lw).eq.0 .and. l .le.ldown) then
         if (rich(nwet+l).ge.rimax) then
         lob(lw) = l-1
         mask(lw) = 1
      else
         tieo(lw) = tieo(lw)+dsch(nwet+l)
      end if
      end if
         ENDDO
         ENDDO
         ENDDO
#endif

#ifdef MPI
      do kk = icej1,icej2
         do lw = lb0(kk),le0(kk)
#else
      do kk = 2,nz
         do lw = lb(kk),le(kk)
#endif
            nwet = indwet(lw)
            ldown = lazc(lw)
            tieo(lw) = amin1(tieo(lw),tiemax)
            tieu(lw) = dsch(nwet+ldown)
            lun(lw) = ldown+1
            mask(lw) = 0
         end do
      end do
#ifndef MPI
      do j = ilo,1,-1
!CDIR VECTOR NODEP
         do ih = lind(j)+1,lind(j+1)
            ll = lev(ih)
            lw = llw(ll)
            if (mask(lw).eq.0.and.j.gt.lob(lw)) then
               if (rich(ll).ge.rimax) then
                  mask(lw) = 1
                  lun(lw) = j+1
               else
                  tieu(lw) = tieu(lw)+dsch(ll-1)
               end if
            end if
         end do
      end do
#else
       DO J=J1,J2
        DO lw=lb0(j),le0(j)
         nwet=indwet(lw)
         ldown=lazc(lw)
         DO l=ldown,1,-1
         ll=nwet+l
         if (mask(lw).eq.0.and.l.gt.lob(lw)) then
               if (rich(ll).ge.rimax) then
                  mask(lw) = 1
                  lun(lw) = l+1
               else
                  tieu(lw) = tieu(lw)+dsch(ll-1)
               end if
           end if
         ENDDO
         ENDDO
         ENDDO
#endif
c
#ifdef MPI
      do kk = icej1,icej2
      do lw = lb0(kk),le0(kk)
#else
      do kk = 2,nz
      do lw = lb(kk),le(kk)
#endif
      ldown = lazc(lw)
      tieu(lw)=amin1(tieu(lw),tiemid)
      if (lob(lw) .eq. ldown) lun(lw) = ldown+1
      end do
      end do

#ifdef MPI
      do kk  = icej1,icej2
      lwa = lb0(kk)
      lwe = le0(kk)
#else
      do kk  = 2,nz
      lwa = lb(kk)
      lwe = le(kk)
#endif
      if (lwa.le.lwe) then
         llb = indwet(lwa)+1
         lle = indwet(lwe)+lazc(lwe)
         do ll = llb,lle
            lw = llw(ll)
            j = ll-indwet(lw)

      if (j.ne.1) then
         if (j.le.lob(lw)) then
            ctioq = (cde*tieo(lw))**2
            radi = amax1(0.,dudz(ll)-szahl(ll)*bvfr(ll))
            avc(ll) = ctioq*sqrt(radi)
         else if (j.lt.lun(lw)) then
            tieum=dsch(ll)
            tiemix=min(tieum,tiemid)
            ctimq = (cde*tiemix)**2
            radi = amax1(0.,dudz(ll)-szahl(ll)*bvfr(ll))
            avc(ll) = ctimq*sqrt(radi)
         else
            ctiuq = (cde*tieu(lw))**2
            radi = amax1(0.,dudz(ll)-szahl(ll)*bvfr(ll))
            avc(ll) = ctiuq*sqrt(radi)
         end if
      end if
         end do
      end if
      end do

#ifdef MPI
      do kk = icej1,icej2 
      lwa = lb0(kk)
      lwe = le0(kk)
#else
      do kk= 2,nz
      lwa = lb(kk)
      lwe = le(kk)
#endif
      if (lwa.le.lwe) then
         llb = indwet(lwa)+1
         lle = indwet(lwe)+lazc(lwe)
         do ll = llb,lle
            j = ll-indwet(llw(ll))
            if (j.eq.1) then
               avc(ll) = 0.
            else
               avc(ll) = amin1(avmax,avc(ll))
               avc(ll) = amax1(avmin,avc(ll))
            end if
         end do
      end if
      end do
c
      return
      end


*��3 konv
c-----------------------------------------------------------------------
      subroutine konv(dt,dz,dx,dy,zz,dz1,dzbod,
     & ukonv,vkonv,uc,vc)
c-----------------------------------------------------------------------
c     calculation of convective terms after arakawa und lamb(1977), J7 scheme
c     temporal discretisation after matsuno (2-steps scheme) 
c     after  Mesinger and Arakawa (1976)
c-----------------------------------------------------------------------
      include 'C_model.f'
      parameter(khor1=khor+1)
      dimension uc(ndrei),vc(ndrei),dx(m,n),dy(m,n)
      dimension ukonv(ndrei),vkonv(ndrei) 
      dimension dz(ilo),zz(m,n),dzbod(m,n)
      dimension ua(m,n,ilo),va(m,n,ilo),ddz(ilo) 
      dimension delu(m,n,ilo),delv(m,n,ilo)
      real fst(m,n,ilo),gst(m,n,ilo) 
      real ug(m,n,ilo),vg(m,n,ilo) 
#ifndef MPI 
      common /geo1/ ltief(m,n),izet(m,n)
#endif
      dimension schmu(m,n,ilo),schmv(m,n,ilo)
#ifdef MPI
      include 'C_index.f'
      include 'C_mpi.f'
      include 'C_parallel.f'
      common /ind/ iwet(khor1),ldep(khor),lazc(khor),indend(n),
     *             isornr(n),isorsr(n),islab(n)
#endif
c------------------------------------------------------------------------------
c    epsilon minimum layer depth ---###### should go to main
c------------------------------------------------------------------------------

      epsilon=1.

        epsil=1.e-7
       dth=dt*0.5
      zweidr=1./12. 
      sechsl=1./24.
      ninter=33
      do 399 i=2,ilo
        ddz(i)=dz(i)   !pd!
 399  continue
      ddz(1)=dz1       !pd!

      call deco1d3d(ua,uc,ndrei)
      call deco1d3d(va,vc,ndrei)

c------------------------------------------------------------------------
c     fst    = f* array after arakawa und lamb (1977) 
c     gst    = g* array after arakawa und lamb (1977) 
c     computation of f* and g* arrays after 
c     arakawa und lamb (1977) 
c-----------------------------------------------------------------------
#ifdef MPI
      do 1999 j=JB1,JB2
      do 1999 i=IB1,IB2
#else
      do 1999 i=1,m
      do 1999 j=1,n
#endif
      do 1999 k=1,ilo
          fst(i,j,k)=0.
          gst(i,j,k)=0.
          delu(i,j,k)=0.
          delv(i,j,k)=0.
          ug(i,j,k)=0.
          vg(i,j,k)=0.
 1999 continue
c
#ifdef MPI
      do 2000 j=JA1,JA2 
      do 2000 i=IA1,IA2
#else
      do 2000 j=1,n 
      do 2000 i=1,m
#endif
         jp1=min(j+1,n)
         ip1=min(i+1,m)
	  jm1=max(j-1,1)
        im1=max(i-1,1) 

      do 2200 k=1,ltief(i,j)
      fst (i,j,k) = 0.5*(ua(i,jm1,k) + ua(i,j  ,k)) 
      gst (i,j,k) = 0.5*(va(im1,j,k) + va(i,j,k))
       schmu(i,j,k)=dz(k)
       schmv(i,j,k)=dz(k)
 2200 continue 
c
      if (ltief(i,j).gt.0) then
         k=ltief(i,j)
c
         if(ltief(i,j).eq.ltief(i,jp1)) then
            schmu(i,j,k)=(dzbod(i,j)+dzbod(i,jp1))/2.
            else if (ltief(i,j).lt.ltief(i,jp1)) then 
            schmu(i,j,k)=(dzbod(i,j)+ddz(k))/2.
         else
            schmu(i,j,k)=dzbod(i,j)/2.
         end if
c
         if(ltief(i,j).eq.ltief(ip1,j)) then
            schmv(i,j,k)=(dzbod(i,j)+dzbod(ip1,j))/2.
         else if(ltief(i,j).lt.ltief(ip1,j)) then
            schmv(i,j,k)=(dzbod(i,j)+ddz(k))/2.
         else
            schmv(i,j,k)=dzbod(i,j)/2.
         end if 
c 
         schmu(i,j,1)=schmu(i,j,1)+(zz(i,j)+zz(i,jp1))/2.
         schmv(i,j,1)=schmv(i,j,1)+(zz(i,j)+zz(ip1,j))/2.
       end if

 2000 continue
#ifdef MPI
      do 11 j=JA1,JA2
      do 11 i=IA1,IA2
#else      
      do 11 i=1,m
      do 11 j=1,n
#endif
      do 11 k=1,ltief(i,j)
	schmu(i,j,k)=max(epsilon,schmu(i,j,k))
	schmv(i,j,k)=max(epsilon,schmv(i,j,k))
  11   continue


      do 2211 izeit=1,2
#ifdef MPI
      do 1 j=JA1,JA2
      do 1 i=IA1,IA2
#else 
      do 1 i=1,m
      do 1 j=1,n
#endif
      do 1 k=1,ltief(i,j)
         ug(i,j,k)=ua(i,j,k)+delu(i,j,k)
         ug(i,j,k)=ug(i,j,k)/schmu(i,j,k) 
         vg(i,j,k)=va(i,j,k)+delv(i,j,k)
         vg(i,j,k)=vg(i,j,k)/schmv(i,j,k)
  1   continue

c--------------------------------------------------------------------------------------
c     arakawa j-7 for u-equation 
c-------------------------------------------------------------------------------------
      do 10 k=1,ilo
#ifdef MPI
      do 10 j=icej1,icej2
      do 10 i=icei1,icei2
#else
      do 10 i=2,m-1 
      do 10 j=2,n-1 
#endif 
      if (k .gt. ltief(i,j)) goto 10
c      
      dl2=dx(i,j)*dx(i,j)+dy(i,j)*dy(i,j) 
      wuzdl=-dth/sqrt(dl2) 
      dxr=-dth/dx(i,j) 
      dyr=-dth/dy(i,j) 
c 
      fuuo = zweidr * (ug(i,j,k) + ug(i,j+1,k)) 
     &       * (fst(i-1,j+1,k) + 2.*fst(i,j+1,k) + fst(i+1,j+1,k)) 
      fuuw = zweidr * (ug(i,j,k) + ug(i,j-1,k)) 
     &       * (fst(i-1,j,k) + 2.*fst(i,j,k) + fst(i+1,j,k)) 
      guun = zweidr * (ug(i,j,k) + ug(i-1,j,k)) 
     &       * (gst(i,j  ,k) + gst(i-1,j  ,k) 
     &        + gst(i,j+1,k) + gst(i-1,j+1,k)) 
      guus = zweidr * (ug(i,j,k) + ug(i+1,j,k)) 
     &       * (gst(i,j  ,k) + gst(i+1,j  ,k) 
     &        + gst(i,j+1,k) + gst(i+1,j+1,k)) 
      fsuno = sechsl * (ug(i,j,k) + ug(i-1,j+1,k)) 
     &       * (gst(i,j+1,k) + gst(i-1,j+1,k) 
     &        + fst(i,j+1,k) + fst(i-1,j+1,k)) 
      fsusw = sechsl * (ug(i,j,k) + ug(i+1,j-1,k)) 
     &       * (gst(i,j,k) + gst(i+1,j,k) 
     &        + fst(i,j,k) + fst(i+1,j,k)) 
      gsunw = sechsl* (ug(i,j,k) + ug(i-1,j-1,k)) 
     &       * (gst(i,j,k) + gst(i-1,j,k) 
     &        - fst(i,j,k) - fst(i-1,j,k)) 
      gsuso = sechsl * (ug(i,j,k) + ug(i+1,j+1,k)) 
     &       * (gst(i,j+1,k) + gst(i+1,j+1,k) 
     &        - fst(i,j+1,k) - fst(i+1,j+1,k)) 
c 
      delu(i,j,k)=   dxr*(fuuo-fuuw) + dyr*(guun-guus) 
     &      + wuzdl*(fsuno-fsusw+gsunw-gsuso) 
c 

c------------------------------------------------------------------------------- 
c     arakawa j-7 for v-equation 
c-------------------------------------------------------------------------------
c 
      fvvo = zweidr * (vg(i,j,k) + vg(i,j+1,k)) 
     &      * (fst(i,j  ,k) + fst(i,j+1  ,k) 
     &       + fst(i+1,j,k) + fst(i+1,j+1,k)) 
      fvvw = zweidr * (vg(i,j,k) + vg(i,j-1,k)) 
     &      * (fst(i,j  ,k) + fst(i,j-1  ,k) 
     &       + fst(i+1,j,k) + fst(i+1,j-1,k)) 
      gvvn = zweidr * (vg(i,j,k) + vg(i-1,j,k)) 
     &      * (gst(i,j-1,k) + 2.*gst(i,j,k) + gst(i,j+1,k)) 
      gvvs = zweidr * (vg(i,j,k) + vg(i+1,j,k)) 
     &      * (gst(i+1,j-1,k) + 2.*gst(i+1,j,k) + gst(i+1,j+1,k)) 
      fsvno = sechsl * (vg(i,j,k) + vg(i-1,j+1,k)) 
     &       * (gst(i,j,k) + gst(i,j+1,k) 
     &        + fst(i,j,k) + fst(i,j+1,k)) 
      fsvsw = sechsl * (vg(i,j,k) + vg(i+1,j-1,k)) 
     &       * (gst(i+1,j,k) + gst(i+1,j-1,k) 
     &        + fst(i+1,j,k) + fst(i+1,j-1,k)) 
      gsvnw = sechsl * (vg(i,j,k) + vg(i-1,j-1,k)) 
     &       * (gst(i,j,k) + gst(i,j-1,k) 
     &        - fst(i,j,k) - fst(i,j-1,k)) 
      gsvso = sechsl * (vg(i,j,k) + vg(i+1,j+1,k)) 
     &       * (gst(i+1,j+1,k) + gst(i+1,j  ,k) 
     &        - fst(i+1,j  ,k) - fst(i+1,j+1,k)) 
c
      delv(i,j,k) =  dxr * (fvvo - fvvw) + dyr * (gvvn - gvvs) 
     &      + wuzdl * (fsvno - fsvsw + gsvnw - gsvso) 
c 
 10   continue
#ifdef MPI
      IF(izeit .EQ. 1) THEN
      call INTI_COMM_VAR

      call com_vertical(ctemp,iel1,delu,ilo,1,1,0,3,0,lazc)
      call com_vertical(ctemp,iel1,delv,ilo,1,1,0,3,1,lazc)

      call com_horizontal(ctemp,iel1,delu,ilo,1,1,0,3,0,icdepb1,
     *idepb1,icdept1
     *,idept1,ncsendb1,nsendb1,ncsendt1,nsendt1,nbhlayer1,nthlayer1)
      call com_horizontal(ctemp,iel1,delv,ilo,1,1,0,3,1,icdepb1,
     *idepb1,icdept1
     *,idept1,ncsendb1,nsendb1,ncsendt1,nsendt1,nbhlayer1,nthlayer1)
      
      call MPI_Waitall(nreq1, ir_r_b, status, ierr)

      INUM1=0
      INUM2=0

      call VER_UNPACK(ctemp,iel1,delu,ilo,1,1,0,3,lazc)
      call VER_UNPACK(ctemp,iel1,delv,ilo,1,1,0,3,lazc)

      INUM3=0
      INUM4=0
      call HOR_UNPACK(ctemp,iel1,delu,ilo,1,1,nbhlayer1,nthlayer1,
     &idepb1,idept1,icdepb1,icdept1,0,3,iwet)
      call HOR_UNPACK(ctemp,iel1,delv,ilo,1,1,nbhlayer1,nthlayer1,
     &idepb1,idept1,icdepb1,icdept1,0,3,iwet)
      ENDIF
#endif
 2211 continue
      call comp3d1d(delu,ukonv)
      call comp3d1d(delv,vkonv)
      return 
      end


*��4-5 deco1d3d and comp3d1d
c-----------------------------------------------------------------------
      subroutine deco1d3d (u,ucomp,ntot) 
c-----------------------------------------------------------------------
c     expanding 1-d-array ->  3-d-arry 
c-----------------------------------------------------------------------
      include 'C_model.f'
      include 'C_index.f'
      parameter (khor1=khor+1) 

      dimension u(m,n,ilo) 
      dimension ucomp(ntot) 
      common /ind/ iwet(khor1),ldep(khor),lazc(khor), 
     * indend(n),isornr(n),isorsr(n),islab(n) 
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
            do j=1,ilo 
               u(i,k,j) = 0.0 
            enddo
         enddo
      enddo
#ifdef MPI
      do k=JB1,JB2 
         lwa = lb2(k)
         lwe = le2(k)
#else
      do k=1,n 
         lwa = lb(k)
         lwe = le(k)
#endif
         if (lwa.le.lwe) then
            llb = indwet(lwa)+1
            lle = indwet(lwe)+lazc(lwe)
            do ll=llb,lle
               j = indver(ll)
               i = indi(ll)
               u(i,k,j) = ucomp(ll) 
            end do
         end if
      end do
c
      return 
      end

c-----------------------------------------------------------------------
      subroutine comp3d1d (d3,d1) 
c-----------------------------------------------------------------------
c    compressing 3-d arrays to compressed 1-d arrays 
c-----------------------------------------------------------------------

      include 'C_model.f'
      include 'C_index.f'
      parameter (khor1=khor+1) 

      dimension d3(m,n,ilo),d1(ndrei) 
      common /ind/ iwet(khor1),ldep(khor),lazc(khor), 
     *  indend(n),isornr(n),isorsr(n),islab(n)
#ifdef MPI 
      common /cord/ I1,I2,J1,J2,IA1,IA2,JA1,JA2,IB1,IB2,JB1,JB2
     1,ID1,ID2,JD1,JD2
      common /lbkhor/ khorl,lzet(khor),lb0(n),le0(n),
     * lb1(n),le1(n),lb2(n),le2(n)
#endif
#ifdef MPI
      do k=J1,J2
       lwa = lb0(k)
         lwe = le0(k)
#else
      do k=1,n 
         lwa = lb(k)
         lwe = le(k)
#endif
         if (lwa.le.lwe) then
            llb = indwet(lwa)+1
            lle = indwet(lwe)+lazc(lwe)
            do ll=llb,lle
               j = indver(ll)
               i = indi(ll)
               d1(ll) = d3(i,k,j) 
!               write(*,*) 'kurt',d1(ll)
            end do
         end if
      end do

      return 
      end




*��6 motmit
c-----------------------------------------------------------------------
      subroutine motmit (htmit,vtmit,vtmic,slip,botdic,ukonv,vkonv,
     * nndrei,ht,pah_fac)
c-----------------------------------------------------------------------
c       calculating the explicit terms in the equation of motion  
c       and integration of equation of motion for the calculation of zeta
c       we consider a trick-term for bottom friction to enable this.
c       See  Backhaus (1983) for details and explanations  
c-----------------------------------------------------------------------

      include 'C_model.f'
      include 'C_index.f'
      parameter(ilop1=ilo+1) 
      parameter(khor1=khor+1) 

      dimension ukonv(nndrei),vkonv(nndrei)
      dimension ht(m,n)
      double precision wo,sum,suma,sumt,sumdep,dpxy 
 
      common /ind/ iwet(khor1),ldep(khor),lazc(khor), 
     *indend(n),isornr(n),isorsr(n),islab(n) 
 
      common uc(ndrei),vc(ndrei),stc(ndrei),avc(ndrei),z(m,n) 
      common zac(khor),wobc(khor),stuvc(khor),fricu(khor),cxc(khor) 
      common cyc(khor),pac(khor),txc(khor),tyc(khor) 
      common stpc(ndrei),sac(ndrei),tec(ndrei) 
      common pres(ilo),wc(ndrei),fricv(khor) 

      common /dreh/ sinfu(m),sinfv(m),cosfu(m),cosfv(m),sincx(m), 
     ssincy(m),bx(m),by(m),pxu(m),pyv(m),pyu(m),pxv(m) 
c 
      common /gitter/ dt,r,g,dl,dlr,dlrh,dln(m),rdln(m),rad,dth, 
     ddlvo(m),dlvu(m),gh,rdt,dt2,r4,coru(m),corv(m),dtrdln(m),dtdlr 
      common /met/ windx,windy,wlam,pnull,rhoq(ilo),refrho(ilo),stress 


c     ------------ length of  common /aux/ = m*(ilo*19+10) 
      common /aux/ u(ilo,m,3),v(ilo,m,3), 
     1 pp_(ilo,m,3),st(ilo,m,3),w(ilo,m,3),ru(ilo,m),avv(ilo,m,2), 
     2 rv(ilo,m),sttu(m,2),za(m,3), 
     3 sttv(m),cx(m),cy(m)
      dimension tausx(khor),tausy(khor),pp(0:ndrei)
c     --------------- length of  common /auxint/ = m*(ilo*6+3) 
      common /auxint/ jc(ilo,m,3),lay(ilo,m,3),laz(m,3) 
c!!!!!!!!! 
      common/num/dz(ilo),av(ilo),ah(ilo),dh(ilo),pd(ilo), 
     pprd(ilo),pr2d(ilo)
      dimension r2d(ndrei),tkw(ndrei),tau(ndrei), d(ndrei), 
     qqa(ndrei),qbet(ilop1),qn(ndrei),rd(ilop1) 
      dimension fr(khor),sttuu(0:khor),sttvv(0:khor),jcc(m,3)
      logical init
      common /vec3/ init,iisu(khor),iisv(khor),lw2(-1:1,-1:1,khor),
     1 lww(-1:1,-1:1,ndrei),lev(ndrei),lind(ilo+1),jc2(khor),ley
      data init /.true./
      dimension taub(khor),sum(khor),sumtau(khor),sumt(khor),suma(khor),
     1 vc0(0:ndrei),tauo(khor),kl(ilo),ubot(khor),vbot(khor),slay(khor),
     2 zero(khor),dif(khor),uc0(0:ndrei),ley(ndrei)
      dimension pah_fac(ndrei)
c--------for calculating bottom shear stress in subroutine schubsp
      common /taubot/ taubu(khor), taubv(khor)
#ifdef MPI
      include 'C_mpi.f'
      include 'C_parallel.f'
#endif
c
      mz = m-1 
      nz = n-1 
c
      if (init) then
         init = .false.
c
      do k = 1,ilo
         kl(k) = 0
      end do
      do lw = 1,khor
         jc2(lw) = 0
      end do
c
      do j = 1,n
         do lw = lb(j),le(j)
            ldown = lazc(lw)
            do k = 1,ldown
               kl(k) = kl(k)+1
            end do
            if (ldown.gt.1) jc2(lw) = 1
         end do
      end do
c
      lind(1) = 0
      do k = 1,ilo
         lind(k+1) = lind(k)+kl(k)
      end do
c
      do k = 1,ilo
         kl(k) = 0
      end do
c
      do j = 1,n
         do lw = lb(j),le(j)
            nwet = indwet(lw)
            ldown = lazc(lw)
            do k = 1,ldown
               kl(k) = kl(k)+1
               lev(lind(k)+kl(k)) = nwet+k
               ley(nwet+k) = pd(k) 
            end do
         end do
      end do
c
         do kk = 1,2
            k = kk+1
            do lw = lb(kk),le(kk)
               i = iwet(lw)
               jcc(i,k) = lw
            end do
         end do
         do kk = 2,nz
            do l = 2,3
               do i = 1,m
                  jcc(i,l-1) = jcc(i,l)
               end do
            end do
            do i = 1,m
               jcc(i,3) = 0
            end do
            do lw = lb(kk+1),le(kk+1)
               i = iwet(lw)
               jcc(i,3) = lw
            end do
            do lw = lb(kk),le(kk)
               i = iwet(lw)
               do jj = -1,1
                  do ii = -1,1
                     lw2(ii,jj,lw) = jcc(i+ii,2+jj)
                  end do
               end do
            end do
         end do
c
         do kk = 1,2
            k = kk+1
            do i = 1,m
               do j = 1,ilo
                  jc(j,i,k) = 0
               end do
            end do
            do lw = lb(kk),le(kk)
               i = iwet(lw)
               ldown = lazc(lw)
               nwet = indwet(lw)
               do j = 1,ldown
                  jc(j,i,k) = nwet+j
               end do
            end do
         end do
         do kk = 2,nz
            do i = 1,m
               do j = 1,ilo
                  jc(j,i,1) = jc(j,i,2)
                  jc(j,i,2) = jc(j,i,3)
               end do
            end do
            do i = 1,m
               do j = 1,ilo
                  jc(j,i,3) = 0
               end do
            end do
            do lw = lb(kk+1),le(kk+1)
               i = iwet(lw)
               ldown = lazc(lw)
               nwet = indwet(lw)
               do j = 1,ldown
                  jc(j,i,3) = nwet+j
               end do
            end do
            do lw = lb(kk),le(kk)
               i = iwet(lw)
               ldown = lazc(lw)
               nwet = indwet(lw)
               do j = 1,ldown
                  do ii = -1,1
                     do jj = -1,1
                        lww(ii,jj,nwet+j) = jc(j,i+ii,2+jj)
                     end do
                  end do
               end do
            end do
         end do
c
         do lw = 1,khor
            iisu(lw) = 0
            iisv(lw) = 0
         end do
         do kk = 1,2 
            k = kk+1 
            do i = 1,m
               do j = 1,ilo
                  jc(j,i,k) = 0
               end do
               laz(i,k) = 0
            end do
            lwa = lb(kk)
            lwe = le(kk)
            if (lwa.le.lwe) then
               llb = indwet(lwa)+1
               lle = indwet(lwe)+lazc(lwe)
               do ll = llb,lle
                  i = indi(ll)
                  j = indver(ll)
                  jc(j,i,k) = 1 
               end do
            end if
            do lw = lb(kk),le(kk)
               i = iwet(lw) 
               laz(i,k) = lazc(lw) 
            end do
         end do
c
         do kk = 2,nz 
            do l = 2,3 
               do i = 1,m 
                  do j = 1,ilo 
                     jc(j,i,l-1) = jc(j,i,l) 
                  end do
                  laz(i,l-1) = laz(i,l) 
               end do
            end do
            do i = 1,m 
               laz(i,3) = 0 
            end do
            do i = 1,m 
               do j = 1,ilo 
                  jc(j,i,3) = 0 
               end do
            end do
            do lw = lb(kk+1),le (kk+1)
               i = iwet(lw) 
               ldown = lazc(lw) 
               do j=1,ldown 
                  jc(j,i,3) = 1 
               end do
               laz(i,3) = ldown 
            end do
            do lw = lb(kk),le(kk)
               i = iwet(lw) 
               ic = min0 (i-1,1) 
               icv = min0 (kk-2,1) 
               icu = min0 (i-2,1) 
               icu = max0 (icu,0) 
               jpu = jc(1,i,2)*jc(1,i,3)*icu 
               if (jpu.eq.1) iisu(lw) = min0(laz(i,2),laz(i,3)) 
               jpv = jc(1,i,2)*jc(1,i+1,2)*ic*icv 
               if (jpv.eq.1) iisv(lw) = min0(laz(i,2),laz(i+1,2)) 
            end do
         end do
      end if

      do j = 1,n
         do lw = lb(j),le(j)
            ldown = lazc(lw)
            nwet = indwet(lw)
            ley(nwet+ldown) = ldep(lw) 
         end do
      end do
c
      r4dt = 4.*r4*dt 
      rdyq = dlr*dlr 

c----------------------------------------------------------------
c    epsilon minimum layer thickness ##### should go to main
c---------------------------------------------------------------

      epsilon=0.5  ![m]

c 
      pp(0)=0
#ifdef MPI
      do kk=JA1,JA2
       do lw = lb1(kk),le1(kk)
#else
      do kk = 1,n
         do lw = lb(kk),le(kk)
#endif
            nwet = indwet(lw)
            ldown = lazc(lw) 
            du = 0.0 
            ps=pac(lw) 		 !1
            do j = 1,ldown 
c-------------------------------------------------------------------
c     pressure integral
c------------------------------------------------------------------

               do = du 
               du = dh(j)*stc(nwet+j) 
               ps = do+du+ps 
               pp(nwet+j) = ps 
            end do
         end do
      end do
c 
c+++++++++++++++++++++++++++++++++++++++++++++++ 
c     equation of motion, column loop 
c+++++++++++++++++++++++++++++++++++++++++++++++ 
c 
      uc0(0) = 0.
      vc0(0) = 0.
#ifdef MPI
      do J=JB1,JB2
      lwa=lb2(J)
      lwe=le2(J)
      if(lwe .ge. lwa) then
      nwet1=indwet(lwa)+1
      nwet2=indwet(lwe)+lazc(lwe)
      do l=nwet1,nwet2
         uc0(l)=uc(l)
         vc0(l)=vc(l)
      end do
      endif 
      end do
#else
      uc0(1:ndrei) = uc
      vc0(1:ndrei) = vc
#endif
      sttuu = 0.
      sttvv = 0.
c
c------------------------------------------- 
c     u - component 
c------------------------------------------- 

c+++++++++++++++++++++++++++++++++++++++++++++++++++ 
c     equation of motion, row loop 
c+++++++++++++++++++++++++++++++++++++++++++++++++++ 

c------------------------------------------------------------------
c  wind stress and vertical exchange
c------------------------------------------------------------------
#ifdef MPI
      do kk  = icej1,icej2
      lwa = lb0(kk)
      lwe = le0(kk)
#else
      do kk  = 2,nz
      lwe = le(kk)
      lwa = lb(kk)
#endif
      if (lwa.le.lwe) then
         llb = indwet(lwa)+1
         lle = indwet(lwe)+lazc(lwe)
!CDIR VECTOR NODEP
         do ll = llb,lle
            lw = llw(ll)
            j = ll-indwet(lw)
      if (iisu(lw).ne.0) then
      is = iisu(lw)
      if (j.eq.1) then
         
         d(ll) = pd(1) 
         r2d(ll) = 1. 
      else if (j.le.is) then
         d(ll) = pd(j) 
         qa(ll) = uc0(ll)/d(ll) 
c         honc=0.25 
         j00 = min(1,lww(0,0,ll))
         j0p = min(1,lww(0,1,ll))
         jm0 = min(1,lww(-1,0,ll))
         jmp = min(1,lww(-1,1,ll))
         hanc = real(j00+j0p+jm0+jmp) 
         honc = 1./hanc 
         qn(ll) = (vc0(ll)+vc0(lww(0,1,ll))
     1    +vc0(lww(-1,0,ll))+vc0(lww(-1,1,ll)))*honc 
         r2d(ll) = pr2d(j) 
      end if
      end if
      end do
      end if
      end do

   76 continue 
#ifndef MPI
      do j = 1,ilo
         do ih = lind(j)+1,lind(j+1)
            ll = lev(ih)
            lw = llw(ll)
      if (iisu(lw).ne.0) then
      is = iisu(lw)
      if (j.eq.1) then
         suma(lw) = uc0(ll) 
      else if (j.le.is) then
         suma(lw) = suma(lw)+uc0(ll) 
      end if
      end if
         end do
      end do
#else
      DO J=J1,J2
        DO lw=lb0(j),le0(j)
         nwet=indwet(lw)
         is = iisu(lw)
         if (is.ne.0) then
         do l=1,ilo
         if (l.eq.1) then
         suma(lw) = uc0(nwet+l)
         else if(l.le.is) then
         suma(lw) = suma(lw)+uc0(nwet+l)
         endif
      enddo
      endif
      enddo 
      enddo
#endif
#ifdef MPI
      do kk = icej1,icej2
      do lw = lb0(kk),le0(kk)
#else
      do kk = 2,nz
      do lw = lb(kk),le(kk)
#endif
      nwet = indwet(lw)
      rdich = stress/(rhoq(1)+stc(nwet+1)) 
c     rdich = stress 
      tausx(lw) = txc(lw)*rdich 
      if (iisu(lw).ne.0) then
      is = iisu(lw)
      ism = max0(is-1,1) 
      dreid = jc2(lw)*jc2(lw2(0,1,lw)) 
      zweid = 1.-dreid 
      tau(nwet+1) = tausx(lw) 

      d(nwet+is) = 0.5*(ley(nwet+is)+ley(lww(0,1,nwet+is)))              !Irina	 
      d(nwet+1) = max(epsilon,d(nwet+1)+0.5*(zac(lw)+zac(lw2(0,1,lw))))   
c      d(nwet+is) = amin0(ley(nwet+is),ley(lww(0,1,nwet+is)))
c      d(nwet+1) = max(epsilon,d(nwet+1)+amin1(zac(lw),zac(lw2(0,1,lw))))   
 
 
      qa(nwet+is) = uc0(nwet+is)/d(nwet+is)
      qa(nwet+1)  = uc0(nwet+1) /d(nwet+1)    
c      honc=0.25 
      j00 = min(1,lww(0,0,nwet+1))
      j0p = min(1,lww(0,1,nwet+1))
      jm0 = min(1,lww(-1,0,nwet+1))
      jmp = min(1,lww(-1,1,nwet+1))
      hanc = real(j00+j0p+jm0+jmp) 
      honc = 1./hanc 
      qn(nwet+1) = (vc0(nwet+1)+vc0(lww(0,1,nwet+1))
     1 +vc0(lww(-1,0,nwet+1))+vc0(lww(-1,1,nwet+1)))*honc 
      if (jc2(lw).eq.1) r2d(nwet+2) = 1./(d(nwet+1)+d(nwet+2)) 
      r2d(nwet+is) = dreid/(d(nwet+ism)+d(nwet+is))  + zweid 
      end if
      end do
      end do

#ifdef MPI
      do kk = icej1,icej2
      do lw = lb0(kk),le0(kk)
#else
      do kk = 2,nz
      do lw = lb(kk),le(kk)
#endif
      nwet = indwet(lw)
      if (iisu(lw).gt.1) then

      uvm = qa(nwet+1)
      uvu = qa(nwet+2) 
      wstu = wc(nwet+2)+wc(lww(0,1,nwet+2))

      tkw(nwet+1) = -0.25*wstu*(uvm+uvu)

      if (jc2(lw).eq.1) tau(nwet+2) = (avc(nwet+2)
     1 +avc(lww(0,1,nwet+2)))*(uvm-uvu) 

      end if
      end do
      end do

#ifdef MPI
      do kk  = ICEJ1,ICEJ2
      lwa = lb0(kk)
      lwe = le0(kk)
#else
      do kk  = 2,nz
      lwa = lb(kk)
      lwe = le(kk)
#endif
      if (lwa.le.lwe) then
         llb = indwet(lwa)+1
         lle = indwet(lwe)+lazc(lwe)
         do ll = llb,lle
            lw = llw(ll)
            j = ll-indwet(lw)
      nwet = indwet(lw)
      ldown = lazc(lw) 
      if (iisu(lw).ne.0) then
      is = iisu(lw)

c-------------------------------------------------------------
c    vertical shear (non-linear terms)
c-------------------------------------------------------------

      if (j.ge.3.and.j.le.is ) then
      uvo = qa(ll-2)
      uvm = qa(ll-1)
      uvu = qa(ll) 
      wsto = wc(ll-1)+wc(lww(0,1,ll-1))
cc      wgsto=wgstu 
      wstu = wc(ll)+wc(lww(0,1,ll))
cc      wgstu=wstu 
      tkw(ll-1) = 0.25*(wsto*(uvo+uvm)-wstu*(uvm+uvu)) 
cc     --0.5*uvm*(wgsto - wgstu) 
 1212 tau(ll) = (avc(ll)+avc(lww(0,1,ll)))*(uvm-uvu) 
      end if
      end if
      end do
      end if
      end do

#ifdef MPI
      do kk = ICEJ1,ICEJ2
      do lw = lb0(kk),le0(kk)
#else
      do kk = 2,nz
      do lw = lb(kk),le(kk)
#endif
      if (iisu(lw).ne.0) then
      is = iisu(lw)
      ism = max(1,is-1)
      nwet = indwet(lw)
      uvo = qa(nwet+ism)
      uvm = qa(nwet+is)
      if (is.gt.1) then
      wstu = wc(nwet+is)+wc(lww(0,1,nwet+is))
      else
      wstu = 0.
      end if
      wsto = wstu 
ccn      wgstu=wboden(i,kk)+wboden(i,kk+1)
      tkw(nwet+is) = 0.25*(wsto*(uvo+uvm)) 
cc     --0.5*uvm*(wsto-wgstu) 
      end if
      end do
      end do

#ifdef MPI
      do kk = icej1,icej2
      do lw = lb0(kk),le0(kk)
#else
      do kk = 2,nz
      do lw = lb(kk),le(kk)
#endif
      if (iisu(lw).ne.0) then
      is = iisu(lw)

c---------------------------------------------------------------------- 
c     bottom friction (tau(is+1)) 
c     compute  pre-defined (botdic) bottom layer for bottom friction 
c     and vertically intergate the flow over the (fictive) layer botdic 
c----------------------------------------------------------------------
      ubot(lw) = 0. 
      vbot(lw) = 0. 
      slay(lw) = 0. 
      zero(lw) = 1. 
      dif(lw) = botdic 
      end if
      end do
      end do
#ifndef MPI
      do j = ilo,1,-1 
!CDIR VECTOR NODEP
         do ih = lind(j)+1,lind(j+1)
            ll = lev(ih)
            lw = llw(ll)
      if (iisu(lw).ne.0) then
      is = iisu(lw)

c----------------------------------------------------------------------- 
c     bottom friction (tau(is+1)) 
c     compute  pre-defined (botdic) bottom layer for bottom friction 
c     and vertically intergate the flow over the (fictive) layer botdic 
c-----------------------------------------------------------------------

      if (j.le.is) then
      dif(lw) = dif(lw)-d(ll) 
      eins = 0.5*(1.+sign(1.,dif(lw))) 
      red = (eins+(1.-eins)*(1.+dif(lw)/d(ll)))*zero(lw) 
      ubot(lw) = ubot(lw) + uc0(ll)*red 
      vbot(lw) = vbot(lw) + qn(ll)   *red 
      slay(lw) = slay(lw)+d(ll)*red 
      zero(lw) = eins 
      end if
      end if
      end do
      end do
#else
        DO J=J1,J2
        DO lw=lb0(j),le0(j)
         nwet=indwet(lw)
         is = iisu(lw)
         if (is.ne.0) then
         do l=ilo,1,-1
         if(l.le.is) then
         k=nwet+l
         dif(lw) = dif(lw)-d(k) 
         eins = 0.5*(1.+sign(1.,dif(lw))) 
         red = (eins+(1.-eins)*(1.+dif(lw)/d(k)))*zero(lw) 
         ubot(lw) = ubot(lw) + uc0(k)*red 
         vbot(lw) = vbot(lw) + qn(k)   *red 
         slay(lw) = slay(lw)+d(k)*red 
         zero(lw) = eins
         endif
      enddo
      endif
      enddo 
      enddo
#endif
#ifdef MPI
      do kk = icej1,icej2
      do lw = lb0(kk),le0(kk)
#else
      do kk = 2,nz
      do lw = lb(kk),le(kk)
#endif
      nwet = indwet(lw)
      if (iisu(lw).ne.0) then
      is = iisu(lw)

c-----------------------------------------------------------------------------
c     frac = rdt*sqrt(ubot**2+vbot**2)/slay**2 
c     scaling factor for transport normalised to botdic 
c     is slay/d(nwet+is), giving the transport for a fictive 
c     bottom layer botdic with direction of u(is) 
c     scaling factor times frac gives / 
c     slay * rdt..... /(slay**2*d(nwet+is)) , one slay cancels 
c     and hence division by (slay*d(nwet+is)) : 
      frac = rdt*sqrt(ubot(lw)**2+vbot(lw)**2)/(slay(lw)*d(nwet+is)) 
c 
c 
c      frac = r4dt*sqrt((u(is,i,2))**2+(qn(nwet+is))**2)/d(nwet+is)**2 
c-------------------------------------------------------------------------------
      rimp = 1./(1.+frac) 
      fricu(lw) = rimp 
      fr(lw) = rimp 
      end if
      end do
      end do

#ifdef MPI
      do kk = icej1,icej2
      do lw = lb0(kk),le0(kk)
#else
      do kk = 2,nz
      do lw = lb(kk),le(kk)
#endif
      i = iwet(lw) 
      nwet = indwet(lw)
      if (iisu(lw).ne.0) then
      is = iisu(lw)
      ism = max0(is-1,1) 
      dreid = jc2(lw)*jc2(lw2(0,1,lw)) 
      zweid = 1.-dreid 
      tauo(lw) = 0.0 
      taub(lw) = tau(nwet+is)*r2d(nwet+is)*(1.-fricu(lw)) * dreid 
c----------keeping values
      taubu(lw)=taub(lw)
      sumtau(lw) = dt*(tau(nwet+1)-taub(lw))  * (fr(lw)*zweid+dreid) 
      sum(lw) = 0. 
      cxc(lw) = 0.d0 
      end if
      end do
      end do

#ifndef MPI
      do j = ilo,1,-1
!CDIR VECTOR NODEP
         do ih = lind(j)+1,lind(j+1)
            ll = lev(ih)
            lw = llw(ll)
      i = iwet(lw) 
      if (iisu(lw).ne.0) then
      is = iisu(lw)
      ism = max0(is-1,1) 
      cosfuc = cosfu(i) 
      sinfuc = sinfu(i) 
      px = pxu(i) 
      py = pyu(i) 
      rdxq = rdln(i)**2 
      if (j.le.is) then
      tauu = tauo(lw)
      tauo(lw) = tau(ll)*r2d(ll) 
c----------------------------------------------------------------------
c   pressure gradient term with primary and secodary pressure terms, 
c    after wais (referenced by Backhaus (1983)
c----------------------------------------------------------------------

      rost = d(ll) 
     //(rhoq(j)+0.5*(stc(ll)+stc(lww(0,1,ll)))) 
      cxc(lw) = cxc(lw)+rost*fr(lw) 
      jm0 = min(1,lww(-1,0,ll))
      jmp = min(1,lww(-1,1,ll))
      j00 = min(1,lww(0,0,ll))
      j0p = min(1,lww(0,1,ll))
      jp0 = min(1,lww(1,0,ll))
      jpp = min(1,lww(1,1,ll))
      gne=4./real(jm0+jmp+j00+j0p) 
      gse=4./real(jp0+jpp+j00+j0p) 
      dpxy = rost*(dprod(px,pp(lww(0,1,ll)))-dprod(px,pp(ll)) 
     1 +dprod(py,gne*(pp(lww(-1,1,ll))+pp(lww(-1,0,ll))
     2 +pp(ll)+pp(lww(0,1,ll))) 
     3 -gse*(pp(lww(1,1,ll))+pp(lww(1,0,ll))+pp(ll)
     4 +pp(lww(0,1,ll))))) 
c--------------------------------------------------------------------------
c   horizontal turbulence 
c--------------------------------------------------------------------------

      diefe = amin1(d(ll),600.) 
      ruv = ah(j) * diefe*pah_fac(ll)*
     1 (rdxq*(uc0(lww(0,-1,ll))-2.*uc0(ll)+
     2 uc0(lww(0,1,ll)))+rdyq*(uc0(lww(-1,0,ll))
     3 +uc0(lww(1,0,ll))- 
     f (2.-slip*(2.-real(jm0+jp0)))*uc0(ll))) 

c----------------------------------------------------------------------------
c    new advection scheme arakawa J7 replaced vektor-upstream, 
c    non-linear terms calculated in subroutine konv, called from main
c    program; nonlinear terms on ukonv and vkonv fields
c---------------------------------------------------------------------------
 
      uik = uc0(ll) 
      vst = qn(ll) 

c----------------------------------------------------------------
c   all terms which are explicitly forward are combined in sumt 
c----------------------------------------------------------------
c 
      sumt(lw) = dprod(uik,cosfuc)+dprod(vst,sinfuc)-dpxy 
     ++dprod(dt,ruv)-dprod(dt,tkw(ll)) +ukonv(ll) 

      sum(lw) = sum(lw)+sumt(lw)*fr(lw) 
      fr(lw) = 1. 
      uc(ll) = sumt(lw)+vtmic*dt*(tauo(lw)-tauu) 
      end if
      end if
      end do
      end do
#else
        DO J=J1,J2
        DO lw=lb0(j),le0(j)
         nwet=indwet(lw)
         is = iisu(lw)
         i = iwet(lw)
        if (is.ne.0) then
          ism = max0(is-1,1)
         cosfuc = cosfu(i) 
      sinfuc = sinfu(i) 
      px = pxu(i) 
      py = pyu(i) 
      rdxq = rdln(i)**2
         do l=ilo,1,-1
         ll=nwet+l     
      if (l.le.is) then
      tauu = tauo(lw)
      tauo(lw) = tau(ll)*r2d(ll) 
c----------------------------------------------------------------------
c   pressure gradient term with primary and secodary pressure terms, 
c    after wais (referenced by Backhaus (1983)
c----------------------------------------------------------------------
      rost = d(ll) 
     //(rhoq(l)+0.5*(stc(ll)+stc(lww(0,1,ll)))) 
      cxc(lw) = cxc(lw)+rost*fr(lw) 
      jm0 = min(1,lww(-1,0,ll))
      jmp = min(1,lww(-1,1,ll))
      j00 = min(1,lww(0,0,ll))
      j0p = min(1,lww(0,1,ll))
      jp0 = min(1,lww(1,0,ll))
      jpp = min(1,lww(1,1,ll))
      gne=4./real(jm0+jmp+j00+j0p) 
      gse=4./real(jp0+jpp+j00+j0p) 
      dpxy = rost*(dprod(px,pp(lww(0,1,ll)))-dprod(px,pp(ll)) 
     1 +dprod(py,gne*(pp(lww(-1,1,ll))+pp(lww(-1,0,ll))
     2 +pp(ll)+pp(lww(0,1,ll))) 
     3 -gse*(pp(lww(1,1,ll))+pp(lww(1,0,ll))+pp(ll)
     4 +pp(lww(0,1,ll))))) 
c--------------------------------------------------------------------------
c   horizontal turbulence 
c--------------------------------------------------------------------------
      diefe = amin1(d(ll),600.) 
      ruv = ah(l) * diefe*pah_fac(ll)*
     1 (rdxq*(uc0(lww(0,-1,ll))-2.*uc0(ll)+
     2 uc0(lww(0,1,ll)))+rdyq*(uc0(lww(-1,0,ll))
     3 +uc0(lww(1,0,ll))- 
     f (2.-slip*(2.-real(jm0+jp0)))*uc0(ll))) 
c----------------------------------------------------------------------------
c    new advection scheme arakawa J7 replaced vektor-upstream, 
c    non-linear terms calculated in subroutine konv, called from main
c---------------------------------------------------------------------------
      uik = uc0(ll) 
      vst = qn(ll) 
c----------------------------------------------------------------
c   all terms which are explicitly forward are combined in sumt 
c----------------------------------------------------------------
c 
      sumt(lw) = dprod(uik,cosfuc)+dprod(vst,sinfuc)-dpxy 
     ++dprod(dt,ruv)-dprod(dt,tkw(ll)) +ukonv(ll) 
      sum(lw) = sum(lw)+sumt(lw)*fr(lw) 
      fr(lw) = 1. 
      uc(ll) = sumt(lw)+vtmic*dt*(tauo(lw)-tauu) 
      end if
      end do
      end if
      end do
      end do
#endif
#ifdef MPI
      do kk = icej1,icej2
      do lw = lb0(kk),le0(kk)
#else
      do kk = 2,nz
      do lw = lb(kk),le(kk)
#endif
      if (iisu(lw).ne.0) then
      is = iisu(lw)
      ism = max0(is-1,1) 
      nwet = indwet(lw)
      uc(nwet+ism) = uc(nwet+ism)-vtmit*dt*taub(lw)
      uc(nwet+1) = uc(nwet+1)+vtmit*dt*tausx(lw) 
      sttuu(lw) = htmit*(sum(lw)+sumtau(lw)-suma(lw))+suma(lw) 
      end if
      end do
      end do
#ifdef MPI      
c     communication sttuu and cxc
      call INTI_COMM_VAR
      call com_horizontal(sttuu(1),khor,dtemp,iel,0,1,1,2,0,icdepb0,
     *idepb0,icdept0,
     *idept0,ncsendb0,nsendb0,ncsendt0,nsendt0,nbhlayer0,nthlayer0)
      call com_horizontal(cxc,khor,dtemp,iel,0,1,1,2,1,icdepb0,
     *idepb0,icdept0,
     *idept0,ncsendb0,nsendb0,ncsendt0,nsendt0,nbhlayer0,nthlayer0)
#endif
c------------------------------------------------------------------- 
c     v - component 
c-------------------------------------------------------------------
 
#ifdef MPI
      do kk  = icej1,icej2
      lwa = lb0(kk)
      lwe = le0(kk)
#else
      do kk=2,nz
      lwa = lb(kk)
      lwe = le(kk)
#endif
      if (lwa.le.lwe) then
         llb = indwet(lwa)+1
         lle = indwet(lwe)+lazc(lwe)
         do ll = llb,lle
            lw = llw(ll)
            j = ll-indwet(lw)
      i = iwet(lw) 
      nwet = indwet(lw)
      ldown = lazc(lw) 
      if (iisv(lw).ne.0) then
      is = iisv(lw)
      ism = max0(is-1,1) 
      if (j.eq.1) then
         d(nwet+1) = pd(1) 
         r2d(nwet+1) = 1. 
      else if (j.le.is) then
         d(nwet+j) = pd(j) 
         qa(nwet+j) = vc0(nwet+j)/d(nwet+j) 
c         honc=0.25 
         j0m = min(1,lww(0,-1,nwet+j))
         j00 = min(1,lww(0,0,nwet+j))
         jpm = min(1,lww(1,-1,nwet+j))
         jp0 = min(1,lww(1,0,nwet+j))
         hanc = real(j0m+j00+jpm+jp0) 
         honc = 1./hanc 
         qn(nwet+j) = (uc0(nwet+j)+uc0(lww(1,0,nwet+j))
     1    +uc0(lww(0,-1,nwet+j))+uc0(lww(1,-1,nwet+j)))*honc 
         r2d(nwet+j) = pr2d(j) 
      end if
      end if
      end do
      end if
      end do
#ifndef MPI 
      do j = 1,ilo
!CDIR VECTOR NODEP
         do ih = lind(j)+1,lind(j+1)
            ll = lev(ih)
            lw = llw(ll)
      if (iisv(lw).ne.0) then
      is = iisv(lw)
      if (j.eq.1) then
         suma(lw) = vc0(ll) 
      else if (j.le.is) then
         suma(lw) = suma(lw)+vc0(ll) 
      end if
      end if
         end do
      end do
#else
        DO J=J1,J2
        DO lw=lb0(j),le0(j)
         nwet=indwet(lw)
         is = iisv(lw)
         if (is.ne.0) then
         do l=1,ilo
         if (l.eq.1) then
         suma(lw) = vc0(nwet+l)
         else if(l.le.is) then
         suma(lw) = suma(lw)+vc0(nwet+l)
         endif
         enddo
      endif
      enddo 
      enddo
#endif
#ifdef MPI
      do kk = icej1,icej2
      do lw = lb0(kk),le0(kk)
#else
      do kk = 2,nz
      do lw = lb(kk),le(kk)
#endif
      nwet = indwet(lw)
      ldown = lazc(lw) 
      rdich = stress/(rhoq(1)+stc(nwet+1)) 
c     rdich = stress 
      tau(nwet+1) = tyc(lw)*rdich 
      if (iisv(lw).ne.0) then
      is = iisv(lw)
      ism = max0(is-1,1) 
      dreid = jc2(lw)*jc2(lw2(1,0,lw)) 
      zweid = 1.-dreid 
c      honc=0.25 
      j0m = min(1,lww(0,-1,nwet+1))
      j00 = min(1,lww(0,0,nwet+1))
      jpm = min(1,lww(1,-1,nwet+1))
      jp0 = min(1,lww(1,0,nwet+1))
      hanc = real(j0m+j00+jpm+jp0) 
      honc = 1./hanc 
      qn(nwet+1) = (uc0(nwet+1)+uc0(lww(1,0,nwet+1))
     1 +uc0(lww(0,-1,nwet+1))+uc0(lww(1,-1,nwet+1)))*honc 

      d(nwet+is) = 0.5*(ley(nwet+is)+ley(lww(1,0,nwet+is)))            !Irina
      d(nwet+1) =max(epsilon,d(nwet+1)+0.5*(zac(lw)+zac(lw2(1,0,lw)))) 
c      d(nwet+is) = amin0(ley(nwet+is),ley(lww(1,0,nwet+is))) 
c      d(nwet+1) =max(epsilon,d(nwet+1)+amin1(zac(lw),zac(lw2(1,0,lw)))) !irina

      qa(nwet+is) = vc0(nwet+is)/d(nwet+is) 
      qa(nwet+1) = vc0(nwet+1) / d(nwet+1 ) 
      if (jc2(lw).eq.1) r2d(nwet+2) = 1./(d(nwet+1)+d(nwet+2)) 
      r2d(nwet+is) = dreid/(d(nwet+ism)+d(nwet+is))  + zweid 
      end if
      end do
      end do


#ifdef MPI
      do kk = icej1,icej2 
      do lw = lb0(kk),le0(kk)
#else
      do kk = 2,nz
      do lw = lb(kk),le(kk)
#endif
      nwet = indwet(lw)
      if (iisv(lw).gt.1) then
      uvm = qa(nwet+1)
      uvu = qa(nwet+2) 
      wstu = wc(nwet+2)+wc(lww(1,0,nwet+2)) 
      tkw(nwet+1) = -0.25*wstu*(uvm+uvu)  !Val
c      tkw(nwet+1) =   0.25*wstu*(uvm+uvu)  !Irina
      if (jc2(lw).eq.1) tau(nwet+2) = (avc(nwet+2)
     1 +avc(lww(1,0,nwet+2)))*(uvm-uvu) 
      end if
      end do
      end do

#ifdef MPI
      do kk  = icej1,icej2 
      lwa = lb0(kk)
      lwe = le0(kk)
#else
      do kk = 2,nz
      lwa= lb(kk)
      lwe= le(kk)
#endif

      if (lwa.le.lwe) then
         llb = indwet(lwa)+1
         lle = indwet(lwe)+lazc(lwe)
         do ll = llb,lle
            lw = llw(ll)
            j = ll-indwet(lw)
      nwet = indwet(lw)
      ldown = lazc(lw) 
      if (iisv(lw).ne.0) then
      is = iisv(lw)
c     ------ nichtlin. term und vertical shear 
      if (j.ge.3.and.j.le.is ) then
      uvo = qa(ll-2)
      uvm = qa(ll-1)
      uvu = qa(ll) 
      wsto = wc(ll-1)+wc(lww(1,0,ll-1))
c      wgsto=wgstu 
      wstu = wc(ll)+wc(lww(1,0,ll)) 
c      wgstu=wstu
      tkw(ll-1) = 0.25*(wsto*(uvo+uvm)-wstu*(uvm+uvu)) 
cc     --0.5*uvm*(wgsto - wgstu) 
      tau(ll) = (avc(ll)+avc(lww(1,0,ll)))*(uvm-uvu) 
      end if 
      end if
      end do
      end if
      end do

#ifdef MPI
      do kk = icej1,icej2 
      do lw = lb0(kk),le0(kk)
#else
      do kk = 2,nz
      do lw = lb(kk),le(kk)
#endif
      if (iisv(lw).ne.0) then
      is = iisv(lw)
      ism = max(1,is-1)
      nwet = indwet(lw)
      uvo = qa(nwet+ism)
      uvm = qa(nwet+is)
      if (is.gt.1) then
      wstu = wc(nwet+is)+wc(lww(1,0,nwet+is)) 
      else
      wstu = 0.
      end if
      wsto = wstu
cc      wgstu=wboden(i,kk)+wboden(i+1,kk) 
      tkw(nwet+is) = 0.25*(wsto*(uvo+uvm)) 
cc     --0.5*uvm*(wsto-wgstu) 
      end if
      end do
      end do

#ifdef MPI
      do kk =  icej1,icej2
!CDIR VECTOR NODEP
      do lw = lb0(kk),le0(kk)
#else
!CDIR VECTOR NODEP
      do kk = 2,nz
      do lw = lb(kk),le(kk)
#endif
      if (iisv(lw).ne.0) then
c-------------------------------------------------------------------- 
c     bottom friction (tau(is+1)) 
c     compute  pre-defined (botdic) bottom layer for bottom friction 
c     and vertically intergate the flow over the (fictive) layer botdic 
c---------------------------------------------------------------------

      ubot(lw) = 0. 
      vbot(lw) = 0. 
      slay(lw) = 0. 
      zero(lw) = 1. 
      dif(lw) = botdic 
      end if
      end do
      end do
#ifndef MPI 
      do j = ilo,1,-1 
!CDIR VECTOR NODEP
         do ih = lind(j)+1,lind(j+1)
            ll = lev(ih)
            lw = llw(ll)
      if (iisv(lw).ne.0) then
      is = iisv(lw)
      if (j.le.is) then
      dif(lw) = dif(lw)-d(ll) 
      eins = 0.5*(1.+sign(1.,dif(lw))) 
      red = (eins+(1.-eins)*(1.+dif(lw)/d(ll)))*zero(lw) 
      vbot(lw) = vbot(lw) + vc0(ll)*red 
      ubot(lw) = ubot(lw) + qn(ll)   *red 
      slay(lw) = slay(lw)+d(ll)*red 
      zero(lw) = eins 
      end if
      end if
      end do
      end do
#else
        DO J=J1,J2
        DO lw=lb0(j),le0(j)
         nwet=indwet(lw)
         is = iisv(lw)
         if (is.ne.0) then
         do l=ilo,1,-1
         if(l.le.is) then
         ll=nwet+l
         dif(lw) = dif(lw)-d(ll) 
      eins = 0.5*(1.+sign(1.,dif(lw))) 
      red = (eins+(1.-eins)*(1.+dif(lw)/d(ll)))*zero(lw) 
      vbot(lw) = vbot(lw) + vc0(ll)*red 
      ubot(lw) = ubot(lw) + qn(ll)   *red 
      slay(lw) = slay(lw)+d(ll)*red 
      zero(lw) = eins
         endif
         enddo
         endif
      enddo 
      enddo
#endif
#ifdef MPI
      do kk =  icej1,icej2
!CDIR VECTOR
      do lw = lb0(kk),le0(kk)
#else
!CDIR VECTOR
      do kk = 2,nz
      do lw = lb(kk),le(kk)
#endif
      nwet = indwet(lw)
      if (iisv(lw).ne.0) then
      is = iisv(lw)
c---------------------------------------------------------------------------
c     frac = rdt*sqrt(ubot**2+vbot**2)/slay**2 
c     scaling factor for transport normalised to botdic 
c     is slay/d(nwet+is), giving the transport for a fictive 
c     bottom layer botdic with direction of u(is) 
c     scaling factor times frac gives / 
c     slay * rdt..... /(slay**2*d(nwet+is)) , one slay cancels 
c     and hence division by (slay*d(nwet+is)) : 
c---------------------------------------------------------------------------
      frac = rdt*sqrt(ubot(lw)**2+vbot(lw)**2)/(slay(lw)*d(nwet+is)) 
c 
c      frac = r4dt*sqrt((v(is,i,2))**2+(qn(nwet+is))**2)/d(nwet+is)**2 
      rimp = 1./(1.+frac) 
      fricv(lw) = rimp 
      fr(lw) = rimp 
      end if
      end do
      end do

#ifdef MPI
      do kk = icej1,icej2       
      do lw = lb0(kk),le0(kk)
#else
!CDIR VECTOR NODEP
      do kk = 2,nz
      do lw = lb(kk),le(kk)
#endif
      nwet = indwet(lw)
      ldown = lazc(lw) 
      if (iisv(lw).ne.0) then
      is = iisv(lw)
      dreid = jc2(lw)*jc2(lw2(1,0,lw)) 
      zweid = 1.-dreid 
      tauo(lw) = 0.0 
      taub(lw) = tau(nwet+is)*r2d(nwet+is)*(1.-fricv(lw)) * dreid 
c----------keeping values
      taubv(lw)=taub(lw)
      sumtau(lw) = dt*(tau(nwet+1)-taub(lw)) * (fr(lw)*zweid+dreid) 
      sum(lw) = 0. 
      cyc(lw) = 0.d0 
      end if
      end do
      end do
#ifndef MPI 
      do j = ilo,1,-1
!CDIR VECTOR NODEP
         do ih = lind(j)+1,lind(j+1)
            ll = lev(ih)
            lw = llw(ll)
      i = iwet(lw) 
      ldown = lazc(lw) 
      if (iisv(lw).ne.0) then
      is = iisv(lw)
      cosfuc = cosfv(i) 
      sinfuc = sinfv(i) 
      py = pyv(i) 
      px = pxv(i) 
      rdxq = rdln(i)**2 
      if (j.le.is) then
      tauu = tauo(lw) 
      tauo(lw) = tau(ll)*r2d(ll) 
c
c------------------------------------------------------------------------
c    pressure gradient terms with primary and secondary pressure terms 
c    according to wais (referred by Backhaus, 1983)
c-----------------------------------------------------------------------

      rost = d(ll) 
     //(rhoq(j)+0.5*(stc(ll)+stc(lww(1,0,ll)))) 
      cyc(lw) = cyc(lw)+rost*fr(lw) 
c 
      jp0 = min(1,lww(1,0,ll))
      jpm = min(1,lww(1,-1,ll))
      j00 = min(1,lww(0,0,ll))
      j0m = min(1,lww(0,-1,ll))
      jpp = min(1,lww(1,1,ll))
      j0p = min(1,lww(0,1,ll))
c 
      gsw=4./real(jp0+jpm+j00+j0m) 
      gse=4./real(jp0+jpp+j00+j0p) 
      dpxy = rost*(dprod(py,pp(ll))-dprod(py,pp(lww(1,0,ll))) 
     1 +dprod(px,gse*(pp(lww(0,1,ll))+pp(ll)
     2 +pp(lww(1,1,ll))+pp(lww(1,0,ll))) 
     3 -gsw*(pp(lww(0,-1,ll))+pp(ll)+pp(lww(1,-1,ll))
     4 +pp(lww(1,0,ll))))) 
c     horizontaler austausch 
      diefe = amin1(d(ll),600.) 
      ruv = ah(j) * diefe*pah_fac(ll)*
     1 (rdyq*(vc0(lww(-1,0,ll))-2.*vc0(ll)
     2 +vc0(lww(1,0,ll)))+ 
     f rdxq*(vc0(lww(0,-1,ll))+vc0(lww(0,1,ll))- 
     f (2.-slip*(2.-real(j0m+j0p)))*vc0(ll))) 

c----------------------------------------------------------------------
c   new advection scheme arakawa J7 replaced vektor upstream
c   nonlinear terms are calculated in subroutine konv and stored on
c   ukonv and vkonv fields, call of subroutine konv from main program
c-----------------------------------------------------------------------

 
      ust = qn(ll) 
      vik = vc0(ll) 
c
c-------------------------------------------------------------------
c   summarize explicit terms in sumt 
c-------------------------------------------------------------------
c 
      sumt(lw) = dprod(vik,cosfuc)-dprod(ust,sinfuc)-dpxy 
     ++dprod(dt,ruv)-dprod(dt,tkw(ll)) +vkonv(ll) 

      sum(lw) = sum(lw)+ sumt(lw)*fr(lw) 
      fr(lw) = 1. 
      vc(ll) = sumt(lw)+vtmic*dt*(tauo(lw)-tauu) 
      end if
      end if
      end do
      end do
#else
        DO J=J1,J2
        DO lw=lb0(j),le0(j)
         nwet=indwet(lw)
         is = iisv(lw)
         i = iwet(lw) 
         if (is.ne.0) then
         cosfuc = cosfv(i) 
      sinfuc = sinfv(i) 
      py = pyv(i) 
      px = pxv(i) 
      rdxq = rdln(i)**2
         do l=ilo,1,-1
          if(l.le.is) then
         ll=nwet+l
       tauu = tauo(lw) 
      tauo(lw) = tau(ll)*r2d(ll) 
c
c------------------------------------------------------------------------
c    pressure gradient terms with primary and secondary pressure terms 
c    according to wais (referred by Backhaus, 1983)
c-----------------------------------------------------------------------
      rost = d(ll) 
     //(rhoq(l)+0.5*(stc(ll)+stc(lww(1,0,ll)))) 
      cyc(lw) = cyc(lw)+rost*fr(lw) 
c 
      jp0 = min(1,lww(1,0,ll))
      jpm = min(1,lww(1,-1,ll))
      j00 = min(1,lww(0,0,ll))
      j0m = min(1,lww(0,-1,ll))
      jpp = min(1,lww(1,1,ll))
      j0p = min(1,lww(0,1,ll))
c 
      gsw=4./real(jp0+jpm+j00+j0m) 
      gse=4./real(jp0+jpp+j00+j0p) 
      dpxy = rost*(dprod(py,pp(ll))-dprod(py,pp(lww(1,0,ll))) 
     1 +dprod(px,gse*(pp(lww(0,1,ll))+pp(ll)
     2 +pp(lww(1,1,ll))+pp(lww(1,0,ll))) 
     3 -gsw*(pp(lww(0,-1,ll))+pp(ll)+pp(lww(1,-1,ll))
     4 +pp(lww(1,0,ll))))) 
c     horizontaler austausch 
      diefe = amin1(d(ll),600.) 
      ruv = ah(l) * diefe*pah_fac(ll)*
     1 (rdyq*(vc0(lww(-1,0,ll))-2.*vc0(ll)
     2 +vc0(lww(1,0,ll)))+ 
     f rdxq*(vc0(lww(0,-1,ll))+vc0(lww(0,1,ll))- 
     f (2.-slip*(2.-real(j0m+j0p)))*vc0(ll))) 
c----------------------------------------------------------------------
c   new advection scheme arakawa J7 replaced vektor upstream
c   nonlinear terms are calculated in subroutine konv and stored on
c   ukonv and vkonv fields, call of subroutine konv from main program
c-----------------------------------------------------------------------
      ust = qn(ll) 
      vik = vc0(ll) 
c
c-------------------------------------------------------------------
c   summarize explicit terms in sumt 
c-------------------------------------------------------------------
c 
      sumt(lw) = dprod(vik,cosfuc)-dprod(ust,sinfuc)-dpxy 
     ++dprod(dt,ruv)-dprod(dt,tkw(ll)) +vkonv(ll) 

      sum(lw) = sum(lw)+ sumt(lw)*fr(lw) 
      fr(lw) = 1. 
      vc(ll) = sumt(lw)+vtmic*dt*(tauo(lw)-tauu) 
      end if
      enddo
      endif
      enddo 
      enddo
#endif
#ifdef MPI
      do kk = icej1,icej2
!CDIR VECTOR NODEP
      do lw = lb0(kk),le0(kk)
#else
!CDIR VECTOR NODEP
      do kk = 2,nz
      do lw = lb(kk),le(kk)
#endif
      if (iisv(lw).ne.0) then
      is = iisv(lw) 
      ism = max0(is-1,1) 
      i = iwet(lw)
      nwet = indwet(lw)
      rdich = stress/(rhoq(1)+stc(nwet+1)) 
c     rdich = stress 
      tausy(lw) = tyc(lw)*rdich 
      vc(nwet+1) = vc(nwet+1)+vtmit*dt*tausy(lw) 
      vc(nwet+ism) = vc(nwet+ism)-vtmit*dt*taub(lw)
      sttvv(lw) = htmit*(sum(lw)+sumtau(lw)-suma(lw))+suma(lw) 
      end if
      end do
      end do
#ifdef MPI      
c     communication sttvv and cyc
      call com_vertical(sttvv(1),khor,dtemp,1,1,0,1,2,0,lazc)
      call com_vertical(cyc,khor,dtemp,1,1,0,1,2,1,lazc)

      call  MPI_Waitall(nreq1, ir_r_b, status, ierr)
      INUM3=0
      INUM4=0     
      call HOR_UNPACK(sttuu(1),khor,dtemp,1,0,1,nbhlayer0,nthlayer0,
     &idepb0,idept0,icdepb0,icdept0,1,2,iwet)
      call HOR_UNPACK(cxc,khor,dtemp,1,0,1,nbhlayer0,nthlayer0,
     &idepb0,idept0,icdepb0,icdept0,1,2,iwet)
      INUM1=0
      INUM2=0
      call VER_UNPACK(sttvv(1),khor,dtemp,1,1,0,1,2,lazc)
      call VER_UNPACK(cyc,khor,dtemp,1,1,0,1,2,lazc)
#endif

#ifdef MPI
      do kk =icej1,icej2         
         do lw = lb0(kk),le0(kk)
#else
      do kk = 2,nz
         do lw = lb(kk),le(kk)
#endif
            i = iwet(lw) 
            dvo = dlvo(i)*dlr 
            dvu = dlvu(i)*dlr 
            stuvc(lw) 
     1      = dt *(dprod(rdln(i),sttuu(lw2(0,-1,lw)))
     2       -dprod(rdln(i),sttuu(lw)) 
     3       +dprod(dvu,sttvv(lw))-dprod(dvo,sttvv(lw2(-1,0,lw))))
         end do
      end do
      return 
      end 



*��7 sor
c-----------------------------------------------------------------------
      subroutine sor(htmit,nosor,eps,ip) 
c-----------------------------------------------------------------------
c     solving the surface elevation equation with red-balack solver  
c-----------------------------------------------------------------------
      include 'C_model.f'
      include 'C_index.f'
c
      parameter(ilop1=ilo+1) 
      parameter(khor1=khor+1) 
c 
      double precision dom,hz 
c 
      common /ind/ iwet(khor1),ldep(khor),lazc(khor), 
     cindend(n),isornr(n),isorsr(n),islab(n) 
c 
      integer nord,sued,ost,west 
      integer*2 isor,ksor 
c 
      common /intit/ itel(1440),jc(m,3),ikk(kasor),isor(kasor), 
     cksor(kasor),nofrez,lone,itende,nfreez 
c 
      common /iterat/ omega(kasor),za(m,3),cx(m,3),cy(m,3),rhoz(m,3), 
     ctcc(8) 
c 
      common /uirand/ ilwrnd(lrp),islrnd(lrp),ldownr(lrp),irangr(iranz), 
     ,              nord,sued,ost,west,kna,kne,ksa,kse,iea,iee,iwa,iwe 
c 
      common /urand/       zday(lrp),zvar(lrp),zalt(lrp)

      common /aux/ cc(kasor,9) 

      common uc(ndrei),vc(ndrei),stc(ndrei),avc(ndrei),z(m,n) 
      common zac(khor),wobc(khor),stuvc(khor),fricu(khor),cxc(khor) 
      common cyc(khor),pac(khor),txc(khor),tyc(khor) 
      common stpc(ndrei),sac(ndrei),tec(ndrei) 
      common pres(ilo),wc(ndrei),fricv(khor) 
c 
      common /dreh/ sinfu(m),sinfv(m),cosfu(m),cosfv(m),sincx(m), 
     ssincy(m),bx(m),by(m),pxu(m),pyv(m),pyu(m),pxv(m) 
c 
      common /gitter/ dt,r,g,dl,dlr,dlrh,dln(m),rdln(m),rad,dth, 
     ddlvo(m),dlvu(m),gh,rdt,dt2,r4,coru(m),corv(m),dtrdln(m),dtdlr 
c 
      common/num/dc(ilo),av(ilo),ad(ilo),dh(ilo),pd(ilo), 
     pprd(ilo),pr2d(ilo),r2d(ilo),tkw(ilo),tau(ilo),dd(ilo), 
     qqa(ilop1),qbet(ilop1),qn(ilop1),rd(ilop1) 
c 
      common /met/ windx,windy,wlam,pnull,rhoq(ilo),refrho(ilo),stress 
      parameter (stride = 4)
#ifdef MPI
      include 'C_mpi.f'
      include 'C_parallel.f'
#endif   
c
      nende = 500 !nosor/20 
      iter = 0
c
   77 continue
#ifndef MPI
      iter = iter+1
      cxcmax = 0.

c-----------------------------------------------------------------------
c red-black ordering for the 9-point stencil
c-----------------------------------------------------------------------
      do 200 kodd = 1,2
      do 190 iodd = 1,2
!CDIR VECTOR NODEP
      do 180 ikrb = 1,nrbsor(iodd,kodd)

c----------------------------------------------------------------------- 
c     index-administration........ 
c ----------------------------------------------------------------------

      ik = irbsor(ikrb,iodd,kodd)
      i = isor(ik) 
      k = ksor(ik) 
c
      om = omega(ik) 
      omc = 1.-om 
      dom = dble(om) 
      hz = dprod(omc,z(i,k)) + dom*( dble(stuvc(ik)) 
     1 +dprod(cc(ik,1),z(i-1,k-1))+dprod(cc(ik,2),z(i-1,k)) 
     2 +dprod(cc(ik,3),z(i-1,k+1)) 
     3 +dprod(cc(ik,4),z(i,k-1))+dprod(cc(ik,5),z(i,k+1)) 
     4 +dprod(cc(ik,6),z(i+1,k-1)) 
     5 +dprod(cc(ik,7),z(i+1,k))+dprod(cc(ik,8),z(i+1,k+1))) 

c-------------------------------------------------------------------------
c   error vector 
c-------------------------------------------------------------------------

      cxc(ik) = abs(z(i,k)-hz) 
      z(i,k) = hz 
  180 continue
  190 continue
  200 continue

c-------------------------------------------------------------------------
c evaluate the first norm of the deviation z(iter)-z(iter-1)
c-------------------------------------------------------------------------

!CDIR VECTOR
      do ik = 1,nosor
         cxcmax = max (cxcmax, abs (cxc(ik)))
      enddo

c###############check this please ????????????????

cDEBUG      res = 0.
cDEBUG      do ik = 1, nosor
cDEBUG         i = isor(ik) 
cDEBUG         k = ksor(ik) 
cDEBUG         tmp = z(i,k)
cDEBUG     1    -dprod(cc(ik,1),z(i-1,k-1))-dprod(cc(ik,2),z(i-1,k)) 
cDEBUG     2    -dprod(cc(ik,3),z(i-1,k+1)) 
cDEBUG     3    -dprod(cc(ik,4),z(i,k-1))-dprod(cc(ik,5),z(i,k+1)) 
cDEBUG     4    -dprod(cc(ik,6),z(i+1,k-1)) 
cDEBUG     5    -dprod(cc(ik,7),z(i+1,k))-dprod(cc(ik,8),z(i+1,k+1)) 
cDEBUG     6    -dble(stuvc(ik))
cDEBUG         res = res + tmp * tmp
cDEBUG      end do
cDEBUG      res = sqrt (res)
cDEBUG      print *, '   iter=', iter, ', max=', cxcmax, ', |res|=', res
c
c stopping criterion 
c
      if (cxcmax .gt. eps .and. iter .lt. nende) goto 77
c
c end of iterations
c
      if (iter .eq. nende) then
         write(*,*) iter,eps,cxcmax
         print *, 'sor: iterations exausted'
         stop 'in sor'
      else
c         print *, 'sor: converged on iter', iter
      end if
c 
cval      itel(ip) = last*1000+liter 
#else
       call cgsolver(cc,stuvc,z)
#endif 
c 
c-----------------------------------------------------------------------------
c    time averaging of surface elevation and multiplication with g and 
c    surface density
c    barotropic pressure gradients  
c     ps = luftdruckanomalie, bezogen auf pnull = 1013 mb 
c-----------------------------------------------------------------------------

      lwe = 0 
      ps = 0. 
c
#ifdef MPI
      do 75 k = J1,J2
      lwa = lb0(k)
      lwe = le0(k)
#else
      do 75 k = 1,n 
      lwa = lwe+1 
      lwe = indend(k) 
#endif
!CDIR VECTOR
      do 75 lw = lwa,lwe 
      i = iwet(lw) 
      ii = indwet(lw)
c      ps = pac(lw)-pnull 

c---------------------------------------------------------------------------------------
c    time weighted averaging with tmit/htmit 
c---------------------------------------------------------------------------------------  
      stuvc(lw) = ps + g*(zac(lw)+htmit*z(i,k))*(rhoq(1)+stc(ii+1)) 
c      write(373,*)stuvc(lw),zac(lw),z(i,k),stc(ii+1),lw
c-------------------------------------------------------------------- 
  75  continue
#ifdef MPI
      call INTI_COMM_VAR
      call com_horizontal(stuvc(1),khor,dtemp,iel,1,1,1,2,1,icdepb1,
     *idepb1,icdept1,
     *idept1,ncsendb1,nsendb1,ncsendt1,nsendt1,nbhlayer1,nthlayer1)
      call com_vertical(stuvc(1),khor,dtemp,iel,1,1,1,2,1,lazc)

      call MPI_Waitall(nreq1, ir_r_b, status, ierr)
      INUM3=0
      INUM4=0
      call HOR_UNPACK(stuvc(1),khor,dtemp,iel,1,1,nbhlayer1,nthlayer1,
     &idepb1,idept1,icdepb1,icdept1,1,2,iwet)
      INUM1=0
      INUM2=0
      call VER_UNPACK(stuvc(1),khor,dtemp,iel,1,1,1,2,lazc)
#endif
      return 
      end 



*��8
c-----------------------------------------------------------------------
      subroutine sorcof(htmit,htmitq) 
c-----------------------------------------------------------------------
c     belegung der koeffizientenmatrix fuer sor 
c     berechnung der rechten seite der zeta-gleichung 
c-----------------------------------------------------------------------
      include 'C_model.f'
      parameter(ilop1=ilo+1) 
      parameter(khor1=khor+1) 
c 
      common /ind/ iwet(khor1),ldep(khor),lazc(khor), 
     cindend(n),isornr(n),isorsr(n),islab(n) 
c 
      integer nord,sued,ost,west 
      integer*2 isor,ksor 
c 
      common /intit/ itel(1440),jc(m,3),ikk(kasor),isor(kasor), 
     cksor(kasor),nofrez,lone,itende,nfreez 
c 
      common /iterat/ omega(kasor),za(m,3),cx(m,3),cy(m,3),rhoz(m,3), 
     ctcc(8) 
c 
      common /uirand/ ilwrnd(lrp),islrnd(lrp),ldownr(lrp),irangr(iranz),
     ,              nord,sued,ost,west,kna,kne,ksa,kse,iea,iee,iwa,iwe 
c 
      common /urand/       zday(lrp),zvar(lrp),zalt(lrp)

c 
ccc   common /tide/ costom,sintom,cosinc,sininc,hcosg(lrp),hsing(lrp) 
c 
      common /aux/ cc(kasor,9) 
      common uc(ndrei),vc(ndrei),stc(ndrei),avc(ndrei),z(m,n) 
      common zac(khor),wobc(khor),stuvc(khor),fricu(khor),cxc(khor) 
      common cyc(khor),pac(khor),txc(khor),tyc(khor) 
      common stpc(ndrei),sac(ndrei),tec(ndrei) 
      common pres(ilo),wc(ndrei),fricv(khor) 
c 
      common /dreh/ sinfu(m),sinfv(m),cosfu(m),cosfv(m),sincx(m), 
     ssincy(m),bx(m),by(m),pxu(m),pyv(m),pyu(m),pxv(m) 
c 
      common /gitter/ dt,r,g,dl,dlr,dlrh,dln(m),rdln(m),rad,dth, 
     ddlvo(m),dlvu(m),gh,rdt,dt2,r4,coru(m),corv(m),dtrdln(m),dtdlr 
c 
      common/num/dc(ilo),av(ilo),ad(ilo),dh(ilo),pd(ilo), 
     pprd(ilo),pr2d(ilo),r2d(ilo),tkw(ilo),tau(ilo),dd(ilo), 
     qqa(ilop1),qbet(ilop1),qn(ilop1),rd(ilop1) 
c 
      common /met/ windx,windy,wlam,pnull,rhoq(ilo),refrho(ilo),stress 

c 
c     berechnung der rechten seite des gleichungssystems in zeta 
c     belegen der koeffizientenmatrix 
c 
#ifdef MPI
       include 'C_mpi.f'
       include 'C_index.f'
       common /cord/ I1,I2,J1,J2,IA1,IA2,JA1,JA2,IB1,IB2,JB1,JB2
     1,ID1,ID2,JD1,JD2
       common /icecor/ ICEI1,ICEI2,ICEJ1,ICEJ2,ICEIA1,ICEIA2,ICEJA1,
     * ICEJA2,ICEX1,ICEX2,ICEY1,ICEY2,ICEP1,ICEP2,ICEQ1,ICEQ2,ICEX3,
     * ICEX4,ICEY3,ICEY4,ICEX5,ICEX6,ICEY5,ICEY6
      common /isor/ isornr1(n),isorsr1(n),isornr2(n),isorsr2(n)
     1,ik1(khor),norsor1
      common /surface/ izets(m,n),izete(m,n)
#endif     
      nz = n-1 
c 
#ifdef MPI
      k1=ICEJA1
      do 172 k=2,3
      k1=k1+k-2 
      lwa=izets(IA1,k1)
      lwe=izete(IA2,k1)
      nwet=indwet(lwa)
      do 150  i=IA1,IA2
#else
      nwet = islab(1) 
      lwe = indend(1) 
      do 172 k = 2,3 
      lwa = lwe+1
      lwe = indend(k)
      do 150  i=1,m 
#endif 
      rhoz(i,k) = 0.0 
      za(i,k) = 0.0 
      cx(i,k) = 0.0 
      cy(i,k) = 0.0 
  150 jc(i,k) = 0 
      do 171 lw = lwa,lwe 
      i = iwet(lw) 
c     -------------------------- 
      cx(i,k) = cxc(lw) 
      cy(i,k) = cyc(lw) 
      za(i,k) = zac(lw) 
      jc(i,k) = 1 
      rhoz(i,k) = rhoq(1)+stc(nwet+1) 
  171 nwet = nwet+lazc(lw) 
  172 continue 
c     ----------------------- koeffizienten im slab kk 
      ik = 0 
      itt=0
#ifdef MPI
      it=0
      JL1=J1
      JL2=ICEJ2
      IF(myid .lt. nPh) JL1=3
      do 175 kk = JL1,JL2
      lwa=izets(IA1,kk+1)
      lwe=izete(IA2,kk+1)
      if(lwe .ge. lwa) then
      nwet=indwet(lwa)
      else
      nwet=0
      endif
      do 173 i = IA1,IA2
      do 1731 l = 2,3
#else
      do 175 kk = 3,nz 
      lwa = indend(kk)+1
      lwe = indend(kk+1)
      nwet = islab(kk) 
c     ----------------------- west-shift der slabs 
      do 173 i = 1,m 
      do 1731 l = 2,3 
#endif 
      cx(i,l-1) = cx(i,l) 
      cy(i,l-1) = cy(i,l) 
      za(i,l-1) = za(i,l) 
      jc(i,l-1) = jc(i,l) 
 1731 rhoz(i,l-1) = rhoz(i,l) 
      cx(i,3) = 0.0 
      cy(i,3) = 0.0 
      za(i,3) = 0.0 
      jc(i,3) = 0 
  173 rhoz(i,3) = 0.0 
c     ------------------------ slab kk+1 fuellen 
      do 174 lw = lwa,lwe 
      i = iwet(lw) 
c     ----------------------------- 
      cx(i,3) = cxc(lw) 
      cy(i,3) = cyc(lw) 
      za(i,3) = zac(lw) 
      jc(i,3) = 1 
      rhoz(i,3) = rhoq(1)+stc(nwet+1) 
  174 nwet = nwet+lazc(lw) 
c 
c 
c     ---------------------------------------------------- 
c     ----------------- loop fuer ik feuchte innere punkte 
c     ---------------------------------------------------- 
c 
c     ------------------ i-sor nordrand (nr) u. suedrand (sr) 
#ifndef MPI
      lwa = isornr (kk) 
      lwe = isorsr (kk) 
#else
      lwa = isornr1 (kk)
      lwe = isorsr1 (kk)
#endif 
c 
#ifdef MPI
      do 175 lw = lwa,lwe 
      itt=itt+1
      ik=ik1(itt)
#else
      do 175 lw = lwa,lwe
      ik=ik+1
#endif 
      i = iwet(lw) 
      gnw=4./(jc(i-1,2)+jc(i-1,1)+jc(i,2)+jc(i,1)) 
      gne=4./(jc(i-1,2)+jc(i-1,3)+jc(i,2)+jc(i,3)) 
      gse=4./(jc(i+1,2)+jc(i+1,3)+jc(i,2)+jc(i,3)) 
      gsw=4./(jc(i+1,2)+jc(i+1,1)+jc(i,2)+jc(i,1)) 
c     ------------------------------------- ost - koeffizienten 
      ce = sincx(i)*cx(i,2) 
      cee = ce*bx(i) 
c     ------------------------------------ west - koeffizienten 
      cw = sincx(i)*cx(i,1) 
      cww = cw*bx(i) 
c     ------------------------------------ nord - koeffizienten 
      cn = dlvo(i)*sincy(i-1)*cy(i-1,2) 
      cnn = cn*by(i-1) 
c     ------------------------------------ sued - koeffizienten 
      cs = dlvu(i)*sincy(i)*cy(i,2) 
      css = cs*by(i) 
      ccc = rhoz(i,2)* 
     f (gnw*(cnn+cww)-gsw*(cww+css)-gne*(cee+cnn)+gse*(cee+css) 
     f +ce+cs+cw+cn) 
c     -------- koeffizienten f. acht umliegende druckpunkte 
      tcc(1) =   -rhoz(i-1,1)*(cww+cnn)*gnw 
      tcc(2) =   -rhoz(i-1,2)*(gnw*(cww+cnn)-gne*(cee+cnn)-cn) 
      tcc(3) =    rhoz(i-1,3)*(cee+cnn)*gne 
      tcc(4) =   -rhoz(i,1)  *(gnw*(cnn+cww)-gsw*(cww+css)-cw) 
      tcc(5) =   -rhoz(i,3)  *(gse*(cee+css)-gne*(cee+cnn)-ce) 
      tcc(6) =    rhoz(i+1,1)*(cww+css)*gsw 
      tcc(7) =   -rhoz(i+1,2)*(gse*(cee+css)-gsw*(css+cww)-cs) 
      tcc(8) =   -rhoz(i+1,3)*(cee+css)*gse 
c     ........................................... 
      cmc = 1./(1.+ htmitq * ccc) 
      cmtc = htmitq * cmc 
c     rublid..................................... 
      stuvc(ik) = (dble(stuvc(lw))+htmit*(-dprod(ccc,za(i,2)) 
     ++dprod(tcc(1),za(i-1,1))+dprod(tcc(2),za(i-1,2)) 
     ++dprod(tcc(3),za(i-1,3))+dprod(tcc(4),za(i,1)) 
     ++dprod(tcc(5),za(i,3))+dprod(tcc(6),za(i+1,1)) 
     ++dprod(tcc(7),za(i+1,2))+dprod(tcc(8),za(i+1,3))))*cmc 
c     ........................................... 
c     ........................................... 
      do 1175 l = 1,8 
 1175 cc(ik,l) = tcc(l)*cmtc 
  175 continue 
c 
c     ------------ start fuer sor aus w(1)*dt berechnen ------ 
c 
#ifdef MPI
      JL1=ICEJA1
      JL2=ICEJA2
      IF(myid .lt. nPh)JL1=3
      do 176 k = JL1,JL2
      lwa = isornr2 (k)
      lwe = isorsr2 (k)
#else
      do 176 k = 3,nz 
      lwa = isornr (k) 
      lwe = isorsr (k) 
#endif 
      do 176 lw = lwa,lwe 
      i = iwet(lw) 
  176 z(i,k) =  dt*wobc(lw) 

      return 
      end 

*��9 konti
      subroutine konti(ht,dzbod,wtest,wsur)
c-----------------------------------------------------------------------
c     solving of continuity equation, calcualting w-component  
c-----------------------------------------------------------------------
      include 'C_model.f'

      parameter(ilop1=ilo+1) 
      parameter(khor1=khor+1) 

      double precision wo,wsur(khor)
      real ht(m,n),dzbod(m,n),wtest(m,n)
 
      common /ind/ iwet(khor1),ldep(khor),lazc(khor), 
     cindend(n),isornr(n),isorsr(n),islab(n) 
 
      common uc(ndrei),vc(ndrei),stc(ndrei),avc(ndrei),z(m,n) 
      common zac(khor),wobc(khor),stuvc(khor),fricu(khor),cxc(khor) 
      common cyc(khor),pac(khor),txc(khor),tyc(khor) 
      common stpc(ndrei),sac(ndrei),tec(ndrei) 
      common pres(ilo),wc(ndrei),fricv(khor) 
 
      common /dreh/ sinfu(m),sinfv(m),cosfu(m),cosfv(m),sincx(m), 
     ssincy(m),bx(m),by(m),pxu(m),pyv(m),pyu(m),pxv(m) 
 
      common /gitter/ dt,r,g,dl,dlr,dlrh,dln(m),rdln(m),rad,dth, 
     ddlvo(m),dlvu(m),gh,rdt,dt2,r4,coru(m),corv(m),dtrdln(m),dtdlr 
 
      common/num/dc(ilo),av(ilo),ad(ilo),dh(ilo),pd(ilo), 
     pprd(ilo),pr2d(ilo),r2d(ilo),tkw(ilo),tau(ilo),dd(ilo), 
     qqa(ilop1),qbet(ilop1),qn(ilop1),rd(ilop1) 
 
      common /met/ windx,windy,wlam,pnull,rhoq(ilo),refrho(ilo),stress 
 
      common /aux/ u(ilo,m,2),v(ilo,m) 
#ifdef MPI
      include 'C_index.f'
      include 'C_mpi.f'
      common /cord/ I1,I2,J1,J2,IA1,IA2,JA1,JA2,IB1,IB2,JB1,JB2
     1,ID1,ID2,JD1,JD2
      common /icecor/ ICEI1,ICEI2,ICEJ1,ICEJ2,ICEIA1,ICEIA2,ICEJA1,
     * ICEJA2,ICEX1,ICEX2,ICEY1,ICEY2,ICEP1,ICEP2,ICEQ1,ICEQ2,ICEX3,
     * ICEX4,ICEY3,ICEY4,ICEX5,ICEX6,ICEY5,ICEY6
      common /surface/ izets(m,n),izete(m,n)
#endif
      mz = m-1 
      nz = n-1 
      u=0
c 
c     prepare moving  arrays for  kk = 1 
c 
#ifndef MPI
      nwet = 0 
      lwa = 1 
      lwe = indend(1) 
#else
      lwa=izets(IB1,JB1)
      lwe=izete(IB2,JB1)
      nwet=indwet(lwa)
#endif  

      do 4 lw = lwa,lwe 
      i = iwet(lw) 
      ldown = lazc(lw) 
      do 3 j = 1,ldown 
      u(j,i,2) = uc(nwet+j) 
   3  continue 
      nwet = nwet+ldown 
   4  continue 
 
#ifdef MPI
      do 76 kk = ICEJA1,ICEJA2
      lwa=izets(IB1,kk)
      lwe=izete(IB2,kk)
      if(lwe .ge. lwa) then
      nwet=indwet(lwa)
      else
      nwet=0
      endif
!      lwa = indend(kk-1)+1
!      lwe = indend(kk)
!      nwet = islab(kk-1)
      do 16 i  = IB1,IB2
#else 
      do 76 kk = 2,nz 
      lwa = indend(kk-1)+1
      lwe = indend(kk)
      nwet = islab(kk-1)
c     --------------------- shifting slabs west
      do 16 i = 1,m 
#endif
      do 16 j = 1,ilo 
      u(j,i,1) = u(j,i,2) 
      u(j,i,2) = 0.0 
      v(j,i) = 0.0 
   16 continue 
c 
c     --------------------- fill w for  slab kk 
      do 76 lw = lwa,lwe 
      i = iwet(lw) 
      im1 = max0(i-1,1) 
      ldown = lazc(lw) 
      dvo = dlvo(i)*dlr 
      dvu = dlvu(i)*dlr 
      wo = 0.d0
c     ---------------- up-down-loop ----------------------- 
c     ---------------------------- (u,v,w  up,  rest down)- 
      do 97 jj = ldown,1,-1 
      u(jj,i,2) = uc (nwet+jj) 
      v(jj,i) = vc (nwet+jj) 
c     ------------------ continuity equation for  k + 1 
      wo = wo+dprod(rdln(i),u(jj,i,1))-dprod(rdln(i),u(jj,i,2)) 
     ++dprod(dvu,v(jj,i))-dprod(dvo,v(jj,im1)) 
   97 wc(nwet+jj) = wo 
c     ------- w at sea surface (compressed) -------- 
      wobc(lw) = wo
      wsur(lw) = wo
      nwet = nwet+ldown 
   76 continue 


      return 
      end 




*Corrections by Irina
*Line ~220:  if(lazc(lw).eq.1) zzz=dfloat(ldep)+dble(zac(lw))   !Irina  Correction#1
*Line ~310:c      if (j.lt.n-2 .and. i.gt.5 .and. i.lt.m-2) then  !Irina  Correction#2 (now commented)
*Line ~442:*----- add the joint river and P-E fluxes of the characteristics ibio=1,nbio-----Correction#4
*Line ~569:c      end if !if (j.lt.n-2 .and. i.gt.5 .and. i.lt.m-2) then  !Irina  Correction#2 (now commented)


c      subroutine stromTVD(zinc,vtmit,vtmic,avmax,ftac,szahl,
c     &sh_wave,nday,time,wgesch)
      subroutine stromTVD(zinc,vtmit,vtmic,avmax,ftac,szahl)
c-----------------------------------------------------------------------
c     transport of salinity and temperature
c     TVD total variation diminishing scheme with superbee
c     
c-----------------------------------------------------------------------
#ifdef _FABM_
      use hamsom_fabm , only: vertical_movement
      include 'C_nbio_fabm.f'
#else
      include 'C_nbio.f'
#endif
      include 'C_model.f'
      include 'C_index.f'
      
      parameter(ilop1=ilo+1) 
      parameter(khor1=khor+1)
      parameter(nbio2=nbio) !take fish out of advection 
c      include 'C_rbio.f'
c      include 'C_bio.f'
      common/T_4d/ Tc (ndrei,nbio)
c      common/T_4d/ Tc (ndrei,nbio),Tc_mit (ndrei,nbio+nsed) 
c      common /T_4d/ Tc(ndrei,nbio),Tc_mit(ndrei,nbio+nsed)
c     & ,Tsed(ndrei,nsed)
      dimension Tc_flux_all(khor,nbio)  !river and P-E ![m**s/sec*timestep]

c parameters used for the TVD procedure:
      integer NSIG
c      double precision COURUL, COURUR, COURVO, COURVU, COURWO, COURWU
c      double precision D1, D2, DDI, DDIZ, DENUM
c      double precision EPSIL
c      double precision PHIR, PHIR1
      double precision RATIO, RSUL, RSUR, RSVO, RSVU, RSWO, RSWU
      double precision SADV
      
      double precision dwdz, divcour
      double precision zzz
      double precision  s,sn,sd,sbet,tsa,tsb,tsc,tsal 
      double precision tsbb

      double precision  XYadv,svert

c      double precision uot,uwt,vnt,vst
c      double precision  wo, wu
     
      dimension XYadv(ndrei,nbio),svert(ndrei,nbio)
      dimension Tl (m,ilo,nbio,3),Tc_k(ilo,nbio),Tl_k(ilo,nbio)
      dimension Tl_bio(nbio)
      dimension s(ndrei,nbio),sn(ndrei,nbio),sd(ndrei,nbio)
      dimension sbet(ndrei,nbio),ftac(khor),szahl(ndrei) 
      dimension tsa(ndrei,nbio),tsb(ndrei,nbio),tsc(ndrei,nbio) 
      dimension tsal(ndrei,nbio)
      dimension wgesch(m,n)
      double precision zinc(m,n)
      logical obwas,advekt,rfals 
c 
c************************  neue vertikaladvektion **********************
       dimension e(ilo),taa(ilo),tac(ilo),tab(ilo),sad(ilo),tad(ilo)
c***********************************************************************     
      common /ind/ iwet(khor1),ldep(khor),lazc(khor), 
     *indend(n),isornr(n),isorsr(n),islab(n) 
c 
      common uc(ndrei),vc(ndrei),stc(ndrei),avc(ndrei),z(m,n) 
      common zac(khor),wobc(khor),stuvc(khor),fricu(khor),cxc(khor) 
      common cyc(khor),pac(khor),txc(khor),tyc(khor) 
      common stpc(ndrei),sac(ndrei),tec(ndrei) 
      common pres(ilo),wc(ndrei),fricv(khor) 
c 
      common /dreh/ sinfu(m),sinfv(m),cosfu(m),cosfv(m),sincx(m), 
     *sincy(m),bx(m),by(m),pxu(m),pyv(m),pyu(m),pxv(m) 
c 
      common /gitter/ dt,r,g,dl,dlr,dlrh,dln(m),rdln(m),rad,dth, 
     *dlvo(m),dlvu(m),gh,rdt,dt2,r4,coru(m),corv(m),dtrdln(m),dtdlr 
      common /met/ windx,windy,wlam,pnull,rhoq(ilo),refrho(ilo),stress 
c 
c     ------------ laenge von common /aux/ = m*(ilo*3*10+12) 
      common /aux/ u(ilo,m,2),v(ilo,m,2),sa(ilo,m,3),te(ilo,m,3) 
     *,za(m,2),su(ilo,m,2),sv(ilo,m),wd(ilo,m,2)
      dimension avd(ndrei,nbio) 
c     ---------------------- laenge von common /auxint/ = m*(ilo*6+3) 
      common /auxint/ jc(ilo,m,3),lay(ilo,m,2),laz(m,2) 
c     new array
      integer jc5(m,n,ilo)
c      common /rivflow/ rivflux(m,n)
 
      common/num/dz(ilo),av(ilo),ah(ilo),dh(ilo),pd(ilo), 
     *pprd(ilo),pr2d(ilo),r2d(ilo),tkw(ilo),tau(ilo),d_(ilo), 
     *qa(ilop1),qbet(ilop1),qn(ilop1),rd(ilop1) 
      logical init
      common /vec44/ init,ljumm(khor),vol0,dtdx(m),dvo(m),dvu(m),
     1 lww(-2:2,-2:2,ndrei),llay(0:khor,ilo)
      data init /.true./
      dimension d(ilo,m),Tc0(0:ndrei,nbio),dd(0:ndrei),r_dt_dd(ndrei)
      dimension dd_dd(ndrei)
      dimension Dvoltj(m,ilo,n)
#ifdef MPI
      include 'C_mpi.f'
      include 'C_parallel.f'
#endif      
c--------  set biology parameters ----------
      ibio0=1  ! ibio= 3,nbio
c      call Rzero(Bmass,nbio+1) !intagrated biomass, output each time step
      lob = 1

      mz = m-1 
      nz = n-1 

ccc  epsilon ist die minimale schichtdicke      
      epsilon=0.5
CCC  EPSIL is a minimum value used by the TVD procedure
      if (init) then
#ifdef MPI
	do I=1,khorl
	lw=lzet(I)
#else    
	do lw = 1,khor
#endif
            ljumm(lw) = 0
         end do

#ifdef MPI
         do j=icej1,icej2
          do lw = lb0(j),le0(j)
#else
        do j=2,nz
          do lw = lb(j),le(j)
#endif
               i = iwet(lw) 
               ljum = min0(i-2,1) * min0(j-2,1) 
               ljum = max0(0,ljum) 
               ljumm(lw) = ljum
            end do
         end do

         dd(0)=0.
         vol0=0.
	   nwet=0

c#ifdef MPI
c      do j=J1,J2
c         do lw=lb0(j),le0(j)
c#else  
      do j=1,n
         do lw=lb(j),le(j)
c#endif
               i = iwet(lw) 
               ldown = lazc(lw) 
               do k = 1,ldown 
	         nwet=nwet+1
                  if (k.ne.1) then
                     if(k.eq.ldown) then
c                        zzz=float(ldep(lw))
                     zzz=dfloat(ldep(lw))
                     else
                        zzz=dble(pd(k))
                     end if
c                     Dvoltj(i,k,j)=zzz*dl*dln(i)	![m**3]
                     Dvoltj(i,k,j)=zzz*dprod(dl,dln(i))
                     vol0=vol0+Dvoltj(i,k,j)          
                     dd(nwet)=zzz
	               dd_dd(nwet)=dprod(dd(nwet),dd(nwet))
                  endif
               end do
            end do
         end do

#ifdef MPI
      do i=I1,I2
#else  
      do i=1,m
#endif
            dtdx(i) = dt/  dln(i)		!dln(j)  = dlam*cosfiu 
            dvu(i)=dln(i)*dlvu(i)
            dvo(i)=dln(i)*dlvo(i)
         end do
            dt_lr=dt*dlr

#ifdef MPI
      do j=J1,J2
      do i=I1,I2

#else  
      do i=1,m
      do j=1,n
#endif
           do k=1,ilo
         jc5(i,j,k)=0
           enddo
          enddo
         enddo

      nwet = 0 
      lwe = 0

      do j = 1,n 
         lwa = lwe+1 
         lwe = indend(j) 
         do lw = lwa,lwe 
          i = iwet(lw)
          j5= jwet(lw)
            do k = 1, lazc(lw)
            jc5(i,j5,k)=nwet+k
            end do
            nwet = nwet+lazc(lw) 
         end do
      end do




       do ii = -2,2
       do jj = -2,2
#ifdef MPI

      do ij=JB1,JB2
      lwa=lb2(ij)
      lwe=le2(ij)
      if(lwe .ge. lwa) then
      nwet1=indwet(lwa)+1
      nwet2=indwet(lwe)+lazc(lwe)
      do  nwet=nwet1,nwet2
           lw=llw(nwet)
c       since there are 3 additional land lines at west and north, no limits
c       for i=iwet(..)-2 and j=jwet(..)-2 :
        i = min(m,iwet(lw)+ii)
        j = min(n,jwet(lw)+jj)
        k=indver(nwet)
        lww(ii,jj,nwet) = jc5(i,j,k)
      enddo
      endif
      enddo           
#else
	do nwet = 1,ndrei
            lw=llw(nwet)
c       since there are 3 additional land lines at west and north, no limits
c       for i=iwet(..)-2 and j=jwet(..)-2 :
        i = min(m,iwet(lw)+ii)
        j = min(n,jwet(lw)+jj)
        k=indver(nwet)
        lww(ii,jj,nwet) = jc5(i,j,k)
        end do
#endif
        enddo
        enddo



         do k = 1,ilo
            llay(0,k) = 0
         end do
#ifdef MPI
         do j=icej1,icej2
          do lw = lb0(j),le0(j)
#else
         do j=2,nz
          do lw = lb(j),le(j)
#endif
           if (ljumm(lw).eq.1) then
            nwet = indwet(lw)
            ldown = lazc(lw) 
            ld  = ldown-1 
            lein= min0(1,ld) 
            if (lein.eq.1) then
               do k = ld,lob,-1 
                  llay(0,k) = llay(0,k)+1
                  llay(llay(0,k),k) = nwet+k
               end do
            end if
           end if
          end do
         end do

         init = .true.
      end if
c





*------------Volume of the 1st layer---------------
      volbio=0.
      nwet=0
#ifdef MPI
      do j = JB1,JB2
      do lw = lb2(j),le2(j)
#else
      do j = 1,n 
      do lw = lb(j),le(j) 
#endif 
       nwet=indwet(lw)
            i = iwet(lw) 
                        zzz=dble(pd(1))+dble(zac(lw))
      if(lazc(lw).eq.1) zzz=dfloat(ldep(lw))+dble(zac(lw))   !Irina  Correction#1
           dd(nwet+1) =  max(epsilon,zzz)
c           Dvoltj(i,1,j)=zzz*dprod(dl,dln(i))
           Dvoltj(i,1,j)=dd(nwet+1)*dprod(dl,dln(i))
           volbio=volbio+Dvoltj(i,1,j)          
           dd_dd(nwet+1)=dprod(dd(nwet+1),dd(nwet+1))
      end do
      end do
#ifdef MPI
      volbio=0.
      volbio1=0.
      do j = J1,J2
      do lw = lb0(j),le0(j)
      i = iwet(lw)
      volbio1=volbio1+Dvoltj(i,1,j)
	 enddo
      enddo
      CALL MPI_ALLREDUCE(volbio1,volbio,1,MPI_REAL,
     & MPI_SUM,MPI_COMM_WORLD,IERR)
       volbio = volbio+vol0
#else      
      volbio = volbio+vol0
#endif
c 
       wmax=0.
       umax=0.
       icount=0
       admax=avmax 
       dtdt = 2.*dt
ccc  epsilon ist die minimale schichtdicke      
c      epsilon=0.5
CCC  EPSIL is a minimum value used by the TVD procedure
      EPSIL=1.0d-10
c 
      avmin = 0.134e-6 
      stfak = 1./1.35 

*-------- T,S ----------------
#ifdef MPI
      do 51 J=JB1,JB2
      lwa=lb2(J)
      lwe=le2(J)
      if(lwe .ge. lwa) then
      nwet1=indwet(lwa)+1
      nwet2=indwet(lwe)+lazc(lwe)
      do  k=nwet1,nwet2
         Tc(k,1)=tec(k)
         Tc(k,2)=sac(k)
         r_dt_dd(k)=dt/dd(k)
      enddo
      endif
 51   continue
 
      do 52 ibio = 1,nbio
      do 52 J=JB1,JB2
      lwa=lb2(J)
      lwe=le2(J)
      if(lwe .ge. lwa) then
      nwet1=indwet(lwa)+1
      nwet2=indwet(lwe)+lazc(lwe)
      do k=nwet1,nwet2
      Tc0(k,ibio) = Tc(k,ibio)
      enddo
      endif
 52   continue
#else
      do nwet=1,ndrei
         Tc(nwet,1)=tec(nwet)
         Tc(nwet,2)=sac(nwet)
         r_dt_dd(nwet)=dt/dd(nwet)
      enddo
*-----------------------------
      do ibio = 1,nbio
         do ll = 1,ndrei
            Tc0(ll,ibio) = Tc(ll,ibio)
         end do
      end do
#endif

      do ibio = 1,nbio
         Tc0(0,ibio) =0. !-100.
      end do
c--------------------------------------------


c 
*------------------- transport cycle -------------------------------------
#ifdef MPI
      do j = icej1,icej2
      lwa = lb0(j)
      lwe = le0(j)
#else
      do j = 2,nz 
      lwa = lb(j)
      lwe = le(j)
#endif
      if (lwa.le.lwe) then
         llb = indwet(lwa)+1
         lle = indwet(lwe)+lazc(lwe)
         do ll = llb,lle
            i = iwet(llw(ll))
            k = ll-indwet(llw(ll))

            if (ljumm(llw(ll)).eq.1) then
c      if (j.lt.n-2 .and. i.gt.5 .and. i.lt.m-2) then  !Irina  Correction#2 (now commented)
*------------ horizontal Advektion terms, TVD -----------

        uo = uc(ll)*min(lww(0,1,ll),1)
     &  /( 0.5* ( dd(ll)+dd(lww(0,1,ll)) ) )
        uot=uc(ll)* min(lww(0,1,ll),1)

       COURUR=abs(dtdx(i)*uo)

       lu = lww(0,-1,ll)
       if (lu.eq.0) then
          uw = 0.
          uwt= 0.0d+0
       else
          uw = uc(lu)/( 0.5* ( dd(ll)+dd(lww(0,-1,ll)) ) )
          uwt= uc(lu)
       end if
     
        COURUL=abs(dtdx(i)*uw)

        lv = lww(-1,0,ll)
        if (lv.eq.0) then
           vn = 0.
           vnt= 0.0d+0
        else
           vn=vc(lv)/( 0.5* ( dd(ll)+dd(lww(-1,0,ll)) ) )
           vnt= vc(lv)
        end if

         COURVO=abs(dt_lr*vn)

        vs=vc(ll)* min(lww(1,0,ll),1)
     &  /( 0.5* ( dd(ll)+dd(lww(1,0,ll)) ) )
        vst=vc(ll)*min(lww(1,0,ll),1)
 
         COURVU=abs(dt_lr*vs)
*------------ vertical advection terms ------------------------
c     SinkD_o =-BioC(23)*r_dt_dd(ll)    !BioC(23) [m/s]
c     SinkD_u = SinkD_o

      if (k.eq.1) then
         llo = ll
      else
         llo = ll-1
      endif
        
      if (k.eq.lazc(llw(ll))) then
         llu = ll
      else
         llu = ll+1
      endif
		
      wo = wc(ll)	      !	 k
      wu = wc(llu)	  !	 k+1

c      if (k.eq.1) then
c	      wo = 0.
c	      SinkD_o = 0.
c      endif
        
c      if (k.eq.lazc(llw(ll))) then
c	      wu=0.
c	      SinkD_u=0.
c      endif

      dwdz=(wo-wu)*r_dt_dd(ll)


      do ibio= ibio0,nbio2
      RSUL=0.0d+0
      do irsu=1,min0(jc5(i,j-1,k),1)
      DENUM=Tc0(ll,ibio)-Tc0(lww(0,-1,ll),ibio) 
      if (abs(DENUM).lt.EPSIL) THEN
       PHIR=0.0d+0
      else
       if (uwt.ge.0.0d+0) NSIG=-1
       if (uwt.lt.0.0d+0) NSIG=1
       if ((NSIG.eq.-1 .and. jc5(i,j-2,k).eq.0) .or. (NSIG.eq.1 .and.
     & jc5(i,j+1,k).eq.0)) then
       PHIR=0.
       else
       RATIO=(Tc0(lww(0,NSIG,ll),ibio)-Tc0(lww(0,NSIG-1,ll),ibio))/DENUM
       PHIR1=max(min(1.,2.*RATIO),min(2.,RATIO))
	 PHIR=max(0.,PHIR1)
       endif
      endif
       D1=PHIR*DENUM
       D2=0.5*abs(uwt)*(1.-COURUL)*D1
       RSUL=D2+0.5*(dprod((uwt+abs(uwt)),Tc0(lww(0,-1,ll),ibio))+
     &              dprod((uwt-abs(uwt)),Tc0(ll,ibio)))
       RSUL=RSUL*dble(dl)
       enddo

c east
      RSUR=0.0d+0
      do irsu=1,min0(jc5(i,j+1,k),1)
      DENUM=(Tc0(lww(0,1,ll),ibio)-Tc0(ll,ibio) )
      if (abs(DENUM).lt.EPSIL) THEN
       PHIR=0.0d+0
      else
       if (uot.ge.0.0d+0) NSIG=-1
       if (uot.lt.0.0d+0) NSIG=1
       if ((NSIG.eq.-1 .and. jc5(i,j-1,k).eq.0) .or. (NSIG.eq.1 .and.
     & jc5(i,j+2,k).eq.0)) then
       PHIR=0.0d+0
       else
c
       RATIO=(Tc0(lww(0,NSIG+1,ll),ibio)-Tc0(lww(0,NSIG,ll),ibio))
     & /DENUM
       PHIR1=max(min(1.,2.*RATIO),min(2.,RATIO))
	 PHIR=max(0.,PHIR1)
       endif
      endif
       D1=PHIR*DENUM
       D2=0.5*abs(uot)*(1.-COURUR)*D1
       RSUR=D2+0.5*(dprod((uot+abs(uot)),Tc0(ll,ibio))+
     &              dprod((uot-abs(uot)),Tc0(lww(0,1,ll),ibio)))
       RSUR=RSUR*dble(dl)
       enddo
c north
      RSVO=0.0d+0
      do irsu=1,min0(jc5(i-1,j,k),1)
      DENUM=(Tc0(lww(-1,0,ll),ibio)-Tc0(ll,ibio) )
      if (abs(DENUM).lt.EPSIL) THEN
       PHIR=0.
      else
       if (vnt.ge.0.0d+0) NSIG=1
       if (vnt.lt.0.0d+0) NSIG=-1
       if ((NSIG.eq.-1 .and. jc5(i-2,j,k).eq.0) .or. (NSIG.eq.1 .and.
     & jc5(i+1,j,k).eq.0)) then
       PHIR=0.
       else
       RATIO=(Tc0(lww(NSIG-1,0,ll),ibio)-Tc0(lww(NSIG,0,ll),ibio))
     & /DENUM
       PHIR1=max(min(1.,2.*RATIO),min(2.,RATIO))
	 PHIR=max(0.,PHIR1)
       endif
      endif
       D1=PHIR*DENUM
       D2=0.5*abs(vnt)*(1.-COURVO)*D1
       RSVO=D2+0.5*(dprod((vnt+abs(vnt)),Tc0(ll,ibio))+
     &              dprod((vnt-abs(vnt)),Tc0(lww(-1,0,ll),ibio)))
       RSVO=RSVO*dvo(i)
       enddo
c south
      RSVU=0.0d+0
       do irsu=1,min0(jc5(i+1,j,k),1)
      DENUM=(Tc0(ll,ibio)-Tc0(lww(1,0,ll),ibio))
      if (abs(DENUM).lt.EPSIL) THEN
       PHIR=0.
      else
       if (vst.ge.0.0d+0) NSIG=1
       if (vst.lt.0.0d+0) NSIG=-1
      if ((NSIG.eq.-1 .and. jc5(i-1,j,k).eq.0) .or. (NSIG.eq.1 .and.
     & jc5(i+2,j,k).eq.0)) then
       PHIR=0.
       else
       RATIO=(Tc0(lww(NSIG,0,ll),ibio)-Tc0(lww(NSIG+1,0,ll),ibio))/DENUM
       PHIR1=max(min(1.,2.*RATIO),min(2.,RATIO))
	 PHIR=max(0.,PHIR1)
       endif
      endif
       D1=PHIR*DENUM
       D2=0.5*abs(vst)*(1.-COURVU)*D1
       RSVU=D2+0.5*(dprod((vst-abs(vst)),Tc0(ll,ibio))+
     &              dprod((vst+abs(vst)),Tc0(lww(1,0,ll),ibio)))
       RSVU=RSVU*dvu(i)
       enddo

*-------- the horizontal advection part to the new value -------------------
      XYadv(ll,ibio)=(RSUL-RSUR+RSVU-RSVO)*dt/Dvoltj(i,k,j)

c#ifndef ECO2
*Correction#4
*----- add the joint river and P-E fluxes of the characteristics ibio=1,nbio--------
c       if (k.eq.1) then
c       XYadv(ll,ibio)=XYadv(ll,ibio)+Tc_flux_all(llw(ll),ibio)/dd(ll)
c       endif

c#endif

        SinkD_o=0.
        SinkD_u=0.
        if(ibio.gt.2)then !biological sinking velocities 
      if (k.eq.1) then
       wo = 0.
        SinkD_o = 0.
!KB       SinkD_u =vertical_movement(i,j,k+1,ibio-2)*r_dt_dd(ll+1)    !BioC(23) [m/s]
      elseif (k.eq.lazc(llw(ll))) then
       wu=0.
       SinkD_u=0.
       SinkD_o =vertical_movement(i,j,k,ibio-2)*r_dt_dd(ll)    !BioC(23) [m/s]
      else
       SinkD_o =vertical_movement(i,j,k,ibio-2)*r_dt_dd(ll)    !BioC(23) [m/s]
       SinkD_u =vertical_movement(i,j,k+1,ibio-2)*r_dt_dd(ll+1)    !BioC(23) [m/s]

       endif
         endif
         D_o2=SinkD_o		 
         D_u2=SinkD_u


c     TVD scheme
      SADV=0.0d+0
c     from above
      RSWO=0.0d+0
       if (k.gt.1) then
      DDI=(dd(ll)+dd(llo))*0.5
      COURWO=abs(wo)*dt/DDI
      DENUM=( Tc0(llo,ibio)-Tc0(ll,ibio))/DDI
      if (abs(DENUM).lt.EPSIL) THEN
       PHIR=0.
      else
      if (wo.ge.0.) NSIG=1
      if (wo.lt.0.) NSIG=-1
      if ((k.eq.lazc(llw(ll)) .and. NSIG.eq.1) .or. (k.eq.2 .and. NSIG
     &         .eq.-1)) then
      PHIR=0.
      else
      DDIZ=(dd(ll+NSIG)+dd(llo+NSIG))*0.5
      RATIO=( Tc0(llo+NSIG,ibio)-Tc0(ll+NSIG,ibio))/(DDIZ*DENUM)
       PHIR1=max(min(1.,2.*RATIO),min(2.,RATIO))
	 PHIR=max(0.,PHIR1)
      endif
      endif
      D1=PHIR*DENUM*DDI
      D2=0.5*abs(wo)*(1.-COURWO)*D1
      RSWO=D2+0.5*(dprod((wo+abs(wo)),Tc0(ll,ibio))+
     &             dprod((wo-abs(wo)),Tc0(llo,ibio)))
      endif
c     from below
      RSWU=0.0d+0
       if (k.lt.lazc(llw(ll))) then
      DDI=(dd(ll)+dd(llu))*0.5 
      COURWU=abs(wu)*(dt)/DDI
      DENUM=( Tc0(ll,ibio)-Tc0(llu,ibio))/DDI
      if (abs(DENUM).lt.EPSIL) THEN
       PHIR=0.
      else
      if (wu.ge.0.) NSIG=1
      if (wu.lt.0.) NSIG=-1
      if ((k.eq.lazc(llw(ll))-1 .and. NSIG.eq.1) .or. (k.eq.1 .and. NSIG
     &         .eq.-1)) then
      PHIR=0.
      else
      DDIZ=(dd(ll+NSIG)+dd(llu+NSIG))*0.5
      RATIO=( Tc0(ll+NSIG,ibio)-Tc0(llu+NSIG,ibio))/(DDIZ*DENUM)
       PHIR1=max(min(1.,2.*RATIO),min(2.,RATIO))
	 PHIR=max(0.,PHIR1)
      endif
      endif
      D1=PHIR*DENUM*DDI
      D2=0.5*abs(wu)*(1.-COURWU)*D1
      RSWU=D2+0.5*(dprod((wu+abs(wu)),Tc0(llu,ibio))+
     &             dprod((wu-abs(wu)),Tc0(ll,ibio)))
      endif
c vertical advection contribution

      SADV=(RSWU-RSWO)*r_dt_dd(ll)


cc----sinking vel is added by upstream---  NEW
	 w_up   = max(0.,D_o2*Tc0(ll,ibio))  +min(0.,D_o2*Tc0(llo,ibio))
	 w_down = max(0.,D_u2*Tc0(llu,ibio)) +min(0.,D_u2*Tc0(ll,ibio))
cc---------------------------------------


      svert(ll,ibio)= SADV
     &       -w_up+w_down           !----sinking---
      if (k.eq.1) then
      svert(ll,ibio) = svert(ll,ibio) - zinc(i,j)*Tc0(ll,ibio)
     & /dd(ll)

      endif ! if (k.eq.1) then
      end do      ! ib = 3, nbio  

      endif !if (ljumm(llw(ll)).eq.1) then

c            end if !if (j.lt.n-2 .and. i.gt.5 .and. i.lt.m-2) then  !Irina  Correction#2 (now commented)
         end do  !i
      end if
c#ifdef MPI
      end do     !j
c#endif

222   format(1x,3(e9.3,2x))
c 
c -----------vertikal diffusion ---------------- 
c 
#ifdef MPI
      do ibio=ibio0,nbio2
      do j = icej1,icej2 
      lwa = lb0(j)
      lwe = le0(j)
#else
      do ibio=ibio0,nbio2
      do j = 2,nz 
      lwa = lb(j)
      lwe = le(j)
#endif
      if (lwa.le.lwe) then
         llb = indwet(lwa)+1
         lle = indwet(lwe)+lazc(lwe)
         do ll = llb,lle
            if (ljumm(llw(ll)).eq.1) then
                  salref = Tc0(ll,ibio)
                  s (ll,ibio)  = Tc0(ll,ibio)
                  sn(ll,ibio)  = Tc0(ll,ibio)
                  stfak=szahl(ll)
                  avd(ll,ibio) = amin1(stfak*avc(ll),admax)
                  if(avd(ll,ibio).lt.0.) avd(ll,ibio)=admax 
            end if
         end do
      end if
      end do
      end do
#ifdef MPI
      do ibio=ibio0,nbio2
      do j = icej1,icej2 
      lwa = lb0(j)
      lwe = le0(j)
#else
      do ibio=ibio0,nbio2
      do j = 2,nz 
      lwa = lb(j)
      lwe = le(j)
#endif
      if (lwa.le.lwe) then
         llb = indwet(lwa)+1
         lle = indwet(lwe)+lazc(lwe)
         do ll = llb,lle
            k = ll-indwet(llw(ll))
            if (ljumm(llw(ll)).eq.1) then
               ldown = lazc(llw(ll))
               if (ldown.gt.1.and.k.gt.1.and.k.le.ldown-1) then
         tsa(ll,ibio) = dprod(dtdt,avd(ll,ibio))  
     1    /(dprod(dd(ll-1),dd(ll))+dd_dd(ll)) 
         tsc(ll,ibio) = dprod(dtdt,avd(ll+1,ibio)) 
     1    /(dprod(dd(ll+1),dd(ll))+dd_dd(ll)) 
         tsb(ll,ibio) = 1.d0+vtmit*tsa(ll,ibio)+vtmit*tsc(ll,ibio) 
         tsbb   = 1.d0-vtmic*tsa(ll,ibio)-vtmic*tsc(ll,ibio) 
         sd(ll,ibio) = vtmic*tsa(ll,ibio)*s(ll-1,ibio)+tsbb
     1    *s(ll,ibio)+vtmic*tsc(ll,ibio)*s(ll+1,ibio) 
               end if
            end if
         end do
      end if
      end do
      end do
c
c --------- ende vertikal diffusion 
c
#ifdef MPI
       do ibio=ibio0,nbio2
      do j = icej1,icej2
      do lw = lb0(j),le0(j)
#else
       do ibio=ibio0,nbio2
      do j = 2,nz 
      do lw = lb(j),le(j)
#endif
         nwet = indwet(lw)
         ldown = lazc(lw) 
         if (ljumm(lw).eq.1) then
               ld  = ldown-1 
               lein= min0(1,ld) 
               if (lein.eq.1) then

      tsa(nwet+1,ibio) = 0.d0 
      tsc(nwet+1,ibio) = dprod(dtdt,avd(nwet+2,ibio)) 
     1 /(dprod(dd(nwet+2),dd(nwet+1))+dd_dd(nwet+1) ) 
      tsb(nwet+1,ibio) = 1.d0+vtmit*tsc(nwet+1,ibio) 
      tsbb   = 1.d0-vtmic*tsc(nwet+1,ibio) 
      sd(nwet+1,ibio) = tsbb*s(nwet+1,ibio)+vtmic*tsc(nwet+1,ibio)
     1 *s(nwet+2,ibio) 

      tsa(nwet+ldown,ibio) = dprod(dtdt,avd(nwet+ldown,ibio)) 
     1 /(dprod(dd(nwet+ld),dd(nwet+ldown))+dd_dd(nwet+ldown) ) 
      tsc(nwet+ldown,ibio) = 0.d0 
      tsb(nwet+ldown,ibio) = 1.d0+vtmit*tsa(nwet+ldown,ibio) 
      tsbb       = 1.d0-vtmic*tsa(nwet+ldown,ibio) 
      sd(nwet+ldown,ibio)  = vtmic*tsa(nwet+ldown,ibio)*s(nwet+ld,ibio)
     1 +tsbb*s(nwet+ldown,ibio) 

               end if
         end if
      end do
      end do
            end do

#ifdef MPI
      do ibio=ibio0,nbio2
      do j = icej1,icej2 
      do lw = lb0(j),le0(j)
#else
      do ibio=ibio0,nbio2
      do j = 2,nz 
      do lw = lb(j),le(j)
#endif
         if (ljumm(lw).eq.1) then
            nwet = indwet(lw)
            ldown = lazc(lw) 
            ld  = ldown-1 
            lein= min0(1,ld) 
            if (lein.eq.1) then
                  tsal(nwet+lob,ibio) =tsb(nwet+lob,ibio) ! lob = 1 (or 2)
                  sbet(nwet+lob,ibio) = sd(nwet+lob,ibio) 
            end if
         end if
      end do
      end do
       end do
c
      do ibio=ibio0,nbio2
      do k = 1,ilo
!CDIR VECTOR NODEP
         do lh = 1,llay(0,k)
            ll = llay(lh,k)
                     tsal(ll+1,ibio) = tsb(ll+1,ibio)-(vtmit
     1                *tsa(ll+1,ibio)/tsal(ll,ibio))*vtmit
     2                *tsc(ll,ibio) 
                     sbet(ll+1,ibio) = sd(ll+1,ibio)+(vtmit
     1                *tsa(ll+1,ibio)
     1                /tsal(ll,ibio))*sbet(ll,ibio) 
         end do
      end do
               end do
c
#ifdef MPI
      do ibio=ibio0,nbio2
      do j = icej1,icej2 
      do lw = lb0(j),le0(j)
#else
      do ibio=ibio0,nbio2
      do j = 2,nz 
      do lw = lb(j),le(j)
#endif
         if (ljumm(lw).eq.1) then
            nwet = indwet(lw)
            ldown = lazc(lw) 
            ld  = ldown-1 
            lein= min0(1,ld) 
            if (lein.eq.1) then
                  sn(nwet+ldown,ibio) = sbet(nwet+ldown,ibio)
     1             /tsal(nwet+ldown,ibio) 
            end if
         end if
      end do
      end do
               end do
      
       do ibio=ibio0,nbio2
      do k = ilo,1,-1
         do lh = 1,llay(0,k)
            ll = llay(lh,k)
               sn(ll,ibio) = (sbet(ll,ibio)+vtmit
     1          *tsc(ll,ibio)*sn(ll+1,ibio))
     2          /tsal(ll,ibio) 
         end do
      end do
            end do
#ifdef MPI
      do ibio = ibio0,nbio2
      do j = icej1,icej2
         lwa = lb0(j)
         lwe = le0(j)
#else
      do ibio = ibio0,nbio2
      do j = 2,nz
         lwa = lb(j)
         lwe = le(j)
#endif
         if (lwa.le.lwe) then
            llb = indwet(lwa)+1
            lle = indwet(lwe)+lazc(lwe)
            do ll = llb,lle
               if (ljumm(llw(ll)).eq.1) then
      Tc(ll,ibio) = sn(ll,ibio)+XYadv(ll,ibio)+svert(ll,ibio)
               end if
            end do
         end if
      end do
      end do

*________ end of   find biology  _______________________

*-------- T,S ----------------
#ifdef MPI
      do J=J1,J2
      lwa=lb0(J)
      lwe=le0(J)
      if(lwe .ge. lwa) then
      nwet1=indwet(lwa)+1
      nwet2=indwet(lwe)+lazc(lwe)
      do k=nwet1,nwet2
         tec(k)=Tc(k,1)
         sac(k)=Tc(k,2)
      end do
      endif
      end do
#else
      do nwet=1,ndrei
         tec(nwet)=Tc(nwet,1)
         sac(nwet)=Tc(nwet,2)
      end do
#endif

c      write(94,'(4(2(e10.4,1x),3(f8.4,1x)))') XYadv(6754,1),
c     & svert(6754,1),sn(6754,1),
c     & Tc0(6754,1),tec(6754),XYadv(540,1),svert(540,1),sn(540,1),
c     & Tc0(540,1),tec(540),
c     & XYadv(6754,2),svert(6754,2),sn(6754,2),Tc0(6754,2),sac(6754),
c     & XYadv(540,2),svert(540,2),sn(540,2),Tc0(540,2),sac(540)

c      print*,'stromTVD finished'
      close(95)
      return 
      end 
