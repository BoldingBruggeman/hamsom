*§§1 advect
*§§2 diffus
*§§3 dythxsum
*§§4 fcuv
*§§5 felld
*§§6 felld1
*§§7 first
*§§8 form
*§§9 full
*§§10 icemod 
*§§11 icevel
*§§12 initial
*§§13 plast
*§§14 relax
*§§15 setice
*§§16 stressw
*§§17 strese  !removed not in use
*§§18 stressa
*§§19 ustern  !remove because not in use
*§§20 xmaxm
*§§21 xsum

*§§1 advect
      subroutine advect(uice1,vice1,hisi,diff1,lad)
c-----------------------------------------------------------------------
c     advect : solving continuity equation
c-----------------------------------------------------------------------
      include 'C_model.f'

      parameter (nx=m,ny=n) 

      dimension hisi(m,n,3),uice1(nx,ny,3),vice1(nx,ny,3),
     * uice(nx,ny),vice(nx,ny) 

      common/step/delice,error,deltt
      common /gitter/ dt,r,g,dl,dlr,dlrh,dln(m),rdln(m),rad,dth, 
     * dlvo(m),dlvu(m),gh,rdt,dt2,r4,coru(m),corv(m),dtrdln(m),dtdlr 
      common/array/mhc(m,n),muv(nx,ny),maphc(m,n) 
      common/fast/muuv(nx,ny)
#ifdef MPI
      include 'C_mpi.f'
      common /cord/ I1,I2,J1,J2,IA1,IA2,JA1,JA2,IB1,IB2,JB1,JB2
     1,ID1,ID2,JD1,JD2
      common /icecor/ ICEI1,ICEI2,ICEJ1,ICEJ2,ICEIA1,ICEIA2,ICEJA1,
     * ICEJA2,ICEX1,ICEX2,ICEY1,ICEY2,ICEP1,ICEP2,ICEQ1,ICEQ2,ICEX3,
     * ICEX4,ICEY3,ICEY4,ICEX5,ICEX6,ICEY5,ICEY6
      integer icheck
#endif
c     now decide if backward euler or leapfrog 
      ll = lad 
      goto (101,102,102) ll 

c     leapfrog 
 101  deltt = delice*2.0 
      k3 = 3 
      goto 103 

c     backward euler 
 102  deltt = delice 
      k3 = 2 

 103  continue 
      k2 = 2 

c     rearrange hs 
#ifdef MPI
      do j=icey5,icey6
      do i=icex5,icex6
#else
      do j=1,ny-1 
      do i=2,nx
#endif 
c        uice(i,j) = uice1(i,j,1) *muuv(nx,ny)
c        vice(i,j) = vice1(i,j,1) *muuv(nx,ny)
         uice(i,j) = uice1(i,j,1) *muuv(i,j)
         vice(i,j) = vice1(i,j,1) *muuv(i,j)
      enddo
      enddo 
#ifdef MPI
      icheck=0
      do j=JB1,JB2
      do i=IB1,IB2
#else
      do j=1,n 
      do i=1,m 
#endif
        hisi(i,j,3) = hisi(i,j,2) 
        hisi(i,j,2) = hisi(i,j,1) 
      enddo
      enddo 

 115  continue 
c     go through standard conservative advection 
c     central difference 

      delty = deltt/(8.0*dl) 
#ifndef MPI
      do j=1,ny-2 
        do i=3,nx 
        deltx = deltt/(8.0*dln(i))
        do im=1,muv(i,j) 
          hisi(i-1,j+1,1) = hisi(i-1,j+1,k3) 
     *     -deltx*(hisi(i-1,j+2,2)-hisi(i-1,j,2)) 
     *     *(uice(i-1,j+1)+uice(i,j+1)+uice(i-1,j)+uice(i,j)) 
     *     -delty*(hisi(i-2,j+1,2)-hisi(i,j+1,2)) 
     *     *(vice(i-1,j+1)+vice(i,j+1)+vice(i-1,j)+vice(i,j)) 
        enddo
        enddo
      enddo 
#else
       if(icheck .EQ. 1)then
       J_1=ICEJ1
       J_2=ICEJ2
       I_1=ICEI1
       I_2=ICEI2
       else
       J_1=ICEJA1
       J_2=ICEJA2
       I_1=ICEIA1
       I_2=ICEIA2
       endif
      do j=J_1,J_2 
      do i=I_1,I_2 
        deltx = deltt/(8.0*dln(i+1))
        do im=1,muv(i+1,j-1) 
          hisi(i,j,1) = hisi(i,j,k3) 
     *     -deltx*(hisi(i,j+1,2)-hisi(i,j-1,2)) 
     *     *(uice(i,j)+uice(i+1,j)+uice(i,j-1)+uice(i+1,j-1)) 
     *     -delty*(hisi(i-1,j,2)-hisi(i+1,j,2)) 
     *     *(vice(i,j)+vice(i+1,j)+vice(i,j-1)+vice(i+1,j-1)) 
        enddo
        enddo
      enddo 
#endif
c     now decide if done 
      goto (131,133,136) ll 
 131  goto 139 
 133  continue 

c     do backward euler correction 
#ifdef MPI
      icheck=1
      do j=JB1,JB2
      do i=IB1,IB2 
#else
      do j=1,n 
      do i=1,m
#endif
      do im=1,mhc(i,j) 
        hisi(i,j,3) = hisi(i,j,2) 
        hisi(i,j,2) = 0.5*(hisi(i,j,1)+hisi(i,j,2)) 
      enddo
      enddo
      enddo 
      ll = 3 
      k3 = 3
 
      goto 115 
 136  continue 

c     fix up h(i,j,2) 
#ifdef MPI
      icheck=1
      do j=JB1,JB2
      do i=IB1,IB2
#else
      do j=1,n 
      do i=1,m 
#endif
        hisi(i,j,2) = hisi(i,j,3) 
      enddo
      enddo 

 139  continue 
c     now do diffusion on h(i,j,k3) 
      do kd=1,2 
        goto (141,142),kd 
 141    continue 
        call diffus(uice,vice,hisi,diff1,kd) 
        goto 143 
 142    continue 
        imitt = m/2
        diff2 = -(dln(imitt)**2)/deltt 
        call diffus(uice,vice,hisi,diff2,kd) 
 143    continue 
#ifdef MPI
        do j=j1,j2
        do i=i1,i2
#else
        do j=1,n 
        do i=1,m
#endif
          hisi(i,j,1) = (hisi(i,j,1)+hisi(i,j,3))*mhc(i,j) 
        enddo
        enddo

      enddo 


      return
      end 



*§§2 diffus
c-----------------------------------------------------------------------
      subroutine diffus(uice,vice,hisi,diff1,kd) 
c-----------------------------------------------------------------------
c     diffus hisi,multiplies by delt,and puts results in hisi(3)
c-----------------------------------------------------------------------      
      include 'C_model.f'

      parameter (nx=m,ny=n) 

      dimension hisi(m,n,3),uice(nx,ny),vice(nx,ny),hisi1(m,n) 

      common /step/ delice,error,deltt
      common /gitter/ dt,r,g,dl,dlr,dlrh,dln(m),rdln(m),rad,dth, 
     * dlvo(m),dlvu(m),gh,rdt,dt2,r4,coru(m),corv(m),dtrdln(m),dtdlr
      common /array/ mhc(m,n),muv(nx,ny),maphc(m,n) 
#ifdef MPI
      common /cord/ I1,I2,J1,J2,IA1,IA2,JA1,JA2,IB1,IB2,JB1,JB2
     1,ID1,ID2,JD1,JD2
      common /icecor/ ICEI1,ICEI2,ICEJ1,ICEJ2,ICEIA1,ICEIA2,ICEJA1,
     * ICEJA2,ICEX1,ICEX2,ICEY1,ICEY2,ICEP1,ICEP2,ICEQ1,ICEQ2,ICEX3,
     * ICEX4,ICEY3,ICEY4,ICEX5,ICEX6,ICEY5,ICEY6
#endif
c     presetting
#ifdef MPI
      do j=JB1,JB2
      do i=IB1,IB2 
#else
      do j=1,n 
      do i=1,m
#endif 
        hisi1(i,j)=0.0 
      enddo
      enddo

c     do diffusion 
      deltyy=deltt*diff1/(dl**2) 
#ifdef MPI
       if(kd .EQ. 2)then
       J_1=ICEJ1
       J_2=ICEJ2
       I_1=ICEI1
       I_2=ICEI2
       else
       J_1=ICEJA1
       J_2=ICEJA2
       I_1=ICEIA1
       I_2=ICEIA2
       endif
      do j=J_1,J_2
      do i=I_1,I_2    
#else
      do j=2,ny-1 
        do i=2,nx-1 
#endif
        deltxx = deltt*diff1/(dln(i)**2)
        do im=1,mhc(i,j) 
          hisi1(i,j) = deltxx*
     *     ((hisi(i,j+1,3)-hisi(i,j,3))*mhc(i,j+1) 
     *     -(hisi(i,j,3)-hisi(i,j-1,3))*mhc(i,j-1)) 
     *     +deltyy*
     *     ((hisi(i-1,j,3)-hisi(i,j,3))*mhc(i-1,j) 
     *     -(hisi(i,j,3)-hisi(i+1,j,3))*mhc(i+1,j)) 
        enddo
        enddo
      enddo

c     result
#ifdef MPI
       if(kd .EQ. 2)then
       J_1=J1
       J_2=J2
       I_1=I1
       I_2=I2
       else
       J_1=JA1
       J_2=JA2
       I_1=IA1
       I_2=IA2
       endif
      icheck1=1
      do j=J_1,J_2
      do i=I_1,I_2
#else
      do j=1,n 
      do i=1,m
#endif 
        hisi(i,j,3) = hisi1(i,j) 
      enddo
      enddo


      return 
      end 


*§§3 dyth
c-----------------------------------------------------------------------
      subroutine dyth(uice,vice,hisi,hisir,fricei) 
c-----------------------------------------------------------------------
c     satisfy dynamic restraint
c-----------------------------------------------------------------------
      include 'C_model.f'

      parameter (nx=m,ny=n) 

      dimension uice(nx,ny,3),vice(nx,ny,3),hisi(m,n,3) 
     * ,hisir(m,n,3),fricei(m,n,3)
 
      common /step/ delice,error,deltt 
      common /gitter/ dt,r,g,dl,dlr,dlrh,dln(m),rdln(m),rad,dth, 
     * dlvo(m),dlvu(m),gh,rdt,dt2,r4,coru(m),corv(m),dtrdln(m),dtdlr
      common /array/ mhc(m,n),muv(nx,ny),maphc(m,n)  
#ifdef MPI
      common /cord/ I1,I2,J1,J2,IA1,IA2,JA1,JA2,IB1,IB2,JB1
     1,ID1,ID2,JD1,JD2
      common /icecor/ ICEI1,ICEI2,ICEJ1,ICEJ2,ICEIA1,ICEIA2,ICEJA1,
     * ICEJA2,ICEX1,ICEX2,ICEY1,ICEY2,ICEP1,ICEP2,ICEQ1,ICEQ2,ICEX3,
     * ICEX4,ICEY3,ICEY4,ICEX5,ICEX6,ICEY5,ICEY6
#endif
      hcr=0.1
#ifdef MPI
      do 2 j=ICEJ1,ICEJ2 
      do 2 i=ICEI1,ICEI2
#else
      do 2 j=2,ny-1 
      do 2 i=2,nx-1 
#endif
      do 2 im=1,mhc(i,j) 
        h = hisi(i,j,1)+hisir(i,j,1) 
        div = (uice(i,j,1)+uice(i+1,j,1)-uice(i,j-1,1)-uice(i+1,j-1,1))/
     *   (2.0*dln(i+1)) 
     *       +(vice(i,j,1)+vice(i,j-1,1)-vice(i+1,j,1)-vice(i+1,j-1,1))/
     *   (2.0*dl) 
        if(div.lt.0.0.and.fricei(i,j,1).ge.0.999) goto 902 
        fricei(i,j,1) = fricei(i,j,1)*(1.-div*delice) 
        goto 905 
 902    continue 
        if(h.lt.hcr) goto 903 
        hisir(i,j,1) = hisir(i,j,1)-h*div*delice 
        goto 905 
 903    continue 
        hisi(i,j,1)  = hisi(i,j,1)-h*div*delice 
 905    continue 
 2    continue

#ifdef MPI
      do j=J1,J2 
      do i=I1,I2 
#else
      do j=1,n 
      do i=1,m
#endif 
        hisi(i,j,1)   = max(0.0,hisi(i,j,1)) 
        hisir(i,j,1)  = max(0.0,hisir(i,j,1)) 
        fricei(i,j,1) = min(1.0,fricei(i,j,1)) 
        fricei(i,j,1) = max(0.0,fricei(i,j,1)) 
      enddo
      enddo
      return 
      end 


*§§4 fcuv
c-----------------------------------------------------------------------
      subroutine fcuv(fricei1) 
c-----------------------------------------------------------------------
c     determine the u, v and h information field
c-----------------------------------------------------------------------
      include 'C_model.f'

      parameter (nx=m,ny=n) 
      dimension fricei1(m,n) 
      common/array/mhc(m,n),muv(nx,ny),maphc(m,n) 
#ifdef MPI
      include 'C_mpi.f'
      include 'C_parallel.f'
      include 'C_index.f'
      parameter(khor1=khor+1)
      common /ind/ iwet(khor1),ldep(khor),lazc(khor),indend(n),
     *             isornr(n),isorsr(n),islab(n) 
#endif
#ifdef MPI
      do J=JA1,JA2
      do I=IA1,IA2
#else
      do j=1,n 
      do i=1,m
#endif
        mhc(i,j)=maphc(i,j) 
        if(maphc(i,j).le.0) mhc(i,j)=9 
      enddo
      enddo
#ifdef MPI
      do j=ICEJA1,ICEJA2
      do i=ICEIA1,ICEIA2
#else
      do j=2,n-1 
      do i=2,m-1
#endif 
        a0=fricei1(i-1,j)+fricei1(i+1,j)
     *    +fricei1(i,j-1)+fricei1(i,j+1) 
        if(mhc(i,j).gt.8) goto 101
        if(a0.lt.1.0e-5) mhc(i,j) = 0 
 101    continue 
      enddo
      enddo
#ifdef MPI
      do j=icey1,icey2
       do i=icex1,icex2
#else
      do j=1,ny-1 
      do i=2,nx
#endif 
        mr1 = mhc(i,j)+mhc(i,j+1)+mhc(i-1,j)+mhc(i-1,j+1) 
        if(mr1.ge.9.or.mr1.eq.0)then 
          muv(i,j) = 0 
        else 
          muv(i,j) = 1 
        endif 
      enddo
      enddo

#ifdef MPI 
      do J=JA1,JA2
      do I=IA1,IA2
#else
      do j=1,n 
      do i=1,m
#endif 
        if(mhc(i,j).eq.9) mhc(i,j) = 0 
      enddo
      enddo
#ifdef MPI
      call INTI_COMM_VAR
      call com_ivertical(0,mhc,2,2,0,0,lazc)
      call com_ivertical(0,muv,2,2,0,1,lazc)
      call com_ihorizontal(0,mhc,2,2,0,0,icdepb2,idepb2,icdept2
     *,idept2,ncsendb2,nsendb2,ncsendt2,nsendt2,nbhlayer2,nthlayer2)
      call com_ihorizontal(0,muv,2,2,0,1,icdepb2,idepb2,icdept2
     *,idept2,ncsendb2,nsendb2,ncsendt2,nsendt2,nbhlayer2,nthlayer2)
      call MPI_Waitall(nreq1, ir_r_b, status, ierr)
      INUM1=0
      INUM2=0
      call VER_IUNPACK(0,mhc,2,2,0,lazc)
      call VER_IUNPACK(0,muv,2,2,0,lazc)
      INUM3=0
      INUM4=0
      call HOR_IUNPACK(0,mhc,2,2,nbhlayer2,nthlayer2,
     &idepb2,idept2,icdepb2,icdept2,0,iwet)
      call HOR_IUNPACK(0,muv,2,2,nbhlayer2,nthlayer2,
     &idepb2,idept2,icdepb2,icdept2,0,iwet)
#endif
      return 
      end 



*§§5 felld
c-----------------------------------------------------------------------
      subroutine felld(uice,vice,eta,f,i,j,k,feta) 
c-----------------------------------------------------------------------
      include 'C_model.f'

      parameter (nx=m,ny=n) 

      dimension uice(nx,ny,3),vice(nx,ny,3),eta(m,n),f(4),feta(4,nx,ny)
      common /dinv/ s1(m),s2,s3(m) 

      f(1) = feta(1,i,j)+s1(i)*(+uice(i,j-1,k)*(eta(i-1,j)+eta(i,j))) 
      f(2) = feta(2,i,j)+s2   *(+uice(i+1,j,k)*(eta(i,j)  +eta(i,j+1))) 
      f(3) = feta(3,i,j)+0.5*s3(i)*(vice(i+1,j-1,k)*eta(i,j)+ 
     * vice(i+1,j,k)*(-eta(i,j) 
     * +eta(i,j+1))-vice(i+1,j+1,k)*eta(i,j+1)+vice(i,j-1,k)*(-eta(i-1,j
     * )+eta(i,j)) 
     * -vice(i-1,j-1,k)*eta(i-1,j)) 
      f(4) = feta(4,i,j)+s3(i)*0.5* 
     * (vice(i+1,j-1,k)*eta(i,j)+vice(i+1,j,k)*(-eta(i,j+1) 
     * +eta(i,j))-vice(i+1,j+1,k)*eta(i,j+1)+vice(i,j-1,k)*(eta(i-1,j) 
     * -eta(i,j)) 
     * -vice(i-1,j-1,k)*eta(i-1,j)) 


      return 
      end 


*§§6 felld1
      subroutine felld1(uice,vice,eta,f,k) 
c-----------------------------------------------------------------------
c     felld1,felld : calculating ice stress i
c-----------------------------------------------------------------------
      include 'C_model.f'

      parameter (nx=m,ny=n) 

      dimension uice(nx,ny,3),vice(nx,ny,3),eta(m,n),f(4,nx,ny) 
      common/dinv/s1(m),s2,s3(m) 
#ifdef MPI
      include 'C_mpi.f'
      common /icecor/ ICEI1,ICEI2,ICEJ1,ICEJ2,ICEIA1,ICEIA2,ICEJA1,
     * ICEJA2,ICEX1,ICEX2,ICEY1,ICEY2,ICEP1,ICEP2,ICEQ1,ICEQ2,ICEX3,
     * ICEX4,ICEY3,ICEY4,ICEX5,ICEX6,ICEY5,ICEY6
#endif
#ifdef MPI
      do j=ICEQ1,ICEQ2 
      do i=ICEP1,ICEP2
#else
      do j=2,ny-2 
      do i=3,nx-1
#endif 
        f(1,i,j) = s1(i)*(uice(i,j+1,k)*(eta(i-1,j+1)+eta(i,j+1))) 
        f(2,i,j) = s2   *(uice(i-1,j,k)*(eta(i-1,j+1)+eta(i-1,j))) 
        f(3,i,j) = s3(i)*
     *   (vice(i,j,k)*  (-eta(i,j)-eta(i-1,j+1)+eta(i,j+1)+eta(i-1,j))
     *   +vice(i,j+1,k)*(-eta(i,j+1)+eta(i-1,j+1)) 
     *   +vice(i-1,j,k)*(-eta(i-1,j+1)+eta(i-1,j)) 
     *   +vice(i-1,j+1,k)*eta(i-1,j+1)) 
        f(4,i,j) = s3(i)*
     *   (vice(i,j,k)*  (eta(i-1,j)+eta(i,j+1)-eta(i,j)-eta(i-1,j+1))
     *   +vice(i,j+1,k)*(eta(i,j+1)-eta(i-1,j+1)) 
     *   +vice(i-1,j,k)*(eta(i-1,j+1)-eta(i-1,j)) 
     *   +vice(i-1,j+1,k)*eta(i-1,j+1)) 
        f(3,i,j) = 0.5*f(3,i,j) 
        f(4,i,j) = 0.5*f(4,i,j) 
      enddo
      enddo
      return 
      end 


*§§7 first   
      subroutine first(iizet,mm,nn,iffir,isice) 
c-----------------------------------------------------------------------
      include 'C_model.f'
#ifdef MPI
      include 'C_mpi.f'
      include 'C_parallel.f'
      include 'C_index.f'
#endif
      parameter (nx=m,ny=n,n4=nx*ny,n3=n*m) 
      dimension iizet(mm,nn)
      character*6 date,hour*10,aaaa 
      common /step/ delice,error,deltt
      common /dinv/ s1(m),s2,s3(m) 
      common /gitter/ dt,r,g,dl,dlr,dlrh,dln(m),rdln(m),rad,dth, 
     * dlvo(m),dlvu(m),gh,rdt,dt2,r4,coru(m),corv(m),dtrdln(m),dtdlr
      common /array/ mhc(m,n),muv(nx,ny),maphc(m,n)  
      common /stresi/ stresi(m,n,3) 
      common /more/ diff1,lad,icount,nst6
      common /ice/ cwa,ccw,tci,sice,rois,cp,tf,epsis,roil,dtroil,
     *  his(m,n),frice(m,n),tis(m,n),tair(m,n),qoc(m,n),icloud,
     *  solcon,frsw,frsi,time,xlat(m),ui(m,n),vi(m,n),us(m,n),vs(m,n)
      common /icevelo/ uice(nx,ny,3),vice(nx,ny,3) 
     *  ,uicec(nx,ny),vicec(nx,ny) 
     *  ,uerr(nx,ny),verr(nx,ny)
     *  ,gairx(nx,ny),gairy(nx,ny),gwatx(nx,ny),gwaty(nx,ny) 
     *  ,drags(nx,ny),draga(nx,ny),amass(nx,ny) 
     *  ,eta(m,n),zeta(m,n),hisir(m,n,3) 
     *  ,hisi(m,n,3),fricei(m,n,3) 
     *  ,forcex(nx,ny),forcey(nx,ny)     
      common /hh/ hisr(m,n)
      dimension u1(nx,ny,3),v1(nx,ny,3) 
c      dimension gairxx(nx,ny),gairyy(nx,ny)

c     u1,v1 : used for output 
c      equivalence(uice,u1) 
c     equivalence(vice,v1) 
      dimension fricei1(m,n),hisi1(m,n) 
c     fricei1,hisi1 : used for output 
      equivalence(fricei,fricei1) 


c     time  steps 
      deltt  = 4.*60.*60. 
      delice = deltt 
      error  = 0.001 

c     diff1 : diffuse coefficient 
      diff1  = 0.004*dl 
      lad    = 2 
      icount = 0 
  
c     read initial fields and find ice-covered extent 
      if(iffir.eq.0) then
         call initial(fricei,hisi,hisir,iizet,m,n)
      end if

      call fcuvs(fricei1) 

      do j=1,n 
      do i=1,m 
        hisi1(i,j)=hisi(i,j,1)+hisir(i,j,1) 
      enddo
      enddo

      call xsum(hisi1,fricei1,summass,sumfricei,sumhisi) 

c      write(6,1111) icount,summass,sumfricei,sumhisi
 1111 format(1x,'time step=',i3,'  mass sum=',e10.4,' fricei sum=',
     *       e12.4,' thick sum=',e12.4)
 
c     get first value of u and v 
cc------------------------------------
cc  ute 01082017 test intiatlize zeta & eta eq 0
cc------------------------------------      
      do j=1,ny-1 
      do i=2,nx 
        zeta(i,j)= 0. 
        eta(i,j)= 0. 
      enddo
      enddo
cc--------------------------------------------
c      print*,'  icemodell 2: ',maxval(eta),maxval(zeta)
c      print*,'  icemodell 3: ',maxval(gwatx),maxval(gwaty)
      if(isice.eq.0)then
      do j=1,ny-1 
      do i=2,nx 
        uice(i,j,1)= 0.02*gairx(i,j)*muv(i,j) 
        uice(i,j,2)=uice(i,j,1) 
        uice(i,j,3)=uice(i,j,1) 
        vice(i,j,1)=0.02*gairy(i,j)*muv(i,j) 
        vice(i,j,2)=vice(i,j,1) 
        vice(i,j,3)=vice(i,j,1) 
        uicec(i,j)=uice(i,j,1) 
        vicec(i,j)=vice(i,j,1) 
      enddo
      enddo
      theta=1.0 
      call form(uice,vice,eta,zeta,drags,draga,gairx,gairy,gwatx 
     *   ,gwaty,forcex,forcey,hisi1,amass,fricei1) 

c     set amass=0 and define initial viscosity 
      do j=1,ny-1 
      do i=2,nx 
        amass(i,j)=0.0 
        zeta(i,j)=hisi1(i,j)*(1.0e+11) 
        eta(i,j)=zeta(i,j)/4.0 
      enddo
      enddo
      call relax(uice,vice,eta,zeta,drags,draga,amass,forcex,forcey 
     *   ,theta,uicec,vicec) 
      elseif(isice.eq.1)then
      do j=1,ny-1 
      do i=2,nx 
        uice(i,j,1)=ui(i,j) 
        uice(i,j,2)=uice(i,j,1) 
        uice(i,j,3)=uice(i,j,1) 
        vice(i,j,1)=vi(i,j) 
        vice(i,j,2)=vice(i,j,1) 
        vice(i,j,3)=vice(i,j,1) 
        uicec(i,j)=uice(i,j,1) 
        vicec(i,j)=vice(i,j,1) 
      enddo
      enddo
      do j=1,ny-1 
      do i=2,nx 
        amass(i,j)=0.0 
        zeta(i,j)=hisi1(i,j)*(1.0e+11) 
        eta(i,j)=zeta(i,j)/4.0 
      enddo
      enddo
      endif
      return
      end



*§§8 form
      subroutine form(uice,vice,eta,zeta,drags,draga,gairx,gairy, 
     * gwatx,gwaty,forcex,forcey,hisi1,amass,fricei1) 
c-----------------------------------------------------------------------
c     form : calculate parameters for relaxation
c-----------------------------------------------------------------------
      include 'C_model.f'

      parameter (nx=m,ny=n) 
      dimension uice(nx,ny,3),vice(nx,ny,3) 
     &,eta(m,n),zeta(m,n) 
     &,drags(nx,ny),draga(nx,ny),dairn(nx,ny),dwatn(nx,ny) 
     &,gairx(nx,ny),gairy(nx,ny),gwatx(nx,ny),gwaty(nx,ny) 
     &,cor(nx,ny),forcex(nx,ny),forcey(nx,ny),amass(nx,ny) 
     &,hisi1(m,n),fricei1(m,n),pressi(m,n) 
      common/step/delice,error,deltt
      common /gitter/ dt,r,g,dl,dlr,dlrh,dln(m),rdln(m),rad,dth, 
     ddlvo(m),dlvu(m),gh,rdt,dt2,r4,coru(m),corv(m),dtrdln(m),dtdlr 
      common/array/mhc(m,n),muv(nx,ny),maphc(m,n) 
       common/fast/muuv(nx,ny)
#ifdef MPI
      include 'C_mpi.f'
      common /cord/ I1,I2,J1,J2,IA1,IA2,JA1,JA2,IB1,IB2,JB1,JB2
     1,ID1,ID2,JD1,JD2
      common /icecor/ ICEI1,ICEI2,ICEJ1,ICEJ2,ICEIA1,ICEIA2,ICEJA1,
     * ICEJA2,ICEX1,ICEX2,ICEY1,ICEY2,ICEP1,ICEP2,ICEQ1,ICEQ2,ICEX3,
     * ICEX4,ICEY3,ICEY4,ICEX5,ICEX6,ICEY5,ICEY6
#endif
c----- set up basic constants 
      pai1=3.14159265/180.0 
      dwat=0.0055 
      dair=0.0018 
      rhowat=1010.0 
      rhoair=1.29 
      rhoice=910.0 
      fcor=1.46e-4*sin(60.0*pai1) 
      sinwin=0.0000 
      coswin=1.0000 
c  Achtung von 16. degree auf 0. wg. Umsetzung auf eine reale nicht
c  geostrophische Wassergeschwindigkeit
      sinwat=sin(16.*pai1) 
      sinwat=0.
      coswat=cos(16.0*pai1)
      coswat=1.
      
cc      g=9.81 
*20 deg gives sin equal to .34202 
*now set up mass per unit fricei and coriolis term 
#ifdef MPI
      do 101 j=icey1,icey2
      do 101 i=icex1,icex2
#else
      do 101 j=1,ny-1 
      do 101 i=2,nx
#endif 
      cor(i,j)=0. 
      do 101 im=1,muv(i,j) 
      amass(i,j)=rhoice*0.25* 
     & (hisi1(i,j)+hisi1(i,j+1)+hisi1(i-1,j)+hisi1(i-1,j+1)) 
      cor(i,j)=amass(i,j)*fcor 
 101  continue 
c----- set up non linear,wind and water drag 
#ifdef MPI
      do 99 j=icey1,icey2
      do 99 i=icex1,icex2
#else
      do 99 j=1,ny-1 
      do 99 i=2,nx
#endif
      dwatn(i,j)=0. 
      do 99 im=1,muv(i,j) 
    
      dwatn(i,j)=rhowat*dwat*sqrt((uice(i,j,1)-gwatx(i,j))**2 
     & +(vice(i,j,1)-gwaty(i,j))**2) 
c----- set up symmettric drag 
      drags(i,j)=dwatn(i,j)*coswat 
c----- set up anti symmettric drag plus coriolis 
      draga(i,j)=dwatn(i,j)*sinwat+cor(i,j) 
 99   continue 
c----- set up forcing field 
#ifdef MPI
      do 104 j=icey1,icey2
      do 104 i=icex1,icex2
#else
      do 104 j=1,ny-1 
      do 104 i=2,nx
#endif 
      do 104 im=1,muv(i,j) 
    
c----- add in current and wind force 
      forcex(i,j)=gairx(i,j)+dwatn(i,j)* 
     & (coswat*gwatx(i,j)-sinwat*gwaty(i,j)) 
      forcey(i,j)=gairy(i,j)+dwatn(i,j)* 
     & (sinwat*gwatx(i,j)+coswat*gwaty(i,j)) 
 104  continue 
c----- set up ice pressiure strangth coefficient 
      pstar=1.0e4 
c----- set up values 
#ifdef MPI
       if(myid.eq.0)print*,'press',maxval(pressi),minval(pressi)
      do 105 j=ja1,ja2
      do 105 i=ia1,ia2
#else
      do 105 j=1,n 
      do 105 i=1,m
#endif 
c15102017 Ute die Belegung von pressi war nicht korrekt
        pressi(i,j)=0. 
      do 105 im=1,mhc(i,j) 
       pressi(i,j)=pstar*hisi1(i,j)*exp(-20.0*(1.0-fricei1(i,j))) 
 105  continue 
       if(myid.eq.0)print*,'press',maxval(pressi),minval(pressi)
c----- calculate stresi, strain, zeta and eta %ca6 not here
      call plast(uice,vice,pressi,eta,zeta) 
           
c----- calculate pressiure force and add to external force 
#ifdef MPI
      do 107 j=icey1,icey2
      do 107 i=icex1,icex2
#else
      do 107 j=1,ny-1 
      do 107 i=2,nx
#endif 
      do 107 im=1,muv(i,j) 
      forcex(i,j)=forcex(i,j)-(0.25/dln(i))* 
     &  (pressi(i,j+1)+pressi(i-1,j+1)-pressi(i,j)-pressi(i-1,j)) 
      forcey(i,j)=forcey(i,j)-(0.25/dl)* 
     &  (pressi(i-1,j)+pressi(i-1,j+1)-pressi(i,j)-pressi(i,j+1)) 
 107  continue 
      return 
      end 

*§§9 full
      SUBROUTINE full (m,n,S)
c-----------------------------------------------------------------------
c    fill 2-D arrays (land points)
c-----------------------------------------------------------------------

      dimension S(M,N)

      DIMENSION SA(M,N)

      do i=1,m
      do j=1,n
        sa(i,j) = -990.
      enddo
      enddo

      DO IOFT=1,1000

        do i=1,m
        do j=1,n
          sa(i,j) = s(i,j)
        enddo
        enddo
      
        DO  J=2,N-1
        DO  I=2,M-1
          IF(Sa(I,J).gE.-900.) THEN
            IF(S(I,J-1).lT.-900.) S(I,J-1) = Sa(I,J)
            IF(S(I-1,J).lT.-900.) S(I-1,J) = Sa(I,J)
            IF(S(I+1,J).lT.-900.) S(I+1,J) = Sa(I,J)
            IF(S(I,J+1).lT.-900.) S(I,J+1) = Sa(I,J)
c            IF(S(I-1,J+1).lT.-900.) S(I-1,J+1) = Sa(I,J)
c            IF(S(I+1,J+1).lT.-900.) S(I+1,J+1) = Sa(I,J)
c            IF(S(I-1,J-1).lT.-900.) S(I-1,J-1) = Sa(I,J)
c            IF(S(I+1,J-1).lT.-900.) S(I+1,J-1) = Sa(I,J)
          ENDIF
        enddo
        enddo

        JOJO = 0
        DO J=2,N-1
        DO I=2,M-1
          if(s(i,j).lt.-900.) jojo = 1
        enddo
        enddo

        IF(JOJO.EQ.0) GOTO 2
      enddo

      PRINT*,'after 1000 iterations not fully filled'
  2   CONTINUE

      do i=1,m
        s(i,1) = s(i,2)
        s(i,n) = s(i,n-1)
      enddo

      do j=1,n
        s(1,j) = s(2,j)
        s(m,j) = s(m-1,j)
      enddo
      
      s(1,1) = s(2,2)
      s(m,n) = s(m-1,n-1)
      s(1,n) = s(2,n-1)
      s(m,1) = s(m-1,2)

      RETURN
      END





*§§10 icemod 
      subroutine icemod 
c-----------------------------------------------------------------------
c     baltic sea ice model by zhanhai zhang; after lepparanta
c     this is a program for viscous plastic sea ice model 
c     h : 38*48    u,v : 37*47 
c-----------------------------------------------------------------------
c     hisi : the level ice thickness 
c     hisir : the ridged ice thickness 
c     fricei : the ice compactness 
c     uice : the ice velocity in x_direction 
c     vice : the ice velocity in y_direction 
c     (nx,ny) : the number of grid points of velocity 
c     (n,m) : the number of grid points of thickness etc. 
c     uicec,vicec : ice velocity used for calculating advective term 
c     eta : shear viscosity coefficient 
c     zeta : bulk viscosity coefficient 
c     gairx,gairy : geostraphical wind speed 
c     gwatx,gwaty : geostraphical current 
c     amass : the ice mass 
c     drags,graga : the water & air coefficient 
c     delice : the time step 
c     dln,dl : the grid distance in x_ & y_ direction 
c     mhc,muv : the thickness & velocity information field 
c-----------------------------------------------------------------------
      include 'C_model.f'

      parameter (nx=m,ny=n,n4=nx*ny,n3=n*m)  

      common/icevelo/ uice(nx,ny,3),vice(nx,ny,3) 
     &,uicec(nx,ny),vicec(nx,ny)
     &,uerr(nx,ny),verr(nx,ny)
     &,gairx(nx,ny),gairy(nx,ny),gwatx(nx,ny),gwaty(nx,ny) 
     &,drags(nx,ny),draga(nx,ny),amass(nx,ny) 
     &,eta(m,n),zeta(m,n),hisir(m,n,3) 
     &,hisi(m,n,3),fricei(m,n,3)
     &,forcex(nx,ny),forcey(nx,ny) 
      character*6 date,hour*10,aaaa 
      common/step/delice,error,deltt
      common/dinv/s1(m),s2,s3(m) 
      common /gitter/ dt,r,g,dl,dlr,dlrh,dln(m),rdln(m),rad,dth, 
     ddlvo(m),dlvu(m),gh,rdt,dt2,r4,coru(m),corv(m),dtrdln(m),dtdlr 
      common/more/ diff1,lad,icount,nst6
      common/array/mhc(m,n),muv(nx,ny),maphc(m,n) 
       common/fast/muuv(nx,ny)
      common/stresi/stresi(m,n,3) 
      dimension u1(nx,ny),v1(nx,ny) 
c----- u1,v1 : used for output 
      equivalence(uice,u1) 
      equivalence(vice,v1) 
c----- fricei1,hisi1 : used for output 
       dimension fricei1(m,n),hisi1(m,n) 
      equivalence(fricei,fricei1) 
c     equivalence(hisi,hisi1) 
c----- uerr,verr : the velocity different 
       dimension uerrv(nx,ny),verrv(nx,ny) 
      equivalence (uerr,uerrv),(verr,verrv)  
#ifdef MPI
      include 'C_mpi.f'
      include 'C_parallel.f'
      include 'C_index.f'
      parameter(khor1=khor+1)
      common /ind/ iwet(khor1),ldep(khor),lazc(khor),indend(n),
     *             isornr(n),isorsr(n),islab(n)
#endif

#ifdef MPI
      do 121 j=icey5,icey6
      do 121 i=icex5,icex6
#else
      do 121 j=1,ny-1 
      do 121 i=2,nx
#endif 
      uice(i,j,3)=uice(i,j,1) 
      vice(i,j,3)=vice(i,j,1) 
      uicec(i,j)=uice(i,j,1) 
      vicec(i,j)=vice(i,j,1) 
 121  continue

#ifdef MPI
      do 122 j=JB1,JB2 
      do 122 i=IB1,IB2
#else 
      do 122 j=1,n 
      do 122 i=1,m
#endif 
      hisi1(i,j)=hisi(i,j,1)+hisir(i,j,1) 
c      fricei1(i,j)=fricei(i,j,1) 
 122  continue 
      theta=1.0 
      delice=deltt/2.0 
      
       call fcuv(fricei1) 
       call form(uice,vice,eta,zeta,drags,draga,gairx,gairy,gwatx 
     &,gwaty,forcex,forcey,hisi1,amass,fricei1) 
       
      call relax(uice,vice,eta,zeta,drags,draga,amass,forcex,forcey 
     &,theta,uicec,vicec) 

c----- do regular time step 
c----- doing backwards time step 
      theta=1.0 
      delice=deltt 
!----ca9 hier passiert was ------
      call form(uice,vice,eta,zeta,drags,draga,gairx,gairy,gwatx
     &,gwaty,forcex,forcey,hisi1,amass,fricei1)
c----- set u(1)=u(2) and same for v 
#ifdef MPI
      do 111 j=icey5,icey6
      do 111 i=icex5,icex6
#else
      do 111 j=1,ny-1 
      do 111 i=2,nx
#endif 
      uice(i,j,3)=uice(i,j,1) 
      vice(i,j,3)=vice(i,j,1) 
      uicec(i,j)=uice(i,j,1) 
      vicec(i,j)=vice(i,j,1)
      uice(i,j,1)=uice(i,j,2)
      vice(i,j,1)=vice(i,j,2) 
 111  continue 
       call relax(uice,vice,eta,zeta,drags,draga,amass,forcex,forcey 
     &,theta,uicec,vicec)
c--- solving contuity equation 
       call advect(uice,vice,hisi,diff1,lad) 
       call advect(uice,vice,fricei,diff1,lad)  

       call advect(uice,vice,hisir,diff1,lad) 
       call dyth(uice,vice,hisi,hisir,fricei) 

c-----find ice-covered extent and calculate ice mass 
c      call fcuv(fricei1) 
c-----print out every 24 hours 
       icount=16
c       print*,'ice icount',icount
       if((icount-icount/4*4).ne.0)goto 96 
 1111 format(3x,'mass sum=',e10.4,' fricei sum=', 
     &       e12.4,' thick sum=',e12.4) 
c      call outha(fricei1,icount,'compa',100.,date) 
#ifdef MPI
      do 502 j=j1,j2
      do 502 i=i1,i2
#else
      do 502 j=1,n 
      do 502 i=1,m 
#endif
      hisi1(i,j)=hisi(i,j,1)+hisir(i,j,1) 
c      if(muv(i,j) .gt. 0 .or. mhc(i,j) .gt. 0) then
c      write(*,*)muv(i,j),mhc(i,j),myid,izet(i,j)
c      endif
 502  continue 
      call xsum(hisi1,fricei1,summass,sumfricei,sumhisi) 
c----- print out square velocity and velocity diference and max change 
      sq1=0.0 
      sq=0.0 
#ifdef MPI
      do 130 j=icey1,icey2
      do 130 i=icex1,icex2 
#else
      do 130 j=1,ny-1 
      do 130 i=2,nx
#endif 
      sq=sq+uice(i,j,1)**2+vice(i,j,1)**2 
      uerr(i,j)=uice(i,j,1)-uice(i,j,2) 
      verr(i,j)=vice(i,j,1)-vice(i,j,2) 
      sq1=sq1+uerr(i,j)**2+verr(i,j)**2 
 130  continue 
#ifdef MPI
      pq=0.0
      pq1=0.0
      call MPI_ALLREDUCE(sq,pq,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,IERR)
      CALL MPI_ALLREDUCE(sq1,pq1,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,IERR)
      sq1=pq1
      sq=pq
#endif
      call xmaxm(uerrv,smu) 
      call xmaxm(verrv,smv) 
      
      sm=max(smu,smv)
#ifdef MPI 
      CALL MPI_REDUCE(sm,sm1,1,MPI_REAL,MPI_MAX,0,MPI_COMM_WORLD,IERR)
      if(myid .eq. 0) sm=sm1
#endif
      IF(myid .eq. 0) write(7,5)sq,sq1,sm 
c      write(6,5)sq,sq1,sm 
 5    format(3x,'square vel, sq. vel diff, max change' 
     &/1x,3e18.10) 
 96   continue 
      end 


*§§11 icevel
      subroutine icevel
c-----------------------------------------------------------------------
c     compute ice velocities
c-----------------------------------------------------------------------
      include 'C_model.f'

      parameter(nx=m,ny=n)
      parameter(khor1=khor+1)    
      include 'C_ice.f'
      common/fast/muuv(nx,ny)
#ifdef MPI
      include 'C_mpi.f'
      common /cord/ I1,I2,J1,J2,IA1,IA2,JA1,JA2,IB1,IB2,JB1,JB2
     1,ID1,ID2,JD1,JD2
c      common /icecor/ ICEI1,ICEI2,ICEJ1,ICEJ2,ICEIA1,ICEIA2,ICEJA1,
c     * ICEJA2,ICEX1,ICEX2,ICEY1,ICEY2,ICEP1,ICEP2,ICEQ1,ICEQ2,ICEX3,
c     * ICEX4,ICEY3,ICEY4,ICEX5,ICEX6,ICEY5,ICEY6
#endif
#ifdef MPI
      do k=j1,j2
      do i=i1,i2
#else
      do k=1,n
      do i=1,m
#endif
c        ui(i,k) = 0.
c        vi(i,k) = 0.
        im = i+1
        jm =k-1         
        im = max(im,3)
        jm = max(jm,1)
        im = min(nx,im)
        jm = min(ny-2,jm)
        vi(i,k) = (vice(im,jm,1)*real(muuv(im,jm))+
     *   vice(im,jm+1,1)*real(muuv(im,jm+1)))/2.
        ui(i,k) = (uice(im,jm+1,1)*real(muuv(im,jm+1)) +
     *   uice(im-1,jm+1,1)*real(muuv(im-1,jm+1)))/2. 
      enddo
      enddo
      return 
      end


*§§12 initial
      subroutine initial(fricei,hisi,hisir,iizet,mm,nn) 
c-----------------------------------------------------------------------
c     initial : read in initial field of hisi & fricei
c-----------------------------------------------------------------------
      include 'C_model.f'

      parameter (nx=m,ny=n) 

      dimension fricei(mm,nn,3),hisi(mm,nn,3),hisir(mm,nn,3) 
      dimension ia(m,n),ih(m,n) 
      dimension iizet(mm,nn)
      character*6 date,dddd 

      common /array/ mhc(m,n),muv(nx,ny),maphc(m,n) 
      common /fast/ muuv(nx,ny)
      common /dinv/ s1(m),s2,s3(m) 
      common /step/ delice,error,deltt 

      common /gitter/ dt,r,g,dl,dlr,dlrh,dln(m),rdln(m),rad,dth, 
     *  dlvo(m),dlvu(m),gh,rdt,dt2,r4,
     *  coru(m),corv(m),dtrdln(m),dtdlr
   
      common /ice/ cwa,ccw,tci,sice,rois,cp
      common /ice/ tf,epsis,roil,dtroil
      common /ice/ his(m,n),frice(m,n),tis(m,n)
      common /ice/ tair(m,n),qoc(m,n)
      common /ice/ icloud,solcon,frsw,frsi,time,xlat(m)
      common /ice/ ui(m,n),vi(m,n),us(m,n),vs(m,n) 
      common /hh/ hisr(m,n)
  
      epsis=1.e-7

      do j=1,n 
      do i=1,m
        maphc(i,j) = 1
      enddo
      enddo

      do j=1,n
      do i=1,m
        if(iizet(i,j).le.0) maphc(i,j) = 0
      enddo
      enddo
 
      do j=1,n 
      do i=1,m 
c        jm = m-j+1
        hisi(i,j,1) = his(i,j)*maphc(i,j) 
        hisi(i,j,2) = 0.0 
        hisi(i,j,3) = 0.0 
        fricei(i,j,1) = frice(i,j)*maphc(i,j) 
        fricei(i,j,2) = fricei(i,j,1) 
        fricei(i,j,3) = fricei(i,j,1) 
        hisir(i,j,1) = hisr(i,j)*maphc(i,j) 
        hisir(i,j,2) = 0.0 
        hisir(i,j,3) = 0.0 
      enddo
      enddo
c       if(myid.eq.0)print*,'inside initial muuv'
c     fast ice
      do j=1,ny-1
      do i=2,nx 
        muuv(i,j) = maphc(i-1,j+1)
        if(iizet(i-1,j+1).le.10) muuv(i,j) = 0
      enddo
      enddo
      
c     set up delta-x arrays s1, s3 and dln, and s2 
      do j=1,m 
c        jm = m-j+1
        s1(j) = 0.5000/(dln(j)*dln(j)) 
        s3(j) = 0.5000/(dl*dln(j)) 
      enddo
      s2 = 0.5000/(dl*dl) 
      return 
      end 



*§§13 plast
      subroutine plast(uice,vice,pressi,eta,zeta) 
c-----------------------------------------------------------------------
c     plast : calculates strain rates and viscous parameters
c-----------------------------------------------------------------------
      include 'C_model.f'

      parameter (nx=m,ny=n) 

      dimension uice(nx,ny,3),vice(nx,ny,3),pressi(m,n) 
     * ,eta(m,n),zeta(m,n),e11(m,n),e22(m,n),e12(m,n)
 
      common/step/delice,error,deltt 
      common /gitter/ dt,r,g,dl,dlr,dlrh,dln(m),rdln(m),rad,dth, 
     * dlvo(m),dlvu(m),gh,rdt,dt2,r4,coru(m),corv(m),dtrdln(m),dtdlr 
      common/stresi/stresi(m,n,3) 
      common/array/mhc(m,n),muv(nx,ny),maphc(m,n) 
#ifdef MPI
      include'C_mpi.f'
      common /cord/ I1,I2,J1,J2,IA1,IA2,JA1,JA2,IB1,IB2,JB1,JB2
     1,ID1,ID2,JD1,JD2
      common /icecor/ ICEI1,ICEI2,ICEJ1,ICEJ2,ICEIA1,ICEIA2,ICEJA1,
     * ICEJA2,ICEX1,ICEX2,ICEY1,ICEY2,ICEP1,ICEP2,ICEQ1,ICEQ2,ICEX3,
     * ICEX4,ICEY3,ICEY4,ICEX5,ICEX6,ICEY5,ICEY6
      common /lbkhor/ khorl,lzet(khor),lb0(n),le0(n),
     * lb1(n),le1(n),lb2(n),le2(n)
#endif
      eccen = 2.0 
      ecm2  = 1.0/(eccen**2) 
      gmin  = 2.0e-7 
      zmin  = 1.0e8 
#ifdef MPI
      do  j=JA1,JA2
      do  i=IA1,IA2
#else
      do j=1,n
      do i=1,m
#endif
      zeta(i,j) = zmin
      eta(i,j)  = ecm2*zeta(i,j)
      enddo
      enddo

c---- evaluate strain rates 
#ifdef MPI
      do 101 j=ICEJA1,ICEJA2
      do 101 i=ICEIA1,ICEIA2
#else
      do 101 j=2,ny-1 
      do 101 i=2,nx-1
#endif 
      do 101 im=1,mhc(i,j) 
      e11(i,j)=(0.5/dln(i))*(uice(i,j,1)+uice(i+1,j,1) 
     & -uice(i,j-1,1)-uice(i+1,j-1,1)) 
      e22(i,j)=(0.5/dl)*(vice(i,j,1)+vice(i,j-1,1) 
     & -vice(i+1,j,1)-vice(i+1,j-1,1)) 
      e12(i,j)=(0.25/dl)*(uice(i,j,1)+uice(i,j-1,1) 
     & -uice(i+1,j,1)-uice(i+1,j-1,1)) 
     & +(0.25/dln(i))*(vice(i,j,1)+vice(i+1,j,1) 
     & -vice(i,j-1,1)-vice(i+1,j-1,1)) 
c---- evaluate viscosities 
      delt=(e11(i,j)**2+e22(i,j)**2)*(1.0+ecm2)+4.0*ecm2*e12(i,j)**2 
     & +2.0*e11(i,j)*e22(i,j)*(1.0-ecm2) 
      delt1=sqrt(delt) 
c---- put min and max viscosities in 
      delt1=max(gmin,delt1) 
      zeta(i,j)=0.5*pressi(i,j)/delt1 
      zeta(i,j)=max(zmin,zeta(i,j)) 
      eta(i,j)=ecm2*zeta(i,j) 
 101  continue 

c---- calculate ice stresi 
#ifdef MPI
      do 102 j=JA1,JA2
      do 102 i=IA1,IA2
#else
      do 102 j=1,n 
      do 102 i=1,m
#endif 
      do 102 im=1,mhc(i,j) 
        ss11=(zeta(i,j)-eta(i,j))*(e11(i,j)+e22(i,j))-pressi(i,j)*0.5 
        stresi(i,j,1) = 2.0*eta(i,j)*e11(i,j)+ss11 
        stresi(i,j,2) = 2.0*eta(i,j)*e22(i,j)+ss11 
        stresi(i,j,3) = 2.0*eta(i,j)*e12(i,j) 
 102  continue 
      return
      end 



*§§14 relax
      subroutine relax(uice,vice,eta,zeta,drags,draga,amass,
     *  forcex,forcey,theta,uicec,vicec)
c-----------------------------------------------------------------------
c     relax : iteration to obtain ice velocity
c-----------------------------------------------------------------------
      include 'C_model.f'

      parameter (nx=m,ny=n,n4=nx*ny) 
      dimension uice(nx,ny,3),vice(nx,ny,3),eta(m,n),zeta(m,n) 
     &,drags(nx,ny),draga(nx,ny),forcex(nx,ny),forcey(nx,ny),fxeta(4) 
     &,fxzeta(4),fyeta(4),fyzeta(4),uerr(nx,ny),uerrv(nx,ny) 
     &,verr(nx,ny),verrv(nx,ny),amass(nx,ny) 
     &,coefix(nx,ny),coefiy(nx,ny),flx(nx,ny),fly(nx,ny)
     &,uicec(nx,ny),vicec(nx,ny),fls(nx,ny) 
     &,fxe(4,nx,ny),fye(4,nx,ny),fxz(4,nx,ny),fyz(4,nx,ny) 
      common/dinv/s1(m),s2,s3(m) 
      common/step/delice,error,deltt
         common /gitter/ dt,r,g,dl,dlr,dlrh,dln(m),rdln(m),rad,dth, 
     ddlvo(m),dlvu(m),gh,rdt,dt2,r4,coru(m),corv(m),dtrdln(m),dtdlr 
      common/array/mhc(m,n),muv(nx,ny),maphc(m,n) 
c      equivalence (uerr,uerrv),(verr,verrv) 
#ifdef MPI
      include 'C_index.f'
      include 'C_mpi.f'
      include 'C_parallel.f'
      parameter(khor1=khor+1)
      common /ind/ iwet(khor1),ldep(khor),lazc(khor),indend(n),
     *             isornr(n),isorsr(n),islab(n)
      character*12 konv11
#endif
      error=1e-3
      inumber=0
      number=0 
      wfa=1.75 
      k=1 
c----- must update hisi before calling relax
c----- first set u(2)=u(1) 
#ifdef MPI
      do 99 j=icey5,icey6
      do 99 i=icex5,icex6
#else
      do 99 j=1,ny-1 
      do 99 i=2,nx
#endif 
c----- make sure bdry pts are equal to zero 
      uice(i,j,2)=uice(i,j,1) 
      vice(i,j,2)=vice(i,j,1) 
      uice(i,j,1)=uice(i,j,3)*real(muv(i,j)) 
      vice(i,j,1)=vice(i,j,3)*real(muv(i,j)) 
 99   continue 
c----- set up coefficients of diagonal components 
      the2=2.0*theta 
#ifdef MPI
      do 102 j=icey1,icey2
      do 102 i=icex1,icex2
#else
      do 102 j=1,ny-1 
      do 102 i=2,nx
#endif 
      coefix(i,j)=0. 
      coefiy(i,j)=0. 
      do 102 im=1,muv(i,j) 
      drs5=0.5*drags(i,j) 
      dra5=0.5*draga(i,j) 
      et1=0.5*(eta(i,j)+eta(i,j+1)+eta(i-1,j)+eta(i-1,j+1)) 
      zt1=0.5*(zeta(i,j)+zeta(i,j+1)+zeta(i-1,j)+zeta(i-1,j+1)) 
      coefix(i,j)=the2*(drs5+(et1+zt1)*s1(i)+et1*s2) 
      coefiy(i,j)=the2*(drs5+(et1+zt1)*s2+et1*s1(i)) 
c
      if(abs(coefix(i,j)).lt.1.e-7)  then
       write(7,*)'coefix ist zu klein',coefix(i,j),i,j,drs5,dra5
       write(6,*)'coefix ist zu klein',coefix(i,j),i,j,drs5,dra5,
     + muv(i,j)      
       coefix(i,j)=coefix(i,j)+2.e-7
       end if
       
      if(abs(coefiy(i,j)).lt.1.e-7) then
        write(7,*)'coefiy ist zu klein',coefiy(i,j),i,j,drs5,dra5
        write(6,*)'coefiy ist zu klein',coefiy(i,j),i,j,drs5,dra5,
     +  muv(i,j)
        coefiy(i,j)=coefiy(i,j)+2.e-7
       end if
      
      coefix(i,j)=1.0/coefix(i,j) 
      coefiy(i,j)=1.0/coefiy(i,j) 
      fly(i,j)=the2*dra5*coefiy(i,j) 
      flx(i,j)=the2*dra5*coefix(i,j) 
      aax=(1.0+flx(i,j)*draga(i,j)*coefiy(i,j))
      if(abs(aax).lt.1.e-7) write(7,*)'aax ist zu klein',aax     
cc       aax=max(1.e-7,aax)
      fls(i,j)=(1.0)/aax 
 102  continue 
c----- set u(3)=u(1) 
 100  continue
#ifdef MPI
      do 121 j=icey3,icey4
      do 121 i=icex3,icex4
#else
      do 121 j=1,ny-1
      do 121 i=2,nx
#endif 
      uice(i,j,3)=uice(i,j,1)
      vice(i,j,3)=vice(i,j,1)
 121  continue

c----- begin sweep 
      call felld1(uice,vice,eta,fxe,1) 
      call felld1(uice,vice,zeta,fxz,1) 
      call felld1(vice,uice,eta,fye,1) 
      call felld1(vice,uice,zeta,fyz,1)
#ifdef MPI
c      do 103 j=ICEQ1,ICEQ2
c      do 103 i=ICEP1,ICEP2
      do 103 i=ICEP2,ICEP1,-1
      do 103 j=ICEQ2,ICEQ1,-1
#else      
      do 103 i=3,nx-1
      do 103 j=2,ny-2
#endif 
      do 103 im=1,muv(i,j) 
      k=1 
      i7=i 
      j7=j 
      call felld(uice,vice,eta,fxeta,i7,j7,k,fxe) 
      call felld(vice,uice,eta,fyeta,i7,j7,k,fye) 
      call felld(uice,vice,zeta,fxzeta,i7,j7,k,fxz) 
      call felld(vice,uice,zeta,fyzeta,i7,j7,k,fyz)
cccccccccccccccccccccccccccccccc      
      fx3=0.5*(fxeta(1)+fxzeta(1)+fxeta(2)+fxeta(3)+fxzeta(4)-fxeta(4)) 
      fx3=fx3*2.0*theta 
      fy3=0.5*(fyeta(1)+fyeta(2)+fyzeta(2)+fyzeta(3)-fyeta(3)+fyeta(4)) 
      fy3=fy3*2.0*theta
ccccccccccccccccccccccccccccc
      fxcp=amass(i,j)*theta* 
     & (uicec(i,j)*(uice(i,j+1,1)-uice(i,j-1,1))/(2.0*dln(i)) 
     & +vicec(i,j)*(uice(i-1,j,1)-uice(i+1,j,1))/(2.0*dl)) 
      fx3=fx3-fxcp 
      fycp=amass(i,j)*theta* 
     & (uicec(i,j)*(vice(i,j+1,1)-vice(i,j-1,1))/(2.0*dln(i)) 
     & +vicec(i,j)*(vice(i-1,j,1)-vice(i+1,j,1))/(2.0*dl)) 
      fy3=fy3-fycp 
ccccccccccccccccccc
      f11=(forcex(i,j)+fx3)*coefix(i,j) 
      f22=(forcey(i,j)+fy3)*coefiy(i,j) 
      uicor=(f11+flx(i,j)*f22)*fls(i,j) 
      vicor=(f22-fly(i,j)*f11)*fls(i,j) 
      uice(i,j,1)=uice(i,j,1)+wfa*(uicor-uice(i,j,1)) 
      vice(i,j,1)=vice(i,j,1)+wfa*(vicor-vice(i,j,1))
 103  continue
      number=number+1 
      
      if(number.gt.100) go to 201 
      if(number.gt.20) go to 202 
      if(number.gt.10) go to 203 
 202  continue 
      wfa=1.00 
      go to 204 
 203  continue 
      wfa=1.25 
 204  continue 
c----- check max error 
c----- form error matrix 
#ifdef MPI
      call INTI_COMM_VAR
      call com_vertical(ctemp,iel,uice(:,:,1),iel1,1,1,0,2,0,lazc)
      call com_vertical(ctemp,iel,vice(:,:,1),iel1,1,1,0,2,1,lazc)
      call com_horizontal(ctemp,iel,uice(:,:,1),iel1,1,1,0,2,0,icdepb1,
     *idepb1,icdept1
     *,idept1,ncsendb1,nsendb1,ncsendt1,nsendt1,nbhlayer1,nthlayer1)
      call com_horizontal(ctemp,iel,vice(:,:,1),iel1,1,1,0,2,1,icdepb1,
     * idepb1,icdept1
     *,idept1,ncsendb1,nsendb1,ncsendt1,nsendt1,nbhlayer1,nthlayer1)
      call MPI_Waitall(nreq1, ir_r_b, status, ierr)
      INUM1=0
      INUM2=0
      call VER_UNPACK(ctemp,iel,uice(:,:,1),iel1,1,1,0,2,lazc)
      call VER_UNPACK(ctemp,iel,vice(:,:,1),iel1,1,1,0,2,lazc)
      INUM3=0
      INUM4=0
      call HOR_UNPACK(ctemp,iel,uice(:,:,1),iel1,1,1,nbhlayer1,
     &nthlayer1, idepb1,idept1,icdepb1,icdept1,0,2,iwet)
      call HOR_UNPACK(ctemp,iel,vice(:,:,1),iel1,1,1,nbhlayer1,
     &nthlayer1,idepb1,idept1,icdepb1,icdept1,0,2,iwet)

      do 104 j=icey3,icey4
      do 104 i=icex3,icex4
      uerr(i,j)=uice(i,j,1)-uice(i,j,3)
      verr(i,j)=vice(i,j,1)-vice(i,j,3)
 104  continue
c----- find error 
      call xmaxm2(uerr,t1)
      call xmaxm2(verr,t2)
      t1=max(t1,t2)
      CALL MPI_ALLREDUCE(t1,t2,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,IERR)
      t1=t2
      if(t1.lt.error) go to 200
      go to 100
      if(myid .eq. 0) write(*,*) 'relax do not converge'
 201  continue
      number=0
      inumber=inumber+1
      do 219 j=icey5,icey6
      do 219 i=icex5,icex6
      uice(i,j,1)=uice(i,j,2)*real(muv(i,j))
      vice(i,j,1)=vice(i,j,2)*real(muv(i,j))
 219  continue
      error=1e-3*inumber
      goto 100
      
#else
      do 104 j=1,ny-1 
      do 104 i=2,nx 
      uerr(i,j)=uice(i,j,1)-uice(i,j,3)
      verr(i,j)=vice(i,j,1)-vice(i,j,3)
 104  continue
      call xmaxm(uerr,t1)
      call xmaxm(verr,t2)
      t1=max(t1,t2)
      if(t1.lt.error) go to 200
      go to 100 
 201  continue
#endif
 200  continue
#ifdef MPI
      call INTI_COMM_VAR
      call com_vertical(ctemp,iel,uice(:,:,1),iel1,2,2,0,2,0,lazc)
      call com_vertical(ctemp,iel,vice(:,:,1),iel1,2,2,0,2,1,lazc)
      call com_horizontal(ctemp,iel,uice(:,:,1),iel1,2,2,0,2,0,icdepb2,
     *idepb2,icdept2
     *,idept2,ncsendb2,nsendb2,ncsendt2,nsendt2,nbhlayer2,nthlayer2)
      call com_horizontal(ctemp,iel,vice(:,:,1),iel1,2,2,0,2,1,icdepb2,
     *idepb2,icdept2
     *,idept2,ncsendb2,nsendb2,ncsendt2,nsendt2,nbhlayer2,nthlayer2)
      call MPI_Waitall(nreq1, ir_r_b, status, ierr)
      INUM1=0
      INUM2=0
      call VER_UNPACK(ctemp,iel,uice(:,:,1),iel1,2,2,0,2,lazc)
      call VER_UNPACK(ctemp,iel,vice(:,:,1),iel1,2,2,0,2,lazc)
      INUM3=0
      INUM4=0
      call HOR_UNPACK(ctemp,iel,uice(:,:,1),iel1,2,2,nbhlayer2,
     &nthlayer2,idepb2,idept2,icdepb2,icdept2,0,2,iwet)
      call HOR_UNPACK(ctemp,iel,vice(:,:,1),iel1,2,2,nbhlayer2,
     &nthlayer2,idepb2,idept2,icdepb2,icdept2,0,2,iwet)
      do 119 j=icey5,icey6
      do 119 i=icex5,icex6
      uice(i,j,3)=uice(i,j,1)
      vice(i,j,3)=vice(i,j,1)
 119  continue
#endif
 12   format(1x,'max error and u and v power',3f12.5) 
 1    format(1x,'no of iterations are',1i5) 
      return 
      end 


*§§15 setice
      subroutine setice(iindex,izet,mm,nn)
c-----------------------------------------------------------------------
c     list of constant for ice model which must be given
c-----------------------------------------------------------------------

      include 'C_model.f'
      parameter(nx=m,ny=n)
      parameter(khor1=khor+1)   
      dimension iindex(mm,nn),izet(mm,nn)

      include 'C_ice.f' 
      dtsec=dt
c     drag coefficient for ice water stress 
      cwa = 0.0025
c     volumetric specific heat of sea-water cw (j/m**3 k)) 
c      ccw = 4.07e6         
c     latent heat of fusion of ice (j/kg) 
      hlatis = 3.36e5 
c     density of sea ice (kg/m**3) 
      rois = 930. 
c     thermal conductivity of ice (j/(m s k)) 
      tci = 2.033 
c     ice-salinity (psu) 
c      sice = 5. 
c     minimum ice-thickness 
      epsis = 1.e-7 
c     solar constant 
      solcon  = 1365.0 
c     cloud cover between 1(min) and 8(max) 
c      icloud = 5 
c     air temperature (if no data) 
       airtem = -15.0 
c     start (in hours) for shortwave radiation 
c      start = 30.*5.*24.0 * 3600. 

c-----derived constants  for ice-model 
c----- no changes here
 
c     freezing temperature of seaice 
c      tf = tfreez(sice,0.0) 
c     volumetric latent heat of fusion of ice , (j/m**3) 
      roil = rois*hlatis 

c     dt divided by latent heat of fusion of ice 
      dtroil = dtsec/roil 


      return 
      end 



*§§16 stressw
      subroutine stressw(iindex,izet,jc,zz,mm,nn,dz,iilo,tx,ty)
c-----------------------------------------------------------------------
c   calculate water stress for sea ice 
c----------------------------------------------------------------------
      include 'C_model.f'
      include 'C_index.f'
      parameter(nx=m,ny=n)
      parameter(khor1=khor+1)    

      dimension iindex(mm,nn),izet(mm,nn),jc(mm,nn)
      dimension zz(mm,nn),dz(iilo)
      dimension tx(mm,nn),ty(mm,nn)
      dimension su(m,n),sv(m,n)
#ifdef MPI
      include 'C_mpi.f'
      common /cord/ I1,I2,J1,J2,IA1,IA2,JA1,JA2,IB1,IB2,JB1
     1,ID1,ID2,JD1,JD2
      common /icecor/ ICEI1,ICEI2,ICEJ1,ICEJ2,ICEIA1,ICEIA2,ICEJA1,
     * ICEJA2,ICEX1,ICEX2,ICEY1,ICEY2,ICEP1,ICEP2,ICEQ1,ICEQ2,ICEX3,
     * ICEX4,ICEY3,ICEY4,ICEX5,ICEX6,ICEY5,ICEY6
      common /lbkhor/ khorl,lzet(khor),lb0(n),le0(n),
     * lb1(n),le1(n),lb2(n),le2(n)
#endif               
      include 'C_ice.f'
#ifdef MPI
      do k=JA1,JA2
      do i=IA1,IA2
#else
      do k=1,n
      do i=1,m
#endif
        su(i,k) = 0.
        sv(i,k) = 0.
      end do
      end do
      

#ifdef MPI
      do k = JA1,JA2       
         do lw = lb1(k),le1(k)
#else
       do k = 1,n
         do lw = lb(k),le(k)
#endif
            i = iwet(lw)
         
        kp1 = min(n,k+1)
        ip1 = min(m,i+1)     
      
        dfu = dz(1)+0.5*(zz(i,k)+zz(i,kp1))
        dfv = dz(1)+0.5*(zz(i,k)+zz(ip1,k))
       
        su(i,k) = uc(iindex(i,k)+1)/dfu
        sv(i,k) = vc(iindex(i,k)+1)/dfv
      
        if((i.lt.m).and.(k.gt.1)) then
c          ii = m-i
          gairx(i+1,k-1) = (tx(i,k-1)+tx(i+1,k-1))/2.
          gairy(i+1,k-1) = (ty(i,k)+ty(i,k-1))/2.
        endif
         end do
      end do
      
      do j=1,ny
      do i=1,nx
        gwatx(i,j)=0.
        gwaty(i,j)=0.
      enddo
      enddo

#ifdef MPI
      do k = JA1,JA2       
      do lw = lb1(k),le1(k)
#else
        do k = 1,n
         do lw = lb(k),le(k)
#endif
            i = iwet(lw)
            nwet = indwet(lw)
        kp1 = min(n,k+1)
        ip1 = min(m,i+1)
      
        dfu = dz(1)+0.5*(zz(i,k)+zz(i,kp1))
        dfv = dz(1)+0.5*(zz(i,k)+zz(ip1,k))
      
        uwi = ui(i,k)-su(i,k)
        vwi = vi(i,k)-sv(i,k)
        rdich = (stc(nwet+1)+rhoq(1))*1000.
 
        if((i.lt.m).and.(k.gt.1)) then
          uu = (su(i,k-1)+su(i+1,k-1))/2.
          vv = (sv(i,k)+sv(i,k-1))/2.
c	  ii = m-i
          gwatx(i+1,k-1) = uu
          gwaty(i+1,k-1) = vv
        endif
     
        txi = rdich*cwa*sqrt(uwi**2+vwi**2)*uwi
        tyi = rdich*cwa*sqrt(uwi**2+vwi**2)*vwi

        txc(lw) = frice(i,k)*txi+(1.0-frice(i,k))*txc(lw)
        tyc(lw) = frice(i,k)*tyi+(1.0-frice(i,k))*tyc(lw)
c-------------test test Ute 180917 das wars nicht
c        txc(lw) = txc(lw)
c        tyc(lw) = tyc(lw)

         end do
      end do
      return 
      end


*§§18 stressa
      subroutine stressa(iindex,izet,jc,zz,mm,nn,dz,iilo,tx,ty)
c-----------------------------------------------------------------------
      include 'C_model.f'

      parameter(nx=m,ny=n)
      parameter(khor1=khor+1)    

      dimension iindex(mm,nn),izet(mm,nn),jc(mm,nn)
      dimension zz(mm,nn),dz(iilo)
      dimension tx(mm,nn),ty(mm,nn)
      dimension su(m,n),sv(m,n)
               
      include 'C_ice.f'

      do 323 k=1,n
      do 323 i=1,m
        su(i,k) = 0.
        sv(i,k) = 0.
        do 323 l=1,jc(i,k)
         
          kp1 = min(n,k+1)
          ip1 = min(m,i+1)     
      
          dfu = dz(1)+0.5*(zz(i,k)+zz(i,kp1))
          dfv = dz(1)+0.5*(zz(i,k)+zz(ip1,k))
       
          su(i,k) = uc(iindex(i,k)+1)/dfu
          sv(i,k) = vc(iindex(i,k)+1)/dfv
     
          if((i.lt.m).and.(k.gt.1)) then
c            ii = m-i
            gairx(i+1,k-1) = (tx(i,k-1)+tx(i+1,k-1))/2.
            gairy(i+1,k-1) = (ty(i,k)+ty(i,k-1))/2.
          endif
 323  continue      
      
      do j=1,ny
      do i=1,nx
        gwatx(i,j) = 0.
        gwaty(i,j) = 0.
      enddo
      enddo
       
      do 324 k=1,n
      do 324 i=1,m
      do 324 l=1,jc(i,k)

        kp1 = min(n,k+1)
        ip1 = min(m,i+1)
      
        dfu = dz(1)+0.5*(zz(i,k)+zz(i,kp1))
        dfv = dz(1)+0.5*(zz(i,k)+zz(ip1,k))
       
        uwi = ui(i,k) - su(i,k)
        vwi = vi(i,k) - sv(i,k)
        rdich = (stc(iindex(i,k)+1)+rhoq(1))*1000.
 
        if((i.lt.m).and.(k.gt.1)) then
          uu = (su(i,k-1)+su(i+1,k-1))/2.
          vv = (sv(i,k)+sv(i,k-1))/2.
c	  ii = m-i
        gwatx(i+1,k-1) = uu
        gwaty(i+1,k-1) = vv
        endif
        txi = rdich*cwa*abs(uwi)*uwi
        tyi = rdich*cwa*abs(vwi)*vwi
  324 continue
      return 
      end


*§§20 xmaxm
      subroutine xmaxm(uerrv,x) 
c-----------------------------------------------------------------------
c     find maximum
c-----------------------------------------------------------------------
      include 'C_model.f'

      parameter (nx=m,ny=n,n4=nx*ny) 
      dimension uerrv(nx,ny) 
#ifdef MPI
      common /icecor/ ICEI1,ICEI2,ICEJ1,ICEJ2,ICEIA1,ICEIA2,ICEJA1,
     * ICEJA2,ICEX1,ICEX2,ICEY1,ICEY2,ICEP1,ICEP2,ICEQ1,ICEQ2,ICEX3,
     * ICEX4,ICEY3,ICEY4,ICEX5,ICEX6,ICEY5,ICEY6
#endif
      x = 0.0 
#ifdef MPI
      do 100 j=ICEY1,ICEY2
      do 100 i=ICEX1,ICEX2
#else
      do 100 j=1,ny-1 
      do 100 i=2,nx
#endif
      x = max(abs(uerrv(i,j)),x) 
 100  continue 

      return 
      end                                


*§§21 xsum
c-----------------------------------------------------------------------
      subroutine xsum(hisi1,fricei1,ss1,ss2,ss3) 
c-----------------------------------------------------------------------
c     sum up a variable
c-----------------------------------------------------------------------
      include 'C_model.f'

      parameter (nx=m,ny=n,n3=n*m) 
      dimension hisi1(nx,ny),fricei1(nx,ny) 

      ss1 = 0.0 
      ss2 = 0.0 
      ss3 = 0.0 
      do 100 j=1,ny 
      do 100 i=1,nx
        ss1 = ss1+hisi1(i,j)*fricei1(i,j) 
        ss2 = ss2+fricei1(i,j) 
        ss3 = ss3+hisi1(i,j) 
 100  continue 
      return 
      end 


*§§4 fcuv
c-----------------------------------------------------------------------
      subroutine fcuvs(fricei1) 
c-----------------------------------------------------------------------
c     determine the u, v and h information field
c-----------------------------------------------------------------------
      include 'C_model.f'
      parameter (nx=m,ny=n) 
      dimension fricei1(m,n) 
      common/array/mhc(m,n),muv(nx,ny),maphc(m,n)
       
      do j=1,n 
      do i=1,m
        mhc(i,j)=maphc(i,j) 
        if(maphc(i,j).le.0) mhc(i,j)=9 
      enddo
      enddo

      do j=2,n-1 
      do i=2,m-1
        a0=fricei1(i-1,j)+fricei1(i+1,j)
     *    +fricei1(i,j-1)+fricei1(i,j+1) 
        if(mhc(i,j).gt.8) goto 101
        if(a0.lt.1.0e-5) mhc(i,j) = 0 
 101    continue 
      enddo
      enddo

      do j=1,ny-1 
      do i=2,nx
        mr1 = mhc(i,j)+mhc(i,j+1)+mhc(i-1,j)+mhc(i-1,j+1) 
        if(mr1.ge.9.or.mr1.eq.0)then 
          muv(i,j) = 0 
        else 
          muv(i,j) = 1 
        endif 
      enddo
      enddo

      do j=1,n 
      do i=1,m
        if(mhc(i,j).eq.9) mhc(i,j) = 0 
      enddo
      enddo
      return 
      end 
*§§20 xmaxm
      subroutine xmaxm2(uerrv,x)
c-----------------------------------------------------------------------
c     find maximum
c-----------------------------------------------------------------------
      include 'C_model.f'

      parameter (nx=m,ny=n,n4=nx*ny) 
      dimension uerrv(nx,ny)
#ifdef MPI
      common /icecor/ ICEI1,ICEI2,ICEJ1,ICEJ2,ICEIA1,ICEIA2,ICEJA1,
     * ICEJA2,ICEX1,ICEX2,ICEY1,ICEY2,ICEP1,ICEP2,ICEQ1,ICEQ2,ICEX3,
     * ICEX4,ICEY3,ICEY4,ICEX5,ICEX6,ICEY5,ICEY6
#endif
      x = 0.0 
#ifdef MPI
      do 100 j=ICEY3,ICEY4
      do 100 i=ICEX3,ICEX4
#else
      do 100 j=1,ny-1
      do 100 i=2,nx
#endif
      x = max(abs(uerrv(i,j)),x) 
 100  continue 

      return 
      end        




