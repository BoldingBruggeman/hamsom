c-----------------------------------------------------------------------
c     somice declaration of common-blocks used in ice modell
c-----------------------------------------------------------------------
      common /gitter/ dt,r,g,dl,dlr,dlrh,dln(m),rdln(m),rad,dth,
     *  dlvo(m),dlvu(m),gh,rdt,dt2,r4,coru(m),corv(m),dtrdln(m),dtdlr

      common /ind/ iwet(khor1),ldep(khor),lazc(khor),indend(n),
     *             isornr(n),isorsr(n),islab(n)

      common uc(ndrei),vc(ndrei),stc(ndrei),avc(ndrei),z(m,n) 
      common zac(khor),wobc(khor),stuvc(khor),fricu(khor),cxc(khor) 
      common cyc(khor),pac(khor),txc(khor),tyc(khor) 
      common stpc(ndrei),sac(ndrei),tec(ndrei) 
      common pres(ilo),wc(ndrei),fricv(khor) 

      common /ice/ cwa,ccw,tci,sice,rois,cp,tf,epsis,roil,dtroil,
     *  his(m,n),frice(m,n),tis(m,n),tair(m,n),qoc(m,n),icloud,
     *  solcon,frsw,frsi,time,xlat(m),ui(m,n),vi(m,n),us(m,n),vs(m,n)
     
      common /icevelo/ uice(nx,ny,3),vice(nx,ny,3),
     *  uicec(nx,ny),vicec(nx,ny),uerr(nx,ny),verr(nx,ny),
     *  gairx(nx,ny),gairy(nx,ny),gwatx(nx,ny),gwaty(nx,ny),
     *  drags(nx,ny),draga(nx,ny),amass(nx,ny),
     *  eta(m,n),zeta(m,n),hisir(m,n,3),hisi(m,n,3),fricei(m,n,3),
     *  forcex(nx,ny),forcey(nx,ny)

      common /stresi/ stresi(m,n,3)
      common /more/ diff1,lad,icount,nst6
      common /step/ delice,error,deltt
      common /array/ mhc(m,n),muv(nx,ny),maphc(m,n)
      common /met/ windx,windy,wlam,pnull,rhoq(ilo),refrho(ilo),stress
      common /hh/ hisr(m,n)

 
