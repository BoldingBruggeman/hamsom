*===============  subroutines to be changed ================================
*��1 kotief               !specify open boundary lines
*��2 gatmit
*��3 readzdv
*��4 read_bconditions
*��6 liq_boundary         !set open boundary lines  
*��7 tsrneu               !specify open boundary lines
*��9 uvrand               !specify open boundary lines
*��10 bforce              !specify open boundary lines
*��11 zetadb              !specify open boundary lines
*��12 readtsin
*��13 topocorr




*��1 kotief
c-----------------------------------------------------------------------
      subroutine kotief (dz,tz,nhor,ktot) 
c-----------------------------------------------------------------------
c     estimation of topography dependent arrays 
c-----------------------------------------------------------------------
      include 'C_model.f'
      include 'C_index.f'
      include 'C_mpi.f' 
      parameter(khor1=khor+1) 
      integer tz(m,n) 
      dimension dz(ilo) 
 
      common /ind/ iwet(khor1),ldep(khor),lazc(khor), 
     *  indend(n),isornr(n),isorsr(n),islab(n) 
 
      idm = 1 
      nhor = 0 
      ktot = 0 
      do 6 k = 1,n 
      do 5 i = 1,m 
      it = tz(i,k) 
      itt = it 
      if (it.le.0) goto 5 
      do 2 j = 1,ilo 
      jj = j 
      if (it.le.dz(j)) goto 3 
    2 continue 
      it = dz(ilo) 
      tz(i,k) = it 
      IF(myid .EQ. 0)THEN
      write (7,603) tz(i,k),itt,i,k 
  603 format (2x,'gitief, depth reduction to  dz(ilo),', 
     2'  tneu, talt, (i,k)/',2x,2i8,4x,2i4)
      endif 
    3 iz = 0 
      j = jj 
      if (j.gt.1) iz = dz(j-1) 
      idep = it-iz 
      if(idep-idm+1) 45,44,4 
   44 tz(i,k) = tz(i,k)+1 
      it = tz(i,k)
      IF(myid .EQ. 0)THEN 
      write (6,602) tz(i,k),itt,i,k
      endif 
      goto 3 
   45 j = j-1 
      jj = max  (j,0) 
      tz(i,k) = iz 
      it = iz 
      IF(myid .EQ. 0)THEN
      write (7,602) tz(i,k),itt,i,k 
  602 format (2x,'gitief, depth changed, tneu, talt, (i,k)/', 
     22x,2i8,4x,2i3)
      endif 
      if (jj-1) 5,3,3 
    4 nhor = nhor+1 
      ktot = ktot+jj 
      ldep(nhor) = idep 
      iwet(nhor) = i 
      lazc(nhor) = jj 
    5 continue 
      indend(k) = nhor 
      islab (k) = ktot 
    6 continue 
c     -----------------beginning and end- indicees for sor-slabs 


      lwe = 0 
      isor=3 
      do 9 k = 1,n 
*     if(k.gt.100) isor=2 
      lwa = lwe+1 
      lwe = indend(k) 
      lxan=-2 
      do 7 lw = lwa,lwe 
      i = iwet(lw) 
      lxan=lw 
c--------------------------------------------------------------------
c      Attention!!!!  this has to be adapted to new topography
c      exclude wet open boundary lines from sor iteration
c  N and W boundaries:inner First and Second lines to be excluded
c  S and E boundaries:inner First lines to be excluded
c--------------------------------------------------------------------
      if((k.le.17).and.(i.gt.150)) goto 7 ! has to be changed 
      if((i.le.67).and.(k.le.57 )) goto 7 ! has to be changed

      if (i.ge.isor) goto 88 
    7 continue 
   88 isornr(k) = lxan 
      lwp = lwe+lwa 
      lxen=-4 
      do 8 ll = lwa,lwe 
      lxen=lwp-ll 
      lw = lwp-ll 
      i = iwet(lw) 
c--------------------------------------------------------------
c      Attention!!!!  this has to be adapted to new topography
c      exclude wet open boundary lines from sor iteration
c      the same as above 
c--------------------------------------------------------------
      if((k.le.17).and.(i.gt.150)) goto 8  !has to be changed
      if((i.le.67).and.(k.le.57 )) goto 8  !has to be changed


      if (i.le.m-1) goto 99 
    8 continue 
   99 isorsr(k) = lxen 
    9 continue 
c
      nwet = 0 
      lwe = 0
      do j = 1,n 
         lwa = lwe+1 
         lwe = indend(j) 
         lb(j) = lwa
         le(j) = lwe
         do lw = lwa,lwe 
            indwet(lw) = nwet
            jwet(lw) = j
            do k = 1, lazc(lw)
               indi(nwet+k) = iwet(lw) 
               indver(nwet+k) = k
               llw(nwet+k)=lw
            end do
            nwet = nwet+lazc(lw) 
         end do
      end do
      IF(myid .EQ. 0)THEN
      write (7,600) nhor,ktot 
      write (7,'(a2,15(1x,i5))')'nr',(isornr(j),j=1,n) 
      write (7,'(a2,15(1x,i5))')'sr',(isorsr(j),j=1,n) 
 600  format (/1x,'number of wet grid points/cells (2d,3d) /',2i8) 
      endif
		k2=0
		do k=1,n
		do i=isornr(k),isorsr(k)
		k2=k2+1
		  enddo
		  enddo
      IF(myid .EQ. 0)THEN
      write (7,*)'//////' 
      write (7,*)' KSOR  ',k2 
      write (*,*)' KSOR  ',k2 
      endif
      return 
      end

*��2 gatmit
c-----------------------------------------------------------------------
      subroutine gatmit (deltat) 
c-----------------------------------------------------------------------
c     time step dependent and grid dependent variables  
c----------------------------------------------------------------------- 
      include 'C_model.f'

      parameter(ilop1=ilo+1) 
      parameter(khor1=khor+1) 
 
      common uc(ndrei),vc(ndrei),stc(ndrei),avc(ndrei),z(m,n) 
      common zac(khor),wobc(khor),stuvc(khor)
      common fricu(khor),cxc(khor) 
      common cyc(khor),pac(khor),txc(khor),tyc(khor) 
      common stpc(ndrei),sac(ndrei),tec(ndrei) 
      common pres(ilo),wc(ndrei),fricv(khor) 
      common /dreh/ sinfu(m),sinfv(m),cosfu(m),cosfv(m),sincx(m), 
     *  sincy(m),bx(m),by(m),pxu(m),pyv(m),pyu(m),pxv(m) 
      common /gitter/ dt,r,g,dl,dlr,dlrh,dln(m),rdln(m),rad,dth, 
     *  dlvo(m),dlvu(m),gh,rdt,dt2,r4,
     *  coru(m),corv(m),dtrdln(m),dtdlr 
      common /gitter2/ fianf,yambdanf,dphi,dlambda,rearth
      common /num/ dc(ilo),av(ilo),ad(ilo),dh(ilo),pd(ilo), 
     *  prd(ilo),pr2d(ilo),r2d(ilo),tkw(ilo),tau(ilo),dd(ilo), 
     *  qa(ilop1),qbet(ilop1),qn(ilop1),rd(ilop1) 
      common /met/ windx,windy,wlam,pnull,
     *  rhoq(ilo),refrho(ilo),stress 

      include 'C_xy.f'

      g = 9.81   ! defined in program main
      rearth=6350000.


*---------------------------------------------------------------------
*		set XY grid for your model domain           
*         print out in control file #7:  pathpc//pathout//control
*---------------------------------------------------------------------

*____model north_b _________________________________
c---- set  first point  i=1, j=1 --------------
c     northern boundary  phi= 65 grad 56. min  N 
       f = 65.+59./60.  

c     western boundary  lam= -4 grad 10. min (W)
       y = -4.0-10./60.  ! 


c---- set  XY step im minutes
		phimin= 6.					 ! phi
		yammin=10.					 ! lam
*--------------------------------------------------------------
*--------------------------------------------------------------


               fianf = f 
             yambdanf= y

*		XY step in degree 
		phigrad  = phimin/60.			
		yamgrad  = yammin/60.



          do i=1,m
          xt(i) = fianf-real(i-1)*phigrad			  !phi = 1 min
	    cosxt(i)=cos(xt(i)*rad)
          enddo

        do j=1,n
        yt(j) = yambdanf+real(j-1)*yamgrad            !lam = 2 min
        enddo

	do j=1,n
        yu(j) = yambdanf+real(j-1+0.5)*yamgrad            !lam = 2 min
	enddo

	do i=1,m
          xv(i) = fianf-real(i-1+0.5)*phigrad			  !phi = 1 min
	enddo

c     gezeit - und gitterkonstanten 
      dt = deltat 
      dth = dt/3600. 
      rdt = r*dt 
      r4 = 0.25*r 
      dt2 = dt*0.5 
      gh = 0.5*g 
      gdt = g*dt 
      pi = 4.*atan(1.) 
      rad = pi/180. 
      aeqgrd = 2.*pi*6.37104e6/360.   ! 111200 [m/grad]

          do i=1,m
	    cosxt(i)=cos(xt(i)*rad)
          enddo


      dl     = phigrad *aeqgrd   !step phi in [m] 
      dlam   = yamgrad *aeqgrd   !step lam in [m]
      dphi   =phigrad * rad   !step phi in rad
      dlambda=yamgrad * rad   !step lam in rad


      dlr = 1./dl 
      dtdlr = dt*dlr 
      dlrh = 0.5*dlr 
      omto = 4.*pi/86164. 
 
 
      fiu = fianf * rad + dphi 
      fiv = fiu-.5*dphi 
        write(7,*)' j,     xv(j),     dlvo(j),     dlvu(j)'

      do j=1,m 
        fiu = fiu-dphi 		        ! 1=fianf,2,3,4   u-point
        fivo = fiv 					! u(i-1)=fianf-1/2 North
        fiv = fiv-dphi                ! u (i) =fianf+1/2 South
        cosfiu = cos(fiu) 
        dlvu(j) = cos(fiv)/cosfiu  !u (i)
        dlvo(j) = cos(fivo)/cosfiu !u(i-1)
        write(7,*)j,xv(j),dlvo(j),dlvu(j)
        dln(j) = dlam*cosfiu            ! at T(j) point
        rdln(j) = 1./dln(j) 
        dtrdln(j) = dt*rdln(j) 
        dlnv = dlam*cosfiu*dlvu(j) 
        fu = omto*sin(fiu) 
        fv = omto*sin(fiv) 
        coru(j) = fu 
        corv(j) = fv 
        fudt = fu*dt 
        fvdt = fv*dt 
        sinfu(j) = sin(fudt) 
        sinfv(j) = sin(fvdt) 
        cosfu(j) = cos(fudt) 
        cosfv(j) = cos(fvdt) 

c       sincx and sincy for variable time averageing in s.o.r. 
        sincx(j) = gdt*sinfu(j)*rdln(j)**2 / fu 
        sincy(j) = gdt*sinfv(j)*dlr**2 / fv 
        bx(j) = dln(j)*(sin(fudt*0.5))**2/(2.*dl*sinfu(j)) 
        by(j) = -dl*(sin(fvdt*0.5))**2/(2.*dlnv *sinfv(j)) 
 
c       corriolis rotation after Roland Wais
c       application on pressure gradient terms 
 
        pxu(j) = sinfu(j)*rdln(j)/fu 
        pyu(j) = (sin(fudt*0.5))**2*dlr/(2.*fu) 
        pyv(j) = sinfv(j)*dlr/fv 
        pxv(j) = -(sin(fvdt*0.5))**2/(2.*dlnv*fv) 
      enddo

      return 
      end 


*��3 readzdv
c-----------------------------------------------------------------------
      subroutine readzdv(zval,lrp,nzac)
c-----------------------------------------------------------------------
c     reads 1-h tidal boundary conditions for sea surface elevations
c-----------------------------------------------------------------------
       dimension zval(lrp)

       read(nzac,'(8(1x,e12.6))') zval

c      read(nzac,2004,err=123) zvar
c2004  format(8(1x,e12.6))

c     read(nzac)zvar
      goto 1
 123  write(*,*)'error in reading boundary values z_v'
      pause

1     return
      end     



*��3 readzdv
c-----------------------------------------------------------------------
      subroutine readzdv_30min(zvar,lrp,nzac)! dhanya 14.12.09
c-----------------------------------------------------------------------
c     reads 1-h tidal boundary conditions for sea surface elevations
c-----------------------------------------------------------------------
       dimension zvar(lrp)
*ASCII
c      read(nzac,2004,err=123) zvar
c2004  format(10f8.5)

*Binary
c      read(nzac)zvar

      goto 1
 123  write(*,*)'error in reading boundary values z_v'
      pause

1     return
      end     



*��4 read_bconditions
c-----------------------------------------------------------------------
      subroutine read_bconditions (nfile,nzactp,nday)
c-----------------------------------------------------------------------
c     boundary conditions for zd (SSE daily variations); and monthly T and S
c-----------------------------------------------------------------------
      include 'C_model.f'
      include 'C_mpi.f' 
      common /tscdb/ sarum(lvrp),tarum(lvrp)    !,bio_rum(lvrp,nbio)
      common /tsc/ sarum_old(lvrp),tarum_old(lvrp)
      common /urand/       zday(lrp),zvar(lrp),zalt(lrp)
      common /uranddb/ zday_old(lrp),zvar_old(lrp)
      common /liq_b/lq_b(iranz,5),lq_b2D(iranz*m),lq_b3D(iranz*m*ilo),
     &  num_b(7),nul_p(6)



c*     Z_D
c		do l=1,lrp
c		zday_old(l)=zday(l)
c		enddo
c		do l=1,lvrp
c		sarum_old(l)=sarum(l)
c		tarum_old(l)=tarum(l)
c		enddo

      read(nzactp,1011,err=123) nwjah,nwmon,nwday,(zday(l),l=1,lrp)
 1011   format( i4,2i2,39(/15f8.5) )
      IF(myid .eq. 0)
     &write(*,*)' boundary values z_d  were read, nday=', nwday,
     &maxval(zday)
         if(nday.ne.nwday)goto 123

*     T  S      !north_b
      if(nday.eq.1)then

      read(nfile,*,err=124) nwday
      do l=1,lvrp
      read(nfile,'( 2(1x, e12.6) )',err=124) tarum(l), sarum(l)
      enddo
      IF(myid .eq. 0)
     &print*,' boundary values t,s  were read, nday=', nwday
     &,maxval(sarum),maxval(tarum)
         if(nday.ne.nwday)goto 123
      endif
      goto 1
      
 123  write(*,*)'error in reading boundary values z_d, nday=',nday,nwday
      pause
	goto 1

 124  write(*,*)'error in reading boundary values T_S, nday=',nday,nwday
      pause
      
1     return
      end    


*��6 liq_boundary
      subroutine liq_boundary(nmodel)
! lq_b(iranz,5) open boundary i,j positions; lq_b(num,i) num=1,5; i=i1,i2,j1,j2,dirflow(+/-1) 
! f2D lq_b2D(lrp)     lw  -index for 2D array of liquid boundaries i,j positons
! f3D lq_b3D(lvrp)    nwet-index for 3D array of liquid boundaries i,j,k positons
 
      include 'C_model.f'
      include 'C_paths.f'
      parameter (khor1=khor+1)
      
	  common /ind/ iwet(khor1),ldep(khor),lazc(khor),indend(n),
     *             isornr(n),isorsr(n),islab(n)

      common /gitter2/ fianf,yambdanf,dphi,dlambda,rearth
      common /jwet/jwet(khor)	

      common /liq_b/lq_b(iranz,5),lq_b2D(iranz*m),lq_b3D(iranz*m*ilo),
     &  num_b(7),nul_p(6)

       integer num_n(7)         ! num_n is the number of points in each l- boundary 
	                           ! num_b is the total number from boundary to boundary
 

      character nmodel*7
      include 'C_xy.f'

      nul_p(1)=4
      nul_p(2)=5
      nul_p(3)=m-2
	  nul_p(4)=4
	  nul_p(5)=5
      nul_p(6)=n-2


	      !n2D=   585  n3D=  5789
		num=iranz
		do i=1,m*num
		lq_b2D(i)=0
		enddo

		do j=1,m*ilo*num
		lq_b3D(j)=0
		enddo

c-------------------------------------------------------------
c Attention!!!!  Specify your boundary lines (num=1,2,..,iranz)
c  N and W boundaries:inner Second lines to be specified
c  S and E boundaries:inner First lines to be specified
c-------------------------------------------------------------
*------ north_b--------------------------
*1 western boundary: indexes for external line: i=156:164 and j=16:16 
		num=1
		lq_b(num,1)=  156  !  i1
		lq_b(num,2)=  164  !  i2
		lq_b(num,3)=17     !  j1 inner second line to be specified
		lq_b(num,4)=17     !  j2 inner second line to be specified
		lq_b(num,5)=+1     !  direction
*2 northern boundary: indexes for external line: i=66:66 and j=12:57
		num=2
		lq_b(num,1)=67    !  i1 inner second line to be specified
		lq_b(num,2)=67    !  i2 inner second line to be specified
		lq_b(num,3)=12	  !  j1
		lq_b(num,4)=57    !  j2
		lq_b(num,5)=+1    !  direction
*-----------------------------------------

	do iii=1,7
	   num_b(iii)=0
	   num_n(iii)=0
      enddo

  	        n2D=0
	        n3D=0
	        nn2D=0
   	        DO num=1,iranz
      lwe = 0 
      nwet = 0 
 	  do 3 j=1,n
        lwa = lwe+1 
        lwe = indend(j) 
          do 4 lw=lwa,lwe 
	    jwet(lw)=j
          lump = lazc(lw) 
          i = iwet(lw) 
			 do inum=lq_b(num,1),lq_b(num,2)	  
               do jnum=lq_b(num,3),lq_b(num,4)

	 		  if(inum.eq.i.and.jnum.eq.j)then
                n2D=n2D+1 
                lq_b2D(n2D)= lw
			  nwet2=nwet
			  do k=1,lump 
                nwet2= nwet2+1 
	            n3D=n3D+1
                lq_b3D(n3D)= nwet2             
			  enddo
		        nn2D=nn2D+1
			  endif

   		     enddo
               enddo
	   do k=1,lump 
	   nwet=nwet+1
	   enddo
4        continue  ! j,i cycles
3        continue  ! j,i cycles
	   num_b(num)=nn2D
	        ENDDO	!	 do num=1,iranz
	   num_b(6)=n2D
	   num_b(7)=n3D

		GOTO 10  !bln prepering and liq_b output
	    open(121,file=pathpc2//nmodel//'/coastline/lq_check.dat')
	    open(122,file=pathpc2//nmodel//'/coastline/lq_ij.dat')
	    open(123,file=pathpc2//nmodel//'/coastline/lq_mline.bln')


		lnum=num_b(1)
		in1=1
	
	   	        DO num=1,iranz
	    lnum=0
		in2=num_b(num)
		num_n(num) =in2-in1+1
		write(123,*)num_n(num)
		do lll=in1,in2
		lnum=lnum+1
		lw=lq_b2D(lll)
		i = iwet(lw) 
		j=  jwet(lw)
		 lump = lazc(lw) 

       write(122,'( 4(1x,i4),2(1x,f11.6) )')lll,i,j,lump,yt(j),xt(i)

				enddo   !do lll=1,num_b(num)
c          	write(123,*)yt(jwet(lq_b2Din1)),xt(iwet(f2D(in1))
			in1=num_b(num)+1

	        ENDDO	!	 do num=1,iranz

		write(121,'(3(a5,i6,1x))')' n2D=',n2d,' n3D=',n3D,'nwet=',nwet
		write(121,'(a42)')'    1     2     3     4     5    S2D   S3D'  
		write(121,'(7(1x,i5))') (num_b(num),num=1,7)
		write(121,'(7(1x,i5))') (num_n(num),num=1,7)
          write(121,'(122(1x,i6))')(lq_b3D(lw),lw=1,num_b(7))
	  stop

10	  continue
      return
      end


*7 tsrneu
c-----------------------------------------------------------------------
        subroutine tsrneu (iindex,izet,time,dz,delsal)
c-----------------------------------------------------------------------
c    boundary conditions for temperature and salinity
c-----------------------------------------------------------------------

        include 'C_model.f'

		  parameter(ilop1=ilo+1) 
		  parameter(khor1=khor+1) 
        dimension f(m,n,ilo)
        dimension dz(ilo)
        dimension iindex(m,n),izet(m,n)
	  

      common /tscdb/ sarum(lvrp),tarum(lvrp) !,bio_rum(lvrp,nbio)
 	common /tsc/ sarum_old(lvrp),tarum_old(lvrp)

        common /liq_b/lq_b(iranz,5),lq_b2D(iranz*m),lq_b3D(iranz*m*ilo),
     &  num_b(7),nul_p(6)

      common /gitter/ dt,r,g,dl,dlr,dlrh,dln(m),rdln(m),rad,dth,
     *  dlvo(m),dlvu(m),gh,rdt,dt2,r4,coru(m),corv(m),dtrdln(m),dtdlr

      common /ind/ iwet(khor1),ldep(khor),lazc(khor),indend(n),
     *             isornr(n),isorsr(n),islab(n)
      common uc(ndrei),vc(ndrei),stc(ndrei),avc(ndrei),z(m,n) 
      common zac(khor),wobc(khor),stuvc(khor),fricu(khor),cxc(khor) 
      common cyc(khor),pac(khor),txc(khor),tyc(khor) 
      common stpc(ndrei),sac(ndrei),tec(ndrei) 
      common pres(ilo),wc(ndrei),fricv(khor) 
      real delsal(6)
      integer isa
c--------------------------------------------------------------------
c    setting of new boundary values
c--------------------------------------------------------------------
* S and T values are prescribed to i,j indices of boundary lines defined in lq_b3D(l)
* with respect to their specifications in 'subroutine liq_boundary'
         do l=1,lvrp   
	   sac(lq_b3D(l))=sarum(l) 
c	   tec(lq_b3D(l))=tarum(l)
         enddo

c--------------------------------------------------------------------
c    setting of equal boundary lines
c--------------------------------------------------------------------

c     Westrand
      do i=156,164
        l   = iindex(i,16)  !First line
        ll  = iindex(i,17)  !Second line
c----changed by Ute Daewel 26.4.2010
	lll = iindex(i,18)  !third line to right T from inside to boundary

        lw  = izet(i,16)
	               lump=0
        if(lw.gt.0)lump = lazc(lw)
        do k=1,lump
c----------temperature from inside to boundary
	  tec (ll+k)  = tec (lll+k)	 
          tec (l +k)  = tec (ll+k)
c---if transport is directed outside the domain S is set from the inside
          if(uc(ll+k).lt.0.)then
          sac (ll+k)  = sac (lll+k)
          endif
          sac (l +k)  = sac (ll+k)
        enddo
      enddo

c     Nordrand 
      do j=12,57
        l   = iindex(66,j)  !First line
        ll  = iindex(67,j)  !Second line
	lll = iindex(68,j)  !third line to right T from inside to boundary

        lw  = izet(66,j)
	             lump=0
        if(lw.gt.0)lump = lazc(lw)
        do k=1,lump
c----incorporate anual salinity anomaly Ute Daewel
	if (k.le.9)then
	if (j.ge.12.and.j.lt.35)then
	isa=1
	elseif(j.ge.35.and.j.lt.47)then
	isa=2
	elseif(j.ge.47.and.j.lt.58)then
	isa=3
	endif
	elseif (k.gt.9)then
        if (j.ge.12.and.j.lt.35)then
        isa=4
        elseif(j.ge.35.and.j.lt.47)then
        isa=5
        elseif(j.ge.47.and.j.lt.58)then
        isa=6
        endif
	endif


c----------temperature from inside to boundary
          tec (ll+k)  = tec (lll+k)     
          tec (l +k)  = tec(ll+k)
c-----if v is directed northward, the salinity from the inside is written to the boundary points
          if (vc(ll+k).gt.0.)then
          sac(ll+k)=sac(lll+k)
          else
	  sac (ll+k)  = sac(ll+k)+delsal(isa)  !Ute annual variability from ICES data
          endif
c	if(sac(ll+k).gt.35.5)then
c	print*,'rand sac',ll+k,sac(ll+k)
c	endif
c	if(sac(ll+k).gt.35.5)then
c	print*,'rand sac',ll+k,sac(ll+k)
c	endif
c	  sac (ll+k)  = sac(ll+k)
          sac (l +k)  = sac (ll+k)
       enddo
      enddo
 


	    goto 5 		  ! output
	     lwe = 0
         nwet=0
	    do j=1,n
        lwa = lwe+1
        lwe = indend(j)
        do lw=lwa,lwe
        i = iwet(lw)
               lump = lazc(lw)
        do k=1,lump
        nwet = nwet+1
         f(i,j,k) = tec (nwet)
		enddo
        enddo
        enddo

		open(43,file='test_T_tsr')
        if(mod(i-1,10).eq.0)write(43,'(3i3, 9(1x,i6))')
     &	    1,k, i,( (j), j=4,n-2)
		do k=1,ilo
	    write(43,*)'k=',k,'////, timestep=',time/(60*10)+1
	    write(43,'(3i3, 9(1x,i6))')
     &	    1,k, i,( (j), j=4,n-2)
		do i=4,m-2
		write(43,'(3i3, 9(1x,f6.2))')
     &	    1, k, i,( f(i,j,k), j=4,n-2)
        if(mod(i-1,11).eq.0)write(43,'(3i3, 9(1x,i6))')
     &	    1,k, i,( (j), j=4,n-2)
	enddo
	enddo
		close(43)
		write(*,*)'TSRneu output'
	
c		pause
5      return 
       end 

#ifdef _FABM_

*��5 read_bio_bconditions
c-----------------------------------------------------------------------
      subroutine read_bio_bconditions (nfile,nfile2,nday,ijulu,month)
c-----------------------------------------------------------------------
c     boundary values for bio-variables 
c-----------------------------------------------------------------------
      include 'C_model.f'
      include 'C_nbio_fabm.f'
      include 'C_mpi.f'
      parameter(khor1=khor+1)
      common /ind/ iwet(khor1),ldep(khor),lazc(khor), 
     *  indend(n),isornr(n),isorsr(n),islab(n) 
      
      common /biobd/ bio_rum(lvrp,nbio)
c      common /tscdb/ sarum(lvrp),tarum(lvrp),bio_rum(lvrp,nbio)
c	common /tsc/ sarum_old(lvrp),tarum_old(lvrp)
      common /urand/       zday(lrp),zvar(lrp),zalt(lrp)
      common /uranddb/ zday_old(lrp),zvar_old(lrp)
      common /liq_b/lq_b(iranz,5),lq_b2D(iranz*m),lq_b3D(iranz*m*ilo),
     &  num_b(7),nul_p(6)
      dimension pl(lrp,2),ps(lrp,2)


      if(nday.eq.1)then
      read(nfile,*,err=124) nwmonht
      do l=1,lvrp
      read(nfile,'( 11(1x, e12.6) )',err=124) (bio_rum(l,i),i=3,13)
       bio_rum(l,14)=bio_rum(l,12)/20.  !...opal....
      enddo
      IF(myid .eq. 0)
     &print*,' boundary values bio_bd  were read, month=',nwmonht
      call change_units_bd(bio_rum,1)

      goto 1
      

 124  write(*,*)'error in reading boundary values bio_bd, month=',
     & nwmonht
      pause
      endif

1     continue
      goto 2
      if(mod(ijulu-1,5).eq.0.and.ijulu.ne.366)then
	write(*,*)'read bio boundary day,ijulu ',nday,ijulu
      do i=10,lrp
      read(nfile2,'(5(1x,f10.4))')ff,pl(i,1),pl(i,2),ps(i,1),ps(i,2)
      enddo      
	 
      nwet=0

	do lww=1,lrp
      lump = lazc(lq_b2D(lww))
      do k=1,lump
      nwet=nwet+1
	if(lww.ge.10)then
	if(k.le.5)then
	bio_rum(nwet,3)=ps(lww,1)
	bio_rum(nwet,4)=pl(lww,1)
	else
	bio_rum(nwet,3)=ps(lww,2)
	bio_rum(nwet,4)=pl(lww,2)
      endif
	if(k.gt.12)then
	bio_rum(nwet,3)=0.0001
	bio_rum(nwet,4)=0.0001
      endif
      endif
      enddo
	enddo	
	if(month.eq.1.or.month.eq.12)then
	do nwet=1,lvrp
      bio_rum(nwet,3)=min(bio_rum(nwet,3),1.6)
      bio_rum(nwet,4)=min(bio_rum(nwet,4),1.6)
      enddo
	endif
	write(*,*)'read bio boundary were read'
      endif
      
2     continue	
c	do nwet=1,lvrp
c	write(*,*)'nio_bc 10', bio_rum(nwet,10),nwet
cc	pause
c       enddo
      return
      end

cUD-------------------------------------------------
      subroutine change_units_bd(f,n_choice)
      include 'C_model.f'
      include 'C_nbio_fabm.f'
       parameter (ips=3,ipl=4,ibg=15,izs=5,izl=6,idet=7,inh4=8,
     &idom=9,ino3=10,ipo4=11,isio=12,io2=13,iop=14,ifi1=16,
     &ised1=1,ised2=2,ised3=3,imb1=4)
       real REDF(10)
      data
     &      REDF(1)/6.625/,  !C_N    PrmREDF(1) mol_C/mol_N=    6.625
     &      REDF(2)/106.0/,  !C_P    PrmREDF(2) mol_C/mol_P=  106.000
     &      REDF(3)/6.625/,  !C_SiO  PrmREDF(3) mol_C/mol_Si=    6.625
     &      REDF(4)/16.0/,  !N_P    PrmREDF(4) mol_N/mol_P=   16.000
     &      REDF(5)/1.0/,   !N_SiO  PrmREDF(5) mol_N/mol_Si=    1.000
     &      REDF(6)/12.010/,   !C_Cmg  PrmREDF(6) mg_C/mmol_C=   12.010
     &      REDF(7)/44.6608009/,   !mm_ml  PrmREDF(7)(mmolO2/m**3)/(mlO2/l)   = 44.6608009
     &      REDF(8)/14.007/,   !  mmolN/mgN   
     &      REDF(9)/30.97/,   ! mmolP/mgP
     &      REDF(10)/28.09/   ! mmolSi/mgSii
      dimension f(lvrp,nbio)

      if(n_choice.eq.1)then

*------change units
         do l=1,lvrp

c         do ibio=3,7  !Ps, Pl, Zs, Zl, De

c         f(l,ibio)=f(l,ibio)        !*REDF(6)  !input in   mgC /m**3 to [mmolC/m**3]

c         enddo

        do ibio=8,10 ! NH4, NO2, NO3
         f(l,ibio)=f(l,ibio)*REDF(1)*REDF(6)  !input in mmolN /m**3 to [mmolC/m**3]
         enddo
         f(l,ipo4)  =f(l,ipo4)  *REDF(2)*REDF(6)  !input in mmolP /m**3 to [mmolC/m**3]
         f(l,isio)  =f(l,isio)  *REDF(3)*REDF(6)  !input in mmolSi/m**3 to [mmolC/m**3]
         f(l,iop)  =f(l,iop)  *REDF(3)*REDF(6)  !input in mmolSi/m**3 to [mmolC/m**3]
         f(l,io2)  =f(l,io2)  *REDF(7)          !input in mmolO2/m**3 to milliliters O2/liter 
         enddo !l=1,nrei

      else
*------change units
      endif
      return
      end
cuD---------------------------------------------------------------------------------------


*��8 bioneu
c-----------------------------------------------------------------------
c      subroutine bioneu (iindex,izet,time,ivar,ibio)
      subroutine bioneu (iindex,izet,time,irum)
c-----------------------------------------------------------------------
      include 'C_model.f'
      include 'C_nbio_fabm.f'
      include 'C_mpi.f'
c        include 'C_rbio.f'
		  parameter(ilop1=ilo+1) 
		  parameter(khor1=khor+1) 

        dimension f(m,n,ilo)
        dimension iindex(m,n),izet(m,n)
      common /gitter/ dt,r,g,dl,dlr,dlrh,dln(m),rdln(m),rad,dth,
     *  dlvo(m),dlvu(m),gh,rdt,dt2,r4,coru(m),corv(m),dtrdln(m),dtdlr
      common /ind/ iwet(khor1),ldep(khor),lazc(khor),indend(n),
     *             isornr(n),isorsr(n),islab(n)
      common /liq_b/lq_b(iranz,5),lq_b2D(iranz*m),lq_b3D(iranz*m*ilo),
     &  num_b(7),nul_p(6)
      common /T_4d/ Tc(ndrei,nbio),dTc(ndrei,nbio),
     &Tc_mit(ndrei,nbio+nsed)
c      common /tscdb/ sarum(lvrp),tarum(lvrp),bio_rum(lvrp,nbio)
      common /biobd/ bio_rum(lvrp,nbio)
      integer ivar,ibio,irum(nbio)

**** conditions for NO2, NO3, P, Si, O2 from WAO 01**********
c      if(ibio>0)then
       do ivar=3,nbio
       if(irum(ivar).eq.1)then
       do l=1,lvrp
       Tc (lq_b3D(l),ivar)=bio_rum(l,ivar)
       enddo
       endif
! Zs,Zl,D,NH4 are not set, no external condition used


c     Westrand
      do i=156,164
        l   = iindex(i,16)
        ll  = iindex(i,17)
        lll = iindex(i,18)
c        if(ibio.ge.10)lll = iindex(i,17)
c        if(ibio.eq.13)lll = iindex(i,18)
        if(irum(ivar).eq.1)lll = iindex(i,17)

        lw  = izet(i,16)
                   lump=0
        if(lw.gt.0)lump = lazc(lw)
        do k=1,lump
          Tc(l +k,ivar)  = Tc(lll+k,ivar)
          Tc(ll+k,ivar)  = Tc(lll+k,ivar)
        enddo
      enddo


c     Nordrand 
      do j=12,57
	  l   = iindex(66,j)
        ll  = iindex(67,j)
        lll = iindex(68,j)
c        if(ibio.ge.10)lll = iindex(67,j)
c        if(ibio.eq.13)lll = iindex(68,j)
        if(irum(ivar).eq.1)lll = iindex(67,j)
        lw  = izet(66,j)
	             lump=0
        if(lw.gt.0)lump = lazc(lw)
        do k=1,lump
          Tc(l +k,ivar)  = Tc(lll+k,ivar)
          Tc(ll+k,ivar)  = Tc(lll+k,ivar)
       enddo
      enddo
       enddo

       end 
#endif

*��9 uvrand
c-----------------------------------------------------------------------
      subroutine uvrand (izet,iindex,time)
c-----------------------------------------------------------------------
c     boundary values for u, v, w and Av  
c-----------------------------------------------------------------------

      include 'C_model.f'
      parameter(ilop1=ilo+1) 
      parameter(khor1=khor+1) 
 
      dimension izet(m,n)
      integer iindex(m,n)
c      integer*2 isor,ksor
      dimension f(m,n,ilo),f1(m,n,ilo)
      
      common /ind/ iwet(khor1),ldep(khor),lazc(khor), 
     * indend(n),isornr(n),isorsr(n),islab(n) 
      
      common uc(ndrei),vc(ndrei),stc(ndrei),avc(ndrei),z(m,n)
      common zac(khor),wobc(khor),stuvc(khor),fricu(khor),cxc(khor)
      common cyc(khor),pac(khor),txc(khor),tyc(khor)
      common stpc(ndrei),sac(ndrei),tec(ndrei)
      common pres(ilo),wc(ndrei),fricv(khor)

*=======================================================
c Attention!!!!!! Specify your open boundary lines------
*-1-----Northern boundary 
*     u(1)=u(2)=u(3)
*     v(1)=v(2)
*     w(1)=w(2)=w(3)
*     avc(1)=avc(2)=avc(3)
*-2-----Western boundary
*     u(1)=u(2)
*     v(1)=v(2)=v(3)
*     w(1)=w(2)=w(3)
*     avc(1)=avc(2)=avc(3)
*-3,4---S and E boundaries
*     u(1)=u(2)
*     v(1)=v(2)
*     w(1)=w(2)
*     avc(1)=avc(2)
*=======================================================

cNorthern boundary
      do j=12,57
	  l   = iindex(66,j)  !First line
        ll  = iindex(67,j)  !Second line
        lll = iindex(68,j)  !Third line

        lw  = izet(66,j)
	             lump=0
        if(lw.gt.0)lump = lazc(lw)
        do k=1,lump
          uc (l +k)  = uc (lll+k)
          uc (ll+k)  = uc (lll+k)
          vc (l +k)  = vc (ll +k)

          wc (l +k)  = wc (lll+k)
          wc (ll+k)  = wc (lll+k)
          avc(l +k)  = avc(lll+k)
          avc(ll+k)  = avc(lll+k)
        enddo
      enddo


cWestern boundary
      do i=156,164
        l   = iindex(i,16) !First line
        ll  = iindex(i,17) !Second line
        lll = iindex(i,18) !Third line

        lw  = izet(i,16)
	             lump=0
        if(lw.gt.0)lump = lazc(lw)
        do k=1,lump
          vc (l +k)  = vc (lll+k)
          vc (ll+k)  = vc (lll+k)
          uc (l +k)  = uc (ll +k)
 
          wc (l +k)  = wc (lll+k)
          wc (ll+k)  = wc (lll+k)
          avc(l +k)  = avc(lll+k)
          avc(ll+k)  = avc(lll+k)
        enddo
      enddo

	    goto 5    !output
		 
		open(43,file='test_uvrand_u')
		open(44,file='test_uvrand_v')

	   lwe = 0
         nwet=0
	    do j=1,n
        lwa = lwe+1
        lwe = indend(j)
        do lw=lwa,lwe
        i = iwet(lw)
               lump = lazc(lw)
        do k=1,lump
        nwet = nwet+1
         f1(i,j,k) = uc (nwet)
         f (i,j,k) = vc (nwet)
		enddo
        enddo
        enddo

		do k=1,ilo

!u
      write(43,*)'k=',k,'////, timestep=',time/(60*2)+1
      write(43,'(3i3, 9(1x,i6))')
     &	    1,k, i,( (j), j=4,n-2)

		do i=4,m-2
		write(43,'(3i3, 9(1x,f6.2))')
     &	    1, k, i,( f1(i,j,k), j=4,n-2)
  	if(mod(i-1,11).eq.0)write(43,'(3i3, 9(1x,i6))')
     &	    1,k, i,( (j), j=4,n-2)
        enddo


!v
      write(44,*)'k=',k,'////, timestep=',time/(60*2)+1
	write(44,'(3i3, 9(1x,i6))')
     &	    1,k, i,( (j), j=4,n-2)

		do i=4,m-2
		write(44,'(3i3, 9(1x,f6.2))')
     &	    1, k, i,( f(i,j,k), j=4,n-2)
	if(mod(i-1,11).eq.0)write(44,'(3i3, 9(1x,i6))')
     &	    1,k, i,( (j), j=4,n-2)
          enddo


		enddo !k
		close(43)
		close(44)
		write(*,*)'uvrand otput'
c		pause

5      return 
      end 



*��10 bforce
      subroutine bforce (f)
c-------------------------------------------------------------------
c   boundary mirror points for atmospheric forcing
c-------------------------------------------------------------------
      include 'C_model.f'
      dimension f(m,n)

      do j=1,57
          f(66,j)=f(68,j)
          f(67,j)=f(68,j)
        enddo
c
       do i=155,165
            f(i,16)=f(i,18)
            f(i,17)=f(i,18)
          enddo


       return
       end


*��11 zetadb
      subroutine zetadb(nday,iminu  ,iindex,nzac,time)
c-----------------------------------------------------------------------
c    boundary values for sea surface elevation 
c-----------------------------------------------------------------------
      include 'C_model.f'
      parameter(ilop1=ilo+1) 
      parameter(khor1=khor+1) 

      dimension iindex(m,n),zmis(m,n)
      common /urand/       zday(lrp),zvar(lrp),zalt(lrp)
      common /uranddb/ zday_old(lrp),zvar_old(lrp)
       
      common uc(ndrei),vc(ndrei),stc(ndrei),avc(ndrei),z(m,n)
      common zac(khor),wobc(khor),stuvc(khor),fricu(khor),cxc(khor)
      common cyc(khor),pac(khor),txc(khor),tyc(khor)
      common stpc(ndrei),sac(ndrei),tec(ndrei)
      common pres(ilo),wc(ndrei),fricv(khor)

 
      common /liq_b/lq_b(iranz,5),lq_b2D(iranz*m),lq_b3D(iranz*m*ilo),
     &  num_b(7),nul_p(6)
 
c----Attention!!!, must be changed with respect to time step
       if(mod(iminu,20).eq.0) then  
       call readzdv(zvar,lrp,nzac) !z_v every 20min, tidal SSE
      endif 
          nn2D=0
          DO num=1,iranz
c       print*,'bound',num,lq_b(num,:)
          do inum=lq_b(num,1),lq_b(num,2)
             do jnum=lq_b(num,3),lq_b(num,4)
              nn2D=nn2D+1

           zmis(inum,jnum) = zday(nn2D) +zvar(nn2D)-zalt(nn2D)   
           zalt(nn2D)      = zday(nn2D) +zvar(nn2D)

           enddo
           enddo
             ENDDO	!	 do num=1,iranz

*---------north_b----------------
           do i=156,164  
           do j=16,17   
          z(i,j)=zmis(i,17) 
            enddo
            enddo

            do j=12,57   
           do i=66,67   
           z(i,j)=zmis(67,j) 
           enddo
            enddo 
*---------north_b----------------


         goto 5 
         k=1

         open(43,file='test_zetadb')
         write(43,*)nday,time,'////',time/(60*2)+1,'iminu=',iminu
         write(43,'(3i3, 9(1x,i9))')
     &    1,k, i,( (j), j=4,n-2)
        do i=4,m-2
        write(43,'(3i3, 9(1x,f9.5))')
     &     1, k, i,( z(i,j), j=4,n-2)
        if(mod(i-1,11).eq.0)write(43,'(3i3, 9(1x,i9))')
     &    1,k, i,( (j), j=4,n-2)
        enddo
        close(43)
        write(*,*)'zetadb otput'
c	pause

 5      return 
        end


*��12 readtsin
c-----------------------------------------------------------------------
      subroutine readtsini(nsanf,ntanf,sac,tec,n3d)
c-----------------------------------------------------------------------
c     read temperature and salinity initial conditions  
c-----------------------------------------------------------------------
      include 'C_model.f'
     
      dimension t(m,n,ilo),s(m,n,ilo)
      real sac(ndrei),tec(ndrei)

      read(nsanf,1000) s
      read(ntanf,1000) t
 1000 format(20f8.3)
#ifdef MPI
      call comp3d1d_s (s,sac)
      call comp3d1d_s (t,tec)
#else      
      call comp3d1d (s,sac)
      call comp3d1d (t,tec)
#endif
   
      return
      end



*��13 topocorr
c-----------------------------------------------------------------------
      subroutine topocorr(izet,iizet,ht,m,n,imax)
c-----------------------------------------------------------------------
c     topography corrections
c-----------------------------------------------------------------------
      include 'C_mpi.f' 
      integer imax
      dimension izet(m,n),iizet(m,n),ht(m,n)

*===========================================================
c------ Attention!!!!! ----------------------------
c     set corrections of your topography
c-------------------------------------------------

c-*1-- set minimal thickness of water column --------------
      do i=1,m
      do k=1,n
      if(izet(i,k).gt.0.and.izet(i,k).lt.4) izet(i,k)=4 ! [m]
      enddo
      enddo
      
      izet(97,10)=izet(97,10)+1
      izet(127,28)=izet(127,28)+1

	

c-*2-- correct topography -------------

      izet(168,25) = 7

      do i=156,164
        izet(i,16) = izet(i,18)
        izet(i,17) = izet(i,18)
      enddo

      do i=1,m
      do j=1,n
        if((i.le.67).and.(j.le.58)) izet(i,j) = 0
        if((i.le.68).and.(j.le.10)) izet(i,j) = 0
        if((i.le.69).and.(j.le.10)) izet(i,j) = 0
        if((i.le.70).and.(j.le.10)) izet(i,j) = 0
        if((i.le.71).and.(j.le.10)) izet(i,j) = 0
        if((i.le.72).and.(j.le.8)) izet(i,j)  = 0
        if((i.le.73).and.(j.le.8)) izet(i,j)  = 0
        if((i.le.74).and.(j.le.8)) izet(i,j)  = 0
      enddo
      enddo

      do i=66,67
      do k=11,57
        izet(i,k) = izet(68,k)
      enddo
      enddo

      izet(116,94) = 20
      izet(112,91) = 20
      izet(111,92) = 20
      izet(102,92) = 20
      izet(101,92) = 20


c     small belt

      izet(103,91) = 30
      izet(104,91) = 30
      izet(105,91) = 30
      izet(106,91) = 30
      izet(107,91) = 30
      izet(108,91) = 30
      izet(109,91) = 30
      izet(111,91) = 30
      izet(112,91) = 30
      izet(113,91) = 30
      izet(114,91) = 30
      izet(115,91) = 30
      izet(99,92)  = 30
      izet(100,92) = 30
      izet(101,92) = 30
      izet(102,92) = 30
      izet(103,92) = 30
      izet(106,92) = 30
      izet(107,92) = 30
      izet(108,92) = 30
      izet(109,92) = 30
      izet(110,92) = 30
      izet(111,92) = 30
      izet(115,92) = 30
      izet(97,93)  = 30
      izet(98,93)  = 30
      izet(99,93)  = 30
      izet(100,93) = 30
      izet(101,93) = 30
      izet(108,93) = 30
      izet(109,93) = 30
      izet(115,93) = 30
      izet(116,93) = 30
      izet(97,94)  = 30
      izet(115,94) = 30
      izet(116,94) = 30
      izet(95,95)  = 30
      izet(96,95)  = 30
      izet(97,95)  = 30
      izet(98,95)  = 30
      izet(116,95) = 30
      izet(117,95) = 30
      izet(95,96)  = 30
      izet(96,96)  = 30
      izet(97,96)  = 30
      izet(98,96)  = 30
      izet(116,96) = 30
      izet(117,96) = 30
      izet(118,96) = 30
      izet(96,97)  = 30
      izet(97,97)  = 30
      izet(117,97) = 30
      izet(118,97) = 30
      izet(117,98) = 30
      izet(118,98) = 30
      izet(116,99) = 30
      izet(117,99) = 30
c!!!
      izet(116,100) = 30
      izet(108,102) = 30
      izet(109,102) = 30
      izet(110,102) = 30
      izet(111,102) = 30
      izet(112,102) = 30
c!!!
      izet(109,103) = 30
      izet(110,103) = 30
      izet(111,103) = 30
      izet(109,104) = 30
      izet(110,104) = 30
      izet(111,104) = 30
      izet(110,105) = 30
      izet(111,105) = 30


      izet(86,147) = 220
      izet(86,148) = 220
      izet(87,146) = 218
      izet(87,147) = 250
      izet(87,148) = 220
      izet(88,146) = 230
      izet(88,147) = 235
      izet(88,148) = 210
      izet(89,146) = 210
      izet(89,147) = 232
      izet(90,146) = 212
      izet(90,147) = 215
      izet(91,146) = 217
      izet(91,147) = 196

c     Stolpe Rinne Danziger Tief
      izet(108,125) = 56
      izet(107,125) = 57
      izet(109,125) = 56
      izet(108,126) = 75
      izet(108,127) = 76
      izet(108,128) = 75
      izet(108,129) = 77
      izet(108,130) = 85
      izet(108,131) = 92
      izet(108,132) = 95
      izet(108,133) = 95
      izet(108,134) = 82
      izet(108,135) = 85
      izet(108,136) = 89
      izet(108,137) = 89
      izet(59,141)  = 202
      izet(60,141)  = 210
      izet(60,142)  = 208
      izet(59,142)  = 190
      izet(60,143)  = 177
      izet(61,143)  = 187
      izet(105,117) = 62
      izet(106,117) = 64
      izet(107,117) = 63
      izet(108,117) = 60
      izet(105,118) = 64
      izet(106,118) = 65
      izet(107,118) = 65
      izet(108,118) = 61

c     oeresund
      izet(99,102)  = 24
      izet(103,102) = 12
      izet(104,103) = 18

c     Skagerrak
      izet(84,78)=248

c----------------------------------------------------------
c End of set corrections of your topography
*===========================================================    
      imax=0
      do k=1,n
      do i=1,m
        if(izet(i,k).gt.imax)imax = izet(i,k)
        iizet(i,k) = izet(i,k)
        ht(i,k) = real(izet(i,k))
      enddo
      enddo
      IF(myid .EQ. 0)THEN 
      write(7,1000)imax
 1000 format(' maximum depth = ',i5)
      endif
      return
      end


*===============END OF  suboroutines to be changed ================================
