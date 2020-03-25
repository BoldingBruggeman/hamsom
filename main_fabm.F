c-----------------------------------------------------------------------
c
c     ECOSMO main program
c     North Sea/ Baltic Sea version  last changes: Oct.2007
c     North Sea/ Baltic Sea version  last changes: Oct.2018
c     model code describing for the physical part of ECOSMO     
c     ecosmo e2e is removed 20.03.2020
c references (Schrum et al. 2006; Daewel and Schrum 2013)
c-----------------------------------------------------------------------

c----------------------------------------------------------------------
c     Include files:
c		model dimension
c		index
c		input-output path information.f
c----------------------------------------------------------------------
#ifdef COASTDAT
      use ncio
#endif
#ifdef _FABM_
      use hamsom_fabm
#endif
      include 'C_model.f'
      include 'C_index.f'
      include 'C_paths.f'
      parameter (nx=m,ny=n,ilop1=ilo+1,khor1=khor+1)
      parameter (nlar=m*(ilo*6+3))

c----------------------------------------------------------------------
c     Parameter settings
c	  numeintr: number of runoff sources considered, either particular river or
c		 runoff in a larger region
c	  nanz:	 number of runoff source points on the model grid
c	  numeintr: number of river runoff values considered
c		    either single river sources or integrated runoff areas
c	  ianzpkt
c	  nceplam: array size for atmospheric input, longitude index
c	  ncephi:  array size for atmospheric input, latitude index
c	  nera:	   number of atmospheric input points
c
c
c---------------------------------------------------------------------
      parameter (nanz=222,numeintr=nanz)
      parameter (ianzpkt=9)
#ifdef NCEP
      parameter (nceplam =30,ncephi  =12,nera =nceplam *ncephi)
      parameter (ncep2lam=18,ncep2phi =10,nera2=ncep2lam*ncep2phi)
#endif
#ifdef COASTDAT
      parameter (ncdlam =145,ncdphi  =107,nera =ncdlam *ncdphi)
#endif 
      parameter (iradia=9)

c-----------------------------------------------------------------------
c     Parameter setttings for partcle tracking module
c     parameter(istat=1) !! particles for movie_
c-----------------------------------------------------------------------

      parameter (ndrei1=ndrei+1)!! particles
      parameter(istat=1)

c-----------------------------------------------------------------------
c     include files for
c	    sea ice dynamics
c	    xy
c	    NPZD ecosystem module
c	    output
c-----------------------------------------------------------------------
      include 'C_ice.f'
      include 'C_xy.f'
      include 'C_daily_physics.f'
#ifdef MPI
      include 'C_mpi.f'
      include 'C_parallel.f'
#endif
c----------------------------------------------------------------------
c    control parameter
c---------------------------------------------------------------------------
      logical neust1,backup,mitout,progno,withice,loxstat,loxstat1
      logical river_on
      logical radinp
      logical if_wet
c---------------------------------------------------------------------------
c     control arrays for slab formulation--
c	     --only wet grid points are adrssed and calculated:
c     iindex(m,n):start index
c     izet(m,n): input array for topography field, later used as index array
c     iizet(m,n):integer topography array
c     jjc(m,n): land mask array
c     jjjc(m,n,ilo): dry grid point mask array
c     id3sur(khor): corresponding grid index of a 2-d compressed array in a 3-d compressed array
c     indm, indn, rivlat, rivlon index arrays for river input
c     itagan: control array for number of days per month
c--------------------------------------------------------------------------

      integer iindex(m,n),iizet(m,n),jjc(m,n),jjjc(m,n,ilo)
      integer id3sur(khor)
      integer indm(nanz),indn(nanz),rivnum(nanz) !,nfuriv(nbio)
      real rivlon(nanz),rivlat(nanz),rivinp(nanz)
      
      integer nord,sued,ost,west
      integer*2 isor,ksor
      integer itagan(12)
c-----------------------------------------------------------------------------
c     control arrays for output statistics
c------------------------------------------------------------------------------
      integer idudo2(m,n)
      dimension dudo3(m,n,ilo)
      dimension dudo4(m,n,ilo),dudo5(m,n,ilo),dudo6(m,n,ilo)
      dimension dudo2(m,n)

      integer	nrmess(m,n,ianzpkt),nrmessp(m,n,ianzpkt)
      dimension weight(m,n,ianzpkt),weightp(m,n,ianzpkt)
      dimension bwind(nera)
#ifdef NCEP
c-----------------------------------------------------------------------------
c     arrays needed for interpolation of atmospheric forcing data
c----------------------------------------------------------------------------
      dimension rlamu(nera),rphiu(nera),rlamup(nera2),rphiup(nera2)

      common /era/ erau10 (nera),erav10 (nera),erap  (nera2),
     *		   eratl  (nera),erash	(nera),eracl (nera),
     *		   erapri (nera),
     *		   eraussr(nera),eradssr(nera),
     *		   eraulwr(nera),eradlwr(nera)
c------------------------------------------------------------------------------
c     common block era:
c		   input arrays for atmospheric forcing
c     intermediate arrays used for atmospheric input and radiation
c     atmospheric forcing arrays on the model grid
c     arrays for climatic temperature corrections
c------------------------------------------------------------------------------

      dimension erassr(nera),erastr(nera)
      character*12 fileu10,filev10,filep,filetl,filesh,filecl,filepri
      character*12 fileussr,filedssr,fileulwr,filedlwr
#endif
#ifdef COASTDAT
      integer metind1(nxra,4)
      real cdlam(nera),cdphi(nera)
       real dump(ncdlam,ncdphi)
      real erau10 (nera),erav10 (nera),erap  (nera),
     *		   eratl  (nera),erash	(nera),eracl (nera),
     *		   erapri (nera)
       real eratp(nera),eraqv(nera),erassr(nera)
       character metappend*24
       character metpost*7
      character(len=20) file_outx
      character(len=11) fileuv10,fileddsr
      character(len=7) fileu10,filev10
      character(len=6) filep,filetl,filecl
      character(len=10) filepri
      character(len=8) fileussr,filedssr,fileulwr,filedlwr
      character(len=7) filesh,filetd
#endif
      character(len=10) carg
      dimension pat(m,n),bwindn(m,n),txx(m,n),tyy(m,n),diru(m,n)
      dimension taair(khor),shum(khor),bew(khor),precip(khor)
      dimension dswr1(khor)
      dimension einstr(m,n),einstrice(m,n),dumm(m,n),dumm1(m,n)
      dimension paac2(khor),txac2(khor),tyac2(khor)
      dimension paac(khor),txac(khor),tyac(khor)
      dimension wgesch(m,n),qwt(khor),evap(khor)
      dimension qw(khor),taup(khor),wspeed(khor)
      dimension q10(khor),q10a(khor),tl10(khor),tl10a(khor)
      dimension cddu(khor)

      dimension tempkorr(m,n),ssrkorr(m,n),strkorr(m,n)
!      dimension testt(10,10)
c-----------------------------------------------------------------------------
c     mradi:, nradi:, iradi
c     szahl:
c     teca: climatic temperature distribution
c     saca: climatic salinity distribution
c     z2d:
c     qoi:
c     qii:
c     exbio:
c-----------------------------------------------------------------------------
      dimension mradi(iradia),nradi(iradia),iradi(iradia)

      dimension szahl(ndrei),teca(ndrei),saca(ndrei)
      dimension zmalt(khor)
      dimension z2d(m,n)

      dimension qoi(m,n),qii(m,n)
c      dimension exbio(ndrei)

c------------------------------------------------------------------------------
c     grid specifications and exchange
c------------------------------------------------------------------------------

      dimension ht(m,n),dzbod(m,n)
      dimension dz(ilo),pav(ilo),pah(ilo),pah_fac(ndrei)
      dimension zdif(m,n),zic(khor)
      dimension dx(m,n),dy(m,n)
#ifndef MPI
      common /geo1/ ltief(m,n),izet(m,n)
#endif
      dimension AAA(0:khor)
      common /zeit/ ijahr
c--------------------------------------------------------------------------------
c     ukonv, vkonv: help arrays for non-linear terms
c     xpro
c     eintr:
c     ftac:
c     zpe:
c     zriv:
c---------------------------------------------------------------------------------
      dimension ukonv(ndrei),vkonv(ndrei)
      dimension xpro(nanz),eintr(nanz)
      dimension ftac(khor)
      dimension zpe(m,n),zriv(m,n)
c------------------------------------------------------------------------------------
c    u22t, v22t, w22t checking for time step stability critera
c------------------------------------------------------------------------------------
      dimension u22t(ndrei),v22t(ndrei) ,w22t(ndrei)


      dimension txac_2(khor),tyac_2(khor),wspeed_2(khor)
      dimension qwt_2(khor),q10_2(khor),tl10_2(khor)
      dimension fqs_2(khor),fql_2(khor),evap_2(khor)

c------Ute test
      real  n2init(m,n,ilo),delsal(6)
      integer spinup,noyears
c      character*11 filewNH4,filewNO3  !for wet deposition scenarios
c      dimension wetNH4(m,n),wetNO3(m,n)
c      dimension wetNH41d(khor),wetNO31d(khor),dumm1d(khor)
#ifdef WADDEN
       real bound_wadden(43,2)
       common /wadden/ bound_wadden
#endif 
*----tecplot------------
      dimension f(m,n,ilo,nbio)

*----tecplot------------

      dimension pd2(m,n,0:ilo+2),dzt(m,n,ilo)

c--------------------------------------------------------------------------------------
c     common block aux and auxint; internal common blocks, help arrays
c	      common block lengths aux = m*(ilo*19+10),auxint = m*(ilo*6+3)
c	      attention !!!!!!!
c--------------------------------------------------------------------------------------

      common /aux/ gross(ngro)
      common /auxint/ large(nlar)
c-------------------------------------------------------------------------------------
c     common blocks
c		init:counter, help and iteration
c		iterat: iteration, surface elevation
c		uirand: boundary arrarys
c		dreh: corriolis rotation term
c		gitter2: grid
c		num: help slab arrays
c		radi:radiation
c		radis:radiation in case of sea ice
c		wetpos:index arrays
c		urand, uranddb: new and old boundary sea level
c		tsc, tscdb: external temperature and salinity open boundary values
c		fast: ice model
c		dinv
c		atm
c		liq_b: open boundary control
c-------------------------------------------------------------------------------------
      common /intit/ itel(1440),jc(m,3),ikk(kasor),isor(kasor),
     *	     ksor(kasor),nofrez,lone,itende,nfreez
      common /iterat/ omega(kasor),za(m,3),cx(m,3),cy(m,3),rhoz(m,3),
     *	     tcc(8)
      common /uirand/ilwrnd(lrp),islrnd(lrp),ldownr(lrp),
     *	     irangr(iranz),iea,iee,iwa,iwe,
     *	     nord,sued,ost,west,kna,kne,ksa,kse
      common /dreh/ sinfu(m),sinfv(m),cosfu(m),cosfv(m),sincx(m),
     *	     sincy(m),bx(m),by(m),pxu(m),pyv(m),pyu(m),pxv(m)
      common /gitter2/ fianf,yambdanf,dphi,dlambda,rearth
      common /num/ dc(ilo),av(ilo),ad(ilo),dh(ilo),pd(ilo),
     *	     prd(ilo),pr2d(ilo),r2d(ilo),tkw(ilo),tau(ilo),
     *	     dd(ilo),qa(ilop1),qbet(ilop1),qn(ilop1),rd(ilop1)

      common /radi/fqg(m,n),fqr(khor),fqs(khor),fql(khor)
      common /radis/fqgmit(m,n),fqrmit(khor),fqsmit(khor),fqlmit(khor)

c      common /wetpos/ indx(nstat),indy(nstat)

      common /urand/	   zday(lrp),zvar(lrp),zalt(lrp)
      common /uranddb/ zday_old(lrp),zvar_old(lrp)
c--- bio_rum bio boundary
c      common /tscdb/ sarum(lvrp),tarum(lvrp),bio_rum(lvrp,nbio)
      common /tscdb/ sarum(lvrp),tarum(lvrp)
      common /tsc/ sarum_old(lvrp),tarum_old(lvrp)

      common /fast/ muuv(nx,ny)
      common /dinv/ s1(m),s2,s3(m)
      common /atm/ rnil,istop

      common /liq_b/lq_b(iranz,5),lq_b2D(iranz*m),lq_b3D(iranz*m*ilo),
     &	num_b(7),nul_p(6)

c------------------------------------------------------------------------------
c     common blocks
c		T_4d:variables NPZD module+ temperature and salinity
c		bio_riv: nutrient loads from rivers
c------------------------------------------------------------------------------

      common /T_4d/ Tc(ndrei,nbio) !Tc_mit(ndrei,nbio+nsed)
c      common /bio_riv/riv_bio(nbio,nanz),
c     &		     zriv_bio(nbio,m,n)


      double precision zinc(m,n),zinc1(m,n)
      double precision wsur(khor)
      dimension zaczac(khor)
c--------------------------------------------------------------------
      dimension Tc_flux_r(khor,nbio),Tc_flux_e(khor,nbio), !river and P-E fulxes
     & Tc_flux_all(khor,nbio)


c----------------------------------------------------------------------------------
c     common locus: particle tracking module
c----------------------------------------------------------------------------------

      common /locus/ xfort(10,istat),yfort(10,istat)
     &,zfort(10,istat),lac(ndrei1),ito(m,n,ilo)
c--------------------------------------------------------------------------------
c     Input output character variables
c--------------------------------------------------------------------------------
#ifdef MPI
      character outfile*10, startfile*10,konv11*6,outfile1*10
     * ,outfile2*10
#else
      character outfile*8, startfile*10
#endif

      character*12 metein
      character*10 metloesch
      character*6 mrunnext
      character*5 ausnext

      character*7  topofile
      character*16 rivinput
      character*12 rivpoint

c--------------------------------------------------------------------------------
c     Input character variables atmospheric forcing
c--------------------------------------------------------------------------------

      character*81 cdum
      character*14 filetkor
      character*19 metinde2

      character*1 readt
      integer ihz, iclock0, iclock1
      real time



c------------------------------------------------------------------------------
c ! T_S boudary conditions   flie=t&s_YYYY, nfutsb1 = 65
c	 z_d,z_v:  zet boundary daily means, and zet daily variations
c------------------------------------------------------------------------------
      character*4 ijahr_1
      character*12 filradou
c      character*16 rivinp_bio

      character*12 klimtex1,klimtex2,sal_aa
      character*10 T_S_bd !bio_bd
      character*10 z_d,z_v

      character text1*9,text2*45,text3*22,jtext*2,ltext*3

      character*8 timefig
      character*20 filel
#ifdef MPI
      character gridinf*10,fbackup*6,control*13,ppp*3
#else
      character gridinf*10,fbackup*6,control*10,ppp*3
#endif
      character*40 text_m

      character*7  nmodel
      character*19 pathtsb,pathzeb,pathout,pathtsb2,pathout2
      character*17 pathtop,pathset,pathriv,pathmet
#ifdef MPI
      common /cord1/ icord(nprocs,4)
#endif

#ifdef _FABM_
!KB hardcoding npel to 5 !!!
      call init_hamsom_fabm(m,n,ilo,ndrei,5,sac,tec)
#if 1      
      write(*,*) indend
      stop 'kaj'
#endif      
!KB      stop 'init_hamsom_fabm()'
#endif


#ifdef MPI
      call init_parallel(myid,nprocs,ierr)
      if(nprocs .le. 1) then
      write(*,*) 'increase number of processors else use sequential'
      stop
      endif
#else
      myid=0
      nprocs=0
#endif

c---------------------------------------------------------------------------------
c     set input/output paths
c---------------------------------------------------------------------------------
      nmodel='north_b'	!model name
      pathtsb=nmodel//'/input_YYYY/'	! pathpc//pathtsb
      pathzeb=nmodel//'/input_YYYY/'	! pathpc//pathzeb
      pathout=nmodel//'/f_out_YYYY/'	! pathpc2//pathout
      pathout2=nmodel//'/n_out_YYYY/'	! pathpc2//pathout

      pathtop= 'h_input/topinput/'
      pathset= 'h_input/setinput/'
      pathriv= 'h_input/rivinput/'
#ifdef NCEP
      pathmet= 'h_input/metinput/'
#endif
#ifdef COASTDAT
      pathmet= 'h_input/metinCD2/'
#endif


c-------------------------------------------------------------------
c     read in from the screen, star year for calculation
c-------------------------------------------------------------------
c      read (*,*) ijahr

c-------------------------------------------------------------------
c      ppp:    characteristic character for the run (input/output names)
c      avmax:	   max. vertical exchange coefficient
c      read_pah:   max. horizontal exchange coefficient
c-----------------------------i------------------------------------------
       ppp='ttt'
c       open(42,file='input_'//ppp)
cc	open(42,file='input_t03')
c       read (42,*) ijahr_start
c       read (42,*) noyears
      call getarg(1,carg)
      read(carg,*)ijahr_start 
      noyears=1 
c      close(42)
      ioutput=0
      moutput=0
      avmax=1000.e-4  !from 1990

       read_pah=50.  !do 1993 pah=90, after 50
c       read_pah=90.  !do 1993 pah=90, after 50
c      read_pah=60.  !do 1993 pah=90, after 50
c      read_pah=45.  !do 1993 pah=90, after 50

	    avmax1=avmax
      IF(myid .EQ. 0)THEN
      print*,'jahr_start',ijahr_start
      write(*,'(a7,e16.5)')'avmax',avmax
      write(*,*) 'ppp= ', ppp
      write(*,'(a7,f 6.1)')'pah',read_pah
      ENDIF

c-----------------------------------------------------------------------
c	gridinf: grid info file
c	control: control output
c	outfile: monthly file containing daily averaged output
c	 YYMMDDHH-->
c	       YY: year
c	       MM: month
c	       DD: day
c	       HH: hour
c	filradou:
c-----------------------------------------------------------------------

      gridinf  = 'pppgrdinfo'
      control  = 'pppcontrol'
      outfile  = 'pppYYMMb'
      filradou = 'pppradou.dat'

c-----------------------------------------------------------------------
c     set input file names
c	 topofile: bathymtry
c	 rivpoint: river input coordinates
c	 metinde2: atmospheric input points
c	 atmospheric input files, separate parameter
c	      u10. 10 m wind E-W direction
c	      v10. 10 m wind N-S direction
c	      p	   sea surface pressure
c	      tl   2m air temperature
c	      sh   2m specific humidity
c	      cl   cloudiness
c	      pri  precipitation rate
c	      ussr short wave radiation upward
c	      dssr short wave radiation downward
c	      ulwr  long wave radition upward
c	      dlwr  long wave radiation downward
c
c    find information about input units and corresponding fields behind
c-----------------------------------------------------------------------

      topofile	= 'modtopo'

      rivpoint	= 'rivpoint.dat'
#ifdef NCEP
      metinde2	= 'positions_ncep.nsbs'
      fileu10	= 'UUGRD_YYYYMM'  ! 1 10m wind x   67	era u10	   txx(i,j)
      filev10	= 'VVGRD_YYYYMM'  ! 2   68	era v10	   tyy (i,j)
      filep	= 'SSSLP_YYYYMM'  ! 3   69	era p	   pat (i,j)
      filetl	= 'TTTMP_YYYYMM'  ! 4 2m air temp  70	era tl	   taair(l)
      filesh	= 'SPFHD_YYYYMM'  ! 5 spec hum	   71   taup(l)
      filecl	= 'TCTOT_YYYYMM'  ! 6   72	era cl	   bew(l)
      filepri	= 'PRATE_YYYYMM'  ! 7   73

      fileussr	= 'USWRF_YYYYMM'  ! 8 radiation	   74
      filedssr	= 'DSWRF_YYYYMM'  ! 9	   75
      fileulwr	= 'ULWRF_YYYYMM'  ! 10	   36
      filedlwr	= 'DLWRF_YYYYMM'  ! 11	   37
#endif  
#ifdef COASTDAT 
      metappend ='1948_01.ECOSMO_domain.nc'
      metpost   ='cDII.00'  
      fileuv10  = '.UVlat_10M.'  ! 1 10m wind x   67 era u10
      fileu10   = '.U_10M.'  ! 1 10m wind x   67 era u10
      filev10   = '.V_10M.'  ! 2                    68   era v10
      filep     = '.PMSL.'  ! 3 ai oressure at sea lavel Pa
      filetl    = '.T_2M.'  ! 4 2m air temp  Kelvin 70   era tl
      filesh    = '.QV_2M.'  ! 5 spec hum   kg/kg  71
      filetd    = '.TD_2M.'  ! 5 dew point T Kelvin     71
      filecl    = '.CLCT.'  ! 6             72   era cl
      filepri   = '.TOT_PREC.'  ! precipitation in kg m-2 h-1
c      fileussr = 'USWRF_YYYYMM'  ! 8 radiation    74
      filedssr  = '.ASOB_S.'  ! 9averaged surface net downward sw
      fileddsr  = '.ASWDIFD_S.'  ! 9averaged surface net downward sw
      fileulwr  = '.ALWU_S.'  ! 10                  36
      filedlwr  = '.ALWD_S.'  ! 11                  37
#endif



c-------------------------------------------------------------------------
c     input file names for external open boundary values;
c	find corresponding input units behind
c-------------------------------------------------------------------------

      T_S_bd	= 't_s_YYYYMM'	!  nfutsb   = 65
c----Ute Daewel 26.4.2010
      sal_aa	= 'sal_aayy.dat'!  anual anomaly for boundary salinty from ICES data
c      bio_bd	= 'bio_1984MM'	!  nfutsb   = 65
      z_d	= 'z_d_YYYYMM'	!  nzactp   = 18
      z_v	= 'z_v_YYYYMM'	!  nzac	    = 19

c-----------------------------------------------------------------------
c	begin model setup
c-----------------------------------------------------------------------
      IF(myid .EQ. 0)THEN
      print*,' '
      print*,'start modell setup'
      print*,' '
      ENDIF
c---------------------------------------------------------------------------------
c     control of time performance
c-------------------------------------------------------------------------------
      call ftrace_region_begin ('rosh')
      call system_clock (count_rate=ihz)
      call system_clock(count=iclock0)
#ifdef WADDEN
      nwa=261
      open(unit=nwa,file= pathpc//'h_input/wadden_ECOSMO.dat',
     &form='formatted')
      do iip=1,43
      read(nwa,'(2(e16.7))')bound_wadden(iip,1),bound_wadden(iip,2)
      enddo
      close(nwa)
      if(myid.eq.0.) print*,'test w1',bound_wadden(:,1)
#endif
c-----------------------------------------------------------------------
c	begin year-loop
c	model currently used for one-year runs
c-----------------------------------------------------------------------
c--------+x years spinup----------------------------------------------
c-------if spinup =0 start from climatologi and WOCA else from start up file
c--------------------------------------------------------------------
      spinup=-1
c      noyears=noyears+spinup
      DO 5000 iyears=1,noyears
      icccount=0
      if (iyears.ge.spinup.and.spinup.gt.0)then
      ijahr=ijahr_start+iyears-spinup
      elseif (iyears.ge.spinup.and.spinup.le.0)then
      ijahr=ijahr_start+iyears-1
      else
      ijahr=ijahr_start
      endif
      imetind=0
c-----------------------------------------------------------------------
c     starttime settings
c     mona-mone	  --> month loop
c     itaga-itage --> day loop
c     mjar counter for years
c     print of beginning and end of model run
c-----------------------------------------------------------------------
      itaga  =	 1
      itanf  =	 1
      mona   =	 1

      itage  =31
      mone   =12

c-----------------------------------------------------------------------
c    control  of specification of  setup file
c    ibup=0 set up from daily mean (as stored in monthly files)
c    ibup=1 set up from start up file==> last time step, last year
c-----------------------------------------------------------------------
      ibup =1
      if(mona.gt.1)ibup =0


      write(ijahr_1(1:4),'(i4)')ijahr-1
      mjar   = ijahr-1900
      if(ijahr.ge.2000) mjar = ijahr-2000
      mjar_1 = mjar-1
      if(ijahr.eq.2000)mjar_1 = 99

      IF(myid .EQ. 0)THEN
      write(6,3406)itanf,mona,mjar
      write(6,3407)itage,mone,mjar

 3406 format(3x,'start :',3i3)
 3407 format(3x,'end   :',3i3)
      print*,' '
      ENDIF
c-----------------------------------------------------------------------
c     input/output control
c     neust1 = .true. ==> cold start
c     backup = .true.  ==> write backup to a start file in the end of month
c     mitout = .true.  ==> write daily mean output in monthly files
c     stat= .true.  ==> write statistics to std
c-----------------------------------------------------------------------

      neust1 = .false.
      backup = .false.
      mitout = .true.
      loxstat=	.false.
      loxstat1= .false.

c-----------------------------------------------------------------------
c     prognostic/diagnostic and sea dynamics
c     progno = .false. ==> no SR strom
c     withice= .false. ==> no icemodel
c     river_on= .true. ==> runoff considered
c     radinp=.true. radiation from file, alternative: radiation calculated from
c	      implemented parameterisations
c-----------------------------------------------------------------------

      progno = .true.
      withice  = .true.
      river_on = .true.
      radinp   = .true.
c----------------------------------------------------------------------
c     dz(ilo) vector contains the respective lower bound
c     of the respective model layer
c-----------------------------------------------------------------------

      data dz/5.,10.,15.,20.,25.,30.,35.,40.,48.,56.,64.,72.,
     * 80.,88.,100.,125.,150.,200.,400.,630./

c----------------------------------------------------------------------
c     timestep organisation
c     iviert =	loops per day
c     iper   = time steps per 12/h
c     nsteps = time steps per day
c----time stepping for output
c     tout = deltat for output
c     r1step = timesteps per tout
c     rstep =1/nstep
c--------------------------------
c     deltat = time step in seconds
c     iet= time step for reading met-input
c     i_ice= time step for the dynamic sea ice model
c-----------------------------------------------------------------------
      iviert =  2
c      min_time_step=60*30	!30m
c      min_time_step=60*20      !20m
      min_time_step=60*10      !cd1 1990 onwards
c      min_time_step=60*5      !10m fur b20
      if(ijahr.eq.2010.and.lmon.eq.2)then 
      min_time_step=60*5      !10m fur b20
      endif
      deltat = real(min_time_step)
        tout=24*60*60          !output alle 60 min

       IF(myid .EQ. 0)THEN
       write(*,*)'time steps in seconds	 =',min_time_step,
     &	     ' s, deltat=', deltat
	ENDIF

	nsteps = 24 * 60  *60 /min_time_step
        iper=nsteps/iviert

	IF(myid .EQ. 0)THEN
       write(*,*)'time steps per   DATE,   iper=', iper
	ENDIF
	r1steps = tout /min_time_step
	toutsteps = 24 * 60  *60 /tout

       rsteps =	1./(real(r1steps))
       IF(myid .EQ. 0)THEN
       write(*,*)'time steps per    DAY, nsteps=', nsteps
       write(*,*)'output deltat, tout,r1steps=', tout,r1steps,toutsteps
       ENDIF
c time step for meteorological input
c       6 hours for NCEP
c       1 hour for COASTDAT
c...............................................
#ifdef COASTDAT
       imet=1
#endif
#ifdef NCEP
       imet=6
#endif
       imetdep=6 * 60 *60/ min_time_step   !  time step for nutrient deposition  
       it_init=0
       ireadmet = imet * 60 *60/ min_time_step   !    
c--------------------------------------------------------       
       i_ice = 1 * 60 *60 / min_time_step   !every 4 hours

#ifdef MPI
       nreadmet=24/6
#endif
c-----------------------------------------------------------------------
c     constants:
c	 g= gravitation
c	 r= bottom friction
c	 botidc= thickness of imagary bottom
c	 cdeav=
c	 avmax=max. vertical exchange coefficient, now since critical this is
c	       read from file
c
c	 slip conditions
c     -1 == noslip , 0 == semislip , 1 == fullslip
c
c	 stress = scalig of wind stress tau/rho
c	 wind stress in Pa requires stress = 0.001  since rho=O(1.)
c-----------------------------------------------------------------------

      g = 9.81
      r	     =	 0.005
      botdic =	 10.
      cdeav  =	 0.05

c      avmax  = 1000.e-4  !red from input file

      slip =1
      stress = 1.e-3


c-----------------------------------------------------------------------
c     s.o.r.   weights for time averaging in SOR
c     tmit from interval 0-1, (0 = explicit, 1 = strongly implicit)
c     tmit is weigthing factor for the time level n+1
c-----------------------------------------------------------------------

      vtmit  = 0.85
      vtmic  = 1.-vtmit
      htmit  = 0.85
      htmitq = htmit*htmit


c-----------------------------------------------------------------------
c     iteration constants
c-----------------------------------------------------------------------

      eps    = 1.e-7

c-----------------------------------------------------------------------
c     organisation for zeta boundary values, routine zetadb
c-----------------------------------------------------------------------

      ilom  = ilo-1

c-----------------------------------------------------------------------
c      begin default settings
c-----------------------------------------------------------------------

      nz = n-1
      mz = m-1

      time   = 0.0
      ijulu  = 0
      ihouu  = 0
      iminu  = 0
      isecu  = 0
      jday=0
      ihour  = 0
      igstep = 0
      igstep1 = 0
      nstart = 0
      imal   = 0
      smdriv = 0.0
      tzei   = 0.0


c-----------------------------------------------------------------------
c     set days per month, number of skipped days and start time estimatn
c     start time parameter time in seconds
c-----------------------------------------------------------------------
c     leap-year check

      iss = 0
      if(mod(ijahr,4).eq.0.and.mod(ijahr,100).ne.0) then
       iss = 1
      endif
      if (mod(ijahr,400).eq.0) then
       iss = 1
      endif
      data (itagan(k),k=3,12) /31,30,31,30,31,31,30,31,30,31/
      itagan(1) = 31
      itagan(2) = 28+iss

      iskipd = 0

      if(itanf.eq.1.and.mona.eq.1)then
	iskipd = 0
      else
	do lstmon=1,mona
	  if(lstmon.eq.mona)then
	    itagvo = itanf-1
	  else
	    itagvo = itagan(lstmon)
	  endif
	  iskipd = iskipd+itagvo
	enddo
      endif

      time = 0.0
      do i=1,iskipd
	time = time+86400.0
      enddo

c-----------------------------------------------------------------------
c      begin   zero initialisation
c-----------------------------------------------------------------------

      pnull = 101.3

      do l=1,nlar
	large(l) = 0
      enddo
      do l=1,ngro
	gross(l) = 0.0
      enddo


      do l=1,iranz
	irangr(l) = 0
      enddo

      do l=1,lrp
	itel(l)	  = 0
	zalt (l)  = 0.0
	zday(l)	 = 0.0
	zvar(l)	 = 0.0
	zday_old(l)  = 0.0
	zvar_old(l)  = 0.0
	ldownr(l) = 0
	ilwrnd(l) = 0
	islrnd(l) = 0
      enddo

      do i=1,lvrp
      tarum(lvrp)=0.0
      sarum(lvrp)=0.0
      tarum_old(lvrp)=0.0
      sarum_old(lvrp)=0.0
      enddo

      do i=1,ndrei
	uc(i)	 = 0.0
	vc(i)	 = 0.0
	ukonv(i) = 0.0
	vkonv(i) = 0.0
	wc(i)	 = 0.0

	u22t(i)	   = 0.0	! for checking dt criteria (after call motmit)
	v22t(i)	   = 0.0	!
	w22t(i)	   = 0.0	!

	umit(i)	 = 0.0
	vmit(i)	 = 0.0
	wcmit(i) = 0.0
	scmit(i) = 0.0
	szmit(i) = 0.0
	acmit(i) = 0.0
	tcmit(i) = 0.0
	sac(i)	 = 0.0
	tec(i)	 = 0.0
	saca(i)	 = 0.0
	teca(i)	 = 0.0
	avc(i)	 = 0.0
	stpc(i)	 = 0.0
	stc(i)	 = 0.0
      enddo

      do j=1,n
	islab(j)  = 0
	isornr(j) = 0
	isorsr(j) = 0
	indend(j) = 0
c	 rlamg(j)  = 0.0
	do i=1,m
c	   rphig(i)	  = 0.0
	  frice(i,j)	 = 0.0
	  his(i,j)	 = 0.0
	  ui(i,j)	 = 0.0
	  vi(i,j)	 = 0.0
	  tis(i,j)	 = 0.0
	  izet(i,j)	 = 0
	  jjc(i,j)	 = 0
	  iindex(i,j)	 = 0
	  frimit(i,j)	 = 0.0
	  hismit(i,j)	 = 0.0
	  hisrmit(i,j)	 = 0.0
	  tismit(i,j)	 = 0.0
	  uimit(i,j)	 = 0.0
	  vimit(i,j)	 = 0.0
	  z(i,j)	 = 0.0
	  tair(i,j)	 = 0.0
	  dumm(i,j)	 = 0.0
	  einstrice(i,j) = 0.0
	  einstr(i,j)	 = 0.0
	  zdif(i,j)	 = 1.2
c-- test ice array
           drags(i,j)	 = 0.0
           draga(i,j)	 = 0.0
           gwatx(i,j)	 = 0.0
           gwaty(i,j)	 = 0.0
           amass(i,j)	 = 0.0
           muv(i,j)	 = 0.0
           mhc(i,j)	 = 0.0
           maphc(i,j)	 = 0.0
           forcex(i,j)	 = 0.0
           forcey(i,j)	 = 0.0
           uerr(i,j)	 = 0.0
           verr(i,j)	 = 0.0
           hisr(i,j)	 = 0.0
c--------------------------           
           z2d(i,j)	 = 0.0
	  txx(i,j)	 = 0.0
	  tyy(i,j)	 = 0.0
	  fqg(i,j)	 = 0.0
	  fqgmit(i,j)	 = 0.0
	  qois(i,j)	 = 0.0
	  qoi(i,j)	 = 0.0
	  qiis(i,j)	 = 0.0
	  qii(i,j)	 = 0.0
	  zriv(i,j)	 = 0.
	  zpe(i,j)	 = 0.
	  zinc(i,j)	 = 0.0d+0
	  zinc1(i,j)	  = 0.0d+0
	  do jj = 1,ilo
	    jjjc(i,j,jj) = 0
	  enddo
	enddo
      enddo

      do i=1,khor
	zic(i)	 = 0.0
	zac(i)	 = 0.0
	zaczac(i)= 0.0
	zmit(i)	 = 0.0
	zmalt(i) = 0.0
	fricu(i) = 0.0
	fqr(i) = 0.0
	fricv(i) = 0.0
	wobc(i)	 = 0.0
	stuvc(i) = 0.0
	pac(i)	 = pnull
	paac2(i) = pnull
	txac(i) = 0.0
	tyac(i) = 0.0
	txac2(i) = 0.0
	tyac2(i) = 0.0
	txc(i)	 = 0.0
	tyc(i)	 = 0.0
	cxc(i)	 = 0.0
	cyc(i)	 = 0.0
	ldep(i)	 = 0
	lazc(i)	 = 0
	iwet(i)	 = 0
	qw(i)	  = 0.0
	qwt(i)	   = 0.0
	bew(i)	  = 0.0
	taup(i)	  = 0.0
        dswr1(i)   = 0.0
	taair(i)  = 0.0
	shum(i)	 = 0.0
	wspeed(i) = 0.0
	id3sur(i) = 0
      enddo


	 txac_2=0.;tyac_2=0.;wspeed_2=0.
	 qwt_2=0.;q10_2=0.;tl10_2=0.;
	 fqs_2=0.;fql_2=0.;evap_2=0.;

#ifdef _FABM_
#ifdef MPI
      call deco1d3d_s (dudo3,sac,ndrei)
#else
      call deco1d3d (dudo3,sac,ndrei)
#endif
#endif

c-------------------------------------------------------------------------------
c   zero-initialisation of biological variables
c-------------------------------------------------------------------------------

	 call Rzero(Tc  ,ndrei*nbio)      !Biologycal 3D fields

c-----------------------------------------------------------------------
c      begin change time dependent file names
c-----------------------------------------------------------------------
      IF(myid .EQ. 0)THEN
      print*,'	change time dependent file names'
      ENDIF
      write(pathtsb(15:18),'(i4.4)') 1948
      write(pathzeb(15:18),'(i4.4)') ijahr
      write(pathout(15:18),'(i4.4)') ijahr
      write(pathout2(15:18),'(i4.4)') ijahr

      IF(myid .EQ. 0)THEN
	write(*,*)'pathout=',pathout
	write(*,*)'//////'
      ENDIF

#ifdef NCEP
      write(fileu10 ( 7:10),'(i4.4)') ijahr	!1     67
      write(filev10 ( 7:10),'(i4.4)') ijahr	!2     68
      write(filep   ( 7:10),'(i4.4)') ijahr	!3     69
      write(filetl  ( 7:10),'(i4.4)') ijahr	!4     70
      write(filesh  ( 7:10),'(i4.4)') ijahr	!5     71
      write(filecl  ( 7:10),'(i4.4)') ijahr	!6     72
      write(filepri ( 7:10),'(i4.4)') ijahr	!7     73

      write(fileussr( 7:10),'(i4.4)') ijahr	!8     74
      write(filedssr( 7:10),'(i4.4)') ijahr	!9     75
      write(fileulwr( 7:10),'(i4.4)') ijahr	!10    36
      write(filedlwr( 7:10),'(i4.4)') ijahr	!11    37
#endif
#ifdef COASTDAT
      write(metappend ( 1:4),'(i4.4)') ijahr    !1     67
#endif
c----------------------------------------------------------------
c    initial conditions and runoff input names,
c    attention! climatic conditions defined after 1998
c---------------------------------------------------------------

	T_S_bd	(5:8)='clim'
c	 if (ijahr.ge.1980.and.ijahr.le.1988)then
c	 T_S_bd(1:3)='T_S'
c	 else
c	 T_S_bd(1:3)='t_s'
c	 endif
c	write(T_S_bd  (5:8),'(i4)') ijahr
       write(sal_aa  (7:8),'(i2.2)')mjar
       write(z_d     (5:8),'(i4)') ijahr
       write(z_v     (5:8),'(i4)') ijahr

      rivinput	= 'rivinpyy_dis.dat'
      write(rivinput(7:8),'(i2.2)')  mjar

c      if(ijahr.ge.1998) rivinput = 'rivinp_k.dat'  ! for 1998-.... climatic means

c-----------------------------------------------------------------
c    Output file names
c-----------------------------------------------------------------

      write(outfile(4:5),'(i2.2)')    mjar
      write(control(1:3),'(a3)')   ppp
#ifdef MPI
      write(outfile1(4:5),'(i2.2)')    mjar
      write(control(11:13),'(i3.3)') myid
#endif
      write(jtext(1:2),'(i2.2)')      mjar

      write(outfile(1:3),'(a3)')   ppp
      write(gridinf(1:3),'(a3)')   ppp
      write(filradou(1:3),'(a3)')  ppp

c-----------------------------------------------
c     open control file, writing general info
c-----------------------------------------------
      nfucon = 7
      if(myid .eq. 0) write(*,*) pathpc2//pathout//control
      open(unit=nfucon,
     &file=pathpc2//pathout//control,form='formatted')
c-------------------------------------------------------

c------------Define model output file units------------------------------
c  daily means of model output parameters will be stored in Monthly files 'pppYYMM'
c  these Monthly output files will be opened in the monthly loop
c-----------------------------------------------------------------------
      nout2    = 285
      nout1    = 185
      nout     = 85	     !pppYYMM  physical variables,
      nout_b1  = 186
      nout_b   = 86		   !pppYYMMb biomass/concentration variables
      nout_f1  = 187
      nout_f   = 87		   !pppYYMMf production(bio fluxes) variables
      nout_s   = 88		   !pppYYMMs space bio fluxes variables


c-----------------------------------------------------------------------
c      print*,'	 open input files'
c-----------------------------------------------------------------------
      nfutop = 10
      open(unit=nfutop,file= pathpc//pathtop//topofile,form='formatted')

      IF(neust1)THEN
      nfuins = 13
      open(nfuins,file= pathpc//pathset//'inisalt',form='formatted')

      nfuint = 14
      open(nfuint,file= pathpc//pathset//'initemp',form='formatted')
      ENDIF !IF(neust1)THEN
#ifdef NCEP
      nfumet = 15
      open(unit=nfumet,file= pathpc//pathmet//metinde2,form='formatted'
     &,status='old')
#endif
      nfuriv_points = 16
      open(unit=nfuriv_points,file= pathpc//pathriv//rivpoint,
     &form='formatted')

      write(rivinput(10:12),'(a3)')  'dis'
      nfurin = 17
      open(unit=nfurin,file= pathpc//pathriv//rivinput,form='formatted')

      nfusa=76
      open(unit=nfusa,file= pathpc//pathzeb//sal_aa)  !salinity anomality Ute Daewel


c-----------------------------------------------------------------------
      IF(myid .EQ. 0)THEN
      print*,'	read salinityanomality		 file : ',sal_aa
      endif
c-----------------------------------------------------------------------
	read(nfusa,'(6(1x,f12.6))')delsal
	close(nfusa)
c-----------------------------------------------------------------------^M
      IF(myid .EQ. 0)THEN
	print*,delsal
      endif
c-----------------------------------------------------------------------
      IF(myid .EQ. 0)THEN
      print*,'	read modelltopo		  file : ',topofile
      endif
c-----------------------------------------------------------------------

      read(nfutop,2000) izet
 2000 format(18i4)
      close(nfutop)

c-----------------------------------------------------------------------
c     correct topography if necessary, blind out areas, dig channels
c-----------------------------------------------------------------------

      call topocorr(izet,iizet,ht,m,n,imax)

c-----------------------------------------------------------------------
c     compute topo-dependent grid values and
c     time-dependent grid values with SR kotief and SR gatmit
c-----------------------------------------------------------------------

      call kotief (dz,izet,nhor,ntot)
      call gatmit (deltat)
      call gittop (ltief,dzbod)
      call liq_boundary(nmodel)
      print*,'ntot',ntot
c---------------------------------------------------------------------
c   compile index for distributing coastdat fields to the model grid
c  coastdat is on 2km resolution
c if model grid is finer interpolation has to be installed
c----------------------------------------------------------------------
#ifdef COASTDAT
      call nc_read(pathpc3//pathmet//metpost//filetl//metappend,
     &"lon",
     &dump(:,:),start=[1,1],count=[ncdlam,ncdphi])
      ixc=0
       do icn=1,ncdlam
       do icm=1,ncdphi
       ixc=ixc+1
       cdlam(ixc)=dump(icn,icm)
        enddo
       enddo
      call nc_read(pathpc3//pathmet//metpost//filetl//metappend,
     &"lat",
     &dump(:,:),start=[1,1],count=[ncdlam,ncdphi])
      ixc=0
       do icn=1,ncdlam
       do icm=1,ncdphi
       ixc=ixc+1
       cdphi(ixc)=dump(icn,icm)
       enddo
       enddo
      ixc=0
c       call convert_metind(metind1,cdphi,cdlam,ncdlam,ncdphi,nera)
c-----------------------------------------------------  
      ifunc  = 5
      dmax   = 44000.
      pardum = 0.
      call intpork(weight,nrmess,cdlam,cdphi,
     *	      ifunc,pardum,dmax,ianzpkt,nera)

#endif

c-----------------------------------------------------------------------
      IF(myid .EQ. 0)THEN
      print*,'	read riverpoint data	  file : ',rivpoint
      ENDIF
c-----------------------------------------------------------------------

      if(river_on) then

      do ii=1,nanz
      read(nfuriv_points,2001) rrivnum,rindm,rindn,rivlat(ii),rivlon(ii)
        rivnum(ii)=nint(rrivnum)
        indm(ii)=nint(rindm) 
        indn(ii)=nint(rindn)
 2001	format(5e16.7)
c-----------------------------------------------------------------
      enddo
      close(nfuriv_points)
       endif   !	 IF(river_on) then
 

#ifdef NCEP
c-----------------------------------------------------------------------
      IF(myid .EQ. 0)THEN
      print*,'	interpolate met. data to model grid positions ***1'
      ENDIF
c-----------------------------------------------------------------------
	  do i=1,nera
	   read (nfumet,1001) rlamu(i),rphiu(i)
	  enddo
 1001 format(2(1X,f10.5))
      close(nfumet)
      ifunc  = 5
      dmax   = 170000.
      pardum = 0.
      call intpork(weight,nrmess,rlamu,rphiu,
     *	      ifunc,pardum,dmax,ianzpkt,nera)
     
      IF(myid .EQ. 0)THEN
      write(*,*) pathpc//pathmet//'position2_ncep.nsbs'
      ENDIF
      open(unit=nfumet,file= pathpc//pathmet//'position2_ncep.nsbs',
     &form='formatted',status='old')

c-----------------------------------------------------------------------
      IF(myid .EQ. 0)THEN
      print*,'	interpolate met. data to model grid positions ***2'
      ENDIF
c     different grid for air pressure
c-----------------------------------------------------------------------
      do i=1,nera2
      read (nfumet,1001) rlamup(i),rphiup(i)
      enddo
      close(nfumet)

      ifunc=5
      dmax=350000.
      dmax=290000. !Irina
      pardum=0.
      call intpork(weightp,nrmessp,rlamup,rphiup,
     *      ifunc,pardum,dmax,ianzpkt,nera2)
#endif

c-----------------------------------------------------------------------
c     set grid index arrays
c-----------------------------------------------------------------------
      lwe = 0
      nwet=0
       izet=0
      do k=1,n
	lwa = lwe+1
	lwe = indend(k)
	do lw=lwa,lwe
	  i = iwet(lw)
	  lump = lazc(lw)
	  jjc(i,k) = 1
	  iindex(i,k) = nwet
	  id3sur(lw) = nwet+1
	  do jj=1,lump
	    nwet = nwet+1
	  enddo
	  izet(i,k) = lw
	enddo
      enddo

#ifdef MPI
       CALL domain_decomp(lazc,iwet)
       if(iyears .eq. 1)then
       CALL INDEX_MPI(icord)
       endif
       CALL temp_izet(izet)
       CALL local_izet(izet,izets,izete,icord)
       if(nPh .gt. 1) CALL ver_comm_cord(lazc,izet,icord)
       if(nPv .gt. 1) then
       CALL hor_comm_b22(icord,izet,lazc)
       CALL hor_comm_t22(icord,izet,lazc)
       CALL hor_comm_b11(icord,izet,lazc)
       CALL hor_comm_t11(icord,izet,lazc)
       CALL hor_comm_b00(icord,izet,lazc)
       CALL hor_comm_t00(icord,izet,lazc)
       ctemp=0.0
       ictemp=0
      iel=1
      iel1=1
      do j=1,n
      do i=1,m
      dtemp(i,j)=0.0
      idtemp(i,j)=0
      enddo
      enddo
      endif
!      IF(myid .EQ. 0)THEN
!      open(unit=144,file='icord')
!      DO I=1,nprocs
!      write(144,*)icord(I,1),icord(I,2),icord(I,3),icord(I,4)
!      ENDDO
!      close(144)
!      write(*,*)'ldep(izet(97,10)):',ldep(izet(97,10))
!      write(*,*)'ldep(izet(127,28)):',ldep(izet(127,28))
!      ENDIF
#endif
c-----------------------------------------------------------------------
c     exchange coefficientis, eddy viscosity
c	   pav= vertical coefficient
c	   pah= horizontal coefficient
c-----------------------------------------------------------------------

	 data pav/ilo*50.e-4/
	  do i=1,ilo
	  pah(i)=read_pah
	  enddo

      pah_fac=1.

c-----------------------------------------------------------------------
c     set vertical constants and exchange coefficients
c-----------------------------------------------------------------------

      do j=1,ilo
	dc(j)  = dz(j)
	rhoq(j)= 0.027
	r2d(j) = 0.0
	ad(j)  = pah(j)
	pah(j) = ad(j)
	av(j)  = pav(j)
      enddo

      gh      = 0.5*g
      dd(1)   = dz(1)
      dz1     = dz(1)
      dh(1)   = gh*dd(1)
      pd(1)   = dd(1)
      prd(1)  = 1./dd(1)
      pr2d(1) = 1.0

      do j=2,ilo
	dd(j)	= dz(j)-dz(j-1)
	dh(j)	= gh*dd(j)
	pd(j)	= dd(j)
	prd(j)	= 1./dd(j)
	pr2d(j) = 1./(dd(j-1)+dd(j))
      enddo

      do j = 1,n
      do i = 1,m
	dx(i,j)=dln(i)
	dy(i,j)=dl
      enddo
      enddo

c-----------------------------------------------------------------------
c     compute reference density and pressure
c-----------------------------------------------------------------------

      dzz = 0.0
      sal = 33.0
      tem = 12.
      do k=1,ilo
	pres(k)	  = (dzz+0.5*pd(k))*0.1	      !pres(k)0=>dz(k-1)&dz(k)
	refrho(k) = sigma(sal,tem,pres(k))
	rhoq(k)	  = 0.999842594+refrho(k)
	dzz	  = dz(k)
      enddo

      do l = 1,ntot
	sac(l) = sal
	tec(l) = tem
      enddo


c-----------------------------------------------------------------------
c     compute space depending omega
c-----------------------------------------------------------------------

      ahaz    = 500.*dt*dlr*dlr
      omi = 1.3
      oma = 1.65
      tomi = 300.
      delom = (oma-omi)/(dz(ilo)-tomi)

      ik = 0
      do k=3,nz
	lwa = isornr(k)
	lwe = isorsr(k)
	do lw=lwa,lwe
	  is = lazc(lw)
	  tif = ldep(lw)+dz(is)-dd(is)-tomi
	  tif = amax1 (tif,0.)
	  i = iwet(lw)
	  ik = ik+1
	  isor(ik) = i
	  ksor(ik) = k
	  omega(ik) = omi+tif*delom
	  zdif(i,k) = omega(ik)
	enddo
      enddo

      nosor = ik
      IF(myid .EQ. 0)THEN
       write(*,*)nosor,'nosor'
      write(7,7028)nosor
 7028 format(' anzahl der zu iter. gitterpunkte (kasor): ',i8)
      endif
c--------------------------------------------------------------------------
c    Red-Black ordering of SOR-points
c--------------------------------------------------------------------------

      do 180 kodd = 1,2
      do 180 iodd = 1,2
      nrbsor(iodd,kodd) = 0
      do 180 ik = 1,nosor
      i = isor(ik)
      k = ksor(ik)
      if (mod (i,2).ne.iodd-1 .or. mod (k,2).ne.kodd-1) goto 180
      nrbsor(iodd,kodd) = nrbsor(iodd,kodd)+1
      irbsor(nrbsor(iodd,kodd),iodd,kodd) = ik
  180 continue
      if (nrbsor(1,1)+nrbsor(1,2)+nrbsor(2,1)+nrbsor(2,2).ne.nosor)
     1 stop 'error: wrong red-black ordering of sor-points'

      do nena=1,2
      do ik=1,nosor
	i = isor(ik)
	k = ksor(ik)
	zdif(i,k) = omega(ik)+ahaz*
     *	 (zdif(i-1,k)+zdif(i+1,k)+zdif(i,k-1)+zdif(i,k+1)-4.*zdif(i,k))
	omega(ik) = zdif(i,k)
      enddo
      enddo
#ifdef MPI
      call NEW_ISORNR(isornr,isorsr,iwet)
#endif
       nout_start=89

c--------------------------------------------------------------------
c==============	 READING START FILES ================================
c--------------------------------------------------------------------

	if(ibup.eq.0)then !start from end of month, unformatted BINARY
*1===================READ PHYSICS=================================
	  nout = 85	     ! physical	    variables Monthly output
	  outfile(1:3)=ppp
!	   if(ijahr.eq.2001)outfile(1:3)='y02'

	  write(outfile(4:5),'(i2.2)') mjar
	  write(outfile(6:7),'(i2.2)') mona-1

	  open(unit=nout,
     &	  file=pathpc2//pathout//outfile(1:7),form='unformatted')
	  IF(myid .EQ. 0)THEN
	  print*,'read pyhsics monthly file open: ',outfile(1:7)
	  ENDIF
	  do nday=1,itagan(mona-1)
	  call IO_daily_physics(nout,mjar,lmon,nday,ivier,imal,2) ! nnn write=1 or read=2
	  enddo
	  call IO_start_physics(nout,mjar,lmon,nday,ivier,imal,2) ! nnn write=1 or read=2
      close(nout)

c-----------------------------------------------------------------------------
c     zero initialisation
c----------------------------------------------------------------------------
	do i=1,ndrei
	  umit(i)  = 0.0
	  vmit(i)  = 0.0
	  wcmit(i) = 0.0
	  acmit(i) = 0.0
	  szmit(i) = 0.0
	  tcmit(i) = 0.0
	  scmit(i) = 0.0
	enddo

	do i=1,m
	do j=1,n
	  frimit(i,j)  = 0.0
	  hismit(i,j)  = 0.0
	  hisrmit(i,j) = 0.0
	  tismit(i,j)  = 0.0
	  uimit(i,j)   = 0.0
	  vimit(i,j)   = 0.0
	  fqgmit(i,j)  = 0.0
	  qois(i,j)    = 0.0
	  qiis(i,j)    = 0.0
	enddo
	enddo

	do i=1,khor
	  zmit(i)    = 0.0
	  fqrmit(i)  = 0.0
	  fqsmit(i)  = 0.0
	  fqlmit(i)  = 0.0
	enddo


      else	  !ibup=1  read backup/start files
      print*,'read start files'
c      if(iyears.eq.1)then
      if(ijahr.eq.1990)then
c      if(ijahr.eq.1979)then
#ifdef SASCII
      print*,'start from ascii file'
      isice=0
       goto 6574 ! start from ASCII physics start file
#else
      startfile(1:3)='cd3'
c      startfile(1:3)=ppp
      isice=1
#endif
      else
       startfile(1:3)=ppp
      isice=1
      endif

*21===================READ PHYSICS=================================
c------------------------------------------------------------------------------
c	read binary physics start files
c-----------------------------------------------------------------------------
      if(iyears.lt.spinup+1)then
      write(startfile(4:5),'(i2.2)') mjar
      startfile(6:9)	      ='1231'
      open(nout_start,file=
     &pathpc2//pathout//'/'//startfile(1:9)//'_start'
     &,form='unformatted')
      else
      write(startfile(4:5),'(i2.2)') mjar_1
      startfile(6:9)='1231'

      open(nout_start,file=
     &pathpc2//pathout(1:14)//ijahr_1//'/'//startfile(1:9)//'_start'
     &,form='unformatted',status='old',access='stream')
      endif
c      IF(myid .EQ. 0)THEN
      write(*,*)
     &pathpc2//pathout(1:14)//ijahr_1//'/'//startfile(1:9)//'_start'
c      ENDIF
      call IO_start_physics(nout_start,mj,ilm,indl,iivier,iimal,2) ! nnn write=1 or read=2
      close(nout_start)

6574  continue
c-----------------------------------------------------------------------------
c      start from ASCII physics start file
c-----------------------------------------------------------------------------
       nfuset=127
       startfile(1:9)='clmYY1231'

	write(startfile(4:5),'(i2.2)') 47
	IF(myid .EQ. 0)THEN
       write(*,*)'!!!ibup=1 file era= ', pathpc//pathtsb//startfile(1:9)
	ENDIF

      open(unit=nfuset,file=pathpc//pathtsb//startfile(1:9),
     &status='old')
      read (nfuset,2003,err=9876,end=9874) imald,mjard,monatd,ltagd
      read (nfuset,2004,err=9876,end=9874) zalt,z,zac
      read (nfuset,2004,err=9876,end=9874) txc,tyc,pac
      read (nfuset,2006,err=9876,end=9874) (uc(l),l=1,ntot)
      read (nfuset,2006,err=9876,end=9874) (vc(l),l=1,ntot)
      read (nfuset,2004,err=9876,end=9874) (wc(l),l=1,ntot)
      read (nfuset,2004,err=9876,end=9874) (avc(l),l=1,ntot)
      read (nfuset,2004,err=9876,end=9874) (sac(l),l=1,ntot)
      read (nfuset,2004,err=9876,end=9874) (tec(l),l=1,ntot)
      read (nfuset,2004,err=9876,end=9874) frice,his,hisr,tis
      read (nfuset,2004,err=9876,end=9874) ui,vi
      close(unit=nfuset)

 2003	format(4i5)
 2004	format(8(1x,e12.6))
 2005	format(5i5)
      endif	  !ibup=1

      goto 9875

 9876 print*,'check: backup-file error',nfuset
      stop'abnormal end'
 9874 print*,'check: backup-file end',nfuset
      stop'abnormal end'
 9875 continue

#ifdef MPI
      call deco1d2d_s(z2d,zac,iwet,indend,1.)
#else
      call deco1d2d(z2d,zac,iwet,indend,1.)
#endif
!      call estate
      goto 8974
 8973 continue
c-----------------------------------------------------------------------
      IF(myid .EQ. 0)THEN
      print*,'	cold start		  mmdd : ',mona,itanf
c-----------------------------------------------------------------------
      write(7,8757)
 8757 format(1x,'cold start !!!!',/)
      ENDIF
 8974 continue



c------------------------------------------------------------------
c    vertical grid and corresponding help arrays
c-----------------------------------------------------------------

       ilo2=ilo+2
      lwe = 0
	do j=1,n
	lwa = lwe+1
	lwe = indend(j)
	do lw=lwa,lwe
	  lump = lazc(lw)
	  i = iwet(lw)
	    pd2(i,j,0)=-zac(lw)			! w1-point vertikal positions
	    pd2(i,j,1)=dz(1)+zac(lw)		! t1-points vertikal positions
	    dzt(i,j,1)=(dz(1)-zac(lw))*0.5	! w-points vertikal positions
	do k=2,lump-1
	  pd2(i,j,k)=dz(k)-dz(k-1)			   !tk-points vertikal positions
	  dzt(i,j,k)=dz(k-1)+(dz(k)-dz(k-1))*0.5

	enddo  !k
	    pd2(i,j,lump)=dzbod(i,j)			   !bottom layer thickness
	    if(lump.eq.1)then
		pd2(i,j,lump+1)=dzbod(i,j)		   !real bottom depth
		dzt(i,j,lump)=dzbod(i,j)*0.5
		else
		pd2(i,j,lump+1)=dz(lump-1)+dzbod(i,j)	   !real bottom depth
		dzt(i,j,lump)=dz(lump-1)+dzbod(i,j)*0.5
	    endif
		pd2(i,j,ilo2)  =lump			   ! lump


	enddo  !i
       enddo   !j


      IF(neust1)THEN
c========== cold start ==================================
c--------------------------------------------------------
c     read temperature and salinity initial conditions
c--------------------------------------------------------
c-----------------------------------------------------------------------
      IF(myid .EQ. 0)THEN
      print*,'read setup t,s file:','inisalt',' ','initemp'
      ENDIF
c-----------------------------------------------------------------------

	call readtsini (nfuins,nfuint,sac,tec,ndrei)

	close(unit=nfuins)
	close(unit=nfuint)
c-----------------------------------------------------------------------
c     check if the T,S-array is correctly filled
c-----------------------------------------------------------------------
      icho = 0
#ifdef MPI
      call deco1d3d_s (dudo3,sac,ndrei)
#else
      call deco1d3d (dudo3,sac,ndrei)
#endif
      call deco1d2di(idudo2,lazc,iwet,indend,1)

      do i=2,m
      do j=2,n
      do k=1,idudo2(i,j)
      if(dudo3(i,j,k).eq.-99.0)then
      IF(myid .EQ. 0)THEN
       print*,'i,j,k, salt(i,j,k) ',i,j,k,dudo3(i,j,k)
      ENDIF
	    if(k.ne.1)then
		   dudo3(i,j,k) = dudo3(i,j,(k-1))
	    else
		   dudo3(i,j,k) = dudo3(i,j-1,k)
		   if(dudo3(i,j-1,k).le.0.)dudo3(i,j,k) = dudo3(i-1,j,k)
	    endif
      endif
      enddo
      enddo
      enddo

      IF(myid .EQ. 0) print*,'holes ',icho
#ifdef MPI
      call comp3d1d_s(dudo3,sac)
#else
      call comp3d1d (dudo3,sac)
#endif
	IF(myid .EQ. 0)THEN
      do iu=1,ndrei
	if(sac(iu).gt.35.5) print*,'salt too much'
	enddo
	endif
c ###	  print*,'salt bottom gotl ',dudo3(87,147,(idudo2(87,147)))

	do i=1,m
	do j=1,n
	do k=1,ilo
	  dudo3(i,j,k) = -99.9
	enddo
	enddo
	enddo
#ifdef MPI
	call deco1d3d_s(dudo3,tec,ndrei)
#else
	call deco1d3d (dudo3,tec,ndrei)
#endif
	do j=2,n
	do i=2,m
	do k=1,idudo2(i,j)
	 if(dudo3(i,j,k).eq.-99.0)then
          if(k.ne.1)then
          dudo3(i,j,k) = dudo3(i,j,(k-1))
           else
          dudo3(i,j,k) = dudo3(i,j-1,k)
          if(dudo3(i,j-1,k).eq.-99.0)dudo3(i,j,k) = dudo3(i-1,j,k)
           endif
        endif
        enddo
        enddo
        enddo
        IF(myid .EQ. 0)THEN
        do i=1,m
        do j=1,n
         do k=1,idudo2(i,j)
	    if(dudo3(i,j,k).le.-1.)
     &	  print*,'catastrophy in initial:T= new ',dudo3(i,j,k),i,j,k
	 enddo
	enddo
	enddo
	endif
#ifdef MPI
	call comp3d1d_s (dudo3,tec)
#else
	call comp3d1d (dudo3,tec)
#endif
!	 call estate


c=====End of cold start ==================================
      ENDIF !	    IF(neust1)THEN

!c-----------------------------------------------------------------------
	IF(myid .EQ. 0)THEN
	print*,'  setup dynamic ice model  '
	endif
!c-----------------------------------------------------------------------

	call setice(iindex,izet,m,n)

	do j=1,ny
	do i=1,nx
          gairx(i,j) = 0.0
          gairy(i,j) = 0.0
        enddo
        enddo
c        goto 1992
c        if(iyears.gt.1)then
c        call estate
        do i=1,m
        do j=1,n
           do ii=1,jjc(i,j)
        txx(i,j)=txc(izet(i,j))
        tyy(i,j)=tyc(izet(i,j))
        enddo
        enddo
        enddo
cc--------------------------------------------------------------
cc		wind stress for ice
cc--------------------------------------------------------------
c#ifdef MPI
c	do j=JB1,JB2
c	do i=IB1,IB2
c#else
c	do i=1,m
c	do j=1,n
c#endif
c	   wgesch(i,j) = sqrt(txx(i,j)**2+tyy(i,j)**2)
c	   cd = 1.28*(0.8+0.065*wgesch(i,j))
c	   cd = cd*1.e-3
c	   txx(i,j) = cd*wgesch(i,j)*txx(i,j)
c	   tyy(i,j) = cd*wgesch(i,j)*tyy(i,j)
c	enddo
c	enddo
c        endif
c        goto 1992
        
         lwa = lb(1)
         lwe = le(n)
         if (lwa.le.lwe) then
            llb = indwet(lwa)+1
            lle = indwet(lwe)+lazc(lwe)
            do ll = llb,lle
               k = indver(ll)
               stc(ll) = sigma(sac(ll),tec(ll),pres(k))-refrho(k)
               if (k.gt.1) stpc(ll) = sigma(sac(ll),tec(ll),pres(k-1))-
     1          refrho(k-1)
            end do
         end if
        call stressa(iindex,izet,jjc,z2d,m,n,dz,ilo,txx,tyy)
c----------------------------------------------
	iffir = 0
	call first (iizet,m,n,iffir,isice)
	iffir = 1
        if(.not. withice)then
	do j=1,n
	 do i=1,m
	   ui(i,j)    = 0.0
	   vi(i,j)    = 0.0
	enddo
	enddo
        endif


c-----------------------------------------------------------------------
      IF(myid .EQ. 0)THEN
      print*,'	write grid info file : ',gridinf
      endif
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
       IF(myid .EQ. 0)THEN
       print*,'  skip riv. input if start month ne 1 file : ',rivinput
       endif
c-----------------------------------------------------------------------

         if(river_on)then
        do i=1,mona-1
        itagend = itagan(i)
        do ii=1,itagend
          call read_riv_input(nfurin,rivinp,nanz)
        enddo
        IF(myid .EQ. 0)THEN
        write(*,*) 'fresh runoff skipped YY MM',ijahr,i
        endif
        enddo
	endif !	       if(river_on)then

      write(outfile(4:5),'(i2.2)') mjar

#ifndef MPI
	  outfile(1:3)=ppp
	  open(121,
     &	  file=pathpc2//pathout//'Qh_'//outfile(1:5),form='unformatted')
#else
	  outfile(1:3)=ppp
	  write(outfile(4:5),'(i2.2)') mjar
	  if(myid .eq. nprocs-1) open(121,
     &	  file=pathpc2//pathout//'Qh_'//outfile(1:5),form='unformatted')
	  outfile1(1:3)=ppp
	  write(outfile1(4:5),'(i2.2)') mjar
#endif
c***********************************************************************
c  begin monthly loop
      IF(myid .EQ. 0)THEN
      print*,' '
      print*,'Here we go: starting with time loop'
      print*,' '
      endif
c***********************************************************************


cKB      do 2111 lmon=mona,mone
      do 2111 lmon=1,1
      igstep1 = 0
      it_init = 0
      monat = lmon
      avmax=1000.e-4 !from 2002
c      if(lmon.eq.1)avmax=1400.e-4 !from 2002  y m d h m s jul	2002  1 28  7 30  0  28
c      if(lmon.eq.1)avmax=600.e-4 !from 2002  y m d h m s jul  2002  1 28  7 30	 0  28

c------------------------------------------------------------------------
c    Output units, daily means in monthly files, unformatted BINARY
c    Open output files
c-----------------------------------------------------------------------

      nout     = 85		   ! physical	  variables Monthly output
      nout_b   = 86		   ! biomass	  variables Monthly output
      nout_f   = 87		   ! production	  variables Monthly output
      nout_s   = 88		   ! space fluxes variables Monthly output
#ifndef MPI
	  outfile(1:3)=ppp
	  write(outfile(6:7),'(i2.2)') monat
	  open(unit=nout,
     &	  file=pathpc2//pathout//outfile(1:7),form='unformatted')
	  print*,'output pyhsics monthly file open: ',outfile(1:7)

#ifdef _FABM_
          outfile(8:8)='b'
          open(unit=nout_b,
     &    file=pathpc2//pathout//outfile,form='unformatted')
          print*,'output biomass monthly file open: ',outfile
#endif

#else
	 write(outfile1(6:8),'(i3.3)') myid
	 open(122,
     &	 file=pathpc2//pathout//'Qh_'//outfile1(1:8),form='unformatted')
	 outfile(1:3)=ppp

	 write(outfile(6:7),'(i2.2)') monat


	  if(myid .eq. 0) then
	  open(unit=nout,
     &	  file=pathpc2//pathout//outfile(1:7),form='unformatted')
	  endif

	  write(outfile(6:8),'(i3.3)') myid

	  open(unit=nout1,
     &	  file=pathpc2//pathout//outfile(1:8),form='unformatted')

#ifdef _FABM_
          if(myid .eq. 1)then
	  write(outfile(6:7),'(i2.2)') monat
          outfile(8:8)='b'
          open(unit=nout_b,
     &    file=pathpc2//pathout//outfile(1:8),form='unformatted')
          endif

          write(outfile(6:8),'(i3.3)') myid

          outfile(9:9)='b'
          open(unit=nout_b1,
     &    file=pathpc2//pathout//outfile(1:9),form='unformatted')
#endif

#endif
c-----------------------------------------------------------
c     open files for boundary condtions
c     z_d, z_v	and  t,s boundary condition  unit 18, 19, 65
c-----------------------------------------------------------
      write(z_d	   (9:10),'(i2.2)') monat
      write(z_v	   (9:10),'(i2.2)') monat
      write(T_S_bd (9:10),'(i2.2)') monat
c      write(bio_bd (9:10),'(i2.2)') monat

	nzactp =  18
      open(unit=nzactp,file=pathpc//pathzeb//z_d,form='formatted',
     & status='old')
	nzac =	  19
      open(unit=nzac,file= pathpc//pathzeb//z_v,form='formatted',
     & status='old')!20m tidal variation

      if(iskipd.gt.0)then
	     IF(myid .EQ. 0)THEN
	     write(*,*)iskipd
	     endif
       iread=(itanf-1)*24*3
       do itagvo=1,iread ! every 20 min
	      IF(myid .EQ. 0)THEN
              print*,'read z_v', iskipd, itanf, iread
	      write(*,*)itagvo
	      endif
	  call readzdv(zvar,lrp,nzac) !z_v  every hour
         
       enddo
	  iskipd=0
	   IF(myid .EQ. 0)THEN
	   write(*,*)'readzv is skipped, day=',itanf
	   endif
      endif
	nfutsb1 = 65
      open(65,	     file= pathpc//nmodel//'/input_1980/'//T_S_bd) !if clim

	  IF(myid .EQ. 0)THEN
	  write(*,*)pathpc//pathzeb//z_d
	  write(*,*)pathpc//pathzeb//z_v
	  write(*,*)pathpc//nmodel//'/input_1980/'//T_S_bd
	  ENDIF
#ifdef COASTDAT
      write(metappend ( 6:7),'(i2.2)') monat    !1     67
#endif
#ifdef NCEP
c-----------------------------------------------------------------------
c     open monthly meteorologic files
c-----------------------------------------------------------------------

      write(fileu10 (11:12),'(i2.2)') monat	!1     67
      write(filev10 (11:12),'(i2.2)') monat	!2     68
      write(filep   (11:12),'(i2.2)') monat	!3     69
      write(filetl  (11:12),'(i2.2)') monat	!4     70
      write(filesh  (11:12),'(i2.2)') monat	!5     71
      write(filecl  (11:12),'(i2.2)') monat	!6     72
      write(filepri (11:12),'(i2.2)') monat	!7     73

      write(fileussr(11:12),'(i2.2)') monat	!8     74
      write(filedssr(11:12),'(i2.2)') monat	!9     75
      write(fileulwr(11:12),'(i2.2)') monat	!10	36
      write(filedlwr(11:12),'(i2.2)') monat	!11    37
      IF(myid .EQ. 0)THEN
      write(*,*)'open meteo file  : ', pathpc//pathmet//fileu10
      ENDIF
      open(67,file=pathpc//pathmet//fileu10 ,status='old',readonly)
      open(68,file=pathpc//pathmet//filev10 ,status='old',readonly)
      open(69,file=pathpc//pathmet//filep   ,status='old',readonly)
      open(70,file=pathpc//pathmet//filetl  ,status='old',readonly)
      open(71,file=pathpc//pathmet//filesh  ,status='old',readonly)
      open(72,file=pathpc//pathmet//filecl  ,status='old',readonly)
      open(73,file=pathpc//pathmet//filepri ,status='old',readonly)
c      do ii=67,73
c      read(ii,'(a81)') cdum
c      enddo
      if (radinp) then
      open(74,file=pathpc//pathmet//fileussr ,status='old',readonly)
      open(75,file=pathpc//pathmet//filedssr ,status='old',readonly)

      open(36,file=pathpc//pathmet//fileulwr ,status='old',readonly)
      open(37,file=pathpc//pathmet//filedlwr ,status='old',readonly)
c      do ii=36,37
c      read(ii,'(a81)') cdum
c      read(ii+38,'(a81)') cdum
c      enddo
      endif !if (radinp) then

c-------------------------------------------------------------
c    read climatic corrections for the atmospheric forcing
c-------------------------------------------------------------

      if (radinp) then
      filetkor='2tnc_ec_MM.dat'
      write(filetkor(9:10),'(i2.2)') monat

      filetkor(1:8)  = '2tnc_ec_MM.dat'
      open (97,file=pathpc//pathmet//filetkor,status='old')
      read(97,'(a50)') cdum
      do jm=1,m
	do jn=1,n/5
	  read(97,2020) (tempkorr(jm,kk),kk=(jn-1)*5+1,jn*5)
	enddo
	if (mod(n,5).ne.0) then
	  ihilf=n/5
	  read(97,2020) (tempkorr(jm,kk),kk=ihilf*5+1,n)
	endif
      enddo
2020  format(5(1X,E15.8))
      close(97)

      filetkor(1:8)  = 'ssrnc_ecMM.dat'
      open (97,file=pathpc//pathmet//filetkor,status='old')
      read(97,'(a50)') cdum
      do jm=1,m
	do jn=1,n/5
	  read(97,2020) (ssrkorr(jm,kk),kk=(jn-1)*5+1,jn*5)
	enddo
	if (mod(n,5).ne.0) then
	  ihilf=n/5
	  read(97,2020) (ssrkorr(jm,kk),kk=ihilf*5+1,n)
	endif
      enddo
      close(97)

      filetkor(1:8)  = 'strnc_ecMM.dat'

      open (97,file=pathpc//pathmet//filetkor,status='old')
      read(97,'(a50)') cdum
      do jm=1,m
	do jn=1,n/5
	  read(97,2020) (strkorr(jm,kk),kk=(jn-1)*5+1,jn*5)
	enddo
	if (mod(n,5).ne.0) then
	  ihilf=n/5
	  read(97,2020) (strkorr(jm,kk),kk=ihilf*5+1,n)
	endif
      enddo
      close(97)


      endif !if (radinp) then
#endif  !ncep

      

c-----------------------------------------------------------------------
c     day loop (days per month)
c-----------------------------------------------------------------------

      if(lmon.eq.mone)then
	itagend = itage
      else
	itagend = itagan(monat)
      endif

      do 1111 nday=itaga,itagend
      ltag = nday
      jday=jday+1

c-----------------------------------------------------------------------------
c   read values of boundary conditions and river nutrient loads
c-----------------------------------------------------------------------------
c   1. for daily z_d and monthly S and T:
      call read_bconditions (nfutsb1,nzactp,nday) !=>zday(lrp),sarum(lvrp),tarum(lvrp)
c     read river discharge daily      
      if(river_on) then
          call Rzero(rivinp,nanz)
          call read_riv_input(nfurin,rivinp,nanz)
          eintr=rivinp
      if(myid.eq.0)print*,'eintr ',eintr
      endif         

c   3. set necessary equal boundary lines for T and S
      call tsrneu(iindex,izet,time,dz,delsal)  !set equal boundary lines for T,S

      call estate
c-----------------------------------------------------------------------
c     to monitor model CPU time
c-----------------------------------------------------------------------
      call ftrace_region_end ('rosh')
      call ftrace_region_begin ('time_loop')

c-----------------------------------------------------------------------
c     date loop (dates per day)
c-----------------------------------------------------------------------

      do 1010 ivier=1,iviert

      do l = 1,ngro
	gross(l) = 0.0
      enddo

c-----------------------------------------------------------------------
c     dt-loop (dt per date)
c-----------------------------------------------------------------------

      do 1000 ip=1,iper

      icccount=icccount+1
c#ifdef MPI
c      if(progno)then
c      call INTI_COMM_VAR
c      do I=3,nbio-1
c      call com_vertical(Tc(:,I),ndrei,dtemp,iel,2,2,1,3,0,lazc)
c      enddo
c      I=nbio
c      call com_vertical(Tc(:,I),ndrei,dtemp,iel,2,2,1,3,1,lazc)
c      do I=3,nbio-1
c      call com_horizontal(Tc(:,I),ndrei,dtemp,iel,2,2,1,3,0,icdepb0,
c     * idepb0,icdept0
c     *,idept0,ncsendb0,nsendb0,ncsendt0,nsendt0,nbhlayer0,nthlayer0)
c      enddo
c      I=nbio
c      call com_horizontal(Tc(:,I),ndrei,dtemp,iel,2,2,1,3,1,icdepb0,
c     *idepb0,icdept0
c     *,idept0,ncsendb0,nsendb0,ncsendt0,nsendt0,nbhlayer0,nthlayer0)
c      endif
c#endif
      imal = imal+1

      call settime(time,ijulu,ihouu,iminu,isecu)

c-----------------------------------------------------------------------
*     time output on screen
c-----------------------------------------------------------------------
       if(myid .eq. 0) THEN
	if(mod(ip-1,1).eq.0)
     & write(*,2829)ijahr,lmon,nday,ihouu,iminu,isecu,ijulu
       endif
 2829 format('y m d h m s jul ',i5,5i3,i4)


c-----------------------------------------------------------------------^M
c	read and process meteorological data, every 6 hours
c-----------------------------------------------------------------------^M
	 if  (mod(ip+ireadmet-1,ireadmet).eq.0) then
        imetind=imetind+1
c-----------------------------------------------------------------------
	if(myid .EQ. 0)
     &	print*,'  read met. data	    file : '
c-----------------------------------------------------------------------
#ifdef NCEP
        call readmet (erau10,nera ,67,ijahr)
	call readmet (erav10,nera ,68,ijahr)
	call readmet (erap  ,nera2,69,ijahr)
	call readmet (eratl ,nera ,70,ijahr)
	call readmet (erash ,nera ,71,ijahr)
	call readmet (eracl ,nera ,72,ijahr)
	call readmet (erapri,nera ,73,ijahr)

       if (radinp) then
	call readmet (eraussr,nera,74,ijahr)
	call readmet (eradssr,nera,75,ijahr)
	call readmet (eraulwr,nera,36,ijahr)
	call readmet (eradlwr,nera,37,ijahr)
	  do i=1,nera
	  erassr(i)=eradssr(i)-eraussr(i)
	  erastr(i)=eradlwr(i)-eraulwr(i)
	  enddo
       endif ! if (radinp) then
        if (myid.eq.0)then
         print*,'test readmet temp max min:',maxval(eratl),minval(eratl)
       endif  

	   if((lmon.eq.mona).and.(nday.lt.itanf)) goto 1000 ! skip dates if another start than 1.1


c--------------------------------------------------------------------------------------------------
c	 interpolation of atmospheric data to the model grid and removal of boundary gradients
c	 in case of nday lt itanf no calculation, but just reading until start time of calculation
c--------------------------------------------------------------------------------------------------
c----------------decreas wind in case of model crush---------------------------

       
c-----------------------------------------------------------------------
         if(myid .eq. 0) print*,'interpolation of met. data'
c-----------------------------------------------------------------------

        do ib=1,nera
c----------------decreas wind in case of model crush---------------------------
         if(ijahr.eq.1989.and.monat.eq.2)then
          erau10(ib)=0.9*erau10(ib)  !max in data -8.2  21.12
          erav10(ib)=0.9*erav10(ib)  !max in data -8.2  21.12
         endif
c------------------------------------------

	  bwind(ib)=sqrt(erau10(ib)**2+erav10(ib)**2)
	enddo
	call interpol2(bwindn,bwind,weight,nrmess,ianzpkt,nera)
	call bforce(bwindn)
	call interpol2(tyy   ,erav10,weight,nrmess,ianzpkt,nera)
	call bforce (tyy)
	call interpol2(txx   ,erau10,weight,nrmess,ianzpkt,nera)
	call bforce (txx)
	call interpol2(pat,erap ,weightp,nrmessp,ianzpkt,nera2) ! pressure
	call bforce (pat)

c----------------------------------------------------------------
c     correction of wind from wind speed and direction after interpolation on model grid
c----------------------------------------------------------------
#ifdef MPI
	do inu=JB1,JB2
	do imu=IB1,IB2
#else
	do imu =1,m
	do inu=1,n
#endif
	  call uv2sr(txx(imu,inu),tyy(imu,inu),speedu,dirou,0,1)
	  diru(imu,inu)=dirou
	  call uv2sr(txxu,tyyu,bwindn(imu,inu),diru(imu,inu),1,1)
	  txx(imu,inu)=txxu
	  tyy(imu,inu)=tyyu
	enddo
	enddo

c----------------------------------------------------------------
c     temperature and spec. humidity, precipitation
c----------------------------------------------------------------

        call interpol2(dumm,eratl,weight,nrmess,ianzpkt,nera)
#ifdef MPI
	do j=JD1,JD2
	do i=ID1,ID2
#else
        do i=1,m
        do j=1,n
#endif
	    dumm(i,j)=dumm(i,j)+tempkorr(i,j)
	enddo
	enddo
	   call bforce (dumm)
#ifdef MPI
	do j=JB1,JB2
	do i=IB1,IB2
#else
	do i=1,m
	do j=1,n
#endif
	do ill=1,jjc(i,j)
	  taair(izet(i,j)) = dumm(i,j)
cccccc b08 sensitivity airtemperature for use in IBM
c          taair(izet(i,j)) =taair(izet(i,j))+1.2
c--------------------------------------------------- 

	enddo
	enddo
	enddo

	call interpol2(dumm,erash,weight,nrmess,ianzpkt,nera)
	   call bforce (dumm)

#ifdef MPI
	do j=JB1,JB2
	do i=IB1,IB2
#else
	do i=1,m
	do j=1,n
#endif
	do ill=1,jjc(i,j)
	  shum(izet(i,j)) = dumm(i,j)
	enddo
	enddo
	enddo
#ifdef MPI
	do j=JB1,JB2
	do i=IB1,IB2
#else
	do i=1,m
	do j=1,n
#endif
	do ill=1,jjc(i,j)
	  call shumtaup(shum(izet(i,j)),pat(i,j),taup(izet(i,j)))
	   taup(izet(i,j))=taup(izet(i,j))+tempkorr(i,j)
	enddo
	enddo
	enddo
*xt4,xc4---------------------------
c	   do lw=1,khor
c	   taair(lw)=taair(lw)+3. !Tair increase of 3deg ppp='xt0'
c	   taup(lw)=taup(lw)+3.
c	   enddo
*-----------------------------------



	call interpol2(dumm,erapri,weight,nrmess,ianzpkt,nera)
	 call bforce (dumm)
#ifdef MPI
	do j=JB1,JB2
	do i=IB1,IB2
#else
	do i=1,m
	do j=1,n
#endif
	do ill=1,jjc(i,j)

	  precip(izet(i,j)) = dumm(i,j)	  ![kg/m**2/s] = [mm/s]
	 if (precip(izet(i,j)).lt.0.0) precip(izet(i,j))=0.
	enddo
	enddo
	enddo

c-----------------------------------------------------------------------
c			 Radiation
c-----------------------------------------------------------------------
	if (radinp) then
c------------------------------------------------------------------------
c	 shortwave radiation taken from ncep
c------------------------------------------------------------------------
	  call interpol2(dumm,erassr,weight,nrmess,ianzpkt,nera)
#ifdef MPI
	do j=JD1,JD2
	do i=ID1,ID2
#else
	do i=1,m
	do j=1,n
#endif
	    dumm(i,j)=dumm(i,j) +ssrkorr(i,j)
	enddo
	enddo
	   call bforce (dumm)
#ifdef MPI
	do j=JB1,JB2
	do i=IB1,IB2
#else
	  do i=1,m
	  do j=1,n
#endif
	  do ill=1,jjc(i,j)
           fqg(i,j) = dumm(i,j)  ! [W/m**2]
            dswr1(izet(i,j))=dumm(i,j)
	  enddo
	  enddo
	  enddo
#ifdef MPI
	do j=JB1,JB2
	do i=IB1,IB2
#else
	  do i=1,m
	  do j=1,n
#endif
	  do ill=1,jjc(i,j)
c	   if(frice(i,j).gt.0.0)then
	   if(frice(i,j).gt.epsis)then
	      einstrice(i,j) = fqg(i,j)
c	       einstrice(i,j) = fqg(i,j)*(1.0-0.65)/(1.0-0.07)
	   else
	      einstr(i,j) = fqg(i,j)
	   endif
	  enddo
	  enddo
	  enddo
c------------------------------------------------------------------------
c	 longwave radiation taken from ncep
c------------------------------------------------------------------------
	  call interpol2(dumm,erastr,weight,nrmess,ianzpkt,nera)
#ifdef MPI
	do j=JD1,JD2
	do i=ID1,ID2
#else
	do i=1,m
	do j=1,n
#endif
	    dumm(i,j)=dumm(i,j)+strkorr(i,j)
	enddo
	enddo
	   call bforce (dumm)
#ifdef MPI
	do j=JB1,JB2
	do i=IB1,IB2
#else
	  do i=1,m
	  do j=1,n
#endif
	  do ill=1,jjc(i,j)
	    fqr(izet(i,j)) =-dumm(i,j)
	  enddo
	  enddo
	  enddo

	endif !	 if (radinp) then

c--------------------------------------------------------------------------------------
c      cloudiness
c-------------------------------------------------------------------------------------

	 call interpol2(dumm,eracl,weight,nrmess,ianzpkt,nera)
	 call bforce (dumm)

	 dumm=dumm/100.	 !in data cloudiness in %
#ifdef MPI
	do j=JB1,JB2
	do i=IB1,IB2
#else
	do i=1,m
	do j=1,n
#endif
	do ill=1,jjc(i,j)
	  bew(izet(i,j)) = max(0.,dumm(i,j)) !there are negative values in NCEP dta
	enddo
	  dumm(i,j) = max(0.,dumm(i,j))*8.   ! used further for calculation if(.not.radinput)
	enddo
	enddo
#endif   !NCEP
#ifdef COASTDAT
       it_init=it_init+1
c---------------------------------------------------------------
c   read and sort metinput on mode grid using metind1 from
convert_metind
c-------------------------------------------------------------------
      if(myid.eq.0)print*,'read coastdat meteorology from CLM',it_init
      call nc_read(pathpc3//pathmet//metpost//fileuv10//metappend,
     &"U_10M",
     &dump(:,:),start=[1,1,1,it_init],count=[ncdlam,ncdphi,1,1])
      if(myid.eq.0)print*,'read u10'
       ixc=0
       do icn=1,ncdlam
       do icm=1,ncdphi
       ixc=ixc+1
       erau10(ixc)=dump(icn,icm)
       enddo
       enddo
c-------------------------------------------------------------
      call nc_read(pathpc3//pathmet//metpost//fileuv10//metappend,
     &"V_10M",
     &dump(:,:),start=[1,1,1,it_init],count=[ncdlam,ncdphi,1,1])
       ixc=0
       do icn=1,ncdlam
       do icm=1,ncdphi
       ixc=ixc+1
       erav10(ixc)=dump(icn,icm)
       enddo
       enddo
c-----------------------------------------------------------------------
         if(myid .eq. 0) print*,'interpolation of met. data'
c-----------------------------------------------------------------------

        do ib=1,nera
	  bwind(ib)=sqrt(erau10(ib)**2+erav10(ib)**2)
	enddo
        
c         if(ijahr.eq.1979.and.monat.eq.2)then
c         do ib=1,nera
c          erau10(ib)=0.8*erau10(ib)
c          erav10(ib)=0.8*erav10(ib)
c	  bwind(ib)=sqrt(erau10(ib)**2+erav10(ib)**2)
c	enddo
c        if(myid.eq.0)print*,'maxval wind', maxval(bwind,1)
cc        endif
c       endif
c         if(ijahr.eq.1987.and.monat.eq.1)then
c        do ib=1,nera
c          erau10(ib)=0.9*erau10(ib)
c          erav10(ib)=0.9*erav10(ib)
c	  bwind(ib)=sqrt(erau10(ib)**2+erav10(ib)**2)
c	enddo
c        if(myid.eq.0)print*,'maxval wind', maxval(bwind,1)
c       endif
c         if(ijahr.eq.1996.and.monat.eq.2)then
c        do ib=1,nera
c          erau10(ib)=0.9*erau10(ib)
c          erav10(ib)=0.9*erav10(ib)
c	  bwind(ib)=sqrt(erau10(ib)**2+erav10(ib)**2)
c	enddo
c        if(myid.eq.0)print*,'maxval wind', maxval(bwind,1)
c       endif
c         if(ijahr.eq.2010.and.monat.eq.2)then
c        do ib=1,nera
c          erau10(ib)=0.9*erau10(ib)
c          erav10(ib)=0.9*erav10(ib)
c	  bwind(ib)=sqrt(erau10(ib)**2+erav10(ib)**2)
c	enddo
c        if(myid.eq.0)print*,'maxval wind', maxval(bwind,1)
c       endif
       
 	call interpol2(bwindn,bwind,weight,nrmess,ianzpkt,nera)
	call bforce(bwindn)
	call interpol2(tyy   ,erav10,weight,nrmess,ianzpkt,nera)
	call bforce (tyy)
	call interpol2(txx   ,erau10,weight,nrmess,ianzpkt,nera)
	call bforce (txx)

c----------------------------------------------------------------
c     correction of wind from wind speed and direction after interpolation on model grid
c----------------------------------------------------------------
#ifdef MPI
	do inu=JB1,JB2
	do imu=IB1,IB2
#else
	do imu =1,m
	do inu=1,n
#endif
	  call uv2sr(txx(imu,inu),tyy(imu,inu),speedu,dirou,0,1)
	  diru(imu,inu)=dirou
	  call uv2sr(txxu,tyyu,bwindn(imu,inu),diru(imu,inu),1,1)
	  txx(imu,inu)=txxu
	  tyy(imu,inu)=tyyu
	enddo
	enddo
       
c       endif !itestxxx
c-------------------------------------------------------------
      call nc_read(pathpc3//pathmet//metpost//filetl//metappend,
     &"T_2M",
     &dump(:,:),start=[1,1,1,it_init],count=[ncdlam,ncdphi,1,1])
      if(myid.eq.0)print*,'read t_2m'
       ixc=0
       do icn=1,ncdlam
       do icm=1,ncdphi
       ixc=ixc+1
       eratl(ixc)=dump(icn,icm)
        enddo
        enddo
        call interpol2(dumm,eratl,weight,nrmess,ianzpkt,nera)
         call bforce (dumm)
#ifdef MPI
        do j=JB1,JB2
        do i=IB1,IB2
#else
        do i=1,m
        do j=1,n
#endif
        do ill=1,jjc(i,j)
          taair(izet(i,j)) = dumm(i,j)
        enddo
        enddo
        enddo
c       if(myid.eq.-1)then
c
c       file_outx='test_dump_tttt_mm.nc'
c       write(file_outx(11:14),'(i4.4)')icccount
c      write(file_outx(16:17),'(i2.2)')myid
c       if(myid.eq.0) print*,'write test ',file_outx 
c       call nc_create(file_outx,overwrite=.TRUE.
c     &,netcdf4=.TRUE.,author="Ute Daewel")
c      call nc_write_dim(file_outx,"lat",x=65.+59./60.,dx=-6./60.,nx=m,
c     &units="degrees_north")
c      call nc_write_dim(file_outx,"lon",x=-4.-10./60.,dx=10./60.,nx=n,
c     &units="degrees_east")
c      call nc_write(file_outx,"dumm",
c     &dumm(:,:),dim1="lat",dim2="lon",
c     &start=[1,1],count=[m,n])
c       endif
c-------------------------------------------------------------
      call nc_read(pathpc3//pathmet//metpost//filep//metappend,"PMSL",
     &dump(:,:),start=[1,1,it_init],count=[ncdlam,ncdphi,1])
       ixc=0
       do icn=1,ncdlam
       do icm=1,ncdphi
       ixc=ixc+1
       erap(ixc)=dump(icn,icm)
        enddo
        enddo
       call interpol2(pat,erap ,weight,nrmess,ianzpkt,nera) ! pressure 	
       call bforce (pat)
c-------------------------------------------------------------i

      call nc_read(pathpc3//pathmet//metpost//filecl//metappend,"CLCT",
     &dump(:,:),start=[1,1,it_init],count=[ncdlam,ncdphi,1])
      if(myid.eq.0)print*,'read clct'
       ixc=0
       do icn=1,ncdlam
       do icm=1,ncdphi
       ixc=ixc+1
       eracl(ixc)=dump(icn,icm)
        enddo
        enddo
       call interpol2(dumm,eracl ,weight,nrmess,ianzpkt,nera) ! pressure 	
       call bforce (dumm)
#ifdef MPI
        do j=JB1,JB2
        do i=IB1,IB2
#else
        do i=1,m
        do j=1,n
#endif
        do ill=1,jjc(i,j)
          bew(izet(i,j)) = max(0.,dumm(i,j)) !there are negative values
        enddo
        enddo
        enddo
c------------------------------------------------------------
      call nc_read(pathpc3//pathmet//metpost//filepri//metappend,
     &"TOT_PREC",
     &dump(:,:),start=[1,1,it_init],count=[ncdlam,ncdphi,1])
       ixc=0
       do icn=1,ncdlam
       do icm=1,ncdphi
       ixc=ixc+1
       erapri(ixc)=dump(icn,icm)
        enddo
        enddo
       call interpol2(dumm,erapri ,weight,nrmess,ianzpkt,nera) ! pressure 	
       call bforce (dumm)
       if(myid.eq.-1)then

       file_outx='test_dump_tttt_mm.nc'
       write(file_outx(11:14),'(i4.4)')icccount
      write(file_outx(16:17),'(i2.2)')myid
       if(myid.eq.0) print*,'write test ',file_outx 
       call nc_create(file_outx,overwrite=.TRUE.
     &,netcdf4=.TRUE.,author="Ute Daewel")
      call nc_write_dim(file_outx,"lat",x=65.+59./60.,dx=-6./60.,nx=m,
     &units="degrees_north")
      call nc_write_dim(file_outx,"lon",x=-4.-10./60.,dx=10./60.,nx=n,
     &units="degrees_east")
      call nc_write(file_outx,"dumm",
     &dumm(:,:),dim1="lat",dim2="lon",
     &start=[1,1],count=[m,n])
       endif
      if(myid.eq.0)print*,'read precip'
#ifdef MPI
        do j=JB1,JB2
        do i=IB1,IB2
#else
        do i=1,m
        do j=1,n
#endif
        do ill=1,jjc(i,j)
         precip(izet(i,j)) = dumm(i,j)/(real(imet*60*60)) !0.8 for CD2 correctionconvert to kg/m**2/sec
	 if (precip(izet(i,j)).lt.0.0) precip(izet(i,j))=0.
        enddo
        enddo
        enddo
c      print*,'read met p',maxval(precip),minval(precip)
c-------------------------------------------------------------
      call nc_read(pathpc3//pathmet//metpost//filetd//metappend,"TD_2M",
     &dump(:,:),start=[1,1,1,it_init],count=[ncdlam,ncdphi,1,1])
       ixc=0
       do icn=1,ncdlam
       do icm=1,ncdphi
       ixc=ixc+1
       eratp(ixc)=dump(icn,icm)
        enddo
        enddo
       call interpol2(dumm,eratp ,weight,nrmess,ianzpkt,nera) ! pressure 	
       call bforce (dumm)
#ifdef MPI
        do j=JB1,JB2
       do i=IB1,IB2
#else
        do i=1,m
        do j=1,n
#endif
        do ill=1,jjc(i,j)
          taup(izet(i,j)) = dumm(i,j) !there are negative values in NCEP
        enddo
        enddo
        enddo
c      if(myid.eq.0)print*,'read td'

c-------------------------------------------------------------
      call nc_read(pathpc3//pathmet//metpost//filesh//metappend,"QV_2M",
     &dump(:,:),start=[1,1,1,it_init],count=[ncdlam,ncdphi,1,1])
       print*,'read qv'
       ixc=0
       do icn=1,ncdlam
       do icm=1,ncdphi
       ixc=ixc+1
       eraqv(ixc)=dump(icn,icm)
        enddo
        enddo
       call interpol2(dumm,eraqv ,weight,nrmess,ianzpkt,nera) ! pressure 	
       call bforce (dumm)
       if(myid.eq.0)print*,'read shum', maxval(shum)
#ifdef MPI
        do j=JB1,JB2
        do i=IB1,IB2
#else
        do i=1,m
        do j=1,n
#endif
        do ill=1,jjc(i,j)
          shum(izet(i,j)) = dumm(i,j) !there are negative values in NCEP
        enddo
        enddo
        enddo
c#ifdef MPI
c	do j=JB1,JB2
c	do i=IB1,IB2
c#else
c	do i=1,m
c	do j=1,n
c#endif
c	do ill=1,jjc(i,j)
c	  call shumtaup(shum(izet(i,j)),pat(i,j),taup(izet(i,j)))
c	   taup(izet(i,j))=taup(izet(i,j))
c	enddo
c	enddo
c	enddo
c------------------------------------------------------------
      call nc_read(pathpc3//pathmet//metpost//filedssr//metappend,
     &"ASOB_S",
     &dump(:,:),start=[1,1,it_init],count=[ncdlam,ncdphi,1])
       ixc=0
       do icn=1,ncdlam
       do icm=1,ncdphi
       ixc=ixc+1
       erassr(ixc)=dump(icn,icm)
        enddo
        enddo
       call interpol2(dumm,erassr ,weight,nrmess,ianzpkt,nera) ! pressure 	
       call bforce (dumm)
#ifdef MPI
        do j=JB1,JB2
        do i=IB1,IB2
#else
          do i=1,m
          do j=1,n
#endif
          do ill=1,jjc(i,j)
           fqg(i,j) = dumm(i,j)  ! [W/m**2]
           dswr1(izet(i,j))=dumm(i,j)
          enddo
          enddo
          enddo
#ifdef MPI
        do j=JB1,JB2
        do i=IB1,IB2
#else
          do i=1,m
          do j=1,n
#endif
          do ill=1,jjc(i,j)
           if(frice(i,j).gt.epsis)then
              einstrice(i,j) = fqg(i,j)
           else
              einstr(i,j) = fqg(i,j)
           endif
          enddo
          enddo
          enddo

c--------------------------------------------------------------
#endif !COASTDAT
c        print*,'metinp txx tyy',maxval(txx),maxval(tyy)
c        print*,'metinp pat taair',maxval(pat),maxval(taair)
c        print*,'metinp shum taup',maxval(shum),maxval(taup)
c        print*,'metinp bew pre',maxval(bew),maxval(precip)
c        print*,'metinp fqr fqg',maxval(fqr),maxval(fqg)
c        print*,'metinp fqr fqg min',minval(fqr),minval(fqg)

c---------------------------------------------------------------------
c-----------------------------------------------------------------------
	 print*,'  set met values on slab'
c-----------------------------------------------------------------------
#ifdef MPI
        do j=JB1,JB2
        do i=IB1,IB2
#else
        do i=1,m
        do j=1,n
#endif
           pat(i,j) = pat(i,j)*0.001
           do ii=1,jjc(i,j)
            taair(izet(i,j)) = taair(izet(i,j))-273.15
           taup(izet(i,j))  = taup(izet(i,j))-273.15
#ifdef COASTDAT
           bew(izet(i,j))   = bew(izet(i,j))*8.
#else
           bew(izet(i,j))   = bew(izet(i,j))*8.
#endif
           paac2(izet(i,j)) = paac(izet(i,j))
           txac2(izet(i,j)) = txac(izet(i,j))
           tyac2(izet(i,j)) = tyac(izet(i,j))
           txac(izet(i,j))  = txx(i,j)
           tyac(izet(i,j))  = tyy(i,j)
           paac(izet(i,j))  = pat(i,j)
c           q10a(izet(i,j))  = q10(izet(i,j))	 !was used with aaa and bbb, now unsed
c           tl10a(izet(i,j)) = tl10(izet(i,j))
         enddo
       enddo
       enddo

c-----------------------------------------------------------------------
	 print*,'  compute heat fluxes'
c-----------------------------------------------------------------------
	 call fluxes(izet,taair,shum,bew,paac,txac,tyac,
     *	 uc,vc,tec,qwt,id3sur,wspeed,q10,tl10,
     *	 fqs,fql,cddu,evap)
c#ifdef NCEP
#ifdef MPI
	do j=1,khorl
	i=lzet(j)
#else
	do i=1,khor
#endif
c-----------------------------------------------------------------------
c	 long wave backward radiation (rbackrad)

c	 long wave radiation after Idso and Jackson (1969)
	 rbackrad = rbacknn(tl10(i),tec(id3sur(i)),bew(i))

c	 long wave radiation after Maykut (1986)
c	 rbackrad = rback(tl10a(i),tec(id3sur(i),cn1,sigem)

c	 long wave radiation after  Bignami et al. (1995)
c	 rbackrad =  rbackn((tl10a(i),tec(id3sur(i),cn1,eaa,sigem)
c-----------------------------------------------------------------------
         fqr(i) = rbackrad   !air - STT

	enddo
c#endif
*=========output========================================================

#ifndef MPI
	 write(121) txac,tyac,wspeed
        write(121) taair,dswr1,tl10
        write(121) precip
#else
	 DO i=1,khorl
	 k=lzet(i)
       WRITE(122)txac(k),tyac(k),wspeed(k),taair(k),
     &   dswr1(k),tl10(k),precip(k)
       ENDDO
#endif
*=========output========================================================


c#ifdef MPI
c      if(myid.eq.0)then
c      write(*,*)'Before W Wind 30% increase:'
c      write(*,*)'evaU==',maxval(erau10),minval(erau10)
c      write(*,*)'evaV== ',maxval(erav10),minval(erav10)
c      write(*,*)'txac==',maxval(txac),minval(txac)
c      write(*,*)'tyac== ',maxval(tyac),minval(tyac)
c      write(*,*)'After Wind 30% increase:'
c      write(*,*)'txac==',maxval(txac),minval(txac)
c      write(*,*)'tyac== ',maxval(tyac),minval(tyac)
c       endif
c#endif
c--------------------------------------------------------------
c		wind stress for ice
c--------------------------------------------------------------
#ifdef MPI
	do j=JB1,JB2
	do i=IB1,IB2
#else
	do i=1,m
	do j=1,n
#endif
           do ii=1,jjc(i,j)
	 
           wgesch(i,j) = sqrt(txx(i,j)**2+tyy(i,j)**2)
	   cd = 1.28*(0.8+0.065*wgesch(i,j))
	   cd = cd*1.e-3
	   txx(i,j) = cd*wgesch(i,j)*txx(i,j)
	   tyy(i,j) = cd*wgesch(i,j)*tyy(i,j)
	enddo
	enddo
	enddo


	if(nstart.eq.0) then
	  nstart = 1
c----------------------------------------------------------------------------
c   10m reference level for temperature and humidty after subroutine fluxes
c---------------------------------------------------------------------------
#ifdef MPI
	  do i=1,khorl
	  ijij=lzet(i)
#else
	  do ijij=1,khor
#endif
	    q10a(ijij)	= q10(ijij)
	    tl10a(ijij) = tl10(ijij)
	  enddo
	endif  !if(nstart.eq.0) then

      else
	if((lmon.eq.mona).and.(nday.lt.itanf)) goto 1000
      endif   ! metinread

c#ifdef MPI
c      if(progno)then
c      call MPI_Waitall(nreq1,ir_r_b, status, ierr)
c      INUM1=0
c      INUM2=0
c      do I=3,nbio
c      call VER_UNPACK(Tc(:,I),ndrei,dtemp,iel,2,2,1,3,lazc)
c      enddo
c      INUM3=0
c      INUM4=0
c      do I=3,nbio
c      call HOR_UNPACK(Tc(:,I),ndrei,dtemp,iel,2,2,nbhlayer0,nthlayer0,
c     &idepb0,idept0,icdepb0,icdept0,1,3,iwet)
c      enddo
c      endif
c#endif

c-----------------------------------------------------------------------
c     if temporal interpolation of wind and pressure is required in between
c     met time steps, it is done here
c-----------------------------------------------------------------------

      atime = 1.
      btime = 0.
#ifdef MPI
	  do j=ja1,ja2
	  do i=ia1,ia2
	  do k=1,jjc(i,j)
	  l=izet(I,J)
#else
      do l=1,khor
#endif
	 pac(l) = atime*paac(l)+btime*paac2(l)
	 txc(l) = atime*txac(l)+btime*txac2(l)
	 tyc(l) = atime*tyac(l)+btime*tyac2(l)
      enddo
#ifdef MPI
       enddo
       enddo
#endif


	if (.not. radinp) then
c-----------------------------------------------------------------------
c	print*,'  compute longwave radiation'
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c      print*,'	 compute shortwave radiation, subroutine feld'
c-----------------------------------------------------------------------

	call feld(dumm,einstr,einstrice, !===> incoming radiation [W/m**2]
     *mjar,time,
     *itest,jjc)



      endif  !	if (.not. radinp) then


#ifdef MPI
	  do l=1,khorl
	  i=lzet(l)
#else
	  do i= 1,khor
#endif
	  qw(i) = qwt(i)-fqr(i)
	  enddo

c-----------------------------------------------------------------------
c      print*,'	 compute water ice stress'
c-----------------------------------------------------------------------
      call estate
      call stressw(iindex,izet,jjc,z2d,m,n,dz,ilo,txx,tyy)

      if(withice.and.(mod(ip+i_ice-1,i_ice).eq.0)) then
c-----------------------------------------------------------------------
c     sea ice dynamics
      IF(MYID .EQ. 0) print*,'	compute	 sea ice dynamics '
c-----------------------------------------------------------------------
#ifdef MPI
	do j=JB1,JB2
	do i=IB1,IB2
#else
	do j=1,n
	do i=1,m
#endif
	   fricei(i,j,1) = frice(i,j)
	   hisi(i,j,1)	 = his(i,j)
	   hisir(i,j,1)	 = hisr(i,j)
	enddo
	enddo
	call icemod
	call icevel

#ifdef MPI
	do j=J1,J2
	do i=I1,I2
#else
	do j=1,n
	do i=1,m
#endif
	   frice(i,j) = fricei(i,j,1)
	   his(i,j)   = hisi(i,j,1)
	   hisr(i,j)  = hisir(i,j,1)
	enddo
	enddo

	endif !if(withice.and.(mod(ip+i_ice-1,i_ice).eq.0)) then

c-----------------------------------------------------------------------
c     dummy arrays set to zero
c-----------------------------------------------------------------------
      call settoz(nlar,ngro,large,gross)


c-----------------------------------------------------------------------
c      print*,' call  motmit'
c     subroutine konv:
c	     preparation for solving the non-linear terms in the equation
c	     of motion; sheme used here is the Arakawa J7 Operator, energy
c	     and enstrophy conserving sheme
c     subroutine motmit:
c	     solves the explicit part of the equation of motion (u,v)
c	     implicit terms are barotropic pressure gradients and
c			  vertical exchange coefficients
c     subroutine uvrand:
c	     applies boundary conditions for the transport field
c-----------------------------------------------------------------------
      IF(MYID .EQ. 0) print*,' motmit, konv, uvrand '
      
      call deco1d2d(z2d,zac,iwet,indend,1.)
      call konv (dt,pd,dx,dy,z2d,dz1,dzbod,
     *		 ukonv,vkonv,uc,vc)
      call motmit (htmit,vtmit,vtmic,slip,botdic,ukonv,vkonv,
     *	ndrei,ht,pah_fac)
#ifndef MPI
      call uvrand(izet,iindex,time)
#else
      call uvrand_par(uc,izet,iindex,0,1)
      call uvrand_par(vc,izet,iindex,1,0)
      call uvrand_par(wc,izet,iindex,0,0)  !Ute 24.10.
      call uvrand_par(avc,izet,iindex,0,0)
#endif
c-----------------------------------------------------------------------
c     dummy arrays set to zero
c-----------------------------------------------------------------------
      call settoz(nlar,ngro,large,gross)

c-----------------------------------------------------------------------
c      print*,'	 compute sor'
c      subroutine sorcof:
c	   prepares coefficient matrix for iteration
c      subroutine zetadb:
c	   boundary values in zeta
c      subroutine sor: solving the implicit equation for the
c	    surface elevation, red-black scheme is used
c-----------------------------------------------------------------------
      call sorcof(htmit,htmitq)
      call zetadb(nday,iminu,iindex,nzac,time)
      call sor(htmit,nosor,eps,ip)

c-----------------------------------------------------------------------
c     dummy arrays set to zero
c-----------------------------------------------------------------------
      call settoz(nlar,ngro,large,gross)
c-----------------------------------------------------------------------
c      print*,'	 compute drux '
c      subroutine druxav: adds time averaged barotropic pressure gradients
c	   to the equation of motion (implicit terms) and adds implicit
c	   vertical turbulent exchange (time dependent av).
c	   final solving of equation of motion
c	   turbulence closure with analytical k-e approach:
c	   after Kochergin, 1987 and Pohlmann (1996) and as updated by Schrum (1997);
c-----------------------------------------------------------------------
      call druxav(vtmit,avmax,cdeav,szahl,ukonv,vkonv,ndrei)
c-----------------------------------------------------------------------
c      consideration of river run-off for sea surface elevation
c      and consideration of nutrient loads
c-----------------------------------------------------------------------
	IF(river_on)THEN

      do ii=1,nanz
      zriv(indm(ii),indn(ii))=0.
      enddo

      do ii=1,nanz
      zriv(indm(ii),indn(ii)) = zriv(indm(ii),indn(ii))
     &+deltat*eintr(ii)/(dl*dln(indm(ii)))
      enddo
      IF(MYID .EQ. 0) print*,' end river on'

      ENDIF !IF(river_on)THEN

c-----------------------------------------------------------------------
c	Influence of precipitation and evaporation on surface elevation
c-----------------------------------------------------------------------

	evafac	= deltat/1000.0
c      write(*,*)'3 evapMAIN=',maxval(evap),minval(evap)
#ifdef MPI
      do j=J1,J2
	  lwa=isornr1(j)
	  lwe=isorsr1(j)
#else
      do j=1,n
	  lwa=isornr(j)
	  lwe=isorsr(j)
#endif
      do lw=lwa,lwe
	 i = iwet(lw)
	   zpe(i,j)=(precip(lw)-evap(lw) )*evafac   ![m]
      enddo
      enddo
c       print*,'test precip evap',maxval(precip),maxval(evap)
c---------------------------------------------------------
c	set zinc1 aus konti gleichung statt aus implicit
c---------------------------------------------------------
      lwe = 0
#ifdef MPI
      do j=J1,J2
	lwa = lb0(j)
	lwe = le0(j)
#else
      do j=1,n
	lwa = lwe+1
	lwe = indend(j)
#endif
	do lw=lwa,lwe
	  i = iwet(lw)
	  zac(lw) =zac(lw)+ z(i,j)

cc	     +zpe(i,j)+zriv(i,j)   !! CS_rivers
	enddo
      enddo
c-------------------------------------------------------------------------
c     check and preventing of dry grid boxes for extreme situations
c-------------------------------------------------------------------------
cc	do j=1,-1
#ifdef MPI
	do j=J1,J2
	  lwa=isornr1(j)
	  lwe=isorsr1(j)
#else
	do j=1,n
	  lwa=isornr(j)
	  lwe=isorsr(j)
#endif
      do lw=lwa,lwe
	  i = iwet(lw);lump = lazc(lw)
	    if(lump.eq.1.)then
	    ddd=real(ldep(lw))
	    if(ddd+zac(lw).lt.0.5)zac(lw)=0.5-ddd
	    else
	    if(dz(1)+zac(lw).lt.0.5)zac(lw)=0.5-dz(1)
	    endif
	enddo
      enddo
c--------------------------------------------
c     output statistics on screen
c--------------------------------------------
#ifdef MPI
      call uvrand_par(uc,izet,iindex,0,1)
      call uvrand_par(vc,izet,iindex,1,0)
      call uvrand_par(wc,izet,iindex,0,0)
      call uvrand_par(avc,izet,iindex,0,0)
      call INTI_COMM_VAR

      call com_vertical(uc,ndrei,dtemp,iel,2,2,1,3,0,lazc)
      call com_vertical(vc,ndrei,dtemp,iel,2,2,1,3,0,lazc)
      call com_vertical(zac,khor,dtemp,iel,2,2,1,2,1,lazc)
      call com_horizontal(uc,ndrei,dtemp,iel,2,2,1,3,0,icdepb2,
     *idepb2,icdept2
     *,idept2,ncsendb2,nsendb2,ncsendt2,nsendt2,nbhlayer2,nthlayer2)
      call com_horizontal(vc,ndrei,dtemp,iel,2,2,1,3,0,icdepb2,
     *idepb2,icdept2
     *,idept2,ncsendb2,nsendb2,ncsendt2,nsendt2,nbhlayer2,nthlayer2)
      call com_horizontal(zac,khor,dtemp,iel,2,2,1,2,1,icdepb2,
     *idepb2,icdept2
     *,idept2,ncsendb2,nsendb2,ncsendt2,nsendt2,nbhlayer2,nthlayer2)
#endif

c-----------------------------------------------------------------------
c     boundary conditions for transports
c-----------------------------------------------------------------------
c set necessary equal boundary lines for u,v,w,avc:
#ifndef MPI
      call uvrand (izet,iindex,time)
#else
      call MPI_Waitall(nreq1, ir_r_b, status, ierr)
      INUM1=0
      INUM2=0
      call VER_UNPACK(uc,ndrei,dtemp,iel,2,2,1,3,lazc)
      call VER_UNPACK(vc,ndrei,dtemp,iel,2,2,1,3,lazc)
      call VER_UNPACK(zac,khor,dtemp,iel,2,2,1,2,lazc)
      INUM3=0
      INUM4=0
      call HOR_UNPACK(uc,ndrei,dtemp,iel,2,2,nbhlayer2,nthlayer2,
     &idepb2,idept2,icdepb2,icdept2,1,3,iwet)
      call HOR_UNPACK(vc,ndrei,dtemp,iel,2,2,nbhlayer2,nthlayer2,
     &idepb2,idept2,icdepb2,icdept2,1,3,iwet)
      call HOR_UNPACK(zac,khor,dtemp,iel,2,2,nbhlayer2,nthlayer2,
     &idepb2,idept2,icdepb2,icdept2,1,2,iwet)
#endif
c-----------------------------------------------------------------------
c     dummy arrays set to zero
c-----------------------------------------------------------------------
      call settoz(nlar,ngro,large,gross)
c-----------------------------------------------------------------------
c     subroutine konti: equation of continuity => compute w
      if (myid.eq.0)print*,' konti: equation of continuity => compute w'
      if(myid.eq.0)print*,'before konti',maxval(uc),minval(uc)
c-----------------------------------------------------------------------
      call konti(ht,dzbod,wtest,wsur)
c set necessary equal boundary lines for u,v,w,avc:
#ifndef MPI
      call uvrand(izet,iindex,time)
#else
      call uvrand_par(uc,izet,iindex,0,1)
      call uvrand_par(vc,izet,iindex,1,0)
      call uvrand_par(wc,izet,iindex,0,0)
      call uvrand_par(avc,izet,iindex,0,0)
      call INTI_COMM_VAR
      call com_vertical(wc,ndrei,dtemp,iel,1,1,1,3,0,lazc)
      call com_vertical(avc,ndrei,dtemp,iel,1,1,1,3,1,lazc)
      call com_horizontal(wc,ndrei,dtemp,iel,1,1,1,3,0,icdepb1,
     *idepb1,icdept1
     *,idept1,ncsendb1,nsendb1,ncsendt1,nsendt1,nbhlayer1,nthlayer1)
      call com_horizontal(avc,ndrei,dtemp,iel,1,1,1,3,1,icdepb1,
     *idepb1,icdept1
     *,idept1,ncsendb1,nsendb1,ncsendt1,nsendt1,nbhlayer1,nthlayer1)

      call MPI_Waitall(nreq1, ir_r_b, status, ierr)

      INUM1=0
      INUM2=0

      call VER_UNPACK(wc,ndrei,dtemp,iel,1,1,1,3,lazc)
      call VER_UNPACK(avc,ndrei,dtemp,iel,1,1,1,3,lazc)

      INUM3=0
      INUM4=0
      call HOR_UNPACK(wc,ndrei,dtemp,iel,1,1,nbhlayer1,nthlayer1,
     &idepb1,idept1,icdepb1,icdept1,1,3,iwet)
      call HOR_UNPACK(avc,ndrei,dtemp,iel,1,1,nbhlayer1,nthlayer1,
     &idepb1,idept1,icdepb1,icdept1,1,3,iwet)
#endif
c-----------------------------------------------------------------------
c     dummy arrays set to zero
c-----------------------------------------------------------------------
      call settoz(nlar,ngro,large,gross)
      if(progno) then
c-----------------------------------------------------------------------
c	prognostic calculation of temperature and salinity
c	advection and vertical diffusion
c-----------------------------------------------------------------------

c---------------------------------------------------------
c	set zinc1 aus konti gleichung statt aus implicit
c---------------------------------------------------------

      lwe = 0
#ifdef MPI
      do j=J1,J2
	lwa = lb0(j)
	lwe = le0(j)
#else
      do j=1,n
	lwa = lwe+1
	lwe = indend(j)
#endif
	do lw=lwa,lwe
	  i = iwet(lw)
	  zinc1(i,j)=wsur(lw)*dble(deltat)
c     &		 +dble(zpe(i,j))+dble(zriv(i,j))
c----------------------------------------------
c	set open boundary values for zinc1
c----------------------------------------------
      if(i.lt.68.and.j.lt.60)then
      zinc1(i,j)=zinc1(68,j)
      elseif(i.gt.150.and.j.lt.18)then
      zinc1(i,j)=zinc1(i,18)
      endif
c-----------------------------------------------
      zinc(i,j)=dble(z(i,j))
     &+dble(zpe(i,j))+dble(zriv(i,j))
       enddo
       enddo
c-----------------------------------------------------------------------
c     subroutine strom3
c	advection and vertical diffusion of T,S and state variables of
c	radiation into lower water layers is considered in strom3 (from 2nd layer)
c-----------------------------------------------------------------------
     
       call stromTVD(zinc1,
     &  vtmit,vtmic,avmax,ftac,szahl) !54 file outtest
c       call stromTVD(zinc1,
c     &  vtmit,vtmic,avmax,ftac,szahl,
c     &  einstr,nday,time,54,wgesch) !54 file outtest
       call estate
c-----------------------------------------------------
c  calculate biology 
c---------------------------------------------------
#ifdef _FABM_
!KB hardcoding npel to 5 !!!
      call do_hamsom_fabm(ndrei,5)
!KB      stop 'do_hamsom_fabm()'
#endif

c-----par rausschreiben--------------------------------------
c	if(rlicht(1).ne.0.)then
c	l1count=l1count+1
c	write(l2count(1:4),'(i4.4)')l1count
c	nlicht=104
c      open(nlicht,
c     & file=pathpc2//'licht_'//l2count//'.dat',form='formatted',
c     & status='unknown')
c
c	 print*,'licht2 ',rlicht,nlicht
c	write(nlicht,'(20(1x,e12.4))') (rlicht(il),il=1,20)
c	close(nlicht)
c	endif

c---------------------------------------------------------
c      add rivers and e-p to zac
c---------------------------------------------------------
      lwe = 0
#ifdef MPI
      do j=J1,J2
	lwa = lb0(j)
	lwe = le0(j)
#else
      do j=1,n
	lwa = lwe+1
	lwe = indend(j)
#endif
	do lw=lwa,lwe
	  i = iwet(lw)
	  zac(lw) =zac(lw) +zpe(i,j)+zriv(i,j)
	enddo
      enddo
c-------------------------------------------------------------------------
c     check and preventing of dry grid boxes for extreme situations
c-------------------------------------------------------------------------
cc	do j=1,-1
#ifdef MPI
	do j=J1,J2
	  lwa=isornr1(j)
	  lwe=isorsr1(j)
#else
	do j=1,n
	  lwa=isornr(j)
	  lwe=isorsr(j)
#endif
      do lw=lwa,lwe
       i = iwet(lw);lump = lazc(lw)
       if(lump.eq.1.)then
       ddd=real(ldep(lw))
       if((ddd+zac(lw)).lt.0.5)zac(lw)=0.5-ddd
        else
        if((dz(1)+zac(lw)).lt.0.5)zac(lw)=0.5-dz(1)
        endif
      enddo
      enddo

*=============================================================================
*---------Fluxes due to river and preci-evap [X*m*time_step]
      nwet=1
#ifdef	MPI
      do k=1,khorl
       lw=lzet(k)
#else
      do lw=1,khor
#endif
	nwet=id3sur(lw)
	i=iwet(lw);j=jwet(lw);
*T
	Tc_flux_r(lw,1)=tec(nwet)*zriv(i,j)	  ![DegC*m*time step]
	Tc_flux_e(lw,1)=tec(nwet)*zpe(i,j)	  ![DegC*m*time step]
*S
	Tc_flux_r(lw,2)=sac(nwet)*zriv(i,j)	  ![DegC*m*time step]
	Tc_flux_e(lw,2)=sac(nwet)*zpe(i,j)	  ![DegC*m*time step]
*Tc
c      do ibio=3,nbio
c	Tc_flux_r(lw,ibio)=zriv_bio(ibio,i,j)	 ![mmolX*m/time step]
c        if(i.eq.indm(153).and.j.eq.indn(153).and.ibio.eq.ino3)then
c       lwx=lw
c       endif

c	Tc_flux_e(lw,ibio)=0.			       ![mmolX*m/time step]
c	if (ibio.eq.ino3)then
c	Tc_flux_e(lw,ibio)=wetNH4(i,j)		![mmolX*m/time step]
c	elseif (ibio.eq.inh4)then
c	Tc_flux_e(lw,ibio)=wetNO3(i,j)		![mmolX*m/time step]
c	else
c	Tc_flux_e(lw,ibio)=0.			       ![mmolX*m/time step]
c	endif

	enddo

#ifdef MPI
       do k=1,khorl
       lw=lzet(k)
       do ibio=1,nbio
       Tc_flux_all(lw,ibio)=Tc_flux_r(lw,ibio)+Tc_flux_e(lw,ibio)
c       Tc_flux_all(lw,ibio)=Tc_flux_e(lw,ibio) !test without rivers
       enddo
       enddo
#else
	Tc_flux_all=Tc_flux_r+Tc_flux_e
#endif
*=============================================================================
*-----------------------------------------------------------------------
*  New salinity due to 1) the fresh voulume flux(d_vol)
*		   and 2) conservation of salt mass due to density changes(d_rho)
*-----------------------------------------------------------------------
c      do ibio=3,nbio
c      if(myid.eq.0)then
c       print*,'test riv ',ibio,
c     &  maxval(Tc(id3sur,ibio)),minval(Tc(id3sur,ibio)) 
c       endif
c      enddo
      foult=0
#ifdef MPI
      do j=j1,j2
      lwa=isornr1(j)
      lwe=isorsr1(j)
#else
      do j=1,n
      lwa=isornr(j)
      lwe=isorsr(j)
#endif
      do lw=lwa,lwe
      vol_new=pd(1)
      if(lazc(lw).eq.1)vol_new=ldep(lw)
      i = iwet(lw)

      vol_new=vol_new+zac(lw)
      dzet=zpe(i,j)+zriv(i,j)    !the total volume
      sfluss=0.0                 !salinity of the total volume
      rhofri=sigma(sfluss,tec(id3sur(lw)),pres(1))-refrho(1)+rhoq(1)
      rho_old=stc(id3sur(lw))+rhoq(1)
      rho_new=(rho_old*(vol_new-dzet) +rhofri*dzet)/(vol_new-dzet)
      d_rho=rho_old/rho_new
*New salinity:
      sac(id3sur(lw))=sac(id3sur(lw))*d_rho
c----------new nutrient concentrations ----------------------------------

c      do ibio=3,nbio
c      Tc(id3sur(lw),ibio)=Tc(id3sur(lw),ibio)*(vol_new-dzet)/(vol_new)+
c     &Tc_flux_all(lw,ibio)/(vol_new)
c      enddo
      enddo
      enddo
      
      if(foult.eq.1) then
      write(*,*)'vol_new.lt.0.... see file=98'
      pause
      endif
*--------End of *New salinity due -----------------


      call tsrneu(iindex,izet,time,dz,delsal)  !set equal boundary lines for T,S
      call estate
c	print*,'computing strom done'
      IF(MYID .EQ. 0) print*,'tsrneu & estate'

c--------------------------------------------------------------
c     check stability criteria
c     -->> remember all loops over the 3-d grid ndrei grid points
c	   in the time loop
c	   are time consuming and should be used with care
c--------------------------------------------------------------
c      goto 11199
#ifdef MPIP
	if(mod(ip-1,1).eq.0) then

      lwe = 0
      nwet=0
      umax  = -9999999.9
      vmax  = -9999999.9
      wmax  = -9999999.9
#ifdef MPI
	do J=J1,J2
	lwa = lb0(J)
	lwe = le0(J)
#else
	do j=1,n
	lwa = lwe+1
	lwe = indend(j)
#endif
	do lw=lwa,lwe
	  i = iwet(lw)
	    lump = lazc(lw)
       nwet = indwet(lw)
	do k=1,lump

      nwet = nwet+1
      dzz=  pd(k)
      if(k.eq.1)   dzz= (pd(k)+zac(lw))
      if(k.eq.lump)then
      dzz=float( ldep(lw) )
      if(k.eq.1)dzz=float( ldep(lw) )+zac(lw)
      endif
      u22t(nwet)=uc(nwet)/dzz
      v22t(nwet)=vc(nwet)/dzz
      w22t(nwet)=wc(nwet)


	if(abs(u22t(nwet)).gt.umax)then
	umax=abs(u22t(nwet))
	iumax=i
	jumax=j
	kumax=k
	numax=nwet
	zumax=dzz
	lumpu=lump
	endif

	if(abs(v22t(nwet)).gt.vmax)then
	vmax=abs(v22t(nwet))
	ivmax=i
	jvmax=j
	kvmax=k
	nvmax=nwet
	zvmax=dzz
	lumpv=lump
	endif

	if(abs(w22t(nwet)).gt.wmax)then
	wmax=abs(w22t(nwet))
	iwmax=i
	jwmax=j
	kwmax=k
	nwmax=nwet
	zwmax=dzz
	lumpW=lump
	 endif

	enddo !k=1,lump
      enddo
      enddo
c      print*,'test4',umax,uminm,vmax,vmin,wmax,wmin
c	print*,'computing strom done'
#ifdef MPI
      CALL MPI_ALLGATHER(umax,1,MPI_REAL,rumax,1,MPI_REAL,
     &MPI_COMM_WORLD,IERR)
      CALL MPI_ALLGATHER(vmax,1,MPI_REAL,rvmax,1,MPI_REAL,
     &MPI_COMM_WORLD,IERR)
      CALL MPI_ALLGATHER(wmax,1,MPI_REAL,rwmax,1,MPI_REAL,
     &MPI_COMM_WORLD,IERR)
      umax  = -9999999.9
      vmax  = -9999999.9
      wmax  = -9999999.9
      do I=1,nPh*nPv
      if(rumax(I).gt.umax)then
      umax=rumax(I)
      itmax1=I-1
      endif
      if(rvmax(I).gt.vmax)then
      vmax=rvmax(I)
      itmax2=I-1
      endif
      if(rwmax(I).gt.wmax)then
      wmax=rwmax(I)
      itmax3=I-1
      endif
      enddo
      If(umax.ne.0. .AND. myid .eq. ITMAX1 ) write(*,155)
     &'Umax:i,j,k',u22t(numax),iumax,jumax,kumax,zumax,'lump=',lumpu,
     &					       dl/(umax*min_time_step)
       if(vmax.ne.0. .AND. myid .eq. ITMAX2)write(*,155)
     &'Vmax:i,j,k',v22t(nvmax),ivmax,jvmax,kvmax,zvmax,'lump=',lumpv,
     &					       dl/(vmax*min_time_step)
       if(wmax.ne.0. .AND. myid .eq. ITMAX3)write(*,155)
     &'Wmax:i,j,k',w22t(nwmax),iwmax,jwmax,kwmax,zwmax,'lump=',lumpw,
     &					   zwmax/(wmax*min_time_step)
#else
       if(umax.ne.0.)write(*,155)
     &'Umax:i,j,k',u22t(numax),iumax,jumax,kumax,zumax,'lump=',lumpu,
     &					       dl/(umax*min_time_step)
      if(vmax.ne.0.)write(*,155)
     &'Vmax:i,j,k',v22t(nvmax),ivmax,jvmax,kvmax,zvmax,'lump=',lumpv,
     &					       dl/(vmax*min_time_step)
	if(wmax.ne.0.)write(*,155)
     &'Wmax:i,j,k',w22t(nwmax),iwmax,jwmax,kwmax,zwmax,'lump=',lumpw,
     &					   zwmax/(wmax*min_time_step)
#endif
155	format(a10,f10.4,3(1x,i3),1x,f10.5,1x,a5,i3,1x,f10.5)
c      print*,'test4c'

      if(myid .EQ. 0)write(*,*)'========T&S==============='
*T
      lwe = 0
      nwet=0
      umax  = -9999999.9
      umin  =  9999999.9
#ifdef MPI
	do J=J1,J2
	lwa = lb0(j)
	lwe = le0(j)
#else
	do j=1,n
	lwa = lwe+1
	lwe = indend(j)
#endif
	do lw=lwa,lwe
	  i = iwet(lw)
	    lump = lazc(lw)
	nwet=indwet(lw)
	do k=1,lump
	nwet = nwet+1
		      dzz=  pd(k)
	 if(k.eq.1)   dzz= (pd(k)+zac(lw))
	 if(k.eq.lump)dzz=float( ldep(lw) )

	if(abs(tec(nwet)).gt.umax)then
	     umax=abs(tec(nwet))
	iumax=i
	jumax=j
	kumax=k
c	 zumax=dzz
	zumax=zac(lw)
	    endif
	if(abs(tec(nwet)).le.umin)then
	     umin=abs(tec(nwet))
	iumin=i
	jumin=j
	kumin=k
	zumin=dzz
	    endif
      enddo
      enddo
      enddo
#ifdef MPI
      CALL MPI_ALLGATHER(umax,1,MPI_REAL,rumax,1,MPI_REAL,
     &MPI_COMM_WORLD,IERR)
      CALL MPI_ALLGATHER(umin,1,MPI_REAL,rvmax,1,MPI_REAL,
     &MPI_COMM_WORLD,IERR)
      umax  = -9999999.9
      umin  = 9999999.9
      do I=1,nPh*nPv
      if(rumax(I).gt.umax)then
      umax=rumax(I)
      itmax1=I-1
      endif
      if(rvmax(I).le.umin)then
      umin=rvmax(I)
      itmax2=I-1
      endif
      enddo
	IF(myid .EQ. itmax1) write(*,'(a14,f10.4,3(1x,i3),1(1x,f10.5))')
     & '    Tmax:i,j,k',umax,iumax,jumax,kumax,zumax
	IF(myid .EQ. itmax2) write(*,'(a14,f10.4,3(1x,i3),1(1x,f10.5))')
     & '    Tmin:i,j,k',umin,iumin,jumin,kumin,zumin
#else
	  write(*,'(a14,f10.4,3(1x,i3),1(1x,f10.5))')
     & '    Tmax:i,j,k',umax,iumax,jumax,kumax,zumax
	  write(*,'(a14,f10.4,3(1x,i3),1(1x,f10.5))')
     & '    Tmin:i,j,k',umin,iumin,jumin,kumin,zumin
#endif
      lwe = 0
      nwet=0
      umax  = -9999999.9
      umin  =  9999999.9
#ifdef MPI
	do J=J1,J2
	lwa = lb0(J)
	lwe = le0(J)
#else
	do j=1,n
	lwa = lwe+1
	lwe = indend(j)
#endif
	do lw=lwa,lwe
	  i = iwet(lw)
	    lump = lazc(lw)
	nwet=indwet(lw)
	do k=1,lump
	nwet = nwet+1
		      dzz=  pd(k)
	 if(k.eq.1)   dzz= (pd(k)+zac(lw))
	 if(k.eq.lump)dzz=float( ldep(lw) )

	if(abs(sac(nwet)).gt.umax)then
	     umax=abs(sac(nwet))
	iumax=i
	jumax=j
	kumax=k
	zumax=dzz
	    endif
	if(abs(sac(nwet)).le.umin)then
	     umin=abs(sac(nwet))
	iumin=i
	jumin=j
	kumin=k
	zumin=dzz
	    endif
      enddo
      enddo
      enddo
#ifdef MPI
      CALL MPI_ALLGATHER(umax,1,MPI_REAL,rumax,1,MPI_REAL,
     &MPI_COMM_WORLD,IERR)
      CALL MPI_ALLGATHER(umin,1,MPI_REAL,rvmax,1,MPI_REAL,
     &MPI_COMM_WORLD,IERR)
      umax  = -9999999.9
      umin  = 9999999.9
      do I=1,nPh*nPv
      if(rumax(I).gt.umax)then
      umax=rumax(I)
      itmax1=I-1
      endif
      if(rvmax(I).le.umin)then
      umin=rvmax(I)
      itmax2=I-1
      endif
      enddo
	IF(myid .EQ.itmax1)  write(*,'(a14,f10.4,3(1x,i3),1(1x,f10.5))')
     & '    Smax:i,j,k',umax,iumax,jumax,kumax,zumax
	IF(myid .EQ. itmax2)  write(*,'(a14,f10.4,3(1x,i3),1(1x,f10.5))')
     & '    Smin:i,j,k',umin,iumin,jumin,kumin,zumin
#else
	  write(*,'(a14,f10.4,3(1x,i3),1(1x,f10.5))')
     & '    Smax:i,j,k',umax,iumax,jumax,kumax,zumax
	  write(*,'(a14,f10.4,3(1x,i3),1(1x,f10.5))')
     & '    Smin:i,j,k',umin,iumin,jumin,kumin,zumin
#endif
	endif  !	if(mod(ip,10).eq.0)
#endif

c11199	  continue

c-----------------------------------------------------------------------
c      subroutine trmice:
c	  compute sea ice dynamics and temperature change of the model
c	  surface layer due to air sea fluxes
c-----------------------------------------------------------------------
	call trmice(iindex,izet,jjc,m,n,wspeed,bew,taup,
     *	  einstrice,einstr,taair,dt,dz,ilo,ltief,qw,khor,
     *	  qoi,qii)
c-----------------------------------------------------------------------
c     boundary conditions for T,S and biology
c-----------------------------------------------------------------------
c     set necessary equal boundary lines for T and S
      call tsrneu(iindex,izet,time,dz,delsal)  !set equal boundary lines for T,S
#ifdef MPI

      call INTI_COMM_VAR

      if((withice .and. (mod(ip+i_ice,i_ice).eq.0)) .or.
     &(mod(ip+ireadmet,ireadmet).eq.0)) then
c     &(mod(ip+ireadmet-1,ireadmet).eq.0)) then
      call com_vertical(ctemp,iel1,frice(1,1),iel,2,2,0,2,0,lazc)
      call com_vertical(ctemp,iel1,his(1,1),iel,2,2,0,2,0,lazc)
      call com_vertical(ctemp,iel1,hisr(1,1),iel,2,2,0,2,0,lazc)
      endif

      call com_vertical(tec,ndrei,dtemp,iel,2,2,1,3,0,lazc)
      call com_vertical(sac,ndrei,dtemp,iel,2,2,1,3,1,lazc)

      if((withice .and. (mod(ip+i_ice,i_ice).eq.0)) .or.
     &(mod(ip+ireadmet,ireadmet).eq. 0)) then
c     &(mod(ip+ireadmet-1,ireadmet).eq. 0)) then
      call com_horizontal(ctemp,iel1,frice(1,1),iel,2,2,0,2,0,
     *icdepb2,idepb2,icdept2
     *,idept2,ncsendb2,nsendb2,ncsendt2,nsendt2,nbhlayer2,nthlayer2)
      call com_horizontal(ctemp,iel1,his(1,1),iel,2,2,0,2,0,
     *icdepb2,idepb2,icdept2
     *,idept2,ncsendb2,nsendb2,ncsendt2,nsendt2,nbhlayer2,nthlayer2)
      call com_horizontal(ctemp,iel1,hisr(1,1),iel,2,2,0,2,0,
     *icdepb2,idepb2,icdept2
     *,idept2,ncsendb2,nsendb2,ncsendt2,nsendt2,nbhlayer2,nthlayer2)
      endif

      call com_horizontal(tec,ndrei,dtemp,iel,2,2,1,3,0,icdepb2,
     *idepb2,icdept2
     *,idept2,ncsendb2,nsendb2,ncsendt2,nsendt2,nbhlayer2,nthlayer2)
      call com_horizontal(sac,ndrei,dtemp,iel,2,2,1,3,1,icdepb2,
     *idepb2,icdept2
     *,idept2,ncsendb2,nsendb2,ncsendt2,nsendt2,nbhlayer2,nthlayer2)
#endif
#ifdef MPI
      call MPI_Waitall(nreq1, ir_r_b, status, ierr)
      INUM1=0
      INUM2=0
      if((withice .and. (mod(ip+i_ice,i_ice).eq.0)) .or.
     &(mod(ip+ireadmet,ireadmet).eq. 0))then
c     &(mod(ip+ireadmet-1,ireadmet).eq. 0))then
      call VER_UNPACK(ctemp,iel1,frice(1,1),iel,2,2,0,2,lazc)
      call VER_UNPACK(ctemp,iel1,his(1,1),iel,2,2,0,2,lazc)
      call VER_UNPACK(ctemp,iel1,hisr(1,1),iel,2,2,0,2,lazc)
      endif
      call VER_UNPACK(tec,ndrei,dtemp,iel1,2,2,1,3,lazc)
      call VER_UNPACK(sac,ndrei,dtemp,iel1,2,2,1,3,lazc)
      INUM3=0
      INUM4=0
      if((withice .and. (mod(ip+i_ice,i_ice).eq.0)) .or.
     &(mod(ip+ireadmet,ireadmet).eq. 0))then
c     &(mod(ip+ireadmet-1,ireadmet).eq. 0))then
      call HOR_UNPACK(ctemp,iel1,frice(1,1),iel,2,2,nbhlayer2,nthlayer2,
     &idepb2,idept2,icdepb2,icdept2,0,2,iwet)
      call HOR_UNPACK(ctemp,iel1,his(1,1),iel,2,2,nbhlayer2,nthlayer2,
     &idepb2,idept2,icdepb2,icdept2,0,2,iwet)
      call HOR_UNPACK(ctemp,iel1,hisr(1,1),iel,2,2,nbhlayer2,nthlayer2,
     &idepb2,idept2,icdepb2,icdept2,0,2,iwet)
      endif
      call HOR_UNPACK(tec,ndrei,dtemp,iel1,2,2,nbhlayer2,nthlayer2,
     &idepb2,idept2,icdepb2,icdept2,1,3,iwet)
      call HOR_UNPACK(sac,ndrei,dtemp,iel1,2,2,nbhlayer2,nthlayer2,
     &idepb2,idept2,icdepb2,icdept2,1,3,iwet)
#endif

      call estate
c-----------------------------------------------
      endif  !	if(progno) then
#ifdef MPI
      if(.not. progno) then
      call INTI_COMM_VAR
      if((withice .and. (mod(ip+i_ice,i_ice).eq.0)) .or.
     &(mod(ip+ireadmet,ireadmet).eq.0)) then
c     &(mod(ip+ireadmet-1,ireadmet).eq. 0))then
      call INTI_COMM_VAR
      call com_vertical(ctemp,iel1,frice(1,1),iel,2,2,0,2,0,lazc)
      call com_vertical(ctemp,iel1,his(1,1),iel,2,2,0,2,0,lazc)
      call com_vertical(ctemp,iel1,hisr(1,1),iel,2,2,0,2,1,lazc)
      call com_horizontal(ctemp,iel1,frice(1,1),iel,2,2,0,2,0,icdepb2,
     *idepb2,icdept2
     *,idept2,ncsendb2,nsendb2,ncsendt2,nsendt2,nbhlayer2,nthlayer2)
      call com_horizontal(ctemp,iel1,his(1,1),iel,2,2,0,2,0,icdepb2,
     *idepb2,icdept2
     *,idept2,ncsendb2,nsendb2,ncsendt2,nsendt2,nbhlayer2,nthlayer2)
      call com_horizontal(ctemp,iel1,hisr(1,1),iel,2,2,0,2,1,icdepb2,
     *idepb2,icdept2
     *,idept2,ncsendb2,nsendb2,ncsendt2,nsendt2,nbhlayer2,nthlayer2)
      endif
      endif
#endif

c-----------------------------------------------------------------------
c     prepare daily mean output: accumulate daily mean
c-----------------------------------------------------------------------
#ifdef MPI
      do I=1,khorl
      lw=lzet(I)
#else
      do lw=1,khor
#endif
	zmit(lw) = zmit(lw)+zac(lw)
      enddo
#ifdef MPI

      do J=J1,J2
      do i=I1,I2
#else
      do i=1,m
      do j=1,n
#endif
	frimit(i,j)  = frimit(i,j)+frice(i,j)
	hismit(i,j)  = hismit(i,j)+his(i,j)
	hisrmit(i,j) = hisrmit(i,j)+hisr(i,j)
	tismit(i,j)  = tismit(i,j)+tis(i,j)
	uimit(i,j)   = uimit(i,j)+ui(i,j)
	vimit(i,j)   = vimit(i,j)+vi(i,j)
      enddo
      enddo

#ifndef MPI
      do lw=1,ndrei
	umit(lw)  = umit(lw)+uc(lw)
	vmit(lw)  = vmit(lw)+vc(lw)
	scmit(lw) = scmit(lw)+sac(lw)
	wcmit(lw) = wcmit(lw)+wc(lw)
	szmit(lw) = szmit(lw)+szahl(lw)
	acmit(lw) = acmit(lw)+avc(lw)
	tcmit(lw) = tcmit(lw)+tec(lw)
      enddo
#else
      do J=J1,J2
      lwa=lb0(J)
      lwe=le0(J)
      if(lwe .ge. lwa) then
      nwet1=indwet(lwa)+1
      nwet2=indwet(lwe)+lazc(lwe)
      do k=nwet1,nwet2
	umit(k)	 = umit(k)+uc(k)
	vmit(k)	 = vmit(k)+vc(k)
	scmit(k) = scmit(k)+sac(k)
	wcmit(k) = wcmit(k)+wc(k)
	szmit(k) = szmit(k)+szahl(k)
	acmit(k) = acmit(k)+avc(k)
	tcmit(k) = tcmit(k)+tec(k)
      enddo
      endif
      enddo
#endif
c-----------------------------------------------------------------------
c    daily mean radiation
c-----------------------------------------------------------------------
!	     write(nout) fqgmit	 !net short wave radiation
!	     write(nout) fqrmit	 !net long wave radiation
!	     write(nout) fqsmit	 !turbulent heat flux
!	     write(nout) fqlmit	 !latent  neat flux
!	   write(nout) qois	 !ice-ocean heat flux
!	   write(nout) qiis	 !heat conduction through the ice
C	qwt(i)	= -evap(i)*xlat- hsens(i) SENSIBLE
C	fqs(i) = hsens(i)*(-1.0)
C	fql(i) = evap(i)*(-1.0)*xlat

#ifdef MPI
       do J=J1,J2
       do i=I1,I2
#else
	do i=1,m
	do j=1,n
#endif
	  fqgmit(i,j) = fqgmit(i,j)+fqg(i,j) !radiation
	  qois(i,j) = qois(i,j)+qoi(i,j)
	  qiis(i,j) = qiis(i,j)+qii(i,j)
	enddo
	enddo
#ifdef MPI
	do I=1,khorl
	lw=lzet(I)
#else
	do lw=1,khor
#endif
	  fqrmit(lw) = fqrmit(lw)+fqr(lw)
	  fqsmit(lw) = fqsmit(lw)+fqs(lw)
	  fqlmit(lw) = fqlmit(lw)+fql(lw)
	enddo

#ifdef MPI
      if(.not. progno) then
      if((withice .and. (mod(ip+i_ice,i_ice).eq.0)) .or.
     &(mod(ip+ireadmet,ireadmet).eq.0)) then
c     &(mod(ip+ireadmet-1,ireadmet).eq. 0))then

      call MPI_Waitall(nreq1, ir_r_b, status, ierr)
      call VER_UNPACK(ctemp,iel,frice(1,1),iel1,2,2,0,2,lazc)
      call VER_UNPACK(ctemp,iel,his(1,1),iel1,2,2,0,2,lazc)
      call VER_UNPACK(ctemp,iel,hisr(1,1),iel1,2,2,0,2,lazc)

      INUM3=0
      INUM4=0
      call HOR_UNPACK(ctemp,iel,frice(1,1),iel1,2,2,nbhlayer2,nthlayer2,
     &idepb2,idept2,icdepb2,icdept2,0,2,iwet)
      call HOR_UNPACK(ctemp,iel,his(1,1),iel1,2,2,nbhlayer2,nthlayer2,
     &idepb2,idept2,icdepb2,icdept2,0,2,iwet)
      call HOR_UNPACK(ctemp,iel,hisr(1,1),iel1,2,2,nbhlayer2,nthlayer2,
     &idepb2,idept2,icdepb2,icdept2,0,2,iwet)
      endif
      endif
#endif
       igstep = igstep+1
c-----------------------------------------------------------------------
c     check end of the day
c-----------------------------------------------------------------------

      if(igstep.eq.r1steps)then   !timstep for output
      igstep1=igstep1+1
      igstep = 0

c-----------------------------------------------------------------------
c	actual ice-mass esum
c-----------------------------------------------------------------------

c	esum=0.
c#ifdef MPI
c	esum1=0.
c	do j=J1,J2
c	do i=I1,I2
c#else
c	do i=1,m
c	do j=1,n
c#endif
c	  esum = esum+(his(i,j)+hisr(i,j))*frice(i,j)
c	enddo
c	enddo
#ifdef MPI
c       CALL MPI_ALLREDUCE(esum,esum1,1,MPI_REAL,MPI_SUM,
c     & MPI_COMM_WORLD,IERR)
c	esum=esum1
c-----------------------------------------------------------------------
c	compute daily mean
c-----------------------------------------------------------------------
	do I=1,khorl
	lw=lzet(I)
#else
	do lw=1,khor
#endif
	  zmit(lw)  = zmit(lw)*rsteps
	  zmalt(lw) = zmit(lw)
	  fqrmit(lw)  = fqrmit(lw)*rsteps
	  fqsmit(lw)  = fqsmit(lw)*rsteps
	  fqlmit(lw)  = fqlmit(lw)*rsteps
	enddo
#ifdef MPI
	do j=J1,J2
	do i=I1,I2
#else
	do i=1,m
	do j=1,n
#endif
	  fqgmit(i,j)	= fqgmit(i,j)*rsteps
	  qois(i,j)   = qois(i,j)*rsteps
	  qiis(i,j)   = qiis(i,j)*rsteps
	  hismit(i,j)  = hismit(i,j)*rsteps
	  frimit(i,j)  = frimit(i,j)*rsteps
	  hisrmit(i,j) = hisrmit(i,j)*rsteps
	  tismit(i,j)  = tismit(i,j)*rsteps
	  uimit(i,j)   = uimit(i,j)*rsteps
	  vimit(i,j)   = vimit(i,j)*rsteps
	enddo
	enddo
#ifndef MPI
	do lw=1,ndrei
	  umit(lw)  = umit(lw)*rsteps
	  vmit(lw)  = vmit(lw)*rsteps
	  wcmit(lw) = wcmit(lw)*rsteps*1000.0
	  szmit(lw) = szmit(lw)*rsteps
	  acmit(lw) = acmit(lw)*rsteps*10000.0
	  scmit(lw) = scmit(lw)*rsteps
	  tcmit(lw) = tcmit(lw)*rsteps
	enddo

#else
      do J=J1,J2
      lwa=lb0(J)
      lwe=le0(J)
      if(lwe .ge. lwa) then
      nwet1=indwet(lwa)+1
      nwet2=indwet(lwe)+lazc(lwe)
      do k=nwet1,nwet2
	  umit(k)  = umit(k)*rsteps
	  vmit(k)  = vmit(k)*rsteps
	  wcmit(k) = wcmit(k)*rsteps*1000.0
	  szmit(k) = szmit(k)*rsteps
	  acmit(k) = acmit(k)*rsteps*10000.0
	  scmit(k) = scmit(k)*rsteps
	  tcmit(k) = tcmit(k)*rsteps
	enddo
	endif
	enddo

#endif


	call ftrace_region_begin ('mitout')

c-------------------------------------------------------------------------------
c    subroutine IO_daily_physics: daily mean output
c    subroutine IO_start_physics: creates a backup file every month
c
c-------------------------------------------------------------------------------
	IF(mitout) THEN
#ifndef MPI
	  call IO_daily_physics(nout,mjar,lmon,nday,ivier,imal,1) ! nnn write=1 or read=2
c	if(nday.eq.itagend)	 ! monthly backup
	if(igstep1.eq.itagend*toutsteps)	 ! monthly backup
     &	  call IO_start_physics(nout,mjar,lmon,nday,ivier,imal,1) ! nnn write=1 or read=2
#else
       call IO_daily_physics_p(nout1,mjar,lmon,nday,ivier,ip,imal)
      IF(myid .EQ. 0)THEN
      print*,'txc,tyc,pac writestart ',minval(txc),maxval(txc),
     &maxval(tyc),minval(pac),maxval(pac)
      endif
       
       if(igstep1.eq.itagend*toutsteps) ! monthly backup
     & call IO_start_physics_p(nout1,mjar,lmon,nday,ivier,ip,imal)
#endif
#ifdef _FABM_
       call ftrace_region_begin ('biological_outpu')
#ifndef MPI
       write (nout_b) cc1d
#else
       DO ib=1,5  ! npel = 5 here hardcoded needs to be changed 
       DO J=J1,J2
       lwa=lb0(J)
       lwe=le0(J)
       if(lwe .ge. lwa) then
       nwet1=indwet(lwa)+1
       nwet2=indwet(lwe)+lazc(lwe)
       write (nout_b1) (cc1d(k,ib),k=nwet1,nwet2)
       endif
       enddo
       enddo
#endif
#endif
      ENDIF!	    if(miaout) then
	call ftrace_region_end ('mitout')

 2006	format(8(1x,e18.12))

	if(myid .eq. 0) write(7,*)'OUTPUT-FILE fuer den',mjar,monat,nday

c-----------------------------------------------------------------------
c	set daily mean to zero
c-----------------------------------------------------------------------

	do lw = 1,khor
	  zmit(lw)   = 0.0
	  fqrmit(lw) = 0.0
	  fqsmit(lw) = 0.0
	  fqlmit(lw) = 0.0
	enddo

	do lw = 1,ndrei
	  umit(lw)  = 0.0
	  vmit(lw)  = 0.0
	  tcmit(lw) = 0.0
	  szmit(lw) = 0.0
	  scmit(lw) = 0.0
	  wcmit(lw) = 0.0
	  acmit(lw) = 0.0
	enddo

      do i=1,m
      do j=1,n
       fqgmit(i,j)  = 0.0
       qois(i,j)    = 0.0
       qiis(i,j)    = 0.0
       hismit(i,j)  = 0.0
       tismit(i,j)  = 0.0
       hisrmit(i,j) = 0.0
       uimit(i,j)   = 0.0
       vimit(i,j)   = 0.0
       frimit(i,j)  = 0.0
      enddo
      enddo
      endif !(end output igstep .eq. r1step)
      time = time+dt
 1000 continue
c-----------------------------------------------------------------------
c     end of dt loop
c-----------------------------------------------------------------------

      if(myid .eq. 0)write(7,606) ivier,mjar,lmon,nday,(itel(i),i=1,12)
  606 format (2x,'termin',i2,3x,'j,m,d, /',3i4,3x,
     * 'steps / iterationen = ',12i5)

 1010 continue
c----------------------------------------------------------------------
c     end of date loop
c-----------------------------------------------------------------------
      call ftrace_region_end ('time_loop')

c-------------------------------------------------------------------------
c    create a start-up file
c------------------------------------------------------------------------
cftrace	     call ftrace_region_begin ('daily_backup')
      if(backup) then ! _start file

	  startfile(1:3)=ppp
	  write(startfile(4:5),'(i2.2)') mjar
	  write(startfile(6:7),'(i2.2)') monat
	  write(startfile(8:9),'(i2.2)') nday
       open(nout_start,
     & file=pathpc2//pathout//startfile(1:9)//'_start'
     &,form='unformatted')
       call IO_start_physics(nout_start,mjar,lmon,nday,ivier,ip,1) ! nnn write=1 or read=2
       close(nout_start)

      endif !  if(backup) then
cftrace	     call ftrace_region_end ('daily_backup')

 1111 continue
c-----------------------------------------------------------------------
c     end of day loop
c-----------------------------------------------------------------------

c	metinput
	close(unit=67)
	close(unit=68)
	close(unit=69)
	close(unit=70)
	close(unit=71)
	close(unit=72)
	close(unit=73)
	close(unit=74)
	close(unit=75)

**z_v  file monthly to close
c	 close(unit=99)
c	openboundary files :monthly
	close(unit=18)
	close(unit=19)
	close(unit=65)
	close(unit=15)
c-----------------------------------------------------------------------
c     end of month loop
c     close Output daily means in monthly files, unformatted BINARY
c-----------------------------------------------------------------------

#ifdef MPI
      close(122)
      close(nout1)
      close(nout_b1)
      close(nout_f1)

      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
      nmetpermon=nreadmet*itagend
       if(myid .eq. nprocs-1)
     & call combine_meteo_data(nanz,outfile1,pathout,nmetpermon)

       if(myid .eq. 0) then
       if(lmon .eq. 12) then
       startfile(1:3)=ppp
       write(startfile(4:5),'(i2.2)') mjar
       write(startfile(6:7),'(i2.2)') monat
       write(startfile(8:9),'(i2.2)') nday-1
       open(nout_start,
     & file=pathpc2//pathout//startfile(1:9)//'_start'
     &,form='unformatted',status='unknown',access='stream',
     &action='write')
       endif
       call combine_output_daily(nout,pathout,monat,ppp,
     & mjar,icord,itagend,nout_start,lmon,toutsteps)
       if(lmon .eq. 12) close(nout_start)
       endif

       if(myid .eq. 1) then
       if(lmon .eq. 12) then
       startfile(1:3)=ppp
       write(startfile(4:5),'(i2.2)') mjar
       write(startfile(6:7),'(i2.2)') monat
       write(startfile(8:9),'(i2.2)') nday-1
       startfile(10:10)='b'
       open(nout_start,
     & file=pathpc2//pathout//startfile//'_start'
     &,form='unformatted')
       endif

       if(lmon .eq. 12) close(nout_start)

       call combine_output_bio(nout_b,pathout,monat,ppp,
     & mjar,icord,itagend,nout_start,lmon,toutsteps,5)  !npel=5 hardcoded
       endif
       CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
#endif
      close(nout  ) !	= 85		   ! physical	  variables Monthly output
      close(nout_b) !	= 86		   ! production	  variables Monthly output
      close(nout_f) !	= 87		   ! production	  variables Monthly output
      close(nout_s) !	= 88		   ! space fluxes variables Monthly output


      write(7,6999) mjar,monat,itaga,ltag
 6999 format (/3x,'jahr, monat /',2i5,'	 von tag',i4,
     * '  bis tag',i4,'	 simuliert')

 2111 continue
c-----------------------------------------------------------------------
c     end of month loop
c     close Output daily means in monthly files, unformatted BINARY
c-----------------------------------------------------------------------

c	radi, evap

c	 close(unit=77)
c	 close(unit=78)
        close(nfurin)
        close(20)
        close(21)

c----------------------------------------------------------------------------
c     write start file for 31.12.yyyy
c----------------------------------------------------------------------------
#ifndef MPI
	  startfile(1:3)=ppp
	  write(startfile(4:5),'(i2.2)') mjar
	  write(startfile(6:7),'(i2.2)') monat
	  write(startfile(8:9),'(i2.2)') nday-1
       open(nout_start,
     & file=pathpc2//pathout//startfile(1:9)//'_start'
     &,form='unformatted')

       call IO_start_physics(nout_start,mjar,lmon,nday,ivier,ip,1) ! nnn write=1 or read=2
       close(nout_start)
#endif
       close(121)

 5000	continue ! DO years
c	 close(unit=19) nzac

	close(7)  !control file
       close(98)
      call system_clock(count=iclock1)
       t = (iclock1 - iclock0)
       time = t/real(ihz)

       open(111,file='INTEL_cputime')
c	open(111,file='PGI_cputime')
       do iii=1,999
       read(111,*,end=129)readt
       enddo

 129   continue
c	write(111,*)'YEAR: FFLAGS = -fpp -fast:', time, time/60,'[min]'	     !1
c	write(111,*)'FFLAGS = -fpp -fast -mtune=pentium4:', time, time/60,'[min]'    !2
c	write(111,*)'FFLAGS = -fpp -fast -mtune=pentium4 -parallel:', time, time/60,'[min]'    !3

	write(111,*) ppp,ijahr,' file 98'
c	 write(111,*)'-O2 -ipo -no-prec-div -static -xW -mtune=pentium4'


c	 write(111,*)'YEAR:p79,498:'
	  write(111,*)'-fpp -O3'


c	write(111,*)'YEAR: f74 FFLAGS = -fast:',time,time/60,'[min]' !1
c	write(111,*)'YEAR e75: FFLAGS = -g:',time,time/60,'[min]' !1
c	write(111,*)'FFLAGS = -fastsse:', time	!2
c	write(111,*)'FFLAGS = -fastsse –Mipa=fast:', time  !3
c	write(111,*)'FFLAGS = -Mipa=fast:', time   !4
c	write(111,*)'FFLAGS = –Mvect=sse:', time   !5
c	write(111,*)'FFLAGS = –Mconcur=dist:block:', time   !6
c	write(111,*)'FFLAGS =  -fast –Mipa=inline:', time   !8

	write(111,*)time/60.,'[min]',time/60./60.,'[h]'	   !

       close(111)
c-----------------------------------------------------------------------
      stop'normaler programm **** s t o p ****'
c-----------------------------------------------------------------------
	call MPI_FINALIZE(ierr)
	end

#ifndef FTRACE
      subroutine ftrace_region_begin(str)
      character*(*) str
      return
      end
      subroutine ftrace_region_end(str)
      character*(*) str
      return
      end
#endif
