module hamsom_fabm

#ifdef _FABM_

use fabm
use fabm_types
#define DEBUG
#ifdef DEBUG
use fabm_debug
#endif
implicit NONE

private

   EXTERNAL comp2d1d, comp3d1d, deco1d2d, deco1d3d, deco1d2di

   class (type_fabm_model), pointer :: model
   integer :: ny,nx,nz ! size of the FABM 3D model
   integer, dimension(:), allocatable :: variable_order

   integer :: fabmunit=120
   integer :: stat
   integer :: m,n,ilo
   integer :: nbio,nsed,khor,ndrei

   real, allocatable, dimension(:,:), target :: bottom_stress 
   real, allocatable, dimension(:,:,:,:), target :: pelagic 
   real, allocatable, dimension(:,:,:), target :: h 
   real, allocatable, dimension(:,:,:), target :: sediment 
   logical, allocatable, dimension(:,:,:) :: mask 
   integer, allocatable, dimension(:,:), target :: bindx 
   real, allocatable, dimension(:,:,:), target :: surf_flux 
   real, allocatable, dimension(:,:,:), target :: surf_sms 
   real, allocatable, dimension(:,:,:,:), target :: interior_sms
   real, allocatable, dimension(:,:,:), target :: bott_flux 
   real, allocatable, dimension(:,:,:), target :: bott_sms

#ifndef MPI
   integer :: myid=0
#endif

   integer :: il,ih,jl,jh

   public configure_fabm, initialize_fabm, update_fabm, clean_fabm
   public fabm_var_index

contains

#ifndef MPI
subroutine configure_fabm(m_,n_,ilo_,khor_,ndrei_,nbio_,nsed_)
#else
subroutine configure_fabm(myid,m_,n_,ilo_,khor_,ndrei_,nbio_,nsed_)
   integer, intent(in) :: myid
#endif
   integer, intent(in) :: m_,n_,ilo_ ! y (m -> from north), x (n -> from west), z (ilo -> from top)
   integer, intent(in) :: khor_
   integer, intent(in) :: ndrei_
   integer, intent(in) :: nbio_
   integer, intent(in) :: nsed_

   integer :: ivar,nvar

   m=m_
   n=n_
   ilo=ilo_
   ndrei=ndrei_
   khor=khor_
   nbio=nbio_
   nsed=nsed_

   if (myid .eq. 0) then
      write(fabmunit,*) 'configure_fabm()'
      write(fabmunit,*) '  m,n,ilo= ',m,n,ilo
      write(fabmunit,*) '  ndrei=   ',ndrei
      write(fabmunit,*) '  khor=    ',khor
      write(fabmunit,*) '  nbio=    ',nbio
      write(fabmunit,*) '  nsed=    ',nsed
   end if

   model => fabm_create_model('fabm.yaml')

   call allocate_fabm(myid,m,n,ilo,nbio,nsed)

   if (myid .eq. 0) then
      do ivar=1,size(model%interior_state_variables)
         write(fabmunit,*) ivar,trim(model%interior_state_variables(ivar)%name)
      end do
   end if
   if (myid .eq. 0) then
      do ivar=1,size(model%bottom_state_variables)
         write(fabmunit,*) ivar,trim(model%bottom_state_variables(ivar)%name)
      end do
   end if
   return
end subroutine configure_fabm

!-----------------------------------------------------------------------

#ifndef MPI
subroutine initialize_fabm(lazc,dt,iwet,indend,icord,pd2,Tc,Tsed,einstr,taubot)
#else
subroutine initialize_fabm(myid,dt,lazc,iwet,indend,icord,pd2,Tc,Tsed,einstr,taubot)
   integer, intent(in) :: myid
#endif
   real, intent(in) :: dt
   integer, intent(in) :: lazc(:),iwet(:),indend(:),icord(:,:)
   real, intent(in) :: pd2(:,:,:)
   real, intent(in) :: Tc(:,:)
   real, intent(in) :: Tsed(:,:)
   real, intent(in) :: einstr(:,:)
   real, intent(in) :: taubot(:,:)

   integer :: i,j,k,l
   integer :: ivar

   write(fabmunit,*) 'initialize_fabm()',myid

   if(myid .eq. 0) then
       write(fabmunit,*) size(lazc),shape(lazc),rank(lazc)
       write(fabmunit,*) lbound(lazc),ubound(lazc)
       write(fabmunit,*) size(iwet),shape(iwet),rank(iwet)
       write(fabmunit,*) lbound(iwet),ubound(iwet)
       write(fabmunit,*) size(indend),shape(indend),rank(indend)
       write(fabmunit,*) lbound(indend),ubound(indend)
       write(fabmunit,*) size(Tc),shape(Tc),rank(Tc)
       write(fabmunit,*) lbound(Tc),ubound(Tc)
       write(fabmunit,*) size(Tsed),shape(Tsed),rank(Tsed)
       write(fabmunit,*) lbound(Tsed),ubound(Tsed)
   end if

   ! FABM pelagics is being 'initialized'
   pelagic = -10.
   do l=1,nbio
#ifdef MPI
      call deco1d3d_s(pelagic(:,:,:,l),Tc(:,l),ndrei)
#else
      call deco1d3d(pelagic(:,:,:,l),Tc(:,l),ndrei)
#endif
   end do

   do l=1,nsed
#ifdef MPI
      call deco1d2d_s(sediment(:,:,l),Tsed(:,l),iwet,indend,1)
#else
      call deco1d2d(sediment(:,:,l),Tsed(:,l),iwet,indend,1)
#endif
   end do

   ! here we initialize the 3D calculation mask
   bindx = 0
   ! is serial and parallel the same?
   call deco1d2di(bindx,lazc,iwet,indend,1)
   mask = .false.
   do i=1,n
      do j=1,m
         if(bindx(j,i) .gt. 0) mask(j,i,1:bindx(j,i)) = .true.
      end do
   end do
   call print_mask(myid,m,n,icord)

   ! to save some typing
   jl = icord(myid+1,1)
   jh = icord(myid+1,2)
   il = icord(myid+1,3)
   ih = icord(myid+1,4)

   call model%set_domain(m,n,ilo,dt)
   call model%set_domain_start(jl,il,1)
   call model%set_domain_stop(jh,ih,ilo) 
   call model%set_bottom_index(bindx(:,:))
   call model%set_mask(mask,mask(:,:,1))

#if 0
   write(*,*) jl,jh
   write(*,*) il,ih
   write(*,*) 1,ilo
   write(*,*) size(pd2(1:m,1:n,1:ilo))
#endif
   call model%link_interior_data(fabm_standard_variables%cell_thickness, pd2(1:m,1:n,1:ilo))

   ! set pointers to environmental forcing
   call model%link_horizontal_data(fabm_standard_variables%surface_downwelling_shortwave_flux, einstr(1:m,1:n))
!   call model%link_horizontal_data(standard_variables%wind_speed,wspd_fabm(1:ii, 1:jj))
!   call model%link_horizontal_data(standard_variables%mole_fraction_of_carbon_dioxide_in_air,atmco2_fabm(1:ii,1:jj))
   call model%link_horizontal_data(fabm_standard_variables%bottom_stress, taubot(1:m,1:n))

   call model%link_interior_data(fabm_standard_variables%temperature, pelagic(:,:,:,1))
   call model%link_interior_data(fabm_standard_variables%practical_salinity, pelagic(:,:,:,2))

   ! link to pelagic state variables
   do ivar = 1,size(model%interior_state_variables)
      call model%link_interior_state_data(ivar,pelagic(:,:,:,ivar+2))
   end do
   ! link to benthic state variables
   do ivar = 1,size(model%bottom_state_variables)
      call model%link_bottom_state_data(ivar,sediment(:,:,ivar))
   end do

   call model%start()
!KB - check this - only if not restart
#if 1
   do k=1,ilo
      do i=il,ih
         call model%initialize_interior_state(jl,jh,i,k)
      end do      
   end do      
!   call model%initialize_bottom_state()
#endif   

   return
end subroutine initialize_fabm

!-----------------------------------------------------------------------

subroutine update_fabm(myid,imal,dt,iwet,indend,Tc,Tsed)
   integer, intent(in) :: myid
   integer, intent(in) :: imal
   real, intent(in) :: dt
   integer, intent(in) :: iwet(:),indend(:)
   real, dimension(:,:), intent(inout) :: Tc
   real, dimension(:,:), intent(inout) :: Tsed

   integer :: i,j,k,l

   if (myid .eq. 0) then
      write(*,*) 'update_fabm()'
   end if

   do l=1,nbio
#ifdef MPI
      call deco1d3d_s(pelagic(:,:,:,l),Tc(:,l),ndrei)
#else
      call deco1d3d(pelagic(:,:,:,l),Tc(:,l),ndrei)
#endif
   end do
   do l=1,nsed
#ifdef MPI
      call deco1d2d_s(sediment(:,:,l),Tsed(:,l),iwet,indend,1)
#else
      call deco1d2d(sediment(:,:,l),Tsed(:,l),iwet,indend,1)
#endif
   end do

   call model%prepare_inputs(t=imal*dt)

! get_vertical_move - maybe send individually

   ! here the surface is updated
   surf_flux = 0.
   do i=il,ih
      call model%get_surface_sources(jl,jh,i,surf_flux(jl:jh,i,:))
   end do

   ! here the pelagic is updated
   interior_sms = 0.
   do k=1,ilo
      do i=il,ih
         call model%get_interior_sources(jl,jh,i,k,interior_sms(jl:jh,i,k,:))
      end do
   end do

   ! here the bottom is updated
   bott_flux = 0.
   bott_sms = 0.
   do i=il,ih
      call model%get_bottom_sources(jl,jh,i,bott_flux(jl:jh,i,:),bott_sms(jl:jh,i,:))
   end do

   call model%finalize_outputs()

   do l=3,nbio
#ifdef MPI
      call comp3d1d_s(pelagic(:,:,:,l),Tc(:,l))
#else
      call comp3d1d(pelagic(:,:,:,l),Tc(:,l))
#endif
   end do
   do l=1,nsed
#ifdef MPI
      call comp2d1d(sediment(:,:,l),Tsed(:,l))
!KB      call comp2d1d_s(sediment(:,:,l),Tsed(:,l))
#else
      call comp2d1d(sediment(:,:,l),Tsed(:,l))
#endif
   end do

   return
end subroutine update_fabm

!-----------------------------------------------------------------------

subroutine clean_fabm(myid)
   integer, intent(in) :: myid
   if (myid .eq. 0) then
      write(fabmunit,*) 'clean_fabm()'
   end if
end subroutine clean_fabm

!-----------------------------------------------------------------------

subroutine allocate_fabm(myid,m,n,ilo,nbio,nsed)
   integer, intent(in) :: myid
   integer, intent(in) :: m,n,ilo ! y (m -> from north), x (n -> from west), z (ilo -> from top)
   integer, intent(in) :: nbio
   integer, intent(in) :: nsed

   if (myid .eq. 0) then
      write(fabmunit,*) 'allocate_fabm()'
   end if

   allocate(bottom_stress(m,n),stat=stat)
   if (stat /= 0) stop 'allocate_fabm(): Error allocating memory (bottom_stress)'

   allocate(h(m,n,ilo),stat=stat)
   if (stat /= 0) stop 'allocate_fabm(): Error allocating memory (h)'

!   allocate(variable_order(3:nbio),stat=stat)
!   if (stat /= 0) stop 'allocate_fabm(): Error allocating memory (variable_order)'

   allocate(pelagic(m,n,ilo,nbio),stat=stat)
   if (stat /= 0) stop 'allocate_fabm(): Error allocating memory (pelagic)'

   allocate(sediment(m,n,nsed),stat=stat)
   if (stat /= 0) stop 'allocate_fabm(): Error allocating memory (sediment)'

   allocate(mask(m,n,ilo),stat=stat)
   if (stat /= 0) stop 'allocate_fabm(): Error allocating memory (mask)'

   allocate(bindx(m,n),stat=stat)
   if (stat /= 0) stop 'allocate_fabm(): Error allocating memory (bindx)'

   allocate(surf_flux(m,n,nbio-2),stat=stat)
   if (stat /= 0) stop 'allocate_fabm(): Error allocating memory (surf_flux)'

!KB   allocate(interior_sms(m,n,ilo,3:nbio),stat=stat)
   allocate(interior_sms(m,n,ilo,nbio-2),stat=stat)
   if (stat /= 0) stop 'allocate_fabm(): Error allocating memory (interior_sms)'

   allocate(bott_flux(m,n,nbio-2),stat=stat)
   if (stat /= 0) stop 'allocate_fabm(): Error allocating memory (bott_flux)'

   allocate(bott_sms(m,n,nsed),stat=stat)
   if (stat /= 0) stop 'allocate_fabm(): Error allocating memory (bott_sms)'

   return
end subroutine allocate_fabm

function fabm_var_index(varname) result(nvar)
   character(len=*), intent(in) :: varname
   integer :: nvar

   integer :: idx,n

   do n=1,size(model%interior_state_variables)
      idx = index(trim(model%interior_state_variables(n)%name),trim(varname),back=.true.)
      ! a pelagic state variable has been found
      if (idx > 0) then
         nvar = n
         return
      end if
   end do

   do n=1,size(model%bottom_state_variables)
      idx = index(trim(model%bottom_state_variables(n)%name),trim(varname),back=.true.)
      ! a bottom state variable has been found
      if (idx > 0) then
         nvar = n
         return
      end if
   end do
   nvar = 0
   return
end function fabm_var_index

!-----------------------------------------------------------------------

subroutine print_mask(myid,m,n,icord)
   integer, intent(in) :: myid
   integer, intent(in) :: m,n
   integer, intent(in) :: icord(:,:)
   integer :: i,j,l

   if (myid .eq. 0) then
      write(fabmunit,*) 'global mask'
      do j=1,m
         write(fabmunit,'(5000(L1))') (mask(j,i,1), i=1,n)
      end do
      do l=1,18
         write(fabmunit,*) 'mask for domain# ',l
         write(fabmunit,*) icord(l,1),icord(l,2),icord(l,3),icord(l,4)
         do j=icord(l,1),icord(l,2)
            write(fabmunit,'(5000(L1))') (mask(j,i,1), i=icord(l,3),icord(l,4))
         end do
      end do
   end if
end subroutine print_mask
#endif

!-----------------------------------------------------------------------

end module hamsom_fabm
