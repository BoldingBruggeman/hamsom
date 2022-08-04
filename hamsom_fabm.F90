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

   EXTERNAL deco1d2di

   class (type_fabm_model), pointer :: model
   integer :: ny,nx,nz ! size of the FABM 3D model 
                       ! y (ny -> from north), x (nx -> from west), z (nz -> from top)

   integer :: fabmunit=120
   integer :: stat
   integer :: npel,nsed

   real, public, allocatable, dimension(:,:,:,:), target :: pelagic 
   real, public, allocatable, dimension(:,:,:), target :: sediment 
   real, public, allocatable, dimension(:,:,:,:), target :: interior_sms
   real, public, allocatable, dimension(:,:,:), target :: bottom_sms

   real, allocatable, dimension(:,:,:,:), target :: vertical_movement

   real, allocatable, dimension(:,:,:), target :: surf_flux 
   real, allocatable, dimension(:,:,:), target :: surf_sms 
   real, allocatable, dimension(:,:,:), target :: bott_flux 
   real, allocatable, dimension(:,:), target :: bottom_stress 
   real, allocatable, dimension(:,:,:), target :: h 
   logical, allocatable, dimension(:,:,:) :: fabm_mask 
   integer, allocatable, dimension(:,:), target :: bindx 

#ifndef MPI
   integer :: myid=0
#endif

   integer :: il,ih,jl,jh

   public configure_fabm, initialize_fabm, update_fabm, clean_fabm
   public fabm_var_index

contains

#ifndef MPI
subroutine configure_fabm(m,n,ilo,nbio,nsed_)
#else
subroutine configure_fabm(myid,m,n,ilo,nbio,nsed_)
   integer, intent(in) :: myid
#endif
   integer, intent(in) :: m,n,ilo ! y (m -> from north), x (n -> from west), z (ilo -> from top)
   integer, intent(in) :: nbio
   integer, intent(in) :: nsed_

   integer :: ivar,nvar

   ny=m
   nx=n
   nz=ilo
   npel=nbio
   nsed=nsed_

   if (myid .eq. 0) then
      write(fabmunit,*) 'configure_fabm()'
      write(fabmunit,*) '  ny,nx,nz= ',ny,nx,nz
      write(fabmunit,*) '  npel=     ',npel
      write(fabmunit,*) '  nsed=     ',nsed
   end if

   model => fabm_create_model('fabm.yaml')

   call allocate_fabm(myid,npel,nsed)

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
end subroutine configure_fabm

!-----------------------------------------------------------------------

#ifndef MPI
subroutine initialize_fabm(lazc,dt,iwet,indend,icord,pd2,einstr,taubot)
#else
subroutine initialize_fabm(myid,dt,lazc,iwet,indend,icord,pd2,einstr,taubot)
   integer, intent(in) :: myid
#endif
   real, intent(in) :: dt
   integer, intent(in) :: lazc(:),iwet(:),indend(:),icord(:,:)
   real, intent(in) :: pd2(:,:,:)
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
   end if

   ! FABM pelagic and sediment is being 'initialized' in the calling program

   ! here we initialize the 3D calculation mask
   bindx = 0
   ! is serial and parallel the same?
   call deco1d2di(bindx,lazc,iwet,indend,1)
   fabm_mask = .false.
   do i=1,nx
      do j=1,ny
         if(bindx(j,i) .gt. 0) fabm_mask(j,i,1:bindx(j,i)) = .true.
      end do
   end do
   call print_mask(myid,ny,nx,icord)

   ! to save some typing
   jl = icord(myid+1,1)
   jh = icord(myid+1,2)
   il = icord(myid+1,3)
   ih = icord(myid+1,4)

   call model%set_domain(ny,nx,nz,dt)
   call model%set_domain_start(jl,il,1)
   call model%set_domain_stop(jh,ih,nz) 
   call model%set_bottom_index(bindx(:,:))
   call model%set_mask(fabm_mask,fabm_mask(:,:,1))

   call model%link_interior_data(fabm_standard_variables%cell_thickness, pd2(1:ny,1:nx,1:nz))

   ! set pointers to environmental forcing
   call model%link_horizontal_data(fabm_standard_variables%surface_downwelling_shortwave_flux, einstr(1:ny,1:nx))
!   call model%link_horizontal_data(standard_variables%wind_speed,wspd_fabm(1:ii, 1:jj))
!   call model%link_horizontal_data(standard_variables%mole_fraction_of_carbon_dioxide_in_air,atmco2_fabm(1:ii,1:jj))
   call model%link_horizontal_data(fabm_standard_variables%bottom_stress, taubot(1:ny,1:nx))

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
   do k=1,nz
      do i=il,ih
         call model%initialize_interior_state(jl,jh,i,k)
      end do      
   end do      
   do i=il,ih
      call model%initialize_bottom_state(jl,jh,i)
   end do      
#endif   
end subroutine initialize_fabm

!-----------------------------------------------------------------------

subroutine update_fabm(myid,imal,ltief,pd2)
   integer, intent(in) :: myid
   integer, intent(in) :: imal
   integer, intent(in) :: ltief(:,:)
   real, dimension(:,:,:), intent(inout) :: pd2

   integer :: i,j,k,l

   if (myid .eq. 0) then
      write(*,*) 'update_fabm()'
   end if

   call model%prepare_inputs(t=real(imal,rk))

   ! here the surface is updated
   surf_flux = 0.
   do i=il,ih
      call model%get_surface_sources(jl,jh,i,surf_flux(jl:jh,i,:))
   end do

   ! here the pelagic is updated
   interior_sms = 0.
   do k=1,nz
      do i=il,ih
         call model%get_interior_sources(jl,jh,i,k,interior_sms(jl:jh,i,k,:))
      end do
   end do

   ! here the bottom is updated
   bott_flux = 0.
   bottom_sms = 0.
   do i=il,ih
      call model%get_bottom_sources(jl,jh,i,bott_flux(jl:jh,i,:),bottom_sms(jl:jh,i,:))
   end do

   ! fold the surface and bottom flux terms
   do i=il,ih
      do j=jl,jh
         if (ltief(j,i) /= 0) then
            ! surface
            k=1
            interior_sms(j,i,k,:)=interior_sms(j,i,k,:)+surf_flux(j,i,:)/pd2(j,i,k)
            ! bottom
            k=ltief(j,i)
            interior_sms(j,i,k,:)=interior_sms(j,i,k,:)+bott_flux(j,i,:)/pd2(j,i,k)
         end if
      end do
   end do

   ! vertical velocities
   do k=1,nz
      do i=il,ih
         call model%get_vertical_movement(jl,jh,i,k,vertical_movement(jl:jh,i,k,:))
      end do
   end do

   call model%finalize_outputs()
end subroutine update_fabm

!-----------------------------------------------------------------------

subroutine clean_fabm(myid)
   integer, intent(in) :: myid
   if (myid .eq. 0) then
      write(fabmunit,*) 'clean_fabm()'
   end if
end subroutine clean_fabm

!-----------------------------------------------------------------------

subroutine allocate_fabm(myid,npel,nsed)
   integer, intent(in) :: myid
   integer, intent(in) :: npel
   integer, intent(in) :: nsed

   if (myid .eq. 0) then
      write(fabmunit,*) 'allocate_fabm()'
   end if

   allocate(bottom_stress(ny,nx),stat=stat)
   if (stat /= 0) stop 'allocate_fabm(): Error allocating memory (bottom_stress)'

!KB   allocate(h(ny,nx,nz),stat=stat)
!KB   if (stat /= 0) stop 'allocate_fabm(): Error allocating memory (h)'

!KB   allocate(variable_order(3:npel),stat=stat)
!KB   if (stat /= 0) stop 'allocate_fabm(): Error allocating memory (variable_order)'

   allocate(pelagic(ny,nx,nz,npel),stat=stat)
   if (stat /= 0) stop 'allocate_fabm(): Error allocating memory (pelagic)'

   allocate(sediment(ny,nx,nsed),stat=stat)
   if (stat /= 0) stop 'allocate_fabm(): Error allocating memory (sediment)'

   allocate(fabm_mask(ny,nx,nz),stat=stat)
   if (stat /= 0) stop 'allocate_fabm(): Error allocating memory (fabm_mask)'

   allocate(bindx(ny,nx),stat=stat)
   if (stat /= 0) stop 'allocate_fabm(): Error allocating memory (bindx)'

   allocate(surf_flux(ny,nx,npel-2),stat=stat)
   if (stat /= 0) stop 'allocate_fabm(): Error allocating memory (surf_flux)'

   allocate(interior_sms(ny,nx,nz,npel-2),stat=stat)
   if (stat /= 0) stop 'allocate_fabm(): Error allocating memory (interior_sms)'

   allocate(bott_flux(ny,nx,npel-2),stat=stat)
   if (stat /= 0) stop 'allocate_fabm(): Error allocating memory (bott_flux)'

   allocate(bottom_sms(ny,nx,nsed),stat=stat)
   if (stat /= 0) stop 'allocate_fabm(): Error allocating memory (bott_sms)'

   allocate(vertical_movement(ny,nx,nz,npel-2),stat=stat)
   if (stat /= 0) stop 'allocate_fabm(): Error allocating memory (vertical_movement)'
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
end function fabm_var_index

subroutine print_mask(myid,ny,nx,icord)
   integer, intent(in) :: myid
   integer, intent(in) :: ny,nx
   integer, intent(in) :: icord(:,:)
   integer :: i,j,l,ic

   ic=size(icord,1)
   if (myid .eq. 0) then
      write(fabmunit,*) 'global mask'
      do j=1,ny
         write(fabmunit,'(5000(L1))') (fabm_mask(j,i,1), i=1,nx)
      end do
      do l=1,ic
         write(fabmunit,*) 'fabm_mask for domain# ',l
         write(fabmunit,*) icord(l,1),icord(l,2),icord(l,3),icord(l,4)
         do j=icord(l,1),icord(l,2)
            write(fabmunit,'(5000(L1))') (fabm_mask(j,i,1), i=icord(l,3),icord(l,4))
         end do
      end do
   end if
end subroutine print_mask
#endif

!-----------------------------------------------------------------------

end module hamsom_fabm
