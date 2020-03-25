module hamsom_fabm

#ifdef _FABM_

use fabm
use fabm_config
use fabm_types, only: attribute_length, output_none
use fabm_standard_variables, only: type_global_standard_variable
implicit NONE

private

   EXTERNAL comp3d1d, deco1d3d

   integer :: stat

   real, allocatable, dimension(:,:), public :: cc1d
   real, allocatable, dimension(:,:,:,:) :: cc 
   logical, allocatable, dimension(:,:,:), target :: mask 
   real, allocatable, dimension(:,:,:), target :: S, T
   ! use pointers to HAMSOM state variables

   public init_hamsom_fabm, do_hamsom_fabm, clean_hamsom_fabm

contains

subroutine init_hamsom_fabm(m,n,ilo,ndrei,npel,sac,tec)
   integer, intent(in) :: m,n,ilo
   integer, intent(in) :: ndrei
   integer, intent(in) :: npel
   real, intent(in) :: sac(:),tec(:)

   integer :: k,l

   write(*,*) 'Inside init_hamsom_fabm()'
   write(*,*) '  m,n,ilo= ',m,n,ilo
   write(*,*) '  ndrei=   ',ndrei
   write(*,*) '  npel=    ',npel

   call allocate_hamsom_fabm(m,n,ilo,ndrei,npel)

   do l=1,npel
      do k=1,ndrei
         cc1d(k,l) = k+l/10.
      end do
   end do

   cc = -10.
   mask = .false.

   ! FABM pelagics is being 'initialized'
   do l=1,npel
      call deco1d3d(cc(:,:,:,l),cc1d(:,l),ndrei)
   end do

!   call deco1d3d(S,sac,ndrei)
!   call deco1d3d(T,tec,ndrei)

   ! here we initialize the 3D mask - Ute - check
!   where (S(:,:,:) .gt. 0.) mask = 1
   where (cc(:,:,:,1) .gt. 0.) mask = .true.

   ! call model%set_bottom_index with a 2D array 
   ! call model%set_bottom_index(indend) ! ltief, lacz !!!

#if 0
   write(*,*) 'cc1d: '
   write(*,*) cc1d(:,1)
   write(*,*) 'cc: '
   write(*,*) cc(:,:,1,1)
   write(*,*) 'mask: '
   write(*,*) mask(:,:,1)
   stop
#endif

   return
end subroutine init_hamsom_fabm


subroutine allocate_hamsom_fabm(m,n,ilo,ndrei,npel)
   integer, intent(in) :: m,n,ilo ! y (m -> from north), x (n -> from west), z (ilo -> from top)
   integer, intent(in) :: ndrei
   integer, intent(in) :: npel

   allocate(cc1d(ndrei,npel),stat=stat)
   if (stat /= 0) stop 'init_hamsom_fabm(): Error allocating memory (cc1d)'

   allocate(cc(m,n,ilo,npel),stat=stat)
   if (stat /= 0) stop 'init_hamsom_fabm(): Error allocating memory (cc)'

   allocate(mask(m,n,ilo),stat=stat)
   if (stat /= 0) stop 'init_hamsom_fabm(): Error allocating memory (mask)'

   allocate(S(m,n,ilo),stat=stat)
   if (stat /= 0) stop 'init_hamsom_fabm(): Error allocating memory (S)'

   allocate(T(m,n,ilo),stat=stat)
   if (stat /= 0) stop 'init_hamsom_fabm(): Error allocating memory (T)'

   return
end subroutine allocate_hamsom_fabm


subroutine do_hamsom_fabm(ndrei,npel)
   integer, intent(in) :: ndrei,npel
   integer :: l
   write(*,*) 'Inside do_hamsom_fabm()'

   cc1d = cc1d+0.1 ! effect of advection/diffusion etc.

   do l=1,npel
      call deco1d3d(cc(:,:,:,l),cc1d(:,l),ndrei)
   end do

   ! here cc should be updated

   do l=1,npel
      call comp3d1d(cc(:,:,:,l),cc1d(:,l))
   end do

   return
end subroutine do_hamsom_fabm

subroutine clean_hamsom_fabm()
   write(*,*) 'Inside clean_hamsom_fabm()'
end subroutine clean_hamsom_fabm

#endif
end module hamsom_fabm
