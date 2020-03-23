module hamsom_fabm

#ifdef _FABM_

use fabm
use fabm_config
use fabm_types, only: attribute_length, output_none
use fabm_standard_variables, only: type_global_standard_variable
implicit NONE

private

   EXTERNAL comp3d1d

   integer :: stat

   real, allocatable :: cc(:,:,:,:), cc1d(:,:)

   public init_hamsom_fabm, do_hamsom_fabm, clean_hamsom_fabm

contains

subroutine init_hamsom_fabm(m,n,ilo,ndrei,npel)
   integer, intent(in) :: m,n,ilo
   integer, intent(in) :: ndrei
   integer, intent(in) :: npel

   integer :: i,j,k,l

   write(*,*) 'Inside init_hamsom_fabm()'
   write(*,*) '  m,n,ilo= ',m,n,ilo
   write(*,*) '  ndrei=   ',ndrei
   write(*,*) '  npel=    ',npel

   allocate(cc(m,n,ilo,npel),stat=stat)
   if (stat /= 0) stop 'init_hamsom_fabm(): Error allocating memory (cc)'
   do l=1,npel
      do k=1,ilo
         do j=1,n
            do i=1,m
               cc(i,j,k,l) = 100000*i+10000*j+100*k+l
            end do
         end do
      end do
   end do
   if (stat /= 0) stop 'init_hamsom_fabm(): Error allocating memory (cc1d)'
   allocate(cc1d(ndrei,npel),stat=stat)
   do l=1,npel
      call comp3d1d(cc(:,:,:,l),cc1d(:,l))
   end do
#if 0
   write(*,*) 'cc', cc
   write(*,*) 'cc1d',cc1d
   stop 'kaj'
#endif
   return
end subroutine init_hamsom_fabm

subroutine do_hamsom_fabm()
   write(*,*) 'Inside do_hamsom_fabm()'
end subroutine do_hamsom_fabm

subroutine clean_hamsom_fabm()
   write(*,*) 'Inside clean_hamsom_fabm()'
end subroutine clean_hamsom_fabm

#endif
end module hamsom_fabm
