module hamsom_fabm

#ifdef _FABM_

use fabm
use fabm_config
use fabm_types, only: attribute_length, output_none
use fabm_standard_variables, only: type_global_standard_variable
implicit NONE

private
   integer :: stat

   real, allocatable :: cc(:,:,:,:), cc1d(:,:)

   public init_hamsom_fabm, do_hamsom_fabm, clean_hamsom_fabm

contains

subroutine init_hamsom_fabm(m,n,ilo,ndrei,npel)
   integer, intent(in) :: m,n,ilo
   integer, intent(in) :: ndrei
   integer, intent(in) :: npel
   write(*,*) 'Inside init_hamsom_fabm()'
   write(*,*) '  ndrei= ',ndrei
   allocate(cc(m,n,ilo,npel),stat=stat)
   allocate(cc1d(ndrei,npel),stat=stat)
end subroutine init_hamsom_fabm

subroutine do_hamsom_fabm()
   write(*,*) 'Inside do_hamsom_fabm()'
end subroutine do_hamsom_fabm

subroutine clean_hamsom_fabm()
   write(*,*) 'Inside clean_hamsom_fabm()'
end subroutine clean_hamsom_fabm

#endif
end module hamsom_fabm
