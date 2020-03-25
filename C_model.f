*     c_model.f

c______ Model north_b _________________
      parameter (m=177,n=207,ilo=20,nxra=m*n)                  ! i,j,k model size
      parameter (khor=8216,ndrei=82108,kasor=8106)    ! 2d->1d, 3d->1d, 2d->1d for SOR
      parameter (lrp=55,iranz=2,lvrp=840)             ! liquid boudaries size
c      parameter (nbio=23)  !biological parameters, see C_bio.f
c      parameter (nbio=15)  !biological parameters, see C_bio.f
      parameter (nbio=2)  !biological parameters, see C_bio.f
c Attention!!! ngro must be recalculated depending on model domain geometry by
c              ngro = (max(m*(ilo*20+10),(kasor*8)))
      parameter (ngro=m*(ilo*20+10)) ! north_b
c_____   north_b _________________
