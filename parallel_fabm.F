C******************* INITIALIZE MPI *******************************
      subroutine init_parallel(myid,nprocs,ierr)
      include "mpif.h"
      integer :: nprocs_test
      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, myid,ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs_test, ierr)
      if (nprocs_test .ne. nprocs) then
        if (myid .eq. 0) write(*,'(/a80,i5/)') "ERROR: Number of "//
     &     "runtime PEs does not fit nprocs set in C_mpi.f: nprocs =",
     &      nprocs
        call MPI_ABORT(MPI_COMM_WORLD, ierr)
      endif
      return
      end

C******************* Load Balancing ********************************      

      subroutine domain_decomp(lazc,iwet)
      include 'C_model.f'
      include 'C_index.f'
      include 'C_mpi.f'
      PARAMETER (khor1=khor+1, nprocs1=nprocs+1)
      common /cord1/  icord(nprocs,4)
      dimension Nblk(nprocs1),lazc(khor),iwet(khor1)
      dimension lwet(n),ntwet(m,n),nvcell(nprocs1),nhcell(m)
      integer iim
      common /bcor/ JWB1(2),JWB2(2),IWB1(2),IWB2(2),NWB
     &,JNB1(2),JNB2(2),INB1(2),INB2(2),NNB
      icord=0
      lwet=0
      ntwet=0
      nvcell=0
      nhcell=0
      Nblk=0
      
      NNB=1
c     Northen boundary      
      INB1(1)=66
      INB2(1)=68
      JNB1(1)=12
      JNB2(1)=57

      NWB=1
c     Eastern Boundary 
      IWB1(1)=156
      IWB2(1)=164
      JWB1(1)=16
      JWB2(1)=18

      isur_con=0
      iim=1
      nblk_vproc=(ndrei+isur_con*khor)/nPv
      n_blk=nblk_vproc
      np=1
      Nblk(np)=0
      nvcell(np)=0
      np=2

      DO J=1,n
 8    continue
      lwe=le(J)
      if(lwe .gt. 0) then
      lwet(J)=indwet(lwe)+lazc(lwe)+lwe*isur_con
      IF(lwet(J) .GE. n_blk)THEN 
      IF( JWB1(iim) .GE. Nblk(np)+1 .AND. JWB1(iim) .LE. J) THEN
      icol1=J-JWB1(iim)
      IF(icol1 .GE. 2) THEN
      Nblk(np)=J
      nvcell(np)=indwet(lwe)+lazc(lwe)
      np=np+1
      n_blk=n_blk+nblk_vproc
      iim=iim+1
      ELSE IF(icol1 .EQ. 0) THEN
      Nblk(np)=J-1
      lwe1=le(j-1)
      nvcell(np)=indwet(lwe1)+lazc(lwe1)
      np=np+1
      n_blk=n_blk+nblk_vproc
      iim=iim+1
      goto 8
      ENDIF
      ELSE
      nextra1=abs(lwet(J)-n_blk)
      nextra2=abs(lwet(J-1)-n_blk)
      if(nextra1 .le. nextra2)then
      Nblk(np)=J
      nvcell(np)=indwet(lwe)+lazc(lwe)
      np=np+1
      n_blk=lwet(J)+nblk_vproc
      else
      Nblk(np)=J-1
      lwe1=le(j-1)
      nvcell(np)=indwet(lwe1)+lazc(lwe1)
      np=np+1
      n_blk=lwet(J-1)+nblk_vproc
      goto 8
      ENDIF
      ENDIF
      ENDIF
      ENDIF
      ENDDO
     
      Nblk(nPv+1)=n
      nvcell(nPv+1)=ndrei
      ntwet=0
      DO j=1,n
      lwa=lb(j)
      lwe=le(j)
      DO lw=lwa,lwe
      i=iwet(lw)
      ntwet(i,j)=lazc(lw)
      ENDDO
      ENDDO
      
      DO nv=1,nPv
      nblk_cell=(nvcell(nv+1)-nvcell(nv))/nPh
      ncell=0
      DO I=1,M
      DO J=Nblk(nv)+1,Nblk(nv+1)
      ncell=ncell+ntwet(I,J)
      ENDDO
      nhcell(I)=ncell
      ENDDO
      
      DO nh=1,nPh
      icord(nh+(nv-1)*nPh,3)=Nblk(nv)+1         !J1
      icord(nh+(nv-1)*nPh,4)=Nblk(nv+1)       !J2
      ENDDO
      
      icell=0
      nh=1
      nblk_cell1=nblk_cell
      DO I=1,M
      icell=icell+nhcell(I)
 11   continue
      If(nhcell(I) .GE. nblk_cell1 )THEN
      iextra1=abs(nblk_cell1-nhcell(I))
      iextra2=abs(nblk_cell1-nhcell(I-1))
      if(iextra1 .le. iextra2)then
      icord(nh+(nv-1)*nPh,2)=I       !I2
c      nblk_cell1=nblk_cell1+nblk_cell
      nblk_cell1=nhcell(I)+nblk_cell
      nh=nh+1
      else
      icord(nh+(nv-1)*nPh,2)=I-1
c      nblk_cell1=nblk_cell1+nblk_cell
      nblk_cell1=nhcell(I-1)+nblk_cell
      nh=nh+1
      goto 11
      endif
      ENDIF
      ENDDO
      
      icord(nv*nPh,2)=M
      icord((nv-1)*nPh+1,1)=1
      ENDDO
      
      DO nv=1,nPv
      DO nh=2,nPh
      icord(nh+(nv-1)*nPh,1)=icord(nh-1+(nv-1)*nPh,2)+1
      ENDDO
      ENDDO

      DO I=1,NNB
      DO J=1,nPh*nPv
      JK1=INB1(I)
      JK2=INB2(I)
      JK3=JNB1(I)
      JK4=JNB2(I)
      IF(ICORD(J,3) .GT. JK4 .OR. ICORD(J,4) .LT. JK3 
     &.OR. ICORD(J,2) .LT. JK1 .OR. ICORD(J,1) .GT. JK2) goto 7

      IF(JK1 .GE. ICORD(J,1) .AND. JK1 .LE. ICORD(J,2)) then
      imt=ICORD(J,2)-JK1
      if(imt .eq. 1) then
      ICORD(J,2)=ICORD(J,2)+1
      ICORD(J+1,1)=ICORD(J+1,1)+1
      else if (imt .eq. 0) then
      ICORD(J,2)=ICORD(J,2)-1
      ICORD(J+1,1)=ICORD(J+1,1)-1
      endif
      ENDIF
 7    continue       
      ENDDO
      ENDDO
      DO I=1,nPh*nPv
      msub=ICORD(I,2)-ICORD(I,1)+1
      nsub=ICORD(I,4)-ICORD(I,3)+1
      if(msub .lt. 10 .and. nsub .lt. 10) then
      write(*,*)'Reduce nPh and nPv'
c      stop
      else if(msub .lt. 10) then
      write(*,*)'Reduce nPh'
c      stop
      else if(nsub .lt. 8) then
      write(*,*)'Reduce nPv' 
c      stop
      endif
      ENDDO
      return
      end 
C**********Defining INDEX for Each Processor*****************************************
      subroutine INDEX_MPI(icord)
      include 'C_mpi.f'
      include 'C_model.f'
      dimension icord(nprocs,4)
       common /cord/ I1,I2,J1,J2,IA1,IA2,JA1,JA2,IB1,IB2,JB1,JB2
     1,ID1,ID2,JD1,JD2
      common /icecor/ ICEI1,ICEI2,ICEJ1,ICEJ2,ICEIA1,ICEIA2,ICEJA1,
     * ICEJA2,ICEX1,ICEX2,ICEY1,ICEY2,ICEP1,ICEP2,ICEQ1,ICEQ2,ICEX3,
     * ICEX4,ICEY3,ICEY4,ICEX5,ICEX6,ICEY5,ICEY6
      I1=icord(myid+1,1)
      I2=icord(myid+1,2)
      J1=icord(myid+1,3)
      J2=icord(myid+1,4)
      if(mod(myid,nPh) .EQ. 0)THEN
      IA1=I1
      IA2=I2+1
      IB1=I1
      IB2=I2+2 
      ICEI1=2
      ICEI2=I2   
      ICEIA1=2
      ICEIA2=I2+1
      ICEX1=2
      ICEX2=I2
      ICEX3=2
      ICEX4=I2+1
      ICEX5=2
      ICEX6=I2+2
      ICEP1=3
      ICEP2=I2
      ID1=1
      ID2=I2+4
      ELSE IF(mod(myid+1,nPh) .EQ. 0)THEN
      IA1=I1-1
      IA2=I2
      IB1=I1-2
      IB2=I2
      ICEI1=I1
      ICEI2=I2-1
      ICEIA1=I1-1
      ICEIA2=I2-1
      ID1=I1-4
      ID2=I2
      ICEX1=I1
      ICEX2=I2
      ICEX3=I1-1
      ICEX4=I2
      ICEX5=I1-2
      ICEX6=I2
      ICEP1=I1
      ICEP2=I2-1
      ELSE
      IA1=I1-1
      IA2=I2+1
      IB1=I1-2
      IB2=I2+2
      ICEI1=I1
      ICEI2=I2
      ID1=I1-4
      ID2=I2+4
      ICEIA1=I1-1
      ICEIA2=I2+1
      ICEX1=I1
      ICEX2=I2
      ICEX3=I1-1
      ICEX4=I2+1
      ICEX5=I1-2
      ICEX6=I2+2
      ICEP1=I1
      ICEP2=I2
      ENDIF
      
      if(myid .LT. nPh) THEN
      JA1=J1
      JA2=J2+1
      JB1=J1
      JB2=J2+2
      ICEJ1=2
      ICEJ2=J2
      ICEJA1=2
      ICEJA2=J2+1
      JD1=J1
      JD2=J2+4
      ICEY1=1
      ICEY2=J2
      ICEY3=1
      ICEY4=J2+1
      ICEY5=1
      ICEY6=J2+2
      ICEQ1=2
      ICEQ2=J2
      else if(myid .GT. nPh*(nPv-1)-1) THEN
      JA1=J1-1
      JA2=J2
      JB1=J1-2
      JB2=J2
      ICEJ1=J1
      ICEJ2=J2-1
      ICEJA1=J1-1
      ICEJA2=J2-1
      JD1=J1-4
      JD2=J2
      ICEY1=J1
      ICEY2=J2-1
      ICEY3=J1-1
      ICEY4=J2-1
      ICEY5=J1-2
      ICEY6=J2-1
      ICEQ1=J1
      ICEQ2=J2-2
      else
      JA1=J1-1
      JA2=J2+1
      JB1=J1-2
      JB2=J2+2
      ICEJ1=J1
      ICEJ2=J2
      ICEJA1=J1-1
      ICEJA2=J2+1
      JD1=J1-4
      JD2=J2+4
      ICEY1=J1
      ICEY2=J2
      ICEY3=J1-1
      ICEY4=J2+1
      ICEY5=J1-2
      ICEY6=J2+2
      ICEQ1=J1
      ICEQ2=J2
      endif
cctest fic IA2=178 when nproc =1 Ute Daewel
!      if(nprocs.eq.1)then
!      IA2=min(IA2,m)
!      IB2=min(IB2,m)
!      JB2=min(JB2,n)
!      endif
c      WRITE(*,*) I1,I2,J1,J2,IA1,IA2,JA1,JA2,IB1,IB2,JB1,JB2,IN1,IN2,
c     * JN1,JN2,IM1,IM2,JM1,JM2,IT1,IT2,JT1,JT2,IQ1,IQ2,JQ1,JQ2,ICEI1,
c     * ICEI2,ICEJ1,ICEJ2
      return
      end  

**************************** VERTICAL COMMUNICATION ELEMENTS************************************************************************************************
      subroutine ver_comm_cord(lazc,izet,icord)
      include 'C_model.f'
      include 'C_mpi.f'
      common /vcord/ nvrlayer(4),nvllayer(4),nvrvar(4*n,3),nvlvar(4*n,3)
     *,nvrdep(4),nvldep(4)
      dimension icord(nprocs,4),izet(m,n),lazc(khor)      
      nvrlayer=0
      nvllayer=0
      nvrvar=0
      nvlvar=0
      nvrdep=0
      nvldep=0 
      I1=icord(myid+1,1)
      I2=icord(myid+1,2)
      J1=icord(myid+1,3)
      J2=icord(myid+1,4)
c      WRITE(*,*)I1,I2,J1,J2 
      if(mod((myid+1),nPh) .ne. 0)THEN
      it=0
      it1=0
      DO I=I2-1,I2+2
      k=0
      kdep=0
      DO J=J1,J2
      if(izet(I,J) .ne. 0) THEN
      k=k+1
      nvrvar(it1+k,1)=izet(I,J)
      nvrvar(it1+k,2)=I
      nvrvar(it1+k,3)=J
      kdep=kdep+lazc(izet(I,J))
      ENDIF
      ENDDO
      it=it+1
      it1=it1+k
      nvrlayer(it)=k
      nvrdep(it)=kdep
      ENDDO
      ENDIF

      if(mod(myid,nPh) .ne. 0)THEN
      it=0
      it1=0
      DO I=I1-2,I1+1
      k=0
      kdep=0
      DO J=J1,J2
      if(izet(I,J) .ne. 0) THEN
      k=k+1
      nvlvar(it1+k,1)=izet(I,J)
      nvlvar(it1+k,2)=I
      nvlvar(it1+k,3)=J
      kdep=kdep+lazc(izet(I,J))
      ENDIF
      ENDDO
      it=it+1
      it1=it1+k
      nvllayer(it)=k
      nvldep(it)=kdep
      ENDDO
      ENDIF
c      DO I=1,4
c      WRITE(*,*) nvldep(I),nvrdep(I),nvllayer(I),nvrlayer(I),myid
c      ENDDO
      return
      end
*********************Exact START and END cell for the subdomain******************************************
      subroutine temp_izet(izet)
      include 'C_model.f'
      dimension izet(m,n)
      common /surface/ izets(m,n),izete(m,n)
      izets=0
      izete=0
      it=1
      DO 6 J=1,n
      DO 6 I=1,m
      izets(i,j)=it
      if(izet(i,j).ne.0) it=it+1
      izete(i,j)=it-1
 6    continue
      return
      end
*********************************Inter connectivity of negihbours domain X direction****************************************************************
c where 'b' represents datas are send from bottom to top layer
c '22' represents if domain is varying from I1 to I2 
c then we are looking for connectivity I1-2 to I1+2 
      subroutine hor_comm_b22(icord,izet,lazc)
      include 'C_model.f'
      include 'C_index.f'
      include 'C_mpi.f'
      dimension ifst(nprocs),isec(nprocs),izet(m,n)
      dimension lazc(khor),icord(nprocs,4)
      common /surface/ izets(m,n),izete(m,n)
      common /hcordb2/ idepb2(nprocs),icdepb2(nprocs),
     & nbhlayer2(nprocs*6,7),nsendb2(nprocs),ncsendb2(nprocs)
      integer i,j,k,I1,I2,I3,I4,id,id1,ibcomm,it,nt,kt,ikt,JS1,JS2,
     &IS1,IS2,lw1,lw2,lw3,lw4,np1,np2 
      idepb2=0
      icdepb2=0
      nbhlayer2=0
      nsendb2=0
      ncsendb2=0
      ifst=0
      isec=0
      k=1
      DO J=1,nPv-1
      DO I=1,nPh
      id=(i-1)+(J-1)*nPh
      I1=icord(id+1,1)
      I2=icord(id+1,2)

      IF(mod(id,nPh) .EQ. 0)THEN
      I1=I1
      I2=I2+2
      ELSE IF(mod(id+1,nPh) .EQ. 0) THEN
      I1=I1-2
      I2=I2
      ELSE
      I1=I1-2
      I2=I2+2
      ENDIF

      k=1
 6    id1=(k-1)+nPh*J
      I3=icord(id1+1,1)
      I4=icord(id1+1,2)
      if(I1 .GE. I3 .AND. I1 .LE. I4 ) THEN
      ifst(id+1)=id1
      else
      k=k+1
      goto 6
      endif

c      k=k-1
c      if(k .le. 0) k=k+1
 7    id1=(k-1)+nPh*J
      I3=icord(id1+1,1)
      I4=icord(id1+1,2)
      IF(I2 .GE. I3 .AND. I2 .LE. I4 ) THEN
      isec(id+1)=id1
      ELSE
      k=k+1
      goto 7
      ENDIF
      ENDDO
      ENDDO

      if(nPh .eq. 1) then
      DO I=1,nPh*(nPv-1)
      ifst(I)=I
      isec(I)=I
      ENDDO
      endif

      ibcomm=0
      DO I=1,nPh*(nPv-1)
      idepb2(I)=isec(I)-ifst(I)+1
      ibcomm=ibcomm+idepb2(I)
      ENDDO
      icdepb2(1)=1
      DO I=2,nPh*(nPv-1)
      icdepb2(I)=idepb2(I-1)+icdepb2(I-1)
      ENDDO 

      it=1
      DO I=1,nPh*(nPv-1)
      DO K=1,idepb2(I)
      nbhlayer2(it,1)=ifst(I)+K-1
      nbhlayer2(it,7)=I-1
      if(K .EQ. 1 .AND. mod(I-1,nPh) .EQ. 0) THEN
      nbhlayer2(it,2)=icord(I,1)
      else IF (K .EQ. 1)THEN
      nbhlayer2(it,2)=icord(I,1)-2
      else
      nbhlayer2(it,2)=icord(ifst(I)+K,1)
      endif
      if(K .EQ. idepb2(I) .AND. mod(I,nPh) .EQ. 0) THEN
      nbhlayer2(it,3)=icord(I,2)
      ELSE if(K .EQ. idepb2(I))THEN
      nbhlayer2(it,3)=icord(I,2)+2
      ELSE
      nbhlayer2(it,3)=icord(ifst(I)+K,2)
      endif
      it=it+1
      ENDDO
      ENDDO

      ikt=1 
      ncsendb2=0
      nsendb2=0
      DO I=1,nPh*(nPv-1)
      j=nPh+I-1
      nt=0
      kt=0
      DO k=1,ibcomm
      IF(nbhlayer2(k,1) .eq. j) THEN
      nbhlayer2(ikt,4)=k
      ikt=ikt+1
      nt=nt+1
      ENDIF
      ENDDO
      nsendb2(I)=nt
      ENDDO

      ncsendb2(1)=1
      DO I=2,nPh*(nPv-1)
      ncsendb2(I)=nsendb2(I-1)+ncsendb2(I-1)
      ENDDO
      DO I=1,ibcomm
      JS1=icord(nbhlayer2(I,1)+1,3)
      JS2=JS1+1
      IS1=nbhlayer2(I,2)
      IS2=nbhlayer2(I,3)
      lw1=izets(IS1,JS1)
      lw2=izete(IS2,JS1)
      lw3=izets(IS1,JS2)
      lw4=izete(IS2,JS2)
      np1=0
      np2=0
      
      if(lw2 .ge. lw1) np1=indwet(lw2)+lazc(lw2)-indwet(lw1)
      if(lw4 .ge. lw3) np2=indwet(lw4)+lazc(lw4)-indwet(lw3)
      nbhlayer2(I,5)=np1
      nbhlayer2(I,6)=np2
      ENDDO
      return
      end
************************************************************************************************************************************************************
      subroutine hor_comm_t22(icord,izet,lazc)
      include 'C_model.f'
      include 'C_index.f' 
      include 'C_mpi.f'
      dimension ifst(nprocs),isec(nprocs),izet(m,n)
      dimension lazc(khor),icord(nprocs,4)
      common /surface/ izets(m,n),izete(m,n)
      common /hcordt2/ idept2(nprocs),icdept2(nprocs),
     &nsendt2(nprocs),ncsendt2(nprocs),nthlayer2(nprocs*6,7)
      integer i,j,k,I1,I2,I3,I4,id,id1,ibcomm,it,nt,kt,itk,JS1,JS2,
     &IS1,IS2,lw1,lw2,lw3,lw4,np1,np2

      idept2=0
      icdept2=0
      nthlayer2=0
      nsendt2=0
      ncsendt2=0

      ifst=0
      isec=0     

      k=1

      DO J=1,nPv-1
      DO I=1,nPh
      id=(i-1)+J*nPh
      I1=icord(id+1,1)
      I2=icord(id+1,2)

      IF(mod(id,nPh) .EQ. 0)THEN
      I1=I1
      I2=I2+2
      ELSE IF(mod(id+1,nPh) .EQ. 0) THEN
      I1=I1-2
      I2=I2
      ELSE
      I1=I1-2
      I2=I2+2
      ENDIF

      k=1
 6    id1=k-1+(J-1)*nPh
      I3=icord(id1+1,1)
      I4=icord(id1+1,2)
      if(I1 .GE. I3 .AND. I1 .LE. I4 ) THEN
      ifst(id+1-nPh)=id1
      else
      k=k+1
      goto 6
      ENDIF
      k=1
 7    id1=k-1+(J-1)*nPh
      I3=icord(id1+1,1)
      I4=icord(id1+1,2)
      if(I2 .GE. I3 .AND. I2 .LE. I4 ) THEN
      isec(id+1-nPh)=id1
      else
      k=k+1
      goto 7
      ENDIF
      ENDDO
      ENDDO
 
      if(nPh .eq. 1) then
      DO I=1,nPh*(nPv-1)
      ifst(I)=I-1
      isec(I)=I-1
      ENDDO
      endif

      ibcomm=0
      DO I=1,nPh*(nPv-1)
      idept2(I)=isec(I)-ifst(I)+1
      ibcomm=ibcomm+idept2(I)
      ENDDO
      icdept2(1)=1
      DO I=2,nPh*(nPv-1)
      icdept2(I)=idept2(I-1)+icdept2(I-1)
      ENDDO 

      it=1
      DO I=nPh+1,nPh*nPv
      
      DO K=1,idept2(I-nPh)
      nthlayer2(it,1)=ifst(I-nPh)+K-1
      nthlayer2(it,7)=I-1
      if(K .EQ. 1 .AND. mod(I-1,nPh) .EQ. 0) THEN
      nthlayer2(it,2)=icord(I,1)
      
      else IF (K .EQ. 1)THEN
      nthlayer2(it,2)=icord(I,1)-2
      else
      nthlayer2(it,2)=icord(ifst(I-nPh)+K,1)
      endif
      if(K .EQ. idept2(I-nPh) .AND. mod(I,nPh) .EQ. 0) THEN
      nthlayer2(it,3)=icord(I,2)
      ELSE if(K .EQ. idept2(I-nPh))THEN
      nthlayer2(it,3)=icord(I,2)+2
      ELSE
      nthlayer2(it,3)=icord(ifst(I-nPh)+K,2)
      endif
      it=it+1
      ENDDO
      ENDDO

      itk=1
      DO I=1,nPh*(nPv-1)
      J=I-1
      nt=0
      DO k=1,ibcomm
      IF(nthlayer2(k,1) .EQ. J) THEN
      nthlayer2(itk,4)=k
      itk=itk+1
      nt=nt+1
      ENDIF
      ENDDO
      nsendt2(I)=nt
      ENDDO

      ncsendt2(1)=1
      DO I=2,nPh*(nPv-1)
      ncsendt2(I)=nsendt2(I-1)+ncsendt2(I-1)
      ENDDO
      DO I=1,ibcomm
      JS1=icord(nthlayer2(I,1)+1,4)
      JS2=JS1-1
      IS1=nthlayer2(I,2)
      IS2=nthlayer2(I,3)
      lw1=izets(IS1,JS1)
      lw2=izete(IS2,JS1)
      lw3=izets(IS1,JS2)
      lw4=izete(IS2,JS2)
      np1=0
      np2=0
      if(lw2 .ge. lw1) np1=indwet(lw2)+lazc(lw2)-indwet(lw1)
      if(lw4 .ge. lw3) np2=indwet(lw4)+lazc(lw4)-indwet(lw3)
      nthlayer2(I,5)=np1
      nthlayer2(I,6)=np2
      ENDDO
      return
      end
************************************************************************************************************************************************************
      subroutine hor_comm_b11(icord,izet,lazc)
      include 'C_model.f'
      include 'C_index.f'
      include 'C_mpi.f'
      dimension ifst(nprocs),isec(nprocs),izet(m,n)
      dimension lazc(khor),icord(nprocs,4)
      common /surface/ izets(m,n),izete(m,n)
      common /hcordb1/ idepb1(nprocs),icdepb1(nprocs),
     &nbhlayer1(nprocs*6,7),nsendb1(nprocs),ncsendb1(nprocs)

      idepb1=0
      icdepb1=0
      nbhlayer1=0
      nsendb1=0
      ncsendb1=0
      ifst=0
      isec=0     
      k=1
      DO J=1,nPv-1
      DO I=1,nPh
      id=(i-1)+(J-1)*nPh
      I1=icord(id+1,1)
      I2=icord(id+1,2)
      IF(mod(id,nPh) .EQ. 0)THEN
      I1=I1
      I2=I2+1
      ELSE IF(mod(id+1,nPh) .EQ. 0) THEN
      I1=I1-1
      I2=I2
      ELSE
      I1=I1-1
      I2=I2+1
      ENDIF
      k=1
c      if(k .le. 0) k=k+1
 6    id1=(k-1)+nPh*J
      I3=icord(id1+1,1)
      I4=icord(id1+1,2)
      if(I1 .GE. I3 .AND. I1 .LE. I4 ) THEN
      ifst(id+1)=id1
      else
      k=k+1
      goto 6
      ENDIF
      k=k-1
      if(k .le. 0) k=1
 7    id1=(k-1)+nPh*J
      I3=icord(id1+1,1)
      I4=icord(id1+1,2)
      if(I2 .GE. I3 .AND. I2 .LE. I4 ) THEN
      isec(id+1)=id1
      else
      k=k+1
      goto 7
      ENDIF
      ENDDO
      ENDDO
      if(nPh .eq. 1) then
      DO I=1,nPh*(nPv-1)
      ifst(I)=I
      isec(I)=I
      ENDDO
      endif

      ibcomm=0
      DO I=1,nPh*(nPv-1)
      idepb1(I)=isec(I)-ifst(I)+1
      ibcomm=ibcomm+idepb1(I)
      ENDDO
      icdepb1(1)=1
      DO I=2,nPh*(nPv-1)
      icdepb1(I)=idepb1(I-1)+icdepb1(I-1)
      ENDDO 

      it=1

      DO I=1,nPh*(nPv-1)
      
      DO K=1,idepb1(I)
      nbhlayer1(it,1)=ifst(I)+K-1
      nbhlayer1(it,7)=I-1
      if(K .EQ. 1 .AND. mod(I-1,nPh) .EQ. 0) THEN
      nbhlayer1(it,2)=icord(I,1)
      
      else IF (K .EQ. 1)THEN
      nbhlayer1(it,2)=icord(I,1)-1
      else
      nbhlayer1(it,2)=icord(ifst(I)+K,1)
      endif
      if(K .EQ. idepb1(I) .AND. mod(I,nPh) .EQ. 0) THEN
      nbhlayer1(it,3)=icord(I,2)
      ELSE if(K .EQ. idepb1(I))THEN
      nbhlayer1(it,3)=icord(I,2)+1
      ELSE
      nbhlayer1(it,3)=icord(ifst(I)+K,2)
      endif
      it=it+1
      ENDDO
      ENDDO
      ncsendb1=0
      nsendb1=0 
      it=1
      DO I=1,nPh*(nPv-1)
      J=nPh+I-1
      nt=0
      kt=0
      DO k=1,ibcomm
      IF(nbhlayer1(k,1) .EQ. J) THEN
      nbhlayer1(it,4)=k
      it=it+1
      nt=nt+1
c      ELSE
c      kt=kt+1
c      if(kt .GT. 2) goto 4
      ENDIF
      ENDDO
 4    continue
      nsendb1(I)=nt
      ENDDO
      ncsendb1(1)=1
      DO I=2,nPh*(nPv-1)
      ncsendb1(I)=nsendb1(I-1)+ncsendb1(I-1)
      ENDDO
      
      DO I=1,ibcomm
      JS1=icord(nbhlayer1(I,1)+1,3)
      JS2=JS1+1
      IS1=nbhlayer1(I,2)
      IS2=nbhlayer1(I,3)
      lw1=izets(IS1,JS1)
      lw2=izete(IS2,JS1)
      lw3=izets(IS1,JS2)
      lw4=izete(IS2,JS2)
      np1=indwet(lw2)+lazc(lw2)-indwet(lw1)
      np2=indwet(lw4)+lazc(lw4)-indwet(lw3)
      if(np1 .lt. 0)np1=0
      if(np2 .lt. 0)np2=0
      nbhlayer1(I,5)=np1
      nbhlayer1(I,6)=np2
c      IF(myid .EQ. 0)  write(*,*)np1,np2,I,IS1,IS2,Js1,Js2
      ENDDO
      
c      write(*,*)ibcomm,'nbhlayer1',myid
      return
      end
************************************************************************************************************************************************************
      subroutine hor_comm_t11(icord,izet,lazc)
      include 'C_model.f'
      include 'C_index.f'
      include 'C_mpi.f'
      dimension ifst(nprocs),isec(nprocs),izet(m,n)
      dimension lazc(khor),icord(nprocs,4)
      common /surface/ izets(m,n),izete(m,n)
      common /hcordt1/ idept1(nprocs),icdept1(nprocs),
     &nthlayer1(nprocs*6,7),nsendt1(nprocs),ncsendt1(nprocs)
      idept1=0
      icdept1=0
      nthlayer1=0
      nsendt1=0
      ncsendt1=0
      ifst=0
      isec=0
      k=1
      DO J=1,nPv-1
      DO I=1,nPh
      id=(i-1)+J*nPh
      I1=icord(id+1,1)
      I2=icord(id+1,2)
      IF(mod(id,nPh) .EQ. 0)THEN
      I1=I1
      I2=I2+1
      ELSE IF(mod(id+1,nPh) .EQ. 0) THEN
      I1=I1-1
      I2=I2
      ELSE
      I1=I1-1
      I2=I2+1
      ENDIF
      k=1
c      if(k .le. 0) k=k+1
 6    id1=k-1+nPh*(J-1)
      I3=icord(id1+1,1)
      I4=icord(id1+1,2)
      if(I1 .GE. I3 .AND. I1 .LE. I4 ) THEN
      ifst(id+1-nPh)=id1
      else
      k=k+1
      goto 6
      ENDIF
      if(k .le. 0) k=1
 7    id1=k-1+nPh*(J-1)
      I3=icord(id1+1,1)
      I4=icord(id1+1,2)
      if(I2 .GE. I3 .AND. I2 .LE. I4 ) THEN
      isec(id+1-nPh)=id1
      else
      k=k+1
      goto 7
      ENDIF
      ENDDO
      ENDDO
      if(nPh .eq. 1) then
      DO I=1,nPh*(nPv-1)
      ifst(I)=I-1
      isec(I)=I-1
      ENDDO
      endif
      
      ibcomm=0
      DO I=1,nPh*(nPv-1)
      idept1(I)=isec(I)-ifst(I)+1
      ibcomm=ibcomm+idept1(I)
      ENDDO
      icdept1(1)=1
      DO I=2,nPh*(nPv-1)
      icdept1(I)=idept1(I-1)+icdept1(I-1)
      ENDDO 
      it=1
      DO I=nPh+1,nPh*nPv
      
      DO K=1,idept1(I-nPh)
      nthlayer1(it,1)=ifst(I-nPh)+K-1
      nthlayer1(it,7)=I-1
      if(K .EQ. 1 .AND. mod(I-1,nPh) .EQ. 0) THEN
      nthlayer1(it,2)=icord(I,1)
      
      else IF (K .EQ. 1)THEN
      nthlayer1(it,2)=icord(I,1)-1
      else
      nthlayer1(it,2)=icord(ifst(I-nPh)+K,1)
      endif
      if(K .EQ. idept1(I-nPh) .AND. mod(I,nPh) .EQ. 0) THEN
      nthlayer1(it,3)=icord(I,2)
      ELSE if(K .EQ. idept1(I-nPh))THEN
      nthlayer1(it,3)=icord(I,2)+1
      ELSE
      nthlayer1(it,3)=icord(ifst(I-nPh)+K,2)
      endif
      it=it+1
      ENDDO
      ENDDO
      ncsendt1=0
      nsendt1=0
      it=1
      DO I=1,nPh*(nPv-1)
      J=I-1
      nt=0
      kt=0
      DO k=1,ibcomm
      IF(nthlayer1(k,1) .EQ. J) THEN
      nthlayer1(it,4)=k
      it=it+1
      nt=nt+1
c      ELSE
c      kt=kt+1
c      if(kt .GT. 2) goto 4
      ENDIF
      ENDDO
 4    continue
      nsendt1(I)=nt
      ENDDO
      ncsendt1(1)=1
      DO I=2,nPh*(nPv-1)
      ncsendt1(I)=nsendt1(I-1)+ncsendt1(I-1)
      ENDDO
      
      DO I=1,ibcomm
      JS1=icord(nthlayer1(I,1)+1,4)
      JS2=JS1-1
      IS1=nthlayer1(I,2)
      IS2=nthlayer1(I,3)
      lw1=izets(IS1,JS1)
      lw2=izete(IS2,JS1)
      lw3=izets(IS1,JS2)
      lw4=izete(IS2,JS2)
      np1=indwet(lw2)+lazc(lw2)-indwet(lw1)
      np2=indwet(lw4)+lazc(lw4)-indwet(lw3)
      if(np1 .lt. 0)np1=0
      if(np2 .lt. 0)np2=0
      nthlayer1(I,5)=np1
      nthlayer1(I,6)=np2
c      IF(myid .EQ. 0)  write(*,*)np1,np2,I,IS1,IS2,Js1,Js2
      ENDDO
c      write(*,*)ibcomm,'nthlayer1',myid
      return
      end
************************************************************************************************************************************************************
      subroutine hor_comm_b00(icord,izet,lazc)
      include 'C_model.f'
      include 'C_index.f'
      include 'C_mpi.f'
      dimension ifst(nprocs),isec(nprocs),izet(m,n)
      dimension lazc(khor),icord(nprocs,4)
      common /surface/ izets(m,n),izete(m,n)
      common /hcordb0/ idepb0(nprocs),icdepb0(nprocs),
     &nbhlayer0(nprocs*6,7),nsendb0(nprocs),ncsendb0(nprocs)

      idepb0=0
      icdepb0=0
      nbhlayer0=0
      nsendb0=0
      ncsendb0=0
      ifst=0
      isec=0     
      k=1
      DO J=1,nPv-1
      DO I=1,nPh
      id=(i-1)+(J-1)*nPh
      I1=icord(id+1,1)
      I2=icord(id+1,2)
      k=1
c      if(k .le. 0) k=k+1
 6    id1=(k-1)+nPh*J
      I3=icord(id1+1,1)
      I4=icord(id1+1,2)
      if(I1 .GE. I3 .AND. I1 .LE. I4 ) THEN
      ifst(id+1)=id1
      else
      k=k+1
      goto 6
      ENDIF
      k=k-1
      if(k .le. 0) k=1
 7    id1=(k-1)+J*nPh
      I3=icord(id1+1,1)
      I4=icord(id1+1,2)
      if(I2 .GE. I3 .AND. I2 .LE. I4 ) THEN
      isec(id+1)=id1
      else
      k=k+1
      goto 7
      ENDIF
      ENDDO
      ENDDO
      if(nPh .eq. 1) then
      DO I=1,nPh*(nPv-1)
      ifst(I)=I
      isec(I)=I
      ENDDO
      endif

      ibcomm=0
      DO I=1,nPh*(nPv-1)
      idepb0(I)=isec(I)-ifst(I)+1
      ibcomm=ibcomm+idepb0(I)
      ENDDO
      icdepb0(1)=1
      DO I=2,nPh*(nPv-1)
      icdepb0(I)=idepb0(I-1)+icdepb0(I-1)
      ENDDO 

      it=1

      DO I=1,nPh*(nPv-1)
      
      DO K=1,idepb0(I)
      nbhlayer0(it,1)=ifst(I)+K-1
      nbhlayer0(it,7)=I-1
      if(K .EQ. 1 .AND. mod(I-1,nPh) .EQ. 0) THEN
      nbhlayer0(it,2)=icord(I,1)
      
      else IF (K .EQ. 1)THEN
      nbhlayer0(it,2)=icord(I,1)
      else
      nbhlayer0(it,2)=icord(ifst(I)+K,1)
      endif
      if(K .EQ. idepb0(I) .AND. mod(I,nPh) .EQ. 0) THEN
      nbhlayer0(it,3)=icord(I,2)
      ELSE if(K .EQ. idepb0(I))THEN
      nbhlayer0(it,3)=icord(I,2)
      ELSE
      nbhlayer0(it,3)=icord(ifst(I)+K,2)
      endif
      it=it+1
      ENDDO
      ENDDO
      ncsendb0=0
      nsendb0=0
      it=1
      DO I=1,nPh*(nPv-1)
      J=nPh+I-1
      nt=0
      kt=0
      DO k=1,ibcomm
      IF(nbhlayer0(k,1) .EQ. J) THEN
      nbhlayer0(it,4)=k
      it=it+1
      nt=nt+1
      ENDIF
      ENDDO
 4    continue
      nsendb0(I)=nt
      ENDDO
      ncsendb0(1)=1
      DO I=2,nPh*(nPv-1)
      ncsendb0(I)=nsendb0(I-1)+ncsendb0(I-1)
      ENDDO
      
      DO I=1,ibcomm
      JS1=icord(nbhlayer0(I,1)+1,3)
      JS2=JS1+1
      IS1=nbhlayer0(I,2)
      IS2=nbhlayer0(I,3)
      lw1=izets(IS1,JS1)
      lw2=izete(IS2,JS1)
      lw3=izets(IS1,JS2)
      lw4=izete(IS2,JS2)
      np1=indwet(lw2)+lazc(lw2)-indwet(lw1)
      np2=indwet(lw4)+lazc(lw4)-indwet(lw3)
      if(np1 .lt. 0)np1=0
      if(np2 .lt. 0)np2=0
      nbhlayer0(I,5)=np1
      nbhlayer0(I,6)=np2
      ENDDO
c      write(*,*)ibcomm,'nbhlayer0',myid
      return
      end
************************************************************************************************************************************************************
      subroutine hor_comm_t00(icord,izet,lazc)
      include 'C_model.f'
      include 'C_index.f'
      include 'C_mpi.f'
      dimension ifst(nprocs),isec(nprocs),izet(m,n)
      dimension lazc(khor),icord(nprocs,4)
      common /surface/ izets(m,n),izete(m,n)
      common /hcordt0/ idept0(nprocs),icdept0(nprocs),
     &nthlayer0(nprocs*6,7),nsendt0(nprocs),ncsendt0(nprocs)
      idept0=0
      icdept0=0
      nthlayer0=0
      nsendt0=0
      ncsendt0=0
      ifst=0
      isec=0     
      k=1
      DO J=1,nPv-1
      DO I=1,nPh
      id=(i-1)+J*nPh
      I1=icord(id+1,1)
      I2=icord(id+1,2)
      k=1
      if(k .le. 0) k=1
 6    id1=nPh*(J-1)+k-1
      I3=icord(id1+1,1)
      I4=icord(id1+1,2)
      if(I1 .GE. I3 .AND. I1 .LE. I4 ) THEN
      ifst(id+1-nPh)=id1
      else
      k=k+1
      goto 6
      ENDIF
      if(k .le. 0) k=1
 7    id1=k-1+nPh*(J-1)
      I3=icord(id1+1,1)
      I4=icord(id1+1,2)
      if(I2 .GE. I3 .AND. I2 .LE. I4 ) THEN
      isec(id+1-nPh)=id1
      else
      k=k+1
      goto 7
      ENDIF
      ENDDO
      ENDDO
      if(nPh .eq. 1) then
      DO I=1,nPh*(nPv-1)
      ifst(I)=I-1
      isec(I)=I-1
      ENDDO
      endif
      
      ibcomm=0
      DO I=1,nPh*(nPv-1)
      idept0(I)=isec(I)-ifst(I)+1
      ibcomm=ibcomm+idept0(I)
      ENDDO
      icdept0(1)=1
      DO I=2,nPh*(nPv-1)
      icdept0(I)=idept0(I-1)+icdept0(I-1)
      ENDDO 
      
      it=1
      DO I=nPh+1,nPh*nPv
      
      DO K=1,idept0(I-nPh)
      nthlayer0(it,1)=ifst(I-nPh)+K-1
      nthlayer0(it,7)=I-1
      if(K .EQ. 1 .AND. mod(I-1,nPh) .EQ. 0) THEN
      nthlayer0(it,2)=icord(I,1)
      
      else IF (K .EQ. 1)THEN
      nthlayer0(it,2)=icord(I,1)
      else
      nthlayer0(it,2)=icord(ifst(I-nPh)+K,1)
      endif
      if(K .EQ. idept0(I-nPh) .AND. mod(I,nPh) .EQ. 0) THEN
      nthlayer0(it,3)=icord(I,2)
      ELSE if(K .EQ. idept0(I-nPh))THEN
      nthlayer0(it,3)=icord(I,2)
      ELSE
      nthlayer0(it,3)=icord(ifst(I-nPh)+K,2)
      endif
      it=it+1
      ENDDO
      ENDDO
      ncsendt0=0
      nsendt0=0
      it=1
      DO I=1,nPh*(nPv-1)
      J=I-1
      nt=0
      kt=0
      DO k=1,ibcomm
      IF(nthlayer0(k,1) .EQ. J) THEN
      nthlayer0(it,4)=k
      it=it+1
      nt=nt+1
c      ELSE
c      kt=kt+1
c      if(kt .GT. 2) goto 4
      ENDIF
      ENDDO
 4    continue
      nsendt0(I)=nt
      ENDDO
      ncsendt0(1)=1
      DO I=2,nPh*(nPv-1)
      ncsendt0(I)=nsendt0(I-1)+ncsendt0(I-1)
      ENDDO
      
      DO I=1,ibcomm
      JS1=icord(nthlayer0(I,1)+1,4)
      JS2=JS1-1
      IS1=nthlayer0(I,2)
      IS2=nthlayer0(I,3)
      lw1=izets(IS1,JS1)
      lw2=izete(IS2,JS1)
      lw3=izets(IS1,JS2)
      lw4=izete(IS2,JS2)
      np1=indwet(lw2)+lazc(lw2)-indwet(lw1)
      np2=indwet(lw4)+lazc(lw4)-indwet(lw3)
      if(np1 .lt. 0)np1=0
      if(np2 .lt. 0)np2=0
      nthlayer0(I,5)=np1
      nthlayer0(I,6)=np2
      ENDDO
c      write(*,*)ibcomm,'nthlayer0',myid
      return
      end
*************************************************************************************************************************************************************
      subroutine deco1d2d_s (zf,zc,iwet,indend,fac) 
c-------------------------------------------------------------------
c     creates  2-d arrays from 1-d arrays  
c-----------------------------------------------------------------------izets
      include 'C_model.f'
      parameter (khor1=khor+1)
      dimension iwet(khor1),indend(n) 
      dimension zf(m,n),zc(khor) 
      
      do k=1,n 
      do i=1,m 
         zf(i,k) = 0.0 
      enddo
      enddo

      lwe = 0 
      do k=1,n 
        lwa = lwe+1
        lwe = indend(k)
        do lw=lwa,lwe 
          i  = iwet(lw) 
          zz = zc(lw) 
          zf(i,k) = zz*fac 
        enddo
      enddo
      return 
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
*§§4-5 deco1d3d and comp3d1d

c-----------------------------------------------------------------------
      subroutine deco1d3d_s (u,ucomp,ntot) 
c-----------------------------------------------------------------------
c     expanding 1-d-array ->  3-d-arry 
c-----------------------------------------------------------------------

      include 'C_model.f'
      include 'C_index.f'
      parameter (khor1=khor+1) 
      dimension u(m,n,ilo) 
      dimension ucomp(ntot) 
      common /ind/ iwet(khor1),ldep(khor),lazc(khor), 
     * indend(n),isornr(n),isorsr(n),islab(n) 
      do k=1,n 
         do i=1,m 
            do j=1,ilo 
               u(i,k,j) = -99.0 
            enddo
         enddo
      enddo

      do k=1,n 
         lwa = lb(k)
         lwe = le(k)
         if (lwa.le.lwe) then
            llb = indwet(lwa)+1
            lle = indwet(lwe)+lazc(lwe)
            do ll=llb,lle
               j = indver(ll)
               i = indi(ll)
               u(i,k,j) = ucomp(ll)
            end do
         end if
      end do
      return 
      end

c-----------------------------------------------------------------------
      subroutine comp3d1d_s (d3,d1) 
c-----------------------------------------------------------------------
c    compressing 3-d arrays to compressed 1-d arrays 
c-----------------------------------------------------------------------

      include 'C_model.f'
      include 'C_index.f'
      parameter (khor1=khor+1) 
      dimension d3(m,n,ilo),d1(ndrei) 
      common /ind/ iwet(khor1),ldep(khor),lazc(khor), 
     *  indend(n),isornr(n),isorsr(n),islab(n)
        do k=1,n 
         lwa = lb(k)
         lwe = le(k)
         if (lwa.le.lwe) then
            llb = indwet(lwa)+1
            lle = indwet(lwe)+lazc(lwe)
            do ll=llb,lle
               j = indver(ll)
               i = indi(ll)
               d1(ll) = d3(i,k,j) 
            end do
         end if
      end do
      return 
      end

**************storing start and end wet grid points for each layer of subdoamin*************
c if I1 & I2 are the coordinates of subdomain
c where '2' represents each subdoamin layer varing from I1-2 to I2+2
      subroutine local_izet(izet,izets,izete,icord)
      include 'C_model.f'
      include 'C_mpi.f'
      common /lbkhor/ khorl,lzet(khor),lb0(n),le0(n),
     * lb1(n),le1(n),lb2(n),le2(n)
      common /cord/ I1,I2,J1,J2,IA1,IA2,JA1,JA2,IB1,IB2,JB1,JB2
     1,ID1,ID2,JD1,JD2
      dimension icord(nprocs,4),izets(m,n),izete(m,n),izet(m,n)
      common /out11/ mzet(khor,nprocs),lkhor(nprocs)
      mzet=0
      lkhor=0
      lzet=0
      khorl=0
      lb0=0
      le0=0
      lb1=0
      le1=0
      lb2=0
      le2=0
      print*,'test icord',icord,IA2
      do nn=1,nprocs
      k=0
      DO J=ICORD(nn,3),ICORD(nn,4)
      DO I=ICORD(nn,1),ICORD(nn,2)
      if(izet(I,J) .ne. 0) THEN
      k=k+1
      mzet(k,nn)=izet(I,J)
      ENDIF
      ENDDO
      ENDDO
      lkhor(nn)=k
      enddo



      k=0
      DO J=J1,J2
      DO I=I1,I2
      if(izet(I,J) .ne. 0) THEN
      k=k+1
      lzet(k)=izet(I,J)
      ENDIF
      ENDDO
      ENDDO
      khorl=k

      DO J=JB1,JB2

c      IF(J .GE. J1 .AND. J .LE. J2) THEN
       lb0(J)=izets(I1,J)
       le0(J)=izete(I2,J)
c      ENDIF

c      IF(J .GE. JA1 .AND. J .LE. JA2) THEN
      lb1(J)=izets(IA1,J)
      
      le1(J)=izete(IA2,J)
c      ENDIF

      lb2(J)=izets(IB1,J)
      le2(J)=izete(IB2,J)

      ENDDO 
      return
      end
*************************Initialization of some variables for coomunication************************************
      subroutine INTI_COMM_VAR
      common /initc/ INUM1,INUM2,INUM3(6),INUM4(6)
     &,IRNUM1,IRNUM2,IRNUM3(6),IRNUM4(6),nreq1,nreq2,
     & nreq3,nreq4
       integer I
c     pointer for vertical communication  
      INUM1=0
      INUM2=0
      IRNUM1=0
      IRNUM2=0
c     pointer for horizontal communication
      DO I=1,6
      INUM3(I)=0
      INUM4(I)=0
      IRNUM3(I)=0
      IRNUM4(I)=0
      enddo
c     communication checking pointers for horizontal communication 
      nreq1=0
      nreq2=0
      nreq3=0
      nreq4=0
      return
      end
************************Comunication for Y direction layer************************************
      subroutine com_vertical(A,indr,B,iilo,NR,NL,ICON,NDIM,ISEND,lazc)
      include 'C_mpi.f'
      include 'C_model.f'
      include 'C_index.f'
      dimension A(indr),B(m,n,iilo),lazc(khor)
      common /verexc/ temp1(10*N*ilo),temp2(10*N*ilo),temp3(10*N*ilo),
     &temp4(10*N*ilo)
      common /cord/ I1,I2,J1,J2,IA1,IA2,JA1,JA2,IB1,IB2,JB1,JB2
     1,ID1,ID2,JD1,JD2
      common /initc/ INUM1,INUM2,INUM3(6),INUM4(6)
     &,IRNUM1,IRNUM2,IRNUM3(6),IRNUM4(6),nreq1,nreq2,
     & nreq3,nreq4 
      common /vcord/nvrlayer(4),nvllayer(4),nvrvar(4*n,3),nvlvar(4*n,3)
     *,nvrdep(4),nvldep(4)
      common /topo1/ izet(m,n),ltief(m,n)
      common /hreq/ir_r_b(50),ir_s_b(6),ir_r_t(6),ir_s_t(6)

      IF(mod(myid,nPh) .ne. 0 .and. nPh .gt. 1)THEN
      IF(NR .EQ. 0) goto 4
      IF(NDIM .EQ. 2) THEN
      IF(NR .EQ. 1) THEN
      nvar3=nvllayer(2)
      ELSE IF(NR .EQ. 2) THEN
      nvar3=nvllayer(1)+nvllayer(2)
      ENDIF
      ELSE
      IF(NR .EQ. 1) THEN
      nvar3=nvldep(2)
      ELSE IF(NR .EQ. 2) THEN
      nvar3=nvldep(1)+nvldep(2)
      ENDIF
      ENDIF
      IRNUM2=IRNUM2+nvar3
      IF(IRNUM2 .GT. 0 .AND. ISEND .EQ. 1) THEN
      nreq1=nreq1+1
      CALL MPI_IRECV(temp2,IRNUM2,MPI_REAL,myid-1,0
     3,MPI_COMM_WORLD,ir_r_b(nreq1),ierr)
      ENDIF

 4    continue
      ENDIF

      IF(mod((myid+1),nPh) .ne. 0 .and. nPh .gt. 1)THEN
      IF(NL .EQ. 0) goto 2
      IF(NDIM .EQ. 2) THEN
      IF(NL .EQ. 1) THEN
      nvar=nvrlayer(3)
      ELSE
      nvar=nvrlayer(3)+nvrlayer(4)
      ENDIF
      ELSE
      IF(NL .EQ. 1) THEN
      nvar=nvrdep(3)
      ELSE
      nvar=nvrdep(3)+nvrdep(4)
      ENDIF
      ENDIF
      IRNUM1=IRNUM1+nvar
      IF(IRNUM1 .GT. 0 .AND. ISEND .EQ. 1) THEN
      nreq1=nreq1+1
      CALL MPI_IRECV(temp4,IRNUM1,MPI_REAL,myid+1,1,MPI_COMM_WORLD
     3,ir_r_b(nreq1),ierr)
      ENDIF

 2    continue
      ENDIF

      IF(mod((myid+1),nPh) .ne. 0 .and. nPh .gt. 1)THEN
      IF(NR .EQ. 0) goto 1

      IF(NDIM .EQ. 2)THEN
      IF(NR .EQ. 1) THEN
      nvar1=nvrlayer(2)
      
      IF(ICON .EQ. 1) THEN
      npos=nvrlayer(1)    
      DO k=1,nvar1
      lw=nvrvar(npos+k,1)
      INUM1=INUM1+1
      temp1(INUM1)=A(lw)
      ENDDO
      ELSE         !IF(ICON .EQ. 1) 
      I=I2
      DO J=J1,J2
      IF(IZET(I,J) .GT. 0)THEN
      INUM1=INUM1+1
      temp1(INUM1)=B(I,J,1)
      ENDIF
      ENDDO
      ENDIF       !IF(ICON .EQ. 1) 
      ELSE        !IF(NR .EQ. 1)
      nvar1=nvrlayer(2)+nvrlayer(1)
      IF(ICON .EQ. 1) THEN
      DO k=1,nvar1
      lw=nvrvar(k,1)
      INUM1=INUM1+1
      temp1(INUM1)=A(lw)
      ENDDO
      ELSE         !IF(ICON .EQ. 1)
      DO I=I2-1,I2
      DO J=J1,J2
      IF(IZET(I,J) .GT. 0)THEN
      INUM1=INUM1+1
      temp1(INUM1)=B(I,J,1)
      ENDIF
      ENDDO
      ENDDO

      ENDIF         !IF(ICON .EQ. 1)
      ENDIF         !IF(NR .EQ. 1)
      ELSE

      IF(NR .EQ. 1) THEN
      nvar0=nvrlayer(2)
      npos=nvrlayer(1)
      IP1=I2
      IP2=I2
      ELSE
      nvar0=nvrlayer(1)+nvrlayer(2)
      npos=0
      IP1=I2-1
      IP2=I2
      ENDIF

      
      IF(ICON .GT. 0) THEN
      DO k=1,nvar0
      lw=nvrvar(npos+k,1)
      nwet=indwet(lw)
      DO J=1,lazc(lw)
      INUM1=INUM1+1
      temp1(INUM1)=A(nwet+J)
      ENDDO
      ENDDO
      ELSE
      DO 9 I=IP1,IP2
      DO 9 J=J1,J2
      ldown=ltief(I,J)
      DO 9 K=1,ldown
      INUM1=INUM1+1
      temp1(INUM1)=B(I,J,K)
 9    continue
      ENDIF

      ENDIF
     
      IF(INUM1 .GT. 0 .AND. ISEND .EQ. 1) THEN      
      nreq1=nreq1+1
      CALL MPI_ISEND(temp1,INUM1,MPI_REAL,myid+1,0,MPI_COMM_WORLD
     3,ir_r_b(nreq1), ierr)
      ENDIF

 1    continue      
      ENDIF 

      IF(mod(myid,nPh) .ne. 0 .and. nPh .gt. 1)THEN
      IF(NL .EQ. 0) goto 3
      IF(NDIM .EQ. 2)THEN
      IF(NL .EQ. 1) THEN
      nvar2=nvllayer(3)
      IF(ICON .EQ. 1) THEN
      npos=nvllayer(1)+ nvllayer(2)   
      
      DO k=1,nvar2
      lw=nvlvar(npos+k,1)
      INUM2=INUM2+1
      temp3(INUM2)=A(lw)
      ENDDO
      ELSE
      I=I1
      DO J=J1,J2
      IF(IZET(I,J) .GT. 0)THEN
      INUM2=INUM2+1
      temp3(INUM2)=B(I,J,1)
      ENDIF
      ENDDO
      ENDIF
      ELSE
      
      nvar2=nvllayer(3)+nvllayer(4)
      IF(ICON .EQ. 1) THEN
      npos=nvllayer(1)+ nvllayer(2)
      DO k=1,nvar2
      lw=nvlvar(npos+k,1)
      INUM2=INUM2+1
      temp3(INUM2)=A(lw)
      ENDDO
      ELSE
      DO I=I1,I1+1
      DO J=J1,J2
      IF(IZET(I,J) .GT. 0)THEN
      INUM2=INUM2+1
      temp3(INUM2)=B(I,J,1)
      ENDIF
      ENDDO
      ENDDO
      ENDIF
      ENDIF
      ELSE
      
      IF(NL .EQ. 1) THEN
      nvar0=nvllayer(3)
      npos=nvllayer(1)+nvllayer(2)
      IP1=I1
      IP2=I1
      ELSE
      nvar0=nvllayer(3)+nvllayer(4)
      npos=nvllayer(1)+nvllayer(2)
      IP1=I1
      IP2=I1+1
      ENDIF

      IF(ICON .GT. 0) THEN
      DO k=1,nvar0
      lw=nvlvar(npos+k,1)
      nwet=indwet(lw)
      DO J=1,lazc(lw)
      INUM2=INUM2+1
      temp3(INUM2)=A(nwet+J)
      ENDDO
      ENDDO
      ELSE
      DO 8 I=IP1,IP2
      DO 8 J=J1,J2
      ldown=ltief(I,J)
      DO 8 K=1,ldown
      INUM2=INUM2+1
      temp3(INUM2)=B(I,J,K)
 8    continue
      
      ENDIF
      ENDIF
      IF(INUM2 .GT. 0 .AND. ISEND .EQ. 1) THEN      
      nreq1=nreq1+1
      CALL MPI_ISEND(temp3,INUM2,MPI_REAL,myid-1,1,MPI_COMM_WORLD
     3,ir_r_b(nreq1),ierr)
      ENDIF

 3    continue      
      ENDIF
      return
      end

**********************Communicatind data between neighobrs in X-direction **************************************************
      subroutine com_horizontal(A,indr,B,iilo,NB,NT,ICON,NDIM,ISEND,
     *icdepb,idepb,icdept,idept,ncsendb,nsendb,ncsendt,nsendt,
     *nbhlayer,nthlayer)
      include 'C_mpi.f'
      include 'C_model.f'
      include 'C_index.f'
      PARAMETER(khor1=khor+1)
      dimension A(indr),B(M,N,iilo)
      dimension icdepb(nprocs),idepb(nprocs),icdept(nprocs)
     & ,idept(nprocs),nsendb(nprocs),ncsendb(nprocs),nsendt(nprocs),
     & ncsendt(nprocs),nbhlayer(nprocs*6,7),nthlayer(nprocs*6,7)
      common /surface/izets(m,n),izete(m,n)
      common /cord/ I1,I2,J1,J2,IA1,IA2,JA1,JA2,IB1,IB2,JB1,JB2
     1,ID1,ID2,JD1,JD2
      common /horexc/t_recv_b(m*ilo*nPh*2,6),t_send_b(m*ilo*nPh*2,6),
     2 t_recv_t(m*ilo*nPh*2,6),t_send_t(m*ilo*nPh*2,6)
      common /hreq/ir_r_b(50),ir_s_b(6),ir_r_t(6),ir_s_t(6)
      common /initc/ INUM1,INUM2,INUM3(6),INUM4(6)
     &,IRNUM1,IRNUM2,IRNUM3(6),IRNUM4(6),nreq1,nreq2,
     & nreq3,nreq4       
      common /topo1/izet(m,n),ltief(m,n)
      common /ind/ iwet(khor1),ldep(khor),lazc(khor),indend(n),
     *             isornr(n),isorsr(n),islab(n) 

      If(myid .GT. nPh-1 .and. nPv .gt. 1) THEN
      IF(NT .EQ. 0) goto 13
      npos=icdept(myid+1-nPh)-1
      DO 8 I=1,IDEPT(myid+1-nPh)
      id=nthlayer(npos+I,1)
      IF(NDIM .EQ. 2)THEN
      IX=nthlayer(npos+I,2)
      IY=nthlayer(npos+I,3)
      lw1=izets(IX,J1-1)
      lw2=izete(IY,J1-1)

      ndata4=lw2-lw1+1
      IF(ndata4 .LT. 0) ndata4=0
      IF(NT .EQ. 2) THEN
      lw3=izets(IX,J1-2)
      lw4=izete(IY,J1-2)
      ndata1=lw4-lw3+1
      IF(ndata1 .LT. 0) ndata1=0
      ndata4=ndata4+ndata1
      ENDIF
      ELSE
      ndata4=nthlayer(npos+I,5)
      IF(NT .EQ. 2) ndata4=ndata4+nthlayer(npos+I,6)
      ENDIF
      IRNUM4(I)=IRNUM4(I)+ndata4
      IF(IRNUM4(I) .GT. 0 .AND. ISEND .EQ. 1) THEN
      nreq1=nreq1+1
      call MPI_IRECV(t_recv_t(1,I),IRNUM4(I),MPI_REAL,id,0,
     5 MPI_COMM_WORLD, ir_r_b(nreq1), ierr)
      ENDIF
 8    continue
 13   continue
      endif
      If(myid .LT. nprocs-nPh .and. nPv .gt. 1) THEN           
      IF(NB .EQ. 0) goto 11
      npos=icdepb(myid+1)-1
      DO 6 I=1,idepb(myid+1)
      id=nbhlayer(npos+I,1)
      IF(NDIM .EQ. 2)THEN
      IX=nbhlayer(npos+I,2)
      IY=nbhlayer(npos+I,3)
      J=J2+1
      lw1=izets(IX,J)
      lw2=izete(IY,J)
      ndata=lw2-lw1+1 
      IF(ndata .LT. 0) ndata=0    
      IF(NB .EQ. 2) THEN
      lw3=izets(IX,J+1)
      lw4=izete(IY,J+1)
      ndata1=lw4-lw3+1
      IF(ndata1 .LT. 0) ndata1=0
      ndata=ndata+ndata1
      ENDIF
      ELSE
      ndata=nbhlayer(npos+I,5)
      IF(NB .EQ. 2) THEN
      ndata=ndata+nbhlayer(npos+I,6)
      ENDIF
      ENDIF
      IRNUM3(I)=IRNUM3(I)+ndata
      IF( IRNUM3(I).GT. 0 .AND. ISEND .EQ. 1) THEN
      nreq1=nreq1+1
      call MPI_IRECV(t_recv_b(1,I),IRNUM3(I),MPI_REAL,id,1,
     &MPI_COMM_WORLD, ir_r_b(nreq1), ierr)
      ENDIF
 6    continue
 11   continue
      endif
      If(myid .LT. nprocs-nPh .and. nPv .gt. 1) THEN
      IF(NT .EQ. 0) goto 12
      npos=ncsendt(myid+1)-1
      DO 7 I=1,nsendt(myid+1)
      kk=nthlayer(npos+I,4)
      id=nthlayer(kk,7)
      IX=nthlayer(kk,2)
      IY=nthlayer(kk,3)
      lw1=izets(IX,J2)
      lw2=izete(IY,J2)

      IF(NDIM .EQ. 2) THEN
      ndata3=lw2-lw1+1
      IF(ndata3 .LT. 0) ndata3=0

      IF(ICON .EQ. 1) THEN          
      DO K=1,ndata3
      INUM3(I)=INUM3(I)+1
      t_send_b(INUM3(I),I)=A(lw1+K-1)
      enddo
      ELSE
      DO K=1,ndata3
      II=iwet(lw1+K-1)
      INUM3(I)=INUM3(I)+1
      t_send_b(INUM3(I),I)=B(II,J2,1)
      ENDDO
      ENDIF
      
      
      IF(NT .EQ. 2) THEN
      lw3=izets(IX,J2-1)
      lw4=izete(IY,J2-1)
      ndata1=lw4-lw3+1
      IF(ndata1 .LT. 0) ndata1=0
      
      IF(ICON .EQ. 1) THEN
      DO K=1,ndata1
      INUM3(I)=INUM3(I)+1
      t_send_b(INUM3(I),I)=A(lw3+K-1)
      enddo
      ELSE
      DO K=1,ndata1
      II=iwet(lw3+K-1)
      INUM3(I)=INUM3(I)+1
      t_send_b(INUM3(I),I)=B(II,J2-1,1)
      ENDDO
      ENDIF
      ENDIF
      ELSE
      ndata3=nthlayer(kk,5)
      nwet1=indwet(lw1)
      IF(ICON .EQ. 1) THEN
      DO K=1,ndata3
      INUM3(I)=INUM3(I)+1
      t_send_b(INUM3(I),I)=A(nwet1+K)
      enddo
      ELSE
      n1=lw2-lw1+1
      IF(n1 .LT. 0) n1=0
      DO K1=1,n1
      II=iwet(lw1+K1-1)
      ldown=ltief(II,J2)
      DO K=1,ldown
      INUM3(I)=INUM3(I)+1
      t_send_b(INUM3(I),I)=B(II,J2,K)
      ENDDO
      ENDDO
      ENDIF
      IF(NT .EQ. 2) THEN
      lw3=izets(IX,J2-1)
      lw4=izete(IY,J2-1)
      n2=lw4-lw3+1
      IF(n2 .LT. 0) n2=0
      n1=nthlayer(kk,6)
      nwet1=indwet(lw3)
      IF(ICON .EQ. 1) THEN
      DO K=1,n1
      INUM3(I)=INUM3(I)+1
      t_send_b(INUM3(I),I)=A(nwet1+K)
      enddo
      ELSE
      DO K1=1,n1
      II=iwet(lw3+K1-1)
      ldown=ltief(II,J2-1)
      DO K=1,ldown
      INUM3(I)=INUM3(I)+1
      t_send_b(INUM3(I),I)=B(II,J2-1,K)
      ENDDO
      ENDDO
      ENDIF
      ENDIF
      ENDIF
      IF(INUM3(I) .GT. 0 .AND. ISEND .GT. 0) THEN
      nreq1=nreq1+1
      call MPI_ISEND(t_send_b(1,I),INUM3(I),MPI_REAL,id,0,MPI_COMM_WORLD
     &,ir_r_b(nreq1), ierr)
      ENDIF      
 7    continue
 12   continue
         
      ENDIF

      If(myid .GT. nPh-1 .and. nPv .gt. 1) THEN
      IF(NB .EQ. 0) goto 14
      npos=ncsendb(myid+1-nPh)-1
      J=J1
      DO I=1,nsendb(myid+1-nPh)
      kk=nbhlayer(npos+I,4)
      id=nbhlayer(kk,7)
      IX=nbhlayer(kk,2)
      IY=nbhlayer(kk,3)
      lw1=izets(IX,J1)
      lw2=izete(IY,J1)
      IF(NDIM .EQ. 2)THEN
      ndata2=lw2-lw1+1
      IF(ndata2 .LT. 0) ndata2=0

      IF(ICON .EQ. 1) THEN
      DO K=1,ndata2
      INUM4(I)=INUM4(I)+1
      t_send_t(INUM4(I),I)=A(lw1+K-1)
      enddo
      ELSE
      DO K=1,ndata2
      II=iwet(lw1+K-1)
      INUM4(I)=INUM4(I)+1
      t_send_t(INUM4(I),I)=B(II,J1,1)
      ENDDO
      ENDIF

      IF(NB .EQ. 2) THEN
      lw3=izets(IX,J1+1)
      lw4=izete(IY,J1+1)
      ndata1=lw4-lw3+1
      IF(ndata1 .LT. 0) ndata1=0
     
      IF(ICON .EQ. 1) THEN
      DO K=1,ndata1
      INUM4(I)=INUM4(I)+1
      t_send_t(INUM4(I),I)=A(lw3+K-1)
      enddo
      ELSE
      DO K=1,ndata1
      II=iwet(lw3+K-1)
      INUM4(I)=INUM4(I)+1
      t_send_t(INUM4(I),I)=B(II,J1+1,1)
      ENDDO      
      ENDIF
      ENDIF
      ELSE
      ndata2=nbhlayer(KK,5)
      nwet1=indwet(lw1)
      IF(ICON .EQ. 1) THEN
      DO K=1,ndata2
      INUM4(I)=INUM4(I)+1
c     if (INUM4(I).ge.10620)then
c      print*,INUM4(I),K,nwet1+K,KK
c     endif
      t_send_t(INUM4(I),I)=A(nwet1+K)
      enddo
      ELSE
      n1=lw2-lw1+1
      IF(n1 .LT. 0) n1=0
      DO K1=1,n1
      II=iwet(lw1+K1-1)
      ldown=ltief(II,J1)
      DO K=1,ldown
      INUM4(I)=INUM4(I)+1
      t_send_t(INUM4(I),I)=B(II,J1,K)
      ENDDO
      ENDDO
      ENDIF
      IF(NB .EQ. 2) THEN
      lw3=izets(IX,J1+1)
      lw4=izete(IY,J1+1)
      n2=lw4-lw3+1
      IF(n2 .LT. 0) n2=0
      n1=nbhlayer(KK,6)
      nwet1=indwet(lw3)

      IF(ICON .EQ. 1) THEN
      DO K=1,n1
      INUM4(I)=INUM4(I)+1
      t_send_t(INUM4(I),I)=A(nwet1+K)
      enddo
      ELSE
      DO K1=1,n2
      II=iwet(lw3+K1-1)
      ldown=ltief(II,J1+1)
      DO K=1,ldown
      INUM4(I)=INUM4(I)+1
      t_send_t(INUM4(I),I)=B(II,J1+1,K)
      ENDDO
      ENDDO
      ENDIF
      ENDIF
      ENDIF
      IF(INUM4(I) .GT. 0 .AND. ISEND .EQ. 1) THEN
      nreq1=nreq1+1
      call MPI_ISEND(t_send_t(1,I),INUM4(I),MPI_REAL,id,1,MPI_COMM_WORLD
     &,ir_r_b(nreq1), ierr)
      ENDIF
      ENDDO
      ENDIF
 14   continue           
      return
      end
**************************distributing data on neighbors in Y direction********************************
      subroutine VER_UNPACK(A,indr,B,iilo,NR,NL,ICON,NDIM,lazc)
      include 'C_model.f'
      include 'C_index.f'
      include 'C_mpi.f'
 
      dimension A(indr),B(m,n,iilo),lazc(khor)
      common /topo1/ izet(m,n),ltief(m,n)
      common /initc/ INUM1,INUM2,INUM3(6),INUM4(6)
     &,IRNUM1,IRNUM2,IRNUM3(6),IRNUM4(6),nreq1,nreq2,
     & nreq3,nreq4
       common /verexc/ temp1(10*N*ilo),temp2(10*N*ilo),temp3(10*N*ilo),
     &temp4(10*N*ilo)
      common /vcord/ nvrlayer(4),nvllayer(4),nvrvar(4*n,3),nvlvar(4*n,3)
     *,nvrdep(4),nvldep(4)
      common /cord/ I1,I2,J1,J2,IA1,IA2,JA1,JA2,IB1,IB2,JB1,JB2
     1,ID1,ID2,JD1,JD2
 
      IF(mod((myid+1),nPh) .ne. 0 .and. nPh .gt. 1) THEN
      IF(NL .EQ. 0) goto 1

      IF(NDIM .EQ. 2)THEN

      IF(NL .EQ. 1) THEN
      nvar=nvrlayer(3)
      II1=I2+1
      II2=I2+1
      ELSE
      nvar=nvrlayer(3)+nvrlayer(4)
      II1=I2+1
      II2=I2+2
      ENDIF

      IF(ICON .EQ. 1) THEN
      npos=nvrlayer(1)+nvrlayer(2)
      DO I=1,nvar
      lw=nvrvar(npos+I,1)
      INUM1=INUM1+1
      A(lw)=temp4(INUM1)
      ENDDO
      ELSE
c      WRITE(*,*)II1,II2,J1,J2,myid,NL
      DO 3 I=II1,II2
      DO 3 J=J1,J2
      IF(IZET(I,J) .GT. 0)THEN
      INUM1=INUM1+1
      B(I,J,1)=temp4(INUM1)
      ENDIF
 3    continue
      ENDIF
      ELSE
      IF(NL .EQ. 1) THEN
      nvar=nvrlayer(3)
      II1=I2+1
      II2=I2+1
      ELSE
      nvar=nvrlayer(3)+nvrlayer(4)
      II1=I2+1
      II2=I2+2
      ENDIF
      IF(ICON .EQ. 1) THEN
      npos=nvrlayer(1)+nvrlayer(2)
     
      DO 4 I=1,nvar
      lw=nvrvar(npos+I,1)
      nwet=indwet(lw)
      DO 4 J=1,lazc(lw)
      INUM1=INUM1+1
      A(nwet+J)=temp4(INUM1)
 4    continue
      ELSE
      DO 5 I=II1,II2
      DO 5 J=J1,J2
      ldown=ltief(I,J)
      DO 5 K=1,ldown
      INUM1=INUM1+1
      B(I,J,K)=temp4(INUM1)
 5    continue
      ENDIF
      ENDIF
 1    continue
      ENDIF

      IF(mod(myid,nPh) .ne. 0 .and. nPh .gt. 1) THEN
      IF(NR .EQ. 0) goto 2
      IF(NDIM .EQ. 2)THEN
      IF(NR .EQ. 1) THEN
      nvar=nvllayer(2)
      II1=I1-1
      II2=I1-1
      npos=nvllayer(1)
      ELSE
      nvar=nvllayer(1)+nvllayer(2)
      II1=I1-2
      II2=I1-1
      npos=0
      ENDIF
      IF(ICON .EQ. 1) THEN
      DO I=1,nvar
      lw=nvlvar(npos+I,1)
      INUM2=INUM2+1
      A(lw)=temp2(INUM2)
      ENDDO
      ELSE
      DO 6 I=II1,II2
      DO 6 J=J1,J2
      IF(IZET(I,J) .GT. 0)THEN
      INUM2=INUM2+1
      B(I,J,1)=temp2(INUM2)
      ENDIF
 6    continue
      ENDIF
      ELSE
      IF(NR .EQ. 1) THEN
      nvar=nvllayer(2)
      II1=I1-1
      II2=I1-1
      npos=nvllayer(1)
      ELSE
      nvar=nvllayer(1)+nvllayer(2)
      II1=I1-2
      II2=I1-1
      npos=0
      ENDIF
      IF(ICON .EQ. 1) THEN
      DO 7 I=1,nvar
      lw=nvlvar(npos+I,1)
      nwet=indwet(lw)
      DO 7 J=1,lazc(lw)
      INUM2=INUM2+1
      A(nwet+J)=temp2(INUM2)
 7    continue
      ELSE
      DO 8 I=II1,II2
      DO 8 J=J1,J2
      ldown=ltief(I,J)
      DO 8 K=1,ldown
      INUM2=INUM2+1
      B(I,J,K)=temp2(INUM2)
 8    continue
      ENDIF
      ENDIF
 2    continue
      ENDIF 
      return
      end
**************************Distributing data in X-direction on neighbors****************************
      subroutine HOR_UNPACK(A,indr,B,iilo,NB,NT,nbhlayer,nthlayer,
     &idepb,idept,icdepb,icdept,ICON,NDIM,iwet)

      include 'C_model.f'
      include 'C_mpi.f'
      include 'C_index.f'
      parameter(khor1=khor+1)
      dimension A(indr),B(m,n,iilo),iwet(khor1)
      dimension icdepb(nprocs),idepb(nprocs),icdept(nprocs)
     &,idept(nprocs)
      dimension nbhlayer(nprocs*6,7),nthlayer(nprocs*6,7)
      common /horexc/ t_recv_b(m*ilo*nPh*2,6),t_send_b(m*ilo*nPh*2,6),
     2 t_recv_t(m*ilo*nPh*2,6),t_send_t(m*ilo*nPh*2,6)
      common /initc/ INUM1,INUM2,INUM3(6),INUM4(6)
     &,IRNUM1,IRNUM2,IRNUM3(6),IRNUM4(6),nreq1,nreq2,
     & nreq3,nreq4
      common /surface/ izets(m,n),izete(m,n)
      common /topo1/ izet(m,n),ltief(m,n)
      common /cord/ I1,I2,J1,J2,IA1,IA2,JA1,JA2,IB1,IB2,JB1,JB2
     1,ID1,ID2,JD1,JD2
      If(myid .LT. nprocs-nPh .and. nPv .gt. 1) THEN           
      IF(NB .EQ. 0) goto 1
      npos=icdepb(myid+1)-1
      DO I=1,IDEPB(myid+1)
      IX=nbhlayer(npos+I,2)
      IY=nbhlayer(npos+I,3)
      J=J2+1
      lw1=izets(IX,J)
      lw2=izete(IY,J)
      n1=lw2-lw1+1
      it=INUM3(I)
      if(n1 .LT. 0) n1=0
      IF(NB. EQ. 2) THEN
      lw3=izets(IX,J+1)
      lw4=izete(IY,J+1)
      n2=lw4-lw3+1
      if(n2 .LT. 0) n2=0
      ENDIF
      IF(NDIM .EQ. 2) THEN
      INUM3(I)=INUM3(I)+n1
      IF(ICON .EQ. 1) THEN
      
      DO k=1,n1 
      A(lw1+k-1)=t_recv_b(it+k,I)
      enddo
      
      ELSE
      DO k=1,n1
      II=iwet(lw1+k-1)
      B(II,J,1)=t_recv_b(it+k,I)
      ENDDO
      
      ENDIF
      IF(NB .EQ. 2) THEN
      
      it=INUM3(I)
      INUM3(I)=INUM3(I)+n2
      IF(ICON .EQ. 1) THEN
      
      DO k=1,n2 
      A(lw3+k-1)=t_recv_b(it+k,I)
      enddo
      
      ELSE
      
      DO k=1,n2
      II=iwet(lw3+k-1)
      B(II,J+1,1)=t_recv_b(it+k,I)
      ENDDO
      ENDIF
      ENDIF
      ELSE IF(NDIM .EQ. 3) THEN
      num=0
      nwet1=indwet(lw1)
      n3=nbhlayer(npos+I,5)
      INUM3(I)=INUM3(I)+n3
      IF(ICON .EQ. 1) THEN
      DO k=1,n3
      A(nwet1+k)=t_recv_b(it+k,I)
      enddo
      ELSE
c      write(*,*)n3,myid,INUM3(I),I,'n3'
      DO k1=1,n1
      II=iwet(lw1+k1-1)
      ldown=ltief(II,J)
      DO K=1,ldown
      num=num+1
      B(II,J,K)=t_recv_b(it+num,I)
      ENDDO
      ENDDO
      ENDIF

      IF(NB .EQ. 2) THEN
      it=INUM3(I)
      n4=nbhlayer(npos+I,6)
      INUM3(I)=INUM3(I)+n4
      nwet2=indwet(lw3)
      IF(ICON .EQ. 1) THEN
      DO k=1,n4
      A(nwet2+k)=t_recv_b(it+k,I)
      enddo
      ELSE
      num=0
      DO k1=1,n2
      II=iwet(lw3+k1-1)
      ldown=ltief(II,J+1)
      DO K=1,ldown
      num=num+1
      B(II,J+1,K)=t_recv_b(it+num,I)
      ENDDO
      ENDDO
      ENDIF
      ENDIF
      ELSE
      WRITE(*,*) 'WRONG DIMENSION'
      STOP
      ENDIF
      ENDDO
      ENDIF
 1    continue

      If(myid .GT. nPh-1 .and. nPv .gt. 1) THEN           
      IF(NT .EQ. 0) goto 2 
      npos=icdept(myid+1-nPh)-1
c      if(myid .EQ. 7) WRITE(*,*)npos,IDEPT(myid+1-nPh)
      DO I=1,IDEPT(myid+1-nPh)
      IX=nthlayer(npos+I,2)
      IY=nthlayer(npos+I,3)
      J=J1-1
      lw1=izets(IX,J)
      lw2=izete(IY,J)
      n1=lw2-lw1+1
      it=INUM4(I)
      if(n1 .LT. 0) n1=0
      IF(NT. EQ. 2) THEN
      lw3=izets(IX,J-1)
      lw4=izete(IY,J-1)
      n2=lw4-lw3+1
      if(n2 .LT. 0) n2=0
      ENDIF
      IF(NDIM .EQ. 2) THEN
      INUM4(I)=INUM4(I)+n1
      IF(ICON .EQ. 1) THEN
      
      DO k=1,n1 
      A(lw1+k-1)=t_recv_t(it+k,I)
      enddo
      
      ELSE
      DO k=1,n1
      II=iwet(lw1+k-1)
      B(II,J,1)=t_recv_t(it+k,I)
c      if(myid .EQ. 7) WRITE(*,*)lw1+k-1,t_recv_t(it+k,I),Ii,J,it+k,'HH'
      ENDDO
      
      ENDIF
      IF(NT .EQ. 2) THEN
      
      it=INUM4(I)
      INUM4(I)=INUM4(I)+n2
      IF(ICON .EQ. 1) THEN
      
      DO k=1,n2 
      A(lw3+k-1)=t_recv_t(it+k,I)
      enddo
      
      ELSE
      
      DO k=1,n2
      II=iwet(lw3+k-1)
      B(II,J-1,1)=t_recv_t(it+k,I)
      ENDDO
      
      
      ENDIF
      ENDIF
      ELSE IF(NDIM .EQ. 3) THEN
      num=0
      nwet1=indwet(lw1)
      n3=nthlayer(npos+I,5)
      INUM4(I)=INUM4(I)+n3
      IF(ICON .EQ. 1) THEN
      DO k=1,n3
      A(nwet1+k)=t_recv_t(it+k,I)
      enddo
      ELSE
c      write(*,*)n3,myid,INUM4(I),I,'n3'

      DO k1=1,n1
      II=iwet(lw1+k1-1)
      ldown=ltief(II,J)
      DO K=1,ldown
      num=num+1
      B(II,J,K)=t_recv_t(it+num,I)
      ENDDO
      ENDDO
      ENDIF

      IF(NT .EQ. 2) THEN
      it=INUM4(I)
      n4=nthlayer(npos+I,6)
      INUM4(I)=INUM4(I)+n4
      nwet2=indwet(lw3)
      IF(ICON .EQ. 1) THEN
      DO k=1,n4
      A(nwet2+k)=t_recv_t(it+k,I)
      enddo
      ELSE
      num=0
      DO k1=1,n2
      II=iwet(lw3+k1-1)
      ldown=ltief(II,J-1)
      DO K=1,ldown
      num=num+1
      B(II,J-1,K)=t_recv_t(it+num,I)
      ENDDO
      ENDDO
      ENDIF
      ENDIF
      ELSE
      WRITE(*,*) 'WRONG DIMENSION'
      STOP
      ENDIF
      ENDDO
      ENDIF
 2    continue      
      return
      end

**************subdoamin: starting and end points for calculation of surface level************************
      subroutine  NEW_ISORNR(isornr,isorsr,iwet)
      include 'C_model.f'
      include 'C_mpi.f'
      parameter(khor1=khor+1)
      dimension isornr(n),isorsr(n),iwet(khor1)
      common /isor/ isornr1(n),isorsr1(n),isornr2(n),isorsr2(n)
     1,ik1(khor),norsor1
      common /cord/ I1,I2,J1,J2,IA1,IA2,JA1,JA2,IB1,IB2,JB1,JB2
     1,ID1,ID2,JD1,JD2

      isornr1=0      
      isorsr1=-1
      norsor1=0
      isornr2=0
      isorsr2=-1
      ik1=0
      ik = 0
      nz=n-1
      it1=0

      do k=1,n
        lwa = isornr(k)
        lwe = isorsr(k)
        it=0
        it1=0
        do lw=lwa,lwe
      i = iwet(lw)
      if(k .ge. 3 .and. k .le. nz) ik = ik+1
      if(k .ge. J1 .and. k .le. J2 .AND. i .ge. I1 .and. i .le. I2)then 
          if(k .ge. 3 .and. k .le. nz)then 
          norsor1=norsor1+1
          ik1(norsor1)=ik
          endif
          it=it+1
          if(it .EQ. 1)isornr1(k)=lw
          isorsr1(k)=lw
       endif 
      if(k .ge. JA1 .and. k .le. JA2 .and. i .ge. 
     3 IA1 .and. i .le. IA2)then 
          it1=it1+1
          if(it1 .EQ. 1)isornr2(k)=lw
          isorsr2(k)=lw
      endif
        enddo
      enddo
      return
      end
C************* same like as uvrand except one matrix at a time*******************************
      subroutine uvrand_par (A,izet,iindex,ITESTX,ITESTY)
      include 'C_model.f'
      include 'C_mpi.f'
      parameter(khor1=khor+1)
      dimension A(ndrei)
      dimension izet(m,n),iindex(m,n)
      integer ITESTX,ITESTY
      common /ind/ iwet(khor1),ldep(khor),lazc(khor), 
     * indend(n),isornr(n),isorsr(n),islab(n)
      common /bcor/ JWB1(2),JWB2(2),IWB1(2),IWB2(2),NWB
     &,JNB1(2),JNB2(2),INB1(2),INB2(2),NNB 

c      if(myid .EQ. 0) WRITE(*,*) INB1(1),INB2(1),JNB1(1),JNB2(1)
c    &,IWB1(1),IWB2(1),JWB1(1),JWB2(1) 
        do j=JNB1(1),JNB2(1)
        l   = iindex(INB1(1),j)  !First line
        ll  = iindex(INB1(1)+1,j)  !Second line
        lll = iindex(INB2(1),j)  !Third line
        IF(ITESTX .EQ. 1) THEN
        lll=ll      
        ENDIF
        lw  = izet(INB1(1),j)
        lump=0
        if(lw.gt.0)lump = lazc(lw)
        do k=1,lump
          A (l +k)  = A (lll+k)
          A (ll+k)  = A (lll+k)
         enddo
      enddo
      
     

      do i=IWB1(1),IWB2(1)
        l   = iindex(i,JWB1(1)) !First line
        ll  = iindex(i,JWB1(1)+1) !Second line
        lll = iindex(i,JWB2(1)) !Third line
        IF(ITESTY .EQ. 1) THEN
        lll=ll      
        ENDIF
        lw  = izet(i,JWB1(1))
             lump=0
        if(lw.gt.0)lump = lazc(lw)
        do k=1,lump
          A (l +k)  = A (lll+k)
          A (ll+k)  = A (lll+k)
        enddo
      enddo
      return
      end
***********************same like as com_vertical except Integer datas are commnicated*************************
      subroutine com_ivertical(IA,IB,NR,NL,ICON,ISEND,lazc)
      include 'C_mpi.f'
      include 'C_model.f'
      include 'C_index.f'
      dimension IA(1),IB(m,n,1),lazc(khor)
      common /iverexc/ itemp1(2*N),itemp2(2*N),itemp3(2*N),
     &itemp4(2*N)
      common /cord/ I1,I2,J1,J2,IA1,IA2,JA1,JA2,IB1,IB2,JB1,JB2
     1,ID1,ID2,JD1,JD2
      common /initc/ INUM1,INUM2,INUM3(6),INUM4(6)
     &,IRNUM1,IRNUM2,IRNUM3(6),IRNUM4(6),nreq1,nreq2,
     & nreq3,nreq4 
      common /vcord/nvrlayer(4),nvllayer(4),nvrvar(4*n,3),nvlvar(4*n,3)
     *,nvrdep(4),nvldep(4)
      common /topo1/ izet(m,n),ltief(m,n)
      common /hreq/ ir_r_b(50),ir_s_b(6),ir_r_t(6),ir_s_t(6)

      IF(mod(myid,nPh) .ne. 0)THEN
      IF(NR .EQ. 0) goto 4
      IF(NR .EQ. 1) THEN
      nvar3=nvllayer(2)
      ELSE IF(NR .EQ. 2) THEN
      nvar3=nvllayer(1)+nvllayer(2)
      ENDIF
      IRNUM2=IRNUM2+nvar3
      IF(IRNUM2 .GT. 0 .AND. ISEND .EQ. 1) THEN
      nreq1=nreq1+1
      CALL MPI_IRECV(itemp2,IRNUM2,MPI_INTEGER,myid-1,0,
     3 MPI_COMM_WORLD,ir_r_b(nreq1), ierr)
      ENDIF

 4    continue
      ENDIF

      IF(mod((myid+1),nPh) .ne. 0)THEN

      IF(NL .EQ. 0) goto 2
      IF(NL .EQ. 1) THEN
      nvar=nvrlayer(3)
      ELSE
      nvar=nvrlayer(3)+nvrlayer(4)
      ENDIF
      IRNUM1=IRNUM1+nvar
      IF(IRNUM1 .GT. 0 .AND. ISEND .EQ. 1) THEN
      nreq1=nreq1+1
      CALL MPI_IRECV(itemp4,IRNUM1,MPI_INTEGER,myid+1,1,MPI_COMM_WORLD,
     3ir_r_b(nreq1),ierr)
      ENDIF
 2    continue
      ENDIF

      IF(mod((myid+1),nPh) .ne. 0)THEN
      IF(NR .EQ. 0) goto 1
      IF(NR .EQ. 1) THEN
      nvar1=nvrlayer(2)
      IF(ICON .EQ. 1) THEN
      npos=nvrlayer(1)    
      DO k=1,nvar1
      lw=nvrvar(npos+k,1)
      INUM1=INUM1+1
      itemp1(INUM1)=IA(lw)
      ENDDO
      ELSE         !IF(ICON .EQ. 1) 
      I=I2
      DO J=J1,J2
      IF(IZET(I,J) .GT. 0)THEN
      INUM1=INUM1+1
      itemp1(INUM1)=IB(I,J,1)
      ENDIF
      ENDDO
      ENDIF       !IF(ICON .EQ. 1) 
      ELSE        !IF(NR .EQ. 1)
      nvar1=nvrlayer(2)+nvrlayer(1)
      IF(ICON .EQ. 1) THEN
      DO k=1,nvar1
      lw=nvrvar(k,1)
      INUM1=INUM1+1
      itemp1(INUM1)=IA(lw)
      ENDDO
      ELSE         !IF(ICON .EQ. 1)
      DO I=I2-1,I2
      DO J=J1,J2
      IF(IZET(I,J) .GT. 0)THEN
      INUM1=INUM1+1
      itemp1(INUM1)=IB(I,J,1)
      ENDIF
      ENDDO
      ENDDO
      ENDIF         !IF(ICON .EQ. 1)
      ENDIF         !IF(NR .EQ. 1)
      IF(INUM1 .GT. 0 .AND. ISEND .EQ. 1) THEN      
      nreq1=nreq1+1
      CALL MPI_ISEND(itemp1,INUM1,MPI_INTEGER,myid+1,0,MPI_COMM_WORLD,
     3ir_r_b(nreq1), ierr)
      ENDIF
 1    continue  
      ENDIF 
      

      IF(mod(myid,nPh) .ne. 0)THEN
      IF(NL .EQ. 0) goto 3

      IF(NL .EQ. 1) THEN
      nvar2=nvllayer(3)
      IF(ICON .EQ. 1) THEN
      npos=nvllayer(1)+ nvllayer(2)   
      DO k=1,nvar2
      lw=nvlvar(npos+k,1)
      INUM2=INUM2+1
      itemp3(INUM2)=IA(lw)
      ENDDO
      ELSE
      I=I1
      DO J=J1,J2
      IF(IZET(I,J) .GT. 0)THEN
      INUM2=INUM2+1
      itemp3(INUM2)=IB(I,J,1)
      ENDIF
      ENDDO
      ENDIF
      ELSE
      nvar2=nvllayer(3)+nvllayer(4)
      IF(ICON .EQ. 1) THEN
      npos=nvllayer(1)+ nvllayer(2)
      DO k=1,nvar2
      lw=nvlvar(npos+k,1)
      INUM2=INUM2+1
      itemp3(INUM2)=IA(lw)
      ENDDO
      ELSE
      DO I=I1,I1+1
      DO J=J1,J2
      IF(IZET(I,J) .GT. 0)THEN
      INUM2=INUM2+1
      itemp3(INUM2)=IB(I,J,1)
      ENDIF
      ENDDO
      ENDDO
      ENDIF
      ENDIF
      
      IF(INUM2 .GT. 0 .AND. ISEND .EQ. 1) THEN      
      nreq1=nreq1+1
      CALL MPI_ISEND(itemp3,INUM2,MPI_INTEGER,myid-1,1,MPI_COMM_WORLD,
     3 ir_r_b(nreq1),ierr)
      ENDIF

 3    continue      
      ENDIF
      return
      end

*************************************************************************************************************************************************************
      subroutine com_ihorizontal(IA,IB,NB,NT,ICON,ISEND,icdepb,idepb,
     *icdept,idept,ncsendb,nsendb,ncsendt,nsendt,nbhlayer,nthlayer)
      include 'C_mpi.f'
      include 'C_model.f'
      include 'C_index.f'
      parameter(khor1=khor+1) 
      dimension IA(1),IB(M,N,1)
      dimension icdepb(nprocs),idepb(nprocs),icdept(nprocs)
     & ,idept(nprocs)
      dimension nbhlayer(nprocs*6,7),nthlayer(nprocs*6,7)
      dimension nsendb(nprocs),ncsendb(nprocs),nsendt(nprocs),
     & ncsendt(nprocs)
      common /hreq/ ir_r_b(50),ir_s_b(6),ir_r_t(6),ir_s_t(6)
      common /ihorexc/ it_recv_b(m*nPh,6),it_send_b(m*nPh,6),
     2 it_recv_t(m*nPh,6),it_send_t(m*nPh,6)
      common /initc/ INUM1,INUM2,INUM3(6),INUM4(6)
     &,IRNUM1,IRNUM2,IRNUM3(6),IRNUM4(6),nreq1,nreq2,
     & nreq3,nreq4
      common /surface/ izets(m,n),izete(m,n)
      common /topo1/ izet(m,n),ltief(m,n)
      common /cord/ I1,I2,J1,J2,IA1,IA2,JA1,JA2,IB1,IB2,JB1,JB2
     1,ID1,ID2,JD1,JD2
      common /ind/ iwet(khor1),ldep(khor),lazc(khor),indend(n),
     *             isornr(n),isorsr(n),islab(n)

      If(myid .GT. nPh-1) THEN
      IF(NT .EQ. 0) goto 13
      npos=icdept(myid+1-nPh)-1
      DO 8 I=1,IDEPT(myid+1-nPh)
      id=nthlayer(npos+I,1)
      IX=nthlayer(npos+I,2)
      IY=nthlayer(npos+I,3)
      lw1=izets(IX,J1-1)
      lw2=izete(IY,J1-1)
      ndata4=lw2-lw1+1
      IF(ndata4 .LT. 0) ndata4=0
      IF(NT .EQ. 2) THEN
      lw3=izets(IX,J1-2)
      lw4=izete(IY,J1-2)
      ndata1=lw4-lw3+1
      IF(ndata1 .LT. 0) ndata1=0
      ndata4=ndata4+ndata1
      ENDIF
      IRNUM4(I)=IRNUM4(I)+ndata4
      IF(IRNUM4(I) .GT. 0 .AND. ISEND .EQ. 1) THEN
      nreq1=nreq1+1
      call MPI_IRECV(it_recv_t(1,I),IRNUM4(I),MPI_INTEGER,id,0,
     5 MPI_COMM_WORLD, ir_r_b(nreq1), ierr)
      ENDIF
 8    continue
 13   continue
      ENDIF

      If(myid .LT. nprocs-nPh) THEN           
      IF(NB .EQ. 0) goto 11
      npos=icdepb(myid+1)-1
      DO 6 I=1,IDEPB(myid+1)
      id=nbhlayer(npos+I,1)
      IX=nbhlayer(npos+I,2)
      IY=nbhlayer(npos+I,3)
      J=J2+1
      lw1=izets(IX,J)
      lw2=izete(IY,J)
      ndata=lw2-lw1+1 
      IF(ndata .LT. 0) ndata=0    
      IF(NB .EQ. 2) THEN
      lw3=izets(IX,J+1)
      lw4=izete(IY,J+1)
      ndata1=lw4-lw3+1
      IF(ndata1 .LT. 0) ndata1=0
      ndata=ndata+ndata1
      ENDIF
      
      IRNUM3(I)=IRNUM3(I)+ndata
      IF( IRNUM3(I).GT. 0 .AND. ISEND .EQ. 1) THEN
      nreq1=nreq1+1
      call MPI_IRECV(it_recv_b(1,I),IRNUM3(I),MPI_INTEGER,id,1,
     &MPI_COMM_WORLD, ir_r_b(nreq1), ierr)
      ENDIF
 6    continue
 11   continue
      ENDIF
      
      If(myid .LT. nprocs-nPh) THEN
      IF(NT .EQ. 0) goto 12
      npos=ncsendt(myid+1)-1
      DO 7 I=1,nsendt(myid+1)
      kk=nthlayer(npos+I,4)
      id=nthlayer(kk,7)
      IX=nthlayer(kk,2)
      IY=nthlayer(kk,3)
      lw1=izets(IX,J2)
      lw2=izete(IY,J2)

      
      ndata3=lw2-lw1+1
      IF(ndata3 .LT. 0) ndata3=0

      IF(ICON .EQ. 1) THEN          
      DO K=1,ndata3
      INUM3(I)=INUM3(I)+1
      it_send_b(INUM3(I),I)=IA(lw1+K-1)
      enddo
      ELSE
      DO K=1,ndata3
      II=iwet(lw1+K-1)
      INUM3(I)=INUM3(I)+1
      it_send_b(INUM3(I),I)=IB(II,J2,1)
      ENDDO
      ENDIF
      
      
      IF(NT .EQ. 2) THEN
      lw3=izets(IX,J2-1)
      lw4=izete(IY,J2-1)
      ndata1=lw4-lw3+1
      IF(ndata1 .LT. 0) ndata1=0
      
      IF(ICON .EQ. 1) THEN
      DO K=1,ndata1
      INUM3(I)=INUM3(I)+1
      it_send_b(INUM3(I),I)=IA(lw3+K-1)
      enddo
      ELSE
      DO K=1,ndata1
      II=iwet(lw3+K-1)
      INUM3(I)=INUM3(I)+1
      it_send_b(INUM3(I),I)=IB(II,J2-1,1)
      ENDDO
      ENDIF
      ENDIF
      
      IF(INUM3(I) .GT. 0 .AND. ISEND .GT. 0) THEN
      nreq1=nreq1+1
      call MPI_ISEND(it_send_b(1,I),INUM3(I),MPI_INTEGER,id,0,
     &MPI_COMM_WORLD,ir_r_b(nreq1), ierr)
      ENDIF      
 7    continue
 12   continue
         
      ENDIF
      
      If(myid .GT. nPh-1) THEN
      IF(NB .EQ. 0) goto 14
      npos=ncsendb(myid+1-nPh)-1
      J=J1
      DO I=1,nsendb(myid+1-nPh)
      kk=nbhlayer(npos+I,4)
      id=nbhlayer(KK,7)
      IX=nbhlayer(KK,2)
      IY=nbhlayer(KK,3)
      lw1=izets(IX,J1)
      lw2=izete(IY,J1)
      
      ndata2=lw2-lw1+1
      IF(ndata2 .LT. 0) ndata2=0
      IF(ICON .EQ. 1) THEN
      DO K=1,ndata2
      INUM4(I)=INUM4(I)+1
      it_send_t(INUM4(I),I)=IA(lw1+K-1)
      enddo
      ELSE
      DO K=1,ndata2
      II=iwet(lw1+K-1)
      INUM4(I)=INUM4(I)+1
      it_send_t(INUM4(I),I)=IB(II,J1,1)
      ENDDO
      ENDIF
      IF(NB .EQ. 2) THEN
      lw3=izets(IX,J1+1)
      lw4=izete(IY,J1+1)
      ndata1=lw4-lw3+1
      IF(ndata1 .LT. 0) ndata1=0
      IF(ICON .EQ. 1) THEN
      DO K=1,ndata1
      INUM4(I)=INUM4(I)+1
      it_send_t(INUM4(I),I)=IA(lw3+K-1)
      enddo
      ELSE
      DO K=1,ndata1
      II=iwet(lw3+K-1)
      INUM4(I)=INUM4(I)+1
      it_send_t(INUM4(I),I)=IB(II,J1+1,1)
      ENDDO      
      ENDIF
      ENDIF
      IF(INUM4(I) .GT. 0 .AND. ISEND .EQ. 1) THEN
      nreq1=nreq1+1
      call MPI_ISEND(it_send_t(1,I),INUM4(I),MPI_INTEGER,id,1,
     &MPI_COMM_WORLD,ir_r_b(nreq1), ierr)
      ENDIF
      ENDDO
      ENDIF
 14   continue           
      return
      end
************************************************************************************************************************************************************
      subroutine VER_IUNPACK(IA,IB,NR,NL,ICON,lazc)
      include 'C_model.f'
      include 'C_index.f'
      include 'C_mpi.f'
 
      dimension IA(1),IB(m,n,1),lazc(khor)
      common /topo1/ izet(m,n),ltief(m,n)
      common /initc/ INUM1,INUM2,INUM3(6),INUM4(6)
     &,IRNUM1,IRNUM2,IRNUM3(6),IRNUM4(6),nreq1,nreq2,
     & nreq3,nreq4
       common /iverexc/ itemp1(2*N),itemp2(2*N),itemp3(2*N),
     &itemp4(2*N)
      common /vcord/ nvrlayer(4),nvllayer(4),nvrvar(4*n,3),nvlvar(4*n,3)
     *,nvrdep(4),nvldep(4)
      common /cord/ I1,I2,J1,J2,IA1,IA2,JA1,JA2,IB1,IB2,JB1,JB2
     1,ID1,ID2,JD1,JD2
 
      IF(mod((myid+1),nPh) .ne. 0) THEN
      IF(NL .EQ. 0) goto 1
      IF(NL .EQ. 1) THEN
      nvar=nvrlayer(3)
      II1=I2+1
      II2=I2+1
      ELSE
      nvar=nvrlayer(3)+nvrlayer(4)
      II1=I2+1
      II2=I2+2
      ENDIF
      IF(ICON .EQ. 1) THEN
      npos=nvrlayer(1)+nvrlayer(2)
      DO I=1,nvar
      lw=nvrvar(npos+I,1)
      INUM1=INUM1+1
      IA(lw)=itemp4(INUM1)
      ENDDO
      ELSE
      DO 3 I=II1,II2
      DO 3 J=J1,J2
      IF(IZET(I,J) .GT. 0)THEN
      INUM1=INUM1+1
      IB(I,J,1)=itemp4(INUM1)
      ENDIF
 3    continue
      ENDIF
 1    continue
      ENDIF

      IF(mod(myid,nPh) .ne. 0) THEN
      IF(NR .EQ. 0) goto 2
      IF(NR .EQ. 1) THEN
      nvar=nvllayer(2)
      II1=I1-1
      II2=I1-1
      npos=nvllayer(1)
      ELSE
      nvar=nvllayer(1)+nvllayer(2)
      II1=I1-2
      II2=I1-1
      npos=0
      ENDIF
      IF(ICON .EQ. 1) THEN
      DO I=1,nvar
      lw=nvlvar(npos+I,1)
      INUM2=INUM2+1
      IA(lw)=itemp2(INUM2)
      ENDDO
      ELSE
      DO 6 I=II1,II2
      DO 6 J=J1,J2
      IF(IZET(I,J) .GT. 0)THEN
      INUM2=INUM2+1
      IB(I,J,1)=itemp2(INUM2)
      ENDIF
 6    continue
      ENDIF
 2    continue
      ENDIF 
      return
      end
************************************************************************************************************************************************************
      subroutine HOR_IUNPACK(IA,IB,NB,NT,nbhlayer,nthlayer,idepb,idept,
     &icdepb,icdept,ICON,iwet)
      include 'C_model.f'
      include 'C_mpi.f'
      include 'C_index.f'
      parameter(khor1=khor+1)
      dimension IA(1),IB(m,n,1),iwet(khor1)
      dimension icdepb(nprocs),idepb(nprocs),icdept(nprocs)
     &,idept(nprocs)
      dimension nbhlayer(nprocs*6,7),nthlayer(nprocs*6,7)
      common /ihorexc/ it_recv_b(m*nPh,6),it_send_b(m*nPh,6),
     2 it_recv_t(m*nPh,6),it_send_t(m*nPh,6)
      common /initc/ INUM1,INUM2,INUM3(6),INUM4(6)
     &,IRNUM1,IRNUM2,IRNUM3(6),IRNUM4(6),nreq1,nreq2,
     & nreq3,nreq4
      common /surface/ izets(m,n),izete(m,n)
      common /topo1/ izet(m,n),ltief(m,n)
      common /cord/ I1,I2,J1,J2,IA1,IA2,JA1,JA2,IB1,IB2,JB1,JB2
     1,ID1,ID2,JD1,JD2

      If(myid .LT. nprocs-nPh) THEN           
      IF(NB .EQ. 0) goto 1
      npos=icdepb(myid+1)-1
      DO I=1,IDEPB(myid+1)
      IX=nbhlayer(npos+I,2)
      IY=nbhlayer(npos+I,3)
      J=J2+1
      lw1=izets(IX,J)
      lw2=izete(IY,J)
      n1=lw2-lw1+1
      it=INUM3(I)
      if(n1 .LT. 0) n1=0
      IF(NB. EQ. 2) THEN
      lw3=izets(IX,J+1)
      lw4=izete(IY,J+1)
      n2=lw4-lw3+1
      if(n2 .LT. 0) n2=0
      ENDIF
      
      INUM3(I)=INUM3(I)+n1
      IF(ICON .EQ. 1) THEN
      DO k=1,n1 
      IA(lw1+k-1)=it_recv_b(it+k,I)
      enddo
      ELSE
      DO k=1,n1
      II=iwet(lw1+k-1)
      IB(II,J,1)=it_recv_b(it+k,I)
      ENDDO
      ENDIF
      IF(NB .EQ. 2) THEN
      it=INUM3(I)
      INUM3(I)=INUM3(I)+n2
      IF(ICON .EQ. 1) THEN
      DO k=1,n2 
      IA(lw3+k-1)=it_recv_b(it+k,I)
      enddo
      ELSE
      DO k=1,n2
      II=iwet(lw3+k-1)
      IB(II,J+1,1)=it_recv_b(it+k,I)
      ENDDO
      ENDIF
      ENDIF
      ENDDO
      ENDIF 
 1    continue

      If(myid .GT. nPh-1) THEN           
      IF(NT .EQ. 0) goto 1
      npos=icdept(myid+1-nPh)-1
      DO I=1,IDEPT(myid+1-nPh)
      IX=nthlayer(npos+I,2)
      IY=nthlayer(npos+I,3)
      J=J1-1
      lw1=izets(IX,J)
      lw2=izete(IY,J)
      n1=lw2-lw1+1
      it=INUM4(I)
      if(n1 .LT. 0) n1=0
      IF(NT. EQ. 2) THEN
      lw3=izets(IX,J-1)
      lw4=izete(IY,J-1)
      n2=lw4-lw3+1
      if(n2 .LT. 0) n2=0
      ENDIF
      
      INUM4(I)=INUM4(I)+n1
      IF(ICON .EQ. 1) THEN
      DO k=1,n1 
      IA(lw1+k-1)=it_recv_t(it+k,I)
      enddo
      ELSE
      DO k=1,n1
      II=iwet(lw1+k-1)
      IB(II,J,1)=it_recv_t(it+k,I)
      ENDDO
      ENDIF
      IF(NT .EQ. 2) THEN
      it=INUM4(I)
      INUM4(I)=INUM4(I)+n2
      IF(ICON .EQ. 1) THEN
      DO k=1,n2 
      IA(lw3+k-1)=it_recv_t(it+k,I)
      enddo
      ELSE
      DO k=1,n2
      II=iwet(lw3+k-1)
      IB(II,J-1,1)=it_recv_t(it+k,I)
      ENDDO
      ENDIF
      ENDIF
      ENDDO
      ENDIF 
 2    continue      
      return
      end

******************************************************************************************************************************************************
      subroutine IO_daily_physics_p(nout,mjar,lmon,nday,ivier,ip,imal) 
      include 'C_model.f'
      include 'C_daily_physics.f'
      include 'C_mpi.f'
      include 'C_index.f'
      parameter (khor1=khor+1)
      common /lbkhor/ khorl,lzet(khor),lb0(n),le0(n),
     * lb1(n),le1(n),lb2(n),le2(n)
      common /radis/fqgmit(m,n),fqrmit(khor),fqsmit(khor),fqlmit(khor)
      common /cord/ I1,I2,J1,J2,IA1,IA2,JA1,JA2,IB1,IB2,JB1,JB2
     1,ID1,ID2,JD1,JD2
      common /ind/ iwet(khor1),ldep(khor),lazc(khor),indend(n),
     *             isornr(n),isorsr(n),islab(n)

      if(myid .eq. 0)write (nout) mjar,lmon,nday,ivier,ip

      DO I=1,khorl
      lw=lzet(I)
      write (nout) zmit(lw),fqrmit(lw),fqsmit(lw),fqlmit(lw)
      enddo

      do J=J1,J2
      lwa=lb0(J)
      lwe=le0(J)
      if(lwe .ge. lwa)then
      nwet1=indwet(lwa)+1
      nwet2=indwet(lwe)+lazc(lwe)
      do k=nwet1,nwet2
      write (nout) umit(k),vmit(k),wcmit(k),acmit(k),
     &szmit(k),tcmit(k),scmit(k)
      enddo
      endif
      enddo

      do j=J1,J2
      do i=I1,I2
      write(nout) fqgmit(i,j),frimit(i,j),hismit(i,j),
     & hisrmit(i,j),tismit(i,j),uimit(i,j),vimit(i,j)
     & ,qois(i,j),qiis(i,j)
      enddo
      enddo
      return
      end
*************************************************************************************************************************************************
C remaining routines are for parallel output they are written by different procerssors in various files and commbined in one
*************************************************************************************************************************************************
      subroutine IO_start_physics_p(nout,mjar,lmon,nday,ivier,ip,imal)
     
      include 'C_model.f'
      parameter (nx=m,ny=n,ilop1=ilo+1,khor1=khor+1)
      include 'C_ice.f'
      common /urand/ zday(lrp),zvar(lrp),zalt(lrp)
      include 'C_mpi.f'
      include 'C_index.f'
      common /lbkhor/ khorl,lzet(khor),lb0(n),le0(n),
     * lb1(n),le1(n),lb2(n),le2(n)
      common /cord/ I1,I2,J1,J2,IA1,IA2,JA1,JA2,IB1,IB2,JB1,JB2
     1,ID1,ID2,JD1,JD2
          if(myid .eq. 0) write (nout) ip,mjar,lmon,nday
          if(myid .eq. 0) write (nout) zalt
          DO I=1,khorl
          lw=lzet(I)
          write (nout) zac(lw),txc(lw),tyc(lw),pac(lw)
          ENDDO
          DO 1 J=J1,J2
          lwa=lb0(J)
          lwe=le0(J)
          if(lwe .ge. lwa)then
          nwet1=indwet(lwa)+1
          nwet2=indwet(lwe)+lazc(lwe)
          do k=nwet1,nwet2
          write (nout) uc(k),vc(k),wc(k),avc(k),tec(k),sac(k)
          enddo
          endif
 1        continue
          do 2 j=J1,J2
          do 2 i=I1,I2
          write (nout) z(i,j),ui(i,j),vi(i,j),frice(i,j),
     &    his(i,j),hisr(i,j),tis(i,j)
 2        continue         
       return
       end
*************************************************************************************************************************************************
      subroutine combine_output_daily(nout,pathout,monat,ppp,
     & mjar,icord,itagend,nout_start,lmon,toutsteps)
      include 'C_model.f'
      include 'C_paths.f'
      include 'C_mpi.f'
      parameter (khor1=khor+1)
      dimension A(khor),B(khor),C(khor),D(khor),E(lrp)
      dimension A1(m,n),A2(m,n),A3(m,n),A4(m,n),A5(m,n),A6(m,n),
     &          A7(m,n),A8(m,n),A9(m,n)
      dimension B1(ndrei),B2(ndrei),B3(ndrei),B4(ndrei),B5(ndrei),
     &          B6(ndrei),B7(ndrei)
      character outfile*10
      character ppp*3
      character*19 pathout
      include 'C_index.f'
      common /lbkhor/ khorl,lzet(khor),lb0(n),le0(n),
     * lb1(n),le1(n),lb2(n),le2(n)
      common /ind/ iwet(khor1),ldep(khor),lazc(khor),indend(n),
     *             isornr(n),isorsr(n),islab(n)
      common /surface/ izets(m,n),izete(m,n)
      dimension icord(nprocs,4)
      common /out11/ mzet(khor,nprocs),lkhor(nprocs)

      outfile(1:3)=ppp
      write(outfile(4:5),'(i2.2)') mjar
c      write(outfile(6:7),'(i2.2)') monat

      do nn=1,nprocs
      write(outfile(6:8),'(i3.3)') nn-1
      open(unit=600+nn,file=pathpc2//pathout//outfile(1:8)
     & ,form='unformatted')
      enddo
      do im=1,itagend*toutsteps
      DO nn=1,nprocs
      if( nn .eq. 1) read(600+nn)mjar1,lmon1,nday1,ivier1,ip1 

      do I=1,lkhor(nn)
      lw=mzet(I,nn)
      read (600+nn) A(lw),B(lw),C(lw),D(lw)
      enddo
c      write(*,*)im,nn,lkhor(nn)

      I11=icord(nn,1)
      I22=icord(nn,2)
      J11=icord(nn,3)
      J22=icord(nn,4)
c      write(*,*)I11,I22,J11,J22
      do J=J11,J22
      lwa=izets(I11,J)
      lwe=izete(I22,J)
      if(lwe .ge. lwa)then
      nwet1=indwet(lwa)+1
      nwet2=indwet(lwe)+lazc(lwe)
      do k=nwet1,nwet2
      read (600+nn) B1(k),B2(k),B3(k),B4(k),
     &B5(k),B6(k),B7(k)
      enddo
      endif
      enddo

      do j=J11,J22
      do i=I11,I22
      read(600+nn) A1(i,j),A2(i,j),A3(i,j),
     & A4(i,j),A5(i,j),A6(i,j),A7(i,j)
     & ,A8(i,j),A9(i,j)
      enddo
      enddo
      enddo
      write (nout) mjar1,lmon1,nday1,ivier1,ip1
      write (nout) A 
      write (nout) B1
      write (nout) B2
      write (nout) B3
      write (nout) B4
      write (nout) B5 
      write (nout) B6
      write (nout) B7
      write(nout)  A2
      write(nout)  A3
      write(nout)  A4
      write(nout)  A5
      write(nout)  A6
      write(nout)  A7
      write(nout)  A1
      write(nout)  B 
      write(nout)  C 
      write(nout)  D
      write(nout)  A8
      write(nout)  A9
      enddo

      DO nn=1,nprocs
      if(nn .eq. 1) then
      read (600+nn) ip1,mjar1,lmon1,nday1
      read (600+nn) E
      endif
      do I=1,lkhor(nn)
      lw=mzet(I,nn)
      read (600+nn) A(lw),B(lw),C(lw),D(lw)
      enddo

      I11=icord(nn,1)
      I22=icord(nn,2)
      J11=icord(nn,3)
      J22=icord(nn,4)
      do 1 J=J11,J22
      lwa=izets(I11,J)
      lwe=izete(I22,J)
      if(lwe .ge. lwa)then
      nwet1=indwet(lwa)+1
      nwet2=indwet(lwe)+lazc(lwe)
      do k=nwet1,nwet2
      read (600+nn)B1(k),B2(k),B3(k),B4(k),B5(k),B6(k)
      enddo
      endif
 1    continue

      do 2 j=J11,J22
      do 2 i=I11,I22
      read (600+nn) A1(i,j),A2(i,j),A3(i,j),A4(i,j),
     &A5(i,j),A6(i,j),A7(i,j)
 2    continue
      enddo
      write (nout) ip1,mjar1,lmon1,nday1
      write (nout) E,A1,A
      write (nout) B,C,D
      write (nout) B1
      write (nout) B2
      write (nout) B3
      write (nout) B4
      write (nout) B5
      write (nout) B6
      write (nout) A4,A5,A6,A7
      write (nout) A2,A3
      if(lmon .eq. 12)then
      write (nout_start) ip1,mjar1,lmon1,nday1
      write (nout_start) E,A1,A
      write (nout_start) B,C,D
      write (nout_start) B1
      write (nout_start) B2
      write (nout_start) B3
      write (nout_start) B4
      write (nout_start) B5
      write (nout_start) B6
      write (nout_start) A4,A5,A6,A7
      write (nout_start) A2,A3
c      INQUIRE(UNIT=nout_start, POS=mypos)
c      PRINT *, " write start files position end", mypos
      endif
      do nn=1,nprocs
      close(600+nn)
      enddo

!      do j=1,n
!      do i=1,m
!      A1(i,j)=1.0
!      enddo
!      enddo
!      lwe=0
!      do j=1,n
!      lwa=lwe+1
!      lwe=indend(j)
!      do lw=lwa,lwe
!      nwet=indwet(lw)
!      i=iwet(lw)
!      A1(i,j)=0.0
!      do k=1,lazc(lw)
!      A1(i,j)=A1(i,j)+B3(nwet+k)
!      enddo
!      enddo
!      enddo
!      do j=1,n
!      write(71,*)(A1(i,j),i=1,m)
!      enddo
!      write(71) A1
      return
      end

***********************************************************************************************************************************************
      subroutine combine_output_start(nout,pathout,monat,ppp,
     & mjar,icordn)
      include 'C_model.f'
      include 'C_paths.f'
      include 'C_mpi.f'
      include 'C_index.f'
      parameter (khor1=khor+1)
      common /lbkhor/ khorl,lzet(khor),lb0(n),le0(n),
     * lb1(n),le1(n),lb2(n),le2(n)
      dimension A(khor),B(khor),C(khor),D(khor),E(lrp)
      dimension A1(ndrei),A2(ndrei),A3(ndrei),A4(ndrei),
     *A5(ndrei),A6(ndrei)
      dimension B1(m,n),B2(m,n),B3(m,n),B4(m,n),B5(m,n),B6(m,n),B7(m,n)
      character ppp*3,outfile*10
      character*19 pathout,command
      common /ind/ iwet(khor1),ldep(khor),lazc(khor),indend(n),
     *             isornr(n),isorsr(n),islab(n)
      common /surface/ izets(m,n),izete(m,n)
      dimension icord(nprocs,4)
      common /out11/ mzet(khor,nprocs),lkhor(nprocs)

      outfile(1:3)=ppp
      write(outfile(4:5),'(i2.2)') mjar
c      write(outfile(6:7),'(i2.2)') monat

      DO nn=1,nprocs
      write(outfile(6:8),'(i3.3)') nn-1
      open(unit=600+nn,file=pathpc2//pathout//outfile(1:8)
     &,form='unformatted')
      ENDDO

      DO nn=1,nprocs
      if(nn .eq. 1) then
      read (600+nn) ip1,mjar1,lmon1,nday1
      read (600+nn) E
      endif
      do I=1,lkhor(nn)
      lw=mzet(I,nn)
      read (600+nn) A(lw),B(lw),C(lw),D(lw)
      enddo 

      I11=icord(nn,1)
      I22=icord(nn,2)
      J11=icord(nn,3)
      J22=icord(nn,4)
      do 1 J=J11,J22
      lwa=izets(I11,J)
      lwe=izete(I22,J)
      if(lwe .ge. lwa)then
      nwet1=indwet(lwa)+1
      nwet2=indwet(lwe)+lazc(lwe)
      do k=nwet1,nwet2
      read (600+nn)A1(k),A2(k),A3(k),A4(k),A5(k),A6(k)
      enddo
      endif
 1    continue

      do 2 j=J11,J22
      do 2 i=I11,I22
      read (600+nn) B1(i,j),B2(i,j),B3(i,j),B4(i,j),
     &B5(i,j),B6(i,j),B7(i,j)
 2    continue
      enddo
      write (nout) ip1,mjar1,lmon1,nday1
      write (nout) E,B1,A
      write (nout) B,C,D
      write (nout) A1
      write (nout) A2
      write (nout) A3
      write (nout) A4
      write (nout) A5
      write (nout) A6
      write (nout) B4,B5,B6,B7
      write (nout) B2,B3
      do nn=1,nprocs
      close(600+nn)
c      write(outfile(8:9),'(i2.2)') nn-1
c      command= pathpc2//pathout//outfile(1:9)
c      write(*,*)command
c      call system(rm command)
      enddo
      return
      end

**************************************************************************************************************************************************
       subroutine combine_output_bio(nout,pathout,monat,ppp,
     & mjar,icord,itagend,nout_start,lmon,toutsteps)
      include 'C_model.f'
      include 'C_paths.f'
      include 'C_mpi.f'
      include 'C_index.f'
      include 'C_nbio_fabm.f'
      parameter (khor1=khor+1)
      common /lbkhor/ khorl,lzet(khor),lb0(n),le0(n),
     * lb1(n),le1(n),lb2(n),le2(n)
      dimension A(ndrei,3:nbio+nsed)  !dimension need to be npel make sure to change here hardcoded
      character ppp*3,outfile*10
      character*19 pathout
      common /ind/ iwet(khor1),ldep(khor),lazc(khor),indend(n),
     *             isornr(n),isorsr(n),islab(n)
      common /surface/ izets(m,n),izete(m,n)
      dimension icord(nprocs,4)
      common /out11/ mzet(khor,nprocs),lkhor(nprocs)

      outfile(1:3)=ppp
      write(outfile(4:5),'(i2.2)') mjar
c      write(outfile(6:7),'(i2.2)') monat 
      outfile(9:9)='b'

      DO nn=1,nprocs
      write(outfile(6:8),'(i3.3)') nn-1
      open(unit=600+nn,file=pathpc2//pathout//outfile
     &,form='unformatted')
      ENDDO

      do I=1,itagend*toutsteps
      do nn=1,nprocs
      I11=icord(nn,1)
      I22=icord(nn,2)
      J11=icord(nn,3)
      J22=icord(nn,4)
      do ib=3,nbio+nsed
      do J=J11,J22
      lwa=izets(I11,J)
      lwe=izete(I22,J)
      if(lwe .ge. lwa) then
      nwet1=indwet(lwa)+1
      nwet2=indwet(lwe)+lazc(lwe)
      read(600+nn)(A(k,ib),k=nwet1,nwet2)
      endif
      enddo
      enddo
      enddo
      write (nout) A
      enddo
      
      do nn=1,nprocs
      I11=icord(nn,1)
      I22=icord(nn,2)
      J11=icord(nn,3)
      J22=icord(nn,4)
      do ib=3,nbio+nsed
      do J=J11,J22
      lwa=izets(I11,J)
      lwe=izete(I22,J)
      if(lwe .ge. lwa) then
      nwet1=indwet(lwa)+1
      nwet2=indwet(lwe)+lazc(lwe)
      read(600+nn)(A(k,ib),k=nwet1,nwet2)
      endif
      enddo
      enddo
      enddo

      write (nout)A
      if(lmon .eq. 12) then
      write (nout_start)A
      endif
      do nn=1,nprocs
      close(600+nn)
      enddo

      return
      end

*************************************************************************************************************************************************
      subroutine combine_meteo_data(nanz,outfile1,pathout,nloop)
      include 'C_mpi.f'
      include 'C_model.f'
      include 'C_paths.f'
      character outfile1*10
      character*19 pathout
      dimension A(khor),B(khor),C(khor),D(khor),E(khor),F(khor),G(khor)
      common /out11/ mzet(khor,nprocs),lkhor(nprocs)
      do nn=1,nprocs
      write(outfile1(6:8),'(i3.3)') nn-1
      open(unit=600+nn,
     & file=pathpc2//pathout//'Qh_'//outfile1(1:8),form='unformatted')
      enddo

      DO J=1,nloop
      do nn=1,nprocs
      do I=1,lkhor(nn)
      k=mzet(I,nn)
      read(600+nn) A(k),B(k),C(k),D(k),
     &E(k),F(k),G(k)
      enddo
      enddo

      write(nanz) A,B,C
      write(nanz) D,E,F
      write(nanz) G
      ENDDO

      do nn=1,nprocs
      close(600+nn)
      enddo

      return
      end


**************************************************************************
