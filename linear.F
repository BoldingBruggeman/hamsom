      subroutine cgsolver(cc,stuvc,z)
      include 'C_model.f'
      parameter(khor1=khor+1)
      real limit,sumr0,sumr1,alpha,beta
      integer niter,norsor
      dimension z(m,n),stuvc(khor),z1(m,n),zd(kasor)
      dimension cc(kasor,9),r0(kasor),pr0(m,n),par0(m,n)
      integer*2 isor,ksor
      common /intit/ itel(1440),jc(m,3),ikk(kasor),isor(kasor),
     cksor(kasor),nofrez,lone,itende,nfreez
      common /ind/ iwet(khor1),ldep(khor),lazc(khor),indend(n),
     *             isornr(n),isorsr(n),islab(n)
      include 'C_mpi.f'
      include 'C_parallel.f'
      include 'C_index.f'

      eps    = 1.e-7
      sumr0=0.0           

       do it=1,norsor1
        ik=ik1(it)
        i = isor(ik) 
        k = ksor(ik)
        r0(ik)=dble(stuvc(ik)) -z(i,k)
     1 +dprod(cc(ik,1),z(i-1,k-1))+dprod(cc(ik,2),z(i-1,k)) 
     2 +dprod(cc(ik,3),z(i-1,k+1)) 
     3 +dprod(cc(ik,4),z(i,k-1))+dprod(cc(ik,5),z(i,k+1)) 
     4 +dprod(cc(ik,6),z(i+1,k-1)) 
     5 +dprod(cc(ik,7),z(i+1,k))+dprod(cc(ik,8),z(i+1,k+1))
        pr0(i,k)=r0(ik)
        sumr0=sumr0+r0(ik)*r0(ik)
       z1(i,k)=z(i,k)
      enddo
       sumr1=0.0
       call MPI_ALLREDUCE(sumr0,sumr1,1,MPI_REAL,MPI_SUM,
     3 MPI_COMM_WORLD,IERR)

      call INTI_COMM_VAR
      call com_vertical(ctemp,iel1,pr0,iel,1,1,0,2,1,lazc)
      call com_horizontal(ctemp,iel1,pr0,iel,1,1,0,2,1,icdepb1,
     *idepb1,icdept1
     *,idept1,ncsendb1,nsendb1,ncsendt1,nsendt1,nbhlayer1,nthlayer1)

      
!      call MPI_Waitall(nreq1, ir_r_b, status, ierr)
!      INUM1=0
!      INUM2=0

!      call VER_UNPACK(ctemp,iel1,pr0,iel,1,1,0,2,lazc)

!      INUM3=0
!      INUM4=0
!      call HOR_UNPACK(ctemp,iel1,pr0,iel,1,1,nbhlayer1,nthlayer1,
!     &idepb1,idept1,icdepb1,icdept1,0,2,iwet)

       DO iter=1,500
       INUM1=0
       INUM2=0
       INUM3=0
       INUM4=0
       cxcmax1=0.0 
       sumpar2=0.0
       sumpar0=0.0
       call MPI_Waitall(nreq1, ir_r_b, status, ierr)
       call VER_UNPACK(ctemp,iel1,pr0,iel,1,1,0,2,lazc)
       call HOR_UNPACK(ctemp,iel1,pr0,iel,1,1,nbhlayer1,nthlayer1,
     &idepb1,idept1,icdepb1,icdept1,0,2,iwet)

       do it=1,norsor1
       ik=ik1(it) 
       i = isor(ik) 
       k = ksor(ik) 
       par0(i,k)=-pr0(i,k)
     1 +dprod(cc(ik,1),pr0(i-1,k-1))+dprod(cc(ik,2),pr0(i-1,k)) 
     2 +dprod(cc(ik,3),pr0(i-1,k+1)) 
     3 +dprod(cc(ik,4),pr0(i,k-1))+dprod(cc(ik,5),pr0(i,k+1)) 
     4 +dprod(cc(ik,6),pr0(i+1,k-1)) 
     5 +dprod(cc(ik,7),pr0(i+1,k))+dprod(cc(ik,8),pr0(i+1,k+1))
       sumpar2=sumpar2+par0(i,k)*pr0(i,k)
       enddo           

       call MPI_ALLREDUCE(sumpar2,sumpar0,1,MPI_REAL,MPI_SUM,
     3 MPI_COMM_WORLD,IERR) 

       alpha=sumr1/sumpar0
       sumr0=0.0
       sumr2=0.0

       do it=1,norsor1
       ik=ik1(it)
       i = isor(ik) 
       k = ksor(ik)
       r0(ik)=r0(ik)-alpha*par0(i,k)
       sumr2=sumr2+r0(ik)*r0(ik) 
       z(i,k)=z(i,k)-alpha*pr0(i,k)     
       zd(ik)=z(i,k)-z1(i,k)
       cxcmax1 = max (cxcmax1, abs (zd(ik)))
       enddo

       call MPI_ALLREDUCE(sumr2,sumr0,1,MPI_REAL,MPI_SUM,
     3 MPI_COMM_WORLD,IERR)
       call MPI_ALLREDUCE(cxcmax1,cxcmax,1,MPI_REAL,MPI_MAX,
     3 MPI_COMM_WORLD,IERR)
       beta=sumr0/sumr1
c       if(myid .eq. 0) write(*,*)cxcmax,eps
       if (cxcmax .lt. eps) goto 1
       sumr1=sumr0

       do it=1,norsor1
       ik=ik1(it)
       i = isor(ik) 
       k = ksor(ik)
       pr0(i,k)=r0(ik)+beta*pr0(i,k)
       z1(i,k)=z(i,k)
       enddo

      call INTI_COMM_VAR
      call com_vertical(ctemp,iel1,pr0,iel,1,1,0,2,1,lazc)
      call com_horizontal(ctemp,iel1,pr0,iel,1,1,0,2,1,icdepb1,
     *idepb1,icdept1
     *,idept1,ncsendb1,nsendb1,ncsendt1,nsendt1,nbhlayer1,nthlayer1)
       ENDDO

 1     continue

       if( iter .ge. 500) then
       write(*,*) 'sor not converging'
       write(*,*)cxcmax,eps
       STOP
       endif

       return 
       END
