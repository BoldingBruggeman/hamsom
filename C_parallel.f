      common /topo1/ izet(m,n),ltief(m,n)
      common /bcor/ JWB1(2),JWB2(2),IWB1(2),IWB2(2),NWB
     &,JNB1(2),JNB2(2),INB1(2),INB2(2),NNB
      common /cord/ I1,I2,J1,J2,IA1,IA2,JA1,JA2,IB1,IB2,JB1,JB2
     1,ID1,ID2,JD1,JD2
      common /icecor/ ICEI1,ICEI2,ICEJ1,ICEJ2,ICEIA1,ICEIA2,ICEJA1,
     * ICEJA2,ICEX1,ICEX2,ICEY1,ICEY2,ICEP1,ICEP2,ICEQ1,ICEQ2,ICEX3,
     * ICEX4,ICEY3,ICEY4,ICEX5,ICEX6,ICEY5,ICEY6
      common /vcord/ nvrlayer(4),nvllayer(4),nvrvar(4*n,3),nvlvar(4*n,3)
     *,nvrdep(4),nvldep(4)
      common /surface/ izets(m,n),izete(m,n)
      common /hcordb2/ idepb2(nprocs),icdepb2(nprocs),
     & nbhlayer2(nprocs*6,7),nsendb2(nprocs),ncsendb2(nprocs)
      common /hcordt2/ idept2(nprocs),icdept2(nprocs),
     &nsendt2(nprocs),ncsendt2(nprocs),nthlayer2(nprocs*6,7)
      common /hcordb1/ idepb1(nprocs),icdepb1(nprocs),
     &nbhlayer1(nprocs*6,7),nsendb1(nprocs),ncsendb1(nprocs)
      common /hcordt1/ idept1(nprocs),icdept1(nprocs),
     &nthlayer1(nprocs*6,7),nsendt1(nprocs),ncsendt1(nprocs)
      common /hcordb0/ idepb0(nprocs),icdepb0(nprocs),
     &nbhlayer0(nprocs*6,7),nsendb0(nprocs),ncsendb0(nprocs)
      common /hcordt0/ idept0(nprocs),icdept0(nprocs),
     &nthlayer0(nprocs*6,7),nsendt0(nprocs),ncsendt0(nprocs)
      common /lbkhor/ khorl,lzet(khor),lb0(n),le0(n),
     * lb1(n),le1(n),lb2(n),le2(n)
      common /verexc/ temp1(10*N*ilo),temp2(10*N*ilo),temp3(10*N*ilo),
     &temp4(10*N*ilo)
      common /initc/ INUM1,INUM2,INUM3(6),INUM4(6)
     &,IRNUM1,IRNUM2,IRNUM3(6),IRNUM4(6),nreq1,nreq2,
     & nreq3,nreq4  
      common /horexc/t_recv_b(m*ilo*nPh,6),t_send_b(m*ilo*nPh,6),
     2 t_recv_t(m*ilo*nPh,6),t_send_t(m*ilo*nPh,6)
      common /hreq/ir_r_b(50),ir_s_b(6),ir_r_t(6),ir_s_t(6)
      common /isor/ isornr1(n),isorsr1(n),isornr2(n),isorsr2(n)
     1,ik1(khor),norsor1
      common /iverexc/ itemp1(2*N),itemp2(2*N),itemp3(2*N),
     &itemp4(2*N)
      common /ihorexc/ it_recv_b(m*nPh,6),it_send_b(m*nPh,6),
     2 it_recv_t(m*nPh,6),it_send_t(m*nPh,6)
      common /extra/ ctemp,ictemp,dtemp(m,n),idtemp(m,n),iel,iel1
      common /MAX11/ rumax(nprocs),rvmax(nprocs),rwmax(nprocs)
      common /out11/ mzet(khor,nprocs),lkhor(nprocs) 
