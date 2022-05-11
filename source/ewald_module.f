      module ewald_module
      
c***********************************************************************
c     
c     dl_poly module for defining ewald sum arrays
c     copyright - daresbury laboratory
c     author    - w. smith    sep 2003
c     adapted for solvation, free energy and excitation
c     - p.-a. cazade oct 2007
c     
c***********************************************************************
      
      use config_module
      use error_module
      use exclude_module
      use metafreeze_module
      use pair_module
      use property_module
      use setup_module
      use solvation_module
      
      implicit none
      
      real(8), allocatable :: ckc(:),cks(:),clm(:),slm(:)
      real(8), allocatable :: elc(:,:),els(:,:)
      real(8), allocatable :: emc(:,:),ems(:,:)
      real(8), allocatable :: enc(:,:),ens(:,:)
      real(8), allocatable :: ewlbuf(:),erc(:),fer(:)
      
      save ckc,cks,clm,slm,elc,emc,enc,els,ems,ens,erc,fer,ewlbuf
      
      contains
      
      subroutine alloc_ewald_arrays(idnode,mxnode)
      
      implicit none
      
      integer, parameter :: nnn=6
      
      logical safe
      integer i,fail,idnode,mxnode
      dimension fail(nnn)
      
      safe=.true.

c     allocate arrays
      
      fail(:)=0
      
      allocate (ckc(mxewld),cks(mxewld),stat=fail(1))
      allocate (clm(mxewld),slm(mxewld),stat=fail(2))
      allocate (elc(mxewld,0:1),els(mxewld,0:1),stat=fail(3))
      allocate (emc(mxewld,0:kmaxb),ems(mxewld,0:kmaxb),stat=fail(4))
      allocate (enc(mxewld,0:kmaxc),ens(mxewld,0:kmaxc),stat=fail(5))
      allocate (ewlbuf(mxebuf),erc(mxegrd),fer(mxegrd),stat=fail(6))

      if(any(fail.gt.0))safe=.false.      
      if(mxnode.gt.1)call gstate(safe)    
      if(.not.safe)call error(idnode,1240)
      
      end subroutine alloc_ewald_arrays
      
      subroutine erfcgen(alpha,drewd,rcut)
      
c***********************************************************************
c     
c     dlpoly routine for generating interpolation tables for 
c     erfc and its derivative - for use with ewald sum.
c     
c     copyright daresbury laboratory 1994
c     author t.forester dec 1994
c     
c***********************************************************************
      
      implicit none
      
      integer i
      real(8) alpha,drewd,rcut,a1,a2,a3,a4,a5,pp,rrr,rsq,tt,exp1
      
      data a1,a2,a3/0.254829592d0,-0.284496736d0,1.421413741d0/
      data a4,a5,pp/-1.453152027d0,1.061405429d0,0.3275911d0/
      
c     look-up tables for real space part of ewald sum
      
      drewd=rcut/dble(mxegrd-4)
      
      do i=1,mxegrd
        
        rrr=dble(i)*drewd
        rsq=rrr*rrr
        tt=1.d0/(1.d0+pp*alpha*rrr)
        exp1=exp(-(alpha*rrr)**2)
        erc(i)=tt*(a1+tt*(a2+tt*(a3+tt*(a4+tt*a5))))*exp1/rrr
        fer(i)=(erc(i)+2.d0*(alpha/sqrpi)*exp1)/rsq
        
      enddo
      
      return
      end subroutine erfcgen
      
      subroutine ewald1
     x  (lsolva,llsolva,lfree,lghost,idnode,mxnode,natms,imcon,
     x  kmax1,kmax2,kmax3,engcpe,vircpe,alpha,volm,epsq)
      
c***********************************************************************
c     
c     dl_poly subroutine for calculating coulombic forces in a
c     periodic system using ewald's method
c     
c     parallel replicated data version (part 1)
c     
c     copyright - daresbury laboratory
c     author    - w. smith march 1992.
c     modified  - t. forester april 1993: t3d adaptation
c     modified  - p.-a. cazade oct 2007: solvation, free energy etc
c     
c     part 1 - reciprocal space terms (fourier part)
c     
c     note - in loop over all k vectors k=2pi(ll/cl,mm/cl,nn/cl)
c     the values of ll,mm and nn are selected so that the symmetry of
c     reciprocal lattice is taken into account i.e. the following
c     rules apply.
c     
c     ll ranges over the values 0 to kmax1 only.
c     
c     mm ranges over 0 to kmax2 when ll=0 and over
c     -kmax2 to kmax2 otherwise.
c     nn ranges over 1 to kmax3 when ll=mm=0 and over
c     -kmax3 to kmax3 otherwise.
c     
c     hence the result of the summation must be doubled at the end.
c     
c     stress tensor added t.forester may 1994
c     
c***********************************************************************
      
      implicit none
      
      logical newjob,lconsw,safe,leven,lsolva,llsolva,lfree,lghost
      integer isol,jsol,ksol,iisol,jjsol,kksol,kstep
      integer idnode,mxnode,natms,imcon,kmax1,kmax2,kmax3
      integer iatm0,iatm1,i,j,limit,l,npass,ipass,kkk,nmin
      integer mmin,ll,m,mm,n,nn
      real(8) engcpe,vircpe,alpha,volm,epsq,omg,qchg,qfix,qforce
      real(8) twopi,rvolm,ralph,det,rcpcut,rcpct2,engsic,ssx
      real(8) ssy,ssz,rkx1,rky1,rkz1,cs,rkx2,rky2,rkz2,eng1
      real(8) rkx3,rky3,rkz3,rksq,ckcs,ckss,rrksq,akk,bkk,akv
      real(8) scal1,scale,virprs,ckc1s,cks1s,ckc2s,cks2s,fkk
      real(8) term1a,term2a,term1b,term2b
      
      dimension omg(9)
      
      save newjob,engsic,qchg
      
      data newjob/.true./,lconsw/.true./,safe/.true./,leven/.true./
      
      twopi=2.d0*pi
      
      kstep=2
      if(lfree.or.lghost)kstep=6
      
      if(alpha.lt.1.d-6)return
      if(mxewld.ne.msatms)call error(idnode,330)
      
c     set up atoms numbers for nodes
      
      iatm0=(idnode*natms)/mxnode+1
      iatm1=((idnode+1)*natms)/mxnode
      
c     initialise coulombic potential energy
      
      engcpe=0.d0
      vircpe=0.d0
      
c     initalize stress tensor working arrays
      
      do i=1,9
        omg(i)=0.d0
      enddo
      
c     set working parameters
      
      rvolm=twopi/volm
      ralph=-0.25d0/alpha**2
      
c     set switch for TO, RD and HP boundary conditions
      
      if(imcon.eq.4.or.imcon.eq.5.or.imcon.eq.7)then
        
        lconsw=.false.
        rvolm=0.5d0*rvolm
        
      endif
      
c     construct reciprocal lattice vectors and set k vector range
      
      call invert(cell,rcell,det)
      if(abs(det).lt.1.d-6)call error(idnode,120)
      call dcell(rcell,buffer)
      
      rcpcut=min(dble(kmax1)*buffer(7),dble(kmax2)*buffer(8),
     x  dble(kmax3)*buffer(9))
      rcpcut=rcpcut*1.05d0*twopi
      rcpct2=rcpcut**2
      
      if(newjob)then
        
c     calculate self interaction correction (sic) and net system charge
        
        qchg=0.d0
        qchg0=0.d0
        qchg1=0.d0
        qchg2=0.d0
        engsic=0.d0
        engsic0=0.d0
        engsic2=0.d0
        
c     set solvation charge correction variables
        
        if(lsolva)then
          
          qfix_sol(:)=0.d0
          cou_sol_sic(:)=0.d0
          if(lghost)then
            
            qfix_exc(:)=0.d0
            cou_exc_sic(:)=0.d0
            
          endif
          
        endif
        
        if(lghost)then
          
c     set excitation sic and charge correction variables
          
          do i=iatm0,iatm1
            
            if(atm_fre(i).ne.2)then
              qchg0=qchg0+chge(i)
              engsic0=engsic0+chge(i)**2
            endif
            if(atm_fre(i).ne.1)then
              qchg2=qchg2+chge(i)
              engsic2=engsic2+chge(i)**2
            endif
            
          enddo
          
        elseif(lfree)then
          
c     set free energy sic and charge correction variables
          
          do i=iatm0,iatm1
            
            if(atm_fre(i).eq.0)then
              qchg0=qchg0+chge(i)
              engsic0=engsic0+chge(i)**2
            elseif(atm_fre(i).eq.1)then
              qchg1=qchg1+chge(i)
              engsic0=engsic0+lambda1*chge(i)**2
              engsic2=engsic2-chge(i)**2
            elseif(atm_fre(i).eq.2)then
              qchg2=qchg2+chge(i)
              engsic0=engsic0+lambda2*chge(i)**2
              engsic2=engsic2+chge(i)**2
            endif
            
          enddo
          
        else
          
c     set normal sic and charge correction variables
          
          do i=iatm0,iatm1
            
            qchg=qchg+chge(i)
            engsic=engsic+chge(i)**2
            
          enddo
          
        endif
        
        if(lsolva)then
          
          if(lghost)then
            
c     set excitation sic and charge correction arrays
            
            do i=iatm0,iatm1
              
              kkk=loc2(atmolt(i),atmolt(i))
              if(atm_fre(i).ne.2)then
                cou_sol_sic(kkk)=cou_sol_sic(kkk)+chge(i)**2
                qfix_sol(atmolt(i))=qfix_sol(atmolt(i))+chge(i)
              endif
              if(atm_fre(i).ne.1)then
                cou_exc_sic(kkk)=cou_exc_sic(kkk)+chge(i)**2
                qfix_exc(atmolt(i))=qfix_exc(atmolt(i))+chge(i)
              endif
              
            enddo
            
          else
            
c     set solvation sic and charge correction arrays
            
            do i=iatm0,iatm1
              
              kkk=loc2(atmolt(i),atmolt(i))
              cou_sol_sic(kkk)=cou_sol_sic(kkk)+chge(i)**2
              qfix_sol(atmolt(i))=qfix_sol(atmolt(i))+chge(i)
              
            enddo
            
          endif
          
        endif
        
c     calculate global values for correction variables and arrays
        
        if(mxnode.gt.1)then
          
          buffer(1)=qchg
          buffer(2)=qchg0
          buffer(3)=qchg1
          buffer(4)=qchg2
          call gdsum(buffer(1),4,buffer(5))
          qchg=buffer(1)
          qchg0=buffer(2)
          qchg1=buffer(3)
          qchg2=buffer(4)
          if(lsolva)then
            
            call gdsum(qfix_sol(1),mxtmls,buffer(1))
            if(lghost)call gdsum(qfix_exc(1),mxtmls,buffer(1))
            
          endif
          
        endif
        
c     store self interaction correction terms
        
        engsic=-r4pie0/epsq*alpha*engsic/sqrpi
        engsic0=-r4pie0/epsq*alpha*engsic0/sqrpi
        engsic2=-r4pie0/epsq*alpha*engsic2/sqrpi
        
        if(lsolva)then
          
          cou_sol_sic(:)=-r4pie0/epsq*alpha*cou_sol_sic(:)/sqrpi
          if(lghost)cou_exc_sic(:)=-r4pie0/epsq*alpha*
     x      cou_exc_sic(:)/sqrpi
          
        endif
        
        newjob=.false.
        
      endif
      
      if(lfree.or.lghost)then
        
        qchg=qchg0
        engsic=engsic0
        
      endif
      
c     calculate and store exponential factors
c     convert real to reciprocal space coordinates
      
      i=0
      
      do j=iatm0,iatm1
        
        i=i+1
        elc(i,0)=1.d0
        emc(i,0)=1.d0
        enc(i,0)=1.d0
        els(i,0)=0.d0
        ems(i,0)=0.d0
        ens(i,0)=0.d0
        ssx=rcell(1)*xxx(j)+rcell(4)*yyy(j)+rcell(7)*zzz(j)
        ssy=rcell(2)*xxx(j)+rcell(5)*yyy(j)+rcell(8)*zzz(j)
        ssz=rcell(3)*xxx(j)+rcell(6)*yyy(j)+rcell(9)*zzz(j)
        elc(i,1)=cos(twopi*ssx)
        emc(i,1)=cos(twopi*ssy)
        enc(i,1)=cos(twopi*ssz)
        els(i,1)=sin(twopi*ssx)
        ems(i,1)=sin(twopi*ssy)
        ens(i,1)=sin(twopi*ssz)
        
      enddo
      
      limit=i
      
      do l=2,kmax2
        
        do i=1,limit
          
          emc(i,l)=emc(i,l-1)*emc(i,1)-ems(i,l-1)*ems(i,1)
          ems(i,l)=ems(i,l-1)*emc(i,1)+emc(i,l-1)*ems(i,1)
          
        enddo
        
      enddo
      
      do l=2,kmax3
        
        do i=1,limit
          
          enc(i,l)=enc(i,l-1)*enc(i,1)-ens(i,l-1)*ens(i,1)
          ens(i,l)=ens(i,l-1)*enc(i,1)+enc(i,l-1)*ens(i,1)
          
        enddo
        
      enddo
      
c     start of main loop over k vectors
      
      npass=1
      if(mxnode.gt.16)npass=2
      if((mxnode.gt.1).and.(mxebuf.gt.5000))npass=2
      
      do ipass=1,npass
        
        kkk=0
        mmin=0
        nmin=1
        if(llsolva)kksol=0
        
        do ll=0,kmax1
          
          l=ll
          rkx1=twopi*dble(ll)*rcell(1)
          rky1=twopi*dble(ll)*rcell(4)
          rkz1=twopi*dble(ll)*rcell(7)
          
c     put cos(i,L) terms into cos(i,0) array
          
          if(l.eq.1)then
            
            do i=1,limit
              
              elc(i,0)=elc(i,1)
              els(i,0)=els(i,1)
              
            enddo
            
          elseif(l.gt.1)then
            
            do i=1,limit
              
              cs=elc(i,0)
              elc(i,0)=cs*elc(i,1)-els(i,0)*els(i,1)
              els(i,0)=els(i,0)*elc(i,1)+cs*els(i,1)
              
            enddo
            
          endif
          
          do mm=mmin,kmax2
            
            m=iabs(mm)
            rkx2=rkx1+twopi*dble(mm)*rcell(2)
            rky2=rky1+twopi*dble(mm)*rcell(5)
            rkz2=rkz1+twopi*dble(mm)*rcell(8)
            
c     set temporary products of exponential terms
            
            if(mm.ge.0)then
              
              do i=1,limit
                
                clm(i)=elc(i,0)*emc(i,m)-els(i,0)*ems(i,m)
                slm(i)=els(i,0)*emc(i,m)+ems(i,m)*elc(i,0)
                
              enddo
              
            else
              
              do i=1,limit
                
                clm(i)=elc(i,0)*emc(i,m)+els(i,0)*ems(i,m)
                slm(i)=els(i,0)*emc(i,m)-ems(i,m)*elc(i,0)
                
              enddo
              
            endif
            
            do nn=nmin,kmax3
              
              n=iabs(nn)
              
              if(.not.lconsw)then
                
                if(imcon.eq.7)then
                  
                  leven=(mod(l+m,2).eq.0)
                  
                else
                  
                  leven=(mod(l+m+n,2).eq.0)
                  
                endif
                
              endif
              
              if(lconsw.or.leven)then
                
                rkx3=rkx2+twopi*dble(nn)*rcell(3)
                rky3=rky2+twopi*dble(nn)*rcell(6)
                rkz3=rkz2+twopi*dble(nn)*rcell(9)
                
c     test on magnitude of k vector
                
                rksq=rkx3*rkx3+rky3*rky3+rkz3*rkz3
                
                if(rksq.le.rcpct2)then
                  
c     calculate exp(ikr) terms and product with charges
                  
                  i=0
                  
                  if(nn.ge.0)then
                    
                    if(lfree.or.lghost)then
                      
                      do j=iatm0,iatm1
                        
                        i=i+1
                        if(atm_fre(j).eq.1)then
                          ckc1(i)=chge(j)*(clm(i)*enc(i,n)-
     x                      slm(i)*ens(i,n))
                          cks1(i)=chge(j)*(slm(i)*enc(i,n)+
     x                      clm(i)*ens(i,n))
                        elseif(atm_fre(j).eq.2)then
                          ckc2(i)=chge(j)*(clm(i)*enc(i,n)-
     x                      slm(i)*ens(i,n))
                          cks2(i)=chge(j)*(slm(i)*enc(i,n)+
     x                      clm(i)*ens(i,n))
                        else
                          ckc(i)=chge(j)*(clm(i)*enc(i,n)-
     x                      slm(i)*ens(i,n))
                          cks(i)=chge(j)*(slm(i)*enc(i,n)+
     x                      clm(i)*ens(i,n))
                        endif
                        
                      enddo
                      
                    else
                      
                      do j=iatm0,iatm1
                        
                        i=i+1
                        ckc(i)=chge(j)*(clm(i)*enc(i,n)-slm(i)*ens(i,n))
                        cks(i)=chge(j)*(slm(i)*enc(i,n)+clm(i)*ens(i,n))
                        
                      enddo
                      
                    endif
                    
                  else
                    
                    if(lfree.or.lghost)then
                      
                      do j=iatm0,iatm1
                        
                        i=i+1
                        if(atm_fre(j).eq.1)then
                          ckc1(i)=chge(j)*(clm(i)*enc(i,n)+
     x                      slm(i)*ens(i,n))
                          cks1(i)=chge(j)*(slm(i)*enc(i,n)-
     x                      clm(i)*ens(i,n))
                        elseif(atm_fre(j).eq.2)then
                          ckc2(i)=chge(j)*(clm(i)*enc(i,n)+
     x                      slm(i)*ens(i,n))
                          cks2(i)=chge(j)*(slm(i)*enc(i,n)-
     x                      clm(i)*ens(i,n))
                        else
                          ckc(i)=chge(j)*(clm(i)*enc(i,n)+
     x                      slm(i)*ens(i,n))
                          cks(i)=chge(j)*(slm(i)*enc(i,n)-
     x                      clm(i)*ens(i,n))
                        endif
                        
                      enddo
                      
                    else
                      
                      do j=iatm0,iatm1
                        
                        i=i+1
                        ckc(i)=chge(j)*(clm(i)*enc(i,n)+slm(i)*ens(i,n))
                        cks(i)=chge(j)*(slm(i)*enc(i,n)-clm(i)*ens(i,n))
                        
                      enddo
                      
                    endif
                    
                  endif
                  
                  if(ipass.eq.1)then
                    
c     calculate vector sums
                    
                    ckcs=0.d0
                    ckss=0.d0
                    
                    if(lfree.or.lghost)then
                      
                      ckc1s=0.d0
                      cks1s=0.d0
                      ckc2s=0.d0
                      cks2s=0.d0
                      
                    endif
                    
                    if(llsolva)then
                      
                      ckc_sol_sum(:)=0.d0
                      cks_sol_sum(:)=0.d0
                      
                      if(lghost)then
                        ckc_fre_sum(:)=0.d0
                        cks_fre_sum(:)=0.d0
                      endif
                      
                    endif
                    
                    if(lfree.or.lghost)then
                      
                      i=0
                      do j=iatm0,iatm1
                        
                        i=i+1
                        if(atm_fre(j).eq.1)then
                          ckc1s=ckc1s+ckc1(i)
                          cks1s=cks1s+cks1(i)
                        elseif(atm_fre(j).eq.2)then
                          ckc2s=ckc2s+ckc2(i)
                          cks2s=cks2s+cks2(i)
                        else
                          ckcs=ckcs+ckc(i)
                          ckss=ckss+cks(i)
                        endif
                        
                      enddo
                      
                    else
                      
                      do i=1,limit
                        
                        ckcs=ckcs+ckc(i)
                        ckss=ckss+cks(i)
                        
                      enddo
                      
                    endif
                    
                    if(llsolva)then
                      
                      i=0
                      if(lghost)then
                        
                        do j=iatm0,iatm1
                          
                          i=i+1
                          if(atm_fre(j).eq.1)then
                            ckc_sol_sum(atmolt(j))=
     x                        ckc_sol_sum(atmolt(j))+ckc1(i)
                            cks_sol_sum(atmolt(j))=
     x                        cks_sol_sum(atmolt(j))+cks1(i)
                          elseif(atm_fre(j).eq.2)then
                            ckc_fre_sum(atmolt(j))=
     x                        ckc_fre_sum(atmolt(j))+ckc2(i)
                            cks_fre_sum(atmolt(j))=
     x                        cks_fre_sum(atmolt(j))+cks2(i)
                          else
                            ckc_sol_sum(atmolt(j))=
     x                        ckc_sol_sum(atmolt(j))+ckc(i)
                            cks_sol_sum(atmolt(j))=
     x                        cks_sol_sum(atmolt(j))+cks(i)
                          endif
                          
                        enddo
                        
                      else
                        
                        do j=iatm0,iatm1
                          
                          i=i+1
                          ckc_sol_sum(atmolt(j))=ckc_sol_sum(atmolt(j))+
     x                      ckc(i)
                          cks_sol_sum(atmolt(j))=cks_sol_sum(atmolt(j))+
     x                      cks(i)
                          
                        enddo
                        
                      endif
                      
                    endif
                    
c     perform global summation of exp(ikr) terms or store if npass=2
                    
                    if(npass.eq.2)then
                      
                      if(kkk+kstep.le.mxebuf)then
                        
                        ewlbuf(kkk+1)=ckcs
                        ewlbuf(kkk+2)=ckss
                        
                        if(lfree.or.lghost)then
                          
                          ewlbuf(kkk+3)=ckc1s
                          ewlbuf(kkk+4)=cks1s
                          ewlbuf(kkk+5)=ckc2s
                          ewlbuf(kkk+6)=cks2s
                          
                        endif
                        
                        if(llsolva)then
                          
                          do isol=1,mxtmls
                            
                            ebuf_sol1(kksol+isol)=ckc_sol_sum(isol)
                            ebuf_sol2(kksol+isol)=cks_sol_sum(isol)
                            
                            if(lghost)then
                              
                              ebuf_exc1(kksol+isol)=ckc_fre_sum(isol)
                              ebuf_exc2(kksol+isol)=cks_fre_sum(isol)
                              
                            endif
                            
                          enddo
                          
                        endif
                        
                      else
                        
                        safe=.false.
                        
                      endif
                      
                    elseif(mxnode.gt.1)then
                      
                      buffer(1)=ckcs
                      buffer(2)=ckss
                      call gdsum(buffer(1),2,buffer(3))
                      ckcs=buffer(1)
                      ckss=buffer(2)
                      
                      if(lfree.or.lghost)then
                        
                        buffer(1)=ckc1s
                        buffer(2)=cks1s
                        buffer(3)=ckc2s
                        buffer(4)=cks2s
                        call gdsum(buffer(1),4,buffer(5))
                        ckc1s=buffer(1)
                        cks1s=buffer(2)
                        ckc2s=buffer(3)
                        cks2s=buffer(4)
                        
                      endif
                      
                      if(llsolva)then
                        
                        call gdsum(ckc_sol_sum(1),mxtmls,buffer(1))
                        call gdsum(cks_sol_sum(1),mxtmls,buffer(1))
                        
                        if(lghost)then
                          
                          call gdsum(ckc_fre_sum(1),mxtmls,buffer(1))
                          call gdsum(cks_fre_sum(1),mxtmls,buffer(1))
                          
                        endif
                        
                      endif
                      
                    endif
                    
                  endif
                  
                  if(ipass.eq.npass)then
                    
                    if(npass.eq.2)then
                      
                      ckcs=ewlbuf(kkk+1)
                      ckss=ewlbuf(kkk+2)
                      
                      if(lfree.or.lghost)then
                        
                        ckc1s=ewlbuf(kkk+3)
                        cks1s=ewlbuf(kkk+4)
                        ckc2s=ewlbuf(kkk+5)
                        cks2s=ewlbuf(kkk+6)
                        
                      endif
                      
                      if(llsolva)then
                        
                        do isol=1,mxtmls
                          
                          ckc_sol_sum(isol)=ebuf_sol1(kksol+isol)
                          cks_sol_sum(isol)=ebuf_sol2(kksol+isol)
                          
                          if(lghost)then
                            
                            ckc_fre_sum(isol)=ebuf_exc1(kksol+isol)
                            cks_fre_sum(isol)=ebuf_exc2(kksol+isol)
                            
                          endif
                          
                        enddo
                        
                      endif
                      
                    endif
                    
c     calculate akk coefficients
                    
                    rrksq=1.d0/rksq
                    if(lconsw)then
                      akk=exp(ralph*rksq)*rrksq
                    else
                      akk=4.0d0*exp(ralph*rksq)*rrksq
                    endif
                    bkk=akk
                    akv=2.d0*akk*(rrksq-ralph)
                    
c     accumulate potential energy and virial terms
                    
                    if(lghost)then
                      
                      engcpe=engcpe+akk*((ckcs*ckcs+ckss*ckss)
     x                  +(ckc1s*ckc1s+cks1s*cks1s)+2.d0*
     x                  (ckc1s*ckcs+cks1s*ckss))
                      virprs=akv*((ckcs*ckcs+ckss*ckss)
     x                  +(ckc1s*ckc1s+cks1s*cks1s)+2.d0*
     x                  (ckc1s*ckcs+cks1s*ckss))
                      
                    elseif(lfree)then
                      
                      term1a=(ckc1s*ckc1s+cks1s*cks1s)
                      term1b=2.d0*(ckc1s*ckcs+cks1s*ckss)
                      term2a=(ckc2s*ckc2s+cks2s*cks2s)
                      term2b=2.d0*(ckc2s*ckcs+cks2s*ckss)
                      
                      engcpe=engcpe+akk*
     x                  ((ckcs*ckcs+ckss*ckss)+
     x                  lambda1*(term1a+term1b)+
     x                  lambda2*(term2a+term2b))
                      
                      cou_fre=cou_fre+akk*
     x                  ((term2a+term2b)-(term1a+term1b))
                      
                      virprs=akv*
     x                  ((ckcs*ckcs+ckss*ckss)+
     x                  lambda1*(term1a+term1b)+
     x                  lambda2*(term2a+term2b))
                      
                      cou_vir=cou_vir+akv*rksq*
     x                  ((term2a+term2b)-(term1a+term1b))
                      
                    else
                      
                      engcpe=engcpe+akk*(ckcs*ckcs+ckss*ckss)
                      virprs=akv*(ckcs*ckcs+ckss*ckss)
                      
                    endif
                    
                    if(llsolva)then
                      
                      ksol=0
                      do isol=1,mxtmls
                        
                        fkk=2.d0*akk
                        do jsol=1,isol
                          
                          ksol=ksol+1
                          if(isol.eq.jsol)fkk=akk
                          
                          cou_sol(ksol)=cou_sol(ksol)+
     x                      fkk*(ckc_sol_sum(isol)*ckc_sol_sum(jsol)+
     x                      cks_sol_sum(isol)*cks_sol_sum(jsol))
                          
                          if(lghost)then
                            
                            cou_exc(ksol)=cou_exc(ksol)+
     x                        akk*(ckc_fre_sum(isol)*ckc_fre_sum(jsol)+
     x                        cks_fre_sum(isol)*cks_fre_sum(jsol))
                            
                          endif
                          
                        enddo
                        
                      enddo
                      
                    endif
                    
c     contributions to stress tensor
                    
                    omg(1)=omg(1)-virprs*rkx3*rkx3
                    omg(5)=omg(5)-virprs*rky3*rky3
                    omg(9)=omg(9)-virprs*rkz3*rkz3
                    omg(2)=omg(2)-virprs*rkx3*rky3
                    omg(3)=omg(3)-virprs*rkx3*rkz3
                    omg(6)=omg(6)-virprs*rky3*rkz3
                    
c     calculate force on each site
                    
                    i=0
                    if(lghost)then
                      
                      do j=iatm0,iatm1
                        
                        i=i+1
                        if(atm_fre(j).eq.2)then
                          qforce=bkk*(cks2(i)*(ckcs+ckc2s)-
     x                      ckc2(i)*(ckss+cks2s))
                        elseif(atm_fre(j).eq.1)then
                          qforce=bkk*(cks1(i)*(ckcs+ckc1s)-
     x                      ckc1(i)*(ckss+cks1s))
                        else
                          qforce=bkk*(cks(i)*(ckcs+ckc1s)-
     x                      ckc(i)*(ckss+cks1s))
                        endif
                        
                        fxx(j)=fxx(j)+rkx3*qforce
                        fyy(j)=fyy(j)+rky3*qforce
                        fzz(j)=fzz(j)+rkz3*qforce
                        
                      enddo
                      
                    elseif(lfree)then
                      
                      do j=iatm0,iatm1
                        
                        i=i+1
                        if(atm_fre(j).eq.0)then
                          qforce=bkk*(cks(i)*(ckcs+lambda1*ckc1s+
     x                      lambda2*ckc2s)-ckc(i)*(ckss+lambda1*
     x                      cks1s+lambda2*cks2s))
                        elseif(atm_fre(j).eq.1)then
                          qforce=lambda1*bkk*(cks1(i)*(ckcs+ckc1s)
     x                      -ckc1(i)*(ckss+cks1s))
                        elseif(atm_fre(j).eq.2)then
                          qforce=lambda2*bkk*(cks2(i)*(ckcs+ckc2s)
     x                      -ckc2(i)*(ckss+cks2s))
                        endif
                        
                        fxx(j)=fxx(j)+rkx3*qforce
                        fyy(j)=fyy(j)+rky3*qforce
                        fzz(j)=fzz(j)+rkz3*qforce
                        
                      enddo
                      
                    else
                      
                      do j=iatm0,iatm1
                        
                        i=i+1
                        qforce=bkk*(cks(i)*ckcs-ckc(i)*ckss)
                        fxx(j)=fxx(j)+rkx3*qforce
                        fyy(j)=fyy(j)+rky3*qforce
                        fzz(j)=fzz(j)+rkz3*qforce
                        
                      enddo
                      
                    endif
                    
c     end vector loop
                    
                  endif
                  
                  kkk=kkk+kstep
                  
                  if(llsolva)kksol=kksol+mxtmls
                  
                endif
                
              endif
              
            enddo
            
            nmin=-kmax3
            
          enddo
          
          mmin=-kmax2
          
        enddo
        
c     delayed global sum of exp(ikr) terms for npass=2 case
        
        if(ipass.eq.1.and.npass.eq.2)then
          
          if(safe)then
            
            call gdsum(ewlbuf,kkk,buffer)
            
            if(llsolva)then
              
              call gdsum(ebuf_sol1,kksol,buffer)
              call gdsum(ebuf_sol2,kksol,buffer)
              
              if(lghost)then
                
                call gdsum(ebuf_exc1,kksol,buffer)
                call gdsum(ebuf_exc2,kksol,buffer)
                
              endif
              
            endif
            
            do i=1,limit
              
              elc(i,0)=1.d0
              els(i,0)=0.d0
              
            enddo
            
          else
            
            if(idnode.eq.0)then
              
              write(nrite,'(a,i10)')
     x          'dimension of ewlbuf array required ',kkk
              write(nrite,'(a,i10)')
     x          'dimension of current  ewlbuf array ',mxebuf
              
            endif
            
            call error(idnode,46)
            
          endif
          
        endif
        
      enddo
      
c     reduce sums by factor mxnode for global summation
      
      engcpe=engcpe/dble(mxnode)
      
      if(lfree)then
        cou_fre=cou_fre/dble(mxnode)
        cou_vir=cou_vir/dble(mxnode)
      endif
      
      if(llsolva)then
        cou_sol(:)=cou_sol(:)/dble(mxnode)
        if(lghost)cou_exc(:)=cou_exc(:)/dble(mxnode)
      endif
      
      do i=1,9
        omg(i)=omg(i)/dble(mxnode)
      enddo
      
c     correction for charged systems
      
      if(lfree)then
        
        qfix=-(0.5d0*pi*r4pie0/epsq)*(((qchg/alpha)**2+lambda1*
     x    (qchg1/alpha)**2+lambda2*(qchg2/alpha)**2+2.d0*lambda1*
     x    (qchg1/alpha)*(qchg/alpha)+2.d0*lambda2*(qchg2/alpha)*
     x    (qchg/alpha))/volm)/dble(mxnode)
        
        qfix_fre=-(0.5d0*pi*r4pie0/epsq)*(((qchg2/alpha)**2-
     x    (qchg1/alpha)**2+2.d0*(qchg2/alpha)*(qchg/alpha)
     x    -2.d0*(qchg1/alpha)*(qchg/alpha))/volm)/dble(mxnode)
        
      else
        
        qfix=-(0.5d0*pi*r4pie0/epsq)*((qchg/alpha)**2/volm)/
     x    dble(mxnode)
        
      endif
      
c     add self interaction correction to potential
      
      if(lconsw)then
        
        eng1=engcpe
        engcpe=2.d0*rvolm*r4pie0*engcpe/epsq+engsic+qfix
        
        if(lfree)then
          
          cou_vir=2.d0*rvolm*r4pie0*(cou_vir-3.d0*cou_fre)/epsq-
     x      3.d0*qfix_fre
          cou_fre=2.d0*rvolm*r4pie0*cou_fre/epsq+qfix_fre+engsic2
          
        endif
        
        if(llsolva)then
          
          cou_sol(:)=2.d0*rvolm*r4pie0*cou_sol(:)/epsq+cou_sol_sic(:)
          if(lghost)cou_exc(:)=2.d0*rvolm*r4pie0*cou_exc(:)/epsq+
     x      cou_exc_sic(:)
          
        endif        
        
        scal1=2.d0*rvolm*r4pie0/epsq
        scale=4.d0*rvolm*r4pie0/epsq
        
      else
        
        eng1=engcpe
        engcpe=rvolm*r4pie0*engcpe/epsq+engsic+qfix
        
        if(lfree)then
          
          cou_vir=rvolm*r4pie0*(cou_vir-3.d0*cou_fre)/epsq-
     x      3.d0*qfix_fre
          cou_fre=rvolm*r4pie0*cou_fre/epsq+qfix_fre+engsic2
          
        endif
        
        if(llsolva)then
          
          cou_sol(:)=rvolm*r4pie0*cou_sol(:)/epsq+cou_sol_sic(:)
          if(lghost)cou_exc(:)=rvolm*r4pie0*cou_exc(:)/epsq+
     x      cou_exc_sic(:)
          
        endif        
        
        scal1=rvolm*r4pie0/epsq
        scale=2.d0*rvolm*r4pie0/epsq
        
      endif
      
      if(llsolva)then
        
        ksol=0
        do isol=1,mxtmls
          
          fkk=1.d0
          do jsol=1,isol
            
            ksol=ksol+1
            if(isol.eq.jsol)fkk=0.5d0
            
            cou_sol(ksol)=cou_sol(ksol)-
     x        ((fkk*pi*r4pie0/epsq)*qfix_sol(isol)*
     x        qfix_sol(jsol)/(alpha*alpha*volm*
     x        dble(mxnode)))
              
            if(lghost)then
              
              cou_exc(ksol)=cou_exc(ksol)-
     x          ((fkk*pi*r4pie0/epsq)*qfix_exc(isol)*
     x          qfix_exc(jsol)/(alpha*alpha*volm*
     x          dble(mxnode)))
              
            endif
            
          enddo
          
        enddo
        
      endif
      
c     calculate final forces
      
      do i=iatm0,iatm1
        
        fxx(i)=scale*fxx(i)
        fyy(i)=scale*fyy(i)
        fzz(i)=scale*fzz(i)
        
      enddo
      
c     calculate stress tensor (symmetrical)
      
      stress(1)=stress(1)+scal1*(omg(1)+eng1)+qfix
      stress(2)=stress(2)+scal1*omg(2)
      stress(3)=stress(3)+scal1*omg(3)
      stress(4)=stress(4)+scal1*omg(2)
      stress(5)=stress(5)+scal1*(omg(5)+eng1)+qfix
      stress(6)=stress(6)+scal1*omg(6)
      stress(7)=stress(7)+scal1*omg(3)
      stress(8)=stress(8)+scal1*omg(6)
      stress(9)=stress(9)+scal1*(omg(9)+eng1)+qfix
      
c     virial term
      
      vircpe=-scal1*(omg(1)+omg(5)+omg(9)+3.d0*eng1)-3.d0*qfix
      
      cou_fre=0.d0
      cou_vir=0.d0
      return
      end subroutine ewald1
      
      subroutine ewald2
     x  (lsolva,lfree,lghost,iatm,ik,engcpe,vircpe,drewd,rcut,epsq)
      
c***********************************************************************
c     
c     dl_poly subroutine for calculating coulombic forces in a
c     periodic system using ewald's method
c     
c     parallel replicated data version (part 2)
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith march 1992.
c     T3d optimised. t.forester july 1994
c     modified  - p.-a. cazade oct 2007: solvation, free energy etcc     
c     part 2 - real space terms. 
c     
c     Tabulated potential in r space
c     3pt interpolation
c     
c     t. forester March 1993
c     {stress tensor : t.forester june 1994}
c     
c***********************************************************************
      
      implicit none
      
      logical lsolva,lfree,lghost,lselect,lskip,idrive,jdrive
      integer m,ik,iatm,jatm,ll,l1,l2,kkk
      real(8) engcpe,vircpe,drewd,rcut,epsq
      real(8) chgprd,rsq,rrr,ppp,vk0,vk1,vk2,t1,t2,erfcr,egamma
      real(8) rcsq,rdrewd,chgea,fx,fy,fz,fi
      real(8) strs(6),strs_loc(6)
      
      dimension fi(3)
      
CDIR$ CACHE_ALIGN fi
      
      lskip=(lfree.or.lghost)
      if(lmetadyn)idrive=driven(ltype(iatm))
      
c     set cutoff condition for pair forces
      
      rcsq=rcut**2
      
c     reciprocal of interpolation interval
      
      rdrewd=1.d0/drewd
      
c     initialise stress tensor accumulators
      
      strs(:)=0.d0
      strs_loc(:)=0.d0
      
c     initialise potential energy and virial
      
      engcpe=0.d0
      vircpe=0.d0
      
c     start of primary loop for forces evaluation
      
      chgea=chge(iatm)/epsq*r4pie0
      
      if(abs(chgea).gt.1.d-10)then
        
        fi(1)=fxx(iatm)
        fi(2)=fyy(iatm)
        fi(3)=fzz(iatm)
        
        do m=1,ik
          
c     atomic index and charge product
          
          jatm=ilist(m)
          if(lmetadyn)jdrive=driven(ltype(jatm))          
          
          if(lskip)then
            if(atm_fre(iatm)*atm_fre(jatm).eq.2)cycle
          endif
          
          chgprd=chgea*chge(jatm)
          
c     Ignore interaction if product of charges is zero
          
          if(abs(chgprd).gt.1.d-10)then
            
c     calculate interatomic distance
            
            rsq=rsqdf(m)
            
c     apply truncation of potential
            
            if(rcsq.gt.rsq)then
              
              rrr=sqrt(rsq)              
              
              ll=int(rrr*rdrewd)
              l1=ll+1
              l2=ll+2
              ppp=rrr*rdrewd-dble(ll)
              
c     calculate interaction energy using 3-point interpolation
              
              vk0=erc(ll)
              vk1=erc(l1)
              vk2=erc(l2)
              t1=vk0+(vk1-vk0)*ppp
              t2=vk1+(vk2-vk1)*(ppp-1.0d0)
              erfcr=(t1+(t2-t1)*ppp*0.5d0)*chgprd
              
c     calculate forces using 3pt interpolation
              
              vk0=fer(ll)
              vk1=fer(l1)
              vk2=fer(l2)
              t1=vk0+(vk1-vk0)*ppp
              t2=vk1+(vk2-vk1)*(ppp-1.0d0)
              egamma=(t1+(t2-t1)*ppp*0.5d0)*chgprd
              
c     set selection control
              
              lselect=.true.
              
c     set double index
              
              if(lsolva)kkk=loc2(atmolt(iatm),atmolt(jatm))
              
              if(lghost)then
                
c     selected excitation option
                
                if((atm_fre(iatm).ne.1).and.(atm_fre(jatm).ne.1))then
                  
c     reset selection control
                  
                  lselect=(atm_fre(iatm)+atm_fre(jatm).eq.0)
                  
c     calculate solvation energy
                  
                  if(lsolva)cou_exc(kkk)=cou_exc(kkk)+erfcr
                  
                endif
                
              elseif(lfree)then
                
c     selected free energy option
                
                if((atm_fre(iatm).eq.1).or.(atm_fre(jatm).eq.1))then
                  
c     set hamiltonian mixing parameter
                  
                  cou_fre=cou_fre-erfcr
                  cou_vir=cou_vir+egamma*rsq
                  erfcr=lambda1*erfcr
                  egamma=lambda1*egamma
                  
                elseif((atm_fre(iatm).eq.2).or.(atm_fre(jatm).eq.2))then
                  
c     set hamiltonian mixing parameter
                  
                  cou_fre=cou_fre+erfcr
                  cou_vir=cou_vir-egamma*rsq
                  erfcr=lambda2*erfcr
                  egamma=lambda2*egamma
                  
                endif
                
              endif
              
              if(lselect)then
                
c     calculate potential energy and virial
                
                engcpe=engcpe+erfcr
                vircpe=vircpe-egamma*rsq
                
c     calculate solvation energy
                
                if(lsolva)cou_sol(kkk)=cou_sol(kkk)+erfcr
                
c     calculate forces
                
                fx=egamma*xdf(m)
                fy=egamma*ydf(m)
                fz=egamma*zdf(m)
                
                fi(1)=fi(1)+fx
                fi(2)=fi(2)+fy
                fi(3)=fi(3)+fz
                
                fxx(jatm)=fxx(jatm)-fx
                fyy(jatm)=fyy(jatm)-fy
                fzz(jatm)=fzz(jatm)-fz
                
c     calculate stress tensor
                
                strs(1)=strs(1)+xdf(m)*fx
                strs(2)=strs(2)+xdf(m)*fy
                strs(3)=strs(3)+xdf(m)*fz
                strs(4)=strs(4)+ydf(m)*fy
                strs(5)=strs(5)+ydf(m)*fz
                strs(6)=strs(6)+zdf(m)*fz
                
              endif
              
c     metadynamics local parameters
              
              if(lmetadyn.and.(idrive.or.jdrive))then
                
c     local energy and virial
                
                eng_loc=eng_loc+erfcr
                vir_loc=vir_loc-egamma*rsq
                
c     local forces
                
                fxx_loc(iatm)=fxx_loc(iatm)+fx
                fyy_loc(iatm)=fyy_loc(iatm)+fy
                fzz_loc(iatm)=fzz_loc(iatm)+fz
                
                fxx_loc(jatm)=fxx_loc(jatm)-fx
                fyy_loc(jatm)=fyy_loc(jatm)-fy
                fzz_loc(jatm)=fzz_loc(jatm)-fz
                
c     local stress tensor
                
                strs_loc(1)=strs_loc(1)+xdf(m)*fx
                strs_loc(2)=strs_loc(2)+xdf(m)*fy
                strs_loc(3)=strs_loc(3)+xdf(m)*fz
                strs_loc(4)=strs_loc(4)+ydf(m)*fy
                strs_loc(5)=strs_loc(5)+ydf(m)*fz
                strs_loc(6)=strs_loc(6)+zdf(m)*fz
                
              endif
              
            endif
            
          endif
          
        enddo
      
c     load temps back to fxx(iatm) etc
        
        fxx(iatm)=fi(1)
        fyy(iatm)=fi(2)
        fzz(iatm)=fi(3)
        
c     complete stress tensor
        
        stress(1)=stress(1)+strs(1)
        stress(2)=stress(2)+strs(2)
        stress(3)=stress(3)+strs(3)
        stress(4)=stress(4)+strs(2)
        stress(5)=stress(5)+strs(4)
        stress(6)=stress(6)+strs(5)
        stress(7)=stress(7)+strs(3)
        stress(8)=stress(8)+strs(5)
        stress(9)=stress(9)+strs(6)

        if(lmetadyn)then
          
          stress_loc(1)=stress_loc(1)+strs_loc(1)
          stress_loc(2)=stress_loc(2)+strs_loc(2)
          stress_loc(3)=stress_loc(3)+strs_loc(3)
          stress_loc(4)=stress_loc(4)+strs_loc(2)
          stress_loc(5)=stress_loc(5)+strs_loc(4)
          stress_loc(6)=stress_loc(6)+strs_loc(5)
          stress_loc(7)=stress_loc(7)+strs_loc(3)
          stress_loc(8)=stress_loc(8)+strs_loc(5)
          stress_loc(9)=stress_loc(9)+strs_loc(6)
          
        endif
        
      endif
      
      return
      end subroutine ewald2
      
      subroutine ewald3
     x  (lsolva,lfree,lghost,iatm,ilst,engcpe,vircpe,alpha,epsq)
      
c***********************************************************************
c     
c     dl_poly subroutine for calculating coulombic forces in a
c     periodic system using ewald's method
c     
c     parallel replicated data version (part 3)
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith dec 1992.
c     modified  - p.-a. cazade oct 2007: solvation, free energy etcc     
c     
c     stress stensor added t.forester may 1994
c     
c***********************************************************************
      implicit none
      
      logical lsolva,lfree,lghost,lselect,lskip,idrive,jdrive
      integer iatm,jatm,ilst,m,kkk
      real(8) engcpe,vircpe,alpha,epsq,a1,a2,a3
      real(8) a4,a5,pp,rr3,r10,r42,r216,chgea,chgprd,rrr,rsq,alpr
      real(8) alpr2,erfr,egamma,tt,exp1,fx,fy,fz
      real(8) strs(6),strs_loc(6)
      
      data a1,a2,a3/0.254829592d0,-0.284496736d0,1.421413741d0/
      data a4,a5,pp/-1.453152027d0,1.061405429d0,0.3275911d0/
      data rr3/0.333333333333d0/,r10/0.1d0/,r42/0.02380952381d0/
      data r216/4.62962962963d-3/
      
      lskip=(lfree.or.lghost)
      if(lmetadyn)idrive=driven(ltype(iatm))
      
c     initialise potential energy and virial
      
      engcpe=0.d0
      vircpe=0.d0
      
c     initialise stress accumulators
      
      strs(:)=0.d0
      strs_loc(:)=0.d0
      
c     start of primary loop for forces evaluation
      
      chgea=chge(iatm)/epsq*r4pie0
      
      do m=1,nexatm(ilst)
        
c     atomic index and charge product
        
        jatm=lexatm(ilst,m)
        if(lmetadyn)jdrive=driven(ltype(jatm))
        
        if(lskip)then
          if(atm_fre(iatm)*atm_fre(jatm).eq. 2)cycle
        endif
        
        chgprd=chgea*chge(jatm)
        
c     calculate interatomic distance
        
        rsq=xdf(m)**2+ydf(m)**2+zdf(m)**2
        
        rrr=sqrt(rsq)
        alpr=rrr*alpha
        alpr2=alpr*alpr
        
c     calculate error function and derivative
        
        if(alpr.lt.1.d-2)then
          
          erfr=2.d0*chgprd*(alpha/sqrpi)*
     x      (1.d0+alpr2*(-rr3+alpr2*(r10+alpr2*(-r42+alpr2*r216))))
          
          egamma=-4.d0*chgprd*(alpha**3/sqrpi)*
     x      (rr3+alpr2*(-2.d0*r10+alpr2*(3.d0*r42-4.d0*alpr2*r216)))
          
        else
          
          tt=1.d0/(1.d0+pp*alpha*rrr)
          exp1=exp(-(alpha*rrr)**2)
          erfr=(1.d0-tt*(a1+tt*(a2+tt*(a3+tt*(a4+tt*a5))))*exp1)*
     x      chgprd/rrr
          egamma=-(erfr-2.d0*chgprd*(alpha/sqrpi)*exp1)/rsq
          
        endif
        
c     set selection control
        
        lselect=.true.
        
c     set double index
        
        if(lsolva)kkk=loc2(atmolt(iatm),atmolt(jatm))
        
        if(lghost)then
          
c     selected excitation option
          
          if((atm_fre(iatm).ne.1).and.(atm_fre(jatm).ne.1))then
            
c     reset selection control
            
            lselect=(atm_fre(iatm)+atm_fre(jatm).eq.0)
            
c     calculate solvation energy
                  
            if(lsolva)cou_exc(kkk)=cou_exc(kkk)-erfr
            
          endif
          
        elseif(lfree)then
          
c     selected free energy option
          
          if((atm_fre(iatm).eq.1).or.(atm_fre(jatm).eq.1))then
            
c     set hamiltonian mixing parameter
            
            cou_fre=cou_fre+erfr
            cou_vir=cou_vir+egamma*rsq
            erfr=lambda1*erfr
            egamma=lambda1*egamma
            
          elseif((atm_fre(iatm).eq.2).or.(atm_fre(jatm).eq.2))then
            
c     set hamiltonian mixing parameter
            
            cou_fre=cou_fre-erfr
            cou_vir=cou_vir-egamma*rsq
            erfr=lambda2*erfr
            egamma=lambda2*egamma
            
          endif
          
        endif
        
        if(lselect)then
          
c     calculate potential energy and virial

          engcpe=engcpe-erfr
          vircpe=vircpe-egamma*rsq
        
c     calculate solvation energy
          
          if(lsolva)cou_sol(kkk)=cou_sol(kkk)-erfr
            
c     calculate forces
          
          fx=egamma*xdf(m)
          fy=egamma*ydf(m)
          fz=egamma*zdf(m)
          
          fxx(iatm)=fxx(iatm)+fx
          fyy(iatm)=fyy(iatm)+fy
          fzz(iatm)=fzz(iatm)+fz
          
          fxx(jatm)=fxx(jatm)-fx
          fyy(jatm)=fyy(jatm)-fy
          fzz(jatm)=fzz(jatm)-fz
          
c     calculate stress tensor
          
          strs(1)=strs(1)+xdf(m)*fx
          strs(2)=strs(2)+xdf(m)*fy
          strs(3)=strs(3)+xdf(m)*fz
          strs(4)=strs(4)+ydf(m)*fy
          strs(5)=strs(5)+ydf(m)*fz
          strs(6)=strs(6)+zdf(m)*fz
          
        endif
        
c     metadynamics local parameters
        
        if(lmetadyn.and.(idrive.or.jdrive))then
          
c     local energy and virial

          eng_loc=eng_loc-erfr
          vir_loc=vir_loc-egamma*rsq
        
c     local forces
          
          fxx_loc(iatm)=fxx_loc(iatm)+fx
          fyy_loc(iatm)=fyy_loc(iatm)+fy
          fzz_loc(iatm)=fzz_loc(iatm)+fz
          
          fxx_loc(jatm)=fxx_loc(jatm)-fx
          fyy_loc(jatm)=fyy_loc(jatm)-fy
          fzz_loc(jatm)=fzz_loc(jatm)-fz
          
c     local stress tensor
          
          strs_loc(1)=strs_loc(1)+xdf(m)*fx
          strs_loc(2)=strs_loc(2)+xdf(m)*fy
          strs_loc(3)=strs_loc(3)+xdf(m)*fz
          strs_loc(4)=strs_loc(4)+ydf(m)*fy
          strs_loc(5)=strs_loc(5)+ydf(m)*fz
          strs_loc(6)=strs_loc(6)+zdf(m)*fz
          
        endif
        
      enddo
      
c     complete stress tensor
      
      stress(1)=stress(1)+strs(1)
      stress(2)=stress(2)+strs(2)
      stress(3)=stress(3)+strs(3)
      stress(4)=stress(4)+strs(2)
      stress(5)=stress(5)+strs(4)
      stress(6)=stress(6)+strs(5)
      stress(7)=stress(7)+strs(3)
      stress(8)=stress(8)+strs(5)
      stress(9)=stress(9)+strs(6)

      if(lmetadyn)then
        
        stress_loc(1)=stress_loc(1)+strs_loc(1)
        stress_loc(2)=stress_loc(2)+strs_loc(2)
        stress_loc(3)=stress_loc(3)+strs_loc(3)
        stress_loc(4)=stress_loc(4)+strs_loc(2)
        stress_loc(5)=stress_loc(5)+strs_loc(4)
        stress_loc(6)=stress_loc(6)+strs_loc(5)
        stress_loc(7)=stress_loc(7)+strs_loc(3)
        stress_loc(8)=stress_loc(8)+strs_loc(5)
        stress_loc(9)=stress_loc(9)+strs_loc(6)
        
      endif
      
      return
      end subroutine ewald3
      
      subroutine ewald4(lsolva,lfree,lghost,iatm,ik,engcpe,vircpe,
     x  engcpl,vircpl,drewd,rcut,epsq)
      
c***********************************************************************
c     
c     dl_poly subroutine for calculating coulombic forces in a
c     periodic system using ewald's method
c     
c     modified to allow direct calculation of primary (short-range)
c     interactions for multiple-time step corrections
c      
c     primary neighbours are taken out of the Ewald sum
c     electrostatics are evaluated directly instead
c     
c     parallel replicated data version (part 2)
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith march 1992.
c     
c     modified  - t. forester february 1993
c     
c     part 2 - real space terms
c     
c***********************************************************************
      
      implicit none
      
      logical lsolva,lfree,lghost,lselect,lskip
      integer iatm,ik,m,jatm,ll,i,kkk
      real(8) engcpe,vircpe,engcpl,vircpl,drewd,rcut,epsq,rrr
      real(8) rcsq,rdrewd,strs,strl,chgea,fi,fli,rsq,chgprd,coul
      real(8) vk0,vk1,vk2,t1,t2,erfcr,fcoul,egamma,fx,fy,fz,ppp
      
      dimension fi(3),fli(3),strs(6),strl(6)
      
CDIR$ CACHE_ALIGN fi
CDIR$ CACHE_ALIGN fli
      
      lskip=(lfree.or.lghost)
      
c     set cutoff condition for pair forces
      
      rcsq=rcut**2
      
c     reciprocal of interpolation interval
      
      rdrewd=1.d0/drewd
      
c     initialise potential energy and virial
      
      engcpe=0.d0
      vircpe=0.d0
      
      engcpl=0.d0
      vircpl=0.d0
      
c     initialise stress tensor accumulators
      
      do i=1,6
        
        strs(i)=0.d0
        strl(i)=0.d0
        
      enddo
      
c     start of primary loop for forces evaluation
      
      chgea=chge(iatm)/epsq*r4pie0
      if(abs(chgea).gt.1.d-10)then
        
c     temporary arrays for cache aligning
        
        fi(1)=fxx(iatm)
        fi(2)=fyy(iatm)
        fi(3)=fzz(iatm)
        
        fli(1)=flx(iatm)
        fli(2)=fly(iatm)
        fli(3)=flz(iatm)
        
        do m=1,ik
          
c     atomic index and charge product
          
          jatm=ilist(m)
          
          if(lskip)then
            if(atm_fre(iatm)*atm_fre(jatm).eq.2)cycle
          endif
          
          chgprd=chgea*chge(jatm)
          
          if(abs(chgprd).gt.1.d-10)then
            
c     calculate interatomic distance
            
            rsq=rsqdf(m)
            
            if(rcsq.gt.rsq)then
              
c     coulombic energy and coulombic force
              
              rrr=sqrt(rsq)
              coul=chgprd/rrr
              fcoul=coul/rsq
              
c     calculate Ewald term using 3-point interpolation
              
              ll=int(rrr*rdrewd)
              ppp=rrr*rdrewd-dble(ll)
              
              vk0=erc(ll)
              vk1=erc(ll+1)
              vk2=erc(ll+2)
              t1=vk0+(vk1-vk0)*ppp
              t2=vk1+(vk2-vk1)*(ppp-1.0d0)
              erfcr=(t1+(t2-t1)*ppp*0.5d0)*chgprd              
              
c     calculate Ewald forces using 3pt interpolation
              
              vk0=fer(ll)
              vk1=fer(ll+1)
              vk2=fer(ll+2)
              t1=vk0+(vk1-vk0)*ppp
              t2=vk1+(vk2-vk1)*(ppp-1.0d0)
              egamma=(t1+(t2-t1)*ppp*0.5d0)*chgprd
              egamma=egamma-fcoul
              
c     set selection control
              
              lselect=.true.
              
c     set double index
              
              if(lsolva)kkk=loc2(atmolt(iatm),atmolt(jatm))
              
              if(lghost)then
                
c     selected excitation option
                
                if((atm_fre(iatm).ne.1).and.(atm_fre(jatm).ne.1))then
                  
c     reset selection control
                  
                  lselect=(atm_fre(iatm)+atm_fre(jatm).eq.0)
                  
c     calculate solvation energy
                  
                  if(lsolva)then
                    cou_exc(kkk)=cou_exc(kkk)+coul
                    cou_exc_lng(kkk)=cou_exc_lng(kkk)+erfcr-coul
                  endif
                  
                endif
                
              elseif(lfree)then
                
c     selected free energy option
                
                if((atm_fre(iatm).eq.1).or.(atm_fre(jatm).eq.1))then
                  
c     set hamiltonian mixing parameter
                  
                  cou_fre=cou_fre-coul
                  cou_vir=cou_vir+fcoul*rsq
                  coul=lambda1*coul
                  erfcr=lambda1*erfcr
                  egamma=lambda1*egamma
                  
                elseif((atm_fre(iatm).eq.2).or.(atm_fre(jatm).eq.2))then
                  
c     set hamiltonian mixing parameter

                  cou_fre=cou_fre+coul
                  cou_vir=cou_vir-fcoul*rsq
                  coul=lambda2*coul
                  erfcr=lambda2*erfcr
                  egamma=lambda2*egamma
                  
                endif
                
              endif
              
              if(lselect)then
                
c     calculate potential energy and virial
                
                engcpe=engcpe+coul
                engcpl=engcpl+erfcr-coul
                vircpe=vircpe-fcoul*rsq
                vircpl=vircpl-egamma*rsq
                
c     calculate solvation energy
                
                if(lsolva)then
                  cou_sol(kkk)=cou_sol(kkk)+coul
                  cou_sol_lng(kkk)=cou_sol_lng(kkk)+erfcr-coul
                endif
                
c     add in contributions to the long-range force
                
                fx=egamma*xdf(m)
                fy=egamma*ydf(m)
                fz=egamma*zdf(m)
                
                fli(1)=fli(1)+fx
                fli(2)=fli(2)+fy
                fli(3)=fli(3)+fz
                
                flx(jatm)=flx(jatm)-fx
                fly(jatm)=fly(jatm)-fy
                flz(jatm)=flz(jatm)-fz
                
c     calculate stress tensor
                
                strl(1)=strl(1)+xdf(m)*fx
                strl(2)=strl(2)+xdf(m)*fy
                strl(3)=strl(3)+xdf(m)*fz
                strl(4)=strl(4)+ydf(m)*fy
                strl(5)=strl(5)+ydf(m)*fz
                strl(6)=strl(6)+zdf(m)*fz
              
c     add in contributions to instantaneous force
                
                fx=fcoul*xdf(m)
                fy=fcoul*ydf(m)
                fz=fcoul*zdf(m)
                
                fi(1)=fi(1)+fx
                fi(2)=fi(2)+fy
                fi(3)=fi(3)+fz
                
                fxx(jatm)=fxx(jatm)-fx
                fyy(jatm)=fyy(jatm)-fy
                fzz(jatm)=fzz(jatm)-fz
                
c     calculate stress tensor
                
                strs(1)=strs(1)+xdf(m)*fx
                strs(2)=strs(2)+xdf(m)*fy
                strs(3)=strs(3)+xdf(m)*fz
                strs(4)=strs(4)+ydf(m)*fy
                strs(5)=strs(5)+ydf(m)*fz
                strs(6)=strs(6)+zdf(m)*fz
                
              endif
              
            endif
            
          endif
          
        enddo
        
c     copy back temporaries
        
        fxx(iatm)=fi(1)
        fyy(iatm)=fi(2)
        fzz(iatm)=fi(3)
        
        flx(iatm)=fli(1)
        fly(iatm)=fli(2)
        flz(iatm)=fli(3)
        
c     complete stress tensorS
        
        stresl(1)=stresl(1)+strl(1)
        stresl(2)=stresl(2)+strl(2)
        stresl(3)=stresl(3)+strl(3)
        stresl(4)=stresl(4)+strl(2)
        stresl(5)=stresl(5)+strl(4)
        stresl(6)=stresl(6)+strl(5)
        stresl(7)=stresl(7)+strl(3)
        stresl(8)=stresl(8)+strl(5)
        stresl(9)=stresl(9)+strl(6)
        
        stress(1)=stress(1)+strs(1)
        stress(2)=stress(2)+strs(2)
        stress(3)=stress(3)+strs(3)
        stress(4)=stress(4)+strs(2)
        stress(5)=stress(5)+strs(4)
        stress(6)=stress(6)+strs(5)
        stress(7)=stress(7)+strs(3)
        stress(8)=stress(8)+strs(5)
        stress(9)=stress(9)+strs(6)
        
      endif
      
      return
      end subroutine ewald4
      
      end module ewald_module
