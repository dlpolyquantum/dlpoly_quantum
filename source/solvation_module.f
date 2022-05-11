      module solvation_module
      
c***********************************************************************
c     
c     dl_poly module for defining decomposition of energy arrays
c     to calculate solvation energies
c     copyright - daresbury laboratory
c     author    - p.-a. cazade oct 2007
c     adapted   - w. smith     aug 2008
c     
c***********************************************************************
      
      use setup_module
      use error_module
      use config_module
      
      implicit none
      
      logical lcomp(9)
      
      integer mxtmls_sol2,mxtmls_sol3,mxtmls_sol4
      integer ntcons_ghost,ngrp_ghost,nswitch,niswitch
      integer mxtmls_exc2,mxtmls_exc3,mxtmls_exc4
      integer nfrn,ifrn,mfree,kfree,ind_fre(4)
      
      real(8) pfree,lambda1,lambda2,dlambda,eng_kin_fre
      real(8) elrc2,virlrc2,engsic0,engsic2,elrc_sav,qfix_fre
      real(8) virlrc_sav,volm_sav,elrc_fre,elrc_fre_sav,vlrc_fre
      real(8) qchg0,qchg1,qchg2,ang_fre,bnd_fre,dih_fre
      real(8) inv_fre,tbp_fre,fbp_fre,cou_fre,elrc2_sav
      real(8) vdw_fre,eng_cfg_fre,shl_fre,virlrc2_sav,vlrc_fre_sav
      real(8) qchg_sav,vir_cfg_fre,ang_vir,bnd_vir,dih_vir
      real(8) inv_vir,tbp_vir,fbp_vir,cou_vir,vdw_vir,shl_vir
      
      integer, allocatable :: atm_fre(:)
      integer, allocatable :: atmolt(:),rigid_sol(:),const_sol(:)
      integer, allocatable :: lstgot_sol(:),natm_sol(:)
      
      real(8), allocatable :: elrc_sol(:),elrc_sol_sav(:),shl_sol(:)
      real(8), allocatable :: cou_sol(:),vdw_sol(:),bnd_sol(:)
      real(8), allocatable :: ckc_sol_sum(:),cks_sol_sum(:)
      real(8), allocatable :: cou_sol_sic(:),ebuf_sol1(:),inv_sol(:)
      real(8), allocatable :: ang_sol(:),dih_sol(:),en3_sol(:)
      real(8), allocatable :: qfix_sol(:),ebuf_sol2(:),en4_sol(:)
      real(8), allocatable :: vdw_sol_lng(:),cou_sol_lng(:)
      real(8), allocatable :: degfre_sol(:),degrot_sol(:),temp_sol(:)
      real(8), allocatable :: vxo_sol(:),vyo_sol(:),vzo_sol(:)
      real(8), allocatable :: ckc1(:),cks1(:),ckc2(:),cks2(:)
      real(8), allocatable :: ckc_fre_sum(:),cks_fre_sum(:)
      real(8), allocatable :: ebuf_exc1(:),ebuf_exc2(:)
      real(8), allocatable :: vxo_fre(:),vyo_fre(:),vzo_fre(:)
      real(8), allocatable :: elrc_exc(:),elrc_exc_sav(:)
      real(8), allocatable :: cou_exc(:),vdw_exc(:),bnd_exc(:)
      real(8), allocatable :: ang_exc(:),dih_exc(:),en4_exc(:)
      real(8), allocatable :: vdw_exc_lng(:),cou_exc_lng(:)
      real(8), allocatable :: shl_exc(:),en3_exc(:),inv_exc(:)
      real(8), allocatable :: qfix_exc(:),cou_exc_sic(:),weight_sav(:)
      
      save atmolt,rigid_sol,const_sol,lstgot_sol,natm_sol,elrc_sol
      save elrc_sol_sav,cou_sol,vdw_sol,bnd_sol,ckc_sol_sum,cks_sol_sum
      save cou_sol_sic,ebuf_sol1,en4_sol,ang_sol,dih_sol,en3_sol
      save qfix_sol,ebuf_sol2,shl_sol,vdw_sol_lng,cou_sol_lng,degfre_sol
      save degrot_sol,inv_sol,temp_sol,vxo_sol,vyo_sol,vzo_sol
      save mxtmls_sol2,mxtmls_sol3,mxtmls_sol4,lcomp
      
      save nfrn,ifrn,mfree,kfree,pfree,lambda1,lambda2,dlambda
      save eng_kin_fre,elrc2,virlrc2,engsic0,engsic2,cks_fre_sum
      save virlrc_sav,volm_sav,elrc_fre,elrc_fre_sav,qchg1,qchg2
      save ang_fre,bnd_fre,dih_fre,inv_fre,tbp_fre,fbp_fre,qchg0
      save cou_fre,vdw_fre,eng_cfg_fre,elrc2_sav,elrc_sav,vlrc_fre
      save ang_vir,bnd_vir,dih_vir,inv_vir,tbp_vir,fbp_vir,cou_vir
      save vdw_vir,shl_vir,vir_cfg_fre,qfix_fre,virlrc2_sav
      save ind_fre,atm_fre,ckc1,cks1,ckc2,cks2,ckc_fre_sum
      save ebuf_exc1,ebuf_exc2,vxo_fre,vyo_fre,vzo_fre,vlrc_fre_sav
      save weight_sav
      
      save ntcons_ghost,ngrp_ghost,qchg_sav,nswitch,niswitch
      save mxtmls_exc2,mxtmls_exc3,mxtmls_exc4
      save cou_exc,vdw_exc,bnd_exc,ang_exc,dih_exc,en4_exc,vdw_exc_lng
      save cou_exc_lng,shl_exc,en3_exc,inv_exc,elrc_exc,elrc_exc_sav
      save qfix_exc,cou_exc_sic
      
      contains
      
      subroutine alloc_sol_arrays(idnode,mxnode)
      
c***********************************************************************
c     
c     dl_poly routine for allocating solvation module arrays
c     copyright - daresbury laboratory
c     author    - p.-a. cazade oct 2007
c     adapted   - w. smith     aug 2008
c     
c***********************************************************************
      
      implicit none
      
      integer, parameter :: nnn=37
      
      logical safe
      integer i,fail,idnode,mxnode
      dimension fail(nnn)
      
      safe=.true.
      
      mxtmls_sol2=((mxtmls_sol+1)*mxtmls_sol)/2
      mxtmls_sol3=(((mxtmls_sol+3)*mxtmls_sol+2)*mxtmls_sol)/6
      mxtmls_sol4=((((mxtmls_sol+6)*mxtmls_sol+11)*mxtmls_sol+6)*
     x  mxtmls_sol)/24
      
c     allocate arrays
      
      fail(:)=0
      
      allocate (cou_sol(mxtmls_sol2),stat=fail(1))
      allocate (vdw_sol(mxtmls_sol2),stat=fail(2))
      allocate (shl_sol(mxtmls_sol),stat=fail(3))
      allocate (ebuf_sol1(mxebuf_sol),stat=fail(4))
      allocate (cou_sol_sic(mxtmls_sol2),stat=fail(5))
      allocate (cks_sol_sum(mxtmls_sol),stat=fail(6))
      allocate (ckc_sol_sum(mxtmls_sol),stat=fail(7))
      allocate (bnd_sol(mxtmls_sol),stat=fail(8))
      allocate (ang_sol(mxtmls_sol),stat=fail(9))
      allocate (dih_sol(mxtmls_sol),stat=fail(10))
      allocate (atmolt(mxatms_sol),stat=fail(11))
      allocate (en3_sol(mxtmls_sol3),stat=fail(12))
      allocate (en4_sol(mxtmls_sol4),stat=fail(13))
      allocate (qfix_sol(mxtmls_sol),stat=fail(14))
      allocate (elrc_sol(mxtmls_sol2),stat=fail(15))
      allocate (elrc_sol_sav(mxtmls_sol2),stat=fail(16))
      allocate (ebuf_sol2(mxebuf_sol),stat=fail(23))
      allocate (rigid_sol(mxtmls_sol),stat=fail(24))
      allocate (const_sol(mxtmls_sol),stat=fail(25))
      allocate (degfre_sol(mxtmls_sol),stat=fail(26))
      allocate (degrot_sol(mxtmls_sol),stat=fail(27))
      allocate (natm_sol(mxtmls_sol),stat=fail(28))
      allocate (lstgot_sol(mxatms_sol),stat=fail(29))
      allocate (temp_sol(mxtmls_sol),stat=fail(30))
      allocate (vxo_sol(mxatms_sol),stat=fail(31))
      allocate (vyo_sol(mxatms_sol),stat=fail(32))
      allocate (vzo_sol(mxatms_sol),stat=fail(33))
      allocate (vdw_sol_lng(mxtmls_sol2),stat=fail(34))
      allocate (cou_sol_lng(mxtmls_sol2),stat=fail(35))
      allocate (inv_sol(mxtmls_sol),stat=fail(36))
      allocate (weight_sav(mxatms_fre),stat=fail(37))
      
      if(any(fail.gt.0))safe=.false.      
      if(mxnode.gt.1)call gstate(safe)    
      if(.not.safe)call error(idnode,540)
      
c     initialise accumulators
      
      lcomp(:)=.false.
      elrc_sol(:)=0.d0
      cou_sol(:)=0.d0
      vdw_sol(:)=0.d0
      en3_sol(:)=0.d0
      en4_sol(:)=0.d0
      bnd_sol(:)=0.d0
      ang_sol(:)=0.d0
      dih_sol(:)=0.d0
      inv_sol(:)=0.d0
      
      return
      end subroutine alloc_sol_arrays
      
      subroutine solva_temp(idnode,mxnode,natms,keyver)
      
c***********************************************************************
c     
c     dl_poly routine for solvation module
c     copyright - daresbury laboratory
c     author    - p.-a. cazade oct 2007
c     adapted   - w. smith     aug 2008
c     
c***********************************************************************
      
      implicit none
      
      integer i,ii,idnode,mxnode,natms,keyver
      real(8) vvx,vvy,vvz
      
      temp_sol(:)=0.d0
      
      do i=idnode+1,natms,mxnode
        
        if(keyver.eq.0)then
          
          vvx=0.5d0*(vxx(i)+vxo_sol(i))
          vvy=0.5d0*(vyy(i)+vyo_sol(i))
          vvz=0.5d0*(vzz(i)+vzo_sol(i))
          
        else
          
          vvx=vxx(i)
          vvy=vyy(i)
          vvz=vzz(i)
          
        endif
        
        if(degfre_sol(atmolt(i)).ge.1.d0)then
          
          temp_sol(atmolt(i))=temp_sol(atmolt(i))+weight(i)*
     x      (vvx*vvx+vvy*vvy+vvz*vvz)/(boltz*degfre_sol(atmolt(i)))
          
        endif
        
      enddo
      
c     global sum
      
      if(mxnode.gt.1)call gdsum(temp_sol,mxtmls_sol,buffer)
        
      return
      end subroutine solva_temp
      
      subroutine alloc_free_arrays(idnode,mxnode)
      
c***********************************************************************
c     
c     dl_poly routine to allocate free energy arrays
c     copyright - daresbury laboratory
c     author    - p.-a. cazade oct 2007
c     adapted   - w.smith      aug 2008
c     
c***********************************************************************
      
      implicit none
      
      integer, parameter :: nnn=12
      logical safe
      integer i,fail,idnode,mxnode
      dimension fail(nnn)

      safe=.true.

c     allocate arrays

      fail(:)=0
      
      allocate (ebuf_exc1(mxebuf_fre),stat=fail(1))
      allocate (ckc1(mxewld_fre),stat=fail(2))
      allocate (cks1(mxewld_fre),stat=fail(3))
      allocate (ckc2(mxewld_fre),stat=fail(4))
      allocate (cks2(mxewld_fre),stat=fail(5))
      allocate (cks_fre_sum(mxtmls_fre),stat=fail(6))
      allocate (ckc_fre_sum(mxtmls_fre),stat=fail(7))
      allocate (atm_fre(mxatms_fre),stat=fail(8))
      allocate (ebuf_exc2(mxebuf_fre),stat=fail(9))
      allocate (vxo_fre(mxatms_fre),stat=fail(10))
      allocate (vyo_fre(mxatms_fre),stat=fail(11))
      allocate (vzo_fre(mxatms_fre),stat=fail(12))
      
      if(any(fail.gt.0))safe=.false.      
      if(mxnode.gt.1)call gstate(safe)    
      if(.not.safe)call error(idnode,541)
      
      return
      
      end subroutine alloc_free_arrays
      
      subroutine lrcorrect_fre(lfree,volm,elrc,virlrc)
      
c***********************************************************************
c     
c     dl_poly routine for free energy module
c     copyright - daresbury laboratory
c     author    - p.-a. cazade oct 2007
c     adapted   - w.smith      aug 2008
c     
c***********************************************************************
      
      logical lfree
      real(8) volm,elrc,virlrc
      
      elrc=elrc_sav*(volm_sav/volm)
      elrc2=elrc2_sav*(volm_sav/volm)
      virlrc=virlrc_sav*(volm_sav/volm)
      virlrc2=virlrc2_sav*(volm_sav/volm)
      if(lfree)then
        elrc_fre=elrc_fre_sav*(volm_sav/volm)
        vlrc_fre=vlrc_fre_sav*(volm_sav/volm)
      endif
      
      end subroutine lrcorrect_fre
      
      subroutine free_kinetic(lfrmas,idnode,mxnode,keyver)
      
c***********************************************************************
c     
c     dl_poly routine for free energy module
c     calculate kinetic energy difference between states
c     
c     copyright - daresbury laboratory
c     author    - p.-a. cazade oct 2007
c     adapted   - w. smith     aug 2008: parallel version
c     
c***********************************************************************
      
      implicit none
      
      logical lfrmas
      integer i,idnode,mxnode,keyver
      real(8) fac
      
      fac=dlambda
      eng_kin_fre=0.d0
      
      if(keyver.eq.0)then
        
        if(lfrmas)fac=dlambda/lambda1
        do i=ind_fre(1)+idnode,ind_fre(2),mxnode
          
          eng_kin_fre=eng_kin_fre-
     x      fac*weight(i)*((vxx(i)+vxo_fre(i))**2+
     x      (vyy(i)+vyo_fre(i))**2+(vzz(i)+vzo_fre(i))**2)
          
        enddo
        
        if(lfrmas)fac=dlambda/lambda2
        do i=ind_fre(3)+idnode,ind_fre(4),mxnode
          
          eng_kin_fre=eng_kin_fre+
     x      fac*weight(i)*((vxx(i)+vxo_fre(i))**2+
     x      (vyy(i)+vyo_fre(i))**2+(vzz(i)+vzo_fre(i))**2)
          
        enddo
        
        eng_kin_fre=eng_kin_fre/8.d0

      else
        
        if(lfrmas)fac=dlambda/lambda1
        do i=ind_fre(1)+idnode,ind_fre(2),mxnode
          
          eng_kin_fre=eng_kin_fre-fac*weight(i)*(vxx(i)**2+
     x      vyy(i)**2+vzz(i)**2)
          
        enddo
        
        if(lfrmas)fac=dlambda/lambda2
        do i=ind_fre(3)+idnode,ind_fre(4),mxnode
          
          eng_kin_fre=eng_kin_fre+fac*weight(i)*(vxx(i)**2+
     x      vyy(i)**2+vzz(i)**2)
          
        enddo
        
        eng_kin_fre=eng_kin_fre/2.d0
        
      endif
      
c     global sum
      
      if(mxnode.gt.1)then
        
        buffer(1)=eng_kin_fre
        call gdsum(buffer(1),1,buffer(2))
        eng_kin_fre=buffer(1)
        
      endif
      
      return
      end subroutine free_kinetic
      
      subroutine freegen()
      
c***********************************************************************
c     
c     dl_poly routine for free energy module: select mixing scheme
c     copyright - daresbury laboratory
c     author    - p.-a. cazade oct 2007
c     adapted   - w.smith      aug 2008
c     
c***********************************************************************
      
      implicit none
      
      integer i,j,fac1,fac2
      real(8) sigma1,sigma2,acc,arg,gss,tt,pp,a1,a2,a3,a4,a5,err
      
      data a1,a2,a3/0.254829592d0,-0.284496736d0,1.421413741d0/
      data a4,a5,pp/-1.453152027d0,1.061405429d0,0.3275911d0/

      if(mfree.eq.1)then
        
c     linear mixing
        
        lambda1=(1.d0-pfree)
        lambda2=pfree
        dlambda=1.d0
        
      elseif(mfree.eq.2)then
        
c     nonlinear mixing
        
        lambda1=(1.d0-pfree)**kfree
        lambda2=(1.d0-(1.d0-pfree)**kfree)
        dlambda=dble(kfree)*(1.d0-pfree)**(kfree-1)
        
      elseif(mfree.eq.3)then
        
c     trigonmetric mixing
        
        lambda2=0.5d0*(1.d0+sin(pi*(pfree-0.5d0)))
        lambda1=1.d0-lambda2
        dlambda=0.5d0*pi*cos(pi*(pfree-0.5d0))
        
      elseif(mfree.eq.4)then
        
c     error function mixing
        
        acc=12.d0
        arg=2.d0*sqrt(2.302585093*acc)
        gss=exp(-(arg*(pfree-0.5d0))**2)
        tt=1.d0/(1.d0+pp*arg*abs(pfree-0.5d0))
        err=1.d0-tt*(a1+tt*(a2+tt*(a3+tt*(a4+tt*a5))))*gss
        lambda2=0.5d0*(1.d0+sign(err,(pfree-0.5d0)))
        lambda1=1.d0-lambda2
        dlambda=gss*arg/sqrpi
        
      elseif(mfree.eq.5)then
        
c     polynomial mixing
        
        sigma1=0.d0
        
        do i=0,kfree-1
          
          fac1=1
          fac2=1
          
          do j=0,i-1
            
            fac1=fac1*(kfree-1+i-j)
            fac2=fac2*(i-j)
            
          enddo
          
          sigma1=sigma1+(dble(fac1/fac2))*pfree**i
          
        enddo
        
        lambda1=sigma1*(1.d0-pfree)**kfree
        lambda2=1.d0-lambda1
        dlambda=sigma1*kfree*(1.d0-pfree)**(kfree-1)
        
        sigma2=0.d0
        
        do i=1,kfree-1
          
          fac1=1
          fac2=1
          
          do j=0,i-1
            
            fac1=fac1*(kfree-1+i-j)
            fac2=fac2*(i-j)
            
          enddo
          
          sigma2=sigma2+(dble(fac1*i/fac2))*pfree**(i-1)
          
        enddo
        
        dlambda=dlambda-sigma2*(1.d0-pfree)**kfree
        
      else
        
c     spline kernel mixing
        
        arg=pfree-0.5d0
        lambda2=2.d0*pfree-8.d0*arg**3*(1.d0-abs(arg))-0.5d0
        lambda1=1.d0-lambda2
        dlambda=2.d0+arg**2*(32.d0*abs(arg)-24.d0)
        
      endif
      
      return
      end subroutine freegen
      
      subroutine free_energy_write(idnode,nstep,engunit)
      
c***********************************************************************
c     
c     dl_poly subroutine for writing free energy file at selected
c     intervals in simulation
c     
c     copyright - daresbury  laboratory
c     author    - p.-a. cazade dec 2007
c     adapted   - w. smith     aug 2008
c     
c***********************************************************************
      
      implicit none
      
      logical newjob
      integer idnode,natms,nstep
      real(8) engunit
      
      save newjob
      data newjob/.true./
      
      if(idnode.eq.0)then
          
c     open the FREENG file if new job or file closed
        
        if(newjob)then
          
          newjob = .false.
          open(nfrnwr,file='FREENG',position='append')
          
        endif
        
        if(nstep.eq.nfrn.or.nstep.eq.ifrn)then
          
          write(nfrnwr,'(80a1)')cfgname
          
          if(abs(engunit-9648.530821d0).le.1.d-10) write(nfrnwr,
     x      "(' ENERGY UNITS=electron Volts ')")
          if(abs(engunit-418.4d0).le.1.d-10)       write(nfrnwr,
     x      "(' ENERGY UNITS=kcal/mol ')")
          if(abs(engunit-1.d2).le.1.d-10)          write(nfrnwr,
     x      "(' ENERGY UNITS=kjoule/mol ')")
          if(abs(engunit-boltz).lt.1.d-10)         write(nfrnwr,
     x      "(' ENERGY UNITS=kelvin ')")
          if(abs(engunit-1.d0).lt.1.d-10)          write(nfrnwr,
     x      "(' ENERGY UNITS=DL_POLY Internal UNITS ')")
          
          write(nfrnwr,'(1p,4e16.8)')pfree,lambda1,lambda2,dlambda
          
        endif
        
        if(mod(nstep-nfrn,ifrn).eq.0)then
          
          write(nfrnwr,"(i10,1p,2e16.8)")
     x      nstep,eng_cfg_fre/engunit,vir_cfg_fre/engunit
          
        endif
        
      endif
      
      return
      end subroutine free_energy_write
      
      subroutine alloc_exi_arrays(idnode,mxnode)
      
c***********************************************************************
c     
c     dl_poly routine to allocate excited state arrays
c     copyright - daresbury laboratory
c     author    - p.-a. cazade oct 2007
c     adapted   - w. smith     aug 2008
c     
c***********************************************************************
      
      implicit none
      
      logical safe
      integer, parameter :: nnn=15
      integer i,fail,idnode,mxnode
      dimension fail(nnn)

      safe=.true.
      mxtmls_exc2=((mxtmls_exc+1)*mxtmls_exc)/2
      mxtmls_exc3=(((mxtmls_exc+3)*mxtmls_exc+2)*mxtmls_exc)/6
      mxtmls_exc4=((((mxtmls_exc+6)*mxtmls_exc+11)*mxtmls_exc+6)*
     x  mxtmls_exc)/24

c     allocate arrays
      
      fail(:)=0
      
      allocate (cou_exc(mxtmls_exc2),stat=fail(1))
      allocate (vdw_exc(mxtmls_exc2),stat=fail(2))
      allocate (bnd_exc(mxtmls_exc),stat=fail(3))
      allocate (ang_exc(mxtmls_exc),stat=fail(4))
      allocate (dih_exc(mxtmls_exc),stat=fail(5))
      allocate (en3_exc(mxtmls_exc3),stat=fail(6))
      allocate (en4_exc(mxtmls_exc4),stat=fail(7))
      allocate (shl_exc(mxtmls_exc),stat=fail(8))
      allocate (vdw_exc_lng(mxtmls_exc2),stat=fail(9))
      allocate (cou_exc_lng(mxtmls_exc2),stat=fail(10))
      allocate (inv_exc(mxtmls_exc),stat=fail(11))
      allocate (elrc_exc(mxtmls_exc2),stat=fail(12))
      allocate (elrc_exc_sav(mxtmls_exc2),stat=fail(13))
      allocate (qfix_exc(mxtmls_exc),stat=fail(14))
      allocate (cou_exc_sic(mxtmls_exc2),stat=fail(15))
      
      if(any(fail.gt.0))safe=.false.      
      if(mxnode.gt.1)call gstate(safe)    
      if(.not.safe)call error(idnode,542)
      
c     initialise accumulators
      
      bnd_exc(:)=0.d0
      ang_exc(:)=0.d0
      dih_exc(:)=0.d0
      inv_exc(:)=0.d0
      en3_exc(:)=0.d0
      en4_exc(:)=0.d0
      elrc_exc(:)=0.d0
      cou_exc(:)=0.d0
      vdw_exc(:)=0.d0
      
      return
      end subroutine alloc_exi_arrays
      
      subroutine update_ghost(idnode,mxnode)
      
c***********************************************************************
c     
c     dl_poly routine for excited state module
c     update the positions of ghost atoms
c     copyright - daresbury laboratory
c     author    - p.-a. cazade oct 2007
c     adapted   - w. smith     aug 2008
c     
c***********************************************************************
      
      implicit none
      
      integer i,ii,idnode,mxnode
      
      if(ind_fre(1).lt.ind_fre(3))then
        
        ii=ind_fre(3)
        
        do i=ind_fre(1),ind_fre(2)
          
          xxx(ii)=xxx(i)
          yyy(ii)=yyy(i)
          zzz(ii)=zzz(i)
          
          ii=ii+1
          
        enddo
        
      else
        
        ii=ind_fre(1)
        
        do i=ind_fre(3),ind_fre(4)
          
          xxx(ii)=xxx(i)
          yyy(ii)=yyy(i)
          zzz(ii)=zzz(i)
          
          ii=ii+1
          
        enddo
        
      endif
      
      return
      end subroutine update_ghost
      
      subroutine copy_force(idnode,mxnode)
      
c***********************************************************************
c     
c     dl_poly routine for excited state module
c     copyright - daresbury laboratory
c     author    - p.-a. cazade oct 2007
c     adapted   - w. smith     aug 2008
c     
c***********************************************************************
      
      implicit none
      
      integer i,ii,idnode,mxnode
      
      ii=ind_fre(3)
      
      do i=ind_fre(1),ind_fre(2)
        
        fxx(ii)=fxx(i)
        fyy(ii)=fyy(i)
        fzz(ii)=fzz(i)
        
        ii=ii+1
        
      enddo
      
      return
      end subroutine copy_force
      
      subroutine switch(elrc,virlrc)
c***********************************************************************
c     
c     dl_poly routine for switching system in excitation simulation
c     copyright - daresbury laboratory
c     author    - w. smith   sep 2008
c     adapted from p.-a. cazade oct 2007
c     
c***********************************************************************
      
      real(8) :: swap,elrc,virlrc
      real(8), allocatable :: cou_sic_swp(:),qfix_swp(:)
      real(8), allocatable :: elrc_swp(:)
      
      allocate(cou_sic_swp(mxtmls_exc2),qfix_swp(mxtmls_exc))
      allocate(elrc_swp(mxtmls_exc2))
      
      swap=elrc
      elrc=elrc2
      elrc2=swap
      
      swap=engsic0
      engsic0=engsic2
      engsic2=swap
      
      swap=virlrc
      virlrc=virlrc2
      virlrc2=swap
      
      swap=elrc_sav
      elrc_sav=elrc2_sav
      elrc2_sav=swap

      swap=virlrc_sav
      virlrc_sav=virlrc2_sav
      virlrc2_sav=swap

      swap=qchg0
      qchg0=qchg2
      qchg2=swap
      
      cou_sic_swp(:)=cou_sol_sic(:)
      cou_sol_sic(:)=cou_exc_sic(:)
      cou_exc_sic(:)=cou_sic_swp(:)
      
      qfix_swp(:)=qfix_sol(:)
      qfix_sol(:)=qfix_exc(:)
      qfix_exc(:)=qfix_swp(:)
      
      elrc_swp(:)=elrc_sol(:)
      elrc_sol(:)=elrc_exc(:)
      elrc_exc(:)=elrc_swp(:)
      
      elrc_swp(:)=elrc_sol_sav(:)
      elrc_sol_sav(:)=elrc_exc_sav(:)
      elrc_exc_sav(:)=elrc_swp(:)
      
      deallocate(elrc_swp,cou_sic_swp,qfix_swp)
      
      return
      
      end subroutine switch
      
      subroutine lrcorrect_sol(lghost,volm)
      
c***********************************************************************
c     
c     dl_poly routine for excited state module
c     copyright - daresbury laboratory
c     author    - p.-a. cazade oct 2007
c     adapted   - w. smith     aug 2008
c     
c***********************************************************************
      
      logical :: lghost
      real(8) :: volm 
      
      elrc_sol(:)=elrc_sol_sav(:)*(volm_sav/volm)
      if(lghost)elrc_exc(:)=elrc_exc_sav(:)*(volm_sav/volm)
      
      return
      end subroutine lrcorrect_sol
      
      subroutine switch_atm(lfrmas)
      
c***********************************************************************
c     
c     dl_poly routine for excitation module
c     copyright - daresbury laboratory
c     author    - p.-a. cazade oct 2007
c     adapted   - w.smith      aug 2008
c     
c***********************************************************************
      
      logical lfrmas
      integer i,at1_swap,at2_swap
      
      at1_swap=ind_fre(1)
      at2_swap=ind_fre(2)
      ind_fre(1)=ind_fre(3)
      ind_fre(2)=ind_fre(4)
      ind_fre(3)=at1_swap
      ind_fre(4)=at2_swap
      
      if(lfrmas)then
        
        do i=ind_fre(1),ind_fre(2)
          
          atm_fre(i)=1
          weight(i)=lambda1*weight_sav(i)
          
        enddo
        
        do i=ind_fre(3),ind_fre(4)
          
          atm_fre(i)=2
          weight(i)=lambda2*weight_sav(i)
          
        enddo
        
      else
        
        do i=ind_fre(1),ind_fre(2)
          atm_fre(i)=1
        enddo
        
        do i=ind_fre(3),ind_fre(4)
          atm_fre(i)=2
        enddo
        
      endif
      
      return
      end subroutine switch_atm
      
      subroutine solvation_write
     x  (lexcite,lswitch,idnode,natms,nstep,nsolva,isolva,
     x  tstep,engunit,elrc)
      
c***********************************************************************
c     
c     dl_poly subroutine for writing solva file at selected
c     intervals in simulation
c     
c     copyright - daresbury  laboratory
c     author    - p.-a. cazade jun 2007
c     adapted   - w. smith     aug 2008
c     
c***********************************************************************
      
      implicit none
      
      logical newjob,lexcite,lswitch
      integer idnode,natms,nstep,nsolva,isolva,i,j,k
      integer mxtmls2,mxtmls3,mxtmls4
      real(8) tstep,engunit,elrc
      character*80 aa,bb,cc,dd
      
      save newjob
      data newjob/.true./
      
      mxtmls2=((mxtmls+1)*mxtmls)/2
      mxtmls3=(((mxtmls+3)*mxtmls+2)*mxtmls)/6
      mxtmls4=((((mxtmls+6)*mxtmls+11)*mxtmls+6)*mxtmls)/24
      
      if(idnode.eq.0)then
        
c     open the SOLVAT file if new job or file closed
        
        if(newjob)then
          
          newjob=.false.
          open(nsolwr,file='SOLVAT',position='append')
          
        endif
        
c     write file header block
        
        if(nstep.eq.nsolva.or.nstep.eq.isolva)then
          
          write(nsolwr,'(80a1)')cfgname
          
          if(abs(engunit-9648.530821d0).le.1.d-10)write(nsolwr,
     x      "('ENERGY UNITS=electron Volts ')")
          if(abs(engunit-418.4d0).le.1.d-10)      write(nsolwr,
     x      "('ENERGY UNITS=kcal/mol ')")
          if(abs(engunit-1.d2).le.1.d-10)         write(nsolwr,
     x      "('ENERGY UNITS=kjoule/mol ')")
          if(abs(engunit-boltz).lt.1.d-10)        write(nsolwr,
     x      "('ENERGY UNITS=kelvin ')")
          if(abs(engunit-1.d0).lt.1.d-10)         write(nsolwr,
     x      "('ENERGY UNITS=DL_POLY Internal UNITS ')")
          
          write(nsolwr,'(2i10)')natms,mxtmls
          write(nsolwr,'(1x,11a4)')' lex','lsw',' bnd',' ang',
     x      ' dih',' inv',' shl',' cou',' vdw',' 3bd',' 4bd'
          write(nsolwr,'(11l4)')lexcite,lswitch,lcomp
          
        endif
        
c     write out periodic data
        
        if(mod(nstep-nsolva,isolva).eq.0)then
          
c     mark start of time step data
          
          if(lexcite)then
            write(nsolwr,'("timestep",i10,f12.5,1p,2e14.6)')
     x        nstep,tstep,elrc/engunit,elrc2/engunit
          else
            write(nsolwr,'("timestep",i10,f12.5,1p,e14.6)')
     x        nstep,tstep,elrc/engunit
          endif
          
c     write intramolecular data
          
          write(nsolwr,'(1p,5e14.6)')temp_sol
          
          if(lcomp(1))then
            if(lexcite)then
              write(nsolwr,'(1p,5e14.6)')bnd_sol(:)/engunit
              write(nsolwr,'(1p,5e14.6)')bnd_exc(:)/engunit
            else
              write(nsolwr,'(1p,5e14.6)')bnd_sol(:)/engunit
            endif
          endif
          if(lcomp(2))then
            if(lexcite)then
              write(nsolwr,'(1p,5e14.6)')ang_sol(:)/engunit
              write(nsolwr,'(1p,5e14.6)')ang_exc(:)/engunit
            else
              write(nsolwr,'(1p,5e14.6)')ang_sol(:)/engunit
            endif
          endif
          if(lcomp(3))then
            if(lexcite)then
              write(nsolwr,'(1p,5e14.6)')dih_sol(:)/engunit
              write(nsolwr,'(1p,5e14.6)')dih_exc(:)/engunit
            else
              write(nsolwr,'(1p,5e14.6)')dih_sol(:)/engunit
            endif
          endif
          if(lcomp(4))then
            if(lexcite)then
              write(nsolwr,'(1p,5e14.6)')inv_sol(:)/engunit
              write(nsolwr,'(1p,5e14.6)')inv_exc(:)/engunit
            else
              write(nsolwr,'(1p,5e14.6)')inv_sol(:)/engunit
            endif
          endif
          
c     write core-shell data
          
          if(lcomp(5))then
            if(lexcite)then
              write(nsolwr,'(1p,5e14.6)')shl_sol(:)/engunit
              write(nsolwr,'(1p,5e14.6)')shl_exc(:)/engunit
            else
              write(nsolwr,'(1p,5e14.6)')shl_sol(:)/engunit
            endif
          endif
          
c     write coulombic data
          
          if(lcomp(6))then
            if(lexcite)then
              write(nsolwr,'(1p,5e14.6)')cou_sol(:)/engunit
              write(nsolwr,'(1p,5e14.6)')cou_exc(:)/engunit
            else
              write(nsolwr,'(1p,5e14.6)')cou_sol(:)/engunit
            endif
          endif
          
c     write vdw data
          
          if(lcomp(7))then
            if(lexcite)then
              write(nsolwr,'(1p,5e14.6)')vdw_sol(:)/engunit
              write(nsolwr,'(1p,5e14.6)')vdw_exc(:)/engunit
            else
              write(nsolwr,'(1p,5e14.6)')vdw_sol(:)/engunit
            endif
          endif
          
c     write 3-body data
          
          if(lcomp(8))then
            if(lexcite)then
              write(nsolwr,'(1p,5e14.6)')en3_sol(:)/engunit
              write(nsolwr,'(1p,5e14.6)')en3_exc(:)/engunit
            else
              write(nsolwr,'(1p,5e14.6)')en3_sol(:)/engunit
            endif
          endif
          
c     write 4-body data
          
          if(lcomp(9))then
            if(lexcite)then
              write(nsolwr,'(1p,5e14.6)')en4_sol(:)/engunit
              write(nsolwr,'(1p,5e14.6)')en4_exc(:)/engunit
            else
              write(nsolwr,'(1p,5e14.6)')en4_sol(:)/engunit
            endif
          endif
        
        endif
        
c     close SOLVAT file at regular intervals
        
        if(.not.newjob.and.mod(nstep,ndump).eq.0)then
          
          close(nsolwr)
          newjob=.true.
          
        endif
        
      endif
      
      return
      end subroutine solvation_write
      
      end module solvation_module
