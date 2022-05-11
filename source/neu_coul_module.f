      module neu_coul_module
      
c***********************************************************************
c     
c     dl_poly module for defining neutral group coulomb terms
c     copyright - daresbury laboratory
c     
c     author  - w. smith     aug 2006
c     adapted - p.-a. cazade oct 2007 : solvation, free energy etc
c     adapted - w. smith     jan 2011 : metadynamics
c     
c***********************************************************************
      
      use config_module
      use error_module
      use ewald_module
      use metafreeze_module
      use pair_module
      use setup_module
      use solvation_module
      
      contains
      
      subroutine coul0neu
     x  (lsolva,lfree,lexcite,ik,engcpe,vircpe,epsq)
      
c***********************************************************************
c     
c     dl_poly subroutine for calculating coulombic force.
c     1/r potential, no truncation or damping.
c     neutral group implementation
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1994
c     author    - t. forester    march 1994
c     adapted for solvation, free energy and excitation
c     - p.-a. cazade oct 2007
c     
c***********************************************************************
      
      implicit none
      
      logical lsolva,lfree,lexcite,lselect,lskip,idrive,jdrive
      integer ik,m,iatm,jatm,kkk
      real(8) engcpe,vircpe,epsq,chgprd,rsq,rrr,coul,reps,fcoul
      real(8) fx,fy,fz
      real(8) strs(6),strs_loc(6)
      
      lskip=(lfree.or.lexcite)
      
c     initialise stress tensor accumulators
      
      strs(:)=0.d0
      strs_loc(:)=0.d0
      
c     initialise potential energy and virial
      
      engcpe=0.d0
      vircpe=0.d0
      
c     start of primary loop for forces evaluation
      
      reps=r4pie0/epsq
      do m=1,ik
        
c     atomic index and charge product
        
        iatm=ilist(m)
        jatm=jlist(m)
        
c     metadynamics local definitions

        if(lmetadyn)then

          idrive=driven(ltype(iatm))
          jdrive=driven(ltype(jatm))

        endif

        if(lskip)then
          if(atm_fre(iatm)*atm_fre(jatm).eq.2)cycle
        endif
        
        chgprd=chge(jatm)*chge(iatm)*reps
        
        if(abs(chgprd).gt.1.d-10)then
          
          rsq=rsqdf(m)
          rrr=sqrt(rsq)
          
c     calculate coulomb energy and force
          
          coul=chgprd/rrr
          fcoul=coul/rsq
          
c     set selection control
          
          lselect=.true.
          
c     set double index
          
          if(lsolva)kkk=loc2(atmolt(iatm),atmolt(jatm))
          
          if(lexcite)then
            
c     selected excitation option
            
            if((atm_fre(iatm).ne.1).and.(atm_fre(jatm).ne.1))then
              
c     reset selection control
              
              lselect=(atm_fre(iatm)+atm_fre(jatm).eq.0)
              
c     calculate solvation energy
                  
              if(lsolva)cou_exc(kkk)=cou_exc(kkk)+coul
              
            endif
            
          elseif(lfree)then
            
c     selected free energy option
            
            if((atm_fre(iatm).eq.1).or.(atm_fre(jatm).eq.1))then
              
c     set hamiltonian mixing parameter
              
              cou_fre=cou_fre-coul
              cou_vir=cou_vir+coul
              coul=lambda1*coul
              fcoul=lambda1*fcoul
              
            elseif((atm_fre(iatm).eq.2).or.(atm_fre(jatm).eq.2))then
              
c     set hamiltonian mixing parameter
              
              cou_fre=cou_fre+coul
              cou_vir=cou_vir-coul
              coul=lambda2*coul
              fcoul=lambda2*fcoul
              
            endif
            
          endif
          
          if(lselect)then
            
c     calculate potential energy and virial
          
            engcpe=engcpe+coul
            vircpe=vircpe-coul
            
c     calculate solvation energy
          
            if(lsolva)cou_sol(kkk)=cou_sol(kkk)+coul
            
c     calculate forces
            
            fx=fcoul*xdf(m)
            fy=fcoul*ydf(m)
            fz=fcoul*zdf(m)
            
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
          
c     metadynamics local definitions

          if(lmetadyn.and.(idrive.or.jdrive))then
            
c     local energy and virial
          
            eng_loc=eng_loc+coul
            vir_loc=vir_loc-coul
            
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
      end subroutine coul0neu
      
      subroutine coul2neu
     x  (lsolva,lfree,lexcite,ik,engcpe,vircpe,epsq)
      
c***********************************************************************
c     
c     dl_poly subroutine for calculating coulombic forces 
c     assuming a distance dependant dielectric `constant'.
c     neutral group implementation
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1994
c     author    - t. forester    march 1994
c     adapted for solvation, free energy and excitation
c     - p.-a. cazade oct 2007
c     
c***********************************************************************
      
      implicit none
      
      logical lsolva,lfree,lexcite,lselect,lskip,idrive,jdrive
      integer ik,m,iatm,jatm,kkk
      real(8) engcpe,vircpe,epsq,fx,fy,fz,chgprd,rrsq,coul,egamma
      real(8) strs(6),strs_loc(6)
      
      lskip=(lfree.or.lexcite)
      
c     initialise stress tensor accumulators
      
      strs(:)=0.d0
      strs_loc(:)=0.d0
      
c     initialise potential energy and virial
      
      engcpe=0.d0
      vircpe=0.d0
      
c     start of primary loop for forces evaluation
      
      do m=1,ik
        
c     atomic index and charge product
        
        iatm=ilist(m)
        jatm=jlist(m)
        
c     metadynamics local definitions

        if(lmetadyn)then

          idrive=driven(ltype(iatm))
          jdrive=driven(ltype(jatm))

        endif

        if(lskip)then
          if(atm_fre(iatm)*atm_fre(jatm).eq.2)cycle
        endif
        
        chgprd=chge(iatm)*chge(jatm)*r4pie0/epsq
        
        if(abs(chgprd).gt.1.d-10)then
          
c     calculate potential energy
          
          rrsq=1.d0/rsqdf(m)
          coul=chgprd*rrsq          
          egamma=2.d0*coul*rrsq
          
c     set selection control
          
          lselect=.true.
          
c     set double index
          
          if(lsolva)kkk=loc2(atmolt(iatm),atmolt(jatm))
          
          if(lexcite)then
            
c     selected excitation option
            
            if((atm_fre(iatm).ne.1).and.(atm_fre(jatm).ne.1))then
              
c     reset selection control
              
              lselect=(atm_fre(iatm)+atm_fre(jatm).eq.0)
              
c     calculate solvation energy
                  
              if(lsolva)cou_exc(kkk)=cou_exc(kkk)+coul
              
            endif
            
          elseif(lfree)then
            
c     selected free energy option
            
            if((atm_fre(iatm).eq.1).or.(atm_fre(jatm).eq.1))then
              
c     set hamiltonian mixing parameter
              
              cou_fre=cou_fre-coul
              cou_vir=cou_vir+2.d0*coul
              coul=lambda1*coul
              egamma=lambda1*egamma
              
            elseif((atm_fre(iatm).eq.2).or.(atm_fre(jatm).eq.2))then
              
c     set hamiltonian mixing parameter
              
              cou_fre=cou_fre+coul
              cou_vir=cou_vir-2.d0*coul
              coul=lambda2*coul
              egamma=lambda2*egamma
              
            endif
            
          endif
          
          if(lselect)then
            
c     calculate potential energy and Virial
            
            engcpe=engcpe+coul
            vircpe=vircpe-2.d0*coul
          
c     calculate solvation energy
          
            if(lsolva)cou_sol(kkk)=cou_sol(kkk)+coul
            
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
          
c     metadynamics local definitions

          if(lmetadyn.and.(idrive.or.jdrive))then
            
c     local energy and virial

            eng_loc=eng_loc+coul
            vir_loc=vir_loc-2.d0*coul
          
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
      end subroutine coul2neu
      
      subroutine coul3neu
     x  (lsolva,lfree,lexcite,ik,engcpe,vircpe,epsq,rcut,alpha)
      
c***********************************************************************
c     
c     dl_poly subroutine for calculating coulombic force.
c     reaction field  potential
c     Ref: M Neumann, J Chem Phys, 82, 5633, (1985)
c     adapted for fennell-gezelter coulombic model
c     by w.smith june 2007
c     Ref: CJ Fennell and JD Gezelter, J Chem Phys, 
c     124, 234104, (2006)
c     adapted for solvation, free energy and excitation
c               - p.-a. cazade oct 2007
c     
c     neutral group implementation
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1995
c     author    - t. forester february 1995
c     stress tensor - t.forester   feb 1995
c     
c***********************************************************************
      
      implicit none
      
      logical lsolva,lfree,lexcite,lselect,lskip,newjob,idrive,jdrive
      integer ik,m,iatm,jatm,kkk
      real(8) engcpe,vircpe,rcut,epsq,vcon,fcon,rdr,ppp,erc1,fer1
      real(8) rcsq,b0,rfld0,rfld1,rfld2,chgprd,rsq,coul,omega,fcoul
      real(8) fx,fy,fz,rrr,alpha,a1,a2,a3,a4,a5,pp,tt,exp1
      real(8) strs(6),strs_loc(6)

      save newjob,b0,rfld0,rfld1,rfld2,vcon,fcon

      data a1,a2,a3/0.254829592d0,-0.284496736d0,1.421413741d0/
      data a4,a5,pp/-1.453152027d0,1.061405429d0,0.3275911d0/
      data newjob/.true./
      
      if(newjob)then
        
c     reaction field terms
        
        b0=2.d0*(epsq-1.d0)/(2.d0*epsq+1.d0)
        rfld0=b0/rcut**3
        rfld1=(1.d0+b0*0.5d0)/rcut
        rfld2=rfld0*0.5d0
        
c     screened coulomb terms
        
        tt=1.d0/(1.d0+pp*alpha*rcut)
        exp1=exp(-(alpha*rcut)**2)
        erc1=tt*(a1+tt*(a2+tt*(a3+tt*(a4+tt*a5))))*exp1/rcut
        fer1=(erc1+2.d0*(alpha/sqrpi)*exp1)/(rcut*rcut)
        vcon=erc1+rfld2*rcut**2-rfld1
        fcon=rcut*fer1-rfld0*rcut
        
      endif
      
      lskip=(lfree.or.lexcite)
      
c     initialise stress tensor accumulators
      
      strs(:)=0.d0
      strs_loc(:)=0.d0
      
c     initialise potential energy and virial
      
      engcpe=0.d0
      vircpe=0.d0     
      
c     start of primary loop for forces evaluation
      
      do m=1,ik
        
c     atomic index and charge product
        
        iatm=ilist(m)
        jatm=jlist(m)
        
c     metadynamics local definitions

        if(lmetadyn)then

          idrive=driven(ltype(iatm))
          jdrive=driven(ltype(jatm))

        endif

        if(lskip)then
          if(atm_fre(iatm)*atm_fre(jatm).eq.2)cycle
        endif
        
        chgprd=chge(jatm)*chge(iatm)*r4pie0
        if(abs(chgprd).gt.1.d-10)then
          
c     calculate interatomic distance
          
          rsq=rsqdf(m)
          rrr=sqrt(rsq)
          
c     error function terms
          
          tt=1.d0/(1.d0+pp*alpha*rrr)
          exp1=exp(-(alpha*rrr)**2)
          erc1=tt*(a1+tt*(a2+tt*(a3+tt*(a4+tt*a5))))*exp1/rrr
          fer1=(erc1+2.d0*(alpha/sqrpi)*exp1)/rsq
          
c     calculate potential energy
          
          omega=erc1-vcon+fcon*(rrr-rcut)
          coul=chgprd*(omega+rfld2*rsq-rfld1)
          
c     calculate forces
          
          fcoul=chgprd*(fer1-fcon/rrr-rfld0)
          
c     set selection control
          
          lselect=.true.
          
c     set double index
          
          if(lsolva)kkk=loc2(atmolt(iatm),atmolt(jatm))
          
          if(lexcite)then
            
c     selected excitation option
            
            if((atm_fre(iatm).ne.1).and.(atm_fre(jatm).ne.1))then
              
c     reset selection control
              
              lselect=(atm_fre(iatm)+atm_fre(jatm).eq.0)
              
c     calculate solvation energy
                  
              if(lsolva)cou_exc(kkk)=cou_exc(kkk)+coul
              
            endif
            
          elseif(lfree)then
            
c     selected free energy option
            
            if((atm_fre(iatm).eq.1).or.(atm_fre(jatm).eq.1))then
              
c     set hamiltonian mixing parameter
              
              cou_fre=cou_fre-coul
              cou_vir=cou_vir+fcoul*rsq
              coul=lambda1*coul
              fcoul=lambda1*fcoul
              
            elseif((atm_fre(iatm).eq.2).or.(atm_fre(jatm).eq.2))then
              
c     set hamiltonian mixing parameter
              
              cou_fre=cou_fre+coul
              cou_vir=cou_vir-fcoul*rsq
              coul=lambda2*coul
              fcoul=lambda2*fcoul
              
            endif
            
          endif
          
          if(lselect)then
            
c     calculate coulombic energy and virial
            
            engcpe=engcpe+coul
            vircpe=vircpe-fcoul*rsq
            
c     calculate solvation energy
            
            if(lsolva)cou_sol(kkk)=cou_sol(kkk)+coul
            
c     calculate coulombic force
            
            fx=fcoul*xdf(m)
            fy=fcoul*ydf(m)
            fz=fcoul*zdf(m)
            
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
          
c     metadynamics local definitions

          if(lmetadyn.and.(idrive.or.jdrive))then
            
c     local energy and virial

            eng_loc=eng_loc+coul
            vir_loc=vir_loc-fcoul*rsq
          
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
      end subroutine coul3neu
      
      end module neu_coul_module
