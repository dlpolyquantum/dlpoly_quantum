      module bonds_module
      
c***********************************************************************
c     
c     dl_poly module for defining bond potential arrays
c     copyright - daresbury laboratory
c     
c     author    - w. smith     sep 2003
c     modified  - p.-a.cazade  oct 2007 : solvation etc.
c     modified  - d. quigley       2010 : metadynamics
c     
c     A quartic Morse bond potential is added
c
c     authors   - M.R.Momeni and F.A.Shakib 
c     copyright - M.R.Momeni and F.A. Shakib 2021
c
c     Method Development and Materials Simulation Laboratory
c
c***********************************************************************
      
      use config_module
      use error_module
      use metafreeze_module
      use parse_module
      use setup_module
      use site_module
      use solvation_module
      use utility_module
      
      implicit none
      
      real(8), allocatable :: prmbnd(:,:)
      integer, allocatable :: listbnd(:,:)
      integer, allocatable :: numbonds(:),keybnd(:),lstbnd(:,:)
      
      save prmbnd,listbnd,numbonds,keybnd,lstbnd
      
      contains
      
      subroutine alloc_bnd_arrays(idnode,mxnode)
      
      implicit none

      logical safe
      integer i,fail,idnode,mxnode
      dimension fail(5)

      safe=.true.
      
      do i=1,5
        fail(i)=0
      enddo
      
      allocate (prmbnd(mxtbnd,mxpbnd),stat=fail(1))
      allocate (numbonds(mxtmls),stat=fail(2))
      allocate (keybnd(mxtbnd),stat=fail(3))
      allocate (lstbnd(mxtbnd,3),stat=fail(4))
      allocate (listbnd(mxbond,4),stat=fail(5))
      
      if(any(fail.gt.0))safe=.false.
      if(mxnode.gt.1)call gstate(safe)
      if(.not.safe)call error(idnode,1030)
      
c     initialse numbonds array
      
      do i=1,mxtmls
        numbonds(i)=0
      enddo
      
      end subroutine alloc_bnd_arrays
      
      subroutine define_bonds
     x  (safe,idnode,itmols,nbonds,nsite,engunit)
      
c***********************************************************************
c     
c     dl_poly subroutine for defining bonds
c     
c     copyright - daresbury laboratory
c     author    - w. smith    aug 2003
c     
c***********************************************************************
      
      implicit none
      
      logical safe
      character*8 keyword
      character*1 message(80)
      integer idnode,itmols,nbonds,nsite,ntmp,ibond,keytmp
      integer iatm1,iatm2,isite1,isite2,idum,j
      real(8) engunit,parpot(mxpbnd)
      
      ntmp=intstr(record,lenrec,idum)
      numbonds(itmols)=numbonds(itmols)+ntmp

      if(idnode.eq.0)then
        write(nrite,"(/,1x,'number of chemical bonds',
     x    7x,i10)")numbonds(itmols)
        write(nrite,"(/,/,1x,'chemical bond details:',
     x    /,/,18x,'unit',5x,'key',5x,'index',5x,'index',28x,
     x    'parameters', /)")
      endif
      
      do ibond=1,ntmp

        nbonds=nbonds+1
        if(nbonds.gt.mxtbnd)call error(idnode,30)

c     read bond potential parameters

        call getrec(safe,idnode,nfield)
        if(.not.safe)return
        
        call copystring(record,message,80)
        call lowcase(record,4)
        call getword(keyword,record,4,lenrec)

        if(keyword(1:4).eq.'harm')then
          keytmp=1
        elseif(keyword(1:4).eq.'-hrm')then
          keytmp=-1
        elseif(keyword(1:4).eq.'mors')then
          keytmp=2
        elseif(keyword(1:4).eq.'-mrs')then
          keytmp=-2
        elseif(keyword(1:4).eq.'12-6')then
          keytmp=3
        elseif(keyword(1:4).eq.'-126')then
          keytmp=-3
        elseif(keyword(1:4).eq.'rhrm')then
          keytmp=4
        elseif(keyword(1:4).eq.'-rhm')then
          keytmp=-4
        elseif(keyword(1:4).eq.'quar')then
          keytmp=5
        elseif(keyword(1:4).eq.'-qur')then
          keytmp=-5
        elseif(keyword(1:4).eq.'buck')then
          keytmp=6
        elseif(keyword(1:4).eq.'-bck')then
          keytmp=-6
        elseif(keyword(1:4).eq.'fene')then
          keytmp=7
        elseif(keyword(1:4).eq.'-fen')then
          keytmp=-7
        elseif(keyword(1:4).eq.'coul')then
          keytmp=8
        elseif(keyword(1:4).eq.'-cou')then
          keytmp=-8
        elseif(keyword(1:2).eq.'lj')then
          keytmp=9
        elseif(keyword(1:3).eq.'-lj')then
          keytmp=-9
c *******************************************************************      
c     M.R.Momeni & F.A.Shakib
c     Method Development and Materials Simulation Laboratory
        elseif(keyword(1:4).eq.'qmor')then
          keytmp=10
        elseif(keyword(1:4).eq.'-qmr')then
          keytmp=-10
c *******************************************************************      
c *******************************************************************      
c     Dil Limbu
c     Method Development and Materials Simulation Laboratory
        elseif(keyword(1:4).eq.'mmst')then
          keytmp=11
        elseif(keyword(1:4).eq.'-mst')then
          keytmp=-11
c *******************************************************************      
        else
          if(idnode.eq.0)write(nrite,*)message
          call error(idnode,444)
        endif

        iatm1=intstr(record,lenrec,idum)
        iatm2=intstr(record,lenrec,idum)
        parpot(1)=dblstr(record,lenrec,idum)
        parpot(2)=dblstr(record,lenrec,idum)
        parpot(3)=dblstr(record,lenrec,idum)
        parpot(4)=dblstr(record,lenrec,idum)
        
        isite1=nsite-numsit(itmols)+iatm1
        isite2=nsite-numsit(itmols)+iatm2
        
c     test for frozen atom pairs
        
        if(lfzsit(isite1)*lfzsit(isite2).ne.0)then
          
          numbonds(itmols)=numbonds(itmols)-1
          if(idnode.eq.0)
     x      write(nrite,"(4x,a8,i10,4x,a4,2i10,2x,10f15.6)")
     x      '*frozen*',ibond,keyword(1:4),iatm1,iatm2,
     x      (parpot(j),j=1,mxpbnd)
          
        else
          
          if(idnode.eq.0)
     x      write(nrite,"(12x,i10,4x,a4,2i10,2x,10f15.6)")
     x      ibond,keyword(1:4),iatm1,iatm2,(parpot(j),j=1,mxpbnd)
          
        endif
        
        keybnd(nbonds)=keytmp
        lstbnd(nbonds,1)=iatm1
        lstbnd(nbonds,2)=iatm2
        prmbnd(nbonds,1)=parpot(1)
        prmbnd(nbonds,2)=parpot(2)
        prmbnd(nbonds,3)=parpot(3)
        prmbnd(nbonds,4)=parpot(4)
        
c     convert energy units to internal units
        
        if(abs(keytmp).eq.3)then
          prmbnd(nbonds,2)=prmbnd(nbonds,2)*engunit
        endif
        if(abs(keytmp).eq.5)then
          prmbnd(nbonds,3)=prmbnd(nbonds,3)*engunit
          prmbnd(nbonds,4)=prmbnd(nbonds,4)*engunit
        endif
        if(abs(keytmp).eq.6)then
          prmbnd(nbonds,3)=prmbnd(nbonds,3)*engunit
        endif
        if(abs(keytmp).ne.8)then
          prmbnd(nbonds,1)=prmbnd(nbonds,1)*engunit
        endif

      enddo
      
      return
      end subroutine define_bonds
      
      subroutine bndfrc
     x  (lsolva,lfree,lexcite,idnode,imcon,mxnode,ntbond,epsq,
     x  engbnd,virbnd)
      
c***********************************************************************
c     
c     dl_poly subroutine for calculating chemical bond energy and 
c     force terms in molecular dynamics.
c     
c     replicated data - blocked  data version
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith        july 1992
c     modified  - t. forester    march 1993 
c     modified  - t. forester    march 1994 
c     modified  - t. forester    may   1994 
c     modified  - t. forester    nov   1994 
c     modified  - w. smith       nov   2006
c     modified  - p.-a. cazade   oct   2007, solvation etc.
c     modified  - d. quigley           2007, metdynamics
c     
c     Potentials and forces of quartic Morse bond are added
c     copyright - M.R.Momeni and F.A. Shakib
c     authors   - M.R.Momeni and F.A.Shakib 2021
c
c***********************************************************************
      
      implicit none
      
      logical safe,lsolva,lfree,lexcite,lselect
      logical idrive,jdrive
      integer i,fail,ibnd1,ibnd2,idnode,mxnode,ii,ia,ib,imcon
      integer keyb,kk,ntbond
      real(8) strs(6),strs_loc(6)
      real(8) rab,rrab,omega,gamma,fx,fy,fz,engbnd,virbnd,epsq
      real(8), allocatable :: xdab(:),ydab(:),zdab(:)
      real(8) :: dr 
      
      safe=.true.
      
c     allocate work arrays
      
      fail=0
      allocate (xdab(msbad),ydab(msbad),zdab(msbad),stat=fail)
      if(fail.ne.0)safe=.false.
      if(mxnode.gt.1)call gstate(safe)
      if(.not.safe)call error(idnode,1040)
      
c     check size of work arrays
      
      if((ntbond-mxnode+1)/mxnode.gt.msbad)call error(idnode,418)
      
c     block indices
      
      ibnd1=(idnode*ntbond)/mxnode+1
      ibnd2=((idnode+1)*ntbond)/mxnode
      
c     initialise accumulators
      
      engbnd=0.d0
      virbnd=0.d0
      bnd_fre=0.d0
      bnd_vir=0.d0
      strs(:)=0.d0
      strs_loc(:)=0.d0

      if(lsolva)then
        
        lcomp(1)=.true.
        bnd_sol(:)=0.d0
        if(lexcite)bnd_exc(:)=0.d0
        
      endif
      
c     calculate atom separation vectors
      
      ii=0
      do i=ibnd1,ibnd2
        
        ii=ii+1
        
c     indices of bonded atoms
        
        ia=listbnd(ii,2)
        ib=listbnd(ii,3)
        
c     metadynamics local definitions
        
        if(lmetadyn)then
          
          idrive=driven(ltype(ia))
          jdrive=driven(ltype(ib))
          
        endif

c     components of bond vector
        
        xdab(ii)=xxx(ia)-xxx(ib)
        ydab(ii)=yyy(ia)-yyy(ib)
        zdab(ii)=zzz(ia)-zzz(ib)
        
      enddo
      
c     periodic boundary condition
      
      call images(imcon,0,1,ii,cell,xdab,ydab,zdab)
      
c     loop over all specified chemical bond potentials
      
      ii=0
      do i=ibnd1,ibnd2
        
        ii=ii+1
        
c     define components of bond vector
        
        rrab=0.d0
        rab=sqrt(xdab(ii)**2+ydab(ii)**2+zdab(ii)**2)
        if(rab.gt.1.d-6)rrab=1.d0/rab
        
c     index of potential function parameters
        
        kk=listbnd(ii,1)
        keyb=abs(keybnd(kk))
        
c     calculate scalar constant terms
        
        if(keyb.eq.0)then
          
c     null interaction
          
          omega=0.d0
          gamma=0.d0
          
        elseif(keyb.eq.1)then
          
c     harmonic potential
          
          omega=0.5d0*prmbnd(kk,1)*(rab-prmbnd(kk,2))**2
          gamma=prmbnd(kk,1)*(rab-prmbnd(kk,2))*rrab
          
        else if(keyb.eq.2)then
          
c     morse potential
          
          omega=prmbnd(kk,1)*((1.d0-exp(-prmbnd(kk,3)*
     x      (rab-prmbnd(kk,2))))**2-1.d0)
          gamma=2.d0*prmbnd(kk,1)*prmbnd(kk,3)*(1.d0-
     x      exp(-prmbnd(kk,3)*(rab-prmbnd(kk,2))))*
     x      exp(-prmbnd(kk,3)*(rab-prmbnd(kk,2)))*rrab
          
        else if(keyb.eq.3)then
          
c     12-6 potential
          
          omega=(prmbnd(kk,1)*rrab**6-prmbnd(kk,2))*rrab**6
          gamma=(6.d0*prmbnd(kk,2)-12.d0*prmbnd(kk,1)*rrab**6)*
     x      rrab**8
          
        elseif(keyb.eq.4)then
          
c     restrained harmonic
          
          rab=rab-prmbnd(kk,2)
          omega=0.5d0*prmbnd(kk,1)*(min(abs(rab),prmbnd(kk,3)))**2
     x      +prmbnd(kk,1)*prmbnd(kk,3)*max(abs(rab)-prmbnd(kk,3),0.d0)
          gamma=rrab*prmbnd(kk,1)*(sign(min(abs(rab),prmbnd(kk,3)),rab))
          
        elseif(keyb.eq.5)then
          
c     quartic potential
          
          omega=0.5d0*prmbnd(kk,1)*(rab-prmbnd(kk,2))**2+
     x      1.d0/3.d0*prmbnd(kk,3)*(rab-prmbnd(kk,2))**3+
     x      0.25d0*prmbnd(kk,4)*(rab-prmbnd(kk,2))**4
          gamma=rrab*(prmbnd(kk,1)*(rab-prmbnd(kk,2))+
     x      prmbnd(kk,3)*(rab-prmbnd(kk,2))**2+
     x      prmbnd(kk,4)*(rab-prmbnd(kk,2))**3)
          
        else if(keyb.eq.6)then
          
c     buckingham exp-6 potential
          
          omega=prmbnd(kk,1)*exp(-rab/prmbnd(kk,2))-prmbnd(kk,3)*
     x      rrab**6
          gamma=-rrab*prmbnd(kk,1)*exp(-rab/prmbnd(kk,2))/prmbnd(kk,2)+
     x      6.d0*prmbnd(kk,3)*rrab**8
          
        else if(keyb.eq.7)then
          
c     FENE bond potential
          
          omega=-0.5d0*prmbnd(kk,1)*prmbnd(kk,2)**2*log(1.d0-
     x      ((rab-prmbnd(kk,3))/prmbnd(kk,2))**2)
          gamma=rrab*prmbnd(kk,1)*(rab-prmbnd(kk,3))/
     x      (1.d0-((rab-prmbnd(kk,3))/prmbnd(kk,2))**2)
          
        else if(keyb.eq.8)then
          
c     coulomb bond potential
          
          omega=prmbnd(kk,1)*prmbnd(kk,2)*rrab*r4pie0/epsq
          gamma=-omega*rrab*rrab
          
        else if(keyb.eq.9)then
                    
c     lennard-jones potential
          
          omega=4.d0*prmbnd(kk,1)*(prmbnd(kk,2)/rab)**6*
     x      ((prmbnd(kk,2)/rab)**6-1.d0)
          gamma=-24.d0*prmbnd(kk,1)*(prmbnd(kk,2)/rab)**6*
     x      (2.d0*(prmbnd(kk,2)/rab)**6-1.d0)/rab**2

c *******************************************************************      
c     M.R.Momeni & F.A.Shakib
c     Method Development and Materials Simulation Laboratory

        else if(keyb.eq.10)then

c     quartic morse potential

         omega=prmbnd(kk,1)*((prmbnd(kk,3)**2*(rab-prmbnd(kk,2))**2)
     x         -(prmbnd(kk,3)**3*(rab-prmbnd(kk,2))**3)
     x         +(0.58333d0*prmbnd(kk,3)**4*(rab-prmbnd(kk,2))**4))
         gamma=(prmbnd(kk,1)*((2.d0*prmbnd(kk,3)**2*(rab-prmbnd(kk,2)))
     x         -(3.d0*prmbnd(kk,3)**3*(rab-prmbnd(kk,2))**2)
     x         +(2.33333d0*prmbnd(kk,3)**4*(rab-prmbnd(kk,2))**3)))
     x         *rrab

c *******************************************************************     
c     Dil Limbu 
c     Method Development and Materials Simulation Laboratory

        else if(keyb.eq.11)then

c     MM3 bond stretch potential
c      k0 is in mdyn unit 
         dr = (rab-prmbnd(kk,2))
         omega=71.94d0*prmbnd(kk,1)*dr**2*(1.d0-2.55d0*dr
     x         +(7d0/12.d0)*(2.55d0*dr)**2)

         gamma=71.94d0*prmbnd(kk,1)*dr*(2.d0-3.d0*2.55d0*dr
     x         +4.d0*(7.d0/12.d0)*(2.55d0*dr)**2)*rrab

c *******************************************************************     


        else
          
c     undefined potential
          
          omega=0.d0
          gamma=0.d0
          safe=.false.
          
        endif
        
c     indices of bonded atoms
        
        ia=listbnd(ii,2)
        ib=listbnd(ii,3)
        
c     set selection control
        
        lselect=.true.
        
        if(lexcite)then
          
c     selected excitation option
          
          if((atm_fre(ia).ne.1).and.(atm_fre(ib).ne.1))then
            
c     reset selection control
            
            lselect=(atm_fre(ia)+atm_fre(ib).eq.0)
            
            if(lsolva)then
              bnd_exc(atmolt(ia))=bnd_exc(atmolt(ia))+omega
            endif
            
          endif
          
        elseif(lfree)then
          
c     selected free energy option
          
          if((atm_fre(ia).eq.1).or.(atm_fre(ib).eq.1))then
            
c     set hamiltonian mixing parameter
            
            bnd_fre=bnd_fre-omega
            bnd_vir=bnd_vir-gamma*rab*rab
            omega=lambda1*omega
            gamma=lambda1*gamma
            
          elseif((atm_fre(ia).eq.2).or.(atm_fre(ib).eq.2))then
            
c     set hamiltonian mixing parameter
            
            bnd_fre=bnd_fre+omega
            bnd_vir=bnd_vir+gamma*rab*rab
            omega=lambda2*omega
            gamma=lambda2*gamma
            
          endif
          
        endif
        
        if(lselect)then
          
c     calculate bond energy and virial
        
          engbnd=engbnd+omega
          virbnd=virbnd+gamma*rab*rab
          
c     calculate solvation energy
        
          if(lsolva)then
            bnd_sol(atmolt(ia))=bnd_sol(atmolt(ia))+omega
          endif
          
c     calculate forces
          
          fx=-gamma*xdab(ii)
          fy=-gamma*ydab(ii)
          fz=-gamma*zdab(ii)
          
          fxx(ia)=fxx(ia)+fx
          fyy(ia)=fyy(ia)+fy
          fzz(ia)=fzz(ia)+fz
          
          fxx(ib)=fxx(ib)-fx
          fyy(ib)=fyy(ib)-fy
          fzz(ib)=fzz(ib)-fz
          
c     calculate stress tensor
        
          strs(1)=strs(1)+xdab(ii)*fx
          strs(2)=strs(2)+xdab(ii)*fy
          strs(3)=strs(3)+xdab(ii)*fz
          strs(4)=strs(4)+ydab(ii)*fy
          strs(5)=strs(5)+ydab(ii)*fz
          strs(6)=strs(6)+zdab(ii)*fz
          
        endif
        
c     metadynamics local parameters
        
        if(lmetadyn.and.(idrive.or.jdrive))then
          
c     local energy and virial
          
          eng_loc=eng_loc+omega
          vir_loc=vir_loc+gamma*rab*rab
          
c     local forces
          
          fxx_loc(ia)=fxx_loc(ia)+fx
          fyy_loc(ia)=fyy_loc(ia)+fy
          fzz_loc(ia)=fzz_loc(ia)+fz
          
          fxx_loc(ib)=fxx_loc(ib)-fx
          fyy_loc(ib)=fyy_loc(ib)-fy
          fzz_loc(ib)=fzz_loc(ib)-fz
          
c     local stress tensor
          
          strs_loc(1)=strs_loc(1)+xdab(ii)*fx
          strs_loc(2)=strs_loc(2)+xdab(ii)*fy
          strs_loc(3)=strs_loc(3)+xdab(ii)*fz
          strs_loc(4)=strs_loc(4)+ydab(ii)*fy
          strs_loc(5)=strs_loc(5)+ydab(ii)*fz
          strs_loc(6)=strs_loc(6)+zdab(ii)*fz
          
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
      
c     check for undefined potentials
      
      if(mxnode.gt.1)call gstate(safe)
      if(.not.safe)call error(idnode,444)
      
c     sum contributions to potential and virial
      
      if(mxnode.gt.1)then
        
        buffer(1)=engbnd
        buffer(2)=virbnd
        buffer(3)=bnd_fre
        buffer(4)=bnd_vir
        call gdsum(buffer(1),4,buffer(5))        
        engbnd=buffer(1)
        virbnd=buffer(2)
        bnd_fre=buffer(3)
        bnd_vir=buffer(4)
        
        if(lsolva)then
          
          call gdsum(bnd_sol(1),mxtmls,buffer(1))
          if(lexcite)call gdsum(bnd_exc(1),mxtmls,buffer(1))
          
        endif
        
      endif
      
      deallocate (xdab,ydab,zdab,stat=fail)
      
      return
      end subroutine bndfrc
      
      end module bonds_module
