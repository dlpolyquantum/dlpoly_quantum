      module external_field_module
      
c***********************************************************************
c     
c     dl_poly module for defining external field potential arrays
c     copyright - daresbury laboratory
c     author    - w. smith    oct 2003
c     
c     Two repulsive Lennard-Jones walls at zmin and max are added
c
c     copyright - M.R.Momeni and F.A. Shakib
c     authors   - M.R.Momeni and F.A.Shakib 2021
c
c     Method Development and Materials Simulation Laboratory
c
c***********************************************************************
      
      use config_module
      use error_module
      use parse_module
      use setup_module
      use utility_module
      
      implicit none
      
      real(8), allocatable :: prmfld(:)
      
      save prmfld
      
      contains
      
      subroutine alloc_fld_arrays(idnode,mxnode)
      
      implicit none
      
      logical safe
      integer fail,idnode,mxnode
      
      safe=.true.
      
c     allocate arrays
      
      fail=0
      
      allocate (prmfld(mxfld),stat=fail)
      if(fail.gt.0 )safe=.false.      
      if(mxnode.gt.1)call gstate(safe)    
      if(.not.safe)call error(idnode,1200)
      
      end subroutine alloc_fld_arrays
      
      subroutine define_external_field
     x  (safe,lunits,idnode,keyfld,engunit)
      
c***********************************************************************
c     
c     dl_poly subroutine to define external fields
c     
c     copyright - daresbury laboratory
c     author    - w. smith    aug 2003
c     amended   - p.-l. chau  jun 2009 z-restraint option
c     
c***********************************************************************
      
      implicit none
      
      logical safe,lunits
      character*8 keyword
      character*1 message(80)
      integer idnode,keyfld,nfld,i,k,idum
      real(8) engunit
      
      call getrec(safe,idnode,nfield)
      if(.not.safe)return
      
      call strip(record,lenrec)
      call lowcase(record,lenrec)
      call copystring(record,message,80)
      call getword(keyword,record,4,lenrec)
      
      if(keyword(1:4).eq.'elec') then
        keyfld=1 
      elseif(keyword(1:4).eq.'oshr') then
        keyfld=2
      elseif(keyword(1:4).eq.'shrx') then
        keyfld=3
      elseif(keyword(1:4).eq.'grav') then
        keyfld=4
      elseif(keyword(1:4).eq.'magn') then
        keyfld=5
      elseif(keyword(1:4).eq.'sphr') then
        keyfld=6
      elseif(keyword(1:4).eq.'zbnd') then
        keyfld=7
      elseif(keyword(1:4).eq.'zres') then
        keyfld=9
c *******************************************************************      
c     M.R.Momeni & F.A.Shakib
c     Method Development and Materials Simulation Laboratory
      elseif(keyword(1:4).eq.'ljwl') then
        keyfld=10
c *******************************************************************      
      else
        if(idnode.eq.0) write(nrite,*) message
        call error(idnode,454)
      endif
      
      do i=1,mxfld
        prmfld(i)=0.d0
      enddo
      
      nfld=intstr(record,lenrec,idum)
      if(nfld.eq.0)nfld=5
      call getrec(safe,idnode,nfield)
      if(.not.safe)return
      do k=1,nfld
        
        prmfld(k)=dblstr(record,lenrec,idum)
        if(idum.gt.lenrec.and.k.lt.nfld)then
          
          call getrec(safe,idnode,nfield)
          if(.not.safe)return
          
        endif
        
      enddo
      
      if(idnode.eq.0) then
        
        write(nrite,"(/,/,1x,'external field key ',13x,a4,
     x    /,/,30x,'external field parameters')") keyword(1:4)
        write(nrite,"(2(/,1x,1p,5e15.5))") prmfld
        
      endif      
      
c     convert to internal units
      
      if(keyfld.eq.1.or.keyfld.eq.4.or.keyfld.eq.5) then
        
        if(.not.lunits)call error(idnode,6)
        
        do i=1,3
          prmfld(i)=prmfld(i)*engunit
        enddo
        
      elseif(keyfld.eq.2.or.keyfld.eq.6.or.keyfld.eq.7.or.
     x      keyfld.eq.10) then
        
        prmfld(1)=prmfld(1)*engunit
        
      elseif(keyfld.eq.9) then
        
        prmfld(3)=prmfld(3)*engunit
        
      endif
      
      return
      end subroutine define_external_field
      
      subroutine extnfld
     x  (idnode,imcon,keyfld,mxnode,natms,engfld,virfld)
      
c***********************************************************************
c     
c     dl_poly routine for application of an external field
c     
c     replicated data version / block data
c     
c     copyright daresbury laboratory 1993
c     author  -    t.forester october 1993
c     amended -    t.forester dec 1994
c     amended -    p.-l. chau jun 2009 z-restraint option
c     
c     Two 12-6 Lennard-Jones walls at zmin and zmax are added
c     copyright - M.R.Momeni and F.A. Shakib
c     authors    - M.R.Momeni and F.A.Shakib 2021
c
c***********************************************************************
      
      implicit none
      
      integer idnode,imcon,keyfld,mxnode,natms,iatm1,iatm2,i,j
      integer istart,ifinish,numresat
      real(8) engfld,virfld,rz,rrr,gamma,zdif,totwgt
      real(8) com(3)
      real(8) zdif_min,ljf_min,ljp_min 
      real(8) zdif_max,ljf_max,ljp_max 

c     energy and virial accumulators 
      
      engfld=0.d0
      virfld=0.d0
      
c     block indices
      
      iatm1=(idnode*natms)/mxnode+1
      iatm2=((idnode+1)*natms)/mxnode
      
      if(keyfld.eq.1) then
        
c     electric field: prmfld(1-3) are field components
        
        do i=iatm1,iatm2
          
          fxx(i)=fxx(i)+chge(i)*prmfld(1)
          fyy(i)=fyy(i)+chge(i)*prmfld(2)
          fzz(i)=fzz(i)+chge(i)*prmfld(3)
          
        enddo
        
      elseif(keyfld.eq.2) then
        
c     oscillating shear: orthorhombic box:  Fx=a*cos(b.2.pi.z/L)
        
        rz=2.d0*pi/cell(9)
        
        do i=iatm1,iatm2
          
          fxx(i)=fxx(i)+prmfld(1)*cos(prmfld(2)*zzz(i)*rz)
          
        enddo
        
      elseif(keyfld.eq.3.and.imcon.eq.6) then
        
c     continuous shear of walls : 2D periodic box (imcon=6)
c     shear rate=prmfld(1) angstrom per ps for atoms at
c     abs(z) > prmfld(2)
        
        do i=iatm1,iatm2
          
          if(abs(zzz(i)).gt.prmfld(2)) then
            
            vxx(i)=0.5d0*sign(prmfld(1),zzz(i))
            
          endif
          
        enddo
        
      elseif(keyfld.eq.4) then
        
c     gravitational field: field components given by prmfld(1-3)
        
        do i=iatm1,iatm2
          
          fxx(i)=fxx(i)+prmfld(1)*weight(i)
          fyy(i)=fyy(i)+prmfld(2)*weight(i)
          fzz(i)=fzz(i)+prmfld(3)*weight(i)
          
        enddo
        
      elseif(keyfld.eq.5) then
        
c     magnetic field: field components given by prmfld(1-3)
        
        do i=iatm1,iatm2
          
          fxx(i)=fxx(i)+(vyy(i)*prmfld(3)-vzz(i)*prmfld(2))
     x      *chge(i)
          fyy(i)=fyy(i)+(vzz(i)*prmfld(1)-vxx(i)*prmfld(3))
     x      *chge(i)
          fzz(i)=fzz(i)+(vxx(i)*prmfld(2)-vyy(i)*prmfld(1))
     x      *chge(i)
          
        enddo
        
      elseif(keyfld.eq.6) then
        
c     containing sphere : r^(-n) potential
        
        do i=iatm1,iatm2
          
          rrr=sqrt(xxx(i)**2+yyy(i)**2+zzz(i)**2)
          if(rrr.gt.prmfld(4)) then
            rrr=prmfld(2)-rrr
            if(rrr.lt.0.d0) rrr=0.1d0
            
            gamma =prmfld(1)*rrr**(-prmfld(3))
            engfld=engfld+gamma
            
            gamma=-prmfld(3)*gamma/((prmfld(2)-rrr)*rrr)
            
            fxx(i)=fxx(i)+gamma*xxx(i)
            fyy(i)=fyy(i)+gamma*yyy(i)
            fzz(i)=fzz(i)+gamma*zzz(i)
            
          endif
          
        enddo
        
      elseif(keyfld.eq.7) then
        
c     repulsive wall (harmonic) starting at z0
        
        do i=iatm1,iatm2
          
          if(prmfld(3)*zzz(i).gt.prmfld(3)*prmfld(2)) then
            
            zdif=zzz(i)-prmfld(2)
            gamma=-prmfld(1)*zdif
            
            fzz(i)=fzz(i)+gamma
            engfld=engfld-gamma*zdif/2.
            
          endif
          
        enddo
        
      elseif(keyfld.eq.9) then
        
c     keyfld=9. restrain molecule z-position
c     prmfld(1) is number of first atom of restrained molecule
c     prmfld(2) is number of last atom of restrained molecule
c     prmfld(3) is the restraining constant
c     prmfld(4) is z-min
c     prmfld(5) is z-max
        
        istart=nint(prmfld(1))
        ifinish=nint(prmfld(2))
        numresat=ifinish-istart+1
        
c     calculate the centre of mass of the molecule
        
        call getcom_mol(istart,ifinish,imcon,idnode,mxnode,totwgt,com)
        
c     apply restraint force according to location
        
        if(com(3).lt.prmfld(4))then
          
c     if centre of mass is below z-min, activate restraining force
          
          do i=istart,ifinish
            
            fzz(i)=fzz(i)-prmfld(3)*(weight(i)/totwgt)*
     x        (com(3)-prmfld(4))/mxnode
            
          enddo
          
        elseif(com(3).gt.prmfld(5))then
          
c     if centre of mass if above z-max, activate restraining force
          
          do i=istart,ifinish
            
            fzz(i)=fzz(i)-prmfld(3)*(weight(i)/totwgt)*
     x        (com(3)-prmfld(5))/mxnode
            
          enddo
          
        endif
        
c **********************************************************************      
c     M.R.Momeni & F.A.Shakib
c     Method Development and Materials Simulation Laboratory

      elseif(keyfld.eq.10) then

c     Two repulsive walls (Lennard-Jones) starting at zmin and zmax
        
        do i=iatm1,iatm2
          
          zdif_min=zzz(i)-prmfld(3)
          zdif_max=zzz(i)-prmfld(4)
          
          ljf_min=4.d0*prmfld(1)*((-12.d0*(prmfld(2)**12)
     x           /(zdif_min**13))+(6.d0*(prmfld(2)**6)/(zdif_min**7)))
          ljf_max=4.d0*prmfld(1)*((-12.d0*(prmfld(2)**12)
     x           /(zdif_max**13))+(6.d0*(prmfld(2)**6)/(zdif_max**7)))

          fzz(i)=fzz(i)-ljf_min-ljf_max

          ljp_min=4.d0*prmfld(1)*((prmfld(2)/zdif_min)**12
     x           -(prmfld(2)/zdif_min)**6)
          ljp_max=4.d0*prmfld(1)*((prmfld(2)/zdif_max)**12
     x           -(prmfld(2)/zdif_max)**6)

        enddo
        
      else
 
c **********************************************************************      

c     unidentified field potential error exit
        
        call error(idnode,454)
        
      endif
      
c     global sum of external field potential and virial
      
      if(mxnode.gt.1)then
        
        buffer(1)=engfld
        buffer(2)=virfld
        call gdsum(buffer(1),2,buffer(3))
        engfld=buffer(1)
        virfld=buffer(2)
        
      endif
      
      return
      end subroutine extnfld
      
      end module external_field_module
