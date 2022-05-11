      module inversion_module

c***********************************************************************
c     
c     dl_poly module for defining inversion potential arrays
c     copyright - daresbury laboratory
c     
c     author  - w. smith  sep 2003
c     adapted - w. smith  jun 2008 : solvation, free energy, excitation
c     adapted - w. smith  jan 2011 : metadynamics
c     
c***********************************************************************

      use config_module
      use error_module
      use metafreeze_module
      use parse_module
      use property_module
      use setup_module
      use site_module
      use solvation_module

      implicit none

      real(8), allocatable :: prminv(:,:)
      integer, allocatable :: listinv(:,:)
      integer, allocatable :: numinv(:),keyinv(:),lstinv(:,:)

      save prminv,listinv,numinv,keyinv,lstinv

      contains

      subroutine alloc_inv_arrays(idnode,mxnode)
      
      implicit none
      
      logical safe
      integer i,fail,idnode,mxnode
      dimension fail(5)

      safe=.true.

c     allocate arrays

      fail(:)=0

      allocate (prminv(mxtinv,mxpinv),stat=fail(1))
      allocate (numinv(mxtmls),stat=fail(2))
      allocate (keyinv(mxtinv),stat=fail(3))
      allocate (lstinv(mxtinv,4),stat=fail(4))
      allocate (listinv(mxinv,5),stat=fail(5))

      if(any(fail.gt.0))safe=.false.      
      if(mxnode.gt.1)call gstate(safe)    
      if(.not.safe)call error(idnode,1120)
      
c     initialise numin array
      
      do i=1,mxtmls
         numinv(i)=0
      enddo

      end subroutine alloc_inv_arrays

      subroutine define_inversions
     x  (safe,idnode,itmols,ninver,nsite,engunit)

c***********************************************************************
c     
c     dl_poly subroutine for defining inversion angle potentials
c     note: includes calcite planar potential which is not strictly
c     an inversion potential and is different in implementation
c     
c     copyright - daresbury laboratory
c     author    - w. smith    aug 2003
c     
c***********************************************************************

      implicit none

      logical safe
      character*8 keyword
      character*1 message(80)
      integer idnode,itmols,ninver,nsite,ntmp,inv,j,keytmp
      integer iatm1,iatm2,iatm3,iatm4,isite1,isite2,isite3,isite4
      integer ia,ja,idum
      real(8) engunit,parpot(mxpinv)

      ntmp=intstr(record,lenrec,idum)
      numinv(itmols)=numinv(itmols)+ntmp

      if(idnode.eq.0)then
        write(nrite,"(/,1x,'number of inversion terms',
     x    6x,i10)")numinv(itmols)
        Write(nrite,"(/,1x,'inversion angle details:',
     x       /,/,12x,'unit',5x,'key',5x,'index',5x,'index',5x,'index',
     x       5x,'index',7x,'f-const',8x,'angle',8x,'factor',/)")
      endif
      
      do inv=1,ntmp
        
        ninver=ninver+1
        if(ninver.gt.mxtinv)call error(idnode,73)
        
c     read inversion potential parameters
        
        call getrec(safe,idnode,nfield)
        if(.not.safe)return
        
        call copystring(record,message,80)
        call lowcase(record,4)
        call getword(keyword,record,4,lenrec)
        
        if(keyword(1:4).eq.'harm')then
          keytmp=1
        elseif(keyword(1:4).eq.'hcos')then
          keytmp=2
        elseif(keyword(1:4).eq.'plan')then
          keytmp=3
        elseif(keyword(1:4).eq.'calc')then
          keytmp=4
        else
          if(idnode.eq.0)write(nrite,*)record
          call error(idnode,449)
        endif
        
        iatm1=intstr(record,lenrec,idum)
        iatm2=intstr(record,lenrec,idum)
        iatm3=intstr(record,lenrec,idum)
        iatm4=intstr(record,lenrec,idum)
        parpot(1)=dblstr(record,lenrec,idum)
        parpot(2)=dblstr(record,lenrec,idum)
        
        isite1=nsite-numsit(itmols)+iatm1
        isite2=nsite-numsit(itmols)+iatm2
        isite3=nsite-numsit(itmols)+iatm3
        isite4=nsite-numsit(itmols)+iatm4
        
c     test for frozen atom pairs
        
        if(lfzsit(isite1)*lfzsit(isite2)*
     x    lfzsit(isite3)*lfzsit(isite4).ne.0)then
          
          numinv(itmols)=numinv(itmols)-1
          if(idnode.eq.0)
     x      write(nrite,"(4x,a8,i10,4x,a4,4i10,10f15.6)")
     x      '*frozen*',inv,keyword(1:4),iatm1,iatm2,iatm3,iatm4,
     x      (parpot(j),j=1,mxpinv)
          
        else
          
          if(idnode.eq.0)
     x      write(nrite,"(12x,i10,4x,a4,4i10,10f15.6)")
     x      inv,keyword(:4),iatm1,iatm2,iatm3,iatm4,
     x      (parpot(j),j=1,mxpinv)
          
        endif
        
        keyinv(ninver)=keytmp
        lstinv(ninver,1)=iatm1
        lstinv(ninver,2)=iatm2
        lstinv(ninver,3)=iatm3
        lstinv(ninver,4)=iatm4
        prminv(ninver,1)=parpot(1)
        prminv(ninver,2)=parpot(2)

c     convert energies to internal units and angles to radians
        
        prminv(ninver,1)=prminv(ninver,1)*engunit
        
        if(keytmp.eq.2)then
          
          prminv(ninver,2)=cos(prminv(ninver,2)*(pi/180.d0))
          
        elseif(keytmp.eq.4)then
          
          prminv(ninver,2)=prminv(ninver,2)*engunit
          
        endif
        
      enddo
      
      return
      end subroutine define_inversions

      subroutine invfrc
     x  (lsolva,lfree,lexcite,idnode,imcon,mxnode,ntinv,enginv,virinv)

c***********************************************************************
c     
c     dl_poly subroutine for calculating inversion energy and force 
c     terms in molecular dynamics.
c     
c     copyright - daresbury laboratory 1996
c     author    - w. smith       may   1996
c     modified  - w.smith        jan   2011 : metadynamics
c     
c***********************************************************************
            
      implicit none

      logical safe,lsolva,lfree,lexcite,lselect,lcalcite
      logical idrive,jdrive,kdrive,ldrive
      integer idnode,imcon,mxnode,ntinv,fail1,fail2
      integer fail3,inv1,inv2,i,ii,ia,id,kk,ib,ic
      real(8) strs(6),strs_loc(6)
      real(8) xab,yab,zab,rab2,rrab,xac,yac,zac,rac2,rrac,xad,yad
      real(8) zad,rad2,rrad,rbc,rcd,rdb,ubx,uby,ubz,ubn,rub,vbx
      real(8) vby,vbz,rvb,wwb,ucx,ucy,ucz,ucn,vcx,vcy,vcz,rvc,wwc
      real(8) udx,udy,udz,udn,vdx,vdy,vdz,vdn,rvd,wwd,cosb,cosc
      real(8) cosd,thb,thc,thd,gamb,gamc,gamd,rubc,rubd,rucd,rucb
      real(8) rudb,rudc,rvbc,rvbd,rvcd,rvcb,rvdb,rvdc,fax,fay,faz
      real(8) fbx,fby,fbz,fcx,fcy,fcz,fdx,fdy,fdz,vbn,vcn,ruc,rud
      real(8) uuu,uu2,uun,uux,uuy,uuz,enginv,virinv,omega,gamma
      real(8) gamvir

      real(8), allocatable :: xdab(:),ydab(:),zdab(:)
      real(8), allocatable :: xdac(:),ydac(:),zdac(:)
      real(8), allocatable :: xdad(:),ydad(:),zdad(:)
      
      data fail1,fail2,fail3/0,0,0/

      allocate (xdab(msbad),ydab(msbad),zdab(msbad),stat=fail1)
      allocate (xdac(msbad),ydac(msbad),zdac(msbad),stat=fail2)
      allocate (xdad(msbad),ydad(msbad),zdad(msbad),stat=fail3)
      if(fail1.ne.0.or.fail2.ne.0.or.fail3.ne.0)
     x     call error(idnode,1130)

c     check size of work arrays
      
      if((ntinv-mxnode+1)/mxnode.gt.msbad)call error(idnode,427)

c     block indices
      
      inv1=(idnode*ntinv)/mxnode+1
      inv2=((idnode+1)*ntinv)/mxnode
      
      safe=.true.
      
c     initialise accumulators
      
      enginv=0.d0
      virinv=0.d0
      inv_fre=0.d0
      inv_vir=0.d0
      strs(:)=0.d0
      strs_loc(:)=0.d0

      if(lsolva)then
        
        lcomp(4)=.true.
        inv_sol(:)=0.d0
        if(lexcite)inv_exc(:)=0.d0
        
      endif
      
c     calculate bond vectors
      
      ii=0
      do i=inv1,inv2
        
        ii=ii+1
        
c     potential energy function type
        
        kk=listinv(ii,1)
        
c     indices of bonded atoms
        
        ia=listinv(ii,2)
        ib=listinv(ii,3)
        ic=listinv(ii,4)
        id=listinv(ii,5)
        
c     define components of bond vectors
        
        xdab(ii)=xxx(ib)-xxx(ia)
        ydab(ii)=yyy(ib)-yyy(ia)
        zdab(ii)=zzz(ib)-zzz(ia)
        
        if(keyinv(kk).eq.4)then

          xdac(ii)=xxx(ic)-xxx(ib)
          ydac(ii)=yyy(ic)-yyy(ib)
          zdac(ii)=zzz(ic)-zzz(ib)
          
          xdad(ii)=xxx(id)-xxx(ib)
          ydad(ii)=yyy(id)-yyy(ib)
          zdad(ii)=zzz(id)-zzz(ib)
          
        else
          
          xdac(ii)=xxx(ic)-xxx(ia)
          ydac(ii)=yyy(ic)-yyy(ia)
          zdac(ii)=zzz(ic)-zzz(ia)
          
          xdad(ii)=xxx(id)-xxx(ia)
          ydad(ii)=yyy(id)-yyy(ia)
          zdad(ii)=zzz(id)-zzz(ia)
          
        endif   
        
      enddo
      
c     periodic boundary condition
      
      call images(imcon,0,1,ii,cell,xdab,ydab,zdab)
      call images(imcon,0,1,ii,cell,xdac,ydac,zdac)
      call images(imcon,0,1,ii,cell,xdad,ydad,zdad)
      
c     loop over all specified inversions
      
      ii=0
      do i=inv1,inv2
        
        ii=ii+1
        
c     select potential energy function type
        
        kk=listinv(ii,1)
        lcalcite=(keyinv(kk).eq.4)
        
c     define components of bond vectors
        
        xab=xdab(ii)
        yab=ydab(ii)
        zab=zdab(ii)
        rab2=xab*xab+yab*yab+zab*zab
        rrab=1.d0/sqrt(rab2)
        
        xac=xdac(ii)
        yac=ydac(ii)
        zac=zdac(ii)
        rac2=xac*xac+yac*yac+zac*zac
        rrac=1.d0/sqrt(rac2)
        
        xad=xdad(ii)
        yad=ydad(ii)
        zad=zdad(ii)
        rad2=xad*xad+yad*yad+zad*zad
        rrad=1.d0/sqrt(rad2)
        
        if(lcalcite)then
          
c     calculate vector normal to plane
        
          uux=yac*zad-zac*yad
          uuy=zac*xad-xac*zad
          uuz=xac*yad-yac*xad
          uun=1.d0/sqrt(uux**2+uuy**2+uuz**2)
          uux=uun*uux
          uuy=uun*uuy
          uuz=uun*uuz
          uuu=xab*uux+yab*uuy+zab*uuz
          
        else
          
c     scalar products of bond vectors
        
          rbc=xab*xac+yab*yac+zab*zac
          rcd=xac*xad+yac*yad+zac*zad
          rdb=xad*xab+yad*yab+zad*zab
        
c     calculate bond-angle-plane vectors
        
          ubx=xac*rrac+xad*rrad
          uby=yac*rrac+yad*rrad
          ubz=zac*rrac+zad*rrad
          ubn=1.d0/sqrt(ubx**2+uby**2+ubz**2)
          ubx=ubn*ubx
          uby=ubn*uby
          ubz=ubn*ubz
          rub=xab*ubx+yab*uby+zab*ubz
          
          vbx=xac*rrac-xad*rrad
          vby=yac*rrac-yad*rrad
          vbz=zac*rrac-zad*rrad
          vbn=1.d0/sqrt(vbx**2+vby**2+vbz**2)
          vbx=vbn*vbx
          vby=vbn*vby
          vbz=vbn*vbz
          rvb=xab*vbx+yab*vby+zab*vbz
          wwb=sqrt(rub**2+rvb**2)
          
          ucx=xad*rrad+xab*rrab
          ucy=yad*rrad+yab*rrab
          ucz=zad*rrad+zab*rrab
          ucn=1.d0/sqrt(ucx**2+ucy**2+ucz**2)
          ucx=ucn*ucx
          ucy=ucn*ucy
          ucz=ucn*ucz
          ruc=xac*ucx+yac*ucy+zac*ucz
          
          vcx=xad*rrad-xab*rrab
          vcy=yad*rrad-yab*rrab
          vcz=zad*rrad-zab*rrab
          vcn=1.d0/sqrt(vcx**2+vcy**2+vcz**2)
          vcx=vcn*vcx
          vcy=vcn*vcy
          vcz=vcn*vcz
          rvc=xac*vcx+yac*vcy+zac*vcz
          wwc=sqrt(ruc**2+rvc**2)
          
          udx=xab*rrab+xac*rrac
          udy=yab*rrab+yac*rrac
          udz=zab*rrab+zac*rrac
          udn=1.d0/sqrt(udx**2+udy**2+udz**2)
          udx=udn*udx
          udy=udn*udy
          udz=udn*udz
          rud=xad*udx+yad*udy+zad*udz
          
          vdx=xab*rrab-xac*rrac
          vdy=yab*rrab-yac*rrac
          vdz=zab*rrab-zac*rrac
          vdn=1.d0/sqrt(vdx**2+vdy**2+vdz**2)
          vdx=vdn*vdx
          vdy=vdn*vdy
          vdz=vdn*vdz
          rvd=xad*vdx+yad*vdy+zad*vdz
          wwd=sqrt(rud**2+rvd**2)
          
c     calculate inversion angle cosines
          
          cosb=wwb*rrab
          cosc=wwc*rrac
          cosd=wwd*rrad
          
        endif
        
c     calculate potential energy and scalar force term
        
        if(keyinv(kk).eq.1)then
          
c     key=1 for harmonic inversion potential
          
          thb=acos(cosb)
          thc=acos(cosc)
          thd=acos(cosd)
          omega=0.5d0*prminv(kk,1)*((thb-prminv(kk,2))**2+
     x      (thc-prminv(kk,2))**2+(thd-prminv(kk,2))**2)/3.d0
          gamvir=0.d0
          
          gamb=0.d0
          if(abs(thb).gt.1.d-12)
     x      gamb=prminv(kk,1)*(thb-prminv(kk,2))/(3.d0*sin(thb))
          gamc=0.d0
          if(abs(thc).gt.1.d-12)
     x      gamc=prminv(kk,1)*(thc-prminv(kk,2))/(3.d0*sin(thc))
          gamd=0.d0
          if(abs(thd).gt.1.d-12)
     x      gamd=prminv(kk,1)*(thd-prminv(kk,2))/(3.d0*sin(thd))
          
        else if(keyinv(kk).eq.2)then
          
c     key=2 for harmonic cosine inversion potential
          
          omega=0.5d0*prminv(kk,1)*((cosb-prminv(kk,2))**2+
     x      (cosc-prminv(kk,2))**2+(cosb-prminv(kk,2))**2)/3.d0
          gamvir=0.d0

          gamb=-prminv(kk,1)*(cosb-prminv(kk,2))/3.d0
          gamc=-prminv(kk,1)*(cosc-prminv(kk,2))/3.d0
          gamd=-prminv(kk,1)*(cosd-prminv(kk,2))/3.d0
          
        else if(keyinv(kk).eq.3)then

c     key=3 for planar inversion potentials
          
          omega=prminv(kk,1)*((1.d0-cosb)+(1.d0-cosc)+(1.d0-cosd))/3.d0
          gamvir=0.d0
          
          gamb=prminv(kk,1)/3.d0
          gamc=prminv(kk,1)/3.d0
          gamd=prminv(kk,1)/3.d0

        else if(keyinv(kk).eq.4)then
          
c     key=4 for planar calcite potential
          
          uu2=uuu*uuu
          omega=uu2*(prminv(kk,1)+prminv(kk,2)*uu2)
          gamvir=uu2*(2.d0*prminv(kk,1)+4.d0*prminv(kk,2)*uu2)
          gamma=-uuu*(2.d0*prminv(kk,1)+4.d0*prminv(kk,2)*uu2)

        else

c     undefined potential
          
          safe=.false.
          gamvir=0.d0
          gamb=0.d0
          gamc=0.d0
          gamd=0.d0
          
        endif
        
c     indices of bonded atoms
        
        ia=listinv(ii,2)
        ib=listinv(ii,3)
        ic=listinv(ii,4)
        id=listinv(ii,5)

c     metadynamics local definitions
        
        if(lmetadyn)then
          
          idrive=driven(ltype(ia))
          jdrive=driven(ltype(ib))
          kdrive=driven(ltype(ic))
          ldrive=driven(ltype(id))
          
        endif
        
c     set selection control
        
        lselect=.true.
        
        if(lexcite)then
          
c     selected excitation option
        
          if((atm_fre(ia).ne.1).and.(atm_fre(ib).ne.1).and.
     x      (atm_fre(ic).ne.1).and.(atm_fre(id).ne.1))then
            
c     reset selection control
            
            lselect=(atm_fre(ia)+atm_fre(ib)+atm_fre(ic)+
     x        atm_fre(id).eq.0)
            
            if(lsolva)then
              inv_exc(atmolt(ia))=inv_exc(atmolt(ia))+omega
            endif
            
          endif
          
        elseif(lfree)then
          
c     selected free energy option
          
          if((atm_fre(ia).eq.1).or.(atm_fre(ib).eq.1).or.
     x      (atm_fre(ic).eq.1).or.(atm_fre(id).eq.1))then
            
c     set hamiltonian mixing parameter

            inv_fre=inv_fre-omega
            omega=lambda1*omega
            gamb=lambda1*gamb
            gamc=lambda1*gamc
            gamd=lambda1*gamd
            
          elseif((atm_fre(ia).eq.2).or.(atm_fre(ib).eq.2).or.
     x        (atm_fre(ic).eq.2).or.(atm_fre(id).eq.2))then
            
c     set hamiltonian mixing parameter

            inv_fre=inv_fre+omega
            omega=lambda2*omega
            gamb=lambda2*gamb
            gamc=lambda2*gamc
            gamd=lambda2*gamd
                        
          endif
          
        endif
        
c     calculate bond and u,v scalar products
        
        if(.not.lcalcite)then
          
          rubc=xab*ucx+yab*ucy+zab*ucz
          rubd=xab*udx+yab*udy+zab*udz
          rucd=xac*udx+yac*udy+zac*udz
          rucb=xac*ubx+yac*uby+zac*ubz
          rudb=xad*ubx+yad*uby+zad*ubz
          rudc=xad*ucx+yad*ucy+zad*ucz
          
          rvbc=xab*vcx+yab*vcy+zab*vcz
          rvbd=xab*vdx+yab*vdy+zab*vdz
          rvcd=xac*vdx+yac*vdy+zac*vdz
          rvcb=xac*vbx+yac*vby+zac*vbz
          rvdb=xad*vbx+yad*vby+zad*vbz
          rvdc=xad*vcx+yad*vcy+zad*vcz
          
        endif
        
        if(lselect)then
          
c     calculate potential energy and virial
          
          enginv=enginv+omega
          virinv=virinv+gamvir
          
c     calculate solvation energy
          
          if(lsolva)then
            inv_sol(atmolt(ia))=inv_sol(atmolt(ia))+omega
          endif
          
c     calculate atomic forces
          
          if(lcalcite)then
            
            fax=-gamma*uux
            fay=-gamma*uuy
            faz=-gamma*uuz
            
            fcx=gamma*uun*((yad*zab-zad*yab)-uuu*(yad*uuz-zad*uuy))
            fcy=gamma*uun*((zad*xab-xad*zab)-uuu*(zad*uux-xad*uuz))
            fcz=gamma*uun*((xad*yab-yad*xab)-uuu*(xad*uuy-yad*uux))
            
            fdx=gamma*uun*((yab*zac-zab*yac)-uuu*(zac*uuy-yac*uuz))
            fdy=gamma*uun*((zab*xac-xab*zac)-uuu*(xac*uuz-zac*uux))
            fdz=gamma*uun*((xab*yac-yab*xac)-uuu*(yac*uux-xac*uuy))
            
            fbx=-(fax+fcx+fdx)
            fby=-(fay+fcy+fdy)
            fbz=-(faz+fcz+fdz)
            
c     stress tensor calculation for calcite terms
            
            strs(1)=strs(1)+uuu*gamma*uux*uux
            strs(2)=strs(2)+uuu*gamma*uux*uuy
            strs(3)=strs(3)+uuu*gamma*uux*uuz
            strs(4)=strs(4)+uuu*gamma*uuy*uuy
            strs(5)=strs(5)+uuu*gamma*uuy*uuz
            strs(6)=strs(6)+uuu*gamma*uuz*uuz
          
          else
            
            fbx=gamb*(-cosb*xab*rrab**2+rrab*(rub*ubx+rvb*vbx)/wwb)
     x      +(ruc*ucn*rrab*(xac-ruc*ucx-(rbc-ruc*rubc)*xab*rrab**2)
     x      -rvc*vcn*rrab*(xac-rvc*vcx-(rbc-rvc*rvbc)*xab*rrab**2))
     x      *gamc*rrac/wwc
     x      +(rud*udn*rrab*(xad-rud*udx-(rdb-rud*rubd)*xab*rrab**2)
     x      +rvd*vdn*rrab*(xad-rvd*vdx-(rdb-rvd*rvbd)*xab*rrab**2))
     x      *gamd*rrad/wwd
          
            fby=gamb*(-cosb*yab*rrab**2+rrab*(rub*uby+rvb*vby)/wwb)
     x      +(ruc*ucn*rrab*(yac-ruc*ucy-(rbc-ruc*rubc)*yab*rrab**2)
     x      -rvc*vcn*rrab*(yac-rvc*vcy-(rbc-rvc*rvbc)*yab*rrab**2))
     x      *gamc*rrac/wwc
     x      +(rud*udn*rrab*(yad-rud*udy-(rdb-rud*rubd)*yab*rrab**2)
     x      +rvd*vdn*rrab*(yad-rvd*vdy-(rdb-rvd*rvbd)*yab*rrab**2))
     x      *gamd*rrad/wwd
          
            fbz=gamb*(-cosb*zab*rrab**2+rrab*(rub*ubz+rvb*vbz)/wwb)
     x      +(ruc*ucn*rrab*(zac-ruc*ucz-(rbc-ruc*rubc)*zab*rrab**2)
     x      -rvc*vcn*rrab*(zac-rvc*vcz-(rbc-rvc*rvbc)*zab*rrab**2))
     x      *gamc*rrac/wwc
     x      +(rud*udn*rrab*(zad-rud*udz-(rdb-rud*rubd)*zab*rrab**2)
     x      +rvd*vdn*rrab*(zad-rvd*vdz-(rdb-rvd*rvbd)*zab*rrab**2))
     x      *gamd*rrad/wwd
          
            fcx=gamc*(-cosc*xac*rrac**2+rrac*(ruc*ucx+rvc*vcx)/wwc)
     x      +(rud*udn*rrac*(xad-rud*udx-(rcd-rud*rucd)*xac*rrac**2)
     x      -rvd*vdn*rrac*(xad-rvd*vdx-(rcd-rvd*rvcd)*xac*rrac**2))
     x      *gamd*rrad/wwd
     x      +(rub*ubn*rrac*(xab-rub*ubx-(rbc-rub*rucb)*xac*rrac**2)
     x      +rvb*vbn*rrac*(xab-rvb*vbx-(rbc-rvb*rvcb)*xac*rrac**2))
     x      *gamb*rrab/wwb
          
            fcy=gamc*(-cosc*yac*rrac**2+rrac*(ruc*ucy+rvc*vcy)/wwc)
     x      +(rud*udn*rrac*(yad-rud*udy-(rcd-rud*rucd)*yac*rrac**2)
     x      -rvd*vdn*rrac*(yad-rvd*vdy-(rcd-rvd*rvcd)*yac*rrac**2))
     x      *gamd*rrad/wwd
     x      +(rub*ubn*rrac*(yab-rub*uby-(rbc-rub*rucb)*yac*rrac**2)
     x      +rvb*vbn*rrac*(yab-rvb*vby-(rbc-rvb*rvcb)*yac*rrac**2))
     x      *gamb*rrab/wwb
          
            fcz=gamc*(-cosc*zac*rrac**2+rrac*(ruc*ucz+rvc*vcz)/wwc)
     x      +(rud*udn*rrac*(zad-rud*udz-(rcd-rud*rucd)*zac*rrac**2)
     x      -rvd*vdn*rrac*(zad-rvd*vdz-(rcd-rvd*rvcd)*zac*rrac**2))
     x      *gamd*rrad/wwd
     x      +(rub*ubn*rrac*(zab-rub*ubz-(rbc-rub*rucb)*zac*rrac**2)
     x      +rvb*vbn*rrac*(zab-rvb*vbz-(rbc-rvb*rvcb)*zac*rrac**2))
     x      *gamb*rrab/wwb
          
            fdx=gamd*(-cosd*xad*rrad**2+rrad*(rud*udx+rvd*vdx)/wwd)
     x      +(rub*ubn*rrad*(xab-rub*ubx-(rdb-rub*rudb)*xad*rrad**2)
     x      -rvb*vbn*rrad*(xab-rvb*vbx-(rdb-rvb*rvdb)*xad*rrad**2))
     x      *gamb*rrab/wwb
     x      +(ruc*ucn*rrad*(xac-ruc*ucx-(rcd-ruc*rudc)*xad*rrad**2)
     x      +rvc*vcn*rrad*(xac-rvc*vcx-(rcd-rvc*rvdc)*xad*rrad**2))
     x      *gamc*rrac/wwc
          
            fdy=gamd*(-cosd*yad*rrad**2+rrad*(rud*udy+rvd*vdy)/wwd)
     x      +(rub*ubn*rrad*(yab-rub*uby-(rdb-rub*rudb)*yad*rrad**2)
     x      -rvb*vbn*rrad*(yab-rvb*vby-(rdb-rvb*rvdb)*yad*rrad**2))
     x      *gamb*rrab/wwb
     x      +(ruc*ucn*rrad*(yac-ruc*ucy-(rcd-ruc*rudc)*yad*rrad**2)
     x      +rvc*vcn*rrad*(yac-rvc*vcy-(rcd-rvc*rvdc)*yad*rrad**2))
     x      *gamc*rrac/wwc
          
            fdz=gamd*(-cosd*zad*rrad**2+rrad*(rud*udz+rvd*vdz)/wwd)
     x      +(rub*ubn*rrad*(zab-rub*ubz-(rdb-rub*rudb)*zad*rrad**2)
     x      -rvb*vbn*rrad*(zab-rvb*vbz-(rdb-rvb*rvdb)*zad*rrad**2))
     x      *gamb*rrab/wwb
     x      +(ruc*ucn*rrad*(zac-ruc*ucz-(rcd-ruc*rudc)*zad*rrad**2)
     x      +rvc*vcn*rrad*(zac-rvc*vcz-(rcd-rvc*rvdc)*zad*rrad**2))
     x      *gamc*rrac/wwc
          
            fax=-(fbx+fcx+fdx)
            fay=-(fby+fcy+fdy)
            faz=-(fbz+fcz+fdz)
            
c     stress tensor calculation for inversion terms
            
            strs(1)=strs(1)+xab*fbx+xac*fcx+xad*fdx 
            strs(2)=strs(2)+yab*fbx+yac*fcx+yad*fdx 
            strs(3)=strs(3)+zab*fbx+zac*fcx+zad*fdx 
            strs(4)=strs(4)+yab*fby+yac*fcy+yad*fdy 
            strs(5)=strs(5)+yab*fbz+yac*fcz+yad*fdz 
            strs(6)=strs(6)+zab*fbz+zac*fcz+zad*fdz 
          
          endif
          
          fxx(ia)=fxx(ia)+fax
          fyy(ia)=fyy(ia)+fay
          fzz(ia)=fzz(ia)+faz
          
          fxx(ib)=fxx(ib)+fbx
          fyy(ib)=fyy(ib)+fby
          fzz(ib)=fzz(ib)+fbz
        
          fxx(ic)=fxx(ic)+fcx
          fyy(ic)=fyy(ic)+fcy
          fzz(ic)=fzz(ic)+fcz
          
          fxx(id)=fxx(id)+fdx
          fyy(id)=fyy(id)+fdy
          fzz(id)=fzz(id)+fdz
          
        endif
        
c     metadynamics local parameters
        
        if(lmetadyn.and.(idrive.or.jdrive.or.kdrive.or.ldrive))then
          
c     local energy and virial
          
          eng_loc=eng_loc+omega
          vir_loc=vir_loc+gamvir
          
c     local forces
          
          fxx_loc(ia)=fxx_loc(ia)+fax
          fyy_loc(ia)=fyy_loc(ia)+fay
          fzz_loc(ia)=fzz_loc(ia)+faz
          
          fxx_loc(ib)=fxx_loc(ib)+fbx
          fyy_loc(ib)=fyy_loc(ib)+fby
          fzz_loc(ib)=fzz_loc(ib)+fbz
          
          fxx_loc(ic)=fxx_loc(ic)+fcx
          fyy_loc(ic)=fyy_loc(ic)+fcy
          fzz_loc(ic)=fzz_loc(ic)+fcz
          
          fxx_loc(id)=fxx_loc(id)+fdx
          fyy_loc(id)=fyy_loc(id)+fdy
          fzz_loc(id)=fzz_loc(id)+fdz
          
c     local stress tensor
          
          strs_loc(1)=strs_loc(1)+xab*fbx+xac*fcx+xad*fdx 
          strs_loc(2)=strs_loc(2)+yab*fbx+yac*fcx+yad*fdx 
          strs_loc(3)=strs_loc(3)+zab*fbx+zac*fcx+zad*fdx 
          strs_loc(4)=strs_loc(4)+yab*fby+yac*fcy+yad*fdy 
          strs_loc(5)=strs_loc(5)+yab*fbz+yac*fcz+yad*fdz 
          strs_loc(6)=strs_loc(6)+zab*fbz+zac*fcz+zad*fdz 
          
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
      if(.not.safe)call error(idnode,449)
      
c     sum contributions over all nodes
      
      if(mxnode.gt.1)then
        
        buffer(1)=enginv
        buffer(2)=virinv
        buffer(3)=inv_fre
        buffer(4)=inv_vir
        call gdsum(buffer(1),4,buffer(5))
        enginv=buffer(1)
        virinv=buffer(2)
        inv_fre=buffer(3)
        inv_vir=buffer(4)
        
        if(lsolva)then
          
          call gdsum(inv_sol(1),mxtmls,buffer(1))
          if(lexcite)call gdsum(inv_exc(1),mxtmls,buffer(1))
          
        endif

      endif
      
      deallocate (xdab,ydab,zdab,stat=fail1)
      deallocate (xdac,ydac,zdac,stat=fail2)
      deallocate (xdad,ydad,zdad,stat=fail3)
      
      return
      end subroutine invfrc
      
      end module inversion_module
