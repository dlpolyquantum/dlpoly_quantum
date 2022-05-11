      module four_body_module

c***********************************************************************
c     
c     dl_poly module for defining four-body potential arrays
c     copyright - daresbury laboratory
c     author  - w. smith  sep 2003
c     adapted - w. smith  aug 2008 : solvation, free energy, excitation 
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
      use utility_module

      implicit none

      logical, allocatable :: filter(:)
      real(8), allocatable :: prmfbp(:,:),rcut4b(:)
      integer, allocatable :: lstfbp(:),ltpfbp(:),latfbp(:)

      save prmfbp,rcut4b,lstfbp,ltpfbp,latfbp,filter

      contains
      
      subroutine alloc_fbp_arrays(idnode,mxnode)

      implicit none

      logical safe
      integer i,fail,idnode,mxnode
      dimension fail(6)
      
      safe=.true.
      
c     allocate arrays
      
      fail(:)=0
      
      allocate (prmfbp(mxfbp,mxpfbp),stat=fail(1))
      allocate (rcut4b(mxfbp),stat=fail(2))
      allocate (lstfbp(mxfbp),stat=fail(3))
      allocate (ltpfbp(mxfbp),stat=fail(4))
      allocate (latfbp(mxatms),stat=fail(5))
      allocate (filter(mxsite),stat=fail(6))
      
      if(any(fail.gt.0))safe=.false.      
      if(mxnode.gt.1)call gstate(safe)    
      if(.not.safe)call error(idnode,1140)
      
      end subroutine alloc_fbp_arrays

      subroutine define_four_body
     x  (safe,lunits,lmols,idnode,ntpfbp,ntpatm,
     x  rcutfb,engunit)

c***********************************************************************
c     
c     dl_poly subroutine for defining four body potentials
c     
c     copyright - daresbury laboratory
c     author    - w. smith    aug 2003
c     
c***********************************************************************

      implicit none

      logical safe,lunits,lmols
      character*8 keyword
      character*8 atom0,atom1,atom2,atom3
      character*1 message(80)
      integer idnode,ntpfbp,ntpatm,ifbp,itpfbp,keypot,katm0
      integer i,katm1,katm2,katm3,jtpatm,ka1,ka2,ka3,keyfbp,kfbp
      integer j,fail,idum
      real(8) rcutfb,engunit
      real(8), allocatable :: parpot(:)

      data fail/0/

      allocate (parpot(mxpfbp),stat=fail)
      if(fail.ne.0)call error(idnode,1150)

      ntpfbp=intstr(record,lenrec,idum)
      
      if(idnode.eq.0) then
        
        write(nrite,"(/,/,1x,'number of specified four ',
     x    'body potentials',i10)") ntpfbp
        write(nrite,"(/,/,16x,'atom 1  ','atom 2  ','atom 3  ',
     x    'atom 4  ',3x,' key',30x,'parameters'/,/)")
        
      endif      
      if(ntpfbp.gt.mxfbp) call error(idnode,89)
      if(.not.lunits) call error(idnode,6)
      if(.not.lmols) call error(idnode,13)
      
      do i=1,mxsite
        filter(i)=.false.
      enddo

      do ifbp=1,mxfbp
        lstfbp(ifbp)=0
      enddo
      
      do ifbp=1,mxfbp,mx3fbp
        lstfbp(ifbp)=-1
      enddo
      
      rcutfb=0.d0
      
      do itpfbp=1,ntpfbp
        
        do i=1,mxpfbp
          parpot(i)=0.d0
        enddo
        
        call getrec(safe,idnode,nfield)
        if(.not.safe)return

c     Note the order!! atom0 is the central atom

        call copystring(record,message,80)
        call getword(atom0,record,8,lenrec)
        call getword(atom1,record,8,lenrec)
        call getword(atom2,record,8,lenrec)
        call getword(atom3,record,8,lenrec)
        call lowcase(record,4)
        call getword(keyword,record,4,lenrec)

        if(keyword(1:4).eq.'harm') then
          keypot=1
        elseif(keyword(1:4).eq.'hcos') then
          keypot=2
        elseif(keyword(1:4).eq.'plan') then
          keypot=3
        else
          if(idnode.eq.0) write(nrite,*) message
          call error(idnode,443)
        endif

        parpot(1)=dblstr(record,lenrec,idum)
        parpot(2)=dblstr(record,lenrec,idum)
        parpot(3)=dblstr(record,lenrec,idum)
        
        if(idnode.eq.0) 
     x    write(nrite,"(16x,4a8,4x,a4,1x,1p,9e13.5)") 
     x    atom0,atom1,atom2,atom3,keyword(1:4),(parpot(j),j=1,mxpfbp)
        
        katm0=0
        katm1=0
        katm2=0
        katm3=0
        
        do jtpatm=1,ntpatm
          
          if(atom0.eq.unqatm(jtpatm))katm0=jtpatm
          if(atom1.eq.unqatm(jtpatm))katm1=jtpatm
          if(atom2.eq.unqatm(jtpatm))katm2=jtpatm
          if(atom3.eq.unqatm(jtpatm))katm3=jtpatm
          
        enddo
        
        if(katm0.eq.0.or.katm1.eq.0.or.katm2.eq.0.or.
     x    katm3.eq.0) call error(idnode,91)
        
        filter(katm0)=.true.
        filter(katm1)=.true.
        filter(katm2)=.true.
        filter(katm3)=.true.
        
        ka1=max(katm1,katm2,katm3)
        ka3=min(katm1,katm2,katm3)
        ka2=katm1+katm2+katm3-ka1-ka3
        keyfbp=ka3+(ka2*(ka2-1))/2+(ka1*(ka1**2-1))/6+
     x    (katm0-1)*mx3fbp

        if(keyfbp.gt.mxfbp) call error(idnode,101)

c     convert parameters to internal units
        
        parpot(1)=parpot(1)*engunit
        parpot(2)=parpot(2)*(pi/180.d0)

        if(keypot.eq.2)then

          parpot(2)=cos(parpot(2))

        endif

        if(lstfbp(keyfbp).gt.0) call error(idnode,19)
        lstfbp(keyfbp)=itpfbp
        ltpfbp(itpfbp)=keypot
        kfbp=mx3fbp*((keyfbp-1)/mx3fbp)+1
        if(lstfbp(kfbp).lt.0)lstfbp(kfbp)=0

c     calculate max four body cutoff
        
        rcutfb=max(rcutfb,parpot(3))
        rcut4b(itpfbp)=parpot(3)

c     store four body potential parameters
        
        do i=1,mxpfbp
          prmfbp(itpfbp,i)=parpot(i)
        enddo

      enddo

      if(rcutfb.lt.1.d-6)call error(idnode,453)
      
      deallocate (parpot,stat=fail)

      return
      end subroutine define_four_body

      subroutine fbpfrc
     x  (lsolva,lfree,lexcite,idnode,mxnode,natms,imcon,rcutfb,
     x  engfbp,virfbp)
      
c***********************************************************************
c     
c     dl_poly subroutine for calculating four body inversion forces
c     arising from the inversion angle between three atoms around a
c     nominated central atom
c     
c     Note: the subroutine converts coordinates to reduced units
c     to avoid a call to images.f. The link cell algorithm used
c     here necessitates a parallelepiped cell geometry
c     
c     copyright - daresbury laboratory 1996
c     author   - w.smith july 1996
c     adapted   - w.smith aug 2008 solvation, free energy etc
c     
c***********************************************************************
      
      implicit none
      
      logical safe,lsolva,lfree,lexcite,lselect,lskip
      logical idrive,jdrive,kdrive,ldrive
      integer idnode,mxnode,natms,imcon,nix,niy,niz
      integer i,j,k,nbx,nby,nbz,ncells,ix,iy,iz,icell,jx,jy
      integer jz,jj,kk,ia,ib,ifbp,jfbp,kfbp,jklbd,kkfbp,ktyp,ii
      integer ic,ll,id,lfbp,l,limit,jcell,kkk
      real(8) rcutfb,engfbp,virfbp,vbn,vcn,pterm,xm,ym,zm
      real(8) strs,cprp,det,xdc,ydc,zdc,sxx,syy,szz,sxab,strs_loc
      real(8) syab,szab,xab,yab,zab,rab2,sxac,syac,szac,xac,yac
      real(8) zac,rac2,sxad,syad,szad,xad,yad,zad,rad2,rrab,rrac
      real(8) rrad,rbc,rcd,rdb,ubx,uby,ubz,ubn,rub,vbx,vby,vbz
      real(8) rvb,wwb,ucx,ucy,ucz,ucn,ruc,vcx,vcy,vcz,rvc,wwc
      real(8) udx,udy,udz,udn,rud,vdx,vdy,vdz,vdn,rvd,wwd,cosb
      real(8) cosc,cosd,thb,thc,thd,gamb,gamc,gamd,rubc,rubd
      real(8) rucd,rucb,rudb,rudc,rvbc,rvbd,rvcd,rvcb,rvdb,rvdc
      real(8) fax,fay,faz,fbx,fby,fbz,fcx,fcy,fcz,fdx,fdy,fdz
      dimension cprp(10),strs(6),nix(27),niy(27),niz(27),strs_loc(6)
      
      data nix/ 0,-1,-1,-1, 0, 0,-1, 1,-1, 0, 1,-1, 0, 1,
     x  1, 1, 1, 0, 0, 1,-1, 1, 0,-1, 1, 0,-1/
      data niy/ 0, 0,-1, 1, 1, 0, 0, 0,-1,-1,-1, 1, 1, 1,
     x  0, 1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1/
      data niz/ 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     x  0, 0, 0, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1/
      
      lskip=(lfree.or.lexcite)
      
c     flag for undefined potentials
      
      safe=.true.
      
c     initialise accumulators
      
      engfbp=0.d0
      virfbp=0.d0
      fbp_fre=0.d0
      fbp_vir=0.d0
      strs(:)=0.d0
      strs_loc(:)=0.d0
      
      if(lsolva)then
        
        lcomp(9)=.true.
        en4_sol(:)=0.d0
        if(lexcite)en4_exc(:)=0.d0
        
      endif
      
c     create mock cell vectors for non-periodic system
      
      if(imcon.eq.0) then
        
        xm=0.d0
        ym=0.d0
        zm=0.d0
        do i=1,natms
          xm=max(xm,abs(xxx(i)))
          ym=max(ym,abs(yyy(i)))
          zm=max(zm,abs(zzz(i)))
        enddo
        
        cell(1)=2.d0*xm+rcutfb
        cell(2)=0.d0
        cell(3)=0.d0
        cell(4)=0.d0
        cell(5)=2.d0*ym+rcutfb
        cell(6)=0.d0
        cell(7)=0.d0
        cell(8)=0.d0
        cell(9)=2.d0*zm+rcutfb
      
      endif
      
c     check for appropriate boundary conditions
      
      if(imcon.gt.3)call error(idnode,79)
      call invert(cell,rcell,det)
      call dcell(cell,cprp)
      
c     calculate link cell numbers
      
      nbx=int(cprp(7)/(rcutfb+1.d-6))
      nby=int(cprp(8)/(rcutfb+1.d-6))
      nbz=int(cprp(9)/(rcutfb+1.d-6))
      ncells=nbx*nby*nbz
      if(ncells.gt.mxcell) then
        
        if(idnode.eq.0) write(nrite,'(a,i6)')
     x    'number of required link cells in routine fbpfrc is ',ncells
        write(nrite,'(a,i6)')
     x    'number of default link cells in routine fbpfrc is ',mxcell
        call error(idnode,87)
        
      endif
      
c     transform atomic coordinates and construct link cells
      
      do l=1,ncells
        
        lct(l)=0
        lst(l)=0
        
      enddo
      
      xdc=dble(nbx)
      ydc=dble(nby)
      zdc=dble(nbz)
      
      do i=1,natms
        
        if(filter(ltype(i)))then
          
          sxx=rcell(1)*xxx(i)+rcell(4)*yyy(i)+rcell(7)*zzz(i)
          syy=rcell(2)*xxx(i)+rcell(5)*yyy(i)+rcell(8)*zzz(i)
          szz=rcell(3)*xxx(i)+rcell(6)*yyy(i)+rcell(9)*zzz(i)
          
          xxx(i)=sxx
          yyy(i)=syy
          zzz(i)=szz
          
          ix=min(int(xdc*(sxx+0.5d0)),nbx-1)
          iy=min(int(ydc*(syy+0.5d0)),nby-1)
          iz=min(int(zdc*(szz+0.5d0)),nbz-1)
          k=1+ix+nbx*(iy+nby*iz)
          lst(k)=lst(k)+1
          link(i)=lct(k)
          lct(k)=i
          
        endif
        
      enddo
      
c     loop over central atoms of inversion
      
      ix=0
      iy=1
      iz=1
      do icell=1,ncells
        
        ix=ix+1
        if(ix.gt.nbx)then
          ix=1
          iy=iy+1
          if(iy.gt.nby)then
            iy=1
            iz=iz+1
          endif
        endif
        
        k=0
        do kk=1,27
          
          jx=ix+nix(kk)
          jy=iy+niy(kk)
          jz=iz+niz(kk)
          
          if(jx.gt.nbx)jx=1
          if(jy.gt.nby)jy=1
          if(jz.gt.nbz)jz=1
          if(jx.lt.1)jx=jx+nbx
          if(jy.lt.1)jy=jy+nby
          if(jz.lt.1)jz=jz+nbz
          
          jcell=jx+nbx*(jy-1+nby*(jz-1))
          j=lct(jcell)
          
          do ii=1,lst(jcell)
            
            k=k+1
            latfbp(k)=j
            j=link(j)
            
          enddo
          
        enddo
        
        limit=k
        
        do ii=1,lst(icell)
          
          ia=latfbp(ii)
          if(lmetadyn)idrive=driven(ltype(ia))
          ifbp=mx3fbp*(ltype(ia)-1)
          if(mod(ia,mxnode).eq.idnode.and.lstfbp(ifbp+1).ge.0)then
            
          do jj=1,limit-2
          
          ib=latfbp(jj)
          if(lmetadyn)jdrive=driven(ltype(ib))
          
          do kk=jj+1,limit-1
          
          ic=latfbp(kk)
          if(lmetadyn)kdrive=driven(ltype(ic))
                
          do ll=kk+1,limit
          
          id=latfbp(ll)
          if(lmetadyn)ldrive=driven(ltype(id))
                  
          if(lskip)then
            
            if((atm_fre(ia).eq.1.or.atm_fre(ib).eq.1.or.
     x        atm_fre(ic).eq.1.or.atm_fre(id).eq.1).and.
     x        (atm_fre(ia).eq.2.or.atm_fre(ib).eq.2.or.
     x        atm_fre(ic).eq.2.or.atm_fre(id).eq.2))cycle
            
          endif
          
          jfbp=max(ltype(ib),ltype(ic),ltype(id))
          lfbp=min(ltype(ib),ltype(ic),ltype(id))
          kfbp=ltype(ib)+ltype(ic)+ltype(id)-jfbp-lfbp
          jklbd=ifbp+lfbp+(kfbp*(kfbp-1))/2+(jfbp*(jfbp**2-1))/6
          kkfbp=lstfbp(jklbd)
          if(kkfbp.gt.0)then
                    
          sxab=xxx(ib)-xxx(ia)
          sxab=sxab-nint(sxab)
          syab=yyy(ib)-yyy(ia)
          syab=syab-nint(syab)
          szab=zzz(ib)-zzz(ia)
          szab=szab-nint(szab)
          
          xab=cell(1)*sxab+cell(4)*syab+cell(7)*szab
          if(abs(xab).lt.rcutfb)then
          
          yab=cell(2)*sxab+cell(5)*syab+cell(8)*szab
          if(abs(yab).lt.rcutfb)then
          
          zab=cell(3)*sxab+cell(6)*syab+cell(9)*szab
          if(abs(zab).lt.rcutfb)then
          
          rab2=xab*xab+yab*yab+zab*zab
          
          sxac=xxx(ic)-xxx(ia)
          sxac=sxac-nint(sxac)
          syac=yyy(ic)-yyy(ia)
          syac=syac-nint(syac)
          szac=zzz(ic)-zzz(ia)
          szac=szac-nint(szac)
          
          xac=cell(1)*sxac+cell(4)*syac+cell(7)*szac
          if(abs(xac).lt.rcutfb)then
                            
          yac=cell(2)*sxac+cell(5)*syac+cell(8)*szac
          if(abs(yac).lt.rcutfb)then
          
          zac=cell(3)*sxac+cell(6)*syac+cell(9)*szac
          if(abs(zac).lt.rcutfb)then
          
          rac2=xac*xac+yac*yac+zac*zac
          
          sxad=xxx(id)-xxx(ia)
          sxad=sxad-nint(sxad)
          syad=yyy(id)-yyy(ia)
          syad=syad-nint(syad)
          szad=zzz(id)-zzz(ia)
          szad=szad-nint(szad)
          
          xad=cell(1)*sxad+cell(4)*syad+cell(7)*szad
          if(abs(xad).lt.rcutfb)then
          
          yad=cell(2)*sxad+cell(5)*syad+cell(8)*szad
          if(abs(yad).lt.rcutfb)then
          
          zad=cell(3)*sxad+cell(6)*syad+cell(9)*szad
          if(abs(zad).lt.rcutfb)then
          
          rad2=xad*xad+yad*yad+zad*zad
          
          if(rcut4b(kkfbp)**2.ge.max(rab2,rac2,rad2))then
          
          rrab=1.d0/sqrt(rab2)
          rrac=1.d0/sqrt(rac2)
          rrad=1.d0/sqrt(rad2)
          
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
          if(abs(cosb).gt.1.d0)cosb=sign(1.d0,cosb)
          if(abs(cosc).gt.1.d0)cosc=sign(1.d0,cosc)
          if(abs(cosd).gt.1.d0)cosd=sign(1.d0,cosd)
          
c     select potential energy function type
          
          ktyp=ltpfbp(kkfbp)
          
c     calculate potential energy and scalar force term
          
          if(ktyp.eq.1)then
            
c     key=1 for harmonic potential
            
            thb=acos(cosb)
            thc=acos(cosc)
            thd=acos(cosd)
            
            pterm=0.5d0*prmfbp(kkfbp,1)*
     x        ((thb-prmfbp(kkfbp,2))**2+
     x        (thc-prmfbp(kkfbp,2))**2+
     x        (thd-prmfbp(kkfbp,2))**2)/3.d0
            
            gamb=0.d0
            if(abs(thb).gt.1.d-12)gamb=prmfbp(kkfbp,1)*
     x        (thb-prmfbp(kkfbp,2))/(3.d0*sin(thb))
            gamc=0.d0
            if(abs(thc).gt.1.d-12)gamc=prmfbp(kkfbp,1)*
     x        (thc-prmfbp(kkfbp,2))/(3.d0*sin(thc))
            gamd=0.d0
            if(abs(thd).gt.1.d-12)gamd=prmfbp(kkfbp,1)*
     x        (thd-prmfbp(kkfbp,2))/(3.d0*sin(thd))
            
          else if(ktyp.eq.2)then
            
c     key=2 for harmonic cosine inversion potential
            
            pterm=0.5d0*prmfbp(kkfbp,1)*
     x        ((cosb-prmfbp(kkfbp,2))**2+
     x        (cosc-prmfbp(kkfbp,2))**2+
     x        (cosd-prmfbp(kkfbp,2))**2)/3.d0
            
            gamb=-prmfbp(kkfbp,1)*(cosb-prmfbp(kkfbp,2))/3.d0
            gamc=-prmfbp(kkfbp,1)*(cosc-prmfbp(kkfbp,2))/3.d0
            gamd=-prmfbp(kkfbp,1)*(cosd-prmfbp(kkfbp,2))/3.d0
            
          else if(ktyp.eq.3)then
            
c     key=3 for planar inversion potentials
            
            pterm=prmfbp(kkfbp,1)*
     x        (3.d0-cosb-cosc-cosd)/3.d0
            
            gamb=-prmfbp(kkfbp,1)/3.d0
            gamc=-prmfbp(kkfbp,1)/3.d0
            gamd=-prmfbp(kkfbp,1)/3.d0
            
          else
            
c     undefined potential
            
            safe=.false.
            pterm=0.d0
            gamb=0.d0
            gamc=0.d0
            gamd=0.d0
            
          endif
          
c     set selection control
          
          lselect=.true.
          
c     set quadruple index
          
          if(lsolva)
     x      kkk=loc4(atmolt(ia),atmolt(ib),atmolt(ic),atmolt(id))
          
          if(lexcite)then
            
c     selected excitation option
            
            if((atm_fre(ia).ne.1).and.(atm_fre(ib).ne.1).and.
     x        (atm_fre(ic).ne.1).and.(atm_fre(id).ne.1))then
              
c     reset selection control
              
              lselect=(atm_fre(ia)+atm_fre(ib)+atm_fre(ic)+
     x          atm_fre(id).eq.0)
              
              if(lsolva)en4_exc(kkk)=en4_exc(kkk)+pterm
              
            endif
            
          elseif(lfree)then
            
c     selected free energy option
            
            if((atm_fre(ia).eq.1).or.(atm_fre(ib).eq.1).or.
     x        (atm_fre(ic).eq.1).or.(atm_fre(id).eq.1))then
              
c     set hamiltonian mixing parameter
              
              fbp_fre=fbp_fre-pterm
              pterm=lambda1*pterm
              gamb=lambda1*gamb
              gamc=lambda1*gamc
              gamd=lambda1*gamd
              
            elseif((atm_fre(ia).eq.2).or.(atm_fre(ib).eq.2).or.
     x          (atm_fre(ic).eq.2).or.(atm_fre(id).eq.2))then
              
c     set hamiltonian mixing parameter
              
              fbp_fre=fbp_fre+pterm
              pterm=lambda2*pterm
              gamb=lambda2*gamb
              gamc=lambda2*gamc
              gamd=lambda2*gamd
              
            endif
            
          endif
          
          if(lselect)then
            
c     calculate potential
            
            engfbp=engfbp+pterm
            
            if(lsolva)en4_sol(kkk)=en4_sol(kkk)+pterm
            
c     calculate bond and u,v scalar products
            
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
            
c     calculate atomic forces
            
            fbx=gamb*(-cosb*xab*rrab**2+rrab*(rub*ubx+rvb*vbx)/wwb)
     x        +(ruc*ucn*rrab*(xac-ruc*ucx-(rbc-ruc*rubc)*xab*rrab**2)
     x        - rvc*vcn*rrab*(xac-rvc*vcx-(rbc-rvc*rvbc)*xab*rrab**2))
     x        * gamc*rrac/wwc
     x        +(rud*udn*rrab*(xad-rud*udx-(rdb-rud*rubd)*xab*rrab**2)
     x        + rvd*vdn*rrab*(xad-rvd*vdx-(rdb-rvd*rvbd)*xab*rrab**2))
     x        * gamd*rrad/wwd
            
            fby=gamb*(-cosb*yab*rrab**2+rrab*(rub*uby+rvb*vby)/wwb)
     x        +(ruc*ucn*rrab*(yac-ruc*ucy-(rbc-ruc*rubc)*yab*rrab**2)
     x        - rvc*vcn*rrab*(yac-rvc*vcy-(rbc-rvc*rvbc)*yab*rrab**2))
     x        * gamc*rrac/wwc
     x        +(rud*udn*rrab*(yad-rud*udy-(rdb-rud*rubd)*yab*rrab**2)
     x        + rvd*vdn*rrab*(yad-rvd*vdy-(rdb-rvd*rvbd)*yab*rrab**2))
     x        * gamd*rrad/wwd
            
            fbz=gamb*(-cosb*zab*rrab**2+rrab*(rub*ubz+rvb*vbz)/wwb)
     x        +(ruc*ucn*rrab*(zac-ruc*ucz-(rbc-ruc*rubc)*zab*rrab**2)
     x        - rvc*vcn*rrab*(zac-rvc*vcz-(rbc-rvc*rvbc)*zab*rrab**2))
     x        * gamc*rrac/wwc
     x        +(rud*udn*rrab*(zad-rud*udz-(rdb-rud*rubd)*zab*rrab**2)
     x        + rvd*vdn*rrab*(zad-rvd*vdz-(rdb-rvd*rvbd)*zab*rrab**2))
     x        * gamd*rrad/wwd
            
            fcx=gamc*(-cosc*xac*rrac**2+rrac*(ruc*ucx+rvc*vcx)/wwc)
     x        +(rud*udn*rrac*(xad-rud*udx-(rcd-rud*rucd)*xac*rrac**2)
     x        - rvd*vdn*rrac*(xad-rvd*vdx-(rcd-rvd*rvcd)*xac*rrac**2))
     x        * gamd*rrad/wwd
     x        +(rub*ubn*rrac*(xab-rub*ubx-(rbc-rub*rucb)*xac*rrac**2)
     x        + rvb*vbn*rrac*(xab-rvb*vbx-(rbc-rvb*rvcb)*xac*rrac**2))
     x        * gamb*rrab/wwb
            
            fcy=gamc*(-cosc*yac*rrac**2+rrac*(ruc*ucy+rvc*vcy)/wwc)
     x        +(rud*udn*rrac*(yad-rud*udy-(rcd-rud*rucd)*yac*rrac**2)
     x        - rvd*vdn*rrac*(yad-rvd*vdy-(rcd-rvd*rvcd)*yac*rrac**2))
     x        * gamd*rrad/wwd
     x        +(rub*ubn*rrac*(yab-rub*uby-(rbc-rub*rucb)*yac*rrac**2)
     x        + rvb*vbn*rrac*(yab-rvb*vby-(rbc-rvb*rvcb)*yac*rrac**2))
     x        * gamb*rrab/wwb
            
            fcz=gamc*(-cosc*zac*rrac**2+rrac*(ruc*ucz+rvc*vcz)/wwc)
     x        +(rud*udn*rrac*(zad-rud*udz-(rcd-rud*rucd)*zac*rrac**2)
     x        - rvd*vdn*rrac*(zad-rvd*vdz-(rcd-rvd*rvcd)*zac*rrac**2))
     x        * gamd*rrad/wwd
     x        +(rub*ubn*rrac*(zab-rub*ubz-(rbc-rub*rucb)*zac*rrac**2)
     x        + rvb*vbn*rrac*(zab-rvb*vbz-(rbc-rvb*rvcb)*zac*rrac**2))
     x        * gamb*rrab/wwb
            
            fdx=gamd*(-cosd*xad*rrad**2+rrad*(rud*udx+rvd*vdx)/wwd)
     x        +(rub*ubn*rrad*(xab-rub*ubx-(rdb-rub*rudb)*xad*rrad**2)
     x        - rvb*vbn*rrad*(xab-rvb*vbx-(rdb-rvb*rvdb)*xad*rrad**2))
     x        * gamb*rrab/wwb
     x        +(ruc*ucn*rrad*(xac-ruc*ucx-(rcd-ruc*rudc)*xad*rrad**2)
     x        + rvc*vcn*rrad*(xac-rvc*vcx-(rcd-rvc*rvdc)*xad*rrad**2))
     x        * gamc*rrac/wwc
            
            fdy=gamd*(-cosd*yad*rrad**2+rrad*(rud*udy+rvd*vdy)/wwd)
     x        +(rub*ubn*rrad*(yab-rub*uby-(rdb-rub*rudb)*yad*rrad**2)
     x        - rvb*vbn*rrad*(yab-rvb*vby-(rdb-rvb*rvdb)*yad*rrad**2))
     x        * gamb*rrab/wwb
     x        +(ruc*ucn*rrad*(yac-ruc*ucy-(rcd-ruc*rudc)*yad*rrad**2)
     x        + rvc*vcn*rrad*(yac-rvc*vcy-(rcd-rvc*rvdc)*yad*rrad**2))
     x        * gamc*rrac/wwc
            
            fdz=gamd*(-cosd*zad*rrad**2+rrad*(rud*udz+rvd*vdz)/wwd)
     x        +(rub*ubn*rrad*(zab-rub*ubz-(rdb-rub*rudb)*zad*rrad**2)
     x        - rvb*vbn*rrad*(zab-rvb*vbz-(rdb-rvb*rvdb)*zad*rrad**2))
     x        * gamb*rrab/wwb
     x        +(ruc*ucn*rrad*(zac-ruc*ucz-(rcd-ruc*rudc)*zad*rrad**2)
     x        + rvc*vcn*rrad*(zac-rvc*vcz-(rcd-rvc*rvdc)*zad*rrad**2))
     x        * gamc*rrac/wwc
            
            fax=-(fbx+fcx+fdx)
            fay=-(fby+fcy+fdy)
            faz=-(fbz+fcz+fdz)
            
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
            
c     stress tensor calculation for inversion terms
            
            strs(1)=strs(1)+xab*fbx+xac*fcx+xad*fdx 
            strs(2)=strs(2)+yab*fbx+yac*fcx+yad*fdx 
            strs(3)=strs(3)+zab*fbx+zac*fcx+zad*fdx 
            strs(4)=strs(4)+yab*fby+yac*fcy+yad*fdy 
            strs(5)=strs(5)+yab*fbz+yac*fcz+yad*fdz 
            strs(6)=strs(6)+zab*fbz+zac*fcz+zad*fdz 
          
          endif
          
c     metadynamics local parameters
          
          if(lmetadyn.and.(idrive.or.jdrive.or.kdrive.or.ldrive))then
            
c     local energy (no virial)

            eng_loc=eng_loc+pterm

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
          
          endif
          endif
          endif
          endif
          endif
          endif
          endif
          endif
          endif
          endif
          endif

          enddo
          enddo
          enddo
          
          endif

        enddo

      enddo

c     check for undefined potentials

      if(mxnode.gt.1) call gstate(safe)
      if(.not.safe) call error(idnode,443)

c     global sum of four body potential and virial

      buffer(1)=engfbp
      buffer(2)=virfbp
      buffer(3)=fbp_fre
      buffer(4)=fbp_vir
      call gdsum(buffer(1),4,buffer(5))
      engfbp=buffer(1)
      virfbp=buffer(2)
      fbp_fre=buffer(3)
      fbp_vir=buffer(4)

c     sum up solvation energies
        
        if(lsolva)then

          call gdsum(en4_sol,mxtmls_sol4,buffer(1))
          if(lexcite)call gdsum(en4_exc,mxtmls_exc4,buffer(1))
          
        endif
        
c     restore coordinate array to original representation

      do i=1,natms

        if(filter(ltype(i)))then
          
          sxx=xxx(i)
          syy=yyy(i)
          szz=zzz(i)
          
          xxx(i)=cell(1)*sxx+cell(4)*syy+cell(7)*szz
          yyy(i)=cell(2)*sxx+cell(5)*syy+cell(8)*szz
          zzz(i)=cell(3)*sxx+cell(6)*syy+cell(9)*szz
          
        endif
        
      enddo

c     restore cell vector

      if(imcon.eq.0) then
        cell(1)=0.d0
        cell(5)=0.d0
        cell(9)=0.d0
      endif

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
      end subroutine fbpfrc
      
      end module four_body_module
