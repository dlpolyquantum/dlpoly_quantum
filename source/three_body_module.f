      module three_body_module

c***********************************************************************
c     
c     dl_poly module for defining three-body potential arrays
c     copyright - daresbury laboratory
c     author  - w. smith  sep 2003
c     adapted - w. smith  aug 2008 : solvation, free energy excitation
c     adapted - w. smith  jan 2011 : metadynamics
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

      logical, allocatable :: filter(:)
      real(8), allocatable :: prmtbp(:,:),rcut3b(:)
      integer, allocatable :: lsttbp(:),ltptbp(:),lattbp(:)

      save prmtbp,rcut3b,lsttbp,ltptbp,lattbp,filter

      contains
      
      subroutine alloc_tbp_arrays(idnode,mxnode)
      
      implicit none

      logical safe
      integer i,fail,idnode,mxnode
      dimension fail(6)
      
      safe=.true.

c     allocate arrays

      fail(:)=0
      
      allocate (prmtbp(mxtbp,mxptbp),stat=fail(1))
      allocate (rcut3b(mxtbp),stat=fail(2))
      allocate (lsttbp(mxtbp),stat=fail(3))
      allocate (ltptbp(mxtbp),stat=fail(4))
      allocate (lattbp(mxatms),stat=fail(5))
      allocate (filter(mxsite),stat=fail(6))
      
      if(any(fail.gt.0))safe=.false.      
      if(mxnode.gt.1)call gstate(safe)    
      if(.not.safe)call error(idnode,1170)
      
      end subroutine alloc_tbp_arrays

      subroutine define_three_body
     x  (safe,lunits,lmols,idnode,ntptbp,ntpatm,rcuttb,engunit)

c***********************************************************************
c     
c     dl_poly subroutine for defining three body potentials
c     
c     copyright - daresbury laboratory
c     author    - w. smith    aug 2003
c     
c***********************************************************************

      implicit none

      logical safe,lunits,lmols
      character*8 keyword
      character*8 atom0,atom1,atom2
      character*1 message(80)
      integer idnode,ntptbp,ntpatm,fail,i,itbp,itptbp,keypot
      integer idum,katm1,katm2,katm0,j,keytbp,ktbp,jtpatm
      real(8) rcuttb,engunit
      real(8), allocatable :: parpot(:)

      data fail/0/

      allocate (parpot(mxptbp),stat=fail)
      if(fail.ne.0)call error(idnode,1180)

      ntptbp=intstr(record,lenrec,idum)
      
      if(idnode.eq.0) then
        
        write(nrite,"(/,/,1x,'number of specified three ',
     x    'body potentials',i10)") ntptbp
        write(nrite,"(/,/,16x,'atom 1  ','atom 2  ','atom 3  ',
     x    3x,' key',30x,'parameters'/,/)")
        
      endif      
      if(ntptbp.gt.mxtbp) call error(idnode,83)
      if(.not.lunits) call error(idnode,6)
      if(.not.lmols) call error(idnode,13)
      
      do i=1,mxsite
        filter(i)=.false.
      enddo
      
      do itbp=1,mxtbp
        lsttbp(itbp)=0
      enddo
      
      do itbp=1,mxtbp,mx2tbp
        lsttbp(itbp)=-1
      enddo
      
      rcuttb=0.d0
      
      do itptbp=1,ntptbp
        
        do i=1,mxptbp
          parpot(i)=0.d0
        enddo
        
        call getrec(safe,idnode,nfield)
        if(.not.safe)return
        
c     Note the order!! atom0 is central 
        
        call copystring(record,message,80)
        call getword(atom1,record,8,lenrec)
        call getword(atom0,record,8,lenrec)
        call getword(atom2,record,8,lenrec)
        call lowcase(record,lenrec)
        call getword(keyword,record,4,lenrec)
        
        if(keyword(1:4).eq.'harm') then
          keypot=0
        elseif(keyword(1:4).eq.'thrm') then
          keypot=1
        elseif(keyword(1:4).eq.'shrm') then
          keypot=2
        elseif(keyword(1:4).eq.'bvs1') then
          keypot=3
        elseif(keyword(1:4).eq.'bvs2') then
          keypot=4
        elseif(keyword(1:4).eq.'hbnd') then
          keypot=5
        else
          if(idnode.eq.0) write(nrite,*) message
          call error(idnode,442)
        endif
        
        parpot(1)=dblstr(record,lenrec,idum)
        parpot(2)=dblstr(record,lenrec,idum)
        parpot(3)=dblstr(record,lenrec,idum)
        parpot(4)=dblstr(record,lenrec,idum)
        parpot(5)=dblstr(record,lenrec,idum)
        
        if(idnode.eq.0) 
     x    write(nrite,"(16x,3a8,4x,a4,1x,1p,9e13.5)") 
     x    atom1,atom0,atom2,keyword(1:4),(parpot(j),j=1,mxptbp)
        
        katm0=0
        katm1=0
        katm2=0
        
        do jtpatm=1,ntpatm
          
          if(atom0.eq.unqatm(jtpatm))katm0=jtpatm
          if(atom1.eq.unqatm(jtpatm))katm1=jtpatm
          if(atom2.eq.unqatm(jtpatm))katm2=jtpatm
          
        enddo
        
        if(katm0.eq.0.or.katm1.eq.0.or.katm2.eq.0) 
     x    call error(idnode,84)
        
        filter(katm0)=.true.
        filter(katm1)=.true.
        filter(katm2)=.true.
        
        keytbp=(max(katm1,katm2)*(max(katm1,katm2)-1))/2+
     x    min(katm1,katm2)+(katm0-1)*mx2tbp
        
        if(keytbp.gt.mxtbp) call error(idnode,86)

c     convert parameters to internal units
        
        parpot(1)=parpot(1)*engunit
        if(keypot.ne.5)parpot(2)=parpot(2)*(pi/180.d0)
        
        if(lsttbp(keytbp).gt.0) call error(idnode,18)
        lsttbp(keytbp)=itptbp
        ltptbp(itptbp)=keypot
        ktbp=mx2tbp*((keytbp-1)/mx2tbp)+1
        if(lsttbp(ktbp).lt.0)lsttbp(ktbp)=0

c     calculate max three body cutoff
        
        rcuttb=max(rcuttb,parpot(5))
        rcut3b(itptbp)=parpot(5)

c     store three body potential parameters
        
        do i=1,4
          prmtbp(itptbp,i)=parpot(i)
        enddo
        if(mxptbp.ge.6) then
          do i=6,mxptbp
            prmtbp(itptbp,i-1)=parpot(i-1)
          enddo
        endif
      enddo

      if(rcuttb.lt.1.d-6)call error(idnode,451)
      
      deallocate (parpot,stat=fail)

      return
      end subroutine define_three_body
      
      subroutine thbfrc
     x  (lsolva,lfree,lexcite,idnode,mxnode,natms,imcon,rcuttb,
     x  engtbp,virtbp)

c***********************************************************************
c     
c     dl_poly subroutine for calculating three body forces arising 
c     from the included angle between three atoms
c
c     Note: the subroutine converts coordinates to reduced units
c     to avoid a call to images.f. The link cell algorithm used
c     here necessitates a parallelepiped cell geometry
c     
c     copyright - daresbury laboratory 1994
c     author    - w.smith mar 1994 
c     adapted   - w.smith aug 2008 solvation, free energy etc
c     
c***********************************************************************

      implicit none

      logical safe,lsolva,lfree,lexcite,lselect,lskip
      logical idrive,jdrive,kdrive
      integer idnode,mxnode,natms,imcon,nix,niy,niz
      integer i,nbx,nby,nbz,ncells,l,ix,iy,iz,k,icell,kk,jx,jy,jz
      integer j,jcell,ii,itbp,limit,last,ktbp,jtbp,jktbp,kktbp
      integer ia,ib,ic,ktyp,jj,jk,kkk
      real(8) rcuttb,engtbp,virtbp,tterm,uterm,xm,ym,zm,cprp,det
      real(8) xdc,ydc,zdc,sxx,syy,szz,sxab,syab,szab,xab,yab,zab
      real(8) rab,sxbc,sybc,szbc,xbc,ybc,zbc,rbc,xac,yac,zac,rac
      real(8) rrab,rrbc,rrac,cost,sint,theta,pterm,gamma,vterm
      real(8) gamsa,gamsb,gamsc,scrn,fxa,fya,fza,fxc,fyc,fzc,strs
      real(8) strs_loc
      
      dimension nix(27),niy(27),niz(27),cprp(10),strs(6),strs_loc(9)
      
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
      
      engtbp=0.d0
      virtbp=0.d0
      tbp_fre=0.d0
      tbp_vir=0.d0
      strs(:)=0.d0
      strs_loc(:)=0.d0
      
      if(lsolva)then
        
        lcomp(8)=.true.
        en3_sol(:)=0.d0
        if(lexcite)en3_exc(:)=0.d0
        
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
        
        cell(1)=2.d0*xm+rcuttb
        cell(2)=0.d0
        cell(3)=0.d0
        cell(4)=0.d0
        cell(5)=2.d0*ym+rcuttb
        cell(6)=0.d0
        cell(7)=0.d0
        cell(8)=0.d0
        cell(9)=2.d0*zm+rcuttb
        
      endif
      
c     check for appropriate boundary conditions
      
      if(imcon.gt.3)call error(idnode,67)
      call invert(cell,rcell,det)
      call dcell(cell,cprp)
      
c     calculate link cell numbers
      
      nbx=int(cprp(7)/(rcuttb+1.d-6))
      nby=int(cprp(8)/(rcuttb+1.d-6))
      nbz=int(cprp(9)/(rcuttb+1.d-6))
      ncells=nbx*nby*nbz
      if(ncells.gt.mxcell) then
        
        if(idnode.eq.0) then
          
          write(nrite,'(a,i6)')
     x      'number of required link cells in routine thbfrc is ',ncells
          write(nrite,'(a,i6)')
     x      'number of default link cells in routine thbfrc is ',mxcell
          call error(idnode,69)
          
        endif
        
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
      
c     loop over central atoms of angles
      
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
        
c     construct mini-list of neighbour cell contents
        
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
            lattbp(k)=j
            j=link(j)
            
          enddo
          
        enddo
        
        limit=k
        
        do ii=1,lst(icell)
          
          i=lattbp(ii)
          itbp=mx2tbp*(ltype(i)-1)
          if(mod(i,mxnode).eq.idnode.and.lsttbp(itbp+1).ge.0)then
          
          last=limit
          
          do kk=1,limit/2
              
          if(kk.gt.(limit-1)/2)last=limit/2
          
          do jj=1,last
          
          j=lattbp(jj)
          jk=jj+kk
          if(jk.gt.limit)jk=jk-limit
          k=lattbp(jk)
          if(i.ne.j.and.i.ne.k)then
                  
          jtbp=max(ltype(j),ltype(k))
          ktbp=min(ltype(j),ltype(k))
          jktbp=itbp+(jtbp*(jtbp-1))/2+ktbp
          kktbp=lsttbp(jktbp)
          if(kktbp.gt.0)then
          
c     make labels etc consistent with angfrc.f
          
          ia=j
          ib=i
          ic=k
          
          if(lmetadyn)then
            
            idrive=driven(ltype(ia))
            jdrive=driven(ltype(ib))
            kdrive=driven(ltype(ic))
            
          endif
          
          if(lskip)then
            
            if((atm_fre(ia).eq.1.or.atm_fre(ib).eq.1.or.
     x        atm_fre(ic).eq.1).and.(atm_fre(ia).eq.2.or.
     x        atm_fre(ib).eq.2.or.atm_fre(ic).eq.2))cycle
            
          endif
          
          sxab=xxx(ia)-xxx(ib)
          sxab=sxab-nint(sxab)
          syab=yyy(ia)-yyy(ib)
          syab=syab-nint(syab)
          szab=zzz(ia)-zzz(ib)
          szab=szab-nint(szab)
          
          xab=cell(1)*sxab+cell(4)*syab+cell(7)*szab
          if(abs(xab).lt.rcuttb)then
          
          yab=cell(2)*sxab+cell(5)*syab+cell(8)*szab
          if(abs(yab).lt.rcuttb)then
          
          zab=cell(3)*sxab+cell(6)*syab+cell(9)*szab
          if(abs(zab).lt.rcuttb)then
                          
          sxbc=xxx(ic)-xxx(ib)
          sxbc=sxbc-nint(sxbc)
          sybc=yyy(ic)-yyy(ib)
          sybc=sybc-nint(sybc)
          szbc=zzz(ic)-zzz(ib)
          szbc=szbc-nint(szbc)
          
          xbc=cell(1)*sxbc+cell(4)*sybc+cell(7)*szbc
          if(abs(xbc).lt.rcuttb)then
                            
          ybc=cell(2)*sxbc+cell(5)*sybc+cell(8)*szbc
          if(abs(ybc).lt.rcuttb)then
          
          zbc=cell(3)*sxbc+cell(6)*sybc+cell(9)*szbc
          if(abs(zbc).lt.rcuttb)then
          
          ktyp=ltptbp(kktbp)
          rab=sqrt(xab*xab+yab*yab+zab*zab)
          rbc=sqrt(xbc*xbc+ybc*ybc+zbc*zbc)
          
          if(rcut3b(kktbp).ge.max(rab,rbc))then
          
          xac=xab-xbc
          yac=yab-ybc
          zac=zab-zbc
          rac=sqrt(xac*xac+yac*yac+zac*zac)
          
          rrab=1.d0/rab
          rrbc=1.d0/rbc
          rrac=1.d0/rac
          
c     normalise direction vectors
          
          xab=xab*rrab
          yab=yab*rrab
          zab=zab*rrab
          
          xbc=xbc*rrbc
          ybc=ybc*rrbc
          zbc=zbc*rrbc
          
          xac=xac*rrac
          yac=yac*rrac
          zac=zac*rrac
          
          cost=(xab*xbc+yab*ybc+zab*zbc)
          if(abs(cost).gt.1.d0)cost=sign(1.d0,cost)
          if(ktyp.ne.5)then
            
            sint=max(1.d-8,sqrt(1.d0-cost*cost))
            theta=acos(cost)
            
          endif
          
          if(ktyp.eq.0)then
            
c     harmonic angle potential
            
            pterm=0.5d0*prmtbp(kktbp,1)*(theta-prmtbp(kktbp,2))**2
            gamma=prmtbp(kktbp,1)*(theta-prmtbp(kktbp,2))/sint
            
            vterm=0.d0
            gamsa=0.d0
            gamsc=0.d0
            gamsb=0.d0
            
          elseif(ktyp.eq.1)then
            
c     truncated harmonic valence angle potential
            
            scrn=exp(-(rab**8+rbc**8)/prmtbp(kktbp,3)**8)
            pterm=scrn*0.5d0*prmtbp(kktbp,1)*(theta-prmtbp(kktbp,2))**2
            vterm=-8.d0*pterm*(rab**8+rbc**8)/prmtbp(kktbp,3)**8
            gamma=scrn*prmtbp(kktbp,1)*(theta-prmtbp(kktbp,2))/sint
            gamsa=(8.d0*pterm/prmtbp(kktbp,3)**8)*rab**7
            gamsc=(8.d0*pterm/prmtbp(kktbp,3)**8)*rbc**7
            gamsb=0.d0
            
          elseif(ktyp.eq.2)then
            
c     screened harmonic valence angle potential
            
            scrn=exp(-(rab/prmtbp(kktbp,3)+rbc/prmtbp(kktbp,4)))
            pterm=scrn*0.5d0*prmtbp(kktbp,1)*(theta-prmtbp(kktbp,2))**2
            vterm=-pterm*(rab/prmtbp(kktbp,3)+rbc/prmtbp(kktbp,4))
            gamma=scrn*prmtbp(kktbp,1)*(theta-prmtbp(kktbp,2))/sint
            gamsa=(pterm/prmtbp(kktbp,3))
            gamsc=(pterm/prmtbp(kktbp,4))
            gamsb=0.d0
            
          elseif(ktyp.eq.3)then
            
c     screened vessal potential type 1
            
            scrn=exp(-(rab/prmtbp(kktbp,3)+rbc/prmtbp(kktbp,4)))
            pterm=scrn*prmtbp(kktbp,1)/
     x        (8.d0*(prmtbp(kktbp,2)-pi)**2)*
     x        ((prmtbp(kktbp,2)-pi)**2-(theta-pi)**2)**2
            vterm=-pterm*(rab/prmtbp(kktbp,3)+rbc/prmtbp(kktbp,4))
            gamma=scrn*prmtbp(kktbp,1)/
     x        (2.d0*(prmtbp(kktbp,2)-pi)**2)*
     x        ((prmtbp(kktbp,2)-pi)**2-(theta-pi)**2)*(theta-pi)/sint
            gamsa=(pterm/prmtbp(kktbp,3))
            gamsc=(pterm/prmtbp(kktbp,4))
            gamsb=0.d0
            
          elseif(ktyp.eq.4)then
            
c     truncated vessal potential type 2 - use with sw1
            
            scrn=exp(-(rab**8+rbc**8)/prmtbp(kktbp,4)**8)
            pterm=scrn*prmtbp(kktbp,1)*(theta**prmtbp(kktbp,3)*
     x        (theta-prmtbp(kktbp,2))**2*(theta+prmtbp(kktbp,2)-
     x        2.d0*pi)**2-0.5d0*prmtbp(kktbp,3)*pi**(prmtbp(kktbp,3)
     x        -1.d0)*(theta-prmtbp(kktbp,2))**2*(pi-
     x        prmtbp(kktbp,2))**3)
            vterm=-8.d0*pterm*(rab**8+rbc**8)/prmtbp(kktbp,4)**8
            gamma=scrn*prmtbp(kktbp,1)*(theta**(prmtbp(kktbp,3)-1.d0)*
     x        (theta-prmtbp(kktbp,2))*(theta+prmtbp(kktbp,2)-
     x        2.d0*pi)*((prmtbp(kktbp,3)+4.d0)*theta**2-2.d0*pi*
     x        (prmtbp(kktbp,3)+2.d0)*theta+prmtbp(kktbp,3)*
     x        prmtbp(kktbp,2)*(2.d0*pi-prmtbp(kktbp,2)))-
     x        prmtbp(kktbp,3)*pi**(prmtbp(kktbp,3)-1.d0)*
     x        (theta-prmtbp(kktbp,2))*(pi-prmtbp(kktbp,2))**3)/sint
            gamsa=(8.d0*pterm/prmtbp(kktbp,4)**8)*rab**7
            gamsc=(8.d0*pterm/prmtbp(kktbp,4)**8)*rbc**7
            gamsb=0.d0
            
          elseif(ktyp.eq.5)then
            
            if(min(rab,rbc).lt.1.5d0)then
              
              scrn=(5.d0*(prmtbp(kktbp,2)/rac)**2-6.d0)*
     x          (prmtbp(kktbp,2)/rac)**10
              tterm=prmtbp(kktbp,1)*cost**4
              pterm=scrn*tterm
              uterm=60.d0*((prmtbp(kktbp,2)/rac)**2-1.d0)*
     x          (prmtbp(kktbp,2)/rac)**10
              vterm=tterm*uterm
              gamma=scrn*4.d0*prmtbp(kktbp,1)*cost**3
              gamsb=tterm*uterm/rac
              gamsa=0.d0
              gamsc=0.d0
              
            endif
            
          else
            
            safe=.false.
            pterm=0.d0
            vterm=0.d0
            gamma=0.d0
            gamsa=0.d0
            gamsb=0.d0
            gamsc=0.d0
            
          endif
          
c     set selection control
          
          lselect=.true.
          
c     set triple index
          
          if(lsolva)kkk=loc3(atmolt(ia),atmolt(ib),atmolt(ic))
          
          if(lexcite)then
            
c     selected excitation option
            
            if((atm_fre(ia).ne.1).and.(atm_fre(ib).ne.1)
     x        .and.(atm_fre(ic).ne.1))then
              
c     reset selection control
              
              lselect=(atm_fre(ia)+atm_fre(ib)+atm_fre(ic).eq.0)
              
              if(lsolva)en3_exc(kkk)=en3_exc(kkk)+pterm
              
            endif
            
          elseif(lfree)then
            
c     selected free energy option
            
            if((atm_fre(ia).eq.1).or.(atm_fre(ib).eq.1)
     x        .or.(atm_fre(ic).eq.1))then
              
c     set hamiltonian mixing parameter
              
              tbp_fre=tbp_fre-pterm
              tbp_vir=tbp_vir-vterm
              pterm=lambda1*pterm
              vterm=lambda1*vterm
              gamma=lambda1*gamma
              gamsa=lambda1*gamsa
              gamsb=lambda1*gamsb
              gamsc=lambda1*gamsc
              
            elseif((atm_fre(ia).eq.2).or.(atm_fre(ib).eq.2)
     x          .or.(atm_fre(ic).eq.2))then
              
c     set hamiltonian mixing parameter
              
              tbp_fre=tbp_fre+pterm
              tbp_vir=tbp_vir+vterm
              pterm=lambda2*pterm
              vterm=lambda2*vterm
              gamma=lambda2*gamma
              gamsa=lambda2*gamsa
              gamsb=lambda2*gamsb
              gamsc=lambda2*gamsc
              
            endif
            
          endif
          
          if(lselect)then
            
c     calculate potential and virial
            
            engtbp=engtbp+pterm
            virtbp=virtbp+vterm
            
            if(lsolva)en3_sol(kkk)=en3_sol(kkk)+pterm
            
c     calculate atomic forces
            
            fxa=gamma*(xbc-xab*cost)*rrab+gamsa*xab+gamsb*xac
            fya=gamma*(ybc-yab*cost)*rrab+gamsa*yab+gamsb*yac
            fza=gamma*(zbc-zab*cost)*rrab+gamsa*zab+gamsb*zac
            
            fxc=gamma*(xab-xbc*cost)*rrbc+gamsc*xbc-gamsb*xac
            fyc=gamma*(yab-ybc*cost)*rrbc+gamsc*ybc-gamsb*yac
            fzc=gamma*(zab-zbc*cost)*rrbc+gamsc*zbc-gamsb*zac
            
            fxx(ia)=fxx(ia)+fxa
            fyy(ia)=fyy(ia)+fya
            fzz(ia)=fzz(ia)+fza
            
            fxx(ib)=fxx(ib)-fxa-fxc
            fyy(ib)=fyy(ib)-fya-fyc
            fzz(ib)=fzz(ib)-fza-fzc
            
            fxx(ic)=fxx(ic)+fxc
            fyy(ic)=fyy(ic)+fyc
            fzz(ic)=fzz(ic)+fzc
            
c     calculate stress tensor
            
            strs(1)=strs(1)+rab*xab*fxa+rbc*xbc*fxc
            strs(2)=strs(2)+rab*xab*fya+rbc*xbc*fyc
            strs(3)=strs(3)+rab*xab*fza+rbc*xbc*fzc
            strs(4)=strs(4)+rab*yab*fya+rbc*ybc*fyc
            strs(5)=strs(5)+rab*yab*fza+rbc*ybc*fzc
            strs(6)=strs(6)+rab*zab*fza+rbc*zbc*fzc
            
          endif

c     metadynamics local parameters
        
          if(lmetadyn)then
            
c     local energy and virial
            
            eng_loc=eng_loc+pterm
            vir_loc=vir_loc+vterm
            
c     local forces
            
            fxx_loc(ia)=fxx_loc(ia)+fxa
            fyy_loc(ia)=fyy_loc(ia)+fya
            fzz_loc(ia)=fzz_loc(ia)+fza
            
            fxx_loc(ib)=fxx_loc(ib)-fxa-fxc
            fyy_loc(ib)=fyy_loc(ib)-fya-fyc
            fzz_loc(ib)=fzz_loc(ib)-fza-fzc
            
            fxx_loc(ic)=fxx_loc(ic)+fxc
            fyy_loc(ic)=fyy_loc(ic)+fyc
            fzz_loc(ic)=fzz_loc(ic)+fzc
            
c     local stress tensor
            
            strs_loc(1)=strs_loc(1)+rab*xab*fxa+rbc*xbc*fxc
            strs_loc(2)=strs_loc(2)+rab*xab*fya+rbc*xbc*fyc
            strs_loc(3)=strs_loc(3)+rab*xab*fza+rbc*xbc*fzc
            strs_loc(4)=strs_loc(4)+rab*yab*fya+rbc*ybc*fyc
            strs_loc(5)=strs_loc(5)+rab*yab*fza+rbc*ybc*fzc
            strs_loc(6)=strs_loc(6)+rab*zab*fza+rbc*zbc*fzc
            
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
          
          endif

        enddo

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

      if(mxnode.gt.1) call gstate(safe)
      if(.not.safe) call error(idnode,442)

c     global sum of three body potential and virial

      if(mxnode.gt.1)then

        buffer(1)=engtbp
        buffer(2)=virtbp
        buffer(3)=tbp_fre
        buffer(4)=tbp_vir
        call gdsum(buffer(1),4,buffer(5))
        engtbp=buffer(1)
        virtbp=buffer(2)
        tbp_fre=buffer(3)
        tbp_vir=buffer(4)

c     sum up solvation energies
        
        if(lsolva)then

          call gdsum(en3_sol,mxtmls_sol3,buffer(1))
          if(lexcite)call gdsum(en3_exc,mxtmls_exc3,buffer(1))
          
        endif
        
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
      
      return

      end subroutine thbfrc
      
      end module three_body_module
