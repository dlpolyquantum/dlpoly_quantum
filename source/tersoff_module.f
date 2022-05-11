      module tersoff_module

c***********************************************************************
c     
c     dl_poly module for defining tersoff potential arrays
c     
c     copyright - daresbury laboratory
c     author    - w. smith    dec 2003
c     
c***********************************************************************

      use config_module
      use error_module
      use parse_module
      use setup_module
      use site_module
      use utility_module

      implicit none

      logical, allocatable :: filter(:)
      integer, allocatable :: lstter(:),ltpter(:),lattsf(:)
      real(8), allocatable :: prmter(:,:),prmter2(:,:)
      real(8), allocatable :: vmbp(:,:,:),gmbp(:,:,:)
      real(8), allocatable :: xtf(:),ytf(:),ztf(:),rtf(:)
      real(8), allocatable :: ert(:),eat(:),grt(:),gat(:)
      real(8), allocatable :: scr(:),gcr(:),gam(:),gvr(:)

      save xtf,ytf,ztf,rtf,ert,eat,grt,gat,scr,gcr,gam,filter
      save prmter,prmter2,lstter,ltpter,lattsf,vmbp,gmbp

      contains
      
      subroutine alloc_ter_arrays(idnode,mxnode)
      
      implicit none
      
      integer, parameter :: nnn=20

      logical safe
      integer i,fail,idnode,mxnode,npairs
      dimension fail(nnn)

      safe=.true.

c     allocate arrays

      fail(:)=0
      
      npairs=(mxter*(mxter+1))/2
      allocate (prmter(mxter,mxpter),stat=fail(1))
      allocate (prmter2(2,npairs),stat=fail(2))
      allocate (lstter(mxter),stat=fail(3))
      allocate (ltpter(mxter),stat=fail(4))
      allocate (lattsf(mxatms),stat=fail(5))
      allocate (xtf(mxatms),stat=fail(6))
      allocate (ytf(mxatms),stat=fail(7))
      allocate (ztf(mxatms),stat=fail(8))
      allocate (rtf(mxatms),stat=fail(9))
      allocate (ert(mxatms),stat=fail(10))
      allocate (eat(mxatms),stat=fail(11))
      allocate (grt(mxatms),stat=fail(12))
      allocate (gat(mxatms),stat=fail(13))
      allocate (scr(mxatms),stat=fail(14))
      allocate (gcr(mxatms),stat=fail(15))
      allocate (gam(mxatms),stat=fail(16))
      allocate (gvr(mxatms),stat=fail(17))
      allocate (vmbp(mxgrid,npairs,3),stat=fail(18))
      allocate (gmbp(mxgrid,npairs,3),stat=fail(19))
      allocate (filter(mxsite),stat=fail(20))
      
      if(any(fail.gt.0))safe=.false.      
      if(mxnode.gt.1)call gstate(safe)    
      if(.not.safe)call error(idnode,1945)
      
      end subroutine alloc_ter_arrays

      subroutine define_tersoff
     x  (safe,lunits,lmols,idnode,ntpter,ntpatm,rctter,engunit)

c***********************************************************************
c     
c     dl_poly subroutine for defining tersoff potentials
c
c     based on potential form defined in:
c     J. Tersoff, Phys. Rev. B 39 (1989) 5566.
c     
c     copyright - daresbury laboratory
c     author    - w. smith    sep 2003
c     
c***********************************************************************

      implicit none

      logical safe,lunits,lmols
      character*8 keyword
      character*8 atom0,atom1
      character*1 message(80)
      integer fail,idnode,ntpter,ntpatm,i,idum,j,itpter
      integer keypot,jtpatm,k,katm0,katm1,icross,npairs,ktyp
      real(8) rctter,engunit
      real(8), allocatable :: parpot(:)
      data fail/0/

      allocate (parpot(mxpter),stat=fail)
      if(fail.ne.0)call error(idnode,1955)

      ntpter=intstr(record,lenrec,idum)
      
      if(idnode.eq.0) then
        
        write(nrite,"(/,/,1x,'number of specified tersoff ',
     x    'atom potentials',i10)") ntpter
        write(nrite,"(/,/,16x,'atom    ',3x,' key',30x,
     x    'parameters'/,/)")
        
      endif      
      if(ntpter.gt.mxter) call error(idnode,88)
      if(.not.lunits) call error(idnode,6)
      if(.not.lmols) call error(idnode,13)
      
      rctter=0.d0

      do i=1,mxter
        lstter(i)=-1
      enddo

      do i=1,mxsite
        filter(i)=.false.
      enddo

      k=0
      do i=1,mxter
        do j=1,i

          k=k+1
          prmter2(1,k)=0.d0
          prmter2(2,k)=0.d0

        enddo
      enddo
      
      do itpter=1,ntpter
        
        do i=1,mxpter
          parpot(i)=0.d0
        enddo
        
        call getrec(safe,idnode,nfield)
        if(.not.safe)return

        call copystring(record,message,80)
        call getword(atom0,record,8,lenrec)
        call lowcase(record,4)
        call getword(keyword,record,4,lenrec)

        if(keyword(1:4).eq.'ters') then

          keypot=1

        else

          if(idnode.eq.0) write(nrite,*) message
          call error(idnode,1972)

        endif

        parpot(1)=dblstr(record,lenrec,idum)  ! A_i
        parpot(2)=dblstr(record,lenrec,idum)  ! a_i
        parpot(3)=dblstr(record,lenrec,idum)  ! B_i
        parpot(4)=dblstr(record,lenrec,idum)  ! b_i
        parpot(5)=dblstr(record,lenrec,idum)  ! R_i

        call getrec(safe,idnode,nfield)
        if(.not.safe)return

        parpot(6)=dblstr(record,lenrec,idum)    ! S_i
        parpot(7)=dblstr(record,lenrec,idum)    ! beta_i
        parpot(8)=dblstr(record,lenrec,idum)    ! eta_i
        parpot(9)=dblstr(record,lenrec,idum)    ! c_i
        parpot(10)=dblstr(record,lenrec,idum)   ! d_i
        parpot(11)=dblstr(record,lenrec,idum)   ! h_i
        
        if(idnode.eq.0) then

          write(nrite,"(16x,a8,2x,a4,2x,1p,5e13.5)") 
     x      atom0,keyword(1:4),(parpot(j),j=1,5)
          write(nrite,"(32x,1p,5e13.5)")(parpot(j),j=6,mxpter)

        endif
        
        katm0=0
        
        do jtpatm=1,ntpatm
          if(atom0.eq.unqatm(jtpatm))katm0=jtpatm
        enddo
        
        if(katm0.eq.0)call error(idnode,92)

        filter(katm0)=.true.
        
c     convert parameters to internal units
        
        if(keypot.eq.1)then
          
          parpot(1)=parpot(1)*engunit
          parpot(3)=parpot(3)*engunit

        endif
        if(lstter(katm0).gt.-1) call error(idnode,21)
        lstter(katm0)=itpter
        ltpter(itpter)=keypot

c     calculate max tersoff cutoff
        
        rctter=max(rctter,parpot(6))
        
c     store tersoff single atom potential parameters
        
        do i=1,mxpter
          prmter(itpter,i)=parpot(i)
        enddo

      enddo

      if(rctter.lt.1.d-6)call error(idnode,1953)
      
c     start processing double atom potential parameters

      npairs=(ntpter*(ntpter+1))/2

      if(idnode.eq.0)then

        write(nrite,"(/,/,1x,'number of tersoff ',
     x    'cross terms',i10)") npairs
        write(nrite,"(/,/,16x,'atom    ','atom    ',10x,
     x    'parameters'/,/)")

      endif

      do icross=1,npairs

        call getrec(safe,idnode,nfield)
        if(.not.safe)return

        call getword(atom0,record,8,lenrec)
        call getword(atom1,record,8,lenrec)

        parpot(1)=dblstr(record,lenrec,idum)  ! chi_ij
        parpot(2)=dblstr(record,lenrec,idum)  ! omega_ij

        katm0=0
        katm1=0
        
        do jtpatm=1,ntpatm
          if(atom0.eq.unqatm(jtpatm))katm0=jtpatm
          if(atom1.eq.unqatm(jtpatm))katm1=jtpatm
        enddo
        
        if(katm0.eq.0.or.katm1.eq.0)call error(idnode,92)
        
        filter(katm0)=.true.
        filter(katm1)=.true.
        
        ktyp=loc2(lstter(katm0),lstter(katm1))
        prmter2(1,ktyp)=parpot(1)
        prmter2(2,ktyp)=parpot(2)

        if(idnode.eq.0)write(nrite,"(16x,a8,a8,1p,2e13.5)") 
     x     atom0,atom1,(parpot(j),j=1,2)

      enddo

c     generate tersoff interpolation arrays

      call tergen(ntpatm,rctter)

      return
      end subroutine define_tersoff

      subroutine tergen(ntpatm,rctter)

c***********************************************************************
c     
c     dl_poly subroutine for generating potential energy and 
c     force arrays for tersoff forces only
c
c     based on potential form defined in:
c     J. Tersoff, Phys. Rev. B 39 (1989) 5566.
c     
c     copyright - daresbury laboratory
c     author    - w. smith dec 2003
c     
c***********************************************************************
      
      implicit none

      integer ntpatm,katm0,katm1,ipt,jpt,kpt,i
      real(8) dlrpot,rctter,baij,saij,bbij,sbij,rij,sij,att,arg
      real(8) rrr,rep
      
c     define grid resolution for potential arrays
      
      dlrpot=rctter/dble(mxgrid-4)

c     construct arrays for all types of short ranged  potential
      
      do katm0=1,ntpatm

        if(filter(katm0))then
          
          ipt=lstter(katm0)
          
          do katm1=1,katm0
            
            if(filter(katm1))then
              
              jpt=lstter(katm1)
              
              if((ltpter(ipt).eq.1).and.(ltpter(jpt).eq.1))then
                
                kpt=loc2(ipt,jpt)
                
c     define tersoff parameters
                
                baij=sqrt(prmter(ipt,1)*prmter(jpt,1))
                saij=0.5d0*(prmter(ipt,2)+prmter(jpt,2))
                bbij=sqrt(prmter(ipt,3)*prmter(jpt,3))
                sbij=0.5d0*(prmter(ipt,4)+prmter(jpt,4))
                rij=sqrt(prmter(ipt,5)*prmter(jpt,5))
                sij=sqrt(prmter(ipt,6)*prmter(jpt,6))
                
c     store potential cutoff
                
                vmbp(1,kpt,1)=sij
                
c     calculate screening function
                
                do i=2,mxgrid
                  
                  rrr=dble(i)*dlrpot
                  
                  if(rrr.le.rij)then
                    
                    vmbp(i,kpt,1)=1.d0
                    gmbp(i,kpt,1)=0.d0
                    
                  else
                    
                    arg=pi*(rrr-rij)/(sij-rij)
                    vmbp(i,kpt,1)=0.5d0*(1.d0+cos(arg))
                    gmbp(i,kpt,1)=0.5d0*pi*rrr*sin(arg)/(sij-rij)
                    
                  endif
                  
                enddo
                
c     calculate screened repulsion function
                
                do i=2,mxgrid
                  
                  rrr=dble(i)*dlrpot
                  
                  rep=baij*exp(-saij*rrr)
                  vmbp(i,kpt,2)=rep*vmbp(i,kpt,1)
                  gmbp(i,kpt,2)=rep*(gmbp(i,kpt,1)+
     x              saij*rrr*vmbp(i,kpt,1))
                  
                enddo
                
c     calculate screened attraction function
                
                do i=2,mxgrid
                  
                  rrr=dble(i)*dlrpot
                  
                  att=bbij*exp(-sbij*rrr)
                  vmbp(i,kpt,3)=att*vmbp(i,kpt,1)
                  gmbp(i,kpt,3)=att*(gmbp(i,kpt,1)+
     x              sbij*rrr*vmbp(i,kpt,1))
                  
                enddo
                
              endif
              
            endif
            
          enddo
          
        endif
        
      enddo
          
      return
      end subroutine tergen

      subroutine tersoff
     x  (idnode,mxnode,natms,imcon,rctter,engter,virter)

c***********************************************************************
c     
c     dl_poly subroutine for calculating potential and forces 
c     due to a tersoff potential
c
c     Note: the subroutine converts coordinates to reduced units
c     to avoid a call to images.f. The link cell algorithm used
c     here necessitates a parallelepiped, cubic or orthorhombic
c     cell geometry
c
c     based on potential form defined in:
c     J. Tersoff, Phys. Rev. B 39 (1989) 5566.
c     
c     copyright - daresbury laboratory
c     author    - w.smith dec 2003
c     
c***********************************************************************

      implicit none

      logical safe
      integer idnode,mxnode,natms,imcon,nix,niy,niz,i,j,nbx
      integer nby,nbz,ncells,l,ix,iy,iz,k,icell,kk,jx,jy,jz,jcell,iter
      integer limit,ii,iatm
      real(8) rctter,engter,virter,xm,ym,zm,det,cprp
      real(8) xdc,ydc,zdc,sxx,syy,szz,strs

      dimension nix(27),niy(27),niz(27),cprp(10),strs(6)

      data nix/ 0,-1,-1,-1, 0, 0,-1, 1,-1, 0, 1,-1, 0, 1,
     x  1, 1, 1, 0, 0, 1,-1, 1, 0,-1, 1, 0,-1/
      data niy/ 0, 0,-1, 1, 1, 0, 0, 0,-1,-1,-1, 1, 1, 1,
     x  0, 1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1/
      data niz/ 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     x  0, 0, 0, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1/
      
c     flag for undefined potentials

      safe=.true.

c     initialise potential energy and virial

      engter=0.d0
      virter=0.d0

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

        cell(1)=2.d0*xm+rctter
        cell(2)=0.d0
        cell(3)=0.d0
        cell(4)=0.d0
        cell(5)=2.d0*ym+rctter
        cell(6)=0.d0
        cell(7)=0.d0
        cell(8)=0.d0
        cell(9)=2.d0*zm+rctter
        
      endif

c     initialise stress tensor accumulators

      strs(:)=0.d0

c     check for appropriate boundary conditions
      
      if(imcon.gt.3)call error(idnode,1974)
      call invert(cell,rcell,det)
      call dcell(cell,cprp)

c     calculate link cell numbers
      
      nbx=int(cprp(7)/(rctter+1.d-6))
      nby=int(cprp(8)/(rctter+1.d-6))
      nbz=int(cprp(9)/(rctter+1.d-6))
      if(nbx.lt.3.or.nby.lt.3.or.nbz.lt.3)then
        
        if(idnode.eq.0)then
          
          write(nrite,'(a,3i5)')
     x      'tersoff link cell decomposition is',nbx,nby,nbz
          
        endif
        
        call error(idnode,1977)
        
      endif
      ncells=nbx*nby*nbz

      if(ncells.gt.mxcell) then
        
        if(idnode.eq.0) then

          write(nrite,'(a,i6)')
     x      'number of required link cells in tersoff.f is ',ncells
          write(nrite,'(a,i6)')
     x      'number of  default link cells in tersoff.f is ',mxcell
          call error(idnode,1976)

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
            lattsf(k)=j
            j=link(j)
            
          enddo
          
        enddo

        limit=k
        
        do ii=1,lst(icell)
          
          iatm=lattsf(ii)
          iter=lstter(ltype(iatm))

          if(mod(iatm,mxnode).eq.idnode.and.iter.ge.0)then

c     construct working arrays by interpolation

            call terint(iatm,limit,rctter)

c     calculate three body (attractive) terms

            call tersoff3(ii,limit,engter,virter,strs)

          endif
          
        enddo

      enddo
      
c     calculate stress tensor
              
      stress(1)=stress(1)+strs(1)
      stress(2)=stress(2)+strs(2)
      stress(3)=stress(3)+strs(3)
      stress(4)=stress(4)+strs(2)
      stress(5)=stress(5)+strs(4)
      stress(6)=stress(6)+strs(5)
      stress(7)=stress(7)+strs(3)
      stress(8)=stress(8)+strs(5)
      stress(9)=stress(9)+strs(6)

c     check for undefined potentials

      if(mxnode.gt.1) call gstate(safe)
      if(.not.safe) call error(idnode,1978)

c     global sum of three body potential and virial

      if(mxnode.gt.1)then

        buffer(1)=engter
        buffer(2)=virter
        call gdsum(buffer(1),2,buffer(3))
        engter=buffer(1)
        virter=buffer(2)

      endif

c    remove effect of double counting

      engter=0.5d0*engter
      virter=0.5d0*virter

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
      end subroutine tersoff

      subroutine terint(iatm,limit,rctter)

c***********************************************************************
c     
c     dl_poly subroutine for constructing working arrays for
c     a tersoff potential using interpolation
c
c     based on potential form defined in:
c     J. Tersoff, Phys. Rev. B 39 (1989) 5566.
c     
c     copyright - daresbury laboratory
c     author    - w.smith dec 2003
c     
c***********************************************************************
      
      implicit none

      integer iatm,jatm,jj,limit,iter,jter,jjter,ll
      real(8) rctter,sxij,syij,szij,ppp,t1,t2,rdr
      real(8) vk0,vk1,vk2,gk0,gk1,gk2
     
      rdr=dble(mxgrid-4)/rctter
      iter=lstter(ltype(iatm))

c    initialise working arrays
            
      do jj=1,limit
        
        xtf(jj)=0.d0
        ytf(jj)=0.d0
        ztf(jj)=0.d0
        rtf(jj)=0.d0
        ert(jj)=0.d0
        eat(jj)=0.d0
        grt(jj)=0.d0
        gat(jj)=0.d0
        scr(jj)=0.d0
        gcr(jj)=0.d0
        
      enddo

c    construct working arrays
      
      do jj=1,limit
        
        jatm=lattsf(jj)
        jter=lstter(ltype(jatm))
        
        if(jatm.ne.iatm.and.jter.ge.0)then
          
          sxij=xxx(jatm)-xxx(iatm)
          sxij=sxij-nint(sxij)
          syij=yyy(jatm)-yyy(iatm)
          syij=syij-nint(syij)
          szij=zzz(jatm)-zzz(iatm)
          szij=szij-nint(szij)
          
          xtf(jj)=cell(1)*sxij+cell(4)*syij+cell(7)*szij
          ytf(jj)=cell(2)*sxij+cell(5)*syij+cell(8)*szij
          ztf(jj)=cell(3)*sxij+cell(6)*syij+cell(9)*szij
          rtf(jj)=sqrt(xtf(jj)**2+ytf(jj)**2+ztf(jj)**2)
          xtf(jj)=xtf(jj)/rtf(jj)
          ytf(jj)=ytf(jj)/rtf(jj)
          ztf(jj)=ztf(jj)/rtf(jj)
            
          jjter=loc2(iter,jter)
          if(rtf(jj).le.vmbp(1,jjter,1))then
            
            ll=int(rdr*rtf(jj))
            ppp=rtf(jj)*rdr-dble(ll)
            
c    interpolate screening function
            
            vk0=vmbp(ll,jjter,1)
            vk1=vmbp(ll+1,jjter,1)
            vk2=vmbp(ll+2,jjter,1)
            t1=vk0+(vk1-vk0)*ppp
            t2=vk1+(vk2-vk1)*(ppp-1.0d0)
            scr(jj)=t1+(t2-t1)*ppp*0.5d0
            
c    interpolate derivative of screening function
            
            gk0=gmbp(ll,jjter,1)
            gk1=gmbp(ll+1,jjter,1)
            gk2=gmbp(ll+2,jjter,1)
            t1=gk0+(gk1-gk0)*ppp
            t2=gk1+(gk2-gk1)*(ppp-1.0d0)
            gcr(jj)=-(t1+(t2-t1)*ppp*0.5d0)/rtf(jj)
            
c    interpolate repulsive component of energy
            
            vk0=vmbp(ll,jjter,2)
            vk1=vmbp(ll+1,jjter,2)
            vk2=vmbp(ll+2,jjter,2)
            t1=vk0+(vk1-vk0)*ppp
            t2=vk1+(vk2-vk1)*(ppp-1.0d0)
            ert(jj)=t1+(t2-t1)*ppp*0.5d0
            
c    interpolate derivative of repulsive function
            
            gk0=gmbp(ll,jjter,2)
            gk1=gmbp(ll+1,jjter,2)
            gk2=gmbp(ll+2,jjter,2)
            t1=gk0+(gk1-gk0)*ppp
            t2=gk1+(gk2-gk1)*(ppp-1.0d0)
            grt(jj)=-(t1+(t2-t1)*ppp*0.5d0)/rtf(jj)
            
c    interpolate attractive component of energy
            
            vk0=vmbp(ll,jjter,3)
            vk1=vmbp(ll+1,jjter,3)
            vk2=vmbp(ll+2,jjter,3)
            t1=vk0+(vk1-vk0)*ppp
            t2=vk1+(vk2-vk1)*(ppp-1.0d0)
            eat(jj)=t1+(t2-t1)*ppp*0.5d0
            
c    interpolate derivative of attractive function
            
            gk0=gmbp(ll,jjter,3)
            gk1=gmbp(ll+1,jjter,3)
            gk2=gmbp(ll+2,jjter,3)
            t1=gk0+(gk1-gk0)*ppp
            t2=gk1+(gk2-gk1)*(ppp-1.0d0)
            gat(jj)=-(t1+(t2-t1)*ppp*0.5d0)/rtf(jj)
            
          endif
          
        endif
        
      enddo
      
      return
      end subroutine terint

      subroutine tersoff3(ii,limit,engter,virter,strs)

c***********************************************************************
c     
c     dl_poly subroutine for calculating three body contributions
c     to a tersoff potential and tersoff potential, virial and
c     atomic forces
c
c     based on potential form defined in:
c     J. Tersoff, Phys. Rev. B 39 (1989) 5566.
c     
c     copyright - daresbury laboratory
c     author    - w.smith dec 2003
c     
c***********************************************************************

      implicit none

      logical flag
      integer iatm,jatm,katm,limit,ii,jj,kk
      integer iter,jter,kter,jjter,kkter
      real(8) cost,gtheta,ci,di,hi,gamma,gam_j,gam_k,eterm,gterm,strs
      real(8) fxa,fya,fza,fxc,fyc,fzc,engter,virter,vterm,gam_ij,bi,ei

      dimension strs(6)
      
      iatm=lattsf(ii)
      iter=lstter(ltype(iatm))

      bi=prmter(iter,7)
      ei=prmter(iter,8)
      ci=prmter(iter,9)
      di=prmter(iter,10)
      hi=prmter(iter,11)

      do jj=1,limit
        
        jatm=lattsf(jj)
        jter=lstter(ltype(jatm))
        if(jter.ge.0.and.iatm.ne.jatm)then
          
          jjter=loc2(iter,jter)
          if(rtf(jj).le.vmbp(1,jjter,1))then
            
            flag=.false.
            
c     potential energy and virial terms
            
            vterm=0.d0
            eterm=0.d0
            
c     initialise work  arrays
            
            do kk=1,limit

              gam(kk)=0.d0
              gvr(kk)=0.d0
              
            enddo

c     calculate bond factor

            do kk=1,limit
              
              katm=lattsf(kk)
              kter=lstter(ltype(katm))
              if(kter.ge.0.and.iatm.ne.katm.and.jatm.ne.katm)then
                
                kkter=loc2(iter,kter)
                
                if(rtf(kk).le.vmbp(1,kkter,1))then
                
                  cost=(xtf(jj)*xtf(kk)+ytf(jj)*ytf(kk)+ztf(jj)*ztf(kk))
                  gtheta=1.d0+(ci/di)**2-ci**2/(di**2+(hi-cost)**2)
                  eterm=eterm+gtheta*prmter2(2,kkter)*scr(kk)
                  vterm=vterm+gtheta*prmter2(2,kkter)*gcr(kk)*rtf(kk)
                  gvr(kk)=2.d0*ci**2*(hi-cost)/(di**2+(hi-cost)**2)**2
                  gam(kk)=gtheta
                  flag=.true.

                endif
                
              endif
              
            enddo
            
            if(flag)then
              
c     tersoff energy and virial

              gam_ij=prmter2(1,jjter)*(1.d0+(bi*eterm)**ei)**(-0.5d0/ei)
              gamma=0.5d0*prmter2(1,jjter)*bi*(bi*eterm)**(ei-1.d0)*
     x          eat(jj)*(1.d0+(bi*eterm)**ei)**(-0.5d0/ei-1.d0)
              engter=engter+ert(jj)-gam_ij*eat(jj)
              virter=virter+gamma*vterm+(grt(jj)-gam_ij*gat(jj))*rtf(jj)

c     calculate 3-body forces
              
              do kk=1,limit
                
                katm=lattsf(kk)
                kter=lstter(ltype(katm))
                if(kter.ge.0.and.iatm.ne.katm.and.jatm.ne.katm)then
                  
                  kkter=loc2(iter,kter)
                  if(rtf(kk).le.vmbp(1,kkter,1))then
                    
                    gam_j=0.5d0*gamma*prmter2(2,kkter)*scr(kk)*gvr(kk)
                    gam_k=0.5d0*gamma*prmter2(2,kkter)*gcr(kk)*gam(kk)
                    
c     calculate contribution to atomic forces
                    
                    cost=(xtf(jj)*xtf(kk)+ytf(jj)*ytf(kk)+
     x                ztf(jj)*ztf(kk))
                    
                    fxa=gam_j*(xtf(kk)-xtf(jj)*cost)/rtf(jj)
                    fya=gam_j*(ytf(kk)-ytf(jj)*cost)/rtf(jj)
                    fza=gam_j*(ztf(kk)-ztf(jj)*cost)/rtf(jj)
                    
                    fxc=gam_j*(xtf(jj)-xtf(kk)*cost)/rtf(kk)-
     x                gam_k*xtf(kk)
                    fyc=gam_j*(ytf(jj)-ytf(kk)*cost)/rtf(kk)-
     x                gam_k*ytf(kk)
                    fzc=gam_j*(ztf(jj)-ztf(kk)*cost)/rtf(kk)-
     x                gam_k*ztf(kk)
                    
                    fxx(jatm)=fxx(jatm)+fxa
                    fyy(jatm)=fyy(jatm)+fya
                    fzz(jatm)=fzz(jatm)+fza
                    
                    fxx(iatm)=fxx(iatm)-(fxa+fxc)
                    fyy(iatm)=fyy(iatm)-(fya+fyc)
                    fzz(iatm)=fzz(iatm)-(fza+fzc)
                    
                    fxx(katm)=fxx(katm)+fxc
                    fyy(katm)=fyy(katm)+fyc
                    fzz(katm)=fzz(katm)+fzc
                    
c     calculate contributions to stress tensor
                    
                    strs(1)=strs(1)+(fxa*xtf(jj)*rtf(jj)+
     x                fxc*xtf(kk)*rtf(kk))
                    strs(2)=strs(2)+(fxa*ytf(jj)*rtf(jj)+
     x                fxc*ytf(kk)*rtf(kk))
                    strs(3)=strs(3)+(fxa*ztf(jj)*rtf(jj)+
     x                fxc*ztf(kk)*rtf(kk))
                    strs(4)=strs(4)+(fya*ytf(jj)*rtf(jj)+
     x                fyc*ytf(kk)*rtf(kk))
                    strs(5)=strs(5)+(fya*ztf(jj)*rtf(jj)+
     x                fyc*ztf(kk)*rtf(kk))
                    strs(6)=strs(6)+(fza*ztf(jj)*rtf(jj)+
     x                fzc*ztf(kk)*rtf(kk))

                  endif
                  
                endif
                
              enddo
            
            else
              
              gam_ij=prmter2(1,jjter)
              engter=engter+ert(jj)-gam_ij*eat(jj)
              virter=virter+(grt(jj)-gam_ij*gat(jj))*rtf(jj)

            endif

c     calculate two body force terms
            
            gterm=0.5d0*(grt(jj)-gam_ij*gat(jj))
            fxx(iatm)=fxx(iatm)+xtf(jj)*gterm
            fyy(iatm)=fyy(iatm)+ytf(jj)*gterm
            fzz(iatm)=fzz(iatm)+ztf(jj)*gterm
            fxx(jatm)=fxx(jatm)-xtf(jj)*gterm
            fyy(jatm)=fyy(jatm)-ytf(jj)*gterm
            fzz(jatm)=fzz(jatm)-ztf(jj)*gterm
            
c     calculate contributions to stress tensor
            
            strs(1)=strs(1)-gterm*rtf(jj)*xtf(jj)*xtf(jj)
            strs(2)=strs(2)-gterm*rtf(jj)*xtf(jj)*ytf(jj)
            strs(3)=strs(3)-gterm*rtf(jj)*xtf(jj)*ztf(jj)
            strs(4)=strs(4)-gterm*rtf(jj)*ytf(jj)*ytf(jj)
            strs(5)=strs(5)-gterm*rtf(jj)*ytf(jj)*ztf(jj)
            strs(6)=strs(6)-gterm*rtf(jj)*ztf(jj)*ztf(jj)

          endif
          
        endif
        
      enddo
      
      return
      end subroutine tersoff3
      
      end module tersoff_module
