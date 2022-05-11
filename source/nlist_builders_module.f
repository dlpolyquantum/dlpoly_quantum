      module nlist_builders_module
      
c***********************************************************************
c     
c     dl_poly module for defining neighbourlist builder routines
c     copyright - daresbury laboratory
c     author    - w. smith    aug 2006
c     pimd adaptation - w. smith jun 2016
c     
c***********************************************************************
      
      use config_module
      use error_module
      use exclude_module
      use pair_module
      use setup_module
      use utility_module
      
      contains
      
      subroutine nlist_driver
     x  (newlst,lneut,lnsq,loglnk,ltad,natms,nbeads,idnode,mxnode,
     x  imcon,nneut,keyfce,rcut,delr,tstep)
      
c*********************************************************************
c     
c     dl_poly subroutine to select and implement neighbour list 
c     builders for short range force calculations
c     
c     copyright - daresbury laboratory
c     author    - w. smith june 2006
c     
c*********************************************************************
      
      implicit none
      
      logical newlst,lneut,lnsq,loglnk,newjob,ltad
      integer natms,idnode,mxnode,imcon,nneut,keyfce,nbeads
      real(8) rcut,delr,tstep
      
      save newjob
      
      data newjob/.true./

c     skip if no pair force calculations required
      
      if(keyfce.gt.0)then
        
c     test for updating the Verlet list
        
        if(ltad)then
          
          call vertest2(newlst,idnode,mxnode,natms,imcon,delr,tstep)
          
        else
          
          call vertest(newlst,idnode,mxnode,natms,nbeads,delr,tstep)
          
        endif
        
c     set up nonbonded interaction (verlet) list
        
        newlst=(newjob.or.newlst)
        
        if(newlst)then
          
          if(.not.lneut)then
            
            if(lnsq)then 
              
c     calculate distant interactions explicitly
              
              call parlst_nsq
     x          (newlst,natms,nbeads,idnode,mxnode,imcon,rcut)
              
            elseif(loglnk)then
              
c     ignore real space distant interactions
              
              call parlink
     x          (newlst,natms,nbeads,idnode,mxnode,imcon,rcut,delr)
              
            else
              
              call parlst
     x          (newlst,natms,nbeads,idnode,mxnode,imcon,rcut,delr)
              
            endif
            
          else
            
            if(.not.loglnk)then 
              
              call parneulst
     x          (newlst,lneut,nneut,idnode,mxnode,imcon,rcut,delr)
              
            else
              
              call parlinkneu
     x          (newlst,natms,nneut,idnode,mxnode,imcon,rcut,delr)
              
            endif
            
          endif
          
        endif
        
      endif
      
      newjob=.false.
      
      return
      end subroutine nlist_driver
      
      subroutine parlst
     x  (newlst,natms,nbeads,idnode,mxnode,imcon,rcut,delr)
      
c***********************************************************************
c     
c     dl_poly subroutine for constructing the verlet neighbour
c     list based on the brode-ahlrichs scheme
c     frozen atoms taken into account
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1992
c     author    - w.smith    march 1992
c     modified  - t.forester october 1993
c     
c***********************************************************************
      
      implicit none
      
      logical lchk,newlst,ldo,lfrzi
      integer natms,nbeads,idnode,mxnode,imcon,ibig,last,nsatm
      integer mpm2,npm2,idum,i,j,k,m,n,ib,jb,nb,kbase,nbase
      real(8) rcut,delr,rclim,rsq
      
      if(newlst)then
        
        ibig=0
        
c     check size of work array
        
        if(mxxdf.lt.(natms+1)/2)then
          
          if(idnode.eq.0) write(nrite,*) 'mxxdf must be at least ',
     x      (natms+1)/2
          call  error(idnode,474)
          
        endif
        
c     set control variables
        
        kbase=0
        nbase=0
        lchk=.true.
        mpm2=natms/2
        npm2=(natms-1)/2
        nsatm=(natms-idnode-1)/mxnode+1
        
c     set cutoff radius
        
        rclim=(rcut+delr)**2
        
c     initialise neighbour list entry count
        
        do i=1,mslist
          
          lentry(i)=0
          
        enddo
        
c     loop over pimd beads
        
        do k=1,nbeads
          
          last=natms
          
c     initialise excluded atom target counter
          
          do i=1,msatms
            
            noxatm(i)=1
            
          enddo
        
c     outer loop over atoms
          
          do m=1,mpm2
            
            if(m.gt.npm2)last=mpm2
            
c     inner loop over atoms
            
            n=0
            
            do i=idnode+1,last,mxnode
              
c     calculate atom indices
              
              j=i+m
              if(j.gt.natms)j=j-natms
              ib=i+kbase
              jb=j+kbase
              
c     calculate interatomic displacements
              
              n=n+1
              xdf(n)=xxx(ib)-xxx(jb)
              ydf(n)=yyy(ib)-yyy(jb)
              zdf(n)=zzz(ib)-zzz(jb)
              
            enddo
            
c     apply minimum image convention
            
            call images(imcon,0,1,n,cell,xdf,ydf,zdf)
            
c     allocate atoms to neighbour list
            
            n=0
            
            do i=idnode+1,last,mxnode
              
c     calculate atom indices
              
              n=n+1
              j=i+m
              if(j.gt.natms)j=j-natms
              ib=i+kbase
              jb=j+kbase
              nb=n+nbase

c     identify frozen atom
              
              lfrzi=(lstfrz(ib).ne.0)
              
c     reject atoms in excluded pair list
              
              if((nexatm(n).gt.0).and.(lexatm(n,noxatm(n)).eq.j))
     x          then
                
                noxatm(n)=min(noxatm(n)+1,nexatm(n))
                
              else
                
c     reject frozen atom pairs
                
                ldo=.true.
                if(lfrzi)ldo=(lstfrz(jb).eq.0)
                
                if(ldo)then
                  
c     calculate interatomic distance
                  
                  if(imcon.eq.6)then
                    
                    rsq=xdf(n)*xdf(n)+ydf(n)*ydf(n)
                    
                  else
                    
                    rsq=xdf(n)*xdf(n)+ydf(n)*ydf(n)+zdf(n)*zdf(n)
                    
                  endif
                  
c     running check of neighbour list array capacity
                  
                  if(rsq.lt.rclim)then
                    
                    lentry(nb)=lentry(nb)+1
                    
                    if(lentry(nb).gt.mxlist)then
                      
                      lchk=.false.
                      ibig=max(lentry(nb),ibig)
                      
                    endif
                    
c     compile neighbour list array
                    
                    if(lchk)then
                      
                      list(nb,lentry(nb))=jb
                      
                    endif
                    
                  endif
                  
                endif
                
              endif
              
            enddo
            
          enddo
          
          nbase=nbase+nsatm
          kbase=kbase+natms
          
c     terminate job if neighbour list array exceeded
          
          if(mxnode.gt.1) call gstate(lchk)
          
          if(.not.lchk)then
            
            call gimax(ibig,1,idum)
            
            if(idnode.eq.0)then
              write(nrite,*) ' mxlist must be at least  ',ibig
              write(nrite,*) ' mxlist is currently ',mxlist
            endif
            
            call error(idnode,110)
            
          endif
          
        enddo
                
c     check all excluded atoms are accounted for
        
        do i=1,nsatm
          
          if(nexatm(i).gt.0.and.noxatm(i).ne.nexatm(i))lchk=.false.
          
        enddo
        
        if(mxnode.gt.1) call gstate(lchk)
        if(.not.lchk) call error(idnode,160)
        
      endif
      
      return
      end subroutine parlst
      
      subroutine parlink
     x  (newlst,natms,nbeads,idnode,mxnode,imcon,rcut,delr)
      
c***********************************************************************
c     
c     dl_poly subroutine for constructing the verlet neighbour
c     list based on link-cell method.
c     frozen atoms taken into account
c     
c     to be used with the link version of exclude :exclude_link
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1993
c     author    - t. forester september 1993.
c     
c***********************************************************************
      
      implicit none
      
      logical lchk,newlst,linc,newjob,lfrzi,ldo
      integer natms,nbeads,idnode,mxnode,imcon,idum,nix,niy,niz
      integer i,ibig,nsbcll,ilx,ily,ilz,ncells,ix,iy,iz,j,icell
      integer ic,n,kc,ik,jx,jy,jz,jc,ixl,numatm,fail,nbase,ibase
      integer nsatm,k,irat
      real(8) rcut,delr,rcsq,xm,ym,zm,det,xdc,ydc,zdc,tx,ty,tz
      real(8) cx,cy,cz,sxd,syd,szd,xd,yd,zd,rsq
      
      real(8), allocatable :: uxx(:),uyy(:),uzz(:)
      
      dimension nix(508),niy(508),niz(508)
      
      save newjob
      data newjob/.true./
      
      data nix/0,1,0,0,-1,1,0,-1,1,0,-1,1,-1,1,2,0,0,-2,2,-1,1,0,-2,2,0,
     x 0,-1,1,0,-1,1,-2,2,-2,2,-1,1,-1,1,-1,1,-2,2,0,-2,2,0,-2,2,-2,2,
     x -1,1,-2,2,-2,2,-1,1,-2,2,-2,2,3,0,0,-3,3,-1,1,0,-3,3,0,0,-1,1,0,
     x -1,1,-3,3,-3,3,-1,1,-1,1,-1,1,-3,3,-2,2,0,-3,3,0,0,-2,2,0,-2,2,
     x -3,3,-3,3,-2,2,-1,1,-3,3,-3,3,-1,1,-1,1,-2,2,-2,2,-1,1,-2,2,-3,3,
     x -3,3,-2,2,-2,2,-2,2,-3,3,0,-3,3,0,-3,3,-3,3,-1,1,-3,3,-3,3,-1,1,
     x -3,3,-3,3,-2,2,-3,3,-3,3,-2,2,-3,3,-3,3,4,0,0,-4,4,-1,1,0,-4,4,0,
     x 0,-1,1,0,-1,1,-4,4,-4,4,-1,1,-1,1,-1,1,-4,4,-2,2,0,-4,4,0,0,-2,2,
     x 0,-2,2,-4,4,-4,4,-2,2,-1,1,-4,4,-4,4,-1,1,-1,1,-2,2,-2,2,-1,1,-2,
     x 2,-4,4,-4,4,-2,2,-2,2,-2,2,-4,4,-3,3,0,-4,4,0,0,-3,3,0,-3,3,-4,4,
     x -4,4,-3,3,-1,1,-4,4,-4,4,-1,1,-1,1,-3,3,-3,3,-1,1,-3,3,-4,4,-4,4,
     x -3,3,-2,2,-4,4,-4,4,-2,2,-2,2,-3,3,-3,3,-2,2,-3,3,-4,4,-4,4,-3,3,
     x -3,3,-3,3,-4,4,0,-4,4,0,-4,4,-4,4,-1,1,-4,4,-4,4,-1,1,-4,4,-4,4,
     x -2,2,-4,4,-4,4,-2,2,-4,4,-4,4,-3,3,-4,4,-4,4,-3,3,5,0,0,-5,5,-1,
     x 1,0,-5,5,0,0,-1,1,0,-1,1,-5,5,-5,5,-1,1,-1,1,-1,1,-5,5,-2,2,0,-5,
     x 5,0,0,-2,2,0,-2,2,-5,5,-5,5,-2,2,-1,1,-5,5,-5,5,-1,1,-1,1,-2,2,
     x -2,2,-1,1,-2,2,-5,5,-5,5,-2,2,-2,2,-2,2,-5,5,-3,3,0,-5,5,0,0,-3,
     x 3,0,-3,3,-5,5,-5,5,-3,3,-1,1,-5,5,-5,5,-1,1,-1,1,-3,3,-3,3,-1,1,
     x -3,3,-5,5,-5,5,-3,3,-2,2,-5,5,-5,5,-2,2,-2,2,-3,3,-3,3,-2,2,-3,3,
     x -5,5,-5,5,-3,3,-3,3,-3,3/
      data niy/  0,0,1,0,1,1,-1,0,0,1,-1,-1,1,1,0,2,0,1,1,2,2,-2,0,0,2,
     x -1,0,0,1,-2,-2,-1,-1,1,1,2,2,-1,-1,1,1,2,2,-2,0,0,2,-2,-2,2,2,-2,
     x -2,-1,-1,1,1,2,2,-2,-2,2,2,0,3,0,1,1,3,3,-3,0,0,3,-1,0,0,1,-3,-3,
     x -1,-1,1,1,3,3,-1,-1,1,1,2,2,3,3,-3,0,0,3,-2,0,0,2,-3,-3,-2,-2,2,
     x 2,3,3,-3,-3,-1,-1,1,1,3,3,-2,-2,-1,-1,1,1,2,2,-3,-3,-2,-2,2,2,3,
     x 3,-2,-2,2,2,3,3,-3,0,0,3,-3,-3,3,3,-3,-3,-1,-1,1,1,3,3,-3,-3,3,3,
     x -3,-3,-2,-2,2,2,3,3,-3,-3,3,3,0,4,0,1,1,4,4,-4,0,0,4,-1,0,0,1,-4,
     x -4,-1,-1,1,1,4,4,-1,-1,1,1,2,2,4,4,-4,0,0,4,-2,0,0,2,-4,-4,-2,-2,
     x 2,2,4,4,-4,-4,-1,-1,1,1,4,4,-2,-2,-1,-1,1,1,2,2,-4,-4,-2,-2,2,2,
     x 4,4,-2,-2,2,2,3,3,4,4,-4,0,0,4,-3,0,0,3,-4,-4,-3,-3,3,3,4,4,-4,
     x -4,-1,-1,1,1,4,4,-3,-3,-1,-1,1,1,3,3,-4,-4,-3,-3,3,3,4,4,-4,-4,
     x -2,-2,2,2,4,4,-3,-3,-2,-2,2,2,3,3,-4,-4,-3,-3,3,3,4,4,-3,-3,3,3,
     x 4,4,-4,0,0,4,-4,-4,4,4,-4,-4,-1,-1,1,1,4,4,-4,-4,4,4,-4,-4,-2,-2,
     x 2,2,4,4,-4,-4,4,4,-4,-4,-3,-3,3,3,4,4,0,5,0,1,1,5,5,-5,0,0,5,-1,
     x 0,0,1,-5,-5,-1,-1,1,1,5,5,-1,-1,1,1,2,2,5,5,-5,0,0,5,-2,0,0,2,-5,
     x -5,-2,-2,2,2,5,5,-5,-5,-1,-1,1,1,5,5,-2,-2,-1,-1,1,1,2,2,-5,-5,
     x -2,-2,2,2,5,5,-2,-2,2,2,3,3,5,5,-5,0,0,5,-3,0,0,3,-5,-5,-3,-3,3,
     x 3,5,5,-5,-5,-1,-1,1,1,5,5,-3,-3,-1,-1,1,1,3,3,-5,-5,-3,-3,3,3,5,
     x 5,-5,-5,-2,-2,2,2,5,5,-3,-3,-2,-2,2,2,3,3,-5,-5,-3,-3,3,3,5,5,-3,
     x -3,3,3/
      data niz/0,0,0,1,0,0,1,1,1,1,1,1,1,1,0,0,2,0,0,0,0,1,1,1,1,2,2,2,
     x 2,1,1,1,1,1,1,1,1,2,2,2,2,0,0,2,2,2,2,1,1,1,1,2,2,2,2,2,2,2,2,2,
     x 2,2,2,0,0,3,0,0,0,0,1,1,1,1,3,3,3,3,1,1,1,1,1,1,1,1,3,3,3,3,0,0,
     x 0,0,2,2,2,2,3,3,3,3,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,3,
     x 3,3,2,2,2,2,2,2,2,2,3,3,3,3,0,0,3,3,3,3,1,1,1,1,3,3,3,3,3,3,3,3,
     x 2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,0,0,4,0,0,0,0,1,1,1,1,4,4,4,4,1,
     x 1,1,1,1,1,1,1,4,4,4,4,0,0,0,0,2,2,2,2,4,4,4,4,1,1,1,1,1,1,1,1,2,
     x 2,2,2,2,2,2,2,4,4,4,4,4,4,4,4,2,2,2,2,2,2,2,2,4,4,4,4,0,0,0,0,3,
     x 3,3,3,4,4,4,4,1,1,1,1,1,1,1,1,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,2,
     x 2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,3,3,3,3,3,3,3,3,4,
     x 4,4,4,0,0,4,4,4,4,1,1,1,1,4,4,4,4,4,4,4,4,2,2,2,2,4,4,4,4,4,4,4,
     x 4,3,3,3,3,4,4,4,4,4,4,4,4,0,0,5,0,0,0,0,1,1,1,1,5,5,5,5,1,1,1,1,
     x 1,1,1,1,5,5,5,5,0,0,0,0,2,2,2,2,5,5,5,5,1,1,1,1,1,1,1,1,2,2,2,2,
     x 2,2,2,2,5,5,5,5,5,5,5,5,2,2,2,2,2,2,2,2,5,5,5,5,0,0,0,0,3,3,3,3,
     x 5,5,5,5,1,1,1,1,1,1,1,1,3,3,3,3,3,3,3,3,5,5,5,5,5,5,5,5,2,2,2,2,
     x 2,2,2,2,3,3,3,3,3,3,3,3,5,5,5,5,5,5,5,5,3,3,3,3,3,3,3,3,5,5,5,5/
      
      data fail/0/
      
      numatm=nbeads*natms
      nsatm=(natms-idnode-1)/mxnode+1
      
      if(newlst)then
        
        if(imcon.eq.4.or.imcon.eq.5.or.imcon.eq.7)
     x    call error(idnode,300)
        lchk=.true.
        ibig=0

c     allocate work arrays

        allocate (uxx(natms),uyy(natms),uzz(natms),stat=fail)
        if(fail.ne.0)call error(idnode,1890)
        
c     real space cut off 
        
        rcsq=(rcut+delr)**2
        
c     create mock cell vector for non-periodic system
        
        if(imcon.eq.0.or.imcon.eq.6)then
          
c     find maximum x,y,z postions
          
          xm=0.d0
          ym=0.d0
          zm=0.d0
          
          do i=1,numatm
            
            xm=max(xm,abs(xxx(i)))
            ym=max(ym,abs(yyy(i)))
            zm=max(zm,abs(zzz(i)))
            
          enddo
          
          if(imcon.eq.0)then
            
            cell(1)=max(2.d0*xm+rcut+delr,3.d0*(rcut+delr))
            cell(5)=max(2.d0*ym+rcut+delr,3.d0*(rcut+delr))
            cell(2)=0.d0
            cell(3)=0.d0
            cell(4)=0.d0
            cell(6)=0.d0
            cell(7)=0.d0
            cell(8)=0.d0
            
          endif
          
          cell(9)=max(2.d0*zm+rcut+delr,3.d0*(rcut+delr),cell(9))
          
        endif
        
        call dcell(cell,celprp)
        call invert(cell,rcell,det)
        
c     number of neighbour subcells
        
        irat=1
        if (irat.eq.1)then 
          nsbcll=14
        elseif(irat.eq.2)then
          nsbcll=63
        elseif(irat.eq.3)then
          nsbcll=156
        elseif(irat.eq.4)then
          nsbcll=307
        elseif(irat.eq.5)then
          nsbcll=508
        endif
        
c     link cell numbers in principal directions

        ilx=int(celprp(7)*dble(irat)/(rcut+delr))
        ily=int(celprp(8)*dble(irat)/(rcut+delr))
        ilz=int(celprp(9)*dble(irat)/(rcut+delr))
        
c     check there are enough link cells
        
        linc=.true.
        if(ilx.lt.2*irat+1)linc=.false.
        if(ily.lt.2*irat+1)linc=.false.
        if(ilz.lt.2*irat+1)linc=.false.
        if(.not.linc) call error(idnode,305)
        ncells=ilx*ily*ilz
        if(ncells.gt.mxcell) call error(idnode,392)
        
c     link-cell cutoff for reduced space
        
        xdc=dble(ilx)
        ydc=dble(ily)
        zdc=dble(ilz)
        
c     ensure all atoms in MD cell on first pass
        
        if(newjob)then
          
          newjob=.false.
          call images(imcon,idnode,mxnode,numatm,cell,xxx,yyy,zzz)
          
          if(mxnode.gt.1)  call merge
     x      (idnode,mxnode,numatm,mxbuff,xxx,yyy,zzz,buffer)
          
        endif

c     loop over pimd beads

        do k=1,nbeads

          ibase=(k-1)*natms
          nbase=(k-1)*nsatm
          
c     zero last entry in neighbour list
          
          do i=1,mslist
            lentry(i)=0
          enddo
          
c     zero head of link cell chain
          
          do i=1,ncells
            lct(i)=0
          enddo
          
c     zero link array
          
          do i=1,natms
            link(i)=0
          enddo
          
c     calculate reduced space coordinates
          
          do i=1,natms
            
            tx=xxx(i+ibase)
            ty=yyy(i+ibase)
            tz=zzz(i+ibase)
            
            uxx(i)=(rcell(1)*tx+rcell(4)*ty+rcell(7)*tz)+0.5d0
            uyy(i)=(rcell(2)*tx+rcell(5)*ty+rcell(8)*tz)+0.5d0
            uzz(i)=(rcell(3)*tx+rcell(6)*ty+rcell(9)*tz)+0.5d0
            
          enddo
          
c     link neighbours 
          
          do i=1,natms
            
            ix=min(int(xdc*uxx(i)),ilx-1)
            iy=min(int(ydc*uyy(i)),ily-1)
            iz=min(int(zdc*uzz(i)),ilz-1)
            
            icell=1+ix+ilx*(iy+ily*iz)
            
            j=lct(icell)
            lct(icell)=i
            link(i)=j
            
          enddo
          
c     set control variables for loop over subcells
          
          ix=1
          iy=1
          iz=1
          
c     primary loop over subcells
          
          do ic=1,ncells
            
            n=lct(ic)
            if(n.gt.0)then
              
c     secondary loop over subcells
              
              ik=0
              
              do kc=1,nsbcll
                
                i=n
                
                cx=0.d0
                cy=0.d0
                cz=0.d0
                jx=ix+nix(kc)
                jy=iy+niy(kc)
                jz=iz+niz(kc)
                
c     minimum image convention
                
                if(jx.gt.ilx)then
                  
                  jx=jx-ilx
                  cx=1.d0
                  
                elseif(jx.lt.1)then
                  
                  jx=jx+ilx
                  cx=-1.d0
                  
                endif
                
                if(jy.gt.ily)then
                  
                  jy=jy-ily
                  cy=1.d0
                  
                elseif(jy.lt.1)then
                  
                  jy=jy+ily
                  cy=-1.d0
                  
                endif
                
                if(jz.gt.ilz)then
                  
                  jz=jz-ilz
                  cz=1.d0
                  
                elseif(jz.lt.1)then
                  
                  jz=jz+ilz
                  cz=-1.d0
                  
                endif
                
c     index of neighbouring cell
                
                jc=jx+ilx*((jy-1)+ily*(jz-1))
                j=lct(jc)
                
c     ignore if empty
                
                if(j.gt.0)then
                  
                  do while(i.ne.0)
                    
c     test if site is of interest to this node
                    
                    if(mod(i-1,mxnode).eq.idnode)then
                      
c     i's index for this processor
                      
                      ik=(i-1)/mxnode+nbase+1
                      
c     test if i is a frozen atom
                      
                      lfrzi=(lstfrz(i).ne.0)
                      
                      if(ic.eq.jc)j=link(i)
                      if(j.gt.0)then
                        
                        do while(j.ne.0)
                          
c     test of frozen atom pairs
                          
                          ldo=.true.
                          if(lfrzi)ldo=(lstfrz(j).eq.0)
                          
                          if(ldo)then
                            
c     distance in real space : minimum image applied
                            
                            sxd=uxx(j)-uxx(i)+cx
                            syd=uyy(j)-uyy(i)+cy
                            szd=uzz(j)-uzz(i)+cz
                            
                            xd=cell(1)*sxd+cell(4)*syd+cell(7)*szd
                            yd=cell(2)*sxd+cell(5)*syd+cell(8)*szd
                            zd=cell(3)*sxd+cell(6)*syd+cell(9)*szd
                            
                            if(imcon.eq.6)then
                              
                              rsq=xd**2+yd**2
                              
                            else
                              
                              rsq=xd**2+yd**2+zd**2
                              
                            endif
                            
c     cut-off test
                            
                            if(rcsq.gt.rsq)then
                              
c     excluded atom test
                              
                              linc=.true.
                              do ixl=1,nexatm(ik)
                                
                                if(lexatm(ik,ixl).eq.j) linc=.false.
                                
                              enddo
                              
                              if(linc)then
                                
                                lentry(ik)=lentry(ik)+1
                                
                                if(lentry(ik).gt.mxlist)then
                                  
                                  ibig=max(ibig,lentry(ik))
                                  lchk=.false.
                                  
                                else
                                  
                                  list(ik,lentry(ik))=j
                                  
                                endif
                                
                              endif
                              
                            endif
                            
                          endif
                          
                          j=link(j)
                          
                        enddo
                        
                      endif
                      
                    endif
                    
                    j=lct(jc)
                    i=link(i)
                    
                  enddo
                  
                endif
                
              enddo
              
            endif
            
            ix=ix+1
            if(ix.gt.ilx)then
              
              ix=1
              iy=iy+1
              
              if(iy.gt.ily)then
                
                iy=1
                iz=iz+1
                
              endif
              
            endif
            
          enddo
          
c     terminate job if neighbour list array exceeded
          
          if(mxnode.gt.1) call gstate(lchk)
          
          if(.not.lchk)then
            
            call gimax(ibig,1,idum)
            if(idnode.eq.0)then
              write(nrite,*) ' mxlist must be >=  ',ibig
              write(nrite,*) ' mxlist is currenty ',mxlist
            endif
            call error(idnode,106)
            
          endif

c     end of loop over pimd beads
          
        enddo
        
c     deallocate work arrays
        
        deallocate (uxx,uyy,uzz,stat=fail)
        
      endif
      
      return
      end subroutine parlink
      
      subroutine parneulst
     x  (newlst,lneut,nneut,idnode,mxnode,imcon,rcut,delr)
      
c***********************************************************************
c     
c     dl_poly subroutine for constructing the verlet neighbour
c     list based on the brode-ahlrichs scheme
c     frozen atoms taken into account
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1992
c     author    - t.forester   april   1994
c     
c***********************************************************************
      
      implicit none
      
      logical lchk,newlst,lneut,safe,lfrzi
      integer nneut,idnode,mxnode,imcon,idum,fail,mpm2,npm2,ibig
      integer i,ill,ia,im,jmlast,jmwrap,nuei1,nuei2,ii,jm1,jm,jj0,jj2
      integer j,ii1
      real(8) rcut,delr,rclim,fi,rrr,rcl1
      
      logical, allocatable :: lms(:)
      dimension fi(3)
      
      data fail/0/
      
C     DIR$ CACHE_ALIGN fi
      
c     allocate work arrays
      
      allocate (lms(mxneut),stat=fail)
      if(fail.ne.0)call error(idnode,1910)
      
      if(newlst.and.lneut)then
        
c     set control variables
        
        safe=.true. 
        lchk= .true.
        mpm2=(nneut+2)/2
        npm2=(nneut+1)/2 
        
c     set cutoff radius
        
        rcl1=(rcut+delr)
        rclim=(rcut+delr)**2
        ibig=0
        ill=0
        
c     construct pair force neighbour list: neutral groups
        
        do i=1,mslist
          
          lentry(i)=0
          
        enddo
        
c     outer loop over groups
        
        ia=0
        
        do im=idnode+1,nneut,mxnode
          
          ia=ia+1 
          if(im.ge.mpm2) mpm2=npm2
          
          lms(1)=.false.
          do j=2,mpm2
            lms(j)=.true.
          enddo
          
          jmlast=min(nneut,im+mpm2-1)
          jmwrap=max(0,im+mpm2-1-nneut)
          
c     loop over atomic pairs
          
          nuei1=neulst(im)
          nuei2=neulst(im+1)-1
          
          do i=nuei1,nuei2
            
            fi(1)=xxx(i)
            fi(2)=yyy(i)
            fi(3)=zzz(i)
            lfrzi=(lstfrz(i).eq.0)
            
            ii=0
            jm1=1
            do jm=im+1,jmlast
              
              jm1=jm1+1
              if(lms(jm1))then
                
                jj0=neulst(jm)
                jj2=neulst(jm+1)-1
                
                do j=jj0,jj2
                  
                  ii=ii+1
                  if(ii.le.mxxdf)then
                    
                    xdf(ii)=fi(1)-xxx(j)
                    ydf(ii)=fi(2)-yyy(j)
                    zdf(ii)=fi(3)-zzz(j)
                    
                  else
                    
                    ibig=max(ibig,ii)
                    safe=.false.
                    
                  endif
                  
                enddo
                
              endif
              
            enddo
            
            do jm=1,jmwrap
              
              jm1=jm1+1
              if(lms(jm1))then
                
                jj0=neulst(jm)
                jj2=neulst(jm+1)-1
                
                do j=jj0,jj2
                  
                  ii=ii+1
                  if(ii.le.mxxdf)then
                    
                    xdf(ii)=fi(1)-xxx(j)
                    ydf(ii)=fi(2)-yyy(j)
                    zdf(ii)=fi(3)-zzz(j)
                    
                  else
                    
                    safe=.false.
                    ibig=max(ibig,ii)
                    
                  endif
                  
                enddo
                
              endif
              
            enddo
            
c     apply minimum image convention
            
            ii1=min(ii,mxxdf)
            call images(imcon,0,1,ii1,cell,xdf,ydf,zdf)
            
c     search for those in cutoff
            
            ii=0
            jm1=1
            do jm=im+1,jmlast
              
              jm1=jm1+1
              if(lms(jm1))then
                
                jj0=neulst(jm)
                jj2=neulst(jm+1)-1
                
                do j=jj0,jj2
                  
                  ii=ii+1
                  if(ii.le.mxxdf)then
                    
                    if(lms(jm1))then
                      
                      if(lfrzi)then
                        if(abs(zdf(ii)).lt.rcl1)then
                          if(abs(ydf(ii)).lt.rcl1)then
                            if(abs(xdf(ii)).lt.rcl1)then
                              rrr=xdf(ii)**2+ydf(ii)**2+zdf(ii)**2
                              if(rrr.lt.rclim) lms(jm1)=.false.
                            endif
                          endif
                        endif
                        
                      elseif(lstfrz(j).eq.0)then 
                        if(abs(zdf(ii)).lt.rcl1)then
                          if(abs(ydf(ii)).lt.rcl1)then
                            if(abs(xdf(ii)).lt.rcl1)then
                              rrr=xdf(ii)**2+ydf(ii)**2+zdf(ii)**2
                              if(rrr.lt.rclim) lms(jm1)=.false.
                            endif
                          endif
                        endif
                        
                      endif
                      
                    endif
                    
                  endif
                  
                enddo
                
              endif
              
            enddo
            
            do jm=1,jmwrap
              
              jm1=jm1+1
              if(lms(jm1))then
                
                jj0=neulst(jm)
                jj2=neulst(jm+1)-1
                
                do j=jj0,jj2
                  
                  ii=ii+1
                  if(ii.le.mxxdf)then
                    
                    if(lms(jm1))then
                      
                      if(lfrzi)then
                        
                        if(abs(zdf(ii)).lt.rcl1)then
                          if(abs(ydf(ii)).lt.rcl1)then
                            if(abs(xdf(ii)).lt.rcl1)then
                              rrr=xdf(ii)**2+ydf(ii)**2+zdf(ii)**2
                              if(rrr.lt.rclim) lms(jm1)=.false.
                            endif
                          endif
                        endif
                        
                      elseif(lstfrz(j).eq.0)then
                        
                        if(abs(zdf(ii)).lt.rcl1)then
                          if(abs(ydf(ii)).lt.rcl1)then
                            if(abs(xdf(ii)).lt.rcl1)then
                              rrr=xdf(ii)**2+ydf(ii)**2+zdf(ii)**2
                              if(rrr.lt.rclim) lms(jm1)=.false.
                            endif
                          endif
                        endif
                        
                      endif
                      
                    endif
                    
                  endif
                  
                enddo
                
              endif
              
            enddo
            
          enddo
          
c     compile neighbour list for ia
c     with running check of neighbour list array capacity
          
          jm1=0
          do jm=im,jmlast
            
            jm1=jm1+1
            if(.not.lms(jm1))then
              
              lentry(ia)=lentry(ia)+1
              if(lentry(ia).le.mxlist)then
                
                list(ia,lentry(ia))=jm
                
              else
                
                ill=max(ill,lentry(ia))
                lchk=.false.
                
              endif
              
            endif
            
          enddo
          
          do jm=1,jmwrap
            
            jm1=jm1+1
            if(.not.lms(jm1))then
              
              lentry(ia)=lentry(ia)+1
              if(lentry(ia).le.mxlist)then
                
                list(ia,lentry(ia))=jm
                
              else
                
                ill=max(ill,lentry(ia))
                lchk=.false.
                
              endif
              
            endif
            
          enddo
          
        enddo
        
c     terminate job if neighbour list array exceeded
        
        if(mxnode.gt.1) call gstate(lchk)
        if(.not.lchk)then
          
          call gimax(ill,1,idum)
          if(idnode.eq.0)then
            write(nrite,*) ' mxlist must be at least  ',ill
            write(nrite,*) ' mxlist is currently ',mxlist
          endif
          call error(idnode,108)
          
        endif   
        
c     terminate job if work arrays exceeded
        
        if(mxnode.gt.1) call gstate(safe)
        if(.not.safe)then
          call gimax(ibig,1,idum)
          if(idnode.eq.0)then
            write(nrite,*)'mxxdf must be at least ',ibig
            write(nrite,*)'mxxdf is currently ',mxxdf
          endif
          call  error(idnode,476)
        endif
        
      endif
      
c     deallocate work arrays
      
      deallocate(lms,stat=fail)
      
      return
      end subroutine parneulst
      
      subroutine parlinkneu
     x  (newlst,natms,nneut,idnode,mxnode,imcon,rcut,delr)
      
c***********************************************************************
c     
c     dl_poly subroutine for constructing the verlet neighbour
c     list based on link-cell method with neutral groups
c     frozen atoms taken into account
c     
c     to be used with the link version of exclude :excludeneu_link
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1996
c     author    - t. forester january 1996.
c     
c***********************************************************************
      
      implicit none
      
      logical lchk,newlst,linc,newjob,lfrzi,ldo,swop
      integer natms,nneut,idnode,mxnode,imcon,idum,fail,ibig
      integer nix,niy,niz,i,nsbcll,ilx,ily,ilz,ncells,ix,iy,iz
      integer icell,j,ic,ii,kc,jx,jy,jz,jc,ineu,ik,jneu,ineua,jneua
      integer ika,jneua1,i1,j1
      real(8) rcut,delr,rcsq,xm,ym,zm,det,xdc,ydc,zdc,tx,ty,tz
      real(8) cx,cy,cz,sxd,syd,szd,rsq,xd,yd,zd
      
      real(8), allocatable :: uxx(:),uyy(:),uzz(:)
      
      dimension nix(14),niy(14),niz(14)
      
      save newjob
      
      data newjob/.true./
      data nix/0,1,0,0,-1,1,0,-1,1,0,-1,1,-1,1/
      data niy/ 0,0,1,0,1,1,-1,0,0,1,-1,-1,1,1/
      data niz/0,0,0,1,0,0,1,1,1,1,1,1,1,1/
      data fail/0/
      
c     allocate work arrays
      
      allocate (uxx(mxatms),uyy(mxatms),uzz(mxatms),stat=fail)
      if(fail.ne.0)call error(idnode,1900)
      
      lchk=.true.
      ibig=0
      if(newlst)then
        
        if(imcon.eq.4.or.imcon.eq.5.or.imcon.eq.7)
     x    call error(idnode,300)
        
c     zero link arrays
        
        do i=1,natms
          
          link(i)=0
          
        enddo
        
c     construct pair force neighbour list
        
        do i=1,mslist
          
          lentry(i)=0
          
        enddo
        
c     real space cut off 
        
        rcsq=(rcut+delr)**2
c     
c     create mock cell vector for non-periodic system
        
        if(imcon.eq.0.or.imcon.eq.6)then
          
c     find maximum x,y,z postions
          
          xm=0.d0
          ym=0.d0
          zm=0.d0
          
          do i=1,natms
            
            xm=max(xm,abs(xxx(i)))
            ym=max(ym,abs(yyy(i)))
            zm=max(zm,abs(zzz(i)))
            
          enddo
          
          if(imcon.eq.0)then
            
            cell(1)=max(2.d0*xm+rcut+delr,3.d0*(rcut+delr))
            cell(5)=max(2.d0*ym+rcut+delr,3.d0*(rcut+delr))
            cell(2)=0.d0
            cell(3)=0.d0
            cell(4)=0.d0
            cell(6)=0.d0
            cell(7)=0.d0
            cell(8)=0.d0
            
          endif
          
          cell(9)=max(2.d0*zm+rcut+delr,3.d0*(rcut+delr),cell(9))
          
        endif
        
        call dcell(cell,celprp)
        call invert(cell,rcell,det)
        
c     number of subcells
        
        nsbcll=14
        
        ilx=int(celprp(7)/(rcut+delr))
        ily=int(celprp(8)/(rcut+delr))
        ilz=int(celprp(9)/(rcut+delr))
c     
c     check there are enough link cells
        
        linc=.false.
        if(ilx.lt.3)linc=.true.
        if(ily.lt.3)linc=.true.
        if(ilz.lt.3)linc=.true.
        if(linc) call error(idnode,305)
        
        ncells=ilx*ily*ilz
        if(ncells.gt.mxcell)then
          
          if(idnode.eq.0) write(nrite,*) 'mxcell must be >= ',ncells
          call  error(idnode,392)
          
        endif
        
c     calculate link cell indices
        
        do i=1,ncells
          
          lct(i)=0
          
        enddo
        
c     link-cell cutoff for reduced space
        
        xdc=dble(ilx)
        ydc=dble(ily)
        zdc=dble(ilz)
        
c     reduced space coordinates
        
        if(newjob)then
          
          newjob=.false.
          call images(imcon,idnode,mxnode,natms,cell,xxx,yyy,zzz)
          
          if(mxnode.gt.1)  call merge
     x      (idnode,mxnode,natms,mxbuff,xxx,yyy,zzz,buffer)
          
        endif
        
        do i=1,natms
          
          tx=xxx(i)
          ty=yyy(i)
          tz=zzz(i)
          
          uxx(i)=(rcell(1)*tx+rcell(4)*ty+rcell(7)*tz)+0.5d0
          uyy(i)=(rcell(2)*tx+rcell(5)*ty+rcell(8)*tz)+0.5d0
          uzz(i)=(rcell(3)*tx+rcell(6)*ty+rcell(9)*tz)+0.5d0
          
        enddo
        
c     link neighbours 
        
        do i=1,natms
          
          ix=min(int(xdc*uxx(i)),ilx-1)
          iy=min(int(ydc*uyy(i)),ily-1)
          iz=min(int(zdc*uzz(i)),ilz-1)
          
          icell=1+ix+ilx*(iy+ily*iz)
          
          j=lct(icell)
          lct(icell)=i
          link(i)=j
          
        enddo
        
c     set control variables for loop over subcells
        
        ix=1
        iy=1
        iz=1
        
c     primary loop over subcells
        
        do ic=1,ncells
          
          ii=lct(ic)
          if(ii.gt.0)then
            
c     secondary loop over subcells
            
            do kc=1,nsbcll
              
              i=ii
              
              cx=0.d0
              cy=0.d0
              cz=0.d0
              jx=ix+nix(kc)
              jy=iy+niy(kc)
              jz=iz+niz(kc)
              
c     minimum image convention
              
              if(jx.gt.ilx)then
                
                jx=jx-ilx
                cx=1.d0
                
              elseif(jx.lt.1)then
                
                jx=jx+ilx
                cx=-1.d0
                
              endif
              
              if(jy.gt.ily)then
                
                jy=jy-ily
                cy=1.d0
                
              elseif(jy.lt.1)then
                
                jy=jy+ily
                cy=-1.d0
                
              endif
              
              if(jz.gt.ilz)then
                
                jz=jz-ilz
                cz=1.d0
                
              elseif(jz.lt.1)then
                
                jz=jz+ilz
                cz=-1.d0
                
              endif
              
c     index of neighbouring cell
              
              jc=jx+ilx*((jy-1)+ily*(jz-1))
              j=lct(jc)
              
c     ignore if empty
              
              if(j.gt.0)then
                
                do while(i.ne.0)
                  
c     test if site is of interest to this node
                  
                  ineu=lstneu(i)
                  ik=0
                  
c     i's  group index for this processor
                  
                  if(mod(ineu-1,mxnode).eq.idnode)
     x              ik=((ineu-1)/mxnode)+1
                  
c     test if i is a frozen atom
                  
                  lfrzi=(lstfrz(i).ne.0)
                  
                  if(ic.eq.jc) j=link(i)
                  if(j.gt.0)then
                    
                    do while(j.ne.0)
                      
                      jneu=lstneu(j)
                      
c     swop tests for switching of group indices,
c     ldo for 'doing' interaction
                      
                      swop=.false.
                      ldo=(ik.gt.0)
                      jneua=jneu
                      ineua=ineu
                      ika=ik
                      
c     keep only Brode-Ahlrichs pairs
                      
                      if(jneua.ge.ineua)then
                        
                        if(jneua-ineua.gt.nneut/2)then 
                          
                          swop=(mod(jneu-1,mxnode).eq.idnode)
                          if(swop)then 
                            ldo=((nneut+ineua-jneua).le.(nneut-1)/2)
                          else
                            ldo=.false.
                          endif
                          
                        endif
                        
                      elseif(nneut+jneua-ineua.gt.(nneut-1)/2)then
                        
                        swop=(mod(jneu-1,mxnode).eq.idnode)
                        if(swop)then
                          ldo=((ineua-jneua).le.nneut/2)
                        else
                          ldo=.false.
                        endif
                        
                      endif
                      
                      if(swop.and.ldo)then
                        jneua=ineu
                        ineua=jneu
                        ika=((jneu-1)/mxnode)+1
                      endif
                      
c     test of frozen atom pairs
                      
                      if(lfrzi.and.ldo)ldo=(lstfrz(j).eq.0)
                      
c     check we haven't already included this group in the list ...
                      
                      jneua1=0
                      do while
     x                  (ldo.and.jneua1.lt.min(lentry(ika),mxlist))
                        
                        jneua1=jneua1+1
                        if(list(ika,jneua1).eq.jneua) ldo=.false.
                        
                      enddo
                      
                      if(ldo)then
                        
c     distance in real space : minimum image applied
                        
                        sxd=uxx(j)-uxx(i)+cx
                        syd=uyy(j)-uyy(i)+cy
                        szd=uzz(j)-uzz(i)+cz
                        
                        xd=cell(1)*sxd+cell(4)*syd+cell(7)*szd
                        yd=cell(2)*sxd+cell(5)*syd+cell(8)*szd
                        zd=cell(3)*sxd+cell(6)*syd+cell(9)*szd
                        
                        rsq=xd*xd+yd*yd+zd*zd
                        
c     test of distance
                        
                        if(rsq.lt.rcsq)then
                          
                          lentry(ika)=lentry(ika)+1
                          if(lentry(ika).gt.mxlist)then
                            
                            ibig=max(ibig,lentry(ika))
                            lchk=.false.
                            
                          else
                            
                            list(ika,lentry(ika))=jneua
                            
                          endif
                          
                        endif
                        
                      endif
                      
                      j=link(j)
                      
                    enddo
                    
                  endif
                  
                  j=lct(jc)
                  i=link(i)
                  
                enddo
                
              endif
              
            enddo
            
          endif
          
          ix=ix+1
          if(ix.gt.ilx)then
            
            ix=1
            iy=iy+1
            
            if(iy.gt.ily)then
              
              iy=1
              iz=iz+1
              
            endif
            
          endif
          
        enddo
        
c     terminate job if neighbour list array exceeded
        
        if(mxnode.gt.1) call gstate(lchk)
        
        if(.not.lchk)then
          
          call gimax(ibig,1,idum)
          if(idnode.eq.0)then
            write(nrite,*)'mxlist must be at least ',ibig
            write(nrite,*)'mxlist is currently ',mxlist
          endif
          call error(idnode,107)
          
        endif
        
c     sort list into order ..
c     use link as a work array
        
        ik=0
        do i=1+idnode,nneut,mxnode
          
          ik=ik+1
          do j=1,lentry(ik)
            
            link(j)=list(ik,j)
            
          enddo
          call shellsort(lentry(ik),link)
          
c     ensure Brode-Ahlrichs ordering
          
          i1=lentry(ik)+1
          j1=0
          do j=1,lentry(ik)
            
            if(link(j).ge.i)then
              
              j1=j1+1
              list(ik,j1)=link(j)
              i1=min(i1,j)
              
            endif
            
          enddo
          
          do j=1,i1-1
            
            j1=j1+1
            list(ik,j1)=link(j)
            
          enddo
          
        enddo
        
      endif
      
c     deallocate work arrays 
      
      deallocate (uxx,uyy,uzz,stat=fail)
      
      return
      end subroutine parlinkneu
      
      subroutine parlst_nsq
     x  (newlst,natms,nbeads,idnode,mxnode,imcon,rcut)
      
c***********************************************************************
c     
c     dl_poly subroutine for constructing the verlet neighbour
c     list based on the brode-ahlrichs scheme
c     frozen atom option included
c     
c     to be used with multiple_nsq
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1994
c     author    - t. forester march 1994     
c     adapted   - w.smith aug 2008 - solvation, excitation etc
c     
c     stress tensor : t.forester may 1994
c     
c***********************************************************************
      
      implicit none
      
      logical lchk,newlst
      integer natms,nbeads,idnode,mxnode,imcon,ibig,last,nsatm
      integer mpm2,npm2,idum,i,j,k,m,n,ib,jb,nb,kbase,nbase
      real(8) rcut,rclim,rsq

      if(newlst)then
        
        ibig=0
        
c     check size of work array
        
        if(mxxdf.lt.(natms+1)/2)then
          
          if(idnode.eq.0) write(nrite,*) 'mxxdf must be greater than ',
     x      (natms+1)/2
          call  error(idnode,475)
          
        endif
        
c     set control variables
        
        kbase=0
        nbase=0
        lchk=.true.
        mpm2=natms/2
        npm2=(natms-1)/2
        nsatm=(natms-idnode-1)/mxnode+1
        
c     set cutoff radius - ignore border width
        
        rclim=rcut**2
        
c     construct pair force neighbour list
        
        do i=1,mslist
          
          lentry(i)=0
          
        enddo
        
c     loop over pimd beads

        do k=1,nbeads
          
          last=natms

c     initialise excluded atom target counter

          do i=1,msatms
            
            noxatm(i)=1

          enddo

c     outer loop over atoms
          
          do m=1,mpm2
            
            if(m.gt.npm2)last=mpm2
            
c     inner loop over atoms
            
            n=0
            
            do i=idnode+1,last,mxnode
              
c     calculate atom indices
              
              j=i+m
              if(j.gt.natms)j=j-natms
              ib=i+kbase
              jb=j+kbase
              
c     calculate interatomic displacements
              
              n=n+1
              xdf(n)=xxx(ib)-xxx(jb)
              ydf(n)=yyy(ib)-yyy(jb)
              zdf(n)=zzz(ib)-zzz(jb)
              
            enddo
            
c     apply minimum image convention
            
            call images(imcon,0,1,n,cell,xdf,ydf,zdf)
            
c     allocate atoms to neighbour list
            
            n=0
            
            do i=idnode+1,last,mxnode
              
c     calculate atom indices
              
              n=n+1
              j=i+m
              if(j.gt.natms)j=j-natms
              ib=i+kbase
              jb=j+kbase
              nb=n+nbase
              
c     reject atoms in excluded pair list
              
              if((nexatm(n).gt.0).and.(lexatm(n,noxatm(n)).eq.j))
     x          then
                
                noxatm(n)=min(noxatm(n)+1,nexatm(n))
                
c     reject frozen atom pairs
                
              elseif(lstfrz(ib).eq.0.or.lstfrz(jb).eq.0)then
                
c     calculate interatomic distance
                
                if(imcon.eq.6)then
                  
                  rsq=xdf(n)**2+ydf(n)**2
                  
                else
                  
                  rsq=xdf(n)**2+ydf(n)**2+zdf(n)**2
                  
                endif
                
c     running check of neighbour list array capacity
                
                if(rsq.lt.rclim)then
                  
                  lentry(nb)=lentry(nb)+1
                  
                  if(lentry(nb).gt.mxlist)then
                    
                    ibig=max(ibig,lentry(nb))
                    lchk=.false.
                    
                  endif
                  
c     compile neighbour list array
                  
                  if(lchk)then
                    
                    list(nb,lentry(nb))=jb
                    
                  endif
                  
                endif
                
              endif
              
            enddo
            
          enddo
          
          nbase=nbase+nsatm
          kbase=kbase+natms
          
c     terminate job if neighbour list array exceeded
          
          if(mxnode.gt.1)call gstate(lchk)
          
          if(.not.lchk)then
            
            call gimax(ibig,1,idum)
            
            if(idnode.eq.0)then
              write(nrite,*) ' mxlist must be >=  ',ibig
              write(nrite,*) ' mxlist is currenty ',mxlist
            endif
            
            call error(idnode,109)
            
          endif
          
        enddo
        
c     check all excluded atoms are accounted for
        
        do i=1,nsatm
          
          if(nexatm(i).gt.0.and.noxatm(i).ne.nexatm(i))lchk=.false.
          
        enddo
        
        if(mxnode.gt.1)call gstate(lchk)
        if(.not.lchk) call error(idnode,160)
        
      endif  
      
      return
      end subroutine parlst_nsq
      
      subroutine primlst(idnode,mxnode,natms,imcon,rprim)
      
c*********************************************************************
c     
c     dlpoly routine to split interaction list into primary and 
c     secondary neighbours for use with multiple timestep method
c     
c     copyright daresbury laboratory
c     author - t. forester february 1993
c     
c*********************************************************************
      
      implicit none
      
      integer idnode,mxnode,natms,imcon,i,j,k,n
      real(8) rprim,rprim2,rsq
      
      rprim2=rprim*rprim

      n=0
      
      do i=1+idnode,natms,mxnode
        
        n=n+1
        
        do k=1,lentry(n)
          
          j=iabs(list(n,k))
          xdf(k)=xxx(i)-xxx(j)
          ydf(k)=yyy(i)-yyy(j)
          zdf(k)=zzz(i)-zzz(j)
          
        enddo           
        
c     apply minimum image convention
        
        call images(imcon,0,1,lentry(n),cell,xdf,ydf,zdf)
        
c     assign atoms as primary or secondary
        
        do j=1,lentry(n)
          
c     calculate interatomic distance

          if(imcon.eq.6)then
          
            rsq=xdf(j)**2+ydf(j)**2

          else

            rsq=xdf(j)**2+ydf(j)**2+zdf(j)**2

          endif
          
          if(rsq.lt.rprim2)then
            
c     compile primary neighbour list array  : -ve indices
            
            list(n,j)=-iabs(list(n,j))
            
          else
            
c     compile secondary neighbour list array : +ve indices
            
            list(n,j)=iabs(list(n,j))
            
          endif
          
        enddo
        
      enddo
      
      return
      end subroutine primlst
      
      subroutine prneulst(newlst,imcon,idnode,mxnode,nneut,rprim)
      
c***********************************************************************
c     
c     dlpoly routine to partition neutral group list into 
c     primary and secondary groups
c     loops over group ineu
c     
c     replicated data version
c     
c     copyright daresbury laboratory 1994
c     author t.forester april 1994
c     
c***********************************************************************
      
      implicit none
      
      logical newlst,lchk,ldo
      integer imcon,idnode,mxnode,nneut,ineu,ia,jj,ibig,n
      integer jj0,jneu,j,i,idum
      real(8) rprim,rclim,xi,yi,zi,rrr
      
      lchk=.true.
      
      if(newlst)then
        
c     set primary cutoff limit
        
        rclim=rprim*rprim
        
c     set list to negative - signal for seconary shell
        
        ia=0
        do ineu=idnode+1,nneut,mxnode
          
          ia=ia+1
          
          do jj=1,lentry(ia)
            
            list(ia,jj)=-abs(list(ia,jj))
            
          enddo
          
        enddo
        
c     loop over neutral group ineu sites
        
        lchk=.true.
        ibig=0
        
        ia=0
        do ineu=idnode+1,nneut,mxnode
          
          ia=ia+1
          
          n=0
          do i=neulst(ineu),neulst(ineu+1)-1
            
            xi=xxx(i)
            yi=yyy(i)
            zi=zzz(i)
            
            do jj=1,lentry(ia)
              
              jneu=-list(ia,jj)
              jj0=neulst(jneu)
              
              if(ineu.eq.jneu) jj0=i+1
              
c     loop over jneu sites
              
              do j=jj0,neulst(jneu+1)-1
                
                n=n+1
                if(n.le.mxxdf)then
                  xdf(n)=xi-xxx(j)
                  ydf(n)=yi-yyy(j)
                  zdf(n)=zi-zzz(j)
                else
                  lchk=.false.
                  ibig=n
                endif
                
              enddo
              
            enddo
            
          enddo
          
c     apply minimum image convention
          
          n=min(n,mxxdf)
          call images(imcon,0,1,n,cell,xdf,ydf,zdf)
          
c     allocate groups to primary or secondary shell
c     on basis of closest atom-atom interactions
          
          n=0
          do i=neulst(ineu),neulst(ineu+1)-1
            
            do jj=1,lentry(ia)
              
              jneu=list(ia,jj)
              ldo=(jneu.lt.0)
              if(ldo)then
                jneu=-jneu
                jj0=neulst(jneu)
                if(ineu.eq.jneu) jj0=i+1
                
c     loop over jneu sites
                
                do j=jj0,neulst(jneu+1)-1
                  
                  if(ldo)then 
                    
                    n=min(n+1,mxxdf)
                    
                    if(abs(xdf(n)).lt.rprim)then
                      if(abs(ydf(n)).lt.rprim)then
                        if(abs(zdf(n)).lt.rprim)then
                          
c     calculate interatomic distance
                          
                          rrr=xdf(n)**2+ydf(n)**2+zdf(n)**2
                          
c     put in primary list if found any interaction close enough
                          
                          if(rrr.le.rclim)then
                            ldo=.false.
                            list(ia,jj)=jneu
                          endif
                          
                        endif
                      endif
                    endif
                    
                  endif
                  
                enddo
                
              endif
              
            enddo
            
          enddo
          
        enddo
        
        if(mxnode.gt.1) call gstate(lchk)
        if(.not.lchk)then
          
          call gimax(ibig,1,idum)
          if(idnode.eq.0)then
            write(nrite,*) 'mxxdf must be at least ',ibig
            write(nrite,*) 'mxxdf is currently     ',mxxdf
          endif
          call  error(idnode,477)
          
        endif
        
      endif
      
      return
      end subroutine prneulst
      
      subroutine vertest(newlst,idnode,mxnode,natms,nbeads,delr,tstep)
      
c*********************************************************************
c     
c     DL_POLY subroutime to test for updating of Verlet list
c     replicated data version
c     
c     copyright daresbury laboratory 1993
c     author -       t. forester may 1993
c     
c*********************************************************************
      
      implicit none
      
      logical newlst,newjob
      integer idnode,mxnode,natms,nbeads,numatm,i,j,k,moved,ibuff,fail
      real(8) rmax,dr,delr,tstep
      
      real(8), allocatable :: xold(:),yold(:),zold(:)
      
      save newjob,xold,yold,zold
      
      data newjob/.true./,fail/0/
      
      numatm=nbeads*natms
      
      if((numatm-1)/mxnode+1.gt.mslist) call error(idnode,112)
      
      if(newjob)then
        
c     set up initial arrays 
        
        allocate (xold(mslist),yold(mslist),zold(mslist),stat=fail)
        if(fail.ne.0)call error(idnode,1930)
        
        j=0
        do i=idnode+1,numatm,mxnode
          
          j=j+1
          xold(j)=0.d0
          yold(j)=0.d0
          zold(j)=0.d0
          
        enddo
        
        newjob=.false.
        newlst=.true.
        
      else
        
c     integrate velocities 
        
        j=0
        do i=idnode+1,numatm,mxnode
          
          j=j+1
          xold(j)=xold(j)+vxx(i)
          yold(j)=yold(j)+vyy(i)
          zold(j)=zold(j)+vzz(i)
          
        enddo
        
c     maximum displacement 
        
        rmax=(delr/2.d0)**2
        
c     test atomic displacements
        
        moved=0
        
        do k=1,j
          
          dr=tstep**2*(xold(k)**2+yold(k)**2+zold(k)**2)
          if(dr.gt.rmax) moved=moved+1
          
        enddo
        
c     global sum of moved atoms
        
        if(mxnode.gt.1) call gisum(moved,1,ibuff)
        
c     test for new verlet list
        
        newlst=(moved.ge.2)
        
c     update stored positions
        
        if(newlst)then
          
          do k=1,j
            
            xold(k)=0.d0
            yold(k)=0.d0
            zold(k)=0.d0
            
          enddo
          
        endif
        
      endif
      
      return
      end subroutine vertest
      
      subroutine vertest2(newlst,idnode,mxnode,natms,imcon,delr,tstep)
      
c*********************************************************************
c     
c     DL_POLY subroutime to test for updating of Verlet list
c     replicated data version (version 2)
c     
c     copyright daresbury laboratory
c     author -       w.smith    2007
c     
c*********************************************************************
      
      implicit none
      
      logical newlst,newjob
      integer idnode,mxnode,natms,imcon,i,j,k,moved,ibuff,fail
      real(8) rmax,dr,delr,tstep
      
      real(8), allocatable :: xold(:),yold(:),zold(:)
      real(8), allocatable :: xdif(:),ydif(:),zdif(:)
      
      save newjob,xold,yold,zold
      
      data newjob/.true./,fail/0/
      
      if((natms+mxnode-1)/mxnode.gt.msatms) call error(idnode,112)
      
c     set up initial arrays 
      
      allocate (xdif(msatms),ydif(msatms),zdif(msatms),stat=fail)
      if(fail.ne.0)call error(idnode,1930)
      
      if(newjob)then
        
        allocate (xold(msatms),yold(msatms),zold(msatms),stat=fail)
        if(fail.ne.0)call error(idnode,1930)
        
        j=0
        do i=idnode+1,natms,mxnode
          
          j=j+1
          xold(j)=xxx(i)
          yold(j)=yyy(i)
          zold(j)=zzz(i)
          
        enddo
        
        newjob=.false.
        newlst=.true.
        
      else
        
c     calculate atomic shifts
        
        j=0
        do i=idnode+1,natms,mxnode
          
          j=j+1
          xdif(j)=xxx(i)-xold(j)
          ydif(j)=yyy(i)-yold(j)
          zdif(j)=zzz(i)-zold(j)
          
        enddo
        
c     minimum image calculation
        
        call images(imcon,0,1,j,cell,xdif,ydif,zdif)
        
c     maximum displacement 
        
        rmax=(delr/2.d0)**2
        
c     test atomic displacements
        
        moved=0
        
        do k=1,j
          
          dr=(xdif(k)**2+ydif(k)**2+zdif(k)**2)
          if(dr.gt.rmax)moved=moved+1
          
        enddo
        
c     global sum of moved atoms
        
        if(mxnode.gt.1) call gisum(moved,1,ibuff)
        
c     test for new verlet list
        
        newlst=(moved.ge.2)
        
c     update stored positions
        
        if(newlst)then
          
          j=0
          do i=idnode+1,natms,mxnode
            
            j=j+1
            xold(j)=xxx(i)
            yold(j)=yyy(i)
            zold(j)=zzz(i)
            
          enddo
          
        endif
        
      endif
      
c     deallocate arrays
      
      deallocate(xdif,ydif,zdif,stat=fail)
      
      return
      end subroutine vertest2
      
      end module nlist_builders_module
