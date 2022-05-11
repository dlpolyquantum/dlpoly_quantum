      module tether_module

c***********************************************************************
c     
c     dl_poly module for defining tether potential arrays
c     
c     copyright - daresbury laboratory
c     author    - w. smith    oct 2003
c     
c***********************************************************************

      use config_module
      use error_module
      use parse_module
      use setup_module
      use site_module
      use utility_module

      implicit none

      real(8), allocatable :: prmtet(:,:)
      integer, allocatable :: listtet(:,:)
      integer, allocatable :: numteth(:),keytet(:),lsttet(:)
      real(8), allocatable :: xxs(:),yys(:),zzs(:)

      save prmtet,lsttet,listtet,numteth,keytet,xxs,yys,zzs

      contains
      
      subroutine alloc_tet_arrays(idnode,mxnode)
      
      implicit none
      
      integer, parameter :: nnn=8

      logical safe
      integer i,fail,idnode,mxnode
      dimension fail(nnn)

      safe=.true.

c     allocate arrays

      fail(:)=0
      
      allocate (prmtet(mxteth,mxpbnd),stat=fail(1))
      allocate (numteth(mxtmls),stat=fail(2))
      allocate (keytet(mxteth),stat=fail(3))
      allocate (lsttet(mxteth),stat=fail(4))
      allocate (listtet(msteth,2),stat=fail(5))
      allocate (xxs(mxatms),stat=fail(6))
      allocate (yys(mxatms),stat=fail(7))
      allocate (zzs(mxatms),stat=fail(8))
      
      if(any(fail.gt.0))safe=.false.      
      if(mxnode.gt.1)call gstate(safe)    
      if(.not.safe)call error(idnode,1260)

c     initialise numteth array
      
      do i=1,mxtmls
         numteth(i)=0
      enddo
      
      end subroutine alloc_tet_arrays

      subroutine define_tethers
     x  (safe,idnode,itmols,nteth,nsite,engunit)

c***********************************************************************
c     
c     dl_poly subroutine for defining tether bonds
c     
c     copyright - daresbury laboratory
c     author    - w. smith    aug 2003
c     
c***********************************************************************

      implicit none

      logical safe
      character*8 keyword
      character*1 message(80)
      real(8) engunit, parpot(mxpbnd)
      integer idnode,itmols,nteth,nsite,ntmp,iteth,idum
      integer iatm1,isite1,j,keytmp

      ntmp=intstr(record,lenrec,idum)
      numteth(itmols)=numteth(itmols)+ntmp

      if(idnode.eq.0)then
        write(nrite,"(/,1x,'number of tethered atoms ',
     x    6x,i10)")numteth(itmols)
        write(nrite,"(/,' tethered atom details:',/,/,
     x    12x,'unit',5x,'key',6x,'atom',19x,'parameters',/) ")
      endif

      do iteth=1,ntmp

        nteth=nteth+1
        if(nteth.gt.mxteth)call error(idnode,62)

        call getrec(safe,idnode,nfield)
        if(.not.safe)return
        call strip(record,lenrec)
        call copystring(record,message,80)
        call lowcase(record,4)
        call getword(keyword,record,4,lenrec)

        if(keyword(1:4).eq.'harm')then
          keytmp=1
        elseif(keyword(1:4).eq.'rhrm')then
          keytmp=2
        elseif(keyword(1:4).eq.'quar')then
          keytmp=3
        else
          if(idnode.eq.0)write(nrite,*)message
          call error(idnode,450)
        endif

        iatm1=intstr(record,lenrec,idum)
        parpot(1)=dblstr(record,lenrec,idum)
        parpot(2)=dblstr(record,lenrec,idum)
        parpot(3)=dblstr(record,lenrec,idum)
        
        isite1=nsite-numsit(itmols)+iatm1
        
c     test for frozen atom 
        
        if(lfzsit(isite1).ne.0)then
          
          numteth(itmols)=numteth(itmols)-1
          if(idnode.eq.0)
     x      write(nrite,"(4x,a8,i10,4x,a4,i10,2x,10f15.6)")
     x      '*frozen*',iteth,keyword(1:4),iatm1,(parpot(j),j=1,mxpbnd)
          
        else
          
          if(idnode.eq.0)
     x      write(nrite,"(12x,i10,4x,a4,i10,2x,10f15.6)")
     x      iteth,keyword(1:4),iatm1,(parpot(j),j=1,mxpbnd)
          
        endif           

c     store parameters
        
        keytet(nteth)=keytmp
        lsttet(nteth)=iatm1
        prmtet(nteth,1)=parpot(1)
        prmtet(nteth,2)=parpot(2)
        prmtet(nteth,3)=parpot(3)

c     convert energy units to internal units
        
        if(abs(keytmp).eq.1)then
          prmtet(nteth,1)=prmtet(nteth,1)*engunit
        elseif(abs(keytmp).eq.2)then
          prmtet(nteth,1)=prmtet(nteth,1)*engunit
        elseif(abs(keytmp).eq.3)then
          prmtet(nteth,1)=prmtet(nteth,1)*engunit
          prmtet(nteth,2)=prmtet(nteth,2)*engunit
          prmtet(nteth,3)=prmtet(nteth,3)*engunit
        endif

      enddo

      return
      end subroutine define_tethers

      subroutine tethfrc
     x  (idnode,mxnode,imcon,natms,nstep,ntteth,engtet,virtet)

c***********************************************************************
c     
c     dl_poly routine to tether atoms to initial positions
c     includes stress tensor
c     
c     replicated data version : block data
c     
c     copyright daresbury laboratory 1994
c     author     t.forester feb 1994
c     amended    t.forester dec 1994 : block data
c     
c***********************************************************************

      implicit none

      logical safe
      integer idnode,mxnode,imcon,natms,nstep,ntteth,i,ii,ia,kk
      integer itet1,itet2,fail

      real(8) engtet,virtet,rab
      real(8) rrab,omega,gamma

      real(8), allocatable :: xdab(:),ydab(:),zdab(:)

      data safe/.true./
      
      allocate (xdab(msbad),ydab(msbad),zdab(msbad),stat=fail)
      if(fail.ne.0)call error(idnode,1270)

c     set up reference positions at start of job

      if(nstep.le.1)then

        do i=1,natms

          xxs(i)=xxx(i)
          yys(i)=yyy(i)
          zzs(i)=zzz(i)

        enddo

      endif

c     check size of work arrays

      if((ntteth-mxnode+1)/mxnode.gt.msbad) call error(idnode,420)

c     block indices

      itet1=(idnode*ntteth)/mxnode+1
      itet2=((idnode+1)*ntteth)/mxnode
      
      ii=0
      do i=itet1,itet2

        ii=ii+1

c     atomic index

        ia= listtet(ii,2)

c     tether vector

        xdab(ii)=xxx(ia)-xxs(ia)
        ydab(ii)=yyy(ia)-yys(ia)
        zdab(ii)=zzz(ia)-zzs(ia)

      enddo

c     ignore  periodic boundary condition
      
      call images(imcon,0,1,ii,cell,xdab,ydab,zdab)

c     zero tether energy and virial accumulators
      
      engtet=0.d0
      virtet=0.d0

c     loop over all specified tethered atoms

      ii=0
      do i=itet1,itet2
        
        ii=ii+1

c     define components of bond vector
        
        rab=sqrt(xdab(ii)**2+ydab(ii)**2+zdab(ii)**2)

c     check for possible zero length vector

        if(rab.lt.1.d-10)then

          rrab =0.d0

        else

          rrab=1.d0/rab

        endif

c     index of potential function parameters

        kk=listtet(ii,1)

c     calculate scalar constant terms

        if(keytet(kk).eq.1)then

c     harmonic function

          omega=0.5d0*prmtet(kk,1)*rab**2
          gamma=prmtet(kk,1)

        elseif(keytet(kk).eq.2)then

c     restrained harmonic: 

          omega=0.5d0*prmtet(kk,1)*(min(rab,prmtet(kk,2)))**2
     x      +prmtet(kk,1)*prmtet(kk,2)*
     x      (sign(max(rab-prmtet(kk,2),0.d0),rab))
          gamma=prmtet(kk,1)*(sign(min(rab,prmtet(kk,2)),rab))*rrab

        elseif(keytet(kk).eq.3)then

c     quartic potential

          omega=0.5d0*prmtet(kk,1)*rab**2 +
     x      1.d0/3.d0*prmtet(kk,2)*rab**3+
     x      0.25d0*prmtet(kk,3)*rab**4
          gamma=(prmtet(kk,1)*rab +
     x      prmtet(kk,2)*rab**2 +
     x      prmtet(kk,3)*rab**3)*rrab

        else
          safe=.false.
          omega=0.d0
          gamma=0.d0
        endif
        
        gamma=-gamma

c     calculate tether energy and virial

        engtet=engtet+omega
        virtet=virtet-gamma*rab*rab
        
c     index of atom
        
        ia=listtet(ii,2)

c     calculate atomic forces
        
        fxx(ia)=fxx(ia)+gamma*xdab(ii)
        fyy(ia)=fyy(ia)+gamma*ydab(ii)
        fzz(ia)=fzz(ia)+gamma*zdab(ii)

c     stress tensor 

        stress(1)=stress(1)+xdab(ii)*gamma*xdab(ii)
        stress(2)=stress(2)+xdab(ii)*gamma*ydab(ii)
        stress(3)=stress(3)+xdab(ii)*gamma*zdab(ii)
        stress(4)=stress(4)+ydab(ii)*gamma*xdab(ii)
        stress(5)=stress(5)+ydab(ii)*gamma*ydab(ii)
        stress(6)=stress(6)+ydab(ii)*gamma*zdab(ii)
        stress(7)=stress(7)+zdab(ii)*gamma*xdab(ii)
        stress(8)=stress(8)+zdab(ii)*gamma*ydab(ii)
        stress(9)=stress(9)+zdab(ii)*gamma*zdab(ii)

      enddo

c     check for undefined potentials

      if(mxnode.gt.1)call gstate(safe)
      if(.not.safe)call error(idnode,450)

c     sum contributions to potential and virial

      if(mxnode.gt.1)then

        buffer(3)=engtet
        buffer(4)=virtet

        call gdsum(buffer(3),2,buffer(1))

        engtet=buffer(3)
        virtet=buffer(4)

      endif

      deallocate (xdab,ydab,zdab,stat=fail)
      
      return
      end subroutine tethfrc
      
      subroutine xscale(idnode,mxnode,natms,keyens,imcon,tstep)

c***********************************************************************
c     
c     dl_poly routine to scale positions with change in box shape
c     
c     parallel replicated data version
c     
c     copyright daresbury laboratory 1995
c     author t.forester      october 1995
c     
c***********************************************************************

      implicit none
      
      integer idnode,mxnode,natms,keyens,imcon,iatm0,iatm1,i
      real(8) tstep,xa,ya,za,totmas,xcmo,ycmo,zcmo
      
c     assign block of atoms to processor

      iatm0=(idnode*natms)/mxnode+1
      iatm1=((idnode+1)*natms)/mxnode

      if((keyens.eq.4).or.(keyens.eq.6))then

c     berendsen npt/nst

        do i=iatm0,iatm1

          xa=eta(1)*xxs(i)+eta(2)*yys(i)+eta(3)*zzs(i)
          ya=eta(4)*xxs(i)+eta(5)*yys(i)+eta(6)*zzs(i)
          za=eta(7)*xxs(i)+eta(8)*yys(i)+eta(9)*zzs(i)

          xxs(i)=xa
          yys(i)=ya
          zzs(i)=za

        enddo

      elseif(keyens.eq.5.or.keyens.eq.7)then

c     hoover npt/nst
        
        totmas=0.d0
        do i=1,natms
          if(rmass(i).gt.0.d0)totmas=totmas+weight(i)
        enddo
        
        xcmo=0.d0
        ycmo=0.d0
        zcmo=0.d0

        do i=1,natms

          if(rmass(i).gt.0.d0)then

            xcmo=xcmo+weight(i)*xxs(i)
            ycmo=ycmo+weight(i)*yys(i)
            zcmo=zcmo+weight(i)*zzs(i)

          endif

        enddo
        xcmo=xcmo/totmas
        ycmo=ycmo/totmas
        zcmo=zcmo/totmas

        do i=iatm0,iatm1

          xa=xxs(i)-xcmo
          ya=yys(i)-ycmo
          za=zzs(i)-zcmo

          xxs(i)=xxs(i)+tstep*(eta(1)*xa+eta(2)*ya+eta(3)*za)
          yys(i)=yys(i)+tstep*(eta(2)*xa+eta(5)*ya+eta(6)*za)
          zzs(i)=zzs(i)+tstep*(eta(3)*xa+eta(6)*ya+eta(9)*za)

        enddo

        call images(imcon,idnode,mxnode,natms,cell,xxs,yys,zzs)

      endif

      if(mxnode.gt.1)
     x  call merge(idnode,mxnode,natms,mxbuff,xxs,yys,zzs,buffer)
      
      return 
      end subroutine xscale
      
      end module tether_module

