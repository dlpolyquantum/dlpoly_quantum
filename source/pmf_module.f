      module pmf_module
      
c***********************************************************************
c     
c     dl_poly module for defining potential of mean force arrays
c     copyright - daresbury laboratory
c     author    - w. smith    oct 2003
c     
c***********************************************************************
      
      use config_module
      use error_module
      use ensemble_tools_module
      use lf_motion_module
      use lf_rotation1_module
      use parse_module
      use property_module
      use setup_module
      use shake_module
      use vv_motion_module
      use utility_module
      
      implicit none
      
      integer npmf
      real(8) prmpmf,pmfnrm
      real(8), allocatable :: pmfwght(:)
      integer, allocatable :: numpmf(:)
      integer, allocatable :: indpmf(:)
      integer, allocatable :: listpm(:)
      integer, allocatable :: lstpmt(:)
      integer, allocatable :: lstpmf(:,:)
      
      dimension npmf(2),pmfnrm(2)
      
      save npmf,prmpmf,pmfnrm,pmfwght,numpmf,indpmf,listpm
      save lstpmt,lstpmf
      
      contains
      
      subroutine alloc_pmf_arrays(idnode,mxnode)
      
c***********************************************************************
c     
c     dl_poly subroutine for allocating pmf arrays
c     copyright - daresbury laboratory
c     author    - w. smith    oct 2003
c     
c***********************************************************************
      
      implicit none
      
      integer, parameter :: nnn=6
      
      logical safe
      integer i,idnode,mxnode,fail
      dimension fail(nnn)

      safe=.true.

c     allocate arrays

      fail(:)=0
      
      allocate (pmfwght(mxspmf),stat=fail(1))
      allocate (indpmf(mxspmf),stat=fail(2))
      allocate (numpmf(mxtmls),stat=fail(3))
      allocate (listpm(mxatms),stat=fail(4))
      allocate (lstpmt(mxatms),stat=fail(5))
      allocate (lstpmf(mxspmf,mspmf),stat=fail(6))

      if(any(fail.gt.0))safe=.false.            
      if(mxnode.gt.1)call gstate(safe)    
      if(.not.safe)call error(idnode,1210)

c     initialise array numpmf
      
      do i=1,mxtmls
        numpmf(i)=0
      enddo
      
      end subroutine alloc_pmf_arrays
      
      subroutine define_pmf(safe,idnode,itmols,nspmf)
      
c***********************************************************************
c     
c     dl_poly subroutine for defining pmf units
c     
c     copyright - daresbury laboratory
c     author    - w. smith    aug 2003
c     
c***********************************************************************
      
      implicit none
      
      logical safe
      integer idnode,itmols,nspmf,ipmf,jpmf,iatm1,idum
      real(8) wght
      
      numpmf(itmols)=1
      prmpmf=dblstr(record,lenrec,idum)
      
      if(idnode.eq.0) then
        write(nrite,"(/,1x,' PMF      bondlength :',
     x    5x,f20.10)") prmpmf
        write(nrite,
     x    "(/,/,12x,'unit, site and weight details:'
     x    ,/,/,16x,'unit',6x,'index',5x,'weight')")
      endif
      
      do ipmf=1,2
        
        call getrec(safe,idnode,nfield)
        if(.not.safe)return
        call strip(record,lenrec)
        call lowcase(record,lenrec)
        if(.not.findstring('pmf unit',record,idum))
     x    call error(idnode,462)
        npmf(ipmf)=intstr(record,lenrec,idum)
        
        do jpmf=1,npmf(ipmf)
          
          call getrec(safe,idnode,nfield)
          if(.not.safe)return
          
          iatm1=intstr(record,lenrec,idum)
          wght=dblstr(record,lenrec,idum)
          if(wght.le.1.d-10) wght=1.d0
          
          nspmf=nspmf+1
          
          if(nspmf.gt.mxspmf) call error(idnode,460)
          
          indpmf(nspmf)=iatm1
          pmfwght(nspmf)=wght
          
          if(idnode.eq.0) then
            
            if(jpmf.eq.1) then
              write(nrite,"(16x,i5,i10,f12.6)")
     x          ipmf,indpmf(nspmf),pmfwght(nspmf)
            else
              write(nrite,"(21x,i10,f12.6)")
     x          indpmf(nspmf),pmfwght(nspmf)
            endif
            
          endif
          
        enddo
        
      enddo
      
      return
      end subroutine define_pmf
      
      subroutine pmf_vectors
     x  (img,nspmf,imcon,cell,xxx,yyy,zzz,xxt,yyt,zzt,xxa,yya,zza,
     x  dxp,dyp,dzp)
      
c***********************************************************************
c     
c     dl_poly subroutine for constructing vectors for PMF calculations
c     
c     copyright - daresbury laboratory
c     adapted by w.smith october 2005
c     original by t.forester aug 1995
c     
c     set variable img true for PBC shifted vectors
c     
c***********************************************************************
      
      implicit none
      
      logical img
      integer nspmf,imcon,k,jj,kk,ipmf,i,i1,i2
      
      real(8) xxx(mxatms),yyy(mxatms),zzz(mxatms)
      real(8) xxt(mxatms),yyt(mxatms),zzt(mxatms)
      real(8) xxa(2,mspmf),yya(2,mspmf),zza(2,mspmf)
      real(8) dxp(mspmf),dyp(mspmf),dzp(mspmf),cell(9)
      
      do k=1,nspmf
        
        jj=0
        kk=0
        
c     calculate difference vectors
        
        do ipmf=1,2
          
          i1=lstpmf(jj+1,k)
          
c     position difference vectors
          
          do i=1,npmf(ipmf)
            
            jj=jj+1
            i2=lstpmf(jj,k)
            xxt(i)=xxx(i2)-xxx(i1)
            yyt(i)=yyy(i2)-yyy(i1)
            zzt(i)=zzz(i2)-zzz(i1)
            
          enddo
          
c     correct for periodic images - assume less than half box length
          
          if(img)call images(imcon,0,1,npmf(ipmf),cell,xxt,yyt,zzt)
          
c     create weighted coordinate
          
          xxa(ipmf,k)=0.d0
          yya(ipmf,k)=0.d0
          zza(ipmf,k)=0.d0
          
          do i=1,npmf(ipmf)
            
            kk=kk+1
            xxa(ipmf,k)=xxa(ipmf,k)+pmfwght(kk)*xxt(i)
            yya(ipmf,k)=yya(ipmf,k)+pmfwght(kk)*yyt(i)
            zza(ipmf,k)=zza(ipmf,k)+pmfwght(kk)*zzt(i)
            
          enddo
          
          xxa(ipmf,k)=xxa(ipmf,k)/pmfnrm(ipmf)+xxx(i1)
          yya(ipmf,k)=yya(ipmf,k)/pmfnrm(ipmf)+yyy(i1)
          zza(ipmf,k)=zza(ipmf,k)/pmfnrm(ipmf)+zzz(i1)
          
        enddo
        
        dxp(k)=xxa(2,k)-xxa(1,k)
        dyp(k)=yya(2,k)-yya(1,k)
        dzp(k)=zza(2,k)-zza(1,k)
        
      enddo
      
c     periodic boundary condition for pmf vectors
      
      if(img)call images(imcon,0,1,nspmf,cell,dxp,dyp,dzp)
      
      return
      end subroutine pmf_vectors
      
      subroutine pmflf
     x  (safe,safep,lshmov,idnode,imcon,mxnode,natms,nscons,
     x  ntcons,nspmf,ntpmf,engke,tolnce,tstep,vircon,virpmf)
      
c***********************************************************************
c     
c     dl_poly subroutine for integrating newtonian equations of
c     motion in molecular dynamics. Verlet leapfrog With RD-SHAKE
c     and PMF_SHAKE - for potential of mean force calculations.
c     
c     parallel replicated data version : block data
c     adapted from dl_poly routine nve_1.f
c     
c     copyright - daresbury laboratory 1995
c     author  - t.forester aug 1995
c     
c***********************************************************************
      
      implicit none
      
      integer, parameter :: nnn=8
      
      logical safe,lshmov,safep,unsafe,img
      integer idnode,imcon,mxnode,natms,nscons,ntcons,nspmf,ntpmf
      integer fail,iatm0,iatm1,i,j,k,jj,ii,ipmf,icyc
      real(8) engke,tolnce,tstep,vircon,virpmf,strpmf,summas
      real(8) rstep,viracc,strkin,strcon
      
      real(8), allocatable :: dxx(:),dyy(:),dzz(:)
      real(8), allocatable :: xxt(:),yyt(:),zzt(:)
      real(8), allocatable :: xxo(:),yyo(:),zzo(:)
      real(8), allocatable :: uxx(:),uyy(:),uzz(:)
      real(8), allocatable :: xxa(:,:),yya(:,:),zza(:,:)
      real(8), allocatable :: dxp(:),dyp(:),dzp(:)
      real(8), allocatable :: dxt(:),dyt(:),dzt(:)
      real(8), allocatable :: txx(:),tyy(:),tzz(:)
      dimension strpmf(9),strcon(9),strkin(9),summas(2),fail(nnn)
      
      do i=1,nnn
        fail(i)=0
      enddo
      allocate(dxx(mxcons),dyy(mxcons),dzz(mxcons),stat=fail(1))
      allocate(xxt(mxatms),yyt(mxatms),zzt(mxatms),stat=fail(2))
      allocate(xxo(msatms),yyo(msatms),zzo(msatms),stat=fail(3))
      allocate(uxx(mxatms),uyy(mxatms),uzz(mxatms),stat=fail(4))
      allocate(xxa(2,mspmf),yya(2,mspmf),zza(2,mspmf),stat=fail(5))
      allocate(dxp(mspmf),dyp(mspmf),dzp(mspmf),stat=fail(6))
      allocate(dxt(mxcons),dyt(mxcons),dzt(mxcons),stat=fail(7))
      allocate(txx(mxatms),tyy(mxatms),tzz(mxatms),stat=fail(8))
      do i=1,nnn
        if(fail(i).ne.0)call error(idnode,1220)
      enddo
      
c     block indices
      
      iatm0=(idnode*natms)/mxnode+1
      iatm1=((idnode+1)*natms)/mxnode
      
c     constraint virials
      
      vircon=0.d0
      virpmf=0.d0
      
c     temporary stress tensor accumulators
      
      do i=1,9
        
        strcns(i)=0.d0
        strpmf(i)=0.d0
        
      enddo
      
c     store initial values of position
      
      j=0
      do i=iatm0,iatm1
        
        j=j+1
        xxo(j)=xxx(i)
        yyo(j)=yyy(i)
        zzo(j)=zzz(i)
        
      enddo
      
c     construct current bond vectors
      
      do k=1,nscons
        
c     indices of atoms in bond
        
        i=listcon(k,2)
        j=listcon(k,3)
        
c     calculate current bond vector
        
        dxx(k)=xxx(i)-xxx(j)
        dyy(k)=yyy(i)-yyy(j)
        dzz(k)=zzz(i)-zzz(j)
        
      enddo
      
c     periodic boundary condition for bond vectors
      
      call images(imcon,0,1,nscons,cell,dxx,dyy,dzz)
      
c     calculate mass terms for PMF units
      
      jj=0
      do ipmf=1,2
        
        summas(ipmf)=0.d0
        pmfnrm(ipmf)=0.d0

        if(nspmf.gt.0)then
          
          do i=1,npmf(ipmf)
            
            jj=jj+1
            ii=lstpmf(jj,1)
            summas(ipmf)=summas(ipmf)+weight(ii)
            pmfnrm(ipmf)=pmfnrm(ipmf)+pmfwght(jj)
            
          enddo

        endif
        
      enddo
      
c     calculate PMF bond constraints and store initial positions
      
      img=.true.
      call pmf_vectors
     x  (img,nspmf,imcon,cell,xxx,yyy,zzz,xxt,yyt,zzt,
     x  xxa,yya,zza,dxp,dyp,dzp)
      
c     move atoms by leapfrog algorithm
      
      safe=(ntcons.eq.0)
      safep=(ntpmf.eq.0)
      
      j=0
      do i=iatm0,iatm1
        
        j=j+1
        
c     update velocities
        
        uxx(i)=vxx(i)+tstep*rmass(i)*fxx(i)
        uyy(i)=vyy(i)+tstep*rmass(i)*fyy(i)
        uzz(i)=vzz(i)+tstep*rmass(i)*fzz(i)
        
c     update positions
        
        xxx(i)=xxo(j)+tstep*uxx(i)
        yyy(i)=yyo(j)+tstep*uyy(i)
        zzz(i)=zzo(j)+tstep*uzz(i)
        
      enddo
      
c     RDSHAKE procedure 
      
      if(ntcons.gt.0.or.ntpmf.gt.0) then
        
c     global exchange of configuration data
        
        if(mxnode.gt.1)call merge
     x    (idnode,mxnode,natms,mxbuff,xxx,yyy,zzz,buffer)
        
c     apply constraint corrections - iteratively
        
        icyc=0
        unsafe=.true.
        
        do while(unsafe.and.icyc.lt.mxshak)
          
          icyc=icyc+1
          
c     apply bond constraints
          
          viracc=0.d0
          if(ntcons.gt.0)call rdshake_1
     x      (safe,lshmov,idnode,imcon,mxnode,natms,nscons,
     x      tolnce,tstep,viracc,dxx,dyy,dzz,dxt,dyt,dzt,
     x      txx,tyy,tzz,xxt,yyt,zzt,strcon)
          
          vircon=vircon+viracc
          do i=1,9
            strcns(i)=strcns(i)+strcon(i)
          enddo
          
c     apply pmf constraints
          
          viracc=0.d0
          if(ntpmf.gt.0)call pmf_shake
     x      (safep,idnode,mxnode,imcon,natms,nspmf,tstep,tolnce,
     x      viracc,xxt,yyt,zzt,strcon,summas,dxp,dyp,dzp,
     x      xxa,yya,zza)
          
          virpmf=virpmf+viracc
          do i=1,9
            strpmf(i)=strpmf(i)+strcon(i)
          enddo
          
          unsafe=(.not.(safe.and.safep.and.abs(viracc).lt.1.d-10))
          
        enddo
        
        safep=.not.unsafe
        
c     calculate velocity correction
        
        j=0
        rstep=1.d0/tstep
        do i=iatm0,iatm1
          
c     update corrected velocity
          
          j=j+1
          uxx(i)=(xxx(i)-xxo(j))*rstep
          uyy(i)=(yyy(i)-yyo(j))*rstep
          uzz(i)=(zzz(i)-zzo(j))*rstep
          
c     calculate the corrected forces
          
          fxx(i)=(uxx(i)-vxx(i))*weight(i)*rstep
          fyy(i)=(uyy(i)-vyy(i))*weight(i)*rstep
          fzz(i)=(uzz(i)-vzz(i))*weight(i)*rstep
          
        enddo
        
      endif
      
c     calculate velocity at full time step
      
      do i=iatm0,iatm1
        
        vxx(i)=0.5d0*(vxx(i)+uxx(i))
        vyy(i)=0.5d0*(vyy(i)+uyy(i))
        vzz(i)=0.5d0*(vzz(i)+uzz(i))
        
      enddo
      
c     calculate kinetic energy
      
      engke=getkin(natms,idnode,mxnode)
      
c     kinetic contribution to stress tensor
      
      call kinstress(natms,idnode,mxnode,strkin)
      
c     total contributions to stress tensor
      
      do i=1,9
        stress(i)=stress(i)+strpmf(i)+strcns(i)+strkin(i)
      enddo
      
c     add pmf and constraint virials
      
      vircon=vircon+virpmf
      
c     restore half step velocity
      
      do i=iatm0,iatm1
        
        vxx(i)=uxx(i)
        vyy(i)=uyy(i)
        vzz(i)=uzz(i)
        
      enddo
      
c     periodic boundary condition
      
      call images(imcon,idnode,mxnode,natms,cell,xxx,yyy,zzz)
      
c     global exchange of configuration data
      
      if(mxnode.gt.1)then
        
        call merge(idnode,mxnode,natms,mxbuff,xxx,yyy,zzz,buffer)
        call merge(idnode,mxnode,natms,mxbuff,vxx,vyy,vzz,buffer)
        if(ntcons.gt.0)call merge
     x    (idnode,mxnode,natms,mxbuff,fxx,fyy,fzz,buffer)
        
      endif
      
c     deallocate work arrays
      
      deallocate(dxx,dyy,dzz,xxt,yyt,zzt,stat=fail(1))
      deallocate(uxx,uyy,uzz,dxp,dyp,dzp,stat=fail(2))
      deallocate(txx,tyy,tzz,xxo,yyo,zzo,stat=fail(3))
      deallocate(dxt,dyt,dzt,xxa,yya,zza,stat=fail(4))
      
      return
      end subroutine pmflf
      
      subroutine pmflfq_1
     x  (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x  ntcons,ntfree,engke,engrot,quattol,tolnce,tstep,vircom,
     x  vircon,safep,nspmf,ntpmf,virpmf)
      
c***********************************************************************
c     
c     dlpoly routine to integrate rigid body equations of motion
c     using implicit leapfrog quaternion algorithm
c     plus bond constraints - PROVIDED rigid body sites
c     and constraint sites do not coincide.
c     
c     parallel replicated data version : block data
c     
c     tqx,tqy,tqz = torque in lab fixed frame
c     omx,omy,omz = angular velocity in body fixed frame (principal axis)
c     rotinx,y,z  = rotational inertia in body fixed frame
c     
c     copyright daresbury laboratory 1993
c     author      t.forester october 1993
c     amended     t.forester dec 1994 : block data
c     amended     w.smith sep 1999 : euler equation
c     
c**********************************************************************
      
      implicit none
      
      integer, parameter :: nnn=14
      real(8), parameter :: pt5=0.5d0
      
      logical safe,safeq,lshmov,newjob,safep,unsafe,img
      integer imcon,idnode,mxnode,natms,ngrp,nscons,ntcons,ntfree
      integer fail,i,igrp,igrp1,igrp2,idum,ifre1,ifre2,j,ifre
      integer jg,ig,k,id,jr,nspmf,ntpmf,jj,ii,ipmf,icyc
      real(8) engke,engrot,quattol,tolnce,tstep,vircom,vircon
      real(8) strkin,rot,rstep,rtsq,engtrn,vaa,vbb,vcc,virpmf
      real(8) trx,try,trz,delx,dely,delz,engfke,viracc
      real(8) strgrp,tqx,tqy,tqz,fmx,fmy,fmz
      real(8) strpmf,strcon,summas
      
      real(8), allocatable :: opx(:),opy(:),opz(:)
      real(8), allocatable :: oqx(:),oqy(:),oqz(:)
      real(8), allocatable :: dtx(:),dty(:),dtz(:)
      real(8), allocatable :: dxx(:),dyy(:),dzz(:)
      real(8), allocatable :: xxt(:),yyt(:),zzt(:)
      real(8), allocatable :: uxx(:),uyy(:),uzz(:)
      real(8), allocatable :: txx(:),tyy(:),tzz(:)
      real(8), allocatable :: dxt(:),dyt(:),dzt(:)
      real(8), allocatable :: xxo(:),yyo(:),zzo(:)
      real(8), allocatable :: xx1(:),yy1(:),zz1(:)
      real(8), allocatable :: vxo(:),vyo(:),vzo(:)
      real(8), allocatable :: gcxo(:),gcyo(:),gczo(:)
      
      
      real(8), allocatable :: xxa(:,:),yya(:,:),zza(:,:)
      real(8), allocatable :: dxp(:),dyp(:),dzp(:)
      
      dimension strkin(9),strgrp(9),rot(9),fail(nnn)
      
      dimension strpmf(9),strcon(9),summas(2)
      
      save igrp1,igrp2,ifre1,ifre2
      
      data newjob/.true./
      
c     allocate working arrays
      
      do i=1,nnn
        fail(i)=0
      enddo
      allocate (opx(msgrp),opy(msgrp),opz(msgrp),stat=fail(1))
      allocate (oqx(msgrp),oqy(msgrp),oqz(msgrp),stat=fail(2))
      allocate (dtx(mxatms),dty(mxatms),dtz(mxatms),stat=fail(3))
      allocate (dxx(mxcons),dyy(mxcons),dzz(mxcons),stat=fail(4))
      allocate (uxx(mxatms),uyy(mxatms),uzz(mxatms),stat=fail(5))
      allocate (xxt(mxatms),yyt(mxatms),zzt(mxatms),stat=fail(6))
      allocate (txx(mxatms),tyy(mxatms),tzz(mxatms),stat=fail(7))
      allocate (dxt(mxcons),dyt(mxcons),dzt(mxcons),stat=fail(8))
      allocate (xxo(msatms),yyo(msatms),zzo(msatms),stat=fail(9))
      allocate (xx1(msatms),yy1(msatms),zz1(msatms),stat=fail(10))
      allocate (vxo(msatms),vyo(msatms),vzo(msatms),stat=fail(11))
      allocate (gcxo(msgrp),gcyo(msgrp),gczo(msgrp),stat=fail(12))
      
      allocate(xxa(2,mspmf),yya(2,mspmf),zza(2,mspmf),stat=fail(13))
      allocate(dxp(mspmf),dyp(mspmf),dzp(mspmf),stat=fail(14))
      
      do i=1,nnn
        if(fail(i).ne.0)call error(idnode,1500)
      enddo
      
      if(newjob)then
        
c     group block indices
        
        igrp1=(idnode*ngrp)/mxnode+1
        igrp2=((idnode+1)*ngrp)/mxnode
        
c     free atom block indices
        
        ifre1=(idnode*ntfree)/mxnode+1
        ifre2=((idnode+1)*ntfree)/mxnode
        
c     check work arrays are large enough
        
        safe=(igrp2-igrp1+1.le.msgrp) 
        if(mxnode.gt.1) call gstate(safe)
        if(.not.safe) then 
          igrp=igrp2-igrp1+1
          call gimax(igrp,1,idum)
          if(idnode.eq.0) write(nrite,*) ' make msgrp >=',igrp
          call  error(idnode,506)
        endif
        
        newjob=.false.
        
      endif
      
      safe=.false.
      
c     constraint virials
      
      vircon=0.d0
      virpmf=0.d0
      
c     temporary stress tensor accumulators
      
      do i=1,9
        
        strcns(i)=0.d0
        strpmf(i)=0.d0
        
      enddo
      
c     store initial values of position and velocity
      
      j=0
      do ifre=ifre1,ifre2
        
        j=j+1
        i=lstfre(ifre)
        xxo(j)=xxx(i)
        yyo(j)=yyy(i)
        zzo(j)=zzz(i)
        vxo(j)=vxx(i)
        vyo(j)=vyy(i)
        vzo(j)=vzz(i)
        
      enddo
      
      jg=0
      do ig=igrp1,igrp2
        
        jg=jg+1
        gcxo(jg)=gcmx(ig)
        gcyo(jg)=gcmy(ig)
        gczo(jg)=gcmz(ig)
        
      enddo
      
c     construct current bond vectors
      
      do k=1,nscons
        
        i=listcon(k,2)
        j=listcon(k,3)
        
        dxx(k)=xxx(i)-xxx(j)
        dyy(k)=yyy(i)-yyy(j)
        dzz(k)=zzz(i)-zzz(j)
        
      enddo
      
c     periodic boundary condition for bond vectors
      
      call images(imcon,0,1,nscons,cell,dxx,dyy,dzz)
      
c     calculate mass terms for PMF units
      
      jj=0
      do ipmf=1,2
        
        summas(ipmf)=0.d0
        pmfnrm(ipmf)=0.d0
        
        do i=1,npmf(ipmf)
          
          jj=jj+1
          ii=lstpmf(jj,1)
          summas(ipmf)=summas(ipmf)+weight(ii)
          pmfnrm(ipmf)=pmfnrm(ipmf)+pmfwght(jj)
          
        enddo
        
      enddo
      
c     calculate PMF bond constraints and store initial positions
      
      img=.true.
      call pmf_vectors
     x  (img,nspmf,imcon,cell,xxx,yyy,zzz,xxt,yyt,zzt,
     x  xxa,yya,zza,dxp,dyp,dzp)
      
c     move atoms by leapfrog algorithm
      
      safe=(ntcons.eq.0)
      safep=(ntpmf.eq.0)
      
c     calculate atom displacements from rigid body com
      
      jg=0
      jr=0
      do ig=igrp1,igrp2
        
        jg=jg+1
        id=lstgtp(ig)
        
        do j=1,numgsit(id)
          
          jr=jr+1
          i=lstrgd(jr)
          
          dtx(jr)=xxx(i)-gcxo(jg)
          dty(jr)=yyy(i)-gcyo(jg)
          dtz(jr)=zzz(i)-gczo(jg)
          
        enddo
        
      enddo
      
c     minimum images
      
      call images(imcon,0,1,jr,cell,dtx,dty,dtz)
      
c     integrate 'free' particles
      
      j=0
      do ifre=ifre1,ifre2
        
        j=j+1
        i=lstfre(ifre)
        
c     advance velocity by leapfrog
        
        uxx(i)=vxo(j)+tstep*rmass(i)*fxx(i)
        uyy(i)=vyo(j)+tstep*rmass(i)*fyy(i)
        uzz(i)=vzo(j)+tstep*rmass(i)*fzz(i)
        
c     advance position by leapfrog
        
        xxx(i)=xxo(j)+tstep*uxx(i)
        yyy(i)=yyo(j)+tstep*uyy(i)
        zzz(i)=zzo(j)+tstep*uzz(i)
        
      enddo
      
      if(ntcons.gt.0.or.ntpmf.gt.0) then
        
c     store integrated positions
        
        j=0
        do ifre=ifre1,ifre2
          
          i=lstfre(ifre)
          j=j+1
          xx1(j)=xxx(i)
          yy1(j)=yyy(i)
          zz1(j)=zzz(i)
          
        enddo
        
c     global exchange of configuration data
        
        if(mxnode.gt.1)call merge1
     x    (idnode,mxnode,natms,lstme,xxx,yyy,zzz,buffer)
        
c     apply constraint corrections - iteratively
        
        icyc=0
        unsafe=.true.
        
        do while(unsafe.and.icyc.lt.mxshak)
          
          icyc=icyc+1
          
c     apply bond constraints
          
          viracc=0.d0
          if(ntcons.gt.0)call rdshake_1
     x      (safe,lshmov,idnode,imcon,mxnode,natms,nscons,
     x      tolnce,tstep,viracc,dxx,dyy,dzz,dxt,dyt,dzt,
     x      txx,tyy,tzz,xxt,yyt,zzt,strcon)
          
          vircon=vircon+viracc
          do i=1,9
            strcns(i)=strcns(i)+strcon(i)
          enddo
          
c     apply pmf constraints
          
          viracc=0.d0
          if(ntpmf.gt.0)call pmf_shake
     x      (safep,idnode,mxnode,imcon,natms,nspmf,tstep,tolnce,
     x      viracc,xxt,yyt,zzt,strcon,summas,dxp,dyp,dzp,
     x      xxa,yya,zza)
          
          virpmf=virpmf+viracc
          do i=1,9
            strpmf(i)=strpmf(i)+strcon(i)
          enddo
          
          unsafe=(.not.(safe.and.safep.and.abs(viracc).lt.1.d-10))
          
        enddo
        
        safep=.not.unsafe
        
c     calculate force and velocity corrections
        
        j=0
        rstep=1.d0/tstep
        rtsq=1.d0/tstep**2
        do ifre=ifre1,ifre2
          
          j=j+1
          i=lstfre(ifre)
          
c     calculate force correction
          
          fxx(i)=fxx(i)+(xxx(i)-xx1(j))*weight(i)*rtsq
          fyy(i)=fyy(i)+(yyy(i)-yy1(j))*weight(i)*rtsq
          fzz(i)=fzz(i)+(zzz(i)-zz1(j))*weight(i)*rtsq
          
c     calculate velocity correction
          
          uxx(i)=uxx(i)+(xxx(i)-xx1(j))*rstep
          uyy(i)=uyy(i)+(yyy(i)-yy1(j))*rstep
          uzz(i)=uzz(i)+(zzz(i)-zz1(j))*rstep
          
        enddo
        
c     end of shake corrections
        
      endif
      
c     estimate full step velocity
      
      j=0
      do ifre=ifre1,ifre2
        
        j=j+1
        i=lstfre(ifre)
        
        vxx(i)=0.5d0*(uxx(i)+vxo(j))
        vyy(i)=0.5d0*(uyy(i)+vyo(j))
        vzz(i)=0.5d0*(uzz(i)+vzo(j))
        
      enddo
      
c     calculate new kinetic energy at current timestep
      
      engfke=getkinf(ntfree,idnode,mxnode)
      
c     kinetic contribution to stress tensor
      
      call kinstressf(ntfree,idnode,mxnode,strkin)
      
c     restore free atom half step velocity
      
      do ifre=ifre1,ifre2
        
        i=lstfre(ifre)
        vxx(i)=uxx(i)
        vyy(i)=uyy(i)
        vzz(i)=uzz(i)
        
      enddo
      
c     *************  Rigid body motion ****************************
      
c     translational rigid body motion
      
      jg=0
      jr=0
      do ig=igrp1,igrp2
        
        jg=jg+1
        id=lstgtp(ig)
        
c     calculate net force on rigid body
        
        fmx=0.d0
        fmy=0.d0
        fmz=0.d0
        do j=1,numgsit(id)
          
          jr=jr+1
          i=lstrgd(jr)
          fmx=fmx+fxx(i)
          fmy=fmy+fyy(i)
          fmz=fmz+fzz(i)
          
        enddo
        
c     advance velocity by leapfrog
        
        uxx(ig)=gvxx(ig)+fmx*tstep/gmass(id)
        uyy(ig)=gvyy(ig)+fmy*tstep/gmass(id)
        uzz(ig)=gvzz(ig)+fmz*tstep/gmass(id)
        
c     advance position by leapfrog
        
        gcmx(ig)=gcmx(ig)+tstep*uxx(ig)
        gcmy(ig)=gcmy(ig)+tstep*uyy(ig)
        gcmz(ig)=gcmz(ig)+tstep*uzz(ig)
        
c     estimate velocity at full time step
        
        gvxx(ig)=0.5d0*(gvxx(ig)+uxx(ig))
        gvyy(ig)=0.5d0*(gvyy(ig)+uyy(ig))
        gvzz(ig)=0.5d0*(gvzz(ig)+uzz(ig))
        
      enddo
      
c     calculate rigid body translational kinetic energy
      
      engtrn=getkint(ngrp,idnode,mxnode)
      
c     total translational kinetic energy
      
      engke=engtrn+engfke
      
c     calculate ridid body kinetic stress tensor
      
      call kinstressg(ngrp,idnode,mxnode,strgrp)
      
c     restore rigid body half timestep velocity
      
      do ig=igrp1,igrp2
        
        gvxx(ig)=uxx(ig)
        gvyy(ig)=uyy(ig)
        gvzz(ig)=uzz(ig)
        
      enddo
      
c     calculate rigid body contribution to stress tensor
      
      call bodystress(idnode,mxnode,ngrp,vircom,strbod,dtx,dty,dtz)
      
c     calculate torques in lab frame
      
      jr=0
      jg=0
      do ig=igrp1,igrp2
        
        jg=jg+1
        id=lstgtp(ig)
        
        tqx=0.d0
        tqy=0.d0
        tqz=0.d0
        do j=1,numgsit(id)
          
          jr=jr+1
          i=lstrgd(jr)
          
          tqx=tqx+dty(jr)*fzz(i)-dtz(jr)*fyy(i)
          tqy=tqy+dtz(jr)*fxx(i)-dtx(jr)*fzz(i)
          tqz=tqz+dtx(jr)*fyy(i)-dty(jr)*fxx(i)
          
        enddo
        
c     current rotational matrix 
        
        call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)
        
c     store current angular velocity
        
        opx(jg)=omx(ig)
        opy(jg)=omy(ig)
        opz(jg)=omz(ig)
        
c     iterate angular velocity for time step n (e. yezdimer)
        
        do i=1,5
          
          trx=(tqx*rot(1)+tqy*rot(4)+tqz*rot(7))*rotinx(id,2)
     x      +(rotiny(id,1)-rotinz(id,1))*opy(jg)*opz(jg)*rotinx(id,2)
          try=(tqx*rot(2)+tqy*rot(5)+tqz*rot(8))*rotiny(id,2)
     x      +(rotinz(id,1)-rotinx(id,1))*opz(jg)*opx(jg)*rotiny(id,2)
          trz=(tqx*rot(3)+tqy*rot(6)+tqz*rot(9))*rotinz(id,2)
     x      +(rotinx(id,1)-rotiny(id,1))*opx(jg)*opy(jg)*rotinz(id,2)
          
          delx=tstep*trx
          dely=tstep*try
          delz=tstep*trz
          
c     improved angular velocity at time step n
          
          opx(jg)=omx(ig)+delx*pt5
          opy(jg)=omy(ig)+dely*pt5
          opz(jg)=omz(ig)+delz*pt5
          
        enddo
        
c     angular velocity at time step n+1/2
        
        uxx(ig)=omx(ig)+delx
        uyy(ig)=omy(ig)+dely
        uzz(ig)=omz(ig)+delz
        
c     angular velocity at time step n+1  (needed for quat algorithm)
        
        oqx(jg)=omx(ig)+delx*1.5d0
        oqy(jg)=omy(ig)+dely*1.5d0
        oqz(jg)=omz(ig)+delz*1.5d0
        
c     angular velocity at timestep n
        
        omx(ig)=omx(ig)+pt5*delx
        omy(ig)=omy(ig)+pt5*dely
        omz(ig)=omz(ig)+pt5*delz
        
      enddo
      
c     rotational kinetic energy
      
      engrot=getkinr(ngrp,idnode,mxnode)
      
c     restore half step angular velocity
      
      jg=0
      do ig=igrp1,igrp2
        
        jg=jg+1
        opx(jg)=omx(ig)
        opy(jg)=omy(ig)
        opz(jg)=omz(ig)
        omx(ig)=uxx(ig)
        omy(ig)=uyy(ig)
        omz(ig)=uzz(ig)
        
      enddo
      
c     assign new quaternions
      
      call update_quaternions
     x  (safeq,igrp1,igrp2,tstep,quattol,opx,opy,opz,oqx,oqy,oqz)
      
c     complete stress tensor
      
      do i=1,9
        stress(i)=stress(i)+strkin(i)+strgrp(i)+strcns(i)+strbod(i)+
     x    strpmf(i)
      enddo
      
c     add pmf and constraint virials
      
      vircon=vircon+virpmf
      
c     minimum images of group positions and particle positions
      
      call images(imcon,0,1,natms,cell,xxx,yyy,zzz)
      call images(imcon,idnode,mxnode,ngrp,cell,gcmx,gcmy,gcmz)
      
c     new atomic positions for atoms in rigid bodies - relative to c.o.m
      
      jr=0
      do ig=igrp1,igrp2
        
        id=lstgtp(ig)
        
c     new rotational matrix
        
        call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)
        
        do j=1,numgsit(id)
          
          jr=jr+1
          i=lstrgd(jr)
          
          xxx(i)=rot(1)*gxx(id,j)+rot(2)*gyy(id,j)+rot(3)*gzz(id,j)
     x      +gcmx(ig)
          yyy(i)=rot(4)*gxx(id,j)+rot(5)*gyy(id,j)+rot(6)*gzz(id,j)
     x      +gcmy(ig)
          zzz(i)=rot(7)*gxx(id,j)+rot(8)*gyy(id,j)+rot(9)*gzz(id,j)
     x      +gcmz(ig)
          
c     new atomic velocites in body frame
          
          vaa=omy(ig)*gzz(id,j)-omz(ig)*gyy(id,j)
          vbb=omz(ig)*gxx(id,j)-omx(ig)*gzz(id,j)
          vcc=omx(ig)*gyy(id,j)-omy(ig)*gxx(id,j)
          
c     new atomic velocites in lab frame
          
          vxx(i)=rot(1)*vaa+rot(2)*vbb+rot(3)*vcc+gvxx(ig)
          vyy(i)=rot(4)*vaa+rot(5)*vbb+rot(6)*vcc+gvyy(ig)
          vzz(i)=rot(7)*vaa+rot(8)*vbb+rot(9)*vcc+gvzz(ig)
          
        enddo
        
      enddo
      
      if(mxnode.gt.1) then
        
c     merge new group coordinates and velocities
        
        call merge(idnode,mxnode,ngrp,mxbuff,gcmx,gcmy,gcmz,buffer)
        call merge(idnode,mxnode,ngrp,mxbuff,gvxx,gvyy,gvzz,buffer)
        
c     merge new atomic coordinates and velocities
        
        call merge1(idnode,mxnode,natms,lstme,xxx,yyy,zzz,buffer)
        call merge1(idnode,mxnode,natms,lstme,vxx,vyy,vzz,buffer)
        
      endif
      
c     ensure all atoms are within cell boundaries
      
      call images(imcon,0,1,natms,cell,xxx,yyy,zzz)
      
c     deallocate work arrays
      
      deallocate(opx,opy,opz,xxt,yyt,zzt,stat=fail(1))
      deallocate(oqx,oqy,oqz,dtx,dty,dtz,stat=fail(2))
      deallocate(dxx,dyy,dzz,uxx,uyy,uzz,stat=fail(3))
      deallocate(txx,tyy,tzz,dxt,dyt,dzt,stat=fail(4))
      deallocate(xxo,yyo,zzo,xx1,yy1,zz1,stat=fail(5))
      deallocate(vxo,vyo,vzo,gcxo,gcyo,gczo,stat=fail(6))
      deallocate(dxp,dyp,dzp,xxa,yya,zza,stat=fail(7))
      
      return
      end subroutine pmflfq_1
      
      subroutine pmf_shake
     x  (safep,idnode,mxnode,imcon,natms,nspmf,tstep,tolnce,
     x  virpmf,xxt,yyt,zzt,strpmf,summas,dxp,dyp,dzp,
     x  xxa,yya,zza)
      
c***********************************************************************
c     
c     dlpoly constraint subroutine for potential of mean force calc.
c     accummulates constraint force to maintain reaction coordinate
c     
c     assume bond vectors dxp,dyp,dzp are input
c     dxp=(sum) wght*xxx(i,1) - (sum) wght*xxx(j,2) etc
c     
c     copyright daresbury laboratory 1995
c     author t.forester august 1995
c     
c***********************************************************************
      
      implicit none
      
      logical safep,img
      integer idnode,mxnode,imcon,natms,nspmf,fail,icyc,k,jj
      integer ii,i,ipmf
      real(8) tstep,tolnce,virpmf,xxt,yyt,zzt,strpmf,summas
      real(8) dxp,dyp,dzp,xxa,yya,zza,amt,strs1,strs2,strs3,strs5
      real(8) strs6,strs9,tstep2,dis,omega2,eps,gamma,gammi
      
      dimension dxp(mspmf),dyp(mspmf),dzp(mspmf)
      dimension xxt(mxatms),yyt(mxatms),zzt(mxatms)
      dimension xxa(2,mspmf),yya(2,mspmf),zza(2,mspmf)
      dimension amt(2),strpmf(9),summas(2)
      
      real(8), allocatable :: dxt(:),dyt(:),dzt(:),dsq(:)
      data fail/0/
      
c     allocate work arrays
      
      allocate (dxt(mxcons),dyt(mxcons),dzt(mxcons),dsq(mspmf),
     x  stat=fail)
      if(fail.ne.0)call error(idnode,1230)
      
      if(mxcons.lt.nspmf) call error(idnode,492)
      if(mspmf .lt.nspmf) call error(idnode,458)
      
c     timestep squared
      
      tstep2=tstep*tstep
      
c     accumulators for strpmf tensor
      
      virpmf=0.d0
      strs1=0.d0
      strs2=0.d0
      strs3=0.d0
      strs5=0.d0
      strs6=0.d0
      strs9=0.d0
      
c     application of constraint (shake) algorithm
      
      icyc=0
      safep=.false.
      
      do while(.not.safep.and.icyc.lt.mxshak)
        
        icyc=icyc+1
        
c     calculate bond constraint length
        
        img=.true.
        call pmf_vectors
     x    (img,nspmf,imcon,cell,xxx,yyy,zzz,xxt,yyt,zzt,
     x    xxa,yya,zza,dxt,dyt,dzt)
        
        amt(1)= tstep2/summas(1)
        amt(2)=-tstep2/summas(2)
        
        eps=0.d0
        dis=prmpmf
        omega2=dis*dis
        
        do k=1,nspmf
          
          dsq(k)=dxt(k)**2+dyt(k)**2+dzt(k)**2
          eps=max(eps,abs((omega2-dsq(k))/(2.0d0*dis)))
          
        enddo
        
c     check convergence condition
        
        safep=(eps.le.tolnce)
        
c     bypass calculations if converged
        
        if(.not.safep)then
          
          do k=1,nspmf
            
            gamma=(omega2-dsq(k))/(-2.d0*(amt(2)-amt(1))*
     x        (dxp(k)*dxt(k)+dyp(k)*dyt(k)+dzp(k)*dzt(k)))
            
c     accumulate pmf virial
            
            virpmf=virpmf+gamma*(dxp(k)**2+dyp(k)**2+dzp(k)**2)
            
            strs1=strs1-gamma*dxp(k)*dxp(k)
            strs2=strs2-gamma*dxp(k)*dyp(k)
            strs3=strs3-gamma*dxp(k)*dzp(k)
            strs5=strs5-gamma*dyp(k)*dyp(k)
            strs6=strs6-gamma*dyp(k)*dzp(k)
            strs9=strs9-gamma*dzp(k)*dzp(k)
            
c     improve approximate atomic positions
            
            jj=0
            do ipmf=1,2
              
              gammi=-gamma*amt(ipmf)
              
              do ii=1,npmf(ipmf)
                
                jj=jj+1
                i=lstpmf(jj,k)
                
                xxx(i)=xxx(i)+dxp(k)*gammi
                yyy(i)=yyy(i)+dyp(k)*gammi
                zzz(i)=zzz(i)+dzp(k)*gammi
                
              enddo
              
            enddo
            
          enddo
          
        endif
        
      enddo
      
c     complete strpmf tensor
      
      strpmf(1)=strs1
      strpmf(2)=strs2
      strpmf(3)=strs3
      strpmf(4)=strs2
      strpmf(5)=strs5
      strpmf(6)=strs6
      strpmf(7)=strs3
      strpmf(8)=strs6
      strpmf(9)=strs9
      
c     splice coordinate arrays across nodes
      
      if(mxnode.gt.1)then
        
        buffer(1)=virpmf
        call gdsum(buffer(1),1,buffer(2))
        virpmf=buffer(1)
        call gdsum(strpmf,9,buffer)
        call splice 
     x    (idnode,natms,listpm,lstpmt,xxx,yyy,zzz,buffer)
        
      endif
      
c     deallocate work arrays
      
      deallocate(dxt,dyt,dzt,dsq,stat=fail)
      
      return
      end subroutine pmf_shake
      
      subroutine pmfvv
     x  (safe,safep,lshmov,isw,idnode,mxnode,imcon,natms,nscons,
     x  ntcons,nspmf,ntpmf,engke,tolnce,tstep,vircon,virpmf)
      
c***********************************************************************
c     
c     dl_poly subroutine for integrating newtonian equations of
c     motion in molecular dynamics. Velocity Verlet With PMF_RATTLE
c     for potential of mean force calculations.
c     
c     copyright - daresbury laboratory
c     adapted by w.smith october 2005
c     original by t.forester aug 1995
c     
c***********************************************************************
      
      implicit none
      
      integer, parameter :: nnn=7
      
      logical safe,lshmov,safep,unsafe,newjob,img
      integer idnode,imcon,mxnode,natms,nscons,ntcons,nspmf,ntpmf
      integer isw,mxtop
      integer fail,iatm0,iatm1,i,j,k,jj,i1,ipmf,icyc
      real(8) engke,tolnce,tstep,vircon,virpmf,strcon,summas
      real(8) viracc,strpmf,strkin
      
      real(8), allocatable :: dxx(:),dyy(:),dzz(:)
      real(8), allocatable :: xxt(:),yyt(:),zzt(:)
      real(8), allocatable :: dxp(:),dyp(:),dzp(:)
      real(8), allocatable :: dxt(:),dyt(:),dzt(:)
      real(8), allocatable :: txx(:),tyy(:),tzz(:)
      real(8), allocatable :: xxa(:,:),yya(:,:),zza(:,:)
      
      dimension strcon(9),strpmf(9),strkin(9),summas(2),fail(nnn)
      save summas,newjob,strpmf
      data newjob/.true./
      
c     set default safety flags
      
      safe=(ntcons.eq.0)
      safep=(ntpmf.eq.0)
      if(safe.and.safep)return
      
      do i=1,nnn
        fail(i)=0
      enddo
      mxtop=max(mxcons,mspmf)
      allocate(dxx(mxcons),dyy(mxcons),dzz(mxcons),stat=fail(1))
      allocate(xxt(mxatms),yyt(mxatms),zzt(mxatms),stat=fail(2))
      allocate(dxp(mspmf),dyp(mspmf),dzp(mspmf),stat=fail(3))
      allocate(dxt(mxtop),dyt(mxtop),dzt(mxtop),stat=fail(4))
      allocate(txx(mxatms),tyy(mxatms),tzz(mxatms),stat=fail(5))
      allocate(xxa(2,mspmf),yya(2,mspmf),zza(2,mspmf),stat=fail(6))
      do i=1,nnn
        if(fail(i).ne.0)call error(idnode,1220)
      enddo
      
c     block indices
      
      iatm0=(idnode*natms)/mxnode+1
      iatm1=((idnode+1)*natms)/mxnode
      
c     calculate mass terms for PMF units
      
      if(newjob)then
        
        jj=0
        do ipmf=1,2
          
          summas(ipmf)=0.d0
          pmfnrm(ipmf)=0.d0
          
          if(nspmf.gt.0)then
            
            do i=1,npmf(ipmf)
              
              jj=jj+1
              i1=lstpmf(jj,1)
              summas(ipmf)=summas(ipmf)+weight(i1)
              pmfnrm(ipmf)=pmfnrm(ipmf)+pmfwght(jj)
              
            enddo

          endif
          
        enddo
        
        newjob=.false.
        
      endif
      
c     construct current bond vectors
      
      do k=1,nscons
        
c     indices of atoms in bond
        
        i=listcon(k,2)
        j=listcon(k,3)
        
c     calculate current bond vector
        
        dxx(k)=xxx(i)-xxx(j)
        dyy(k)=yyy(i)-yyy(j)
        dzz(k)=zzz(i)-zzz(j)
        
      enddo
      
c     periodic boundary condition for bond vectors
      
      call images(imcon,0,1,nscons,cell,dxx,dyy,dzz)
      
c     calculate PMF bond constraints and store initial positions
      
      img=.true.
      call pmf_vectors
     x  (img,nspmf,imcon,cell,xxx,yyy,zzz,xxt,yyt,zzt,xxa,yya,zza,
     x  dxp,dyp,dzp)
      
c     update velocities
      
      do i=iatm0,iatm1
        
        vxx(i)=vxx(i)+0.5d0*tstep*rmass(i)*fxx(i)
        vyy(i)=vyy(i)+0.5d0*tstep*rmass(i)*fyy(i)
        vzz(i)=vzz(i)+0.5d0*tstep*rmass(i)*fzz(i)
        
      enddo
      
      if(mxnode.gt.1)
     x  call merge(idnode,mxnode,natms,mxbuff,vxx,vyy,vzz,buffer)
      
      if(isw.eq.1)then
        
c     constraint virials
        
        vircon=0.d0
        virpmf=0.d0
        
c     temporary stress tensor accumulators
        
        do i=1,9
          
          strcns(i)=0.d0
          strpmf(i)=0.d0
          
        enddo
        
c     update positions
        
        do i=iatm0,iatm1
          
          xxx(i)=xxx(i)+tstep*vxx(i)
          yyy(i)=yyy(i)+tstep*vyy(i)
          zzz(i)=zzz(i)+tstep*vzz(i)
          
        enddo
        
        if(mxnode.gt.1)call merge
     x    (idnode,mxnode,natms,mxbuff,xxx,yyy,zzz,buffer)
        
      endif
      
c     apply constraint corrections - iteratively
      
      icyc=0
      unsafe=.true.
      
      do while(unsafe.and.icyc.lt.mxshak)
        
        icyc=icyc+1
        
        if(isw.eq.1)then
          
c     apply bond constraints
          
          if(ntcons.gt.0)then
            
            safe=.false.
            call rdrattle_r
     x        (safe,lshmov,idnode,imcon,mxnode,natms,nscons,
     x        tolnce,tstep,viracc,dxx,dyy,dzz,dxt,dyt,dzt,
     x        txx,tyy,tzz,xxt,yyt,zzt,strcon)
            
            vircon=vircon+viracc
            do i=1,9
              strcns(i)=strcns(i)+strcon(i)
            enddo
            
          endif
          
c     apply pmf constraints
          
          if(ntpmf.gt.0)then
            
            safep=.false.
            call pmf_rattle_r
     x        (safep,idnode,mxnode,imcon,natms,nspmf,tstep,tolnce,
     x        viracc,summas,xxt,yyt,zzt,xxa,yya,zza,dxp,dyp,dzp,
     x        dxt,dyt,dzt,strcon)
            
            virpmf=virpmf+viracc
            do i=1,9
              strpmf(i)=strpmf(i)+strcon(i)
            enddo
            
          endif
          
          unsafe=(.not.(safe.and.safep.and.abs(viracc).le.1.d-10))
          
        endif
        
        if(isw.eq.2)then
          
c     apply rattle velocity constraints
          
          if(ntcons.gt.0)then
            
            safe=.false.
            call rdrattle_v
     x        (safe,idnode,mxnode,natms,nscons,tolnce,tstep,
     x        dxx,dyy,dzz,txx,tyy,tzz,xxt,yyt,zzt)
            
          endif
          
c     apply pmf velocity constraints
          
          if(ntpmf.gt.0)then
            
            safep=.false.
            call pmf_rattle_v
     x        (safep,idnode,mxnode,imcon,natms,nspmf,tstep,tolnce,
     x        summas,dxp,dyp,dzp,xxt,yyt,zzt,xxa,yya,zza,dxt,dyt,dzt)
            
          endif
          
          unsafe=(.not.(safe.and.safep))
          
        endif
        
      enddo
      
      safep=(.not.unsafe)
      
c     periodic boundary condition
      
      if(isw.eq.1)call images
     x  (imcon,idnode,mxnode,natms,cell,xxx,yyy,zzz)
      
      if(isw.eq.2)then
        
c     calculate kinetic energy
        
        engke=getkin(natms,idnode,mxnode)
        
c     kinetic contribution to stress tensor
        
        call kinstress(natms,idnode,mxnode,strkin)
        
c     total contributions to stress tensor
        
        do i=1,9
          stress(i)=stress(i)+strpmf(i)+strcns(i)+strkin(i)
        enddo
        
c     add pmf and constraint virials
        
        vircon=vircon+virpmf
        
      endif
      
c     global exchange of configuration data
      
      if(mxnode.gt.1)then
        
        call merge(idnode,mxnode,natms,mxbuff,vxx,vyy,vzz,buffer)
        if(isw.eq.1)call merge
     x    (idnode,mxnode,natms,mxbuff,xxx,yyy,zzz,buffer)
        if(ntcons.gt.0)call merge
     x    (idnode,mxnode,natms,mxbuff,fxx,fyy,fzz,buffer)
        
      endif
      
c     deallocate work arrays
      
      deallocate(dxx,dyy,dzz,xxt,yyt,zzt,stat=fail(1))
      deallocate(dxp,dyp,dzp,dxt,dyt,dzt,stat=fail(2))
      deallocate(txx,tyy,tzz,xxa,yya,zza,stat=fail(3))
      
      return
      end subroutine pmfvv
      
      subroutine pmf_rattle_r
     x  (safep,idnode,mxnode,imcon,natms,nspmf,tstep,tolnce,
     x  virpmf,summas,xxt,yyt,zzt,xxa,yya,zza,dxp,dyp,dzp,
     x  dxt,dyt,dzt,strpmf)
      
c***********************************************************************
c     
c     dlpoly constraint subroutine for potential of mean force calc.
c     accumulates constraint force to maintain reaction coordinate.
c     velocity verlet adaptation
c     
c     assume bond vectors dxp,dyp,dzp are input
c     dxp=(sum) wght*xxx(i,1) - (sum) wght*xxx(j,2) etc
c     
c     copyright daresbury laboratory
c     adapted by w.smith october 2005
c     original t.forester august 1995
c     
c***********************************************************************
      
      implicit none
      
      logical safep,img
      integer idnode,mxnode,imcon,natms,nspmf,icyc,k,jj
      integer i1,i,ipmf
      real(8) tstep,tolnce,virpmf,xxt,yyt,zzt,strpmf,summas,gamma
      real(8) dxp,dyp,dzp,xxa,yya,zza,amt,tstep2,dis,omega2,eps,gammi
      real(8) strs1,strs2,strs3,strs5,strs6,strs9,dxt,dyt,dzt,dsq
      
      dimension dxp(mspmf),dyp(mspmf),dzp(mspmf)
      dimension dxt(mxcons),dyt(mxcons),dzt(mxcons)
      dimension xxt(mxatms),yyt(mxatms),zzt(mxatms)
      dimension xxa(2,mspmf),yya(2,mspmf),zza(2,mspmf)
      dimension amt(2),strpmf(9),summas(2)
      
c     timestep squared
      
      tstep2=tstep*tstep
      
c     pmf virial
      
      virpmf=0.d0
      
c     accumulators for stress tensor
      
      strs1=0.d0
      strs2=0.d0
      strs3=0.d0
      strs5=0.d0
      strs6=0.d0
      strs9=0.d0
      
c     array bound check
      
      if(mxcons.lt.nspmf) call error(idnode,492)
      if(mspmf .lt.nspmf) call error(idnode,458)
      
c     application of constraint (shake) algorithm
      
      icyc=0
      img=.true.
      safep=.false.
      dis=prmpmf
      omega2=dis*dis
      amt(1)= tstep2/summas(1)
      amt(2)=-tstep2/summas(2)
      
      do while(.not.safep.and.icyc.lt.mxshak)
        
        icyc=icyc+1
        
        call pmf_vectors
     x    (img,nspmf,imcon,cell,xxx,yyy,zzz,xxt,yyt,zzt,xxa,yya,zza,
     x    dxt,dyt,dzt)
        
c     check convergence
        
        eps=0.d0
        
        do k=1,nspmf
          
          dsq=dxt(k)**2+dyt(k)**2+dzt(k)**2
          eps=max(eps,abs((omega2-dsq)/dis))
          
        enddo
        
        eps=eps*0.5d0
        
c     verification of convergence
        
        safep=(eps.lt.tolnce)
        
c     bypass calculations if converged
        
        if(.not.safep)then
          
          do k=1,nspmf
            
            dsq=dxt(k)**2+dyt(k)**2+dzt(k)**2
            gamma=(omega2-dsq)/(-(amt(2)-amt(1))*
     x        (dxp(k)*dxt(k)+dyp(k)*dyt(k)+dzp(k)*dzt(k)))
            
c     accumulate pmf virial
            
            virpmf=virpmf+gamma*(dxp(k)**2+dyp(k)**2+dzp(k)**2)
            
            strs1=strs1-gamma*dxp(k)*dxp(k)
            strs2=strs2-gamma*dxp(k)*dyp(k)
            strs3=strs3-gamma*dxp(k)*dzp(k)
            strs5=strs5-gamma*dyp(k)*dyp(k)
            strs6=strs6-gamma*dyp(k)*dzp(k)
            strs9=strs9-gamma*dzp(k)*dzp(k)
            
c     improve approximate atomic positions and velocities
            
            jj=0
            do ipmf=1,2
              
              gammi=-0.5d0*gamma*amt(ipmf)
              
              do i1=1,npmf(ipmf)
                
                jj=jj+1
                i=lstpmf(jj,k)
                
                xxx(i)=xxx(i)+dxp(k)*gammi
                yyy(i)=yyy(i)+dyp(k)*gammi
                zzz(i)=zzz(i)+dzp(k)*gammi
                vxx(i)=vxx(i)+dxp(k)*gammi/tstep
                vyy(i)=vyy(i)+dyp(k)*gammi/tstep
                vzz(i)=vzz(i)+dzp(k)*gammi/tstep
                
              enddo
              
            enddo
            
          enddo
          
        endif
        
      enddo
      
c     complete stress tensor
      
      strpmf(1)=strs1
      strpmf(2)=strs2
      strpmf(3)=strs3
      strpmf(4)=strs2
      strpmf(5)=strs5
      strpmf(6)=strs6
      strpmf(7)=strs3
      strpmf(8)=strs6
      strpmf(9)=strs9
      
c     splice coordinate and velocity arrays across nodes
      
      if(mxnode.gt.1)then
        
        buffer(1)=virpmf
        call gdsum(buffer(1),1,buffer(2))
        virpmf=buffer(1)
        call gdsum(strpmf,9,buffer)
        call splice 
     x    (idnode,natms,listpm,lstpmt,xxx,yyy,zzz,buffer)
        call splice 
     x    (idnode,natms,listpm,lstpmt,vxx,vyy,vzz,buffer)
        
      endif
      
      return
      end subroutine pmf_rattle_r
      
      subroutine pmf_rattle_v
     x  (safep,idnode,mxnode,imcon,natms,nspmf,tstep,tolnce,summas,
     x  dxp,dyp,dzp,vxt,vyt,vzt,vxa,vya,vza,vxp,vyp,vzp)
      
c***********************************************************************
c     
c     dlpoly constraint subroutine for potential of mean force calc.
c     accumulates velocity correction for second constraint condition
c     velocity verlet adaptation
c     
c     assume bond vectors dxp,dyp,dzp are input
c     dxp=(sum) wght*xxx(i,1) - (sum) wght*xxx(j,2) etc
c     
c     copyright daresbury laboratory
c     author w.smith october 2005
c     
c***********************************************************************
      
      implicit none
      
      logical safep,img
      integer idnode,mxnode,imcon,natms,nspmf,icyc,k,jj
      integer i1,i,ipmf
      real(8) tstep,tolnce,summas,gamma,vxt,vyt,vzt,vxa,vya,vza
      real(8) vxp,vyp,vzp,dxp,dyp,dzp,amt,omega,eps,gammi,tolvel
      
      dimension dxp(mspmf),dyp(mspmf),dzp(mspmf)
      dimension vxp(mspmf),vyp(mspmf),vzp(mspmf)
      dimension vxt(mxatms),vyt(mxatms),vzt(mxatms)
      dimension vxa(2,mspmf),vya(2,mspmf),vza(2,mspmf)
      dimension amt(2),summas(2)
      
c     constraint convergence tolerance
      
      tolvel=tolnce/tstep
      
c     array bound check
      
      if(mxcons.lt.nspmf) call error(idnode,492)
      if(mspmf .lt.nspmf) call error(idnode,458)
      
c     application of constraint (shake) algorithm
      
      icyc=0
      img=.false.
      safep=.false.
      amt(1)= 0.5d0*tstep/summas(1)
      amt(2)=-0.5d0*tstep/summas(2)
      
      do while(.not.safep.and.icyc.lt.mxshak)
        
        icyc=icyc+1
        
        call pmf_vectors
     x    (img,nspmf,imcon,cell,vxx,vyy,vzz,vxt,vyt,vzt,vxa,vya,vza,
     x    vxp,vyp,vzp)
        
c     check convergence
        
        eps=0.d0
        do k=1,nspmf
          
          omega=dxp(k)*vxp(k)+dyp(k)*vyp(k)+dzp(k)*vzp(k)
          eps=max(eps,abs(omega)/prmpmf)
          
        enddo
        
c     verification of convergence
        
        safep=(eps.lt.tolvel)
        
c     bypass calculations if converged
        
        if(.not.safep)then
          
          do k=1,nspmf
            
            omega=dxp(k)*vxp(k)+dyp(k)*vyp(k)+dzp(k)*vzp(k)
            gamma=omega/((amt(2)-amt(1))*
     x        (dxp(k)**2+dyp(k)**2+dzp(k)**2))
            
c     improve approximate atomic velocities
            
            jj=0
            do ipmf=1,2
              
              gammi=-gamma*amt(ipmf)
              
              do i1=1,npmf(ipmf)
                
                jj=jj+1
                i=lstpmf(jj,k)
                
                vxx(i)=vxx(i)+dxp(k)*gammi
                vyy(i)=vyy(i)+dyp(k)*gammi
                vzz(i)=vzz(i)+dzp(k)*gammi
                
              enddo
              
            enddo
            
          enddo
          
        endif
        
      enddo
      
c     splice velocity arrays across nodes
      
      if(mxnode.gt.1)call splice 
     x  (idnode,natms,listpm,lstpmt,vxx,vyy,vzz,buffer)
      
      return
      end subroutine pmf_rattle_v
      
      end module pmf_module
