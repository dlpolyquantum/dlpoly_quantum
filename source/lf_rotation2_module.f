      module lf_rotation2_module
      
c***********************************************************************
c     
c     dl_poly module 2 for verlet leap frog rotational integration 
c     schemes
c     copyright - daresbury laboratory
c     author    - w. smith    aug 2006
c     
c***********************************************************************
      
      use config_module
      use ensemble_tools_module
      use error_module
      use lf_rotation1_module
      use property_module
      use rigid_body_module
      use setup_module
      use shake_module
      use site_module
      use utility_module
      
      contains
      
      subroutine nveq_2
     x  (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x  ntcons,ntfree,engke,engrot,quattol,tolnce,tstep,vircom,
     x  vircon)
      
c***********************************************************************
c     
c     dlpoly routine to integrate rigid body equations of motion
c     using implicit leapfrog quaternion algorithm
c     plus bond constraints - rigid body sites and constraint sites 
c     may coincide.
c     
c     parallel replicated data version : block data
c     
c     tqx,tqy,tqz = torque in lab fixed frame (input)
c     omx,omy,omz = angular velocity in body fixed frame (principal axis)
c     rotinx,y,z  = rotational inertia in body fixed frame
c     
c     copyright daresbury laboratory 1995
c     author      t.forester june  1995
c     amended     w.smith nov 2005
c     
c**********************************************************************
      
      implicit none
      
      integer, parameter :: nnn=16
      real(8), parameter :: pt5=0.5d0
      
      logical safe,safeq,lshmov,newstep,newjob,cycle
      integer fail,imcon,idnode,mxnode,natms,ngrp,nscons
      integer ntcons,ntfree,igrp,igrp1,igrp2,idum,ifre,ifre1,ifre2
      integer i,j,k,jg,ig,jr,id,mxshak1,icyc
      real(8) engke,engrot,quattol,tolnce,tstep,vircom,vircon
      real(8) rot,strkin,strcon,strgrp,engfke,engtrn
      real(8) delx,dely,delz,trx,try,trz,vaa,vbb,vcc,viracc
      
      real(8), allocatable :: opx(:),opy(:),opz(:)
      real(8), allocatable :: oqx(:),oqy(:),oqz(:)
      real(8), allocatable :: dtx(:),dty(:),dtz(:)
      real(8), allocatable :: dxx(:),dyy(:),dzz(:)
      real(8), allocatable :: uxx(:),uyy(:),uzz(:)
      real(8), allocatable :: txx(:),tyy(:),tzz(:)
      real(8), allocatable :: xxt(:),yyt(:),zzt(:)
      real(8), allocatable :: dxt(:),dyt(:),dzt(:)
      real(8), allocatable :: xxo(:),yyo(:),zzo(:)
      real(8), allocatable :: vxo(:),vyo(:),vzo(:)
      real(8), allocatable :: fmx(:),fmy(:),fmz(:)
      real(8), allocatable :: tqx(:),tqy(:),tqz(:)
      real(8), allocatable :: omxo(:),omyo(:),omzo(:)
      real(8), allocatable :: gvxo(:),gvyo(:),gvzo(:)
      real(8), allocatable :: gcxo(:),gcyo(:),gczo(:)
      real(8), allocatable :: qn0(:),qn1(:),qn2(:),qn3(:)
      
      dimension rot(9),strkin(9),strcon(9),strgrp(9),fail(nnn)
      
      save newjob,igrp1,igrp2,ifre1,ifre2
      
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
      allocate (txx(mxatms),tyy(mxatms),tzz(mxatms),stat=fail(6))
      allocate (xxt(mxatms),yyt(mxatms),zzt(mxatms),stat=fail(7))
      allocate (dxt(mxcons),dyt(mxcons),dzt(mxcons),stat=fail(8))
      allocate (xxo(mxatms),yyo(mxatms),zzo(mxatms),stat=fail(9))
      allocate (vxo(msatms),vyo(msatms),vzo(msatms),stat=fail(10))
      allocate (fmx(msgrp),fmy(msgrp),fmz(msgrp),stat=fail(11))
      allocate (tqx(msgrp),tqy(msgrp),tqz(msgrp),stat=fail(12))
      allocate (omxo(msgrp),omyo(msgrp),omzo(msgrp),stat=fail(13))
      allocate (gvxo(msgrp),gvyo(msgrp),gvzo(msgrp),stat=fail(14))
      allocate (gcxo(msgrp),gcyo(msgrp),gczo(msgrp),stat=fail(15))
      allocate (qn0(msgrp),qn1(msgrp),qn2(msgrp),qn3(msgrp),
     x  stat=fail(16))
      do i=1,nnn
        if(fail(i).ne.0)call error(idnode,1620)
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
      cycle=.true.
      
c     store initial values of position and velocity
      
      do i=1,natms
        
        xxo(i)=xxx(i)
        yyo(i)=yyy(i)
        zzo(i)=zzz(i)
        
      enddo
      
      j=0
      do ifre=ifre1,ifre2
        
        j=j+1
        i=lstfre(ifre)
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
        gvxo(jg)=gvxx(ig)
        gvyo(jg)=gvyy(ig)
        gvzo(jg)=gvzz(ig)
        omxo(jg)=omx(ig)
        omyo(jg)=omy(ig)
        omzo(jg)=omz(ig)
        qn0(jg)=q0(ig)
        qn1(jg)=q1(ig)
        qn2(jg)=q2(ig)
        qn3(jg)=q3(ig)
        
      enddo
      
      if(ntcons.gt.0)then
        
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
        
      endif
      
c     calculate atom displacement from rigid body com
      
      jg=0
      jr=0
      do ig=igrp1,igrp2
        
        jg=jg+1
        id=lstgtp(ig)
        
        do j=1,numgsit(id)
          
          jr=jr+1
          i=lstrgd(jr)
          
          dtx(jr)=xxo(i)-gcxo(jg)
          dty(jr)=yyo(i)-gcyo(jg)
          dtz(jr)=zzo(i)-gczo(jg)
          
        enddo
        
      enddo
      
c     minimum images
      
      call images(imcon,0,1,jr,cell,dtx,dty,dtz)
      
c     accumulators for constraint and virial stress tensor
      
      vircon=0.d0
      do i=1,9
        strcns(i)=0.d0
      enddo
      
c     start of shake cycle
      
      icyc=0
      mxshak1=mxshak
      if(ntcons.eq.0)mxshak1=1
      do while(cycle.and.icyc.le.mxshak1)
        
        icyc=icyc+1
        
c     restore original quaternions for this step

        jg=0
        do ig=igrp1,igrp2

          jg=jg+1
          q0(ig)=qn0(jg)
          q1(ig)=qn1(jg)
          q2(ig)=qn2(jg)
          q3(ig)=qn3(jg)
          
        enddo

c     integrate 'free' particles
        
        j=0
        do ifre=ifre1,ifre2
          
          j=j+1
          i=lstfre(ifre)
          
c     advance velocity by leapfrog
          
          uxx(j)=vxo(j)+tstep*rmass(i)*fxx(i)
          uyy(j)=vyo(j)+tstep*rmass(i)*fyy(i)
          uzz(j)=vzo(j)+tstep*rmass(i)*fzz(i)
          
c     advance position by leapfrog
          
          xxx(i)=xxo(i)+tstep*uxx(j)
          yyy(i)=yyo(i)+tstep*uyy(j)
          zzz(i)=zzo(i)+tstep*uzz(j)
          
c     estimate full step velocities
          
          vxx(i)=pt5*(uxx(j)+vxo(j))
          vyy(i)=pt5*(uyy(j)+vyo(j))
          vzz(i)=pt5*(uzz(j)+vzo(j))
          
        enddo
        
c     calculate new kinetic energy at current timestep
        
        engfke=getkinf(ntfree,idnode,mxnode)
        
c     kinetic contribution to stress tensor
        
        call kinstressf(ntfree,idnode,mxnode,strkin)
        
c     restore half step velocities
        
        j=0
        do ifre=ifre1,ifre2
          
          j=j+1
          i=lstfre(ifre)
          vxx(i)=uxx(j)
          vyy(i)=uyy(j)
          vzz(i)=uzz(j)
          
        enddo
        
c     *************  Rigid body motion ****************************
        
        jg=0
        jr=0
        do ig=igrp1,igrp2
          
          jg=jg+1
          id=lstgtp(ig)
          
c     calculate centre of mass forces
          
          fmx(jg)=0.d0
          fmy(jg)=0.d0
          fmz(jg)=0.d0
          do j=1,numgsit(id)
            
            jr=jr+1
            i=lstrgd(jr)
            
            fmx(jg)=fmx(jg)+fxx(i)
            fmy(jg)=fmy(jg)+fyy(i)
            fmz(jg)=fmz(jg)+fzz(i)
            
          enddo
          
c     advance velocity by leapfrog
          
          uxx(jg)=gvxo(jg)+tstep/gmass(id)*fmx(jg)
          uyy(jg)=gvyo(jg)+tstep/gmass(id)*fmy(jg)
          uzz(jg)=gvzo(jg)+tstep/gmass(id)*fmz(jg)
          
c     advance position by leapfrog
          
          gcmx(ig)=gcxo(jg)+tstep*uxx(jg)
          gcmy(ig)=gcyo(jg)+tstep*uyy(jg)
          gcmz(ig)=gczo(jg)+tstep*uzz(jg)
          
c     centre of mass velocities at full step
          
          gvxx(ig)=pt5*(uxx(jg)+gvxo(jg))
          gvyy(ig)=pt5*(uyy(jg)+gvyo(jg))
          gvzz(ig)=pt5*(uzz(jg)+gvzo(jg))

        enddo
        
c     translational kinetic energy 
        
        engtrn=getkint(ngrp,idnode,mxnode)
        
c     kinetic contribution to stress tensor
        
        call kinstressg(ngrp,idnode,mxnode,strgrp)
        
c     restore half step velocity
        
        jg=0
        do ig=igrp1,igrp2
          
          jg=jg+1
          gvxx(ig)=uxx(jg)
          gvyy(ig)=uyy(jg)
          gvzz(ig)=uzz(jg)
          
        enddo
        
c     calculate rigid body contribution to stress tensor
        
        call bodystress(idnode,mxnode,ngrp,vircom,strbod,dtx,dty,dtz)
        
c     torques in lab frame
        
        jr=0
        jg=0
        do ig=igrp1,igrp2
          
          jg=jg+1
          id=lstgtp(ig)
          
          tqx(jg)=0.d0
          tqy(jg)=0.d0
          tqz(jg)=0.d0
          do j=1,numgsit(id)
            
            jr=jr+1
            i=lstrgd(jr)
            
            tqx(jg)=tqx(jg)+dty(jr)*fzz(i)-dtz(jr)*fyy(i)
            tqy(jg)=tqy(jg)+dtz(jr)*fxx(i)-dtx(jr)*fzz(i)
            tqz(jg)=tqz(jg)+dtx(jr)*fyy(i)-dty(jr)*fxx(i)
            
          enddo
          
c     current rotational matrix 
          
          call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)
          
c     angular velocity at time step n (first guess)
          
          opx(jg)=omxo(jg)
          opy(jg)=omyo(jg)
          opz(jg)=omzo(jg)
          
c     iterate angular velocity for time step n (e. yezdimer)
          
          do i=1,5
            
            trx=(tqx(jg)*rot(1)+tqy(jg)*rot(4)+tqz(jg)*rot(7))*
     x        rotinx(id,2)
     x        +(rotiny(id,1)-rotinz(id,1))*opy(jg)*opz(jg)*
     x        rotinx(id,2)
            try=(tqx(jg)*rot(2)+tqy(jg)*rot(5)+tqz(jg)*rot(8))*
     x        rotiny(id,2)
     x        +(rotinz(id,1)-rotinx(id,1))*opz(jg)*opx(jg)*
     x        rotiny(id,2)
            trz=(tqx(jg)*rot(3)+tqy(jg)*rot(6)+tqz(jg)*rot(9))*
     x        rotinz(id,2)
     x        +(rotinx(id,1)-rotiny(id,1))*opx(jg)*opy(jg)*
     x        rotinz(id,2)
            
c     improved angular velocity at time step n
            
            delx=tstep*trx
            dely=tstep*try
            delz=tstep*trz
            
            opx(jg)=omxo(jg)+delx*pt5
            opy(jg)=omyo(jg)+dely*pt5
            opz(jg)=omzo(jg)+delz*pt5
            
          enddo
          
c     angular velocity at time step n+1/2
          
          uxx(jg)=omxo(jg)+delx
          uyy(jg)=omyo(jg)+dely
          uzz(jg)=omzo(jg)+delz
          
c     angular velocity at time step n+1 (needed for quat algorithm)
          
          oqx(jg)=omxo(jg)+delx*1.5d0
          oqy(jg)=omyo(jg)+dely*1.5d0
          oqz(jg)=omzo(jg)+delz*1.5d0
          
c     angular velocity at time step n
          
          omx(ig)=opx(jg)
          omy(ig)=opy(jg)
          omz(ig)=opz(jg)
          
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
          omx(ig)=uxx(jg)
          omy(ig)=uyy(jg)
          omz(ig)=uzz(jg)
          
        enddo
        
c     assign new quaternions
        
        call update_quaternions
     x    (safeq,igrp1,igrp2,tstep,quattol,opx,opy,opz,oqx,oqy,oqz)
        
c     minimum images of group positions and particle positions
        
        call images(imcon,0,1,natms,cell,xxx,yyy,zzz)
        call images(imcon,idnode,mxnode,ngrp,cell,gcmx,gcmy,gcmz)
        
c     new atomic positions for atoms in rigid bodies
        
        jr=0
        jg=0
        do ig=igrp1,igrp2
          
          jg=jg+1
          id=lstgtp(ig)
          
c     new rotational matrix
          
          call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)
          
          do j=1,numgsit(id)
            
            jr=jr+1
            i=lstrgd(jr)
            
            xxx(i)=rot(1)*gxx(id,j)+rot(2)*gyy(id,j)+rot(3)*gzz(id,j)
     x        +gcmx(ig)
            yyy(i)=rot(4)*gxx(id,j)+rot(5)*gyy(id,j)+rot(6)*gzz(id,j)
     x        +gcmy(ig)
            zzz(i)=rot(7)*gxx(id,j)+rot(8)*gyy(id,j)+rot(9)*gzz(id,j)
     x        +gcmz(ig)
            
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
        
c     merge new atomic coordinates 
        
        if(mxnode.gt.1)call merge1
     x    (idnode,mxnode,natms,lstme,xxx,yyy,zzz,buffer)
        
        if(ntcons.gt.0) then
          
c     apply constraint correction
          
          newstep=.false.
          if(icyc.eq.1)newstep=.true.
          
          call qshake
     x      (newstep,safe,lshmov,idnode,imcon,mxnode,natms,
     x      nscons,tolnce,tstep,viracc,dxx,dyy,dzz,dxt,dyt,
     x      dzt,txx,tyy,tzz,xxt,yyt,zzt,strcon)
          
          if(abs(viracc).le.1.d-10)cycle=.false.
          
          vircon=vircon+viracc
          do i=1,9
            strcns(i)=strcns(i)+strcon(i)
          enddo
          
c     end of shake corrections
          
        endif
        
      enddo
      
      if(mxnode.gt.1) then
        
c     merge new group coordinates and velocities
        
        call merge(idnode,mxnode,ngrp,mxbuff,gcmx,gcmy,gcmz,buffer)
        call merge(idnode,mxnode,ngrp,mxbuff,gvxx,gvyy,gvzz,buffer)
        
c     merge new atomic velocities
        
        call merge1(idnode,mxnode,natms,lstme,xxx,yyy,zzz,buffer)
        call merge1(idnode,mxnode,natms,lstme,vxx,vyy,vzz,buffer)
        
c     merge new quaternions
        
        call merge4(idnode,mxnode,ngrp,mxbuff,q0,q1,q2,q3,buffer)
        
      endif
      
c     total kinetic energy
      
      engke=engfke+engtrn
      
c     ensure all atoms are within cell boundaries
      
      call images(imcon,0,1,natms,cell,xxx,yyy,zzz)
      
c     complete stress tensor
      
      do i=1,9
        stress(i)=stress(i)+strkin(i)+strgrp(i)+strcns(i)+strbod(i)
      enddo
      
c     deallocate work arrays
      
      deallocate (opx,opy,opz,oqx,oqy,oqz,stat=fail(1))
      deallocate (dtx,dty,dtz,dxx,dyy,dzz,stat=fail(2))
      deallocate (uxx,uyy,uzz,txx,tyy,tzz,stat=fail(3))
      deallocate (dxt,dyt,dzt,xxo,yyo,zzo,stat=fail(4))
      deallocate (vxo,vyo,vzo,fmx,fmy,fmz,stat=fail(5))
      deallocate (tqx,tqy,tqz,omxo,omyo,omzo,stat=fail(6))
      deallocate (gvxo,gvyo,gvzo,gcxo,gcyo,gczo,stat=fail(7))
      deallocate (xxt,yyt,zzt,qn0,qn1,qn2,qn3,stat=fail(8))
      
      return
      end subroutine nveq_2
      
      subroutine nvtq_b2
     x  (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x  ntcons,ntfree,engke,engrot,quattol,sigma,taut,tolnce,
     x  tstep,vircom,vircon)
      
c***********************************************************************
c     
c     dlpoly routine to integrate rigid body equations of motion
c     using implicit leapfrog quaternion algorithm
c     plus bond constraints - rigid body sites and constraint sites
c     may coincide.
c     
c     verlet leapfrog with Berendsen thermostat.
c
c     parallel replicated data version : block data
c     
c     tqx,tqy,tqz=torque in lab fixed frame (input)
c     omx,omy,omz=angular velocity in body fixed frame (principal axis)
c     rotinx,y,z =rotational inertia in body fixed frame
c     
c     copyright daresbury laboratory 1995
c     author      t.forester june 1995
c     amended     w.smith nov 2005
c     
c**********************************************************************
      
      implicit none
      
      integer, parameter :: nnn=16
      real(8), parameter :: pt5=0.5d0
      
      logical safe,safeq,lshmov,newstep,newjob,cycle
      integer fail,imcon,idnode,mxnode,natms,ngrp,nscons
      integer ntcons,ntfree,igrp,igrp1,igrp2,idum,ifre,ifre1,ifre2
      integer i,j,k,jg,ig,jr,id,mxshak1,icyc
      real(8) engke,engrot,quattol,tolnce,tstep,vircom,vircon,engfke
      real(8) rot,strkin,strgrp,strcon,engtrn,trx,try,trz
      real(8) delx,dely,delz,vaa,vbb,vcc,viracc,engtot,chit0
      real(8) sigma,taut
      
      real(8), allocatable :: opx(:),opy(:),opz(:)
      real(8), allocatable :: oqx(:),oqy(:),oqz(:)
      real(8), allocatable :: dtx(:),dty(:),dtz(:)
      real(8), allocatable :: dxx(:),dyy(:),dzz(:)
      real(8), allocatable :: uxx(:),uyy(:),uzz(:)
      real(8), allocatable :: txx(:),tyy(:),tzz(:)
      real(8), allocatable :: xxt(:),yyt(:),zzt(:)
      real(8), allocatable :: dxt(:),dyt(:),dzt(:)
      real(8), allocatable :: xxo(:),yyo(:),zzo(:)
      real(8), allocatable :: vxo(:),vyo(:),vzo(:)
      real(8), allocatable :: fmx(:),fmy(:),fmz(:)
      real(8), allocatable :: tqx(:),tqy(:),tqz(:)
      real(8), allocatable :: omxo(:),omyo(:),omzo(:)
      real(8), allocatable :: gvxo(:),gvyo(:),gvzo(:)
      real(8), allocatable :: gcxo(:),gcyo(:),gczo(:)
      real(8), allocatable :: qn0(:),qn1(:),qn2(:),qn3(:)
      
      dimension rot(9),strkin(9),strcon(9),strgrp(9),fail(nnn)
      
      save chit0,igrp1,igrp2,ifre1,ifre2,newjob
      
      data chit0/1.d0/,newjob/.true./
      
c     allocate working arrays
      
      do i=1,nnn
        fail(i)=0
      enddo
      allocate (opx(msgrp),opy(msgrp),opz(msgrp),stat=fail(1))
      allocate (oqx(msgrp),oqy(msgrp),oqz(msgrp),stat=fail(2))
      allocate (dtx(mxatms),dty(mxatms),dtz(mxatms),stat=fail(3))
      allocate (dxx(mxcons),dyy(mxcons),dzz(mxcons),stat=fail(4))
      allocate (uxx(mxatms),uyy(mxatms),uzz(mxatms),stat=fail(5))
      allocate (txx(mxatms),tyy(mxatms),tzz(mxatms),stat=fail(6))
      allocate (xxt(mxatms),yyt(mxatms),zzt(mxatms),stat=fail(7))
      allocate (dxt(mxcons),dyt(mxcons),dzt(mxcons),stat=fail(8))
      allocate (xxo(mxatms),yyo(mxatms),zzo(mxatms),stat=fail(9))
      allocate (vxo(msatms),vyo(msatms),vzo(msatms),stat=fail(10))
      allocate (fmx(msgrp),fmy(msgrp),fmz(msgrp),stat=fail(11))
      allocate (tqx(msgrp),tqy(msgrp),tqz(msgrp),stat=fail(12))
      allocate (omxo(msgrp),omyo(msgrp),omzo(msgrp),stat=fail(13))
      allocate (gvxo(msgrp),gvyo(msgrp),gvzo(msgrp),stat=fail(14))
      allocate (gcxo(msgrp),gcyo(msgrp),gczo(msgrp),stat=fail(15))
      allocate (qn0(msgrp),qn1(msgrp),qn2(msgrp),qn3(msgrp),
     x  stat=fail(16))
      do i=1,nnn
        if(fail(i).ne.0)call error(idnode,1630)
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
      cycle=.true.
      
c     store initial values of position and velocity
      
      do i=1,natms
        
        xxo(i)=xxx(i)
        yyo(i)=yyy(i)
        zzo(i)=zzz(i)
        
      enddo
      
      j=0
      do ifre=ifre1,ifre2
        
        j=j+1
        i=lstfre(ifre)
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
        gvxo(jg)=gvxx(ig)
        gvyo(jg)=gvyy(ig)
        gvzo(jg)=gvzz(ig)
        omxo(jg)=omx(ig)
        omyo(jg)=omy(ig)
        omzo(jg)=omz(ig)
        qn0(jg)=q0(ig)
        qn1(jg)=q1(ig)
        qn2(jg)=q2(ig)
        qn3(jg)=q3(ig)
        
      enddo
      
c     construct current bond vectors
      
      if(ntcons.gt.0)then
        
        do k=1,nscons
          
          i=listcon(k,2)
          j=listcon(k,3)
          
          dxx(k)=xxx(i)-xxx(j)
          dyy(k)=yyy(i)-yyy(j)
          dzz(k)=zzz(i)-zzz(j)
          
        enddo
        
c     periodic boundary condition for bond vectors
        
        call images(imcon,0,1,nscons,cell,dxx,dyy,dzz)
        
      endif
      
c     calculate atom displacement from rigid body com
      
      jg=0
      jr=0
      do ig=igrp1,igrp2
        
        jg=jg+1
        id=lstgtp(ig)
        
        do j=1,numgsit(id)
          
          jr=jr+1
          i=lstrgd(jr)
          
          dtx(jr)=xxo(i)-gcxo(jg)
          dty(jr)=yyo(i)-gcyo(jg)
          dtz(jr)=zzo(i)-gczo(jg)
          
        enddo
        
      enddo
      
c     minimum images
        
        call images(imcon,0,1,jr,cell,dtx,dty,dtz)
        
c     accumulators for constraint stress and virial
      
      vircon=0.d0
      do i=1,9
        strcns(i)=0.d0
      enddo
      
c     shake and thermostat iterations start here
      
      icyc=0
      mxshak1=mxshak
      if(ntcons.eq.0)mxshak1=2
      do while(cycle.and.icyc.le.mxshak1)
        
        icyc=icyc+1
        
c     restore original quaternions for this step

        jg=0
        do ig=igrp1,igrp2

          jg=jg+1
          q0(ig)=qn0(jg)
          q1(ig)=qn1(jg)
          q2(ig)=qn2(jg)
          q3(ig)=qn3(jg)
          
        enddo

c     integrate 'free' particles
        
        j=0
        do ifre=ifre1,ifre2
          
          j=j+1
          i=lstfre(ifre)
          
c     advance velocity by leapfrog
          
          uxx(j)=(vxo(j)+tstep*rmass(i)*fxx(i))*chit0
          uyy(j)=(vyo(j)+tstep*rmass(i)*fyy(i))*chit0
          uzz(j)=(vzo(j)+tstep*rmass(i)*fzz(i))*chit0
          
c     advance position by leapfrog
          
          xxx(i)=xxo(i)+tstep*uxx(j)
          yyy(i)=yyo(i)+tstep*uyy(j)
          zzz(i)=zzo(i)+tstep*uzz(j)
          
c     calculate full time velocity
          
          vxx(i)=pt5*(uxx(j)+vxo(j))
          vyy(i)=pt5*(uyy(j)+vyo(j))
          vzz(i)=pt5*(uzz(j)+vzo(j))
          
        enddo
        
c     calculate kinetic energy at current timestep
        
        engfke=getkinf(ntfree,idnode,mxnode)
        
c     calculate kinetic stress tensor
        
        call kinstressf(ntfree,idnode,mxnode,strkin)
        
c     restore half step velocity
        
        j=0
        do ifre=ifre1,ifre2
          
          j=j+1
          i=lstfre(ifre)
          vxx(i)=uxx(j)
          vyy(i)=uyy(j)
          vzz(i)=uzz(j)
          
        enddo
        
c     *************  Rigid body motion ****************************
        
c     translational kinetic energy
        
        jg=0
        jr=0
        do ig=igrp1,igrp2
          
          jg=jg+1
          id=lstgtp(ig)
          
c     calculate centre of mass forces
          
          fmx(jg)=0.d0
          fmy(jg)=0.d0
          fmz(jg)=0.d0
          do j=1,numgsit(id)
            
            jr=jr+1
            i=lstrgd(jr)
            
            fmx(jg)=fmx(jg)+fxx(i)
            fmy(jg)=fmy(jg)+fyy(i)
            fmz(jg)=fmz(jg)+fzz(i)
            
          enddo
          
c     update centre of mass velocities
          
          gvxx(ig)=gvxo(jg)+tstep/gmass(id)*fmx(jg)
          gvyy(ig)=gvyo(jg)+tstep/gmass(id)*fmy(jg)
          gvzz(ig)=gvzo(jg)+tstep/gmass(id)*fmz(jg)
          
        enddo
        
c     translational kinetic energy 
        
        engtrn=getkint(ngrp,idnode,mxnode)
        
c     kinetic stress tensor
        
        call kinstressg(ngrp,idnode,mxnode,strgrp)
        
        jg=0
        do ig=igrp1,igrp2
          
          jg=jg+1
          
c     advance velocity by leapfrog
          
          gvxx(ig)=(gvxo(jg)+tstep/gmass(id)*fmx(jg))*chit0
          gvyy(ig)=(gvyo(jg)+tstep/gmass(id)*fmy(jg))*chit0
          gvzz(ig)=(gvzo(jg)+tstep/gmass(id)*fmz(jg))*chit0
          
c     advance position by leapfrog
          
          gcmx(ig)=gcxo(jg)+tstep*gvxx(ig)
          gcmy(ig)=gcyo(jg)+tstep*gvyy(ig)
          gcmz(ig)=gczo(jg)+tstep*gvzz(ig)
          
        enddo
        
c     calculate rigid body stress tensor
        
        call bodystress(idnode,mxnode,ngrp,vircom,strbod,dtx,dty,dtz)
        
c     torques in lab frame
        
        jr=0
        jg=0
        do ig=igrp1,igrp2
          
          jg=jg+1
          id=lstgtp(ig)
          
          tqx(jg)=0.d0
          tqy(jg)=0.d0
          tqz(jg)=0.d0
          do j=1,numgsit(id)
            
            jr=jr+1
            i=lstrgd(jr)
            
            tqx(jg)=tqx(jg)+dty(jr)*fzz(i)-dtz(jr)*fyy(i)
            tqy(jg)=tqy(jg)+dtz(jr)*fxx(i)-dtx(jr)*fzz(i)
            tqz(jg)=tqz(jg)+dtx(jr)*fyy(i)-dty(jr)*fxx(i)
            
          enddo
          
c     current rotational matrix 
          
          call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)
          
c     angular velocity at time step n (first guess)
          
          opx(jg)=omxo(jg)
          opy(jg)=omyo(jg)
          opz(jg)=omzo(jg)
          
c     iterate angular velocity for time step n (e. yezdimer)
          
          do i=1,5
            
            trx=(tqx(jg)*rot(1)+tqy(jg)*rot(4)+tqz(jg)*rot(7))*
     x        rotinx(id,2)
     x        +(rotiny(id,1)-rotinz(id,1))*opy(jg)*opz(jg)*
     x        rotinx(id,2)
            try=(tqx(jg)*rot(2)+tqy(jg)*rot(5)+tqz(jg)*rot(8))*
     x        rotiny(id,2)
     x        +(rotinz(id,1)-rotinx(id,1))*opz(jg)*opx(jg)*
     x        rotiny(id,2)
            trz=(tqx(jg)*rot(3)+tqy(jg)*rot(6)+tqz(jg)*rot(9))*
     x        rotinz(id,2)
     x        +(rotinx(id,1)-rotiny(id,1))*opx(jg)*opy(jg)*
     x        rotinz(id,2)
            
c     improved angular velocity at time step n
            
            delx=tstep*trx
            dely=tstep*try
            delz=tstep*trz
            
            opx(jg)=(omxo(jg)+delx*pt5)
            opy(jg)=(omyo(jg)+dely*pt5)
            opz(jg)=(omzo(jg)+delz*pt5)
            
          enddo
          
c     scaled angular velocity at timestep n
          
          omx(ig)=opx(jg)*chit0
          omy(ig)=opy(jg)*chit0
          omz(ig)=opz(jg)*chit0
          
c     angular velocity at time step n+1 (needed for quat algorithm)
          
          oqx(jg)=(omxo(jg)+delx*1.5d0)*chit0
          oqy(jg)=(omyo(jg)+dely*1.5d0)*chit0
          oqz(jg)=(omzo(jg)+delz*1.5d0)*chit0
          
c     angular velocity at full time step
          
          uxx(jg)=(omxo(jg)+delx)*chit0
          uyy(jg)=(omyo(jg)+dely)*chit0
          uzz(jg)=(omzo(jg)+delz)*chit0
          
        enddo
        
c     rotational kinetic energy
        
        engrot=getkinr(ngrp,idnode,mxnode)
        
c     restore half step angular velocities
        
        jg=0
        do ig=igrp1,igrp2
          
          jg=jg+1
          opx(jg)=omx(ig)
          opy(jg)=omy(ig)
          opz(jg)=omz(ig)
          omx(ig)=uxx(jg)
          omy(ig)=uyy(jg)
          omz(ig)=uzz(jg)
          
        enddo
        
c     assign new quaternions
        
        call update_quaternions
     x    (safeq,igrp1,igrp2,tstep,quattol,opx,opy,opz,oqx,oqy,oqz)
        
c     new atomic positions for atoms in rigid bodies
        
        jg=0
        jr=0
        do ig=igrp1,igrp2
          
          jg=jg+1
          id=lstgtp(ig)
          
c     new rotational matrix
          
          call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)
          
          do j=1,numgsit(id)
            
            jr=jr+1
            i=lstrgd(jr)
            
            xxx(i)=rot(1)*gxx(id,j)+rot(2)*gyy(id,j)+rot(3)*gzz(id,j)
     x        +gcmx(ig)
            yyy(i)=rot(4)*gxx(id,j)+rot(5)*gyy(id,j)+rot(6)*gzz(id,j)
     x        +gcmy(ig)
            zzz(i)=rot(7)*gxx(id,j)+rot(8)*gyy(id,j)+rot(9)*gzz(id,j)
     x        +gcmz(ig)
            
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
        
c     merge new atomic coordinates 
        
        if(mxnode.gt.1)call merge1
     x    (idnode,mxnode,natms,lstme,xxx,yyy,zzz,buffer)
        
c     new estimate of chit
        
        engke=engfke+engtrn
        engtot=engke+engrot 
        chit0=sqrt(1.d0+tstep/taut*(sigma/engtot-1.d0))
        
        if(ntcons.gt.0) then
          
c     apply constraint correction
          
          newstep=.false.
          if(icyc.eq.1)newstep=.true.
          
          call qshake
     x      (newstep,safe,lshmov,idnode,imcon,mxnode,natms,
     x      nscons,tolnce,tstep,viracc,dxx,dyy,dzz,dxt,dyt,
     x      dzt,txx,tyy,tzz,xxt,yyt,zzt,strcon)
          
          if(abs(viracc).le.1.d-10)cycle=.false.
          
          vircon=vircon+viracc
          do i=1,9
            strcns(i)=strcns(i)+strcon(i)
          enddo
          
c     end of shake corrections
          
        endif
        
      enddo
      
c     minimum images of group positions and particle positions
      
      call images(imcon,0,1,natms,cell,xxx,yyy,zzz)
      call images(imcon,idnode,mxnode,ngrp,cell,gcmx,gcmy,gcmz)
      
c     new atomic positions for atoms in rigid bodies
      
      jr=0
      jg=0
      do ig=igrp1,igrp2
        
        jg=jg+1
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
      
      
      if(mxnode.gt.1)then
        
c     merge new group coordinates and velocities
        
        call merge(idnode,mxnode,ngrp,mxbuff,gcmx,gcmy,gcmz,buffer)
        call merge(idnode,mxnode,ngrp,mxbuff,gvxx,gvyy,gvzz,buffer)
        
c     merge new atomic velocities
        
        call merge1(idnode,mxnode,natms,lstme,vxx,vyy,vzz,buffer)
        
c     merge new quaternions
        
        call merge4(idnode,mxnode,ngrp,mxbuff,q0,q1,q2,q3,buffer)
        
      endif
      
c     ensure all atoms are within cell boundaries
      
      call images(imcon,0,1,natms,cell,xxx,yyy,zzz)
      
c     complete stress tensor
      
      do i=1,9
        stress(i)=stress(i)+strcns(i)+strkin(i)+strgrp(i)+strbod(i)
      enddo
      
c     deallocate work arrays
      
      deallocate (opx,opy,opz,oqx,oqy,oqz,stat=fail(1))
      deallocate (dtx,dty,dtz,dxx,dyy,dzz,stat=fail(2))
      deallocate (uxx,uyy,uzz,txx,tyy,tzz,stat=fail(3))
      deallocate (dxt,dyt,dzt,xxo,yyo,zzo,stat=fail(4))
      deallocate (vxo,vyo,vzo,tqx,tqy,tqz,stat=fail(5))
      deallocate (fmx,fmy,fmz,omxo,omyo,omzo,stat=fail(6))
      deallocate (gvxo,gvyo,gvzo,gcxo,gcyo,gczo,stat=fail(7))
      deallocate (xxt,yyt,zzt,qn0,qn1,qn2,qn3,stat=fail(8))
      
      return
      end subroutine nvtq_b2
      
      subroutine nvtq_h2
     x  (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x  ntcons,ntfree,conint,consv,chit,engke,engrot,quattol,
     x  sigma,taut,tolnce,tstep,vircom,vircon)
      
c***********************************************************************
c     
c     dlpoly routine to integrate rigid body equations of motion
c     using implicit leapfrog quaternion algorithm
c     plus bond constraints - rigid body sites and constraint sites
c     may coincide.
c     
c     verlet leapfrog with Hoover thermostat.
c
c     parallel replicated data version : block data
c     
c     tqx,tqy,tqz = torque in lab fixed frame (input)
c     omx,omy,omz = angular velocity in body fixed frame (principal axis)
c     rotinx,y,z  = rotational inertia in body fixed frame
c     
c     copyright daresbury laboratory 1995
c     author      t.forester june 1995
c     amended     w.smith nov 2005
c     
c**********************************************************************
      
      implicit none
      
      integer, parameter :: nnn=16
      real(8), parameter :: pt5=0.5d0
      
      logical safe,safeq,lshmov,newstep,newjob,cycle
      integer fail,imcon,idnode,mxnode,natms,ngrp,nscons
      integer ntcons,ntfree,igrp,igrp1,igrp2,idum,ifre,ifre1,ifre2
      integer i,j,k,jg,ig,jr,id,mxshak1,icyc
      real(8) engke,engrot,quattol,tolnce,tstep,vircom,vircon,engtot
      real(8) rot,strkin,strcon,strgrp,engtrn,engfke,trx,try,trz
      real(8) delx,dely,delz,vaa,vbb,vcc,viracc,sigma,taut,chit0
      real(8) chitnew,chitp,conint,consv,chit,qmass
      
      real(8), allocatable :: opx(:),opy(:),opz(:)
      real(8), allocatable :: oqx(:),oqy(:),oqz(:)
      real(8), allocatable :: dtx(:),dty(:),dtz(:)
      real(8), allocatable :: dxx(:),dyy(:),dzz(:)
      real(8), allocatable :: uxx(:),uyy(:),uzz(:)
      real(8), allocatable :: txx(:),tyy(:),tzz(:)
      real(8), allocatable :: xxt(:),yyt(:),zzt(:)
      real(8), allocatable :: dxt(:),dyt(:),dzt(:)
      real(8), allocatable :: xxo(:),yyo(:),zzo(:)
      real(8), allocatable :: vxo(:),vyo(:),vzo(:)
      real(8), allocatable :: fmx(:),fmy(:),fmz(:)
      real(8), allocatable :: tqx(:),tqy(:),tqz(:)
      real(8), allocatable :: omxo(:),omyo(:),omzo(:)
      real(8), allocatable :: gvxo(:),gvyo(:),gvzo(:)
      real(8), allocatable :: gcxo(:),gcyo(:),gczo(:)
      real(8), allocatable :: qn0(:),qn1(:),qn2(:),qn3(:)
      
      dimension rot(9),strkin(9),strcon(9),strgrp(9),fail(nnn)
      
      save newjob,igrp1,igrp2,ifre1,ifre2,qmass
      
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
      allocate (txx(mxatms),tyy(mxatms),tzz(mxatms),stat=fail(6))
      allocate (xxt(mxatms),yyt(mxatms),zzt(mxatms),stat=fail(7))
      allocate (dxt(mxcons),dyt(mxcons),dzt(mxcons),stat=fail(8))
      allocate (xxo(mxatms),yyo(mxatms),zzo(mxatms),stat=fail(9))
      allocate (vxo(msatms),vyo(msatms),vzo(msatms),stat=fail(10))
      allocate (fmx(msgrp),fmy(msgrp),fmz(msgrp),stat=fail(11))
      allocate (tqx(msgrp),tqy(msgrp),tqz(msgrp),stat=fail(12))
      allocate (omxo(msgrp),omyo(msgrp),omzo(msgrp),stat=fail(13))
      allocate (gvxo(msgrp),gvyo(msgrp),gvzo(msgrp),stat=fail(14))
      allocate (gcxo(msgrp),gcyo(msgrp),gczo(msgrp),stat=fail(15))
      allocate (qn0(msgrp),qn1(msgrp),qn2(msgrp),qn3(msgrp),
     x  stat=fail(16))
      do i=1,nnn
        if(fail(i).ne.0)call error(idnode,1640)
      enddo
      
      if(newjob)then
        
c     mass parameters for thermostat
        
        qmass=2.d0*sigma*taut**2
        
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
        
      endif
      
      safe=.false.
      cycle=.true.
      
c     store initial values of position and velocity
      
      do i=1,natms
        
        xxo(i)=xxx(i)
        yyo(i)=yyy(i)
        zzo(i)=zzz(i)
        
      enddo
      
      j=0
      do ifre=ifre1,ifre2
        
        j=j+1
        i=lstfre(ifre)
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
        gvxo(jg)=gvxx(ig)
        gvyo(jg)=gvyy(ig)
        gvzo(jg)=gvzz(ig)
        omxo(jg)=omx(ig)
        omyo(jg)=omy(ig)
        omzo(jg)=omz(ig)
        qn0(jg)=q0(ig)
        qn1(jg)=q1(ig)
        qn2(jg)=q2(ig)
        qn3(jg)=q3(ig)
        
      enddo
      
c     construct current bond vectors
      
      if(ntcons.gt.0)then
        
        do k=1,nscons
          
          i=listcon(k,2)
          j=listcon(k,3)
          
          dxx(k)=xxx(i)-xxx(j)
          dyy(k)=yyy(i)-yyy(j)
          dzz(k)=zzz(i)-zzz(j)
          
        enddo
        
c     periodic boundary condition for bond vectors
        
        call images(imcon,0,1,nscons,cell,dxx,dyy,dzz)
        
      endif
      
c     calculate atoms displacement from rigid body com
      
      jg=0
      jr=0
      do ig=igrp1,igrp2
        
        jg=jg+1
        id=lstgtp(ig)
        
        do j=1,numgsit(id)
          
          jr=jr+1
          i=lstrgd(jr)
          
          dtx(jr)=xxo(i)-gcxo(jg)
          dty(jr)=yyo(i)-gcyo(jg)
          dtz(jr)=zzo(i)-gczo(jg)
          
        enddo
        
      enddo
      
c     minimum images
      
      call images(imcon,0,1,jr,cell,dtx,dty,dtz)
      
c     initial thermostat variable
      
      chit0=chit
      
c     accumulators for constraint stress and virial
      
      vircon=0.d0
      do i=1,9
        strcns(i)=0.d0
      enddo
      
c     shake iterations and thermostat iterations start here
      
      icyc=0
      mxshak1=mxshak
      if(ntcons.eq.0)mxshak1=3
      do while(cycle.and.icyc.le.mxshak1)
        
        icyc=icyc+1
        
c     restore original quaternions for this step

        jg=0
        do ig=igrp1,igrp2

          jg=jg+1
          q0(ig)=qn0(jg)
          q1(ig)=qn1(jg)
          q2(ig)=qn2(jg)
          q3(ig)=qn3(jg)
          
        enddo

c     integrate 'free' particles
        
        j=0
        do ifre=ifre1,ifre2
          
          j=j+1
          i=lstfre(ifre)
          
c     advance velocity by leapfrog
          
          uxx(j)=vxo(j)+tstep*(fxx(i)*rmass(i)-(chit0)*
     x      pt5*(vxx(i)+vxo(j)))
          uyy(j)=vyo(j)+tstep*(fyy(i)*rmass(i)-(chit0)*
     x      pt5*(vyy(i)+vyo(j)))
          uzz(j)=vzo(j)+tstep*(fzz(i)*rmass(i)-(chit0)*
     x      pt5*(vzz(i)+vzo(j)))
          
c     advance position by leapfrog
          
          xxx(i)=xxo(i)+tstep*uxx(j)
          yyy(i)=yyo(i)+tstep*uyy(j)
          zzz(i)=zzo(i)+tstep*uzz(j)
          
c     estimate full step velocities
          
          vxx(i)=pt5*(vxo(j)+uxx(j))
          vyy(i)=pt5*(vyo(j)+uyy(j))
          vzz(i)=pt5*(vzo(j)+uzz(j))
          
        enddo
        
c     calculate new kinetic energy at current timestep
        
        engfke=getkinf(ntfree,idnode,mxnode)
        
c     calculate kinetic stress tensor
        
        call kinstressf(ntfree,idnode,mxnode,strkin)
        
c     restore half step velocity
        
        j=0
        do ifre=ifre1,ifre2
          
          j=j+1
          i=lstfre(ifre)
          vxx(i)=uxx(j)
          vyy(i)=uyy(j)
          vzz(i)=uzz(j)
          
        enddo
        
c     *************  Rigid body motion ****************************
        
c     translational kinetic energy
        
        jg=0
        jr=0
        do ig=igrp1,igrp2
          
          jg=jg+1
          id=lstgtp(ig)
          
c     calculate centre of mass forces
          
          fmx(jg)=0.d0
          fmy(jg)=0.d0
          fmz(jg)=0.d0
          do j=1,numgsit(id)
            
            jr=jr+1
            i=lstrgd(jr)
            
            fmx(jg)=fmx(jg)+fxx(i)
            fmy(jg)=fmy(jg)+fyy(i)
            fmz(jg)=fmz(jg)+fzz(i)
            
          enddo
          
c     advance velocity by leapfrog
          
          uxx(jg)=gvxo(jg)+tstep*(fmx(jg)/gmass(id)-chit0*
     x      pt5*(gvxx(ig)+gvxo(jg)))
          uyy(jg)=gvyo(jg)+tstep*(fmy(jg)/gmass(id)-chit0*
     x      pt5*(gvyy(ig)+gvyo(jg)))
          uzz(jg)=gvzo(jg)+tstep*(fmz(jg)/gmass(id)-chit0*
     x      pt5*(gvzz(ig)+gvzo(jg)))
          
c     advance position by leapfrog
          
          gcmx(ig)=gcxo(jg)+tstep*uxx(jg)
          gcmy(ig)=gcyo(jg)+tstep*uyy(jg)
          gcmz(ig)=gczo(jg)+tstep*uzz(jg)
          
c     centre of mass velocities at half-step
          
          gvxx(ig)=pt5*(gvxo(jg)+uxx(jg))
          gvyy(ig)=pt5*(gvyo(jg)+uyy(jg))
          gvzz(ig)=pt5*(gvzo(jg)+uzz(jg))
          
        enddo
        
c     translational kinetic energy 
        
        engtrn=getkint(ngrp,idnode,mxnode)
        
c     kinetic contribution to stress tensor
        
        call kinstressg(ngrp,idnode,mxnode,strgrp)
        
c     restore half step velocity
        
        jg=0
        do ig=igrp1,igrp2
          
          jg=jg+1
          gvxx(ig)=uxx(jg)
          gvyy(ig)=uyy(jg)
          gvzz(ig)=uzz(jg)
          
        enddo
        
c     calculate rigid body stress tensor
        
        call bodystress(idnode,mxnode,ngrp,vircom,strbod,dtx,dty,dtz)
        
c     torques in lab frame
        
        jr=0
        jg=0
        do ig=igrp1,igrp2
          
          jg=jg+1
          id=lstgtp(ig)
          
          tqx(jg)=0.d0
          tqy(jg)=0.d0
          tqz(jg)=0.d0
          do j=1,numgsit(id)
            
            jr=jr+1
            i=lstrgd(jr)
            
            tqx(jg)=tqx(jg)+dty(jr)*fzz(i)-dtz(jr)*fyy(i)
            tqy(jg)=tqy(jg)+dtz(jr)*fxx(i)-dtx(jr)*fzz(i)
            tqz(jg)=tqz(jg)+dtx(jr)*fyy(i)-dty(jr)*fxx(i)
            
          enddo
          
c     current rotational matrix 
          
          call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)
          
c     angular velocity at time step n (first guess)
          
          opx(jg)=omxo(jg)
          opy(jg)=omyo(jg)
          opz(jg)=omzo(jg)
          
c     iterate angular velocity for time step n (e. yezdimer)
          
          do i=1,5
            
            trx=(tqx(jg)*rot(1)+tqy(jg)*rot(4)+tqz(jg)*rot(7))*
     x        rotinx(id,2)
     x        +(rotiny(id,1)-rotinz(id,1))*opy(jg)*opz(jg)*
     x        rotinx(id,2)
            try=(tqx(jg)*rot(2)+tqy(jg)*rot(5)+tqz(jg)*rot(8))*
     x        rotiny(id,2)
     x        +(rotinz(id,1)-rotinx(id,1))*opz(jg)*opx(jg)*
     x        rotiny(id,2)
            trz=(tqx(jg)*rot(3)+tqy(jg)*rot(6)+tqz(jg)*rot(9))*
     x        rotinz(id,2)
     x        +(rotinx(id,1)-rotiny(id,1))*opx(jg)*opy(jg)*
     x        rotinz(id,2)
            
c     improved angular velocity at time step n
            
            opx(jg)=omxo(jg)+pt5*tstep*trx
            opy(jg)=omyo(jg)+pt5*tstep*try
            opz(jg)=omzo(jg)+pt5*tstep*trz
            
          enddo
          
c     correction due to thermostat
          
          delx=tstep*(trx-chit0*pt5*(omx(ig)+omxo(jg)))
          dely=tstep*(try-chit0*pt5*(omy(ig)+omyo(jg)))
          delz=tstep*(trz-chit0*pt5*(omz(ig)+omzo(jg)))
          
c     angular velocity at time step n
          
          omx(ig)=omxo(jg)+delx*pt5
          omy(ig)=omyo(jg)+dely*pt5
          omz(ig)=omzo(jg)+delz*pt5
          
c     angular velocity at time step n+1/2
          
          uxx(jg)=omxo(jg)+delx
          uyy(jg)=omyo(jg)+dely
          uzz(jg)=omzo(jg)+delz
          
c     angular velocity at time step n+1 (needed for quat algorithm)
          
          oqx(jg)=omxo(jg)+delx*1.5d0
          oqy(jg)=omyo(jg)+dely*1.5d0
          oqz(jg)=omzo(jg)+delz*1.5d0
          
        enddo
        
c     rotational kinetic energy
        
        engrot=getkinr(ngrp,idnode,mxnode)
        
c     restore half step angular velocities
        
        jg=0
        do ig=igrp1,igrp2
          
          jg=jg+1
          opx(jg)=omx(ig)
          opy(jg)=omy(ig)
          opz(jg)=omz(ig)
          omx(ig)=uxx(jg)
          omy(ig)=uyy(jg)
          omz(ig)=uzz(jg)
          
        enddo
        
c     assign new quaternions
        
        call update_quaternions
     x    (safeq,igrp1,igrp2,tstep,quattol,opx,opy,opz,oqx,oqy,oqz)
        
c     new atomic positions for atoms in rigid bodies
        
        jg=0
        jr=0
        do ig=igrp1,igrp2
          
          jg=jg+1
          id=lstgtp(ig)
          
c     new rotational matrix
          
          call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)
          
          do j=1,numgsit(id)
            
            jr=jr+1
            i=lstrgd(jr)
            
            xxx(i)=rot(1)*gxx(id,j)+rot(2)*gyy(id,j)+rot(3)*gzz(id,j)
     x        +gcmx(ig)
            yyy(i)=rot(4)*gxx(id,j)+rot(5)*gyy(id,j)+rot(6)*gzz(id,j)
     x        +gcmy(ig)
            zzz(i)=rot(7)*gxx(id,j)+rot(8)*gyy(id,j)+rot(9)*gzz(id,j)
     x        +gcmz(ig)
            
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
        
c     merge new atomic coordinates 
        
        if(mxnode.gt.1)call merge1
     x    (idnode,mxnode,natms,lstme,xxx,yyy,zzz,buffer)
        
c     new estimate of chit
        
        engke=engfke+engtrn
        engtot=engke+engrot
        chitp=2.d0*(engtot-sigma)/qmass
        chitnew=chit+tstep*chitp
        chit0=pt5*(chit+chitnew)
        
        if(ntcons.gt.0) then
          
c     apply constraint correction
          
          newstep=.false.
          if(icyc.eq.1)newstep=.true.
          
          call qshake
     x      (newstep,safe,lshmov,idnode,imcon,mxnode,natms,
     x      nscons,tolnce,tstep,viracc,dxx,dyy,dzz,dxt,dyt,
     x      dzt,txx,tyy,tzz,xxt,yyt,zzt,strcon)
          
          if(abs(viracc).le.1.d-10.and.icyc.gt.3)cycle=.false.
          
          vircon=vircon+viracc
          do i=1,9
            strcns(i)=strcns(i)+strcon(i)
          enddo
          
c     end of shake corrections
          
        endif
        
      enddo
      
c     minimum images of group positions and particle positions
      
      call images(imcon,0,1,natms,cell,xxx,yyy,zzz)
      call images(imcon,idnode,mxnode,ngrp,cell,gcmx,gcmy,gcmz)
      
c     new atomic positions for atoms in rigid bodies
      
      jr=0
      jg=0
      do ig=igrp1,igrp2
        
        jg=jg+1
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
        
c     merge new atomic velocities
        
        call merge1(idnode,mxnode,natms,lstme,vxx,vyy,vzz,buffer)
        
c     merge new quaternions
        
        call merge4(idnode,mxnode,ngrp,mxbuff,q0,q1,q2,q3,buffer)
        
      endif
      
c     ensure all atoms are within cell boundaries
      
      call images(imcon,0,1,natms,cell,xxx,yyy,zzz)
      
c     complete stress tensor
      
      do i=1,9
        stress(i)=stress(i)+strcns(i)+strkin(i)+strgrp(i)+strbod(i)
      enddo
      
c     update thermostat variable
      
      chit=chitnew
      
c     conserved quantity less kinetic and potential energy terms
      
      conint=conint+tstep*chit0*qmass/taut**2
      consv=conint+pt5*qmass*chit0**2
      
c     deallocate work arrays
      
      deallocate (opx,opy,opz,oqx,oqy,oqz,stat=fail(1))
      deallocate (dtx,dty,dtz,dxx,dyy,dzz,stat=fail(2))
      deallocate (uxx,uyy,uzz,txx,tyy,tzz,stat=fail(3))
      deallocate (dxt,dyt,dzt,xxo,yyo,zzo,stat=fail(4))
      deallocate (vxo,vyo,vzo,tqx,tqy,tqz,stat=fail(5))
      deallocate (fmx,fmy,fmz,omxo,omyo,omzo,stat=fail(6))
      deallocate (gvxo,gvyo,gvzo,gcxo,gcyo,gczo,stat=fail(7))
      deallocate (xxt,yyt,zzt,qn0,qn1,qn2,qn3,stat=fail(8))
      
      return
      end subroutine nvtq_h2
      
      subroutine nptq_b2
     x  (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x  ntcons,ntfree,ntpatm,elrc,engke,engrot,virlrc,press,
     x  quattol,sigma,taup,taut,tolnce,tstep,vircom,vircon,
     x  virtot,volm)
      
c***********************************************************************
c     
c     dlpoly routine to integrate rigid body equations of motion
c     using implicit leapfrog quaternion algorithm
c     plus bond constraints - rigid body sites and constraint sites
c     may coincide.
c     
c     verlet leapfrog with Berendsen thermostat and barostat.
c     (cell may change volume)
c
c     parallel replicated data version : block data
c     
c     tqx,tqy,tqz=torque in lab fixed frame (input)
c     omx,omy,omz=angular velocity in body fixed frame (principal axis)
c     rotinx,y,z =rotational inertia in body fixed frame
c     
c     copyright daresbury laboratory 1995
c     author      t.forester june 1995
c     amended     w.smith     sep 1999
c     amended     w.smith     nov 2005
c     
c**********************************************************************
      
      implicit none
      
      integer, parameter :: nnn=16
      real(8), parameter :: pt5=0.5d0
      
      logical safe,safeq,lshmov,newstep,newjob,cycle
      integer fail,imcon,idnode,mxnode,natms,ngrp,nscons
      integer ntcons,ntfree,igrp,igrp1,igrp2,idum,ifre,ifre1,ifre2
      integer i,j,k,jg,ig,jr,id,mxshak1,icyc,ntpatm
      real(8) engke,engrot,quattol,tolnce,tstep,vircom,vircon
      real(8) rot,strkin,strgrp,strcon,engtrn,engfke,trx,try,trz
      real(8) delx,dely,delz,czero
      real(8) vaa,vbb,vcc,viracc,beta,elrc,virlrc,press,engtot
      real(8) sigma,taup,taut,virtot,volm,cell0,elrc0,virlrc0
      real(8) chit0,volm0,chip0,psyst,scale
      
      real(8), allocatable :: dens0(:)
      real(8), allocatable :: opx(:),opy(:),opz(:)
      real(8), allocatable :: oqx(:),oqy(:),oqz(:)
      real(8), allocatable :: dtx(:),dty(:),dtz(:)
      real(8), allocatable :: dxx(:),dyy(:),dzz(:)
      real(8), allocatable :: uxx(:),uyy(:),uzz(:)
      real(8), allocatable :: txx(:),tyy(:),tzz(:)
      real(8), allocatable :: xxt(:),yyt(:),zzt(:)
      real(8), allocatable :: dxt(:),dyt(:),dzt(:)
      real(8), allocatable :: xxo(:),yyo(:),zzo(:)
      real(8), allocatable :: vxo(:),vyo(:),vzo(:)
      real(8), allocatable :: fmx(:),fmy(:),fmz(:)
      real(8), allocatable :: tqx(:),tqy(:),tqz(:)
      real(8), allocatable :: omxo(:),omyo(:),omzo(:)
      real(8), allocatable :: gvxo(:),gvyo(:),gvzo(:)
      real(8), allocatable :: gcxo(:),gcyo(:),gczo(:)
      real(8), allocatable :: qn0(:),qn1(:),qn2(:),qn3(:)
      
      dimension rot(9),strkin(9),strcon(9),strgrp(9),fail(nnn)
      dimension cell0(9),czero(9)
      
      save newjob,volm0,elrc0,virlrc0,czero,chit0,chip0,dens0
      save igrp1,igrp2,ifre1,ifre2
      
      data newjob/.true./
      data beta/7.3728d-3/
      
c     allocate working arrays
      
      do i=1,nnn
        fail(i)=0
      enddo
      allocate (opx(msgrp),opy(msgrp),opz(msgrp),stat=fail(1))
      allocate (oqx(msgrp),oqy(msgrp),oqz(msgrp),stat=fail(2))
      allocate (dtx(mxatms),dty(mxatms),dtz(mxatms),stat=fail(3))
      allocate (dxx(mxcons),dyy(mxcons),dzz(mxcons),stat=fail(4))
      allocate (uxx(mxatms),uyy(mxatms),uzz(mxatms),stat=fail(5))
      allocate (txx(mxatms),tyy(mxatms),tzz(mxatms),stat=fail(6))
      allocate (xxt(mxatms),yyt(mxatms),zzt(mxatms),stat=fail(7))
      allocate (dxt(mxcons),dyt(mxcons),dzt(mxcons),stat=fail(8))
      allocate (xxo(mxatms),yyo(mxatms),zzo(mxatms),stat=fail(9))
      allocate (vxo(msatms),vyo(msatms),vzo(msatms),stat=fail(10))
      allocate (fmx(msgrp),fmy(msgrp),fmz(msgrp),stat=fail(11))
      allocate (tqx(msgrp),tqy(msgrp),tqz(msgrp),stat=fail(12))
      allocate (omxo(msgrp),omyo(msgrp),omzo(msgrp),stat=fail(13))
      allocate (gvxo(msgrp),gvyo(msgrp),gvzo(msgrp),stat=fail(14))
      allocate (gcxo(msgrp),gcyo(msgrp),gczo(msgrp),stat=fail(15))
      allocate (qn0(msgrp),qn1(msgrp),qn2(msgrp),qn3(msgrp),
     x  stat=fail(16))
      do i=1,nnn
        if(fail(i).ne.0)call error(idnode,1670)
      enddo
      
c     store initial values of volume and long range corrections
      
      if(newjob) then
        
        volm0=volm
        elrc0=elrc
        virlrc0=virlrc
        fail(1)=0
        allocate (dens0(mxatyp),stat=fail(1))
        if(fail(1).ne.0)call error(idnode,1660)
        do i=1,ntpatm
          dens0(i)=dens(i)
        enddo
        do i=1,9
          czero(i)=cell(i)
        enddo
        newjob=.false.
        
        chit0=1.d0
        chip0=1.d0
        
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
        
      endif
      
      safe=.false.
      cycle=.true.
      
c     set constraint stress and virial and inital cell vectors
      
      vircon=0.d0
      do i=1,9
        
        strcns(i)=0.d0
        cell0(i)=cell(i)
        
      enddo
      
c     store initial values of position and velocity
      
      do i=1,natms
        
        xxo(i)=xxx(i)
        yyo(i)=yyy(i)
        zzo(i)=zzz(i)
        
      enddo
      
      j=0
      do ifre=ifre1,ifre2
        
        j=j+1
        i=lstfre(ifre)
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
        gvxo(jg)=gvxx(ig)
        gvyo(jg)=gvyy(ig)
        gvzo(jg)=gvzz(ig)
        omxo(jg)=omx(ig)
        omyo(jg)=omy(ig)
        omzo(jg)=omz(ig)
        qn0(jg)=q0(ig)
        qn1(jg)=q1(ig)
        qn2(jg)=q2(ig)
        qn3(jg)=q3(ig)
        
      enddo
      
c     construct current bond vectors
      
      if(ntcons.gt.0)then
        
        do k=1,nscons
          
          i=listcon(k,2)
          j=listcon(k,3)
          
          dxx(k)=xxx(i)-xxx(j)
          dyy(k)=yyy(i)-yyy(j)
          dzz(k)=zzz(i)-zzz(j)
          
        enddo
        
c     periodic boundary condition for bond vectors
        
        call images(imcon,0,1,nscons,cell,dxx,dyy,dzz)
        
      endif
      
c     calculate atom displacement from rigid body com
      
      jg=0
      jr=0
      do ig=igrp1,igrp2
        
        jg=jg+1
        id=lstgtp(ig)
        
        do j=1,numgsit(id)
          
          jr=jr+1
          i=lstrgd(jr)
          
          dtx(jr)=xxo(i)-gcxo(jg)
          dty(jr)=yyo(i)-gcyo(jg)
          dtz(jr)=zzo(i)-gczo(jg)
          
        enddo
        
      enddo
      
c     minimum images
      
      call images(imcon,0,1,jr,cell,dtx,dty,dtz)
      
c     length scaling parameter
      
      scale=chip0**(1.d0/3.d0)
      
c     shake iterations and thermostat iterations start here
      
      icyc=0
      mxshak1=mxshak
      if(ntcons.eq.0)mxshak1=3
      do while(cycle.and.icyc.le.mxshak1)
        
        icyc=icyc+1
        
c     restore cell vectors
        
        do i=1,9
          cell(i)=cell0(i)
        enddo

c     restore original quaternions for this step

        jg=0
        do ig=igrp1,igrp2

          jg=jg+1
          q0(ig)=qn0(jg)
          q1(ig)=qn1(jg)
          q2(ig)=qn2(jg)
          q3(ig)=qn3(jg)
          
        enddo

c     integrate 'free' particles
        
        j=0
        do ifre=ifre1,ifre2
          
          j=j+1
          i=lstfre(ifre)
          
c     advance velocity by leapfrog
          
          uxx(j)=(vxo(j)+tstep*rmass(i)*fxx(i))*chit0
          uyy(j)=(vyo(j)+tstep*rmass(i)*fyy(i))*chit0
          uzz(j)=(vzo(j)+tstep*rmass(i)*fzz(i))*chit0
          
c     advance position by leapfrog
          
          xxx(i)=xxo(i)*scale+tstep*uxx(j)
          yyy(i)=yyo(i)*scale+tstep*uyy(j)
          zzz(i)=zzo(i)*scale+tstep*uzz(j)
          
c     estimate full step velocity
          
          vxx(i)=pt5*(vxo(j)+uxx(j))
          vyy(i)=pt5*(vyo(j)+uyy(j))
          vzz(i)=pt5*(vzo(j)+uzz(j))
          
        enddo
        
c     calculate kinetic energy
        
        engfke=getkinf(ntfree,idnode,mxnode)
        
c     calculate kinetic stress tensor
        
        call kinstressf(ntfree,idnode,mxnode,strkin)
        
c     restore half step velocity
        
        j=0
        do ifre=ifre1,ifre2
          
          j=j+1
          i=lstfre(ifre)
          vxx(i)=uxx(j)
          vyy(i)=uyy(j)
          vzz(i)=uzz(j)
          
        enddo
        
c     *************  Rigid body motion ****************************
        
c     translational kinetic energy
        
        jg=0
        jr=0
        do ig=igrp1,igrp2
          
          jg=jg+1
          id=lstgtp(ig)
          
c     calculate centre of mass forces
          
          fmx(jg)=0.d0
          fmy(jg)=0.d0
          fmz(jg)=0.d0
          do j=1,numgsit(id)
            
            jr=jr+1
            i=lstrgd(jr)
            
            fmx(jg)=fmx(jg)+fxx(i)
            fmy(jg)=fmy(jg)+fyy(i)
            fmz(jg)=fmz(jg)+fzz(i)
            
          enddo
          
c     advance velocity by leapfrog
          
          uxx(jg)=(gvxo(jg)+tstep/gmass(id)*fmx(jg))*chit0
          uyy(jg)=(gvyo(jg)+tstep/gmass(id)*fmy(jg))*chit0
          uzz(jg)=(gvzo(jg)+tstep/gmass(id)*fmz(jg))*chit0
          
c     advance position by leapfrog
          
          gcmx(ig)=gcxo(jg)*scale+tstep*uxx(jg)
          gcmy(ig)=gcyo(jg)*scale+tstep*uyy(jg)
          gcmz(ig)=gczo(jg)*scale+tstep*uzz(jg)
          
c     centre of mass velocities at full step
          
          gvxx(ig)=pt5*(gvxo(jg)+uxx(jg))
          gvyy(ig)=pt5*(gvyo(jg)+uyy(jg))
          gvzz(ig)=pt5*(gvzo(jg)+uzz(jg))
          
        enddo
        
c     translational kinetic energy 
        
        engtrn=getkint(ngrp,idnode,mxnode)
        
c     kinetic contribution to stress tensor
        
        call kinstressg(ngrp,idnode,mxnode,strgrp)
        
c     restore half step velocity
        
        jg=0
        do ig=igrp1,igrp2
          
          jg=jg+1
          gvxx(ig)=uxx(jg)
          gvyy(ig)=uyy(jg)
          gvzz(ig)=uzz(jg)
          
        enddo
        
c     calculate rigid body stress tensor
        
        call bodystress(idnode,mxnode,ngrp,vircom,strbod,dtx,dty,dtz)
        
c     torques in lab frame
        
        jr=0
        jg=0
        do ig=igrp1,igrp2
          
          jg=jg+1
          id=lstgtp(ig)
          
          tqx(jg)=0.d0
          tqy(jg)=0.d0
          tqz(jg)=0.d0
          do j=1,numgsit(id)
            
            jr=jr+1
            i=lstrgd(jr)
            
            tqx(jg)=tqx(jg)+dty(jr)*fzz(i)-dtz(jr)*fyy(i)
            tqy(jg)=tqy(jg)+dtz(jr)*fxx(i)-dtx(jr)*fzz(i)
            tqz(jg)=tqz(jg)+dtx(jr)*fyy(i)-dty(jr)*fxx(i)
            
          enddo
          
c     current rotational matrix 
          
          call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)
          
c     angular velocity at time step n (first guess)
          
          opx(jg)=omxo(jg)
          opy(jg)=omyo(jg)
          opz(jg)=omzo(jg)
          
c     iterate angular velocity for time step n (e. yezdimer)
          
          do i=1,5
            
            trx=(tqx(jg)*rot(1)+tqy(jg)*rot(4)+tqz(jg)*rot(7))*
     x        rotinx(id,2)
     x        +(rotiny(id,1)-rotinz(id,1))*opy(jg)*opz(jg)*
     x        rotinx(id,2)
            try=(tqx(jg)*rot(2)+tqy(jg)*rot(5)+tqz(jg)*rot(8))*
     x        rotiny(id,2)
     x        +(rotinz(id,1)-rotinx(id,1))*opz(jg)*opx(jg)*
     x        rotiny(id,2)
            trz=(tqx(jg)*rot(3)+tqy(jg)*rot(6)+tqz(jg)*rot(9))*
     x        rotinz(id,2)
     x        +(rotinx(id,1)-rotiny(id,1))*opx(jg)*opy(jg)*
     x        rotinz(id,2)
            
            delx=tstep*trx
            dely=tstep*try
            delz=tstep*trz
            
c     improved angular velocity at time step n
            
            opx(jg)=(omxo(jg)+pt5*delx)
            opy(jg)=(omyo(jg)+pt5*dely)
            opz(jg)=(omzo(jg)+pt5*delz)
            
          enddo
          
c     scaled angular velocity at time step n
          
          omx(ig)=opx(jg)*chit0
          omy(ig)=opy(jg)*chit0
          omz(ig)=opz(jg)*chit0
          
c     angular velocity at time step n+1  (needed for quat algorithm)
          
          oqx(jg)=(omxo(jg)+delx*1.5d0)*chit0
          oqy(jg)=(omyo(jg)+dely*1.5d0)*chit0
          oqz(jg)=(omzo(jg)+delz*1.5d0)*chit0
          
c     angular velocity at time step n+1/2
          
          uxx(jg)=(omxo(jg)+delx)*chit0
          uyy(jg)=(omyo(jg)+dely)*chit0
          uzz(jg)=(omzo(jg)+delz)*chit0
          
        enddo
        
c     rotational kinetic energy
        
        engrot=getkinr(ngrp,idnode,mxnode)
        
c     restore half step velocities
        
        jg=0
        do ig=igrp1,igrp2
          
          jg=jg+1
          opx(jg)=omx(ig)
          opy(jg)=omy(ig)
          opz(jg)=omz(ig)
          omx(ig)=uxx(jg)
          omy(ig)=uyy(jg)
          omz(ig)=uzz(jg)
          
        enddo
        
c     assign new quaternions
        
        call update_quaternions
     x    (safeq,igrp1,igrp2,tstep,quattol,opx,opy,opz,oqx,oqy,oqz)
        
c     minimum images of group positions and particle positions
        
        call images(imcon,0,1,natms,cell,xxx,yyy,zzz)
        call images(imcon,idnode,mxnode,ngrp,cell,gcmx,gcmy,gcmz)
        
c     new atomic positions for atoms in rigid bodies
        
        jg=0
        jr=0
        do ig=igrp1,igrp2
          
          jg=jg+1
          id=lstgtp(ig)
          
c     new rotational matrix
          
          call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)
          
          do j=1,numgsit(id)
            
            jr=jr+1
            i=lstrgd(jr)
            
            xxx(i)=rot(1)*gxx(id,j)+rot(2)*gyy(id,j)+rot(3)*gzz(id,j)
     x        +gcmx(ig)
            yyy(i)=rot(4)*gxx(id,j)+rot(5)*gyy(id,j)+rot(6)*gzz(id,j)
     x        +gcmy(ig)
            zzz(i)=rot(7)*gxx(id,j)+rot(8)*gyy(id,j)+rot(9)*gzz(id,j)
     x        +gcmz(ig)
            
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
        
c     merge new atomic coordinates 
        
        if(mxnode.gt.1)call merge1
     x    (idnode,mxnode,natms,lstme,xxx,yyy,zzz,buffer)
        
c     pressure control variable
        
        engke=engfke+engtrn
        psyst=(2.d0*engke-virtot-vircon-vircom)/(3.d0*volm)
        chip0=1.d0+beta*tstep*(psyst-press)/taup
        scale=chip0**(1.d0/3.d0)
        
c     new estimate of chit
        
        engtot=engke+engrot 
        chit0=sqrt(1.d0+tstep/taut*(sigma/engtot-1.d0))
        
        if(ntcons.gt.0) then
          
c     apply constraint correction
          
          newstep=.false.
          if(icyc.eq.1)newstep=.true.
          
c     new cell parameters 
          
          do i=1,9
            cell(i)=scale*cell0(i)
          enddo
          
          call qshake
     x      (newstep,safe,lshmov,idnode,imcon,mxnode,natms,
     x      nscons,tolnce,tstep,viracc,dxx,dyy,dzz,dxt,dyt,
     x      dzt,txx,tyy,tzz,xxt,yyt,zzt,strcon)
          
          if(abs(viracc).le.1.d-10)cycle=.false.
          
          vircon=vircon+viracc
          do i=1,9
            strcns(i)=strcns(i)+strcon(i)
          enddo
          
c     end of shake corrections
          
        endif
        
      enddo
      
c     scale cell vectors
      
      scale=((chip0*volm)/volm0)**(1.d0/3.d0)
      do i=1,9
        cell(i)=scale*czero(i)
      enddo
      
c     construct scaling tensor (for later!)
      
      do i=2,8
        eta(i)=0.d0
      enddo
      eta(1)=scale
      eta(5)=scale
      eta(9)=scale
      
c     adjust long range corrections and number density
      
      elrc=elrc0*(volm0/volm)
      virlrc=virlrc0*(volm0/volm)
      
      do k=1,ntpatm
        dens(k)=dens0(k)*(volm0/volm)
      enddo
      
      if(mxnode.gt.1) then
        
c     merge new group coordinates and velocities
        
        call merge(idnode,mxnode,ngrp,mxbuff,gcmx,gcmy,gcmz,buffer)
        call merge(idnode,mxnode,ngrp,mxbuff,gvxx,gvyy,gvzz,buffer)
        
c     merge new atomic velocities
        
        call merge1(idnode,mxnode,natms,lstme,xxx,yyy,zzz,buffer)
        call merge1(idnode,mxnode,natms,lstme,vxx,vyy,vzz,buffer)
        
c     merge new quaternions
        
        call merge4(idnode,mxnode,ngrp,mxbuff,q0,q1,q2,q3,buffer)
        
      endif
      
c     ensure all atoms are within cell boundaries
      
      call images(imcon,0,1,natms,cell,xxx,yyy,zzz)
      
c     complete stress tensor
      
      do i=1,9
        stress(i)=stress(i)+strcns(i)+strkin(i)+strgrp(i)+strbod(i)
      enddo

c     deallocate work arrays
      
      deallocate (opx,opy,opz,oqx,oqy,oqz,stat=fail(1))
      deallocate (dtx,dty,dtz,dxx,dyy,dzz,stat=fail(2))
      deallocate (uxx,uyy,uzz,txx,tyy,tzz,stat=fail(3))
      deallocate (dxt,dyt,dzt,xxo,yyo,zzo,stat=fail(4))
      deallocate (vxo,vyo,vzo,tqx,tqy,tqz,stat=fail(5))
      deallocate (fmx,fmy,fmz,omxo,omyo,omzo,stat=fail(6))
      deallocate (gvxo,gvyo,gvzo,gcxo,gcyo,gczo,stat=fail(7))
      deallocate (xxt,yyt,zzt,qn0,qn1,qn2,qn3,stat=fail(8))
      
      return
      end subroutine nptq_b2
      
      subroutine nptq_h2
     x  (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x  ntcons,ntfree,ntpatm,chip,chit,consv,conint,elrc,engke,
     x  engrot,virlrc,press,quattol,sigma,taup,taut,temp,tolnce,
     x  tstep,vircom,vircon,virtot,volm)
      
c***********************************************************************
c     
c     dlpoly routine to integrate rigid body equations of motion
c     using implicit leapfrog quaternion algorithm
c     plus bond constraints - rigid body sites and constraint sites
c     may coincide.
c     
c     verlet leapfrog with Hoover thermostat and barostat. 
c     (cell may change volume)
c     
c     parallel replicated data version : block data
c     
c     tqx,tqy,tqz=torque in lab fixed frame (input)
c     omx,omy,omz=angular velocity in body fixed frame (principal axis)
c     rotinx,y,z =rotational inertia in body fixed frame
c     
c     copyright daresbury laboratory 1995
c     author      t.forester june 1995
c     amended     w.smith sep 1999 : euler equation
c     
c**********************************************************************
      
      implicit none
      
      integer, parameter :: nnn=16
      real(8), parameter :: pt5=0.5d0
      
      logical safe,safeq,lshmov,newstep,newjob,cycle
      integer fail,imcon,idnode,mxnode,natms,ngrp,nscons
      integer ntcons,ntfree,igrp,igrp1,igrp2,idum,ifre,ifre1,ifre2
      integer i,j,k,jg,ig,jr,id,mxshak1,icyc,ntpatm
      real(8) engke,engrot,quattol,tolnce,tstep,vircom,vircon,com
      real(8) rot,strkin,strcon,strgrp,vom,engtrn,trx,try,trz
      real(8) delx,dely,delz,engfke
      real(8) vaa,vbb,vcc,viracc,pmass,qmass,totmas,czero
      real(8) chip,chit,consv,conint,elrc,virlrc,press,sigma,taup,taut
      real(8) temp,virtot,cell0,volm0,elrc0,virlrc0
      real(8) chit0,chip0,chipnew,chipp,engtot,chitnew,chitp,volnew
      real(8) scale,volm,vold,cons1,cons2,cons3
      
      real(8), allocatable :: dens0(:)
      real(8), allocatable :: opx(:),opy(:),opz(:)
      real(8), allocatable :: oqx(:),oqy(:),oqz(:)
      real(8), allocatable :: dtx(:),dty(:),dtz(:)
      real(8), allocatable :: dxx(:),dyy(:),dzz(:)
      real(8), allocatable :: uxx(:),uyy(:),uzz(:)
      real(8), allocatable :: txx(:),tyy(:),tzz(:)
      real(8), allocatable :: xxt(:),yyt(:),zzt(:)
      real(8), allocatable :: dxt(:),dyt(:),dzt(:)
      real(8), allocatable :: xxo(:),yyo(:),zzo(:)
      real(8), allocatable :: vxo(:),vyo(:),vzo(:)
      real(8), allocatable :: fmx(:),fmy(:),fmz(:)
      real(8), allocatable :: tqx(:),tqy(:),tqz(:)
      real(8), allocatable :: omxo(:),omyo(:),omzo(:)
      real(8), allocatable :: gvxo(:),gvyo(:),gvzo(:)
      real(8), allocatable :: gcxo(:),gcyo(:),gczo(:)
      real(8), allocatable :: qn0(:),qn1(:),qn2(:),qn3(:)
      
      dimension rot(9),strkin(9),strcon(9),strgrp(9),fail(nnn)
      dimension czero(9),cell0(9),com(3),vom(3)
      
      save newjob,volm0,elrc0,virlrc0,czero,dens0,pmass,qmass
      save igrp1,igrp2,ifre1,ifre2,totmas
      
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
      allocate (txx(mxatms),tyy(mxatms),tzz(mxatms),stat=fail(6))
      allocate (xxt(mxatms),yyt(mxatms),zzt(mxatms),stat=fail(7))
      allocate (dxt(mxcons),dyt(mxcons),dzt(mxcons),stat=fail(8))
      allocate (xxo(mxatms),yyo(mxatms),zzo(mxatms),stat=fail(9))
      allocate (vxo(msatms),vyo(msatms),vzo(msatms),stat=fail(10))
      allocate (fmx(msgrp),fmy(msgrp),fmz(msgrp),stat=fail(11))
      allocate (tqx(msgrp),tqy(msgrp),tqz(msgrp),stat=fail(12))
      allocate (omxo(msgrp),omyo(msgrp),omzo(msgrp),stat=fail(13))
      allocate (gvxo(msgrp),gvyo(msgrp),gvzo(msgrp),stat=fail(14))
      allocate (gcxo(msgrp),gcyo(msgrp),gczo(msgrp),stat=fail(15))
      allocate (qn0(msgrp),qn1(msgrp),qn2(msgrp),qn3(msgrp),
     x  stat=fail(16))
      do i=1,nnn
        if(fail(i).ne.0)call error(idnode,1670)
      enddo
      
      if(newjob) then
        
c     inertia parameter for Nose-Hoover thermostat
        
        qmass=2.0d0*sigma*taut**2
        pmass=2.0d0*sigma*taup**2
        
c     store initial values of volume and long range corrections
        
        volm0=volm
        elrc0=elrc
        virlrc0=virlrc
        fail(1)=0
        allocate (dens0(mxatyp),stat=fail(1))
        if(fail(1).ne.0)call error(idnode,1680)
        do i=1,ntpatm
          dens0(i)=dens(i)
        enddo
        do i=1,9
          czero(i)=cell(i)
        enddo
        
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
        
c     calculate system mass
        
        totmas=getmass(natms,idnode,mxnode)
        
        newjob=.false.
        
      endif
      
      safe=.false.
      cycle=.true.
      
c     ensure total momentum is zero
      
      call getvom(natms,idnode,mxnode,totmas,vom)
      
      do ifre=ifre1,ifre2
        
        i=lstfre(ifre)
        vxx(i)=vxx(i)-vom(1)
        vyy(i)=vyy(i)-vom(2)
        vzz(i)=vzz(i)-vom(3)
        
      enddo
      
      do ig=igrp1,igrp2
        
        gvxx(ig)=gvxx(ig)-vom(1)
        gvyy(ig)=gvyy(ig)-vom(2)
        gvzz(ig)=gvzz(ig)-vom(3)
        
      enddo
      
c     store initial values of position and velocity
      
      do i=1,natms
        
        xxo(i)=xxx(i)
        yyo(i)=yyy(i)
        zzo(i)=zzz(i)
        
      enddo
      
      j=0
      do ifre=ifre1,ifre2
        
        j=j+1
        i=lstfre(ifre)
        vxo(j)=vxx(i)
        vyo(j)=vyy(i)
        vzo(j)=vzz(i)
        
      enddo
      
      jg=0
      do ig=igrp1,igrp2
        
        jg=jg+1
        omxo(jg)=omx(ig)
        omyo(jg)=omy(ig)
        omzo(jg)=omz(ig)
        gcxo(jg)=gcmx(ig)
        gcyo(jg)=gcmy(ig)
        gczo(jg)=gcmz(ig)
        gvxo(jg)=gvxx(ig)
        gvyo(jg)=gvyy(ig)
        gvzo(jg)=gvzz(ig)
        qn0(jg)=q0(ig)
        qn1(jg)=q1(ig)
        qn2(jg)=q2(ig)
        qn3(jg)=q3(ig)
        
      enddo
      
c     calculate centre of mass
      
      call getcom(natms,idnode,mxnode,totmas,com)
      
c     construct current bond vectors
      
      if(ntcons.gt.0)then
        
        do k=1,nscons
          
          i=listcon(k,2)
          j=listcon(k,3)
          
          dxx(k)=xxx(i)-xxx(j)
          dyy(k)=yyy(i)-yyy(j)
          dzz(k)=zzz(i)-zzz(j)
          
        enddo
        
c     periodic boundary condition for bond vectors
        
        call images(imcon,0,1,nscons,cell,dxx,dyy,dzz)
        
      endif
      
c     calculate atom displacement from rigid body com
      
      jg=0
      jr=0
      do ig=igrp1,igrp2
        
        jg=jg+1
        id=lstgtp(ig)
        
        do j=1,numgsit(id)
          
          jr=jr+1
          i=lstrgd(jr)
          
          dtx(jr)=xxo(i)-gcxo(jg)
          dty(jr)=yyo(i)-gcyo(jg)
          dtz(jr)=zzo(i)-gczo(jg)
          
        enddo
        
      enddo
      
c     minimum images
      
      call images(imcon,0,1,jr,cell,dtx,dty,dtz)
      
c     initial thermostat and barostat variables
      
      chit0=chit
      chip0=chip
      chipnew=chip
      
c     initialise constraint stress and virial
      
      vircon=0.d0
      do i=1,9
        
        strcns(i)=0.d0
        cell0(i)=cell(i)
        
      enddo
      
c     shake iterations and thermostat iterations start here
      
      icyc=0
      mxshak1=mxshak
      if(ntcons.eq.0)mxshak1=4
      do while(cycle.and.icyc.le.mxshak1)
        
        icyc=icyc+1
        
c     restore cell vectors
        
        do i=1,9
          cell(i)=cell0(i)
        enddo
        
c     restore original quaternions for this step

        jg=0
        do ig=igrp1,igrp2

          jg=jg+1
          q0(ig)=qn0(jg)
          q1(ig)=qn1(jg)
          q2(ig)=qn2(jg)
          q3(ig)=qn3(jg)
          
        enddo

c     integrate unconstrained new positions
        
        j=0
        do ifre=ifre1,ifre2
          
          j=j+1
          i=lstfre(ifre)
          
c     advance velocity using leapfrog
          
          uxx(j)=vxo(j)+tstep*(fxx(i)*rmass(i)-(chit0+chip0)*
     x      pt5*(vxx(i)+vxo(j)))
          uyy(j)=vyo(j)+tstep*(fyy(i)*rmass(i)-(chit0+chip0)*
     x      pt5*(vyy(i)+vyo(j)))
          uzz(j)=vzo(j)+tstep*(fzz(i)*rmass(i)-(chit0+chip0)*
     x      pt5*(vzz(i)+vzo(j)))
          
c     advance position using leapfrog
          
          xxx(i)=xxo(i)+tstep*(uxx(j)+
     x      chipnew*((xxx(i)+xxo(i))*pt5-com(1)))
          yyy(i)=yyo(i)+tstep*(uyy(j)+
     x      chipnew*((yyy(i)+yyo(i))*pt5-com(2)))
          zzz(i)=zzo(i)+tstep*(uzz(j)+
     x      chipnew*((zzz(i)+zzo(i))*pt5-com(3)))
          
c     estimate full step velocity
          
          vxx(i)=pt5*(uxx(j)+vxo(j))
          vyy(i)=pt5*(uyy(j)+vyo(j))
          vzz(i)=pt5*(uzz(j)+vzo(j))
          
        enddo
        
c     calculate new kinetic energy at current timestep
        
        engfke=getkinf(ntfree,idnode,mxnode)
        
c     calculate kinetic stress tensor
        
        call kinstressf(ntfree,idnode,mxnode,strkin)
        
c     restore half step velocity
        
        j=0
        do ifre=ifre1,ifre2
          
          j=j+1
          i=lstfre(ifre)
          vxx(i)=uxx(j)
          vyy(i)=uyy(j)
          vzz(i)=uzz(j)
          
        enddo
        
c     *************  Rigid body motion ****************************
        
c     translational kinetic energy
        
        jg=0
        jr=0
        do ig=igrp1,igrp2
          
          jg=jg+1
          id=lstgtp(ig)
          
c     calculate centre of mass forces
          
          fmx(jg)=0.d0
          fmy(jg)=0.d0
          fmz(jg)=0.d0
          do j=1,numgsit(id)
            
            jr=jr+1
            i=lstrgd(jr)
            
            fmx(jg)=fmx(jg)+fxx(i)
            fmy(jg)=fmy(jg)+fyy(i)
            fmz(jg)=fmz(jg)+fzz(i)
            
          enddo
          
c     advance velocity by leapfrog
          
          uxx(jg)=gvxo(jg)+tstep*(fmx(jg)/gmass(id)-(chit0+chip0)*
     x      pt5*(gvxx(ig)+gvxo(jg)))
          uyy(jg)=gvyo(jg)+tstep*(fmy(jg)/gmass(id)-(chit0+chip0)*
     x      pt5*(gvyy(ig)+gvyo(jg)))
          uzz(jg)=gvzo(jg)+tstep*(fmz(jg)/gmass(id)-(chit0+chip0)*
     x      pt5*(gvzz(ig)+gvzo(jg)))
          
c     advance position by leapfrog
          
          gcmx(ig)=gcxo(jg)+tstep*(uxx(jg)+
     x      chipnew*((gcxo(jg)+gcmx(ig))*pt5-com(1)))
          gcmy(ig)=gcyo(jg)+tstep*(uyy(jg)+
     x      chipnew*((gcyo(jg)+gcmy(ig))*pt5-com(2)))
          gcmz(ig)=gczo(jg)+tstep*(uzz(jg)+
     x      chipnew*((gczo(jg)+gcmz(ig))*pt5-com(3)))
          
c     estimate full step velocities
          
          gvxx(ig)=pt5*(gvxo(jg)+uxx(jg))
          gvyy(ig)=pt5*(gvyo(jg)+uyy(jg))
          gvzz(ig)=pt5*(gvzo(jg)+uzz(jg))
          
        enddo
        
c     translational kinetic energy 
        
        engtrn=getkint(ngrp,idnode,mxnode)
        
c     kinetic contribution to stress tensor
        
        call kinstressg(ngrp,idnode,mxnode,strgrp)
        
c     restore half step velocities
        
        jg=0
        do ig=igrp1,igrp2
          
          jg=jg+1
          gvxx(ig)=uxx(jg)
          gvyy(ig)=uyy(jg)
          gvzz(ig)=uzz(jg)
          
        enddo
        
c     calculate rigid body stress tensor
        
        call bodystress(idnode,mxnode,ngrp,vircom,strbod,dtx,dty,dtz)
        
c     torques in lab frame
        
        jr=0
        jg=0
        do ig=igrp1,igrp2
          
          jg=jg+1
          id=lstgtp(ig)
          
          tqx(jg)=0.d0
          tqy(jg)=0.d0
          tqz(jg)=0.d0
          do j=1,numgsit(id)
            
            jr=jr+1
            i=lstrgd(jr)
            
            tqx(jg)=tqx(jg)+dty(jr)*fzz(i)-dtz(jr)*fyy(i)
            tqy(jg)=tqy(jg)+dtz(jr)*fxx(i)-dtx(jr)*fzz(i)
            tqz(jg)=tqz(jg)+dtx(jr)*fyy(i)-dty(jr)*fxx(i)
            
          enddo
          
c     current rotational matrix 
          
          call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)
          
c     angular velocity at time step n (first guess)
          
          opx(jg)=omxo(jg)
          opy(jg)=omyo(jg)
          opz(jg)=omzo(jg)
          
c     iterate angular velocity for time step n (e. yezdimer)
          
          do i=1,5
            
            trx=(tqx(jg)*rot(1)+tqy(jg)*rot(4)+tqz(jg)*rot(7))*
     x        rotinx(id,2)
     x        +(rotiny(id,1)-rotinz(id,1))*opy(jg)*opz(jg)*
     x        rotinx(id,2)
            try=(tqx(jg)*rot(2)+tqy(jg)*rot(5)+tqz(jg)*rot(8))*
     x        rotiny(id,2)
     x        +(rotinz(id,1)-rotinx(id,1))*opz(jg)*opx(jg)*
     x        rotiny(id,2)
            trz=(tqx(jg)*rot(3)+tqy(jg)*rot(6)+tqz(jg)*rot(9))*
     x        rotinz(id,2)
     x        +(rotinx(id,1)-rotiny(id,1))*opx(jg)*opy(jg)*
     x        rotinz(id,2)
            
c     improved angular velocity at time step n
            
            opx(jg)=omxo(jg)+pt5*tstep*trx
            opy(jg)=omyo(jg)+pt5*tstep*try
            opz(jg)=omzo(jg)+pt5*tstep*trz
            
          enddo
          
c     correction due to thermostat
          
          delx=tstep*(trx-chit0*pt5*(omx(ig)+omxo(jg)))
          dely=tstep*(try-chit0*pt5*(omy(ig)+omyo(jg)))
          delz=tstep*(trz-chit0*pt5*(omz(ig)+omzo(jg)))
          
c     angular velocity at time step n
          
          omx(ig)=omxo(jg)+delx*pt5
          omy(ig)=omyo(jg)+dely*pt5
          omz(ig)=omzo(jg)+delz*pt5
          
c     angular velocity at time step n+1/2
          
          uxx(jg)=omxo(jg)+delx
          uyy(jg)=omyo(jg)+dely
          uzz(jg)=omzo(jg)+delz
          
c     angular velocity at time step n+1  (needed for quat algorithm)
          
          oqx(jg)=omxo(jg)+delx*1.5d0
          oqy(jg)=omyo(jg)+dely*1.5d0
          oqz(jg)=omzo(jg)+delz*1.5d0
          
        enddo
        
c     rotational kinetic energy
        
        engrot=getkinr(ngrp,idnode,mxnode)
        
c     restore half step velocities
        
        jg=0
        do ig=igrp1,igrp2
          
          jg=jg+1
          opx(jg)=omx(ig)
          opy(jg)=omy(ig)
          opz(jg)=omz(ig)
          omx(ig)=uxx(jg)
          omy(ig)=uyy(jg)
          omz(ig)=uzz(jg)
          
        enddo
        
c     assign new quaternions
        
        call update_quaternions
     x    (safeq,igrp1,igrp2,tstep,quattol,opx,opy,opz,oqx,oqy,oqz)
        
c     minimum images of group positions and particle positions
        
        call images(imcon,0,1,natms,cell,xxx,yyy,zzz)
        call images(imcon,idnode,mxnode,ngrp,cell,gcmx,gcmy,gcmz)
        
c     new atomic positions for atoms in rigid bodies
        
        jg=0
        jr=0
        do ig=igrp1,igrp2
          
          jg=jg+1
          id=lstgtp(ig)
          
c     new rotational matrix
          
          call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)
          
          do j=1,numgsit(id)
            
            jr=jr+1
            i=lstrgd(jr)
            
            xxx(i)=rot(1)*gxx(id,j)+rot(2)*gyy(id,j)+rot(3)*gzz(id,j)
     x        +gcmx(ig)
            yyy(i)=rot(4)*gxx(id,j)+rot(5)*gyy(id,j)+rot(6)*gzz(id,j)
     x        +gcmy(ig)
            zzz(i)=rot(7)*gxx(id,j)+rot(8)*gyy(id,j)+rot(9)*gzz(id,j)
     x        +gcmz(ig)
            
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
        
c     merge new atomic coordinates 
        
        if(mxnode.gt.1)call merge1
     x    (idnode,mxnode,natms,lstme,xxx,yyy,zzz,buffer)
        
c     new estimate of chip and chit
        
        engke=engfke+engtrn
        chipp=(2.d0*engke-virtot-vircom-vircon-3.d0*press*volm)/pmass-
     x    chit0*chip0
        chipnew=chip+tstep*chipp
        chip0=pt5*(chip+chipnew)
        
        engtot=engke+engrot
        chitp=(2.d0*(engtot-sigma)+pmass*chip0**2-boltz*temp)/qmass
        chitnew=chit+tstep*chitp
        chit0=pt5*(chit+chitnew)
        
        if(ntcons.gt.0) then
          
c     apply constraint correction
          
          newstep=.false.
          if(icyc.eq.1)newstep=.true.
          
c     estimate new cell tensor
          
          volnew=volm*exp(3.d0*tstep*chipnew)
          scale=(volnew/volm0)**(1.d0/3.d0)
          do i=1,9
            cell(i)=czero(i)*scale
          enddo
          
          call qshake
     x      (newstep,safe,lshmov,idnode,imcon,mxnode,natms,
     x      nscons,tolnce,tstep,viracc,dxx,dyy,dzz,dxt,dyt,
     x      dzt,txx,tyy,tzz,xxt,yyt,zzt,strcon)
          
          if(abs(viracc).le.1.d-10)cycle=.false.
          
          vircon=vircon+viracc
          do i=1,9
            strcns(i)=strcns(i)+strcon(i)
          enddo
          
c     end of shake corrections
          
        endif
        
      enddo
      
c     update volume
      
      vold=volm
      volm=volm*exp(3.d0*tstep*chipnew)
      
c     scale cell vectors-isotropic
      
      scale=(volm/volm0)**(1.d0/3.d0)
      do i=1,9
        cell(i)=czero(i)*scale
      enddo
      
c     construct scaling tensor (for later!)
      
      do i=2,8
        eta(i)=0.d0
      enddo
      eta(1)=chipnew
      eta(5)=chipnew
      eta(9)=chipnew
      
c     adjust long range corrections and number density
      
      elrc=elrc0*(volm0/volm)
      virlrc=virlrc0*(volm0/volm)
      
      do k=1,ntpatm
        dens(k)=dens0(k)*(volm0/volm)
      enddo
      
c     update thermostat and barostat variables
      
      chit=chitnew
      chip=chipnew
      
c     conserved quantity less kinetic and potential energy terms
      
      conint=conint+tstep*chit0*(qmass/taut**2+boltz*temp)
      cons1=pt5*qmass*chit0**2
      cons2=press*vold
      cons3=pt5*pmass*chip0**2
      consv=conint+cons1+cons2+cons3
      
      if(mxnode.gt.1) then
        
c     merge new group coordinates and velocities
        
        call merge(idnode,mxnode,ngrp,mxbuff,gcmx,gcmy,gcmz,buffer)
        call merge(idnode,mxnode,ngrp,mxbuff,gvxx,gvyy,gvzz,buffer)
        
c     merge new atomic velocities
        
        call merge1(idnode,mxnode,natms,lstme,vxx,vyy,vzz,buffer)
        
c     merge new quaternions
        
        call merge4(idnode,mxnode,ngrp,mxbuff,q0,q1,q2,q3,buffer)
        
      endif
      
c     ensure all atoms are within cell boundaries
      
      call images(imcon,0,1,natms,cell,xxx,yyy,zzz)
      
c     complete stress tensor
      
      do i=1,9
        stress(i)=stress(i)+strcns(i)+strkin(i)+strgrp(i)+strbod(i)
      enddo
      
c     deallocate work arrays
      
      deallocate (opx,opy,opz,oqx,oqy,oqz,stat=fail(1))
      deallocate (dtx,dty,dtz,dxx,dyy,dzz,stat=fail(2))
      deallocate (uxx,uyy,uzz,txx,tyy,tzz,stat=fail(3))
      deallocate (dxt,dyt,dzt,xxo,yyo,zzo,stat=fail(4))
      deallocate (vxo,vyo,vzo,tqx,tqy,tqz,stat=fail(5))
      deallocate (fmx,fmy,fmz,omxo,omyo,omzo,stat=fail(6))
      deallocate (gvxo,gvyo,gvzo,gcxo,gcyo,gczo,stat=fail(7))
      deallocate (xxt,yyt,zzt,qn0,qn1,qn2,qn3,stat=fail(8))
      
      return
      end subroutine nptq_h2
      
      subroutine nstq_b2
     x  (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x  ntcons,ntfree,ntpatm,mode,elrc,engke,engrot,virlrc,press,
     x  quattol,sigma,taup,taut,tolnce,tstep,vircom,vircon,volm)
      
c***********************************************************************
c     
c     dlpoly routine to integrate rigid body equations of motion
c     using implicit leapfrog quaternion algorithm
c     plus bond constraints. Rigid body sites and constraint sites may
c     coincide. 
c     
c     verlet leapfrog with Berendsen thermostat and barostat. 
c     (cell may change shape)
c     
c     parallel replicated data version : block data
c     
c     for systems using bond constraints
c     
c     tqx,tqy,tqz=torque in lab fixed frame (input)
c     omx,omy,omz=angular velocity in body fixed frame (principal axis)
c     rotinx,y,z =rotational inertia in body fixed frame
c     
c     copyright daresbury laboratory 1995
c     author      t.forester june  1995
c     amended     w.smith nov 2005
c     
c**********************************************************************
      
      implicit none
      
      integer, parameter :: nnn=16
      real(8), parameter :: pt5=0.5d0
      
      logical safe,lshmov,newjob,safeq,newstep,cycle
      integer fail,imcon,idnode,mxnode,natms,ngrp,nscons,mode
      integer ntcons,ntfree,igrp,igrp1,igrp2,idum,ifre,ifre1,ifre2
      integer i,j,k,jg,ig,jr,id,mxshak1,icyc,ntpatm
      real(8) engke,engrot,quattol,tolnce,tstep,vircom,vircon
      real(8) rot,strkin,strcon,strgrp,engtrn,trx,try,trz
      real(8) delx,dely,delz,engfke
      real(8) vaa,vbb,vcc,viracc,elrc,virlrc,press,sigma
      real(8) taup,taut,volm,cell0,volm0,elrc0,chit0,uni
      real(8) beta,stres0,engtot,virlrc0
      
      real(8), allocatable :: dens0(:)
      real(8), allocatable :: opx(:),opy(:),opz(:)
      real(8), allocatable :: oqx(:),oqy(:),oqz(:)
      real(8), allocatable :: dtx(:),dty(:),dtz(:)
      real(8), allocatable :: dxx(:),dyy(:),dzz(:)
      real(8), allocatable :: uxx(:),uyy(:),uzz(:)
      real(8), allocatable :: txx(:),tyy(:),tzz(:)
      real(8), allocatable :: xxt(:),yyt(:),zzt(:)
      real(8), allocatable :: dxt(:),dyt(:),dzt(:)
      real(8), allocatable :: xxo(:),yyo(:),zzo(:)
      real(8), allocatable :: vxo(:),vyo(:),vzo(:)
      real(8), allocatable :: fmx(:),fmy(:),fmz(:)
      real(8), allocatable :: tqx(:),tqy(:),tqz(:)
      real(8), allocatable :: omxo(:),omyo(:),omzo(:)
      real(8), allocatable :: gvxo(:),gvyo(:),gvzo(:)
      real(8), allocatable :: gcxo(:),gcyo(:),gczo(:)
      real(8), allocatable :: qn0(:),qn1(:),qn2(:),qn3(:)
      
      dimension rot(9),strkin(9),strcon(9),strgrp(9),fail(nnn)
      dimension cell0(9),uni(9),stres0(9)
      
      save newjob,volm0,elrc0,virlrc0,chit0,dens0
      save igrp1,igrp2,ifre1,ifre2
      
      data uni/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/
      data newjob/.true./
      data beta/7.3728d-3/
      
c     allocate working arrays
      
      do i=1,nnn
        fail(i)=0
      enddo
      allocate (opx(msgrp),opy(msgrp),opz(msgrp),stat=fail(1))
      allocate (oqx(msgrp),oqy(msgrp),oqz(msgrp),stat=fail(2))
      allocate (dtx(mxatms),dty(mxatms),dtz(mxatms),stat=fail(3))
      allocate (dxx(mxcons),dyy(mxcons),dzz(mxcons),stat=fail(4))
      allocate (uxx(mxatms),uyy(mxatms),uzz(mxatms),stat=fail(5))
      allocate (txx(mxatms),tyy(mxatms),tzz(mxatms),stat=fail(6))
      allocate (xxt(mxatms),yyt(mxatms),zzt(mxatms),stat=fail(7))
      allocate (dxt(mxcons),dyt(mxcons),dzt(mxcons),stat=fail(8))
      allocate (xxo(mxatms),yyo(mxatms),zzo(mxatms),stat=fail(9))
      allocate (vxo(msatms),vyo(msatms),vzo(msatms),stat=fail(10))
      allocate (fmx(msgrp),fmy(msgrp),fmz(msgrp),stat=fail(11))
      allocate (tqx(msgrp),tqy(msgrp),tqz(msgrp),stat=fail(12))
      allocate (omxo(msgrp),omyo(msgrp),omzo(msgrp),stat=fail(13))
      allocate (gvxo(msgrp),gvyo(msgrp),gvzo(msgrp),stat=fail(14))
      allocate (gcxo(msgrp),gcyo(msgrp),gczo(msgrp),stat=fail(15))
      allocate (qn0(msgrp),qn1(msgrp),qn2(msgrp),qn3(msgrp),
     x  stat=fail(16))
      do i=1,nnn
        if(fail(i).ne.0)call error(idnode,1690)
      enddo
      
c     store initial values of volume, long range corrections etc
      
      if(newjob) then
        
        chit0=1.d0
        volm0=volm
        elrc0=elrc
        virlrc0=virlrc
        fail(1)=0
        allocate (dens0(mxatyp),stat=fail(1))
        if(fail(1).ne.0)call error(idnode,1700)
        do i=1,ntpatm
          dens0(i)=dens(i)
        enddo
        
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
      cycle=.true.
      
c     set virial, strain and stress tensors
      
      vircon=0.d0
      do i=1,9
        
        eta(i)=uni(i)
        strcns(i)=0.d0
        cell0(i)=cell(i)
        stres0(i)=stress(i)
        
      enddo
      
c     store initial values of position and velocity
      
      do i=1,natms
        
        xxo(i)=xxx(i)
        yyo(i)=yyy(i)
        zzo(i)=zzz(i)
        
      enddo
      
      j=0
      do ifre=ifre1,ifre2
        
        j=j+1
        i=lstfre(ifre)
        vxo(j)=vxx(i)
        vyo(j)=vyy(i)
        vzo(j)=vzz(i)
        
      enddo
      
      jg=0
      do ig=igrp1,igrp2
        
        jg=jg+1
        omxo(jg)=omx(ig)
        omyo(jg)=omy(ig)
        omzo(jg)=omz(ig)
        gcxo(jg)=gcmx(ig)
        gcyo(jg)=gcmy(ig)
        gczo(jg)=gcmz(ig)
        gvxo(jg)=gvxx(ig)
        gvyo(jg)=gvyy(ig)
        gvzo(jg)=gvzz(ig)
        qn0(jg)=q0(ig)
        qn1(jg)=q1(ig)
        qn2(jg)=q2(ig)
        qn3(jg)=q3(ig)
        
      enddo
      
c     construct current bond vectors
      
      if(ntcons.gt.0)then
        
        do k=1,nscons
          
          i=listcon(k,2)
          j=listcon(k,3)
          
          dxx(k)=xxx(i)-xxx(j)
          dyy(k)=yyy(i)-yyy(j)
          dzz(k)=zzz(i)-zzz(j)
          
        enddo
        
c     periodic boundary condition for bond vectors
        
        call images(imcon,0,1,nscons,cell,dxx,dyy,dzz)
        
      endif
      
c     calculate atom displacement from rigid body com
      
      jg=0
      jr=0
      do ig=igrp1,igrp2
        
        jg=jg+1
        id=lstgtp(ig)
        
        do j=1,numgsit(id)
          
          jr=jr+1
          i=lstrgd(jr)
          
          dtx(jr)=xxo(i)-gcxo(jg)
          dty(jr)=yyo(i)-gcyo(jg)
          dtz(jr)=zzo(i)-gczo(jg)
          
        enddo
        
      enddo
      
c     minimum images
      
      call images(imcon,0,1,jr,cell,dtx,dty,dtz)
      
c     shake iterations and thermostat iterations start here
      
      icyc=0
      mxshak1=mxshak
      if(ntcons.eq.0)mxshak1=4
      do while(cycle.and.icyc.le.mxshak1)
        
        icyc=icyc+1
        
c     restore original quaternions for this step

        jg=0
        do ig=igrp1,igrp2

          jg=jg+1
          q0(ig)=qn0(jg)
          q1(ig)=qn1(jg)
          q2(ig)=qn2(jg)
          q3(ig)=qn3(jg)
          
        enddo

c     unconstrained new positions
        
        j=0
        do ifre=ifre1,ifre2
          
          j=j+1
          i=lstfre(ifre)
          
c     advance velocity using leapfrog
          
          uxx(j)=(vxo(j)+tstep*rmass(i)*fxx(i))*chit0
          uyy(j)=(vyo(j)+tstep*rmass(i)*fyy(i))*chit0
          uzz(j)=(vzo(j)+tstep*rmass(i)*fzz(i))*chit0
          
c     update positions
          
          xxx(i)=tstep*uxx(j)+eta(1)*xxo(i)+eta(4)*yyo(i)+eta(7)*zzo(i)
          yyy(i)=tstep*uyy(j)+eta(2)*xxo(i)+eta(5)*yyo(i)+eta(8)*zzo(i)
          zzz(i)=tstep*uzz(j)+eta(3)*xxo(i)+eta(6)*yyo(i)+eta(9)*zzo(i)
          
c     calculate velocity at full time step
          
          vxx(i)=pt5*(uxx(j)+vxo(j))
          vyy(i)=pt5*(uyy(j)+vyo(j))
          vzz(i)=pt5*(uzz(j)+vzo(j))
          
        enddo
        
c     calculate kinetic energy
        
        engfke=getkinf(ntfree,idnode,mxnode)
        
c     kinetic stress tensor
        
        call kinstressf(ntfree,idnode,mxnode,strkin)
        
c     restore half step velocities
        
        j=0
        do ifre=ifre1,ifre2
          
          j=j+1
          i=lstfre(ifre)
          vxx(i)=uxx(j)
          vyy(i)=uyy(j)
          vzz(i)=uzz(j)
          
        enddo
        
c     ********: rigid body motion - thermostated  :***********
        
c     translational kinetic energy
        
        jg=0
        jr=0
        do ig=igrp1,igrp2
          
          jg=jg+1
          id=lstgtp(ig)
          
c     calculate centre of mass forces
          
          fmx(jg)=0.d0
          fmy(jg)=0.d0
          fmz(jg)=0.d0
          do j=1,numgsit(id)
            
            jr=jr+1
            i=lstrgd(jr)
            
            fmx(jg)=fmx(jg)+fxx(i)
            fmy(jg)=fmy(jg)+fyy(i)
            fmz(jg)=fmz(jg)+fzz(i)
            
          enddo
          
c     calculate thermostated velocities
          
          uxx(jg)=(gvxo(jg)+tstep/gmass(id)*fmx(jg))*chit0
          uyy(jg)=(gvyo(jg)+tstep/gmass(id)*fmy(jg))*chit0
          uzz(jg)=(gvzo(jg)+tstep/gmass(id)*fmz(jg))*chit0
          
c     update positions
          
          gcmx(ig)=tstep*uxx(jg)+
     x      eta(1)*gcxo(jg)+eta(4)*gcyo(jg)+eta(7)*gczo(jg)
          gcmy(ig)=tstep*uyy(jg)+
     x      eta(2)*gcxo(jg)+eta(5)*gcyo(jg)+eta(8)*gczo(jg)
          gcmz(ig)=tstep*uzz(jg)+
     x      eta(3)*gcxo(jg)+eta(6)*gcyo(jg)+eta(9)*gczo(jg)
          
c     centre of mass velocities at full step
          
          gvxx(ig)=pt5*(gvxo(jg)+uxx(jg))
          gvyy(ig)=pt5*(gvyo(jg)+uyy(jg))
          gvzz(ig)=pt5*(gvzo(jg)+uzz(jg))
          
        enddo
        
c     calculate kinetic energy
        
        engtrn=getkint(ngrp,idnode,mxnode)
        
c     kinetic stress tensor
        
        call kinstressg(ngrp,idnode,mxnode,strgrp)
        
c     restore half step velocity
        
        jg=0
        do ig=igrp1,igrp2
          
          jg=jg+1
          gvxx(ig)=uxx(jg)
          gvyy(ig)=uyy(jg)
          gvzz(ig)=uzz(jg)
          
        enddo
        
c     calculate rigid body stress tensor
        
        call bodystress(idnode,mxnode,ngrp,vircom,strbod,dtx,dty,dtz)
        
c     torques in lab frame
        
        jr=0
        jg=0
        do ig=igrp1,igrp2
          
          jg=jg+1
          id=lstgtp(ig)
          
          tqx(jg)=0.d0
          tqy(jg)=0.d0
          tqz(jg)=0.d0
          do j=1,numgsit(id)
            
            jr=jr+1
            i=lstrgd(jr)
            
            tqx(jg)=tqx(jg)+dty(jr)*fzz(i)-dtz(jr)*fyy(i)
            tqy(jg)=tqy(jg)+dtz(jr)*fxx(i)-dtx(jr)*fzz(i)
            tqz(jg)=tqz(jg)+dtx(jr)*fyy(i)-dty(jr)*fxx(i)
            
          enddo
          
c     current rotational matrix 
          
          call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)
          
c     angular velocity at time step n (first guess)
          
          opx(jg)=omxo(jg)
          opy(jg)=omyo(jg)
          opz(jg)=omzo(jg)
          
c     iterate angular velocity for time step n (e. yezdimer)
          
          do i=1,5
            
            trx=(tqx(jg)*rot(1)+tqy(jg)*rot(4)+tqz(jg)*rot(7))*
     x        rotinx(id,2)
     x        +(rotiny(id,1)-rotinz(id,1))*opy(jg)*opz(jg)*
     x        rotinx(id,2)
            try=(tqx(jg)*rot(2)+tqy(jg)*rot(5)+tqz(jg)*rot(8))*
     x        rotiny(id,2)
     x        +(rotinz(id,1)-rotinx(id,1))*opz(jg)*opx(jg)*
     x        rotiny(id,2)
            trz=(tqx(jg)*rot(3)+tqy(jg)*rot(6)+tqz(jg)*rot(9))*
     x        rotinz(id,2)
     x        +(rotinx(id,1)-rotiny(id,1))*opx(jg)*opy(jg)*
     x        rotinz(id,2)
            
            delx=tstep*trx
            dely=tstep*try
            delz=tstep*trz
            
c     improved angular velocity at time step n
            
            opx(jg)=(omxo(jg)+pt5*delx)
            opy(jg)=(omyo(jg)+pt5*dely)
            opz(jg)=(omzo(jg)+pt5*delz)
            
          enddo
          
c     scaled angular velocity at time step n
          
          omx(ig)=opx(jg)*chit0
          omy(ig)=opy(jg)*chit0
          omz(ig)=opz(jg)*chit0
          
c     angular velocity at time step n+1 (needed for quat algorithm)
          
          oqx(jg)=(omxo(jg)+delx*1.5d0)*chit0
          oqy(jg)=(omyo(jg)+dely*1.5d0)*chit0
          oqz(jg)=(omzo(jg)+delz*1.5d0)*chit0
          
c     angular velocity at time step n+1/2
          
          uxx(jg)=(omxo(jg)+delx)*chit0
          uyy(jg)=(omyo(jg)+dely)*chit0
          uzz(jg)=(omzo(jg)+delz)*chit0
          
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
          omx(ig)=uxx(jg)
          omy(ig)=uyy(jg)
          omz(ig)=uzz(jg)
          
        enddo
        
c     assign new quaternions
        
        call update_quaternions
     x    (safeq,igrp1,igrp2,tstep,quattol,opx,opy,opz,oqx,oqy,oqz)
        
c     minimum images of group positions and particle positions
        
        call images(imcon,0,1,natms,cell,xxx,yyy,zzz)
        call images(imcon,idnode,mxnode,ngrp,cell,gcmx,gcmy,gcmz)
        
c     new atomic positions for atoms in rigid bodies
        
        jg=0
        jr=0
        do ig=igrp1,igrp2
          
          jg=jg+1
          id=lstgtp(ig)
          
c     new rotational matrix
          
          call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)
          
          do j=1,numgsit(id)
            
            jr=jr+1
            i=lstrgd(jr)
            
            xxx(i)=rot(1)*gxx(id,j)+rot(2)*gyy(id,j)+rot(3)*gzz(id,j)
     x        +gcmx(ig)
            yyy(i)=rot(4)*gxx(id,j)+rot(5)*gyy(id,j)+rot(6)*gzz(id,j)
     x        +gcmy(ig)
            zzz(i)=rot(7)*gxx(id,j)+rot(8)*gyy(id,j)+rot(9)*gzz(id,j)
     x        +gcmz(ig)
            
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
        
c     merge new atomic coordinates 
        
        if(mxnode.gt.1)call merge1
     x    (idnode,mxnode,natms,lstme,xxx,yyy,zzz,buffer)
        
c     calculate total stress tensor
        
        do i=1,9
          stress(i)=stres0(i)+strcns(i)+strkin(i)+strgrp(i)+strbod(i)
        enddo
        
c     calculate new cell tensor
        
        call mat_mul(eta,cell0,cell)
        
c     calculate eta tensor
        
        do i=1,9
          eta(i)=beta*tstep/taup*(stress(i)/volm-press*uni(i))+uni(i)
        enddo
        if(mode.gt.0)then
          eta(3)=0.d0
          eta(6)=0.d0
          eta(7)=0.d0
          eta(8)=0.d0
          if(mode.lt.3)then
            eta(2)=0.d0
            eta(4)=0.d0
            if(mode.eq.2)then
              eta(1)=0.5d0*(eta(1)+eta(5))
              eta(5)=eta(1)
            endif
          endif
        endif
        
c     new estimate of chit
        
        engke=engfke+engtrn
        engtot=engke+engrot
        chit0=sqrt(1.d0+tstep/taut*(sigma/engtot-1.d0))
        
        if(ntcons.gt.0) then
          
c     apply constraint correction
          
          newstep=.false.
          if(icyc.eq.1)newstep=.true.
          
          call qshake
     x      (newstep,safe,lshmov,idnode,imcon,mxnode,natms,
     x      nscons,tolnce,tstep,viracc,dxx,dyy,dzz,dxt,dyt,
     x      dzt,txx,tyy,tzz,xxt,yyt,zzt,strcon)
          
          if(abs(viracc).le.1.d-10)cycle=.false.
          
          vircon=vircon+viracc
          do i=1,9
            strcns(i)=strcns(i)+strcon(i)
          enddo
          
c     end of shake corrections
          
        endif
        
      enddo
      
c     update volume
      
      volm=volm*eta(1)*eta(5)*eta(9)
      
c     adjust long range corrections and number density
      
      elrc=elrc0*(volm0/volm)
      virlrc=virlrc0*(volm0/volm)
      
      do k=1,ntpatm
        dens(k)=dens0(k)*(volm0/volm)
      enddo
      
      if(mxnode.gt.1) then
        
c     merge new group coordinates and velocities
        
        call merge(idnode,mxnode,ngrp,mxbuff,gcmx,gcmy,gcmz,buffer)
        call merge(idnode,mxnode,ngrp,mxbuff,gvxx,gvyy,gvzz,buffer)
        
c     merge new atomic velocities
        
        call merge1(idnode,mxnode,natms,lstme,xxx,yyy,zzz,buffer)
        call merge1(idnode,mxnode,natms,lstme,vxx,vyy,vzz,buffer)
        
c     merge new quaternions
        
        call merge4(idnode,mxnode,ngrp,mxbuff,q0,q1,q2,q3,buffer)
        
      endif
      
c     ensure all atoms are within cell boundaries
      
      call images(imcon,0,1,natms,cell,xxx,yyy,zzz)
      
c     deallocate work arrays
      
      deallocate (opx,opy,opz,oqx,oqy,oqz,stat=fail(1))
      deallocate (dtx,dty,dtz,dxx,dyy,dzz,stat=fail(2))
      deallocate (uxx,uyy,uzz,txx,tyy,tzz,stat=fail(3))
      deallocate (dxt,dyt,dzt,xxo,yyo,zzo,stat=fail(4))
      deallocate (vxo,vyo,vzo,tqx,tqy,tqz,stat=fail(5))
      deallocate (fmx,fmy,fmz,omxo,omyo,omzo,stat=fail(6))
      deallocate (gvxo,gvyo,gvzo,gcxo,gcyo,gczo,stat=fail(7))
      deallocate (xxt,yyt,zzt,qn0,qn1,qn2,qn3,stat=fail(8))
      
      return
      end subroutine nstq_b2
      
      subroutine nstq_h2
     x  (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x  ntcons,ntfree,ntpatm,mode,chit,conint,consv,elrc,engke,engrot,
     x  virlrc,press,quattol,sigma,taup,taut,temp,tolnce,tstep,
     x  vircom,vircon,volm)
      
c***********************************************************************
c     
c     dlpoly routine to integrate rigid body equations of motion
c     using implicit leapfrog quaternion algorithm
c     plus bond constraints- rigid body sites and constraint sites 
c     may coincide. 
c     
c     verlet leapfrog with Hoover like thermostat and barostat. 
c     (cell may change shape)
c     
c     parallel replicated data version : block data
c     
c     for systems using bond constraints
c     
c     tqx,tqy,tqz=torque in lab fixed frame (input)
c     omx,omy,omz=angular velocity in body fixed frame (principal axis)
c     rotinx,y,z =rotational inertia in body fixed frame
c     
c     copyright daresbury laboratory 1995
c     author      t.forester june  1995
c     amended     w.smith sep 1999 : euler equation
c     
c**********************************************************************
      
      implicit none
      
      integer, parameter :: nnn=16
      real(8), parameter :: pt5=0.5d0
      
      logical safe,safeq,lshmov,newstep,newjob,cycle
      integer fail,imcon,idnode,mxnode,natms,ngrp,nscons
      integer ntcons,ntfree,igrp,igrp1,igrp2,idum,ifre,ifre1,ifre2
      integer i,j,k,jg,ig,jr,id,mxshak1,icyc,ntpatm,mode
      real(8) engke,engrot,quattol,tolnce,tstep,vircom,vircon
      real(8) rot,strkin,strcon,strgrp,engtrn,vxt,vyt,vzt,trx,try,trz
      real(8) delx,dely,delz,vaa,vbb,vcc,viracc,com,fac,etadot
      real(8) chit,conint,consv,virlrc,elrc,press,sigma,taut,taup
      real(8) volm,eta0,etanew,cell0,volm0,elrc0,virlrc0,pmass,qmass
      real(8) totmas,chit0,xxa,yya,zza,chip,chitp,vom,engfke
      real(8) chitnew,vold,cons1,cons2,cons3,temp,uni
      real(8) stres0,engtot
      
      real(8), allocatable :: dens0(:)
      real(8), allocatable :: opx(:),opy(:),opz(:)
      real(8), allocatable :: oqx(:),oqy(:),oqz(:)
      real(8), allocatable :: dtx(:),dty(:),dtz(:)
      real(8), allocatable :: dxx(:),dyy(:),dzz(:)
      real(8), allocatable :: uxx(:),uyy(:),uzz(:)
      real(8), allocatable :: txx(:),tyy(:),tzz(:)
      real(8), allocatable :: xxt(:),yyt(:),zzt(:)
      real(8), allocatable :: dxt(:),dyt(:),dzt(:)
      real(8), allocatable :: xxo(:),yyo(:),zzo(:)
      real(8), allocatable :: vxo(:),vyo(:),vzo(:)
      real(8), allocatable :: fmx(:),fmy(:),fmz(:)
      real(8), allocatable :: tqx(:),tqy(:),tqz(:)
      real(8), allocatable :: omxo(:),omyo(:),omzo(:)
      real(8), allocatable :: gvxo(:),gvyo(:),gvzo(:)
      real(8), allocatable :: gcxo(:),gcyo(:),gczo(:)
      real(8), allocatable :: qn0(:),qn1(:),qn2(:),qn3(:)
      
      dimension rot(9),strkin(9),strcon(9),strgrp(9),fail(nnn),vom(3)
      dimension cell0(9),eta0(9),etanew(9),stres0(9),uni(9),com(3)
      
      save newjob,volm0,elrc0,virlrc0,dens0,pmass,qmass
      save igrp1,igrp2,ifre1,ifre2,totmas
      
      data uni/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/
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
      allocate (txx(mxatms),tyy(mxatms),tzz(mxatms),stat=fail(6))
      allocate (xxt(mxatms),yyt(mxatms),zzt(mxatms),stat=fail(7))
      allocate (dxt(mxcons),dyt(mxcons),dzt(mxcons),stat=fail(8))
      allocate (xxo(mxatms),yyo(mxatms),zzo(mxatms),stat=fail(9))
      allocate (vxo(msatms),vyo(msatms),vzo(msatms),stat=fail(10))
      allocate (fmx(msgrp),fmy(msgrp),fmz(msgrp),stat=fail(11))
      allocate (tqx(msgrp),tqy(msgrp),tqz(msgrp),stat=fail(12))
      allocate (omxo(msgrp),omyo(msgrp),omzo(msgrp),stat=fail(13))
      allocate (gvxo(msgrp),gvyo(msgrp),gvzo(msgrp),stat=fail(14))
      allocate (gcxo(msgrp),gcyo(msgrp),gczo(msgrp),stat=fail(15))
      allocate (qn0(msgrp),qn1(msgrp),qn2(msgrp),qn3(msgrp),
     x  stat=fail(16))
      do i=1,nnn
        if(fail(i).ne.0)call error(idnode,1710)
      enddo
      
      if(newjob) then
        
c     inertia parameter for Nose-Hoover thermostat
        
        qmass=2.0d0*sigma*taut**2
        pmass=2.0d0*sigma*taup**2
        
c     store initial values of volume, long range corrections etc
        
        volm0=volm
        elrc0=elrc
        virlrc0=virlrc
        fail(1)=0
        allocate (dens0(mxatyp),stat=fail(1))
        if(fail(1).ne.0)call error(idnode,1720)
        do i=1,ntpatm
          dens0(i)=dens(i)
        enddo
        
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
        
c     system total mass
        
        totmas=getmass(natms,idnode,mxnode)
        
        newjob=.false.
        
      endif
      
      safe=.false.
      cycle=.true.
      
c     ensure total momentum is zero
      
      call getvom(natms,idnode,mxnode,totmas,vom)
      
c     correction to velocities
      
      do ifre=ifre1,ifre2
        
        i=lstfre(ifre)
        vxx(i)=vxx(i)-vom(1)
        vyy(i)=vyy(i)-vom(2)
        vzz(i)=vzz(i)-vom(3)
        
      enddo
      
      do ig=igrp1,igrp2
        
        gvxx(ig)=gvxx(ig)-vom(1)
        gvyy(ig)=gvyy(ig)-vom(2)
        gvzz(ig)=gvzz(ig)-vom(3)
        
      enddo
      
c     store initial values of position and velocity
      
      do i=1,natms
        
        xxo(i)=xxx(i)
        yyo(i)=yyy(i)
        zzo(i)=zzz(i)
        
      enddo
      
      j=0
      do ifre=ifre1,ifre2
        
        j=j+1
        i=lstfre(ifre)
        vxo(j)=vxx(i)
        vyo(j)=vyy(i)
        vzo(j)=vzz(i)
        
      enddo
      
      jg=0
      do ig=igrp1,igrp2
        
        jg=jg+1
        omxo(jg)=omx(ig)
        omyo(jg)=omy(ig)
        omzo(jg)=omz(ig)
        gcxo(jg)=gcmx(ig)
        gcyo(jg)=gcmy(ig)
        gczo(jg)=gcmz(ig)
        gvxo(jg)=gvxx(ig)
        gvyo(jg)=gvyy(ig)
        gvzo(jg)=gvzz(ig)
        qn0(jg)=q0(ig)
        qn1(jg)=q1(ig)
        qn2(jg)=q2(ig)
        qn3(jg)=q3(ig)
        
      enddo
      
c     calculate centre of mass
      
      call getcom(natms,idnode,mxnode,totmas,com)
      
c     construct current bond vectors
      
      if(ntcons.gt.0)then
        
        do k=1,nscons
          
          i=listcon(k,2)
          j=listcon(k,3)
          
          dxx(k)=xxx(i)-xxx(j)
          dyy(k)=yyy(i)-yyy(j)
          dzz(k)=zzz(i)-zzz(j)
          
        enddo
        
c     periodic boundary condition for bond vectors
        
        call images(imcon,0,1,nscons,cell,dxx,dyy,dzz)
        
      endif
      
c     calculate atom displacement from rigid body com
      
      jg=0
      jr=0
      do ig=igrp1,igrp2
        
        jg=jg+1
        id=lstgtp(ig)
        
        do j=1,numgsit(id)
          
          jr=jr+1
          i=lstrgd(jr)
          
          dtx(jr)=xxo(i)-gcxo(jg)
          dty(jr)=yyo(i)-gcyo(jg)
          dtz(jr)=zzo(i)-gczo(jg)
          
        enddo
        
      enddo
      
c     minimum images
      
      call images(imcon,0,1,jr,cell,dtx,dty,dtz)
      
c     initial thermostat and barostat variables and new cell
      
      chit0=chit
      do i=1,9
        
        eta0(i)=eta(i)
        cell0(i)=cell(i)
        
      enddo
      
c     initialise constraint stress and virial
      
      vircon=0.d0
      do i=1,9
        
        strcns(i)=0.d0
        stres0(i)=stress(i)
        
      enddo
      
c     shake and thermostat iterations start here
      
      icyc=0
      mxshak1=mxshak
      if(ntcons.eq.0)mxshak1=4
      do while(cycle.and.icyc.le.mxshak1)
        
        icyc=icyc+1
        
c     restore cell vectors
        
        do i=1,9
          cell(i)=cell0(i)
        enddo
        
c     restore original quaternions for this step

        jg=0
        do ig=igrp1,igrp2

          jg=jg+1
          q0(ig)=qn0(jg)
          q1(ig)=qn1(jg)
          q2(ig)=qn2(jg)
          q3(ig)=qn3(jg)
          
        enddo

c     integrate unconstrained new positions
        
        j=0
        do ifre=ifre1,ifre2
          
          j=j+1
          i=lstfre(ifre)
          
c     advance velocity using leapfrog
          
          vxt=pt5*(vxx(i)+vxo(j))
          vyt=pt5*(vyy(i)+vyo(j))
          vzt=pt5*(vzz(i)+vzo(j))
          
          uxx(j)=vxo(j)+tstep*(fxx(i)*rmass(i)-
     x      (eta0(1)+chit0)*vxt-eta0(4)*vyt-eta0(7)*vzt)
          uyy(j)=vyo(j)+tstep*(fyy(i)*rmass(i)-
     x      eta0(2)*vxt-(eta0(5)+chit0)*vyt-eta0(8)*vzt)
          uzz(j)=vzo(j)+tstep*(fzz(i)*rmass(i)-
     x      eta0(3)*vxt-eta0(6)*vyt-(eta0(9)+chit0)*vzt)
          
c     advance positions using leapfrog
          
          xxa=(xxx(i)+xxo(i))*pt5-com(1)
          yya=(yyy(i)+yyo(i))*pt5-com(2)
          zza=(zzz(i)+zzo(i))*pt5-com(3)
          
          xxx(i)=xxo(i)+tstep*(uxx(j)+
     x      eta0(1)*xxa+eta0(4)*yya+eta0(7)*zza)
          yyy(i)=yyo(i)+tstep*(uyy(j)+
     x      eta0(2)*xxa+eta0(5)*yya+eta0(8)*zza)
          zzz(i)=zzo(i)+tstep*(uzz(j)+
     x      eta0(3)*xxa+eta0(6)*yya+eta0(9)*zza)
          
c     estimate full step  velocities
          
          vxx(i)=pt5*(uxx(j)+vxo(j))
          vyy(i)=pt5*(uyy(j)+vyo(j))
          vzz(i)=pt5*(uzz(j)+vzo(j))
          
        enddo
        
c     kinetic energy at current timestep
        
        engfke=getkinf(ntfree,idnode,mxnode)
        
c     kinetic contribution stress tensor
        
        call kinstressf(ntfree,idnode,mxnode,strkin)
        
c     restore half step velocities
        
        j=0
        do ifre=ifre1,ifre2
          
          j=j+1
          i=lstfre(ifre)
          vxx(i)=uxx(j)
          vyy(i)=uyy(j)
          vzz(i)=uzz(j)
          
        enddo
        
c     ********: rigid body motion - thermostated  :************
        
c     translational kinetic energy
        
        jg=0
        jr=0
        do ig=igrp1,igrp2
          
          jg=jg+1
          id=lstgtp(ig)
          
c     calculate centre of mass forces
          
          fmx(jg)=0.d0
          fmy(jg)=0.d0
          fmz(jg)=0.d0
          do j=1,numgsit(id)
            
            jr=jr+1
            i=lstrgd(jr)
            
            fmx(jg)=fmx(jg)+fxx(i)
            fmy(jg)=fmy(jg)+fyy(i)
            fmz(jg)=fmz(jg)+fzz(i)
            
          enddo
          
c     calculate thermostated velocities
          
          vxt=pt5*(gvxx(ig)+gvxo(jg))
          vyt=pt5*(gvyy(ig)+gvyo(jg))
          vzt=pt5*(gvzz(ig)+gvzo(jg))
          
          uxx(jg)=gvxo(jg)+tstep*(fmx(jg)/gmass(id)-
     x      (chit0+eta0(1))*vxt-eta0(4)*vyt-eta0(7)*vzt)
          uyy(jg)=gvyo(jg)+tstep*(fmy(jg)/gmass(id)-
     x      eta0(2)*vxt-(eta0(5)+chit0)*vyt-eta0(8)*vzt)
          uzz(jg)=gvzo(jg)+tstep*(fmz(jg)/gmass(id)-
     x      eta0(3)*vxt-eta0(6)*vyt-(eta0(9)+chit0)*vzt)
          
c     advance positions using leapfrog
          
          xxa=(gcmx(ig)+gcxo(jg))*pt5-com(1)
          yya=(gcmy(ig)+gcyo(jg))*pt5-com(2)
          zza=(gcmz(ig)+gczo(jg))*pt5-com(3)
          
          gcmx(ig)=gcxo(jg)+tstep*(uxx(jg)+
     x      eta0(1)*xxa+eta0(4)*yya+eta0(7)*zza)
          gcmy(ig)=gcyo(jg)+tstep*(uyy(jg)+
     x      eta0(2)*xxa+eta0(5)*yya+eta0(8)*zza)
          gcmz(ig)=gczo(jg)+tstep*(uzz(jg)+
     x      eta0(3)*xxa+eta0(6)*yya+eta0(9)*zza)
          
c     estimate full step velocities
          
          gvxx(ig)=pt5*(gvxo(jg)+uxx(jg))
          gvyy(ig)=pt5*(gvyo(jg)+uyy(jg))
          gvzz(ig)=pt5*(gvzo(jg)+uzz(jg))
          
        enddo
        
c     calculate kinetic energy
        
        engtrn=getkint(ngrp,idnode,mxnode)
        
c     kinetic contribution to stress tensor
        
        call kinstressg(ngrp,idnode,mxnode,strgrp)
        
c     restore half step velocities
        
        jg=0
        do ig=igrp1,igrp2
          
          jg=jg+1
          gvxx(ig)=uxx(jg)
          gvyy(ig)=uyy(jg)
          gvzz(ig)=uzz(jg)
          
        enddo
        
c     calculate rigid body stress tensor
        
        call bodystress(idnode,mxnode,ngrp,vircom,strbod,dtx,dty,dtz)
        
c     torques in lab frame
        
        jr=0
        jg=0
        do ig=igrp1,igrp2
          
          jg=jg+1
          id=lstgtp(ig)
          
          tqx(jg)=0.d0
          tqy(jg)=0.d0
          tqz(jg)=0.d0
          do j=1,numgsit(id)
            
            jr=jr+1
            i=lstrgd(jr)
            
            tqx(jg)=tqx(jg)+dty(jr)*fzz(i)-dtz(jr)*fyy(i)
            tqy(jg)=tqy(jg)+dtz(jr)*fxx(i)-dtx(jr)*fzz(i)
            tqz(jg)=tqz(jg)+dtx(jr)*fyy(i)-dty(jr)*fxx(i)
            
          enddo
          
c     current rotational matrix 
          
          call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)
          
c     angular velocity at time step n (first guess)
          
          opx(jg)=omxo(jg)
          opy(jg)=omyo(jg)
          opz(jg)=omzo(jg)
          
c     iterate angular velocity for time step n (e. yezdimer)
          
          do i=1,5
            
            trx=(tqx(jg)*rot(1)+tqy(jg)*rot(4)+tqz(jg)*rot(7))*
     x        rotinx(id,2)
     x        +(rotiny(id,1)-rotinz(id,1))*opy(jg)*opz(jg)*
     x        rotinx(id,2)
            try=(tqx(jg)*rot(2)+tqy(jg)*rot(5)+tqz(jg)*rot(8))*
     x        rotiny(id,2)
     x        +(rotinz(id,1)-rotinx(id,1))*opz(jg)*opx(jg)*
     x        rotiny(id,2)
            trz=(tqx(jg)*rot(3)+tqy(jg)*rot(6)+tqz(jg)*rot(9))*
     x        rotinz(id,2)
     x        +(rotinx(id,1)-rotiny(id,1))*opx(jg)*opy(jg)*
     x        rotinz(id,2)
            
c     improved angular velocity at time step n
            
            opx(jg)=omxo(jg)+pt5*tstep*trx
            opy(jg)=omyo(jg)+pt5*tstep*try
            opz(jg)=omzo(jg)+pt5*tstep*trz
            
          enddo
          
c     correction due to thermostat
          
          delx=tstep*(trx-chit0*pt5*(omx(ig)+omxo(jg)))
          dely=tstep*(try-chit0*pt5*(omy(ig)+omyo(jg)))
          delz=tstep*(trz-chit0*pt5*(omz(ig)+omzo(jg)))
          
c     angular velocity at time step n
          
          omx(ig)=omxo(jg)+delx*pt5
          omy(ig)=omyo(jg)+dely*pt5
          omz(ig)=omzo(jg)+delz*pt5
          
c     angular velocity at time step n+1/2
          
          uxx(jg)=omxo(jg)+delx
          uyy(jg)=omyo(jg)+dely
          uzz(jg)=omzo(jg)+delz
          
c     angular velocity at time step n+1 (needed for quat algorithm)
          
          oqx(jg)=omxo(jg)+delx*1.5d0
          oqy(jg)=omyo(jg)+dely*1.5d0
          oqz(jg)=omzo(jg)+delz*1.5d0
          
        enddo
        
c     rotational kinetic energy
        
        engrot=getkinr(ngrp,idnode,mxnode)
        
c     restore half step velocities
        
        jg=0
        do ig=igrp1,igrp2
          
          jg=jg+1
          opx(jg)=omx(ig)
          opy(jg)=omy(ig)
          opz(jg)=omz(ig)
          omx(ig)=uxx(jg)
          omy(ig)=uyy(jg)
          omz(ig)=uzz(jg)
          
        enddo
        
c     assign new quaternions
        
        call update_quaternions
     x    (safeq,igrp1,igrp2,tstep,quattol,opx,opy,opz,oqx,oqy,oqz)
        
c     minimum images of group positions and particle positions
        
        call images(imcon,0,1,natms,cell,xxx,yyy,zzz)
        call images(imcon,idnode,mxnode,ngrp,cell,gcmx,gcmy,gcmz)
        
c     new atomic positions for atoms in rigid bodies
        
        jg=0
        jr=0
        do ig=igrp1,igrp2
          
          jg=jg+1
          id=lstgtp(ig)
          
c     new rotational matrix
          
          call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)
          
          do j=1,numgsit(id)
            
            jr=jr+1
            i=lstrgd(jr)
            
            xxx(i)=rot(1)*gxx(id,j)+rot(2)*gyy(id,j)+rot(3)*gzz(id,j)
     x        +gcmx(ig)
            yyy(i)=rot(4)*gxx(id,j)+rot(5)*gyy(id,j)+rot(6)*gzz(id,j)
     x        +gcmy(ig)
            zzz(i)=rot(7)*gxx(id,j)+rot(8)*gyy(id,j)+rot(9)*gzz(id,j)
     x        +gcmz(ig)
            
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
        
c     merge new atomic coordinates 
        
        if(mxnode.gt.1)call merge1
     x    (idnode,mxnode,natms,lstme,xxx,yyy,zzz,buffer)
        
c     calculate total stress tensor
        
        do i=1,9
          stress(i)=stres0(i)+strcns(i)+strkin(i)+strgrp(i)+strbod(i)
        enddo
        
c     propagate eta
        
        fac=9.d0
        do i=1,9
          etanew(i)=eta(i)+tstep*((stress(i)-press*volm*uni(i))/pmass-
     x      chit0*eta0(i))
        enddo
        if(mode.gt.0)then
          fac=5.d0
          etanew(3)=0.d0
          etanew(6)=0.d0
          etanew(7)=0.d0
          etanew(8)=0.d0
          if(mode.lt.3)then
            fac=3.d0
            etanew(2)=0.d0
            etanew(4)=0.d0
            if(mode.eq.2)then
              fac=2.d0
              etanew(1)=0.5d0*(etanew(1)+etanew(5))
              etanew(5)=etanew(1)
            endif
          endif
        endif
        do i=1,9
          eta0(i)=pt5*(etanew(i)+eta(i)) 
        enddo
        
c     propagate chit
        
        etadot=sdot0(9,eta0,eta0)
        if(mode.eq.2)etadot=etadot-eta0(1)**2
        engke=engfke+engtrn
        engtot=engke+engrot
        chitp=(2.d0*(engtot-sigma)+pmass*etadot-fac*boltz*temp)/qmass
        chitnew=chit+tstep*chitp
        chit0=pt5*(chit+chitnew)
        
c     estimate new cell parameters
        
        call cell_propagate(tstep,cell,etanew)
        
        if(ntcons.gt.0) then
          
c     apply constraint correction
          
          newstep=.false.
          if(icyc.eq.1)newstep=.true.
          
          call qshake
     x      (newstep,safe,lshmov,idnode,imcon,mxnode,natms,
     x      nscons,tolnce,tstep,viracc,dxx,dyy,dzz,dxt,dyt,
     x      dzt,txx,tyy,tzz,xxt,yyt,zzt,strcon)
          
          if(abs(viracc).le.1.d-10)cycle=.false.
          
          vircon=vircon+viracc
          do i=1,9
            strcns(i)=strcns(i)+strcon(i)
          enddo
          
c     end of shake corrections
          
        endif
        
      enddo
      
c     update thermostat and barostat variables
      
      chit=chitnew
      do i=1,9
        eta(i)=etanew(i)
      enddo
      
c     update volume
      
      chip=eta(1)+eta(5)+eta(9)
      vold=volm
      volm=volm*exp(tstep*chip)
      
c     adjust long range corrections and number density
      
      elrc=elrc0*(volm0/volm)
      virlrc=virlrc0*(volm0/volm)
      
      do k=1,ntpatm
        dens(k)=dens0(k)*(volm0/volm)
      enddo
      
c     conserved quantity less kinetic and potential energy
      
      conint=conint+tstep*chit0*(qmass/taut**2+fac*boltz*temp)
      cons1=pt5*qmass*chit0**2
      cons2=press*vold
      etadot=sdot0(9,eta0,eta0)
      if(mode.eq.2)etadot=etadot-eta0(1)**2
      cons3=pt5*pmass*etadot
      consv=conint+cons1+cons2+cons3
      
      if(mxnode.gt.1) then
        
c     merge new group coordinates and velocities
        
        call merge(idnode,mxnode,ngrp,mxbuff,gcmx,gcmy,gcmz,buffer)
        call merge(idnode,mxnode,ngrp,mxbuff,gvxx,gvyy,gvzz,buffer)
        
c     merge new atomic velocities
        
        call merge1(idnode,mxnode,natms,lstme,vxx,vyy,vzz,buffer)
        
c     merge new quaternions
        
        call merge4(idnode,mxnode,ngrp,mxbuff,q0,q1,q2,q3,buffer)
        
      endif
      
c     ensure all atoms are within cell boundaries
      
      call images(imcon,0,1,natms,cell,xxx,yyy,zzz)
      
c     deallocate work arrays
      
      deallocate (opx,opy,opz,oqx,oqy,oqz,stat=fail(1))
      deallocate (dtx,dty,dtz,dxx,dyy,dzz,stat=fail(2))
      deallocate (uxx,uyy,uzz,txx,tyy,tzz,stat=fail(3))
      deallocate (dxt,dyt,dzt,xxo,yyo,zzo,stat=fail(4))
      deallocate (vxo,vyo,vzo,tqx,tqy,tqz,stat=fail(5))
      deallocate (fmx,fmy,fmz,omxo,omyo,omzo,stat=fail(6))
      deallocate (gvxo,gvyo,gvzo,gcxo,gcyo,gczo,stat=fail(7))
      deallocate (xxt,yyt,zzt,qn0,qn1,qn2,qn3,stat=fail(8))
      
      return
      end subroutine nstq_h2
      
      subroutine qshake
     x  (newstep,safe,lshmov,idnode,imcon,mxnode,natms,
     x  nscons,tolnce,tstep,vircon,dxx,dyy,dzz,dxt,dyt,
     x  dzt,txx,tyy,tzz,xxt,yyt,zzt,stresh)
      
c***********************************************************************
c     
c     dl_poly subroutine for appling bond constraint corrections after
c     atomic integration. Assumes rigid bodies connected by constraints
c     If this is not so use rdshake_1 instead
c     Must be used in conjunction with leapfrog integration algorithms
c     
c     copyright - daresbury laboratory 1995
c     author    - t. forester june 1995
c     
c***********************************************************************
      
      implicit none
      
      logical safe,lshmov,newstep,newjob
      integer fail,idnode,imcon,mxnode,natms,nscons,i,j,k
      integer ik,ig,id,jj
      real(8) tolnce,tstep,vircon,dxx,dyy,dzz,dxt,dyt,dzt,txx,tyy
      real(8) tzz,xxt,yyt,zzt,stresh,tstep2,esig,dis,dis2
      real(8) xxa,yya,zza,tax,tay,taz,doti,amti,amtj
      real(8) trx,try,trz,vix,viy,viz,vxi,vyi,vzi
      real(8) vjx,vjy,vjz,vxj,vyj,vzj,gamma,dli,dlj,rot
      
      real(8), allocatable :: redmass(:),esig1(:)
      
      dimension xxt(mxatms),yyt(mxatms),zzt(mxatms)
      dimension txx(mxatms),tyy(mxatms),tzz(mxatms)
      dimension dxx(mxcons),dyy(mxcons),dzz(mxcons)
      dimension dxt(mxcons),dyt(mxcons),dzt(mxcons)
      dimension stresh(9),rot(9)
      
      save newjob,esig1,redmass
      
      data newjob/.true./,fail/0/
      
      if(newjob)then
        
        allocate (redmass(mxcons),esig1(mxcons),stat=fail)
        if(fail.ne.0)call error(idnode,1610)
        newjob=.false.
        
      endif
      
c     constraint virial
      
      vircon=0.d0
      
c     accumulators for stress tensor
      
      do i=1,9
        stresh(i)=0.d0
      enddo
      
c     timestep squared
      
      tstep2=tstep*tstep
      
c     one iteration of constraint (shake) algorithm
      
      do k=1,nscons
        
        i=listcon(k,2)
        j=listcon(k,3)
        
c     calculate temporary bond vector
        
        dxt(k)=xxx(i)-xxx(j)
        dyt(k)=yyy(i)-yyy(j)
        dzt(k)=zzz(i)-zzz(j)
        
      enddo
      
c     periodic boundary condition
      
      call images(imcon,0,1,nscons,cell,dxt,dyt,dzt)
      
c     calculate maximum error in bondlength
      
      esig=0.d0
      do k=1,nscons
        
c     set bond parameter
        
        dis=prmcon(listcon(k,1))
        dis2=dis*dis
        esig1(k)=0.5d0*(dis2-(dxt(k)**2+dyt(k)**2+dzt(k)**2))/dis2
        esig=max(esig,abs(esig1(k)))

      enddo
      
c     global verification of convergence
      
      safe=(esig.lt.tolnce)
      
      if(mxnode.gt.1)call gstate(safe)
      
c     terminate iteration if all tolerances satisfied 
      
      if (.not.safe) then
        
c     initialise force increment arrays
        
        do i=1,natms
          
          xxt(i)=0.d0
          yyt(i)=0.d0
          zzt(i)=0.d0
          
        enddo
        
c     calculate constraint forces
        
        ik=0
        do k=1,nscons
          
          i=listcon(k,2)
          j=listcon(k,3)
          
c     assign effective reduced mass
          
          if(newstep) then
            
            ig=lstbod(i)
            
            if(ig.eq.0) then
              
              amti=rmass(i)
              
            else
              
              ik=ik+1
              id=lstgtp(ig)
              
              call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)
              
              jj=lstcsit(ik)
              
c     site to com in lab frame
              
              xxa=(gxx(id,jj)*rot(1)+gyy(id,jj)*rot(2)+
     x          gzz(id,jj)*rot(3))
              yya=(gxx(id,jj)*rot(4)+gyy(id,jj)*rot(5)+
     x          gzz(id,jj)*rot(6))
              zza=(gxx(id,jj)*rot(7)+gyy(id,jj)*rot(8)+
     x          gzz(id,jj)*rot(9))
              
c     find cross product between interatomic vector and vector to com
              
              tax=yya*dzz(k)-zza*dyy(k)
              tay=zza*dxx(k)-xxa*dzz(k)
              taz=xxa*dyy(k)-yya*dxx(k)
              
c     transform to body fixed frame
              
              trx=(tax*rot(1)+tay*rot(4)+taz*rot(7))*rotinx(id,2)
              try=(tax*rot(2)+tay*rot(5)+taz*rot(8))*rotiny(id,2)
              trz=(tax*rot(3)+tay*rot(6)+taz*rot(9))*rotinz(id,2)
              
c     direction of induced velocites in body frame
              
              vix=try*gzz(id,jj)-trz*gyy(id,jj)
              viy=trz*gxx(id,jj)-trx*gzz(id,jj)
              viz=trx*gyy(id,jj)-try*gxx(id,jj)
              
c     transform to lab frame
              
              vxi=vix*rot(1)+viy*rot(2)+viz*rot(3)
              vyi=vix*rot(4)+viy*rot(5)+viz*rot(6)
              vzi=vix*rot(7)+viy*rot(8)+viz*rot(9)
              
c     find dot product between induced translational and rotational velocities
              
              doti=abs(vxi*dxx(k)+vyi*dyy(k)+vzi*dzz(k))
              doti=doti/dis2
              
              amti=(1.d0/gmass(id)+doti)
              
            endif
            
            ig=lstbod(j)
            if(ig.eq.0) then
              
              amtj=rmass(j)
              
            else
              
              ik=ik+1
              id=lstgtp(ig)
              
              call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)
              
              jj=lstcsit(ik)
              
c     site to com in lab frame
              
              xxa=(gxx(id,jj)*rot(1)+gyy(id,jj)*rot(2)+
     x          gzz(id,jj)*rot(3))
              yya=(gxx(id,jj)*rot(4)+gyy(id,jj)*rot(5)+
     x          gzz(id,jj)*rot(6))
              zza=(gxx(id,jj)*rot(7)+gyy(id,jj)*rot(8)+
     x          gzz(id,jj)*rot(9))
              
c     find cross product between interatomic vector and vector to com
              
              tax =yya*dzz(k)-zza*dyy(k)
              tay =zza*dxx(k)-xxa*dzz(k)
              taz =xxa*dyy(k)-yya*dxx(k)
              
c     transform to body fixed frame
              
              trx=(tax*rot(1)+tay*rot(4)+taz*rot(7))*rotinx(id,2)
              try=(tax*rot(2)+tay*rot(5)+taz*rot(8))*rotiny(id,2)
              trz=(tax*rot(3)+tay*rot(6)+taz*rot(9))*rotinz(id,2)
              
c     direction of induced velocites in body frame
              
              vjx=try*gzz(id,jj)-trz*gyy(id,jj)
              vjy=trz*gxx(id,jj)-trx*gzz(id,jj)
              vjz=trx*gyy(id,jj)-try*gxx(id,jj)
              
c     transform to lab frame
              
              vxj=vjx*rot(1)+vjy*rot(2)+vjz*rot(3)
              vyj=vjx*rot(4)+vjy*rot(5)+vjz*rot(6)
              vzj=vjx*rot(7)+vjy*rot(8)+vjz*rot(9)
              
c     find dot product between induced translational and rotational velocities
              
              doti=abs(vxj*dxx(k)+vyj*dyy(k)+vzj*dzz(k))
              doti=doti/dis2
              
              amtj=(1.d0/gmass(id)+doti)
              
            endif
            
            redmass(k)=1.d0/(amti+amtj)/tstep2
            
          endif
          
c     constraint force parameter 
          
          gamma=esig1(k)*redmass(k)
          
c     accumulate bond virial
          
          vircon=vircon-gamma*(dxx(k)**2+dyy(k)**2+dzz(k)**2)
          
          stresh(1)=stresh(1)+gamma*dxx(k)*dxx(k)
          stresh(2)=stresh(2)+gamma*dxx(k)*dyy(k)
          stresh(3)=stresh(3)+gamma*dxx(k)*dzz(k)
          stresh(5)=stresh(5)+gamma*dyy(k)*dyy(k)
          stresh(6)=stresh(6)+gamma*dyy(k)*dzz(k)
          stresh(9)=stresh(9)+gamma*dzz(k)*dzz(k)
          
c     improved atomic force
          
          xxt(i)=xxt(i)+dxx(k)*gamma
          yyt(i)=yyt(i)+dyy(k)*gamma
          zzt(i)=zzt(i)+dzz(k)*gamma
          
          xxt(j)=xxt(j)-dxx(k)*gamma
          yyt(j)=yyt(j)-dyy(k)*gamma
          zzt(j)=zzt(j)-dzz(k)*gamma
          
        enddo
        
c     transport temporary positions to other nodes
        
        if(mxnode.gt.1)then
          
          if(lshmov) call shmove
     x      (idnode,mxnode,natms,lashap,lishap,xxt,yyt,zzt,
     x      txx,tyy,tzz,buffer)
          
        endif
        
        do k=1,nscons
          
          i=listcon(k,2)
          j=listcon(k,3)
          
          dli=1.d0/dble(listme(i))
          dlj=1.d0/dble(listme(j))
          
          fxx(i)=fxx(i)+xxt(i)*dli
          fyy(i)=fyy(i)+yyt(i)*dli
          fzz(i)=fzz(i)+zzt(i)*dli
          fxx(j)=fxx(j)+xxt(j)*dlj
          fyy(j)=fyy(j)+yyt(j)*dlj
          fzz(j)=fzz(j)+zzt(j)*dlj
          
        enddo
        
c     complete (symmetrical) stress tensor
        
        stresh(4)=stresh(2)
        stresh(7)=stresh(3)
        stresh(8)=stresh(6)
        
c     splice force arrays across nodes
        
        if(mxnode.gt.1)then
          
          buffer(1)=vircon
          call gdsum(buffer(1),1,buffer(2))
          vircon=buffer(1)
          call gdsum(stresh,9,buffer)
          call splice 
     x      (idnode,natms,listme,listot,fxx,fyy,fzz,buffer)
          
        endif
        
      endif
      
      return
      end subroutine qshake
      
      end module lf_rotation2_module
