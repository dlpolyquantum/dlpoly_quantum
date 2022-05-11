      module lf_rotation1_module
      
c***********************************************************************
c     
c     dl_poly module 1 for verlet leap frog rotational integration 
c     schemes
c     copyright - daresbury laboratory
c     author    - w. smith    aug 2006
c     
c***********************************************************************
      
      use config_module
      use ensemble_tools_module
      use error_module
      use lf_motion_module
      use property_module
      use rigid_body_module
      use setup_module
      use shake_module
      use site_module
      use utility_module
      
      contains
      
      subroutine update_quaternions
     x  (safeq,igrp1,igrp2,tstep,quattol,opx,opy,opz,oqx,oqy,oqz)

c***********************************************************************
c     
c     dlpoly routine to update the quaternion arrays as part of
c     the leapfrog algorithm
c     
c     copyright daresbury laboratory
c     author   -  w.smith october 2005
c     based on -  t.forester oct. 1993
c     
c**********************************************************************

      implicit none

      real(8), parameter :: pt5=0.5d0

      logical safeq
      integer igrp1,igrp2,jg,itq,ig
      real(8) qn0,qn1,qn2,qn3,qn0a,qn1a,qn2a,qn3a,qn0b,qn1b,qn2b,qn3b
      real(8) rnorm,tstep,quattol,eps
      real(8) opx(msgrp),opy(msgrp),opz(msgrp)
      real(8) oqx(msgrp),oqy(msgrp),oqz(msgrp)

      jg=0
      safeq=.true.

      do ig=igrp1,igrp2

        jg=jg+1

c     first iteration of new quaternions (lab fixed)
        
        qn0=q0(ig)+(-q1(ig)*opx(jg)-q2(ig)*opy(jg)-q3(ig)*opz(jg))
     x    *tstep*pt5
        qn1=q1(ig)+( q0(ig)*opx(jg)-q3(ig)*opy(jg)+q2(ig)*opz(jg))
     x    *tstep*pt5
        qn2=q2(ig)+( q3(ig)*opx(jg)+q0(ig)*opy(jg)-q1(ig)*opz(jg))
     x    *tstep*pt5
        qn3=q3(ig)+(-q2(ig)*opx(jg)+q1(ig)*opy(jg)+q0(ig)*opz(jg))
     x    *tstep*pt5
        
        qn0b=0.d0
        qn1b=0.d0
        qn2b=0.d0
        qn3b=0.d0
        
        itq=0
        eps=1.0d9
        do while((itq.lt.mxquat).and.(eps.gt.quattol))
          
          itq=itq+1
          
          qn0a=pt5*(-q1(ig)*opx(jg)-q2(ig)*opy(jg)-q3(ig)*opz(jg))
     x      +pt5*(-qn1*oqx(jg)-qn2*oqy(jg)-qn3*oqz(jg))
          qn1a=pt5*(  q0(ig)*opx(jg)-q3(ig)*opy(jg)+q2(ig)*opz(jg))
     x      +   pt5*( qn0*oqx(jg)-qn3*oqy(jg)+qn2*oqz(jg))
          qn2a=pt5*(  q3(ig)*opx(jg)+q0(ig)*opy(jg)-q1(ig)*opz(jg))
     x      +   pt5*( qn3*oqx(jg)+qn0*oqy(jg)-qn1*oqz(jg))
          qn3a=pt5*(-q2(ig)*opx(jg)+q1(ig)*opy(jg)+q0(ig)*opz(jg))
     x      +   pt5*(-qn2*oqx(jg)+qn1*oqy(jg)+qn0*oqz(jg))
          
          qn0=q0(ig)+pt5*qn0a*tstep
          qn1=q1(ig)+pt5*qn1a*tstep
          qn2=q2(ig)+pt5*qn2a*tstep
          qn3=q3(ig)+pt5*qn3a*tstep
          
          rnorm=1.d0/sqrt(qn0**2+qn1**2+qn2**2+qn3**2)
          qn0=qn0*rnorm
          qn1=qn1*rnorm
          qn2=qn2*rnorm
          qn3=qn3*rnorm
          
c     convergence test 
          
          eps=sqrt(((qn0a-qn0b)**2+(qn1a-qn1b)**2+(qn2a-qn2b)**2
     x      +(qn3a-qn3b)**2)*tstep**2)
          
          qn0b=qn0a
          qn1b=qn1a
          qn2b=qn2a
          qn3b=qn3a
          
        enddo
        
        if(itq.ge.mxquat) safeq=.false.
        
c     store new quaternions
        
        q0(ig)=qn0
        q1(ig)=qn1
        q2(ig)=qn2
        q3(ig)=qn3
        
      enddo

      return
      end subroutine update_quaternions
      
      subroutine nveq_1
     x  (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x  ntcons,ntfree,engke,engrot,quattol,tolnce,tstep,vircom,
     x  vircon)
      
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

      integer, parameter :: nnn=12
      real(8), parameter :: pt5=0.5d0

      logical safe,safeq,lshmov,newjob
      integer imcon,idnode,mxnode,natms,ngrp,nscons,ntcons,ntfree
      integer fail,i,igrp,igrp1,igrp2,idum,ifre1,ifre2,j,ifre
      integer jg,ig,k,id,jr
      real(8) engke,engrot,quattol,tolnce,tstep,vircom,vircon
      real(8) strkin,rot,rstep,rtsq,engtrn,vaa,vbb,vcc
      real(8) trx,try,trz,delx,dely,delz,engfke
      real(8) strgrp,tqx,tqy,tqz,fmx,fmy,fmz

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

      dimension strkin(9),strgrp(9),rot(9),fail(nnn)
      
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

      if(ntcons.eq.0) safe=.true.
      if(ntcons.gt.0) then

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

c     apply constraint correction
        
        call rdshake_1
     x    (safe,lshmov,idnode,imcon,mxnode,natms,nscons,
     x    tolnce,tstep,vircon,dxx,dyy,dzz,dxt,dyt,dzt,
     x    txx,tyy,tzz,xxt,yyt,zzt,strcns)
        
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

c     calculate rigid body kinetic stress tensor
      
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

c       iterate angular velocity for time step n (e. yezdimer)

        do i=1,5
        
          trx=(tqx*rot(1)+tqy*rot(4)+tqz*rot(7))*rotinx(id,2)
     x     +(rotiny(id,1)-rotinz(id,1))*opy(jg)*opz(jg)*rotinx(id,2)
          try=(tqx*rot(2)+tqy*rot(5)+tqz*rot(8))*rotiny(id,2)
     x     +(rotinz(id,1)-rotinx(id,1))*opz(jg)*opx(jg)*rotiny(id,2)
          trz=(tqx*rot(3)+tqy*rot(6)+tqz*rot(9))*rotinz(id,2)
     x     +(rotinx(id,1)-rotiny(id,1))*opx(jg)*opy(jg)*rotinz(id,2)

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
        stress(i)=stress(i)+strkin(i)+strgrp(i)+strcns(i)+strbod(i)
      enddo
      
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
     x     +gcmx(ig)
          yyy(i)=rot(4)*gxx(id,j)+rot(5)*gyy(id,j)+rot(6)*gzz(id,j)
     x     +gcmy(ig)
          zzz(i)=rot(7)*gxx(id,j)+rot(8)*gyy(id,j)+rot(9)*gzz(id,j)
     x     +gcmz(ig)

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
        if(ntcons.gt.0)
     x    call merge1(idnode,mxnode,natms,lstme,fxx,fyy,fzz,buffer)
        
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
      
      return
      end subroutine nveq_1

      subroutine nvtq_b1
     x  (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x  ntcons,ntfree,engke,engrot,quattol,sigma,taut,tolnce,
     x  tstep,vircom,vircon)
      
c***********************************************************************
c     
c     dlpoly routine to integrate rigid body equations of motion
c     using implicit leapfrog quaternion algorithm
c     plus bond constraints - PROVIDED rigid body sites
c     and constraint sites do not coincide.
c     
c     dl_poly subroutine for integrating newtonian equations of
c     motion in molecular dynamics - verlet leapfrog with Berendsen
c     thermostat.
c     
c     parallel replicated data version : block data
c     
c     for systems using bond constraints
c     
c     parallel replicated data version : block data
c     
c     tqx,tqy,tqz = torque in lab fixed frame
c     omx,omy,omz = angular velocity in body fixed frame (principle axis)
c     rotinx,y,z = rotational inertia in body fixed frame
c     
c     copyright daresbury laboratory 1993
c     author      t.forester october 1993
c     amended     t.forester dec 1994 : block data
c     amended     w.smith sep 1999 : euler equation
c     
c**********************************************************************
      
      implicit none
      
      integer, parameter :: nnn=16
      real(8), parameter :: pt5=0.5d0

      logical safe,lshmov,safeq,newjob
      integer imcon,idnode,mxnode,natms,ngrp,nscons,ntcons,ntfree
      integer i,fail,igrp,igrp1,igrp2,ifre,ifre1,ifre2,jg,ig
      integer j,k,jr,id,mxiter,iter,idum
      real(8) engke,engrot,quattol,sigma,taut,tolnce,tstep,vircom
      real(8) vircon,strkin,strgrp,rot,rstep,rtsq
      real(8) engtrn,trx,try,trz,chit0,rgmas,engfke
      real(8) vaa,vbb,vcc,engtot,viracc,strcon
      
      real(8), allocatable :: opx(:),opy(:),opz(:)
      real(8), allocatable :: oqx(:),oqy(:),oqz(:)
      real(8), allocatable :: dtx(:),dty(:),dtz(:)
      real(8), allocatable :: xxt(:),yyt(:),zzt(:)
      real(8), allocatable :: dxx(:),dyy(:),dzz(:)
      real(8), allocatable :: uxx(:),uyy(:),uzz(:)
      real(8), allocatable :: txx(:),tyy(:),tzz(:)
      real(8), allocatable :: dxt(:),dyt(:),dzt(:)
      real(8), allocatable :: xxo(:),yyo(:),zzo(:)
      real(8), allocatable :: xx1(:),yy1(:),zz1(:)
      real(8), allocatable :: vxo(:),vyo(:),vzo(:)
      real(8), allocatable :: fmx(:),fmy(:),fmz(:)
      real(8), allocatable :: tqx(:),tqy(:),tqz(:)
      real(8), allocatable :: gcxo(:),gcyo(:),gczo(:)
      real(8), allocatable :: omxo(:),omyo(:),omzo(:)
      real(8), allocatable :: gvxo(:),gvyo(:),gvzo(:)

      dimension strkin(9),strgrp(9),strcon(9),rot(9),fail(nnn)
      
      save igrp1,igrp2,ifre1,ifre2,newjob
      
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
      allocate (fmx(msgrp),fmy(msgrp),fmz(msgrp),stat=fail(12))
      allocate (tqx(msgrp),tqy(msgrp),tqz(msgrp),stat=fail(13))
      allocate (omxo(msgrp),omyo(msgrp),omzo(msgrp),stat=fail(14))
      allocate (gvxo(msgrp),gvyo(msgrp),gvzo(msgrp),stat=fail(15))
      allocate (gcxo(msgrp),gcyo(msgrp),gczo(msgrp),stat=fail(16))
      do i=1,nnn
        if(fail(i).ne.0)call error(idnode,1510)
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

c     initialise constraint virial terms
      
      vircon=0.d0
      do i=1,9
        strcns(i)=0.d0
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
        omxo(jg)=omx(ig)
        omyo(jg)=omy(ig)
        omzo(jg)=omz(ig)
        gcxo(jg)=gcmx(ig)
        gcyo(jg)=gcmy(ig)
        gczo(jg)=gcmz(ig)
        gvxo(jg)=gvxx(ig)
        gvyo(jg)=gvyy(ig)
        gvzo(jg)=gvzz(ig)

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

c     periodic boundary condition for displacements
      
      call images(imcon,0,1,jr,cell,dtx,dty,dtz)

c     estimate velocity and temperature at half-time step
      
      j=0
      do ifre=ifre1,ifre2
        
        j=j+1
        i=lstfre(ifre)

        vxx(i)=vxo(j)+pt5*tstep*rmass(i)*fxx(i)
        vyy(i)=vyo(j)+pt5*tstep*rmass(i)*fyy(i)
        vzz(i)=vzo(j)+pt5*tstep*rmass(i)*fzz(i)

      enddo

c     kinetic energy of free atoms

      engfke=getkinf(ntfree,idnode,mxnode)

c     estimate kinetic energy of rigid body com
      
      jg=0
      jr=0
      do ig=igrp1,igrp2
        
        jg=jg+1
        id=lstgtp(ig)

c     calculate net force on rigid body
        
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

c     centre of mass velocities at half-step
        
        gvxx(ig)=gvxo(jg)+pt5*tstep*fmx(jg)/gmass(id)
        gvyy(ig)=gvyo(jg)+pt5*tstep*fmy(jg)/gmass(id)
        gvzz(ig)=gvzo(jg)+pt5*tstep*fmz(jg)/gmass(id)

      enddo

c     translation kinetic energy of rigid bodies

      engtrn=getkint(ngrp,idnode,mxnode)
      
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
     x      rotinx(id,2)
     x     +(rotiny(id,1)-rotinz(id,1))*opy(jg)*opz(jg)*rotinx(id,2)
          try=(tqx(jg)*rot(2)+tqy(jg)*rot(5)+tqz(jg)*rot(8))*
     x      rotiny(id,2)
     x     +(rotinz(id,1)-rotinx(id,1))*opz(jg)*opx(jg)*rotiny(id,2)
          trz=(tqx(jg)*rot(3)+tqy(jg)*rot(6)+tqz(jg)*rot(9))*
     x      rotinz(id,2)
     x     +(rotinx(id,1)-rotiny(id,1))*opx(jg)*opy(jg)*rotinz(id,2)

c     improved angular velocity at time step n
          
          opx(jg)=omx(ig)+pt5*tstep*trx
          opy(jg)=omy(ig)+pt5*tstep*try
          opz(jg)=omz(ig)+pt5*tstep*trz
        
        enddo
        
c     store angular velocity at timestep n
        
        omx(ig)=opx(jg)
        omy(ig)=opy(jg)
        omz(ig)=opz(jg)
        
      enddo
      
c     rigid body rotational kinetic energy

      engrot=getkinr(ngrp,idnode,mxnode)

c     temperature scaling  coefficient - taut is the relaxation time
      
      engtot=engfke+engrot+engtrn
      chit0= sqrt(1.d0+tstep/taut*(sigma/engtot-1.d0))

c     begin iterations !!-----------------------------------------------

      mxiter=2
      if(ntcons.eq.0) mxiter=mxiter-1
      
      do iter=1,mxiter

c     unconstrained new positions
        
        j=0
        do ifre=ifre1,ifre2

          j=j+1
          i=lstfre(ifre)

c     advance velocity using leapfrog

          uxx(i)=(vxo(j)+tstep*rmass(i)*fxx(i))*chit0
          uyy(i)=(vyo(j)+tstep*rmass(i)*fyy(i))*chit0
          uzz(i)=(vzo(j)+tstep*rmass(i)*fzz(i))*chit0

c     update positions 
          
          xxx(i)=tstep*uxx(i)+xxo(j)
          yyy(i)=tstep*vyy(i)+yyo(j)
          zzz(i)=tstep*uzz(i)+zzo(j)

        enddo

        if(ntcons.eq.0) safe=.true.
        if(ntcons.gt.0.and.iter.eq.1) then

c     store integrated positions

          j=0
          do ifre=ifre1,ifre2

            j=j+1
            i=lstfre(ifre)
            xx1(j)=xxx(i)
            yy1(j)=yyy(i)
            zz1(j)=zzz(i)

          enddo

c     global exchange of configuration data
          
          if(mxnode.gt.1)call merge1
     x      (idnode,mxnode,natms,lstme,xxx,yyy,zzz,buffer)

c     apply constraint correction
          
          call rdshake_1
     x      (safe,lshmov,idnode,imcon,mxnode,natms,nscons,
     x      tolnce,tstep,viracc,dxx,dyy,dzz,dxt,dyt,dzt,
     x      txx,tyy,tzz,xxt,yyt,zzt,strcon)
          
c     contribution to constraint virial 
          
          vircon=vircon+viracc
          do i=1,9
            strcns(i)=strcns(i)+strcon(i)
          enddo
        
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
        
        if(iter.eq.mxiter)call kinstressf(ntfree,idnode,mxnode,strkin)
        
c     restore free atom half step velocity
        
        do ifre=ifre1,ifre2
          
          i=lstfre(ifre)
          vxx(i)=uxx(i)
          vyy(i)=uyy(i)
          vzz(i)=uzz(i)
          
        enddo
        
c     ********: rigid body motion - thermostated  :************
        
c     ***** step 1 : integrate centre of mass motion *********
        
        jg =0
        do ig=igrp1,igrp2
          
          jg=jg+1
          id=lstgtp(ig)

c     calculate thermostated velocities
          
          rgmas=1.d0/gmass(id)
          uxx(ig)=(gvxo(jg)+tstep*(fmx(jg)*rgmas))*chit0
          uyy(ig)=(gvyo(jg)+tstep*(fmy(jg)*rgmas))*chit0
          uzz(ig)=(gvzo(jg)+tstep*(fmz(jg)*rgmas))*chit0

c     update positions 
          
          gcmx(ig)=gcxo(jg)+tstep*uxx(ig)
          gcmy(ig)=gcyo(jg)+tstep*uyy(ig)
          gcmz(ig)=gczo(jg)+tstep*uzz(ig)
          
c     calculate full step velocities
          
          gvxx(ig)=0.5d0*(gvxo(jg)+uxx(ig))
          gvyy(ig)=0.5d0*(gvyo(jg)+uyy(ig))
          gvzz(ig)=0.5d0*(gvzo(jg)+uzz(ig))
          
        enddo
        
c     calculate kinetic energy
        
        engtrn=getkint(ngrp,idnode,mxnode)
        
        if(iter.eq.mxiter)call kinstressg(ngrp,idnode,mxnode,strgrp)
        
c     restore half step velocities 

        do ig=igrp1,igrp2

          gvxx(ig)=uxx(ig)
          gvyy(ig)=uyy(ig)
          gvzz(ig)=uzz(ig)
          
        enddo
        
c     ****** step 2 : integrate rotational motion **********
        
        jg=0
        safeq=.true.
        do ig=igrp1,igrp2
          
          jg=jg+1
          id=lstgtp(ig)

c     current rotational matrix 
          
          call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)
        
          trx=(tqx(jg)*rot(1)+tqy(jg)*rot(4)+tqz(jg)*rot(7))*
     x      rotinx(id,2)
     x     +(rotiny(id,1)-rotinz(id,1))*opy(jg)*opz(jg)*
     x      rotinx(id,2)
          try=(tqx(jg)*rot(2)+tqy(jg)*rot(5)+tqz(jg)*rot(8))*
     x      rotiny(id,2)
     x     +(rotinz(id,1)-rotinx(id,1))*opz(jg)*opx(jg)*
     x      rotiny(id,2)
          trz=(tqx(jg)*rot(3)+tqy(jg)*rot(6)+tqz(jg)*rot(9))*
     x      rotinz(id,2)
     x     +(rotinx(id,1)-rotiny(id,1))*opx(jg)*opy(jg)*
     x      rotinz(id,2)

c     scaled angular velocity at time step n
            
          omx(ig)=(omxo(jg)+pt5*tstep*trx)*chit0
          omy(ig)=(omyo(jg)+pt5*tstep*try)*chit0
          omz(ig)=(omzo(jg)+pt5*tstep*trz)*chit0

c     angular velocity at time step n+1/2
          
          uxx(ig)=(omxo(jg)+tstep*trx)*chit0
          uyy(ig)=(omyo(jg)+tstep*try)*chit0
          uzz(ig)=(omzo(jg)+tstep*trz)*chit0

c     angular velocity at time step n+1 (needed for quat algorithm)
          
          oqx(jg)=(omxo(jg)+1.5d0*tstep*trx)*chit0
          oqy(jg)=(omyo(jg)+1.5d0*tstep*try)*chit0
          oqz(jg)=(omzo(jg)+1.5d0*tstep*trz)*chit0

        enddo
        
c     calculate rigid body rotational energy
        
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

c     total translational kinetic energy
      
        engke=engtrn+engfke
      
c     total kinetic energy
        
        engtot=engke+engrot
        
c     improved prediction of chit

        chit0=sqrt(1.d0+tstep/taut*(sigma/engtot-1.d0))

c     end of thermostat/barostat iterations

      enddo

c     assign new quaternions

      call update_quaternions
     x  (safeq,igrp1,igrp2,tstep,quattol,opx,opy,opz,oqx,oqy,oqz)

c     minimum images 
      
      call images(imcon,0,1,natms,cell,xxx,yyy,zzz)
      call images(imcon,idnode,mxnode,ngrp,cell,gcmx,gcmy,gcmz)

c     new atomic positions for atoms in rigid bodies
      
      jr=0
      do ig=igrp1,igrp2

        id=lstgtp(ig)

c     new rotational matrix
        
        call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)
        
        do j=1,numgsit(id)
          
          jr=jr+1
          i=lstrgd(jr)
          
          xxx(i)=rot(1)*gxx(id,j)+rot(2)*gyy(id,j)+rot(3)*gzz(id,j)
     x     +gcmx(ig)
          yyy(i)=rot(4)*gxx(id,j)+rot(5)*gyy(id,j)+rot(6)*gzz(id,j)
     x     +gcmy(ig)
          zzz(i)=rot(7)*gxx(id,j)+rot(8)*gyy(id,j)+rot(9)*gzz(id,j)
     x     +gcmz(ig)

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
        if(ntcons.gt.0)
     x    call merge1(idnode,mxnode,natms,lstme,fxx,fyy,fzz,buffer)
      
      endif

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
      deallocate (xx1,yy1,zz1,vxo,vyo,vzo,stat=fail(5))
      deallocate (fmx,fmy,fmz,tqx,tqy,tqz,stat=fail(6))
      deallocate (omxo,omyo,omzo,gcxo,gcyo,gczo,stat=fail(7))
      deallocate (xxt,yyt,zzt,gvxo,gvyo,gvzo,stat=fail(8))
      
      return
      end subroutine nvtq_b1

      subroutine nvtq_h1
     x  (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x  ntcons,ntfree,chit,consv,conint,engke,engrot,quattol,
     x  sigma,taut,tolnce,tstep,vircom,vircon)
      
c***********************************************************************
c     
c     dlpoly routine to integrate rigid body equations of motion
c     using implicit leapfrog quaternion algorithm
c     plus bond constraints - PROVIDED rigid body sites
c     and constraint sites do not coincide.
c     
c     dl_poly subroutine for integrating newtonian equations of
c     motion in molecular dynamics - verlet leapfrog with Hoover 
c     thermostat.
c     
c     for systems using bond constraints
c     
c     parallel replicated data version : block data
c     
c     tqx,tqy,tqz = torque in lab fixed frame
c     omx,omy,omz = angular velocity in body fixed frame (principle axis)
c     rotinx,y,z  = rotational inertia in body fixed frame
c     
c     copyright daresbury laboratory 1993
c     author      t.forester october 1993
c     amended     t.forester dec 1994 : block data
c     amended     w.smith nov 2005
c     
c**********************************************************************
      
      implicit none
      
      integer, parameter :: nnn=16
      real(8), parameter :: pt5=0.5d0
      
      logical safe,lshmov,safeq,newjob
      integer imcon,idnode,mxnode,natms,ngrp,nscons,ntcons,ntfree
      integer fail,i,idum,igrp,igrp1,igrp2,ifre,ifre1,ifre2
      integer j,k,ig,jg,jr,id,iter,mxiter
      real(8) chit,consv,conint,engke,engrot,quattol,sigma,taut
      real(8) tolnce,tstep,vircom,vircon,strkin,strgrp,strcon,rot
      real(8) rstep,rtsq,qmass,engtrn,cons1,engtot,vaa,vbb,vcc
      real(8) chit0,chitnew,chitp,viracc,rgmas,trx,try,trz,delx
      real(8) dely,delz,engfke

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
      real(8), allocatable :: fmx(:),fmy(:),fmz(:)
      real(8), allocatable :: tqx(:),tqy(:),tqz(:)
      real(8), allocatable :: omxo(:),omyo(:),omzo(:)
      real(8), allocatable :: gvxo(:),gvyo(:),gvzo(:)
      real(8), allocatable :: gcxo(:),gcyo(:),gczo(:)

      dimension strkin(9),strgrp(9),strcon(9),rot(9),fail(nnn)
      
      save igrp1,igrp2,ifre1,ifre2,qmass,newjob
      
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
      allocate (fmx(msgrp),fmy(msgrp),fmz(msgrp),stat=fail(12))
      allocate (tqx(msgrp),tqy(msgrp),tqz(msgrp),stat=fail(13))
      allocate (omxo(msgrp),omyo(msgrp),omzo(msgrp),stat=fail(14))
      allocate (gvxo(msgrp),gvyo(msgrp),gvzo(msgrp),stat=fail(15))
      allocate (gcxo(msgrp),gcyo(msgrp),gczo(msgrp),stat=fail(16))
      do i=1,nnn
        if(fail(i).ne.0)call error(idnode,1520)
      enddo

      if(newjob)then
        
c     inertia parameter for Nose-Hoover thermostat
        
        qmass=2.0d0*sigma*taut**2
        
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

c     initialise constraint virial terms
      
      vircon=0.d0
      do i=1,9
        strcns(i)=0.d0
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
        omxo(jg)=omx(ig)
        omyo(jg)=omy(ig)
        omzo(jg)=omz(ig)
        gcxo(jg)=gcmx(ig)
        gcyo(jg)=gcmy(ig)
        gczo(jg)=gcmz(ig)
        gvxo(jg)=gvxx(ig)
        gvyo(jg)=gvyy(ig)
        gvzo(jg)=gvzz(ig)

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

c     periodic boundary condition for displacements
      
      call images(imcon,0,1,jr,cell,dtx,dty,dtz)

c     estimate velocity and temperature at half-time step
      
      j=0
      do ifre=ifre1,ifre2
        
        j=j+1
        i=lstfre(ifre)

        vxx(i)=vxo(j)+pt5*tstep*rmass(i)*fxx(i)
        vyy(i)=vyo(j)+pt5*tstep*rmass(i)*fyy(i)
        vzz(i)=vzo(j)+pt5*tstep*rmass(i)*fzz(i)

      enddo

c     kinetic energy of free atoms

      engfke=getkinf(ntfree,idnode,mxnode)
      
c     estimate kinetic energy of rigid body com
      
      jg=0
      jr=0
      do ig=igrp1,igrp2
        
        jg=jg+1

c     calculate net force on rigid body
        
        fmx(jg)=0.d0
        fmy(jg)=0.d0
        fmz(jg)=0.d0
        id=lstgtp(ig)
        do j=1,numgsit(id)
          
          jr=jr+1
          i=lstrgd(jr)
          fmx(jg)=fmx(jg)+fxx(i)
          fmy(jg)=fmy(jg)+fyy(i)
          fmz(jg)=fmz(jg)+fzz(i)
          
        enddo

c     centre of mass velocities at half-step
        
        gvxx(ig)=gvxo(jg)+pt5*tstep*fmx(jg)/gmass(id)
        gvyy(ig)=gvyo(jg)+pt5*tstep*fmy(jg)/gmass(id)
        gvzz(ig)=gvzo(jg)+pt5*tstep*fmz(jg)/gmass(id)

      enddo

c     translation kinetic energy of rigid bodies

      engtrn=getkint(ngrp,idnode,mxnode)

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
     x      rotinx(id,2)
     x     +(rotiny(id,1)-rotinz(id,1))*opy(jg)*opz(jg)*
     x      rotinx(id,2)
          try=(tqx(jg)*rot(2)+tqy(jg)*rot(5)+tqz(jg)*rot(8))*
     x      rotiny(id,2)
     x     +(rotinz(id,1)-rotinx(id,1))*opz(jg)*opx(jg)*
     x      rotiny(id,2)
          trz=(tqx(jg)*rot(3)+tqy(jg)*rot(6)+tqz(jg)*rot(9))*
     x      rotinz(id,2)
     x     +(rotinx(id,1)-rotiny(id,1))*opx(jg)*opy(jg)*
     x      rotinz(id,2)

c     improved angular velocity at time step n
          
          opx(jg)=omx(ig)+pt5*tstep*trx
          opy(jg)=omy(ig)+pt5*tstep*try
          opz(jg)=omz(ig)+pt5*tstep*trz
        
        enddo

c     store angular velocity at timestep n
        
        omx(ig)=opx(jg)
        omy(ig)=opy(jg)
        omz(ig)=opz(jg)
        
      enddo
      
c     rigid body rotational kinetic energy

      engrot=getkinr(ngrp,idnode,mxnode)

c     propagate chit

      engke=engfke+engtrn
      engtot=engke+engrot
      chitp=2.d0*(engtot-sigma)/qmass
      chitnew=chit+tstep*chitp
      chit0=0.5d0*(chit+chitnew)

c     begin iterations !!-----------------------------------------------

      mxiter=4
      if(ntcons.eq.0) mxiter=mxiter-1
      
      do iter=1,mxiter

c     unconstrained new positions
        
        j=0
        do ifre=ifre1,ifre2

          j=j+1
          i=lstfre(ifre)

c     advance velocity using leapfrog

          uxx(i)=vxo(j)+tstep*(fxx(i)*rmass(i)-chit0*vxx(i))
          uyy(i)=vyo(j)+tstep*(fyy(i)*rmass(i)-chit0*vyy(i))
          uzz(i)=vzo(j)+tstep*(fzz(i)*rmass(i)-chit0*vzz(i))

c     advance positions using leapfrog

          xxx(i)=xxo(j)+tstep*uxx(i)
          yyy(i)=yyo(j)+tstep*uyy(i)
          zzz(i)=zzo(j)+tstep*uzz(i)

        enddo

        if(ntcons.eq.0) safe=.true.
        if(ntcons.gt.0.and.iter.eq.1) then

c     store integrated positions

          j=0
          do ifre=ifre1,ifre2
            
            j=j+1
            i=lstfre(ifre)
            xx1(j)=xxx(i)
            yy1(j)=yyy(i)
            zz1(j)=zzz(i)
            
          enddo

c     global exchange of configuration data
          
          if(mxnode.gt.1)call merge1
     x      (idnode,mxnode,natms,lstme,xxx,yyy,zzz,buffer)

c     apply constraint correction
          
          call rdshake_1
     x      (safe,lshmov,idnode,imcon,mxnode,natms,nscons,
     x      tolnce,tstep,viracc,dxx,dyy,dzz,dxt,dyt,dzt,
     x      txx,tyy,tzz,xxt,yyt,zzt,strcon)

c     contribution to constraint virial 
          
          vircon=vircon+viracc
          do i=1,9
            strcns(i)=strcns(i)+strcon(i)
          enddo

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
        
        if(iter.eq.mxiter)then
          
c     kinetic contribution to stress tensor
          
          call kinstressf(ntfree,idnode,mxnode,strkin)
          
c     restore free atom half step velocity
          
          do ifre=ifre1,ifre2
            
            i=lstfre(ifre)
            vxx(i)=uxx(i)
            vyy(i)=uyy(i)
            vzz(i)=uzz(i)
            
          enddo
          
        endif
        
c     ********: rigid body motion - thermostated  :************
        
c     ***** step 1 : integrate centre of mass motion *********
        
        jg =0
        do ig=igrp1,igrp2

          jg=jg+1
          id=lstgtp(ig)

c     calculate thermostated velocities
          
          rgmas=1.d0/gmass(id)
          uxx(ig)=gvxo(jg)+tstep*(fmx(jg)*rgmas-chit0*gvxx(ig))
          uyy(ig)=gvyo(jg)+tstep*(fmy(jg)*rgmas-chit0*gvyy(ig))
          uzz(ig)=gvzo(jg)+tstep*(fmz(jg)*rgmas-chit0*gvzz(ig))

c     update positions

          gcmx(ig)=gcxo(jg)+tstep*uxx(ig)
          gcmy(ig)=gcyo(jg)+tstep*uyy(ig)
          gcmz(ig)=gczo(jg)+tstep*uzz(ig)

c     calculate full step velocities
          
          gvxx(ig)=0.5d0*(gvxo(jg)+uxx(ig))
          gvyy(ig)=0.5d0*(gvyo(jg)+uyy(ig))
          gvzz(ig)=0.5d0*(gvzo(jg)+uzz(ig))
          
        enddo

c     calculate kinetic energy
        
        engtrn=getkint(ngrp,idnode,mxnode)
        
        if(iter.eq.mxiter)then

          call kinstressg(ngrp,idnode,mxnode,strgrp)
          
c     restore half step velocities 
          
          do ig=igrp1,igrp2
            
            gvxx(ig)=uxx(ig)
            gvyy(ig)=uyy(ig)
            gvzz(ig)=uzz(ig)
            
          enddo
        
        endif
        
c     ****** step 2 : integrate rotational motion **********
        
        safeq=.true.
        engrot=0.d0
        
        jg=0
        do ig=igrp1,igrp2
          
          jg=jg+1
          id=lstgtp(ig)

c     current rotational matrix 
          
          call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)

          trx=(tqx(jg)*rot(1)+tqy(jg)*rot(4)+tqz(jg)*rot(7))*
     x      rotinx(id,2)
     x     +(rotiny(id,1)-rotinz(id,1))*opy(jg)*opz(jg)*
     x      rotinx(id,2)
          try=(tqx(jg)*rot(2)+tqy(jg)*rot(5)+tqz(jg)*rot(8))*
     x      rotiny(id,2)
     x     +(rotinz(id,1)-rotinx(id,1))*opz(jg)*opx(jg)*
     x      rotiny(id,2)
          trz=(tqx(jg)*rot(3)+tqy(jg)*rot(6)+tqz(jg)*rot(9))*
     x      rotinz(id,2)
     x     +(rotinx(id,1)-rotiny(id,1))*opx(jg)*opy(jg)*
     x      rotinz(id,2)

c     correction due to thermostat

          delx=tstep*(trx-chit0*pt5*(omx(ig)+omxo(jg)))
          dely=tstep*(try-chit0*pt5*(omy(ig)+omyo(jg)))
          delz=tstep*(trz-chit0*pt5*(omz(ig)+omzo(jg)))

c     angular velocity at time step n
          
          omx(ig)=omxo(jg)+delx*pt5
          omy(ig)=omyo(jg)+dely*pt5
          omz(ig)=omzo(jg)+delz*pt5

c     angular velocity at time step n+1/2
          
          uxx(ig)=omxo(jg)+delx
          uyy(ig)=omyo(jg)+dely
          uzz(ig)=omzo(jg)+delz

c     angular velocity at time step n+1 (needed for quat algorithm)
          
          oqx(jg)=omxo(jg)+delx*1.5d0
          oqy(jg)=omyo(jg)+dely*1.5d0
          oqz(jg)=omzo(jg)+delz*1.5d0

        enddo

c     calculate rigid body rotational energy
        
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

c     improved prediction of chit 

        engke=engfke+engtrn
        engtot=engke+engrot
        chitp=2.d0*(engtot-sigma)/qmass
        chitnew=chit+tstep*chitp
        chit0=0.5d0*(chit+chitnew)

c     end of thermostat iterations

      enddo

c     assign new quaternions

      call update_quaternions
     x  (safeq,igrp1,igrp2,tstep,quattol,opx,opy,opz,oqx,oqy,oqz)

c     minimum images 
      
      call images(imcon,0,1,natms,cell,xxx,yyy,zzz)
      call images(imcon,idnode,mxnode,ngrp,cell,gcmx,gcmy,gcmz)

c     new atomic positions for atoms in rigid bodies
      
      jr=0
      do ig=igrp1,igrp2

        id=lstgtp(ig)

c     new rotational matrix
        
        call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)
        
        do j=1,numgsit(id)
          
          jr=jr+1
          i=lstrgd(jr)
          
          xxx(i)=rot(1)*gxx(id,j)+rot(2)*gyy(id,j)+rot(3)*gzz(id,j)
     x     +gcmx(ig)
          yyy(i)=rot(4)*gxx(id,j)+rot(5)*gyy(id,j)+rot(6)*gzz(id,j)
     x     +gcmy(ig)
          zzz(i)=rot(7)*gxx(id,j)+rot(8)*gyy(id,j)+rot(9)*gzz(id,j)
     x     +gcmz(ig)

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

c     update thermostat variable

      chit=chitnew

c     conserved quantity less kinetic and potential energy terms
      
      conint=conint+tstep*chit0*qmass/taut**2
      cons1=0.5d0*qmass*chit0**2
      consv=conint+cons1

      if(mxnode.gt.1) then

c     merge new group coordinates and velocities

        call merge(idnode,mxnode,ngrp,mxbuff,gcmx,gcmy,gcmz,buffer)
        call merge(idnode,mxnode,ngrp,mxbuff,gvxx,gvyy,gvzz,buffer)
      
c     merge new atomic coordinates and velocities
      
        call merge1(idnode,mxnode,natms,lstme,xxx,yyy,zzz,buffer)
        call merge1(idnode,mxnode,natms,lstme,vxx,vyy,vzz,buffer)
        if(ntcons.gt.0)
     x    call merge1(idnode,mxnode,natms,lstme,fxx,fyy,fzz,buffer)
      
      endif

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
      deallocate (xx1,yy1,zz1,vxo,vyo,vzo,stat=fail(5))
      deallocate (fmx,fmy,fmz,tqx,tqy,tqz,stat=fail(6))
      deallocate (omxo,omyo,omzo,gcxo,gcyo,gczo,stat=fail(7))
      deallocate (xxt,yyt,zzt,gvxo,gvyo,gvzo,stat=fail(8))
      
      return
      end subroutine nvtq_h1

      subroutine nptq_b1
     x  (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x  ntcons,ntfree,ntpatm,elrc,engke,engrot,virlrc,press,
     x  quattol,sigma,taup,taut,tolnce,tstep,virtot,vircom,
     x  vircon,volm)
      
c***********************************************************************
c     
c     dlpoly routine to integrate rigid body equations of motion
c     using implicit leapfrog quaternion algorithm
c     plus bond constraints - PROVIDED rigid body sites
c     and constraint sites do not coincide.
c     
c     dl_poly subroutine for integrating newtonian equations of
c     motion in molecular dynamics - verlet leapfrog with Berendsen
c     thermostat and barostat.
c     isothermal compressibility (beta) set to that of liquid water
c     = 0.007372 dlpoly units
c     
c     parallel replicated data version : block data
c     
c     for systems using bond constraints
c     
c     parallel replicated data version : block data
c     
c     tqx,tqy,tqz = torque in lab fixed frame (input)
c     omx,omy,omz = angular velocity in body fixed frame (principl axis)
c     rotinx,y,z  = rotational inertia in body fixed frame
c     
c     copyright daresbury laboratory 1993
c     author      t.forester october 1993
c     amended     t.forester dec 1994 : block data
c     amended     w.smith nov 2005
c     
c**********************************************************************
      
      implicit none

      integer, parameter :: nnn=16
      real(8), parameter :: pt5=0.5d0

      logical safe,lshmov,newjob,safeq
      integer imcon,idnode,mxnode,natms,ngrp,nscons,ntcons,ntfree
      integer ntpatm,i,fail,igrp,igrp1,igrp2,idum,ifre,ifre1
      integer ifre2,j,jg,ig,jr,k,id,iter,mxiter
      real(8) elrc,engke,engrot,virlrc,press,quattol,sigma,taup,taut
      real(8) tolnce,tstep,virtot,vircom,vircon,volm,rot,engfke,uni
      real(8) cell0,beta,volm0,elrc0,virlrc0,rstep,rtsq
      real(8) engtrn,trx,try,trz,chip0,scale,engtot,chit0,viracc,czero
      real(8) rgmas,vaa,vbb,vcc,strkin,strcon,strgrp,psyst
      
      real(8), allocatable :: dens0(:)
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
      real(8), allocatable :: fmx(:),fmy(:),fmz(:)
      real(8), allocatable :: tqx(:),tqy(:),tqz(:)
      real(8), allocatable :: omxo(:),omyo(:),omzo(:)
      real(8), allocatable :: gvxo(:),gvyo(:),gvzo(:)
      real(8), allocatable :: gcxo(:),gcyo(:),gczo(:)

      dimension fail(nnn),rot(9),cell0(9),czero(9),uni(9)
      dimension strcon(9),strgrp(9),strkin(9)

      save newjob,volm0,elrc0,virlrc0,cell0,dens0,igrp1,igrp2
      save ifre1,ifre2
      
      data newjob/.true./,beta/7.3728d-3/
      data uni/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/
      
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
      allocate (fmx(msgrp),fmy(msgrp),fmz(msgrp),stat=fail(12))
      allocate (tqx(msgrp),tqy(msgrp),tqz(msgrp),stat=fail(13))
      allocate (omxo(msgrp),omyo(msgrp),omzo(msgrp),stat=fail(14))
      allocate (gvxo(msgrp),gvyo(msgrp),gvzo(msgrp),stat=fail(15))
      allocate (gcxo(msgrp),gcyo(msgrp),gczo(msgrp),stat=fail(16))
      do i=1,nnn
        if(fail(i).ne.0)call error(idnode,1530)
      enddo
     
      if(newjob) then

c     store initial values of volume and long range corrections

        volm0=volm
        elrc0=elrc
        virlrc0=virlrc
        fail(1)=0
        allocate (dens0(mxatyp),stat=fail(1))
        if(fail(1).ne.0)call error(idnode,1540)
        do i=1,ntpatm
          dens0(i)=dens(i)
        enddo
        do i=1,9
          cell0(i)=cell(i)
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

c     constraint stress tensor accumulators

      vircon=0.d0
      do i=1,9
        strcns(i)=0.d0
      enddo
      
c     current cell vectors
      
      do i=1,9
        czero(i)=cell(i)
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
        omxo(jg)=omx(ig)
        omyo(jg)=omy(ig)
        omzo(jg)=omz(ig)
        gcxo(jg)=gcmx(ig)
        gcyo(jg)=gcmy(ig)
        gczo(jg)=gcmz(ig)
        gvxo(jg)=gvxx(ig)
        gvyo(jg)=gvyy(ig)
        gvzo(jg)=gvzz(ig)

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

c     periodic boundary condition for displacements
      
      call images(imcon,0,1,jr,cell,dtx,dty,dtz)

c     estimate velocity and temperature at half-time step
      
      j=0
      do ifre=ifre1,ifre2
        
        j=j+1
        i=lstfre(ifre)

        vxx(i)=vxo(j)+pt5*tstep*rmass(i)*fxx(i)
        vyy(i)=vyo(j)+pt5*tstep*rmass(i)*fyy(i)
        vzz(i)=vzo(j)+pt5*tstep*rmass(i)*fzz(i)

      enddo

c     kinetic energy of free atoms

      engfke=getkinf(ntfree,idnode,mxnode)

c     estimate kinetic energy of rigid body com
      
      jg=0
      jr=0
      do ig=igrp1,igrp2
        
        jg=jg+1
        id=lstgtp(ig)

c     calculate net force on rigid body
        
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

c     centre of mass velocities at half-step
        
        gvxx(ig)=gvxo(jg)+pt5*tstep*fmx(jg)/gmass(id)
        gvyy(ig)=gvyo(jg)+pt5*tstep*fmy(jg)/gmass(id)
        gvzz(ig)=gvzo(jg)+pt5*tstep*fmz(jg)/gmass(id)

      enddo

c     translation kinetic energy of rigid bodies

      engtrn=getkint(ngrp,idnode,mxnode)
      
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

        opx(jg)=omx(ig)
        opy(jg)=omy(ig)
        opz(jg)=omz(ig)

c       iterate angular velocity for time step n (e. yezdimer)

        do i=1,5
        
          trx=(tqx(jg)*rot(1)+tqy(jg)*rot(4)+tqz(jg)*rot(7))*
     x      rotinx(id,2)
     x     +(rotiny(id,1)-rotinz(id,1))*opy(jg)*opz(jg)*
     x      rotinx(id,2)
          try=(tqx(jg)*rot(2)+tqy(jg)*rot(5)+tqz(jg)*rot(8))*
     x      rotiny(id,2)
     x     +(rotinz(id,1)-rotinx(id,1))*opz(jg)*opx(jg)*
     x      rotiny(id,2)
          trz=(tqx(jg)*rot(3)+tqy(jg)*rot(6)+tqz(jg)*rot(9))*
     x      rotinz(id,2)
     x     +(rotinx(id,1)-rotiny(id,1))*opx(jg)*opy(jg)*
     x      rotinz(id,2)

c     improved angular velocity at time step n

          opx(jg)=omx(ig)+pt5*tstep*trx
          opy(jg)=omy(ig)+pt5*tstep*try
          opz(jg)=omz(ig)+pt5*tstep*trz

        enddo

c     store angular velocity at timestep n
        
        omx(ig)=opx(jg)
        omy(ig)=opy(jg)
        omz(ig)=opz(jg)
        
      enddo
      
c     rigid body rotational kinetic energy

      engrot=getkinr(ngrp,idnode,mxnode)

c     pressure control variable - taup is the relaxation time

      engke=engfke+engtrn
      psyst=(2.d0*engke-virtot-vircon-vircom)/(3.d0*volm)
      chip0=1.d0+beta*tstep*(psyst-press)/taup
      scale=chip0**(1.d0/3.d0)
      
c     temperature scaling  coefficient - taut is the relaxation time
      
      engtot=engke+engrot
      chit0= sqrt(1.d0+tstep/taut*(sigma/engtot-1.d0))

c     begin iterations !!-----------------------------------------------

      mxiter=5
      if(ntcons.eq.0) mxiter=mxiter-1
      
      do iter=1,mxiter

c     unconstrained new positions
        
        j=0
        do ifre=ifre1,ifre2

          j=j+1
          i=lstfre(ifre)

c     advance velocity using leapfrog

          uxx(i)=(vxo(j)+tstep*rmass(i)*fxx(i))*chit0
          uyy(i)=(vyo(j)+tstep*rmass(i)*fyy(i))*chit0
          uzz(i)=(vzo(j)+tstep*rmass(i)*fzz(i))*chit0

c     update positions
          
          xxx(i)=tstep*uxx(i)+scale*xxo(j)
          yyy(i)=tstep*uyy(i)+scale*yyo(j)
          zzz(i)=tstep*uzz(i)+scale*zzo(j)

        enddo

c     estimate new cell tensor

        do i=1,9
          cell(i)=scale*czero(i)
        enddo

        if(ntcons.eq.0) safe=.true.
        if(ntcons.gt.0) then

c     store integrated positions

          j=0
          do ifre=ifre1,ifre2

            j=j+1
            i=lstfre(ifre)
            xx1(j)=xxx(i)
            yy1(j)=yyy(i)
            zz1(j)=zzz(i)

          enddo

c     global exchange of configuration data
          
          if(mxnode.gt.1)call merge1
     x      (idnode,mxnode,natms,lstme,xxx,yyy,zzz,buffer)

c     apply constraint correction
          
          call rdshake_1
     x      (safe,lshmov,idnode,imcon,mxnode,natms,nscons,
     x      tolnce,tstep,viracc,dxx,dyy,dzz,dxt,dyt,dzt,
     x      txx,tyy,tzz,xxt,yyt,zzt,strcon)

c     contribution to constraint virial 
          
          vircon=vircon+viracc
          do i=1,9
            strcns(i)=strcns(i)+strcon(i)
          enddo

c     calculate force correction
          
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
        
        if(iter.eq.mxiter)call kinstressf(ntfree,idnode,mxnode,strkin)
        
c     restore free atom half step velocity
        
        do ifre=ifre1,ifre2
          
          i=lstfre(ifre)
          vxx(i)=uxx(i)
          vyy(i)=uyy(i)
          vzz(i)=uzz(i)
          
        enddo
        
c     ********: rigid body motion - thermostated  :************

c     ***** step 1 : integrate centre of mass motion *********
        
        jg =0
        do ig=igrp1,igrp2

          jg=jg+1
          id=lstgtp(ig)

c     calculate thermostated velocities
          
          rgmas=1.d0/gmass(id)
          uxx(ig)=(gvxo(jg)+tstep*(fmx(jg)*rgmas))*chit0
          uyy(ig)=(gvyo(jg)+tstep*(fmy(jg)*rgmas))*chit0
          uzz(ig)=(gvzo(jg)+tstep*(fmz(jg)*rgmas))*chit0

c     update positions : 
          
          gcmx(ig)=scale*gcxo(jg)+tstep*uxx(ig)
          gcmy(ig)=scale*gcyo(jg)+tstep*uyy(ig)
          gcmz(ig)=scale*gczo(jg)+tstep*uzz(ig)
          
c     calculate full step velocities
          
          gvxx(ig)=0.5d0*(gvxo(jg)+uxx(ig))
          gvyy(ig)=0.5d0*(gvyo(jg)+uyy(ig))
          gvzz(ig)=0.5d0*(gvzo(jg)+uzz(ig))
          
        enddo
        
c     calculate kinetic energy
        
        engtrn=getkint(ngrp,idnode,mxnode)
        
        if(iter.eq.mxiter)call kinstressg(ngrp,idnode,mxnode,strgrp)
        
c     restore half step velocities 

        do ig=igrp1,igrp2

          gvxx(ig)=uxx(ig)
          gvyy(ig)=uyy(ig)
          gvzz(ig)=uzz(ig)
          
        enddo
        
c     ****** step 2 : integrate rotational motion **********
        
        jg=0
        safeq=.true.
        do ig=igrp1,igrp2
          
          jg=jg+1
          id=lstgtp(ig)

c     current rotational matrix 
          
          call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)
        
          trx=(tqx(jg)*rot(1)+tqy(jg)*rot(4)+tqz(jg)*rot(7))*
     x      rotinx(id,2)
     x     +(rotiny(id,1)-rotinz(id,1))*opy(jg)*opz(jg)*
     x      rotinx(id,2)
          try=(tqx(jg)*rot(2)+tqy(jg)*rot(5)+tqz(jg)*rot(8))*
     x      rotiny(id,2)
     x     +(rotinz(id,1)-rotinx(id,1))*opz(jg)*opx(jg)*
     x      rotiny(id,2)
          trz=(tqx(jg)*rot(3)+tqy(jg)*rot(6)+tqz(jg)*rot(9))*
     x      rotinz(id,2)
     x     +(rotinx(id,1)-rotiny(id,1))*opx(jg)*opy(jg)*
     x      rotinz(id,2)

c     scaled angular velocity at time step n
          
          omx(ig)=(omxo(jg)+pt5*tstep*trx)*chit0
          omy(ig)=(omyo(jg)+pt5*tstep*try)*chit0
          omz(ig)=(omzo(jg)+pt5*tstep*trz)*chit0

c     angular velocity at time step n+1/2
          
          uxx(ig)=(omxo(jg)+tstep*trx)*chit0
          uyy(ig)=(omyo(jg)+tstep*try)*chit0
          uzz(ig)=(omzo(jg)+tstep*trz)*chit0

c     angular velocity at time step n+1 (needed for quat algorithm)
          
          oqx(jg)=(omxo(jg)+1.5d0*tstep*trx)*chit0
          oqy(jg)=(omyo(jg)+1.5d0*tstep*try)*chit0
          oqz(jg)=(omzo(jg)+1.5d0*tstep*trz)*chit0

        enddo
        
c     calculate rigid body rotational energy
        
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
        
c     total translational kinetic energy
      
        engke=engtrn+engfke
      
c     total kinetic energy
        
        engtot=engke+engrot
        
c     improved pressure control variable

        psyst=(2.d0*engke-virtot-vircon-vircom)/(3.d0*volm)
        chip0=1.d0+beta*tstep*(psyst-press)/taup
        scale=chip0**(1.d0/3.d0)
        
c     improved temperature scaling  coefficient
        
        chit0=sqrt(1.d0+tstep/taut*(sigma/engtot-1.d0))
        
c     end of thermostat/barostat iterations

      enddo

c     scale cell vectors

      scale=((chip0*volm)/volm0)**(1.d0/3.d0)

      do i=1,9
        cell(i)=scale*cell0(i)
      enddo

c     construct scaling tensor (for later!)

      do i=1,9
        eta(i)=scale*uni(i)
      enddo

c     assign new quaternions

      call update_quaternions
     x  (safeq,igrp1,igrp2,tstep,quattol,opx,opy,opz,oqx,oqy,oqz)

c     minimum images 
      
      call images(imcon,0,1,natms,cell,xxx,yyy,zzz)
      call images(imcon,idnode,mxnode,ngrp,cell,gcmx,gcmy,gcmz)

c     new atomic positions for atoms in rigid bodies
      
      jr=0
      do ig=igrp1,igrp2

        id=lstgtp(ig)

c     new rotational matrix
        
        call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)
        
        do j=1,numgsit(id)
          
          jr=jr+1
          i=lstrgd(jr)
          
          xxx(i)=rot(1)*gxx(id,j)+rot(2)*gyy(id,j)+rot(3)*gzz(id,j)
     x     +gcmx(ig)
          yyy(i)=rot(4)*gxx(id,j)+rot(5)*gyy(id,j)+rot(6)*gzz(id,j)
     x     +gcmy(ig)
          zzz(i)=rot(7)*gxx(id,j)+rot(8)*gyy(id,j)+rot(9)*gzz(id,j)
     x     +gcmz(ig)

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
        if(ntcons.gt.0)
     x    call merge1(idnode,mxnode,natms,lstme,fxx,fyy,fzz,buffer)
      
      endif

c     adjust long range corrections and number density
      
      elrc=elrc0*(volm0/volm)
      virlrc=virlrc0*(volm0/volm)
      
      do k=1,ntpatm
        dens(k)=dens0(k)*(volm0/volm)
      enddo

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
      deallocate (xx1,yy1,zz1,vxo,vyo,vzo,stat=fail(5))
      deallocate (fmx,fmy,fmz,tqx,tqy,tqz,stat=fail(6))
      deallocate (omxo,omyo,omzo,gcxo,gcyo,gczo,stat=fail(7))
      deallocate (xxt,yyt,zzt,gvxo,gvyo,gvzo,stat=fail(8))
      
      return
      end subroutine nptq_b1

      subroutine nptq_h1
     x (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x  ntcons,ntfree,ntpatm,chip,chit,consv,conint,elrc,engke,
     x  engrot,virlrc,press,quattol,sigma,taup,taut,temp,tolnce,
     x  tstep,virtot,vircom,vircon,volm)
      
c***********************************************************************
c     
c     dlpoly routine to integrate rigid body equations of motion
c     using implicit leapfrog quaternion algorithm
c     plus bond constraints - PROVIDED rigid body sites
c     and constraint sites do not coincide.
c     
c     dl_poly subroutine for integrating newtonian equations of
c     motion in molecular dynamics - verlet leapfrog with Hoover
c     thermostat and barostat (Melchionna et al variant)
c     
c     parallel replicated data version : block data
c     
c     for systems using bond constraints
c     
c     parallel replicated data version : block data
c     
c     tqx,tqy,tqz = torque in lab fixed frame (input)
c     omx,omy,omz = angular velocity in body fixed frame (principl axis)
c     rotinx,y,z  = rotational inertia in body fixed frame
c     
c     copyright daresbury laboratory 1993
c     author      t.forester october 1993
c     amended     t.forester dec 1994 : block data
c     amended     w.smith sep 1999 : euler equation
c     
c**********************************************************************
      
      implicit none
      
      integer, parameter :: nnn=16
      real(8), parameter :: pt5=0.5d0

      logical safe,safeq,lshmov,newjob
      integer imcon,idnode,mxnode,natms,ngrp,nscons,ntcons,ntfree
      integer ntpatm,fail,i,igrp,igrp1,igrp2,ifre,ifre1,ifre2
      integer j,k,ig,jg,jr,id,iter,mxiter,idum
      real(8) chip,chit,consv,conint,elrc,engke,engrot,virlrc,press
      real(8) quattol,sigma,taup,taut,temp,tolnce,tstep,virtot,vircom
      real(8) vircon,volm,cell0,rot,volm0,elrc0,rtsq,uni
      real(8) virlrc0,strkin,rstep,qmass,pmass,strgrp,strcon
      real(8) trx,try,trz,chipp,chipnew,chip0,engtot,chitp
      real(8) chitnew,chit0,volnew,scale,viracc,rgmas
      real(8) vaa,vbb,vcc,vold,cons1,cons2,cons3,delx,dely,delz
      real(8) engtrn,totmas,com,vom,engfke

      real(8), allocatable :: dens0(:)
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
      real(8), allocatable :: fmx(:),fmy(:),fmz(:)
      real(8), allocatable :: tqx(:),tqy(:),tqz(:)
      real(8), allocatable :: omxo(:),omyo(:),omzo(:)
      real(8), allocatable :: gvxo(:),gvyo(:),gvzo(:)
      real(8), allocatable :: gcxo(:),gcyo(:),gczo(:)

      dimension fail(nnn),rot(9),cell0(9),uni(9)
      dimension strkin(9),strcon(9),strgrp(9),com(3),vom(3)

      save newjob,volm0,elrc0,virlrc0,cell0,dens0,pmass,qmass
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
      allocate (xxt(mxatms),yyt(mxatms),zzt(mxatms),stat=fail(6))
      allocate (txx(mxatms),tyy(mxatms),tzz(mxatms),stat=fail(7))
      allocate (dxt(mxcons),dyt(mxcons),dzt(mxcons),stat=fail(8))
      allocate (xxo(msatms),yyo(msatms),zzo(msatms),stat=fail(9))
      allocate (xx1(msatms),yy1(msatms),zz1(msatms),stat=fail(10))
      allocate (vxo(msatms),vyo(msatms),vzo(msatms),stat=fail(11))
      allocate (fmx(msgrp),fmy(msgrp),fmz(msgrp),stat=fail(12))
      allocate (tqx(msgrp),tqy(msgrp),tqz(msgrp),stat=fail(13))
      allocate (omxo(msgrp),omyo(msgrp),omzo(msgrp),stat=fail(14))
      allocate (gvxo(msgrp),gvyo(msgrp),gvzo(msgrp),stat=fail(15))
      allocate (gcxo(msgrp),gcyo(msgrp),gczo(msgrp),stat=fail(16))
      do i=1,nnn
        if(fail(i).ne.0)call error(idnode,1550)
      enddo

      if(newjob) then

c     store initial values of volume and long range corrections

        volm0=volm
        elrc0=elrc
        virlrc0=virlrc
        fail(1)=0
        allocate (dens0(mxatyp),stat=fail(1))
        if(fail(1).ne.0)call error(idnode,1560)

        do i=1,ntpatm
          dens0(i)=dens(i)
        enddo

        do i=1,9
          cell0(i)=cell(i)
        enddo

c     inertia parameter for Nose-Hoover thermostat
        
        qmass=2.0d0*sigma*taut**2
        pmass=2.0d0*sigma*taup**2
        
c     total system mass
        
        totmas=getmass(natms,idnode,mxnode)
        
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

c     temporary stress tensor accumulators

      vircon=0.d0
      do i=1,9
        strcns(i)=0.d0
      enddo

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
        omxo(jg)=omx(ig)
        omyo(jg)=omy(ig)
        omzo(jg)=omz(ig)
        gcxo(jg)=gcmx(ig)
        gcyo(jg)=gcmy(ig)
        gczo(jg)=gcmz(ig)
        gvxo(jg)=gvxx(ig)
        gvyo(jg)=gvyy(ig)
        gvzo(jg)=gvzz(ig)

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

c     periodic boundary condition for displacements
      
      call images(imcon,0,1,jr,cell,dtx,dty,dtz)

c     calculate centre of mass

      call getcom(natms,idnode,mxnode,totmas,com)
      
c     estimate velocity at full step
      
      j=0
      do ifre=ifre1,ifre2
        
        j=j+1
        i=lstfre(ifre)

        vxx(i)=vxo(j)+pt5*tstep*rmass(i)*fxx(i)
        vyy(i)=vyo(j)+pt5*tstep*rmass(i)*fyy(i)
        vzz(i)=vzo(j)+pt5*tstep*rmass(i)*fzz(i)

      enddo

c     kinetic energy of free atoms

      engfke=getkinf(ntfree,idnode,mxnode)

c     estimate kinetic energy of rigid body com
      
      jg=0
      jr=0
      do ig=igrp1,igrp2
        
        jg=jg+1
        id=lstgtp(ig)

c     calculate net force on rigid body
        
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

c     centre of mass velocities at half-step
        
        gvxx(ig)=gvxo(jg)+pt5*tstep*fmx(jg)/gmass(id)
        gvyy(ig)=gvyo(jg)+pt5*tstep*fmy(jg)/gmass(id)
        gvzz(ig)=gvzo(jg)+pt5*tstep*fmz(jg)/gmass(id)

      enddo

c     translation kinetic energy of rigid bodies

      engtrn=getkint(ngrp,idnode,mxnode)
      
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

        opx(jg)=omx(ig)
        opy(jg)=omy(ig)
        opz(jg)=omz(ig)

c       iterate angular velocity for time step n (e. yezdimer)

        do i=1,5
        
          trx=(tqx(jg)*rot(1)+tqy(jg)*rot(4)+tqz(jg)*rot(7))*
     x      rotinx(id,2)
     x     +(rotiny(id,1)-rotinz(id,1))*opy(jg)*opz(jg)*
     x      rotinx(id,2)
          try=(tqx(jg)*rot(2)+tqy(jg)*rot(5)+tqz(jg)*rot(8))*
     x      rotiny(id,2)
     x     +(rotinz(id,1)-rotinx(id,1))*opz(jg)*opx(jg)*
     x      rotiny(id,2)
          trz=(tqx(jg)*rot(3)+tqy(jg)*rot(6)+tqz(jg)*rot(9))*
     x      rotinz(id,2)
     x     +(rotinx(id,1)-rotiny(id,1))*opx(jg)*opy(jg)*
     x      rotinz(id,2)

c     improved angular velocity at time step n

          opx(jg)=omx(ig)+pt5*tstep*trx
          opy(jg)=omy(ig)+pt5*tstep*try
          opz(jg)=omz(ig)+pt5*tstep*trz
          
        enddo

c     store angular velocity at timestep n
        
        omx(ig)=opx(jg)
        omy(ig)=opy(jg)
        omz(ig)=opz(jg)
        
      enddo
      
c     rigid body rotational kinetic energy

      engrot=getkinr(ngrp,idnode,mxnode)

c     propagate chip

      engke=engfke+engtrn
      chipp=(2.d0*engke-virtot-vircon-vircom-3.d0*press*volm)/pmass-
     x  chit*chip
      chipnew=chip+tstep*chipp
      chip0=0.5d0*(chip+chipnew)

c     propagate chit

      engtot=engke+engrot
      chitp=(2.d0*(engtot-sigma)+pmass*chip0**2-boltz*temp)/qmass
      chitnew=chit+tstep*chitp
      chit0=0.5d0*(chit+chitnew)

c     begin iterations !!-----------------------------------------------

      mxiter=5
      if(ntcons.eq.0) mxiter=mxiter-1
      
      do iter=1,mxiter

c     unconstrained new positions
        
        j=0
        do ifre=ifre1,ifre2

          j=j+1
          i=lstfre(ifre)

c     advance velocity using leapfrog

          uxx(i)=vxo(j)+tstep*(fxx(i)*rmass(i)-(chit0+chip0)*vxx(i))
          uyy(i)=vyo(j)+tstep*(fyy(i)*rmass(i)-(chit0+chip0)*vyy(i))
          uzz(i)=vzo(j)+tstep*(fzz(i)*rmass(i)-(chit0+chip0)*vzz(i))

c     advance positions using leapfrog

          xxx(i)=xxo(j)+tstep*(uxx(i)+
     x      chipnew*((xxx(i)+xxo(j))*0.5d0-com(1)))
          yyy(i)=yyo(j)+tstep*(uyy(i)+
     x      chipnew*((yyy(i)+yyo(j))*0.5d0-com(2)))
          zzz(i)=zzo(j)+tstep*(uzz(i)+
     x      chipnew*((zzz(i)+zzo(j))*0.5d0-com(3)))

        enddo

c     estimate new cell parameters

        volnew=volm*exp(3.d0*tstep*chipnew)
        scale=(volnew/volm0)**(1.d0/3.d0)
        do i=1,9
          cell(i)=cell0(i)*scale
        enddo
        
        if(ntcons.eq.0) safe=.true.
        if(ntcons.gt.0) then

c     store integrated positions

          j=0
          do ifre=ifre1,ifre2

            j=j+1
            i=lstfre(ifre)
            xx1(j)=xxx(i)
            yy1(j)=yyy(i)
            zz1(j)=zzz(i)

          enddo

c     global exchange of configuration data
          
          if(mxnode.gt.1)call merge1
     x      (idnode,mxnode,natms,lstme,xxx,yyy,zzz,buffer)

c     apply constraint correction
          
          call rdshake_1
     x      (safe,lshmov,idnode,imcon,mxnode,natms,nscons,
     x      tolnce,tstep,viracc,dxx,dyy,dzz,dxt,dyt,dzt,
     x      txx,tyy,tzz,xxt,yyt,zzt,strcon)

c     contribution to constraint virial 
          
          vircon=vircon+viracc
          do i=1,9
            strcns(i)=strcns(i)+strcon(i)
          enddo
          
c     calculate force correction
          
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
        
        if(iter.eq.mxiter)then
          
          call kinstressf(ntfree,idnode,mxnode,strkin)
          
c     restore free atom half step velocity
          
          do ifre=ifre1,ifre2
            
            i=lstfre(ifre)
            vxx(i)=uxx(i)
            vyy(i)=uyy(i)
            vzz(i)=uzz(i)
            
          enddo
        
        endif
        
c     ********: rigid body motion - thermostated  :***********

c     ***** step 1 : integrate centre of mass motion *********
        
        jg =0
        do ig=igrp1,igrp2

          jg=jg+1
          id=lstgtp(ig)

c     calculate thermostated velocities
          
          rgmas=1.d0/gmass(id)
          uxx(ig)=gvxo(jg)+tstep*(fmx(jg)*rgmas-(chit0+chip0)*
     x      gvxx(ig))
          uyy(ig)=gvyo(jg)+tstep*(fmy(jg)*rgmas-(chit0+chip0)*
     x      gvyy(ig))
          uzz(ig)=gvzo(jg)+tstep*(fmz(jg)*rgmas-(chit0+chip0)*
     x      gvzz(ig))

c     advance positions using leapfrog

          gcmx(ig)=gcxo(jg)+tstep*(uxx(ig)+
     x      chipnew*((gcmx(ig)+gcxo(jg))*0.5d0-com(1)))
          gcmy(ig)=gcyo(jg)+tstep*(uyy(ig)+
     x      chipnew*((gcmy(ig)+gcyo(jg))*0.5d0-com(2)))
          gcmz(ig)=gczo(jg)+tstep*(uzz(ig)+
     x      chipnew*((gcmz(ig)+gczo(jg))*0.5d0-com(3)))
          
c     calculate full step velocities
          
          gvxx(ig)=0.5d0*(gvxo(jg)+uxx(ig))
          gvyy(ig)=0.5d0*(gvyo(jg)+uyy(ig))
          gvzz(ig)=0.5d0*(gvzo(jg)+uzz(ig))
          
        enddo

c     calculate kinetic energy and stress tensor
        
        engtrn=getkint(ngrp,idnode,mxnode)
        
        if(iter.eq.mxiter)call kinstressg(ngrp,idnode,mxnode,strgrp)
        
c     restore half step velocities 

        do ig=igrp1,igrp2

          gvxx(ig)=uxx(ig)
          gvyy(ig)=uyy(ig)
          gvzz(ig)=uzz(ig)
          
        enddo
        
c     ****** step 2 : integrate rotational motion **********
        
        jg=0
        safeq=.true.
        do ig=igrp1,igrp2
          
          jg=jg+1
          id=lstgtp(ig)

c     current rotational matrix 
          
          call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)

          trx=(tqx(jg)*rot(1)+tqy(jg)*rot(4)+tqz(jg)*rot(7))*
     x      rotinx(id,2)
     x     +(rotiny(id,1)-rotinz(id,1))*opy(jg)*opz(jg)*
     x      rotinx(id,2)
          try=(tqx(jg)*rot(2)+tqy(jg)*rot(5)+tqz(jg)*rot(8))*
     x      rotiny(id,2)
     x     +(rotinz(id,1)-rotinx(id,1))*opz(jg)*opx(jg)*
     x      rotiny(id,2)
          trz=(tqx(jg)*rot(3)+tqy(jg)*rot(6)+tqz(jg)*rot(9))*
     x      rotinz(id,2)
     x     +(rotinx(id,1)-rotiny(id,1))*opx(jg)*opy(jg)*
     x      rotinz(id,2)
          
c     correction due to thermostat

          delx=tstep*(trx-chit0*pt5*(omx(ig)+omxo(jg)))
          dely=tstep*(try-chit0*pt5*(omy(ig)+omyo(jg)))
          delz=tstep*(trz-chit0*pt5*(omz(ig)+omzo(jg)))

c     angular velocity at time step n
          
          omx(ig)=omxo(jg)+delx*pt5
          omy(ig)=omyo(jg)+dely*pt5
          omz(ig)=omzo(jg)+delz*pt5

c     angular velocity at time step n+1/2
          
          uxx(ig)=omxo(jg)+delx
          uyy(ig)=omyo(jg)+dely
          uzz(ig)=omzo(jg)+delz

c     angular velocity at time step n+1 (needed for quat algorithm)
          
          oqx(jg)=omxo(jg)+delx*1.5d0
          oqy(jg)=omyo(jg)+dely*1.5d0
          oqz(jg)=omzo(jg)+delz*1.5d0

        enddo

c     calculate rigid body rotational energy
        
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
        
c     improved prediction of chip

        engke=engfke+engtrn
        chipp=(2.d0*engke-virtot-vircom-vircon-3.d0*press*volm)/pmass-
     x    chit0*chip0
        chipnew=chip+tstep*chipp
        chip0=0.5d0*(chip+chipnew)

c     improved prediction of chit
        
        engtot=engke+engrot
        chitp=(2.d0*(engtot-sigma)+pmass*chip0**2-boltz*temp)/qmass
        chitnew=chit+tstep*chitp
        chit0=0.5d0*(chit+chitnew)

c     end of thermostat iterations

      enddo

c     assign new quaternions

      call update_quaternions
     x  (safeq,igrp1,igrp2,tstep,quattol,opx,opy,opz,oqx,oqy,oqz)

c     minimum images
      
      call images(imcon,0,1,natms,cell,xxx,yyy,zzz)
      call images(imcon,idnode,mxnode,ngrp,cell,gcmx,gcmy,gcmz)

c     new atomic positions for atoms in rigid bodies
      
      jr=0
      do ig=igrp1,igrp2

        id=lstgtp(ig)

c     new rotational matrix
        
        call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)
        
        do j=1,numgsit(id)
          
          jr=jr+1
          i=lstrgd(jr)
          
          xxx(i)=rot(1)*gxx(id,j)+rot(2)*gyy(id,j)+rot(3)*gzz(id,j)
     x     +gcmx(ig)
          yyy(i)=rot(4)*gxx(id,j)+rot(5)*gyy(id,j)+rot(6)*gzz(id,j)
     x     +gcmy(ig)
          zzz(i)=rot(7)*gxx(id,j)+rot(8)*gyy(id,j)+rot(9)*gzz(id,j)
     x     +gcmz(ig)

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

c     update volume

      vold=volm
      volm=volm*exp(3.d0*tstep*chipnew)

c     scale cell vectors - isotropic

      scale=(volm/volm0)**(1.d0/3.d0)
      do i=1,9
        cell(i)=cell0(i)*scale
      enddo

c     construct scaling tensor (for later!)

      do i=1,9
        eta(i)=chipnew*uni(i)
      enddo

c     update thermostat and barostat variables

      chit=chitnew
      chip=chipnew

c     adjust long range corrections and number density
      
      elrc=elrc0*(volm0/volm)
      virlrc=virlrc0*(volm0/volm)
      
      do k=1,ntpatm
        dens(k)=dens0(k)*(volm0/volm)
      enddo

c     conserved quantity less kinetic and potential energy terms
      
      conint=conint+tstep*chit0*(qmass/taut**2+boltz*temp)
      cons1=0.5d0*qmass*chit0**2
      cons2=press*vold
      cons3=0.5d0*pmass*chip0**2
      consv=conint+cons1+cons2+cons3

      if(mxnode.gt.1) then

c     merge new group coordinates and velocities

        call merge(idnode,mxnode,ngrp,mxbuff,gcmx,gcmy,gcmz,buffer)
        call merge(idnode,mxnode,ngrp,mxbuff,gvxx,gvyy,gvzz,buffer)
      
c     merge new atomic coordinates and velocities
      
        call merge1(idnode,mxnode,natms,lstme,xxx,yyy,zzz,buffer)
        call merge1(idnode,mxnode,natms,lstme,vxx,vyy,vzz,buffer)
        if(ntcons.gt.0)
     x    call merge1(idnode,mxnode,natms,lstme,fxx,fyy,fzz,buffer)
      
      endif

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
      deallocate (xx1,yy1,zz1,vxo,vyo,vzo,stat=fail(5))
      deallocate (fmx,fmy,fmz,tqx,tqy,tqz,stat=fail(6))
      deallocate (omxo,omyo,omzo,gcxo,gcyo,gczo,stat=fail(7))
      deallocate (xxt,yyt,zzt,gvxo,gvyo,gvzo,stat=fail(8))
      
      return
      end subroutine nptq_h1

      subroutine nstq_b1
     x  (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x  ntcons,ntfree,ntpatm,mode,elrc,engke,engrot,virlrc,press,
     x  quattol,sigma,taup,taut,tolnce,tstep,vircom,vircon,volm)
      
c***********************************************************************
c     
c     dlpoly routine to integrate rigid body equations of motion
c     using implicit leapfrog quaternion algorithm
c     plus bond constraints - PROVIDED rigid body sites
c     and constraint sites do not coincide.
c     
c     dl_poly subroutine for integrating newtonian equations of
c     motion in molecular dynamics - verlet leapfrog with Berendsen
c     thermostat and barostat. (cell may change shape).
c     isothermal compressibility (beta) set to that of liquid water
c     = 0.007372 dlpoly units
c     
c     parallel replicated data version : block data
c     
c     for systems using bond constraints
c     
c     parallel replicated data version : block data
c     
c     tqx,tqy,tqz = torque in lab fixed frame (input)
c     omx,omy,omz = angular velocity in body fixed frame (principl axis)
c     rotinx,y,z  = rotational inertia in body fixed frame
c     
c     copyright daresbury laboratory 1993
c     author      t.forester october 1993
c     amended     t.forester dec 1994 : block data
c     amended     w.smith sep 1999 : euler equation
c     
c**********************************************************************
      
      implicit none
      
      integer, parameter :: nnn=16
      real(8), parameter :: pt5=0.5d0
      
      logical safe,lshmov,newjob,safeq
      integer imcon,idnode,mxnode,natms,ngrp,nscons,ntcons,ntfree
      integer ntpatm,fail,i,idum,igrp,igrp1,igrp2,ifre,ifre1,mode
      integer ifre2,j,k,ig,jg,jr,id,iter,mxiter
      real(8) elrc,engke,engrot,virlrc,press,quattol,sigma,taup,taut
      real(8) tolnce,tstep,vircom,vircon,volm,beta,uni,cell0
      real(8) volm0,elrc0,virlrc0,rot,rstep,rtsq,engfke
      real(8) engtrn,trx,try,trz,engtot,chit0,rgmas
      real(8) vaa,vbb,vcc,viracc,strkin,strcon,strgrp,stres0

      real(8), allocatable :: dens0(:)
      real(8), allocatable :: opx(:),opy(:),opz(:)
      real(8), allocatable :: oqx(:),oqy(:),oqz(:)
      real(8), allocatable :: dtx(:),dty(:),dtz(:)
      real(8), allocatable :: dxx(:),dyy(:),dzz(:)
      real(8), allocatable :: uxx(:),uyy(:),uzz(:)
      real(8), allocatable :: xxt(:),yyt(:),zzt(:)
      real(8), allocatable :: txx(:),tyy(:),tzz(:)
      real(8), allocatable :: dxt(:),dyt(:),dzt(:)
      real(8), allocatable :: xxo(:),yyo(:),zzo(:)
      real(8), allocatable :: xx1(:),yy1(:),zz1(:)
      real(8), allocatable :: vxo(:),vyo(:),vzo(:)
      real(8), allocatable :: fmx(:),fmy(:),fmz(:)
      real(8), allocatable :: tqx(:),tqy(:),tqz(:)
      real(8), allocatable :: omxo(:),omyo(:),omzo(:)
      real(8), allocatable :: gvxo(:),gvyo(:),gvzo(:)
      real(8), allocatable :: gcxo(:),gcyo(:),gczo(:)

      dimension strkin(9),strcon(9),strgrp(9),stres0(9),rot(9)
      dimension fail(nnn),uni(9),cell0(9)

      save newjob,volm0,elrc0,virlrc0,dens0
      save igrp1,igrp2,ifre1,ifre2

      data newjob/.true./, beta/7.3728d-3/
      data uni/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/
      
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
      allocate (xxo(mxxdf),yyo(mxxdf),zzo(mxxdf),stat=fail(9))
      allocate (xx1(mxxdf),yy1(mxxdf),zz1(mxxdf),stat=fail(10))
      allocate (vxo(msatms),vyo(msatms),vzo(msatms),stat=fail(11))
      allocate (fmx(msgrp),fmy(msgrp),fmz(msgrp),stat=fail(12))
      allocate (tqx(msgrp),tqy(msgrp),tqz(msgrp),stat=fail(13))
      allocate (omxo(msgrp),omyo(msgrp),omzo(msgrp),stat=fail(14))
      allocate (gvxo(msgrp),gvyo(msgrp),gvzo(msgrp),stat=fail(15))
      allocate (gcxo(msgrp),gcyo(msgrp),gczo(msgrp),stat=fail(16))
      do i=1,nnn
        if(fail(i).ne.0)call error(idnode,1570)
      enddo

      if(newjob) then

c     store initial values of volume, long range corrections etc

        volm0=volm
        elrc0=elrc
        virlrc0=virlrc
        fail(1)=0
        allocate (dens0(mxatyp),stat=fail(1))
        if(fail(1).ne.0)call error(idnode,1580)
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
        
        safe= (igrp2-igrp1+1.le.msgrp) 
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
      
c     temporary stress tensor accumulators and new cell

      vircon=0.d0
      do i=1,9

        strcns(i)=0.d0
        cell0(i)=cell(i)
        stres0(i)=stress(i)

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
        omxo(jg)=omx(ig)
        omyo(jg)=omy(ig)
        omzo(jg)=omz(ig)
        gcxo(jg)=gcmx(ig)
        gcyo(jg)=gcmy(ig)
        gczo(jg)=gcmz(ig)
        gvxo(jg)=gvxx(ig)
        gvyo(jg)=gvyy(ig)
        gvzo(jg)=gvzz(ig)

      enddo

c     construct current bond vectors - required by shake
      
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
      
c     calculate atom displacement from coms
      
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

c     estimate velocity at half-time step
      
      j=0
      do ifre=ifre1,ifre2
        
        j=j+1
        i=lstfre(ifre)

        vxx(i)=vxo(j)+pt5*tstep*rmass(i)*fxx(i)
        vyy(i)=vyo(j)+pt5*tstep*rmass(i)*fyy(i)
        vzz(i)=vzo(j)+pt5*tstep*rmass(i)*fzz(i)

      enddo

c     estimate kinetic energy of free atoms

      engfke=getkinf(ntfree,idnode,mxnode)

c     stress tensor of free atoms

      call kinstressf(ntfree,idnode,mxnode,strkin)

c     estimate translational kinetic energy of rigid bodies
      
      jg=0
      jr=0
      do ig=igrp1,igrp2
        
        jg=jg+1
        id=lstgtp(ig)
        
c     calculate net force on rigid body
        
        fmx(jg)=0.d0
        fmy(jg)=0.d0
        fmz(jg)=0.d0
        
        do j=1,numgsit(id)
          
          jr=jr+1
          i=lstrgd(jr)

c     forces on centre of mass
          
          fmx(jg)=fmx(jg)+fxx(i)
          fmy(jg)=fmy(jg)+fyy(i)
          fmz(jg)=fmz(jg)+fzz(i)
          
        enddo

c     centre of mass velocities at half-step
        
        gvxx(ig)=gvxo(jg)+pt5*tstep*fmx(jg)/gmass(id)
        gvyy(ig)=gvyo(jg)+pt5*tstep*fmy(jg)/gmass(id)
        gvzz(ig)=gvzo(jg)+pt5*tstep*fmz(jg)/gmass(id)
        
      enddo

c     translational kinetic energy of rigid body
      
      engtrn=getkint(ngrp,idnode,mxnode)
      
c     stress tensor of rigid body
      
      call kinstressg(ngrp,idnode,mxnode,strgrp)
      
c     calculate rigid body contribution to stress tensor

      call bodystress(idnode,mxnode,ngrp,vircom,strbod,dtx,dty,dtz)

c     calculate torques in lab frame
      
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

        opx(jg)=omx(ig)
        opy(jg)=omy(ig)
        opz(jg)=omz(ig)

c       iterate angular velocity for time step n (e. yezdimer)

        do i=1,5
        
          trx=(tqx(jg)*rot(1)+tqy(jg)*rot(4)+tqz(jg)*rot(7))*
     x      rotinx(id,2)
     x     +(rotiny(id,1)-rotinz(id,1))*opy(jg)*opz(jg)*
     x      rotinx(id,2)
          try=(tqx(jg)*rot(2)+tqy(jg)*rot(5)+tqz(jg)*rot(8))*
     x      rotiny(id,2)
     x     +(rotinz(id,1)-rotinx(id,1))*opz(jg)*opx(jg)*
     x      rotiny(id,2)
          trz=(tqx(jg)*rot(3)+tqy(jg)*rot(6)+tqz(jg)*rot(9))*
     x      rotinz(id,2)
     x     +(rotinx(id,1)-rotiny(id,1))*opx(jg)*opy(jg)*
     x      rotinz(id,2)

c     improved angular velocity at time step n

          opx(jg)=omx(ig)+pt5*tstep*trx
          opy(jg)=omy(ig)+pt5*tstep*try
          opz(jg)=omz(ig)+pt5*tstep*trz

        enddo

c     store angular velocity at timestep n
        
        omx(ig)=opx(jg)
        omy(ig)=opy(jg)
        omz(ig)=opz(jg)
        
      enddo

c     calculate rotational kinetic energy of rigid bodies

      engrot=getkinr(ngrp,idnode,mxnode)

c     complete stress tensor

      do i=1,9
        stress(i)=stres0(i)+strkin(i)+strgrp(i)+strbod(i)
      enddo
      
c     find eta - taup is the relaxation time

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

c     temperature scaling  coefficient - taut is the relaxation time
      
      engke=engfke+engtrn
      engtot=engke+engrot
      chit0=sqrt(1.d0+tstep/taut*(sigma/engtot-1.d0))
      
c     begin iterations !!-----------------------------------------------

      mxiter=5
      if(ntcons.eq.0) mxiter=mxiter-1
      
      do iter=1,mxiter

c     unconstrained new positions
        
        j=0
        do ifre=ifre1,ifre2

          j=j+1
          i=lstfre(ifre)

c     advance velocity using leapfrog

          uxx(i)=(vxo(j)+tstep*rmass(i)*fxx(i))*chit0
          uyy(i)=(vyo(j)+tstep*rmass(i)*fyy(i))*chit0
          uzz(i)=(vzo(j)+tstep*rmass(i)*fzz(i))*chit0

c     update positions
          
          xxx(i)=tstep*uxx(i)+eta(1)*xxo(j)+eta(4)*yyo(j)+eta(7)*zzo(j)
          yyy(i)=tstep*uyy(i)+eta(2)*xxo(j)+eta(5)*yyo(j)+eta(8)*zzo(j)
          zzz(i)=tstep*uzz(i)+eta(3)*xxo(j)+eta(6)*yyo(j)+eta(9)*zzo(j)

        enddo

c     estimate new cell tensor

        call mat_mul(eta,cell0,cell)

        if(ntcons.eq.0) safe=.true.
        if(ntcons.gt.0) then

c     store integrated positions

          j=0
          do ifre=ifre1,ifre2

            j=j+1
            i=lstfre(ifre)
            xx1(j)=xxx(i)
            yy1(j)=yyy(i)
            zz1(j)=zzz(i)

          enddo

c     global exchange of configuration data
          
          if(mxnode.gt.1)call merge1
     x      (idnode,mxnode,natms,lstme,xxx,yyy,zzz,buffer)
          
c     apply constraint correction
          
          call rdshake_1
     x      (safe,lshmov,idnode,imcon,mxnode,natms,nscons,
     x      tolnce,tstep,viracc,dxx,dyy,dzz,dxt,dyt,dzt,
     x      txx,tyy,tzz,xxt,yyt,zzt,strcon)

c     contribution to constraint virial 
          
          vircon=vircon+viracc
          do i=1,9
            strcns(i)=strcns(i)+strcon(i)
          enddo

c     calculate other constraint corrections
          
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
        
c     calculate new kinetic energy 
        
        engfke=getkinf(ntfree,idnode,mxnode)
        
c     calculate current stress tensor

        call kinstressf(ntfree,idnode,mxnode,strkin)

c     restore half step velocity
          
        if(iter.eq.mxiter)then
          
          do ifre=ifre1,ifre2
            
            i=lstfre(ifre)
            vxx(i)=uxx(i)
            vyy(i)=uyy(i)
            vzz(i)=uzz(i)
            
          enddo
          
        endif
        
c     ********: rigid body motion - thermostated  :***********

c     ***** step 1 : integrate centre of mass motion *********
        
        jg =0
        do ig=igrp1,igrp2

          jg=jg+1
          id=lstgtp(ig)

c     calculate thermostated velocities
          
          rgmas=1.d0/gmass(id)
          uxx(ig)=(gvxo(jg)+tstep*(fmx(jg)*rgmas))*chit0
          uyy(ig)=(gvyo(jg)+tstep*(fmy(jg)*rgmas))*chit0
          uzz(ig)=(gvzo(jg)+tstep*(fmz(jg)*rgmas))*chit0

c     update positions : 

          gcmx(ig)=tstep*uxx(ig)+eta(1)*gcxo(jg)+eta(4)*gcyo(jg)+
     x      eta(7)*gczo(jg)
          gcmy(ig)=tstep*uyy(ig)+eta(2)*gcxo(jg)+eta(5)*gcyo(jg)+
     x      eta(8)*gczo(jg)
          gcmz(ig)=tstep*uzz(ig)+eta(3)*gcxo(jg)+eta(6)*gcyo(jg)+
     x      eta(9)*gczo(jg)

c     full step com velocity

          gvxx(ig)=0.5d0*(gvxo(jg)+uxx(ig))
          gvyy(ig)=0.5d0*(gvyo(jg)+uyy(ig))
          gvzz(ig)=0.5d0*(gvzo(jg)+uzz(ig))
          
        enddo

c     calculate rigid body kinetic energy and stress tensor
        
        engtrn=getkint(ngrp,idnode,mxnode)

        call kinstressg(ngrp,idnode,mxnode,strgrp)

c     restore half step velocity

        do ig=igrp1,igrp2

          gvxx(ig)=uxx(ig)
          gvyy(ig)=uyy(ig)
          gvzz(ig)=uzz(ig)

        enddo

c     ****** step 2 : integrate rotational motion **********
        
        jg=0
        safeq=.true.
        do ig=igrp1,igrp2
          
          jg=jg+1
          id=lstgtp(ig)

c     current rotational matrix 
          
          call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)
        
          trx=(tqx(jg)*rot(1)+tqy(jg)*rot(4)+tqz(jg)*rot(7))*
     x      rotinx(id,2)
     x     +(rotiny(id,1)-rotinz(id,1))*opy(jg)*opz(jg)*
     x      rotinx(id,2)
          try=(tqx(jg)*rot(2)+tqy(jg)*rot(5)+tqz(jg)*rot(8))*
     x      rotiny(id,2)
     x     +(rotinz(id,1)-rotinx(id,1))*opz(jg)*opx(jg)*
     x      rotiny(id,2)
          trz=(tqx(jg)*rot(3)+tqy(jg)*rot(6)+tqz(jg)*rot(9))*
     x      rotinz(id,2)
     x     +(rotinx(id,1)-rotiny(id,1))*opx(jg)*opy(jg)*
     x      rotinz(id,2)

c     scaled angular velocity at time step n
            
          omx(ig)=(omxo(jg)+pt5*tstep*trx)*chit0
          omy(ig)=(omyo(jg)+pt5*tstep*try)*chit0
          omz(ig)=(omzo(jg)+pt5*tstep*trz)*chit0

c     angular velocity at time step n+1/2
          
          uxx(ig)=(omxo(jg)+tstep*trx)*chit0
          uyy(ig)=(omyo(jg)+tstep*try)*chit0
          uzz(ig)=(omzo(jg)+tstep*trz)*chit0

c     angular velocity at time step n+1 (needed for quat algorithm)
          
          oqx(jg)=(omxo(jg)+1.5d0*tstep*trx)*chit0
          oqy(jg)=(omyo(jg)+1.5d0*tstep*try)*chit0
          oqz(jg)=(omzo(jg)+1.5d0*tstep*trz)*chit0

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
        
c     complete stress tensor - add all contributions

        do i=1,9
          stress(i)=stres0(i)+strkin(i)+strgrp(i)+strbod(i)+strcns(i)
        enddo

c     improved prediction of eta and chit 

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

        engke=engfke+engtrn
        engtot=engke+engrot
        chit0=sqrt(1.d0+tstep/taut*(sigma/engtot-1.d0))

c     end of thermostat/barostat iterations

      enddo

c     update cell vectors

      call mat_mul(eta,cell0,cell)

c     update volume
      
      volm=volm*eta(1)*eta(5)*eta(9)

c     adjust long range corrections and number density
      
      elrc=elrc0*(volm0/volm)
      virlrc=virlrc0*(volm0/volm)
      
      do k=1,ntpatm
        dens(k)=dens0(k)*(volm0/volm)
      enddo

c     assign new quaternions

      call update_quaternions
     x  (safeq,igrp1,igrp2,tstep,quattol,opx,opy,opz,oqx,oqy,oqz)

c     minimum images 
      
      call images(imcon,0,1,natms,cell,xxx,yyy,zzz)
      call images(imcon,idnode,mxnode,ngrp,cell,gcmx,gcmy,gcmz)

c     new atomic positions for atoms in rigid bodies
      
      jr=0
      do ig=igrp1,igrp2

        id=lstgtp(ig)

c     new rotational matrix
        
        call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)
        
        do j=1,numgsit(id)
          
          jr=jr+1
          i=lstrgd(jr)
          
          xxx(i)=rot(1)*gxx(id,j)+rot(2)*gyy(id,j)+rot(3)*gzz(id,j)
     x     +gcmx(ig)
          yyy(i)=rot(4)*gxx(id,j)+rot(5)*gyy(id,j)+rot(6)*gzz(id,j)
     x     +gcmy(ig)
          zzz(i)=rot(7)*gxx(id,j)+rot(8)*gyy(id,j)+rot(9)*gzz(id,j)
     x     +gcmz(ig)

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
        if(ntcons.gt.0)
     x    call merge1(idnode,mxnode,natms,lstme,fxx,fyy,fzz,buffer)
      
      endif

c     ensure all atoms are within cell boundaries

      call images(imcon,0,1,natms,cell,xxx,yyy,zzz)

c     deallocate work arrays

      deallocate (opx,opy,opz,oqx,oqy,oqz,stat=fail(1))
      deallocate (dtx,dty,dtz,dxx,dyy,dzz,stat=fail(2))
      deallocate (uxx,uyy,uzz,txx,tyy,tzz,stat=fail(3))
      deallocate (dxt,dyt,dzt,xxo,yyo,zzo,stat=fail(4))
      deallocate (xx1,yy1,zz1,vxo,vyo,vzo,stat=fail(5))
      deallocate (fmx,fmy,fmz,tqx,tqy,tqz,stat=fail(6))
      deallocate (omxo,omyo,omzo,gcxo,gcyo,gczo,stat=fail(7))
      deallocate (xxt,yyt,zzt,gvxo,gvyo,gvzo,stat=fail(8))
      
      return
      end subroutine nstq_b1

      subroutine nstq_h1
     x  (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x  ntcons,ntfree,ntpatm,mode,chit,conint,consv,elrc,engke,engrot,
     x  virlrc,press,quattol,sigma,taup,taut,temp,tolnce,tstep,
     x  vircom,vircon,volm)
      
c***********************************************************************
c     
c     dlpoly routine to integrate rigid body equations of motion
c     using implicit leapfrog quaternion algorithm
c     plus bond constraints - PROVIDED rigid body sites
c     and constraint sites do not coincide.
c     
c     dl_poly subroutine for integrating newtonian equations of
c     motion in molecular dynamics - verlet leapfrog with Hoover like
c     thermostat and barostat. (cell may change shape).
c     
c     parallel replicated data version : block data
c     
c     for systems using bond constraints
c     
c     tqx,tqy,tqz = torque in lab fixed frame (input)
c     omx,omy,omz = angular velocity in body fixed frame (principle axis)
c     rotinx,y,z  = rotational inertia in body fixed frame
c     
c     copyright daresbury laboratory 1995
c     author      t.forester june  1995
c     amended     w.smith sep 1999 : euler equation
c     
c**********************************************************************
      
      implicit none
      
      integer, parameter :: nnn=16
      real(8), parameter :: pt5=0.5d0

      logical safe,lshmov,newjob,safeq
      integer imcon,idnode,mxnode,natms,ngrp,nscons,ntcons,ntfree
      integer ntpatm,fail,i,igrp,igrp1,igrp2,ifre,ifre1,ifre2
      integer j,k,ig,jg,jr,id,iter,mxiter,idum,mode
      real(8) chit,conint,consv,elrc,engke,engrot,virlrc,press,quattol
      real(8) sigma,taup,taut,temp,tolnce,tstep,vircom,vircon,volm
      real(8) strkin,strcon,strgrp,eta0,etanew,rot,cell0,volm0,stres0
      real(8) elrc0,virlrc0,rstep,rtsq,qmass,pmass,totmas
      real(8) engtrn,trx,try,trz,engtot,engfke,fac,etadot
      real(8) chitp,chitnew,chit0,xxa,yya,zza,viracc,rgmas,uni
      real(8) delx,dely,delz,vold,cons1,cons2,cons3
      real(8) vaa,vbb,vcc,chip,com,vom

      real(8), allocatable :: dens0(:)
      real(8), allocatable :: opx(:),opy(:),opz(:)
      real(8), allocatable :: oqx(:),oqy(:),oqz(:)
      real(8), allocatable :: dtx(:),dty(:),dtz(:)
      real(8), allocatable :: dxx(:),dyy(:),dzz(:)
      real(8), allocatable :: uxx(:),uyy(:),uzz(:)
      real(8), allocatable :: xxt(:),yyt(:),zzt(:)
      real(8), allocatable :: txx(:),tyy(:),tzz(:)
      real(8), allocatable :: dxt(:),dyt(:),dzt(:)
      real(8), allocatable :: xxo(:),yyo(:),zzo(:)
      real(8), allocatable :: xx1(:),yy1(:),zz1(:)
      real(8), allocatable :: vxo(:),vyo(:),vzo(:)
      real(8), allocatable :: fmx(:),fmy(:),fmz(:)
      real(8), allocatable :: tqx(:),tqy(:),tqz(:)
      real(8), allocatable :: omxo(:),omyo(:),omzo(:)
      real(8), allocatable :: gvxo(:),gvyo(:),gvzo(:)
      real(8), allocatable :: gcxo(:),gcyo(:),gczo(:)

      dimension eta0(9),etanew(9),rot(9),cell0(9),fail(nnn),uni(9)
      dimension strkin(9),strcon(9),strgrp(9),stres0(9),com(3),vom(3)
      
      save newjob,volm0,elrc0,virlrc0,dens0,pmass,qmass,totmas
      save igrp1,igrp2,ifre1,ifre2
      
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
      allocate (xxt(mxatms),yyt(mxatms),zzt(mxatms),stat=fail(6))
      allocate (txx(mxatms),tyy(mxatms),tzz(mxatms),stat=fail(7))
      allocate (dxt(mxcons),dyt(mxcons),dzt(mxcons),stat=fail(8))
      allocate (xxo(msatms),yyo(msatms),zzo(msatms),stat=fail(9))
      allocate (xx1(msatms),yy1(msatms),zz1(msatms),stat=fail(10))
      allocate (vxo(msatms),vyo(msatms),vzo(msatms),stat=fail(11))
      allocate (fmx(msgrp),fmy(msgrp),fmz(msgrp),stat=fail(12))
      allocate (tqx(msgrp),tqy(msgrp),tqz(msgrp),stat=fail(13))
      allocate (omxo(msgrp),omyo(msgrp),omzo(msgrp),stat=fail(14))
      allocate (gvxo(msgrp),gvyo(msgrp),gvzo(msgrp),stat=fail(15))
      allocate (gcxo(msgrp),gcyo(msgrp),gczo(msgrp),stat=fail(16))
      do i=1,nnn
        if(fail(i).ne.0)call error(idnode,1590)
      enddo

      if(newjob) then
        
c     store initial values of volume, long range corrections etc
      
        volm0=volm
        elrc0=elrc
        virlrc0=virlrc
        fail(1)=0
        allocate (dens0(mxatyp),stat=fail(1))
        if(fail(1).ne.0)call error(idnode,1600)

        do i=1,ntpatm
          dens0(i)=dens(i)
        enddo
        
c     inertia parameter for Nose-Hoover thermostat
        
        qmass=2.0d0*sigma*taut**2
        pmass=2.0d0*sigma*taup**2

c     calculate total system mass 

        totmas=getmass(natms,idnode,mxnode)
        
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
      
c     initialise stress tensor accumulators

      vircon=0.d0

      do i=1,9

        strcns(i)=0.d0
        cell0(i)=cell(i)
        stres0(i)=stress(i)
        
      enddo

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
        omxo(jg)=omx(ig)
        omyo(jg)=omy(ig)
        omzo(jg)=omz(ig)
        gcxo(jg)=gcmx(ig)
        gcyo(jg)=gcmy(ig)
        gczo(jg)=gcmz(ig)
        gvxo(jg)=gvxx(ig)
        gvyo(jg)=gvyy(ig)
        gvzo(jg)=gvzz(ig)
        
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

c     calculate centre of mass
      
      call getcom(natms,idnode,mxnode,totmas,com)
      
c     estimate velocity at full step
      
      j=0
      do ifre=ifre1,ifre2
        
        j=j+1
        i=lstfre(ifre)

        vxx(i)=vxo(j)+pt5*tstep*rmass(i)*fxx(i)
        vyy(i)=vyo(j)+pt5*tstep*rmass(i)*fyy(i)
        vzz(i)=vzo(j)+pt5*tstep*rmass(i)*fzz(i)
        
      enddo
      
c     estimate kinetic energy
        
      engfke=getkinf(ntfree,idnode,mxnode)

c     estimate stress tensor
      
      call kinstressf(ntfree,idnode,mxnode,strkin)
        
c     estimate rigid body translational kinetic energy
      
      jg=0
      jr=0
      do ig=igrp1,igrp2
        
        jg=jg+1
        id=lstgtp(ig)

c     calculate net force on rigid body
        
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

c     centre of mass velocities at half-step
        
        gvxx(ig)=gvxo(jg)+pt5*tstep/gmass(id)*fmx(jg)
        gvyy(ig)=gvyo(jg)+pt5*tstep/gmass(id)*fmy(jg)
        gvzz(ig)=gvzo(jg)+pt5*tstep/gmass(id)*fmz(jg)

      enddo
      
c     rigid body translational kinetic energy
        
      engtrn=getkint(ngrp,idnode,mxnode)

c     rigid body stress tensor
        
      call kinstressg(ngrp,idnode,mxnode,strgrp)
      
c     calculate rgid body contribution to stress tensor
      
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

        opx(jg)=omx(ig)
        opy(jg)=omy(ig)
        opz(jg)=omz(ig)

c       iterate angular velocity for time step n (e. yezdimer)

        do i=1,5
        
          trx=(tqx(jg)*rot(1)+tqy(jg)*rot(4)+tqz(jg)*rot(7))*
     x      rotinx(id,2)
     x     +(rotiny(id,1)-rotinz(id,1))*opy(jg)*opz(jg)*
     x      rotinx(id,2)
          try=(tqx(jg)*rot(2)+tqy(jg)*rot(5)+tqz(jg)*rot(8))*
     x      rotiny(id,2)
     x     +(rotinz(id,1)-rotinx(id,1))*opz(jg)*opx(jg)*
     x      rotiny(id,2)
          trz=(tqx(jg)*rot(3)+tqy(jg)*rot(6)+tqz(jg)*rot(9))*
     x      rotinz(id,2)
     x     +(rotinx(id,1)-rotiny(id,1))*opx(jg)*opy(jg)*
     x      rotinz(id,2)

c     improved angular velocity at time step n

          opx(jg)=omx(ig)+pt5*tstep*trx
          opy(jg)=omy(ig)+pt5*tstep*try
          opz(jg)=omz(ig)+pt5*tstep*trz

        enddo

c     store angular velcoity at timestep n
        
        omx(ig)=opx(jg)
        omy(ig)=opy(jg)
        omz(ig)=opz(jg)
        
      enddo

c     calculate rotational kinetic energy of rigid bodies

      engrot=getkinr(ngrp,idnode,mxnode)

c     complete stress tensor

      do i=1,9
        stress(i)=stres0(i)+strkin(i)+strgrp(i)+strbod(i)
      enddo
      
c     propagate eta
      
      fac=9.d0
      do i=1,9
        etanew(i)=eta(i)+tstep*((stress(i)-press*volm*uni(i))/pmass-
     x    chit*eta(i))
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
        eta0(i)=0.5d0*(etanew(i)+eta(i)) 
      enddo

c     propagate chit
      
      engke=engfke+engtrn
      engtot=engke+engrot
      etadot=sdot0(9,eta0,eta0)
      if(mode.eq.2)etadot=etadot-eta0(1)**2
      chitp=(2.d0*(engtot-sigma)+pmass*etadot-fac*boltz*temp)/qmass
      chitnew=chit+tstep*chitp
      chit0=0.5d0*(chit+chitnew)

c     begin iterations !!-----------------------------------------------

      mxiter=5
      if(ntcons.eq.0) mxiter=mxiter-1
      
      do iter=1,mxiter

c     unconstrained new positions
        
        j=0
        do ifre=ifre1,ifre2

          j=j+1
          i=lstfre(ifre)

c     advance velocity using leapfrog

          uxx(i)=vxo(j)+tstep*(fxx(i)*rmass(i)-
     x      (chit0+eta0(1))*vxx(i)-eta0(4)*vyy(i)-eta0(7)*vzz(i))
          uyy(i)=vyo(j)+tstep*(fyy(i)*rmass(i)-
     x      eta0(2)*vxx(i)-(eta0(5)+chit0)*vyy(i)-eta0(8)*vzz(i))
          uzz(i)=vzo(j)+tstep*(fzz(i)*rmass(i)-
     x      eta0(3)*vxx(i)-eta0(6)*vyy(i)-(eta0(9)+chit0)*vzz(i))

c     advance positions using leapfrog

          xxa=(xxx(i)+xxo(j))*0.5d0-com(1)
          yya=(yyy(i)+yyo(j))*0.5d0-com(2)
          zza=(zzz(i)+zzo(j))*0.5d0-com(3)

          xxx(i)=xxo(j)+tstep*(uxx(i)+
     x      eta0(1)*xxa+eta0(4)*yya+eta0(7)*zza)
          yyy(i)=yyo(j)+tstep*(uyy(i)+
     x      eta0(2)*xxa+eta0(5)*yya+eta0(8)*zza)
          zzz(i)=zzo(j)+tstep*(uzz(i)+
     x      eta0(3)*xxa+eta0(6)*yya+eta0(9)*zza)
          
        enddo

c     estimate new cell parameters
        
        do i=1,9
          cell(i)=cell0(i)
        enddo
        call cell_propagate(tstep,cell,etanew)

        if(ntcons.eq.0) safe=.true.
        if(ntcons.gt.0) then

c     store integrated positions

          j=0
          do ifre=ifre1,ifre2

            j=j+1
            i=lstfre(ifre)
            xx1(j)=xxx(i)
            yy1(j)=yyy(i)
            zz1(j)=zzz(i)

          enddo

c     global exchange of configuration data
          
          if(mxnode.gt.1)call merge1
     x      (idnode,mxnode,natms,lstme,xxx,yyy,zzz,buffer)

c     apply constraint correction
          
          call rdshake_1
     x      (safe,lshmov,idnode,imcon,mxnode,natms,nscons,
     x      tolnce,tstep,viracc,dxx,dyy,dzz,dxt,dyt,dzt,
     x      txx,tyy,tzz,xxt,yyt,zzt,strcon)

c     contribution to constraint virial 
          
          vircon=vircon+viracc
          do i=1,9
            strcns(i)=strcns(i)+strcon(i)
          enddo
          
c     calculate other constraint corrections
          
          j=0
          rstep=1.d0/tstep
          rtsq=1.d0/tstep**2
          do ifre=ifre1,ifre2

            j=j+1
            i=lstfre(ifre)

c     calculate corrected velocity
            
            uxx(i)=uxx(i)+(xxx(i)-xx1(j))*rstep
            uyy(i)=uyy(i)+(yyy(i)-yy1(j))*rstep
            uzz(i)=uzz(i)+(zzz(i)-zz1(j))*rstep

c     calculate the corrected forces
            
            fxx(i)=fxx(i)+(xxx(i)-xx1(j))*weight(i)*rtsq
            fyy(i)=fyy(i)+(yyy(i)-yy1(j))*weight(i)*rtsq
            fzz(i)=fzz(i)+(zzz(i)-zz1(j))*weight(i)*rtsq
            
          enddo

c     end of shake corrections
          
        endif
        
c     estimate velocity at the full step
        
        j=0
        do ifre=ifre1,ifre2
          
          j=j+1
          i=lstfre(ifre)
          vxx(i)=0.5d0*(uxx(i)+vxo(j))
          vyy(i)=0.5d0*(uyy(i)+vyo(j))
          vzz(i)=0.5d0*(uzz(i)+vzo(j))
          
        enddo
        
c     calculate kinetic energy
        
        engfke=getkinf(ntfree,idnode,mxnode)
        
c     kinetic contribution to stress tensor
        
        call kinstressf(ntfree,idnode,mxnode,strkin)
        
c     restore free half step velocity
          
        if(iter.eq.mxiter)then
          
          do ifre=ifre1,ifre2
            
            i=lstfre(ifre)
            vxx(i)=uxx(i)
            vyy(i)=uyy(i)
            vzz(i)=uzz(i)
            
          enddo
          
        endif
        
c     ********: rigid body motion - thermostated  :***********
        
c     ***** step 1 : integrate centre of mass motion *********
        
        jg =0
        do ig=igrp1,igrp2

          jg=jg+1
          id=lstgtp(ig)

c     calculate thermostated velocities
          
          rgmas=1.d0/gmass(id)
          uxx(ig)=gvxo(jg)+tstep*(fmx(jg)*rgmas-
     x      (chit0+eta0(1))*gvxx(ig)-eta0(4)*gvyy(ig)-eta0(7)*gvzz(ig))
          uyy(ig)=gvyo(jg)+tstep*(fmy(jg)*rgmas-
     x      eta0(2)*gvxx(ig)-(eta0(5)+chit0)*gvyy(ig)-eta0(8)*gvzz(ig))
          uzz(ig)=gvzo(jg)+tstep*(fmz(jg)*rgmas-
     x      eta0(3)*gvxx(ig)-eta0(6)*gvyy(ig)-(eta0(9)+chit0)*gvzz(ig))

c     advance positions using leapfrog

          xxa=(gcmx(ig)+gcxo(jg))*0.5d0-com(1)
          yya=(gcmy(ig)+gcyo(jg))*0.5d0-com(2)
          zza=(gcmz(ig)+gczo(jg))*0.5d0-com(3)

          gcmx(ig)=gcxo(jg)+tstep*(uxx(ig)+
     x      eta0(1)*xxa+eta0(4)*yya+eta0(7)*zza)
          gcmy(ig)=gcyo(jg)+tstep*(uyy(ig)+
     x      eta0(2)*xxa+eta0(5)*yya+eta0(8)*zza)
          gcmz(ig)=gczo(jg)+tstep*(uzz(ig)+
     x      eta0(3)*xxa+eta0(6)*yya+eta0(9)*zza)

c     full step com velocity

          gvxx(ig)=0.5d0*(gvxo(jg)+uxx(ig))
          gvyy(ig)=0.5d0*(gvyo(jg)+uyy(ig))
          gvzz(ig)=0.5d0*(gvzo(jg)+uzz(ig))
          
        enddo

c     calculate rigid body kinetic energy and stress tensor
        
        engtrn=getkint(ngrp,idnode,mxnode)

        call kinstressg(ngrp,idnode,mxnode,strgrp)

c     restore half step velocity

        if(iter.eq.mxiter)then
          
          do ig=igrp1,igrp2
            
            gvxx(ig)=uxx(ig)
            gvyy(ig)=uyy(ig)
            gvzz(ig)=uzz(ig)
            
          enddo
          
        endif
        
c     ****** step 2 : integrate rotational motion **********
        
        jg=0
        safeq=.true.
        do ig=igrp1,igrp2
          
          jg=jg+1
          id=lstgtp(ig)

c     current rotational matrix 
          
          call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)
        
          trx=(tqx(jg)*rot(1)+tqy(jg)*rot(4)+tqz(jg)*rot(7))*
     x      rotinx(id,2)
     x     +(rotiny(id,1)-rotinz(id,1))*opy(jg)*opz(jg)*
     x      rotinx(id,2)
          try=(tqx(jg)*rot(2)+tqy(jg)*rot(5)+tqz(jg)*rot(8))*
     x      rotiny(id,2)
     x     +(rotinz(id,1)-rotinx(id,1))*opz(jg)*opx(jg)*
     x      rotiny(id,2)
          trz=(tqx(jg)*rot(3)+tqy(jg)*rot(6)+tqz(jg)*rot(9))*
     x      rotinz(id,2)
     x     +(rotinx(id,1)-rotiny(id,1))*opx(jg)*opy(jg)*
     x      rotinz(id,2)

c     correction due to thermostat

          delx=tstep*(trx-chit0*pt5*(omx(ig)+omxo(jg)))
          dely=tstep*(try-chit0*pt5*(omy(ig)+omyo(jg)))
          delz=tstep*(trz-chit0*pt5*(omz(ig)+omzo(jg)))

c     angular velocity at time step n
          
          omx(ig)=omxo(jg)+delx*pt5
          omy(ig)=omyo(jg)+dely*pt5
          omz(ig)=omzo(jg)+delz*pt5

c     angular velocity at time step n+1/2
          
          uxx(ig)=omxo(jg)+delx
          uyy(ig)=omyo(jg)+dely
          uzz(ig)=omzo(jg)+delz

c     angular velocity at time step n+1  (needed for quat algorithm)
          
          oqx(jg)=omxo(jg)+delx*1.5d0
          oqy(jg)=omyo(jg)+dely*1.5d0
          oqz(jg)=omzo(jg)+delz*1.5d0

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
        
c     complete stress tensor - add all contributions

        do i=1,9
          stress(i)=stres0(i)+strkin(i)+strgrp(i)+strbod(i)+strcns(i)
        enddo

c     improved prediction of eta and chit 

        do i=1,9
          etanew(i)=eta(i)+tstep*((stress(i)-press*volm*uni(i))/pmass-
     x      chit0*eta0(i))
        enddo
        if(mode.gt.0)then
          etanew(3)=0.d0
          etanew(6)=0.d0
          etanew(7)=0.d0
          etanew(8)=0.d0
          if(mode.lt.3)then
            etanew(2)=0.d0
            etanew(4)=0.d0
            if(mode.eq.2)then
              etanew(1)=0.5d0*(etanew(1)+etanew(5))
              etanew(5)=etanew(1)
            endif
          endif
        endif
        do i=1,9
          eta0(i)=0.5d0*(etanew(i)+eta(i)) 
        enddo

        engke=engfke+engtrn
        engtot=engke+engrot
        etadot=sdot0(9,eta0,eta0)
        if(mode.eq.2)etadot=etadot-eta0(1)**2
        chitp=(2.d0*(engtot-sigma)+pmass*etadot-fac*boltz*temp)/qmass
        chitnew=chit+tstep*chitp
        chit0=0.5d0*(chit+chitnew)

c     end of thermostat/barostat iterations

      enddo

c     assign new quaternions

      call update_quaternions
     x  (safeq,igrp1,igrp2,tstep,quattol,opx,opy,opz,oqx,oqy,oqz)

c     minimum images 
      
      call images(imcon,0,1,natms,cell,xxx,yyy,zzz)
      call images(imcon,idnode,mxnode,ngrp,cell,gcmx,gcmy,gcmz)

c     new atomic positions for atoms in rigid bodies
      
      jr=0
      do ig=igrp1,igrp2

        id=lstgtp(ig)

c     new rotational matrix
        
        call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)
        
        do j=1,numgsit(id)
          
          jr=jr+1
          i=lstrgd(jr)
          
          xxx(i)=rot(1)*gxx(id,j)+rot(2)*gyy(id,j)+rot(3)*gzz(id,j)
     x     +gcmx(ig)
          yyy(i)=rot(4)*gxx(id,j)+rot(5)*gyy(id,j)+rot(6)*gzz(id,j)
     x     +gcmy(ig)
          zzz(i)=rot(7)*gxx(id,j)+rot(8)*gyy(id,j)+rot(9)*gzz(id,j)
     x     +gcmz(ig)

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
      cons1=0.5d0*qmass*chit0**2
      cons2=press*vold
      etadot=sdot0(9,eta0,eta0)
      if(mode.eq.2)etadot=etadot-eta0(1)**2
      cons3=0.5d0*pmass*etadot
      consv=conint+cons1+cons2+cons3
      
      if(mxnode.gt.1) then

c     merge new group coordinates and velocities
      
        call merge(idnode,mxnode,ngrp,mxbuff,gcmx,gcmy,gcmz,buffer)
        call merge(idnode,mxnode,ngrp,mxbuff,gvxx,gvyy,gvzz,buffer)
      
c     merge new atomic coordinates and velocities
      
        call merge1(idnode,mxnode,natms,lstme,xxx,yyy,zzz,buffer)
        call merge1(idnode,mxnode,natms,lstme,vxx,vyy,vzz,buffer)
        if(ntcons.gt.0)
     x    call merge1(idnode,mxnode,natms,lstme,fxx,fyy,fzz,buffer)
      
      endif

c     ensure all atoms are within cell boundaries

      call images(imcon,0,1,natms,cell,xxx,yyy,zzz)

c     deallocate work arrays

      deallocate (opx,opy,opz,oqx,oqy,oqz,stat=fail(1))
      deallocate (dtx,dty,dtz,dxx,dyy,dzz,stat=fail(2))
      deallocate (uxx,uyy,uzz,txx,tyy,tzz,stat=fail(3))
      deallocate (dxt,dyt,dzt,xxo,yyo,zzo,stat=fail(4))
      deallocate (xx1,yy1,zz1,vxo,vyo,vzo,stat=fail(5))
      deallocate (fmx,fmy,fmz,tqx,tqy,tqz,stat=fail(6))
      deallocate (omxo,omyo,omzo,gcxo,gcyo,gczo,stat=fail(7))
      deallocate (xxt,yyt,zzt,gvxo,gvyo,gvzo,stat=fail(8))
      
      return
      end subroutine nstq_h1
      
      end module lf_rotation1_module
