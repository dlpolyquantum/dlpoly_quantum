      module vv_rotation1_module
      
c***********************************************************************
c     
c     dl_poly module 1 for velocity verlet rotational integration 
c     schemes
c     copyright - daresbury laboratory
c     author    - w. smith    aug 2006
c     adapted   - d. quigley - metadynamics
c     
c***********************************************************************
      
      use config_module
      use ensemble_tools_module
      use error_module
      use metafreeze_module, only : lmetadyn
      use property_module
      use rigid_body_module
      use setup_module
      use shake_module
      use site_module
      use vv_motion_module
      use utility_module
      
      contains
      
      subroutine rotate_omega
     x  (idnode,mxnode,ngrp,tstep,p0,p1,p2,p3,dtx,dty,dtz)

c***********************************************************************
c     
c     dl_poly subroutine for updating the angular velocity and momentum
c     for rigid bodies
c     
c     copyright - daresbury laboratory
c     author    - w. smith  sept 2005
c     
c***********************************************************************

      implicit none

      integer idnode,mxnode,ngrp,i,j,jr,jrs,ig,igrp1,igrp2,id
      real(8) ftx,fty,ftz,fmx,fmy,fmz,tstep,tqx,tqy,tqz,tq0,tq1,tq2,tq3

      real(8) p0(mxgrp),p1(mxgrp),p2(mxgrp),p3(mxgrp)
      real(8) dtx(mxatms),dty(mxatms),dtz(mxatms),rot(9)

c     group block indices
        
      igrp1=(idnode*ngrp)/mxnode+1
      igrp2=((idnode+1)*ngrp)/mxnode

      jr=0
      do ig=igrp1,igrp2
        
c     fmx,fmy,fmz represent force on c.o.m.
        
        jrs=jr
        fmx=0.d0
        fmy=0.d0
        fmz=0.d0
        id=lstgtp(ig)
        do j=1,numgsit(id)
          
          jr=jr+1
          i=lstrgd(jr)
          fmx=fmx+fxx(i)
          fmy=fmy+fyy(i)
          fmz=fmz+fzz(i)
          
        enddo
        
c     current rotational matrix
          
        call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)

c     calculate torque in principal frame
          
        jr=jrs
        ftx=0.d0
        fty=0.d0
        ftz=0.d0
        do j=1,numgsit(id)
          
          jr=jr+1
          i=lstrgd(jr)
          ftx=ftx+dty(jr)*fzz(i)-dtz(jr)*fyy(i)
          fty=fty+dtz(jr)*fxx(i)-dtx(jr)*fzz(i)
          ftz=ftz+dtx(jr)*fyy(i)-dty(jr)*fxx(i)
          
        enddo
        tqx=ftx*rot(1)+fty*rot(4)+ftz*rot(7)
        tqy=ftx*rot(2)+fty*rot(5)+ftz*rot(8)
        tqz=ftx*rot(3)+fty*rot(6)+ftz*rot(9)

c     calculate quaternion torques

        tq0=2.0d0*(-q1(ig)*tqx-q2(ig)*tqy-q3(ig)*tqz)
        tq1=2.0d0*( q0(ig)*tqx-q3(ig)*tqy+q2(ig)*tqz)
        tq2=2.0d0*( q3(ig)*tqx+q0(ig)*tqy-q1(ig)*tqz)
        tq3=2.0d0*(-q2(ig)*tqx+q1(ig)*tqy+q0(ig)*tqz)
        
c     update quaternion momentum by half timestep

        p0(ig)=p0(ig)+tq0*tstep
        p1(ig)=p1(ig)+tq1*tstep
        p2(ig)=p2(ig)+tq2*tstep
        p3(ig)=p3(ig)+tq3*tstep

c     update centre of mass velocity by half timestep

        gvxx(ig)=gvxx(ig)+fmx*tstep/gmass(id)
        gvyy(ig)=gvyy(ig)+fmy*tstep/gmass(id)
        gvzz(ig)=gvzz(ig)+fmz*tstep/gmass(id)
        
      enddo

      return
      end subroutine rotate_omega

      subroutine nosquish(ig,tstep,qq0,qq1,qq2,qq3,pp0,pp1,pp2,pp3)

c***********************************************************************
c     
c     dlpoly routine to implement the symplectic no_squish quaternion 
c     algorithm of miller et al j.chem.phys 116 (2002) 8649
c     
c     copyright daresbury laboratory
c     author      m. leslie jan 2004
c     amended     w.smith   mar 2005
c     
c**********************************************************************

      implicit none

      integer m,ig,id
      real(8) zetax,zetay,zetaz,tstep,cs,sn,trstep

      integer, parameter :: mrot=10
      real(8), parameter :: ov4=0.25d0
      real(8), parameter :: ov8=0.125d0

      real(8) qq0(*),qq1(*),qq2(*),qq3(*)
      real(8) pp0(*),pp1(*),pp2(*),pp3(*)

      real(8) qn1(0:3),pq2(0:3),qn2(0:3),pq3(0:3)
      real(8) qn3(0:3),pq4(0:3)

c     rotational time step

      trstep=tstep/dble(mrot)

c     rotation: iterate over mrot rotational time steps

      id=lstgtp(ig)

      do m=1,mrot
        
        zetaz=ov8*rotinz(id,2)*trstep*
     x    (-pp0(ig)*qq3(ig)+pp1(ig)*qq2(ig)-
     x    pp2(ig)*qq1(ig)+pp3(ig)*qq0(ig))
        cs=cos(zetaz)
        sn=sin(zetaz)
        qn1(0)=cs*qq0(ig)-sn*qq3(ig)
        qn1(1)=cs*qq1(ig)+sn*qq2(ig)
        qn1(2)=cs*qq2(ig)-sn*qq1(ig)
        qn1(3)=cs*qq3(ig)+sn*qq0(ig)
        pq2(0)=cs*pp0(ig)-sn*pp3(ig)
        pq2(1)=cs*pp1(ig)+sn*pp2(ig)
        pq2(2)=cs*pp2(ig)-sn*pp1(ig)
        pq2(3)=cs*pp3(ig)+sn*pp0(ig)
        
        zetay=ov8*rotiny(id,2)*trstep*
     x    (-pq2(0)*qn1(2)-pq2(1)*qn1(3)+
     x    pq2(2)*qn1(0)+pq2(3)*qn1(1))
        cs=cos(zetay)
        sn=sin(zetay)
        qn2(0)=cs*qn1(0)-sn*qn1(2)
        qn2(1)=cs*qn1(1)-sn*qn1(3)
        qn2(2)=cs*qn1(2)+sn*qn1(0)
        qn2(3)=cs*qn1(3)+sn*qn1(1)
        pq3(0)=cs*pq2(0)-sn*pq2(2)
        pq3(1)=cs*pq2(1)-sn*pq2(3)
        pq3(2)=cs*pq2(2)+sn*pq2(0)
        pq3(3)=cs*pq2(3)+sn*pq2(1)
        
        zetax=ov4*rotinx(id,2)*trstep*
     x    (-pq3(0)*qn2(1)+pq3(1)*qn2(0)+
     x    pq3(2)*qn2(3)-pq3(3)*qn2(2))
        cs=cos(zetax)
        sn=sin(zetax)
        qn3(0)=cs*qn2(0)-sn*qn2(1)
        qn3(1)=cs*qn2(1)+sn*qn2(0)
        qn3(2)=cs*qn2(2)+sn*qn2(3)
        qn3(3)=cs*qn2(3)-sn*qn2(2)
        pq4(0)=cs*pq3(0)-sn*pq3(1)
        pq4(1)=cs*pq3(1)+sn*pq3(0)
        pq4(2)=cs*pq3(2)+sn*pq3(3)
        pq4(3)=cs*pq3(3)-sn*pq3(2)
        
        zetay=ov8*rotiny(id,2)*trstep*
     x    (-pq4(0)*qn3(2)-pq4(1)*qn3(3)+
     x    pq4(2)*qn3(0)+pq4(3)*qn3(1))
        cs=cos(zetay)
        sn=sin(zetay)
        qn2(0)=cs*qn3(0)-sn*qn3(2)
        qn2(1)=cs*qn3(1)-sn*qn3(3)
        qn2(2)=cs*qn3(2)+sn*qn3(0)
        qn2(3)=cs*qn3(3)+sn*qn3(1)
        pq3(0)=cs*pq4(0)-sn*pq4(2)
        pq3(1)=cs*pq4(1)-sn*pq4(3)
        pq3(2)=cs*pq4(2)+sn*pq4(0)
        pq3(3)=cs*pq4(3)+sn*pq4(1)
        
        zetaz=ov8*rotinz(id,2)*trstep*
     x    (-pq3(0)*qn2(3)+pq3(1)*qn2(2)-
     x    pq3(2)*qn2(1)+pq3(3)*qn2(0))
        cs=cos(zetaz)
        sn=sin(zetaz)
        qq0(ig)=cs*qn2(0)-sn*qn2(3)
        qq1(ig)=cs*qn2(1)+sn*qn2(2)
        qq2(ig)=cs*qn2(2)-sn*qn2(1)
        qq3(ig)=cs*qn2(3)+sn*qn2(0)
        pp0(ig)=cs*pq3(0)-sn*pq3(3)
        pp1(ig)=cs*pq3(1)+sn*pq3(2)
        pp2(ig)=cs*pq3(2)-sn*pq3(1)
        pp3(ig)=cs*pq3(3)+sn*pq3(0)
        
      enddo

      return
      end subroutine nosquish

      subroutine nveqvv_1
     x  (safe,lshmov,isw,imcon,idnode,mxnode,natms,ngrp,nscons,
     x  ntcons,ntfree,engke,engrot,tolnce,tstep,vircom,vircon)
      
c***********************************************************************
c     
c     dlpoly routine to integrate rigid body equations of motion
c     using the symplectic no_squish quaternion algorithm of 
c     miller et al j.chem.phys 116 (2002) 8649
c     plus bond constraints - provided rigid body sites
c     and constraint sites do not coincide.
c     
c     parallel replicated data version : block data
c     
c     omx,omy,omz=angular velocity in body fixed frame (principal axes)
c     rotinx,y,z =rotational inertia in body fixed frame
c     
c     copyright daresbury laboratory
c     author      m. leslie jan 2004
c     amended     w.smith   jan 2005: f90 conversion
c     
c**********************************************************************

      implicit none

      logical safe,lshmov,newjob

      integer isw,imcon,idnode,mxnode,natms,ngrp,nscons,ntcons
      integer ntfree,i,j,k,igrp1,igrp2,igrp,ifre1,ifre2,jr
      integer id,ifre,jrs,idum,ig

      real(8) engke,engrot,tolnce,tstep,vircom,vircon
      real(8) tqx,tqy,tqz,ftx,fty,ftz
      real(8) vaa,vbb,vcc,engtrn,fmx,fmy,fmz
      real(8) qt0,qt1,qt2,qt3,opx,opy,opz,engfke

      integer, parameter :: nnn=6
      real(8), parameter :: pt5=0.5d0

      integer fail(nnn)
      real(8) rot(9),strkin(9),strgrp(9)

      real(8), allocatable :: xxt(:),yyt(:),zzt(:)
      real(8), allocatable :: dxx(:),dyy(:),dzz(:)
      real(8), allocatable :: txx(:),tyy(:),tzz(:)
      real(8), allocatable :: dxt(:),dyt(:),dzt(:)
      real(8), allocatable :: dtx(:),dty(:),dtz(:)
      real(8), allocatable :: p0(:),p1(:),p2(:),p3(:)
      
      save newjob,p0,p1,p2,p3
      
      data newjob/.true./
      
c     free atom block indices
      
      ifre1=(idnode*ntfree)/mxnode+1
      ifre2=((idnode+1)*ntfree)/mxnode
      
c     group block indices

      igrp1=(idnode*ngrp)/mxnode+1
      igrp2=((idnode+1)*ngrp)/mxnode

c     check work arrays are large enough

      safe=(igrp2-igrp1+1.le.msgrp) 
      if(mxnode.gt.1) call gstate(safe)
      if(.not.safe)then 
        igrp=igrp2-igrp1+1
        call gimax(igrp,1,idum)
        if(idnode.eq.0) write(nrite,*) ' make msgrp >=',igrp
        call  error(idnode,506)
      endif
      if(ntcons.gt.0)safe=.false.

c     allocate working arrays

      do i=1,nnn
        fail(i)=0
      enddo

      allocate(dtx(mxatms),dty(mxatms),dtz(mxatms),stat=fail(1))
      if(ntcons.gt.0)then

        allocate(dxx(mxcons),dyy(mxcons),dzz(mxcons),stat=fail(2))
        allocate(txx(mxatms),tyy(mxatms),tzz(mxatms),stat=fail(3))
        allocate(xxt(mxatms),yyt(mxatms),zzt(mxatms),stat=fail(4))
        allocate(dxt(mxcons),dyt(mxcons),dzt(mxcons),stat=fail(5))

      endif
      if(newjob)then
        allocate(p0(mxgrp),p1(mxgrp),p2(mxgrp),p3(mxgrp),stat=fail(6))
      endif
      do i=1,nnn
        if(fail(i).gt.0)call error(idnode,2100)
      enddo

      newjob=.false.
      
      if(ntcons.gt.0)then

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

      endif

c     atom displacement from rigid body centre of mass

      jr=0
      do ig=igrp1,igrp2
        
        id=lstgtp(ig)
        
        do j=1,numgsit(id)
          
          jr=jr+1
          i=lstrgd(jr)
          
          dtx(jr)=xxx(i)-gcmx(ig)
          dty(jr)=yyy(i)-gcmy(ig)
          dtz(jr)=zzz(i)-gcmz(ig)
          
        enddo
        
      enddo

c     periodic boundary condition for displacement vectors

      call images(imcon,0,1,jr,cell,dtx,dty,dtz)
      
c     calculate quaternion momenta at start of time step

      if(isw.eq.1)then

        do ig=igrp1,igrp2
          
          id=lstgtp(ig)
          opx=omx(ig)*rotinx(id,1)
          opy=omy(ig)*rotiny(id,1)
          opz=omz(ig)*rotinz(id,1)
          p0(ig)=2.0d0*(-q1(ig)*opx-q2(ig)*opy-q3(ig)*opz)
          p1(ig)=2.0d0*( q0(ig)*opx-q3(ig)*opy+q2(ig)*opz)
          p2(ig)=2.0d0*( q3(ig)*opx+q0(ig)*opy-q1(ig)*opz)
          p3(ig)=2.0d0*(-q2(ig)*opx+q1(ig)*opy+q0(ig)*opz)
          
        enddo

      endif

c     update free atom velocities 1/2 time step first and second stages

      do ifre=ifre1,ifre2

        i=lstfre(ifre)
        vxx(i)=vxx(i)+pt5*tstep*rmass(i)*fxx(i)
        vyy(i)=vyy(i)+pt5*tstep*rmass(i)*fyy(i)
        vzz(i)=vzz(i)+pt5*tstep*rmass(i)*fzz(i)
        
      enddo

c     *************  Rigid body motion ****************************
c     operations common to first and second stages

      jr=0
      do ig=igrp1,igrp2
        
c     calculate com force arrays 
        
        jrs=jr
        fmx=0.d0
        fmy=0.d0
        fmz=0.d0
        id=lstgtp(ig)
        do j=1,numgsit(id)
          
          jr=jr+1
          i=lstrgd(jr)
          fmx=fmx+fxx(i)
          fmy=fmy+fyy(i)
          fmz=fmz+fzz(i)

        enddo

c     current rotation matrix 
        
        call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)

c     calculate torque in principal frame

        jr=jrs
        ftx=0.d0
        fty=0.d0
        ftz=0.d0
        do j=1,numgsit(id)
          jr=jr+1
          i=lstrgd(jr)
          ftx=ftx+dty(jr)*fzz(i)-dtz(jr)*fyy(i)
          fty=fty+dtz(jr)*fxx(i)-dtx(jr)*fzz(i)
          ftz=ftz+dtx(jr)*fyy(i)-dty(jr)*fxx(i)
        enddo

        tqx=ftx*rot(1)+fty*rot(4)+ftz*rot(7)
        tqy=ftx*rot(2)+fty*rot(5)+ftz*rot(8)
        tqz=ftx*rot(3)+fty*rot(6)+ftz*rot(9)

c     calculate quaternion torques

        qt0=2.0d0*(-q1(ig)*tqx-q2(ig)*tqy-q3(ig)*tqz)
        qt1=2.0d0*( q0(ig)*tqx-q3(ig)*tqy+q2(ig)*tqz)
        qt2=2.0d0*( q3(ig)*tqx+q0(ig)*tqy-q1(ig)*tqz)
        qt3=2.0d0*(-q2(ig)*tqx+q1(ig)*tqy+q0(ig)*tqz)

c     update quaternion momenta by 1/2 time step
        
        p0(ig)=p0(ig)+qt0*pt5*tstep
        p1(ig)=p1(ig)+qt1*pt5*tstep
        p2(ig)=p2(ig)+qt2*pt5*tstep
        p3(ig)=p3(ig)+qt3*pt5*tstep

c     update centre of mass velocity by 1/2 time step

        gvxx(ig)=gvxx(ig)+fmx*pt5*tstep/gmass(id)
        gvyy(ig)=gvyy(ig)+fmy*pt5*tstep/gmass(id)
        gvzz(ig)=gvzz(ig)+fmz*pt5*tstep/gmass(id)

      enddo

c     merge centre of mass velocities from all nodes

      if(mxnode.gt.1)call merge
     x  (idnode,mxnode,ngrp,mxbuff,gvxx,gvyy,gvzz,buffer)

c     first stage of velocity verlet algorithm

      if(isw.eq.1)then

c     move centre of mass by full time step

        do ig=igrp1,igrp2

          gcmx(ig)=gcmx(ig)+tstep*gvxx(ig)
          gcmy(ig)=gcmy(ig)+tstep*gvyy(ig)
          gcmz(ig)=gcmz(ig)+tstep*gvzz(ig)

        enddo

c     merge centre of mass position from all nodes

        if(mxnode.gt.1)call merge
     x    (idnode,mxnode,ngrp,mxbuff,gcmx,gcmy,gcmz,buffer)

c     rotate rigid groups: nosquish algorithm

        jr=0
        do ig=igrp1,igrp2
          
          call nosquish(ig,tstep,q0,q1,q2,q3,p0,p1,p2,p3)
          
        enddo

c     new atomic positions for atoms in rigid bodies-relative to c.o.m
        
        k=0
        do ig=igrp1,igrp2

          call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)

          id=lstgtp(ig)

          do j=1,numgsit(id)
            
            k=k+1
            i=lstme(k)
            xxx(i)=rot(1)*gxx(id,j)+rot(2)*gyy(id,j)+rot(3)*gzz(id,j)
     x        +gcmx(ig)
            yyy(i)=rot(4)*gxx(id,j)+rot(5)*gyy(id,j)+rot(6)*gzz(id,j)
     x        +gcmy(ig)
            zzz(i)=rot(7)*gxx(id,j)+rot(8)*gyy(id,j)+rot(9)*gzz(id,j)
     x        +gcmz(ig)

          enddo
          
        enddo

c     update positions of free particles
        
        do ifre=ifre1,ifre2
          
          i=lstfre(ifre)
          xxx(i)=xxx(i)+tstep*vxx(i)
          yyy(i)=yyy(i)+tstep*vyy(i)
          zzz(i)=zzz(i)+tstep*vzz(i)
          
        enddo

c     merge atom positions

        if(mxnode.gt.1)call merge1
     x    (idnode,mxnode,natms,lstme,xxx,yyy,zzz,buffer)

c     apply shake corrections to bond constraints

        if(ntcons.gt.0)then
          
          if(mxnode.gt.1)call merge1
     x      (idnode,mxnode,natms,lstme,vxx,vyy,vzz,buffer)

          safe=.false.
          call rdrattle_r
     x      (safe,lshmov,idnode,imcon,mxnode,natms,nscons,
     x      tolnce,tstep,vircon,dxx,dyy,dzz,dxt,dyt,dzt,
     x      txx,tyy,tzz,xxt,yyt,zzt,strcns)
          if(.not.safe)return

        endif

c     end of first stage of velocity verlet algorithm

      else

c     second stage of velocity verlet algorithm

        jr=0
        do ig=igrp1,igrp2

c     new angular momenta and velocities
          
          opx=pt5*(-q1(ig)*p0(ig)+q0(ig)*p1(ig)+
     x      q3(ig)*p2(ig)-q2(ig)*p3(ig))
          opy=pt5*(-q2(ig)*p0(ig)-q3(ig)*p1(ig)+
     x      q0(ig)*p2(ig)+q1(ig)*p3(ig))
          opz=pt5*(-q3(ig)*p0(ig)+q2(ig)*p1(ig)-
     x      q1(ig)*p2(ig)+q0(ig)*p3(ig))

          id=lstgtp(ig)

          omx(ig)=opx*rotinx(id,2)
          omy(ig)=opy*rotiny(id,2)
          omz(ig)=opz*rotinz(id,2)
          
          call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)

          do j=1,numgsit(id)
            
            jr=jr+1
            i=lstrgd(jr)
            
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

c     merge velocities from all nodes

        if(mxnode.gt.1)call merge1
     x    (idnode,mxnode,natms,lstme,vxx,vyy,vzz,buffer)

c     correct constraint bond velocities using rattle

        if(ntcons.gt.0)then

          safe=.false.
          call rdrattle_v
     x      (safe,idnode,mxnode,natms,nscons,tolnce,tstep,
     x      dxx,dyy,dzz,txx,tyy,tzz,xxt,yyt,zzt)
          if(.not.safe)return

        endif

c     calculate rigid body contribution to stress tensor

        call bodystress(idnode,mxnode,ngrp,vircom,strbod,dtx,dty,dtz)

c     end of second stage of velocity verlet algorithm
        
      endif

c     periodic boundary condition
          
      if(isw.eq.2)then
        
        call images(imcon,idnode,mxnode,natms,cell,xxx,yyy,zzz)
        call images(imcon,idnode,mxnode,ngrp,cell,gcmx,gcmy,gcmz)
        
c     merge position data
        
        if(mxnode.gt.1)then
          
          call merge(idnode,mxnode,natms,mxbuff,xxx,yyy,zzz,buffer)
          call merge(idnode,mxnode,ngrp,mxbuff,gcmx,gcmy,gcmz,buffer)
          
        endif
        
      endif
      
c     calculate kinetic energy
      
      if(isw.eq.2)then

        engfke=getkinf(ntfree,idnode,mxnode)
        call getking(ngrp,idnode,mxnode,engtrn,engrot)
        engke=engfke+engtrn

c     kinetic contribution to stress tensor
        
        call kinstressf(ntfree,idnode,mxnode,strkin)        
        call kinstressg(ngrp,idnode,mxnode,strgrp)
        
c     add contributions to stress tensor
        
        do i=1,9
          stress(i)=stress(i)+strbod(i)+strcns(i)+strkin(i)+strgrp(i)
        enddo
        
      endif

c     deallocate working arrays

      deallocate(dtx,dty,dtz,stat=fail(1))
      if(ntcons.gt.0)then

        deallocate(dxx,dyy,dzz,txx,tyy,tzz,stat=fail(2))
        deallocate(xxt,yyt,zzt,dxt,dyt,dzt,stat=fail(3))

      endif
      
      return
      end subroutine nveqvv_1

      subroutine nvtqvv_b1
     x  (safe,lshmov,isw,imcon,idnode,mxnode,natms,ngrp,nscons,
     x  ntcons,ntfree,engke,engrot,taut,sigma,tolnce,tstep,
     x  vircom,vircon)
      
c***********************************************************************
c     
c     dlpoly routine to integrate rigid body equations of motion
c     using the symplectic no_squish quaternion algorithm of 
c     miller et al j.chem.phys 116 (2002) 8649
c     plus bond constraints - provided rigid body sites
c     and constraint sites do not coincide.
c     
c     nvt ensemble - Berendsen thermostat (n.b. not symplectic)
c     
c     parallel replicated data version : block data
c     
c     omx,omy,omz=angular velocity in body fixed frame (principal axes)
c     rotinx,y,z =rotational inertia in body fixed frame
c     
c     copyright daresbury laboratory
c     author      w.smith   mar 2005
c     
c**********************************************************************

      implicit none

      logical safe,lshmov,newjob

      integer isw,imcon,idnode,mxnode,natms,ngrp,nscons,ntcons
      integer ntfree,i,j,k,igrp1,igrp2,igrp,ifre1,ifre2,jr
      integer id,ifre,jrs,idum,ig,iatm0,iatm1

      real(8) engke,engrot,tolnce,tstep,vircom,vircon
      real(8) tqx,tqy,tqz,ftx,fty,ftz
      real(8) vaa,vbb,vcc,engtrn,fmx,fmy,fmz
      real(8) qt0,qt1,qt2,qt3,opx,opy,opz,taut,sigma,engtke
      real(8) chit,engfke

      integer, parameter :: nnn=6
      real(8), parameter :: pt5=0.5d0

      integer fail(nnn)
      real(8) rot(9),strkin(9),strgrp(9)

      real(8), allocatable :: xxt(:),yyt(:),zzt(:)
      real(8), allocatable :: dxx(:),dyy(:),dzz(:)
      real(8), allocatable :: txx(:),tyy(:),tzz(:)
      real(8), allocatable :: dxt(:),dyt(:),dzt(:)
      real(8), allocatable :: dtx(:),dty(:),dtz(:)
      real(8), allocatable :: p0(:),p1(:),p2(:),p3(:)
      
      save newjob,p0,p1,p2,p3

      data newjob/.true./
      
c     atom block indices
      
      iatm0=(idnode*natms)/mxnode+1
      iatm1=((idnode+1)*natms)/mxnode
      
c     free atom block indices
      
      ifre1=(idnode*ntfree)/mxnode+1
      ifre2=((idnode+1)*ntfree)/mxnode
      
c     group block indices

      igrp1=(idnode*ngrp)/mxnode+1
      igrp2=((idnode+1)*ngrp)/mxnode

c     check work arrays are large enough

      safe=(igrp2-igrp1+1.le.msgrp) 
      if(mxnode.gt.1) call gstate(safe)
      if(.not.safe)then 
        igrp=igrp2-igrp1+1
        call gimax(igrp,1,idum)
        if(idnode.eq.0) write(nrite,*) ' make msgrp >=',igrp
        call  error(idnode,506)
      endif
      if(ntcons.gt.0)safe=.false.

c     allocate working arrays

      do i=1,nnn
        fail(i)=0
      enddo
      allocate(dtx(mxatms),dty(mxatms),dtz(mxatms),stat=fail(1))
      if(ntcons.gt.0)then

        allocate(dxx(mxcons),dyy(mxcons),dzz(mxcons),stat=fail(2))
        allocate(txx(mxatms),tyy(mxatms),tzz(mxatms),stat=fail(3))
        allocate(xxt(mxatms),yyt(mxatms),zzt(mxatms),stat=fail(4))
        allocate(dxt(mxcons),dyt(mxcons),dzt(mxcons),stat=fail(5))

      endif
      if(newjob)then
        allocate(p0(mxgrp),p1(mxgrp),p2(mxgrp),p3(mxgrp),stat=fail(6))
      endif
      do i=1,nnn
        if(fail(i).gt.0)call error(idnode,2120)
      enddo

      newjob=.false.

      if(ntcons.gt.0)then

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

      endif

c     atom displacement from rigid body centre of mass

      jr=0
      do ig=igrp1,igrp2
        
        id=lstgtp(ig)
        
        do j=1,numgsit(id)
          
          jr=jr+1
          i=lstrgd(jr)
          
          dtx(jr)=xxx(i)-gcmx(ig)
          dty(jr)=yyy(i)-gcmy(ig)
          dtz(jr)=zzz(i)-gcmz(ig)
          
        enddo
        
      enddo

c     periodic boundary condition for displacement vectors

      call images(imcon,0,1,jr,cell,dtx,dty,dtz)
      
c     calculate quaternion momenta at start of time step

      if(isw.eq.1)then

        do ig=igrp1,igrp2
          
          id=lstgtp(ig)
          opx=omx(ig)*rotinx(id,1)
          opy=omy(ig)*rotiny(id,1)
          opz=omz(ig)*rotinz(id,1)
          p0(ig)=2.0d0*(-q1(ig)*opx-q2(ig)*opy-q3(ig)*opz)
          p1(ig)=2.0d0*( q0(ig)*opx-q3(ig)*opy+q2(ig)*opz)
          p2(ig)=2.0d0*( q3(ig)*opx+q0(ig)*opy-q1(ig)*opz)
          p3(ig)=2.0d0*(-q2(ig)*opx+q1(ig)*opy+q0(ig)*opz)
          
        enddo

      endif

c     update free atom velocities 1/2 time step first and second stages

      do ifre=ifre1,ifre2

        i=lstfre(ifre)
        vxx(i)=vxx(i)+pt5*tstep*rmass(i)*fxx(i)
        vyy(i)=vyy(i)+pt5*tstep*rmass(i)*fyy(i)
        vzz(i)=vzz(i)+pt5*tstep*rmass(i)*fzz(i)
        
      enddo

c     *************  Rigid body motion ****************************
c     operations common to first and second stages

      jr=0
      do ig=igrp1,igrp2
        
c     calculate com force arrays 
        
        jrs=jr
        fmx=0.d0
        fmy=0.d0
        fmz=0.d0
        id=lstgtp(ig)
        do j=1,numgsit(id)
          
          jr=jr+1
          i=lstrgd(jr)
          fmx=fmx+fxx(i)
          fmy=fmy+fyy(i)
          fmz=fmz+fzz(i)

        enddo

c     current rotation matrix 
        
        call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)

c     calculate torque in principal frame

        jr=jrs
        ftx=0.d0
        fty=0.d0
        ftz=0.d0
        do j=1,numgsit(id)
          jr=jr+1
          i=lstrgd(jr)
          ftx=ftx+dty(jr)*fzz(i)-dtz(jr)*fyy(i)
          fty=fty+dtz(jr)*fxx(i)-dtx(jr)*fzz(i)
          ftz=ftz+dtx(jr)*fyy(i)-dty(jr)*fxx(i)
        enddo
        tqx=ftx*rot(1)+fty*rot(4)+ftz*rot(7)
        tqy=ftx*rot(2)+fty*rot(5)+ftz*rot(8)
        tqz=ftx*rot(3)+fty*rot(6)+ftz*rot(9)

c     calculate quaternion torques

        qt0=2.0d0*(-q1(ig)*tqx-q2(ig)*tqy-q3(ig)*tqz)
        qt1=2.0d0*( q0(ig)*tqx-q3(ig)*tqy+q2(ig)*tqz)
        qt2=2.0d0*( q3(ig)*tqx+q0(ig)*tqy-q1(ig)*tqz)
        qt3=2.0d0*(-q2(ig)*tqx+q1(ig)*tqy+q0(ig)*tqz)

c     update quaternion momenta by 1/2 time step

        p0(ig)=p0(ig)+qt0*pt5*tstep
        p1(ig)=p1(ig)+qt1*pt5*tstep
        p2(ig)=p2(ig)+qt2*pt5*tstep
        p3(ig)=p3(ig)+qt3*pt5*tstep

c     update centre of mass velocity by 1/2 time step

        gvxx(ig)=gvxx(ig)+fmx*pt5*tstep/gmass(id)
        gvyy(ig)=gvyy(ig)+fmy*pt5*tstep/gmass(id)
        gvzz(ig)=gvzz(ig)+fmz*pt5*tstep/gmass(id)

      enddo

c     first stage of velocity verlet algorithm

      if(isw.eq.1)then

c     move centre of mass by full time step

        do ig=igrp1,igrp2

          gcmx(ig)=gcmx(ig)+tstep*gvxx(ig)
          gcmy(ig)=gcmy(ig)+tstep*gvyy(ig)
          gcmz(ig)=gcmz(ig)+tstep*gvzz(ig)

        enddo

c     merge group coms from all nodes

        if(mxnode.gt.1)call merge
     x    (idnode,mxnode,ngrp,mxbuff,gcmx,gcmy,gcmz,buffer)

c     rotate rigid groups: nosquish algoritm

        jr=0
        do ig=igrp1,igrp2
          call nosquish(ig,tstep,q0,q1,q2,q3,p0,p1,p2,p3)
        enddo

c     new atomic positions for atoms in rigid bodies-relative to c.o.m
        
        k=0
        do ig=igrp1,igrp2

          call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)
          
          id=lstgtp(ig)

          do j=1,numgsit(id)
            
            k=k+1
            i=lstme(k)
            xxx(i)=rot(1)*gxx(id,j)+rot(2)*gyy(id,j)+rot(3)*gzz(id,j)
     x        +gcmx(ig)
            yyy(i)=rot(4)*gxx(id,j)+rot(5)*gyy(id,j)+rot(6)*gzz(id,j)
     x        +gcmy(ig)
            zzz(i)=rot(7)*gxx(id,j)+rot(8)*gyy(id,j)+rot(9)*gzz(id,j)
     x        +gcmz(ig)

          enddo
          
        enddo

c     update positions of free particles
        
        do ifre=ifre1,ifre2
          
          i=lstfre(ifre)
          xxx(i)=xxx(i)+tstep*vxx(i)
          yyy(i)=yyy(i)+tstep*vyy(i)
          zzz(i)=zzz(i)+tstep*vzz(i)
          
        enddo

c     merge atom positions

        if(mxnode.gt.1)call merge1
     x    (idnode,mxnode,natms,lstme,xxx,yyy,zzz,buffer)

c     apply shake corrections to bond constraints

        if(ntcons.gt.0)then
          
          if(mxnode.gt.1)call merge1
     x      (idnode,mxnode,natms,lstme,vxx,vyy,vzz,buffer)

          safe=.false.
          call rdrattle_r
     x      (safe,lshmov,idnode,imcon,mxnode,natms,nscons,
     x      tolnce,tstep,vircon,dxx,dyy,dzz,dxt,dyt,dzt,
     x      txx,tyy,tzz,xxt,yyt,zzt,strcns)
          if(.not.safe)return

        endif

c     end of first stage of velocity verlet algorithm

      else

c     second stage of velocity verlet algorithm

        jr=0
        do ig=igrp1,igrp2

c     new angular momenta and velocities
          
          opx=pt5*(-q1(ig)*p0(ig)+q0(ig)*p1(ig)+
     x      q3(ig)*p2(ig)-q2(ig)*p3(ig))
          opy=pt5*(-q2(ig)*p0(ig)-q3(ig)*p1(ig)+
     x      q0(ig)*p2(ig)+q1(ig)*p3(ig))
          opz=pt5*(-q3(ig)*p0(ig)+q2(ig)*p1(ig)-
     x      q1(ig)*p2(ig)+q0(ig)*p3(ig))

          id=lstgtp(ig)

          omx(ig)=opx*rotinx(id,2)
          omy(ig)=opy*rotiny(id,2)
          omz(ig)=opz*rotinz(id,2)
          
c     new rotation matrix
          
          call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)

          do j=1,numgsit(id)
            
            jr=jr+1
            i=lstrgd(jr)
            
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

c     merge velocities from all nodes

        if(mxnode.gt.1)call merge1
     x    (idnode,mxnode,natms,lstme,vxx,vyy,vzz,buffer)

c     correct constraint bond velocities using rattle

        if(ntcons.gt.0)then

          safe=.false.
          call rdrattle_v
     x      (safe,idnode,mxnode,natms,nscons,tolnce,tstep,
     x      dxx,dyy,dzz,txx,tyy,tzz,xxt,yyt,zzt)
          if(.not.safe)return

        endif

c     calculate rigid body contribution to stress tensor

        call bodystress(idnode,mxnode,ngrp,vircom,strbod,dtx,dty,dtz)

c     end of second stage of velocity verlet algorithm
        
      endif

c     calculate kinetic energy
      
      if(isw.eq.2)then

        engfke=getkinf(ntfree,idnode,mxnode)
        call getking(ngrp,idnode,mxnode,engtrn,engrot)
        engtke=engfke+engtrn+engrot
        engke=engfke+engtrn

c     apply Berendsen thermostat - taut is the relaxation time
        
        chit=sqrt(1.d0+tstep/taut*(sigma/engtke-1.d0))
        
        engke=engke*chit**2
        engtrn=engtrn*chit**2
        engrot=engrot*chit**2
        
c     thermostat velocities 
        
        do i=iatm0,iatm1
          
          if(lstfrz(i).eq.0)then

            vxx(i)=chit*vxx(i)
            vyy(i)=chit*vyy(i)
            vzz(i)=chit*vzz(i)
            
          endif
          
        enddo

c     merge velocities from all nodes

        if(mxnode.gt.1)call merge
     x    (idnode,mxnode,natms,mxbuff,vxx,vyy,vzz,buffer)

c     thermostat rigid body velocities

        do ig=igrp1,igrp2
          
          omx(ig)=chit*omx(ig)
          omy(ig)=chit*omy(ig)
          omz(ig)=chit*omz(ig)
          gvxx(ig)=chit*gvxx(ig)
          gvyy(ig)=chit*gvyy(ig)
          gvzz(ig)=chit*gvzz(ig)

        enddo

c     kinetic contribution to stress tensor
        
        call kinstressf(ntfree,idnode,mxnode,strkin)        
        call kinstressg(ngrp,idnode,mxnode,strgrp)
        
c     add contributions to stress tensor
        
        do i=1,9
          stress(i)=stress(i)+strbod(i)+strcns(i)+strkin(i)+strgrp(i)
        enddo
        
      endif

c     periodic boundary condition
          
      if(isw.eq.2)then
        
        call images(imcon,idnode,mxnode,natms,cell,xxx,yyy,zzz)
        call images(imcon,idnode,mxnode,ngrp,cell,gcmx,gcmy,gcmz)
        
c     exchange position data
        
        if(mxnode.gt.1)then
          
          call merge(idnode,mxnode,natms,mxbuff,xxx,yyy,zzz,buffer)
          call merge(idnode,mxnode,ngrp,mxbuff,gcmx,gcmy,gcmz,buffer)
          
        endif
        
      endif
      
c     merge group velocities from all processors

      if(mxnode.gt.1)call merge
     x  (idnode,mxnode,ngrp,mxbuff,gvxx,gvyy,gvzz,buffer)

c     deallocate working arrays

      deallocate(dtx,dty,dtz,stat=fail(1))
      if(ntcons.gt.0)then

        deallocate(dxx,dyy,dzz,txx,tyy,tzz,stat=fail(2))
        deallocate(xxt,yyt,zzt,dxt,dyt,dzt,stat=fail(3))

      endif
      
      return
      end subroutine nvtqvv_b1

      subroutine nvtqvv_h1
     x  (safe,lshmov,isw,imcon,idnode,mxnode,natms,ngrp,nscons,
     x  ntcons,ntfree,ntshl,keyshl,chit,consv,conint,engke,engrot,
     x  taut,sigma,tolnce,tstep,vircom,vircon,chit_shl,sigma_shl)
      
c***********************************************************************
c     
c     dlpoly routine to integrate rigid body equations of motion
c     using the symplectic no_squish quaternion algorithm of 
c     miller et al j.chem.phys 116 (2002) 8649
c     plus bond constraints - provided rigid body sites
c     and constraint sites do not coincide.
c     
c     nvt ensemble - nose-hoover thermostat Molec Phys 87 (1996) 1117
c     
c     parallel replicated data version : block data
c     
c     omx,omy,omz=angular velocity in body fixed frame (principal axes)
c     rotinx,y,z =rotational inertia in body fixed frame
c     
c     copyright daresbury laboratory
c     author      w.smith april 2005
c     adapted     d.quigley : metadynamics
c     
c**********************************************************************

      implicit none

      integer, parameter :: nnn=6
      real(8), parameter :: pt5=0.5d0

      logical safe,lshmov,newjob
      integer isw,imcon,idnode,mxnode,natms,ngrp,nscons,ntcons,ntfree
      integer ifre1,ifre2,igrp1,igrp2,igrp,i,j,k,jr
      integer id,ig,ifre,jrs,idum
      real(8) chit,consv,conint,engke,engrot,taut,sigma,tolnce,tstep
      real(8) vircom,vircon,hstep,qmass,opx,opy,opz,fmx,fmy,fmz,engtrn
      real(8) ftx,fty,ftz,tqx,tqy,tqz,qt0,qt1,qt2,qt3,vaa,vbb,vcc
      real(8) engfke

      integer fail(nnn)
      real(8) rot(9),strkin(9),strgrp(9)

      real(8), allocatable :: xxt(:),yyt(:),zzt(:)
      real(8), allocatable :: dxx(:),dyy(:),dzz(:)
      real(8), allocatable :: txx(:),tyy(:),tzz(:)
      real(8), allocatable :: dxt(:),dyt(:),dzt(:)
      real(8), allocatable :: dtx(:),dty(:),dtz(:)
      real(8), allocatable :: p0(:),p1(:),p2(:),p3(:)
      
c     metadynamics shell thermostat variables

      integer ntshl,keyshl
      real(8) sigma_shl

      logical,save :: lfirst=.true.
      real(8)      :: chit_shl  
      real(8),save :: qmass_shl
      real(8)      :: shlke

c     end metadynamics shell thermostat variables

      save newjob,p0,p1,p2,p3

      data newjob/.true./
      
c     timestep parameters

      hstep=pt5*tstep

c     nose-hoover inertia parameter

      qmass=2.d0*sigma*taut**2

c     free atom block indices
      
      ifre1=(idnode*ntfree)/mxnode+1
      ifre2=((idnode+1)*ntfree)/mxnode
      
c     group block indices

      igrp1=(idnode*ngrp)/mxnode+1
      igrp2=((idnode+1)*ngrp)/mxnode

c     check work arrays are large enough

      safe=(igrp2-igrp1+1.le.msgrp) 
      if(mxnode.gt.1) call gstate(safe)
      if(.not.safe)then 
        igrp=igrp2-igrp1+1
        call gimax(igrp,1,idum)
        if(idnode.eq.0) write(nrite,*) ' make msgrp >=',igrp
        call  error(idnode,506)
      endif
      if(ntcons.gt.0)safe=.false.

c     allocate working arrays

      do i=1,nnn
        fail(i)=0
      enddo

      allocate(dtx(mxatms),dty(mxatms),dtz(mxatms),stat=fail(1))
      if(ntcons.gt.0)then

        allocate(dxx(mxcons),dyy(mxcons),dzz(mxcons),stat=fail(2))
        allocate(txx(mxatms),tyy(mxatms),tzz(mxatms),stat=fail(3))
        allocate(xxt(mxatms),yyt(mxatms),zzt(mxatms),stat=fail(4))
        allocate(dxt(mxcons),dyt(mxcons),dzt(mxcons),stat=fail(5))

      endif
      if(newjob)then
        allocate(p0(mxgrp),p1(mxgrp),p2(mxgrp),p3(mxgrp),stat=fail(6))
      endif
      do i=1,nnn
        if(fail(i).gt.0)call error(idnode,2140)
      enddo

      newjob=.false.

      if(lmetadyn.and.lfirst.and.(ntshl>0))then
        if(idnode.eq.0)then
          write(*,*)"Warning - Metadynamics Modification"
          write(*,*)"========================="
          write(*,*)"Coupling core-shell motion thermostat at 1 K"
        endif
        lfirst=.false.
c     use same relaxation time for global and core-shell?
        qmass_shl=2.d0*sigma_shl*taut**2
      endif

      if(ntcons.gt.0)then

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

      endif

c     atom displacement from rigid body centre of mass

      jr=0
      do ig=igrp1,igrp2
        
        id=lstgtp(ig)
        
        do j=1,numgsit(id)
          
          jr=jr+1
          i=lstrgd(jr)
          
          dtx(jr)=xxx(i)-gcmx(ig)
          dty(jr)=yyy(i)-gcmy(ig)
          dtz(jr)=zzz(i)-gcmz(ig)
          
        enddo
        
      enddo

c     periodic boundary condition for displacement vectors

      call images(imcon,0,1,jr,cell,dtx,dty,dtz)
      
c     apply thermostat for first stage

      if(isw.eq.1)then
        
        call nvtqscl
     x    (idnode,mxnode,ntfree,ngrp,engfke,engtrn,engrot,sigma,
     x    hstep,qmass,taut,chit,conint)
        
c     metadynamics shell thermostat
c     ====================================================
c     Must first merge update velocities as the core-shell
c     velocites are not distributed according to the same
c     rules.
c     ====================================================
        if(lmetadyn.and.keyshl.eq.1)then
          if(mxnode.gt.1)call merge
     x      (idnode,mxnode,natms,mxbuff,vxx,vyy,vzz,buffer)
          call nvtscale_shl
     x      (idnode,mxnode,ntshl,shlke,sigma_shl,hstep,qmass_shl,
     x      taut,chit_shl,conint)      
        endif
        
      endif
      
c     calculate quaternion momenta at start of time step
      
      if(isw.eq.1)then
        
        do ig=igrp1,igrp2
          
          id=lstgtp(ig)
          opx=omx(ig)*rotinx(id,1)
          opy=omy(ig)*rotiny(id,1)
          opz=omz(ig)*rotinz(id,1)
          p0(ig)=2.0d0*(-q1(ig)*opx-q2(ig)*opy-q3(ig)*opz)
          p1(ig)=2.0d0*( q0(ig)*opx-q3(ig)*opy+q2(ig)*opz)
          p2(ig)=2.0d0*( q3(ig)*opx+q0(ig)*opy-q1(ig)*opz)
          p3(ig)=2.0d0*(-q2(ig)*opx+q1(ig)*opy+q0(ig)*opz)
          
        enddo

      endif

c     update free atom velocities 1/2 time step first and second stages

      do ifre=ifre1,ifre2

        i=lstfre(ifre)
        vxx(i)=vxx(i)+hstep*rmass(i)*fxx(i)
        vyy(i)=vyy(i)+hstep*rmass(i)*fyy(i)
        vzz(i)=vzz(i)+hstep*rmass(i)*fzz(i)
        
      enddo

c     *************  Rigid body motion ****************************
c     operations common to first and second stages

      jr=0
      do ig=igrp1,igrp2
        
c     calculate com force arrays 
        
        jrs=jr
        fmx=0.d0
        fmy=0.d0
        fmz=0.d0
        id=lstgtp(ig)
        do j=1,numgsit(id)
          
          jr=jr+1
          i=lstrgd(jr)
          fmx=fmx+fxx(i)
          fmy=fmy+fyy(i)
          fmz=fmz+fzz(i)

        enddo

c     current rotation matrix 
        
        call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)

c     calculate torque in principal frame

        jr=jrs
        ftx=0.d0
        fty=0.d0
        ftz=0.d0
        do j=1,numgsit(id)
          jr=jr+1
          i=lstrgd(jr)
          ftx=ftx+dty(jr)*fzz(i)-dtz(jr)*fyy(i)
          fty=fty+dtz(jr)*fxx(i)-dtx(jr)*fzz(i)
          ftz=ftz+dtx(jr)*fyy(i)-dty(jr)*fxx(i)
        enddo
        tqx=ftx*rot(1)+fty*rot(4)+ftz*rot(7)
        tqy=ftx*rot(2)+fty*rot(5)+ftz*rot(8)
        tqz=ftx*rot(3)+fty*rot(6)+ftz*rot(9)

c     calculate quaternion torques

        qt0=2.0d0*(-q1(ig)*tqx-q2(ig)*tqy-q3(ig)*tqz)
        qt1=2.0d0*( q0(ig)*tqx-q3(ig)*tqy+q2(ig)*tqz)
        qt2=2.0d0*( q3(ig)*tqx+q0(ig)*tqy-q1(ig)*tqz)
        qt3=2.0d0*(-q2(ig)*tqx+q1(ig)*tqy+q0(ig)*tqz)

c     update quaternion momenta by 1/2 time step

        p0(ig)=p0(ig)+qt0*hstep
        p1(ig)=p1(ig)+qt1*hstep
        p2(ig)=p2(ig)+qt2*hstep
        p3(ig)=p3(ig)+qt3*hstep

c     update centre of mass velocity by 1/2 time step

        gvxx(ig)=gvxx(ig)+fmx*hstep/gmass(id)
        gvyy(ig)=gvyy(ig)+fmy*hstep/gmass(id)
        gvzz(ig)=gvzz(ig)+fmz*hstep/gmass(id)

      enddo

c     first stage of velocity verlet algorithm

      if(isw.eq.1)then

c     move centre of mass by full time step

        do ig=igrp1,igrp2

          gcmx(ig)=gcmx(ig)+tstep*gvxx(ig)
          gcmy(ig)=gcmy(ig)+tstep*gvyy(ig)
          gcmz(ig)=gcmz(ig)+tstep*gvzz(ig)

        enddo

c     merge group coms from all nodes

        if(mxnode.gt.1)call merge
     x    (idnode,mxnode,ngrp,mxbuff,gcmx,gcmy,gcmz,buffer)

c     rotate rigid groups: nosquish algoritm

        jr=0
        do ig=igrp1,igrp2
          call nosquish(ig,tstep,q0,q1,q2,q3,p0,p1,p2,p3)
        enddo

c     new atomic positions for atoms in rigid bodies-relative to c.o.m
        
        k=0
        do ig=igrp1,igrp2

          call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)
          
          id=lstgtp(ig)

          do j=1,numgsit(id)
            
            k=k+1
            i=lstme(k)
            xxx(i)=rot(1)*gxx(id,j)+rot(2)*gyy(id,j)+rot(3)*gzz(id,j)
     x        +gcmx(ig)
            yyy(i)=rot(4)*gxx(id,j)+rot(5)*gyy(id,j)+rot(6)*gzz(id,j)
     x        +gcmy(ig)
            zzz(i)=rot(7)*gxx(id,j)+rot(8)*gyy(id,j)+rot(9)*gzz(id,j)
     x        +gcmz(ig)

          enddo
          
        enddo

c     update positions of free particles
        
        do ifre=ifre1,ifre2
          
          i=lstfre(ifre)
          xxx(i)=xxx(i)+tstep*vxx(i)
          yyy(i)=yyy(i)+tstep*vyy(i)
          zzz(i)=zzz(i)+tstep*vzz(i)
          
        enddo

c     merge atom positions

        if(mxnode.gt.1)call merge1
     x    (idnode,mxnode,natms,lstme,xxx,yyy,zzz,buffer)
      
c     apply shake corrections to bond constraints

        if(ntcons.gt.0)then
          
          if(mxnode.gt.1)call merge1
     x      (idnode,mxnode,natms,lstme,vxx,vyy,vzz,buffer)

          safe=.false.
          call rdrattle_r
     x      (safe,lshmov,idnode,imcon,mxnode,natms,nscons,
     x      tolnce,tstep,vircon,dxx,dyy,dzz,dxt,dyt,dzt,
     x      txx,tyy,tzz,xxt,yyt,zzt,strcns)
          if(.not.safe)return

        endif

c     end of first stage of velocity verlet algorithm

      else

c     second stage of velocity verlet algorithm

        jr=0
        do ig=igrp1,igrp2

c     new angular momenta and velocities
          
          opx=pt5*(-q1(ig)*p0(ig)+q0(ig)*p1(ig)+
     x      q3(ig)*p2(ig)-q2(ig)*p3(ig))
          opy=pt5*(-q2(ig)*p0(ig)-q3(ig)*p1(ig)+
     x      q0(ig)*p2(ig)+q1(ig)*p3(ig))
          opz=pt5*(-q3(ig)*p0(ig)+q2(ig)*p1(ig)-
     x      q1(ig)*p2(ig)+q0(ig)*p3(ig))

          id=lstgtp(ig)

          omx(ig)=opx*rotinx(id,2)
          omy(ig)=opy*rotiny(id,2)
          omz(ig)=opz*rotinz(id,2)
          
c     new rotation matrix
          
          call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)

          do j=1,numgsit(id)
            
            jr=jr+1
            i=lstrgd(jr)
            
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

c     merge velocities from all nodes

        if(mxnode.gt.1)call merge1
     x    (idnode,mxnode,natms,lstme,vxx,vyy,vzz,buffer)

c     correct constraint bond velocities using rattle

        if(ntcons.gt.0)then

          safe=.false.
          call rdrattle_v
     x      (safe,idnode,mxnode,natms,nscons,tolnce,tstep,
     x      dxx,dyy,dzz,txx,tyy,tzz,xxt,yyt,zzt)
          if(.not.safe)return

        endif

c     calculate rigid body contribution to stress tensor

        call bodystress(idnode,mxnode,ngrp,vircom,strbod,dtx,dty,dtz)

c     apply thermostat for second stage and calculate kinetic energy
      
        call nvtqscl
     x    (idnode,mxnode,ntfree,ngrp,engfke,engtrn,engrot,sigma,
     x    hstep,qmass,taut,chit,conint)

        engke=engfke+engtrn

c     merge velocities from all nodes

        if(mxnode.gt.1)call merge1
     x    (idnode,mxnode,natms,lstme,vxx,vyy,vzz,buffer)

c     metadynamics shell thermostat
c     ====================================================
c     Must first merge update velocities as the core-shell
c     velocites are not distributed according to the same
c     rules.
c     ====================================================
        if(lmetadyn.and.keyshl.eq.1)then
          call nvtscale_shl
     x      (idnode,mxnode,ntshl,shlke,sigma_shl,hstep,qmass_shl,
     x      taut,chit_shl,conint)      
        endif
        
c     conserved quantity less kinetic and potential energy terms
        
        consv=conint+0.5d0*qmass*chit**2

c     metadynamics shell thermostat

        if(lmetadyn.and.keyshl.eq.1)then
          consv=consv+0.5d0*qmass_shl*chit_shl**2
        endif

c     kinetic contribution to stress tensor
        
        call kinstressf(ntfree,idnode,mxnode,strkin)        
        call kinstressg(ngrp,idnode,mxnode,strgrp)
        
c     add contributions to stress tensor
        
        do i=1,9
          stress(i)=stress(i)+strbod(i)+strcns(i)+strkin(i)+strgrp(i)
        enddo
        
c     end of second stage of velocity verlet algorithm
        
      endif

c     periodic boundary condition
          
      if(isw.eq.2)then
        
        call images(imcon,idnode,mxnode,natms,cell,xxx,yyy,zzz)
        call images(imcon,idnode,mxnode,ngrp,cell,gcmx,gcmy,gcmz)
        
c     merge position data
        
        if(mxnode.gt.1)then
          
          call merge(idnode,mxnode,natms,mxbuff,xxx,yyy,zzz,buffer)
          call merge(idnode,mxnode,ngrp,mxbuff,gcmx,gcmy,gcmz,buffer)
          
        endif
        
      endif
      
c     merge group velocities from all processors

      if(mxnode.gt.1)call merge
     x  (idnode,mxnode,ngrp,mxbuff,gvxx,gvyy,gvzz,buffer)

c     deallocate working arrays

      deallocate(dtx,dty,dtz,stat=fail(1))
      if(ntcons.gt.0)then

        deallocate(dxx,dyy,dzz,txx,tyy,tzz,stat=fail(2))
        deallocate(xxt,yyt,zzt,dxt,dyt,dzt,stat=fail(3))

      endif
      
      return
      end subroutine nvtqvv_h1

      subroutine nptqvv_b1
     x  (safe,lshmov,isw,imcon,idnode,mxnode,natms,ngrp,nscons,
     x  ntcons,ntfree,ntpatm,engke,engrot,press,taut,taup,sigma,
     x  tolnce,tstep,vircom,vircon,elrc,virlrc,virtot,volm)
      
c***********************************************************************
c     
c     dlpoly routine to integrate rigid body equations of motion
c     using the symplectic no_squish quaternion algorithm of 
c     miller et al j.chem.phys 116 (2002) 8649
c     plus bond constraints - provided rigid body sites
c     and constraint sites do not coincide.
c     
c     npt ensemble - Berendsen thermostat and barostat 
c     (n.b. not symplectic)
c     
c     isothermal compressibility (beta) set to that of liquid water
c     = 0.007372 dlpoly units
c     
c     parallel replicated data version : block data
c     
c     omx,omy,omz=angular velocity in body fixed frame (principal axes)
c     rotinx,y,z =rotational inertia in body fixed frame
c     
c     copyright daresbury laboratory
c     author      w.smith april 2005
c     
c**********************************************************************

      implicit none

      logical safe,lshmov,newjob

      integer isw,imcon,idnode,mxnode,natms,ngrp,nscons,ntcons
      integer ntfree,i,j,k,igrp1,igrp2,igrp,ifre1,ifre2,jr,kk
      integer id,ifre,jrs,idum,ig,iatm0,iatm1,ntpatm,iter,mxiter

      real(8) engke,engrot,tolnce,tstep,vircom,vircon
      real(8) tqx,tqy,tqz,ftx,fty,ftz
      real(8) vaa,vbb,vcc,engtrn,fmx,fmy,fmz
      real(8) qt0,qt1,qt2,qt3,opx,opy,opz,taut,sigma,engtke
      real(8) chit,chip,beta,volm,volm0,elrc,elrc0,virlrc,virlrc0
      real(8) virtot,psyst,press,taup,scale,engfke

      integer, parameter :: nnn=11
      real(8), parameter :: pt5=0.5d0

      integer fail(nnn)
      real(8) rot(9),strkin(9),strgrp(9),uni(9)

      real(8), allocatable :: dens0(:)
      real(8), allocatable :: xxo(:),yyo(:),zzo(:)
      real(8), allocatable :: vxo(:),vyo(:),vzo(:)
      real(8), allocatable :: xxt(:),yyt(:),zzt(:)
      real(8), allocatable :: gxo(:),gyo(:),gzo(:)
      real(8), allocatable :: dxx(:),dyy(:),dzz(:)
      real(8), allocatable :: txx(:),tyy(:),tzz(:)
      real(8), allocatable :: dxt(:),dyt(:),dzt(:)
      real(8), allocatable :: dtx(:),dty(:),dtz(:)
      real(8), allocatable :: p0(:),p1(:),p2(:),p3(:)
      real(8), allocatable :: b0(:),b1(:),b2(:),b3(:)
      real(8), allocatable :: c0(:),c1(:),c2(:),c3(:)
      
      save newjob,volm0,elrc0,virlrc0,iatm0,iatm1,dens0
      save p0,p1,p2,p3,ifre1,ifre2,igrp1,igrp2

      data uni/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/
      data newjob/.true./,beta/7.3728d-3/
      
      safe=.true.

      if(newjob)then
        
        fail(1)=0
        allocate(dens0(mxatyp),stat=fail(1))
        if(fail(1).gt.0)call error(idnode,2160)

c     store initial values of volume and long range corrections
        
        volm0=volm
        elrc0=elrc
        virlrc0=virlrc
        do i=1,ntpatm
          dens0(i)=dens(i)
        enddo
        
c     atom block indices
        
        iatm0=(idnode*natms)/mxnode+1
        iatm1=((idnode+1)*natms)/mxnode
        
c     free atom block indices
        
        ifre1=(idnode*ntfree)/mxnode+1
        ifre2=((idnode+1)*ntfree)/mxnode
        
c     group block indices
        
        igrp1=(idnode*ngrp)/mxnode+1
        igrp2=((idnode+1)*ngrp)/mxnode
        
c     check work arrays are large enough
        
        safe=(igrp2-igrp1+1.le.msgrp) 
        if(mxnode.gt.1) call gstate(safe)
        if(.not.safe)then 
          igrp=igrp2-igrp1+1
          call gimax(igrp,1,idum)
          if(idnode.eq.0) write(nrite,*) ' make msgrp >=',igrp
          call  error(idnode,506)
        endif

      endif

c     allocate working arrays

      do i=1,nnn
        fail(i)=0
      enddo
      allocate(dtx(mxatms),dty(mxatms),dtz(mxatms),stat=fail(1))
      if(ntcons.gt.0)then

        allocate(dxx(mxcons),dyy(mxcons),dzz(mxcons),stat=fail(2))
        allocate(txx(mxatms),tyy(mxatms),tzz(mxatms),stat=fail(3))
        allocate(xxt(mxatms),yyt(mxatms),zzt(mxatms),stat=fail(4))
        allocate(dxt(mxcons),dyt(mxcons),dzt(mxcons),stat=fail(5))

        if(isw.eq.1)then

          allocate(gxo(mxgrp),gyo(mxgrp),gzo(mxgrp),stat=fail(6))
          allocate(b0(mxgrp),b1(mxgrp),b2(mxgrp),b3(mxgrp),stat=fail(7))
          allocate(c0(mxgrp),c1(mxgrp),c2(mxgrp),c3(mxgrp),stat=fail(8))
          allocate(xxo(mxatms),yyo(mxatms),zzo(mxatms),stat=fail(9))
          allocate(vxo(mxatms),vyo(mxatms),vzo(mxatms),stat=fail(10))

        endif

      endif
      if(newjob)then         
        allocate(p0(mxgrp),p1(mxgrp),p2(mxgrp),p3(mxgrp),stat=fail(11))
      endif
      do i=1,nnn
        if(fail(i).gt.0)call error(idnode,2170)
      enddo

      newjob=.false.
      if(ntcons.gt.0)safe=.false.

c     atom displacement from rigid body centre of mass

      jr=0
      do ig=igrp1,igrp2
        
        id=lstgtp(ig)
        
        do j=1,numgsit(id)
          
          jr=jr+1
          i=lstrgd(jr)
          
          dtx(jr)=xxx(i)-gcmx(ig)
          dty(jr)=yyy(i)-gcmy(ig)
          dtz(jr)=zzz(i)-gcmz(ig)
          
        enddo
        
      enddo

c     periodic boundary condition for displacement vectors

      call images(imcon,0,1,jr,cell,dtx,dty,dtz)
      
c     construct current bond vectors
      
      if(ntcons.gt.0)then

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

      endif

      if(isw.eq.1)then

c     calculate kinetic energy

        engfke=getkinf(ntfree,idnode,mxnode)
        call getking(ngrp,idnode,mxnode,engtrn,engrot)
        engke=engfke+engtrn

c     calculate quaternion momenta at start of time step

        do ig=igrp1,igrp2
          
          id=lstgtp(ig)
          opx=omx(ig)*rotinx(id,1)
          opy=omy(ig)*rotiny(id,1)
          opz=omz(ig)*rotinz(id,1)
          p0(ig)=2.0d0*(-q1(ig)*opx-q2(ig)*opy-q3(ig)*opz)
          p1(ig)=2.0d0*( q0(ig)*opx-q3(ig)*opy+q2(ig)*opz)
          p2(ig)=2.0d0*( q3(ig)*opx+q0(ig)*opy-q1(ig)*opz)
          p3(ig)=2.0d0*(-q2(ig)*opx+q1(ig)*opy+q0(ig)*opz)
          
        enddo

      endif

c     update free atom velocities 1/2 time step first and second stages

      do ifre=ifre1,ifre2

        i=lstfre(ifre)
        vxx(i)=vxx(i)+pt5*tstep*rmass(i)*fxx(i)
        vyy(i)=vyy(i)+pt5*tstep*rmass(i)*fyy(i)
        vzz(i)=vzz(i)+pt5*tstep*rmass(i)*fzz(i)
        
      enddo

c     rigid body motion for first and second stages

      jr=0
      do ig=igrp1,igrp2
        
c     calculate com force arrays 
        
        jrs=jr
        fmx=0.d0
        fmy=0.d0
        fmz=0.d0
        id=lstgtp(ig)
        do j=1,numgsit(id)
          
          jr=jr+1
          i=lstrgd(jr)
          fmx=fmx+fxx(i)
          fmy=fmy+fyy(i)
          fmz=fmz+fzz(i)

        enddo

c     current rotation matrix 
        
        call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)

c     calculate torque in principal frame

        jr=jrs
        ftx=0.d0
        fty=0.d0
        ftz=0.d0
        do j=1,numgsit(id)
          jr=jr+1
          i=lstrgd(jr)
          ftx=ftx+dty(jr)*fzz(i)-dtz(jr)*fyy(i)
          fty=fty+dtz(jr)*fxx(i)-dtx(jr)*fzz(i)
          ftz=ftz+dtx(jr)*fyy(i)-dty(jr)*fxx(i)
        enddo
        tqx=ftx*rot(1)+fty*rot(4)+ftz*rot(7)
        tqy=ftx*rot(2)+fty*rot(5)+ftz*rot(8)
        tqz=ftx*rot(3)+fty*rot(6)+ftz*rot(9)

c     calculate quaternion torques

        qt0=2.0d0*(-q1(ig)*tqx-q2(ig)*tqy-q3(ig)*tqz)
        qt1=2.0d0*( q0(ig)*tqx-q3(ig)*tqy+q2(ig)*tqz)
        qt2=2.0d0*( q3(ig)*tqx+q0(ig)*tqy-q1(ig)*tqz)
        qt3=2.0d0*(-q2(ig)*tqx+q1(ig)*tqy+q0(ig)*tqz)

c     update quaternion momenta by 1/2 time step

        p0(ig)=p0(ig)+qt0*pt5*tstep
        p1(ig)=p1(ig)+qt1*pt5*tstep
        p2(ig)=p2(ig)+qt2*pt5*tstep
        p3(ig)=p3(ig)+qt3*pt5*tstep

c     update centre of mass velocity by 1/2 time step

        gvxx(ig)=gvxx(ig)+fmx*pt5*tstep/gmass(id)
        gvyy(ig)=gvyy(ig)+fmy*pt5*tstep/gmass(id)
        gvzz(ig)=gvzz(ig)+fmz*pt5*tstep/gmass(id)

      enddo

c     first stage of velocity verlet algorithm

      if(isw.eq.1)then

c     store current integration variables
        
        if(ntcons.gt.0)then

          do i=1,natms
            
            xxo(i)=xxx(i)
            yyo(i)=yyy(i)
            zzo(i)=zzz(i)
            vxo(i)=vxx(i)
            vyo(i)=vyy(i)
            vzo(i)=vzz(i)

          enddo
          do ig=1,ngrp
            
            b0(ig)=q0(ig)
            b1(ig)=q1(ig)
            b2(ig)=q2(ig)
            b3(ig)=q3(ig)
            c0(ig)=p0(ig)
            c1(ig)=p1(ig)
            c2(ig)=p2(ig)
            c3(ig)=p3(ig)
            gxo(ig)=gcmx(ig)
            gyo(ig)=gcmy(ig)
            gzo(ig)=gcmz(ig)
            
          enddo

        endif

c     iteration required if ntcons > 0

        mxiter=1
        if(ntcons.gt.0)mxiter=2
        do iter=1,mxiter

          scale=1.d0

          if(iter.eq.mxiter)then
            
c     calculate system pressure
            
            psyst=(2.d0*engke-virtot-vircon-vircom)/(3.d0*volm)
            
c     apply Berendsen barostat
            
            chip=1.d0+beta*tstep*(psyst-press)/taup
            scale=chip**(1.d0/3.d0)
            volm=chip*volm

c     reset cell parameters for new volume
            
            do i=1,9
              cell(i)=scale*cell(i)
            enddo

          endif

c     update centre of mass position

          do ig=igrp1,igrp2

            gcmx(ig)=scale*gcmx(ig)+tstep*gvxx(ig)
            gcmy(ig)=scale*gcmy(ig)+tstep*gvyy(ig)
            gcmz(ig)=scale*gcmz(ig)+tstep*gvzz(ig)

          enddo

c     merge group coms from all nodes

          if(mxnode.gt.1)call merge
     x      (idnode,mxnode,ngrp,mxbuff,gcmx,gcmy,gcmz,buffer)

c     rotate rigid groups: nosquish algorithm

          jr=0
          do ig=igrp1,igrp2
            call nosquish(ig,tstep,q0,q1,q2,q3,p0,p1,p2,p3)
          enddo

c     new atomic positions for atoms in rigid bodies-relative to c.o.m
          
          k=0
          do ig=igrp1,igrp2
            
            call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)
            
            id=lstgtp(ig)
            
            do j=1,numgsit(id)
              
              k=k+1
              i=lstme(k)
              xxx(i)=gcmx(ig)+
     x          rot(1)*gxx(id,j)+rot(2)*gyy(id,j)+rot(3)*gzz(id,j)
              yyy(i)=gcmy(ig)+
     x          rot(4)*gxx(id,j)+rot(5)*gyy(id,j)+rot(6)*gzz(id,j)
              zzz(i)=gcmz(ig)+
     x          rot(7)*gxx(id,j)+rot(8)*gyy(id,j)+rot(9)*gzz(id,j)
              
            enddo
            
          enddo
          
c     update positions of free particles
          
          do ifre=ifre1,ifre2
            
            i=lstfre(ifre)
            xxx(i)=scale*xxx(i)+tstep*vxx(i)
            yyy(i)=scale*yyy(i)+tstep*vyy(i)
            zzz(i)=scale*zzz(i)+tstep*vzz(i)
            
          enddo
          
c     merge position data

          if(mxnode.gt.1)call merge1
     x      (idnode,mxnode,natms,lstme,xxx,yyy,zzz,buffer)
          
c     apply shake corrections to bond constraints
          
          if(ntcons.gt.0)then
            
            if(mxnode.gt.1)call merge1
     x        (idnode,mxnode,natms,lstme,vxx,vyy,vzz,buffer)

            safe=.false.
            call rdrattle_r
     x        (safe,lshmov,idnode,imcon,mxnode,natms,nscons,
     x        tolnce,tstep,vircon,dxx,dyy,dzz,dxt,dyt,dzt,
     x        txx,tyy,tzz,xxt,yyt,zzt,strcns)
            if(.not.safe)return

          endif

          if(iter.lt.mxiter)then

            do i=1,natms

              xxx(i)=xxo(i)
              yyy(i)=yyo(i)
              zzz(i)=zzo(i)
              vxx(i)=vxo(i)
              vyy(i)=vyo(i)
              vzz(i)=vzo(i)

            enddo
            do ig=1,ngrp

              q0(ig)=b0(ig)
              q1(ig)=b1(ig)
              q2(ig)=b2(ig)
              q3(ig)=b3(ig)
              p0(ig)=c0(ig)
              p1(ig)=c1(ig)
              p2(ig)=c2(ig)
              p3(ig)=c3(ig)
              gcmx(ig)=gxo(ig)
              gcmy(ig)=gyo(ig)
              gcmz(ig)=gzo(ig)

            enddo

          endif

        enddo

c     adjust long range corrections and number density
        
        elrc=elrc0*(volm0/volm)
        virlrc=virlrc0*(volm0/volm)
        
        do kk=1,ntpatm
          dens(kk)=dens0(kk)*(volm0/volm)
        enddo

c     construct scaling tensor for tethered bonds

        do i=1,9
          eta(i)=scale*uni(i)
        enddo

c     end of first stage of velocity verlet algorithm

      else

c     second stage of velocity verlet algorithm

        jr=0
        do ig=igrp1,igrp2

c     new angular momenta and velocities
          
          opx=pt5*(-q1(ig)*p0(ig)+q0(ig)*p1(ig)+
     x      q3(ig)*p2(ig)-q2(ig)*p3(ig))
          opy=pt5*(-q2(ig)*p0(ig)-q3(ig)*p1(ig)+
     x      q0(ig)*p2(ig)+q1(ig)*p3(ig))
          opz=pt5*(-q3(ig)*p0(ig)+q2(ig)*p1(ig)-
     x      q1(ig)*p2(ig)+q0(ig)*p3(ig))

          id=lstgtp(ig)

          omx(ig)=opx*rotinx(id,2)
          omy(ig)=opy*rotiny(id,2)
          omz(ig)=opz*rotinz(id,2)
          
c     new rotation matrix
          
          call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)

          do j=1,numgsit(id)
            
            jr=jr+1
            i=lstrgd(jr)
            
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

c     merge velocities and forces from all nodes

        if(mxnode.gt.1)call merge1
     x    (idnode,mxnode,natms,lstme,vxx,vyy,vzz,buffer)

c     correct constraint bond velocities using rattle

        if(ntcons.gt.0)then

          safe=.false.
          call rdrattle_v
     x      (safe,idnode,mxnode,natms,nscons,tolnce,tstep,
     x      dxx,dyy,dzz,txx,tyy,tzz,xxt,yyt,zzt)
          if(.not.safe)return

        endif

c     calculate rigid body contribution to stress tensor

        call bodystress(idnode,mxnode,ngrp,vircom,strbod,dtx,dty,dtz)

c     calculate kinetic energy
        
        if(isw.eq.2)then

          engfke=getkinf(ntfree,idnode,mxnode)
          call getking(ngrp,idnode,mxnode,engtrn,engrot)
          engtke=engfke+engtrn+engrot
          engke=engfke+engtrn

c     apply Berendsen thermostat - taut is the relaxation time
          
          chit=sqrt(1.d0+tstep/taut*(sigma/engtke-1.d0))
          
          engke=engke*chit**2
          engtrn=engtrn*chit**2
          engrot=engrot*chit**2
          
c     thermostat velocities 
          
          do i=iatm0,iatm1
            
            if(lstfrz(i).eq.0)then

              vxx(i)=chit*vxx(i)
              vyy(i)=chit*vyy(i)
              vzz(i)=chit*vzz(i)
              
            endif
            
          enddo

c     merge velocities from all nodes

          if(mxnode.gt.1)call merge
     x      (idnode,mxnode,natms,mxbuff,vxx,vyy,vzz,buffer)

c     thermostat rigid body velocities

          do ig=igrp1,igrp2
            
            omx(ig)=chit*omx(ig)
            omy(ig)=chit*omy(ig)
            omz(ig)=chit*omz(ig)
            gvxx(ig)=chit*gvxx(ig)
            gvyy(ig)=chit*gvyy(ig)
            gvzz(ig)=chit*gvzz(ig)

          enddo

c     merge group velocities from all processors

          if(mxnode.gt.1)call merge
     x      (idnode,mxnode,ngrp,mxbuff,gvxx,gvyy,gvzz,buffer)

c     kinetic contribution to stress tensor
          
          call kinstressf(ntfree,idnode,mxnode,strkin)        
          call kinstressg(ngrp,idnode,mxnode,strgrp)
          
c     add contributions to stress tensor
          
          do i=1,9
            stress(i)=stress(i)+strbod(i)+strcns(i)+strkin(i)+strgrp(i)
          enddo
          
        endif
        
c     end of second stage of velocity verlet algorithm
        
      endif

c     periodic boundary condition
          
      if(isw.eq.2)then
        
        call images(imcon,idnode,mxnode,natms,cell,xxx,yyy,zzz)
        call images(imcon,idnode,mxnode,ngrp,cell,gcmx,gcmy,gcmz)
        
c     merge position data
        
        if(mxnode.gt.1)then
          
          call merge(idnode,mxnode,natms,mxbuff,xxx,yyy,zzz,buffer)
          call merge(idnode,mxnode,ngrp,mxbuff,gcmx,gcmy,gcmz,buffer)
          
        endif
        
      endif
      
c     deallocate working arrays

      deallocate(dtx,dty,dtz,stat=fail(1))
      if(ntcons.gt.0)then

        deallocate(dxx,dyy,dzz,txx,tyy,tzz,stat=fail(2))
        deallocate(xxt,yyt,zzt,dxt,dyt,dzt,stat=fail(3))
        
        if(isw.eq.1)then
          
          deallocate(vxo,vyo,vzo,b0,b1,b2,b3,stat=fail(4))
          deallocate(xxo,yyo,zzo,c0,c1,c2,c3,stat=fail(5))
          deallocate(gxo,gyo,gzo,stat=fail(6))
          
        endif
        
      endif
      
      return
      end subroutine nptqvv_b1

      subroutine nptqvv_h1
     x  (safe,lshmov,isw,idnode,mxnode,natms,imcon,ngrp,nscons,
     x  ntcons,ntpatm,ntfree,ntshl,keyshl,tstep,taut,taup,sigma,
     x  temp,chip,chit,consv,conint,engke,engrot,elrc,tolnce,
     x  vircon,virtot,virlrc,vircom,volm,press,chit_shl,sigma_shl)
      
c***********************************************************************
c     
c     dlpoly routine to integrate rigid body equations of motion
c     using the symplectic no_squish quaternion algorithm of 
c     miller et al j.chem.phys 116 (2002) 8649
c     plus bond constraints - provided rigid body sites
c     and constraint sites do not coincide.
c     
c     npt ensemble - nose-hoover thermostat Molec Phys 87 (1996) 1117
c     
c     parallel replicated data version : block data
c     
c     omx,omy,omz=angular velocity in body fixed frame (principal axes)
c     rotinx,y,z =rotational inertia in body fixed frame
c     
c     copyright daresbury laboratory
c     author      w.smith may 2005
c     adapted     d.quigley : metadynamics
c     
c**********************************************************************

      implicit none

      integer, parameter :: nnn=12
      integer, parameter :: ncyc=5
      real(8), parameter :: pt5=0.5d0

      logical safe,lshmov,newjob
      integer isw,imcon,idnode,mxnode,natms,ngrp,nscons,ntcons,ntfree
      integer ifre1,ifre2,igrp1,igrp2,igrp,i,j,k,jr
      integer id,ig,ifre,jrs,idum,mxiter,iter,ntpatm,kk,icyc
      real(8) chit,consv,conint,engke,engrot,taut,sigma,tolnce,tstep
      real(8) vircom,vircon,hstep,qmass,opx,opy,opz,fmx,fmy,fmz,engtrn
      real(8) ftx,fty,ftz,tqx,tqy,tqz,qt0,qt1,qt2,qt3,vaa,vbb,vcc
      real(8) taup,temp,press,virtot,vzero,chit0,chip0,cons0
      real(8) chip,volm,elrc,elrc0,virlrc,virlrc0,qstep,pmass,totmas
      real(8) volm0,scale,cxx,cyy,czz,engfke,fstep

      integer fail(nnn)
      real(8) rot(9),strkin(9),strgrp(9),cell0(9),com(3),vom(3),uni(9)

      real(8), allocatable :: dens0(:)
      real(8), allocatable :: xxt(:),yyt(:),zzt(:)
      real(8), allocatable :: dxx(:),dyy(:),dzz(:)
      real(8), allocatable :: txx(:),tyy(:),tzz(:)
      real(8), allocatable :: dxt(:),dyt(:),dzt(:)
      real(8), allocatable :: dtx(:),dty(:),dtz(:)
      real(8), allocatable :: xxo(:),yyo(:),zzo(:)
      real(8), allocatable :: vxo(:),vyo(:),vzo(:)
      real(8), allocatable :: gxo(:),gyo(:),gzo(:)
      real(8), allocatable :: oxo(:),oyo(:),ozo(:)
      real(8), allocatable :: p0(:),p1(:),p2(:),p3(:)
      real(8), allocatable :: b0(:),b1(:),b2(:),b3(:)
      real(8), allocatable :: gvxo(:),gvyo(:),gvzo(:)
      
c     metadynamics shell thermostat variables

      integer ntshl,keyshl
      real(8) sigma_shl

      logical,save :: lfirst=.true.
      real(8)      :: chit_shl  
      real(8),save :: qmass_shl
      real(8)      :: shlke

c     end metadynamics shell thermostat variables

      save newjob,hstep,qstep,fstep,pmass,qmass
      save p0,p1,p2,p3,ifre1,ifre2,igrp1,igrp2,volm0,elrc0,virlrc0
      save totmas,dens0,cell0

      data uni/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/
      data newjob/.true./
      
      safe=.true.
      
      if(newjob)then

c     allocate density storage array

        fail(1)=0
        allocate(dens0(mxatyp),stat=fail(1))
        if(fail(1).gt.0)call error(idnode,2200)

c     store intitial parameters
        
        volm0=volm
        elrc0=elrc
        virlrc0=virlrc
        hstep=0.5d0*tstep
        fstep=0.5d0*tstep/dble(ncyc)
        qstep=0.25d0*tstep/dble(ncyc)
        
        do i=1,ntpatm
          dens0(i)=dens(i)
        enddo
        
        do i=1,9
          cell0(i)=cell(i)
        enddo

c     total system mass
        
        totmas=getmass(natms,idnode,mxnode)

c     nose-hoover thermostat and barostat inertia parameter
        
        qmass=2.d0*sigma*taut**2
        pmass=2.d0*sigma*taup**2
        
c     free atom block indices
        
        ifre1=(idnode*ntfree)/mxnode+1
        ifre2=((idnode+1)*ntfree)/mxnode
        
c     group block indices
        
        igrp1=(idnode*ngrp)/mxnode+1
        igrp2=((idnode+1)*ngrp)/mxnode
        
c     check work arrays are large enough
        
        safe=(igrp2-igrp1+1.le.msgrp) 
        if(mxnode.gt.1) call gstate(safe)
        if(.not.safe)then 
          igrp=igrp2-igrp1+1
          call gimax(igrp,1,idum)
          if(idnode.eq.0) write(nrite,*) ' make msgrp >=',igrp
          call error(idnode,506)
        endif

      endif

      if(ntcons.gt.0)safe=.false.
      
c     allocate working arrays
      
      do i=1,nnn
        fail(i)=0
      enddo
      allocate(dtx(mxatms),dty(mxatms),dtz(mxatms),stat=fail(1))
      
      if(ntcons.gt.0)then
        
        allocate(gxo(mxgrp),gyo(mxgrp),gzo(mxgrp),stat=fail(2))
        allocate(oxo(mxgrp),oyo(mxgrp),ozo(mxgrp),stat=fail(3))
        allocate(xxo(mxatms),yyo(mxatms),zzo(mxatms),stat=fail(4))
        allocate(vxo(mxatms),vyo(mxatms),vzo(mxatms),stat=fail(5))
        allocate(dxx(mxcons),dyy(mxcons),dzz(mxcons),stat=fail(6))
        allocate(txx(mxatms),tyy(mxatms),tzz(mxatms),stat=fail(7))
        allocate(xxt(mxatms),yyt(mxatms),zzt(mxatms),stat=fail(8))
        allocate(dxt(mxcons),dyt(mxcons),dzt(mxcons),stat=fail(9))
        allocate(b0(mxgrp),b1(mxgrp),b2(mxgrp),b3(mxgrp),stat=fail(10))
        allocate(gvxo(mxgrp),gvyo(mxgrp),gvzo(mxgrp),stat=fail(11))
        
      endif
      if(newjob)then
        allocate(p0(mxgrp),p1(mxgrp),p2(mxgrp),p3(mxgrp),stat=fail(12))
      endif
      do i=1,nnn
        if(fail(i).gt.0)call error(idnode,2210)
      enddo

      newjob=.false.

      if(lmetadyn.and.lfirst.and.(ntshl>0))then
        if(idnode.eq.0)then
          write(*,*)"Warning - Metadynamics Modification"
          write(*,*)"========================="
          write(*,*)"Coupling core-shell motion thermostat at 1 K"
        endif
        lfirst=.false.
c     use same relaxation time for global and core-shell?
        qmass_shl=2.d0*sigma_shl*taut**2
      endif
      
      if(ntcons.gt.0)then

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

      endif

c     atom displacement from rigid body centre of mass

      jr=0
      do ig=igrp1,igrp2
        
        id=lstgtp(ig)
        
        do j=1,numgsit(id)
          
          jr=jr+1
          i=lstrgd(jr)
          
          dtx(jr)=xxx(i)-gcmx(ig)
          dty(jr)=yyy(i)-gcmy(ig)
          dtz(jr)=zzz(i)-gcmz(ig)
          
        enddo
        
      enddo

c     periodic boundary condition for displacement vectors

      call images(imcon,0,1,jr,cell,dtx,dty,dtz)
      
c     first stage of velocity verlet algorithm

      if(isw.eq.1)then

c     remove system centre of mass velocity
        
        call getvom(natms,idnode,mxnode,totmas,vom)
        
        do i=1,natms
          
          vxx(i)=vxx(i)-vom(1)
          vyy(i)=vyy(i)-vom(2)
          vzz(i)=vzz(i)-vom(3)
          
        enddo
        
        do ig=1,ngrp
          
          gvxx(ig)=gvxx(ig)-vom(1)
          gvyy(ig)=gvyy(ig)-vom(2)
          gvzz(ig)=gvzz(ig)-vom(3)
          
        enddo

c     store current integration variables if ntcons > 0
        
        if(ntcons.gt.0)then
          
          vzero=volm
          chit0=chit
          chip0=chip
          cons0=conint
          do i=1,natms
            
            xxo(i)=xxx(i)
            yyo(i)=yyy(i)
            zzo(i)=zzz(i)
            vxo(i)=vxx(i)
            vyo(i)=vyy(i)
            vzo(i)=vzz(i)
            
          enddo
          do ig=1,ngrp
            
            b0(ig)=q0(ig)
            b1(ig)=q1(ig)
            b2(ig)=q2(ig)
            b3(ig)=q3(ig)
            oxo(ig)=omx(ig)
            oyo(ig)=omy(ig)
            ozo(ig)=omz(ig)
            gxo(ig)=gcmx(ig)
            gyo(ig)=gcmy(ig)
            gzo(ig)=gcmz(ig)
            gvxo(ig)=gvxx(ig)
            gvyo(ig)=gvyy(ig)
            gvzo(ig)=gvzz(ig)
            
          enddo
          
        endif
        
      endif

c     iteration necessary if ntcons > 0 and isw=1

      mxiter=1
      if(isw.eq.1.and.ntcons.gt.0)mxiter=3
      do iter=1,mxiter
        
        if(isw.eq.1)then
          
          do icyc=1,ncyc
            
c     integrate and apply npt thermostat
            
            call nptqscl_t
     x        (idnode,mxnode,ntfree,ngrp,engfke,engtrn,engrot,temp,
     x        sigma,qstep,pmass,qmass,taut,chip,chit,conint)
            
c     metadynamics shell thermostat
c     ====================================================
c     Must first merge update velocities as the core-shell
c     velocites are not distributed according to the same
c     rules.
c     ====================================================
            if(lmetadyn.and.keyshl.eq.1)then
              if(mxnode.gt.1)call merge
     x          (idnode,mxnode,natms,mxbuff,vxx,vyy,vzz,buffer)
              call nvtscale_shl
     x          (idnode,mxnode,ntshl,shlke,sigma_shl,hstep,qmass_shl,
     x          taut,chit_shl,conint)      
            endif           
            
c     integrate and apply npt barostat
            
            call nptqscl_p
     x        (idnode,mxnode,ntfree,ngrp,engfke,engtrn,fstep,pmass,
     x        chip,chit,volm,press,vircon,virtot,vircom)
            
c     integrate and apply npt thermostat
            
            call nptqscl_t
     x        (idnode,mxnode,ntfree,ngrp,engfke,engtrn,engrot,temp,
     x        sigma,qstep,pmass,qmass,taut,chip,chit,conint)
            
c     metadynamics shell thermostat
c     ====================================================
c     Must first merge update velocities as the core-shell
c     velocites are not distributed according to the same
c     rules.
c     ====================================================
            if(lmetadyn.and.keyshl.eq.1)then
              if(mxnode.gt.1)call merge
     x          (idnode,mxnode,natms,mxbuff,vxx,vyy,vzz,buffer)
              call nvtscale_shl
     x          (idnode,mxnode,ntshl,shlke,sigma_shl,hstep,qmass_shl,
     x          taut,chit_shl,conint)      
            endif           
            
          enddo
          
          engke=engfke+engtrn
          
c     scale cell vectors - isotropic
          
          scale=(volm/volm0)**(1.d0/3.d0)
          do i=1,9
            cell(i)=cell0(i)*scale
          enddo
          
c     calculate quaternion momenta
          
          do ig=igrp1,igrp2
            
            id=lstgtp(ig)
            opx=omx(ig)*rotinx(id,1)
            opy=omy(ig)*rotiny(id,1)
            opz=omz(ig)*rotinz(id,1)
            p0(ig)=2.0d0*(-q1(ig)*opx-q2(ig)*opy-q3(ig)*opz)
            p1(ig)=2.0d0*( q0(ig)*opx-q3(ig)*opy+q2(ig)*opz)
            p2(ig)=2.0d0*( q3(ig)*opx+q0(ig)*opy-q1(ig)*opz)
            p3(ig)=2.0d0*(-q2(ig)*opx+q1(ig)*opy+q0(ig)*opz)
            
          enddo
          
        endif

c     update free atom velocities 

        do ifre=ifre1,ifre2

          i=lstfre(ifre)
          vxx(i)=vxx(i)+hstep*rmass(i)*fxx(i)
          vyy(i)=vyy(i)+hstep*rmass(i)*fyy(i)
          vzz(i)=vzz(i)+hstep*rmass(i)*fzz(i)
          
        enddo

c     *************  Rigid body motion ****************************

        jr=0
        do ig=igrp1,igrp2
          
c     calculate com force arrays 
          
          jrs=jr
          fmx=0.d0
          fmy=0.d0
          fmz=0.d0
          id=lstgtp(ig)
          do j=1,numgsit(id)
            
            jr=jr+1
            i=lstrgd(jr)
            fmx=fmx+fxx(i)
            fmy=fmy+fyy(i)
            fmz=fmz+fzz(i)
            
          enddo
          
c     current rotation matrix 
          
          call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)
          
c     calculate torque in principal frame
          
          jr=jrs
          ftx=0.d0
          fty=0.d0
          ftz=0.d0
          do j=1,numgsit(id)
            
            jr=jr+1
            i=lstrgd(jr)
            ftx=ftx+dty(jr)*fzz(i)-dtz(jr)*fyy(i)
            fty=fty+dtz(jr)*fxx(i)-dtx(jr)*fzz(i)
            ftz=ftz+dtx(jr)*fyy(i)-dty(jr)*fxx(i)
            
          enddo
          tqx=ftx*rot(1)+fty*rot(4)+ftz*rot(7)
          tqy=ftx*rot(2)+fty*rot(5)+ftz*rot(8)
          tqz=ftx*rot(3)+fty*rot(6)+ftz*rot(9)
          
c     calculate quaternion torques
          
          qt0=2.0d0*(-q1(ig)*tqx-q2(ig)*tqy-q3(ig)*tqz)
          qt1=2.0d0*( q0(ig)*tqx-q3(ig)*tqy+q2(ig)*tqz)
          qt2=2.0d0*( q3(ig)*tqx+q0(ig)*tqy-q1(ig)*tqz)
          qt3=2.0d0*(-q2(ig)*tqx+q1(ig)*tqy+q0(ig)*tqz)
          
c     update quaternion momenta by 1/2 time step
          
          p0(ig)=p0(ig)+qt0*hstep
          p1(ig)=p1(ig)+qt1*hstep
          p2(ig)=p2(ig)+qt2*hstep
          p3(ig)=p3(ig)+qt3*hstep
          
c     update centre of mass velocity by 1/2 time step
          
          gvxx(ig)=gvxx(ig)+fmx*hstep/gmass(id)
          gvyy(ig)=gvyy(ig)+fmy*hstep/gmass(id)
          gvzz(ig)=gvzz(ig)+fmz*hstep/gmass(id)
          
        enddo
        
        if(isw.eq.1)then
          
c     calculate system centre of mass
          
          call getcom(natms,idnode,mxnode,totmas,com)
          
c     move centre of mass by full time step
          
          do ig=igrp1,igrp2
            
            cxx=gcmx(ig)-com(1)
            cyy=gcmy(ig)-com(2)
            czz=gcmz(ig)-com(3)
            gcmx(ig)=gcmx(ig)+tstep*(gvxx(ig)+chip*cxx)
            gcmy(ig)=gcmy(ig)+tstep*(gvyy(ig)+chip*cyy)
            gcmz(ig)=gcmz(ig)+tstep*(gvzz(ig)+chip*czz)
            
          enddo
          
c     merge group coms from all nodes
          
          if(mxnode.gt.1)call merge
     x      (idnode,mxnode,ngrp,mxbuff,gcmx,gcmy,gcmz,buffer)
          
c     rotate rigid groups: nosquish algorithm
          
          do ig=igrp1,igrp2
            call nosquish(ig,tstep,q0,q1,q2,q3,p0,p1,p2,p3)
          enddo
          
c     new atomic positions for atoms in rigid bodies-relative to c.o.m
          
          k=0
          do ig=igrp1,igrp2
            
            call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)
            
            id=lstgtp(ig)
            
            do j=1,numgsit(id)
              
              k=k+1
              i=lstme(k)
              xxx(i)=rot(1)*gxx(id,j)+rot(2)*gyy(id,j)+rot(3)*gzz(id,j)
     x          +gcmx(ig)
              yyy(i)=rot(4)*gxx(id,j)+rot(5)*gyy(id,j)+rot(6)*gzz(id,j)
     x          +gcmy(ig)
              zzz(i)=rot(7)*gxx(id,j)+rot(8)*gyy(id,j)+rot(9)*gzz(id,j)
     x          +gcmz(ig)
              
            enddo
            
          enddo
          
c     update positions of free particles
          
          do ifre=ifre1,ifre2
            
            i=lstfre(ifre)
            cxx=xxx(i)-com(1)
            cyy=yyy(i)-com(2)
            czz=zzz(i)-com(3)
            xxx(i)=xxx(i)+tstep*(vxx(i)+chip*cxx)
            yyy(i)=yyy(i)+tstep*(vyy(i)+chip*cyy)
            zzz(i)=zzz(i)+tstep*(vzz(i)+chip*czz)
            
          enddo
          
c     merge position data
          
          if(mxnode.gt.1)call merge1
     x      (idnode,mxnode,natms,lstme,xxx,yyy,zzz,buffer)
          
c     apply shake corrections to bond constraints
          
          if(ntcons.gt.0)then
            
            if(mxnode.gt.1)call merge1
     x        (idnode,mxnode,natms,lstme,vxx,vyy,vzz,buffer)
            
            safe=.false.
            call rdrattle_r
     x        (safe,lshmov,idnode,imcon,mxnode,natms,nscons,
     x        tolnce,tstep,vircon,dxx,dyy,dzz,dxt,dyt,dzt,
     x        txx,tyy,tzz,xxt,yyt,zzt,strcns)
            if(.not.safe)return
            
          endif
          
c     restore original integration parameters if iter < mxiter
          
          if(iter.lt.mxiter)then
            
            volm=vzero
            chit=chit0
            chip=chip0
            conint=cons0
            do i=1,natms
              
              xxx(i)=xxo(i)
              yyy(i)=yyo(i)
              zzz(i)=zzo(i)
              vxx(i)=vxo(i)
              vyy(i)=vyo(i)
              vzz(i)=vzo(i)
              
            enddo
            do ig=1,ngrp
              
              q0(ig)=b0(ig)
              q1(ig)=b1(ig)
              q2(ig)=b2(ig)
              q3(ig)=b3(ig)
              omx(ig)=oxo(ig)
              omy(ig)=oyo(ig)
              omz(ig)=ozo(ig)
              gcmx(ig)=gxo(ig)
              gcmy(ig)=gyo(ig)
              gcmz(ig)=gzo(ig)
              gvxx(ig)=gvxo(ig)
              gvyy(ig)=gvyo(ig)
              gvzz(ig)=gvzo(ig)
              
            enddo
            
          endif
          
        endif

c     operations for second stage only

        if(isw.eq.2)then
          
          jr=0
          do ig=igrp1,igrp2
            
c     new angular momenta and velocities
            
            opx=pt5*(-q1(ig)*p0(ig)+q0(ig)*p1(ig)+
     x        q3(ig)*p2(ig)-q2(ig)*p3(ig))
            opy=pt5*(-q2(ig)*p0(ig)-q3(ig)*p1(ig)+
     x        q0(ig)*p2(ig)+q1(ig)*p3(ig))
            opz=pt5*(-q3(ig)*p0(ig)+q2(ig)*p1(ig)-
     x        q1(ig)*p2(ig)+q0(ig)*p3(ig))
            
            id=lstgtp(ig)
            
            omx(ig)=opx*rotinx(id,2)
            omy(ig)=opy*rotiny(id,2)
            omz(ig)=opz*rotinz(id,2)
            
c     new rotation matrix
            
            call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)
            
            do j=1,numgsit(id)
              
              jr=jr+1
              i=lstrgd(jr)
              
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
          
c     merge velocities from all nodes
          
          if(mxnode.gt.1)call merge1
     x      (idnode,mxnode,natms,lstme,vxx,vyy,vzz,buffer)
          
c     correct constraint bond velocities using rattle
          
          if(ntcons.gt.0)then
            
            safe=.false.
            call rdrattle_v
     x        (safe,idnode,mxnode,natms,nscons,tolnce,tstep,
     x        dxx,dyy,dzz,txx,tyy,tzz,xxt,yyt,zzt)
            if(.not.safe)return
            
          endif
          
c     kinetic contribution to stress tensor
          
          call kinstressf(ntfree,idnode,mxnode,strkin)        
          call kinstressg(ngrp,idnode,mxnode,strgrp)
          
c     calculate rigid body contribution to stress tensor

          call bodystress(idnode,mxnode,ngrp,vircom,strbod,dtx,dty,dtz)

c     add contributions to stress tensor
          
          do i=1,9
            stress(i)=stress(i)+strbod(i)+strcns(i)+strkin(i)+strgrp(i)
          enddo
          
          do icyc=1,ncyc
            
c     integrate and apply npt thermostat
            
            call nptqscl_t
     x        (idnode,mxnode,ntfree,ngrp,engfke,engtrn,engrot,temp,
     x        sigma,qstep,pmass,qmass,taut,chip,chit,conint)
            
c     metadynamics shell thermostat
c     ====================================================
c     Must first merge update velocities as the core-shell
c     velocites are not distributed according to the same
c     rules.
c     ====================================================
            if(lmetadyn.and.keyshl.eq.1)then
              if(mxnode.gt.1)call merge
     x          (idnode,mxnode,natms,mxbuff,vxx,vyy,vzz,buffer)
              call nvtscale_shl
     x          (idnode,mxnode,ntshl,shlke,sigma_shl,hstep,qmass_shl,
     x          taut,chit_shl,conint)      
            endif           
            
c     integrate and apply npt barostat
            
            call nptqscl_p
     x        (idnode,mxnode,ntfree,ngrp,engfke,engtrn,fstep,pmass,
     x        chip,chit,volm,press,vircon,virtot,vircom)
            
c     integrate and apply npt thermostat
            
            call nptqscl_t
     x        (idnode,mxnode,ntfree,ngrp,engfke,engtrn,engrot,temp,
     x        sigma,qstep,pmass,qmass,taut,chip,chit,conint)
            
c     metadynamics shell thermostat
c     ====================================================
c     Must first merge update velocities as the core-shell
c     velocites are not distributed according to the same
c     rules.
c     ====================================================
            if(lmetadyn.and.keyshl.eq.1)then
              if(mxnode.gt.1)call merge
     x          (idnode,mxnode,natms,mxbuff,vxx,vyy,vzz,buffer)
              call nvtscale_shl
     x          (idnode,mxnode,ntshl,shlke,sigma_shl,hstep,qmass_shl,
     x          taut,chit_shl,conint)      
            endif           
            
          enddo

          engke=engfke+engtrn
          
          if(mxnode.gt.1)call merge1
     x      (idnode,mxnode,natms,lstme,vxx,vyy,vzz,buffer)

c     scale cell vectors - isotropic
          
          scale=(volm/volm0)**(1.d0/3.d0)
          do i=1,9
            cell(i)=cell0(i)*scale
          enddo
          
c     calculate conserved variable

          consv=conint+0.5d0*qmass*chit**2+press*volm
     x      +0.5d0*pmass*chip**2
          
c     metadynamics shell thermostat
          
          if(lmetadyn.and.keyshl.eq.1)then
            consv=consv+0.5d0*qmass_shl*chit_shl**2
          endif
          
c     end of second stage of velocity verlet algorithm
          
        endif
        
c     end of iteration cycle

      enddo
          
c     periodic boundary condition
          
      if(isw.eq.2)then
        
        call images(imcon,idnode,mxnode,natms,cell,xxx,yyy,zzz)
        call images(imcon,idnode,mxnode,ngrp,cell,gcmx,gcmy,gcmz)
        
c     merge position data
        
        if(mxnode.gt.1)then
          
          call merge(idnode,mxnode,natms,mxbuff,xxx,yyy,zzz,buffer)
          call merge(idnode,mxnode,ngrp,mxbuff,gcmx,gcmy,gcmz,buffer)
          
        endif
        
      endif
      
c     adjust long range corrections and number density
      
      elrc=elrc0*(volm0/volm)
      virlrc=virlrc0*(volm0/volm)
      do kk=1,ntpatm
        dens(kk)=dens0(kk)*(volm0/volm)
      enddo
      
c     construct scaling tensor (for tethered atoms)
      
      do i=1,9
        eta(i)=chip*uni(i)
      enddo
      
c     deallocate working arrays
      
      deallocate(dtx,dty,dtz,stat=fail(1))
      
      if(ntcons.gt.0)then

        deallocate(dxx,dyy,dzz,txx,tyy,tzz,stat=fail(3))
        deallocate(xxt,yyt,zzt,dxt,dyt,dzt,stat=fail(4))
        deallocate(xxo,yyo,zzo,oxo,oyo,ozo,stat=fail(5))
        deallocate(vxo,vyo,vzo,b0,b1,b2,b3,stat=fail(6))
        deallocate(gxo,gyo,gzo,gvxo,gvyo,gvzo,stat=fail(7))

      endif
      
      return
      end subroutine nptqvv_h1

      subroutine nstqvv_b1
     x  (safe,lshmov,isw,imcon,idnode,mxnode,natms,ngrp,nscons,
     x  ntcons,ntfree,ntpatm,mode,engke,engrot,press,taut,taup,sigma,
     x  tolnce,tstep,vircom,vircon,elrc,virlrc,volm)
      
c***********************************************************************
c     
c     dlpoly routine to integrate rigid body equations of motion
c     using the symplectic no_squish quaternion algorithm of 
c     miller et al j.chem.phys 116 (2002) 8649
c     plus bond constraints - provided rigid body sites
c     and constraint sites do not coincide.
c     
c     nst ensemble - Berendsen thermostat and barostat 
c     (n.b. not symplectic)
c     
c     isothermal compressibility (beta) set to that of liquid water
c     = 0.007372 dlpoly units
c     
c     parallel replicated data version : block data
c     
c     omx,omy,omz=angular velocity in body fixed frame (principal axes)
c     rotinx,y,z =rotational inertia in body fixed frame
c     
c     copyright daresbury laboratory
c     author      w.smith may 2005
c     
c**********************************************************************

      implicit none

      logical safe,lshmov,newjob

      integer isw,imcon,idnode,mxnode,natms,ngrp,nscons,ntcons
      integer ntfree,i,j,k,igrp1,igrp2,igrp,ifre1,ifre2,jr,kk,mode
      integer id,ifre,jrs,idum,ig,iatm0,iatm1,ntpatm,iter,mxiter

      real(8) engke,engrot,tolnce,tstep,vircom,vircon,engfke
      real(8) tqx,tqy,tqz,ftx,fty,ftz
      real(8) vaa,vbb,vcc,engtrn,fmx,fmy,fmz
      real(8) qt0,qt1,qt2,qt3,opx,opy,opz,taut,sigma,engtke
      real(8) chit,beta,volm,volm0,elrc,elrc0,virlrc,virlrc0
      real(8) press,taup,xtmp,ytmp,ztmp

      integer, parameter :: nnn=11
      real(8), parameter :: pt5=0.5d0

      integer fail(nnn)
      real(8) rot(9),strkin(9),strgrp(9)
      real(8) celp(10),uni(9)

      real(8), allocatable :: dens0(:)
      real(8), allocatable :: xxo(:),yyo(:),zzo(:)
      real(8), allocatable :: vxo(:),vyo(:),vzo(:)
      real(8), allocatable :: xxt(:),yyt(:),zzt(:)
      real(8), allocatable :: gxo(:),gyo(:),gzo(:)
      real(8), allocatable :: dxx(:),dyy(:),dzz(:)
      real(8), allocatable :: txx(:),tyy(:),tzz(:)
      real(8), allocatable :: dxt(:),dyt(:),dzt(:)
      real(8), allocatable :: dtx(:),dty(:),dtz(:)
      real(8), allocatable :: p0(:),p1(:),p2(:),p3(:)
      real(8), allocatable :: b0(:),b1(:),b2(:),b3(:)
      real(8), allocatable :: c0(:),c1(:),c2(:),c3(:)
      
      save newjob,volm0,elrc0,virlrc0,iatm0,iatm1,dens0
      save p0,p1,p2,p3,ifre1,ifre2,igrp1,igrp2

      data newjob/.true./,beta/7.3728d-3/
      data uni/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/
      
      safe=.true.

      if(newjob)then
        
        fail(1)=0
        allocate(dens0(mxatyp),stat=fail(1))
        if(fail(1).gt.0)call error(idnode,2240)

c     store initial values of volume and long range corrections
        
        volm0=volm
        elrc0=elrc
        virlrc0=virlrc
        do i=1,ntpatm
          dens0(i)=dens(i)
        enddo
        
c     atom block indices
        
        iatm0=(idnode*natms)/mxnode+1
        iatm1=((idnode+1)*natms)/mxnode
        
c     free atom block indices
        
        ifre1=(idnode*ntfree)/mxnode+1
        ifre2=((idnode+1)*ntfree)/mxnode
        
c     group block indices
        
        igrp1=(idnode*ngrp)/mxnode+1
        igrp2=((idnode+1)*ngrp)/mxnode
        
c     check work arrays are large enough
        
        safe=(igrp2-igrp1+1.le.msgrp) 
        if(mxnode.gt.1) call gstate(safe)
        if(.not.safe)then 
          igrp=igrp2-igrp1+1
          call gimax(igrp,1,idum)
          if(idnode.eq.0) write(nrite,*) ' make msgrp >=',igrp
          call  error(idnode,506)
        endif

      endif

c     allocate working arrays

      do i=1,nnn
        fail(i)=0
      enddo
      allocate(dtx(mxatms),dty(mxatms),dtz(mxatms),stat=fail(1))
      if(ntcons.gt.0)then

        allocate(dxx(mxcons),dyy(mxcons),dzz(mxcons),stat=fail(2))
        allocate(txx(mxatms),tyy(mxatms),tzz(mxatms),stat=fail(3))
        allocate(xxt(mxatms),yyt(mxatms),zzt(mxatms),stat=fail(4))
        allocate(dxt(mxcons),dyt(mxcons),dzt(mxcons),stat=fail(5))

        if(isw.eq.1)then

          allocate(gxo(mxgrp),gyo(mxgrp),gzo(mxgrp),stat=fail(6))
          allocate(b0(mxgrp),b1(mxgrp),b2(mxgrp),b3(mxgrp),stat=fail(7))
          allocate(c0(mxgrp),c1(mxgrp),c2(mxgrp),c3(mxgrp),stat=fail(8))
          allocate(xxo(mxatms),yyo(mxatms),zzo(mxatms),stat=fail(9))
          allocate(vxo(mxatms),vyo(mxatms),vzo(mxatms),stat=fail(10))

        endif

      endif
      if(newjob)then         
        allocate(p0(mxgrp),p1(mxgrp),p2(mxgrp),p3(mxgrp),stat=fail(11))
      endif
      do i=1,nnn
        if(fail(i).gt.0)call error(idnode,2250)
      enddo

      newjob=.false.
      if(ntcons.gt.0)safe=.false.

c     atom displacement from rigid body centre of mass

      jr=0
      do ig=igrp1,igrp2
        
        id=lstgtp(ig)
        
        do j=1,numgsit(id)
          
          jr=jr+1
          i=lstrgd(jr)
          
          dtx(jr)=xxx(i)-gcmx(ig)
          dty(jr)=yyy(i)-gcmy(ig)
          dtz(jr)=zzz(i)-gcmz(ig)
          
        enddo
        
      enddo

c     periodic boundary condition for displacement vectors

      call images(imcon,0,1,jr,cell,dtx,dty,dtz)
      
c     construct current bond vectors
      
      if(ntcons.gt.0)then

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

      endif

      if(isw.eq.1)then

c     calculate kinetic energy

        engfke=getkinf(ntfree,idnode,mxnode)
        call getking(ngrp,idnode,mxnode,engtrn,engrot)
        engke=engfke+engtrn

c     calculate quaternion momenta at start of time step

        do ig=igrp1,igrp2
          
          id=lstgtp(ig)
          opx=omx(ig)*rotinx(id,1)
          opy=omy(ig)*rotiny(id,1)
          opz=omz(ig)*rotinz(id,1)
          p0(ig)=2.0d0*(-q1(ig)*opx-q2(ig)*opy-q3(ig)*opz)
          p1(ig)=2.0d0*( q0(ig)*opx-q3(ig)*opy+q2(ig)*opz)
          p2(ig)=2.0d0*( q3(ig)*opx+q0(ig)*opy-q1(ig)*opz)
          p3(ig)=2.0d0*(-q2(ig)*opx+q1(ig)*opy+q0(ig)*opz)
          
        enddo

      endif

c     update free atom velocities 1/2 time step first and second stages

      do ifre=ifre1,ifre2

        i=lstfre(ifre)
        vxx(i)=vxx(i)+pt5*tstep*rmass(i)*fxx(i)
        vyy(i)=vyy(i)+pt5*tstep*rmass(i)*fyy(i)
        vzz(i)=vzz(i)+pt5*tstep*rmass(i)*fzz(i)
        
      enddo

c     rigid body motion for first and second stages

      jr=0
      do ig=igrp1,igrp2
        
c     calculate com force arrays 
        
        jrs=jr
        fmx=0.d0
        fmy=0.d0
        fmz=0.d0
        id=lstgtp(ig)
        do j=1,numgsit(id)
          
          jr=jr+1
          i=lstrgd(jr)
          fmx=fmx+fxx(i)
          fmy=fmy+fyy(i)
          fmz=fmz+fzz(i)

        enddo

c     current rotation matrix 
        
        call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)

c     calculate torque in principal frame

        jr=jrs
        ftx=0.d0
        fty=0.d0
        ftz=0.d0
        do j=1,numgsit(id)
          jr=jr+1
          i=lstrgd(jr)
          ftx=ftx+dty(jr)*fzz(i)-dtz(jr)*fyy(i)
          fty=fty+dtz(jr)*fxx(i)-dtx(jr)*fzz(i)
          ftz=ftz+dtx(jr)*fyy(i)-dty(jr)*fxx(i)
        enddo
        tqx=ftx*rot(1)+fty*rot(4)+ftz*rot(7)
        tqy=ftx*rot(2)+fty*rot(5)+ftz*rot(8)
        tqz=ftx*rot(3)+fty*rot(6)+ftz*rot(9)

c     calculate quaternion torques

        qt0=2.0d0*(-q1(ig)*tqx-q2(ig)*tqy-q3(ig)*tqz)
        qt1=2.0d0*( q0(ig)*tqx-q3(ig)*tqy+q2(ig)*tqz)
        qt2=2.0d0*( q3(ig)*tqx+q0(ig)*tqy-q1(ig)*tqz)
        qt3=2.0d0*(-q2(ig)*tqx+q1(ig)*tqy+q0(ig)*tqz)

c     update quaternion momenta by 1/2 time step

        p0(ig)=p0(ig)+qt0*pt5*tstep
        p1(ig)=p1(ig)+qt1*pt5*tstep
        p2(ig)=p2(ig)+qt2*pt5*tstep
        p3(ig)=p3(ig)+qt3*pt5*tstep

c     update centre of mass velocity by 1/2 time step

        gvxx(ig)=gvxx(ig)+fmx*pt5*tstep/gmass(id)
        gvyy(ig)=gvyy(ig)+fmy*pt5*tstep/gmass(id)
        gvzz(ig)=gvzz(ig)+fmz*pt5*tstep/gmass(id)

      enddo

c     first stage of velocity verlet algorithm

      if(isw.eq.1)then

c     store current integration variables
        
        if(ntcons.gt.0)then

          do i=1,natms
            
            xxo(i)=xxx(i)
            yyo(i)=yyy(i)
            zzo(i)=zzz(i)
            vxo(i)=vxx(i)
            vyo(i)=vyy(i)
            vzo(i)=vzz(i)

          enddo
          do ig=1,ngrp
            
            b0(ig)=q0(ig)
            b1(ig)=q1(ig)
            b2(ig)=q2(ig)
            b3(ig)=q3(ig)
            c0(ig)=p0(ig)
            c1(ig)=p1(ig)
            c2(ig)=p2(ig)
            c3(ig)=p3(ig)
            gxo(ig)=gcmx(ig)
            gyo(ig)=gcmy(ig)
            gzo(ig)=gcmz(ig)
            
          enddo

        endif

c     extract previous constraint terms from stress tensor

      if(isw.eq.1)then          

        do i=1,9
          stress(i)=stress(i)-strcns(i)
        enddo

      endif

c     iteration required if ntcons > 0

        mxiter=1
        if(ntcons.gt.0)mxiter=2
        do iter=1,mxiter

c     zero scaling matrix

          do i=1,9
            eta(i)=uni(i)
          enddo

          if(iter.eq.mxiter)then
            
c     calculate Berendsen barostat
            
            do i=1,9
              eta(i)=tstep*beta*(stress(i)+strcns(i)-
     x          press*volm*uni(i))/(taup*volm)+uni(i)
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

c     reset cell parameters for new volume
          
            call mat_mul(eta,cell,cell)
            
c     calculate new volume
            
            call dcell(cell,celp)
            volm=celp(10)

          endif

c     update centre of mass position

          do ig=igrp1,igrp2

            xtmp=eta(1)*gcmx(ig)+eta(4)*gcmy(ig)+eta(7)*gcmz(ig)
            ytmp=eta(2)*gcmx(ig)+eta(5)*gcmy(ig)+eta(8)*gcmz(ig)
            ztmp=eta(3)*gcmx(ig)+eta(6)*gcmy(ig)+eta(9)*gcmz(ig)
            gcmx(ig)=tstep*gvxx(ig)+xtmp
            gcmy(ig)=tstep*gvyy(ig)+ytmp
            gcmz(ig)=tstep*gvzz(ig)+ztmp

          enddo

c     merge group coms from all nodes

          if(mxnode.gt.1)call merge
     x      (idnode,mxnode,ngrp,mxbuff,gcmx,gcmy,gcmz,buffer)

c     rotate rigid groups: nosquish algorithm

          jr=0
          do ig=igrp1,igrp2
            call nosquish(ig,tstep,q0,q1,q2,q3,p0,p1,p2,p3)
          enddo

c     new atomic positions for atoms in rigid bodies-relative to c.o.m
          
          k=0
          do ig=igrp1,igrp2
            
            call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)
            
            id=lstgtp(ig)
            
            do j=1,numgsit(id)
              
              k=k+1
              i=lstme(k)
              xxx(i)=gcmx(ig)+
     x          rot(1)*gxx(id,j)+rot(2)*gyy(id,j)+rot(3)*gzz(id,j)
              yyy(i)=gcmy(ig)+
     x          rot(4)*gxx(id,j)+rot(5)*gyy(id,j)+rot(6)*gzz(id,j)
              zzz(i)=gcmz(ig)+
     x          rot(7)*gxx(id,j)+rot(8)*gyy(id,j)+rot(9)*gzz(id,j)
              
            enddo
            
          enddo
          
c     update positions of free particles
          
          do ifre=ifre1,ifre2
            
            i=lstfre(ifre)
            xxx(i)=tstep*vxx(i)+
     x        eta(1)*xxx(i)+eta(4)*yyy(i)+eta(7)*zzz(i)
            yyy(i)=tstep*vyy(i)+
     x        eta(2)*xxx(i)+eta(5)*yyy(i)+eta(8)*zzz(i)
            zzz(i)=tstep*vzz(i)+
     x        eta(3)*xxx(i)+eta(6)*yyy(i)+eta(9)*zzz(i)
            
          enddo
          
c     merge position data

          if(mxnode.gt.1)call merge1
     x      (idnode,mxnode,natms,lstme,xxx,yyy,zzz,buffer)
          
c     apply shake corrections to bond constraints
          
          if(ntcons.gt.0)then
            
            if(mxnode.gt.1)call merge1
     x        (idnode,mxnode,natms,lstme,vxx,vyy,vzz,buffer)

            safe=.false.
            call rdrattle_r
     x        (safe,lshmov,idnode,imcon,mxnode,natms,nscons,
     x        tolnce,tstep,vircon,dxx,dyy,dzz,dxt,dyt,dzt,
     x        txx,tyy,tzz,xxt,yyt,zzt,strcns)
            if(.not.safe)return

          endif

          if(iter.lt.mxiter)then

            do i=1,natms

              xxx(i)=xxo(i)
              yyy(i)=yyo(i)
              zzz(i)=zzo(i)
              vxx(i)=vxo(i)
              vyy(i)=vyo(i)
              vzz(i)=vzo(i)

            enddo
            do ig=1,ngrp

              q0(ig)=b0(ig)
              q1(ig)=b1(ig)
              q2(ig)=b2(ig)
              q3(ig)=b3(ig)
              p0(ig)=c0(ig)
              p1(ig)=c1(ig)
              p2(ig)=c2(ig)
              p3(ig)=c3(ig)
              gcmx(ig)=gxo(ig)
              gcmy(ig)=gyo(ig)
              gcmz(ig)=gzo(ig)

            enddo

          endif

        enddo

c     adjust long range corrections and number density
        
        elrc=elrc0*(volm0/volm)
        virlrc=virlrc0*(volm0/volm)
        
        do kk=1,ntpatm
          dens(kk)=dens0(kk)*(volm0/volm)
        enddo

c     end of first stage of velocity verlet algorithm

      else

c     second stage of velocity verlet algorithm

        jr=0
        do ig=igrp1,igrp2

c     new angular momenta and velocities
          
          opx=pt5*(-q1(ig)*p0(ig)+q0(ig)*p1(ig)+
     x      q3(ig)*p2(ig)-q2(ig)*p3(ig))
          opy=pt5*(-q2(ig)*p0(ig)-q3(ig)*p1(ig)+
     x      q0(ig)*p2(ig)+q1(ig)*p3(ig))
          opz=pt5*(-q3(ig)*p0(ig)+q2(ig)*p1(ig)-
     x      q1(ig)*p2(ig)+q0(ig)*p3(ig))

          id=lstgtp(ig)

          omx(ig)=opx*rotinx(id,2)
          omy(ig)=opy*rotiny(id,2)
          omz(ig)=opz*rotinz(id,2)
          
c     new rotation matrix
          
          call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)

          do j=1,numgsit(id)
            
            jr=jr+1
            i=lstrgd(jr)
            
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

c     merge velocities and forces from all nodes

        if(mxnode.gt.1)call merge1
     x    (idnode,mxnode,natms,lstme,vxx,vyy,vzz,buffer)

c     correct constraint bond velocities using rattle

        if(ntcons.gt.0)then

          safe=.false.
          call rdrattle_v
     x      (safe,idnode,mxnode,natms,nscons,tolnce,tstep,
     x      dxx,dyy,dzz,txx,tyy,tzz,xxt,yyt,zzt)
          if(.not.safe)return

        endif

c     calculate rigid body contribution to stress tensor

        call bodystress(idnode,mxnode,ngrp,vircom,strbod,dtx,dty,dtz)

c     calculate kinetic energy
        
        if(isw.eq.2)then

          engfke=getkinf(ntfree,idnode,mxnode)
          call getking(ngrp,idnode,mxnode,engtrn,engrot)
          engtke=engfke+engtrn+engrot
          engke=engfke+engtrn

c     apply Berendsen thermostat - taut is the relaxation time
          
          chit=sqrt(1.d0+tstep/taut*(sigma/engtke-1.d0))
          
          engke=engke*chit**2
          engtrn=engtrn*chit**2
          engrot=engrot*chit**2
          
c     thermostat velocities 
          
          do i=iatm0,iatm1
            
            if(lstfrz(i).eq.0)then

              vxx(i)=chit*vxx(i)
              vyy(i)=chit*vyy(i)
              vzz(i)=chit*vzz(i)
              
            endif
            
          enddo

c     merge velocities from all nodes

          if(mxnode.gt.1)call merge
     x      (idnode,mxnode,natms,mxbuff,vxx,vyy,vzz,buffer)

c     thermostat rigid body velocities

          do ig=igrp1,igrp2
            
            omx(ig)=chit*omx(ig)
            omy(ig)=chit*omy(ig)
            omz(ig)=chit*omz(ig)
            gvxx(ig)=chit*gvxx(ig)
            gvyy(ig)=chit*gvyy(ig)
            gvzz(ig)=chit*gvzz(ig)

          enddo

c     merge group velocities from all processors

          if(mxnode.gt.1)call merge
     x      (idnode,mxnode,ngrp,mxbuff,gvxx,gvyy,gvzz,buffer)

c     kinetic contribution to stress tensor
          
          call kinstressf(ntfree,idnode,mxnode,strkin)        
          call kinstressg(ngrp,idnode,mxnode,strgrp)
          
c     add contributions to stress tensor
          
          do i=1,9
            stress(i)=stress(i)+strbod(i)+strcns(i)+strkin(i)+strgrp(i)
          enddo
          
        endif

c     end of second stage of velocity verlet algorithm
        
      endif
      
c     periodic boundary condition
          
      if(isw.eq.2)then
        
        call images(imcon,idnode,mxnode,natms,cell,xxx,yyy,zzz)
        call images(imcon,idnode,mxnode,ngrp,cell,gcmx,gcmy,gcmz)
        
c     merge position data
        
        if(mxnode.gt.1)then
          
          call merge(idnode,mxnode,natms,mxbuff,xxx,yyy,zzz,buffer)
          call merge(idnode,mxnode,ngrp,mxbuff,gcmx,gcmy,gcmz,buffer)
          
        endif
        
      endif

c     deallocate working arrays

      deallocate(dtx,dty,dtz,stat=fail(1))
      if(ntcons.gt.0)then

        deallocate(dxx,dyy,dzz,txx,tyy,tzz,stat=fail(2))
        deallocate(xxt,yyt,zzt,dxt,dyt,dzt,stat=fail(3))
        
        if(isw.eq.1)then
          
          deallocate(vxo,vyo,vzo,b0,b1,b2,b3,stat=fail(4))
          deallocate(xxo,yyo,zzo,c0,c1,c2,c3,stat=fail(5))
          deallocate(gxo,gyo,gzo,stat=fail(6))
          
        endif
        
      endif
      
      return
      end subroutine nstqvv_b1

      subroutine nstqvv_h1
     x  (safe,lshmov,isw,idnode,mxnode,natms,imcon,ngrp,nscons,
     x  ntcons,ntpatm,ntfree,mode,ntshl,keyshl,tstep,taut,taup,
     x  sigma,temp,chit,consv,conint,engke,engrot,elrc,tolnce,
     x  vircon,virlrc,vircom,volm,press,chit_shl,sigma_shl)
      
c***********************************************************************
c     
c     dlpoly routine to integrate rigid body equations of motion
c     using the symplectic no_squish quaternion algorithm of 
c     miller et al j.chem.phys 116 (2002) 8649
c     plus bond constraints - provided rigid body sites
c     and constraint sites do not coincide.
c     
c     nst ensemble - nose-hoover thermostat Molec Phys 87 (1996) 1117
c     
c     parallel replicated data version : block data
c     
c     omx,omy,omz=angular velocity in body fixed frame (principal axes)
c     rotinx,y,z =rotational inertia in body fixed frame
c     
c     copyright daresbury laboratory
c     author      w.smith may 2005
c     adapted     d.quigley : metadynamics
c     
c**********************************************************************

      implicit none

      integer, parameter :: nnn=12
      integer, parameter :: ncyc=5
      real(8), parameter :: pt5=0.5d0

      logical safe,lshmov,newjob
      integer isw,imcon,idnode,mxnode,natms,ngrp,nscons,ntcons,ntfree
      integer ifre1,ifre2,igrp1,igrp2,igrp,i,j,k,jr,mode
      integer id,ig,ifre,jrs,idum,mxiter,iter,ntpatm,kk,icyc
      real(8) chit,consv,conint,engke,engrot,taut,sigma,tolnce,tstep
      real(8) vircom,vircon,hstep,qmass,opx,opy,opz,fmx,fmy,fmz,engtrn
      real(8) ftx,fty,ftz,tqx,tqy,tqz,qt0,qt1,qt2,qt3,vaa,vbb,vcc
      real(8) taup,temp,press,vzero,chit0,cons0
      real(8) chip2,volm,elrc,elrc0,virlrc,virlrc0,qstep,pmass,totmas
      real(8) volm0,cxx,cyy,czz,engfke,fstep

      integer fail(nnn)
      real(8) rot(9),strkin(9),strgrp(9),com(3),vom(3)
      real(8) czero(9),eta0(9)

      real(8), allocatable :: dens0(:)
      real(8), allocatable :: xxt(:),yyt(:),zzt(:)
      real(8), allocatable :: dxx(:),dyy(:),dzz(:)
      real(8), allocatable :: txx(:),tyy(:),tzz(:)
      real(8), allocatable :: dxt(:),dyt(:),dzt(:)
      real(8), allocatable :: dtx(:),dty(:),dtz(:)
      real(8), allocatable :: xxo(:),yyo(:),zzo(:)
      real(8), allocatable :: vxo(:),vyo(:),vzo(:)
      real(8), allocatable :: gxo(:),gyo(:),gzo(:)
      real(8), allocatable :: oxo(:),oyo(:),ozo(:)
      real(8), allocatable :: p0(:),p1(:),p2(:),p3(:)
      real(8), allocatable :: b0(:),b1(:),b2(:),b3(:)
      real(8), allocatable :: gvxo(:),gvyo(:),gvzo(:)
      
c     metadynamics shell thermostat variables

      integer ntshl,keyshl
      real(8) sigma_shl

      logical,save :: lfirst=.true.
      real(8)      :: chit_shl  
      real(8),save :: qmass_shl
      real(8)      :: shlke

c     end metadynamics shell thermostat variables
      
      save newjob,hstep,qstep,fstep,pmass,qmass
      save p0,p1,p2,p3,ifre1,ifre2,igrp1,igrp2,volm0,elrc0,virlrc0
      save totmas,dens0

      data newjob/.true./
      
      safe=.true.
      if(newjob)then

c     allocate density storage array

        fail(1)=0
        allocate(dens0(mxatyp),stat=fail(1))
        if(fail(1).gt.0)call error(idnode,2280)

c     store intitial parameters
        
        volm0=volm
        elrc0=elrc
        virlrc0=virlrc
        hstep=0.5d0*tstep
        fstep=0.5d0*tstep/dble(ncyc)
        qstep=0.25d0*tstep/dble(ncyc)
        
        do i=1,ntpatm
          dens0(i)=dens(i)
        enddo
        
c     total system mass
        
        totmas=getmass(natms,idnode,mxnode)

c     nose-hoover thermostat and barostat inertia parameter
        
        qmass=2.d0*sigma*taut**2
        pmass=2.d0*sigma*taup**2
        
c     free atom block indices
        
        ifre1=(idnode*ntfree)/mxnode+1
        ifre2=((idnode+1)*ntfree)/mxnode
        
c     group block indices
        
        igrp1=(idnode*ngrp)/mxnode+1
        igrp2=((idnode+1)*ngrp)/mxnode
        
c     check work arrays are large enough
        
        safe=(igrp2-igrp1+1.le.msgrp) 
        if(mxnode.gt.1) call gstate(safe)
        if(.not.safe)then 
          igrp=igrp2-igrp1+1
          call gimax(igrp,1,idum)
          if(idnode.eq.0) write(nrite,*) ' make msgrp >=',igrp
          call error(idnode,506)
        endif

      endif

      if(ntcons.gt.0)safe=.false.
      
c     allocate working arrays
      
      do i=1,nnn
        fail(i)=0
      enddo
      allocate(dtx(mxatms),dty(mxatms),dtz(mxatms),stat=fail(1))
      if(ntcons.gt.0)then

        allocate(gxo(mxgrp),gyo(mxgrp),gzo(mxgrp),stat=fail(2))
        allocate(oxo(mxgrp),oyo(mxgrp),ozo(mxgrp),stat=fail(3))
        allocate(xxo(mxatms),yyo(mxatms),zzo(mxatms),stat=fail(4))
        allocate(vxo(mxatms),vyo(mxatms),vzo(mxatms),stat=fail(5))
        allocate(dxx(mxcons),dyy(mxcons),dzz(mxcons),stat=fail(6))
        allocate(txx(mxatms),tyy(mxatms),tzz(mxatms),stat=fail(7))
        allocate(xxt(mxatms),yyt(mxatms),zzt(mxatms),stat=fail(8))
        allocate(dxt(mxcons),dyt(mxcons),dzt(mxcons),stat=fail(9))
        allocate(b0(mxgrp),b1(mxgrp),b2(mxgrp),b3(mxgrp),stat=fail(10))
        allocate(gvxo(mxgrp),gvyo(mxgrp),gvzo(mxgrp),stat=fail(11))

      endif
      if(newjob)then
        allocate(p0(mxgrp),p1(mxgrp),p2(mxgrp),p3(mxgrp),stat=fail(12))
      endif
      do i=1,nnn
        if(fail(i).gt.0)call error(idnode,2290)
      enddo

      newjob=.false.

      if(lmetadyn.and.lfirst.and.(ntshl>0))then
        if(idnode.eq.0)then
          write(*,*)"Warning - Metadynamics Modification"
          write(*,*)"========================="
          write(*,*)"Coupling core-shell motion thermostat at 1 K"
        endif
        lfirst=.false.
c     use same relaxation time for global and core-shell?
        qmass_shl=2.d0*sigma_shl*taut**2
      endif
      
      if(ntcons.gt.0)then

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

      endif

c     atom displacement from rigid body centre of mass

      jr=0
      do ig=igrp1,igrp2
        
        id=lstgtp(ig)
        
        do j=1,numgsit(id)
          
          jr=jr+1
          i=lstrgd(jr)
          
          dtx(jr)=xxx(i)-gcmx(ig)
          dty(jr)=yyy(i)-gcmy(ig)
          dtz(jr)=zzz(i)-gcmz(ig)
          
        enddo
        
      enddo

c     periodic boundary condition for displacement vectors

      call images(imcon,0,1,jr,cell,dtx,dty,dtz)
      
c     first stage of velocity verlet algorithm

      if(isw.eq.1)then

c     remove system centre of mass velocity
        
        call getvom(natms,idnode,mxnode,totmas,vom)
        
        do i=1,natms
          
          vxx(i)=vxx(i)-vom(1)
          vyy(i)=vyy(i)-vom(2)
          vzz(i)=vzz(i)-vom(3)
          
        enddo
        
        do ig=1,ngrp
          
          gvxx(ig)=gvxx(ig)-vom(1)
          gvyy(ig)=gvyy(ig)-vom(2)
          gvzz(ig)=gvzz(ig)-vom(3)
          
        enddo

c     store current integration variables if ntcons > 0
        
        if(ntcons.gt.0)then
          
          vzero=volm
          chit0=chit
          cons0=conint
          do i=1,9

            eta0(i)=eta(i)
            czero(i)=cell(i)

          enddo
          do i=1,natms
            
            xxo(i)=xxx(i)
            yyo(i)=yyy(i)
            zzo(i)=zzz(i)
            vxo(i)=vxx(i)
            vyo(i)=vyy(i)
            vzo(i)=vzz(i)
            
          enddo
          do ig=1,ngrp
            
            b0(ig)=q0(ig)
            b1(ig)=q1(ig)
            b2(ig)=q2(ig)
            b3(ig)=q3(ig)
            oxo(ig)=omx(ig)
            oyo(ig)=omy(ig)
            ozo(ig)=omz(ig)
            gxo(ig)=gcmx(ig)
            gyo(ig)=gcmy(ig)
            gzo(ig)=gcmz(ig)
            gvxo(ig)=gvxx(ig)
            gvyo(ig)=gvyy(ig)
            gvzo(ig)=gvzz(ig)
            
          enddo
          
        endif
        
      endif

c     iteration necessary if ntcons > 0 and isw=1

      mxiter=1
      if(isw.eq.1.and.ntcons.gt.0)mxiter=3
      do iter=1,mxiter

        if(isw.eq.1)then

          do icyc=1,ncyc
            
c     integrate and apply nst thermostat
            
            call nstqscl_t
     x        (idnode,mxnode,ntfree,ngrp,mode,engfke,engtrn,engrot,
     x        temp,sigma,qstep,pmass,qmass,taut,chit,conint)
            
c     metadynamics shell thermostat
c     ====================================================
c     Must first merge update velocities as the core-shell
c     velocites are not distributed according to the same
c     rules.
c     ====================================================
            if(lmetadyn.and.keyshl.eq.1)then
              if(mxnode.gt.1)call merge
     x          (idnode,mxnode,natms,mxbuff,vxx,vyy,vzz,buffer)
              call nvtscale_shl
     x          (idnode,mxnode,ntshl,shlke,sigma_shl,hstep,qmass_shl,
     x          taut,chit_shl,conint)      
            endif  
            
c     integrate and apply nst barostat
            
            call nstqscl_p
     x        (idnode,mxnode,ntfree,ngrp,mode,fstep,pmass,chit,
     x        press,volm)
            
c     integrate and apply nst thermostat
            
            call nstqscl_t
     x        (idnode,mxnode,ntfree,ngrp,mode,engfke,engtrn,engrot,
     x        temp,sigma,qstep,pmass,qmass,taut,chit,conint)
          
c     metadynamics shell thermostat
c     ====================================================
c     Must first merge update velocities as the core-shell
c     velocites are not distributed according to the same
c     rules.
c     ====================================================
            if(lmetadyn.and.keyshl.eq.1)then
              if(mxnode.gt.1)call merge
     x          (idnode,mxnode,natms,mxbuff,vxx,vyy,vzz,buffer)
              call nvtscale_shl
     x          (idnode,mxnode,ntshl,shlke,sigma_shl,hstep,qmass_shl,
     x          taut,chit_shl,conint)      
            endif  
            
          enddo
          
          engke=engfke+engtrn
          
c     calculate quaternion momenta
          
          do ig=igrp1,igrp2
            
            id=lstgtp(ig)
            opx=omx(ig)*rotinx(id,1)
            opy=omy(ig)*rotiny(id,1)
            opz=omz(ig)*rotinz(id,1)
            p0(ig)=2.0d0*(-q1(ig)*opx-q2(ig)*opy-q3(ig)*opz)
            p1(ig)=2.0d0*( q0(ig)*opx-q3(ig)*opy+q2(ig)*opz)
            p2(ig)=2.0d0*( q3(ig)*opx+q0(ig)*opy-q1(ig)*opz)
            p3(ig)=2.0d0*(-q2(ig)*opx+q1(ig)*opy+q0(ig)*opz)
            
          enddo
          
        endif

c     update free atom velocities 

        do ifre=ifre1,ifre2

          i=lstfre(ifre)
          vxx(i)=vxx(i)+hstep*rmass(i)*fxx(i)
          vyy(i)=vyy(i)+hstep*rmass(i)*fyy(i)
          vzz(i)=vzz(i)+hstep*rmass(i)*fzz(i)
          
        enddo

c     *************  Rigid body motion ****************************

        jr=0
        do ig=igrp1,igrp2
          
c     calculate com force arrays 
          
          jrs=jr
          fmx=0.d0
          fmy=0.d0
          fmz=0.d0
          id=lstgtp(ig)
          do j=1,numgsit(id)
            
            jr=jr+1
            i=lstrgd(jr)
            fmx=fmx+fxx(i)
            fmy=fmy+fyy(i)
            fmz=fmz+fzz(i)

          enddo

c     current rotation matrix 
          
          call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)

c     calculate torque in principal frame

          jr=jrs
          ftx=0.d0
          fty=0.d0
          ftz=0.d0
          do j=1,numgsit(id)

            jr=jr+1
            i=lstrgd(jr)
            ftx=ftx+dty(jr)*fzz(i)-dtz(jr)*fyy(i)
            fty=fty+dtz(jr)*fxx(i)-dtx(jr)*fzz(i)
            ftz=ftz+dtx(jr)*fyy(i)-dty(jr)*fxx(i)

          enddo
          tqx=ftx*rot(1)+fty*rot(4)+ftz*rot(7)
          tqy=ftx*rot(2)+fty*rot(5)+ftz*rot(8)
          tqz=ftx*rot(3)+fty*rot(6)+ftz*rot(9)

c     calculate quaternion torques

          qt0=2.0d0*(-q1(ig)*tqx-q2(ig)*tqy-q3(ig)*tqz)
          qt1=2.0d0*( q0(ig)*tqx-q3(ig)*tqy+q2(ig)*tqz)
          qt2=2.0d0*( q3(ig)*tqx+q0(ig)*tqy-q1(ig)*tqz)
          qt3=2.0d0*(-q2(ig)*tqx+q1(ig)*tqy+q0(ig)*tqz)

c     update quaternion momenta by 1/2 time step

          p0(ig)=p0(ig)+qt0*hstep
          p1(ig)=p1(ig)+qt1*hstep
          p2(ig)=p2(ig)+qt2*hstep
          p3(ig)=p3(ig)+qt3*hstep

c     update centre of mass velocity by 1/2 time step

          gvxx(ig)=gvxx(ig)+fmx*hstep/gmass(id)
          gvyy(ig)=gvyy(ig)+fmy*hstep/gmass(id)
          gvzz(ig)=gvzz(ig)+fmz*hstep/gmass(id)

        enddo

        if(isw.eq.1)then

c     calculate system centre of mass
          
          call getcom(natms,idnode,mxnode,totmas,com)

c     move centre of mass by full time step

          do ig=igrp1,igrp2

            cxx=gcmx(ig)-com(1)
            cyy=gcmy(ig)-com(2)
            czz=gcmz(ig)-com(3)
            gcmx(ig)=gcmx(ig)+
     x        tstep*(gvxx(ig)+eta(1)*cxx+eta(4)*cyy+eta(7)*czz)
            gcmy(ig)=gcmy(ig)+
     x        tstep*(gvyy(ig)+eta(2)*cxx+eta(5)*cyy+eta(8)*czz)
            gcmz(ig)=gcmz(ig)+
     x        tstep*(gvzz(ig)+eta(3)*cxx+eta(6)*cyy+eta(9)*czz)

          enddo

c     merge group coms from all nodes

          if(mxnode.gt.1)call merge
     x      (idnode,mxnode,ngrp,mxbuff,gcmx,gcmy,gcmz,buffer)

c     rotate rigid groups: nosquish algorithm

          do ig=igrp1,igrp2
            call nosquish(ig,tstep,q0,q1,q2,q3,p0,p1,p2,p3)
          enddo

c     new atomic positions for atoms in rigid bodies-relative to c.o.m
          
          k=0
          do ig=igrp1,igrp2

            call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)
            
            id=lstgtp(ig)

            do j=1,numgsit(id)
              
              k=k+1
              i=lstme(k)
              xxx(i)=rot(1)*gxx(id,j)+rot(2)*gyy(id,j)+rot(3)*gzz(id,j)
     x          +gcmx(ig)
              yyy(i)=rot(4)*gxx(id,j)+rot(5)*gyy(id,j)+rot(6)*gzz(id,j)
     x          +gcmy(ig)
              zzz(i)=rot(7)*gxx(id,j)+rot(8)*gyy(id,j)+rot(9)*gzz(id,j)
     x          +gcmz(ig)

            enddo
            
          enddo

c     update positions of free particles
          
          do ifre=ifre1,ifre2
            
            k=k+1
            i=lstfre(ifre)
            cxx=xxx(i)-com(1)
            cyy=yyy(i)-com(2)
            czz=zzz(i)-com(3)
            xxx(i)=xxx(i)+
     x        tstep*(vxx(i)+eta(1)*cxx+eta(4)*cyy+eta(7)*czz)
            yyy(i)=yyy(i)+
     x        tstep*(vyy(i)+eta(2)*cxx+eta(5)*cyy+eta(8)*czz)
            zzz(i)=zzz(i)+
     x        tstep*(vzz(i)+eta(3)*cxx+eta(6)*cyy+eta(9)*czz)
            
          enddo
          
c     merge position data

          if(mxnode.gt.1)call merge1
     x      (idnode,mxnode,natms,lstme,xxx,yyy,zzz,buffer)

c     apply shake corrections to bond constraints

          if(ntcons.gt.0)then
            
            if(mxnode.gt.1)call merge1
     x        (idnode,mxnode,natms,lstme,vxx,vyy,vzz,buffer)

c     subtract old constraint terms from stress tensor
            
            do i=1,9
              stress(i)=stress(i)-strcns(i)
            enddo
            
c     correct constraint bonds using rattle

            safe=.false.
            call rdrattle_r
     x        (safe,lshmov,idnode,imcon,mxnode,natms,nscons,
     x        tolnce,tstep,vircon,dxx,dyy,dzz,dxt,dyt,dzt,
     x        txx,tyy,tzz,xxt,yyt,zzt,strcns)
            if(.not.safe)return

c     add new constraint terms to stress tensor

            do i=1,9
              stress(i)=stress(i)+strcns(i)
            enddo

          endif

c     restore original integration parameters if iter < mxiter

          if(iter.lt.mxiter)then
            
            volm=vzero
            chit=chit0
            conint=cons0
            do i=1,9
              
              eta(i)=eta0(i)
              cell(i)=czero(i)
              
            enddo
            do i=1,natms
              
              xxx(i)=xxo(i)
              yyy(i)=yyo(i)
              zzz(i)=zzo(i)
              vxx(i)=vxo(i)
              vyy(i)=vyo(i)
              vzz(i)=vzo(i)
              
            enddo
            do ig=1,ngrp
              
              q0(ig)=b0(ig)
              q1(ig)=b1(ig)
              q2(ig)=b2(ig)
              q3(ig)=b3(ig)
              omx(ig)=oxo(ig)
              omy(ig)=oyo(ig)
              omz(ig)=ozo(ig)
              gcmx(ig)=gxo(ig)
              gcmy(ig)=gyo(ig)
              gcmz(ig)=gzo(ig)
              gvxx(ig)=gvxo(ig)
              gvyy(ig)=gvyo(ig)
              gvzz(ig)=gvzo(ig)
              
            enddo

          endif

        endif

c     operations for second stage only

        if(isw.eq.2)then
          
          jr=0
          do ig=igrp1,igrp2
            
c     new angular momenta and velocities
            
            opx=pt5*(-q1(ig)*p0(ig)+q0(ig)*p1(ig)+
     x        q3(ig)*p2(ig)-q2(ig)*p3(ig))
            opy=pt5*(-q2(ig)*p0(ig)-q3(ig)*p1(ig)+
     x        q0(ig)*p2(ig)+q1(ig)*p3(ig))
            opz=pt5*(-q3(ig)*p0(ig)+q2(ig)*p1(ig)-
     x        q1(ig)*p2(ig)+q0(ig)*p3(ig))
            
            id=lstgtp(ig)
            
            omx(ig)=opx*rotinx(id,2)
            omy(ig)=opy*rotiny(id,2)
            omz(ig)=opz*rotinz(id,2)
            
c     new rotation matrix
            
            call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)
            
            do j=1,numgsit(id)
              
              jr=jr+1
              i=lstrgd(jr)
              
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
          
c     merge velocities from all nodes
          
          if(mxnode.gt.1)call merge1
     x      (idnode,mxnode,natms,lstme,vxx,vyy,vzz,buffer)
          
c     correct constraint bond velocities using rattle
          
          if(ntcons.gt.0)then
            
            safe=.false.
            call rdrattle_v
     x        (safe,idnode,mxnode,natms,nscons,tolnce,tstep,
     x        dxx,dyy,dzz,txx,tyy,tzz,xxt,yyt,zzt)
            if(.not.safe)return
            
c     add constraint terms to stress tensor

            do i=1,9
              stress(i)=stress(i)+strcns(i)
            enddo

          endif
          
c     kinetic terms for stress tensor
        
          call kinstressf(ntfree,idnode,mxnode,strkin)        
          call kinstressg(ngrp,idnode,mxnode,strgrp)
        
c     calculate rigid body contribution to stress tensor

          call bodystress(idnode,mxnode,ngrp,vircom,strbod,dtx,dty,dtz)

c     add contributions to stress tensor
          
          do i=1,9
            stress(i)=stress(i)+strkin(i)+strgrp(i)+strbod(i)
          enddo
          
          do icyc=1,ncyc
            
c     integrate and apply nst thermostat
            
            call nstqscl_t
     x        (idnode,mxnode,ntfree,ngrp,mode,engfke,engtrn,engrot,
     x        temp,sigma,qstep,pmass,qmass,taut,chit,conint)
            
c     metadynamics shell thermostat
c     ====================================================
c     Must first merge update velocities as the core-shell
c     velocites are not distributed according to the same
c     rules.
c     ====================================================
            if(lmetadyn.and.keyshl.eq.1)then
              if(mxnode.gt.1)call merge
     x          (idnode,mxnode,natms,mxbuff,vxx,vyy,vzz,buffer)
              call nvtscale_shl
     x          (idnode,mxnode,ntshl,shlke,sigma_shl,hstep,qmass_shl,
     x          taut,chit_shl,conint)      
            endif  
            
c     integrate and apply nst barostat
            
            call nstqscl_p
     x        (idnode,mxnode,ntfree,ngrp,mode,fstep,pmass,chit,
     x        press,volm)
            
c     integrate and apply nst thermostat
            
            call nstqscl_t
     x        (idnode,mxnode,ntfree,ngrp,mode,engfke,engtrn,engrot,
     x        temp,sigma,qstep,pmass,qmass,taut,chit,conint)
          
c     metadynamics shell thermostat
c     ====================================================
c     Must first merge update velocities as the core-shell
c     velocites are not distributed according to the same
c     rules.
c     ====================================================
            if(lmetadyn.and.keyshl.eq.1)then
              if(mxnode.gt.1)call merge
     x          (idnode,mxnode,natms,mxbuff,vxx,vyy,vzz,buffer)
              call nvtscale_shl
     x          (idnode,mxnode,ntshl,shlke,sigma_shl,hstep,qmass_shl,
     x          taut,chit_shl,conint)      
            endif  
            
          enddo
          
          engke=engfke+engtrn
          
          if(mxnode.gt.1)call merge1
     x      (idnode,mxnode,natms,lstme,vxx,vyy,vzz,buffer)

c     conserved quantity less kinetic and potential energy terms

          chip2=sdot0(9,eta,eta)
          if(mode.eq.2)chip2=chip2-eta(1)**2
          consv=conint+0.5d0*qmass*chit**2+0.5d0*pmass*chip2+press*volm
          
c     metadynamics shell thermostat

          if(lmetadyn.and.keyshl.eq.1)then
            consv=consv+0.5d0*qmass_shl*chit_shl**2
          endif
          
c     end of second stage of velocity verlet algorithm
          
        endif
        
c     end of iteration cycle

      enddo

c     periodic boundary condition
          
      if(isw.eq.2)then
        
        call images(imcon,idnode,mxnode,natms,cell,xxx,yyy,zzz)
        call images(imcon,idnode,mxnode,ngrp,cell,gcmx,gcmy,gcmz)
        
c     merge position data
        
        if(mxnode.gt.1)then
          
          call merge(idnode,mxnode,natms,mxbuff,xxx,yyy,zzz,buffer)
          call merge(idnode,mxnode,ngrp,mxbuff,gcmx,gcmy,gcmz,buffer)
          
        endif
        
      endif
      
c     adjust long range corrections and number density
      
      elrc=elrc0*(volm0/volm)
      virlrc=virlrc0*(volm0/volm)
      do kk=1,ntpatm
        dens(kk)=dens0(kk)*(volm0/volm)
      enddo
      
c     deallocate working arrays
      
      deallocate(dtx,dty,dtz,stat=fail(1))
      
      if(ntcons.gt.0)then

        deallocate(dxx,dyy,dzz,txx,tyy,tzz,stat=fail(3))
        deallocate(xxt,yyt,zzt,dxt,dyt,dzt,stat=fail(4))
        deallocate(xxo,yyo,zzo,oxo,oyo,ozo,stat=fail(5))
        deallocate(vxo,vyo,vzo,b0,b1,b2,b3,stat=fail(6))
        deallocate(gxo,gyo,gzo,gvxo,gvyo,gvzo,stat=fail(7))

      endif
      
      return
      end subroutine nstqvv_h1
      
      end module vv_rotation1_module

