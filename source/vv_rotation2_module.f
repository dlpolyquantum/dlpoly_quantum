      module vv_rotation2_module
      
c***********************************************************************
c     
c     dl_poly module 2 for velocity verlet rotational integration 
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
      use vv_rotation1_module
      use utility_module
      
      contains
      
      subroutine qrattle_r
     x  (newstep,safe,lshmov,idnode,imcon,mxnode,natms,
     x  nscons,tolnce,tstep,vircon,dxx,dyy,dzz,dxt,dyt,
     x  dzt,txx,tyy,tzz,xxt,yyt,zzt,stresh)

c***********************************************************************
c     
c     dl_poly subroutine for applying bond constraint corrections after
c     atomic integration in the velocity verlet scheme. assumes rigid
c     bodies connected by constraints.  must be used in conjunction with
c     velocity verlet integration algorithm. note the iteration is
c     handled by the calling routine.
c     
c     copyright - daresbury laboratory
c     author    - w. smith february 2005
c     
c***********************************************************************

      implicit none
      
      logical safe,lshmov,newstep,newjob
      integer fail,idnode,imcon,mxnode,natms,nscons,i,j,k,ik
      real(8) tolnce,tstep,vircon,dxx,dyy,dzz,dxt,dyt,dzt,txx,tyy
      real(8) tzz,xxt,yyt,zzt,stresh,strs1,strs2,strs3,strs5,strs6
      real(8) strs9,tstep2,esig,dis2,tqa,tqb,gamma,dli,dlj

      real(8), allocatable :: esig1(:),ggx(:),ggy(:),ggz(:)

      dimension xxt(mxatms),yyt(mxatms),zzt(mxatms)
      dimension txx(mxatms),tyy(mxatms),tzz(mxatms)
      dimension dxx(mxcons),dyy(mxcons),dzz(mxcons)
      dimension dxt(mxcons),dyt(mxcons),dzt(mxcons)
      dimension stresh(9),tqa(3),tqb(3)

      save newjob,esig1,ggx,ggy,ggz

      data newjob/.true./,fail/0/
      
      if(newjob)then
        
        allocate (esig1(mxcons),ggx(mxcons),ggy(mxcons),ggz(mxcons),
     x    stat=fail)
        if(fail.ne.0)call error(idnode,1615)
        newjob=.false.

      endif

c     accumulators for stress tensor

      strs1=0.d0
      strs2=0.d0
      strs3=0.d0
      strs5=0.d0
      strs6=0.d0
      strs9=0.d0

c     timestep squared
      
      tstep2=tstep*tstep

c     constraint bond vectors are dxx,dyy,dzz (input)
      

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

        dis2=prmcon(listcon(k,1))**2
        esig1(k)=(dis2-(dxt(k)**2+dyt(k)**2+dzt(k)**2))/dis2
        esig=max(esig,abs(esig1(k)))

      enddo
      
c     global verification of convergence

      safe=(esig.lt.tolnce)
      
      if(mxnode.gt.1)call gstate(safe)

c     continue if any tolerances unsatisfied 
      
      if(.not.safe)then

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
          dis2=prmcon(listcon(k,1))**2

          if(newstep)then

            call pivot(1,i,k,ik,tqa,dxx,dyy,dzz)
            call pivot(1,j,k,ik,tqb,dxx,dyy,dzz)

            ggx(k)=tstep2*(tqa(1)+tqb(1))/dis2
            ggy(k)=tstep2*(tqa(2)+tqb(2))/dis2
            ggz(k)=tstep2*(tqa(3)+tqb(3))/dis2

          endif

c     constraint force parameter 
          
          gamma=esig1(k)/(dxt(k)*ggx(k)+dyt(k)*ggy(k)+dzt(k)*ggz(k))
          
c     accumulate bond virial
          
          vircon=vircon-gamma*(dxx(k)**2+dyy(k)**2+dzz(k)**2)

          strs1=strs1+gamma*dxx(k)*dxx(k)
          strs2=strs2+gamma*dxx(k)*dyy(k)
          strs3=strs3+gamma*dxx(k)*dzz(k)
          strs5=strs5+gamma*dyy(k)*dyy(k)
          strs6=strs6+gamma*dyy(k)*dzz(k)
          strs9=strs9+gamma*dzz(k)*dzz(k)

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

c     complete stress tensor
        
        stresh(1)=stresh(1)+strs1
        stresh(2)=stresh(2)+strs2
        stresh(3)=stresh(3)+strs3
        stresh(4)=stresh(4)+strs2
        stresh(5)=stresh(5)+strs5
        stresh(6)=stresh(6)+strs6
        stresh(7)=stresh(7)+strs3
        stresh(8)=stresh(8)+strs6
        stresh(9)=stresh(9)+strs9
        
c     splice force arrays across nodes

        if(mxnode.gt.1)then

          call splice 
     x      (idnode,natms,listme,listot,fxx,fyy,fzz,buffer)
          
        endif
        

      endif
      
      return
      end subroutine qrattle_r

      subroutine qrattle_v
     x  (newstep,safe,lshmov,idnode,mxnode,natms,
     x  nscons,tolnce,tstep,dxx,dyy,dzz,txx,tyy,tzz,
     x  xxt,yyt,zzt)

c***********************************************************************
c     
c     dl_poly subroutine for applying bond constraint corrections after
c     atomic integration in the velocity verlet scheme. assumes rigid
c     bodies connected by constraints.  must be used in conjunction with
c     velocity verlet integration algorithm. note the iteration is
c     handled by the calling routine.
c     
c     copyright - daresbury laboratory
c     author    - w. smith february 2005
c     
c***********************************************************************

      implicit none
      
      logical safe,lshmov,newstep,newjob
      integer fail,idnode,mxnode,natms,nscons,i,j,k,ik
      real(8) tolnce,tstep,dxx,dyy,dzz,txx,tyy,tzz,tqa,tqb
      real(8) xxt,yyt,zzt,tstep2,esig,gamma,dli,dlj
      real(8) tolvel

      real(8), allocatable :: esig2(:),hhx(:),hhy(:),hhz(:)

      dimension xxt(mxatms),yyt(mxatms),zzt(mxatms)
      dimension txx(mxatms),tyy(mxatms),tzz(mxatms)
      dimension dxx(mxcons),dyy(mxcons),dzz(mxcons)
      dimension tqa(3),tqb(3)

      save newjob,esig2,hhx,hhy,hhz

      data newjob/.true./,fail/0/
      
      if(newjob)then
        
        allocate (esig2(mxcons),hhx(mxcons),hhy(mxcons),hhz(mxcons),
     x    stat=fail)
        if(fail.ne.0)call error(idnode,1625)
        newjob=.false.

      endif

c     constraint bond vectors are dxx,dyy,dzz (input)
      
c     half timestep
      
      tstep2=tstep/2.d0

c     tolerance for velocity convergence

      tolvel=tolnce/tstep

c     calculate maximum error in constraint
      
      esig=0.d0
      do k=1,nscons

        i=listcon(k,2)
        j=listcon(k,3)
        esig2(k)=(dxx(k)*(vxx(i)-vxx(j))+dyy(k)*(vyy(i)-vyy(j))+
     x    dzz(k)*(vzz(i)-vzz(j)))
        esig=max(esig,abs(esig2(k)))

      enddo
      
c     global verification of convergence

      safe=(esig.lt.tolvel)
      
      if(mxnode.gt.1)then
        call gstate(safe)
      endif

c     continue if all tolerances satisfied else return to calling routine 
      
      if(.not.safe)then

c     initialise velocity correction arrays
        
        do i=1,natms
          
          xxt(i)=0.d0
          yyt(i)=0.d0
          zzt(i)=0.d0
          
        enddo

c     calculate constraint correction
        
        ik=0
        do k=1,nscons
          
          i=listcon(k,2)
          j=listcon(k,3)

          if(newstep)then

            call pivot(2,i,k,ik,tqa,dxx,dyy,dzz)
            call pivot(2,j,k,ik,tqb,dxx,dyy,dzz)

            hhx(k)=tstep2*(tqa(1)+tqb(1))
            hhy(k)=tstep2*(tqa(2)+tqb(2))
            hhz(k)=tstep2*(tqa(3)+tqb(3))

          endif

c     constraint force parameter 
          
          gamma=esig2(k)/(dxx(k)*hhx(k)+dyy(k)*hhy(k)+dzz(k)*hhz(k))

c     improved atomic force
          
          xxt(i)=xxt(i)-dxx(k)*gamma
          yyt(i)=yyt(i)-dyy(k)*gamma
          zzt(i)=zzt(i)-dzz(k)*gamma
          xxt(j)=xxt(j)+dxx(k)*gamma
          yyt(j)=yyt(j)+dyy(k)*gamma
          zzt(j)=zzt(j)+dzz(k)*gamma
          
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

          fxx(i)=fxx(i)+dli*xxt(i)
          fyy(i)=fyy(i)+dli*yyt(i)
          fzz(i)=fzz(i)+dli*zzt(i)
          fxx(j)=fxx(j)+dlj*xxt(j)
          fyy(j)=fyy(j)+dlj*yyt(j)
          fzz(j)=fzz(j)+dlj*zzt(j)
          
        enddo

c     splice force arrays across nodes
        
        if(mxnode.gt.1)then

          call splice 
     x      (idnode,natms,listme,listot,fxx,fyy,fzz,buffer)
          
        endif

      endif
      
      return
      end subroutine qrattle_v

      subroutine pivot(k,i,kk,ik,tqq,dxx,dyy,dzz)

c***********************************************************************
c     
c     dl_poly subroutine for computing pivot vector for velocity
c     corrections to bonds between rigid bodies
c     must be used in conjunction with qrattle routines:
c     if k=1 - use with qrattle_r
c     if k=2 - use with qrattle_v
c     
c     copyright - daresbury laboratory
c     author    - w. smith february 2005
c     
c***********************************************************************

      implicit none

      integer k,i,kk,ik,ig,id,jj
      real(8) xxa,yya,zza,tax,tay,taz,trx,try,trz,vix,viy,viz
      real(8) rot(9),tqq(3),dxx(mxcons),dyy(mxcons),dzz(mxcons)

      ig=lstbod(i)

      if(ig.eq.0)then

c     atoms in constraint bonds
        
        tqq(1)=dxx(kk)*rmass(i)
        tqq(2)=dyy(kk)*rmass(i)
        tqq(3)=dzz(kk)*rmass(i)
        
      else

c     terms for rigid body atoms

        ik=ik+1
        id=lstgtp(ig)

        tqq(1)=dxx(kk)/gmass(id)
        tqq(2)=dyy(kk)/gmass(id)
        tqq(3)=dzz(kk)/gmass(id)

        call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)
        
        jj=lstcsit(ik)

c     com-site vector in lab frame

        xxa=(gxx(id,jj)*rot(1)+gyy(id,jj)*rot(2)+gzz(id,jj)*rot(3))
        yya=(gxx(id,jj)*rot(4)+gyy(id,jj)*rot(5)+gzz(id,jj)*rot(6))
        zza=(gxx(id,jj)*rot(7)+gyy(id,jj)*rot(8)+gzz(id,jj)*rot(9))

c     cross product of com-site vector and interatomic vector
        
        tax=yya*dzz(kk)-zza*dyy(kk)
        tay=zza*dxx(kk)-xxa*dzz(kk)
        taz=xxa*dyy(kk)-yya*dxx(kk)

c     transform to body fixed frame
        
        trx=(tax*rot(1)+tay*rot(4)+taz*rot(7))*rotinx(id,2)
        try=(tax*rot(2)+tay*rot(5)+taz*rot(8))*rotiny(id,2)
        trz=(tax*rot(3)+tay*rot(6)+taz*rot(9))*rotinz(id,2)

        if(k.eq.1)then

c     direction of induced velocites in body frame
          
          vix=try*gzz(id,jj)-trz*gyy(id,jj)
          viy=trz*gxx(id,jj)-trx*gzz(id,jj)
          viz=trx*gyy(id,jj)-try*gxx(id,jj)
          
c     transform to lab frame
          
          tqq(1)=tqq(1)+vix*rot(1)+viy*rot(2)+viz*rot(3)
          tqq(2)=tqq(2)+vix*rot(4)+viy*rot(5)+viz*rot(6)
          tqq(3)=tqq(3)+vix*rot(7)+viy*rot(8)+viz*rot(9)

        elseif(k.eq.2)then

c     transform to lab frame

          tqq(1)=tqq(1)+trx*rot(1)+try*rot(2)+trz*rot(3)
          tqq(2)=tqq(2)+trx*rot(4)+try*rot(5)+trz*rot(6)
          tqq(3)=tqq(3)+trx*rot(7)+try*rot(8)+trz*rot(9)
          
        endif

      endif

      return
      end subroutine pivot

      subroutine nveqvv_2
     x  (safe,lshmov,isw,imcon,idnode,mxnode,natms,ngrp,nscons,
     x  ntcons,ntfree,engke,engrot,tolnce,tstep,vircom,vircon)
      
c***********************************************************************
c     
c     dlpoly subroutine to integrate rigid body equations of motion
c     using the symplectic no_squish quaternion algorithm of 
c     miller et al j.chem.phys 116 (2002) 8649
c     plus bond constraints -  including rigid body sites linked
c     by constraint sites (qrattle algorithm)
c     
c     parallel replicated data version : block data
c     
c     omx,omy,omz=angular velocity in body fixed frame (principal axes)
c     rotinx,y,z =rotational inertia in body fixed frame
c     
c     copyright daresbury laboratory
c     author      w.smith   jan 2005
c     amended     w.smith   feb 2005: qrattle added
c     
c**********************************************************************

      implicit none

      logical newstep,safe,lshmov,newjob

      integer isw,imcon,idnode,mxnode,natms,ngrp,nscons,ntcons
      integer ntfree,i,j,k,igrp1,igrp2,igrp,ifre1,ifre2,jg,jr
      integer id,ifre,jrs,icyc,mxshk,idum,ig

      real(8) engke,engrot,tolnce,tstep,vircom,vircon
      real(8) engtrn
      real(8) vaa,vbb,vcc,opx,opy,opz,ftx,fty,ftz
      real(8) fmx,fmy,fmz,tqx,tqy,tqz,tq0,tq1,tq2,tq3

      integer, parameter :: nnn=13
      real(8), parameter :: pt5=0.5d0

      integer fail(nnn)
      real(8) rot(9),strkin(9),strgrp(9)

      real(8), allocatable :: gvxo(:),gvyo(:),gvzo(:)
      real(8), allocatable :: vxo(:),vyo(:),vzo(:)
      real(8), allocatable :: b0(:),b1(:),b2(:),b3(:)
      real(8), allocatable :: c0(:),c1(:),c2(:),c3(:)
      real(8), allocatable :: p0(:),p1(:),p2(:),p3(:)
      real(8), allocatable :: xxo(:),yyo(:),zzo(:)
      real(8), allocatable :: gxo(:),gyo(:),gzo(:)
      real(8), allocatable :: dxx(:),dyy(:),dzz(:)
      real(8), allocatable :: txx(:),tyy(:),tzz(:)
      real(8), allocatable :: xxt(:),yyt(:),zzt(:)
      real(8), allocatable :: dxt(:),dyt(:),dzt(:)
      real(8), allocatable :: dtx(:),dty(:),dtz(:)
      real(8), allocatable :: fxo(:),fyo(:),fzo(:)

      save newstep,newjob,p0,p1,p2,p3,ifre1,ifre2,igrp1,igrp2

      data newjob/.true./
      
c     set array allocation error flags

      do i=1,nnn
        fail(i)=0
      enddo

c     assign initial parameters

      if(newjob)then
        
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

        allocate(p0(mxgrp),p1(mxgrp),p2(mxgrp),p3(mxgrp),stat=fail(1))
        newjob=.false.

      endif

c     allocate working arrays

      allocate(gvxo(msgrp),gvyo(msgrp),gvzo(msgrp),stat=fail(2))
      allocate(vxo(msatms),vyo(msatms),vzo(msatms),stat=fail(3))
      allocate(b0(msgrp),b1(msgrp),b2(msgrp),b3(msgrp),
     x  stat=fail(4))
      allocate(c0(msgrp),c1(msgrp),c2(msgrp),c3(msgrp),
     x  stat=fail(5))
      if(isw.eq.1)then

        allocate(xxo(mxatms),yyo(mxatms),zzo(mxatms),stat=fail(6))
        allocate(gxo(msgrp),gyo(msgrp),gzo(msgrp),stat=fail(7))

      endif
      if(ntcons.gt.0)then

        allocate(dxx(mxcons),dyy(mxcons),dzz(mxcons),stat=fail(8))
        allocate(txx(mxatms),tyy(mxatms),tzz(mxatms),stat=fail(9))
        allocate(xxt(mxatms),yyt(mxatms),zzt(mxatms),stat=fail(10))
        allocate(dxt(mxcons),dyt(mxcons),dzt(mxcons),stat=fail(11))

      endif
      allocate(dtx(mxatms),dty(mxatms),dtz(mxatms),stat=fail(12))
      if(isw.eq.2)then
        allocate(fxo(mxatms),fyo(mxatms),fzo(mxatms),stat=fail(13))
      endif

c     check array allocation error flags

      do i=1,nnn
        if(fail(i).gt.0)call error(idnode,2110)
      enddo

c     initialise constraint virial

      if(isw.eq.1)then

        vircon=0.d0

        do i=1,9
          strcns(i)=0.d0
        enddo

      endif

c     construct current bond vectors
      
      if(ntcons.gt.0)then

        do k=1,nscons
          
          i=listcon(k,2)
          j=listcon(k,3)
          dxx(k)=xxx(i)-xxx(j)
          dyy(k)=yyy(i)-yyy(j)
          dzz(k)=zzz(i)-zzz(j)
          
        enddo

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

c     store key config data at start of time step

      if(isw.eq.1)then

c     atom positions

        do i=1,natms
          
          xxo(i)=xxx(i)
          yyo(i)=yyy(i)
          zzo(i)=zzz(i)

        enddo

c     rigid body positions

        j=0
        do i=igrp1,igrp2
          
          j=j+1
          gxo(j)=gcmx(i)
          gyo(j)=gcmy(i)
          gzo(j)=gcmz(i)

        enddo

      endif

c     store free atom velocities

      j=0
      do ifre=ifre1,ifre2
        
        j=j+1
        i=lstfre(ifre)
        vxo(j)=vxx(i)
        vyo(j)=vyy(i)
        vzo(j)=vzz(i)
        
      enddo
      
c     store rigid body quaternions, momenta and cartesian velocities

      j=0
      do i=igrp1,igrp2

        j=j+1
        b0(j)=q0(i)
        b1(j)=q1(i)
        b2(j)=q2(i)
        b3(j)=q3(i)
        c0(j)=p0(i)
        c1(j)=p1(i)
        c2(j)=p2(i)
        c3(j)=p3(i)
        gvxo(j)=gvxx(i)
        gvyo(j)=gvyy(i)
        gvzo(j)=gvzz(i)

      enddo

c     store forces if isw = 2
      
      if(isw.eq.2)then

        do i=1,natms

          fxo(i)=fxx(i)
          fyo(i)=fyy(i)
          fzo(i)=fzz(i)

        enddo
      
      endif

c     -------------- start of shake iteration cycle -------------------
      
      icyc=0
      mxshk=1
      safe=.false.
      newstep=.true.
      if(ntcons.gt.0)mxshk=mxshak
      do while(.not.safe.and.icyc.lt.mxshk)
        
        icyc=icyc+1

c     update velocities of free atoms 1/2 timestep

        j=0
        do ifre=ifre1,ifre2
          
          j=j+1
          i=lstfre(ifre)
          vxx(i)=vxo(j)+(pt5*tstep*rmass(i))*fxx(i)
          vyy(i)=vyo(j)+(pt5*tstep*rmass(i))*fyy(i)
          vzz(i)=vzo(j)+(pt5*tstep*rmass(i))*fzz(i)
          
        enddo
        
c     *************  rigid body motion ****************************

c     operations common to first and second stages

        jg=0
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
          
          jg=jg+1
          call getrotmat(b0(jg),b1(jg),b2(jg),b3(jg),rot)

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

          tq0=2.0d0*(-b1(jg)*tqx-b2(jg)*tqy-b3(jg)*tqz)
          tq1=2.0d0*( b0(jg)*tqx-b3(jg)*tqy+b2(jg)*tqz)
          tq2=2.0d0*( b3(jg)*tqx+b0(jg)*tqy-b1(jg)*tqz)
          tq3=2.0d0*(-b2(jg)*tqx+b1(jg)*tqy+b0(jg)*tqz)

c     update quaternion momentum by half timestep

          p0(ig)=c0(jg)+tq0*pt5*tstep
          p1(ig)=c1(jg)+tq1*pt5*tstep
          p2(ig)=c2(jg)+tq2*pt5*tstep
          p3(ig)=c3(jg)+tq3*pt5*tstep

c     update centre of mass velocity by half timestep

          gvxx(ig)=gvxo(jg)+fmx*pt5*tstep/gmass(id)
          gvyy(ig)=gvyo(jg)+fmy*pt5*tstep/gmass(id)
          gvzz(ig)=gvzo(jg)+fmz*pt5*tstep/gmass(id)
          
        enddo

c     first stage of velocity verlet algorithm

        if(isw.eq.1)then

          jg=0
          do ig=igrp1,igrp2

            jg=jg+1

c     update centre of mass position by full time step
            
            gcmx(ig)=gxo(jg)+tstep*gvxx(ig)
            gcmy(ig)=gyo(jg)+tstep*gvyy(ig)
            gcmz(ig)=gzo(jg)+tstep*gvzz(ig)
            
c     calculate rotation of rigid groups: nosquish algorithm
            
            q0(ig)=b0(jg)
            q1(ig)=b1(jg)
            q2(ig)=b2(jg)
            q3(ig)=b3(jg)
            call nosquish(ig,tstep,q0,q1,q2,q3,p0,p1,p2,p3)
            
          enddo

c     new atomic positions for atoms in rigid bodies - relative to com
          
          k=0
          do ig=igrp1,igrp2
            
c     new rotational matrix
            
            call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)
            
            id=lstgtp(ig)
            do j=1,numgsit(id)
              
              k=k+1
              i=lstme(k)
              xxx(i)=rot(1)*gxx(id,j)+rot(2)*gyy(id,j)+
     x          rot(3)*gzz(id,j)+gcmx(ig)
              yyy(i)=rot(4)*gxx(id,j)+rot(5)*gyy(id,j)+
     x          rot(6)*gzz(id,j)+gcmy(ig)
              zzz(i)=rot(7)*gxx(id,j)+rot(8)*gyy(id,j)+
     x          rot(9)*gzz(id,j)+gcmz(ig)
              
            enddo
            
          enddo
          
c     update positions of free particles to full time step
          
          do ifre=ifre1,ifre2
            
            i=lstfre(ifre)
            xxx(i)=xxo(i)+tstep*vxx(i)
            yyy(i)=yyo(i)+tstep*vyy(i)
            zzz(i)=zzo(i)+tstep*vzz(i)
            
          enddo

c     merge free atom positions

          if(mxnode.gt.1)call merge1
     x      (idnode,mxnode,natms,lstme,xxx,yyy,zzz,buffer)

c     apply rattle corrections to bond constraints
          
          if(ntcons.gt.0)then
            
            if(mxnode.gt.1)call merge1
     x        (idnode,mxnode,natms,lstme,vxx,vyy,vzz,buffer)

            call qrattle_r
     x        (newstep,safe,lshmov,idnode,imcon,mxnode,natms,
     x        nscons,tolnce,tstep,vircon,dxx,dyy,dzz,dxt,dyt,
     x        dzt,txx,tyy,tzz,xxt,yyt,zzt,strcns)

          endif
          
c     end of first stage 

        endif

        if(isw.eq.2)then
          
c     second stage of velocity verlet algorithm
          
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

c     new rotational matrix
            
            call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)

c     new atomic velocites
            
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
          
          if(ntcons.gt.0)then
            
            if(mxnode.gt.1)then

              call merge1
     x          (idnode,mxnode,natms,lstme,vxx,vyy,vzz,buffer)
              call merge1
     x          (idnode,mxnode,natms,lstme,fxx,fyy,fzz,buffer)
              
            endif

c     correct constraint bond velocities using rattle

            call qrattle_v
     x        (newstep,safe,lshmov,idnode,mxnode,natms,
     x        nscons,tolnce,tstep,dxx,dyy,dzz,txx,tyy,tzz,
     x        xxt,yyt,zzt)

          endif
          
c     end of second stage
          
        endif
        
        newstep=.false.
        
      enddo

c     check shake convergence

      if(.not.safe)call error(idnode,105)

c     sum constraint virial and stress across processors
      
      if(mxnode.gt.1.and.isw.eq.1)then
        
        buffer(1)=vircon
        call gdsum(buffer(1),1,buffer(2))
        vircon=buffer(1)
        call gdsum(strcns,9,buffer)
        
      endif

c     -------------- end of shake iteration cycle -------------------
      
c     calculate kinetic energy
      
      if(isw.eq.2)then

        engke=getkinf(ntfree,idnode,mxnode)      
        call getking(ngrp,idnode,mxnode,engtrn,engrot)

        engke=engke+engtrn
        
c     kinetic contribution to stress tensor
        
        call kinstressf(ntfree,idnode,mxnode,strkin)
        call kinstressg(ngrp,idnode,mxnode,strgrp)
        
c     rigid body contribution to stress tensor

        call bodystress(idnode,mxnode,ngrp,vircom,strbod,dtx,dty,dtz)

c     add contributions to stress tensor
        
        do i=1,9
          stress(i)=stress(i)+strbod(i)+strcns(i)+strkin(i)+strgrp(i)
        enddo

      endif
      
      if(mxnode.gt.1)then

c     merge new group coordinates and velocities

        if(isw.eq.1)
     x    call merge(idnode,mxnode,ngrp,mxbuff,gcmx,gcmy,gcmz,buffer)
        call merge(idnode,mxnode,ngrp,mxbuff,gvxx,gvyy,gvzz,buffer)
        
c     merge new atomic coordinates and velocities
        
        if(isw.eq.1)
     x    call merge1(idnode,mxnode,natms,lstme,xxx,yyy,zzz,buffer)
        call merge1(idnode,mxnode,natms,lstme,vxx,vyy,vzz,buffer)
        
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
      
c     restore forces if isw = 2

      if(isw.eq.2)then

        do i=1,natms

          fxx(i)=fxo(i)
          fyy(i)=fyo(i)
          fzz(i)=fzo(i)

        enddo
      
      endif

c     deallocate working arrays

      deallocate(gvxo,gvyo,gvzo,vxo,vyo,vzo,stat=fail(1))
      deallocate(b0,b1,b2,b3,c0,c1,c2,c3,stat=fail(2))
      deallocate(dtx,dty,dtz,stat=fail(3))
      if(isw.eq.1)then
        deallocate(xxo,yyo,zzo,gxo,gyo,gzo,stat=fail(4))
      endif
      if(ntcons.gt.0)then

        deallocate(dxx,dyy,dzz,txx,tyy,tzz,stat=fail(5))
        deallocate(dxt,dyt,dzt,xxt,yyt,zzt,stat=fail(6))

      endif
      if(isw.eq.2)then
        deallocate(fxo,fyo,fzo,stat=fail(7))
      endif
      
      return
      end subroutine nveqvv_2

      subroutine nvtqvv_b2
     x  (safe,lshmov,isw,imcon,idnode,mxnode,natms,ngrp,nscons,
     x  ntcons,ntfree,engke,engrot,taut,sigma,tolnce,tstep,
     x  vircom,vircon)
      
c***********************************************************************
c     
c     dlpoly subroutine to integrate rigid body equations of motion
c     using the symplectic no_squish quaternion algorithm of 
c     miller et al j.chem.phys 116 (2002) 8649
c     plus bond constraints -  including rigid body sites linked
c     by constraint sites (qrattle algorithm)
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

      logical newstep,safe,lshmov,newjob

      integer isw,imcon,idnode,mxnode,natms,ngrp,nscons,ntcons
      integer ntfree,i,j,k,igrp1,igrp2,igrp,ifre1,ifre2,jg,jr
      integer id,ifre,jrs,icyc,mxshk,idum,ig

      real(8) engke,engrot,tolnce,tstep,vircom,vircon,engtke
      real(8) engtrn
      real(8) vaa,vbb,vcc,opx,opy,opz,ftx,fty,ftz
      real(8) fmx,fmy,fmz,tqx,tqy,tqz,tq0,tq1,tq2,tq3,taut,sigma,chit

      integer, parameter :: nnn=13
      real(8), parameter :: pt5=0.5d0

      integer fail(nnn)
      real(8) rot(9),strkin(9),strgrp(9)

      real(8), allocatable :: gvxo(:),gvyo(:),gvzo(:)
      real(8), allocatable :: vxo(:),vyo(:),vzo(:)
      real(8), allocatable :: b0(:),b1(:),b2(:),b3(:)
      real(8), allocatable :: c0(:),c1(:),c2(:),c3(:)
      real(8), allocatable :: p0(:),p1(:),p2(:),p3(:)
      real(8), allocatable :: xxo(:),yyo(:),zzo(:)
      real(8), allocatable :: gxo(:),gyo(:),gzo(:)
      real(8), allocatable :: dxx(:),dyy(:),dzz(:)
      real(8), allocatable :: txx(:),tyy(:),tzz(:)
      real(8), allocatable :: xxt(:),yyt(:),zzt(:)
      real(8), allocatable :: dxt(:),dyt(:),dzt(:)
      real(8), allocatable :: dtx(:),dty(:),dtz(:)
      real(8), allocatable :: fxo(:),fyo(:),fzo(:)

      save newstep,newjob,p0,p1,p2,p3,ifre1,ifre2,igrp1,igrp2

      data newjob/.true./
      
c     set array alocation error flags

      do i=1,nnn
        fail(i)=0
      enddo

c     assign initial parameters

      if(newjob)then

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

        allocate(p0(mxgrp),p1(mxgrp),p2(mxgrp),p3(mxgrp),stat=fail(1))
        newjob=.false.

      endif

c     allocate working arrays

      allocate(gvxo(msgrp),gvyo(msgrp),gvzo(msgrp),stat=fail(2))
      allocate(vxo(msatms),vyo(msatms),vzo(mxatms),stat=fail(3))
      allocate(b0(msgrp),b1(msgrp),b2(msgrp),b3(msgrp),
     x  stat=fail(4))
      allocate(c0(msgrp),c1(msgrp),c2(msgrp),c3(msgrp),
     x  stat=fail(5))
      if(isw.eq.1)then

        allocate(xxo(mxatms),yyo(mxatms),zzo(mxatms),stat=fail(6))
        allocate(gxo(msgrp),gyo(msgrp),gzo(msgrp),stat=fail(7))

      endif
      if(ntcons.gt.0)then

        allocate(dxx(mxcons),dyy(mxcons),dzz(mxcons),stat=fail(8))
        allocate(txx(mxatms),tyy(mxatms),tzz(mxatms),stat=fail(9))
        allocate(xxt(mxatms),yyt(mxatms),zzt(mxatms),stat=fail(10))
        allocate(dxt(mxcons),dyt(mxcons),dzt(mxcons),stat=fail(11))

      endif
      allocate(dtx(mxatms),dty(mxatms),dtz(mxatms),stat=fail(12))
      if(isw.eq.2)then
        allocate(fxo(mxatms),fyo(mxatms),fzo(mxatms),stat=fail(13))
      endif

c     check array allocation error flags

      do i=1,nnn
        if(fail(i).gt.0)call error(idnode,2130)
      enddo

c     initialise constraint virial

      if(isw.eq.1)then

        vircon=0.d0

        do i=1,9
          strcns(i)=0.d0
        enddo

      endif

c     construct current bond vectors
      
      if(ntcons.gt.0)then

        do k=1,nscons
          
          i=listcon(k,2)
          j=listcon(k,3)
          dxx(k)=xxx(i)-xxx(j)
          dyy(k)=yyy(i)-yyy(j)
          dzz(k)=zzz(i)-zzz(j)
          
        enddo

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

c     store key config data at start of time step

      if(isw.eq.1)then

c     atom positions

        do i=1,natms
          
          xxo(i)=xxx(i)
          yyo(i)=yyy(i)
          zzo(i)=zzz(i)
          
        enddo

c     rigid body positions

        j=0
        do i=igrp1,igrp2
          
          j=j+1
          gxo(j)=gcmx(i)
          gyo(j)=gcmy(i)
          gzo(j)=gcmz(i)

        enddo

      endif

c     store free atom atom velocities

      j=0
      do ifre=ifre1,ifre2
        
        j=j+1
        i=lstfre(ifre)
        vxo(j)=vxx(i)
        vyo(j)=vyy(i)
        vzo(j)=vzz(i)
        
      enddo

c     store rigid body quaternions, momenta and cartesian velocities

      j=0
      do i=igrp1,igrp2

        j=j+1
        b0(j)=q0(i)
        b1(j)=q1(i)
        b2(j)=q2(i)
        b3(j)=q3(i)
        c0(j)=p0(i)
        c1(j)=p1(i)
        c2(j)=p2(i)
        c3(j)=p3(i)
        gvxo(j)=gvxx(i)
        gvyo(j)=gvyy(i)
        gvzo(j)=gvzz(i)

      enddo

c     store forces if isw = 2

      if(isw.eq.2)then

        do i=1,natms

          fxo(i)=fxx(i)
          fyo(i)=fyy(i)
          fzo(i)=fzz(i)

        enddo
      
      endif

c     -------------- start of shake iteration cycle -------------------
      
      icyc=0
      mxshk=1
      safe=.false.
      newstep=.true.
      if(ntcons.gt.0)mxshk=mxshak
      do while(.not.safe.and.icyc.lt.mxshk)
        
        icyc=icyc+1
        
c     update velocities of free atoms 1/2 timestep

        j=0
        do ifre=ifre1,ifre2

          j=j+1
          i=lstfre(ifre)
          vxx(i)=vxo(j)+(pt5*tstep*rmass(i))*fxx(i)
          vyy(i)=vyo(j)+(pt5*tstep*rmass(i))*fyy(i)
          vzz(i)=vzo(j)+(pt5*tstep*rmass(i))*fzz(i)
          
        enddo

c     *************  rigid body motion ****************************

c     operations common to first and second stages

        jg=0
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
          
          jg=jg+1
          call getrotmat(b0(jg),b1(jg),b2(jg),b3(jg),rot)

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

          tq0=2.0d0*(-b1(jg)*tqx-b2(jg)*tqy-b3(jg)*tqz)
          tq1=2.0d0*( b0(jg)*tqx-b3(jg)*tqy+b2(jg)*tqz)
          tq2=2.0d0*( b3(jg)*tqx+b0(jg)*tqy-b1(jg)*tqz)
          tq3=2.0d0*(-b2(jg)*tqx+b1(jg)*tqy+b0(jg)*tqz)

c     update quaternion momentum by half timestep

          p0(ig)=c0(jg)+tq0*pt5*tstep
          p1(ig)=c1(jg)+tq1*pt5*tstep
          p2(ig)=c2(jg)+tq2*pt5*tstep
          p3(ig)=c3(jg)+tq3*pt5*tstep

c     update centre of mass velocity by half timestep

          gvxx(ig)=gvxo(jg)+fmx*pt5*tstep/gmass(id)
          gvyy(ig)=gvyo(jg)+fmy*pt5*tstep/gmass(id)
          gvzz(ig)=gvzo(jg)+fmz*pt5*tstep/gmass(id)
          
        enddo

c     first stage of velocity verlet algorithm

        if(isw.eq.1)then

          jg=0
          do ig=igrp1,igrp2

            jg=jg+1

c     update centre of mass position by full time step
            
            gcmx(ig)=gxo(jg)+tstep*gvxx(ig)
            gcmy(ig)=gyo(jg)+tstep*gvyy(ig)
            gcmz(ig)=gzo(jg)+tstep*gvzz(ig)
            
c     calculate rotation of rigid groups: nosquish algorithm
            
            q0(ig)=b0(jg)
            q1(ig)=b1(jg)
            q2(ig)=b2(jg)
            q3(ig)=b3(jg)
            call nosquish(ig,tstep,q0,q1,q2,q3,p0,p1,p2,p3)
            
          enddo

c     new atomic positions for atoms in rigid bodies - relative to com
          
          k=0
          do ig=igrp1,igrp2
            
            call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)

            id=lstgtp(ig)
            do j=1,numgsit(id)
              
              k=k+1
              i=lstme(k)
              xxx(i)=rot(1)*gxx(id,j)+rot(2)*gyy(id,j)+
     x          rot(3)*gzz(id,j)+gcmx(ig)
              yyy(i)=rot(4)*gxx(id,j)+rot(5)*gyy(id,j)+
     x          rot(6)*gzz(id,j)+gcmy(ig)
              zzz(i)=rot(7)*gxx(id,j)+rot(8)*gyy(id,j)+
     x          rot(9)*gzz(id,j)+gcmz(ig)
              
            enddo
            
          enddo
          
c     update positions of free particles to full time step
          
          do ifre=ifre1,ifre2
            
            i=lstfre(ifre)
            xxx(i)=xxo(i)+tstep*vxx(i)
            yyy(i)=yyo(i)+tstep*vyy(i)
            zzz(i)=zzo(i)+tstep*vzz(i)
            
          enddo
          
c     merge atom positions

          if(mxnode.gt.1)call merge1
     x      (idnode,mxnode,natms,lstme,xxx,yyy,zzz,buffer)

c     apply rattle corrections to bond constraints
          
          if(ntcons.gt.0)then
            
            if(mxnode.gt.1)call merge1
     x          (idnode,mxnode,natms,lstme,vxx,vyy,vzz,buffer)

            call qrattle_r
     x        (newstep,safe,lshmov,idnode,imcon,mxnode,natms,
     x        nscons,tolnce,tstep,vircon,dxx,dyy,dzz,dxt,dyt,
     x        dzt,txx,tyy,tzz,xxt,yyt,zzt,strcns)

          endif
          
c     end of first stage 

        else
          
c     second stage of velocity verlet algorithm
          
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

c     new rotational matrix
            
            call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)

c     new atomic velocites
            
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
          
          if(ntcons.gt.0)then
            
            if(mxnode.gt.1)then

              call merge1
     x          (idnode,mxnode,natms,lstme,vxx,vyy,vzz,buffer)
              call merge1
     x          (idnode,mxnode,natms,lstme,fxx,fyy,fzz,buffer)
              
            endif

c     correct constraint bond velocities using rattle

            call qrattle_v
     x        (newstep,safe,lshmov,idnode,mxnode,natms,
     x        nscons,tolnce,tstep,dxx,dyy,dzz,txx,tyy,tzz,
     x        xxt,yyt,zzt)

          endif
          
c     end of second stage
          
        endif
        
        newstep=.false.
        
      enddo

c     check shake convergence

      if(.not.safe)call error(idnode,105)

c     sum constraint virial and stress across processors
      
      if(mxnode.gt.1.and.isw.eq.1)then
        
        buffer(1)=vircon
        call gdsum(buffer(1),1,buffer(2))
        vircon=buffer(1)
        call gdsum(strcns,9,buffer)
        
      endif

c     -------------- end of shake iteration cycle -------------------
      
c     calculate kinetic energy
      
      if(isw.eq.2)then

        engke=getkinf(ntfree,idnode,mxnode)      
        call getking(ngrp,idnode,mxnode,engtrn,engrot)
        engtke=engke+engtrn+engrot
        engke=engke+engtrn
        
c     apply Berendsen thermostat - taut is the relaxation time
        
        chit=sqrt(1.d0+tstep/taut*(sigma/engtke-1.d0))

        engke=engke*chit**2
        engtrn=engtrn*chit**2
        engrot=engrot*chit**2

c     thermostat velocities of free particles
        
        do ifre=ifre1,ifre2
          
          i=lstfre(ifre)
          if(lstfrz(i).ne.0)then
            
            vxx(i)=chit*vxx(i)
            vyy(i)=chit*vyy(i)
            vzz(i)=chit*vzz(i)
            
          endif

        enddo

c     thermostat rigid body velocities

        do ig=igrp1,igrp2
          
          id=lstgtp(ig)
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
        
c     rigid body contribution to stress tensor

        call bodystress(idnode,mxnode,ngrp,vircom,strbod,dtx,dty,dtz)

c     add contributions to stress tensor
        
        do i=1,9
          stress(i)=stress(i)+strbod(i)+strcns(i)+strkin(i)+strgrp(i)
        enddo

      endif
      
      if(mxnode.gt.1)then

c     merge new group coordinates and velocities

        if(isw.eq.1)
     x    call merge(idnode,mxnode,ngrp,mxbuff,gcmx,gcmy,gcmz,buffer)
        call merge(idnode,mxnode,ngrp,mxbuff,gvxx,gvyy,gvzz,buffer)
        
c     merge new atomic coordinates and velocities
        
        if(isw.eq.1)
     x    call merge1(idnode,mxnode,natms,lstme,xxx,yyy,zzz,buffer)
        call merge1(idnode,mxnode,natms,lstme,vxx,vyy,vzz,buffer)
        
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
      
c     restore forces if isw = 2

      if(isw.eq.2)then

        do i=1,natms

          fxx(i)=fxo(i)
          fyy(i)=fyo(i)
          fzz(i)=fzo(i)

        enddo
      
      endif

c     deallocate working arrays

      deallocate(gvxo,gvyo,gvzo,vxo,vyo,vzo,stat=fail(1))
      deallocate(b0,b1,b2,b3,c0,c1,c2,c3,stat=fail(2))
      deallocate(dtx,dty,dtz,stat=fail(3))
      if(isw.eq.1)then
        deallocate(xxo,yyo,zzo,gxo,gyo,gzo,stat=fail(4))
      endif
      if(ntcons.gt.0)then

        deallocate(dxx,dyy,dzz,txx,tyy,tzz,stat=fail(5))
        deallocate(dxt,dyt,dzt,xxt,yyt,zzt,stat=fail(6))

      endif
      if(isw.eq.2)then
        deallocate(fxo,fyo,fzo,stat=fail(7))
      endif
      
      return
      end subroutine nvtqvv_b2

      subroutine nvtqvv_h2
     x  (safe,lshmov,isw,imcon,idnode,mxnode,natms,ngrp,nscons,
     x  ntcons,ntfree,ntshl,keyshl,chit,consv,conint,engke,engrot,
     x  taut,sigma,tolnce,tstep,vircom,vircon,chit_shl,sigma_shl)
      
c***********************************************************************
c     
c     dlpoly subroutine to integrate rigid body equations of motion
c     using the symplectic no_squish quaternion algorithm of 
c     miller et al j.chem.phys 116 (2002) 8649
c     plus bond constraints -  including rigid body sites linked
c     by constraint sites (qrattle algorithm)
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

      logical newstep,safe,lshmov,newjob

      integer isw,imcon,idnode,mxnode,natms,ngrp,nscons,ntcons
      integer ntfree,i,j,k,igrp1,igrp2,igrp,ifre1,ifre2,jg,jr
      integer id,ifre,jrs,icyc,mxshk,idum,ig

      real(8) engke,engrot,tolnce,tstep,vircom,vircon
      real(8) engtrn,vaa,vbb,vcc,opx,opy,opz,ftx,fty,ftz
      real(8) fmx,fmy,fmz,tqx,tqy,tqz,tq0,tq1,tq2,tq3,engfke
      real(8) taut,sigma,chit,hstep,qmass,conint,consv

      integer, parameter :: nnn=13
      real(8), parameter :: pt5=0.5d0

      integer fail(nnn)
      real(8) rot(9),strkin(9),strgrp(9)

      real(8), allocatable :: gvxo(:),gvyo(:),gvzo(:)
      real(8), allocatable :: vxo(:),vyo(:),vzo(:)
      real(8), allocatable :: b0(:),b1(:),b2(:),b3(:)
      real(8), allocatable :: c0(:),c1(:),c2(:),c3(:)
      real(8), allocatable :: p0(:),p1(:),p2(:),p3(:)
      real(8), allocatable :: xxo(:),yyo(:),zzo(:)
      real(8), allocatable :: gxo(:),gyo(:),gzo(:)
      real(8), allocatable :: dxx(:),dyy(:),dzz(:)
      real(8), allocatable :: txx(:),tyy(:),tzz(:)
      real(8), allocatable :: xxt(:),yyt(:),zzt(:)
      real(8), allocatable :: dxt(:),dyt(:),dzt(:)
      real(8), allocatable :: dtx(:),dty(:),dtz(:)
      real(8), allocatable :: fxo(:),fyo(:),fzo(:)

c     metadynamics shell thermostat variables

      integer ntshl,keyshl
      real(8) sigma_shl

      logical,save :: lfirst=.true.
      real(8)      :: chit_shl  
      real(8),save :: qmass_shl
      real(8)      :: shlke

c     end metadynamics shell thermostat variables
      
      save newstep,newjob,p0,p1,p2,p3,hstep,qmass,ifre1,ifre2
      save igrp1,igrp2

      data newjob/.true./
      
c     set array allocation error flags

      do i=1,nnn
        fail(i)=0
      enddo

c     assign initial parameters

      if(newjob)then
        
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

        allocate(p0(mxgrp),p1(mxgrp),p2(mxgrp),p3(mxgrp),stat=fail(1))
        newjob=.false.

      endif

c     allocate working arrays

      allocate(gvxo(msgrp),gvyo(msgrp),gvzo(msgrp),stat=fail(2))
      allocate(vxo(msatms),vyo(msatms),vzo(msatms),stat=fail(3))
      allocate(b0(msgrp),b1(msgrp),b2(msgrp),b3(msgrp),
     x  stat=fail(4))
      allocate(c0(msgrp),c1(msgrp),c2(msgrp),c3(msgrp),
     x  stat=fail(5))
      if(isw.eq.1)then

        allocate(xxo(mxatms),yyo(mxatms),zzo(mxatms),stat=fail(6))
        allocate(gxo(msgrp),gyo(msgrp),gzo(msgrp),stat=fail(7))

      endif
      if(ntcons.gt.0)then

        allocate(dxx(mxcons),dyy(mxcons),dzz(mxcons),stat=fail(8))
        allocate(txx(mxatms),tyy(mxatms),tzz(mxatms),stat=fail(9))
        allocate(xxt(mxatms),yyt(mxatms),zzt(mxatms),stat=fail(10))
        allocate(dxt(mxcons),dyt(mxcons),dzt(mxcons),stat=fail(11))

      endif
      allocate(dtx(mxatms),dty(mxatms),dtz(mxatms),stat=fail(12))
      if(isw.eq.2)then
        allocate(fxo(mxatms),fyo(mxatms),fzo(mxatms),stat=fail(13))
      endif

      do i=1,nnn
        if(fail(i).gt.0)call error(idnode,2150)
      enddo

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
      
c     initialise constraint virial

      if(isw.eq.1)then

        vircon=0.d0

        do i=1,9
          strcns(i)=0.d0
        enddo

      endif

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

c     construct current bond vectors
      
      if(ntcons.gt.0)then

        do k=1,nscons
          
          i=listcon(k,2)
          j=listcon(k,3)
          dxx(k)=xxx(i)-xxx(j)
          dyy(k)=yyy(i)-yyy(j)
          dzz(k)=zzz(i)-zzz(j)
          
        enddo

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

c     store key config data at start of time step

      if(isw.eq.1)then

c     atom positions

        do i=1,natms
          
          xxo(i)=xxx(i)
          yyo(i)=yyy(i)
          zzo(i)=zzz(i)
          
        enddo

c     rigid body positions

        j=0
        do i=igrp1,igrp2
          
          j=j+1
          gxo(j)=gcmx(i)
          gyo(j)=gcmy(i)
          gzo(j)=gcmz(i)

        enddo

      endif

c     store free atom atom velocities
      
      j=0
      do ifre=ifre1,ifre2
        
        j=j+1
        i=lstfre(ifre)
        vxo(j)=vxx(i)
        vyo(j)=vyy(i)
        vzo(j)=vzz(i)
        
      enddo

c     store rigid body quaternions, momenta and cartesian velocities

      j=0
      do i=igrp1,igrp2

        j=j+1
        b0(j)=q0(i)
        b1(j)=q1(i)
        b2(j)=q2(i)
        b3(j)=q3(i)
        c0(j)=p0(i)
        c1(j)=p1(i)
        c2(j)=p2(i)
        c3(j)=p3(i)
        gvxo(j)=gvxx(i)
        gvyo(j)=gvyy(i)
        gvzo(j)=gvzz(i)

      enddo

c     store forces if isw = 2

      if(isw.eq.2)then

        do i=1,natms

          fxo(i)=fxx(i)
          fyo(i)=fyy(i)
          fzo(i)=fzz(i)

        enddo
      
      endif

c     -------------- start of shake iteration cycle -------------------
      
      icyc=0
      mxshk=1
      safe=.false.
      newstep=.true.
      if(ntcons.gt.0)mxshk=mxshak

      do while(.not.safe.and.icyc.lt.mxshk)
        
        icyc=icyc+1
        
c     update velocities of free atoms 1/2 timestep

        j=0
        do ifre=ifre1,ifre2

          j=j+1
          i=lstfre(ifre)
          vxx(i)=vxo(j)+(pt5*tstep*rmass(i))*fxx(i)
          vyy(i)=vyo(j)+(pt5*tstep*rmass(i))*fyy(i)
          vzz(i)=vzo(j)+(pt5*tstep*rmass(i))*fzz(i)
          
        enddo

c     *************  rigid body motion ****************************

c     operations common to first and second stages

        jg=0
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
          
          jg=jg+1
          call getrotmat(b0(jg),b1(jg),b2(jg),b3(jg),rot)

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

          tq0=2.0d0*(-b1(jg)*tqx-b2(jg)*tqy-b3(jg)*tqz)
          tq1=2.0d0*( b0(jg)*tqx-b3(jg)*tqy+b2(jg)*tqz)
          tq2=2.0d0*( b3(jg)*tqx+b0(jg)*tqy-b1(jg)*tqz)
          tq3=2.0d0*(-b2(jg)*tqx+b1(jg)*tqy+b0(jg)*tqz)

c     update quaternion momentum by half timestep

          p0(ig)=c0(jg)+tq0*pt5*tstep
          p1(ig)=c1(jg)+tq1*pt5*tstep
          p2(ig)=c2(jg)+tq2*pt5*tstep
          p3(ig)=c3(jg)+tq3*pt5*tstep

c     update centre of mass velocity by half timestep

          gvxx(ig)=gvxo(jg)+fmx*pt5*tstep/gmass(id)
          gvyy(ig)=gvyo(jg)+fmy*pt5*tstep/gmass(id)
          gvzz(ig)=gvzo(jg)+fmz*pt5*tstep/gmass(id)
          
        enddo

c     first stage of velocity verlet algorithm

        if(isw.eq.1)then

          jg=0
          do ig=igrp1,igrp2

            jg=jg+1

c     update centre of mass position by full time step
            
            gcmx(ig)=gxo(jg)+tstep*gvxx(ig)
            gcmy(ig)=gyo(jg)+tstep*gvyy(ig)
            gcmz(ig)=gzo(jg)+tstep*gvzz(ig)
            
c     calculate rotation of rigid groups: nosquish algorithm
            
            q0(ig)=b0(jg)
            q1(ig)=b1(jg)
            q2(ig)=b2(jg)
            q3(ig)=b3(jg)
            call nosquish(ig,tstep,q0,q1,q2,q3,p0,p1,p2,p3)
            
          enddo

c     new atomic positions for atoms in rigid bodies - relative to com
          
          k=0
          do ig=igrp1,igrp2
            
            call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)

            id=lstgtp(ig)
            do j=1,numgsit(id)
              
              k=k+1
              i=lstme(k)
              xxx(i)=rot(1)*gxx(id,j)+rot(2)*gyy(id,j)+
     x          rot(3)*gzz(id,j)+gcmx(ig)
              yyy(i)=rot(4)*gxx(id,j)+rot(5)*gyy(id,j)+
     x          rot(6)*gzz(id,j)+gcmy(ig)
              zzz(i)=rot(7)*gxx(id,j)+rot(8)*gyy(id,j)+
     x          rot(9)*gzz(id,j)+gcmz(ig)
              
            enddo
            
          enddo
          
c     update positions of free particles to full time step
          
          do ifre=ifre1,ifre2
            
            i=lstfre(ifre)
            xxx(i)=xxo(i)+tstep*vxx(i)
            yyy(i)=yyo(i)+tstep*vyy(i)
            zzz(i)=zzo(i)+tstep*vzz(i)
            
          enddo

c     merge atom positions

          if(mxnode.gt.1)call merge1
     x      (idnode,mxnode,natms,lstme,xxx,yyy,zzz,buffer)

c     apply rattle corrections to bond constraints
          
          if(ntcons.gt.0)then
            
            if(mxnode.gt.1)call merge1
     x        (idnode,mxnode,natms,lstme,vxx,vyy,vzz,buffer)

            call qrattle_r
     x        (newstep,safe,lshmov,idnode,imcon,mxnode,natms,
     x        nscons,tolnce,tstep,vircon,dxx,dyy,dzz,dxt,dyt,
     x        dzt,txx,tyy,tzz,xxt,yyt,zzt,strcns)

          endif
          
c     end of first stage 

        else
          
c     second stage of velocity verlet algorithm
          
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

c     new rotational matrix
            
            call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)

c     new atomic velocites
            
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
          
          if(ntcons.gt.0)then
            
            if(mxnode.gt.1)then

              call merge1
     x          (idnode,mxnode,natms,lstme,vxx,vyy,vzz,buffer)
              call merge1
     x          (idnode,mxnode,natms,lstme,fxx,fyy,fzz,buffer)
              
            endif

c     correct constraint bond velocities using rattle

            call qrattle_v
     x        (newstep,safe,lshmov,idnode,mxnode,natms,
     x        nscons,tolnce,tstep,dxx,dyy,dzz,txx,tyy,tzz,
     x        xxt,yyt,zzt)

          endif
          
c     end of second stage
          
        endif
        
        newstep=.false.
        
      enddo

c     check shake convergence
      
      if(.not.safe)call error(idnode,105)

c     sum constraint virial and stress across processors
      
      if(mxnode.gt.1.and.isw.eq.1)then
        
        buffer(1)=vircon
        call gdsum(buffer(1),1,buffer(2))
        vircon=buffer(1)
        call gdsum(strcns,9,buffer)
        
      endif

c     -------------- end of shake iteration cycle -------------------
      
c     apply thermostat for second stage and calculate kinetic energy
      
      if(isw.eq.2)then

c     rigid body contribution to stress tensor

        call bodystress(idnode,mxnode,ngrp,vircom,strbod,dtx,dty,dtz)

c     apply thermostat for second stage and calculate kinetic energy
      
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
        
        engke=engfke+engtrn

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
        
      endif

      if(mxnode.gt.1)then

c     merge new group coordinates and velocities

        if(isw.eq.1)
     x    call merge(idnode,mxnode,ngrp,mxbuff,gcmx,gcmy,gcmz,buffer)
        call merge(idnode,mxnode,ngrp,mxbuff,gvxx,gvyy,gvzz,buffer)
        
c     merge new atomic coordinates and velocities
        
        if(isw.eq.1)
     x    call merge1(idnode,mxnode,natms,lstme,xxx,yyy,zzz,buffer)
        call merge1(idnode,mxnode,natms,lstme,vxx,vyy,vzz,buffer)
        
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
      
c     restore forces if isw = 2

      if(isw.eq.2)then

        do i=1,natms

          fxx(i)=fxo(i)
          fyy(i)=fyo(i)
          fzz(i)=fzo(i)

        enddo
      
      endif

c     deallocate working arrays

      deallocate(gvxo,gvyo,gvzo,vxo,vyo,vzo,stat=fail(1))
      deallocate(b0,b1,b2,b3,c0,c1,c2,c3,stat=fail(2))
      deallocate(dtx,dty,dtz,stat=fail(3))
      if(isw.eq.1)then
        deallocate(xxo,yyo,zzo,gxo,gyo,gzo,stat=fail(4))
      endif
      if(ntcons.gt.0)then

        deallocate(dxx,dyy,dzz,txx,tyy,tzz,stat=fail(5))
        deallocate(dxt,dyt,dzt,xxt,yyt,zzt,stat=fail(6))

      endif
      if(isw.eq.2)then
        deallocate(fxo,fyo,fzo,stat=fail(7))
      endif
      
      return
      end subroutine nvtqvv_h2

      subroutine nptqvv_b2
     x  (safe,lshmov,isw,imcon,idnode,mxnode,natms,ngrp,nscons,
     x  ntcons,ntfree,ntpatm,engke,engrot,press,taut,taup,sigma,
     x  tolnce,tstep,vircom,vircon,elrc,virlrc,virtot,volm)
      
c***********************************************************************
c     
c     dlpoly subroutine to integrate rigid body equations of motion
c     using the symplectic no_squish quaternion algorithm of 
c     miller et al j.chem.phys 116 (2002) 8649
c     plus bond constraints -  including rigid body sites linked
c     by constraint sites (qrattle algorithm)
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
c     author      w.smith   sep 2005
c     
c**********************************************************************

      implicit none

      logical newstep,safe,lshmov,newjob

      integer isw,imcon,idnode,mxnode,natms,ngrp,nscons,ntcons
      integer ntfree,i,j,k,igrp1,igrp2,igrp,ifre1,ifre2,jg,jr
      integer id,ifre,jrs,icyc,mxshk,idum,ig,ntpatm
      integer iter,mxiter

      real(8) engke,engrot,tolnce,tstep,vircom,vircon,engtke,engtrn
      real(8) vaa,vbb,vcc,opx,opy,opz,ftx,fty,ftz,volm0
      real(8) fmx,fmy,fmz,tqx,tqy,tqz,tq0,tq1,tq2,tq3,taut,sigma,chit
      real(8) volm,elrc0,elrc,virlrc0,virlrc,scale,psyst,virtot,chip
      real(8) beta,press,taup,engfke,vzero

      integer, parameter :: nnn=13
      real(8), parameter :: pt5=0.5d0

      integer fail(nnn)
      real(8) rot(9),strkin(9),strgrp(9),cell0(9),uni(9)

      real(8), allocatable :: dens0(:)
      real(8), allocatable :: xxo(:),yyo(:),zzo(:)
      real(8), allocatable :: vxo(:),vyo(:),vzo(:)
      real(8), allocatable :: gxo(:),gyo(:),gzo(:)
      real(8), allocatable :: dxx(:),dyy(:),dzz(:)
      real(8), allocatable :: txx(:),tyy(:),tzz(:)
      real(8), allocatable :: xxt(:),yyt(:),zzt(:)
      real(8), allocatable :: dxt(:),dyt(:),dzt(:)
      real(8), allocatable :: dtx(:),dty(:),dtz(:)
      real(8), allocatable :: fxo(:),fyo(:),fzo(:)
      real(8), allocatable :: b0(:),b1(:),b2(:),b3(:)
      real(8), allocatable :: c0(:),c1(:),c2(:),c3(:)
      real(8), allocatable :: p0(:),p1(:),p2(:),p3(:)
      real(8), allocatable :: gvxo(:),gvyo(:),gvzo(:)

      save newstep,newjob,volm0,elrc0,virlrc0,dens0
      save p0,p1,p2,p3,ifre1,ifre2,igrp1,igrp2

      data uni/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/
      data newjob/.true./,beta/7.3728d-3/
      
c     set array alocation error flags
      
      do i=1,nnn
        fail(i)=0
      enddo
      
c     assign initial parameters

      if(newjob)then

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
        
        allocate(dens0(mxatyp),stat=fail(1))
        if(fail(1).gt.0)call error(idnode,2180)
        
c     store initial values of volume and long range corrections
        
        volm0=volm
        elrc0=elrc
        virlrc0=virlrc
        do i=1,ntpatm
          dens0(i)=dens(i)
        enddo
        
        allocate(p0(mxgrp),p1(mxgrp),p2(mxgrp),p3(mxgrp),stat=fail(1))
        newjob=.false.

      endif
      
c     allocate working arrays
      
      if(isw.eq.1)then

        allocate(xxo(mxatms),yyo(mxatms),zzo(mxatms),stat=fail(2))
        allocate(gxo(msgrp),gyo(msgrp),gzo(msgrp),stat=fail(3))

      endif
      if(ntcons.gt.0)then

        allocate(dxx(mxcons),dyy(mxcons),dzz(mxcons),stat=fail(4))
        allocate(txx(mxatms),tyy(mxatms),tzz(mxatms),stat=fail(5))
        allocate(xxt(mxatms),yyt(mxatms),zzt(mxatms),stat=fail(6))
        allocate(dxt(mxcons),dyt(mxcons),dzt(mxcons),stat=fail(7))

      endif
      allocate(fxo(mxatms),fyo(mxatms),fzo(mxatms),stat=fail(8))
      allocate(dtx(mxatms),dty(mxatms),dtz(mxatms),stat=fail(9))
      allocate(vxo(msatms),vyo(msatms),vzo(msatms),stat=fail(10))
      allocate(gvxo(msgrp),gvyo(msgrp),gvzo(msgrp),stat=fail(11))
      allocate(b0(msgrp),b1(msgrp),b2(msgrp),b3(msgrp),
     x  stat=fail(12))
      allocate(c0(msgrp),c1(msgrp),c2(msgrp),c3(msgrp),
     x  stat=fail(13))

c     check array allocation error flags

      do i=1,nnn
        if(fail(i).gt.0)call error(idnode,2190)
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
      
c     store key config data at start of time step
      
      if(isw.eq.1)then
        
c     cell parameters
        
        vzero=volm
        do i=1,9
          cell0(i)=cell(i)
        enddo
        
c     atom positions

        do i=1,natms
          
          xxo(i)=xxx(i)
          yyo(i)=yyy(i)
          zzo(i)=zzz(i)
          
        enddo

c     rigid body positions

        j=0
        do i=igrp1,igrp2
          
          j=j+1
          gxo(j)=gcmx(i)
          gyo(j)=gcmy(i)
          gzo(j)=gcmz(i)

        enddo
        
      endif

c     store free atom velocities
      
      j=0
      do ifre=ifre1,ifre2
        
        j=j+1
        i=lstfre(ifre)
        vxo(j)=vxx(i)
        vyo(j)=vyy(i)
        vzo(j)=vzz(i)
        
      enddo

c     store rigid body quaternions, momenta and cartesian velocities

      j=0
      do i=igrp1,igrp2

        j=j+1
        b0(j)=q0(i)
        b1(j)=q1(i)
        b2(j)=q2(i)
        b3(j)=q3(i)
        c0(j)=p0(i)
        c1(j)=p1(i)
        c2(j)=p2(i)
        c3(j)=p3(i)
        gvxo(j)=gvxx(i)
        gvyo(j)=gvyy(i)
        gvzo(j)=gvzz(i)

      enddo

c     store forces

      do i=1,natms
        
        fxo(i)=fxx(i)
        fyo(i)=fyy(i)
        fzo(i)=fzz(i)
        
      enddo

      if(isw.eq.1)then
        
c     calculate kinetic energy
        
        engfke=getkinf(ntfree,idnode,mxnode)
        call getking(ngrp,idnode,mxnode,engtrn,engrot)
        engke=engfke+engtrn
        
      endif
      
c     -------------- start of barostat iteration cycle -----------------
      
      mxiter=1
      if(isw.eq.1.and.ntcons.gt.0)mxiter=3
      do iter=1,mxiter
        
        if(isw.eq.1)then
          
c     restore cell parameters
          
          volm=vzero
          do i=1,9
            cell(i)=cell0(i)
          enddo
          
c     calculate system pressure
          
          vircon=-(strcns(1)+strcns(5)+strcns(9))
          psyst=(2.d0*engke-virtot-vircon-vircom)/(3.d0*volm)
          
c     apply Berendsen barostat
          
          chip=1.d0+beta*tstep*(psyst-press)/taup
          chip=1.d0
          scale=chip**(1.d0/3.d0)
          volm=chip*volm
          
c     reset cell parameters for new volume
          
          do i=1,9
            cell(i)=scale*cell(i)
          enddo
          
c     reset constraint virial
          
          vircon=0.d0
          do i=1,9
            strcns(i)=0.d0
          enddo
          
        endif

c     -------------- start of shake iteration cycle -------------------
        
        icyc=0
        mxshk=1
        safe=.false.
        newstep=.true.
        if(ntcons.gt.0)mxshk=mxshak
        do while(.not.safe.and.icyc.lt.mxshk)
          
          icyc=icyc+1
          
c     update velocities of free atoms 1/2 timestep
          
          j=0
          do ifre=ifre1,ifre2
            
            j=j+1
            i=lstfre(ifre)
            vxx(i)=vxo(j)+(pt5*tstep*rmass(i))*fxx(i)
            vyy(i)=vyo(j)+(pt5*tstep*rmass(i))*fyy(i)
            vzz(i)=vzo(j)+(pt5*tstep*rmass(i))*fzz(i)
            
          enddo

c     *************  rigid body motion ****************************

c     operations common to first and second stages

          jg=0
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
            
            jg=jg+1
            call getrotmat(b0(jg),b1(jg),b2(jg),b3(jg),rot)

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

            tq0=2.0d0*(-b1(jg)*tqx-b2(jg)*tqy-b3(jg)*tqz)
            tq1=2.0d0*( b0(jg)*tqx-b3(jg)*tqy+b2(jg)*tqz)
            tq2=2.0d0*( b3(jg)*tqx+b0(jg)*tqy-b1(jg)*tqz)
            tq3=2.0d0*(-b2(jg)*tqx+b1(jg)*tqy+b0(jg)*tqz)

c     update quaternion momentum by half timestep

            p0(ig)=c0(jg)+tq0*pt5*tstep
            p1(ig)=c1(jg)+tq1*pt5*tstep
            p2(ig)=c2(jg)+tq2*pt5*tstep
            p3(ig)=c3(jg)+tq3*pt5*tstep

c     update centre of mass velocity by half timestep

            gvxx(ig)=gvxo(jg)+fmx*pt5*tstep/gmass(id)
            gvyy(ig)=gvyo(jg)+fmy*pt5*tstep/gmass(id)
            gvzz(ig)=gvzo(jg)+fmz*pt5*tstep/gmass(id)
            
          enddo
          
c     first stage of velocity verlet algorithm

          if(isw.eq.1)then

            jg=0
            do ig=igrp1,igrp2

              jg=jg+1

c     update centre of mass position by full time step
              
              gcmx(ig)=scale*gxo(jg)+tstep*gvxx(ig)
              gcmy(ig)=scale*gyo(jg)+tstep*gvyy(ig)
              gcmz(ig)=scale*gzo(jg)+tstep*gvzz(ig)
              
c     calculate rotation of rigid groups: nosquish algorithm
              
              q0(ig)=b0(jg)
              q1(ig)=b1(jg)
              q2(ig)=b2(jg)
              q3(ig)=b3(jg)
              call nosquish(ig,tstep,q0,q1,q2,q3,p0,p1,p2,p3)
              
            enddo

c     new atomic positions for atoms in rigid bodies - relative to com
            
            k=0
            do ig=igrp1,igrp2
              
              call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)

              id=lstgtp(ig)
              do j=1,numgsit(id)
                
                k=k+1
                i=lstme(k)
                xxx(i)=rot(1)*gxx(id,j)+rot(2)*gyy(id,j)+
     x            rot(3)*gzz(id,j)+gcmx(ig)
                yyy(i)=rot(4)*gxx(id,j)+rot(5)*gyy(id,j)+
     x            rot(6)*gzz(id,j)+gcmy(ig)
                zzz(i)=rot(7)*gxx(id,j)+rot(8)*gyy(id,j)+
     x            rot(9)*gzz(id,j)+gcmz(ig)
                
              enddo
              
            enddo

c     update positions of free particles to full time step
            
            do ifre=ifre1,ifre2
              
              i=lstfre(ifre)
              xxx(i)=scale*xxo(i)+tstep*vxx(i)
              yyy(i)=scale*yyo(i)+tstep*vyy(i)
              zzz(i)=scale*zzo(i)+tstep*vzz(i)
              
            enddo

c     merge atom positions
            
            if(mxnode.gt.1)call merge1
     x        (idnode,mxnode,natms,lstme,xxx,yyy,zzz,buffer)

c     apply rattle corrections to bond constraints
            
            if(ntcons.gt.0)then
              
              if(mxnode.gt.1)call merge1
     x          (idnode,mxnode,natms,lstme,vxx,vyy,vzz,buffer)
              
              call qrattle_r
     x          (newstep,safe,lshmov,idnode,imcon,mxnode,natms,
     x          nscons,tolnce,tstep,vircon,dxx,dyy,dzz,dxt,dyt,
     x          dzt,txx,tyy,tzz,xxt,yyt,zzt,strcns)
              
            endif
            
c     end of first stage 
            
          endif

c     second stage of velocity verlet algorithm
          
          if(isw.eq.2)then

            jr=0
            do ig=igrp1,igrp2
              
c     new angular momenta and velocities

              opx=pt5*(-q1(ig)*p0(ig)+q0(ig)*p1(ig)+
     x          q3(ig)*p2(ig)-q2(ig)*p3(ig))
              opy=pt5*(-q2(ig)*p0(ig)-q3(ig)*p1(ig)+
     x          q0(ig)*p2(ig)+q1(ig)*p3(ig))
              opz=pt5*(-q3(ig)*p0(ig)+q2(ig)*p1(ig)-
     x          q1(ig)*p2(ig)+q0(ig)*p3(ig))

              id=lstgtp(ig)

              omx(ig)=opx*rotinx(id,2)
              omy(ig)=opy*rotiny(id,2)
              omz(ig)=opz*rotinz(id,2)

c     new rotational matrix
              
              call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)

c     new atomic velocites
              
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
            
            if(ntcons.gt.0)then
              
              if(mxnode.gt.1)then
                
                call merge1
     x            (idnode,mxnode,natms,lstme,vxx,vyy,vzz,buffer)
                call merge1
     x            (idnode,mxnode,natms,lstme,fxx,fyy,fzz,buffer)
                
              endif

c     correct constraint bond velocities using rattle

              call qrattle_v
     x          (newstep,safe,lshmov,idnode,mxnode,natms,
     x          nscons,tolnce,tstep,dxx,dyy,dzz,txx,tyy,tzz,
     x          xxt,yyt,zzt)

            endif
            
c     end of second stage
            
          endif
          
          newstep=.false.
          
        enddo

c     check shake convergence
        
        if(.not.safe)call error(idnode,105)
        
c     sum constraint virial and stress across processors
      
        if(mxnode.gt.1.and.isw.eq.1)then
          
          buffer(1)=vircon
          call gdsum(buffer(1),1,buffer(2))
          vircon=buffer(1)
          call gdsum(strcns,9,buffer)
          
        endif

c     -------------- end of shake iteration cycle -------------------

c     rigid body contribution to stress tensor

        if(isw.eq.2)call bodystress
     x    (idnode,mxnode,ngrp,vircom,strbod,dtx,dty,dtz)

c     restore  forces
          
        do i=1,natms
          
          fxx(i)=fxo(i)
          fyy(i)=fyo(i)
          fzz(i)=fzo(i)
          
        enddo

c     -------------- end of barostat iteration cycle ----------------
        
      enddo
      
      if(isw.eq.1)then

c     adjust long range corrections and number density
        
        elrc=elrc0*(volm0/volm)
        virlrc=virlrc0*(volm0/volm)
        
        do k=1,ntpatm
          dens(k)=dens0(k)*(volm0/volm)
        enddo
        
c     construct scaling tensor for tethered bonds
        
        do i=1,9
          eta(i)=uni(i)*scale
        enddo
        
      endif
      
      if(isw.eq.2)then
        
c     calculate kinetic energy
      
        engfke=getkinf(ntfree,idnode,mxnode)      
        call getking(ngrp,idnode,mxnode,engtrn,engrot)
        engtke=engfke+engtrn+engrot
        engke=engfke+engtrn
        
c     apply Berendsen thermostat - taut is the relaxation time
        
        chit=sqrt(1.d0+tstep/taut*(sigma/engtke-1.d0))
        
        engke=engke*chit**2
        engtrn=engtrn*chit**2
        engrot=engrot*chit**2

c     thermostat velocities of free particles
        
        do ifre=ifre1,ifre2
          
          i=lstfre(ifre)
          if(lstfrz(i).ne.0)then
            
            vxx(i)=chit*vxx(i)
            vyy(i)=chit*vyy(i)
            vzz(i)=chit*vzz(i)
            
          endif
          
        enddo
        
c     thermostat rigid body velocities
        
        do ig=igrp1,igrp2
          
          id=lstgtp(ig)
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
      
      if(mxnode.gt.1)then

c     merge new group coordinates and velocities

        if(isw.eq.1)
     x    call merge(idnode,mxnode,ngrp,mxbuff,gcmx,gcmy,gcmz,buffer)
        call merge(idnode,mxnode,ngrp,mxbuff,gvxx,gvyy,gvzz,buffer)
        
c     merge new atomic coordinates and velocities
        
        if(isw.eq.1)
     x    call merge1(idnode,mxnode,natms,lstme,xxx,yyy,zzz,buffer)
        call merge1(idnode,mxnode,natms,lstme,vxx,vyy,vzz,buffer)
        
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

      deallocate(gvxo,gvyo,gvzo,vxo,vyo,vzo,stat=fail(1))
      deallocate(b0,b1,b2,b3,c0,c1,c2,c3,stat=fail(2))
      deallocate(dtx,dty,dtz,fxo,fyo,fzo,stat=fail(3))
      if(isw.eq.1)then
        deallocate(xxo,yyo,zzo,gxo,gyo,gzo,stat=fail(4))
      endif
      if(ntcons.gt.0)then

        deallocate(dxx,dyy,dzz,txx,tyy,tzz,stat=fail(5))
        deallocate(dxt,dyt,dzt,xxt,yyt,zzt,stat=fail(6))

      endif
      
      return
      end subroutine nptqvv_b2

      subroutine nptqvv_h2
     x  (safe,lshmov,isw,idnode,mxnode,natms,imcon,ngrp,nscons,
     x  ntcons,ntpatm,ntfree,ntshl,keyshl,tstep,taut,taup,sigma,
     x  temp,chip,chit,consv,conint,engke,engrot,elrc,tolnce,
     x  vircom,vircon,virtot,virlrc,volm,press,chit_shl,sigma_shl)
      
c***********************************************************************
c     
c     dlpoly subroutine to integrate rigid body equations of motion
c     using the symplectic no_squish quaternion algorithm of 
c     miller et al j.chem.phys 116 (2002) 8649
c     plus bond constraints -  including rigid body sites linked
c     by constraint sites (qrattle algorithm)
c     
c     npt ensemble - nose-hoover thermostat Molec Phys 87 (1996) 1117
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

      logical newstep,safe,lshmov,newjob

      integer isw,imcon,idnode,mxnode,natms,ngrp,nscons,ntcons
      integer ntfree,i,j,k,igrp1,igrp2,igrp,ifre1,ifre2,jr
      integer id,ifre,icyc,mxshk,idum,ig,ntpatm
      integer jcyc,iter,mxiter

      real(8) engke,engrot,tolnce,tstep,vircom,vircon
      real(8) engtrn,vaa,vbb,vcc,opx,opy,opz,engfke
      real(8) taut,taup,sigma,chit,hstep,qmass,conint,consv
      real(8) cxx,cyy,czz,scale,virtot,press,chip,temp
      real(8) volm,pmass,totmas,qstep,fstep,volm0,elrc
      real(8) virlrc,elrc0,virlrc0,chit0,chip0,vzero,cons0

      integer, parameter :: nnn=16
      integer, parameter :: ncyc=5
      real(8), parameter :: pt5=0.5d0

      integer fail(nnn)
      real(8) rot(9),strkin(9),strgrp(9),cell0(9),com(3),vom(3),uni(9)

      real(8), allocatable :: dens0(:)
      real(8), allocatable :: xxo(:),yyo(:),zzo(:)
      real(8), allocatable :: vxo(:),vyo(:),vzo(:)
      real(8), allocatable :: vx1(:),vy1(:),vz1(:)
      real(8), allocatable :: fxo(:),fyo(:),fzo(:)
      real(8), allocatable :: gxo(:),gyo(:),gzo(:)
      real(8), allocatable :: oxo(:),oyo(:),ozo(:)
      real(8), allocatable :: dxx(:),dyy(:),dzz(:)
      real(8), allocatable :: txx(:),tyy(:),tzz(:)
      real(8), allocatable :: xxt(:),yyt(:),zzt(:)
      real(8), allocatable :: dxt(:),dyt(:),dzt(:)
      real(8), allocatable :: dtx(:),dty(:),dtz(:)
      real(8), allocatable :: b0(:),b1(:),b2(:),b3(:)
      real(8), allocatable :: c0(:),c1(:),c2(:),c3(:)
      real(8), allocatable :: p0(:),p1(:),p2(:),p3(:)
      real(8), allocatable :: gvxo(:),gvyo(:),gvzo(:)
      real(8), allocatable :: gvx1(:),gvy1(:),gvz1(:)

c     metadynamics shell thermostat variables

      integer ntshl,keyshl
      real(8) sigma_shl

      logical,save :: lfirst=.true.
      real(8)      :: chit_shl  
      real(8),save :: qmass_shl
      real(8)      :: shlke

c     end metadynamics shell thermostat variables
      
      save newstep,newjob,p0,p1,p2,p3,hstep,fstep,qmass,ifre1,ifre2
      save igrp1,igrp2,volm0,elrc0,virlrc0,qstep,dens0,totmas
      save pmass

      data uni/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/
      data newjob/.true./
      
c     set array allocation error flags

      do i=1,nnn
        fail(i)=0
      enddo

c     assign initial parameters

      if(newjob)then
        
c     store intitial parameters
        
        volm0=volm
        elrc0=elrc
        virlrc0=virlrc
        hstep=0.5d0*tstep
        fstep=0.5d0*tstep/dble(ncyc)
        qstep=0.25d0*tstep/dble(ncyc)
        
c     allocate density storage array

        fail(1)=0
        allocate(dens0(mxatyp),stat=fail(1))
        if(fail(1).gt.0)call error(idnode,2220)

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
          call  error(idnode,506)
        endif

        allocate(p0(mxgrp),p1(mxgrp),p2(mxgrp),p3(mxgrp),stat=fail(1))
        newjob=.false.

      endif

      if(ntcons.gt.0)safe=.false.

c     allocate working arrays

      if(isw.eq.1)then

        allocate(xxo(mxatms),yyo(mxatms),zzo(mxatms),stat=fail(2))
        allocate(gxo(mxgrp),gyo(mxgrp),gzo(mxgrp),stat=fail(3))

      endif
      if(ntcons.gt.0)then
        
        allocate(dxx(mxcons),dyy(mxcons),dzz(mxcons),stat=fail(4))
        allocate(txx(mxatms),tyy(mxatms),tzz(mxatms),stat=fail(5))
        allocate(xxt(mxatms),yyt(mxatms),zzt(mxatms),stat=fail(6))
        allocate(dxt(mxcons),dyt(mxcons),dzt(mxcons),stat=fail(7))
        
      endif
      allocate(dtx(mxatms),dty(mxatms),dtz(mxatms),stat=fail(8))
      allocate(vxo(mxatms),vyo(mxatms),vzo(mxatms),stat=fail(9))
      allocate(vx1(mxatms),vy1(mxatms),vz1(mxatms),stat=fail(10))
      allocate(fxo(mxatms),fyo(mxatms),fzo(mxatms),stat=fail(11))
      allocate(oxo(mxatms),oyo(mxatms),ozo(mxatms),stat=fail(12))
      allocate(b0(mxgrp),b1(mxgrp),b2(mxgrp),b3(mxgrp),
     x  stat=fail(13))
      allocate(c0(mxgrp),c1(mxgrp),c2(mxgrp),c3(mxgrp),
     x  stat=fail(14))
      allocate(gvxo(mxgrp),gvyo(mxgrp),gvzo(mxgrp),stat=fail(15))
      allocate(gvx1(mxgrp),gvy1(mxgrp),gvz1(mxgrp),stat=fail(16))

      do i=1,nnn
        if(fail(i).gt.0)call error(idnode,2230)
      enddo
      
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
      
c     construct current bond vectors
      
      if(ntcons.gt.0)then

        do k=1,nscons
          
          i=listcon(k,2)
          j=listcon(k,3)
          dxx(k)=xxx(i)-xxx(j)
          dyy(k)=yyy(i)-yyy(j)
          dzz(k)=zzz(i)-zzz(j)
          
        enddo

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
      
c     store key config data at start of timestep

      vzero=volm
      chit0=chit
      chip0=chip
      cons0=conint
      do i=1,9
        cell0(i)=cell(i)
      enddo

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

c     store atom positions
        
        do i=1,natms
          
          xxo(i)=xxx(i)
          yyo(i)=yyy(i)
          zzo(i)=zzz(i)
          
        enddo
        
c     store rigid body positions
        
        do ig=1,ngrp
          
          gxo(ig)=gcmx(ig)
          gyo(ig)=gcmy(ig)
          gzo(ig)=gcmz(ig)
          
        enddo
        
      endif
      
c     store free atom velocities

      do i=1,natms
        
        vxo(i)=vxx(i)
        vyo(i)=vyy(i)
        vzo(i)=vzz(i)
        
      enddo

c     store forces

      do i=1,natms
        
        fxo(i)=fxx(i)
        fyo(i)=fyy(i)
        fzo(i)=fzz(i)
        
      enddo
      
c     store rigid body quaternions, angular and cartesian velocities

      do ig=1,ngrp
        
        b0(ig)=q0(ig)
        b1(ig)=q1(ig)
        b2(ig)=q2(ig)
        b3(ig)=q3(ig)
        oxo(ig)=omx(ig)
        oyo(ig)=omy(ig)
        ozo(ig)=omz(ig)
        gvxo(ig)=gvxx(ig)
        gvyo(ig)=gvyy(ig)
        gvzo(ig)=gvzz(ig)

      enddo
      
c     iteration necessary if ntcons > 0 and isw=1

      mxiter=1
      if(isw.eq.1.and.ntcons.gt.0)mxiter=3
      do iter=1,mxiter

c     integration of barostat and thermostat (part 1)
        
        if(isw.eq.1)then
          
c     restore cell parameters

          volm=vzero
          chit=chit0
          chip=chip0
          conint=cons0
          do i=1,9
            cell(i)=cell0(i)
          enddo

c     restore free atom velocities
          
          do i=1,natms
            
            vxx(i)=vxo(i)
            vyy(i)=vyo(i)
            vzz(i)=vzo(i)
            
          enddo

c     restore rigid body quaternions angular and cartesian velocities

          do ig=1,ngrp
            
            omx(ig)=oxo(ig)
            omy(ig)=oyo(ig)
            omz(ig)=ozo(ig)
            gvxx(ig)=gvxo(ig)
            gvyy(ig)=gvyo(ig)
            gvzz(ig)=gvzo(ig)
            
          enddo
          
c     current constraint virial

          vircon=-(strcns(1)+strcns(5)+strcns(9))

          do jcyc=1,ncyc
            
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
          
c     translational kinetic energy

          engke=engfke+engtrn
          
c     scale cell vectors - isotropic
          
          scale=(volm/vzero)**(1.d0/3.d0)
          do i=1,9
            cell(i)=cell0(i)*scale
          enddo
          
c     reset constraint virial and stress
          
          vircon=0.d0
          do i=1,9
            strcns(i)=0.d0
          enddo
          
c     calculate quaternion momenta
          
          do ig=igrp1,igrp2
            
            id=lstgtp(ig)
            opx=omx(ig)*rotinx(id,1)
            opy=omy(ig)*rotiny(id,1)
            opz=omz(ig)*rotinz(id,1)
            p0(ig)=2.0d0*(-b1(ig)*opx-b2(ig)*opy-b3(ig)*opz)
            p1(ig)=2.0d0*( b0(ig)*opx-b3(ig)*opy+b2(ig)*opz)
            p2(ig)=2.0d0*( b3(ig)*opx+b0(ig)*opy-b1(ig)*opz)
            p3(ig)=2.0d0*(-b2(ig)*opx+b1(ig)*opy+b0(ig)*opz)
            
          enddo
          
        endif
        
c     store intermediate velocities

        do ifre=ifre1,ifre2
          
          i=lstfre(ifre)
          vx1(i)=vxx(i)
          vy1(i)=vyy(i)
          vz1(i)=vzz(i)
          
        enddo
        do ig=igrp1,igrp2
          
          c0(ig)=p0(ig)
          c1(ig)=p1(ig)
          c2(ig)=p2(ig)
          c3(ig)=p3(ig)
          gvx1(ig)=gvxx(ig)
          gvy1(ig)=gvyy(ig)
          gvz1(ig)=gvzz(ig)
          
        enddo

c     -------------- start of shake iteration cycle -------------------
        
        icyc=0
        mxshk=1
        safe=.false.
        newstep=.true.
        if(ntcons.gt.0)mxshk=mxshak
        do while(.not.safe.and.icyc.lt.mxshk)
          
          icyc=icyc+1
          
c     update velocities of free atoms 1/2 timestep
          
          do ifre=ifre1,ifre2
            
            i=lstfre(ifre)
            vxx(i)=vx1(i)+hstep*rmass(i)*fxx(i)
            vyy(i)=vy1(i)+hstep*rmass(i)*fyy(i)
            vzz(i)=vz1(i)+hstep*rmass(i)*fzz(i)
            
          enddo
          
c     *************  rigid body motion ****************************
          
c     restore rigid body quaternions, angular momenta and velocities

          do ig=igrp1,igrp2
            
            q0(ig)=b0(ig)
            q1(ig)=b1(ig)
            q2(ig)=b2(ig)
            q3(ig)=b3(ig)
            p0(ig)=c0(ig)
            p1(ig)=c1(ig)
            p2(ig)=c2(ig)
            p3(ig)=c3(ig)
            gvxx(ig)=gvx1(ig)
            gvyy(ig)=gvy1(ig)
            gvzz(ig)=gvz1(ig)
            
          enddo

c     calculate new rigid body velocities

          call rotate_omega
     x      (idnode,mxnode,ngrp,hstep,p0,p1,p2,p3,dtx,dty,dtz)

c     first stage of velocity verlet algorithm

          if(isw.eq.1)then
            
c     calculate system centre of mass
            
            call getcom(natms,idnode,mxnode,totmas,com)
            
c     update centre of mass position by full time step
      
            do ig=igrp1,igrp2
              
              cxx=gxo(ig)-com(1)
              cyy=gyo(ig)-com(2)
              czz=gzo(ig)-com(3)
              gcmx(ig)=gxo(ig)+tstep*(gvxx(ig)+chip*cxx)
              gcmy(ig)=gyo(ig)+tstep*(gvyy(ig)+chip*cyy)
              gcmz(ig)=gzo(ig)+tstep*(gvzz(ig)+chip*czz)
              
c     calculate rotation of rigid groups: nosquish algorithm
              
              call nosquish(ig,tstep,q0,q1,q2,q3,p0,p1,p2,p3)
              
            enddo
            
c     merge group coms from all nodes
            
            if(mxnode.gt.1)call merge
     x        (idnode,mxnode,ngrp,mxbuff,gcmx,gcmy,gcmz,buffer)
            
c     new atomic positions for atoms in rigid bodies - relative to com
            
            k=0
            do ig=igrp1,igrp2
              
              call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)

              id=lstgtp(ig)
              do j=1,numgsit(id)
                
                k=k+1
                i=lstme(k)
                xxx(i)=rot(1)*gxx(id,j)+rot(2)*gyy(id,j)+
     x            rot(3)*gzz(id,j)+gcmx(ig)
                yyy(i)=rot(4)*gxx(id,j)+rot(5)*gyy(id,j)+
     x            rot(6)*gzz(id,j)+gcmy(ig)
                zzz(i)=rot(7)*gxx(id,j)+rot(8)*gyy(id,j)+
     x            rot(9)*gzz(id,j)+gcmz(ig)
                
              enddo
              
            enddo
            
c     update positions of free particles to full time step
            
            do ifre=ifre1,ifre2
              
              i=lstfre(ifre)
              cxx=xxo(i)-com(1)
              cyy=yyo(i)-com(2)
              czz=zzo(i)-com(3)
              xxx(i)=xxo(i)+tstep*(vxx(i)+chip*cxx)
              yyy(i)=yyo(i)+tstep*(vyy(i)+chip*cyy)
              zzz(i)=zzo(i)+tstep*(vzz(i)+chip*czz)
              
            enddo

c     merge atom positions

            if(mxnode.gt.1)call merge1
     x        (idnode,mxnode,natms,lstme,xxx,yyy,zzz,buffer)

c     apply rattle corrections to bond constraints
            
            if(ntcons.gt.0)then
              
              if(mxnode.gt.1)call merge1
     x          (idnode,mxnode,natms,lstme,vxx,vyy,vzz,buffer)
              
              call qrattle_r
     x          (newstep,safe,lshmov,idnode,imcon,mxnode,natms,
     x          nscons,tolnce,tstep,vircon,dxx,dyy,dzz,dxt,dyt,
     x          dzt,txx,tyy,tzz,xxt,yyt,zzt,strcns)
              
            endif

c     end of first stage 
            
          endif
          
c     second stage of velocity verlet algorithm
          
          if(isw.eq.2)then
            
            jr=0
            do ig=igrp1,igrp2
              
c     new angular momenta and velocities

              opx=pt5*(-q1(ig)*p0(ig)+q0(ig)*p1(ig)+
     x          q3(ig)*p2(ig)-q2(ig)*p3(ig))
              opy=pt5*(-q2(ig)*p0(ig)-q3(ig)*p1(ig)+
     x          q0(ig)*p2(ig)+q1(ig)*p3(ig))
              opz=pt5*(-q3(ig)*p0(ig)+q2(ig)*p1(ig)-
     x          q1(ig)*p2(ig)+q0(ig)*p3(ig))

              id=lstgtp(ig)

              omx(ig)=opx*rotinx(id,2)
              omy(ig)=opy*rotiny(id,2)
              omz(ig)=opz*rotinz(id,2)
              
c     new rotational matrix
              
              call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)

c     new atomic velocites
              
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
            
            if(ntcons.gt.0)then
              
              if(mxnode.gt.1)then

                call merge1
     x            (idnode,mxnode,natms,lstme,vxx,vyy,vzz,buffer)
                call merge1
     x            (idnode,mxnode,natms,lstme,fxx,fyy,fzz,buffer)
                
              endif

c     correct constraint bond velocities using rattle

              call qrattle_v
     x          (newstep,safe,lshmov,idnode,mxnode,natms,
     x          nscons,tolnce,tstep,dxx,dyy,dzz,txx,tyy,tzz,
     x          xxt,yyt,zzt)

            endif
            
c     end of second stage
            
          endif
          
          newstep=.false.
          
        enddo

c     check shake convergence
        
        if(.not.safe)call error(idnode,105)

c     sum constraint virial and stress across processors
      
        if(mxnode.gt.1.and.isw.eq.1)then
          
          buffer(1)=vircon
          call gdsum(buffer(1),1,buffer(2))
          vircon=buffer(1)
          call gdsum(strcns,9,buffer)
          
        endif

c     -------------- end of shake iteration cycle -------------------
        
c     rigid body contribution to stress tensor
        
        call bodystress(idnode,mxnode,ngrp,vircom,strbod,dtx,dty,dtz)
        
c     integration of barostat and thermostat (part 2)
        
        if(isw.eq.2)then

c     current constraint virial

          vircon=-(strcns(1)+strcns(5)+strcns(9))

          do jcyc=1,ncyc
            
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
          
c     translational kinetic energy

          engke=engfke+engtrn

c     scale cell vectors - isotropic
          
          scale=(volm/vzero)**(1.d0/3.d0)
          do i=1,9
            cell(i)=cell0(i)*scale
          enddo
          
        endif

c     restore  forces
        
        do i=1,natms
          
          fxx(i)=fxo(i)
          fyy(i)=fyo(i)
          fzz(i)=fzo(i)
          
        enddo

c     -------------- end of barostat iteration cycle ----------------

      enddo

      if(isw.eq.2)then

c     calculate conserved variable

        consv=conint+0.5d0*qmass*chit**2+press*volm
     x    +0.5d0*pmass*chip**2
        
c     metadynamics shell thermostat

        if(lmetadyn.and.keyshl.eq.1)then
          consv=consv+0.5d0*qmass_shl*chit_shl**2
        endif
        
c     merge velocity arrays
        
        if(mxnode.gt.1)call merge1
     x    (idnode,mxnode,natms,lstme,vxx,vyy,vzz,buffer)
          
c     kinetic contribution to stress tensor
        
        call kinstressf(ntfree,idnode,mxnode,strkin)        
        call kinstressg(ngrp,idnode,mxnode,strgrp)
        
c     add contributions to stress tensor
        
        do i=1,9
          stress(i)=stress(i)+strbod(i)+strcns(i)+strkin(i)+strgrp(i)
        enddo
        
      endif

c     adjust long range corrections and number density
      
      elrc=elrc0*(volm0/volm)
      virlrc=virlrc0*(volm0/volm)
      do k=1,ntpatm
        dens(k)=dens0(k)*(volm0/volm)
      enddo
      
c     construct scaling tensor (for tethered atoms)
      
      do i=1,9
        eta(i)=chip*uni(i)
      enddo
      
      if(mxnode.gt.1)then

c     merge new group coordinates and velocities

        if(isw.eq.1)
     x    call merge(idnode,mxnode,ngrp,mxbuff,gcmx,gcmy,gcmz,buffer)
        call merge(idnode,mxnode,ngrp,mxbuff,gvxx,gvyy,gvzz,buffer)
        
c     merge new atomic coordinates and velocities
        
        if(isw.eq.1)
     x    call merge1(idnode,mxnode,natms,lstme,xxx,yyy,zzz,buffer)
        call merge1(idnode,mxnode,natms,lstme,vxx,vyy,vzz,buffer)
        
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
      
c     deallocate working arrays

      deallocate(gvxo,gvyo,gvzo,vxo,vyo,vzo,stat=fail(1))
      deallocate(b0,b1,b2,b3,c0,c1,c2,c3,stat=fail(2))
      deallocate(oxo,oyo,ozo,dtx,dty,dtz,stat=fail(3))
      if(isw.eq.1)then
        deallocate(xxo,yyo,zzo,gxo,gyo,gzo,stat=fail(4))
      endif
      if(ntcons.gt.0)then

        deallocate(dxx,dyy,dzz,txx,tyy,tzz,stat=fail(5))
        deallocate(dxt,dyt,dzt,xxt,yyt,zzt,stat=fail(6))

      endif
      deallocate(fxo,fyo,fzo,stat=fail(7))
      
      return
      end subroutine nptqvv_h2

      subroutine nstqvv_b2
     x  (safe,lshmov,isw,imcon,idnode,mxnode,natms,ngrp,nscons,
     x  ntcons,ntfree,ntpatm,mode,engke,engrot,press,taut,taup,sigma,
     x  tolnce,tstep,vircom,vircon,elrc,virlrc,volm)
      
c***********************************************************************
c     
c     dlpoly subroutine to integrate rigid body equations of motion
c     using the symplectic no_squish quaternion algorithm of 
c     miller et al j.chem.phys 116 (2002) 8649
c     plus bond constraints -  including rigid body sites linked
c     by constraint sites (qrattle algorithm)
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
c     author      w.smith   sep 2005
c     
c**********************************************************************

      implicit none

      logical newstep,safe,lshmov,newjob

      integer isw,imcon,idnode,mxnode,natms,ngrp,nscons,ntcons
      integer ntfree,i,j,k,igrp1,igrp2,igrp,ifre1,ifre2,jg,jr
      integer id,ifre,icyc,mxshk,idum,ig,ntpatm,mode
      integer iter,mxiter,jrs

      real(8) engke,engrot,tolnce,tstep,vircom,vircon,engtke
      real(8) engtrn,taut,sigma,chit
      real(8) vaa,vbb,vcc,opx,opy,opz,ftx,fty,ftz,volm0
      real(8) fmx,fmy,fmz,tqx,tqy,tqz,tq0,tq1,tq2,tq3
      real(8) volm,elrc0,elrc,virlrc0,virlrc
      real(8) beta,press,taup,engfke,hstep,vzero

      integer, parameter :: nnn=13
      real(8), parameter :: pt5=0.5d0

      integer fail(nnn)
      real(8) rot(9),strkin(9),strgrp(9),uni(9),celp(10)
      real(8) cell0(9)

      real(8), allocatable :: dens0(:)
      real(8), allocatable :: xxo(:),yyo(:),zzo(:)
      real(8), allocatable :: vxo(:),vyo(:),vzo(:)
      real(8), allocatable :: gxo(:),gyo(:),gzo(:)
      real(8), allocatable :: dxx(:),dyy(:),dzz(:)
      real(8), allocatable :: txx(:),tyy(:),tzz(:)
      real(8), allocatable :: xxt(:),yyt(:),zzt(:)
      real(8), allocatable :: dxt(:),dyt(:),dzt(:)
      real(8), allocatable :: dtx(:),dty(:),dtz(:)
      real(8), allocatable :: fxo(:),fyo(:),fzo(:)
      real(8), allocatable :: b0(:),b1(:),b2(:),b3(:)
      real(8), allocatable :: c0(:),c1(:),c2(:),c3(:)
      real(8), allocatable :: p0(:),p1(:),p2(:),p3(:)
      real(8), allocatable :: gvxo(:),gvyo(:),gvzo(:)

      save newstep,newjob,volm0,elrc0,virlrc0,dens0
      save p0,p1,p2,p3,ifre1,ifre2,igrp1,igrp2,hstep

      data newjob/.true./,beta/7.3728d-3/
      data uni/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/
      
c     set array alocation error flags
      
      do i=1,nnn
        fail(i)=0
      enddo
      
c     assign initial parameters

      if(newjob)then

c     timestep parameters
        
        hstep=pt5*tstep
        
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
        
        allocate(dens0(mxatyp),stat=fail(1))
        if(fail(1).gt.0)call error(idnode,2260)
        
c     store initial values of volume and long range corrections
        
        volm0=volm
        elrc0=elrc
        virlrc0=virlrc
        do i=1,ntpatm
          dens0(i)=dens(i)
        enddo
        
        allocate(p0(mxgrp),p1(mxgrp),p2(mxgrp),p3(mxgrp),stat=fail(1))

        newjob=.false.

      endif
      
c     allocate working arrays
      
      if(isw.eq.1)then

        allocate(xxo(mxatms),yyo(mxatms),zzo(mxatms),stat=fail(2))
        allocate(gxo(msgrp),gyo(msgrp),gzo(msgrp),stat=fail(3))

      endif
      if(ntcons.gt.0)then

        allocate(dxx(mxcons),dyy(mxcons),dzz(mxcons),stat=fail(4))
        allocate(txx(mxatms),tyy(mxatms),tzz(mxatms),stat=fail(5))
        allocate(xxt(mxatms),yyt(mxatms),zzt(mxatms),stat=fail(6))
        allocate(dxt(mxcons),dyt(mxcons),dzt(mxcons),stat=fail(7))

      endif
      allocate(fxo(mxatms),fyo(mxatms),fzo(mxatms),stat=fail(8))
      allocate(dtx(mxatms),dty(mxatms),dtz(mxatms),stat=fail(9))
      allocate(vxo(msatms),vyo(msatms),vzo(msatms),stat=fail(10))
      allocate(gvxo(msgrp),gvyo(msgrp),gvzo(msgrp),stat=fail(11))
      allocate(b0(msgrp),b1(msgrp),b2(msgrp),b3(msgrp),
     x  stat=fail(12))
      allocate(c0(msgrp),c1(msgrp),c2(msgrp),c3(msgrp),
     x  stat=fail(13))

c     check array allocation error flags

      do i=1,nnn
        if(fail(i).gt.0)call error(idnode,2270)
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
      
c     store key config data at start of time step
      
      if(isw.eq.1)then

c     cell parameters
        
        vzero=volm
        do i=1,9
          cell0(i)=cell(i)
        enddo
        
c     atom positions

        do i=1,natms
          
          xxo(i)=xxx(i)
          yyo(i)=yyy(i)
          zzo(i)=zzz(i)
          
        enddo

c     rigid body positions

        j=0
        do i=igrp1,igrp2
          
          j=j+1
          gxo(j)=gcmx(i)
          gyo(j)=gcmy(i)
          gzo(j)=gcmz(i)

        enddo

      endif

c     store free atom velocities
      
      j=0
      do ifre=ifre1,ifre2
        
        j=j+1
        i=lstfre(ifre)
        vxo(j)=vxx(i)
        vyo(j)=vyy(i)
        vzo(j)=vzz(i)
        
      enddo

c     store rigid body quaternions, momenta and cartesian velocities

      j=0
      do i=igrp1,igrp2

        j=j+1
        b0(j)=q0(i)
        b1(j)=q1(i)
        b2(j)=q2(i)
        b3(j)=q3(i)
        c0(j)=p0(i)
        c1(j)=p1(i)
        c2(j)=p2(i)
        c3(j)=p3(i)
        gvxo(j)=gvxx(i)
        gvyo(j)=gvyy(i)
        gvzo(j)=gvzz(i)

      enddo

c     store forces

      do i=1,natms
        
        fxo(i)=fxx(i)
        fyo(i)=fyy(i)
        fzo(i)=fzz(i)
        
      enddo
      
c     extract previous constraint terms from stress tensor

      if(isw.eq.1)then
        
        do i=1,9
          stress(i)=stress(i)-strcns(i)
        enddo

      endif

c     -------------- start of barostat iteration cycle -----------------
      
      mxiter=1
      if(isw.eq.1.and.ntcons.gt.0)mxiter=3
      do iter=1,mxiter
        
        do i=1,9
          eta(i)=uni(i)
        enddo
        
        if(isw.eq.1)then
          
c     restore cell parameters
          
          volm=vzero
          do i=1,9
            cell(i)=cell0(i)
          enddo
          
c     calculate Berendsen barostat
          
          do i=1,9
            eta(i)=tstep*beta*((stress(i)+strcns(i))/volm-
     x        press*uni(i))/taup+uni(i)
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

c     reset constraint virial
        
          vircon=0.d0
          
          do i=1,9
            strcns(i)=0.d0
          enddo
          
        endif
        
c     -------------- start of shake iteration cycle -------------------
        
        icyc=0
        mxshk=1
        safe=.false.
        newstep=.true.
        if(ntcons.gt.0)mxshk=mxshak
        do while(.not.safe.and.icyc.lt.mxshk)
          
          icyc=icyc+1
          
c     update velocities of free atoms 1/2 timestep
          
          j=0
          do ifre=ifre1,ifre2
            
            j=j+1
            i=lstfre(ifre)
            vxx(i)=vxo(j)+hstep*rmass(i)*fxx(i)
            vyy(i)=vyo(j)+hstep*rmass(i)*fyy(i)
            vzz(i)=vzo(j)+hstep*rmass(i)*fzz(i)
            
          enddo
          
c     *************  rigid body motion ****************************
          
c     operations common to first and second stages
          
          jg=0
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
            
            jg=jg+1
            call getrotmat(b0(jg),b1(jg),b2(jg),b3(jg),rot)
            
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
            
            tq0=2.0d0*(-b1(jg)*tqx-b2(jg)*tqy-b3(jg)*tqz)
            tq1=2.0d0*( b0(jg)*tqx-b3(jg)*tqy+b2(jg)*tqz)
            tq2=2.0d0*( b3(jg)*tqx+b0(jg)*tqy-b1(jg)*tqz)
            tq3=2.0d0*(-b2(jg)*tqx+b1(jg)*tqy+b0(jg)*tqz)
            
c     update quaternion momentum by half timestep
            
            p0(ig)=c0(jg)+tq0*hstep
            p1(ig)=c1(jg)+tq1*hstep
            p2(ig)=c2(jg)+tq2*hstep
            p3(ig)=c3(jg)+tq3*hstep
            
c     update centre of mass velocity by half timestep
            
            gvxx(ig)=gvxo(jg)+fmx*hstep/gmass(id)
            gvyy(ig)=gvyo(jg)+fmy*hstep/gmass(id)
            gvzz(ig)=gvzo(jg)+fmz*hstep/gmass(id)
            
          enddo
          
c     first stage of velocity verlet algorithm
          
          if(isw.eq.1)then
            
            jg=0
            do ig=igrp1,igrp2
              
              jg=jg+1
              
c     update centre of mass position by full time step
              
              gcmx(ig)=tstep*gvxx(ig)+
     x          eta(1)*gxo(jg)+eta(4)*gyo(jg)+eta(7)*gzo(jg)
              gcmy(ig)=tstep*gvyy(ig)+
     x          eta(2)*gxo(jg)+eta(5)*gyo(jg)+eta(8)*gzo(jg)
              gcmz(ig)=tstep*gvzz(ig)+
     x          eta(3)*gxo(jg)+eta(6)*gyo(jg)+eta(9)*gzo(jg)
              
c     calculate rotation of rigid groups: nosquish algorithm
              
              q0(ig)=b0(jg)
              q1(ig)=b1(jg)
              q2(ig)=b2(jg)
              q3(ig)=b3(jg)
              call nosquish(ig,tstep,q0,q1,q2,q3,p0,p1,p2,p3)
              
            enddo

c     new atomic positions for atoms in rigid bodies - relative to com
          
            k=0
            do ig=igrp1,igrp2
              
              call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)

              id=lstgtp(ig)
              do j=1,numgsit(id)
                
                k=k+1
                i=lstme(k)
                xxx(i)=rot(1)*gxx(id,j)+rot(2)*gyy(id,j)+
     x            rot(3)*gzz(id,j)+gcmx(ig)
                yyy(i)=rot(4)*gxx(id,j)+rot(5)*gyy(id,j)+
     x            rot(6)*gzz(id,j)+gcmy(ig)
                zzz(i)=rot(7)*gxx(id,j)+rot(8)*gyy(id,j)+
     x            rot(9)*gzz(id,j)+gcmz(ig)

              enddo
              
            enddo
            
c     update positions of free particles to full time step
            
            do ifre=ifre1,ifre2
              
              i=lstfre(ifre)
              xxx(i)=tstep*vxx(i)+
     x          eta(1)*xxo(i)+eta(4)*yyo(i)+eta(7)*zzo(i)
              yyy(i)=tstep*vyy(i)+
     x          eta(2)*xxo(i)+eta(5)*yyo(i)+eta(8)*zzo(i)
              zzz(i)=tstep*vzz(i)+
     x          eta(3)*xxo(i)+eta(6)*yyo(i)+eta(9)*zzo(i)
              
            enddo
            
c     merge atom positions
            
            if(mxnode.gt.1)call merge1
     x        (idnode,mxnode,natms,lstme,xxx,yyy,zzz,buffer)
            
c     apply rattle corrections to bond constraints
            
            if(ntcons.gt.0)then
              
              if(mxnode.gt.1)call merge1
     x            (idnode,mxnode,natms,lstme,vxx,vyy,vzz,buffer)

              call qrattle_r
     x          (newstep,safe,lshmov,idnode,imcon,mxnode,natms,
     x          nscons,tolnce,tstep,vircon,dxx,dyy,dzz,dxt,dyt,
     x          dzt,txx,tyy,tzz,xxt,yyt,zzt,strcns)
              
            endif
            
c     end of first stage 
            
          endif

c     second stage of velocity verlet algorithm
          
          if(isw.eq.2)then
            
            jr=0
            do ig=igrp1,igrp2
              
c     new angular momenta and velocities
              
              opx=pt5*(-q1(ig)*p0(ig)+q0(ig)*p1(ig)+
     x          q3(ig)*p2(ig)-q2(ig)*p3(ig))
              opy=pt5*(-q2(ig)*p0(ig)-q3(ig)*p1(ig)+
     x          q0(ig)*p2(ig)+q1(ig)*p3(ig))
              opz=pt5*(-q3(ig)*p0(ig)+q2(ig)*p1(ig)-
     x          q1(ig)*p2(ig)+q0(ig)*p3(ig))
              
              id=lstgtp(ig)
              
              omx(ig)=opx*rotinx(id,2)
              omy(ig)=opy*rotiny(id,2)
              omz(ig)=opz*rotinz(id,2)
              
c     new rotational matrix
              
              call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)
              
c     new atomic velocites
              
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
            
            if(ntcons.gt.0)then
            
              if(mxnode.gt.1)then
                
                call merge1
     x            (idnode,mxnode,natms,lstme,vxx,vyy,vzz,buffer)
                call merge1
     x            (idnode,mxnode,natms,lstme,fxx,fyy,fzz,buffer)
                
              endif
              
c     correct constraint bond velocities using rattle
              
              call qrattle_v
     x          (newstep,safe,lshmov,idnode,mxnode,natms,
     x          nscons,tolnce,tstep,dxx,dyy,dzz,txx,tyy,tzz,
     x          xxt,yyt,zzt)
              
            endif
            
c     end of second stage
            
          endif
          
          newstep=.false.
          
        enddo

c     check shake convergence
        
        if(.not.safe)call error(idnode,105)
        
c     sum constraint virial and stress across processors
        
        if(mxnode.gt.1.and.isw.eq.1)then
          
          buffer(1)=vircon
          call gdsum(buffer(1),1,buffer(2))
          vircon=buffer(1)
          call gdsum(strcns,9,buffer)
          
        endif

c     -------------- end of shake iteration cycle -------------------
      
c     rigid body contribution to stress tensor

        if(isw.eq.2)call bodystress
     x    (idnode,mxnode,ngrp,vircom,strbod,dtx,dty,dtz)

c     restore  forces
        
        do i=1,natms
          
          fxx(i)=fxo(i)
          fyy(i)=fyo(i)
          fzz(i)=fzo(i)
          
        enddo

c     -------------- end of barostat iteration cycle ----------------
        
      enddo
      
      if(isw.eq.1)then

c     adjust long range corrections and number density
        
        elrc=elrc0*(volm0/volm)
        virlrc=virlrc0*(volm0/volm)
        
        do k=1,ntpatm
          dens(k)=dens0(k)*(volm0/volm)
        enddo
        
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

c     thermostat velocities of free particles
        
        do ifre=ifre1,ifre2
          
          i=lstfre(ifre)
          if(lstfrz(i).ne.0)then
            
            vxx(i)=chit*vxx(i)
            vyy(i)=chit*vyy(i)
            vzz(i)=chit*vzz(i)
            
          endif
          
        enddo
        
c     thermostat rigid body velocities
        
        do ig=igrp1,igrp2
          
          id=lstgtp(ig)
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
      
      if(mxnode.gt.1)then

c     merge new group coordinates and velocities

        if(isw.eq.1)
     x    call merge(idnode,mxnode,ngrp,mxbuff,gcmx,gcmy,gcmz,buffer)
        call merge(idnode,mxnode,ngrp,mxbuff,gvxx,gvyy,gvzz,buffer)
        
c     merge new atomic coordinates and velocities
        
        if(isw.eq.1)
     x    call merge1(idnode,mxnode,natms,lstme,xxx,yyy,zzz,buffer)
        call merge1(idnode,mxnode,natms,lstme,vxx,vyy,vzz,buffer)
        
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

      deallocate(gvxo,gvyo,gvzo,vxo,vyo,vzo,stat=fail(1))
      deallocate(b0,b1,b2,b3,c0,c1,c2,c3,stat=fail(2))
      deallocate(dtx,dty,dtz,fxo,fyo,fzo,stat=fail(3))
      if(isw.eq.1)then
        deallocate(xxo,yyo,zzo,gxo,gyo,gzo,stat=fail(4))
      endif
      if(ntcons.gt.0)then

        deallocate(dxx,dyy,dzz,txx,tyy,tzz,stat=fail(5))
        deallocate(dxt,dyt,dzt,xxt,yyt,zzt,stat=fail(6))

      endif
      
      return
      end subroutine nstqvv_b2

      subroutine nstqvv_h2
     x  (safe,lshmov,isw,idnode,mxnode,natms,imcon,ngrp,nscons,
     x  ntcons,ntpatm,ntfree,mode,ntshl,keyshl,tstep,taut,taup,
     x  sigma,temp,chit,consv,conint,engke,engrot,elrc,tolnce,
     x  vircom,vircon,virlrc,volm,press,chit_shl,sigma_shl)
      
c***********************************************************************
c     
c     dlpoly subroutine to integrate rigid body equations of motion
c     using the symplectic no_squish quaternion algorithm of 
c     miller et al j.chem.phys 116 (2002) 8649
c     plus bond constraints -  including rigid body sites linked
c     by constraint sites (qrattle algorithm)
c     
c     nst ensemble - nose-hoover thermostat Molec Phys 87 (1996) 1117
c     
c     parallel replicated data version : block data
c     
c     omx,omy,omz=angular velocity in body fixed frame (principal axes)
c     rotinx,y,z =rotational inertia in body fixed frame
c     
c     copyright daresbury laboratory
c     author      w.smith sept 2005
c     adapted     d. quigley : metadynamics
c     
c**********************************************************************

      implicit none

      logical newstep,safe,lshmov,newjob

      integer isw,imcon,idnode,mxnode,natms,ngrp,nscons,ntcons
      integer ntfree,i,j,k,igrp1,igrp2,igrp,ifre1,ifre2,jr
      integer id,ifre,icyc,mxshk,idum,ig,ntpatm,mode
      integer jcyc,iter,mxiter

      real(8) engke,engrot,tolnce,tstep,vircom,vircon
      real(8) engtrn,engfke
      real(8) vaa,vbb,vcc,opx,opy,opz
      real(8) taut,taup,sigma,chit,hstep,qmass,conint,consv
      real(8) cxx,cyy,czz,press,chip2,temp
      real(8) volm,pmass,totmas,qstep,fstep,volm0,elrc
      real(8) virlrc,elrc0,virlrc0,chit0,vzero,cons0

      integer, parameter :: nnn=16
      integer, parameter :: ncyc=5
      real(8), parameter :: pt5=0.5d0

      integer fail(nnn)
      real(8) rot(9),strkin(9),strgrp(9),com(3),vom(3)
      real(8) cell0(9),eta0(9),stres0(9)

      real(8), allocatable :: dens0(:)
      real(8), allocatable :: xxo(:),yyo(:),zzo(:)
      real(8), allocatable :: vxo(:),vyo(:),vzo(:)
      real(8), allocatable :: vx1(:),vy1(:),vz1(:)
      real(8), allocatable :: fxo(:),fyo(:),fzo(:)
      real(8), allocatable :: gxo(:),gyo(:),gzo(:)
      real(8), allocatable :: oxo(:),oyo(:),ozo(:)
      real(8), allocatable :: dxx(:),dyy(:),dzz(:)
      real(8), allocatable :: txx(:),tyy(:),tzz(:)
      real(8), allocatable :: xxt(:),yyt(:),zzt(:)
      real(8), allocatable :: dxt(:),dyt(:),dzt(:)
      real(8), allocatable :: dtx(:),dty(:),dtz(:)
      real(8), allocatable :: b0(:),b1(:),b2(:),b3(:)
      real(8), allocatable :: c0(:),c1(:),c2(:),c3(:)
      real(8), allocatable :: p0(:),p1(:),p2(:),p3(:)
      real(8), allocatable :: gvxo(:),gvyo(:),gvzo(:)
      real(8), allocatable :: gvx1(:),gvy1(:),gvz1(:)
      
c     metadynamics shell thermostat variables
      
      integer ntshl,keyshl
      real(8) sigma_shl
      
      logical,save :: lfirst=.true.
      real(8)      :: chit_shl  
      real(8),save :: qmass_shl
      real(8)      :: shlke
      
c     end metadynamics shell thermostat variables
      
      save newstep,newjob,p0,p1,p2,p3,hstep,fstep,qmass,ifre1,ifre2
      save igrp1,igrp2,volm0,elrc0,virlrc0,qstep,dens0,totmas
      save pmass

      data newjob/.true./
      
c     set array allocation error flags

      do i=1,nnn
        fail(i)=0
      enddo

c     assign initial parameters

      if(newjob)then
        
c     store intitial parameters
        
        volm0=volm
        elrc0=elrc
        virlrc0=virlrc
        hstep=0.5d0*tstep
        fstep=0.5d0*tstep/dble(ncyc)
        qstep=0.25d0*tstep/dble(ncyc)
        
c     allocate density storage array

        fail(1)=0
        allocate(dens0(mxatyp),stat=fail(1))
        if(fail(1).gt.0)call error(idnode,2220)

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
          call  error(idnode,506)
        endif

        allocate(p0(mxgrp),p1(mxgrp),p2(mxgrp),p3(mxgrp),stat=fail(1))
        newjob=.false.

      endif

      if(ntcons.gt.0)safe=.false.

c     allocate working arrays

      if(isw.eq.1)then

        allocate(xxo(mxatms),yyo(mxatms),zzo(mxatms),stat=fail(2))
        allocate(gxo(mxgrp),gyo(mxgrp),gzo(mxgrp),stat=fail(3))

      endif
      if(ntcons.gt.0)then
        
        allocate(dxx(mxcons),dyy(mxcons),dzz(mxcons),stat=fail(4))
        allocate(txx(mxatms),tyy(mxatms),tzz(mxatms),stat=fail(5))
        allocate(xxt(mxatms),yyt(mxatms),zzt(mxatms),stat=fail(6))
        allocate(dxt(mxcons),dyt(mxcons),dzt(mxcons),stat=fail(7))
        
      endif
      allocate(dtx(mxatms),dty(mxatms),dtz(mxatms),stat=fail(8))
      allocate(vxo(mxatms),vyo(mxatms),vzo(mxatms),stat=fail(9))
      allocate(vx1(mxatms),vy1(mxatms),vz1(mxatms),stat=fail(10))
      allocate(fxo(mxatms),fyo(mxatms),fzo(mxatms),stat=fail(11))
      allocate(oxo(mxatms),oyo(mxatms),ozo(mxatms),stat=fail(12))
      allocate(b0(mxgrp),b1(mxgrp),b2(mxgrp),b3(mxgrp),
     x  stat=fail(13))
      allocate(c0(mxgrp),c1(mxgrp),c2(mxgrp),c3(mxgrp),
     x  stat=fail(14))
      allocate(gvxo(mxgrp),gvyo(mxgrp),gvzo(mxgrp),stat=fail(15))
      allocate(gvx1(mxgrp),gvy1(mxgrp),gvz1(mxgrp),stat=fail(16))

      do i=1,nnn
        if(fail(i).gt.0)call error(idnode,2230)
      enddo

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
      
c     construct current bond vectors
      
      if(ntcons.gt.0)then

        do k=1,nscons
          
          i=listcon(k,2)
          j=listcon(k,3)
          dxx(k)=xxx(i)-xxx(j)
          dyy(k)=yyy(i)-yyy(j)
          dzz(k)=zzz(i)-zzz(j)
          
        enddo

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
      
c     store key config data at start of timestep

      vzero=volm
      chit0=chit
      cons0=conint
      do i=1,9

        cell0(i)=cell(i)
        eta0(i)=eta(i)
        stres0(i)=stress(i)
        
      enddo

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

c     store atom positions
        
        do i=1,natms
          
          xxo(i)=xxx(i)
          yyo(i)=yyy(i)
          zzo(i)=zzz(i)
          
        enddo
        
c     store rigid body positions
        
        do ig=1,ngrp
          
          gxo(ig)=gcmx(ig)
          gyo(ig)=gcmy(ig)
          gzo(ig)=gcmz(ig)
          
        enddo
        
      endif
      
c     store free atom velocities

      do i=1,natms
        
        vxo(i)=vxx(i)
        vyo(i)=vyy(i)
        vzo(i)=vzz(i)
        
      enddo

c     store forces

      do i=1,natms
        
        fxo(i)=fxx(i)
        fyo(i)=fyy(i)
        fzo(i)=fzz(i)
        
      enddo
      
c     store rigid body quaternions, angular and cartesian velocities

      do ig=1,ngrp
        
        b0(ig)=q0(ig)
        b1(ig)=q1(ig)
        b2(ig)=q2(ig)
        b3(ig)=q3(ig)
        oxo(ig)=omx(ig)
        oyo(ig)=omy(ig)
        ozo(ig)=omz(ig)
        gvxo(ig)=gvxx(ig)
        gvyo(ig)=gvyy(ig)
        gvzo(ig)=gvzz(ig)

      enddo
      
c     iteration necessary if ntcons > 0 and isw=1
      
      mxiter=1
      if(isw.eq.1.and.ntcons.gt.0)mxiter=3
      do iter=1,mxiter

c     integration of barostat and thermostat (part 1)

        if(isw.eq.1)then

c     restore cell parameters

          volm=vzero
          chit=chit0
          conint=cons0
          do i=1,9

            cell(i)=cell0(i)
            eta(i)=eta0(i)
            
          enddo

c     restore free atom velocities
          
          do i=1,natms
            
            vxx(i)=vxo(i)
            vyy(i)=vyo(i)
            vzz(i)=vzo(i)
            
          enddo

c     restore rigid body quaternions angular and cartesian velocities

          do ig=1,ngrp
            
            omx(ig)=oxo(ig)
            omy(ig)=oyo(ig)
            omz(ig)=ozo(ig)
            gvxx(ig)=gvxo(ig)
            gvyy(ig)=gvyo(ig)
            gvzz(ig)=gvzo(ig)
            
          enddo
          
c     kinetic contributions to stress tensor
          
          call kinstressf(ntfree,idnode,mxnode,strkin)        
          call kinstressg(ngrp,idnode,mxnode,strgrp)

          do jcyc=1,ncyc
            
c     integrate and apply nst thermostat
            
            call nstqscl_t2
     x        (idnode,mxnode,ntfree,ngrp,mode,engfke,engtrn,engrot,temp,
     x        sigma,qstep,pmass,qmass,taut,chit,conint,strkin,strgrp)
            
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
            
            call nstqscl_p2
     x        (idnode,mxnode,ntfree,ngrp,mode,fstep,pmass,chit,press,
     x        volm,strkin,strgrp)
            
c     integrate and apply nst thermostat
            
            call nstqscl_t2
     x        (idnode,mxnode,ntfree,ngrp,mode,engfke,engtrn,engrot,temp,
     x        sigma,qstep,pmass,qmass,taut,chit,conint,strkin,strgrp)
            
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
          
c     translational kinetic energy

          engke=engfke+engtrn
          
c     reset constraint virial and stress

          vircon=0.d0
          do i=1,9

            stress(i)=stress(i)-strcns(i)-strbod(i)
            strcns(i)=0.d0

          enddo
          
c     calculate quaternion momenta
          
          do ig=igrp1,igrp2
            
            id=lstgtp(ig)
            opx=omx(ig)*rotinx(id,1)
            opy=omy(ig)*rotiny(id,1)
            opz=omz(ig)*rotinz(id,1)
            p0(ig)=2.0d0*(-b1(ig)*opx-b2(ig)*opy-b3(ig)*opz)
            p1(ig)=2.0d0*( b0(ig)*opx-b3(ig)*opy+b2(ig)*opz)
            p2(ig)=2.0d0*( b3(ig)*opx+b0(ig)*opy-b1(ig)*opz)
            p3(ig)=2.0d0*(-b2(ig)*opx+b1(ig)*opy+b0(ig)*opz)
            
          enddo
          
        endif
        
c     store intermediate velocities

        do ifre=ifre1,ifre2
          
          i=lstfre(ifre)
          vx1(i)=vxx(i)
          vy1(i)=vyy(i)
          vz1(i)=vzz(i)
          
        enddo
        do ig=igrp1,igrp2
          
          c0(ig)=p0(ig)
          c1(ig)=p1(ig)
          c2(ig)=p2(ig)
          c3(ig)=p3(ig)
          gvx1(ig)=gvxx(ig)
          gvy1(ig)=gvyy(ig)
          gvz1(ig)=gvzz(ig)
          
        enddo

c     -------------- start of shake iteration cycle -------------------
        
        icyc=0
        mxshk=1
        safe=.false.
        newstep=.true.
        if(ntcons.gt.0)mxshk=mxshak
        
        do while(.not.safe.and.icyc.lt.mxshk)
          
          icyc=icyc+1
          
c     update velocities of free atoms 1/2 timestep
          
          do ifre=ifre1,ifre2
            
            i=lstfre(ifre)
            vxx(i)=vx1(i)+hstep*rmass(i)*fxx(i)
            vyy(i)=vy1(i)+hstep*rmass(i)*fyy(i)
            vzz(i)=vz1(i)+hstep*rmass(i)*fzz(i)
            
          enddo
          
c     *************  rigid body motion ****************************
          
c     restore rigid body quaternions, angular momenta and velocities

          do ig=igrp1,igrp2
            
            q0(ig)=b0(ig)
            q1(ig)=b1(ig)
            q2(ig)=b2(ig)
            q3(ig)=b3(ig)
            p0(ig)=c0(ig)
            p1(ig)=c1(ig)
            p2(ig)=c2(ig)
            p3(ig)=c3(ig)
            gvxx(ig)=gvx1(ig)
            gvyy(ig)=gvy1(ig)
            gvzz(ig)=gvz1(ig)
            
          enddo
          
c     calculate new rigid body velocities

          call rotate_omega
     x      (idnode,mxnode,ngrp,hstep,p0,p1,p2,p3,dtx,dty,dtz)

c     first stage of velocity verlet algorithm

          if(isw.eq.1)then
            
c     calculate system centre of mass
            
            call getcom(natms,idnode,mxnode,totmas,com)

c     update centre of mass position by full time step
            
            do ig=igrp1,igrp2

              cxx=gxo(ig)-com(1)
              cyy=gyo(ig)-com(2)
              czz=gzo(ig)-com(3)
              gcmx(ig)=gxo(ig)+
     x          tstep*(gvxx(ig)+eta(1)*cxx+eta(4)*cyy+eta(7)*czz)
              gcmy(ig)=gyo(ig)+
     x          tstep*(gvyy(ig)+eta(2)*cxx+eta(5)*cyy+eta(8)*czz)
              gcmz(ig)=gzo(ig)+
     x          tstep*(gvzz(ig)+eta(3)*cxx+eta(6)*cyy+eta(9)*czz)

c     calculate rotation of rigid groups: nosquish algorithm
              
              call nosquish(ig,tstep,q0,q1,q2,q3,p0,p1,p2,p3)
              
            enddo
            
c     merge group coms from all nodes

            if(mxnode.gt.1)call merge
     x        (idnode,mxnode,ngrp,mxbuff,gcmx,gcmy,gcmz,buffer)

c     new atomic positions for atoms in rigid bodies - relative to com
            
            k=0
            do ig=igrp1,igrp2
              
              call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)
              
              id=lstgtp(ig)
              do j=1,numgsit(id)
                
                k=k+1
                i=lstme(k)
                xxx(i)=rot(1)*gxx(id,j)+rot(2)*gyy(id,j)+
     x            rot(3)*gzz(id,j)+gcmx(ig)
                yyy(i)=rot(4)*gxx(id,j)+rot(5)*gyy(id,j)+
     x            rot(6)*gzz(id,j)+gcmy(ig)
                zzz(i)=rot(7)*gxx(id,j)+rot(8)*gyy(id,j)+
     x            rot(9)*gzz(id,j)+gcmz(ig)
                
              enddo
              
            enddo
            
c     update positions of free particles to full time step
            
            do ifre=ifre1,ifre2
              
              i=lstfre(ifre)
              cxx=xxo(i)-com(1)
              cyy=yyo(i)-com(2)
              czz=zzo(i)-com(3)
              xxx(i)=xxo(i)+
     x          tstep*(vxx(i)+eta(1)*cxx+eta(4)*cyy+eta(7)*czz)
              yyy(i)=yyo(i)+
     x          tstep*(vyy(i)+eta(2)*cxx+eta(5)*cyy+eta(8)*czz)
              zzz(i)=zzo(i)+
     x          tstep*(vzz(i)+eta(3)*cxx+eta(6)*cyy+eta(9)*czz)
              
            enddo

c     merge atom positions

            if(mxnode.gt.1)call merge1
     x        (idnode,mxnode,natms,lstme,xxx,yyy,zzz,buffer)

c     apply rattle corrections to bond constraints
            
            if(ntcons.gt.0)then
              
              if(mxnode.gt.1)call merge1
     x          (idnode,mxnode,natms,lstme,vxx,vyy,vzz,buffer)

              call qrattle_r
     x          (newstep,safe,lshmov,idnode,imcon,mxnode,natms,
     x          nscons,tolnce,tstep,vircon,dxx,dyy,dzz,dxt,dyt,
     x          dzt,txx,tyy,tzz,xxt,yyt,zzt,strcns)
              
            endif

c     end of first stage 

          endif

c     second stage of velocity verlet algorithm
          
          if(isw.eq.2)then
            
            jr=0
            do ig=igrp1,igrp2
              
c     new angular momenta and velocities

              opx=pt5*(-q1(ig)*p0(ig)+q0(ig)*p1(ig)+
     x          q3(ig)*p2(ig)-q2(ig)*p3(ig))
              opy=pt5*(-q2(ig)*p0(ig)-q3(ig)*p1(ig)+
     x          q0(ig)*p2(ig)+q1(ig)*p3(ig))
              opz=pt5*(-q3(ig)*p0(ig)+q2(ig)*p1(ig)-
     x          q1(ig)*p2(ig)+q0(ig)*p3(ig))

              id=lstgtp(ig)

              omx(ig)=opx*rotinx(id,2)
              omy(ig)=opy*rotiny(id,2)
              omz(ig)=opz*rotinz(id,2)
              
c     new rotational matrix
              
              call getrotmat(q0(ig),q1(ig),q2(ig),q3(ig),rot)

c     new atomic velocites
              
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
            
            if(ntcons.gt.0)then
              
              if(mxnode.gt.1)then

                call merge1
     x            (idnode,mxnode,natms,lstme,vxx,vyy,vzz,buffer)
                call merge1
     x            (idnode,mxnode,natms,lstme,fxx,fyy,fzz,buffer)
                
              endif

c     correct constraint bond velocities using rattle

              call qrattle_v
     x          (newstep,safe,lshmov,idnode,mxnode,natms,
     x          nscons,tolnce,tstep,dxx,dyy,dzz,txx,tyy,tzz,
     x          xxt,yyt,zzt)

            endif
            
c     end of second stage
            
          endif
          
          newstep=.false.
          
        enddo

c     check shake convergence
        
        if(.not.safe)call error(idnode,105)
        
c     sum constraint virial and stress across processors
        
        if(mxnode.gt.1.and.isw.eq.1)then
          
          buffer(1)=vircon
          call gdsum(buffer(1),1,buffer(2))
          vircon=buffer(1)
          call gdsum(strcns,9,buffer)
          
        endif
        
c     -------------- end of shake iteration cycle -------------------

c     calculate rigid body contribution to stress tensor
          
        call bodystress(idnode,mxnode,ngrp,vircom,strbod,dtx,dty,dtz)
          
c     add new constraint terms to stress tensor
          
        if(isw.eq.1)then

          do i=1,9
            stress(i)=stress(i)+strcns(i)+strbod(i)
          enddo
        
        endif

c     integration of barostat and thermostat (part 2)
        
        if(isw.eq.2)then
          
c     kinetic contributions to stress tensor
          
          call kinstressf(ntfree,idnode,mxnode,strkin)        
          call kinstressg(ngrp,idnode,mxnode,strgrp)
          
c     add kinetic and body contributions to stress tensor
          
          do i=1,9
            stress(i)=stres0(i)+strkin(i)+strgrp(i)+strbod(i)+strcns(i)
          enddo

          do jcyc=1,ncyc
            
c     integrate and apply nst thermostat
            
            call nstqscl_t2
     x        (idnode,mxnode,ntfree,ngrp,mode,engfke,engtrn,engrot,temp,
     x        sigma,qstep,pmass,qmass,taut,chit,conint,strkin,strgrp)
            
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
            
            call nstqscl_p2
     x        (idnode,mxnode,ntfree,ngrp,mode,fstep,pmass,chit,press,
     x        volm,strkin,strgrp)
            
c     integrate and apply nst thermostat
            
            call nstqscl_t2
     x        (idnode,mxnode,ntfree,ngrp,mode,engfke,engtrn,engrot,temp,
     x        sigma,qstep,pmass,qmass,taut,chit,conint,strkin,strgrp)
            
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
          
c     translational kinetic energy
          
          engke=engfke+engtrn

c     sum up all contributions to stress tensor
          
          do i=1,9
            stress(i)=stres0(i)+strkin(i)+strgrp(i)+strbod(i)+strcns(i)
          enddo
          
        endif

c     restore forces

        do i=1,natms
          
          fxx(i)=fxo(i)
          fyy(i)=fyo(i)
          fzz(i)=fzo(i)
          
        enddo
        
c     -------------- end of barostat iteration cycle ----------------

      enddo

      if(isw.eq.2)then

c     calculate conserved variable

        chip2=sdot0(9,eta,eta)
        if(mode.eq.2)chip2=chip2-eta(1)**2
        consv=conint+0.5d0*qmass*chit**2+0.5d0*pmass*chip2+press*volm
        
c     metadynamics shell thermostat

        if(lmetadyn.and.keyshl.eq.1)then
          consv=consv+0.5d0*qmass_shl*chit_shl**2
        endif

c     merge velocity arrays
        
        if(mxnode.gt.1)call merge1
     x    (idnode,mxnode,natms,lstme,vxx,vyy,vzz,buffer)
        
      endif

c     adjust long range corrections and number density
      
      elrc=elrc0*(volm0/volm)
      virlrc=virlrc0*(volm0/volm)
      do k=1,ntpatm
        dens(k)=dens0(k)*(volm0/volm)
      enddo
      
      if(mxnode.gt.1)then

c     merge new group coordinates and velocities

        if(isw.eq.1)
     x    call merge(idnode,mxnode,ngrp,mxbuff,gcmx,gcmy,gcmz,buffer)
        call merge(idnode,mxnode,ngrp,mxbuff,gvxx,gvyy,gvzz,buffer)
        
c     merge new atomic coordinates and velocities
        
        if(isw.eq.1)
     x    call merge1(idnode,mxnode,natms,lstme,xxx,yyy,zzz,buffer)
        call merge1(idnode,mxnode,natms,lstme,vxx,vyy,vzz,buffer)
        
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

      deallocate(gvxo,gvyo,gvzo,vxo,vyo,vzo,stat=fail(1))
      deallocate(b0,b1,b2,b3,c0,c1,c2,c3,stat=fail(2))
      deallocate(dtx,dty,dtz,stat=fail(3))
      if(isw.eq.1)then
        deallocate(xxo,yyo,zzo,gxo,gyo,gzo,stat=fail(4))
      endif
      if(ntcons.gt.0)then

        deallocate(dxx,dyy,dzz,txx,tyy,tzz,stat=fail(5))
        deallocate(dxt,dyt,dzt,xxt,yyt,zzt,stat=fail(6))

      endif
      if(isw.eq.2)then
        deallocate(fxo,fyo,fzo,stat=fail(7))
      endif
      
      return
      end subroutine nstqvv_h2
      
      end module vv_rotation2_module
