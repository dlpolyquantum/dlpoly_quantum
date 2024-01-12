      module vv_motion_module

c***********************************************************************
c     
c     dl_poly module for velocity verlet integration schemes
c     copyright - daresbury laboratory
c     author    - w. smith    aug 2006
c     adapted   - d. quigley - metadynamics
c     
c     nvtvv_nhc and nptvv_nhc subroutines are added for NVT and
c     NPT ensembles together with Nose-Hoover Chain thermostat/barostat
c
c     copyright - M.R.Momeni and F.A. Shakib
c     authors   - M.R.Momeni and F.A.Shakib 2021
c
c     Method Development and Materials Simulation Laboratory
c
c***********************************************************************
      
      use config_module
      use ensemble_tools_module
      use error_module
      use metafreeze_module,only : lmetadyn
      use property_module
      use setup_module
      use shake_module
      use site_module
      use utility_module
      use nhc_module
      use water_module

      contains
      
      subroutine rdrattle_r
     x  (safe,lshmov,idnode,imcon,mxnode,natms,nscons,
     x  tolnce,tstep,vircon,dxx,dyy,dzz,dxt,dyt,dzt,
     x  txx,tyy,tzz,xxt,yyt,zzt,stresh)

c***********************************************************************
c     
c     dl_poly subroutine for applying bond constraint corrections after
c     atomic integration. rattle algorithm
c     must be used in conjunction with integration algorithms
c     
c     copyright - daresbury laboratory
c     author    - w. smith october 2002
c     amended   - w. smith january 2005 : f90 conversion
c     
c***********************************************************************

      implicit none
      
      logical check,safe,lshmov
      integer idnode,imcon,mxnode,natms,nscons,icyc,i,j,k
      real(8) tolnce,tstep,vircon,dxx,dyy,dzz,dxt,dyt,dzt,esig
      real(8) txx,tyy,tzz,xxt,yyt,zzt,stresh,dx,dy,dz,dis,omega2
      real(8) strs1,strs2,strs3,strs5,strs6,strs9,amti,amtj,gamma
      real(8) gammi,gammj,dli,dlj
      
      dimension stresh(9)
      dimension xxt(mxatms),yyt(mxatms),zzt(mxatms)
      dimension txx(mxatms),tyy(mxatms),tzz(mxatms)
      dimension dxx(mxcons),dyy(mxcons),dzz(mxcons)
      dimension dxt(mxcons),dyt(mxcons),dzt(mxcons)
      
c     constraint virial

      vircon=0.d0

c     accumulators for stress tensor

      strs1=0.d0
      strs2=0.d0
      strs3=0.d0
      strs5=0.d0
      strs6=0.d0
      strs9=0.d0

c     test size of work arrays

      check=.true.
      if(mxxdf.lt.nscons)check=.false.
      if(mxnode.gt.1)call gstate(check)
      if(.not.check)call error(idnode,412)

c     application of constraint (shake) algorithm
      
      icyc=0
      safe=.false.

      do while(.not.safe.and.icyc.lt.mxshak)

        icyc=icyc+1

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
          
          dx=dxt(k)
          dy=dyt(k)
          dz=dzt(k)
          dis=prmcon(listcon(k,1))
          esig=max(esig,abs(dx*dx+dy*dy+dz*dz-dis*dis)/dis)
          
        enddo
        
        esig=esig*0.5d0

c     global verification of convergence
        
        safe=(esig.lt.tolnce)
        
        if(mxnode.gt.1)call gstate(safe)

c     bypass calculations if all tolerances satisfied 
        
        if(.not.safe)then
          
c     initialise increment arrays
          
          do i=1,natms
            
            xxt(i)=0.d0
            yyt(i)=0.d0
            zzt(i)=0.d0
            
          enddo
          
c     calculate constraint forces
          
          do k=1,nscons
            
            i=listcon(k,2)
            j=listcon(k,3)
            
c     set constraint parameters
            
            dis=prmcon(listcon(k,1))
            omega2= dis*dis
            amti= tstep*rmass(i)
            amtj=-tstep*rmass(j)
            
c     constraint force parameter
            
            dx=dxt(k)
            dy=dyt(k)
            dz=dzt(k)
            
            gamma=(omega2-(dx*dx+dy*dy+dz*dz))/
     x        (-tstep*(amti-amtj)*(dxx(k)*dx+dyy(k)*dy+dzz(k)*dz))
            
c     accumulate bond virial
            
            vircon=vircon+gamma*(dxx(k)**2+dyy(k)**2+dzz(k)**2)
            
            strs1=strs1-gamma*dxx(k)*dxx(k)
            strs2=strs2-gamma*dxx(k)*dyy(k)
            strs3=strs3-gamma*dxx(k)*dzz(k)
            strs5=strs5-gamma*dyy(k)*dyy(k)
            strs6=strs6-gamma*dyy(k)*dzz(k)
            strs9=strs9-gamma*dzz(k)*dzz(k)
            
c     improve approximate constraint force
            
            gammi=-0.5d0*gamma*amti
            xxt(i)=xxt(i)+dxx(k)*gammi
            yyt(i)=yyt(i)+dyy(k)*gammi
            zzt(i)=zzt(i)+dzz(k)*gammi
            
            gammj=-0.5d0*gamma*amtj
            xxt(j)=xxt(j)+dxx(k)*gammj
            yyt(j)=yyt(j)+dyy(k)*gammj
            zzt(j)=zzt(j)+dzz(k)*gammj
            
          enddo
          
c     sum up constraint forces across nodes
          
          if(mxnode.gt.1)then
            
            if(lshmov)call shmove
     x        (idnode,mxnode,natms,lashap,lishap,xxt,yyt,zzt,
     x        txx,tyy,tzz,buffer)
            
          endif
          
          do k=1,nscons
            
            i=listcon(k,2)
            j=listcon(k,3)
            
            dli=1.0d0/dble(listme(i))
            dlj=1.0d0/dble(listme(j))
            
            xxx(i)=xxx(i)+tstep*dli*xxt(i)
            yyy(i)=yyy(i)+tstep*dli*yyt(i)
            zzz(i)=zzz(i)+tstep*dli*zzt(i)
            xxx(j)=xxx(j)+tstep*dlj*xxt(j)
            yyy(j)=yyy(j)+tstep*dlj*yyt(j)
            zzz(j)=zzz(j)+tstep*dlj*zzt(j)
            
            vxx(i)=vxx(i)+dli*xxt(i)
            vzz(i)=vzz(i)+dli*zzt(i)
            vyy(i)=vyy(i)+dli*yyt(i)
            vxx(j)=vxx(j)+dlj*xxt(j)
            vyy(j)=vyy(j)+dlj*yyt(j)
            vzz(j)=vzz(j)+dlj*zzt(j)
            
          enddo
          
        endif
        
      enddo

c     error exit for non-convergence

      if(.not.safe)return
      
c     complete stress tensor
      
      stresh(1)=strs1
      stresh(2)=strs2
      stresh(3)=strs3
      stresh(4)=strs2
      stresh(5)=strs5
      stresh(6)=strs6
      stresh(7)=strs3
      stresh(8)=strs6
      stresh(9)=strs9
      
c     splice coordinate arrays across nodes

      if(mxnode.gt.1)then

        buffer(1)=vircon
        call gdsum(buffer(1),1,buffer(2))
        vircon=buffer(1)
        call gdsum(stresh,9,buffer)
        call splice 
     x    (idnode,natms,listme,listot,xxx,yyy,zzz,buffer)
        call splice 
     x    (idnode,natms,listme,listot,vxx,vyy,vzz,buffer)

      endif
      
      return
      end subroutine rdrattle_r

      subroutine rdrattle_v
     x  (safe,idnode,mxnode,natms,nscons,tolnce,tstep,
     x  dxx,dyy,dzz,txx,tyy,tzz,xxt,yyt,zzt)

c*********************************************************************
c     
c     dl_poly subroutine for applying a bond constraints to 
c     the velocities of the constrained atoms using the rattle
c     procedure (replicated data version)
c     
c     copyright - daresbury laboratory
c     author w.smith october 2002
c     
c*********************************************************************

      implicit none

      logical safe
      integer idnode,mxnode,natms,nscons,icyc
      integer i,j,k
      real(8) tolnce,tstep,dxx,dyy,dzz,txx,tyy,tzz,xxt,yyt,zzt
      real(8) esig,amti,amtj,gamma,gammi,gammj,dli,dlj,tolvel

      dimension txx(mxatms),tyy(mxatms),tzz(mxatms)
      dimension xxt(mxatms),yyt(mxatms),zzt(mxatms)
      dimension dxx(mxcons),dyy(mxcons),dzz(mxcons)
      
c     constraint convergence tolerance

      tolvel=tolnce/tstep

c     start of rattle cycle
      
      icyc=0
      safe=.false.

      do while(.not.safe.and.icyc.lt.mxshak)

        icyc=icyc+1

c     initialise velocity correction arrays
        
        do i=1,natms
          
          xxt(i)=0.d0
          yyt(i)=0.d0
          zzt(i)=0.d0

        enddo

c     calculate velocity constraint corrections
        
        esig=0.d0

        do k=1,nscons
          
          i=listcon(k,2)
          j=listcon(k,3)
          
          amti= 0.5d0*tstep*rmass(i)
          amtj=-0.5d0*tstep*rmass(j)

c     calculate constraint force parameter

          gamma=(dxx(k)*(vxx(i)-vxx(j))+dyy(k)*(vyy(i)-vyy(j))+
     x      dzz(k)*(vzz(i)-vzz(j)))/((amti-amtj)*
     x      (dxx(k)**2+dyy(k)**2+dzz(k)**2))
          esig=max(esig,abs(gamma))

c     improve approximate constraint force

          gammi=-gamma*amti
          xxt(i)=xxt(i)+gammi*dxx(k)
          yyt(i)=yyt(i)+gammi*dyy(k)
          zzt(i)=zzt(i)+gammi*dzz(k)
          gammj=-gamma*amtj
          xxt(j)=xxt(j)+gammj*dxx(k)
          yyt(j)=yyt(j)+gammj*dyy(k)
          zzt(j)=zzt(j)+gammj*dzz(k)
          
        enddo

c     global verification of convergence
        
        safe=(esig.lt.tolvel)
        if(mxnode.gt.1)then
          
          call gstate(safe)
          
        endif

c     terminate iteration if constraints satisfied
        
        if(.not.safe)then
          
c     transport velocity adjustments to other nodes
          
          if(mxnode.gt.1)then
            
            call shmove
     x        (idnode,mxnode,natms,lashap,lishap,xxt,yyt,zzt,
     x        txx,tyy,tzz,buffer)
            
          endif
          
c     update velocities
          
          do k=1,nscons
            
            i=listcon(k,2)
            j=listcon(k,3)
            
            dli=1.d0/dble(listme(i))
            vxx(i)=vxx(i)+dli*xxt(i)
            vyy(i)=vyy(i)+dli*yyt(i)
            vzz(i)=vzz(i)+dli*zzt(i)
            dlj=1.d0/dble(listme(j))
            vxx(j)=vxx(j)+dlj*xxt(j)
            vyy(j)=vyy(j)+dlj*yyt(j)
            vzz(j)=vzz(j)+dlj*zzt(j)
            
          enddo
        
        endif
        
      enddo

c     error exit if rattle fails

      if(.not.safe)return
      
c     splice velocity arrays across nodes
      
      if(mxnode.gt.1)then

        call splice
     x    (idnode,natms,listme,listot,vxx,vyy,vzz,buffer)

      endif
      
      return
      end subroutine rdrattle_v

      subroutine nvevv_1
     x  (safe,lshmov,lmsite,isw,idnode,mxnode,natms,imcon,
     x  nscons,ntcons,ntpmls,tstep,engke,tolnce,vircon,g_qt4f)

c***********************************************************************
c     
c     dl_poly subroutine for integrating newtonian equations of
c     motion in molecular dynamics - velocity verlet incorporating
c     bond constraints via the shake/rattle algorithm
c     
c     nve ensemble
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 
c     author    - w. smith october 2002
c     amended   - w. smith january 2005 : f90 conversion
c
c     updated to account for qtip4p/f water model
c     Velocity/Verlet algorithm rewritten in simplectic form
c
c     copyright - M.R.Momeni and F.A.Shakib
c     authors   - M.R.Momeni and F.A.Shakib 2021
c
c     Method Development and Materials Simulation Laboratory
c     
c***********************************************************************

      implicit none
      
      integer, parameter :: nnn=4

      logical safe,lshmov,lmsite
      integer isw,idnode,mxnode,natms,imcon,nscons,ntcons,i,j,k
      integer iatm0,iatm1,ntpmls
      real(8) tstep,engke,tolnce,vircon
      integer fail(nnn)
      real(8) strkin(9)
      real(8) g_qt4f

      real(8), allocatable :: xxt(:),yyt(:),zzt(:)
      real(8), allocatable :: txx(:),tyy(:),tzz(:)
      real(8), allocatable :: dxx(:),dyy(:),dzz(:)
      real(8), allocatable :: dxt(:),dyt(:),dzt(:)
      
      safe=.true.

c     allocate working arrays

      if(ntcons.gt.0)then
        
        do i=1,nnn
          fail(i)=0
        enddo
        
        allocate(xxt(mxatms),yyt(mxatms),zzt(mxatms),stat=fail(1))
        allocate(txx(mxatms),tyy(mxatms),tzz(mxatms),stat=fail(2))
        allocate(dxx(mxcons),dyy(mxcons),dzz(mxcons),stat=fail(3))
        allocate(dxt(mxcons),dyt(mxcons),dzt(mxcons),stat=fail(4))
        
        do i=1,nnn
          if(fail(i).gt.0)call error(idnode,1980)
        enddo

      endif

c     block indices
      
      iatm0=(idnode*natms)/mxnode+1
      iatm1=((idnode+1)*natms)/mxnode

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
      
      if(ntcons.gt.0)
     x  call images(imcon,0,1,nscons,cell,dxx,dyy,dzz)

c     first stage of velocity verlet algorithm

      if(isw.eq.1)then

        do i=iatm0,iatm1

          vxx(i)=vxx(i)+0.5d0*tstep*rmass(i)*fxx(i)
          vyy(i)=vyy(i)+0.5d0*tstep*rmass(i)*fyy(i)
          vzz(i)=vzz(i)+0.5d0*tstep*rmass(i)*fzz(i)

        enddo
   
c     update positions
        
        do i=iatm0,iatm1
          
          xxx(i)=xxx(i)+tstep*vxx(i)
          yyy(i)=yyy(i)+tstep*vyy(i)
          zzz(i)=zzz(i)+tstep*vzz(i)
          
        enddo

c  update the position of M-site if water model qtip4p/f requested

!        if(lmsite)then

!          call qtip4pf(idnode,mxnode,imcon,nbeads,ntpmls,g_qt4f)

!        endif

c     merge position data
        
        if(mxnode.gt.1)
     x    call merge(idnode,mxnode,natms,mxbuff,xxx,yyy,zzz,buffer)
        
c  update the position of M-site if water model qtip4p/f requested

        if(lmsite)then

          call qtip4pf(idnode,mxnode,imcon,nbeads,ntpmls,g_qt4f)

        endif

c     apply shake corrections to bond constraints

        if(ntcons.gt.0)then
          
          if(mxnode.gt.1)
     x      call merge(idnode,mxnode,natms,mxbuff,vxx,vyy,vzz,buffer)

          safe=.false.
          call rdrattle_r
     x      (safe,lshmov,idnode,imcon,mxnode,natms,nscons,
     x      tolnce,tstep,vircon,dxx,dyy,dzz,dxt,dyt,dzt,
     x      txx,tyy,tzz,xxt,yyt,zzt,strcns)
          if(.not.safe)return

        endif

c     second stage of velocity verlet algorithm
        
      else
        
c  Redistribute the M-site forces if water model qtip4p/f requested

        if(lmsite)then

          call qt4_force_redist(idnode,mxnode,nbeads,ntpmls,g_qt4f)

        endif

        do i=iatm0,iatm1

          vxx(i)=vxx(i)+0.5d0*tstep*rmass(i)*fxx(i)
          vyy(i)=vyy(i)+0.5d0*tstep*rmass(i)*fyy(i)
          vzz(i)=vzz(i)+0.5d0*tstep*rmass(i)*fzz(i)

        enddo

c     merge velocity data

        if(mxnode.gt.1)
     x    call merge(idnode,mxnode,natms,mxbuff,vxx,vyy,vzz,buffer)

c     correct constraint bond velocities using rattle

        if(ntcons.gt.0)then

          safe=.false.
          call rdrattle_v
     x      (safe,idnode,mxnode,natms,nscons,tolnce,tstep,
     x      dxx,dyy,dzz,txx,tyy,tzz,xxt,yyt,zzt)
          if(.not.safe)return

        endif

c     calculate kinetic energy
        
        engke=getkin(natms,idnode,mxnode)

c     kinetic contribution to stress tensor
        
        call kinstress(natms,idnode,mxnode,strkin)

c     add contributions to stress tensor
        
        do i=1,9
          stress(i)=stress(i)+strcns(i)+strkin(i)
        enddo

      endif
      
c     periodic boundary condition
      
      if(isw.eq.2)then
        
        call images(imcon,idnode,mxnode,natms,cell,xxx,yyy,zzz)
        if(mxnode.gt.1)
     x    call merge(idnode,mxnode,natms,mxbuff,xxx,yyy,zzz,buffer)
        
      endif
      
c     deallocate working arrays

      if(ntcons.gt.0)then

        deallocate(xxt,yyt,zzt,txx,tyy,tzz,stat=fail(1))
        deallocate(dxx,dyy,dzz,dxt,dyt,dzt,stat=fail(2))

      endif
      
      return
      end subroutine nvevv_1

      subroutine nvtvv_b1
     x  (safe,lshmov,isw,idnode,mxnode,natms,imcon,nscons,ntcons,
     x  tstep,taut,sigma,engke,tolnce,vircon)

c***********************************************************************
c     
c     dl_poly subroutine for integrating newtonian equations of
c     motion in molecular dynamics - velocity verlet incorporating
c     bond constraints via the shake/rattle algorithm
c     
c     nvt ensemble - Berendsen thermostat (n.b. not symplectic)
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 
c     author    - w. smith october 2002
c     amended   - w. smith january 2005 : f90 conversion
c     
c***********************************************************************

      implicit none
      
      integer, parameter :: nnn=4

      logical safe,lshmov
      integer isw,idnode,mxnode,natms,imcon,nscons,ntcons
      integer i,j,k,iatm0,iatm1
      real(8) tstep,taut,sigma,engke,tolnce,vircon,chit
      integer fail(nnn)
      real(8) strkin(9)
      
      real(8), allocatable :: xxt(:),yyt(:),zzt(:)
      real(8), allocatable :: txx(:),tyy(:),tzz(:)
      real(8), allocatable :: dxx(:),dyy(:),dzz(:)
      real(8), allocatable :: dxt(:),dyt(:),dzt(:)
      
      safe=.true.

c     allocate working arrays

      if(ntcons.gt.0)then

        do i=1,nnn
          fail(i)=0
        enddo
        
        allocate(xxt(mxatms),yyt(mxatms),zzt(mxatms),stat=fail(1))
        allocate(txx(mxatms),tyy(mxatms),tzz(mxatms),stat=fail(2))
        allocate(dxx(mxcons),dyy(mxcons),dzz(mxcons),stat=fail(3))
        allocate(dxt(mxcons),dyt(mxcons),dzt(mxcons),stat=fail(4))
        
        do i=1,nnn
          if(fail(i).gt.0)call error(idnode,1990)
        enddo

      endif

c     block indices
      
      iatm0=(idnode*natms)/mxnode+1
      iatm1=((idnode+1)*natms)/mxnode

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
      
      if(ntcons.gt.0)
     x  call images(imcon,0,1,nscons,cell,dxx,dyy,dzz)

c     update velocities
      
      do i=iatm0,iatm1
        
        vxx(i)=vxx(i)+0.5d0*tstep*rmass(i)*fxx(i)
        vyy(i)=vyy(i)+0.5d0*tstep*rmass(i)*fyy(i)
        vzz(i)=vzz(i)+0.5d0*tstep*rmass(i)*fzz(i)
        
      enddo

c     first pass of velocity verlet algorithm
      
      if(isw.eq.1)then

c     update positions
        
        do i=iatm0,iatm1
          
          xxx(i)=xxx(i)+tstep*vxx(i)
          yyy(i)=yyy(i)+tstep*vyy(i)
          zzz(i)=zzz(i)+tstep*vzz(i)
          
        enddo

c     merge position data
        
        if(mxnode.gt.1)
     x    call merge(idnode,mxnode,natms,mxbuff,xxx,yyy,zzz,buffer)

c     apply shake corrections to bond constraints
        
        if(ntcons.gt.0)then
          
          if(mxnode.gt.1)
     x      call merge(idnode,mxnode,natms,mxbuff,vxx,vyy,vzz,buffer)

          safe=.false.
          call rdrattle_r
     x      (safe,lshmov,idnode,imcon,mxnode,natms,nscons,
     x      tolnce,tstep,vircon,dxx,dyy,dzz,dxt,dyt,dzt,
     x      txx,tyy,tzz,xxt,yyt,zzt,strcns)
          if(.not.safe)return

        endif

c     second pass of velocity verlet algorithm
        
      else

c     merge velocity data
        
        if(mxnode.gt.1)
     x    call merge(idnode,mxnode,natms,mxbuff,vxx,vyy,vzz,buffer)

c     correct constraint bond velocities using rattle
        
        if(ntcons.gt.0)then
          
          safe=.false.
          call rdrattle_v
     x      (safe,idnode,mxnode,natms,nscons,tolnce,tstep,
     x      dxx,dyy,dzz,txx,tyy,tzz,xxt,yyt,zzt)
          if(.not.safe)return
          
        endif

c     calculate kinetic energy
        
        engke=getkin(natms,idnode,mxnode)

c     apply Berendsen thermostat - taut is the relaxation time
        
        chit=sqrt(1.d0+tstep/taut*(sigma/engke-1.d0))
        
        do i=iatm0,iatm1
          
          vxx(i)=chit*vxx(i)
          vyy(i)=chit*vyy(i)
          vzz(i)=chit*vzz(i)
          
        enddo
        
        engke=engke*chit**2

c     merge velocity data
        
        if(mxnode.gt.1)
     x    call merge(idnode,mxnode,natms,mxbuff,vxx,vyy,vzz,
     x    buffer)

c     kinetic contribution to stress tensor
        
        call kinstress(natms,idnode,mxnode,strkin)

c     add contributions to stress tensor
        
        do i=1,9
          stress(i)=stress(i)+strkin(i)+strcns(i)
        enddo

      endif

c     periodic boundary condition
      
      if(isw.eq.2)then
        
        call images(imcon,idnode,mxnode,natms,cell,xxx,yyy,zzz)
        if(mxnode.gt.1)
     x    call merge(idnode,mxnode,natms,mxbuff,xxx,yyy,zzz,buffer)
        
      endif
      
c     deallocate working arrays

      if(ntcons.gt.0)then

        deallocate(xxt,yyt,zzt,txx,tyy,tzz,stat=fail(1))
        deallocate(dxx,dyy,dzz,dxt,dyt,dzt,stat=fail(2))

      endif
      
      return
      end subroutine nvtvv_b1

      subroutine nvtvv_e1
     x  (safe,lshmov,isw,idnode,mxnode,natms,imcon,nscons,ntcons,
     x  tstep,engke,tolnce,vircon)

c***********************************************************************
c     
c     dl_poly subroutine for integrating newtonian equations of
c     motion in molecular dynamics - velocity verlet incorporating
c     bond constraints via the shake/rattle algorithm
c     
c     nvt ensemble - evans thermostat
c     Comp. Phys. reports 1, 299, (1984)
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 
c     author    - w. smith october 2002
c     amended   - w. smith january 2005 : f90 conversion
c     
c***********************************************************************

      implicit none

      integer, parameter :: nnn=4

      logical safe,lshmov
      integer isw,idnode,mxnode,natms,imcon,nscons,ntcons
      integer i,j,k,iatm0,iatm1
      real(8) tstep,engke,tolnce,vircon,vdotf,scale,chit
      integer fail(nnn)
      real(8) strkin(9)

      real(8), allocatable :: xxt(:),yyt(:),zzt(:)
      real(8), allocatable :: txx(:),tyy(:),tzz(:)
      real(8), allocatable :: dxx(:),dyy(:),dzz(:)
      real(8), allocatable :: dxt(:),dyt(:),dzt(:)
      
      safe=.true.

c     allocate working arrays

      if(ntcons.gt.0)then

        do i=1,nnn
          fail(i)=0
        enddo
        
        allocate(xxt(mxatms),yyt(mxatms),zzt(mxatms),stat=fail(1))
        allocate(txx(mxatms),tyy(mxatms),tzz(mxatms),stat=fail(2))
        allocate(dxx(mxcons),dyy(mxcons),dzz(mxcons),stat=fail(3))
        allocate(dxt(mxcons),dyt(mxcons),dzt(mxcons),stat=fail(4))
        
        do i=1,nnn
          if(fail(i).gt.0)call error(idnode,2000)
        enddo

      endif

      if(ntcons.eq.0)safe=.true.

c     block indices
      
      iatm0=(idnode*natms)/mxnode+1
      iatm1=((idnode+1)*natms)/mxnode

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
      
      if(ntcons.gt.0)
     x  call images(imcon,0,1,nscons,cell,dxx,dyy,dzz)

c     first stage of velocity verlet algorithm
      
      if(isw.eq.1)then

c     calculate kinetic energy and evans thermostat parameter
        
        engke=0.d0
        vdotf=0.d0
        do i=iatm0,iatm1
          
          engke=engke+weight(i)*(vxx(i)**2+vyy(i)**2+vzz(i)**2)
          vdotf=vdotf+vxx(i)*fxx(i)+vyy(i)*fyy(i)+vzz(i)*fzz(i)
          
        enddo
        
        if(mxnode.gt.1)then
          
          buffer(1)=engke
          buffer(2)=vdotf
          call gdsum(buffer(1),2,buffer(3))
          engke=buffer(1)
          vdotf=buffer(2)
          
        endif
        engke=0.5d0*engke
        chit=0.5d0*vdotf/engke

c     thermostat the velocities
        
        scale=(1.d0-0.5d0*tstep*chit)
        do i=iatm0,iatm1
          
          vxx(i)=scale*vxx(i)
          vyy(i)=scale*vyy(i)
          vzz(i)=scale*vzz(i)
          
        enddo

c     update velocities
        
        do i=iatm0,iatm1
          
          vxx(i)=vxx(i)+0.5d0*tstep*rmass(i)*fxx(i)
          vyy(i)=vyy(i)+0.5d0*tstep*rmass(i)*fyy(i)
          vzz(i)=vzz(i)+0.5d0*tstep*rmass(i)*fzz(i)
          
        enddo

c     update positions
        
        do i=iatm0,iatm1
          
          xxx(i)=xxx(i)+tstep*vxx(i)
          yyy(i)=yyy(i)+tstep*vyy(i)
          zzz(i)=zzz(i)+tstep*vzz(i)
          
        enddo

c     merge position data
        
        if(mxnode.gt.1)
     x    call merge(idnode,mxnode,natms,mxbuff,xxx,yyy,zzz,buffer)

c     apply shake corrections to bond constraints
        
        if(ntcons.gt.0)then
          
          if(mxnode.gt.1)
     x      call merge(idnode,mxnode,natms,mxbuff,vxx,vyy,vzz,buffer)

          safe=.false.
          call rdrattle_r
     x      (safe,lshmov,idnode,imcon,mxnode,natms,nscons,
     x      tolnce,tstep,vircon,dxx,dyy,dzz,dxt,dyt,dzt,
     x      txx,tyy,tzz,xxt,yyt,zzt,strcns)
          if(.not.safe)return

        endif

c     second stage of velocity verlet algorithm
        
      else


c     update velocities
        
        do i=iatm0,iatm1
          
          vxx(i)=vxx(i)+0.5d0*tstep*rmass(i)*fxx(i)
          vyy(i)=vyy(i)+0.5d0*tstep*rmass(i)*fyy(i)
          vzz(i)=vzz(i)+0.5d0*tstep*rmass(i)*fzz(i)
          
        enddo

c     merge velocity data
        
        if(mxnode.gt.1)
     x    call merge(idnode,mxnode,natms,mxbuff,vxx,vyy,vzz,buffer)

c     correct constraint bond velocities using rattle
        
        if(ntcons.gt.0)then
          
          safe=.false.
          call rdrattle_v
     x      (safe,idnode,mxnode,natms,nscons,tolnce,tstep,
     x      dxx,dyy,dzz,txx,tyy,tzz,xxt,yyt,zzt)
          if(.not.safe)return
          
        endif

c     calculate kinetic energy and evans thermostat parameter
        
        engke=0.d0
        vdotf=0.d0
        do i=iatm0,iatm1
          
          engke=engke+weight(i)*(vxx(i)**2+vyy(i)**2+vzz(i)**2)
          vdotf=vdotf+vxx(i)*fxx(i)+vyy(i)*fyy(i)+vzz(i)*fzz(i)
          
        enddo
        
        if(mxnode.gt.1)then
          
          buffer(1)=engke
          buffer(2)=vdotf
          call gdsum(buffer(1),2,buffer(3))
          engke=buffer(1)
          vdotf=buffer(2)
          
        endif
        engke=0.5d0*engke
        chit=0.5d0*vdotf/engke
        scale=(1.d0-0.5d0*tstep*chit)

c     scale velocities
        
        do i=iatm0,iatm1
          
          vxx(i)=scale*vxx(i)
          vyy(i)=scale*vyy(i)
          vzz(i)=scale*vzz(i)
          
        enddo
        
        engke=engke*scale**2

c     merge velocity data
        
        if(mxnode.gt.1)
     x    call merge(idnode,mxnode,natms,mxbuff,vxx,vyy,vzz,buffer)

c     kinetic contribution to stress tensor
        
        call kinstress(natms,idnode,mxnode,strkin)

c     add contributions to stress tensor
        
        do i=1,9
          stress(i)=stress(i)+strkin(i)+strcns(i)
        enddo

      endif

c     periodic boundary condition
      
      if(isw.eq.2)then
        
        call images(imcon,idnode,mxnode,natms,cell,xxx,yyy,zzz)
        if(mxnode.gt.1)
     x    call merge(idnode,mxnode,natms,mxbuff,xxx,yyy,zzz,buffer)
        
      endif
      
c     deallocate working arrays

      if(ntcons.gt.0)then

        deallocate(xxt,yyt,zzt,txx,tyy,tzz,stat=fail(1))
        deallocate(dxx,dyy,dzz,dxt,dyt,dzt,stat=fail(2))

      endif
      
      return
      end subroutine nvtvv_e1

      subroutine nvtvv_h1
     x  (safe,lshmov,isw,idnode,mxnode,natms,imcon,nscons,ntcons,
     x  ntshl,keyshl,tstep,taut,sigma,chit,consv,conint,engke,
     x  tolnce,vircon,chit_shl,sigma_shl,
     x  lmsite,g_qt4f,ntpmls)

c***********************************************************************
c     
c     dl_poly subroutine for integrating newtonian equations of
c     motion in molecular dynamics - velocity verlet incorporating
c     bond constraints via the shake/rattle algorithm
c     
c     nvt ensemble - nose-hoover thermostat
c     Molecular Physics 87 (1996) 1117
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 
c     author    - w. smith october 2002
c     amended   - w. smith january 2005 : f90 conversion
c     adapted   - d. quigley - metadynamics
c     
c***********************************************************************

      implicit none

      integer, parameter :: nnn=4
      
      logical safe,lshmov,lmsite
      integer isw,idnode,mxnode,natms,imcon,nscons,ntcons
      integer i,j,k,iatm0,iatm1
      real(8) tstep,taut,sigma,chit,consv,conint,engke,tolnce,vircon
      real(8) hstep,qmass
      integer fail(nnn)
      real(8) strkin(9)

      integer ntpmls
      real(8) g_qt4f

      real(8), allocatable :: xxt(:),yyt(:),zzt(:)
      real(8), allocatable :: txx(:),tyy(:),tzz(:)
      real(8), allocatable :: dxx(:),dyy(:),dzz(:)
      real(8), allocatable :: dxt(:),dyt(:),dzt(:)
      
c     metadynamics shell thermostat variables

      integer ntshl,keyshl
      real(8) sigma_shl

      logical,save :: lfirst=.true.
      real(8)      :: chit_shl  
      real(8),save :: qmass_shl
      real(8)      :: shlke

c     end metadynamics shell thermostat variables
      
      safe=.true.

c     allocate working arrays

      if(ntcons.gt.0)then

        do i=1,nnn
          fail(i)=0
        enddo
        
        allocate(xxt(mxatms),yyt(mxatms),zzt(mxatms),stat=fail(1))
        allocate(txx(mxatms),tyy(mxatms),tzz(mxatms),stat=fail(2))
        allocate(dxx(mxcons),dyy(mxcons),dzz(mxcons),stat=fail(3))
        allocate(dxt(mxcons),dyt(mxcons),dzt(mxcons),stat=fail(4))
        
        do i=1,nnn
          if(fail(i).gt.0)call error(idnode,2010)
        enddo

      endif

c     inertia parameter for Nose-Hoover thermostat
      
      hstep=0.5d0*tstep
      qmass=2.d0*sigma*taut**2

c     block indices
      
      iatm0=(idnode*natms)/mxnode+1
      iatm1=((idnode+1)*natms)/mxnode

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
      
      if(ntcons.gt.0)
     x  call images(imcon,0,1,nscons,cell,dxx,dyy,dzz)

c     first stage of velocity verlet algorithm
      
      if(isw.eq.1)then

c     integrate and apply nvt thermostat

        call nvtscale
     x    (idnode,mxnode,natms,engke,sigma,hstep,qmass,taut,
     x    chit,conint)

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
        
c     update velocities
        
        do i=iatm0,iatm1
 
          vxx(i)=vxx(i)+hstep*rmass(i)*fxx(i)
          vyy(i)=vyy(i)+hstep*rmass(i)*fyy(i)
          vzz(i)=vzz(i)+hstep*rmass(i)*fzz(i)

        enddo

c     update positions
        
        do i=iatm0,iatm1
          
          xxx(i)=xxx(i)+tstep*vxx(i)
          yyy(i)=yyy(i)+tstep*vyy(i)
          zzz(i)=zzz(i)+tstep*vzz(i)
          
        enddo

c     merge position data
        
        if(mxnode.gt.1)
     x    call merge(idnode,mxnode,natms,mxbuff,xxx,yyy,zzz,buffer)

c     update the position of M-site if water model qtip4p/f requested

        if(lmsite)then

          call qtip4pf(idnode,mxnode,imcon,nbeads,ntpmls,g_qt4f)

        endif

c     apply shake corrections to bond constraints
        
        if(ntcons.gt.0)then
          
          if(mxnode.gt.1)
     x      call merge(idnode,mxnode,natms,mxbuff,vxx,vyy,vzz,buffer)

          safe=.false.
          call rdrattle_r
     x      (safe,lshmov,idnode,imcon,mxnode,natms,nscons,
     x      tolnce,tstep,vircon,dxx,dyy,dzz,dxt,dyt,dzt,
     x      txx,tyy,tzz,xxt,yyt,zzt,strcns)
          if(.not.safe)return

        endif

c     second stage of velocity verlet algorithm
        
      else

c  Redistribute the M-site forces if water model qtip4p/f requested

        if(lmsite)then

          call qt4_force_redist(idnode,mxnode,nbeads,ntpmls,g_qt4f)

        endif

c     update velocities
        
        do i=iatm0,iatm1
 
          vxx(i)=vxx(i)+hstep*rmass(i)*fxx(i)
          vyy(i)=vyy(i)+hstep*rmass(i)*fyy(i)
          vzz(i)=vzz(i)+hstep*rmass(i)*fzz(i)
          
        enddo
        
        if(ntcons.gt.0)then

c     merge velocity data
          
          if(mxnode.gt.1)
     x      call merge(idnode,mxnode,natms,mxbuff,vxx,vyy,vzz,buffer)

c     correct constraint bond velocities using rattle
          
          safe=.false.
          call rdrattle_v
     x      (safe,idnode,mxnode,natms,nscons,tolnce,tstep,
     x      dxx,dyy,dzz,txx,tyy,tzz,xxt,yyt,zzt)
          if(.not.safe)return
          
        endif

c     integrate and apply nvt thermostat
        
        call nvtscale
     x    (idnode,mxnode,natms,engke,sigma,hstep,qmass,taut,
     x    chit,conint)
        
c     merge velocity data
        
        if(mxnode.gt.1)
     x    call merge(idnode,mxnode,natms,mxbuff,vxx,vyy,vzz,buffer)

c     metdynamics shell thermostat
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
        
        call kinstress(natms,idnode,mxnode,strkin)

c     add contributions to stress tensor
        
        do i=1,9
          stress(i)=stress(i)+strkin(i)+strcns(i)
        enddo

      endif

c     periodic boundary condition
      
      if(isw.eq.2)then
        
        call images(imcon,idnode,mxnode,natms,cell,xxx,yyy,zzz)
        if(mxnode.gt.1)
     x    call merge(idnode,mxnode,natms,mxbuff,xxx,yyy,zzz,buffer)
        
      endif
      
c     deallocate working arrays

      if(ntcons.gt.0)then

        deallocate(xxt,yyt,zzt,txx,tyy,tzz,stat=fail(1))
        deallocate(dxx,dyy,dzz,dxt,dyt,dzt,stat=fail(2))

      endif
      
      return
      end subroutine nvtvv_h1

      subroutine nvtvv_nhc
     x      (safe,lshmov,lmsite,isw,idnode,mxnode,natms,imcon,nscons,
     x      ntcons,nchain,nrespa,ntpmls,         
     x      ntshl,keyshl,tstep,taut,sigma,sigma_nhc,g_qt4f,chit,
     x      consv,conint,engke,tolnce,vircon,chit_shl,sigma_shl)

c***********************************************************************
c     
c     dl_poly quantum subroutine for integrating newtonian equations of
c     motion in molecular dynamics - velocity verlet 
c     
c     nvt ensemble - nose-hoover-chain (NHC) thermostat
c     with Suzuki-Yoshida scheme and RESPA algorithm
c     
c     copyright - M.R.Momeni and F.A.Shakib
c     author    - M.R.Momeni and F.A.Shakib 2021
c
c     Method Development and Materials Simulation Laboratory
c
c***********************************************************************

      implicit none

      integer, parameter :: nnn=4
      
      logical safe,lmsite,lshmov
      integer isw,idnode,mxnode,natms,imcon,nscons,ntcons
      integer i,j,k,iatm0,iatm1
      integer nchain,nrespa,ntpmls
      real(8) tstep,taut,sigma,chit,consv,conint,engke,tolnce,vircon
      real(8) hstep,qmass_t,heta_nhc,hpeta_nhc
      real(8) heta_1,heta_rest,hpeta_1,hpeta_rest
      real(8) sigma_nhc,qmass_part,g_qt4f
      integer fail(nnn)
      real(8) strkin(9)

      real(8), allocatable :: xxt(:),yyt(:),zzt(:)
      real(8), allocatable :: txx(:),tyy(:),tzz(:)
      real(8), allocatable :: dxx(:),dyy(:),dzz(:)
      real(8), allocatable :: dxt(:),dyt(:),dzt(:)

c     metadynamics shell thermostat variables

      integer ntshl,keyshl
      real(8) sigma_shl

      logical,save :: lfirst=.true.
      real(8)      :: chit_shl  
      real(8),save :: qmass_shl
      real(8)      :: shlke

c     end metadynamics shell thermostat variables
      
      safe=.true.

c cccccccc shake/rattle testing cccccccc

c     allocate working arrays

      if(ntcons.gt.0)then

        do i=1,nnn
          fail(i)=0
        enddo
        
        allocate(xxt(mxatms),yyt(mxatms),zzt(mxatms),stat=fail(1))
        allocate(txx(mxatms),tyy(mxatms),tzz(mxatms),stat=fail(2))
        allocate(dxx(mxcons),dyy(mxcons),dzz(mxcons),stat=fail(3))
        allocate(dxt(mxcons),dyt(mxcons),dzt(mxcons),stat=fail(4))
        
        do i=1,nnn
          if(fail(i).gt.0)call error(idnode,2010)
        enddo

      endif

c cccccccccccccccccccccccccccccccccccccc

c     inertia parameter for Nose-Hoover Chain thermostat
      
      hstep=0.5d0*tstep
      qmass_t=2.d0*sigma*taut**2
      qmass_part=2.d0*sigma_nhc*taut**2

c     block indices
      
      iatm0=(idnode*natms)/mxnode+1
      iatm1=((idnode+1)*natms)/mxnode

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

c cccccccc shake/rattle testing cccccccc

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

      if(ntcons.gt.0)
     x  call images(imcon,0,1,nscons,cell,dxx,dyy,dzz)

c cccccccccccccccccccccccccccccccccccccc


c     first stage of velocity verlet algorithm
      
      if(isw.eq.1)then

c     integrate and apply nvt thermostat

        call nhc_part
     x     (idnode,mxnode,natms,nchain,nrespa,engke,sigma,sigma_nhc,
     x     hstep,qmass_t,qmass_part,taut)

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
        
c     update velocities
        
        do i=iatm0,iatm1
 
          vxx(i)=vxx(i)+hstep*rmass(i)*fxx(i)
          vyy(i)=vyy(i)+hstep*rmass(i)*fyy(i)
          vzz(i)=vzz(i)+hstep*rmass(i)*fzz(i)

        enddo

c     update positions
        
        do i=iatm0,iatm1
          
          xxx(i)=xxx(i)+tstep*vxx(i)
          yyy(i)=yyy(i)+tstep*vyy(i)
          zzz(i)=zzz(i)+tstep*vzz(i)
          
        enddo

c     merge position data, i.e. merge coordinate arrays across
c     a number of processors        

        if(mxnode.gt.1)
     x    call merge(idnode,mxnode,natms,mxbuff,xxx,yyy,zzz,buffer)

c     update the position of M-site if water model qtip4p/f requested

        if(lmsite)then

          call qtip4pf(idnode,mxnode,imcon,nbeads,ntpmls,g_qt4f)

        endif

c cccccccc shake/rattle testing cccccccc

c     apply shake corrections to bond constraints
        
        if(ntcons.gt.0)then
          
          if(mxnode.gt.1)
     x      call merge(idnode,mxnode,natms,mxbuff,vxx,vyy,vzz,buffer)

          safe=.false.
          call rdrattle_r
     x      (safe,lshmov,idnode,imcon,mxnode,natms,nscons,
     x      tolnce,tstep,vircon,dxx,dyy,dzz,dxt,dyt,dzt,
     x      txx,tyy,tzz,xxt,yyt,zzt,strcns)
          if(.not.safe)return

        endif

c cccccccccccccccccccccccccccccccccccccc

c     second stage of velocity verlet algorithm
        
      else

c  Redistribute the M-site forces if water model qtip4p/f requested

        if(lmsite)then

          call qt4_force_redist(idnode,mxnode,nbeads,ntpmls,g_qt4f)

        endif

c     update velocities
        
        do i=iatm0,iatm1
 
          vxx(i)=vxx(i)+hstep*rmass(i)*fxx(i)
          vyy(i)=vyy(i)+hstep*rmass(i)*fyy(i)
          vzz(i)=vzz(i)+hstep*rmass(i)*fzz(i)
          
        enddo
        
c cccccccc shake/rattle testing cccccccc

      if(ntcons.gt.0)then

c     merge velocity data
          
          if(mxnode.gt.1)
     x      call merge(idnode,mxnode,natms,mxbuff,vxx,vyy,vzz,buffer)

c     correct constraint bond velocities using rattle
          
          safe=.false.
          call rdrattle_v
     x      (safe,idnode,mxnode,natms,nscons,tolnce,tstep,
     x      dxx,dyy,dzz,txx,tyy,tzz,xxt,yyt,zzt)
          if(.not.safe)return
          
        endif

c cccccccccccccccccccccccccccccccccccccc

c     integrate and apply nvt thermostat
        
        call nhc_part
     x     (idnode,mxnode,natms,nchain,nrespa,engke,sigma,sigma_nhc,
     x     hstep,qmass_t,qmass_part,taut)
        
c     merge velocity data
        
        if(mxnode.gt.1)
     x    call merge(idnode,mxnode,natms,mxbuff,vxx,vyy,vzz,buffer)

c     metdynamics shell thermostat
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

       heta_1=2.d0*sigma*eta_nhc(1)

       heta_rest=0.d0
         do j=2,nchain
            heta_rest=heta_rest+eta_nhc(j)
         enddo

       heta_nhc=(2.d0*sigma_nhc*heta_rest)+heta_1

       hpeta_1=0.5d0*(peta(1)**2/qmass_t)

       hpeta_rest=0.d0
         do j=2,nchain
            hpeta_rest=hpeta_rest+0.5d0*(peta(j)**2/qmass_part)
         enddo

       hpeta_nhc=hpeta_1+hpeta_rest

       consv=heta_nhc+hpeta_nhc

c     metadynamics shell thermostat

        if(lmetadyn.and.keyshl.eq.1)then
           consv=consv+0.5d0*qmass_shl*chit_shl**2
        endif

c     kinetic contribution to stress tensor
        
        call kinstress(natms,idnode,mxnode,strkin)

c     add contributions to stress tensor
        
        do i=1,9
          stress(i)=stress(i)+strkin(i)+strcns(i)
        enddo

      endif

c     periodic boundary condition
      
      if(isw.eq.2)then
        
        call images(imcon,idnode,mxnode,natms,cell,xxx,yyy,zzz)
        if(mxnode.gt.1)
     x    call merge(idnode,mxnode,natms,mxbuff,xxx,yyy,zzz,buffer)
        
      endif
      
c cccccccc shake/rattle testing cccccccc

c     deallocate working arrays

      if(ntcons.gt.0)then

        deallocate(xxt,yyt,zzt,txx,tyy,tzz,stat=fail(1))
        deallocate(dxx,dyy,dzz,dxt,dyt,dzt,stat=fail(2))

      endif

c cccccccccccccccccccccccccccccccccccccc

      return
      end subroutine nvtvv_nhc

      subroutine nptvv_b1
     x  (safe,lshmov,isw,idnode,mxnode,natms,imcon,nscons,
     x  ntcons,ntpatm,tstep,taut,taup,sigma,engke,press,elrc,
     x  virlrc,tolnce,virtot,vircon,volm)
      
c***********************************************************************
c     
c     dl_poly subroutine for integrating newtonian equations of
c     motion in molecular dynamics - velocity verlet incorporating
c     bond constraints via the shake/rattle algorithm
c     
c     npt ensemble - Berendsen thermostat and barostat 
c     (n.b. not symplectic)
c     
c     isothermal compressibility (beta) set to that of liquid water
c     = 0.007372 dlpoly units
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 
c     author    - w. smith december 2002
c     amended   - w. smith  january 2005 : f90 conversion
c     
c***********************************************************************
      
      implicit none

      integer, parameter :: nnn=7

      logical newjob,safe,lshmov
      integer isw,idnode,mxnode,natms,imcon,nscons,ntcons,ntpatm
      integer i,j,k,iatm0,iatm1,mxiter,iter,kk
      real(8) tstep,taut,taup,sigma,engke,press,elrc,virlrc,tolnce
      real(8) virtot,vircon,volm,volm0,elrc0,virlrc0,psyst
      real(8) chit,chip,scale,beta

      integer fail(nnn)
      real(8) strkin(9),uni(9)

      real(8), allocatable :: dens0(:)
      real(8), allocatable :: xxo(:),yyo(:),zzo(:)
      real(8), allocatable :: xxt(:),yyt(:),zzt(:)
      real(8), allocatable :: txx(:),tyy(:),tzz(:)
      real(8), allocatable :: dxx(:),dyy(:),dzz(:)
      real(8), allocatable :: dxt(:),dyt(:),dzt(:)
      real(8), allocatable :: vxo(:),vyo(:),vzo(:)
      
      save newjob,volm0,elrc0,virlrc0,iatm0,iatm1,dens0
      
      data uni/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/
      data newjob/.true./,beta/7.3728d-3/
      
      safe=.true.
      
      if(newjob)then

c     block indices
        
        iatm0=(idnode*natms)/mxnode+1
        iatm1=((idnode+1)*natms)/mxnode

c     allocate density storage array

        fail(1)=0
        allocate(dens0(mxatyp),stat=fail(1))
        if(fail(1).gt.0)call error(idnode,2020)

c     store initial values of volume and long range corrections
        
        volm0=volm
        elrc0=elrc
        virlrc0=virlrc
        do i=1,ntpatm
          dens0(i)=dens(i)
        enddo
        
        newjob=.false.
        
      endif

c     allocate working arrays

      if(ntcons.gt.0)then

        do i=1,nnn
          fail(i)=0
        enddo
        allocate(xxo(mxatms),yyo(mxatms),zzo(mxatms),stat=fail(2))
        allocate(xxt(mxatms),yyt(mxatms),zzt(mxatms),stat=fail(3))
        allocate(txx(mxatms),tyy(mxatms),tzz(mxatms),stat=fail(4))
        allocate(dxx(mxcons),dyy(mxcons),dzz(mxcons),stat=fail(5))
        allocate(dxt(mxcons),dyt(mxcons),dzt(mxcons),stat=fail(6))
        allocate(vxo(mxatms),vyo(mxatms),vzo(mxatms),stat=fail(7))
        do i=1,nnn
          if(fail(i).gt.0)call error(idnode,2030)
        enddo

      endif

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

c     first stage of velocity verlet algorithm
      
      if(isw.eq.1)then          

c     calculate kinetic energy
        
        engke=getkin(natms,idnode,mxnode)

c     update velocities
        
        do i=iatm0,iatm1
 
          vxx(i)=vxx(i)+0.5d0*tstep*rmass(i)*fxx(i)
          vyy(i)=vyy(i)+0.5d0*tstep*rmass(i)*fyy(i)
          vzz(i)=vzz(i)+0.5d0*tstep*rmass(i)*fzz(i)
          
        enddo
          
c     store current integration variables
        
        if(ntcons.gt.0)then

          do i=iatm0,iatm1
            
            xxo(i)=xxx(i)
            yyo(i)=yyy(i)
            zzo(i)=zzz(i)
            vxo(i)=vxx(i)
            vyo(i)=vyy(i)
            vzo(i)=vzz(i)
            
          enddo
        
        endif

c     iteration required if ntcons > 0 and isw = 1

        mxiter=1
        if(isw.eq.1.and.ntcons.gt.0)mxiter=2
        do iter=1,mxiter
          
          scale=1.d0
          
          if(iter.eq.mxiter)then

c     calculate system pressure
            
            psyst=(2.d0*engke-virtot-vircon)/(3.d0*volm)

c     apply Berendsen barostat taup is relaxation time
            
            chip=1.d0+beta*tstep*(psyst-press)/taup
            scale=chip**(1.d0/3.d0)
            volm=chip*volm
            
c     reset cell parameters for new volume
          
            do i=1,9
              cell(i)=scale*cell(i)
            enddo
          
          endif

c     update positions
          
          do i=iatm0,iatm1
            
            xxx(i)=scale*xxx(i)+tstep*vxx(i)
            yyy(i)=scale*yyy(i)+tstep*vyy(i)
            zzz(i)=scale*zzz(i)+tstep*vzz(i)
            
          enddo
          
c     merge position data
          
          if(mxnode.gt.1)
     x      call merge(idnode,mxnode,natms,mxbuff,xxx,yyy,zzz,buffer)
          
c     apply shake corrections to bond constraints
          
          if(ntcons.gt.0)then
            
            if(mxnode.gt.1)
     x        call merge(idnode,mxnode,natms,mxbuff,vxx,vyy,vzz,buffer)
            
            safe=.false.
            call rdrattle_r
     x        (safe,lshmov,idnode,imcon,mxnode,natms,nscons,
     x        tolnce,tstep,vircon,dxx,dyy,dzz,dxt,dyt,dzt,
     x        txx,tyy,tzz,xxt,yyt,zzt,strcns)
            if(.not.safe)return
            
          endif
          
          if(iter.lt.mxiter)then
            
            do i=iatm0,iatm1

              xxx(i)=xxo(i)
              yyy(i)=yyo(i)
              zzz(i)=zzo(i)
              vxx(i)=vxo(i)
              vyy(i)=vyo(i)
              vzz(i)=vzo(i)

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

c     second pass of velocity verlet algorithm
        
      else

c     update velocities
        
        do i=iatm0,iatm1
 
          vxx(i)=vxx(i)+0.5d0*tstep*rmass(i)*fxx(i)
          vyy(i)=vyy(i)+0.5d0*tstep*rmass(i)*fyy(i)
          vzz(i)=vzz(i)+0.5d0*tstep*rmass(i)*fzz(i)
          
        enddo

c     apply Berendsen thermostat taut is relaxation time
        
        engke=getkin(natms,idnode,mxnode)
        chit=sqrt(1.d0+tstep/taut*(sigma/engke-1.d0))
        
        do i=iatm0,iatm1
          
          vxx(i)=chit*vxx(i)
          vyy(i)=chit*vyy(i)
          vzz(i)=chit*vzz(i)
          
        enddo

c     merge velocity data
        
        if(mxnode.gt.1)
     x    call merge(idnode,mxnode,natms,mxbuff,vxx,vyy,vzz,buffer)
          
c     correct constraint bond velocities using rattle
        
        if(ntcons.gt.0)then

          safe=.false.
          call rdrattle_v
     x      (safe,idnode,mxnode,natms,nscons,tolnce,tstep,
     x      dxx,dyy,dzz,txx,tyy,tzz,xxt,yyt,zzt)
          if(.not.safe)return
          
        endif

c     kinetic contribution to stress tensor
        
        call kinstress(natms,idnode,mxnode,strkin)
        engke=0.5d0*(strkin(1)+strkin(5)+strkin(9))

c     add contributions to stress tensor
        
        do i=1,9
          stress(i)=stress(i)+strkin(i)+strcns(i)
        enddo

      endif

c     periodic boundary condition
      
      if(isw.eq.2)then
        
        call images(imcon,idnode,mxnode,natms,cell,xxx,yyy,zzz)
        if(mxnode.gt.1)
     x    call merge(idnode,mxnode,natms,mxbuff,xxx,yyy,zzz,buffer)
        
      endif
      
c     deallocate working arrays

      if(ntcons.gt.0)then

        deallocate(xxo,yyo,zzo,vxo,vyo,vzo,stat=fail(1))
        deallocate(txx,tyy,tzz,dxx,dyy,dzz,stat=fail(2))
        deallocate(dxt,dyt,dzt,xxt,yyt,zzt,stat=fail(3))

      endif
      
      return
      end subroutine nptvv_b1

      subroutine nptvv_h1
     x  (safe,lshmov,isw,idnode,mxnode,natms,imcon,nscons,ntcons,
     x  ntpatm,ntshl,keyshl,tstep,taut,taup,sigma,temp,chip,chit,
     x  consv,conint,engke,elrc,tolnce,vircon,virtot,virlrc,volm,
     x  press,chit_shl,sigma_shl)

c***********************************************************************
c     
c     dl_poly subroutine for integrating newtonian equations of
c     motion in molecular dynamics - velocity verlet incorporating
c     bond constraints via the shake/rattle algorithm
c     
c     npt ensemble - Melchionna, Ciccotti and Holian
c     Molecular Physics 78 (1993) 533
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 
c     author    - w. smith november 2002
c     amended   - w. smith january 2005: f90 conversion
c     adapted   - d. quigley - metadynamics
c     
c***********************************************************************

      implicit none

      integer, parameter :: nnn=7
      integer, parameter :: ncyc=5
      
      logical safe,lshmov,newjob
      integer isw,idnode,mxnode,natms,imcon,nscons,ntcons,ntpatm
      integer i,j,k,iatm0,iatm1,mxiter,iter,kk,icyc
      real(8) tstep,taup,taut,sigma,temp,chip,chit,consv,conint
      real(8) engke,elrc,tolnce,vircon,virtot,virlrc,volm,press
      real(8) volm0,elrc0,virlrc0,hstep,qstep,totmas,qmass,pmass
      real(8) vzero,chit0,chip0,cons0,scale,fstep

      real(8), allocatable :: dens0(:)
      real(8), allocatable :: xxo(:),yyo(:),zzo(:)
      real(8), allocatable :: xxt(:),yyt(:),zzt(:)
      real(8), allocatable :: txx(:),tyy(:),tzz(:)
      real(8), allocatable :: dxx(:),dyy(:),dzz(:)
      real(8), allocatable :: dxt(:),dyt(:),dzt(:)
      real(8), allocatable :: vxo(:),vyo(:),vzo(:)

      integer fail(nnn)
      real(8) cell0(9),com(3),vom(3),strkin(9),uni(9)

      save newjob,totmas,volm0,elrc0,virlrc0,dens0
      save cell0,iatm0,iatm1,hstep,qstep,fstep,pmass,qmass

c     metadynamics shell thermostat variables

      integer ntshl,keyshl
      real(8) sigma_shl

      logical,save :: lfirst=.true.
      real(8)      :: chit_shl 
      real(8),save :: qmass_shl
      real(8)      :: shlke

c     end metadynamics shell thermostat variables

      data uni/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/
      data newjob/.true./
      
      safe=.true.
      
      if(newjob)then

c     block indices
        
        iatm0=(idnode*natms)/mxnode+1
        iatm1=((idnode+1)*natms)/mxnode

c     allocate density storage array

        fail(1)=0
        allocate(dens0(mxatyp),stat=fail(1))
        if(fail(1).gt.0)call error(idnode,2040)

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

c     inertia parameter for Nose-Hoover thermostat and barostat
        
        qmass=2.d0*sigma*taut**2
        pmass=2.d0*sigma*taup**2
        
        newjob=.false.
        
      endif

c     allocate working arrays

      if(ntcons.gt.0)then

        do i=1,nnn
          fail(i)=0
        enddo
        allocate(xxo(mxatms),yyo(mxatms),zzo(mxatms),stat=fail(2))
        allocate(vxo(mxatms),vyo(mxatms),vzo(mxatms),stat=fail(3))
        allocate(xxt(mxatms),yyt(mxatms),zzt(mxatms),stat=fail(4))
        allocate(txx(mxatms),tyy(mxatms),tzz(mxatms),stat=fail(5))
        allocate(dxx(mxcons),dyy(mxcons),dzz(mxcons),stat=fail(6))
        allocate(dxt(mxcons),dyt(mxcons),dzt(mxcons),stat=fail(7))
        do i=1,nnn
          if(fail(i).gt.0)call error(idnode,2050)
        enddo
        
      endif
      
      if(lmetadyn.and.lfirst.and.(ntshl>0))then
        if(idnode.eq.0)then
          write(*,*)"Warning - Metdynamics Modification"
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
      
c     first stage of velocity verlet algorithm
      
      if(isw.eq.1)then
        
c     store current integration variables if ntcons > 0
        
        if(ntcons.gt.0)then
          
          vzero=volm
          chit0=chit
          chip0=chip
          cons0=conint
          do i=iatm0,iatm1
            
            xxo(i)=xxx(i)
            yyo(i)=yyy(i)
            zzo(i)=zzz(i)
            vxo(i)=vxx(i)
            vyo(i)=vyy(i)
            vzo(i)=vzz(i)
            
          enddo
          
        endif
        
c     iteration necessary if ntcons > 0 and isw = 1
        
        mxiter=1
        if(isw.eq.1.and.ntcons.gt.0)mxiter=2
        do iter=1,mxiter
          
c     volume integration parameter
          
          do icyc=1,ncyc
            
c     integrate and apply npt thermostat
            
            call nptscale_t
     x        (idnode,mxnode,natms,engke,temp,sigma,qstep,pmass,qmass,
     x        taut,chip,chit,conint)
            
c     metdynamics shell thermostat
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
            
            call nptscale_p
     x        (idnode,mxnode,natms,engke,fstep,pmass,chip,chit,
     x        volm,press,vircon,virtot)
            
c     integrate and apply npt thermostat
            
            call nptscale_t
     x        (idnode,mxnode,natms,engke,temp,sigma,qstep,pmass,qmass,
     x        taut,chip,chit,conint)
            
c     metdynamics shell thermostat
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
          
c     scale cell vectors - isotropic
          
          scale=(volm/volm0)**(1.d0/3.d0)
          do i=1,9
            cell(i)=cell0(i)*scale
          enddo
          
c     update velocities
          
          do i=iatm0,iatm1
            
            vxx(i)=vxx(i)+hstep*rmass(i)*fxx(i)
            vyy(i)=vyy(i)+hstep*rmass(i)*fyy(i)
            vzz(i)=vzz(i)+hstep*rmass(i)*fzz(i)
            
          enddo

c     calculate system centre of mass
          
          call getcom(natms,idnode,mxnode,totmas,com)

c     update positions
          
          scale=exp(tstep*chip)
          do i=iatm0,iatm1
            
            xxx(i)=scale*(xxx(i)-com(1))+tstep*vxx(i)+com(1)
            yyy(i)=scale*(yyy(i)-com(2))+tstep*vyy(i)+com(2)
            zzz(i)=scale*(zzz(i)-com(3))+tstep*vzz(i)+com(3)
            
          enddo

c     merge position data
          
          if(mxnode.gt.1)
     x      call merge(idnode,mxnode,natms,mxbuff,xxx,yyy,zzz,buffer)

c     apply shake corrections to bond constraints
          
          if(ntcons.gt.0)then
            
            if(mxnode.gt.1)
     x        call merge(idnode,mxnode,natms,mxbuff,vxx,vyy,vzz,buffer)
            
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

            do i=iatm0,iatm1
              
              xxx(i)=xxo(i)
              yyy(i)=yyo(i)
              zzz(i)=zzo(i)
              vxx(i)=vxo(i)
              vyy(i)=vyo(i)
              vzz(i)=vzo(i)
              
            enddo
            
          endif
          
        enddo

c     second stage of velocity verlet algorithm
        
      else

c     update velocities
        
        do i=iatm0,iatm1
          
          vxx(i)=vxx(i)+hstep*rmass(i)*fxx(i)
          vyy(i)=vyy(i)+hstep*rmass(i)*fyy(i)
          vzz(i)=vzz(i)+hstep*rmass(i)*fzz(i)
          
        enddo
        
        if(ntcons.gt.0)then

c     merge velocity data
          
          if(mxnode.gt.1)
     x      call merge(idnode,mxnode,natms,mxbuff,vxx,vyy,vzz,buffer)

c     correct constraint bond velocities using rattle
          
          safe=.false.
          call rdrattle_v
     x      (safe,idnode,mxnode,natms,nscons,tolnce,tstep,
     x      dxx,dyy,dzz,txx,tyy,tzz,xxt,yyt,zzt)
          if(.not.safe)return
          
        endif

        do icyc=1,ncyc
          
c     integrate and apply npt thermostat
          
          call nptscale_t
     x      (idnode,mxnode,natms,engke,temp,sigma,qstep,pmass,qmass,
     x      taut,chip,chit,conint)
          
c     metdynamics shell thermostat
c     ====================================================
c     Must first merge update velocities as the core-shell
c     velocites are not distributed according to the same
c     rules.
c     ====================================================
          if(lmetadyn.and.keyshl.eq.1)then
            if(mxnode.gt.1)call merge
     x        (idnode,mxnode,natms,mxbuff,vxx,vyy,vzz,buffer)
            call nvtscale_shl
     x        (idnode,mxnode,ntshl,shlke,sigma_shl,hstep,qmass_shl,
     x        taut,chit_shl,conint)      
          endif  
          
c     integrate and apply npt barostat
          
          call nptscale_p
     x      (idnode,mxnode,natms,engke,fstep,pmass,chip,chit,
     x      volm,press,vircon,virtot)
          
c     integrate and apply npt thermostat
          
          call nptscale_t
     x      (idnode,mxnode,natms,engke,temp,sigma,qstep,pmass,qmass,
     x      taut,chip,chit,conint)

c     metadynamics shell thermostat
c     ====================================================
c     Must first merge update velocities as the core-shell
c     velocites are not distributed according to the same
c     rules.
c     ====================================================
          if(lmetadyn.and.keyshl.eq.1)then
            if(mxnode.gt.1)call merge
     x        (idnode,mxnode,natms,mxbuff,vxx,vyy,vzz,buffer)
            call nvtscale_shl
     x        (idnode,mxnode,ntshl,shlke,sigma_shl,hstep,qmass_shl,
     x        taut,chit_shl,conint)      
          endif  
          
        enddo

c     remove system centre of mass velocity
        
        call getvom(natms,idnode,mxnode,totmas,vom)
        
        do i=iatm0,iatm1
          
          vxx(i)=vxx(i)-vom(1)
          vyy(i)=vyy(i)-vom(2)
          vzz(i)=vzz(i)-vom(3)
          
        enddo

c     merge velocity data
        
        if(mxnode.gt.1)
     x    call merge(idnode,mxnode,natms,mxbuff,vxx,vyy,vzz,buffer)

c     scale cell vectors - isotropic
        
        scale=(volm/volm0)**(1.d0/3.d0)
        do i=1,9
          cell(i)=cell0(i)*scale
        enddo

c     conserved quantity less kinetic and potential energy terms
        
        consv=conint+0.5d0*qmass*chit**2+
     x    0.5d0*pmass*chip**2+press*volm

c     kinetic contribution to stress tensor
        
        call kinstress(natms,idnode,mxnode,strkin)
        engke=0.5d0*(strkin(1)+strkin(5)+strkin(9))

c     metadynamics shell thermostat

        if(lmetadyn.and.keyshl.eq.1)then
          consv=consv+0.5d0*qmass_shl*chit_shl**2
        endif

c     add contributions to stress tensor
        
        do i=1,9
          stress(i)=stress(i)+strkin(i)+strcns(i)
        enddo

      endif

c     periodic boundary condition
      
      if(isw.eq.2)then
        
        call images(imcon,idnode,mxnode,natms,cell,xxx,yyy,zzz)
        if(mxnode.gt.1)
     x    call merge(idnode,mxnode,natms,mxbuff,xxx,yyy,zzz,buffer)
        
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

      if(ntcons.gt.0)then

        deallocate(xxo,yyo,zzo,vxo,vyo,vzo,stat=fail(1))
        deallocate(txx,tyy,tzz,dxx,dyy,dzz,stat=fail(2))
        deallocate(dxt,dyt,dzt,xxt,yyt,zzt,stat=fail(3))

      endif
      
      return
      end subroutine nptvv_h1

      subroutine nptvv_nhc
     x    (safe,lmsite,isw,idnode,mxnode,natms,imcon,
     x    nchain,nrespa,ntpmls,        
     x    ntpatm,ntshl,keyshl,tstep,taut,taup,sigma,
     x    sigma_nhc,sigma_volm,alpha_volm,virtot,vircon,virlrc,
     x    g_qt4f,press,volm,chit,consv,
     x    conint,engke,elrc,chit_shl,sigma_shl)

c***********************************************************************
c     
c     dl_poly quantum subroutine for integrating newtonian equations of
c     motion in molecular dynamics - velocity verlet 
c     
c     npt ensemble - nose-hoover-chain (NHC) thermostat/barostat
c     with isotropic cell fluctuations, Suzuki-Yoshida scheme, and
c     Martyna-Tobias-Klein algorithm
c     
c     copyright - M.R.Momeni and F.A.Shakib
c     author    - M.R.Momeni and F.A.Shakib 2021
c
c     Method Development and Materials Simulation Laboratory
c     
c***********************************************************************

      implicit none

      integer, parameter :: nnn=7
      
      logical safe,newjob,lmsite
      integer isw,idnode,mxnode,natms,imcon,ntpatm
      integer i,j,k,iatm0,iatm1,kk
      integer nchain,nrespa,ntpmls
      real(8) tstep,taut,taup,chit,consv,conint,engke,elrc
      real(8) hstep,sigma,qmass_t
      real(8) sigma_nhc,qmass_baro,qmass_part
      real(8) volm_mass,sigma_volm,alpha_volm,g_qt4f
      real(8) press,volm,virtot,vircon,virlrc
      real(8) heta_nhc,hpeta_nhc,heta_1,heta_rest,hpeta_1,hpeta_rest
      real(8) hepsilon,hksi,hpksi
      real(8) volm0,elrc0,virlrc0,scale
      real(8) a_2n(6),sinh_v,sinh_r

      real(8), allocatable :: dens0(:)

      integer fail(nnn)
      real(8) cell0(9),com(3),vom(3),strkin(9),uni(9)

      save newjob,volm0,elrc0,virlrc0,dens0
      save cell0,iatm0,iatm1,hstep
      save qmass_t,qmass_baro,qmass_part,volm_mass

c     metadynamics shell thermostat variables

      integer ntshl,keyshl
      real(8) sigma_shl

      logical,save :: lfirst=.true.
      real(8)      :: chit_shl  
      real(8),save :: qmass_shl
      real(8)      :: shlke

c     end metadynamics shell thermostat variables

      real(8) :: Pint

      data uni/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/
      data newjob/.true./

      safe=.true.

      if(newjob)then

c     block indices

        iatm0=(idnode*natms)/mxnode+1
        iatm1=((idnode+1)*natms)/mxnode

c     allocate density storage array

        fail(1)=0
        allocate(dens0(mxatyp),stat=fail(1))
        if(fail(1).gt.0)call error(idnode,2040)

c     store intitial parameters

        volm0=volm
        elrc0=elrc
        virlrc0=virlrc
        hstep=0.5d0*tstep

        do i=1,ntpatm
          dens0(i)=dens(i)
        enddo

        do i=1,9
          cell0(i)=cell(i)
        enddo

c     inertia parameter for Nose-Hoover thermostat and barostat

        qmass_t=2.d0*sigma*taut**2
        qmass_part=2.d0*sigma_nhc*taut**2
        qmass_baro=2.d0*sigma_nhc*taup**2
        volm_mass=sigma_volm*taup**2

        newjob=.false.

      endif
        
c     a_2n parameters for NHC barostat

        a_2n(1)=1.d0
        a_2n(2)=1.d0/6.d0
        a_2n(3)=1.d0/120.d0
        a_2n(4)=1.d0/5040.d0
        a_2n(5)=1.d0/362880.d0
        a_2n(6)=1.d0/39916800.d0

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

c     first stage of velocity verlet algorithm
      
      if(isw.eq.1)then

c     integrate and apply nhc barostat

        call nhc_baro
     x     (idnode,mxnode,natms,nchain,nrespa,sigma_nhc,
     x     hstep,volm_mass,qmass_baro,taup,v_epsilon)

c     integrate and apply nvt thermostat

        call nhc_part
     x     (idnode,mxnode,natms,nchain,nrespa,engke,sigma,sigma_nhc,
     x     hstep,qmass_t,qmass_part,taut)

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

c     update the momentum variable of the volume

        v_epsilon=v_epsilon+hstep*(1.d0/volm_mass)*
     x     ((alpha_volm*2.0*engke)-virtot-vircon-3.d0*(press*volm))   

c     first integrate sinh(alpha*v_epsilon*tstep/4)

        sinh_v=0.d0 
        do i=1,6
          sinh_v=sinh_v+a_2n(i)
     x          *(alpha_volm*v_epsilon*hstep*0.5d0)**(2*(i-1))
        enddo         
          
c     update velocities

        do i=iatm0,iatm1
 
          vxx(i)=vxx(i)*exp(-alpha_volm*v_epsilon*hstep)+hstep*rmass(i)
     x    *fxx(i)*exp(-alpha_volm*v_epsilon*hstep*0.5d0)*sinh_v
          vyy(i)=vyy(i)*exp(-alpha_volm*v_epsilon*hstep)+hstep*rmass(i)
     x    *fyy(i)*exp(-alpha_volm*v_epsilon*hstep*0.5d0)*sinh_v
          vzz(i)=vzz(i)*exp(-alpha_volm*v_epsilon*hstep)+hstep*rmass(i)
     x    *fzz(i)*exp(-alpha_volm*v_epsilon*hstep*0.5d0)*sinh_v

        enddo

c     first integrate sinh(v_epsilon*tstep/2)

        sinh_r=0.d0 
        do i=1,6
          sinh_r=sinh_r+a_2n(i)*(v_epsilon*hstep)**(2*(i-1))
        enddo         
          
c     update positions
        
        do i=iatm0,iatm1
          
          xxx(i)=xxx(i)*exp(v_epsilon*tstep)+tstep*vxx(i)
     x    *exp(v_epsilon*hstep)*sinh_r           
          yyy(i)=yyy(i)*exp(v_epsilon*tstep)+tstep*vyy(i)
     x    *exp(v_epsilon*hstep)*sinh_r           
          zzz(i)=zzz(i)*exp(v_epsilon*tstep)+tstep*vzz(i)
     x    *exp(v_epsilon*hstep)*sinh_r           
          
        enddo

c     update the volume based on v_epsilon

        volm=volm*exp(3.d0*tstep*v_epsilon)

c     scale cell vectors - isotropic

          scale=(volm/volm0)**(1.d0/3.d0)
c          scale=exp(tstep*v_epsilon)
          do i=1,9
            cell(i)=cell0(i)*scale
c            cell(i)=cell(i)*scale
          enddo

c     merge position data
        
        if(mxnode.gt.1)
     x    call merge(idnode,mxnode,natms,mxbuff,xxx,yyy,zzz,buffer)

c     update the position of M-site if water model qtip4p/f requested

        if(lmsite)then

          call qtip4pf(idnode,mxnode,imcon,nbeads,ntpmls,g_qt4f)

        endif

c     second stage of velocity verlet algorithm
        
      else

c  Redistribute the M-site forces if water model qtip4p/f requested

        if(lmsite)then

          call qt4_force_redist(idnode,mxnode,nbeads,ntpmls,g_qt4f)

        endif

c     update velocities
        
        do i=iatm0,iatm1
 
          vxx(i)=vxx(i)*exp(-alpha_volm*v_epsilon*hstep)+hstep*rmass(i)
     x    *fxx(i)*exp(-alpha_volm*v_epsilon*hstep*0.5d0)*sinh_v
          vyy(i)=vyy(i)*exp(-alpha_volm*v_epsilon*hstep)+hstep*rmass(i)
     x    *fyy(i)*exp(-alpha_volm*v_epsilon*hstep*0.5d0)*sinh_v
          vzz(i)=vzz(i)*exp(-alpha_volm*v_epsilon*hstep)+hstep*rmass(i)
     x    *fzz(i)*exp(-alpha_volm*v_epsilon*hstep*0.5d0)*sinh_v

        enddo

c     update the momentum variable of the volume

        v_epsilon=v_epsilon+hstep*(1.d0/volm_mass)*
     x     ((alpha_volm*2.0*engke)-virtot-vircon-3.d0*(press*volm))   

c     integrate and apply nvt thermostat
        
        call nhc_part
     x     (idnode,mxnode,natms,nchain,nrespa,engke,sigma,sigma_nhc,
     x     hstep,qmass_t,qmass_part,taut)
        
c     metdynamics shell thermostat
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
        
c     integrate and apply nhc barostat

        call nhc_baro
     x     (idnode,mxnode,natms,nchain,nrespa,sigma_nhc,
     x     hstep,volm_mass,qmass_baro,taup,v_epsilon)

c     merge velocity data
        
        if(mxnode.gt.1)
     x    call merge(idnode,mxnode,natms,mxbuff,vxx,vyy,vzz,buffer)

c     conserved quantity less kinetic and potential energy terms

       hepsilon=0.5d0*volm_mass*(v_epsilon**2)

       hksi=0.d0
         do j=1,nchain
            hksi=hksi+ksi(j)
         enddo

       hksi=2.d0*sigma_nhc*hksi

       hpksi=0.d0
         do j=1,nchain
            hpksi=hpksi+0.5d0*(pksi(j)**2/qmass_baro)
         enddo

       heta_1=2.d0*sigma*eta_nhc(1)

       heta_rest=0.d0
         do j=2,nchain
            heta_rest=heta_rest+eta_nhc(j)
         enddo

       heta_nhc=(2.d0*sigma_nhc*heta_rest)+heta_1

       hpeta_1=0.5d0*(peta(1)**2/qmass_t)

       hpeta_rest=0.d0
         do j=2,nchain
            hpeta_rest=hpeta_rest+0.5d0*(peta(j)**2/qmass_part)
         enddo

       hpeta_nhc=hpeta_1+hpeta_rest

       consv=hepsilon+hksi+hpksi+heta_nhc+hpeta_nhc+(press*volm)

c     kinetic contribution to stress tensor
        
        call kinstress(natms,idnode,mxnode,strkin)
        engke=0.5d0*(strkin(1)+strkin(5)+strkin(9))

c     metadynamics shell thermostat

        if(lmetadyn.and.keyshl.eq.1)then
           consv=consv+0.5d0*qmass_shl*chit_shl**2
        endif

c     add contributions to stress tensor
        
        do i=1,9
          stress(i)=stress(i)+strkin(i)+strcns(i)
        enddo

      endif

c     periodic boundary condition
      
      if(isw.eq.2)then
        
        call images(imcon,idnode,mxnode,natms,cell,xxx,yyy,zzz)
        if(mxnode.gt.1)
     x    call merge(idnode,mxnode,natms,mxbuff,xxx,yyy,zzz,buffer)
        
      endif

c     adjust long range corrections and number density

      elrc=elrc0*(volm0/volm)
      virlrc=virlrc0*(volm0/volm)
      do kk=1,ntpatm
        dens(kk)=dens0(kk)*(volm0/volm)
      enddo

      return
      end subroutine nptvv_nhc

      subroutine nstvv_b1
     x  (safe,lshmov,isw,idnode,mxnode,natms,imcon,nscons,
     x  ntcons,ntpatm,mode,tstep,taut,taup,sigma,engke,press,elrc,
     x  virlrc,tolnce,vircon,volm)
      
c***********************************************************************
c     
c     dl_poly subroutine for integrating newtonian equations of
c     motion in molecular dynamics - velocity verlet incorporating
c     bond constraints via the shake/rattle algorithm
c     
c     anisotropic npt ensemble - Berendsen thermostat and barostat 
c     (n.b. not symplectic)
c     
c     isothermal compressibility (beta) set to that of liquid water
c     = 0.007372 dlpoly units
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 
c     author    - w. smith december 2002
c     amended   - w. smith january 2005 : f90 conversion
c     
c***********************************************************************
      
      implicit none

      integer, parameter :: nnn=7

      logical newjob,safe,lshmov
      integer isw,idnode,mxnode,natms,imcon,nscons,ntcons,ntpatm
      integer i,j,k,iatm0,iatm1,mxiter,iter,kk,mode
      real(8) tstep,taut,taup,sigma,engke,press,elrc,virlrc,tolnce,beta
      real(8) vircon,volm,volm0,elrc0,virlrc0,chit
      real(8) xtmp,ytmp,ztmp

      real(8), allocatable :: dens0(:)
      real(8), allocatable :: xxo(:),yyo(:),zzo(:)
      real(8), allocatable :: xxt(:),yyt(:),zzt(:)
      real(8), allocatable :: txx(:),tyy(:),tzz(:)
      real(8), allocatable :: dxx(:),dyy(:),dzz(:)
      real(8), allocatable :: dxt(:),dyt(:),dzt(:)
      real(8), allocatable :: vxo(:),vyo(:),vzo(:)

      integer fail(nnn)
      real(8) uni(9),strkin(9),celp(10)

      save newjob,volm0,elrc0,virlrc0,iatm0,iatm1,dens0
      
      data uni/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/
      data newjob/.true./,beta/7.3728d-3/
      
      safe=.true.
      
      if(newjob)then

c     block indices
        
        iatm0=(idnode*natms)/mxnode+1
        iatm1=((idnode+1)*natms)/mxnode

c     allocate density storage array

        fail(1)=0
        allocate(dens0(mxatyp),stat=fail(1))
        if(fail(1).gt.0)call error(idnode,2060)

c     store initial values of volume and long range corrections
        
        volm0=volm
        elrc0=elrc
        virlrc0=virlrc
        do i=1,ntpatm
          dens0(i)=dens(i)
        enddo
        
        newjob=.false.
        
      endif

c     allocate working arrays

      if(ntcons.gt.0)then

        do i=1,nnn
          fail(i)=0
        enddo
        allocate(xxo(mxatms),yyo(mxatms),zzo(mxatms),stat=fail(2))
        allocate(xxt(mxatms),yyt(mxatms),zzt(mxatms),stat=fail(3))
        allocate(txx(mxatms),tyy(mxatms),tzz(mxatms),stat=fail(4))
        allocate(dxx(mxcons),dyy(mxcons),dzz(mxcons),stat=fail(5))
        allocate(dxt(mxcons),dyt(mxcons),dzt(mxcons),stat=fail(6))
        allocate(vxo(mxatms),vyo(mxatms),vzo(mxatms),stat=fail(7))
        do i=1,nnn
          if(fail(i).gt.0)call error(idnode,2070)
        enddo
        
      endif

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

c     first stage of velocity verlet algorithm
      
      if(isw.eq.1)then          

c     extract previous constraint terms from stress tensor

        do i=1,9
          stress(i)=stress(i)-strcns(i)
        enddo

c     update velocities
        
        do i=iatm0,iatm1
          
          vxx(i)=vxx(i)+0.5d0*tstep*rmass(i)*fxx(i)
          vyy(i)=vyy(i)+0.5d0*tstep*rmass(i)*fyy(i)
          vzz(i)=vzz(i)+0.5d0*tstep*rmass(i)*fzz(i)
          
        enddo

c     store current integration variables
        
        if(ntcons.gt.0)then
          
          do i=iatm0,iatm1
            
            xxo(i)=xxx(i)
            yyo(i)=yyy(i)
            zzo(i)=zzz(i)
            vxo(i)=vxx(i)
            vyo(i)=vyy(i)
            vzo(i)=vzz(i)
            
          enddo
          
        endif
        
c     iteration required if ntcons > 0 and isw = 1

        mxiter=1
        if(isw.eq.1.and.ntcons.gt.0)mxiter=2
        do iter=1,mxiter

c     zero scaling matrix

          do i=1,9
            eta(i)=uni(i)
          enddo

          if(iter.eq.mxiter)then

c     calculate Berendsen barostat - taup is relaxation time
            
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

c     update positions
          
          do i=iatm0,iatm1
            
            xtmp=eta(1)*xxx(i)+eta(4)*yyy(i)+eta(7)*zzz(i)
            ytmp=eta(2)*xxx(i)+eta(5)*yyy(i)+eta(8)*zzz(i)
            ztmp=eta(3)*xxx(i)+eta(6)*yyy(i)+eta(9)*zzz(i)
            xxx(i)=tstep*vxx(i)+xtmp
            yyy(i)=tstep*vyy(i)+ytmp
            zzz(i)=tstep*vzz(i)+ztmp
            
          enddo

c     merge position data
          
          if(mxnode.gt.1)
     x      call merge(idnode,mxnode,natms,mxbuff,xxx,yyy,zzz,buffer)

c     apply shake corrections to bond constraints
          
          if(ntcons.gt.0)then
            
            if(mxnode.gt.1)
     x        call merge(idnode,mxnode,natms,mxbuff,vxx,vyy,vzz,buffer)

            safe=.false.
            call rdrattle_r
     x        (safe,lshmov,idnode,imcon,mxnode,natms,nscons,
     x        tolnce,tstep,vircon,dxx,dyy,dzz,dxt,dyt,dzt,
     x        txx,tyy,tzz,xxt,yyt,zzt,strcns)
            if(.not.safe)return
            
          endif
          
          if(iter.lt.mxiter)then

            do i=iatm0,iatm1

              xxx(i)=xxo(i)
              yyy(i)=yyo(i)
              zzz(i)=zzo(i)
              vxx(i)=vxo(i)
              vyy(i)=vyo(i)
              vzz(i)=vzo(i)

            enddo
            
          endif

        enddo

c     adjust long range corrections and number density
        
        elrc=elrc0*(volm0/volm)
        virlrc=virlrc0*(volm0/volm)
        
        do kk=1,ntpatm
          dens(kk)=dens0(kk)*(volm0/volm)
        enddo
        
c     second pass of velocity verlet algorithm
        
      else
        
c     update velocities
        
        do i=iatm0,iatm1
 
          vxx(i)=vxx(i)+0.5d0*tstep*rmass(i)*fxx(i)
          vyy(i)=vyy(i)+0.5d0*tstep*rmass(i)*fyy(i)
          vzz(i)=vzz(i)+0.5d0*tstep*rmass(i)*fzz(i)
          
        enddo
        
c     apply Berendsen thermostat taut is relaxation time
        
        engke=getkin(natms,idnode,mxnode)
        chit=sqrt(1.d0+tstep/taut*(sigma/engke-1.d0))
        
        do i=iatm0,iatm1
          
          vxx(i)=chit*vxx(i)
          vyy(i)=chit*vyy(i)
          vzz(i)=chit*vzz(i)
          
        enddo
        
c     merge velocity data
          
        if(mxnode.gt.1)
     x    call merge(idnode,mxnode,natms,mxbuff,vxx,vyy,vzz,buffer)
          
c     correct constraint bond velocities using rattle
        
        if(ntcons.gt.0)then

          safe=.false.
          call rdrattle_v
     x      (safe,idnode,mxnode,natms,nscons,tolnce,tstep,
     x      dxx,dyy,dzz,txx,tyy,tzz,xxt,yyt,zzt)
          if(.not.safe)return
          
        endif

c     kinetic contribution to stress tensor
        
        call kinstress(natms,idnode,mxnode,strkin)
        engke=0.5d0*(strkin(1)+strkin(5)+strkin(9))
        
c     add contributions to stress tensor
        
        do i=1,9
          stress(i)=stress(i)+strkin(i)+strcns(i)
        enddo
        
      endif

c     periodic boundary condition
      
      if(isw.eq.2)then
        
        call images(imcon,idnode,mxnode,natms,cell,xxx,yyy,zzz)
        if(mxnode.gt.1)
     x    call merge(idnode,mxnode,natms,mxbuff,xxx,yyy,zzz,buffer)
        
      endif
      
c     deallocate working arrays

      if(ntcons.gt.0)then

        deallocate(xxo,yyo,zzo,stat=fail(1))
        deallocate(txx,tyy,tzz,dxx,dyy,dzz,stat=fail(2))
        deallocate(dxt,dyt,dzt,xxt,yyt,zzt,stat=fail(3))

      endif
      
      return
      end subroutine nstvv_b1

      subroutine nstvv_h1
     x  (safe,lshmov,isw,idnode,mxnode,natms,imcon,nscons,ntcons,
     x  ntpatm,mode,ntshl,keyshl,tstep,taut,taup,sigma,temp,chit,
     x  consv,conint,engke,elrc,tolnce,vircon,virlrc,volm,press,
     x  chit_shl,sigma_shl)

c***********************************************************************
c     
c     dl_poly subroutine for integrating newtonian equations of
c     motion in molecular dynamics - velocity verlet incorporating
c     bond constraints via the shake/rattle algorithm
c     
c     anisotropic npt ensemble - Melchionna, Ciccotti and Holian
c     Molecular Physics 78 (1993) 533
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 
c     author    - w. smith november 2002
c     amended   - w. smith january 2005 : f90 conversion
c     adapted   - d. quigley - metadynamics
c     
c***********************************************************************

      implicit none

      integer, parameter :: nnn=7
      integer, parameter :: ncyc=5
      
      logical safe,lshmov,newjob
      integer isw,idnode,mxnode,natms,imcon,nscons,ntcons,ntpatm
      integer i,j,k,iatm0,iatm1,mxiter,iter,kk,icyc,mode
      real(8) tstep,taup,taut,sigma,temp,chit,consv,conint,chit0
      real(8) engke,elrc,tolnce,vircon,virlrc,volm,press,volm0
      real(8) elrc0,virlrc0,hstep,qstep,totmas,qmass,pmass
      real(8) cons0,cxx,cyy,czz,chip2,fstep

      real(8), allocatable :: dens0(:)
      real(8), allocatable :: xxo(:),yyo(:),zzo(:)
      real(8), allocatable :: xxt(:),yyt(:),zzt(:)
      real(8), allocatable :: txx(:),tyy(:),tzz(:)
      real(8), allocatable :: dxx(:),dyy(:),dzz(:)
      real(8), allocatable :: dxt(:),dyt(:),dzt(:)
      real(8), allocatable :: vxo(:),vyo(:),vzo(:)

      integer fail(nnn)
      real(8) com(3),vom(3),czero(9),strkin(9),eta0(9),celp(10)
      
c     metadynamics shell thermostat variables

      integer ntshl,keyshl
      real(8) sigma_shl

      logical,save :: lfirst=.true.
      real(8)      :: chit_shl  
      real(8),save :: qmass_shl
      real(8)      :: shlke

c     end metdynamics shell thermostat variables
      
      data newjob/.true./

      save newjob,totmas,volm0,elrc0,virlrc0,dens0
      save iatm0,iatm1,hstep,qstep,pmass,qmass,fstep
      
      safe=.true.

      if(newjob)then

c     block indices
        
        iatm0=(idnode*natms)/mxnode+1
        iatm1=((idnode+1)*natms)/mxnode

c     allocate density storage array

        fail(1)=0
        allocate(dens0(mxatyp),stat=fail(1))
        if(fail(1).gt.0)call error(idnode,2080)

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

c     system total mass

        totmas=getmass(natms,idnode,mxnode)

c     inertia parameter for Nose-Hoover thermostat and barostat
        
        qmass=2.d0*sigma*taut**2
        pmass=2.d0*sigma*taup**2

        newjob=.false.
        
      endif

c     allocate working arrays

      if(ntcons.gt.0)then

        do i=1,nnn
          fail(i)=0
        enddo
        
        allocate(xxo(mxatms),yyo(mxatms),zzo(mxatms),stat=fail(2))
        allocate(vxo(mxatms),vyo(mxatms),vzo(mxatms),stat=fail(3))
        allocate(xxt(mxatms),yyt(mxatms),zzt(mxatms),stat=fail(4))
        allocate(txx(mxatms),tyy(mxatms),tzz(mxatms),stat=fail(5))
        allocate(dxx(mxcons),dyy(mxcons),dzz(mxcons),stat=fail(6))
        allocate(dxt(mxcons),dyt(mxcons),dzt(mxcons),stat=fail(7))

        do i=1,nnn
          if(fail(i).gt.0)call error(idnode,2090)
        enddo

      endif

      if(lmetadyn.and.lfirst.and.(ntshl>0))then
        if(idnode.eq.0)then
          write(*,*)"Warning - Metdynamics Modification"
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

c     first stage of velocity verlet algorithm

      if(isw.eq.1)then

c     store current integration variables
        
        if(ntcons.gt.0)then
          
          chit0=chit
          cons0=conint
          do i=1,9

            czero(i)=cell(i)
            eta0(i)=eta(i)

          enddo
          do i=iatm0,iatm1
            
            xxo(i)=xxx(i)
            yyo(i)=yyy(i)
            zzo(i)=zzz(i)
            vxo(i)=vxx(i)
            vyo(i)=vyy(i)
            vzo(i)=vzz(i)
            
          enddo
          
        endif

c     subtract kinetic terms from stress tensor
        
        call kinstress(natms,idnode,mxnode,strkin)
        do i=1,9
          stress(i)=stress(i)-strkin(i)
        enddo

c     iteration necessary if ntcons > 0 and isw = 1

        mxiter=1
        if(isw.eq.1.and.ntcons.gt.0)mxiter=2
        do iter=1,mxiter

c     calculate current volume

          call dcell(cell,celp)
          volm=celp(10)

          do icyc=1,ncyc
            
c     integrate and apply nst thermostat
            
            call nstscale_t
     x        (idnode,mxnode,natms,mode,engke,temp,sigma,qstep,
     x        pmass,qmass,taut,chit,conint)
            
c     metdynamics shell thermostat
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
            
            call nstscale_p
     x        (idnode,mxnode,natms,mode,fstep,pmass,chit,press,volm)
            
c     integrate and apply nst thermostat
            
            call nstscale_t
     x        (idnode,mxnode,natms,mode,engke,temp,sigma,qstep,
     x        pmass,qmass,taut,chit,conint)

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

c     update velocities
          
          do i=iatm0,iatm1
   
            vxx(i)=vxx(i)+hstep*rmass(i)*fxx(i)
            vyy(i)=vyy(i)+hstep*rmass(i)*fyy(i)
            vzz(i)=vzz(i)+hstep*rmass(i)*fzz(i)
            
          enddo

c     calculate system centre of mass
          
          call getcom(natms,idnode,mxnode,totmas,com)

c     update positions
          
          do i=iatm0,iatm1
            
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
          
          if(mxnode.gt.1)
     x      call merge(idnode,mxnode,natms,mxbuff,xxx,yyy,zzz,buffer)

c     apply shake corrections to bond constraints
          
          if(ntcons.gt.0)then
            
            if(mxnode.gt.1)
     x        call merge(idnode,mxnode,natms,mxbuff,vxx,vyy,vzz,buffer)
            
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

            chit=chit0
            conint=cons0
            do i=1,9

              cell(i)=czero(i)
              eta(i)=eta0(i)

            enddo
            do i=iatm0,iatm1
              
              xxx(i)=xxo(i)
              yyy(i)=yyo(i)
              zzz(i)=zzo(i)
              vxx(i)=vxo(i)
              vyy(i)=vyo(i)
              vzz(i)=vzo(i)
              
            enddo

          endif

        enddo

c     second stage of velocity verlet algorithm
        
      else

c     update velocities
        
        do i=iatm0,iatm1
 
          vxx(i)=vxx(i)+hstep*rmass(i)*fxx(i)
          vyy(i)=vyy(i)+hstep*rmass(i)*fyy(i)
          vzz(i)=vzz(i)+hstep*rmass(i)*fzz(i)
          
        enddo
        
        if(ntcons.gt.0)then

c     merge velocity data
          
          if(mxnode.gt.1)
     x      call merge(idnode,mxnode,natms,mxbuff,vxx,vyy,vzz,buffer)

c     add constraint terms to stress tensor

            do i=1,9
              stress(i)=stress(i)+strcns(i)
            enddo

c     correct constraint bond velocities using rattle
          
          safe=.false.
          call rdrattle_v
     x      (safe,idnode,mxnode,natms,nscons,tolnce,tstep,
     x      dxx,dyy,dzz,txx,tyy,tzz,xxt,yyt,zzt)
          if(.not.safe)return
          
        endif

        do icyc=1,ncyc
          
c     integrate and apply nst thermostat
          
          call nstscale_t
     x      (idnode,mxnode,natms,mode,engke,temp,sigma,qstep,
     x      pmass,qmass,taut,chit,conint)
          
c     metadynamics shell thermostat
c     ====================================================
c     Must first merge update velocities as the core-shell
c     velocites are not distributed according to the same
c     rules.
c     ====================================================
          if(lmetadyn.and.keyshl.eq.1)then
            if(mxnode.gt.1)call merge
     x        (idnode,mxnode,natms,mxbuff,vxx,vyy,vzz,buffer)
            call nvtscale_shl
     x        (idnode,mxnode,ntshl,shlke,sigma_shl,hstep,qmass_shl,
     x        taut,chit_shl,conint)      
          endif  
          
c     integrate and apply nst barostat
          
          call nstscale_p
     x      (idnode,mxnode,natms,mode,fstep,pmass,chit,press,volm)
          
c     integrate and apply nst thermostat
          
          call nstscale_t
     x      (idnode,mxnode,natms,mode,engke,temp,sigma,qstep,
     x      pmass,qmass,taut,chit,conint)
          
c     metadynamics shell thermostat
c     ====================================================
c     Must first merge update velocities as the core-shell
c     velocites are not distributed according to the same
c     rules.
c     ====================================================
          if(lmetadyn.and.keyshl.eq.1)then
            if(mxnode.gt.1)call merge
     x        (idnode,mxnode,natms,mxbuff,vxx,vyy,vzz,buffer)
            call nvtscale_shl
     x        (idnode,mxnode,ntshl,shlke,sigma_shl,hstep,qmass_shl,
     x        taut,chit_shl,conint)      
          endif  
          
        enddo
        
c     remove system centre of mass velocity
        
        call getvom(natms,idnode,mxnode,totmas,vom)
        
        do i=iatm0,iatm1
          
          vxx(i)=vxx(i)-vom(1)
          vyy(i)=vyy(i)-vom(2)
          vzz(i)=vzz(i)-vom(3)
          
        enddo

c     merge velocity data
        
        if(mxnode.gt.1)
     x    call merge(idnode,mxnode,natms,mxbuff,vxx,vyy,vzz,buffer)

c     conserved quantity less kinetic and potential energy terms
        
        chip2=sdot0(9,eta,eta)
        if(mode.eq.2)chip2=chip2-eta(1)**2
        consv=conint+0.5d0*qmass*chit**2+0.5d0*pmass*chip2+press*volm

c     metadynamics shell thermostat

        if(lmetadyn.and.keyshl.eq.1)then
          consv=consv+0.5d0*qmass_shl*chit_shl**2
        endif
        
c     kinetic contribution to stress tensor
        
        call kinstress(natms,idnode,mxnode,strkin)

c     add contributions to stress tensor
        
        do i=1,9
          stress(i)=stress(i)+strkin(i)
        enddo

      endif

c     periodic boundary condition
      
      if(isw.eq.2)then
        
        call images(imcon,idnode,mxnode,natms,cell,xxx,yyy,zzz)
        if(mxnode.gt.1)
     x    call merge(idnode,mxnode,natms,mxbuff,xxx,yyy,zzz,buffer)
        
      endif
      
c     adjust long range corrections and number density
        
        elrc=elrc0*(volm0/volm)
        virlrc=virlrc0*(volm0/volm)
        do kk=1,ntpatm
          dens(kk)=dens0(kk)*(volm0/volm)
        enddo

c     deallocate working arrays

      if(ntcons.gt.0)then

        deallocate(xxo,yyo,zzo,vxo,vyo,vzo,stat=fail(1))
        deallocate(txx,tyy,tzz,dxx,dyy,dzz,stat=fail(2))
        deallocate(dxt,dyt,dzt,xxt,yyt,zzt,stat=fail(3))

      endif
      
      return
      end subroutine nstvv_h1
      
      end module vv_motion_module
