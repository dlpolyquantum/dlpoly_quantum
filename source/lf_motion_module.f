      module lf_motion_module

c***********************************************************************
c     
c     dl_poly module for verlet leap frog integration schemes
c     copyright - daresbury laboratory
c     author    - w. smith    aug 2006
c     
c***********************************************************************
      
      use config_module
      use ensemble_tools_module
      use error_module
      use property_module
      use setup_module
      use shake_module
      use site_module
      use utility_module
      
      contains

      subroutine rdshake_1
     x  (safe,lshmov,idnode,imcon,mxnode,natms,nscons,
     x  tolnce,tstep,vircon,dxx,dyy,dzz,dxt,dyt,dzt,
     x  txx,tyy,tzz,xxt,yyt,zzt,stresh)

c***********************************************************************
c     
c     dl_poly subroutine for applying bond constraint corrections after
c     atomic integration.
c     Must be used in conjunction with integration algorithms
c     
c     assume bond vectors dxx,dyy,dzz are input
c     dxx =xxx(i) - xxx(j) etc
c
c     copyright - daresbury laboratory 1994
c     author    - w. smith august 1992.
c     amended   - t. forester march 1994.
c     
c***********************************************************************
      
      implicit none
      
      logical safe,lshmov
      integer idnode,imcon,mxnode,natms,nscons,icyc,i,j,k
      real(8) tolnce,tstep,vircon,stresh,dxx,dyy,dzz,strs1,strs2
      real(8) strs3,strs5,strs6,strs9,tstep2,esig,esig1
      real(8) dis,amti,amtj,omega2,gamma,gammi,gammj,dli,dlj
      real(8) dxt,dyt,dzt,txx,tyy,tzz,xxt,yyt,zzt

      dimension stresh(9)
      dimension dxx(mxcons),dyy(mxcons),dzz(mxcons)
      dimension dxt(mxcons),dyt(mxcons),dzt(mxcons)
      dimension xxt(mxatms),yyt(mxatms),zzt(mxatms)
      dimension txx(mxatms),tyy(mxatms),tzz(mxatms)
      
c     test size of work arrays

      safe=.true.
      if(mxxdf.lt.nscons)safe=.false.

      if(mxnode.gt.1) call gstate(safe)
      if(.not.safe) call error(idnode,412)

c     timestep squared

      tstep2=tstep*tstep

c     accumulators for stress tensor

      vircon=0.d0
      strs1=0.d0
      strs2=0.d0
      strs3=0.d0
      strs5=0.d0
      strs6=0.d0
      strs9=0.d0

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
          
          dis=prmcon(listcon(k,1))
          esig1=abs(dxt(k)**2+dyt(k)**2+dzt(k)**2-dis**2)/dis
          esig=max(esig,esig1)
          
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
            omega2=dis*dis
            amti= tstep2/weight(i)
            amtj=-tstep2/weight(j)
            
            if(lstfrz(i).ne.0) amti=0.d0
            if(lstfrz(j).ne.0) amtj=0.d0
            
c     constraint force parameter
            
            gamma=(omega2-(dxt(k)**2+dyt(k)**2+dzt(k)**2))/
     x        (-2.d0*(amti-amtj)*
     x        (dxx(k)*dxt(k)+dyy(k)*dyt(k)+dzz(k)*dzt(k)))
            
c     accumulate bond virial
            
            vircon=vircon+gamma*(dxx(k)**2+dyy(k)**2+dzz(k)**2)
            
            strs1=strs1-gamma*dxx(k)*dxx(k)
            strs2=strs2-gamma*dxx(k)*dyy(k)
            strs3=strs3-gamma*dxx(k)*dzz(k)
            strs5=strs5-gamma*dyy(k)*dyy(k)
            strs6=strs6-gamma*dyy(k)*dzz(k)
            strs9=strs9-gamma*dzz(k)*dzz(k)
            
c     improve approximate atomic positions
            
            gammi=-gamma*amti
            xxt(i)=xxt(i)+dxx(k)*gammi
            yyt(i)=yyt(i)+dyy(k)*gammi
            zzt(i)=zzt(i)+dzz(k)*gammi
            
            gammj=-gamma*amtj
            xxt(j)=xxt(j)+dxx(k)*gammj
            yyt(j)=yyt(j)+dyy(k)*gammj
            zzt(j)=zzt(j)+dzz(k)*gammj
            
          enddo
          
c     transport temporary positions to other nodes
          
          if(mxnode.gt.1)then
            
            if(lshmov) call shmove
     x        (idnode,mxnode,natms,lashap,lishap,xxt,yyt,zzt,
     x        txx,tyy,tzz,buffer)
            
          endif
          
          do k=1,nscons
            
            i=listcon(k,2)
            j=listcon(k,3)
            
            dli=1.d0/dble(listme(i))
            dlj=1.d0/dble(listme(j))
            
            xxx(i)=xxx(i)+xxt(i)*dli
            yyy(i)=yyy(i)+yyt(i)*dli
            zzz(i)=zzz(i)+zzt(i)*dli
            xxx(j)=xxx(j)+xxt(j)*dlj
            yyy(j)=yyy(j)+yyt(j)*dlj
            zzz(j)=zzz(j)+zzt(j)*dlj
            
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

      endif

      return
      end subroutine rdshake_1

      subroutine nve_1
     x  (safe,lshmov,idnode,imcon,mxnode,natms,nscons,ntcons,
     x  engke,tolnce,tstep,vircon)

c***********************************************************************
c     
c     dl_poly subroutine for integrating newtonian equations of
c     motion in molecular dynamics. Verlet leapfrog With RD-SHAKE
c     
c     parallel replicated data version : block data
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith august 1992.
c     amended   - t.forester sept 1994
c     amended   - t.forester dec  1994 : block data
c     amended   - w.smith    oct  2005
c     
c***********************************************************************
      
      implicit none

      logical safe,lshmov
      integer idnode,imcon,mxnode,natms,nscons,ntcons,fail,iatm0
      integer iatm1,i,j,k
      real(8) engke,tolnce,tstep,vircon,strkin,rstep

      real(8), allocatable :: xxt(:),yyt(:),zzt(:)
      real(8), allocatable :: txx(:),tyy(:),tzz(:)
      real(8), allocatable :: uxx(:),uyy(:),uzz(:)
      real(8), allocatable :: dxx(:),dyy(:),dzz(:)
      real(8), allocatable :: dxt(:),dyt(:),dzt(:)
      real(8), allocatable :: xxo(:),yyo(:),zzo(:)
      real(8), allocatable :: vxo(:),vyo(:),vzo(:)
      
      dimension strkin(9),fail(7)
      
c     allocate working arrays

      do i=1,7
         fail(i)=0
      enddo
      allocate (xxt(mxatms),yyt(mxatms),zzt(mxatms),stat=fail(1))
      allocate (txx(mxatms),tyy(mxatms),tzz(mxatms),stat=fail(2))
      allocate (uxx(mxatms),uyy(mxatms),uzz(mxatms),stat=fail(3))
      allocate (dxx(mxcons),dyy(mxcons),dzz(mxcons),stat=fail(4))
      allocate (dxt(mxcons),dyt(mxcons),dzt(mxcons),stat=fail(5))
      allocate (xxo(msatms),yyo(msatms),zzo(msatms),stat=fail(6))
      allocate (vxo(msatms),vyo(msatms),vzo(msatms),stat=fail(7))
      do i=1,7
         if(fail(i).ne.0)call error(idnode,1380)
      enddo

      safe=.false.
      
c     block indices

      iatm0=(idnode*natms)/mxnode+1
      iatm1=((idnode+1)*natms)/mxnode

c     store initial values of position and velocity
      
      j=0
      do i=iatm0,iatm1
        
        j=j+1
        xxo(j)=xxx(i)
        yyo(j)=yyy(i)
        zzo(j)=zzz(i)
        vxo(j)=vxx(i)
        vyo(j)=vyy(i)
        vzo(j)=vzz(i)
        
      enddo

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
      
c     move atoms by leapfrog algorithm
      
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
      
c     start of bond constraint procedures

      if(ntcons.eq.0)safe=.true.
      if(ntcons.gt.0)then

c     global exchange of configuration data
        
        if(mxnode.gt.1)call merge
     x    (idnode,mxnode,natms,mxbuff,xxx,yyy,zzz,buffer)

c     apply constraint correction
        
        call rdshake_1
     x    (safe,lshmov,idnode,imcon,mxnode,natms,nscons,
     x    tolnce,tstep,vircon,dxx,dyy,dzz,dxt,dyt,dzt,
     x    txx,tyy,tzz,xxt,yyt,zzt,strcns)

c     calculate velocity correction
        
        j=0
        rstep=1.d0/tstep
        do i=iatm0,iatm1
          
          j=j+1

c     calculate corrected velocity
          
          uxx(i)=(xxx(i)-xxo(j))*rstep
          uyy(i)=(yyy(i)-yyo(j))*rstep
          uzz(i)=(zzz(i)-zzo(j))*rstep
          
c     calculate the corrected forces
          
          fxx(i)=(uxx(i)-vxo(j))*weight(i)*rstep
          fyy(i)=(uyy(i)-vyo(j))*weight(i)*rstep
          fzz(i)=(uzz(i)-vzo(j))*weight(i)*rstep
          
        enddo
        
      endif

c     calculate full timestep velocity

      do i=iatm0,iatm1

        vxx(i)=0.5d0*(vxx(i)+uxx(i))
        vyy(i)=0.5d0*(vyy(i)+uyy(i))
        vzz(i)=0.5d0*(vzz(i)+uzz(i))

      enddo
      
c     calculate kinetic energy
      
      engke=getkin(natms,idnode,mxnode)
      
c     kinetic contribution to stress tensor
      
      call kinstress(natms,idnode,mxnode,strkin)

c     add contributions to stress tensor
      
      do i=1,9
        stress(i)=stress(i)+strcns(i)+strkin(i)
      enddo

c     periodic boundary condition
      
      call images(imcon,idnode,mxnode,natms,cell,xxx,yyy,zzz)
      
c     updated velocity
      
      do i=iatm0,iatm1

        vxx(i)=uxx(i)
        vyy(i)=uyy(i)
        vzz(i)=uzz(i)

      enddo

c     global exchange of configuration data
      
      if(mxnode.gt.1)then

        call merge(idnode,mxnode,natms,mxbuff,xxx,yyy,zzz,buffer)
        call merge(idnode,mxnode,natms,mxbuff,vxx,vyy,vzz,buffer)
        if(ntcons.gt.0)call merge
     x    (idnode,mxnode,natms,mxbuff,fxx,fyy,fzz,buffer)
        
      endif

c     deallocate work arrays

      deallocate (xxt,yyt,zzt,txx,tyy,tzz,stat=fail(1))
      deallocate (uxx,uyy,uzz,dxx,dyy,dzz,stat=fail(2))
      deallocate (dxt,dyt,dzt,xxo,yyo,zzo,stat=fail(3))
      deallocate (vxo,vyo,vzo,stat=fail(4))
      
      return
      end subroutine nve_1

      subroutine nvt_e1
     x  (safe,lshmov,idnode,imcon,mxnode,natms,nscons,ntcons,
     x  engke,tolnce,tstep,vircon)
      
c***********************************************************************
c     
c     dl_poly subroutine for integrating newtonian equations of
c     motion in molecular dynamics - verlet leapfrog with Evans
c     thermostat.
c     Comp. Phys. reports 1, 299, (1984)
c     
c     parallel replicated data version : block data
c     
c     for systems using bond CONSTRAINTS.
c     
c     copyright - daresbury laboratory
c     author    - t forester july 1993
c     amended   - w.smith october 2005
c     
c***********************************************************************

      implicit none

      logical safe,lshmov
      integer idnode,imcon,mxnode,natms,nscons,ntcons,fail,iatm0
      integer iatm1,i,j,k,iter,mxiter
      real(8) engke,tolnce,tstep,vircon,strkin
      real(8) rstep,chit,viracc,strcon,vdotf

      real(8), allocatable :: xxt(:),yyt(:),zzt(:)
      real(8), allocatable :: txx(:),tyy(:),tzz(:)
      real(8), allocatable :: uxx(:),uyy(:),uzz(:)
      real(8), allocatable :: dxx(:),dyy(:),dzz(:)
      real(8), allocatable :: dxt(:),dyt(:),dzt(:)
      real(8), allocatable :: xxo(:),yyo(:),zzo(:)
      real(8), allocatable :: vxo(:),vyo(:),vzo(:)

      dimension strkin(9),strcon(9),fail(7)
      
c     allocate working arrays

      do i=1,7
         fail(i)=0
      enddo
      allocate (xxt(mxatms),yyt(mxatms),zzt(mxatms),stat=fail(1))
      allocate (txx(mxatms),tyy(mxatms),tzz(mxatms),stat=fail(2))
      allocate (uxx(mxatms),uyy(mxatms),uzz(mxatms),stat=fail(3))
      allocate (dxx(mxcons),dyy(mxcons),dzz(mxcons),stat=fail(4))
      allocate (dxt(mxcons),dyt(mxcons),dzt(mxcons),stat=fail(5))
      allocate (xxo(msatms),yyo(msatms),zzo(msatms),stat=fail(6))
      allocate (vxo(msatms),vyo(msatms),vzo(msatms),stat=fail(7))
      do i=1,7
         if(fail(i).ne.0)call error(idnode,1390)
      enddo

      safe=.false.

c     block indices

      iatm0=(idnode*natms)/mxnode+1
      iatm1=((idnode+1)*natms)/mxnode
      
c     initialise constraint virial accumulators

      vircon=0.d0
      do i=1,9
        strcns(i)=0.d0
      enddo

c     store initial positions and velocities
      
      j=0
      do i=iatm0,iatm1
        
        j=j+1
        xxo(j)=xxx(i)
        yyo(j)=yyy(i)
        zzo(j)=zzz(i)
        vxo(j)=vxx(i)
        vyo(j)=vyy(i)
        vzo(j)=vzz(i)
        
      enddo

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
        
c     begin temperature control iteration
      
      mxiter=3
      if(ntcons.eq.0)mxiter=2
      chit=0.d0
      
      do iter=1,mxiter
        
c     move atoms by leapfrog algorithm
        
        j=0
        
        do i=iatm0,iatm1

          j=j+1
          
c     update velocities
          
          uxx(i)=vxo(j)+tstep*(rmass(i)*fxx(i)-chit*vxx(i))
          uyy(i)=vyo(j)+tstep*(rmass(i)*fyy(i)-chit*vyy(i))
          uzz(i)=vzo(j)+tstep*(rmass(i)*fzz(i)-chit*vzz(i))
          
c     update positions
          
          xxx(i)=xxo(j)+tstep*uxx(i)
          yyy(i)=yyo(j)+tstep*uyy(i)
          zzz(i)=zzo(j)+tstep*uzz(i)
          
        enddo
        
c     start of bond constraint procedures
        
        if(ntcons.eq.0)safe=.true.
        if(ntcons.gt.0)then
          
c     merge configuration data
          
          if(mxnode.gt.1)call merge
     x      (idnode,mxnode,natms,mxbuff,xxx,yyy,zzz,buffer)
          
c     apply constraint correction
          
          call rdshake_1
     x      (safe,lshmov,idnode,imcon,mxnode,natms,nscons,
     x      tolnce,tstep,viracc,dxx,dyy,dzz,dxt,dyt,dzt,
     x      txx,tyy,tzz,xxt,yyt,zzt,strcon)
          
c     accumulate constraint virial terms

          vircon=vircon+viracc
          do i=1,9
            strcns(i)=strcns(i)+strcon(i)
          enddo
          
c     calculate other constraint corrections
          
          j=0
          rstep=1.d0/tstep
          do i=iatm0,iatm1
            
            j=j+1
            
c     calculate corrected velocity
            
            uxx(i)=(xxx(i)-xxo(j))*rstep
            uyy(i)=(yyy(i)-yyo(j))*rstep
            uzz(i)=(zzz(i)-zzo(j))*rstep
            
c     calculate the corrected forces
            
            fxx(i)=(uxx(i)-vxo(j))*weight(i)*rstep
            fyy(i)=(uyy(i)-vyo(j))*weight(i)*rstep
            fzz(i)=(uzz(i)-vzo(j))*weight(i)*rstep
            
          enddo
          
c     end of shake corrections

        endif

c     estimate velocity at the full step

        j=0
        do i=iatm0,iatm1
          
          j=j+1

          vxx(i)=0.5d0*(uxx(i)+vxo(j))
          vyy(i)=0.5d0*(uyy(i)+vyo(j))
          vzz(i)=0.5d0*(uzz(i)+vzo(j))

        enddo

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
        chit=vdotf/engke
        engke=0.5d0*engke
        
c     end of thermal constraint iteration
        
      enddo
      
c     kinetic contribution to stress tensor

      call kinstress(natms,idnode,mxnode,strkin)

c     add contributions to stress tensor
      
      do i=1,9
        stress(i)=stress(i)+strcns(i)+strkin(i)
      enddo

c     periodic boundary condition
      
      call images(imcon,idnode,mxnode,natms,cell,xxx,yyy,zzz)

c     updated velocity

      do i=iatm0,iatm1

        vxx(i)=uxx(i)
        vyy(i)=uyy(i)
        vzz(i)=uzz(i)

      enddo

c     global exchange of configuration data
      
      if(mxnode.gt.1)then
        
        call merge(idnode,mxnode,natms,mxbuff,xxx,yyy,zzz,buffer)
        call merge(idnode,mxnode,natms,mxbuff,vxx,vyy,vzz,buffer)
        if(ntcons.gt.0)call merge
     x    (idnode,mxnode,natms,mxbuff,fxx,fyy,fzz,buffer)
        
      endif

c     deallocate work arrays

      deallocate (xxt,yyt,zzt,txx,tyy,tzz,stat=fail(1))
      deallocate (uxx,uyy,uzz,dxx,dyy,dzz,stat=fail(2))
      deallocate (dxt,dyt,dzt,xxo,yyo,zzo,stat=fail(3))
      deallocate (vxo,vyo,vzo,stat=fail(4))
      
      return
      end subroutine nvt_e1

      subroutine nvt_b1
     x  (safe,lshmov,idnode,imcon,mxnode,natms,nscons,ntcons,
     x  engke,taut,sigma,tolnce,tstep,vircon)
      
c***********************************************************************
c     
c     dl_poly subroutine for integrating newtonian equations of
c     motion in molecular dynamics - verlet leapfrog with Berendsen
c     thermostat.
c     
c     parallel replicated data version : block data
c     
c     for systems using bond constraints
c     
c     copyright - daresbury laboratory 1993
c     author    -    t. forester   may 1993
c     amended : t.forester sept 1994
c     amended : t.forester  dec 1994 : block data
c     amended   - w.smith   oct 2005
c     
c***********************************************************************

      implicit none

      logical safe,lshmov
      integer idnode,imcon,mxnode,natms,nscons,ntcons,fail
      integer iatm0,iatm1,i,j,k,maxit,iter
      real(8) engke,taut,sigma,tolnce,tstep,vircon,strkin,viracc
      real(8) rstep,rtsq,chit0,strcon

      real(8), allocatable :: xxt(:),yyt(:),zzt(:)
      real(8), allocatable :: txx(:),tyy(:),tzz(:)
      real(8), allocatable :: uxx(:),uyy(:),uzz(:)
      real(8), allocatable :: dxx(:),dyy(:),dzz(:)
      real(8), allocatable :: dxt(:),dyt(:),dzt(:)
      real(8), allocatable :: xxo(:),yyo(:),zzo(:)
      real(8), allocatable :: xx1(:),yy1(:),zz1(:)
      real(8), allocatable :: vxo(:),vyo(:),vzo(:)

      dimension strkin(9),fail(8),strcon(9)
      
c     allocate working arrays

      do i=1,8
         fail(i)=0
      enddo
      allocate (xxt(mxatms),yyt(mxatms),zzt(mxatms),stat=fail(1))
      allocate (txx(mxatms),tyy(mxatms),tzz(mxatms),stat=fail(2))
      allocate (uxx(mxatms),uyy(mxatms),uzz(mxatms),stat=fail(3))
      allocate (dxx(mxcons),dyy(mxcons),dzz(mxcons),stat=fail(4))
      allocate (dxt(mxcons),dyt(mxcons),dzt(mxcons),stat=fail(5))
      allocate (xxo(msatms),yyo(msatms),zzo(msatms),stat=fail(6))
      allocate (xx1(msatms),yy1(msatms),zz1(msatms),stat=fail(7))
      allocate (vxo(msatms),vyo(msatms),vzo(msatms),stat=fail(8))
      do i=1,8
         if(fail(i).ne.0)call error(idnode,1400)
      enddo

      safe=.false.

c     set up block indices

      iatm0=(idnode*natms)/mxnode+1
      iatm1=((idnode+1)*natms)/mxnode
      
c     initialise constraint virial accumulators

      vircon=0.d0
      do i=1,9
        strcns(i)=0.d0
      enddo

c     store initial values of position and velocity

      j=0
      do i=iatm0,iatm1
        
        j=j+1
        xxo(j)=xxx(i)
        yyo(j)=yyy(i)
        zzo(j)=zzz(i)
        vxo(j)=vxx(i)
        vyo(j)=vyy(i)
        vzo(j)=vzz(i)

      enddo
      
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
      
c     estimate kinetic energy at full timestep

      j=0
      do i=iatm0,iatm1

        j=j+1

        vxx(i)=vxo(j)+0.5d0*tstep*rmass(i)*fxx(i)
        vyy(i)=vyo(j)+0.5d0*tstep*rmass(i)*fyy(i)
        vzz(i)=vzo(j)+0.5d0*tstep*rmass(i)*fzz(i)

      enddo

      engke=getkin(natms,idnode,mxnode)
      
c     begin iterations !!-----------------------------------------------

      maxit=3
      if(ntcons.eq.0) maxit=maxit-1
      do iter=1,maxit

c     temperature scaling  coefficient - taut is the decay constant
        
        chit0=sqrt(1.d0+tstep/taut*(sigma/engke-1.d0))

c     unconstrained new positions with thermostat
        
        j=0
        do i=iatm0,iatm1

          j=j+1

c     advance velocity using leapfrog

          uxx(i)=(vxo(j)+tstep*rmass(i)*fxx(i))*chit0
          uyy(i)=(vyo(j)+tstep*rmass(i)*fyy(i))*chit0
          uzz(i)=(vzo(j)+tstep*rmass(i)*fzz(i))*chit0

c     advance positions using leapfrog

          xxx(i)=xxo(j)+tstep*uxx(i)
          yyy(i)=yyo(j)+tstep*uyy(i)
          zzz(i)=zzo(j)+tstep*uzz(i)

c     store uncorrected positions

          xx1(j)=xxx(i)
          yy1(j)=yyy(i)
          zz1(j)=zzz(i)

        enddo

c     start of bond constraint procedures

        if(ntcons.eq.0)safe=.true.
        if(ntcons.gt.0.and.iter.eq.1)then

c     global exchange of configuration data
          
          if(mxnode.gt.1)call merge
     x      (idnode,mxnode,natms,mxbuff,xxx,yyy,zzz,buffer)

c     apply constraint correction
          
          call rdshake_1
     x      (safe,lshmov,idnode,imcon,mxnode,natms,nscons,
     x      tolnce,tstep,viracc,dxx,dyy,dzz,dxt,dyt,dzt,
     x      txx,tyy,tzz,xxt,yyt,zzt,strcon)

c     accumulate constraint virial terms

          vircon=vircon+viracc
          do i=1,9
            strcns(i)=strcns(i)+strcon(i)
          enddo

c     calculate other constraint corrections
          
          j=0
          rstep=1.d0/tstep
          rtsq=1.d0/tstep**2 
          do i=iatm0,iatm1
            
            j=j+1

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

c     calculate kinetic energy
        
        j=0
        do i=iatm0,iatm1
          
          j=j+1

c     estimate velocity at the full step

          vxx(i)=0.5d0*(uxx(i)+vxo(j))
          vyy(i)=0.5d0*(uyy(i)+vyo(j))
          vzz(i)=0.5d0*(uzz(i)+vzo(j))

        enddo

        engke=getkin(natms,idnode,mxnode)
        
c     end of thermostat iterations

      enddo

c     kinetic contribution to stress tensor

      call kinstress(natms,idnode,mxnode,strkin)

c     add contributions to stress tensor
      
      do i=1,9
        stress(i)=stress(i)+strcns(i)+strkin(i)
      enddo

c     periodic boundary condition
      
      call images(imcon,idnode,mxnode,natms,cell,xxx,yyy,zzz)

c     updated velocity

      do i=iatm0,iatm1

        vxx(i)=uxx(i)
        vyy(i)=uyy(i)
        vzz(i)=uzz(i)

      enddo

c     global exchange of configuration data
      
      if(mxnode.gt.1)then
        
        call merge(idnode,mxnode,natms,mxbuff,xxx,yyy,zzz,buffer)
        call merge(idnode,mxnode,natms,mxbuff,vxx,vyy,vzz,buffer)
        if(ntcons.gt.0)call merge
     x    (idnode,mxnode,natms,mxbuff,fxx,fyy,fzz,buffer)
        
      endif

c     deallocate work arrays

      deallocate (xxt,yyt,zzt,txx,tyy,tzz,stat=fail(1))
      deallocate (uxx,uyy,uzz,dxx,dyy,dzz,stat=fail(2))
      deallocate (dxt,dyt,dzt,xxo,yyo,zzo,stat=fail(3))
      deallocate (xx1,yy1,zz1,vxo,vyo,vzo,stat=fail(4))
      
      return
      end subroutine nvt_b1

      subroutine nvt_h1
     x  (safe,lshmov,idnode,imcon,mxnode,natms,nscons,ntcons,
     x  chit,consv,conint,engke,taut,sigma,tolnce,tstep,vircon)
      
c***********************************************************************
c     
c     dl_poly subroutine for integrating newtonian equations of
c     motion in molecular dynamics - verlet leapfrog with Hoover 
c     thermostat.
c     
c     parallel replicated data version : block data
c     
c     for systems using bond constraints
c     
c     copyright - daresbury laboratory
c     author    - t. forester may 1993
c     amended   - w.smith october 2005
c     
c***********************************************************************

      implicit none

      logical safe,lshmov
      integer idnode,imcon,mxnode,natms,nscons,ntcons,fail,i,j,k
      integer iatm0,iatm1,maxit,iter
      real(8) chit,consv,conint,engke,taut,sigma,tolnce,tstep,vircon
      real(8) strkin,rstep,rtsq,qmass,chitp,chit0,viracc
      real(8) chitnew,strcon

      real(8), allocatable :: xxt(:),yyt(:),zzt(:)
      real(8), allocatable :: txx(:),tyy(:),tzz(:)
      real(8), allocatable :: uxx(:),uyy(:),uzz(:)
      real(8), allocatable :: dxx(:),dyy(:),dzz(:)
      real(8), allocatable :: dxt(:),dyt(:),dzt(:)
      real(8), allocatable :: xxo(:),yyo(:),zzo(:)
      real(8), allocatable :: xx1(:),yy1(:),zz1(:)
      real(8), allocatable :: vxo(:),vyo(:),vzo(:)

      dimension strkin(9),fail(8),strcon(9)
      
c     allocate working arrays

      do i=1,8
         fail(i)=0
      enddo
      allocate (xxt(mxatms),yyt(mxatms),zzt(mxatms),stat=fail(1))
      allocate (txx(mxatms),tyy(mxatms),tzz(mxatms),stat=fail(2))
      allocate (uxx(mxatms),uyy(mxatms),uzz(mxatms),stat=fail(3))
      allocate (dxx(mxcons),dyy(mxcons),dzz(mxcons),stat=fail(4))
      allocate (dxt(mxcons),dyt(mxcons),dzt(mxcons),stat=fail(5))
      allocate (xxo(msatms),yyo(msatms),zzo(msatms),stat=fail(6))
      allocate (xx1(msatms),yy1(msatms),zz1(msatms),stat=fail(7))
      allocate (vxo(msatms),vyo(msatms),vzo(msatms),stat=fail(8))
      do i=1,8
         if(fail(i).ne.0)call error(idnode,1410)
      enddo

      safe=.false.

c     set up block indices

      iatm0=(idnode*natms)/mxnode+1
      iatm1=((idnode+1)*natms)/mxnode

c     inertia parameter for Nose-Hoover thermostat
      
      qmass=2.0d0*sigma*taut**2

c     initialise constraint virial accumulators

      vircon=0.d0
      do i=1,9
        strcns(i)=0.d0
      enddo

c     store initial values of position and velocity

      j=0
      do i=iatm0,iatm1
        
        j=j+1
        xxo(j)=xxx(i)
        yyo(j)=yyy(i)
        zzo(j)=zzz(i)
        vxo(j)=vxx(i)
        vyo(j)=vyy(i)
        vzo(j)=vzz(i)

      enddo
      
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
      
c     estimate velocities at full time step
      
      j=0
      do i=iatm0,iatm1
        
        j=j+1
        vxx(i)=vxo(j)+0.5d0*tstep*rmass(i)*fxx(i)
        vyy(i)=vyo(j)+0.5d0*tstep*rmass(i)*fyy(i)
        vzz(i)=vzo(j)+0.5d0*tstep*rmass(i)*fzz(i)
        
      enddo
      
c     kinetic energy at full time step
      
      engke=getkin(natms,idnode,mxnode)

c     propagate chit

      chitp=2.d0*(engke-sigma)/qmass
      chitnew=chit+tstep*chitp
      chit0=0.5d0*(chit+chitnew)

c     begin iterations !!-----------------------------------------------

      maxit=4
      if(ntcons.eq.0) maxit=maxit-1
      
      do iter=1,maxit

c     unconstrained new positions
        
        j=0
        do i=iatm0,iatm1

          j=j+1

c     advance velocity using leapfrog

          uxx(i)=vxo(j)+tstep*(fxx(i)*rmass(i)-chit0*vxx(i))
          uyy(i)=vyo(j)+tstep*(fyy(i)*rmass(i)-chit0*vyy(i))
          uzz(i)=vzo(j)+tstep*(fzz(i)*rmass(i)-chit0*vzz(i))

c     advance positions using leapfrog

          xxx(i)=xxo(j)+tstep*uxx(i)
          yyy(i)=yyo(j)+tstep*uyy(i)
          zzz(i)=zzo(j)+tstep*uzz(i)

c     store uncorrected positions

          xx1(j)=xxx(i)
          yy1(j)=yyy(i)
          zz1(j)=zzz(i)

        enddo

c     start of bond constraint procedures

        if(ntcons.eq.0)safe=.true.
        if(ntcons.gt.0.and.iter.eq.1)then

c     global exchange of configuration data
          
          if(mxnode.gt.1)call merge
     x      (idnode,mxnode,natms,mxbuff,xxx,yyy,zzz,buffer)
          
c     apply constraint correction
          
          call rdshake_1
     x      (safe,lshmov,idnode,imcon,mxnode,natms,nscons,
     x      tolnce,tstep,viracc,dxx,dyy,dzz,dxt,dyt,dzt,
     x      txx,tyy,tzz,xxt,yyt,zzt,strcon)

c     accumulate constraint virial terms

          vircon=vircon+viracc
          do i=1,9
            strcns(i)=strcns(i)+strcon(i)
          enddo
          
c     calculate other constraint corrections
          
          j=0
          rstep=1.d0/tstep
          rtsq=1.d0/tstep**2 
          do i=iatm0,iatm1
            
            j=j+1

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
        do i=iatm0,iatm1
          
          j=j+1
          vxx(i)=0.5d0*(uxx(i)+vxo(j))
          vyy(i)=0.5d0*(uyy(i)+vyo(j))
          vzz(i)=0.5d0*(uzz(i)+vzo(j))

        enddo

c     calculate kinetic energy

        engke=getkin(natms,idnode,mxnode)
        
c     improved prediction of chit 

        chitp=2.d0*(engke-sigma)/qmass
        chitnew=chit+tstep*chitp
        chit0=0.5d0*(chit+chitnew)

c     end of thermostat iterations

      enddo

c     kinetic contribution to stress tensor

      call kinstress(natms,idnode,mxnode,strkin)

c     add contributions to stress tensor
      
      do i=1,9
        stress(i)=stress(i)+strcns(i)+strkin(i)
      enddo

c     update thermostat

      chit=chitnew

c     conserved quantity less kinetic and potential energy terms
      
      conint=conint+tstep*chit0*qmass/taut**2
      consv=conint+0.5d0*qmass*chit0**2

c     periodic boundary condition
      
      call images(imcon,idnode,mxnode,natms,cell,xxx,yyy,zzz)

c     updated velocity

      do i=iatm0,iatm1

        vxx(i)=uxx(i)
        vyy(i)=uyy(i)
        vzz(i)=uzz(i)

      enddo

c     global exchange of configuration data
      
      if(mxnode.gt.1)then
        
        call merge(idnode,mxnode,natms,mxbuff,xxx,yyy,zzz,buffer)
        call merge(idnode,mxnode,natms,mxbuff,vxx,vyy,vzz,buffer)
        if(ntcons.gt.0)call merge
     x    (idnode,mxnode,natms,mxbuff,fxx,fyy,fzz,buffer)
        
      endif

c     deallocate work arrays

      deallocate (xxt,yyt,zzt,txx,tyy,tzz,stat=fail(1))
      deallocate (uxx,uyy,uzz,dxx,dyy,dzz,stat=fail(2))
      deallocate (dxt,dyt,dzt,xxo,yyo,zzo,stat=fail(3))
      deallocate (xx1,yy1,zz1,vxo,vyo,vzo,stat=fail(4))
      
      return
      end subroutine nvt_h1

      subroutine npt_b1
     x  (safe,lshmov,idnode,imcon,mxnode,natms,ntpatm,nscons,
     x  ntcons,elrc,engke,virlrc,press,taup,taut,sigma,tolnce,
     x  tstep,virtot,vircon,volm)
      
c***********************************************************************
c     
c     dl_poly subroutine for integrating newtonian equations of
c     motion in molecular dynamics - verlet leapfrog with Berendsen
c     thermostat and isotropic pressure control
c     isothermal compressibility (beta) set to that of liquid water
c     = 0.007372 dlpoly units
c     
c     parallel replicated data version
c     
c     for systems using bond CONSTRAINTS. Frozen atoms feb 1994
c     
c     copyright - daresbury laboratory 1993
c     author    - t. forester dec 1993
c     amended   - w.smith     oct 2005
c     
c***********************************************************************

      implicit none
      
      logical safe,lshmov,newjob
      integer idnode,imcon,mxnode,natms,ntpatm,nscons,ntcons
      integer fail,i,j,k,iatm0,iatm1,maxit,iter
      real(8) elrc,engke,virlrc,press,taup,taut,sigma,tolnce,tstep
      real(8) virtot,vircon,volm,strkin,beta,volm0,cell0
      real(8) elrc0,virlrc0,rstep,rtsq,psyst,chip0,scale
      real(8) chit0,viracc,strcon
      
      real(8), allocatable :: dens0(:)
      real(8), allocatable :: xxt(:),yyt(:),zzt(:)
      real(8), allocatable :: txx(:),tyy(:),tzz(:)
      real(8), allocatable :: uxx(:),uyy(:),uzz(:)
      real(8), allocatable :: dxx(:),dyy(:),dzz(:)
      real(8), allocatable :: dxt(:),dyt(:),dzt(:)
      real(8), allocatable :: xxo(:),yyo(:),zzo(:)
      real(8), allocatable :: xx1(:),yy1(:),zz1(:)
      real(8), allocatable :: vxo(:),vyo(:),vzo(:)

      dimension strkin(9),cell0(9),fail(8),strcon(9)

      save newjob,volm0,elrc0,virlrc0,dens0

      data newjob/.true./
      data beta/7.3728d-3/
      
c     allocate working arrays

      do i=1,8
         fail(i)=0
      enddo
      allocate (xxt(mxatms),yyt(mxatms),zzt(mxatms),stat=fail(1))
      allocate (txx(mxatms),tyy(mxatms),tzz(mxatms),stat=fail(2))
      allocate (uxx(mxatms),uyy(mxatms),uzz(mxatms),stat=fail(3))
      allocate (dxx(mxcons),dyy(mxcons),dzz(mxcons),stat=fail(4))
      allocate (dxt(mxcons),dyt(mxcons),dzt(mxcons),stat=fail(5))
      allocate (xxo(msatms),yyo(msatms),zzo(msatms),stat=fail(6))
      allocate (xx1(msatms),yy1(msatms),zz1(msatms),stat=fail(7))
      allocate (vxo(msatms),vyo(msatms),vzo(msatms),stat=fail(8))
      do i=1,8
         if(fail(i).ne.0)call error(idnode,1420)
      enddo

      safe=.false.

c     store initial values of volume and long range corrections

      if(newjob)then

        volm0=volm
        elrc0=elrc
        virlrc0=virlrc
        allocate (dens0(mxatyp),stat=fail(1))
        if(fail(1).ne.0)call error(idnode,1430)
        do i=1,ntpatm
          dens0(i)=dens(i)
        enddo
        newjob=.false.

      endif

c     set up block indices

      iatm0=(idnode*natms)/mxnode+1
      iatm1=((idnode+1)*natms)/mxnode
      
c     initialise constraint virial terms

      vircon=0.d0
      do i=1,9
        strcns(i)=0.d0
      enddo

c     store initial cell vectors

      do i=1,9
        cell0(i)=cell(i)
      enddo

c     store initial values of position and velocity
      
      j=0
      do i=iatm0,iatm1
        
        j=j+1
        xxo(j)=xxx(i)
        yyo(j)=yyy(i)
        zzo(j)=zzz(i)
        vxo(j)=vxx(i)
        vyo(j)=vyy(i)
        vzo(j)=vzz(i)

      enddo
      
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
      
c     estimate velocity at full timestep

      j=0
      do i=iatm0,iatm1

        j=j+1
        vxx(i)=vxo(j)+0.5d0*tstep*rmass(i)*fxx(i)
        vyy(i)=vyo(j)+0.5d0*tstep*rmass(i)*fyy(i)
        vzz(i)=vzo(j)+0.5d0*tstep*rmass(i)*fzz(i)

      enddo

c     kinetic energy at current timestep

      engke=getkin(natms,idnode,mxnode)

c     pressure control variable - taup is pressure relaxation time

      psyst=(2.d0*engke-virtot-vircon)/(3.d0*volm)
      chip0=1.d0+beta*tstep*(psyst-press)/taup
      scale=chip0**(1.d0/3.d0)
      
c     temperature scaling  coefficient - taut is temperature relaxation time

      chit0=sqrt(1.d0+tstep/taut*(sigma/engke-1.d0))

c     begin iterations !!-----------------------------------------------

      maxit=5
      if(ntcons.eq.0)maxit=maxit-1
      
      do iter=1,maxit

c     unconstrained new positions
        
        j=0
        do i=iatm0,iatm1

          j=j+1

c     advance velocity using leapfrog

          uxx(i)=(vxo(j)+tstep*rmass(i)*fxx(i))*chit0
          uyy(i)=(vyo(j)+tstep*rmass(i)*fyy(i))*chit0
          uzz(i)=(vzo(j)+tstep*rmass(i)*fzz(i))*chit0

c     update positions
          
          xxx(i)=tstep*uxx(i)+scale*xxo(j)
          yyy(i)=tstep*uyy(i)+scale*yyo(j)
          zzz(i)=tstep*uzz(i)+scale*zzo(j)

        enddo
        
c     start of bond constraint procedures

        if(ntcons.eq.0)safe=.true.
        if(ntcons.gt.0)then

c     store integrated positions

          j=0
          do i=iatm0,iatm1

            j=j+1
            xx1(j)=xxx(i)
            yy1(j)=yyy(i)
            zz1(j)=zzz(i)

          enddo

c     global exchange of configuration data
          
          if(mxnode.gt.1)call merge
     x      (idnode,mxnode,natms,mxbuff,xxx,yyy,zzz,buffer)
          
c     estimate new cell tensor
          
          do i=1,9
            cell(i)=cell0(i)*scale
          enddo
          
c     apply constraint correction
          
          call rdshake_1
     x      (safe,lshmov,idnode,imcon,mxnode,natms,nscons,
     x      tolnce,tstep,viracc,dxx,dyy,dzz,dxt,dyt,dzt,
     x      txx,tyy,tzz,xxt,yyt,zzt,strcon)

c     accumulate constraint virial terms

          vircon=vircon+viracc
          do i=1,9
            strcns(i)=strcns(i)+strcon(i)
          enddo
          
c     calculate other constraint corrections
          
          j=0
          rstep=1.d0/tstep
          rtsq=1.d0/tstep**2 
          do i=iatm0,iatm1
            
            j=j+1

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
        do i=iatm0,iatm1
          
          j=j+1
          vxx(i)=0.5d0*(uxx(i)+vxo(j))
          vyy(i)=0.5d0*(uyy(i)+vyo(j))
          vzz(i)=0.5d0*(uzz(i)+vzo(j))

        enddo

c     calculate kinetic energy

        engke=getkin(natms,idnode,mxnode)
        
c     improved prediction of chip and chit

        psyst=(2.d0*engke-virtot-vircon)/(3.d0*volm)
        chip0=1.d0+beta*tstep*(psyst-press)/taup
        scale=chip0**(1.d0/3.d0)
        chit0=sqrt(1.d0+tstep/taut*(sigma/engke-1.d0))
        
c     end of thermostat iterations

      enddo

c     kinetic contribution to stress tensor

      call kinstress(natms,idnode,mxnode,strkin)

c     add contributions to stress tensor
      
      do i=1,9
        stress(i)=stress(i)+strcns(i)+strkin(i)
      enddo

c     update volume

      volm=volm*chip0

c     scale cell vectors - isotropic

      do i=1,9
        cell(i)=cell0(i)*scale
      enddo

c     construct scaling tensor (for use with tethers)

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

c     periodic boundary condition
      
      call images(imcon,idnode,mxnode,natms,cell,xxx,yyy,zzz)

c     updated velocity

      do i=iatm0,iatm1

        vxx(i)=uxx(i)
        vyy(i)=uyy(i)
        vzz(i)=uzz(i)

      enddo

c     global exchange of configuration data
      
      if(mxnode.gt.1)then
        
        call merge(idnode,mxnode,natms,mxbuff,xxx,yyy,zzz,buffer)
        call merge(idnode,mxnode,natms,mxbuff,vxx,vyy,vzz,buffer)
        if(ntcons.gt.0)call merge
     x    (idnode,mxnode,natms,mxbuff,fxx,fyy,fzz,buffer)
        
      endif

c     deallocate work arrays

      deallocate (xxt,yyt,zzt,txx,tyy,tzz,stat=fail(1))
      deallocate (uxx,uyy,uzz,dxx,dyy,dzz,stat=fail(2))
      deallocate (dxt,dyt,dzt,xxo,yyo,zzo,stat=fail(3))
      deallocate (xx1,yy1,zz1,vxo,vyo,vzo,stat=fail(4))
      
      return
      end subroutine npt_b1

      subroutine npt_h1
     x  (safe,lshmov,idnode,imcon,mxnode,natms,ntpatm,nscons,
     x  ntcons,chip,chit,conint,consv,elrc,engke,virlrc,press,
     x  taup,taut,sigma,temp,tolnce,tstep,virtot,vircon,volm)
      
c***********************************************************************
c     
c     dl_poly subroutine for integrating newtonian equations of
c     motion in molecular dynamics - verlet leapfrog with Hoover 
c     thermostat+piston.
c     
c     reference: Melchionna, Ciccotti and Holian,
c     Mol Phys 1993, 78, p533
c     
c     parallel replicated data version
c     
c     for systems using bond constraints (using atomic pressure)
c     
c     copyright daresbury laboratory 1995
c     author    -    s. melchionna   april 1995
c     and       -    t. forester     april 1995
c     amended   -    w. smith     october  2005
c     
c***********************************************************************

      implicit none
      
      logical safe,lshmov,newjob
      integer idnode,imcon,mxnode,natms,ntpatm,nscons,ntcons
      integer i,j,k,iatm0,iatm1,fail,maxit,iter
      real(8) chip,chit,conint,consv,elrc,engke,virlrc,press
      real(8) taup,taut,sigma,temp,tolnce,tstep,virtot,vircon,volm
      real(8) strcon,volm0,elrc0,virlrc0,rstep,rtsq,qmass
      real(8) chipnew,chitp,chitnew,chit0,volnew,scale,viracc,vold
      real(8) cons1,cons2,cons3,strkin,cell0
      real(8) pmass,totmas,chipp,chip0,com,vom

      real(8), allocatable :: dens0(:)
      real(8), allocatable :: xxt(:),yyt(:),zzt(:)
      real(8), allocatable :: txx(:),tyy(:),tzz(:)
      real(8), allocatable :: uxx(:),uyy(:),uzz(:)
      real(8), allocatable :: dxx(:),dyy(:),dzz(:)
      real(8), allocatable :: dxt(:),dyt(:),dzt(:)
      real(8), allocatable :: xxo(:),yyo(:),zzo(:)
      real(8), allocatable :: xx1(:),yy1(:),zz1(:)
      real(8), allocatable :: vxo(:),vyo(:),vzo(:)

      dimension strcon(9),fail(8),strkin(9),com(3),vom(3),cell0(9)
      
      save newjob,volm0,elrc0,virlrc0,cell0,dens0

      data newjob/.true./

c     allocate working arrays

      do i=1,8
         fail(i)=0
      enddo
      allocate (xxt(mxatms),yyt(mxatms),zzt(mxatms),stat=fail(1))
      allocate (txx(mxatms),tyy(mxatms),tzz(mxatms),stat=fail(2))
      allocate (uxx(mxatms),uyy(mxatms),uzz(mxatms),stat=fail(3))
      allocate (dxx(mxcons),dyy(mxcons),dzz(mxcons),stat=fail(4))
      allocate (dxt(mxcons),dyt(mxcons),dzt(mxcons),stat=fail(5))
      allocate (xxo(msatms),yyo(msatms),zzo(msatms),stat=fail(6))
      allocate (xx1(msatms),yy1(msatms),zz1(msatms),stat=fail(7))
      allocate (vxo(msatms),vyo(msatms),vzo(msatms),stat=fail(8))
      do i=1,8
         if(fail(i).ne.0)call error(idnode,1440)
      enddo

      safe=.false.

c     store initial values of volume, long range corrections etc

      if(newjob)then

        volm0=volm
        elrc0=elrc
        virlrc0=virlrc
        allocate (dens0(mxatyp),stat=fail(1))
        if(fail(1).ne.0)call error(idnode,1450)
        do i=1,ntpatm
          dens0(i)=dens(i)
        enddo
        do i=1,9
          cell0(i)=cell(i)
        enddo
        newjob=.false.

      endif

c     set up block indices

      iatm0=(idnode*natms)/mxnode+1
      iatm1=((idnode+1)*natms)/mxnode
      
c     inertia parameter for Nose-Hoover thermostat
      
      qmass=2.0d0*sigma*taut**2
      pmass=2.0d0*sigma*taup**2

c     initialise constraint virial accumulators

      vircon=0.d0
      do i=1,9
        strcns(i)=0.d0
      enddo

c     store initial values of position and velocity
      
      j=0
      do i=iatm0,iatm1
        
        j=j+1
        xxo(j)=xxx(i)
        yyo(j)=yyy(i)
        zzo(j)=zzz(i)
        vxo(j)=vxx(i)
        vyo(j)=vyy(i)
        vzo(j)=vzz(i)

      enddo

c     total system mass

      totmas=getmass(natms,idnode,mxnode)

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
      
c     calculate centre of mass

      call getcom(natms,idnode,mxnode,totmas,com)

c     estimate kinetic energy at current timestep

      j=0
      do i=iatm0,iatm1

        j=j+1

c     estimate velocity at the full step

        vxx(i)=vxo(j)+0.5d0*tstep*rmass(i)*fxx(i)
        vyy(i)=vyo(j)+0.5d0*tstep*rmass(i)*fyy(i)
        vzz(i)=vzo(j)+0.5d0*tstep*rmass(i)*fzz(i)

      enddo

c     calculate kinetic energy

      engke=getkin(natms,idnode,mxnode)
        
c     propagate chip

      chipp=(2.d0*engke-virtot-vircon-3.d0*press*volm)/pmass-
     x  chit*chip
      chipnew=chip+tstep*chipp
      chip0=0.5d0*(chip+chipnew)

c     propagate chit

      chitp=(2.d0*(engke-sigma)+pmass*chip**2-boltz*temp)/qmass
      chitnew=chit+tstep*chitp
      chit0=0.5d0*(chit+chitnew)

c     begin iterations !!-----------------------------------------------

      maxit=5
      if(ntcons.eq.0) maxit=maxit-1
      
      do iter=1,maxit

c     unconstrained new positions
        
        j=0
        do i=iatm0,iatm1

          j=j+1

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
        
c     start of bond constraint procedures

        if(ntcons.eq.0)safe=.true.
        if(ntcons.gt.0)then

c     store integrated positions

          j=0
          do i=iatm0,iatm1

            j=j+1
            xx1(j)=xxx(i)
            yy1(j)=yyy(i)
            zz1(j)=zzz(i)

          enddo

c     global exchange of configuration data
          
          if(mxnode.gt.1)call merge
     x      (idnode,mxnode,natms,mxbuff,xxx,yyy,zzz,buffer)

c     estimate new cell tensor

          volnew=volm*exp(3.d0*tstep*chipnew)
          scale=(volnew/volm0)**(1.d0/3.d0)
          do i=1,9
            cell(i)=cell0(i)*scale
          enddo

c     apply constraint correction
          
          call rdshake_1
     x      (safe,lshmov,idnode,imcon,mxnode,natms,nscons,
     x      tolnce,tstep,viracc,dxx,dyy,dzz,dxt,dyt,dzt,
     x      txx,tyy,tzz,xxt,yyt,zzt,strcon)

c     accumulate constraint virial terms

          vircon=vircon+viracc
          do i=1,9
            strcns(i)=strcns(i)+strcon(i)
          enddo

c     calculate other constraint corrections
          
          j=0
          rstep=1.d0/tstep
          rtsq=1.d0/tstep**2
          do i=iatm0,iatm1
            
            j=j+1

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
        do i=iatm0,iatm1
          
          j=j+1

          vxx(i)=0.5d0*(uxx(i)+vxo(j))
          vyy(i)=0.5d0*(uyy(i)+vyo(j))
          vzz(i)=0.5d0*(uzz(i)+vzo(j))

        enddo

c     calculate kinetic energy

        engke=getkin(natms,idnode,mxnode)
        
c     improved prediction of chip and chit 

        chipp=(2.d0*engke-virtot-vircon-3.d0*press*volm)/pmass-
     x    chit0*chip0
        chipnew=chip+tstep*chipp
        chip0=0.5d0*(chip+chipnew)

        chitp=(2.d0*(engke-sigma)+pmass*chip0**2-boltz*temp)/qmass
        chitnew=chit+tstep*chitp
        chit0=0.5d0*(chit+chitnew)

c     end of thermostat iterations

      enddo

c     kinetic contribution to stress tensor

      call kinstress(natms,idnode,mxnode,strkin)

c     add contributions to stress tensor
      
      do i=1,9
        stress(i)=stress(i)+strcns(i)+strkin(i)
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
      cons1=0.5d0*qmass*chit0**2
      cons2=press*vold
      cons3=0.5d0*pmass*chip0**2
      consv=conint+cons1+cons2+cons3

c     periodic boundary condition
      
      call images(imcon,idnode,mxnode,natms,cell,xxx,yyy,zzz)
      
c     updated velocity
      
      do i=iatm0,iatm1

        vxx(i)=uxx(i)
        vyy(i)=uyy(i)
        vzz(i)=uzz(i)

      enddo

c     ensure total momentum is zero

      call getvom(natms,idnode,mxnode,totmas,vom)

      do i=iatm0,iatm1
        
        vxx(i)=vxx(i)-vom(1)
        vyy(i)=vyy(i)-vom(2)
        vzz(i)=vzz(i)-vom(3)
        
      enddo

c     global exchange of configuration data
      
      if(mxnode.gt.1)then
        
        call merge(idnode,mxnode,natms,mxbuff,xxx,yyy,zzz,buffer)
        call merge(idnode,mxnode,natms,mxbuff,vxx,vyy,vzz,buffer)
        if(ntcons.gt.0)call merge
     x    (idnode,mxnode,natms,mxbuff,fxx,fyy,fzz,buffer)
        
      endif

c     deallocate work arrays

      deallocate (xxt,yyt,zzt,txx,tyy,tzz,stat=fail(1))
      deallocate (uxx,uyy,uzz,dxx,dyy,dzz,stat=fail(2))
      deallocate (dxt,dyt,dzt,xxo,yyo,zzo,stat=fail(3))
      deallocate (xx1,yy1,zz1,vxo,vyo,vzo,stat=fail(4))
      
      return
      end subroutine npt_h1

      subroutine nst_b1
     x  (safe,lshmov,idnode,imcon,mxnode,natms,ntpatm,nscons,
     x  ntcons,mode,elrc,engke,virlrc,press,taup,taut,sigma,tolnce,
     x  tstep,vircon,volm)
      
c***********************************************************************
c     
c     dl_poly subroutine for integrating newtonian equations of
c     motion in molecular dynamics - verlet leapfrog with Berendsen
c     thermostat and anisotropic pressure control
c     isothermal compressibility (beta) set to that of liquid water
c     = 0.007372 dlpoly units
c     
c     parallel replicated data version
c     
c     for systems using bond CONSTRAINTS. Frozen atoms feb 1994
c     
c     copyright - daresbury laboratory 1993
c     author    - t. forester december 1993
c     amended   - w. smith  october    2005
c     
c***********************************************************************

      implicit none

      logical safe,lshmov,newjob
      integer idnode,imcon,mxnode,natms,ntpatm,nscons,ntcons,mode
      integer fail,i,j,k,iatm0,iatm1,maxit,iter
      real(8) elrc,engke,virlrc,press,taup,taut,sigma,tolnce,tstep
      real(8) vircon,volm,beta,volm0,elrc0,virlrc0,rstep,rtsq,chit0
      real(8) viracc,strkin,strcon,cell0,stres0,uni
      
      real(8), allocatable :: dens0(:)
      real(8), allocatable :: xxt(:),yyt(:),zzt(:)
      real(8), allocatable :: txx(:),tyy(:),tzz(:)
      real(8), allocatable :: uxx(:),uyy(:),uzz(:)
      real(8), allocatable :: dxx(:),dyy(:),dzz(:)
      real(8), allocatable :: dxt(:),dyt(:),dzt(:)
      real(8), allocatable :: xxo(:),yyo(:),zzo(:)
      real(8), allocatable :: xx1(:),yy1(:),zz1(:)
      real(8), allocatable :: vxo(:),vyo(:),vzo(:)

      dimension strkin(9),strcon(9),cell0(9),fail(8),stres0(9),uni(9)

      save newjob,volm0,elrc0,virlrc0,dens0

      data uni/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/
      data beta/7.3728d-3/
      data newjob/.true./
      
c     allocate working arrays

      do i=1,8
         fail(i)=0
      enddo
      allocate (xxt(mxatms),yyt(mxatms),zzt(mxatms),stat=fail(1))
      allocate (txx(mxatms),tyy(mxatms),tzz(mxatms),stat=fail(2))
      allocate (uxx(mxatms),uyy(mxatms),uzz(mxatms),stat=fail(3))
      allocate (dxx(mxcons),dyy(mxcons),dzz(mxcons),stat=fail(4))
      allocate (dxt(mxcons),dyt(mxcons),dzt(mxcons),stat=fail(5))
      allocate (xxo(msatms),yyo(msatms),zzo(msatms),stat=fail(6))
      allocate (xx1(msatms),yy1(msatms),zz1(msatms),stat=fail(7))
      allocate (vxo(msatms),vyo(msatms),vzo(msatms),stat=fail(8))
      do i=1,8
         if(fail(i).ne.0)call error(idnode,1460)
      enddo

      safe=.false.

c     store initial values of volume, long range corrections etc

      if(newjob)then

        volm0=volm
        elrc0=elrc
        virlrc0=virlrc
        allocate (dens0(mxatyp),stat=fail(1))
        if(fail(1).ne.0)call error(idnode,1470)
        do i=1,ntpatm
          dens0(i)=dens(i)
        enddo
        newjob=.false.

      endif

c     set up block indices

      iatm0=(idnode*natms)/mxnode+1
      iatm1=((idnode+1)*natms)/mxnode
      
c     initialise constraint virial accumulators

      vircon=0.d0
      do i=1,9
        strcns(i)=0.d0
      enddo

c     store original cell vectors and stress tensor

      do i=1,9
        
        cell0(i)=cell(i)
        stres0(i)=stress(i)
        
      enddo

c     store initial values of position and velocity

      j=0
      do i=iatm0,iatm1
        
        j=j+1
        xxo(j)=xxx(i)
        yyo(j)=yyy(i)
        zzo(j)=zzz(i)
        vxo(j)=vxx(i)
        vyo(j)=vyy(i)
        vzo(j)=vzz(i)

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
      
c     estimate kinetic energy at current timestep

      j=0
      do i=iatm0,iatm1

        j=j+1
        vxx(i)=vxo(j)+0.5d0*tstep*rmass(i)*fxx(i)
        vyy(i)=vyo(j)+0.5d0*tstep*rmass(i)*fyy(i)
        vzz(i)=vzo(j)+0.5d0*tstep*rmass(i)*fzz(i)

      enddo

c     calculate kinetic energy

      engke=getkin(natms,idnode,mxnode)
        
c     kinetic contribution to stress tensor

      call kinstress(natms,idnode,mxnode,strkin)

c     current estimate of stres tensor

      do i=1,9
        stress(i)=stres0(i)+strkin(i)
      enddo
      
c     initial estimate of eta matrix and chit 

      chit0=sqrt(1.d0+tstep/taut*(sigma/engke-1.d0))
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
      
c     begin iterations !!-----------------------------------------------

      maxit=5
      if(ntcons.eq.0)maxit=maxit-1
      
      do iter=1,maxit

c     unconstrained new positions
        
        j=0
        do i=iatm0,iatm1

          j=j+1

c     advance velocity using leapfrog

          uxx(i)=(vxo(j)+tstep*rmass(i)*fxx(i))*chit0
          uyy(i)=(vyo(j)+tstep*rmass(i)*fyy(i))*chit0
          uzz(i)=(vzo(j)+tstep*rmass(i)*fzz(i))*chit0

c     update positions
          
          xxx(i)=tstep*uxx(i)+
     x      eta(1)*xxo(j)+eta(4)*yyo(j)+eta(7)*zzo(j)
          yyy(i)=tstep*uyy(i)+
     x      eta(2)*xxo(j)+eta(5)*yyo(j)+eta(8)*zzo(j)
          zzz(i)=tstep*uzz(i)+
     x      eta(3)*xxo(j)+eta(6)*yyo(j)+eta(9)*zzo(j)

        enddo
        
c     start of bond constraint procedures
        
        if(ntcons.eq.0)safe=.true.
        if(ntcons.gt.0)then

c     store integrated positions

          j=0
          do i=iatm0,iatm1

            j=j+1
            xx1(j)=xxx(i)
            yy1(j)=yyy(i)
            zz1(j)=zzz(i)

          enddo

c     global exchange of configuration data
          
          if(mxnode.gt.1)call merge
     x      (idnode,mxnode,natms,mxbuff,xxx,yyy,zzz,buffer)
          
c     estimate new cell parameters
      
          call mat_mul(eta,cell0,cell)

c     apply constraint correction
          
          call rdshake_1
     x      (safe,lshmov,idnode,imcon,mxnode,natms,nscons,
     x      tolnce,tstep,viracc,dxx,dyy,dzz,dxt,dyt,dzt,
     x      txx,tyy,tzz,xxt,yyt,zzt,strcon)

c     accumulate constraint virial terms

          vircon=vircon+viracc
          do i=1,9
            strcns(i)=strcns(i)+strcon(i)
          enddo

c     calculate other constraint corrections
          
          j=0
          rstep=1.d0/tstep
          rtsq=1.d0/tstep**2 
          do i=iatm0,iatm1
            
            j=j+1

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
        do i=iatm0,iatm1
          
          j=j+1

          vxx(i)=0.5d0*(uxx(i)+vxo(j))
          vyy(i)=0.5d0*(uyy(i)+vyo(j))
          vzz(i)=0.5d0*(uzz(i)+vzo(j))

        enddo

c     calculate kinetic energy

        engke=getkin(natms,idnode,mxnode)
        
c     kinetic contribution to stress tensor
        
        call kinstress(natms,idnode,mxnode,strkin)
        
c     add contributions to stress tensor
        
        do i=1,9
          stress(i)=stres0(i)+strcns(i)+strkin(i)
        enddo

c     improved calculation of eta matrix and chit 

        chit0=sqrt(1.d0+tstep/taut*(sigma/engke-1.d0))
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
        
c     end of thermostat and barostat iterations

      enddo

c     update volume
      
      volm=volm*eta(1)*eta(5)*eta(9)

c     adjust cell vectors - anisotropic

      call mat_mul(eta,cell0,cell)

c     adjust long range corrections and number density
      
      elrc=elrc0*(volm0/volm)
      virlrc=virlrc0*(volm0/volm)
      
      do k=1,ntpatm
        dens(k)=dens0(k)*(volm0/volm)
      enddo

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

      deallocate (xxt,yyt,zzt,txx,tyy,tzz,stat=fail(1))
      deallocate (uxx,uyy,uzz,dxx,dyy,dzz,stat=fail(2))
      deallocate (dxt,dyt,dzt,xxo,yyo,zzo,stat=fail(3))
      deallocate (xx1,yy1,zz1,vxo,vyo,vzo,stat=fail(4))
      
      return
      end subroutine nst_b1

      subroutine nst_h1
     x  (safe,lshmov,idnode,imcon,mxnode,natms,ntpatm,nscons,
     x  ntcons,mode,chit,conint,consv,elrc,engke,virlrc,press,
     x  taup,taut,sigma,temp,tolnce,tstep,vircon,volm)
      
c***********************************************************************
c     
c     dl_poly subroutine for integrating newtonian equations of
c     motion in molecular dynamics - verlet leapfrog with Hoover 
c     thermostat+piston.
c     
c     Parrinello - Rahman type : changing cell shape.
c     
c     reference: Melchionna, Ciccotti and Holian,
c     Mol Phys 1993, 78, p533
c     
c     parallel replicated data version
c     
c     for systems using bond constraints (using atomic pressure)
c     
c     copyright daresbury laboratory 1995
c     author    -    t. forester     june  1995
c     
c***********************************************************************

      implicit none
      
      logical safe,lshmov,newjob
      integer idnode,imcon,mxnode,natms,ntpatm,nscons,ntcons
      integer fail,i,j,k,iatm0,iatm1,maxit,iter,mode
      real(8) chip,chit,conint,consv,elrc,engke,virlrc,press
      real(8) taup,taut,sigma,temp,tolnce,tstep,vircon,volm
      real(8) strcon,strkin,etanew,eta0,cell0,volm0,elrc0,virlrc0
      real(8) rstep,rtsq,pmass,qmass,totmas,com,vom,uni,fac
      real(8) chitp,chitnew,chit0,xxa,yya,zza,etadot
      real(8) viracc,cons1,cons2,cons3,vold,stres0

      real(8), allocatable :: dens0(:)
      real(8), allocatable :: xxt(:),yyt(:),zzt(:)
      real(8), allocatable :: txx(:),tyy(:),tzz(:)
      real(8), allocatable :: uxx(:),uyy(:),uzz(:)
      real(8), allocatable :: dxx(:),dyy(:),dzz(:)
      real(8), allocatable :: dxt(:),dyt(:),dzt(:)
      real(8), allocatable :: xxo(:),yyo(:),zzo(:)
      real(8), allocatable :: xx1(:),yy1(:),zz1(:)
      real(8), allocatable :: vxo(:),vyo(:),vzo(:)

      dimension strkin(9),strcon(9),fail(8),com(3),vom(3)
      dimension etanew(9),eta0(9),cell0(9),stres0(9),uni(9)

      save newjob,volm0,elrc0,virlrc0,dens0

      data newjob/.true./,uni/1.d0,3*0.d0,1.d0,3*0.d0,1.d0/
      
c     allocate working arrays
      
      do i=1,8
         fail(i)=0
      enddo
      allocate (xxt(mxatms),yyt(mxatms),zzt(mxatms),stat=fail(1))
      allocate (txx(mxatms),tyy(mxatms),tzz(mxatms),stat=fail(2))
      allocate (uxx(mxatms),uyy(mxatms),uzz(mxatms),stat=fail(3))
      allocate (dxx(mxcons),dyy(mxcons),dzz(mxcons),stat=fail(4))
      allocate (dxt(mxcons),dyt(mxcons),dzt(mxcons),stat=fail(5))
      allocate (xxo(msatms),yyo(msatms),zzo(msatms),stat=fail(6))
      allocate (xx1(msatms),yy1(msatms),zz1(msatms),stat=fail(7))
      allocate (vxo(msatms),vyo(msatms),vzo(msatms),stat=fail(8))
      do i=1,8
         if(fail(i).ne.0)call error(idnode,1480)
      enddo

      safe=.false.

c     store initial values of volume, long range corrections etc

      if(newjob)then

        volm0=volm
        elrc0=elrc
        virlrc0=virlrc
        allocate (dens0(mxatyp),stat=fail(1))
        if(fail(1).ne.0)call error(idnode,1490)
        do i=1,ntpatm
          dens0(i)=dens(i)
        enddo
        newjob=.false.

      endif

c     set up block indices

      iatm0=(idnode*natms)/mxnode+1
      iatm1=((idnode+1)*natms)/mxnode
      
c     inertia parameter for Nose-Hoover thermostat
      
      qmass=2.0d0*sigma*taut**2
      pmass=2.0d0*sigma*taup**2

c     initialise constraint virial accumulators

      vircon=0.d0
      do i=1,9
        strcns(i)=0.d0
      enddo

c     store original cell vectors and stress tensor

      do i=1,9
        
        cell0(i)=cell(i)
        stres0(i)=stress(i)
        
      enddo

c     store initial values of position and velocity
      
      j=0
      do i=iatm0,iatm1
        
        j=j+1
        xxo(j)=xxx(i)
        yyo(j)=yyy(i)
        zzo(j)=zzz(i)
        vxo(j)=vxx(i)
        vyo(j)=vyy(i)
        vzo(j)=vzz(i)

      enddo

c     total system mass

      totmas=getmass(natms,idnode,mxnode)

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
      
c     calculate system centre of mass

      call getcom(natms,idnode,mxnode,totmas,com)

c     estimate kinetic energy at current timestep

      j=0
      do i=iatm0,iatm1

        j=j+1

c     estimate position at current time step

        xxx(i)=xxo(j)+tstep*(vxo(j)+tstep*rmass(i)*fxx(i))
        yyy(i)=yyo(j)+tstep*(vyo(j)+tstep*rmass(i)*fyy(i))
        zzz(i)=zzo(j)+tstep*(vzo(j)+tstep*rmass(i)*fzz(i))
        
c     estimate velocity at the full step

        vxx(i)=vxo(j)+0.5d0*tstep*rmass(i)*fxx(i)
        vyy(i)=vyo(j)+0.5d0*tstep*rmass(i)*fyy(i)
        vzz(i)=vzo(j)+0.5d0*tstep*rmass(i)*fzz(i)

      enddo

c     calculate kinetic energy

      engke=getkin(natms,idnode,mxnode)
        
c     kinetic contribution to stress tensor
        
      call kinstress(natms,idnode,mxnode,strkin)

c     initial estimate of stress tensor

      do i=1,9
        stress(i)=stres0(i)+strkin(i)
      enddo

c     propagation of eta

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

      etadot=sdot0(9,eta0,eta0)
      if(mode.eq.2)etadot=etadot-eta0(1)**2
      chitp=(2.d0*(engke-sigma)+pmass*etadot-fac*boltz*temp)/qmass
      chitnew=chit+tstep*chitp
      chit0=0.5d0*(chit+chitnew)

c     begin iterations !!-----------------------------------------------

      maxit=5
      if(ntcons.eq.0) maxit=maxit-1
      do iter=1,maxit

c     unconstrained new positions
        
        j=0
        do i=iatm0,iatm1

          j=j+1

c     advance velocity using leapfrog

          uxx(i)=vxo(j)+tstep*(fxx(i)*rmass(i)-(eta0(1)+chit0)*vxx(i)-
     x      eta0(4)*vyy(i)-eta0(7)*vzz(i))
          uyy(i)=vyo(j)+tstep*(fyy(i)*rmass(i)-(eta0(5)+chit0)*vyy(i)-
     x      eta0(2)*vxx(i)-eta0(8)*vzz(i))
          uzz(i)=vzo(j)+tstep*(fzz(i)*rmass(i)-(eta0(9)+chit0)*vzz(i)-
     x      eta0(3)*vxx(i)-eta0(6)*vyy(i))

c     advance positions using leapfrog

          xxa=(xxx(i)+xxo(j))*0.5d0-com(1)
          yya=(yyy(i)+yyo(j))*0.5d0-com(2)
          zza=(zzz(i)+zzo(j))*0.5d0-com(3)

          xxx(i)=xxo(j)+tstep*(uxx(i)+
     x      etanew(1)*xxa+etanew(4)*yya+etanew(7)*zza)
          yyy(i)=yyo(j)+tstep*(uyy(i)+
     x      etanew(2)*xxa+etanew(5)*yya+etanew(8)*zza)
          zzz(i)=zzo(j)+tstep*(uzz(i)+
     x      etanew(3)*xxa+etanew(6)*yya+etanew(9)*zza)

        enddo
        
c     start of bond constraint procedures
        
        if(ntcons.eq.0)safe=.true.
        if(ntcons.gt.0)then

c     store integrated positions

          j=0
          do i=iatm0,iatm1

            j=j+1
            xx1(j)=xxx(i)
            yy1(j)=yyy(i)
            zz1(j)=zzz(i)

          enddo

c     global exchange of configuration data
          
          if(mxnode.gt.1)call merge
     x      (idnode,mxnode,natms,mxbuff,xxx,yyy,zzz,buffer)

c     estimate new cell parameters
        
        do i=1,9
          cell(i)=cell0(i)
        enddo
        call cell_propagate(tstep,cell,etanew)

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
          do i=iatm0,iatm1
            
            j=j+1

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
        do i=iatm0,iatm1
          
          j=j+1

          vxx(i)=0.5d0*(uxx(i)+vxo(j))
          vyy(i)=0.5d0*(uyy(i)+vyo(j))
          vzz(i)=0.5d0*(uzz(i)+vzo(j))

        enddo

c     calculate kinetic energy

        engke=getkin(natms,idnode,mxnode)
        
c     kinetic contribution to stress tensor
        
        call kinstress(natms,idnode,mxnode,strkin)
        
c     add contributions to stress tensor
        
        do i=1,9
          stress(i)=stres0(i)+strcns(i)+strkin(i)
        enddo

c     improved prediction of eta

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

c     improved prediction of chit

        etadot=sdot0(9,eta0,eta0)
        if(mode.eq.2)etadot=etadot-eta0(1)**2
        chitp=(2.d0*(engke-sigma)+pmass*etadot-fac*boltz*temp)/qmass
        chitnew=chit+tstep*chitp
        chit0=0.5d0*(chit+chitnew)

c     end of thermostat iterations

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

c     adjust cell vectors - anisotropic

      do i=1,9
        cell(i)=cell0(i)
      enddo
      call cell_propagate(tstep,cell,eta)
      
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

c     restore half step velocity

      do i=iatm0,iatm1
        
        vxx(i)=uxx(i)
        vyy(i)=uyy(i)
        vzz(i)=uzz(i)
        
      enddo

c     ensure total momentum is zero

      call getvom(natms,idnode,mxnode,totmas,vom)

      do i=iatm0,iatm1
        
        vxx(i)=vxx(i)-vom(1)
        vyy(i)=vyy(i)-vom(2)
        vzz(i)=vzz(i)-vom(3)
        
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

      deallocate (xxt,yyt,zzt,txx,tyy,tzz,stat=fail(1))
      deallocate (uxx,uyy,uzz,dxx,dyy,dzz,stat=fail(2))
      deallocate (dxt,dyt,dzt,xxo,yyo,zzo,stat=fail(3))
      deallocate (xx1,yy1,zz1,vxo,vyo,vzo,stat=fail(4))
      
      return
      end subroutine nst_h1
      
      end module lf_motion_module
