      module optimiser_module
      
c***********************************************************************
c     
c     dl_poly module for defining structural optimiser routines
c     copyright - daresbury laboratory
c     author    - w. smith    aug 2006
c     
c***********************************************************************
      
      use config_module
      use rigid_body_module
      use setup_module
      use shake_module
      use site_module
      use utility_module
      
      real(8), allocatable :: xxt(:),yyt(:),zzt(:)
      real(8), allocatable :: ggx(:),ggy(:),ggz(:)
      real(8), allocatable :: hhx(:),hhy(:),hhz(:)
      real(8), allocatable :: oxx(:),oyy(:),ozz(:)
      real(8), allocatable :: txx(:),tyy(:),tzz(:)
      real(8), allocatable :: uxx(:),uyy(:),uzz(:)
      real(8), allocatable :: dtx(:),dty(:),dtz(:)
      real(8), allocatable :: dxx(:),dyy(:),dzz(:)
      
      save hhx,hhy,hhz
      
      contains
      
      subroutine optimisation_selector
     x  (loptim,stropt,lzero,idnode,mxnode,natms,imcon,ntcons,
     x  nscons,ngrp,ntfree,keystr,keytol,engcfg,tstep,opttol)
      
c***********************************************************************
c     
c     dl_poly routine for selecting and running a nominated 
c     structure optimisation algorithm using energy minimisation
c     
c     copyright - daresbury laboratory
c     author    - w. smith june 2006
c     
c***********************************************************************
      
      implicit none
      
      logical loptim,stropt,lzero
      integer idnode,mxnode,natms,imcon,nscons,ngrp,ntfree,keystr
      integer keytol,ntcons
      real(8) engcfg,tstep,opttol,hnorm,grad0,grad1,ff1,sgn
      
      save grad0,grad1,ff1,sgn,hnorm

      stropt=.false.
      
      if(loptim)then
        
c     conjugate gradient structure optimisation
        
        call strucopt
     x    (stropt,keystr,keytol,idnode,mxnode,natms,ntcons,nscons,
     x    imcon,ngrp,ntfree,tstep,opttol,engcfg,hnorm,grad0,grad1,
     x    ff1,sgn)
        
      else if(lzero)then
        
c     zero kelvin structure optimisation 
        
        call zero_kelvin
     x    (stropt,idnode,mxnode,imcon,natms,ngrp,ntfree,opttol)
        
      endif
      
      return
      end subroutine optimisation_selector

      subroutine zero_kelvin
     x  (stropt,idnode,mxnode,imcon,natms,ngrp,ntfree,opttol)
      
c***********************************************************************
c     
c     dl_poly routine for zero Kelvin temperature optimization
c     if velocity.Force < 0 then velocity is set to zero in 
c     preparation for integration of equations of motion
c     
c     parallel replicated data version : block data
c     
c     copyright daresbury laboratory 1994
c     author t.forester     march 1994
c     amended t.forester    dec 1994 : block data
c     
c***********************************************************************
      
      implicit none

      logical stropt
      integer idnode,mxnode,imcon,natms,ngrp,ntfree,fail,i
      integer iatm0,iatm1,igrp1,igrp2,ifre1,ifre2,jr,ig,j,id
      real(8) dot,fsq,fcomx,fcomy,fcomz,trx,try,trz,tax,tay,taz
      real(8) rot,ggg,opttol
      
      dimension rot(9)

      data fail/0/

c     allocate work arrays

      allocate (xxt(mxatms),yyt(mxatms),zzt(mxatms),stat=fail)
      if(fail.ne.0)call error(idnode,1920)
      
c     block indices
      
      iatm0=(idnode*natms)/mxnode+1
      iatm1=((idnode+1)*natms)/mxnode
      
c     check for convergence of forces 
      
      ggg=0.d0
      do i=iatm0,iatm1
        ggg=ggg+fxx(i)**2+fyy(i)**2+fzz(i)**2
      enddo
      
      if(mxnode.gt.1)then
        buffer(1)=ggg
        call gdsum(buffer(1),1,buffer(2))
        ggg=buffer(1)
      endif
      
c     check convergence condition for forces
      
      if(opttol.ge.abs(ggg)/dble(natms))then
        
        stropt=.true.
        return
        
      endif

      if(ngrp.eq.0) then

c     take component of velocity in direction of force
        
        do i=iatm0,iatm1

          dot=vxx(i)*fxx(i)+vyy(i)*fyy(i)+vzz(i)*fzz(i)
          if(dot.lt.0.d0) then

            vxx(i)=0.d0
            vyy(i)=0.d0
            vzz(i)=0.d0

          else

c     take component of velocity in direction of force

            fsq=(fxx(i)**2+fyy(i)**2+fzz(i)**2)
            fsq=dot/max(1.d-10,fsq)
            vxx(i)=fxx(i)*fsq
            vyy(i)=fyy(i)*fsq
            vzz(i)=fzz(i)*fsq

          endif
          
        enddo
        
      else

c     block indices for groups and free atoms

        igrp1=(idnode*ngrp)/mxnode+1
        igrp2=((idnode+1)*ngrp)/mxnode

        ifre1=(idnode*ntfree)/mxnode+1
        ifre2=((idnode+1)*ntfree)/mxnode

        do j=ifre1,ifre2

c     reset atomic velocities 
          
          i=lstfre(j)
          
          dot=vxx(i)*fxx(i)+vyy(i)*fyy(i)+vzz(i)*fzz(i)
          if(dot.lt.0.d0) then

            vxx(i)=0.d0
            vyy(i)=0.d0
            vzz(i)=0.d0

          else

c     take component of velocity in direction of force

            fsq=(fxx(i)**2+fyy(i)**2+fzz(i)**2)
            fsq=dot/max(1.d-10,fsq)
            vxx(i)=fxx(i)*fsq
            vyy(i)=fyy(i)*fsq
            vzz(i)=fzz(i)*fsq
            
          endif
          
        enddo

        jr=0
        do ig=igrp1,igrp2

c     reset rigid body velocites (linear and angular)
          
          fcomx=0.d0
          fcomy=0.d0
          fcomz=0.d0

          id=lstgtp(ig)
          do j=1,numgsit(id)
            
            jr=jr+1
            i=lstrgd(jr)

c     forces on com
            
            fcomx=fcomx+fxx(i)
            fcomy=fcomy+fyy(i)
            fcomz=fcomz+fzz(i)
            
          enddo
          
          dot=gvxx(ig)*fcomx+gvyy(ig)*fcomy+gvzz(ig)*fcomz
          if(dot.lt.0.d0) then

            gvxx(ig)=0.d0
            gvyy(ig)=0.d0
            gvzz(ig)=0.d0

          else

c     take component of velocity in direction of force

            fsq=(fcomx**2+fcomy**2+fcomz**2)
            fsq=dot/max(1.d-10,fsq)
            gvxx(ig)=fcomx*fsq
            gvyy(ig)=fcomy*fsq
            gvzz(ig)=fcomz*fsq

          endif

        enddo

c     site to com distances
        
        jr=0
        do ig=igrp1,igrp2
          
          do j=1,numgsit(lstgtp(ig))
            
            jr=jr+1
            i=lstrgd(jr)
            
            xxt(jr)=xxx(i)-gcmx(ig)
            yyt(jr)=yyy(i)-gcmy(ig)
            zzt(jr)=zzz(i)-gcmz(ig)
            
          enddo
          
        enddo

c     minimum images
        
        call images(imcon,0,1,jr,cell,xxt,yyt,zzt)

c     calculate torques in lab frame
        
        jr=0
        do ig=igrp1,igrp2
          
          trx=0.d0
          try=0.d0
          trz=0.d0
          
          id=lstgtp(ig)
          do j=1,numgsit(id)
            
            jr=jr+1
            i=lstrgd(jr)
            
            trx=trx+yyt(jr)*fzz(i)-zzt(jr)*fyy(i)
            try=try+zzt(jr)*fxx(i)-xxt(jr)*fzz(i)
            trz=trz+xxt(jr)*fyy(i)-yyt(jr)*fxx(i)
            
          enddo
          
          rot(1)=q0(ig)**2+q1(ig)**2-q2(ig)**2-q3(ig)**2
          rot(2)=2.d0*(q1(ig)*q2(ig)-q0(ig)*q3(ig))
          rot(3)=2.d0*(q1(ig)*q3(ig)+q0(ig)*q2(ig))
          rot(4)=2.d0*(q1(ig)*q2(ig)+q0(ig)*q3(ig))
          rot(5)=q0(ig)**2-q1(ig)**2+q2(ig)**2-q3(ig)**2
          rot(6)=2.d0*(q2(ig)*q3(ig)-q0(ig)*q1(ig))
          rot(7)=2.d0*(q1(ig)*q3(ig)-q0(ig)*q2(ig))
          rot(8)=2.d0*(q2(ig)*q3(ig)+q0(ig)*q1(ig))
          rot(9)=q0(ig)**2-q1(ig)**2-q2(ig)**2+q3(ig)**2

c     transform to body fixed frame
          
          tax=(trx*rot(1)+try*rot(4)+trz*rot(7))*rotinx(id,2)
          tay=(trx*rot(2)+try*rot(5)+trz*rot(8))*rotiny(id,2)
          taz=(trx*rot(3)+try*rot(6)+trz*rot(9))*rotinz(id,2)
          
          dot=omx(ig)*tax+omy(ig)*tay+omz(ig)*taz
          if(dot.le.0.d0) then

            omx(ig)=0.d0
            omy(ig)=0.d0
            omz(ig)=0.d0

          else

c     take component of velocity in direction of torque
            
            fsq=(tax**2+tay**2+taz**2)
            fsq=dot/max(1.d-10,fsq)
            omx(ig)=tax*fsq
            omy(ig)=tay*fsq
            omz(ig)=taz*fsq
            
          endif
          
        enddo
        
      endif
      
c     deallocate work arrays

      deallocate (xxt,yyt,zzt,stat=fail)
      
      return
      end subroutine zero_kelvin

      subroutine strucopt
     x  (stropt,keystr,keytol,idnode,mxnode,natms,ntcons,nscons,
     x  imcon,ngrp,ntfree,tstep,opttol,fnew,hnorm,grad0,grad1,
     x  ff1,sgn)
      
c***********************************************************************
c     
c     dl_poly subroutine for optimising molecular structures
c     based on conjugate gradient method
c     
c     copyright - daresbury laboratory
c     author    - w. smith    dec 2005
c     
c     note. basis of minimisation criterion :
c           keytol=0 : absolute force
c           keytol=1 : absolute energy
c           keytol=2 : absolute displacement
c     
c***********************************************************************
      
      implicit none
      
      integer, parameter :: nnn=8
      
      logical stropt,newjob,engchk
      integer keystr,keytol,idnode,mxnode,natms,i,j,k
      integer iatm0,iatm1,fail,ngrp,ntcons,nscons,jr
      integer imcon,ig,jf,id,ntfree,igrp1,igrp2,ifre1,ifre2
      real(8) hnorm,grad0,grad1,grad2,ff1,stride,tstep,step
      real(8) ggg,fnew,fff,gam2,sgn,opttol,dischk
            
      dimension fail(nnn)
            
      save iatm0,iatm1,igrp1,igrp2,engchk,ifre1,ifre2,newjob
      
      data newjob/.true./,engchk/.false./
      
c     define initial data
      
      do i=1,nnn
        fail(i)=0
      enddo
      if(newjob)then
        allocate(hhx(mxatms),hhy(mxatms),hhz(mxatms),stat=fail(1))
      endif        
      allocate(ggx(mxatms),ggy(mxatms),ggz(mxatms),stat=fail(2))
      allocate(oxx(mxatms),oyy(mxatms),ozz(mxatms),stat=fail(3))
      allocate(dtx(mxatms),dty(mxatms),dtz(mxatms),stat=fail(4))
      if(ngrp.gt.0)then
        
        allocate(txx(mxatms),tyy(mxatms),tzz(mxatms),stat=fail(5))
        allocate(uxx(mxatms),uyy(mxatms),uzz(mxatms),stat=fail(6))
        
      endif
      if(ntcons.gt.0)then
        
        if(ngrp.eq.0)
     x    allocate(txx(mxatms),tyy(mxatms),tzz(mxatms),stat=fail(5))
        allocate(dxx(mxatms),dyy(mxatms),dzz(mxatms),stat=fail(7))
        allocate(xxt(mxatms),yyt(mxatms),zzt(mxatms),stat=fail(8))
        
      endif
      
      do i=1,nnn
        if(fail(i).ne.0)call error(idnode,1925)
      enddo
        
      if(newjob)then
        
c     define atoms for this node
        
        iatm0=(idnode*natms)/mxnode+1
        iatm1=((idnode+1)*natms)/mxnode
        
c     group block indices
        
        igrp1=(idnode*ngrp)/mxnode+1
        igrp2=((idnode+1)*ngrp)/mxnode
        
c     free atom block indices
        
        ifre1=(idnode*ntfree)/mxnode+1
        ifre2=((idnode+1)*ntfree)/mxnode
        
        newjob=.false.
        
      endif
      
c     step length for relaxation
      
      if(ntcons.gt.0)then
        step=tstep**2
      else
        step=10.d0*tstep**2
      endif
      
c     current system configuration energy
      
      fff=fnew
      
c     initialise conjugate gradient position arrays
      
      do i=1,natms
        
        oxx(i)=xxx(i)
        oyy(i)=yyy(i)
        ozz(i)=zzz(i)
        ggx(i)=0.d0
        ggy(i)=0.d0
        ggz(i)=0.d0
        
      enddo
      
c     define constraint bonds
      
      if(ntcons.gt.0)then
        
c     calculate constraint bond vector

        do k=1,nscons
          
          i=listcon(k,2)
          j=listcon(k,3)
          
          dxx(k)=xxx(i)-xxx(j)
          dyy(k)=yyy(i)-yyy(j)
          dzz(k)=zzz(i)-zzz(j)
          
        enddo
        
c     periodic boundary condition
        
        call images(imcon,0,1,nscons,cell,dxx,dyy,dzz)
        
c     calculate pseudo forces for constraint bonds
        
        call pseudo_shake(nscons,natms,mxnode,fff)
        
        do i=1,natms
          
          ggx(i)=fxx(i)+ggx(i)
          ggy(i)=fyy(i)+ggy(i)
          ggz(i)=fzz(i)+ggz(i)
          
        enddo
      
      else
        
        do i=1,natms
          
          ggx(i)=fxx(i)
          ggy(i)=fyy(i)
          ggz(i)=fzz(i)
          
        enddo
        
      endif
      
c     calculate pseudo forces for rigid bodies
      
      if(ngrp.gt.0)call torque_split
     x  (ngrp,idnode,mxnode,imcon,ggx,ggy,ggz,txx,tyy,tzz,
     x  uxx,uyy,uzz,dtx,dty,dtz) 

c     determine magnitude of 3N force vector
      
      ggg=0.d0
      
      if(ngrp.eq.0)then
        
        do i=iatm0,iatm1
          ggg=ggg+ggx(i)**2+ggy(i)**2+ggz(i)**2
        enddo
        
      else
        
        do jf=ifre1,ifre2
          
          i=lstfre(jf)
          ggg=ggg+ggx(i)**2+ggy(i)**2+ggz(i)**2
          
        enddo
        
        jr=0
        do ig=igrp1,igrp2
          
          id=lstgtp(ig)
          do j=1,numgsit(id)
            
            jr=jr+1
            i=lstrgd(jr)
            ggg=ggg+ggx(i)**2+ggy(i)**2+ggz(i)**2
            
          enddo
        
        enddo
        
      endif
      
      if(mxnode.gt.1)then
        buffer(1)=ggg
        call gdsum(buffer(1),1,buffer(2))
        ggg=buffer(1)
      endif
      ggg=sqrt(ggg)
      
c     check convergence condition for forces
      
      if(keytol.eq.0.and.opttol.ge.abs(ggg)/dble(natms))stropt=.true.
      
      if(keystr.eq.0) then
        
c     set original search direction
        
        ff1=fff
        hnorm=ggg
        grad0=ggg
        grad1=ggg
        
        if(ngrp.eq.0)then
          
          do i=iatm0,iatm1
            
            hhx(i)=ggx(i)
            hhy(i)=ggy(i)
            hhz(i)=ggz(i)
            oxx(i)=oxx(i)+step*hhx(i)
            oyy(i)=oyy(i)+step*hhy(i)
            ozz(i)=ozz(i)+step*hhz(i)
            
          enddo
          
        else
          
          do jf=ifre1,ifre2
            
            i=lstfre(jf)
            hhx(i)=ggx(i)
            hhy(i)=ggy(i)
            hhz(i)=ggz(i)
            
          enddo
          
          jr=0
          do ig=igrp1,igrp2
            
            id=lstgtp(ig)
            do j=1,numgsit(id)
              
              jr=jr+1
              i=lstrgd(jr)
              hhx(i)=ggx(i)
              hhy(i)=ggy(i)
              hhz(i)=ggz(i)
              
            enddo
            
          enddo
          
          call turn_rigid_body
     x      (igrp1,igrp2,ifre1,ifre2,step,hhx,hhy,hhz,
     x      uxx,uyy,uzz,txx,tyy,tzz,oxx,oyy,ozz)
          
        endif
        
        keystr=1
        sgn=1.d0
        
      elseif(keystr.eq.1)then
        
c     check convergence condition for energy

        if(engchk.and.keytol.eq.1.and.
     x    opttol.ge.abs(fff-ff1))stropt=.true.
        engchk=.false.
        
c     line search along chosen direction
        
        ff1=fff
        grad2=grad1
        grad1=0.d0
        do i=iatm0,iatm1
          grad1=grad1+hhx(i)*ggx(i)+hhy(i)*ggy(i)+hhz(i)*ggz(i)
        enddo
        if(mxnode.gt.1)then
          buffer(1)=grad1
          call gdsum(buffer(1),1,buffer(2))
          grad1=buffer(1)
        endif
        grad1=sgn*grad1/hnorm
        
c     linear extrapolation to minimum
        
        stride=sgn*step
        if(grad1.lt.0.d0)then
          
          keystr=2
          stride=sgn*step*grad1/(grad2-grad1)
          
        endif
        
        if(ngrp.eq.0)then
          
          do i=iatm0,iatm1
            
            oxx(i)=oxx(i)+stride*hhx(i)
            oyy(i)=oyy(i)+stride*hhy(i)
            ozz(i)=ozz(i)+stride*hhz(i)
            
          enddo
          
        else

          call turn_rigid_body
     x      (igrp1,igrp2,ifre1,ifre2,stride,hhx,hhy,hhz,
     x      uxx,uyy,uzz,txx,tyy,tzz,oxx,oyy,ozz)
          
        endif
        
      elseif(keystr.eq.2)then
        
c     construct conjugate search vector
        
        ff1=fff
        gam2=(ggg/grad0)**2
        hnorm=0.d0
        grad0=ggg
        grad1=0.d0
        do i=iatm0,iatm1
          
          hhx(i)=ggx(i)+gam2*hhx(i)
          hhy(i)=ggy(i)+gam2*hhy(i)
          hhz(i)=ggz(i)+gam2*hhz(i)
          hnorm=hnorm+hhx(i)**2+hhy(i)**2+hhz(i)**2
          grad1=grad1+hhx(i)*ggx(i)+hhy(i)*ggy(i)+hhz(i)*ggz(i)
          
        enddo
        if(mxnode.gt.1)then
          
          buffer(1)=hnorm
          buffer(2)=grad1
          call gdsum(buffer(1),2,buffer(3))
          hnorm=buffer(1)
          grad1=buffer(2)
          
        endif
        hnorm=sqrt(hnorm)
        grad1=grad1/hnorm
        sgn=sign(1.d0,grad1)
        grad1=sgn*grad1
        stride=sgn*step

        if(ngrp.eq.0)then
          
          do i=iatm0,iatm1
            
            oxx(i)=oxx(i)+stride*hhx(i)
            oyy(i)=oyy(i)+stride*hhy(i)
            ozz(i)=ozz(i)+stride*hhz(i)
            
          enddo
        
        else

          call turn_rigid_body
     x      (igrp1,igrp2,ifre1,ifre2,stride,hhx,hhy,hhz,
     x      uxx,uyy,uzz,txx,tyy,tzz,oxx,oyy,ozz)
          
        endif

        engchk=.true.
        keystr=1
        
      endif 
      
c     merge coordinate arrays
      
      if(mxnode.gt.1)then
        
        if(ngrp.eq.0)then
          
          call merge
     x      (idnode,mxnode,natms,mxbuff,oxx,oyy,ozz,buffer)
          
        else
          
          call merge1
     x      (idnode,mxnode,natms,lstme,oxx,oyy,ozz,buffer)
          
        endif
        
      endif
      
c     reassign atomic positions and calculate max displacement
      
      dischk=0.d0
      do i=1,natms
        
        dischk=max(dischk,(xxx(i)-oxx(i))**2+
     x    (yyy(i)-oyy(i))**2+(zzz(i)-ozz(i))**2)

        xxx(i)=oxx(i)
        yyy(i)=oyy(i)
        zzz(i)=ozz(i)
        
      enddo
      
c     check convergence condition for position

      if(keytol.eq.2.and.keystr.gt.0.and.
     x  opttol.ge.sqrt(dischk))stropt=.true.
      
c     deallocate working arrays
      
      deallocate(ggx,ggy,ggz,dtx,dty,dtz,oxx,oyy,ozz,stat=fail(1))
      if(ngrp.gt.0)then
        deallocate(txx,tyy,tzz,uxx,uyy,uzz,stat=fail(2))
      endif
      if(ntcons.gt.0)then
        
        deallocate(dxx,dyy,dzz,xxt,yyt,zzt,stat=fail(3))
        if(ngrp.eq.0)deallocate(txx,tyy,tzz,stat=fail(4))
        
      endif
      
      return
      end subroutine strucopt
      
      subroutine pseudo_shake(nscons,natms,mxnode,fff)

c***********************************************************************
c     
c     dl_poly subroutine treating rigid bonds as stiff harmonic bonds
c     suitable for conjugate gradient minimisation
c     
c     copyright - daresbury laboratory
c     author    - w. smith    may 2006
c     
c***********************************************************************
      
      implicit none

      real(8), parameter :: harm=1.d6

      integer i,j,k,natms,nscons,mxnode
      real(8) fff,engbnd,dis,rrr,gamma
      
c     calculate energy and force
      
      engbnd=0.d0
      
      do k=1,nscons
        
        i=listcon(k,2)
        j=listcon(k,3)
        
        dis=prmcon(listcon(k,1))
        rrr=sqrt(dxx(k)**2+dyy(k)**2+dzz(k)**2)
        engbnd=engbnd+0.5d0*harm*(rrr-dis)**2
        gamma=harm*(rrr-dis)/rrr
        ggx(i)=ggx(i)-dxx(k)*gamma
        ggy(i)=ggy(i)-dyy(k)*gamma
        ggz(i)=ggz(i)-dzz(k)*gamma
        
        ggx(j)=ggx(j)+dxx(k)*gamma
        ggy(j)=ggy(j)+dyy(k)*gamma
        ggz(j)=ggz(j)+dzz(k)*gamma
        
      enddo
      
c     global sum of pseudo forces

      call global_sum_forces(natms,mxnode,ggx,ggy,ggz)
      if(mxnode.gt.1)then
        buffer(1)=engbnd
        call gdsum(buffer(1),1,buffer(2))
        engbnd=buffer(1)
      endif
      fff=fff+engbnd
      
      return
      end subroutine pseudo_shake

      subroutine torque_split
     x  (ngrp,idnode,mxnode,imcon,ggx,ggy,ggz,txx,tyy,tzz,
     x  uxx,uyy,uzz,dtx,dty,dtz) 

c***********************************************************************
c     
c     dl_poly subroutine for resolving torques into equivalent atomic
c     forces suitable for conjugate gradient minimisation
c     
c     copyright - daresbury laboratory
c     author    - w. smith    may 2006
c     
c***********************************************************************
      
      implicit none
      
      integer i,j,ig,id,jr,jrs,ngrp,igrp1,igrp2,idnode,imcon,mxnode

      real(8) fmx,fmy,fmz,tqx,tqy,tqz,trq,txx,tyy,tzz
      real(8) ggx,ggy,ggz,tmp,taq,scale
      real(8) uxx,uyy,uzz,dtx,dty,dtz

      dimension ggx(mxatms),ggy(mxatms),ggz(mxatms)
      dimension txx(mxatms),tyy(mxatms),tzz(mxatms)
      dimension uxx(mxatms),uyy(mxatms),uzz(mxatms)
      dimension dtx(mxatms),dty(mxatms),dtz(mxatms)
      
c     group block indices
        
      igrp1=(idnode*ngrp)/mxnode+1
      igrp2=((idnode+1)*ngrp)/mxnode
      
c     calculate centres of mass of rigid bodies
      
      jr=0
      do ig=igrp1,igrp2
        
c     working com is first site in group
        
        i=lstrgd(jr+1)
        txx(ig)=xxx(i)
        tyy(ig)=yyy(i)
        tzz(ig)=zzz(i)
        
        id=lstgtp(ig)
        
        do j=1,numgsit(id)
          
          jr=jr+1
          i=lstrgd(jr)
          uxx(jr)=xxx(i)-txx(ig)
          uyy(jr)=yyy(i)-tyy(ig)
          uzz(jr)=zzz(i)-tzz(ig)
          
        enddo
        
      enddo
      
c     minimum image from working com
      
      call images(imcon,0,1,jr,cell,uxx,uyy,uzz)
      
      jr=0
      do ig=igrp1,igrp2
        
        gcmx(ig)=0.d0
        gcmy(ig)=0.d0
        gcmz(ig)=0.d0
        
        id=lstgtp(ig)
        
        do j=1,numgsit(id)
          
          jr=jr+1
          i=lstrgd(jr)
          
          gcmx(ig)=gcmx(ig)+weight(i)*uxx(jr)
          gcmy(ig)=gcmy(ig)+weight(i)*uyy(jr)
          gcmz(ig)=gcmz(ig)+weight(i)*uzz(jr)
          
        enddo
        
c     final centre of mass
        
        gcmx(ig)=gcmx(ig)/gmass(id)+txx(ig)
        gcmy(ig)=gcmy(ig)/gmass(id)+tyy(ig)
        gcmz(ig)=gcmz(ig)/gmass(id)+tzz(ig)
        
      enddo
      
c     calculate atom displacements from rigid body com

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
      
c     minimum images
      
      call images(imcon,0,1,jr,cell,dtx,dty,dtz)
        
c     resolve rigid body forces and torques to orthogonal atomic basis 
      
      jr=0
      do ig=igrp1,igrp2
        
        id=lstgtp(ig)
        
c     calculate net force on rigid body          
        
        jrs=jr
        fmx=0.d0
        fmy=0.d0
        fmz=0.d0
        do j=1,numgsit(id)
          
          jr=jr+1
          i=lstrgd(jr)
          fmx=fmx+ggx(i)
          fmy=fmy+ggy(i)
          fmz=fmz+ggz(i)
          
        enddo
        fmx=fmx/dble(numgsit(id))
        fmy=fmy/dble(numgsit(id))
        fmz=fmz/dble(numgsit(id))

c     calculate torque on rigid body
        
        jr=jrs
        tqx=0.d0
        tqy=0.d0
        tqz=0.d0
        do j=1,numgsit(id)
          
          jr=jr+1
          i=lstrgd(jr)
          
          tqx=tqx+dty(jr)*ggz(i)-dtz(jr)*ggy(i)
          tqy=tqy+dtz(jr)*ggx(i)-dtx(jr)*ggz(i)
          tqz=tqz+dtx(jr)*ggy(i)-dty(jr)*ggx(i)
          
        enddo
        
c     magnitude of torque
        
        trq=sqrt(tqx**2+tqy**2+tqz**2)
        
c     construct unit vectors for new site forces
        
        jr=jrs
        do j=1,numgsit(id)
          
          jr=jr+1
          i=lstrgd(jr)
          
          txx(i)=dty(jr)*tqz-tqy*dtz(jr)
          tyy(i)=dtz(jr)*tqx-tqz*dtx(jr)
          tzz(i)=dtx(jr)*tqy-tqx*dty(jr)
          tmp=sqrt(txx(i)**2+tyy(i)**2+tzz(i)**2)
          if(tmp.gt.1.d-10)then
            
            txx(i)=txx(i)/tmp
            tyy(i)=tyy(i)/tmp
            tzz(i)=tzz(i)/tmp
          
          else
            
            txx(i)=0.d0
            tyy(i)=0.d0
            tzz(i)=0.d0
            
          endif
          
        enddo
        
c     construct unit vectors for site location
        
        jr=jrs
        do j=1,numgsit(id)
          
          jr=jr+1
          i=lstrgd(jr)
          uxx(i)=(tyy(i)*tqz-tqy*tzz(i))/trq
          uyy(i)=(tzz(i)*tqx-tqz*txx(i))/trq
          uzz(i)=(txx(i)*tqy-tqx*tyy(i))/trq
          
        enddo
        
c     scale unit vectors to working lengths
        
        jr=jrs
        taq=0.d0
        do j=1,numgsit(id)
          
          jr=jr+1
          i=lstrgd(jr)
          tmp=dtx(jr)*uxx(i)+dty(jr)*uyy(i)+dtz(jr)*uzz(i)
          taq=taq+tmp**2
          txx(i)=tmp*txx(i)
          tyy(i)=tmp*tyy(i)
          tzz(i)=tmp*tzz(i)
          uxx(i)=tmp*uxx(i)
          uyy(i)=tmp*uyy(i)
          uzz(i)=tmp*uzz(i)
          
        enddo
        
c     calculate force scale factor
        
        scale=trq/taq
        
c     final site forces
        
        jr=jrs
        do j=1,numgsit(id)
          
          jr=jr+1
          i=lstrgd(jr)
          
          txx(i)=scale*txx(i)
          tyy(i)=scale*tyy(i)
          tzz(i)=scale*tzz(i)
          ggx(i)=fmx
          ggy(i)=fmy
          ggz(i)=fmz
          
        enddo
        
      enddo
      
      return
      end subroutine torque_split

      subroutine turn_rigid_body
     x  (igrp1,igrp2,ifre1,ifre2,step,hhx,hhy,hhz,
     x  uxx,uyy,uzz,txx,tyy,tzz,oxx,oyy,ozz)
      
c***********************************************************************
c     
c     dl_poly routine for updating positions of atoms in a rigid body
c     during a conjugate gradient minimisation
c     
c     copyright daresbury laboratory
c     author w.smith       may  2006
c     
c     note: coz=cos(theta)-1
c           zin=sin(theta)/theta
c     
c***********************************************************************
      
      implicit none
      
      integer i,j,jr,jf,ig,id,igrp1,igrp2,ifre1,ifre2
      real(8) step,hhx,hhy,hhz,uxx,uyy,uzz,txx,tyy,tzz
      real(8) oxx,oyy,ozz,uuu,ttt,the2,coz,zin

      dimension hhx(mxatms),hhy(mxatms),hhz(mxatms)
      dimension oxx(mxatms),oyy(mxatms),ozz(mxatms)
      dimension txx(mxatms),tyy(mxatms),tzz(mxatms)
      dimension uxx(mxatms),uyy(mxatms),uzz(mxatms)
      
c     update free atom positions

      do jf=ifre1,ifre2
        
        i=lstfre(jf)
        oxx(i)=oxx(i)+step*hhx(i)
        oyy(i)=oyy(i)+step*hhy(i)
        ozz(i)=ozz(i)+step*hhz(i)
        
      enddo
      
c     update rigid body atoms
      
      jr=0
      do ig=igrp1,igrp2
        
        id=lstgtp(ig)
        do j=1,numgsit(id)
          
          jr=jr+1
          i=lstrgd(jr)
          uuu=uxx(i)**2+uyy(i)**2+uzz(i)**2
          if(uuu.gt.1.d-10)then
            
            ttt=txx(i)**2+tyy(i)**2+tzz(i)**2
            the2=(ttt/uuu)*step**2
            
            coz=-the2*(1.d0-the2*(1.d0-the2*(1.d0-the2*(1.d0-the2*(1.d0-
     x        the2/132.d0)/90.d0)/56.d0)/30.d0)/12.d0)/2.d0
            zin=-the2*(1.d0-the2*(1.d0-the2*(1.d0-the2*(1.d0-the2*(1.d0-
     x        the2/156.d0)/110.d0)/72.d0)/42.d0)/20.d0)/6.d0+1.d0
            
            oxx(i)=oxx(i)+coz*uxx(i)+step*(hhx(i)+zin*txx(i))
            oyy(i)=oyy(i)+coz*uyy(i)+step*(hhy(i)+zin*tyy(i))
            ozz(i)=ozz(i)+coz*uzz(i)+step*(hhz(i)+zin*tzz(i))
            
          else
            
            oxx(i)=oxx(i)+step*hhx(i)
            oyy(i)=oyy(i)+step*hhy(i)
            ozz(i)=ozz(i)+step*hhz(i)
            
          endif
          
        enddo
        
      enddo

      return
      end subroutine turn_rigid_body
      
      end module optimiser_module
