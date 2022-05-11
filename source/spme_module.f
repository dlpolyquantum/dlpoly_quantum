      module spme_module

c***********************************************************************
c     
c     dl_poly module for defining spme arrays
c     
c     copyright - daresbury laboratory
c     author    - w. smith    sep 2003
c     
c***********************************************************************
      
      use config_module
      use error_module
      use setup_module
      use utility_module
      
      implicit none
      
      real(8), allocatable :: csp(:),qqc(:,:,:),ffttable(:)
      real(8), allocatable :: bspx(:,:), bspy(:,:), bspz(:,:)
      real(8), allocatable :: bsdx(:,:), bsdy(:,:), bsdz(:,:)
      integer, allocatable :: key1(:),key2(:),key3(:)
      complex(8), allocatable :: ww1(:), ww2(:), ww3(:)
      complex(8), allocatable :: qqq(:,:,:)
      complex(8), allocatable :: bscx(:), bscy(:),bscz(:)
CFFTW      pointer, save :: fplan, bplan
      
      save csp,qqc,qqq,ww1,ww2,ww3,bscx,bscy,bscz,ffttable
      save bspx,bspy,bspz,bsdx,bsdy,bsdz,key1,key2,key3
      
      contains
      
      subroutine alloc_spme_arrays(idnode,mxnode)
      
      implicit none
      
      integer, parameter :: nnn=9
      
      logical safe
      integer i,fail,idnode,mxnode
      dimension fail(nnn)

      safe=.true.

c     allocate arrays

      fail(:)=0
      
      allocate (csp(mxspl),ffttable(mxftab),stat=fail(1))
      allocate (bspx(mxspme,mxspl),bspy(mxspme,mxspl),stat=fail(2))
      allocate (bspz(mxspme,mxspl),bsdx(mxspme,mxspl),stat=fail(3))
      allocate (bsdy(mxspme,mxspl),bsdz(mxspme,mxspl),stat=fail(4))
      allocate (bscx(kmaxd),bscy(kmaxe),bscz(kmaxf),stat=fail(5))
      allocate (key1(kmaxd),key2(kmaxe),key3(kmaxf),stat=fail(6))
      allocate (ww1(kmaxd),ww2(kmaxe),ww3(kmaxf),stat=fail(7))
      allocate (qqc(kmaxd,kmaxe,kmaxf),stat=fail(8))
      allocate (qqq(kmaxd,kmaxe,kmaxf),stat=fail(9))
      
      if(any(fail.gt.0))safe=.false.      
      if(mxnode.gt.1)call gstate(safe)    
      if(.not.safe)call error(idnode,1750)
      
      end subroutine alloc_spme_arrays

      subroutine bspcoe(nospl,kmax1,kmax2,kmax3)

c**********************************************************************
c     
c     dl_poly subroutine to calculate B-spline coefficients for 
c     Euler exponential splines.
c     
c     copyright - daresbury laboratory 1998
c     author    - w. smith july 1998
c     
c***********************************************************************

      implicit none

      integer nospl,kmax1,kmax2,kmax3,k,i,j
      complex(8) ccc

c     calculate B-splines at knots

        csp(1)=0.d0
        csp(2)=1.d0
        
        do k=3,nospl
          
          csp(k)=0.d0
          
          do j=k,2,-1
            
            csp(j)=(dble(j-1)*csp(j)+dble(k-j+1)*csp(j-1))/dble(k-1)
            
          enddo
          
        enddo
        
c     calculate B-spline coefficients

      do i=0,kmax1-1

        ccc=(0.d0,0.d0)

        do k=0,nospl-2

          ccc=ccc+csp(k+2)*ww1(mod(i*k,kmax1)+1)

        enddo

        bscx(i+1)=ww1(mod(i*(nospl-1),kmax1)+1)/ccc

      enddo

      do i=0,kmax2-1

        ccc=(0.d0,0.d0)

        do k=0,nospl-2

          ccc=ccc+csp(k+2)*ww2(mod(i*k,kmax2)+1)

        enddo

        bscy(i+1)=ww2(mod(i*(nospl-1),kmax2)+1)/ccc

      enddo

      do i=0,kmax3-1

        ccc=(0.d0,0.d0)

        do k=0,nospl-2

          ccc=ccc+csp(k+2)*ww3(mod(i*k,kmax3)+1)

        enddo

        bscz(i+1)=ww3(mod(i*(nospl-1),kmax3)+1)/ccc

      enddo

      return
      end subroutine bspcoe

      subroutine bspgen(idnode,mxnode,natms,nospl,txx,tyy,tzz)

c***********************************************************************
c
c     dl_poly subroutine to calculate B-splines for SPME method
c
c     copyright - daresbury laboratory 1998
c     author    - w. smith july 1998
c     
c***********************************************************************

      implicit none

      integer nospl,natms,idnode,mxnode,iatm0,iatm1,i,j,k
      real(8) aaa,bbb,ccc,txx,tyy,tzz
      dimension txx(mxatms),tyy(mxatms),tzz(mxatms)
      
c     set up atoms numbers for nodes
      
      iatm0=(idnode*natms)/mxnode+1
      iatm1=((idnode+1)*natms)/mxnode
      
c     construct B-splines

      do i=iatm0,iatm1

        bsdx(i,1)=1.d0
        bsdy(i,1)=1.d0
        bsdz(i,1)=1.d0
        bsdx(i,2)=-1.d0
        bsdy(i,2)=-1.d0
        bsdz(i,2)=-1.d0
        bspx(i,1)=txx(i)-int(txx(i))
        bspy(i,1)=tyy(i)-int(tyy(i))
        bspz(i,1)=tzz(i)-int(tzz(i))
        bspx(i,2)=1.d0-txx(i)+int(txx(i))
        bspy(i,2)=1.d0-tyy(i)+int(tyy(i))
        bspz(i,2)=1.d0-tzz(i)+int(tzz(i))

      enddo
      
      do k=3,nospl
        
        do i=iatm0,iatm1

          bspx(i,k)=0.d0
          bspy(i,k)=0.d0
          bspz(i,k)=0.d0

        enddo
        
        do j=k,2,-1

          if(k.eq.nospl)then
            
            do i=iatm0,iatm1
              
              bsdx(i,j)=bspx(i,j)-bspx(i,j-1)
              bsdy(i,j)=bspy(i,j)-bspy(i,j-1)
              bsdz(i,j)=bspz(i,j)-bspz(i,j-1)
              
            enddo
            
          endif
          
          do i=iatm0,iatm1
            
            aaa=txx(i)+dble(j-1)-int(txx(i))
            bbb=tyy(i)+dble(j-1)-int(tyy(i))
            ccc=tzz(i)+dble(j-1)-int(tzz(i))
            bspx(i,j)=(aaa*bspx(i,j)+(dble(k)-aaa)*bspx(i,j-1))/
     x        dble(k-1)
            bspy(i,j)=(bbb*bspy(i,j)+(dble(k)-bbb)*bspy(i,j-1))/
     x        dble(k-1)
            bspz(i,j)=(ccc*bspz(i,j)+(dble(k)-ccc)*bspz(i,j-1))/
     x        dble(k-1)

          enddo
          
        enddo
        
        if(k.eq.nospl)then
          
          do i=iatm0,iatm1
            
            bsdx(i,1)=bspx(i,1)
            bsdy(i,1)=bspy(i,1)
            bsdz(i,1)=bspz(i,1)
            
          enddo
          
        endif
        
        do i=iatm0,iatm1

          bspx(i,1)=(txx(i)-int(txx(i)))*bspx(i,1)/dble(k-1)
          bspy(i,1)=(tyy(i)-int(tyy(i)))*bspy(i,1)/dble(k-1)
          bspz(i,1)=(tzz(i)-int(tzz(i)))*bspz(i,1)/dble(k-1)

        enddo
        
      enddo
      
      return
      end subroutine bspgen

      subroutine ewald_spme
     x (idnode,mxnode,natms,imcon,kmax1,kmax2,kmax3,nospl,
     x  engcpe,vircpe,alpha,volm,epsq)

c***********************************************************************
c     
c     dl_poly subroutine for calculating coulombic forces in a
c     periodic system using the smoothed particle mesh ewald method
c     due to Essmann et al J. Chem. Phys. 103 (1995) 8577.
c     
c     parallel replicated data version (part 1)
c     
c     copyright - daresbury laboratory 1998
c     author    - w. smith july 1998
c     additional FFT code - j. geronowicz sept 1999
c     
c     part 1 - reciprocal space terms (fourier part)
c     
c***********************************************************************

      implicit none

      logical newjob,lconsw
      
      integer idnode,mxnode,natms,imcon,kmax1,kmax2,kmax3,nospl
      integer npass,i,nnn,ipass,l,ll,k,kk,j,jj,fail,iatm0,iatm1
      real(8) engcpe,vircpe,alpha,volm,epsq,omg,bb1,bb2,bb3,qchg
      real(8) twopi,engsic,rvolm,ralph,shiftx,shifty,shiftz,det,qfix
      real(8) tx,ty,tz,rcpcut,rcpct2,rkx1,rky1,rkz1,rkx2,rky2,rkz2
      real(8) rkx3,rky3,rkz3,rksq,akv,eng1,den,scal1,tmp,rclprp
      real(8), allocatable :: txx(:),tyy(:),tzz(:)
CESSL      integer inc2,inc3
CSGIC      real(8) nauxfft(4)

      dimension omg(9),rclprp(10)
      complex(8) cpetot,vterm
      save newjob,engsic,qchg,iatm0,iatm1
      
      data newjob/.true./,fail/0/
CSGIC      data nauxfft/3,0,0,0/

c     allocate temporary arrays

      allocate (txx(mxatms),tyy(mxatms),tzz(mxatms),stat=fail)
      if(fail.ne.0)call error(idnode,1760)

      npass=1
      lconsw=.true.
      twopi=2.d0*pi
      
      if(newjob)then
        
        newjob=.false.
        
c     set up atoms numbers for nodes
        
        iatm0=(idnode*natms)/mxnode+1
        iatm1=((idnode+1)*natms)/mxnode
        
c     calculate self interaction correction and net system charge
        
        qchg=0.d0
        engsic=0.d0
        
        do i=iatm0,iatm1
          
          qchg=qchg+chge(i)
          engsic=engsic+chge(i)**2
          
        enddo
        
        if(mxnode.gt.1)then
          
          buffer(1)=qchg
          buffer(2)=engsic
          call gdsum(buffer(1),2,buffer(3))
          qchg  =buffer(1)
          engsic=buffer(2)
          
        endif
      
        engsic=-r4pie0/epsq*alpha*engsic/sqrpi

c     initialise the complex exponential arrays

CCRAY        call spl_cexp(kmax1,kmax2,kmax3,ww1,ww2,ww3)
CESSL        call spl_cexp(kmax1,kmax2,kmax3,ww1,ww2,ww3)
CFFTW        call spl_cexp(kmax1,kmax2,kmax3,ww1,ww2,ww3)
CSGIC        call spl_cexp(kmax1,kmax2,kmax3,ww1,ww2,ww3)

c     initialise the default fft routine

      call dlpfft3(1,1,kmax1,kmax2,kmax3,key1,key2,key3,
     x  ww1,ww2,ww3,qqq)
      
c     calculate B-spline coefficients

        call bspcoe(nospl,kmax1,kmax2,kmax3)
        
      endif
      
c     initialise coulombic potential energy
      
      engcpe=0.d0
      vircpe=0.d0

c     initalize stress tensor working arrays

      do i = 1,9
        omg(i) = 0.d0
      enddo

c     set working parameters

      rvolm=twopi/volm
      ralph=-0.25d0/alpha**2

c     set switch for TO, RD and HP boundary conditions

      if(imcon.eq.4.or.imcon.eq.5.or.imcon.eq.7) then

        npass=2
        lconsw=.false.
        rvolm=0.5d0*rvolm
        shiftx=0.5d0*dble(kmax1)
        shifty=0.5d0*dble(kmax2)
        shiftz=0.5d0*dble(kmax3)
        if(imcon.eq.7)shiftz=0.d0

      endif

c     convert cell coordinates to fractional coordinates

      call invert(cell,rcell,det)
      if(abs(det).lt.1.d-6)call error(idnode,120)
      
      do i=iatm0,iatm1
        
        txx(i)=dble(kmax1)*(rcell(1)*xxx(i)+rcell(4)*yyy(i)+
     x    rcell(7)*zzz(i)+0.5d0)
        tyy(i)=dble(kmax2)*(rcell(2)*xxx(i)+rcell(5)*yyy(i)+
     x    rcell(8)*zzz(i)+0.5d0)
        tzz(i)=dble(kmax3)*(rcell(3)*xxx(i)+rcell(6)*yyy(i)+
     x    rcell(9)*zzz(i)+0.5d0)

      enddo
      
c     construct B-splines for atoms
      
      call bspgen(idnode,mxnode,natms,nospl,txx,tyy,tzz)
      
c     zero 3D charge array
      
      nnn=kmaxd*kmaxe*kmaxf
      call set_block(nnn,0.d0,qqc)
      
c     construct 3D charge array
      
      do ipass=1,npass
        
        do i=iatm0,iatm1
          
          do l=1,nospl
            
            ll=int(tzz(i))-l+2
            if(ll.gt.kmax3)ll=1
            if(ll.lt.1)ll=ll+kmax3
            do k=1,nospl
              
              kk=int(tyy(i))-k+2
              if(kk.gt.kmax2)kk=1
              if(kk.lt.1)kk=kk+kmax2
              
              do j=1,nospl
                
                jj=int(txx(i))-j+2
                if(jj.gt.kmax1)jj=1
                if(jj.lt.1)jj=jj+kmax1
                
                qqc(jj,kk,ll)=qqc(jj,kk,ll)+
     x            chge(i)*bspx(i,j)*bspy(i,k)*bspz(i,l)
                
              enddo
              
            enddo
            
          enddo
          
        enddo
        
        if(.not.lconsw)then

          do i=iatm0,iatm1

            tx=txx(i)-shiftx
            ty=tyy(i)-shifty
            tz=tzz(i)-shiftz
            txx(i)=txx(i)-sign(shiftx,tx)
            tyy(i)=tyy(i)-sign(shifty,ty)
            tzz(i)=tzz(i)-sign(shiftz,tz)

          enddo

        endif
        
      enddo

c     global sum of charge array
        
      if(mxnode.gt.1) call gdsum(qqc,nnn,buffer)
        
c     load charge array into complex array for FFT

      call cpy_rtc(nnn,qqc,qqq)

c     calculate inverse 3D FFT of charge array (in place).
      
CFFTW      call fftwnd_f77_one(fplan,qqq,0)

CESSL      inc2=kmaxd
CESSL      inc3=kmaxd*kmaxe
CESSL      call dcft3(qqq,inc2,inc3,qqq,inc2,inc3,kmax1,kmax2,kmax3,
CESSL     x  -1,1.d0,buffer,mxbuff)

CSGIC      call zzfft3d( -1,kmax1,kmax2,kmax3,1.d0,qqq,kmaxd,kmaxe,
CSGIC     x              qqq,kmaxd,kmaxe,ffttable,buffer,nauxfft )

CCRAY      call ccfft3d( -1,kmax1,kmax2,kmax3,1.d0,qqq,kmaxd,kmaxe,
CCRAY     x              qqq,kmaxd,kmaxe,ffttable,buffer,nauxfft )
      
      call dlpfft3(0,1,kmax1,kmax2,kmax3,key1,key2,key3,
     x  ww1,ww2,ww3,qqq)
      
c     set reciprocal space cutoff

      call dcell(rcell,rclprp)

      rcpcut=0.5d0*min(dble(kmax1)*rclprp(7),dble(kmax2)*rclprp(8),
     x  dble(kmax3)*rclprp(9))
      rcpcut=rcpcut*1.05d0*twopi
      rcpct2=rcpcut**2
      
c     calculate convolution of charge array with gaussian function

      do l=1,kmax3

        ll=l-1
        if(l.gt.kmax3/2)ll=l-kmax3-1
        tmp=twopi*dble(ll)
        rkx1=tmp*rcell(3)
        rky1=tmp*rcell(6)
        rkz1=tmp*rcell(9)
        bb3=real(bscz(l)*conjg(bscz(l)))

        do k=1,kmax2

          kk=k-1
          if(k.gt.kmax2/2)kk=k-kmax2-1
          tmp=twopi*dble(kk)
          rkx2=rkx1+tmp*rcell(2)
          rky2=rky1+tmp*rcell(5)
          rkz2=rkz1+tmp*rcell(8)
          bb2=bb3*real(bscy(k)*conjg(bscy(k)))
          
          do j=1,kmax1

            jj=j-1
            if(j.gt.kmax1/2)jj=j-kmax1-1
            tmp=twopi*dble(jj)
            rkx3=rkx2+tmp*rcell(1)
            rky3=rky2+tmp*rcell(4)
            rkz3=rkz2+tmp*rcell(7)
            bb1=bb2*real(bscx(j)*conjg(bscx(j)))
                
            rksq=rkx3*rkx3+rky3*rky3+rkz3*rkz3

            if(rksq.gt.1.d-6.and.rksq.le.rcpct2)then

              vterm=bb1*exp(ralph*rksq)/rksq*qqq(j,k,l)
              akv=2.d0*(1.d0/rksq-ralph)*real(vterm*conjg(qqq(j,k,l)))
              omg(1)=omg(1)-rkx3*rkx3*akv
              omg(5)=omg(5)-rky3*rky3*akv
              omg(9)=omg(9)-rkz3*rkz3*akv
              omg(2)=omg(2)-rkx3*rky3*akv
              omg(3)=omg(3)-rkx3*rkz3*akv
              omg(6)=omg(6)-rky3*rkz3*akv
              qqq(j,k,l)=vterm

            else

              qqq(j,k,l)=(0.d0,0.d0)

            endif

          enddo

        enddo

      enddo
      
CFFTW      call fftwnd_f77_one(bplan,qqq,0)
CESSL      call dcft3(qqq,inc2,inc3,qqq,inc2,inc3,kmax1,kmax2,kmax3,
CESSL     x  1,1.d0,buffer,mxbuff)

CSGIC      call zzfft3d( 1,kmax1,kmax2,kmax3,1.d0,qqq,kmaxd,kmaxe,
CSGIC     x              qqq,kmaxd,kmaxe,ffttable,buffer,nauxfft )

CCRAY      call ccfft3d( 1,kmax1,kmax2,kmax3,1.d0,qqq,kmaxd,kmaxe,
CCRAY     x              qqq,kmaxd,kmaxe,ffttable,buffer,nauxfft )

      call dlpfft3(0,-1,kmax1,kmax2,kmax3,key1,key2,key3,
     x  ww1,ww2,ww3,qqq)

c     calculate atomic forces
      
      call spme_for
     x  (idnode,mxnode,nospl,natms,kmax1,kmax2,kmax3,rvolm,
     x  epsq,txx,tyy,tzz)

c     complete product of charge array and its gaussian convolution

      call ele_prd(nnn,qqq,qqc,qqq)
      
c     correction for charged systems

      qfix=-(0.5d0*pi*r4pie0/epsq)*((qchg/alpha)**2/volm)/
     x  dble(mxnode)
      
c     calculate total energy

      call scl_csum(nnn,cpetot,qqq)

      eng1=real(cpetot)
      den=1.d0/dble(npass)
      engcpe=engcpe+(den*rvolm*r4pie0*eng1/epsq+engsic)/
     x   dble(mxnode)+qfix

c     calculate stress tensor (symmetrical)

      scal1=den*rvolm*r4pie0/(epsq*dble(mxnode))
      stress(1) = stress(1)+scal1*(omg(1)+eng1)+qfix
      stress(2) = stress(2)+scal1*omg(2)
      stress(3) = stress(3)+scal1*omg(3)
      stress(4) = stress(4)+scal1*omg(2)
      stress(5) = stress(5)+scal1*(omg(5)+eng1)+qfix
      stress(6) = stress(6)+scal1*omg(6)
      stress(7) = stress(7)+scal1*omg(3)
      stress(8) = stress(8)+scal1*omg(6)
      stress(9) = stress(9)+scal1*(omg(9)+eng1)+qfix

c     virial term

      vircpe=vircpe-scal1*(omg(1)+omg(5)+omg(9)+3.d0*eng1)-3.d0*qfix
      
c     deallocate temporary arrays

      deallocate (txx,tyy,tzz,stat=fail)
      
      return
      end subroutine ewald_spme

      subroutine spme_for
     x  (idnode,mxnode,nospl,natms,kmax1,kmax2,kmax3,rvolm,
     x  epsq,txx,tyy,tzz)

c***********************************************************************
c     
c     dl_poly subroutine for calculating coulombic forces in a
c     periodic system using smoothed particle mesh ewald method
c     
c     parallel replicated data version (part 1)
c     
c     copyright - daresbury laboratory 1998
c     author    - w. smith oct 1998
c     
c     part 1 - reciprocal space terms (fourier part)
c
c***********************************************************************

      implicit none

      integer idnode,mxnode,nospl,natms,kmax1,kmax2,kmax3,i,ll
      integer iatm0,iatm1,kk,k,j,jj,l
      real(8) rvolm,epsq,txx,tyy,tzz,fff,fac,bdx,bdy,bdz
      real(8) det,qsum

      dimension txx(mxatms),tyy(mxatms),tzz(mxatms)
      dimension fff(3)

      fac=-2.d0*rvolm*r4pie0/epsq
      call invert(cell,rcell,det)

c     set up atom numbers for nodes

      iatm0 = (idnode*natms)/mxnode + 1
      iatm1 = ((idnode+1)*natms)/mxnode

c     calculate forces

      do i=iatm0,iatm1

        fxx(i)=0.d0
        fyy(i)=0.d0
        fzz(i)=0.d0

        do l=1,nospl
          
          ll=int(tzz(i))-l+2
          if(ll.gt.kmax3)ll=1
          if(ll.lt.1)ll=ll+kmax3
          
          do k=1,nospl
            
            kk=int(tyy(i))-k+2
            if(kk.gt.kmax2)kk=1
            if(kk.lt.1)kk=kk+kmax2
            
            do j=1,nospl
              
              jj=int(txx(i))-j+2
              if(jj.gt.kmax1)jj=1
              if(jj.lt.1)jj=jj+kmax1
              
              qsum=real(qqq(jj,kk,ll))
              bdx=qsum*bsdx(i,j)*bspy(i,k)*bspz(i,l)*dble(kmax1)
              bdy=qsum*bspx(i,j)*bsdy(i,k)*bspz(i,l)*dble(kmax2)
              bdz=qsum*bspx(i,j)*bspy(i,k)*bsdz(i,l)*dble(kmax3)
              
              fxx(i)=fxx(i)+fac*chge(i)*(bdx*rcell(1)+bdy*rcell(2)+
     x          bdz*rcell(3))
              fyy(i)=fyy(i)+fac*chge(i)*(bdx*rcell(4)+bdy*rcell(5)+
     x          bdz*rcell(6))
              fzz(i)=fzz(i)+fac*chge(i)*(bdx*rcell(7)+bdy*rcell(8)+
     x          bdz*rcell(9))
              
            enddo
            
          enddo
          
        enddo

      enddo

c     remove COM drift arising from SPME approximations

      fff(1)=0.d0
      fff(2)=0.d0
      fff(3)=0.d0

      do i=iatm0,iatm1

        fff(1)=fff(1)+fxx(i)
        fff(2)=fff(2)+fyy(i)
        fff(3)=fff(3)+fzz(i)

      enddo

      if(mxnode.gt.1)call gdsum(fff,3,buffer)

      fff(1)=fff(1)/dble(natms)
      fff(2)=fff(2)/dble(natms)
      fff(3)=fff(3)/dble(natms)

      do i=iatm0,iatm1

        fxx(i)=fxx(i)-fff(1)
        fyy(i)=fyy(i)-fff(2)
        fzz(i)=fzz(i)-fff(3)

      enddo
      
      return
      end subroutine spme_for

      subroutine dlpfft3
     x  (ind,isw,ndiv1,ndiv2,ndiv3,key1,key2,key3,ww1,ww2,ww3,aaa)

c***********************************************************************
c     
c     dl-poly 3D fast fourier transform routine (in place)
c     
c     copyright daresbury laboratory 1998
c     author w smith july 1998
c     
c***********************************************************************
      
      implicit none
      
      logical lkx,lky,lkz
      integer ind,isw,ndiv1,ndiv2,ndiv3,key1,key2,key3,i,idm,kkk
      integer nu1,nu2,nu3,iii,jjj,j,jj2,num,l,kk1,k12,k
      real(8) tpi,arg

      dimension key1(ndiv1),key2(ndiv2),key3(ndiv3)
      complex(8) ww1(ndiv1),ww2(ndiv2),ww3(ndiv3)
      complex(8) ttt,aaa(ndiv1,ndiv2,ndiv3)
      save nu1,nu2,nu3

      data tpi/6.283185307179586d0/
      
      if(ind.gt.0)then

c     check FFT array dimensions

        idm=1
        lkx=.true.
        lky=.true.
        lkz=.true.

        do i=1,30
          
          idm=2*idm

          if(idm.eq.ndiv1)then

            lkx=.false.
            nu1=i

          endif
          if(idm.eq.ndiv2)then

            lky=.false.
            nu2=i

          endif
          if(idm.eq.ndiv3)then

            lkz=.false.
            nu3=i

          endif
          
        enddo
        
        if(lkx.or.lky.or.lkz)then
          
          write(*,*)'error - FFT array not 2**N'
          stop
          
        endif
        
c     set reverse bit address arrays
        
        do kkk=1,ndiv1

          iii=0
          jjj=kkk-1

          do j=1,nu1

            jj2=jjj/2
            iii=2*(iii-jj2)+jjj
            jjj=jj2

          enddo

          key1(kkk)=iii+1

        enddo
        
        do kkk=1,ndiv2

          iii=0
          jjj=kkk-1

          do j=1,nu2

            jj2=jjj/2
            iii=2*(iii-jj2)+jjj
            jjj=jj2

          enddo

          key2(kkk)=iii+1

        enddo
        
        do kkk=1,ndiv3

          iii=0
          jjj=kkk-1

          do j=1,nu3

            jj2=jjj/2
            iii=2*(iii-jj2)+jjj
            jjj=jj2

          enddo

          key3(kkk)=iii+1

        enddo
        
c     initialise complex exponential factors
        
        ww1(1)=(1.d0,0.d0)

        do i=1,ndiv1/2

          arg=(tpi/dble(ndiv1))*dble(i)
          ww1(i+1)=cmplx(cos(arg),sin(arg),kind=8)
          ww1(ndiv1+1-i)=conjg(ww1(i+1))

        enddo
        
        ww2(1)=(1.d0,0.d0)

        do i=1,ndiv2/2

          arg=(tpi/dble(ndiv2))*dble(i)
          ww2(i+1)=cmplx(cos(arg),sin(arg),kind=8)
          ww2(ndiv2+1-i)=conjg(ww2(i+1))

        enddo
        
        ww3(1)=(1.d0,0.d0)

        do i=1,ndiv3/2

          arg=(tpi/dble(ndiv3))*dble(i)
          ww3(i+1)=cmplx(cos(arg),sin(arg),kind=8)
          ww3(ndiv3+1-i)=conjg(ww3(i+1))

        enddo
        
        return

      endif

c     take conjugate of exponentials if required
      
      if(isw.lt.0)then
        
        do i=1,ndiv1

          ww1(i)=conjg(ww1(i))

        enddo

        do i=1,ndiv2

          ww2(i)=conjg(ww2(i))

        enddo

        do i=1,ndiv3

          ww3(i)=conjg(ww3(i))

        enddo
        
      endif
      
c     perform fourier transform in X direction
      
      kkk=0
      num=ndiv1/2
      
      do l=1,nu1

        do while(kkk.lt.ndiv1)

          do i=1,num
            
            iii=key1(kkk/num+1)
            kk1=kkk+1
            k12=kk1+num
            
            do j=1,ndiv2
              
              do k=1,ndiv3
                
                ttt=aaa(k12,j,k)*ww1(iii)
                aaa(k12,j,k)=aaa(kk1,j,k)-ttt
                aaa(kk1,j,k)=aaa(kk1,j,k)+ttt
                
              enddo
              
            enddo
            
            kkk=kkk+1
            
          enddo
          
          kkk=kkk+num
          
        enddo
        
        kkk=0
        num=num/2

      enddo
      
c     unscramble the fft using bit address array
      
      do kkk=1,ndiv1

        iii=key1(kkk)

        if(iii.gt.kkk)then

          do j=1,ndiv2

            do k=1,ndiv3

              ttt=aaa(kkk,j,k)
              aaa(kkk,j,k)=aaa(iii,j,k)
              aaa(iii,j,k)=ttt

            enddo

          enddo

        endif

      enddo
      
c     perform fourier transform in Y direction
      
      kkk=0
      num=ndiv2/2

      do l=1,nu2

        do while(kkk.lt.ndiv2)

          do i=1,num
            
            iii=key2(kkk/num+1)
            kk1=kkk+1
            k12=kk1+num
            
            do j=1,ndiv1
              
              do k=1,ndiv3
                
                ttt=aaa(j,k12,k)*ww2(iii)
                aaa(j,k12,k)=aaa(j,kk1,k)-ttt
                aaa(j,kk1,k)=aaa(j,kk1,k)+ttt
                
              enddo
              
            enddo
            
            kkk=kkk+1
            
          enddo
          
          kkk=kkk+num
          
        enddo

        kkk=0
        num=num/2

      enddo
      
c     unscramble the fft using bit address array
      
      do kkk=1,ndiv2

        iii=key2(kkk)

        if(iii.gt.kkk)then

          do j=1,ndiv1

            do k=1,ndiv3

              ttt=aaa(j,kkk,k)
              aaa(j,kkk,k)=aaa(j,iii,k)
              aaa(j,iii,k)=ttt

            enddo

          enddo

        endif

      enddo
      
c     perform fourier transform in Z direction
      
      kkk=0
      num=ndiv3/2

      do l=1,nu3

        do while(kkk.lt.ndiv3)

          do i=1,num

            iii=key3(kkk/num+1)
            kk1=kkk+1
            k12=kk1+num
            
            do j=1,ndiv1
              
              do k=1,ndiv2
                
                ttt=aaa(j,k,k12)*ww3(iii)
                aaa(j,k,k12)=aaa(j,k,kk1)-ttt
                aaa(j,k,kk1)=aaa(j,k,kk1)+ttt
                
              enddo
              
            enddo
            
            kkk=kkk+1
            
          enddo
          
          kkk=kkk+num
          
        enddo

        kkk=0
        num=num/2

      enddo
      
c     unscramble the fft using bit address array
      
      do kkk=1,ndiv3

        iii=key3(kkk)

        if(iii.gt.kkk)then

          do j=1,ndiv1

            do k=1,ndiv2

              ttt=aaa(j,k,kkk)
              aaa(j,k,kkk)=aaa(j,k,iii)
              aaa(j,k,iii)=ttt

            enddo

          enddo

        endif

      enddo
      
c     restore exponentials to unconjugated values if necessary
      
      if(isw.lt.0)then
        
        do i=1,ndiv1

          ww1(i)=conjg(ww1(i))

        enddo
        
        do i=1,ndiv2

          ww2(i)=conjg(ww2(i))

        enddo
        
        do i=1,ndiv3

          ww3(i)=conjg(ww3(i))

        enddo
        
      endif
      
      return
      end subroutine dlpfft3

      subroutine spl_cexp(ndiv1,ndiv2,ndiv3,ww1,ww2,ww3)

c***********************************************************************
c     
c     dl-poly routine to create complex exponential arrays for
c     b-splines
c     
c     copyright daresbury laboratory 1998
c     author w smith oct 1998
c     
c***********************************************************************
      
      implicit none
      
      integer ndiv1,ndiv2,ndiv3,i
      real(8) tpi,arg
      complex(8) ww1(ndiv1),ww2(ndiv2),ww3(ndiv3)

      data tpi/6.283185307179586d0/

c     initialise complex exponential factors
      
      ww1(1)=(1.d0,0.d0)

      do i=1,ndiv1/2

        arg=(tpi/dble(ndiv1))*dble(i)
        ww1(i+1)=cmplx(cos(arg),sin(arg),kind=8)
        ww1(ndiv1+1-i)=conjg(ww1(i+1))

      enddo
      
      ww2(1)=(1.d0,0.d0)

      do i=1,ndiv2/2

        arg=(tpi/dble(ndiv2))*dble(i)
        ww2(i+1)=cmplx(cos(arg),sin(arg),kind=8)
        ww2(ndiv2+1-i)=conjg(ww2(i+1))

      enddo
      
      ww3(1)=(1.d0,0.d0)

      do i=1,ndiv3/2

        arg=(tpi/dble(ndiv3))*dble(i)
        ww3(i+1)=cmplx(cos(arg),sin(arg),kind=8)
        ww3(ndiv3+1-i)=conjg(ww3(i+1))

      enddo
      
      return
      end subroutine spl_cexp
      
      end module spme_module
