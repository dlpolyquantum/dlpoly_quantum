      module utility_module

c***********************************************************************
c     
c     dl_poly module for utility subroutines and functions
c     
c     copyright - daresbury laboratory
c     author    - w. smith    aug 2006
c     
c***********************************************************************
      
      implicit none

      contains
      
      subroutine global_sum_forces(natms,mxnode,gxx,gyy,gzz)
      
c***********************************************************************
c     
c     dl_poly subroutine to perform global sum of atomic forces as
c     requred by replicated data strategy
c     
c     copyright - daresbury laboratory 
c     author    - w.smith december 2005
c     
c***********************************************************************
      
      use config_module
      
      implicit none
      
      integer natms,mxnode,i,j
      real(8) gxx(*),gyy(*),gzz(*)

      if(mxnode.gt.1) then
        
        j=0
        do i=1,natms
          
          buffer(j+1)=gxx(i)
          buffer(j+2)=gyy(i)
          buffer(j+3)=gzz(i)
          j=j+3
          
        enddo
        call gdsum(buffer(1),3*natms,buffer(3*natms+1))
        j=0
        do i=1,natms
          
          gxx(i)=buffer(j+1)
          gyy(i)=buffer(j+2)
          gzz(i)=buffer(j+3)
          j=j+3
          
        enddo
        
      endif
      
      return
      end subroutine global_sum_forces
      
      subroutine images
     x  (imcon,idnode,mxnode,natms,cell,xxx,yyy,zzz)
      
c***********************************************************************
c     
c     dl_poly subroutine for calculating the minimum image
c     of atom pairs within a specified MD cell
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith march 1992.
c     T3D optimised version. t.forester july 1994
c     
c     for
c     imcon=0 no boundary conditions apply
c     imcon=1 standard cubic boundaries apply
c     imcon=2 orthorhombic boundaries apply
c     imcon=3 parallelepiped boundaries apply
c     imcon=4 truncated octahedron boundaries apply
c     imcon=5 rhombic dodecahedron boundaries apply
c     imcon=6 x-y parallelogram boundary conditions : no periodicity in z
c     imcon=7 hexagonal prism boundaries apply
c     
c     note: in all cases the centre of the cell is at (0,0,0)
c     warning - replicated data version: does not re-merge 
c     coordinate arrays
c     
c***********************************************************************
      
      use error_module
      
      implicit none

      integer imcon,idnode,mxnode,natms,iatm1,iatm2,i
      real(8) cell,xxx,yyy,zzz,aaa,bbb,ccc,det,rt2,rt3,ssx
      real(8) ssy,ssz,ddd,xss,yss,zss,rcell

      dimension xxx(*),yyy(*),zzz(*)
      dimension cell(9),rcell(9)

      data rt2/1.41421356623d0/,rt3/1.7320508075d0/

      if(imcon.gt.0) then

c     block indices

        iatm1 = (idnode*natms)/mxnode+1
        iatm2 = ((idnode+1)*natms)/mxnode

      endif
      
      if(imcon.eq.1)then

c     standard cubic boundary conditions
        
        
        aaa=1.d0/cell(1)

        do i=iatm1,iatm2
          xxx(i)=xxx(i)-cell(1)*nint(aaa*xxx(i))
          yyy(i)=yyy(i)-cell(1)*nint(aaa*yyy(i))
          zzz(i)=zzz(i)-cell(1)*nint(aaa*zzz(i))
        enddo
        
      else if(imcon.eq.2)then

c     rectangular (slab) boundary conditions
        
        aaa=1.d0/cell(1)
        bbb=1.d0/cell(5)
        ccc=1.d0/cell(9)
        
        do i=iatm1,iatm2
          
          xxx(i)=xxx(i)-cell(1)*nint(aaa*xxx(i))
          yyy(i)=yyy(i)-cell(5)*nint(bbb*yyy(i))
          zzz(i)=zzz(i)-cell(9)*nint(ccc*zzz(i))
          
        enddo
        
      else if(imcon.eq.3)then

c     parallelepiped boundary conditions
        
        call invert(cell,rcell,det)
        
        do i=iatm1,iatm2
          
          ssx=(rcell(1)*xxx(i)+rcell(4)*yyy(i)+rcell(7)*zzz(i))
          ssy=(rcell(2)*xxx(i)+rcell(5)*yyy(i)+rcell(8)*zzz(i))
          ssz=(rcell(3)*xxx(i)+rcell(6)*yyy(i)+rcell(9)*zzz(i))
          
          xss=ssx-nint(ssx)
          yss=ssy-nint(ssy)
          zss=ssz-nint(ssz)
          
          xxx(i)=(cell(1)*xss+cell(4)*yss+cell(7)*zss)
          yyy(i)=(cell(2)*xss+cell(5)*yss+cell(8)*zss)
          zzz(i)=(cell(3)*xss+cell(6)*yss+cell(9)*zss)
          
        enddo
        
      else if(imcon.eq.4)then

c     truncated octahedral boundary conditions
        
        if(.not.(abs(cell(1)-cell(5)).lt.1.d-6.and.
     x    abs(cell(5)-cell(9)).lt.1.d-6)) call error(idnode,130)
        
        aaa=1.d0/cell(1)
        
        do i=iatm1,iatm2
          
          xxx(i)=xxx(i)-cell(1)*nint(aaa*xxx(i))
          yyy(i)=yyy(i)-cell(1)*nint(aaa*yyy(i))
          zzz(i)=zzz(i)-cell(1)*nint(aaa*zzz(i))
          
          if((abs(xxx(i))+abs(yyy(i))+abs(zzz(i))).ge.
     x      (0.75d0*cell(1)))then
            
            xxx(i)=xxx(i)-0.5d0*sign(cell(1),xxx(i))
            yyy(i)=yyy(i)-0.5d0*sign(cell(1),yyy(i))
            zzz(i)=zzz(i)-0.5d0*sign(cell(1),zzz(i))
            
          endif
          
        enddo
        
      else if(imcon.eq.5)then

c     rhombic dodecahedral boundary conditions
        
        if(.not.(abs(cell(1)-cell(5)).lt.1.d-6.and.
     x    abs(cell(9)-cell(1)*rt2).lt.1.d-6)) 
     x    call error(idnode,140)
        
        aaa=1.d0/cell(1)
        bbb=1.d0/cell(9)
        
        do i=iatm1,iatm2
          
          xxx(i)=xxx(i)-cell(1)*nint(aaa*xxx(i))
          yyy(i)=yyy(i)-cell(1)*nint(aaa*yyy(i))
          zzz(i)=zzz(i)-cell(9)*nint(bbb*zzz(i))
          
          if((abs(xxx(i))+abs(yyy(i))+abs(rt2*zzz(i))).ge.
     x      cell(1))then
            
            xxx(i)=xxx(i)-0.5d0*sign(cell(1),xxx(i))
            yyy(i)=yyy(i)-0.5d0*sign(cell(1),yyy(i))
            zzz(i)=zzz(i)-0.5d0*sign(cell(9),zzz(i))
            
          endif
          
        enddo
        
      else if(imcon.eq.6) then

c     x-y boundary conditions 

        det = cell(1)*cell(5) - cell(2)*cell(4)

        if(abs(det).lt.1.d-6)call error(idnode,120)
        
        det = 1.d0/det

        rcell(1) =  det*cell(5)
        rcell(2) = -det*cell(2)
        rcell(4) = -det*cell(4)
        rcell(5) =  det*cell(1)
        
        do i=iatm1,iatm2

          ssx = rcell(1)*xxx(i) + rcell(4)*yyy(i)
          ssy = rcell(2)*xxx(i) + rcell(5)*yyy(i)

          xss = ssx - nint(ssx)
          yss = ssy - nint(ssy)

          xxx(i)=cell(1)*xss + cell(4)*yss
          yyy(i)=cell(2)*xss + cell(5)*yss

        enddo

      else if(imcon.eq.7) then

c     hexagonal prism boundary conditions
        
        if(abs(cell(1)-rt3*cell(5)).ge.1.d-6)
     x    call error(idnode,135)
        
        aaa=cell(1)/(rt3*2.d0)
        bbb=cell(1)/rt3
        ccc=rt3/cell(1)
        ddd=1.d0/cell(9)
        
        do i=iatm1,iatm2
          
          yyy(i)=yyy(i)-bbb*nint(ccc*yyy(i))
          zzz(i)=zzz(i)-cell(9)*nint(ddd*zzz(i))
          
          if((abs(yyy(i))+abs(rt3*xxx(i))).ge.bbb)then
            
            xxx(i)=xxx(i)-rt3*sign(aaa,xxx(i))
            yyy(i)=yyy(i)-sign(aaa,yyy(i))
            
          endif
          
        enddo
        
      endif
      
      return
      end subroutine images

      subroutine imagesrev
     x  (imcon,cell,xxb,yyb,zzb,xxc,yyc,zzc)
      
c***********************************************************************
c     
c     dl_poly subroutine for calculating the minimum image
c     of bead atoms back from a specified MD cell to supercell
c     for HISTORYc in PIMD
c
c     for
c     imcon=1 standard cubic boundaries apply
c     imcon=3 parallelepiped boundaries apply
c     
c     note: in all cases the centre of the beads averages
c     
c     Dil Limbu
c     Method Development and Materials Simulation Laboratory
c***********************************************************************
      
      use error_module
      
      implicit none

      integer imcon,i
      real(8) cell,rcell,xxb,yyb,zzb
      real(8) aaa,bbb,ccc,det
      real(8) ssx,ssy,ssz,ddd,xss,yss,zss
      real(8) xxc,yyc,zzc

      dimension xxb(nbeads),yyb(nbeads),zzb(nbeads)
      dimension cell(9),rcell(9)


      xxc=sum(xxb(1:nbeads))/nbeads
      yyc=sum(yyb(1:nbeads))/nbeads
      zzc=sum(zzb(1:nbeads))/nbeads

      if(imcon.eq.1)then

c     standard cubic boundary conditions
        
        aaa=1.d0/cell(1)

c       do i=1,nbeads
        xxb(:)=xxb(:)-cell(1)*nint(aaa*(xxb(:)-xxc))
        yyb(:)=yyb(:)-cell(1)*nint(aaa*(yyb(:)-yyc))
        zzb(:)=zzb(:)-cell(1)*nint(aaa*(zzb(:)-zzc))
c       enddo

      else if(imcon.eq.3)then

c     parallelepiped boundary conditions
        
        call invert(cell,rcell,det)
        
        do i=1,nbeads
          
          ssx=(rcell(1)*xxb(i)+rcell(4)*yyb(i)+rcell(7)*zzb(i))
          ssy=(rcell(2)*xxb(i)+rcell(5)*yyb(i)+rcell(8)*zzb(i))
          ssz=(rcell(3)*xxb(i)+rcell(6)*yyb(i)+rcell(9)*zzb(i))
          
          xss=ssx-nint(ssx-xxc*rcell(1))
          yss=ssy-nint(ssy-yyc*rcell(5))
          zss=ssz-nint(ssz-zzc*rcell(9))
          
          xxb(i)=(cell(1)*xss+cell(4)*yss+cell(7)*zss)
          yyb(i)=(cell(2)*xss+cell(5)*yss+cell(8)*zss)
          zzb(i)=(cell(3)*xss+cell(6)*yss+cell(9)*zss)
          
        enddo

      endif

      xxc=sum(xxb(1:nbeads))/nbeads
      yyc=sum(yyb(1:nbeads))/nbeads
      zzc=sum(zzb(1:nbeads))/nbeads

      return
      end subroutine imagesrev


      subroutine config_write(fname,levcfg,imcon,natms,engcfg)
      
c***********************************************************************
c     
c     dl_poly subroutine for writing CONFIG files
c     
c     copyright - daresbury laboratory 
c     author    - w. smith aug 2007
c     
c***********************************************************************
      
      use config_module
      use setup_module
      
      implicit none
      
      character*6 fname
      
      integer i,natms,levcfg,imcon,nstep
      real(8) engcfg

      open(nconf,file=fname,form='formatted')
      
      write(nconf,'(80a1)') cfgname
      write(nconf,'(3i10,1p,g20.12)') levcfg,imcon,natms,engcfg
      if(imcon.gt.0) write(nconf,'(3f20.12)') cell
      
      do i=1,natms
        
        write(nconf,'(a8,i10)') atmnam(i),i
        write(nconf,'(3g20.10)') xxx(i),yyy(i),zzz(i)
        if(levcfg.gt.0)write(nconf,'(3g20.12)')
     x    vxx(i),vyy(i),vzz(i)
        if(levcfg.gt.1)write(nconf,'(3g20.12)') 
     x    fxx(i),fyy(i),fzz(i)
        
      enddo
      
      close (nconf)
      
      return
      end subroutine config_write
      
      subroutine bomb(idnode,nyr,nmn,ndy)

c***********************************************************************
c     
c     dl_poly subroutine to set an expiry date in a compiled program
c     
c     copyright - daresbury laboratory 
c     author    - w. smith    oct 2002
c
c***********************************************************************

      use setup_module

      implicit none

      logical safe
      integer info(8)
      character*12 dat,tim,zon
      integer idnode,nyr,nmn,ndy

      safe=.true.

      call date_and_time(dat,tim,zon,info)
      
      if(info(1).gt.nyr)then

        safe=.false.

      else if(info(1).eq.nyr)then

        if(info(2).gt.nmn)then

          safe=.false.

        else if(info(2).eq.nmn)then

          if(info(3).ge.ndy)safe=.false.

        endif

      endif

      if(.not.safe)then

        if(idnode.eq.0)write(nrite,'(a,/,a)')
     x    'THE EXPIRY DATE OF THIS EXECUTABLE HAS PASSED.',
     X    'PLEASE CONTACT W.SMITH@DL.AC.UK FOR A NEW LICENCE'

        call exitcomms()

      endif

      return
      end subroutine bomb

      subroutine cpy_rtc(nnn,aaa,bbb)

c**********************************************************************
c
c     dl_poly subroutine for copying a real array into a complex array
c     of the same dimension
c
c     copyright daresbury laboratory 1998
c     author w.smith oct 1998
c
c**********************************************************************

      implicit none

      integer i,nnn
      real(8) aaa(*)
      complex(8) bbb(*)

      do i=1,nnn

        bbb(i)=cmplx(aaa(i),0.d0,kind=8)

      enddo

      return
      end subroutine cpy_rtc

      function duni()

c*********************************************************************
c     
c     dl_poly random number generator based on the universal
c     random number generator of marsaglia, zaman and tsang
c     (stats and prob. lett. 8 (1990) 35-39.) it must be
c     called once to initialise parameters u,c,cd,cm
c     
c     copyright daresbury laboratory 1992
c     author -  w.smith         july 1992
c     
c*********************************************************************

      implicit none

      logical new
      integer ir,jr,i,j,k,l,m,ii,jj
      real(4) s,t,u,c,cd,cm,uni
      real(8) duni
      dimension u(97)
      save u,c,cd,cm,uni,ir,jr,new
      data new/.true./

      if(new)then

c     initial values of i,j,k must be in range 1 to 178 (not all 1)
c     initial value of l must be in range 0 to 168.

        i=12
        j=34
        k=56
        l=78
c     
        ir=97
        jr=33
        new=.false.

        do 200 ii=1,97
          s=0.0
          t=0.5
          do 100 jj=1,24
            m=mod(mod(i*j,179)*k,179)
            i=j
            j=k
            k=m
            l=mod(53*l+1,169)
            if(mod(l*m,64).ge.32)s=s+t
            t=0.5*t
  100     continue
          u(ii)=s
  200   continue
        c =  362436.0/16777216.0
        cd= 7654321.0/16777216.0
        cm=16777213.0/16777216.0
      else

c     calculate random number
        uni=u(ir)-u(jr)
        if(uni.lt.0.0)uni=uni+1.0
        u(ir)=uni
        ir=ir-1
        if(ir.eq.0)ir=97
        jr=jr-1
        if(jr.eq.0)jr=97
        c=c-cd
        if(c.lt.0.0)c=c+cm
        uni=uni-c
        if(uni.lt.0.0)uni=uni+1.0
        duni=dble(uni)
      endif
      
      return
      end function duni
      
      function puni(key,uuu)

c*********************************************************************
c     
c     dl_poly random number generator based on the universal
c     random number generator of marsaglia, zaman and tsang
c     (stats and prob. lett. 8 (1990) 35-39.) it must be
c     called once to initialise parameters u,c,cd,cm
c     
c     copyright daresbury laboratory 1992
c     author -  w.smith         july 1992
c     
c     parallel version - w.smith sep 2016
c
c     set key=0 to return a random number
c     set key=1 to initialise random number generator
c     set key=2 to reload control variables at restart
c     set key=3 to return control variables for saving
c     
c*********************************************************************
      
      implicit none
      
      integer ir,jr,i,j,k,l,m,ii,jj,key,idnode,mynode
      real(4) s,t,c,cd,cm,uni
      real(8) puni,uuu(102)
      real(4) u(97)
      save u,c,cd,cm,uni,ir,jr
      
      if(key.eq.1)then
        
c     initialise random number generator
c     initial values of i,j,k must be in range 1 to 178 (not all 1)
c     initial value of l must be in range 0 to 168.
        
        idnode=mynode()
        i=mod(11+idnode,178)+1
        j=mod(33+idnode,178)+1
        k=mod(55+idnode,178)+1
        l=mod(77+idnode,168)+1
        if(i.eq.1.and.j.eq.1.and.k.eq.1)then
          j=13
          k=131
        endif
        
        ir=97
        jr=33
        
        do ii=1,97
          s=0.0
          t=0.5
          do jj=1,24
            m=mod(mod(i*j,179)*k,179)
            i=j
            j=k
            k=m
            l=mod(53*l+1,169)
            if(mod(l*m,64).ge.32)s=s+t
            t=0.5*t
          enddo
          u(ii)=s
        enddo
        
        c =  362436.0/16777216.0
        cd= 7654321.0/16777216.0
        cm=16777213.0/16777216.0
        
      else if(key.eq.2)then
        
c     reload control variables at restart
        
        do ii=1,97
          u(ii)=uuu(ii)
        enddo
        c=uuu(98)
        cd=uuu(99)
        cm=uuu(100)
        ir=nint(uuu(101))
        jr=nint(uuu(102))
        
      else if(key.eq.3)then
        
c     return control variables for saving
        
        do ii=1,97
          uuu(ii)=u(ii)
        enddo
        uuu(98)=c
        uuu(99)=cd
        uuu(100)=cm
        uuu(101)=dble(ir)
        uuu(102)=dble(jr)
        
      else
        
c     generate random number
        
        uni=u(ir)-u(jr)
        if(uni.lt.0.0)uni=uni+1.0
        u(ir)=uni
        ir=ir-1
        if(ir.eq.0)ir=97
        jr=jr-1
        if(jr.eq.0)jr=97
        c=c-cd
        if(c.lt.0.0)c=c+cm
        uni=uni-c
        if(uni.lt.0.0)uni=uni+1.0
        puni=dble(uni)
        
      endif
      
      return
      end function puni
      
      subroutine ele_prd(nnn,aaa,bbb,ccc)

c**********************************************************************
c
c     dl_poly subroutine for element by element product of
c     a real array (bbb) and a complex array (ccc)
c
c     copyright daresbury laboratory 1998
c     author w.smith july 1998
c
c**********************************************************************

      implicit none

      integer i,nnn
      real(8) bbb(*)
      complex(8) aaa(*),ccc(*)

      do i=1,nnn

        aaa(i)=bbb(i)*ccc(i)

      enddo

      return
      end subroutine ele_prd

      subroutine gauss(natms,vxx,vyy,vzz)

c*********************************************************************
c     
c     dl_poly subroutine for constructing velocity arrays
c     with a gaussian distribution of unit variance.
c     
c     based on the Box-Muller method
c     
c     note - this version uses a universal random number 
c     generator, which generates pseudo-random numbers between
c     0 and 1. it is based on the algorithm of marsaglia, zaman
c     and tsang in: stats and prob. lett. 8 (1990) 35-39.
c     
c     copyright daresbury laboratory 2007
c     author - w. smith         nov  2007
c     
c*********************************************************************
      
      use setup_module
      
      implicit none

      integer natms,i
      real(8) vxx,vyy,vzz,rrr,rr1,rr2
      
      dimension vxx(natms),vyy(natms),vzz(natms)
      
c     initialise random number generator
      
      rrr=duni()
      
c     calculate gaussian random numbers
      
      do i=1,2*(natms/2),2
        
        rr1=sqrt(-2.d0*log(duni()))
        rr2=2.d0*pi*duni()
        vxx(i)=rr1*cos(rr2)
        vxx(i+1)=rr1*sin(rr2)

        rr1=sqrt(-2.d0*log(duni()))
        rr2=2.d0*pi*duni()
        vyy(i)=rr1*cos(rr2)
        vyy(i+1)=rr1*sin(rr2)

        rr1=sqrt(-2.d0*log(duni()))
        rr2=2.d0*pi*duni()
        vzz(i)=rr1*cos(rr2)
        vzz(i+1)=rr1*sin(rr2)
        
      enddo
      if(mod(natms,2).ne.0)then
        
        rr1=sqrt(-2.d0*log(duni()))
        rr2=2.d0*pi*duni()
        vxx(natms)=rr1*cos(rr2)
        vyy(natms)=rr1*sin(rr2)
        rr1=sqrt(-2.d0*log(duni()))
        rr2=2.d0*pi*duni()
        vzz(natms)=rr1*cos(rr2)
        
      endif
      
      return
      end subroutine gauss

      subroutine invert(a,b,d)

c***********************************************************************
c     
c     dl_poly subroutine to invert a 3 * 3 matrix using cofactors
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith       april 1992
c     
c***********************************************************************

      implicit none

      real(8) a,b,d,r

      dimension a(9),b(9)

c     calculate adjoint matrix
      b(1)=a(5)*a(9)-a(6)*a(8)
      b(2)=a(3)*a(8)-a(2)*a(9)
      b(3)=a(2)*a(6)-a(3)*a(5)
      b(4)=a(6)*a(7)-a(4)*a(9)
      b(5)=a(1)*a(9)-a(3)*a(7)
      b(6)=a(3)*a(4)-a(1)*a(6)
      b(7)=a(4)*a(8)-a(5)*a(7)
      b(8)=a(2)*a(7)-a(1)*a(8)
      b(9)=a(1)*a(5)-a(2)*a(4)

c     calculate determinant
      d=a(1)*b(1)+a(4)*b(2)+a(7)*b(3)
      r=0.d0
      if(abs(d).gt.0.d0)r=1.d0/d

c     complete inverse matrix
      b(1)=r*b(1)
      b(2)=r*b(2)
      b(3)=r*b(3)
      b(4)=r*b(4)
      b(5)=r*b(5)
      b(6)=r*b(6)
      b(7)=r*b(7)
      b(8)=r*b(8)
      b(9)=r*b(9)

      return
      end subroutine invert

      subroutine jacobi(a,v,n)

c***********************************************************************
c     
c     diagonalisation of real symmetric matices by jacobi method
c     
c     input parameters:
c     
c     a(n,n) is the matrix to be diagonalised
c     v(n,n) is the eigenvector matrix
c     n   is the dimension of the matrices
c     
c     jacobi processes lower triangle only (upper triangle unchanged)
c     
c     variable rho sets absolute tolerance on convergence
c     variable tes is a moving tolerance that diminishes
c     on each pass until at true convergence tes<rho
c     
c     author w.smith 1993
c     
c***********************************************************************

      implicit none

      logical pass
      integer n,i,j,k
      real(8) a,v,rho,tes,scl,v1,v2,v3,u,omg,s,c,tem

      dimension a(n,n),v(n,n)

      rho=1.0d-16
      tes=0.0d0
      scl=0.0d0

c     initialize eigenvectors

      do i=1,n
        do j=1,n
          v(i,j)=0.0d0
        enddo
        v(i,i)=1.0d0
      enddo

c     rescale matrix for optimal accuracy

      do i=1,n
        if(abs(a(i,i)).gt.scl)scl=abs(a(i,i))
      enddo
      do i=1,n
        do j=1,i
          a(i,j)=a(i,j)/scl
        enddo
      enddo

c     set initial value of moving tolerance

      do i=2,n
        do j=1,i-1
          tes=tes+2.0d0*a(i,j)*a(i,j)
        enddo
      enddo
      tes=sqrt(tes)

c     recycle until absolute tolerance satisfied

      do while(tes.gt.rho)

        tes=tes/dble(n)
        if(tes.lt.rho)tes=rho
        
c     jacobi diagonalisation
        
        pass=.true.
        
c     recycle until moving tolerance satisfied
        
        do while(pass)
          
          pass=.false.
          
          do i=2,n
            
            do j=1,i-1
              
              if(abs(a(i,j)).ge.tes)then
                pass=.true.
                v1=a(j,j)
                v2=a(i,j)
                v3=a(i,i)
                u=0.5d0*(v1-v3)
                if(abs(u).lt.rho)then
                  omg=-1.0d0
                else
                  omg=-v2/sqrt(v2*v2+u*u)
                  if(u.lt.0.0d0)omg=-omg
                endif
                s=omg/sqrt(2.0d0*(1.0d0+sqrt(1.0d0-omg*omg)))
                c=sqrt(1.0d0-s*s)
                do k=1,n
                  if(k.ge.i)then
                    tem=a(k,j)*c-a(k,i)*s
                    a(k,i)=a(k,j)*s+a(k,i)*c
                    a(k,j)=tem
                  else if(k.lt.j)then
                    tem=a(j,k)*c-a(i,k)*s
                    a(i,k)=a(j,k)*s+a(i,k)*c
                    a(j,k)=tem
                  else
                    tem=a(k,j)*c-a(i,k)*s
                    a(i,k)=a(k,j)*s+a(i,k)*c
                    a(k,j)=tem
                  endif
                  tem=v(k,j)*c-v(k,i)*s
                  v(k,i)=v(k,j)*s+v(k,i)*c
                  v(k,j)=tem
                enddo
                a(j,j)=v1*c*c+v3*s*s-2.0d0*v2*s*c
                a(i,i)=v1*s*s+v3*c*c+2.0d0*v2*s*c
                a(i,j)=(v1-v3)*s*c+v2*(c*c-s*s)
              endif
              
            enddo
            
          enddo
          
        enddo

      enddo

c     rescale matrix

      do i=1,n
        do j=1,i
          a(i,j)=scl*a(i,j)
        enddo
      enddo

      return
      end subroutine jacobi

      subroutine scl_csum(nnn,tot,aaa)

c**********************************************************************
c
c     dl_poly subroutine to calculate the scalar sum of the elements
c     of a complex array
c
c     copyright daresbury laboratory 1998
c     author w.smith july 1998
c
c**********************************************************************

      implicit none

      integer i,nnn
      complex(8) aaa(*),tot

      tot=(0.d0,0.d0)

      do i=1,nnn

        tot=tot+aaa(i)

      enddo

      return
      end subroutine scl_csum

      subroutine set_block(nnn,ccc,aaa)

c**********************************************************************
c
c     dl_poly subroutine to initialise an array to a single value
c
c     copyright daresbury laboratory 1998
c     author w.smith july 1998
c
c**********************************************************************

      implicit none

      integer i,nnn
      real(8) ccc,aaa(nnn)

      do i=1,nnn,2

        aaa(i)=ccc
        aaa(i+1)=ccc

      enddo
      
      return
      end subroutine set_block

      subroutine shellsort(n,list)

c***********************************************************************
c     
c     dlpoly shell sort routine. 
c     Sorts an array of integers into ascending order
c     
c     copyright daresbury laboratory 1993
c     author - t.forester   november 1993
c     
c***********************************************************************

      implicit none

      integer n,list,nn,nl,i,j,ix,imax

      dimension list(*)

c     set up sort

      if(n.gt.1) then

c     number of lists

        nl = n/2

c     iterate shell sort

        do while(nl.gt.0)

          do nn = 1,nl
            
c     begin insertion sort on nnth list
            
            do i = nn+nl,n,nl
              
              imax = list(i)
              ix = i
              
c     find location for insertion
              
              j = i
              do while(j.ge.nl+1)
                
                j = j-nl
                if (list(j).gt.imax) then
                  ix = j
                else
                  j =1
                endif
                
              enddo
              
c     insert in index array
              
              do j = i,ix+nl,-nl
                list(j) = list(j-nl)
              enddo
              
              list(ix) = imax
              
            enddo
            
          enddo
        
          nl = nl/2

        enddo
        
      endif

      return
      end subroutine shellsort

      subroutine fcap(lfcap,natms,fmax,temp)
      
c*********************************************************************
c     
c     DLPOLY routinue for limiting the absolute magnitude of
c     forces. Used in equilibration period only
c     
c     copyright daresbury laboratory 1993
c     
c     author -     t. forester march 1993
c     amended-     t. forester  sept 1994
c     
c*********************************************************************

      use config_module
      
      implicit none

      logical lfcap
      integer natms,i
      real(8) fmax,temp,fmax1,fmax2,fxc,fyc,fzc,fmod,fscale
      
      if(lfcap) then

c     maximum force permitted
        
        fmax1 = boltz*fmax*temp
        fmax2 = fmax1*fmax1

c     cap forces and conserve linear momentum
        
        fxc = 0.d0
        fyc = 0.d0
        fzc = 0.d0
        
        do i = 1,natms
          
          fmod = fxx(i)**2 + fyy(i)**2 + fzz(i)**2
          
          if(fmod.gt.fmax2) then
            
            fscale = sqrt(fmax2/fmod)
            
            fxx(i) = fxx(i)*fscale
            fyy(i) = fyy(i)*fscale
            fzz(i) = fzz(i)*fscale
            
          endif

c     accummulate forces - to check on momentum conservation
          
          fxc = fxc + fxx(i)
          fyc = fyc + fyy(i)
          fzc = fzc + fzz(i)
          
        enddo

c     ensure net forces sum to zero
        
        fxc = -fxc/dble(natms)
        fyc = -fyc/dble(natms)
        fzc = -fzc/dble(natms)

c     conserve momentum
        
        do i = 1,natms
          
          fxx(i) = fxx(i) + fxc
          fyy(i) = fyy(i) + fyc
          fzz(i) = fzz(i) + fzc
          
        enddo
        
      endif
      
      return
      end subroutine fcap

      subroutine freeze(natms)

c***********************************************************************
c     
c     dlpoly routine to quench forces and velocities on 'frozen' atoms
c     replicated data version - blocked data
c     
c     copyright daresbury laboratory 1994
c     author t.forester nov 1994
c     
c***********************************************************************

      use config_module

      implicit none

      integer natms,i

      do i = 1,natms
        
        if(lstfrz(i).ne.0) then
          
          vxx(i) = 0.d0
          vyy(i) = 0.d0
          vzz(i) = 0.d0
          fxx(i) = 0.d0
          fyy(i) = 0.d0
          fzz(i) = 0.d0
          
        endif
        
      enddo
      
      return
      end subroutine freeze

      subroutine mat_mul(aaa,bbb,ccc)

c***********************************************************************
c     
c     dlpoly utility to multiply 3x3 matrices
c
c     copyright daresbury laboratory
c     author      w.smith  oct  2005
c     
c**********************************************************************

      implicit none

      integer i
      real(8) aaa(9),bbb(9),ccc(9),tmp(9)

      tmp(1)=aaa(1)*bbb(1)+aaa(4)*bbb(2)+aaa(7)*bbb(3)
      tmp(2)=aaa(2)*bbb(1)+aaa(5)*bbb(2)+aaa(8)*bbb(3)
      tmp(3)=aaa(3)*bbb(1)+aaa(6)*bbb(2)+aaa(9)*bbb(3)

      tmp(4)=aaa(1)*bbb(4)+aaa(4)*bbb(5)+aaa(7)*bbb(6)
      tmp(5)=aaa(2)*bbb(4)+aaa(5)*bbb(5)+aaa(8)*bbb(6)
      tmp(6)=aaa(3)*bbb(4)+aaa(6)*bbb(5)+aaa(9)*bbb(6)

      tmp(7)=aaa(1)*bbb(7)+aaa(4)*bbb(8)+aaa(7)*bbb(9)
      tmp(8)=aaa(2)*bbb(7)+aaa(5)*bbb(8)+aaa(8)*bbb(9)
      tmp(9)=aaa(3)*bbb(7)+aaa(6)*bbb(8)+aaa(9)*bbb(9)
      
      do i=1,9
        ccc(i)=tmp(i)
      enddo
      
      return
      end subroutine mat_mul

      subroutine getrotmat(q0,q1,q2,q3,rot)
      
c***********************************************************************
c     
c     dlpoly utility to  construct rotation matrix
c     from quaternions using x convention for euler angles
c
c     copyright daresbury laboratory
c     author      w.smith   mar 2005
c     
c**********************************************************************

      implicit none
      
      real(8) q0,q1,q2,q3,rot(9)
      
      rot(1)=q0**2+q1**2-q2**2-q3**2
      rot(2)=2.d0*(q1*q2-q0*q3)
      rot(3)=2.d0*(q1*q3+q0*q2)
      rot(4)=2.d0*(q1*q2+q0*q3)
      rot(5)=q0**2-q1**2+q2**2-q3**2
      rot(6)=2.d0*(q2*q3-q0*q1)
      rot(7)=2.d0*(q1*q3-q0*q2)
      rot(8)=2.d0*(q2*q3+q0*q1)
      rot(9)=q0**2-q1**2-q2**2+q3**2
      
      return
      end subroutine getrotmat

      function sdot0(n,aaa,bbb)

c***********************************************************************
c     
c     dlpoly utility to calculate scalar product of two arrays
c
c     copyright daresbury laboratory
c     author      w.smith  july 2005
c     
c**********************************************************************

      implicit none

      integer n,i
      real(8) sdot0,aaa,bbb

      dimension aaa(*),bbb(*)

      sdot0=0.d0

      do i=1,n
        sdot0=sdot0+aaa(i)*bbb(i)
      enddo

      return
      end function sdot0

      function sdot1(natms,idnode,mxnode,aaa,bbb)

c***********************************************************************
c     
c     dlpoly utility to calculate scalar product of two arrays
c     distributed version
c     
c     copyright daresbury laboratory
c     author      w.smith  july 2005
c     
c**********************************************************************
      
      use config_module
      
      implicit none
      
      integer natms,idnode,mxnode,i,iatm0,iatm1
      real(8) sdot1,aaa,bbb
      
      dimension aaa(*),bbb(*)
      
c     assign block of atoms to processor

      iatm0=(idnode*natms)/mxnode + 1
      iatm1=((idnode+1)*natms)/mxnode
      
      sdot1=0.d0
      
      do i=iatm0,iatm1
        sdot1=sdot1+aaa(i)*bbb(i)
      enddo
      
      if(mxnode.gt.1)then
        buffer(1)=sdot1
        call gdsum(buffer(1),1,buffer(2))
        sdot1=buffer(1)
      endif
      
      return
      end function sdot1

      integer function loc2(i,j)

c*********************************************************************
c
c     calculates double index array minimum reference
c
c     copyright daresbury laboratory
c     author w.smith november 2005
c
c*********************************************************************
      
      integer i,j
      
      loc2=(max(i,j)*(max(i,j)-1))/2+min(i,j)
      
      return
      end function loc2

      integer function loc3(i,j,k)

c*********************************************************************
c
c     calculates triple index array minimum reference
c
c     copyright daresbury laboratory
c     author w.smith september 2008
c
c*********************************************************************
      
      integer i,j,k,a,b,c,u,v,w
      
      a=max(i,j)
      b=min(a,k)
      c=min(i,j)
      u=max(a,k)
      v=max(b,c)
      w=min(b,c)
      loc3=(u*(u*u-1))/6+(v*(v-1))/2+w
      
      return
      end function loc3

      integer function loc4(i,j,k,l)

c*********************************************************************
c
c     calculates quaduple index array minimum reference
c
c     copyright daresbury laboratory
c     author w.smith september 2008
c
c*********************************************************************
      
      integer i,j,k,l,a,b,c,d,e,f,t,u,v,w
      
      a=max(i,j)
      b=max(k,l)
      c=min(i,j)
      d=min(k,l)
      e=max(c,d)
      f=min(a,b)
      t=max(a,b)
      u=max(e,f)
      v=min(e,f)
      w=min(c,d)
      loc4=((((t+2)*t-1)*t-2)*t)/24+(u*(u*u-1))/6+(v*(v-1))/2+w
      
      return
      end function loc4

      character*3 function intstr3(nnn)

c*********************************************************************
c
c     converts a 3 digit integer to a string "001" etc.
c
c     copyright daresbury laboratory
c     author w.smith november 2005
c
c*********************************************************************

      implicit none

      integer nnn

      write(intstr3,'(i3.3)')nnn

      return
      end function intstr3
      

      subroutine traject
     x      (ltraj,idnode,imcon,istraj,keytrj,natms,
     x      nstraj,nstep,tstep,lpimd)
          
c***********************************************************************
c     
c     dl_poly subroutine for writing history file at selected
c     intervals in simulation
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith dec 1992.
c
c***********************************************************************

      use setup_module
      use config_module

      implicit none
      
      logical newjob,ltraj,lpimd
      integer idnode,imcon,istraj,keytrj,natms,nstraj,nstep,i
      real(8) tstep
      
      real(8) xxc,yyc,zzc
      integer :: natmc
      integer,parameter :: nhistc = 22

      save newjob
      data newjob/.true./
      
      natmc = natms/nbeads

      if(ltraj.and.idnode.eq.0)then
        
c     open the history file if new job or file closed
        
        if(newjob)then
          
          newjob = .false.

          open(nhist,file='HISTORY',position='append')
c    if(lpimd == .true.) then write CENTROID position of HISTORYc file
          if(lpimd)open(nhistc,file='HISTORYc',position='append')

        endif
        
        if(nstep.eq.nstraj.or.nstep.eq.istraj)then
          
          write(nhist,'(80a1)') cfgname
          write(nhist,'(3i10)') keytrj,imcon,natms

          if(lpimd)then
            write(nhistc,'(80a1)') cfgname
            write(nhistc,'(3i10)') keytrj,imcon,natmc
          endif

        endif
        
        if(mod(nstep-nstraj,istraj).eq.0)then
          
          write(nhist,'(a8,4i10,f12.6)') 'timestep',
     x         nstep,natms,keytrj,imcon,tstep

          if(imcon.gt.0) write(nhist,'(3g12.4)') cell
 
          do i = 1,natms

            write(nhist,'(a8,i10,2f12.6)')
     x        atmnam(i),i,weight(i),chge(i)
            write(nhist,'(1p,3e12.4)') xxx(i),yyy(i),zzz(i)
            if(keytrj.ge.1)then
              write(nhist,'(1p,3e12.4)') vxx(i),vyy(i),vzz(i)
            endif
            if(keytrj.ge.2)then
              write(nhist,'(1p,3e12.4)') fxx(i),fyy(i),fzz(i)
            endif

          enddo

          if(lpimd)then
            write(nhistc,'(a8,4i10,f12.6)') 'timestep',
     x         nstep,natmc,keytrj,imcon,tstep
 
            if(imcon.gt.0) write(nhistc,'(3g12.4)') cell
          
            do i = 1,natmc

c               call imagesrev
c     x         (imcon,cell,xxx(i:natms:natmc),yyy(i:natms:natmc),
c     x          zzz(i:natms:natmc),xxc,yyc,zzc)

              write(nhistc,'(a8,i10,2f12.6)')
     x          atmnam(i),i,weight(i),chge(i)
              write(nhistc,'(1p,3e12.4)')sum(xxx(i:natms:natmc))/nbeads,
     x                 sum(yyy(i:natms:natmc))/nbeads,
     x                 sum(zzz(i:natms:natmc))/nbeads
c              write(nhistc,'(1p,3e12.4)')xxc,yyc,zzc
              if(keytrj.ge.1)then
              write(nhistc,'(1p,3e12.4)')sum(vxx(i:natms:natmc))/nbeads,
     x                 sum(vyy(i:natms:natmc))/nbeads,
     x                 sum(vzz(i:natms:natmc))/nbeads
              endif
              if(keytrj.ge.2)then
              write(nhistc,'(1p,3e12.4)')sum(fxx(i:natms:natmc))/nbeads,
     x                 sum(fyy(i:natms:natmc))/nbeads,
     x                 sum(fzz(i:natms:natmc))/nbeads
              endif
            enddo
          endif

        endif

c     close history file at regular intervals
        
        if(.not.newjob.and.mod(nstep,ndump).eq.0)then
          
          close (nhist)
          if(lpimd)close (nhistc)
          newjob = .true.
          
        endif
        
      endif
      
      return
      end subroutine traject
      
      subroutine traject_u
     x     (ltraj,idnode,imcon,istraj,keytrj,natms,nstraj,nstep,tstep)
      
c***********************************************************************
c     
c     dl_poly subroutine for writing history file at selected
c     intervals in simulation
c     
c     Unformatted, double precision version
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith dec 1992.
c     
c***********************************************************************
      
      use setup_module
      use config_module
      
      implicit none
      
      logical newjob,ltraj
      integer idnode,imcon,istraj,keytrj,natms,nstraj,nstep,i
      real(8) tstep
      
      save newjob
      data newjob/.true./
      
      if(ltraj.and.idnode.eq.0)then
        
c     open the history file if new job or file closed
        
        if(newjob)  then
          
          newjob = .false.
          open(nhist,file='HISTORY',form='unformatted',
     x      position='append')
          
        endif
        if(nstep.eq.nstraj.or.nstep.eq.istraj)then
          
          write(nhist) cfgname
          write(nhist) dble(natms)
          write(nhist) (atmnam(i),i=1,natms)
          write(nhist) (weight(i),i=1,natms)
          write(nhist) (chge(i),i=1,natms)
          
        endif
        
        if(mod(nstep-nstraj,istraj).eq.0)then
          
          write(nhist)dble(nstep),dble(natms),dble(keytrj),
     x      dble(imcon),tstep
          
          if(imcon.gt.0) write(nhist) cell
          
          write(nhist) (xxx(i),i = 1,natms)
          write(nhist) (yyy(i),i = 1,natms)
          write(nhist) (zzz(i),i = 1,natms)
          
          if(keytrj.ge.1)then
            write(nhist) (vxx(i),i = 1,natms)
            write(nhist) (vyy(i),i = 1,natms)
            write(nhist) (vzz(i),i = 1,natms)
          endif
          if(keytrj.ge.2)then
            write(nhist) (fxx(i),i = 1,natms)
            write(nhist) (fyy(i),i = 1,natms)
            write(nhist) (fzz(i),i = 1,natms)
          endif
          
        endif
        
c     close history file at regular intervals
        
        if(.not.newjob.and.mod(nstep,ndump).eq.0)then
          
          close (nhist)
          newjob=.true.
          
        endif
        
      endif
      
      return
      end subroutine traject_u
      
      subroutine getcom_mol(ibeg,iend,imcon,idnode,mxnode,molmas,com)
      
c*********************************************************************
c     
c     dl_poly routine to calculate centre of mass of a molecule
c     specified between two atomic indices of the configuration
c     
c     copyright daresbury laboratory
c     author - w.smith june 2009
c     
c*********************************************************************
      
      use config_module
      
      implicit none
      
      integer ibeg,iend,imcon,i,idnode,mxnode,iatm0,iatm1,nmol
      integer j,fail
      real(8) molmas
      real(8) com(3)
      real(8), allocatable :: disx(:),disy(:),disz(:)
      data fail/0/
      
      nmol=iend-ibeg+1
      if(mxnode.eq.1.or.nmol.lt.mxnode)then
        
        iatm0=ibeg
        iatm1=iend
        
      else
        
        iatm0=(idnode*nmol)/mxnode+ibeg
        iatm1=((idnode+1)*nmol)/mxnode+ibeg-1
        nmol=(iend-ibeg+1)/mxnode+1
        
      endif
      
      allocate (disx(nmol),disy(nmol),disz(nmol),stat=fail)
      
      com(1)=0.d0
      com(2)=0.d0
      com(3)=0.d0
      molmas=0.d0
      
      j=0
      do i=iatm0,iatm1
        
        j=j+1
        disx(j)=xxx(i)-xxx(ibeg)
        disy(j)=yyy(i)-yyy(ibeg)
        disz(j)=zzz(i)-zzz(ibeg)
        
      enddo
      nmol=j
      
      call images(imcon,0,1,nmol,cell,disx,disy,disz)
      
      j=0
      do i=iatm0,iatm1
        
        j=j+1
        molmas=molmas+weight(i)
        com(1)=com(1)+weight(i)*(disx(j)+xxx(ibeg))
        com(2)=com(2)+weight(i)*(disy(j)+yyy(ibeg))
        com(3)=com(3)+weight(i)*(disz(j)+zzz(ibeg))

      enddo

      nmol=iend-ibeg+1
      if(mxnode.gt.1.and.nmol.ge.mxnode)then
        
        buffer(1)=com(1)
        buffer(2)=com(2)
        buffer(3)=com(3)
        buffer(4)=molmas
        call gdsum(buffer(1),4,buffer(5))
        com(1)=buffer(1)
        com(2)=buffer(2)
        com(3)=buffer(3)
        molmas=buffer(4)
        
      endif
      
      com(1)=com(1)/molmas
      com(2)=com(2)/molmas
      com(3)=com(3)/molmas

      call images(imcon,0,1,1,cell,com(1),com(2),com(3))

      deallocate(disx,disy,disz)
      
      return
      end subroutine getcom_mol

      subroutine timchk(ktim,time)
      
c***********************************************************************
c     
c     dlpoly timing routine for time elapsed in seconds
c     
c     copyright daresbury laboratory
c     author w.smith nov 2003
c
c***********************************************************************

      use setup_module

      implicit none

      logical init
      character*12 dat,tim,zon
      integer idnode,mynode,ktim,day
      real(8) time,told,tsum,tnow
      integer info(8)

      save init,idnode,told,tsum,day

      data init/.true./

   10 format(/,' time elapsed since job start = ',f15.3,' seconds',/)

      call date_and_time(dat,tim,zon,info)
      
      if(init)then

         tsum=0.d0
         time=0.d0
         day=info(3)
         idnode=mynode()
         told=3600.d0*dble(info(5))+60.d0*dble(info(6))+
     x         dble(info(7))+0.001d0*dble(info(8))
         init=.false.

      else

         tnow=3600.d0*dble(info(5))+60.d0*dble(info(6))+
     x         dble(info(7))+0.001d0*dble(info(8))
         if(day.ne.info(3))then
           told=told-86400.d0
           day=info(3)
         endif
         tsum=tsum+tnow-told
         told=tnow
         time=tsum

      endif

      if(ktim.gt.0.and.idnode.eq.0) write(nrite,10)time

      return
      end subroutine timchk
      
      subroutine get_prntime(hms,timelp,prntim)
      
c***********************************************************************
c     
c     dlpoly routine for casting cpu elapsed time into days, hours, 
c     minutes and seconds for printing (input timelp in seconds)
c     copyright daresbury laboratory
c     author w.smith may 2009
c
c***********************************************************************
      
      implicit none
      
      character*1 hms
      real(8) timelp,prntim
      
      if(timelp.ge.8.64d4)then
        hms='d'
        prntim=timelp/8.64d4
      elseif(timelp.ge.3.6d3)then
        hms='h'
        prntim=timelp/3.6d3
      elseif(timelp.ge.6.0d1)then
        hms='m'
        prntim=timelp/6.0d1
      else
        hms='s'
        prntim=timelp
      endif

      return
      end subroutine get_prntime
      
      subroutine get_simtime(dec,nstep,tstep,simtim)
      
c***********************************************************************
c     
c     dlpoly routine for casting simulation time into microseconds,
c     nanoseconds, picoseconds and femtoseconds for printing (input
c     tstep in picoseconds)
c     copyright daresbury laboratory
c     author w.smith may 2009
c     
c***********************************************************************
      
      implicit none
      
      character*1 dec
      integer nstep
      real(8) tmptim,simtim,tstep
      
      tmptim=tstep*dble(nstep)

      if(tmptim.ge.1.0d6)then
        dec='m'
        simtim=tmptim*1.0d-6
      elseif(tmptim.ge.1.0d3)then
        dec='n'
        simtim=tmptim*1.0d-3
      elseif(tmptim.ge.1.0d0)then
        dec='p'
        simtim=tmptim
      else
        dec='f'
        simtim=tmptim*1.0d3
      endif

      return
      end subroutine get_simtime

      end module utility_module
