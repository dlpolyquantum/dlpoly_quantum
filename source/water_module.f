      module water_module
      
c***********************************************************************
c     
c     dl_poly quantum module for defining water model qtip4p/f
c     
c     authors    - M.R. Momeni & F.A. Shakib     
c     copyright  - M.R. Momeni & F.A. Shakib 2021
c
c     Method Development and Materials Simulation Laboratory
c     New Jersey Institute of Technology
c
c***********************************************************************
     
      use config_module
      use parse_module 
      use error_module
      use setup_module
      use site_module
      use utility_module
      
      implicit none

      integer nqt4
      integer, allocatable :: qt4index(:)      
      integer, allocatable :: numvoids(:),keyvod(:),lstvod(:,:)
      
      save numvoids,keyvod,lstvod,qt4index
      
      contains

      subroutine alloc_water_arrays(idnode,mxnode)

c***********************************************************************
c     
c     dl_poly quantum subroutine for defining water arrays
c
c     authors    - M.R. Momeni & F.A. Shakib     
c     copyright  - M.R. Momeni & F.A. Shakib 2021
c
c     Method Development and Materials Simulation Laboratory
c     
c***********************************************************************
      
      
      implicit none
      
      integer, parameter :: nnn=4
      
      logical safe
      integer i,fail,idnode,mxnode
      integer numatm
      dimension fail(nnn)
      
      safe=.true.
      numatm=mxatms*nbeads

c     allocate arrays

      fail(:)=0
     
c     allocation is based on arrays that are scanned in setup_module

      allocate (numvoids(mxtmls),stat=fail(1))
      allocate (keyvod(mxtvod),stat=fail(2))
      allocate (lstvod(mxtvod,3),stat=fail(3))
      allocate (qt4index(numatm),stat=fail(4))
      
      if(any(fail.gt.0))safe=.false.      
      if(mxnode.gt.1)call gstate(safe)    
      if(.not.safe)call error(idnode,1091)

c     initialise void arrays in case a molecular type does not have
c     void interactions.
     
      do i=1,mxtmls
         numvoids(i)=0
      enddo

c     initialise qt4index array for water model qtip4p/f
      
      do i=1,mxatms
         qt4index(i)=0
      enddo

      end subroutine alloc_water_arrays

      subroutine define_voids(safe,idnode,mxnode,itmols,nvoids)

c***********************************************************************
c     
c     dl_poly quantum subroutine for defining intramolecular interactions in
c     water molecules due to M-site to be excluded in ewald summation.
c     
c     authors    - M.R. Momeni & F.A. Shakib     
c     copyright  - M.R. Momeni & F.A. Shakib 2021
c
c     Method Development and Materials Simulation Laboratory
c     
c***********************************************************************

      implicit none
      
      logical safe
      character*8 keyword
      character*1 message(80)
      integer idnode,mxnode,itmols,nvoids,ntmp,ivoid,keytmp
      integer iatm1,iatm2,idum
      
      ntmp=intstr(record,lenrec,idum)
      numvoids(itmols)=numvoids(itmols)+ntmp

      if(idnode.eq.0)then
        write(nrite,"(/,1x,'number of void intramolecular interactions',
     x    7x,i10)")numvoids(itmols)
      endif
      
      do ivoid=1,ntmp

        nvoids=nvoids+1
        if(nvoids.gt.mxtvod)call error(idnode,3003)

c     read void parameters

        call getrec(safe,idnode,nfield)
        if(.not.safe)return
        
        call copystring(record,message,80)
        call lowcase(record,4)
        call getword(keyword,record,4,lenrec)

        if(keyword(1:4).eq.'qtp4')then
          keytmp=1
        else
          if(idnode.eq.0)write(nrite,*)message
          call error(idnode,444)
        endif

        iatm1=intstr(record,lenrec,idum)
        iatm2=intstr(record,lenrec,idum)
        
        keyvod(nvoids)=keytmp
        lstvod(nvoids,1)=iatm1
        lstvod(nvoids,2)=iatm2
        
      enddo
      
      return
      
      end subroutine define_voids

      subroutine water_index(idnode,mxnode,nbeads,ntpmls,iqt4)

c***********************************************************************
c     
c     dl_poly quantum subroutine for assigning inidces to O and H
c     atoms of water to be used for updating the coordinates
c     of the corresponding M-sites in water model qtip4p/f.
c     
c     authors    - M.R. Momeni & F.A. Shakib     
c     copyright  - M.R. Momeni & F.A. Shakib 2021
c
c     Method Development and Materials Simulation Laboratory
c     
c***********************************************************************

      implicit none
     
      integer i,j,n,itmols,idnode,mxnode
      integer, intent(in) :: nbeads
      integer, intent(in) :: ntpmls
      integer, intent(out) :: iqt4


c     check for errors

      if (numsit(ntpmls).ne.4) then
        write(nrite,*)'M-site coordinates should be entered as the
     x  4th atom in each water molecule.'
        stop
      endif

c     qtip4p/f water model indices

      iqt4=0
      nqt4=0
      do n=1,nbeads
       do itmols=1,ntpmls
        do i=1,nummols(itmols)
          do j=1,numsit(itmols)
            nqt4=nqt4+1
            if (itmols.eq.ntpmls.and.j.eq.numsit(ntpmls)) then
              iqt4=iqt4+1
              qt4index(iqt4)=nqt4
            endif
          enddo
        enddo
       enddo
      enddo  

      return

      end subroutine water_index

      subroutine qtip4pf(idnode,mxnode,imcon,nbeads,ntpmls,g_qt4f)

c***********************************************************************
c     
c     dl_poly quantum subroutine for updating the coordinates of M-sites 
c     in water model qtip4p/f based on the coordinates of the 
c     corresponding O and H atoms.
c     
c     authors    - M.R. Momeni & F.A. Shakib     
c     copyright  - M.R. Momeni & F.A. Shakib 2021
c
c     Method Development and Materials Simulation Laboratory
c     
c***********************************************************************
      
      implicit none
     
      integer i,n,m,mqt4,idnode,mxnode,imcon
      integer oxi,h1i,h2i
      integer, intent(in) :: nbeads
      integer, intent(in) :: ntpmls
      real(8) xho(2),yho(2),zho(2),g_qt4f_1
      real(8), intent(in) :: g_qt4f

      g_qt4f_1=0.5d0*(1.d0-g_qt4f)

      do n=1,nbeads

        m=n-1

        do i=1,nummols(ntpmls)

          mqt4=qt4index(i+(nummols(ntpmls)*m))
          oxi=qt4index(i+(nummols(ntpmls)*m))-3
          h1i=qt4index(i+(nummols(ntpmls)*m))-2
          h2i=qt4index(i+(nummols(ntpmls)*m))-1

          xho(1)=xxx(h1i)-xxx(oxi)
          yho(1)=yyy(h1i)-yyy(oxi)
          zho(1)=zzz(h1i)-zzz(oxi)
          xho(2)=xxx(h2i)-xxx(oxi)
          yho(2)=yyy(h2i)-yyy(oxi)
          zho(2)=zzz(h2i)-zzz(oxi)
          call images(imcon,0,1,2,cell,xho,yho,zho)
          xxx(mqt4)=xxx(oxi)+g_qt4f_1*(xho(1)+xho(2))
          yyy(mqt4)=yyy(oxi)+g_qt4f_1*(yho(1)+yho(2))
          zzz(mqt4)=zzz(oxi)+g_qt4f_1*(zho(1)+zho(2))

        enddo  
      enddo

      return

      end subroutine qtip4pf

      subroutine qt4_force_redist(idnode,mxnode,nbeads,ntpmls,g_qt4f)

c***********************************************************************
c     
c     dl_poly quantum subroutine for distributing the forces on M-sites in water
c     model qtip4p/f to the O and H atoms.
c     
c     authors    - M.R. Momeni & F.A. Shakib     
c     copyright  - M.R. Momeni & F.A. Shakib 2021
c
c     Method Development and Materials Simulation Laboratory
c     
c***********************************************************************
      
      implicit none
     
      integer i,m,n,mqt4,idnode,mxnode
      integer oxi,h1i,h2i
      integer, intent(in) :: nbeads
      integer, intent(in) :: ntpmls
      real(8), intent(in) :: g_qt4f
      real(8) g_qt4f_1

      g_qt4f_1=0.5d0*(1.d0-g_qt4f)

c     redistribute the force on the m-site to oxygen and hydrogen

      do n=1,nbeads

        m=n-1

        do i=1,nummols(ntpmls)
          mqt4=qt4index(i+(nummols(ntpmls)*m))
          oxi=qt4index(i+(nummols(ntpmls)*m))-3
          h1i=qt4index(i+(nummols(ntpmls)*m))-2
          h2i=qt4index(i+(nummols(ntpmls)*m))-1
          fxx(oxi)=fxx(oxi)+g_qt4f*fxx(mqt4)
          fyy(oxi)=fyy(oxi)+g_qt4f*fyy(mqt4)
          fzz(oxi)=fzz(oxi)+g_qt4f*fzz(mqt4)
          fxx(h1i)=fxx(h1i)+g_qt4f_1*fxx(mqt4)
          fyy(h1i)=fyy(h1i)+g_qt4f_1*fyy(mqt4)
          fzz(h1i)=fzz(h1i)+g_qt4f_1*fzz(mqt4)
          fxx(h2i)=fxx(h2i)+g_qt4f_1*fxx(mqt4)
          fyy(h2i)=fyy(h2i)+g_qt4f_1*fyy(mqt4)
          fzz(h2i)=fzz(h2i)+g_qt4f_1*fzz(mqt4)
        enddo

      enddo  

      return

      end subroutine qt4_force_redist

      end module water_module
