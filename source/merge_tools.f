      subroutine merge(idnode,mxnode,natms,nbuff,xxx,yyy,zzz,buffer)

c*********************************************************************
c     
c     dl_poly subroutine for merging coordinate arrays across
c     a number of processors
c     
c     parallel replicated data version
c
c     copyright - daresbury laboratory 1992
c     author    - w. smith november 1992.
c     MPI version - t. forester may 1995
c     CPP version - w.smith may 1995
c
c*********************************************************************

      use error_module

      implicit none

      integer idnode,mxnode,natms,nbuff,nsize,ierr,iatm1,iatm2
      integer j,i,k,jdnode,kdnode,katm0,katm1
      real*8 xxx(natms),yyy(natms),zzz(natms),buffer(nbuff)

      include "comms.inc"

      integer status(MPI_STATUS_SIZE), request

CMPIU      define MPI_SEND MPI_SEND_
CMPIU      define MPI_IRECV MPI_IRECV_
CMPIU      define MPI_WAIT MPI_WAIT_

c     check that buffer is large enough

      nsize=(natms+mxnode-1)/mxnode
      if(nbuff.lt.6*nsize)call error(idnode,47)

c     load initial transfer buffer

      j=0

c     set up this nodes atoms

      iatm1 = (idnode*natms)/mxnode + 1
      iatm2 = ((idnode+1)*natms)/mxnode

      do i=iatm1,iatm2

         buffer(j+1)=xxx(i)
         buffer(j+2)=yyy(i)
         buffer(j+3)=zzz(i)
         j=j+3

      enddo


      call gsync()

c     identity of neighbour node for systolic transfer

      jdnode=mod(idnode+1,mxnode)

      do k=1,mxnode-1

c     identity of node of origin of incoming data

         kdnode=mod(idnode+mxnode-k,mxnode)

c     identity of incoming  atoms

         katm0 = (kdnode*natms)/mxnode + 1
         katm1 = ((kdnode+1)*natms)/mxnode

c     systolic data pulse to transfer data

         call MPI_IRECV(buffer(3*nsize+1),3*nsize,MPI_DOUBLE_PRECISION,
     x        MPI_ANY_SOURCE,Merge_tag+k,MPI_COMM_WORLD,request,ierr)

         call MPI_SEND(buffer(1),3*nsize,MPI_DOUBLE_PRECISION,jdnode,
     x        Merge_tag+k,MPI_COMM_WORLD,ierr)

         call MPI_WAIT(request,status,ierr)

c     merge the incoming data into current arrays

         j=3*nsize

         do i=katm0,katm1

            xxx(i)=buffer(j+1)
            yyy(i)=buffer(j+2)
            zzz(i)=buffer(j+3)
            j=j+3

         enddo

c     shift new data to start of buffer

         do i=1,3*nsize

            buffer(i)=buffer(3*nsize+i)

         enddo

      enddo
      
      return
      end

      subroutine merge1(idnode,mxnode,natms,lstme,xxx,yyy,zzz,buffer)

c*********************************************************************
c     
c     dl_poly subroutine for merging together coordinate arrays
c     across a number of processors during rigid body algorithm
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1993
c     author    - t.forester  november 1993
c     systolic pulse version. T3D t.forester sept 1994
c     MPI version - t.forester may 1995
c     CPP version - w.smith may 1995
c     
c*********************************************************************
      
      use setup_module
      use error_module

      implicit none

      include "comms.inc"
      
      logical safe
      integer idnode,mxnode,natms,ierr,i,j,k,l,mxb,nout,nin
      integer nin1,jdnode,j1,kdnode
      integer lstme(natms)
      real*8 xxx(mxatms),yyy(mxatms),zzz(mxatms),buffer(mxbuff)
      
      integer status(MPI_STATUS_SIZE), request

CMPIU      define MPI_SEND MPI_SEND_
CMPIU      define MPI_IRECV MPI_IRECV_
CMPIU      define MPI_WAIT MPI_WAIT_

      
      safe =.true.
      
c     load up buffers
      
      j=1
      l=1
      do while(lstme(l).gt.0.and.l.le.natms)
        
        i=lstme(l)
        buffer(j+1)=dble(i)
        buffer(j+2)=xxx(i)
        buffer(j+3)=yyy(i)
        buffer(j+4)=zzz(i)
        j=j+4
        l=l+1
        
      enddo
      
c     length of message
      
      buffer(1) = dble(j)
      
c     array position for incoming messages
      
      mxb = mxbuff/2
      
c     load initial transfer buffer
      
      call gsync()
      
c     identity of neighbour node for systolic transfer
      
      jdnode=mod(idnode+1,mxnode)
      
      do k=1,mxnode-1
        
c     identity of node of origin of incoming data
        
        kdnode=mod(idnode+mxnode-k,mxnode)
        
c     out going message size
        
        nout = nint(buffer(1))
        
        call MPI_IRECV(nin,1,MPI_INTEGER,
     x    MPI_ANY_SOURCE,Merge1_tag+k,MPI_COMM_WORLD,request,ierr)
        
        call MPI_SEND(nout,1,MPI_INTEGER,jdnode,
     x    Merge1_tag+k,MPI_COMM_WORLD,ierr)
        
        call MPI_WAIT(request,status,ierr)
        
        call MPI_IRECV(buffer(mxb),nin,MPI_DOUBLE_PRECISION,
     x    MPI_ANY_SOURCE,Merge1_tag+k,MPI_COMM_WORLD,request,ierr)
        
        call MPI_SEND(buffer(1),nout,MPI_DOUBLE_PRECISION,jdnode,
     x    Merge1_tag+k,MPI_COMM_WORLD,ierr)
        
        call MPI_WAIT(request,status,ierr)
        
c     check buffer array not exceeded
        
        if(nin.gt.mxbuff-mxb) safe =.false.
        
c     position of first data element in incoming array
        
        nin1 = (nin-1)/4
        j = mxb+1
        
        do j1=1,nin1
          
          i = nint(buffer(j))
          xxx(i)=buffer(j+1)
          yyy(i)=buffer(j+2)
          zzz(i)=buffer(j+3)
          j=j+4
          
        enddo
        
c     shift new data to start of buffer
        
        do i=1,nin
          
          buffer(i)=buffer(mxb-1+i)
          
        enddo
        
      enddo
      
c     global check 
      
      if(mxnode.gt.1) call gstate(safe)
      if(.not.safe) call error(idnode,47)
      
      return
      end

      subroutine merge4(idnode,mxnode,ngrp,nbuff,q0,q1,q2,q3,buffer)

c*********************************************************************
c     
c     dl_poly subroutine for merging coordinate arrays across
c     a number of processors
c     
c     parallel replicated data version
c
c     copyright - daresbury laboratory 1994
c     author    - t.forester  february 1994
c     T3D version - sept 1994 t.forester
c     MPI version - t.forester may 1995
c     CPP version - w.smith may 1995
c
c*********************************************************************

      use error_module
      
      implicit none
      
      integer idnode,mxnode,ngrp,nbuff,ierr,nsize,i,j,igrp1,igrp2
      integer k,jdnode,kdnode,kgrp1,kgrp2
      real*8 q0(ngrp),q1(ngrp),q2(ngrp),q3(ngrp),buffer(nbuff)

      include "comms.inc"

      integer status(MPI_STATUS_SIZE), request

CMPIU      define MPI_SEND MPI_SEND_
CMPIU      define MPI_IRECV MPI_IRECV_
CMPIU      define MPI_WAIT MPI_WAIT_


c     check that buffer is large enough

      nsize=(ngrp+mxnode-1)/mxnode
      if(nbuff.lt.8*nsize)call error(idnode,47)

c     load initial transfer buffer

      j=0

      igrp1 = (idnode*ngrp)/mxnode+1
      igrp2 = ((idnode+1)*ngrp)/mxnode

      do i=igrp1,igrp2

         buffer(j+1)=q0(i)
         buffer(j+2)=q1(i)
         buffer(j+3)=q2(i)
         buffer(j+4)=q3(i)
         j=j+4

      enddo

      call gsync()

c     identity of neighbour node for systolic transfer

      jdnode=mod(idnode+1,mxnode)
     
      do k=1,mxnode-1

c     identity of node of origin of incoming data

         kdnode=mod(idnode+mxnode-k,mxnode)

c	identity of incoming groups 

         kgrp1 = (kdnode*ngrp)/mxnode+1
         kgrp2 = ((kdnode+1)*ngrp)/mxnode

         call MPI_IRECV(buffer(4*nsize+1),4*nsize,MPI_DOUBLE_PRECISION,
     x        MPI_ANY_SOURCE,Merge4_tag+k,MPI_COMM_WORLD,request,ierr)

         call MPI_SEND(buffer(1),4*nsize,MPI_DOUBLE_PRECISION,jdnode,
     x        Merge4_tag+k,MPI_COMM_WORLD,ierr)

         call MPI_WAIT(request,status,ierr)

c     merge the incoming data into current arrays

         j=4*nsize

         do i=kgrp1,kgrp2

            q0(i)=buffer(j+1)
            q1(i)=buffer(j+2)
            q2(i)=buffer(j+3)
            q3(i)=buffer(j+4)
            j=j+4

         enddo

c     shift new data to start of buffer

         do i=1,4*nsize

            buffer(i)=buffer(4*nsize+i)

         enddo

      enddo
      
      return
      end

      subroutine shlmerge(idnode,mxnode,ntshl)
      
c***********************************************************************
c     
c     dl_poly subroutine for merging core-shell velocity data
c     to restore data replication on all nodes
c     
c     copyright - daresbury laboratory 1993
c     author    - w. smith february 1993
c     MPI version - w. smith june 1995
c     CPP version - w. smith june 1995
c     
c***********************************************************************
     
      use setup_module
      use config_module
      use core_shell_module
      use error_module

      implicit none

      include "comms.inc"
      
      integer idnode,mxnode,ntshl,ierr,i,j,k,n,m,ishl1,ishl2,nsize
      integer jdnode,kshl1,kshl2,kdnode

      integer status(MPI_STATUS_SIZE), request

CMPIU      define MPI_SEND MPI_SEND_
CMPIU      define MPI_IRECV MPI_IRECV_
CMPIU      define MPI_WAIT MPI_WAIT_

c     check that buffer is large enough
      
      nsize=8*((ntshl+mxnode-1)/mxnode)
      
      if(mxbuff.lt.2*nsize)call error(idnode,425)

c     block indices

      ishl1 = (idnode*ntshl)/mxnode+1
      ishl2 = ((idnode+1)*ntshl)/mxnode

c     load initial transfer buffer
      
      n=0
      m=0
      
      do k=ishl1,ishl2

        m=m+1
        
c     indices of core and shell
        
        i=listshl(m,2)
        j=listshl(m,3)
        buffer(n+1)=dble(i)
        buffer(n+2)=dble(j)
        buffer(n+3)=vxx(i)
        buffer(n+4)=vyy(i)
        buffer(n+5)=vzz(i)
        buffer(n+6)=vxx(j)
        buffer(n+7)=vyy(j)
        buffer(n+8)=vzz(j)
        n=n+8

      enddo
      
      call gsync()
      
c     identity of neighbour node for systolic transfer
      
      jdnode=mod(idnode+1,mxnode)
      
      do k=1,mxnode-1

c     identity of node of origin of incoming data
        
        kdnode=mod(idnode+mxnode-k,mxnode)
        
c     systolic data pulse to transfer data
        
        call MPI_IRECV(buffer(nsize+1),nsize,MPI_DOUBLE_PRECISION,
     x    MPI_ANY_SOURCE,Shell_tag+k,MPI_COMM_WORLD,request,ierr)

        call MPI_SEND(buffer(1),nsize,MPI_DOUBLE_PRECISION,jdnode,
     x    Shell_tag+k,MPI_COMM_WORLD,ierr)

         call MPI_WAIT(request,status,ierr)

c     merge the incoming data into current arrays
        
        n=nsize

c     block indices

        kshl1 = (kdnode*ntshl)/mxnode+1
        kshl2 = ((kdnode+1)*ntshl)/mxnode

        do m=kshl1,kshl2
          
          i=nint(buffer(n+1))
          j=nint(buffer(n+2))
          
          vxx(i)=buffer(n+3)
          vyy(i)=buffer(n+4)
          vzz(i)=buffer(n+5)
          vxx(j)=buffer(n+6)
          vyy(j)=buffer(n+7)
          vzz(j)=buffer(n+8)

          n=n+8
          
        enddo
        
c     shift new data to start of buffer
        
        do i=1,nsize
          
          buffer(i)=buffer(nsize+i)
          
        enddo
        
      enddo
      
      return
      end

      subroutine shmove
     x     (idnode,mxnode,natms,lashap,lishap,xxt,yyt,zzt,
     x      txx,tyy,tzz,buffer)

c***********************************************************************
c     
c     dl_poly subroutine for passing coordinate updates between
c     nodes during the shake iteration cycle
c     
c     parallel replicated data algorithm 
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith august 1992.
c     MPI version - t.forester may 1995
c     CPP version - w.smith may 1995
c     
c***********************************************************************

      use setup_module

      implicit none

      include "comms.inc"

      integer idnode, mxnode, natms,ierr,i,k,j0,k0,n,jdnode,j
      integer lishap(mxlshp),lashap(mxnode)

      integer status(MPI_STATUS_SIZE), request

CMPIU      define MPI_SEND MPI_SEND_
CMPIU      define MPI_IRECV MPI_IRECV_
CMPIU      define MPI_WAIT MPI_WAIT_

      real*8 xxt(mxatms),yyt(mxatms),zzt(mxatms)
      real*8 txx(mxatms),tyy(mxatms),tzz(mxatms)
      real*8 buffer(mxbuff)

c     store increments to be transferred

      do i=1,natms

         txx(i)=xxt(i)
         tyy(i)=yyt(i)
         tzz(i)=zzt(i)

      enddo

c     transfer coordinate data to all nodes

      call gsync()

      do k=1,mxnode-1

         i=0
         j0=0
         if(k.gt.1)j0=lashap(k-1)

         do j=j0+1,lashap(k)

            buffer(i+1)=txx(lishap(j))
            buffer(i+2)=tyy(lishap(j))
            buffer(i+3)=tzz(lishap(j))
            i=i+3                             
                                               
         enddo

c     inter node communication

         k0=0

         if(k+1.lt.mxnode)k0=lashap(mxnode-k-1)
         n=3*(lashap(mxnode-k)-k0)
         jdnode=mod(idnode+k,mxnode)

c     check for zero length messages

         if(n.gt.0) call MPI_IRECV(buffer(i+1),n,MPI_DOUBLE_PRECISION,
     x        MPI_ANY_SOURCE,Shmove_tag+k,MPI_COMM_WORLD,request,ierr)

         if(i.gt.0) call MPI_SEND(buffer(1),i,MPI_DOUBLE_PRECISION,
     x        jdnode,Shmove_tag+k,MPI_COMM_WORLD,ierr)

         if(n.gt.0) call MPI_WAIT(request,status,ierr)

c     consolidate transferred data

         do j=k0+1,lashap(mxnode-k)

            xxt(lishap(j))=xxt(lishap(j))+buffer(i+1)
            yyt(lishap(j))=yyt(lishap(j))+buffer(i+2)
            zzt(lishap(j))=zzt(lishap(j))+buffer(i+3)
            i=i+3

         enddo
         
      enddo
      
      return
      end

      subroutine splice
     x      (idnode,natms,listme,listot,xxx,yyy,zzz,buffer)

c*********************************************************************
c     
c     dl_poly subroutine for splicing together coordinate arrays
c     across a number of processors during shake algorithm
c     
c     parallel replicated data version
c
c     copyright - daresbury laboratory 1993
c     author    - w. smith       march 1993
c
c     second version of splice
c
c*********************************************************************

      use setup_module
      use error_module
      
      implicit none

      integer idnode,natms,listme,listot,j,n3,i,lastot
      real*8 xxx,yyy,zzz,buffer

      dimension listme(mxatms),listot(mxatms)
      dimension xxx(natms),yyy(natms),zzz(natms)
      dimension buffer(mxbuff)

c     check buffer size

      if(mxbuff.lt.6*natms) call error(idnode,190)

c     load initial transfer buffers

      j=3*natms
      n3=3*natms

      do i=1,natms

         if(listot(i).gt.0)then

            if(listme(i).gt.0)then

               buffer(j+1)=xxx(i)
               buffer(j+2)=yyy(i)
               buffer(j+3)=zzz(i)

            else

               buffer(j+1)=0.d0
               buffer(j+2)=0.d0
               buffer(j+3)=0.d0

            endif

            j=j+3

         endif

      enddo

      lastot=j-n3

c     splice constraint coordinates

      if(lastot.gt.0) call gdsum(buffer(n3+1),lastot,buffer(1))

c     reconstitute coordinate arrays

      j=n3

      do i=1,natms

         if(listot(i).gt.0)then

            xxx(i)=buffer(j+1)/dble(listot(i))
            yyy(i)=buffer(j+2)/dble(listot(i))
            zzz(i)=buffer(j+3)/dble(listot(i))

            j=j+3

         endif

      enddo
      
      return
      end
