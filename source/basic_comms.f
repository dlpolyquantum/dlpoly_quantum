      subroutine initcomms()
      
c*********************************************************************
c     
c     communication harness initialisation
c     
c     copyright - daresbury laboratory
c     MPI version - t.forester may 1995
c     CPP version - w.smith may 1995
c     
c*********************************************************************
      
      implicit none
      
      include "comms.inc"
      
      integer ierr

CMPIU      define MPI_init MPI_init_

      call MPI_init(ierr)

      return
      end

      subroutine machine(idnode,mxnode)

c*********************************************************************
c     
c     dl_poly subroutine for obtaining charcteristics of
c     the computer on which the program is being run
c     
c     copyright daresbury laboratory 1992
c     author - w.smith july 1992
c     
c     MPI version - t.forester may 1995
c
c*********************************************************************

      implicit none

      integer idnode,mxnode,mynode,numnodes

      idnode=mynode()
      mxnode=numnodes()

      return
      end

      integer function mynode()

c*********************************************************************
c
c     routine to determine identity of processing node 
c
c     copyright - daresbury laboratory
c     MPI version - t.forester may 1995
c
c*********************************************************************

      implicit none

      include "comms.inc"

      integer ierr

CMPIU define MPI_COMM_RANK MPI_COMM_RANK_

      call MPI_COMM_RANK(MPI_COMM_WORLD, mynode ,ierr)

      return
      end

      integer function nodedim()

c*********************************************************************
c
c     calculate dimension of hypercube
c
c     copyright - daresbury laboratory
c     MPI version - t.forester may 1995
c
c*********************************************************************

      implicit none

      include "comms.inc"

      integer i,n,ierr,mxnode

CMPIU      define MPI_COMM_SIZE MPI_COMM_SIZE_

      call MPI_COMM_SIZE(MPI_COMM_WORLD, mxnode ,ierr)
      n=1
      nodedim = -1
      do i=0,16

         if(n.eq.mxnode)nodedim=i
         n=2*n

      enddo

      return
      end

      integer function numnodes()

c*********************************************************************
c
c     calculate number of nodes
c
c     copyright - daresbury laboratory
c     MPI version - t.forester may 1995
c
c*********************************************************************

      implicit none

      include "comms.inc"

      integer ierr

CMPIU      define MPI_COMM_SIZE MPI_COMM_SIZE_

      call MPI_COMM_SIZE(MPI_COMM_WORLD, numnodes, ierr)

      return
      end

      subroutine csend(tagmsg,buf,length,pe,idum)

c*********************************************************************
c
c     Intel-like  csend (double precision)
c
c     copyright - daresbury laboratory
c     MPI version - t.forester may 1995
c     CPP version - w.smith may 1995
c
c*********************************************************************

      implicit none

      include "comms.inc"

      integer tagmsg,length,pe,idum

      integer ierr
      real(8) buf(*)

CMPIU      define MPI_send MPI_send_

      call MPI_send(buf,length,MPI_DOUBLE_PRECISION,pe,tagmsg,
     x     MPI_COMM_WORLD,ierr)

      return
      end

      subroutine crecv(tagmsg,buf,length)

c*********************************************************************
c
c     Intel-like  crecv (double precision)
c
c     copyright - daresbury laboratory
c     MPI version - t.forester may 1995
c     CPP version - w.smith may 1995
c
c*********************************************************************

      implicit none

      include "comms.inc"

      integer tagmsg,length

      integer ierr
      integer status(MPI_STATUS_SIZE)
      real(8) buf(*)

CMPIU      define MPI_RECV MPI_RECV_

      call MPI_RECV(buf,length,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,
     x     tagmsg,MPI_COMM_WORLD,status,ierr)

      return 
      end

      subroutine gisum(aaa,nnn,bbb)

c***********************************************************************
c     
c     dl_poly global summation subroutine for hypercube - MPI version
c     integer version
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith march 1992.
c     MPI version - t.forester may 1995
c     CPP version - w.smith may 1995
c     
c***********************************************************************
      
      use setup_module

      implicit none

      integer nnn,i,ierror,iii,kk,k,k0,k1,k2,msg1,msg2
      integer aaa(nnn),bbb(nnn)

      include "comms.inc"

      integer status(MPI_STATUS_SIZE)

CMPIU      define MPI_allreduce MPI_allreduce_

      call MPI_allreduce(aaa,bbb,nnn,MPI_INTEGER,
     x  MPI_SUM,MPI_COMM_WORLD,ierror)

      do i = 1,nnn
        aaa(i) = bbb(i)
      enddo

      return
      end

      subroutine gdsum(aaa,nnn,bbb)

c***********************************************************************
c     
c     dl_poly global summation subroutine for MPI - hypercube assumed
c     double precision version
c     
c     copyright - daresbury laboratory 1995
c     author    - w. smith march 1992.
c     MPI version - t.forester may 1995
c     CPP version - w.smith may 1995
c     
c***********************************************************************

      implicit none

      integer nnn,i,iii,kk,k1,k2,ierror
      real(8) aaa(nnn),bbb(nnn)

      include "comms.inc"

      integer status(MPI_STATUS_SIZE)

CMPIU      define MPI_allreduce MPI_allreduce_

      call MPI_allreduce(aaa,bbb,nnn,MPI_DOUBLE_PRECISION,
     x  MPI_SUM,MPI_COMM_WORLD,ierror)

        do i = 1,nnn
          aaa(i) = bbb(i)
        enddo

      return
      end

      subroutine gimax(aaa,nnn,bbb)

c***********************************************************************
c     
c     dl_poly global maximum subroutine for hypercube - MPI version
c     integer version
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith march 1992.
c     MPI version - t.forester may 1995
c     CPP version - w.smith may 1995
c     
c***********************************************************************
      
      use setup_module

      implicit none

      integer nnn,i,iii,kk,k1,k2,k,k0msg1,msg2,ierror
      integer aaa(nnn),bbb(nnn)

      include "comms.inc"

      integer status(MPI_STATUS_SIZE)
CMPIU      define MPI_allreduce MPI_allreduce_
      
      call MPI_allreduce(aaa,bbb,nnn,MPI_INTEGER,
     x   MPI_MAX,MPI_COMM_WORLD,ierror)
      
      do i = 1,nnn
        aaa(i) = bbb(i)
      enddo

      return
      end

      subroutine gstate(check)

c***********************************************************************
c     
c     dl_poly global status subroutine : gisum version
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith       march 1992
c     MPI version -  t. forester may 1995
c     
c***********************************************************************


      implicit none

      logical check
      integer i(1),j(1)

      i(1) = 0
      if(.not.check) i(1) = 1

      call gisum(i,1,j)
      
      check = (i(1).eq.0)

      return
      end

      subroutine gsync()

c*********************************************************************
c     
c     barrier / synchronization routine
c
c     copyright - daresbury laboratory
c     MPI version - t.forester may 1995
c     CPP version - w.smith
c
c*********************************************************************

      implicit none

      integer ierr

      include "comms.inc"

CMPIU      define MPI_BARRIER MPI_BARRIER_

      call  MPI_BARRIER(MPI_COMM_WORLD,ierr)

      return
      end

      subroutine exitcomms()

c*********************************************************************
c
c     exitcomms: exit from communication harness
c
c     copyright - daresbury laboratory
c     MPI version - t.forester may 1995
c     CPP version - w.smith may 1995
c
c*********************************************************************

      implicit none

      include "comms.inc"

      integer ierr
CMPIU      define MPI_FINALIZE MPI_FINALIZE_

      call MPI_FINALIZE(ierr)
      call exit(0)

      return
      end
