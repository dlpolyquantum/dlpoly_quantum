      subroutine passcon
     x     (lshmov,idnode,mxnode,natms,nscons,lashap,lishap,listme,
     x     listin,listot,listcon,lstfrz)
      
c*********************************************************************
c     
c     dl_poly subroutine for passing information about bond 
c     constraints between nodes
c     
c     parallel replicated data version assuming direct node-node
c     connection (i.e. this version may be intel specific)
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith august 1992.
c     MPI version t.forester may 1995
c     CPP version - w.smith may 1995
c     
c***********************************************************************

      use setup_module
      use error_module

      implicit none

      include "comms.inc"

      logical safe,lshmov
      integer idnode,mxnode,natms,nscons,lashap,lishap,listme,ierr
      integer listin,listot,listcon,lstfrz,i,j,k,jdnode,idum

      dimension listme(mxatms),listin(mxatms),listot(mxatms)
      dimension lishap(mxlshp),lashap(mxnode),listcon(mxcons,3)
      dimension lstfrz(mxatms)

      integer status(MPI_STATUS_SIZE), request

CMPIU define MPI_SEND MPI_SEND_
CMPIU define MPI_IRECV MPI_IRECV_
CMPIU define MPI_WAIT MPI_WAIT_

      safe=.true.

      do i=1,natms
         
         listme(i)=0
         
      enddo
      
      do k=1,nscons
         
         i=listcon(k,2)
         j=listcon(k,3)
         listme(i)=listme(i)+1
         listme(j)=listme(j)+1
         
      enddo
      
      if(mxnode.gt.1)then
         
         j=0
         call gsync()
         do k=1,mxnode-1
            
            jdnode=mod(idnode+mxnode-k,mxnode)

            call MPI_IRECV(listin,natms,MPI_INTEGER,
     x        MPI_ANY_SOURCE,Passcon_tag+k,MPI_COMM_WORLD,request,ierr)
            
            call MPI_SEND(listme,natms,MPI_INTEGER,jdnode,
     x           Passcon_tag+k,MPI_COMM_WORLD,ierr)

            call MPI_WAIT(request,status,ierr)

            do i=1,natms
               
               if((listme(i).gt.0).and.(listin(i).gt.0.and.
     x              lstfrz(i).eq.0))then
                  
                  j=j+1
                  if(j.gt.mxlshp)then

                     safe=.false.

                  else

                     lishap(j)=i

                  endif
                  
               endif
               
            enddo
            
            lashap(k)=j
            
         enddo
         
      endif

c     check for global error condition

      if(mxnode.gt.1) call gstate(safe)

      if(.not.safe)call error(idnode,103)

      if(mxnode.gt.1) then
         call gisum(j,1,idum)
         if(idnode.eq.0) write(nrite,'(/,a,14x,i10)')
     x     ' shared atoms from passcon',j/2
         lshmov = (j.gt.0)
      endif

c     keep record of all atoms subject to constraints
      
      do i=1,natms
         
         if(listme(i).gt.0)then
            
            listot(i)=1
            
         else
            
            listot(i)=0
            
         endif
         
      enddo
      
      if(mxnode.gt.1)call gisum(listot,natms,listin)
      
      return
      end

      subroutine passpmf
     x  (idnode,mxnode,natms,nspmf,listpm,listin,lstpmt,lstpmf,npmf)

c*********************************************************************
c     
c     dl_poly subroutine for passing information about PMF
c     constraints between nodes
c     
c     parallel replicated data version assuming direct node-node
c     connection (i.e. this version may be intel specific)
c     
c     copyright - daresbury laboratory 1995
c     author    - t.forester august 1995.
c     
c***********************************************************************

      use setup_module
      use error_module

      implicit none

      integer idnode,mxnode,natms,nspmf,listpm,listin,lstpmt,lstpmf
      integer npmf,i,j,k

      dimension listpm(mxatms),listin(mxatms),lstpmt(mxatms)
      dimension lstpmf(mxspmf,mspmf),npmf(2)

      do i=1,mxatms
        
        listpm(i)=0
        lstpmt(i)=0
        
      enddo
      
      do k=1,nspmf
        
        do j=1,npmf(1)+npmf(2)

          i=lstpmf(j,k)
          listpm(i)=1
          lstpmt(i)=1
          
        enddo

      enddo
      
c     keep record of all atoms subject to pmf constraints
      
      if(mxnode.gt.1)call gisum(lstpmt,natms,listin)
      
      return
      end

      subroutine passquat
     x  (lcnb,idnode,mxnode,natms,ngrp,nscons,ntpmls,listin,
     x  listcon,lstrgd,lstout,lstcsit,lstgtp,nummols,numgrp,numgsit)

c*********************************************************************
c     
c     dl_poly subroutine for passing information about rigid body 
c     atoms involved in bond constraints between nodes
c     
c     parallel replicated data version assuming direct node-node
c     connection
c     
c     copyright - daresbury laboratory 1995
c     author    - t. forester december 1995.
c     
c***********************************************************************
      
      use setup_module
      use error_module

      implicit none

      include "comms.inc"

      logical lcnb,safe
      integer idnode,mxnode,natms,ngrp,nscons,ntpmls,listin
      integer listcon,lstrgd,lstout,lstcsit,lstgtp,nummols,numgrp
      integer numgsit,igrp1,igrp2,i,jr,igrp,itmols,imols,lgrp,id
      integer jj,ik,j,k
      
      dimension listin(mxatms)
      dimension listcon(mxcons,3),lstcsit(2*mxcons)
      dimension lstout(mxatms),lstrgd(mxgatm)
      dimension nummols(mxtmls),numgrp(mxtmls),numgsit(mxungp)
      dimension lstgtp(mxgrp)

      integer status(MPI_STATUS_SIZE)
      
c     block indices for groups
      
      igrp1 = (idnode*ngrp)/mxnode + 1
      igrp2 = ((idnode+1)*ngrp)/mxnode
      
c     locate site indices of atoms in constraints

      do i = 1,natms
        listin(i) = 0
      enddo

c     loop over molecule types

      jr = 0 
      igrp = 0
      do itmols=1,ntpmls

c     loop over molecules in system
        
        do imols=1,nummols(itmols)

c     construct rigid body site list: each processor has a different copy
          
          do lgrp=1,numgrp(itmols)
            
            igrp=igrp+1
            
            if((igrp.ge.igrp1).and.(igrp.le.igrp2)) then
                
              id = lstgtp(igrp)
              do jj = 1,numgsit(id)
                  
                jr = jr +1
                i = lstrgd(jr)
                listin(i) = jj

              enddo
            endif
          enddo
        enddo
      enddo

      if(mxnode.gt.1) call gisum(listin,natms,lstout)

      safe = .true.
      ik = 0
      do k=1,nscons
        
        i=listcon(k,2)
        j=listcon(k,3)

        if(listin(i).ne.0) then
          ik = ik + 1
          lstcsit(ik) = listin(i)
          safe = .false.
        endif

        if(listin(j).ne.0) then
          ik = ik + 1
          lstcsit(ik) = listin(j)
          safe = .false.
        endif

      enddo

c     lcnb flags bodies connected by constraints

      if(mxnode.gt.1) call gstate(safe)
      lcnb = (.not.safe)
      
      return
      end



