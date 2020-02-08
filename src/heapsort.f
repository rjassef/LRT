c     Heap sort routine. Addapted from the generic code example in
c     wikipedia, http://en.wikipedia.org/wiki/Heapsort . Not very
c     polished nor optimized.

      subroutine heapsort(a,icount)
      implicit real*8 (a-h,o-z)

      real*8 a(*)

      call heapify(a,icount)

      iend = icount
      if(iend.le.1) goto 110
 100  continue
         call swap(a,iend,1)
         iend = iend-1
         call siftdown(a,1,iend)
         if(iend.le.1) goto 110
         goto 100
 110  continue

      return
      end

      subroutine heapify(a,icount)
      implicit real*8 (a-h,o-z)

      real*8 a(*)

      istart = icount/2
      
      if(istart.lt.1) goto 110
 100  continue
         call siftdown(a,istart,icount)
         istart = istart - 1
         if(istart.lt.1) goto 110
         goto 100
 110  continue

      return 
      end

      subroutine siftdown(a,istart,iend)
      implicit real*8 (a-h,o-z)

      real*8 a(*)

      iroot = istart

      if(iroot*2.gt.iend) goto 110
 100  continue
         ichild = iroot*2
         if(ichild.lt.iend.and.a(ichild).lt.a(ichild+1)) then
            ichild = ichild + 1
         endif
         if(a(iroot).lt.a(ichild)) then
            call swap(a,iroot,ichild)
            iroot = ichild
         else
            return
         endif
         if(iroot*2.gt.iend) goto 110
         goto 100
 110  continue

      return
      end

      subroutine swap(a,i,j)
      implicit real*8(a-h,o-z)
      
      real*8 a(*)

      aux = a(i)
      a(i) = a(j)
      a(j) = aux
      
      return
      end

