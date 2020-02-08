cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Program to calculate photometric redshifts for testing the LRT
c     libraries.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      program main
      implicit real*8 (a-h,o-z)

      real*8 jy(11),ejy(11),comp(4)
      integer jyuse(11)

      character*200 line

c     Initialize the subroutines.
      call pzinit('bandmag.dat',1,1,0,0.d0,1.5d0,0.01d0,1)

c     Open the files used.
      open(unit=11,file='gal_phot.dat',status='old')
      open(unit=12,file='zphot.dat',status='unknown')

 200  read(11,'(a)')line
         if(line(1:1).eq.'#') goto 200
      continue
      read(line,*)ntarg,nchan
         

      print*
      print*,'Calculating photometric redshifts for ',ntarg,' objects'


c     Start the main cycle.
      do i=1,ntarg

c     Read the 6'' aperture Photometry File
         read(11,*)z
         read(11,*)(jy(j),j=1,nchan)
         read(11,*)(ejy(j),j=1,nchan)
         read(11,*)(jyuse(j),j=1,nchan)

c     Calculate number of usable bands
         m = 0
         do j = 1,nchan
            m = m + jyuse(j)
         enddo

c     Calculate the photometric redshifts.
         call pza(jy,ejy,jyuse,zp,chigal,chinop,0,i)

         write(12,100)i,m,zp,z,chigal,chinop
 100     format(2i10,4E20.6)
        
      enddo

      close(11)
      close(12)

      print*,'Done'
      print*


      stop 
      end

