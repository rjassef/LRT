cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Program to calculate K corrections redshifts for testing the LRT
c     libraries.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      program main
      implicit real*8 (a-h,o-z)
      parameter (NCMAX=32)

      real*8 jy(11),ejy(11)
      integer jyuse(11)
      real*8 jymod(11),jycorr(11),comp(4)
      real*8 mag_kcorr(11),mag(11)

      real*8 jyzero(NCMAX),con(NCMAX),lbar(NCMAX)
      common /cal1/jyzero,con,lbar

      real*8 ebv,igm
      real*8 cov(4,4)

      character*200 line


c     Initialize the subroutines.
      call kcinit('bandmag.dat',1,1,1,1)

c     Open the files used.
      open(unit=11,file='gal_phot.dat',status='old')
      open(unit=12,file='kcorr.dat',status='unknown')

 200  read(11,'(a)')line
         if(line(1:1).eq.'#') goto 200
      continue
      read(line,*)ntarg,nchan

      print*
      print*,'Calculating K-Corrections for ',ntarg,' objects'


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

c     Calculate the photometric redshifts
         z0 = 0.d0

         call kca(jy,ejy,jyuse,z,z0,jymod,jycorr,comp,
     *        cov,ebv,igm,chi2,0)

c     Calculate the K-corrected and observed magntiudes.
         do j = 1,nchan
            if(jyuse(j).eq.1) then
               mag(j) = -2.5*log10(jy(j)/jyzero(j))
               mag_kcorr(j) = -2.5*log10(jy(j)/jyzero(j)/jycorr(j))
            else
               mag(j) = -99.
               mag_kcorr(j) = -99.
            endif
         enddo

         write(12,100)i,m,chi2,(mag(j),mag_kcorr(j),j=1,nchan),ebv,igm
 100     format(2i10,25E20.6)
        
      enddo

      close(11)
      close(12)

      print*,'Done'
      print*


      stop 
      end
