cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Program to calculate K corrections redshifts for testing the LRT
c     libraries.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


        program main
        implicit real*8 (a-h,o-z)
        parameter (NCMAX=40)

        real*8 jy(11),ejy(11)
        integer jyuse(11)
        real*8 jymod(11),jycorr(11),comp(4)
        real*8 mag_kcorr(11),mag(11)

        real*8 jyzero(NCMAX),con(NCMAX),lbar(NCMAX)
        common /cal1/jyzero,con,lbar

        real*8 ebv,igm
        real*8 cov(4,4)

        character*200 line


c     Initialize the subroutines. Will only fit as Main Sequence for this test.
        call starinit('bandmag_stars.dat',1,1)

c     Open the files used.
        open(unit=11,file='fake_stellar_catalog.dat',status='old')
        open(unit=12,file='star_model.dat',status='unknown')

200     read(11,'(a)')line
        if(line(1:1).eq.'#') goto 200
            continue
        read(line,*)ntarg,nchan

        print*
        print*,'Calculating stellar fits for ',ntarg,' objects'


c     Start the main cycle.
        do i=1,ntarg

c     Read the 6'' aperture Photometry File
            read(11,*)(jy(j),j=1,nchan)

            do j=1,nchan
                ejy(j) = 0.1d0 !*jy(j)
                if(jy(j).gt.0.d0) then
                    jyuse(j) = 1
                else
                    jyuse(j) = 0
                endif
            enddo

c     Calculate number of usable bands
            m = 0
            do j = 1,nchan
                m = m + jyuse(j)
            enddo

            call stf(jy,ejy,jyuse,jymod,comp,ns_best,chi2,0)

            write(12,100)i,m,ns_best,chi2,(comp(l),l=1,2)
100         format(3i10,25E20.6)
      
        enddo

        close(11)
        close(12)

        print*,'Done'
        print*


        stop 
        end
