c     This function sets the bands by reading the bands names and the
c     zeropoints of the filters.  The arguments are as follow:
c
c     filtname = string containing the name of the file on which the
c                bands and zero points are listed (character*100)
c
c     Note that in the filters file, the bands have to be listed in the
c     same order as the magnitudes are supplied. This is the first
c     function to be called.
c
      subroutine setfilt(filtname)
      implicit real*8(a-h,o-z)
      parameter (NCMAX=32,NWMAX=350,NSMAX=4,NTMAX=4)

      character filtname*(*)
      character*200 bname(NCMAX)
      common /bands/bname
      integer inorm2(NCMAX)
      common /normal/inorm2
      real*8 jyzero(NCMAX),con(NCMAX),lbar(NCMAX)
      common /cal1/jyzero,con,lbar
      character*100 line
      real*8 magzero(NCMAX)
      real*8 jy(NCMAX),ejy(NCMAX)
      integer nchan
      common /data1b/jy,ejy,nchan
      integer verbose
      common /verb/verbose
      
      if(verbose.eq.1) then
         print*
         print*,
     *    '*********************************************************'
         print*,'Reading Filters and Setting Normalizations...'
      endif
      
      
      jchan = 0
      open(unit=14,file=filtname,status='old')

 100  read(14,'(a)',end=101)line
        if(line(1:1).eq.'#'.or.line.eq.' ') goto 100
        jchan = jchan + 1
        read(line,*)bname(jchan),inorm2(jchan),jyzero(jchan)
        goto 100
 101  continue

      close(14)
      nchan = jchan
      call setconstants(nchan)

      if(verbose.eq.1) then      
         print*,
     *    '*********************************************************'
         print*
      endif
      
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     This function calculates all the necessary constants for the
c     filters calibrations.
c      
      subroutine setconstants(nchan)
      implicit real*8 (a-h,o-z)
      parameter (NCMAX=32,NGMAX=17000,NWMAX=350,NSMAX=4,NTMAX=4)

      character*200 bname(NCMAX)
      common /bands/bname

      real*8 jyzero(NCMAX),con(NCMAX),lbar(NCMAX)
      common /cal1/jyzero,con,lbar

      integer nfilt(NCMAX)
      real*8 lamf(NCMAX,10000),eff(NCMAX,10000)
      common /filt1/lamf,eff,nfilt

      integer inorm2(NCMAX)
      common /normal/inorm2

      integer nvega
      real*8 lvega(11000),vegaraw(11000),vinterp(11000)
      character*200 fname,path,auxname
      integer VEGA,IRAC,AB
      character*20 tname

      integer verbose
      common /verb/verbose

c     Set normalization meanings.
      VEGA = 1
      IRAC = 2
      AB   = 3
      if(verbose.eq.1) then
         print*,'Setting up ',nchan,' bands '
      endif
c
c     Get Path
      call getenv('LRTPATH',path)
      call noblank(path,ip1,ip2)

      if(ip1.eq.101) then
         write(path,*)"."
         ip1 = 2
         ip2 = 2
      endif
c
c     Read in the Vega spectrum (flambda), convert angstroms to microns
c     and flambda to fnu.
c
      open(unit=13,file=path(ip1:ip2)//'/specs/Vega.sed',status='old')
      read(13,*)nvega
      do i=1,nvega
         read(13,*)lvega(i),vegaraw(i)
         lvega(i)   = lvega(i)/10000.0
         vegaraw(i) = vegaraw(i)*lvega(i)*lvega(i)
      enddo
      close(unit=13)
c
c     Read the Filters
c
      do jchan=1,nchan
         call noblank(bname(jchan),ib1,ib2)
         inorm = inorm2(jchan)
         open(unit=13,file=path(ip1:ip2)//'/Filters/'//bname(jchan
     *        )(ib1:ib2)//'.filter',status='old')
         read(13,*)nfilt(jchan)
         do i=1,nfilt(jchan)
            read(13,*)lamf(jchan,i),eff(jchan,i)
            lamf(jchan,i) = lamf(jchan,i)/10000.0
            if ((lamf(jchan,i).gt.100).or.(lamf(jchan,i).lt.0.01)) then 
               write(0,*)'you probably have an angstron/micron conversion 
     *              problem '
               write(0,*)'in setconstants for filter ',fname
               stop
            endif
         enddo
         close(unit=13)
c
c     Interpolate the vega spectrum onto this grid 
c
         if (inorm.eq.VEGA) then
            do kwave=1,nfilt(jchan)
               do ll=2,nvega
                  if ((lvega(ll-1).le.lamf(jchan,kwave)).and.
     *                 (lvega(ll).gt.lamf(jchan,kwave))) then
                     slope          = (lamf(jchan,kwave)-lvega(ll-1))/
     *                    (lvega(ll)-lvega(ll-1))
                     vinterp(kwave) = vegaraw(ll-1) + 
     *                    slope*(vegaraw(ll)-vegaraw(ll-1))
                  endif
               enddo
            enddo
         endif
c
c     Construct the normalization constants
c     For IRAC spectrum is  1/nu = lamf
c     For optical spectrum is Vega
c     For AB is the AB flat fnu source.
c
         a1 = 0.0
         a2 = 0.0
         a3 = 0.0
         do kwave=2,nfilt(jchan)
            iok = 0
            if (inorm.eq.IRAC) then
               slow = lamf(jchan,kwave-1)
               shig = lamf(jchan,kwave)
               iok  = 1
            endif
            if (inorm.eq.VEGA) then
               slow = vinterp(kwave-1)
               shig = vinterp(kwave)
               iok  = 1
            endif
            if (inorm.eq.AB) then
               slow = 1.0
               shig = 1.0
               iok  = 1
            endif
            if (iok.eq.0) then
               write(0,*),'unknown normalizing spectrum ',inorm,
     *              ' for band ',bname(jchan)
               stop
            endif

            del  = lamf(jchan,kwave)-lamf(jchan,kwave-1)

            vlow = eff(jchan,kwave-1)/lamf(jchan,kwave-1)**2
            vhig = eff(jchan,kwave  )/lamf(jchan,kwave)**2
            a1 = a1 + 0.5*del*(vlow+vhig)

            vlow = eff(jchan,kwave-1)*slow/lamf(jchan,kwave-1)**2
            vhig = eff(jchan,kwave  )*shig/lamf(jchan,kwave  )**2
            a2 = a2 + 0.5*del*(vlow+vhig)

            vlow = eff(jchan,kwave-1)*slow/lamf(jchan,kwave-1)
            vhig = eff(jchan,kwave  )*shig/lamf(jchan,kwave  )
            a3 = a3 + 0.5*del*(vlow+vhig)

         enddo
         con(jchan) = a2/a1/a3   
         lbar(jchan) = a3/a2    
      enddo
      if(verbose.eq.1) then
         print*,'Done'
      endif

      return 
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     This routine sets the templates of Assef, Kochanek et al. 2007 (in
c     prep.).  Arguments are as follows:
c
c     num = number of templates to define the model. Can only be 3 or
c           4. (integer)
c
c     Note that this is the second function to be called, right after
c     setfilt.
c
      subroutine settemp(num)
      implicit real*8 (a-h,o-z)
      parameter (NCMAX=32,NWMAX=350,NSMAX=4,NTMAX=4)
      
      integer num

      real*8 bedge(NWMAX)
      real*8 bcen(NWMAX)
      common /wavegrid/bedge,bcen,nwave

      real*8 spec(NSMAX,NWMAX),specuse(NSMAX,NWMAX)
      common /specmod1/spec,specuse,nspec

      real*8 alpha_norm(NSMAX)
      common /alphanorm/alpha_norm

      real*8 wgt(NCMAX,NWMAX)
      real*8 c(NCMAX)
      common /weights1/wgt,c
      real*8 jy(NCMAX),ejy(NCMAX)
      integer nchan
      common /data1b/jy,ejy,nchan

      character*100 specname
      character*100 channame
      character*200 path

      integer verbose
      common /verb/verbose

      integer ivaryobj(NSMAX)
      common /ivary/ivaryobj

      integer use_red_igm_prior
      common /red_igm_prior/use_red_igm_prior

c     Get Path
      call getenv('LRTPATH',path)
      call noblank(path,ip1,ip2)

c     Set the Path to the current directory if it has not been set
      if(ip1.eq.101) then
         write(path,*)"."
         ip1 = 2
         ip2 = 2
      endif

c     Set the name of the file with the spectra.
      if(num.eq.1) then
         write(specname,*)path(ip1:ip2)
     *        ,'/specs/seds.dat'
      else if(num.eq.2) then
         write(specname,*)path(ip1:ip2)
     *        ,'/specs/richards_seds.dat'
      else
         write(0,*)'Invalid value for inum: ',num
         write(0,*)'Exiting program.'
         stop
      endif

c     Set the range of wavelengths uniformly spaced in log space. The
c     templates that will be read will be compared to this grid. If you
c     which to change the templates remember to modify this part of the
c     code or comply to the grid.
      nwave     =  300
      bminl     =  dlog10(0.03d0)
      bmaxl     =  dlog10(30.0d0)
      dbl       = (bmaxl-bminl)/float(nwave)
      do kwave=1,nwave+1
         bedge(kwave)  =  10.d0**(bminl+dbl*float(kwave-1))
      enddo
      do kwave=1,nwave
         bcen(kwave) = 0.5d0*(bedge(kwave)+bedge(kwave+1))
      enddo

c     Read the templates and check the wavelenghts matches with the
c     grid.
      if(verbose.eq.1) then
         print*
         print*,
     *   '*********************************************************'
         print*,'Reading Templates... '
      endif

      call noblank(specname,is1,is2)
      open(unit=13,file=specname(is1:is2),form='formatted',status='old')
      read(13,*)nwt,nspec
      if(verbose.eq.1) then
         print*,'Setting up ',nspec,' Templates'
      endif
      if (nwt.ne.nwave) then
         write(0,*)'startup file has too few wavelengths ',nwt,nwave
         stop
      endif
      do kwave=1,nwave
         read(13,*)bt,tt,(spec(ll,kwave),ll=1,nspec)
         if (dabs(bt-bcen(kwave))/bt.gt.0.01d0) then
            write(0,*)'wavelength mismatch in startup file ',bt,
     *           bcen(kwave)
            stop
         endif
      enddo
      close(unit=13)  

c     Normalize the templates to 10**10 L_sun at 10 pc. This is only
c     important for determining bolometric luminosities. Note that
c     during most of the running of the codes, this normalization is
c     neglected and is only used when we want to get physical units for
c     the contribution of each component to the SED.  
c     Note that for the AGN template we only integrate longwards of
c     Lyalpha.
      do l = 1,nspec
         alpha_norm(l) = 0.d0
         do kwave=1,nwave-1
            if(l.ne.1.or.bedge(kwave).gt.0.1216d0) then
               alpha_norm(l) = alpha_norm(l) + 0.5d0*(spec(l,kwave+1)
     *              /bedge(kwave+1)**2+ spec(l,kwave)/bedge(kwave)**2)
     *              *(bedge(kwave+1)-bedge(kwave))
            endif
         enddo
         alpha_norm(l) = 3.826d5/(4*3.14159d0*3.086d0**2*alpha_norm(l))
      enddo

c     Initialize the ivaryobj files. The default is to use all templates.
      do l=1,nspec
         ivaryobj(l) = 1
      enddo

c     Initialize the use of the reddening prior. 
      use_red_igm_prior = 1

      if(verbose.eq.1) then
         print*,'Done'
         print*,
     c   '*********************************************************'
         print*
      endif

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     This function gets the weights needed to fit the galaxies. It
c     calculates the integral of the filter over two continuous steps in
c     lambda. See the paper for a clearer explanation.
c
      function getweight(zuse,jchan,kwave)
      implicit real*8 (a-h,o-z)
      parameter (NCMAX=32,NWMAX=350,NSMAX=4,NTMAX=4)
      
      real*8 bedge(NWMAX)
      real*8 bcen(NWMAX)
      common /wavegrid/bedge,bcen,nwave
      real*8 fterp(20),wterp(20)
      
      real*8 jyzero(NCMAX),con(NCMAX),lbar(NCMAX)
      common /cal1/jyzero,con,lbar
      
      integer nfilt(NCMAX)
      real*8 lamf(NCMAX,10000),eff(NCMAX,10000)
      common /filt1/lamf,eff,nfilt
      
      wemin = bedge(kwave  )*(1.0+zuse)
      wemax = bedge(kwave+1)*(1.0+zuse)

      nterp = 20
      jj    = 1
      dwe   = (wemax-wemin)/float(nterp-1)
      do kk=1,nterp 
         wterp(kk) = wemin + dwe*float(kk-1)
         fterp(kk) = 0.d0
         wedge     = wterp(kk)
         if ((wedge.ge.lamf(jchan,1)).and.
     *        (wedge.lt.lamf(jchan,nfilt(jchan)))) then
 100        if (wedge.ge.lamf(jchan,jj+1)) then
               jj = jj + 1
               go to 100
            endif
            if (jj+1.gt.nfilt(jchan)) then
               write(0,*)'bracket went past edge of filter '
               stop
            endif
            if ((wedge.lt.lamf(jchan,jj)).or.
     *           (wedge.gt.lamf(jchan,jj+1))) then
               write(0,*)'error in bracket ',
     *              lamf(jchan,jj),wedge,lamf(jchan,jj+1)
               write(0,*)jj,nfilt(jchan)
               stop
            endif
            slope     = (wedge-lamf(jchan,jj))/(lamf(jchan,jj+1)
     *           -lamf(jchan,jj))
            fterp(kk) = eff(jchan,jj) + slope*(eff(jchan,jj+1)-eff(jchan
     *           ,jj))
         endif
      enddo

c Now do the integral
      getweight = 0.d0
      do kk=1,nterp-1 
         wlow      = wterp(kk)
         whig      = wterp(kk+1)
         vlow      = fterp(kk  )/wlow
         vhig      = fterp(kk+1)/whig
         getweight = getweight + 0.5d0*(whig-wlow)*(vhig+vlow)
      enddo

      getweight = con(jchan)*getweight

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     This subroutine calculates the wavelength ranges on which the
c     filters are really contributing to the overall flux. This is set
c     for optimization of the algorithm but it is not terribly
c     important.
c
      subroutine getrange(jchan)
      implicit real*8 (a-h,o-z)
      parameter (NCMAX=32,NGMAX=17000,NWMAX=350,NSMAX=4,NTMAX=4)
      
      real*8 wgt(NCMAX,NWMAX)
      real*8 c(NCMAX)
      common /weights1/wgt,c
      integer jwmin(NCMAX),jwmax(NCMAX)
      common /weights2/jwmin,jwmax
      
      real*8 bedge(NWMAX)
      real*8 bcen(NWMAX)
      common /wavegrid/bedge,bcen,nwave
      
      real*8 wint(NWMAX)

      integer jyuse(NCMAX)
      common /data2/jyuse

      integer verbose
      common /verb/verbose


c     Work out the integral
      wint(1) = wgt(jchan,1)
      do kwave=2,nwave
         wint(kwave) = wint(kwave-1) + wgt(jchan,kwave)
      enddo


c     Now work out the contributing range
      wnorm=wint(nwave)


c     Trap worst case error. Do not use the channel for fitting if this
c     is the case.
      if (wnorm.le.0) then
         jwmin(jchan) = 1
         jwmax(jchan) = 1
         jyuse(jchan) = 0
c         write(0,*)'no contribution to weights ',jchan,wnorm
         return
      endif

c     Get the wavelength ranges.
      do kwave=1,nwave
         wint(kwave) = wint(kwave)/wnorm
      enddo
      do kwave=1,nwave
         if (wint(kwave).lt.0.999) jwmax(jchan) = kwave 
      enddo
      do kwave=nwave,1,-1
         if (wint(kwave).gt.0.001) jwmin(jchan) = kwave 
      enddo

c     Now needs to extend the coverage 1 resolution element to really
c     include 99.9% of the flux in the band.  
c     If the filter has a limit at the end of the spectrum, eliminate
c     the channel from the fit.
      if(jwmin(jchan).gt.1) then
         jwmin(jchan) = jwmin(jchan)-1
      else
         jyuse(jchan) = 0
         if(verbose.eq.1) then
c            write(0,*)'Channel ',jchan
c     *           ,' is partially outside of the SEDs wavelength range.',
c     *           'Channel is eliminated from the fit.'
         endif
      endif
      if(jwmax(jchan).lt.nwave) then
         jwmax(jchan) = jwmax(jchan)+1
      else
         jyuse(jchan) = 0
         if(verbose.eq.1) then
            write(0,*)'Channel ',jchan
     *           ,' is partially outside of the SEDs wavelength range.',
     *           'Channel is eliminated from the fit.'
         endif
      endif         

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c     This subroutine sets the luminosity distances for all redshifts in
c     the table.
c     The arguments are
c     
c     zin = vector containing the redshifts of the table (real*8)
c 
c     nin = number of redshifts in vector zin (integer)
c
c     This function does not need to be called by the user at any
c     time. 
c
      subroutine setdist(zin,nin)
      parameter (NCMAX=32,NWMAX=350,NSMAX=4)
      implicit real*8 (a-h,o-z)
      
      real*8 zin(*)
      
      real*8 dcom(5000),dlum(5000),dmod(5000),dvol(5000),rstar(5000)
      common /distance/dcom,dlum,dmod,dvol,rstar

      real*8 zin2(5000)      
      common /distance2/zin2,dz,zinmin,zinmax,nin2
      
      real*8 om,ol,ok,H0,DH
      common /cosmo/om,ol,ok,H0,DH

      real*8 DM(10000)
      
      real*8 aux(4)
      character*50 line

      integer verbose
      common /verb/verbose

      nin2 = nin
      do i=1,nin2
         zin2(i) = zin(i)
      enddo

      zinmin = zin2(1)
      zinmax = zin2(nin)
      dz = (zinmax-zinmin)/float(nin-1)

c     Set up the cosmology
      if(verbose.eq.1) then
         print*
         print*,
     c   '*********************************************************'
         print*,'Setting Cosmology'
      endif


c     Look for the cosmology file cosmo.dat . If not found, set the
c     defaults.
      jj=1
      icheck=0
      open(unit=20,file='cosmo.dat',status='old',err=101)
 100  read(20,*,end=102)line 
        if(icheck.eq.0.and.verbose.eq.1) then
           print*,'Found file cosmo.dat'
           print*,'Setting up User Cosmology'
           print*
           icheck=1
        endif           
        if(line(1:1).eq.'#'.or.line.eq.' ') goto 100
        read(line,*)aux(jj)
        jj = jj+1
        goto 100

 101  continue
      if(verbose.eq.1) then
         print*, 'No Cosmology File'
         print*, 'Setting Defaults'
         print*
      endif
      aux(1) =  0.3d0         !Omega_Matter
      aux(2) =  0.7d0         !Omega_Lambda
      aux(3) =  0.0d0         !Omega_k   
      aux(4) = 70.0d0         !Hubble Parameter Today in km/s/Mpc

 102  continue
      close(20)

      om = aux(1)
      ol = aux(2)
      ok = aux(3)
      H0 = aux(4)

      if(verbose.eq.1) then
         print*,'Omega_Matter = ',om
         print*,'Omega_Lambda = ',ol
         print*,'Omega_K      = ',ok
         print*,'H0           = ',H0
         print*
      endif


c     Set the Hubble Distance
      DH = 3d5/H0

c     Calculate the Comoving Distance. Calculations follow the
c     formalism of Hogg, D., 1999, astro-ph/9905116
      npt     = 100
      dcom(1) = 0.d0
      do i=1,nin-1
         zmin = zin(i)
         zmax = zin(i+1)
         delz = (zmax-zmin)/float(npt-1) 
         dadd = 0.d0
         do k=1,npt
            ztemp = zmin + delz*float(k-1)
            x   = 1+ztemp
            val = 1.d0/sqrt(om*x**3+ok*x**2+ol)
            if (k.ne.1) dadd = dadd + 0.5d0*delz*(val+vold) 
            vold  = val
         enddo
         dcom(i+1) = dcom(i) + DH*dadd 
      enddo

c     Determine the transverse comoving distance
      if(ok.eq.0.d0) then 
         do i=1,nin
            DM(i) = dcom(i)
         enddo
      else if(ok.gt.0.d0) then
         do i=1,nin
            DM(i) = DH/dsqrt(ok) * dsinh(dsqrt(ok)*dcom(i)/DH)
         enddo
      else
         do i=1,nin
            DM(i) = DH/dsqrt(-ok) * dsin(dsqrt(-ok)*dcom(i)/DH)
         enddo
      endif

c     Calculate now the Luminosity Distance, the Distance Modulus
c     and the Comoving Volume.
      do i=1,nin
         dlum(i) = DM(i)*(1.d0+zin(i))
         x       = 1.d0+zin(i)
         dvol(i) = 2.d0*DH*DM(i)*DM(i)/dsqrt(om*x**3+ok*x**2+ol)
         if (i.ne.1) then
            dmod(i) = 5.d0*dlog10(1.0d6*dlum(i)/10.d0)
         endif
      enddo
      dmod(1) = dmod(2)

      if(verbose.eq.1) then
         print*,'Done'
         print*,
     c   '*********************************************************'
         print*
      endif

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine closemag(jchan,nchan,jclose)
      implicit real*8 (a-h,o-z)
      parameter (NCMAX=32)

      real*8 jyzero(NCMAX),con(NCMAX),lbar(NCMAX)
      common /cal1/jyzero,con,lbar

      integer jyuse(NCMAX)
      common /data2/jyuse

      real*8 dist,dist_min

      dist_min = 1.d10

      do j=1,nchan
         if(jyuse(j).eq.1) then
            dist = abs(lbar(jchan)-lbar(j))
            if(dist.lt.dist_min) then
               dist_min = dist
               jclose = j
            endif
         endif
      enddo
      
      return
      end    

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine total_corr(z,jchan,nchan,jcorr,magtotuse)
      implicit real*8 (a-h,o-z)
      parameter (NCMAX=32)

      real*8 jyzero(NCMAX),con(NCMAX),lbar(NCMAX)
      common /cal1/jyzero,con,lbar

      integer jyuse(NCMAX)
      common /data2/jyuse

      integer magtotuse(NCMAX)



      dist_min = 1d20

      do j=1,nchan
         if(jyuse(j).eq.1.and.magtotuse(j).eq.1) then
            dist = abs(lbar(jchan)-lbar(j)/(1+z))
            if(dist.lt.dist_min) then
               dist_min = dist
               jcorr = j
            endif
         endif
      enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     the transmission bluewards of Lyman alpha
c
      function transmit(alamrest,zqso,scale)
      implicit real*8 (a-h,o-z)
      
c     alamrest = rest wavelength
      alamobs = alamrest*(1.0+zqso)
      
c     wavelength of lyman alpha in microns
      wavelya = 0.121567d0
      wavelyb = 0.102518d0
      if (alamrest.ge.wavelya) then
         transmit = 1.d0
         return
      endif
      
c     this model is from Fan et al. AJ 132 117 2006
c     Lyman alpha absorbers
      onepzabs = alamobs/wavelya
      tau      = 0.85d0*(0.2d0*onepzabs)**4.3d0
c     Lyman beta absorbers
      if (alamrest.le.wavelyb) then
         onepzabs = alamobs/wavelyb
         tau      = tau + 0.38d0*(0.2d0*onepzabs)**4.3d0
      endif
      
c     Lyman limit absorbers -- uses number of lyman limit systems from
c     Stengler-Larrea et al. ApJ 444 64 1995 assumes optical depth of
c     tau0 for an absorber counts up mean number of absorbers to the
c     quasar Lyman limit multiplies by tau0 and uses that as the optical
c     depth to shorter wavelengths wavelength of lyman limit in microns
      wavelylim = 0.091127d0
      if (alamrest.le.wavelylim) then
         tau0     = 2.d0
         onepzabs = 1.d0 + zqso
         tau      = tau + 0.1d0*tau0*(onepzabs**2.5d0-1.d0)
      endif
      
      transmit = exp(-scale*tau)
      
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Extinction curves. Details in Assef et al. (2008)
c
      function rl(x,rv)
      implicit real*8 (a-h,o-z)

      if (x.lt.1.1) then
         a =  0.574*x**1.61
         b = -0.527*x**1.61
      else
         if (x.lt.3.3) then
            y = x-1.82
            a = 1.0+0.17699*y-0.50447*y*y-0.02427*y**3+0.72085*y**4 +
     1           0.01979*y**5-0.77530*y**6+0.32999*y**7
            b = 1.41338*y+2.28305*y**2+1.07233*y**3-5.38434*y**4-
     1           0.62251*y**5+5.30260*y**6-2.09002*y**7
         else
            if (x.lt.8.0) then
               if (x.lt.5.9) then
                  fa = 0.0
                  fb = 0.0
               else
                  fa = -0.04473*(x-5.9)**2 - 0.009779*(x-5.9)**3
                  fb =  0.21300*(x-5.9)**2 + 0.120700*(x-5.9)**3
               endif
               a =  1.752 - 0.316*x - 0.104/((x-4.67)**2+0.341) + fa
               b = -3.090 + 1.825*x + 1.206/((x-4.67)**2+0.263) + fb
            else
               y = x-8.0
               a = -1.073-0.628*y+0.137*y**2-0.070*y**3
               b = 13.670+4.257*y-0.420*y**2+0.374*y**3
            endif
         endif
      endif
      
      rgal = a+b/rv
      
c     compute the SMC curve
      c1 = -5.68
      c2 =  2.53
      c3 =  0.76
      x0 =  4.56
      g0 =  1.68
      c4 =  0.60
      
      df = x*x/((x*x-x0*x0)**2+(x*g0)**2)
      if (x.gt.5.9) then
         ff = 0.5392*(x-5.9)**2+0.05644*(x-5.9)**3
      else
         ff = 0.0
      endif
      asmc  = 1.0
      bsmc  = c1 + c2*x + c3*df + c4*ff
      
      rsmc = asmc+bsmc/rv
      
      if (x.lt.4.0) then
         rl = max(rgal,rsmc)
      else
         rl = rsmc
      endif
      
      rl = rv*rl
      
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Subroutine to initialize the reddening absorption curves.
c
      subroutine set_red(ired)
      implicit real*8 (a-h,o-z)
      parameter(NWMAX=350)

      real*8 tau(NWMAX),ebv,igm
      common /dust/tau,ebv,igm

      real*8 bedge(NWMAX)
      real*8 bcen(NWMAX)
      common /wavegrid/bedge,bcen,nwave

      real*8 emin,emax,de
      real*8 gmin,gmax,dg
      integer ne,ng
      common /redpars/emin,emax,de,ne
      common /igmpars/gmin,gmax,dg,ng

      if(ired.ne.0.and.ired.ne.1) then
         write(0,*)'Invalid value for parameter ired: ',ired
         write(0,*)'Exiting program'
      endif


      if(ired.eq.1) then
         rv = 3.1d0
         do kwave=1,nwave
            xarg       = 1.d0/bcen(kwave)
            tau(kwave) = rl(xarg,rv)
         enddo
         emin = -2.0d0
         emax =  1.5d0
         ne   =  36
         de   = (emax-emin)/(1.d0*float(max(ne-1,1)))
         ne   = ne + 1 !needs to add 1 so that the first value is 0.
      else
         emin = 0.d0
         emax = 0.d0
         ne   = 1
         de   = 0.d0
      endif

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Subroutine to initialize the IGM absorption.
c
      subroutine set_igm(iigm)
      implicit real*8 (a-h,o-z)
      parameter(NWMAX=350)

      real*8 tau(NWMAX),ebv,igm
      common /dust/tau,ebv,igm

      real*8 bedge(NWMAX)
      real*8 bcen(NWMAX)
      common /wavegrid/bedge,bcen,nwave

      real*8 emin,emax,de
      real*8 gmin,gmax,dg
      integer ne,ng
      common /redpars/emin,emax,de,ne
      common /igmpars/gmin,gmax,dg,ng


      if(iigm.ne.0.and.iigm.ne.1) then
         write(0,*)'Invalid value for parameter iigm: ',iigm
         write(0,*)'Exiting program'
      endif

      if(iigm.eq.1) then
         gmin = 0.0d0
         gmax = 1.4d0
         ng   = 8
         dg   = (gmax-gmin)/float(max(ng-1,1))         
      else
         gmin = 1.d0
         gmax = 1.d0
         ng   = 1
         dg   = 0.d0
      endif

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Subroutine to read the zero point corrections. It must be called
c     after setfilt and before settemp.
c
      subroutine read_zpc
      implicit real*8 (a-h,o-z)
      parameter (NCMAX=32,NWMAX=350,NSMAX=4,NTMAX=4,NGMAX=17000)

      real*8 wgt(NCMAX,NWMAX)
      real*8 c(NCMAX)
      common /weights1/wgt,c

      real*8 jy(NCMAX),ejy(NCMAX)
      integer nchan
      common /data1b/jy,ejy,nchan

      character*500 line

      integer verbose
      common /verb/verbose


c     See if the file with the zero point corrections is present.
      jj = 1
      icheck = 0
      open(unit=16,file='channel.zpc',status='old',err=101)
 100  read(16,'(a)',end=102)line
         if(icheck.eq.0.and.verbose.eq.1) then
            print*
            print*,
     *    '*********************************************************'
            print*,'Found file channel.zpc'
            print*,'Reading zero point corrections'
            icheck=1
         endif
         if(line(1:1).eq.'#'.or.line.eq.' ') goto 100
         read(line,*)c(jj)
         jj = jj + 1
         goto 100
 101  continue
c     Only gets here if there is no file with zero point corrections.
         do j = 1,nchan
            c(j) = 1.d0
         enddo
 102  continue

c     Read the channel zero point corrections.
      if(verbose.eq.1.and.icheck.eq.1) then
         print*,'Done'
         print*,
     *    '*********************************************************'
         print*
      endif

      return
      end

