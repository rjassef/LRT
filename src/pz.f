c     Subroutine to initialize the K-Corrections program. This
c     subroutine will initialize the bands, the templates and the
c     luminosity priors.
c
      subroutine pzinit(filtname,inum,ired,iigm,zmin,zmax,dz,verb_flag)
      implicit real*8 (a-h,o-z)
      parameter (NCMAX=40,NWMAX=350,NSMAX=4,NTMAX=4)

      character filtname*(*)

      integer verbose,verb_flag
      common /verb/verbose

      integer pzon
      common /regen/zmax2,zmin2,dz2,pzon

      verbose = verb_flag
      pzon = 1
      zmax2 = zmax
      zmin2 = zmin
      dz2 = dz

c     Initialize the filters.
      call setfilt(filtname)
c     Read the zero point normalizations or otherwise initialize them
      call read_zpc
c     Initialize the templates.
      call settemp(inum)
c     Initialize the reddening.
      call set_red(ired)
c     Initialize the IGM absorption
      call set_igm(iigm)
c     Initialize Luminosity priors
      call setlumprior
c     Create the photoz grid.
      call photoz_grid(zmin,zmax,dz)
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Subroutine to set the luminosity priors. The function will look
c     for a file named 'prior.dat' . If it finds it, it will read
c     expecting the following order for the parameters (lines starting
c     with # will be skipped):
c
c     Use_prior
c     Channel Number
c     Mstar
c     Alpha
c
c     If Use_prior = 0, the priors will not be applied and if its 1,
c     they will be used. Any other number will make the program crash.
c     Channel number refers to the number of the band for which the
c     luminosity prior will be applied. The number is the position on
c     which the band was declared on the magnitudes file.  Mstar and
c     alpha are the Schechter Function (Schechter, 1976) parameters.
c
      subroutine setlumprior
      implicit real*8 (a-h,o-z)
      parameter (NCMAX=40,NWMAX=350,NSMAX=4,NTMAX=4)

      real*8 mstar,alpha
      integer uselump,lumchan
      common /lumprior/mstar,alpha,uselump,lumchan

      real*8 prior(5)
      character*100 line

      integer verbose
      common /verb/verbose

c     Set up the luminosity priors.
      if(verbose.eq.1) then
         print*
         print*,
     *    '*********************************************************'
         print*,'Setting the Luminosity priors'
      endif


c     Try to open the file prior.dat and read the parameters. If the
c     file is not found, go to the default prior.
      jj = 1
      icheck=0
      open(unit=16,file='prior.dat',status='old',err=201)
 200  read(16,'(a)',end=202)line
        if(icheck.eq.0.and.verbose.eq.1) then
           print*,'Found file prior.dat'
           print*,'Setting up User Luminosity Prior'
           print*
           icheck=1
        endif
        if(line(1:1).eq.'#'.or.line.eq.' ') goto 200
        read(line,*)prior(jj)
        jj = jj+1
        goto 200

c     Jump here if the prior.dat file doesn't exist.
 201  continue
      if(verbose.eq.1) then
         print*,'No priors file found'
         print*,'Assuming Default Behaviour'
         print*
      endif
      prior(1) = 0d0
      prior(2) = 2d0
      prior(3) = -21.4d0
      prior(4) = -0.7d0

 202  continue
      close(16)

      uselump = int(prior(1))
      lumchan = int(prior(2))
      mstar   = prior(3)
      alpha   = prior(4)

c     Print settings.
      if(uselump.eq.1) then
         if(verbose.eq.1) then
            print*,'Luminosity Priors = On'
            print*,'Band   = ',lumchan
            print*,'M_star = ',mstar
            print*,'alpha  = ',alpha
         endif
      else if(uselump.eq.0) then
         if(verbose.eq.1) then
            print*,'Luminosity Priors = Off'
         endif
      else
         write(0,*)'Invalid number for using luminosity priors.'
         write(0,*)'Recieved ',uselump,' when expecting 0 or 1.'
         write(0,*)'Aborting Program.'
         stop
      endif

      if(verbose.eq.1) then
         print*
         print*,'Done'
         print*,
     *    '*********************************************************'
         print*
      endif

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     This function writes the photoz tables.
c     It should only be called after setfilt and settemp.
c
      subroutine photoz_grid(zmin,zmax,dz)
      implicit real*8 (a-h,o-z)
      parameter (NCMAX=40,NGMAX=17000,NSMAX=4,NWMAX=350,NXMAX=5000)
      parameter (NEMAX=37,NIMAX=8)


      real*8 zbin(NGMAX)

      real*8 wgt(NCMAX,NWMAX)
      real*8 c(NCMAX)
      common /weights1/wgt,c
      integer jwmin(NCMAX),jwmax(NCMAX)
      common /weights2/jwmin,jwmax

      real*8 vec(NSMAX)
      real*8 jymod(NSMAX,NCMAX)
      real*8 jymodtot(NCMAX)
      common /models/jymod,jymodtot,vec

      real*8 z
      real*8 jy(NCMAX),ejy(NCMAX)
      integer nchan
      common /data1b/jy,ejy,nchan
      common /data1/z

      real*8 jyzero(NCMAX),con(NCMAX),lbar(NCMAX)
      common /cal1/jyzero,con,lbar

      real*8 spec(NSMAX,NWMAX),specuse(NSMAX,NWMAX)
      common /specmod1/spec,specuse,nspec

      real*8 bedge(NWMAX)
      real*8 bcen(NWMAX)
      common /wavegrid/bedge,bcen,nwave

      integer verbose
      common /verb/verbose

      real*8 tau(NWMAX),ebv,igm
      common /dust/tau,ebv,igm

      real*8 priorjy(1000,NSMAX,NCMAX)
      common /priortemp/priorjy

      real*8 emin,emax,de
      real*8 gmin,gmax,dg
      integer ne,ng
      common /redpars/emin,emax,de,ne
      common /igmpars/gmin,gmax,dg,ng

      real*8 tigm(NWMAX)

      real*8 zgal(NXMAX),galjy(NXMAX,NSMAX,NCMAX,NEMAX,NIMAX)
      common /galtemp/zgal,galjy,ngalz,ngalt

      character*200 path, line


c     Open the table file and write the header.
      nz   = (zmax-zmin)/dz + 1

      ngalz = nz
      ngalt = nspec

c     Start the iteration.
      if(verbose.eq.1) then
         print*
         print*,
     *   '*********************************************************'
         print*,'Constructing Photo-z Grid'
      endif
      do iz=1,nz
         zbin(iz) = zmin + dz*float(iz-1)
         zgal(iz) = zbin(iz)
         zuse = zbin(iz)
c     Build the weights.
         do jchan=1,nchan
            do kwave=1,nwave
               wgt(jchan,kwave) = getweight(zuse,jchan,kwave)
            enddo
            call getrange(jchan)
         enddo
         do ie=1,ne
            if(ie.eq.1) then
               euse = 0.d0
            else
               euse = emin + de*float(ie-2)
               euse = 10.d0**euse
            endif
            do ig=1,ng
               guse = gmin + dg*float(ig-1)

               do kwave=1,nwave
                  tigm(kwave) = transmit(bcen(kwave),zuse,guse)
               enddo

               do jchan=1,nchan
c     Build the flux models only if it changes anything.
                  if(jwmin(jchan).le.0.121567d0.and.
     *                 ig.gt.1) goto 200

                  do l=1,nspec
                     galjy(iz,l,jchan,ie,ig) = 0.d0
                     priorjy(iz,l,jchan) = 0.d0
                     do kwave=jwmin(jchan),jwmax(jchan)
                        if(l.ne.1) then
                           dust = 1.d0
                        else
                           dust = 10.d0**(-0.4d0*tau(kwave)*euse)
                        endif
                        galjy(iz,l,jchan,ie,ig) =
     *                       galjy(iz,l,jchan,ie,ig) + c(jchan)*
     *                       (1.d0+zbin(iz))*spec(l,kwave)*
     *                       wgt(jchan,kwave)*dust*tigm(kwave)
                        priorjy(iz,l,jchan) = priorjy(iz,l,jchan) +
     *                       c(jchan)*(1.d0+zbin(iz))*spec(l,kwave)*
     *                       wgt(jchan,kwave)
                     enddo
                  enddo
 200              continue
               enddo

            enddo
         enddo
      enddo

      if(verbose.eq.1) then
         print*,'Done'
         print*,
     *   '*********************************************************'
         print*
      endif

      call setdist(zgal,nz)

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     This function calculates photometric redshifts for galaxies based
c     on the algorithm and templates described by Assef & Kochanek et
c     al., 2007 (in preparation).  The function arguments are:
c
c
c     mag    =  NCMAX dimensionl vector with the observed magnitudes
c               of the galaxies (real*8)
c
c     emag   =  NCMAX dimensional vector with the errors in the
c               observed magnitudes (real*8)
c
c     maguse =  NCMAX dimensional vector which holds which mags will
c               be used for fitting the templates. mag(j) = 1 means that
c               band j will be used and 0 that it will not (integer)
c
c     num    =  number of bands (integer)
c
c     zobj   =  variable to return the best fitted redshift (real*8)
c
c     chigal =  variable to return the chi^2 of the best fit (real*8)
c
c     chinop =  variable to return the chi^2 of the best fit without the
c               prior contribution (real*8)
c
c     op     =  1 if mag in magnitudes. 0 if mag in Jy. (integer)
c
c     chi2zop=  1 to write the chi-squared distribution to the fort.90
c               file. 0 to not do it.
c
c
c     This function should only be called after the photoz tables have
c     been set with settable. Notice that this and all the subroutines
c     to follow do not necessarily follow the 80 character per line
c     fortran 77 convention.
c
      subroutine pza(mag,emag,maguse,zobj,chigal,chinop,op,chi2zop)
      implicit real*8(a-h,o-z)
      parameter (NCMAX=40,NWMAX=350,NSMAX=4,NTMAX=4)

      integer op, chi2zop

      real*8 mag(*),emag(*)
      integer maguse(*)

      real*8 jy(NCMAX),ejy(NCMAX)
      integer nchan
      common /data1b/jy,ejy,nchan

      integer jyuse(NCMAX)
      common /data2/jyuse

      real*8 jyzero(NCMAX),con(NCMAX),lbar(NCMAX)
      common /cal1/jyzero,con,lbar

      real*8 mstar,alpha
      integer uselump,lumchan
      common /lumprior/mstar,alpha,uselump,lumchan

      integer opchi2z
      common /chiz/opchi2z

      real*8 spec(NSMAX,NWMAX),specuse(NSMAX,NWMAX)
      common /specmod1/spec,specuse,nspec

Cf2py intent(inout) zobj,chigal,chinop

      opchi2z = chi2zop


c     See which bands will be used for fitting.
      m = 0
      do j = 1,nchan
         if(maguse(j).ne.0.and.maguse(j).ne.1.and.
     *        maguse(j).ne.2) then
            write(0,*)'maguse(',j,'=',maguse(j),') not equal to 2, 1
     *        or 0 in function pza.'
            write(0,*)'Aborting program'
            stop
         else
            jyuse(j) = maguse(j)
            if(jyuse(j).ge.1) then
               m = m + 1
            endif
         endif
      enddo

c     Only continue if there are at least as many magnitudes as chosen
c     number of templates.
      if(m.lt.nspec) then
         write(0,*)'Too few magnitudes to calculate photo-z (<'
     *        ,nspec,').'
         zobj = -1.d0
         chigal = -1.d0
         chinop = -1.d0
         goto 500
      endif

      if(op.eq.1) then
c     Transform from magnitudes to fluxes if op==1.
         do jchan = 1,nchan
            if(jyuse(jchan).eq.1) then
               jy(jchan)  = jyzero(jchan)*10.0**(-0.4*mag(jchan))
               ejy(jchan) = 0.4*log(10.0)*emag(jchan)*jy(jchan)
            else if(jyuse(jchan).eq.2) then
               jy(jchan)  = 0.d0
               ejy(jchan) = jyzero(jchan)*10.0**(-0.4*emag(jchan))
            else
               jy(jchan)  = 0.d0
               ejy(jchan) = 0.d0
            endif
            ejy(jchan) = ejy(jchan)**2
         enddo
      else if(op.eq.0) then
         do jchan = 1,nchan
            jy(jchan)  = mag(jchan)
            if(jyuse(jchan).eq.2) jy(jchan) = 0.d0
            ejy(jchan) = emag(jchan)**2
         enddo
      else
         write(0,*)'Wrong Value for op in function pza.'
         write(0,*)'Recieved ',op,' when expecting 0 or 1.'
         write(0,*)'Aborting Program.'
      endif

c     Fit for the redshift.
      call pza_fitgal(nmag,zobj,chigal,chinop)

 500  continue

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     This function fits the redshift to the galaxy.  The arguments are:
c
c
c     nmag   = variable to return the number of used bands (integer)
c
c     zbest  =  variable on which the best fit redshift will
c               be returned (real*8)
c
c     chigal =  variable on which the chi^2 of the best fit
c               is returned. When considering priors, this can
c               be negative (real*8)
c
c     chinop =  chi-squared of the best fit without the prior.
c
c
c     This fuction does not need to be called by the user at any
c     time. It is embbeded in the pza function.
c
      subroutine pza_fitgal(nmag,zbest,chigal,chinop)
      implicit real*8 (a-h,o-z)
      parameter (NCMAX=40,NWMAX=350,NSMAX=4,NXMAX=5000)
      parameter (NEMAX=37,NIMAX=8)

      real*8 zgal(NXMAX),galjy(NXMAX,NSMAX,NCMAX,NEMAX,NIMAX)
      common /galtemp/zgal,galjy,ngalz,ngalt

      real*8 jy(NCMAX),ejy(NCMAX)
      integer nchan
      common /data1b/jy,ejy,nchan

      integer jyuse(NCMAX)
      common /data2/jyuse

      real*8 jyzero(NCMAX),con(NCMAX),lbar(NCMAX)
      common /cal1/jyzero,con,lbar

      real*8 dcom(5000),dlum(5000),dmod(5000),dvol(5000),rstar(5000)
      common /distance/dcom,dlum,dmod,dvol,rstar

      real*8 a(NSMAX,NSMAX),b(NSMAX)
      real*8 asave(NSMAX,NSMAX),bsave(NSMAX)
      real*8 tempw(NSMAX),tempv(NSMAX,NSMAX),tempz(NSMAX)
      integer itempv(NSMAX)
      real*8 temps(NSMAX)

      real*8 vfit(NSMAX)

      real*8 mstar,alpha
      integer uselump,lumchan
      common /lumprior/mstar,alpha,uselump,lumchan

      real*8 vec(NSMAX)
      common /galtype/vec

      real*8 priorjy(1000,NSMAX,NCMAX)
      common /priortemp/priorjy

      real*8 tau(NWMAX),ebv,igm
      common /dust/tau,ebv,igm

      real*8 emin,emax,de
      real*8 gmin,gmax,dg
      integer ne,ng
      common /redpars/emin,emax,de,ne
      common /igmpars/gmin,gmax,dg,ng

      real*8 jymodx(NSMAX,NCMAX)
      common /modelsx/jymodx

      integer opchi2z
      common /chiz/opchi2z


c     Initialize values for return in case of failure.
      chigal    = 1.0d32
      chigalcut = 1.0d32
      zbest     = -1.d0

      if(opchi2z.ne.0) write(90,*)opchi2z,ngalz-1

c     Figure out the number of bands to use and exit if none can be used.
      nmag = 0
      do jchan=1,nchan
         if (jyuse(jchan).ge.1) nmag = nmag + 1
      enddo

c     Initialize the parameters for the iteration.
      chigal    = -1.d0
      zbest     = -1.d0
      istart    = 1
      nfitz = ngalz
      nfitt = ngalt

c     Start Main Cycle.
      do k=2,nfitz
         zval = zgal(k)

         do ie=1,ne
            if(ie.eq.1) then
               euse = 0.d0
            else
               euse = emin + de*float(ie-2)
               euse = 10.d0**euse
            endif
            do ig=1,ng
               guse = gmin + dg*float(ig-1)
               chival = (euse/0.5d0)**2 + ((guse-1.d0)/0.5d0)**2

c     Fill the matrices
               maxdim = NSMAX
               call clearmat(asave,bsave,maxdim,ngalt)
               do jchan=1,nchan
                  if (jyuse(jchan).ge.1) then
                     do l1=1,nfitt
                        vfit(l1) = galjy(k,l1,jchan,ie,ig)
                     enddo
                     do l1=1,nfitt
                        bsave(l1) = bsave(l1) + jy(jchan)*vfit(l1)/
     *                       ejy(jchan)
                        do l2=l1,nfitt
                           asave(l1,l2) =   asave(l1,l2) +
     *                          vfit(l1)*vfit(l2)/ejy(jchan)
                        enddo
                     enddo
                  endif
               enddo

c     Include the bright part of the galaxy luminosity prior if priors
c     are used.
               if(uselump.eq.1) then
                  jchan = lumchan
                  fstar = jyzero(jchan)*10.d0**(-0.4d0*(mstar+dmod(k)))
                  do l=2,nfitt
                     bsave(l) = bsave(l) -
     *                    priorjy(1,l,jchan)/fstar
                  enddo
               endif

               call symmat(asave,bsave,maxdim,nfitt)

               do l1=1,nfitt
                  b(l1) = bsave(l1)
                  do l2=1,nfitt
                     a(l1,l2) = asave(l1,l2)
                  enddo
               enddo

c     Solve assuming only positive coefficients. If convergence fails,
c     revert to the slower version going through all possible
c     combinations.
               call my_nnls_2(a,maxdim,nfitt,nfitt,b,temps,MODE,its,0)
               if(MODE.eq.3) then
                  do l1=1,nfitt
                     b(l1) = bsave(l1)
                     do l2=1,nfitt
                        a(l1,l2) = asave(l1,l2)
                     enddo
                  enddo
                  do l=1,nfitt
                     do j=1,nchan
                        jymodx(l,j) = galjy(k,l,j,ie,ig)
                     enddo
                  enddo
                  call ANNLS(a,maxdim,nfitt,nfitt,b,temps)
               endif

c               print*,temps,nfitt

c     Compute the goodness of fit
               do jchan=1,nchan
                  estjy = 0.d0
                  do l=1,nfitt
                     estjy = estjy + temps(l)*galjy(k,l,jchan,ie,ig)
                  enddo
                  if (jyuse(jchan).ge.1) then
                     chival = chival + (jy(jchan)-estjy)**2/ejy(jchan)
                  endif
               enddo

c     Apply luminosity priors if any.
               prior = 0.d0
               if(uselump.eq.1) then
                  estjy = 0.d0
                  jchan = lumchan
                  do l=2,nfitt
                     estjy = estjy + temps(l)*priorjy(1,l,jchan)
                  enddo
                  restr  = estjy
                  ineg = 0
                  if (estjy.le.0.d0) then
                     ineg = 1
                  else
                     rmag   = -2.5d0*dlog10(estjy/jyzero(jchan))
                     rabs   =  rmag-dmod(k)
                     xi     = -0.4d0*(rabs-mstar)
                     prior  = -2.0d0*(log(dvol(k))-10.d0**xi+(1.d0+alpha)*
     *                    xi*dlog(10.d0))
                  endif
                  chival    = chival+prior
               endif

c     Keep result if chi^2 is minimum or if its the first iteration.
               if ((chival.lt.chigal).or.(istart.eq.1)) then
                  istart     = 0
                  chigal     = chival
                  if(uselump.eq.1.and.ineg.eq.0) then
                     chinop = chival - prior
                  else
                     chinop  = chival
                  endif
                  zbest      = zval
                  do l=1,NSMAX
                     vec(l) = temps(l)
                  enddo
                  ebv = euse
                  igm = guse
                  chinop = chinop -
     *                 ((ebv/0.5d0)**2 + ((igm-1.d0)/0.5d0)**2)
               endif

            enddo
         enddo

c     If warranted, print the chi2 distribution.
         if(opchi2z.ne.0) then
            write(90,100)zval,chival-prior,chival
 100        format(3E20.6)
         endif

      enddo

      return
      end
