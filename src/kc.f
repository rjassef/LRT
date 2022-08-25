c     Subroutine to initialize the K-Corrections program. This
c     subroutine will initialize the bands, the templates and the
c     luminosity distance grid (not used by the main function).

c     inum = 1 -> Standard AGN SED
c     inum = 2 -> Richards et al. (2006) converged SED
c
      subroutine kcinit(filtname,inum,ired,iigm,verb_flag)
      implicit real*8 (a-h,o-z)
      parameter (NCMAX=40,NWMAX=350,NSMAX=4,NTMAX=4)

      character filtname*(*)
      real*8 zin(5000)

      integer verbose,verb_flag
      common /verb/verbose

      integer pzon
      common /regen/zmax2,zmin2,dz2,pzon

      verbose = verb_flag
      pzon = 0

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

c     First set redshift ranges to call setdist.
      zinmin = 0.0d0
      zinmax = 6.0d0
      nin  = 601
      dz = (zinmax-zinmin)/float(nin-1)
      do i = 1,nin
         zin(i) = zinmin + dz*float(i-1)
      enddo
c     Initialize luminosity distance grid
      call setdist(zin,nin)

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     The KCA function calculates K corrections based on the templates
c     generated by Assef & Kochanek et al. 2007 (in preparation).
c     The arguments are:
c
c
c     num     = number of bands on which the photometry is being
c               supplied (integer)
c
c     mag     = vector of n components containing the observed
c               magnitudes (real*8)
c
c     emag    = vector of n components containing the errors in the
c               magnitudes (real*8)
c
c     maguse  = vector of n components containing which magnitudes will
c               be used to fit the galaxy to the templates. 1 means to
c               be used and 0 to not be used.  For bands not used in
c               the fit, a kcorrection will also be calculated as well
c               as a model flux. (integer)
c
c     zobj    = redshift of the galaxy (real*8)
c
c     z0      = redshift to which you want to k correct (real*8)
c
c     magmod  = vector of n components on which the modeled magnitudes
c               will be returned (real*8)
c
c     magcorr = vector of n components on which the kcorrections in
c               magnitudes are returned (real*8)
c
c     comp    = vector with at least dimensions equal to the number
c               of templates used (real*8)
c
c     op      = 1 if mag is in magnitudes and 0 if it is in flux.
c
c
c     Note that the following functions do not necessarily follow the 80
c     character per line fortran 77 convention. Compilation with g77 has
c     to include the -ffixed-line-length-none flag. Compilation with
c     intel fortran compiler, ifort, has to include the -extend_source
c     flag. Also note that min.f should be compiled simultaneously.
c
c     This function should only be called after calling setfilt and
c     settemp.
c
      subroutine kca(mag,emag,maguse,zobj,z0,magmod,magcorr,comp,
     *     covx,xebv,xigm,chi2,op)
      implicit real*8(a-h,o-z)
      parameter (NCMAX=40,NWMAX=350,NSMAX=4,NTMAX=4)

      real*8 mag(*),emag(*),magmod(*),magcorr(*),comp(*)
      real*8 covx(NSMAX,NSMAX)
      real*8 xebv,xigm
      integer maguse(*),op

      real*8 z
      real*8 jy(NCMAX),ejy(NCMAX)
      integer nchan
      common /data1b/jy,ejy,nchan
      common /data1/z

      integer jyuse(NCMAX)
      common /data2/jyuse

      real*8 vec(NSMAX)
      real*8 jymod(NSMAX,NCMAX)
      real*8 jymodtot(NCMAX)
      common /models/jymod,jymodtot,vec

      real*8 jyzero(NCMAX),con(NCMAX),lbar(NCMAX)
      common /cal1/jyzero,con,lbar

      real*8 bedge(NWMAX)
      real*8 bcen(NWMAX)
      common /wavegrid/bedge,bcen,nwave

      real*8 wgt(NCMAX,NWMAX)
      real*8 c(NCMAX)
      common /weights1/wgt,c
      integer jwmin(NCMAX),jwmax(NCMAX)
      common /weights2/jwmin,jwmax

      integer ivaryobj(NSMAX)
      common /ivary/ivaryobj

      real*8 spec(NSMAX,NWMAX),specuse(NSMAX,NWMAX)
      common /specmod1/spec,specuse,nspec
      common /specnorm/bminnorm,bmaxnorm

      real*8 alpha_norm(NSMAX)
      common /alphanorm/alpha_norm

      real*8 wgt0(NCMAX,NWMAX)
      real*8 jymodz0(NSMAX,NCMAX),jymodz0tot(NCMAX),jycorr(NCMAX)

      real*8 tau(NWMAX),ebv,igm
      common /dust/tau,ebv,igm

      real*8 tigm(NWMAX)

      real*8 cov(NSMAX,NSMAX)
      common /vecerr/cov

      real*8 chimin
      common /minchi/chimin

Cf2py intent(inout) xebv,xigm,chi2

      z = zobj

c     See which bands will be used for fitting.
      m = 0
      do j = 1,nchan
         if(maguse(j).ne.0.and.maguse(j).ne.1.and.
     *        maguse(j).ne.2) then
            write(0,*)'maguse(',j,')=',maguse(j),' not equal to 1 or 0
     *           in function kca.'
            write(0,*)'Aborting program'
            stop
         else
            jyuse(j) = maguse(j)
            if(jyuse(j).ge.1) m = m + 1
         endif
      enddo

c     Only continue if there are at least as many magnitudes as chosen
c     number of templates.
      nspec_fit = 0
      do l=1,nspec
         nspec_fit = nspec_fit + ivaryobj(l)
      enddo
      if(m.lt.nspec_fit) then
         write(0,*)'Too few magnitudes to calculate K-Corrections (<'
     *        ,nspec_fit,').'
         do j=1,nchan
            magmod(j)  = -99.
            magcorr(j) = -99.
         enddo
         do l=1,nspec
            comp(l) = -1
         enddo
         goto 500
      endif

      if(op.eq.1) then
c     If op = 1, transform from magnitudes to fluxes.
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
c     If op = 0, data is in flux.
         do jchan = 1,nchan
            jy(jchan)  = mag(jchan)
            if(jyuse(jchan).eq.2) jy(jchan) = 0.d0
            ejy(jchan) = emag(jchan)**2
         enddo
      else
         write(0,*)'Not a valid value for op in function kca.'
         write(0,*)'Input 0 if data is in Jy and 1 if data is in
     *        magnitudes'
         stop
      endif

c     Build the weights
      do jchan=1,nchan
         do kwave=1,nwave
            wgt(jchan,kwave) = getweight(z,jchan,kwave)
         enddo
         call getrange(jchan)
      enddo

c     Now Fit the Model Fluxes to the spectra. Allow for fixed
c     components with a provided amplitude.

      vecfac = DL(z)**2*1d10*3d-9/(1.d0+z)
      do l = 1,nspec
         if(ivaryobj(l).eq.1) then
            vec(l) = 0.d0
         else
            vec(l) = comp(l)*alpha_norm(l)/vecfac
         endif
      enddo
c      do l=1,nspec
c         vec(l) = 0.d0
c      enddo
c      do l = 1,nspec
c         ivaryobj(l) = 1
c      enddo
      call kca_fitgal

c     Calculate the models at redshift z and redshift 0
      do jchan=1,nchan
         do kwave=1,nwave
            wgt(jchan,kwave) = getweight(z0,jchan,kwave)
         enddo
         call getrange(jchan)
      enddo

      do k = 1,nwave
         tigm(k) = transmit(bcen(k),z,igm)
      enddo
      do l=1,nspec
         do j=1,nchan
            jymodz0(l,j) = 0.d0
            do k=jwmin(j),jwmax(j)
               if(l.ne.1) then
                  dust = 1.d0
               else
                  dust = 10.d0**(-0.4d0*tau(k)*ebv)
               endif
               jymodz0(l,j) = jymodz0(l,j) + c(j)*spec(l,k)*
     *              wgt(j,k)*dust*tigm(k)
            enddo
         enddo
      enddo


      do j=1,nchan
         jymodtot(j)   = 0.d0
         jymodz0tot(j) = 0.d0
         do l = 1,nspec
            jymodtot(j)   = jymodtot(j) + vec(l)*jymod(l,j)
            jymodz0tot(j) = jymodz0tot(j) + vec(l)*jymodz0(l,j)
         enddo
         jycorr(j) = jymodtot(j)/jymodz0tot(j)*(1.d0+z)/(1.d0+z0)
      enddo


      if(op.eq.1) then
c     Write the results in magnitudes if op==1.
         do j = 1,nchan
            magmod(j)  = -2.5*log10(jymodtot(j)/jyzero(j))
            magcorr(j) = -2.5*log10(jycorr(j))
         enddo
      else if(op.eq.0) then
         do j = 1,nchan
            magmod(j)  = jymodtot(j)
            magcorr(j) = jycorr(j)
         enddo
      endif


c     Copy the components of vector vec to the output comp and scale
c     them to specific luminosities units.
      vecfac = DL(z)**2*1d10*3d-9/(1.d0+z)
      do l = 1,nspec
         comp(l) = vec(l)*vecfac/alpha_norm(l)
      enddo
      do l1=1,nspec
         do l2=1,nspec
            covx(l1,l2) = cov(l1,l2)*vecfac**2/
     *           (alpha_norm(l1)*alpha_norm(l2))
         enddo
      enddo

      xebv = ebv
      xigm = igm
      chi2 = chimin

 500  continue

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Determine the fit coefficients for individual galaxies.
      subroutine kca_fitgal
      implicit real*8 (a-h,o-z)
      parameter (NCMAX=40,NGMAX=17000,NWMAX=350,NSMAX=4,NTMAX=4)

      real*8 z
      real*8 jy(NCMAX),ejy(NCMAX)
      integer nchan
      common /data1b/jy,ejy,nchan
      common /data1/z

      integer jyuse(NCMAX)
      common /data2/jyuse

      real*8 vec(NSMAX)
      real*8 jymod(NSMAX,NCMAX)
      real*8 jymodtot(NCMAX)
      common /models/jymod,jymodtot,vec

      real*8 spec(NSMAX,NWMAX),specuse(NSMAX,NWMAX)
      common /specmod1/spec,specuse,nspec
      common /specnorm/bminnorm,bmaxnorm

      real*8 wgt(NCMAX,NWMAX)
      real*8 c(NCMAX)
      common /weights1/wgt,c
      integer jwmin(NCMAX),jwmax(NCMAX)
      common /weights2/jwmin,jwmax

      integer ivaryobj(NSMAX)
      common /ivary/ivaryobj

      real*8 bedge(NWMAX)
      real*8 bcen(NWMAX)
      common /wavegrid/bedge,bcen,nwave

      real*8 a(10,10),b(10)
      real*8 atemp(10,10),btemp(10)
      real*8 asave(10,10),bsave(10)
      real*8 tempw(10),tempv(10,10),tempz(10)
      integer itempv(10)
      real*8 temps(10),temps2(10)

      real*8 chitabebv(100),chitabigm(100)
      real*8 tabebv(100),tabigm(100)

      real*8 jymodx(NSMAX,NCMAX)
      common /modelsx/jymodx

      real*8 vecbest(NSMAX)
      real*8 jymodbest(NSMAX,NCMAX)
      real*8 jymodtotbest(NCMAX)

      real*8 tau(NWMAX),ebv,igm
      common /dust/tau,ebv,igm

      real*8 emin,emax,de
      real*8 gmin,gmax,dg
      integer ne,ng
      common /redpars/emin,emax,de,ne
      common /igmpars/gmin,gmax,dg,ng

      real*8 cov(NSMAX,NSMAX)
      common /vecerr/cov

      real*8 tigm(NWMAX)

      real*8 chimin
      common /minchi/chimin

      real*8 work(500)
      integer ipiv(10)

      integer use_red_igm_prior
      common /red_igm_prior/use_red_igm_prior


      chimin = 1.d32

c     Set some needed reddening parameters.
      ng0 = ng
      do j=1,ne
         chitabebv(j) = 1.d32
      enddo
      do j=1,ng
         chitabigm(j) = 1.d32
      enddo

c     Start the main cycle.
      do ie = 1,ne+1
         ngx = 1
         if(ie.le.ne) ngx = ng0
         do ig = 1,ngx
c     If not on the last step, go through the grid. Otherwise, fit for
c     the best value of E(B-V) and IGM.
            if(ie.le.ne) then
               if(ie.eq.1) then
                  euse = 0.d0
               else
                  euse = emin + de*float(ie-2)
                  euse = 10.d0**euse
               endif
               guse = gmin + dg*float(ig-1)
            else
               euse = tabebv(iebst)
               if((iebst.ne.1).and.(iebst.ne.ne)) then
                  y1 = chitabebv(iebst-1)
                  y2 = chitabebv(iebst  )
                  y3 = chitabebv(iebst+1)
                  x1 = tabebv(iebst-1)
                  x2 = tabebv(iebst  )
                  x3 = tabebv(iebst+1)
                  aa = ((y3-y2)*(x2-x1) - (y2-y1)*(x3-x2))/
     *                 ((x3**2-x2**2)*(x2-x1)-(x2**2-x1**2)*(x3-x2))
                  bb = ((y3-y2)-aa*(x3**2-x2**2))/(x3-x2)
                  euse = -bb/(2.d0*aa)
                  if((euse.lt.tabebv(iebst-1)).or.(euse.gt.tabebv(iebst+1)))
     *                 euse = tabebv(iebst)
               endif
               guse = tabigm(igbst)
               if((igbst.ne.1).and.(igbst.ne.ng)) then
c     For the IGM it is OK to do the quadratic interpolation this way,
c     as it is uniformly sampled in linear space.
                  den = 2.d0*chitabigm(igbst)-chitabigm(igbst-1)-
     *                 chitabigm(igbst+1)
                  if(den.lt.0.d0) then
                     guse = tabigm(igbst) + 0.5*dg*
     *                    (chitabigm(igbst+1)-chitabigm(igbst-1))/den
                     if((guse.lt.tabigm(igbst-1)).or.
     *                    (guse.gt.tabigm(igbst+1)))
     *                    guse = tabigm(igbst)
                  endif
               endif
            endif

            tabebv(ie) = euse
            tabigm(ig) = guse
c     This prior makes the IGM and reddening be as little as possible
c     without altering the fit.
            if(use_red_igm_prior.eq.1) then
               chi        = (euse/0.5d0)**2 + ((guse-1.d0)/0.5d0)**2
            else
               chi = 0.d0
            endif

c     Work out the contribution from each template to the object
            do k=1,nwave
               tigm(k) = transmit(bcen(k),z,guse)
            enddo
            do l=1,nspec
               do j=1,nchan
                  jymod(l,j) = 0.d0
                  do k=jwmin(j),jwmax(j)
                     if(l.ne.1) then
                        dust = 1.d0
                     else
                        dust = 10.d0**(-0.4d0*tau(k)*euse)
                     endif
                     jymod(l,j) = jymod(l,j) + c(j)*spec(l,k)*
     *                    wgt(j,k)*dust*tigm(k)
                  enddo
               enddo
            enddo

c     Compute the present model
            maxdim = 10
            call clearmat(atemp,btemp,maxdim,nspec)
            do j=1,nchan
               if (jyuse(j).ge.1) then
                  do l1=1,nspec
                     btemp(l1) = btemp(l1) + jy(j)*jymod(l1,j)/ejy(j)
                     do l2=l1,nspec
                        atemp(l1,l2) = atemp(l1,l2) +
     *                       jymod(l1,j)*jymod(l2,j)/ejy(j)
                     enddo
                  enddo
               endif
            enddo
            call symmat(atemp,btemp,maxdim,nspec)


c     Having built the matrix assuming everything is varying, rearrange
c     the equations for when some are held fixed
            nm1 = 0
            do l1=1,nspec
               if (ivaryobj(l1).eq.1) then
                  nm1    = nm1 + 1
                  nm2    = 0
                  b(nm1) = btemp(l1)
                  do l2=1,nspec
                     if (ivaryobj(l2).eq.1) then
                        nm2        = nm2 + 1
                        a(nm1,nm2) = atemp(l1,l2)
                     else
                        b(nm1) = b(nm1) - vec(l2)*atemp(l1,l2)
                     endif
                  enddo
               endif
            enddo

c     Save the matrices
            do l1=1,nm1
               bsave(l1) = b(l1)
               do l2=1,nm1
                  asave(l1,l2) = a(l1,l2)
               enddo
            enddo

c     Solve assuming only positive coefficients. If convergence fails,
c     revert to the slower version going through all possible
c     combinations.
            call my_nnls_2(a,maxdim,nm1,nm1,b,temps,MODE,its,0)
            if(MODE.eq.3) then
               do l1=1,nm1
                  b(l1) = bsave(l1)
                  do l2=1,nm2
                     a(l1,l2) = asave(l1,l2)
                  enddo
               enddo
               nm3 = 0
               do l=1,nspec
                  if(ivaryobj(l).eq.1) then
                     nm3 = nm3 + 1
                     do j=1,nchan
                        jymodx(nm3,j) = jymod(l,j)
                     enddo
                  endif
               enddo
               call ANNLS(a,maxdim,nm1,nm1,b,temps)
            endif

c     Copy solution out into final vector
            nm1  = 0
            do l1=1,nspec
               if (ivaryobj(l1).eq.1) then
                  nm1     = nm1 + 1
                  vec(l1) = temps(nm1)
               endif
            enddo

c     Calculate the chi-square of the fit.
            do j=1,nchan
               jymodtot(j) = 0.d0
               do l=1,nspec
                  jymodtot(j) = jymodtot(j) + vec(l)*jymod(l,j)
               enddo
               if(jyuse(j).ge.1) then
                  diff = jy(j)-jymodtot(j)
                  chi  = chi + diff*diff/ejy(j)
               endif
            enddo

            if(chi.le.chitabebv(ie)) chitabebv(ie) = chi
            if(chi.le.chitabigm(ig)) chitabigm(ig) = chi
            if(chi.le.chimin) then
               chimin = chi
               ebv    = euse
               igm    = guse
               iebst  = ie
               igbst  = ig
               do l=1,nspec
                  vecbest(l) = vec(l)
                  do j=1,nchan
                     jymodbest(l,j) = jymod(l,j)
                  enddo
               enddo
               do j=1,nchan
                  jymodtotbest(j) = jymodtot(j)
               enddo
c     Estimate the template errors.
               nm1 = 0
               do l=1,nspec
                  nm1 = nm1 + ivaryobj(l)
               enddo
               lwork = 500
               call dgetrf(nm1,nm1,asave,maxdim,ipiv,INFO)
               call dgetri(nm1,asave,maxdim,ipiv,work,lwork,INFO)
               nm1 = 0.d0
               do l1=1,nspec
                  nm2 = 0
                  if(ivaryobj(l1).eq.1) then
                     nm1 = nm1 + 1
                     do l2=1,nspec
                        if(ivaryobj(l2).eq.1) then
                           nm2 = nm2 + 1
                           cov(l1,l2) = asave(nm1,nm2)
                        else
                           cov(l1,l2) = 0.d0
                        endif
                     enddo
                  else
                     do l2=1,nspec
                        cov(l1,l2) = 0.d0
                     enddo
                  endif
               enddo
            endif

         enddo
      enddo

c     Now that we've finished the main cycle, get the best fit
c     values. Notice that this is not immediate from the last cycle, as
c     it is possible for the interpolation scheme to give a larger chi2
c     as the surface is not necessarily a good paraboloid.
      do l=1,nspec
         vec(l) = vecbest(l)
         do j=1,nchan
            jymod(l,j) =jymodbest(l,j)
         enddo
      enddo
      do j=1,nchan
         jymodtot(j) = jymodtotbest(j)
      enddo

c     Finally remove the prior term from the chi2.
      if(use_red_igm_prior.eq.1) then
         chimin = chimin - ((ebv/0.5d0)**2 + ((igm-1.d0)/0.5d0)**2)
      endif

      return
      end
