c     This function returns the luminosity distance given a
c     redshift. This function should only be called after initializing
c     KCA or PZA.
c
      real*8 function DL(z)
      implicit real*8 (a-h,o-z)

      real*8 dcom(5000),dlum(5000),dmod(5000),dvol(5000),rstar(5000)
      common /distance/dcom,dlum,dmod,dvol,rstar

      real*8 zin2(5000)
      common /distance2/zin2,dz,zinmin,zinmax,nin2

      real*8 om,ol,ok,H0,DH
      common /cosmo/om,ol,ok,H0,DH

      real*8 DC,DM

c     Calculate the position of zin that holds the closest, but lower,
c     value of z. Also check that everything is ok.
      i = (z-zinmin)/dz+1
      if(i.gt.nin2) i=nin2

c     Complete the integration up to the redshift needed.
      npt = 1000
      zmin = zin2(i)
      zmax = z
      delz = (zmax-zmin)/float(npt-1)
      dadd = 0.d0
      do k=1,npt
         ztemp = zmin + delz*float(k-1)
         x   = ztemp
         val = 1.d0/sqrt(om*(1.d0+x)**3+ok*(1.d0+x)**2+ol)
         if (k.ne.1) dadd = dadd + 0.5d0*delz*(val+vold)
         vold  = val
      enddo
      DC = dcom(i) + DH*dadd

c     Calculate DM and DL.
      if(ok.gt.0.d0) then
         DM = DH*sinh(sqrt(ok)*DC/DH)/sqrt(ok)
      else if(ok.lt.0d0) then
         DM = DH*sin(sqrt(-ok)*DC/DH)/sqrt(-ok)
      else
         DM = DC
      endif
      DL = (1.d0+z)*DM

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     This function returns the comoving volume given a redshift.
c
      real*8 function vc(z)
      implicit real*8 (a-h,o-z)

c     If argument is negative, return -1.
      if(z.lt.0.d0) then
         vc = -1.d0
      else
         vc = (4.d0*3.14159d0/3.d0)*DL(z)**3/(1.d0+z)**3
      endif

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Given a certain vector comp, holding the specific luminosities of
c     each component (just like the one returned by the kca function),
c     this subroutine returns the fluxes you would expect at a given
c     redshift. This function if for creating mock galaxies using the
c     templates.
      subroutine get_mags(comp,euse,guse,z,jymodtot)
      implicit real*8 (a-h,o-z)
      parameter (NCMAX=40,NWMAX=350,NSMAX=4)

      real*8 comp(*),jymodtot(*)
      real*8 euse,guse
      real*8 vec(NSMAX)

      real*8 alpha_norm(NSMAX)
      common /alphanorm/alpha_norm

      real*8 wgt(NCMAX,NWMAX)
      real*8 c(NCMAX)
      common /weights1/wgt,c
      integer jwmin(NCMAX),jwmax(NCMAX)
      common /weights2/jwmin,jwmax

      real*8 bedge(NWMAX)
      real*8 bcen(NWMAX)
      common /wavegrid/bedge,bcen,nwave

      real*8 spec(NSMAX,NWMAX),specuse(NSMAX,NWMAX)
      common /specmod1/spec,specuse,nspec

      real*8 jyzero(NCMAX),sat(NCMAX),con(NCMAX),lbar(NCMAX)
      common /cal1/jyzero,sat,con,lbar

      real*8 jy(NCMAX),ejy(NCMAX)
      integer nchan
      common /data1b/jy,ejy,nchan

      real*8 jymod(NSMAX,NCMAX)

      real*8 tau(NWMAX),ebv,igm
      common /dust/tau,ebv,igm

      real*8 tigm(NWMAX)


c     Transform the component vectors back to the natural units of the
c     code. Note DL(z) is in Mpc.
      if(z.gt.0.d0) then
         DL_use = DL(z)/(1.d0+z)
      else if(z.eq.-1.d0) then
c     Assume we are getting an absolute magnitude if z=-1. So DL=10pc=1e-5Mpc
         DL_use = 1.d-5
         z=0
      else
         write(6,*)'Redshift must be positive.'
      endif
      vecfac = DL_use**2*1d10*3d-9
      do l = 1,nspec
         vec(l) = comp(l)*alpha_norm(l)/vecfac
      enddo

      do jchan=1,nchan
         do kwave=1,nwave
            wgt(jchan,kwave) = getweight(z,jchan,kwave)
         enddo
         call getrange(jchan)
      enddo

      do kwave=1,nwave
         tigm(kwave) = transmit(bcen(kwave),z,guse)
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
               jymod(l,j) = jymod(l,j) + c(j)*spec(l,k)*wgt(j,k)*
     *              dust*tigm(k)
            enddo
         enddo
      enddo

      do j=1,nchan
         jymodtot(j)   = 0.d0
         do l = 1,nspec
            jymodtot(j)   = jymodtot(j) + vec(l)*jymod(l,j)
         enddo
      enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Clear the matrix
      subroutine clearmat(a,b,maxdim,nwave)
      real*8 b(maxdim),a(maxdim,maxdim)

      do k1=1,nwave
         b(k1) = 0.0
         do k2=1,nwave
            a(k1,k2) = 0.0
         enddo
      enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Symmetrize the matrix
      subroutine symmat(a,b,maxdim,nwave)
      real*8 b(maxdim),a(maxdim,maxdim)

      do k1=1,nwave
         do k2=1,k1-1
            a(k1,k2) = a(k2,k1)
         enddo
      enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Get chi2 of the fit using fluxes.
c
      subroutine get_chi_jy(jy,ejy,jymod,jyuse,nchan,chi2)
      implicit real*8 (a-h,o-z)

      real*8 jy(*),jymod(*),ejy(*)
      integer jyuse(*)


      chi2 = 0.d0

      do j=1,nchan
         if(jyuse(j).ge.1) then
            chi2 = chi2 + ((jy(j)-jymod(j))/ejy(j))**2
         endif
      enddo

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c     Get chi2 of the fit using mags.
c
      subroutine get_chimag(mag,emag,magmod,maguse,nchan,chi2)
      implicit real*8 (a-h,o-z)
      parameter (NCMAX=40)

      real*8 mag(*),magmod(*),emag(*)
      integer maguse(*)

      real*8 jy(NCMAX),ejy(NCMAX),jymod(NCMAX)

      real*8 jyzero(NCMAX),sat(NCMAX),con(NCMAX),lbar(NCMAX)
      common /cal1/jyzero,sat,con,lbar


      do jchan=1,nchan
         jy(jchan)   = jyzero(jchan)*10.d0**(-0.4d0*mag(jchan))
         ejy(jchan)  = jy(jchan)*emag(jchan)*2.5d0/dlog(10.d0)
         jymod(jchan)= jyzero(jchan)*10**(-0.4*magmod(jchan))
      enddo


      chi2 = 0.d0

      do j=1,nchan
         if(maguse(j).eq.1) then
            chi2 = chi2 + ((jy(j)-jymod(j))/ejy(j))**2
         endif
      enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine noblank(fname,i,j)
      implicit real*8 (a-h,o-z)

      character fname*(*)

      do i=1,100
         if(fname(i:i).ne.' ') goto 100
      enddo
 100  continue
      do j=100,1,-1
         if(fname(j:j).ne.' ') goto 110
      enddo
 110  continue

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Get the AGN Luminosity density at a specific wavelength in um

      real*8 function get_Lnu(comp,z,xlam)
      implicit real*8(a-h,o-z)
      parameter (NCMAX=40,NWMAX=350,NSMAX=4,NTMAX=4)

      real*8 comp(*)

      real*8 bedge(NWMAX)
      real*8 bcen(NWMAX)
      common /wavegrid/bedge,bcen,nwave

      real*8 spec(NSMAX,NWMAX),specuse(NSMAX,NWMAX)
      common /specmod1/spec,specuse,nspec
      common /specnorm/bminnorm,bmaxnorm

      real*8 vec(NSMAX)

      real*8 alpha_norm(NSMAX)
      common /alphanorm/alpha_norm

c     Comp to vec
      vecfac = DL(z)**2*1d10*3d-9/(1.d0+z)
      do l = 1,nspec
         vec(l) = comp(l)*alpha_norm(l)/vecfac
      enddo

      pi = 4.d0*datan(1.d0)
      flux = 0.d0
      dflux = 0.d0
      ksave = -1
      do k=1,nwave
         if(bcen(k).ge.xlam) then
            ksave = k
            goto 70
         endif
      enddo
 70   continue
      suse = (spec(1,ksave)-spec(1,ksave-1))/(bcen(ksave)-bcen(ksave-1))
     *     * (xlam-bcen(ksave)) + spec(1,ksave)
      flux  = suse*vec(1)
      flux = flux*1.d-23
      get_Lnu = 4.d0*pi*(DL(z)*3.086d24)**2 * flux/(1.d0+z)

      return
      end

cccccccccccc
c
c     Mass to Light ratios estimated from g-r color and K band Luminosity
c     calibrations c of Bell et al. (2003). A correction of -0.15dex is taken
c     to change the IMF from "diet" Salpeter to Kroupa, as the latter should
c     be a description of galaxies.
c
      real*8 function ML_Bell03_II(comp,z)
      implicit real*8 (a-h,o-z)
      parameter(NWMAX=350,NCMAX=40,NSMAX=4)
      parameter(pi=3.14159d0)

      real*8 M
      real*8 comp(*)
      real*8 vec(NSMAX)

      real*8 alpha_norm(NSMAX)
      common /alphanorm/alpha_norm

c     Precomputed. Holds templates convolved with Ks band.
      real*8 jymod_ml_kband(NSMAX)
      data jymod_ml_kband/0.d0,1.7223013079080471d0,
     *     1.0651971527188366d0,0.91153915287120157d0/
      real*8 jyzero_kband

      real*8 spec(NSMAX,NWMAX),specuse(NSMAX,NWMAX)
      common /specmod1/spec,specuse,nspec
      common /specnorm/bminnorm,bmaxnorm

c     Precomputed. Holds the K-band mass to light ratio for each
c     template.
      real*8 ML(NSMAX)
      data ML/0.d0,0.62351559837353399d0,
     *     0.59180425337814679d0,0.47329597624195591d0/

c     Ks-band zero point.
      jyzero_kband = 666.7d0


c     Comp to vec
c      vecfac = DL(z)**2*1d10*3d-9/(1.d0+z)
      vecfac = 1.d10*3d-9
      do l = 1,nspec
         vec(l) = comp(l)*alpha_norm(l)/vecfac
      enddo

      M = 0.d0
      do l=2,nspec
         flux = vec(l)*jymod_ml_kband(l)
         xmag = -2.5d0*dlog10(flux/jyzero_kband) - 25.d0
         flum = 10.d0**(-0.4d0*(xmag-3.32d0))
         M = M + flum*ML(l)
      enddo

      ML_Bell03_II = M
      return
      end
