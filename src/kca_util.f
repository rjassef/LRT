c     This function calculates the maximum volume on which a galaxy
c     could have been found. Area is the area of the survey area in
c     square degrees.
c
      real*8 function vmax(comp,ebv,igm,jchan,mlim,zlim,area)
      implicit real*8 (a-h,o-z)
      parameter (NCMAX=32,NWMAX=350,NSMAX=4)
      
      real*8 comp(NSMAX)
      real*8 mlim
      real*8 ebv,igm

c     Find the maximum redshift to which this galaxy could have been
c     found.
      call findzlim(comp,ebv,igm,jchan,mlim,zlim,0)
c     Now calculate vmax.
      domega = area/41253.d0
      vmax   = domega*vc(zlim)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     This Function determines the maximum redshift up to which a galaxy
c     could be found, given its best fit SED. Important for VMAX
c     calculations. This function is called directly by the vmax
c     function. Note that the maximum redshift to which this algorithm
c     runs is 3.0. Please modify the value of zlimhig if you need higher
c     redshifts. This function works by bracketing the maximum
c     redshift. If the maximum redshift is higher than three or if the
c     algorithm doesn't converge, the function will return zlim=4.
c
      subroutine findzlim(comp,ebvx,igmx,jchan,mlim,zlim,idbg)
      implicit real*8 (a-h,o-z)
      parameter (NCMAX=32,NWMAX=350,NSMAX=4)

      real*8 mlim,magl
      real*8 ebvx,igmx
      real*8 vec(NSMAX),comp(*)

      real*8 alpha(NSMAX)
      common /alphanorm/alpha

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

      real*8 jyzero(NCMAX),con(NCMAX),lbar(NCMAX)
      common /cal1/jyzero,con,lbar

      real*8 jymodlim(NSMAX),jymodlimtot

      real*8 tau(NWMAX),ebv,igm
      common /dust/tau,ebv,igm

      real*8 tigm(NWMAX)

      real*8 emin,emax,de
      real*8 gmin,gmax,dg
      integer ne,ng
      common /redpars/emin,emax,de,ne     
      common /igmpars/gmin,gmax,dg,ng


c     Initialize the range on which a galaxy can be found or the value
c     returned in case of failure.
      zlimlow      = 0.0d0
      zlimhig      = 7.0d0
      zlim_default = -4.0d0
      icont = 0

c     Initiate main cycle.
 300  continue
      zlim    = 0.5*(zlimhig+zlimlow)

c     Put the vectors in the internal units of the program.
      vecfac = DL(zlim)**2*1d10*3d-9/(1.d0+zlim)
      do l = 1,nspec
         vec(l) = comp(l)*alpha(l)/vecfac
      enddo

c     Build the weights.
      do kwave=1,nwave
         wgt(jchan,kwave) = getweight(zlim,jchan,kwave)
      enddo
      call getrange(jchan)

c     Get the model flux.
      do kwave=1,nwave
         tigm(kwave) = transmit(bcen(kwave),zlim,igm)
      enddo
      do l=1,nspec
         jymodlim(l) = 0.d0
         do k=jwmin(jchan),jwmax(jchan)
            if(l.ne.1) then
               dust = 1.d0
            else
               dust = 10.d0**(-0.4d0*tau(k)*euse)
            endif
            jymodlim(l) = jymodlim(l) + 
     *           c(jchan)*spec(l,k)*wgt(jchan,k)*dust*tigm(k)
         enddo
      enddo
      jymodlimtot = 0.d0
      do l = 1,nspec
         jymodlimtot = jymodlimtot + 
     *        vec(l)*jymodlim(l)
      enddo

c     Get the magnitude, compare it to the limit magnitude and decide
c     how to bracket.
      magl = -2.5*log10(jymodlimtot/jyzero(jchan))
      if(idbg.eq.1) then
         print*,zlim,zlimhig,zlimlow,magl,mlim
         pause
      endif
      if(magl.le.mlim) zlimlow = zlim
      if(magl.gt.mlim) zlimhig = zlim

c     Convergence criteria.
      if(abs(magl-mlim).lt.0.002d0) goto 301
      if(abs(zlimhig-zlimlow).lt.0.00001d0.and.zlimhig.eq.3.0d0) then
         zlim = zlim_default
         goto 301
      endif

c     If it cannot converge, print this statement and exit the
c     subroutine.
      icont = icont + 1
      if(icont.gt.100) then
         write(0,*)'Could not converge in 100 iterations in function 
     *        findzlim. Exiting subroutine and setting zlim to -4.'
c         write(0,*)zlim,zlimhig,zlimlow,magl,mlim
         pause
         zlim = zlim_default
         goto 301
      endif

      goto 300
c     Continue here when exiting the cycle.
 301  continue

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This function returns the bolometric luminosity of a galaxy. This
c     is calculated just by adding the specific luminosities of each
c     component.
c
      real*8 function bol_lum(vec)
      implicit real*8 (a-h,o-z)
      parameter (NCMAX=32,NWMAX=350,NSMAX=4)

      real*8 vec(*)

      real*8 spec(NSMAX,NWMAX),specuse(NSMAX,NWMAX)
      common /specmod1/spec,specuse,nspec

      bol_lum = 0.0d0
      do l = 1,nspec
         bol_lum = bol_lum + vec(l)
      enddo

      return 
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      real*8 function err_lum(cov)
      implicit real*8 (a-h,o-z)
      parameter(NSMAX=4,NWMAX=350)

      real*8 cov(NSMAX,NSMAX)

      real*8 spec(NSMAX,NWMAX),specuse(NSMAX,NWMAX)
      common /specmod1/spec,specuse,nspec


      err_lum = 0.d0
      do l1=1,nspec
         err_lum = err_lum + cov(l1,l1)
         do l2=l1+1,nspec
            err_lum = err_lum + 2.d0*cov(l1,l2)
         enddo
      enddo
      err_lum = dsqrt(err_lum)

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      real*8 function err_lum_host(cov)
      implicit real*8 (a-h,o-z)
      parameter(NSMAX=4,NWMAX=350)

      real*8 cov(NSMAX,NSMAX)

      real*8 spec(NSMAX,NWMAX),specuse(NSMAX,NWMAX)
      common /specmod1/spec,specuse,nspec

      real*8 lum


      err_lum = 0.d0
      do l1=2,nspec
         err_lum = err_lum + cov(l1,l1)
         do l2=l1+1,nspec
            err_lum = err_lum + 2.d0*cov(l1,l2)
         enddo
      enddo
      err_lum_host = sqrt(err_lum)

      return
      end

