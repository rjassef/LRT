c     Subroutine to estimate zero point corrections given a training
c     set. The training set must be file having in each line the
c     redshift of a galaxy, all magnitudes (or fluxes), all errors and
c     all use flags (1 to use the band and zero to not). See the README
c     file or the online manual
c     (http://www.astronomy.ohio-state.edu/~rjassef/instructions.html)
c     for a more thorough explanation.

      subroutine fitzero(filename,niter2,chifrac,corr,op)
      implicit real*8 (a-h,o-z)
      parameter (NCMAX=32,NWMAX=350,NSMAX=4,NTMAX=4,NGMAX=17000)

      real*8 corr(*)
      integer op   
      character filename*(*)
      character*1000 line
      real*8 jya(NGMAX,NCMAX),ejya(NGMAX,NCMAX)
      integer jyusea(NGMAX,NCMAX)
      real*8 z(NGMAX)      
      real*8 wgta(NGMAX,NCMAX,NWMAX)
      real*8 jwmina(NGMAX,NCMAX),jwmaxa(NGMAX,NCMAX)
      real*8 dchi(NGMAX),chitemp(NGMAX)
      integer dchiuse(NGMAX)
      real*8 termu(NCMAX),terml(NCMAX)
      real*8 vecold(NGMAX,NSMAX)
      integer verbose
      common /verb/verbose
      integer pzon
      common /regen/zmax,zmin,dz,pzon

ccccccccccccccccccccccccccccccccccccccc
c     Common blocks, usually just used for compatibility issues in this
c     routine.
      real*8 bedge(NWMAX)
      real*8 bcen(NWMAX)
      common /wavegrid/bedge,bcen,nwave

      real*8 wgt(NCMAX,NWMAX)
      real*8 c(NCMAX)
      common /weights1/wgt,c
      integer jwmin(NCMAX),jwmax(NCMAX)
      common /weights2/jwmin,jwmax 

      real*8 jy(NCMAX),ejy(NCMAX)
      integer nchan
      common /data1b/jy,ejy,nchan

      integer jyuse(NCMAX)
      common /data2/jyuse

      real*8 vec(NSMAX)
      real*8 jymod(NSMAX,NCMAX)
      real*8 jymodtot(NCMAX)
      common /models/jymod,jymodtot,vec

      real*8 spec(NSMAX,NWMAX),specuse(NSMAX,NWMAX)
      common /specmod1/spec,specuse,nspec
      common /specnorm/bminnorm,bmaxnorm

      real*8 jyzero(NCMAX),con(NCMAX),lbar(NCMAX)
      common /cal1/jyzero,con,lbar

      integer ivaryobj(NSMAX)
      common /ivary/ivaryobj     

      real*8 tau(NWMAX),ebv,igm
      common /dust/tau,ebv,igm

cccccccccccccccccccccccccccccccccccccccc

      niter = niter2


c     Read in the training set.
      if(verbose.eq.1) then
         print*
         print*,
     *    '*********************************************************'
         print*,'Estimating zero point corrections'
         print*,'Reading training set...'
      endif

      i = 1
      open(unit=11,file=filename,status='old')
 100  read(11,'(a)',end=101)line
         if(line(1:1).eq.'#'.or.line.eq.' ') goto 100
         read(line,*)z(i),(jya(i,j),j=1,nchan),
     *        (ejya(i,j),j=1,nchan),(jyusea(i,j),j=1,nchan)

c     Do not consider galaxies with less or equal number of magnitudes
c     than the number of templates + 2 (IGM + reddening).
         m = 0
         do j=1,nchan
            if(jyusea(i,j).ge.1) m = m + 1
         enddo       
         if(m.gt.nspec+2) i = i+1

         if(i.gt.NGMAX) then
            write(0,*)
            write(0,*)'Aborting program in subroutine fitzero.'
            write(0,*)
            write(0,*)'Too many galaxies in training set.'
            write(0,*)'A maximum of ',NGMAX,' is currently allowed.'
            write(0,*)'Please remove extra galaxies from the' 
            write(0,*)'training set or edit the source file fitzero.f'
            write(0,*)'and adjust the value of NGMAX where needed.'
            write(0,*)
            stop
         endif

         goto 100
 101  continue
      close(11)
      ntarg = i-1

c     Move everything to fluxes if in magnitudes
      if(op.eq.1) then
         do i=1,ntarg
            do j = 1,nchan
               jya(i,j)  = jyzero(j)*10.d0**(-0.4d0*jya(i,j))
               ejya(i,j) = 0.4d0*log(10.d0)*ejya(i,j)*jya(i,j)
            enddo
         enddo
      else if(op.ne.0) then
         write(0,*)'Not a valid value for op in function fitzero.'
         write(0,*)'Input 0 if data is in Jy and 1 if data is in magnitudes
     *        '
         stop
      endif

c     Calculate the weights of each object.
      if(verbose.eq.1) then
         print*,'Building weights...'
      endif
      do i=1,ntarg
         do j=1,nchan
             do k=1,nwave
               wgt(j,k) = getweight(z(i),j,k)
               wgta(i,j,k) = wgt(j,k)
            enddo
            call getrange(j)         
            jwmaxa(i,j) = jwmax(j)
            jwmina(i,j) = jwmin(j)
         enddo
      enddo

c     Set all dchiuse(i) to 1 for first run
      do i=1,ntarg
         dchiuse(i) = 1
         dchi(i)    = 0.d0
      enddo

 333  continue

c     Start the main cycle
      niter = niter + 1
      do kk = 1,niter
         kkk = kk/10.d0
         if(kk.eq.kkk*10.d0.and.verbose.eq.1) print*,'Iteration ',kk

c     Restart all termu and terml
         do j=1,nchan
            termu(j) = 0.d0
            terml(j) = 0.d0
         enddo
         chi2 = 0.d0

         do i=1,ntarg

c     Set all the variables needed for compatibility
            do j=1,nchan
               jy(j)  = jya(i,j)
               ejy(j) = ejya(i,j)**2
               jyuse(j) = jyusea(i,j)
               do k = 1,nwave
                  wgt(j,k) = wgta(i,j,k)
               enddo
               jwmax(j) = jwmaxa(i,j)
               jwmin(j) = jwmina(i,j)
            enddo
            
c     Now Fit the Model Fluxes to the spectra
            do l=1,nspec
               vec(l) = 0.d0
            enddo
            do l = 1,nspec
               ivaryobj(l) = 1
            enddo
            chiold = dchi(i)
            call fitzero_fitgal(kk,isuccess,z(i))

c     Add the flux terms.
            if(dchiuse(i).eq.1) then
               do j=1,nchan
                  if(jyuse(j).ge.1) then
                     termu(j) = termu(j) + jymodtot(j)*jy(j)/ejy(j)
                     terml(j) = terml(j) + jymodtot(j)**2/ejy(j)
                  endif
               enddo
            endif
            
c     Calculate the chi2 for each galaxy
            dchi(i) = 0.d0
            do j=1,nchan
               if(jyuse(j).ge.1) then
                  dchi(i) = dchi(i) + (jy(j)-jymodtot(j))**2/ejy(j)
               endif
            enddo
            if(dchiuse(i).eq.1) chi2 = chi2 + dchi(i)
         enddo
      
c     Calculate the corrections on each band, keeping the first one
c     fixed
         do j=1,nchan
            if(j.ne.1.and.terml(j).gt.0.d0) then
               c(j) = (termu(j)/terml(j))*c(j)
            endif
         enddo

c     Determine which galaxies to use on next run.
         if(kk.lt.niter) then
            do i=1,ntarg
               chitemp(i) = dchi(i)
            enddo
            call heapsort(chitemp,ntarg)
            nlim   = int(chifrac*float(ntarg))
            nlim   = max(1,min(ntarg,nlim))
            chilim = chitemp(nlim)
            do i=1,ntarg
               dchiuse(i) = 1
               if (dchi(i).ge.chilim) dchiuse(i) = 0
            enddo      
         endif

c     If this is the first run, reset the constants
         if(kk.eq.1) then
            do j=1,nchan
               c(j)=1.d0
            enddo
         endif

         write(12,200)(c(j),j=1,nchan),chi2
 200     format(12E20.6)
         
      enddo

      if(verbose.eq.1) then
         print*,'Estimated channel corrections:'
         print*,(c(j),j=1,nchan)
      endif

      do j=1,nchan
         corr(j) = c(j)
      enddo

c     If estimating photometric redshifts, regenerate photometric
c     redshift tables.
      if(pzon.eq.1) then
         if(verbose.eq.1) then
            print*, 'Resetting photometric redshift table'
            print*,
     *       '*********************************************************'
            print*
         endif
         call photoz_grid(zmin,zmax,dz)
      endif

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Determine the fit coefficients for individual galaxies.
      subroutine fitzero_fitgal(niter,isuccess,z)
      implicit real*8 (a-h,o-z)
      parameter (NCMAX=32,NGMAX=17000,NWMAX=350,NSMAX=4,NTMAX=4)
      
      real*8 jy(NCMAX),ejy(NCMAX)
      common /data1b/jy,ejy,nchan
      
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

      real*8 vecold(NSMAX)

      real*8 chitabebv(100),chitabigm(100)
      real*8 tabebv(100),tabigm(100)

      real*8 jymodx(NSMAX,NCMAX)
      common /modelsx/jymodx

      real*8 vecbest(NSMAX)
      real*8 jymodbest(NSMAX,NCMAX)
      real*8 jymodtotbest(NCMAX)

      real*8 tau(NWMAX),ebv,igm
      common /dust/tau,ebv,igm

      real*8 tigm(NWMAX)

      real*8 emin,emax,de
      real*8 gmin,gmax,dg
      integer ne,ng
      common /redpars/emin,emax,de,ne     
      common /igmpars/gmin,gmax,dg,ng


      chimin = 1.d32
      
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
            chi        = (euse/0.5d0)**2 + ((guse-1.d0)/0.5d0)**2

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
                        atemp(l1,l2) = atemp(l1,l2) + jymod(l1,j)*
     *                       jymod(l2,j)/ejy(j) 
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
            nm1 = 0
            do l1=1,nspec
               if(ivaryobj(l1).eq.1) then
                  nm1     = nm1 + 1
                  vec(l1) = temps(nm1)
               endif
            enddo
                  

c     Calculate the chi2 of the solution.
            do j=1,nchan
               jymodtot(j) = 0.d0
               do l=1,nspec
                  jymodtot(j) = jymodtot(j) + vec(l)*jymod(l,j)
               enddo
            enddo
            do j=1,nchan
               if(jyuse(j).ge.1) then
                  chi = chi + (jy(j)-jymodtot(j))**2/ejy(j)
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
               nm1 = 0
               do l=1,nspec
                  nm1 = nm1 + ivaryobj(l)
               enddo
               do j=1,nchan
                  jymodtotbest(j) = jymodtot(j)
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
      chimin = chimin - ((ebv/0.5d0)**2 + ((igm-1.d0)/0.5d0)**2)

      return
      end

