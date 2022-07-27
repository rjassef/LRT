c     Subroutine to initialize the star-fitting routines. This
c     subroutine will initialize the bands and the templates.
c
        subroutine starinit(filtname,inum,verb_flag)
        implicit real*8 (a-h,o-z)

        character filtname*(*)

        integer verbose,verb_flag
        common /verb/verbose

        verbose = verb_flag

c     Initialize the filters.
        call setfilt(filtname)
c     Read the zero point normalizations or otherwise initialize them
        call read_zpc
c     Initialize the templates.
        call settemp_star(inum)

        return
        end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine settemp_star(num)
        implicit real*8 (a-h,o-z)
        parameter (NWMAX=350,NSTMAX=54)

        integer num

        real*8 bedge(NWMAX)
        real*8 bcen(NWMAX)
        common /wavegrid/bedge,bcen,nwave

        real*8 spec_stars(NSTMAX,NWMAX)
        common /specmod_stars/spec_stars,nspec_stars

        integer ivaryobj(4)
        common /ivary/ivaryobj

        character*100 specname
        character*200 path

        integer verbose
        common /verb/verbose

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
c   1 is Main Sequence
c   2 is Giant Stars
c   3 is Super Giant Stars
        if(num.eq.1) then
            write(specname,*)path(ip1:ip2), '/specs/lrt_kc04_MS.dat'
        else if(num.eq.2) then
            write(specname,*)path(ip1:ip2), '/specs/lrt_kc04_GS.dat'
        else if(num.eq.3) then
            write(specname,*)path(ip1:ip2), '/specs/lrt_kc04_SGS.dat'
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
        print*,'*********************************************************'
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
            read(13,*)bt,tt,(spec_stars(ll,kwave),ll=1,nspec)
            if (dabs(bt-bcen(kwave))/bt.gt.0.01d0) then
                write(0,*)'wavelength mismatch in startup file ',bt,bcen(kwave)
                stop
            endif
        enddo
        close(unit=13)
        nspec_stars = nspec

c     Initialize the ivaryobj files. The default is to use all templates.
        do l=1,4
            ivaryobj(l) = 1
        enddo

        if(verbose.eq.1) then
            print*,'Done'
            print*,'*********************************************************'
            print*
        endif

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
        subroutine stf(mag,emag,maguse,magmod,comp,nstar_best,chi2,op)
        implicit real*8(a-h,o-z)
        parameter (NCMAX=40,NWMAX=350,NSTMAX=54)

        real*8 mag(NCMAX),emag(NCMAX),magmod(NCMAX),comp(NCMAX)
        integer maguse(NCMAX),op

        real*8 z
        real*8 jy(NCMAX),ejy(NCMAX)
        integer nchan
        common /data1b/jy,ejy,nchan
        common /data1/z

        integer jyuse(NCMAX)
        common /data2/jyuse

        real*8 vec(2)
        real*8 jymod(NSTMAX,NCMAX)
        real*8 jymodtot(NCMAX)
        common /models_stars/jymod,jymodtot,vec,ns_best

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

        real*8 chimin
        common /minchi/chimin

c   We are fitting stars, so z=0.
        z = 0.d0

c     See which bands will be used for fitting.
        m = 0
        do j = 1,nchan
            if(maguse(j).ne.0.and.maguse(j).ne.1.and.maguse(j).ne.2) then
                write(0,*)'maguse(',j,')=',maguse(j),' not equal to 1 or 0 in function kca.'
                write(0,*)'Aborting program'
                stop
            else
                jyuse(j) = maguse(j)
                if(jyuse(j).ge.1) m = m + 1
            endif
        enddo

c     Since we are fitting two templates at a time, make sure there are at least 2 magnitudes/fluxes.
        nspec_fit = 2
        if(m.lt.nspec_fit) then
            write(0,*)'Too few magnitudes to fit stellar models (<',nspec_fit,').'
            do j=1,nchan
                magmod(j)  = -99.
            enddo
            do l=1,2
                comp(l) = -1.d0
            enddo
            nstar_best = -1
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
            write(0,*)'Input 0 if data is in Jy and 1 if data is in magnitudes'
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
        call run_stellar_fit

        if(op.eq.1) then
c     Write the results in magnitudes if op==1.
            do j = 1,nchan
                magmod(j)  = -2.5*log10(jymodtot(j)/jyzero(j))
            enddo
        else if(op.eq.0) then
            do j = 1,nchan
                magmod(j)  = jymodtot(j)
            enddo
        endif


c     Copy the components of vector vec to the output comp and scale
c     them to specific luminosities units. Assume a distance of 10pc.
        do l = 1,2
            comp(l) = vec(l)
        enddo
        nstar_best = ns_best
        chi2=chimin

500     continue

        return
        end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Determine the fit coefficients for individual galaxies.
        subroutine run_stellar_fit
        implicit real*8 (a-h,o-z)
        parameter (NCMAX=40,NWMAX=350,NSTMAX=54)
        
        real*8 z
        real*8 jy(NCMAX),ejy(NCMAX)
        integer nchan
        common /data1b/jy,ejy,nchan
        common /data1/z      

        integer jyuse(NCMAX)
        common /data2/jyuse

        real*8 vec(2)
        real*8 jymod(NSTMAX,NCMAX)
        real*8 jymodtot(NCMAX)
        common /models_stars/jymod,jymodtot,vec,ns_best

        real*8 jymodx(4,NCMAX)
        common /modelsx/jymodx

        real*8 spec_stars(NSTMAX,NWMAX)
        common /specmod_stars/spec_stars,nspec_stars

        real*8 wgt(NCMAX,NWMAX)
        real*8 c(NCMAX)
        common /weights1/wgt,c
        integer jwmin(NCMAX),jwmax(NCMAX)
        common /weights2/jwmin,jwmax
        
        real*8 bedge(NWMAX)
        real*8 bcen(NWMAX)
        common /wavegrid/bedge,bcen,nwave
        
        real*8 a(10,10),b(10)
        real*8 asave(10,10),bsave(10)
        real*8 temps(10)

        real*8 vecbest(2)
        real*8 jymodtotbest(NCMAX)

        real*8 chimin
        common /minchi/chimin

c     Set a large value for chi2min to start.
        chimin = 1.d32

c     For compatibility+simplicity.
        nspec = nspec_stars

c     Work out the contribution from each template to the object
        do l=1,nspec
            do j=1,nchan
                jymod(l,j) = 0.d0
                do k=jwmin(j),jwmax(j)
                    jymod(l,j) = jymod(l,j) + c(j)*spec_stars(l,k)*
     *  wgt(j,k)
                enddo
            enddo
        enddo

c     Compute the present model. Here, we run in pairs of templates.
        maxdim = 10
        do ns=1,nspec-1

c     Build the matrices.
            nm1 = 2
            call clearmat(a,b,maxdim,nm2)
            do j=1,nchan
                if (jyuse(j).ge.1) then
                    do l1=1,2
                        ll1 = l1+ns-1
                        b(l1) = b(l1) + jy(j)*jymod(ll1,j)/ejy(j) 
                        do l2=l1,2
                            ll2=l2+ns-1
                            a(l1,l2) = a(l1,l2) + jymod(ll1,j)*
     *                                  jymod(ll2,j)/ejy(j) 
                        enddo
                    enddo
                endif
            enddo
            call symmat(a,b,maxdim,nm1)
                
c     Save the matrices
            nm1 = 2
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
                do l=ns,ns+1
                    nm3 = nm3 + 1
                    do j=1,nchan
                        jymodx(nm3,j) = jymod(l,j)
                    enddo
                enddo
                call ANNLS(a,maxdim,nm1,nm1,b,temps)
            endif

c     Copy solution out into final vector
            do l=1,2
                vec(l) = temps(l)
            enddo

c     Calculate the chi-square of the fit.
            chi = 0.d0
            do j=1,nchan
                jymodtot(j) = 0.d0
                do l=1,2
                    jymodtot(j) = jymodtot(j) + vec(l)*jymod(l+ns-1,j)
                enddo
                if(jyuse(j).ge.1) then
                    diff = jy(j)-jymodtot(j)
                    chi  = chi + diff*diff/ejy(j)
                endif
            enddo

            if(chi.le.chimin) then
                chimin = chi
                ns_best = ns
                do l=1,2
                    vecbest(l) = vec(l)
                enddo
                do j=1,nchan
                    jymodtotbest(j) = jymodtot(j)
                enddo           
            endif 
        enddo

c     Now that we've finished the main cycle, get the best fit
c     values. 
        do l=1,2
            vec(l) = vecbest(l)
        enddo
        do j=1,nchan
            jymodtot(j) = jymodtotbest(j)
        enddo

        return
        end

c     Given a certain vector comp, holding the specific luminosities of
c     each component (just like the one returned by the kca function),
c     this subroutine returns the fluxes you would expect at a given
c     redshift. This function if for creating mock galaxies using the
c     templates.
        subroutine get_fluxes_stars(comp,ns_best,jymodtot)
        implicit real*8 (a-h,o-z)
        parameter (NCMAX=40,NWMAX=350,NSTMAX=54)

        real*8 comp(*),jymodtot(*)
        real*8 vec(2)

        real*8 wgt(NCMAX,NWMAX)
        real*8 c(NCMAX)
        common /weights1/wgt,c
        integer jwmin(NCMAX),jwmax(NCMAX)
        common /weights2/jwmin,jwmax

        real*8 bedge(NWMAX)
        real*8 bcen(NWMAX)
        common /wavegrid/bedge,bcen,nwave

        real*8 spec_stars(NSTMAX,NWMAX)
        common /specmod_stars/spec_stars,nspec_stars 


        do l=1,2
            vec(l) = comp(l)
        enddo

        z=0.d0
        do jchan=1,nchan
            do kwave=1,nwave
                wgt(jchan,kwave) = getweight(z,jchan,kwave)
            enddo
            call getrange(jchan)
        enddo

        do j=1,nchan
            jymodtot(j) = 0.d0
            do l=1,2
                ll = l+ns_best-1
                do k=jwmin(j),jwmax(j)
                    jymodtot(j) = jymodtot(j) + vec(l) * 
     *                          c(j)*spec_stars(ll,k)*wgt(j,k)
                enddo
            enddo
        enddo

        return
        end




