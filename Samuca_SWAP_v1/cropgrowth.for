! File VersionID:
!   $Id: cropgrowth.for 197 2011-02-18 10:42:47Z kroes006 $
! ----------------------------------------------------------------------
      subroutine CropGrowth(task) 
! ----------------------------------------------------------------------
!     Date               : Aug 2004   
!     Purpose            : Call proper crop routines for initialization,
!                          calculation of rate/state variables and output
! ----------------------------------------------------------------------

      use variables
      implicit none

      integer task,numcrop
      character messag*200

      goto (1000, 2000, 3000) task

1000  continue

! === initialization =========================================================

      if (flbaresoil) then
! --- set bare soil condition  -----------------------------------------------
        call nocrop (rd,lai,cf,ch,albedo,rsc)
      else
! ---   fixed crop development -----------------------------------------------
        if (croptype(icrop) .eq. 1) call CropFixed(1)
! ---   detailed crop growth -------------------------------------------------
        if (croptype(icrop) .eq. 2) call Wofost(1)
! ---   detailed grass growth  -----------------------------------------------
        if (croptype(icrop) .eq. 3) call Grass(1)
! ---   detailed sugarcane growth --------------------------------------------
        if (croptype(icrop) .eq. 4) call Samuca(1)
      endif

!     error message in case of wofost crop growht and initial data from file 
      if(.not.flbaresoil) then
        if(swinco.eq.3 .and. croptype(icrop).ge.2) then
          write(messag,'(2a)')
     &     'Warning: Wofost crop growth in combination with ',
     &     'initial input from file is not implemented yet !'
          call warn ('CropGrowth',messag,logf,swscre)
        endif
      endif

      return

2000  continue
      
! === calculation of potential crop rate and state variables =================

! --- number of crop
      if (flcropend) then
        numcrop = icrop-1
      else
        numcrop = icrop
      endif

! --- detailed crop growth -------------------------------------------------
      if(numcrop.gt.0) then
        if (croptype(numcrop) .eq. 2) call Wofost(2)
        if (croptype(numcrop) .eq. 3) call Grass(2)
        if (croptype(numcrop) .eq. 4) call Samuca(2)
      endif

      return

3000  continue

! === calculation of actual crop rate and state variables =================

! --- number of crop
      if (flcropend) then
        numcrop = icrop-1
      else
        numcrop = icrop
      endif

      if(numcrop.gt.0) then
! ---   fixed crop development -----------------------------------------------
        if (croptype(numcrop) .eq. 1) call CropFixed(3)
! ---   detailed crop growth -------------------------------------------------
        if (croptype(numcrop) .eq. 2) call Wofost(3)
! ---   detailed grass growth  -----------------------------------------------
        if (croptype(numcrop) .eq. 3) call Grass(3)
! ---   detailed sugarcane growth --------------------------------------------
        if (croptype(icrop) .eq. 4) call Samuca(3)
      endif

      return
      end 

! ----------------------------------------------------------------------
      subroutine astro 
     &   (logf,swscre,daynr,lat,dayl,daylp,sinld,cosld,dsinb,dsinbe,dso)
! ----------------------------------------------------------------------
! Subroutine astro (daynr,lat,dayl,daylp,sinld,cosld)
! Authors: this routine Astro is based on Sastro (Daniel van Kraalingen)
! Date   : 28-November-2005 
! Purpose: This subroutine calculates solar constant, daily
!          extraterrestrial radiation, daylength and some intermediate
!          variables required by other routines. The routine has been
!          rewritten such that latitudes from pole to pole can be used.
!
! Formal parameters:  (I=input,O=output,C=control,IN=init,T=time)
! name   type meaning                                     units  class
! ----   ---- -------                                     -----  -----
! logf    I4  Internal number of logbook output file *.LOG   -      I
! swscre  I4  Switch of screen display:  0 = no display;     -      I
!             1 = summary water balance; 2 = daynumber
! daynr   I4  Day of year (Jan 1st = 1)                      d      I  
! lat     R8  Latitude of the site                       degrees    I  
! dayl    R8  Astronomical daylength (base = 0 degrees)      h      O  
! daylp   R8  Photoperiodic daylength (base = -4 degrees)    h      O  
! sinld   R8  Intermediate variable for other subroutine     -      O  
! cosld   R8  Intermediate variable for other subroutine     -      O  
! dsinb   R8  Daily total of sine of solar height            s      O  
! dsinbe  R8  Daily integral of sine of solar height         s      O  
!             corrected for lower transmission at low                  
!             elevation                                                
! sc      R8  Solar constant at day=daynr                   W/m2     O  
! dso     R8  Daily extraterrestrial radiation            J/m2/d    O  
!                                                                      
! Fatal error checks (on input): lat > 90, lat < -90
! Warnings          : lat above polar circle, lat within polar circle  
! Subprograms called: Warning
! File usage        : none
!----------------------------------------------------------------------
      implicit none
 
!     formal parameters
      integer logf,swscre,daynr
      real*8  lat,dayl,daylp,sinld,cosld,dsinb,dsinbe,dso

!     local parameters
      real*8  angle,aob,dec,pi,rad,zza,zzcos,zzsin,sc,help1
      character messag*200

      data    pi /3.1415926d0/,angle /-4.0d0/
! ----------------------------------------------------------------------
! --- declination of the sun as a function of daynr
!     (see ref.manual: Radiation term: 23.45*rad=0.409 en (90-10)*rad=1.39)
      rad = pi/180.d0
      dec = -asin(sin(23.45d0*rad)*cos(2.d0*pi*dble(daynr+10)/365.0d0))
! --- some intermediate variables
      sinld = sin(rad*lat)*sin(dec)
      cosld = cos(rad*lat)*cos(dec)
      aob = sinld/cosld
! --- calculation of daylenght and photoperiodic daylength
!     solution for polar circle altutude adopted from 
!     Daniel van Kraalingen (routine Sastro, dd 12-june-1996,version 1.1)
      if (aob.lt.-1.0d0) then
        messag = 'Warning: latitude above polar circle, daylength= 0hrs'
        call warn ('Astro',messag,logf,swscre)
        dayl = 0.0d0
        zzcos =  0.0d0
        zzsin =  1.0d0
      else if (aob.gt.1.0d0) then
        messag = 'Warning: latitude within polar circle,daylength=24hrs'
        call warn ('Astro',messag,logf,swscre)
        dayl = 24.0d0
        zzcos =  0.0d0
        zzsin = -1.0d0
      else
        dayl  = 12.0d0*(1.0d0+2.0d0*asin(aob)/pi)
        help1 = (-sin(angle*rad)+sinld)/cosld
        if (help1.gt.1.0d0) then
          daylp = 24.0d0
        else
          daylp = 12.0d0*(1.0d0+2.0d0*asin(help1)/pi)
        endif
!        write(logf,*) 'help1=',help1,'daylp=',daylp
        zza   = pi*(12.0d0+dayl)/24.0d0
        zzcos = cos (zza)
        zzsin = sin (zza)
      endif

!     Daily integral of sine of solar height (DSINB) with a
!     correction for lower atmospheric transmission at lower solar
!     elevations (DSINBE)
      dsinb  = 2.0d0*3600.0d0*(dayl*0.50d0*sinld-12.0d0*cosld*zzcos/pi)
      dsinbe = 2.0d0*3600.0d0*(dayl*(0.50d0*sinld+0.20d0*sinld**2.0d0+
     &      0.10d0*cosld**2.0d0)-(12.0d0*cosld*zzcos+
     &      9.6d0*sinld*cosld*zzcos+2.4d0*cosld**2.0d0*zzcos*zzsin)/pi)

!     Solar constant and daily extraterrestrial radiation
      sc = 1370.0d0*(1.0d0+0.033d0*cos (2.0d0*pi*daynr/365.d0))
      dso  = sc*dsinb

      return
      end
! ----------------------------------------------------------------------
      subroutine cropfixed (task)
! ----------------------------------------------------------------------
!     date               : august 2004                           
!     purpose            : simple crop growth routine for swap 
! ----------------------------------------------------------------------
      use variables
      implicit none

! --- local variables
      integer stepnr,i,icgs,task,lcc,node
      
      real*8  gctb(2*magrs),rdtb(2*magrs),kytb(2*magrs),nihil,
     &        afgen,cptr0,ctr0,cptr(magrs),ctr(magrs),crt(magrs),
     &        rely(magrs),help(magrs),dtsum,dvr,phead                   ! NwRootExtr

      parameter (nihil=1.0d-7)

      save
! ----------------------------------------------------------------------

      goto (1000,2000,3000) task

1000  continue

! === initialization ===================================================

! --- read crop data
      gctb = 0.0d0
      cftb = 0.0d0
      chtb = 0.0d0
      rdtb = 0.0d0
      call readcropfixed (cropfil(icrop),pathcrop,idev,lcc,tsumea,
     & tsumam,tbase,kdif,kdir,gctb,swgc,cftb,swcf,rdtb,rdctb,hlim1,
     & hlim2u,hlim2l,hlim3h,hlim3l,hlim4,rsc,adcrh,adcrl,kytb,
     & cofab,logf,schedule,swinter,pfreetb,pstemtb,scanopytb,
     & avprectb,avevaptb,cumdens,chtb,albedo,swetr,
     & flsolute,ecmax,ecslop,c2eca,c2ecb,c2ecf,numlay,
     & alphacrit,swroottyp,wiltpoint,rootradius,rootcoefa,rsw)           ! NwRootExtr

! --- development stage
      dvs = 0.0d0

! --- initial lai or sc
      lai = afgen (gctb,(2*magrs),dvs)
      if (swgc.eq.2) then
        gc = lai
        lai = lai*3.0d0
      endif

! --- initial crop factor or crop height
      cf = afgen (cftb,(2*magrs),dvs)
      ch = afgen (chtb,(2*magrs),dvs)

! --- actual rooting depth [cm]
      rd = min (rds,afgen (rdtb,(2*magrs),dvs))

! --- initial summation variables of the crop
      cptr0 = 0.
      ctr0 = 0.

! --- init arrays with cum. pot. and act. transpiration 
      cptr = 0.0d0
      ctr = 0.0d0

! --- initialize matric flux potential                                  ! NwRootExtr
      phead = 0.d0   ! dummy                                            !
      if (swroottyp .eq. 2) then                                        !
        node = 1                                                        ! dummy node nr
        call MatricFlux(1,phead,node)                                   !
! ---   output of matric flux potential                                 !
!        if (swmfp.eq.1) call outmatricflux(2,mfp,numnod,tcum,          !
!     &   mflux,z,outfil,pathwork,project,ptra,h)                       !
      endif                                                             ! NwRootExtr

      return          

2000  continue

! === calculate potential rate and state variables ======================

3000  continue

! === calculate actual rate and state variables ======================

! --- increase in temperature sum
      dtsum = max (0.0d0,tav-tbase)

! --- development rate
      if (idev.eq.1) then
        dvr = 2.0/lcc
      elseif (idev.eq.2) then
        if (dvs.lt.1.0d0) then
          dvr = dtsum/tsumea
        else
          dvr = dtsum/tsumam
        endif
      endif

! --- determination of current growing stage
      do i = 1,magrs
        help(i) = kytb(2*i-1)
      end do
      icgs = stepnr(help,magrs,dvs)

! --- water stress
      if(abs(ptra).lt.nihil) then
        reltr = 1.0d0
      else
        reltr = max(min(tra/ptra,1.0d0),0.0d0)
      endif

! ----integrals of the crop --------------------------------------------

! --- phenological development stage
      dvs = dvs+dvr

! --- leaf area index or soil cover fraction    
      lai = afgen (gctb,(2*magrs),dvs)
      if (swgc.eq.2) then
        gc = lai
        lai = lai*3.0d0
      endif

! --- crop factor or crop height
      cf = afgen (cftb,(2*magrs),dvs)
      ch = afgen (chtb,(2*magrs),dvs)

! --- rooting depth [cm]
      rd = min (rds,afgen (rdtb,(2*magrs),dvs))

! --- cumulative relative transpiration, total growing season 
      cptr0 = cptr0 + ptra  
      ctr0 = ctr0  + tra
      if (cptr0.le.nihil) then
        crt0 = 0.0d0
      else
        crt0 = max(min(ctr0/cptr0,1.0d0),0.0d0)
      endif

! --- cumulative relative transpiration, current growing stage
      cptr(icgs) = cptr(icgs) + ptra
      ctr(icgs) = ctr(icgs)  + tra
      if (cptr(icgs).le.nihil) then
        crt(icgs) = 0.0d0
      else
        crt(icgs) = max(min(ctr(icgs)/cptr(icgs),1.0d0),0.0d0)
      endif

! --- relative yield per growing stage and cumulated
      crely = 1.0d0 
      do i = 1,icgs
        rely(i) = 1.0d0-((1.0d0-crt(i))*kytb(2*i))
        crely = crely*rely(i)
      end do

      return
      end

! ----------------------------------------------------------------------
      subroutine cropoutput(task) 
! ----------------------------------------------------------------------
!     Date               : Aug 2004   
!     Purpose            : open and write crop output files 
! ----------------------------------------------------------------------

      use variables
      implicit none

! --- local variables ------------------
      integer task,getun,numcrop
      character messag*200
      character*160 filnam,filtext

      goto (1000, 2000, 3000) task

1000  continue

! === open crop output file and write headers =====================
      
! --- open crop output file
      if (flopencropoutput) then
! ---   open crop output file and write general header
        if (trim(outfil).eq.trim(cropfil(1))) then
          Messag = 'The name of the input crop-file 
     &    ('//trim(cropfil(icrop))//') cannot be equal to the name of '
     &   //'the output crop-file '//trim(outfil)//' Adjust a filename !'
          call fatalerr ('crops',messag)
        endif
        filnam = trim(pathwork)//trim(outfil)//'.crp'
        crp = getun (20,90)
        call fopens(crp,filnam,'new','del')
        filtext = 'output data of simple or detailed crop growth model'
        call writehead (crp,1,filnam,filtext,project)

! ---   write header fixed crop growth
        if (croptype(icrop) .eq. 1) 
     &    call OutCropFixed(1,date,t,daycrop,dvs,lai,cf,rd,
     &          crt0,crely,crp,ch)
! ---   write header detailed crop growth 
        if (croptype(icrop) .eq. 2) 
     &    call OutWofost(1,date,daycrop,crp,t,dvs,lai,cf,rd,ch,
     &                        crt0,crely,crt1,cwdmpot,cwdm,wsopot,wso)
! ---   write header detailed grass growth
        if (croptype(icrop) .eq. 3) 
     &    call OutGrass(1,date,daycrop,crp,t,lai,rd,dvs,cf,ch,
     &                    crt0,crt1,tagppot,tagp,tagptpot,tagpt)

          flopencropoutput = .false.

      else
! ---   header for second and subsequent crops

! ---   write header fixed crop growth
        if (croptype(icrop).eq.1 .and. swheader.eq.1) 
     &    call OutCropFixed(1,date,t,daycrop,dvs,lai,cf,rd,
     &           crt0,crely,crp,ch)

! ---   write header detailed crop growth 
        if (croptype(icrop).eq.2 .and. swheader.eq.1) 
     &    call OutWofost(1,date,daycrop,crp,t,dvs,lai,cf,rd,ch,
     &                        crt0,crely,crt1,cwdmpot,cwdm,wsopot,wso)
! ---   write header detailed grass growth
        if (croptype(icrop).eq.3 .and. swheader.eq.1) 
     &    call OutGrass(1,date,daycrop,crp,t,lai,rd,dvs,cf,ch,
     &                    crt0,crt1,tagppot,tagp,tagptpot,tagpt)

      endif

      return

2000  continue

! --- write actual data ----------------------------------------------------

! --- number of crop
      if (flcropend) then
        numcrop = icrop-1
      else
        numcrop = icrop
      endif

! --- fixed crop file
      if (croptype(numcrop) .eq. 1) 
     &  call OutCropFixed(2,date,t,daycrop,dvs,lai,cf,rd,
     &           crt0,crely,crp,ch)

! --- detailed crop growth 
      if (croptype(numcrop) .eq. 2) 
     &  call OutWofost(2,date,daycrop,crp,t,dvs,lai,cf,rd,ch,
     &                        crt0,crely,crt1,cwdmpot,cwdm,wsopot,wso)

! --- detailed grass growth
      if (croptype(numcrop) .eq. 3) 
     &  call OutGrass(2,date,daycrop,crp,t,lai,rd,dvs,cf,ch,
     &                    crt0,crt1,tagppot,tagp,tagptpot,tagpt)

      return

3000  continue
! --- close crop output file ------------------------------------------------

      close (crp)

      return
      end 

! ----------------------------------------------------------------------
      subroutine nocrop (rd,lai,cf,ch,albedo,rsc)
! ----------------------------------------------------------------------
      real*8  rd,lai,cf,ch,albedo,rsc
! ----------------------------------------------------------------------
      rd = 0.0d0
      lai = 0.0d0
      cf = 0.0d0
      ch = 12.d0
      albedo = 0.23d0
      rsc = 70.d0

      return
      end

! ----------------------------------------------------------------------
      subroutine readcropfixed (crpfil,pathcrop,idev,lcc,tsumea,
     & tsumam,tbase,kdif,kdir,gctb,swgc,cftb,swcf,rdtb,rdctb,hlim1,
     & hlim2u,hlim2l,hlim3h,hlim3l,hlim4,rsc,adcrh,adcrl,kytb,
     & cofab,logf,schedule,swinter,pfreetb,pstemtb,scanopytb,
     & avprectb,avevaptb,cumdens,chtb,albedo,swetr,
     & flsolute,ecmax,ecslop,c2eca,c2ecb,c2ecf,numlay,
     & alphacrit,swroottyp,wiltpoint,rootradius,rootcoefa,rsw)              ! NwRootExtr
! ----------------------------------------------------------------------
!     Update             : July 2009
!     date               : July 2002             
!     purpose            : get crop parameters from cropfile
! ----------------------------------------------------------------------
      implicit  none
      include  'arrays.fi'
      
      integer   crp,idev,i,lcc,swgc,swcf,logf,ifnd,getun2,schedule
      integer   swinter,swetr,numlay,swroottyp                              ! NwRootExtr
      logical   flsolute 
      real*8    gctb(2*magrs),cftb(2*magrs),rdtb(2*magrs),kytb(2*magrs)
      real*8    adcrh,adcrl,tbase,tsumam,tsumea,rdctb(22)
      real*8    hlim1,hlim2u,hlim2l,hlim3h,hlim3l,hlim4,rsc
      real*8    sum,afgen,kdif,kdir,cofab,chtb(2*magrs)
      real*8    tinter(magrs),pfree(magrs),pstem(magrs)
      real*8    scanopy(magrs),avprec(magrs),avevap(magrs)
      real*8    pfreetb(2*magrs),pstemtb(2*magrs),scanopytb(2*magrs)
      real*8    avprectb(2*magrs),avevaptb(2*magrs)
      real*8    depth,rootdis(202),cumdens(202)
      real*8    dvsinput(magrs),cfinput(magrs),chinput(magrs)
      real*8    ecmax,ecslop,c2eca,c2ecb,c2ecf(maho),albedo,alphacrit
      real*8    wiltpoint,rootradius,rootcoefa,rsw                          ! NwRootExtr
      logical   rdinqr
      character crpfil*(*),pathcrop*(*)
! locals
      integer   swc2ecf
      character message*200,filnam*200
! ----------------------------------------------------------------------

! --- initialise and start reading
      filnam = trim(pathcrop)//trim(crpfil)//'.crp'
      crp = getun2 (10,90,2)
      call rdinit(crp,logf,filnam)

! --- phenology
      call rdsinr ('idev',1,2,idev)
      if (idev.eq.1) then
        call rdsinr ('lcc',1,366,lcc)
      elseif (idev.eq.2) then
        call rdsdor ('tsumea',0.0d0,10000.0d0,tsumea)
        call rdsdor ('tsumam',0.0d0,10000.0d0,tsumam)
        call rdsdor ('tbase',-10.0d0, 30.0d0,tbase)
      endif

! --- assimilation                        
      call rdsdor ('kdif',0.0d0,2.0d0,kdif)
      call rdsdor ('kdir',0.0d0,2.0d0,kdir)
     
! --- LAI or soil cover fraction 
      call rdsinr ('swgc',1,2,swgc)
      if (swgc.eq.1) then
        call rdador ('gctb',0.0d0,12.0d0,gctb,(2*magrs),ifnd)
      elseif (swgc.eq.2) then
        call rdador ('gctb',0.0d0,2.0d0,gctb,(2*magrs),ifnd)
      endif

! --- Crop factor or crop height
      call rdsinr ('swcf',1,2,swcf)

! --- check use of crop factors in case of ETref
      if (swetr.eq.1 .and. swcf.eq.2) then
        message = 'If ETref is used (SWETR = 1), always define crop '//
     &           'factors (SWCF = 1)' 
        call fatalerr ('ReadCropFixed',message)
      endif

      if (swcf.eq.1) then
! ---   crop factor is input
        call rdador ('dvs',0.0d0,2.0d0,dvsinput,(magrs),ifnd)
        call rdfdor ('cf',0.0d0,2.0d0,cfinput,(magrs),ifnd)
! ---   store values in cftb
        do i = 1,ifnd
          cftb(i*2) = cfinput(i) 
          cftb(i*2-1) = dvsinput(i)
        enddo
        chtb = -99.99d0
      else
! ---   crop height is input
        call rdador ('dvs',0.0d0,2.0d0,dvsinput,(magrs),ifnd)
        call rdfdor ('ch',0.0d0,1.0d4,chinput,(magrs),ifnd)
! ---   store values in chtb
        do i = 1,ifnd
          chtb(i*2) = chinput(i) 
          chtb(i*2-1) = dvsinput(i)
        enddo
        cftb = -99.99d0
      endif

! --- reflection coefficient and crop resistance
      if (swcf.eq.1) then
! ---   use standard values for ETref
        albedo = 0.23d0
        rsc = 70.0d0
        rsw = 0.0d0
      else
! ---   use crop specific values
        call rdsdor ('albedo',0.0d0,1.0d0,albedo)
        call rdsdor ('rsc',0.0d0,1.0d6,rsc)
        call rdsdor ('rsw',0.0d0,1.0d6,rsw)
      endif

! --- rooting depth
      call rdador ('rdtb',0.0d0,1000.0d0,rdtb,(2*magrs),ifnd)

! --- yield response
      call rdador ('kytb',0.0d0,5.0d0,kytb,(2*magrs),ifnd)

! --- water use
      swroottyp = 1
      if(rdinqr('swroottyp')) then
        call rdsinr ('swroottyp',1,2,swroottyp)                         ! NwRootExtr
      endif
      if (swroottyp.eq.1) then                                          !
         call rdsdor ('hlim1' ,-100.0d0,100.0d0,hlim1)                  !
         call rdsdor ('hlim2u',-1000.0d0,100.0d0,hlim2u)                !
         call rdsdor ('hlim2l',-1000.0d0,100.0d0,hlim2l)                !
         call rdsdor ('hlim3h',-10000.0d0,100.0d0,hlim3h)               !
         call rdsdor ('hlim3l',-10000.0d0,100.0d0,hlim3l)               !
         call rdsdor ('hlim4' ,-16000.0d0,100.0d0,hlim4)                !
         call rdsdor ('adcrh',0.0d0,5.0d0,adcrh)                        !
         call rdsdor ('adcrl',0.0d0,5.0d0,adcrl)                        !
!       Criticial stress index for compensation of root water uptake (-)
        alphacrit = 1.0d0
        if(rdinqr('alphacrit')) then
          call rdsdor ('alphacrit',0.2d0,1.0d0,alphacrit)
        endif
      else                                                              !
        call rdsdor ('wiltpoint',-1.0d6,-1.0d2,wiltpoint)               !
        call rdsdor ('rootradius',0.0001d0,1.0d0,rootradius)            !
        call rdsdor ('rootcoefa',0.0d0,1.0d0,rootcoefa)                 !
      endif                                                             ! NwRootExtr


! --- salt stress
      if (flsolute) then
        call rdsdor ('ecmax', 0.0d0,20.0d0,ecmax)
        call rdsdor ('ecslop',0.0d0,40.0d0,ecslop)
        call rdsdor ('c2eca', 0.0d0,1000.0d0,c2eca)
        call rdsdor ('c2ecb', 0.0d0,10.0d0,c2ecb)
        call rdsinr ('swc2ecf',1,2,swc2ecf)
        if (swc2ecf.eq.1) then
          call rdsdor ('c2ecf', 0.d0, 10.d0, c2ecf(1))
          if(numlay.ge.2) then
            do i = 2,numlay
              c2ecf(i) = c2ecf(1)
            enddo
          endif
        else if (swc2ecf.eq.2) then
          call rdfdor ('c2ecf', 0.d0, 10.d0, c2ecf,maho,numlay)
        endif
      endif

! --- interception
      call rdsinr ('swinter',0,2,swinter)
      if (swinter .eq. 1) then
        call rdsdor ('cofab',0.0d0,1.0d0,cofab)
      else if (swinter .eq. 2) then
        call rdador ('t',0.d0,366.d0,tinter,(magrs),ifnd)
        call rdfdor ('pfree',0.d0,1.d0,pfree,(magrs),ifnd)
        call rdfdor ('pstem',0.d0,1.d0,pstem,(magrs),ifnd)
        call rdfdor ('scanopy',0.d0,10.d0,scanopy,(magrs),ifnd)
        call rdfdor ('avprec',0.d0,100.d0,avprec,(magrs),ifnd)
        call rdfdor ('avevap',0.d0,10.d0,avevap,(magrs),ifnd)
        do i = 1, ifnd
          pfreetb(i*2) = pfree(i)
          pfreetb(i*2-1) = tinter(i)
          pstemtb(i*2) = pstem(i)
          pstemtb(i*2-1) = tinter(i)
          scanopytb(i*2) = scanopy(i)
          scanopytb(i*2-1) = tinter(i)
          avprectb(i*2) = avprec(i)
          avprectb(i*2-1) = tinter(i)
          avevaptb(i*2) = avevap(i)
          avevaptb(i*2-1) = tinter(i)
        end do
      endif

! --- read table with root distribution coefficients
      call rdador ('rdctb',0.0d0,100.0d0,rdctb,22,ifnd)
      
! --- determine whether irrigation scheduling is applied
      call rdsinr ('schedule',0,1,schedule)

! --- close file with crop data
      close (crp)


! --- CALCULATE NORMALIZED CUMULATIVE ROOT DENSITY FUNCTION
      if (swroottyp .eq. 1) then                                        !
! ---   root water extraction according to Feddes function              ! NwRootExtr
              

! ---   specify array ROOTDIS with root density distribution
        do i = 0,100
          depth = 0.01d0 * dble(i)
          rootdis(i*2+1) = depth
          rootdis(i*2+2) = afgen(rdctb,22,depth)
        enddo
        
! ---   calculate cumulative root density function
        do i = 1,202,2
! ---     relative depths
          cumdens(i) = rootdis(i)
        enddo
        sum = 0.d0
        cumdens(2) = 0.d0
        do i = 4,202,2
! ---     cumulative root density
          sum = sum + (rootdis(i-2)+rootdis(i)) * 0.5d0
     &               * (cumdens(i-1)-cumdens(i-3))
          cumdens(i) = sum
        enddo

! ---   normalize cumulative root density function to one
        do i = 2,202,2
          cumdens(i) = cumdens(i) / sum
        enddo
        endif                                                             ! NwRootExtr
      
      return
      end
! ----------------------------------------------------------------------
      subroutine wofost(task)
! ----------------------------------------------------------------------
!     date               : october 2004
!     purpose            : detailed crop growth routine
! ----------------------------------------------------------------------
      use variables
      implicit none
 
      integer   i1,ilvold,task,ilvoldpot,node

      real*8    lv(366),lvage(366),sla(366)
      real*8    lvpot(366),lvagepot(366),slapot(366)
      real*8    wlv,wrt,wst
      real*8    admi,afgen,amax,asrc,ccheck,cosld,cvf
      real*8    laicr,laiexp,laimax,lasum,mres,mrest,rdm
      real*8    dalv,dayl,delt,dmi,drlv,cptr0,ctr0,cptr1,ctr1
      real*8    drrt,drst,dslv,dslv1,dslv2,dslvt,dteff,dtga,dtsum,dvr
      real*8    dvred,dwlv,dwrt,dwst,fcheck,fl,fo,fr,fs
      real*8    fysdel,gass,gasst,gla,glaiex,glasol,grlv,grrt,grst
      real*8    gwrt,gwso,gwst,pgass,rest,rmres,rr
      real*8    sinld,slat,tadw,teff,twlv,twst
      real*8    wlvpot,lasumpot,wrtpot,wstpot
      real*8    dwrtpot,dwlvpot,dwstpot,dtgapot,tadwpot
      real*8    pgasspot,gasspot,rmrespot,mrespot,asrcpot,dmipot,rrpot
      real*8    admipot,grrtpot,drrtpot,gwrtpot,grlvpot
      real*8    dslvpot,restpot,dalvpot,drlvpot,gwsopot
      real*8    glasolpot,slatpot,glapot,grstpot,drstpot,gwstpot
      real*8    dslvtpot,twlvpot,twstpot
      real*8    dsinb,dsinbe,dso
      real*8    nihil,phead,dvspast1
      character tmp*11,messag*200
d     real*8    co2rootfix,co2rootloss,co2shootfix
d     character komma*1

      parameter (nihil=1.0d-7)
      parameter (delt=1.0d0)

      save
! ----------------------------------------------------------------------
      goto (1000,2000,3000) task

1000  continue

! === initialization ====================================================

! --- read crop data
      call readwofost (cropfil(icrop),pathcrop,swcf,cftb,idsl,dlo,dlc,
     &  tsumea,tsumam,dtsmtb,dvsend,tdwi,laiem,rgrlai,slatb,spa,
     &  ssa,span,tbase,kdif,kdir,eff,amaxtb,tmpftb,tmnftb,cvl,cvo,cvr,
     &  cvs,q10,rml,rmo,rmr,rms,rfsetb,frtb,fltb,fstb,fotb,perdl,rdrrtb,
     &  rdrstb,hlim1,hlim2u,hlim2l,hlim3h,hlim3l,hlim4,rsc,adcrh,adcrl,
     &  cofab,rdi,rri,rdc,rdctb,logf,schedule,cumdens,chtb,albedo,swetr,
     &  flsolute,ecmax,ecslop,c2eca,c2ecb,c2ecf,numlay,relni,cfet,
     &  alphacrit,swroottyp,wiltpoint,rootradius,rootcoefa,rsw)          ! NwRootExtr


d     open(unit=87,file='crop_CO2.csv',status='unknown')
d     write(87,*) 'date,co2rootfix,co2rootloss,co2shootfix'
d     open(unit=88,file='crop_dm.csv',status='unknown')
d     write(88,*) 'date,wrtpot,wrt,wstpot,wst,rdpot,rd,laipot,lai,
d    &gwrt,gwst,drrt,drlv,drst'

! --- maximum rooting depth & actual rooting depth
      rdm = min(rds,rdc)
      rd = min(rdi,rdm)
      rdpot = min(rdi,rdm)

! --- initial values of crop parameters
      swinter = 1
      dvs = 0.
      fr = afgen (frtb,30,dvs)
      fl = afgen (fltb,30,dvs)
      fs = afgen (fstb,30,dvs)
      fo = afgen (fotb,30,dvs)
      sla(1) = afgen (slatb,30,dvs)
      lvage(1) = 0.0d0
      ilvold = 1
      slapot(1) = afgen (slatb,30,dvs)
      lvagepot(1) = 0.0d0
      ilvoldpot = 1

! --- initial state variables of the crop
      wrt = fr*tdwi
      wrtpot = wrt
      tadw = (1.0d0-fr)*tdwi
      tadwpot = tadw
      wst = fs*tadw
      wstpot = wst
      wso = fo*tadw
      wsopot = wso
      wlv = fl*tadw
      wlvpot = wlv
!      laiem = wlv*sla(1)  is input !
      lv(1) = wlv
      lvpot(1) = wlv
      lasum = laiem     
      lasumpot = laiem     
      laiexp = laiem     
      glaiex = 0.0d0
      laimax = laiem
      lai = lasum+ssa*wst+spa*wso 
      laipot = lai 
      dwrt = 0.0d0
      dwrtpot = 0.0d0
      dwlv = 0.0d0
      dwlvpot = 0.0d0
      dwst = 0.0d0
      dwstpot = 0.0d0
      cf = afgen (cftb,(2*magrs),dvs)
      ch = afgen (chtb,(2*magrs),dvs)

! --- initial summation variables of the crop
      gasst = 0.0d0
      gasstpot = 0.0d0
      mrest = 0.0d0 
      mrestpot = 0.0d0 
      cptr0 = 0.0d0
      ctr0 = 0.0d0
      cptr1 = 0.0d0
      ctr1 = 0.0d0
      cwdm = 0.0d0
      cwdmpot = 0.0d0

! --- initialize matric flux potential                                  ! NwRootExtr
      phead = 0.d0   ! dummy                                            !
      if (swroottyp .eq. 2) then                                        !
        node = 1                                                        ! dummy node nr
        call MatricFlux(1,phead,node)                                   !
! ---   output of matric flux potential                                 !
!        if (swmfp.eq.1) call outmatricflux(2,mfp,numnod,tcum,          !
!     &   mflux,z,outfil,pathwork,project,ptra,h)                       !
      endif                                                             ! NwRootExtr

      return

2000  continue

! === calculate potential rate and state variables =====================

! --- rates of change of the crop variables ----------------------------

! --- phenological development rate  
      call astro(logf,swscre,daymeteo,lat,
     &           dayl,daylp,sinld,cosld,dsinb,dsinbe,dso)

! --- increase in temperature sum
      dtsum = afgen (dtsmtb,30,tav)

      if (dvs.lt.1.0d0) then     
! --- development during vegetative phase
        dvred = 1.0d0
        if (idsl.ge.1) 
     &          dvred = max(0.0d0,min(1.0d0,(daylp-dlc)/(dlo-dlc)))
        dvr = dvred*dtsum/tsumea
      else
! --- development during generative phase
        dvr = dtsum/tsumam
      endif    

! --- Correction for transition to values above 1.0 
!     (as suggested by PVWalsum, 20090821)
      if (dvs .lt. 1.0d0 .and. (dvs+dvr*delt) .gt. 1.0d0) then
        dvspast1 = (dvs+dvr*delt) - 1.0d0
        dvspast1 = dvspast1*(tsumea/dvred)/tsumam
        dvr      = ((1.0d0 + dvspast1) - dvs)/delt
      endif


! == = daily dry matter production 

! --- gross assimilation

      amax = afgen (amaxtb,30,dvs)
! --- correction for sub-optimum average daytemperature
      amax = amax * afgen (tmpftb,30,tavd)
      call totass (daynr,dayl,amax,eff,laipot,kdif,rad,sinld,cosld,
     &             dtgapot)
! --- correction for low minimum temperature
      dtgapot = dtgapot * afgen (tmnftb,30,tmnr)
! --- potential assimilation in kg ch2o per ha
      pgasspot = dtgapot * 30.0d0/44.0

! --- water stress reduction of pgass to gass (not for potential conditions)
      reltr = 1.0d0
      gasspot = pgasspot * reltr

! --- respiration and partitioning of carbohydrates between growth and
! --- maintenance respiration
      rmrespot = (rmr*wrtpot+rml*wlvpot+rms*wstpot+rmo*wsopot)*
     &            afgen(rfsetb,30,dvs)
      teff = q10**((tav-25.0d0)/10.0d0)
      mrespot = dmin1(gasspot,rmrespot*teff)
      asrcpot = gasspot - mrespot

! --- partitioning factors
      fr = afgen(frtb,30,dvs)
      fl = afgen(fltb,30,dvs)
      fs = afgen(fstb,30,dvs)
      fo = afgen(fotb,30,dvs)
! --- check on partitioning
      fcheck = fr+(fl+fs+fo)*(1.0d0-fr) - 1.0d0
      if (abs(fcheck).gt.0.0001d0) then
        write(tmp,'(f6.3)') dvs
        tmp = adjustl (tmp)
        Messag ='The sum of partitioning factors for leaves, stems'//
     &    ' and storage organs is not equal to one at development stage'
     &    //trim(tmp)//'.'
        call fatalerr ('cropd',messag)
      endif

! --- dry matter increase
      cvf = 1.0d0/((fl/cvl+fs/cvs+fo/cvo)*(1.0d0-fr)+fr/cvr)
      dmipot = cvf*asrcpot
! --- check on carbon balance
      ccheck = (gasspot-mrespot-(fr+(fl+fs+fo)*(1.0d0-fr))*dmipot/cvf)
     &         /max(0.0001d0,gasspot)      
      if (abs(ccheck).gt.0.0001d0) then
        Messag ='The carbon balance is not correct'
        call fatalerr ('cropd',messag)
      endif

! == = growth rate by plant organ

! --- root extension
      rrpot = min (rdm-rdpot,rri)
      if (fr.le.0.0d0.or.pgasspot.lt.1.0d0) rrpot = 0.0d0

! --- growth rate roots and aerial parts
      admipot = (1.0d0-fr)*dmipot
      grrtpot = fr*dmipot
      drrtpot = wrtpot*afgen (rdrrtb,30,dvs)
      gwrtpot = grrtpot - drrtpot

! --- weight of new leaves
      grlvpot = fl*admipot

! --- death of leaves due to water stress or high lai
      laicr = 3.2d0/kdif
      dslvpot = wlvpot*max(0.0d0,min(0.03d0,0.03d0*
     &           (laipot-laicr)/laicr))

! --- death of leaves due to exceeding life span:

! --- first: leaf death due to water stress or high lai is imposed 
! ---        on array until no more leaves have to die or all leaves
! ---        are gone

      restpot = dslvpot*delt
      i1 = ilvoldpot

      do while (restpot.gt.lvpot(i1).and.i1.ge.1)
        restpot = restpot - lvpot(i1) 
        i1 = i1-1
      enddo

! --- then: check if some of the remaining leaves are older than span,
! ---       sum their weights

      dalvpot = 0.0d0
      if (lvagepot(i1).gt.span .and. restpot.gt.0.0d0 .and.i1.ge.1) then
        dalvpot = lvpot(i1) - restpot
        restpot = 0.0d0
        i1 = i1-1
      endif

      do while (i1.ge.1.and.lvagepot(i1).gt.span)
        dalvpot = dalvpot+lvpot(i1)
        i1 = i1-1
      enddo

      dalvpot = dalvpot/delt

! --- finally: calculate total death rate leaves
      drlvpot = dslvpot + dalvpot

! --- physiologic ageing of leaves per time step
      fysdel = max (0.0d0,(tav-tbase)/(35.0d0-tbase))

! --- specific leaf area valid for current timestep
      slatpot = afgen (slatb,30,dvs)

! --- calculation of specific leaf area in case of exponential growth:
! --- leaf area not to exceed exponential growth curve
      if (laiexp.lt.6.0d0) then
        dteff = max (0.0d0,tav-tbase)
! ---   increase in leaf area during exponential growth
        glaiex = laiexp*rgrlai*dteff
! ---   source-limited increase in leaf area
        glasolpot = grlvpot*slatpot
! ---   actual increase is determined by lowest value
        glapot = min (glaiex,glasolpot)
! ---   slat will be modified in case gla equals glaiex
        if (grlvpot.gt.0.0d0) slatpot = glapot/grlvpot
      endif  

! --- growth rate stems
      grstpot = fs*admipot
! --- death rate stems
      drstpot = afgen (rdrstb,30,dvs)*wstpot
! --- net growth rate stems
      gwstpot = grstpot - drstpot

! --- growth rate storage organs
      gwsopot = fo*admipot

! ----integrals of the crop --------------------------------------------

! --- leaf death (due to water stress or high lai) is imposed on array 
! --- untill no more leaves have to die or all leaves are gone

      dslvtpot = dslvpot*delt
      i1 = ilvoldpot
      do while (dslvtpot.gt.0.and.i1.ge.1)
        if (dslvtpot.ge.lvpot(i1)) then
          dslvtpot = dslvtpot-lvpot(i1)
          lvpot(i1) = 0.0d0
          i1 = i1-1
        else
          lvpot(i1) = lvpot(i1)-dslvtpot
          dslvtpot = 0.0d0
        endif
      enddo

! --- leaves older than span die
      do while (lvagepot(i1).ge.span.and.i1.ge.1)
        lvpot(i1) = 0.0d0
        i1 = i1-1
      enddo

! --- oldest class with leaves
      ilvoldpot = i1

! --- shifting of contents, updating of physiological age
      do i1 = ilvoldpot,1,-1
        lvpot(i1+1) = lvpot(i1)
        slapot(i1+1) = slapot(i1)
        lvagepot(i1+1) = lvagepot(i1)+fysdel*delt
      enddo
      ilvoldpot = ilvoldpot + 1

! --- new leaves in class 1
      lvpot(1) = grlvpot*delt
      slapot(1) = slatpot
      lvagepot(1) = 0.0d0 

! --- calculation of new leaf area and weight
      lasumpot = 0.0d0
      wlvpot = 0.0d0
      do i1 = 1,ilvoldpot
        lasumpot = lasumpot + lvpot(i1)*slapot(i1)
        wlvpot = wlvpot + lvpot(i1)
      enddo

! --- leaf area index in case of exponential growth
      laiexp = laiexp+glaiex*delt

! --- dry weight of living plant organs
      wrtpot = wrtpot + gwrtpot*delt
      wstpot = wstpot + gwstpot*delt
      wsopot = wsopot + gwsopot*delt

! --- total above ground biomass
      tadwpot = wlvpot + wstpot + wsopot
      tadwpot = tadwpot ! for Forcheck

! --- dry weight of dead plant organs (roots,leaves & stems)
      dwrtpot = dwrtpot + drrtpot*delt
      dwlvpot = dwlvpot + drlvpot*delt
      dwstpot = dwstpot + drstpot*delt

! --- dry weight of dead and living plant organs
      twlvpot = wlvpot + dwlvpot
      twstpot = wstpot + dwstpot
      cwdmpot = twlvpot + twstpot + wsopot

! --- total gross assimilation and maintenance respiration
      gasstpot = gasspot + gasstpot
      mrestpot = mrespot + mrestpot

! --- leaf area index
      laipot = lasumpot + ssa*wstpot + spa*wsopot

! --- rooting depth
      rdpot = rdpot + rrpot

      return

3000  continue

! === calculate actual rate and state variables =====================

! --- rates of change of the crop variables ----------------------------

!     correction of potential transpiration in relation to reference crop
!     (default = 1.0, range = 0.8 - 1.2)
      ptra = cfet*ptra

! --- gross assimilation

      call totass (daynr,dayl,amax,eff,lai,kdif,rad,sinld,cosld,dtga)
! --- correction for low minimum temperature
      dtga = dtga * afgen (tmnftb,30,tmnr)
! --- potential assimilation in kg ch2o per ha
      pgass = dtga * 30.0d0/44.0

! --- water stress reduction of pgass to gass
      if(abs(ptra).lt.nihil) then
        reltr = 1.0d0
      else
        reltr = max(0.0d0,min(1.0d0,tra/ptra))
      endif
      gass = pgass * reltr
! --- management factor 
!     (nitrogen and other forms of stress, not accounted for)
      gass = gass * relni

! --- respiration and partitioning of carbohydrates between growth and
! --- maintenance respiration
      rmres = (rmr*wrt+rml*wlv+rms*wst+rmo*wso)*afgen(rfsetb,30,dvs)
      mres = dmin1(gass,rmres*teff)
      asrc = gass-mres

! --- dry matter increase
      dmi = cvf*asrc
! --- check on carbon balance
      ccheck = (gass-mres-(fr+(fl+fs+fo)*(1.0d0-fr))*dmi/cvf)
     &         /max(0.0001d0,gass)      
      if (abs(ccheck).gt.0.0001d0) then
        Messag ='The carbon balance is not correct'
        call fatalerr ('cropd',messag)
      endif

! --- growth rate by plant organ

! --- root extension
      rr = min (rdm-rd,rri)
      if (fr.le.0.0d0.or.pgass.lt.1.0d0) rr = 0.0d0

! --- growth rate roots and aerial parts
      admi = (1.0d0-fr)*dmi
      grrt = fr*dmi
      drrt = wrt*afgen (rdrrtb,30,dvs)
      gwrt = grrt-drrt

!       CO2 fixation and loss
d       co2rootfix = grrt*44.0d0/33.0d0
d       co2rootloss = drrt*44.0d0/33.0d0


! --- weight of new leaves
      grlv = fl*admi

! --- death of leaves due to water stress or high lai
      if(abs(ptra).lt.nihil) then
        dslv1 = 0.0d0
      else
        dslv1 = wlv*(1.0d0-tra/ptra)*perdl
      endif
      laicr = 3.2d0/kdif
      dslv2 = wlv*max(0.0d0,min(0.03d0,0.03d0*(lai-laicr)/laicr))
      dslv = max (dslv1,dslv2) 

! --- death of leaves due to exceeding life span:

! --- first: leaf death due to water stress or high lai is imposed on array
! ---        until no more leaves have to die or all leaves are gone

      rest = dslv*delt
      i1 = ilvold

      do while (rest.gt.lv(i1).and.i1.ge.1)
        rest = rest-lv(i1) 
        i1 = i1-1
      enddo

! --- then: check if some of the remaining leaves are older than span,
! ---       sum their weights

      dalv = 0.0d0
      if (lvage(i1).gt.span .and. rest.gt.0.0d0 .and.i1.ge.1) then
        dalv = lv(i1)-rest
        rest = 0.0d0
        i1 = i1-1
      endif

      do while (i1.ge.1.and.lvage(i1).gt.span)
        dalv = dalv+lv(i1)
        i1 = i1-1
      enddo

      dalv = dalv/delt

! --- finally: calculate total death rate leaves
      drlv = dslv+dalv

! --- specific leaf area valid for current timestep
      slat = afgen (slatb,30,dvs)

! --- calculation of specific leaf area in case of exponential growth:
! --- leaf area not to exceed exponential growth curve
      if (laiexp.lt.6.0d0) then
! ---   source-limited increase in leaf area
        glasol = grlv*slat
! ---   actual increase is determined by lowest value
        gla = min (glaiex,glasol)
! ---   slat will be modified in case gla equals glaiex
        if (grlv.gt.0.0d0) slat = gla/grlv
      endif  

! --- growth rate stems
      grst = fs*admi
! --- death rate stems
      drst = afgen (rdrstb,30,dvs)*wst
! --- net growth rate stems
      gwst = grst-drst

! --- growth rate storage organs
      gwso = fo*admi

! ----integrals of the crop --------------------------------------------

! --- phenological development stage
      dvs = dvs+dvr*delt

! --- leaf death (due to water stress or high lai) is imposed on array 
! --- untill no more leaves have to die or all leaves are gone

      dslvt = dslv*delt
      i1 = ilvold
      do while (dslvt.gt.0.and.i1.ge.1)
        if (dslvt.ge.lv(i1)) then
          dslvt = dslvt-lv(i1)
          lv(i1) = 0.0d0
          i1 = i1-1
        else
          lv(i1) = lv(i1)-dslvt
          dslvt = 0.0d0
        endif
      enddo

! --- leaves older than span die
      do while (lvage(i1).ge.span.and.i1.ge.1)
        lv(i1) = 0.0d0
        i1 = i1-1
      enddo

! --- oldest class with leaves
      ilvold = i1

! --- shifting of contents, updating of physiological age
      do i1 = ilvold,1,-1
        lv(i1+1) = lv(i1)
        sla(i1+1) = sla(i1)
        lvage(i1+1) = lvage(i1)+fysdel*delt
      enddo
      ilvold = ilvold+1

! --- new leaves in class 1
      lv(1) = grlv*delt
      sla(1) = slat
      lvage(1) = 0.0d0 

! --- calculation of new leaf area and weight
      lasum = 0.0d0
      wlv = 0.0d0
      do i1 = 1,ilvold
        lasum = lasum+lv(i1)*sla(i1)
        wlv = wlv+lv(i1)
      enddo

! --- dry weight of living plant organs
      wrt = wrt+gwrt*delt
      wst = wst+gwst*delt
      wso = wso+gwso*delt

! --- total above ground biomass
      tadw = wlv+wst+wso

! --- dry weight of dead plant organs (roots,leaves & stems)
      dwrt = dwrt+drrt*delt
      dwlv = dwlv+drlv*delt
      dwst = dwst+drst*delt

! --- dry weight of dead and living plant organs
!     twrt = wrt+dwrt
      twlv = wlv+dwlv
      twst = wst+dwst
      cwdm = twlv+twst+wso

! --- total gross assimilation and maintenance respiration
      gasst = gass + gasst
      mrest = mres + mrest

! --- leaf area index
      lai = lasum+ssa*wst+spa*wso
! --- determine maximum lai
      laimax = max (lai,laimax)

! --- rooting depth
      rd = rd+rr

! --- crop factor or crop height
      cf = afgen (cftb,(2*magrs),dvs)
      ch = afgen (chtb,(2*magrs),dvs)

!     CO2 fixation and root,shoot developm
d     co2shootfix = tagp * 44.0d0/33.0d0
d     komma = ","
d     write(87,'(a12,1x,3(a,f12.4))') 
d    &  date,komma,co2rootfix,komma,co2rootloss,komma,co2shootfix
d     write(88,'(a12,1x,20(a,f12.4))') 
d    &  date,komma,wrtpot,komma,wrt,komma,wstpot,komma,wst,
d    &  komma,rdpot,komma,rd,komma,laipot,komma,lai,
d    &  komma,gwrt,komma,gwst,komma,drrt,komma,drlv,komma,drst

! --- cumulative relative transpiration
      cptr0 = cptr0 + ptra  
      ctr0 = ctr0  + tra
      if (cptr0.le.nihil) then
        crt0 = 0.0d0
      else 
        crt0 = max(0.0d0,min(1.0d0,ctr0/cptr0))
      endif

      if (dvs.ge.1.0d0) then
        cptr1 = cptr1+ptra
        ctr1 = ctr1 + tra
        if (cptr1.le.nihil) then
          crt1 = 0.0d0
        else 
          crt1 = max(0.0d0,min(1.0d0,ctr1/cptr1))
        endif
      else
        crt1 = 1.0d0
      endif

! --- crop finish conditions based on dvs or lai
      if ( ((dvs.ge.dvsend) .or. (lai.le.0.002d0.and.dvs.gt.0.5d0)) 
     &     .and. (.not. flCropEnd) ) then
        flCropOutput =.true.
        flCropEnd = .true.
        icrop = icrop + 1
        flBareSoil = .true.
      endif

      return
      end

! ----------------------------------------------------------------------
      subroutine totass (daynr,dayl,amax,eff,lai,kdif,avrad,
     $sinld,cosld,dtga)
! ----------------------------------------------------------------------
! --- author: daniel van kraalingen, 1986
! --- calculates daily total gross assimilation (dtga) by performing
! --- a gaussian integration over time. at three different times of 
! --- the day, irradiance is computed and used to calculate the instan- 
! --- taneous canopy assimilation, whereafter integration takes place.
! --- more information on this routine is given by spitters et al./1988
! --- subroutines and functions called: assim, radiat
! ----------------------------------------------------------------------
      implicit none

      integer i,daynr

      real*8  lai,kdif
      real*8  amax,avrad,cosld,dayl,dtga,eff,fgros,gausr,hour,pardif
      real*8  pardir,sinb,sinld

      data    gausr /0.3872983d0/
! ----------------------------------------------------------------------
! --- three point gaussian integration over day
      dtga = 0.
      if (amax.lt.1.0d-10) return
      do 10 i=1,3
        hour = 12.0d0+dayl*0.5*(0.5d0+(i-2)*gausr)
! --- at a specified hour, diffuse and direct irradiance is computed
        call radiat (daynr,hour,dayl,sinld,cosld,avrad,sinb,pardir,
     $              pardif)
! --- irradiance and crop properties determine assimilation
        call assim (amax,eff,lai,kdif,sinb,pardir,pardif,fgros)
        if(i.eq.2) fgros=fgros*1.6
        dtga = dtga+fgros
10    continue
      dtga =dtga*dayl/3.6

      return
      end

! ----------------------------------------------------------------------
      subroutine radiat (daynr,hour,dayl,sinld,cosld,avrad,sinb,
     $                   pardir,pardif)
! ----------------------------------------------------------------------
! --- author: daniel van kraalingen, 1986
! --- calculates the fluxes of diffuse and direct photosynthetically
! --- active radiation from the total daily shortwave radiation actually
! --- received (avrad) for a given day of the year and hour of the day.
! --- the input variables dayl, sinld and cosld are calculated in astro.
! --- for more information: see spitters et al. (1988).
! ----------------------------------------------------------------------
      implicit none

      integer daynr

      real*8  aob,atmtr,avrad,cosld,dayl,dsinb,dsinbe,dso,frdif,hour
      real*8  par,pardif,pardir,pi,sc,sinb,sinld

      data    pi /3.1415926d0/
! ----------------------------------------------------------------------
! --- calculations on solar elevation
! --- sine of solar elevation sinb
      aob = sinld/cosld
      sinb = max (0.0d0,sinld+cosld*cos(2.0*pi*(hour+12.0d0)/24.0))
! --- integral of sinb
      dsinb = 3600.*(dayl*sinld+24.*cosld*sqrt(1.0d0-aob*aob)/pi)
! --- integral of sinb, corrected for lower atmospheric transmission
! --- at low solar elevations
      dsinbe = 3600.*(dayl*(sinld+0.4*(sinld*sinld+cosld*cosld*0.5))+
     $         12.0*cosld*(2.0d0+3.0*0.4*sinld)*sqrt(1.0d0-aob*aob)/pi)

! --- solar constant and daily extraterrestrial radiation
      sc = 1370.*(1.0d0+0.033*cos(2.0*pi*daynr/365.))
      dso = sc*dsinb

! --- diffuse light fraction from atmospheric transmission
      atmtr = avrad/dso
      if (atmtr.gt.0.75d0) frdif = 0.23d0
      if (atmtr.le.0.75d0.and.atmtr.gt.0.35d0) frdif = 1.33d0-1.46*atmtr
      if (atmtr.le.0.35d0.and.atmtr.gt.0.07d0) 
     $ frdif = 1.0d0-2.3*(atmtr-0.07d0)**2
      if (atmtr.le.0.07d0) frdif = 1.0d0

! --- photosynthetic active radiation, diffuse and direct
      par = 0.5*avrad*sinb*(1.0d0+0.4*sinb)/dsinbe
      pardif = min (par,sinb*frdif*atmtr*0.5*sc)
      pardir = par-pardif

      return
      end

! ----------------------------------------------------------------------
      subroutine assim (amax,eff,lai,kdif,sinb,pardir,pardif,fgros)
! ----------------------------------------------------------------------
!     author: daniel van kraalingen, 1986
!     calculates the gross co2 assimilation rate of the whole crop, 
!     fgros, by performing a gaussian integration over depth in the 
!     crop canopy. at three different depths in the canopy, i.e. for
!     different values of lai, the assimilation rate is computed for
!     given fluxes of photosynthetically active radiation, whereafter
!     integration over depth takes place. for more information: see 
!     spitters et al. (1988). the input variables sinb, pardir and 
!     pardif are calculated in radiat.
! ----------------------------------------------------------------------
      implicit none

      integer i

      real*8  lai,laic,kdif,kdirbl,kdirt
      real*8  amax,eff,fgl,fgros,fgrsh,fgrsun,fslla,gausr,pardif,pardir
      real*8  refh,refs,scv,sinb,visd,visdf,vispp,visshd,vist

      data    gausr /0.3872983d0/
! ----------------------------------------------------------------------
! --- extinction coefficients kdif,kdirbl,kdirt
      scv = 0.2d0
      refh = (1.0d0-sqrt(1.0d0-scv))/(1.0d0+sqrt(1.0d0-scv))
      refs = refh*2.0/(1.0d0+1.6*sinb)
      kdirbl = (0.5/sinb)*kdif/(0.8*sqrt(1.0d0-scv))
      kdirt = kdirbl*sqrt(1.0d0-scv)

! --- three point gaussian integration over lai
      fgros = 0.
      do 10 i = 1,3
        laic = 0.5*lai+gausr*(i-2)*lai
! --- absorbed diffuse radiation (vidf),light from direct
! --- origine (vist) and direct light(visd)
        visdf = (1.0d0-refs)*pardif*kdif  *exp(-kdif  *laic)
        vist = (1.0d0-refs)*pardir*kdirt *exp(-kdirt *laic)
        visd = (1.0d0-scv) *pardir*kdirbl*exp(-kdirbl*laic)
! --- absorbed flux in w/m2 for shaded leaves and assimilation
        visshd = visdf+vist-visd
        fgrsh = amax*(1.0d0-exp(-visshd*eff/amax))
! --- direct light absorbed by leaves perpendicular on direct
! --- beam and assimilation of sunlit leaf area
        vispp = (1.0d0-scv)*pardir/sinb
        if (vispp.le.0.0d0) fgrsun = fgrsh
        if (vispp.gt.0.0d0) fgrsun = amax*(1.0d0-
     $    (amax-fgrsh)*(1.0d0-exp(-vispp*eff/amax))/ (eff*vispp))
! --- fraction of sunlit leaf area (fslla) and local
! --- assimilation rate (fgl)
        fslla = exp(-kdirbl*laic)
        fgl = fslla*fgrsun+(1.0d0-fslla)*fgrsh
! --- integration
        if (i.eq.2) fgl = fgl*1.6
        fgros = fgros+fgl
10    continue
      fgros = fgros*lai/3.6

      return
      end


! ----------------------------------------------------------------------
      subroutine grass(task)
! ----------------------------------------------------------------------
!     Date               : November 2004
!     Purpose            : detailed grass growth routine 
! ----------------------------------------------------------------------
      use variables
      implicit none
 
      integer   i1,ilvold,task,ilvoldpot,idregr,idregrpot
      integer   idelaypot,idelay,node

      real*8    laicr,laiexp,laimax,lasum,mres,rid
!     real*8    rdm,rr,rrpot
      real*8    lv(366),lvage(366),sla(366)
      real*8    admi,afgen,amax,asrc,ccheck,cosld,cvf,rlwtb(22)
      real*8    dalv,dayl,delt,dmi,drlv,cptr0,ctr0
      real*8    drrt,drst,dslv,dslv1,dslv2,dslvt,dteff,dtga
      real*8    dwlv,dwrt,dwst,fcheck,fl,fr,fs,laiexppot
      real*8    fysdel,gass,gla,glaiex,glasol,grlv,grrt,grst
      real*8    gwrt,gwst,pgass,rest,rmres
      real*8    sinld,slat,teff,twlv,twst,wlv,wrt,wst,glaiexpot
      real*8    lvpot(366),lvagepot(366),slapot(366)
      real*8    wlvpot,lasumpot,wrtpot,wstpot,drst1,drst2
      real*8    dwrtpot,dwlvpot,dwstpot,dtgapot,drst1pot,drst2pot
      real*8    pgasspot,gasspot,rmrespot,mrespot,asrcpot,dmipot
      real*8    admipot,grrtpot,drrtpot,gwrtpot,grlvpot,dslv1pot
      real*8    dslv2pot,dslvpot,restpot,dalvpot,drlvpot
      real*8    glasolpot,slatpot,glapot,grstpot,drstpot,gwstpot
      real*8    dslvtpot,twlvpot,twstpot,tagpspot,tagps
      real*8    dsinb,dsinbe,dso
      real*8    nihil,phead
      real*8    dmharvest1,dmharvest2,dmlastharvest,dateharvest(999)
      real*8    wrtmax,rdm
      real*8    grazingfactor
      real*8    nsuptab(magrs),dmfac(magrs),relnitab(2*magrs),nsupply        ! Nwgrassland
      integer   daylastharvest,iharvest,iharvestpot,swharvest,swgrazing
      logical   flharvestpot,flharvest,flgrazing,flgrazingpot
      character tmp*11,messag*200
d     real*8    co2rootfix,co2rootloss,co2shootfix
d     character komma*1

      parameter (nihil=1.0d-7)
      parameter (delt=1.0d0)

      save


! ----------------------------------------------------------------------
      goto (1000,2000,3000) task

1000  continue

! === initialization ====================================================

! --- read grass input data
      call readgrass (cropfil(icrop),pathcrop,tdwi,laiem,rgrlai,slatb,
     &  ssa,span,tbase,kdif,kdir,eff,amaxtb,tmpftb,tmnftb,cvl,cvr,cvs,
     &  q10,rml,rmr,rms,rfsetb,frtb,fltb,fstb,perdl,rdrrtb,
     &  rdrstb,hlim1,hlim2u,hlim2l,hlim3h,hlim3l,hlim4,rsc,adcrh,adcrl,
     &  cofab,rdi,rri,rdc,rdctb,rlwtb,logf,schedule,cumdens,
     &  flsolute,ecmax,ecslop,c2eca,c2ecb,c2ecf,numlay,dateharvest,
     &  swharvest,dmharvest1,dmharvest2,swgrazing,grazingfactor,
     &  daylastharvest,dmlastharvest,wrtmax,
     &  nsuptab,dmfac,relnitab,nsupply,
     &  swcf,swetr,cftb,chtb,cfet,alphacrit,
     &  swroottyp,wiltpoint,rootradius,rootcoefa,rsw)      

         
! --- initial values
      swinter = 1
      iharvest = 1
      iharvestpot = 1
      flharvest = .false.
      flharvestpot = .false.
      
      flgrazing = .false.
      flgrazingpot = .false.

d     open(unit=87,file='grass_CO2.csv',status='unknown')
d     write(87,*) 'date,co2rootfix,co2rootloss,co2shootfix'
d     open(unit=88,file='grass_dm.csv',status='unknown')
d     write(88,*) 'date,wrtpot,wrt,wstpot,wst,rdpot,rd,laipot,lai,
d    &gwrt,gwst,drrt,drlv,drst'

! --- development stage (not used by Grassland, instead Daynrs are used)
      dvs = -99.99d0

! --- maximum rooting depth & actual rooting depth
      rdm = min(rds,rdc)
      rd = min(rdi,rdm)
      rdpot = min(rdi,rdm)

! --- initial values of crop parameters
      rid = 1.0d0
      fr = afgen (frtb,30,rid)
      fl = afgen (fltb,30,rid)
      fs = afgen (fstb,30,rid)
      sla(1) = afgen (slatb,30,rid)
      lvage(1) = 0.d0
      ilvold = 1
      idregr = 0
      slapot(1) = afgen (slatb,30,rid)
      lvagepot(1) = 0.d0
      ilvoldpot = 1
      idregrpot = 0

! --- initial state variables of the crop
      wrt = fr*tdwi
      wrtpot = wrt
      wst = fs*(1.0d0-fr)*tdwi
      wstpot = wst
      wlv = laiem/sla(1)
      wlvpot = wlv
      lv(1) = wlv
      lvpot(1) = lv(1)
      lasum = laiem
      lasumpot = lasum     
      glaiex = 0.0d0
      laiexp = laiem
      laiexppot = laiem
      laimax = laiem
      lai = lasum+ssa*wst
      laipot = lai
      dwrt = 0.d0
      dwrtpot = dwrt
      dwlv = 0.d0
      dwlvpot = dwlv
      dwst = 0.d0
      dwstpot = dwst
      rid = dble(daycrop)
      cf = afgen (cftb,(2*magrs),rid)
      ch = afgen (chtb,(2*magrs),rid)

! --- initial summation variables of the crop
      tagp = wlv+wst
      tagppot = tagp
      cptr0 = 0.0d0
      ctr0 = 0.0d0
      tagpt = 0.0d0
      tagptpot = tagpt

! --- initialize matric flux potential                                  ! NwRootExtr
      phead = 0.d0   ! dummy                                            !
      if (swroottyp .eq. 2) then                                        !
        node = 1                                                        ! dummy node nr
        call MatricFlux(1,phead,node)                                   !
! ---   output of matric flux potential                                 !
!        if (swmfp.eq.1) call outmatricflux(2,mfp,numnod,tcum,          !
!     &   mflux,z,outfil,pathwork,project,ptra,h)                       !
      endif                                                             ! NwRootExtr

      return

2000  continue

! === calculate potential rate and state variables ======================================

! --- rates of change of the grass variables ---------------------------------------------

      rid = dble(daycrop)
      cf = afgen (cftb,(2*magrs),rid)
      ch = afgen (chtb,(2*magrs),rid)

! --- skip in case of regrowth
!      if (daycrop.ne.0.and.daycrop.lt.idregrpot) goto 2100
      if (daycrop.eq.0 .or.daycrop.ge.idregrpot) then

! ===   daily dry matter production ===

! ---   gross assimilation
        amax = afgen (amaxtb,30,rid)
! ---   correction for sub-optimum average daytemperature
        amax = amax * afgen (tmpftb,30,tavd)
        call astro(logf,swscre,daymeteo,lat,
     &           dayl,daylp,sinld,cosld,dsinb,dsinbe,dso)
        call totass (daynr,dayl,amax,eff,laipot,kdif,rad,sinld,cosld,
     &             dtgapot)
! ---   correction for low minimum temperature
        dtgapot = dtgapot * afgen (tmnftb,30,tmnr)
! ---   potential assimilation in kg ch2o per ha
        pgasspot = dtgapot * 30.0d0/44.0d0

! --- water stress reduction of pgass to gass (not for potential conditions)
        reltr = 1.0d0
        gasspot = pgasspot * reltr

! ---   respiration and partitioning of carbohydrates between growth and
! ---   maintenance respiration
        rmrespot=(rmr*wrtpot+rml*wlvpot+rms*wstpot)*afgen(rfsetb,30,rid)
        teff = q10**((tav-25.0d0)/10.0d0)
        mrespot = min (gasspot,rmrespot*teff)
        asrcpot = gasspot-mrespot

! ---   partitioning factors
        fr = afgen(frtb,30,rid)
        fl = afgen(fltb,30,rid)
        fs = afgen(fstb,30,rid)
! ---   check on partitioning
        fcheck = fr+(fl+fs)*(1.0d0-fr) - 1.0d0
        if (abs(fcheck).gt.0.0001d0) then
          write(tmp,'(f6.3)') rid
          tmp = adjustl (tmp)
          Messag ='The sum of partitioning factors for leaves, stems'//
     &    ' and storage organs is not equal to one at time '
     &    //trim(tmp)//'.'
          call fatalerr ('cropd',messag)
        endif

! ---   dry matter increase
        cvf = 1.0d0/((fl/cvl+fs/cvs)*(1.0d0-fr)+fr/cvr)
        dmipot = cvf*asrcpot
! ---   check on carbon balance
        ccheck = (gasspot-mrespot-(fr+(fl+fs)*(1.0d0-fr))*dmipot/cvf)
     &         /max(0.0001d0,gasspot)      
        if (abs(ccheck).gt.0.0001d0) then
          Messag ='The carbon balance is not correct'
          call fatalerr ('cropd',messag)
        endif


! ===   growth rate by plant organ ===

! ---   root length (not used, because Rooting depth is (for grassland) 
!       dependent on available root biomass (weight)
!        rrpot = min (rdm-rdpot,rri)
!        if (fr.le.0.or.pgasspot.lt.1.0d0) rrpot = 0.0d0

! ---   growth rate roots and aerial parts
! ---   after reaching a live weight of wrtmax (default 2500 kg), the
! ---   growth of the roots is balanced by the death of root tissue
        grrtpot = fr*dmipot
        if (wrtpot.gt.wrtmax) then
          drrtpot = grrtpot
        else
          drrtpot = wrtpot*afgen (rdrrtb,30,rid)
        endif
        gwrtpot = grrtpot-drrtpot

! ---   growth rate leaves

! ---   weight of new leaves
        admipot = (1.0d0-fr)*dmipot
        grlvpot = fl*admipot

! ---   death of leaves due to water stress or high lai
        dslv1pot = 0.0d0
        laicr = 3.2d0/kdif
        dslv2pot=wlvpot*max(0.0d0,
     &                  min(0.03d0,0.03d0*(laipot-laicr)/laicr))
        dslvpot = max (dslv1pot,dslv2pot) 

! ---   death of leaves due to exceeding life span;
! ---   leaf death is imposed on array until no more leaves have
! ---   to die or all leaves are gone

        restpot = dslvpot*delt
        i1 = ilvoldpot

        do while (restpot.gt.lvpot(i1).and.i1.ge.1)
          restpot = restpot-lvpot(i1) 
          i1 = i1-1
        enddo

! ---   check if some of the remaining leaves are older than span,
! ---   sum their weights

        dalvpot = 0.0d0
        if (lvagepot(i1).gt.span.and.restpot.gt.0.and.i1.ge.1) then
          dalvpot = lvpot(i1)-restpot
          restpot = 0.0d0
          i1 = i1-1
        endif

        do while (i1.ge.1.and.lvagepot(i1).gt.span)
          dalvpot = dalvpot+lvpot(i1)
          i1 = i1-1
        enddo

        dalvpot = dalvpot/delt

! ---   death rate leaves and growth rate living leaves
        drlvpot = dslvpot+dalvpot

! ---   physiologic ageing of leaves per time step
        fysdel = max (0.0d0,(tav-tbase)/(35.0d0-tbase))

! ---   leaf area not to exceed exponential growth curve
        slatpot = afgen (slatb,30,rid)
        if (laiexppot.lt.6.0d0) then
          dteff = max (0.0d0,tav-tbase)
          glaiexpot = laiexppot*rgrlai*dteff
! ---   source-limited increase in leaf area
          glasolpot = grlvpot*slatpot
          glapot = min (glaiexpot,glasolpot)
! ---   adjustment of specific leaf area of youngest leaf class
          if (grlvpot.gt.0.0d0) slatpot = glapot/grlvpot
        endif  

! ---   growth rate stems
        grstpot = fs*admipot
! ---   death of stems due to water stress is zero in case of potential growth
        drst1pot = 0.0d0
! ---   death of stems due to ageing
        drst2pot = afgen (rdrstb,30,rid)*wstpot
        drstpot = (drst1pot+drst2pot)/delt 
        gwstpot = grstpot-drstpot

! ----  integrals of the crop --------------------------------------------
!       after cutting, growth is initialized again and the weight of the sward is stored
! ---   harvest criteria (open to personal choice of user) 
!  INPUT    dmharvest, daylastharvest, dmlastharvest
!      if (tagppot.gt.4200.0d0 .or. 
!     &             (daycrop.gt.210 .and. tagppot.gt.3700.0d0)) then
        if (swharvest.eq.1) then
          if(tagppot.gt.dmharvest1 .or. (daycrop.gt.daylastharvest .and. 
     &                                 tagppot.gt.dmlastharvest).or. 
     &                                 flgrazingpot .eqv. .true.) then
            if (swgrazing.eq.1) then        ! grazing
                flgrazingpot = .true.
                      
                grlvpot = 0.0d0                           ! growing rate leaves adaption
                gwstpot = -1.0d0*grazingfactor* wstpot    ! growing rate stems adaption 
                 i1 = ilvoldpot
                 do while (i1.ge.1)
                    lvpot(i1) = (1.0d0-grazingfactor)*lvpot(i1)
                    i1 = i1 -1
                 end do
                 if(tagppot.lt.800.0d00) then            ! assumed last grazing day
                    flgrazingpot = .false.
                    flharvestpot = .true.
                 endif 
            else if (swgrazing.eq.2) then         ! growth continues followed by mowing
                if (tagppot.gt.dmharvest2.or. 
     &                             (daycrop.gt.daylastharvest
     &                              .and.tagppot.gt.dmlastharvest)) then
                    flharvestpot = .true.
                else
                    flharvestpot = .false.
                endif
            end if
          else
            flharvestpot = .false.
          endif
        endif
        if (swharvest.eq.2) then                         ! mowing using mowing dates
          if(t1900.gt.dateharvest(iharvestpot)) then
            iharvestpot = iharvestpot + 1
            flharvestpot = .true.
          else
            flharvestpot = .false.
          endif
        endif
        if (flharvestpot) then
          flharvestpot = .false.
          lasumpot = laiem
          slapot(1) = afgen (slatb,30,rid)
          wlvpot = lasumpot/slapot(1)
          fl = afgen (fltb,30,rid)
          fs = afgen (fstb,30,rid)
          wstpot = fs/fl*wlvpot
          dwlvpot = 0.0d0
          dwstpot = 0.0d0
          lvagepot(1) = 0.0d0
          ilvoldpot = 1
          laiexppot = laiem
          lvpot(1) = wlvpot

          gwstpot = 0.0d0
          gwrtpot = 0.0d0
          drlvpot = 0.0d0
          drstpot = 0.0d0
          drrtpot = 0.0d0

          tagpspot =max(0.0d0,(tagppot-(wlvpot+dwlvpot+wstpot+dwstpot)))
          tagptpot = tagptpot + tagpspot

! ---     regrowth delay after handbook p.r.

          if (tagpspot.lt.2000.0d0) idelaypot=1
          if (tagpspot.ge.2000.0d0.and.tagpspot.lt.2500.0d0) idelaypot=2
          if (tagpspot.ge.2500.0d0.and.tagpspot.lt.3000.0d0) idelaypot=3
          if (tagpspot.ge.3000.0d0.and.tagpspot.lt.3500.0d0) idelaypot=4
          if (tagpspot.ge.3500.0d0.and.tagpspot.lt.4000.0d0) idelaypot=5
          if (tagpspot.ge.4000.0d0) idelaypot = 6

          idregrpot = daycrop + idelaypot + 3

        endif

        if (daycrop.ge.idregrpot) then

! ---     leaf death is imposed on array untill no more leaves have to die or all leaves are gone

          dslvtpot = dslvpot*delt
          i1 = ilvoldpot
          do while (dslvtpot.gt.0.and.i1.ge.1)
            if (dslvtpot.ge.lvpot(i1)) then
              dslvtpot = dslvtpot-lvpot(i1)
              lvpot(i1) = 0.0d0
              i1 = i1-1
            else
              lvpot(i1) = lvpot(i1)-dslvtpot
              dslvtpot = 0.0d0
            endif
          enddo

          do while (lvagepot(i1).ge.span.and.i1.ge.1)
            lvpot(i1) = 0.0d0
            i1 = i1-1
          enddo

          ilvoldpot = i1

! ---     shifting of contents, integration of physiological age
          do i1 = ilvoldpot,1,-1
            lvpot(i1+1) = lvpot(i1)
            slapot(i1+1) = slapot(i1)
            lvagepot(i1+1) = lvagepot(i1)+fysdel*delt
          enddo
          ilvoldpot = ilvoldpot+1

! ---     new leaves in class 1
          lvpot(1) = grlvpot*delt
          slapot(1) = slatpot
          lvagepot(1) = 0.d0 

! ---     calculation of new leaf area and weight
          lasumpot = 0.d0
          wlvpot = 0.d0
          do i1 = 1,ilvoldpot
            lasumpot = lasumpot+lvpot(i1)*slapot(i1)
            wlvpot = wlvpot+lvpot(i1)
          enddo

          laiexppot = laiexppot+glaiexpot*delt

        endif
      endif

! --- dry weight of living plant organs
      wrtpot = wrtpot+gwrtpot*delt
      wstpot = wstpot+gwstpot*delt

! --- dry weight of dead plant organs (roots,leaves & stems)
      dwrtpot = dwrtpot+drrtpot*delt
      dwlvpot = dwlvpot+drlvpot*delt
      dwstpot = dwstpot+drstpot*delt

! --- dry weight of dead and living plant organs
      twlvpot = wlvpot+dwlvpot
      twstpot = wstpot+dwstpot
      tagppot = twlvpot+twstpot

! --- leaf area index
      laipot = lasumpot+ssa*wstpot

! --- rooting depth as function of available root weight (not root rate!)
!      rdpot = rdpot+rrpot
      rdpot = afgen (rlwtb,22,wrtpot)
      rdpot = min(rdpot,rdm)

      return

3000  continue

! === calculate actual rate and state variables ======================================

! --- rates of change of the crop variables ---------------------------------------------

!     correction of potential transpiration in relation to reference crop
!     (default = 1.0, range = 0.8 - 1.2)
      ptra = cfet*ptra

! --- skip in case of regrowth
      if (daycrop.eq.0.or.daycrop.ge.idregr) then

! ===   daily dry matter production ===

! ---   gross assimilation
        amax = afgen (amaxtb,30,rid)
! ---   correction for sub-optimum average daytemperature
        amax = amax * afgen (tmpftb,30,tavd)
! ---   gross assimilation
        call astro(logf,swscre,daymeteo,lat, 
     &           dayl,daylp,sinld,cosld,dsinb,dsinbe,dso)
        call totass (daynr,dayl,amax,eff,lai,kdif,rad,sinld,cosld,dtga)
! ---   correction for low minimum temperature
        dtga = dtga * afgen (tmnftb,30,tmnr)
! ---   potential assimilation in kg ch2o per ha
        pgass = dtga * 30.0d0/44.0d0

! ---   water stress reduction of pgass to gass
        if(abs(ptra).lt.nihil) then
          reltr = 1.0d0
        else
          reltr = max(0.0d0,min(1.0d0,tra/ptra))
        endif
        gass = pgass * reltr
! ---   nitrogen stress reduction of pgass to gass
        relni = afgen(relnitab,magrs,nsupply)
        gass = gass * relni

! ---   respiration and partitioning of carbohydrates between growth and
! ---   maintenance respiration
        rmres = (rmr*wrt+rml*wlv+rms*wst)*afgen(rfsetb,30,rid)
        teff = q10**((tav-25.0d0)/10.0d0)
        mres = min (gass,rmres*teff)
        asrc = gass-mres

! ---   dry matter increase
        cvf = 1.0d0/((fl/cvl+fs/cvs)*(1.0d0-fr)+fr/cvr)
        dmi = cvf*asrc
! ---   check on carbon balance
        ccheck = (gass-mres-(fr+(fl+fs)*(1.0d0-fr))*dmi/cvf)
     &         /max(0.0001d0,gass)      
        if (abs(ccheck).gt.0.0001d0) then
          Messag ='The carbon balance is not correct'
          call fatalerr ('cropd',messag)
        endif

! ===   growth rate by plant organ ===

! ---   root length
!        rr = min (rdm-rd,rri)
!        if (fr.le.0.or.pgass.lt.1.0d0) rr = 0.0d0
!        rr = 0.0d0

! ---   growth rate roots and aerial parts
! ---   after reaching a live weight of wrtmax (default 2500 kg), the
! ---   growth of the roots is balanced by the death of root tissue
        grrt = fr*dmi
        if (wrt.gt.wrtmax) then
!original drrt = wrt - wrtmax
          drrt = grrt
!         CO2 loss
d         co2rootloss = drrt*44.0d0/33.0d0
!original  grrt = 0.0d0
        else
          drrt = wrt*afgen (rdrrtb,30,rid)
        endif
        admi = (1.0d0-fr)*dmi
!original drrt = 0.0d0
        gwrt = grrt-drrt

!       CO2 fixation and loss
d       co2rootfix = grrt*44.0d0/33.0d0
d       co2rootloss = co2rootloss + drrt*44.0d0/33.0d0

! ---   growth rate leaves

! ---   weight of new leaves
        grlv = fl*admi

! ---   death of leaves due to water stress or high lai
        if(abs(ptra).lt.nihil) then
          dslv1 = 0.0d0
        else
          dslv1 = wlv*(1.0d0-tra/ptra)*perdl
        endif
        laicr = 3.2d0/kdif
        dslv2 = wlv*max(0.0d0,min(0.03d0,0.03d0*(lai-laicr)/laicr))
        dslv = max (dslv1,dslv2) 

! ---   death of leaves due to exceeding life span;
! ---   leaf death is imposed on array until no more leaves have
! ---   to die or all leaves are gone

        rest = dslv*delt
        i1 = ilvold

        do while (rest.gt.lv(i1).and.i1.ge.1)
          rest = rest-lv(i1) 
          i1 = i1-1
        enddo

! ---   check if some of the remaining leaves are older than span,
! ---   sum their weights

        dalv = 0.0d0
        if (lvage(i1).gt.span.and.rest.gt.0.and.i1.ge.1) then
          dalv = lv(i1)-rest
          rest = 0.0d0
          i1 = i1-1
        endif

        do while (i1.ge.1.and.lvage(i1).gt.span)
          dalv = dalv+lv(i1)
          i1 = i1-1
        enddo

        dalv = dalv/delt

! ---   death rate leaves and growth rate living leaves
        drlv = dslv+dalv

! ---   physiologic ageing of leaves per time step
        slat = afgen (slatb,30,rid)

! ---   leaf area not to exceed exponential growth curve
        if (laiexp.lt.6.0d0) then
          dteff = max (0.0d0,tav-tbase)
          glaiex = laiexp*rgrlai*dteff
! ---     source-limited increase in leaf area
          glasol = grlv*slat
          gla = min (glaiex,glasol)
! ---     adjustment of specific leaf area of youngest leaf class
          if (grlv.gt.0.0d0) slat = gla/grlv
        endif  

! ---   growth rate stems
        grst = fs*admi
! ---   death of stems due to water stress
        if(abs(ptra).lt.nihil) then
          drst1 = 0.0d0
        else
          drst1 = wst*(1.0d0-tra/ptra)*perdl
        endif
! ---   death of stems due to ageing
        drst2 = afgen (rdrstb,30,rid)*wst
        drst = (drst1+drst2)/delt 
        gwst = grst-drst

! ----  integrals of the crop --------------------------------------------

!       after cutting, growth is initialized again and the weight of the sward is stored

        if (swharvest.eq.1) then
          if(tagp.gt.dmharvest1 .or. (daycrop.gt.daylastharvest .and. 
     &                              tagp.gt.dmlastharvest).or. 
     &                              flgrazing .eqv. .true.) then
            if (swgrazing.eq.1) then               ! grazing
                flgrazing = .true.
                grlv = 0.0d0                       ! growing rate leaves adaption
                gwst = -1.0d0*grazingfactor* wst   ! growing rate stems adaption
                 i1 = ilvold
                 do while (i1.ge.1)
                    lv(i1) = (1.0d0-grazingfactor)*lv(i1) ! reduction of leaf weight
                    i1 = i1 -1
                 end do
                 if(tagp.lt.800.0d00) then  ! assumed last grazing day
                    flgrazing = .false.
                    flharvest = .true.
                 endif 
            else if (swgrazing.eq.2) then   ! growth continues followed by harvest
                if (tagp.gt.dmharvest2.or. 
     &                         (daycrop.gt.daylastharvest .and. 
     &                                   tagp.gt.dmlastharvest)) then
                    flharvest = .true.
                else
                    flharvest = .false.
                endif
            end if
          else
            flharvest = .false.
          endif
        endif
        if (swharvest.eq.2) then            ! mowing using mowing dates
          if(t1900.gt.dateharvest(iharvest)) then
            iharvest = iharvest + 1
            flharvest = .true.
          else
            flharvest = .false.
          endif
        endif
        if (flharvest) then
          flharvest = .false.
          lasum = laiem
          sla(1) = afgen (slatb,30,rid)
          wlv = lasum/sla(1)
          wst = fs/fl*wlv
          dwlv = 0.0d0
          dwst = 0.0d0
          lvage(1) = 0.0d0
          ilvold = 1
          laiexp = laiem
          lv(1) = wlv

          gwst = 0.0d0
          gwrt = 0.0d0
          drlv = 0.0d0
          drst = 0.0d0
          drrt = 0.0d0

          tagps = max (0.0d0,(tagp-(wlv+dwlv+wst+dwst)))
          tagpt = tagpt + tagps

! ---     regrowth delay after handbook p.r.

          if (tagps.lt.2000.0d0) idelay = 1           
          if (tagps.ge.2000.0d0.and.tagps.lt.2500.0d0) idelay = 2
          if (tagps.ge.2500.0d0.and.tagps.lt.3000.0d0) idelay = 3
          if (tagps.ge.3000.0d0.and.tagps.lt.3500.0d0) idelay = 4
          if (tagps.ge.3500.0d0.and.tagps.lt.4000.0d0) idelay = 5
          if (tagps.ge.4000.0d0) idelay = 6

          idregr = daycrop + idelay + 3

        endif

        if (daycrop.ge.idregr) then

! ---     physiologic ageing of leaves per time step
          fysdel = max (0.0d0,(tav-tbase)/(35.0d0-tbase))

! ---     leaf death is imposed on array untill no more leaves have to die or all leaves are gone

          dslvt = dslv*delt
          i1 = ilvold
          do while (dslvt.gt.0.and.i1.ge.1)
            if (dslvt.ge.lv(i1)) then
              dslvt = dslvt-lv(i1)
              lv(i1) = 0.0d0
              i1 = i1-1
            else
              lv(i1) = lv(i1)-dslvt
              dslvt = 0.0d0
            endif
          enddo

          do while (lvage(i1).ge.span.and.i1.ge.1)
            lv(i1) = 0.0d0
            i1 = i1-1
          enddo

          ilvold = i1

! ---     shifting of contents, integration of physiological age
          do i1 = ilvold,1,-1
            lv(i1+1) = lv(i1)
            sla(i1+1) = sla(i1)
            lvage(i1+1) = lvage(i1)+fysdel*delt
          enddo
          ilvold = ilvold+1

! ---     new leaves in class 1
          lv(1) = grlv*delt
          sla(1) = slat
          lvage(1) = 0.d0 

! ---     calculation of new leaf area and weight
          lasum = 0.d0
          wlv = 0.d0
          do i1 = 1,ilvold
            lasum = lasum+lv(i1)*sla(i1)
            wlv = wlv+lv(i1)
          enddo

          laiexp = laiexp+glaiex*delt

        endif
      endif

! --- dry weight of living plant organs
      wrt = wrt+gwrt*delt
      wst = wst+gwst*delt

! --- dry weight of dead plant organs (roots,leaves & stems)
      dwrt = dwrt+drrt*delt
      dwlv = dwlv+drlv*delt
      dwst = dwst+drst*delt

! --- dry weight of dead and living plant organs
!     twrt = wrt+dwrt
      twlv = wlv+dwlv
      twst = wst+dwst
      tagp = twlv+twst

! --- leaf area index
      lai = lasum+ssa*wst
      laimax = max (lai,laimax)

! --- rooting depth as function of root weight
!      rd = rd+rr
      rd = afgen (rlwtb,22,wrt)
      rd= min(rd,rdm)

!     CO2 fixation and root,shoot developm
d     co2shootfix = tagp * 44.0d0/33.0d0
d     komma = ","
d     write(87,'(a12,1x,3(a,f12.4))') 
d    &  date,komma,co2rootfix,komma,co2rootloss,komma,co2shootfix
d     write(88,'(a12,1x,20(a,f12.4))') 
d    &  date,komma,wrtpot,komma,wrt,komma,wstpot,komma,wst,
d    &  komma,rdpot,komma,rd,komma,laipot,komma,lai,
d    &  komma,gwrt,komma,gwst,komma,drrt,komma,drlv,komma,drst

! --- cumulative relative transpiration
      cptr0 = cptr0 + ptra  
      ctr0 = ctr0  + tra
      if (cptr0.le.nihil) then
        crt0 = 0.0d0
      else
        crt0 = max(0.0d0,min(1.0d0,ctr0/cptr0))
      endif

      return
      end

! NwRootExtr New Subroutine for new Rootextraction ! NwRootExtr 

! ----------------------------------------------------------------------
      subroutine MatricFlux(task,phead,node) 
! ----------------------------------------------------------------------
!     Date               : October 2006   
!     Purpose            : Initialize and calculate matric flux potential
! ----------------------------------------------------------------------

      use variables
      implicit none

      integer task,lay,count,start,node
      real*8  phead1,phead2,wcontent,conduc1,conduc2,watcon,hconduc
      real*8  phead,logphead

      goto (1000, 2000) task

1000  continue

! === initialization =========================================================

      do lay = 1,numlay
        do count = 1,600
          mfluxtable(lay,count) = 0.0d0
        enddo
      enddo

      start = int(100.d0*log10(-wiltpoint))
      do lay = 1,numlay
        phead1 = -10.d0**(dble(start)/100.d0)

!       find first Node of the Layer
        Node = nod1lay(lay)

        wcontent = watcon(Node,phead1,cofgen(1,Node),
     &                    swsophy,numtab,sptab)
        conduc1 = hconduc (Node,wcontent,cofgen(1,Node),swfrost,10.d0,
     &                     swsophy,numtab,sptab)

        do count = start-1,1,-1
          phead2 = -10.d0**(dble(count)/100.d0)
          wcontent = watcon(Node,phead2,cofgen(1,Node),
     &                      swsophy,numtab,sptab)
          conduc2 = hconduc (Node,wcontent,cofgen(1,Node),swfrost,10.d0,
     &                       swsophy,numtab,sptab)
          mfluxtable(lay,count) = mfluxtable(lay,count+1) + 
     &                 0.5d0 * (conduc1 + conduc2) * (phead2 - phead1) 
          phead1 = phead2
          conduc1 = conduc2
        enddo
      enddo

      return

2000  continue

! === calculation of matric flux potential ===================================

      lay = layer(node)
      if (phead .lt. wiltpoint) then
! ---   very dry range 
         mflux(node) = 0.0d0
      elseif (phead .gt. -1.023293d0) then
! ---   very wet range (> -10^0.01)
         mflux(node) = mfluxtable(lay,1)
      else  
! ---   direct access table, with linear interpolation
        logphead = 100.d0*log10(-phead)
        count = int(logphead)
        mflux(node) = (logphead-dble(count))*mfluxtable(lay,count+1) +
     &                (dble(count+1)-logphead)*mfluxtable(lay,count)
      endif

      return
      end 

        
      
      subroutine Samuca(task)
      
      !************************************************************************
      !*     SUBROUTINE SAMUCA
      !************************************************************************
      !*  SAMuCA - AGRONOMIC MODULAR SIMULATOR FOR SUGARCANE     
      !*  Written in Microsoft Visual Studio FORTRAN for PC-compatible machines              
	!*  Author: FABIO R. MARIN Date: 11/10/2010
	!*  Last time edited: 09-08-2015 by Murilo dos S. Vianna 
	!************************************************************************
	!*     This subroutine simulates the growth of the plant using pre-determined
	!*     conditions.Hourly values of temperature and photosyntetically active
	!*     radiation come from WEATHER subroutine and daily values of availability
	!*     of water in the soil come from SW subroutine. This subroutine supplies
	!*     the SW subroutine with daily values of leaf area index (LAI) (FM, Oct 2010).
	!C****************************************************************************
	!
	!*                  LIST OF VARIABLES 
	!-----------------------------------------------------------------------  
	
      use Variables      
      
	implicit none	
      
      character*4  outfile          !4-character name
      character*90 path_out         !path to the output file
      
      integer      task             !Controler of call ordering from the main program
      integer      i                !Counter	
      integer      j                !counter
      integer      co2              !CO2 concentration at atmosphere [ppm]
      integer      rdenshape        !User can choose to use a shapefactor for root density (=1) or not (=0)     
      integer      node
      
      logical      flemerged        !Crop emergence flag
          
      real         agefactor        !Age factor to reduce photosynthesis and PER
      real         chudec           !Termal time for tillering deccay
      real         chuem            !Cumulative heat units to emergence (Degree-days)
      real         chumat           !Termal time for tillering maturation stage
      real         chupeak          !Termal time for tillering peak
      real         chustk           !Cumulative heat units for stalk emergence (Degree-days)
      real         cumla(150)       !Cumulative Leaf Area cm
      real         ddeadlm          !
      real         ddealla          !
      real         deadln           !Dead leaf number
      real         di               !daily accumulated temperature above TB (degree days)
      real         diac             !Cumulative Degree-Days
      real         diaclf           !
      real         diam             !average stem diameter (cm)
      real         dileaf           !
      real         dla              !
      real         dleafdm          !Daily incremental leaf dry mass
      real         dlfdma           !
      real         dnleaf           !Incremental leaf number
      real         dnstk            !
      real         dpercoeff        !
      real         dsuc             !
      real         dw               !incremental total plant dry matter weight (kg m-2)
      real         dwa              !incremental canopy dry matter weight (kg m-2)
      real         dwl              !incremental leaf dry matter weight (kg m-2)
      real         dwr              !incremental root dry matter weight (kg m-2)
      real         dws              !incremental stalk dry matter weight (kg m-2)
      real         dwsuc            !incremental sucrose mass (kg m-2)
      real         dwwater          !incremental stalk fresh matter weight (kg m-2)
      real         e                !conversion efficiency of CH2O to plant tissue (g g-1)
      real         epp              !
      real         esred            !
      real         esw              !
      real         etp              !
      real         extcoef          !Extintion Coeficcient
      real         hour(24)         !hour counter
      real         internode(100,4) !internode array
      real         la               !leaf area
      real         la_stk           !leaf area stalk
      real         lgpf             !Light gross photosynthesis
      real         li               !Light Interception - non-dimensional factor to be multoplied be PAR
      real         ln               !Number of Green Leaves (# Leaves/stalk)
      real         ln_stk           !
      real         lntotal          !number of leaves (green + senesced)
      real         maxgl            !maximum number of leaves
      real         mla              !Maximum leaf Area (cm2) as a function of leaf number
      real         noden            !
      real         nstk             !
      real         nstkzero         !Printing (fake) variable to write the number of stalks before emergence. In real, it is 0.1 even before to allow plant start to grow.
      real         par              !
      real         perday           !plant elongation rate per day
      real         pg               !canopy gross photosynthesis rate (t ha-2 day-1)
      real         pho(4)           !
      real         phyloc           !
      real         plantdepth       !
      real         pleng            !
      real         pol              !ratio between WSUC and WSFRESH (%)
      real         popmat           !
      real         poppeak          !
      real         rdepth           !
      real         resp             !
      real         rgpf             !
      real         rootdmzero       !
      real         rowsp            !row spacing
      real         rue              !
      real         sgpf             !
      real         shootdepth       !Shoot depth before emergence
      real         sla              !specific leaf area (m2 kg-1)
      real         srad             !Daily solar radiation (MJ m-2)
      real         srl              !
      real         stalkgpf         !
      real         stkdmc           !Dry matter fraction in the stalk fresh mass
      real         stress_k         !
      real         stress_n         !
      real         stress_pho       !
      real         suc_stk          !
      real         sucmax           !
      real         swface           !soil water excess stress factor, 1
      real         swfacp           !soil water deficit stress factor, 1
      real         tb               !base temperature above which growth occurs (oC)
      real         tbleaf           !
      real         thour(24)        !
      real         tmax             !Daily maximum temperature (oC)
      real         tmin             !Daily manimum temperature (oC)
      real         t_mean            !Daily mean temperature (oC)
      real         trwu             !
      real         tstress          !
      real         w                !total plant dry matter weight (ton ha-1)
      real         w_stk            !
      real         wa               !aerial dry matter weight (ton ha-1)
      real         wa_stk           !
      real         wl               !
      real         wr               !root dry matter weight (ton ha-1)
      real         wr_stk           !
      real         wsdm             !stalk dry matter weight (ton ha-1)
      real         wsfresh          !stalk fresh matter weight (ton ha-1)
      real         wstkwat          !
      real         wsuc             !
      real         wsuc_old         !
      real         qrotts           !
      real         rwuep1
      real         rwuep2
      real         rue_modifier
      real         rdm
      real         ptrats
      real         epp_cm
      real*8       afgen
      real*8       watcon
      real*8       depth
      real*8       rootdis(202)
      real*8       sum
      real         trwuafter
      real         rootshape        !Factor related to the root density shape in soil profile
      real*8       top
      real*8       bot
      real         prwu(macp)
      real         tqropot
      real         rld(macp)
      real         rgf(macp)
      real         wpp(macp)
      real         fcp(macp)
      real         stp(macp)
      real         rootsene
      real         swcon1
      real         swcon2(macp)
      real         swcon3
      real         laythi(macp)
      real         prwulay(macp)
      real         wuf
      real        eoratio
      real        kc_min
      real        lai_max
      
      real*8, parameter :: hwpp = -15000.
      real*8, parameter :: hfcp = -330.
      real*8, parameter :: hstp = -1.d-5
              
      save      
             
      GOTO (1000,2000,3000) task

1000  CONTINUE
      
      !Use the below call To read parameters from a file.crp and initialize variables 
     !! call readsamuca (crpfil,pathcrop,tdwi,laiem,rgrlai,slatb,
     !!&  ssa,span,tbase,kdif,kdir,eff,amaxtb,tmpftb,tmnftb,cvl,cvr,cvs,
     !!&  q10,rml,rmr,rms,rfsetb,frtb,fltb,fstb,perdl,rdrrtb,
     !!&  rdrstb,hlim1,hlim2u,hlim2l,hlim3h,hlim3l,hlim4,rsc,adcrh,adcrl,
     !!&  cofab,rdi,rri,rdc,rdctb,rlwtb,logf,schedule,cumdens,
     !!&  flsolute,ecmax,ecslop,c2eca,c2ecb,c2ecf,numlay,dateharvest,
     !!&  swharvest,dmharvest1,dmharvest2,swgrazing,grazingfactor,
     !!&  daylastharvest,dmlastharvest,wrtmax,                            ! Nwgrassland
     !!&  nsuptab,dmfac,relnitab,nsupply,                                 ! Nwgrassland
     !!&  swcf,swetr,cftb,chtb,cfet,alphacrit,
     !!&  swroottyp,wiltpoint,rootradius,rootcoefa,rsw)                   ! NwRootExtr

                         
      !Initializing crop state variables
      
          agefactor  = 1.
		lntotal    = 0.0
		ln         = 0.0
		deadln     = 0.
		la         = 0.0047
		lai        = 0.
		w          = 0.
		wr         = 0.1
		wa         = 0.
		wsdm       = 0.
		wsuc       = 0.
		wsuc_old   = 0.
		noden      = 0.
		nstk       = 0.1
		wsfresh    = 0.
		stkdmc     = 0.
		pol        = 0.
		pleng      = 0.
		dwsuc      = 0.
		dwwater    = 0.
		wstkwat    = 0.
		stkdmc     = 0.
		wsfresh    = 0.
		wstkwat    = 0.
		shootdepth = 0.
          rdi        = 0.
          diac       = 0.
          di         = 0.
          diaclf     = 0.
          rootsene   = 0.
          swfacp = 1.0d0  
          swface = 1.0d0
          wuf = 1.
          
          
          flemerged = .FALSE.
          
          do i = 1,macp
              rgf(i) = 0.
              rld(i) = 0.
              wpp(i) = 0.
              fcp(i) = 0.
              stp(i) = 0.
              swcon2(i) = 0.
              laythi(macp) = 0.
              prwulay(macp) = 0.
          end do
          
      !-------------------------------------------------
      !-------INPUTS VARIABLES FOR THE MODEL------------
      !---IT SHOULD BE READABLE FORM A FILE IN FUTURE---
      !-------------------------------------------------
          !following the crop file for grass
          !part1
          idev = 1
          !lcc = 366
          
          !part2
          kdif = 0.75
          kdir = 0.75
          
          !part3
          swgc = 1
          
          !part4
          swcf = 2
          albedo = .23
          rsc = 70.
                    
          !part5
          rdc = 200.
          swroottyp = 1
          alphacrit = 1.0d0
          
          !Part6
          
          !Part7 ---> soil water extraction by plant roots (Feddes approach)
          
          HLIM1  =     -0.01  ! No water extraction at higher pressure heads, [-100..100 cm, R]
          HLIM2U =     -1.0  ! h below which optimum water extr. starts for top layer, [-1000..100 cm, R]
          HLIM2L =     -1.0  ! h below which optimum water extr. starts for sub layer, [-1000..100 cm, R]
          HLIM3H =    -200.0  ! h below which water uptake red. starts at high Tpot, [-10000..100 cm, R]
          HLIM3L =    -800.0  ! h below which water uptake red. starts at low Tpot, [-10000..100 cm, R]
          HLIM4  =   -15000.0 ! No water extraction at lower pressure heads, [-16000..100 cm, R]
          ADCRH  =       0.5  ! Level of high atmospheric demand, [0..5 cm/d, R]     
          ADCRL  =       0.1  ! Level of low atmospheric demand,  [0..5 cm/d, R]     
          
          
          
          !Part8 - Salt stress simulation
          !part9
          
          swinter = 1
          cofab = 0.25
          
          !part10 - root density as function of root depth        
          
                    
          !Samuca Inputs
          
          rowsp =         140. !row-spacing [cm]
          plantdepth =     25. !plant depth on planting
          co2 = 390 ! ppm - Current Scenario 
		          
          !Samuca crop parameters
          
          maxgl       =    5.500
          tb          =   11.053
          rue         =    2.207
          sla         =   50.000
          extcoef     =    0.700
          sgpf        =    0.880
          dpercoeff   =    0.323
          sucmax      =    0.950
          srl         =   18.000
          chustk      = 1100.000
          chupeak     = 1200.000
          chudec      = 1500.000
          chumat      = 2500.000
          popmat      =   11.000
          poppeak     =   20.000
          phyloc      =  169.000
          mla         =  450.000
          kc_min      =   0.75
          eoratio     =   1.15
          lai_max     =   mla * poppeak * maxgl / 10000
          
      !-------------------------------------------------    
      !------------END OF INPUTS------------------------    
      !-------------------------------------------------    
            
          
      ! --- maximum rooting depth & actual rooting depth
          rdi = plantdepth
          rdm = min(rds,rdc) !rdc is set on crop file and rds in master file
          rd = min(rdi,rdm)
          rdpot = min(rdi,rdm)
          
                 
      ! --- CALCULATE NORMALIZED CUMULATIVE ROOT DENSITY FUNCTION (Murilo Vianna)
      !User can choose the root density shape type
      !rdenshape = 1 use shape factor
      !rdenshape = 0 provide the table with relative root density as function of relative root depth
          
      rdenshape = 1
      if (rdenshape .eq. 1) then                                        ! NwRootExtr
                
          do i = 1,22,2
          rdctb(i) = real(i)/22.                    
          end do
          
          rootshape = 2.
          do i = 2,22,2
          rdctb(i) = (1-rdctb(i-1))**rootshape
          end do
          
      else
      
          rdctb(1) = 0.
          rdctb(2) = 1.
          
          rdctb(3) = 0.25
          rdctb(4) = 0.9
          
          rdctb(5) = 0.5
          rdctb(6) = 0.7
          
          rdctb(7) = 0.75
          rdctb(8) = 0.4
          
          rdctb(9) = 1.
          rdctb(10)= 0.01
          
      end if
          
      ! ---   specify array ROOTDIS with root density distribution
      do i = 0,100
          
          depth = 0.01d0 * dble(i)
          rootdis(i*2+1) = depth
          rootdis(i*2+2) = afgen(rdctb,22,depth)
                   
      enddo     
      
      ! ---   calculate cumulative root density function
      do i = 1,202,2
      ! ---     relative depths
          cumdens(i) = rootdis(i)
      enddo
      
        sum = 0.d0
        cumdens(2) = 0.d0
        do i = 4,202,2
      ! ---     cumulative root density
          sum = sum + (rootdis(i-2)+rootdis(i)) * 0.5d0
     &               * (cumdens(i-1)-cumdens(i-3))
          cumdens(i) = sum
      enddo
      
      ! ---   normalize cumulative root density function to one
        do i = 2,202,2
          cumdens(i) = cumdens(i) / sum
        enddo         
	          
		do i = 1,100
			internode(i,1) = 0.
			internode(i,2) = 0.
			internode(i,3) = 0.
			internode(i,4) = 0.
		end do
			
		do i = 1,150
			cumla(i) = 0.
          end do
          
      !Potential Root water uptake (RWU) 
          
      ! --- calculate swc in wilting point(hlim4), field capacity, saturated point and swcon2 for potential rwu    
          do node = 1, macp
          wpp(node) = watcon(node,hwpp,cofgen(1,node),swsophy,
     &    numtab,sptab)          
          fcp(node) = watcon(node,hfcp,cofgen(1,node),swsophy,
     &    numtab,sptab)    
          stp(node) = watcon(node,hstp,cofgen(1,node),swsophy,
     &    numtab,sptab)    
          
              if (wpp(node) .GT. 0.30) then
              swcon2(node) = 45.0
              end if
          
          swcon2(node) = 120. - 250. * wpp(node)
          
          !swcon2(node) = 62.
          
          laythi(node) = abs(abs((z(node)+ 0.5*dz(node))) - 
     &     abs((z(node) - 0.5*dz(node))))
                      
          end do
       
          swcon1 = 1.32E-3 !L/?/?
          swcon3 = 7.01    !L/?/?          
          
          !swcon1 = 2.67E-3 !L/?/?
          !swcon3 = 6.68    !L/?/? 
          
                    
      !--- Linking with SWAP variables    
          
      tmax = tmx                  !(oC)
      tmin = tmn                  !(oC)
      t_mean = .5 * (tmax + tmin) !(oC)
      par = (rad * .5)/0.1E7      !(MJ/m2.d
      epp = ptra                  !(cm/d)
      ch = pleng * 100.           !(cm) Used for Penmam-M. method for et0
          
      !Opening output file...
          
      path_out = 'C:\Murilo\DOUTORADO\WUR\SWAP_Sugarcanev1\Data'
      
      OPEN (45,FILE='C:\Murilo\DOUTORADO\WUR\SWAP_Sugarcanev1\Data\
     &Plant_'//trim(project)//'.OUT',STATUS='REPLACE',RECL=1180)
      
      OPEN (44,FILE='C:\Murilo\DOUTORADO\WUR\SWAP_Sugarcanev1\Data\
     &Plant_TEST.OUT',STATUS='REPLACE',RECL=1180)
          
          !Writing the outputfile header
            
		WRITE(45,11) '#Simulating for ', project
		WRITE(45,13) 
		WRITE(45,14)
		WRITE(45,15)
		WRITE(45,16)
		WRITE(45,17)
	
11        FORMAT (2A,/) 		
12        FORMAT (2A,/)                                                  
13        FORMAT('Results of plant growth simulation in daily output:')				
14        FORMAT('        Day   Days   Days  Accum    Plant   Stalk',
     &'    Root   Stalk  Sucros                          ',
     &' Plant   Green   Stress  Stress         ')!    Stalk')
15        FORMAT('         of  after  after  D.Day   Weight  Weight',
     &'  Weight  Weight  Weight     POL     LAI  Tiller  Height',
     &'  Leaves   Factor  Factor        ')!   Diamet')      
16        FORMAT('Year   Year  Simul  Plant   oC.d   DMt/ha  DMt/ha',
     &'  DMt/ha  FMt/ha  DMt/ha       %   m2/m2    #/m2       m',
     &'   #/Stk   Expans  Photos       ')!       cm')
17        FORMAT('----  -----  -----  -----  ------ -------  ------',
     &'  ------  ------  ------  ------  ------  ------  ------',
     &'  -------  ------  ------        ')       
	
       
	RETURN      
      
2000      CONTINUE
          
      !Linking with SWAP Variables
          
      tmax = tmx                  !(oC)
      tmin = tmn                  !(oC)
      t_mean = .5 * (tmax + tmin) !(oC)
      par = (rad * .5)/0.1E7      !(MJ/m2.d)
      epp = ptra                  !(cm/d)
      ch = pleng * 100.           !(cm) Used for Penmam-M. method for et0
      
      
	!************************************************************************
	!     EMERGENCE OF PLANT - Calculates the thermal time for emergence
	!       and trigger the germination when Plant Shoot reaches the soil surface 
	!************************************************************************

	!Calculating degree-days using a unique tb
          
		if (t_mean .GE. tb .and. t_mean .LE. 42) THEN
			di = (t_mean-tb)     ! Degree-Days for today
		    ELSE 
			    di = 0.0
                  !t_mean = tb
		    ENDIF
		
		diac = diac + di         !Cumulative Degree-Days
          
          !Crop coeficient (Kc) - Based on Field Experiment data for RB867515          
          cf = kc_min + (eoratio - kc_min) * lai / lai_max
          cf = max(0.,cf)
		      
          IF (.not. flemerged) THEN
          	
		shootdepth = plantdepth - diac * 0.08 !Assuming the primary shoot growth rate of 0.08 cm/DG (Keating et. al. 1999)
		
		    IF (shootdepth .LT. 1.d-8) THEN
			    flemerged = .TRUE.          !Emerged!
			    shootdepth = 0.
              END IF
              
          END IF	    
			
      !--------------LEAF DEGREE-DAYS------------- 
		
		tbleaf = tb !Could different tb for leaf fenology
						
		IF (t_mean .LT. tbleaf) THEN  !Never will be possible with above condition...
			!t_mean = tbleaf
			dileaf = 0.! Degree-Days
		ELSE 
			dileaf = (t_mean-tbleaf) 
          ENDIF		

		diaclf = diaclf + dileaf	     
         
      !*****************************************************************************
	!*     Subroutine PGS
	!*     Calculates the canopy gross photosysntesis rate (PG)
	!*     Input: SWFACP,PAR,LAI,DIAC,TMN,W,EXTCOEF)
	!*     Output: PG,RESP
	!******************************************************************************
	!-----------------------------------------------------------------------  
            
      
	!-----------------------------------------------------------------------
	!     Calculate light interception of the canopy
	!-----------------------------------------------------------------------
		li    = 1. - EXP(-extcoef * lai)
          
	!Age factor to reduce the photosynthesis due to age of crop
		agefactor = EXP(-0.000401*(diac-4000)) ! Changed to reduce the photosynthesis rates after a while
		agefactor = MIN(1.,agefactor)      
			
		!Computing the Temperature Stress on Photosynthesis
		!########Aqui é necessário computar o numero de horas em as temperaturas estiveram acima desses valores.
		!Importante utilizar método de Thoraria de Parton & Logan
		
          CALL THOURS(tmax,tmin,thour,tstress)!check if this thresholds are correct for gross photsynthesis
          		              
      !Computing gross photosynthesis without stress during the first days
      
		swfacp = 1. !No water stress (potential photosynthesis)
		rue_modifier =  ((0.0001282051  * co2) + 0.95)
		pg = rue * rue_modifier * li * par  * tstress   * swfacp  !* AGEFACTOR
		!  (g /MJ * dimensionless * dimensionless * MJ/m2 [soil] d-1  * m2 [leaf] / m2 [soil])
		pg = pg  / 100.              ! (converting to t / ha d)
		
                 
      !*****************************************************************************
	!*     Subroutine PARTS
	!*     Calculates the Partitioning rules for stalk, root and leaf
	!******************************************************************************
		
		IF (diac .LT. chustk) THEN  !.AND. DIAC > CHUEM
				rgpf = .3			
				lgpf = .1
				!stalkgpf = .9		
		ELSEIF (diac .GE. chustk) THEN
				stalkgpf = sgpf
				lgpf = (1-(sgpf)) 
				rgpf = .2				
          END IF
                                  
      
      !Calculating root depth as function of the cumulative degree-days
      !PLANTDEPTH = 25.0   ! PLANTING DEPTH
      !Laclau and Laclau (2009) rate of deepening of the root front (0.53 cm day-1 or 0.048 cm oC-1 day-1) over the first 4 months after planting, and an increase thereafter to 1.75 cm day-1 (0.22 cm oC-1 day-1) in the irrigated crop and 1.86 cm day-1 (0.24 cm oC-1 day-1)
            
      IF (diac .LT. 1000) THEN   ! A variation in RDEPTH should be computed dur the water stress - the higher WS, the deeper the root goes
        rdepth = (0.08 * diac) + plantdepth
        !RDEPTH = 0.048 * DIAC !- Original results from Laclau & Laclau (2009)
      ELSE
        rdepth = (0.08 * diac) + plantdepth
        !RDEPTH = .22 * DIAC !- Original results from Laclau & Laclau (2009)
      END IF
            
      rd = min(rdepth,rdm)
                     
      RETURN
      
3000    CONTINUE    
          
      IF(flemerged) THEN
          
        swfacp = 1.0d0  !Soil water stress factor on photosynthesis
        swface = 1.0d0  !soil Water Stress factor on growth
        rwuep1 = 2.0d0  !Root water uptake transpiration ratio to start growth stress
        rwuep2 = 1.0d0  !Root water uptake transpiration ratio to start photosynthesis
        
        !!Linking the soil water balance with SWAP
        !        
        !IF (epp .GT. 0.0001) THEN
        !      IF ((trwu/epp_cm) .LT. rwuep1) THEN
        !          swface =  MIN(1.,(1/rwuep2*trwuafter / epp_cm))
        !      END IF
        !      IF ((trwu/epp_cm) .LT. rwuep2) THEN
        !          swfacp = MIN(1.,(1/rwuep1*trwuafter / epp_cm))
        !      END IF
        !ENDIF
        
        ! --- reset root water extraction array
        do node = 1,numnod
          prwu(node) = 0.0d0
        end do
        tqropot = 0.
        
        ! --- determine lowest compartment containing roots
          node = 1
          do while ((z(node) - dz(node)*0.5d0) .gt. (-rd + 1.d-8))
              node = node + 1
          end do
          noddrz = node
        
        ! ---  root senescence
        rootsene = 1- ((.15 * wr) / (wr + 1)) 
        
        ! ---   calculate root length density
        do node = 1,noddrz
          top = abs((z(node) + 0.5*dz(node)) / rd)
          bot = abs((z(node) - 0.5*dz(node)) / rd)
      
          rgf(node) = 
     &        (afgen(cumdens,202,bot)-afgen(cumdens,202,top))
                              
          rld(node) = wr * srl * rgf(node) * rootsene/ laythi(node)                 
        enddo
          
       ! --- potential RWU 
       do node = 1, noddrz
          prwu(node) = max(0.,min(0.07, 
     &  swcon1*EXP(MIN((swcon2(node)*(theta(node)-wpp(node))),40.))/
     &     (swcon3-LOG(rld(node)))))
          
        prwulay(node) = prwu(node) * laythi(node) * rld(node)
          
          tqropot = tqropot + prwulay(node)
       enddo
      
      ! --- water stress factors in growth and photosynthesis 
      if (ptra .le. 1.d-5 .OR. tqropot .le. 1.d-5) then
        wuf = 0.
        swfacp = 1.
        swface = 1.
      else
        wuf = max(tqropot/ptra,0.)
        
          if (wuf .lt. rwuep1) then          
          swface = max(0.,min((1./rwuep1) * wuf,1.))      
          else
          swface = 1.
          endif
          
          if (wuf .lt. rwuep2) then          
          swfacp = max(0.,min((1./rwuep2) * wuf,1.))      
          else
          swfacp = 1.
          endif         
        
      endif          
      
          pg = pg * swfacp
          
         IF (diac .LT. chustk) THEN			! Growing primary without tillers
			dnstk = 0.
			ln = .1                     !In doubt because the energy retained in the stalk
			lntotal = ln			                  
	
			dw =  pg/1.270      
			!dW =  PG/1.420    !conveting from (kg [Glucose] m-2 ground) to (kg [Plant Material] m-2)
			dw =  MAX(dw, 0.)
              
	
              CALL LAIS(daynr,dileaf,diaclf,mla,chustk,dnstk,nstk,ln,
     & dnleaf,maxgl,swface,sla,dla,dleafdm,wl,dw,lgpf,rgpf,sgpf,
     & ddeadlm,deadln,lntotal,phyloc,li,ddealla,cumla)
			
			perday = 0.
			dwa = (1. - rgpf) * dw  ! kg/m2ground
			dwr = rgpf * dw        ! kg/m2ground 
			dws = 0.
			dwsuc = 0.
			
		ELSE	!Growing prymary stalk plus tillers 
			
			CALL DIAMPERS(daynr,swface,noden,thour,nstk,wsdm,dnleaf,
     & wsuc,lntotal,ln,pleng,diac,chustk,perday,diam,internode,
     & sucmax)
			! Tillers Rate calculated in Stalks/m2 - That's why they were divided by ROWSP
			
			IF (diac .GT. chustk .AND. diac .LT. chupeak) THEN 
				dnstk = (poppeak/(chupeak-chustk)) * di !/ ROWSP				   ! initial tiller grow
			ELSEIF (diac .GT. chupeak .AND. diac .LT. chudec) THEN
				dnstk = 0.								                           ! peak tiller 
			ELSEIF (diac .GT. chudec .AND. diac .LT. chumat) THEN
				dnstk = (-(poppeak - popmat) / (chumat-chudec)) * di !/ ROWSP	   ! reduction phase to mature
			ELSEIF (diac .GT. (chumat)) THEN
				dnstk = 0.		            						               ! late stable tiller population
			END IF
	
			dw =  pg/1.270    !See Scarpare et al. (2012) e Machado (1981)  
			!dW =  PG/1.420      !conveting from (kg [Glucose] m-2 ground) to (kg [Plant Material] m-2). See sheet Plant Composition.xls
			dw =  MAX(dw, 0.)
              
              			
			CALL LAIS(daynr,dileaf,diaclf,mla,chustk,dnstk,nstk,ln,dnleaf,maxgl,
     & swface,sla,dla,dleafdm,wl,dw,lgpf,rgpf,sgpf,ddeadlm,
     & deadln,lntotal,phyloc,li,ddealla,cumla)
			
			dwa = (1.- rgpf) * dw
			dwr = rgpf * dw
			dws = stalkgpf * dw 
			dwsuc = wsuc - wsuc_old
		ENDIF
		
		IF (wsdm .GT. 0.) THEN
	
			IF (swfacp .LT. 0.9 .AND. dws .LT. 0.5) THEN 
				
				IF ((wstkwat/wsfresh) .GT. 0.75 ) THEN 
					
					dwwater = -( 0.003 * wstkwat) 
					!write (*,*) "umido"
				
				ELSEIF ((wstkwat/wsfresh) .GT. 0.7 ) THEN 
					
					dwwater = -( 0.003 * wstkwat) 
					!write (*,*) "umido"
					
				ELSEIF  ((wstkwat/wsfresh) .GT. 0.6 ) THEN 
					dwwater = -( 0.002 * wstkwat) 
					!write (*,*) "medio"
					
				ELSEIF  ((wstkwat/wsfresh) .GT. 0.45 ) THEN 
					dwwater = -( 0.001 * wstkwat) 
					!write (*,*) "medio"
	
				ELSEIF  ((wstkwat/wsfresh) .LT. 0.45 ) THEN 
					dwwater = 0.0
				!write (*,*) "seco"
				ENDIF
			
			ELSE 
				dwwater = (3.607 * dws) - (3.078 * dwsuc)
			END IF
	
		ELSE
			dwwater = 0.
			stkdmc = .3    !Conteúdo de massa seca na planta quando ainda não há colmo.
		ENDIF    
	
		! integration of daily state variable
		
		nstk      = nstk    + dnstk			         !Number of stalks
		ln        = ln      + dnleaf                    !Number of green leaves per Stem
		la        = la      + dla - ddealla                     
		lai       = la      * nstk
		w         = w       + dw
		wa        = wa      + dwa
		wr        = wr      + dwr
		wsdm      = wsdm    + dws
		wstkwat   = wstkwat + dwwater
		wl        = wl      + dleafdm
		stkdmc    = wsdm    / (wsdm + wstkwat)
		wsfresh   = wsdm    + wstkwat
		
		IF (noden .GT. 1) THEN
			wsuc = 0.
			pleng = 0.
			DO i = 1,noden
				wsuc = wsuc + internode(i,3)
				pleng = pleng + internode(i,4)
			END DO
		ELSE
			wsuc = 0.
          END IF         
                 
          
		!------------------
		! Armazendando o valor de Wsuc para que no próximo dia (Wsuc_Old) possa ser calculada a taxa de incremento
		! e com ela a taxa de variação do conteúdo de água no colmo (dWSuc)
	
	
		wsuc_old = wsuc
	
		la  = MAX(la,0.0)
		w   = MAX(w, 0.0)
		wa  = MAX(wa, 0.0)
		wr  = MAX(wr,0.0)
		wsdm  = MAX(wsdm,0.0)
		wsfresh = MAX(wsfresh,0.0)
		wl  = MAX(wl,0.0)
		pleng = MAX(pleng,0.0)
		
		
		!****************************
		! Calculating the Fresh Mass and POL% of Cane based on Cane DM and Suc Content
		! after Martines SASTA 2001. Apud Canegro (2008). This should be included in a Suroutine for clarity of the software.
		! POL has been calculated by fresh cane mass
		! Calculates the water loss when stalk dry mass rate is zero due to water stress
		! Daily rate of water loss when dWS is zero: (1-0,998523) - From Ap. Taboado Experiment. (FM - March of 2014)
	
		IF (wsuc==0. .OR.  wsfresh==0.) THEN
			pol = 0.
		ELSE
			pol = wsuc / wsfresh * 100
		END IF
	     
          
        ENDIF
	!************************************************************************
	!     OUTPUT
	!************************************************************************
	        
			IF (nstk == .1) THEN
				nstkzero = 0.
				rootdmzero = 0.
				WRITE(45,9) iyear,daynr,daycum,daycrop,diac,w,wsdm,rootdmzero,
     & wsfresh,wsuc,pol,lai,nstkzero,pleng,ln,swface,swfacp
			ELSE
				WRITE(45,9) iyear,daynr,daycum,daycrop,diac,w,wsdm,wr,wsfresh,wsuc,
     & pol,lai,nstk,pleng,ln,swface,swfacp
9           FORMAT(I4,4X,I3,4X,I3,4X,I3,6F8.1,4F8.2,9F8.2)
              ENDIF
              
              
      
      IF (flcropend) THEN
          CLOSE(45)
          CLOSE(44)
      END IF
      
	
	    RETURN
      
		END SUBROUTINE Samuca
	!************************************************************************
	
		
	SUBROUTINE PGS(swfacp,chustk,par,lai,diac,t_mean,w,extcoef,pg,rue,li)
	IMPLICIT NONE
	
	!*****************************************************************************
	!*     Subroutine PGS
	!*     Calculates the canopy gross photosysntesis rate (PG)
	!*     Input: SWFACP,PAR,LAI,DIAC,TMN,W,EXTCOEF)
	!*     Output: PG,RESP
	!******************************************************************************
	!-----------------------------------------------------------------------  
		

          
		REAL par
		REAL lai
		REAL diac
		REAL li                      !Canopy Light interception
		REAL pg
		REAL e_fac
		REAL pratio
		REAL ct
		REAL swface
		REAL swfacp
		REAL rowspc
		REAL w
		REAL agefactor
		REAL extcoef
		REAL chudec
		REAL chustk
		REAL t_mean
		REAL rue
		REAL tstress
		REAL co2
		REAL rue_modifier
          REAL epp
          REAL trwu
          REAL rwuep2
		REAL rwuep1
		REAL rdm
		REAL qr
		REAL trwuafter
          REAL ptrats
          REAL epp_cm
          
		
		
	co2 = 390 ! ppm - Current Scenario ALTERAR TAMBÉM OS VALORES QUE CONSTAM NO MODULO SOIL WATER - EPs
	!CO2 = 550 ! ppm - Scenario 4 - ALTERAR TAMBÉM OS VALORES QUE CONSTAM NO MODULO SOIL WATER - EPs
	!CO2 = 450 ! ppm - Scenario 3 - ALTERAR TAMBÉM OS VALORES QUE CONSTAM NO MODULO SOIL WATER - EPs
	!-----------------------------------------------------------------------
	!     Calculate maximum photosynthesis as function of PAR, g CH2O/m2
	!-----------------------------------------------------------------------
		par = par !* 3.969354729                                   !The multiplier 3.9 for PAR is to convert Mj to micromol. In the current version RUE is informed in terms of g/MJ
																	!For CropGro PhotoMode, see Sheet c:\canemodel\Photosynthesis_Cropgro.xls
	!-----------------------------------------------------------------------
	!     Calculate reduction in photosynthesis due to incomplete canopy.
	!-----------------------------------------------------------------------
		li    = 1. - EXP(-extcoef * lai) 
          
          
        !Block pasted from the soilwater balance subroutine of SAMUCA
        !********************************************************************  
        swfacp = 1.0  !Soil water stress factor on photosynthesis
        swface = 1.0  !soil Water Stress factor on growth
        rwuep1 = 1.00 !Root water uptake transpiration ratio to start growth stress
        rwuep2 = 2.00 !Root water uptake transpiration ratio to start photosynthesis
        
        !Linking
        
        epp = ptrats
        
             
        IF (epp .GT. 0.0001) THEN
              IF ((trwu/epp_cm) .LT. rwuep1) THEN
                  swface =  MIN(1.,(1/rwuep2*trwuafter / epp_cm))
              END IF
              IF ((trwu/epp_cm) .LT. rwuep2) THEN
                  swfacp = MIN(1.,(1/rwuep1*trwuafter / epp_cm))
              END IF
        ENDIF
        !********************************************************************   
        
	
	!Compute daily gross photosynthesis (g CH2O/m2 leaf/d)
		agefactor = EXP(-0.000401*(diac-4000)) ! Changed to reduce the photosynthesis rates after a while
		agefactor = MIN(1.,agefactor)
	
		!E_FAC - Effect of any stress on canopy photosynthesis (0-1). KEEP 1 UNTIL ADVANCES
		IF (swfacp .LT. .99) THEN
				e_fac = agefactor * swfacp !Changed in 11/5/2010 to reduce the photosynthesis rates after a while
			ELSE
				e_fac = agefactor 
		ENDIF
	
	!###############AgeFactor must be included as soon as possible to explain the late decrease in leaf photosynhtesis due the decrase in leaf N
		
		!Computing the Temperature Stress on Photosynthesis
		!########Aqui é necessário computar o numero de horas em as temperaturas estiveram acima desses valores.
		!Importante utilizar método de Thoraria de Parton & Logan
		
		IF (t_mean .LT. 15) THEN
			tstress = -0.50 + 0.10*t_mean
		ELSEIF (t_mean .GE. 42.) THEN
			tstress = -0.0666667*t_mean + 3.33333
		ELSE 
			tstress = 1.
		ENDIF
		!--------------    
		!Computing gross photosynthesis without stress during the first days
		
		rue_modifier =  ((0.0001282051  * co2) + 0.95)
		pg = rue * rue_modifier * li * par  * tstress   * swfacp  !* AGEFACTOR
		!  (g /MJ * dimensionless * dimensionless * MJ/m2 [soil] d-1  * m2 [leaf] / m2 [soil])
		
		pg = pg  / 100              ! (converting to t / ha d)
		!write (*,*) 'Simulating for Atmosphere in Plant Module', CO2,  RUE_MODIFIER
	
		RETURN
		END SUBROUTINE PGS
	
	
	
	SUBROUTINE PARTS(diac,chustk,chuem,lgpf,sgpf,rgpf,stalkgpf)
	
      Use Variables
      
      IMPLICIT NONE 
	!*****************************************************************************
	!*     Subroutine PARTS
	!*     Calculates the Partitioning rules for stalk, root and leaf
	!******************************************************************************
	!-----------------------------------------------------------------------  
		
          
		REAL diac
		REAL chustk
		REAL chuem
		REAL lgpf
		REAL sgpf
		REAL rgpf
		REAL contador
		REAL minrgpf
		REAL maxlgpf
		REAL stalkgpf
		
          save
          
		IF (diac .LT. chustk) THEN  !.AND. DIAC > CHUEM
				rgpf = .3			
				lgpf = .1
				!stalkgpf = .9
		
		ELSEIF (diac .GE. chustk) THEN
				stalkgpf = sgpf
				lgpf = (1-(sgpf)) 
				rgpf = .2
				
	!#######################retirar ao finalizar
				!if ((rgpf + lgpf + stalkgpf)==1.) then
				!    continue
				!else
				!    pause
				!endif
				
		END IF
		RETURN
		
		END SUBROUTINE PARTS
	
	
	
	SUBROUTINE LAIS(daynr,dileaf,diaclf,mla,chustk,dnstk,nstk,ln,dnleaf,
     & maxgl,swface,sla,dla,dleafdm,wl,dw,lgpf,rgpf,sgpf,ddeadlm,deadln,
     & lntotal,phyloc,li,ddealla,cumla)
	IMPLICIT NONE
	!-----------------------------------------------------------------------  
		
		REAL ln
		REAL dla
		REAL swfacp
		REAL dnleaf
		REAL sla
		REAL la
		REAL agefactor
		REAL chustk
		REAL e_fac
		REAL cumla(150)
		REAL dnstk
		REAL dleafdm
		REAL dw
		REAL lgpf
		REAL maxgl
		REAL sgpf
		REAL wl
		REAL phyloc
		REAL li
		REAL swface
		REAL di
		REAL nstk
		REAL ddeadlm
		REAL ddealla
		REAL deadln
		REAL lntotal
		REAL diac
		REAL diaclf !acumulated degree-days based on leaf tb temperature
		REAL dileaf !degree-days based on leaf tb temperature
		REAL tbleaf !base temperature for leaf development
		REAL rgpf
		REAL mla
		INTEGER daynr
		INTEGER i
		
	!-----------------------------------------------------------------------  
		!Calculating the Age effect on leaf production
		!####### include a function or explian the function below 
		agefactor = EXP(-0.000401*(diaclf-(4000))) ! it'll be kept equal 1 untill N Stress algorithm had been finished
		agefactor = MIN(1.,agefactor)
		!IF (SWFACE < .99) THEN
		e_fac = swface !* AGEFACTOR
	! ELSE
		!   E_FAC = SWFACE
		!ENDIF		
		
	
		!Calculating the daily leaf number increment 
		!#leaves/m2 = leaf/oC.d  *  oC.d * 
		dnleaf = 1/phyloc * dileaf
	!------------
	! Limiting the maximum number of leaves to MAXGL
		IF (ln .GT. maxgl) THEN     !* NSTK)
			
			ln = ln - 1	   	  
			deadln  = deadln + 1 !lntotal will be computed in the integration block 	      
			ddealla = cumla(int(deadln))     !(m2)
				
		ELSE
			
			ddealla = 0.
		
		ENDIF
		
		lntotal = ln + deadln  
	
	!---------------
	!Senescing leaves because of Water Stress
	!	IF (SWFACE < 0.7) THEN
	!	   	  LN = LN - (1 - SWFACE)  
	!	   	  DEADLN  = DEADLN  + (1 - SWFACE) 
	!    ENDIF
	
	!	IF (LI > 0.85) THEN
	!	   	  LN = LN - 1  
	!	   	  DEADLN  = DEADLN  + 1 
	!    ENDIF
	
	!----------------------
	! Calculatin the leaf area increment
	! We are not assuming that the PHYLOC would be affected by the Water Stress. We are just including a 
	! restriction for the expasion of leaves out of the expansion algorithm. In other words, 
	! in spite of our believe still being that Water Stress affects area expansion rather than leaf production
	! our algorthim does not allow to include the water stress easily. 
	!--
		!m2l/m2ground       #leaves/stalk           m2/leaf      
			dla       =        dnleaf           *   (mla/10000)  *   e_fac  !SWFACE-Included by Fabio Marin in 10/22/2010
	
	!---------------------
	!Checking the Needed mass for reach the potential leaf area increment (kg/m2 SOLO))
	!kgleaf/stalk	m2leaf/stk * g/mm2leaf *  (kg/g)   * mm2/m
		
		dleafdm   =  dla     *  (1/sla)  
	!write (*,*) "efac_folha", e_fac, agefactor
		IF (dleafdm .GT. (lgpf * dw)/nstk) THEN	        !checking there is enough biomass for leaf growth 
			dla = ((lgpf * dw)/nstk) * (sla)          !reducing the leaf growth due the lack of photosynthates
		END IF                                              
	
		cumla(int(lntotal+1)) = cumla(int(lntotal+1)) + dla
		dleafdm = dla*(1/sla)                                	   	! Recalculating dLEAFDM after the dLAI has been reduced
	
	! IMPORTANT: remaining to include an stress control to die leaves under severe conditions 
	
	!-----------------------------------------------------------------------  
	RETURN
	END SUBROUTINE LAIS
	!************************************************************************
	
	
	SUBROUTINE THOURS(TMAX,TMIN,THOUR,TSTRESS)
      USE Variables
	IMPLICIT NONE
	!*****************************************************************************
	!*     Subroutine THOURS
	!*     Calculates the Hourly temperature and daily integration of dPER 
	!******************************************************************************
	!-----------------------------------------------------------------------  
		
	
		REAL thour(24)	                           !Temperature related variables !oC
		REAL tsunset                               !Temperature related variables !oC
		REAL tmax                                  !Temperature related variables !oC
		REAL tmin                                  !Temperature related variables !oC
		REAL :: a=2.                               !original constants from Parton and Logan paper
		REAL :: b=2.2                              !original constants from Parton and Logan paper
		REAL :: c=-0.17                            !original constants from Parton and Logan paper
		REAL decsol                                ! astronomic variables
		REAL x                                     ! astronomic variables
		REAL ahn                                   ! astronomic variables
		REAL timnight                              ! time related variables
		REAL timday                                ! time related variables
		REAL tim                                   ! time related variables
		REAL tstress                               !Temperature stress indicator
		REAL ts(24)                                ! Temperature stress indicator
		REAL dj(365) 
		REAL sunset(365)
		REAL sunrise(365)
		REAL qo(365)
		REAL photop(365)
		REAL nigthp(365)    
		REAL, Parameter :: pi = 3.14159265         !trigonometric variables
		REAL, Parameter :: d_2_r = 3.14159265/180. !trigonometric variables
		REAL, Parameter :: r_2_d = 180./3.14159265 !trigonometric variables
		INTEGER n
		          
          save
          
		DO tim=1,24
			ts(tim) = 0.
			tstress = 0.
		ENDDO
	
		! initializing arrays
			DO n = 1, 365
				dj(n) = n
				sunset(n) = 0
				sunrise(n) = 0
				qo(n) = 0
				photop(n) = 0
			ENDDO
		
		!calculating photoperiod
			DO n = 1,365 
				decsol = 23.45 * SIN((d_2_r*(dj(n)-81.)*(360./365.)))
				ahn  = ACOS((-TAN(d_2_r*lat)*TAN(d_2_r*decsol)))
				photop (n) = 2 * (r_2_d*ahn) / 15
				nigthp (n) = 24 - photop(n)
				sunrise(n) = 12 - photop(n)/2
				sunset(n) = 12 + photop(n)/2
              ENDDO		
		
              !
		!Initial Conditions
			DO tim=1,24
				
				IF (daynr .GT. 365) THEN
					n = daynr - 365
				ELSE
					n = daynr
				END IF
			
				!Calculating SUNSET temperature follow Parton & Logan (1981)
				
				tsunset = (tmax - tmin) * sin((pi*photop(n)-c)/(photop(n)+2*a))
     &           + tmin
				! Estimating daylight temperatures
				IF (tim .GT. sunrise(n) .AND. tim .LT. sunset(n)) THEN
					timday=tim-sunrise(n)
					thour(tim) = tmin + (tmax - tmin) * SIN((pi*timday)/(photop(n)
     &                 + 2*a))
				
				ELSE ! nocturnal temperatures
					IF (tim .GT. sunset(n) .AND. tim .LT. 24) THEN
						timnight = tim - sunset(n)
					ELSE
						timnight = (24.-sunset(n)) + tim
					ENDIF
					thour(tim) = tmin + (tsunset-tmin) * EXP(-b * timnight)/(nigthp(n)/
     &                (2*a))
				ENDIF
						
			ENDDO
		
		DO tim=1,24
	
		!Computing the Temperature Stress on Photosynthesis
		!########Aqui é necessário computar o numero de horas em as temperaturas estiveram acima desses valores.
			IF (thour(tim) .LT. 15) THEN
				ts(tim) = -0.50 + 0.10*ts(tim)
			ELSEIF (thour(tim) .GT. 35) THEN
				ts(tim)= -0.0666667*ts(tim) + 3.33333
			ELSE
				ts(tim)= 1.
			ENDIF
			tstress = tstress + ts(tim)
		ENDDO
		
		tstress = max(0.,min(1.,tstress / 24.))
	
	RETURN
	END SUBROUTINE THOURS
	!************************************************************************
	
	
	SUBROUTINE DIAMPERS(daynr,swface,noden,thour,nstk,wsdm,
     & dnleaf,wsuc,lntotal,ln,pleng,diac,chustk,perday,diam,internode,
     & sucmax)
	!-----------------------------------------------------------------------  
	IMPLICIT NONE   !Inserted on 9/8/2010 - finalized in 9/9/2010.
	!*****************************************************************************
	!*     Subroutine DIAMPERS
	!*     Calculates the Hourly dPER and integrates in a daily value PERDAY, 
	!*     Calculates the Mean Diameter of Stalks
	!*     Calculates the Sucrose Mass
	!******************************************************************************
		
		REAL thour(24)	 
		REAL INTERNODE(100,4)            !Declaring the internode array to store: 1)internode number,2)diameter,3)lenght,4)sucrose
		REAL sucpf(200)
		REAL perday
		REAL dper
		REAL agefactor
		REAL chustk
		REAL diac
		REAL diam_stk
		REAL diam
		REAL pleng
		REAL ileng
		REAL nstk
		REAL swface
		REAL dpercoeff
		REAL dnleaf
		REAL ln
		REAL wsuc
		REAL noden
		REAL nodencurrent
		REAL lntotal
		REAL sucpfacum
		REAL fiberpf
		REAL sucmax
		REAL i
		REAL wsdm
		INTEGER h
		!INTEGER year
		INTEGER daynr
		!INTEGER das
		!INTEGER dap
			
		
		!Maximum sucrose content in stalk - Check!!!!!!!!!!!!!!!!
		
		perday = 0.	
		dpercoeff = .3 
		DO h=1,24
			IF (thour(h) .GE. 35) THEN                                       !The 2 equation system is not in use. The same relation has being used in two temp ranges. See PER.xlsx
				dper = (2.95 * thour(h) - 0.8516) * dpercoeff * swface    !.16 means the fraction of plant elongation due to stalk elongation - from Canegro
			ELSEIF (thour(h) .GE. 33 .AND. thour(h) .LT. 35) THEN            !The 2 equation system is not in use. The same relation has being used in two temp ranges. See PER.xlsx
				dper = (1.9625 * thour(h) - 0.8516) * dpercoeff * swface  !.16 means the fraction of plant elongation due to stalk elongation - from Canegro
			ELSEIF (thour(h) .GE. 30 .AND. thour(h) .LT. 33) THEN            !The 2 equation system is not in use. The same relation has being used in two temp ranges. See PER.xlsx
					dper = (0.4875 * thour(h) - 0.8516) * dpercoeff * swface !.16 means the fraction of plant elongation due to stalk elongation - from Canegro
				!  dPER = (0.1625 * THOUR(H) - 0.8516) * .3 * SWFACE !.16 means the fraction of plant elongation due to stalk elongation - from Canegro
			ELSEIF (thour(h) .GE. 25 .AND. thour(h) .LT. 30) THEN 
					dper = (0.325 * thour(h) - 0.8516) * dpercoeff * swface  !.16 means the fraction of plant elongation due to stalk elongation - from Canegro
				! 	dPER = (0.1625 * THOUR(H) - 0.8516) * .3 * SWFACE !Originally .16 !The value .3 was added to put the values closer to Brazilian values.
			ELSEIF (thour(h) .LT. 25) THEN 
					dper = (0.10 * thour(h) - 0.8516) * dpercoeff * swface   !Originally .16 !The value .3 was added to put the values closer to Brazilian values.
			ENDIF
	
				dper = (MAX(0.0, dper))
				agefactor = EXP(-0.000401*(4000)) ! Factor is an age reduction factor for dPER, based on N. G. Inman-Bamber et al. 2008 Australian Journal of Agricultural Research, Fig.3
				agefactor = MIN(1.,agefactor)
				dper = dper * agefactor  !mm/h
				perday = perday + dper   !mm/day
	
		END DO
		
		IF (nstk .LT. 9) THEN
			diam_stk = -.077 * nstk + 3.0443
		ELSE
			diam_stk = -.0256 * nstk**2 + .4206 *  nstk + .7763
		END IF
	
		!DIAM_STK is a factor base on number of stalk per square meter used to calculate stem diameter       
		!DIAM is the average stem diameter in cm
	
		diam = .17409 + .1803853 * LOG(pleng*100)+.542597*diam_stk
		
		diam = (MAX(0.0, diam))
	
		nodencurrent = (lntotal)
		nodencurrent = INT(nodencurrent)
	
		IF (noden .LT. 1) THEN
			noden = 1
		END IF
		
		IF (nodencurrent .LT. noden) THEN       
			nodencurrent = noden - 1     
		END IF
	
			!INTERNODE ARRAY:
			! 1 - internode #
			! 2 - average plant diameter (cm) - the numbers recorded in each internode doesn't means the diameter of each internode - to be fixed
			! 3 - Sucrose Mass (kg/m2) 
			! 4 - Internode lenght (m)
	
		IF (noden == nodencurrent) THEN
				
				internode(noden,4) = internode(noden,4)  + perday/1000  !INTERNODE(NODEN+1,5) = INTERNODE(NODEN+1,5) + PERDAY/1000
			
		ELSE  
				
				internode(noden ,1) = INT(lntotal)               ! INTERNODE(NODEN+1,1) = INT(LNTOTAL)
				internode(noden ,2) = diam                       ! INTERNODE(NODEN+1,3) = DIAM
				
			!Calculating the sucrose mass in each internode
				DO i = 1, noden
					internode(i,3) = 0.       !Mass of Sucrose (kg/m2)
				END DO
	
				sucpf = 0.
	
				DO i = 1,noden
					sucpfacum = sucpfacum + ((sucmax/(1+EXP(-i+6))))
				END DO
			
				DO i =  1,noden
					sucpf(i) = ((sucmax/(1+EXP(-i+6))))
				END DO
	
				DO i = 1,noden
	
					internode(i,3) = sucpf(noden-i+1) * wsdm / noden     !Mass of Sucrose (kg/m2)
	
				END DO
				
			END IF            
		noden = nodencurrent 
		
	RETURN
      END SUBROUTINE DIAMPERS
	
	!************************************************************************

      
      

      ! ----------------------------------------------------------------------
      subroutine readsamuca (crpfil,pathcrop,tdwi,laiem,rgrlai,slatb,
     &  ssa,span,tbase,kdif,kdir,eff,amaxtb,tmpftb,tmnftb,cvl,cvr,cvs,
     &  q10,rml,rmr,rms,rfsetb,frtb,fltb,fstb,perdl,rdrrtb,
     &  rdrstb,hlim1,hlim2u,hlim2l,hlim3h,hlim3l,hlim4,rsc,adcrh,adcrl,
     &  cofab,rdi,rri,rdc,rdctb,rlwtb,logf,schedule,cumdens,
     &  flsolute,ecmax,ecslop,c2eca,c2ecb,c2ecf,numlay,dateharvest,
     &  swharvest,dmharvest1,dmharvest2,swgrazing,grazingfactor,
     &  daylastharvest,dmlastharvest,wrtmax,                            ! Nwgrassland
     &  nsuptab,dmfac,relnitab,nsupply,                                 ! Nwgrassland
     &  swcf,swetr,cftb,chtb,cfet,alphacrit,
     &  swroottyp,wiltpoint,rootradius,rootcoefa,rsw)                   ! NwRootExtr

! ----------------------------------------------------------------------
!     date               : november 2004   
!     purpose            : read parameters for grass growth routine
! ----------------------------------------------------------------------
      implicit none
      include  'arrays.fi'

      integer crp,i,logf,ifnd,getun2,schedule,numlay,swroottyp
      integer swcf,swetr
      logical flsolute,rdinqr
      real*8 slatb(30),amaxtb(30),tmpftb(30),depth,rootdis(202)
      real*8 tmnftb(30),rfsetb(30),frtb(30),fltb(30),rdrrtb(30)
      real*8 rdrstb(30),kdif,kdir,laiem,cofab,fstb(30)
      real*8 hlim1,hlim2u,hlim2l,hlim3h,hlim3l,hlim4,rdctb(22),rlwtb(22)
      real*8 sum,adcrh,adcrl,afgen,cvl,cvr,cvs,eff,rmr
      real*8 perdl,q10,rdc,rdi,rgrlai,rml,rms,rri,span,ssa,tbase,tdwi
      real*8 cumdens(202)
      real*8 ecmax,ecslop,c2eca,c2ecb,c2ecf(maho)
      real*8 cftb(2*magrs),chtb(2*magrs)
      real*8 albedo,rsc,rsw,cfet,alphacrit
      real*8 wiltpoint,rootradius,rootcoefa                           ! NwRootExtr 
      integer daylastharvest,swharvest,swgrazing                      ! Nwgrassland
      real*8 dmharvest1,dmharvest2, dmlastharvest,dateharvest(999)    ! Nwgrassland
      real*8 nsuptab(magrs),dmfac(magrs),relnitab(2*magrs),nsupply    ! Nwgrassland
      real*8 wrtmax,grazingfactor                                     ! Nwgrassland
      character crpfil*(*),pathcrop*(*)
! locals
      integer   swc2ecf
      real*8    dnrinput(magrs),cfinput(magrs),chinput(magrs)
      character message*200,filnam*200
! ----------------------------------------------------------------------


! --- initialise and start reading
      filnam = trim(pathcrop)//trim(crpfil)//'.crp'
      crp = getun2 (10,90,2)
      call rdinit(crp,logf,filnam)

! --- ET related params  ---------

! --- crop factor or crop height
      call rdsinr ('swcf',1,2,swcf)

! --- check use of crop factors in case of ETref
      if (swetr.eq.1 .and. swcf.eq.2) then
        message = 'If ETref is used (SWETR = 1), always define crop '//
     &           'factors (SWCF = 1)' 
        call fatalerr ('ReadGrass',message)
      endif

      if (swcf.eq.1) then
! ---   crop factor is input
        call rdador ('dnr',0.0d0,366.0d0,dnrinput,(magrs),ifnd)
        call rdfdor ('cf',0.0d0,2.0d0,cfinput,(magrs),ifnd)
! ---   store values in cftb
        do i = 1,ifnd
          cftb(i*2) = cfinput(i) 
          cftb(i*2-1) = dnrinput(i)
        enddo
        chtb = -99.99d0
      else
! ---   crop height is input
        call rdador ('dnr',0.0d0,366.0d0,dnrinput,(magrs),ifnd)
        call rdfdor ('ch',0.0d0,1.0d4,chinput,(magrs),ifnd)
! ---   store values in chtb
        do i = 1,ifnd
          chtb(i*2) = chinput(i) 
          chtb(i*2-1) = dnrinput(i)
        enddo
        cftb = -99.99d0
      endif

! --- reflection coefficient and crop resistance
      if (swcf.eq.1) then
! ---   use standard values for ETref
        albedo = 0.23d0
        rsc = 70.0d0
        rsw = 0.0d0
      else
! ---   use crop specific values
        call rdsdor ('albedo',0.0d0,1.0d0,albedo)
        call rdsdor ('rsc',0.0d0,1.0d6,rsc)
        call rdsdor ('rsw',0.0d0,1.0d6,rsw)
      endif


! --- crop growth related params ---------

! --- initial
      call rdsdor ('tdwi',0.0d0,10000.0d0,tdwi)
      call rdsdor ('laiem',0.0d0,10.0d0,laiem)
      call rdsdor ('rgrlai',0.0d0,1.0d0,rgrlai)

! --- green area
      call rdador ('slatb',0.0d0,366.0d0,slatb,30,ifnd)
      call rdsdor ('ssa',0.0d0,1.0d0,ssa)
      call rdsdor ('span',0.0d0,366.0d0,span)
      call rdsdor ('tbase',-10.0d0,30.0d0,tbase)

! --- assimilation
      call rdsdor ('kdif',0.0d0,2.0d0,kdif)
      call rdsdor ('kdir',0.0d0,2.0d0,kdir)
      call rdsdor ('eff',0.0d0,10.0d0,eff)
      call rdador ('amaxtb',0.0d0,366.0d0,amaxtb,30,ifnd)
      call rdador ('tmpftb',-10.0d0,50.0d0,tmpftb,30,ifnd)
      call rdador ('tmnftb',-10.0d0,50.0d0,tmnftb,30,ifnd)

! --- conversion of assimilates into biomass
      call rdsdor ('cvl',0.0d0,1.0d0,cvl)
      call rdsdor ('cvr',0.0d0,1.0d0,cvr)
      call rdsdor ('cvs',0.0d0,1.0d0,cvs)

! --- maintenance respiration
      call rdsdor ('q10',0.0d0,5.0d0,q10)
      call rdsdor ('rml',0.0d0,1.0d0,rml)
      call rdsdor ('rmr',0.0d0,1.0d0,rmr)
      call rdsdor ('rms',0.0d0,1.0d0,rms)
      call rdador ('rfsetb',0.0d0,366.0d0,rfsetb,30,ifnd)

! --- partitioning
      call rdador ('frtb',0.0d0,366.0d0,frtb,30,ifnd)
      call rdador ('fltb',0.0d0,366.0d0,fltb,30,ifnd)
      call rdador ('fstb',0.0d0,366.0d0,fstb,30,ifnd)

! --- death rates
      call rdsdor ('perdl',0.0d0,3.0d0,perdl)
      call rdador ('rdrrtb',0.0d0,366.0d0,rdrrtb,30,ifnd)
      call rdador ('rdrstb',0.0d0,366.0d0,rdrstb,30,ifnd)

! --- water use
      swroottyp = 1
      if(rdinqr('swroottyp')) then
        call rdsinr ('swroottyp',1,2,swroottyp)                         ! NwRootExtr
      endif
      if (swroottyp.eq.1) then                                          ! NwRootExtr
        call rdsdor ('hlim1' ,-100.0d0,100.0d0,hlim1)                   !
        call rdsdor ('hlim2u',-1000.0d0,100.0d0,hlim2u)                 !
        call rdsdor ('hlim2l',-1000.0d0,100.0d0,hlim2l)                 !
        call rdsdor ('hlim3h',-10000.0d0,100.0d0,hlim3h)                !
        call rdsdor ('hlim3l',-10000.0d0,100.0d0,hlim3l)                !
        call rdsdor ('hlim4' ,-16000.0d0,100.0d0,hlim4)                 !
                                                                        !
        call rdsdor ('adcrh',0.0d0,5.0d0,adcrh)                         !
        call rdsdor ('adcrl',0.0d0,5.0d0,adcrl)                         !
!       Criticial stress index for compensation of root water uptake (-)
        alphacrit = 1.0d0
        if(rdinqr('alphacrit')) then
          call rdsdor ('alphacrit',0.2d0,1.0d0,alphacrit)
        endif

      else                                                              !
        call rdsdor ('wiltpoint',-1.0d6,-1.0d3,wiltpoint)               !
        call rdsdor ('rootradius',0.0001d0,1.0d0,rootradius)            !
        call rdsdor ('rootcoefa',0.0d0,1.0d0,rootcoefa)                 !
      endif                                                             ! NwRootExtr

!     correction factor to relate potential transpiration to reference crop
!     (default = 1.0, range = 0.8 - 1.2)
      cfet = 1.0d0
      if(rdinqr('cfet')) then
        call rdsdor ('cfet',0.5d0,1.5d0,cfet)
      endif


! --- salt stress
      if (flsolute) then
        call rdsdor ('ecmax', 0.0d0,20.0d0,ecmax)
        call rdsdor ('ecslop',0.0d0,40.0d0,ecslop)
        call rdsdor ('c2eca', 0.0d0,1000.0d0,c2eca)
        call rdsdor ('c2ecb', 0.0d0,10.0d0,c2ecb)
        call rdsinr ('swc2ecf',1,2,swc2ecf)
        if (swc2ecf.eq.1) then
          call rdsdor ('c2ecf', 0.d0, 10.d0, c2ecf(1))
          if(numlay.ge.2) then
            do i = 2,numlay
              c2ecf(i) = c2ecf(1)
            enddo
          endif
        else if (swc2ecf.eq.2) then
          call rdfdor ('c2ecf', 0.d0, 10.d0, c2ecf,maho,numlay)
        endif
      endif

! --- interception
      call rdsdor ('cofab',0.0d0,1.0d0,cofab)

! --- rooting
      call rdador ('rdctb',0.0d0,100.0d0,rdctb,22,ifnd)
      call rdador ('rlwtb',0.0d0,5000.0d0,rlwtb,22,ifnd)
      call rdsdor ('rdi',0.0d0,1000.0d0,rdi)
      call rdsdor ('rri',0.0d0,100.0d0,rri)
      call rdsdor ('rdc',0.0d0,1000.0d0,rdc)

! --- nutrient stress
      call rdador ('nsuptab',0.0d0,1000.0d0,nsuptab,magrs,ifnd)
      call rdfdor ('dmfac',0.0d0,1.0d0,dmfac,magrs,ifnd)
!     store values in relni
      do i = 1,ifnd
        relnitab(i*2) = dmfac(i)
        relnitab(i*2-1) = nsuptab(i)
      enddo
      call rdsdor ('Nsupply',0.0d0,1000.0d0,Nsupply)


! --- harvest: dm or dates
      call rdsinr ('swharvest',1,2,swharvest)
      if(swharvest.eq.1) then
        call rdsdor ('dmharvest1',0.0d0,100000.0d0,dmharvest1)
        call rdsinr ('swgrazing',1,2,swgrazing)
        if(swgrazing.eq.2) then
            call rdsdor ('dmharvest2',0.0d0,100000.0d0,dmharvest2)
        else if (swgrazing.eq.1) then
            call rdsdor ('grazingfactor',0.0d0,1.0d0,grazingfactor)
        end if
        call rdsinr ('daylastharvest',1,366,daylastharvest)
        call rdsdor ('dmlastharvest',0.0d0,100000.0d0,dmlastharvest)
      elseif(swharvest.eq.2) then
        call rdatim ('dateharvest',dateharvest,999,ifnd)
      endif

! --- maximum weight of roots (kg/ha dm)
      call rdsdor ('wrtmax',0.0d0,100000.0d0,wrtmax)

! --- determine whether irrigation scheduling is applied
      call rdsinr ('schedule',0,1,schedule)

      close (crp)

! --- CALCULATE NORMALIZED CUMULATIVE ROOT DENSITY FUNCTION
      
      if (swroottyp .eq. 1) then                                        ! NwRootExtr

! ---   specify array ROOTDIS with root density distribution
        do i = 0,100
          depth = 0.01d0 * dble(i)
          rootdis(i*2+1) = depth
          rootdis(i*2+2) = afgen(rdctb,22,depth)
        enddo

! ---   calculate cumulative root density function
        do i = 1,202,2
! ---     relative depths
          cumdens(i) = rootdis(i)
        enddo
        sum = 0.d0
        cumdens(2) = 0.d0
        do i = 4,202,2
! ---     cumulative root density
          sum = sum + (rootdis(i-2)+rootdis(i)) * 0.5d0
     &               * (cumdens(i-1)-cumdens(i-3))
          cumdens(i) = sum
        enddo

! ---   normalize cumulative root density function to one
        do i = 2,202,2
          cumdens(i) = cumdens(i) / sum
        enddo

      endif                                                             ! NwRootExtr


      return
      end