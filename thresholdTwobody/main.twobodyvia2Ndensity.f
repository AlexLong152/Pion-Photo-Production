c     hgrie Aug 2020: v1.0 fewbody-Compton
c     new Aug 2020, based on 3He density codes with the following datings/changes:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     hgrie May 2017, revised May 2018 (see below)
c     based on Andreas Nogga's "template" files common-densities/2Ndensity-module/testcompdens.F90 in May 2017/2018.
c                                           and common-densities/2Ndensity-module/CompDens.F90 in May 2017/2018.
c
c     This file adapted from main.twobody.f, plus changes.
c          added read of density matrix,
c          eliminated spectator(3)-integrations and calls, and calls to wave function
c          densityFileName now used to set name of input density file
c          split setquad setting up angular integrations into separate routines
c                 for (12) and spectator (3) system
c          here, (12) integration only -- spectator (3) integration provided by 2Ndensity
c     j12max can now be specified in input file: default is j12max=2 for onebody and j12max=1 for twobody,
c     as in previous runs (where they were hardwired in code).
c     These values give MEs which are converged to better than 0.7% -- see documentation/.
c
c     twoSmax/twoMz dependence: only via array size of ResultAB
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     TO do:
c     Does output of multiple angles & energies to same output file work? May need adjusting output file name.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CHANGES:
c
c     hgrie Sep 2020: in readinputTwobody():
c           removed NP12p = 3rd variable in integration-grid
c           line of input file. It is actually never used for anything!
c           That also reduced number of arguments of readinputTwobody()!
c     hgrie Aug/Sep 2020: rewrote makedensityfilename() to deal with extracting densities from a .gz or downloading from server
c     hgrie June 2018: renamed "parity" to "symmetry -- see notes in usesymmetry+*.f
c       
c     hgrie May 2018: decluttered files: remove obsolete variables, unify look, documentation,...
c     hgrie May 2018: new subroutines to read input, weeded out unused variables, new input.dat format. 
c     hgrie May 2018: rewritten to accommodate hdf5 format for input files of 2N density rho
c     hgrie May 2018:
c           All quantum numbers which start with "two" run over integers
c                  Examples: 
c                     twoMz,twoMzp: magnetic quantum numbers of in/out target nucleus, times 2.
c                     twoSnucl: 2 x spin of target nucleus
c           For all such variables, run over 2xQM values.
c                  Examples:
c                     in do-loops: twoMz runs from +twoSnucl to -twoSnucl, in steps of -2
c                     in array: Resultxx(twoMzp,twoMz) runs over   (-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)
c
c     This numbering agrees with Andreas' assignments of magnetic quantum numbers.
c     In ARRAYs, this means a slight waste of space since we leave many array entries unused.
c     Example S=1/2: 2x2=4 entries needed, but array is (-1:1,-1:1): 3x3=9  : 225% of needed size
c     Example S=1:   3x3=9 entries needed, but array is (-2:2,-2:2): 5x5=25 : 280% of needed size
c
c     Still, our arrays are for small nuclei -- we do not really waste a lot. Not a time/storage issue.
c
c     hgrie May 2018: outsourced symmetry+output into sub routine outputroutine(), identical for onebody and twobody
c      
c     Implemented symmetry for arbitrary nucleon spin:
c     Use Mzp>=0, and for Mzp=0, run only over Mz>=0
c     -- that's still 2 more than necessary since ME(+0->+0) = ME(-0->-0) and ME(+0->-0) = ME(-0->+0)
c     but it's good enough, saving lots of CPU time.
c     see manuscript "Compton Densities Approach" pp11-12
c      
c     In May 2018, Andreas replaced fkltt() by a more sophisticated and parallel routine initclebsch() in
c     common-densities/2Ndensity-module/clebsch.F .
c      
c     hgrie June 2017: implemented run over several energies and angles into same output file
c                      implemented: when input file contains a 2Ndensity filename which contains "XXX"
c                                   and "YYY" strings, then these are automatically replaced by
c                                   XXX => energy of run (in numeric format used by Andreas)
c                                   YYY => angle of run (in numeric format used by Andreas)
c                      That reduces error-proneness.
c      
c    modified by hgrie June 2014:
c     calculate nucleon amplitudes outside fewbody loops
c     use symmetry to calculate only amplitudes with 3He out-spin +1/2
c     modified hgrie 20 June 2014: add option to use LebedevLaikov or Gaussian
c     integration for theta & phi separately,
c     for solid angle integral in (12) system 
c     modifed by hgrie Oct 2014:
c     use spherical or cartesian basis, and symmetry or not in either.
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   

      PROGRAM twobodydensitymain
      
      USE CompDens              ! needs module CompDens.mod
      
      IMPLICIT NONE
c**********************************************************************
      include '../common-densities/params.def'
      include '../common-densities/calctype.def'
      include '../common-densities/constants.def'
c**********************************************************************
c     
c     program argument management
      integer narg              ! number of arguments
      character*500 inputfile   ! argument string is written to this
c     
c*********************************************************************
      integer NP12A,NP12B,NP12
      real*8  P12A,P12B,P12C
      real*8 P12MAG(Npmax),AP12MAG(Npmax)
c     
c**********************************************************************
c     
c     Input information:
c     
c     inUnitno-unit where input information is stored
c     outUnitno-unit where output is written
c     kgamma-photon momentum in cm frame in units of MeV
c     Nangles-number of angles at which calculation is to be done
c     thetaL-photon scattering angle in lab frame, in degrees
c     thetacm-photon scattering angle in cm frame, in radians
c     calctype-which calculation to do
c     
c     frame-which frame the results are to be presented for
c     1=lab. frame
c     2=c.m. frame
c     outfile-name of output file
c     
      integer inUnitno,outUnitno
      
      real*8 Egamma,kgamma,thetaL,thetacm,Elow,Ehigh,Einterval
      
      real*8 thetaLow,thetaHigh,thetaInterval
      integer calctype,frame,Nangles,Nenergy,ienergy,j ! number of energies/angles; index for energies/angles
      character*200 descriptors  ! additional descriptors of calculation

      integer,parameter :: variedA = 0 ! no variedA since no 1N amplitude, but define here for makeoutputfilename()
      character*3 nucleus ! name of nucleus to be considered, for output file name
      integer Anucl             ! target nucleus mass number
      integer twoSnucl          ! 2 x target nucleus spin
      real*8 Mnucl               ! mass of target nucleus
      character*500 outfile
      character*500 densityFileName,originaldensityFileName ! second for multiple energies or angles
      
c*********************************************************************
c     Quadrature variables:
c     
c     Nordth,Nthbins-number of quadratures/bin and number of bins for 
c     theta integration
c     Nordphi,Nphibins-number of quadratures/bin and number of bins for phi integration
c     Nth-total number of theta quadratures
c     Nphi12-total number of phi quadratures
c     pq,wp-radial quadratures and weights, set up on [0,infty]
c     thq,wth-theta quadratures and weights, set up on [0,PI]
c     phi12,wphi-phi quadratures and weights, set up on [0,2 PI]
c     
      integer Nordth12,Nordphi12,Nthbins12,Nphibins12

      integer Nth12,Nphi12
      
      real*8 th12(Nangmax),phi12(Nangmax)
      
      integer AngularType12,Nanggrid12
      real*8 angweight12(Nangmax,Nangmax)
c     
c----------------------------------------------------------------------
c     TODO: redo this     
c     Momentum variables:
c     
c     pp-magnitude of final-state relative three-momentum vector,
c     in IA=p + 1/2(k - k') as a vector sum. In IA the
c     kinematics are 
c     
c                 \    /
c                  \  /   
c     k' + k/2 + p  \/        -k/2 + p
c     ---------------------------------------
c     
c     
c     
c     ---------------------------------------
c     - k/2 - p
c     
      real*8 k,kth,kphi,kp,kpth,kpphi,Qk,Qkth,Qkphi
      real*8 t,omega
c     
c**********************************************************************
c     
      integer m12,mt12 ! projections of total ang mom & isospin of (12) subsystem: automatically integers
      
      integer i,ip12
      integer l12,s12,j12,t12 ! orb ang mom, spin, total ang mom, isospin of (12): automatically integers
      
      integer j12max            ! max total ang mom in (12) subsystem -- =1 suffices for 1% accuracy.
      
c     projections of target nucleus' in-spin, out-spin
      integer twoMz,twoMzp      !  -- these two not use right now
      
c     hgrie June 2014: added variable twoMzplimit; changed by flag "nosymmetry" in input file.
c     Value twoMzplimit = 0 calculates half of the amplitudes, the other amps then from symmetry
c     Value twoMzplimit = -twoSnucl calculates all amplitudes
      integer twoMzplimit,twoMzlimit ! latter only for verbose output

      real*8 frac

      complex*16, allocatable :: Resultx(:,:),Resulty(:,:) ! twoMz from -twoSnucl to twoSnucl, stepsize 2; rest blank.
c     That means arrays are less than 2^2=4 times bigger than need be, but that's ok since quite small anyway. 
      
      logical cartesian         ! hgrie Oct 2014: for output in Cartesian basis of photon polarisations
      
      integer verbosity         ! verbosity index for stdout hgrie June 2014

      integer test              ! a generic integer for testing i/o
      logical testtf            ! a generic logical for testing i/o
c     if density file generated from a .gz, delete that temporary file after each energy/angle
c     if downloaded and .gz, also delete the download.
c     0: do not delete; 1: delete un-gz'd file; 2: delete downloaded and un-gz'd file 
      integer rmDensityFileLater 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     end OF VARIABLE DECLARATIONS, BEGINNING OF CODING
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      write(*,*) 
      write(*,*) "================================================================================"
      write(*,*) "Twobody Contributions to Few-Nucleon Compton Scattering Calculated Via 2N-Density Matrix"
      write(*,*) "================================================================================"
      write(*,*) "   Version 1.0"
      write(*,*) "      D. Phillips/A. Nogga/hgrie starting August 2020   "
      write(*,*) "      based on 3He codes: D. Phillips/A. Nogga/hgrie starting May 2018"
      write(*,*) "                          D. Phillips/B. Strasberg/hgrie starting June 2014"
      write(*,*) "                          with A. Margaryan 2016-17, modifying codes by D. Shukla 2007/8"
      write(*,*)
c**********************************************************************
c     Reading the input file from command line
c**********************************************************************

c     get the number of arguments
      narg=command_argument_count()

c     if you have 1 argument, write it to inputfile, otherwise stop 
      if (narg.eq.1) then
         call get_command_argument(1, inputfile)
      else
         write(*,*) "*** ERROR: Pass one input file as argument!"
         stop
      end if
c     
c**********************************************************************
c     Reading in data from the input file
c**********************************************************************
      inUnitno=13
      outUnitno=10

      open(unit=inUnitno, file= inputfile, status= 'OLD',iostat=test)
      if (test .ne. 0) stop "*** ERROR: Could not open input file!!! Aborting."
      
      call ReadinputCommon(Elow,Ehigh,Einterval,frame,
     &     thetaLow,thetaHigh,thetaInterval,
     &     outfile,descriptors,densityFileName,inUnitno,
     &     nucleus,Anucl,twoSnucl,Mnucl,
     &     twoMzplimit,cartesian,verbosity)
c      
      call ReadinputTwobody(inUnitno,calctype,descriptors,
c---- Variables to control radial quadrature settings------------------------
     &     NP12A,NP12B,P12A,P12B,P12C,
c---- Variables to control angular quadrature settings------------------------
     &     AngularType12,Nanggrid12,
     &     Nordth12,Nordphi12,
     &     NthBins12,NphiBins12,
     &     j12max)
      
      call ReadinputCommonComments(descriptors,inUnitno,verbosity)
      
      close(unit=inUnitno,iostat=test)
      if (test .ne. 0) stop "*** ERROR: Could not close input file!!! Aborting."
c
      call makeoutputfilename(outfile,calctype,nucleus,descriptors,densityFileName,variedA,
     &     Elow,Ehigh,Einterval,thetaLow,thetaHigh,thetaInterval,verbosity)
      
c**********************************************************************
c     FINAL PRELIMINARIES
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c     hgrie June 2017: keep original filename: needed for replacements of energy & angle later 
      originaldensityFileName = densityFileName

      thetaL=0.d0
      thetacm=0.d0
c     BS: new Nangles algebra due to changed input form
      Nangles=int((thetaHigh-thetaLow)/thetaInterval)+1
      Nenergy=int((Ehigh-Elow)/Einterval)+1
      
c**********************************************************************
c     (12) integration set-up
c          spectator (3) integration is provided by 2Ndensity, so nothing to do. 
c     define total number of integration points for (12) mom magnitude
      NP12 = NP12A+NP12B
c
c     Set up radial quadratures for (12) integration
      call TRNS(NP12A,NP12B,NP12,P12A,P12B,P12C,P12MAG,AP12MAG)
      
c     Set up angular quadratures for (12) integration
      call Setquad12(th12,Nth12,phi12,Nphi12,
     &     Nordth12,Nthbins12,Nordphi12,Nphibins12,
     &     AngularType12,angweight12,Nanggrid12,1 !verbosity
     &     )

      write (*,*) "***************************** END OF INITIALISATION *****************************"
      
c**********************************************************************
      open(unit=outUnitno, file=outfile,iostat=test)
      if (test .ne. 0) stop "*** ERROR: Could not open output file!!! Aborting."
c**********************************************************************
c     Loop over Energies
c**********************************************************************
      do j=1,Nenergy
         Egamma=Elow+Einterval*(j-1)
         ienergy=int(Egamma)
         write(*,*) "cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc"
         write(*,*) "cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc"
         write(*,*)
         write(*,*) "Incoming photon energy (rounded) = ", ienergy, " MeV"
         write(*,*) "Number of Angles Nangles         = ", Nangles
c**********************************************************************
c     Loop over angles
c**********************************************************************
         do i=1,Nangles
            if (frame.eq.lab) then
               kgamma=Egamma/sqrt(1.d0 + 2.0d0*Egamma/Mnucl)
c     BS: new theta algebra due to changed input
c     hgrie Sep 2014: if thetaL is ZERO degrees, actual calculated at 1 Degree  
               thetaL=(thetaLow+real(i-1)*thetaInterval)*Pi/180.d0
               if (thetaL.eq.0.0d0) then
                  thetaL=1.0d0*Pi/180.d0
                  write(*,*) "   Replaced input angle 0 deg with 1 deg."
               end if   
               frac=(Mnucl + (Mnucl + Egamma)*(dcos(thetaL) - 1.0d0))/
     &              (Mnucl + Egamma*(1.d0 - dcos(thetaL)))
               thetacm=dacos(frac)
            else
               thetacm=(thetaLow+real(i-1)*thetaInterval)*Pi/180.d0
               if (thetacm.eq.0.0d0) then
                  thetacm=1.0d0*Pi/180.d0
                  write(*,*) "   Replaced input angle 0 deg with 1 deg."
               end if   
               kgamma=Egamma      
            end if
            if (frame.eq.cm) then
               write (outUnitno,*) "cm ","omega = ",Egamma,"thetacm = ",thetacm*180.0/Pi
            else if (frame.eq.lab) then
               write (outUnitno,*) "lab ","omega = ",Egamma,"thetalab = ",thetaL*180.0/Pi
            end if   
            write(*,*)
            write(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
            write(*,*) 'Calculating amps: theta =',thetacm*180.0/Pi,'deg: angle #',i
c**********************************************************************
            call calcphotonmomenta(k,kth,kphi,t,kp,kpth,kpphi,omega,
     &           Qk,Qkth,Qkphi,kgamma,thetacm,verbosity)
c           write(*,*) "In main"
c           write(*,*) "Egamma=", Egamma
c           write(*,*) "k=", k
c**********************************************************************
c     be a good boy and initialise everything to 0, overwriting entries from previous ω/θ
c**********************************************************************
            allocate(Resultx(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl))
            allocate(Resulty(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl))
c           allocate(Resultxx(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl))
c           allocate(Resultxy(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl))
c           allocate(Resultyx(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl))
c           allocate(Resultyy(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl))
            Resultx=c0
            Resulty=c0
c           Resultxx=c0
c           Resultxy=c0
c           Resultyx=c0
c           Resultyy=c0
c**********************************************************************
c     hgrie June 2017: create name of 1Ndensity file for given energy and angle, unpack it
c     define correct formats for energy and angle
c     hgrie May 2018: outsourced into subroutine common-densities/makedensityfilename.f
            densityFileName = originaldensityFileName
            call makedensityfilename(densityFileName,Egamma,thetacm,rmDensityFileLater,verbosity)
c**********************************************************************
c     hgrie May 2017: read 2Ndensity
            call read2Ndensity(densityFileName,Anucl,omega,thetacm,j12max,P12MAG,AP12MAG,NP12,verbosity)
c**********************************************************************      
c     hgrie Aug/Sep 2020: delete the local .h5 file if one was generated from .gz
            if (rmDensityFileLater.gt.0) then
               call EXECUTE_COMMAND_LINE("rm "//densityFileName, WAIT=.True., EXITSTAT=test )
               if (test.ne.0) stop "*** ERROR: Could not remove .h5 file created from .gz"
               write(*,*) "   Removed .h5 file unzipped from .gz."
               if (rmDensityFileLater.ge.2) then
                  INQUIRE(FILE=TRIM(densityFileName)//".gz", EXIST=testtf)
                  if ( testtf ) then
                     call EXECUTE_COMMAND_LINE("rm "//TRIM(densityFileName)//".gz", WAIT=.True., EXITSTAT=test )
                     if (test.ne.0) stop "*** ERROR: Could not remove .h5 file created from .gz"
                     write(*,*) "   Removed .h5.gz file downloaded."
                  end if
               end if
            end if
c**********************************************************************
            write(*,*) "*********Now convoluting 2N helicity amplitudes with 2N density matrix.*********"  
            do mt12=0,0                ! only charged pion exchange at OQ4 => (12) subsystem is (pn) 
               do j12=0,j12max                      ! total ang mom (12); usually Jmax=1 for 1% convergence
                  do s12=0,1                        ! spin (12)
                     do l12=abs(j12-s12),j12+s12    ! angular mom. (12)
                        t12=(1-(-1)**(l12+s12+1))/2 ! isospin (12)
                        do m12=-j12,j12             ! spin projection (12)
                           do ip12=1,NP12           ! integration over momentum magnitude (12)
                              call twobodyfinalstatesumsvia2Ndensity(
     &                             Resultx,Resulty,
     &                             Anucl,twoSnucl,twoMzplimit,j12,m12,l12,s12,t12,mt12,
     &                             k,thetacm,
     &                             ip12,P12MAG(ip12),AP12MAG(ip12),     
     &                             p12MAG,AP12MAG,NP12,                                 
     &                             th12,phi12,Nth12,Nphi12,j12max,
     &                             AngularType12,angweight12,calctype, Mnucl, verbosity)
                           end do !ip12
                        end do    !m12
                     end do       !l12
                  end do          !s12
               end do             !j12
            end do                !mt12
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc            
            if (verbosity.ge.3) then
               do twoMzp=twoSnucl,twoMzplimit,-2
c         for Mzp=0, run only over Mz>=0 -- that's still 2 more than necessary, but good enough -- see cure below
                  if ((twoMzp.eq.0).and.(twoMzplimit.eq.0)) then
                     twoMzlimit = 0
                  else
                     twoMzlimit = -twoSnucl
                  end if   
                  do twoMz=twoSnucl,twoMzlimit,-2
c                    write (*,*) "Resultx(twoMzp=",twoMzp,", twoMz=",twoMz,"): ",Resultx(twoMzp,twoMz)
c hgrie Aug 2020: now cure: for Mzp=Mz=0, only calculate xx and yy, since xy and yx must be zero               
                     if ((twoMzplimit.eq.0).and.(twoMzp.eq.0).and.(twoMz.eq.0)) then
                        continue
                     else
                        write (*,*) "Resultx(twoMzp=",twoMzp,", twoMz=",twoMz,"): ",Resultx(twoMzp,twoMz)
                        write (*,*) "Resulty(twoMzp=",twoMzp,", twoMz=",twoMz,"): ",Resulty(twoMzp,twoMz)
                     end if
c      end cure
c                    write (*,*) "Resultyy(twoMzp=",twoMzp,", twoMz=",twoMz,"): ",Resultyy(twoMzp,twoMz)
                  end do
               end do
            end if
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      found in usesymmetry+writeoutput-densities.f
            call outputPiPhoto(outUnitno,cartesian,twoSnucl,twoMzplimit,
     &           Resultx,Resulty,verbosity) 
            
c     be a good boy and deallocate arrays. Compilers do that automatically for simple programs. Better safe than sorry.
            deallocate (Resultx,Resulty, STAT=test ) ! test becomes nonzero if this fails
            if (test .ne. 0) stop "*** ERROR: Arrays ResultAB: Deallocation error. Abort."
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         end do                 !Nangles
         write(*,*) "cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc"
         write(*,*) "cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc"
      end do                    !Nenergy
      close(outUnitno,iostat=test)
      if (test .ne. 0) stop "*** ERROR: Could not close output file!!!"
      
      write (*,*) '*** Wrote output to file: ',TRIM(outfile)

      stop
      end PROGRAM

      
