c     hgrie Aug 2020: v1.0 fewbody-Compton
c     hgrie Aug 2020: for usesymmetry.and.Mzp=0.andMz=0, Resultyx and Resultyx not calculated
c             They must be zero by symmetry, see manu-script "Compton Densities Approach" p.53
c     new Aug 2020, based on 3He density codes with the following datings/changes:
c     hgrie May 2018: produce onebody amplitudes from 1Ndensities.
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     TO do:
c     Does output of multiple angles & energies to same output file work? May need adjusting output file name.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CHANGES:
c
c     hgrie Aug/Sep 2020: rewrote makedensityfilename() to deal with extracting densities from a .gz or downloading from server
c     hgrie June 2018: renamed "parity" to "symmetry -- see notes in usesymmetry+*.f
c       
c     hgrie May 2018: version 1 based on traditional main.onebody.f
c      
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   
      PROGRAM onebodydensitymain
      
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
      
      real*8 FT_sPlus, FT_sMinus, FL_sPlus, FL_sMinus
      real*8 E_prot, E_neut, L_prot, L_neut, mPionPlus

      real*8 thetaLow,thetaHigh,thetaInterval
      integer calctype,frame,Nangles,Nenergy,ienergy,j ! number of energies/angles; index for energies/angles
      character*200 descriptors  ! additional descriptors of calculation

      character*3 nucleus ! name of nucleus to be considered, for output file name
      integer Anucl             ! target nucleus mass number
      integer twoSnucl          ! 2 x target nucleus spin
      real*8 Mnucl               ! mass of target nucleus
      
      character*500 outfile
      character*500 densityFileName,originaldensityFileName ! second for multiple energies or angles
      
c*********************************************************************
c     Quadrature variables: 
c     onebody knows only about 1N amplitude's Feynman parameter integration
c     
      integer Nx                ! grid size 
      real*8 xq(Nxmax),wx(Nxmax)! points & weights
      
c     
c----------------------------------------------------------------------
c     
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
c     2* projections of single-nucleon's isospin & spin, both in & out
c     NB: this is the nucleon struck by the photons
      integer twomt1N,twomt1Np,twom1N,twom1Np
      
      integer rindx
      integer i
      integer a 
      integer aa
c     projections of target nucleus' in-spin, out-spin
      integer twoMz,twoMzp

c     hgrie June 2014: added variable twoMzplimit; changed by flag "nosymmetry" in input file.
c     Value twoMzplimit = 0 calculates half of the amplitudes, the other amps then from symmetry
c     Value twoMzplimit = -twoSnucl calculates all amplitudes
      integer twoMzplimit
      
c     integer twoMzlimit ! for symmetry calculation: Mzp>=0 *and* for Mzp=0, only Mz>=0, else Mz between +Snucl and -Snucl
      real*8 frac

c     complex*16, allocatable :: Resultxx(:,:),Resultxy(:,:) ! twoMz from -twoSnucl to twoSnucl, stepsize 2; rest blank.
c     complex*16, allocatable :: Resultyx(:,:),Resultyy(:,:) ! twoMz from -twoSnucl to twoSnucl, stepsize 2; rest blank.
      complex*16, allocatable :: SResultx(:,:),SResulty(:,:), SResultz(:,:) ! twoMz from -twoSnucl to
      complex*16, allocatable :: VResultx(:,:),VResulty(:,:), VResultz(:,:) ! twoMz from -twoSnucl to
      complex*16, allocatable :: LResultPlus(:,:), LResultMinus(:,:)
c     complex*16 :: nucS(3)
c     That means arrays are less than 2^2=4 times bigger than need be, but that's ok since quite small anyway. 

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c     added by hgrie May 2018: arrays which hold the 1N amplitude in the basis (twomt1N,twom1Np,twom1N,L1N,ML1N)
      
c     NB: these are the quantum numbers of the nucleon struck by the photons
c         (in traditional approach, index "3" was used for spectator) -- see also docu to read1Ndensity()
c     oneNspinbasisXY, where
c     XY: photon helicities in Cartesian basis X: out; Y: in
c     index 1: twomt1N     conserved isospin of nucleon (1: proton, -1: neutron)
c     index 2: twom1Np     1N out-spin                  (1: up, -1: down)
c     index 3: twom1N      1N in-spin                   (1: up, -1: down)
c     index 4: L1N         angular momentum of 1N op.   (integer starts at 0, up to L1Nmax; "K" in Andreas' notes)
c     index 5: ML1N        mag. quantum of L1Nop        (-L1N to L1N: 2L1N+1 integers; "κ" in Andreas' notes)
      
c     At the moment, L1N & ML1N are meaningless (L=ML=0), but they are implemented here as stump already. 

      integer,parameter :: L1Nmax=0    
      integer L1N, ML1N
c     complex*8 :: epsL_DOTSigma(-1:1,-1:1)  
      complex*16 :: tmpPlus(-1:1, -1:1)
      complex*16 tmpMinus(-1:1, -1:1)
c     integer*8 :: test2(-1:1,-1:1)  
c     complex*16,allocatable :: oneNspinbasisxx(:,:,:,:,:) ! (twomt1N,twom1Np,twom1N,L1N,ML1N)
c     complex*16,allocatable :: oneNspinbasisxy(:,:,:,:,:) ! (twomt1N,twom1Np,twom1N,L1N,ML1N)
c     complex*16,allocatable :: oneNspinbasisyx(:,:,:,:,:) ! (twomt1N,twom1Np,twom1N,L1N,ML1N)
c     complex*16,allocatable :: oneNspinbasisyy(:,:,:,:,:) ! (twomt1N,twom1Np,twom1N,L1N,ML1N)
      
c     1N amps of proton or neutron outside fewbody loops: array with 1: proton (twomt1N=+1); -1: neutron (twomt1N=-1)
c     real*8 A1(-1:1),A2(-1:1),A3(-1:1),A4(-1:1),A5(-1:1),A6(-1:1)
c     real*8 A1p,A2p,A3p,A4p,A5p,A6p
c     real*8 A1n,A2n,A3n,A4n,A5n,A6n
      
      integer variedA           !BS: integer variedA to indicate which A is varied by calctype=VaryA
      logical cartesian         !hgrie Oct 2014: for output in Cartesian basis of photon polarisations
      
      integer verbosity         !verbosity index for stdout hgrie June 2014
      
      integer test              ! a generic integer for testing i/o
      logical testtf            ! a generic logical for testing i/o
c     if density file generated from a .gz, delete that temporary file after each energy/angle
c     if downloaded and .gz, also delete the download.
c     0: do not delete; 1: delete un-gz'd file; 2: delete downloaded and un-gz'd file 
      integer rmDensityFileLater  

      real*8 dummy

c     for calculating magnetic-moment insertions on the way: (twoMzp,twoMzp,twomt1N)
c     real*8,allocatable     :: FSinsertion0(:,:,:),FSinsertionz(:,:,:),FSinsertiony(:,:,:),FSinsertionx(:,:,:)
c     real*8,allocatable     :: FVinsertion0(:,:,:),FVinsertionz(:,:,:),FVinsertiony(:,:,:),FVinsertionx(:,:,:)
c     complex*8,allocatable     :: FS_L(:,:,:,:)
c     complex*8  :: FS_T(-1:1, -1:1)
c     complex*8  :: FV_T(-1:1, -1:1)
c     real*8  :: FPlus(-1:1, -1:1)
c     real*8  :: FMinus(-1:1, -1:1)
c     complex F1
c     real*8 FST
c     complex*16 FST1
c     complex*16 FST2
c     complex*16 FST3
c     complex*16 FST4
c     complex*16 FST5
c     complex*16 FST6
c     real*8 FVT
c     complex*16 FVT1
c     complex*16 FVT2
c     complex*16 FVT3
c     complex*16 FVT4
c     complex*16 FVT5
c     complex*16 FVT6
      real*8 Paper3HeFTSPlusV
      real*8 Paper3HeFTSMinusV
      integer mOne
      integer mTwo
      complex F2
      real*8 :: Iden(-1:1,-1:1)  ! (ms3p,ms3): sigma-0=unit matrix
      complex*16 :: sigmax(-1:1,-1:1)  ! (ms3p,ms3): sigma-x
c     real*8 :: isigmay(-1:1,-1:1) ! (ms3p,ms3): I times sigma-y !!!!
      complex*16 :: sigmay(-1:1,-1:1) ! (ms3p,ms3): sigma-y
      complex*16 :: sigmaz(-1:1,-1:1)  ! (ms3p,ms3): sigma-z
      real*8,parameter ::  munucleon(-1:1) = (/kappan,0.d0,kappap+1.d0/)! indices of entries: (-1,0,+1)!!!!
c     real*8 :: eps(3)   ! polarization vector
c     real*8 :: epsMat(3,3)   ! polarization vector
c     real*8 :: epsL(3)  ! eps_L as defined in paper
c     real*8 :: epsT(3)  ! eps_T as defined in paper
c     complex*8 :: epsTemp(3) 
      real*8 :: kVec(3)  ! Momentum vector of k
      real*8 :: kHat(3)  ! normalized kVec
      real*8 temp
      integer mOneArr(2)
      integer mTwoArr(2)
      integer ii
      integer jj
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     end OF VARIABLE DECLARATIONS, BEGINNING OF CODING
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Reading the input file from command line
c**********************************************************************
c     get the number of arguments
      narg=command_argument_count()
      inputfile='input.dat'
c     if you have 1 argument, write it to inputfile, otherwise stop 
      ! if (narg.eq.1) then
      !    call get_command_argument(1, inputfile)
      ! else
      !    write(*,*) "*** ERROR: Pass one input file as argument!"
      !    stop
      ! end if
c     
c     
c**********************************************************************
c     Reading in data from the input file
c**********************************************************************
      inUnitno=13
      outUnitno=10
      verbosity=0
      open(unit=inUnitno, file=inputfile, status= 'OLD',iostat=test)
      if (test .ne. 0) stop "*** ERROR: Could not open input file!!! Aborting."
      call ReadinputCommon(Elow,Ehigh,Einterval,frame,
     &     thetaLow,thetaHigh,thetaInterval,
     &     outfile,descriptors,densityFileName,inUnitno,
     &     nucleus,Anucl,twoSnucl,Mnucl,
     &     twoMzplimit,cartesian,verbosity)
      write(*,*) "Remember to change the value of the energy to match the density "
c      
      call ReadinputOnebody(inUnitno,calctype,variedA,descriptors,
c---- Variable to control Feynman quadrature settings------------------------
     &     Nx,
     &     verbosity)

      call ReadinputCommonComments(descriptors,inUnitno,verbosity)
      
      close(unit=inUnitno,iostat=test)
      if (test .ne. 0) stop "*** ERROR: Could not close input file!!! Aborting."
c
c     call makeoutputfilename(outfile,calctype,nucleus,descriptors,densityFileName,variedA,
c    &     Elow,Ehigh,Einterval,thetaLow,thetaHigh,thetaInterval,verbosity)
c**********************************************************************
c     FINAL PRELIMINARIES
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
c     hgrie June 2017: keep original filename: needed for replacements of energy & angle later 
      originaldensityFileName = densityFileName

      thetaL=0.d0
      thetacm=0.d0
c**********************************************************************
c     Setting up quadratures for the Feynman integrals
      call AnglePtsWts(Nx,1,Nxmax,0.d0,1.0d0,xq,wx,Nx,verbosity)
c**********************************************************************
c     BS: new Nangles algebra due to changed input form
      Nangles=int((thetaHigh-thetaLow)/thetaInterval)+1
      Nenergy=int((Ehigh-Elow)/Einterval)+1
     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Nov 30 2022 Hard coded values from https://arxiv.org/abs/1103.3400v2 for 3He
c     3He values

          mPionPlus = 139.57

          FT_sPlus=0.017
          FT_sMinus=1.480
          FL_sPlus=-0.079
          FL_sMinus=1.479
          
          E_prot = -1.16E-3/mPionPlus
          E_neut = 2.13E-3/mPionPlus
          L_prot = -1.35E-3/mPionPlus
          L_neut = -2.41E-3/mPionPlus
c       Epsilon stuff. WLOG photon comes from the z direction
c**********************************************************************
      open(unit=outUnitno, file=outfile,iostat=test)
      if (test .ne. 0) stop "*** ERROR: Could not open output file!!! Aborting."
c**********************************************************************
c     Loop over Energies
c**********************************************************************
      do j=1,Nenergy

         Egamma=Elow+Einterval*(j-1)
         ienergy=int(Egamma)

          kVec = (/0.0,0.0, real(Egamma)/)
          kHat = kVec/SQRT(DOT_PRODUCT(kVec, kVec))
          write(*,*) "kVec=", kVec
          write(*,*) "kHat=", kHat
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
c**********************************************************************
            call calcphotonmomenta(k,kth,kphi,t,kp,kpth,kpphi,omega,
     &           Qk,Qkth,Qkphi,kgamma,thetacm,verbosity)
            
c**********************************************************************
c      be a good boy and initialise everything to 0, overwriting entries from previous ω/θ
c**********************************************************************
c           allocate (oneNspinbasisxx(-1:1,-1:1,-1:1,0:L1Nmax,-L1Nmax:L1Nmax))
c           allocate (oneNspinbasisxy(-1:1,-1:1,-1:1,0:L1Nmax,-L1Nmax:L1Nmax))
c           allocate (oneNspinbasisyx(-1:1,-1:1,-1:1,0:L1Nmax,-L1Nmax:L1Nmax))
c           allocate (oneNspinbasisyy(-1:1,-1:1,-1:1,0:L1Nmax,-L1Nmax:L1Nmax))
c           oneNspinbasisxx = c0
c           oneNspinbasisxy = c0
c           oneNspinbasisyx = c0
c           oneNspinbasisyy = c0

            allocate(VResultx(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl))
            allocate(VResulty(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl))
            allocate(VResultz(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl))
            VResultx=c0
            VResulty=c0
            VResultz=c0

            allocate(LResultPlus(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl))
            allocate(LResultMinus(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl))
            allocate(SResultx(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl))
            allocate(SResultz(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl))
            allocate(SResulty(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl))
            LResultPlus=c0
            LResultMinus=c0
            SResultx=c0
            SResultz=c0
            SResulty=c0
c     for calculating electric FF on the way
            dummy=0.d0
c     for calculating magnetic moment insertions on the way

c           allocate(FS_L(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl,-1:1))
c           allocate(FSinsertion0(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl,-1:1))
c           allocate(FSinsertionx(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl,-1:1))
c           allocate(FSinsertiony(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl,-1:1))
c           allocate(FSinsertionz(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl,-1:1))
c           FSinsertion0 =0.0
c           FSinsertionx =0.0
c           FSinsertiony =0.0
c           FSinsertionz =0.0


c           allocate(FVinsertion0(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl,-1:1))
c           allocate(FVinsertionx(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl,-1:1))
c           allocate(FVinsertiony(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl,-1:1))
c           allocate(FVinsertionz(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl,-1:1))
c           FVinsertion0 =0.0
c           FVinsertionx =0.0
c           FVinsertiony =0.0
c           FVinsertionz =0.0

            sigmax=0.d0
            sigmay=0.d0
            sigmaz=0.d0
c     
            sigmax(1,-1)=dcmplx(1.d0,0)
            sigmax(-1,1)=dcmplx(1.d0,0)
c
            sigmay(1,-1)=dcmplx(0, -1.d0)
            sigmay(-1,1)=dcmplx (0, 1.d0) 
c
            sigmaz(1,1)=dcmplx(1.d0,0)
            sigmaz(-1,-1)=dcmplx(-1.d0,0)

            Iden=0.d0
            Iden(-1,-1)=1
            Iden(1,1)=1
c**********************************************************************
c     hgrie June 2017: create name of 1Ndensity file for given energy and angle, unpack it
c     define correct formats for energy and angle
c     hgrie May 2018: outsourced into subroutine common-densities/makedensityfilename.f
            densityFileName = originaldensityFileName
            call makedensityfilename(densityFileName,Egamma,thetacm,rmDensityFileLater,verbosity)
c**********************************************************************
c     hgrie May 2018: read 1N density
            call read1Ndensity(densityFileName,Anucl,twoSnucl,omega,thetacm,verbosity)
c           write(*,*) "#########################################################################################"
c           write(*,*) "#########################################################################################"
c           eps = (/1,0,0/)
c           epsT = eps - DOT_PRODUCT(eps, kHat)*kHat
c           epsL = DOT_PRODUCT(eps,kHat)*kHat
c           epsTemp=CMPLX(epsT, KIND=8)! cast to complex
            tmpPlus = Iden+sigmaz
            tmpMinus = Iden-sigmaz
c           write(*,*) "tmpPlus(-1,-1)=",Int(tmpPlus(-1,-1))
c           write(*,*) "tmpMinus(1,1)=",Int(tmpMinus(1,1))

c           LResultPlus=c0
c           LResultMinus=c0
c           do rindx=1,maxrho1bindex
c                   CALL get1Nqnnum(rindx,twom1N,twomt1N,twoMz,twom1Np,twomt1Np,twoMzp,L1N,ML1N)
c                   If ((L1N.eq.0).and.((twomt1N.eq.1).and.(twomt1Np.eq.1))) then
c                       LResultPlus(twoMzp,twoMz)=LResultPlus(twoMzp,twoMz)+Anucl*rho1b(rindx)*sigmax(twom1Np,twom1N)
c    &                          *tmpPlus(twomt1Np,twomt1N)/2
c                       LResultMinus(twoMzp,twoMz)=LResultMinus(twoMzp,twoMz)+Anucl*rho1b(rindx)*sigmax(twom1Np,twom1N)
c    &                           *tmpMinus(twomt1Np,twomt1N)/2
c                   end if ! L1N
c           end do              !rindx   
c           LResultMinus= 2* LResultMinus! S has factor of 1/2 with it
c           LResultPlus = 2* LResultPlus
c           write(*,*) ""
c           write(*,*) "My results for sigmax (twomt1N.eq.1).and.(twomt1Np.eq.1))"
c           write(*,*) "##############################################################################"
c           write(*,*) Real(LResultPlus)
c           write(*,*) Real(LResultMinus)
c           write(*,*) "(twoMzp,twoMz)    =  |  (1, -1)                  |   (-1,1)"
c           write(*,*) "---------------------|---------------------------|----------------------------"
c           write(*,*) "Real(LResultPlus) =  |", Real(LResultPlus(1,-1)),"|",Real(LResultPlus(-1,1))
c           write(*,*) "Real(LResultMinus)=  |", Real(LResultMinus(1,-1)),"|", Real(LResultMinus(-1,1))
c           write(*,*) ""
c           write(*,*) ""
c           LResultPlus=c0
c           LResultMinus=c0
c           do rindx=1,maxrho1bindex
c                   CALL get1Nqnnum(rindx,twom1N,twomt1N,twoMz,twom1Np,twomt1Np,twoMzp,L1N,ML1N)
c                   If ((L1N.eq.0).and.((twomt1N.eq.1).and.(twomt1Np.eq.1))) then
c                       LResultPlus(twoMzp,twoMz)=LResultPlus(twoMzp,twoMz)+Anucl*rho1b(rindx)*sigmax(twom1Np,twom1N)
c    &                          *tmpPlus(twomt1Np,twomt1N)/2
c                       LResultMinus(twoMzp,twoMz)=LResultMinus(twoMzp,twoMz)+Anucl*rho1b(rindx)*sigmax(twom1Np,twom1N)
c    &                           *tmpMinus(twomt1Np,twomt1N)/2
c                   end if ! L1N
c           end do              !rindx   
c           LResultMinus= 2* LResultMinus! S has factor of 1/2 with it
c           LResultPlus = 2* LResultPlus
c           write(*,*) "My results for sigmax (twomt1N.eq.-1).and.(twomt1Np.eq.-1)"
c           write(*,*) "##############################################################################"
c           write(*,*) "(twoMzp,twoMz)    =  |  (1, -1)                  |   (-1,1)"
c           write(*,*) "---------------------|---------------------------|----------------------------"
c           write(*,*) "Real(LResultPlus) =  |", Real(LResultPlus(1,-1)),"|",Real(LResultPlus(-1,1))
c           write(*,*) "Real(LResultMinus)=  |", Real(LResultMinus(1,-1)),"|", Real(LResultMinus(-1,1))
c           write(*,*) ""
c           write(*,*) ""
            LResultPlus=c0
            LResultMinus=c0
            do rindx=1,maxrho1bindex
                    CALL get1Nqnnum(rindx,twom1N,twomt1N,twoMz,twom1Np,twomt1Np,twoMzp,L1N,ML1N)
                    If (L1N.eq.0) then
                        LResultPlus(twoMzp,twoMz)=LResultPlus(twoMzp,twoMz)+Anucl*rho1b(rindx)*sigmax(twom1Np,twom1N)
     &                          *tmpPlus(twomt1Np,twomt1N)/2
                        LResultMinus(twoMzp,twoMz)=LResultMinus(twoMzp,twoMz)+Anucl*rho1b(rindx)*sigmax(twom1Np,twom1N)
     &                           *tmpMinus(twomt1Np,twomt1N)/2
                    end if ! L1N
            end do              !rindx   
            LResultMinus= 2* LResultMinus! S has factor of 1/2 with it
            LResultPlus = 2* LResultPlus
            write(*,*) "My results for sigmax"
            write(*,*) "------------------------------------------------------------------------------"
c           write(*,*) Real(LResultPlus)
c           write(*,*) Real(LResultMinus)
            write(*,*) "(twoMzp,twoMz)    =  |  (1, -1)                  |   (-1,1)"
            write(*,*) "---------------------|---------------------------|----------------------------"
            write(*,*) "Real(LResultPlus) =  |", Real(LResultPlus(1,-1)),"|",Real(LResultPlus(-1,1))
            write(*,*) "Real(LResultMinus)=  |", Real(LResultMinus(1,-1)),"|", Real(LResultMinus(-1,1))
            write(*,*) ""
            write(*,*) ""
            LResultPlus=c0
            LResultMinus=c0
            do rindx=1,maxrho1bindex
                    CALL get1Nqnnum(rindx,twom1N,twomt1N,twoMz,twom1Np,twomt1Np,twoMzp,L1N,ML1N)
                    If (L1N.eq.0) then
                        LResultPlus(twoMzp,twoMz)=LResultPlus(twoMzp,twoMz)+Anucl*rho1b(rindx)*sigmay(twom1Np,twom1N)
     &                          *tmpPlus(twomt1Np,twomt1N)/2
                        LResultMinus(twoMzp,twoMz)=LResultMinus(twoMzp,twoMz)+Anucl*rho1b(rindx)*sigmay(twom1Np,twom1N)
     &                           *tmpMinus(twomt1Np,twomt1N)/2
                    end if ! L1N
            end do              !rindx   
            LResultMinus= 2* LResultMinus! S has factor of 1/2 with it
            LResultPlus = 2* LResultPlus
            LResultMinus = LResultMinus/sigmay
            LResultPlus = LResultPlus/sigmay
            write(*,*) ""
            write(*,*) "My results for sigmay"
            write(*,*) "-----------------------------------------------------------------------------"
c           write(*,*) "Imag(LResultPlus)=", AImag(LResultPlus(1,-1)),AImag(LResultPlus(-1,1))
c           write(*,*) "Imag(LResultMinus)=",AImag(LResultMinus(1,-1)),AImag(LResultMinus(-1,1))
            write(*,*) "LResultPlus=", Real(LResultPlus(1,-1)),Real(LResultPlus(-1,1))
            write(*,*) "LResultMinus=",Real(LResultMinus(1,-1)),Real(LResultMinus(-1,1))
            write(*,*) ""

            LResultPlus=c0
            LResultMinus=c0
            do rindx=1,maxrho1bindex
                    CALL get1Nqnnum(rindx,twom1N,twomt1N,twoMz,twom1Np,twomt1Np,twoMzp,L1N,ML1N)
c                   If ((L1N.eq.0).and.(twoMzp.eq.twoMz))then
                    If (L1N.eq.0) then
                        LResultPlus(twoMzp,twoMz)=LResultPlus(twoMzp,twoMz)+Anucl*rho1b(rindx)*sigmaz(twom1Np,twom1N)
     &                          *tmpPlus(twomt1Np,twomt1N)/2
                        LResultMinus(twoMzp,twoMz)=LResultMinus(twoMzp,twoMz)+Anucl*rho1b(rindx)*sigmaz(twom1Np,twom1N)
     &                           *tmpMinus(twomt1Np,twomt1N)/2
                    end if ! L1N
            end do              !rindx   
            LResultMinus= 2* LResultMinus! S has factor of 1/2 with it
            LResultPlus = 2* LResultPlus
            write(*,*) ""
            write(*,*) "My results for sigmaz "
            write(*,*) "-----------------------------------------------------------------------------"
            write(*,*) "Real(LResultPlus)=", Real(LResultPlus(-1,-1)/sigmaz(-1,-1)),Real(LResultPlus(1,1)/sigmaz(1,1))
            write(*,*) "Real(LResultMinus)=", Real(LResultMinus(-1,-1)/sigmaz(-1,-1)),Real(LResultMinus(1,1)/sigmaz(1,1))
            write(*,*) ""
            write(*,*) ""
            write(*,*) "####################################################################"
            write(*,*) "###################       End of My results      ###################"
            write(*,*) "####################################################################"
c           write(*,*) "Lenkewitz Thesis 3He result:"
c           write(*,*) "F^S+V=", -0.065, "and F^S-V=",1.801
c           write(*,*) ""
            write(*,*) "Lenkewitz 2011 paper 3He result"
            write(*,*) "F_T^S+V=0.017 and F_T^S-V=1.480"
            write(*,*) "F_L^S+V=-0.079 and F_L^S-V=1.479"
            write(*,*) ""
c           write(*,*) "Identity check", dcmplx(0,-1.d0)*MATMUL(MATMUL(sigmax, sigmay),sigmaz)
c           write(*,*) "Identity check", Real(MATMUL(sigmax, sigmax))
c           write(*,*) "Identity check", Real(MATMUL(sigmay, sigmay))
c           write(*,*) "Identity check", Real(MATMUL(sigmaz, sigmaz))
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     hgrie May 2018: symmetry and output
c            call outputroutine(outUnitno,cartesian,twoSnucl,twoMzplimit,
c     &           Resultxx,Resultxy,Resultyx,Resultyy,verbosity)
            
c     be a good boy and deallocate arrays. Compilers do that automatically for simple programs. Better safe than sorry.
            deallocate (LResultMinus, LResultPlus, SResultx,SResulty,VResulty, VResultx, STAT=test ) ! test becomes nonzero if this fails
c           if (test .ne. 0) stop "*** ERROR: Arrays ResulyAB: Deallocation error. Abort."
c          deallocate (oneNspinbasisxx,oneNspinbasisxy,oneNspinbasisyx,oneNspinbasisyy, STAT=test ) ! test nonzero if fails
c           if (test .ne. 0) stop "*** ERROR: Arrays insertion: Deallocation error. Abort."
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     hgrie Aug/Sep 2020: delete the local .dat file if one was generated from .gz
            if (rmDensityFileLater.gt.0) then
               call EXECUTE_COMMAND_LINE("rm "//densityFileName, WAIT=.True., EXITSTAT=test )
               if (test.ne.0) stop "*** ERROR: Could not remove .dat file created from .gz"
               write(*,*) "   Removed .dat file unzipped from .gz."
               if (rmDensityFileLater.ge.2) then
                  INQUIRE(FILE=TRIM(densityFileName)//".gz", EXIST=testtf)
                  if ( testtf ) then
                     call EXECUTE_COMMAND_LINE("rm "//TRIM(densityFileName)//".gz", WAIT=.True., EXITSTAT=test )
                     if (test.ne.0) stop "*** ERROR: Could not remove .dat file created from .gz"
                     write(*,*) "   Removed .dat.gz file downloaded."
                  end if
               end if
            end if
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         end do                 !Nangles
      end do                    !Nenergy
      close(outUnitno,iostat=test)
      if (test .ne. 0) stop "*** ERROR: Could not close output file!!!"
     
      write (*,*) '*** Wrote output to file: ',TRIM(outfile)
      
      write(*,*) "Remember to change the value of the energy to match the density "
      stop
      
c20   format(' ',A,I6,A,8I8,A,E24.15,SP,E25.15," I")
c30   format(' ',A,5I4,A,F20.13,SP,F21.13," I")
c40   format(A,2F18.13)
      
      end PROGRAM
