
import os.path

from stat import S_IRWXG,S_IRWXU
from shutil import rmtree
import warnings

import numpy as np
from numpy import sqrt
from numpy import array as npa
from numpy import power,absolute,logical_and,column_stack,zeros,empty,append,\
sqrt,absolute,recarray


from math import pow,exp,log10,pi

try:
    from rootpy.plotting import Canvas,Hist,Hist2D,Graph
    from rootpy.plotting.style import set_style
    from rootpy.io import root_open


    warnings.simplefilter("ignore")
except:
    print "Could not load in root_numpy or rootpy, they are required to run this module."

defaultValues  = [3,2500,2805.,\
10026.35,10026.35,3080.0,6.35,1000.,\
0.043, 0.133,0.002,\
10.,0.2,0.25,0.28,0.35,1.7,\
'merged_ntuple_watchman','merged_ntuple_watchman','null', 'processed_watchman.root',\
10.,2.0, 100.0, 9, 0.65,0.1,\
'day','boulby', \
1.0]

docstring = """
    Usage: watchmakers.py [options]

    Arguments:

    Options:

    ## System options

    --docker            Is watchmakers running with docker, is so changes jobs files
    --singularity       Will change job submission script to reflect use of singularity
    --simg=<simg>       What singularity image will you be using [Default: /usr/gapps/adg/geant4/rat_pac_and_dependency/AITWATCHMAN/watch.simg]
    --newVers=<nv>      Major revision to Watchmakers. By default on, but can use old vers [Default: 1]
    --force             Forcing the recreation of the root_file,bonsai_root_file and log folders
    --noRoot            Allows to generate scripts without loading in any ROOT module
    -v                  Verbose. Allow print out of additional information.


    ## Create macros and job scripts for a user defined detector configuration

    -m                  Generate rat-pac macro files
    -j                  Create rat-pac/bonsai submision scripts for option above. Can be run with -m.

    -N=<N>                 Number of rat-pac macros per unique process [Default: %d]
    -e=<runBeamEntry>      Number of entries per macro (U/Th event x5) [Default: %d]
    --depth=<depthD>       Depth of detector (for fast neutron spectra) [Default: %f]
    --tankRadius=<TR>      Total radius of tank (mm) [Default: %f]
    --halfHeight=<HH>      Half height of tank (mm) [Default: %f]
    --shieldThick=<ST>     Steel->PMT distance (mm) [Default: %f]
    --vetoThickR=<VTR>     Steel->PMT radius distance (mm)
    --vetoThickZ=<VTZ>     Steel->PMT height distance (mm)
    --steelThick=<StT>     Steel Thickness (mm)     [Default: %f]
    --fidThick=<fT>        Fiducial volume-> PMT Thickness (mm) [Default: %f]
    --detectMedia=<_dM>    Detector media (doped_water,...)
    --collectionEff=<CE>   Collection efficiency (e.g.: 0.85,0.67,0.475)
    --pmtModel=<_PMTM>     PMT Model (r7081_lqe/r7081_hqe for 10inch or r11780_lqe/r11780_hqe for 12inch) 
    --photocath =<_PC>     PMT photocathode (R7081HQE)
    --pmtCtrPoint          Point inner PMTs of Watchman geometry to detector center
    -C=<cov>               Pick a single detector Configuration (25pct, SuperK,...). If not config, all
                           are generated
    --ipc=<_ipc>	   Inner PMTs photocoverage. If set WILL OVERIDE DEFAULT CONFIGURATION COVERAGE
    --vpc=<_vpc>           Veto PMTs photocoverage.
    --vetoModel=<_vPMTM>   Veto PMT Model (ETEL, r7081_lqe/r7081_hqe for 10inch or r11780_lqe/r11780_hqe for 12inch) 


    ## Analysis - efficiency and sensitivity evaluation. Once jobs have run, do these steps

    -M                  Merge result files from trial ntuples. Step one.
    --histograms        Apply n9-position cuts and create efficiency histograms. Step two.
    --PMTAccHists       Apply n9/good_pos/good_dir cuts and Look at Acceptance for WaterVolume events in PMT volume. Optional.
    --evalRate          Evaluate the rates and combine to efficiency histrograms Step three.

    --U238_PPM=<_Uppm>  Concentration of U-238 in glass [Default: %f]
    --Th232_PPM=<_Thp>  Concentration of Th-232 in glass [Default: %f]
    --K_PPM=<_K>        Concentration of K-40 in glass [Default: 16.0]
    --Rn222=<_Rn>       Radon activity in water SK 2x10^-3 Bq/m^3 [Default: %f]
    
    --U238_vPPM=<_vUppm>  Concentration of U-238 in veto pmt glass (ULB 0.0431 STD 0.341)
    --Th232_vPPM=<_vThp>  Concentration of Th-232 in veto pmt glass (ULB 0.133 STD 1.33)
    --K_vPPM=<_vK>        Concentration of K-40 in veto pmt glass (ULB 36.0 STD 260.0) 
    
    --U238_Gd=<_U238Gd>    Activity of U238 in Gd (upper) mBq/kg [Default: %f ]
    --Th232_Gd=<_Th232Gd>     Activity of Th232 in Gd (upper) mBq/kg [Default: %f]
    --U235_Gd=<_U235Gd>    Activity of U235 in Gd (upper) mBq/kg [Default: %f]
    --U238_Gd_l=<_U238Gd_l>    Activity of U238 in Gd (lower) mBq/kg [Default: %f]
    --Th232_Gd_l=<_Th232Gd_l>    Activity of Th232 in Gd (lower) mBq/kg [Default: %f]
    --U235_Gd_l=<_U235Gd_l>    Activity of U235 in Gd (lower) mBq/kg [Default: %f]

    --totRate		     ALTERNATIVE. Flag to use total rate from radiopurity google spreadsheet
    --U238_PMT=<_Upmt>  Total concentration of U-238 in glass [Default: 20.]
    --Th232_PMT=<_Tpmt> Total concentration of Th-232 in glass [Default: 20.]
    --K_PMT=<_Kpmt>     Total concentration of K-40 in glass [Default: 16.0]
    --Rn222_WV=<_RnWV>  Total Radon activity in water SK 2x10^-3 Bq/m^3 [Default: 20.]

    ## Misc options. Are in development or are historical, may not function or work as intended

    --sensitivity       Read analyzed results and evaluate sensitivity
    --efficiency        Read merged files and perform analysis
    --findRate          Find rate as a function of input flags.
    -n                  generate ntuple from single rat-pac root files
    --extractNtup       generate ntuple from all rat-pac root files
    -f=<ifile>          Input file [Default: %s]
    --ft=<ifile2>       Input file type [Default: %s]
    --ntupleout=<outN>  Name of ntuple out [Default: %s]
    -o=<outputfile>     Efficiency output file [Default: %s]
    --PDFs              Extract PDFs with series of cuts
    -r=<rate>           rate of accidentals in hz [Default: %f]
    -d=<distance>       Maximal distance between two events (m) [Default: %f]
    -t=<time>           Maximal time between two events (micro) [Default: %f]
    -T=<tubes>          Minimal number of tubes hit [Default: %d]
    --minN9=<_MPE>      Minimal number of photoelectron [Default: 9.]
    -g=<goodness>       Bonsai position goodness parameter [Default: %f]
    -G=<Goodness>       Bonsai direction goodness parameter [Default: %f]
    --RNRedux=<_RNR>    Reduction due to time/spatial veto arround (>0.90) [Default: 0.9]
    -P=<proc>           Pick a single physics process to analyis/merge (used for ntup)
    -L=<loc>            Pick a single physics location to analyis/merge (used for ntup)
    --timeScale=<_ts>   Integration period (sec,day,month,year) [Default: %s]
    --site=<_site>      Site of the experiment (boulby,fairport) [Default: %s]
    --OnOff=<_OOratio>  Ratio of reactor on to reactor off [Default: %d]
    --cores=<_cores>    Number of cores to discover [Default: 1]
    --bkgdSys=<_BSys>   Systematic background percentage [Default: 0.2]
    --40MWth            Option to sensitivity to do case scenarios
    --40MWthSig=<_SL>   Sigma discovery [Default: 3.0]


    """ % (defaultValues[0],defaultValues[1],defaultValues[2],defaultValues[3],defaultValues[4],\
           defaultValues[5],defaultValues[6],defaultValues[7],defaultValues[8],\
           defaultValues[9],defaultValues[10],defaultValues[11],defaultValues[12],\
           defaultValues[13],defaultValues[14],defaultValues[15],defaultValues[16],\
           defaultValues[17],defaultValues[18],defaultValues[19],defaultValues[20],\
           defaultValues[21],defaultValues[22],defaultValues[23],defaultValues[24],\
           defaultValues[25],defaultValues[26],defaultValues[27],defaultValues[28],\
	   defaultValues[29])

try:
    import docopt
    arguments = docopt.docopt(docstring)
    print '\nUsing docopt as the user control interface\n'
except ImportError:
    print 'docopt is not a recognized module, it is required to run this module'



if (arguments['--vetoThickZ'] and arguments['--vetoThickR']):
        print "Argument provided for veto radius and height"
else:
        print "No argument provided for veto radius and height, assumes shieldThick value for both"
        arguments['--vetoThickR'] = arguments['--shieldThick']
        arguments['--vetoThickZ'] = arguments['--shieldThick']

# Unless specified the default glass is standard in the vetos
if not (arguments['--U238_vPPM']):
        arguments['--U238_vPPM'] = '0.341' 
if not (arguments['--Th232_vPPM']):
        arguments['--Th232_vPPM'] = '1.33'
if not (arguments['--K_vPPM']):
        arguments['--K_vPPM'] = '260.0'
        
if (arguments['--vpc']):
     print "Arguments provided for the veto covereage"
else:
     arguments['--vpc'] = "0.002"
    
if arguments['--noRoot']:
    print 'Not loading any ROOT modules. usefull for generating files on oslic'
    
else:
    from ROOT import TRandom3
    from ROOT import TChain,TGraph,TGraphErrors,gSystem,gROOT,TH1D,TH2D,TFile,TCanvas,TF1
    from ROOT import THStack,Double
    from ROOT import kRed,kBlue,kGreen,kCyan,kOrange

    from ROOT import kOrange as kO,kBlue as kB,kGreen as kG
    from ROOT import kMagenta as kM,kAzure as kA,kRed as kR
    from ROOT import TCanvas,TLine, TLatex
    from ROOT import gStyle,gPad,TPaletteAxis

    gSystem.Load("$RATROOT/lib/libRATEvent")
    gSystem.AddIncludePath(" -I$RATROOT/include")

    # gROOT.LoadMacro("$WATCHENV/watchmakers/goldenFileExtractor.C")
    # from ROOT import goldenFileExtractor

    gStyle.SetOptStat(1112211)




def loadSimulationParametersNew():
    #Chain and subsequent isotopes
    d = {}

    d['CHAIN_238U_NA'] =['234Pa','214Pb','214Bi','210Bi','210Tl']
    d['CHAIN_232Th_NA'] = ['228Ac','212Pb','212Bi','208Tl']
    d['CHAIN_235U_NA'] = ['231Th','223Fr','211Pb','211Bi','207Tl']
    d['40K_NA']         = ['40K']
    d['CHAIN_222Rn_NA'] = ['214Pb','214Bi','210Bi','210Tl']
    d['STEEL_ACTIVITY'] = ['60Co','137Cs']

    d['ibd_p'] = ['promptPositron']
    d['ibd_n'] = ['delayedNeutron']
    d['pn_ibd']   = ['promptDelayedPair']
    # Fast neutron contamination
    d['FN'] = ['QGSP_BERT_EMV','QGSP_BERT_EMX','QGSP_BERT','QGSP_BIC','QBBC',\
    'QBBC_EMZ','FTFP_BERT','QGSP_FTFP_BERT']
    A,Z = ['9','11'],['3','3']
    ZA = A
    for i in range(len(A)):
        ZA[i] = str(int(A[i])*1000 +int(Z[i]))
    d['A_Z'] =  ZA

    process = {'40K_NA':['WaterVolume','PMT','VETO','CONCRETE','TANK','ROCK','SHIELD'],\
    'CHAIN_238U_NA':['PMT','VETO','CONCRETE','TANK','ROCK','GD','SHIELD'],\
    'CHAIN_232Th_NA':['PMT','VETO','CONCRETE','TANK','ROCK','GD','SHIELD'],\
    'CHAIN_235U_NA':['TANK','GD'],\
    'CHAIN_222Rn_NA':['WaterVolume'],\
    'STEEL_ACTIVITY':['TANK','CABLE'],\
    'FN':['ROCK','CONCRETE'],\
    'A_Z':['WaterVolume'],\
    'ibd_p':['WaterVolume'],\
    'ibd_n':['WaterVolume'],\
    'pn_ibd':['WaterVolume']}

    #Photocoverage selected
    # coverage = ['15pct','20pct','25pct','30pct','35pct','SuperK','WatchmanSphere']
    coverage = ['25pct']

    if arguments['-C']:
        coverage = [arguments['-C']]

    return d,process,coverage



def loadActivity():

    d,process,coverage = loadSimulationParametersNew()

    ##Evaluate the total mass of PMT glass in kg
    mass, diameter = 1.4, 10.0/0.039 #inch_per_mm # from Hamamatsu tech details
    areaPerPMT = pi*diameter*diameter/4.
    pmtRadius = float(arguments['--tankRadius'])-float(arguments['--steelThick'])-float(arguments['--vetoThickR'])
    pmtHeight = float(arguments['--halfHeight'])-float(arguments['--steelThick'])-float(arguments['--vetoThickZ'])
    psupArea = (2*pmtHeight)*2*pi*pmtRadius + 2.*(pi*pmtRadius**2)
    numPMTs = psupArea/areaPerPMT

    cPMTs = [float(s.strip('pct'))/100.*numPMTs for s in coverage]
    mPMTs = [s*mass for s in cPMTs]

    # print "Num pmts", cPMTs
    # print "Total Mass",mPMTs

    M_U238,Lambda_U238,Abund_U238 = 3.953e-25,4.916e-18,0.992745
    PPM_U238    = float(arguments["--U238_PPM"])
    ActivityU238= Lambda_U238*PPM_U238/M_U238/1e6
    mPMTsU238 = [s*ActivityU238 for s in mPMTs]
    print 'U238',mPMTsU238, ', PPM:',PPM_U238

    M_Th232,Lambda_Th232,Abund_Th232 = 3.853145e-25, 1.57e-18,1.0
    PPM_Th232    = float(arguments["--Th232_PPM"])
    ActivityTh232= Lambda_Th232*PPM_Th232/M_Th232/1e6
    mPMTsTh232 = [s*ActivityTh232 for s in mPMTs]
    print 'Th232',mPMTsTh232, ', PPM:',PPM_Th232

    M_K,Lambda_K,Abund_K = 6.636286e-26,1.842e-18,0.00117
    PPM_K    = float(arguments["--K_PPM"])
    ActivityK= Lambda_K*PPM_K/M_K/1e6
    mPMTsK = [s*ActivityK for s in mPMTs]
    print 'K',mPMTsK, ', PPM:',PPM_K

    print

    return mPMTs




def loadPMTActivity():

    d,process,coverage = loadSimulationParametersNew()

    ##Evaluate the total mass of PMT glass in kg
    mass, diameter = 1.4, 10.0/0.039 #inch_per_mm # from Hamamatsu tech details
    areaPerPMT = pi*diameter*diameter/4.
    pmtRadius = float(arguments['--tankRadius'])-float(arguments['--steelThick'])-float(arguments['--vetoThickR'])
    pmtHeight = float(arguments['--halfHeight'])-float(arguments['--steelThick'])-float(arguments['--vetoThickZ'])
    psupArea = (2*pmtHeight)*2*pi*pmtRadius + 2.*(pi*pmtRadius**2)
    numPMTs = psupArea/areaPerPMT
    cPMTs = [float(s.strip('pct'))/100.*numPMTs for s in coverage]
    mPMTs = [s*mass for s in cPMTs]

    print "Num pmts", cPMTs
    print "Total Mass",mPMTs

    M_U238,Lambda_U238,Abund_U238 = 3.953e-25,4.916e-18,0.992745
    PPM_U238    = float(arguments["--U238_PPM"])
    ActivityU238= Lambda_U238*PPM_U238/M_U238/1e6
    mPMTsU238 = [s*ActivityU238 for s in mPMTs]
    print 'U238',mPMTsU238, ', PPM:',PPM_U238,'activity per PMT:', ActivityU238*mass,'Bq per PMT per isotope in chain'

    M_Th232,Lambda_Th232,Abund_Th232 = 3.853145e-25, 1.57e-18,1.0
    PPM_Th232    = float(arguments["--Th232_PPM"])
    ActivityTh232= Lambda_Th232*PPM_Th232/M_Th232/1e6
    mPMTsTh232 = [s*ActivityTh232 for s in mPMTs]
    print 'Th232',mPMTsTh232, ', PPM:',PPM_Th232, 'activity per PMT:', ActivityTh232*mass,'Bq per PMT per isotope in chain'

    M_K,Lambda_K,Abund_K = 6.636286e-26,1.842e-18,0.00117
    PPM_K    = float(arguments["--K_PPM"])
    ActivityK= Lambda_K*PPM_K/M_K/1e6*Abund_K
    mPMTsK = [s*ActivityK for s in mPMTs]
    print 'K',mPMTsK, ', PPM:',PPM_K, 'activity per PMT:', ActivityK*mass,'Bq per PMT per isotope in chain'

    print

    return mPMTs,mPMTsU238,mPMTsTh232,mPMTsK

def loadVETOActivity():
    d,process,coverage = loadSimulationParametersNew()
    if (arguments['--vetoModel']=="r11780_hqe" or arguments['--vetoModel']=="r11780_lqe"):
        mass, diameter = 2.0, 12.0/0.039 #inch_per_mm
    elif (arguments['--vetoModel']=="r7081_hqe" or arguments['--vetoModel']=="r7081_lqe"):
        mass, diameter = 1.4, 10.0/0.039 #inch_per_mm
    else: #assuming this would be the ETL pmt
        mass, diameter = 2.0, 11.0/0.039
    areaPerPMT = pi*diameter*diameter/4.0
    pmtRadius = float(arguments['--tankRadius'])-float(arguments['--steelThick'])-float(arguments['--vetoThickR'])+700.0#still need to find dim
    pmtHeight =float(arguments['--halfHeight'])-float(arguments['--steelThick'])-float(arguments['--vetoThickZ'])+700.0#still need to find dim
    psupArea = (2*pmtHeight)*2*pi*pmtRadius + 2.*(pi*pmtRadius**2)
    numPMTs = psupArea/areaPerPMT
    cVETOs = float(arguments["--vpc"])*numPMTs
    mVETOs = mass*cVETOs

    print "Number of Veto", cVETOs
    print "Total Mass", mVETOs

    M_U238,Lambda_U238,Abund_U238 = 3.953e-25,4.916e-18,0.992745
    PPM_U238    = float(arguments["--U238_vPPM"])
    ActivityU238= Lambda_U238*PPM_U238/M_U238/1e6
    mVETOsU238 = mVETOs*ActivityU238
    print 'U238',mVETOsU238, ', PPM:',PPM_U238,'activity per PMT:', ActivityU238*mass,'Bq per PMT per isotope in chain'

    M_Th232,Lambda_Th232,Abund_Th232 = 3.853145e-25, 1.57e-18,1.0
    PPM_Th232    = float(arguments["--Th232_vPPM"])
    ActivityTh232= Lambda_Th232*PPM_Th232/M_Th232/1e6
    mVETOsTh232 = mVETOs*ActivityTh232
    print 'Th232',mVETOsTh232, ', PPM:',PPM_Th232, 'activity per PMT:', ActivityTh232*mass,'Bq per PMT per isotope in chain'


    M_K,Lambda_K,Abund_K = 6.636286e-26,1.842e-18,0.00117
    PPM_K    = float(arguments["--K_PPM"])
    ActivityK= Lambda_K*PPM_K/M_K/1e6*Abund_K
    mVETOsK = mVETOs*ActivityK
    print 'K',mVETOsK, ', PPM:',PPM_K, 'activity per PMT:', ActivityK*mass,'Bq per PMT per isotope in chain'

    print

    return mVETOs,mVETOsU238,mVETOsTh232,mVETOsK

def loadTankActivity():                                         ##added by Leah: activity from steel in tank
    ##MASS OF STEEL USED IN KG -- assuming use of steel grade 304
    density = 8000                                              ##kg/m^3
    r1 = (float(arguments["--tankRadius"]))/1000.               ##outer radius inc steel thick in m
    h1 = 2 * (float(arguments["--halfHeight"]))/1000.           ##outer height inc steel thick in m
    V1 = pi * h1 * r1**2                                        ##outer volume inc steel thick in m^3

    r2 = r1 - (float(arguments["--steelThick"]))/1000.         ##inner radius in m
    h2 = h1 - 2*(float(arguments["--steelThick"]))/1000.       ##inner height in m
    V2 = pi * h2 * r2**2                                        ##inner volume in m^3

    tankvol = V1 - V2                                           ##hollow cylinder m^3
    tankmass = tankvol*density

    print "Total steel mass",tankmass,"kg"

    ##STAINLESS STEEL CONTAINS 60CO, 137CS
    act_60co = 19e-3                                           ##mean value from "Measurements of extremely low radioactivity in stainless steel" paper Bq per kg
    tankact_60co = tankmass * act_60co                          ##activity in Bq
    print "60Co in tank steel:\n activity per kilogram:", act_60co,"Bq/kg,\n total activity:",tankact_60co,"Bq"

    act_137cs = 0.77e-3                                         ##mean value from "Measurements of extremely low radioactivity in stainless steel" paper Bq per kg
    tankact_137cs = tankmass * act_137cs                        ##activity in Bq
    print "137Cs in tank steel:\n activity per kilogram:", act_137cs,"Bq/kg,\n total activity:",tankact_137cs,"Bq"

    return tankmass,tankact_60co,tankact_137cs


def loadConcreteActivity():                                     ##added by Leah: activity from concrete
    ##MASS OF CONCRETE USED IN KG -- assuming normal-weight concrete (NWC)
    density = 2300                                              ##kg/m^3
    thickness = 0.5                                             ## slab thickness in metres - alter for desired value
    concvol = 25.5*(pi*pow(13.,2)-pi*pow(12.5,2)) + 0.5*pi*(pow(13.,2)) #0.5m thick, 25m high concrete 'tube': outer diameter 26m, 
									#inner diameter 25m, plus 0.5m base of diameter 26m
    concmass = concvol * density

    print "Total concrete slab mass",concmass,"kg"

##CONCRETE CONTAINS 238U, 232TH, 40K

    act_238u = 61                                              ##mean UK value from "Natural radioactivity in building materials in the European Union" paper Bq per kg
    concact_238u = concmass * act_238u                         ##activity in Bq
    print "238U in concrete slab:\n activity per kilogram:", act_238u,"Bq/kg,\n total activity:",concact_238u,"Bq"

    act_232th = 30                                              ##mean UK value from "Natural radioactivity in building materials in the European Union" paper Bq per kg
    concact_232th = concmass * act_232th                        ##activity in Bq
    print "232Th in concrete slab:\n activity per kilogram:", act_232th,"Bq/kg,\n total activity:",concact_232th,"Bq"

    act_40k = 493                                               ##mean UK value from "Natural radioactivity in building materials in the European Union" paper Bq per kg
    concact_40k = concmass * act_40k                            ##activity in Bq
    print "40K in concrete slab:\n activity per kilogram:", act_40k,"Bq/kg,\n total activity:",concact_40k,"Bq"

    return concmass,concact_238u,concact_232th,concact_40k
###gravel under concrete? 24in


def loadShotcreteActivity():
    ##MASS OF SHOTCRETE USED IN KG -- assuming same shotcrete as used in snolab (Shotcrete Application in SNOLAB paper)
    density = 2400                                              ##kg/m^3 -- approx from "Shotcrete in Tunnel Construction - Introduction to the basic technology of sprayed concrete"
    thickness = 0.1                                             ##thickness in metres
    r1 = 25                                                     ##outer radius in m
    h1 = 25                                                     ##outer height in m
    V1 = pi * h1 * r1**2                                        ##outer volume in m^3

    r2 = r1 - thickness                                         ##inner radius in m
    h2 = h1 - 2*thickness                                       ##inner height in m
    V2 = pi * h2 * r2**2                                        ##inner volume in m^3

    shotvol = V1 - V2                                           ##hollow cylinder m^3
    shotmass = shotvol*density

    print "Total shotcrete mass",shotmass,"kg"

    mass_238u = 3.953e-5
    lambda_238u = 4.916e-18
    abund_238u = 0.992745
    ppm_238u = 2.6
    act_238u =( (lambda_238u * ppm_238u)/mass_238u/1e6)*shotmass
    print "238U in shotcrete coating:\n ppm:", ppm_238u,"\n total activity:",act_238u,"Bq"

    mass_232th = 3.853145e-25
    lambda_232th = 1.57e-18
    abund_232th = 1.0
    ppm_232th = 14
    act_232th = ((lambda_232th * ppm_232th)/mass_232th/1e6)*shotmass
    print "232Th in shotcrete coating:\n ppm:", ppm_232th,"\n total activity:",act_232th,"Bq"

    mass_40k = 6.636286e-26
    lambda_40k = 1.842e-18
    abund_40k = 0.00117
    ppm_40k = 16000           ##1.6%
    act_40k = ((lambda_40k * ppm_40k)/mass_40k/1e6)*shotmass
    print "40K in shotcrete coating:\n ppm:", ppm_40k,"\n total activity:",act_40k,"Bq"

    return shotmass,act_238u,act_232th,act_40k



def loadRockActivity():
    ##NB MOVED TO ROCK SALT CAVERN RATHER THAN ROCK...
    ##5M THICKNESS
    ##MASS OF SALT IN WALL IN KG -- assuming pure rock salt walls (Overview of the European Underground Facilities paper)
    density = 2165                                              ##kg/m^3 -- approx from "Physical Properties Data for Rock Salt"
    thickness = 5                                               ##thickness in metres
    r1 = 13 + thickness                                         ##outer radius in m
    h1 = 25.5 + 2*thickness                                     ##outer height in m
    V1 = pi * h1 * r1**2                                        ##outer volume in m^3

    r2 = 13                                                     ##inner radius in m
    h2 = 25.5                                                   ##inner height in m
    V2 = pi * h2 * r2**2                                        ##inner volume in m^3

    rockvol = V1 - V2                                           ##hollow cylinder m^3 outside 25m cylindrical cavern 
								##plus 0.5m concrete layer on walls and floor 
    rockmass = rockvol*density

    print "Total relevant rock mass",rockmass,"kg"

    mass_238u = 3.953e-5
    lambda_238u = 4.916e-18
    abund_238u = 0.992745
    ppm_238u = 0.067
    act_238u =( (lambda_238u * ppm_238u)/mass_238u/1e6)*rockmass
    print "238U in shotcrete coating:\n ppm:", ppm_238u,"\n total activity:",act_238u,"Bq"

    mass_232th = 3.853145e-25
    lambda_232th = 1.57e-18
    abund_232th = 1.0
    ppm_232th = 0.125
    act_232th = ((lambda_232th * ppm_232th)/mass_232th/1e6)*rockmass
    print "232Th in shotcrete coating:\n ppm:", ppm_232th,"\n total activity:",act_232th,"Bq"

    mass_40k = 6.636286e-26
    lambda_40k = 1.842e-18
    abund_40k = 0.00117
    ppm_40k = 1130
    act_40k = ((lambda_40k * ppm_40k)/mass_40k/1e6)*rockmass
    print "40K in shotcrete coating:\n ppm:", ppm_40k,"\n total activity:",act_40k,"Bq"

    return rockmass,act_238u,act_232th,act_40k


def loadGdActivity():

    d,process,coverage = loadSimulationParametersNew()
    tankRadius = float(arguments['--tankRadius']) - float(arguments['--steelThick'])
    halfHeight = float(arguments['--halfHeight']) - float(arguments['--steelThick'])
    nKiloTons = pi*pow(tankRadius/1000.,2)*(2.*halfHeight/1000.)/1000.
    GdU238    = float(arguments["--U238_Gd"]) / 1000. * nKiloTons * 1e6 * 0.002 # bq/kg * kg of water * Gd(SO4)3 concentration
    GdTh232   = float(arguments["--Th232_Gd"])/ 1000. * nKiloTons * 1e6 * 0.002 # bq/kg * kg of water * Gd(SO4)3 concentration
    GdU235    = float(arguments["--U235_Gd"]) / 1000. * nKiloTons * 1e6 * 0.002  #bq/kg * kg of water * Gd(SO4)3 concentration
    GdU238_l    = float(arguments["--U238_Gd_l"]) / 1000. * nKiloTons * 1e6 * 0.002 # bq/kg * kg of water * Gd(SO4)3 concentration
    GdTh232_l   = float(arguments["--Th232_Gd_l"])/ 1000. * nKiloTons * 1e6 * 0.002 # bq/kg * kg of water * Gd(SO4)3 concentration
    GdU235_l    = float(arguments["--U235_Gd_l"]) / 1000. * nKiloTons * 1e6 * 0.002  #bq/kg * kg of water * Gd(SO4)3 concentration


    return GdU238,GdTh232,GdU235,GdU238_l,GdTh232_l,GdU235_l





