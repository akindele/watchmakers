from .load import *
from ROOT import TChain, gROOT, TFile



# The purpose of this class is to handle the input/ouput operations of
# WATCHMAKERS (WM). This include creating directories and files for the different
# operations of WM, such as creating macros, jobs and so forth.

def testCreateDirectory(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
    else:
        rmtree(directory)
        os.makedirs(directory)

def deleteDirectory(directory):
    if os.path.exists(directory):
        rmtree(directory)

def testCreateDirectoryIfNotExist(directory):

    if os.path.exists(directory):
        print('''There is already a directory here. %s
No new directory has been made. Bad idea. Consider saving current files
and using --force.        \n'''%(directory))
    if not os.path.exists(directory):
        os.makedirs(directory)

def macroGeneratorNew(percentage,location,element,_dict,runs,events,dirOpt):

    covPCT = {'10pct':0.1,'15pct':0.15,'20pct':0.2,\
    '25pct':0.25,'30pct':0.30,'35pct':0.35,'40pct':0.40,'SuperK':0.39,'WatchmanSphere':0.20}

    additionalString,additionalCommands,additionalMacStr,additionalMacOpt = testEnabledCondition(arguments)
    d,process,coverage = loadSimulationParametersNew()
    # loadActivity()

    #Part of the macro that is the same for all jobs
    dir = os.getcwd()

    depth = float(arguments["--depth"])
    rate = 1.0
    if 'pct' in percentage:
        header = '''
/glg4debug/glg4param omit_muon_processes  0.0
/glg4debug/glg4param omit_hadronic_processes  0.0

/rat/db/set DETECTOR experiment "Watchman"
/rat/db/set DETECTOR detector_factory "Watchman"
/rat/db/set WATCHMAN_PARAMS photocathode_coverage %4.2f
%s

/run/initialize

# BEGIN EVENT LOOP
/rat/proc lesssimpledaq
# /rat/proc fitbonsai
# /rat/proc fitcentroid
# /rat/proc fitpath
/rat/proc count
/rat/procset update 1000

# Use IO.default_output_filename
/rat/proclast outroot
# /rat/procset file "%s/root_files%s%sfile_%d.root" Now set with -o object
#END EVENT LOOP
''' %(covPCT[percentage],additionalMacOpt,dir,additionalMacStr,dirOpt,runs)
    else:
        if location == 'GUNITE' or location == 'CONCRETE':
            return ''
        else:
            header = '''
/glg4debug/glg4param omit_muon_processes  0.0
/glg4debug/glg4param omit_hadronic_processes  0.0

/rat/db/set DETECTOR experiment "%s"
/rat/db/set DETECTOR geo_file "%s/%s.geo"

/run/initialize

# BEGIN EVENT LOOP
/rat/proc lesssimpledaq
# /rat/proc fitbonsai
# /rat/proc fitcentroid
# /rat/proc fitpath
/rat/proc count
/rat/procset update 1000

# Use IO.default_output_filename
/rat/proclast outroot
# /rat/procset file "%s/root_files%s%sfile_%d.root" Now set with -o object
#END EVENT LOOP
''' %(percentage,percentage,percentage,dir,additionalMacStr,dirOpt,runs)

    if element in d['FN']:
        line1 = '''
/generator/add combo fastneutron:regexfill:poisson
/generator/pos/set rock_[0-9]+
/generator/vtx/set 0 0 0
/generator/fastneutron/depth %f
/generator/fastneutron/enthresh 10.0
/generator/fastneutron/sidewalls 1.0
/generator/rate/set %f
/run/beamOn %d'''%(depth,rate,events)

    elif element in d['CHAIN_238U_NA'] or element in d['CHAIN_232Th_NA'] or element in d['40K_NA'] or element in d['STEEL_ACTIVITY'] or element in d['CHAIN_235U_NA']:
        if location == 'PMT':
            line1 = '''
/generator/add decaychain %s:regexfill:poisson
/generator/pos/set inner_pmts[0-9]+
/generator/rate/set %f
/run/beamOn %d''' %(element,rate,events*2)
        elif location == 'VETO':
            line1 = '''
/generator/add decaychain %s:regexfill:poisson
/generator/pos/set veto_pmts[0-9]+
/generator/rate/set %f
/run/beamOn %d''' %(element,rate,events*2)
        else:
            locat = location.lower()
            if locat == 'watervolume' or locat == 'gd':
                locat = 'detector'
                xTimes = 2
            else:
                xTimes = 10
            line1 = '''
/generator/add decaychain %s:regexfill:poisson
/generator/pos/set %s+
/generator/rate/set %f
/run/beamOn %d''' %(element,locat,rate,events*xTimes)

    elif element in d['ibd_p']:
        line1 ='''
/generator/add combo spectrum:fill:poisson
/generator/vtx/set e+ %s
/generator/pos/set 0 0 0
/generator/rate/set %f
/run/beamOn %d'''%(element,rate,events/3)

    elif element in d['ibd_n']:
        line1 = '''
/generator/add combo gun2:fill:poisson
/generator/vtx/set %s  0 0 0 0 0.001 0.20
/generator/pos/set 0 0 0
/generator/rate/set %f
/run/beamOn %d'''%('neutron',rate,events/3)

    elif element in d['pn_ibd']:
        line1 = '''
/generator/add combo ibd:fill:poisson
/generator/vtx/set %s  1 0 0
/generator/pos/set 0 0 0
/generator/rate/set %f
/run/beamOn %d'''%(element,rate,events/3)

    elif element in d['A_Z']:
        A =  int(int(element)/1000)
        Z = int(element) - A*1000
        line1 = '''
/generator/add combo isotope:fill:poisson
/generator/vtx/set GenericIon  1 0 0
/generator/isotope/A %s.0
/generator/isotope/Z %s.0
/generator/isotope/E 0.0
/generator/pos/set 0 0 0
/generator/rate/set %f
/run/beamOn %d'''%(A,Z,rate,events/4)


    else:
        print('Could not find ',element,location.lower())
        line1 = ''
    return header+line1


def jobString(percentage,j,runs,models,arguments):
#    directory = "/p/lscratche/adg/Watchboy/simplifiedData/rp_sim/wm"
    directory   = os.getcwd()
    softDir     = "/usr/gapps/adg/geant4/rat_pac_and_dependency"
    ratDir      = os.environ['RATROOT']
    rootDir     = os.environ['ROOTSYS']
    g4Dir       =  os.environ['G4INSTALL']
    watchmakersDir = os.environ['WATCHENV']
    try:
        sheffield   = os.environ['SHEFFIELD']
        # print 'Running on sheffield cluster'
    except:
        # print 'Not running on sheffield cluster'
        sheffield =  0

    software    = "%s/bin/rat" %(ratDir)
    # d,iso,loc,coverage,coveragePCT = loadSimulationParameters()

    if arguments['--newVers']:
        d,proc,coverage = loadSimulationParameters()
    else:
        d,iso,loc,coverage,coveragePCT = loadSimulationParameters()

    ele =  d["%s"%(iso[j])]
    location = loc[j]

    goodness     = float(arguments['-g'])
    case = int(arguments['-j'])

    additionalString,additionalCommands,additionalMacStr,additionalMacOpt = testEnabledCondition(arguments)


    N            = int(arguments["-N"])
    rate         = float(arguments["-r"])
    timemask     = float(arguments['-t'])*1000.0
    distancemask = float(arguments['-d'])
    goodness     = float(arguments['-g'])
    dirGoodness  = float(arguments['-G'])
    minNHIT      = float(arguments['-T'])
    fIn          = arguments["-f"]
    fidR         = (float(arguments["--tankRadius"])-float(arguments["--vetoThickR"])-float(arguments["--steelThick"])-float(arguments["--fidThick"]))/1000.
    fidZ         = (float(arguments["--halfHeight"])-float(arguments["--vetoThickZ"])-float(arguments["--steelThick"])-float(arguments["--fidThick"]))/1000.
    pmtR         = (float(arguments["--tankRadius"])-float(arguments["--steelThick"])-float(arguments["--vetoThickR"]))/1000.
    pmtZ         = (float(arguments["--halfHeight"])-float(arguments["--steelThick"])-float(arguments["--vetoThickZ"]))/1000.
    tankR        = float(arguments["--tankRadius"])/1000.
    tankZ        = float(arguments["--halfHeight"])/1000.
    outF         = arguments["--ntupleout"]
    superNova    = arguments["--supernovaFormat"]

    if sheffield:

        line1 = """#!/bin/sh
#MSUB -N WM_%s_%s_%d_%s    #name of job
#MSUB -A ait         # sets bank account
#MSUB -l nodes=1:ppn=1,walltime=73:59:59,partition=borax  # uses 1 node
#MSUB -q pbatch         #pool
#MSUB -o %s/log_case%s%s/wmpc_%s_%s_%d.log
#MSUB -e %s/log_case%s%s/wmpc_%s_%s_%d.err
#MSUB -d %s  # directory to run from
#MSUB -V
#MSUB                     # no more psub commands

source %s/bin/thisroot.sh
source /cvmfs/sft.cern.ch/lcg/releases/LCG_85swan2/gcc/4.9.3/x86_64-slc6/setup.sh
source /usr/local/geant4/setup.sh 10.2
source %s/geant4make.sh
source %s/env.sh
source %s/env_wm.sh
export G4NEUTRONHP_USE_ONLY_PHOTONEVAPORATION=1
export SHEFFIELD=1\n
    """ %(percentage,location,runs,additionalMacStr,\
    directory,case,additionalMacStr,percentage,location,runs,\
    directory,case,additionalMacStr,percentage,location,runs,\
    directory,\
    rootDir,g4Dir,ratDir,watchmakersDir)

    elif arguments["--docker"]:
        line1 = """#!/bin/sh
#MSUB -N WM_%s_%s_%d_%s    #name of job
#MSUB -A ait         # sets bank account
#MSUB -l nodes=1:ppn=1,walltime=73:59:59,partition=borax  # uses 1 node
#MSUB -q pbatch         #pool
#MSUB -o %s/log_case%s%s/wmpc_%s_%s_%d.log
#MSUB -e %s/log_case%s%s/wmpc_%s_%s_%d.err
#MSUB -d %s  # directory to run from
#MSUB -V
#MSUB                     # no more psub commands

#source %s/bin/thisroot.sh
#source %s/../../../bin/geant4.sh
#source %s/geant4make.sh
#source %s/env.sh
#source %s/env_wm.sh
source $HOME/software/docker_watchman/env.sh
export G4NEUTRONHP_USE_ONLY_PHOTONEVAPORATION=1\n
    """ %(percentage,location,runs,additionalMacStr,\
    directory,case,additionalMacStr,percentage,location,runs,\
    directory,case,additionalMacStr,percentage,location,runs,\
    directory,\
    rootDir,g4Dir,g4Dir,ratDir,watchmakersDir)

    elif  arguments["--singularity"]:
        line1 = """#!/bin/sh
#MSUB -N WM_%s_%s_%d_%s    #name of job
#MSUB -A ait         # sets bank account
#MSUB -l nodes=1:ppn=1,walltime=73:59:59,partition=borax  # uses 1 node
#MSUB -q pbatch         #pool
#MSUB -o %s/log_case%s%s/wmpc_%s_%s_%d.log
#MSUB -e %s/log_case%s%s/wmpc_%s_%s_%d.err
#MSUB -d %s  # directory to run from
#MSUB -V
#MSUB                     # no more psub commands

#source %s/bin/thisroot.sh
#source %s/../../../bin/geant4.sh
#source %s/geant4make.sh
#source %s/env.sh
#source %s/env_wm.sh

export G4NEUTRONHP_USE_ONLY_PHOTONEVAPORATION=1\n
    """ %(percentage,location,runs,additionalMacStr,\
    directory,case,additionalMacStr,percentage,location,runs,\
    directory,case,additionalMacStr,percentage,location,runs,\
    directory,\
    rootDir,g4Dir,g4Dir,ratDir,watchmakersDir)

    else:
        line1 = """#!/bin/sh
#MSUB -N WM_%s_%s_%d_%s    #name of job
#MSUB -A ait         # sets bank account
#MSUB -l nodes=1:ppn=1,walltime=73:59:59,partition=borax  # uses 1 node
#MSUB -q pbatch         #pool
#MSUB -o %s/log_case%s%s/wmpc_%s_%s_%d.log
#MSUB -e %s/log_case%s%s/wmpc_%s_%s_%d.err
#MSUB -d %s  # directory to run from
#MSUB -V
#MSUB                     # no more psub commands

source %s/bin/thisroot.sh
source %s/../../../bin/geant4.sh
source %s/geant4make.sh
source %s/env.sh
source %s/env_wm.sh
export G4NEUTRONHP_USE_ONLY_PHOTONEVAPORATION=1\n
    """ %(percentage,location,runs,additionalMacStr,\
    directory,case,additionalMacStr,percentage,location,runs,\
    directory,case,additionalMacStr,percentage,location,runs,\
    directory,\
    rootDir,g4Dir,g4Dir,ratDir,watchmakersDir)

    if arguments['--newVers']:
        print()
    else:
        for mods in models:
            if location == "FN":
                line1 += "export PHYSLIST=%s\n" %(mods)
            if case == 1 or case == 2 or case == 4:
                _log = "log_case%s%s/%s/%s/rat.%s_%s_%s_%d.log" %(case,additionalMacStr,mods,percentage,percentage,mods,location,runs)
                _mac = "%s/macro%s/%s/%s/run%s_%s_%d.mac" %(directory,additionalMacStr,mods,percentage,mods,location,runs)
                line1 += "%s -l %s %s\n" %(software,_log,_mac)
            if case >= 3:
                fileN = "root_files%s/%s/%s/watchman_%s_%s_%s_%d.root" %(additionalString,mods,percentage,mods,percentage,location,runs)
                if additionalString != "":
                    fileNO = "ntuple_root_files%s/%s/%s/watchman_%s_%s_%s%s_%d.root" %(additionalString,mods,percentage,mods,percentage,location,additionalString,runs)
                    if sheffield:
                        line1 += "root -b -l -q %s/watchmakers/\'goldenFileExtractor.C(\"%s\",\"%s\",%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\")\'\n" %(watchmakersDir,fileN,fileNO,minNHIT,goodness,dirGoodness,timemask,\
                                        rate,distancemask,fidR,fidZ,pmtR,pmtZ,tankR,tankZ)
                        ###int goldenFileExtractor(const char *file, const char *outfile = "null", double nhit_min =3., double goodness_min = 0.1, double goodness_dir = 0.1, double timeWindow_ns = 100000, double rate = 10.0, double maxDistance = 2.0, double fidBoundR = 5.42, double fidBoundZ = 5.42, double pmtBoundR = 6.42, double pmtBoundZ = 6.42, double tankBoundR = 8.02635, double tankBoundZ = 8.02635 ) {
                    else:
                        line1 += "watch -n %s -f %s --ntupleout %s\n" %(additionalCommands,fileN,fileNO)
    return line1,case



def generateMacros(N,e):
    if arguments['--newVers']:
        d,proc,coverage = loadSimulationParameters()
    else:
        d,iso,loc,coverage,coveragePCT = loadSimulationParameters()

    additionalString,additionalCommands,additionalMacStr,additionalMacOpt = testEnabledCondition(arguments)

    print(additionalMacOpt)
    print(N,e)

    ##Clean or create macro directories
    for j in range(len(iso)):
        for ii in d["%s"%(iso[int(j)])]:
            for idx,cover in enumerate(coverage):
                dir = "macro%s/%s/%s" %(additionalMacStr,ii,cover)
                testCreateDirectory(dir)

    for j in range(len(iso)):
        for ii in d["%s"%(iso[int(j)])]:
            for idx,cover in enumerate(coverage):
                for val in range(N):
                    line = macroGenerator(cover,ii,loc[j],val,e )
                    dir = "macro%s/%s/%s" %(additionalMacStr,ii,cover)

                    outfile = open("%s/run%s_%s_%d.mac" %(dir,ii,\
                    loc[j],val),"w+")
                    outfile.writelines(line)
                    outfile.close
    return 0


def generateMacrosNew(N,e):

    d,proc,coverage = loadSimulationParametersNew()
    additionalString,additionalCommands,additionalMacStr,additionalMacOpt = testEnabledCondition(arguments)

    print(additionalMacOpt)
    print(N,e)

    cnt = 0

    for _p in proc:
        for _loc in proc[_p]:
            for idx,_cover in enumerate(coverage):
                for _element in d[_p]:
                    cond1 = ('pct' in _cover)==0 and (_loc == 'GUNITE' or _loc == 'CONCRETE')!=1
                    cond2 = ('pct' in _cover)==1
                    print(_cover,_loc,cond1,cond2, cond1 or cond2)
                    if cond1 or cond2:
                        for i in range(int(N/10+1)):
                            dir = "macro%s/%s/%s/%s/%s/run%08d/"%(additionalMacStr,_cover,_loc,_element,_p,i*10)
                            testCreateDirectory(dir)
                            cnt+=1
                        for val in range(N):
                            i = int(val/10)
                            print(_p," ",i," ",val)
                            dir = "%s/%s/%s/%s/%s/run%08d/"%(additionalMacStr,_cover,_loc,_element,_p,i*10)
                            outfile = open("macro%s/run_%08d.mac" %(dir,val),"w+")
                            line = macroGeneratorNew(_cover,_loc,_element,_p,val,e,dir)
                            outfile.writelines(line)
                            outfile.close


    return 0



def generateJobsNew(N,arguments):

    d,proc,coverage = loadSimulationParametersNew()

    case = arguments["-j"]

    additionalString,additionalCommands,additionalMacStr,additionalMacOpt = testEnabledCondition(arguments)

    try:
        sheffield   = os.environ['SHEFFIELD']
        # print 'Running on sheffield cluster'
    except:
        # print 'Not running on sheffield cluster'
        sheffield =  0

    '''Find wheter the jobs folder exist: if not create, if yes clean and recreate'''


    for _p in proc:
        for _loc in proc[_p]:
            for idx,_cover in enumerate(coverage):
                for _element in d[_p]:
                    # print cnt,_p,element,_loc,cover
                    print(_p,_loc,_cover,_element)
                    for i in range(int(N/10+1)):
                        dir = "root_files%s/%s/%s/%s/%s/run%08d"%(additionalMacStr,_cover,_loc,_element,_p,i*10)
                        if arguments['--force']:
                            print('Using force to recreate dir:',dir)
                            testCreateDirectory(dir)
                        else:
                            testCreateDirectoryIfNotExist(dir)
                        dir = "bonsai_root_files%s/%s/%s/%s/%s/run%08d"%(additionalMacStr,_cover,_loc,_element,_p,i*10)
                        if arguments['--force']:
                            print('Using force to recreate dir:',dir)
                            testCreateDirectory(dir)
                        else:
                            testCreateDirectoryIfNotExist(dir)
                        dir = "log%s/%s/%s/%s/%s/run%08d"%(additionalMacStr,_cover,_loc,_element,_p,i*10)
                        if arguments['--force']:
                            print('Using force to recreate dir:',dir)
                            testCreateDirectory(dir)
                        else:
                            testCreateDirectoryIfNotExist(dir)
                        #dir = "jobs%s/%s/%s/%s/%s/run%08d"%(additionalMacStr,_cover,_loc,_element,_p,i*10)
                        #if arguments['--force']:
                        #    print 'Using force to recreate dir:',dir
                        #    testCreateDirectory(dir)
                        #else:
                        #    testCreateDirectoryIfNotExist(dir)


    '''Make sure that the softlink are correct for Bonsai input'''

    ratDir      = os.environ['RATROOT']
    if arguments["--docker"]:

        src = '%s/fit_param.dat'%(ratDir)
        dst = os.getcwd()+'/fit_param.dat'
        if not os.path.exists(dst):
            os.system("rsync -avz %s %s "%(src,dst))

        src = '%s/like.bin'%(ratDir)
        dst = os.getcwd()+'/like.bin'
        if not os.path.exists(dst):
            os.system("rsync -avz %s %s "%(src,dst))
    elif arguments["--singularity"]:
        src = '/src/rat-pac/fit_param.dat'
        dst = '%s/fit_param.dat'%(os.environ['PWD'])
        if not os.path.exists(dst):
            os.system("singularity exec --bind $PWD %s /usr/bin/rsync -avz %s %s "%(arguments['--simg'],src,dst))

        src = '/src/rat-pac/like.bin'
        dst = '%s/like.bin'%(os.environ['PWD'])
        if not os.path.exists(dst):
            os.system("singularity exec --bind $PWD %s /usr/bin/rsync -avz %s %s "%(arguments['--simg'],src,dst))

    else:
        src = ratDir+'/fit_param.dat'
        dst = os.getcwd()+'/fit_param.dat'
        if not os.path.exists(dst):
            os.symlink(src,dst)

        src = ratDir+'/like.bin'
        dst = os.getcwd()+'/like.bin'
        if not os.path.exists(dst):
            os.symlink(src,dst)



    rootDir     = os.environ['ROOTSYS']
    softDir     = os.environ['RATROOT']
    rootDir     = os.environ['ROOTSYS']
    g4Dir       =  os.environ['G4INSTALL']
    watchmakersDir = os.environ['WATCHENV']
    directory   = os.getcwd()

    if arguments["--singularity"]:
        srat    = "singularity exec --bind $PWD %s /src/rat-pac/bin/rat" %(arguments["--simg"])
        sbonsai = "singularity exec --bind $PWD %s /src/rat-pac/tools/bonsai/bonsai"%(arguments["--simg"])


    outfile_jobs = open('sub_job_%s'%(additionalMacStr),"w+")
    outfile_jobs.writelines("#!/bin/sh\n")
    for _p in proc:
        for _loc in proc[_p]:
            for idx,_cover in enumerate(coverage):
                for _element in d[_p]:
                    # print cnt,_p,element,_loc,cover
                    cond1 = ('pct' in _cover)==0 and (_loc == 'GUNITE' or _loc == 'CONCRETE')!=1
                    cond2 = ('pct' in _cover)==1
                    if cond1 or cond2:
                        dir = "jobs%s/%s/%s/%s/%s"%(additionalMacStr,_cover,_loc,_element,_p)
                        if arguments['--force']:
                            print('Using force to recreate dir:',dir)
                            testCreateDirectory(dir)
                        else:
                            testCreateDirectoryIfNotExist(dir)
                        if sheffield:
                            outfile_jobs.writelines('condor_qsub %s\n'%(dir+'/job%08d.sh'%(0)))
                        else:
                            outfile_jobs.writelines('msub %s\n'%(dir+'/job%08d.sh'%(0)))
                        log = "log%s/%s/%s/%s/%s/log"%(additionalMacStr,_cover,_loc,_element,_p)
                        for i in range(int(N/10+1)):
                            dir = "jobs%s/%s/%s/%s/%s"%(additionalMacStr,_cover,_loc,_element,_p)
                            outfile = open(dir+'/job%08d.sh'%(i*10),"w+")
                            #os.chmod(dir+'/job%08d.sh'%(i*10),S_IRWXG)
                            #os.chmod(dir+'/job%08d.sh'%(i*10),S_IRWXU)
                            job_line = "%s_%s_%s_%s_%s_%s"%(additionalMacStr,_cover,_loc,_element,_p,i*10)
                            if arguments["--docker"]:
                                outfile.writelines("""#!/bin/sh
#MSUB -N job_%s    #name of job
#MSUB -A ait         # sets bank account
#MSUB -l nodes=1:ppn=1,walltime=73:59:59,partition=borax  # uses 1 node
#MSUB -q pbatch         #pool
#MSUB -o %s
#MSUB -e %s
#MSUB -d %s  # directory to run from
#MSUB -V
#MSUB                     # no more psub commands

source $HOME/software/docker_watchman/env.sh

"""%(job_line,log+'.out',log+'.err',directory))
                            elif arguments["--singularity"]:
                                outfile.writelines("""#!/bin/sh
#MSUB -N job_%s    #name of job
#MSUB -A ait         # sets bank account
#MSUB -l nodes=1:ppn=1,walltime=73:59:59,partition=borax  # uses 1 node
#MSUB -q pbatch         #pool
#MSUB -o %s
#MSUB -e %s
#MSUB -d %s  # directory to run from
#MSUB -V
#MSUB                     # no more psub commands



"""%(job_line,log+'.out',log+'.err',directory))
                            else:
                                outfile.writelines("""#!/bin/sh
#MSUB -N job_%s    #name of job
#MSUB -A ait         # sets bank account
#MSUB -l nodes=1:ppn=1,walltime=73:59:59,partition=borax  # uses 1 node
#MSUB -q pbatch         #pool
#MSUB -o %s
#MSUB -e %s
#MSUB -d %s  # directory to run from
#MSUB -V
#MSUB                     # no more psub commands

source %s/bin/thisroot.sh
source %s/../../../bin/geant4.sh
source %s/geant4make.sh
source %s/env.sh
source %s/env_wm.sh\n\n"""%(job_line,log+'.out',log+'.err',directory,\
                            rootDir,g4Dir,g4Dir,ratDir,watchmakersDir))
                            for _j in range(10):
                                if i*10+_j < N:
                                    mac = "macro%s/%s/%s/%s/%s/run%08d/run_%08d.mac"%(additionalMacStr,_cover,_loc,_element,_p,i*10,i*10+_j)
                                    r_outfile = "root_files%s/%s/%s/%s/%s/run%08d/run_%08d.root"%(additionalMacStr,_cover,_loc,_element,_p,i*10,i*10+_j)
                                    l_outfile = "log%s/%s/%s/%s/%s/run%08d/run_%08d.log"%(additionalMacStr,_cover,_loc,_element,_p,i*10,i*10+_j)
                                    b_outfile = "bonsai_root_files%s/%s/%s/%s/%s/run%08d/run_%08d.root"%(additionalMacStr,_cover,_loc,_element,_p,i*10,i*10+_j)
                                    if arguments["--docker"]:
                                        lines = '''drat %s %s %s
(dbonsai %s %s || dbonsai %s %s ||dbonsai %s %s ||dbonsai %s %s || echo \"Could not run bonsai after 4 tries.\")>> %s\n\n'''%(mac,r_outfile,l_outfile,r_outfile,b_outfile,r_outfile,b_outfile,r_outfile,b_outfile,r_outfile,b_outfile,l_outfile)
                                    elif arguments["--singularity"]:
                                        lines = '''%s %s -o %s -l %s
(%s %s %s || %s %s %s ||%s %s %s ||%s %s %s || echo \"Could not run bonsai after 4 tries.\")>> %s\n\n'''%(srat,mac,r_outfile,l_outfile,sbonsai,r_outfile,b_outfile,sbonsai,r_outfile,b_outfile,sbonsai,r_outfile,b_outfile,sbonsai,r_outfile,b_outfile,l_outfile)
                                    else:
                                        lines = '''rat %s -o %s -l %s
(bonsai %s %s || bonsai %s %s ||bonsai %s %s ||bonsai %s %s || echo \"Could not run bonsai after 4 tries.\")>> %s\n\n'''%(mac,r_outfile,l_outfile,r_outfile,b_outfile,r_outfile,b_outfile,r_outfile,b_outfile,r_outfile,b_outfile,l_outfile)
                                    outfile.writelines(lines)
                            if i*10+10  < N:
                                outfile.writelines("(msub %s || ./%s)"%(dir+'/job%08d.sh'%((i+1)*10),dir+'/job%08d.sh'%((i+1)*10)))
                            outfile.close()





    return 0



def deleteAllWorkDirectories():
    if arguments['--newVers']:
        d,proc,coverage = loadSimulationParameters()
    else:
        d,iso,loc,coverage,coveragePCT = loadSimulationParameters()

    dir = "log"
    if os.path.exists(dir):
        rmtree(dir)

    dir = "jobs"
    if os.path.exists(dir):
        rmtree(dir)

    for idx,cover in enumerate(coverage):
        dir = "macro"
        if os.path.exists(dir):
            rmtree(dir)

    if os.path.exists('fit_param.dat'):
        os.remove('fit_param.dat')

    if os.path.exists('like.bin'):
        os.remove('like.bin')

    if os.path.exists('sub_jobs'):
        os.remove('sub_jobs')



def testEnabledCondition(arguments):

    additionalString      = ""
    additionalCommands    = ""

    additionalMacStr = ""
    additionalMacOpt = ""

    # Commands required for root_file


    if (arguments['--detectMedia']):
        additionalMacOpt += "/rat/db/set GEO[detector] material \"%s\"\n" %(arguments['--detectMedia'])
        additionalMacStr += "_detectorMedia_%s" %(arguments['--detectMedia'])
        additionalString += "_detectorMedia_%s" %(arguments['--detectMedia'])

    if (arguments['--collectionEff']):
        additionalMacOpt += "/rat/db/set GEO[inner_pmts] efficiency_correction %f\n" %(float(arguments['--collectionEff']))
        additionalMacStr += "_collectionEfficiency_%f" %(float(arguments['--collectionEff']))
        additionalString += "_collectionEfficiency_%f" %(float(arguments['--collectionEff']))

    if (arguments['--pmtModel']):
        additionalMacOpt += "/rat/db/set GEO[inner_pmts] pmt_model \"%s\"\n" %((arguments['--pmtModel']))
        additionalMacStr += "_pmtModel_%s" %((arguments['--pmtModel']))
        additionalString += "_pmtModel_%s" %((arguments['--pmtModel']))


    if (arguments['--ipc']):
        additionalMacOpt += "/rat/db/set WATCHMAN_PARAMS photocathode_coverage %s \n" %(arguments['--ipc'])
        additionalMacStr += "_inner_coverage_%s" %(arguments['--ipc'])
        additionalString += "_inner_coverage_%s" %(arguments['--ipc'])

    if (arguments['--vpc']):
        additionalMacOpt += "/rat/db/set WATCHMAN_PARAMS veto_coverage %s \n" %(arguments['--vpc'])
        additionalMacStr += "_veto_coverage_%s" %(arguments['--vpc'])
        additionalString += "_veto_coverage_%s" %(arguments['--vpc'])

    if (arguments['--vetoModel']):
        additionalMacOpt += "/rat/db/set GEO[veto_pmts] pmt_model \"%s\"\n" %((arguments['--vetoModel']))
        additionalMacStr += "_vetoModel_%s" %((arguments['--vetoModel']))
        additionalString += "_vetoModel_%s" %((arguments['--vetoModel']))


    baseValue = 7

    if float(arguments['--tankRadius']) != defaultValues[3]:
        additionalMacOpt += "/rat/db/set GEO[tank] r_max %f\n" %(float(arguments['--tankRadius']))
        additionalMacOpt += "/rat/db/set GEO[detector] r_max %f\n" %(float(arguments['--tankRadius'])-1.5875)
        additionalMacOpt += "/rat/db/set GEO[shield] detector_size_d %f\n" %(float(arguments['--tankRadius'])*2)
        additionalMacStr += "_tankRadius_%f" %(float(arguments['--tankRadius']))
        additionalString += "_tankRadius_%f" %(float(arguments['--tankRadius']))
    else:
        additionalMacOpt += "/rat/db/set GEO[tank] r_max %f\n" %(float(defaultValues[3]))
        additionalMacOpt += "/rat/db/set GEO[detector] r_max %f\n" %(float(defaultValues[3])-1.5875)
        additionalMacOpt += "/rat/db/set GEO[shield] detector_size_d %f\n" %(float(defaultValues[3])*2)


    if float(arguments['--halfHeight'])!= defaultValues[4]:
        additionalMacOpt += "/rat/db/set GEO[tank] size_z %f\n" %(float(arguments['--halfHeight']))
        additionalMacOpt += "/rat/db/set GEO[shield] detector_size_z %f\n" %(float(arguments['--halfHeight'])*2)
        additionalMacOpt += "/rat/db/set GEO[detector] size_z %f\n" %(float(arguments['--halfHeight'])-1.5875)
        additionalMacOpt += "/rat/db/set GEO[cables] size_z %f\n" %(float(arguments['--halfHeight'])-1.5875)
        additionalMacStr += "_halfHeight_%f" %(float(arguments['--halfHeight']))
        additionalString += "_halfHeight_%f" %(float(arguments['--halfHeight']))
    else:
        additionalMacOpt += "/rat/db/set GEO[tank] size_z %f\n" %(float(defaultValues[4]))
        additionalMacOpt += "/rat/db/set GEO[shield] detector_size_z %f\n" %(float(defaultValues[4])*2)
        additionalMacOpt += "/rat/db/set GEO[detector] size_z %f\n" %(float(defaultValues[4])-1.5875)
        additionalMacOpt += "/rat/db/set GEO[cables] size_z %f\n" %(float(defaultValues[4])-1.5875)

    if float(arguments['--vetoThickR'])>0:
        additionalMacOpt += "/rat/db/set GEO[shield] veto_thickness_r %f\n" %(float(arguments['--vetoThickR']))
        additionalMacStr += "_vetoThickR_%f" %(float(arguments['--vetoThickR']))
        additionalString += "_vetoThickR_%f" %(float(arguments['--vetoThickR']))

    if float(arguments['--vetoThickZ'])>0:
        additionalMacOpt += "/rat/db/set GEO[shield] veto_thickness_z %f\n" %(float(arguments['--vetoThickZ']))
        additionalMacStr += "_vetoThickZ_%f" %(float(arguments['--vetoThickZ']))
        additionalString += "_vetoThickZ_%f" %(float(arguments['--vetoThickZ']))


    if float(arguments['--steelThick'])!= defaultValues[6]:
        additionalMacOpt += "/rat/db/set GEO[shield] steel_thickness %f\n" %(float(arguments['--steelThick']))
        additionalMacStr += "_steelThickness_%f" %(float(arguments['--steelThick']))
        additionalString += "_steelThickness_%f" %(float(arguments['--steelThick']))
    else:
        additionalMacOpt += "/rat/db/set GEO[shield] steel_thickness %f\n" %(float(defaultValues[6]))


    if float(arguments['--fidThick'])!= defaultValues[7]:
        additionalString += "_fidThickness_%f" %(float(arguments['--fidThick']))
        additionalMacStr += "_fidThickness_%f" %(float (arguments['--fidThick']))
        additionalCommands +=" --fidThick %f" %(float(arguments['--fidThick']))

    if arguments['--pmtCtrPoint']:
        additionalMacOpt += '/rat/db/set GEO[inner_pmts] orientation "point"\n'
        additionalMacOpt += '/rat/db/set GEO[shield] orientation_inner "point"\n'
        additionalMacStr += "_pmtCtrPoint_"
        additionalString += "_pmtCtrPoint_"

    if float(arguments['--U238_PPM'])!= defaultValues[8]:
        additionalString += "_U238_PPM_%f" %(float(arguments['--U238_PPM']))
        additionalCommands +=" --U238_PPM %f" %(float(arguments['--U238_PPM']))

    if float(arguments['--Th232_PPM'])!= defaultValues[9]:
        additionalString += "_Th232_PPM_%f" %(float(arguments['--Th232_PPM']))
        additionalCommands +=" --Th232_PPM %f" %(float(arguments['--Th232_PPM']))

    if float(arguments['--Rn222'])!= defaultValues[10]:
        additionalString += "_Rn222_%f" %(float(arguments['--Rn222']))
        additionalCommands +=" --Rn222 %f" %(float(arguments['--Rn222']))

    if additionalString == "":
        additionalString = "_default"

    if additionalMacStr =="":
        additionalMacStr = "_default"

    return  additionalString,additionalCommands,additionalMacStr,additionalMacOpt


def mergeNtupleFilesNew(arguments):
    # Read external requirements
    #arguments = docopt.docopt(docstring)
    # Load internal requirements
    d,proc,coverage = loadSimulationParametersNew()

    trees = {}

    additionalString,additionalCommands,additionalMacStr,additionalMacOpt = testEnabledCondition(arguments)

    directory   = os.getcwd()

    pathFinal = "bonsai_root_files%s/" %(additionalString)


    cnt = 0



    N = int(arguments['-N'])

    for _p in proc:
        for _loc in proc[_p]:
            for idx,_cover in enumerate(coverage):
                for _element in d[_p]:
                    cond1 = ('pct' in _cover)==0 and (_loc == 'GUNITE' or _loc == 'CONCRETE')!=1
                    cond2 = ('pct' in _cover)==1
                    if cond1 or cond2:
                        _tmp =  "%s_%s_%s_%s"%(_cover,_loc,_element,_p)
                        trees[_tmp] = TChain("data")
                        trees[_tmp+'_RS'] = TChain("runSummary")
                        print(_tmp)
                        for _ii in range(N):# Covers up to 1000 jobs,
                            i = int(_ii/10)
                            dir = "%s/bonsai_root_files%s/%s/%s/%s/%s/run%08d/run_%08d.root"%(directory,additionalMacStr,_cover,_loc,_element,_p,i*10,_ii)
                            try:
                                if os.path.exists(dir):
                                    _ff = TFile(dir)
                                    _data = _ff.Get('data')
                                    _tot = _data.GetEntries()
                                    # There is a tree, and there at least 0 events in it...
                                    if _tot>=0:
                                        trees[_tmp].Add(dir)
                                        trees[_tmp+'_RS'].Add(dir)
                                    else:
                                        print('No tree in file',dir)
                                _ff.Close()
                            except:
                                print('Could not read ',dir)
                        data = gROOT.FindObject('data')
                        fLocation = "%s/bonsai_root_files%s/%s/merged_%s_%s_%s.root"%(directory,additionalMacStr,_cover,_loc,_element,_p)
                        fLocationSum = "%s/bonsai_root_files%s/%s/mergedSumary_%s_%s_%s.root"%(directory,additionalMacStr,_cover,_loc,_element,_p)
                        nEntry = data.GetEntries()
                        runSummary = gROOT.FindObject('runSummary')
                        totalEntries = runSummary.GetEntries()
                        runSummary.GetEntry(0)
                        try:
                            totEvents = totalEntries*runSummary.nEvents
                        except:
                            print('Something went wrong')
                            totEvents = -1
                        dir = "%s/bonsai_root_files%s/%s/%s/%s/%s/run********/run_********.root"%(directory,additionalMacStr,_cover,_loc,_element,_p)
                        print('Merging :\n\t',dir, '\n->\n\t', fLocation, \
                        ';\ntotal entries (trigger/total):', nEntry,'/',totEvents,\
                        ',merged a total of ',totalEntries,'files.')
                        _f = TFile(fLocation,"recreate")
                        data.Write()
                        runSummary.Write()
                        _f.Close()

                        print('done\n')
                        trees[_tmp].Delete()
                        trees[_tmp+'_RS'].Delete()


    return 0
