#import watchmakers as PR

from watchmakers.load import *
from watchmakers.analysis import *

from ROOT import kOrange as kO,kBlue as kB,kGreen as kG
from ROOT import kMagenta as kM,kAzure as kA,kRed as kR
from ROOT import TCanvas,TH2D, gStyle,TFile,gROOT
#from .NeutrinoOscillation import nuOsc


from .io_operations import testEnabledCondition
if arguments['--noRoot']:
    print('Not loading the neutrinoOscillation module')
else:
    print('With root')    
import watchmakers.NeutrinoOscillation as nuOsc

from decimal import *
setcontext(ExtendedContext)

from numpy import max
t = arguments['--timeScale']

fidRadius = float(arguments['--tankRadius'])-float(arguments['--steelThick'])-float(arguments['--vetoThickR'])-float(arguments['--fidThick'])
fidHeight = float(arguments['--halfHeight'])-float(arguments['--steelThick'])-float(arguments['--vetoThickZ'])-float(arguments['--fidThick'])

pmtRadius = float(arguments['--tankRadius'])-float(arguments['--steelThick'])-float(arguments['--vetoThickR'])
pmtHeight = float(arguments['--halfHeight'])-float(arguments['--steelThick'])-float(arguments['--vetoThickZ'])

detectorRadius  = float(arguments['--tankRadius'])-float(arguments['--steelThick'])
detectorHeight  = float(arguments['--halfHeight'])-float(arguments['--steelThick'])

def drange(start, stop, step):
    rii= start
    while rii<stop:
        yield rii
        rii+=step

def sensitivityMapPass2New():

    site = arguments["--site"]

    # Need to fix this for future running

    OnOffRatio = float(arguments["--OnOff"])
    print(site,'with on-off ratio of ',OnOffRatio)

    cores = int(arguments["--cores"])

    if arguments["--RNRedux"]:
        rnRedux = float(arguments["--RNRedux"])
        if rnRedux>1:
            print("Value of reduction of radionuclide greater than 1, setting to 0")
            rnRedux = 0.0
    else:
        rnRedux = 0.0

    if t == 'sec':
        timeAdjustment = 24*3600.
    if t == 'day':
        timeAdjustment = 1.0
    if t == 'month':
        timeAdjustment = 1./31.
    if t == 'year':
        timeAdjustment = 1./365.
    maxTime = 14400.*timeAdjustment

    print('\nEvaluation based on geoneutrinos.org')
    #parameters  = loadAnalysisParametersNew(t)
    rates       = {}#parameters[11]
    rates["boulby_S"]=1.0
    rates["imb_S"]=1.0
    print('Wrong rates for now')
    sizeDetc    = 2.*pi*pow(fidRadius/1000.,2)*fidHeight/1000./1000.
    sizeTank    = 2.*pi*pow(detectorRadius/1000.,2)*detectorHeight/1000./1000.
    FVkTonRatio = (pow(fidRadius,2)*fidHeight)/(pow(detectorRadius,2)*detectorHeight)
    boulbyRate,imbRate = rates["boulby_S"]*FVkTonRatio,rates["imb_S"]*FVkTonRatio
    print(' boulby rates: %4.2e per %s per %4.2f kton; [r: %4.2f m; z: %4.2f m]'\
    %(boulbyRate,t,sizeDetc,fidRadius/1000.,fidHeight/1000.))
    #fast neutrons
    # print 'Debug',rates["boulby_S"],FVkTonRatio
    print('\nEvaluation not based on geoneutrinos.org')
    detectorMedium,detectorMass,reactorPower,reactorStandoff = 1,sizeDetc*1000.,1.5,24.98
    experiment = nuOsc.NeutrinoOscillation(detectorMedium,detectorMass,reactorPower,reactorStandoff)
    preOsc,afterOsc = experiment.FindRate()
    print(' Neutrino rate pre osc: %4.2f; neutrino rate post osc: %4.2f at %4.2f GWth, at %4.2f km, for %4.2f kton' %(preOsc,afterOsc,reactorPower,reactorStandoff,detectorMass/1000.))
    detectorMedium,detectorMass,reactorPower,reactorStandoff = 1,sizeDetc*1000.,1.575,24.98
    experiment = nuOsc.NeutrinoOscillation(detectorMedium,detectorMass,reactorPower,reactorStandoff)
    preOsc,afterOsc = experiment.FindRate()
    print(' Neutrino rate pre osc: %4.2f; neutrino rate post osc: %4.2f at %4.2f GWth, at %4.2f km, for %4.2f kton' %(preOsc,afterOsc,reactorPower,reactorStandoff,detectorMass/1000.))
    print('')

    proc,loca,type,color,lineS,acc,scale   = [],[],[],[],[],[],[]

    proc        += ['QGSP_BERT_EMV','QGSP_BERT_EMX','QGSP_BERT','QGSP_BIC',\
    'QBBC','QBBC_EMZ','FTFP_BERT','QGSP_FTFP_BERT']
    _t          = 'FN%s' % (site)
    loca        += [_t,_t, _t,_t,_t,_t,_t,_t]
    type        += ['si','si','si','si','si','si','si','si']
    acc         += ['corr','corr','corr','corr','corr','corr','corr', 'corr']
    color       += [kO+6,kO+6,kO+6,kO+6,kO+6,kO+6,kO+6,kO+6]
    lineS       += [2,2,2,2,2,2,2,2]
    scale       += [1./8.,1./8.,1./8.,1./8.,1./8.,1./8.,1./8.,1./8.]

    #radionuclides
    proc        += ['9003','11003']
    _t          =  'RN%s' % (site)
    loca        += [_t,_t]
    type        += ['ei','ei']
    acc         += ['di','di']
    color       += [kG+3,kG+3]
    lineS       += [1,1]
    scale       += [1.,1.]
    #ibds
    if site == 'boulby':
        proc    += ['boulby','boulby','neutron']
    else:
        proc    += ['imb','imb','neutron']
    loca        += ['S','S','N%s'%(site)]
    type        += ['ei','ei','ei']
    acc         += ['di','corr','corr']
    color       += [kA+0,kA-0,kA-0]
    lineS       += [1,2,2]
    scale       += [-1.0,0.0,0.0]

    c1 = TCanvas('c1','c1',1618,1000)
    c1.SetRightMargin(0.53)
    c1.Divide(2,2)

    additionalString,additionalCommands,additionalMacStr,additionalMacOpt = testEnabledCondition(arguments)
    if additionalString == "":
        additionalString = "_default"
    location = 'Boulby '

    hist = TH2D('hist','3#sigma discovery phase space -  %s '%(location),31,9.5,40.5,30,9.5,39.5)
    hist.SetXTitle('photocoverage [%]')
    hist.SetYTitle('photoelectron threhsold cut [p.e.]')
    hist.SetZTitle('off-time [%s] for 3 #sigma discovery'%(t))
    hist.GetZaxis().SetTitleOffset(-.55);
    hist.GetZaxis().SetTitleColor(1);
    hist.GetZaxis().CenterTitle();
    gStyle.SetOptStat(0)
    gStyle.SetPalette(55)

    h = {}

    binR,rangeRmin,rangeRmax = 31,0.45,3.55
    binwidthR = (rangeRmax-rangeRmin)/binR
    binN,rangeNmin,rangeNmax = 48,7.5,55.5
    binwidthN = (rangeNmax-rangeNmin)/binN

    _cov = arguments['-C']

    d,proc,coverage = loadSimulationParametersNew()

    minAchieve = 0
    
    for idx,_cover in enumerate(coverage):
        _str = "bonsai_root_files%s/%s/results.root"%(additionalMacStr,_cover)
        outfile = TFile(_str,"RECREATE")
        for _p in proc:
            for _loc in proc[_p]:
                for _element in d[_p]:
                    _tag = "%s_%s_%s_%s"%(_cover,_loc,_element,_p)
                    _file = "bonsai_root_files%s/%s/merged_%s_%s_%s.root"%(additionalMacStr,_cover,_loc,_element,_p)
                    print(_tag)
                    obtainEventEfficiency(_cover,_file,_tag,outfile)

                    print('')
        print('Saving outfile:',_str)
        outfile.Close()
    return 0

def EfficiencyMapInPMTVol():
    '''For the given configuration on initializing watchmakers,
    Generates PMT Volume Efficiency histograms for all merged
    files.'''
    n9 = int(arguments['--minN9'])
    good_pos = float(arguments['-g'])
    good_dir = float(arguments['-G'])
    print("Evaluating sensitivity in PMT Volume for all WaterVolume types. Using given minimum parameter requirements")
    additionalString,additionalCommands,additionalMacStr,additionalMacOpt = testEnabledCondition(arguments)
    if additionalString == "":
        additionalString = "_default"
    d,proc,coverage = loadSimulationParametersNew()
    for _p in proc:
        for _loc in proc[_p]:
            if _loc != "WaterVolume":
                continue
            for idx,_cover in enumerate(coverage):
                for _element in d[_p]:
                    _tag = "%s_%s_%s_%s"%(_cover,_loc,_element,_p)
                    _file = "bonsai_root_files%s/%s/merged_%s_%s_%s.root"%(additionalMacStr,_cover,_loc,_element,_p)
                    print(_tag)
                    obtainEfficiencyInPMTVol(_cover,_file,_tag, _n9=n9, _posGood=good_pos, _dirGood=good_dir)
                    print('')




def readEfficiencyHistogram():

    hist = {}

    d,proc,coverage = loadSimulationParametersNew()
    additionalString,additionalCommands,additionalMacStr,additionalMacOpt = testEnabledCondition(arguments)

    print('Reading in root tree')
    for idx,_cover in enumerate(coverage):
        _str = "bonsai_root_files%s/%s/results.root"%(additionalMacStr,_cover)
        outfile = TFile(_str,"READ")
        for _p in proc:
            for _loc in proc[_p]:
                for _element in d[_p]:
                    _tag = "hist%s_%s_%s_%s"%(_cover,_loc,_element,_p)
                    hist[_tag] =  outfile.Get(_tag)
                    print(hist[_tag].GetEntries())

    print('done processing result file. What is in the directory:')
    gROOT.ProcessLine('.ls')

    print('\nLoading PMT activity:')
    mPMTs,mPMTsU238,mPMTsTh232,mPMTsK40 = loadPMTActivity()
    print('done.')

    print('\nLoading Veto activity:')
    mVETOs,mVETOsU238,mVETOsTh232,mVETOsK40 = loadVETOActivity()
    print('done.')

    print('\nLoading Gd activity:')
    GdU238,GdTh232,GdU235,GdU238_l,GdTh232_l,GdU235_l = loadGdActivity()
    print('done.')

    tankRadius  = float(arguments["--tankRadius"])-float(arguments['--steelThick'])
    tankHeight  = float(arguments["--halfHeight"])-float(arguments['--steelThick'])
    nKiloTons = pi*pow(tankRadius/1000.,2)*(2.*tankHeight/1000.)
    rRn222 = float(arguments["--Rn222"])*nKiloTons
    print('\nLoaded Rn-222 activity of ',rRn222,'Bq per water volume, assumed a rate of %4.3e Bq/m^3'%(float(arguments["--Rn222"])))


    FreeProtons = 0.668559
    TNU         = FreeProtons* nKiloTons /1000./365./24./3600.
    boulbyIBDRate   = 800.*TNU
    print('\nLoaded an IBD rate of ',boulbyIBDRate,' events per water volume per second, assumed a rate of %4.3e TNU'%(boulbyIBDRate/TNU))

    innerRad = 12.5 #meters
    outerRad = 13.5 #meters
    rockMass = (pi*pow(outerRad,2)*(2.*outerRad)-pi*pow(innerRad,2)*(2.*innerRad))*power(100.,3)*2.39
    # volumeR         = power(22.,3)-power(20.,3)# Rock cavern e.g. (22m x 22m x 22m) - (20m x 20m x 20m)
    # density         = 2.39 #from McGrath
    # rockMass        = volumeR*power(100.,3)*density
    #Mass of rock evalyated
    avgMuon         = npa([180.,264.])
    avgMuonNC       = power(avgMuon,0.849)
    avgNFluxMag     = 1e-6
    muonRate        = npa([7.06e-7,4.09e-8]) # mu/cm2/s
    tenMeVRatio     = npa([7.51/34.1,1.11/4.86])
    fastNeutrons    = rockMass*avgMuonNC*avgNFluxMag*muonRate*tenMeVRatio
    FN_boulby = fastNeutrons[1]

    avgRNYieldRC    = power(avgMuon,0.73)
    skRNRate        = 0.5e-7 # 1/mu/g cm2
    avgMuonSK       = power(219.,0.73)
    skMuFlux        = 1.58e-7 #mu/cm2/sec
    radionuclideRate= (skRNRate*avgRNYieldRC/avgMuonSK)*muonRate*nKiloTons*1e9
    RN_boulby        = radionuclideRate[1]
    print('\nLoaded mass of rock %e g. Fast Neutron Yield %e per sec; radionuclide yield %e per sec'%(rockMass,FN_boulby,RN_boulby))


    print('\n What are the maximum efficiency/rate found in each histogram:')
    _sing = 0.0
    lineU238PMT,lineTh232PMT,lineKPMT = '','',''
    lineU238VETO,lineTh232VETO,lineKVETO = '','',''
    lineFNROCK = ''
    lineU238GUN,lineTh232GUN,lineKGUN = '','',''
    lineU238ROCK,lineTh232ROCK,lineKROCK = '','',''
    lineU238CONC,lineTh232CONC,lineKCONC = '','',''
    lineRn222WaterVolume = ''
    lineU238GD,lineU235GD,lineTH232GD = '','',''
    linePromptWaterVolume,lineDelayedWaterVolume,linePromptDelayedWaterVolume = '','',''
    lineELSE = ''

    firstGo = 1
    for _t in hist:
        if firstGo:
            h = hist[_t].Clone()
            h.SetZTitle('singles rate (Hz)')
            h.SetTitle('Singles rate')
            h.SetName('hSinglesRate')
            h.Reset()
            hn = hist[_t].Clone()
            hn.SetZTitle('neutron rate (Hz)')
            hn.SetTitle('rate')
            hn.SetName('hNeutronRate')
            hn.Reset()
            hee = hist[_t].Clone()
            hee.SetZTitle('positron efficiency')
            hee.SetTitle('efficiency')
            hee.SetName('hPositronEfficiency')
            hee.Reset()
            hne = hist[_t].Clone()
            hne.SetZTitle('neutron efficiency')
            hne.SetTitle('efficiency')
            hne.SetName('hNeutronEfficiency')
            hne.Reset()

            firstGo =0
        if 'PMT' in _t and 'CHAIN_238U_NA' in _t:
            if '210Tl' in _t:
                _sing+=hist[_t].GetMaximum()*mPMTsU238[0]*0.002
                h.Add(hist[_t],mPMTsU238[0]*0.002)
                lineU238PMT += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*mPMTsU238[0]*0.002)
            else:
                _sing+=hist[_t].GetMaximum()*mPMTsU238[0]
                h.Add(hist[_t],mPMTsU238[0])
                lineU238PMT+= "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*mPMTsU238[0])
        elif 'PMT' in _t and 'CHAIN_232Th_NA' in _t:
            _sing+=hist[_t].GetMaximum()*mPMTsTh232[0]
            h.Add(hist[_t],mPMTsTh232[0])
            lineTh232PMT += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*mPMTsTh232[0])
        elif 'PMT' in _t and '40K_NA' in _t:
            _sing+=hist[_t].GetMaximum()*mPMTsK40[0]
            h.Add(hist[_t],mPMTsK40[0])
            lineKPMT += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*mPMTsK40[0])

        elif 'VETO' in _t and 'CHAIN_238U_NA' in _t:
            if '210Tl' in _t:
                _sing+=hist[_t].GetMaximum()*mVETOsU238*0.002
                h.Add(hist[_t],mVETOsU238*0.002)
                lineU238VETO += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*mVETOsU238*0.002)
            else:
                _sing+=hist[_t].GetMaximum()*mVETOsU238
                h.Add(hist[_t],mVETOsU238)
                lineU238VETO+= "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*mVETOsU238)
        elif 'VETO' in _t and 'CHAIN_232Th_NA' in _t:
            print(mVETOsTh232)
            _sing+=hist[_t].GetMaximum()*mVETOsTh232
            h.Add(hist[_t],mVETOsTh232)
            lineTh232VETO += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*mVETOsTh232)
        elif 'VETO' in _t and '40K_NA' in _t:
            _sing+=hist[_t].GetMaximum()*mVETOsK40
            h.Add(hist[_t],mVETOsK40)
            lineKVETO += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*mVETOsK40)

        elif 'GUNITE' in _t and 'CHAIN_238U_NA' in _t:
            if '210Tl' in _t:
                _sing+=hist[_t].GetMaximum()*mPMTsU238[0]*0.002
                lineU238GUN += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*mPMTsU238[0]*0.002)
            else:
                _sing+=hist[_t].GetMaximum()*mPMTsU238[0]
                lineU238GUN+= "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*mPMTsU238[0])
        elif 'GUNITE' in _t and 'CHAIN_232Th_NA' in _t:
            _sing+=hist[_t].GetMaximum()*mPMTsTh232[0]
            lineTh232GUN += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*mPMTsTh232[0])
        elif 'GUNITE' in _t and '40K_NA' in _t:
            _sing+=hist[_t].GetMaximum()*mPMTsK40[0]
            lineKGUN += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*mPMTsK40[0])

        elif 'ROCK' in _t and 'CHAIN_238U_NA' in _t:
            if '210Tl' in _t:
                _sing+=hist[_t].GetMaximum()*mPMTsU238[0]*0.002
                lineU238ROCK += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*mPMTsU238[0]*0.002)
            else:
                _sing+=hist[_t].GetMaximum()*mPMTsU238[0]
                lineU238ROCK+= "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*mPMTsU238[0])
        elif 'ROCK' in _t and 'CHAIN_232Th_NA' in _t:
            _sing+=hist[_t].GetMaximum()*mPMTsTh232[0]
            lineTh232ROCK += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*mPMTsTh232[0])
        elif 'ROCK' in _t and '40K_NA' in _t:
            _sing+=hist[_t].GetMaximum()*mPMTsK40[0]
            lineKROCK += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*mPMTsK40[0])
        elif 'ROCK' in _t and '_FN' in _t:
            lineFNROCK += "%50s %e %15.10f (%15.10f per day)\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*FN_boulby,hist[_t].GetMaximum()*FN_boulby*3600.*24.)
        elif 'CONC' in _t and 'CHAIN_238U_NA' in _t:
            if '210Tl' in _t:
                _sing+=hist[_t].GetMaximum()*mPMTsU238[0]*0.002
                lineU238CONC += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*mPMTsU238[0]*0.002)
            else:
                _sing+=hist[_t].GetMaximum()*mPMTsU238[0]
                lineU238CONC+= "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*mPMTsU238[0])
        elif 'CONC' in _t and 'CHAIN_232Th_NA' in _t:
            _sing+=hist[_t].GetMaximum()*mPMTsTh232[0]
            lineTh232CONC += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*mPMTsTh232[0])
        elif 'CONC' in _t and '40K_NA' in _t:
            _sing+=hist[_t].GetMaximum()*mPMTsK40[0]
            lineKCONC += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*mPMTsK40[0])

        elif 'GD' in _t and 'CHAIN_238U_NA' in _t:
            if '210Tl' in _t:
                _sing+=hist[_t].GetMaximum()*GdU238_l*0.002
                h.Add(hist[_t],GdU238_l*0.002)
                lineU238GD += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*GdU238_l*0.002)
            elif '234Pa' in _t:
                _sing+=hist[_t].GetMaximum()*GdU238
                h.Add(hist[_t],GdU238)
                lineU238GD += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*GdU238)
            else:
                _sing+=hist[_t].GetMaximum()*GdU238_l
                h.Add(hist[_t],GdU238_l)
                lineU238GD += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*GdU238_l)
        elif 'GD' in _t and 'CHAIN_232Th_NA' in _t:
            if '228Ac' in _t:
                _sing+=hist[_t].GetMaximum()*GdTh232
                h.Add(hist[_t],GdTh232)
                lineTH232GD += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*GdTh232)
            else:
                _sing+=hist[_t].GetMaximum()*GdTh232_l
                h.Add(hist[_t],GdTh232_l)
                lineTH232GD += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*GdTh232_l)
        elif 'GD' in _t and 'CHAIN_235U_NA' in _t:
            if '231Th' in _t:
                _sing+=hist[_t].GetMaximum()*GdU235
                h.Add(hist[_t],GdU235)
                lineU235GD += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*GdU235)
            else:
                _sing+=hist[_t].GetMaximum()*GdU235_l
                h.Add(hist[_t],GdU235_l)
                lineU235GD += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*GdU235_l)


        elif 'WaterVolume' in _t and 'CHAIN_222Rn_NA' in _t:
            if '210Tl' in _t:
                _sing+=hist[_t].GetMaximum()*rRn222*0.002
                h.Add(hist[_t],rRn222*0.002)
                lineRn222WaterVolume += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*rRn222*0.002)
            else:
                _sing+=hist[_t].GetMaximum()*rRn222
                h.Add(hist[_t],rRn222)
                lineRn222WaterVolume+= "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*rRn222)

        elif 'WaterVolume' in _t and 'promptPositron' in _t:
            _day=hist[_t].GetMaximum()*boulbyIBDRate*3600.*24
            hee.Add(hist[_t],1.)#Only add efficiency for positron
            linePromptWaterVolume += "%50s %e %15.10f per sec (%15.10f per day)\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*boulbyIBDRate,_day)

        elif 'WaterVolume' in _t and 'delayedNeutron' in _t:
            _day=hist[_t].GetMaximum()*boulbyIBDRate*3600.*24
            hn.Add(hist[_t],boulbyIBDRate)
            hne.Add(hist[_t],1.)
            linePromptWaterVolume += "%50s %e %15.10f per sec (%15.10f per day)\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*boulbyIBDRate,_day)

        elif 'WaterVolume' in _t and 'promptDelayedPair' in _t:
            _day=hist[_t].GetMaximum()*boulbyIBDRate*3600.*24
            linePromptWaterVolume += "%50s %e %15.10f per sec (%15.10f per day)\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*boulbyIBDRate,_day)

        else:
            lineELSE += "%50s %e\n"%(_t,hist[_t].GetMaximum())


                    # hist =
                    # print _hist.GetMaximum()
    print('')
    print('PMT U-238 \n', lineU238PMT)
    print('PMT Th-232\n', lineTh232PMT)
    print('PMT K\n', lineKPMT,'\n')
    print('VETO U-238 \n', lineU238VETO)
    print('VETO Th-232 \n', lineTh232VETO)
    print('VETO K \n', lineKVETO,'\n')
    print('Water Rn-222\n',lineRn222WaterVolume,'\n')
    print('Gunite U-238 \n', lineU238GUN)
    print('Gunite Th-232\n', lineTh232GUN)
    print('Gunite K\n', lineKGUN,'\n')
    print('Concrete U-238 \n', lineU238CONC)
    print('Concrete Th-232\n', lineTh232CONC)
    print('Concrete K\n', lineKCONC,'\n')
    print('ROCK U-238 \n', lineU238ROCK)
    print('ROCK Th-232\n', lineTh232ROCK)
    print('ROCK K\n', lineKROCK,'\n')
    print('ROCK Fast Neutron\n',lineFNROCK,'\n')
    print('Else  \n', lineELSE,'\n')
    print('Total singles rate:\t\t\t',_sing,'events per sec at minimum buffer distance of 0.5 m\n')
    print('Signal information')
    print('Prompt positron Water volume \n', linePromptWaterVolume)
    signal = ['WaterVolume_delayedNeutron_ibd_n','WaterVolume_promptPositron_ibd_p']
    _str =  "bonsai_root_files%s/%s/histograms_U238_%4.3fPPM_Th232_%4.3fPPM_K_%4.3fPPM.root"%(additionalMacStr,_cover,float(arguments["--U238_PPM"]),float(arguments["--Th232_PPM"]),float(arguments["--K_PPM"]))


    timeAcc = 0.0001*86400.

    offsets_n9 = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17]  ## bin numbers
    offsets_dtw = [0]       ## bin numbers

    binR,rangeRmin,rangeRmax = 31,0.45,3.55
    binwidthR = (rangeRmax-rangeRmin)/binR
    binN,rangeNmin,rangeNmax = 48,7.5,55.5
    binwidthN = (rangeNmax-rangeNmin)/binN

    _cov = arguments['-C']


    pmtRadius = float(arguments['--tankRadius'])-float(arguments['--steelThick'])-float(arguments['--vetoThickR'])
    pmtHeight = float(arguments['--halfHeight'])-float(arguments['--steelThick'])-float(arguments['--vetoThickZ'])
    detectorRadius  = float(arguments['--tankRadius'])-float(arguments['--steelThick'])
    detectorHeight  = float(arguments['--halfHeight'])-float(arguments['--steelThick'])

    _maxSignal,_maxBkgd,_maxSoverB,_maxOffn9,_maxOff_dtw = -1,-1,-1,-1,-1
    _maxSignal2,_maxBkgd2,_maxSoverB2,_maxOffn92,_maxOff_dtw2,_maxOffset2 = -1,-1,-1,-1,-1,-1


    line,_line,_line2= ("",),"",""
    _histograms = {}
    for offset in offsets_n9:
        for fv_offset in offsets_dtw:
            _histograms["sOverB_%d"%(offset)]= h.Clone()
            _histograms["sOverB_%d"%(offset)].SetZTitle('signal/sqrt(signal+background)')
            _histograms["sOverB_%d"%(offset)].SetYTitle('n9 cut on prompt')
            _histograms["sOverB_%d"%(offset)].SetTitle('hSoB - offset %2d'%(offset))
            _histograms["sOverB_%d"%(offset)].SetName('hSoB%d'%(offset))
            _histograms["sOverB_%d"%(offset)].Reset()

            _histograms["signal_%d"%(offset)]= h.Clone()
            _histograms["signal_%d"%(offset)].SetZTitle('evts/day')
            _histograms["signal_%d"%(offset)].SetYTitle('n9 cut on prompt')
            _histograms["signal_%d"%(offset)].SetTitle('Signal - offset %2d'%(offset))
            _histograms["signal_%d"%(offset)].SetName('hSignal%d'%(offset))
            _histograms["signal_%d"%(offset)].Reset()

            _histograms["background_%d"%(offset)]= h.Clone()
            _histograms["background_%d"%(offset)].SetZTitle('evts/day')
            _histograms["background_%d"%(offset)].SetYTitle('n9 cut on prompt')
            _histograms["background_%d"%(offset)].SetTitle('accidentals - offset %2d'%(offset))
            _histograms["background_%d"%(offset)].SetName('hBackground%d'%(offset))
            _histograms["background_%d"%(offset)].Reset()



            _proc = '_%d_%d_%s'%(offset,fv_offset,_cov)

            _maxSignal,_maxBkgd,_maxSoverB,_maxOffn9,_maxOff_dtw = -1,-1,-1,-1,-1
            for _d in range(binR-fv_offset-1):
                for _n in range(binN-offset-1):
                    _db,_nb=_d+1,_n+1
                    _p_d  = hee.GetXaxis().GetBinCenter(_db)
                    _p_n9 = hee.GetYaxis().GetBinCenter(_nb)
                    _n_d  = hne.GetXaxis().GetBinCenter(_db+fv_offset)
                    _n_n9 = hne.GetYaxis().GetBinCenter(_nb+offset)
                    fidRadius,fidHeight = pmtRadius-_p_d*1e3,pmtHeight-_p_d*1e3
                    ratioScaling = (pow(fidRadius,2)*fidHeight)/(pow(detectorRadius,2)*detectorHeight)
                    _p_v  = hee.GetBinContent(_db,_nb)/ratioScaling
                    _n_v  = hne.GetBinContent(_db+fv_offset,_nb+offset)/ratioScaling
                    _rate_v  = hn.GetBinContent(_db+fv_offset,_nb+offset)

                    _signal = _rate_v*_p_v*86400.

                    _p_d_1  = h.GetXaxis().GetBinCenter(_db)
                    _p_n9_1 = h.GetYaxis().GetBinCenter(_nb)
                    _n_d_1  = h.GetXaxis().GetBinCenter(_db+fv_offset)
                    _n_n9_1 = h.GetYaxis().GetBinCenter(_nb+offset)
                    _p_v_1  = h.GetBinContent(_db,_nb)
                    _n_v_1  = h.GetBinContent(_db+fv_offset,_nb+offset)

                    _background = _p_v_1*_n_v_1*timeAcc*0.05
                    sob = _signal/sqrt(_signal+_background)

                    _histograms["sOverB_%d"%(offset)].SetBinContent(_db,_nb,sob)
                    _histograms["signal_%d"%(offset)].SetBinContent(_db,_nb,_signal)
                    _histograms["background_%d"%(offset)].SetBinContent(_db,_nb,_background)

                    if sob >_maxSoverB:
                        _maxSoverB = _signal/sqrt(_signal+_background)
                        _maxSignal = _signal
                        _maxBkgd   = _background
                        _maxOffn9   = _n_n9
                        _maxOff_dtw = _n_d
                        _line = ("Offset:%3d, beta/n: Wall-dist (%4.1f,%4.1f) m, n9 cut (%d,%d), rel. efficiency (%4.2f,%4.2f), neutron rate :(%4.2f per day), combined eff/rate : %4.2f per day;"\
                            %(offset,_p_d,_n_d\
                            ,_p_n9,_n_n9\
                            ,_p_v,_n_v\
                            ,_rate_v*86400.\
                            ,_rate_v*_p_v*86400.),)
                        _line2 =    ("acc. rate (%5.3f,%5.3f): acc. combined rate: %4.3f per day (pre-prox)"\
                         %(_p_v_1,_n_v_1,_p_v_1*_n_v_1*timeAcc),)

                    if sob >_maxSoverB2:
                        _maxSoverB2 = _signal/sqrt(_signal+_background)
                        _maxSignal2 = _signal
                        _maxBkgd2   = _background
                        _maxOffn92   = _n_n9
                        _maxOff_dtw2 = _n_d
                        _maxOffset2 = offset


            print('Offset:',str(offset).rjust(3,' '),',Found max S/sqrt(S+B)',_maxSoverB,',(S,B,n9,dtw):(',_maxSignal,_maxBkgd,_maxOffn9,_maxOff_dtw,')')
            line += (_line + _line2,)

    print('\n\nMore info on the maximal sensitivity found:')
    # print line

    for _l in line:
        for i in range(len(_l)):
            print(_l[i], end=' ')
        print('')

    cut_signal = 0.9
    sigma = 4.65 # 3-sigma 95% of the time, according to Owen report
    _S = _maxSignal2*cut_signal
    _B = _maxSignal2*cut_signal*1.15 + _maxBkgd2 # Includes a 15% other reactor component
    OnOffRatio = 1.5

    if _S >0:
        T3SIGMA = sigma**2*(_B +(_S+_B)/OnOffRatio)/_S/_S

    else:
        T3SIGMA = 1e999
    metric = T3SIGMA*OnOffRatio + T3SIGMA
    _res = "%s %4.1f %3d %3d %4.3f %4.3f %4.1f %4.1f" %(_cover,_maxOff_dtw2,_maxOffn92,_maxOffn92-_maxOffset2,_maxSignal2,_maxBkgd2,T3SIGMA,metric)
    _strRes = "results_DTW_%dmm_U238_%4.3fPPM_Th232_%4.3fPPM_K_%4.3fPPM.txt"%(float(arguments['--vetoThickR']),float(arguments["--U238_PPM"]),float(arguments["--Th232_PPM"]),float(arguments["--K_PPM"]))
    # _strRes = 'res.txt'
    with open(_strRes,'a') as file:
        file.write(_res+'\n')

    print('\n\nWriting histograms to file',_str)
    f_root = TFile(_str,"recreate")
    h.Write()
    hn.Write()
    hne.Write()
    hee.Write()
    for offset in offsets_n9:
        _histograms["sOverB_%d"%(offset)].Write()
        _histograms["signal_%d"%(offset)].Write()
        _histograms["background_%d"%(offset)].Write()
    f_root.Close()

    print('\n\nall done.')



def findRate():


    d,proc,coverage = loadSimulationParametersNew()
    additionalString,additionalCommands,additionalMacStr,additionalMacOpt = testEnabledCondition(arguments)

    print('\nLoading TANK activity:')
    tankmass,tankact_60co,tankact_137cs = loadTankActivity()
    print('done.')

    print('\nLoading CONCRETE activity:')
    concmass,concact_238u,concact_232th,concact_40k = loadConcreteActivity()
    print('done.')

    print('\nLoading Shotcrete activity:')
    shotmass,act_238u,act_232th,act_40k = loadShotcreteActivity()
    print('done.')

    print('\nLoading ROCK activity:')
    rockmass,act_238u,act_232th,act_40k = loadRockActivity()
    print('done.')


    print('\nLoading PMT activity:')
    mPMTs,mPMTsU238,mPMTsTh232,mPMTsK40 = loadPMTActivity()
    print('done.')

    print('\nLoading VETO activity:')
    mVETOs,mVETOsU238,mVETOsTh232,mVETOsK40 = loadVETOActivity()
    print('done.')


    tankRadius  = float(arguments["--tankRadius"])-float(arguments['--steelThick'])
    tankHeight  = float(arguments["--halfHeight"])-float(arguments['--steelThick'])
    nKiloTons = pi*pow(tankRadius/1000.,2)*(2.*tankHeight/1000.)
    rRn222 = float(arguments["--Rn222"])*nKiloTons
    print('\nLoaded Rn-222 activity of ',rRn222,'Bq per water volume, assumed a rate of %4.3e Bq/m^3'%(float(arguments["--Rn222"])))

    FreeProtons = 0.668559
    TNU         = FreeProtons* nKiloTons /1000./365./24./3600.
    boulbyIBDRate   = 800.*TNU
    print('\nLoaded an IBD rate of ',boulbyIBDRate,' events per water volume per second, assumed a rate of %4.3e TNU'%(boulbyIBDRate/TNU))

    innerRad = 12.5 #meters
    outerRad = 13.5 #meters
    rockMass = (pi*pow(outerRad,2)*(2.*outerRad)-pi*pow(innerRad,2)*(2.*innerRad))*power(100.,3)*2.39
    # volumeR         = power(22.,3)-power(20.,3)# Rock cavern e.g. (22m x 22m x 22m) - (20m x 20m x 20m)
    # density         = 2.39 #from McGrath
    # rockMass        = volumeR*power(100.,3)*density
    #Mass of rock evalyated
    avgMuon         = npa([180.,264.])
    avgMuonNC       = power(avgMuon,0.849)
    avgNFluxMag     = 1e-6
    muonRate        = npa([7.06e-7,4.09e-8]) # mu/cm2/s
    tenMeVRatio     = npa([7.51/34.1,1.11/4.86])
    fastNeutrons    = rockMass*avgMuonNC*avgNFluxMag*muonRate*tenMeVRatio
    FN_boulby = fastNeutrons[1]

    avgRNYieldRC    = power(avgMuon,0.73)
    skRNRate        = 0.5e-7 # 1/mu/g cm2
    avgMuonSK       = power(219.,0.73)
    skMuFlux        = 1.58e-7 #mu/cm2/sec
    radionuclideRate= (skRNRate*avgRNYieldRC/avgMuonSK)*muonRate*nKiloTons*1e9
    RN_boulby        = radionuclideRate[1]
    print('\nLoaded mass of rock %e g. Fast Neutron Yield %e per sec; radionuclide yield %e per sec'%(rockMass,FN_boulby,RN_boulby))
    _strRes = "rate_%dmm_U238_%4.3fPPM_Th232_%4.3fPPM_K_%4.3fPPM_%s.txt"%(float(arguments['--vetoThickR']),float(arguments["--U238_PPM"]),float(arguments["--Th232_PPM"]),float(arguments["--K_PPM"]),arguments["-C"])
    _res = "%e %e %d %e %e %e"%(boulbyIBDRate,rRn222,mPMTs[0],mPMTsU238[0],mPMTsTh232[0],mPMTsK40[0])
    with open(_strRes,'a') as file:
        file.write(_res+'\n')

    return boulbyIBDRate,rRn222,mPMTsU238,mPMTsTh232,mPMTsK40
