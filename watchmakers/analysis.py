from load import *
from io_operations import testEnabledCondition,writeResultsToFile

def drange(start, stop, step):
    rii= start
    while rii<stop:
        yield rii
        rii+=step

def logx_logy_array(nbins = 500,xmin = 1e-2,xmax = 30.,ymin = 1e-9,ymax = 1e3):
    #x-axis
    logxmin = log10(xmin)
    logxmax = log10(xmax)
    xbinwidth= (logxmax-logxmin)/float(nbins)
    xbins	= zeros(nbins,dtype=float)
    xbins[0] = xmin
    #y-axis
    logymin = log10(ymin)
    logymax = log10(ymax)
    ybinwidth= (logymax-logymin)/float(nbins)
    ybins	= zeros(nbins,dtype=float)
    ybins[0] = ymin
    for i in range(nbins):
        xbins[i] = xmin + pow(10,logxmin+i*xbinwidth)
        ybins[i] = ymin + pow(10,logymin+i*ybinwidth)
    return nbins,xbins,ybins


def obtainEventEfficiency(cover,file,_tag,_distance2pmt=1,_n9=8,_dist=30.0,\
_posGood=0.1,_dirGood=0.1,_pe=8,_nhit=8,_itr = 1.5):
    # covPCT  = coveragePCT[cover]
    para = testEnabledCondition(arguments)
    additionalString  = para[0]
    arbre = {}
    arbre["rfile"] = TFile(file)
    print 'Reading', file
    try:
        runSummary = arbre["rfile"].Get('runSummary')
        Entries = runSummary.GetEntries()
        runSummary.GetEntry(Entries-1)
        events = 0
        _eventPerRun = runSummary.nEvents
    except:
        print 'File',file,'did not have run associated with it. Returning empty histogram.'
        binR,rangeRmin,rangeRmax = 31,0.45,3.55
        binwidthR = (rangeRmax-rangeRmin)/binR
        binN,rangeNmin,rangeNmax = 48,7.5,55.5
        binwidthN = (rangeNmax-rangeNmin)/binN
        h = TH2D('hist%s'%(_tag),'EMPTY - Rate of events -  %s '%(_tag),binR,rangeRmin,rangeRmax,binN,rangeNmin,rangeNmax)
        h.SetXTitle('distance from wall [m]')
        h.SetYTitle('n9 cut')
        h.SetZTitle('efficiency')
        h.GetZaxis().SetTitleOffset(-.55);
        h.GetZaxis().SetTitleColor(1);
        h.GetZaxis().CenterTitle();
        h.SaveAs("bonsai_root_files%s/%s/hist%s.C"%(additionalString,cover,_tag))
        return -1

    for i in range(10):
        events+= runSummary.subEventTally[i]
    totalEvents = float(Entries)*_eventPerRun

    arbre["data"]   = arbre["rfile"].Get('data')
    _someEntries = arbre["data"].GetEntries()

    binR,rangeRmin,rangeRmax = 31,0.45,3.55
    binwidthR = (rangeRmax-rangeRmin)/binR
    binN,rangeNmin,rangeNmax = 48,7.5,55.5
    binwidthN = (rangeNmax-rangeNmin)/binN
    h = TH2D('hist%s'%(_tag),'Rate of events -  %s '%(_tag),binR,rangeRmin,rangeRmax,binN,rangeNmin,rangeNmax)
    h.SetXTitle('distance from wall [m]')
    h.SetYTitle('n9 cut')
    h.SetZTitle('efficiency')
    h.GetZaxis().SetTitleOffset(-.55);
    h.GetZaxis().SetTitleColor(1);
    h.GetZaxis().CenterTitle();
    for _d in drange(rangeRmin+binwidthR/2.,rangeRmax,binwidthR):
        minAchieve = 0

        print '\nD:',_d
        # h['hist%s'%(_tag)].Fill(_d,rangeNmin+binwidthN/2.0,eff)
        for _n9 in range(int(rangeNmin+binwidthN/2.0),int(rangeNmax)):
            cond = "closestPMT/1000.>%f"%(_d)
            cond += "&& good_pos>%f " %(_posGood)
            cond += "&& n9 > %f && nhit > %f && pe > %f" %(_n9,_nhit,_pe)
            cond += "&& pe/nhit < %f" %(_itr)
            cond += "&& sqrt(pow(x-mcx,2)+pow(y-mcy,2)+pow(z-mcz,2))/1000.<%f"%(_dist)

            if _someEntries !=0:
                if minAchieve == 0:
                    # _evts,eff,minR,tot = obtainEventEfficiency(_cov,_file,_distance2pmt=_d,_n9=_n)
                    # print cond
                    # print arbre["data"]

                    evts = arbre["data"].Draw("",cond,"goff")
                    eff = evts/totalEvents
                    if eff == 0:
                        eff = 1/totalEvents
                        minAchieve = 1
                    h.Fill(_d,_n9,eff)
                else:
                    h.Fill(_d,_n9,1./totalEvents)
                    eff = 1./totalEvents
            else:
                h.Fill(_d,_n9,1./totalEvents)
                eff =  1./totalEvents
            print '(%2d,%4.2e),'%(_n9,eff),

    h.SaveAs("bonsai_root_files%s/%s/hist%s.C"%(additionalString,cover,_tag))
    arbre["rfile"].Close()
    del arbre
    return eff


def integralCoincidence(R,lowerBound,upperBound):
    low = -exp(-lowerBound) * (1+lowerBound )
    up  = -exp(-upperBound) * (1+upperBound)
    return up - low


def histIntegral(s,f,cut):
    # H = {}
    Hs = f.Get(s)

    if Hs!= None:
        a = Hs.Integral(0,int(cut*10))
        N = Hs.Integral(0,5000)
    else:
        a = 1.0
        N = 1.0

    # print 'histInt',cut,s,Hs,a,N
    if N !=0:
        return (1.0 - a/N )
    else:
        return 0
