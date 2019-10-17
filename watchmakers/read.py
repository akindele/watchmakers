# #import watchmakers as PR
#
# from watchmakers.load import *
# from watchmakers.analysis import *
# from io_operations import testEnabledCondition
# import watchmakers.NeutrinoOscillation as nuOsc
# from decimal import *
# setcontext(ExtendedContext)
#
# from numpy import max
# t = arguments['--timeScale']
#
# fidRadius = float(arguments['--tankRadius'])-float(arguments['--steelThick'])-float(arguments['--vetoThickR'])-float(arguments['--fidThick'])
# fidHeight = float(arguments['--halfHeight'])-float(arguments['--steelThick'])-float(arguments['--vetoThickZ'])-float(arguments['--fidThick'])
# detectorRadius  = float(arguments['--tankRadius'])-float(arguments['--steelThick'])
# detectorHeight  = float(arguments['--halfHeight'])-float(arguments['--steelThick'])
#
# def drange(start, stop, step):
#     rii= start
#     while rii<stop:
#         yield rii
#         rii+=step
#
#
# def accidentalCanvas(Graphs):
#     # Not currently in use
#     proc    = []
#     loca    = []
#     type    = []
#     color   = []
#     lineS   = []
#     acc     = []
#
#     proc    += ['234Pa','214Pb','214Bi','210Bi','210Tl']
#     loca    += ['PMT','PMT','PMT','PMT','PMT']
#     type    += ['si','si','si','si','si']
#     color   += [kO+0,  kO+1, kO+2, kO+3, kO-2]
#
#     proc    +=['232Th','228Ac','212Pb','212Bi','208Tl']
#     loca    +=['PMT','PMT','PMT','PMT','PMT']
#     type    += ['si','si','si','si','si']
#     color   += [kM+0,  kM-6, kM+2, kM-3, kM-7]
#
#     proc    +=['214Pb','214Bi','210Bi','210Tl']
#     loca    +=['FV','FV','FV','FV']
#     type    += ['ei','ei','ei','ei']
#     color   += [kB+0,  kB+1, kB+2, kB+3]
#
#     proc    +=['214Pb','214Bi','210Bi','210Tl']
#     loca    +=['FV','FV','FV','FV']
#     type    += ['si','si','si','si']
#     color   += [kG+0,  kG+1, kG+2, kG+3]
#
#     proc    +=['QGSP_BERT_EMV','QGSP_BERT_EMX','QGSP_BERT','QGSP_BIC',\
#     'QBBC','QBBC_EMZ','FTFP_BERT','QGSP_FTFP_BERT']
#     loca    +=['FN','FN','FN','FN','FN','FN','FN','FN']
#     type    +=['si','si','si','si','si','si','si','si']
#     color   += [kO+2,kO+4,kO+2,kO+4,kO+2,kO+4,kO+2,kO+4]
#
#
# def siteCanvas(Graphs,site,cut,hist):
#     # Historical, not currently in use
#     proc        = []
#     loca        = []
#     type        = []
#     color       = []
#     lineS       = []
#     acc         = []
#     scale       = []
#
#     #fast neutrons
#     proc        += ['QGSP_BERT_EMV','QGSP_BERT_EMX','QGSP_BERT','QGSP_BIC',\
#     'QBBC','QBBC_EMZ','FTFP_BERT','QGSP_FTFP_BERT']
#     _t          = 'FN%s' % (site)
#     loca        += [_t,_t, _t,_t,_t,_t,_t,_t]
#     type        += ['si','si','si','si','si','si','si','si']
#     acc         += ['corr','corr','corr','corr','corr','corr','corr', 'corr']
#     color       += [kO+6,kO+6,kO+6,kO+6,kO+6,kO+6,kO+6,kO+6]
#     lineS       += [2,2,2,2,2,2,2,2]
#     scale       += [1./8.,1./8.,1./8.,1./8.,1./8.,1./8.,1./8.,1./8.]
#
#     #ibds
#     if site == 'boulby':
#         proc    += ['boulby','boulby','neutron']
#     else:
#         proc    += ['imb','imb','neutron']
#     loca        += ['S','S','N%s'%(site)]
#     type        += ['ei','ei','ei']
#     acc         += ['di','corr','corr']
#     color       += [kA+0,kA-0,kA-0]
#     lineS       += [1,2,2]
#     scale       += [-1.0,0.0,0.0]
#
#     proc        += ['9003','11003']
#     _t          =  'RN%s' % (site)
#     loca        += [_t,_t]
#     type        += ['ei','ei']
#     acc         += ['di','di']
#     color       += [kG+3,kG+3]
#
#
#     lineS       += [1,1,1,1,1,1,1,1]
#     scale       += [1.,1.,1.,1.,1.,1.,1.,1.]
#
#
#     for PC in range(10,41):
#         S = 0.
#         B = Graphs['all_cut_accidental'].Eval(PC)
#
#         for i,arr in enumerate(proc):
#             _str = 'scaled_%s_%s_%s_1_abs_%s' %(type[i],proc[i],loca[i],acc[i])
#             if scale[i] > 0.:
#                 B+=Graphs[_str].Eval(PC)*scale[i]
#             elif scale[i] < 0:
#                 S+=Graphs[_str].Eval(PC)
#                 #        print cut,S,B,S/B,S/sqrt(B+S)
#         hist.Fill(PC,cut,S/sqrt(B+S))
#
#     C1 = TCanvas('C1','%s'%(t),1200,800)
#
#     Graphs['all_cut_accidental'].Draw('AC')
#     Graphs['all_cut_accidental'].SetLineColor(2)
#     Graphs['all_cut_accidental'].SetLineWidth(3)
#     Graphs['all_cut_accidental'].SetTitle('doubles rate (all cuts)')
#     Graphs['all_cut_accidental'].GetYaxis().SetTitle('event rate [%s^{-1}]'%(t))
#     Graphs['all_cut_accidental'].GetXaxis().SetTitle('photo-coverage [%]')
#     if site == 'boulby':
#         Graphs['all_cut_accidental'].GetYaxis().SetRangeUser(0,70/31.)
#     else:
#         Graphs['all_cut_accidental'].GetYaxis().SetRangeUser(0,400/31.)
#
#     for i,arr in enumerate(proc):
#         _str = 'scaled_%s_%s_%s_1_abs_%s' %(type[i],proc[i],loca[i],acc[i])
#         Graphs[_str].Draw('same')
#         Graphs[_str].SetLineColor(color[i])
#         Graphs[_str].SetLineStyle(lineS[i])
#         if loca[i] == 'PMT':
#             _descrpt  = '%s %s glass'%(proc[i],loca[i])
#             Graphs[_str].SetTitle(_descrpt)
#         elif loca[i] == 'FV':
#             _descrpt  = '%s %s in water (%s)'%(proc[i],loca[i],type[i])
#             Graphs[_str].SetTitle(_descrpt)
#         else:
#             if proc[i]=='imb':
#                 Graphs[_str].SetTitle('Prompt-positron at Perry')
#             if proc[i]=='boulby':
#                 Graphs[_str].SetTitle('Prompt-positron at Boulby')
#         Graphs[_str].SetMarkerStyle(1)
#
#     t1 = TLatex()
#     if site == 'boulby':
#         line = TLine(7, 4.,42, 4)
#         t1.DrawLatex( 10,65., 'pe > %d.0'%(cut))
#     else:
#         line = TLine(7, 40,42, 40)
#         t1.DrawLatex( 10,363., 'pe > %d.0'%(cut))
#
#     line.SetLineColor(4)
#     line.Draw('same')
#     C1.SetGridy()
#     C1.Update()
#
#     #    C1.SaveAs('gif/pc_%s_%d.gif'%(site,cut))
#     C1.SaveAs('gif/pc_%s.gif+'%(site))
#     return hist
#
#
# def sensitivityMap():
#     # Historical, not currently in use
#     a =1
#
#
# def sensitivityMapNew():
#
#     # Function to find the optimal signal to background as a function
#
#     detectorMedium      = 1
#     detectorMass        = 1000.
#     reactorPower        = 0.04
#     reactorStandoff     = 25.0
#     experiment = nuOsc.NeutrinoOscillation(detectorMedium,detectorMass,reactorPower,reactorStandoff)
#     preOsc,afterOsc = experiment.FindRate()
#     print preOsc,afterOsc
#     site = arguments["--site"]
#     if site != 'boulby':
#         site = ''
#     # Need to fix this for future running
#
#     OnOffRatio = float(arguments["--OnOff"])
#     print site,OnOffRatio
#
#     cores = int(arguments["--cores"])
#
#     if arguments["--RNRedux"]:
#         rnRedux = float(arguments["--RNRedux"])
#         if rnRedux>1:
#             print "Value of reduction of radionuclide greater than 1, setting to 0"
#             rnRedux = 0.0
#     else:
#         rnRedux = 0.0
#
#     if t == 'sec':
#         timeAdjustment = 24*3600.
#
#     if t == 'day':
#         timeAdjustment = 1.0
#
#     if t == 'month':
#         timeAdjustment = 1./31.
#
#     if t == 'year':
#         timeAdjustment = 1./365.
#
#     maxTime = 14400.*timeAdjustment
#
#     proc        = []
#     loca        = []
#     type        = []
#     color       = []
#     lineS       = []
#     acc         = []
#     scale       = []
#
#     parameters  = loadAnalysisParameters(t)
#     rates       = parameters[11]
#     FVkTonRatio = (pow(fidRadius,2)*fidHeight)/(pow(detectorRadius,2)*detectorHeight)
#     boulbyRate,imbRate = rates["boulby_S"]*FVkTonRatio,rates["imb_S"]*FVkTonRatio
#     print 'rates:',imbRate,boulbyRate, ' per ', t
#
#     #fast neutrons
#     proc        += ['QGSP_BERT_EMV','QGSP_BERT_EMX','QGSP_BERT','QGSP_BIC',\
#     'QBBC','QBBC_EMZ','FTFP_BERT','QGSP_FTFP_BERT']
#     _t          = 'FN%s' % (site)
#     loca        += [_t,_t, _t,_t,_t,_t,_t,_t]
#     type        += ['si','si','si','si','si','si','si','si']
#     acc         += ['corr','corr','corr','corr','corr','corr','corr', 'corr']
#     color       += [kO+6,kO+6,kO+6,kO+6,kO+6,kO+6,kO+6,kO+6]
#     lineS       += [2,2,2,2,2,2,2,2]
#     scale       += [1./8.,1./8.,1./8.,1./8.,1./8.,1./8.,1./8.,1./8.]
#     #radionuclides
#     proc        += ['9003','11003']
#     _t          =  'RN%s' % (site)
#     loca        += [_t,_t]
#     type        += ['ei','ei']
#     acc         += ['di','di']
#     color       += [kG+3,kG+3]
#     lineS       += [1,1]
#     scale       += [1.,1.]
#     #ibds
#     if site == 'boulby':
#         proc    += ['boulby','boulby','neutron']
#     else:
#         proc    += ['imb','imb','neutron']
#     loca        += ['S','S','N%s'%(site)]
#     type        += ['ei','ei','ei']
#     acc         += ['di','corr','corr']
#     color       += [kA+0,kA-0,kA-0]
#     lineS       += [1,2,2]
#     scale       += [-1.0,0.0,0.0]
#
#     c1 = TCanvas('c1','c1',1618,1000)
#     c1.SetRightMargin(0.53)
#     c1.Divide(2,2)
#     additionalString,additionalCommands,additionalMacStr,additionalMacOpt = testEnabledCondition(arguments)
#
#     if additionalString == "":
#         additionalString = "_default"
#
#     if site == 'boulby':
#         location = 'Boulby '
#     else:
#         location = 'Fairport '
#
#     hist = TH2D('hist','3#sigma discovery phase space -  %s '%(location),31,9.5,40.5,30,9.5,39.5)
#     hist.SetXTitle('photocoverage [%]')
#     hist.SetYTitle('photoelectron threhsold cut [p.e.]')
#     hist.SetZTitle('off-time [%s] for 3 #sigma discovery'%(t))
#     hist.GetZaxis().SetTitleOffset(-.55);
#     hist.GetZaxis().SetTitleColor(1);
#     hist.GetZaxis().CenterTitle();
#
#     gStyle.SetOptStat(0)
#     gStyle.SetPalette(55)
#
#     fileIN = 'processed_watchman_March%s.root' %(additionalString)
#     fileIN = 'ntuple_root_files%s/processed_watchman.root' %(additionalString)
#
#     print 'Will read in', fileIN
#     print "%s %s %10s %10s %10s %10s %10s %10s %10s %10s %10s" %('PC','cut','S','eff','BAC','BRN','BFN','BRS','B_tot','SSBB','T3SIGMA')
#     _graph = Graph()
#     _graph.SetName('_graph')
#     _graph.SetTitle('off-time [%s] required till a 3#sigma discovery'%(t))
#
#     _graphS = Graph()
#     _graphS.SetName('_graphS')
#     _graphS.SetTitle('IBDs from %d-core at full power'%(cores))
#
#     _graphB = Graph()
#     _graphB.SetName('_graphB')
#     _graphB.SetTitle('total background rate after optimization')
#
#     _graphBFN = Graph()
#     _graphBFN.SetLineColor(kO+6)
#     _graphBFN.SetName('_graphBFN')
#     _graphBFN.SetTitle('total fast neutron component')
#
#     _graphBRN = Graph()
#     _graphBRN.SetLineColor(kG+3)
#     _graphBRN.SetName('_graphBFN')
#     _graphBRN.SetTitle('total radionuclide component')
#
#     _graphBAC = Graph()
#     _graphBAC.SetLineColor(2)
#     _graphBAC.SetName('_graphBFN')
#     _graphBAC.SetTitle('total accidental component')
#
#     _graphBRS = Graph()
#     _graphBRS.SetLineColor(4)
#     _graphBRS.SetName('_graphBRS')
#     _graphBRS.SetTitle('other reactor antineutrinos')
#
#     _graphSIZE = Graph()
#     _graphSIZE.SetLineColor(4)
#     _graphSIZE.SetName('_graphBRS')
#     _graphSIZE.SetTitle('3-month discovery of 40 MWth at 25 km')
#
#     graphDict = {}
#     for pcVal in range(10,41,5):
#         graphDict["%dpct"%(pcVal)] = Graph()
#         graphDict["%dpct"%(pcVal)].SetLineColor(4)
#         graphDict["%dpct"%(pcVal)].SetName('_graph%dpct'%(pcVal))
#         graphDict["%dpct"%(pcVal)].SetTitle('40 MWth standoff Mass relation - signal and background scaled from %dpct values'%(pcVal))
#         graphDict["%dpct_sys10"%(pcVal)] = Graph()
#         graphDict["%dpct_sys10"%(pcVal)].SetLineColor(4)
#         graphDict["%dpct_sys10"%(pcVal)].SetName('_graph%dpct_sys10'%(pcVal))
#         graphDict["%dpct_sys10"%(pcVal)].SetTitle('%dpct values - sys 10 percent'%(pcVal))
#         graphDict["%dpct_sys20"%(pcVal)] = Graph()
#         graphDict["%dpct_sys20"%(pcVal)].SetLineColor(4)
#         graphDict["%dpct_sys20"%(pcVal)].SetName('_graph%dpct_sys20'%(pcVal))
#         graphDict["%dpct_sys20"%(pcVal)].SetTitle('%dpct values - sys 20 percent'%(pcVal))
#
#
#     _graph.SetLineWidth(2)
#     _graphS.SetLineWidth(2)
#     _graphS.SetLineColor(4)
#     _graphB.SetLineWidth(2)
#     _graphBAC.SetLineWidth(2)
#     _graphBFN.SetLineWidth(2)
#     _graphBRN.SetLineWidth(2)
#     _graphBRS.SetLineWidth(2)
#     _graphBAC.SetLineStyle(2)
#     _graphBFN.SetLineStyle(2)
#     _graphBRN.SetLineStyle(2)
#     _graphBRS.SetLineStyle(2)
#
#     minY    =   0
#
#     gCnt = 0
#     for PC in drange(10.,41.,5.0):
#         #        print 'processing photocoverage ',PC,'% :',
#         S_tmp           = 0.
#         B_tmp           = 0.0
#         BFN_tmp           = 0.0
#         BRN_tmp           = 0.0
#         BAC_tmp           = 0.0
#         BRS_tmp           = 0.0
#         SoverS_B        = 0.
#         SoverS_B_sys20  = 0.
#         T3SIGMA_MIN     = 16000000
#         cut_tmp     = 0.
#
#         c1.cd(1)
#         hist.Draw('colz')
#         gPad.SetTicks()
#
#         c1.cd(2)
#         _graph.Draw('AL')
#         _graph.GetXaxis().SetTitle('photocoverage [%]')
#         _graph.GetYaxis().SetTitle('time [%s]'%(t))
#         gPad.SetTicks()
#         gPad.SetGridy()
#
#         c1.cd(3)
#         _graphS.Draw('AL')
#         y = _graphS.GetYaxis().GetXmax()
#         _graphS.GetYaxis().SetRangeUser(0.,_graphS.GetYaxis().GetXmax()*2.)
#         _graphS.GetXaxis().SetTitle('photocoverage [%]')
#         _graphS.GetYaxis().SetTitle('event rate [%s^{-1}]'%(t))
#         c1.cd(3).BuildLegend(0.42, 0.84, 0.85, 0.88)
#         gPad.SetGridy()
#         gPad.SetTicks()
#
#         c1.cd(4)
#         _graphB.Draw('AL')
#         _graphB.GetYaxis().SetRangeUser(0., _graphB.GetYaxis().GetXmax()*2.)
#         _graphB.GetXaxis().SetTitle('photocoverage [%]')
#         _graphB.GetYaxis().SetTitle('event rate [%s^{-1}]'%(t))
#         _graphBAC.Draw('same')
#         _graphBFN.Draw('same')
#         _graphBRN.Draw('same')
#         c1.cd(4).BuildLegend(0.42, 0.72, 0.85, 0.88)
#         gPad.SetGridy()
#         gPad.SetTicks()
#
#         c1.Update()
#
#         for cut in range(4,40):
#             Graphs =obtainAbsoluteEfficiency(fileIN,timeScale=t,cut=cut)
#
#             S       = 0.  # Signal
#             BFN     = 0.  # Fast Neutron backgrounds
#             BRN     = 0.  # Radionuclides backgrounds
#             BAC     = 0.  # Accidental backgrounds
#             BRS     = 0.  # Reactor Signal backgrounds
#             _str    = 'all_cut_accidental'
#             BAC = Graphs[_str].Eval(PC)
#             for i,arr in enumerate(proc):
#                 _str = 'scaled_%s_%s_%s_1_abs_%s' %(type[i],proc[i],loca[i],acc[i])
#                 if scale[i] > 0.:
#                     if loca[i] == 'FN%s'%(site):
# 			#print BFN, _str
# 			#print Graphs
#                         BFN += Graphs[_str].Eval(PC)*scale[i]
#                         # print 'BFN',BFN,Graphs[_str].Eval(PC),scale[i]
#                     elif loca[i] == 'RN%s'%(site):
#                         BRN += Graphs[_str].Eval(PC)*scale[i]*(1.0-rnRedux)
#                         # print 'BRN',BRN,Graphs[_str].Eval(PC),scale[i]
#                 elif scale[i] < 0:
#
#                     S =Graphs[_str].Eval(PC)
#
#                     #Include 29% other reactor
#                     if site == 'boulby' and cores ==1 :
#                         BRS    = S*1.29 # 29.% other reactor
#                         B      = BAC+BRN+BFN+BRS
#                         SSBB   = S/sqrt(B + (S+B)/OnOffRatio)
#                         T3SIGMA = 9.*(B +(S+B)/OnOffRatio)/S/S
#                         nuRate = boulbyRate
#                     elif site == 'boulby' and cores ==2 :
#                         S       *= 2.                 #2-cores
#                         BRS     = S * 0.29
#                         B       = BAC+BRN+BFN+BRS
#                         SSBB    = S/sqrt(S+B)
#                         T3SIGMA = 9.*(S+B)/S/S #1.29% other reactor
#                         nuRate = boulbyRate
#
#                     else:
#                         BRS     = S * 0.05
#                         B       = BAC+BRN+BFN+BRS
#                         SSBB   = S/sqrt(B + (S+B)/OnOffRatio)
#                         T3SIGMA = 9.*(B + (S+B)/OnOffRatio)/S/S#5% other reactor
#                         nuRate = imbRate
#                     # For debugging purposes
#                     # print "%3d %3d %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f" %(PC,cut,S,S/nuRate,BAC,BRN,BFN,BRS,B,SSBB,T3SIGMA)
#
#                     if  T3SIGMA < maxTime:
#                         print "%3d %3d %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f" %(PC,cut,S,S/nuRate,BAC,BRN,BFN,BRS,B,SSBB,T3SIGMA)
#                         # print "%3d %3d %10.4e %10.4f %10.4e %10.4e %10.4e %10.4e %10.4e" %(PC,cut,S,S/nuRate,BAC,BRN,BFN,BRS,B)
#                         hist.Fill(PC,cut,T3SIGMA)
#                     if  T3SIGMA >= maxTime:
#                         hist.Fill(PC,cut,maxTime)
#                     if T3SIGMA <   T3SIGMA_MIN:
#                         S_tmp         = S
#                         B_tmp         = B
#                         BFN_tmp       = BFN
#                         BRN_tmp       = BRN
#                         BAC_tmp       = BAC
#                         BRS_tmp       = BRS
#                         T3SIGMA_MIN   = T3SIGMA
#                         _signal_Eff   = S/nuRate
#
#         _graph.SetPoint(gCnt,PC,T3SIGMA_MIN)
#         _graphS.SetPoint(gCnt,PC,S_tmp)
#         _graphB.SetPoint(gCnt,PC,B_tmp)
#         _graphBFN.SetPoint(gCnt,PC,BFN_tmp)
#         _graphBRN.SetPoint(gCnt,PC,BRN_tmp)
#         _graphBAC.SetPoint(gCnt,PC,BAC_tmp)
#         _graphBRS.SetPoint(gCnt,PC,BRS_tmp)
#
#         # _sOff,_T3sOff = findEquivalentStandoff(S_tmp,BAC_tmp+BFN_tmp,BRN_tmp,_signal_Eff)
#
#         if arguments["--40MWth"]:
#             _size,_T3 = findScaledVolume(S_tmp,BAC_tmp+BFN_tmp,BRN_tmp,1.0)
#             graphCnt = 0
#             sigmaLvl = float(arguments["--40MWthSig"])
#             if 1==1:#PC == 40:
#                 for _detectSize in drange(1.,50.,1.):
#                     _sOff,_sRate,_sRatePostEff = findEquivalentStandoff(S_tmp,\
#                     BAC_tmp+BFN_tmp,BRN_tmp,_signal_Eff,0.,6,_detectSize,sigmaLvl)
#                     graphDict["%dpct"%(PC)].SetPoint(graphCnt,_sOff,_detectSize)
#                     #Do the 10% systematics
#                     _sOff,_sRate,_sRatePostEff = findEquivalentStandoff(S_tmp,\
#                     BAC_tmp+BFN_tmp,BRN_tmp,_signal_Eff,0.05,6,_detectSize,sigmaLvl)
#                     graphDict["%dpct_sys10"%(PC)].SetPoint(graphCnt,_sOff,_detectSize)
#                     #Do the 20% systematics
#                     _sOff,_sRate,_sRatePostEff = findEquivalentStandoff(S_tmp,\
#                     BAC_tmp+BFN_tmp,BRN_tmp,_signal_Eff,0.1,6,_detectSize,sigmaLvl)
#                     graphDict["%dpct_sys20"%(PC)].SetPoint(graphCnt,_sOff,_detectSize)
#                     graphCnt+=1
#
#
#             _graphSIZE.SetPoint(gCnt,PC,_size)
#         # print PC,_size,_T3,_sOff,_T3sOff
#
#         gCnt+=1
#
#     c1.cd(1)
#     hist.Draw('colz')
#     gPad.SetTicks()
#
#     c1.cd(2)
#     _graph.Draw('AL')
#     _graph.GetXaxis().SetTitle('photocoverage [%]')
#     _graph.GetYaxis().SetTitle('time [%s]'%(t))
#     gPad.SetTicks()
#     gPad.SetGridy()
#
#
#     c1.cd(3)
#     _graphS.Draw('AL')
#     y = _graphS.GetYaxis().GetXmax()
#     _graphS.GetYaxis().SetRangeUser(0.,_graphS.GetYaxis().GetXmax()*2.)
#     _graphS.GetXaxis().SetTitle('photocoverage [%]')
#     _graphS.GetYaxis().SetTitle('event rate [%s^{-1}]'%(t))
#     c1.cd(3).BuildLegend(0.42, 0.84, 0.85, 0.88)
#     gPad.SetGridy()
#     gPad.SetTicks()
#
#     c1.cd(4)
#     _graphS.Draw('AL')
#     if site == 'boulby':
#         _graphS.GetYaxis().SetRangeUser(0., _graphS.GetYaxis().GetXmax()*3.)
#     else:
#         _graphS.GetYaxis().SetRangeUser(0., _graphS.GetYaxis().GetXmax()*2.)
#     _graphS.GetXaxis().SetTitle('photocoverage [%]')
#     _graphS.GetYaxis().SetTitle('event rate [%s^{-1}]'%(t))
#     _graphB.Draw('same')
#     _graphBAC.Draw('same')
#     _graphBFN.Draw('same')
#     _graphBRN.Draw('same')
#     _graphBRS.Draw('same')
#     c1.cd(4).BuildLegend(0.42, 0.72, 0.85, 0.88)
#     gPad.SetGridy()
#     gPad.SetTicks()
#
#     c3 = TCanvas('c3','c3',1618,1000)
#     c3.Divide(1,1)
#     c3.cd(1)
#     _graphS.Draw('AL')
#     if site == 'boulby':
#         _graphS.GetYaxis().SetRangeUser(0., _graphS.GetYaxis().GetXmax()*3.)
#     else:
#         _graphS.GetYaxis().SetRangeUser(0., _graphS.GetYaxis().GetXmax()*2.)
#     _graphS.GetXaxis().SetTitle('photocoverage [%]')
#     _graphS.GetYaxis().SetTitle('event rate [%s^{-1}]'%(t))
#     _graphB.Draw('same')
#     _graphBAC.Draw('same')
#     _graphBFN.Draw('same')
#     _graphBRN.Draw('same')
#     _graphBRS.Draw('same')
#     c3.cd(1).BuildLegend(0.42, 0.72, 0.85, 0.88)
#     gPad.SetGridy()
#     gPad.SetTicks()
#     c3.Update()
#
#
#     c2 = TCanvas('c2','c2',1618,1000)
#     # c2.SetRightMargin(0.53)
#     c2.Divide(1,1)
#     c2.cd(1)
#     _graph.Draw('AL')
#     _graph.GetXaxis().SetTitle('photocoverage [%]')
#     _graph.GetYaxis().SetTitle('time [%s]'%(t))
#     gPad.SetTicks()
#     gPad.SetGridy()
#
#     c4 = TCanvas('c4','c4',1618,1000)
#     c4.Divide(1,1)
#     c4.cd(1)
#     _graphSIZE.Draw('AL')
#     _graphSIZE.GetXaxis().SetTitle('photocoverage [%]')
#     _graphSIZE.GetYaxis().SetTitle('Detector FV size [kton]')
#     c4.Update()
#     if arguments["--40MWth"]:
#         _pc = ['10','15','20','25','30','35','40']
#         _pc = ['40']
#         _Canvas = {}
#         for pc in _pc:
#             _Canvas["c%s"%(pc)] = TCanvas("c%s"%(pc),"c%s"%(pc),1618,1000)
#             _Canvas["c%s"%(pc)].Divide(1,1)
#             graphDict["%spct"%(pc)].Draw('ALP')
#             graphDict["%spct_sys10"%(pc)].Draw('LPsame')
#             graphDict["%spct_sys20"%(pc)].Draw('LPsame')
#             _Canvas["c%s"%(pc)].Update()
#             _Canvas["c%s"%(pc)].SaveAs('results/Canvas%s.gif'%(pc))
#             _Canvas["c%s"%(pc)].SaveAs('results/Canvas%s.C'%(pc))
#
#
#     analysisString =""
#     analysisString += arguments["--site"]
#     analysisString += arguments["--cores"]
#     analysisString += arguments["--RNRedux"]
#     analysisString += arguments["--timeScale"]
#     analysisString += arguments["--OnOff"]
#
#     hist.SaveAs('results/SensitivityMap%s%s.C'%(additionalString,analysisString))
#     _graph.SaveAs('results/SensitivityGraph%s%s.C'%(additionalString,analysisString))
#     _graphS.SaveAs('results/SensitivityGraphS%s%s.C'%(additionalString,analysisString))
#     _graphB.SaveAs('results/SensitivityGraphB%s%s.C'%(additionalString,analysisString))
#
#     c1.SaveAs('results/SensitivityCanvas%s%s.C'%(additionalString,analysisString))
#     c1.SaveAs('results/SensitivityCanvas%s%s.gif'%(additionalString,analysisString))
#
#     c2.SaveAs('results/SensitivityTime%s%s.C'%(additionalString,analysisString))
#     c2.SaveAs('results/SensitivityTime%s%s.gif'%(additionalString,analysisString))
#
#     c3.SaveAs('results/SensitivityRate%s%s.C'%(additionalString,analysisString))
#     c3.SaveAs('results/SensitivityRate%s%s.gif'%(additionalString,analysisString))
#
#     c4.SaveAs('results/SensitivitySize%s%s.C'%(additionalString,analysisString))
#     c4.SaveAs('results/SensitivitySize%s%s.gif'%(additionalString,analysisString))
#
#
# def runSensitivity():
#     hBoulby = TH2D('hBoulby','hBoulby',50,0.5,50.5,50,0.5,50.5)
#     print 'Boulby:'
#     for cut in range(4,40):
#         Graphs = obtainAbsoluteEfficiency('processed_data_watchman.root',\
#         timeScale=t,cut=cut)
#         hBoulby =  siteCanvas(Graphs,'boulby',cut,hBoulby)
#     hBoulby.SaveAs('hBoulby.C')
#
#
#     hFairport = TH2D('hFairport','hFairport',50,0.5,50.5,50,0.5,50.5)
#     print 'Fairport:'
#     for cut in range(4,40):
#         Graphs = obtainAbsoluteEfficiency('processed_data_watchman.root',\
#         timeScale=t,cut=cut)
#         hFairport =  siteCanvas(Graphs,'',cut,hFairport)
#     hFairport.SaveAs('hFairport.C')
#
#     return 0
#
# def findScaledVolume(S,B_Surface,B_Volume,OnOffRatio = 1.0):
#     B_Volume += S*0.02 # Reactor background for site of interest 15TNU/747TNU
#     S *= 40./1500. # Transform Boubly rate to 40MWth rate
#
#     for size in drange(1.,500.,0.5):
#         _rad    = pow(size*1000./pi/2.,1./3.)
#         _surf   = pow(_rad,2)/pow(fidRadius,2)
#         _volume = pow(_rad,3)/(pow(fidRadius,2)*fidHeight)
#         _S      = S * _volume
#         _B      = B_Surface * _surf + B_Volume * _volume
#         SSBB   = _S/sqrt(_B + (_S + _B)/OnOffRatio)
#
#         T3SIGMA = 9.*(_B +(_S+_B)/OnOffRatio)/_S/_S
#         # except:
#         #     print 'Could not'
#         # print "%3d %10.2f %10.2f %10.2f %10.2f"%(size, SSBB, T3SIGMA,_S,_B)
#         if T3SIGMA < 3.00001:
#             break
#     return size,T3SIGMA
#
#
#
# def findEquivalentStandoff(S,B_Surface,B_Volume,_signal_Eff,_sigma_b = 0.1,\
# t = 3.0,detectorMass=1.,N_sigma=3):
#     B_Volume += S*0.02 # Reactor background for site of interest 15TNU/747TNU
#
#     # First, find rate on observed signal needed
#     _rad    = pow(detectorMass*1000./pi/2.,1./3.)
#     _surf   = pow(_rad,2)/pow(fidRadius,2)
#     _volume = pow(_rad,3)/(pow(fidRadius,2)*fidHeight)
#     _b = B_Surface*_surf + B_Volume*_volume
#     _sigma_b*=_b
#
#     # Evaluate the required signal to obtain N_sigma over t
#     C_square = N_sigma*N_sigma
#     _s = C_square/(2.*t)*(1+sqrt(1. + 4.*_b*t/C_square+\
#     4.*_sigma_b*_sigma_b*t*t/C_square))
#
#     if _s > _sigma_b:
#         _t = C_square*(_s+_b)/(_s*_s - _sigma_b*_sigma_b*C_square)
#     else:
#         print 'Background systematics greater than signal'
#         return -1,-1,-1
#     print 'd,r,s,v: ',detectorMass,_rad,_surf,_volume,_b,_sigma_b,_s,t,_t,
#
#     # Define the medium of the detector (1:water) and power of test reactor
#     detectorMedium,reactorPower      = 1,0.04 #0.04
#     # Cycle through the standoff and stop when signal is observed
#     experiment = nuOsc.NeutrinoOscillation(detectorMedium,\
#     1.*1000.,reactorPower,25)
#     preOsc,afterOsc = experiment.FindRate()
#     print preOsc,afterOsc,
#     for _standoff in drange(0.05,50,0.05):
#             experiment = nuOsc.NeutrinoOscillation(detectorMedium,\
#             detectorMass*1000.,reactorPower,_standoff)
#             preOsc,afterOsc = experiment.FindRate()
#             preOsc,afterOsc = preOsc*365./12.,afterOsc*365./12.
#             if afterOsc*_signal_Eff < _s:
#                 break
#
#     print 'result:',_standoff,afterOsc*_signal_Eff
#     return _standoff, afterOsc,afterOsc*_signal_Eff
