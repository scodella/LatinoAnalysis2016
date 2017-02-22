import optparse
import numpy
import ROOT
import os.path
import math

from LatinoAnalysis.Gardener.gardening import TreeCloner

#
#  _ _|      |         /     _ _|                    ___|                |            ____|             |                           
#    |    _` |        /        |    __|   _ \      \___ \    __|   _` |  |   _ \      |     _` |   __|  __|   _ \    __|  __| 
#    |   (   |       /         |  \__ \  (   |           |  (     (   |  |   __/      __|  (   |  (     |    (   |  |   \__ \ 
#  ___| \__,_|     _/        ___| ____/ \___/      _____/  \___| \__,_| _| \___|     _|   \__,_| \___| \__| \___/  _|   ____/ 
#                                                                                                                             
#

class IdIsoSFStopFiller(TreeCloner):

    def __init__(self):
        pass

    def __del__(self):
        pass

    def help(self):
        return '''Add a new lepton scale factor weight based on id/isolation scale factors data/MC.'''

    def addOptions(self,parser):
        description = self.help()
        group = optparse.OptionGroup(parser,self.label, description)
        group.add_option('-c', '--cmssw', dest='cmssw', help='cmssw version (naming convention may change)', default='ICHEP2016', type='string')
        group.add_option('-f', '--readfastsim', dest='readfastsim', help='read FastSim SFs', default=0, type='int')
        group.add_option( '--idLepKind' , dest='idLepKind', help='kind of lepton id', default=None) 
        group.add_option( '--BCDEFtoGHRatio', dest='BCDEFtoGHRatio', help='Ratio of BCDEF to GH data used', type='float'  ,    default=0.549763) # 19.72/35.87
        parser.add_option_group(group)
        return group



    def checkOptions(self,opts):

        cmssw_base = os.getenv('CMSSW_BASE')

        if opts.cmssw == "ICHEP2016" :
            
            # https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF#Data_leading_order_FullSim_MC_co
            self.fileMuonId = self._openRootFile(cmssw_base+'/src/LatinoAnalysis/Gardener/python/data/idiso/ICHEP2016/Stop/TnP_MuonID_NUM_MediumID_DENOM_generalTracks_VAR_map_pt_eta.root') 
            self.fileMuonIso = self._openRootFile(cmssw_base+'/src/LatinoAnalysis/Gardener/python/data/idiso/ICHEP2016/Stop/TnP_MuonID_NUM_MultiIsoVT_DENOM_MediumID_VAR_map_pt_eta.root') 
            self.fileMuonIP2D = self._openRootFile(cmssw_base+'/src/LatinoAnalysis/Gardener/python/data/idiso/ICHEP2016/Stop/TnP_MuonID_NUM_TightIP2D_DENOM_MediumID_VAR_map_pt_eta.root')
            self.fileMuonIP3D = self._openRootFile(cmssw_base+'/src/LatinoAnalysis/Gardener/python/data/idiso/ICHEP2016/Stop/TnP_MuonID_NUM_TightIP3D_DENOM_MediumID_VAR_map_pt_eta.root') 

            self.MuId   = self._getRootObj(self.fileMuonId,   'pt_abseta_PLOT_pair_probeMultiplicity_bin0')
            self.MuIso  = self._getRootObj(self.fileMuonIso,  'pt_abseta_PLOT_pair_probeMultiplicity_bin0_&_Medium2016_pass')  
            self.MuIP2D = self._getRootObj(self.fileMuonIP2D, 'pt_abseta_PLOT_pair_probeMultiplicity_bin0_&_Medium2016_pass') 
            self.MuIP3D = self._getRootObj(self.fileMuonIP3D, 'pt_abseta_PLOT_pair_probeMultiplicity_bin0_&_Medium2016_pass')

            # See also https://twiki.cern.ch/twiki/bin/view/CMS/MuonWorkInProgressAndPagResults
            self.fileMuonReco = self._openRootFile(cmssw_base+'/src/LatinoAnalysis/Gardener/python/data/idiso/ICHEP2016/Stop/ratios_MuonTrk.root')
            self.MuReco = self._getRootObj(self.fileMuonReco, 'ratio_eta') 
            
            # https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF#FullSim_FastSim_TTBar_MC_compari
            self.fileFastSimMuonId = self._openRootFile(cmssw_base+'/src/LatinoAnalysis/Gardener/python/data/idiso/ICHEP2016/Stop/sf_mu_medium.root') 
            self.fileFastSimMuonIso = self._openRootFile(cmssw_base+'/src/LatinoAnalysis/Gardener/python/data/idiso/ICHEP2016/Stop/sf_mu_mediumID_multiVT.root') 
            self.fileFastSimMuonIP2D = self._openRootFile(cmssw_base+'/src/LatinoAnalysis/Gardener/python/data/idiso/ICHEP2016/Stop/sf_mu_tightIP2D.root')
            self.fileFastSimMuonIP3D = self._openRootFile(cmssw_base+'/src/LatinoAnalysis/Gardener/python/data/idiso/ICHEP2016/Stop/sf_mu_tightIP3D.root')
            
            self.FastSimMuId   = self._getRootObj(self.fileFastSimMuonId,   'histo2D')
            self.FastSimMuIso  = self._getRootObj(self.fileFastSimMuonIso,  'histo2D') 
            self.FastSimMuIP2D = self._getRootObj(self.fileFastSimMuonIP2D, 'histo2D') 
            self.FastSimMuIP3D = self._getRootObj(self.fileFastSimMuonIP3D, 'histo2D') 

            # https://twiki.cern.ch/twiki/bin/viewauth/CMS/SUSLeptonSF#Data_leading_order_FullSim_M_AN1
            self.fileElectronIdIso = self._openRootFile(cmssw_base+'/src/LatinoAnalysis/Gardener/python/data/idiso/ICHEP2016/Stop/scaleFactors__ElectronIdIso.root')
            self.ElId   = self._getRootObj(self.fileElectronIdIso,  'GsfElectronToTight')
            self.ElIso  = self._getRootObj(self.fileElectronIdIso,  'CutBasedTightElectronToMultiIsoVT')

            # See also https://twiki.cern.ch/twiki/bin/view/CMS/EgammaIDRecipesRun2#Electron_efficiencies_and_scale
            self.fileElectronReco = self._openRootFile(cmssw_base+'/src/LatinoAnalysis/Gardener/python/data/idiso/ICHEP2016/Stop/egammaEffi.txt_SF2D.root')
            self.ElReco = self._getRootObj(self.fileElectronReco,  'EGamma_SF2D')
            
            # https://twiki.cern.ch/twiki/bin/viewauth/CMS/SUSLeptonSF#FullSim_FastSim_TTBar_MC_com_AN1
            self.fileFastSimElectronIdIso = self._openRootFile(cmssw_base+'/src/LatinoAnalysis/Gardener/python/data/idiso/ICHEP2016/Stop/sf_el_tightCB_MultiVT__FastSim.root')
            self.FastSimElIdIso = self._getRootObj(self.fileFastSimElectronIdIso,   'histo2D')
            
        elif opts.cmssw == "Full2016" : 
            
            if opts.idLepKind == "Ghent" : 

                # https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF#Scale_Factors_for_Moriond2017
                self.fileMuonId = self._openRootFile(cmssw_base+'/src/LatinoAnalysis/Gardener/python/data/idiso/Full2016/Stop/TnP_NUM_MediumID_DENOM_generalTracks_VAR_map_pt_eta.root') 
                self.fileMuonIso = self._openRootFile(cmssw_base+'/src/LatinoAnalysis/Gardener/python/data/idiso/Full2016/Stop/TnP_NUM_RelIsoVTight_DENOM_MediumID_VAR_map_pt_eta.root') 
                self.fileMuonIP2D = self._openRootFile(cmssw_base+'/src/LatinoAnalysis/Gardener/python/data/idiso/Full2016/Stop/TnP_NUM_TightIP2D_DENOM_MediumID_VAR_map_pt_eta.root')
                self.fileMuonIP3D = self._openRootFile(cmssw_base+'/src/LatinoAnalysis/Gardener/python/data/idiso/Full2016/Stop/TnP_NUM_TightIP3D_DENOM_MediumID_VAR_map_pt_eta.root') 
                
                self.MuId   = self._getRootObj(self.fileMuonId,   'SF')
                self.MuIso  = self._getRootObj(self.fileMuonIso,  'pt_abseta_ratio')  
                self.MuIP2D = self._getRootObj(self.fileMuonIP2D, 'SF') 
                self.MuIP3D = self._getRootObj(self.fileMuonIP3D, 'SF')
                
                # See also https://twiki.cern.ch/twiki/bin/view/CMS/MuonWorkInProgressAndPagResults
                self.fileMuonReco = self._openRootFile(cmssw_base+'/src/LatinoAnalysis/Gardener/python/data/idiso/Full2016/Stop/Tracking_EfficienciesAndSF_BCDEFGH.root')
                self.MuReco = self._getRootObj(self.fileMuonReco, 'ratio_eff_eta3_dr030e030_corr') 
            
                # https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF#FullSim_FastSim_TTBar_MC_compari  --> TOBEUPDATED
                self.fileFastSimMuonId = self._openRootFile(cmssw_base+'/src/LatinoAnalysis/Gardener/python/data/idiso/ICHEP2016/Stop/sf_mu_medium.root') 
                self.fileFastSimMuonIso = self._openRootFile(cmssw_base+'/src/LatinoAnalysis/Gardener/python/data/idiso/ICHEP2016/Stop/sf_mu_mediumID_multiVT.root') 
                self.fileFastSimMuonIP2D = self._openRootFile(cmssw_base+'/src/LatinoAnalysis/Gardener/python/data/idiso/ICHEP2016/Stop/sf_mu_tightIP2D.root')
                self.fileFastSimMuonIP3D = self._openRootFile(cmssw_base+'/src/LatinoAnalysis/Gardener/python/data/idiso/ICHEP2016/Stop/sf_mu_tightIP3D.root')
                
                self.FastSimMuId   = self._getRootObj(self.fileFastSimMuonId,   'histo2D')
                self.FastSimMuIso  = self._getRootObj(self.fileFastSimMuonIso,  'histo2D') 
                self.FastSimMuIP2D = self._getRootObj(self.fileFastSimMuonIP2D, 'histo2D') 
                self.FastSimMuIP3D = self._getRootObj(self.fileFastSimMuonIP3D, 'histo2D') 
                
                # https://twiki.cern.ch/twiki/bin/viewauth/CMS/SUSLeptonSF#Data_leading_order_FullSim_M_AN1
                self.fileElectronIdIso = self._openRootFile(cmssw_base+'/src/LatinoAnalysis/Gardener/python/data/idiso/Full2016/Stop/scaleFactors_electrons.root')
                self.ElId   = self._getRootObj(self.fileElectronIdIso,  'GsfElectronToCutBasedStopsDilepton')
                self.ElIso  = self._getRootObj(self.fileElectronIdIso,  'CutBasedStopsDileptonToRelIso012')

                # See also https://twiki.cern.ch/twiki/bin/view/CMS/EgammaIDRecipesRun2#Electron_efficiencies_and_scale
                self.fileElectronReco = self._openRootFile(cmssw_base+'/src/LatinoAnalysis/Gardener/python/data/idiso/Full2016/Stop/egammaEffi.txt_EGM2D.root')
                self.ElReco = self._getRootObj(self.fileElectronReco,  'EGamma_SF2D')
            
                # https://twiki.cern.ch/twiki/bin/viewauth/CMS/SUSLeptonSF#FullSim_FastSim_TTBar_MC_com_AN1  --> TOBEUPDATED
                self.fileFastSimElectronIdIso = self._openRootFile(cmssw_base+'/src/LatinoAnalysis/Gardener/python/data/idiso/ICHEP2016/Stop/sf_el_tightCB_MultiVT__FastSim.root')
                self.FastSimElIdIso = self._getRootObj(self.fileFastSimElectronIdIso,   'histo2D')
            
            elif opts.idLepKind == "POG" : 
            
                # 
                self.fileMuonId = self._openRootFile(cmssw_base+'/src/LatinoAnalysis/Gardener/python/data/idiso/Full2016/Stop/EfficienciesAndSF_GH_ID.root') 
                self.fileMuonIso = self._openRootFile(cmssw_base+'/src/LatinoAnalysis/Gardener/python/data/idiso/Full2016/Stop/EfficienciesAndSF_GH_ISO.root') 
                self.fileMuonIdBF = self._openRootFile(cmssw_base+'/src/LatinoAnalysis/Gardener/python/data/idiso/Full2016/Stop/EfficienciesAndSF_BCDEF_ID.root') 
                self.fileMuonIsoBF = self._openRootFile(cmssw_base+'/src/LatinoAnalysis/Gardener/python/data/idiso/Full2016/Stop/EfficienciesAndSF_BCDEF_ISO.root') 
                self.fileMuonIP2D = self._openRootFile(cmssw_base+'/src/LatinoAnalysis/Gardener/python/data/idiso/Full2016/Stop/TnP_NUM_TightIP2D_DENOM_MediumID_VAR_map_pt_eta.root')
                self.fileMuonIP3D = self._openRootFile(cmssw_base+'/src/LatinoAnalysis/Gardener/python/data/idiso/Full2016/Stop/TnP_NUM_TightIP3D_DENOM_MediumID_VAR_map_pt_eta.root') 
                
                self.MuId    = self._getRootObj(self.fileMuonId,    'MC_NUM_MediumID2016_DEN_genTracks_PAR_pt_eta/pt_abseta_ratio')
                self.MuIso   = self._getRootObj(self.fileMuonIso,   'TightISO_MediumID_pt_eta/pt_abseta_ratio')
                self.MuIdBF  = self._getRootObj(self.fileMuonIdBF,  'MC_NUM_MediumID2016_DEN_genTracks_PAR_pt_eta/pt_abseta_ratio')
                self.MuIsoBF = self._getRootObj(self.fileMuonIsoBF, 'TightISO_MediumID_pt_eta/pt_abseta_ratio')  
                self.MuIP2D  = self._getRootObj(self.fileMuonIP2D,  'SF') 
                self.MuIP3D  = self._getRootObj(self.fileMuonIP3D,  'SF')

                # See also https://twiki.cern.ch/twiki/bin/view/CMS/MuonWorkInProgressAndPagResults
                self.fileMuonReco = self._openRootFile(cmssw_base+'/src/LatinoAnalysis/Gardener/python/data/idiso/Full2016/Stop/Tracking_EfficienciesAndSF_BCDEFGH.root')
                self.MuReco = self._getRootObj(self.fileMuonReco, 'ratio_eff_eta3_dr030e030_corr') 
            
                # https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF#FullSim_FastSim_TTBar_MC_compari  --> TOBEUPDATED
                self.fileFastSimMuonId = self._openRootFile(cmssw_base+'/src/LatinoAnalysis/Gardener/python/data/idiso/ICHEP2016/Stop/sf_mu_medium.root') 
                self.fileFastSimMuonIso = self._openRootFile(cmssw_base+'/src/LatinoAnalysis/Gardener/python/data/idiso/ICHEP2016/Stop/sf_mu_mediumID_multiVT.root') 
                self.fileFastSimMuonIP2D = self._openRootFile(cmssw_base+'/src/LatinoAnalysis/Gardener/python/data/idiso/ICHEP2016/Stop/sf_mu_tightIP2D.root')
                self.fileFastSimMuonIP3D = self._openRootFile(cmssw_base+'/src/LatinoAnalysis/Gardener/python/data/idiso/ICHEP2016/Stop/sf_mu_tightIP3D.root')
            
                self.FastSimMuId   = self._getRootObj(self.fileFastSimMuonId,   'histo2D')
                self.FastSimMuIso  = self._getRootObj(self.fileFastSimMuonIso,  'histo2D') 
                self.FastSimMuIP2D = self._getRootObj(self.fileFastSimMuonIP2D, 'histo2D') 
                self.FastSimMuIP3D = self._getRootObj(self.fileFastSimMuonIP3D, 'histo2D') 

                # https://twiki.cern.ch/twiki/bin/viewauth/CMS/SUSLeptonSF#Data_leading_order_FullSim_M_AN1
                self.fileElectronIdIso = self._openRootFile(cmssw_base+'/src/LatinoAnalysis/Gardener/python/data/idiso/Full2016/Stop/egammaEffi.txt_EGM2D_idiso.root')
                self.ElIdIso   = self._getRootObj(self.fileElectronIdIso,  'EGamma_SF2D')

                # See also https://twiki.cern.ch/twiki/bin/view/CMS/EgammaIDRecipesRun2#Electron_efficiencies_and_scale
                self.fileElectronReco = self._openRootFile(cmssw_base+'/src/LatinoAnalysis/Gardener/python/data/idiso/Full2016/Stop/egammaEffi.txt_EGM2D.root')
                self.ElReco = self._getRootObj(self.fileElectronReco,  'EGamma_SF2D')
            
                # https://twiki.cern.ch/twiki/bin/viewauth/CMS/SUSLeptonSF#FullSim_FastSim_TTBar_MC_com_AN1  --> TOBEUPDATED
                self.fileFastSimElectronIdIso = self._openRootFile(cmssw_base+'/src/LatinoAnalysis/Gardener/python/data/idiso/ICHEP2016/Stop/sf_el_tightCB_MultiVT__FastSim.root')
                self.FastSimElIdIso = self._getRootObj(self.fileFastSimElectronIdIso,   'histo2D')
            
        self.cmssw          = opts.cmssw
        self.readfastsim    = opts.readfastsim
        self.idLepKind      = opts.idLepKind
        self.BCDEFtoGHRatio = opts.BCDEFtoGHRatio

        print " cmssw          = ", self.cmssw
        print " readfastsim    = ", self.readfastsim
        print " idLepKind      = ", self.idLepKind
        print " BCDEFtoGHRatio = ", self.BCDEFtoGHRatio

    # X axis: Pt, Y axis: abs(Eta)        
    def _getHistoValuePtAbsEta(self, h2, pt, eta):

        nPtBins = h2.GetNbinsX()
        pt = min(pt, h2.GetXaxis().GetBinCenter(nPtBins))
        pt = max(pt, h2.GetXaxis().GetBinLowEdge(1))

        eta = abs(eta)

        value = h2.GetBinContent(h2.FindBin(pt, eta))
        error = h2.GetBinError  (h2.FindBin(pt, eta))

        return value, error
            

    # X axis: Eta        
    def _getTGraphValueEta(self, h1, eta):

        for point in xrange(h1.GetN()):
            xx = yy = ROOT.Double(0)
            h1.GetPoint(point, xx, yy)
            ex = h1.GetErrorX(point)
            if (eta >= xx-ex and eta < xx+ex):
                value = yy
                error = h1.GetErrorY(point)
                return value, error

        #print "idisoScaleFactorsStop: TGraph value eta not found"
        return 1., 0.
            

    # X axis: Eta, Y axis: Pt        
    def _getHistoValueEtaPt(self, h2, pt, eta):

        nPtBins = h2.GetNbinsY()
        pt = min(pt, h2.GetYaxis().GetBinCenter(nPtBins))

        pt = max(pt, h2.GetYaxis().GetBinLowEdge(1))

        value = h2.GetBinContent(h2.FindBin(eta, pt))
        error = h2.GetBinError  (h2.FindBin(eta, pt))

        return value, error


    def _getWeight (self, kindLep, pt, eta, tight, ScaleFactorLevel, nvtx):

        if kindLep == 'ele' :
            if pt<10. or abs(eta)>=2.5 :
                return 1.0, 1.0, 1.0, 0.0

        if kindLep == 'mu' :
            if pt<10. or abs(eta)>=2.4 :
                return 1.0, 1.0, 1.0, 0.0
            
        # Id*Iso
        if ScaleFactorLevel == 0:
            
            if kindLep == 'ele' :

                scaleFactor = 1.
                relative_error_scaleFactor = 0.

                if self.idLepKind != "POG" :
                    IdSF,  IdSFErr  = self._getHistoValuePtAbsEta(self.ElId,  pt, eta)
                    IsoSF, IsoSFErr = self._getHistoValuePtAbsEta(self.ElIso, pt, eta)
                    scaleFactor = IdSF*IsoSF
                    relative_error_scaleFactor = math.sqrt( (IdSFErr/IdSF)*(IdSFErr/IdSF) + (IsoSFErr/IsoSF)*(IsoSFErr/IsoSF) )
                else :
                    idisoSF, idisoSFErr  = self._getHistoValueEtaPt(self.ElIdIso,  pt, eta)
                    scaleFactor = idisoSF
                    relative_error_scaleFactor = idisoSFErr/idisoSF

                # This includes systematic errors. Put on statistics because of how it's written AnalysisCMS
                scaleFactor_do = scaleFactor*(1. -  relative_error_scaleFactor)
                scaleFactor_up = scaleFactor*(1. +  relative_error_scaleFactor)
          
                # No systematic uncertainty for the time being
                return scaleFactor, scaleFactor_do, scaleFactor_up, 0.0
            
            elif kindLep == 'mu' :

                IdSF,   IdSFErr   = self._getHistoValuePtAbsEta(self.MuId,   pt, eta)
                IsoSF,  IsoSFErr  = self._getHistoValuePtAbsEta(self.MuIso,  pt, eta)
                IP2DSF, IP2DSFErr = self._getHistoValuePtAbsEta(self.MuIP2D, pt, eta)
                IP3DSF, IP3DSFErr = self._getHistoValuePtAbsEta(self.MuIP3D, pt, eta)

                if self.idLepKind == "POG" :
                    toss_a_coin = ROOT.gRandom.Rndm()
                    if toss_a_coin <=  self.BCDEFtoGHRatio :
                        IdSF,   IdSFErr   = self._getHistoValuePtAbsEta(self.MuIdBF,   pt, eta)
                        IsoSF,  IsoSFErr  = self._getHistoValuePtAbsEta(self.MuIsoBF,  pt, eta)

                scaleFactor = IdSF*IsoSF*IP2DSF*IP3DSF # ?

                # This includes systematic errors. Put on statistics because of how it's written AnalysisCMS  
                relative_error_scaleFactor = 0.03
                if self.idLepKind == "POG" : # Preliminary
                    relative_error_scaleFactor = math.sqrt( (IdSFErr/IdSF)*(IdSFErr/IdSF) + 0.01*0.01 +
                                                            (IsoSFErr/IsoSF)*(IsoSFErr/IsoSF) + 0.005*0.005 )
                scaleFactor_do = scaleFactor*(1. -  relative_error_scaleFactor)
                scaleFactor_up = scaleFactor*(1. +  relative_error_scaleFactor)

                # No systematic uncertainty for the time being
                return scaleFactor, scaleFactor_do, scaleFactor_up, 0.0
                
            # not a lepton ... like some default value: and what can it be if not a lepton? ah ah 
            # --> it happens for default values -9999.
            return 1.0, 1.0, 1.0, 0.0
        
        # Reco scale factor
        elif ScaleFactorLevel == 1:
            
            if kindLep == 'ele' : # https://twiki.cern.ch/twiki/bin/view/CMS/EgammaIDRecipesRun2#Electron_efficiencies_and_scale

                relative_error_pt = 0.
                if pt < 20 :
                    if self.cmssw == "ICHEP2016" :
                        relative_error_pt = 0.03
                    elif self.cmssw == "Full2016" :
                        relative_error_pt = 0.01
                    pt = 20
                elif pt > 80 :
                    if self.cmssw == "Full2016" :
                        relative_error_pt = 0.01

                relative_error_pu = 0.
                if self.cmssw == "ICHEP2016" :
                    relative_error_pu = 0.01

                RecoSF, RecoSFErr  = self._getHistoValueEtaPt(self.ElReco,  pt, eta)

                scaleFactor = RecoSF

                # This includes systematic errors. Put on statistics because of how it's written AnalysisCMS
                relative_error_scaleFactor = math.sqrt( (RecoSFErr/RecoSF)*(RecoSFErr/RecoSF) + relative_error_pt*relative_error_pt + relative_error_pu*relative_error_pu )
                scaleFactor_do = scaleFactor*(1. -  relative_error_scaleFactor)
                scaleFactor_up = scaleFactor*(1. +  relative_error_scaleFactor)
          
                # No systematic uncertainty for the time being
                return scaleFactor, scaleFactor_do, scaleFactor_up, 0.0
            
            elif kindLep == 'mu' :

                RecoSF, RecoSFErr  = self._getTGraphValueEta(self.MuReco,  eta)

                scaleFactor = RecoSF

                if self.BCDEFtoGHRatio != 0.549763 :
                    CorrFactor = 0.995 - (1. - self.BCDEFtoGHRatio)*(0.995 - 0.999)/(1. - 0.549763)
                    scaleFactor = scaleFactor*CorrFactor/0.999

                # This includes systematic errors. Put on statistics because of how it's written AnalysisCMS
                relative_error_scaleFactor = RecoSFErr/RecoSF
                scaleFactor_do = scaleFactor*(1. -  relative_error_scaleFactor)
                scaleFactor_up = scaleFactor*(1. +  relative_error_scaleFactor)
                
                # No systematic uncertainty for the time being
                return scaleFactor, scaleFactor_do, scaleFactor_up, 0.0
                
            # not a lepton ... like some default value: and what can it be if not a lepton? ah ah 
            # --> it happens for default values -9999.
            return 1.0, 1.0, 1.0, 0.0
        
        # FastSim scale factor
        elif ScaleFactorLevel == 2:
            
            if kindLep == 'ele' :
                    
                IdIsoSF, IdIsoSFErr  = self._getHistoValuePtAbsEta(self.FastSimElIdIso,  pt, eta)

                scaleFactor = IdIsoSF

                # This includes systematic errors. Put on statistics because of how it's written AnalysisCMS  
                relative_error_scaleFactor = 0.02
                scaleFactor_do = scaleFactor*(1. -  relative_error_scaleFactor)
                scaleFactor_up = scaleFactor*(1. +  relative_error_scaleFactor)
          
                # No systematic uncertainty for the time being
                return scaleFactor, scaleFactor_do, scaleFactor_up, 0.0
            
            elif kindLep == 'mu' :

                IdSF,   IdSFErr   = self._getHistoValuePtAbsEta(self.FastSimMuId,   pt, eta)
                IsoSF,  IsoSFErr  = self._getHistoValuePtAbsEta(self.FastSimMuIso,  pt, eta)
                IP2DSF, IP2DSFErr = self._getHistoValuePtAbsEta(self.FastSimMuIP2D, pt, eta)
                IP3DSF, IP3DSFErr = self._getHistoValuePtAbsEta(self.FastSimMuIP3D, pt, eta)

                scaleFactor = IdSF*IsoSF*IP2DSF*IP3DSF # ?

                # This includes systematic errors. Put on statistics because of how it's written AnalysisCMS  
                relative_error_scaleFactor = 0.02
                scaleFactor_do = scaleFactor*(1. -  relative_error_scaleFactor)
                scaleFactor_up = scaleFactor*(1. +  relative_error_scaleFactor)

                # No systematic uncertainty for the time being
                return scaleFactor, scaleFactor_do, scaleFactor_up, 0.0
                
            # not a lepton ... like some default value: and what can it be if not a lepton? ah ah 
            # --> it happens for default values -9999.
            return 1.0, 1.0, 1.0, 0.0
                
        # Wrong choice of scale factor level
        print "idisoScaleFactorsStop: Wrong choice of scale factor level for idisoScaleFactorsStop"
        return 1.0, 1.0, 1.0, 0.0


    def process(self,**kwargs):
        tree  = kwargs['tree']
        input = kwargs['input']
        output = kwargs['output']

        self.connect(tree,input)
        
        self.namesOldBranchesToBeModifiedVector = [
           'std_vector_lepton_recoW',
           'std_vector_lepton_recoW_Up',
           'std_vector_lepton_recoW_Down',                              
           'std_vector_lepton_idisoW',
           'std_vector_lepton_idisoW_Up',
           'std_vector_lepton_idisoW_Down',                              
           'std_vector_lepton_idisoW_Syst'                              
           ]

        if self.readfastsim == 1 :
            print 'adding fastsim'
            self.namesNewBranchesVector = ['std_vector_lepton_fastsimW',
                                           'std_vector_lepton_fastsimW_Up',
                                           'std_vector_lepton_fastsimW_Down']
        
            self.clone(output,self.namesOldBranchesToBeModifiedVector + self.namesNewBranchesVector)
            
        else :
            self.clone(output,self.namesOldBranchesToBeModifiedVector)
            
        bvector_reco =  ROOT.std.vector(float) ()
        self.otree.Branch('std_vector_lepton_recoW',bvector_reco)
        bvector_reco_Up =  ROOT.std.vector(float) ()
        self.otree.Branch('std_vector_lepton_recoW_Up',bvector_reco_Up)
        bvector_reco_Down =  ROOT.std.vector(float) ()
        self.otree.Branch('std_vector_lepton_recoW_Down',bvector_reco_Down)

        bvector_idiso =  ROOT.std.vector(float) ()
        self.otree.Branch('std_vector_lepton_idisoW',bvector_idiso)
        bvector_idiso_Up =  ROOT.std.vector(float) ()
        self.otree.Branch('std_vector_lepton_idisoW_Up',bvector_idiso_Up)
        bvector_idiso_Down =  ROOT.std.vector(float) ()
        self.otree.Branch('std_vector_lepton_idisoW_Down',bvector_idiso_Down)
        bvector_idiso_Syst =  ROOT.std.vector(float) ()
        self.otree.Branch('std_vector_lepton_idisoW_Syst',bvector_idiso_Syst)

        bvector_fastsim =  ROOT.std.vector(float) ()
        bvector_fastsim_Up =  ROOT.std.vector(float) ()
        bvector_fastsim_Down =  ROOT.std.vector(float) ()
        if self.readfastsim == 1 :
            self.otree.Branch('std_vector_lepton_fastsimW',bvector_fastsim)
            self.otree.Branch('std_vector_lepton_fastsimW_Up',bvector_fastsim_Up)
            self.otree.Branch('std_vector_lepton_fastsimW_Down',bvector_fastsim_Down)
            
        nentries = self.itree.GetEntries()
        print 'Total number of entries: ',nentries 
        savedentries = 0
                
        # avoid dots to go faster
        itree     = self.itree
        otree     = self.otree

        print '- Starting eventloop'
        step = 5000
        for i in xrange(nentries):
            itree.GetEntry(i)

            ## print event count
            if i > 0 and i%step == 0.:
              print i,'events processed.'

            bvector_reco.clear()
            bvector_reco_Up.clear()
            bvector_reco_Down.clear()

            bvector_idiso.clear()
            bvector_idiso_Up.clear()
            bvector_idiso_Down.clear()
            bvector_idiso_Syst.clear()

            if self.readfastsim == 1 :
                bvector_fastsim.clear()
                bvector_fastsim_Up.clear()
                bvector_fastsim_Down.clear()

            for iLep in xrange(len(itree.std_vector_lepton_pt)) :
             
              pt = itree.std_vector_lepton_pt [iLep]
              eta = itree.std_vector_lepton_eta [iLep]
              flavour = itree.std_vector_lepton_flavour [iLep]
              
              kindLep = 'nonlep' # ele or mu
              if abs (flavour) == 11 : 
                kindLep = 'ele'
              elif abs (flavour) == 13 :
                kindLep = 'mu'
 

              # is tight lepton? 1=tight, 0=loose
              # scale factor level: 0=IdIso, 1=Reco, 2=FastSim
              w, w_lo, w_up, w_syst = self._getWeight(kindLep, pt, eta, 1, 0, itree.nvtx)
              
              bvector_idiso.push_back(w)
              bvector_idiso_Up.push_back(w_up)
              bvector_idiso_Down.push_back(w_lo)             
              bvector_idiso_Syst.push_back(w_syst)             

              w, w_lo, w_up, w_syst = self._getWeight(kindLep, pt, eta, 1, 1, itree.nvtx)
             
              bvector_reco.push_back(w)
              bvector_reco_Up.push_back(w_up)
              bvector_reco_Down.push_back(w_lo)                

              if self.readfastsim == 1:

                  w, w_lo, w_up, w_syst = self._getWeight(kindLep, pt, eta, 1, 2, itree.nvtx)
              
                  bvector_fastsim.push_back(w)
                  bvector_fastsim_Up.push_back(w_up)
                  bvector_fastsim_Down.push_back(w_lo)             
                  

            otree.Fill()
            savedentries+=1

        self.disconnect()
        print '- Eventloop completed'
        print '   Saved: ', savedentries, ' events'


