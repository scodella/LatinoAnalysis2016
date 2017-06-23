#!/usr/bin/env python

from LatinoAnalysis.Gardener.gardening         import gardener_cli
from LatinoAnalysis.Gardener.gardening         import ModuleManager,Pruner,Grafter,AliasGrafter,RootWeighter
#from LatinoAnalysis.Gardener.ww                import WWPruner, WWFlagsGrafter
#from LatinoAnalysis.Gardener.efficiencies      import EffLepFiller,EffTrgFiller

# pileup
from LatinoAnalysis.Gardener.variables.pileup  import PUpper

# trigger efficiencies
from LatinoAnalysis.Gardener.variables.efficiencies               import EffTrgFiller
from LatinoAnalysis.Gardener.variables.triggerMaker               import triggerCalculator
from LatinoAnalysis.Gardener.variables.triggerMaker               import triggerMaker

# id/isolation scale factors
from LatinoAnalysis.Gardener.variables.idisoScaleFactors          import IdIsoSFFiller
from LatinoAnalysis.Gardener.variables.multiIdisoScaleFactors     import MultiIdIsoSFFiller
from LatinoAnalysis.Gardener.variables.idisoScaleFactorsStop      import IdIsoSFStopFiller


# selections
from LatinoAnalysis.Gardener.variables.l2Sel                      import L2SelFiller
from LatinoAnalysis.Gardener.variables.l1Sel                      import L1SelFiller
from LatinoAnalysis.Gardener.variables.LeptonSel                  import LeptonSel
from LatinoAnalysis.Gardener.variables.TauCleaning                import TauCleaning 


# kinematic variables
from LatinoAnalysis.Gardener.variables.l2Kin                      import L2KinFiller
from LatinoAnalysis.Gardener.variables.l3Kin                      import L3KinFiller
from LatinoAnalysis.Gardener.variables.l4Kin                      import L4KinFiller

# fake weights adder for W+jet sample
from LatinoAnalysis.Gardener.variables.fakeWeight                 import FakeWeightFiller
from LatinoAnalysis.Gardener.variables.multiFakeWeight            import multiFakeWeightFiller


# lepton pt corrector
from LatinoAnalysis.Gardener.variables.lepScaleCorrector          import LeptonPtCorrector

# qqWW corrections
from LatinoAnalysis.Gardener.variables.wwNLLcorrectionWeight      import wwNLLcorrectionWeightFiller
from LatinoAnalysis.Gardener.variables.qq2vvEWKcorrectionsWeight  import qq2vvEWKcorrectionsWeightFiller

# qqWZ and qqZZ corrections
from LatinoAnalysis.Gardener.variables.qq2wzEWKcorrectionsWeight  import qq2wzEWKcorrectionsWeightFiller
from LatinoAnalysis.Gardener.variables.qq2zzEWKcorrectionsWeight  import qq2zzEWKcorrectionsWeightFiller

# new variables
from LatinoAnalysis.Gardener.variables.WW2jVar                    import WW2jVarFiller
from LatinoAnalysis.Gardener.variables.WWVar                      import WWVarFiller
from LatinoAnalysis.Gardener.variables.ElectronsVar               import ElectronsVarFiller
from LatinoAnalysis.Gardener.variables.DMVar                      import DMVarFiller
from LatinoAnalysis.Gardener.variables.XWWVar                     import XWWVarFiller
from LatinoAnalysis.Gardener.variables.dymvaVar                   import DymvaVarFiller
from LatinoAnalysis.Gardener.variables.dymvaHiggs                 import DymvaHiggsFiller

from LatinoAnalysis.Gardener.variables.chargeFlipWeight           import chargeFlipWeight
# mucca
from LatinoAnalysis.Gardener.variables.muccaMvaVar                import MuccaMvaVarFiller
from LatinoAnalysis.Gardener.variables.muccaMonoHVar              import MuccaMonoHVarFiller
from LatinoAnalysis.Gardener.variables.muccaMonoHFullVar          import MuccaMonoHFullVarFiller

# mrww
from LatinoAnalysis.Gardener.variables.MrWWVar                    import MrWWVarFiller   

# specific variables for MC
from LatinoAnalysis.Gardener.variables.mcWeights                  import mcWeightsFiller
from LatinoAnalysis.Gardener.variables.mcWeightsCount             import mcWeightsCounter
from LatinoAnalysis.Gardener.variables.GenVar                     import genVariablesFiller
# gen lepton matching
from LatinoAnalysis.Gardener.variables.genMatchVar                import GenMatchVarFiller

#generic tools
from LatinoAnalysis.Gardener.variables.TLorentzVectorCreator      import TLorentzVectorCreator

# filter duplicates in data
from LatinoAnalysis.Gardener.variables.filterDuplicates           import FilterDuplicates

# filter using a JSON in data
from LatinoAnalysis.Gardener.variables.filterJson                 import FilterJSON


# JES uncertainty
from LatinoAnalysis.Gardener.variables.jetScaleUncertainty        import JESTreeMaker

# bpog sfale factors
from LatinoAnalysis.Gardener.variables.btagPogScaleFactors        import btagPogScaleFactors
from LatinoAnalysis.Gardener.variables.allBtagPogScaleFactors     import allBtagPogScaleFactors
from LatinoAnalysis.Gardener.variables.allBtagPogScaleFactorsICHEP     import allBtagPogScaleFactorsICHEP
# lepton pT scale uncertainty and resolution
from LatinoAnalysis.Gardener.variables.lepScaleUncertainty        import LeppTScalerTreeMaker
from LatinoAnalysis.Gardener.variables.lepResolutionUncertainty   import LeptonResolutionTreeMaker
# MET uncertainty
from LatinoAnalysis.Gardener.variables.metUncertainty             import MetUncertaintyTreeMaker
from LatinoAnalysis.Gardener.variables.metUnclustered             import MetUnclusteredTreeMaker
from LatinoAnalysis.Gardener.variables.metXYshift                 import MetXYshiftTreeMaker
# QCD uncertainty
from LatinoAnalysis.Gardener.variables.qcdUncertainty             import QcdUncertaintyTreeMaker
# PDF uncertainty
from LatinoAnalysis.Gardener.variables.pdfUncertainty             import PdfUncertaintyTreeMaker
# EWK singlet reweighter
from LatinoAnalysis.Gardener.variables.BWEwkSingletReweighter     import BWEwkSingletReweighter
# PDF and scale uncertainty
from LatinoAnalysis.Gardener.variables.pdfAndScaleUncertainty     import PdfAndScaleUncertaintyTreeMaker
# GenPT for the top
from LatinoAnalysis.Gardener.variables.TopGenPt                   import TopGenPt

# ggH uncertainty LHCXSWG
from LatinoAnalysis.Gardener.variables.ggHUncertainty             import ggHUncertaintyMaker

# generic formula adder
from LatinoAnalysis.Gardener.variables.genericFormulaAdder        import genericFormulaAdder

if __name__ == '__main__':

    print "gardener"
    
    modules = ModuleManager()
    modules['filter']           = Pruner()
    modules['adder']            = Grafter()
    modules['alias']            = AliasGrafter()
    modules['rootweighter']     = RootWeighter()
    #modules['wwfilter']         = WWPruner()
    #modules['wwflagger']        = WWFlagsGrafter()
    modules['puadder']          = PUpper()


# filter duplicates
    modules['filterduplicates']          = FilterDuplicates()

# filter using a json file
    modules['filterjson']          = FilterJSON()


# trigger efficiency
    modules['efftfiller']       = EffTrgFiller()
    modules['trigMaker']        = triggerMaker()

# id/isolation scale factors
    modules['idisofiller'] = IdIsoSFFiller()
    modules['multiidiso']  = MultiIdIsoSFFiller()
    modules['idisostopfiller'] = IdIsoSFStopFiller()

# specific variables for MC

    modules['mcweightscounter'] = mcWeightsCounter()
    modules['mcweightsfiller']  = mcWeightsFiller()

# generator level variables
    modules['genvariablesfiller']  = genVariablesFiller()
# gen lepton matching
    modules['genmatchvarfiller']  = GenMatchVarFiller()




# new variables
    modules['ww2jvarfiller']    = WW2jVarFiller()
    modules['wwvarfiller']      = WWVarFiller()
    modules['electronidfiller'] = ElectronsVarFiller()
    modules['dmvarfiller']      = DMVarFiller()
    modules['xwwvarfiller']     = XWWVarFiller()
    modules['dymvaVarFiller']   = DymvaVarFiller()
    modules['dymvaHiggsFiller']   = DymvaHiggsFiller()

# Charge Flip
    modules['chFlipProba']      = chargeFlipWeight()

# mucca
    modules['muccaMvaVarFiller']       = MuccaMvaVarFiller()
    modules['muccaMonoHVarFiller']     = MuccaMonoHVarFiller()
    modules['muccaMonoHFullVarFiller'] = MuccaMonoHFullVarFiller()

# mrWW
    modules['mrWWvarfiller']   = MrWWVarFiller()



# add nll re-weight for ww
    modules['wwNLLcorrections']      =  wwNLLcorrectionWeightFiller()
# add ewk re-weight for ww
    modules['wwEWKcorrections']      =  qq2vvEWKcorrectionsWeightFiller()

# add ewk re-weight for wz and zz
    modules['wzEWKcorrections']      =  qq2wzEWKcorrectionsWeightFiller()
    modules['zzEWKcorrections']      =  qq2zzEWKcorrectionsWeightFiller()

# add bpog SF
    modules['btagPogScaleFactors']   = btagPogScaleFactors()
    modules['allBtagPogScaleFactors'] = allBtagPogScaleFactors()
    modules['allBtagPogScaleFactorsICHEP'] = allBtagPogScaleFactorsICHEP()

# generic tool
    modules['tlorentzvectorfiller']  = TLorentzVectorCreator()

# apply selections and update variables
    modules['l2selfiller']     = L2SelFiller()
    modules['l1selfiller']     = L1SelFiller()
    modules['lepSel']          = LeptonSel()
    modules['cleanTau']        = TauCleaning()

# update kinematic variables
    modules['l2kinfiller']     = L2KinFiller()
    modules['l3kinfiller']     = L3KinFiller()
    modules['l4kinfiller']     = L4KinFiller()

# Nuisances
    modules['JESTreeMaker']           = JESTreeMaker()
    modules['LeppTScalerTreeMaker']   = LeppTScalerTreeMaker()
    modules['leptonResolution']       = LeptonResolutionTreeMaker()
    modules['metUncertainty']         = MetUncertaintyTreeMaker()
    modules['metUnclustered']         = MetUnclusteredTreeMaker()
    modules['pdfUncertainty']         = PdfUncertaintyTreeMaker()
    modules['qcdUncertainty']         = QcdUncertaintyTreeMaker()
    modules['pdfAndScaleUncertainty'] = PdfAndScaleUncertaintyTreeMaker()
    
    
# fake weights
    modules['fakeWeights']      = FakeWeightFiller()
    modules['multiFakeWeights']      = multiFakeWeightFiller()

# lepton pt corrector
    modules['letPtCorrector']   = LeptonPtCorrector()
# MET corrector
    modules['metXYshift']     = MetXYshiftTreeMaker()

# EWK bw weights
    modules['BWEwkSingletReweighter'] = BWEwkSingletReweighter()

# Top gen pt for reweighting
    modules['TopGenPt'] = TopGenPt()
    
# ggH uncertainty LHCXSWG
    modules['ggHUncertainty'] = ggHUncertaintyMaker()

# generic formula adder
    modules['genericFormulaAdder'] = genericFormulaAdder()



    gardener_cli( modules )


