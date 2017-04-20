
eosProdBase= '/eos/cms/'
eosTargBaseIn = '/eos/user/x/xjanssen/HWW2015/'
eosTargBaseOut= '/eos/user/x/xjanssen/HWW2015/'


# ---- production to run on
# .... .... this is defined by mkGardener in "-p" "--prods" option

Productions= {

#### 74x / 21Oct & 21OctBis tags / miniAOD v2

  '21Oct_25ns_MC'   : {
                        'isData'  : False ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_spring15_miniaodv2_25ns.py' , 
                      #  'dirExt'  : 'LatinoTrees' ,
                        'dirExt'  : 'split' ,
                        'gDocID'  : '1wH73CYA_T4KMkl1Cw-xLTj8YG7OPqayDnP53N-lZwFQ' ,
                        #'puData'  : '/afs/cern.ch/user/x/xjanssen/public/MyDataPileupHistogram.root',
                        'puData'  : '/afs/cern.ch/user/x/xjanssen/public/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_v2_from256630_PileupHistogram.root' ,
                        #'bigSamples': ['DYJetsToLL_M-10to50'] ,
                        'cmssw' : '74x' ,
                      } ,

  '21Oct_Run2015D_05Oct2015' : {
                        'isData'  : True ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_dataD_05Oct2015_25ns.py' ,
                        #'dirExt'  : 'Run2015D_05Oct2015' ,
                        'dirExt'  : 'split' ,
                        'cmssw' : '74x' ,
                      } ,

  '21Oct_Run2015D_PromptReco' : {
                        'isData'  : True ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_dataD_PromptReco_25ns.py' ,
                        #'dirExt'  : 'Run2015D_PromptReco' ,
                        'dirExt'  : 'split' ,
                        'cmssw' : '74x' ,
                      } ,

  '21OctBis_Run2015D_05Oct2015' : {
                        'isData'  : True ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_dataD_05Oct2015_25ns_21OctBis.py',
                        #'dirExt'  : 'Run2015D_PromptReco' ,
                        'dirExt'  : 'split2' ,
                        'cmssw' : '74x' ,
                      } ,

  '21OctBis_Run2015D_PromptReco_0716pb' : {
                        'isData'  : True ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_dataD_PromptReco_25ns_21OctBis_0716pb.py' ,
                        'dirExt'  : 'split2' ,
                        'cmssw' : '74x' ,
                      } ,

  '21OctBis_Run2015D_PromptReco_0851pb' : {
                        'isData'  : True ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_dataD_PromptReco_25ns_21OctBis_0851pb.py' ,
                        'dirExt'  : 'split2' ,
                        'cmssw' : '74x' ,
                      } ,



#### 76x / StarWars tag / miniAOD v1

  '08Jan_25ns_mAODv1_MC'   : {
                        'isData'  : False ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_fall15_miniaod_25ns.py' ,
                        'dirExt'  : 'LatinoTrees' ,
                        'gDocID'  : '1wH73CYA_T4KMkl1Cw-xLTj8YG7OPqayDnP53N-lZwFQ' ,
                        'puData'  : '/afs/cern.ch/user/x/xjanssen/public/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_v2_from256630_PileupHistogram.root' ,
                        'cmssw'   : '763' ,
                      } ,

#### 76x / StarWars tag (v2) / miniAOD v2

  '18Jan_25ns_mAODv2_MC_TEST'   : {
                        'isData'  : False ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_fall15_miniaodv2_25ns.py' ,
                        'dir'     : '/store/group/phys_higgs/cmshww/amassiro/RunII/18Jan/MC/25ns/',
                        'dirExt'  : 'LatinoTrees' ,
                        #'gDocID'  : '1wH73CYA_T4KMkl1Cw-xLTj8YG7OPqayDnP53N-lZwFQ' ,
                        'puData'  : '/afs/cern.ch/user/x/xjanssen/public/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_v2_from256630_PileupHistogram.root' ,
                        'cmssw'   : '763' ,
                      } ,


#### 76x / StarWars tag (v3) / miniAOD v2

  '22Jan_25ns_mAODv2_MC'   : {
                        'isData'  : False ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_fall15_miniaodv2_25ns.py' ,
                        'dirExt'  : 'LatinoTrees' ,
                        #'gDocID'  : '1wH73CYA_T4KMkl1Cw-xLTj8YG7OPqayDnP53N-lZwFQ' ,
                        'puData'  : '/afs/cern.ch/user/x/xjanssen/public/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_v2_from256630_PileupHistogram.root' ,
                        'cmssw'   : '763' ,
                      } ,

  '22Jan_Run2015D_16Dec2015' : {
                        'isData'  : True ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_dataD_16Dec2015-v1_25ns_StarWars.py' ,
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : '763' ,
                        'reName'  : { 
                                       'DoubleEG_Run2015D_25ns-16Dec2015-v2'       : 'Run2015D_16Dec2015_DoubleEG' ,
                                       'DoubleMuon_Run2015D_25ns-16Dec2015-v1'     : 'Run2015D_16Dec2015_DoubleMuon' ,
                                       'MuonEG_Run2015D_25ns-16Dec2015-v1'         : 'Run2015D_16Dec2015_MuonEG' , 
                                       'SingleElectron_Run2015D_25ns-16Dec2015-v1' : 'Run2015D_16Dec2015_SingleElectron' ,
                                       'SingleMuon_Run2015D_25ns-16Dec2015-v1'     : 'Run2015D_16Dec2015_SingleMuon' ,
                                      #'MET_Run2015D_25ns-16Dec2015-v1'            : 'Run2015D_16Dec2015_MET' ,
                                      #'SinglePhoton_Run2015D_25ns-16Dec2015-v1'   : 'Run2015D_16Dec2015_SinglePhoton' ,
                                    }
                      } ,

  '22Jan_Run2015C_16Dec2015' : {
                        'isData'  : True ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_dataC_16Dec2015-v1_25ns_StarWars.py',
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : '763' ,
                        'reName'  : {
                                       'DoubleEG_Run2015C_25ns-16Dec2015-v1'       : 'Run2015C_16Dec2015_DoubleEG' ,
                                       'DoubleMuon_Run2015C_25ns-16Dec2015-v1'     : 'Run2015C_16Dec2015_DoubleMuon' ,
                                       'MuonEG_Run2015C_25ns-16Dec2015-v1'         : 'Run2015C_16Dec2015_MuonEG' ,
                                       'SingleElectron_Run2015C_25ns-16Dec2015-v1' : 'Run2015C_16Dec2015_SingleElectron' ,
                                       'SingleMuon_Run2015C_25ns-16Dec2015-v1'     : 'Run2015C_16Dec2015_SingleMuon' ,
                                   }
                      } ,

  '03Mar_25ns_mAODv2_MC'   : {
                        'isData'  : False ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_fall15_miniaodv2_25ns.py' ,
                        'dir'     : '/store/group/phys_higgs/cmshww/amassiro/RunII/03Mar/MC/25ns/',
                        'dirExt'  : 'LatinoTrees' ,
                        #'gDocID'  : '1wH73CYA_T4KMkl1Cw-xLTj8YG7OPqayDnP53N-lZwFQ' ,
                        'puData'  : '/afs/cern.ch/user/x/xjanssen/public/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_v2_from256630_PileupHistogram.root' ,
                        'cmssw'   : '763' ,
                      } ,


  '03Mar_Run2015D_16Dec2015' : {
                        'isData'  : True ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_dataD_16Dec2015-v1_25ns_StarWars.py' ,
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : '763' ,
                        'reName'  : {
                                       'DoubleEG_Run2015D_25ns-16Dec2015-v2'       : 'Run2015D_16Dec2015_DoubleEG' ,
                                       'DoubleMuon_Run2015D_25ns-16Dec2015-v1'     : 'Run2015D_16Dec2015_DoubleMuon' ,
                                       'MuonEG_Run2015D_25ns-16Dec2015-v1'         : 'Run2015D_16Dec2015_MuonEG' ,
                                       'SingleElectron_Run2015D_25ns-16Dec2015-v1' : 'Run2015D_16Dec2015_SingleElectron' ,
                                       'SingleMuon_Run2015D_25ns-16Dec2015-v1'     : 'Run2015D_16Dec2015_SingleMuon' ,
                                      #'MET_Run2015D_25ns-16Dec2015-v1'            : 'Run2015D_16Dec2015_MET' ,
                                      #'SinglePhoton_Run2015D_25ns-16Dec2015-v1'   : 'Run2015D_16Dec2015_SinglePhoton' ,
                                    }
                      } ,

  '03Mar_Run2015C_16Dec2015' : {
                        'isData'  : True ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_dataC_16Dec2015-v1_25ns_StarWars.py',
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : '763' ,
                        'reName'  : {
                                       'DoubleEG_Run2015C_25ns-16Dec2015-v1'       : 'Run2015C_16Dec2015_DoubleEG' ,
                                       'DoubleMuon_Run2015C_25ns-16Dec2015-v1'     : 'Run2015C_16Dec2015_DoubleMuon' ,
                                       'MuonEG_Run2015C_25ns-16Dec2015-v1'         : 'Run2015C_16Dec2015_MuonEG' ,
                                       'SingleElectron_Run2015C_25ns-16Dec2015-v1' : 'Run2015C_16Dec2015_SingleElectron' ,
                                       'SingleMuon_Run2015C_25ns-16Dec2015-v1'     : 'Run2015C_16Dec2015_SingleMuon' ,
                                   }
                      } ,

# 80X 2016 DAT/MC --> REM: Still using 763 post-processing  --> MoneyMonster tag




  '13May2016_25ns_Spring16_MC'   : {
                        'isData'  : False ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_spring16_miniaod.py' ,
                        'dir'     : '/store/group/phys_higgs/cmshww/amassiro/RunII/2016/May13/MC/25ns/',
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : '763' ,
                        #'gDocID'  : '1wH73CYA_T4KMkl1Cw-xLTj8YG7OPqayDnP53N-lZwFQ' ,
                        #'puData'  : '/afs/cern.ch/user/x/xjanssen/public/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_v2_from256630_PileupHistogram.root' ,
                      } ,

  '20May2016_Run2016B_PromptReco' : {
                        'isData'  : True ,        
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_data_2016_PromptReco.py',
                        'dir'     : '/store/group/phys_higgs/cmshww/amassiro/RunII/2016/May20/data/25ns/',
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : '763' ,
                        'reName'  : {
                                       'DoubleEG_Run2016B-PromptReco-v2'          : 'Run2016B_PromptReco_DoubleEG',
                                       'DoubleMuon_Run2016B-PromptReco-v2'        : 'Run2016B_PromptReco_DoubleMuon',
                                       'MuonEG_Run2016B-PromptReco-v2'            : 'Run2016B_PromptReco_MuonEG',
                                       'SingleElectron_Run2016B-PromptReco-v2'    : 'Run2016B_PromptReco_SingleElectron',   
                                       'SingleMuon_Run2016B-PromptReco-v2'        : 'Run2016B_PromptReco_SingleMuon',
                                    },
                     },


# 80X 2016 DATA/MC: 3June2016_NinjaTurtles_v3 Tag / mAODv2

  '07Jun2016_spring16_mAODv2_0p8fbm1' : {
                        'isData'  : False ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_spring16_miniaod_v2.py' ,
                        'dir'     : '/store/group/phys_higgs/cmshww/amassiro/RunII/2016/Jun07/MC/v2/',
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'ICHEP2016' ,
                        # 0.8 fb-1
                        'puData'  : '/afs/cern.ch/user/x/xjanssen/public/PileupHistogram_805_15Jun.root' , 
                       },


  '07Jun2016_spring16_mAODv2'   : {
                        'isData'  : False ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_spring16_miniaod_v2.py' ,
                        'dir'     : '/store/group/phys_higgs/cmshww/amassiro/RunII/2016/Jun07/MC/v2/',
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'ICHEP2016' ,
                        # 2.6 fb-1
                        'puData'  : '/afs/cern.ch/user/x/xjanssen/public/PileupHistogram_805_28Jun.root' , 
                       },

  '07Jun2016_spring16_mAODv2_4p0fbm1'   : {
                        'isData'  : False ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_spring16_miniaod_v2.py' ,
                        'dir'     : '/store/group/phys_higgs/cmshww/amassiro/RunII/2016/Jun07/MC/v2/',
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'ICHEP2016' ,
                        # 4.0 fb-1
                        'puData'  : '/afs/cern.ch/user/x/xjanssen/public/PileupHistogram_805_5Jul_4fb.root' ,
                       },


  '07Jun2016_spring16_mAODv2_6p3fbm1'   : {
                        'isData'  : False ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_spring16_miniaod_v2.py' ,
                        'dir'     : '/store/group/phys_higgs/cmshww/amassiro/RunII/2016/Jun07/MC/v2/',
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'ICHEP2016' ,
                        # 6.3 fb-1
                        'puData'  : '/afs/cern.ch/user/x/xjanssen/public/PileupHistogram_805_8Jul_6p3fb.root' ,
                       },

  '07Jun2016_spring16_mAODv2_12pXfbm1'   : {
                        'isData'  : False ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_spring16_miniaod_v2.py' ,
                        'dir'     : '/store/group/phys_higgs/cmshww/amassiro/RunII/2016/Jun07/MC/v2/',
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'ICHEP2016' ,
                        # 12.X fb-1
                        'puData'  : '/afs/cern.ch/user/x/xjanssen/public/PileupHistogram_805_271036-276811_63mb_24Aug.root',
                       },

  '07Jun2016_spring16_mAODv2_12pXfbm1_KNU'   : {# To be used @ KNU T3 farm @@@@@@@@@@@@@@@@@@@22
                        'isData'  : False ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_spring16_miniaod_v2.py' ,
                        'dir'     : 'HWW12fb_v2/07Jun2016_spring16_mAODv2_12pXfbm1/',
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'ICHEP2016' ,
                        # 12.X fb-1
                        'puData'  : '/afs/cern.ch/user/x/xjanssen/public/PileupHistogram_805_271036-276811_63mb_24Aug.root',
                       },

  '07Jun2016_spring16_mAODv2_12pXfbm1_repro'   : {
                        'isData'  : False ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_spring16_miniaod_v2.py' ,
                        'dir'     : '/store/group/phys_higgs/cmshww/amassiro/RunII/2016/Jun07/MC/v2/',
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'ICHEP2016' ,
                        # 12.X fb-1
                        'puData'  : '/afs/cern.ch/user/x/xjanssen/public/PileupHistogram_805_271036-276811_63mb_24Aug.root',
                       },

  '07Jun2016_spring16_mAODv2_12pXfbm1_demo'   : {
                        'isData'  : False ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_spring16_miniaod_v2.py' ,
                        'dir'     : '/store/group/phys_higgs/cmshww/amassiro/RunII/2016/Jun07/MC/demo/',
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'ICHEP2016' ,
                        # 12.X fb-1
                        'puData'  : '/afs/cern.ch/user/x/xjanssen/public/PileupHistogram_805_271036-276811_63mb_24Aug.root',
                       },

  '07Jun2016_Run2016B_PromptReco' : {
                        'isData'  : True ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_data_2016_PromptReco.py',
                        'dir'     : '/store/group/phys_higgs/cmshww/amassiro/RunII/2016/Jun07/data/25ns/',
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'ICHEP2016' ,
                        'reName'  : {
                                       'DoubleEG_Run2016B-PromptReco-v2'          : 'Run2016B_PromptReco_DoubleEG',
                                       'DoubleMuon_Run2016B-PromptReco-v2'        : 'Run2016B_PromptReco_DoubleMuon',
                                       'MuonEG_Run2016B-PromptReco-v2'            : 'Run2016B_PromptReco_MuonEG',
                                       'SingleElectron_Run2016B-PromptReco-v2'    : 'Run2016B_PromptReco_SingleElectron',
                                       'SingleMuon_Run2016B-PromptReco-v2'        : 'Run2016B_PromptReco_SingleMuon',
                                    },
                       },
            
  '21Jun2016_Run2016B_PromptReco' : {
                        'isData'  : True ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_data_2016_PromptReco.py',
                        'dir'     : '/store/group/phys_higgs/cmshww/amassiro/RunII/2016/Jun21/data/25ns/',
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'ICHEP2016' ,
                        'reName'  : {
                                       'DoubleEG_Run2016B-PromptReco-v2'          : 'Run2016B_PromptReco_DoubleEG',
                                       'DoubleMuon_Run2016B-PromptReco-v2'        : 'Run2016B_PromptReco_DoubleMuon',
                                       'MuonEG_Run2016B-PromptReco-v2'            : 'Run2016B_PromptReco_MuonEG',
                                       'SingleElectron_Run2016B-PromptReco-v2'    : 'Run2016B_PromptReco_SingleElectron',
                                       'SingleMuon_Run2016B-PromptReco-v2'        : 'Run2016B_PromptReco_SingleMuon',
                                    },
                       },

  '21Jun2016_v2_Run2016B_PromptReco' : {
                        'isData'  : True ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_data_2016_PromptReco.py',
                        'dir'     : '/store/group/phys_higgs/cmshww/amassiro/RunII/2016/Jun21_v2/data/25ns/',
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'ICHEP2016' ,
                        'reName'  : {
                                       'DoubleEG_Run2016B-PromptReco-v2'          : 'Run2016B_PromptReco_DoubleEG',
                                       'DoubleMuon_Run2016B-PromptReco-v2'        : 'Run2016B_PromptReco_DoubleMuon',
                                       'MuonEG_Run2016B-PromptReco-v2'            : 'Run2016B_PromptReco_MuonEG',
                                       'SingleElectron_Run2016B-PromptReco-v2'    : 'Run2016B_PromptReco_SingleElectron',
                                       'SingleMuon_Run2016B-PromptReco-v2'        : 'Run2016B_PromptReco_SingleMuon',
                                    },
                       },
  '21Jun2016_v2_Run2016B_PromptReco_KNU' : {
                        'isData'  : True ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_data_2016_PromptReco.py',
                        'dir'     : 'HWW6p3/21Jun2016_v2_Run2016B_PromptReco/',
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'ICHEP2016' ,
                        'reName'  : {
                                       'DoubleEG_Run2016B-PromptReco-v2'          : 'Run2016B_PromptReco_DoubleEG',
                                       'DoubleMuon_Run2016B-PromptReco-v2'        : 'Run2016B_PromptReco_DoubleMuon',
                                       'MuonEG_Run2016B-PromptReco-v2'            : 'Run2016B_PromptReco_MuonEG',
                                       'SingleElectron_Run2016B-PromptReco-v2'    : 'Run2016B_PromptReco_SingleElectron',
                                       'SingleMuon_Run2016B-PromptReco-v2'        : 'Run2016B_PromptReco_SingleMuon',
                                    },
                       },


  '21Jun2016_v2_Run2016B_PromptReco_repro' : {
                        'isData'  : True ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_data_2016_PromptReco.py',
                        'dir'     : '/store/group/phys_higgs/cmshww/amassiro/RunII/2016/Jun21_v2/data/25ns/',
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'ICHEP2016' ,
                        'reName'  : {
                                       'DoubleEG_Run2016B-PromptReco-v2'          : 'Run2016B_PromptReco_DoubleEG',
                                       'DoubleMuon_Run2016B-PromptReco-v2'        : 'Run2016B_PromptReco_DoubleMuon',
                                       'MuonEG_Run2016B-PromptReco-v2'            : 'Run2016B_PromptReco_MuonEG',
                                       'SingleElectron_Run2016B-PromptReco-v2'    : 'Run2016B_PromptReco_SingleElectron',
                                       'SingleMuon_Run2016B-PromptReco-v2'        : 'Run2016B_PromptReco_SingleMuon',
                                    },
                       },


  '05Jul2016_Run2016B_PromptReco' : {
                        'isData'  : True ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_data_2016_PromptReco.py',
                        'dir'     : '/store/group/phys_higgs/cmshww/amassiro/RunII/2016/Jul05/data/25ns/',
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'ICHEP2016' ,
                        'reName'  : {
                                       'DoubleEG_Run2016B-PromptReco-v2'          : 'Run2016B_PromptReco_DoubleEG',
                                       'DoubleMuon_Run2016B-PromptReco-v2'        : 'Run2016B_PromptReco_DoubleMuon',
                                       'MuonEG_Run2016B-PromptReco-v2'            : 'Run2016B_PromptReco_MuonEG',
                                       'SingleElectron_Run2016B-PromptReco-v2'    : 'Run2016B_PromptReco_SingleElectron',
                                       'SingleMuon_Run2016B-PromptReco-v2'        : 'Run2016B_PromptReco_SingleMuon',
                                    },
                       },
  '05Jul2016_Run2016B_PromptReco_KNU' : {
                        'isData'  : True ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_data_2016_PromptReco.py',
                        'dir'     : 'HWW6p3/05Jul2016_Run2016B_PromptReco/',
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'ICHEP2016' ,
                        'reName'  : {
                                       'DoubleEG_Run2016B-PromptReco-v2'          : 'Run2016B_PromptReco_DoubleEG',
                                       'DoubleMuon_Run2016B-PromptReco-v2'        : 'Run2016B_PromptReco_DoubleMuon',
                                       'MuonEG_Run2016B-PromptReco-v2'            : 'Run2016B_PromptReco_MuonEG',
                                       'SingleElectron_Run2016B-PromptReco-v2'    : 'Run2016B_PromptReco_SingleElectron',
                                       'SingleMuon_Run2016B-PromptReco-v2'        : 'Run2016B_PromptReco_SingleMuon',
                                    },
                       },

  '05Jul2016_Run2016B_PromptReco_repro' : {
                        'isData'  : True ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_data_2016_PromptReco.py',
                        'dir'     : '/store/group/phys_higgs/cmshww/amassiro/RunII/2016/Jul05/data/25ns/',
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'ICHEP2016' ,
                        'reName'  : {
                                       'DoubleEG_Run2016B-PromptReco-v2'          : 'Run2016B_PromptReco_DoubleEG',
                                       'DoubleMuon_Run2016B-PromptReco-v2'        : 'Run2016B_PromptReco_DoubleMuon',
                                       'MuonEG_Run2016B-PromptReco-v2'            : 'Run2016B_PromptReco_MuonEG',
                                       'SingleElectron_Run2016B-PromptReco-v2'    : 'Run2016B_PromptReco_SingleElectron',
                                       'SingleMuon_Run2016B-PromptReco-v2'        : 'Run2016B_PromptReco_SingleMuon',
                                    },
                       },


  '08Jul2016_Run2016B_PromptReco' : {
                        'isData'  : True ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_data_2016_PromptReco.py',
                        'dir'     : '/store/group/phys_higgs/cmshww/amassiro/RunII/2016/Jul08/data/25ns/',
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'ICHEP2016' ,
                        'reName'  : {
                                       # Run2016B
                                       'DoubleEG_Run2016B-PromptReco-v2'          : 'Run2016B_PromptReco_DoubleEG',
                                       'DoubleMuon_Run2016B-PromptReco-v2'        : 'Run2016B_PromptReco_DoubleMuon',
                                       'MuonEG_Run2016B-PromptReco-v2'            : 'Run2016B_PromptReco_MuonEG',
                                       'SingleElectron_Run2016B-PromptReco-v2'    : 'Run2016B_PromptReco_SingleElectron',
                                       'SingleMuon_Run2016B-PromptReco-v2'        : 'Run2016B_PromptReco_SingleMuon',
                                    },
                       },
  '08Jul2016_Run2016B_PromptReco_KNU' : {
                        'isData'  : True ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_data_2016_PromptReco.py',
                        'dir'     : 'HWW6p3/08Jul2016_Run2016B_PromptReco/',
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'ICHEP2016' ,
                        'reName'  : {
                                       # Run2016B
                                       'DoubleEG_Run2016B-PromptReco-v2'          : 'Run2016B_PromptReco_DoubleEG',
                                       'DoubleMuon_Run2016B-PromptReco-v2'        : 'Run2016B_PromptReco_DoubleMuon',
                                       'MuonEG_Run2016B-PromptReco-v2'            : 'Run2016B_PromptReco_MuonEG',
                                       'SingleElectron_Run2016B-PromptReco-v2'    : 'Run2016B_PromptReco_SingleElectron',
                                       'SingleMuon_Run2016B-PromptReco-v2'        : 'Run2016B_PromptReco_SingleMuon',
                                    },
                       },

  '08Jul2016_Run2016B_PromptReco_repro' : {
                        'isData'  : True ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_data_2016_PromptReco.py',
                        'dir'     : '/store/group/phys_higgs/cmshww/amassiro/RunII/2016/Jul08/data/25ns/',
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'ICHEP2016' ,
                        'reName'  : {
                                       # Run2016B
                                       'DoubleEG_Run2016B-PromptReco-v2'          : 'Run2016B_PromptReco_DoubleEG',
                                       'DoubleMuon_Run2016B-PromptReco-v2'        : 'Run2016B_PromptReco_DoubleMuon',
                                       'MuonEG_Run2016B-PromptReco-v2'            : 'Run2016B_PromptReco_MuonEG',
                                       'SingleElectron_Run2016B-PromptReco-v2'    : 'Run2016B_PromptReco_SingleElectron',
                                       'SingleMuon_Run2016B-PromptReco-v2'        : 'Run2016B_PromptReco_SingleMuon',
                                    },
                       },


  '08Jul2016_Run2016C_PromptReco' : {
                        'isData'  : True ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_data_2016_PromptReco.py',
                        'dir'     : '/store/group/phys_higgs/cmshww/amassiro/RunII/2016/Jul08/data/25ns/',
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'ICHEP2016' ,
                        'reName'  : {
                                       # Run2016C
                                       'DoubleEG_Run2016C-PromptReco-v2'          : 'Run2016C_PromptReco_DoubleEG',
                                       'DoubleMuon_Run2016C-PromptReco-v2'        : 'Run2016C_PromptReco_DoubleMuon',
                                       'MuonEG_Run2016C-PromptReco-v2'            : 'Run2016C_PromptReco_MuonEG',
                                       'SingleElectron_Run2016C-PromptReco-v2'    : 'Run2016C_PromptReco_SingleElectron',
                                       'SingleMuon_Run2016C-PromptReco-v2'        : 'Run2016C_PromptReco_SingleMuon',
                                    },
                       },
  '08Jul2016_Run2016C_PromptReco_KNU' : {
                        'isData'  : True ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_data_2016_PromptReco.py',
                        'dir'     : 'HWW6p3/08Jul2016_Run2016C_PromptReco/',
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'ICHEP2016' ,
                        'reName'  : {
                                       # Run2016C
                                       'DoubleEG_Run2016C-PromptReco-v2'          : 'Run2016C_PromptReco_DoubleEG',
                                       'DoubleMuon_Run2016C-PromptReco-v2'        : 'Run2016C_PromptReco_DoubleMuon',
                                       'MuonEG_Run2016C-PromptReco-v2'            : 'Run2016C_PromptReco_MuonEG',
                                       'SingleElectron_Run2016C-PromptReco-v2'    : 'Run2016C_PromptReco_SingleElectron',
                                       'SingleMuon_Run2016C-PromptReco-v2'        : 'Run2016C_PromptReco_SingleMuon',
                                    },
                       },

  '08Jul2016_Run2016C_PromptReco_repro' : {
                        'isData'  : True ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_data_2016_PromptReco.py',
                        'dir'     : '/store/group/phys_higgs/cmshww/amassiro/RunII/2016/Jul08/data/25ns/',
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'ICHEP2016' ,
                        'reName'  : {
                                       # Run2016C
                                       'DoubleEG_Run2016C-PromptReco-v2'          : 'Run2016C_PromptReco_DoubleEG',
                                       'DoubleMuon_Run2016C-PromptReco-v2'        : 'Run2016C_PromptReco_DoubleMuon',
                                       'MuonEG_Run2016C-PromptReco-v2'            : 'Run2016C_PromptReco_MuonEG',
                                       'SingleElectron_Run2016C-PromptReco-v2'    : 'Run2016C_PromptReco_SingleElectron',
                                       'SingleMuon_Run2016C-PromptReco-v2'        : 'Run2016C_PromptReco_SingleMuon',
                                    },
                       },


  '11Jul2016_Run2016C_PromptReco' : {
                        'isData'  : True ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_data_2016_PromptReco.py',
                        'dir'     : '/store/group/phys_higgs/cmshww/amassiro/RunII/2016/Jul11_NoL1T/data/25ns/',
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'ICHEP2016' ,
                        'reName'  : {
                                       # Run2016C
                                       'DoubleEG_Run2016C-PromptReco-v2'          : 'Run2016C_PromptReco_DoubleEG',
                                       'DoubleMuon_Run2016C-PromptReco-v2'        : 'Run2016C_PromptReco_DoubleMuon',
                                       'MuonEG_Run2016C-PromptReco-v2'            : 'Run2016C_PromptReco_MuonEG',
                                       'SingleElectron_Run2016C-PromptReco-v2'    : 'Run2016C_PromptReco_SingleElectron',
                                       'SingleMuon_Run2016C-PromptReco-v2'        : 'Run2016C_PromptReco_SingleMuon',
                                    },
                       },
  '11Jul2016_Run2016C_PromptReco_KNU' : {
                        'isData'  : True ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_data_2016_PromptReco.py',
                        'dir'     : 'HWW12fb/11Jul2016_Run2016C_PromptReco/',
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'ICHEP2016' ,
                        'reName'  : {
                                       # Run2016C
                                       'DoubleEG_Run2016C-PromptReco-v2'          : 'Run2016C_PromptReco_DoubleEG',
                                       'DoubleMuon_Run2016C-PromptReco-v2'        : 'Run2016C_PromptReco_DoubleMuon',
                                       'MuonEG_Run2016C-PromptReco-v2'            : 'Run2016C_PromptReco_MuonEG',
                                       'SingleElectron_Run2016C-PromptReco-v2'    : 'Run2016C_PromptReco_SingleElectron',
                                       'SingleMuon_Run2016C-PromptReco-v2'        : 'Run2016C_PromptReco_SingleMuon',
                                    },
                       },

  '11Jul2016_Run2016C_PromptReco_repro' : {
                        'isData'  : True ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_data_2016_PromptReco.py',
                        'dir'     : '/store/group/phys_higgs/cmshww/amassiro/RunII/2016/Jul11_NoL1T/data/25ns/',
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'ICHEP2016' ,
                        'reName'  : {
                                       # Run2016C
                                       'DoubleEG_Run2016C-PromptReco-v2'          : 'Run2016C_PromptReco_DoubleEG',
                                       'DoubleMuon_Run2016C-PromptReco-v2'        : 'Run2016C_PromptReco_DoubleMuon',
                                       'MuonEG_Run2016C-PromptReco-v2'            : 'Run2016C_PromptReco_MuonEG',
                                       'SingleElectron_Run2016C-PromptReco-v2'    : 'Run2016C_PromptReco_SingleElectron',
                                       'SingleMuon_Run2016C-PromptReco-v2'        : 'Run2016C_PromptReco_SingleMuon',
                                    },
                       },


  '15Jul2016_Run2016C_PromptReco' : {
                        'isData'  : True ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_data_2016_PromptReco.py',
                        'dir'     : '/store/group/phys_higgs/cmshww/amassiro/RunII/2016/Jul15_DCSONLY/data/25ns/',
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'ICHEP2016' ,
                        'reName'  : {
                                       # Run2016C
                                       'DoubleEG_Run2016C-PromptReco-v2'          : 'Run2016C_PromptReco_DoubleEG',
                                       'DoubleMuon_Run2016C-PromptReco-v2'        : 'Run2016C_PromptReco_DoubleMuon',
                                       'MuonEG_Run2016C-PromptReco-v2'            : 'Run2016C_PromptReco_MuonEG',
                                       'SingleElectron_Run2016C-PromptReco-v2'    : 'Run2016C_PromptReco_SingleElectron',
                                       'SingleMuon_Run2016C-PromptReco-v2'        : 'Run2016C_PromptReco_SingleMuon',
                                    },
                       },
  '15Jul2016_Run2016C_PromptReco_KNU' : {
                        'isData'  : True ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_data_2016_PromptReco.py',
                        'dir'     : 'HWW12fb/15Jul2016_Run2016C_PromptReco/',
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'ICHEP2016' ,
                        'reName'  : {
                                       # Run2016C
                                       'DoubleEG_Run2016C-PromptReco-v2'          : 'Run2016C_PromptReco_DoubleEG',
                                       'DoubleMuon_Run2016C-PromptReco-v2'        : 'Run2016C_PromptReco_DoubleMuon',
                                       'MuonEG_Run2016C-PromptReco-v2'            : 'Run2016C_PromptReco_MuonEG',
                                       'SingleElectron_Run2016C-PromptReco-v2'    : 'Run2016C_PromptReco_SingleElectron',
                                       'SingleMuon_Run2016C-PromptReco-v2'        : 'Run2016C_PromptReco_SingleMuon',
                                    },
                       },

  '15Jul2016_Run2016C_PromptReco_repro' : {
                        'isData'  : True ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_data_2016_PromptReco.py',
                        'dir'     : '/store/group/phys_higgs/cmshww/amassiro/RunII/2016/Jul15_DCSONLY/data/25ns/',
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'ICHEP2016' ,
                        'reName'  : {
                                       # Run2016C
                                       'DoubleEG_Run2016C-PromptReco-v2'          : 'Run2016C_PromptReco_DoubleEG',
                                       'DoubleMuon_Run2016C-PromptReco-v2'        : 'Run2016C_PromptReco_DoubleMuon',
                                       'MuonEG_Run2016C-PromptReco-v2'            : 'Run2016C_PromptReco_MuonEG',
                                       'SingleElectron_Run2016C-PromptReco-v2'    : 'Run2016C_PromptReco_SingleElectron',
                                       'SingleMuon_Run2016C-PromptReco-v2'        : 'Run2016C_PromptReco_SingleMuon',
                                    },
                       },


  '15Jul2016_Run2016D_PromptReco' : {
                        'isData'  : True ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_data_2016_PromptReco.py',
                        'dir'     : '/store/group/phys_higgs/cmshww/amassiro/RunII/2016/Jul15_DCSONLY/data/25ns/',
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'ICHEP2016' ,
                        'reName'  : {
                                       # Run2016D
                                       'DoubleEG_Run2016D-PromptReco-v2'          : 'Run2016D_PromptReco_DoubleEG',
                                       'DoubleMuon_Run2016D-PromptReco-v2'        : 'Run2016D_PromptReco_DoubleMuon',
                                       'MuonEG_Run2016D-PromptReco-v2'            : 'Run2016D_PromptReco_MuonEG',
                                       'SingleElectron_Run2016D-PromptReco-v2'    : 'Run2016D_PromptReco_SingleElectron',
                                       'SingleMuon_Run2016D-PromptReco-v2'        : 'Run2016D_PromptReco_SingleMuon',
                                    },
                       },
  '15Jul2016_Run2016D_PromptReco_KNU' : {
                        'isData'  : True ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_data_2016_PromptReco.py',
                        'dir'     : 'HWW12fb/15Jul2016_Run2016D_PromptReco/',
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'ICHEP2016' ,
                        'reName'  : {
                                       # Run2016D
                                       'DoubleEG_Run2016D-PromptReco-v2'          : 'Run2016D_PromptReco_DoubleEG',
                                       'DoubleMuon_Run2016D-PromptReco-v2'        : 'Run2016D_PromptReco_DoubleMuon',
                                       'MuonEG_Run2016D-PromptReco-v2'            : 'Run2016D_PromptReco_MuonEG',
                                       'SingleElectron_Run2016D-PromptReco-v2'    : 'Run2016D_PromptReco_SingleElectron',
                                       'SingleMuon_Run2016D-PromptReco-v2'        : 'Run2016D_PromptReco_SingleMuon',
                                    },
                       },


  '15Jul2016_Run2016D_PromptReco_repro' : {
                        'isData'  : True ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_data_2016_PromptReco.py',
                        'dir'     : '/store/group/phys_higgs/cmshww/amassiro/RunII/2016/Jul15_DCSONLY/data/25ns/',
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'ICHEP2016' ,
                        'reName'  : {
                                       # Run2016D
                                       'DoubleEG_Run2016D-PromptReco-v2'          : 'Run2016D_PromptReco_DoubleEG',
                                       'DoubleMuon_Run2016D-PromptReco-v2'        : 'Run2016D_PromptReco_DoubleMuon',
                                       'MuonEG_Run2016D-PromptReco-v2'            : 'Run2016D_PromptReco_MuonEG',
                                       'SingleElectron_Run2016D-PromptReco-v2'    : 'Run2016D_PromptReco_SingleElectron',
                                       'SingleMuon_Run2016D-PromptReco-v2'        : 'Run2016D_PromptReco_SingleMuon',
                                    },
                       },


  '26Jul2016_Run2016D_PromptReco' : {
                        'isData'  : True ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_data_2016_PromptReco.py',
                        'dir'     : '/store/group/phys_higgs/cmshww/amassiro/RunII/2016/Jul26/data/25ns/',
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'ICHEP2016' ,
                        'reName'  : {
                                       # Run2016D
                                       'DoubleEG_Run2016D-PromptReco-v2'          : 'Run2016D_PromptReco_DoubleEG',
                                       'DoubleMuon_Run2016D-PromptReco-v2'        : 'Run2016D_PromptReco_DoubleMuon',
                                       'MuonEG_Run2016D-PromptReco-v2'            : 'Run2016D_PromptReco_MuonEG',
                                       'SingleElectron_Run2016D-PromptReco-v2'    : 'Run2016D_PromptReco_SingleElectron',
                                       'SingleMuon_Run2016D-PromptReco-v2'        : 'Run2016D_PromptReco_SingleMuon',
                                    },
                       },
  '26Jul2016_Run2016D_PromptReco_KNU' : {
                        'isData'  : True ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_data_2016_PromptReco.py',
                        'dir'     : 'HWW12fb/26Jul2016_Run2016D_PromptReco/',
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'ICHEP2016' ,
                        'reName'  : {
                                       # Run2016D
                                       'DoubleEG_Run2016D-PromptReco-v2'          : 'Run2016D_PromptReco_DoubleEG',
                                       'DoubleMuon_Run2016D-PromptReco-v2'        : 'Run2016D_PromptReco_DoubleMuon',
                                       'MuonEG_Run2016D-PromptReco-v2'            : 'Run2016D_PromptReco_MuonEG',
                                       'SingleElectron_Run2016D-PromptReco-v2'    : 'Run2016D_PromptReco_SingleElectron',
                                       'SingleMuon_Run2016D-PromptReco-v2'        : 'Run2016D_PromptReco_SingleMuon',
                                    },
                       },


  '26Jul2016_Run2016D_PromptReco_repro' : {
                        'isData'  : True ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_data_2016_PromptReco.py',
                        'dir'     : '/store/group/phys_higgs/cmshww/amassiro/RunII/2016/Jul26/data/25ns/',
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'ICHEP2016' ,
                        'reName'  : {
                                       # Run2016D
                                       'DoubleEG_Run2016D-PromptReco-v2'          : 'Run2016D_PromptReco_DoubleEG',
                                       'DoubleMuon_Run2016D-PromptReco-v2'        : 'Run2016D_PromptReco_DoubleMuon',
                                       'MuonEG_Run2016D-PromptReco-v2'            : 'Run2016D_PromptReco_MuonEG',
                                       'SingleElectron_Run2016D-PromptReco-v2'    : 'Run2016D_PromptReco_SingleElectron',
                                       'SingleMuon_Run2016D-PromptReco-v2'        : 'Run2016D_PromptReco_SingleMuon',
                                    },
                       },


  '31Aug2016_Run2016E_PromptReco' : {
                        'isData'  : True ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_data_2016_PromptReco.py',
                        'dir'     : '/store/group/phys_higgs/cmshww/amassiro/RunII/2016/Aug31/data/25ns/',
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'ICHEP2016' ,
                        'reName'  : {
                                       # Run2016D
                                       'DoubleEG_Run2016E-PromptReco-v2'          : 'Run2016E_PromptReco_DoubleEG',
                                       'DoubleMuon_Run2016E-PromptReco-v2'        : 'Run2016E_PromptReco_DoubleMuon',
                                       'MuonEG_Run2016E-PromptReco-v2'            : 'Run2016E_PromptReco_MuonEG',
                                       'SingleElectron_Run2016E-PromptReco-v2'    : 'Run2016E_PromptReco_SingleElectron',
                                       'SingleMuon_Run2016E-PromptReco-v2'        : 'Run2016E_PromptReco_SingleMuon',
                                    },
                       },

  '31Aug2016_Run2016F_PromptReco' : {
                        'isData'  : True ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_data_2016_PromptReco.py',
                        'dir'     : '/store/group/phys_higgs/cmshww/amassiro/RunII/2016/Aug31/data/25ns/',
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'ICHEP2016' ,
                        'reName'  : {
                                       # Run2016D
                                       'DoubleEG_Run2016F-PromptReco-v1'          : 'Run2016F_PromptReco_DoubleEG',
                                       'DoubleMuon_Run2016F-PromptReco-v1'        : 'Run2016F_PromptReco_DoubleMuon',
                                       'MuonEG_Run2016F-PromptReco-v1'            : 'Run2016F_PromptReco_MuonEG',
                                       'SingleElectron_Run2016F-PromptReco-v1'    : 'Run2016F_PromptReco_SingleElectron',
                                       'SingleMuon_Run2016F-PromptReco-v1'        : 'Run2016F_PromptReco_SingleMuon',
                                    },
                       },

#  2016 ReReco DATA: 29October2016_RerecoData Tag

  'Dec2016_Run2016B_ReReco_27p6fbm1'   : {
                        'isData'  : True ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_data_2016_ReReco.py',
                        'dir'     : '/store/group/phys_higgs/cmshww/amassiro/RunII/2016/23Sep2016/data/25ns/',
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'Rereco2016' ,
                        'onlySample' : [
                                       'DoubleEG_Run2016B-23Sep2016-v3', 'DoubleMuon_Run2016B-23Sep2016-v3', 'MuonEG_Run2016B-23Sep2016-v3', 'SingleElectron_Run2016B-23Sep2016-v3', 'SingleMuon_Run2016B-23Sep2016-v3'
                                       ],
                       },

  'Dec2016_Run2016C_ReReco_27p6fbm1'   : {
                        'isData'  : True ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_data_2016_ReReco.py',
                        'dir'     : '/store/group/phys_higgs/cmshww/amassiro/RunII/2016/23Sep2016/data/25ns/',
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'Rereco2016' ,
                        'onlySample' : [
                                       'DoubleEG_Run2016C-23Sep2016-v3', 'DoubleMuon_Run2016C-23Sep2016-v3', 'MuonEG_Run2016C-23Sep2016-v3', 'SingleElectron_Run2016C-23Sep2016-v3', 'SingleMuon_Run2016C-23Sep2016-v3'
                                       ],
                       },

  'Dec2016_Run2016D_ReReco_27p6fbm1'   : {
                        'isData'  : True ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_data_2016_ReReco.py',
                        'dir'     : '/store/group/phys_higgs/cmshww/amassiro/RunII/2016/23Sep2016/data/25ns/',
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'Rereco2016' ,
                        'onlySample' : [
                                       'DoubleEG_Run2016D-23Sep2016-v3', 'DoubleMuon_Run2016D-23Sep2016-v3', 'MuonEG_Run2016D-23Sep2016-v3', 'SingleElectron_Run2016D-23Sep2016-v3', 'SingleMuon_Run2016D-23Sep2016-v3'
                                       ],
                       },

  'Dec2016_Run2016E_ReReco_27p6fbm1'   : {
                        'isData'  : True ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_data_2016_ReReco.py',
                        'dir'     : '/store/group/phys_higgs/cmshww/amassiro/RunII/2016/23Sep2016/data/25ns/',
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'Rereco2016' ,
                        'onlySample' : [
                                       'DoubleEG_Run2016E-23Sep2016-v3', 'DoubleMuon_Run2016E-23Sep2016-v3', 'MuonEG_Run2016E-23Sep2016-v3', 'SingleElectron_Run2016E-23Sep2016-v3', 'SingleMuon_Run2016E-23Sep2016-v3'
                                       ],
                       },

  'Dec2016_Run2016F_ReReco_27p6fbm1'   : {
                        'isData'  : True ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_data_2016_ReReco.py',
                        'dir'     : '/store/group/phys_higgs/cmshww/amassiro/RunII/2016/23Sep2016/data/25ns/',
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'Rereco2016' ,
                        'onlySample' : [
                                       'DoubleEG_Run2016F-23Sep2016-v3', 'DoubleMuon_Run2016F-23Sep2016-v3', 'MuonEG_Run2016F-23Sep2016-v3', 'SingleElectron_Run2016F-23Sep2016-v3', 'SingleMuon_Run2016F-23Sep2016-v3'
                                       ],
                       },

  'Dec2016_Run2016G_ReReco_27p6fbm1'   : {
                        'isData'  : True ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_data_2016_ReReco.py',
                        'dir'     : '/store/group/phys_higgs/cmshww/amassiro/RunII/2016/23Sep2016/data/25ns/',
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'Rereco2016' ,
                        'onlySample' : [
                                       'DoubleEG_Run2016G-23Sep2016-v3', 'DoubleMuon_Run2016G-23Sep2016-v3', 'MuonEG_Run2016G-23Sep2016-v3', 'SingleElectron_Run2016G-23Sep2016-v3', 'SingleMuon_Run2016G-23Sep2016-v3'
                                       ],
                       },

######## 2016 Rereco DATA : 6JAn2017_RogueOne_v6 

  'Dec2016_Run2016B_ReReco' : {
                        'isData'  : True ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_data_2016_ReReco.py' ,
                        'dir'     : '/store/group/phys_higgs/cmshww/amassiro/RunII/2016/Dec2016/data/25ns/' ,
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'Full2016' ,
                        'onlySample' : [
                                        'DoubleEG_Run2016B-23Sep2016-v3',
                                        'DoubleMuon_Run2016B-23Sep2016-v3',
                                        'MuonEG_Run2016B-23Sep2016-v3',
                                        'SingleElectron_Run2016B-23Sep2016-v3',
                                        'SingleMuon_Run2016B-23Sep2016-v3',
#                                       'MET_Run2016B-23Sep2016-v3',   
                                       ],
                       },

  'Dec2016_Run2016C_ReReco' : {
                        'isData'  : True ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_data_2016_ReReco.py' ,
                        'dir'     : '/store/group/phys_higgs/cmshww/amassiro/RunII/2016/Dec2016/data/25ns/' ,
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'Full2016' ,
                        'onlySample' : [
                                        'DoubleEG_Run2016C-23Sep2016-v1',
                                        'DoubleMuon_Run2016C-23Sep2016-v1',
                                        'MuonEG_Run2016C-23Sep2016-v1',
                                        'SingleElectron_Run2016C-23Sep2016-v1',
                                        'SingleMuon_Run2016C-23Sep2016-v1',
#                                       'MET_Run2016C-23Sep2016-v1',   
                                       ],
                       },

  'Dec2016_Run2016D_ReReco' : {
                        'isData'  : True ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_data_2016_ReReco.py' ,
                        'dir'     : '/store/group/phys_higgs/cmshww/amassiro/RunII/2016/Dec2016/data/25ns/' ,
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'Full2016' ,
                        'onlySample' : [
                                        'DoubleEG_Run2016D-23Sep2016-v1',
                                        'DoubleMuon_Run2016D-23Sep2016-v1',
                                        'MuonEG_Run2016D-23Sep2016-v1',
                                        'SingleElectron_Run2016D-23Sep2016-v1',
                                        'SingleMuon_Run2016D-23Sep2016-v1',
#                                       'MET_Run2016D-23Sep2016-v1',   
                                       ],
                       },

  'Dec2016_Run2016E_ReReco' : {
                        'isData'  : True ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_data_2016_ReReco.py' ,
                        'dir'     : '/store/group/phys_higgs/cmshww/amassiro/RunII/2016/Dec2016/data/25ns/' ,
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'Full2016' ,
                        'onlySample' : [
                                        'DoubleEG_Run2016E-23Sep2016-v1',
                                        'DoubleMuon_Run2016E-23Sep2016-v1',
                                        'MuonEG_Run2016E-23Sep2016-v1',
                                        'SingleElectron_Run2016E-23Sep2016-v1',
                                        'SingleMuon_Run2016E-23Sep2016-v1',
#                                       'MET_Run2016E-23Sep2016-v1',   
                                       ],
                       },

  'Dec2016_Run2016F_ReReco' : {
                        'isData'  : True ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_data_2016_ReReco.py' ,
                        'dir'     : '/store/group/phys_higgs/cmshww/amassiro/RunII/2016/Dec2016/data/25ns/' ,
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'Full2016' ,
                        'onlySample' : [
                                        'DoubleEG_Run2016F-23Sep2016-v1',
                                        'DoubleMuon_Run2016F-23Sep2016-v1',
                                        'MuonEG_Run2016F-23Sep2016-v1',
                                        'SingleElectron_Run2016F-23Sep2016-v1',
                                        'SingleMuon_Run2016F-23Sep2016-v1',
#                                       'MET_Run2016F-23Sep2016-v1',   
                                       ],
                       },

  'Dec2016_Run2016G_ReReco' : {
                        'isData'  : True ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_data_2016_ReReco.py' ,
                        'dir'     : '/store/group/phys_higgs/cmshww/amassiro/RunII/2016/Dec2016/data/25ns/' ,
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'Full2016' ,
                        'onlySample' : [
                                        'DoubleEG_Run2016G-23Sep2016-v1',
                                        'DoubleMuon_Run2016G-23Sep2016-v1',
                                        'MuonEG_Run2016G-23Sep2016-v1',
                                        'SingleElectron_Run2016G-23Sep2016-v1',
                                        'SingleMuon_Run2016G-23Sep2016-v1',
#                                       'MET_Run2016G-23Sep2016-v1',   
                                       ],
                       },

  'Dec2016_Run2016H_PromptReco' : {
                        'isData'  : True ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_data_2016_PromptReco_RunH.py' ,
                        'dir'     : '/store/group/phys_higgs/cmshww/amassiro/RunII/2016/Dec2016/data/25ns/' ,
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'Full2016' ,
                        'onlySample' : [
                                        'DoubleEG_Run2016H-PromptReco-v2',
                                        'DoubleMuon_Run2016H-PromptReco-v2',
                                        'MuonEG_Run2016H-PromptReco-v2',
                                        'SingleElectron_Run2016H-PromptReco-v2',
                                        'SingleMuon_Run2016H-PromptReco-v2',
#                                       'MET_Run2016H-PromptReco-v2',
                                        'DoubleEG_Run2016H-PromptReco-v3',
                                        'DoubleMuon_Run2016H-PromptReco-v3',
                                        'MuonEG_Run2016H-PromptReco-v3',
                                        'SingleElectron_Run2016H-PromptReco-v3',
                                        'SingleMuon_Run2016H-PromptReco-v3',
#                                       'MET_Run2016H-PromptReco-v3',
                                       ],
                       },


#### Summer16 MC: Rogueone_v5  

  'Dec2016_summer16_mAODv2'   : {
                        'isData'  : False ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_summer16.py' ,
                        'dir'     : '/store/group/phys_higgs/cmshww/amassiro/RunII/2017/6Jan_RogueOne/MC/',
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'Full2016' ,
                        # 37.X fb-1
                        'puData'  : '/afs/cern.ch/user/x/xjanssen/public/PileupHistogram_Full2016_271036-284044_69p2mb_31Jan17.root',
                       },


######## 2016 ReMiniAOD DATA : Feb2017_TheGreatWall_v5 

  'Feb2017_Run2016B_RemAOD' : {
                        'isData'  : True ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_data_2016_Re-miniAOD.py' ,
                        'dir'     : '/store/group/phys_higgs/cmshww/amassiro/RunII/2016/Feb2017/data/25ns/' ,
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'Full2016' ,
                        'onlySample' : [
                                        'DoubleEG_Run2016B-03Feb2017_ver2-v2',
                                        'DoubleMuon_Run2016B-03Feb2017_ver2-v2',
                                        'MuonEG_Run2016B-03Feb2017_ver2-v2',
                                        'SingleElectron_Run2016B-03Feb2017_ver2-v2',
                                        'SingleMuon_Run2016B-03Feb2017_ver2-v2',
                                        'MET_Run2016B-03Feb2017_ver2-v2',   
                                       ],
                       },

  'Feb2017_Run2016C_RemAOD' : {
                        'isData'  : True ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_data_2016_Re-miniAOD.py' ,
                        'dir'     : '/store/group/phys_higgs/cmshww/amassiro/RunII/2016/Feb2017/data/25ns/' ,
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'Full2016' ,
                        'onlySample' : [
                                        'DoubleEG_Run2016C-03Feb2017-v1',
                                        'DoubleMuon_Run2016C-03Feb2017-v1',
                                        'MuonEG_Run2016C-03Feb2017-v1',
                                        'SingleElectron_Run2016C-03Feb2017-v1',
                                        'SingleMuon_Run2016C-03Feb2017-v1',
                                        'MET_Run2016C-03Feb2017-v1',   
                                       ],
                       },

  'Feb2017_Run2016D_RemAOD' : {
                        'isData'  : True ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_data_2016_Re-miniAOD.py' ,
                        'dir'     : '/store/group/phys_higgs/cmshww/amassiro/RunII/2016/Feb2017/data/25ns/' ,
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'Full2016' ,
                        'onlySample' : [
                                        'DoubleEG_Run2016D-03Feb2017-v1',
                                        'DoubleMuon_Run2016D-03Feb2017-v1',
                                        'MuonEG_Run2016D-03Feb2017-v1',
                                        'SingleElectron_Run2016D-03Feb2017-v1',
                                        'SingleMuon_Run2016D-03Feb2017-v1',
                                        'MET_Run2016D-03Feb2017-v1',   
                                       ],
                       },

  'Feb2017_Run2016E_RemAOD' : {
                        'isData'  : True ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_data_2016_Re-miniAOD.py' ,
                        'dir'     : '/store/group/phys_higgs/cmshww/amassiro/RunII/2016/Feb2017/data/25ns/' ,
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'Full2016' ,
                        'onlySample' : [
                                        'DoubleEG_Run2016E-03Feb2017-v1',
                                        'DoubleMuon_Run2016E-03Feb2017-v1',
                                        'MuonEG_Run2016E-03Feb2017-v1',
                                        'SingleElectron_Run2016E-03Feb2017-v1',
                                        'SingleMuon_Run2016E-03Feb2017-v1',
                                        'MET_Run2016E-03Feb2017-v1',   
                                       ],
                       },

  'Feb2017_Run2016F_RemAOD' : {
                        'isData'  : True ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_data_2016_Re-miniAOD.py' ,
                        'dir'     : '/store/group/phys_higgs/cmshww/amassiro/RunII/2016/Feb2017/data/25ns/' ,
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'Full2016' ,
                        'onlySample' : [
                                        'DoubleEG_Run2016F-03Feb2017-v1',
                                        'DoubleMuon_Run2016F-03Feb2017-v1',
                                        'MuonEG_Run2016F-03Feb2017-v1',
                                        'SingleElectron_Run2016F-03Feb2017-v1',
                                        'SingleMuon_Run2016F-03Feb2017-v1',
                                        'MET_Run2016F-03Feb2017-v1',   
                                       ],
                       },

  'Feb2017_Run2016G_RemAOD' : {
                        'isData'  : True ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_data_2016_Re-miniAOD.py' ,
                        'dir'     : '/store/group/phys_higgs/cmshww/amassiro/RunII/2016/Feb2017/data/25ns/' ,
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'Full2016' ,
                        'onlySample' : [
                                        'DoubleEG_Run2016G-03Feb2017-v1',
                                        'DoubleMuon_Run2016G-03Feb2017-v1',
                                        'MuonEG_Run2016G-03Feb2017-v1',
                                        'SingleElectron_Run2016G-03Feb2017-v1',
                                        'SingleMuon_Run2016G-03Feb2017-v1',
                                        'MET_Run2016G-03Feb2017-v1',   
                                       ],
                       },

  'Feb2017_Run2016G_RemAOD_Dec2016Fix' : {
                        'isData'  : True ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_data_2016_Re-miniAOD.py' ,
                        'dir'     : '/store/group/phys_higgs/cmshww/amassiro/RunII/2016/Feb2017/data/25ns/' ,
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'Full2016' ,
                        'onlySample' : [
                                        'SingleMuon_Run2016G-03Feb2017-v1',
                                       ],
                       },


  'Feb2017_Run2016H_RemAOD' : {
                        'isData'  : True ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_data_2016_Re-miniAOD_RunH.py' ,
                        'dir'     : '/store/group/phys_higgs/cmshww/amassiro/RunII/2016/Feb2017/data/25ns/' ,
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'Full2016' ,
                        'onlySample' : [
                                        'DoubleEG_Run2016H-03Feb2017_ver2-v1',
                                        'DoubleMuon_Run2016H-03Feb2017_ver2-v1',
                                        'MuonEG_Run2016H-03Feb2017_ver2-v1',
                                        'SingleElectron_Run2016H-03Feb2017_ver2-v1',
                                        'SingleMuon_Run2016H-03Feb2017_ver2-v1',
                                        'MET_Run2016H-03Feb2017_ver2-v1',
                                        'DoubleEG_Run2016H-03Feb2017_ver3-v1',
                                        'DoubleMuon_Run2016H-03Feb2017_ver3-v1',
                                        'MuonEG_Run2016H-03Feb2017_ver3-v1',
                                        'SingleElectron_Run2016H-03Feb2017_ver3-v1',
                                        'SingleMuon_Run2016H-03Feb2017_ver3-v1',
                                        'MET_Run2016H-03Feb2017_ver3-v1',
                                       ],
                       },


#### Summer16 MC: Feb2017_TheGreatWall_v5  

  'Feb2017_summer16'   : {
                        'isData'  : False ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_summer16.py' ,
                        'dir'     : '/store/group/phys_higgs/cmshww/amassiro/RunII/2016/Feb2017/MC/',
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'Full2016' ,
                        # 37.X fb-1
                        'puData'  : '/afs/cern.ch/user/x/xjanssen/public/PileupHistogram_Full2016_271036-284044_69p2mb_31Jan17.root',
                       },


#### Stop Preliminary Spring16 MC:   

  'Feb2017_spring16_stop'   : {
                        'isData'  : False ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_spring16_miniaod_v2_Feb2017.py' ,
                        'dir'     : '/gpfs/csic_projects/cms/sluca/SUSY/Stop/Moriond17/LatinoTrees/MC/',
                        'dirExt'  : 'Spring16' ,
                        'cmssw'   : 'Full2016' ,
                        # 37.X fb-1
                        'puData'  : '/gpfs/csic_users/sluca/Stop/CMSSW_8_0_17/src/LatinoAnalysis/Gardener/scripts/PileupHistogram_8020_272007-275376_69.2mb_26Oct.root',
                        'LepKind' : 'Ghent',
                       },  

  'Feb2017_spring16_stop_pog'   : {
                        'isData'  : False ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_spring16_miniaod_v2_Feb2017.py' ,
                        'dir'     : '/gpfs/csic_projects/cms/sluca/SUSY/Stop/Moriond17/LatinoTrees/MC/',
                        'dirExt'  : 'Spring16' ,
                        'cmssw'   : 'Full2016' ,
                        # 37.X fb-1
                        'puData'  : '/gpfs/csic_users/sluca/Stop/CMSSW_8_0_17/src/LatinoAnalysis/Gardener/scripts/PileupHistogram_8020_272007-275376_69.2mb_26Oct.root',
                        'LepKind' : 'POG',
                       }, 

  'Feb2017_spring16_stop_ichep'   : {
                        'isData'  : False ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_spring16_miniaod_v2_Feb2017_WW.py' ,
                        'dir'     : '/gpfs/csic_projects/cms/sluca/SUSY/Stop/Moriond17/LatinoTrees/MC/',
                        'dirExt'  : 'Spring16' ,
                        'cmssw'   : 'ICHEP2016' ,
                        # 37.X fb-1
                        'puData'  : '/gpfs/csic_users/sluca/Stop/CMSSW_8_0_17/src/LatinoAnalysis/Gardener/scripts/PileupHistogram_8020_272007-275376_69.2mb_26Oct.root',
                        'LepKind' : 'None',
                       },

#### Stop Preliminary ReReco

  'Feb2017_Run2016B_ReReco' : {
        'isData'  : True ,
        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_data_2016_ReReco.py',
        'dir'     : '/gpfs/csic_projects/cms/sluca/SUSY/Stop/Moriond17/LatinoTrees/Data/',
        'dirExt'  : 'ReReco' ,
        'cmssw'   : 'Full2016' ,
        'reName'  : {
            # Run2016B
            'DoubleEG_Run2016B-23Sep2016-v3'          : 'Run2016B_23Sep2016_DoubleEG',
            'DoubleMuon_Run2016B-23Sep2016-v3'        : 'Run2016B_23Sep2016_DoubleMuon',
            'MuonEG_Run2016B-23Sep2016-v3'            : 'Run2016B_23Sep2016_MuonEG',
            'SingleElectron_Run2016B-23Sep2016-v3'    : 'Run2016B_23Sep2016_SingleElectron',
            'SingleMuon_Run2016B-23Sep2016-v3'        : 'Run2016B_23Sep2016_SingleMuon',
            },
            'LepKind' : 'Ghent',
        },

  'Feb2017_Run2016B_ReReco_pog' : {
        'isData'  : True ,
        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_data_2016_ReReco.py',
        'dir'     : '/gpfs/csic_projects/cms/sluca/SUSY/Stop/Moriond17/LatinoTrees/Data/',
        'dirExt'  : 'ReReco' ,
        'cmssw'   : 'Full2016' ,
        'reName'  : {
            # Run2016B
            'DoubleEG_Run2016B-23Sep2016-v3'          : 'Run2016B_23Sep2016_DoubleEG',
            'DoubleMuon_Run2016B-23Sep2016-v3'        : 'Run2016B_23Sep2016_DoubleMuon',
            'MuonEG_Run2016B-23Sep2016-v3'            : 'Run2016B_23Sep2016_MuonEG',
            'SingleElectron_Run2016B-23Sep2016-v3'    : 'Run2016B_23Sep2016_SingleElectron',
            'SingleMuon_Run2016B-23Sep2016-v3'        : 'Run2016B_23Sep2016_SingleMuon',
            },
            'LepKind' : 'POG',
        },


#### Stop Preliminary Summer16 MC:   

  'Feb2017_spring16_stop_ghent'   : {
                        'isData'  : False ,
                        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_summer16.py' ,
                        'dir'     : '/gpfs/gaes/cms/store/group/phys_higgs/cmshww/amassiro/RunII/2016/Feb2017/MC/',
                        'dirExt'  : 'LatinoTrees' ,
                        'cmssw'   : 'Full2016' ,
                        # 37.X fb-1
                        'puData'  : '/gpfs/csic_users/sluca/Stop/CMSSW_8_0_17/src/LatinoAnalysis/Gardener/scripts/PileupHistogram_8020_272007-275376_69.2mb_26Oct.root',
                        'LepKind' : 'Ghent',
                       },  

#### Stop ReMiniAOD

  'Feb2017_Run2016B_ReMiniAOD' : {
        'isData'  : True ,
        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_data_2016_Re-miniAOD.py',
        'dir'     : '/gpfs/gaes/cms/store/group/phys_higgs/cmshww/amassiro/RunII/2016/Feb2017/data/25ns/LatinoTrees/',
        'dirExt'  : 'Run2016B' ,
        'cmssw'   : 'Full2016' ,
        'onlySample' : [
            'DoubleEG_Run2016B-03Feb2017_ver2-v2',
            'DoubleMuon_Run2016B-03Feb2017_ver2-v2',
            'MuonEG_Run2016B-03Feb2017_ver2-v2',
            'SingleElectron_Run2016B-03Feb2017_ver2-v2',
            'SingleMuon_Run2016B-03Feb2017_ver2-v2',
            'MET_Run2016B-03Feb2017_ver2-v2',   
            ],
            'LepKind' : 'Ghent',
        },

  'Feb2017_Run2016C_ReMiniAOD' : {
        'isData'  : True ,
        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_data_2016_Re-miniAOD.py',
        'dir'     : '/gpfs/gaes/cms/store/group/phys_higgs/cmshww/amassiro/RunII/2016/Feb2017/data/25ns/LatinoTrees/',
        'dirExt'  : 'Run2016C' ,
        'cmssw'   : 'Full2016' ,
        'onlySample' : [
            'DoubleEG_Run2016C-03Feb2017-v1',
            'DoubleMuon_Run2016C-03Feb2017-v1',
            'MuonEG_Run2016C-03Feb2017-v1',
            'SingleElectron_Run2016C-03Feb2017-v1',
            'SingleMuon_Run2016C-03Feb2017-v1',
            'MET_Run2016C-03Feb2017-v1',   
            ],
            'LepKind' : 'Ghent',
        },

  'Feb2017_Run2016D_ReMiniAOD' : {
        'isData'  : True ,
        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_data_2016_Re-miniAOD.py',
        'dir'     : '/gpfs/gaes/cms/store/group/phys_higgs/cmshww/amassiro/RunII/2016/Feb2017/data/25ns/LatinoTrees/',
        'dirExt'  : 'Run2016D' ,
        'cmssw'   : 'Full2016' ,
        'onlySample' : [
            'DoubleEG_Run2016D-03Feb2017-v1',
            'DoubleMuon_Run2016D-03Feb2017-v1',
            'MuonEG_Run2016D-03Feb2017-v1',
            'SingleElectron_Run2016D-03Feb2017-v1',
            'SingleMuon_Run2016D-03Feb2017-v1',
            'MET_Run2016D-03Feb2017-v1',   
            ],
            'LepKind' : 'Ghent',
        },

  'Feb2017_Run2016E_ReMiniAOD' : {
        'isData'  : True ,
        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_data_2016_Re-miniAOD.py',
        'dir'     : '/gpfs/gaes/cms/store/group/phys_higgs/cmshww/amassiro/RunII/2016/Feb2017/data/25ns/LatinoTrees/',
        'dirExt'  : 'Run2016E' ,
        'cmssw'   : 'Full2016' ,
        'onlySample' : [
            'DoubleEG_Run2016E-03Feb2017-v1',
            'DoubleMuon_Run2016E-03Feb2017-v1',
            'MuonEG_Run2016E-03Feb2017-v1',
            'SingleElectron_Run2016E-03Feb2017-v1',
            'SingleMuon_Run2016E-03Feb2017-v1',
            'MET_Run2016E-03Feb2017-v1',   
            ],
            'LepKind' : 'Ghent',
        },

  'Feb2017_Run2016B_ReMiniAOD_pog' : {
        'isData'  : True ,
        'samples' : 'LatinoTrees/AnalysisStep/test/crab/samples/samples_data_2016_Re-miniAOD.py',
        'dir'     : '/gpfs/gaes/cms/store/group/phys_higgs/cmshww/amassiro/RunII/2016/Feb2017/data/25ns/',
        'dirExt'  : 'LatinoTrees' ,
        'cmssw'   : 'Full2016' ,
        'onlySample' : [
            'DoubleEG_Run2016B-03Feb2017_ver2-v2',
            'DoubleMuon_Run2016B-03Feb2017_ver2-v2',
            'MuonEG_Run2016B-03Feb2017_ver2-v2',
            'SingleElectron_Run2016B-03Feb2017_ver2-v2',
            'SingleMuon_Run2016B-03Feb2017_ver2-v2',
            'MET_Run2016B-03Feb2017_ver2-v2',   
            ],
            'LepKind' : 'POG',
        },


}



# ---- Steps
# .... .... this is defined by mkGardener in "-s" "--steps" option
# .... .... if it is a "chain", it means that the intermediate steps are NOT saved
# .... ....    e.g. 'puadder','baseW','wwNLL' ---> only after all steps the folder will be saved on eos


samples4Syst = [
                 # DY 
                 'DYJetsToLL_M-10to50','DYJetsToLL_M-50','DYJetsToLL_M-10to50ext3',
                  # ... ICHEP16
                 'DY2JetsToLL', 'DY3JetsToLL', 'DY4JetsToLL',
                 'DYJetsToLL_M-50_HT-100to200_MLM' , 'DYJetsToLL_M-50_HT-100to200_MLM_ext' , 
                 'DYJetsToLL_M-50_HT-200to400_MLM' , 
                 'DYJetsToLL_M-50_HT-400to600_MLM' ,
                 'DYJetsToLL_M-50_HT-600toInf_MLM' ,
                  # ... Moriond17
                 'DYJetsToLL_M-5to50_HT-70to100'  , 'DYJetsToLL_M-5to50_HT-100to200' , 'DYJetsToLL_M-5to50_HT-200to400' ,
                 'DYJetsToLL_M-5to50_HT-400to600' , 'DYJetsToLL_M-5to50_HT-600toInf' ,
                 'DYJetsToLL_M-50_HT-70to100'  , 'DYJetsToLL_M-50_HT-100to200' , 'DYJetsToLL_M-50_HT-100to200_ext1' ,
                 'DYJetsToLL_M-50_HT-200to400' , 'DYJetsToLL_M-50_HT-200to400_ext1' , 'DYJetsToLL_M-50_HT-400to600' ,
                 'DYJetsToLL_M-50_HT-600to800' , 'DYJetsToLL_M-50_HT-800to1200' , 'DYJetsToLL_M-50_HT-1200to2500' ,
                 'DYJetsToLL_M-50_HT-2500toInf', 
                 # ... DY -> MuEle
                 'DYJetsToTT_MuEle_M-50' , 'DYJetsToTT_MuEle_M-50_ext1' ,   

                 # Top
                 'TTTo2L2Nu','TT','TTTo2L2Nu_ext1',
                 'ST_s-channel' , 
                 'ST_t-channel_antitop','ST_t-channel_top',
                 'ST_tW_antitop','ST_tW_top',
                 'ST_tW_antitop_noHad' , 'ST_tW_antitop_noHad_ext1' , 'ST_tW_top_noHad' , 'ST_tW_top_noHad_ext1' ,

                 # VV (including WW) 
                 'WWTo2L2Nu','GluGluWWTo2L2Nu_MCFM','GluGluWWTo2L2NuHiggs_MCFM',
                 'WZTo3LNu','WZ','WZTo2L2Q','WZTo3LNu_mllmin01',
                 'ZZ','ZZTo2L2Nu','ZZTo2L2Q','ZZTo4L',
                 'Wg_AMCNLOFXFX','WgStarLNuEE','WgStarLNuMuMu', 'Wg_MADGRAPHMLM', 'Zg' ,

                 # VVV
                 'WZZ','ZZZ','WWZ','WWW','WWG',

                 # Higgs 
                 'GluGluHToTauTau_M125' , 'GluGluHToWWTo2L2Nu_M125','GluGluHToWWTo2L2NuPowheg_M125',
                 'GluGluHToWWTo2L2Nu_alternative_M125',
                 'VBFHToTauTau_M125' , 'VBFHToWWTo2L2Nu_M125','VBFHToWWTo2L2Nu_alternative_M125',
                 'HWminusJ_HToTauTau_M125', 'HWminusJ_HToWW_M125', 'HWminusJ_HToWW_LNu_M125' ,
                 'HWplusJ_HToTauTau_M125' , 'HWplusJ_HToWW_M125',  'HWplusJ_HToWW_LNu_M125' ,
                 'HZJ_HToTauTau_M125'     , 'HZJ_HToWW_M125', 'HZJ_HToWWTo2L2Nu_M125',
                 'bbHToWWTo2L2Nu_M125_yb2', 'bbHToWWTo2L2Nu_M125_ybyt',

                 'ggZH_HToWW_M125', # missing ggZHToTauTau

                 # What ????
                 #'ttHJetToNonbb_M125','TTWJetsToLNu',
                 #'GluGluHToZZTo4L_M125',
                 #'GluGluZH_HToWWTo2L2Nu_M120_noHLT','GluGluZH_HToWWTo2L2Nu_M125_noHLT','GluGluZH_HToWWTo2L2Nu_M130_noHLT',
                 #'HZJ_HToWWTo2L2Nu_M120_noHLT','HZJ_HToWWTo2L2Nu_M125_noHLT','HZJ_HToWWTo2L2Nu_M130_noHLT', 
                 #'HWplusJ_WToLNu_HToWWTo2L2Nu_M125','HWminusJ_WToLNu_HToWWTo2L2Nu_M125',

                 # ... ggH High Mass
                 'GluGluHToWWTo2L2Nu_M130',
                 'GluGluHToWWTo2L2Nu_M135',
                 'GluGluHToWWTo2L2Nu_M140',
                 'GluGluHToWWTo2L2Nu_M145',
                 'GluGluHToWWTo2L2Nu_M150',
                 'GluGluHToWWTo2L2Nu_M155',
                 'GluGluHToWWTo2L2Nu_M160',
                 'GluGluHToWWTo2L2Nu_M165',
                 'GluGluHToWWTo2L2Nu_M170',
                 'GluGluHToWWTo2L2Nu_M175',
                 'GluGluHToWWTo2L2Nu_M180',
                 'GluGluHToWWTo2L2Nu_M190',
                 'GluGluHToWWTo2L2Nu_M200',
                 'GluGluHToWWTo2L2Nu_M210',
                 'GluGluHToWWTo2L2Nu_M230',
                 'GluGluHToWWTo2L2Nu_M250',
                 'GluGluHToWWTo2L2Nu_M270',
                 'GluGluHToWWTo2L2Nu_M300', 'GluGluHToWWTo2L2Nu_JHUGen698_M300' ,
                 'GluGluHToWWTo2L2Nu_M350', 'GluGluHToWWTo2L2Nu_JHUGen698_M350' ,
                 'GluGluHToWWTo2L2Nu_M400', 'GluGluHToWWTo2L2Nu_JHUGen698_M400' ,
                 'GluGluHToWWTo2L2Nu_M450', 'GluGluHToWWTo2L2Nu_JHUGen698_M450' ,
                 'GluGluHToWWTo2L2Nu_M500', 'GluGluHToWWTo2L2Nu_JHUGen698_M500' ,
                 'GluGluHToWWTo2L2Nu_M550', 'GluGluHToWWTo2L2Nu_JHUGen698_M550' ,
                 'GluGluHToWWTo2L2Nu_M600', 'GluGluHToWWTo2L2Nu_JHUGen698_M600' ,
                 'GluGluHToWWTo2L2Nu_M650', 'GluGluHToWWTo2L2Nu_JHUGen698_M650' ,
                 'GluGluHToWWTo2L2Nu_M700', 'GluGluHToWWTo2L2Nu_JHUGen698_M700' ,
                 'GluGluHToWWTo2L2Nu_M750', 'GluGluHToWWTo2L2Nu_JHUGen698_M750' , 'GluGluHToWWTo2L2Nu_M750_NWA',
                 'GluGluHToWWTo2L2Nu_M800', 'GluGluHToWWTo2L2Nu_JHUGen698_M800' ,
                 'GluGluHToWWTo2L2Nu_M900', 'GluGluHToWWTo2L2Nu_JHUGen698_M900' ,
                 'GluGluHToWWTo2L2Nu_M1000','GluGluHToWWTo2L2Nu_JHUGen698_M1000',
                 'GluGluHToWWTo2L2Nu_JHUGen698_M1500' ,
                 'GluGluHToWWTo2L2Nu_JHUGen698_M2000' ,
                 'GluGluHToWWTo2L2Nu_JHUGen698_M2500' ,
                 'GluGluHToWWTo2L2Nu_JHUGen698_M3000' ,
                 # ... VBF High Mass
                 'VBFHToWWTo2L2Nu_M130',
                 'VBFHToWWTo2L2Nu_M135',
                 'VBFHToWWTo2L2Nu_M140',
                 'VBFHToWWTo2L2Nu_M145',
                 'VBFHToWWTo2L2Nu_M150',
                 'VBFHToWWTo2L2Nu_M155',
                 'VBFHToWWTo2L2Nu_M160',
                 'VBFHToWWTo2L2Nu_M165',
                 'VBFHToWWTo2L2Nu_M170',
                 'VBFHToWWTo2L2Nu_M175',
                 'VBFHToWWTo2L2Nu_M180',
                 'VBFHToWWTo2L2Nu_M190',
                 'VBFHToWWTo2L2Nu_M200',
                 'VBFHToWWTo2L2Nu_M210',
                 'VBFHToWWTo2L2Nu_M230',
                 'VBFHToWWTo2L2Nu_M250',
                 'VBFHToWWTo2L2Nu_M270',
                 'VBFHToWWTo2L2Nu_M300', 'VBFHToWWTo2L2Nu_JHUGen698_M300' ,
                 'VBFHToWWTo2L2Nu_M350', 'VBFHToWWTo2L2Nu_JHUGen698_M350' ,
                 'VBFHToWWTo2L2Nu_M400', 'VBFHToWWTo2L2Nu_JHUGen698_M400' ,
                 'VBFHToWWTo2L2Nu_M450', 'VBFHToWWTo2L2Nu_JHUGen698_M450' ,
                 'VBFHToWWTo2L2Nu_M500', 'VBFHToWWTo2L2Nu_JHUGen698_M500' ,
                 'VBFHToWWTo2L2Nu_M550', 'VBFHToWWTo2L2Nu_JHUGen698_M550' ,
                 'VBFHToWWTo2L2Nu_M600', 'VBFHToWWTo2L2Nu_JHUGen698_M600' ,
                 'VBFHToWWTo2L2Nu_M650', 'VBFHToWWTo2L2Nu_JHUGen698_M650' ,
                 'VBFHToWWTo2L2Nu_M700', 'VBFHToWWTo2L2Nu_JHUGen698_M700' ,
                 'VBFHToWWTo2L2Nu_M750', 'VBFHToWWTo2L2Nu_JHUGen698_M750' , 'VBFHToWWTo2L2Nu_M750_NWA',
                 'VBFHToWWTo2L2Nu_M800', 'VBFHToWWTo2L2Nu_JHUGen698_M800' ,
                 'VBFHToWWTo2L2Nu_M900', 'VBFHToWWTo2L2Nu_JHUGen698_M900' ,
                 'VBFHToWWTo2L2Nu_M1000','VBFHToWWTo2L2Nu_JHUGen698_M1000' ,
                 'VBFHToWWTo2L2Nu_JHUGen698_M1500' ,        
                 'VBFHToWWTo2L2Nu_JHUGen698_M2000' ,        
                 'VBFHToWWTo2L2Nu_JHUGen698_M2500' ,        
                 'VBFHToWWTo2L2Nu_JHUGen698_M3000' ,        

                 # VBF 
                 'WpWmJJ_EWK_QCD_noTop','WpWmJJ_QCD_noTop', 'WpWmJJ_EWK_noTop',
                 # VBS
                 'WmWmJJ_EWK_powheg','WpWpJJ_EWK_powheg','WpWpJJ_EWK_aQGC','WpWpJJ_EWK',
                 'WpWpJJ_EWK_QCD','WpWpJJ_QCD','WW_DoubleScattering','WWTo2L2Nu_DoubleScattering','WGJJ','WLLJJToLNu_M-60_EWK_4F',
                 'WLLJJToLNu_M-50_QCD_0Jet','WLLJJToLNu_M-50_QCD_1Jet','WLLJJToLNu_M-50_QCD_2Jet','WLLJJToLNu_M-50_QCD_3Jet',
                 'WLLJJToLNu_M-4To60_EWK_4F','WLLJJToLNu_M-4To50_QCD_0Jet','WLLJJToLNu_M-4To50_QCD_1Jet',
                 'WLLJJToLNu_M-4To50_QCD_2Jet','WLLJJToLNu_M-4To50_QCD_3Jet','tZq_ll','ZZJJTo4L_EWK', 'ZZTo4L',                                   
                 'WZTo3LNu','Wg_AMCNLOFXFX','Wg_MADGRAPHMLM','TTTo2L2Nu','ggZZ4e','ggZZ4m','ggZZ4t','ggZZ2e2m','ggZZ2e2t','ggZZ2m2t',


                 # ttDM
                 'ttDM0001pseudo00010',
                 'ttDM0001pseudo00020',
                 'ttDM0001pseudo00050',
                 'ttDM0001pseudo00100',
                 'ttDM0001pseudo00200',
                 'ttDM0001pseudo00300',
                 'ttDM0001pseudo00500',
                 'ttDM0001scalar00010',
                 'ttDM0001scalar00020',
                 'ttDM0001scalar00050',
                 'ttDM0001scalar00100',
                 'ttDM0001scalar00200',
                 'ttDM0001scalar00300',
                 'ttDM0001scalar00500',

                 'monoH_2HDM_MZp-1000_MA0-300',
                 'monoH_2HDM_MZp-1000_MA0-400',
                 'monoH_2HDM_MZp-1000_MA0-500',
                 'monoH_2HDM_MZp-1000_MA0-600',
                 'monoH_2HDM_MZp-1000_MA0-700',
                 'monoH_2HDM_MZp-1000_MA0-800',
                 'monoH_2HDM_MZp-1200_MA0-400',
                 'monoH_2HDM_MZp-1200_MA0-600',
                 'monoH_2HDM_MZp-1200_MA0-700',
                 'monoH_2HDM_MZp-1200_MA0-800',
                 'monoH_2HDM_MZp-1400_MA0-300',
                 'monoH_2HDM_MZp-1400_MA0-400',
                 'monoH_2HDM_MZp-1400_MA0-500',
                 'monoH_2HDM_MZp-1400_MA0-600',
                 'monoH_2HDM_MZp-1400_MA0-700',
                 'monoH_2HDM_MZp-1400_MA0-800',
                 'monoH_2HDM_MZp-1700_MA0-300',
                 'monoH_2HDM_MZp-1700_MA0-600',
                 'monoH_2HDM_MZp-1700_MA0-700',
                 'monoH_2HDM_MZp-1700_MA0-800',
                 'monoH_2HDM_MZp-2000_MA0-300',
                 'monoH_2HDM_MZp-2000_MA0-600',
                 'monoH_2HDM_MZp-2000_MA0-700',
                 'monoH_2HDM_MZp-2000_MA0-800',
                 'monoH_2HDM_MZp-2500_MA0-400',
                 'monoH_2HDM_MZp-600_MA0-400',
                 'monoH_2HDM_MZp-600_MA0-300',
                 'monoH_2HDM_MZp-800_MA0-300',
                 'monoH_2HDM_MZp-800_MA0-400',
                 'monoH_ZpBaryonic_MZp-10000_MChi-50',
                 'monoH_ZpBaryonic_MZp-10000_MChi-500',
                 'monoH_ZpBaryonic_MZp-1000_MChi-1',
                 'monoH_ZpBaryonic_MZp-1000_MChi-1000',
                 'monoH_ZpBaryonic_MZp-1000_MChi-150',
                 'monoH_ZpBaryonic_MZp-100_MChi-1',
                 'monoH_ZpBaryonic_MZp-100_MChi-10',
                 'monoH_ZpBaryonic_MZp-10_MChi-1',
                 'monoH_ZpBaryonic_MZp-10_MChi-1000',
                 'monoH_ZpBaryonic_MZp-10_MChi-50',
                 'monoH_ZpBaryonic_MZp-10_MChi-500',
                 'monoH_ZpBaryonic_MZp-15_MChi-10',
                 'monoH_ZpBaryonic_MZp-200_MChi-150',
                 'monoH_ZpBaryonic_MZp-300_MChi-1',
                 'monoH_ZpBaryonic_MZp-300_MChi-50',
                 'monoH_ZpBaryonic_MZp-500_MChi-150',
                 'monoH_ZpBaryonic_MZp-500_MChi-500',
                 'monoH_ZpBaryonic_MZp-50_MChi-1',
                 'monoH_ZpBaryonic_MZp-50_MChi-10',
                 'monoH_ZpBaryonic_MZp-50_MChi-50',
                 'monoH_ZpBaryonic_MZp-995_MChi-500',

                 'monoH_2HDM_MZp-2500_MA0-700',
                 'monoH_2HDM_MZp-2500_MA0-600',
                 'monoH_2HDM_MZp-1200_MA0-300',
                 'monoH_2HDM_MZp-2000_MA0-400',
                 'monoH_2HDM_MZp-1200_MA0-500',
                 'monoH_2HDM_MZp-1700_MA0-400',
                 'monoH_2HDM_MZp-800_MA0-600',
                 'monoH_2HDM_MZp-2500_MA0-800',
                 'monoH_2HDM_MZp-800_MA0-500',
                 'monoH_2HDM_MZp-1700_MA0-500',
                 'monoH_2HDM_MZp-2500_MA0-500',

                ] 


Steps= {

# ... Chains

  'MC' :       {
                  'isChain'    : True ,
                  'do4MC'      : True , 
                  'do4Data'    : False,
                  'subTargets' : ['puadder','baseW','wwNLL']
                },

  'MCl2sel' :       {
                  'isChain'    : True ,
                  'do4MC'      : True ,
                  'do4Data'    : False,
                  'subTargets' : ['puadder','baseW','wwNLL','l2sel','genVariables']
                },

  'MCl2loose' :       {
                  'isChain'    : True ,
                  'do4MC'      : True ,
                  'do4Data'    : False,
                  'subTargets' : ['do_l2loose','puadder','baseW','wwNLL','genVariables','genMatchVariables','l2kin','l3kin','l4kin','BWEwkSinglet','BWEwkSinglet_JHUGen698','TopGenPt'],
                },


  'MCl2stop' :       {
                  'isChain'    : True ,
                  'do4MC'      : True ,
                  'do4Data'    : False,
                  'subTargets' : ['do_l2stop','puadder','baseW','wwNLL','genVariables','genMatchVariables','l2kin','l3kin','l4kin','BWEwkSinglet','BWEwkSinglet_JHUGen698','TopGenPt'],
                },

  'MCl2looseCut' :       {
                  'isChain'    : True ,
                  'do4MC'      : True ,
                  'do4Data'    : False,
                  'subTargets' : ['do_l2loose_Cut','puadder','baseW','wwNLL','genVariables','genMatchVariables','l2kin','l3kin','l4kin','do_dymvaHiggs','BWEwkSinglet','wwEWK','wzEWK','zzEWK'] #,'BWEwkSinglet_JHUGen698','TopGenPt'],
                },


  'l2loose'  :       {
                  'isChain'    : True ,
                  'do4MC'      : False ,
                  'do4Data'    : True,
                  'subTargets' : ['do_l2loose','l2kin','l3kin','l4kin'],
                },

  'l2stop'  :       {
                  'isChain'    : True ,
                  'do4MC'      : False ,
                  'do4Data'    : True,
                  'subTargets' : ['do_l2stop','l2kin','l3kin','l4kin'],
                },

  'l2vloose'  :       {
                  'isChain'    : True ,
                  'do4MC'      : False ,
                  'do4Data'    : True,
                  'subTargets' : ['do_l2vloose','l2kin','l3kin','l4kin'],
                },

  'l2looseCut'  :       {
                  'isChain'    : True ,
                  'do4MC'      : False ,
                  'do4Data'    : True,
                  'subTargets' : ['do_l2loose_Cut','l2kin','l3kin','l4kin','do_dymvaHiggs'],
                },

  'l2vlooseCut'  :       {
                  'isChain'    : True ,
                  'do4MC'      : False ,
                  'do4Data'    : True,
                  'subTargets' : ['do_l2vloose_Cut','l2kin','l3kin','l4kin','do_dymvaHiggs'],
                },


# 'l1loose16'  :       {
#                 'isChain'    : True ,
#                 'do4MC'      : False ,
#                 'do4Data'    : True,
#                 'subTargets' : ['l1loose','l2kin'],
#               },



  'MCl2vloose' :       {
                  'isChain'    : True ,
                  'do4MC'      : True ,
                  'do4Data'    : False,
                  'subTargets' : ['do_l2vloose','puadder','baseW','wwNLL','genVariables','genMatchVariables','l2kin','l3kin','l4kin','BWEwkSinglet','TopGenPt'],
		  'onlySample' : [ 
			   	  # VBS
				  'DYJetsToLL_M-50','DYJetsToLL_M-10to50ext3','DYJetsToLL_M-10to50',
                                  'WmWmJJ_EWK_powheg',
 			  	  'WpWpJJ_EWK','WpWpJJ_EWK_QCD','WpWpJJ_QCD','WW_DoubleScattering','WWTo2L2Nu_DoubleScattering','WLLJJToLNu_M-4to60_EWK_QCD','WLLJJToLNu_M-60_EWK_QCD','WGJJ',
                                  'TTToSemiLepton','TTToSemiLeptonic','DY2JetsToLL','DY3JetsToLL','DY4JetsToLL','DYJetsToLL_M-50-LO','Wg_AMCNLOFXFX','Wg_MADGRAPHMLM',
                                  'WpWmJJ_EWK_QCD_noTop','WpWmJJ_QCD_noTop','WpWmJJ_EWK_noTop',
                                  'TTTo2L2Nu_ext1','ST_t-channel_antitop','ST_t-channel_top',
                                  'ST_tW_antitop','ST_tW_top','TTJets', 
                                 ] ,
                },

  'MCl2vlooseCut' :       {
                  'isChain'    : True ,
                  'do4MC'      : True ,
                  'do4Data'    : False,
                  'subTargets' : ['do_l2vloose_Cut','puadder','baseW','wwNLL','genVariables','genMatchVariables','l2kin','l3kin','l4kin','do_dymvaHiggs','BWEwkSinglet','wwEWK','wzEWK','zzEWK'], #,'TopGenPt'],
                  'onlySample' : [
                                  # VBS
                                  'DYJetsToLL_M-50','DYJetsToLL_M-10to50ext3','DYJetsToLL_M-10to50','DYJetsToLL_M-50-LO-ext1',
                                  'WmWmJJ_EWK_powheg','WpWpJJ_EWK_powheg','WpWpJJ_EWK_QCD_aQGC','WpWpJJ_EWK_aQGC',
                                  'WpWpJJ_EWK','WpWpJJ_EWK_QCD','WpWpJJ_QCD','WW_DoubleScattering','WWTo2L2Nu_DoubleScattering','WLLJJToLNu_M-4to60_EWK_QCD','WLLJJToLNu_M-60_EWK_QCD','WGJJ',
                                  'WLLJJToLNu_M-60_EWK_4F','WLLJJToLNu_M-50_QCD_0Jet','WLLJJToLNu_M-50_QCD_1Jet','WLLJJToLNu_M-50_QCD_2Jet','WLLJJToLNu_M-50_QCD_3Jet','WLLJJToLNu_M-4To60_EWK_4F','WLLJJToLNu_M-4To50_QCD_0Jet','WLLJJToLNu_M-4To50_QCD_1Jet','WLLJJToLNu_M-4To50_QCD_2Jet','WLLJJToLNu_M-4To50_QCD_3Jet',
                                  'WZJJ_EWK_QCD','tZq_ll','ZZJJTo4L_EWK', 'ZZTo4L', 'WZTo3LNu',
                                  'TTToSemiLepton','TTToSemiLeptonic','DY2JetsToLL','DY3JetsToLL','DY4JetsToLL','DYJetsToLL_M-50-LO','Wg_AMCNLOFXFX','Wg_MADGRAPHMLM',
                                  'WpWmJJ_EWK_QCD_noTop','WpWmJJ_QCD_noTop','WpWmJJ_EWK_noTop',
                                  'TTTo2L2Nu_ext1','TTTo2L2Nu','ST_t-channel_antitop','ST_t-channel_top',
                                  'ST_tW_antitop','ST_tW_top','TTJets',
                                  'ggZZ4e','ggZZ4m','ggZZ4t','ggZZ2e2m','ggZZ2e2t','ggZZ2m2t',
                                 ] ,
                },


  'MCl1loose' :       {
                  'isChain'    : True ,
                  'do4MC'      : True ,
                  'do4Data'    : False,
                  'subTargets' : ['l1loose','puadder','baseW','wwNLL','genVariables','genMatchVariables','BWEwkSinglet','BWEwkSinglet_JHUGen698','TopGenPt'],
                  'onlySample' : [
                                  #### DY 
                                  'DYJetsToLL_M-10to50','DYJetsToLL_M-50','DYJetsToLL_M-10to50ext3','DYJetsToLL_M-50-LO',
                                  ####
                                  'WJetsToLNu','WJetsToLNu_HT100_200','WJetsToLNu_HT200_400','WJetsToLNu_HT400_600','WJetsToLNu_HT600_800',
                                  'WJetsToLNu_HT800_1200','WJetsToLNu_HT1200_2500','WJetsToLNu_HT2500_inf',
                                  ####
                                  'QCD_Pt-15to20_EMEnriched', 'QCD_Pt-20to30_EMEnriched', 'QCD_Pt-30to50_EMEnriched', 'QCD_Pt-50to80_EMEnriched',
                                  'QCD_Pt-20toInf_MuEnrichedPt15','QCD_Pt-30toInf_DoubleEMEnriched','QCD_Pt-15to20_MuEnrichedPt5',
                                  ####
                                  'QCD_Pt_15to20_bcToE','QCD_Pt_20to30_bcToE','QCD_Pt_30to80_bcToE','QCD_Pt_80to170_bcToE',
                                  'QCD_Pt_170to250_bcToE','QCD_Pt_250toInf_bcToE',
                                  ####
                                  'TT','TTJets',
 
                                 ] ,

                },

  'MCl1looseCut' :       {
                  'isChain'    : True ,
                  'do4MC'      : True ,
                  'do4Data'    : False,
                  'subTargets' : ['l1looseCut','puadder','baseW','wwNLL','genVariables','genMatchVariables'], #,'BWEwkSinglet','BWEwkSinglet_JHUGen698','TopGenPt'],
                  'onlySample' : [
                                  #### DY 
                                  'DYJetsToLL_M-10to50','DYJetsToLL_M-50','DYJetsToLL_M-10to50ext3','DYJetsToLL_M-50-LO','DYJetsToLL_M-50-LO-ext1',
                                  ####
                                  'WJetsToLNu','WJetsToLNu_HT100_200','WJetsToLNu_HT200_400','WJetsToLNu_HT400_600','WJetsToLNu_HT600_800',
                                  'WJetsToLNu_HT800_1200','WJetsToLNu_HT1200_2500','WJetsToLNu_HT2500_inf',
                                  ####
                                  'QCD_Pt-15to20_EMEnriched', 'QCD_Pt-20to30_EMEnriched', 'QCD_Pt-30to50_EMEnriched', 'QCD_Pt-50to80_EMEnriched','QCD_Pt-50to80_EMEnriched_ext1',
                                  'QCD_Pt-20toInf_MuEnrichedPt15','QCD_Pt-30toInf_DoubleEMEnriched','QCD_Pt-15to20_MuEnrichedPt5',
                                  ####
                                  'QCD_Pt_15to20_bcToE','QCD_Pt_20to30_bcToE','QCD_Pt_30to80_bcToE','QCD_Pt_80to170_bcToE',
                                  'QCD_Pt_170to250_bcToE','QCD_Pt_250toInf_bcToE',
                                  ####
                                  'TT','TTJets',

                                 ] ,

                },


  'MCl1vloose' :       {
                  'isChain'    : True ,
                  'do4MC'      : True ,
                  'do4Data'    : False,
                  'subTargets' : ['l1vloose','puadder','baseW','wwNLL','genVariables','genMatchVariables','BWEwkSinglet','BWEwkSinglet_JHUGen698','TopGenPt'],
                  'onlySample' : [
                                  #### DY 
                                  'DYJetsToLL_M-10to50','DYJetsToLL_M-50','DYJetsToLL_M-10to50ext3','DYJetsToLL_M-50-LO',
                                  'DY2JetsToLL','DY3JetsToLL','DY4JetsToLL','DYJetsToLL_M-50-LO',
                                  ####
                                  'WJetsToLNu','WJetsToLNu_HT100_200','WJetsToLNu_HT200_400','WJetsToLNu_HT400_600','WJetsToLNu_HT600_800',
                                  'WJetsToLNu_HT800_1200','WJetsToLNu_HT1200_2500','WJetsToLNu_HT2500_inf',
                                  ####
                                  'QCD_Pt-15to20_EMEnriched', 'QCD_Pt-20to30_EMEnriched', 'QCD_Pt-30to50_EMEnriched', 'QCD_Pt-50to80_EMEnriched',
                                  'QCD_Pt-20toInf_MuEnrichedPt15','QCD_Pt-30toInf_DoubleEMEnriched','QCD_Pt-15to20_MuEnrichedPt5',
                                  ####
                                  'QCD_Pt_15to20_bcToE','QCD_Pt_20to30_bcToE','QCD_Pt_30to80_bcToE','QCD_Pt_80to170_bcToE',
                                  'QCD_Pt_170to250_bcToE','QCD_Pt_250toInf_bcToE',
                                  ####
                                  'TT','TTJets',
                                 ] ,
                },

  'MCl1vlooseCut' :       {
                  'isChain'    : True ,
                  'do4MC'      : True ,
                  'do4Data'    : False,
                  'subTargets' : ['l1vlooseCut','puadder','baseW','wwNLL','genVariables','genMatchVariables'], # ,'BWEwkSinglet','BWEwkSinglet_JHUGen698','TopGenPt'],
                  'onlySample' : [
                                  #### DY 
                                  'DYJetsToLL_M-10to50','DYJetsToLL_M-50','DYJetsToLL_M-10to50ext3','DYJetsToLL_M-50-LO','DYJetsToLL_M-50-LO-ext1',
                                  'DY2JetsToLL','DY3JetsToLL','DY4JetsToLL','DYJetsToLL_M-50-LO',
                                  ####
                                  'WJetsToLNu','WJetsToLNu_HT100_200','WJetsToLNu_HT200_400','WJetsToLNu_HT400_600','WJetsToLNu_HT600_800',
                                  'WJetsToLNu_HT800_1200','WJetsToLNu_HT1200_2500','WJetsToLNu_HT2500_inf',
                                  ####
                                  'QCD_Pt-15to20_EMEnriched', 'QCD_Pt-20to30_EMEnriched', 'QCD_Pt-30to50_EMEnriched', 'QCD_Pt-50to80_EMEnriched','QCD_Pt-50to80_EMEnriched_ext1',
                                  'QCD_Pt-20toInf_MuEnrichedPt15','QCD_Pt-30toInf_DoubleEMEnriched','QCD_Pt-15to20_MuEnrichedPt5',
                                  ####
                                  'QCD_Pt_15to20_bcToE','QCD_Pt_20to30_bcToE','QCD_Pt_30to80_bcToE','QCD_Pt_80to170_bcToE',
                                  'QCD_Pt_170to250_bcToE','QCD_Pt_250toInf_bcToE',
                                  ####
                                  'TT','TTJets',
                                 ] ,
                },


  'MCWgStarsel' : {
                  'isChain'    : True ,
                  'do4MC'      : True ,
                  'do4Data'    : False,
                  'subTargets' : ['do_WgStarsel','puadder','pu2p6','pu4p3','pu6p3','baseW','wwNLL','genVariables','genMatchVariables','wwEWK','wzEWK','zzEWK'],
                  'onlySample' : [
                                   'Wg500','Wg_AMCNLOFXFX','WZTo3LNu','WgStarLNuEE','WgStarLNuMuMu','Wg_MADGRAPHMLM',
                                   'DYJetsToLL_M-10to50','DYJetsToLL_M-50','DYJetsToLL_M-10to50ext3', 
                                   'WZTo2L2Q','WZTo3LNu_mllmin01'  
                                 ]
                },

  'WgStarsel' : {
                  'isChain'    : True ,
                  'do4MC'      : False ,
                  'do4Data'    : True,
                  'subTargets' : ['do_WgStarsel'],
                  'onlySample' : [
                                   'Run2015C_16Dec2015_DoubleMuon' , 'Run2015C_16Dec2015_SingleElectron' , 'Run2015C_16Dec2015_SingleMuon',
                                   'Run2015D_16Dec2015_DoubleMuon' , 'Run2015D_16Dec2015_SingleElectron' , 'Run2015D_16Dec2015_SingleMuon',
                                   'Run2016B_PromptReco_DoubleMuon', 'Run2016B_PromptReco_SingleElectron', 'Run2016B_PromptReco_SingleMuon',
                                   'Run2016C_PromptReco_DoubleMuon', 'Run2016C_PromptReco_SingleElectron', 'Run2016C_PromptReco_SingleMuon',
                                   'Run2016D_PromptReco_DoubleMuon', 'Run2016D_PromptReco_SingleElectron', 'Run2016D_PromptReco_SingleMuon',
                                 ]
                },


  'bSFL2Eff'   :   {
                  'isChain'    : True ,
                  'do4MC'      : True ,
                  'do4Data'    : True ,
                  'subTargets' : ['bPogSF','TrigEff','IdIsoSC'],
                },
      

  'EpTCorr'       :  {
                  'isChain'    : True ,
                  'do4MC'      : True ,
                  'do4Data'    : True ,
                  'subTargets' : ['do_lpTCorrMC','do_lpTCorrData','l2kin','l3kin','l4kin','do_dymvaHiggs'],
                },

  'bSFKinFix'    : {
                  'isChain'    : True ,
                  'do4MC'      : True ,
                  'do4Data'    : False ,
                  'subTargets' : ['bPogSF','genMatchVariables'],
                },

  'KinFix'    : {
                  'isChain'    : True ,
                  'do4MC'      : True ,
                  'do4Data'    : False ,
                  'subTargets' : ['genMatchVariables'],
                },


  'filterjson'   : {
                  'isChain'    : True ,
                  'do4MC'      : False ,
                  'do4Data'    : True ,
                  'subTargets' : ['tagjson', 'selectjson'],
                } ,


   #tagjsonICHEP
  'ICHEPjson' : {
                  'isChain'    : True ,
                  'do4MC'      : False ,
                  'do4Data'    : True ,
                  'subTargets' : ['tagjsonICHEP', 'selectjson'],
                } ,
 


  'bSFLepEff'       :  {
                  'isChain'    : True ,
                  'do4MC'      : True ,
                  'do4Data'    : False ,
                  'subTargets' : ['bPogSF','TrigEff','IdIsoSC']
                    }, 

  'LepSFCut'    :  {
                  'isChain'    : True ,
                  'do4MC'      : True ,
                  'do4Data'    : False ,
                  'subTargets' : ['TrigEff_Cut','IdIsoSC_Cut']
                    },

  'bSFLepEffCut' :  {
                  'isChain'    : True ,
                  'do4MC'      : True ,
                  'do4Data'    : False ,
                  'subTargets' : ['bPogSF','TrigMakerMC','IdIsoSC_Cut']
                    },
  
  'bSFL2pTEffCut' : {
                  'isChain'    : True ,
                  'do4MC'      : True ,
                  'do4Data'    : False ,
                  'subTargets' : ['do_lpTCorrMC','bPogSF','TrigMakerMC','puRunPer','IdIsoSC_Cut','l2kin','l3kin','l4kin','do_dymvaHiggs'],
                  'XonlySample' : [
'DYJetsToLL_M-50','DYJetsToLL_M-10to50','WmWmJJ_EWK_powheg','WpWpJJ_EWK_powheg','WpWpJJ_EWK_aQGC','WpWpJJ_EWK','WpWpJJ_EWK_QCD','WpWpJJ_QCD','WWTo2L2Nu_DoubleScattering','WGJJ','WLLJJToLNu_M-60_EWK_4F','WLLJJToLNu_M-50_QCD_0Jet','WLLJJToLNu_M-50_QCD_1Jet','WLLJJToLNu_M-50_QCD_2Jet','WLLJJToLNu_M-50_QCD_3Jet','WLLJJToLNu_M-4To60_EWK_4F','WLLJJToLNu_M-4To50_QCD_0Jet','WLLJJToLNu_M-4To50_QCD_1Jet','WLLJJToLNu_M-4To50_QCD_2Jet','WLLJJToLNu_M-4To50_QCD_3Jet','tZq_ll','ZZJJTo4L_EWK', 'ZZTo4L', 'WZTo3LNu','Wg_AMCNLOFXFX','Wg_MADGRAPHMLM','TTTo2L2Nu','ST_t-channel_antitop','ST_t-channel_top','ST_tW_antitop','ST_tW_top','ggZZ4e','ggZZ4m','ggZZ4t','ggZZ2e2m','ggZZ2e2t','ggZZ2m2t']
                    },

  'bSFL1pTEffCut' : {
                  'isChain'    : True ,
                  'do4MC'      : True ,
                  'do4Data'    : False ,
                  'subTargets' : ['do_lpTCorrMC','bPogSF','TrigMakerMC','puRunPer','IdIsoSC_Cut','l2kin'],
                    },


  'puextra'      :   {
                  'isChain'    : True ,
                  'do4MC'      : True ,
                  'do4Data'    : False ,
                  'subTargets' : ['puW63mb','puW69mb']
                    }, 

  'bSFL2pTEff'   :   {
                  'isChain'    : True ,
                  'do4MC'      : True ,
                  'do4Data'    : True ,
                  'subTargets' : ['do_lpTCorrMC','do_lpTCorrData','bPogSF','TrigEff','IdIsoSC','l2kin','l3kin','l4kin'],
                  'XonlySample' : [
                                  # DY 
                                  'DYJetsToLL_M-10to50','DYJetsToLL_M-50','DYJetsToLL_M-10to50ext3', 'DYJetsToLL_M-50-LO' ,
                                  'DY2JetsToLL', 'DY3JetsToLL', 'DY4JetsToLL',
                                  'DYJetsToLL_M-50_HT-100to200_MLM' , 'DYJetsToLL_M-50_HT-100to200_MLM_ext' , 
                                  'DYJetsToLL_M-50_HT-200to400_MLM' , 
                                  'DYJetsToLL_M-50_HT-400to600_MLM' ,
                                  'DYJetsToLL_M-50_HT-600toInf_MLM' ,
                                  # WJets
                                  'WJetsToLNu',
                                  # Top
                                  'TTTo2L2Nu','TTWJetsToLNu',
                                  'ST_t-channel_antitop','ST_t-channel_top',
                                  'ST_tW_antitop','ST_tW_top',
                                  'TTToSemiLeptonic','TT',
                                  # VV (including WW) 
                                  'WWTo2L2Nu','GluGluWWTo2L2Nu_MCFM','GluGluWWTo2L2NuHiggs_MCFM',
                                  'WZTo3LNu','WZ',
                                  'ZZ','Zg',
                                  'ZZTo2L2Nu','ZZTo2L2Q','ZZTo4L',
                                  'Wg_AMCNLOFXFX', 
                                  'Wg_MADGRAPHMLM',
                                  'WZTo2L2Q',
                                  # VVV
                                  'WZZ','ZZZ','WWZ','WWW',
                                  # Higgs 
                                  'GluGluHToTauTau_M125', 'GluGluHToWWTo2L2Nu_M125','GluGluHToWWTo2L2NuPowheg_M125',
                                  'GluGluHToWWTo2L2Nu_alternative_M125','VBFHToWWTo2L2Nu_alternative_M125',
                                  'HWminusJ_HToTauTau_M125', 'HWminusJ_HToWW_M125',
                                  'HWplusJ_HToTauTau_M125', 'HWplusJ_HToWW_M125',
                                  'HZJ_HToTauTau_M125', 'HZJ_HToWW_M125',
                                  'VBFHToTauTau_M125', 'VBFHToWWTo2L2Nu_M125',
                                  'ggZH_HToWW_M125', # missing ggZHToTauTau
                                  'ttHJetToNonbb_M125',
                                  'GluGluHToZZTo4L_M125',
                                  'GluGluZH_HToWWTo2L2Nu_M120_noHLT','GluGluZH_HToWWTo2L2Nu_M125_noHLT','GluGluZH_HToWWTo2L2Nu_M130_noHLT',
                                  'HZJ_HToWWTo2L2Nu_M120_noHLT','HZJ_HToWWTo2L2Nu_M125_noHLT','HZJ_HToWWTo2L2Nu_M130_noHLT', 
                                  'HWplusJ_WToLNu_HToWWTo2L2Nu_M125','HWminusJ_WToLNu_HToWWTo2L2Nu_M125',
                                  # ... ggH High Mass
                                  'GluGluHToWWTo2L2Nu_M130',
                                  'GluGluHToWWTo2L2Nu_M135',
                                  'GluGluHToWWTo2L2Nu_M140',
                                  'GluGluHToWWTo2L2Nu_M145',
                                  'GluGluHToWWTo2L2Nu_M150',
                                  'GluGluHToWWTo2L2Nu_M155',
                                  'GluGluHToWWTo2L2Nu_M160',
                                  'GluGluHToWWTo2L2Nu_M165',
                                  'GluGluHToWWTo2L2Nu_M170',
                                  'GluGluHToWWTo2L2Nu_M175',
                                  'GluGluHToWWTo2L2Nu_M180',
                                  'GluGluHToWWTo2L2Nu_M190',
                                  'GluGluHToWWTo2L2Nu_M200',
                                  'GluGluHToWWTo2L2Nu_M210',
                                  'GluGluHToWWTo2L2Nu_M230',
                                  'GluGluHToWWTo2L2Nu_M250',
                                  'GluGluHToWWTo2L2Nu_M270',
                                  'GluGluHToWWTo2L2Nu_M300',
                                  'GluGluHToWWTo2L2Nu_M350',
                                  'GluGluHToWWTo2L2Nu_M400',
                                  'GluGluHToWWTo2L2Nu_M450',
                                  'GluGluHToWWTo2L2Nu_M500',
                                  'GluGluHToWWTo2L2Nu_M550',
                                  'GluGluHToWWTo2L2Nu_M600',
                                  'GluGluHToWWTo2L2Nu_M650',
                                  'GluGluHToWWTo2L2Nu_M700',
                                  'GluGluHToWWTo2L2Nu_M750',
                                  'GluGluHToWWTo2L2Nu_M750_NWA',
                                  'GluGluHToWWTo2L2Nu_M800',
                                  'GluGluHToWWTo2L2Nu_M900',
                                  'GluGluHToWWTo2L2Nu_M1000',
                                  # ... VBF High Mass
                                  'VBFHToWWTo2L2Nu_M130',
                                  'VBFHToWWTo2L2Nu_M135',
                                  'VBFHToWWTo2L2Nu_M140',
                                  'VBFHToWWTo2L2Nu_M145',
                                  'VBFHToWWTo2L2Nu_M150',
                                  'VBFHToWWTo2L2Nu_M155',
                                  'VBFHToWWTo2L2Nu_M160',
                                  'VBFHToWWTo2L2Nu_M165',
                                  'VBFHToWWTo2L2Nu_M170',
                                  'VBFHToWWTo2L2Nu_M175',
                                  'VBFHToWWTo2L2Nu_M180',
                                  'VBFHToWWTo2L2Nu_M190',
                                  'VBFHToWWTo2L2Nu_M200',
                                  'VBFHToWWTo2L2Nu_M210',
                                  'VBFHToWWTo2L2Nu_M230',
                                  'VBFHToWWTo2L2Nu_M250',
                                  'VBFHToWWTo2L2Nu_M270',
                                  'VBFHToWWTo2L2Nu_M300',
                                  'VBFHToWWTo2L2Nu_M350',
                                  'VBFHToWWTo2L2Nu_M400',
                                  'VBFHToWWTo2L2Nu_M450',
                                  'VBFHToWWTo2L2Nu_M500',
                                  'VBFHToWWTo2L2Nu_M550',
                                  'VBFHToWWTo2L2Nu_M600',
                                  'VBFHToWWTo2L2Nu_M650',
                                  'VBFHToWWTo2L2Nu_M700',
                                  'VBFHToWWTo2L2Nu_M750',
                                  'VBFHToWWTo2L2Nu_M750_NWA',
                                  'VBFHToWWTo2L2Nu_M800',
                                  'VBFHToWWTo2L2Nu_M900',
                                  'VBFHToWWTo2L2Nu_M1000',
                                  # PS
                                  'GluGluHToWWTo2L2NuHerwigPS_M125','VBFHToWWTo2L2NuHerwigPS_M125','WWTo2L2NuHerwigPS',
                                  # UE
                                  'GluGluHToWWTo2L2Nu_M125_CUETDown',
                                  'GluGluHToWWTo2L2Nu_M125_CUETUp',
                                  'VBFHToWWTo2L2Nu_M125_CUETDown',
                                  'VBFHToWWTo2L2Nu_M125_CUETUp',
                                  'WWTo2L2Nu_CUETDown',
                                  'WWTo2L2Nu_CUETUp',
 			  # VBS
 			  'WpWpJJ_EWK','WpWpJJ_EWK_QCD','WpWpJJ_QCD','WW_DoubleScattering','WWTo2L2Nu_DoubleScattering','WLLJJToLNu_M-4to60_EWK_QCD','WLLJJToLNu_M-60_EWK_QCD',
                          'WGJJ','EWKZ2Jets','TTToSemiLepton',
                                  # VBF 
                                  'WpWmJJ_EWK_QCD_noTop','WpWmJJ_QCD_noTop', 'WpWmJJ_EWK_noTop',
                                 ] ,
                }, 

  'bSFL2pTEffBStop'   :   {
                  'isChain'    : True ,
                  'do4MC'      : True ,
                  'do4Data'    : True ,
                  'subTargets' : ['do_lpTCorrMC','do_lpTCorrData','bPogSF','TrigEffB','IdIsoSCStop','l2kin','l3kin','l4kin'],
                  },


  'JESup'     :  {
                  'isChain'    : True ,
                  'do4MC'      : True ,
                  'do4Data'    : False,
                  'subTargets' : ['do_JESup','bPogSF','l2kin','l3kin','l4kin'],
                  'onlySample' : [
                                  # DY 
                                  'DYJetsToLL_M-10to50','DYJetsToLL_M-50','DYJetsToLL_M-10to50ext3',
                                  'DY2JetsToLL', 'DY3JetsToLL', 'DY4JetsToLL',
                                  'DYJetsToLL_M-50_HT-100to200_MLM' , 'DYJetsToLL_M-50_HT-100to200_MLM_ext' , 
                                  'DYJetsToLL_M-50_HT-200to400_MLM' , 
                                  'DYJetsToLL_M-50_HT-400to600_MLM' ,
                                  'DYJetsToLL_M-50_HT-600toInf_MLM' ,
                                  # Top
                                  'TTTo2L2Nu','TT','TTTo2L2Nu_ext1',
                                  'ST_t-channel_antitop','ST_t-channel_top',
                                  'ST_tW_antitop','ST_tW_top',
                                  # VV (including WW) 
                                  'WWTo2L2Nu','GluGluWWTo2L2Nu_MCFM','GluGluWWTo2L2NuHiggs_MCFM',
                                  'WZTo3LNu','WZ',
                                  'ZZ','Zg',
                                  'ZZTo2L2Nu','ZZTo2L2Q','ZZTo4L',
                                  'Wg_AMCNLOFXFX','WgStarLNuEE','WgStarLNuMuMu',
                                  'Wg_MADGRAPHMLM',
                                  # VVV
                                  'WZZ','ZZZ','WWZ','WWW',
                                  # Higgs 
                                  'GluGluHToTauTau_M125', 'GluGluHToWWTo2L2Nu_M125','GluGluHToWWTo2L2NuPowheg_M125',
                                  'GluGluHToWWTo2L2Nu_alternative_M125','VBFHToWWTo2L2Nu_alternative_M125',
                                  'HWminusJ_HToTauTau_M125', 'HWminusJ_HToWW_M125',
                                  'HWplusJ_HToTauTau_M125', 'HWplusJ_HToWW_M125',
                                  'HZJ_HToTauTau_M125', 'HZJ_HToWW_M125',
                                  'VBFHToTauTau_M125', 'VBFHToWWTo2L2Nu_M125',
                                  'ggZH_HToWW_M125', # missing ggZHToTauTau
                                  'ttHJetToNonbb_M125','TTWJetsToLNu',
                                  'GluGluHToZZTo4L_M125',
                                  'GluGluZH_HToWWTo2L2Nu_M120_noHLT','GluGluZH_HToWWTo2L2Nu_M125_noHLT','GluGluZH_HToWWTo2L2Nu_M130_noHLT',
                                  'HZJ_HToWWTo2L2Nu_M120_noHLT','HZJ_HToWWTo2L2Nu_M125_noHLT','HZJ_HToWWTo2L2Nu_M130_noHLT', 
                                  'HWplusJ_WToLNu_HToWWTo2L2Nu_M125','HWminusJ_WToLNu_HToWWTo2L2Nu_M125',

                                  # ... ggH High Mass
                                  'GluGluHToWWTo2L2Nu_M130',
                                  'GluGluHToWWTo2L2Nu_M135',
                                  'GluGluHToWWTo2L2Nu_M140',
                                  'GluGluHToWWTo2L2Nu_M145',
                                  'GluGluHToWWTo2L2Nu_M150',
                                  'GluGluHToWWTo2L2Nu_M155',
                                  'GluGluHToWWTo2L2Nu_M160',
                                  'GluGluHToWWTo2L2Nu_M165',
                                  'GluGluHToWWTo2L2Nu_M170',
                                  'GluGluHToWWTo2L2Nu_M175',
                                  'GluGluHToWWTo2L2Nu_M180',
                                  'GluGluHToWWTo2L2Nu_M190',
                                  'GluGluHToWWTo2L2Nu_M200',
                                  'GluGluHToWWTo2L2Nu_M210',
                                  'GluGluHToWWTo2L2Nu_M230',
                                  'GluGluHToWWTo2L2Nu_M250',
                                  'GluGluHToWWTo2L2Nu_M270',
                                  'GluGluHToWWTo2L2Nu_M300',
                                  'GluGluHToWWTo2L2Nu_M350',
                                  'GluGluHToWWTo2L2Nu_M400',
                                  'GluGluHToWWTo2L2Nu_M450',
                                  'GluGluHToWWTo2L2Nu_M500',
                                  'GluGluHToWWTo2L2Nu_M550',
                                  'GluGluHToWWTo2L2Nu_M600',
                                  'GluGluHToWWTo2L2Nu_M650',
                                  'GluGluHToWWTo2L2Nu_M700',
                                  'GluGluHToWWTo2L2Nu_M750',
                                  'GluGluHToWWTo2L2Nu_M750_NWA',
                                  'GluGluHToWWTo2L2Nu_M800',
                                  'GluGluHToWWTo2L2Nu_M900',
                                  'GluGluHToWWTo2L2Nu_M1000',
                                  # ... VBF High Mass
                                  'VBFHToWWTo2L2Nu_M130',
                                  'VBFHToWWTo2L2Nu_M135',
                                  'VBFHToWWTo2L2Nu_M140',
                                  'VBFHToWWTo2L2Nu_M145',
                                  'VBFHToWWTo2L2Nu_M150',
                                  'VBFHToWWTo2L2Nu_M155',
                                  'VBFHToWWTo2L2Nu_M160',
                                  'VBFHToWWTo2L2Nu_M165',
                                  'VBFHToWWTo2L2Nu_M170',
                                  'VBFHToWWTo2L2Nu_M175',
                                  'VBFHToWWTo2L2Nu_M180',
                                  'VBFHToWWTo2L2Nu_M190',
                                  'VBFHToWWTo2L2Nu_M200',
                                  'VBFHToWWTo2L2Nu_M210',
                                  'VBFHToWWTo2L2Nu_M230',
                                  'VBFHToWWTo2L2Nu_M250',
                                  'VBFHToWWTo2L2Nu_M270',
                                  'VBFHToWWTo2L2Nu_M300',
                                  'VBFHToWWTo2L2Nu_M350',
                                  'VBFHToWWTo2L2Nu_M400',
                                  'VBFHToWWTo2L2Nu_M450',
                                  'VBFHToWWTo2L2Nu_M500',
                                  'VBFHToWWTo2L2Nu_M550',
                                  'VBFHToWWTo2L2Nu_M600',
                                  'VBFHToWWTo2L2Nu_M650',
                                  'VBFHToWWTo2L2Nu_M700',
                                  'VBFHToWWTo2L2Nu_M750',
                                  'VBFHToWWTo2L2Nu_M750_NWA',
                                  'VBFHToWWTo2L2Nu_M800',
                                  'VBFHToWWTo2L2Nu_M900',
                                  'VBFHToWWTo2L2Nu_M1000',
                                  # VBF 
                                  'WpWmJJ_EWK_QCD_noTop','WpWmJJ_QCD_noTop', 'WpWmJJ_EWK_noTop',
                                 ] ,
                },

  'JESdo'     :  {
                  'isChain'    : True ,
                  'do4MC'      : True ,
                  'do4Data'    : False,
                  'subTargets' : ['do_JESdo','bPogSF','l2kin','l3kin','l4kin'],
                  'onlySample' : [
                                  # DY 
                                  'DYJetsToLL_M-10to50','DYJetsToLL_M-50','DYJetsToLL_M-10to50ext3',
                                  'DY2JetsToLL', 'DY3JetsToLL', 'DY4JetsToLL',
                                  'DYJetsToLL_M-50_HT-100to200_MLM' , 'DYJetsToLL_M-50_HT-100to200_MLM_ext' , 
                                  'DYJetsToLL_M-50_HT-200to400_MLM' , 
                                  'DYJetsToLL_M-50_HT-400to600_MLM' ,
                                  'DYJetsToLL_M-50_HT-600toInf_MLM' ,
                                  # Top
                                  'TTTo2L2Nu','TT','TTTo2L2Nu_ext1',
                                  'ST_t-channel_antitop','ST_t-channel_top',
                                  'ST_tW_antitop','ST_tW_top',
                                  # VV (including WW) 
                                  'WWTo2L2Nu','GluGluWWTo2L2Nu_MCFM','GluGluWWTo2L2NuHiggs_MCFM',
                                  'WZTo3LNu','WZ',
                                  'ZZ','Zg',
                                  'ZZTo2L2Nu','ZZTo2L2Q','ZZTo4L',
                                  'Wg_AMCNLOFXFX','WgStarLNuEE','WgStarLNuMuMu',
                                  'Wg_MADGRAPHMLM',
                                  # VVV
                                  'WZZ','ZZZ','WWZ','WWW',
                                  # Higgs 
                                  'GluGluHToTauTau_M125', 'GluGluHToWWTo2L2Nu_M125','GluGluHToWWTo2L2NuPowheg_M125',
                                  'GluGluHToWWTo2L2Nu_alternative_M125','VBFHToWWTo2L2Nu_alternative_M125',
                                  'HWminusJ_HToTauTau_M125', 'HWminusJ_HToWW_M125',
                                  'HWplusJ_HToTauTau_M125', 'HWplusJ_HToWW_M125',
                                  'HZJ_HToTauTau_M125', 'HZJ_HToWW_M125',
                                  'VBFHToTauTau_M125', 'VBFHToWWTo2L2Nu_M125',
                                  'ggZH_HToWW_M125', # missing ggZHToTauTau
                                  'ttHJetToNonbb_M125','TTWJetsToLNu',
                                  'GluGluHToZZTo4L_M125',
                                  'GluGluZH_HToWWTo2L2Nu_M120_noHLT','GluGluZH_HToWWTo2L2Nu_M125_noHLT','GluGluZH_HToWWTo2L2Nu_M130_noHLT',
                                  'HZJ_HToWWTo2L2Nu_M120_noHLT','HZJ_HToWWTo2L2Nu_M125_noHLT','HZJ_HToWWTo2L2Nu_M130_noHLT', 
                                  'HWplusJ_WToLNu_HToWWTo2L2Nu_M125','HWminusJ_WToLNu_HToWWTo2L2Nu_M125',

                                  # ... ggH High Mass
                                  'GluGluHToWWTo2L2Nu_M130',
                                  'GluGluHToWWTo2L2Nu_M135',
                                  'GluGluHToWWTo2L2Nu_M140',
                                  'GluGluHToWWTo2L2Nu_M145',
                                  'GluGluHToWWTo2L2Nu_M150',
                                  'GluGluHToWWTo2L2Nu_M155',
                                  'GluGluHToWWTo2L2Nu_M160',
                                  'GluGluHToWWTo2L2Nu_M165',
                                  'GluGluHToWWTo2L2Nu_M170',
                                  'GluGluHToWWTo2L2Nu_M175',
                                  'GluGluHToWWTo2L2Nu_M180',
                                  'GluGluHToWWTo2L2Nu_M190',
                                  'GluGluHToWWTo2L2Nu_M200',
                                  'GluGluHToWWTo2L2Nu_M210',
                                  'GluGluHToWWTo2L2Nu_M230',
                                  'GluGluHToWWTo2L2Nu_M250',
                                  'GluGluHToWWTo2L2Nu_M270',
                                  'GluGluHToWWTo2L2Nu_M300',
                                  'GluGluHToWWTo2L2Nu_M350',
                                  'GluGluHToWWTo2L2Nu_M400',
                                  'GluGluHToWWTo2L2Nu_M450',
                                  'GluGluHToWWTo2L2Nu_M500',
                                  'GluGluHToWWTo2L2Nu_M550',
                                  'GluGluHToWWTo2L2Nu_M600',
                                  'GluGluHToWWTo2L2Nu_M650',
                                  'GluGluHToWWTo2L2Nu_M700',
                                  'GluGluHToWWTo2L2Nu_M750',
                                  'GluGluHToWWTo2L2Nu_M750_NWA',
                                  'GluGluHToWWTo2L2Nu_M800',
                                  'GluGluHToWWTo2L2Nu_M900',
                                  'GluGluHToWWTo2L2Nu_M1000',
                                  # ... VBF High Mass
                                  'VBFHToWWTo2L2Nu_M130',
                                  'VBFHToWWTo2L2Nu_M135',
                                  'VBFHToWWTo2L2Nu_M140',
                                  'VBFHToWWTo2L2Nu_M145',
                                  'VBFHToWWTo2L2Nu_M150',
                                  'VBFHToWWTo2L2Nu_M155',
                                  'VBFHToWWTo2L2Nu_M160',
                                  'VBFHToWWTo2L2Nu_M165',
                                  'VBFHToWWTo2L2Nu_M170',
                                  'VBFHToWWTo2L2Nu_M175',
                                  'VBFHToWWTo2L2Nu_M180',
                                  'VBFHToWWTo2L2Nu_M190',
                                  'VBFHToWWTo2L2Nu_M200',
                                  'VBFHToWWTo2L2Nu_M210',
                                  'VBFHToWWTo2L2Nu_M230',
                                  'VBFHToWWTo2L2Nu_M250',
                                  'VBFHToWWTo2L2Nu_M270',
                                  'VBFHToWWTo2L2Nu_M300',
                                  'VBFHToWWTo2L2Nu_M350',
                                  'VBFHToWWTo2L2Nu_M400',
                                  'VBFHToWWTo2L2Nu_M450',
                                  'VBFHToWWTo2L2Nu_M500',
                                  'VBFHToWWTo2L2Nu_M550',
                                  'VBFHToWWTo2L2Nu_M600',
                                  'VBFHToWWTo2L2Nu_M650',
                                  'VBFHToWWTo2L2Nu_M700',
                                  'VBFHToWWTo2L2Nu_M750',
                                  'VBFHToWWTo2L2Nu_M750_NWA',
                                  'VBFHToWWTo2L2Nu_M800',
                                  'VBFHToWWTo2L2Nu_M900',
                                  'VBFHToWWTo2L2Nu_M1000',
                                  # VBF 
                                  'WpWmJJ_EWK_QCD_noTop','WpWmJJ_QCD_noTop', 'WpWmJJ_EWK_noTop',
                                 ] ,
                },


  'JESMaxup'     :  {
                  'isChain'    : True ,
                  'do4MC'      : True ,
                  'do4Data'    : False,
                  'subTargets' : ['do_JESMaxup','bPogSF','l2kin','l3kin','l4kin'],
                  'onlySample' : [
                                  # DY 
                                  'DYJetsToLL_M-10to50','DYJetsToLL_M-50','DYJetsToLL_M-10to50ext3','DYJetsToLL_M-50-LO','DYJetsToLL_M-50-LO-ext1',
                                  'DY2JetsToLL', 'DY3JetsToLL', 'DY4JetsToLL',
                                  'DYJetsToLL_M-50_HT-100to200_MLM' , 'DYJetsToLL_M-50_HT-100to200_MLM_ext' , 
                                  'DYJetsToLL_M-50_HT-200to400_MLM' , 
                                  'DYJetsToLL_M-50_HT-400to600_MLM' ,
                                  'DYJetsToLL_M-50_HT-600toInf_MLM' ,
                                  'DYJetsToLL_M-50_HT-100to200_ext1','DYJetsToLL_M-50_HT-100to200_ext1','DYJetsToLL_M-50_HT-200to400_ext1',
                                  'DYJetsToLL_M-50_HT-200to400_ext1','DYJetsToLL_M-50_HT-200to400_ext1','DYJetsToLL_M-50_HT-600toInf_ext1',
                                  'DYJetsToLL_M-50_HT-600toInf', 'DYJetsToTT_MuEle_M-50',
                                  # Top
                                  'TTTo2L2Nu','TT','TTTo2L2Nu_ext1',
                                  'ST_t-channel_antitop','ST_t-channel_top',
                                  'ST_tW_antitop','ST_tW_top',
                                  # VV (including WW) 
                                  'WWTo2L2Nu','GluGluWWTo2L2Nu_MCFM','GluGluWWTo2L2NuHiggs_MCFM',
                                  'WZTo3LNu','WZ',
                                  'ZZ','Zg',
                                  'ZZTo2L2Nu','ZZTo2L2Q','ZZTo4L',
                                  'Wg_AMCNLOFXFX','WgStarLNuEE','WgStarLNuMuMu',
                                  'Wg_MADGRAPHMLM','WZTo3LNu_mllmin01',
                                  'ggZZ2e2t', 'ggZZ4t', 'ggZZ4m', 'ggZZ2m2t', 'ggZZ2e2m', 'ggZZ4e',
                                  # VVV
                                  'WZZ','ZZZ','WWZ','WWW',
                                  # Higgs 
                                  'GluGluHToTauTau_M125', 'GluGluHToWWTo2L2Nu_M125','GluGluHToWWTo2L2NuPowheg_M125',
                                  'GluGluHToWWTo2L2Nu_alternative_M125','VBFHToWWTo2L2Nu_alternative_M125',
                                  'HWminusJ_HToTauTau_M125', 'HWminusJ_HToWW_M125',
                                  'HWplusJ_HToTauTau_M125', 'HWplusJ_HToWW_M125',
                                  'HZJ_HToTauTau_M125', 'HZJ_HToWW_M125',
                                  'VBFHToTauTau_M125', 'VBFHToWWTo2L2Nu_M125',
                                  'ggZH_HToWW_M125', # missing ggZHToTauTau
                                  'ttHJetToNonbb_M125','TTWJetsToLNu',
                                  'GluGluHToZZTo4L_M125',
                                  'GluGluZH_HToWWTo2L2Nu_M120_noHLT','GluGluZH_HToWWTo2L2Nu_M125_noHLT','GluGluZH_HToWWTo2L2Nu_M130_noHLT',
                                  'HZJ_HToWWTo2L2Nu_M120_noHLT','HZJ_HToWWTo2L2Nu_M125_noHLT','HZJ_HToWWTo2L2Nu_M130_noHLT', 
                                  'HWplusJ_WToLNu_HToWWTo2L2Nu_M125','HWminusJ_WToLNu_HToWWTo2L2Nu_M125',

                                  # ... ggH High Mass
                                  'GluGluHToWWTo2L2Nu_M130',
                                  'GluGluHToWWTo2L2Nu_M135',
                                  'GluGluHToWWTo2L2Nu_M140',
                                  'GluGluHToWWTo2L2Nu_M145',
                                  'GluGluHToWWTo2L2Nu_M150',
                                  'GluGluHToWWTo2L2Nu_M155',
                                  'GluGluHToWWTo2L2Nu_M160',
                                  'GluGluHToWWTo2L2Nu_M165',
                                  'GluGluHToWWTo2L2Nu_M170',
                                  'GluGluHToWWTo2L2Nu_M175',
                                  'GluGluHToWWTo2L2Nu_M180',
                                  'GluGluHToWWTo2L2Nu_M190',
                                  'GluGluHToWWTo2L2Nu_M200',
                                  'GluGluHToWWTo2L2Nu_M210',
                                  'GluGluHToWWTo2L2Nu_M230',
                                  'GluGluHToWWTo2L2Nu_M250',
                                  'GluGluHToWWTo2L2Nu_M270',
                                  'GluGluHToWWTo2L2Nu_M300',
                                  'GluGluHToWWTo2L2Nu_M350',
                                  'GluGluHToWWTo2L2Nu_M400',
                                  'GluGluHToWWTo2L2Nu_M450',
                                  'GluGluHToWWTo2L2Nu_M500',
                                  'GluGluHToWWTo2L2Nu_M550',
                                  'GluGluHToWWTo2L2Nu_M600',
                                  'GluGluHToWWTo2L2Nu_M650',
                                  'GluGluHToWWTo2L2Nu_M700',
                                  'GluGluHToWWTo2L2Nu_M750',
                                  'GluGluHToWWTo2L2Nu_M750_NWA',
                                  'GluGluHToWWTo2L2Nu_M800',
                                  'GluGluHToWWTo2L2Nu_M900',
                                  'GluGluHToWWTo2L2Nu_M1000',
                                  # ... VBF High Mass
                                  'VBFHToWWTo2L2Nu_M130',
                                  'VBFHToWWTo2L2Nu_M135',
                                  'VBFHToWWTo2L2Nu_M140',
                                  'VBFHToWWTo2L2Nu_M145',
                                  'VBFHToWWTo2L2Nu_M150',
                                  'VBFHToWWTo2L2Nu_M155',
                                  'VBFHToWWTo2L2Nu_M160',
                                  'VBFHToWWTo2L2Nu_M165',
                                  'VBFHToWWTo2L2Nu_M170',
                                  'VBFHToWWTo2L2Nu_M175',
                                  'VBFHToWWTo2L2Nu_M180',
                                  'VBFHToWWTo2L2Nu_M190',
                                  'VBFHToWWTo2L2Nu_M200',
                                  'VBFHToWWTo2L2Nu_M210',
                                  'VBFHToWWTo2L2Nu_M230',
                                  'VBFHToWWTo2L2Nu_M250',
                                  'VBFHToWWTo2L2Nu_M270',
                                  'VBFHToWWTo2L2Nu_M300',
                                  'VBFHToWWTo2L2Nu_M350',
                                  'VBFHToWWTo2L2Nu_M400',
                                  'VBFHToWWTo2L2Nu_M450',
                                  'VBFHToWWTo2L2Nu_M500',
                                  'VBFHToWWTo2L2Nu_M550',
                                  'VBFHToWWTo2L2Nu_M600',
                                  'VBFHToWWTo2L2Nu_M650',
                                  'VBFHToWWTo2L2Nu_M700',
                                  'VBFHToWWTo2L2Nu_M750',
                                  'VBFHToWWTo2L2Nu_M750_NWA',
                                  'VBFHToWWTo2L2Nu_M800',
                                  'VBFHToWWTo2L2Nu_M900',
                                  'VBFHToWWTo2L2Nu_M1000',
                                  # VBF 
                                  'WpWmJJ_EWK_QCD_noTop','WpWmJJ_QCD_noTop', 'WpWmJJ_EWK_noTop',
                                  'WGJJ', 'WLLJJToLNu_M-4to60_EWK_QCD', 'WLLJJToLNu_M-60_EWK_QCD', 'WpWpJJ_EWK_QCD',
                                  'WW_DoubleScattering','WWTo2L2Nu_DoubleScattering', 'WmWmJJ_EWK_powheg', 'WpWpJJ_EWK', 'WpWpJJ_QCD', 
                                  # stop

                                  'latino_T2bW',
                                  'latino_T2tt_dM10to80',
                                  'latino_T2tt_dM10to80_genHT160_genMET80',
                                  'latino_T2tt_mStop-150to250',
                                  'latino_T2tt_mStop-250to350',
                                  'latino_T2tt_mStop-350to400',
                                  'latino_T2tt_mStop-400to1200',
                                  'latino_T2tt_mStop425_mLSP325',
                                  'latino_T2tt_mStop500_mLSP325',
                                  'latino_T2tt_mStop850_mLSP100',

   
                                  # monoH
                                  'monoH_2HDM_MZp-1000_MA0-400',
                                  'monoH_2HDM_MZp-1000_MA0-500',
                                  'monoH_2HDM_MZp-1000_MA0-600',
                                  'monoH_2HDM_MZp-1000_MA0-700',
                                  'monoH_2HDM_MZp-1000_MA0-800',
                                  'monoH_2HDM_MZp-1200_MA0-400',
                                  'monoH_2HDM_MZp-1200_MA0-500',
                                  'monoH_2HDM_MZp-1200_MA0-600',
                                  'monoH_2HDM_MZp-1200_MA0-700',
                                  'monoH_2HDM_MZp-1200_MA0-800',
                                  'monoH_2HDM_MZp-1400_MA0-400',
                                  'monoH_2HDM_MZp-1400_MA0-500',
                                  'monoH_2HDM_MZp-1400_MA0-600',
                                  'monoH_2HDM_MZp-1400_MA0-700',
                                  'monoH_2HDM_MZp-1400_MA0-800',
                                  'monoH_2HDM_MZp-1700_MA0-400',
                                  'monoH_2HDM_MZp-1700_MA0-500',
                                  'monoH_2HDM_MZp-1700_MA0-600',
                                  'monoH_2HDM_MZp-1700_MA0-700',
                                  'monoH_2HDM_MZp-1700_MA0-800',
                                  'monoH_2HDM_MZp-2000_MA0-400',
                                  'monoH_2HDM_MZp-2000_MA0-500',
                                  'monoH_2HDM_MZp-2000_MA0-600',
                                  'monoH_2HDM_MZp-2000_MA0-700',
                                  'monoH_2HDM_MZp-2000_MA0-800',
                                  'monoH_2HDM_MZp-2500_MA0-400',
                                  'monoH_2HDM_MZp-2500_MA0-500',
                                  'monoH_2HDM_MZp-2500_MA0-600',
                                  'monoH_2HDM_MZp-2500_MA0-700',
                                  'monoH_2HDM_MZp-2500_MA0-800',
                                  'monoH_2HDM_MZp-600_MA0-400',
                                  'monoH_2HDM_MZp-800_MA0-400',
                                  'monoH_2HDM_MZp-800_MA0-500',
                                  'monoH_2HDM_MZp-800_MA0-600',
                                  'monoH_2HDM_MZp-1000_MA0-300',
                                  'monoH_2HDM_MZp-1200_MA0-300',
                                  'monoH_2HDM_MZp-1400_MA0-300',
                                  'monoH_2HDM_MZp-1700_MA0-300',
                                  'monoH_2HDM_MZp-2000_MA0-300',
                                  'monoH_2HDM_MZp-2500_MA0-300',
                                  'monoH_2HDM_MZp-600_MA0-300',
                                  'monoH_2HDM_MZp-800_MA0-300',

                                  'monoH_ZpBaryonic_MZp-10000_MChi-1000',
                                  'monoH_ZpBaryonic_MZp-10000_MChi-10',
                                  'monoH_ZpBaryonic_MZp-10000_MChi-150',
                                  'monoH_ZpBaryonic_MZp-10000_MChi-1',
                                  'monoH_ZpBaryonic_MZp-10000_MChi-50',
                                  'monoH_ZpBaryonic_MZp-1000_MChi-1000',
                                  'monoH_ZpBaryonic_MZp-1000_MChi-1',
                                  'monoH_ZpBaryonic_MZp-100_MChi-10',
                                  'monoH_ZpBaryonic_MZp-100_MChi-1',
                                  'monoH_ZpBaryonic_MZp-10_MChi-1000',
                                  'monoH_ZpBaryonic_MZp-10_MChi-10',
                                  'monoH_ZpBaryonic_MZp-10_MChi-150',
                                  'monoH_ZpBaryonic_MZp-10_MChi-1',
                                  'monoH_ZpBaryonic_MZp-10_MChi-500',
                                  'monoH_ZpBaryonic_MZp-10_MChi-50',
                                  'monoH_ZpBaryonic_MZp-15_MChi-10',
                                  'monoH_ZpBaryonic_MZp-2000_MChi-1',
                                  'monoH_ZpBaryonic_MZp-2000_MChi-500',
                                  'monoH_ZpBaryonic_MZp-200_MChi-150',
                                  'monoH_ZpBaryonic_MZp-200_MChi-1',
                                  'monoH_ZpBaryonic_MZp-200_MChi-50',
                                  'monoH_ZpBaryonic_MZp-20_MChi-1',
                                  'monoH_ZpBaryonic_MZp-295_MChi-150',
                                  'monoH_ZpBaryonic_MZp-300_MChi-1',
                                  'monoH_ZpBaryonic_MZp-300_MChi-50',
                                  'monoH_ZpBaryonic_MZp-500_MChi-150',
                                  'monoH_ZpBaryonic_MZp-500_MChi-1',
                                  'monoH_ZpBaryonic_MZp-500_MChi-500',
                                  'monoH_ZpBaryonic_MZp-50_MChi-10',
                                  'monoH_ZpBaryonic_MZp-50_MChi-1',
                                  'monoH_ZpBaryonic_MZp-95_MChi-50',
                                  'monoH_ZpBaryonic_MZp-995_MChi-500',

                               ],               
                },

  'JESMaxdo'     :  {
                  'isChain'    : True ,
                  'do4MC'      : True ,
                  'do4Data'    : False,
                  'subTargets' : ['do_JESMaxdo','bPogSF','l2kin','l3kin','l4kin'],
                  'onlySample' : [
                                  # DY 
                                  'DYJetsToLL_M-10to50','DYJetsToLL_M-50','DYJetsToLL_M-10to50ext3','DYJetsToLL_M-50-LO','DYJetsToLL_M-50-LO-ext1',
                                  'DY2JetsToLL', 'DY3JetsToLL', 'DY4JetsToLL',
                                  'DYJetsToLL_M-50_HT-100to200_MLM' , 'DYJetsToLL_M-50_HT-100to200_MLM_ext' , 
                                  'DYJetsToLL_M-50_HT-200to400_MLM' , 
                                  'DYJetsToLL_M-50_HT-400to600_MLM' ,
                                  'DYJetsToLL_M-50_HT-600toInf_MLM' ,
                                  'DYJetsToLL_M-50_HT-100to200_ext1','DYJetsToLL_M-50_HT-100to200_ext1','DYJetsToLL_M-50_HT-200to400_ext1',
                                  'DYJetsToLL_M-50_HT-200to400_ext1','DYJetsToLL_M-50_HT-200to400_ext1','DYJetsToLL_M-50_HT-600toInf_ext1',
                                  'DYJetsToLL_M-50_HT-600toInf', 'DYJetsToTT_MuEle_M-50',
                                  # Top
                                  'TTTo2L2Nu','TT','TTTo2L2Nu_ext1',
                                  'ST_t-channel_antitop','ST_t-channel_top',
                                  'ST_tW_antitop','ST_tW_top',
                                  # VV (including WW) 
                                  'WWTo2L2Nu','GluGluWWTo2L2Nu_MCFM','GluGluWWTo2L2NuHiggs_MCFM',
                                  'WZTo3LNu','WZ',
                                  'ZZ','Zg',
                                  'ZZTo2L2Nu','ZZTo2L2Q','ZZTo4L',
                                  'Wg_AMCNLOFXFX','WgStarLNuEE','WgStarLNuMuMu',
                                  'Wg_MADGRAPHMLM','WZTo3LNu_mllmin01',
                                  'ggZZ2e2t', 'ggZZ4t', 'ggZZ4m', 'ggZZ2m2t', 'ggZZ2e2m', 'ggZZ4e',
                                  # VVV
                                  'WZZ','ZZZ','WWZ','WWW',
                                  # Higgs 
                                  'GluGluHToTauTau_M125', 'GluGluHToWWTo2L2Nu_M125','GluGluHToWWTo2L2NuPowheg_M125',
                                  'GluGluHToWWTo2L2Nu_alternative_M125','VBFHToWWTo2L2Nu_alternative_M125',
                                  'HWminusJ_HToTauTau_M125', 'HWminusJ_HToWW_M125',
                                  'HWplusJ_HToTauTau_M125', 'HWplusJ_HToWW_M125',
                                  'HZJ_HToTauTau_M125', 'HZJ_HToWW_M125',
                                  'VBFHToTauTau_M125', 'VBFHToWWTo2L2Nu_M125',
                                  'ggZH_HToWW_M125', # missing ggZHToTauTau
                                  'ttHJetToNonbb_M125','TTWJetsToLNu',
                                  'GluGluHToZZTo4L_M125',
                                  'GluGluZH_HToWWTo2L2Nu_M120_noHLT','GluGluZH_HToWWTo2L2Nu_M125_noHLT','GluGluZH_HToWWTo2L2Nu_M130_noHLT',
                                  'HZJ_HToWWTo2L2Nu_M120_noHLT','HZJ_HToWWTo2L2Nu_M125_noHLT','HZJ_HToWWTo2L2Nu_M130_noHLT', 
                                  'HWplusJ_WToLNu_HToWWTo2L2Nu_M125','HWminusJ_WToLNu_HToWWTo2L2Nu_M125',
                                  # ... ggH High Mass
                                  'GluGluHToWWTo2L2Nu_M130',
                                  'GluGluHToWWTo2L2Nu_M135',
                                  'GluGluHToWWTo2L2Nu_M140',
                                  'GluGluHToWWTo2L2Nu_M145',
                                  'GluGluHToWWTo2L2Nu_M150',
                                  'GluGluHToWWTo2L2Nu_M155',
                                  'GluGluHToWWTo2L2Nu_M160',
                                  'GluGluHToWWTo2L2Nu_M165',
                                  'GluGluHToWWTo2L2Nu_M170',
                                  'GluGluHToWWTo2L2Nu_M175',
                                  'GluGluHToWWTo2L2Nu_M180',
                                  'GluGluHToWWTo2L2Nu_M190',
                                  'GluGluHToWWTo2L2Nu_M200',
                                  'GluGluHToWWTo2L2Nu_M210',
                                  'GluGluHToWWTo2L2Nu_M230',
                                  'GluGluHToWWTo2L2Nu_M250',
                                  'GluGluHToWWTo2L2Nu_M270',
                                  'GluGluHToWWTo2L2Nu_M300',
                                  'GluGluHToWWTo2L2Nu_M350',
                                  'GluGluHToWWTo2L2Nu_M400',
                                  'GluGluHToWWTo2L2Nu_M450',
                                  'GluGluHToWWTo2L2Nu_M500',
                                  'GluGluHToWWTo2L2Nu_M550',
                                  'GluGluHToWWTo2L2Nu_M600',
                                  'GluGluHToWWTo2L2Nu_M650',
                                  'GluGluHToWWTo2L2Nu_M700',
                                  'GluGluHToWWTo2L2Nu_M750',
                                  'GluGluHToWWTo2L2Nu_M750_NWA',
                                  'GluGluHToWWTo2L2Nu_M800',
                                  'GluGluHToWWTo2L2Nu_M900',
                                  'GluGluHToWWTo2L2Nu_M1000',
                                  # ... VBF High Mass
                                  'VBFHToWWTo2L2Nu_M130',
                                  'VBFHToWWTo2L2Nu_M135',
                                  'VBFHToWWTo2L2Nu_M140',
                                  'VBFHToWWTo2L2Nu_M145',
                                  'VBFHToWWTo2L2Nu_M150',
                                  'VBFHToWWTo2L2Nu_M155',
                                  'VBFHToWWTo2L2Nu_M160',
                                  'VBFHToWWTo2L2Nu_M165',
                                  'VBFHToWWTo2L2Nu_M170',
                                  'VBFHToWWTo2L2Nu_M175',
                                  'VBFHToWWTo2L2Nu_M180',
                                  'VBFHToWWTo2L2Nu_M190',
                                  'VBFHToWWTo2L2Nu_M200',
                                  'VBFHToWWTo2L2Nu_M210',
                                  'VBFHToWWTo2L2Nu_M230',
                                  'VBFHToWWTo2L2Nu_M250',
                                  'VBFHToWWTo2L2Nu_M270',
                                  'VBFHToWWTo2L2Nu_M300',
                                  'VBFHToWWTo2L2Nu_M350',
                                  'VBFHToWWTo2L2Nu_M400',
                                  'VBFHToWWTo2L2Nu_M450',
                                  'VBFHToWWTo2L2Nu_M500',
                                  'VBFHToWWTo2L2Nu_M550',
                                  'VBFHToWWTo2L2Nu_M600',
                                  'VBFHToWWTo2L2Nu_M650',
                                  'VBFHToWWTo2L2Nu_M700',
                                  'VBFHToWWTo2L2Nu_M750',
                                  'VBFHToWWTo2L2Nu_M750_NWA',
                                  'VBFHToWWTo2L2Nu_M800',
                                  'VBFHToWWTo2L2Nu_M900',
                                  'VBFHToWWTo2L2Nu_M1000',
                                  # VBF 
                                  'WpWmJJ_EWK_QCD_noTop','WpWmJJ_QCD_noTop', 'WpWmJJ_EWK_noTop',
                                  'WGJJ', 'WLLJJToLNu_M-4to60_EWK_QCD', 'WLLJJToLNu_M-60_EWK_QCD', 'WpWpJJ_EWK_QCD',
                                  'WW_DoubleScattering','WWTo2L2Nu_DoubleScattering', 'WmWmJJ_EWK_powheg', 'WpWpJJ_EWK', 'WpWpJJ_QCD', 


                                  # stop

                                  'latino_T2bW',
                                  'latino_T2tt_dM10to80',
                                  'latino_T2tt_dM10to80_genHT160_genMET80',
                                  'latino_T2tt_mStop-150to250',
                                  'latino_T2tt_mStop-250to350',
                                  'latino_T2tt_mStop-350to400',
                                  'latino_T2tt_mStop-400to1200',
                                  'latino_T2tt_mStop425_mLSP325',
                                  'latino_T2tt_mStop500_mLSP325',
                                  'latino_T2tt_mStop850_mLSP100',


                                  # monoH
                                  'monoH_2HDM_MZp-1000_MA0-400',
                                  'monoH_2HDM_MZp-1000_MA0-500',
                                  'monoH_2HDM_MZp-1000_MA0-600',
                                  'monoH_2HDM_MZp-1000_MA0-700',
                                  'monoH_2HDM_MZp-1000_MA0-800',
                                  'monoH_2HDM_MZp-1200_MA0-400',
                                  'monoH_2HDM_MZp-1200_MA0-500',
                                  'monoH_2HDM_MZp-1200_MA0-600',
                                  'monoH_2HDM_MZp-1200_MA0-700',
                                  'monoH_2HDM_MZp-1200_MA0-800',
                                  'monoH_2HDM_MZp-1400_MA0-400',
                                  'monoH_2HDM_MZp-1400_MA0-500',
                                  'monoH_2HDM_MZp-1400_MA0-600',
                                  'monoH_2HDM_MZp-1400_MA0-700',
                                  'monoH_2HDM_MZp-1400_MA0-800',
                                  'monoH_2HDM_MZp-1700_MA0-400',
                                  'monoH_2HDM_MZp-1700_MA0-500',
                                  'monoH_2HDM_MZp-1700_MA0-600',
                                  'monoH_2HDM_MZp-1700_MA0-700',
                                  'monoH_2HDM_MZp-1700_MA0-800',
                                  'monoH_2HDM_MZp-2000_MA0-400',
                                  'monoH_2HDM_MZp-2000_MA0-500',
                                  'monoH_2HDM_MZp-2000_MA0-600',
                                  'monoH_2HDM_MZp-2000_MA0-700',
                                  'monoH_2HDM_MZp-2000_MA0-800',
                                  'monoH_2HDM_MZp-2500_MA0-400',
                                  'monoH_2HDM_MZp-2500_MA0-500',
                                  'monoH_2HDM_MZp-2500_MA0-600',
                                  'monoH_2HDM_MZp-2500_MA0-700',
                                  'monoH_2HDM_MZp-2500_MA0-800',
                                  'monoH_2HDM_MZp-600_MA0-400',
                                  'monoH_2HDM_MZp-800_MA0-400',
                                  'monoH_2HDM_MZp-800_MA0-500',
                                  'monoH_2HDM_MZp-800_MA0-600',
                                  'monoH_2HDM_MZp-1000_MA0-300',
                                  'monoH_2HDM_MZp-1200_MA0-300',
                                  'monoH_2HDM_MZp-1400_MA0-300',
                                  'monoH_2HDM_MZp-1700_MA0-300',
                                  'monoH_2HDM_MZp-2000_MA0-300',
                                  'monoH_2HDM_MZp-2500_MA0-300',
                                  'monoH_2HDM_MZp-600_MA0-300',
                                  'monoH_2HDM_MZp-800_MA0-300',

                                  'monoH_ZpBaryonic_MZp-10000_MChi-1000',
                                  'monoH_ZpBaryonic_MZp-10000_MChi-10',
                                  'monoH_ZpBaryonic_MZp-10000_MChi-150',
                                  'monoH_ZpBaryonic_MZp-10000_MChi-1',
                                  'monoH_ZpBaryonic_MZp-10000_MChi-50',
                                  'monoH_ZpBaryonic_MZp-1000_MChi-1000',
                                  'monoH_ZpBaryonic_MZp-1000_MChi-1',
                                  'monoH_ZpBaryonic_MZp-100_MChi-10',
                                  'monoH_ZpBaryonic_MZp-100_MChi-1',
                                  'monoH_ZpBaryonic_MZp-10_MChi-1000',
                                  'monoH_ZpBaryonic_MZp-10_MChi-10',
                                  'monoH_ZpBaryonic_MZp-10_MChi-150',
                                  'monoH_ZpBaryonic_MZp-10_MChi-1',
                                  'monoH_ZpBaryonic_MZp-10_MChi-500',
                                  'monoH_ZpBaryonic_MZp-10_MChi-50',
                                  'monoH_ZpBaryonic_MZp-15_MChi-10',
                                  'monoH_ZpBaryonic_MZp-2000_MChi-1',
                                  'monoH_ZpBaryonic_MZp-2000_MChi-500',
                                  'monoH_ZpBaryonic_MZp-200_MChi-150',
                                  'monoH_ZpBaryonic_MZp-200_MChi-1',
                                  'monoH_ZpBaryonic_MZp-200_MChi-50',
                                  'monoH_ZpBaryonic_MZp-20_MChi-1',
                                  'monoH_ZpBaryonic_MZp-295_MChi-150',
                                  'monoH_ZpBaryonic_MZp-300_MChi-1',
                                  'monoH_ZpBaryonic_MZp-300_MChi-50',
                                  'monoH_ZpBaryonic_MZp-500_MChi-150',
                                  'monoH_ZpBaryonic_MZp-500_MChi-1',
                                  'monoH_ZpBaryonic_MZp-500_MChi-500',
                                  'monoH_ZpBaryonic_MZp-50_MChi-10',
                                  'monoH_ZpBaryonic_MZp-50_MChi-1',
                                  'monoH_ZpBaryonic_MZp-95_MChi-50',
                                  'monoH_ZpBaryonic_MZp-995_MChi-500',

                                 ] ,
                },


  'LepElepTup':  {
                  'isChain'    : True ,
                  'do4MC'      : True ,
                  'do4Data'    : False,
#                  'subTargets' : ['do_LepElepTup','TrigEff','IdIsoSC','l2kin','l3kin','l4kin'],
                  'subTargets' : ['do_LepElepTup','TrigMakerMC','IdIsoSC_Cut','l2kin','l3kin','l4kin'],
                  'onlySample' : [
                                  # DY 
                                  'DYJetsToLL_M-10to50','DYJetsToLL_M-50','DYJetsToLL_M-10to50ext3','DYJetsToLL_M-50-LO','DYJetsToLL_M-50-LO-ext1',
                                  'DY2JetsToLL', 'DY3JetsToLL', 'DY4JetsToLL',
                                  'DYJetsToLL_M-50_HT-100to200_MLM' , 'DYJetsToLL_M-50_HT-100to200_MLM_ext' , 
                                  'DYJetsToLL_M-50_HT-200to400_MLM' , 
                                  'DYJetsToLL_M-50_HT-400to600_MLM' ,
                                  'DYJetsToLL_M-50_HT-600toInf_MLM' ,
                                  'DYJetsToLL_M-50_HT-100to200_ext1','DYJetsToLL_M-50_HT-100to200_ext1','DYJetsToLL_M-50_HT-200to400_ext1',
                                  'DYJetsToLL_M-50_HT-200to400_ext1','DYJetsToLL_M-50_HT-200to400_ext1','DYJetsToLL_M-50_HT-600toInf_ext1',
                                  'DYJetsToLL_M-50_HT-600toInf', 'DYJetsToTT_MuEle_M-50',
                                  # Top
                                  'TTTo2L2Nu','TT','TTTo2L2Nu_ext1',
                                  'ST_t-channel_antitop','ST_t-channel_top',
                                  'ST_tW_antitop','ST_tW_top',
                                  # VV (including WW) 
                                  'WWTo2L2Nu','GluGluWWTo2L2Nu_MCFM','GluGluWWTo2L2NuHiggs_MCFM','GluGluHToWWTo2L2NuPowheg_M125',
                                  'WZTo3LNu','WZ',
                                  'ZZ','Zg',
                                  'ZZTo2L2Nu','ZZTo2L2Q','ZZTo4L',
                                  'Wg_AMCNLOFXFX','WgStarLNuEE','WgStarLNuMuMu',
                                  'Wg_MADGRAPHMLM','WZTo3LNu_mllmin01',
                                  'ggZZ2e2t', 'ggZZ4t', 'ggZZ4m', 'ggZZ2m2t', 'ggZZ2e2m', 'ggZZ4e',
                                  # VVV
                                  'WZZ','ZZZ','WWZ','WWW',
                                  # Higgs 
                                  'GluGluHToTauTau_M125', 'GluGluHToWWTo2L2Nu_M125','GluGluHToWWTo2L2NuPowheg_M125',
                                  'GluGluHToWWTo2L2Nu_alternative_M125','VBFHToWWTo2L2Nu_alternative_M125',
                                  'HWminusJ_HToTauTau_M125', 'HWminusJ_HToWW_M125',
                                  'HWplusJ_HToTauTau_M125', 'HWplusJ_HToWW_M125',
                                  'HZJ_HToTauTau_M125', 'HZJ_HToWW_M125',
                                  'VBFHToTauTau_M125', 'VBFHToWWTo2L2Nu_M125',
                                  'ggZH_HToWW_M125', # missing ggZHToTauTau
                                  'ttHJetToNonbb_M125','TTWJetsToLNu',
                                  'GluGluHToZZTo4L_M125',
                                  'GluGluZH_HToWWTo2L2Nu_M120_noHLT','GluGluZH_HToWWTo2L2Nu_M125_noHLT','GluGluZH_HToWWTo2L2Nu_M130_noHLT',
                                  'HZJ_HToWWTo2L2Nu_M120_noHLT','HZJ_HToWWTo2L2Nu_M125_noHLT','HZJ_HToWWTo2L2Nu_M130_noHLT', 
                                  'HWplusJ_WToLNu_HToWWTo2L2Nu_M125','HWminusJ_WToLNu_HToWWTo2L2Nu_M125',
                                  # ... ggH High Mass
                                  'GluGluHToWWTo2L2Nu_M130',
                                  'GluGluHToWWTo2L2Nu_M135',
                                  'GluGluHToWWTo2L2Nu_M140',
                                  'GluGluHToWWTo2L2Nu_M145',
                                  'GluGluHToWWTo2L2Nu_M150',
                                  'GluGluHToWWTo2L2Nu_M155',
                                  'GluGluHToWWTo2L2Nu_M160',
                                  'GluGluHToWWTo2L2Nu_M165',
                                  'GluGluHToWWTo2L2Nu_M170',
                                  'GluGluHToWWTo2L2Nu_M175',
                                  'GluGluHToWWTo2L2Nu_M180',
                                  'GluGluHToWWTo2L2Nu_M190',
                                  'GluGluHToWWTo2L2Nu_M200',
                                  'GluGluHToWWTo2L2Nu_M210',
                                  'GluGluHToWWTo2L2Nu_M230',
                                  'GluGluHToWWTo2L2Nu_M250',
                                  'GluGluHToWWTo2L2Nu_M270',
                                  'GluGluHToWWTo2L2Nu_M300',
                                  'GluGluHToWWTo2L2Nu_M350',
                                  'GluGluHToWWTo2L2Nu_M400',
                                  'GluGluHToWWTo2L2Nu_M450',
                                  'GluGluHToWWTo2L2Nu_M500',
                                  'GluGluHToWWTo2L2Nu_M550',
                                  'GluGluHToWWTo2L2Nu_M600',
                                  'GluGluHToWWTo2L2Nu_M650',
                                  'GluGluHToWWTo2L2Nu_M700',
                                  'GluGluHToWWTo2L2Nu_M750',
                                  'GluGluHToWWTo2L2Nu_M750_NWA',
                                  'GluGluHToWWTo2L2Nu_M800',
                                  'GluGluHToWWTo2L2Nu_M900',
                                  'GluGluHToWWTo2L2Nu_M1000',
                                  # ... VBF High Mass
                                  'VBFHToWWTo2L2Nu_M130',
                                  'VBFHToWWTo2L2Nu_M135',
                                  'VBFHToWWTo2L2Nu_M140',
                                  'VBFHToWWTo2L2Nu_M145',
                                  'VBFHToWWTo2L2Nu_M150',
                                  'VBFHToWWTo2L2Nu_M155',
                                  'VBFHToWWTo2L2Nu_M160',
                                  'VBFHToWWTo2L2Nu_M165',
                                  'VBFHToWWTo2L2Nu_M170',
                                  'VBFHToWWTo2L2Nu_M175',
                                  'VBFHToWWTo2L2Nu_M180',
                                  'VBFHToWWTo2L2Nu_M190',
                                  'VBFHToWWTo2L2Nu_M200',
                                  'VBFHToWWTo2L2Nu_M210',
                                  'VBFHToWWTo2L2Nu_M230',
                                  'VBFHToWWTo2L2Nu_M250',
                                  'VBFHToWWTo2L2Nu_M270',
                                  'VBFHToWWTo2L2Nu_M300',
                                  'VBFHToWWTo2L2Nu_M350',
                                  'VBFHToWWTo2L2Nu_M400',
                                  'VBFHToWWTo2L2Nu_M450',
                                  'VBFHToWWTo2L2Nu_M500',
                                  'VBFHToWWTo2L2Nu_M550',
                                  'VBFHToWWTo2L2Nu_M600',
                                  'VBFHToWWTo2L2Nu_M650',
                                  'VBFHToWWTo2L2Nu_M700',
                                  'VBFHToWWTo2L2Nu_M750',
                                  'VBFHToWWTo2L2Nu_M750_NWA',
                                  'VBFHToWWTo2L2Nu_M800',
                                  'VBFHToWWTo2L2Nu_M900',
                                  'VBFHToWWTo2L2Nu_M1000',
                                  # VBF 
                                  'WpWmJJ_EWK_QCD_noTop','WpWmJJ_QCD_noTop', 'WpWmJJ_EWK_noTop',
                                  'WGJJ', 'WLLJJToLNu_M-4to60_EWK_QCD', 'WLLJJToLNu_M-60_EWK_QCD', 'WpWpJJ_EWK_QCD',
                                  'WW_DoubleScattering','WWTo2L2Nu_DoubleScattering', 'WmWmJJ_EWK_powheg', 'WpWpJJ_EWK', 'WpWpJJ_QCD', 

                                  # stop

                                  'latino_T2bW',
                                  'latino_T2tt_dM10to80',
                                  'latino_T2tt_dM10to80_genHT160_genMET80',
                                  'latino_T2tt_mStop-150to250',
                                  'latino_T2tt_mStop-250to350',
                                  'latino_T2tt_mStop-350to400',
                                  'latino_T2tt_mStop-400to1200',
                                  'latino_T2tt_mStop425_mLSP325',
                                  'latino_T2tt_mStop500_mLSP325',
                                  'latino_T2tt_mStop850_mLSP100',


                                  # monoH
                                  'monoH_2HDM_MZp-1000_MA0-400',
                                  'monoH_2HDM_MZp-1000_MA0-500',
                                  'monoH_2HDM_MZp-1000_MA0-600',
                                  'monoH_2HDM_MZp-1000_MA0-700',
                                  'monoH_2HDM_MZp-1000_MA0-800',
                                  'monoH_2HDM_MZp-1200_MA0-400',
                                  'monoH_2HDM_MZp-1200_MA0-500',
                                  'monoH_2HDM_MZp-1200_MA0-600',
                                  'monoH_2HDM_MZp-1200_MA0-700',
                                  'monoH_2HDM_MZp-1200_MA0-800',
                                  'monoH_2HDM_MZp-1400_MA0-400',
                                  'monoH_2HDM_MZp-1400_MA0-500',
                                  'monoH_2HDM_MZp-1400_MA0-600',
                                  'monoH_2HDM_MZp-1400_MA0-700',
                                  'monoH_2HDM_MZp-1400_MA0-800',
                                  'monoH_2HDM_MZp-1700_MA0-400',
                                  'monoH_2HDM_MZp-1700_MA0-500',
                                  'monoH_2HDM_MZp-1700_MA0-600',
                                  'monoH_2HDM_MZp-1700_MA0-700',
                                  'monoH_2HDM_MZp-1700_MA0-800',
                                  'monoH_2HDM_MZp-2000_MA0-400',
                                  'monoH_2HDM_MZp-2000_MA0-500',
                                  'monoH_2HDM_MZp-2000_MA0-600',
                                  'monoH_2HDM_MZp-2000_MA0-700',
                                  'monoH_2HDM_MZp-2000_MA0-800',
                                  'monoH_2HDM_MZp-2500_MA0-400',
                                  'monoH_2HDM_MZp-2500_MA0-500',
                                  'monoH_2HDM_MZp-2500_MA0-600',
                                  'monoH_2HDM_MZp-2500_MA0-700',
                                  'monoH_2HDM_MZp-2500_MA0-800',
                                  'monoH_2HDM_MZp-600_MA0-400',
                                  'monoH_2HDM_MZp-800_MA0-400',
                                  'monoH_2HDM_MZp-800_MA0-500',
                                  'monoH_2HDM_MZp-800_MA0-600',
                                  'monoH_2HDM_MZp-1000_MA0-300',
                                  'monoH_2HDM_MZp-1200_MA0-300',
                                  'monoH_2HDM_MZp-1400_MA0-300',
                                  'monoH_2HDM_MZp-1700_MA0-300',
                                  'monoH_2HDM_MZp-2000_MA0-300',
                                  'monoH_2HDM_MZp-2500_MA0-300',
                                  'monoH_2HDM_MZp-600_MA0-300',
                                  'monoH_2HDM_MZp-800_MA0-300',

                                                                'monoH_ZpBaryonic_MZp-10000_MChi-1000',
                                  'monoH_ZpBaryonic_MZp-10000_MChi-10',
                                  'monoH_ZpBaryonic_MZp-10000_MChi-150',
                                  'monoH_ZpBaryonic_MZp-10000_MChi-1',
                                  'monoH_ZpBaryonic_MZp-10000_MChi-50',
                                  'monoH_ZpBaryonic_MZp-1000_MChi-1000',
                                  'monoH_ZpBaryonic_MZp-1000_MChi-1',
                                  'monoH_ZpBaryonic_MZp-100_MChi-10',
                                  'monoH_ZpBaryonic_MZp-100_MChi-1',
                                  'monoH_ZpBaryonic_MZp-10_MChi-1000',
                                  'monoH_ZpBaryonic_MZp-10_MChi-10',
                                  'monoH_ZpBaryonic_MZp-10_MChi-150',
                                  'monoH_ZpBaryonic_MZp-10_MChi-1',
                                  'monoH_ZpBaryonic_MZp-10_MChi-500',
                                  'monoH_ZpBaryonic_MZp-10_MChi-50',
                                  'monoH_ZpBaryonic_MZp-15_MChi-10',
                                  'monoH_ZpBaryonic_MZp-2000_MChi-1',
                                  'monoH_ZpBaryonic_MZp-2000_MChi-500',
                                  'monoH_ZpBaryonic_MZp-200_MChi-150',
                                  'monoH_ZpBaryonic_MZp-200_MChi-1',
                                  'monoH_ZpBaryonic_MZp-200_MChi-50',
                                  'monoH_ZpBaryonic_MZp-20_MChi-1',
                                  'monoH_ZpBaryonic_MZp-295_MChi-150',
                                  'monoH_ZpBaryonic_MZp-300_MChi-1',
                                  'monoH_ZpBaryonic_MZp-300_MChi-50',
                                  'monoH_ZpBaryonic_MZp-500_MChi-150',
                                  'monoH_ZpBaryonic_MZp-500_MChi-1',
                                  'monoH_ZpBaryonic_MZp-500_MChi-500',
                                  'monoH_ZpBaryonic_MZp-50_MChi-10',
                                  'monoH_ZpBaryonic_MZp-50_MChi-1',
                                  'monoH_ZpBaryonic_MZp-95_MChi-50',
                                  'monoH_ZpBaryonic_MZp-995_MChi-500',
 

                                 ] ,
              },

  'LepElepTdo':  {
                  'isChain'    : True ,
                  'do4MC'      : True ,
                  'do4Data'    : False,
                  #'subTargets' : ['do_LepElepTdo','TrigEff','IdIsoSC','l2kin','l3kin','l4kin'],
                  'subTargets' : ['do_LepElepTdo','TrigMakerMC','IdIsoSC_Cut','l2kin','l3kin','l4kin'],
                  'onlySample' : [
                                  # DY 
                                  'DYJetsToLL_M-10to50','DYJetsToLL_M-50','DYJetsToLL_M-10to50ext3','DYJetsToLL_M-50-LO','DYJetsToLL_M-50-LO-ext1',
                                  'DY2JetsToLL', 'DY3JetsToLL', 'DY4JetsToLL',
                                  'DYJetsToLL_M-50_HT-100to200_MLM' , 'DYJetsToLL_M-50_HT-100to200_MLM_ext' , 
                                  'DYJetsToLL_M-50_HT-200to400_MLM' , 
                                  'DYJetsToLL_M-50_HT-400to600_MLM' ,
                                  'DYJetsToLL_M-50_HT-600toInf_MLM' ,
                                  'DYJetsToLL_M-50_HT-100to200_ext1','DYJetsToLL_M-50_HT-100to200_ext1','DYJetsToLL_M-50_HT-200to400_ext1',
                                  'DYJetsToLL_M-50_HT-200to400_ext1','DYJetsToLL_M-50_HT-200to400_ext1','DYJetsToLL_M-50_HT-600toInf_ext1',
                                  'DYJetsToLL_M-50_HT-600toInf', 'DYJetsToTT_MuEle_M-50',
                                  # Top
                                  'TTTo2L2Nu','TT','TTTo2L2Nu_ext1',
                                  'ST_t-channel_antitop','ST_t-channel_top',
                                  'ST_tW_antitop','ST_tW_top',
                                  # VV (including WW) 
                                  'WWTo2L2Nu','GluGluWWTo2L2Nu_MCFM','GluGluWWTo2L2NuHiggs_MCFM',
                                  'WZTo3LNu','WZ',
                                  'ZZ','Zg',
                                  'ZZTo2L2Nu','ZZTo2L2Q','ZZTo4L',
                                  'Wg_AMCNLOFXFX','WgStarLNuEE','WgStarLNuMuMu',
                                  'Wg_MADGRAPHMLM','WZTo3LNu_mllmin01',
                                  'ggZZ2e2t', 'ggZZ4t', 'ggZZ4m', 'ggZZ2m2t', 'ggZZ2e2m', 'ggZZ4e',
                                  # VVV
                                  'WZZ','ZZZ','WWZ','WWW',
                                  # Higgs 
                                  'GluGluHToTauTau_M125', 'GluGluHToWWTo2L2Nu_M125','GluGluHToWWTo2L2NuPowheg_M125',
                                  'GluGluHToWWTo2L2Nu_alternative_M125','VBFHToWWTo2L2Nu_alternative_M125',
                                  'HWminusJ_HToTauTau_M125', 'HWminusJ_HToWW_M125',
                                  'HWplusJ_HToTauTau_M125', 'HWplusJ_HToWW_M125',
                                  'HZJ_HToTauTau_M125', 'HZJ_HToWW_M125',
                                  'VBFHToTauTau_M125', 'VBFHToWWTo2L2Nu_M125',
                                  'ggZH_HToWW_M125', # missing ggZHToTauTau
                                  'ttHJetToNonbb_M125','TTWJetsToLNu',
                                  'GluGluHToZZTo4L_M125',
                                  'GluGluZH_HToWWTo2L2Nu_M120_noHLT','GluGluZH_HToWWTo2L2Nu_M125_noHLT','GluGluZH_HToWWTo2L2Nu_M130_noHLT',
                                  'HZJ_HToWWTo2L2Nu_M120_noHLT','HZJ_HToWWTo2L2Nu_M125_noHLT','HZJ_HToWWTo2L2Nu_M130_noHLT', 
                                  'HWplusJ_WToLNu_HToWWTo2L2Nu_M125','HWminusJ_WToLNu_HToWWTo2L2Nu_M125',

                                  # ... ggH High Mass
                                  'GluGluHToWWTo2L2Nu_M130',
                                  'GluGluHToWWTo2L2Nu_M135',
                                  'GluGluHToWWTo2L2Nu_M140',
                                  'GluGluHToWWTo2L2Nu_M145',
                                  'GluGluHToWWTo2L2Nu_M150',
                                  'GluGluHToWWTo2L2Nu_M155',
                                  'GluGluHToWWTo2L2Nu_M160',
                                  'GluGluHToWWTo2L2Nu_M165',
                                  'GluGluHToWWTo2L2Nu_M170',
                                  'GluGluHToWWTo2L2Nu_M175',
                                  'GluGluHToWWTo2L2Nu_M180',
                                  'GluGluHToWWTo2L2Nu_M190',
                                  'GluGluHToWWTo2L2Nu_M200',
                                  'GluGluHToWWTo2L2Nu_M210',
                                  'GluGluHToWWTo2L2Nu_M230',
                                  'GluGluHToWWTo2L2Nu_M250',
                                  'GluGluHToWWTo2L2Nu_M270',
                                  'GluGluHToWWTo2L2Nu_M300',
                                  'GluGluHToWWTo2L2Nu_M350',
                                  'GluGluHToWWTo2L2Nu_M400',
                                  'GluGluHToWWTo2L2Nu_M450',
                                  'GluGluHToWWTo2L2Nu_M500',
                                  'GluGluHToWWTo2L2Nu_M550',
                                  'GluGluHToWWTo2L2Nu_M600',
                                  'GluGluHToWWTo2L2Nu_M650',
                                  'GluGluHToWWTo2L2Nu_M700',
                                  'GluGluHToWWTo2L2Nu_M750',
                                  'GluGluHToWWTo2L2Nu_M750_NWA',
                                  'GluGluHToWWTo2L2Nu_M800',
                                  'GluGluHToWWTo2L2Nu_M900',
                                  'GluGluHToWWTo2L2Nu_M1000',
                                  # ... VBF High Mass
                                  'VBFHToWWTo2L2Nu_M130',
                                  'VBFHToWWTo2L2Nu_M135',
                                  'VBFHToWWTo2L2Nu_M140',
                                  'VBFHToWWTo2L2Nu_M145',
                                  'VBFHToWWTo2L2Nu_M150',
                                  'VBFHToWWTo2L2Nu_M155',
                                  'VBFHToWWTo2L2Nu_M160',
                                  'VBFHToWWTo2L2Nu_M165',
                                  'VBFHToWWTo2L2Nu_M170',
                                  'VBFHToWWTo2L2Nu_M175',
                                  'VBFHToWWTo2L2Nu_M180',
                                  'VBFHToWWTo2L2Nu_M190',
                                  'VBFHToWWTo2L2Nu_M200',
                                  'VBFHToWWTo2L2Nu_M210',
                                  'VBFHToWWTo2L2Nu_M230',
                                  'VBFHToWWTo2L2Nu_M250',
                                  'VBFHToWWTo2L2Nu_M270',
                                  'VBFHToWWTo2L2Nu_M300',
                                  'VBFHToWWTo2L2Nu_M350',
                                  'VBFHToWWTo2L2Nu_M400',
                                  'VBFHToWWTo2L2Nu_M450',
                                  'VBFHToWWTo2L2Nu_M500',
                                  'VBFHToWWTo2L2Nu_M550',
                                  'VBFHToWWTo2L2Nu_M600',
                                  'VBFHToWWTo2L2Nu_M650',
                                  'VBFHToWWTo2L2Nu_M700',
                                  'VBFHToWWTo2L2Nu_M750',
                                  'VBFHToWWTo2L2Nu_M750_NWA',
                                  'VBFHToWWTo2L2Nu_M800',
                                  'VBFHToWWTo2L2Nu_M900',
                                  'VBFHToWWTo2L2Nu_M1000',
                                  # VBF 
                                  'WpWmJJ_EWK_QCD_noTop','WpWmJJ_QCD_noTop', 'WpWmJJ_EWK_noTop',
                                  'WGJJ', 'WLLJJToLNu_M-4to60_EWK_QCD', 'WLLJJToLNu_M-60_EWK_QCD', 'WpWpJJ_EWK_QCD',
                                  'WW_DoubleScattering','WWTo2L2Nu_DoubleScattering', 'WmWmJJ_EWK_powheg', 'WpWpJJ_EWK', 'WpWpJJ_QCD', 


                                  # stop

                                  'latino_T2bW',
                                  'latino_T2tt_dM10to80',
                                  'latino_T2tt_dM10to80_genHT160_genMET80',
                                  'latino_T2tt_mStop-150to250',
                                  'latino_T2tt_mStop-250to350',
                                  'latino_T2tt_mStop-350to400',
                                  'latino_T2tt_mStop-400to1200',
                                  'latino_T2tt_mStop425_mLSP325',
                                  'latino_T2tt_mStop500_mLSP325',
                                  'latino_T2tt_mStop850_mLSP100',


                                  # monoH
                                  'monoH_2HDM_MZp-1000_MA0-400',
                                  'monoH_2HDM_MZp-1000_MA0-500',
                                  'monoH_2HDM_MZp-1000_MA0-600',
                                  'monoH_2HDM_MZp-1000_MA0-700',
                                  'monoH_2HDM_MZp-1000_MA0-800',
                                  'monoH_2HDM_MZp-1200_MA0-400',
                                  'monoH_2HDM_MZp-1200_MA0-500',
                                  'monoH_2HDM_MZp-1200_MA0-600',
                                  'monoH_2HDM_MZp-1200_MA0-700',
                                  'monoH_2HDM_MZp-1200_MA0-800',
                                  'monoH_2HDM_MZp-1400_MA0-400',
                                  'monoH_2HDM_MZp-1400_MA0-500',
                                  'monoH_2HDM_MZp-1400_MA0-600',
                                  'monoH_2HDM_MZp-1400_MA0-700',
                                  'monoH_2HDM_MZp-1400_MA0-800',
                                  'monoH_2HDM_MZp-1700_MA0-400',
                                  'monoH_2HDM_MZp-1700_MA0-500',
                                  'monoH_2HDM_MZp-1700_MA0-600',
                                  'monoH_2HDM_MZp-1700_MA0-700',
                                  'monoH_2HDM_MZp-1700_MA0-800',
                                  'monoH_2HDM_MZp-2000_MA0-400',
                                  'monoH_2HDM_MZp-2000_MA0-500',
                                  'monoH_2HDM_MZp-2000_MA0-600',
                                  'monoH_2HDM_MZp-2000_MA0-700',
                                  'monoH_2HDM_MZp-2000_MA0-800',
                                  'monoH_2HDM_MZp-2500_MA0-400',
                                  'monoH_2HDM_MZp-2500_MA0-500',
                                  'monoH_2HDM_MZp-2500_MA0-600',
                                  'monoH_2HDM_MZp-2500_MA0-700',
                                  'monoH_2HDM_MZp-2500_MA0-800',
                                  'monoH_2HDM_MZp-600_MA0-400',
                                  'monoH_2HDM_MZp-800_MA0-400',
                                  'monoH_2HDM_MZp-800_MA0-500',
                                  'monoH_2HDM_MZp-800_MA0-600',
                                  'monoH_2HDM_MZp-1000_MA0-300',
                                  'monoH_2HDM_MZp-1200_MA0-300',
                                  'monoH_2HDM_MZp-1400_MA0-300',
                                  'monoH_2HDM_MZp-1700_MA0-300',
                                  'monoH_2HDM_MZp-2000_MA0-300',
                                  'monoH_2HDM_MZp-2500_MA0-300',
                                  'monoH_2HDM_MZp-600_MA0-300',
                                  'monoH_2HDM_MZp-800_MA0-300',

                                                                'monoH_ZpBaryonic_MZp-10000_MChi-1000',
                                  'monoH_ZpBaryonic_MZp-10000_MChi-10',
                                  'monoH_ZpBaryonic_MZp-10000_MChi-150',
                                  'monoH_ZpBaryonic_MZp-10000_MChi-1',
                                  'monoH_ZpBaryonic_MZp-10000_MChi-50',
                                  'monoH_ZpBaryonic_MZp-1000_MChi-1000',
                                  'monoH_ZpBaryonic_MZp-1000_MChi-1',
                                  'monoH_ZpBaryonic_MZp-100_MChi-10',
                                  'monoH_ZpBaryonic_MZp-100_MChi-1',
                                  'monoH_ZpBaryonic_MZp-10_MChi-1000',
                                  'monoH_ZpBaryonic_MZp-10_MChi-10',
                                  'monoH_ZpBaryonic_MZp-10_MChi-150',
                                  'monoH_ZpBaryonic_MZp-10_MChi-1',
                                  'monoH_ZpBaryonic_MZp-10_MChi-500',
                                  'monoH_ZpBaryonic_MZp-10_MChi-50',
                                  'monoH_ZpBaryonic_MZp-15_MChi-10',
                                  'monoH_ZpBaryonic_MZp-2000_MChi-1',
                                  'monoH_ZpBaryonic_MZp-2000_MChi-500',
                                  'monoH_ZpBaryonic_MZp-200_MChi-150',
                                  'monoH_ZpBaryonic_MZp-200_MChi-1',
                                  'monoH_ZpBaryonic_MZp-200_MChi-50',
                                  'monoH_ZpBaryonic_MZp-20_MChi-1',
                                  'monoH_ZpBaryonic_MZp-295_MChi-150',
                                  'monoH_ZpBaryonic_MZp-300_MChi-1',
                                  'monoH_ZpBaryonic_MZp-300_MChi-50',
                                  'monoH_ZpBaryonic_MZp-500_MChi-150',
                                  'monoH_ZpBaryonic_MZp-500_MChi-1',
                                  'monoH_ZpBaryonic_MZp-500_MChi-500',
                                  'monoH_ZpBaryonic_MZp-50_MChi-10',
                                  'monoH_ZpBaryonic_MZp-50_MChi-1',
                                  'monoH_ZpBaryonic_MZp-95_MChi-50',
                                  'monoH_ZpBaryonic_MZp-995_MChi-500',
 

                                 ] ,
              },


  'LepMupTup':  {
                  'isChain'    : True ,
                  'do4MC'      : True ,
                  'do4Data'    : False,
                  #'subTargets' : ['do_LepMupTup','TrigEff','IdIsoSC','l2kin','l3kin','l4kin'],
                  'subTargets' : ['do_LepMupTup','TrigMakerMC','IdIsoSC_Cut','l2kin','l3kin','l4kin'],
                  'onlySample' : [
                                  # DY 
                                  'DYJetsToLL_M-10to50','DYJetsToLL_M-50','DYJetsToLL_M-10to50ext3','DYJetsToLL_M-50-LO','DYJetsToLL_M-50-LO-ext1',
                                  'DY2JetsToLL', 'DY3JetsToLL', 'DY4JetsToLL',
                                  'DYJetsToLL_M-50_HT-100to200_MLM' , 'DYJetsToLL_M-50_HT-100to200_MLM_ext' , 
                                  'DYJetsToLL_M-50_HT-200to400_MLM' , 
                                  'DYJetsToLL_M-50_HT-400to600_MLM' ,
                                  'DYJetsToLL_M-50_HT-600toInf_MLM' ,
                                  'DYJetsToLL_M-50_HT-100to200_ext1','DYJetsToLL_M-50_HT-100to200_ext1','DYJetsToLL_M-50_HT-200to400_ext1',
                                  'DYJetsToLL_M-50_HT-200to400_ext1','DYJetsToLL_M-50_HT-200to400_ext1','DYJetsToLL_M-50_HT-600toInf_ext1',
                                  'DYJetsToLL_M-50_HT-600toInf', 'DYJetsToTT_MuEle_M-50',
                                  # Top
                                  'TTTo2L2Nu','TT','TTTo2L2Nu_ext1',
                                  'ST_t-channel_antitop','ST_t-channel_top',
                                  'ST_tW_antitop','ST_tW_top',
                                  # VV (including WW) 
                                  'WWTo2L2Nu','GluGluWWTo2L2Nu_MCFM','GluGluWWTo2L2NuHiggs_MCFM',
                                  'WZTo3LNu','WZ',
                                  'ZZ','Zg',
                                  'ZZTo2L2Nu','ZZTo2L2Q','ZZTo4L',
                                  'Wg_AMCNLOFXFX','WgStarLNuEE','WgStarLNuMuMu',
                                  'Wg_MADGRAPHMLM','WZTo3LNu_mllmin01',
                                  'ggZZ2e2t', 'ggZZ4t', 'ggZZ4m', 'ggZZ2m2t', 'ggZZ2e2m', 'ggZZ4e',
                                  # VVV
                                  'WZZ','ZZZ','WWZ','WWW',
                                  # Higgs 
                                  'GluGluHToTauTau_M125', 'GluGluHToWWTo2L2Nu_M125','GluGluHToWWTo2L2NuPowheg_M125',
                                  'GluGluHToWWTo2L2Nu_alternative_M125','VBFHToWWTo2L2Nu_alternative_M125',
                                  'HWminusJ_HToTauTau_M125', 'HWminusJ_HToWW_M125',
                                  'HWplusJ_HToTauTau_M125', 'HWplusJ_HToWW_M125',
                                  'HZJ_HToTauTau_M125', 'HZJ_HToWW_M125',
                                  'VBFHToTauTau_M125', 'VBFHToWWTo2L2Nu_M125',
                                  'ggZH_HToWW_M125', # missing ggZHToTauTau
                                  'ttHJetToNonbb_M125','TTWJetsToLNu',
                                  'GluGluHToZZTo4L_M125',
                                  'GluGluZH_HToWWTo2L2Nu_M120_noHLT','GluGluZH_HToWWTo2L2Nu_M125_noHLT','GluGluZH_HToWWTo2L2Nu_M130_noHLT',
                                  'HZJ_HToWWTo2L2Nu_M120_noHLT','HZJ_HToWWTo2L2Nu_M125_noHLT','HZJ_HToWWTo2L2Nu_M130_noHLT', 
                                  'HWplusJ_WToLNu_HToWWTo2L2Nu_M125','HWminusJ_WToLNu_HToWWTo2L2Nu_M125',

                                  # ... ggH High Mass
                                  'GluGluHToWWTo2L2Nu_M130',
                                  'GluGluHToWWTo2L2Nu_M135',
                                  'GluGluHToWWTo2L2Nu_M140',
                                  'GluGluHToWWTo2L2Nu_M145',
                                  'GluGluHToWWTo2L2Nu_M150',
                                  'GluGluHToWWTo2L2Nu_M155',
                                  'GluGluHToWWTo2L2Nu_M160',
                                  'GluGluHToWWTo2L2Nu_M165',
                                  'GluGluHToWWTo2L2Nu_M170',
                                  'GluGluHToWWTo2L2Nu_M175',
                                  'GluGluHToWWTo2L2Nu_M180',
                                  'GluGluHToWWTo2L2Nu_M190',
                                  'GluGluHToWWTo2L2Nu_M200',
                                  'GluGluHToWWTo2L2Nu_M210',
                                  'GluGluHToWWTo2L2Nu_M230',
                                  'GluGluHToWWTo2L2Nu_M250',
                                  'GluGluHToWWTo2L2Nu_M270',
                                  'GluGluHToWWTo2L2Nu_M300',
                                  'GluGluHToWWTo2L2Nu_M350',
                                  'GluGluHToWWTo2L2Nu_M400',
                                  'GluGluHToWWTo2L2Nu_M450',
                                  'GluGluHToWWTo2L2Nu_M500',
                                  'GluGluHToWWTo2L2Nu_M550',
                                  'GluGluHToWWTo2L2Nu_M600',
                                  'GluGluHToWWTo2L2Nu_M650',
                                  'GluGluHToWWTo2L2Nu_M700',
                                  'GluGluHToWWTo2L2Nu_M750',
                                  'GluGluHToWWTo2L2Nu_M750_NWA',
                                  'GluGluHToWWTo2L2Nu_M800',
                                  'GluGluHToWWTo2L2Nu_M900',
                                  'GluGluHToWWTo2L2Nu_M1000',
                                  # ... VBF High Mass
                                  'VBFHToWWTo2L2Nu_M130',
                                  'VBFHToWWTo2L2Nu_M135',
                                  'VBFHToWWTo2L2Nu_M140',
                                  'VBFHToWWTo2L2Nu_M145',
                                  'VBFHToWWTo2L2Nu_M150',
                                  'VBFHToWWTo2L2Nu_M155',
                                  'VBFHToWWTo2L2Nu_M160',
                                  'VBFHToWWTo2L2Nu_M165',
                                  'VBFHToWWTo2L2Nu_M170',
                                  'VBFHToWWTo2L2Nu_M175',
                                  'VBFHToWWTo2L2Nu_M180',
                                  'VBFHToWWTo2L2Nu_M190',
                                  'VBFHToWWTo2L2Nu_M200',
                                  'VBFHToWWTo2L2Nu_M210',
                                  'VBFHToWWTo2L2Nu_M230',
                                  'VBFHToWWTo2L2Nu_M250',
                                  'VBFHToWWTo2L2Nu_M270',
                                  'VBFHToWWTo2L2Nu_M300',
                                  'VBFHToWWTo2L2Nu_M350',
                                  'VBFHToWWTo2L2Nu_M400',
                                  'VBFHToWWTo2L2Nu_M450',
                                  'VBFHToWWTo2L2Nu_M500',
                                  'VBFHToWWTo2L2Nu_M550',
                                  'VBFHToWWTo2L2Nu_M600',
                                  'VBFHToWWTo2L2Nu_M650',
                                  'VBFHToWWTo2L2Nu_M700',
                                  'VBFHToWWTo2L2Nu_M750',
                                  'VBFHToWWTo2L2Nu_M750_NWA',
                                  'VBFHToWWTo2L2Nu_M800',
                                  'VBFHToWWTo2L2Nu_M900',
                                  'VBFHToWWTo2L2Nu_M1000',
                                  # VBF 
                                  'WpWmJJ_EWK_QCD_noTop','WpWmJJ_QCD_noTop', 'WpWmJJ_EWK_noTop',
                                  'WGJJ', 'WLLJJToLNu_M-4to60_EWK_QCD', 'WLLJJToLNu_M-60_EWK_QCD', 'WpWpJJ_EWK_QCD',
                                  'WW_DoubleScattering','WWTo2L2Nu_DoubleScattering', 'WmWmJJ_EWK_powheg', 'WpWpJJ_EWK', 'WpWpJJ_QCD', 

                                  # stop

                                  'latino_T2bW',
                                  'latino_T2tt_dM10to80',
                                  'latino_T2tt_dM10to80_genHT160_genMET80',
                                  'latino_T2tt_mStop-150to250',
                                  'latino_T2tt_mStop-250to350',
                                  'latino_T2tt_mStop-350to400',
                                  'latino_T2tt_mStop-400to1200',
                                  'latino_T2tt_mStop425_mLSP325',
                                  'latino_T2tt_mStop500_mLSP325',
                                  'latino_T2tt_mStop850_mLSP100',



                                  # monoH
                                  'monoH_2HDM_MZp-1000_MA0-400',
                                  'monoH_2HDM_MZp-1000_MA0-500',
                                  'monoH_2HDM_MZp-1000_MA0-600',
                                  'monoH_2HDM_MZp-1000_MA0-700',
                                  'monoH_2HDM_MZp-1000_MA0-800',
                                  'monoH_2HDM_MZp-1200_MA0-400',
                                  'monoH_2HDM_MZp-1200_MA0-500',
                                  'monoH_2HDM_MZp-1200_MA0-600',
                                  'monoH_2HDM_MZp-1200_MA0-700',
                                  'monoH_2HDM_MZp-1200_MA0-800',
                                  'monoH_2HDM_MZp-1400_MA0-400',
                                  'monoH_2HDM_MZp-1400_MA0-500',
                                  'monoH_2HDM_MZp-1400_MA0-600',
                                  'monoH_2HDM_MZp-1400_MA0-700',
                                  'monoH_2HDM_MZp-1400_MA0-800',
                                  'monoH_2HDM_MZp-1700_MA0-400',
                                  'monoH_2HDM_MZp-1700_MA0-500',
                                  'monoH_2HDM_MZp-1700_MA0-600',
                                  'monoH_2HDM_MZp-1700_MA0-700',
                                  'monoH_2HDM_MZp-1700_MA0-800',
                                  'monoH_2HDM_MZp-2000_MA0-400',
                                  'monoH_2HDM_MZp-2000_MA0-500',
                                  'monoH_2HDM_MZp-2000_MA0-600',
                                  'monoH_2HDM_MZp-2000_MA0-700',
                                  'monoH_2HDM_MZp-2000_MA0-800',
                                  'monoH_2HDM_MZp-2500_MA0-400',
                                  'monoH_2HDM_MZp-2500_MA0-500',
                                  'monoH_2HDM_MZp-2500_MA0-600',
                                  'monoH_2HDM_MZp-2500_MA0-700',
                                  'monoH_2HDM_MZp-2500_MA0-800',
                                  'monoH_2HDM_MZp-600_MA0-400',
                                  'monoH_2HDM_MZp-800_MA0-400',
                                  'monoH_2HDM_MZp-800_MA0-500',
                                  'monoH_2HDM_MZp-800_MA0-600',
                                  'monoH_2HDM_MZp-1000_MA0-300',
                                  'monoH_2HDM_MZp-1200_MA0-300',
                                  'monoH_2HDM_MZp-1400_MA0-300',
                                  'monoH_2HDM_MZp-1700_MA0-300',
                                  'monoH_2HDM_MZp-2000_MA0-300',
                                  'monoH_2HDM_MZp-2500_MA0-300',
                                  'monoH_2HDM_MZp-600_MA0-300',
                                  'monoH_2HDM_MZp-800_MA0-300',

                                                                'monoH_ZpBaryonic_MZp-10000_MChi-1000',
                                  'monoH_ZpBaryonic_MZp-10000_MChi-10',
                                  'monoH_ZpBaryonic_MZp-10000_MChi-150',
                                  'monoH_ZpBaryonic_MZp-10000_MChi-1',
                                  'monoH_ZpBaryonic_MZp-10000_MChi-50',
                                  'monoH_ZpBaryonic_MZp-1000_MChi-1000',
                                  'monoH_ZpBaryonic_MZp-1000_MChi-1',
                                  'monoH_ZpBaryonic_MZp-100_MChi-10',
                                  'monoH_ZpBaryonic_MZp-100_MChi-1',
                                  'monoH_ZpBaryonic_MZp-10_MChi-1000',
                                  'monoH_ZpBaryonic_MZp-10_MChi-10',
                                  'monoH_ZpBaryonic_MZp-10_MChi-150',
                                  'monoH_ZpBaryonic_MZp-10_MChi-1',
                                  'monoH_ZpBaryonic_MZp-10_MChi-500',
                                  'monoH_ZpBaryonic_MZp-10_MChi-50',
                                  'monoH_ZpBaryonic_MZp-15_MChi-10',
                                  'monoH_ZpBaryonic_MZp-2000_MChi-1',
                                  'monoH_ZpBaryonic_MZp-2000_MChi-500',
                                  'monoH_ZpBaryonic_MZp-200_MChi-150',
                                  'monoH_ZpBaryonic_MZp-200_MChi-1',
                                  'monoH_ZpBaryonic_MZp-200_MChi-50',
                                  'monoH_ZpBaryonic_MZp-20_MChi-1',
                                  'monoH_ZpBaryonic_MZp-295_MChi-150',
                                  'monoH_ZpBaryonic_MZp-300_MChi-1',
                                  'monoH_ZpBaryonic_MZp-300_MChi-50',
                                  'monoH_ZpBaryonic_MZp-500_MChi-150',
                                  'monoH_ZpBaryonic_MZp-500_MChi-1',
                                  'monoH_ZpBaryonic_MZp-500_MChi-500',
                                  'monoH_ZpBaryonic_MZp-50_MChi-10',
                                  'monoH_ZpBaryonic_MZp-50_MChi-1',
                                  'monoH_ZpBaryonic_MZp-95_MChi-50',
                                  'monoH_ZpBaryonic_MZp-995_MChi-500',
 

                                 ] ,
              },

  'LepMupTdo':  {
                  'isChain'    : True ,
                  'do4MC'      : True ,
                  'do4Data'    : False,
                  #'subTargets' : ['do_LepMupTdo','TrigEff','IdIsoSC','l2kin','l3kin','l4kin'],
                  'subTargets' : ['do_LepMupTdo','TrigMakerMC','IdIsoSC_Cut','l2kin','l3kin','l4kin'],
                  'onlySample' : [
                                  # DY 
                                  'DYJetsToLL_M-10to50','DYJetsToLL_M-50','DYJetsToLL_M-10to50ext3','DYJetsToLL_M-50-LO','DYJetsToLL_M-50-LO-ext1',
                                  'DY2JetsToLL', 'DY3JetsToLL', 'DY4JetsToLL',
                                  'DYJetsToLL_M-50_HT-100to200_MLM' , 'DYJetsToLL_M-50_HT-100to200_MLM_ext' , 
                                  'DYJetsToLL_M-50_HT-200to400_MLM' , 
                                  'DYJetsToLL_M-50_HT-400to600_MLM' ,
                                  'DYJetsToLL_M-50_HT-600toInf_MLM' ,
                                  'DYJetsToLL_M-50_HT-100to200_ext1','DYJetsToLL_M-50_HT-100to200_ext1','DYJetsToLL_M-50_HT-200to400_ext1',
                                  'DYJetsToLL_M-50_HT-200to400_ext1','DYJetsToLL_M-50_HT-200to400_ext1','DYJetsToLL_M-50_HT-600toInf_ext1',
                                  'DYJetsToLL_M-50_HT-600toInf', 'DYJetsToTT_MuEle_M-50',
                                  # Top
                                  'TTTo2L2Nu','TT','TTTo2L2Nu_ext1',
                                  'ST_t-channel_antitop','ST_t-channel_top',
                                  'ST_tW_antitop','ST_tW_top',
                                  # VV (including WW) 
                                  'WWTo2L2Nu','GluGluWWTo2L2Nu_MCFM','GluGluWWTo2L2NuHiggs_MCFM',
                                  'GluGluHToWWTo2L2Nu_alternative_M125','VBFHToWWTo2L2Nu_alternative_M125',
                                  'WZTo3LNu','WZ',
                                  'ZZ','Zg',
                                  'ZZTo2L2Nu','ZZTo2L2Q','ZZTo4L',
                                  'Wg_AMCNLOFXFX','WgStarLNuEE','WgStarLNuMuMu',
                                  'Wg_MADGRAPHMLM','WZTo3LNu_mllmin01',
                                  'ggZZ2e2t', 'ggZZ4t', 'ggZZ4m', 'ggZZ2m2t', 'ggZZ2e2m', 'ggZZ4e',
                                  # VVV
                                  'WZZ','ZZZ','WWZ','WWW',
                                  # Higgs 
                                  'GluGluHToTauTau_M125', 'GluGluHToWWTo2L2Nu_M125','GluGluHToWWTo2L2NuPowheg_M125',
                                  'GluGluHToWWTo2L2Nu_alternative_M125','VBFHToWWTo2L2Nu_alternative_M125',
                                  'HWminusJ_HToTauTau_M125', 'HWminusJ_HToWW_M125',
                                  'HWplusJ_HToTauTau_M125', 'HWplusJ_HToWW_M125',
                                  'HZJ_HToTauTau_M125', 'HZJ_HToWW_M125',
                                  'VBFHToTauTau_M125', 'VBFHToWWTo2L2Nu_M125',
                                  'ggZH_HToWW_M125', # missing ggZHToTauTau
                                  'ttHJetToNonbb_M125','TTWJetsToLNu',
                                  'GluGluHToZZTo4L_M125',
                                  'GluGluZH_HToWWTo2L2Nu_M120_noHLT','GluGluZH_HToWWTo2L2Nu_M125_noHLT','GluGluZH_HToWWTo2L2Nu_M130_noHLT',
                                  'HZJ_HToWWTo2L2Nu_M120_noHLT','HZJ_HToWWTo2L2Nu_M125_noHLT','HZJ_HToWWTo2L2Nu_M130_noHLT', 
                                  'HWplusJ_WToLNu_HToWWTo2L2Nu_M125','HWminusJ_WToLNu_HToWWTo2L2Nu_M125',
                                  # ... ggH High Mass
                                  'GluGluHToWWTo2L2Nu_M130',
                                  'GluGluHToWWTo2L2Nu_M135',
                                  'GluGluHToWWTo2L2Nu_M140',
                                  'GluGluHToWWTo2L2Nu_M145',
                                  'GluGluHToWWTo2L2Nu_M150',
                                  'GluGluHToWWTo2L2Nu_M155',
                                  'GluGluHToWWTo2L2Nu_M160',
                                  'GluGluHToWWTo2L2Nu_M165',
                                  'GluGluHToWWTo2L2Nu_M170',
                                  'GluGluHToWWTo2L2Nu_M175',
                                  'GluGluHToWWTo2L2Nu_M180',
                                  'GluGluHToWWTo2L2Nu_M190',
                                  'GluGluHToWWTo2L2Nu_M200',
                                  'GluGluHToWWTo2L2Nu_M210',
                                  'GluGluHToWWTo2L2Nu_M230',
                                  'GluGluHToWWTo2L2Nu_M250',
                                  'GluGluHToWWTo2L2Nu_M270',
                                  'GluGluHToWWTo2L2Nu_M300',
                                  'GluGluHToWWTo2L2Nu_M350',
                                  'GluGluHToWWTo2L2Nu_M400',
                                  'GluGluHToWWTo2L2Nu_M450',
                                  'GluGluHToWWTo2L2Nu_M500',
                                  'GluGluHToWWTo2L2Nu_M550',
                                  'GluGluHToWWTo2L2Nu_M600',
                                  'GluGluHToWWTo2L2Nu_M650',
                                  'GluGluHToWWTo2L2Nu_M700',
                                  'GluGluHToWWTo2L2Nu_M750',
                                  'GluGluHToWWTo2L2Nu_M750_NWA',
                                  'GluGluHToWWTo2L2Nu_M800',
                                  'GluGluHToWWTo2L2Nu_M900',
                                  'GluGluHToWWTo2L2Nu_M1000',
                                  # ... VBF High Mass
                                  'VBFHToWWTo2L2Nu_M130',
                                  'VBFHToWWTo2L2Nu_M135',
                                  'VBFHToWWTo2L2Nu_M140',
                                  'VBFHToWWTo2L2Nu_M145',
                                  'VBFHToWWTo2L2Nu_M150',
                                  'VBFHToWWTo2L2Nu_M155',
                                  'VBFHToWWTo2L2Nu_M160',
                                  'VBFHToWWTo2L2Nu_M165',
                                  'VBFHToWWTo2L2Nu_M170',
                                  'VBFHToWWTo2L2Nu_M175',
                                  'VBFHToWWTo2L2Nu_M180',
                                  'VBFHToWWTo2L2Nu_M190',
                                  'VBFHToWWTo2L2Nu_M200',
                                  'VBFHToWWTo2L2Nu_M210',
                                  'VBFHToWWTo2L2Nu_M230',
                                  'VBFHToWWTo2L2Nu_M250',
                                  'VBFHToWWTo2L2Nu_M270',
                                  'VBFHToWWTo2L2Nu_M300',
                                  'VBFHToWWTo2L2Nu_M350',
                                  'VBFHToWWTo2L2Nu_M400',
                                  'VBFHToWWTo2L2Nu_M450',
                                  'VBFHToWWTo2L2Nu_M500',
                                  'VBFHToWWTo2L2Nu_M550',
                                  'VBFHToWWTo2L2Nu_M600',
                                  'VBFHToWWTo2L2Nu_M650',
                                  'VBFHToWWTo2L2Nu_M700',
                                  'VBFHToWWTo2L2Nu_M750',
                                  'VBFHToWWTo2L2Nu_M750_NWA',
                                  'VBFHToWWTo2L2Nu_M800',
                                  'VBFHToWWTo2L2Nu_M900',
                                  'VBFHToWWTo2L2Nu_M1000',
                                  # VBF 
                                  'WpWmJJ_EWK_QCD_noTop','WpWmJJ_QCD_noTop', 'WpWmJJ_EWK_noTop',
                                  'WGJJ', 'WLLJJToLNu_M-4to60_EWK_QCD', 'WLLJJToLNu_M-60_EWK_QCD', 'WpWpJJ_EWK_QCD',
                                  'WW_DoubleScattering','WWTo2L2Nu_DoubleScattering', 'WmWmJJ_EWK_powheg', 'WpWpJJ_EWK', 'WpWpJJ_QCD', 

                                  # stop

                                  'latino_T2bW',
                                  'latino_T2tt_dM10to80',
                                  'latino_T2tt_dM10to80_genHT160_genMET80',
                                  'latino_T2tt_mStop-150to250',
                                  'latino_T2tt_mStop-250to350',
                                  'latino_T2tt_mStop-350to400',
                                  'latino_T2tt_mStop-400to1200',
                                  'latino_T2tt_mStop425_mLSP325',
                                  'latino_T2tt_mStop500_mLSP325',
                                  'latino_T2tt_mStop850_mLSP100',


                                  # monoH
                                  'monoH_2HDM_MZp-1000_MA0-400',
                                  'monoH_2HDM_MZp-1000_MA0-500',
                                  'monoH_2HDM_MZp-1000_MA0-600',
                                  'monoH_2HDM_MZp-1000_MA0-700',
                                  'monoH_2HDM_MZp-1000_MA0-800',
                                  'monoH_2HDM_MZp-1200_MA0-400',
                                  'monoH_2HDM_MZp-1200_MA0-500',
                                  'monoH_2HDM_MZp-1200_MA0-600',
                                  'monoH_2HDM_MZp-1200_MA0-700',
                                  'monoH_2HDM_MZp-1200_MA0-800',
                                  'monoH_2HDM_MZp-1400_MA0-400',
                                  'monoH_2HDM_MZp-1400_MA0-500',
                                  'monoH_2HDM_MZp-1400_MA0-600',
                                  'monoH_2HDM_MZp-1400_MA0-700',
                                  'monoH_2HDM_MZp-1400_MA0-800',
                                  'monoH_2HDM_MZp-1700_MA0-400',
                                  'monoH_2HDM_MZp-1700_MA0-500',
                                  'monoH_2HDM_MZp-1700_MA0-600',
                                  'monoH_2HDM_MZp-1700_MA0-700',
                                  'monoH_2HDM_MZp-1700_MA0-800',
                                  'monoH_2HDM_MZp-2000_MA0-400',
                                  'monoH_2HDM_MZp-2000_MA0-500',
                                  'monoH_2HDM_MZp-2000_MA0-600',
                                  'monoH_2HDM_MZp-2000_MA0-700',
                                  'monoH_2HDM_MZp-2000_MA0-800',
                                  'monoH_2HDM_MZp-2500_MA0-400',
                                  'monoH_2HDM_MZp-2500_MA0-500',
                                  'monoH_2HDM_MZp-2500_MA0-600',
                                  'monoH_2HDM_MZp-2500_MA0-700',
                                  'monoH_2HDM_MZp-2500_MA0-800',
                                  'monoH_2HDM_MZp-600_MA0-400',
                                  'monoH_2HDM_MZp-800_MA0-400',
                                  'monoH_2HDM_MZp-800_MA0-500',
                                  'monoH_2HDM_MZp-800_MA0-600',
                                  'monoH_2HDM_MZp-1000_MA0-300',
                                  'monoH_2HDM_MZp-1200_MA0-300',
                                  'monoH_2HDM_MZp-1400_MA0-300',
                                  'monoH_2HDM_MZp-1700_MA0-300',
                                  'monoH_2HDM_MZp-2000_MA0-300',
                                  'monoH_2HDM_MZp-2500_MA0-300',
                                  'monoH_2HDM_MZp-600_MA0-300',
                                  'monoH_2HDM_MZp-800_MA0-300',

                               
                                  'monoH_ZpBaryonic_MZp-10000_MChi-1000',
                                  'monoH_ZpBaryonic_MZp-10000_MChi-10',
                                  'monoH_ZpBaryonic_MZp-10000_MChi-150',
                                  'monoH_ZpBaryonic_MZp-10000_MChi-1',
                                  'monoH_ZpBaryonic_MZp-10000_MChi-50',
                                  'monoH_ZpBaryonic_MZp-1000_MChi-1000',
                                  'monoH_ZpBaryonic_MZp-1000_MChi-1',
                                  'monoH_ZpBaryonic_MZp-100_MChi-10',
                                  'monoH_ZpBaryonic_MZp-100_MChi-1',
                                  'monoH_ZpBaryonic_MZp-10_MChi-1000',
                                  'monoH_ZpBaryonic_MZp-10_MChi-10',
                                  'monoH_ZpBaryonic_MZp-10_MChi-150',
                                  'monoH_ZpBaryonic_MZp-10_MChi-1',
                                  'monoH_ZpBaryonic_MZp-10_MChi-500',
                                  'monoH_ZpBaryonic_MZp-10_MChi-50',
                                  'monoH_ZpBaryonic_MZp-15_MChi-10',
                                  'monoH_ZpBaryonic_MZp-2000_MChi-1',
                                  'monoH_ZpBaryonic_MZp-2000_MChi-500',
                                  'monoH_ZpBaryonic_MZp-200_MChi-150',
                                  'monoH_ZpBaryonic_MZp-200_MChi-1',
                                  'monoH_ZpBaryonic_MZp-200_MChi-50',
                                  'monoH_ZpBaryonic_MZp-20_MChi-1',
                                  'monoH_ZpBaryonic_MZp-295_MChi-150',
                                  'monoH_ZpBaryonic_MZp-300_MChi-1',
                                  'monoH_ZpBaryonic_MZp-300_MChi-50',
                                  'monoH_ZpBaryonic_MZp-500_MChi-150',
                                  'monoH_ZpBaryonic_MZp-500_MChi-1',
                                  'monoH_ZpBaryonic_MZp-500_MChi-500',
                                  'monoH_ZpBaryonic_MZp-50_MChi-10',
                                  'monoH_ZpBaryonic_MZp-50_MChi-1',
                                  'monoH_ZpBaryonic_MZp-95_MChi-50',
                                  'monoH_ZpBaryonic_MZp-995_MChi-500',

                                 ] ,
              },
                  
                  

   'METup':  {
                  'isChain'    : True ,
                  'do4MC'      : True ,
                  'do4Data'    : False,
                  'subTargets' : ['do_METup','l2kin','l3kin','l4kin'],
                  'onlySample' : [
                                  # DY 
                                  'DYJetsToLL_M-10to50','DYJetsToLL_M-50','DYJetsToLL_M-10to50ext3','DYJetsToLL_M-50-LO','DYJetsToLL_M-50-LO-ext1',
                                  'DY2JetsToLL', 'DY3JetsToLL', 'DY4JetsToLL',
                                  'DYJetsToLL_M-50_HT-100to200_MLM' , 'DYJetsToLL_M-50_HT-100to200_MLM_ext' , 
                                  'DYJetsToLL_M-50_HT-200to400_MLM' , 
                                  'DYJetsToLL_M-50_HT-400to600_MLM' ,
                                  'DYJetsToLL_M-50_HT-600toInf_MLM' ,
                                  'DYJetsToLL_M-50_HT-100to200_ext1','DYJetsToLL_M-50_HT-100to200_ext1','DYJetsToLL_M-50_HT-200to400_ext1',
                                  'DYJetsToLL_M-50_HT-200to400_ext1','DYJetsToLL_M-50_HT-200to400_ext1','DYJetsToLL_M-50_HT-600toInf_ext1',
                                  'DYJetsToLL_M-50_HT-600toInf', 'DYJetsToTT_MuEle_M-50',
                                  # Top
                                  'TTTo2L2Nu','TT','TTTo2L2Nu_ext1',
                                  'ST_t-channel_antitop','ST_t-channel_top',
                                  'ST_tW_antitop','ST_tW_top',
                                  # VV (including WW) 
                                  'WWTo2L2Nu','GluGluWWTo2L2Nu_MCFM','GluGluWWTo2L2NuHiggs_MCFM',
                                  'WZTo3LNu','WZ',
                                  'ZZ','Zg',
                                  'ZZTo2L2Nu','ZZTo2L2Q','ZZTo4L',
                                  'Wg_AMCNLOFXFX','WgStarLNuEE','WgStarLNuMuMu',
                                  'Wg_MADGRAPHMLM','WZTo3LNu_mllmin01',
                                  'ggZZ2e2t', 'ggZZ4t', 'ggZZ4m', 'ggZZ2m2t', 'ggZZ2e2m', 'ggZZ4e',
                                  # VVV
                                  'WZZ','ZZZ','WWZ','WWW',
                                  # Higgs 
                                  'GluGluHToTauTau_M125', 'GluGluHToWWTo2L2Nu_M125','GluGluHToWWTo2L2NuPowheg_M125',
                                  'GluGluHToWWTo2L2Nu_alternative_M125','VBFHToWWTo2L2Nu_alternative_M125',
                                  'HWminusJ_HToTauTau_M125', 'HWminusJ_HToWW_M125',
                                  'HWplusJ_HToTauTau_M125', 'HWplusJ_HToWW_M125',
                                  'HZJ_HToTauTau_M125', 'HZJ_HToWW_M125',
                                  'VBFHToTauTau_M125', 'VBFHToWWTo2L2Nu_M125',
                                  'ggZH_HToWW_M125', # missing ggZHToTauTau
                                  'ttHJetToNonbb_M125','TTWJetsToLNu',
                                  'GluGluHToZZTo4L_M125',
                                  'GluGluZH_HToWWTo2L2Nu_M120_noHLT','GluGluZH_HToWWTo2L2Nu_M125_noHLT','GluGluZH_HToWWTo2L2Nu_M130_noHLT',
                                  'HZJ_HToWWTo2L2Nu_M120_noHLT','HZJ_HToWWTo2L2Nu_M125_noHLT','HZJ_HToWWTo2L2Nu_M130_noHLT', 
                                  'HWplusJ_WToLNu_HToWWTo2L2Nu_M125','HWminusJ_WToLNu_HToWWTo2L2Nu_M125',
                                  # ... ggH High Mass
                                  'GluGluHToWWTo2L2Nu_M130',
                                  'GluGluHToWWTo2L2Nu_M135',
                                  'GluGluHToWWTo2L2Nu_M140',
                                  'GluGluHToWWTo2L2Nu_M145',
                                  'GluGluHToWWTo2L2Nu_M150',
                                  'GluGluHToWWTo2L2Nu_M155',
                                  'GluGluHToWWTo2L2Nu_M160',
                                  'GluGluHToWWTo2L2Nu_M165',
                                  'GluGluHToWWTo2L2Nu_M170',
                                  'GluGluHToWWTo2L2Nu_M175',
                                  'GluGluHToWWTo2L2Nu_M180',
                                  'GluGluHToWWTo2L2Nu_M190',
                                  'GluGluHToWWTo2L2Nu_M200',
                                  'GluGluHToWWTo2L2Nu_M210',
                                  'GluGluHToWWTo2L2Nu_M230',
                                  'GluGluHToWWTo2L2Nu_M250',
                                  'GluGluHToWWTo2L2Nu_M270',
                                  'GluGluHToWWTo2L2Nu_M300',
                                  'GluGluHToWWTo2L2Nu_M350',
                                  'GluGluHToWWTo2L2Nu_M400',
                                  'GluGluHToWWTo2L2Nu_M450',
                                  'GluGluHToWWTo2L2Nu_M500',
                                  'GluGluHToWWTo2L2Nu_M550',
                                  'GluGluHToWWTo2L2Nu_M600',
                                  'GluGluHToWWTo2L2Nu_M650',
                                  'GluGluHToWWTo2L2Nu_M700',
                                  'GluGluHToWWTo2L2Nu_M750',
                                  'GluGluHToWWTo2L2Nu_M750_NWA',
                                  'GluGluHToWWTo2L2Nu_M800',
                                  'GluGluHToWWTo2L2Nu_M900',
                                  'GluGluHToWWTo2L2Nu_M1000',
                                  # ... VBF High Mass
                                  'VBFHToWWTo2L2Nu_M130',
                                  'VBFHToWWTo2L2Nu_M135',
                                  'VBFHToWWTo2L2Nu_M140',
                                  'VBFHToWWTo2L2Nu_M145',
                                  'VBFHToWWTo2L2Nu_M150',
                                  'VBFHToWWTo2L2Nu_M155',
                                  'VBFHToWWTo2L2Nu_M160',
                                  'VBFHToWWTo2L2Nu_M165',
                                  'VBFHToWWTo2L2Nu_M170',
                                  'VBFHToWWTo2L2Nu_M175',
                                  'VBFHToWWTo2L2Nu_M180',
                                  'VBFHToWWTo2L2Nu_M190',
                                  'VBFHToWWTo2L2Nu_M200',
                                  'VBFHToWWTo2L2Nu_M210',
                                  'VBFHToWWTo2L2Nu_M230',
                                  'VBFHToWWTo2L2Nu_M250',
                                  'VBFHToWWTo2L2Nu_M270',
                                  'VBFHToWWTo2L2Nu_M300',
                                  'VBFHToWWTo2L2Nu_M350',
                                  'VBFHToWWTo2L2Nu_M400',
                                  'VBFHToWWTo2L2Nu_M450',
                                  'VBFHToWWTo2L2Nu_M500',
                                  'VBFHToWWTo2L2Nu_M550',
                                  'VBFHToWWTo2L2Nu_M600',
                                  'VBFHToWWTo2L2Nu_M650',
                                  'VBFHToWWTo2L2Nu_M700',
                                  'VBFHToWWTo2L2Nu_M750',
                                  'VBFHToWWTo2L2Nu_M750_NWA',
                                  'VBFHToWWTo2L2Nu_M800',
                                  'VBFHToWWTo2L2Nu_M900',
                                  'VBFHToWWTo2L2Nu_M1000',
                                  # VBF 
                                  'WpWmJJ_EWK_QCD_noTop','WpWmJJ_QCD_noTop', 'WpWmJJ_EWK_noTop',
                                  'WGJJ', 'WLLJJToLNu_M-4to60_EWK_QCD', 'WLLJJToLNu_M-60_EWK_QCD', 'WpWpJJ_EWK_QCD',
                                  'WW_DoubleScattering','WWTo2L2Nu_DoubleScattering', 'WmWmJJ_EWK_powheg', 'WpWpJJ_EWK', 'WpWpJJ_QCD', 

                                  # stop

                                  'latino_T2bW',
                                  'latino_T2tt_dM10to80',
                                  'latino_T2tt_dM10to80_genHT160_genMET80',
                                  'latino_T2tt_mStop-150to250',
                                  'latino_T2tt_mStop-250to350',
                                  'latino_T2tt_mStop-350to400',
                                  'latino_T2tt_mStop-400to1200',
                                  'latino_T2tt_mStop425_mLSP325',
                                  'latino_T2tt_mStop500_mLSP325',
                                  'latino_T2tt_mStop850_mLSP100',


                                  # monoH
                                  'monoH_2HDM_MZp-1000_MA0-400',
                                  'monoH_2HDM_MZp-1000_MA0-500',
                                  'monoH_2HDM_MZp-1000_MA0-600',
                                  'monoH_2HDM_MZp-1000_MA0-700',
                                  'monoH_2HDM_MZp-1000_MA0-800',
                                  'monoH_2HDM_MZp-1200_MA0-400',
                                  'monoH_2HDM_MZp-1200_MA0-500',
                                  'monoH_2HDM_MZp-1200_MA0-600',
                                  'monoH_2HDM_MZp-1200_MA0-700',
                                  'monoH_2HDM_MZp-1200_MA0-800',
                                  'monoH_2HDM_MZp-1400_MA0-400',
                                  'monoH_2HDM_MZp-1400_MA0-500',
                                  'monoH_2HDM_MZp-1400_MA0-600',
                                  'monoH_2HDM_MZp-1400_MA0-700',
                                  'monoH_2HDM_MZp-1400_MA0-800',
                                  'monoH_2HDM_MZp-1700_MA0-400',
                                  'monoH_2HDM_MZp-1700_MA0-500',
                                  'monoH_2HDM_MZp-1700_MA0-600',
                                  'monoH_2HDM_MZp-1700_MA0-700',
                                  'monoH_2HDM_MZp-1700_MA0-800',
                                  'monoH_2HDM_MZp-2000_MA0-400',
                                  'monoH_2HDM_MZp-2000_MA0-500',
                                  'monoH_2HDM_MZp-2000_MA0-600',
                                  'monoH_2HDM_MZp-2000_MA0-700',
                                  'monoH_2HDM_MZp-2000_MA0-800',
                                  'monoH_2HDM_MZp-2500_MA0-400',
                                  'monoH_2HDM_MZp-2500_MA0-500',
                                  'monoH_2HDM_MZp-2500_MA0-600',
                                  'monoH_2HDM_MZp-2500_MA0-700',
                                  'monoH_2HDM_MZp-2500_MA0-800',
                                  'monoH_2HDM_MZp-600_MA0-400',
                                  'monoH_2HDM_MZp-800_MA0-400',
                                  'monoH_2HDM_MZp-800_MA0-500',
                                  'monoH_2HDM_MZp-800_MA0-600',
                                  'monoH_2HDM_MZp-1000_MA0-300',
                                  'monoH_2HDM_MZp-1200_MA0-300',
                                  'monoH_2HDM_MZp-1400_MA0-300',
                                  'monoH_2HDM_MZp-1700_MA0-300',
                                  'monoH_2HDM_MZp-2000_MA0-300',
                                  'monoH_2HDM_MZp-2500_MA0-300',
                                  'monoH_2HDM_MZp-600_MA0-300',
                                  'monoH_2HDM_MZp-800_MA0-300',

                               
                                  'monoH_ZpBaryonic_MZp-10000_MChi-1000',
                                  'monoH_ZpBaryonic_MZp-10000_MChi-10',
                                  'monoH_ZpBaryonic_MZp-10000_MChi-150',
                                  'monoH_ZpBaryonic_MZp-10000_MChi-1',
                                  'monoH_ZpBaryonic_MZp-10000_MChi-50',
                                  'monoH_ZpBaryonic_MZp-1000_MChi-1000',
                                  'monoH_ZpBaryonic_MZp-1000_MChi-1',
                                  'monoH_ZpBaryonic_MZp-100_MChi-10',
                                  'monoH_ZpBaryonic_MZp-100_MChi-1',
                                  'monoH_ZpBaryonic_MZp-10_MChi-1000',
                                  'monoH_ZpBaryonic_MZp-10_MChi-10',
                                  'monoH_ZpBaryonic_MZp-10_MChi-150',
                                  'monoH_ZpBaryonic_MZp-10_MChi-1',
                                  'monoH_ZpBaryonic_MZp-10_MChi-500',
                                  'monoH_ZpBaryonic_MZp-10_MChi-50',
                                  'monoH_ZpBaryonic_MZp-15_MChi-10',
                                  'monoH_ZpBaryonic_MZp-2000_MChi-1',
                                  'monoH_ZpBaryonic_MZp-2000_MChi-500',
                                  'monoH_ZpBaryonic_MZp-200_MChi-150',
                                  'monoH_ZpBaryonic_MZp-200_MChi-1',
                                  'monoH_ZpBaryonic_MZp-200_MChi-50',
                                  'monoH_ZpBaryonic_MZp-20_MChi-1',
                                  'monoH_ZpBaryonic_MZp-295_MChi-150',
                                  'monoH_ZpBaryonic_MZp-300_MChi-1',
                                  'monoH_ZpBaryonic_MZp-300_MChi-50',
                                  'monoH_ZpBaryonic_MZp-500_MChi-150',
                                  'monoH_ZpBaryonic_MZp-500_MChi-1',
                                  'monoH_ZpBaryonic_MZp-500_MChi-500',
                                  'monoH_ZpBaryonic_MZp-50_MChi-10',
                                  'monoH_ZpBaryonic_MZp-50_MChi-1',
                                  'monoH_ZpBaryonic_MZp-95_MChi-50',
                                  'monoH_ZpBaryonic_MZp-995_MChi-500',


                                 ] ,
              },

   'METdo':  {
                  'isChain'    : True ,
                  'do4MC'      : True ,
                  'do4Data'    : False,
                  'subTargets' : ['do_METdo','l2kin','l3kin','l4kin'],
                  'onlySample' : [
                                  # DY 
                                  'DYJetsToLL_M-10to50','DYJetsToLL_M-50','DYJetsToLL_M-10to50ext3','DYJetsToLL_M-50-LO','DYJetsToLL_M-50-LO-ext1',
                                  'DY2JetsToLL', 'DY3JetsToLL', 'DY4JetsToLL',
                                  'DYJetsToLL_M-50_HT-100to200_MLM' , 'DYJetsToLL_M-50_HT-100to200_MLM_ext' , 
                                  'DYJetsToLL_M-50_HT-200to400_MLM' , 
                                  'DYJetsToLL_M-50_HT-400to600_MLM' ,
                                  'DYJetsToLL_M-50_HT-600toInf_MLM' ,
                                  'DYJetsToLL_M-50_HT-100to200_ext1','DYJetsToLL_M-50_HT-100to200_ext1','DYJetsToLL_M-50_HT-200to400_ext1',
                                  'DYJetsToLL_M-50_HT-200to400_ext1','DYJetsToLL_M-50_HT-200to400_ext1','DYJetsToLL_M-50_HT-600toInf_ext1',
                                  'DYJetsToLL_M-50_HT-600toInf', 'DYJetsToTT_MuEle_M-50',
                                  # Top
                                  'TTTo2L2Nu','TT','TTTo2L2Nu_ext1',
=======
>>>>>>> d6dbd9f034fad9f8b783aff7ebd17a6014db64b1
                                  'ST_t-channel_antitop','ST_t-channel_top',
                                  'ST_tW_antitop','ST_tW_top',
                                  'TTToSemiLeptonic','TT',
                                  # VV (including WW) 
                                  'WWTo2L2Nu','GluGluWWTo2L2Nu_MCFM','GluGluWWTo2L2NuHiggs_MCFM',
                                  'WZTo3LNu','WZ',
                                  'ZZ','Zg',
                                  'ZZTo2L2Nu','ZZTo2L2Q','ZZTo4L',
                                  'Wg_AMCNLOFXFX', 
                                  'Wg_MADGRAPHMLM',
                                  'WZTo2L2Q',
                                  # VVV
                                  'WZZ','ZZZ','WWZ','WWW',
                                  # Higgs 
                                  'GluGluHToTauTau_M125', 'GluGluHToWWTo2L2Nu_M125','GluGluHToWWTo2L2NuPowheg_M125',
                                  'GluGluHToWWTo2L2Nu_alternative_M125','VBFHToWWTo2L2Nu_alternative_M125',
                                  'HWminusJ_HToTauTau_M125', 'HWminusJ_HToWW_M125',
                                  'HWplusJ_HToTauTau_M125', 'HWplusJ_HToWW_M125',
                                  'HZJ_HToTauTau_M125', 'HZJ_HToWW_M125',
                                  'VBFHToTauTau_M125', 'VBFHToWWTo2L2Nu_M125',
                                  'ggZH_HToWW_M125', # missing ggZHToTauTau
                                  'ttHJetToNonbb_M125',
                                  'GluGluHToZZTo4L_M125',
                                  'GluGluZH_HToWWTo2L2Nu_M120_noHLT','GluGluZH_HToWWTo2L2Nu_M125_noHLT','GluGluZH_HToWWTo2L2Nu_M130_noHLT',
                                  'HZJ_HToWWTo2L2Nu_M120_noHLT','HZJ_HToWWTo2L2Nu_M125_noHLT','HZJ_HToWWTo2L2Nu_M130_noHLT', 
                                  'HWplusJ_WToLNu_HToWWTo2L2Nu_M125','HWminusJ_WToLNu_HToWWTo2L2Nu_M125',
                                  # ... ggH High Mass
                                  'GluGluHToWWTo2L2Nu_M130',
                                  'GluGluHToWWTo2L2Nu_M135',
                                  'GluGluHToWWTo2L2Nu_M140',
                                  'GluGluHToWWTo2L2Nu_M145',
                                  'GluGluHToWWTo2L2Nu_M150',
                                  'GluGluHToWWTo2L2Nu_M155',
                                  'GluGluHToWWTo2L2Nu_M160',
                                  'GluGluHToWWTo2L2Nu_M165',
                                  'GluGluHToWWTo2L2Nu_M170',
                                  'GluGluHToWWTo2L2Nu_M175',
                                  'GluGluHToWWTo2L2Nu_M180',
                                  'GluGluHToWWTo2L2Nu_M190',
                                  'GluGluHToWWTo2L2Nu_M200',
                                  'GluGluHToWWTo2L2Nu_M210',
                                  'GluGluHToWWTo2L2Nu_M230',
                                  'GluGluHToWWTo2L2Nu_M250',
                                  'GluGluHToWWTo2L2Nu_M270',
                                  'GluGluHToWWTo2L2Nu_M300',
                                  'GluGluHToWWTo2L2Nu_M350',
                                  'GluGluHToWWTo2L2Nu_M400',
                                  'GluGluHToWWTo2L2Nu_M450',
                                  'GluGluHToWWTo2L2Nu_M500',
                                  'GluGluHToWWTo2L2Nu_M550',
                                  'GluGluHToWWTo2L2Nu_M600',
                                  'GluGluHToWWTo2L2Nu_M650',
                                  'GluGluHToWWTo2L2Nu_M700',
                                  'GluGluHToWWTo2L2Nu_M750',
                                  'GluGluHToWWTo2L2Nu_M750_NWA',
                                  'GluGluHToWWTo2L2Nu_M800',
                                  'GluGluHToWWTo2L2Nu_M900',
                                  'GluGluHToWWTo2L2Nu_M1000',
                                  # ... VBF High Mass
                                  'VBFHToWWTo2L2Nu_M130',
                                  'VBFHToWWTo2L2Nu_M135',
                                  'VBFHToWWTo2L2Nu_M140',
                                  'VBFHToWWTo2L2Nu_M145',
                                  'VBFHToWWTo2L2Nu_M150',
                                  'VBFHToWWTo2L2Nu_M155',
                                  'VBFHToWWTo2L2Nu_M160',
                                  'VBFHToWWTo2L2Nu_M165',
                                  'VBFHToWWTo2L2Nu_M170',
                                  'VBFHToWWTo2L2Nu_M175',
                                  'VBFHToWWTo2L2Nu_M180',
                                  'VBFHToWWTo2L2Nu_M190',
                                  'VBFHToWWTo2L2Nu_M200',
                                  'VBFHToWWTo2L2Nu_M210',
                                  'VBFHToWWTo2L2Nu_M230',
                                  'VBFHToWWTo2L2Nu_M250',
                                  'VBFHToWWTo2L2Nu_M270',
                                  'VBFHToWWTo2L2Nu_M300',
                                  'VBFHToWWTo2L2Nu_M350',
                                  'VBFHToWWTo2L2Nu_M400',
                                  'VBFHToWWTo2L2Nu_M450',
                                  'VBFHToWWTo2L2Nu_M500',
                                  'VBFHToWWTo2L2Nu_M550',
                                  'VBFHToWWTo2L2Nu_M600',
                                  'VBFHToWWTo2L2Nu_M650',
                                  'VBFHToWWTo2L2Nu_M700',
                                  'VBFHToWWTo2L2Nu_M750',
                                  'VBFHToWWTo2L2Nu_M750_NWA',
                                  'VBFHToWWTo2L2Nu_M800',
                                  'VBFHToWWTo2L2Nu_M900',
                                  'VBFHToWWTo2L2Nu_M1000',
                                  # PS
                                  'GluGluHToWWTo2L2NuHerwigPS_M125','VBFHToWWTo2L2NuHerwigPS_M125','WWTo2L2NuHerwigPS',
                                  # UE
                                  'GluGluHToWWTo2L2Nu_M125_CUETDown',
                                  'GluGluHToWWTo2L2Nu_M125_CUETUp',
                                  'VBFHToWWTo2L2Nu_M125_CUETDown',
                                  'VBFHToWWTo2L2Nu_M125_CUETUp',
                                  'WWTo2L2Nu_CUETDown',
                                  'WWTo2L2Nu_CUETUp',
 			  # VBS
 			  'WpWpJJ_EWK','WpWpJJ_EWK_QCD','WpWpJJ_QCD','WW_DoubleScattering','WWTo2L2Nu_DoubleScattering','WLLJJToLNu_M-4to60_EWK_QCD','WLLJJToLNu_M-60_EWK_QCD',
                          'WGJJ','EWKZ2Jets','TTToSemiLepton',
                                  # VBF 
                                  'WpWmJJ_EWK_QCD_noTop','WpWmJJ_QCD_noTop', 'WpWmJJ_EWK_noTop',
                                 ] ,
                },


  'JESup'     :  {
                  'isChain'    : True ,
                  'do4MC'      : True ,
                  'do4Data'    : False,
                  'subTargets' : ['do_JESup','bPogSF','l2kin','l3kin','l4kin','do_dymvaHiggs'],
                  'onlySample' : samples4Syst
                },

  'JESdo'     :  {
                  'isChain'    : True ,
                  'do4MC'      : True ,
                  'do4Data'    : False,
                  'subTargets' : ['do_JESdo','bPogSF','l2kin','l3kin','l4kin','do_dymvaHiggs'],
                  'onlySample' : samples4Syst
                },


  'JESMaxup'     :  {
                  'isChain'    : True ,
                  'do4MC'      : True ,
                  'do4Data'    : False,
                  'subTargets' : ['do_JESMaxup','bPogSF','l2kin','l3kin','l4kin'],
                  'onlySample' : samples4Syst
                },

  'JESMaxdo'     :  {
                  'isChain'    : True ,
                  'do4MC'      : True ,
                  'do4Data'    : False,
                  'subTargets' : ['do_JESMaxdo','bPogSF','l2kin','l3kin','l4kin'],
                  'onlySample' : samples4Syst
                },

  'LepElepTCutup':  {
                  'isChain'    : True ,
                  'do4MC'      : True ,
                  'do4Data'    : False,
#                  'subTargets' : ['do_LepElepTup','TrigEff','IdIsoSC','l2kin','l3kin','l4kin'],
                  'subTargets' : ['do_LepElepTup','TrigMakerMCkeepRun','IdIsoSC_Cut','l2kin','l3kin','l4kin','do_dymvaHiggs'],
                  'onlySample' : samples4Syst
              },

  'LepElepTCutdo':  {
                  'isChain'    : True ,
                  'do4MC'      : True ,
                  'do4Data'    : False,
                  #'subTargets' : ['do_LepElepTdo','TrigEff','IdIsoSC','l2kin','l3kin','l4kin'],
                  'subTargets' : ['do_LepElepTdo','TrigMakerMCkeepRun','IdIsoSC_Cut','l2kin','l3kin','l4kin','do_dymvaHiggs'],
                  'onlySample' : samples4Syst
              },


  'LepMupTCutup':  {
                  'isChain'    : True ,
                  'do4MC'      : True ,
                  'do4Data'    : False,
                  #'subTargets' : ['do_LepMupTup','TrigEff','IdIsoSC','l2kin','l3kin','l4kin'],
                  'subTargets' : ['do_LepMupTup','TrigMakerMCkeepRun','IdIsoSC_Cut','l2kin','l3kin','l4kin','do_dymvaHiggs'],
                  'onlySample' : samples4Syst
              },

  'LepMupTCutdo':  {
                  'isChain'    : True ,
                  'do4MC'      : True ,
                  'do4Data'    : False,
                  #'subTargets' : ['do_LepMupTdo','TrigEff','IdIsoSC','l2kin','l3kin','l4kin'],
                  'subTargets' : ['do_LepMupTdo','TrigMakerMCkeepRun','IdIsoSC_Cut','l2kin','l3kin','l4kin','do_dymvaHiggs'],
                  'onlySample' : samples4Syst
              },
                  
                  

   'METup':  {
                  'isChain'    : True ,
                  'do4MC'      : True ,
                  'do4Data'    : False,
                  'subTargets' : ['do_METup','l2kin','l3kin','l4kin','do_dymvaHiggs'],
                  'onlySample' : samples4Syst
              },

   'METdo':  {
                  'isChain'    : True ,
                  'do4MC'      : True ,
                  'do4Data'    : False,
                  'subTargets' : ['do_METdo','l2kin','l3kin','l4kin','do_dymvaHiggs'],
                  'onlySample' : samples4Syst
              },

# ... Individual Steps

  'mcwghtcount':{ 
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False ,
                  'command'    : 'gardener.py mcweightscounter '
                },

# 'mcweights' : {
#                 'isChain'    : False ,
#                 'do4MC'      : True  ,
#                 'do4Data'    : False ,
#                 'command'    : 'gardener.py mcweightsfiller '
#               } ,


  'tagjson'   : {
                  'isChain'    : False ,
                  'do4MC'      : False ,
                  'do4Data'    : True ,
                  #'command'    : 'gardener.py  filterjson --json=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-273450_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt'
                  #'command'    : 'gardener.py  filterjson --json=/user/xjanssen/HWW2015/pudata/Cert_271036-273450_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt'
                  'command'    : 'gardener.py  filterjson --json=/user/xjanssen/HWW2015/pudata/Cert_271036-276811_13TeV_PromptReco_Collisions16_JSON.txt'
                } ,

  'tagjsonICHEP'   : {
                  'isChain'    : False ,
                  'do4MC'      : False ,
                  'do4Data'    : True ,
                  #'command'    : 'gardener.py  filterjson --json=/user/xjanssen/HWW2015/pudata/XXXXXXX.txt'
                  'command'    : 'gardener.py  filterjson --json=/user/xjanssen/HWW2015/pudata/Cert_271036-276811_13TeV_PromptReco_Collisions16_JSON.txt'
                } ,



  'selectjson'   : {
                  'isChain'    : False ,
                  'do4MC'      : False ,
                  'do4Data'    : True ,
                  'command'    : 'gardener.py filter -f \'isJsonOk>0.5\' '
                } ,

  'fakeW12fb'     : {
                  'isChain'    : False ,
                  'do4MC'      : False ,
                  'do4Data'    : True ,
                  'command'    : 'gardener.py fakeWeights ',
                } ,

  'fakeW12fb_v2'     : {
                  'isChain'    : False ,
                  'do4MC'      : False ,
                  'do4Data'    : True ,
                  'command'    : 'gardener.py fakeWeights ',
                } ,

  'fakeWCut' : {
                  'isChain'    : False ,
                  'do4MC'      : False ,
                  'do4Data'    : True ,
                  'command'    : 'gardener.py fakeWeights --cmssw RPLME_CMSSW --idEleKind cut_WP_Tight80X',
                } ,


  'puadder'   : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False ,
                  'command'    : 'gardener.py puadder --data=RPLME_puData --HistName=pileup --branch=puW --kind=trpu '
                } ,

  'puRunPer'   : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False ,
                  'command'    : 'gardener.py puadder --data=RPLME_puData --HistName=pileup --branch=puW --kind=trpu --run --cmssw RPLME_CMSSW'
                } ,

  'pu2p6'   : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False ,
                  'command'    : 'gardener.py puadder --data=/user/xjanssen/HWW2015/pudata/PileupHistogram_805_28Jun.root --HistName=pileup --branch=puW2p6 --kind=trpu '
                }, 

  'pu4p3'   : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False ,
                  'command'    : 'gardener.py puadder --data=/user/xjanssen/HWW2015/pudata/PileupHistogram_805_5Jul_4fb.root --HistName=pileup --branch=puW4p3 --kind=trpu '
                },


  'pu6p3'   : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False ,
                  'command'    : 'gardener.py puadder --data=/user/xjanssen/HWW2015/pudata/PileupHistogram_805_8Jul_6p3fb.root --HistName=pileup --branch=puW6p3 --kind=trpu '
                },

  'puW63mb'  : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False ,
                  'command'    : 'gardener.py puadder --data=/user/xjanssen/HWW2015/pudata/PileupHistogram_805_6.3fb_63mb_26Jul.root  --HistName=pileup --branch=puW63mb --kind=trpu '
                }, 



  'puW69mb'  : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False ,
                  'command'    : 'gardener.py puadder --data=/user/xjanssen/HWW2015/pudata/PileupHistogram_805_6.3fb_69.2mb_26Jul.root  --HistName=pileup --branch=puW69mb --kind=trpu '
                },


  'baseW'     : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False ,
                  'command'    : 'gardener.py adder -v \'baseW/F=RPLME_baseW\' -v \'Xsec/F=RPLME_XSection\' '

                } ,

  'wwNLL'     : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False ,
                  'onlySample' : ['WWTo2L2Nu','WWTo2L2NuHerwigPS','WWTo2L2Nu_CUETUp','WWTo2L2Nu_CUETDown'] ,
                  'command'    : 'gardener.py wwNLLcorrections -m \'powheg\' --cmssw RPLME_CMSSW'
                },

  'wwEWK'     : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False ,
                  'onlySample' : ['WWTo2L2Nu','WWTo2L2NuHerwigPS','WWTo2L2Nu_CUETUp','WWTo2L2Nu_CUETDown'] ,
                  'command'    : 'gardener.py wwEWKcorrections' ,
                } ,

  'wzEWK'     : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False ,
                  'onlySample' : ['WZTo1L1Nu2Q', 'WZTo1L3Nu', 'WZTo2L2Q', 'WZTo3LNu', 'WZTo3LNu_mllmin01', 'WZTo3LNu_mllmin01_ext1', ],
                  'command'    : 'gardener.py wzEWKcorrections' ,
                } ,

  'zzEWK'     : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False ,
                  'onlySample' : ['ZZTo2L2Nu', 'ZZTo2L2Q', ],
                  'command'    : 'gardener.py zzEWKcorrections' ,
                } ,


  'genVariables'  : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False ,
                  'command'    : 'gardener.py genvariablesfiller ',
                  'onlySample' : [
                                  # DY 
                                  'DYJetsToLL_M-5to50-LO',
                                  'DYJetsToLL_M-50-LO-ext1',
                                  'DYJetsToLL_M-10to50-LO' ,
                                  'DYJetsToLL_M-10to50','DYJetsToLL_M-50','DYJetsToLL_M-10to50ext3','DYJetsToLL_M-50-LO',
                                  # ....
                                  'DYJetsToLL_M-50-PSdo' , 'DYJetsToLL_M-50-PSup' , 'DYJetsToLL_M-50-UEdo' , 'DYJetsToLL_M-50-UEup', 
                                  # ... ICHEP16
                                  'DY2JetsToLL', 'DY3JetsToLL', 'DY4JetsToLL','DYJetsToTT_MuEle_M-50',
                                  'DYJetsToLL_M-50_HT-100to200_MLM' , 'DYJetsToLL_M-50_HT-100to200_MLM_ext' , 
                                  'DYJetsToLL_M-50_HT-200to400_MLM' , 
                                  'DYJetsToLL_M-50_HT-400to600_MLM' ,
                                  'DYJetsToLL_M-50_HT-600toInf_MLM' ,
                                  'DYJetsToEE_Pow' ,
                                  'DYJetsToLL_M-50_HT-400to600_ext1' ,
                                  'DYJetsToLL_M-50_HT-600toInf_ext1' ,
                                   # ... Moriond17
                                  'DYJetsToLL_M-5to50_HT-70to100'  , 'DYJetsToLL_M-5to50_HT-100to200' , 'DYJetsToLL_M-5to50_HT-200to400' ,
                                  'DYJetsToLL_M-5to50_HT-400to600' , 'DYJetsToLL_M-5to50_HT-600toInf' ,
                                  'DYJetsToLL_M-50_HT-70to100'  , 'DYJetsToLL_M-50_HT-100to200' , 'DYJetsToLL_M-50_HT-100to200_ext1' ,
                                  'DYJetsToLL_M-50_HT-200to400' , 'DYJetsToLL_M-50_HT-200to400_ext1' , 'DYJetsToLL_M-50_HT-400to600' ,
                                  'DYJetsToLL_M-50_HT-600to800' , 'DYJetsToLL_M-50_HT-800to1200' , 'DYJetsToLL_M-50_HT-1200to2500' ,
                                  'DYJetsToLL_M-50_HT-2500toInf',
                                   # ... DY -> MuEle 
                                   'DYJetsToTT_MuEle_M-50' , 'DYJetsToTT_MuEle_M-50_ext1' ,

                                  # WW ewk
                                  'WpWmJJ_EWK_QCD_noTop','WpWmJJ_QCD_noTop','WpWmJJ_EWK_noTop',
                                  # WW
                                  'WWTo2L2Nu','WWTo2L2NuHerwigPS','WWTo2L2Nu_CUETUp','WWTo2L2Nu_CUETDown' ,
                                  'WWTo2L2Nu_DoubleScattering',

                                  
                 # Higgs 
                 'GluGluHToTauTau_M125' , 'GluGluHToWWTo2L2Nu_M125','GluGluHToWWTo2L2NuPowheg_M125',
                 'GluGluHToWWTo2L2Nu_alternative_M125',
                 'VBFHToTauTau_M125' , 'VBFHToWWTo2L2Nu_M125','VBFHToWWTo2L2Nu_alternative_M125',
                 'HWminusJ_HToTauTau_M125', 'HWminusJ_HToWW_M125', 'HWminusJ_HToWW_LNu_M125' ,
                 'HWplusJ_HToTauTau_M125' , 'HWplusJ_HToWW_M125',  'HWplusJ_HToWW_LNu_M125' ,
                 'HZJ_HToTauTau_M125'     , 'HZJ_HToWW_M125', 'HZJ_HToWWTo2L2Nu_M125',
                 'bbHToWWTo2L2Nu_M125_yb2', 'bbHToWWTo2L2Nu_M125_ybyt',

                 'ggZH_HToWW_M125', # missing ggZHToTauTau

                 # What ????
                 #'ttHJetToNonbb_M125','TTWJetsToLNu',
                 #'GluGluHToZZTo4L_M125',
                 #'GluGluZH_HToWWTo2L2Nu_M120_noHLT','GluGluZH_HToWWTo2L2Nu_M125_noHLT','GluGluZH_HToWWTo2L2Nu_M130_noHLT',
                 #'HZJ_HToWWTo2L2Nu_M120_noHLT','HZJ_HToWWTo2L2Nu_M125_noHLT','HZJ_HToWWTo2L2Nu_M130_noHLT', 
                 #'HWplusJ_WToLNu_HToWWTo2L2Nu_M125','HWminusJ_WToLNu_HToWWTo2L2Nu_M125',

                 # ... ggH High Mass
                 'GluGluHToWWTo2L2Nu_M130',
                 'GluGluHToWWTo2L2Nu_M135',
                 'GluGluHToWWTo2L2Nu_M140',
                 'GluGluHToWWTo2L2Nu_M145',
                 'GluGluHToWWTo2L2Nu_M150',
                 'GluGluHToWWTo2L2Nu_M155',
                 'GluGluHToWWTo2L2Nu_M160',
                 'GluGluHToWWTo2L2Nu_M165',
                 'GluGluHToWWTo2L2Nu_M170',
                 'GluGluHToWWTo2L2Nu_M175',
                 'GluGluHToWWTo2L2Nu_M180',
                 'GluGluHToWWTo2L2Nu_M190',
                 'GluGluHToWWTo2L2Nu_M200',
                 'GluGluHToWWTo2L2Nu_M210',
                 'GluGluHToWWTo2L2Nu_M230',
                 'GluGluHToWWTo2L2Nu_M250',
                 'GluGluHToWWTo2L2Nu_M270',
                 'GluGluHToWWTo2L2Nu_M300', 'GluGluHToWWTo2L2Nu_JHUGen698_M300' ,
                 'GluGluHToWWTo2L2Nu_M350', 'GluGluHToWWTo2L2Nu_JHUGen698_M350' ,
                 'GluGluHToWWTo2L2Nu_M400', 'GluGluHToWWTo2L2Nu_JHUGen698_M400' ,
                 'GluGluHToWWTo2L2Nu_M450', 'GluGluHToWWTo2L2Nu_JHUGen698_M450' ,
                 'GluGluHToWWTo2L2Nu_M500', 'GluGluHToWWTo2L2Nu_JHUGen698_M500' ,
                 'GluGluHToWWTo2L2Nu_M550', 'GluGluHToWWTo2L2Nu_JHUGen698_M550' ,
                 'GluGluHToWWTo2L2Nu_M600', 'GluGluHToWWTo2L2Nu_JHUGen698_M600' ,
                 'GluGluHToWWTo2L2Nu_M650', 'GluGluHToWWTo2L2Nu_JHUGen698_M650' ,
                 'GluGluHToWWTo2L2Nu_M700', 'GluGluHToWWTo2L2Nu_JHUGen698_M700' ,
                 'GluGluHToWWTo2L2Nu_M750', 'GluGluHToWWTo2L2Nu_JHUGen698_M750' , 'GluGluHToWWTo2L2Nu_M750_NWA',
                 'GluGluHToWWTo2L2Nu_M800', 'GluGluHToWWTo2L2Nu_JHUGen698_M800' ,
                 'GluGluHToWWTo2L2Nu_M900', 'GluGluHToWWTo2L2Nu_JHUGen698_M900' ,
                 'GluGluHToWWTo2L2Nu_M1000','GluGluHToWWTo2L2Nu_JHUGen698_M1000',
                 'GluGluHToWWTo2L2Nu_JHUGen698_M1500' ,
                 'GluGluHToWWTo2L2Nu_JHUGen698_M2000' ,
                 'GluGluHToWWTo2L2Nu_JHUGen698_M2500' ,
                 'GluGluHToWWTo2L2Nu_JHUGen698_M3000' ,
                 # ... VBF High Mass
                 'VBFHToWWTo2L2Nu_M130',
                 'VBFHToWWTo2L2Nu_M135',
                 'VBFHToWWTo2L2Nu_M140',
                 'VBFHToWWTo2L2Nu_M145',
                 'VBFHToWWTo2L2Nu_M150',
                 'VBFHToWWTo2L2Nu_M155',
                 'VBFHToWWTo2L2Nu_M160',
                 'VBFHToWWTo2L2Nu_M165',
                 'VBFHToWWTo2L2Nu_M170',
                 'VBFHToWWTo2L2Nu_M175',
                 'VBFHToWWTo2L2Nu_M180',
                 'VBFHToWWTo2L2Nu_M190',
                 'VBFHToWWTo2L2Nu_M200',
                 'VBFHToWWTo2L2Nu_M210',
                 'VBFHToWWTo2L2Nu_M230',
                 'VBFHToWWTo2L2Nu_M250',
                 'VBFHToWWTo2L2Nu_M270',
                 'VBFHToWWTo2L2Nu_M300', 'VBFHToWWTo2L2Nu_JHUGen698_M300' ,
                 'VBFHToWWTo2L2Nu_M350', 'VBFHToWWTo2L2Nu_JHUGen698_M350' ,
                 'VBFHToWWTo2L2Nu_M400', 'VBFHToWWTo2L2Nu_JHUGen698_M400' ,
                 'VBFHToWWTo2L2Nu_M450', 'VBFHToWWTo2L2Nu_JHUGen698_M450' ,
                 'VBFHToWWTo2L2Nu_M500', 'VBFHToWWTo2L2Nu_JHUGen698_M500' ,
                 'VBFHToWWTo2L2Nu_M550', 'VBFHToWWTo2L2Nu_JHUGen698_M550' ,
                 'VBFHToWWTo2L2Nu_M600', 'VBFHToWWTo2L2Nu_JHUGen698_M600' ,
                 'VBFHToWWTo2L2Nu_M650', 'VBFHToWWTo2L2Nu_JHUGen698_M650' ,
                 'VBFHToWWTo2L2Nu_M700', 'VBFHToWWTo2L2Nu_JHUGen698_M700' ,
                 'VBFHToWWTo2L2Nu_M750', 'VBFHToWWTo2L2Nu_JHUGen698_M750' , 'VBFHToWWTo2L2Nu_M750_NWA',
                 'VBFHToWWTo2L2Nu_M800', 'VBFHToWWTo2L2Nu_JHUGen698_M800' ,
                 'VBFHToWWTo2L2Nu_M900', 'VBFHToWWTo2L2Nu_JHUGen698_M900' ,
                 'VBFHToWWTo2L2Nu_M1000','VBFHToWWTo2L2Nu_JHUGen698_M1000' ,
                 'VBFHToWWTo2L2Nu_JHUGen698_M1500' ,        
                 'VBFHToWWTo2L2Nu_JHUGen698_M2000' ,        
                 'VBFHToWWTo2L2Nu_JHUGen698_M2500' ,        
                 'VBFHToWWTo2L2Nu_JHUGen698_M3000' ,        
                                  ]

                },


  'genMatchVariables'  : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False ,
                  'command'    : 'gardener.py genmatchvarfiller ',
                },

  'TopGenPt' :   {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False ,
                  'command'    : 'gardener.py TopGenPt',
                  'onlySample' : [ 
                                   'TTJetsDiLep-LO-ext1',
                                   'TTJets',
                                   'TTTo2L2Nu_alphaS01108',
                                   'TTTo2L2Nu_ext1',
                                   'TTToSemiLepton','TTToSemiLeptonic',
                                   'TTToSemiLeptonic_alphaS01108',
                                 ],
                  },

  'BWEwkSinglet' : {
                  'isChain'    : False ,
                  'do4MC'      : True ,
                  'do4Data'    : False,
      
                  'onlySample' : [
                                  # ... ggH High Mass
                                  'GluGluHToWWTo2L2Nu_M200',
                                  'GluGluHToWWTo2L2Nu_M210',
                                  'GluGluHToWWTo2L2Nu_M230',
                                  'GluGluHToWWTo2L2Nu_M250',
                                  'GluGluHToWWTo2L2Nu_M270',
                                  'GluGluHToWWTo2L2Nu_M300',
                                  'GluGluHToWWTo2L2Nu_M350',
                                  'GluGluHToWWTo2L2Nu_M400',
                                  'GluGluHToWWTo2L2Nu_M450',
                                  'GluGluHToWWTo2L2Nu_M500',
                                  'GluGluHToWWTo2L2Nu_M550',
                                  'GluGluHToWWTo2L2Nu_M600',
                                  'GluGluHToWWTo2L2Nu_M650',
                                  'GluGluHToWWTo2L2Nu_M700',
                                  'GluGluHToWWTo2L2Nu_M750',
                                  #'GluGluHToWWTo2L2Nu_M750_NWA',
                                  'GluGluHToWWTo2L2Nu_M800',
                                  'GluGluHToWWTo2L2Nu_M900',
                                  'GluGluHToWWTo2L2Nu_M1000',
                                  # .... and new JHU samples
                                  'GluGluHToWWTo2L2Nu_JHUGen698_M200',
                                  'GluGluHToWWTo2L2Nu_JHUGen698_M210',
                                  'GluGluHToWWTo2L2Nu_JHUGen698_M230',
                                  'GluGluHToWWTo2L2Nu_JHUGen698_M250',
                                  'GluGluHToWWTo2L2Nu_JHUGen698_M270',
                                  'GluGluHToWWTo2L2Nu_JHUGen698_M300',
                                  'GluGluHToWWTo2L2Nu_JHUGen698_M350',
                                  'GluGluHToWWTo2L2Nu_JHUGen698_M400',
                                  'GluGluHToWWTo2L2Nu_JHUGen698_M450',
                                  'GluGluHToWWTo2L2Nu_JHUGen698_M500',
                                  'GluGluHToWWTo2L2Nu_JHUGen698_M550',
                                  'GluGluHToWWTo2L2Nu_JHUGen698_M600',
                                  'GluGluHToWWTo2L2Nu_JHUGen698_M650',
                                  'GluGluHToWWTo2L2Nu_JHUGen698_M700',
                                  'GluGluHToWWTo2L2Nu_JHUGen698_M750',
                                  'GluGluHToWWTo2L2Nu_JHUGen698_M800',
                                  'GluGluHToWWTo2L2Nu_JHUGen698_M900',
                                  'GluGluHToWWTo2L2Nu_JHUGen698_M1000',
                                  'GluGluHToWWTo2L2Nu_JHUGen698_M1500',
                                  'GluGluHToWWTo2L2Nu_JHUGen698_M2000',
                                  'GluGluHToWWTo2L2Nu_JHUGen698_M2500',
                                  'GluGluHToWWTo2L2Nu_JHUGen698_M3000',


                                  # ... VBF High Mass
                          #       'VBFHToWWTo2L2Nu_M200',
                          #       'VBFHToWWTo2L2Nu_M210',
                          #       'VBFHToWWTo2L2Nu_M230',
                          #       'VBFHToWWTo2L2Nu_M250',
                          #       'VBFHToWWTo2L2Nu_M270',
                          #       'VBFHToWWTo2L2Nu_M300',
                          #       'VBFHToWWTo2L2Nu_M350',
                          #       'VBFHToWWTo2L2Nu_M400',
                          #       'VBFHToWWTo2L2Nu_M450',
                          #       'VBFHToWWTo2L2Nu_M500',
                          #       'VBFHToWWTo2L2Nu_M550',
                          #       'VBFHToWWTo2L2Nu_M600',
                          #       'VBFHToWWTo2L2Nu_M650',
                          #       'VBFHToWWTo2L2Nu_M700',
                          #       'VBFHToWWTo2L2Nu_M750',
                          #       #'VBFHToWWTo2L2Nu_M750_NWA',
                          #       'VBFHToWWTo2L2Nu_M800',
                          #       'VBFHToWWTo2L2Nu_M900',
                          #       'VBFHToWWTo2L2Nu_M1000', 
                                 ],
                  #'command'    : 'gardener.py BWEwkSingletReweighter -p "latino_(GluGlu|VBF)HToWWTo2L2Nu_M([0-9]+)*"',
                  'command'    : 'gardener.py BWEwkSingletReweighter ',
                 },


  'BWEwkSinglet_JHUGen698' : {
                  'isChain'    : False ,
                  'do4MC'      : True ,
                  'do4Data'    : False,
      
                  'onlySample' : [
                                  # ... ggH High Mass
                                  'GluGluHToWWTo2L2Nu_JHUGen698_M200',
                                  'GluGluHToWWTo2L2Nu_JHUGen698_M210',
                                  'GluGluHToWWTo2L2Nu_JHUGen698_M230',
                                  'GluGluHToWWTo2L2Nu_JHUGen698_M250',
                                  'GluGluHToWWTo2L2Nu_JHUGen698_M270',
                                  'GluGluHToWWTo2L2Nu_JHUGen698_M300',
                                  'GluGluHToWWTo2L2Nu_JHUGen698_M350',
                                  'GluGluHToWWTo2L2Nu_JHUGen698_M400',
                                  'GluGluHToWWTo2L2Nu_JHUGen698_M450',
                                  'GluGluHToWWTo2L2Nu_JHUGen698_M500',
                                  'GluGluHToWWTo2L2Nu_JHUGen698_M550',
                                  'GluGluHToWWTo2L2Nu_JHUGen698_M600',
                                  'GluGluHToWWTo2L2Nu_JHUGen698_M650',
                                  'GluGluHToWWTo2L2Nu_JHUGen698_M700',
                                  'GluGluHToWWTo2L2Nu_JHUGen698_M750',
                                  'GluGluHToWWTo2L2Nu_JHUGen698_M800',
                                  'GluGluHToWWTo2L2Nu_JHUGen698_M900',
                                  'GluGluHToWWTo2L2Nu_JHUGen698_M1000',
                                  #'GluGluHToWWTo2L2Nu_JHUGen698_M1500',
                                  #'GluGluHToWWTo2L2Nu_JHUGen698_M2000',
                                  #'GluGluHToWWTo2L2Nu_JHUGen698_M2500',
                                  #'GluGluHToWWTo2L2Nu_JHUGen698_M3000',

                                  # ... VBF High Mass
                                  'VBFHToWWTo2L2Nu_JHUGen698_M200',
                                  'VBFHToWWTo2L2Nu_JHUGen698_M210',
                                  'VBFHToWWTo2L2Nu_JHUGen698_M230',
                                  'VBFHToWWTo2L2Nu_JHUGen698_M250',
                                  'VBFHToWWTo2L2Nu_JHUGen698_M270',
                                  'VBFHToWWTo2L2Nu_JHUGen698_M300',
                                  'VBFHToWWTo2L2Nu_JHUGen698_M350',
                                  'VBFHToWWTo2L2Nu_JHUGen698_M400',
                                  'VBFHToWWTo2L2Nu_JHUGen698_M450',
                                  'VBFHToWWTo2L2Nu_JHUGen698_M500',
                                  'VBFHToWWTo2L2Nu_JHUGen698_M550',
                                  'VBFHToWWTo2L2Nu_JHUGen698_M600',
                                  'VBFHToWWTo2L2Nu_JHUGen698_M650',
                                  'VBFHToWWTo2L2Nu_JHUGen698_M700',
                                  'VBFHToWWTo2L2Nu_JHUGen698_M750',
                                  'VBFHToWWTo2L2Nu_JHUGen698_M800',
                                  'VBFHToWWTo2L2Nu_JHUGen698_M900',
                                  'VBFHToWWTo2L2Nu_JHUGen698_M1000', 
                                  #'VBFHToWWTo2L2Nu_JHUGen698_M1500', 
                                  #'VBFHToWWTo2L2Nu_JHUGen698_M2000', 
                                  #'VBFHToWWTo2L2Nu_JHUGen698_M2500', 
                                  #'VBFHToWWTo2L2Nu_JHUGen698_M3000', 
                                 ],
                  'command'    : 'gardener.py BWEwkSingletReweighter -p "latino_(GluGlu|VBF)HToWWTo2L2Nu_JHUGen698_M([0-9]+)*" --undoCPS=False',
                 },

  'l2sel'    : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : True  ,
                  'command'    : 'gardener.py l2selfiller --kind 1 --cmssw RPLME_CMSSW'
               },

  'l2kin'    : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : True  ,
                  'command'    : 'gardener.py l2kinfiller --cmssw RPLME_CMSSW'
               },

  'l2kin_metXYshift'    : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : True  ,
                  'command'    : 'gardener.py l2kinfiller --cmssw RPLME_CMSSW --met'
               },

  'l3kin'    : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : True  ,
                  'command'    : 'gardener.py l3kinfiller --cmssw RPLME_CMSSW'
               },


  'l4kin'    : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : True  ,
                  'command'    : 'gardener.py l4kinfiller --cmssw RPLME_CMSSW'
               },


  'do_l2loose'    : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : True  ,
                  'command'    : 'gardener.py l2selfiller --kind 2 --cmssw RPLME_CMSSW --selection 1'
               },


  'do_l2stop'    : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : True  ,
                  'command'    : 'gardener.py l2selfiller --kind 5 --cmssw RPLME_CMSSW --selection 1 --idEleKind RPLME_LepKind'
               },

  'do_l2vloose'    : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : True  ,
                  'command'    : 'gardener.py l2selfiller --kind 4 --cmssw RPLME_CMSSW --selection 1'
               },

  'do_l2loose_Cut'    : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : True  ,
                  'command'    : 'gardener.py l2selfiller --kind 2 --cmssw RPLME_CMSSW --selection 1 --idEleKind cut_WP_Tight80X'
               },

  'do_l2vloose'    : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : True  ,
                  'command'    : 'gardener.py l2selfiller --kind 4 --cmssw RPLME_CMSSW --selection 1'
               },

  'do_l2vloose_Cut'    : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : True  ,
                  'command'    : 'gardener.py l2selfiller --kind 4 --cmssw RPLME_CMSSW --selection 1 --idEleKind cut_WP_Tight80X'
               },


  'l2tight'    : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : True  ,
                  'command'    : 'gardener.py filter -f \'std_vector_lepton_isTightLepton[0]>0.5 && std_vector_lepton_isTightLepton[1]>0.5\' '
               },

  'l1loose'    : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : True  ,
                  'command'    : 'gardener.py l1selfiller --kind 2 --cmssw RPLME_CMSSW'
               },

  'l1looseCut'    : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : True  ,
                  'command'    : 'gardener.py l1selfiller --kind 2 --cmssw RPLME_CMSSW --idEleKind cut_WP_Tight80X'
               },

  'l1vloose'    : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : True  ,
                  'command'    : 'gardener.py l1selfiller --kind 4 --cmssw RPLME_CMSSW'
               },

  'l1vlooseCut'    : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : True  ,
                  'command'    : 'gardener.py l1selfiller --kind 4 --cmssw RPLME_CMSSW --idEleKind cut_WP_Tight80X'
               },


  'do_WgStarsel' : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : True  ,
                  'command'    : 'gardener.py l2selfiller --kind 3 --cmssw RPLME_CMSSW'
               },


  'IdIsoSC'  : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False  ,
                  'command'    : 'gardener.py idisofiller  --isoideleAltLumiRatio=0.0404229 --cmssw=RPLME_CMSSW'
               },
               # the number is 0.497/fb / XXX/fb 
               # then 0.497 / 2.6 = 0.19
               # then 0.497 / 4.0 = 0.12
               # then 0.497 / (2.791 + 1.546 + 1.549 + 0.378) = 0.497 / 6.264 = 0.079
               # then 0.497 / (2.791 + 1.546 + 1.549 + 0.378) = 0.497 / 6.264 = 0.079
               # then 0.497 / 12.2950 =   0.0404229              
               #
               # export PATH=$HOME/.local/bin:/afs/cern.ch/cms/lumi/brilconda-1.0.3/bin:$PATH
               # brilcalc lumi --begin  273158 --end 273726 -u /fb -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-274443_13TeV_PromptReco_Collisions16_JSON.txt
               # 0.497 /fb

  'IdIsoSC_Cut'  : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False  ,
                  'command'    : 'gardener.py idisofiller  --isoideleAltLumiRatio=0.0135 --cmssw=RPLME_CMSSW --idEleKind=cut_WP_Tight80X'
               },

  'IdIsoSCStop' : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False  ,
                  'command'    : 'gardener.py idisostopfiller  --readfastsim=1 --cmssw=RPLME_CMSSW --idLepKind RPLME_LepKind --BCDEFtoGHRatio=1.0'
                  }, # Change BCDEFtoGH according to data periods

# Old trigger module


  'TrigEff'  : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False  ,
                  'command'    : 'gardener.py efftfiller  --fixMuonTriggerLumiRatio=0.0483937   --cmssw=RPLME_CMSSW' 
               }, 
               # the number is 0.595/fb / XXX/fb 
               # then 0.595 / 2.6 = 0.23
               # then 0.595 / 4.0 = 0.15
               # then 0.595 / (2.791 + 1.546 + 1.549 + 0.378) = 0.595 / 6.264 = 0.095
               # then 0.595 / 12.2950 =    0.0483937

  'TrigEffB'  : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False  ,
                  'command'    : 'gardener.py efftfiller  --fixMuonTriggerLumiRatio=0.1009844   --cmssw=RPLME_CMSSW'
               }, 


  'TrigEff_Cut'  : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False  ,
                  'command'    : 'gardener.py efftfiller  --fixMuonTriggerLumiRatio=0.0135   --cmssw=RPLME_CMSSW'
               },

# New trigger module

  'TrigMakerMC'    : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False  ,
                  'command'    : 'gardener.py trigMaker  --cmssw=RPLME_CMSSW'
                 },

  'TrigMakerMCkeepRun'    : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False  ,
                  'command'    : 'gardener.py trigMaker  --cmssw=RPLME_CMSSW --keeprun'
                 },

  'TrigMakerData'    : {
                  'isChain'    : False ,
                  'do4MC'      : False  ,
                  'do4Data'    : True  ,
                  'command'    : 'gardener.py trigMaker  --cmssw=RPLME_CMSSW -d'
                 },
 

  'hadd'     : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : True  ,
                  'SizeMethod' : True , 
                  'SizeMax'    : 5e9 , 
                  'bigSamples' : ['DYJetsToLL_M-50','DY2JetsToLL','ZZTo2L2Q','DYJetsToLL_M-50-LO',
                                  'DYJetsToLL_M-50-LO-ext1',
                                  'WZTo2L2Q','TTToSemiLepton','TTToSemiLeptonic','TTTo2L2Nu_ext1','TTJetsDiLep-LO-ext1','TTTo2L2Nu',
                                  'DYJetsToEE_Pow',
                                  'DY1JetsToLL',
                                  #'TTJets',
                                 ],
               },

  'UEPS'     : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False  ,
                  'onlySample' : [ 
                                    'GluGluHToWWTo2L2Nu_M125_CUETDown' , 'VBFHToWWTo2L2Nu_M125_CUETDown' , 'WWTo2L2Nu_CUETDown' ,
                                    'GluGluHToWWTo2L2Nu_M125_CUETUp'   , 'VBFHToWWTo2L2Nu_M125_CUETUp'   , 'WWTo2L2Nu_CUETUp'   ,
                                    'GluGluHToWWTo2L2NuHerwigPS_M125'  , 'VBFHToWWTo2L2NuHerwigPS_M125'  , 'WWTo2L2NuHerwigPS'  ,
                                 ] ,
                  'cpMap' : {
                              'UEdo' : { 
                                          'GluGluHToWWTo2L2Nu_M125_CUETDown' : ['GluGluHToWWTo2L2Nu_M125' ,'GluGluHToWWTo2L2NuPowheg_M125'],
                                          'VBFHToWWTo2L2Nu_M125_CUETDown'    : ['VBFHToWWTo2L2Nu_M125']    ,
                                          'WWTo2L2Nu_CUETDown'               : ['WWTo2L2Nu'] , 
                                       },
                              'UEup' : {
                                          'GluGluHToWWTo2L2Nu_M125_CUETUp'   : ['GluGluHToWWTo2L2Nu_M125' ,'GluGluHToWWTo2L2NuPowheg_M125'],
                                          'VBFHToWWTo2L2Nu_M125_CUETUp'      : ['VBFHToWWTo2L2Nu_M125']    ,
                                          'WWTo2L2Nu_CUETUp'                 : ['WWTo2L2Nu'] ,
                                       },
                              'PS'   : {    
                                          'GluGluHToWWTo2L2NuHerwigPS_M125'  : ['GluGluHToWWTo2L2Nu_M125' ,'GluGluHToWWTo2L2NuPowheg_M125'],
                                          'VBFHToWWTo2L2NuHerwigPS_M125'     : ['VBFHToWWTo2L2Nu_M125'] ,
                                          'WWTo2L2NuHerwigPS'                : ['WWTo2L2Nu'] ,
                                       },
                            },
               },

  'do_lpTCorrMC'  : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False ,
                  'command'    : 'gardener.py letPtCorrector --isData=0 --cmssw=RPLME_CMSSW'
                } ,

  'do_lpTCorrData'  : {
                  'isChain'    : False ,
                  'do4MC'      : False  ,
                  'do4Data'    : True ,
                  'command'    : 'gardener.py letPtCorrector --isData=1 --cmssw=RPLME_CMSSW'
                } ,


  'do_JESup'    : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False ,
                  'command'    : 'gardener.py JESTreeMaker -k 1 --cmssw=RPLME_CMSSW'
                } ,

  'do_JESdo'    : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False ,
                  'command'    : 'gardener.py JESTreeMaker -k -1 --cmssw=RPLME_CMSSW'
                } ,
  
  'do_JESMaxup'    : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False ,
                  'command'    : 'gardener.py JESTreeMaker -k 1 --maxUncertainty --cmssw=RPLME_CMSSW'
                } ,

  'do_JESMaxdo'    : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False ,
                  'command'    : 'gardener.py JESTreeMaker -k -1 --maxUncertainty --cmssw=RPLME_CMSSW'
                } ,

 
  'bPogSF'   :{
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False  ,
                  #'command'    : 'gardener.py btagPogScaleFactors '
                  # --> switch to multiple bTag algo SF:
                  'command'    : 'gardener.py allBtagPogScaleFactors --cmssw=RPLME_CMSSW'
              },

  'do_LepElepTup'    : {
                   'isChain'    : False ,
                   'do4MC'      : True  ,
                   'do4Data'    : False ,
                   'command'    : 'gardener.py LeppTScalerTreeMaker --lepFlavourToChange ele    -v 1.0 --cmssw=RPLME_CMSSW'
                 } ,
  
  'do_LepElepTdo'    : {
                   'isChain'    : False ,
                   'do4MC'      : True  ,
                   'do4Data'    : False ,
                   'command'    : 'gardener.py LeppTScalerTreeMaker --lepFlavourToChange ele   -v -1.0 --cmssw=RPLME_CMSSW'
                 } ,
  
  'do_LepMupTup'    : {
                   'isChain'    : False ,
                   'do4MC'      : True  ,
                   'do4Data'    : False ,
                   'command'    : 'gardener.py LeppTScalerTreeMaker --lepFlavourToChange mu    -v 1.0 --cmssw=RPLME_CMSSW'
                 } ,
  
  'do_LepMupTdo'    : {
                   'isChain'    : False ,
                   'do4MC'      : True  ,
                   'do4Data'    : False ,
                   'command'    : 'gardener.py LeppTScalerTreeMaker --lepFlavourToChange mu   -v -1.0 --cmssw=RPLME_CMSSW'
                 } ,
  
  'do_METup'        : {
                   'isChain'    : False ,
                   'do4MC'      : True  ,
                   'do4Data'    : False ,
                   'command'    : 'gardener.py metUncertainty --kind=Up --cmssw=RPLME_CMSSW --lepton no   --jetresolution no   --unclustered no  '
                 } ,
  
  'do_METdo'        : {
                   'isChain'    : False ,
                   'do4MC'      : True  ,
                   'do4Data'    : False ,
                   'command'    : 'gardener.py metUncertainty --kind=Dn --cmssw=RPLME_CMSSW --lepton no   --jetresolution no   --unclustered no '
                 } ,

   'dymvaHiggs' :   {
                   'isChain'    : True ,
                   'do4MC'      : True ,
                   'do4Data'    : True ,
                   'subTargets' : ['l2kin','do_dymvaHiggs'],
                  },

   'wwvarfiller' : {
                   'isChain'    : False ,
                   'do4MC'      : True ,
                   'do4Data'    : True ,
                   'command'    : 'gardener.py wwvarfiller',
                   },

   'do_dymvaHiggs'  : {
                   'isChain'    : False ,
                   'do4MC'      : True  ,
                   'do4Data'    : True ,
                   'command'    : 'gardener.py dymvaHiggsFiller',
                },
 
   'Mucca'       :  {
                  'isChain'    : True ,
                  'do4MC'      : True ,
                  'do4Data'    : True ,
                  'subTargets' : [
                                  'do_Mucca_1',
                                  'do_Mucca_2',
                                  'do_Mucca_3',
                                  'do_Mucca_4',
                                  'do_Mucca_5'
                                  ],
                },

   'do_Mucca_1'        : {
                   'isChain'    : False ,
                   'do4MC'      : True  ,
                   'do4Data'    : True ,
                   'command'    : 'gardener.py muccaMvaVarFiller --kind 1'
                 } ,
  
   'do_Mucca_2'        : {
                   'isChain'    : False ,
                   'do4MC'      : True  ,
                   'do4Data'    : True ,
                   'command'    : 'gardener.py muccaMvaVarFiller --kind 2'
                 } ,

   'do_Mucca_3'        : {
                   'isChain'    : False ,
                   'do4MC'      : True  ,
                   'do4Data'    : True ,
                   'command'    : 'gardener.py muccaMvaVarFiller --kind 3'
                 } ,

   'do_Mucca_4'        : {
                   'isChain'    : False ,
                   'do4MC'      : True  ,
                   'do4Data'    : True ,
                   'command'    : 'gardener.py muccaMvaVarFiller --kind 4'
                 } ,

   'do_Mucca_5'        : {
                   'isChain'    : False ,
                   'do4MC'      : True  ,
                   'do4Data'    : True ,
                   'command'    : 'gardener.py muccaMvaVarFiller --kind 5'
                 } ,

  # MUCCA for monoH
  # 2HDM model em
  
  'mucca_2HDMadapt_em'      : {
        'isChain'    : False ,
        'do4MC'      : True  ,
        'do4Data'    : True ,
        'command'    : 'gardener.py muccaMonoHVarFiller --signal="2HDMadapt_em" --training="BDT4" --channel="em" --model="2HDM"'
        }  ,
  
  'mucca_2HDMadaptFull_em'    : {
        'isChain'    : False ,
        'do4MC'      : True  ,
        'do4Data'    : True ,
        'command'    : 'gardener.py muccaMonoHFullVarFiller --signal="2HDMadaptFull_em" --training="BDT4" --channel="em" --model="2HDM"'
        } ,
  
  'mucca_2HDMgrad_em'      : {
        'isChain'    : False ,
        'do4MC'      : True  ,
        'do4Data'    : True ,
        'command'    : 'gardener.py muccaMonoHVarFiller --signal="2HDMgrad_em" --training="BDTG4" --channel="em" --model="2HDM"'
        } ,
  
  'mucca_2HDMgradFull_em'    : {
        'isChain'    : False ,
        'do4MC'      : True  ,
        'do4Data'    : True ,
        'command'    : 'gardener.py muccaMonoHFullVarFiller --signal="2HDMgradFull_em" --training="BDTG4" --channel="em" --model="2HDM"'
        } ,
  
  # 2HDM model sf
  
  'mucca_2HDMadapt_sf'      : {
        'isChain'    : False ,
        'do4MC'      : True  ,
        'do4Data'    : True ,
        'command'    : 'gardener.py muccaMonoHVarFiller --signal="2HDMadapt_sf" --training="BDT4" --channel="sf" --model="2HDM"'
        } ,
  
  'mucca_2HDMadaptFull_sf'    : {
        'isChain'    : False ,
        'do4MC'      : True  ,
        'do4Data'    : True ,
        'command'    : 'gardener.py muccaMonoHFullVarFiller --signal="2HDMadaptFull_sf" --training="BDT4" --channel="sf" --model="2HDM"'
        } ,
  
  'mucca_2HDMgrad_sf'      : {
        'isChain'    : False ,
        'do4MC'      : True  ,
        'do4Data'    : True ,
        'command'    : 'gardener.py muccaMonoHVarFiller --signal="2HDMgrad_sf" --training="BDTG4" --channel="sf" --model="2HDM"'
        } ,
  
  'mucca_2HDMgradFull_sf'    : {
        'isChain'    : False ,
        'do4MC'      : True  ,
        'do4Data'    : True ,
        'command'    : 'gardener.py muccaMonoHFullVarFiller --signal="2HDMgradFull_sf" --training="BDTG4" --channel="sf" --model="2HDM"'
        } ,
  
  # Zbar model em
  
  'mucca_Zbaradapt_em'      : {
        'isChain'    : False ,
        'do4MC'      : True  ,
        'do4Data'    : True ,
        'command'    : 'gardener.py muccaMonoHVarFiller --signal="Zbaradapt_em" --training="BDT4" --channel="em" --model="Zbar"'
        } ,
  
  'mucca_ZbaradaptFull_em'    : {
        'isChain'    : False ,
        'do4MC'      : True  ,
        'do4Data'    : True ,
        'command'    : 'gardener.py muccaMonoHFullVarFiller --signal="ZbaradaptFull_em" --training="BDT4" --channel="em" --model="Zbar"'
        } ,
  
  'mucca_Zbargrad_em'      : {
        'isChain'    : False ,
        'do4MC'      : True  ,
        'do4Data'    : True ,
        'command'    : 'gardener.py muccaMonoHVarFiller --signal="Zbargrad_em" --training="BDTG4" --channel="em" --model="Zbar"'
        } ,
  
  'mucca_ZbargradFull_em'    : {
        'isChain'    : False ,
        'do4MC'      : True  ,
        'do4Data'    : True ,
        'command'    : 'gardener.py muccaMonoHFullVarFiller --signal="ZbargradFull_em" --training="BDTG4" --channel="em" --model="Zbar"'
        } ,
  
  # Zbar model sf
  
  'mucca_Zbaradapt_sf'      : {
        'isChain'    : False ,
        'do4MC'      : True  ,
        'do4Data'    : True ,
        'command'    : 'gardener.py muccaMonoHVarFiller --signal="Zbaradapt_sf" --training="BDT4" --channel="sf" --model="Zbar"'
        } ,
  
  'mucca_ZbaradaptFull_sf'    : {
        'isChain'    : False ,
        'do4MC'      : True  ,
        'do4Data'    : True ,
        'command'    : 'gardener.py muccaMonoHFullVarFiller --signal="ZbaradaptFull_sf" --training="BDT4" --channel="sf" --model="Zbar"'
        } ,
  
  'mucca_Zbargrad_sf'      : {
        'isChain'    : False ,
        'do4MC'      : True  ,
        'do4Data'    : True ,
        'command'    : 'gardener.py muccaMonoHVarFiller --signal="Zbargrad_sf" --training="BDTG4" --channel="sf" --model="Zbar"'
        } ,
  
  'mucca_ZbargradFull_sf'    : {
        'isChain'    : False ,
        'do4MC'      : True  ,
        'do4Data'    : True ,
        'command'    : 'gardener.py muccaMonoHFullVarFiller --signal="ZbargradFull_sf" --training="BDTG4" --channel="sf" --model="Zbar"'
        } ,



   # fix datasets names
  'fixdataset_Herwig_nuisance':  {
               'isChain'    : False ,
               'do4MC'      : True ,
               'do4Data'    : False,
               'onlySample' : [
                               #  qqWW
                               'WWTo2L2NuHerwigPS',
                               # ggH
                               'GluGluHToWWTo2L2NuHerwigPS_M125' ,
                               # VBF
                               'VBFHToWWTo2L2NuHerwigPS_M125'
                              ] ,
               'command'    : 'gardener.py adder -v \'dataset/F=42\'  '
           },

  'fakeSel'  : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : True  ,
                  'command'    : 'gardener.py filter -f \'metPfType1 < 20 && mtw1 < 20\' '
           },


  'topSel'  : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : True  ,
                  'command'    : 'gardener.py filter -f \' mll>20 && std_vector_lepton_pt[0]>25 && std_vector_lepton_pt[1]>15 && njet>1 \' '
           },


  'wwSel'    : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : True  ,
                  'command'    : 'gardener.py filter -f \' mll>12 && std_vector_lepton_pt[0]>20 && std_vector_lepton_pt[1]>10 && std_vector_lepton_pt[2]<10 && metPfType1 > 20 && ptll > 30 && (std_vector_lepton_flavour[0] * std_vector_lepton_flavour[1] == -11*13) \' '
           },

  'monohSel'    : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : True  ,
                  'command'    : 'gardener.py filter --keeplist keeplist.txt '
           },



  'vh3lSel'   : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : True  ,
                  'command'    : 'gardener.py filter -f \' std_vector_lepton_isTightLepton[0] > 0.5  && std_vector_lepton_isTightLepton[1] > 0.5  && std_vector_lepton_isTightLepton[2] > 0.5 && std_vector_lepton_pt[0] > 20. && std_vector_lepton_pt[1] > 10. && std_vector_lepton_pt[2] > 10.\' '
           },

  'vh3lFakeSel'   : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : True  ,
                  'command'    : 'gardener.py filter -f \' std_vector_lepton_pt[0] > 20. && std_vector_lepton_pt[1] > 10. && std_vector_lepton_pt[2] > 10.\' ' ,
           },


  'wh2lss1jSel' :  {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : True  ,
                  'command'    : 'gardener.py filter -f \' mll>12 && std_vector_lepton_pt[0]>25 && std_vector_lepton_pt[1]>10 && std_vector_lepton_pt[2]<10 && std_vector_jet_pt[0]>30 && (abs(std_vector_lepton_flavour[1])==13 || std_vector_lepton_pt[1]>13) && ( (std_vector_lepton_flavour[0] * std_vector_lepton_flavour[1] == 11*11) || (std_vector_lepton_flavour[0] * std_vector_lepton_flavour[1] == 13*13) || (std_vector_lepton_flavour[0] * std_vector_lepton_flavour[1] == 11*13)  )  \' '
           },

  'wh2lss1jDYSel' :  {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : True  ,
                  'command'    : 'gardener.py filter -f \' mll>12 && std_vector_lepton_pt[0]>25 && std_vector_lepton_pt[1]>10 && std_vector_lepton_pt[2]<10 && std_vector_jet_pt[0]>30 && (abs(std_vector_lepton_flavour[1])==13 || std_vector_lepton_pt[1]>13) && ( (abs(std_vector_lepton_flavour[0] * std_vector_lepton_flavour[1]) == 11*11) || (std_vector_lepton_flavour[0] * std_vector_lepton_flavour[1] == 13*13) || (std_vector_lepton_flavour[0] * std_vector_lepton_flavour[1] == 11*13)  )  \' ',
                  'onlySample' : ['DYJetsToLL_M-10to50','DYJetsToLL_M-50'],

           },

  'vbsSel'    : {
                  'isChain'    : True  ,
                  'do4MC'      : True  ,
                  'do4Data'    : True  ,
                  'subTargets' : ['do_vbsSel','chFlipProba2j'], 
           },

  'tightVbsSel'    : {
                  'isChain'    : True  ,
                  'do4MC'      : True  ,
                  'do4Data'    : True  ,
                  'subTargets' : ['do_tightVbsSel','chFlipProba2j'],
           },

  'do_vbsSel'    : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : True  ,
                  'command'    : 'gardener.py filter -f \' std_vector_lepton_pt[0]>20 && std_vector_lepton_pt[1]>10 && std_vector_lepton_pt[2]<10 && njet>=2 && mjj>100 \' '
           },

  'do_tightVbsSel'    : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : True  ,
                  'command'    : 'gardener.py filter -f \' std_vector_lepton_pt[0]>25 && std_vector_lepton_pt[1]>20 && std_vector_lepton_pt[2]<10 && njet>=2 && mjj>500 && detajj > 2.5 \' '
           },

  'chFlipProba' : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False  ,
                  'command'    : 'gardener.py chFlipProba',
                  'onlySample' : ['DYJetsToLL_M-10to50','DYJetsToLL_M-50'],
                  },

  'chFlipProba2j' : {
                  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : False  ,
                  'command'    : 'gardener.py chFlipProba --njets 2',
                  'onlySample' : ['DYJetsToLL_M-10to50','DYJetsToLL_M-50','TTTo2L2Nu_ext1','DYJetsToLL_M-50-LO-ext1','TTTo2L2Nu'],
                  },


  'sfSel'    : {  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : True  ,
                  'command'    : 'gardener.py filter -f \'(metTtrk > 20 || metPfType1 > 20) && ptll > 30 && std_vector_lepton_flavour[0] == -std_vector_lepton_flavour[1] && mll > 12 && std_vector_lepton_pt[0] > 20 && std_vector_lepton_pt[1] > 10 \' ',
           },


  'dymvaSel' : {  'isChain'    : False ,
                  'do4MC'      : True  ,
                  'do4Data'    : True  ,
                  'command'    : 'gardener.py filter -f \' mll > 12 && std_vector_lepton_pt[0]>20 && std_vector_lepton_pt[1]>10 && std_vector_lepton_pt[2]<10 && (abs(std_vector_lepton_flavour[1]) == 13 || (std_vector_lepton_pt[0]>25 && std_vector_lepton_pt[1]>13))  && ( (std_vector_lepton_flavour[0] * std_vector_lepton_flavour[1] == -11*11) || (std_vector_lepton_flavour[0] * std_vector_lepton_flavour[1] == -13*13) ) && std_vector_jet_pt[1]<30    && fabs(91.1876 - mll) > 15  && ( std_vector_jet_pt[0] < 20 || std_vector_jet_cmvav2[0] < -0.5884 ) && ( std_vector_jet_pt[1] < 20 || std_vector_jet_cmvav2[1] < -0.5884 ) && ( std_vector_jet_pt[2] < 20 || std_vector_jet_cmvav2[2] < -0.5884 ) && ( std_vector_jet_pt[3] < 20 || std_vector_jet_cmvav2[3] < -0.5884 ) && ( std_vector_jet_pt[4] < 20 || std_vector_jet_cmvav2[4] < -0.5884 ) && ( std_vector_jet_pt[5] < 20 || std_vector_jet_cmvav2[5] < -0.5884 ) && ( std_vector_jet_pt[6] < 20 || std_vector_jet_cmvav2[6] < -0.5884 ) && ( std_vector_jet_pt[7] < 20 || std_vector_jet_cmvav2[7] < -0.5884 ) && ( std_vector_jet_pt[8] < 20 || std_vector_jet_cmvav2[8] < -0.5884 ) && ( std_vector_jet_pt[9] < 20 || std_vector_jet_cmvav2[9] < -0.5884 ) && metTtrk > 20 \' ' ,

                  'onlySample' : [ 'GluGluHToWWTo2L2Nu_alternative_M125' ,
                                   'DYJetsToLL_M-10to50-LO' ,
                                   'DYJetsToLL_M-50-LO-ext1' ,
                                 ],
          },


  'dymvaSel_2j' : {  'isChain'    : False ,
                     'do4MC'      : True  ,
                     'do4Data'    : True  ,
                     'command'    : 'gardener.py filter -f \' mll > 12 && std_vector_lepton_pt[0]>20 && std_vector_lepton_pt[1]>10 && std_vector_lepton_pt[2]<10 && (abs(std_vector_lepton_flavour[1]) == 13 || (std_vector_lepton_pt[0]>25 && std_vector_lepton_pt[1]>13))  && ( (std_vector_lepton_flavour[0] * std_vector_lepton_flavour[1] == -11*11) || (std_vector_lepton_flavour[0] * std_vector_lepton_flavour[1] == -13*13) ) && (std_vector_jet_pt[0]>=30 && std_vector_jet_pt[1]>=30)    && fabs(91.1876 - mll) > 15  && ( std_vector_jet_pt[0] < 20 || std_vector_jet_cmvav2[0] < -0.5884 ) && ( std_vector_jet_pt[1] < 20 || std_vector_jet_cmvav2[1] < -0.5884 ) && ( std_vector_jet_pt[2] < 20 || std_vector_jet_cmvav2[2] < -0.5884 ) && ( std_vector_jet_pt[3] < 20 || std_vector_jet_cmvav2[3] < -0.5884 ) && ( std_vector_jet_pt[4] < 20 || std_vector_jet_cmvav2[4] < -0.5884 ) && ( std_vector_jet_pt[5] < 20 || std_vector_jet_cmvav2[5] < -0.5884 ) && ( std_vector_jet_pt[6] < 20 || std_vector_jet_cmvav2[6] < -0.5884 ) && ( std_vector_jet_pt[7] < 20 || std_vector_jet_cmvav2[7] < -0.5884 ) && ( std_vector_jet_pt[8] < 20 || std_vector_jet_cmvav2[8] < -0.5884 ) && ( std_vector_jet_pt[9] < 20 || std_vector_jet_cmvav2[9] < -0.5884 ) && metTtrk > 20 \' ' ,

                     'onlySample' : [ 'GluGluHToWWTo2L2Nu_alternative_M125' ,
                                      'VBFHToWWTo2L2NuPowheg_M125' ,
                                      'DYJetsToLL_M-10to50-LO' ,
                                      'DYJetsToLL_M-50-LO-ext1' ,
                                    ],
          },


   'SkimSF'   :  {  'isChain'    : False ,
                    'do4MC'      : True  ,
                    'do4Data'    : True  ,
                    'command'    : 'gardener.py filter -f \' mll > 12 && std_vector_lepton_pt[0]>20 && std_vector_lepton_pt[1]>10 && std_vector_lepton_pt[2]<10 && (abs(std_vector_lepton_flavour[1]) == 13 || (std_vector_lepton_pt[0]>25 && std_vector_lepton_pt[1]>13)) && metTtrk > 20 && (dymvaggh>0.6 || dymvavbf>0.6) \' ' ,
                 },


  'metXYshift' : {  'isChain'    : True ,
                    'do4MC'      : True  ,
                    'do4Data'    : True  ,
                    'subTargets' : ['metXYshift_MC','metXYshift_2016B','metXYshift_2016C','metXYshift_2016D'] ,
                 },

  'metXYshift_MC' : {  'isChain'    : False ,
                       'do4MC'      : True  ,
                       'do4Data'    : False  , 
                       'command'    : 'gardener.py metXYshift --cmssw 763 --paraFile metXYshiftPara_MoriondV1_ftnalPar.txt --sample DY' ,
                    },

  'metXYshift_2016B' : {  'isChain'    : False ,
                       'do4MC'      : False  ,
                       'do4Data'    : True   ,
                       'command'    : 'gardener.py metXYshift --cmssw 763 --paraFile metXYshiftPara_MoriondV1_ftnalPar.txt --sample Run2016BDouble' ,
                       'onlySample' : [ 
                                        'Run2016B_PromptReco_DoubleMuon', 'Run2016B_PromptReco_SingleElectron', 'Run2016B_PromptReco_SingleMuon',
                                        'Run2016B_PromptReco_MuonEG' , 'Run2016B_PromptReco_DoubleEG' 
                                      ],
                    },

  'metXYshift_2016C' : {  'isChain'    : False ,
                       'do4MC'      : False  ,
                       'do4Data'    : True   ,
                       'command'    : 'gardener.py metXYshift --cmssw 763 --paraFile metXYshiftPara_MoriondV1_ftnalPar.txt --sample Run2016CDouble' ,
                       'onlySample' : [
                                        'Run2016C_PromptReco_DoubleMuon', 'Run2016C_PromptReco_SingleElectron', 'Run2016C_PromptReco_SingleMuon',
                                        'Run2016C_PromptReco_MuonEG' , 'Run2016C_PromptReco_DoubleEG'
                                      ],
                    },

  'metXYshift_2016D' : {  'isChain'    : False ,
                       'do4MC'      : False  ,
                       'do4Data'    : True   ,
                       'command'    : 'gardener.py metXYshift --cmssw 763 --paraFile metXYshiftPara_MoriondV1_ftnalPar.txt --sample Run2016DDouble' ,
                       'onlySample' : [
                                        'Run2016D_PromptReco_DoubleMuon', 'Run2016D_PromptReco_SingleElectron', 'Run2016D_PromptReco_SingleMuon',
                                        'Run2016D_PromptReco_MuonEG' , 'Run2016D_PromptReco_DoubleEG'
                                      ],
                    },


}
