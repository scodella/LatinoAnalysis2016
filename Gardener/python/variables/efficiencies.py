import optparse
import numpy
import ROOT
import os.path

from LatinoAnalysis.Gardener.gardening import TreeCloner

#from HWWAnalysis.ShapeAnalysis.triggerEffCombiner import TriggerEff

#    ____________            _____ ____       
#   / __/ _/ _/ /  ___ ___  / __(_) / /__ ____
#  / _// _/ _/ /__/ -_) _ \/ _// / / / -_) __/
# /___/_//_//____/\__/ .__/_/ /_/_/_/\__/_/   
#                   /_/                       

class EffLepFiller(TreeCloner):

    def __init__(self):
        pass

    def __del__(self):
        pass

    def help(self):
        return '''Add a new lepton efficiency weight. The source root files for electrons and muons have to be specified'''

    def addOptions(self,parser):
        description = self.help()
        group = optparse.OptionGroup(parser,self.label, description)
        group.add_option('-s', '--isoid', dest='isoidScaleFactorsFile', help='file with scale factors for isolation and id for leptons', default='data/isoidScaleFactors.py')

        parser.add_option_group(group)
        return group



    def checkOptions(self,opts):
       
        # ~~~~
        isoidScaleFactors = {}
        if opts.isoidScaleFactorsFile == None :
          print " Please provide an input file with the scale factors "
           
        elif os.path.exists(opts.isoidScaleFactorsFile) :
          handle = open(opts.isoidScaleFactorsFile,'r')
          exec(handle)
          handle.close()

        #print " isoidScaleFactors = ", isoidScaleFactors
        
        self.isoidScaleFactors = isoidScaleFactors
        
        self.minpt = 0
        self.maxpt = 1000
        
        self.mineta = 0
        self.maxeta = 4.0
        


    def _getWeight (self, kindLep, pt, eta):

        # fix underflow and overflow
        if pt < self.minpt:
          pt = self.minpt
        if pt > self.maxpt:
          pt = self.maxpt
        
        if eta < self.mineta:
          eta = self.mineta
        if eta > self.maxeta:
          eta = self.maxeta
 
        #print " self.isoidScaleFactors = ", self.isoidScaleFactors
        
        if kindLep in self.isoidScaleFactors.keys() : 
          # get the efficiency
          for point in self.isoidScaleFactors[kindLep] : 
            #   pt           eta           down   value  up
            # (( 0.0, 10.0), (0.0, 1.5), ( 0.980, 0.986, 0.999 ) ),

            if ( pt  >= point[0][0] and pt  < point[0][1] and
                 eta >= point[1][0] and eta < point[1][1] ) :
                return point[2][0], point[2][1], point[2][2]

          # default ... it should never happen!
          # print " default ???"
          return 1.0, 1.0, 1.0
 
        # not a lepton ... like some default value
        return 1.0, 1.0, 1.0
   

    def process(self,**kwargs):
        tree  = kwargs['tree']
        input = kwargs['input']
        output = kwargs['output']

        self.connect(tree,input)

        self.namesOldBranchesToBeModifiedVector = [
           'std_vector_lepton_effW',
           'std_vector_lepton_effW_Up',
           'std_vector_lepton_effW_Down'
           ]
        
        self.clone(output,self.namesOldBranchesToBeModifiedVector)


        bvector_eff =  ROOT.std.vector(float) ()
        self.otree.Branch('std_vector_lepton_effW',bvector_eff)
        bvector_eff_Up =  ROOT.std.vector(float) ()
        self.otree.Branch('std_vector_lepton_effW_Up',bvector_eff_Up)
        bvector_eff_Down =  ROOT.std.vector(float) ()
        self.otree.Branch('std_vector_lepton_effW_Down',bvector_eff_Down)
            

        nentries = self.itree.GetEntries()
        print 'Total number of entries: ',nentries 
                
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

            bvector_eff.clear()
            bvector_eff_Up.clear()
            bvector_eff_Down.clear()

            for iLep in xrange(len(itree.std_vector_lepton_pt)) :
             
              pt = itree.std_vector_lepton_pt [iLep]
              eta = itree.std_vector_lepton_pt [iLep]
              flavour = itree.std_vector_lepton_flavour [iLep]
              
              kindLep = 'lep' # ele or mu
              if abs (flavour) == 11 : 
                kindLep = 'ele'
              elif abs (flavour) == 13 :
                kindLep = 'mu'
 
 
              wdo, w, wup = self._getWeight (kindLep, pt, eta)
             
              bvector_eff.push_back(w)
              bvector_eff_Up.push_back(wup)
              bvector_eff_Down.push_back(wdo)             
              
            otree.Fill()

        self.disconnect()
        print '- Eventloop completed'



#    ________________         _____ ____       
#   / __/ _/ _/_  __/______ _/ __(_) / /__ ____
#  / _// _/ _/ / / / __/ _ `/ _// / / / -_) __/
# /___/_//_/  /_/ /_/  \_, /_/ /_/_/_/\__/_/   
#                     /___/                    

class EffTrgFiller(TreeCloner):

    def __init__(self):
        pass

    def help(self):
        return '''Add a new trigger efficiency weight. The source files must be passed as an option'''

    def addOptions(self,parser):
        description = self.help()
        group = optparse.OptionGroup(parser,self.label, description)

        group.add_option('-s', '--effTrig', dest='efficienciesFile', help='file with trigger efficiencies', default='data/triggerEfficiencies.py')


        parser.add_option_group(group)
        return group 

    def checkOptions(self,opts):
       
        # ~~~~
        triggerEfficiencies = {}
        if opts.efficienciesFile == None :
          print " Please provide an input file with the trigger efficiencies "
           
        elif os.path.exists(opts.efficienciesFile) :
          handle = open(opts.efficienciesFile,'r')
          exec(handle)
          handle.close()

        #print " triggerEfficiencies = ", triggerEfficiencies
        
        self.triggerEfficiencies = triggerEfficiencies
        
        self.minpt = 0
        self.maxpt = 1000
        
        self.mineta = 0
        self.maxeta = 4.0
        


    def _getWeight (self, kindLep, pt, eta):

        # fix underflow and overflow
        if pt < self.minpt:
          pt = self.minpt
        if pt > self.maxpt:
          pt = self.maxpt
        
        if eta < self.mineta:
          eta = self.mineta
        if eta > self.maxeta:
          eta = self.maxeta
 
        #print " self.triggerEfficiencies = ", self.triggerEfficiencies
        
        if kindLep in self.triggerEfficiencies.keys() : 
          # get the efficiency
          for point in self.triggerEfficiencies[kindLep] : 
            #   pt           eta           down   value  up
            # (( 0.0, 10.0), (0.0, 1.5), ( 0.980, 0.986, 0.999 ) ),

            if ( pt  >= point[0][0] and pt  < point[0][1] and
                 eta >= point[1][0] and eta < point[1][1] ) :
                return point[2][0], point[2][1], point[2][2]

          # default ... it should never happen!
          # print " default ???"
          return 1.0, 1.0, 1.0
 
        # not a lepton ... like some default value
        return 1.0, 1.0, 1.0
       
       
    def process(self,**kwargs):
        tree  = kwargs['tree']
        input = kwargs['input']
        output = kwargs['output']


        tree  = kwargs['tree']
        input = kwargs['input']
        output = kwargs['output']

        self.connect(tree,input)

        self.namesOldBranchesToBeModifiedVector = [
           'std_vector_lepton_effTrigW',
           'std_vector_lepton_effTrigW_Up',
           'std_vector_lepton_effTrigW_Down'
           ]
        
        self.clone(output,self.namesOldBranchesToBeModifiedVector)


        bvector_eff =  ROOT.std.vector(float) ()
        self.otree.Branch('std_vector_lepton_effTrigW',bvector_eff)
        bvector_eff_Up =  ROOT.std.vector(float) ()
        self.otree.Branch('std_vector_lepton_effTrigW_Up',bvector_eff_Up)
        bvector_eff_Down =  ROOT.std.vector(float) ()
        self.otree.Branch('std_vector_lepton_effTrigW_Down',bvector_eff_Down)
            

        nentries = self.itree.GetEntries()
        print 'Total number of entries: ',nentries 
                
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

            bvector_eff.clear()
            bvector_eff_Up.clear()
            bvector_eff_Down.clear()

            for iLep in xrange(len(itree.std_vector_lepton_pt)) :
             
              pt = itree.std_vector_lepton_pt [iLep]
              eta = itree.std_vector_lepton_pt [iLep]
              flavour = itree.std_vector_lepton_flavour [iLep]
              
              kindLep = 'lep' # ele or mu
              if abs (flavour) == 11 : 
                kindLep = 'ele'
              elif abs (flavour) == 13 :
                kindLep = 'mu'
 
 
              wdo, w, wup = self._getWeight (kindLep, pt, eta)
             
              bvector_eff.push_back(w)
              bvector_eff_Up.push_back(wup)
              bvector_eff_Down.push_back(wdo)             
              
            otree.Fill()

        self.disconnect()
        print '- Eventloop completed'
