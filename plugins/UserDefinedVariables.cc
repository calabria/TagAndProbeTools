#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "TLorentzVector.h"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include <DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h>
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/GeometryVector/interface/VectorUtil.h"

#include <vector>
#include <utility>
#include <map>

using namespace edm;
using namespace std;

class UserDefinedVariables : public edm::EDProducer {
public:
  explicit UserDefinedVariables(const edm::ParameterSet & iConfig);
  virtual ~UserDefinedVariables() ;
  
  virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);

  virtual void beginJob();

private:

  typedef std::vector<double> vdouble;
  typedef std::vector<std::string> vstring;

  edm::InputTag objects_;
  edm::InputTag met_;         
  edm::InputTag triggerResultsTag_;   
  bool isMC_;
  edm::LumiReWeighting LumiWeights2011A_;
  edm::LumiReWeighting LumiWeights2011B_;
  vdouble TrueDist2011A_f_;
  vdouble TrueDist2011B_f_;
  vdouble MCDist_f_;

};

UserDefinedVariables::UserDefinedVariables(const edm::ParameterSet & iConfig) :
  objects_(iConfig.getParameter<edm::InputTag>("objects")),
  met_(iConfig.getParameter<edm::InputTag>("met")),
  triggerResultsTag_(iConfig.getParameter<edm::InputTag>("triggerResults")), 
  isMC_(iConfig.getParameter<bool>("isMC")),
  TrueDist2011A_f_(iConfig.getParameter<vdouble>("TrueDist2011A")),
  TrueDist2011B_f_(iConfig.getParameter<vdouble>("TrueDist2011B")),
  MCDist_f_(iConfig.getParameter<vdouble>("MCDist"))
{
  produces<edm::ValueMap<float> >("Mt");
  produces<edm::ValueMap<float> >("puMCWeightRun2012");
  produces<edm::ValueMap<float> >("puMCWeightRun2012A");
  produces<edm::ValueMap<float> >("triggerBitSingleMu");
  produces<edm::ValueMap<float> >("muTriggerMatching");
  produces<edm::ValueMap<float> >("mvaMET");
  produces<edm::ValueMap<float> >("PFRelIsoDB04v2");
}


UserDefinedVariables::~UserDefinedVariables()
{

}

void UserDefinedVariables::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {

    // read input
    Handle<View<reco::Candidate> > objects;
    iEvent.getByLabel(objects_, objects);

    edm::Handle<pat::METCollection> metHandle;
    iEvent.getByLabel(met_,metHandle);
    const pat::METCollection* met = metHandle.product();

    int nPUVertices = -99;
    //int nOOTPUVertices = -99;
    float mcPUweight2011A = 1;
    float mcPUweight2011B = 1;
    if(isMC_){

	edm::Handle<std::vector< PileupSummaryInfo > >  PupInfo;
	iEvent.getByLabel(edm::InputTag("addPileupInfo"), PupInfo);

	std::vector<PileupSummaryInfo>::const_iterator PVI;

	for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {

	  	nPUVertices = PVI->getTrueNumInteractions();

	}

    	mcPUweight2011A = LumiWeights2011A_.weight(nPUVertices);
    	mcPUweight2011B = LumiWeights2011B_.weight(nPUVertices);

	//std::cout<<"mcPUweight2011A "<<mcPUweight2011A<<std::endl;

    }

    //cout << mcPUweight << " -- " << nPUVertices << endl;

    // prepare vector for output   
    std::vector<float> values;
    std::vector<float> values2;
    std::vector<float> values3;
    std::vector<float> valuesMVAMet;
    std::vector<float> valuesPFRelIsoDB04v2;

    View<reco::Candidate>::const_iterator object; 
    for (object = objects->begin(); object != objects->end(); ++object) {

      float scalarSumPt = (object->p4()).Pt() + ((*met)[0].p4()).Pt();
      float vectorSumPt = (object->p4() + (*met)[0].p4()).Pt() ;
      float Mt = TMath::Sqrt( scalarSumPt*scalarSumPt - vectorSumPt*vectorSumPt ) ;

      values.push_back(Mt);
      values2.push_back(mcPUweight2011A);
      values3.push_back(mcPUweight2011B);
      valuesMVAMet.push_back(metHandle->front().et());

      const pat::Muon* muon = dynamic_cast<const pat::Muon*>(object->clone());
      if(muon) valuesPFRelIsoDB04v2.push_back(muon->userFloat("PFRelIsoDB04v2"));
      else valuesPFRelIsoDB04v2.push_back(-1);

    }

    // convert into ValueMap and store
    std::auto_ptr<ValueMap<float> > valMap(new ValueMap<float>());
    ValueMap<float>::Filler filler(*valMap);
    filler.insert(objects, values.begin(), values.end());
    filler.fill();
    iEvent.put(valMap, "Mt");

    std::auto_ptr<ValueMap<float> > valMap2(new ValueMap<float>());
    ValueMap<float>::Filler filler2(*valMap2);
    filler2.insert(objects, values2.begin(), values2.end());
    filler2.fill();
    iEvent.put(valMap2, "puMCWeightRun2012");

    std::auto_ptr<ValueMap<float> > valMap3(new ValueMap<float>());
    ValueMap<float>::Filler filler3(*valMap3);
    filler3.insert(objects, values3.begin(), values3.end());
    filler3.fill();
    iEvent.put(valMap3, "puMCWeightRun2012A");

    std::auto_ptr<ValueMap<float> > valMapMVAMet(new ValueMap<float>());
    ValueMap<float>::Filler fillerMVA(*valMapMVAMet);
    fillerMVA.insert(objects, valuesMVAMet.begin(), valuesMVAMet.end());
    fillerMVA.fill();
    iEvent.put(valMapMVAMet, "mvaMET");

    std::auto_ptr<ValueMap<float> > valMapIso(new ValueMap<float>());
    ValueMap<float>::Filler fillerIso(*valMapIso);
    fillerIso.insert(objects, valuesPFRelIsoDB04v2.begin(), valuesPFRelIsoDB04v2.end());
    fillerIso.fill();
    iEvent.put(valMapIso, "PFRelIsoDB04v2");

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    string triggerPathSingleMu;
    string HLTfilterSingleMu;

    std::vector<float> muTriggerMatching_;
    std::vector<float> triggerBitSingleMu_;

    int runNumber = iEvent.id().run();

    if(isMC_){

        triggerPathSingleMu = "HLT_IsoMu24_eta2p1_v13";
    	HLTfilterSingleMu = "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15";

    }

    else{

	if(runNumber >= 190456 && runNumber <= 190738){
        	triggerPathSingleMu = "HLT_IsoMu24_eta2p1_v11";
    		HLTfilterSingleMu = "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoFiltered10";
	}
	else if(runNumber >= 191046 && runNumber <= 193621){
        	triggerPathSingleMu = "HLT_IsoMu24_eta2p1_v12";
    		HLTfilterSingleMu = "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15";
	}
	else if(runNumber >= 193834 && runNumber <= 196531){
        	triggerPathSingleMu = "HLT_IsoMu24_eta2p1_v13";
    		HLTfilterSingleMu = "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15";
	}
	else if(runNumber >= 198049 && runNumber <= 199608){
        	triggerPathSingleMu = "HLT_IsoMu24_eta2p1_v14";
    		HLTfilterSingleMu = "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15";
	}
	else if(runNumber >= 199698 && runNumber <= 208686){
        	triggerPathSingleMu = "HLT_IsoMu24_eta2p1_v15";
    		HLTfilterSingleMu = "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15";
	}

    }

    //std::cout<<"triggerPathSingleMu: "<<triggerPathSingleMu.c_str()<<" HLTfilterSingleMu: "<<HLTfilterSingleMu.c_str()<<std::endl;

    edm::Handle<pat::TriggerEvent> triggerHandle;
    iEvent.getByLabel(triggerResultsTag_, triggerHandle);
    const pat::TriggerEvent* trigger = triggerHandle.product();

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    edm::Handle<edm::TriggerResults> hltResults;
    iEvent.getByLabel(edm::InputTag("TriggerResults::HLT"), hltResults);

    const edm::TriggerNames& triggerNames = iEvent.triggerNames(*hltResults);

    //std::cout<<triggerNames.triggerName(1)<<std::endl; 
    //std::cout << "Available trigger Paths:" << std::endl;
    for ( edm::TriggerNames::Strings::const_iterator triggerName = triggerNames.triggerNames().begin(); triggerName != triggerNames.triggerNames().end(); ++triggerName ) {
	unsigned int index = triggerNames.triggerIndex(*triggerName);
	if ( index < triggerNames.size() ) {
		std::string triggerDecision = ( hltResults->accept(index) ) ? "passed" : "failed";
		//std::cout << " triggerName = " << (*triggerName) << " " << triggerDecision << std::endl;
	}
    }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    edm::Handle<pat::TriggerObjectStandAloneCollection > triggerObjsHandle;
    iEvent.getByLabel(edm::InputTag("patTrigger"),triggerObjsHandle);
    const pat::TriggerObjectStandAloneCollection* triggerObjs = triggerObjsHandle.product();
    

    for (object = objects->begin(); object != objects->end(); ++object) {

	if(trigger){

		const pat::TriggerPath *triggerPathSingleMuon =  trigger->path(triggerPathSingleMu);

		//if(triggerPathSingleMuon) std::cout<<"triggerPathSingleMuon "<<std::endl;

		//std::cout<<"triggerPathSingleMuon "<<triggerPathSingleMuon<<" wasRun "<<triggerPathSingleMuon->wasRun()<<" wasAccept "<<triggerPathSingleMuon->wasAccept()<<" prescale "<<triggerPathSingleMuon->prescale()<<std::endl;

	   	if(triggerPathSingleMuon && triggerPathSingleMuon->wasRun() && triggerPathSingleMuon->wasAccept() && triggerPathSingleMuon->prescale() == 1 ) 
			triggerBitSingleMu_.push_back(1);
	    	else if (triggerPathSingleMuon && triggerPathSingleMuon->wasRun() && triggerPathSingleMuon->wasAccept() && triggerPathSingleMuon->prescale() != 1) 
			triggerBitSingleMu_.push_back(2);
	   	else triggerBitSingleMu_.push_back(0);

	}

	bool matchedSingleMu = false;
	for(pat::TriggerObjectStandAloneCollection::const_iterator it = triggerObjs->begin() ; it !=triggerObjs->end() ; it++){
		pat::TriggerObjectStandAlone *aObj = const_cast<pat::TriggerObjectStandAlone*>(&(*it));

	      	if( Geom::deltaR( aObj->triggerObject().p4(), object->p4() ) < 0.3 && aObj->hasFilterLabel(HLTfilterSingleMu) ){
			matchedSingleMu = true;
	      	}
	}

        //std::cout<<"matching: "<<matchedSingleMu<<std::endl;

	if(matchedSingleMu) muTriggerMatching_.push_back(1);
	else muTriggerMatching_.push_back(0);

    }

    //std::cout<<"Size: "<<objects->size()<<" triggerBitSingleMu: "<<triggerBitSingleMu_.size()<<" muTriggerMatching: "<<muTriggerMatching_.size()<<std::endl;

    std::auto_ptr<ValueMap<float> > valMap7(new ValueMap<float>());
    ValueMap<float>::Filler filler7(*valMap7);
    filler7.insert(objects, triggerBitSingleMu_.begin(), triggerBitSingleMu_.end());
    filler7.fill();
    iEvent.put(valMap7, "triggerBitSingleMu");

    std::auto_ptr<ValueMap<float> > valMap8(new ValueMap<float>());
    ValueMap<float>::Filler filler8(*valMap8);
    filler8.insert(objects, muTriggerMatching_.begin(), muTriggerMatching_.end());
    filler8.fill();
    iEvent.put(valMap8, "muTriggerMatching");

}

void
UserDefinedVariables::beginJob()
{

  if(isMC_){

	  std::vector< float > MCDist ;
	  std::vector< float > TrueDist2011A;
	  std::vector< float > TrueDist2011B;

	  int sizeMCDist_f_ = MCDist_f_.size();

	  for( int i=0; i<sizeMCDist_f_; ++i) {
	      TrueDist2011A.push_back(TrueDist2011A_f_[i]);
	      MCDist.push_back(MCDist_f_[i]);
	  }

	  for( int i=0; i<sizeMCDist_f_; ++i) {
	      TrueDist2011B.push_back(TrueDist2011B_f_[i]);
	  }

	  //std::cout<<MCDist_f_.size()<<" "<<TrueDist2011A_f_.size()<<" "<<TrueDist2011B_f_.size()<<std::endl;
	  //std::cout<<MCDist.size()<<" "<<TrueDist2011A.size()<<" "<<TrueDist2011B.size()<<std::endl;

	  LumiWeights2011A_ = edm::LumiReWeighting(MCDist, TrueDist2011A);
	  LumiWeights2011B_ = edm::LumiReWeighting(MCDist, TrueDist2011B);

  }

}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(UserDefinedVariables);

