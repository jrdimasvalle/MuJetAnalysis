#include <memory>
#include "TTree.h"
#include <iomanip>
#include <sstream>
#include <vector>

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "DataFormats/GeometrySurface/interface/Plane.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include <SimDataFormats/Track/interface/SimTrackContainer.h>
#include <SimDataFormats/Vertex/interface/SimVertexContainer.h>
#include <FWCore/Framework/interface/EventSetupRecord.h>
#include "MagneticField/Engine/interface/MagneticField.h"
#include <DataFormats/DetId/interface/DetId.h>
#include "DataFormats/MuonDetId/interface/DTWireId.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "FWCore/Framework/interface/EventSetupRecordImplementation.h"
#include "FWCore/Framework/interface/eventsetuprecord_registration_macro.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "MuonAnalysis/MuonAssociators/interface/L1MuonMatcherAlgo.h"
#include "MuonAnalysis/MuonAssociators/interface/PropagateToMuon.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include <DataFormats/TrackReco/interface/TrackExtra.h>

#include "GEMCode/GEMValidation/interface/SimTrackMatchManager.h"
#include "GEMCode/GEMValidation/interface/Ptassignment.h"


using namespace std;

struct MyTrackRateL1
{
 void init();
 TTree*book(TTree *t, const std::string & name = "l1_particles_");

 Float_t L1_pt;
 Float_t L1_eta;
 Float_t L1_phi;
 Float_t L1_charge;
};

struct MyTrackRateCSC
{
 void init();
 TTree*book(TTree *t, const std::string & name = "trk_eff_csc_");
 


 Int_t lumi;
 Int_t run;
 Float_t eta_SimTrack;
 Float_t phi_SimTrack;
 Float_t pt_SimTrack;
 Float_t vertex_x;
 Float_t vertex_y;
 Float_t vertex_z;
 Int_t csc_station;
 Int_t csc_ring;
 Int_t csc_chamber;
 Int_t nlayerscsc;

 Int_t endcap_st1;
 Int_t endcap_st2;
 Int_t endcap_st3;
 Int_t endcap_st4;
 Float_t Lxy;
 Float_t pzvz;
 Float_t pp_SimTrack;

 Float_t p_SimTrack;
 Float_t p_c_SimTrack;

 Float_t bending_sh_st1;
 Float_t elliptic_value;
 Float_t elliptic_value_st2;

 Float_t csc_gv_eta;
 Float_t csc_gv_phi;
 Float_t csc_gv_pt;
 Float_t csc_gp_r;
 Float_t csc_gp_x;
 Float_t csc_gp_y;
 Float_t csc_gp_z;
 Float_t csc_gp_eta;
 Float_t csc_gp_phi;
 Float_t csc_deltaphi_gp;
 Float_t csc_deltaphi;
 Int_t has_delta_y;
 Float_t delta_y_23_12;
 Float_t delta_x_23_12;
 Float_t delta_y_24_12;
 Float_t delta_x_24_12;
 Int_t has_delta_y_4;
 Float_t Reco_pT_Position;
 Float_t Reco_pT_Direction;
 Float_t Reco_pT_Direction_Smeared;


 Float_t DeltaPhi_Smeared_03_07;
 Float_t DeltaPhi_Smeared_3_7;
 Float_t DeltaPhi_Smeared_6_14;
 Float_t DeltaPhi_Smeared_30_70;


 Float_t Reco_pT_Direction_Smeared_3_7;
 Float_t Reco_pT_Direction_Smeared_6_14;
 Float_t Reco_pT_Direction_Smeared_9_21;
 Float_t Reco_pT_Direction_Smeared_12_28;
 Float_t Reco_pT_Direction_Smeared_15_35;
 Float_t Reco_pT_Direction_Smeared_18_42;
 Float_t Reco_pT_Direction_Smeared_21_49;
 Float_t Reco_pT_Direction_Smeared_24_56;
 Float_t Reco_pT_Direction_Smeared_27_63;
 Float_t Reco_pT_Direction_Smeared_30_70;
 Float_t Reco_pT_Direction_Smeared_03_07;


 Int_t has_pT_Direction;
 Int_t has_pT_Position;
 Int_t csc_chamber_st4;
 Int_t csc_chamber_st3;
 Int_t csc_chamber_st2;
 Int_t nlayers_st2;
 Int_t nlayers_st3;
 Int_t nlayers_st4;
 Float_t csc_gp_second_st3;
 Float_t csc_gp_second_st2;
 Float_t csc_gp_second_st4;
 Float_t csc_bending_angle_12;
 Float_t csc_bending_angle_13;
 Float_t csc_bending_angle_14;


 Float_t delta_y_gp_12;
 Float_t delta_y_gp_14;


 Float_t delta_x_gp_23;
 Float_t delta_x_gp_24;


 Float_t delta_y_gp_34;
 Float_t delta_x_gp_34;

 Float_t csc_deltaphi_gp_12;
 Float_t csc_deltaphi_gp_13;
 Float_t csc_deltaphi_gp_14;
 Float_t csc_deltaphi_gp_23;
 Float_t csc_deltaphi_gp_24;
 Float_t csc_deltaphi_gp_34;


 Int_t has_csc_12;
 Int_t has_csc_13;
 Int_t has_csc_14;
 Int_t has_csc_23;
 Int_t has_csc_24;
 Int_t has_csc_34;


 Float_t dxy; 
 Float_t charge;
 Int_t ntrks;


 Int_t csc_st1_ring;
 Int_t csc_st1_chamber;
 Int_t csc_st1_nlayerscsc;
 Char_t csc_st1_has_csc_sh; // #layers with SimHits > minHitsChamber    bit1: in odd, bit2: even
 Float_t csc_st1_gv_eta;
 Float_t csc_st1_gv_phi;
 Float_t csc_st1_gv_pt;
 Float_t csc_st1_gp_r;
 Float_t csc_st1_gp_x;
 Float_t csc_st1_gp_y;
 Float_t csc_st1_gp_z;
 Float_t csc_st1_gp_eta;
 Float_t csc_st1_gp_phi;
 Float_t csc_st1_bending_sh;
 Float_t csc_st1_deltaphi;
 
 Int_t csc_st2_ring;
 Int_t csc_st2_chamber;
 Int_t csc_st2_nlayerscsc;
 Char_t csc_st2_has_csc_sh; // #layers with SimHits > minHitsChamber    bit1: in odd, bit2: even
 Float_t csc_st2_gv_eta;
 Float_t csc_st2_gv_phi;
 Float_t csc_st2_gv_pt;
 Float_t csc_st2_gp_r;
 Float_t csc_st2_gp_x;
 Float_t csc_st2_gp_y;
 Float_t csc_st2_gp_z;
 Float_t csc_st2_gp_eta;
 Float_t csc_st2_gp_phi;
 Float_t csc_st2_bending_sh;
 Float_t csc_st2_deltaphi;

 Int_t csc_st3_ring;
 Int_t csc_st3_chamber;
 Int_t csc_st3_nlayerscsc;
 Char_t csc_st3_has_csc_sh; // #layers with SimHits > minHitsChamber    bit1: in odd, bit2: even
 Float_t csc_st3_gv_eta;
 Float_t csc_st3_gv_phi;
 Float_t csc_st3_gv_pt;
 Float_t csc_st3_gp_r;
 Float_t csc_st3_gp_x;
 Float_t csc_st3_gp_y;
 Float_t csc_st3_gp_z;
 Float_t csc_st3_gp_eta;
 Float_t csc_st3_gp_phi;
 Float_t csc_st3_bending_sh;
 Float_t csc_st3_deltaphi;

 Int_t csc_st4_ring;
 Int_t csc_st4_chamber;
 Int_t csc_st4_nlayerscsc;
 Char_t csc_st4_has_csc_sh; // #layers with SimHits > minHitsChamber    bit1: in odd, bit2: even
 Float_t csc_st4_gv_eta;
 Float_t csc_st4_gv_phi;
 Float_t csc_st4_gv_pt;
 Float_t csc_st4_gp_r;
 Float_t csc_st4_gp_x;
 Float_t csc_st4_gp_y;
 Float_t csc_st4_gp_z;
 Float_t csc_st4_gp_eta;
 Float_t csc_st4_gp_phi;
 Float_t csc_st4_bending_sh;
 Float_t csc_st4_deltaphi;

 

 Float_t delta_x_gp_12;
 Float_t delta_x_gp_13;
 Float_t delta_x_gp_14;




 

 Float_t csc_p_over_cosh_eta;

 Float_t csc_deltaeta;
 Float_t csc_deltaeta_14;
 Float_t csc_deltaeta_13;
 Float_t csc_deltaeta_12;

 Int_t npar;
 Float_t pt_position_sh,pt_direction_sh;

};


class DisplacedMuonTriggerRateGENSIM : public edm::EDAnalyzer 
{
public:
  explicit DisplacedMuonTriggerRateGENSIM(const edm::ParameterSet&);
  ~DisplacedMuonTriggerRateGENSIM();
  
  virtual void analyze(const edm::Event&, const edm::EventSetup&) ;
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
private:
  
  void analyzeTrackEfficiency(SimTrackMatchManager& match, int trk_no, int st);

  bool isSimTrackGood(const SimTrack &t);
  int detIdToMEStation(int st, int ri);
  int detIdToMBStation(int wh, int st);
  std::vector<string> dtStations_;
  std::set<int> stationsdt_to_use_;
  std::set<int> stationscsc_to_use_;
  std::set<int> l1particles_muons_;
  
  TTree *tree_eff_dt_[56];
  
  TTree *tree_eff_csc_[3];//sim, position based, direction based
  MyTrackRateCSC etrk_csc_[3];

  TTree *tree_eff_l1_[6];
  MyTrackRateL1 etrk_l1_[6];
  float deltaR;
  int does_it_match;

  double vtx_dt;
  double vty_dt;
  double vtz_dt;
  
  edm::ParameterSet cfg_;
  int verbose_;
  int verboseSimTrack_;
  edm::InputTag simInputLabel_;
  double simTrackMinPt_;
  double simTrackMinEta_;
  double simTrackMaxEta_;
  double simTrackOnlyMuon_;
  std::vector<std::pair<int,int> > dtStationsCo_;
  std::vector<std::pair<int,int> > cscStationsCo_;
};

DisplacedMuonTriggerRateGENSIM::DisplacedMuonTriggerRateGENSIM(const edm::ParameterSet& ps)
  : cfg_(ps.getParameterSet("simTrackMatching"))
  , verbose_(ps.getUntrackedParameter<int>("verbose", 0))
{
  auto simTrack = cfg_.getParameter<edm::ParameterSet>("simTrack");
  verboseSimTrack_ = simTrack.getParameter<int>("verbose");
  simInputLabel_ = edm::InputTag("g4SimHits");
  simTrackMinPt_ = simTrack.getParameter<double>("minPt");
  simTrackMinEta_ = simTrack.getParameter<double>("minEta");
  simTrackMaxEta_ = simTrack.getParameter<double>("maxEta");
  simTrackOnlyMuon_ = simTrack.getParameter<bool>("onlyMuon");



  std::vector<int> CSCStationsToUse;
  CSCStationsToUse.push_back(0);  // ALL
  CSCStationsToUse.push_back(1);  // ME11
  CSCStationsToUse.push_back(2);  // ME11a
  CSCStationsToUse.push_back(3);  // ME11b
  CSCStationsToUse.push_back(4);  // ME12
  CSCStationsToUse.push_back(5);  // ME13
  CSCStationsToUse.push_back(6);  // ME21
  CSCStationsToUse.push_back(7);  // ME22
  CSCStationsToUse.push_back(8);  // ME31
  CSCStationsToUse.push_back(9);  // ME32
  CSCStationsToUse.push_back(10);  // ME41
  CSCStationsToUse.push_back(11);  // ME42
  CSCStationsToUse.push_back(12);  // ME1 only
  std::vector<string> stationsDT; 
  stationsDT.push_back("ALL");
  stationsDT.push_back("MB01");
  stationsDT.push_back("MB11");
  stationsDT.push_back("MB21");
  stationsDT.push_back("MB02");
  stationsDT.push_back("MB12");
  stationsDT.push_back("MB22");
  stationsDT.push_back("MB03");
  stationsDT.push_back("MB13");
  stationsDT.push_back("MB23");
  stationsDT.push_back("MB04");
  stationsDT.push_back("MB14");
  stationsDT.push_back("MB24");
  stationsDT.push_back("MB11n");
  stationsDT.push_back("MB21n");
  stationsDT.push_back("MB12n");
  stationsDT.push_back("MB22n");
  stationsDT.push_back("MB13n");
  stationsDT.push_back("MB23n");
  stationsDT.push_back("MB14n");
  stationsDT.push_back("MB24n");
  stationsDT.push_back("STMB2");  

  std::vector<string> L1Ppabc;
  L1Ppabc.push_back("Muon1");
  L1Ppabc.push_back("Muon2");
  L1Ppabc.push_back("Muon3");
  L1Ppabc.push_back("Muon4");
  L1Ppabc.push_back("Muon5");
  L1Ppabc.push_back("Muon6");
  L1Ppabc.push_back("Muon7");
  L1Ppabc.push_back("Muon8");
  L1Ppabc.push_back("Muon9");
  L1Ppabc.push_back("Muon10");
  
  std::vector<int> DtStationsToUse;
  DtStationsToUse.push_back(0);
  DtStationsToUse.push_back(1);
  DtStationsToUse.push_back(2);
  DtStationsToUse.push_back(3);
  DtStationsToUse.push_back(4);
  DtStationsToUse.push_back(5);
  DtStationsToUse.push_back(6);
  DtStationsToUse.push_back(7);
  DtStationsToUse.push_back(8);
  DtStationsToUse.push_back(9);
  DtStationsToUse.push_back(10);
  DtStationsToUse.push_back(11);
  DtStationsToUse.push_back(12);
  DtStationsToUse.push_back(13);
  DtStationsToUse.push_back(14);
  DtStationsToUse.push_back(15);
  DtStationsToUse.push_back(16);
  DtStationsToUse.push_back(17);
  DtStationsToUse.push_back(18);
  DtStationsToUse.push_back(19);
  DtStationsToUse.push_back(20);
  DtStationsToUse.push_back(21);

  std::vector<int> L1Particles;
  L1Particles.push_back(0);
  L1Particles.push_back(1);
  L1Particles.push_back(2);
  L1Particles.push_back(4);
  L1Particles.push_back(5);
  L1Particles.push_back(6);
  L1Particles.push_back(7);
  L1Particles.push_back(8);
  L1Particles.push_back(9);

  copy(CSCStationsToUse.begin(),CSCStationsToUse.end(), inserter(stationscsc_to_use_, stationscsc_to_use_.end()));
  stringstream ss0;
  ss0<<"trk_rate_csc_sim";
  tree_eff_csc_[0] = etrk_csc_[0].book(tree_eff_csc_[0], ss0.str()); 
  stringstream ss1;
  ss1<< "trk_rate_csc_position";
  tree_eff_csc_[1] = etrk_csc_[1].book(tree_eff_csc_[1], ss1.str());
  stringstream ss2;
  ss2<< "trk_rate_csc_direction";
  tree_eff_csc_[2] = etrk_csc_[2].book(tree_eff_csc_[2], ss2.str());


  dtStationsCo_.push_back(std::make_pair(-99,-99));
  dtStationsCo_.push_back(std::make_pair(0,1));
  dtStationsCo_.push_back(std::make_pair(1,1));
  dtStationsCo_.push_back(std::make_pair(2,1));
  dtStationsCo_.push_back(std::make_pair(0,2));
  dtStationsCo_.push_back(std::make_pair(1,2));
  dtStationsCo_.push_back(std::make_pair(2,2));
  dtStationsCo_.push_back(std::make_pair(0,3));
  dtStationsCo_.push_back(std::make_pair(1,3));
  dtStationsCo_.push_back(std::make_pair(2,3));
  dtStationsCo_.push_back(std::make_pair(0,4));
  dtStationsCo_.push_back(std::make_pair(1,4));
  dtStationsCo_.push_back(std::make_pair(2,4));
  dtStationsCo_.push_back(std::make_pair(-1,1));
  dtStationsCo_.push_back(std::make_pair(-2,1));
  dtStationsCo_.push_back(std::make_pair(-1,2));
  dtStationsCo_.push_back(std::make_pair(-2,2));
  dtStationsCo_.push_back(std::make_pair(-1,3));
  dtStationsCo_.push_back(std::make_pair(-2,3));
  dtStationsCo_.push_back(std::make_pair(-1,4));
  dtStationsCo_.push_back(std::make_pair(-2,4));



  cscStationsCo_.push_back(std::make_pair(-99,-99));
  cscStationsCo_.push_back(std::make_pair(1,-99));
  cscStationsCo_.push_back(std::make_pair(1,4));
  cscStationsCo_.push_back(std::make_pair(1,1));
  cscStationsCo_.push_back(std::make_pair(1,2));
  cscStationsCo_.push_back(std::make_pair(1,3));
  cscStationsCo_.push_back(std::make_pair(2,1));
  cscStationsCo_.push_back(std::make_pair(2,2));
  cscStationsCo_.push_back(std::make_pair(3,1));
  cscStationsCo_.push_back(std::make_pair(3,2));
  cscStationsCo_.push_back(std::make_pair(4,1));
  cscStationsCo_.push_back(std::make_pair(4,2));

};


int DisplacedMuonTriggerRateGENSIM::detIdToMEStation(int st, int ri)
{
  auto p(std::make_pair(st, ri));
  return std::find(cscStationsCo_.begin(), cscStationsCo_.end(), p) - cscStationsCo_.begin();
}

int DisplacedMuonTriggerRateGENSIM::detIdToMBStation(int wh,  int st)
{
  auto p(std::make_pair(wh, st));
  return std::find(dtStationsCo_.begin(), dtStationsCo_.end(),p) - dtStationsCo_.begin();
};

DisplacedMuonTriggerRateGENSIM::~DisplacedMuonTriggerRateGENSIM()
{
}

void
DisplacedMuonTriggerRateGENSIM::analyze(const edm::Event& ev, const edm::EventSetup& es)
{
   using namespace edm;

   edm::Handle<edm::SimTrackContainer> sim_tracks;
   ev.getByLabel(simInputLabel_, sim_tracks);
   const edm::SimTrackContainer & sim_track = *sim_tracks.product();

   edm::Handle<edm::SimVertexContainer> sim_vertices;
   ev.getByLabel(simInputLabel_, sim_vertices);
   const edm::SimVertexContainer & sim_vert = *sim_vertices.product();
 for (unsigned int k=0; k<3; k++){
   etrk_csc_[k].init();// sim
   int trk_no=0;
   for (auto& t: *sim_tracks.product()) {
     if(!isSimTrackGood(t)) continue;

     vtx_dt = sim_vert[t.vertIndex()].position().x();
     vty_dt = sim_vert[t.vertIndex()].position().y();
     vtz_dt = sim_vert[t.vertIndex()].position().z();

     SimTrackMatchManager match(t, sim_vert[t.vertIndex()], cfg_, ev, es);
     analyzeTrackEfficiency(match, trk_no, k);
    trk_no = trk_no + 1;
  }

  etrk_csc_[k].ntrks = trk_no;
  tree_eff_csc_[k]->Fill();
  }

  /*
 *    etrk_csc_[0].init();// sim
 *       int trk1_no=0;
 *          for (auto& t: *sim_tracks.product()) {
 *               if(!isSimTrackGood(t)) continue;
 *
 *                    vtx_dt = sim_vert[t.vertIndex()].position().x();
 *                         vty_dt = sim_vert[t.vertIndex()].position().y();
 *                              vtz_dt = sim_vert[t.vertIndex()].position().z();
 *
 *                                   SimTrackMatchManager match(t, sim_vert[t.vertIndex()], cfg_, ev, es);
 *                                        analyzeTrackEfficiency(match, trk1_no, 1);
 *                                            //a, l1_particles, hlt_l2_pp, l2_track, SegmentsDT);
 *
 *                                                trk1_no = trk1_no + 1;
 *                                                  }
 *                                                    etrk_csc_[1].ntrks = trk_no;
 *                                                      tree_eff_csc_[1]->Fill();
 *                                                         */

}

void 
DisplacedMuonTriggerRateGENSIM::analyzeTrackEfficiency(SimTrackMatchManager& match, int trk_no, int st)
{
  const SimHitMatcher& match_sh = match.simhits();
  const SimTrack &t = match_sh.trk();  
  const CSCRecHitMatcher& match_cscrh = match.cscRecHits();
  const HLTTrackMatcher& match_hlt_track = match.hltTracks();
  auto csc_simhits(match_sh.chamberIdsCSC(0));
  float pt_position_tmp=-99;
  float pt_direction_tmp=-99;

  GlobalPoint gp_sh_odd[4];
  GlobalPoint gp_sh_even[4];
  GlobalVector gv_sh_odd[4];
  GlobalVector gv_sh_even[4];
  bool has_csc_sh[4]={false,false,false,false};
  bool odd[4]={false,false,false,false};
  for(auto d: csc_simhits)
  {
    CSCDetId id(d);
    const int cscst(detIdToMEStation(id.station(),id.ring()));
    if (stationscsc_to_use_.count(cscst) == 0) continue;
    int nlayers(match_sh.nLayersWithHitsInSuperChamber(d));

    if (id.station()==1 and (id.ring()==4 or id.ring()==1)){
    int other_ring(id.ring()==4 ? 1 : 4);
    CSCDetId co_id(id.endcap(), id.station(), other_ring, id.chamber());
      auto rawId(co_id.rawId());
      if (csc_simhits.find(rawId) != csc_simhits.end()) {
        nlayers = nlayers+match_sh.nLayersWithHitsInSuperChamber(rawId);

      }
    }

    if (nlayers < 4) continue;
    
    has_csc_sh[id.station()-1]=true;
    if (id.chamber()%2==1){
	 odd[id.station()-1]=true; 
     	 gp_sh_odd[id.station()-1] = match_sh.simHitsMeanPosition(match_sh.hitsInChamber(d));
    	 gv_sh_odd[id.station()-1] = match_sh.simHitsMeanMomentum(match_sh.hitsInChamber(d));
      }else {
	 odd[id.station()-1]=false; 
     	 gp_sh_even[id.station()-1] = match_sh.simHitsMeanPosition(match_sh.hitsInChamber(d));
    	 gv_sh_even[id.station()-1] = match_sh.simHitsMeanMomentum(match_sh.hitsInChamber(d));
      }
       
  }

  int npar=-1;
  if (has_csc_sh[0] and has_csc_sh[1]){
     GlobalPoint gp1,gp2, gp3;
     GlobalVector gv1,gv2;
     if (odd[0] and not(odd[1]) and not(odd[2])){
        gp1=gp_sh_odd[0];
        gp2=gp_sh_even[1];
        gv1=gv_sh_odd[0];
        gv2=gv_sh_even[1];
	npar=0;
	if (has_csc_sh[2] and not(odd[2])) 
        	gp3=gp_sh_even[2];
     }else if (odd[0] and odd[1]) {
        gp1=gp_sh_odd[0];
        gp2=gp_sh_odd[1];
        gv1=gv_sh_odd[0];
        gv2=gv_sh_odd[1];
	npar=1;
	if (has_csc_sh[2] and odd[2]) 
        	gp3=gp_sh_odd[2];
    }else if (not(odd[0]) and not(odd[1])){
        gp1=gp_sh_even[0];
        gp2=gp_sh_even[1];
        gv1=gv_sh_even[0];
        gv2=gv_sh_even[1];
	npar=2;
	if (has_csc_sh[2] and not(odd[2])) 
        	gp3=gp_sh_even[2];
    }else if (not(odd[0]) and odd[1]){
        gp1=gp_sh_even[0];
        gp2=gp_sh_odd[1];
        gv1=gv_sh_odd[0];
        gv2=gv_sh_odd[1];
	npar=3;
	if (has_csc_sh[2] and odd[2]) 
        	gp3=gp_sh_odd[2];
     }
     float csc_bending_angle_12=deltaPhi(gv1.phi(), gv2.phi());
     if (has_csc_sh[2])
	pt_position_tmp=Ptassign_Position_gp(gp1, gp2, gp3, gp2.eta(), npar); //t.momentum().eta() 

     pt_direction_tmp=Ptassign_Direction(csc_bending_angle_12, gp2.eta(), npar);  
     std::cerr <<"case "<< st <<" eta "<< gp2.eta() <<" has csc sh st12 npar "<< npar <<" simpt "<< t.momentum().pt() <<" pt_position "<< pt_position_tmp << " pt_direction "<< pt_direction_tmp <<std::endl;
  
  } 
  etrk_csc_[st].npar = npar;
  if ((st==0 and etrk_csc_[st].pt_SimTrack>t.momentum().pt()) or (st==1 and etrk_csc_[st].pt_position_sh> pt_position_tmp) or (st==2 and etrk_csc_[st].pt_direction_sh>pt_direction_tmp))
	return;

   etrk_csc_[st].run = match.simhits().event().id().run();
    etrk_csc_[st].lumi = match.simhits().event().id().luminosityBlock();
    etrk_csc_[st].pt_position_sh = pt_position_tmp;
    etrk_csc_[st].pt_direction_sh = pt_direction_tmp;
    etrk_csc_[st].pt_SimTrack = t.momentum().pt();
    etrk_csc_[st].phi_SimTrack = t.momentum().phi();
    etrk_csc_[st].eta_SimTrack = t.momentum().eta();
    etrk_csc_[st].vertex_x = vtx_dt;
    etrk_csc_[st].vertex_y = vty_dt;
    etrk_csc_[st].vertex_z = vtz_dt;
    auto pphi = t.momentum().phi();
    etrk_csc_[st].dxy = vtx_dt*sin(pphi) - vty_dt*cos(pphi);
    etrk_csc_[st].pp_SimTrack= t.momentum().z();
    etrk_csc_[st].Lxy = vtx_dt*cos(pphi) + vty_dt*sin(pphi);
    etrk_csc_[st].pzvz = t.momentum().z()*vtz_dt;
    auto totalp = std::sqrt( t.momentum().x()*t.momentum().x() + t.momentum().y()*t.momentum().y() + t.momentum().z()*t.momentum().z());
    etrk_csc_[st].p_SimTrack = totalp;
    etrk_csc_[st].charge = t.charge();
    etrk_csc_[st].p_c_SimTrack = totalp*t.charge();

   

 for(auto d: csc_simhits)
 {
    CSCDetId id(d);
    const int cscst(detIdToMEStation(id.station(),id.ring()));
    if (stationscsc_to_use_.count(cscst) == 0) continue;
    int nlayers(match_sh.nLayersWithHitsInSuperChamber(d));

    if (id.station()==1 and (id.ring()==4 or id.ring()==1)){
    int other_ring(id.ring()==4 ? 1 : 4);
    CSCDetId co_id(id.endcap(), id.station(), other_ring, id.chamber());
      auto rawId(co_id.rawId());
      if (csc_simhits.find(rawId) != csc_simhits.end()) {
        nlayers = nlayers+match_sh.nLayersWithHitsInSuperChamber(rawId);

      }
    }

    if (nlayers < 4) continue;
    if (id.station() == 1){
    	etrk_csc_[st].csc_st1_nlayerscsc = nlayers;
    	etrk_csc_[st].csc_st1_ring = id.ring();
    	etrk_csc_[st].csc_st1_chamber = id.chamber();

    	GlobalPoint hitGp = match_sh.simHitsMeanPosition(match_sh.hitsInChamber(d));
    	etrk_csc_[st].csc_st1_gp_x = hitGp.x();
    	etrk_csc_[st].csc_st1_gp_y = hitGp.y();
    	etrk_csc_[st].csc_st1_gp_z = hitGp.z();
    	etrk_csc_[st].csc_st1_gp_r = hitGp.perp();
    	etrk_csc_[st].csc_st1_gp_eta = hitGp.eta();
    	etrk_csc_[st].csc_st1_gp_phi = hitGp.phi();
    	etrk_csc_[st].csc_st1_bending_sh = match_sh.LocalBendingInChamber(d);
    
    	const bool odd(id.chamber()%2==1);

    	if (odd) etrk_csc_[st].csc_st1_has_csc_sh |= 1;
    	else etrk_csc_[st].csc_st1_has_csc_sh |= 2;

    	GlobalVector ym = match_sh.simHitsMeanMomentum(match_sh.hitsInChamber(d));
    	etrk_csc_[st].csc_st1_gv_eta = ym.eta();
    	etrk_csc_[st].csc_st1_gv_phi = ym.phi();
    	etrk_csc_[st].csc_st1_gv_pt = ym.perp();
    	etrk_csc_[st].csc_st1_deltaphi = deltaPhi(hitGp.phi(), ym.phi());  //Bending Angle Position and Direction


	float local_bending = match_sh.LocalBendingInChamber(d);
	etrk_csc_[st].bending_sh_st1 = local_bending;

        for(auto s_d: csc_simhits)
        {
             CSCDetId s_id(s_d);
             const int s_st(detIdToMEStation(s_id.station(),s_id.ring()));
             if (stationscsc_to_use_.count(s_st) == 0) continue;
             int d_nlayers(match_sh.nLayersWithHitsInSuperChamber(d));
             int s_nlayers(match_sh.nLayersWithHitsInSuperChamber(s_d));

             if(s_nlayers == 0) continue; // Check to have hits in the secondary chamber
             if(d_nlayers == 0) continue; //Check that has hits in previous one
             if(id.station() == s_id.station()) continue; // no double hits in the same station, ME11 included by default.
             if (s_nlayers < 4) continue;

	     GlobalPoint hitGp2 = match_sh.simHitsMeanPosition(match_sh.hitsInChamber(s_d));
	     GlobalVector ym2 = match_sh.simHitsMeanMomentum(match_sh.hitsInChamber(s_d));

             if(s_id.station()==2){

		float local_bending_st2 = match_sh.LocalBendingInChamber(s_d);
                float anglea =   hitGp2.phi();
                //float newxst1 =  hitGp.x()*cos(anglea) + hitGp.y()*sin(anglea);
                float newyst1 = -hitGp.x()*sin(anglea) + hitGp.y()*cos(anglea);
                //float newxst2 =  hitGp2.x()*cos(anglea) + hitGp2.y()*sin(anglea);
                float newyst2 = -hitGp2.x()*sin(anglea) + hitGp2.y()*cos(anglea);
                float csc_bending_angle_variable = deltaPhi(ym.phi(), ym2.phi());
                float delta_y_gp_12_variable = newyst2 - newyst1;

                for(auto t_d: csc_simhits)
                {
                        CSCDetId t_id(t_d);
                        const int t_st(detIdToMEStation(t_id.station(),t_id.ring()));
                        if (stationscsc_to_use_.count(t_st) == 0) continue;
                        int t_nlayers(match_sh.nLayersWithHitsInSuperChamber(t_d));

                        if (t_nlayers < 4) continue;

                        if (t_id.station()==1) continue;
                        if (t_id.station()==2) continue;
                        if (t_id.station()==4) continue;
                        GlobalPoint hitGp3 = match_sh.simHitsMeanPosition(match_sh.hitsInChamber(t_d));
                        float ysst3 = -hitGp3.x()*sin(anglea) + hitGp3.y()*cos(anglea);
                        //float xsst3 = hitGp3.x()*cos(anglea) + hitGp3.y()*sin(anglea);

                        etrk_csc_[st].has_delta_y = 1;
                        etrk_csc_[1].has_delta_y = 1;

                        etrk_csc_[1].delta_y_23_12 =  ysst3 - newyst2;
                        etrk_csc_[st].delta_y_23_12 =  ysst3 - newyst2;


                        float delta_y_gp_23_variable =  ysst3 - newyst2;


                        int parity_case = -1;

                        if (id.chamber()%2 == 1){
                                if (s_id.chamber()%2==0 and t_id.chamber()%2==0) parity_case = 0;
                                if (s_id.chamber()%2==1 and t_id.chamber()%2==1) parity_case = 1;

                        }
                        if (id.chamber()%2 == 0) {
                                if (s_id.chamber()%2==0 and t_id.chamber()%2==0) parity_case = 2;
                                if (s_id.chamber()%2==1 and t_id.chamber()%2==1) parity_case = 3;

                        }




                        if(parity_case == -1){
                                //std::cout<<" Error in the parity assignment -- Continue"<<std::endl;
                                continue;
                                }


                        float etamin = 0.0;


                        if( fabs(hitGp2.eta())>1.6 and fabs(hitGp2.eta()<1.8)) etamin = 1.6;
                        if( fabs(hitGp2.eta())>1.8 and fabs(hitGp2.eta()<2.0)) etamin = 1.8;
                        if( fabs(hitGp2.eta())>2.0 and fabs(hitGp2.eta()<2.2)) etamin = 2.0;
                        if( fabs(hitGp2.eta())>2.2 and fabs(hitGp2.eta()<2.4)) etamin = 2.2;

                        if(etamin == 0.0){
                                continue;
                        }

                       float slope = 0.0;
                       float intercept = 0.0;
                       float prop = 0.0;
                       float sigma= 0.0;

                       float slope_elliptic_Tao = 0.0;
		       float slope_Tao_st2 = 0.0;
                       if (parity_case ==0 ){
                                slope_elliptic_Tao = -15.71;
				slope_Tao_st2 = -27.76;
                                prop = 0.6484;
                                if (etamin == 1.6){
                                    sigma = 0.01319;
                                    slope = 0.05527;
                                    intercept = 0.08944;
                                }
                                if (etamin == 1.8){
                                    sigma = 0.02154;
                                    slope= 0.08295;
                                    intercept = 0.1279;
                                }
                                if (etamin == 2.0){
                                    sigma = 0.03251;
                                    slope = 0.166;
                                    intercept = 0.2158;
                                }
                                if (etamin == 2.2){
                                    sigma = 0.05515;
                                    slope = 0.4952;
                                    intercept = 0.7103;
                                }
                        }
                        if (parity_case ==1 ){
                                slope_elliptic_Tao = -15.16;
				slope_Tao_st2 = -26.17;
                                prop = 0.3542;
                                if (etamin == 1.6){
                                    sigma = 0.01014;
                                    slope = 0.1067;
                                    intercept = 0.1957;
                                }
                                if (etamin == 1.8){
                                    sigma = 0.01997;
                                    slope= 0.1561;
                                    intercept = 0.2654;
                                }
                                if (etamin == 2.0){
                                   sigma = 0.03619;
                                    slope = 0.3156;
                                    intercept = 0.4514;
                                }
                                if (etamin == 2.2){
                                    sigma = 0.05695;
                                    slope = 0.8242;
                                    intercept = 1.071;
                                }
                        }
                        if (parity_case ==2 ){
                                slope_elliptic_Tao = -13.63;
				slope_Tao_st2 = -24.3;
                                prop = 0.5636;
                                if (etamin == 1.6){
                                    sigma = 0.008583;
                                    slope = 0.05624;
                                    intercept = 0.08417;
                                }
                                if (etamin == 1.8){
                                    sigma = 0.02352;
                                    slope= 0.08702;
                                    intercept = 0.1426;
                                }
                                if (etamin == 2.0){
                                    sigma = 0.03006;
                                    slope = 0.1676;
                                    intercept = 0.2198;
                                }
                                if (etamin == 2.2){
                                    sigma = 0.05692;
                                    slope = 0.4953;
                                    intercept = 0.7272;
                                }
                        }
                        if (parity_case ==3 ){
                                slope_elliptic_Tao = -12.92;
				slope_Tao_st2 = -21.68;
                                prop = 0.3217;
                                if (etamin == 1.6){
                                    sigma = 0.006731;
                                    slope = 0.1066;
                                    intercept = 0.2026;
                                }
                                if (etamin == 1.8){
                                    sigma = 0.02152;
                                    slope= 0.1435;
                                    intercept = 0.2118;
                                }
                                if (etamin == 2.0){
                                    sigma = 0.03513;
                                    slope = 0.2874;
                                    intercept = 0.4055;
                                }
                                if (etamin == 2.2){
                                    sigma = 0.05173;
                                    slope = 0.7625;
                                    intercept = 1.075;
                                }

                        }



                        if (csc_bending_angle_variable == 0) continue;
                        if (delta_y_gp_23_variable == 0 ) continue;
                        if (delta_y_gp_12_variable == 0) continue;
                        if (slope == 0) continue;

			if (sigma == 0) continue;

                        etrk_csc_[1].elliptic_value = abs( deltaPhi(hitGp.phi(), hitGp2.phi()) - slope_elliptic_Tao*local_bending);
			etrk_csc_[1].elliptic_value_st2 = abs(deltaPhi(hitGp.phi(), hitGp2.phi()) - slope_Tao_st2*local_bending_st2);

                        etrk_csc_[1].has_pT_Position = 1;
                        etrk_csc_[st].has_pT_Position = 1;
                        float Reco_pT_Position_var =  (( 1/abs((delta_y_gp_23_variable) - prop*(delta_y_gp_12_variable) )  + intercept )/slope);
                        etrk_csc_[st].Reco_pT_Position  =  Reco_pT_Position_var;
                        etrk_csc_[1].Reco_pT_Position  =  Reco_pT_Position_var;


		}
	    }	
	}
     }

    if (id.station() == 2){
    	etrk_csc_[st].csc_st2_nlayerscsc = nlayers;
    	etrk_csc_[st].csc_st2_ring = id.ring();
    	etrk_csc_[st].csc_st2_chamber = id.chamber();

    	GlobalPoint hitGp = match_sh.simHitsMeanPosition(match_sh.hitsInChamber(d));
    	etrk_csc_[st].csc_st2_gp_x = hitGp.x();
    	etrk_csc_[st].csc_st2_gp_y = hitGp.y();
    	etrk_csc_[st].csc_st2_gp_z = hitGp.z();
    	etrk_csc_[st].csc_st2_gp_r = hitGp.perp();
    	etrk_csc_[st].csc_st2_gp_eta = hitGp.eta();
    	etrk_csc_[st].csc_st2_gp_phi = hitGp.phi();
    	etrk_csc_[st].csc_st2_bending_sh = match_sh.LocalBendingInChamber(d);
    
    	const bool odd(id.chamber()%2==1);

    	if (odd) etrk_csc_[st].csc_st2_has_csc_sh |= 1;
    	else etrk_csc_[st].csc_st2_has_csc_sh |= 2;

    	GlobalVector ym = match_sh.simHitsMeanMomentum(match_sh.hitsInChamber(d));
    	etrk_csc_[st].csc_st2_gv_eta = ym.eta();
    	etrk_csc_[st].csc_st2_gv_phi = ym.phi();
    	etrk_csc_[st].csc_st2_gv_pt = ym.perp();
    	etrk_csc_[st].csc_st2_deltaphi = deltaPhi(hitGp.phi(), ym.phi());  //Bending Angle Position and Direction

     }

     
    if (id.station() == 3){
    	etrk_csc_[st].csc_st3_nlayerscsc = nlayers;
    	etrk_csc_[st].csc_st3_ring = id.ring();
    	etrk_csc_[st].csc_st3_chamber = id.chamber();

    	GlobalPoint hitGp = match_sh.simHitsMeanPosition(match_sh.hitsInChamber(d));
    	etrk_csc_[st].csc_st3_gp_x = hitGp.x();
    	etrk_csc_[st].csc_st3_gp_y = hitGp.y();
    	etrk_csc_[st].csc_st3_gp_z = hitGp.z();
    	etrk_csc_[st].csc_st3_gp_r = hitGp.perp();
    	etrk_csc_[st].csc_st3_gp_eta = hitGp.eta();
    	etrk_csc_[st].csc_st3_gp_phi = hitGp.phi();
    	etrk_csc_[st].csc_st3_bending_sh = match_sh.LocalBendingInChamber(d);
   
    	const bool odd(id.chamber()%2==1);

    	if (odd) etrk_csc_[st].csc_st3_has_csc_sh |= 1;
    	else etrk_csc_[st].csc_st3_has_csc_sh |= 2;

    	GlobalVector ym = match_sh.simHitsMeanMomentum(match_sh.hitsInChamber(d));
    	etrk_csc_[st].csc_st3_gv_eta = ym.eta();
    	etrk_csc_[st].csc_st3_gv_phi = ym.phi();
    	etrk_csc_[st].csc_st3_gv_pt = ym.perp();
    	etrk_csc_[st].csc_st3_deltaphi = deltaPhi(hitGp.phi(), ym.phi());  //Bending Angle Position and Direction

     }

    if (id.station() == 4){
    	etrk_csc_[st].csc_st4_nlayerscsc = nlayers;
    	etrk_csc_[st].csc_st4_ring = id.ring();
    	etrk_csc_[st].csc_st4_chamber = id.chamber();

    	GlobalPoint hitGp = match_sh.simHitsMeanPosition(match_sh.hitsInChamber(d));
    	etrk_csc_[st].csc_st4_gp_x = hitGp.x();
    	etrk_csc_[st].csc_st4_gp_y = hitGp.y();
    	etrk_csc_[st].csc_st4_gp_z = hitGp.z();
    	etrk_csc_[st].csc_st4_gp_r = hitGp.perp();
    	etrk_csc_[st].csc_st4_gp_eta = hitGp.eta();
    	etrk_csc_[st].csc_st4_gp_phi = hitGp.phi();
    	etrk_csc_[st].csc_st4_bending_sh = match_sh.LocalBendingInChamber(d);
   
    	const bool odd(id.chamber()%2==1);

    	if (odd) etrk_csc_[st].csc_st4_has_csc_sh |= 1;
    	else etrk_csc_[st].csc_st4_has_csc_sh |= 2;

    	GlobalVector ym = match_sh.simHitsMeanMomentum(match_sh.hitsInChamber(d));
    	etrk_csc_[st].csc_st4_gv_eta = ym.eta();
    	etrk_csc_[st].csc_st4_gv_phi = ym.phi();
    	etrk_csc_[st].csc_st4_gv_pt = ym.perp();
    	etrk_csc_[st].csc_st4_deltaphi = deltaPhi(hitGp.phi(), ym.phi());  //Bending Angle Position and Direction

     }

 } // End of CSC Sim Hits
}

bool 
DisplacedMuonTriggerRateGENSIM::isSimTrackGood(const SimTrack &t)
{
  if (t.noVertex()) return false;
  if (t.noGenpart()) return false;
  if (std::abs(t.type()) != 13 and simTrackOnlyMuon_) return false;
  if (t.momentum().pt() < simTrackMinPt_) return false;

  return true;
}

void MyTrackRateL1::init()
{

 L1_pt = -99.;
 L1_eta = -9.;
 L1_phi = - 99.;
 L1_charge = - 9.;


} 

void MyTrackRateCSC::init()
{

 lumi = - 99;
 endcap_st1 = 0;
 endcap_st2 = 0;
 endcap_st3 = 0;
 endcap_st4 = 0;
 run = - 99;
 pt_SimTrack = - 9;
 phi_SimTrack = - 9;
 eta_SimTrack= - 9;
 vertex_x = - 9;
 vertex_y = - 9;
 vertex_z = - 9;

 dxy = - 9.;
 p_SimTrack = - 9;
 charge = - 9.;
 p_c_SimTrack = - 9;

 Lxy = - 9999.;
 pzvz = - 999.;
 pp_SimTrack = - 99.;

 bending_sh_st1 = - 99.;
 elliptic_value = - 99.;
 elliptic_value_st2 = - 99.;

 csc_gp_y = - 9999;
 csc_gp_x = - 9999;
 csc_gp_r = - 9999;
 csc_gp_z = - 9999;
 nlayerscsc = 0;
 csc_gp_eta = - 9;
 csc_gp_phi = - 9;
 csc_chamber = - 9;
 csc_chamber_st2 = - 9;
 csc_chamber_st3 = - 9;
 csc_chamber_st4 = - 9;
 csc_gv_eta = - 9.;
 csc_gv_phi = - 9.;
 csc_gv_pt = - 9.;

 delta_x_gp_12 = - 999.;
 delta_x_gp_13 = - 999.;
 delta_x_gp_14 = - 999.;
 delta_y_gp_12 = - 999.;


 delta_y_gp_34 = - 999.;

 csc_station = - 99;
 csc_ring = - 99;
 csc_deltaeta = - 99.;
 csc_deltaeta_14 = - 99.;
 csc_deltaeta_13 = - 99.;
 csc_deltaeta_12 = - 99.;

 has_delta_y = 0;
 delta_y_23_12 = - 9999.;
 delta_x_23_12 = - 9999.;


 Reco_pT_Direction_Smeared_3_7 = - 9.0;
 Reco_pT_Direction_Smeared_6_14 = - 9.0;
 Reco_pT_Direction_Smeared_9_21 = - 9.0;
 Reco_pT_Direction_Smeared_12_28 = - 9.0;
 Reco_pT_Direction_Smeared_15_35 = - 9.0;
 Reco_pT_Direction_Smeared_18_42 = - 9.0;
 Reco_pT_Direction_Smeared_21_49 = - 9.0;
 Reco_pT_Direction_Smeared_24_56 = - 9.0;
 Reco_pT_Direction_Smeared_27_63 = - 9.0;
 Reco_pT_Direction_Smeared_30_70 = - 9.0;
 Reco_pT_Direction_Smeared_03_07 = - 9.0;

 Reco_pT_Direction_Smeared = - 9.0;
 Reco_pT_Direction = - 9.;
 Reco_pT_Position = - 9.;

 DeltaPhi_Smeared_03_07 = - 99.;
 DeltaPhi_Smeared_3_7 = - 99.;
 DeltaPhi_Smeared_6_14 = - 99.;
 DeltaPhi_Smeared_30_70 = - 99.;


 has_pT_Direction = 0;
 has_pT_Position = 0;

 has_delta_y_4 = 0;
 delta_y_24_12 = -9999.;
 delta_x_24_12 = - 9999.;
 delta_x_24_12 = - 9999.;
 delta_y_24_12 = - 9999.;
 has_delta_y_4= 0;

 csc_deltaphi_gp_12=-99.;
 csc_deltaphi_gp_13=-99.;
 csc_deltaphi_gp_14=-99.;
 csc_deltaphi_gp_23=-99.;
 csc_deltaphi_gp_24=-99.;
 csc_deltaphi_gp_34=-99.;

 csc_gp_second_st2 = 0.;
 nlayers_st3 = 0;
 nlayers_st4 = 0;
 nlayers_st2 = 0;
 csc_gp_second_st3 = 0.;
 csc_gp_second_st4 = 0.;
 has_csc_12 = 0;
 has_csc_13 = 0;
 has_csc_14 = 0;
 has_csc_23 = 0;
 has_csc_24 = 0;
 has_csc_34 = 0;

 csc_bending_angle_12 = - 99;
 csc_bending_angle_13 = - 99;
 csc_bending_angle_14 = - 99;
 csc_deltaphi_gp = - 99;
 csc_deltaphi = - 99;

 csc_p_over_cosh_eta = - 99.;


 lumi = - 99;
 run = - 99;
 pt_SimTrack = - 9;
 phi_SimTrack = - 9;
 eta_SimTrack= - 9;
 vertex_x = - 9;
 vertex_y = - 9;
 vertex_z = - 9;

 dxy = - 9.;
 charge = - 9.;
 ntrks=0;
 
 Lxy = - 9999.;
 pzvz = - 999.;
 pp_SimTrack = - 99.;

 csc_st1_chamber = - 9; 
 csc_st1_ring = - 99;
 csc_st1_has_csc_sh = 0;
 csc_st1_gp_y = - 9999;
 csc_st1_gp_x = - 9999;
 csc_st1_gp_r = - 9999;
 csc_st1_gp_z = - 9999;
 csc_st1_nlayerscsc = 0;
 csc_st1_gp_eta = - 9;
 csc_st1_gp_phi = - 9;
 csc_st1_gv_eta = - 9.;
 csc_st1_gv_phi = - 9.;
 csc_st1_gv_pt = - 9.;
 csc_st1_bending_sh = -10;
 csc_st1_deltaphi = - 99;

 csc_st2_chamber = - 9; 
 csc_st2_ring = - 99;
 csc_st2_has_csc_sh = 0;
 csc_st2_gp_y = - 9999;
 csc_st2_gp_x = - 9999;
 csc_st2_gp_r = - 9999;
 csc_st2_gp_z = - 9999;
 csc_st2_nlayerscsc = 0;
 csc_st2_gp_eta = - 9;
 csc_st2_gp_phi = - 9;
 csc_st2_gv_eta = - 9.;
 csc_st2_gv_phi = - 9.;
 csc_st2_gv_pt = - 9.;
 csc_st2_bending_sh = -10;
 csc_st2_deltaphi = - 99;

 csc_st3_chamber = - 9; 
 csc_st3_ring = - 99;
 csc_st3_has_csc_sh = 0;
 csc_st3_gp_y = - 9999;
 csc_st3_gp_x = - 9999;
 csc_st3_gp_r = - 9999;
 csc_st3_gp_z = - 9999;
 csc_st3_nlayerscsc = 0;
 csc_st3_gp_eta = - 9;
 csc_st3_gp_phi = - 9;
 csc_st3_gv_eta = - 9.;
 csc_st3_gv_phi = - 9.;
 csc_st3_gv_pt = - 9.;
 csc_st3_bending_sh = -10;
 csc_st3_deltaphi = - 99;

 csc_st4_chamber = - 9; 
 csc_st4_ring = - 99;
 csc_st4_has_csc_sh = 0;
 csc_st4_gp_y = - 9999;
 csc_st4_gp_x = - 9999;
 csc_st4_gp_r = - 9999;
 csc_st4_gp_z = - 9999;
 csc_st4_nlayerscsc = 0;
 csc_st4_gp_eta = - 9;
 csc_st4_gp_phi = - 9;
 csc_st4_gv_eta = - 9.;
 csc_st4_gv_phi = - 9.;
 csc_st4_gv_pt = - 9.;
 csc_st4_bending_sh = -10;
 csc_st4_deltaphi = - 99;

 delta_x_gp_12 = - 999.;
 delta_x_gp_13 = - 999.;
 delta_x_gp_14 = - 999.;
 delta_y_gp_12 = - 999.;
 delta_y_gp_14 = - 999.;


 delta_y_gp_34 = - 999.;

 csc_deltaeta = - 99.;
 csc_deltaeta_14 = - 99.;
 csc_deltaeta_13 = - 99.;
 csc_deltaeta_12 = - 99.;



 csc_deltaphi_gp_12=-99.;
 csc_deltaphi_gp_13=-99.;
 csc_deltaphi_gp_14=-99.;
 csc_deltaphi_gp_23=-99.;
 csc_deltaphi_gp_24=-99.;
 csc_deltaphi_gp_34=-99.;


 csc_bending_angle_12 = - 99;
 csc_bending_angle_13 = - 99;
 csc_bending_angle_14 = - 99;


 npar = -1;
 pt_position_sh=-99;
 pt_direction_sh = -99;
}

TTree*MyTrackRateL1::book(TTree *t, const std::string & name)
{
  edm::Service< TFileService> fs;
  t = fs->make<TTree>(name.c_str(),name.c_str());
 
  t->Branch("L1_pt", &L1_pt);
  t->Branch("L1_eta", &L1_eta);
  t->Branch("L1_charge", &L1_charge);
  t->Branch("L1_phi", &L1_phi);

  return t;
}

TTree*MyTrackRateCSC::book(TTree *t, const std::string & name)
{

  edm::Service< TFileService > fs;
  t = fs->make<TTree>(name.c_str(),name.c_str());

  t->Branch("run", &run);
  t->Branch("lumi", &lumi);
  t->Branch("pt_SimTrack", &pt_SimTrack);
  t->Branch("phi_SimTrack", &phi_SimTrack);
  t->Branch("eta_SimTrack", &eta_SimTrack);
  t->Branch("vertex_x", &vertex_x);
  t->Branch("vertex_y", &vertex_y);
  t->Branch("vertex_z", &vertex_z);

  t->Branch("dxy", &dxy);
  t->Branch("csc_station", &csc_station);
  t->Branch("csc_ring", &csc_ring);
  t->Branch("Lxy", &Lxy);
  t->Branch("pzvz", &Lxy);
  t->Branch("pp_SimTrack", &pp_SimTrack);


  t->Branch("p_SimTrack", &p_SimTrack);
  t->Branch("charge", &charge);
  t->Branch("p_c_SimTrack", &p_c_SimTrack);
  t->Branch("endcap_st1", &endcap_st1);
  t->Branch("endcap_st2", &endcap_st2);
  t->Branch("endcap_st3", &endcap_st3);
  t->Branch("endcap_st4", &endcap_st4);

  t->Branch("bending_sh_st1", &bending_sh_st1);
  t->Branch("elliptic_value", &elliptic_value);
  t->Branch("elliptic_value_st2", &elliptic_value_st2);
  t->Branch("csc_gp_y", &csc_gp_y);
  t->Branch("csc_gp_x", &csc_gp_x);
  t->Branch("csc_gp_r", &csc_gp_r);
  t->Branch("csc_gp_z", &csc_gp_z);
  t->Branch("csc_gp_eta", &csc_gp_eta);
  t->Branch("csc_gp_phi", &csc_gp_phi);
  t->Branch("nlayerscsc", &nlayerscsc);


  t->Branch("delta_x_gp_12", &delta_x_gp_12);
  t->Branch("delta_y_gp_12", &delta_y_gp_12);
  t->Branch("delta_x_gp_13", &delta_x_gp_13);
  t->Branch("delta_y_gp_14", &delta_y_gp_14);
  t->Branch("delta_x_gp_13", &delta_x_gp_13);
  t->Branch("delta_y_gp_14", &delta_y_gp_14);

  t->Branch("delta_x_gp_24", &delta_x_gp_24);

  t->Branch("delta_y_gp_34", &delta_y_gp_34);

  t->Branch("csc_chamber_st4", &csc_chamber_st4);
  t->Branch("csc_chamber_st3", &csc_chamber_st3);
  t->Branch("csc_chamber_st2", &csc_chamber_st2);
  t->Branch("csc_chamber", &csc_chamber);


  t->Branch("csc_gv_eta", &csc_gv_eta);
  t->Branch("csc_gv_pt", &csc_gv_pt);
  t->Branch("csc_gv_phi", &csc_gv_phi);
  t->Branch("csc_deltaphi", &csc_deltaphi);
  t->Branch("csc_deltaphi_gp", &csc_deltaphi_gp);


  t->Branch("nlayers_st2", &nlayers_st2);
  t->Branch("nlayers_st3", &nlayers_st3);
  t->Branch("nlayers_st4", &nlayers_st4);


  t->Branch("has_delta_y", &has_delta_y);
  t->Branch("delta_x_23_from12", &delta_x_23_12);
  t->Branch("delta_y_23_from12", &delta_y_23_12);

  t->Branch("Reco_pT_Direction", &Reco_pT_Direction);
  t->Branch("Reco_pT_Direction_Smeared", &Reco_pT_Direction_Smeared);
  t->Branch("has_pT_Direction", &has_pT_Direction);
  t->Branch("has_pT_Position", &has_pT_Position);
  t->Branch("Reco_pT_Position", &Reco_pT_Position);
  t->Branch("has_delta_y_4", &has_delta_y);
  t->Branch("delta_x_24_from12", &delta_x_24_12);
  t->Branch("delta_y_24_from12", &delta_y_24_12);


  t->Branch("csc_deltaphi_gp_12", &csc_deltaphi_gp_12);
  t->Branch("csc_deltaphi_gp_13", &csc_deltaphi_gp_13);
  t->Branch("csc_deltaphi_gp_14", &csc_deltaphi_gp_14);
  t->Branch("csc_deltaphi_gp_23", &csc_deltaphi_gp_23);
  t->Branch("csc_deltaphi_gp_24", &csc_deltaphi_gp_24);
  t->Branch("csc_deltaphi_gp_34", &csc_deltaphi_gp_34);



  t->Branch("csc_gp_second_st3", &csc_gp_second_st3);
  t->Branch("csc_gp_second_st2", &csc_gp_second_st2);
  t->Branch("csc_gp_second_st4", &csc_gp_second_st4);
  t->Branch("csc_bending_angle_12", &csc_bending_angle_12);
  t->Branch("csc_bending_angle_13", &csc_bending_angle_13);
  t->Branch("csc_bending_angle_14", &csc_bending_angle_14);

  t->Branch("csc_p_over_cosh_eta", &csc_p_over_cosh_eta);

  t->Branch("csc_deltaeta_14", &csc_deltaeta_14);
  t->Branch("csc_deltaeta", &csc_deltaeta);
  t->Branch("csc_deltaeta_13", &csc_deltaeta_13);
  t->Branch("csc_deltaeta_12", &csc_deltaeta_12);
  t->Branch("has_csc_12", &has_csc_12);
  t->Branch("has_csc_13", &has_csc_13);
  t->Branch("has_csc_14", &has_csc_14);
  t->Branch("has_csc_23", &has_csc_23);
  t->Branch("has_csc_24", &has_csc_24);
  t->Branch("has_csc_34", &has_csc_34);





  t->Branch("run", &run);
  t->Branch("lumi", &lumi);
  t->Branch("pt_SimTrack", &pt_SimTrack);
  t->Branch("phi_SimTrack", &phi_SimTrack);
  t->Branch("eta_SimTrack", &eta_SimTrack); 
  t->Branch("vertex_x", &vertex_x);
  t->Branch("vertex_y", &vertex_y);
  t->Branch("vertex_z", &vertex_z);

  t->Branch("dxy", &dxy);
  t->Branch("Lxy", &Lxy);
  t->Branch("pzvz", &Lxy);
  t->Branch("pp_SimTrack", &pp_SimTrack);


  t->Branch("p_SimTrack", &p_SimTrack);
  t->Branch("charge", &charge);
  t->Branch("p_c_SimTrack", &p_c_SimTrack);
  t->Branch("ntrks", &ntrks);


  t->Branch("csc_st1_ring", &csc_st1_ring);
  t->Branch("csc_st1_chamber", &csc_st1_chamber);
  t->Branch("csc_st1_gp_y", &csc_st1_gp_y);
  t->Branch("csc_st1_gp_x", &csc_st1_gp_x);
  t->Branch("csc_st1_gp_r", &csc_st1_gp_r);
  t->Branch("csc_st1_gp_z", &csc_st1_gp_z);
  t->Branch("csc_st1_gp_eta", &csc_st1_gp_eta);
  t->Branch("csc_st1_gp_phi", &csc_st1_gp_phi);
  t->Branch("csc_st1_nlayerscsc", &csc_st1_nlayerscsc);
  t->Branch("csc_st1_gv_eta", &csc_st1_gv_eta);
  t->Branch("csc_st1_gv_pt", &csc_st1_gv_pt);
  t->Branch("csc_st1_gv_phi", &csc_st1_gv_phi);
  t->Branch("csc_st1_deltaphi", &csc_st1_deltaphi);
  t->Branch("csc_st1_bending_sh", &csc_st1_bending_sh);
  t->Branch("csc_st1_has_csc_sh", &csc_st1_has_csc_sh);

  t->Branch("csc_st2_ring", &csc_st2_ring);
  t->Branch("csc_st2_chamber", &csc_st2_chamber);
  t->Branch("csc_st2_gp_y", &csc_st2_gp_y);
  t->Branch("csc_st2_gp_x", &csc_st2_gp_x);
  t->Branch("csc_st2_gp_r", &csc_st2_gp_r);
  t->Branch("csc_st2_gp_z", &csc_st2_gp_z);
  t->Branch("csc_st2_gp_eta", &csc_st2_gp_eta);
  t->Branch("csc_st2_gp_phi", &csc_st2_gp_phi);
  t->Branch("csc_st2_nlayerscsc", &csc_st2_nlayerscsc);
  t->Branch("csc_st2_gv_eta", &csc_st2_gv_eta);
  t->Branch("csc_st2_gv_pt", &csc_st2_gv_pt);
  t->Branch("csc_st2_gv_phi", &csc_st2_gv_phi);
  t->Branch("csc_st2_deltaphi", &csc_st2_deltaphi);
  t->Branch("csc_st2_bending_sh", &csc_st2_bending_sh);
  t->Branch("csc_st2_has_csc_sh", &csc_st2_has_csc_sh);

  t->Branch("csc_st3_ring", &csc_st3_ring);
  t->Branch("csc_st3_chamber", &csc_st3_chamber);
  t->Branch("csc_st3_gp_y", &csc_st3_gp_y);
  t->Branch("csc_st3_gp_x", &csc_st3_gp_x);
  t->Branch("csc_st3_gp_r", &csc_st3_gp_r);
  t->Branch("csc_st3_gp_z", &csc_st3_gp_z);
  t->Branch("csc_st3_gp_eta", &csc_st3_gp_eta);
  t->Branch("csc_st3_gp_phi", &csc_st3_gp_phi);
  t->Branch("csc_st3_nlayerscsc", &csc_st3_nlayerscsc);
  t->Branch("csc_st3_gv_eta", &csc_st3_gv_eta);
  t->Branch("csc_st3_gv_pt", &csc_st3_gv_pt);
  t->Branch("csc_st3_gv_phi", &csc_st3_gv_phi);
  t->Branch("csc_st3_deltaphi", &csc_st3_deltaphi);
  t->Branch("csc_st3_bending_sh", &csc_st3_bending_sh);
  t->Branch("csc_st3_has_csc_sh", &csc_st3_has_csc_sh);

  t->Branch("csc_st4_ring", &csc_st4_ring);
  t->Branch("csc_st4_chamber", &csc_st4_chamber);
  t->Branch("csc_st4_gp_y", &csc_st4_gp_y);
  t->Branch("csc_st4_gp_x", &csc_st4_gp_x);
  t->Branch("csc_st4_gp_r", &csc_st4_gp_r);
  t->Branch("csc_st4_gp_z", &csc_st4_gp_z);
  t->Branch("csc_st4_gp_eta", &csc_st4_gp_eta);
  t->Branch("csc_st4_gp_phi", &csc_st4_gp_phi);
  t->Branch("csc_st4_nlayerscsc", &csc_st4_nlayerscsc);
  t->Branch("csc_st4_gv_eta", &csc_st4_gv_eta);
  t->Branch("csc_st4_gv_pt", &csc_st4_gv_pt);
  t->Branch("csc_st4_gv_phi", &csc_st4_gv_phi);
  t->Branch("csc_st4_deltaphi", &csc_st4_deltaphi);
  t->Branch("csc_st4_bending_sh", &csc_st4_bending_sh);
  t->Branch("csc_st4_has_csc_sh", &csc_st4_has_csc_sh);


  t->Branch("delta_x_gp_12", &delta_x_gp_12);
  t->Branch("delta_y_gp_12", &delta_y_gp_12);
  t->Branch("delta_x_gp_13", &delta_x_gp_13);
  t->Branch("delta_y_gp_14", &delta_y_gp_14);
  t->Branch("delta_x_gp_13", &delta_x_gp_13);
  t->Branch("delta_y_gp_14", &delta_y_gp_14);

  t->Branch("delta_x_gp_24", &delta_x_gp_24);

  t->Branch("delta_x_gp_34", &delta_x_gp_34);
  t->Branch("delta_y_gp_34", &delta_y_gp_34);

  t->Branch("csc_deltaphi_gp_12", &csc_deltaphi_gp_12);
  t->Branch("csc_deltaphi_gp_13", &csc_deltaphi_gp_13);
  t->Branch("csc_deltaphi_gp_14", &csc_deltaphi_gp_14);
  t->Branch("csc_deltaphi_gp_23", &csc_deltaphi_gp_23);
  t->Branch("csc_deltaphi_gp_24", &csc_deltaphi_gp_24);
  t->Branch("csc_deltaphi_gp_34", &csc_deltaphi_gp_34);



  t->Branch("csc_bending_angle_12", &csc_bending_angle_12);
  t->Branch("csc_bending_angle_13", &csc_bending_angle_13);
  t->Branch("csc_bending_angle_14", &csc_bending_angle_14);

  t->Branch("csc_p_over_cosh_eta", &csc_p_over_cosh_eta);

  t->Branch("csc_deltaeta_14", &csc_deltaeta_14);
  t->Branch("csc_deltaeta", &csc_deltaeta);
  t->Branch("csc_deltaeta_13", &csc_deltaeta_13);
  t->Branch("csc_deltaeta_12", &csc_deltaeta_12);

  t->Branch("npar", &npar);
  t->Branch("pt_position_sh", &pt_position_sh);
  t->Branch("pt_direction_sh", &pt_direction_sh);
  return t;
}



void
DisplacedMuonTriggerRateGENSIM::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DisplacedMuonTriggerRateGENSIM);
