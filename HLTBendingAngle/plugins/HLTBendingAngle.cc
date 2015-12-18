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
#include "GEMCode/GEMValidation/interface/SimTrackMatchManager.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include <DataFormats/TrackReco/interface/TrackExtra.h>


using namespace std;

struct MyTrackEffL1
{
 void init();
 TTree*book(TTree *t, const std::string & name = "l1_particles_");

 Float_t L1_pt;
 Float_t L1_eta;
 Float_t L1_phi;
 Float_t L1_charge;
};

struct MyTrackEffCSC
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
 Float_t dxy; 
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
 Float_t charge;


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
 

 Float_t delta_x_gp_12;
 Float_t delta_y_gp_12;
 Float_t delta_y_gp_13;
 Float_t delta_y_gp_14;
 Float_t delta_x_gp_13;
 Float_t delta_x_gp_14;


 Float_t delta_x_gp_23;
 Float_t delta_y_gp_23;
 Float_t delta_x_gp_24;
 Float_t delta_y_gp_24;


 Float_t delta_y_gp_34;
 Float_t delta_x_gp_34;

 
 Float_t csc_deltaphi_gp_12;
 Float_t csc_deltaphi_gp_13;
 Float_t csc_deltaphi_gp_14;
 Float_t csc_deltaphi_gp_23;
 Float_t csc_deltaphi_gp_24;
 Float_t csc_deltaphi_gp_34;


 Float_t csc_p_over_cosh_eta;

 Float_t csc_deltaeta;
 Float_t csc_deltaeta_14;
 Float_t csc_deltaeta_13;
 Float_t csc_deltaeta_12;
 Int_t has_csc_12;
 Int_t has_csc_13;
 Int_t has_csc_14;
 Int_t has_csc_23;
 Int_t has_csc_24;
 Int_t has_csc_34;


};

struct MyTrackEffDT
{
 void init();
 TTree*book(TTree *t, const std::string & name = "trk_eff_dt_");
 Int_t lumi;
 Int_t run;
 Int_t event;
 Char_t charge_dt;

 Float_t deltaphi_dt_rpc_gv;
 Float_t deltaphi_dt_rpc_gp;
 Float_t deltaphi_first_second_gv;
 Float_t deltaphi_first_second_gp;
 Float_t deltaphi_first_third_gv;
 Float_t deltaphi_first_third_gp;
 Float_t deltaphi_first_fourth_gv;
 Float_t deltaphi_first_fourth_gp;

 Float_t wheel_second;
 Float_t eta_gv_second;
 Float_t phi_gv_second;
 Float_t eta_gp_second;
 Float_t phi_gp_second;

 Float_t wheel_third;
 Float_t eta_gv_third;
 Float_t phi_gv_third;
 Float_t phi_gp_third;
 Float_t eta_gp_third;

 Float_t wheel_fourth;
 Float_t eta_gv_fourth;
 Float_t phi_gv_fourth;
 Float_t eta_gp_fourth;
 Float_t phi_gp_fourth;

 Char_t has_second_dtst_hit;
 Char_t has_third_dtst_hit;
 Char_t has_fourth_dtst_hit;

 Float_t pt_calculated_dt;
 Float_t pt_calculated_dt_12;
 Float_t pt_calculated_dt_14;
 Float_t pt_calculated_dt_13;

 Double_t dtvertex_x;
 Double_t dtvertex_y;
 Double_t dtvertex_z;
 Double_t dtvertex_r;
 Float_t dt_dxy;
 Float_t eta_gp;
 Float_t x_gp;
 Float_t y_gp;
 Float_t z_gp;
 Float_t deltaphi_h_g;
 Float_t deltaphi_t_h;
 Float_t deltaphi_t_g;
 Float_t pt_gv;
 Float_t apt_SimTrack_dt;
 Float_t phi_gv;
 Float_t eta_gv;
 Float_t r_gp;
 Float_t phi_gp;

 Float_t pt_SimTrack_dt;
 Float_t eta_SimTrack_dt;
 Float_t phi_SimTrack_dt;
 Char_t has_dt_sh;
 Float_t R_gv;
 Int_t nlayerdt;
 Int_t nslayerdt;
 Float_t Z_gv;
 Float_t X_gv;
 Float_t Y_gv;
 Int_t wheel;
 Int_t station;




};


class HLTBendingAngle : public edm::EDAnalyzer 
{
public:
  explicit HLTBendingAngle(const edm::ParameterSet&);
  ~HLTBendingAngle();
  
  virtual void analyze(const edm::Event&, const edm::EventSetup&) ;
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
private:
  
  void analyzeTrackEfficiency(SimTrackMatchManager& match, int trk_no);
  //, edm::Handle<std::vector<l1extra::L1MuonParticle> > l1p, edm::Handle<std::vector<reco::RecoChargedCandidate> > hlt_l2_pp, edm::Handle<std::vector<reco::TrackExtra> > l2_track, edm::Handle<edm::RangeMap<DTChamberId,edm::OwnVector<DTRecSegment4D,edm::ClonePolicy<DTRecSegment4D> >,edm::ClonePolicy<DTRecSegment4D> > > SegmentsDT);

  bool isSimTrackGood(const SimTrack &t);
  int detIdToMEStation(int st, int ri);
  int detIdToMBStation(int wh, int st);
  std::vector<string> dtStations_;
  std::set<int> stationsdt_to_use_;
  std::set<int> stationscsc_to_use_;
  std::set<int> l1particles_muons_;
  
  TTree *tree_eff_dt_[56];
  MyTrackEffDT etrk_dt_[56];
  
  TTree *tree_eff_csc_[56];
  MyTrackEffCSC etrk_csc_[56];

  TTree *tree_eff_l1_[6];
  MyTrackEffL1 etrk_l1_[6];
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

HLTBendingAngle::HLTBendingAngle(const edm::ParameterSet& ps)
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


  std::vector<string> stationsCSC;
  stationsCSC.push_back("ALL");
  stationsCSC.push_back("ME11");
  stationsCSC.push_back("ME11a");
  stationsCSC.push_back("ME11b");
  stationsCSC.push_back("ME12");
  stationsCSC.push_back("ME13");
  stationsCSC.push_back("ME21");
  stationsCSC.push_back("ME22");
  stationsCSC.push_back("ME31");
  stationsCSC.push_back("ME32");
  stationsCSC.push_back("ME41");
  stationsCSC.push_back("ME42");
  stationsCSC.push_back("ME1");

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

  // auto input = cms.InputTag("g4SimHits","MuonDTHits");
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

  //copy(L1Particles.begin(), L1Particles.end(), inserter(l1particles_muons_, l1particles_muons_.end()));
  //for(auto m: l1particles_muons_)
  //{
    //stringstream ss;
    //ss<<" trk_eff_l1_"<< L1Ppabc[m];
    //tree_eff_l1_[m] = etrk_l1_[m].book(tree_eff_l1_[m], ss.str());    
 // }

  //copy(DtStationsToUse.begin(),DtStationsToUse.end(),inserter(stationsdt_to_use_,stationsdt_to_use_.end()));
  //for (auto m: stationsdt_to_use_)
  //{
   // stringstream ss;
    //ss<< "trk_eff_dt_" << stationsDT[m];
    //tree_eff_dt_[m] = etrk_dt_[m].book(tree_eff_dt_[m], ss.str());    
  //}

  copy(CSCStationsToUse.begin(),CSCStationsToUse.end(), inserter(stationscsc_to_use_, stationscsc_to_use_.end()));
  for (auto m: stationscsc_to_use_)
  {
    stringstream ss;
    ss<< "trk_eff_csc_" << stationsCSC[m];
    tree_eff_csc_[m] = etrk_csc_[m].book(tree_eff_csc_[m], ss.str());

  }



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


int HLTBendingAngle::detIdToMEStation(int st, int ri)
{
  auto p(std::make_pair(st, ri));
  return std::find(cscStationsCo_.begin(), cscStationsCo_.end(), p) - cscStationsCo_.begin();
}

int HLTBendingAngle::detIdToMBStation(int wh,  int st)
{
  auto p(std::make_pair(wh, st));
  return std::find(dtStationsCo_.begin(), dtStationsCo_.end(),p) - dtStationsCo_.begin();
};

HLTBendingAngle::~HLTBendingAngle()
{
}

void
HLTBendingAngle::analyze(const edm::Event& ev, const edm::EventSetup& es)
{
   using namespace edm;

   edm::Handle<edm::SimTrackContainer> sim_tracks;
   ev.getByLabel(simInputLabel_, sim_tracks);
   const edm::SimTrackContainer & sim_track = *sim_tracks.product();

   edm::Handle<edm::SimVertexContainer> sim_vertices;
   ev.getByLabel(simInputLabel_, sim_vertices);
   const edm::SimVertexContainer & sim_vert = *sim_vertices.product();

   if (verboseSimTrack_){
    // std::cout << "Total number of SimTracks in this event: " << sim_track.size() << std::endl;   
    // std::cout << "Total number of SimVertexs in this event: " << sim_vert.size() << std::endl;
   }
   
   //edm::Handle<std::vector<l1extra::L1MuonParticle> > l1_particles;
   //ev.getByLabel("hltL1extraParticles", l1_particles);


   //edm::Handle<std::vector<reco::TrackExtra> > l2_track;
   //ev.getByLabel("hltL2Muons", l2_track);


   //edm::Handle<std::vector<reco::RecoChargedCandidate> > hlt_l2_pp;
   //ev.getByLabel("hltL2MuonCandidatesNoVtx", hlt_l2_pp);


   //edm::Handle<edm::RangeMap<DTChamberId,edm::OwnVector<DTRecSegment4D,edm::ClonePolicy<DTRecSegment4D> >,edm::ClonePolicy<DTRecSegment4D> > > SegmentsDT;
   //ev.getByLabel("hltDt4DSegments", SegmentsDT);

   int trk_no=0;
   for (auto& t: *sim_tracks.product()) {
     if(!isSimTrackGood(t)) continue;
     if (verboseSimTrack_) {
      // std::cout << "Processing SimTrack " << trk_no + 1 << std::endl;      
      // std::cout << "pt(GeV/c) = " << t.momentum().pt() << ", eta = " << t.momentum().eta()  
      //           << ", phi = " << t.momentum().phi() << ", Q = " << t.charge()
      //           << ", vtxIndex = " << t.vertIndex() << std::endl;
     }

     vtx_dt = sim_vert[t.vertIndex()].position().x();
     vty_dt = sim_vert[t.vertIndex()].position().y();
     vtz_dt = sim_vert[t.vertIndex()].position().z();

     SimTrackMatchManager match(t, sim_vert[t.vertIndex()], cfg_, ev, es);
     analyzeTrackEfficiency(match, trk_no);
    //a, l1_particles, hlt_l2_pp, l2_track, SegmentsDT);

    trk_no = trk_no + 1;
  }
}

void 
HLTBendingAngle::analyzeTrackEfficiency(SimTrackMatchManager& match, int trk_no)
{
  const SimHitMatcher& match_sh = match.simhits();
  //const TrackMatcher& match_track = match.tracks();
  const SimTrack &t = match_sh.trk();
  //const SimVertex &vtx = match_sh.vtx();
  const CSCRecHitMatcher& match_cscrh = match.cscRecHits();
  const HLTTrackMatcher& match_hlt_track = match.hltTracks();
  //const SimVertex& vtx = match_sh.vtx();

  for(auto st: stationscsc_to_use_)
  {
    etrk_csc_[st].init();
    etrk_csc_[st].run = match.simhits().event().id().run();
    etrk_csc_[st].lumi = match.simhits().event().id().luminosityBlock();
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
    etrk_csc_[st]. p_c_SimTrack = totalp*t.charge();
  }


 //CSC SimHits Start here
 auto csc_simhits(match_sh.chamberIdsCSC(0));

 for(auto d: csc_simhits)
 {
    CSCDetId id(d);
    const int st(detIdToMEStation(id.station(),id.ring()));
    if (stationscsc_to_use_.count(st) == 0) continue;
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
    etrk_csc_[st].nlayerscsc = nlayers;
    etrk_csc_[st].csc_station = id.station();
    etrk_csc_[1].csc_station = id.station();
    etrk_csc_[st].csc_chamber = id.chamber();
    etrk_csc_[st].csc_ring = id.ring();
    etrk_csc_[1].csc_ring = id.ring();

    GlobalPoint hitGp = match_sh.simHitsMeanPosition(match_sh.hitsInChamber(d));
    etrk_csc_[st].csc_gp_x = hitGp.x();
    etrk_csc_[st].csc_gp_y = hitGp.y();
    etrk_csc_[st].csc_gp_z = hitGp.z();
    etrk_csc_[st].csc_gp_r = hitGp.perp();
    etrk_csc_[st].csc_gp_eta = hitGp.eta();
    etrk_csc_[st].csc_gp_phi = hitGp.phi();
    
    GlobalVector ym = match_sh.simHitsMeanMomentum(match_sh.hitsInChamber(d));
    etrk_csc_[st].csc_gv_eta = ym.eta();
    etrk_csc_[st].csc_gv_phi = ym.phi();
    etrk_csc_[st].csc_gv_pt = ym.perp();
    etrk_csc_[st].csc_deltaphi = deltaPhi(hitGp.phi(), ym.phi());  //Bending Angle Position and Direction


    // Case ME11
    if(id.station()==1){
        etrk_csc_[1].nlayerscsc = nlayers;
        etrk_csc_[1].csc_gp_x = hitGp.x();
        etrk_csc_[1].csc_gp_y = hitGp.y();
        etrk_csc_[1].csc_gp_z = hitGp.z();
        etrk_csc_[1].csc_gp_r = hitGp.perp();
        etrk_csc_[1].csc_gp_eta = hitGp.eta();
        etrk_csc_[1].csc_gp_phi = hitGp.phi();
        etrk_csc_[1].csc_gv_eta = ym.eta();
        etrk_csc_[1].csc_gv_phi = ym.phi();
        etrk_csc_[1].csc_gv_pt = ym.perp();
        etrk_csc_[1].csc_deltaphi = deltaPhi(hitGp.phi(), ym.phi());  //Bending Angle Position and Direction
	etrk_csc_[1].csc_chamber = id.chamber();
	auto totalp = std::sqrt( t.momentum().x()*t.momentum().x() + t.momentum().y()*t.momentum().y() + t.momentum().z()*t.momentum().z());
	etrk_csc_[1].csc_p_over_cosh_eta = totalp/cosh(std::abs(hitGp.eta()));
	etrk_csc_[1].endcap_st1 = id.endcap();
    }



    //Starting look for second hit
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

        GlobalPoint hitGp2 = match_sh.simHitsMeanPosition(match_sh.hitsInChamber(s_d));
        GlobalVector ym2 = match_sh.simHitsMeanMomentum(match_sh.hitsInChamber(s_d));
 
	if (s_nlayers < 4) continue;


	// Special case for ME1
        if(id.station()==1 ){


            if(s_id.station()==2){

		float anglea =   hitGp2.phi();
            	float newxst1 =  hitGp.x()*cos(anglea) + hitGp.y()*sin(anglea);
	        float newyst1 = -hitGp.x()*sin(anglea) + hitGp.y()*cos(anglea);
		float newxst2 =  hitGp2.x()*cos(anglea) + hitGp2.y()*sin(anglea);
                float newyst2 = -hitGp2.x()*sin(anglea) + hitGp2.y()*cos(anglea);
		
		etrk_csc_[1].delta_x_gp_12 = newxst2 - newxst1;
		etrk_csc_[st].delta_x_gp_12 = newxst2 - newxst1;
		etrk_csc_[st].delta_y_gp_12 = newyst2 - newyst1;
		etrk_csc_[1].delta_y_gp_12 = newyst2 - newyst1;

		etrk_csc_[1].endcap_st2 = s_id.endcap();
		etrk_csc_[st].endcap_st2 = s_id.endcap();

                etrk_csc_[1].csc_bending_angle_12 = deltaPhi(ym.phi(), ym2.phi());
                etrk_csc_[st].csc_bending_angle_12 = deltaPhi(ym.phi(), ym2.phi());
                etrk_csc_[1].has_csc_12 = 1;
                etrk_csc_[st].has_csc_12 = 1;                                                     // Ask for this
                etrk_csc_[st].csc_deltaeta_12 = ym.eta() - ym2.eta();
                etrk_csc_[1].csc_deltaeta_12 = ym.eta() - ym2.eta();
		etrk_csc_[1].csc_gp_second_st2 = hitGp2.eta();
		etrk_csc_[st].csc_gp_second_st2 = hitGp2.eta();

		etrk_csc_[st].nlayers_st2 = s_nlayers;
		etrk_csc_[1].nlayers_st2 = s_nlayers;

		etrk_csc_[st].csc_chamber_st2 = s_id.chamber();
		etrk_csc_[1].csc_chamber_st2 = s_id.chamber();
		etrk_csc_[1].csc_deltaphi_gp_12 = deltaPhi(hitGp.phi(),hitGp2.phi());
		etrk_csc_[st].csc_deltaphi_gp_12 = deltaPhi(hitGp.phi(),hitGp2.phi());




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
			float xsst3 = hitGp2.x()*cos(anglea) + hitGp3.y()*sin(anglea);
			
			etrk_csc_[st].has_delta_y = 1;
			etrk_csc_[1].has_delta_y = 1;

			etrk_csc_[1].delta_y_23_12 =  ysst3 - newyst2;
			etrk_csc_[st].delta_y_23_12 =  ysst3 - newyst2;
			etrk_csc_[st].delta_x_23_12 =  xsst3 - newxst2;
			etrk_csc_[1].delta_x_23_12 =  xsst3 - newxst2;
		}



            }
        
            if(s_id.station()==3){
                etrk_csc_[st].csc_bending_angle_13 = deltaPhi(ym.phi(), ym2.phi());
                etrk_csc_[1].csc_bending_angle_13 = deltaPhi(ym.phi(), ym2.phi());
                etrk_csc_[1].has_csc_13 = 1;
                etrk_csc_[st].has_csc_13 = 1;
                etrk_csc_[1].csc_deltaeta_13 = ym.eta() - ym2.eta();
                etrk_csc_[st].csc_deltaeta_13 = ym.eta() - ym2.eta();
		etrk_csc_[1].csc_gp_second_st3 = hitGp2.eta();
		etrk_csc_[st].csc_gp_second_st3 = hitGp2.eta();
		etrk_csc_[st].nlayers_st3 = s_nlayers;
		etrk_csc_[1].nlayers_st3 = s_nlayers;

		etrk_csc_[1].endcap_st3 = s_id.endcap();
		etrk_csc_[st].endcap_st3 = s_id.endcap();
		etrk_csc_[1].csc_chamber_st3 = s_id.chamber();
		etrk_csc_[st].csc_chamber_st3 = s_id.chamber();
		etrk_csc_[1].csc_deltaphi_gp_13 = deltaPhi(hitGp.phi(),hitGp2.phi());
		etrk_csc_[st].csc_deltaphi_gp_13 = deltaPhi(hitGp.phi(),hitGp2.phi());
                

            }

            if(s_id.station()==4){
                etrk_csc_[st].has_csc_14 = 1;
                etrk_csc_[1].has_csc_14 = 1;
                etrk_csc_[1].csc_bending_angle_14 = deltaPhi(ym.phi(), ym2.phi());
                etrk_csc_[st].csc_bending_angle_14 = deltaPhi(ym.phi(), ym2.phi());
                etrk_csc_[st].csc_deltaeta_14 = ym.eta() - ym2.eta();
                etrk_csc_[1].csc_deltaeta_14 = ym.eta() - ym2.eta();
		etrk_csc_[1].csc_gp_second_st4 = hitGp2.eta();
		etrk_csc_[st].csc_gp_second_st4 = hitGp2.eta();
		etrk_csc_[st].nlayers_st4 = s_nlayers;
		etrk_csc_[1].nlayers_st4 = s_nlayers;
		etrk_csc_[1].endcap_st4 = s_id.endcap();
		etrk_csc_[st].endcap_st4 = s_id.endcap();


		etrk_csc_[st].csc_chamber_st4 = s_id.chamber();
		etrk_csc_[1].csc_chamber_st4 = s_id.chamber();
		etrk_csc_[1].csc_deltaphi_gp_14 = deltaPhi(hitGp.phi(),hitGp2.phi());
		etrk_csc_[st].csc_deltaphi_gp_14 = deltaPhi(hitGp.phi(),hitGp2.phi());


            }
        } // End of especial case for ME11



	if (id.station()==2){

		etrk_csc_[1].endcap_st2 = id.endcap();
		if(s_id.station()==3){

			float angleb =   hitGp.phi();
		        float newxst2 =  hitGp.x()*cos(angleb) + hitGp.y()*sin(angleb);
		        float newyst2 = -hitGp.x()*sin(angleb) + hitGp.y()*cos(angleb);
        		float newxst3 =  hitGp2.x()*cos(angleb) + hitGp2.y()*sin(angleb);
        		float newyst3 = -hitGp2.x()*sin(angleb) + hitGp2.y()*cos(angleb);

			etrk_csc_[1].delta_x_gp_23 = newxst3 - newxst2;
			etrk_csc_[st].delta_x_gp_23 = newxst3 - newxst2;
			etrk_csc_[st].delta_y_gp_23 = newyst3 - newyst2;
			etrk_csc_[1].delta_y_gp_23 = newyst3 - newyst2;

              		etrk_csc_[1].has_csc_23 = 1;
              		etrk_csc_[st].has_csc_23 = 1;
			etrk_csc_[st].csc_deltaphi_gp_23 = deltaPhi(hitGp.phi(),hitGp2.phi());
	              	etrk_csc_[1].csc_deltaphi_gp_23 = deltaPhi(hitGp.phi(),hitGp2.phi());

			etrk_csc_[1].endcap_st3 = s_id.endcap();
			etrk_csc_[st].endcap_st3 = s_id.endcap();


		}

		if(s_id.station()==4){
              		etrk_csc_[1].csc_deltaphi_gp_24 = deltaPhi(hitGp.phi(),hitGp2.phi());
        	      	etrk_csc_[st].csc_deltaphi_gp_24 = deltaPhi(hitGp.phi(),hitGp2.phi());
              		etrk_csc_[1].has_csc_24 = 1;
			etrk_csc_[1].endcap_st4 = s_id.endcap();
			etrk_csc_[st].endcap_st4 = s_id.endcap();
              		etrk_csc_[st].has_csc_24 = 1;
		}
		

	}


	if (id.station()==3){


		etrk_csc_[1].endcap_st3 = id.endcap();

		if (s_id.station()==4){
              	  etrk_csc_[st].csc_deltaphi_gp_34 = deltaPhi(hitGp.phi(),hitGp2.phi());
              	  etrk_csc_[1].csc_deltaphi_gp_34 = deltaPhi(hitGp.phi(),hitGp2.phi());
		  etrk_csc_[1].endcap_st4 = s_id.endcap();
		  etrk_csc_[st].endcap_st4 = s_id.endcap();
		  etrk_csc_[1].has_csc_34 = 1;
		  etrk_csc_[st].has_csc_34 = 1;

		}

	}

        //std::cout<<" Second hit in Station ME"<<s_id.station()<<s_id.ring()<<" with nlayers "<<s_nlayers<<std::endl;

    } // End of Second Hit

 } // End of CSC Sim Hits






 //Filling per station CSC

 for(auto st: stationscsc_to_use_)
 {
  tree_eff_csc_[st]->Fill();
 }

}

bool 
HLTBendingAngle::isSimTrackGood(const SimTrack &t)
{
  // select only muon tracks
  if (t.noVertex()) return false;
  if (t.noGenpart()) return false;
  if (std::abs(t.type()) != 13 and simTrackOnlyMuon_) return false;
  if (t.momentum().pt() < simTrackMinPt_) return false;
  //const float eta(std::abs(t.momentum().eta()));
  //if (eta > simTrackMaxEta_ || eta < simTrackMinEta_) return false; 
  return true;
}

void MyTrackEffL1::init()
{

 L1_pt = -99.;
 L1_eta = -9.;
 L1_phi = - 99.;
 L1_charge = - 9.;


} 

void MyTrackEffCSC::init()
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
 delta_y_gp_13 = - 999.;
 delta_y_gp_14 = - 999.;

 delta_x_gp_23 = - 999.;
 delta_y_gp_23 = - 999.;
 delta_x_gp_24 = - 999.;
 delta_y_gp_24 = - 999.;

 delta_x_gp_34 = - 999.;
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
}
void MyTrackEffDT::init()
{
 lumi = -99;
 run= -99;
 event = -99;

 pt_SimTrack_dt = -9.;
 eta_SimTrack_dt=-9.;
 phi_SimTrack_dt=-9.;
 eta_gp = -9.;
 eta_gv = -9.;
 phi_gv= -9.;
 pt_gv= -9.;
 z_gp = -9900.;
 deltaphi_h_g = -9.;
 apt_SimTrack_dt=-999;
 charge_dt = -99;


 deltaphi_first_second_gv=-99.;
 deltaphi_first_second_gp=-99.;
 deltaphi_first_third_gv=-99.;
 deltaphi_first_third_gp=-99.;
 deltaphi_first_fourth_gv=-99.;
 deltaphi_first_fourth_gp=-99.;
 has_second_dtst_hit=0;
 has_third_dtst_hit=0;
 has_fourth_dtst_hit=0;

 wheel_second = -99;
 phi_gp_second= - 99.;
 eta_gp_second = - 99.;
 phi_gv_second = - 99.;
 eta_gv_second = - 99.;

 wheel_third = -99.;
 phi_gp_third =  - 99.;
 eta_gp_third = -99.;
 phi_gv_third = - 99.;
 eta_gv_third = - 99.;

 wheel_fourth = -99.;
 phi_gp_fourth = - 9999.;
 eta_gp_fourth = -99.;
 phi_gv_fourth = - 9999.;
 eta_gv_fourth = -99.;

 pt_calculated_dt= -9;
 pt_calculated_dt_12=-9;
 pt_calculated_dt_13=-9;
 pt_calculated_dt_14=9;
 x_gp = -9900.;
 y_gp = -9900.;
 r_gp = -9900.;
 phi_gp = -99;
 dt_dxy = -9999;
 dtvertex_x=-9999;
 dtvertex_y=-9999;
 dtvertex_z=-9999;
 dtvertex_r=-9999;
 has_dt_sh= 0;
 
 R_gv=-9999.;
 nlayerdt = 0;
 nslayerdt = 0;
 Z_gv=-9999.;
 X_gv=-9999.;
 Y_gv=-9999.;

 wheel = -9;
 station = - 9;


}

TTree*MyTrackEffL1::book(TTree *t, const std::string & name)
{
  edm::Service< TFileService> fs;
  t = fs->make<TTree>(name.c_str(),name.c_str());
 
  t->Branch("L1_pt", &L1_pt);
  t->Branch("L1_eta", &L1_eta);
  t->Branch("L1_charge", &L1_charge);
  t->Branch("L1_phi", &L1_phi);

  return t;
}

TTree*MyTrackEffCSC::book(TTree *t, const std::string & name)
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

  t->Branch("delta_y_gp_23", &delta_y_gp_23);
  t->Branch("delta_x_gp_23", &delta_x_gp_23);
  t->Branch("delta_y_gp_24", &delta_y_gp_24);
  t->Branch("delta_x_gp_24", &delta_x_gp_24);

  t->Branch("delta_x_gp_34", &delta_x_gp_34);
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
  return t;
}

TTree*MyTrackEffDT::book(TTree *t,const std::string & name)
{
  edm::Service< TFileService > fs;
  t = fs->make<TTree>(name.c_str(),name.c_str());

  t->Branch("lumi", &lumi);
  t->Branch("run", &run);
  t->Branch("event", &event);
  t->Branch("eta_SimTrack_dt", &eta_SimTrack_dt);
  t->Branch("pt_SimTrack_dt", &pt_SimTrack_dt);
  t->Branch("eta_gv", &eta_gv);
  t->Branch("deltaphi_first_second_gv", &deltaphi_first_second_gv);
  t->Branch("deltaphi_first_second_gp", &deltaphi_first_second_gp);
  t->Branch("deltaphi_first_third_gv", &deltaphi_first_third_gv);
  t->Branch("deltaphi_first_third_gp", &deltaphi_first_third_gp);
  t->Branch("deltaphi_first_fourth_gv", &deltaphi_first_fourth_gv);
  t->Branch("deltaphi_first_fourth_gp", &deltaphi_first_fourth_gp);

  t->Branch("has_second_dtst_hit", &has_second_dtst_hit);
  t->Branch("has_third_dtst_hit", &has_third_dtst_hit);
  t->Branch("has_fourth_dtst_hit", &has_fourth_dtst_hit);

  t->Branch("wheel", &wheel);
  t->Branch("station", &station);
  t->Branch("wheel_second", &wheel_second);
  t->Branch("eta_gv_second", &eta_gv_second);
  t->Branch("eta_gp_second", &eta_gp_second);
  t->Branch("phi_gv_second", &phi_gv_second);
  t->Branch("eta_gv_second", &eta_gv_second);

  t->Branch("wheel_third", &wheel_third);
  t->Branch("eta_gv_third", &eta_gv_third);
  t->Branch("eta_gp_third", &eta_gp_third);
  t->Branch("phi_gv_third", &phi_gv_third);
  t->Branch("eta_gv_third", &eta_gv_third);

  t->Branch("wheel_fourth", &wheel_fourth);
  t->Branch("eta_gv_fourth", &eta_gv_fourth);
  t->Branch("eta_gp_fourth", &eta_gp_fourth);
  t->Branch("phi_gv_fourth", &phi_gv_fourth);
  t->Branch("eta_gv_fourth", &eta_gv_fourth);


  t->Branch("pt_calculated_dt", &pt_calculated_dt);
  t->Branch("pt_calculated_dt_12", &pt_calculated_dt_12);
  t->Branch("pt_calculated_dt_13", &pt_calculated_dt_13);
  t->Branch("pt_calculated_dt_14", &pt_calculated_dt_14);


  t->Branch("pt_gv", &pt_gv);
  t->Branch("phi_gv", &phi_gv);
  t->Branch("eta_gp", &eta_gp);
  t->Branch("apt_SimTrack_dt", &apt_SimTrack_dt);
  t->Branch("charge_dt", &charge_dt);
  t->Branch("dtvertex_x", &dtvertex_x);
  t->Branch("dtvertex_y", &dtvertex_y);
  t->Branch("dtvertex_z", &dtvertex_z);
  t->Branch("dtvertex_r", &dtvertex_r);
  t->Branch("deltaphi_h_g", &deltaphi_h_g);
  t->Branch("z_gp", &z_gp);
  t->Branch("x_gp", &x_gp);
  t->Branch("y_gp", &y_gp);
  t->Branch("r_gp", &r_gp);
  t->Branch("phi_gp", &phi_gp);
  t->Branch("phi_SimTrack_dt", &phi_SimTrack_dt);
  t->Branch("has_dt_sh", &has_dt_sh);
  t->Branch("nlayerdt", &nlayerdt);
  t->Branch("nslayerdt", &nslayerdt);
  t->Branch("R_gv", &R_gv);
  t->Branch("Z_gv", &Z_gv);
  t->Branch("X_gv", &X_gv);
  t->Branch("Y_gv", &Y_gv);
  t->Branch("dt_dxy", &dt_dxy);



  return t;


}


void
HLTBendingAngle::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HLTBendingAngle);
