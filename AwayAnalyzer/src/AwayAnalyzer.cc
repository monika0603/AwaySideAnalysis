// Description: [Class to analyse the away-side jet medium modifications]
// Original Author:  Monika Sharma
//         Created:  Wed Mar 18 20:36:37 CET 2015


// system include files
#include <memory>
#include <algorithm>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
//HI stuff
//#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "RecoHI/HiCentralityAlgos/interface/CentralityProvider.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include <DataFormats/TrackReco/interface/Track.h>
#include <DataFormats/TrackReco/interface/TrackFwd.h>
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "RecoJets/JetAlgorithms/interface/JetAlgoHelper.h"

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TGraph.h>
#include "TVector3.h"
#include "assert.h"

using namespace std;
using namespace edm;
using namespace reco;

class AwayAnalyzer : public edm::EDAnalyzer {
   public:
      explicit AwayAnalyzer(const edm::ParameterSet&);
      ~AwayAnalyzer();
      static bool vtxSort( const reco::Vertex &  a, const reco::Vertex & b );
      bool TrackQualityCuts(const reco::Track & track, const reco::Vertex & vertexCollectionSelected);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      void initHistos(const edm::Service<TFileService> & fs);

      // ----------member data ---------------------------
      std::map<std::string,TH1F*> evtPerf_;
      std::map<std::string,TH1F*> trkPerf_;
      std::map<std::string,TH1F*> vtxPerf_;
      std::map<std::string,TH2F*> vtxPerf2D_;
      std::map<std::string,TH1F*> hdNdEtaVzBin_;
      std::map<std::string,TH1F*> hEventVzBin_;

      TH1F* events_;
      TH1F* vertices_;
      TH2D* hSignal;
      TH2D* hBackground;

      int nevt_;
      int ntrack_;
      int nvertex_;
      double pi_;
      int tHighPurityTracks_;
      int nVzBins;
      int nHITracks;

      edm::InputTag vertexSrc_;
      edm::InputTag trackSrc_;
      double etaMin_;
      double etaMax_;
      double ptMin_;
      double vertexZMax_;

      std::string qualityString_;
    
      std::vector<double> etaBins_;
      std::vector<double> vzBins_;
      std::vector<double> NptBins_;
      std::vector<std::vector<TVector3>> pVectVect_trg;
      std::vector<std::vector<TVector3>> pVectVect_ass;
      std::vector<double> zvtxVect;
    
      double cutDzErrMax_;
      double cutDxyErrMax_;
      double cutPtErrMax_;
      double cutMultMin_;
      double cutMultMax_;
    
      double etaMinTrg_;
      double etaMaxTrg_;
      double etaMinAsso_;
      double etaMaxAsso_;
      double ptMinTrg_;
      double ptMaxTrg_;
      double ptMinAsso_;
      double ptMaxAsso_;
    
};

AwayAnalyzer::AwayAnalyzer(const edm::ParameterSet& iConfig):
nevt_(0),
ntrack_(0),
nvertex_(0),
pi_(3.1415927),
vertexSrc_(iConfig.getParameter<edm::InputTag>("vertexSrc")),
trackSrc_(iConfig.getParameter<edm::InputTag>("trackSrc")),
etaMin_(iConfig.getParameter<double>("etaMin")),
etaMax_(iConfig.getParameter<double>("etaMax")),
ptMin_(iConfig.getParameter<double>("ptMin")),
vertexZMax_(iConfig.getParameter<double>("vertexZMax")),
qualityString_(iConfig.getParameter<std::string>("qualityString")),
etaBins_(iConfig.getParameter<std::vector<double> >("etaBins")),
vzBins_(iConfig.getParameter<std::vector<double> >("vzBins"))
{
    edm::Service<TFileService> fs;
    initHistos(fs);
    
    nVzBins = vzBins_.size()-1;
    
    NptBins_ = iConfig.getParameter<std::vector<double> >("NptBins");
    cutMultMin_ = iConfig.getParameter<double>("cutMultMin");
    cutMultMax_ = iConfig.getParameter<double>("cutMultMax");
    cutDzErrMax_ = iConfig.getUntrackedParameter<double>("cutDzErrMax", 3.0);
    cutDxyErrMax_ = iConfig.getUntrackedParameter<double>("cutDxyErrMax", 3.0);
    cutPtErrMax_ = iConfig.getUntrackedParameter<double>("cutPtErrMax", 0.1);
    
    etaMinTrg_ = iConfig.getParameter<double>("etaMinTrg");
    etaMaxTrg_ = iConfig.getParameter<double>("etaMaxTrg");
    etaMinAsso_ = iConfig.getParameter<double>("etaMinAsso");
    etaMaxAsso_ = iConfig.getParameter<double>("etaMaxAsso");
    ptMinTrg_ = iConfig.getParameter<double>("ptMinTrg");
    ptMaxTrg_ = iConfig.getParameter<double>("ptMaxTrg");
    ptMinAsso_ = iConfig.getParameter<double>("ptMinAsso");
    ptMaxAsso_ = iConfig.getParameter<double>("ptMaxAsso");
    
}


AwayAnalyzer::~AwayAnalyzer()
{
}

void
AwayAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   tHighPurityTracks_ = 0;
   vector<TVector3> pVect_trg;
   vector<TVector3> pVect_ass;
    
   Handle<reco::TrackCollection> tracks;
   iEvent.getByLabel(trackSrc_, tracks);

   Handle<std::vector<reco::Vertex> > vertex;
   iEvent.getByLabel(vertexSrc_, vertex);

   std::vector<reco::Vertex> vsorted = *vertex;
   // sort the vertcies by number of tracks in descending order                                                                
   // use chi2 as tiebreaker                                                                                                   
   std::sort( vsorted.begin(), vsorted.end(), AwayAnalyzer::vtxSort );
   // skip events with no PV, this should not happen                                                                           
   if( vsorted.size() == 0) return;
   // skip events failing vertex cut                                                                                           
   if( fabs(vsorted[0].z()) > vertexZMax_ ) return;

   nevt_++;
   events_->Fill(0.5);
   evtPerf_["Nvtx"]->Fill(vsorted.size());
   evtPerf_["Ntrk"]->Fill(tracks->size());

   int lumi = iEvent.getLuminosityBlock().luminosityBlock();
   evtPerf_["Lumi"]->Fill(lumi);
   evtPerf_["NvtxLumi"]->Fill(lumi,vsorted.size());

   int vcount = 0;
   for( const auto & vi :  vsorted )
     {
       vtxPerf_["Ntrk"]->Fill(vi.tracksSize());
       vtxPerf_["x"]->Fill(vi.x());
       vtxPerf_["y"]->Fill(vi.y());
       vtxPerf_["z"]->Fill(vi.z());
       vtxPerf2D_["Ntrk2D"]->Fill(vcount,vi.tracksSize());
       vertices_->Fill(0.5);
       vcount++;
       nvertex_++;
     }

   for (unsigned int i =1; i<vsorted.size(); i++)
     {
       double dz = fabs( vsorted[i].z() - vsorted[0].z() );
       double dx = fabs( vsorted[i].x() - vsorted[0].x() );
       double dy = fabs( vsorted[i].y() - vsorted[0].y() );
       double dxy  = sqrt ( dx*dx + dy*dy );
       vtxPerf_["assocVtxDz"]->Fill(dz);
       vtxPerf2D_["assocVtxDzNtrk"]->Fill(dz,vsorted[i].tracksSize() );
       vtxPerf_["assocVtxDxy"]->Fill(dxy);
       vtxPerf2D_["assocVtxDxyNtrk"]->Fill(dxy,vsorted[i].tracksSize() );
     }

   math::XYZPoint vtxPoint(0.0,0.0,0.0);
   double vzErr =0.0, vxErr=0.0, vyErr=0.0;

   // use vertex w most tracks as primary vertex                                                                               
   if( vsorted.size() != 0 )
     {
       vtxPoint=vsorted.begin()->position();
       vzErr=vsorted.begin()->zError();
       vxErr=vsorted.begin()->xError();
       vyErr=vsorted.begin()->yError();
     }
    
    for( const auto & track : *tracks ) //de-referencing the pointer "tracks" to track and auto will automatically know the type of tracks. This is a new way of looping over the tracks. It is guaranteed to run over all the tracks.
    {
        double dxy=0.0, dz=0.0, dxysigma=0.0, dzsigma=0.0;
        dxy = track.dxy(vtxPoint);
        dz = track.dz(vtxPoint);
        dxysigma = sqrt(track.d0Error()*track.d0Error()+vxErr*vyErr);
        dzsigma = sqrt(track.dzError()*track.dzError()+vzErr*vzErr);
        
        if( !TrackQualityCuts(track, vsorted[0])) continue;
        
        tHighPurityTracks_++;
        
        char histoName1[200];
        char histoName2[200];
        for(int kVz=0; kVz<nVzBins; kVz++) {
            if(vsorted[0].z() > vzBins_[kVz] && vsorted[0].z() <= vzBins_[kVz+1])
            {
                sprintf(histoName1, "hdNdEta_VzBin_%d", kVz);
                hdNdEtaVzBin_[histoName1]->Fill(track.eta());
                sprintf(histoName2, "nEventsVzBin_%d", kVz);
                hEventVzBin_[histoName2]->Fill(0.5);
            }
        }
        
        if( track.eta() <= etaMax_ && track.eta() >= etaMin_ && track.pt() > ptMin_)
        {
            trkPerf_["Nhit"]->Fill(track.numberOfValidHits());
            trkPerf_["pt"]->Fill(track.pt());
            trkPerf_["eta"]->Fill( track.eta() );
            trkPerf_["ptHigh"]->Fill(track.pt());
            trkPerf_["phi"]->Fill(track.phi());
            trkPerf_["dxyErr"]->Fill(dxy/dxysigma);
            trkPerf_["dzErr"]->Fill(dz/dzsigma);
            trkPerf_["chi2"]->Fill(track.normalizedChi2());
            trkPerf_["pterr"]->Fill(track.ptError() / track.pt() );
            ntrack_++;
        }
        
    }
    
    if( !(tHighPurityTracks_ >= cutMultMin_ && tHighPurityTracks_ < cutMultMax_)) return;
    evtPerf_["NHPtrk"]->Fill(tHighPurityTracks_);
    
    for( const auto & track : *tracks)
    {
        if( !TrackQualityCuts(track, vsorted[0])) continue;
        
        TVector3 pvector;
        pvector.SetPtEtaPhi(track.pt(),track.eta(),track.phi());
        
        if(track.eta()<=etaMaxAsso_ && track.eta()>=etaMinAsso_
           && track.pt()<ptMaxAsso_ && track.pt()>ptMinAsso_)
        {
            pVect_ass.push_back(pvector);
            trkPerf_["ptAsso"]->Fill(track.pt());
            trkPerf_["etaAsso"]->Fill(track.eta());
            trkPerf_["phiAsso"]->Fill(track.phi());
        }
        
      
        TVector3 pvector1;
        pvector1.SetPtEtaPhi(track.pt(),track.eta(),track.phi());
            
        if(track.eta()<=etaMaxTrg_ && track.eta()>=etaMinTrg_
            && track.pt()<=ptMaxTrg_ && track.pt()>=ptMinTrg_)
        {
            pVect_trg.push_back(pvector);
            trkPerf_["ptTrg"]->Fill(track.pt());
            trkPerf_["etaTrg"]->Fill(track.eta());
            trkPerf_["phiTrg"]->Fill(track.phi());
        }
    }
    
    /////// Calculating the signal for di-hadron correlations ////////////////
    int nMultTrg = (int)pVect_trg.size();
    int nMultAsso = (int)pVect_ass.size();
    
    for(int ntrg=0; ntrg<nMultTrg; ++ntrg)
    {
        TVector3 pvector_trg = (pVect_trg)[ntrg];
        double eta_trg = pvector_trg.Eta();
        double phi_trg = pvector_trg.Phi();
        
        for(int nass=0; nass<nMultAsso; nass++)
        {
            TVector3 pvector_ass = (pVect_ass)[nass];
            double eta_ass = pvector_ass.Eta();
            double phi_ass = pvector_ass.Phi();
            
            double deltaEta = eta_ass - eta_trg;
            double deltaPhi = phi_ass - phi_trg;
            if(deltaPhi > pi_) deltaPhi = deltaPhi - 2*pi_;
            if(deltaPhi < -pi_) deltaPhi = deltaPhi + 2*pi_;
            if(deltaPhi > -pi_ && deltaPhi < -pi_/2.0) deltaPhi = deltaPhi + 2*pi_;
            
            if(deltaEta == 0 && deltaPhi == 0) continue;
            
            hSignal->Fill(deltaEta,deltaPhi,1.0/nMultTrg);
        
        } //Loop over associated particles
    } //Loop over trigger particles
    
    pVectVect_trg.push_back(pVect_trg);
    pVectVect_ass.push_back(pVect_ass);
    zvtxVect.push_back(vsorted[0].z());

}

void
AwayAnalyzer::initHistos(const edm::Service<TFileService> & fs)
{

  events_ = fs->make<TH1F>("events","",1,0,1);
  vertices_ = fs->make<TH1F>("vertices","",1,0,1);

  evtPerf_["Ntrk"] = fs->make<TH1F>("evtNtrk","Tracks per event",100,0,400);
  evtPerf_["NHPtrk"] = fs->make<TH1F>("evtHPNtrk","High purity tracks per event",100,0,400);
  evtPerf_["Nvtx"] = fs->make<TH1F>("evtNvtx","Primary Vertices per event",10,0,10);
  evtPerf_["NvtxLumi"] = fs->make<TH1F>("evtNvtxLumi","Primary Vertices by Lumi",200,0,2000);
  evtPerf_["Lumi"] = fs->make<TH1F>("evtLumi","Events by Lumi",200,0,2000);

  vtxPerf_["Ntrk"] = fs->make<TH1F>("vtxNtrk","Tracks per vertex",50,0,200);
  vtxPerf_["x"] = fs->make<TH1F>("vtxX","Vertex x position",1000,-1,1);
  vtxPerf_["y"] = fs->make<TH1F>("vtxY","Vertex y position",1000,-1,1);
  vtxPerf_["z"] = fs->make<TH1F>("vtxZ","Vertex z position",100,-30,30);
  vtxPerf_["assocVtxDz"] = fs->make<TH1F>("assocVtxDz","Z Distance from first PV; dz (cm)",200,0,50);
  vtxPerf_["assocVtxDxy"] = fs->make<TH1F>("assocVtxDxy","Rho Distance from first PV; dxy (cm)",200,0,4);

  vtxPerf2D_["Ntrk2D"] = fs->make<TH2F>("vtxNtrk2D","Tracks per vertex;vertex (sorted by Ntrk);Ntrk"
					,10,0,10,200,0,200);
  vtxPerf2D_["assocVtxDzNtrk"] = fs->make<TH2F>("assocVtxDzNtrk",
						"Z Distance from first PV vs Ntrk of assoc; dz (cm); Ntrk",
						200,0,50,50,0,200);
  vtxPerf2D_["assocVtxDxyNtrk"] = fs->make<TH2F>("assocVtxDxyNtrk",
						 "Rho Distance from first PV vs Ntrk of assoc; dxy (cm); Ntrk",
						 200,0,4,50,0,200);

  trkPerf_["Nhit"] = fs->make<TH1F>("trkNhit", "Tracks by Number of Valid Hits;N hits",    35,  0,35);
  trkPerf_["pt"] = fs->make<TH1F>("trkPt", "Track p_{T} Distribution;p_{T} [GeV/c]",100,0,50);
  trkPerf_["ptHigh"] = fs->make<TH1F>("trkPtHigh", "Track p_{T} Distribution;p_{T} [GeV/c]",100,0,200);
  trkPerf_["eta"] = fs->make<TH1F>("trkEta", "Track iConfigeudorapidity Distribution;#eta",50,-2.5,2.5);
  trkPerf_["phi"] = fs->make<TH1F>("trkPhi", "Track Azimuthal Distribution;#phi",100,-3.15,3.15);
  trkPerf_["chi2"] = fs->make<TH1F>("trkChi2", "Track Normalized #chi^{2};#chi^{2}/n.d.o.f",60,0,6);
  trkPerf_["pterr"] = fs->make<TH1F>("trkPterr", "Track p_{T} error;#delta p_{T} / p_{T}",50,0,0.2);
  trkPerf_["dxyErr"] = fs->make<TH1F>("trkDxyErr", "Transverse DCA Significance;dxy / #sigma_{dxy}",100,-8,8);
  trkPerf_["dzErr"] = fs->make<TH1F>("trkDzErr", "Longitudinal DCA Significance;dz / #sigma_{dz}",100,-8,8);

  char histoName1[300];
  char histoTitle1[1000];
  char histoName2[200];
  char histoTitle2[200];
    
  for(int kVz=0; kVz<int(vzBins_.size()-1); kVz++) {
      
      sprintf(histoName1, "hdNdEta_VzBin_%d", kVz);
      sprintf(histoTitle1, "dNdEta distribution for %5.2f < V_{z} < %5.2f ", vzBins_[kVz], vzBins_[kVz+1]);
      hdNdEtaVzBin_[histoName1] = fs->make<TH1F>(histoName1, histoTitle1, 51, etaMin_, etaMax_);
      
      sprintf(histoName2, "nEventsVzBin_%d", kVz);
      sprintf(histoTitle2, "No of events for %5.2f < V_{z} < %5.2f ", vzBins_[kVz], vzBins_[kVz+1]);
      hEventVzBin_[histoName2] = fs->make<TH1F>(histoName2, histoTitle2, 1, 0, 1);
    }
    
   hSignal = fs->make<TH2D>("hSignal", ";#Delta#eta;#Delta#phi", 33,-5.0-0.15,5.0+0.15,31,-pi_/2+pi_/32,3*pi_/2-pi_/32);
   hBackground = fs->make<TH2D>("hBackground", ";#Delta#eta;#Delta#phi", 33,-5.0-0.15,5.0+0.15,31,-pi_/2+pi_/32,3*pi_/2-pi_/32);
    
   trkPerf_["ptAsso"] = fs->make<TH1F>("trkPtAsso", "Associated Track p_{T} Distribution;p_{T} [GeV/c]",100,0,10);
   trkPerf_["ptTrg"] = fs->make<TH1F>("trkPtTrg", "Trigger Track p_{T} Distribution;p_{T} [GeV/c]",100,0,10);
   trkPerf_["etaAsso"] = fs->make<TH1F>("trkEtaAsso", "Associated Track pseudorapidity Distribution;#eta",51,-2.5,2.5);
   trkPerf_["etaTrg"] = fs->make<TH1F>("trkEtaTrg", "Trigger Track pseudorapidity Distribution;#eta",51,-2.5,2.5);
   trkPerf_["phiAsso"] = fs->make<TH1F>("trkPhiAsso", "Associated Track Azimuthal Distribution;#phi",100,-3.15,3.15);
   trkPerf_["phiTrg"] = fs->make<TH1F>("trkPhiTrg", "Trigger Track Azimuthal Distribution;#phi",100,-3.15,3.15);
}

bool
AwayAnalyzer::vtxSort( const reco::Vertex &  a, const reco::Vertex & b )
{
  if( a.tracksSize() != b.tracksSize() )
    return  a.tracksSize() > b.tracksSize() ? true : false ;
  else
    return  a.chi2() < b.chi2() ? true : false ;
}

bool
AwayAnalyzer::TrackQualityCuts(const reco::Track & track, const reco::Vertex & vertexCollectionSelected)
{
    
    math::XYZPoint vtxPoint(0.0,0.0,0.0);
    double vzErr =0.0, vxErr=0.0, vyErr=0.0;
    vtxPoint=vertexCollectionSelected.position();
    vzErr=vertexCollectionSelected.zError();
    vxErr=vertexCollectionSelected.xError();
    vyErr=vertexCollectionSelected.yError();
    
    double dxy=0.0, dz=0.0, dxysigma=0.0, dzsigma=0.0;
    dxy = track.dxy(vtxPoint);
    dz = track.dz(vtxPoint);
    dxysigma = sqrt(track.d0Error()*track.d0Error()+vxErr*vyErr);
    dzsigma = sqrt(track.dzError()*track.dzError()+vzErr*vzErr);
    
    if(track.quality(reco::TrackBase::qualityByName(qualityString_)) != 1)
        return false;
    if(fabs(dxy/dxysigma) > cutDxyErrMax_) return false;
    if(fabs(dz/dzsigma) > cutDzErrMax_) return false;
    if(track.ptError() / track.pt() > cutPtErrMax_) return false;
    
    return true;
    
}



// ------------ method called once each job just before starting event loop  ------------
void
AwayAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
AwayAnalyzer::endJob()
{
    ////////////////// Calculating background for pi0-hadron correlations //////////
    int nevttotal_trg = (int)pVectVect_trg.size();
    int nevttotal_ass = (int)pVectVect_ass.size();
    
    unsigned int nBackgroundFill = 0;
    unsigned int countAcceptedTriggerEvents = 0;
    bool *acceptedTriggerEvents = new bool [nevttotal_trg];
        
    // Look at Z-Vertex differences
    for(int nevt_trg=0; nevt_trg<nevttotal_trg; nevt_trg++) {
        acceptedTriggerEvents[nevt_trg] = false;
        for(int nevt_ass=0; nevt_ass<nevttotal_ass; nevt_ass++) {
                
            if(nevt_trg == nevt_ass)
                continue;  // don't use the same event
            if(fabs((zvtxVect)[nevt_trg] - (zvtxVect)[nevt_ass])<=2.0) {
                acceptedTriggerEvents[nevt_trg] = true;
                countAcceptedTriggerEvents++;
                break; // found at least one partner event with a close enough Z vertex
            }
        } // loop over associated events
    } // loop over trigger events
 
    for(int nevt_trg=0; nevt_trg<nevttotal_trg; nevt_trg++) {
        
        if(!acceptedTriggerEvents[nevt_trg])
            continue;  // skip this trigger event which has no partner event close enough in Z
        nBackgroundFill++;
        
        int countGoodAssociated = 0;
        
        for(int nevt_ass=0; nevt_ass<nevttotal_ass; nevt_ass++) {
            if(nevt_trg == nevt_ass) { // check if the random trigger event matches this assocated event
                continue;    // go to the next associated track
            }
            
            if(fabs((zvtxVect)[nevt_trg] - (zvtxVect)[nevt_ass])>2.0){  // check if the Z vertex of the trigger and associated are separated by more than 2 cm
                continue;    // go to the next associated event
            }
            countGoodAssociated++;
        }
        
        if(countGoodAssociated < 1) {
            cout << "\n For nevt_trg " << nevt_trg << " the number of good associated events = " << countGoodAssociated << endl;
            continue;
        }
        
        int takeRandomInterval = 1;
        if(countGoodAssociated > 10)
            takeRandomInterval = countGoodAssociated/10 + 1 ;
        
        int takeAssociated = 0;
        for(int nevt_ass=0; nevt_ass<nevttotal_ass; nevt_ass += takeRandomInterval) {
            
            if(nevt_trg == nevt_ass) { // check if the random trigger event matches this assocated event
                continue;    // go to the next associated track
            }
            
            if(fabs((zvtxVect)[nevt_trg] - (zvtxVect)[nevt_ass])>2.0){  // check if the Z vertex of the trigger and associated are separated by more than 2 cm
                continue;    // go to the next associated event
            }
            
            takeAssociated++;
            if(takeAssociated > 10)
                break;
            
            vector<TVector3> pVectTmp_trg = (pVectVect_trg)[nevt_trg];
            vector<TVector3> pVectTmp_ass = (pVectVect_ass)[nevt_ass];
            int nMult_trg1 = pVectTmp_trg.size();
            int nMult_ass1 = pVectTmp_ass.size();
            
            for(int ntrg=0; ntrg<nMult_trg1; ++ntrg)
            {
                TVector3 pvectorTmp_trg = pVectTmp_trg[ntrg];
                double eta_trg = pvectorTmp_trg.Eta();
                double phi_trg = pvectorTmp_trg.Phi();
                
                for(int nass=0; nass<nMult_ass1; ++nass)
                {
                    
                    TVector3 pvectorTmp_ass = pVectTmp_ass[nass];
                    double eta_ass = pvectorTmp_ass.Eta();
                    double phi_ass = pvectorTmp_ass.Phi();
                    
                    double deltaEta = eta_ass - eta_trg;
                    double deltaPhi = phi_ass - phi_trg;
                    if(deltaPhi > pi_) deltaPhi = deltaPhi - 2*pi_;
                    if(deltaPhi < -pi_) deltaPhi = deltaPhi + 2*pi_;
                    if(deltaPhi > -pi_ && deltaPhi < -pi_/2.0) deltaPhi = deltaPhi + 2*pi_;
                    
                    hBackground->Fill(deltaEta,deltaPhi,1.0/nMult_trg1);
                }
            }
        }
    }
    
    delete [] acceptedTriggerEvents;
}

//define this as a plug-in
DEFINE_FWK_MODULE(AwayAnalyzer);
