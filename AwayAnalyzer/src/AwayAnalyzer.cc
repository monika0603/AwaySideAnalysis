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
      std::map<std::string,TH2F*> trkPerf2D_;
      std::map<std::string,TH3F*> trkPerf3D_;
      std::map<std::string,TH1F*> vtxPerf_;
      std::map<std::string,TH2F*> vtxPerf2D_;
      std::map<std::string,TH1F*> hdNdEtaVzBin_;
      std::map<std::string,TH1F*> hEventVzBin_;

      TH1F* events_;
      TH1F* vertices_;
      TH1F* tracks_;
      TH1F* trackseta_;

      int nevt_;
      int ntrack_;
      int nvertex_;
      int tHighPurityTracks_;
      int nVzBins;

      edm::InputTag vertexSrc_;
      edm::InputTag trackSrc_;
      edm::InputTag jetSrc_;
      double etaMin_;
      double etaMax_;
      double ptMin_;
      double vertexZMax_;

      std::string qualityString_;

      std::vector<double> ptBins_;
      std::vector<double> vzBins_;
      std::vector<double> etaBins_;
      std::vector<double> NptBins_;
    
      double jetEtaMax_;
      double jetEtMin_;
      double jetEtMax_;
      double cutDzErrMax_;
      double cutDxyErrMax_;
      double cutPtErrMax_;
      double cutMultMin_;
      double cutMultMax_;
      double cutMinTrack_;

};

AwayAnalyzer::AwayAnalyzer(const edm::ParameterSet& iConfig):
nevt_(0),
ntrack_(0),
nvertex_(0),
vertexSrc_(iConfig.getParameter<edm::InputTag>("vertexSrc")),
trackSrc_(iConfig.getParameter<edm::InputTag>("trackSrc")),
jetSrc_(iConfig.getParameter<edm::InputTag>("jetSrc")),
etaMin_(iConfig.getParameter<double>("etaMin")),
etaMax_(iConfig.getParameter<double>("etaMax")),
ptMin_(iConfig.getParameter<double>("ptMin")),
vertexZMax_(iConfig.getParameter<double>("vertexZMax")),
qualityString_(iConfig.getParameter<std::string>("qualityString")),
ptBins_(iConfig.getParameter<std::vector<double> >("ptBins")),
etaBins_(iConfig.getParameter<std::vector<double> >("etaBins"))
{
    edm::Service<TFileService> fs;
    initHistos(fs);

    ptBins_ = iConfig.getParameter<std::vector<double> >("ptBins");
    vzBins_ = iConfig.getParameter<std::vector<double> >("vzBins");
    NptBins_ = iConfig.getParameter<std::vector<double> >("NptBins");
    cutMultMin_ = iConfig.getParameter<double>("cutMultMin");
    cutMultMax_ = iConfig.getParameter<double>("cutMultMax");
    cutMinTrack_ = iConfig.getParameter<double>("cutMinTrack");
    cutDzErrMax_ = iConfig.getUntrackedParameter<double>("cutDzErrMax", 3.0);
    cutDxyErrMax_ = iConfig.getUntrackedParameter<double>("cutDxyErrMax", 3.0);
    cutPtErrMax_ = iConfig.getUntrackedParameter<double>("cutPtErrMax", 0.1);
    
    nVzBins = vzBins_.size()-1;
    std::cout<<"Finished initializing the parameters"<<std::endl;

}


AwayAnalyzer::~AwayAnalyzer()
{
}

void
AwayAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   tHighPurityTracks_ = 0;
    
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
        
      /*  char histoName1[200];
        char histoName2[200];
        for(int kVz=0; kVz<nVzBins; kVz++) {
            if(vtxPoint.z() > vzBins_[kVz] && vtxPoint.z() <= vzBins_[kVz+1])
            {
                
                sprintf(histoName1, "hdNdEta_VzBin_%d", kVz);
                hdNdEtaVzBin_[histoName1]->Fill(track.eta());
                
                sprintf(histoName2, "nEventsVzBin_%d", kVz);
                hEventVzBin_[histoName2]->Fill(0.5);
            }
        }*/
        
        if( track.eta() <= etaMax_ && track.eta() >= etaMin_ && track.pt() > ptMin_)
        {
            trackseta_->Fill(0.5);
            trkPerf_["Nhit"]->Fill(track.numberOfValidHits());
            trkPerf_["pt"]->Fill(track.pt());
            trkPerf_["eta"]->Fill( track.eta() );
            trkPerf2D_["etaphi"]->Fill( track.eta(), track.phi() );
            trkPerf_["ptHigh"]->Fill(track.pt());
            trkPerf_["phi"]->Fill(track.phi());
            trkPerf_["dxyErr"]->Fill(dxy/dxysigma);
            trkPerf_["dzErr"]->Fill(dz/dzsigma);
            trkPerf_["chi2"]->Fill(track.normalizedChi2());
            trkPerf_["pterr"]->Fill(track.ptError() / track.pt() );
            trkPerf2D_["etavz"]->Fill( vtxPoint.z(), track.eta() );
            trkPerf3D_["Nhit3D"]->Fill(track.eta(), track.pt(), track.numberOfValidHits());
            trkPerf3D_["phi3D"]->Fill(track.eta(), track.pt(), track.phi());
            trkPerf3D_["dxyErr3D"]->Fill(track.eta(), track.pt(), dxy/dxysigma);
            trkPerf3D_["dzErr3D"]->Fill(track.eta(), track.pt(), dz/dzsigma);
            trkPerf3D_["chi23D"]->Fill(track.eta(), track.pt(), track.normalizedChi2());
            trkPerf3D_["pterr3D"]->Fill(track.eta(), track.pt(), track.ptError() / track.pt() );
            ntrack_++;
        }
        
    }
    
   // tAllTracks_->Fill(nHITracks);
    
    if( !(tHighPurityTracks_ >= cutMultMin_ && tHighPurityTracks_ < cutMultMax_)) return;
   // tHPTracks_->Fill(tHighPurityTracks_);

}

void
AwayAnalyzer::initHistos(const edm::Service<TFileService> & fs)
{

  events_ = fs->make<TH1F>("events","",1,0,1);
  tracks_ = fs->make<TH1F>("tracks","",1,0,1);
  trackseta_ = fs->make<TH1F>("trackseta","",1,0,1);
  vertices_ = fs->make<TH1F>("vertices","",1,0,1);


  std::vector<double> dumBins;
  dumBins.clear();
  for( double i = 0.; i<36.; i += 1.) dumBins.push_back(i);
  trkPerf3D_["Nhit3D"] = fs->make<TH3F>("trkNhit3D", "Tracks by Number of Valid Hits;N hits",
					etaBins_.size()-1, &etaBins_[0],
					ptBins_.size()-1, &ptBins_[0],
					dumBins.size()-1, &dumBins[0]);
  dumBins.clear();
  for( double i = -3.15; i<3.151; i += 6.30/100.) dumBins.push_back(i);
  trkPerf3D_["phi3D"] = fs->make<TH3F>("trkPhi3D", "Track Azimuthal Distribution;#phi",
				       etaBins_.size()-1, &etaBins_[0],
				       ptBins_.size()-1, &ptBins_[0],
				       dumBins.size()-1, &dumBins[0]);
  dumBins.clear();
  for( double i = 0.; i<6.01; i += 6./60.) dumBins.push_back(i);
  trkPerf3D_["chi23D"] = fs->make<TH3F>("trkChi23D", "Track Normalized #chi^{2};#chi^{2}/n.d.o.f",
					etaBins_.size()-1, &etaBins_[0],
					ptBins_.size()-1, &ptBins_[0],
					dumBins.size()-1, &dumBins[0]);
  dumBins.clear();
  for( double i = 0.0; i<0.201; i += 0.2/50.) dumBins.push_back(i);
  trkPerf3D_["pterr3D"] = fs->make<TH3F>("trkPterr3D", "Track p_{T} error;#delta p_{T} / p_{T}",
					 etaBins_.size()-1, &etaBins_[0],
					 ptBins_.size()-1, &ptBins_[0],
					 dumBins.size()-1, &dumBins[0]);
  dumBins.clear();
  for( double i = -8.; i<8.01; i += 16./100.) dumBins.push_back(i);
  trkPerf3D_["dxyErr3D"] = fs->make<TH3F>("trkDxyErr3D", "Transverse DCA Significance;dxy / #sigma_{dxy}",
					  etaBins_.size()-1, &etaBins_[0],
					  ptBins_.size()-1, &ptBins_[0],
					  dumBins.size()-1, &dumBins[0]);
  trkPerf3D_["dzErr3D"] = fs->make<TH3F>("trkDzErr3D", "Longitudinal DCA Significance;dz / #sigma_{dz}",
					 etaBins_.size()-1, &etaBins_[0],
					 ptBins_.size()-1, &ptBins_[0],
					 dumBins.size()-1, &dumBins[0]);
  evtPerf_["Ntrk"] = fs->make<TH1F>("evtNtrk","Tracks per event",100,0,400);
  evtPerf_["Nvtx"] = fs->make<TH1F>("evtNvtx","Primary Vertices per event",10,0,10);
  evtPerf_["ncoll"] = fs->make<TH1F>("ncoll","Event Ncoll from Generator",50,0,50);
  evtPerf_["b"] = fs->make<TH1F>("b","Impact Parameter from Generator",100,0,20);

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

  trkPerf2D_["etaphi"] = fs->make<TH2F>("trkEtaPhi","Track Eta-Phi Map;#eta;#phi",50,-2.5,2.5,100,-3.15,3.15);
  trkPerf2D_["etavz"] = fs->make<TH2F>("trkEtaVz","Track Eta vs Vertex z;Vertex z (cm);#eta",
                                       100,-30,30,100,-3.0,3.0);
 /* for(int kVz=0; kVz<nVzBins; kVz++) {
      char histoName1[200];
      char histoTitle1[200];
      char histoName2[200];
      char histoTitle2[200];
      sprintf(histoName1, "hdNdEta_VzBin_%d", kVz);
      sprintf(histoTitle1, "dNdEta distribution for %5.2f < V_{z} < %5.2f ", vzBins_[kVz], vzBins_[kVz+1]);
      hdNdEtaVzBin_[histoName1] = fs->make<TH1F>(histoName1, histoTitle1, 100, etaMin_, etaMax_);
      
      sprintf(histoName2, "nEventsVzBin_%d", kVz);
      sprintf(histoTitle2, "No of events for %5.2f < V_{z} < %5.2f ", vzBins_[kVz], vzBins_[kVz+1]);
      hEventVzBin_[histoName2] = fs->make<TH1F>(histoName2, histoTitle2, 1, 0, 1);
    }*/
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
}

//define this as a plug-in
DEFINE_FWK_MODULE(AwayAnalyzer);
