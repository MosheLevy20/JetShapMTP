// ################################################################
// Author:  Joel Mazer for the STAR Collaboration
// Affiliation: Rutgers University
//
// ################################################################

#include "StJetShapeMTPMaker.h"
#include "StMemStat.h"

// ROOT includes
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include <THnSparse.h>
#include "TParameter.h"
#include <TProfile.h>
#include "TRandom.h"
#include "TRandom3.h"
#include "TVector3.h"

// STAR includes
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StMaker.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StPicoEvent/StPicoBTowHit.h"
#include "StRoot/StPicoEvent/StPicoBTofHit.h"
#include "StRoot/StPicoEvent/StPicoEmcTrigger.h"
#include "StRoot/StPicoEvent/StPicoMtdTrigger.h"
#include "StRoot/StPicoEvent/StPicoBEmcPidTraits.h"
#include "StRoot/StPicoEvent/StPicoBTofPidTraits.h"
#include "StRoot/StPicoEvent/StPicoMtdPidTraits.h"

// jet-framework includes
#include "StJetFrameworkPicoBase.h"
#include "StRhoParameter.h"
#include "StRho.h"
#include "StJetMakerTask.h"
#include "StFemtoTrack.h"
#include "StEmcPosition2.h"
#include "StCentMaker.h"

// old file kept
#include "StPicoConstants.h"

// centrality includes
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"

//Event Mixing
#include "StEventPoolManager.h"


ClassImp(StJetShapeMTPMaker)


/*To Do

*/


//________________________________________________________________________
StJetShapeMTPMaker::StJetShapeMTPMaker(const char* name, StPicoDstMaker *picoMaker, const char* outName = "", const char* jetMakerName = "", const char* rhoMakerName = "") : StJetFrameworkPicoBase(name) //StMaker(name),
{
  //Important switches
  doUsePrimTracks = kFALSE;
  doppAnalysis = kFALSE;
  fDoEffCorr = kFALSE;

  //Impportant cutoff values 
  fMinPtJet = 15.0;
  fMaxEventTrackPt = 30.0;
  fTrackPtMinCut = 2.0; fTrackPtMaxCut = 30.0; //Track pt cutoffs
  fTrackPhiMinCut = 0.0; fTrackPhiMaxCut = 2.0*TMath::Pi();
  fTrackEtaMinCut = -1.0; fTrackEtaMaxCut = 1.0;

  //All other variables that need to be defined.
  fRequireCentSelection = kFALSE;
  doRejectBadRuns = kFALSE;

  fLeadingJet = 0x0; fSubLeadingJet = 0x0;
  fJets = 0x0 ;
  mPicoDstMaker = 0x0;
  mPicoDst = 0x0;
  mPicoEvent = 0x0;
  fRho = 0x0;
  mCentMaker = 0x0;
  mBaseMaker = 0x0;

  fRunNumber = 0;
  refCorr2 = 0.0;
  JetMaker = 0;
  RhoMaker = 0;
  fRhoVal = 0;
  
  fRunFlag = 0;       // see StJetShapeMTPMaker::fRunFlagEnum
  fCentralityDef = 4; // see StJetFrameworkPicoBase::fCentralityDefEnum
  fCentralitySelectionCut = -99;


  fEventZVtxMinCut = -40.0; fEventZVtxMaxCut = 40.0;
  fCentralityScaled = 0.;
  ref16 = -99; ref9 = -99;
  Bfield = 0.0;
  //mVertex = 0x0;
  zVtx = 0.0;
  fEmcTriggerEventType = 0; fMBEventType = 2;
  
  mOutName = outName;
  fAnalysisMakerName = name;
  fJetMakerName = jetMakerName;
  fRhoMakerName = rhoMakerName;
  for(int i=0; i<8; i++) { fEmcTriggerArr[i] = 0; }

  //Event Mixing
  fPoolMgr = 0x0;
  fDoEventMixing = 0; 
  fMixingTracks = 50000; 
  fNMIXtracks = 5000; 
  fNMIXevents = 1;//this
  fCentBinSize = 5; 
  fCentBinSizeJS = 10; 
  fDoUseMultBins = kFALSE;
  doUseEPBins = kFALSE;
  fnEPBins = 6;

}

//  Function: This method is called every event.
//_____________________________________________________________________________
Int_t StJetShapeMTPMaker::Make() {
  // zero out these global variables
  fCentralityScaled = 0.0, ref9 = 0, ref16 = 0;

  // get PicoDstMaker 
  mPicoDstMaker = static_cast<StPicoDstMaker*>(GetMaker("picoDst"));
  if(!mPicoDstMaker) {
    LOG_WARN << " No PicoDstMaker! Skip! " << endm;
    return kStWarn;
  }

  // construct PicoDst object from maker
  mPicoDst = static_cast<StPicoDst*>(mPicoDstMaker->picoDst());
  if(!mPicoDst) {
    LOG_WARN << " No PicoDst! Skip! " << endm;
    return kStWarn;
  }

  // create pointer to PicoEvent 
  mPicoEvent = static_cast<StPicoEvent*>(mPicoDst->event());
  if(!mPicoEvent) {
    LOG_WARN << " No PicoEvent! Skip! " << endm;
    return kStWarn;
  }

  // get base class pointer
  mBaseMaker = static_cast<StJetFrameworkPicoBase*>(GetMaker("baseClassMaker"));
  if(!mBaseMaker) {
    LOG_WARN << " No baseMaker! Skip! " << endm;
    return kStWarn;
  }

  // get bad run, dead & bad tower lists
  badRuns = mBaseMaker->GetBadRuns();
  deadTowers = mBaseMaker->GetDeadTowers();
  badTowers = mBaseMaker->GetBadTowers();

  // get run number, check bad runs list if desired (kFALSE if bad)
  fRunNumber = mPicoEvent->runId();
  if(doRejectBadRuns) {
    if( !mBaseMaker->IsRunOK(fRunNumber) ) return kStOK;
  }

  //for(int i = 1; i<4801; i++) {  if(!mBaseMaker->IsTowerOK(i))  cout<<"tower: "<<i<<" is not good!!"<<endl;  }
  // can check for bad towers like this (if needed):
  // bool isTowOK = mBaseMaker->IsTowerOK(towerID);
  // ===========================================================================================

  // cut event on max track pt > 30.0 GeV
  if(GetMaxTrackPt() > fMaxEventTrackPt) return kStOK;

  // cut event on max tower Et > 30.0 GeV
  //if(GetMaxTowerEt() > fMaxEventTowerEt) return kStOK;

  // get event B (magnetic) field
  Bfield = mPicoEvent->bField(); 

  // get vertex 3-vector and z-vertex component
  mVertex = mPicoEvent->primaryVertex();
  zVtx = mVertex.z();
  
  // Z-vertex cut: the Aj analysis cut on (-40, 40) for reference
  if((zVtx < fEventZVtxMinCut) || (zVtx > fEventZVtxMaxCut)) return kStOk;


  // ============================ CENTRALITY ============================== //
  // get CentMaker pointer
  mCentMaker = static_cast<StCentMaker*>(GetMaker("CentMaker"));
  if(!mCentMaker) {
    LOG_WARN << " No CenttMaker! Skip! " << endm;
    return kStWarn;
  }

  // centrality variables
  int grefMult = mCentMaker->GetgrefMult(); // see StPicoEvent
  int refMult =  mCentMaker->GetrefMult();  // see StPicoEvent
  ref9 = mCentMaker->GetRef9();   // binning from central -> peripheral
  ref16 = mCentMaker->GetRef16(); // binning from central -> peripheral
  int cent16 = mCentMaker->GetCent16(); // centrality bin from StRefMultCorr (increasing bin corresponds to decreasing cent %) - Don't use except for cut below
  int centbin = mCentMaker->GetRef16();
  double refCorr2 = mCentMaker->GetRefCorr2();
  fCentralityScaled = mCentMaker->GetCentScaled();
  //double refCorr = mCentMaker->GetCorrectedMultiplicity(refMult, zVtx, zdcCoincidenceRate, 0); // example usage
  // for pp analyses:    centbin = 0, cent9 = 0, cent16 = 0, refCorr2 = 0.0, ref9 = 0, ref16 = 0;

  // cut on unset centrality, > 80%
  if(cent16 == -1) return kStOk; // this is for lowest multiplicity events 80%+ centrality, cut on them 

  double kEventActivity = (doppAnalysis) ? (double)grefMult : refCorr2;
  

  // cut on centrality for analysis before doing anything
  if(fRequireCentSelection) { if(!SelectAnalysisCentralityBin(centbin, fCentralitySelectionCut)) return kStOk; }

  // ============================ end of CENTRALITY ============================== //

  // ========================= Trigger Info =============================== //
  // looking at the EMCal triggers - used for QA and deciding on HT triggers
  FillEmcTriggers();

  // get trigger IDs from PicoEvent class and loop over them
  vector<unsigned int> mytriggers = mPicoEvent->triggerIds();
  if(fDebugLevel == StJetFrameworkPicoBase::kDebugEmcTrigger) cout<<"EventTriggers: ";
  for(unsigned int i=0; i<mytriggers.size(); i++) {
    if(fDebugLevel == StJetFrameworkPicoBase::kDebugEmcTrigger) cout<<"i = "<<i<<": "<<mytriggers[i] << ", ";
  }
  if(fDebugLevel == StJetFrameworkPicoBase::kDebugEmcTrigger) cout<<endl;

  // check for MB/HT event
  bool fHaveMBevent = CheckForMB(fRunFlag, fMBEventType);
  bool fHaveMB5event = CheckForMB(fRunFlag, StJetFrameworkPicoBase::kVPDMB5);
  bool fHaveMB30event = CheckForMB(fRunFlag, StJetFrameworkPicoBase::kVPDMB30);
  bool fHaveEmcTrigger = CheckForHT(fRunFlag, fEmcTriggerEventType);
  
  bool fRunForMB = kFALSE;  // used to differentiate pp and AuAu
  if(doppAnalysis)  fRunForMB = (fHaveMBevent) ? kTRUE : kFALSE;
  if(!doppAnalysis) fRunForMB = (fHaveMB5event || fHaveMB30event) ? kTRUE : kFALSE;
  // ======================== end of Triggers ============================= //

  // =========================== JetMaker =============================== //
  // get JetMaker pointer
  JetMaker = static_cast<StJetMakerTask*>(GetMaker(fJetMakerName));
  const char *fJetMakerNameCh = fJetMakerName;
  if(!JetMaker) {
    LOG_WARN << Form(" No %s! Skip! ", fJetMakerNameCh) << endm;
    return kStWarn;
  }


  // ============================= RhoMaker ============================== //
  // get RhoMaker pointer from event: old names "StRho_JetsBG"
  RhoMaker = static_cast<StRho*>(GetMaker(fRhoMakerName));
  const char *fRhoMakerNameCh = fRhoMakerName;
  if(!RhoMaker) {
    LOG_WARN << Form(" No %s! Skip! ", fRhoMakerNameCh) << endm;
    return kStWarn;
  }

  // set rho object, alt fRho = GetRhoFromEvent(fRhoName);
  fRho = static_cast<StRhoParameter*>(RhoMaker->GetRho());
  if(!fRho) {
    LOG_WARN << Form("Couldn't get fRho object! ") << endm;
    return kStWarn;    
  } 
  
  // get rho/area value from rho object     fRho->ls("");
  fRhoVal = fRho->GetVal();


  //=================================================================================
  //COPIED FROM JOEL'S "StMyAnalysisMaker3"; CAN BE CLEANED UP
  // require event mixing
  // 1. First get an event pool corresponding in mult (cent) and
  //    zvertex to the current event. Once initialized, the pool
  //    should contain nMix (reduced) events. This routine does not
  //    pre-scan the chain. The first several events of every chain
  //    will be skipped until the needed pools are filled to the
  //    specified depth. If the pool categories are not too rare, this
  //    should not be a problem. If they are rare, you could lose
  //    statistics.

  // 2. Collect the whole pool's content of tracks into one TObjArray
  //    (bgTracks), which is effectively a single background super-event.

  // 3. The reduced and bgTracks arrays must both be passed into
  //    FillCorrelations() (diff. now). Also nMix should be passed in, so a weight
  //    of 1./nMix can be applied.

  // mix jets from triggered events with tracks from MB events
  // get the trigger bit, need to change trigger bits between different runs

  // declare pool pointer
  StEventPool *pool = 0x0;
  // require event mixing
  if(fDoEventMixing > 0) {
    // convert back to integer bins for mixed event pool - 10% bins (0, 7), 5% bins (0, 15)
    Int_t mixcentbin = TMath::Floor(fCentralityScaled / fCentBinSizeJS);
    //cout<<"fCentralityScaled: "<<fCentralityScaled<<"  fCentBinSizeJS: "<<fCentBinSizeJS<<"  mixcentbin: "<<mixcentbin<<"  zVtx: "<<zVtx<<endl;
    // get the generic event plane angle of the event: using tracks 0.2-2.0 GeV for calculation (option 4)
    // for an angle (0, pi)    
    ///const char *fEventPlaneMakerNameChTemp = fEventPlaneMakerName;
    ///StEventPlaneMaker *EPMaker = static_cast<StEventPlaneMaker*>(GetMaker(Form("%s%i", fEventPlaneMakerNameChTemp, 4)));
    ///double psi2 = (EPMaker) ? (double)EPMaker->GetTPCEP() : -999;
    // initialize event pools - different cases for each dataset
    if(fDoUseMultBins) {
      pool = fPoolMgr->GetEventPool(kEventActivity, zVtx);
      
    } else {
      pool = fPoolMgr->GetEventPool(mixcentbin, zVtx);      
    }
    // check if pool exists
    if(!pool) {
      Form("No pool found for centrality = %.1f, zVtx = %f", (double)mixcentbin, zVtx);
      return kStOK;     
    }
  }
  
  //if(fHaveEmcTrigger) {
  RunJets(pool);  
  //}     // inclusive jet case
  // ==================================================================================================
  // use only tracks from MB events
  //if(fDoEventMixing > 0 && fRunForMB && (!fHaveEmcTrigger)) { // kMB5 or kMB30 - AuAu, kMB - pp (excluding HT)
  if(fDoEventMixing > 0 && fRunForMB){// && pool->GetCurrentNEvents()<2) { // kMB5 or kMB30 - AuAu, kMB - pp (don't exclude HT)
    // update pool: create a list of reduced objects. This speeds up processing and reduces memory consumption for the event pool
    //cout<<"mb event"<<endl;
    pool->UpdatePool(CloneAndReduceTrackList());
    
  } // MB + event mixing 
  //============================================================= 

  return kStOK;
}


void StJetShapeMTPMaker::RunJets(StEventPool *pool)
{
  // need unsubbed tracks for jet shape (but use subbed jet pt) and need to do event mixing (and eta ref) sub on jet shape after
  //the fact
  Int_t njets; 
  
  fJets = static_cast<TClonesArray*>(JetMaker->GetJetsBGsub());
  //fJets = static_cast<TClonesArray*>(JetMaker->GetJets());
  
  njets = fJets->GetEntries();

  //cout<<"njets: "<<njets <<endl;

  for(int ijet = 0; ijet < njets; ijet++) {  // JET LOOP
    // get jet pointer
    
    StJet *jet = static_cast<StJet*>(fJets->At(ijet));
    if(!jet){ 
      cout<<"no jet"<<endl;
      continue;
    }
    
    //cout<<"a jet!" <<endl;
    // get some jet parameters
    double jetarea = jet->Area();
    double jetpt = jet->Pt();
    double jetEta = jet->Eta();
    double jetPhi = TVector2::Phi_0_2pi(jet->Phi());

    //trk = 0, jet = 1
    if (!withinCut(1, jetpt, jetEta)) continue;

    int ptBin = getPtBin(jetpt);
    if (ptBin < 0){
      continue;
    }
    int centBin = getCentBin(fCentralityScaled);

    
    //TObjArray fConstituents = moeGetJetConstituents(jet);
    TObjArray fConstituents = moeGetJetConstituentsSubbed(jet);
    if (fConstituents.GetEntries() < 2) continue;
    fillJetQA(jetpt, jetEta, jetPhi, jetarea, ptBin, centBin);
    fillTrackQA(0, fConstituents, ptBin, centBin);

    double LeSub = calculate_lesub(fConstituents);
    double Girth = calculate_girth(fConstituents, jetpt, jetEta, jetPhi);
    double PtD = calculate_ptd(fConstituents);
    fillJetShapes(0, LeSub, Girth, PtD, ptBin, centBin);
    fConstituents.Clear();
    
    //============================================
    //                  ETA REF
    //============================================
    //ETA REF INDEX = 1 

    double refEta = -jetEta;

    TObjArray fEtaRefTracks = GetEtaRefTracks(refEta, jetPhi);
    fillTrackQA(1, fEtaRefTracks, ptBin, centBin);
    double erLeSub = calculate_lesub(fEtaRefTracks);
    double erGirth = calculate_girth(fEtaRefTracks, jetpt, refEta, jetPhi);
    double erPtD = calculate_ptd(fEtaRefTracks);

    fillJetShapes(1, erLeSub, erGirth, erPtD, ptBin, centBin);
    
    fEtaRefTracks.Clear();
    

    //============================================
    //                  MIXED EVENTS
    //============================================
    //MIXED EVENTS INDEX = 2
    //if(pool->IsReady() || pool->NTracksInPool() > fNMIXtracks || pool->GetCurrentNEvents() >= fNMIXevents) {
      //cout<<"pool->IsReady()"<<pool->IsReady();
    if (fDoEventMixing == 0) continue;
    int nMixedEvents = pool->GetCurrentNEvents();
    if (nMixedEvents >= fNMIXevents) {
      vector<TObjArray> meTracksList = GetMixedEventTracks(pool, jetEta, jetPhi);

      for (int bgEvent = 0; bgEvent < nMixedEvents; ++bgEvent){
          TObjArray meTracks = meTracksList[bgEvent];
          fillTrackQA(2, meTracks, ptBin, centBin);
          double meLeSub = calculate_lesub(meTracks);
          double meGirth = calculate_girth(meTracks, jetpt, refEta, jetPhi);
          double mePtD = calculate_ptd(meTracks);
          fillJetShapes(2, meLeSub, meGirth, mePtD, ptBin, centBin, 1.0/nMixedEvents);
	         meTracks.Clear();
      }
      meTracksList.clear();
    }
    
  }   
}



TObjArray StJetShapeMTPMaker::moeGetJetConstituents(StJet *jet)
{ 
  TObjArray constits(0);
  vector<fastjet::PseudoJet> jetConstituents = jet->GetJetConstituents();
  int nconstituentTracks = jetConstituents.size();
  for(int itrk = 0; itrk< nconstituentTracks; itrk++){
    //get user defined index
    Int_t uid = jetConstituents[itrk].user_index();
    if(uid < 0.) continue; //makes sure we have the track array element
    //if uid = -1 its a ghost particle; otherwise its a tower; might want to use

    StPicoTrack* tempTrk = static_cast<StPicoTrack*>(mPicoDst->track(uid));
    StFemtoTrack *trk = new StFemtoTrack(tempTrk, Bfield, mVertex, doUsePrimTracks);
    if(!trk){ continue; }
    //if(!AcceptTrack(trk, Bfield, mVertex)) { continue; }
   
    double trkpt = trk->Pt();
    double trkphi = trk->Phi_0_2pi();
    double trketa = trk->Eta();

    //if (!(trkpt>fTrackPtMinCut && trkpt<30.0 && trketa>-1 && trketa<1)) continue;
    if (!withinCut(0, trkpt, trketa)) continue; 
    constits.Add(trk);
  }
  return constits;
}


TObjArray StJetShapeMTPMaker::moeGetJetConstituentsSubbed(StJet *jet)
{ 
  TObjArray constits(0);
  vector<fastjet::PseudoJet> jetConstituents = jet->GetJetConstituents();
  int nconstituentTracks = jetConstituents.size();
  for(int itrk = 0; itrk< nconstituentTracks; itrk++){
   
    fastjet::PseudoJet track = jetConstituents.at(itrk);
    StFemtoTrack *trk = new StFemtoTrack(track.px(), track.py(), track.pz());
    if(!trk){ continue; }
    //if(!AcceptTrack(trk, Bfield, mVertex)) { continue; }
   
    double trkpt = trk->Pt();
    double trkphi = trk->Phi_0_2pi();
    //cout<<"femto track phi" << trkphi<<endl;
    double trketa = trk->Eta();

    //if (!(trkpt>fTrackPtMinCut && trkpt<30.0 && trketa>-1 && trketa<1)) continue;
    if (!withinCut(0, trkpt, trketa)) continue; 
    constits.Add(trk);
  }
  return constits;
}



TObjArray StJetShapeMTPMaker::GetEtaRefTracks(Double_t jetEtaBG, Double_t jetPhiBG){
  
  Int_t ntracks = mPicoDst->numberOfTracks();
  TObjArray etaRefTracks(0);
  for(int itrack = 0; itrack < ntracks; itrack++){

    //StFemtoTrack *trk = static_cast<StFemtoTrack*>(mPicoDst->track(itrack));
    StPicoTrack* tempTrk = static_cast<StPicoTrack*>(mPicoDst->track(itrack));
    StFemtoTrack *trk = new StFemtoTrack(tempTrk, Bfield, mVertex, doUsePrimTracks);
    
    if(!trk){ continue; }
    //if(!AcceptTrack(trk, Bfield, mVertex)) { continue; }
    
    double trkphi = trk->Phi_0_2pi();
    double trketa = trk->Eta();
    double trkpt = trk->Pt();

    //if (!(trkpt>fTrackPtMinCut && trkpt<30.0 && trketa>-1 && trketa<1)) continue;
    if (!withinCut(0, trkpt, trketa)) continue; 

    double deltaEtaBG = 1.0*TMath::Abs(jetEtaBG - trketa);
    double deltaPhiBG = 1.0*TMath::Abs(jetPhiBG - trkphi);
    if(deltaPhiBG > 1.0*pi) deltaPhiBG = 2.0*pi - deltaPhiBG;
    double deltaR = 1.0*TMath::Sqrt(deltaEtaBG*deltaEtaBG + deltaPhiBG*deltaPhiBG);

    if (deltaR > jetRadius) continue;
    etaRefTracks.Add(trk);    
  }
  return etaRefTracks;
}


vector<TObjArray> StJetShapeMTPMaker::GetMixedEventTracks(StEventPool *pool, Double_t jetEtaBG, Double_t jetPhiBG){
  int nMix = pool->GetCurrentNEvents();

  vector<TObjArray> bgTracksList;
  TObjArray *bgTracks;
  TObjArray bgFemtoTrks(0);

  for(int jMix = 0; jMix < nMix; jMix++) {
    // get jMix'th event
    bgTracks = pool->GetEvent(jMix);
    const Int_t Nbgtrks = bgTracks->GetEntries();
    for(int ibg = 0; ibg < Nbgtrks; ibg++) {
      StFemtoTrack *trk = static_cast<StFemtoTrack*>(bgTracks->At(ibg));

      if(!trk) continue;
      //if(!AcceptTrack(trk, Bfield, mVertex)) { continue; }
      double trkphi = trk->Phi_0_2pi();
      double trketa = trk->Eta();
      double trkpt = trk->Pt();
      if (!withinCut(0, trkpt, trketa)) continue; 
      
      double deltaEtaBG = 1.0*TMath::Abs(jetEtaBG - trketa);
      double deltaPhiBG = 1.0*TMath::Abs(jetPhiBG - trkphi);
      if(deltaPhiBG > 1.0*pi) deltaPhiBG = 2.0*pi - deltaPhiBG;
      double deltaR = 1.0*TMath::Sqrt(deltaEtaBG*deltaEtaBG + deltaPhiBG*deltaPhiBG);
      if (deltaR > jetRadius) continue;

      bgFemtoTrks.Add(trk);

    }
    bgTracksList.push_back(bgFemtoTrks);
    //cuttoffs?   
  }
  return bgTracksList;
}


double StJetShapeMTPMaker::calculate_lesub(TObjArray arr){
  double le = 0;
  double sub = 0;
  
  for (unsigned j = 0; j < arr.GetEntries(); j++) {
    StFemtoTrack *trk = static_cast<StFemtoTrack*>(arr.At(j));
    double constituent_pt = trk->Pt();
    if (constituent_pt>le){
      sub = le;
      le = constituent_pt;
    }  
    else if (constituent_pt>sub){
      sub = constituent_pt;
    }
  }
  return le - sub;
}

double StJetShapeMTPMaker::calculate_girth(TObjArray arr, double jet_pt, double jet_eta, double jet_phi){
  double g = 0.0;
  for (unsigned j = 0; j < arr.GetEntries(); j++) {
    StFemtoTrack *currTrk = static_cast<StFemtoTrack*>(arr.At(j));
    //cout<<"trk eta:"<<currTrk->Eta()<<" trk phi:"<<currTrk->Phi_0_2pi();
    double r = sqrt(pow(jet_eta - currTrk->Eta(), 2)+ pow(TVector2::Phi_0_2pi(jet_phi) - currTrk->Phi_0_2pi(), 2));
    //cout<<" r:"<<r<<" trk pt:"<<currTrk->Pt()<<endl;
    g += (currTrk->Pt()/jet_pt)*r;
  }

  return g;
}

double StJetShapeMTPMaker::calculate_ptd(TObjArray arr){
  double numerator = 0;
  double denominator = 0;

  for (unsigned j = 0; j < arr.GetEntries(); j++) {
    StFemtoTrack *trk = static_cast<StFemtoTrack*>(arr.At(j));
    numerator += pow(trk->Pt(), 2);

    denominator += trk->Pt();
  }
  double ptD = sqrt(numerator)/denominator;

  return ptD;
}

double StJetShapeMTPMaker::calculate_energy(double px,double py,double pz){
  //mass of pi+ is 0.139570 GeV/c^2
  double m = 0.139570;
  double e = sqrt(pow(m,2)+pow(px,2)+pow(py,2)+pow(pz,2));
  return e;
}


int StJetShapeMTPMaker::getPtBin(float pt){
  if (pt>15.0 && pt<= 25.0){
    return 1;
  }
  else if (pt>25.0){
    return 2;
  }
  else{
    return -1;
  }
}

int StJetShapeMTPMaker::getCentBin(float cent){
  if (cent<=10.0){
    return 1;
  }
  else if (cent>10.0 && cent<=30.0){
    return 2;
  }
  else if (cent>30.0 && cent<=50.0){
    return 3;
  }
  else if (cent>50.0 && cent<=80.0){
    return 4;
  }
  else{
    return 5;
  }
}

bool StJetShapeMTPMaker::withinCut(int type, double pt, double eta){
  //trk = 0, jet = 1
  if (type == 0){
    return (pt>fTrackPtMinCut && pt<30.0 && eta>-1 && eta<1);
  }
  else{
    return (pt>fMinPtJet && TMath::Abs(eta)<fJetEtaCut);// && TMath::Abs(eta)> jetRadius);
  }
}

//
//________________________________________________________________________
StJetShapeMTPMaker::~StJetShapeMTPMaker()
{ /*  */
  // destructor
  for (int i = 0; i < nJetQATypes; ++i)
  {
    for (int j = 0; j < nJetQABGSMethods; ++j)
    {
      for (int k = 0; k < nPtBins; ++k)
      {
        for (int l = 0; l < nCentBins; ++l)
        {

          if (fHistJetQA[i][j][k][l]) delete fHistJetQA[i][j][k][l];        
        }
      }
    }
  }
  //track QA plots
  for (int i = 0; i < nTrkQATypes; ++i)
  {
    for (int j = 0; j < nTrkBGSMethods; ++j)
    {
      for (int k = 0; k < nPtBins; ++k)
      {
        for (int l = 0; l < nCentBins; ++l)
        {
          if (fHistTrackQA[i][j][k][l]) delete fHistTrackQA[i][j][k][l];        
        }
      }
    }
  }

  //Jet Shape plots
  for (int i = 0; i < nJetShapeTypes; ++i)
  {
    for (int j = 0; j < nJetShapeBGSMethods; ++j)
    {
      for (int k = 0; k < nPtBins; ++k)
      {
        for (int l = 0; l < nCentBins; ++l)
        {
          if (fHistJetShapes[i][j][k][l]) delete fHistJetShapes[i][j][k][l];        
        }
      }
    }
  }

  fPoolMgr->ClearPools();
  
}

//
//________________________________________________________________________
Int_t StJetShapeMTPMaker::Init() {
  StJetFrameworkPicoBase::Init();

  // declare histograms
  DeclareHistograms();

  // position object for Emc
  mEmcPosition = new StEmcPosition2();

  // Jet TClonesArray
  fJets = new TClonesArray("StJet"); // will have name correspond to the Maker which made it
  //fJets->SetName(fJetsName);

  return kStOK;
}
//
// Function:  write to output file and close
//________________________________________________________________________
Int_t StJetShapeMTPMaker::Finish() { 
  cout << "StJetShapeMTPMaker::Finish()\n";

  //  Write histos to file and close it.
  if(mOutName!="") {
    TFile *fout = new TFile(mOutName.Data(), "UPDATE");
    fout->cd();
    fout->mkdir(GetName());
    fout->cd(GetName());
    WriteHistograms();
    fout->cd();
    fout->Write();
    fout->Close();
  }

  cout<<"End of StJetShapeMTPMaker::Finish"<<endl;
  StMemStat::PrintMem("End of Finish...");

  return kStOK;
}

//
// initialization of histograms done here + global object setup
//________________________________________________________________________
void StJetShapeMTPMaker::DeclareHistograms() {
  // Jet QA histos
  //types: pt, eta, phi, area
  //BGSubMethods: None, constit, eta ref, eta ref subbed
  //pt bins: all, 15-25, 25-inf
  //cent bins: all, 0-10, 10-30, 30-50, 50-70, 70-100
  //           [type][BGSubMethod][Pt Bin][Cent bin]
  //TH1F           *fHistJetQA[4][3][3][6];

  // Track QA histos
  //types: pt, eta, phi
  //BGSubMethods: None, constit, eta ref
  //pt, cent bins: same as jet QA
  //           [type][BGSubMethod][Pt Bin][Cent bin]


  // Jet Shape histos
  //types: Lesub, Girth, PtD
  //everything else same as jet QA


  //TODO make ranges and bin nums
  
  int QAranges[nJetQATypes][3] = {{50,0,100},{20,-1,1},{20,-1,7},{20,0,1}};
  int shapeRanges[nJetShapeTypes][3] = {{6,0,30},{50,0,1},{10,0,1}};
  
  string typeName;
  string bgsName;
  string ptName;
  string centName;
  
  int nbins;
  int minX;
  int maxX;
  //Jet QA plots
  for (int i = 0; i < nJetQATypes; ++i)
  {
    typeName = QANames[i];
    for (int j = 0; j < nJetQABGSMethods; ++j)
    {
      bgsName = BGSNames[j];
      for (int k = 0; k < nPtBins; ++k)
      {
        ptName = PtBinNames[k];
        for (int l = 0; l < nCentBins; ++l)
        {
          centName = CentBinNames[l];
          nbins = QAranges[i][0];
    minX = QAranges[i][1];
    maxX = QAranges[i][2];
          fHistJetQA[i][j][k][l] = new TH1F(("fHistJet"+typeName+"_BGSubMethod:"+bgsName+"_PtBin:"+ptName+"_CentBin:"+centName).c_str(), ("fHistJet"+typeName+"_BGSubMethod:"+bgsName+"_PtBin:"+ptName+"_CentBin:"+centName).c_str(), nbins, minX, maxX);
        }
      }
    }
  }
  //track QA plots
  for (int i = 0; i < nTrkQATypes; ++i)
  {
    typeName = QANames[i];
    for (int j = 0; j < nTrkBGSMethods; ++j)
    {
      bgsName = BGSNames[j];
      for (int k = 0; k < nPtBins; ++k)
      {
        ptName = PtBinNames[k];
        for (int l = 0; l < nCentBins; ++l)
        {
    nbins = QAranges[i][0];
    minX = QAranges[i][1];
    maxX = QAranges[i][2];
          centName = CentBinNames[l];
          fHistTrackQA[i][j][k][l] = new TH1F(("fHistTrk"+typeName+"_BGSubMethod:"+bgsName+"_PtBin:"+ptName+"_CentBin:"+centName).c_str(), ("fHistTrk"+typeName+"_BGSubMethod:"+bgsName+"_PtBin:"+ptName+"_CentBin:"+centName).c_str(), nbins, minX, maxX);
        }
      }
    }
  }

  //Jet Shape plots
  for (int i = 0; i < nJetShapeTypes; ++i)
  {
    typeName = JetShapeNames[i];
    for (int j = 0; j < nJetShapeBGSMethods; ++j)
    {
      bgsName = BGSNames[j];
      for (int k = 0; k < nPtBins; ++k)
      {
        ptName = PtBinNames[k];
        for (int l = 0; l < nCentBins; ++l)
        { 
    nbins = shapeRanges[i][0];
    minX = shapeRanges[i][1];
    maxX = shapeRanges[i][2];
          centName = CentBinNames[l];
          fHistJetShapes[i][j][k][l] = new TH1F(("fHistJet"+typeName+"_BGSubMethod:"+bgsName+"_PtBin:"+ptName+"_CentBin:"+centName).c_str(), ("fHistJet"+typeName+"_BGSubMethod:"+bgsName+"_PtBin:"+ptName+"_CentBin:"+centName).c_str(), nbins, minX, maxX);
        }
      }
    }
  }

  // =================================================================================================
  //COPIED FROM JOEL'S "StMyAnalysisMaker3"; CAN BE CLEANED UP
  // set up centrality bins for mixed events
  // for pp we need mult bins for event mixing. Create binning here, to also make a histogram from it

  // Setup for Au-Au collisions: cent bin size can only be 5 or 10% bins
  // centrality bins for mixed events
  int nCentralityBinsJS = 100 + 1;
  double multJS = 1.0;
  if(fCentBinSizeJS==5){ // will be most commonly used
    nCentralityBinsJS = 20;
    multJS = 5.0;
  } else if(fCentBinSizeJS==10){
    nCentralityBinsJS = 10;
    multJS = 10.0;
  } else if(fCentBinSizeJS==20){
    nCentralityBinsJS = 5;
    multJS = 20.0;
  }

  // this is temp as the above and various other implementation attempts would not work for both cases
  // need to look into this, but a few hours already is too much.  Don't want to have to have this hard-coded
  Int_t nCentBins = 8;
  Double_t* centralityBins = GenerateFixedBinArray(nCentBins, 0., 8.);

  // for AuAu data
  // +1 to accomodate the fact that we define bins rather than array entries
  const int nMultBinsJS = 26;  // Alt-2 - 27 values, 26 ranges
  Double_t multBinsJS[nMultBinsJS + 1] = {0, 10,15,21,31,42,53,66,   80, 95, 112, 130, 149, 169, 190, 212, 235, 257, 280, 304, 329, 355, 382, 410, 439, 469, 800};
  // TEST August 23, 2019
  //  const int nMultBinsJS = 17; // 18 values, 17 ranges
  //  Double_t multBinsJS[nMultBinsJS + 1] = {0, 10, 14, 21, 29, 40, 54, 71, 91, 115, 143, 176, 214, 257, 307, 364, 430, 800};
  Double_t *multiplicityBinsJS = multBinsJS;

  // cent bins for AuAu data
  Int_t nCentBinsJS = 20;
  Double_t* centralityBinsJSnew = GenerateFixedBinArray(nCentBinsJS, 0., 20.); 

  // for pp data
  const int nMultBinsJSpp = 7;
  //Double_t multBinsJSpp[nMultBinsJSpp + 1] = {0.0, 4., 9, 15, 25, 35, 55, 100.0, 500.0};  // 8 (9)
  Double_t multBinsJSpp[nMultBinsJSpp + 1] = {0.0, 4.0, 6.0, 8.0, 10.0, 13.0, 30., 100.};   // 7 (8)
  Double_t *multiplicityBinsJSpp = multBinsJSpp;

  // z-vertex bins for mixed events
  Int_t nZvBins  = 20; // 4 cm wide, 40 for 2 cm wide
  Double_t* zvBins = GenerateFixedBinArray(nZvBins, -40., 40.); // min/max doesn't matter as data is cut zmin/zmax

  // event plane bins for mixed events
  Int_t nEPBins = 6;   // default 6 from 0-180 degrees, (0 - pi) range
  if(fnEPBins != 6) nEPBins = fnEPBins;
  Double_t* epBins = GenerateFixedBinArray(nEPBins, 0., 1.0*pi);

  // Event Mixing
  Int_t trackDepth = fMixingTracks;
  Int_t poolsize   = 1000;  // Maximum number of events, ignored in the present implementation of AliEventPoolManager
  if(fDoUseMultBins) {
    if(!doppAnalysis && !doUseEPBins) fPoolMgr = new StEventPoolManager(poolsize, trackDepth, nMultBinsJS, (Double_t*)multiplicityBinsJS, nZvBins, (Double_t*)zvBins); // not pp
    if(!doppAnalysis &&  doUseEPBins) fPoolMgr = new StEventPoolManager(poolsize, trackDepth, nMultBinsJS, (Double_t*)multiplicityBinsJS, nZvBins, (Double_t*)zvBins, nEPBins, (Double_t*)epBins); // not pp
    if( doppAnalysis) fPoolMgr = new StEventPoolManager(poolsize, trackDepth, nMultBinsJSpp, (Double_t*)multiplicityBinsJSpp, nZvBins, (Double_t*)zvBins); // is pp

  } else { // centrality binning  
    if(!doUseEPBins) fPoolMgr = new StEventPoolManager(poolsize, trackDepth, nCentBinsJS, (Double_t*)centralityBinsJSnew, nZvBins, (Double_t*)zvBins);
    if( doUseEPBins) fPoolMgr = new StEventPoolManager(poolsize, trackDepth, nCentBinsJS, (Double_t*)centralityBinsJSnew, nZvBins, (Double_t*)zvBins, nEPBins, (Double_t*)epBins);
  }
  

  // ================================================================================================
  SetSumw2();
}
//
void StJetShapeMTPMaker::fillTrackQA(int bgMethod, TObjArray arr, int ptBin, int centBin){
  for (unsigned j = 0; j < arr.GetEntries(); j++) {
    StFemtoTrack *trk = static_cast<StFemtoTrack*>(arr.At(j));
    double trkpt = trk->Pt();
    double trketa = trk->Eta();
    double trkphi = trk->Phi_0_2pi();

    fHistTrackQA[0][bgMethod][0][0]->Fill(trkpt);
    fHistTrackQA[1][bgMethod][0][0]->Fill(trketa);
    fHistTrackQA[2][bgMethod][0][0]->Fill(trkphi);

    fHistTrackQA[0][bgMethod][ptBin][centBin]->Fill(trkpt);
    fHistTrackQA[1][bgMethod][ptBin][centBin]->Fill(trketa);
    fHistTrackQA[2][bgMethod][ptBin][centBin]->Fill(trkphi);
  }
}

void StJetShapeMTPMaker::fillJetQA(double jetpt, double jeteta, double jetphi, double jetarea, int ptBin, int centBin){
  int bgMethod = 0;
  fHistJetQA[0][bgMethod][0][0]->Fill(jetpt);
  fHistJetQA[1][bgMethod][0][0]->Fill(jeteta);
  fHistJetQA[2][bgMethod][0][0]->Fill(jetphi);
  fHistJetQA[3][bgMethod][0][0]->Fill(jetarea);

  fHistJetQA[0][bgMethod][ptBin][centBin]->Fill(jetpt);
  fHistJetQA[1][bgMethod][ptBin][centBin]->Fill(jeteta);
  fHistJetQA[2][bgMethod][ptBin][centBin]->Fill(jetphi);
  fHistJetQA[3][bgMethod][ptBin][centBin]->Fill(jetarea);
}

void StJetShapeMTPMaker::fillJetShapes(int bgMethod, double LeSub, double Girth, double PtD, int ptBin, int centBin, double weight){
  
  fHistJetShapes[0][bgMethod][0][0]->Fill(LeSub, weight);
  fHistJetShapes[1][bgMethod][0][0]->Fill(Girth, weight);
  fHistJetShapes[2][bgMethod][0][0]->Fill(PtD, weight);

  fHistJetShapes[0][bgMethod][ptBin][centBin]->Fill(LeSub, weight);
  fHistJetShapes[1][bgMethod][ptBin][centBin]->Fill(Girth, weight);
  fHistJetShapes[2][bgMethod][ptBin][centBin]->Fill(PtD, weight);
}
// write histograms
//_____________________________________________________________________________
void StJetShapeMTPMaker::WriteHistograms() {
  // writing of histograms done here
  //TODO mk more dirs for organization
  for (int i = 0; i < nJetQATypes; ++i){
    for (int j = 0; j < nJetQABGSMethods; ++j){
      for (int k = 0; k < nPtBins; ++k){
        for (int l = 0; l < nCentBins; ++l){
    if (fHistJetQA[i][j][k][l]->GetEntries() > 0){
          fHistJetQA[i][j][k][l]->Write();  
    }
        }
      }
    }
  }
  //track QA plots
  for (int i = 0; i < nTrkQATypes; ++i){
    for (int j = 0; j < nTrkBGSMethods; ++j){
      for (int k = 0; k < nPtBins; ++k){
        for (int l = 0; l < nCentBins; ++l){
    if (fHistTrackQA[i][j][k][l]->GetEntries() > 0){
          fHistTrackQA[i][j][k][l]->Write();  
    }
        }
      }
    }
  }
  //Jet Shape plots
  for (int i = 0; i < nJetShapeTypes; ++i){
    for (int j = 0; j < nJetShapeBGSMethods; ++j){
      for (int k = 0; k < nPtBins; ++k){
        for (int l = 0; l < nCentBins; ++l){
    //if (j == 3){
      //cout<<"numEntriesRaw: "<<fHistJetShapes[i][j][k][l]->GetEntries()<<"numEntriesBG: "<<fHistJetShapes[i][2][k][l]->GetEntries()<<endl;
      //fHistJetShapes[i][j][k][l]->Add(fHistJetShapes[i][2][k][l],-1);
      //cout<<"j=3"<<"type: "<<i<<"numEntries: "<<fHistJetShapes[i][j][k][l]->GetEntries()<<endl;   }
    if (fHistJetShapes[i][j][k][l]->GetEntries() > 0){
      fHistJetShapes[i][j][k][l]->Write();  
    }
        }
      }
    }
  }
}
//
// OLD user code says: //  Called every event after Make(). 
//_____________________________________________________________________________
void StJetShapeMTPMaker::Clear(Option_t *opt) {
  fJets->Clear();
}

//
//
// __________________________________________________________________________________
void StJetShapeMTPMaker::SetSumw2() {
  for (int i = 0; i < nJetQATypes; ++i){
    for (int j = 0; j < nJetQABGSMethods; ++j){
      for (int k = 0; k < nPtBins; ++k){
        for (int l = 0; l < nCentBins; ++l){
          fHistJetQA[i][j][k][l]->Sumw2();        
        }
      }
    }
  }
  //track QA plots
  for (int i = 0; i < nTrkQATypes; ++i){
    for (int j = 0; j < nTrkBGSMethods; ++j){
      for (int k = 0; k < nPtBins; ++k){
        for (int l = 0; l < nCentBins; ++l){
          fHistTrackQA[i][j][k][l]->Sumw2();        
        }
      }
    }
  }

  //Jet Shape plots
  for (int i = 0; i < nJetShapeTypes; ++i){
    for (int j = 0; j < nJetShapeBGSMethods; ++j){
      for (int k = 0; k < nPtBins; ++k){
        for (int l = 0; l < nCentBins; ++l){
          fHistJetShapes[i][j][k][l]->Sumw2();        
        }
      }
    }
  }
}
//_________________________________________________________________________
void StJetShapeMTPMaker::FillEmcTriggers() {
  // number of Emcal Triggers
  for(int i=0; i<8; i++) { fEmcTriggerArr[i] = 0; }
  Int_t nEmcTrigger = mPicoDst->numberOfEmcTriggers();

  // set kAny true to use of 'all' triggers
  fEmcTriggerArr[StJetFrameworkPicoBase::kAny] = 1;  // always TRUE, so can select on all event (when needed/wanted) 

  // loop over valid EmcalTriggers
  for(int i = 0; i < nEmcTrigger; i++) {
    // get trigger pointer
    StPicoEmcTrigger *emcTrig = static_cast<StPicoEmcTrigger*>(mPicoDst->emcTrigger(i));
    if(!emcTrig) continue;

    // fill for valid triggers
    if(emcTrig->isHT0()) { fEmcTriggerArr[StJetFrameworkPicoBase::kIsHT0] = 1; }
    if(emcTrig->isHT1()) { fEmcTriggerArr[StJetFrameworkPicoBase::kIsHT1] = 1; }
    if(emcTrig->isHT2()) { fEmcTriggerArr[StJetFrameworkPicoBase::kIsHT2] = 1; }
    if(emcTrig->isHT3()) { fEmcTriggerArr[StJetFrameworkPicoBase::kIsHT3] = 1; }
    if(emcTrig->isJP0()) { fEmcTriggerArr[StJetFrameworkPicoBase::kIsJP0] = 1; }
    if(emcTrig->isJP1()) { fEmcTriggerArr[StJetFrameworkPicoBase::kIsJP1] = 1; }
    if(emcTrig->isJP2()) { fEmcTriggerArr[StJetFrameworkPicoBase::kIsJP2] = 1; }
  }
}



TClonesArray* StJetShapeMTPMaker::CloneAndReduceTrackList()
{
  // create array for Femto tracks
  TClonesArray *tracksClone = new TClonesArray("StFemtoTrack");

  // construct variables, get # of tracks
  int nMixTracks = mPicoDst->numberOfTracks();
  int iterTrk = 0;

  // loop over tracks
  for(int i = 0; i < nMixTracks; i++) { 
    // get track pointer
    StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(i));
    if(!trk){ continue; }

    // acceptance and kinematic quality cuts
    //if(!AcceptTrack(trk, Bfield, mVertex)) { continue; }

    // get momentum vector of track - global or primary track
    TVector3 mTrkMom;
    if(doUsePrimTracks) {
      if(!(trk->isPrimary())) continue; // get primary track vector
      mTrkMom = trk->pMom();
    } else {                            // get global track vector
      mTrkMom = trk->gMom(mVertex, Bfield);
    }

    // track variables - used with alt method below
    double pt = mTrkMom.Perp();
    //double phi = mTrkMom.Phi();
    //double eta = mTrkMom.PseudoRapidity();
    //short charge = trk->charge();

    // create StFemtoTracks out of accepted tracks - light-weight object for mixing
    //  StFemtoTrack *t = new StFemtoTrack(pt, eta, phi, charge);
    StFemtoTrack *t = new StFemtoTrack(trk, Bfield, mVertex, doUsePrimTracks);
    if(!t) continue;

    // add light-weight tracks passing cuts to TClonesArray
    ((*tracksClone)[iterTrk]) =  t;

    //delete t;
    ++iterTrk;
  } // end of looping through tracks

  return tracksClone;
}









