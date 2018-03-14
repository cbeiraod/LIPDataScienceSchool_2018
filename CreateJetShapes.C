#include <iostream>
#include "TH1F.h"
#include "TMath.h"
#include "TChain.h"
#include "TSystemDirectory.h"
#include "TList.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TStyle.h"

// This macro computes jet shapes and saves them to a root tree.


// Variables which will contain the jet info

const Int_t kMaxTracks = 500;
const Int_t kMaxTowers = 500;

Int_t ntracks  ;
Int_t ntowers  ;
Float_t jetPt  ;
Float_t jetEta ;
Float_t jetPhi ;
Float_t jetMass;
// Tracks
Float_t trackPt     [kMaxTracks];
Float_t trackEta    [kMaxTracks];
Float_t trackPhi    [kMaxTracks];
Float_t trackCharge [kMaxTracks];
// Towers
Float_t towerE      [kMaxTowers];
Float_t towerEem    [kMaxTowers];
Float_t towerEhad   [kMaxTowers];
Float_t towerEta    [kMaxTowers];
Float_t towerPhi    [kMaxTowers];

void OpenFiles(const char * inputDir);

TChain * tree = 0;

void CreateJetShapes(const char * inputDir , // Loops over all root files in
                                             // this folder
                     const char * fileOut,   // Output file
                     Int_t nentries = -1     // Nunber of jets to be processed
                     ) {

  Float_t mass, shapeLeSub, shapeRadial, shapeDispersion,ntowersLoc;
  Float_t trackChargeSum, trackAveCharge, trackWeightCharge;
  Float_t minTrackPt, maxTrackPt, meanTrackPt, sumTrackPt, weightTrackPt;
  Float_t minTrackDR, maxTrackDR, meanTrackDR;
  Float_t minTowerE, maxTowerE, meanTowerE, sumTowerE, weightTowerE;
  Float_t minTowerEEem, minTowerEEhad, maxTowerEEem, maxTowerEEhad;
  Float_t minTowerEem, maxTowerEem, meanTowerEem;
  Float_t minTowerEhad, maxTowerEhad, meanTowerEhad;
  Float_t towerTrackRatio;

  TFile fOut(fileOut,"recreate");
  TTree *treeOut = new TTree("treeShapes","example jet shapes");

  treeOut->Branch("mass"       ,&mass            , "mass/F");
  treeOut->Branch("ntowers"    ,&ntowersLoc      , "ntowers/F");
  treeOut->Branch("radial"     ,&shapeRadial     , "radial/F");
  treeOut->Branch("dispersion" ,&shapeDispersion , "dispersion/F");
  treeOut->Branch("ntracks"    ,&ntracks         , "ntracks/I");
  treeOut->Branch("trackleadDiff"    , &shapeLeSub       , "trackleadDiff/F");
  treeOut->Branch("trackChargeSum"   , &trackChargeSum   , "trackChargeSum/F");
  treeOut->Branch("trackAveCharge"   , &trackAveCharge   , "trackAveCharge/F");
  treeOut->Branch("trackWeightCharge", &trackWeightCharge, "trackWeightCharge/F");
  treeOut->Branch("minTrackPt"       , &minTrackPt       , "minTrackPt/F");
  treeOut->Branch("maxTrackPt"       , &maxTrackPt       , "maxTrackPt/F");
  treeOut->Branch("meanTrackPt"      , &meanTrackPt      , "meanTrackPt/F");
  treeOut->Branch("sumTrackPt"       , &sumTrackPt       , "sumTrackPt/F");
  treeOut->Branch("weightTrackPt"    , &weightTrackPt    , "weightTrackPt/F");
  treeOut->Branch("minTrackDR"       , &minTrackDR       , "minTrackDR/F");
  treeOut->Branch("maxTrackDR"       , &maxTrackDR       , "maxTrackDR/F");
  treeOut->Branch("meanTrackDR"      , &meanTrackDR      , "meanTrackDR/F");
  treeOut->Branch("minTowerE"        , &minTowerE        , "minTowerE/F");
  treeOut->Branch("maxTowerE"        , &maxTowerE        , "maxTowerE/F");
  treeOut->Branch("meanTowerE"       , &meanTowerE       , "meanTowerE/F");
  treeOut->Branch("sumTowerE"        , &sumTowerE        , "sumTowerE/F");
  treeOut->Branch("weightTowerE"     , &weightTowerE     , "weightTowerE/F");
  treeOut->Branch("minTowerEEem"     , &minTowerEEem     , "minTowerEEem/F");
  treeOut->Branch("minTowerEEhad"    , &minTowerEEhad    , "minTowerEEhad/F");
  treeOut->Branch("maxTowerEEem"     , &maxTowerEEem     , "maxTowerEEem/F");
  treeOut->Branch("maxTowerEEhad"    , &maxTowerEEhad    , "maxTowerEEhad/F");
  treeOut->Branch("minTowerEem"      , &minTowerEem      , "minTowerEem/F");
  treeOut->Branch("maxTowerEem"      , &maxTowerEem      , "maxTowerEem/F");
  treeOut->Branch("meanTowerEem"     , &meanTowerEem     , "meanTowerEem/F");
  treeOut->Branch("minTowerEhad"     , &minTowerEhad     , "minTowerEhad/F");
  treeOut->Branch("maxTowerEhad"     , &maxTowerEhad     , "maxTowerEhad/F");
  treeOut->Branch("meanTowerEhad"    , &meanTowerEhad    , "meanTowerEhad/F");


  treeOut->Branch("towerTrackRatio"  , &towerTrackRatio  , "towerTrackRatio/F");

  OpenFiles(inputDir);

  if (nentries < 0) nentries = tree->GetEntries();
  if (nentries > tree->GetEntries()) {
    std::cout << "Less entries than requested in tree: " << tree->GetEntries() << std::endl;
    nentries = tree->GetEntries();

  }
  if (nentries == 0) {
    std::cout << "nentries == 0? Please check path " << inputDir  << std::endl;

  }
  for(Int_t ientry = 0; ientry < nentries; ientry++){
    tree->GetEntry(ientry);
    if(!(ientry%10000))  {
      printf("\r Processing [%d/%d]",  ientry, nentries);
      fflush(stdout);
    }


    Float_t leadingHadronPt    = -999.;
    Float_t subleadingHadronPt = -999.;
    Float_t jetDispersionSum = 0;
    Float_t jetDispersionSquareSum = 0;
    Int_t numConst = 0;
    ntowersLoc=0;

    shapeRadial = 0.;
    trackChargeSum = 0;
    trackWeightCharge = 0;
    minTrackPt = 1000000000;
    maxTrackPt = 0;
    meanTrackPt = 0;
    sumTrackPt = 0;
    weightTrackPt = 0;
    minTrackDR = 1000000000;
    maxTrackDR = 0;
    Float_t sumTrackDR = 0;
    meanTrackDR = 0;
    minTowerE = 1000000000;
    maxTowerE = 0;
    meanTowerE = 0;
    sumTowerE = 0;
    weightTowerE = 0;
    minTowerEEem = 0;
    minTowerEEhad = 0;
    maxTowerEEem = 0;
    maxTowerEEhad = 0;
    towerTrackRatio = 0;
    minTowerEem = 1000000000;
    maxTowerEem = 0;
    Float_t sumTowerEem = 0;
    meanTowerEem = 0;
    minTowerEhad = 1000000000;
    maxTowerEhad = 0;
    Float_t sumTowerEhad = 0;
    meanTowerEhad = 0;

    for(Int_t itrack = 0; itrack < ntracks; itrack++){
      if (TMath::Abs(trackEta[itrack]) > 20.) continue;

      trackChargeSum += trackCharge[itrack];
      trackWeightCharge += trackPt[itrack]/jetPt * trackCharge[itrack];
      sumTrackPt += trackPt[itrack];
      weightTrackPt += trackPt[itrack]/jetPt * trackPt[itrack];

      // Get leading hadron pt
      if (trackPt[itrack] > leadingHadronPt){
        subleadingHadronPt = leadingHadronPt;
        leadingHadronPt    = trackPt[itrack];
      }
      else if( trackPt[itrack] > subleadingHadronPt){
        subleadingHadronPt = trackPt[itrack];
      }

      if(trackPt[itrack] < minTrackPt)
        minTrackPt = trackPt[itrack];

      Float_t deltaPhi = TMath::Min(Double_t(TMath::Abs(jetPhi-trackPhi[itrack])), Double_t(2*TMath::Pi()- TMath::Abs(jetPhi-trackPhi[itrack])));
      Float_t deltaEta = jetEta-trackEta[itrack];
      Float_t deltaR   = TMath::Sqrt(deltaPhi*deltaPhi + deltaEta*deltaEta);

      if(deltaR < minTrackDR)
        minTrackDR = deltaR;
      if(deltaR > maxTrackDR)
        maxTrackDR = deltaR;
      sumTrackDR += deltaR;

      //Calculate properties important for shape calculation
      jetDispersionSum += trackPt[itrack];
      jetDispersionSquareSum += trackPt[itrack]*trackPt[itrack];
      shapeRadial += trackPt[itrack]/jetPt * deltaR;

      numConst += 1;
    }
    // Calculate the shapes
    if (numConst > 1)
      shapeLeSub = leadingHadronPt - subleadingHadronPt;
    else
      shapeLeSub = 1.;

    if (jetDispersionSum)
      shapeDispersion = TMath::Sqrt(jetDispersionSquareSum)/jetDispersionSum;
    else
      shapeDispersion = 0.;

    if(ntracks > 0)
    {
      trackAveCharge = trackChargeSum/ntracks;
      meanTrackPt = sumTrackPt/ntracks;
      meanTrackDR = sumTrackDR/ntracks;
    }
    else
    {
      trackAveCharge = 0;
      minTrackPt = 0;
      meanTrackPt = 0;
      minTrackDR = 0;
      meanTrackDR = 0;
    }

    for(int itower = 0; itower < ntowers; ++itower)
    {
      if(towerE[itower] < minTowerE)
      {
        minTowerE = towerE[itower];
        minTowerEEem = towerEem[itower];
        minTowerEEhad = towerEhad[itower];
      }
      if(towerE[itower] > maxTowerE)
      {
        maxTowerE = towerE[itower];
        maxTowerEEem = towerEem[itower];
        maxTowerEEhad = towerEhad[itower];
      }

      if(towerEem[itower] < minTowerEem)
        minTowerEem = towerEem[itower];
      if(towerEem[itower] > maxTowerEem)
        maxTowerEem = towerEem[itower];
      if(towerEhad[itower] < minTowerEhad)
        minTowerEhad = towerEhad[itower];
      if(towerEm[itower] > maxTowerEhad)
        maxTowerEhad = towerEhad[itower];

      sumTowerE += towerE[itower];
      weightTowerE += towerE[itower]/jetPt * towerE[itower];
      sumTowerEem += towerEem[itower];
      sumTowerEhad += towerEhad[itower];
    }

    if(ntowers > 0)
    {
      meanTowerE = sumTowerE/ntowers;
      meanTowerEem = sumTowerEem/ntowers;
      meanTowerEhad = sumTowerEhad/ntowers;
    }
    else
    {
      meanTowerE = 0;
      minTowerE = 0;
      meanTowerEem = 0;
      meanTowerEhad = 0;
    }

    mass = jetMass;
    ntowersLoc = ntowers;
    maxTrackPt = leadingHadronPt;

    if(sumTrackPt > 0)
      towerTrackRatio = sumTowerE/sumTrackPt;
    else
      towerTrackRatio = 0;

    // Fill the output tree
    //    std::cout << mass << " " << ntowersLoc <<" " << shapeDispersion << " " << shapeRadial << std::endl;
    treeOut->Fill();

  }
  fOut.cd();
  treeOut->Write();
  printf("\n");
}

 void OpenFiles(const char * inputDir) {
  // Create chains and sets branck addresses
   std::cout << "Input dir: " << inputDir << std::endl;

    tree = new TChain("treeJets");
    TSystemDirectory dir(inputDir, inputDir);
    TList *files = dir.GetListOfFiles();
    if (files) {
      TSystemFile *file;
      TString fname;
      TIter next(files);
      while ((file=(TSystemFile*)next())) {
        fname = file->GetName();
        TString fnameWithPath;
        fnameWithPath.Form("%s/%s", inputDir, fname.Data());

        if (!file->IsDirectory() && fname.EndsWith(".root")) {
          //          std::cout << fnameWithPath.Data() << std::endl;
          tree->AddFile(fnameWithPath);
        }
      }
    }

    tree->SetBranchAddress("jetPt"  ,&jetPt  );
    tree->SetBranchAddress("jetEta" ,&jetEta );
    tree->SetBranchAddress("jetPhi" ,&jetPhi );
    tree->SetBranchAddress("jetMass",&jetMass);

    tree->SetBranchAddress("ntracks",&ntracks);
    tree->SetBranchAddress("ntowers",&ntowers);

    tree->SetBranchAddress("trackPt"     , trackPt    );
    tree->SetBranchAddress("trackEta"    , trackEta   );
    tree->SetBranchAddress("trackPhi"    , trackPhi   );
    tree->SetBranchAddress("trackCharge" , trackCharge);
    tree->SetBranchAddress("towerE"      , towerE     );
    tree->SetBranchAddress("towerEem"    , towerEem   );
    tree->SetBranchAddress("towerEhad"   , towerEhad  );
    tree->SetBranchAddress("towerEta"    , towerEta   );
    tree->SetBranchAddress("towerPhi"    , towerPhi   );


 }
