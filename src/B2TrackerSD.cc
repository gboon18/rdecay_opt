//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: B2TrackerSD.cc 87359 2014-12-01 16:04:27Z gcosmo $
//
/// \file B2TrackerSD.cc
/// \brief Implementation of the B2TrackerSD class

#include "B2TrackerSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

//added 20170611
#include "HistoManager.hh"
//added 20170611

//added 20170613
#include "G4Track.hh"//for G4Track*
//added 20170613

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B2TrackerSD::B2TrackerSD(const G4String& name,
                         const G4String& hitsCollectionName) 
 : G4VSensitiveDetector(name),
   fHitsCollection(NULL)
{
  collectionName.insert(hitsCollectionName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B2TrackerSD::~B2TrackerSD() 
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2TrackerSD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection

  fHitsCollection 
    = new B2TrackerHitsCollection(SensitiveDetectorName, collectionName[0]); 

  // Add this collection in hce

  G4int hcID
    = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection( hcID, fHitsCollection );

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool B2TrackerSD::ProcessHits(G4Step* aStep, 
                                     G4TouchableHistory*)
{  
  // energy deposit
  G4double edep = aStep->GetTotalEnergyDeposit();

  if (edep==0.) return false;

  B2TrackerHit* newHit = new B2TrackerHit();

  newHit->SetTrackID  (aStep->GetTrack()->GetTrackID());
  newHit->SetChamberNb(aStep->GetPreStepPoint()->GetTouchableHandle()
                                               ->GetCopyNumber());
  newHit->SetEdep(edep);
  newHit->SetPos (aStep->GetPostStepPoint()->GetPosition());

  //added 20170613
  G4Track* track = aStep->GetTrack();
  const G4ParticleDefinition* particle = track->GetParticleDefinition();//from TrackingAction.cc
  G4int pid = particle->GetPDGEncoding();
  newHit->SetPID(pid);
  //added 20170613
  
  //added 20170613_1
  G4double ini_energy = track->GetVertexKineticEnergy();
  newHit->SetIniKinEnergy(ini_energy);
  //added 20170613_1

  //added 20170614
   G4String ParticleName = track->GetDynamicParticle()->GetParticleDefinition()->GetParticleName();
   //    if(ParticleName == "opticalphoton")
      // G4cout<<ParticleName<<G4endl;
   newHit->SetOP(ParticleName);
  //added 20170614

  
  fHitsCollection->insert( newHit );
  //G4cout<<"edep : "<<edep<<G4endl;
  //  newHit->Print();//turn this on, B2TrackerHit::Print is operated.

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2TrackerSD::EndOfEvent(G4HCofThisEvent*)
{
  /*useless
  //added 2017612
  B2TrackerHit* newHit = new B2TrackerHit();
  G4int TrackID = newHit->GetTrackID();
  G4int ChamberNb = newHit->GetChamberNb();
  G4double Edep = newHit->GetEdep();G4cout<<"Edep : "<<Edep<<G4endl;
  G4ThreeVector Pos = newHit->GetPos();
  */
  //added 20170612
  
  /*moved to B2TrackerHit.cc
  //added 20170611
  G4int id = 27;
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->FillNtupleDColumn(id,0,Edep);
  analysisManager->FillNtupleDColumn(id,1,Pos.x());
  analysisManager->FillNtupleDColumn(id,2,Pos.y());
  analysisManager->FillNtupleDColumn(id,3,Pos.z());
  analysisManager->FillNtupleDColumn(id,4,TrackID);
  analysisManager->FillNtupleDColumn(id,5,ChamberNb);  
  analysisManager->AddNtupleRow(id);
  //added 20170611
  */
  if ( verboseLevel>1 ) { 
     G4int nofHits = fHitsCollection->entries();
     G4cout << G4endl
            << "-------->Hits Collection: in this event they are " << nofHits 
            << " hits in the tracker chambers: " << G4endl;
     for ( G4int i=0; i<nofHits; i++ ) (*fHitsCollection)[i]->Print();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
