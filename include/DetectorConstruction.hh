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
/// \file DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class
//
// $Id: DetectorConstruction.hh 66586 2012-12-21 10:48:39Z ihrivnac $
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4LogicalVolume;
class G4Material;
class DetectorMessenger;

//added 20170611
class G4GlobalMagFieldMessenger;
//added 20170611
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    DetectorConstruction();
   ~DetectorConstruction();

  public:
  
    virtual G4VPhysicalVolume* Construct();

  //added 20170611
    virtual void ConstructSDandField();
  //added 20170611
  
    void SetTargetLength (G4double value);
    void SetTargetRadius (G4double value);
    void SetTargetMaterial (G4String);
    
    void SetDetectorLength(G4double value);           
    void SetDetectorThickness(G4double value);  
    void SetDetectorMaterial(G4String);               
                   
    void PrintParameters();
    
  public:
      
    G4double GetTargetLength();
    G4double GetTargetRadius();
    G4Material* GetTargetMaterial();       
    G4LogicalVolume* GetLogicTarget();
    
    G4double GetDetectorLength();
    G4double GetDetectorThickness();
    G4Material* GetDetectorMaterial();                 
    G4LogicalVolume* GetLogicDetector();      
                       
  private:
  
    G4double           fTargetLength; 
    G4double           fTargetRadius;
    G4Material*        fTargetMater;
    G4LogicalVolume*   fLogicTarget;
                 
    G4double           fDetectorLength;
    G4double           fDetectorThickness;
    G4Material*        fDetectorMater;
    G4LogicalVolume*   fLogicDetector;

  //20170607
    G4Material*        fReflectorMater;
    G4Material*        fCoverMater;
    G4LogicalVolume*   lReflector;
    G4LogicalVolume*   lCover;
  //20170607

  
    G4LogicalVolume*   lCryst;//added20170128
  ///////////20170128 deleted             
  //   G4double           fWorldLength;
  //   G4double           fWorldRadius;
  /////////////////////
    G4Material*        fWorldMater;     
    G4VPhysicalVolume* fPhysiWorld;
                
    DetectorMessenger* fDetectorMessenger;

  //added 20170614
  G4Material* NaI;
  //added 20170614

  
  //added 20170611
  static G4ThreadLocal G4GlobalMagFieldMessenger* fMagFieldMessenger;
  //added 20170611
  
  private:
    
    void               DefineMaterials();
    G4VPhysicalVolume* ConstructVolumes();     
    G4bool fCheckOverlaps; //added 20170128
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#endif

