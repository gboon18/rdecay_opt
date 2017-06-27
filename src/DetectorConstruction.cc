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
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
// $Id: DetectorConstruction.cc 70755 2013-06-05 12:17:48Z ihrivnac $
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
//added 20170611
#include "G4SDManager.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "B2TrackerSD.hh"
#include "G4AutoDelete.hh"
//added 20170611

#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4RunManager.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4UnitsTable.hh"

//added 20170128
#include "G4Box.hh"
//added 20170607
#include "G4SubtractionSolid.hh"

//added 20170614
#include "G4MaterialTable.hh"
//added 20170614


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//added 20170611
G4ThreadLocal
G4GlobalMagFieldMessenger* DetectorConstruction::fMagFieldMessenger = 0;
//added 20170611


DetectorConstruction::DetectorConstruction()
:G4VUserDetectorConstruction(),
 fTargetMater(0), fLogicTarget(0),
 fDetectorMater(0), fLogicDetector(0),
 fReflectorMater(0),//added20170607
 fCoverMater(0),//added20170607
 lReflector(0),//added20170607
 lCover(0),//added20170607
 lCryst(0),//added20170128 
 fWorldMater(0), fPhysiWorld(0),
 /*
 fCheckOverlaps(true)
{
  DefineMaterials();
  },*/
 fDetectorMessenger(0)
{
  ////////////////We Do not actually use these20170129// 
  fTargetLength      = 4.6*cm; 
  fTargetRadius      = 2.3*cm;
  fDetectorLength    = 5*cm; 
  fDetectorThickness = 2*cm;
  //////////////////////////////////////////////////////
  DefineMaterials();
  fDetectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ delete fDetectorMessenger;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  // build materials
  //
  /*
  fDetectorMater = 
  new G4Material("Germanium", 32, 72.61*g/mole, 5.323*g/cm3);
  */

  
  G4double z, a, fractionmass, density;
  G4String name, symbol;
  G4int ncomponents;//# of components in the material.

  a = 22.99*g/mole;
  G4Element* Na = new G4Element(name = "Sodium", symbol = "Na", z = 11, a);

  a = 126.9*g/mole;
  G4Element* I = new G4Element(name = "Iodine", symbol = "I", z = 53, a);


  density = 3.67*g/cm3;
  //  G4Material* NaI

  NaI = new G4Material(name = "NaI",density,ncomponents=2);
  NaI->AddElement(Na, fractionmass=0.5);
  NaI->AddElement(I, fractionmass=0.5);


  ///////Detector Material/////////////////////////////////////////////////////////////////
  fDetectorMater = NaI;//20170129 NaI as a detector material, NaI is also provided by NIST.
  /////////////////////////////////////////////////////////////////////////////////////////
  G4Element* N  = new G4Element("Nitrogen", "N", 7, 14.01*g/mole);
  G4Element* O  = new G4Element("Oxygen",   "O", 8, 16.00*g/mole);
  G4Element* C = new G4Element("Carbon", "C", 6, 12.01*g/mole);
  G4Element* F = new G4Element("Fluorine", "F", 9, 18.9984*g/mole);
  //

  //20170607//Al, MgO
  density = 2.700*g/cm3;
  a = 26.98*g/mole;
  G4Material* Al = new G4Material(name = "Aluminum", z=13., a, density);
  fCoverMater = Al;
  a = 24.31*g/mole;
  G4Element* Mg = new G4Element(name = "Magnesium", symbol = "Mg", z = 12, a);

  
  density = 2.*g/cm3;
  G4Material* MgO
    = new G4Material(name = "MgO", density, ncomponents = 2);
  MgO->AddElement(Mg, fractionmass = 0.5);
  MgO->AddElement(O, fractionmass = 0.5);
  fReflectorMater = MgO;
  /*
  density = 2.2*g/cm3;
  G4Material* Teflon = new G4Material("Teflon", density, ncomponents = 2);
  Teflon->AddElement(C,2./6.);
  Teflon->AddElement(F,4./6.);
  fReflectorMater = Teflon;//20170608 described as MgO in this script.
  */
  //20170607//
 
  G4Material* Air20 = new G4Material("Air", 1.205*mg/cm3, ncomponents=2,
                      kStateGas, 293.*kelvin, 1.*atmosphere);
    Air20->AddElement(N, fractionmass=0.7);
    Air20->AddElement(O, fractionmass=0.3);
  //
    /////////World Material///////////////////////////////////////////////////////////////
  fWorldMater = Air20;
  /////////////////////////////////////////////////////////////////////////////////////////
  // or use G4 materials data base
  //
  // G4NistManager* man = G4NistManager::Instance();  
  //  fTargetMater = man->FindOrBuildMaterial("G4_CESIUM_IODIDE");

 
  
  /////Target Material////////////////////////////////////////////////////////////////////
  // fTargetMater = new G4Material("Silicon", 14, 28.09*g/mole, 2.3290*g/cm3);//20170608 deleted
  ////////////////////////////////////////////////////////////////////////////////////////////
  //20170608
  G4Material* Si = new G4Material(name = "Silicon", z = 14, 28.09*g/mole, (10./19.)*2.3290*g/cm3);//densidy is reduced by 10/19(Si strip/Air gap raio)
  fTargetMater = Si;
  //20170608
  
 ///G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{
  //added 20170614
  const G4int nEntries = 32;

  G4double PhotonEnergy[nEntries] =
            { 2.034*eV, 2.068*eV, 2.103*eV, 2.139*eV,
              2.177*eV, 2.216*eV, 2.256*eV, 2.298*eV,
              2.341*eV, 2.386*eV, 2.433*eV, 2.481*eV,
              2.532*eV, 2.585*eV, 2.640*eV, 2.697*eV,
              2.757*eV, 2.820*eV, 2.885*eV, 2.954*eV,
              3.026*eV, 3.102*eV, 3.181*eV, 3.265*eV,
              3.353*eV, 3.446*eV, 3.545*eV, 3.649*eV,
              3.760*eV, 3.877*eV, 4.002*eV, 4.136*eV };
 
//
// Water (NaI)
//
/*
   G4double RefractiveIndex1[nEntries]=
            {1.8390,   1.8395,   1.8400,   1.8405,
             1.8410,   1.8415,   1.8420,   1.8425,
             1.8430,   1.8435,   1.8440,   1.8447,
             1.8455,   1.8460,   1.8465,   1.8473,
             1.8477,   1.8485,   1.8490,   1.8495,
             1.8500,   1.8505,   1.8510,   1.8515,
             1.8523,   1.8527,   1.8535,   1.8540,
             1.8545,   1.8550,   1.8555,   1.8563 };
*/
/*
   G4double RefractiveIndex1[nEntries]=
            {	1.766488159,1.7679973633,1.769572262,1.7712157643,
		1.772977708,1.7748162458,1.7767350536,1.7787873582,
		1.7809299411,1.7832187719,1.7856621319,1.7882152532,
		1.7909943174,1.7939575559,1.797116742,1.8004847152,
		1.804137281,1.8080957536,1.8123177871,1.8169595385,
		1.8219866195,1.8275060062,1.8334867877,1.8401328865,
		1.8474304253,1.8555377073,1.864643371,1.8747742158
		,1.8862752549,1.8992321063,1.9141084948,1.9313637927};
    
  G4double Absorption1[nEntries] =
           {3.448*m,  4.082*m,  6.329*m,  9.174*m, 12.346*m, 13.889*m,
           15.152*m, 17.241*m, 18.868*m, 20.000*m, 26.316*m, 35.714*m,
           45.455*m, 47.619*m, 52.632*m, 52.632*m, 55.556*m, 52.632*m,
           52.632*m, 47.619*m, 45.455*m, 41.667*m, 37.037*m, 33.333*m,
           30.000*m, 28.500*m, 27.000*m, 24.500*m, 22.000*m, 19.500*m,
           17.500*m, 14.500*m };
*/
  G4double ScintilFast[nEntries] =
            { 0.000134, 0.004432, 0.053991, 0.241971, 0.398942, 
		0.000134, 0.004432,0.053991, 0.241971, 0.398942,
		0.000134, 0.004432, 0.053991, 0.241971, 0.398942,
		0.000134, 0.004432, 0.053991, 0.241971, 0.398942,
		0.000134, 0.004432, 0.053991, 0.241971, 0.398942,
		0.000134, 0.004432, 0.053991, 0.241971, 0.398942,
		0.000134, 0.004432 };
  G4double ScintilSlow[nEntries] =
            { 0.000010, 0.000020, 0.000030, 0.004000, 0.008000, 0.005000, 0.020000,0.001000,
		0.000010, 0.000020, 0.000030, 0.004000, 0.008000, 0.005000, 0.020000, 0.001000,
		0.000010, 0.000020, 0.000030, 0.004000, 0.008000, 0.005000, 0.020000,0.001000,
		0.000010, 0.000020, 0.000030, 0.004000, 0.008000, 0.005000, 0.020000,0.001000,};

  G4MaterialPropertiesTable* myMPT1 = new G4MaterialPropertiesTable();
  // myMPT1->AddProperty("RINDEX",       PhotonEnergy, RefractiveIndex1,nEntries);
  // myMPT1->AddProperty("ABSLENGTH",    PhotonEnergy, Absorption1,     nEntries);
  myMPT1->AddProperty("FASTCOMPONENT",PhotonEnergy, ScintilFast,     nEntries);
  myMPT1->AddProperty("SLOWCOMPONENT",PhotonEnergy, ScintilSlow,     nEntries);

  myMPT1->AddConstProperty("SCINTILLATIONYIELD",3800./MeV);
  //  myMPT1->AddConstProperty("RESOLUTIONSCALE",1.0);
  myMPT1->AddConstProperty("RESOLUTIONSCALE",2.0);  
  myMPT1->AddConstProperty("FASTTIMECONSTANT", 1.*ns);
  myMPT1->AddConstProperty("SLOWTIMECONSTANT",10.*ns);
  myMPT1->AddConstProperty("YIELDRATIO",0.8);

  NaI->SetMaterialPropertiesTable(myMPT1);
//added 20170614




  // Cleanup old geometry
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

 G4int array = 3; //number of arrays
 //G4double fDetectorX = (150. - 1. - 1.)*mm, fDetectorY = (150. - 1. - 1.)*mm, fDetectorZ = (250. - 1. - 1.)*mm;//size of each detector//20170607, 1mm for the Al cover, 1mm for the MgO reflector.
 G4double fDetectorX = (150.)*mm, fDetectorY = (150.)*mm, fDetectorZ = (250.)*mm;//size of each detector//20170607, 1mm for the Al cover, 1mm for the MgO reflector.//20170608

G4double gap = 0.5*mm; //a gap for wrapping
G4double gapAl = 1.*mm; // Al cover for the NaI//20170607
G4double gapMgO = 1.*mm; // MgO reflector for the NaI//20170607

 G4double dX = fDetectorX + gap + gapAl + gapMgO, dY = fDetectorY + gap + gapAl + gapMgO, dZ = fDetectorZ + gap + gapAl + gapMgO;//20170607
  
G4double  fDetectorWholeX = (dX)*array;//Whole size of the big detector
G4double  fDetectorWholeY = (dY)*array;//Whole size of the big detector
G4double  fDetectorWholeZ = (dZ)*(array-1);//Whole size of the big detector//in z-axis, array is 2 20170131


  
  //     
  // World//from B3a, 20170128
  //
  G4double world_sizeXY = 2.4*fDetectorWholeX;
  G4double world_sizeZ  = 1.2*fDetectorWholeZ;
  
  G4Box* sWorld =
    new G4Box("World",0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ); //its size
 			   
  G4LogicalVolume* lWorld =                         
    new G4LogicalVolume(sWorld,          //its solid
                        fWorldMater,         //its material
                        "World");            //its name
                                   
  fPhysiWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      lWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      fCheckOverlaps);       // checking overlaps 
  ////////////////////////////////////////////////////////////////////
  
  // Target
  //
  //20170227 no target. only the detectors.
  //20170708 AIDA
  G4double fTargetX = 10.*cm, fTargetY = 10.*cm, fTargetZ = 19.*mm;//20170608
  fTargetRadius = (fDetectorX-1*cm)*0.5;//20170608 impotent
  fTargetLength = fDetectorZ-1*cm;//20170608 impotent
  /*
  G4Tubs* 
  sTarget = new G4Tubs("Target",                                   //name
                  0., fTargetRadius, 0.5*fTargetLength, 0.,twopi); //dimensions


  fLogicTarget = new G4LogicalVolume(sTarget,           //shape
                             fTargetMater,              //material
                             "Target");                 //name
                               
           new G4PVPlacement(0,                         //no rotation
                           G4ThreeVector(),             //at (0,0,0)
                           fLogicTarget,                //logical volume
                           "Target",                    //name
                           lWorld,                      //mother  volume//20170128 originally it was lWorld
                           false,                       //no boolean operation
                           0);                          //copy number
  */
  G4Box* 
  sTarget = new G4Box("Target",                                   //name
		      0.5*fTargetX,
		      0.5*fTargetY,
		      0.5*fTargetZ
		      ); //dimensions


  fLogicTarget = new G4LogicalVolume(sTarget,           //shape
                             fTargetMater,              //material
                             "Target");                 //name
                               
           new G4PVPlacement(0,                         //no rotation
                           G4ThreeVector(),             //at (0,0,0)
                           fLogicTarget,                //logical volume
                           "Target",                    //name
                           lWorld,                      //mother  volume//20170128 originally it was lWorld
                           false,                       //no boolean operation
                           0);                          //copy number
    	   
  //Small Each Detector
  //
  G4Box* 
  sCryst = new G4Box("crystal",  
			0.5*fDetectorX,0.5*fDetectorY,0.5*fDetectorZ);


  lCryst = new G4LogicalVolume(sCryst,       //shape
                             fDetectorMater,            //material
                             "crystalLV");               //name

  
  //20170607 MgO reflector
  G4Box*
    rawReflector = new G4Box("rawReflector",
			  0.5*(fDetectorX + gapMgO),0.5*(fDetectorY + gapMgO),0.5*(fDetectorZ + gapMgO));
  G4SubtractionSolid* subReflector = new G4SubtractionSolid("Reflector",
							    rawReflector,sCryst);
  lReflector = new G4LogicalVolume(subReflector,
				  fReflectorMater,
  				  "reflector");
  //20170607
  //20170607 Al cover
  G4Box*
    rawCover = new G4Box("rawCover",
			  0.5*(fDetectorX + gapMgO +gapAl),0.5*(fDetectorY + gapMgO +gapAl),0.5*(fDetectorZ + gapMgO +gapAl));
  G4SubtractionSolid* subCover = new G4SubtractionSolid("Cover",
							    rawCover,rawReflector);
  lCover = new G4LogicalVolume(subCover,
			      fCoverMater,
			      "cover");
  //20170607
  
  //Big Whole Detector 20170128 //20170607 modified for the Al cover, MgO reflector.

  
  G4Box* sDetector = new G4Box("Detector", 0.5*fDetectorWholeX, 0.5*fDetectorWholeY, 0.5*fDetectorWholeZ);
  fLogicDetector = new G4LogicalVolume(sDetector,fWorldMater,"Detector");

  //now, place small detectors in the big whole detector
  G4double fPositionX = -0.5*(fDetectorWholeX+dX);
  G4double fPositionY;
  G4double fPositionZ;
  // G4double fPositionY = -0.5*(fDetectorWholeY+fDetectorY);
  // G4double fPositionZ = -0.5*(fDetectorWholeZ+fDetectorZ);
  //G4int x;
  G4int y;
  G4int z;
  G4int i = 0;// CopyNo.
  
 
        fPositionX = fPositionX + dX;
        fPositionY = -0.5*(fDetectorWholeY+dY);
        fPositionZ = -0.5*(fDetectorWholeZ+dZ);
	for(y = 0 ; y < array ; y++) {
            fPositionY = fPositionY + dY;
            fPositionZ = -0.5*(fDetectorWholeZ+dZ);
	    for(z = 0 ; z < array-1 ; z++) {//array is 2 in z-axis 20170131
		    fPositionZ = fPositionZ + dZ;
	    new G4PVPlacement(0,
		              G4ThreeVector(fPositionX,fPositionY,fPositionZ),
			      lCryst,
			      "s_crystal",
			      lWorld,
			      false,
			      i,//copy number
			      fCheckOverlaps);
	    new G4PVPlacement(0,//20170607 MgO reflector
		              G4ThreeVector(fPositionX,fPositionY,fPositionZ),
			      lReflector,
			      "s_reflector",
			      lWorld,
			      false,
			      i,//copy number
			      fCheckOverlaps);
	    new G4PVPlacement(0,//20170607 Al cover
		              G4ThreeVector(fPositionX,fPositionY,fPositionZ),
			      lCover,
			      "s_cover",
			      lWorld,
			      false,
			      i,//copy number
			      fCheckOverlaps);
	      i++;
	      }
	}
	
	
        fPositionX = -0.5*(fDetectorWholeX+dX);
        fPositionY = -0.5*(fDetectorWholeY+dY);
        fPositionZ = -0.5*(fDetectorWholeZ+dZ);
        fPositionX = fPositionX + 3*dX;
	for(y = 0 ; y < array ; y++) {
            fPositionY = fPositionY + dY;
            fPositionZ = -0.5*(fDetectorWholeZ+dZ);
	    for(z = 0 ; z < array-1 ; z++) {//array is 2 in z-axis 20170131
		    fPositionZ = fPositionZ + dZ;
	    new G4PVPlacement(0,
		              G4ThreeVector(fPositionX,fPositionY,fPositionZ),
			      lCryst,
			      "s_crystal",
			      lWorld,
			      false,
			      i,
			      fCheckOverlaps);
	    new G4PVPlacement(0,//20170607 MgO reflector
		              G4ThreeVector(fPositionX,fPositionY,fPositionZ),
			      lReflector,
			      "s_reflector",
			      lWorld,
			      false,
			      i,//copy number
			      fCheckOverlaps);
	    new G4PVPlacement(0,//20170607 Al cover
		              G4ThreeVector(fPositionX,fPositionY,fPositionZ),
			      lCover,
			      "s_cover",
			      lWorld,
			      false,
			      i,//copy number
			      fCheckOverlaps);	    
	      i++;
	      }
	}

	
	//avoid the target detector area and the beam coming area////////////////////// 

	fPositionX = -0.5*(fDetectorWholeX+dX);
        fPositionY = -0.5*(fDetectorWholeY+dY);
        fPositionZ = -0.5*(fDetectorWholeZ+dZ);
	fPositionX = fPositionX + 2*dX;
	fPositionY = fPositionY + 3*dY;
	for(z = 0 ; z < array - 1 ; z++) {
	  fPositionZ = fPositionZ + dZ;
		  	    new G4PVPlacement(0,
		              G4ThreeVector(fPositionX,fPositionY,fPositionZ),
			      lCryst,
			      "s_crystal",
			      lWorld,
			      false,
			      i,//copy number
			      fCheckOverlaps);
    	                    new G4PVPlacement(0,//20170607 MgO reflector
		              G4ThreeVector(fPositionX,fPositionY,fPositionZ),
			      lReflector,
			      "s_reflector",
			      lWorld,
			      false,
			      i,//copy number
			      fCheckOverlaps);
	                    new G4PVPlacement(0,//20170607 Al cover
		              G4ThreeVector(fPositionX,fPositionY,fPositionZ),
			      lCover,
			      "s_cover",
			      lWorld,
			      false,
			      i,//copy number
			      fCheckOverlaps);			    
		    i++;
	}

	fPositionX = -0.5*(fDetectorWholeX+dX);
        fPositionY = -0.5*(fDetectorWholeY+dY);
        fPositionZ = -0.5*(fDetectorWholeZ+dZ);
	fPositionX = fPositionX + 2*dX;
	fPositionY = fPositionY + dY;
      	for(z = 0 ; z < array - 1 ; z++) {
	  fPositionZ = fPositionZ + dZ;
		  	    new G4PVPlacement(0,
		              G4ThreeVector(fPositionX,fPositionY,fPositionZ),
			      lCryst,
			      "s_crystal",
			      lWorld,
			      false,
			      i,
			      fCheckOverlaps);
	                    new G4PVPlacement(0,//20170607 MgO reflector
		              G4ThreeVector(fPositionX,fPositionY,fPositionZ),
			      lReflector,
			      "s_reflector",
			      lWorld,
			      false,
			      i,//copy number
			      fCheckOverlaps);
	                    new G4PVPlacement(0,//20170607 Al cover
		              G4ThreeVector(fPositionX,fPositionY,fPositionZ),
			      lCover,
			      "s_cover",
			      lWorld,
			      false,
			      i,//copy number
			      fCheckOverlaps);
			    i++;
	}
	  
	/*	   
	   
	/////////////////////////////////////////////////////////////////////////////////////////
  new G4PVPlacement(0,
		    G4ThreeVector(),
	            fLogicDetector,
		    "Detector",
		    lWorld,
		    false,
		    0,
		    fCheckOverlaps);
			      
  //////////////////////////////////////////////////			      
  */

  
  
  PrintParameters();
  
  //always return the root volume
  //
  return fPhysiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//added 20170611
void DetectorConstruction::ConstructSDandField()
{
  // Sensitive detectors

  G4String trackerChamberSDname = "B2/TrackerChamberSD";
  B2TrackerSD* aTrackerSD = new B2TrackerSD(trackerChamberSDname,
                                            "TrackerHitsCollection");
  /*
  G4String trackerChamberSDname = "rdecay02/NaISD";
  B2TrackerSD* aTrackerSD = new B2TrackerSD(trackerChamberSDname,
                                            "NaIHitsCollection");
  */
  G4SDManager::GetSDMpointer()->AddNewDetector(aTrackerSD);
  // Setting aTrackerSD to all logical volumes with the same name                                                       
  // of "Chamber_LV"-->"crystalLV".                                                                                                   
  SetSensitiveDetector("crystalLV", aTrackerSD, true);

  // Create global magnetic field messenger.                                                                            
  // Uniform magnetic field is then created automatically if                                                            
  // the field value is not zero.                                                                                       
  G4ThreeVector fieldValue = G4ThreeVector();
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);

  // Register the field messenger for deleting                                                                          
  G4AutoDelete::Register(fMagFieldMessenger);
}
//added 20170611

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters()
{
  G4cout << "\n Target : Length = " << G4BestUnit(fTargetLength,"Length")
         << " Radius = " << G4BestUnit(fTargetRadius,"Length")  
         << " Material = " << fTargetMater->GetName();
  G4cout << "\n Detector : Length = " << G4BestUnit(fDetectorLength,"Length")
         << " Tickness = " << G4BestUnit(fDetectorThickness,"Length")  
         << " Material = " << fDetectorMater->GetName() << G4endl;          
  G4cout << "\n" << fTargetMater << "\n" << fDetectorMater << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetTargetMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
     G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);   
  
  if (pttoMaterial) { 
    fTargetMater = pttoMaterial;
    if(fLogicTarget) { fLogicTarget->SetMaterial(fTargetMater); }
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetTargetMaterial : "
           << materialChoice << " not found" << G4endl;
  }              
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetDetectorMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
     G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);   
  
  if (pttoMaterial) { 
    fDetectorMater = pttoMaterial;
    if(fLogicDetector) { fLogicDetector->SetMaterial(fDetectorMater); }
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetDetectorMaterial : "
           << materialChoice << " not found" << G4endl;
  }              
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetTargetRadius(G4double value)
{
  fTargetRadius = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetTargetLength(G4double value)
{
  fTargetLength = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetDetectorThickness(G4double value)
{
  fDetectorThickness = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetDetectorLength(G4double value)
{
  fDetectorLength = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetTargetLength()
{
  return fTargetLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetTargetRadius()
{
  return fTargetRadius;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* DetectorConstruction::GetTargetMaterial()
{
  return fTargetMater;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LogicalVolume* DetectorConstruction::GetLogicTarget()
{
  return fLogicTarget;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetDetectorLength()
{
  return fDetectorLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetDetectorThickness()
{
  return fDetectorThickness;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* DetectorConstruction::GetDetectorMaterial()
{
  return fDetectorMater;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LogicalVolume* DetectorConstruction::GetLogicDetector()
{
  // return fLogicDetector;//deleted 20170128
   return lCryst;//20170128
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
