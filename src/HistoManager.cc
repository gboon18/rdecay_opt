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
/// \file HistoManager.cc
/// \brief Implementation of the HistoManager class
//
// $Id: HistoManager.cc 83882 2014-09-22 11:09:30Z maire $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "HistoManager.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager()
  : fFileName("rdecay02")
{
  Book();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager()
{
  delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Book()
{
  // Create or get analysis manager
  // The choice of analysis technology is done via selection of a namespace
  // in HistoManager.hh
  //
  G4AnalysisManager* analysis = G4AnalysisManager::Instance();
  
  analysis->SetFileName(fFileName);
  analysis->SetVerboseLevel(1);
  analysis->SetActivation(true);     //enable inactivation of histos, nTuples
    
  // Default values (to be reset via /analysis/h1/set command)               
  G4int nbins = 100;
  G4double vmin = 0.;
  G4double vmax = 100.;

  // Create all histograms as inactivated 
  // as we have not yet set nbins, vmin, vmax
  //
  ////analysis->SetHistoDirectoryName("histo");  
  ////analysis->SetFirstHistoId(1);
  analysis->SetNtupleMerging(true);//added to merge the outputfile 20170127    
  G4int id = analysis->CreateH1("H10","Energy deposit (MeV) in the target", //histo id=0
                       nbins, vmin, vmax);
  analysis->SetH1Activation(id, false);
    
  id = analysis->CreateH1("H11","Energy deposit (MeV) in the detector",//id=1
                 nbins, vmin, vmax);
  analysis->SetH1Activation(id, false);

  id = analysis->CreateH1("H12","Total energy (MeV) in target and detector",//id=2)
                 nbins, vmin, vmax);
  analysis->SetH1Activation(id, false);

  id = analysis->CreateH1("H13",
			  "Coincidence spectrum (MeV) between the target and detector",//id=3
                 nbins, vmin, vmax);
  analysis->SetH1Activation(id, false);  

  id = analysis->CreateH1("H14",
			  "Anti-coincidence spectrum (MeV) in the traget",//id=4
                 nbins, vmin, vmax);
  analysis->SetH1Activation(id, false);

  id = analysis->CreateH1("H15",
			  "Anti-coincidence spectrum (MeV) in the detector",//id=5
                 nbins, vmin, vmax);
  analysis->SetH1Activation(id, false);  

  id = analysis->CreateH1("H16","Decay emission spectrum (MeV)",//id=6
                 nbins, vmin, vmax);
  analysis->SetH1Activation(id, false);
  //20170223//20170224////////////////////////////////////////////////////////////////////
  for(G4int i = 0 ; i < 16 ; i++){
    // std::string text1 = "0";
    std::string text1 = "H_";
    std::string text2 = "Energy deposit (MeV) in the detector CopyNo ";
    text1 += std::to_string(i);
    text2 += std::to_string(i);
    id = analysis->CreateH1(text1,text2,//histo id=7~
                 nbins, vmin, vmax);
  analysis->SetH1Activation(id, false);
  }
  
  ////////////////////////////////////////////////////////////////////////////////
  // nTuples
  //
  ////analysis->SetNtupleDirectoryName("ntuple");
  ////analysis->SetFirstNtupleId(1);
  //       
  analysis->CreateNtuple("T1", "Emitted Particles");//id=0 <-- this id indication is not sure....20170131
  analysis->CreateNtupleDColumn("PID");       //column 0
  analysis->CreateNtupleDColumn("Energy");    //column 1
  analysis->CreateNtupleDColumn("Time");      //column 2
  analysis->CreateNtupleDColumn("Weight");    //column 3
  /*
  //20170202/////////////////////////////////////////////////
  analysis->CreateNtupleDColumn("Z");         //column 4
  analysis->CreateNtupleDColumn("A");         //column 5    
  ////////////////////////////////////////////////////////////
  */
  analysis->FinishNtuple();
  
  analysis->CreateNtuple("T2", "RadioIsotopes");//id=1
  analysis->CreateNtupleDColumn("PID");       //column 0
  analysis->CreateNtupleDColumn("Time");      //column 1
  analysis->CreateNtupleDColumn("Weight");    //column 2
  analysis->FinishNtuple();
  
  analysis->CreateNtuple("T3", "Energy depositions");//id=2
  analysis->CreateNtupleDColumn("Energy");    //column 0
  analysis->CreateNtupleDColumn("Time");      //column 1
  analysis->CreateNtupleDColumn("Weight");    //column 2
  analysis->FinishNtuple();
  
  analysis->CreateNtuple("RDecayProducts", "All Products of RDecay");//id=3
  analysis->CreateNtupleDColumn("PID");       //column 0
  analysis->CreateNtupleDColumn("Z");         //column 1
  analysis->CreateNtupleDColumn("A");         //column 2    
  analysis->CreateNtupleDColumn("Energy");    //column 3
  analysis->CreateNtupleDColumn("Time");      //column 4
  analysis->CreateNtupleDColumn("Weight");    //column 5
  analysis->FinishNtuple();

  /////////////////////////////////////////////////////20170131 to make each detector work(not really, this saves {almost}everything)
  analysis->CreateNtuple("T101", "EveryStep");//id=4
  analysis->CreateNtupleDColumn("PID");       //column 0
  analysis->CreateNtupleDColumn("Z");         //column 1
  analysis->CreateNtupleDColumn("A");         //column 2    
  analysis->CreateNtupleDColumn("Energy");    //column 3
  analysis->CreateNtupleDColumn("Time");      //column 4
  analysis->CreateNtupleDColumn("x");
  analysis->CreateNtupleDColumn("y");
  analysis->CreateNtupleDColumn("z");
  analysis->CreateNtupleDColumn("px");
  analysis->CreateNtupleDColumn("py");
  analysis->CreateNtupleDColumn("pz");
  analysis->CreateNtupleDColumn("trackID");
  analysis->CreateNtupleDColumn("StepNo");
  //analysis->CreateNtupleDColumn("testc");
  analysis->CreateNtupleDColumn("Weight");    
  /*
  analysis->CreateNtupleDColumn("x0");
  analysis->CreateNtupleDColumn("y0");
  analysis->CreateNtupleDColumn("z0");
  */
  analysis->CreateNtupleDColumn("PatrackID");//20170207 //trackID of the parent particle
  analysis->CreateNtupleDColumn("EnergyDep");   
  analysis->FinishNtuple();
  /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////20170201 to make each detector work
  analysis->CreateNtuple("T001", "Detector0_0_Pre");//id=5, presteppoint
  analysis->CreateNtupleDColumn("PID");       //column 0
  analysis->CreateNtupleDColumn("Z");         //column 1
  analysis->CreateNtupleDColumn("A");         //column 2    
  analysis->CreateNtupleDColumn("Energy");    //column 3
  analysis->CreateNtupleDColumn("Time");      //column 4
  analysis->CreateNtupleDColumn("x");
  analysis->CreateNtupleDColumn("y");
  analysis->CreateNtupleDColumn("z");
  analysis->CreateNtupleDColumn("px");
  analysis->CreateNtupleDColumn("py");
  analysis->CreateNtupleDColumn("pz");
  analysis->CreateNtupleDColumn("trackID");
  analysis->CreateNtupleDColumn("StepNo");
  analysis->CreateNtupleDColumn("Weight");
  analysis->CreateNtupleDColumn("CopyNo");    
  analysis->CreateNtupleDColumn("EnergyDep");    
  analysis->FinishNtuple();
  /////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////20170201 for the Target
  analysis->CreateNtuple("T002", "Target0_0_Pre");//id=6, presteppoint
  analysis->CreateNtupleDColumn("PID");       //column 0
  analysis->CreateNtupleDColumn("Z");         //column 1
  analysis->CreateNtupleDColumn("A");         //column 2    
  analysis->CreateNtupleDColumn("Energy");    //column 3
  analysis->CreateNtupleDColumn("Time");      //column 4
  analysis->CreateNtupleDColumn("x");
  analysis->CreateNtupleDColumn("y");
  analysis->CreateNtupleDColumn("z");
  analysis->CreateNtupleDColumn("px");
  analysis->CreateNtupleDColumn("py");
  analysis->CreateNtupleDColumn("pz");
  analysis->CreateNtupleDColumn("trackID");
  analysis->CreateNtupleDColumn("StepNo");
  analysis->CreateNtupleDColumn("Weight");
  analysis->CreateNtupleDColumn("CopyNo");    
  analysis->CreateNtupleDColumn("EnergyDep");    
  analysis->FinishNtuple();
  /////////////////////////////////////////////////////////////////////////////////////////////
       /////////////////////////////////////////////////////20170201 to make each detector work
  analysis->CreateNtuple("T003", "Detector0_0_Post");//id=7, poststeppoint
  analysis->CreateNtupleDColumn("PID");       //column 0
  analysis->CreateNtupleDColumn("Z");         //column 1
  analysis->CreateNtupleDColumn("A");         //column 2    
  analysis->CreateNtupleDColumn("Energy");    //column 3
  analysis->CreateNtupleDColumn("Time");      //column 4
  analysis->CreateNtupleDColumn("x");
  analysis->CreateNtupleDColumn("y");
  analysis->CreateNtupleDColumn("z");
  analysis->CreateNtupleDColumn("px");
  analysis->CreateNtupleDColumn("py");
  analysis->CreateNtupleDColumn("pz");
  analysis->CreateNtupleDColumn("trackID");
  analysis->CreateNtupleDColumn("StepNo");
  analysis->CreateNtupleDColumn("Weight");
  analysis->CreateNtupleDColumn("CopyNo");    
  analysis->CreateNtupleDColumn("EnergyDep");    
  analysis->FinishNtuple();
   /////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////20170201 for the Target
  analysis->CreateNtuple("T004", "Target0_0_Post");//id=8, poststeppoint
  analysis->CreateNtupleDColumn("PID");       //column 0
  analysis->CreateNtupleDColumn("Z");         //column 1
  analysis->CreateNtupleDColumn("A");         //column 2    
  analysis->CreateNtupleDColumn("Energy");    //column 3
  analysis->CreateNtupleDColumn("Time");      //column 4
  analysis->CreateNtupleDColumn("x");
  analysis->CreateNtupleDColumn("y");
  analysis->CreateNtupleDColumn("z");
  analysis->CreateNtupleDColumn("px");
  analysis->CreateNtupleDColumn("py");
  analysis->CreateNtupleDColumn("pz");
  analysis->CreateNtupleDColumn("trackID");
  analysis->CreateNtupleDColumn("StepNo");
  analysis->CreateNtupleDColumn("Weight");
  analysis->CreateNtupleDColumn("CopyNo");    
  analysis->CreateNtupleDColumn("EnergyDep");    
  analysis->FinishNtuple();
  /////////////////////////////////////////////////////////////////////////////////////////////

  //20170202///////////////////////////////////////////////////
        /////////////////////////////////////////////////////20170202 for initial particle state//in PrimaryGeneratorAction.cc
  //20170206 revived, trajectory.cc
  analysis->CreateNtuple("T000", "Initial Particle state");//id=9
  analysis->CreateNtupleDColumn("PID");       //column 0
  //  analysis->CreateNtupleDColumn("Z");         //column 1
  //  analysis->CreateNtupleDColumn("A");         //column 2    
  analysis->CreateNtupleDColumn("x0");
  analysis->CreateNtupleDColumn("y0");
  analysis->CreateNtupleDColumn("z0");
  analysis->CreateNtupleDColumn("x");
  analysis->CreateNtupleDColumn("y");
  analysis->CreateNtupleDColumn("z");
  analysis->CreateNtupleDColumn("px");
  analysis->CreateNtupleDColumn("py");
  analysis->CreateNtupleDColumn("pz");
  analysis->CreateNtupleDColumn("trackID");
  analysis->CreateNtupleDColumn("PatrackID");
  analysis->CreateNtupleDColumn("traj");
  analysis->FinishNtuple();
  /////////////////////////////////////////////////////////////////////////////////////////////
  analysis->CreateNtuple("T_tot","Energy deposit (MeV) in all the detectors" ); // id=10
  analysis->CreateNtupleDColumn("edep");
  analysis->CreateNtupleDColumn("edep2");//20170512
  analysis->CreateNtupleDColumn("weight");
  //  analysis->CreateNtupleDColumn("copyno");    
  analysis->FinishNtuple();

  //20170223//20170224//////////////////////////////////////////////////////
  for(G4int i = 0 ; i < 16 ; i++){
    std::string text1 = "T_";
    std::string text2 = "Energy deposit (MeV) in the detector CopyNo ";
    text1 += std::to_string(i);
    text2 += std::to_string(i);
    analysis->CreateNtuple(text1, text2); // id=11~
  analysis->CreateNtupleDColumn("edep");
  analysis->CreateNtupleDColumn("edep2");//20170512  
  analysis->CreateNtupleDColumn("weight");
  //  analysis->CreateNtupleDColumn("copyno");    
  analysis->FinishNtuple();
  }
  ////////////////////////////////////////////////////////////////////
  //20170611 Sensitive Detector
  analysis->CreateNtuple("TSD", "SensitiveDetector");//id=27, presteppoint
  analysis->CreateNtupleDColumn("Energy_dep");    //column 0
  analysis->CreateNtupleDColumn("x");         //column 1
  analysis->CreateNtupleDColumn("y");
  analysis->CreateNtupleDColumn("z");
  analysis->CreateNtupleDColumn("trackID");
  analysis->CreateNtupleDColumn("CopyNo");
  analysis->CreateNtupleDColumn("PID");
  analysis->CreateNtupleDColumn("ike"); 
  analysis->CreateNtupleDColumn("L1");
  analysis->CreateNtupleDColumn("L2");
  analysis->CreateNtupleDColumn("light");      
  analysis->FinishNtuple();
  //20170611

  //added 20170613
  analysis->CreateNtuple("TSD_tot", "SensitiveDetector_total_energy");//id=28, presteppoint
  analysis->CreateNtupleDColumn("Energy_tot");    //column 0
  analysis->CreateNtupleDColumn("Energy_tot_gamma");    //column 1
  analysis->CreateNtupleDColumn("Energy_tot_op");    //column 2    
  analysis->CreateNtupleDColumn("light_tot");     //column 3
  analysis->CreateNtupleDColumn("OP");            //column 4      
  analysis->FinishNtuple();
  //added 20170613
  
  analysis->SetNtupleActivation(false);          
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......