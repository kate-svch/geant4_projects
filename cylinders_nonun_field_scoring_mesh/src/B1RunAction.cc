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
// $Id: B1RunAction.cc 99560 2016-09-27 07:03:29Z gcosmo $
//
/// \file B1RunAction.cc
/// \brief Implementation of the B1RunAction class

#include "B1RunAction.hh"
#include "F05PrimaryGeneratorAction.hh"
#include "F05DetectorConstruction.hh"
// #include "B1Run.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "from_B4Analysis.hh"

#include <sstream>


extern G4double q_charge;
//extern G4double distance;
//extern G4double alpha;
extern G4double kinetic_energy;
extern G4double from_start_to_the_detector;
extern G4int number_of_events;
extern G4double cover_thickness;
extern G4double det_sizeZ;
extern G4double det_sizeXY;

//extern G4double from_electrode_to_edge_of_the_world;
//extern G4double from_axis_to_cylinder ;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1RunAction::B1RunAction()
: G4UserRunAction(),
  fEdep(0.),
  fEdep2(0.)
{ 
  // add new units for dose
  // 
  const G4double milligray = 1.e-3*gray;
  const G4double microgray = 1.e-6*gray;
  const G4double nanogray  = 1.e-9*gray;  
  const G4double picogray  = 1.e-12*gray;
   
  new G4UnitDefinition("milligray", "milliGy" , "Dose", milligray);
  new G4UnitDefinition("microgray", "microGy" , "Dose", microgray);
  new G4UnitDefinition("nanogray" , "nanoGy"  , "Dose", nanogray);
  new G4UnitDefinition("picogray" , "picoGy"  , "Dose", picogray); 

  // Register accumulable to the accumulable manager
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->RegisterAccumulable(fEdep);
  accumulableManager->RegisterAccumulable(fEdep2); 


  // set printing event number per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1);

  // Create analysis manager
  // The choice of analysis technology is done via selectin of a namespace in B4Analysis.hh
  auto analysisManager = G4AnalysisManager::Instance();
  G4cout << "Using " << analysisManager->GetType() << G4endl;

  // Create directories
  analysisManager->SetHistoDirectoryName("histograms");
  analysisManager->SetNtupleDirectoryName("ntuple");
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetNtupleMerging(true);

  // Creating histograms
//  analysisManager->CreateH1("fEder","Edep in detector", 100, 0., 5.0*MeV);


//  double num = 2.25149;
//  std::stringstream ss(stringstream::in | stringstream::out);
//  ss << setprecision(5) << num << endl;
//  ss << setprecision(4) << num << endl;


  G4String this_string_is_charge;          // string which will contain the result
  std::ostringstream convert_charge;   // stream used for the conversion
  convert_charge << q_charge / (coulomb);      // insert the textual representation of 'Number' in the characters in the stream
  this_string_is_charge = convert_charge.str();

  G4String this_string_is_distance;          // string which will contain the result
  std::ostringstream convert_distance;   // stream used for the conversion
  convert_distance << from_start_to_the_detector / cm;      // insert the textual representation of 'Number' in the characters in the stream
  this_string_is_distance = convert_distance.str();


//  G4String this_string_is_angle;          // string which will contain the result
//  std::ostringstream convert_angle;   // stream used for the conversion
//  convert_angle << alpha/M_PI*180.0;     // insert the textual representation of 'Number' in the characters in the stream
//  this_string_is_angle = convert_angle.str();

  G4String this_string_is_energy;          // string which will contain the result
  std::ostringstream convert_energy;   // stream used for the conversion
  convert_energy << kinetic_energy/keV;     // insert the textual representation of 'Number' in the characters in the stream
  this_string_is_energy = convert_energy.str();

  G4String this_string_is_number_of_initial_electrons;          // string which will contain the result
  std::ostringstream convert_number_of_init_e;   // stream used for the conversion
  convert_number_of_init_e << number_of_events;     // insert the textual representation of 'Number' in the characters in the stream
  this_string_is_number_of_initial_electrons = convert_number_of_init_e.str();

//  G4String this_string_is_distance_from_start_to_layer;          // string which will contain the result
//  std::ostringstream convert_distance_fstl;   // stream used for the conversion
//  convert_distance_fstl << distance_from_start_to_layer/mm;     // insert the textual representation of 'Number' in the characters in the stream
//  this_string_is_distance_from_start_to_layer = convert_distance_fstl.str();

  G4String this_string_is_cover_sizeZ;          // string which will contain the result
  std::ostringstream convert_cover_sizeZ;    // stream used for the conversion
  convert_cover_sizeZ << cover_thickness/ mm;     // insert the textual representation of 'Number' in the characters in the stream
  this_string_is_cover_sizeZ = convert_cover_sizeZ.str();

  G4String this_string_is_det_sizeZ;          // string which will contain the result
  std::ostringstream convert_det_sizeZ;   // stream used for the conversion
  convert_det_sizeZ << det_sizeZ/ cm;     // insert the textual representation of 'Number' in the characters in the stream
  this_string_is_det_sizeZ = convert_det_sizeZ.str();


  G4String this_string_is_det_sizeXY;          // string which will contain the result
  std::ostringstream convert_det_sizeXY ;   // stream used for the conversion
  convert_det_sizeXY << det_sizeXY/ cm;     // insert the textual representation of 'Number' in the characters in the stream
  this_string_is_det_sizeXY = convert_det_sizeXY.str();


  G4String about_distance_fstl ("Distance from start to the detector = ");
  G4String unit_of_distance_fstl (" mm");
  G4String physical_value_and_charge ("   Edep in detector, q = -");
  G4String unit_of_Energy (" keV/cm");
  G4String about_distance (", distance = ");
  G4String unit_of_distance (" cm");
  G4String about_angle (", angle = ");
  G4String unit_of_angle (" degrees");
  G4String about_energy (", energy = ");
  G4String unit_of_energy (" keV");
  G4String about_number_of_initial_e (", number_of_init_e = ");
  G4String about_cover_sizeZ (", thickness of the cover = ");
  G4String unit_of_cover_sizeZ (" mm");
  G4String about_det_sizes ("Detector: ");
  G4String unit_of_det_sizes (" cm");
  G4String unit_of_charge (" coulomb");


G4String title_preliminary =
        about_distance_fstl +  +  this_string_is_distance + unit_of_distance
       // + this_string_is_distance_from_start_to_layer + unit_of_distance_fstl
        +    about_det_sizes  + this_string_is_det_sizeXY + '*' + this_string_is_det_sizeXY
        +    '*' + this_string_is_det_sizeZ + unit_of_det_sizes +
         about_cover_sizeZ + this_string_is_cover_sizeZ + unit_of_cover_sizeZ +
          physical_value_and_charge + this_string_is_charge + unit_of_charge +
    //    about_distance +  this_string_is_distance + unit_of_distance
       // + about_angle + this_string_is_angle + unit_of_angle
        + about_energy + this_string_is_energy + unit_of_energy +
        about_number_of_initial_e + this_string_is_number_of_initial_electrons;

  const G4String& title = title_preliminary;

//  gr->SetTitle("#splitline{a simple}{graph}");

    analysisManager->CreateH1("fEder", title,
                             100, 0., 5.0*MeV,  "MeV");


// analysisManager->CreateH1("fEder", "Edep in detector",
//                          100, 0., 1.0*MeV,  "MeV");

    // The meaning  of arguments is as follows:
//  G4int CreateH1(const G4String& name, const G4String& title,
//  G4int  nbins, G4double xmin, G4double xmax,
//  const  G4String& unitName = "none",
//  const  G4String& fcnName = "none",
//  const   G4String& binSchemeName = "linear");


//  analysisManager->CreateH1("Egap","Edep in gap", 100, 0., 100*MeV);

  // Creating ntuple
  //
  analysisManager->CreateNtuple("from_B4", "just Edep in detector");
  analysisManager->CreateNtupleDColumn("Edep");
//  analysisManager->CreateNtupleDColumn("Egap");
  analysisManager->FinishNtuple();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1RunAction::~B1RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::BeginOfRunAction(const G4Run*)
{ 
  // inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);

  // reset accumulables to their initial values
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Reset();

   // for ROOT:

   // Get analysis manager
   auto analysisManager = G4AnalysisManager::Instance();

   G4String fileName = "test";

       G4bool fileOpen = analysisManager->OpenFile(fileName);


       if (! fileOpen) {
           G4cerr  << G4endl  << G4endl << "\n----------------------> HistoManager::Book(): cannot open "
                  << analysisManager->GetFileName() << G4endl  << G4endl<< G4endl;
            return;
          }


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::EndOfRunAction(const G4Run* run)
{
  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;

  // Merge accumulables 
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Merge();

  // Compute dose = total energy deposit in a run and its variance
  //
  G4double edep  = fEdep.GetValue();
  G4double edep2 = fEdep2.GetValue();
  
  G4double rms = edep2 - edep*edep/nofEvents;
  if (rms > 0.) rms = std::sqrt(rms); else rms = 0.;  

  const F05DetectorConstruction* detectorConstruction
   = static_cast<const F05DetectorConstruction*>
     (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  G4double mass = detectorConstruction->GetScoringVolume()->GetMass();
  G4double dose = edep/mass;
  G4double rmsDose = rms/mass;

  // Run conditions
  //  note: There is no primary generator action object for "master"
  //        run manager for multi-threaded mode.
  const F05PrimaryGeneratorAction* generatorAction
   = static_cast<const F05PrimaryGeneratorAction*>
     (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  G4String runCondition;
  if (generatorAction)
  {
    const G4ParticleGun* particleGun = generatorAction->GetParticleGun();
    runCondition += particleGun->GetParticleDefinition()->GetParticleName();
    runCondition += " of ";
    G4double particleEnergy = particleGun->GetParticleEnergy();
    runCondition += G4BestUnit(particleEnergy,"Energy");
  }
        
  // Print
  //  
  if (IsMaster()) {
    G4cout
     << G4endl
     << "--------------------End of Global Run-----------------------";
  }
  else {
    G4cout
     << G4endl
     << "--------------------End of Local Run------------------------";
  }
  
  G4cout
     << G4endl
     << " The run consists of " << nofEvents << " "<< runCondition
     << G4endl
     << " Cumulated dose per run, in scoring volume : " 
     << G4BestUnit(dose,"Dose") << " rms = " << G4BestUnit(rmsDose,"Dose")
     << G4endl
     << "------------------------------------------------------------"
     << G4endl
     << G4endl;

  // print histogram statistics
  //
  auto analysisManager = G4AnalysisManager::Instance();
  G4cout << " Edep : mean = "
     << G4BestUnit(analysisManager->GetH1(0)->mean(), "Energy")
     << " rms = "
     << G4BestUnit(analysisManager->GetH1(0)->rms(),  "Energy") << G4endl;

  // save histograms & ntuple
  //
  analysisManager->Write();
  analysisManager->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::AddEdep(G4double edep)
{
  fEdep  += edep;
  fEdep2 += edep*edep;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

