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
// $Id: field05.cc 78053 2013-12-03 08:25:13Z gcosmo $
//
/// \file field/field05/field05.cc
/// \brief Main program of the field/field05 example
//
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef WIN32
#include <unistd.h>
#endif

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "F05SteppingVerbose.hh"
#include "G4RunManager.hh"
#endif

#include "G4ScoringManager.hh"
#include "F05PhysicsList.hh"
#include "F05DetectorConstruction.hh"

#include "F05ActionInitialization.hh"

#include "G4UImanager.hh"

#include "Randomize.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

#include "QBBC.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
    G4Colour orange(1, 0.6, 0.0, 1);
    G4Colour::AddToMap("orange", orange);

  // Choose the Random engine
  //
  G4Random::setTheEngine(new CLHEP::RanecuEngine);

  G4int myseed = 1234;
  if (argc  > 2) myseed = atoi(argv[argc-1]);

  // Construct the default run manager
  //
#ifdef G4MULTITHREADED
  G4MTRunManager * runManager = new G4MTRunManager;
#else
  G4VSteppingVerbose::SetInstance(new F05SteppingVerbose);
  G4RunManager * runManager = new G4RunManager;
#endif

  // Activate UI-command base scorer
  G4ScoringManager * scManager = G4ScoringManager::GetScoringManager();
  scManager->SetVerboseLevel(1);

  G4Random::setTheSeed(myseed);

  // Set mandatory initialization classes
  //
  // Detector construction
  runManager->SetUserInitialization(new F05DetectorConstruction());


  // Physics list
//  runManager->SetUserInitialization(new F05PhysicsList());

  G4VModularPhysicsList* physicsList = new QBBC;
  physicsList->SetVerboseLevel(1);
  runManager->SetUserInitialization(physicsList);

  // User action initialization
  runManager->SetUserInitialization(new F05ActionInitialization());

  // Initialize G4 kernel
  //
  runManager->Initialize();



#ifdef G4VIS_USE
  // Initialize visualization
  //
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();
#endif
  
  // Get the pointer to the User Interface manager
  //
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  
  if (argc!=1)   // batch mode
  {
     G4String command = "/control/execute ";
     G4String fileName = argv[1];
     UImanager->ApplyCommand(command+fileName);    
  }
  else
  {  // interactive mode : define UI session
#ifdef G4UI_USE
     G4UIExecutive* ui = new G4UIExecutive(argc, argv);
#ifdef G4VIS_USE
     UImanager->ApplyCommand("/control/execute vis.mac");
#else
     UImanager->ApplyCommand("/control/execute field05.in");
#endif
     if (ui->IsGUI())
        UImanager->ApplyCommand("/control/execute gui.mac");
     ui->SessionStart();
     delete ui;
#endif
  }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  //                 owned and deleted by the run manager, so they should not
  //                 be deleted in the main() program !

#ifdef G4VIS_USE
  delete visManager;
#endif                
  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
