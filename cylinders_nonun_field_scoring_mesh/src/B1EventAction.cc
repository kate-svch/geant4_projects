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
// $Id: B1EventAction.cc 93886 2015-11-03 08:28:26Z gcosmo $
//
/// \file B1EventAction.cc
/// \brief Implementation of the B1EventAction class

#include "B1EventAction.hh"
#include "B1RunAction.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"

// for ROOT
#include "from_B4Analysis.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"
#include <iomanip>

#include <sstream>


#include <iostream>
#include <fstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//B1EventAction::B1EventAction(B1RunAction* runAction)
//: G4UserEventAction(),
//  fRunAction(runAction),
//  fEdep(0.)
//{}

// This is from B4 and seems to be unneeded
B1EventAction::B1EventAction()
 : G4UserEventAction(),
   fEdep(0.)
{}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1EventAction::~B1EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1EventAction::BeginOfEventAction(const G4Event*)
{    
  fEdep = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


//G4double B1EventAction::get_particle_energy(const G4Event* event )
//{

//    G4PrimaryVertex* primaryVertex = event->GetPrimaryVertex();
//    G4PrimaryParticle* primaryParticle = primaryVertex->GetPrimary();
//    G4double ke = primaryParticle->GetKineticEnergy();


//    std::ofstream myfile;
//    myfile.open ("example.txt");
//    myfile << "Kinetic energy is: " << ke;
//    myfile.close();


//    G4cout << G4endl  << G4endl
//   << "cout:--------------------Kinetic energy is:  " << ke << "----------------"    << G4endl << G4endl;

//          G4cerr << G4endl  << G4endl
//         << "cerr:--------------------Kinetic energy is:  " << ke << "----------------"    << G4endl << G4endl;


//   return ke;

//}

void B1EventAction::EndOfEventAction(const G4Event* event)
{
     // accumulate statistics in run action
    fRunAction->AddEdep(fEdep);

  // Accumulate statistics

  // get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // fill histograms
  if (fEdep != 0) {
  analysisManager->FillH1(0, fEdep);
    }

      G4PrimaryVertex* primaryVertex = event->GetPrimaryVertex();
      G4PrimaryParticle* primaryParticle = primaryVertex->GetPrimary();
      G4double ke = primaryParticle->GetKineticEnergy();

      G4cout << G4endl  << G4endl
       << "cout:--------------------Kinetic energy is:  " << ke << "----------------"    << G4endl << G4endl;



            G4cerr << G4endl  << G4endl
           << "cerr:--------------------Kinetic energy is:  " << ke << "----------------"    << G4endl << G4endl;


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
