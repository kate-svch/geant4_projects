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
//#include "Randomize.hh"
//#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1EventAction::B1EventAction(B1RunAction* runAction)
: G4UserEventAction(),
  fRunAction(runAction),
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

// an old version of void B1EventAction

//void B1EventAction::EndOfEventAction(const G4Event*)
//{
//  // accumulate statistics in run action
//  fRunAction->AddEdep(fEdep);
//}


// a new version of void B1EventAction

//void B1EventAction::EndOfEventAction(const G4Event* event)
void B1EventAction::EndOfEventAction(const G4Event* )
{
  // Accumulate statistics

  // get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // fill histograms
  if (fEdep != 0) {
  analysisManager->FillH1(0, fEdep);
	}
//  analysisManager->FillH1(1, fEnergyGap);

  // fill ntuple
//  analysisManager->FillNtupleDColumn(0, fEdep);
// analysisManager->FillNtupleDColumn(1, fEnergyGap);
//  analysisManager->AddNtupleRow();

 // it's just an output??
  // Print per event (modulo n)
  //
//  auto eventID = event->GetEventID();
//  auto printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
//  if ( ( printModulo > 0 ) && ( eventID % printModulo == 0 ) ) {
//    G4cout << "---> End of event: " << eventID << G4endl;

//    G4cout
//       << "   Absorber: total energy: " << std::setw(7)
//                                        << G4BestUnit(fEnergyAbs,"Energy")
//       << "       total track length: " << std::setw(7)
//                                        << G4BestUnit(fTrackLAbs,"Length")
//       << G4endl
//       << "        Gap: total energy: " << std::setw(7)
//                                        << G4BestUnit(fEnergyGap,"Energy")
//       << "       total track length: " << std::setw(7)
//                                        << G4BestUnit(fTrackLGap,"Length")
//       << G4endl;
//  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
