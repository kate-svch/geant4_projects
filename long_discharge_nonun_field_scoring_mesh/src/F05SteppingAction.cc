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
// $Id: F05SteppingAction.cc 68021 2013-03-13 13:36:07Z gcosmo $
//
/// \file field/field05/src/F05SteppingAction.cc
/// \brief Implementation of the F05SteppingAction class
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "F05SteppingAction.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4ProcessTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"

#include "B1EventAction.hh"
#include "F05DetectorConstruction.hh"

#include "G4UnitsTable.hh"  // for G4BestUnit


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F05SteppingAction::F05SteppingAction(void)
{;}


//F05SteppingAction::F05SteppingAction(B1EventAction* eventAction)
//: G4UserSteppingAction(),
//  fEventAction(eventAction),
//  fScoringVolume(0)
//{}



//F05SteppingAction::F05SteppingAction(
//                      const F05DetectorConstruction* detectorConstruction,
//                      B1EventAction* eventAction)
//  : G4UserSteppingAction(),
//    fDetConstruction(detectorConstruction),
//    fEventAction(eventAction)
//{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F05SteppingAction::~F05SteppingAction()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F05SteppingAction::UserSteppingAction(const G4Step* step)
{

     if (!fScoringVolume) {
         const F05DetectorConstruction* detectorConstruction
           = static_cast<const F05DetectorConstruction*>
             (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
         fScoringVolume = detectorConstruction->GetScoringVolume();
       }

       // get volume of the current step
       G4LogicalVolume* volume
         = step->GetPreStepPoint()->GetTouchableHandle()
           ->GetVolume()->GetLogicalVolume();

       // check if we are in scoring volume
       if (volume != fScoringVolume) return;

       // collect energy deposited in this step
       G4double edepStep = step->GetTotalEnergyDeposit();
       fEventAction->AddEdep(edepStep);
}
