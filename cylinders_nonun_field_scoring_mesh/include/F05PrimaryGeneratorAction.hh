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
// $Id: F05PrimaryGeneratorAction.hh 69563 2013-05-08 12:30:36Z gcosmo $
//
/// \file field/field05/include/F05PrimaryGeneratorAction.hh
/// \brief Definition of the F05PrimaryGeneratorAction class
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef F05PrimaryGeneratorAction_h
#define F05PrimaryGeneratorAction_h 1

#include "globals.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"

class G4Event;
class G4Box;
class G4ParticleGun;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class F05PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:

  F05PrimaryGeneratorAction(void);    
  virtual ~F05PrimaryGeneratorAction();

  virtual void GeneratePrimaries(G4Event*);

  // method to access particle gun
  const G4ParticleGun* GetParticleGun() const { return fParticleGun; }

private:
  G4ParticleGun*  fParticleGun; // pointer a to G4 gun class
  G4Box* fEnvelopeBox;
};


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
