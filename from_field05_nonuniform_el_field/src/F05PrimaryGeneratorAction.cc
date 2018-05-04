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
// $Id: F05PrimaryGeneratorAction.cc 69563 2013-05-08 12:30:36Z gcosmo $
//
/// \file field/field05/src/F05PrimaryGeneratorAction.cc
/// \brief Implementation of the F05PrimaryGeneratorAction class
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "F05PrimaryGeneratorAction.hh"
#include "F05DetectorConstruction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"

#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Box.hh"
//#include "G4RunManager.hh"


extern G4double kinetic_energy;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F05PrimaryGeneratorAction::F05PrimaryGeneratorAction(void)
    : G4VUserPrimaryGeneratorAction(),
      fParticleGun(0),
      fEnvelopeBox(0)
{
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);


  // default particle kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle  = particleTable->FindParticle(particleName="e-");
//  G4ParticleDefinition* particle  = particleTable->FindParticle(particleName="alpha");
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,-1.));
  fParticleGun->SetParticleEnergy(kinetic_energy);
  
//  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
// // G4ParticleDefinition* particle = particleTable->FindParticle("mu+");
//  G4ParticleDefinition* particle = particleTable->FindParticle("gamma");

//  fParticleGun->SetParticleDefinition(particle);
}





//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F05PrimaryGeneratorAction::~F05PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F05PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

  //this function is called at the begining of event
  //
//  G4double Pmu = 517.6*MeV;
//  G4double mu_mass = 105.658*MeV;
//  G4double Emu = std::sqrt(Pmu*Pmu + mu_mass*mu_mass);
//  G4double Kmu = Emu - mu_mass;



  G4double envSizeXY = 0;
  G4double envSizeZ = 0;

  if (!fEnvelopeBox)
  {
    G4LogicalVolume* envLV = G4LogicalVolumeStore::GetInstance()->GetVolume("Envelope");
    if ( envLV ) fEnvelopeBox = dynamic_cast<G4Box*>(envLV->GetSolid());
  }

  if ( fEnvelopeBox ) {
    envSizeXY = fEnvelopeBox->GetXHalfLength()*2.;
    envSizeZ = fEnvelopeBox->GetZHalfLength()*2.;
  }
  else  {
    G4ExceptionDescription msg;
    msg << "Envelope volume of box shape not found.\n";
    msg << "Perhaps you have changed geometry.\n";
    msg << "The gun will be place at the center.";
    G4Exception("B1PrimaryGeneratorAction::GeneratePrimaries()",
     "MyCode0002",JustWarning,msg);
  }

  G4double size = 0.25;
  G4double x0 = size * envSizeXY * (G4UniformRand()-0.5);
  G4double y0 = size * envSizeXY * (G4UniformRand()-0.5);

  G4double z0 = 0.5 * envSizeZ;




//  G4double x0 = 1.*m;
//  G4double y0 =  0.00*m;
//  G4double z0 =  0.00*m;


  fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
  //  fParticleGun->SetParticleEnergy(Kmu);
//  fParticleGun->SetParticlePolarization(G4ThreeVector(0.,1.,0.));
//  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,1.,0.));

  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
