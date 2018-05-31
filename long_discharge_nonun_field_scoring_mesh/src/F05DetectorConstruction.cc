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
// $Id: F05DetectorConstruction.cc 101905 2016-12-07 11:34:39Z gunter $
//
/// \file field/field05/src/F05DetectorConstruction.cc
/// \brief Implementation of the F05DetectorConstruction class
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "F05DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4UserLimits.hh"
#include "G4SystemOfUnits.hh"
#include "F05Field.hh"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"

//#include "G4RepleteEofM.hh"
#include "G4EqEMFieldWithSpin.hh"

#include "G4ClassicalRK4.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4ChordFinder.hh"
#include "G4PropagatorInField.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// The main parameters


// The main parameters
G4double kinetic_energy = 1000*keV;
G4double distance = 400*cm;
G4double alpha =  0;
//G4double alpha =  M_PI/6;
//G4double Electric_field_y = 10.0*kilovolt/cm;   // directed in +y-direction!! (so, it eccelerates elelctrons, which are negative)
G4int number_of_events = 100;

G4double det_sizeZ = 8.*cm;
G4double det_sizeXY = 100.*cm;
G4double cover_sizeZ = 0.02*mm;
G4double space_under_detector = 10*cm;


F05DetectorConstruction::F05DetectorConstruction()
 : fVacuum(0), world_sizeXY(0), world_sizeZ(0),
   solidWorld(0), logicWorld(0), physWorld(0),
    fScoringVolume(0)
{
  // materials
  DefineMaterials();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F05DetectorConstruction::~F05DetectorConstruction()
{
  if (fField) delete fField;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F05DetectorConstruction::DefineMaterials()
{
  G4NistManager* nistMan = G4NistManager::Instance();

  fVacuum = nistMan->FindOrBuildMaterial("G4_Galactic");

  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* F05DetectorConstruction::Construct()
{

        // Get nist material manager
      G4NistManager* nist = G4NistManager::Instance();

      // Envelope, detector and cover parameters



      G4double env_sizeXY = 2.0*distance + det_sizeZ + cover_sizeZ;
      G4double env_sizeZ = distance + det_sizeZ + cover_sizeZ + space_under_detector;
   //   G4Material* env_mat = nist->FindOrBuildMaterial("G4_Galactic");
      G4Material* env_mat = nist->FindOrBuildMaterial("G4_AIR");

      // Option to switch on/off checking of volumes overlaps
      //
      G4bool checkOverlaps = true;

      //
      // World
      //
       world_sizeZ = 1.2*env_sizeZ;
       world_sizeXY  = 1.2*env_sizeXY;
      //G4Material* world_mat = nist->FindOrBuildMaterial("G4_Galactic");
      G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");

  //    G4Box*   // it's already declared as a member of 'F05DetectorConstruction'
      solidWorld =
        new G4Box("World",                       //its name
           0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size

  //    G4LogicalVolume*  // it's already declared as a member of 'F05DetectorConstruction'
              logicWorld =
        new G4LogicalVolume(solidWorld,          //its solid
                            world_mat,           //its material
                            "World");            //its name

    //  G4VPhysicalVolume*  // it's already declared as a member of 'F05DetectorConstruction'
        physWorld =
        new G4PVPlacement(0,                     //no rotation
                          G4ThreeVector(),       //at (0,0,0)
                          logicWorld,            //its logical volume
                          "World",               //its name
                          0,                     //its mother  volume
                          false,                 //no boolean operation
                          0,                     //copy number
                          checkOverlaps);        //overlaps checking

      //
      // Envelope
      //
      G4Box* solidEnv =
        new G4Box("Envelope",                    //its name
            0.5*env_sizeXY, 0.5*env_sizeXY, 0.5*env_sizeZ); //its size

      G4LogicalVolume* logicEnv =
        new G4LogicalVolume(solidEnv,            //its solid
                            env_mat,             //its material
                            "Envelope");         //its name

      new G4PVPlacement(0,                       //no rotation
                       // G4ThreeVector(0, -0.5*env_sizeY, 0),         //at (0,0,0)   // THIS DEFINES THE INTERSECTION OF THE AXIS, REALLY! and the centre of the Envelope volume
                       G4ThreeVector(0, 0, 0),
                        logicEnv,                //its logical volume
                        "Envelope",              //its name
                        logicWorld,              //its mother  volume
                        false,                   //no boolean operation
                        0,                       //copy number
                        checkOverlaps);          //overlaps checking

      //
      // Shape 1 - Detector
      //
    //  G4Material* shape1_mat = nist->FindOrBuildMaterial("G4_Galactic");
    //  G4Material* shape1_mat = nist->FindOrBuildMaterial("G4_AIR");
    //    G4Material* shape1_mat = nist->FindOrBuildMaterial("G4_Pb");


      // ------------------------------------------------------------------------
       // Elements
       // ------------------------------------------------------------------------
        G4double A, Z;
     //    G4Element* elO  = new G4Element("Oxygen",  "O",  Z=8.,  A= 15.9994*g/mole);
         G4Element* elNa = new G4Element("Sodium",  "Na", Z=11., A= 22.989768*g/mole);
         G4Element* elI  = new G4Element("Iodine",  "I",  Z=53., A= 126.90447*g/mole);
      //   G4Element* elMg = new G4Element("Magnesium",  "Mg", Z=12., A= 24.305*g/mole);
         G4Element* elTl  = new G4Element("Thallium",  "Tl",  Z=81., A= 204.38*g/mole);


      // NaI crystal
      G4int natoms, nel;
      G4double density = 3.67 *g/cm3;
      G4Material* NaI = new G4Material("NaI", density, nel= 2);
      NaI-> AddElement(elNa, natoms=1);
      NaI-> AddElement(elI,  natoms=1);


        // Let us make the NaI_Tl alloy

        //  we made the "HUGE MOLECULE" - should we use the mass fractions and define the compound instead?
        // at least, now it works

       density = 3.67 *g/cm3;
       G4Material* NaI_Tl = new G4Material("NaI_Tl", density, nel= 3);
       NaI_Tl-> AddElement(elNa, natoms=500);
       NaI_Tl-> AddElement(elI,  natoms=500);
       NaI_Tl-> AddElement(elTl,  natoms=1);


    // // G4ThreeVector pos1 = G4ThreeVector(1*m, 1*m, 0);
      G4ThreeVector pos1 = G4ThreeVector((distance+det_sizeZ*0.5+cover_sizeZ)*sin(alpha*rad), 0, -(distance+det_sizeZ*0.5+cover_sizeZ)*cos(alpha*rad) + env_sizeZ/2.0);

     //G4Material* shape1_mat = nist->FindOrBuildMaterial("G4_Al");

          G4RotationMatrix* zRot = new G4RotationMatrix;
          zRot->rotateZ(-alpha);                         // Rotates X and Y axes only



//      // Box shape    - DETECTOR

//        G4Box* solidShape1 =
//        new G4Box("Shape1",                       //its name
//           0.5*det_sizeXY , 0.5*det_sizeXY, 0.5*det_sizeZ);     //its size


//      G4LogicalVolume* logicShape1 =
//        new G4LogicalVolume(solidShape1,         //its solid
//        //                    shape1_mat,          //its material
//                            NaI_Tl,
//        //                    NaI,
//                            "Shape1");           //its name



//      new G4PVPlacement(zRot,                       //rotation
//                        pos1,                    //at position
//                        logicShape1,             //its logical volume
//                        "Shape1",                //its name
//                        logicEnv,                //its mother  volume
//                        false,                   //no boolean operation
//                        0,                       //copy number
//                        checkOverlaps);          //overlaps checking






      // Box shape    - aluminium COVER on the DETECTOR

     // G4double cover_sizeY = 0.02*mm;        // is defined earlier

        G4Box* solidCover =
        new G4Box("cover",                       //its name
           0.5*det_sizeXY, 0.5*det_sizeXY, 0.5*cover_sizeZ);     //its size

        G4Material* cover_mat = nist->FindOrBuildMaterial("G4_Al");

      G4LogicalVolume* logicCover =
        new G4LogicalVolume(solidCover,         //its solid
                           cover_mat,          //its material
                            "cover");           //its name

       //   G4ThreeVector posCover = G4ThreeVector((distance+(det_sizeY+cover_sizeY)*0.5)*sin(alpha*rad), -(distance+(det_sizeY+cover_sizeY)*0.5)*cos(alpha*rad) + env_sizeY/2.0, 0);
           G4ThreeVector posCover = G4ThreeVector((distance+(cover_sizeZ)*0.5)*sin(alpha*rad), 0, -(distance+(cover_sizeZ)*0.5)*cos(alpha*rad) + env_sizeZ/2.0);

      new G4PVPlacement(zRot,                        //rotation
                        posCover,                    //at position
                        logicCover,             //its logical volume
                        "Cover",                //its name
                        logicEnv,                //its mother  volume
                        false,                   //no boolean operation
                        0,                       //copy number
                        checkOverlaps);          //overlaps checking







          G4UserLimits* stepLimit;
          stepLimit = new G4UserLimits(5*mm);

          logicWorld->SetUserLimits(stepLimit);

    //     fScoringVolume = logic_detectorTube;

          return physWorld;



}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal F05Field* F05DetectorConstruction::fField = 0;

void F05DetectorConstruction::ConstructSDandField()

{
  if (!fField) {

     fField = new F05Field();

//     G4RepleteEofM* equation = new G4RepleteEofM(fField);
     G4EqEMFieldWithSpin* equation = new G4EqEMFieldWithSpin(fField);
//     equation->SetBField();
//     equation->SetEField();
//     equation->SetSpin();

     G4FieldManager* fieldManager
      = G4TransportationManager::GetTransportationManager()->GetFieldManager();
     fieldManager->SetDetectorField(fField);

     G4MagIntegratorStepper* stepper = new G4ClassicalRK4(equation,12);

     G4double minStep           = 0.01*mm;

     G4ChordFinder* chordFinder =
                    new G4ChordFinder((G4MagneticField*)fField,minStep,stepper);

     // Set accuracy parameters
     G4double deltaChord        = 3.0*mm;
     chordFinder->SetDeltaChord( deltaChord );

     G4double deltaOneStep      = 0.01*mm;
     fieldManager->SetAccuraciesWithDeltaOneStep(deltaOneStep);

     G4double deltaIntersection = 0.1*mm;
     fieldManager->SetDeltaIntersection(deltaIntersection);

     G4TransportationManager* transportManager =
                           G4TransportationManager::GetTransportationManager();

     G4PropagatorInField* fieldPropagator =
                                      transportManager->GetPropagatorInField();

     G4double epsMin            = 2.5e-7*mm;
     G4double epsMax            = 0.05*mm;

     fieldPropagator->SetMinimumEpsilonStep(epsMin);
     fieldPropagator->SetMaximumEpsilonStep(epsMax);

     fieldManager->SetChordFinder(chordFinder);

     G4bool allLocal = true;
     logicWorld->SetFieldManager(fieldManager, allLocal);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
