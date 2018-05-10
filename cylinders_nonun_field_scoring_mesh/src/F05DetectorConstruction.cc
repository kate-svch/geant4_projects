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
//G4double kinetic_energy = 1000.*keV;

G4int number_of_events =  100;    // !!! MUST BE THE SAME VALUE AS IN run1.mac !!!

G4double kinetic_energy = 100.*keV;
G4double from_start_to_the_detector =  4.5*m;     // the parameter - distance from start of the particles to the plane of the detector
G4double distance = from_start_to_the_detector/2.;    // is used in what follows, it's an auxiliary unit!
G4double from_axis_to_cylinder =  3.0*m;
G4double from_electrode_to_edge_of_the_world = 1.0*cm;

G4double det_sizeXY = 3000.*mm; // this is for auxiliary detector and ROOT!!

G4double det_sizeZ = 63.*mm;
G4double det_radius = 63.*mm;
G4double cover_sizeZ = 0.5*mm;
G4double cover_thickness = 0.5*mm;
G4double space_under_detector = 10*cm;
G4double space_near_the_edge = 10*cm;

G4double env_sizeXY = 2*( from_axis_to_cylinder + space_near_the_edge + det_radius + cover_thickness);
G4double env_sizeZ = 2*(distance + det_sizeZ + cover_sizeZ + space_under_detector);

G4double world_coef = 1.;
G4double world_sizeXY = world_coef*env_sizeXY;
G4double world_sizeZ = world_coef*env_sizeZ;   // it's used again later, you should change both blocks simultaneously!!


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

        G4Material* env_mat = nist->FindOrBuildMaterial("G4_Galactic");
  //    G4Material* env_mat = nist->FindOrBuildMaterial("G4_AIR");

      // Option to switch on/off checking of volumes overlaps
      //
      G4bool checkOverlaps = true;

      //
      // World
      //

// Uncomment it, if "World Box size is too small - (0,0,0)"
//      from_start_to_the_detector =  4.5*m;     // the parameter
//      distance = from_start_to_the_detector/2.;    // is used in what follows, it's an auxiliary unit!
//       from_axis_to_cylinder =  3.0*m;
//       from_electrode_to_edge_of_the_world = 10.0*cm;
//                det_sizeZ = 63.*mm;
//       det_radius = 63.*mm;
//          cover_sizeZ = 0.5*mm;
//          cover_thickness = 0.5*mm;
//           space_under_detector = 10*cm;
//       space_near_the_edge = 10*cm;
//      env_sizeXY = 2*( from_axis_to_cylinder + space_near_the_edge + det_radius + cover_thickness);
//    env_sizeZ = 2*(distance + det_sizeZ + cover_sizeZ + space_under_detector);

    world_sizeXY = world_coef*env_sizeXY;
    world_sizeZ = world_coef*env_sizeZ;
      G4Material* world_mat = nist->FindOrBuildMaterial("G4_Galactic");
   //   G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");

      solidWorld =
        new G4Box("World",                       //its name
           0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size

      logicWorld =
        new G4LogicalVolume(solidWorld,          //its solid
                            world_mat,           //its material
                            "World");            //its name

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


       G4Material* cover_mat = nist->FindOrBuildMaterial("G4_Al");


          // Box shape    - DETECTOR, TEMPORARY

//             G4Box* solidShape1 =
//             new G4Box("Shape1",                       //its name
//                0.5*det_sizeXY , 0.5*det_sizeXY, 0.5*det_sizeZ);     //its size


//       G4Box* solidShape1 =
//       new G4Box("Shape1",                       //its name
//          1500*mm, 1500*mm, 500*mm);     //its size

       G4Box* solidShape1 =       new G4Box("Shape1",   0.5*det_sizeXY, 0.5*det_sizeXY, 0.5*det_sizeZ);

           G4LogicalVolume* logicShape1 =
             new G4LogicalVolume(solidShape1,         //its solid
             //                    shape1_mat,          //its material
                                 NaI_Tl,
             //                    NaI,
                                 "Shape1");           //its name

//        G4ThreeVector pos1 = G4ThreeVector(0, 0. , -(distance+ 0.5*det_sizeZ + cover_sizeZ));
                G4ThreeVector pos1 = G4ThreeVector(0, 0. ,  -(distance+ 0.5*det_sizeZ + cover_sizeZ));


           new G4PVPlacement(0,                       //rotation
                             pos1,                    //at position
                             logicShape1,             //its logical volume
                             "Shape1",                //its name
                             logicEnv,                //its mother  volume
                             false,                   //no boolean operation
                             0,                       //copy number
                             checkOverlaps);          //overlaps checking



        // Box shape    - aluminium COVER on the DETECTOR, TEMPORARY

             // G4double cover_sizeY = 0.02*mm;        // is defined earlier

//                G4Box* solidCover =
//                new G4Box("cover",    0.5*det_sizeXY, 0.5*det_sizeXY, 0.5*cover_sizeZ);     //its size



//                G4Box* solidCover =
//                new G4Box("cover",      1500*mm, 1500*mm,  0.5*cover_sizeZ);     //its size


//              G4LogicalVolume* logicCover =
//                new G4LogicalVolume(solidCover,         //its solid
//                                   cover_mat,          //its material
//                                    "cover");           //its name

//             G4ThreeVector posCover = G4ThreeVector(0, 0,  -(distance + 0.5*cover_sizeZ) );



//              new G4PVPlacement(0,                        //rotation
//                                posCover,                    //at position
//                                logicCover,             //its logical volume
//                                "cover",                //its name
//                                logicEnv,                //its mother  volume
//                                false,                   //no boolean operation
//                                0,                       //copy number
//                                checkOverlaps);          //overlaps checking






          G4UserLimits* stepLimit;
          stepLimit = new G4UserLimits(5*mm);

          logicWorld->SetUserLimits(stepLimit);

         fScoringVolume = logicShape1;

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
