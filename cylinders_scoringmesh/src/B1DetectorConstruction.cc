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
// $Id: B1DetectorConstruction.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file B1DetectorConstruction.cc
/// \brief Implementation of the B1DetectorConstruction class

#include "B1DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"


// Let's describe the electric field
#include "G4ElectricField.hh"
#include "G4EqMagElectricField.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4FieldManager.hh"
#include "G4ChordFinder.hh"

#include "G4UniformElectricField.hh"
#include "G4ClassicalRK4.hh"
#include "G4TransportationManager.hh"
#include "G4MagIntegratorDriver.hh"

// The main parameters
G4double kinetic_energy = 100*keV;
G4double from_start_to_the_detector =  4.5*m;     // the parameter
G4double distance = from_start_to_the_detector/2.;    // is used in what follows, it's an auxiliary unit!
G4double from_axis_to_cylinder =  3.0*m;
//G4double alpha =  0;
//G4double alpha =  M_PI/6;
G4double Electric_field_z = 7.0*kilovolt/cm;   // directed in +y-direction!! (so, it eccelerates elelctrons, which are negative)

// It'll be needed for the electric field.

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::B1DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::~B1DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//G4double get_the_field_value() {
//    return Electric_field_y;
//}


G4VPhysicalVolume* B1DetectorConstruction::Construct()
{  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  // Envelope, detector and cover parameters

  G4double det_sizeZ = 63.*mm;
  G4double det_radius = 63.*mm;
  G4double cover_sizeZ = 0.5*mm;
  G4double cover_thickness = 0.5*mm;
  G4double space_under_detector = 10*cm;
  G4double space_near_the_edge = 10*cm;



  // additional detector, for temporary use - from Kostinsky geometry
  G4double det_sizeXY = 100.*cm;


  G4double env_sizeXY = 2*( from_axis_to_cylinder + space_near_the_edge + det_radius + cover_thickness);
  G4double env_sizeZ = 2*(distance + det_sizeZ + cover_sizeZ + space_under_detector);

//  G4double env_sizeXY = 4*m;
//  G4double env_sizeZ =4*m;

  //G4Material* env_mat = nist->FindOrBuildMaterial("G4_Galactic");
  G4Material* env_mat = nist->FindOrBuildMaterial("G4_AIR");
   
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //     
  // World
  //
  G4double world_sizeXY = 1.2*env_sizeXY;
  G4double world_sizeZ = 1.2*env_sizeZ;
  //G4Material* world_mat = nist->FindOrBuildMaterial("G4_Galactic");
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
  
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name
                                   
  G4VPhysicalVolume* physWorld = 
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

 G4ThreeVector pos_NaI = G4ThreeVector(from_axis_to_cylinder, 0. , -(distance+ det_sizeZ + cover_sizeZ));
 // G4ThreeVector pos_NaI = G4ThreeVector(from_axis_to_cylinder, -(distance+ det_sizeY*0.5 + cover_sizeY) + env_sizeY/2.0, 0);

 //G4Material* shape1_mat = nist->FindOrBuildMaterial("G4_Al");


  // Cylinder shape    - DETECTOR
  
   G4Tubs* solid_detectorTube = new G4Tubs("NaI_detector", 0.*m, det_radius, det_sizeZ, 0.*deg, 360.*deg);

  G4LogicalVolume* logic_detectorTube =
    new G4LogicalVolume(solid_detectorTube,         //its solid
                        NaI_Tl,
                        "NaI_detector");           //its name
  new G4PVPlacement(0,                       //no rotation
                    pos_NaI,                    //at position
                    logic_detectorTube,             //its logical volume
                    "NaI_detector",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking



  // Cylinder shape    - aluminium COVER on the DETECTOR

 // it could be used as "cylindric walls around the Nai-detector"
  //   G4Tubs* solid_coverTube = new G4Tubs("Al_cover", det_radius, det_radius + cover_thickness, det_sizeZ, 0.*deg, 360.*deg);

  // this is cover UNDER the detector
  G4Tubs* solid_coverTube = new G4Tubs("Al_cover", 0, det_radius, cover_sizeZ, 0.*deg, 360.*deg);


    G4Material* cover_mat = nist->FindOrBuildMaterial("G4_Al");

  G4LogicalVolume* logic_coverTube =
    new G4LogicalVolume(solid_coverTube,         //its solid
                       cover_mat,          //its material
                        "Al_cover");           //its name

   //  G4ThreeVector pos_cover = G4ThreeVector(from_axis_to_cylinder, -(distance+det_sizeY*0.5+cover_sizeY) + env_sizeY/2.0 + det_sizeY/2., 0);
//   G4ThreeVector pos_cover = G4ThreeVector(from_axis_to_cylinder, 0,  -(distance + 0.5*cover_sizeZ)  );
      G4ThreeVector pos_cover = G4ThreeVector(from_axis_to_cylinder, 0,  -(distance + 0.5*cover_sizeZ) );

  new G4PVPlacement(0,                        //rotation
                    pos_cover,                    //at position
                    logic_coverTube,             //its logical volume
                    "Al_cover",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking




      // Box shape    - DETECTOR, TEMPORARY

         G4Box* solidShape1 =
         new G4Box("Shape1",                       //its name
            0.5*det_sizeXY , 0.5*det_sizeXY, 0.5*det_sizeZ);     //its size


       G4LogicalVolume* logicShape1 =
         new G4LogicalVolume(solidShape1,         //its solid
         //                    shape1_mat,          //its material
                             NaI_Tl,
         //                    NaI,
                             "Shape1");           //its name

    G4ThreeVector pos1 = G4ThreeVector(0, 0. , -(distance+ det_sizeZ + cover_sizeZ));


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

         G4Box* solidCover =
         new G4Box("cover",                       //its name
            0.5*det_sizeXY, 0.5*det_sizeXY, 0.5*cover_sizeZ);     //its size

       G4LogicalVolume* logicCover =
         new G4LogicalVolume(solidCover,         //its solid
                            cover_mat,          //its material
                             "cover");           //its name

      G4ThreeVector posCover = G4ThreeVector(0, 0,  -(distance + 0.5*cover_sizeZ) );


       new G4PVPlacement(0,                        //rotation
                         posCover,                    //at position
                         logicCover,             //its logical volume
                         "Cover",                //its name
                         logicEnv,                //its mother  volume
                         false,                   //no boolean operation
                         0,                       //copy number
                         checkOverlaps);          //overlaps checking


                    
  // Set Shape1 as scoring volume (detector)
  //
 // fScoringVolume = logic_detectorTube;
      fScoringVolume = logicShape1;
  //fScoringVolume = logicEnv;




//  G4ElectricField*  fEMfield;
  G4EqMagElectricField* fEquation;
  G4MagIntegratorStepper* fStepper;
//  G4FieldManager* fFieldManager;
  G4double fMinStep ;
  G4ChordFinder*  fChordFinder ;

  G4MagInt_Driver* fIntgrDriver;  // this one was absent in MSU-description


  // Now we describe the electric field

 //  G4UniformElectricField* entireVolumeField = new G4UniformElectricField(G4ThreeVector(0.0, 500.0*kilovolt/m, 0.0));
  G4UniformElectricField* entireVolumeField = new G4UniformElectricField(G4ThreeVector(0.0, 0.0, Electric_field_z));


   G4FieldManager* entireFieldMgrVolume = new   G4FieldManager(entireVolumeField);
   G4bool includeDaugthers = true;

   // уравнение движения в поле
   fEquation = new G4EqMagElectricField(entireVolumeField);

   G4int nvar = 8; // число переменных
   fStepper = new G4ClassicalRK4( fEquation, nvar );
   entireFieldMgrVolume -> SetDetectorField(entireVolumeField);
   fMinStep = 0.010*mm ; // минимальный шаг 10 микрон
   fIntgrDriver = new G4MagInt_Driver(fMinStep, fStepper, fStepper->GetNumberOfVariables() );
   fChordFinder = new G4ChordFinder(fIntgrDriver);
   entireFieldMgrVolume -> SetChordFinder( fChordFinder );

//   logicEnv ->  SetFieldManager(entireFieldMgrVolume,   includeDaugthers);
   logicWorld ->  SetFieldManager(entireFieldMgrVolume,   includeDaugthers);   // anyway

   // electric field. We had just described it.

               
  //
  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
