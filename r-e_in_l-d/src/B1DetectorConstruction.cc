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

// The main parameters
G4double kinetic_energy = 1000*keV;
G4double distance = 200*cm;
//G4double alpha =  0;
G4double alpha =  M_PI/6;
G4double Electric_field_y = 10.0*kilovolt/cm;   // directed in +y-direction!! (so, it eccelerates elelctrons, which are negative)

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

  G4double det_sizeY = 8.*cm;
  G4double det_sizeXZ = 100.*cm;
  G4double cover_sizeY = 0.02*mm;
  G4double space_under_detector = 10*cm;


  G4double env_sizeXZ = 2.0*distance + det_sizeY + cover_sizeY;
  G4double env_sizeY = distance + det_sizeY + cover_sizeY+space_under_detector;
  //G4Material* env_mat = nist->FindOrBuildMaterial("G4_Galactic");
  G4Material* env_mat = nist->FindOrBuildMaterial("G4_AIR");
   
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //     
  // World
  //
  G4double world_sizeXZ = 1.2*env_sizeXZ;
  G4double world_sizeY  = 1.2*env_sizeY;
  //G4Material* world_mat = nist->FindOrBuildMaterial("G4_Galactic");
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
  
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       0.5*world_sizeXZ, 0.5*world_sizeY, 0.5*world_sizeXZ);     //its size
      
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
        0.5*env_sizeXZ, 0.5*env_sizeY, 0.5*env_sizeXZ); //its size
      
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


 // G4double alpha = M_PI/4.;
 // G4double alpha = 0;



// // G4ThreeVector pos1 = G4ThreeVector(1*m, 1*m, 0);
  G4ThreeVector pos1 = G4ThreeVector((distance+det_sizeY*0.5+cover_sizeY)*sin(alpha*rad), -(distance+det_sizeY*0.5+cover_sizeY)*cos(alpha*rad) + env_sizeY/2.0, 0);

 //G4Material* shape1_mat = nist->FindOrBuildMaterial("G4_Al");


  // Box shape    - DETECTOR
  
    G4Box* solidShape1 =
    new G4Box("Shape1",                       //its name
       0.5*det_sizeXZ , 0.5*det_sizeY, 0.5*det_sizeXZ);     //its size
  
                      
  G4LogicalVolume* logicShape1 =
    new G4LogicalVolume(solidShape1,         //its solid
    //                    shape1_mat,          //its material
                        NaI_Tl,
    //                    NaI,
                        "Shape1");           //its name


  G4RotationMatrix* zRot = new G4RotationMatrix;
  zRot->rotateZ(-alpha);                         // Rotates X and Y axes only

               
  new G4PVPlacement(zRot,                       //rotation
                    pos1,                    //at position
                    logicShape1,             //its logical volume
                    "Shape1",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking


  // Box shape    - aluminium COVER on the DETECTOR

 // G4double cover_sizeY = 0.02*mm;        // is defined earlier

    G4Box* solidCover =
    new G4Box("cover",                       //its name
       0.5*det_sizeXZ, 0.5*cover_sizeY, 0.5*det_sizeXZ);     //its size

    G4Material* cover_mat = nist->FindOrBuildMaterial("G4_Al");

  G4LogicalVolume* logicCover =
    new G4LogicalVolume(solidCover,         //its solid
                       cover_mat,          //its material
                        "cover");           //its name

   //   G4ThreeVector posCover = G4ThreeVector((distance+(det_sizeY+cover_sizeY)*0.5)*sin(alpha*rad), -(distance+(det_sizeY+cover_sizeY)*0.5)*cos(alpha*rad) + env_sizeY/2.0, 0);
       G4ThreeVector posCover = G4ThreeVector((distance+(cover_sizeY)*0.5)*sin(alpha*rad), -(distance+(cover_sizeY)*0.5)*cos(alpha*rad) + env_sizeY/2.0, 0);

  new G4PVPlacement(zRot,                        //rotation
                    posCover,                    //at position
                    logicCover,             //its logical volume
                    "Cover",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking




 //the following is unneeded due to the first two arguments of G4PVPlacement?

//    G4RotationMatrix yRot45deg;   // Rotates X and Y axes only

//    yRot45deg.rotateZ(M_PI/4.*rad);
//    G4ThreeVector  translation(0, 0, 0);
//    G4Box  rotatedCover("rotatedCover", &solidEnv,&solidCover,&yRot45deg,translation);
    // The new coordinate system of the cylinder is translated so that
    // its centre is at +50 on the original Z axis, and it is rotated
    // with its X axis halfway between the original X and Z axes.



                    
  // Set Shape1 as scoring volume (detector)
  //
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
  G4UniformElectricField* entireVolumeField = new G4UniformElectricField(G4ThreeVector(0.0, Electric_field_y, 0.0));
  // G4UniformElectricField* entireVolumeField = new G4UniformElectricField(G4ThreeVector(Electric_field_y, 0.0, 0.0));

   //G4UniformElectricField* entireVolumeField = new G4UniformElectricField(G4ThreeVector(0.0, 0.0, 0.0));
   G4FieldManager* entireFieldMgrVolume = new   G4FieldManager(entireVolumeField);
   G4bool includeDaugthers = true; // или true

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
