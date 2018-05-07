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
// $Id: F05Field.cc 75672 2013-11-05 08:47:41Z gcosmo $
//
/// \file field/field05/src/F05Field.cc
/// \brief Implementation of the F05Field class
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "F05Field.hh"

#include "G4SystemOfUnits.hh"

#include "F05DetectorConstruction.hh"

extern G4double from_electrode_to_edge_of_the_world;
extern G4double world_sizeZ;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F05Field::F05Field() : G4ElectroMagneticField()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F05Field::~F05Field()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F05Field::GetFieldValue( const G4double Point[4], G4double* Bfield ) const
{
  // Point[0],Point[1],Point[2] are x-, y-, z-cordinates, Point[3] is time

//  const G4double Bz = 0.24*tesla;
//  const G4double Er = 2.113987E+6*volt/m;


  const G4double k_coef = 9.*std::pow(10,9)*volt*m/coulomb;
 // const G4double q_charge =10.*std::pow(10,-4)*coulomb;
    const G4double q_charge =1.*std::pow(10,-4)*coulomb;

  G4double Ex,Ey, Ez;
 
  G4double D_distance = 0.5*world_sizeZ + from_electrode_to_edge_of_the_world;
  G4double lateral_R = std::sqrt(std::pow(Point[0],2) + std::pow(Point[1],2));
  G4double r_real = std::sqrt(std::pow(Point[0],2) + std::pow(Point[1],2) + std::pow(D_distance - Point[2], 2));
  G4double r_mirror = std::sqrt(std::pow(Point[0],2) + std::pow(Point[1],2) + std::pow(D_distance + Point[2], 2));

  G4double r_real_minus_three = 1/std::pow(r_real, 3);
  G4double r_mirror_minus_three = 1/std::pow(r_mirror, 3);


  // just one charge, 'point-like', in (0,0, D_distance)
//    Ex = -k_coef*q_charge*Point[0]*(r_real_minus_three);
//    Ey = -k_coef*q_charge*Point[1]*(r_real_minus_three );
//    Ez = k_coef*q_charge*(r_real_minus_three*(D_distance - Point[2]));


  //  one charge in (0,0, D_distance) and it's MIRROR in (0,0, -D_distance)
  Ex = -k_coef*q_charge*Point[0]*(r_real_minus_three - r_mirror_minus_three);
  Ey = -k_coef*q_charge*Point[1]*(r_real_minus_three - r_mirror_minus_three);
  Ez = k_coef*q_charge*(r_real_minus_three*(D_distance - Point[2]) + r_mirror_minus_three*(D_distance + Point[2]));


  // test field
//  Ex = k_coef*q_charge*Point[0]/std::pow(lateral_R, 3);
//  Ey = k_coef*q_charge*Point[1]/std::pow(lateral_R, 3);
//  Ez = k_coef*q_charge/(std::pow(D_distance/2,2) + std::pow(D_distance/2,2))/1000;

//  Ex = 0;
//  Ey = -1.*tesla*Point[2]/D_distance;
//  Ez = 0;


//  if (posR>0){
//     cos_theta = Point[0]/(G4double)posR;
//     sin_theta = Point[1]/(G4double)posR;
//     Ex = -1*Er*cos_theta;//apply radial electric field
//     Ey = -1*Er*sin_theta;
//  }else{
//     Ex=0;
//     Ey=0;
//  }
  
  Bfield[0]=0;
  Bfield[1]=0;
  Bfield[2]=0;

  Bfield[3]=Ex;
  Bfield[4]=Ey;
  Bfield[5]=Ez;

  return;
}
