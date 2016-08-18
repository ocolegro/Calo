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
// * work  make  anxy representation or  warranty, express or implied, *
// * regarding  this  software system or assume anxy liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * anxy work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// Code developed by:
//  S.Larsson and J. Generowicz.
//
//    *************************************
//    *                                   *
//    *    PurgMagTabulatedField3D.cc     *
//    *                                   *
//    *************************************
//
// $Id: PurgMagTabulatedField3D.cc 84477 2014-10-16 08:44:04Z gcosmo $
//

#include "PurgMagTabulatedField3D.hh"
#include "G4SystemOfUnits.hh"
//#include "G4AutoLock.hh"
#include <fstream>
#include <string>
#include <sstream>
#include <iterator>
#include <algorithm>
namespace{
}

PurgMagTabulatedField3D::PurgMagTabulatedField3D(const char* filename, 
						 G4double zOffset )
  :fZoffset(zOffset),invertX(false),invertY(false),invertZ(false)
{    
 
    G4double lenUnit = centimeter;
    G4double fieldUnit = gauss;
  G4cout << "\n-----------------------------------------------------------"
	 << "\n      Magnetic field"
	 << "\n-----------------------------------------------------------";
    
  G4cout << "\n ---> " "Reading the field grid from " << filename << " ... " << endl; 


  ifstream file( filename ); // Open the file for reading.
  
  if (!file.is_open())
    {
      G4ExceptionDescription ed; ed << "Could not open input file " << filename << G4endl;
      ed << "Could not open input file " << filename << G4endl;
      G4Exception("PurgMagTabulatedField3D::PurgMagTabulatedField3D",
		  "pugmag001",FatalException,ed);
    }

    char buffer[512];



    nz = 0;
    std::string line;
    std::istringstream iss;
    std::vector<string>tokens;
    std::string word;

    while(getline(file, line)){
        ++nz;
        iss.str(line);

    }
    nxy = -1;
    while(iss >> word){
        tokens.push_back(word);
        ++nxy;
    }


    tokens.clear();
    iss.clear();
    word.clear();
    tokens.clear();

	xField.resize(nxy);
	yField.resize(nxy);
	zField.resize(nxy);
	for (unsigned ix=0; ix<nxy; ix++) {
		xField[ix].resize(nxy);
		yField[ix].resize(nxy);
		zField[ix].resize(nxy);
		for (unsigned iy=0; iy<nxy; iy++) {
			xField[ix][iy].resize(nz);
			yField[ix][iy].resize(nz);
			zField[ix][iy].resize(nz);
		}
  }
 G4cout << "Resized arrays successfully." <<G4endl;
 G4cout << "Array size [x]: " << yField.size() << G4endl;
 G4cout << "Array size [y]: " << yField[0].size() << G4endl;
 G4cout << "Array size [z]: " << yField[0][0].size() << G4endl;

  file.clear();
  file.seekg(0,ios::beg);
  G4double xval=0.0, yval=0.0, zval=0.0;
    for (unsigned iz = 0; iz < (nz); iz++){
        if(!tokens.empty()){ 
            tokens.clear();
            iss.clear();
        }
        std::getline(file,line);
        iss.str(line);
        while(iss >> word){
            tokens.push_back(word);
        }
        zval = stod(tokens.at(0)) * cm; // Read in the z-coordinate
        for(unsigned ix=0; ix < (nxy-1); ix++){
            double btemp  = stod(tokens.at(ix+1)) ;
            yField[ix][0][iz] = 3*btemp  * fieldUnit;
            xField[ix][0][iz] = 0.0  * fieldUnit;
            zField[ix][0][iz] = 0.0 * fieldUnit;
            for(unsigned iy = 0; iy < nxy-1; iy++){
                yField[ix][iy][iz] = yField[ix][0][iz];
                xField[ix][iy][iz] = xField[ix][0][iz];
                zField[ix][iy][iz] = zField[ix][0][iz];

            }    
         }
  }
  file.close();

  maxxy = nxy/lenUnit;
  maxz = zval;

}

void PurgMagTabulatedField3D::GetFieldValue(const G4double point[4],
				      G4double *Bfield ) const
{

  G4double lenUnit = cm;// milimeter;
  G4double fieldUnit = gauss;
  G4double x = point[0]/lenUnit;
  G4double y = point[1]/lenUnit;
  G4double z = 1e6*lenUnit;
  if (point[2] < 0 && fabs(point[2]) < fZoffset ){
  	z = (point[2] + fZoffset)/lenUnit ;
  }
  if (point[2] < 0 && fabs(point[2]) >= fZoffset ){
        z = fabs((point[2] + fZoffset)/lenUnit) ;
  } 
// G4cout << z << ", " << point[2] << ", " << fZoffset/lenUnit << ", " << maxz/lenUnit << G4endl;
// G4cout << z << ", " << point[2] << ", " << fZoffset << ", " << maxz << G4endl;

 // Check that the point is within the defined region 
  if ( //fabs(x)>=minxy && fabs(x)<=maxxy &&
       //fabs(y)>=minxy && fabs(y)<=maxxy &&
       fabs(z)<=fabs(maxz)  ) {
	G4int xlow  = floor(fabs(x));
	G4int xhigh = xlow + 1;
	G4int ylow  = floor(fabs(y));
	G4int yhigh = ylow + 1;
	G4int zlow =  (floor(2 *fabs(z)));
	G4int zhigh = zlow + 1;

	double xPercL = 1 - (fabs(x) - floor(x));
	double yPercL = 1 - (fabs(y) - floor(y));
	double zPercL = (fabs(2*z) - floor(2*z));


	if (xlow >= nxy-2){
		xlow = nxy - 2;
		xhigh = nxy - 2;
		xPercL = 1;
	}
	if (ylow >= nxy-2){
		ylow = nxy - 2;
		yhigh = nxy - 2;
		yPercL = 1;
	}
	if (zlow >= nz - 1){
		zlow = nz - 1;
		zhigh = nz - 1;
		zPercL = 1;
	}
	//G4cout << zlow << ", " << zhigh << ", " << xlow << ", " << xhigh << ", " << ylow << ", " << yhigh << G4endl;

	Bfield[0] = 0; Bfield[2] = 0;

    Bfield[1] =
      yField[xlow  ][ylow  ][zlow  ] *		xPercL 	  * 	 yPercL    * 	 zPercL +
      yField[xlow  ][ylow  ][zhigh] * 		xPercL 	  * 	 yPercL    *    (1-zPercL)  +
      yField[xlow  ][yhigh][zlow  ] * 		xPercL 	  *    (1-yPercL)  * 	 zPercL +
      yField[xlow  ][yhigh][zhigh] * 		xPercL 	  *    (1-yPercL)  *    (1-zPercL)  +
      yField[xhigh][ylow  ][zlow  ] *    (1-xPercL)  * 	 yPercL    * 	 zPercL +
      yField[xhigh][ylow  ][zhigh] *    (1-xPercL)  * 	 yPercL    *    (1-zPercL)  +
      yField[xhigh][yhigh][zlow  ] *    (1-xPercL)  *    (1-yPercL)  * 	  zPercL +
      yField[xhigh][yhigh][zhigh] *    (1-xPercL)  *    (1-yPercL)  *    (1-zPercL) ;


  } else {
    Bfield[0] = 0.0;
    Bfield[1] = 0.0;
    Bfield[2] = 0.0;
  }
}

std::vector<std::string> PurgMagTabulatedField3D::Split(const std::string &s, char delim) {
	std::stringstream ss(s);
	std::string item;
	std::vector<std::string> tokens;
	while (std::getline(ss, item, delim)) {
		tokens.push_back(item);
	}
	return tokens;
}




