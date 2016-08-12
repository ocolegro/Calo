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
    	G4cout << "iz = " << iz << ", nz = " << nz << G4endl;
        if(!tokens.empty()){ 
            tokens.clear();
            iss.clear();
        }
        std::getline(file,line);
        iss.str(line);
        while(iss >> word){
            tokens.push_back(word);
        }
        //G4cout << G4endl << "Vector size = " << tokens.size() << G4endl;
        nxy = tokens.size();
        zval = stod(tokens.at(0)) * cm; // Read in the z-coordinate
        for(unsigned ix=0; ix < (nxy-2); ix++){
        	//G4cout << "ix = " << ix << ", nxy = " << nxy << G4endl;
            double btemp  = stod(tokens.at(ix+1)) ;
            G4cout << "The field in the y direction at iz = " << iz << ", ix = " << ix << " = " << stod(tokens.at(ix+1)) << G4endl;
            G4cout << "The field in the zval direction at iz = " << iz << "is " << zval/cm << G4endl;

            yField[ix][0][iz] = btemp * fieldUnit;
            G4cout << "The field in the y direction at iz = " << iz << ", ix = " << ix << " = " << yField[ix][0][iz]/gauss << G4endl;

            xField[ix][0][iz] = 0.0  * fieldUnit;
            zField[ix][0][iz] = 0.0 * fieldUnit;
        /* Copy all values along y-axis*/
            for(unsigned iy = 1; iy < nxy-2; iy++){
            	//G4cout << "iy = " << iy << ", nxy = " << nxy << G4endl;
                yField[ix][iy+1][iz] = yField[ix][0][iz];
                xField[ix][iy+1][iz] = xField[ix][0][iz];
                zField[ix][iy+1][iz] = zField[ix][0][iz];

            }    
         }
  }
  file.close();

  maxxy = nxy;
  maxz = zval;
  minz = -1 * maxz;
  minxy = -nxy*lenUnit;
  dxy = maxxy - minxy;
  dz = maxz - minz;

}

void PurgMagTabulatedField3D::GetFieldValue(const G4double point[4],
				      G4double *Bfield ) const
{

  G4double lenUnit = centimeter;
  G4double fieldUnit = gauss;
  G4double x = point[0];
  G4double y = point[1];
  G4double ztrue = point[2] ;//+ fZoffset)/lenUnit ;
  G4double z = ztrue;//150/lenUnit - ztrue;
  G4cout << "The maxz is " << maxz/cm << G4endl;
  G4cout << "The ztrue is " << ztrue/cm << " the z we attempt to read is " << z <<  G4endl;
  bool printField = false;
  if (Bfield[0] == 999){
	  printField = true;
  }
  // Check that the point is within the defined region 
  if ( x>=minxy && x<=maxxy &&
       y>=minxy && y<=maxxy &&
       z>=minz && z<=maxz ) {
    
    // Position of given point within region, normalized to the range
    // [0,1]
    G4double xfraction = (x - minxy) / dxy;
    G4double yfraction = (y - minxy) / dxy;
    G4double zfraction = (z - minz) / dz;

    if (invertX) { xfraction = 1 - xfraction;}
    if (invertY) { yfraction = 1 - yfraction;}
    if (invertZ) { zfraction = 1 - zfraction;}

    // Need addresses of these to pass to modf below.
    // modf uses its second argument as an OUTPUT argument.
    G4double xdindex, ydindex, zdindex;
    
    // Position of the point within the cuboid defined by the
    // nearest surrounding tabulated points
    G4double xlocal = ( std::modf(xfraction*(nxy-1), &xdindex));
    G4double ylocal = ( std::modf(yfraction*(nxy-1), &ydindex));
    G4double zlocal = ( std::modf(zfraction*(nz-1), &zdindex));
    
    // The indices of the nearest tabulated point whose coordinates
    // are all less than those of the given point
    int xindex = static_cast<int>(xdindex);
    int yindex = static_cast<int>(ydindex);
    int zindex = static_cast<int>(zdindex);
    

#ifdef DEBUG_INTERPOLATING_FIELD
    G4cout << "Local x,y,z: " << xlocal << " " << ylocal << " " << zlocal << endl;
    G4cout << "Index x,y,z: " << xindex << " " << yindex << " " << zindex << endl;
    G4double valx0z0, mulx0z0, valx1z0, mulx1z0;
    G4double valx0z1, mulx0z1, valx1z1, mulx1z1;
    valx0z0= table[xindex  ][0][zindex];  mulx0z0=  (1-xlocal) * (1-zlocal);
    valx1z0= table[xindex+1][0][zindex];  mulx1z0=   xlocal    * (1-zlocal);
    valx0z1= table[xindex  ][0][zindex+1]; mulx0z1= (1-xlocal) * zlocal;
    valx1z1= table[xindex+1][0][zindex+1]; mulx1z1=  xlocal    * zlocal;
#endif

        // Full 3-dimensional version
  /*  Bfield[0] =
      xField[xindex  ][yindex  ][zindex  ] * (1-xlocal) * (1-ylocal) * (1-zlocal) +
      xField[xindex  ][yindex  ][zindex+1] * (1-xlocal) * (1-ylocal) *    zlocal  +
      xField[xindex  ][yindex+1][zindex  ] * (1-xlocal) *    ylocal  * (1-zlocal) +
      xField[xindex  ][yindex+1][zindex+1] * (1-xlocal) *    ylocal  *    zlocal  +
      xField[xindex+1][yindex  ][zindex  ] *    xlocal  * (1-ylocal) * (1-zlocal) +
      xField[xindex+1][yindex  ][zindex+1] *    xlocal  * (1-ylocal) *    zlocal  +
      xField[xindex+1][yindex+1][zindex  ] *    xlocal  *    ylocal  * (1-zlocal) +
      xField[xindex+1][yindex+1][zindex+1] *    xlocal  *    ylocal  *    zlocal ;
    Bfield[1] =
      yField[xindex  ][yindex  ][zindex  ] * (1-xlocal) * (1-ylocal) * (1-zlocal) +
      yField[xindex  ][yindex  ][zindex+1] * (1-xlocal) * (1-ylocal) *    zlocal  +
      yField[xindex  ][yindex+1][zindex  ] * (1-xlocal) *    ylocal  * (1-zlocal) +
      yField[xindex  ][yindex+1][zindex+1] * (1-xlocal) *    ylocal  *    zlocal  +
      yField[xindex+1][yindex  ][zindex  ] *    xlocal  * (1-ylocal) * (1-zlocal) +
      yField[xindex+1][yindex  ][zindex+1] *    xlocal  * (1-ylocal) *    zlocal  +
      yField[xindex+1][yindex+1][zindex  ] *    xlocal  *    ylocal  * (1-zlocal) +
      yField[xindex+1][yindex+1][zindex+1] *    xlocal  *    ylocal  *    zlocal ;
    Bfield[2] =
      zField[xindex  ][yindex  ][zindex  ] * (1-xlocal) * (1-ylocal) * (1-zlocal) +
      zField[xindex  ][yindex  ][zindex+1] * (1-xlocal) * (1-ylocal) *    zlocal  +
      zField[xindex  ][yindex+1][zindex  ] * (1-xlocal) *    ylocal  * (1-zlocal) +
      zField[xindex  ][yindex+1][zindex+1] * (1-xlocal) *    ylocal  *    zlocal  +
      zField[xindex+1][yindex  ][zindex  ] *    xlocal  * (1-ylocal) * (1-zlocal) +
      zField[xindex+1][yindex  ][zindex+1] *    xlocal  * (1-ylocal) *    zlocal  +
      zField[xindex+1][yindex+1][zindex  ] *    xlocal  *    ylocal  * (1-zlocal) +
      zField[xindex+1][yindex+1][zindex+1] *    xlocal  *    ylocal  *    zlocal ;*/
    Bfield[0] = 0; Bfield[1] = yField[xindex][yindex][zindex] ; Bfield[2] = 0;
    if (printField){
    	G4cout << "The x,y,z that we are reading in is: " << x/cm << ", " << y/cm << ", " << ztrue/cm << G4endl;
    	G4cout << "The recalled filed, before passing was :  "  << " (" << Bfield[0]/gauss << ", " << Bfield[1]/gauss << ", " << Bfield[2]/gauss << " )" << G4endl;
    }
  } else {
    Bfield[0] = 0.0;
    Bfield[1] = 0.0;
    Bfield[2] = 0.0;
  }
  //std::cout << "bx = " << Bfield[0] << "; by = " << Bfield[1] << "; bz = " << Bfield[2] << std::endl;
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

