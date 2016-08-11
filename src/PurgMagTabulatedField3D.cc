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
  //G4Mutex myPurgMagTabulatedField3DLock = G4MUTEX_INITIALIZER;
}

PurgMagTabulatedField3D::PurgMagTabulatedField3D(const char* filename, 
						 double zOffset ) 
  :fZoffset(zOffset),invertX(false),invertY(false),invertZ(false)
{    
 
  //double lenUnit= meter;
  //double fieldUnit= tesla; 
    double lenUnit = centimeter;
    double fieldUnit = kilogauss;
  G4cout << "\n-----------------------------------------------------------"
	 << "\n      Magnetic field"
	 << "\n-----------------------------------------------------------";
    
  G4cout << "\n ---> " "Reading the field grid from " << filename << " ... " << endl; 

  //
  //This is a thread-local class and we have to avoid that all workers open the 
  //file at the same time
  //G4AutoLock lock(&myPurgMagTabulatedField3DLock);

  ifstream file( filename ); // Open the file for reading.
  
  if (!file.is_open())
    {
      G4ExceptionDescription ed; ed << "Could not open input file " << filename << G4endl;
      ed << "Could not open input file " << filename << G4endl;
      G4Exception("PurgMagTabulatedField3D::PurgMagTabulatedField3D",
		  "pugmag001",FatalException,ed);
    }
  // Ignore first blank line
  //char buffer[256];
  //file.getline(buffer,256);
    char buffer[512];
  // Read table dimensions 
  //file >> nx >> ny >> nz; // Note dodgy order
    // Get number of lines in data file -> nz
    nx = 26;
    ny = nx;
    nz = 0;
    std::string line;
    while(getline(file, line)){
        ++nz;
    }
    std::cout << "Number of lines in data file: " << nz << endl;
/*
 *   G4cout << "  [ Number of values x,y,z: " 
 *      << nx << " " << ny << " " << nz << " ] "
 *      << endl;
 * 
 */
  // Set up storage space for table
  xField.resize(nx*2);
  yField.resize(nx*2);
  zField.resize(nx*2);
  int ix, iy, iz;
  for (ix=0; ix<(2*nx); ix++) {
    xField[ix].resize(ny*2);
    yField[ix].resize(ny*2);
    zField[ix].resize(ny*2);
    for (iy=0; iy<(2*ny); iy++) {
      xField[ix][iy].resize(nz*2);
      yField[ix][iy].resize(nz*2);
      zField[ix][iy].resize(nz*2);
    }
  }
 G4cout << "Resized arrays successfully." <<G4endl;
 G4cout << "Array size [x]: " << yField.size() << G4endl;
 G4cout << "Array size [y]: " << yField[0].size() << G4endl;
 G4cout << "Array size [z]: " << yField[0][0].size() << G4endl;

  // Ignore other header information    
  // The first line whose second character is '0' is considered to
  // be the last line of the header.
  /*
   * do {
   *   file.getline(buffer,256);
   * } while ( buffer[1]!='0');
   */
  // Read in the data
  /*
   * double xval,yval,zval,bx,by,bz;
   * double permeability; // Not used in this example.
   * for (ix=0; ix<nx; ix++) {
   *   for (iy=0; iy<ny; iy++) {
   *     for (iz=0; iz<nz; iz++) {
   *       file >> xval >> yval >> zval >> bx >> by >> bz >> permeability;
   *       if ( ix==0 && iy==0 && iz==0 ) {
   *         minx = xval * lenUnit;
   *         miny = yval * lenUnit;
   *         minz = zval * lenUnit;
   *       }
   *       xField[ix][iy][iz] = bx * fieldUnit;
   *       yField[ix][iy][iz] = by * fieldUnit;
   *       zField[ix][iy][iz] = bz * fieldUnit;
   *     }
   *   }
   * }
   */
  // << nx << G4endl;/ Return to beginning of file
  file.clear();
  file.seekg(0,ios::beg); 
  G4double bval=0.0, xval=0.0, yval=0.0, zval=0.0;
  std::vector<string>tokens; 
  std::istringstream iss;
  std::string word;
    for (iz = 0; iz < (nz-1); iz++){
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
        nx = tokens.size();
        ix =0, iy = 0;
        zval = stod(tokens.at(0)); // Read in the z-coordinate
        //G4cout << "nx = " << nx << G4endl;
        for(ix=0; ix < (nx-2); ix++){
        // Read in all 25 b-field values along the x-axis
            bval = stod(tokens.at(ix+1));
           if ((ix == 0) && (iy == 0) && (iz == 0)){
               minx = -nx*lenUnit;// xval * lenUnit;
               miny = minx; //yval * lenUnit;
               minz = 0;//zval * lenUnit;
           }
           // G4cout << "nx+ix =  "<< nx+ix << " nz+iz = " << nz+iz << G4endl;
            yField[nx+ix][0][nz+iz] = bval * fieldUnit;
           // G4cout << zval << ": yField[" << nx+ix <<"][" << 0 << "][" << nz+iz << "]=" << bval << G4endl;       
            yField[nx-ix][0][nz+iz] = yField[nx+ix][0][nz+iz];
            yField[nx+ix][0][nz-iz] = yField[nx+ix][0][nz+iz];
            yField[nx-ix][0][nz-iz] = yField[nx+ix][0][nz+iz];

            //Apply symmetric considerations along x- axis 
            xField[nx+ix][0][nz+iz] = 0.0  * fieldUnit;
            xField[nx-ix][0][nz+iz] = xField[nx+ix][0][nz+iz];
            xField[nx+ix][0][nz-iz] = xField[nx+ix][0][nz+iz];
            xField[nx-ix][0][nz-iz] = xField[nx+ix][0][nz+iz];
            zField[nx+ix][0][nz+iz] = 0.0 * fieldUnit;
            zField[nx-ix][0][nz+iz] = zField[nx+ix][0][nz+iz];
            zField[nx+ix][0][nz-iz] = zField[nx+ix][0][nz+iz];
            zField[nx-ix][0][nz-iz] = zField[nx+ix][0][nz+iz];
        /* Copy all values along y-axis*/
            for(iy = 1; iy < 2*ny; iy++){
                xField[ix][iy][iz] = xField[ix][0][iz];
                yField[ix][iy][iz] = yField[ix][0][iz];
                zField[ix][iy][iz] = zField[ix][0][iz];
            }    
         }
  }
  file.close();

  //lock.unlock();

  maxx = nx*lenUnit;//xval * lenUnit;
  maxy = maxx;//yval * lenUnit;
  maxz = zval * lenUnit;
  minz = -1 * maxz;
  G4cout << "\n ---> ... done reading " << endl;

  // G4cout << " Read values of field from file " << filename << endl; 
  G4cout << " ---> assumed the order:  x, y, z, Bx, By, Bz "
	 << "\n ---> Min values x,y,z: " 
	 << minx/cm << " " << miny/cm << " " << minz/cm << " cm "
	 << "\n ---> Max values x,y,z: " 
	 << maxx/cm << " " << maxy/cm << " " << maxz/cm << " cm "
	 << "\n ---> The field will be offset by " << zOffset/cm << " cm " << endl;

  // Should really check that the limits are not the wrong way around.
  if (maxx < minx) {swap(maxx,minx); invertX = true;} 
  if (maxy < miny) {swap(maxy,miny); invertY = true;} 
  if (maxz < minz) {swap(maxz,minz); invertZ = true;} 
  G4cout << "\nAfter reordering if neccesary"  
	 << "\n ---> Min values x,y,z: " 
	 << minx/cm << " " << miny/cm << " " << minz/cm << " cm "
	 << " \n ---> Max values x,y,z: " 
	 << maxx/cm << " " << maxy/cm << " " << maxz/cm << " cm ";

  dx = maxx - minx;
  dy = maxy - miny;
  dz = maxz - minz;
  G4cout << "\n ---> Dif values x,y,z (range): " 
	 << dx/cm << " " << dy/cm << " " << dz/cm << " cm in z "
	 << "\n-----------------------------------------------------------" << endl;
}

void PurgMagTabulatedField3D::GetFieldValue(const double point[4],
				      double *Bfield ) const
{
  //std::cout << "In GetFieldValue" << std::endl;
  double x = point[0];
  double y = point[1];
  double z = point[2] + fZoffset;
  
  //std::cout << "B-field vector address: " << Bfield << std::endl;
  //std::cout << "x = " << point[0] << "; y = " << point[1] << "; z = " << point[2] << std::endl;
  //std::cout << "x = " << Bfield[0] << "; y = " << Bfield[1] << "; z = " << Bfield[2] << std::endl;
    
  // Check that the point is within the defined region 
  if ( x>=minx && x<=maxx &&
       y>=miny && y<=maxy && 
       z>=minz && z<=maxz ) {
    
    // Position of given point within region, normalized to the range
    // [0,1]
    double xfraction = (x - minx) / dx;
    double yfraction = (y - miny) / dy; 
    double zfraction = (z - minz) / dz;

    if (invertX) { xfraction = 1 - xfraction;}
    if (invertY) { yfraction = 1 - yfraction;}
    if (invertZ) { zfraction = 1 - zfraction;}

    // Need addresses of these to pass to modf below.
    // modf uses its second argument as an OUTPUT argument.
    double xdindex, ydindex, zdindex;
    
    // Position of the point within the cuboid defined by the
    // nearest surrounding tabulated points
    double xlocal = ( std::modf(xfraction*(nx-1), &xdindex));
    double ylocal = ( std::modf(yfraction*(ny-1), &ydindex));
    double zlocal = ( std::modf(zfraction*(nz-1), &zdindex));
    
    // The indices of the nearest tabulated point whose coordinates
    // are all less than those of the given point
    int xindex = static_cast<int>(xdindex);
    int yindex = static_cast<int>(ydindex);
    int zindex = static_cast<int>(zdindex);
    

#ifdef DEBUG_INTERPOLATING_FIELD
    G4cout << "Local x,y,z: " << xlocal << " " << ylocal << " " << zlocal << endl;
    G4cout << "Index x,y,z: " << xindex << " " << yindex << " " << zindex << endl;
    double valx0z0, mulx0z0, valx1z0, mulx1z0;
    double valx0z1, mulx0z1, valx1z1, mulx1z1;
    valx0z0= table[xindex  ][0][zindex];  mulx0z0=  (1-xlocal) * (1-zlocal);
    valx1z0= table[xindex+1][0][zindex];  mulx1z0=   xlocal    * (1-zlocal);
    valx0z1= table[xindex  ][0][zindex+1]; mulx0z1= (1-xlocal) * zlocal;
    valx1z1= table[xindex+1][0][zindex+1]; mulx1z1=  xlocal    * zlocal;
#endif

        // Full 3-dimensional version
    Bfield[0] =
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
      zField[xindex+1][yindex+1][zindex+1] *    xlocal  *    ylocal  *    zlocal ;
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

