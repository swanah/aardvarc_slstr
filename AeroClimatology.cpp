/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   AeroClimatology.cpp
 * Author: akheckel
 * 
 * Created on 29. März 2017, 12:14
 */

#include <netcdf>
#include <stdexcept>
#include "AeroClimatology.hpp"
#include "Interpolation.hpp"
#include "defs.hpp"

using namespace netCDF;

AeroClimatology::AeroClimatology(const std::string& climName, const int& month) {

    // for now assuming global 1x1° climatology
    nx = 360;
    ny = 180;
    lat_start = 90.;
    lon_start = -180.;
    dlat = -1;
    dlon =  1;
    
    std::vector<size_t> start;
    start.push_back(month-1);
    start.push_back(0);
    start.push_back(0);
    std::vector<size_t> count;
    count.push_back(1);
    count.push_back(ny);
    count.push_back(nx);
    
    NcFile ncF(climName, NcFile::read);
    
    NcVar v = ncF.getVar("AOD550");
    aod = new float[nx * ny];
    v.getVar(start, count, aod);
    v = ncF.getVar("fine_of_total_fraction");
    fineOfTotal = new float[nx * ny];
    v.getVar(start, count, fineOfTotal);
    v = ncF.getVar("lessAbs_of_fine_fraction");
    weakOfFine = new float[nx * ny];
    v.getVar(start, count, weakOfFine);
    v = ncF.getVar("dust_of_coarse_fraction");
    dustOfCoarse = new float[nx * ny];
    v.getVar(start, count, dustOfCoarse);
    
    //ncF.close();
}

AeroClimatology::~AeroClimatology() {
    delete [] aod;
    aod = NULL;
    delete [] dustOfCoarse;
    dustOfCoarse = NULL;
    delete [] weakOfFine;
    weakOfFine = NULL;
    delete [] fineOfTotal;
    fineOfTotal = NULL;
}

void AeroClimatology::getMixPercentages(const GeoPos& gp, float *mixPercentages, float *mixFrac, float *climAod) {
    int ilat = (int) ((gp.lat - lat_start) / dlat);
    int ilon = (int) ((gp.lon - lon_start) / dlon);

    int index;
    float fot, wof, doc;

    if (ilat == ny) ilat--;
    if (ilon == nx) ilon--;
    if (ilat < 0 || ilat > ny || ilon < 0 || ilon > nx) {
        throw std::runtime_error("Error lat/lon index for climatology out of bounds!");
    }

    index = ilat * nx + ilon;
    fot = fineOfTotal[index];
    fot = (fot>0) ? fot : 0.05;
    fot = (fot<1) ? fot : 0.95;
    wof = weakOfFine[index];
    //wof = wof + 0.5 * ( 1 - wof );
    //wof = 0.95;
    
    doc = dustOfCoarse[index];
    //doc = 0.5 * clim->dustOfCoarse[index];
    
    *climAod = aod[index];
    mixFrac[0] = fot; 
    mixFrac[1] = wof; 
    mixFrac[2] = doc;

    mixPercentages[0] = ((1 - fot) * doc)*100; // Dust
    mixPercentages[1] = ((1 - fot) * (1 - doc))*100; // Sea Salt
    mixPercentages[2] = (fot * (1 - wof))*100; // Strong Abs
    mixPercentages[3] = (fot * wof)*100; // Weak Abs
}

void AeroClimatology::getMixPercentagesInt(const GeoPos& gp, float *mixPercentages, float *mixFrac, float *climAod) {
    double x2 = (gp.lon - lon_start) / dlon;
    double y2 = (gp.lat - lat_start) / dlat;
    float fot, wof, doc;
    
    fot = interpol_2d_img(fineOfTotal, x2, y2, nx, ny, false);    
    fot = (fot>0) ? fot : 0.05;
    fot = (fot<1) ? fot : 0.95;

    wof = interpol_2d_img(weakOfFine, x2, y2, nx, ny, false);
    //wof = wof + 0.1 * ( 1 - wof );
    wof = 0.35 * pow((wof + 0.51), 2) + 0.2;
    //wof = 0.95;

    //doc = fineOfTotal[((int)y2) * nx + (int)x2];
    doc = interpol_2d_img(dustOfCoarse, x2, y2, nx, ny, false);
    doc = 2 * sqrt(doc + 1) - 2;
    //doc = 0.8 * doc;

    *climAod = interpol_2d_img(aod, x2, y2, nx, ny, false);

    mixFrac[0] = fot; 
    mixFrac[1] = wof; 
    mixFrac[2] = doc;

    mixPercentages[0] = ((1 - fot) * doc)*100; // Dust
    mixPercentages[1] = ((1 - fot) * (1 - doc))*100; // Sea Salt
    mixPercentages[2] = (fot * (1 - wof))*100; // Strong Abs
    mixPercentages[3] = (fot * wof)*100; // Weak Abs
}
