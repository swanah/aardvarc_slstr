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
#include "AeroClimatology.hpp"

using namespace netCDF;

AeroClimatology::AeroClimatology(const std::string climName, const int month) {

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
    
    ncF.close();
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

