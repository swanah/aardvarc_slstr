/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   S3CciL2Writer.cpp
 * Author: akheckel
 * 
 * Created on 22. November 2017, 08:30
 */

#include <string>
#include <netcdf>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>
#include "S3CciL2Writer.hpp"

using namespace netCDF;


//
// public
//

S3CciL2Writer::S3CciL2Writer(const S3NcdfData& s3NcdfData, const int& numberOfPixels) : s3Data(s3NcdfData), nPix(numberOfPixels) {
    
    initNcdfVars();
    
    int iPix = 0;
    int idx;
    for (int j=0; j<s3Data.sceneOutputHeight; j++) {
        idx = j * s3Data.sceneOutputWidth;
        for (int i=0; i<s3Data.sceneOutputWidth; i++) {
            idx += i;
            if (s3Data.s3AodImgs[0].isValidValue(s3Data.s3AodImgs[0].img[idx])) {
                pixelNumber.vec.push_back(iPix);
                
                for (int iWvl=0; iWvl<N_SLSTR_BANDS; iWvl++){
                    aod[iWvl].vec.push_back(s3Data.s3AodImgs[iWvl].img[idx]);
                    aodUnc[iWvl].vec.push_back(s3Data.s3UncImgs[iWvl].img[idx]);
                }
                
                lat.vec.push_back(s3Data.s3LatImgs[0].img[idx]);
                lon.vec.push_back(s3Data.s3LonImgs[0].img[idx]);
                /**
                for (int iCrn=0; iCrn<N_SLSTR_BANDS; iCrn++){
                    latCrn.vec.push_back(s3Data.s3LatCrnImgs[iCrn].img[idx]);
                    lonCrn.vec.push_back(s3Data.s3LonCrnImgs[iCrn].img[idx]);
                }
                /**/
                
                sza.vec.push_back(s3Data.s3SzaImgs[0].img[idx]);
                vza.vec.push_back(s3Data.s3VzaImgs[0].img[idx]);
                raz.vec.push_back(s3Data.s3RazImgs[0].img[idx]);
                
                iPix++;
            }
        }
    }
    nPix = pixelNumber.vec.size();
}

void S3CciL2Writer::write() {
    printf("writing \n%s\n", s3Data.pars.aodOutName.c_str());
    
    NcFile ncOut(s3Data.pars.aodOutName, NcFile::replace);
    
    NcDim pixDim = ncOut.addDim(pixelNumber.varName, nPix);
    std::vector<NcDim> dimVec;
    dimVec.push_back(pixDim);
    
    addWriteVariable(ncOut, pixelNumber, dimVec);
    for (int iWvl=0; iWvl<N_SLSTR_BANDS; iWvl++){
        addWriteVariable(ncOut, aod[iWvl], dimVec);
    }

    for (int iWvl=0; iWvl<N_SLSTR_BANDS; iWvl++){
        addWriteVariable(ncOut, aodUnc[iWvl], dimVec);
    }

    addWriteVariable(ncOut, lat, dimVec);
    addWriteVariable(ncOut, lon, dimVec);
    for (int iCrn=0; iCrn<N_SLSTR_BANDS; iCrn++){
        //addWriteVariable(ncOut, latCrn[iCrn], dimVec);
        //addWriteVariable(ncOut, lonCrn[iCrn], dimVec);
    }
    
    addWriteVariable(ncOut, sza, dimVec);
    addWriteVariable(ncOut, vza, dimVec);
    addWriteVariable(ncOut, raz, dimVec);
    
    
    
/*    
    NcVar nPixVar = ncOut.addVar(pixelNumber.varName, ncInt, dimVec);
    nPixVar.setCompression(true, true, 5);
    if (!pixelNumber.longName.empty()) {
        nPixVar.putAtt("long_name", pixelNumber.longName);
    }
    if (!pixelNumber.stdName.empty()) {
        nPixVar.putAtt("standard_name", pixelNumber.stdName);
    }
    if (!pixelNumber.unit.empty()) {
        nPixVar.putAtt("units", pixelNumber.unit);
    }
    int valRange[] = {1, nPix};
    nPixVar.putAtt("valid_range", ncInt, 2, &valRange);
    int* data = pixelNumber.vec.data();
    nPixVar.putVar(pixelNumber.vec.data());
    
    NcVar aod550Var = ncOut.addVar(aod550.varName, ncFloat, dimVec);
    aod550Var.setCompression(true, true, 5);
    if (!aod550.longName.empty()) {
        aod550Var.putAtt("long_name", aod550.longName);
    }
    if (!aod550.stdName.empty()) {
        aod550Var.putAtt("standard_name", aod550.stdName);
    }
    if (!aod550.unit.empty()) {
        aod550Var.putAtt("units", aod550.unit);
    }
    float* dataF = aod550.vec.data();
    aod550Var.putVar(aod550.vec.data());
*/    
    writeGlobalAttributes(ncOut);
    //ncOut.close();
}


S3CciL2Writer::~S3CciL2Writer() {
}

//
// private
//

S3CciL2Writer::initNcdfVars() {
    
    pixelNumber.varName = "pixel_number";
    pixelNumber.longName = "pixel index";
    pixelNumber.ncTypeId = ncInt.getId();
    pixelNumber.unit = "1";
    pixelNumber.vec.reserve(nPix);
    
    std::string wlStr[] = {"550", "670", "865", "1600", "2250"};
    for (int i=0; i<N_SLSTR_BANDS; i++){
        aod[i].varName = "AOD" + wlStr[i];
        aod[i].longName = "aerosol optical thickness at " + wlStr[i] + " nm";
        aod[i].stdName = "atmosphere_optical_thickness_due_to_ambient_aerosol";
        aod[i].ncTypeId = ncFloat.getId();
        aod[i].unit = "1";
        aod[i].setFillValue(-999.);
        aod[i].setValidRange(0., 4.);
        aod[i].vec.reserve(nPix);

        aodUnc[i].varName = "AOD" + wlStr[i] + "_uncertainty";
        aodUnc[i].longName = "uncertainty on AOD at " + wlStr[i] + " nm";
        aodUnc[i].stdName = "atmosphere_optical_thickness_due_to_ambient_aerosol";
        aodUnc[i].ncTypeId = ncFloat.getId();
        aodUnc[i].unit = "1";
        aodUnc[i].setValidRange(0., 10.);
        aodUnc[i].setFillValue(-999.);
        aodUnc[i].vec.reserve(nPix);
    }
    /**
    aod[0].varName = "AOD550";
    aod[0].longName = "aerosol optical thickness at 550 nm";
    aod[0].stdName = "atmosphere_optical_thickness_due_to_ambient_aerosol";
    aod[0].ncTypeId = ncFloat.getId();
    aod[0].unit = "1";
    aod[0].setFillValue(-999.);
    aod[0].setValidRange(0., 4.);
    aod[0].vec.reserve(nPix);
    
    aod[1].varName = "AOD670";
    aod[1].longName = "aerosol optical thickness at 659 nm";
    aod[1].stdName = "atmosphere_optical_thickness_due_to_ambient_aerosol";
    aod[1].ncTypeId = ncFloat.getId();
    aod[1].unit = "1";
    aod[1].setValidRange(0., 4.);
    aod[1].setFillValue(-999.);
    aod[1].vec.reserve(nPix);
    
    aod[2].varName = "AOD865";
    aod[2].longName = "aerosol optical thickness at 865 nm";
    aod[2].stdName = "atmosphere_optical_thickness_due_to_ambient_aerosol";
    aod[2].ncTypeId = ncFloat.getId();
    aod[2].unit = "1";
    aod[2].setValidRange(0., 4.);
    aod[2].setFillValue(-999.);
    aod[2].vec.reserve(nPix);
    
    aod[3].varName = "AOD1600";
    aod[3].longName = "aerosol optical thickness at 1610 nm";
    aod[3].stdName = "atmosphere_optical_thickness_due_to_ambient_aerosol";
    aod[3].ncTypeId = ncFloat.getId();
    aod[3].unit = "1";
    aod[3].setValidRange(0., 4.);
    aod[3].setFillValue(-999.);
    aod[3].vec.reserve(nPix);
    
    aod[4].varName = "AOD2250";
    aod[4].longName = "aerosol optical thickness at 2250 nm";
    aod[4].stdName = "atmosphere_optical_thickness_due_to_ambient_aerosol";
    aod[4].ncTypeId = ncFloat.getId();
    aod[4].unit = "1";
    aod[4].setValidRange(0., 4.);
    aod[4].setFillValue(-999.);
    aod[4].vec.reserve(nPix);
    
    aodUnc[0].varName = "AOD550_uncertainty";
    aodUnc[0].longName = "uncertainty on AOT at 550 nm";
    aodUnc[0].stdName = "atmosphere_optical_thickness_due_to_ambient_aerosol";
    aodUnc[0].ncTypeId = ncFloat.getId();
    aodUnc[0].unit = "1";
    aodUnc[0].setValidRange(0., 10.);
    aodUnc[0].setFillValue(-999.);
    aodUnc[0].vec.reserve(nPix);
    
    aodUnc[1].varName = "AOD670_uncertainty";
    aodUnc[1].longName = "uncertainty on AOT at 670 nm";
    aodUnc[1].stdName = "atmosphere_optical_thickness_due_to_ambient_aerosol";
    aodUnc[1].ncTypeId = ncFloat.getId();
    aodUnc[1].unit = "1";
    aodUnc[1].setValidRange(0., 10.);
    aodUnc[1].setFillValue(-999.);
    aodUnc[1].vec.reserve(nPix);
    
    aodUnc[2].varName = "AOD865_uncertainty";
    aodUnc[2].longName = "uncertainty on AOT at 865 nm";
    aodUnc[2].stdName = "atmosphere_optical_thickness_due_to_ambient_aerosol";
    aodUnc[2].ncTypeId = ncFloat.getId();
    aodUnc[2].unit = "1";
    aodUnc[2].setValidRange(0., 10.);
    aodUnc[2].setFillValue(-999.);
    aodUnc[2].vec.reserve(nPix);
    
    aodUnc[3].varName = "AOD1600_uncertainty";
    aodUnc[3].longName = "uncertainty on AOT at 1600 nm";
    aodUnc[3].stdName = "atmosphere_optical_thickness_due_to_ambient_aerosol";
    aodUnc[3].ncTypeId = ncFloat.getId();
    aodUnc[3].unit = "1";
    aodUnc[3].setValidRange(0., 10.);
    aodUnc[3].setFillValue(-999.);
    aodUnc[3].vec.reserve(nPix);
    
    aodUnc[4].varName = "AOD2250_uncertainty";
    aodUnc[4].longName = "uncertainty on AOT at 2250 nm";
    aodUnc[4].stdName = "atmosphere_optical_thickness_due_to_ambient_aerosol";
    aodUnc[4].ncTypeId = ncFloat.getId();
    aodUnc[4].unit = "1";
    aodUnc[4].setValidRange(0., 10.);
    aodUnc[4].setFillValue(-999.);
    aodUnc[4].vec.reserve(nPix);
    /**/
    
    fmAod.varName = "FM_AOD550";
    fmAod.longName = "fine-mode aerosol optical thickness at 550nm";
    fmAod.stdName = "atmosphere_optical_thickness_due_to_ambient_aerosol";
    fmAod.ncTypeId = ncFloat.getId();
    fmAod.unit = "1";
    fmAod.setValidRange(0., 4.);
    fmAod.setFillValue(-999.);
    fmAod.vec.reserve(nPix);
    
    dustAod.varName = "D_AOD550";
    dustAod.longName = "dust aerosol optical thickness at 550nm";
    dustAod.stdName = "atmosphere_optical_thickness_due_to_ambient_aerosol";
    dustAod.ncTypeId = ncFloat.getId();
    dustAod.unit = "1";
    dustAod.setValidRange(0., 4.);
    dustAod.setFillValue(-999.);
    dustAod.vec.reserve(nPix);
    
    absAod.varName = "AAOD550";
    absAod.longName = "aerosol absorption optical thickness at 550nm";
    absAod.stdName = "atmosphere_optical_thickness_due_to_ambient_aerosol";
    absAod.ncTypeId = ncFloat.getId();
    absAod.unit = "1";
    absAod.setValidRange(0., 4.);
    absAod.setFillValue(-999.);
    absAod.vec.reserve(nPix);
    
    ssa.varName = "SSA550";
    ssa.longName = "single scattering albedo at 550nm";
    //ssa.stdName = "";
    ssa.ncTypeId = ncFloat.getId();
    ssa.unit = "1";
    ssa.setValidRange(0., 1.);
    ssa.setFillValue(-999.);
    ssa.vec.reserve(nPix);
    
    lat.varName = "latitude";
    lat.longName = "Latitude at pixel centre";
    lat.stdName = "latitude";
    lat.ncTypeId = ncFloat.getId();
    lat.unit = "degrees_north";
    lat.setValidRange(-90., 90.);
    lat.setFillValue(-999.);
    lat.vec.reserve(nPix);
    
    lon.varName = "longitude";
    lon.longName = "Longitude at pixel centre";
    lon.stdName = "longitude";
    lon.ncTypeId = ncFloat.getId();
    lon.unit = "degrees_east";
    lon.setValidRange(-180., 180.);
    lon.setFillValue(-999.);
    lon.vec.reserve(nPix);

    std::string stmp1[] = {"1", "2", "3", "4"};
    std::string stmp2[] = {"1st", "2nd", "3rd", "4th"};
    for (int i=0; i<4; i++){
        latCrn[i].varName = "pixel_corner_latitude" + stmp1[i];
        latCrn[i].longName = "latitude_" + stmp2[i] + "_corner";
        latCrn[i].stdName = "latitude";
        latCrn[i].ncTypeId = ncFloat.getId();
        latCrn[i].unit = "degrees_north";
        latCrn[i].setValidRange(-90., 90.);
        latCrn[i].setFillValue(-999.);
        latCrn[i].vec.reserve(nPix);

        lonCrn[i].varName = "pixel_corner_longitude" + stmp1[i];
        lonCrn[i].longName = "longitude_" + stmp2[i] + "_corner";
        lonCrn[i].stdName = "longitude";
        lonCrn[i].ncTypeId = ncFloat.getId();
        lonCrn[i].unit = "degrees_east";
        lonCrn[i].setValidRange(-90., 90.);
        lonCrn[i].setFillValue(-999.);
        lonCrn[i].vec.reserve(nPix);
    }

    sza.varName = "sun_zenith_at_center";
    sza.longName = "solar zenith angle";
    sza.stdName = "solar_zenith_angle";
    sza.ncTypeId = ncFloat.getId();
    sza.unit = "degrees";
    sza.setValidRange(0., 90.);
    sza.setFillValue(-999.);
    sza.vec.reserve(nPix);

    vza.varName = "satellite_zenith_at_center";
    vza.longName = "satellite zenith angle";
    vza.stdName = "zenith_angle";
    vza.ncTypeId = ncFloat.getId();
    vza.unit = "degrees";
    vza.setValidRange(0., 90.);
    vza.setFillValue(-999.);
    vza.vec.reserve(nPix);

    raz.varName = "relative_azimuth_at_center";
    raz.longName = "relative azimuth angle";
    raz.ncTypeId = ncFloat.getId();
    raz.unit = "degrees";
    raz.setValidRange(0., 180.);
    raz.setFillValue(-999.);
    raz.vec.reserve(nPix);

    for (int i=0; i<N_SLSTR_BANDS; i++) {
        sdr[i].varName = "surface_reflectance" + wlStr[i];
        sdr[i].longName = "mean bidirectional surface reflectance (nadir)";
        sdr[i].ncTypeId = ncFloat.getId();
        sdr[i].unit = "1";
        sdr[i].setValidRange(0., 2.);
        sdr[i].setFillValue(-999.);
        sdr[i].vec.reserve(nPix);
    }
    
    time.varName = "time";
    time.longName = "time seconds since 1970-01-01 00:00:00 UTC";
    time.stdName = "time";
    time.ncTypeId = ncInt.getId();
    time.unit = "1";
    time.setValidRange(1, 2147483647);
    time.setFillValue(0);
    time.vec.reserve(nPix);
}

S3CciL2Writer::writeGlobalAttributes(netCDF::NcFile& ncFile) {
    
    const S3MetaData& s3MD = s3Data.pars.s3MD;
    const boost::uuids::uuid uuid = boost::uuids::random_generator()();
    const std::string uuidStr = boost::uuids::to_string(uuid);
    
    ncFile.putAtt("Conventions", "CF-1.6");
    ncFile.putAtt("tracking_id", uuidStr);
    ncFile.putAtt("naming_authority", "uk.ac.su.aatsraerosol");
    ncFile.putAtt("title", "AARDVARC CCI aerosol product level 2");
    ncFile.putAtt("product_version", "1.0");
    ncFile.putAtt("summary", "This dataset contains the level-2 aerosol properties products from SLSTR satellite observations. Data are processed by Swansea algorithm");
    ncFile.putAtt("id", ncFile.getName(true));
    ncFile.putAtt("sensor", s3MD.instrument);
    ncFile.putAtt("platform", s3MD.platformFamily + s3MD.platformNumber);
    ncFile.putAtt("resolution", "10km x 10km");
    ncFile.putAtt("projection", "quasi-cartesian instrument grid");
    ncFile.putAtt("cdm_data_type", "Swath");
    ncFile.putAtt("cell", "");
    ncFile.putAtt("inputFileList", s3MD.productName);
    ncFile.putAtt("time_coverage_start", s3MD.startTime);
    ncFile.putAtt("time_coverage_end", s3MD.stopTime);
    ncFile.putAtt("geospatial_lat_min", "");
    ncFile.putAtt("geospatial_lat_max", "");
    ncFile.putAtt("geospatial_lon_min", "");
    ncFile.putAtt("geospatial_lon_max", "");
    ncFile.putAtt("date_created", "");
    ncFile.putAtt("project", "Climate Change Initiative - European Space Agency");
    ncFile.putAtt("references", "http://www.esa-aerosol-cci.org");
    ncFile.putAtt("creator_name", "College of Science, Swansea University");
    ncFile.putAtt("creator_url", "http://www.swan.ac.uk/staff/academic/environmentsociety/geography/northpeter/");
    ncFile.putAtt("creator_email", "p.r.j.north@swansea.ac.uk \n a.heckel@swansea.ac.uk");
    ncFile.putAtt("source", "");
    ncFile.putAtt("keywords", "satellite,observation,atmosphere");
    ncFile.putAtt("keywords_vocabulary", "NASA Global Change Master Directory (GCMD) Science Keywords");
    ncFile.putAtt("standard_name_vocabulary", "NetCDF Climate and Forecast (CF) Metadata Convention version 18");
    ncFile.putAtt("license", "ESA CCI Data Policy: free and open access");
    //ncFile.putAtt("", "");
}
