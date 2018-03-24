/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */ 

/* 
 * File:   S3CciL2Writer.hpp
 * Author: akheckel
 *
 * Created on 22. November 2017, 08:30
 */

#ifndef S3CCIL2WRITER_HPP
#define S3CCIL2WRITER_HPP

#include <netcdf>
#include <vector>
#include "defs.hpp"
#include "S3NcdfData.hpp"


template<class T>
class CciDataVector {
public:
    std::string varName;
    std::string longName;
    std::string stdName;
    std::string unit;
    nc_type ncTypeId;
    
    bool hasValidRange;
    T validRange[2];
    bool hasFillValue;
    T fillValue;
    
    std::vector<T> vec;


    void setValidRange(const T& min, const T& max) {
        hasValidRange = true;
        validRange[0] = min;
        validRange[1] = max;
    }

    void setFillValue(const T& fillVal) {
        hasValidRange = true;
        fillValue = fillVal;
    }
    

private:
            
};



class S3CciL2Writer {
public:
    int nPix;
    const S3NcdfData& s3Data;
    
    CciDataVector<int>   pixelNumber;
    CciDataVector<float> aod[N_SLSTR_BANDS];
    CciDataVector<float> aodUnc[N_SLSTR_BANDS];
    CciDataVector<float> fmAod;
    CciDataVector<float> dustAod;
    CciDataVector<float> absAod;
    CciDataVector<float> ssa;

    CciDataVector<float> lat;
    CciDataVector<float> lon;
    CciDataVector<float> latCrn[4];
    CciDataVector<float> lonCrn[4];

    CciDataVector<float> sza;
    CciDataVector<float> vza;
    CciDataVector<float> raz;

    CciDataVector<float> sdr[N_SLSTR_BANDS];

    CciDataVector<int> time;
    


    S3CciL2Writer(const S3NcdfData& s3NcdfData, const int& numberOfPixels);
    ~S3CciL2Writer();
    
    void write();
    
private:
    S3CciL2Writer();
    S3CciL2Writer(const S3CciL2Writer& orig);  //disabled copy Constructor
    S3CciL2Writer& operator=(const S3CciL2Writer& rhs){ throw std::logic_error("S3NcdfData shoudlnt be copied"); } // disable copy assignment
    
    template<class T>
    void addWriteVariable(netCDF::NcFile& ncFile, CciDataVector<T> varData, std::vector<netCDF::NcDim> dimVec);
    initNcdfVars();
    writeGlobalAttributes(netCDF::NcFile& ncFile);
};


template<class T>
void S3CciL2Writer::addWriteVariable(netCDF::NcFile& ncFile, CciDataVector<T> varData, std::vector<netCDF::NcDim> dimVec){
    long unsigned int nDim = dimVec.size();
    if (nDim <= 0){
        std::cerr << varData.varName << std::endl;
        throw std::runtime_error("adding var to netCDF - number of dimensions is < 1\n");
    }
    if (varData.vec.size() <= 0){
        std::cerr << varData.varName << std::endl;
        throw std::runtime_error("adding var to netCDF - var contains no data elements\n");
    }
    int nDimElem = 1;
    for (int iDim=0; iDim<dimVec.size(); iDim++) {
        nDimElem *= dimVec[iDim].getSize();
    }
    if (nDimElem != varData.vec.size()){
        std::cerr << varData.varName << std::endl;
        throw std::runtime_error("writing var to netCDF - mismatch in number of elements for data and netCDF\n");
    }
    
    netCDF::NcType nct(varData.ncTypeId);
    netCDF::NcVar var = ncFile.addVar(varData.varName, nct, dimVec);
    var.setCompression(true, true, 5);
    if (!varData.longName.empty()) {
        var.putAtt("long_name", varData.longName);
    }
    if (!varData.stdName.empty()) {
        var.putAtt("standard_name", varData.stdName);
    }
    if (!varData.unit.empty()) {
        var.putAtt("units", varData.unit);
    }
    if (varData.hasValidRange) {
        var.putAtt("valid_range", nct, 2, varData.validRange);
    }

    var.putVar(varData.vec.data());
}


#endif /* S3CCIL2WRITER_HPP */

