/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   OceanReflLut.cpp
 * Author: akheckel
 * 
 * Created on 29. März 2017, 13:50
 */

#include <netcdf>
#include "OceanReflLut.hpp"

using namespace netCDF;

OceanReflLut::OceanReflLut(const std::string& lutName) {
    
    NcFile ncF(lutName, NcFile::read);

    readDimVarFloat(&ncF, "wdsp", wdspD);
    readDimVarFloat(&ncF, "wdir", wdirD);
    readDimVarFloat(&ncF, "pigc", pigcD);
    readDimVarByte(&ncF, "model", modelD);
    readDimVarFloat(&ncF, "band", wvlD);
    readDimVarFloat(&ncF, "tau", aodD);
    readDimVarFloat(&ncF, "raz", razD);
    readDimVarFloat(&ncF, "vza", vzaD);
    readDimVarFloat(&ncF, "sza", szaD);
    readLutVar(&ncF, "Rocean");
    
    ncF.close();
}

OceanReflLut::~OceanReflLut() {
    delete [] rsurf;
    rsurf = NULL;
}

/**
 * compute fractional index, includes range check and verifies the returned index is always smaller than dimPar.n-1
 * the last verification avoids testing in the lut interpolation!!! 
 * @param val - actual value along lut dim
 * @param dimPar - corresponding dim parameter
 * @return fractional index in dim axis
 */
float OceanReflLut::getInterPar(const float val, const LutDimParamFloat& dimPar){
    if (val < dimPar.min || val > dimPar.max) throw std::range_error("variable not in range of LUT");
    return (val == dimPar.max) ? (val * 0.9999 - dimPar.min) / dimPar.delta: (val - dimPar.min) / dimPar.delta;
}

float OceanReflLut::getInterPar(const float val, const LutDimParamByte& dimPar){
    if (val < dimPar.min || val > dimPar.max) throw std::range_error("variable not in range of LUT");
    return (val == dimPar.max) ? (val * 0.9999 - dimPar.min) / dimPar.delta: (val - dimPar.min) / dimPar.delta;
}


OceanReflLut::readDimVarByte(NcFile *ncF, std::string varName, LutDimParamByte& dimPar){
    NcVar v;
    v = ncF->getVar(varName);
    dimPar.n = v.getDim(0).getSize();
    char a[dimPar.n];
    v.getVar(&a);
    dimPar.min = a[0];
    dimPar.max = a[dimPar.n-1];
    dimPar.delta = a[1] - a[0]; 
}

OceanReflLut::readDimVarFloat(NcFile *ncF, std::string varName, LutDimParamFloat& dimPar){
    NcVar v;
    v = ncF->getVar(varName);
    dimPar.n = v.getDim(0).getSize();
    float a[dimPar.n];
    v.getVar(&a);
    dimPar.min = a[0];
    dimPar.max = a[dimPar.n-1];
    dimPar.delta = a[1] - a[0]; 
}

OceanReflLut::readLutVar(NcFile *ncF, std::string varName){
    NcVar v;
    v = ncF->getVar(varName);
    nRsurf = 1;
    for (int iDim = v.getDimCount()-1; iDim >= 0; iDim--){
        nRsurf *= v.getDim(iDim).getSize();
    }
    rsurf = new float[nRsurf];
    v.getVar(rsurf);
}

