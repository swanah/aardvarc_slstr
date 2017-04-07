/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   OceanReflLut.hpp
 * Author: akheckel
 *
 * Created on 29. März 2017, 13:50
 */

#ifndef OCEANREFLLUT_HPP
#define OCEANREFLLUT_HPP

#include <netcdf>
#include "defs.hpp"



class OceanReflLut {
public:
    float* rsurf;
    size_t nRsurf;
    LutDimParamFloat szaD;
    LutDimParamFloat vzaD;
    LutDimParamFloat razD;
    LutDimParamFloat aodD;
    LutDimParamFloat wvlD;
    LutDimParamByte  modelD;
    LutDimParamFloat pigcD;
    LutDimParamFloat wdirD;
    LutDimParamFloat wdspD;
    
    OceanReflLut(const std::string& lutName);
    ~OceanReflLut();
    
    float getInterPar(const float val, const LutDimParamFloat& dimPar);
    float getInterPar(const float val, const LutDimParamByte&  dimPar);
    
private:
    OceanReflLut();
    OceanReflLut(const OceanReflLut& orig); //disabled copy constructor
    OceanReflLut& operator=(const OceanReflLut& rhs) { throw std::logic_error("assigning AeroClimatology not implemented"); } // disabled assignment 
    
    readDimVarByte(netCDF::NcFile *ncF, std::string varName, LutDimParamByte& dimPar);
    readDimVarFloat(netCDF::NcFile *ncF, std::string varName, LutDimParamFloat& dimPar);
    readLutVar(netCDF::NcFile *ncF, std::string varName);
};

#endif /* OCEANREFLLUT_HPP */

