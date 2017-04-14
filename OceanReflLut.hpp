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
    LutDimParamFloat wvlD;
    LutDimParamFloat aodD;
    LutDimParamInt  modelD;
    LutDimParamFloat pigcD;
    LutDimParamFloat wdirD;
    LutDimParamFloat wdspD;
    
    OceanReflLut(const std::string& lutName);
    ~OceanReflLut();
    
    float getInterPar(const float val, const LutDimParamFloat& dimPar);
    float getInterPar(const float val, const LutDimParamInt&  dimPar);
    void get_rho_ocean(SlstrPixel* p, float tau);
    void get_rho_ocean_wind(SlstrPixel* p, float tau, float wndspd);
    void get_rho_ocean_sqr_error(float* m_error, SlstrPixel* p, float tau, float d_wind_speed, float d_pigment);
    double interpol_ocean(float szai, float vzai, float razi, int iBand, float taui, float iModel, float pigi, float wdi, float wsi);
    
private:
    OceanReflLut();
    OceanReflLut(const OceanReflLut& orig); //disabled copy constructor
    OceanReflLut& operator=(const OceanReflLut& rhs) { throw std::logic_error("assigning AeroClimatology not implemented"); } // disabled assignment 
    
    readDimVarInt(netCDF::NcFile *ncF, std::string varName, LutDimParamInt& dimPar);
    readDimVarFloat(netCDF::NcFile *ncF, std::string varName, LutDimParamFloat& dimPar);
    readLutVar(netCDF::NcFile *ncF, std::string varName);
};

#endif /* OCEANREFLLUT_HPP */

