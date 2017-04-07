/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   AtmosphericLut.hpp
 * Author: akheckel
 *
 * Created on 30. März 2017, 15:49
 */

#ifndef ATMOSPHERICLUT_HPP
#define ATMOSPHERICLUT_HPP

#include <map>
#include "defs.hpp"

#define CCI_INT 25
#define CCI_NUM  5
#define BEST 0
#define LATEST 1

class AtmosphericLut {
public:
    double* rPath;          // dims: sza, vza, raz, pres, aod, wvl, model
    size_t nRPath;
    double* t;              // dims: sza, pres, aod, wvl, model
    size_t nT;
    double* tGas;           // dims: sza, vza, pres, wvl
    size_t nTgas;
    double* spherA;         // dims: pres, tau, wvl, model
    size_t nSpherA;
    double* difFrac;        // dims: sza, pres, aod, wvl, model
    size_t nDifFrac;
    double* specAodRatio;   // dims: wvl, model
    size_t nSpecAodRatio;
    double* ssa;            // dims: wvl, model
    size_t nSsa;
    
    LutDimParamDouble szaD;
    LutDimParamDouble vzaD;
    LutDimParamDouble razD;
    LutDimParamDouble presD;
    LutDimParamDouble aodD;
    LutDimParamDouble wvlD;
    LutDimParamInt  modelD;

    std::map<int,int> aerModelIdxMap;
    
    AtmosphericLut(const std::string& lutName);
    ~AtmosphericLut();
    
    char getTetrahedronPoints(LutPars *lutpars, int verbose);
    double getTetraVol(int v0, int v1, int v2, int vOrigin, LutPars* lutpars);
    void psInv6s(SlstrPixel* pix, float tau);
    float getInterPar(const float val, const LutDimParamDouble& dimPar);
    double interpolRpath(const float szai, const float vzai, const float razi, const float presi, const float taui, const int iBand, const int iModel);
    
    
private:
    AtmosphericLut();
    AtmosphericLut(const AtmosphericLut& orig); //disabled copy constructor
    AtmosphericLut& operator=(const AtmosphericLut& rhs) { throw std::logic_error("assigning atmospheric LUT not implemented"); } // disabled assignment 

    void readDimVarByte(netCDF::NcFile *ncF, std::string varName, LutDimParamByte& dimPar);
    void readDimVarInt(netCDF::NcFile *ncF, std::string varName, LutDimParamInt& dimPar);
    void readDimVarFloat(netCDF::NcFile *ncF, std::string varName, LutDimParamFloat& dimPar);
    void readDimVarDouble(netCDF::NcFile *ncF, std::string varName, LutDimParamDouble& dimPar);
    void readLutVar(netCDF::NcFile *ncF, std::string varName, double** data, size_t* nData);
    void initIndexMap();
    void getCoefs(LutCoefs* coef, const LutPars* lpars, const float taui, const int iBand, const int iModel);
    double interpolSpherA(const float presi, const float taui, const int iBand, const int iModel);
    double interpolDifFrac(const float szai, const float presi, const float taui, const int iBand, const int iModel);
    double interpolTgas(const float szai, const float vzai, const float presi, const int iBand, const int iModel);
    double interpolTscat(const float szai, const float presi, const float taui, const int iBand, const int iModel);
};

#endif /* ATMOSPHERICLUT_HPP */

