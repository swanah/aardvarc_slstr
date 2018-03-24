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

#include <cmath>
#include <netcdf>
#include <stdexcept>
#include "OceanReflLut.hpp"

using namespace netCDF;

//
// public
//


OceanReflLut::OceanReflLut(const std::string& lutName) {
    
    NcFile ncF(lutName, NcFile::read);

    readDimVarFloat(&ncF, "wdsp", wdspD);
    readDimVarFloat(&ncF, "wdir", wdirD);
    readDimVarFloat(&ncF, "pigc", pigcD);
    readDimVarInt(&ncF, "model", modelD);
    readDimVarFloat(&ncF, "band", wvlD);
    readDimVarFloat(&ncF, "tau", aodD);
    readDimVarFloat(&ncF, "raz", razD);
    readDimVarFloat(&ncF, "vza", vzaD);
    readDimVarFloat(&ncF, "sza", szaD);
    readLutVar(&ncF, "Rocean");
    
    //ncF.close();
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
    if (val < dimPar.min || val > dimPar.max) {
        throw std::range_error("variable not in range of Ocn-LUT");
    }
    return (val == dimPar.max) ? (val * 0.9999 - dimPar.min) / dimPar.delta: (val - dimPar.min) / dimPar.delta;
}

float OceanReflLut::getInterPar(const float val, const LutDimParamInt& dimPar){
    if (val < dimPar.min || val > dimPar.max) {
        throw std::range_error("variable not in range of LUT");
    }
    return (val == dimPar.max) ? (val * 0.9999 - dimPar.min) / dimPar.delta: (val - dimPar.min) / dimPar.delta;
}

/**
 * writes ocean surface reflectance into SlstrPixel.rho_surf[band][view]
 * !! Relies on global index variables being setup previously
 * @param p
 * @param tau
 */
void OceanReflLut::get_rho_ocean(SlstrPixel* p, float tau) { 
    int band;
    float taui;

    taui = getInterPar(tau, aodD); //(tau - TAU_MIN) / T_OCEAN_INT;
    for (band = 0; band < N_SLSTR_BANDS; band++) {
            p->rho_surf[band][0] = interpol_ocean(p->lutpars.szani, p->lutpars.vzani, p->lutpars.ocn_razni, band, taui, p->lutpars.ocn_mi, p->lutpars.ocn_pi, p->lutpars.ocn_wdi, p->lutpars.ocn_wsi);
            p->rho_surf[band][1] = interpol_ocean(p->lutpars.szafi, p->lutpars.vzafi, p->lutpars.ocn_razfi, band, taui, p->lutpars.ocn_mi, p->lutpars.ocn_pi, p->lutpars.ocn_wdi, p->lutpars.ocn_wsi);
    }

}

/**
 * writes ocean surface reflectance into SlstrPixel.rho_surf[band][view]
 * !! Relies on global index variables being setup previously
 * @param p - slstrPixel
 * @param tau - aod
 * @param wndspd - wind speed
 */
void OceanReflLut::get_rho_ocean_wind(SlstrPixel *p, float tau, float wndspd) { // Relies on global index variables being setup previously
    int band;
    float taui, wi;

    taui = getInterPar(tau, aodD); //(tau - TAU_MIN) / T_OCEAN_INT;
    wi = getInterPar(wndspd, wdspD); //(wind - WDSP_MIN) / WDSP_INT;
    for (band = 0; band < N_SLSTR_BANDS; band++) {
            p->rho_surf[band][0] = interpol_ocean(p->lutpars.szani, p->lutpars.vzani, p->lutpars.ocn_razni, band, taui, p->lutpars.ocn_mi, p->lutpars.ocn_pi, p->lutpars.ocn_wdi, wi);
            p->rho_surf[band][1] = interpol_ocean(p->lutpars.szafi, p->lutpars.vzafi, p->lutpars.ocn_razfi, band, taui, p->lutpars.ocn_mi, p->lutpars.ocn_pi, p->lutpars.ocn_wdi, wi);
    }

}

/**
 * 
 * @param m_error - pointer to sqr error to be returned [band][angle]
 * @param p - slstrpixel
 * @param tau - aod
 * @param d_wind_speed - delta wind speed (for error)
 * @param d_pigment - delta pigment conc (for error)
 */
void OceanReflLut::get_rho_ocean_sqr_error(float* m_error, SlstrPixel* p, float tau, float d_wind_speed, float d_pigment) {
    int band;
    float taui, wi, pi;
    float rsurf_w, rsurf_p;

    if (tau < aodD.min) tau = aodD.min;
    if (tau > aodD.max) tau = aodD.max;
    taui = getInterPar(tau, aodD); //(tau - TAU_MIN) / T_OCEAN_INT;
    wi = getInterPar(p->ocn_wind_speed + d_wind_speed, wdspD); //;(p->ocn_wind_speed + d_wind_speed - WDSP_MIN) / WDSP_INT;
    pi = getInterPar(p->ocn_pigment + d_pigment, pigcD); //;(p->ocn_pigment + d_pigment - PIG_MIN) / PIG_INT;
    
    for (band = 0; band < N_SLSTR_BANDS; band++) {
        rsurf_w = interpol_ocean(p->lutpars.szani, p->lutpars.vzani, p->lutpars.ocn_razni, band, taui, p->lutpars.ocn_mi, p->lutpars.ocn_pi, p->lutpars.ocn_wdi, wi);
        rsurf_p = interpol_ocean(p->lutpars.szani, p->lutpars.vzani, p->lutpars.ocn_razni, band, taui, p->lutpars.ocn_mi, pi, p->lutpars.ocn_wdi, p->lutpars.ocn_wsi);
        p->rho_surf[band][0] = interpol_ocean(p->lutpars.szani, p->lutpars.vzani, p->lutpars.ocn_razni, band, taui, p->lutpars.ocn_mi, p->lutpars.ocn_pi, p->lutpars.ocn_wdi, p->lutpars.ocn_wsi);
        m_error[band * N_SLSTR_VIEWS] = pow((rsurf_w - p->rho_surf[band][0]),2);
        m_error[band * N_SLSTR_VIEWS] += pow((rsurf_p - p->rho_surf[band][0]),2);
        
        rsurf_w = interpol_ocean(p->lutpars.szafi, p->lutpars.vzafi, p->lutpars.ocn_razfi, band, taui, p->lutpars.ocn_mi, p->lutpars.ocn_pi, p->lutpars.ocn_wdi, wi);
        rsurf_p = interpol_ocean(p->lutpars.szafi, p->lutpars.vzafi, p->lutpars.ocn_razfi, band, taui, p->lutpars.ocn_mi, pi, p->lutpars.ocn_wdi, p->lutpars.ocn_wsi);
        p->rho_surf[band][1] = interpol_ocean(p->lutpars.szafi, p->lutpars.vzafi, p->lutpars.ocn_razfi, band, taui, p->lutpars.ocn_mi, p->lutpars.ocn_pi, p->lutpars.ocn_wdi, p->lutpars.ocn_wsi);
        m_error[band * N_SLSTR_VIEWS + 1] = pow((rsurf_w - p->rho_surf[band][1]),2);
        m_error[band * N_SLSTR_VIEWS + 1] += pow((rsurf_p - p->rho_surf[band][1]),2);
    }

}


//
// private
//

void OceanReflLut::readDimVarInt(NcFile *ncF, std::string varName, LutDimParamInt& dimPar){
    NcVar v;
    v = ncF->getVar(varName);
    dimPar.n = v.getDim(0).getSize();
    int a[dimPar.n];
    v.getVar(&a);
    dimPar.min = a[0];
    dimPar.max = a[dimPar.n-1];
    dimPar.delta = a[1] - a[0]; 
}

void OceanReflLut::readDimVarFloat(NcFile *ncF, std::string varName, LutDimParamFloat& dimPar){
    NcVar v;
    v = ncF->getVar(varName);
    dimPar.n = v.getDim(0).getSize();
    float a[dimPar.n];
    v.getVar(&a);
    dimPar.min = a[0];
    dimPar.max = a[dimPar.n-1];
    dimPar.delta = a[1] - a[0]; 
}

void OceanReflLut::readLutVar(NcFile *ncF, std::string varName){
    NcVar v;
    v = ncF->getVar(varName);
    nRsurf = 1;
    for (int iDim = v.getDimCount()-1; iDim >= 0; iDim--){
        nRsurf *= v.getDim(iDim).getSize();
    }
    rsurf = new float[nRsurf];
    v.getVar(rsurf);
}

double OceanReflLut::interpol_ocean(float szai, float vzai, float razi, int iBand, float taui, float mi, float pigi, float wdi, float wsi) {

    double t, u, w, z, v, s, r, q;
    int i, j, k, l, m, n, o, p, x1i[2], x2i[2], x3i[2], x4i[2], x5i[2], x6i[2], x7i[2], x8i[2];
    double val = 0.0;

    t = szai - (int) szai;
    u = vzai - (int) vzai;
    w = razi - (int) razi;
    z = taui - (int) taui;
    v = mi   - (int) mi;
    s = pigi - (int) pigi;
    r = wdi  - (int) wdi;
    q = wsi  - (int) wsi;

    for (i = 0; i < 2; i++) {
        x1i[i] = i + (int) szai;
        x2i[i] = i + (int) vzai;
        x3i[i] = i + (int) razi;
        x4i[i] = i + (int) taui;
        x5i[i] = i + (int) mi;
        x6i[i] = i + (int) pigi;
        x7i[i] = i + (int) wdi;
        x8i[i] = i + (int) wsi;
    }

    // Sorting out edge of LUT, so that
    // interpolation does not go over end

    if (x1i[1] == szaD.n) x1i[1] = x1i[0];
    if (x2i[1] == vzaD.n) x2i[1] = x2i[0];
    if (x3i[1] == razD.n) x3i[1] = x3i[0];
    if (x4i[1] == aodD.n) x4i[1] = x4i[0];
    if (x5i[1] == modelD.n) x5i[1] = x5i[0];
    if (x6i[1] == pigcD.n) x6i[1] = x6i[0];
    if (x7i[1] == wdirD.n) x7i[1] = x7i[0];
    if (x8i[1] == wdspD.n) x7i[1] = x8i[0];


    for (i = 0; i < 2; i++)  // szai, x1, t
        for (j = 0; j < 2; j++)  // vzai, x2, u
            for (k = 0; k < 2; k++)  // razi, x3, w
                for (l = 0; l < 2; l++)  // taui, x4, z
                    for (m = 0; m < 2; m++)  // mi, x5, v
                        for (n = 0; n < 2; n++)  // pigi, x6, s
                            for (o = 0; o < 2; o++)  // wdi, x7, r
                                for (p = 0; p < 2; p++)  // wsi, x8, q
                                    val += fabs(1 - i - t) * fabs(1 - j - u) * fabs(1 - k - w) * fabs(1 - l - z) * fabs(1 - m - v) * fabs(1 - n - s) * fabs(1 - o - r) * fabs(1 - p - q)
                                           * rsurf[x8i[p] + wdspD.n * (x7i[o] + wdirD.n * (x6i[n] + pigcD.n * (x5i[m] + modelD.n * (x4i[l] + aodD.n * (iBand + wvlD.n * (x3i[k] + razD.n * (x2i[j] + vzaD.n * x1i[i])))))))];
                                            //* luw[x1i[i]][x2i[j]][x3i[k]][band][x4i[l]][x5i[m]][x6i[n]][x7i[o]];

    return val;
}

