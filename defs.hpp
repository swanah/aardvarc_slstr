/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   defs.hpp
 * Author: akheckel
 *
 * Created on 25. November 2016, 10:43
 */

#ifndef DEFS_HPP
#define DEFS_HPP

#include <cstdio>

#ifndef M_PI
#define M_PI		3.14159265358979323846
#endif

#define N_SLSTR_BANDS 5
#define N_SLSTR_VIEWS 2
#define N_DETECTORS 4

#define N_AER_FRAC 4

#define N_CCI_LUT 35

#define   DEFAULT_PALT 1013.0 // Default surface pressure
#define   DEFAULT_WDSP   3.0 // Default wind speed
#define   DEFAULT_WDIR  90.0 // Wind Direction (relative to SAA???)
#define   DEFAULT_PIG    0.1 // Default ocn pigment concentration (mg/m3)
#define   GLINT_THRS   0.008 // 1.6um threshold on modelled ocean refl. due to wind

#define   SALTCONC  34.3  // Salt concentration (Default = 34.3 ppt)

// AOD retrieval pixel flags (unsigned short = 16 bit)
#define CLR_LAND 1                 // Bit 1
#define CLR_OCEAN_N 2              // Bit 2
#define CLR_OCEAN_F 4              // Bit 3
//#define UNIFORM_0870_FAILED 4      // Bit 3
#define CCI_MIX_FAILED 8           // Bit 4
#define TOA_ARR_ZERO 16            // Bit 5
//#define RR_NEGATIVE 32             // Bit 6
#define AOD_HIGH 32                // Bit 6
#define AOD_ZERO 64                // Bit 7
#define EMOD_PENALTY 128           // Bit 8
#define CURV_NEGATIVE 256          // Bit 9
#define FINE_FRAC_FORCED 512       // Bit 10
#define FIT_AOD_ZERO_CLIM 1024     // Bit 11
#define FIT_AOD_ZERO_FINE 2048     // Bit 12
#define FIT_AOD_ZERO_COARSE 4096   // Bit 13
#define GLINT_NADIR         8192   // Bit 14
#define GLINT_FWARD        16384   // Bit 15



struct ImageProperties {
    int width, height, nPix;
    int xOff, yOff;
    int xRes, yRes;

    bool isBinned;
    int binSize;

    ImageProperties() {
        isBinned = false;
        binSize = 1;
        width = -1;
        height = -1;
        nPix = 0;
        xOff = 0;
        yOff = 0;
        xRes = -1;
        yRes = -1;
    }

    bool operator==(const ImageProperties& rhs) {
        return ( (width == rhs.width)
                &&(height == rhs.height)
                &&(nPix == rhs.nPix)
                &&(xOff == rhs.xOff)
                &&(yOff == rhs.yOff)
                &&(xRes == rhs.xRes)
                &&(yRes == rhs.yRes)
                &&(isBinned == rhs.isBinned)
                &&(binSize == binSize));
    }

    bool operator==(const ImageProperties& rhs) const {
        return ( (width == rhs.width)
                &&(height == rhs.height)
                &&(nPix == rhs.nPix)
                &&(xOff == rhs.xOff)
                &&(yOff == rhs.yOff)
                &&(xRes == rhs.xRes)
                &&(yRes == rhs.yRes)
                &&(isBinned == rhs.isBinned)
                &&(binSize == binSize));
    }
};

struct LutDimParamByte {
    char min, max, delta;
    int n;
};

struct LutDimParamInt {
    int min, max, delta;
    int n;
};

struct LutDimParamFloat {
    float min, max, delta;
    int n;
};

struct LutDimParamDouble {
    double min, max, delta;
    int n;
};

struct LutPars {
    float szani, szafi, vzani, vzafi, razni, razfi;
    float pAlti, o3i, ocn_mi, ocn_wsi, ocn_wdi, ocn_pi;
    float ocn_razni, ocn_razfi;
    float mixing[4];  // Dust, SeaSalt, StrongAbs, WeakAbs
    int tetraP[4][4]; // 
    float mix_frac[3]; // fine_of_total, weak_of_fine, dust_of_coarse
    float climAod;      // climatology AOD
    float climFineMode; // climatology fine mode fraction
};

struct GeoPos {
    float lat, lon;
    char latLonStr[30];

    GeoPos() {}
    GeoPos(float lat, float lon) : lat(lat), lon(lon) {}
    char* toCstr(){
        sprintf(latLonStr, "(%8.4fN/%9.4fE)", lat, lon);
        return latLonStr;
    }

    GeoPos(const GeoPos& other) : lat(other.lat), lon(other.lon) {}
};

struct ViewGeom {
    float for_sat_azim, for_sol_azim, nad_sat_azim, nad_sol_azim;
    float for_sol_zen, for_sat_zen, nad_sol_zen, nad_sat_zen;
    float razn, razf;
};

struct LutCoefs {
    //double fd, cfac, xa, xb, xc;
    double rPath, tDown, tUp, tGas, spherAlb, difFrac;
};

struct SlstrPixel{
    ViewGeom geom;
    GeoPos geo_pos;
    float pAlt;
    float o3;
    float ocn_wind_speed;
    float ocn_wind_dir;
    float ocn_pigment;
    float fclr_lnd;
    float fclr_ocn[N_SLSTR_VIEWS];
    bool  view_clear[N_SLSTR_VIEWS];
    char mixOK, range;
    char dims;
    unsigned short qflag;
    LutPars lutpars;
    LutCoefs inv_coef[N_SLSTR_BANDS][N_SLSTR_VIEWS]; // inversion coefficients [band][angle]
    float DD[N_SLSTR_BANDS][N_SLSTR_VIEWS];       // diffuse fraction per band and angle
    float RR[N_SLSTR_BANDS][N_SLSTR_VIEWS];       // BRF or SDR  per band and angle
    float tarr[N_SLSTR_BANDS][N_SLSTR_VIEWS];     // TOA reflec  per band and angle
    float rho_surf[N_SLSTR_BANDS][N_SLSTR_VIEWS]; // modelled ocean reflec from LUT per band and angle
    float angWeights[N_SLSTR_BANDS][N_SLSTR_VIEWS]; // fit weights in angular model per band and angle
    float aod, fmin, ediff;                       // retrieved aod, fit residual and uncertainty
    float spec_aod_fac[N_SLSTR_BANDS];            // factors for spectral AOD, i.e.(AOD[BAND]/AOD[550])
    float ssa[N_SLSTR_BANDS];                     // single scattering albedo
    float ndvi;                                   // ndvi
    float ndvi_veg_weight;                        // 
    float dust_weight;                            // 
    float rho_glint[N_SLSTR_VIEWS];               // modelled rho ocean for glint testing (R(1.61) at wind speed of 9m/s)
    //float(*fct_emod_tau)(float);                  // fct. pointer to actual aerosol fit routine
    bool dumpPix, dumpRR;                         // switch to dump Pix or just RR to output for debugging
    int x, y;                                     // x/y coord of pixel in image
    float ax, cx;                                 // brent brackets
    float prevFineFrac;                           // previously retrieved fine mode fraction
};


#endif /* DEFS_HPP */

