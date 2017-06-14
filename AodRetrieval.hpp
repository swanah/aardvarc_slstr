/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   AodRetrieval.hpp
 * Author: akheckel
 *
 * Created on 8. April 2017, 14:25
 */

#ifndef AODRETRIEVAL_HPP
#define AODRETRIEVAL_HPP

#include <stdexcept>
#include "defs.hpp"
#include "AtmosphericLut.hpp"
#include "OceanReflLut.hpp"
#include "nr3.h"
#include "mins.h"
#include "mins_ndim.h"


//#define TOL 0.01    // optimisation limit for Brent fits
//#define FTOL 5e-4   // optimisation limit for Powell surface fit


//const float inst_frac_error[] = {.024, .032, .02, .033, .033}; /* Instrument calibration errors adjusted to uncertainty document specs (taken from ORAC)) */
//const float inst_frac_error[] = {.048, .064, .04, .066, .066}; /* Instrument calibration errors adjusted to uncertainty document specs (taken from ORAC)) */
//const float inst_frac_error[] = {.048, .064, .04, .12, .12}; /* Instrument calibration errors adjusted to uncertainty document specs (taken from ORAC)) */
const float inst_frac_error[] = {.048, .064, .04, .066, .12}; /* Instrument calibration errors adjusted to uncertainty document specs (taken from ORAC)) */
//const float m_error_offset[] = {0.005, 0.005, 0.02, 0.01}; /* Land surface model error estimates per channel, derived from MOnte Carlo simulations */
//const float m_error_offset[] = {0.01, 0.01, 0.04, 0.02}; /* Land surface model error estimates per channel, derived from MOnte Carlo simulations */
const float m_error_offset_v[] = {0.01, 0.01, 0.04, 0.02, 0.02}; /* Land surface model error estimates per channel, derived from MOnte Carlo simulations */
const float m_error_offset_b[] = {0.01, 0.01, 0.01, 0.15, 0.15}; /* Land surface model error estimates per channel, derived from MOnte Carlo simulations */
//const float m_error_gain[] = {0.07, 0.07, 0.15, 0.1}; /* Written as sum of absolute ('offset') and relative ('gain') errors, 1sd */
const float m_error_gain_v[] = {0.05, 0.05, 0.10, 0.07, 0.07}; /* Written as sum of absolute ('offset') and relative ('gain') errors, 1sd */
const float m_error_gain_b[] = {0.05, 0.05, 0.05, 0.2, 0.2}; /* Written as sum of absolute ('offset') and relative ('gain') errors, 1sd */
const float rt_abs_error[] = {.006, .006, .006, .006, .006}; /* Approximate uncertainty in surface reflectance from 6S code for known constituents (from Kalashnikova et al) */
const float aer_mod_frac_error[] = {0.05, 0.05, 0.05, 0.05, 0.05}; /* Fractional error per channel in aerosol scattering coefficient (multiple of SSA and phase fn errors */
//float aer_mod_frac_error[] = {0.01, 0.01, 0.01, 0.05}; /* Fractional error per channel in aerosol scattering coefficient (multiple of SSA and phase fn errors */
const float aer_mod_frac_error_d[] = {0.05, 0.05, 0.05, 0.2, 0.2}; /* Fractional error per channel in aerosol scattering coefficient (multiple of SSA and phase fn errors */
const float dToa = 0.01;
const float ocn_d_wind_speed = 6;
const float ocn_d_pigment = 0.1;


class EmodTau {
public:
    Doub fot;
    virtual ~EmodTau(){};
    virtual Doub operator()(const Doub) = 0;
};

class EmodSize {
public:
    SlstrPixel& pix;
    AtmosphericLut& atmLut;
    OceanReflLut& ocnLut;
    EmodTau* emodTau;

    EmodSize(SlstrPixel& pix, AtmosphericLut& atmLut, OceanReflLut& ocnLut, EmodTau* emod_tau) :
    pix(pix), atmLut(atmLut), ocnLut(ocnLut), emodTau(emod_tau) {}

    Doub operator()(const Doub fot) {
        Brent brentAod;
        brentAod.ax = pix.ax;
        brentAod.bx = 1.05 * pix.ax; //(pix.ax + pix.cx) / 2;
        brentAod.cx = pix.cx;
        emodTau->fot = fot;
        float fmin, penalty;
        pix.lutpars.mixing[0] = ((1 - fot) * pix.lutpars.mix_frac[2])*100; // Dust
        pix.lutpars.mixing[1] = ((1 - fot) * (1 - pix.lutpars.mix_frac[2]))*100; // Sea Salt
        pix.lutpars.mixing[2] = (fot * (1 - pix.lutpars.mix_frac[1]))*100; // Strong Abs
        pix.lutpars.mixing[3] = (fot * pix.lutpars.mix_frac[1])*100; // Weak Abs
        atmLut.getTetrahedronPoints(&pix.lutpars, false);
        pix.lutpars.ocn_mi = ocnLut.getInterPar(fot, ocnLut.modelD);
        brentAod.minimize(*emodTau); //brent_d(&ax, bx, &cx, fct_emod_tau, TOL, &pix.aod);
        fmin = brentAod.fmin;
        pix.aod = brentAod.xmin;
        penalty = (fot - pix.lutpars.mix_frac[0]);
        penalty *= penalty;
        //fmin += 5*penalty + 2*fmin*penalty ;
        //fmin += fmin * pow(penalty, 2);
        fmin += 25 * pow(penalty, 4);
        if (pix.prevFineFrac > 0) {
            penalty = (fot - pix.prevFineFrac);
            fmin += 25 * pow(penalty, 4);
        }
        
        return fmin;
    }
};

class EmodUncert {
public:
    SlstrPixel& pix;
    AtmosphericLut& atmLut;
    OceanReflLut& ocnLut;
    double *dRsurf_dToa;

    EmodUncert(SlstrPixel& pix, AtmosphericLut& atmLut, OceanReflLut& ocnLut, double *dRsurf__dToa) :
    pix(pix), atmLut(atmLut), ocnLut(ocnLut), dRsurf_dToa(dRsurf__dToa) {
    }

    Doub operator() (VecDoub_I &p){
        int i, j;
        double DF = 0.3, p6 = 0.35, tot = 0.0;
        double mval[N_SLSTR_BANDS][N_SLSTR_VIEWS];
        double k, dir, dif, g;
        double Rpath, y, m_error, o_error;
        double lim;
        double WGv[] = {1500.0, 1000.0, 1000.0, 1500.0};
        double WGe[] = {1000.0, 1000.0, 1000.0, 1000.0};
        double WGd[] = {1000.0, 1000.0, 1000.0, 0.0};
        double WG[]  = {1000.0, 1000.0, 1000.0, 100.0};


        for (i = 0; i < N_SLSTR_BANDS; i++) {
            for (j = 0; j < N_SLSTR_VIEWS; j++) {

                dir = (1.0 - DF * (double) pix.DD[i][j])*(double) p[5 + j]*(double) p[i];
                g = (1.0 - p6)*(double) p[i];
                dif = (DF * (double) pix.DD[i][j] + g * (1.0 - DF * (double) pix.DD[i][j])) * p6 * (double) p[i] / (1.0 - g);
                mval[i][j] = (dir + dif);
                k = (double) pix.RR[i][j] - mval[i][j];

                /* Model and observation error to calculate chi sq. Could include full error covariane matrix here */
                /* Model errors per channel derived from fit against Monte Carlo 3D model */
                /* Simple error propogation of atmospheric & sensor errors based on net atmospheric effect for solar/view geometry and opt depth */
                m_error = m_error_offset_v[i];
                //m_error = pow(m_error_offset_v[i], 2);
                //switch (i){
                //    case 2  : m_error = pix.ndvi_veg_weight * m_error_offset_v[i] + (1.0 - pix.ndvi_veg_weight) * m_error_offset_b[i]; break;
                //    case 3  : m_error = pix.dust_weight * m_error_offset_v[i] + (1.0 - pix.dust_weight) * m_error_offset_b[i]; break;
                //    default : m_error = m_error_offset_v[i]; break;
                //}
                m_error = m_error*m_error;//pow(m_error, 2);

                //m_error = pix.ndvi_veg_weight * m_error_offset_v[i] + (1.0 - pix.ndvi_veg_weight) * m_error_offset_b[i];
                //m_error = m_error_offset[i] + mval[i][j] * m_error_gain[i];
                //m_error = m_error_offset[i] + pix.RR[3][0] * m_error_gain[i];
                //m_error +=  pix.RR[i][j] * (pix.ndvi_veg_weight * m_error_gain_v[i] + (1.0 - pix.ndvi_veg_weight) * m_error_gain_b[i]);
                //m_error = m_error_offset[i] + pix.RR[i][j]*pix.RR[i][j] * m_error_gain[i];

                // compute o_error view dependent
                y = (pix.tarr[i][j] + dToa - pix.inv_coef[i][j].rPath) / pix.inv_coef[i][j].tGas;
                y /= (pix.inv_coef[i][j].tDown * pix.inv_coef[i][j].tUp);
                *dRsurf_dToa = fabs((y / (1.0 + pix.inv_coef[i][j].spherAlb * y)) - pix.RR[i][j]) / dToa;

                //y = (pix.tarr[i][j] + dToa) * pix.inv_coef[i][j].cfac * pix.inv_coef[i][j].xa - pix.inv_coef[i][j].xb;
                //dRsurf_dToa = fabs((y / (1.0 + y * pix.inv_coef[i][j].xc)) - pix.RR[i][j]) / dToa;

                o_error = (inst_frac_error[i] * pix.RR[i][j] * *dRsurf_dToa);
                //o_error = (0.5 * inst_frac_error[i] * pix.tarr[i][j] * *(dRsurf_dToa + i*N_SLSTR_VIEWS + j));
                o_error = o_error * o_error;

                //Rpath = pix.inv_coef[i][j].xb / (pix.inv_coef[i][j].xa * pix.inv_coef[i][j].cfac);
                //o_error += (aer_mod_frac_error[i] * Rpath);
                //o_error += pow((0.5 * aer_mod_frac_error[i] * Rpath), 2);
                //o_error += 0.5 * (((1 - pix.dust_weight) * aer_mod_frac_error[i] + pix.dust_weight * aer_mod_frac_error_d[i])* Rpath);

                // compute o_error with equal weights for forward and nadir
                /*y = (pix.tarr[i][0] + dToa) * pix.inv_coef[i][0].cfac * pix.inv_coef[i][0].xa - pix.inv_coef[i][0].xb;
                dRsurf_dToa = fabs((y / (1.0 + y * pix.inv_coef[i][0].xc)) - pix.RR[i][0]) / dToa;
                o_error += inst_frac_error[i] /10 * pix.RR[i][0] * dRsurf_dToa;
                //o_error += inst_frac_error[i] * pix.RR[i][j] * pix.RR[i][j] / pix.tarr[i][j];

                Rpath = pix.inv_coef[i][0].xb / (pix.inv_coef[i][0].xa * pix.inv_coef[i][0].cfac);
                o_error += aer_mod_frac_error[i] * Rpath;
                //o_error += aer_mod_frac_error[i] * fabs(pix.tarr[i][j] - pix.RR[i][j]);
                 */

                //o_error = rt_abs_error[i];
                o_error += rt_abs_error[i]*rt_abs_error[i];//pow(rt_abs_error[i],2);


                /* Replace WG to give chi sq. Could still use old WG in addition if normalised to 1 */
                //pix.angWeights[i][j] = 1 / (m_error * m_error + o_error * o_error);
                pix.angWeights[i][j] = 4 / (m_error + o_error);

                tot = tot + pix.angWeights[i][j] * k * k;

            }
        }

        lim = 0.03; if (p[0] < lim) tot = tot + (lim - p[0])*(lim - p[0])*1000.0;
        lim = 0.02; if (p[1] < lim) tot = tot + (lim - p[1])*(lim - p[1])*1000.0;
        lim = 0.01; if (p[2] < lim) tot = tot + (lim - p[2])*(lim - p[2])*1000.0;
        lim = 0.01; if (p[3] < lim) tot = tot + (lim - p[3])*(lim - p[3])*1000.0;
        lim = 0.01; if (p[4] < lim) tot = tot + (lim - p[4])*(lim - p[4])*1000.0;
        lim = 0.49; if (p[5] < lim) tot = tot + (lim - p[5])*(lim - p[5])*1000.0;
        lim = 0.51; if (p[5] > lim) tot = tot + (lim - p[5])*(lim - p[5])*1000.0;
        lim = 0.20; if (p[6] < lim) tot = tot + (lim - p[6])*(lim - p[6])*1000.0;
        //lim = 1.20; if (p[5] > lim) tot = tot + (lim - p[5])*(lim - p[5])*1000.0;
    /*
        lim = 0.01; if (p[0 + 1] < lim) tot = tot + (lim - p[0 + 1])*(lim - p[0 + 1])*1000.0;
        lim = 0.01; if (p[1 + 1] < lim) tot = tot + (lim - p[1 + 1])*(lim - p[1 + 1])*1000.0;
        lim = 0.01; if (p[2 + 1] < lim) tot = tot + (lim - p[2 + 1])*(lim - p[2 + 1])*1000.0;
        lim = 0.01; if (p[3 + 1] < lim) tot = tot + (lim - p[3 + 1])*(lim - p[3 + 1])*1000.0;
        lim = 0.49; if (p[4 + 1] < lim) tot = tot + (lim - p[4 + 1])*(lim - p[4 + 1])*1000.0;
        lim = 0.20; if (p[5 + 1] < lim) tot = tot + (lim - p[5 + 1])*(lim - p[5 + 1])*1000.0;

        lim = 0.70; if (p[0 + 1] > lim) tot = tot + (lim - p[0 + 1])*(lim - p[0 + 1])*1000.0;
        lim = 0.70; if (p[1 + 1] > lim) tot = tot + (lim - p[1 + 1])*(lim - p[1 + 1])*1000.0;
        lim = 0.70; if (p[2 + 1] > lim) tot = tot + (lim - p[2 + 1])*(lim - p[2 + 1])*1000.0;
        lim = 0.70; if (p[3 + 1] > lim) tot = tot + (lim - p[3 + 1])*(lim - p[3 + 1])*1000.0;
        lim = 1.00; if (p[4 + 1] > lim) tot = tot + (lim - p[4 + 1])*(lim - p[4 + 1])*1000.0;
        lim = 1.20; if (p[5 + 1] > lim) tot = tot + (lim - p[5 + 1])*(lim - p[5 + 1])*1000.0;
    */
    /*
        lim = 0.49; if (p[4 + 1] < lim) tot = tot + (lim - p[4 + 1])*(lim - p[4 + 1])*1000.0;
        lim = 0.20; if (p[5 + 1] < lim) tot = tot + (lim - p[5 + 1])*(lim - p[5 + 1])*1000.0;
        lim = 0.51; if (p[4 + 1] > lim) tot = tot + (lim - p[4 + 1])*(lim - p[4 + 1])*1000.0;
    */
        
        //fprintf(stderr, "pix: %d/%d fot: %f aod: %f\n", pix.x, pix.y, pix.lutpars.mixing[3]/pix.lutpars.mix_frac[1], pix.aod);
        //for (int iVec=0; iVec<p.size(); iVec++){
        //    fprintf(stderr, " %g ", p[iVec]);
        //}
        //fprintf(stderr, "\ntot: %g\n\n", tot); fflush(stderr);


        return (float) tot;

    }
};

class EmodOcean {
public:
    SlstrPixel& pix;
    AtmosphericLut& atmLut;
    OceanReflLut& ocnLut;
    double *dRsurf_dToa;

    EmodOcean(SlstrPixel& pix, AtmosphericLut& atmLut, OceanReflLut& ocnLut, double *dRsurf__dToa) :
    pix(pix), atmLut(atmLut), ocnLut(ocnLut), dRsurf_dToa(dRsurf__dToa) {
    }

    Doub operator() (Doub wdsp){
        int i, j;
        float fmin = 0.0;
        float m_error[N_SLSTR_BANDS * N_SLSTR_VIEWS];
        double diff, o_error;

        if (wdsp < ocnLut.wdspD.min) wdsp = ocnLut.wdspD.min;
        if (wdsp > ocnLut.wdspD.max) wdsp = ocnLut.wdspD.max;
        //pixel->ocn_wind_speed = wdsp;
        ocnLut.get_rho_ocean_sqr_error(m_error, &pix, pix.aod, ocn_d_wind_speed, ocn_d_pigment);
        //get_rho_ocean_sqr_error(m_error, pixel, 0.0, ocn_d_wind_speed, ocn_d_pigment);

        for (i = 1; i < N_SLSTR_BANDS; i++) {
            for (j = 0; j < N_SLSTR_VIEWS; j++) {
                if ( ! pix.view_clear[j] ) {
                    pix.angWeights[i][j] = 0;
                }
                else {

                    o_error = (0.5*inst_frac_error[i] * pix.tarr[i][j] * dRsurf_dToa[i*N_SLSTR_VIEWS+j]);
                    o_error *= o_error;
                    
                    o_error += rt_abs_error[i] * rt_abs_error[i];

                    //o_error += aer_mod_frac_error[i] * Rpath;
                    //o_error += pow ( (0.5 * aer_mod_frac_error[i] * Rpath), 2);
                    //o_error += ((1 - pixel->dust_weight) * aer_mod_frac_error[i] + pixel->dust_weight * aer_mod_frac_error_d[i])* Rpath;

                    /* Replace WG to give chi sq. Could still use old WG in addition if normalised to 1 */
                    pix.angWeights[i][j] = 1 / (m_error[i * N_SLSTR_VIEWS + j] + o_error);
                    if (i==1) pix.angWeights[i][j] *= 0.2;
                    
                    diff = pix.RR[i][j] - pix.rho_surf[i][j];
                    fmin = fmin + diff * diff * pix.angWeights[i][j];
                }
            }
        }
        if (!pix.view_clear[0] || !pix.view_clear[1]) fmin *= 2;
        return fmin;
    }    
};

class EmodTauUncert : public EmodTau {
public:
    SlstrPixel& pix;
    AtmosphericLut& atmLut;
    OceanReflLut& ocnLut;

    EmodTauUncert(SlstrPixel& pix, AtmosphericLut& atmLut, OceanReflLut& ocnLut) :
    pix(pix), atmLut(atmLut), ocnLut(ocnLut) {
    }

    Doub operator() (Doub tau){

        int i, j;
        Doub fret = 0.0, fmin = 0.0, y;
        Doub p[] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.5, 0.3};
        double dRsurf_dToa[N_SLSTR_BANDS * N_SLSTR_VIEWS];
        VecDoub pVec(7, &p[0]);

        if (tau < atmLut.aodD.min) tau = atmLut.aodD.min;
        if (tau > atmLut.aodD.max) tau = atmLut.aodD.max;
        pix.aod=tau;
        atmLut.psInv6s(&pix, tau);

        pix.ndvi = (pix.RR[2][0] - pix.RR[1][0]) / (pix.RR[2][0] + pix.RR[1][0]);
        //pix.ndvi = (pix.tarr[2][0] - pix.tarr[1][0]) / (pix.tarr[2][0] + pix.tarr[1][0]);
        pix.ndvi_veg_weight = (1.0 - 0.0) / (0.35 - 0.15) * (pix.ndvi - 0.15);
        if (pix.ndvi_veg_weight < 0) pix.ndvi_veg_weight = 0.0;
        if (pix.ndvi_veg_weight > 1) pix.ndvi_veg_weight = 1.0;
        pix.dust_weight = (0.0 - 1.0) / (0.6 - 0.4) * (pix.tarr[3][0] - 0.6);
        if (pix.dust_weight < 0) pix.dust_weight = 0.0;
        if (pix.dust_weight > 1) pix.dust_weight = 1.0;

        for (i = 0; i < N_SLSTR_BANDS; i++){
            for (j = 0; j < N_SLSTR_VIEWS; j++){
                if (pix.RR[i][j] < 0.001) { //maybe only for i>0
                    fmin = fmin + pow((pix.RR[i][j] - 0.001), 2) * 1000000.0;
                }
                y = (pix.tarr[i][j] + dToa - pix.inv_coef[i][j].rPath) / pix.inv_coef[i][j].tGas;
                y /= (pix.inv_coef[i][j].tDown * pix.inv_coef[i][j].tUp);
                dRsurf_dToa[i*N_SLSTR_VIEWS + j] = fabs((y / (1.0 + pix.inv_coef[i][j].spherAlb * y)) - pix.RR[i][j]) / dToa;
            }
        }
        setBit(&pix.qflag, EMOD_PENALTY, (fmin > 0));

        if (fmin < 10){
            EmodUncert emodUncert(pix, atmLut, ocnLut, dRsurf_dToa);
            Powell<EmodUncert> powell(emodUncert);
            pVec = powell.minimize(pVec);
            fret = powell.fret;
        }
        if (pix.dumpRR){
            //dumpRR(pixel, p, tau);
        }
        return fmin + fret;
    }
};

class EmodTauOcean : public EmodTau {
public:
    SlstrPixel& pix;
    AtmosphericLut& atmLut;
    OceanReflLut& ocnLut;

    EmodTauOcean(SlstrPixel& pix, AtmosphericLut& atmLut, OceanReflLut& ocnLut) :
    pix(pix), atmLut(atmLut), ocnLut(ocnLut) {
    }

    Doub operator() (Doub tau){
        int i, j;
        float y, fmin = 0.0;
        double dRsurf_dToa[N_SLSTR_BANDS*N_SLSTR_VIEWS];
        //float ax = WDSP_MIN;
        //float bx = 6;
        //float cx = WDSP_MAX;

        if (tau < atmLut.aodD.min) tau = atmLut.aodD.min;
        if (tau > atmLut.aodD.max) tau = atmLut.aodD.max;
        pix.aod=tau;
        atmLut.psInv6s(&pix, tau);

        /* penalise -ve reflectances */
        for (i = 0; i < N_SLSTR_BANDS; i++){
            for (j = 0; j < N_SLSTR_VIEWS; j++){
                if (pix.RR[i][j] < -1e-6){
                    fmin = fmin + (pix.RR[i][j] * pix.RR[i][j]) * 10000.0; //1000000.0;
                }
                y = (pix.tarr[i][j] + dToa - pix.inv_coef[i][j].rPath) / pix.inv_coef[i][j].tGas;
                y /= (pix.inv_coef[i][j].tDown * pix.inv_coef[i][j].tUp);
                dRsurf_dToa[i*N_SLSTR_VIEWS + j] = fabs((y / (1.0 + pix.inv_coef[i][j].spherAlb * y)) - pix.RR[i][j]) / dToa;
            }
        }
        setBit(&pix.qflag, EMOD_PENALTY, (fmin > 0));

        //if (fmin < 100.0){
            EmodOcean emodOcean(pix, atmLut, ocnLut, dRsurf_dToa);
            //fmin = fmin + brent_d(&ax, bx, &cx, emod_ocean, TOL, &pixel->ocn_wind_speed);
            fmin = fmin + emodOcean(pix.ocn_wind_speed);
        //}
        return fmin;
    }
};



class AodRetrieval {
public:
    SlstrPixel& pix;
    AtmosphericLut& atmLut;
    OceanReflLut& ocnLut;


    AodRetrieval(SlstrPixel& slstrPixel, AtmosphericLut& aLut, OceanReflLut& oLut);
    ~AodRetrieval();

    void retrieveAodSizeBrent(bool isOcean);

private:
    AodRetrieval();
    AodRetrieval(const AodRetrieval& orig) : pix(orig.pix), atmLut(orig.atmLut), ocnLut(orig.ocnLut) {}; // disable copy construtor
    AodRetrieval& operator=(const AodRetrieval& rhs) { throw std::logic_error("AodRetrieval shoudlnt be copied"); } // disable copy assignment

    float getCurvature(EmodTau *emodTau, char isOcean);
/*
    float emod_size(float fot);
    float emod_tau_ocean_model(float tau);
    float emod_tau_uncer(float tau);
    float emod_uncer(float *);
*/
};



#endif /* AODRETRIEVAL_HPP */

