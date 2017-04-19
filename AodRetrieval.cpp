/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   AodRetrieval.cpp
 * Author: akheckel
 * 
 * Created on 8. April 2017, 14:25
 */

#include <cmath>
#include <stdexcept>
#include "miscUtils.hpp"
#include "AodRetrieval.hpp"



//
// public
//

AodRetrieval::AodRetrieval(SlstrPixel& slstrPixel, AtmosphericLut& aLut, OceanReflLut& oLut)
 : pix(slstrPixel), atmLut(aLut), ocnLut(oLut) {}

AodRetrieval::~AodRetrieval() {}

void AodRetrieval::retrieveAodSizeBrent(bool isOcean){
    float ax = atmLut.aodD.min;
    float cx = atmLut.aodD.max;
    float bx = 0.05;
    float fot, cc;
    pix.lutpars.climFineMode = pix.lutpars.mix_frac[0];
    //float(*fct_emod_tau)(float) = (isOcean) ? emod_tau_ocean_model : emod_tau_uncer;

    // start retrrieval with fine mode fraction from climatology
    // to test amount of aod and determine ax cx limits
    //pix.fmin = brent_d(&ax, bx, &cx, fct_emod_tau, TOL, &pix.aod);
    if ((!isOcean) && (pix.aod > 0.01)){
        pix.ax = atmLut.aodD.min;  //ax - 0.5 * (ax - TAU_MIN);
        pix.cx = cx + 0.02 * (atmLut.aodD.max - cx);
    }
    else {
        pix.ax = atmLut.aodD.min;  //ax - 0.5 * (ax - TAU_MIN);
        pix.cx = atmLut.aodD.max;  //cx + 0.01 * (TAU_MAX - cx);
    }
    EmodTauUncert emod_tau_uncert(pix, atmLut, ocnLut);
    EmodTauOcean  emod_tau_ocean(pix, atmLut, ocnLut);
    EmodSize emod_size(pix, atmLut, ocnLut, &emod_tau_uncert);
    if (isOcean) {
        emod_size.emodTau = &emod_tau_ocean;
    }
/***/
    // start size fitting
    Brent brent;
    brent.ax = 0; brent.cx = 1;
    brent.bx = pix.lutpars.mix_frac[0];
    fot = pix.lutpars.mix_frac[0];
    brent.minimize(emod_size);
    pix.fmin = brent.fmin;
    fot = brent.xmin;

    // finally set pixel values to optimal values
    pix.lutpars.mix_frac[0] = fot;
    pix.lutpars.mixing[0] = ((1 - fot) * pix.lutpars.mix_frac[2])*100; // Dust
    pix.lutpars.mixing[1] = ((1 - fot) * (1 - pix.lutpars.mix_frac[2]))*100; // Sea Salt
    pix.lutpars.mixing[2] = (fot * (1 - pix.lutpars.mix_frac[1]))*100; // Strong Abs
    pix.lutpars.mixing[3] = (fot * pix.lutpars.mix_frac[1])*100; // Weak Abs
    atmLut.getTetrahedronPoints(&pix.lutpars, false);
    pix.lutpars.ocn_mi = ocnLut.getInterPar(fot, ocnLut.modelD); //(fot - MVM_MIN) / MVM_INT;
/***/
    // determine curvature to estimate uncertainty
    // rescaling uncertainty by about a factor of 0.25
    cc = getCurvature(emod_size.emodTau, isOcean);
    if ( cc > 0 ){
        pix.ediff = 1.0 / sqrt(2.0 * cc);
        
        if (!isOcean){
            pix.ediff *= 0.7; // factor 0.25 introduced to temp fix uncertainty overestimation
            if (pix.ediff < 0.02){
                pix.ediff = 0.02 + 0.05 * pix.aod;
            }
        }
    } else {
        setBit(&pix.qflag, CURV_NEGATIVE, (cc <= 0));
        pix.ediff = 0.02 + 0.05 * pix.aod;
    }

}


//
// private
//

    float AodRetrieval::getCurvature(EmodTau *emodTau, char isOcean){
        //float x[] = {-0.03+pixel->aod, pixel->aod, 0.03+pixel->aod};
        //float x[] = {0.75*pixel->aod, pixel->aod, 0.75*pixel->aod+0.5};
        double x[] = {0.7*pix.aod, pix.aod, 0.85*pix.aod};
        double x2[3], y[3];
        double r,cc;
        char i;

        //TODO: fixed constants need adjustment if LUT changes
        /*if (pixel->aod > 1.5) {
            x[2] = 0.825*pixel->aod;
        }
        else*/
        if (pix.aod < 0.05){
            x[0] = 0.002;
        }
        /*if (isOcean) {
            //x[0] = 0.9*pixel->aod;
            x[2] = 0.5*pixel->aod;
        }*/
        for (i=0; i<3; i++){
            x2[i] = x[i] * x[i];
        }
        y[0] = (*emodTau)(x[0]);
        y[2] = (*emodTau)(x[2]);
        if (pix.RR[0][0]<0){
            x[2] = 0.5*pix.aod;
            x2[2] = x[2] * x[2];
            y[2] = (*emodTau)(x[2]);
        }
        y[1] = (*emodTau)(x[1]); //needs to be done last to get correct RR for pixel
        r = (x[2] - x[0]) / (x[1] - x[0]);
        cc = y[2] - y[0] - r * (y[1] - y[0]);
        cc /= (x2[2] - x2[0] - r* (x2[1] - x2[0]));

        return cc;
    }

