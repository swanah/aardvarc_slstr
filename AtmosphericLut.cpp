/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   AtmosphericLut.cpp
 * Author: akheckel
 * 
 * Created on 30. März 2017, 15:49
 */

#include <cmath>
#include <netcdf>
#include <stdexcept>
#include "AtmosphericLut.hpp"

using namespace netCDF;

//
//public
//

AtmosphericLut::AtmosphericLut(const std::string& lutName) {
    NcFile ncF(lutName, NcFile::read);
    readDimVarInt(&ncF, "model", modelD);
    readDimVarDouble(&ncF, "band", wvlD);
    readDimVarDouble(&ncF, "tau", aodD);
    readDimVarDouble(&ncF, "pressure", presD);
    readDimVarDouble(&ncF, "raz", razD);
    readDimVarDouble(&ncF, "vza", vzaD);
    readDimVarDouble(&ncF, "sza", szaD);
    readLutVar(&ncF, "rPath", &rPath, &nRPath);
    readLutVar(&ncF, "T", &t, &nT);
    readLutVar(&ncF, "tGas", &tGas, &nTgas);
    readLutVar(&ncF, "spherAlb", &spherA, &nSpherA);
    readLutVar(&ncF, "D", &difFrac, &nDifFrac);
    readLutVar(&ncF, "spec_aod_ratio", &specAodRatio, &nSpecAodRatio);
    readLutVar(&ncF, "ssa", &ssa, &nSsa);    
    ncF.close();
    initIndexMap();
}

AtmosphericLut::~AtmosphericLut() {
    delete [] rPath;
    rPath = NULL;
    delete [] t;
    t = NULL;
    delete [] spherA;
    spherA = NULL;
    delete [] difFrac;
    difFrac = NULL;
    delete [] specAodRatio;
    specAodRatio = NULL;
    delete [] ssa;
    ssa = NULL;
}

char AtmosphericLut::getTetrahedronPoints(LutPars* lutpars, bool verbose){
    double val_interp = 0.0;
    float fPercent[4];
    int vA = 0, vB = 1, vC = 2, vD = 3, vP = -1;
    int val_a = 1, val_b = 1, val_c = 1, val_d = 1, noRoundingErrors;
    int iT, jT, j, k, l, m, begin, end, iterationsLeft = 5, pointEnclosed = 0, pointNum = 0;
    int botop[4][4]; //bottom, below, above, top
    int points[256][4];
    float score[2] = {1600.0, 0.0};
    int tetraBest[4][4];

    for (iT = 0; iT < 4; iT++) {
        fPercent[iT] = lutpars->mixing[iT];
        botop[iT][1] = CCI_INT * (int) (fPercent[iT] / CCI_INT); //below
        botop[iT][2] = botop[iT][1] + CCI_INT; //above
        botop[iT][0] = botop[iT][1] - CCI_INT; //bottom
        botop[iT][3] = botop[iT][2] + CCI_INT; //top
        if (botop[iT][0] < 0) botop[iT][0] = 0;
        if (botop[iT][3] > 100) botop[iT][3] = 100;
    }

    /***
     // Force the use of 4 pure components
     iterationsLeft=1;
     ***/

    //Try every combination from narrow to broad tetrahedron till enclosure
    while (((float) val_interp != 1.0) && (iterationsLeft > 0)) {
        pointNum = 0;
        switch (iterationsLeft) {
            case 5: // below and above  - creates the smallest Tetrahedron
                begin = 1;
                end = 3;
                break;
            case 4: // below and top
                begin = 1;
                end = 4;
                break;
            case 3: // bottom and above
                begin = 0;
                end = 3;
                break;
            case 2: // bottom and top
                begin = 0;
                end = 4;
                break;
            case 1: // 4 pure components  - creates the largest Tetrahedron
                begin = -1;
                end = 0;
                break;
        }
        if (begin >= 0) {
            for (j = begin; j < end; j++) {
                for (k = begin; k < end; k++) {
                    for (l = begin; l < end; l++) {
                        for (m = begin; m < end; m++) {
                            if ((botop[0][j] + botop[1][k] +
                                    botop[2][l] + botop[3][m]) == 100) {
                                points[pointNum][0] = botop[0][j];
                                points[pointNum][1] = botop[1][k];
                                points[pointNum][2] = botop[2][l];
                                points[pointNum][3] = botop[3][m];
                                pointNum++; // valid point, on to the next ...
                            }
                        }
                    }
                }
            }
        } else { // Use the 4 100% points
            pointNum = 4;
            for (jT = 0; jT < 4; jT++) {
                for (iT = 0; iT < 4; iT++) {
                    if (iT == jT)
                        points[jT][iT] = 100;
                    else
                        points[jT][iT] = 0;
                }
            }
        }
        //Cycle through all in this case to see if any will enclose the mix
        for (j = pointNum - 1; j >= 0; j--) {
            for (k = j - 1; k >= 0; k--) {
                for (l = k - 1; l >= 0; l--) {
                    for (m = l - 1; m >= 0; m--) {

                        for (iT = 0; iT < 4; iT++) {
                            lutpars->tetraP[0][iT] = points[j][iT];
                        }
                        for (iT = 0; iT < 4; iT++) {
                            lutpars->tetraP[1][iT] = points[k][iT];
                        }
                        for (iT = 0; iT < 4; iT++) {
                            lutpars->tetraP[2][iT] = points[l][iT];
                        }
                        for (iT = 0; iT < 4; iT++) {
                            lutpars->tetraP[3][iT] = points[m][iT];
                        }

                        if (j != k && j != l && j != m && k != l && k != m && l != m) {
                            //fprintf(stdout,"Trying %i %i %i %i\n", j,k,l,m);
                            val_interp = (getTetraVol(vA, vB, vC, vP, lutpars) * val_d + getTetraVol(vA, vB, vD, vP, lutpars) * val_c +
                                    getTetraVol(vA, vC, vD, vP, lutpars) * val_b + getTetraVol(vB, vC, vD, vP, lutpars) * val_a) /
                                    getTetraVol(vA, vB, vC, vD, lutpars);
                            //         noRoundingErrors = (int)((val_interp+0.000005) * 100000); //Round to 5 dec places
                            //         val_interp = (double)(noRoundingErrors/100000);
                            noRoundingErrors = (int) ((val_interp + 0.005) * 100); //Round to 2 dec places
                            val_interp = (double) (noRoundingErrors / 100);
                            if ((float) val_interp == 1.0) {
                                //To save processing time while sacrificing some accuracy,
                                //we could just use the 1st one we find
                                /***return(1);        //Comment this out for more accuracy ***/

                                //We have a tetrahedron that is enclosing the mix, but we want to
                                //find the one that is closest to it.
                                score[LATEST] = 0.0;
                                for (jT = 0; jT < 4; jT++) {
                                    for (iT = 0; iT < 4; iT++) {
                                        score[LATEST] = score[LATEST] + fabs((float) lutpars->tetraP[jT][iT] - fPercent[iT]);
                                    }
                                }
                                /*** fprintf(stdout,"Tetra Point 1: %i %i %i %i \n",tetraP[0][0],tetraP[0][1],tetraP[0][2],tetraP[0][3]);
                                     fprintf(stdout,"Tetra Point 2: %i %i %i %i \n",tetraP[1][0],tetraP[1][1],tetraP[1][2],tetraP[1][3]);
                                     fprintf(stdout,"Tetra Point 3: %i %i %i %i \n",tetraP[2][0],tetraP[2][1],tetraP[2][2],tetraP[2][3]);
                                     fprintf(stdout,"Tetra Point 4: %i %i %i %i \n",tetraP[3][0],tetraP[3][1],tetraP[3][2],tetraP[3][3]);
                                     fprintf(stdout,"Score: latest %f best %f \n",score[LATEST],score[BEST]);   ***/
                                if (score[LATEST] < score[BEST]) {
                                    score[BEST] = score[LATEST];
                                    for (jT = 0; jT < 4; jT++) {
                                        for (iT = 0; iT < 4; iT++) {
                                            tetraBest[jT][iT] = lutpars->tetraP[jT][iT];
                                        }
                                    }
                                }
                                pointEnclosed = 1; //We have enclosure on this case
                                iterationsLeft = 0; //No need to widen the tetrahedron
                            }
                        }
                    }
                }
            }
        }
        iterationsLeft--;
    }

    if (pointEnclosed == 1) {
        for (jT = 0; jT < 4; jT++) {
            for (iT = 0; iT < 4; iT++) {
                lutpars->tetraP[jT][iT] = tetraBest[jT][iT];
            }
        }

        /***
         //Swap the vertices to see if it makes a difference
         for(iT=0; iT<4; iT++) {
         tetraP[3][iT]=tetraBest[0][iT];
         }
         for(iT=0; iT<4; iT++) {
         tetraP[2][iT]=tetraBest[1][iT];
         }
         for(iT=0; iT<4; iT++) {
         tetraP[1][iT]=tetraBest[2][iT];
         }
         for(iT=0; iT<4; iT++) {
         tetraP[0][iT]=tetraBest[3][iT];
         }
         ***/
        if (verbose) {
            fprintf(stdout, "Tetra Point 1: %i %i %i %i \n", lutpars->tetraP[0][0], lutpars->tetraP[0][1], lutpars->tetraP[0][2], lutpars->tetraP[0][3]);
            fprintf(stdout, "Tetra Point 2: %i %i %i %i \n", lutpars->tetraP[1][0], lutpars->tetraP[1][1], lutpars->tetraP[1][2], lutpars->tetraP[1][3]);
            fprintf(stdout, "Tetra Point 3: %i %i %i %i \n", lutpars->tetraP[2][0], lutpars->tetraP[2][1], lutpars->tetraP[2][2], lutpars->tetraP[2][3]);
            fprintf(stdout, "Tetra Point 4: %i %i %i %i \n", lutpars->tetraP[3][0], lutpars->tetraP[3][1], lutpars->tetraP[3][2], lutpars->tetraP[3][3]);
        }
        return (1);
    } else
        return (0);
}

double AtmosphericLut::getTetraVol(int v0, int v1, int v2, int vOrigin, LutPars* lutpars) {
    float Xa, Xb, Xc, Ya, Yb, Yc, Za, Zb, Zc;
    double volume;

    if (vOrigin < 0) {
        Xa = (float) lutpars->tetraP[v0][0] - lutpars->mixing[0];
        Ya = (float) lutpars->tetraP[v0][1] - lutpars->mixing[1];
        Za = (float) lutpars->tetraP[v0][2] - lutpars->mixing[2];
        Xb = (float) lutpars->tetraP[v1][0] - lutpars->mixing[0];
        Yb = (float) lutpars->tetraP[v1][1] - lutpars->mixing[1];
        Zb = (float) lutpars->tetraP[v1][2] - lutpars->mixing[2];
        Xc = (float) lutpars->tetraP[v2][0] - lutpars->mixing[0];
        Yc = (float) lutpars->tetraP[v2][1] - lutpars->mixing[1];
        Zc = (float) lutpars->tetraP[v2][2] - lutpars->mixing[2];
    } else {
        Xa = ((float) lutpars->tetraP[v0][0] - lutpars->tetraP[vOrigin][0]);
        Ya = ((float) lutpars->tetraP[v0][1] - lutpars->tetraP[vOrigin][1]);
        Za = ((float) lutpars->tetraP[v0][2] - lutpars->tetraP[vOrigin][2]);
        Xb = ((float) lutpars->tetraP[v1][0] - lutpars->tetraP[vOrigin][0]);
        Yb = ((float) lutpars->tetraP[v1][1] - lutpars->tetraP[vOrigin][1]);
        Zb = ((float) lutpars->tetraP[v1][2] - lutpars->tetraP[vOrigin][2]);
        Xc = ((float) lutpars->tetraP[v2][0] - lutpars->tetraP[vOrigin][0]);
        Yc = ((float) lutpars->tetraP[v2][1] - lutpars->tetraP[vOrigin][1]);
        Zc = ((float) lutpars->tetraP[v2][2] - lutpars->tetraP[vOrigin][2]);
    }

    /***
        if (test==1) {
        double yz, zy, xz, zx, xy, yx, sub1, sub2, sub3, p1, p2, p3;
        yz=(Yb*Zc);
        zy=(Zb*Yc);
        xz=(Xb*Zc);
        zx=(Zb*Xc);
        xy=(Xb*Yc);
        yx=(Yb*Xc);
        sub1=yz-zy;
        sub2=xz-zx;
        sub3=xy-yx;
        p1=Xa*sub1;
        p2=Ya*sub2;
        p3=Za*sub3;
        }
     ***/

    volume = fabs(((Xa * ((Yb * Zc)-(Zb * Yc)))-(Ya * ((Xb * Zc)-(Zb * Xc)))+(Za * ((Xb * Yc)-(Yb * Xc)))) / 6);

    return (volume);
}

void AtmosphericLut::psInv6s(SlstrPixel* p, float tau){
    
    LutCoefs lutcoef[N_SLSTR_VIEWS];
    LutCoefs coeffs[N_SLSTR_VIEWS];
    float spec_aod;
    float ssa_t;
    float taui;
    int band, angle;
    int nsDi, seaSi, strongAi, iAero;
    int vA = 0, vB = 1, vC = 2, vD = 3, vP = -1;
    double vol, yval;
    taui = getInterPar(tau, aodD);//(tau - aodD.min) / aodD.delta;
    
    for (band = 0; band < N_SLSTR_BANDS; band++) {
        for (angle = 0; angle < N_SLSTR_VIEWS; angle++) {
            coeffs[angle].rPath = coeffs[angle].difFrac = coeffs[angle].spherAlb = 0;
            coeffs[angle].tDown = coeffs[angle].tUp = coeffs[angle].tGas = 0;
        }
        spec_aod = 0;
        ssa_t = 0;
        
        nsDi = p->lutpars.tetraP[0][0] / CCI_INT;
        seaSi = p->lutpars.tetraP[0][1] / CCI_INT;
        strongAi = p->lutpars.tetraP[0][2] / CCI_INT;
        iAero = aerModelIdxMap[nsDi*100 + seaSi*10 + strongAi];
        getCoefs(lutcoef, &p->lutpars, taui, band, iAero);
        vol = getTetraVol(vB, vC, vD, vP, &p->lutpars);
        for (angle = 0; angle < N_SLSTR_VIEWS; angle++){
            coeffs[angle].rPath    += vol * lutcoef[angle].rPath;
            coeffs[angle].tDown    += vol * lutcoef[angle].tDown;
            coeffs[angle].tUp      += vol * lutcoef[angle].tUp;
            coeffs[angle].tGas     += vol * lutcoef[angle].tGas;
            coeffs[angle].spherAlb += vol * lutcoef[angle].spherAlb;
            coeffs[angle].difFrac  += vol * lutcoef[angle].difFrac;
        }
        spec_aod += vol * specAodRatio[iAero + modelD.n * band];
        ssa_t += vol * ssa[iAero + modelD.n * band];
        
        nsDi = p->lutpars.tetraP[1][0] / CCI_INT;
        seaSi = p->lutpars.tetraP[1][1] / CCI_INT;
        strongAi = p->lutpars.tetraP[1][2] / CCI_INT;
        iAero = aerModelIdxMap[nsDi*100 + seaSi*10 + strongAi];
        getCoefs(lutcoef, &p->lutpars, taui, band, iAero);
        vol = getTetraVol(vA, vC, vD, vP, &p->lutpars);
        for (angle = 0; angle < N_SLSTR_VIEWS; angle++){
            coeffs[angle].rPath    += vol * lutcoef[angle].rPath;
            coeffs[angle].tDown    += vol * lutcoef[angle].tDown;
            coeffs[angle].tUp      += vol * lutcoef[angle].tUp;
            coeffs[angle].tGas     += vol * lutcoef[angle].tGas;
            coeffs[angle].spherAlb += vol * lutcoef[angle].spherAlb;
            coeffs[angle].difFrac  += vol * lutcoef[angle].difFrac;
        }
        spec_aod += vol * specAodRatio[iAero + modelD.n * band];
        ssa_t += vol * ssa[iAero + modelD.n * band];
        
        nsDi = p->lutpars.tetraP[2][0] / CCI_INT;
        seaSi = p->lutpars.tetraP[2][1] / CCI_INT;
        strongAi = p->lutpars.tetraP[2][2] / CCI_INT;
        iAero = aerModelIdxMap[nsDi*100 + seaSi*10 + strongAi];
        getCoefs(lutcoef, &p->lutpars, taui, band, iAero);
        vol = getTetraVol(vA, vB, vD, vP, &p->lutpars);
        for (angle = 0; angle < N_SLSTR_VIEWS; angle++){
            coeffs[angle].rPath    += vol * lutcoef[angle].rPath;
            coeffs[angle].tDown    += vol * lutcoef[angle].tDown;
            coeffs[angle].tUp      += vol * lutcoef[angle].tUp;
            coeffs[angle].tGas     += vol * lutcoef[angle].tGas;
            coeffs[angle].spherAlb += vol * lutcoef[angle].spherAlb;
            coeffs[angle].difFrac  += vol * lutcoef[angle].difFrac;
        }
        spec_aod += vol * specAodRatio[iAero + modelD.n * band];
        ssa_t += vol * ssa[iAero + modelD.n * band];
        
        nsDi = p->lutpars.tetraP[3][0] / CCI_INT;
        seaSi = p->lutpars.tetraP[3][1] / CCI_INT;
        strongAi = p->lutpars.tetraP[3][2] / CCI_INT;
        iAero = aerModelIdxMap[nsDi*100 + seaSi*10 + strongAi];
        getCoefs(lutcoef, &p->lutpars, taui, band, iAero);
        vol = getTetraVol(vA, vB, vC, vP, &p->lutpars);
        for (angle = 0; angle < N_SLSTR_VIEWS; angle++){
            coeffs[angle].rPath    += vol * lutcoef[angle].rPath;
            coeffs[angle].tDown    += vol * lutcoef[angle].tDown;
            coeffs[angle].tUp      += vol * lutcoef[angle].tUp;
            coeffs[angle].tGas     += vol * lutcoef[angle].tGas;
            coeffs[angle].spherAlb += vol * lutcoef[angle].spherAlb;
            coeffs[angle].difFrac  += vol * lutcoef[angle].difFrac;
        }
        spec_aod += vol * specAodRatio[iAero + modelD.n * band];
        ssa_t += vol * ssa[iAero + modelD.n * band];
        
        vol = getTetraVol(vA, vB, vC, vD, &p->lutpars);
        
        for (angle = 0; angle < N_SLSTR_VIEWS; angle++){
            coeffs[angle].rPath    /= vol;
            coeffs[angle].tDown    /= vol;
            coeffs[angle].tUp      /= vol;
            coeffs[angle].tGas     /= vol;
            coeffs[angle].spherAlb /= vol;
            coeffs[angle].difFrac  /= vol;
            
            if (p->tarr[band][angle] > 1e-6) {
                yval = (p->tarr[band][angle] - coeffs[angle].rPath) / coeffs[angle].tGas;
                yval /= (coeffs[angle].tDown * coeffs[angle].tUp);
                p->RR[band][angle] = yval / (1.0 + coeffs[angle].spherAlb * yval);
                        
                /*yval = coeffs[angle].cfac * p->tarr[band][angle] * coeffs[angle].xa - coeffs[angle].xb;
                p->RR[band][angle] = yval / (1.0 + coeffs[angle].xc * yval);*/
            }
            else {
                p->RR[band][angle] = 0;
            }
            p->DD[band][angle] = coeffs[angle].difFrac;
            
            p->inv_coef[band][angle].rPath = coeffs[angle].rPath;
            p->inv_coef[band][angle].tDown = coeffs[angle].tDown;
            p->inv_coef[band][angle].tUp = coeffs[angle].tUp;
            p->inv_coef[band][angle].tGas = coeffs[angle].tGas;
            p->inv_coef[band][angle].spherAlb = coeffs[angle].spherAlb;
            p->inv_coef[band][angle].difFrac = coeffs[angle].difFrac;
        }
        p->spec_aod_fac[band] = spec_aod / vol;
        p->ssa[band] = ssa_t / vol;
        
    }
}

/**
 * compute fractional index, includes range check and verifies the returned index is always smaller than dimPar.n-1
 * the last verification avoids testing in the lut interpolation!!! 
 * @param val - actual value along lut dim
 * @param dimPar - corresponding dim parameter
 * @return fractional index in dim axis
 */
float AtmosphericLut::getInterPar(const float val, const LutDimParamDouble& dimPar){
    if (val < dimPar.min || val > dimPar.max) throw std::range_error("variable not in range of LUT");
    return (val == dimPar.max) ? (val * 0.9999 - dimPar.min) / dimPar.delta: (val - dimPar.min) / dimPar.delta;
}

//private
    
void AtmosphericLut::readDimVarByte(NcFile *ncF, std::string varName, LutDimParamByte& dimPar){
    NcVar v;
    v = ncF->getVar(varName);
    dimPar.n = v.getDim(0).getSize();
    char a[dimPar.n];
    v.getVar(&a);
    dimPar.min = a[0];
    dimPar.max = a[dimPar.n-1];
    dimPar.delta = a[1] - a[0]; 
}

void AtmosphericLut::readDimVarInt(NcFile *ncF, std::string varName, LutDimParamInt& dimPar){
    NcVar v;
    v = ncF->getVar(varName);
    dimPar.n = v.getDim(0).getSize();
    int a[dimPar.n];
    v.getVar(&a);
    dimPar.min = a[0];
    dimPar.max = a[dimPar.n-1];
    dimPar.delta = a[1] - a[0]; 
}

void AtmosphericLut::readDimVarFloat(NcFile *ncF, std::string varName, LutDimParamFloat& dimPar){
    NcVar v;
    v = ncF->getVar(varName);
    dimPar.n = v.getDim(0).getSize();
    float a[dimPar.n];
    v.getVar(&a);
    dimPar.min = a[0];
    dimPar.max = a[dimPar.n-1];
    dimPar.delta = a[1] - a[0]; 
}

void AtmosphericLut::readDimVarDouble(NcFile *ncF, std::string varName, LutDimParamDouble& dimPar){
    NcVar v;
    v = ncF->getVar(varName);
    dimPar.n = v.getDim(0).getSize();
    double a[dimPar.n];
    v.getVar(&a);
    dimPar.min = a[0];
    dimPar.max = a[dimPar.n-1];
    dimPar.delta = a[1] - a[0]; 
}

void AtmosphericLut::readLutVar(NcFile *ncF, std::string varName, double** data, size_t* nData){
    NcVar v;
    v = ncF->getVar(varName);
    *nData = 1;
    for (int iDim = v.getDimCount()-1; iDim >= 0; iDim--){
        *nData *= v.getDim(iDim).getSize();
    }
    *data = new double[*nData];
    v.getVar(*data);
}

/**
 * Generate a mapping of aerosol mixing indices (iNsD, iSeaS, iSAbs) to aerosol model index (0..34)
 * @param idxMap
 */
void AtmosphericLut::initIndexMap(){
    int iNsD, iSeaS, iStrongA, iWeakA;
    int mix1, mix2, i=0;
    for (iNsD = 0; iNsD < CCI_NUM; iNsD = iNsD + 1) {
        mix1 = CCI_NUM - iNsD; // The mix that we have left to play with
        for (iSeaS = 0; iSeaS < mix1; iSeaS = iSeaS + 1) {
            mix2 = CCI_NUM - iNsD - iSeaS;
            for (iStrongA = 0; iStrongA < mix2; iStrongA = iStrongA + 1) {
                //iWeakA = (CCI_NUM-1) - iNsD - iSeaS - iStrongA;
                aerModelIdxMap[iNsD*100 + iSeaS*10 + iStrongA] = i;
                i++;
            }
        }
    }
}

void AtmosphericLut::getCoefs(LutCoefs* coef, const LutPars* lpars, const float taui, const int iBand, const int iModel) {
    
    coef[0].rPath = interpolRpath(lpars->szani, lpars->vzani, lpars->razni, lpars->pAlti, taui, iBand, iModel);
    coef[0].tDown = interpolTscat(lpars->szani, lpars->pAlti, taui, iBand, iModel);
    coef[0].tUp   = interpolTscat(lpars->vzani, lpars->pAlti, taui, iBand, iModel);
    coef[0].tGas  = interpolTgas(lpars->szani, lpars->vzani, lpars->pAlti, iBand, iModel);
    coef[0].spherAlb = interpolSpherA(lpars->pAlti, taui, iBand, iModel);
    coef[0].difFrac  = interpolDifFrac(lpars->szani, lpars->pAlti, taui, iBand, iModel);

    coef[1].rPath = interpolRpath(lpars->szafi, lpars->vzafi, lpars->razfi, lpars->pAlti, taui, iBand, iModel);
    coef[1].tDown = interpolTscat(lpars->szafi, lpars->pAlti, taui, iBand, iModel);
    coef[1].tUp   = interpolTscat(lpars->vzafi, lpars->pAlti, taui, iBand, iModel);
    coef[1].tGas  = interpolTgas(lpars->szafi, lpars->vzafi, lpars->pAlti, iBand, iModel);
    coef[1].spherAlb = interpolSpherA(lpars->pAlti, taui, iBand, iModel);
    coef[1].difFrac  = interpolDifFrac(lpars->szafi, lpars->pAlti, taui, iBand, iModel);
}


double AtmosphericLut::interpolSpherA(const float presi, const float taui, const int iBand, const int iModel){
    double v[2], w[2];
    int x1i, x2i;
    int bIdx, off[5];
    double val;

    x1i = (int) presi; 
    x2i = (int) taui; 
    
    v[0] = presi - x1i; v[1] = 1 - v[0];
    w[0] = taui  - x2i; w[1] = 1 - w[0];
    
    off[0] = modelD.n * wvlD.n;
    off[1] = off[0] * aodD.n;
    bIdx = iModel + modelD.n * iBand 
           + off[0] * x2i
           + off[1] * x1i;

    val  = v[1] * w[1] * spherA[bIdx                  ]; // x1i,x2i = 0,0 
    val += v[1] * w[0] * spherA[bIdx + off[0]         ]; // x1i,x2i = 0,1 
    val += v[0] * w[1] * spherA[bIdx          + off[1]]; // x1i,x2i = 1,0 
    val += v[0] * w[0] * spherA[bIdx + off[0] + off[1]]; // x1i,x2i = 1,1 
    
    return val;
}

double AtmosphericLut::interpolDifFrac(const float szai, const float presi, const float taui, 
                                   const int iBand, const int iModel){
    double u[2], v[2], w[2];
    int x1i, x2i, x3i;
    int bIdx, off[5];
    double val;

    x1i = (int) szai; 
    x2i = (int) presi; 
    x3i = (int) taui; 
    
    u[0] = szai  - x1i; u[1] = 1 - u[0];
    v[0] = presi - x2i; v[1] = 1 - v[0];
    w[0] = taui  - x3i; w[1] = 1 - w[0];
    
    off[0] = modelD.n * wvlD.n;
    off[1] = off[0] * aodD.n;
    off[2] = off[1] * presD.n;
    bIdx = iModel + modelD.n * iBand 
           + off[0] * x3i
           + off[1] * x2i
           + off[2] * x1i;

    val  = u[1] * v[1] * w[1] * difFrac[bIdx                           ]; // x1i,x2i,x3i = 0,0,0 
    val += u[1] * v[1] * w[0] * difFrac[bIdx + off[0]                  ]; // x1i,x2i,x3i = 0,0,1 
    val += u[1] * v[0] * w[1] * difFrac[bIdx          + off[1]         ]; // x1i,x2i,x3i = 0,1,0 
    val += u[1] * v[0] * w[0] * difFrac[bIdx + off[0] + off[1]         ]; // x1i,x2i,x3i = 0,1,1 
    val += u[0] * v[1] * w[1] * difFrac[bIdx                   + off[2]]; // x1i,x2i,x3i = 1,0,0 
    val += u[0] * v[1] * w[0] * difFrac[bIdx + off[0]          + off[2]]; // x1i,x2i,x3i = 1,0,1 
    val += u[0] * v[0] * w[1] * difFrac[bIdx          + off[1] + off[2]]; // x1i,x2i,x3i = 1,1,0 
    val += u[0] * v[0] * w[0] * difFrac[bIdx + off[0] + off[1] + off[2]]; // x1i,x2i,x3i = 1,1,1 
    
    return val;
}

double AtmosphericLut::interpolTgas(const float szai, const float vzai, const float presi, 
                                   const int iBand, const int iModel){
    double u[2], v[2], w[2];
    int x1i, x2i, x3i;
    int bIdx, off[5];
    double val;

    x1i = (int) szai; 
    x2i = (int) vzai; 
    x3i = (int) presi; 
    
    u[0] = szai  - x1i; u[1] = 1 - u[0];
    v[0] = vzai  - x2i; v[1] = 1 - v[0];
    w[0] = presi - x3i; w[1] = 1 - w[0];
    
    // val  = s[1] * t[1] * u[1] * v[1] * w[1]
    //        * rPath[iModel + modelD.n * (iBand + wvlD.n * (x5i[0] + aodD.n * (x4i[0] + presD.n * (x3i[0] + razD.n * (x2i[0] + vzaD.n * x1i[0])))))]; 
    off[0] = modelD.n * wvlD.n;
    off[1] = off[0] * presD.n;
    off[2] = off[1] * vzaD.n;
    bIdx = iModel + modelD.n * iBand 
           + off[0] * x3i
           + off[1] * x2i
           + off[2] * x1i;

    val  = u[1] * v[1] * w[1] * tGas[bIdx                           ]; // x1i,x2i,x3i = 0,0,0 
    val += u[1] * v[1] * w[0] * tGas[bIdx + off[0]                  ]; // x1i,x2i,x3i = 0,0,1 
    val += u[1] * v[0] * w[1] * tGas[bIdx          + off[1]         ]; // x1i,x2i,x3i = 0,1,0 
    val += u[1] * v[0] * w[0] * tGas[bIdx + off[0] + off[1]         ]; // x1i,x2i,x3i = 0,1,1 
    val += u[0] * v[1] * w[1] * tGas[bIdx                   + off[2]]; // x1i,x2i,x3i = 1,0,0 
    val += u[0] * v[1] * w[0] * tGas[bIdx + off[0]          + off[2]]; // x1i,x2i,x3i = 1,0,1 
    val += u[0] * v[0] * w[1] * tGas[bIdx          + off[1] + off[2]]; // x1i,x2i,x3i = 1,1,0 
    val += u[0] * v[0] * w[0] * tGas[bIdx + off[0] + off[1] + off[2]]; // x1i,x2i,x3i = 1,1,1 
    
    return val;
}

double AtmosphericLut::interpolTscat(const float szai, const float presi, const float taui, 
                                   const int iBand, const int iModel){
    double u[2], v[2], w[2];
    int x1i, x2i, x3i;
    int bIdx, off[5];
    double val;

    x1i = (int) szai; 
    x2i = (int) presi; 
    x3i = (int) taui; 
    
    u[0] = szai  - x1i; u[1] = 1 - u[0];
    v[0] = presi - x2i; v[1] = 1 - v[0];
    w[0] = taui  - x3i; w[1] = 1 - w[0];
    
    // val  = s[1] * t[1] * u[1] * v[1] * w[1]
    //        * rPath[iModel + modelD.n * (iBand + wvlD.n * (x5i[0] + aodD.n * (x4i[0] + presD.n * (x3i[0] + razD.n * (x2i[0] + vzaD.n * x1i[0])))))]; 
    off[0] = modelD.n * wvlD.n;
    off[1] = off[0] * aodD.n;
    off[2] = off[1] * presD.n;
    bIdx = iModel + modelD.n * iBand 
           + off[0] * x3i
           + off[1] * x2i
           + off[2] * x1i;

    val  = u[1] * v[1] * w[1] * t[bIdx                           ]; // x1i,x2i,x3i = 0,0,0 
    val += u[1] * v[1] * w[0] * t[bIdx + off[0]                  ]; // x1i,x2i,x3i = 0,0,1 
    val += u[1] * v[0] * w[1] * t[bIdx          + off[1]         ]; // x1i,x2i,x3i = 0,1,0 
    val += u[1] * v[0] * w[0] * t[bIdx + off[0] + off[1]         ]; // x1i,x2i,x3i = 0,1,1 
    val += u[0] * v[1] * w[1] * t[bIdx                   + off[2]]; // x1i,x2i,x3i = 1,0,0 
    val += u[0] * v[1] * w[0] * t[bIdx + off[0]          + off[2]]; // x1i,x2i,x3i = 1,0,1 
    val += u[0] * v[0] * w[1] * t[bIdx          + off[1] + off[2]]; // x1i,x2i,x3i = 1,1,0 
    val += u[0] * v[0] * w[0] * t[bIdx + off[0] + off[1] + off[2]]; // x1i,x2i,x3i = 1,1,1 
    
    return val;
}

double AtmosphericLut::interpolRpath(const float szai, const float vzai, const float razi, const float presi, const float taui, 
                                     const int iBand, const int iModel){

    double s[2], t[2], u[2], v[2], w[2];
    int x1i, x2i, x3i, x4i, x5i;
    int bIdx, off[5];
    double val;

    x1i = (int) szai; 
    x2i = (int) vzai; 
    x3i = (int) razi; 
    x4i = (int) presi; 
    x5i = (int) taui; 
    
    s[0] = szai  - x1i; s[1] = 1 - s[0];
    t[0] = vzai  - x2i; t[1] = 1 - t[0];
    u[0] = razi  - x3i; u[1] = 1 - u[0];
    v[0] = presi - x4i; v[1] = 1 - v[0];
    w[0] = taui  - x5i; w[1] = 1 - w[0];
    
    // val  = s[1] * t[1] * u[1] * v[1] * w[1]
    //        * rPath[iModel + modelD.n * (iBand + wvlD.n * (x5i[0] + aodD.n * (x4i[0] + presD.n * (x3i[0] + razD.n * (x2i[0] + vzaD.n * x1i[0])))))]; 
    off[0] = modelD.n * wvlD.n;
    off[1] = off[0] * aodD.n;
    off[2] = off[1] * presD.n;
    off[3] = off[2] * razD.n;
    off[4] = off[3] * vzaD.n;
    bIdx = iModel + modelD.n * iBand 
           + off[0] * x5i
           + off[1] * x4i
           + off[2] * x3i
           + off[3] * x2i
           + off[4] * x1i;

    val  = s[1] * t[1] * u[1] * v[1] * w[1] * rPath[bIdx                                             ]; // x1i,x2i,x3i,x4i,x5i = 0,0,0,0,0 
    val += s[1] * t[1] * u[1] * v[1] * w[0] * rPath[bIdx + off[0]                                    ]; // x1i,x2i,x3i,x4i,x5i = 0,0,0,0,1 
    val += s[1] * t[1] * u[1] * v[0] * w[1] * rPath[bIdx          + off[1]                           ]; // x1i,x2i,x3i,x4i,x5i = 0,0,0,1,0 
    val += s[1] * t[1] * u[1] * v[0] * w[0] * rPath[bIdx + off[0] + off[1]                           ]; // x1i,x2i,x3i,x4i,x5i = 0,0,0,1,1 
    val += s[1] * t[1] * u[0] * v[1] * w[1] * rPath[bIdx                   + off[2]                  ]; // x1i,x2i,x3i,x4i,x5i = 0,0,1,0,0 
    val += s[1] * t[1] * u[0] * v[1] * w[0] * rPath[bIdx + off[0]          + off[2]                  ]; // x1i,x2i,x3i,x4i,x5i = 0,0,1,0,1 
    val += s[1] * t[1] * u[0] * v[0] * w[1] * rPath[bIdx          + off[1] + off[2]                  ]; // x1i,x2i,x3i,x4i,x5i = 0,0,1,1,0 
    val += s[1] * t[1] * u[0] * v[0] * w[0] * rPath[bIdx + off[0] + off[1] + off[2]                  ]; // x1i,x2i,x3i,x4i,x5i = 0,0,1,1,1 
    
    val += s[1] * t[0] * u[1] * v[1] * w[1] * rPath[bIdx                            + off[3]         ]; // x1i,x2i,x3i,x4i,x5i = 0,1,0,0,0 
    val += s[1] * t[0] * u[1] * v[1] * w[0] * rPath[bIdx + off[0]                   + off[3]         ]; // x1i,x2i,x3i,x4i,x5i = 0,1,0,0,1 
    val += s[1] * t[0] * u[1] * v[0] * w[1] * rPath[bIdx          + off[1]          + off[3]         ]; // x1i,x2i,x3i,x4i,x5i = 0,1,0,1,0 
    val += s[1] * t[0] * u[1] * v[0] * w[0] * rPath[bIdx + off[0] + off[1]          + off[3]         ]; // x1i,x2i,x3i,x4i,x5i = 0,1,0,1,1 
    val += s[1] * t[0] * u[0] * v[1] * w[1] * rPath[bIdx                   + off[2] + off[3]         ]; // x1i,x2i,x3i,x4i,x5i = 0,1,1,0,0 
    val += s[1] * t[0] * u[0] * v[1] * w[0] * rPath[bIdx + off[0]          + off[2] + off[3]         ]; // x1i,x2i,x3i,x4i,x5i = 0,1,1,0,1 
    val += s[1] * t[0] * u[0] * v[0] * w[1] * rPath[bIdx          + off[1] + off[2] + off[3]         ]; // x1i,x2i,x3i,x4i,x5i = 0,1,1,1,0 
    val += s[1] * t[0] * u[0] * v[0] * w[0] * rPath[bIdx + off[0] + off[1] + off[2] + off[3]         ]; // x1i,x2i,x3i,x4i,x5i = 0,1,1,1,1 

    val += s[0] * t[1] * u[1] * v[1] * w[1] * rPath[bIdx                                     + off[4]]; // x1i,x2i,x3i,x4i,x5i = 1,0,0,0,0 
    val += s[0] * t[1] * u[1] * v[1] * w[0] * rPath[bIdx + off[0]                            + off[4]]; // x1i,x2i,x3i,x4i,x5i = 1,0,0,0,1 
    val += s[0] * t[1] * u[1] * v[0] * w[1] * rPath[bIdx          + off[1]                   + off[4]]; // x1i,x2i,x3i,x4i,x5i = 1,0,0,1,0 
    val += s[0] * t[1] * u[1] * v[0] * w[0] * rPath[bIdx + off[0] + off[1]                   + off[4]]; // x1i,x2i,x3i,x4i,x5i = 1,0,0,1,1 
    val += s[0] * t[1] * u[0] * v[1] * w[1] * rPath[bIdx                   + off[2]          + off[4]]; // x1i,x2i,x3i,x4i,x5i = 1,0,1,0,0 
    val += s[0] * t[1] * u[0] * v[1] * w[0] * rPath[bIdx + off[0]          + off[2]          + off[4]]; // x1i,x2i,x3i,x4i,x5i = 1,0,1,0,1 
    val += s[0] * t[1] * u[0] * v[0] * w[1] * rPath[bIdx          + off[1] + off[2]          + off[4]]; // x1i,x2i,x3i,x4i,x5i = 1,0,1,1,0 
    val += s[0] * t[1] * u[0] * v[0] * w[0] * rPath[bIdx + off[0] + off[1] + off[2]          + off[4]]; // x1i,x2i,x3i,x4i,x5i = 1,0,1,1,1 
    
    val += s[0] * t[0] * u[1] * v[1] * w[1] * rPath[bIdx                            + off[3] + off[4]]; // x1i,x2i,x3i,x4i,x5i = 1,1,0,0,0 
    val += s[0] * t[0] * u[1] * v[1] * w[0] * rPath[bIdx + off[0]                   + off[3] + off[4]]; // x1i,x2i,x3i,x4i,x5i = 1,1,0,0,1 
    val += s[0] * t[0] * u[1] * v[0] * w[1] * rPath[bIdx          + off[1]          + off[3] + off[4]]; // x1i,x2i,x3i,x4i,x5i = 1,1,0,1,0 
    val += s[0] * t[0] * u[1] * v[0] * w[0] * rPath[bIdx + off[0] + off[1]          + off[3] + off[4]]; // x1i,x2i,x3i,x4i,x5i = 1,1,0,1,1 
    val += s[0] * t[0] * u[0] * v[1] * w[1] * rPath[bIdx                   + off[2] + off[3] + off[4]]; // x1i,x2i,x3i,x4i,x5i = 1,1,1,0,0 
    val += s[0] * t[0] * u[0] * v[1] * w[0] * rPath[bIdx + off[0]          + off[2] + off[3] + off[4]]; // x1i,x2i,x3i,x4i,x5i = 1,1,1,0,1 
    val += s[0] * t[0] * u[0] * v[0] * w[1] * rPath[bIdx          + off[1] + off[2] + off[3] + off[4]]; // x1i,x2i,x3i,x4i,x5i = 1,1,1,1,0 
    val += s[0] * t[0] * u[0] * v[0] * w[0] * rPath[bIdx + off[0] + off[1] + off[2] + off[3] + off[4]]; // x1i,x2i,x3i,x4i,x5i = 1,1,1,1,1 

    return val;
}
