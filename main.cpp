/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   main.cpp
 * Author: akheckel
 *
 * Created on 23. November 2016, 09:45
 */

#include <cstdlib>
#include <netcdf>
#include <sys/_intsup.h>
#include <sys/unistd.h>
#include "InputParameter.hpp"
#include "S3MetaData.hpp"
#include "S3NcdfData.hpp"
#include "Images.hpp"
#include "Interpolation.hpp"
#include "miscUtils.hpp"
#include "AeroClimatology.hpp"
#include "OceanReflLut.hpp"
#include "AtmosphericLut.hpp"
#include "AodRetrieval.hpp"

using std::cout;
using std::cerr;
using std::endl;
using namespace netCDF;

void addWriteVar(NcFile* ncOut, const std::vector<NcDim>& dimVec, const S3BasicImage<signed char>& img);
void addWriteVar(NcFile* ncOut, const std::vector<NcDim>& dimVec, const S3BasicImage<short>& img);
void addWriteVar(NcFile* ncOut, const std::vector<NcDim>& dimVec, const S3BasicImage<float>& img);
void addWriteVar(NcFile* ncOut, const std::vector<NcDim>& dimVec, const S3BasicImage<double>& img);
void computeOutputImageProp(ImageProperties* outputImgProp, const ImageProperties& inputImgProp, const int& winSize);


/*
 *
 */
int main(int argc, char** argv) {
    double t1, t2;

    InputParameter pars;
    S3NcdfData s3Data(pars);
    ImageProperties outImgProp;
    computeOutputImageProp(&outImgProp, s3Data.s3MetaData.slstrPInfo.nadirImg0500m, pars.winSize);

    cout << s3Data.s3MetaData.productName << endl;
    try {

        /* reading nadir and oblique radiance images */
        t1 = timer();

        s3Data.readNcdf(outImgProp);

        t2 = timer();
        printf("Runtime %f seconds\n", t2 - t1);

        s3Data.convRad2Refl();
        s3Data.verifyInput();
        s3Data.initResultImgs(outImgProp);
        int& imgWidth  = outImgProp.width;  //s3Data.s3MetaData.slstrPInfo.nadirImg0500m.width;
        int& imgHeight = outImgProp.height; //s3Data.s3MetaData.slstrPInfo.nadirImg0500m.height;

/******/

        t1 = timer();

        // read aerosol climatology
        AeroClimatology aerClim(pars.climFileName, s3Data.s3MetaData.getMonth());

        // read ocean LUT
        OceanReflLut ocnLut(pars.ocnLutFileName);

        // read atmospheric LUT
        AtmosphericLut lut(pars.atmLutFileName);

        t2 = timer();
        printf("Time to read Luts %f seconds\n", t2 - t1);

        SlstrPixel slstrPixel;
        int idx;
        for (slstrPixel.x = 0; slstrPixel.x < imgWidth; slstrPixel.x++){
//        for (slstrPixel.x = 129; slstrPixel.x < 130; slstrPixel.x++){
            fprintf(stdout, "processing %6.2f%%\r", (float)(slstrPixel.x)/imgWidth*100.0); fflush(stdout);
            for (slstrPixel.y = 0; slstrPixel.y < imgHeight; slstrPixel.y++){
//            for (slstrPixel.y = 111; slstrPixel.y < 112; slstrPixel.y++){
//                fprintf(stdout, "processing %5d / %5d\n", slstrPixel.x, slstrPixel.y); fflush(stdout);
                idx = slstrPixel.y * imgWidth + slstrPixel.x;
                if ((slstrPixel.y > 0 ) && (s3Data.isValidPixel(idx))){
                    s3Data.getGeoPos(idx, &slstrPixel.geo_pos);
                    s3Data.getViewGeom(idx, &slstrPixel.geom);

                    slstrPixel.pAlt = DEFAULT_PALT;
                    slstrPixel.ocn_wind_speed = DEFAULT_WDSP;
                    slstrPixel.ocn_wind_dir   = DEFAULT_WDIR;
                    slstrPixel.ocn_pigment    = DEFAULT_PIG;

                    slstrPixel.view_clear[0] = ((s3Data.flags.img[idx] & (4+1)) > 0);
                    slstrPixel.view_clear[1] = ((s3Data.flags.img[idx] & (8+2)) > 0);

                    s3Data.getToaReflec(idx, slstrPixel.tarr);
                    if (slstrPixel.x == 0 && slstrPixel.y == 30){
                        printf("Pixel (%d/%d) is valid and at %s\n", slstrPixel.x, slstrPixel.y, slstrPixel.geo_pos.toCstr());
                        printf("Pixel SZA: %f SAA: %f VZA: %f VAA: %f RAZ: %f\n", slstrPixel.geom.nad_sol_zen, slstrPixel.geom.nad_sol_azim, slstrPixel.geom.nad_sat_zen, slstrPixel.geom.nad_sat_azim, slstrPixel.geom.razn);
                        printf("Pixel SZA: %f SAA: %f VZA: %f VAA: %f RAZ: %f\n", slstrPixel.geom.for_sol_zen, slstrPixel.geom.for_sol_azim, slstrPixel.geom.for_sat_zen, slstrPixel.geom.for_sat_azim, slstrPixel.geom.razf);
                        printf("Pixel TOA-N %f %f %f %f %f\n", slstrPixel.tarr[0][0], slstrPixel.tarr[1][0], slstrPixel.tarr[2][0], slstrPixel.tarr[3][0], slstrPixel.tarr[4][0]);
                        printf("Pixel TOA-O %f %f %f %f %f\n", slstrPixel.tarr[0][1], slstrPixel.tarr[1][1], slstrPixel.tarr[2][1], slstrPixel.tarr[3][1], slstrPixel.tarr[4][1]);
                    }
                    //get mixing paramter
                    //prep lut interpolation
                    aerClim.getMixPercentages(slstrPixel.geo_pos, slstrPixel.lutpars.mixing, slstrPixel.lutpars.mix_frac, &slstrPixel.lutpars.climAod);
                    lut.getTetrahedronPoints(&slstrPixel.lutpars, false);
                    try {
                        slstrPixel.lutpars.razni = lut.getInterPar(slstrPixel.geom.razn, lut.razD);
                        slstrPixel.lutpars.razfi = lut.getInterPar(slstrPixel.geom.razf, lut.razD);

                        slstrPixel.lutpars.szani = lut.getInterPar(slstrPixel.geom.nad_sol_zen, lut.szaD);
                        slstrPixel.lutpars.szafi = lut.getInterPar(slstrPixel.geom.for_sol_zen, lut.szaD);

                        slstrPixel.lutpars.vzani = lut.getInterPar(slstrPixel.geom.nad_sat_zen, lut.vzaD);
                        slstrPixel.lutpars.vzafi = lut.getInterPar(slstrPixel.geom.for_sat_zen, lut.vzaD);

                        slstrPixel.lutpars.pAlti = lut.getInterPar(slstrPixel.pAlt, lut.presD);

                        slstrPixel.lutpars.ocn_razni = ocnLut.getInterPar(slstrPixel.geom.razn, ocnLut.razD);
                        slstrPixel.lutpars.ocn_razfi = ocnLut.getInterPar(slstrPixel.geom.razf, ocnLut.razD);
                        slstrPixel.lutpars.ocn_mi  = ocnLut.getInterPar(slstrPixel.lutpars.mix_frac[0], ocnLut.modelD);
                        slstrPixel.lutpars.ocn_wsi = ocnLut.getInterPar(slstrPixel.ocn_wind_speed, ocnLut.wdspD);
                        slstrPixel.lutpars.ocn_wdi = ocnLut.getInterPar(slstrPixel.ocn_wind_dir, ocnLut.wdirD);
                        slstrPixel.lutpars.ocn_pi  = ocnLut.getInterPar(slstrPixel.ocn_pigment, ocnLut.pigcD);

                        ocnLut.get_rho_ocean_wind(&slstrPixel, 0.05, 9);
                        slstrPixel.rho_glint[0] = slstrPixel.rho_surf[3][0];
                        slstrPixel.rho_glint[1] = slstrPixel.rho_surf[3][1];

                        AodRetrieval retrieval(slstrPixel, lut, ocnLut);
                        if ((s3Data.flags.img[idx] & 15) == 3){
                            //clear land in both views
                            retrieval.retrieveAodSizeBrent(false);
                        }
                        else {
                            //clear ocean in either view
                            retrieval.retrieveAodSizeBrent(true);
                        }
                        s3Data.setRetrievalResults(idx, slstrPixel);
                    }
                    catch (std::range_error){
                        fprintf(stderr, "skipping pixel, is not inside LUT\n");
                        //TODO: set possible flag for outside lut and continue
                    }

                }
                else {
                    //printf("Pixel (%d/%d) is invalid\n", slstrPixel.x, slstrPixel.y);
                }

                s3Data.flags.img[idx] &= slstrPixel.qflag;

            }
        }

        // writing radiance data to ncdf file
        t1 = timer();
        printf("Time to retrieve AOD %f seconds\n", t1 - t2);
/*******/
        cout << "writing: " << s3Data.aodOutName << endl;
        NcFile ncOut(s3Data.aodOutName, NcFile::replace);
        NcDim xDim = ncOut.addDim("columns", imgWidth);
        NcDim yDim = ncOut.addDim("rows", imgHeight);
        std::vector<NcDim> dimVec;
        dimVec.push_back(yDim);
        dimVec.push_back(xDim);

        for (int iView = 0; iView < N_SLSTR_VIEWS; iView++){
            for (int iBand = 0; iBand < N_SLSTR_BANDS; iBand++){
                addWriteVar(&ncOut, dimVec, s3Data.s3RadImgs[iView][iBand]);
                addWriteVar(&ncOut, dimVec, s3Data.s3SdrImgs[iView][iBand]);
                addWriteVar(&ncOut, dimVec, s3Data.s3RPathImgs[iView][iBand]);
                addWriteVar(&ncOut, dimVec, s3Data.s3TDownImgs[iView][iBand]);
                addWriteVar(&ncOut, dimVec, s3Data.s3TUpImgs[iView][iBand]);
                addWriteVar(&ncOut, dimVec, s3Data.s3TGasImgs[iView][iBand]);
                addWriteVar(&ncOut, dimVec, s3Data.s3SpherAImgs[iView][iBand]);
                addWriteVar(&ncOut, dimVec, s3Data.s3DifFracImgs[iView][iBand]);
            }

            addWriteVar(&ncOut, dimVec, s3Data.s3LatImgs[iView]);
            addWriteVar(&ncOut, dimVec, s3Data.s3LonImgs[iView]);
            addWriteVar(&ncOut, dimVec, s3Data.s3SzaImgs[iView]);
            addWriteVar(&ncOut, dimVec, s3Data.s3SaaImgs[iView]);
            addWriteVar(&ncOut, dimVec, s3Data.s3VzaImgs[iView]);
            addWriteVar(&ncOut, dimVec, s3Data.s3VaaImgs[iView]);
        }
        for (int iBand = 0; iBand < N_SLSTR_BANDS; iBand++){
            addWriteVar(&ncOut, dimVec, s3Data.s3AodImgs[iBand]);
        }
        for (int i = 0; i < 3; i++){
            addWriteVar(&ncOut, dimVec, s3Data.s3AerFracImgs[i]);
        }

        addWriteVar(&ncOut, dimVec, s3Data.flags);

        ncOut.close();
        t2 = timer();
        printf("Runtime %f seconds\n", t2 - t1);


        /* writing hi res flag data to ncdf file */
        /*imgWidth  = s3Data.s3MetaData.slstrPInfo.nadirImg0500m.width;
        imgHeight = s3Data.s3MetaData.slstrPInfo.nadirImg0500m.height;
        replaceStringInPlace(s3Data.aodOutName, ".nc", "_flags.nc");
        cout << "writing: " << s3Data.aodOutName << endl;
        NcFile ncOut1(s3Data.aodOutName, NcFile::replace);
        xDim = ncOut1.addDim("columns", imgWidth);
        yDim = ncOut1.addDim("rows", imgHeight);
        dimVec.clear();
        dimVec.push_back(yDim);
        dimVec.push_back(xDim);

        addWriteVar(&ncOut1, dimVec, s3Data.s3LandImg);
        addWriteVar(&ncOut1, dimVec, s3Data.s3CloudImg);
        addWriteVar(&ncOut1, dimVec, s3Data.s3ValidImg);

        ncOut1.close();*/

    }
    catch (exceptions::NcException& e){
        cerr << "unrecoverable error, exiting...\n";
        cerr << e.what();
        return e.errorCode();
    }
    cout << "finished!" << endl;
    return 0;
}

void addWriteVar(NcFile* ncOut, const std::vector<NcDim>& dimVec, const S3BasicImage<signed char>& img){
    NcVar var = ncOut->addVar(img.name, ncByte, dimVec);

    //var.setCompression(true, true, 5);

    if (img.hasOffset){
        var.putAtt("add_offset", ncDouble, img.valOffset);
    }
    if (img.hasScale){
        var.putAtt("scale_factor", ncDouble, img.valScale);
    }
    if (img.hasNoData){
        var.putAtt("_FillValue", ncByte, img.noData);
    }
    if (img.isSpecBand){
        var.putAtt("wavelength", ncFloat, img.wvl);
    }
    //var.putAtt("coordinates", std::string("latitude_an longitude_an"));

    var.putVar(img.img);
}

void addWriteVar(NcFile* ncOut, const std::vector<NcDim>& dimVec, const S3BasicImage<short>& img){
    NcVar var = ncOut->addVar(img.name, ncShort, dimVec);

    var.setCompression(true, true, 5);

    if (img.hasOffset){
        var.putAtt("add_offset", ncDouble, img.valOffset);
    }
    if (img.hasScale){
        var.putAtt("scale_factor", ncDouble, img.valScale);
    }
    if (img.hasNoData){
        var.putAtt("_FillValue", ncShort, img.noData);
    }
    if (img.isSpecBand){
        var.putAtt("wavelength", ncFloat, img.wvl);
    }
    var.putAtt("coordinates", std::string("latitude_an longitude_an"));

    var.putVar(img.img);
}

void addWriteVar(NcFile* ncOut, const std::vector<NcDim>& dimVec, const S3BasicImage<float>& img){
    NcVar var = ncOut->addVar(img.name, ncFloat, dimVec);

    var.setCompression(true, true, 5);

    if (img.hasOffset){
        var.putAtt("add_offset", ncDouble, img.valOffset);
    }
    if (img.hasScale){
        var.putAtt("scale_factor", ncDouble, img.valScale);
    }
    if (img.hasNoData){
        var.putAtt("_FillValue", ncFloat, img.noData);
    }
    if (img.isSpecBand){
        var.putAtt("wavelength", ncFloat, img.wvl);
    }
    var.putVar(img.img);
}

void addWriteVar(NcFile* ncOut, const std::vector<NcDim>& dimVec, const S3BasicImage<double>& img){
    NcVar var = ncOut->addVar(img.name, ncDouble, dimVec);

    var.setCompression(true, true, 5);

    if (img.hasOffset){
        var.putAtt("add_offset", ncDouble, img.valOffset);
    }
    if (img.hasScale){
        var.putAtt("scale_factor", ncDouble, img.valScale);
    }
    if (img.hasNoData){
        var.putAtt("_FillValue", ncDouble, img.noData);
    }
    if (img.isSpecBand){
        var.putAtt("wavelength", ncFloat, img.wvl);
    }
    var.putVar(img.img);
}

void computeOutputImageProp(ImageProperties* outputImgProp, const ImageProperties& inputImgProp, const int& winSize){

    if (winSize < 1) {
        std::string msg("computeOutputImageProp(): binning factor (winSize) < 1!\n");
        throw std::out_of_range(msg);
    }
    else if (winSize == 1){
        *outputImgProp = inputImgProp;
    }
    else {
        outputImgProp->isBinned = true;
        outputImgProp->binSize = winSize;
        outputImgProp->width = inputImgProp.width / winSize;
        if (inputImgProp.width > outputImgProp->width * winSize) outputImgProp->width++;

        outputImgProp->height = inputImgProp.height / winSize;
        if (inputImgProp.height > outputImgProp->height * winSize) outputImgProp->height++;

        outputImgProp->nPix = outputImgProp->width * outputImgProp->height;

        outputImgProp->xRes = inputImgProp.xRes * winSize;
        outputImgProp->yRes = inputImgProp.yRes * winSize;

        outputImgProp->xOff = (int)((double)(inputImgProp.xOff) / winSize + 0.5);  //* inputImgProp.xRes / outputImgProp->xRes;
        outputImgProp->yOff = (int)((double)(inputImgProp.yOff) / winSize + 0.5);  //* inputImgProp.yRes / outputImgProp->yRes;
    }
}

