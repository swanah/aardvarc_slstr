/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   S3NcdfData.cpp
 * Author: akheckel
 * 
 * Created on 25. November 2016, 10:02
 */

#include <limits>
#include <cmath>
#include <netcdf>
#include <vector>
#include "S3NcdfData.hpp"
#include "Images.hpp"
#include "Interpolation.hpp"
#include "miscUtils.hpp"

using std::cout;
using std::cerr;
using std::endl;
using namespace netCDF;

/* const definitions */
    const std::string S3NcdfData::SDR_NAMES[N_SLSTR_VIEWS][5] = {
        {"S1_SDR_n", "S2_SDR_n", "S3_SDR_n", "S5_SDR_n", "S6_SDR_n"},
        {"S1_SDR_o", "S2_SDR_o", "S3_SDR_o", "S5_SDR_o", "S6_SDR_o"}
    };
    const std::string S3NcdfData::RPATH_NAMES[N_SLSTR_VIEWS][5] = {
        {"S1_rPath_n", "S2_rPath_n", "S3_rPath_n", "S5_rPath_n", "S6_rPath_n"},
        {"S1_rPath_o", "S2_rPath_o", "S3_rPath_o", "S5_rPath_o", "S6_rPath_o"}
    };
    const std::string S3NcdfData::TDOWN_NAMES[N_SLSTR_VIEWS][5] = {
        {"S1_tDown_n", "S2_tDown_n", "S3_tDown_n", "S5_tDown_n", "S6_tDown_n"},
        {"S1_tDown_o", "S2_tDown_o", "S3_tDown_o", "S5_tDown_o", "S6_tDown_o"}
    };
    const std::string S3NcdfData::TUP_NAMES[N_SLSTR_VIEWS][5] = {
        {"S1_tUp_n", "S2_tUp_n", "S3_tUp_n", "S5_tUp_n", "S6_tUp_n"},
        {"S1_tUp_o", "S2_tUp_o", "S3_tUp_o", "S5_tUp_o", "S6_tUp_o"}
    };
    const std::string S3NcdfData::TGAS_NAMES[N_SLSTR_VIEWS][5] = {
        {"S1_tGas_n", "S2_tGas_n", "S3_tGas_n", "S5_tGas_n", "S6_tGas_n"},
        {"S1_tGas_o", "S2_tGas_o", "S3_tGas_o", "S5_tGas_o", "S6_tGas_o"}
    };
    const std::string S3NcdfData::SPHERA_NAMES[N_SLSTR_VIEWS][5] = {
        {"S1_SpherA_n", "S2_SpherA_n", "S3_SpherA_n", "S5_SpherA_n", "S6_SpherA_n"},
        {"S1_SpherA_o", "S2_SpherA_o", "S3_SpherA_o", "S5_SpherA_o", "S6_SpherA_o"}
    };
    const std::string S3NcdfData::DIFFRAC_NAMES[N_SLSTR_VIEWS][5] = {
        {"S1_DifFrac_n", "S2_DifFrac_n", "S3_DifFrac_n", "S5_DifFrac_n", "S6_DifFrac_n"},
        {"S1_DifFrac_o", "S2_DifFrac_o", "S3_DifFrac_o", "S5_DifFrac_o", "S6_DifFrac_o"}
    };

    const std::string S3NcdfData::AOD_NAMES[N_SLSTR_BANDS] = {"AOD_0550", "AOD_0659", "AOD_0865", "AOD_1610", "AOD_2250"};
    const std::string S3NcdfData::AER_FRAC_NAMES[N_AER_FRAC] = {"fot_clim", "fineOfTotal", "weakAbsOfFine", "dustOfCoarse"};
    const std::string S3NcdfData::FMIN_NAME     = "fmin";
    const std::string S3NcdfData::SSA_NAME      = "SSA550";
    const std::string S3NcdfData::ABS_AOD_NAME  = "AAOD550";
    const std::string S3NcdfData::D_AOD_NAME    = "D_AOD550";
    const std::string S3NcdfData::FM_AOD_NAME   = "FM_AOD550";
    const std::string S3NcdfData::ANGSTROM_NAME = "ANG550_870";
    const std::string S3NcdfData::UNC_NAME[N_SLSTR_BANDS] = {"AOD_0550_uncertainty", "AOD_0659_uncertainty", "AOD_0865_uncertainty", "AOD_1610_uncertainty", "AOD_2250_uncertainty"};
    const std::string S3NcdfData::RAZ_NAME[N_SLSTR_VIEWS] = {"rel_azimuth_an", "rel_azimuth_ao"};

    const std::string S3NcdfData::LAT_CNR_NAMES[4] = {"pixel_corner_latitude1", "pixel_corner_latitude2", "pixel_corner_latitude3", "pixel_corner_latitude4"};
    const std::string S3NcdfData::LON_CNR_NAMES[4] = {"pixel_corner_longitude1", "pixel_corner_longitude2", "pixel_corner_longitude3", "pixel_corner_longitude4"};

    const std::string S3NcdfData::TIME_NAME = "time";

    const float S3NcdfData::angWvlLog = -1.0 / log( 550. / 865.);
//
//public
//

S3NcdfData::S3NcdfData(const InputParameter& inPar) : pars(inPar), s3DataDir(inPar.slstrProductDir) {
//    : s3DataDir(inPar.slstrProductDir), s3MetaData(inPar.s3MD), szaLim(inPar.szaLimit),
//      doGeoSubset(inPar.doGeoSubset), latLim(inPar.latLim), lonLim(inPar.lonLim) {
    
    isLandMaskAvailable = false;
    isFlagsImgAvailable = false;
}

S3NcdfData::~S3NcdfData() {
}

/**
 * read the input product from all required netcdf files
 * all images are binned, interpolated and shifted according 
 * to the binning parameters and resolution
 * @param outImgProp
 */
void S3NcdfData::readNcdf(const ImageProperties& outImgProp) {
    float wvl[] = {555., 659., 865., 1610., 2250.};
    offCorr[0] = offCorr[1] = 0;
    sceneOutputWidth = outImgProp.width;
    sceneOutputHeight = outImgProp.height;
    readIrrad(s3Irrad);
    createLandMask();
    createNanMask();
    if (pars.useSCloudS3SU) {
        createMyCloudMask();
    }
    else {
        createCloudMask();
    }
    createValidMask();
    flags = S3BasicImage<short>(outImgProp);
    flags.name = "aod_flags";
    flags.initImgArray();
    isFlagsImgAvailable = true;

    s3NPixImg = S3BasicImage<short>(outImgProp);
    s3NPixImg.name = "npix";
    s3NPixImg.setValidLimits(0, pars.winSize*pars.winSize);
    s3NPixImg.initImgArray();
    

    // reading nadir and oblique radiance images
    NcdfImageType imgType;
    std::string ncdfName;
    for (int iView = 0; iView < N_SLSTR_VIEWS; iView++){
        imgType = (iView == NADIR_VIEW) ? Nadir0500 : Obliq0500;

        //read img and bin to output image properties
        for (int iBand = 0; iBand < N_SLSTR_BANDS; iBand++){
            s3RadImgsF[iView][iBand] = S3BasicImage<double>(outImgProp);
            s3RadImgsF[iView][iBand].setValidLimits(0, 10000);
            s3RadImgsF[iView][iBand].setFillVal(-1.);
            s3RadImgsF[iView][iBand].setWvl(wvl[iBand]);
            
            ncdfName = pars.slstrProductDir + "/" + S3MetaData::CHANNEL_RAD_NAME[iView][iBand] + ".nc";
            /**
            if (iView == 1 && iBand < 3){
                offCorr[0] = -1; offCorr[1] = 2;
            }
            else if (iView == 0 && iBand > 2){
                offCorr[0] = 1; offCorr[1] = 0;
            }
            else if (iView == 1 && iBand > 2){
                offCorr[0] = -3; offCorr[1] = 2;
            }
            /**/
            
            readImgBinned(&s3RadImgsF[iView][iBand], ncdfName, S3MetaData::CHANNEL_RAD_NAME[iView][iBand], imgType);

            // read radiances to <short> images
            s3RadImgs[iView][iBand] = S3BasicImage<short>(outImgProp);
            s3RadImgs[iView][iBand].setValidLimits(0, 20000);
            s3RadImgs[iView][iBand].setWvl(wvl[iBand]);
            readImgBinned(&s3RadImgs[iView][iBand], ncdfName, S3MetaData::CHANNEL_RAD_NAME[iView][iBand], imgType);
        }

        ncdfName = pars.slstrProductDir + "/" + S3MetaData::GEODETIC_NAME[iView] + ".nc";
        s3LatImgs[iView] = S3BasicImage<double>(outImgProp);
        s3LatImgs[iView].setValidLimits(-90., 90.);
        s3LatImgs[iView].setFillVal(-999.);
        if (iView == NADIR_VIEW){
            for (int iCnr = 0; iCnr < 4; iCnr++) {
                s3LatCnrImgs[iCnr] = S3BasicImage<double>(outImgProp);
                s3LatCnrImgs[iCnr].setValidLimits(-90., 90.);
                s3LatCnrImgs[iCnr].setFillVal(-999.);
                s3LatCnrImgs[iCnr].name = LAT_CNR_NAMES[iCnr];
            }
        }
        readImgBinned(&s3LatImgs[iView], ncdfName, S3MetaData::LAT_NAME[iView], imgType);

        s3LonImgs[iView] = S3BasicImage<double>(outImgProp);
        s3LonImgs[iView].setValidLimits(-180., 180.);
        s3LonImgs[iView].setFillVal(-999.);
        if (iView == NADIR_VIEW){
            for (int iCnr = 0; iCnr < 4; iCnr++) {
                s3LonCnrImgs[iCnr] = S3BasicImage<double>(outImgProp);
                s3LonCnrImgs[iCnr].setValidLimits(-180., 180.);
                s3LonCnrImgs[iCnr].setFillVal(-999.);
                s3LonCnrImgs[iCnr].name = LON_CNR_NAMES[iCnr];
            }
        }
        readImgBinned(&s3LonImgs[iView], ncdfName, S3MetaData::LON_NAME[iView], imgType);

        imgType = (iView == NADIR_VIEW) ? NadirTpg : ObliqTpg;
        
        ncdfName = pars.slstrProductDir + "/" + S3MetaData::GEOMETRY_NAME[iView] + ".nc";
        s3SzaImgs[iView] = S3BasicImage<double>(outImgProp);
        s3SzaImgs[iView].setValidLimits(0., 90.);
        readImgBinned(&s3SzaImgs[iView], ncdfName, S3MetaData::SZA_NAME[iView], imgType);
        
        s3SaaImgs[iView] = S3BasicImage<double>(outImgProp);
        s3SaaImgs[iView].setValidLimits(0., 360.);
        readImgBinned(&s3SaaImgs[iView], ncdfName, S3MetaData::SAA_NAME[iView], imgType);
        
        s3VzaImgs[iView] = S3BasicImage<double>(outImgProp);
        s3SzaImgs[iView].setValidLimits(0., 90.);
        readImgBinned(&s3VzaImgs[iView], ncdfName, S3MetaData::VZA_NAME[iView], imgType);
        
        s3VaaImgs[iView] = S3BasicImage<double>(outImgProp);
        s3VaaImgs[iView].setValidLimits(0., 360.);
        readImgBinned(&s3VaaImgs[iView], ncdfName, S3MetaData::VAA_NAME[iView], imgType);
    }

    ncdfName = pars.slstrProductDir + "/" + S3MetaData::MET_NAME + ".nc";
    s3PresImg = S3BasicImage<float>(outImgProp);
    s3PresImg.setValidLimits(0., 1100.);
    readImgBinned(&s3PresImg, ncdfName, S3MetaData::PRES_NAME, imgType);
    
    ncdfName = pars.slstrProductDir + "/" + S3MetaData::TIME_NAME[0] + ".nc";
    s3TimeImg = S3BasicImage<unsigned int>(outImgProp);
    s3TimeImg.name = TIME_NAME;
    s3TimeImg.setFillVal((unsigned int)(0));
    s3TimeImg.initImgArray((unsigned int)(0));
    computeAvgTimeImg(&s3TimeImg);

}

/**
 * convert all radiance images to reflectance images
 * using the provided irradiance information in the L1b product
 */
void S3NcdfData::convRad2Refl(){
    double rad;
    const double deg2rad = M_PI / 180;
    for (int iView = 0; iView < N_SLSTR_VIEWS; iView++){
        for (int iBand = 0; iBand < N_SLSTR_BANDS; iBand++){
            double origScale = s3RadImgs[iView][iBand].valScale;
            double origOffset = s3RadImgs[iView][iBand].valOffset;
            short  origNoData = s3RadImgs[iView][iBand].noData;
            s3RadImgs[iView][iBand].valScale = 1e-4;
            s3RadImgs[iView][iBand].valOffset = 0;
            s3RadImgs[iView][iBand].noData = -10000;
            //replaceStringInPlace(s3RadImgs[iView][iBand].name, "radiance", "reflec");
            replaceStringInPlace(s3RadImgsF[iView][iBand].name, "radiance", "reflec");
            for (int i = 0; i < s3RadImgs[iView][iBand].imgP.nPix; i++){
                if (s3RadImgsF[iView][iBand].img[i] != s3RadImgsF[iView][iBand].noData){
                    rad = (double)(s3RadImgs[iView][iBand].img[i]) * origScale + origOffset;
                    rad *= M_PI / (s3Irrad[iView][iBand] * cos(s3SzaImgs[iView].img[i] * deg2rad));
                    //s3RadImgs[iView][iBand].img[i] = (rad - s3RadImgs[iView][iBand].valOffset) / s3RadImgs[iView][iBand].valScale;
                    s3RadImgsF[iView][iBand].img[i] *= M_PI / (s3Irrad[iView][iBand] * cos(s3SzaImgs[iView].img[i] * deg2rad));
                }
                else {
                    //s3RadImgs[iView][iBand].img[i] = s3RadImgs[iView][iBand].noData;
                    s3RadImgsF[iView][iBand].img[i] = s3RadImgsF[iView][iBand].noData;
                }
            }
        }
    }
}

/**
 * apply valid flags (cloud) to falgs and reflectance images
 * testing for consistency
 * i.e. land retrieval requires both nadir and oblique to be valid
 *      ocean retrieval only one view needs to be valid
 */
void S3NcdfData::verifyInput(){
    char valid;
    int iBand;
    
    for (int i = 0; i < s3RadImgsF[0][0].imgP.nPix; i++){
        valid = (flags.img[i] & 15);
        if (valid == 1 || valid == 2){
            // only single view over land is not sufficient 
            // -> setpixel invalid
            flags.img[i] &= ~(3);
            for (iBand = 0; iBand < N_SLSTR_BANDS; iBand++) {
                s3RadImgsF[0][iBand].img[i] = s3RadImgsF[0][iBand].noData;
                s3RadImgsF[1][iBand].img[i] = s3RadImgsF[1][iBand].noData;
            }
        }
        else if (valid == 3){
            // pixel is valid land -> test all reflec are available
            for (iBand = 0; iBand < N_SLSTR_BANDS; iBand++) {
                valid = valid && (s3RadImgsF[0][iBand].img[i] != s3RadImgsF[0][iBand].noData);
                valid = valid && (s3RadImgsF[1][iBand].img[i] != s3RadImgsF[1][iBand].noData);
            }
            if (!valid){
                // at least one reflec missing -> set pixel invalid for land
                flags.img[i] &= ~(3);
                for (iBand = 0; iBand < N_SLSTR_BANDS; iBand++) {
                    s3RadImgsF[0][iBand].img[i] = s3RadImgsF[0][iBand].noData;
                    s3RadImgsF[1][iBand].img[i] = s3RadImgsF[1][iBand].noData;
                }
            }
        }
        else if (valid == 4){
            // test ocean nadir
            for (iBand = 0; iBand < N_SLSTR_BANDS; iBand++) {
                valid = valid && (s3RadImgsF[0][iBand].img[i] != s3RadImgsF[0][iBand].noData);
            }
            if (!valid){
                flags.img[i] &= ~(4);
                for (iBand = 0; iBand < N_SLSTR_BANDS; iBand++) {
                    s3RadImgsF[0][iBand].img[i] = s3RadImgsF[0][iBand].noData;
                }
            }
        }
        else if (valid == 8){
            // test ocean oblique
            for (iBand = 0; iBand < N_SLSTR_BANDS; iBand++) {
                valid = valid && (s3RadImgsF[1][iBand].img[i] != s3RadImgsF[1][iBand].noData);
            }
            if (!valid){
                flags.img[i] &= ~(8);
                for (iBand = 0; iBand < N_SLSTR_BANDS; iBand++) {
                    s3RadImgsF[0][iBand].img[i] = s3RadImgsF[0][iBand].noData;
                }
            }
        }
    }
    
    /**
    for (int i = 0; i < s3RadImgs[0][0].imgP.nPix; i++){
        valid = (flags.img[i] & 15);
        if (valid == 1 || valid == 2){
            // only single view over land is not sufficient 
            // -> setpixel invalid
            flags.img[i] &= ~(3);
            for (iBand = 0; iBand < N_SLSTR_BANDS; iBand++) {
                s3RadImgs[0][iBand].img[i] = s3RadImgs[0][iBand].noData;
                s3RadImgs[1][iBand].img[i] = s3RadImgs[1][iBand].noData;
            }
        }
        else if (valid == 3){
            // pixel is valid land -> test all reflec are available
            for (iBand = 0; iBand < N_SLSTR_BANDS; iBand++) {
                valid = valid && (s3RadImgs[0][iBand].img[i] != s3RadImgs[0][iBand].noData);
                valid = valid && (s3RadImgs[1][iBand].img[i] != s3RadImgs[1][iBand].noData);
            }
            if (!valid){
                // at least one reflec missing -> set pixel invalid for land
                flags.img[i] &= ~(3);
                for (iBand = 0; iBand < N_SLSTR_BANDS; iBand++) {
                    s3RadImgs[0][iBand].img[i] = s3RadImgs[0][iBand].noData;
                    s3RadImgs[1][iBand].img[i] = s3RadImgs[1][iBand].noData;
                }
            }
        }
        else if (valid == 4){
            // test ocean nadir
            for (iBand = 0; iBand < N_SLSTR_BANDS; iBand++) {
                valid = valid && (s3RadImgs[0][iBand].img[i] != s3RadImgs[0][iBand].noData);
            }
            if (!valid){
                flags.img[i] &= ~(4);
                for (iBand = 0; iBand < N_SLSTR_BANDS; iBand++) {
                    s3RadImgs[0][iBand].img[i] = s3RadImgs[0][iBand].noData;
                }
            }
        }
        else if (valid == 8){
            // test ocean oblique
            for (iBand = 0; iBand < N_SLSTR_BANDS; iBand++) {
                valid = valid && (s3RadImgs[1][iBand].img[i] != s3RadImgs[1][iBand].noData);
            }
            if (!valid){
                flags.img[i] &= ~(8);
                for (iBand = 0; iBand < N_SLSTR_BANDS; iBand++) {
                    s3RadImgs[0][iBand].img[i] = s3RadImgs[0][iBand].noData;
                }
            }
        }
    }
    /**/
}

/**
 * initialise all pixels of the results images with fill values
 * @param outImgProp
 */
void S3NcdfData::initResultImgs(const ImageProperties& outImgProp){
    s3RazImgs[0] = S3BasicImage<float>(outImgProp);
    s3RazImgs[0].name = RAZ_NAME[0];
    s3RazImgs[0].setFillVal((float)(-1));
    s3RazImgs[0].initImgArray((float)(-1));
    s3RazImgs[1] = S3BasicImage<float>(outImgProp);
    s3RazImgs[1].name = RAZ_NAME[1];
    s3RazImgs[1].setFillVal((float)(-1));
    s3RazImgs[1].initImgArray((float)(-1));
    for (int iBand = 0; iBand < N_SLSTR_BANDS; iBand++) {
        s3AodImgs[iBand] = S3BasicImage<float>(outImgProp);
        s3AodImgs[iBand].name = AOD_NAMES[iBand];
        s3AodImgs[iBand].setFillVal((float)(-1));
        s3AodImgs[iBand].setValidLimits((float)(0.003), (float)(4));
        s3AodImgs[iBand].initImgArray((float)(-1));

        s3UncImgs[iBand] = S3BasicImage<float>(outImgProp);
        s3UncImgs[iBand].name = UNC_NAME[iBand];
        s3UncImgs[iBand].setFillVal((float)(-1));
        s3UncImgs[iBand].initImgArray((float)(-1));

        for (int iView = 0; iView < N_SLSTR_VIEWS; iView++) {
            s3SdrImgs[iView][iBand] = S3BasicImage<short>(outImgProp);
            s3SdrImgs[iView][iBand].name = SDR_NAMES[iView][iBand];
            s3SdrImgs[iView][iBand].setFillVal((short)(-30000));
            s3SdrImgs[iView][iBand].initImgArray((short)(-30000));
            s3SdrImgs[iView][iBand].setScaleOffset(1e-4, 0);
            s3SdrImgs[iView][iBand].setValidLimits(-10000, 20000);

            /**
            s3RPathImgs[iView][iBand] = S3BasicImage<float>(outImgProp);
            s3RPathImgs[iView][iBand].name = RPATH_NAMES[iView][iBand];
            s3RPathImgs[iView][iBand].setFillVal((float)(-1));
            s3RPathImgs[iView][iBand].initImgArray((float)(-1));
            
            s3TDownImgs[iView][iBand] = S3BasicImage<float>(outImgProp);
            s3TDownImgs[iView][iBand].name = TDOWN_NAMES[iView][iBand];
            s3TDownImgs[iView][iBand].setFillVal((float)(-1));
            s3TDownImgs[iView][iBand].initImgArray((float)(-1));
            
            s3TUpImgs[iView][iBand] = S3BasicImage<float>(outImgProp);
            s3TUpImgs[iView][iBand].name = TUP_NAMES[iView][iBand];
            s3TUpImgs[iView][iBand].setFillVal((float)(-1));
            s3TUpImgs[iView][iBand].initImgArray((float)(-1));
            
            s3TGasImgs[iView][iBand] = S3BasicImage<float>(outImgProp);
            s3TGasImgs[iView][iBand].name = TGAS_NAMES[iView][iBand];
            s3TGasImgs[iView][iBand].setFillVal((float)(-1));
            s3TGasImgs[iView][iBand].initImgArray((float)(-1));
            
            s3SpherAImgs[iView][iBand] = S3BasicImage<float>(outImgProp);
            s3SpherAImgs[iView][iBand].name = SPHERA_NAMES[iView][iBand];
            s3SpherAImgs[iView][iBand].setFillVal((float)(-1));
            s3SpherAImgs[iView][iBand].initImgArray((float)(-1));
            
            s3DifFracImgs[iView][iBand] = S3BasicImage<float>(outImgProp);
            s3DifFracImgs[iView][iBand].name = DIFFRAC_NAMES[iView][iBand];
            s3DifFracImgs[iView][iBand].setFillVal((float)(-1));
            s3DifFracImgs[iView][iBand].initImgArray((float)(-1));
            /**/
        }
    }
    for (int i = 0; i < N_AER_FRAC; i++){
        s3AerFracImgs[i] = S3BasicImage<float>(outImgProp);
        s3AerFracImgs[i].name = AER_FRAC_NAMES[i];
        s3AerFracImgs[i].setFillVal((float)(-1));
        s3AerFracImgs[i].initImgArray((float)(-1));
    }
    s3FminImg = S3BasicImage<float>(outImgProp);
    s3FminImg.name = FMIN_NAME;
    s3FminImg.setFillVal((float)(-1));
    s3FminImg.initImgArray((float)(-1));
    
    s3SsaImg = S3BasicImage<float>(outImgProp);
    s3SsaImg.name = SSA_NAME;
    s3SsaImg.setFillVal((float)(-1));
    s3SsaImg.initImgArray((float)(-1));
    
    s3AbsAodImg = S3BasicImage<float>(outImgProp);
    s3AbsAodImg.name = ABS_AOD_NAME;
    s3AbsAodImg.setFillVal((float)(-1));
    s3AbsAodImg.initImgArray((float)(-1));
    
    s3DustAodImg = S3BasicImage<float>(outImgProp);
    s3DustAodImg.name = D_AOD_NAME;
    s3DustAodImg.setFillVal((float)(-1));
    s3DustAodImg.initImgArray((float)(-1));
    
    s3FmAodImg = S3BasicImage<float>(outImgProp);
    s3FmAodImg.name = FM_AOD_NAME;
    s3FmAodImg.setFillVal((float)(-1));
    s3FmAodImg.initImgArray((float)(-1));
    
    s3AngstromImg = S3BasicImage<float>(outImgProp);
    s3AngstromImg.name = ANGSTROM_NAME;
    s3AngstromImg.setFillVal((float)(-1));
    s3AngstromImg.initImgArray((float)(-1));
    
}

/**
 * test if pixel geolocation and geometry parameters are valid for the retrieval
 * @param idx
 * @return 
 */
bool S3NcdfData::isValidPixel(const int& idx){
    if (s3SzaImgs[0].img[idx] > pars.szaLimit){
        flags.img[idx] &= ~(1 + 4);
    }
    if (s3SzaImgs[1].img[idx] > pars.szaLimit){
        flags.img[idx] &= ~(2 + 8);
    }
    if (pars.doGeoSubset) {
        if ( s3LatImgs[0].img[idx] < pars.latLim[0] || s3LatImgs[0].img[idx] > pars.latLim[1]
             || s3LonImgs[0].img[idx] < pars.lonLim[0] && s3LonImgs[0].img[idx] > pars.lonLim[1] ){
            flags.img[idx] &= ~(1 + 4);
        }
        if ( s3LatImgs[1].img[idx] < pars.latLim[0] || s3LatImgs[1].img[idx] > pars.latLim[1]
             || s3LonImgs[1].img[idx] < pars.lonLim[0] && s3LonImgs[1].img[idx] > pars.lonLim[1] ){
            flags.img[idx] &= ~(2 + 8);
        }
    }

    return ((flags.img[idx] & 15) >= 3);
}

void S3NcdfData::getGeoPos(const int& idx, GeoPos* gp){
    gp->lat = s3LatImgs[0].img[idx];
    gp->lon = s3LonImgs[0].img[idx];
}

void S3NcdfData::getViewGeom(const int& idx, ViewGeom* vg){
    vg->nad_sol_zen  = std::isnan(s3SzaImgs[0].img[idx]) ? 0 : s3SzaImgs[0].img[idx];
    vg->nad_sol_azim = std::isnan(s3SaaImgs[0].img[idx]) ? 0 : s3SaaImgs[0].img[idx];
    vg->nad_sat_zen  = std::isnan(s3VzaImgs[0].img[idx]) ? 0 : s3VzaImgs[0].img[idx];
    vg->nad_sat_azim = std::isnan(s3VaaImgs[0].img[idx]) ? 0 : s3VaaImgs[0].img[idx];
    vg->razn = fabsf(vg->nad_sat_azim - vg->nad_sol_azim);
    if (vg->razn > 180.0) vg->razn = 360.0 - vg->razn;
    s3RazImgs[0].img[idx] = vg->razn;

    vg->for_sol_zen  = std::isnan(s3SzaImgs[1].img[idx]) ? 0 : s3SzaImgs[1].img[idx];
    vg->for_sol_azim = std::isnan(s3SaaImgs[1].img[idx]) ? 0 : s3SaaImgs[1].img[idx];
    vg->for_sat_zen  = std::isnan(s3VzaImgs[1].img[idx]) ? 0 : s3VzaImgs[1].img[idx];
    vg->for_sat_azim = std::isnan(s3VaaImgs[1].img[idx]) ? 0 : s3VaaImgs[1].img[idx];
    vg->razf = fabsf(vg->for_sat_azim - vg->for_sol_azim);
    if (vg->razf > 180.0) vg->razf = 360.0 - vg->razf;
    s3RazImgs[1].img[idx] = vg->razf;
}

void S3NcdfData::getPres(const int& idx, float* pres){
    *pres  = std::isnan(s3PresImg.img[idx]) ? DEFAULT_PALT : s3PresImg.img[idx];
}

void S3NcdfData::getToaReflec(const int& idx, float tarr[][N_SLSTR_VIEWS]){
    for (int iBand = 0; iBand < N_SLSTR_BANDS; iBand++) {
        for (int iView = 0; iView < N_SLSTR_VIEWS; iView++) {
            tarr[iBand][iView] = s3RadImgs[iView][iBand].img[idx] * s3RadImgs[iView][iBand].valScale + s3RadImgs[iView][iBand].valOffset;
        }
    }
}

void S3NcdfData::getToaReflecF(const int& idx, float tarr[][N_SLSTR_VIEWS]){
    for (int iBand = 0; iBand < N_SLSTR_BANDS; iBand++) {
        for (int iView = 0; iView < N_SLSTR_VIEWS; iView++) {
            tarr[iBand][iView] = s3RadImgsF[iView][iBand].img[idx];
        }
    }
}

void S3NcdfData::setRetrievalResults(const int& idx, SlstrPixel& pix){
    setBit(&pix.qflag, AOD_ZERO, (pix.aod < 0.003));
    setBit(&pix.qflag, AOD_HIGH, (pix.aod > 3.003));
    if (pix.aod > 0.003 /*&& pix.aod < 3.003*/){
        for (int iBand = 0; iBand < N_SLSTR_BANDS; iBand++) {
            s3AodImgs[iBand].img[idx] = pix.aod * pix.spec_aod_fac[iBand];
            s3UncImgs[iBand].img[idx] = pix.ediff * pix.spec_aod_fac[iBand];
            for (int iView = 0; iView < N_SLSTR_VIEWS; iView++) {
                if (pix.view_clear[iView]){
                    s3SdrImgs[iView][iBand].img[idx] = (short)((pix.RR[iBand][iView] - s3SdrImgs[iView][iBand].valOffset) / s3SdrImgs[iView][iBand].valScale);
                    /**
                    s3RPathImgs[iView][iBand].img[idx] = pix.inv_coef[iBand][iView].rPath;
                    s3TDownImgs[iView][iBand].img[idx] = pix.inv_coef[iBand][iView].tDown;
                    s3TUpImgs[iView][iBand].img[idx]   = pix.inv_coef[iBand][iView].tUp;
                    s3TGasImgs[iView][iBand].img[idx]  = pix.inv_coef[iBand][iView].tGas;
                    s3SpherAImgs[iView][iBand].img[idx]  = pix.inv_coef[iBand][iView].spherAlb;
                    s3DifFracImgs[iView][iBand].img[idx] = pix.inv_coef[iBand][iView].difFrac;
                    /**/
                }
            }
        }
        s3AerFracImgs[0].img[idx] = pix.lutpars.climFineMode;
        s3AerFracImgs[1].img[idx] = pix.lutpars.mix_frac[0]; // fine mode
        s3AerFracImgs[2].img[idx] = pix.lutpars.mix_frac[1]; // weak of fine
        s3AerFracImgs[3].img[idx] = pix.lutpars.mix_frac[2]; // dust of coarse
        s3FminImg.img[idx] = pix.fmin;
        
        s3SsaImg.img[idx] = pix.ssa[0];
        s3AbsAodImg.img[idx]  = (pix.ssa > 0) ? (1 - pix.ssa[0]) * pix.aod : 0;
        s3DustAodImg.img[idx] = pix.lutpars.mix_frac[2] * (1 - pix.lutpars.mix_frac[0]) * pix.aod;
        s3FmAodImg.img[idx]   = pix.lutpars.mix_frac[0] * pix.aod;
        s3AngstromImg.img[idx] = log(1. / pix.spec_aod_fac[2]) * angWvlLog;
    }
}

void S3NcdfData::computeGeoExtent() {
    int cnrIdx[4] = { 0,
                      (s3LatCnrImgs[0].imgP.height - 1) * s3LatCnrImgs[0].imgP.width,
                      s3LatCnrImgs[0].imgP.height * s3LatCnrImgs[0].imgP.width - 1,
                      s3LatCnrImgs[0].imgP.width - 1
    };

    geoLimits[0] = 999;
    geoLimits[1] = -999;
    geoLimits[2] = 999;
    geoLimits[3] = -999;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            if (geoLimits[0] > s3LatCnrImgs[i].img[cnrIdx[j]]) {
                geoLimits[0] = s3LatCnrImgs[i].img[cnrIdx[j]];
            }
            if (geoLimits[1] < s3LatCnrImgs[i].img[cnrIdx[j]]) {
                geoLimits[1] = s3LatCnrImgs[i].img[cnrIdx[j]];
            }
            if (geoLimits[2] > s3LonCnrImgs[i].img[cnrIdx[j]]) {
                geoLimits[2] = s3LonCnrImgs[i].img[cnrIdx[j]];
            }
            if (geoLimits[3] < s3LonCnrImgs[i].img[cnrIdx[j]]) {
                geoLimits[3] = s3LonCnrImgs[i].img[cnrIdx[j]];
            }
        }
    }
}

double S3NcdfData::getLatMin() {
    if (!isGeoLimitsAvailable) {
        computeGeoExtent();
    }
    return geoLimits[0];
}

double S3NcdfData::getLatMax() {
    if (!isGeoLimitsAvailable) {
        computeGeoExtent();
    }
    return geoLimits[1];
}

double S3NcdfData::getLonMin() {
    if (!isGeoLimitsAvailable) {
        computeGeoExtent();
    }
    return geoLimits[2];
}

double S3NcdfData::getLonMax() {
    if (!isGeoLimitsAvailable) {
        computeGeoExtent();
    }
    return geoLimits[3];
}


//
// private
//
void S3NcdfData::readImgBinned(S3BasicImage<short>* s3Img, const std::string& ncdfName, const std::string varName, 
                               const NcdfImageType& imgType){
    
    //cout << varName << endl;
    // open file
    NcFile ncF(ncdfName, NcFile::read);
    s3Img->name = varName;
    
    NcVar imgVar = ncF.getVar(varName);
    if (imgVar.isNull()) {
        //throw exceptions::NcNotVar("var is null", varName.c_str(), 26);
    }

    // get add_offset, scale_factor and _FillValue attriubutes 
    // from image Variable IF available
    getVarAttSafely(s3Img, imgVar);
    
    if (varName.find("radiance") != std::string::npos){
        if (varName[varName.size()-1] == 'n'){
            getBinRadImg(s3Img, pars.s3MD.slstrPInfo.nadirImg0500m, imgVar);
        }
        else if (varName[varName.size()-1] == 'o'){
            getBinRadImg(s3Img, pars.s3MD.slstrPInfo.obliqImg0500m, imgVar);
        }
    }
}

void S3NcdfData::readImgBinned(S3BasicImage<unsigned short>* s3Img, const std::string& ncdfName, const std::string varName, 
                               const NcdfImageType& imgType){
    
    //cout << varName << endl;
    // open file
    NcFile ncF(ncdfName, NcFile::read);
    s3Img->name = varName;    
    
    NcVar imgVar = ncF.getVar(varName);
    if (imgVar.isNull()) {
        //throw exceptions::NcNotVar("var is null", varName.c_str(), 26);
    }

    // flag image has no attributes, so we dont get them here ;)
    s3Img->setFillVal(32);
    
    if ((varName.find("confid") != std::string::npos)
         || (varName.find("cloud") != std::string::npos)){
        if (varName[varName.size()-1] == 'n'){
            getFlagImg(s3Img, pars.s3MD.slstrPInfo.nadirImg0500m, imgVar);
        }
        else if (varName[varName.size()-1] == 'o'){
            getFlagImg(s3Img, pars.s3MD.slstrPInfo.obliqImg0500m, imgVar);
        }
    }
    
}

void S3NcdfData::readRadImg(S3BasicImage<short>* s3Img, const std::string& ncdfName, const std::string varName, 
                               const NcdfImageType& imgType){
    
    //cout << varName << endl;
    // open file
    NcFile ncF(ncdfName, NcFile::read);
    s3Img->name = varName;    
    
    NcVar imgVar = ncF.getVar(varName);
    if (imgVar.isNull()) {
        //throw exceptions::NcNotVar("var is null", varName.c_str(), 26);
    }

    // get add_offset, scale_factor and _FillValue attriubutes 
    // from image Variable IF available
    getVarAttSafely(s3Img, imgVar);
    
    if (varName[varName.size()-1] == 'n'){
        getRadImg(s3Img, pars.s3MD.slstrPInfo.nadirImg0500m, imgVar);
    }
    else if (varName[varName.size()-1] == 'o'){
        getRadImg(s3Img, pars.s3MD.slstrPInfo.obliqImg0500m, imgVar);
    }
    
}

void S3NcdfData::readSCloudS3SU(S3BasicImage<int>* s3Img, const std::string& ncdfName, const std::string varName){    
    //cout << varName << endl;
    // open file
    NcFile ncF(ncdfName, NcFile::read);
    s3Img->name = varName;    
    
    NcVar imgVar = ncF.getVar(varName);
    if (imgVar.isNull()) {
        //throw exceptions::NcNotVar("var is null", varName.c_str(), 26);
    }
    imgVar.getVar(s3Img->img);
}

void S3NcdfData::readImgBinned(S3BasicImage<float>* s3Img, const std::string& ncdfName, const std::string varName, 
                               const NcdfImageType& imgType){
    
    //cout << varName << endl;
    // open file
    NcFile ncF(ncdfName, NcFile::read);
    s3Img->name = varName;    
    
    NcVar imgVar = ncF.getVar(varName);
    if (imgVar.isNull()) {
        //throw exceptions::NcNotVar("var is null", varName.c_str(), 26);
    }

    // get add_offset, scale_factor and _FillValue attriubutes 
    // from image Variable IF available
    getVarAttSafely(s3Img, imgVar);
    
    if (varName.find("surface_pressure") != std::string::npos){
        getBinPresImg(s3Img, pars.s3MD.slstrPInfo.nadirTpgImg, imgVar);
    }

}

void S3NcdfData::readImgBinned(S3BasicImage<double>* s3Img, const std::string& ncdfName, const std::string varName, 
                               const NcdfImageType& imgType){
    
    // cout << varName << endl;
    // open file
    NcFile ncF(ncdfName, NcFile::read);
    s3Img->name = varName;    
    
    NcVar imgVar = ncF.getVar(varName);
    if (imgVar.isNull()) {
        //throw exceptions::NcNotVar("var is null", varName.c_str(), 26);
    }

    // get add_offset, scale_factor and _FillValue attriubutes 
    // from image Variable IF available
    // getVarAttSafely(s3Img, imgVar); is done in getRadImg for tmp img
    
    if (varName.find("radiance") != std::string::npos){
        if (varName[varName.size()-1] == 'n'){
            getBinRadImgF(s3Img, pars.s3MD.slstrPInfo.nadirImg0500m, imgVar);
        }
        else if (varName[varName.size()-1] == 'o'){
            getBinRadImgF(s3Img, pars.s3MD.slstrPInfo.obliqImg0500m, imgVar);
        }
    }

    if (ncdfName.find("geodetic") != std::string::npos){
        if (varName[varName.size()-1] == 'n'){
            //getBinGeoLocImg(s3Img, pars.s3MD.slstrPInfo.nadirImg0500m, imgVar);
            getCtrGeoLocImg(s3Img, pars.s3MD.slstrPInfo.nadirImg0500m, imgVar);
            getCnrGeoLocImg(pars.s3MD.slstrPInfo.nadirImg0500m, imgVar);
        }
        else if (varName[varName.size()-1] == 'o'){
            //getBinGeoLocImg(s3Img, pars.s3MD.slstrPInfo.obliqImg0500m, imgVar);
            getCtrGeoLocImg(s3Img, pars.s3MD.slstrPInfo.obliqImg0500m, imgVar);
        }
    }
    
    if (ncdfName.find("geometry") != std::string::npos){
        if (varName[varName.size()-1] == 'n'){
            getBinGeomImg(s3Img, pars.s3MD.slstrPInfo.nadirTpgImg, imgVar);
        }
        else if (varName[varName.size()-1] == 'o'){
            getBinGeomImg(s3Img, pars.s3MD.slstrPInfo.obliqTpgImg, imgVar);
        }
    }
    
}

void S3NcdfData::getFlagImg(S3BasicImage<unsigned short>* s3Img, const ImageProperties& imgProp, const NcVar& imgVar){
    char isObliq = (imgProp == pars.s3MD.slstrPInfo.obliqImg0500m) ? 1 : 0;
    // if no binning or shifting is required simply return the img
    if (!s3Img->imgP.isBinned && !isObliq){
        imgVar.getVar(s3Img->img);
        return;
    }
    
    const ImageProperties& origNadirProp = pars.s3MD.slstrPInfo.nadirImg0500m;
    //const int dx = origNadirProp.xOff - ( imgProp.xOff - ((isObliq)?( -1 ) : 0) );
    //const int dy = origNadirProp.yOff - ( imgProp.yOff - ((isObliq)?(  1 ) : 0) );
    const int dx = origNadirProp.xOff - ( imgProp.xOff - offCorr[0] );
    const int dy = origNadirProp.yOff - ( imgProp.yOff - offCorr[1] );

    S3BasicImage<unsigned short> tmp(imgProp);
    imgVar.getVar(tmp.img);
    int i2, j2;
    for (int j1 = 0; j1 < s3Img->imgP.height; j1++){
        for (int i1 = 0; i1 < s3Img->imgP.width; i1++){
            i2 = i1 - dx;
            j2 = j1 - dy;
            if (j2 >= 0 && j2 < tmp.imgP.height && i2 >= 0 && i2 < tmp.imgP.width){
                s3Img->img[ j1 * s3Img->imgP.width + i1 ] = tmp.img[ j2 * tmp.imgP.width + i2 ];
            }
            else {
                s3Img->img[ j1 * s3Img->imgP.width + i1 ] = s3Img->noData;
            }
        }
    }
}

/**
 * get radiance image without binning
 * (implemented for testing valid range of radiance pixels before binning)
 * @param s3Img
 * @param imgProp
 * @param imgVar
 */
void S3NcdfData::getRadImg(S3BasicImage<short>* s3Img, const ImageProperties& imgProp, const NcVar& imgVar){
    char isObliq = (imgProp == pars.s3MD.slstrPInfo.obliqImg0500m) ? 1 : 0;
    // if no binning or shifting is required simply return the img
    if (!s3Img->imgP.isBinned && !isObliq){
        imgVar.getVar(s3Img->img);
        return;
    }
    
    const ImageProperties& origNadirProp = pars.s3MD.slstrPInfo.nadirImg0500m;
    const int dx = origNadirProp.xOff - ( imgProp.xOff - offCorr[0] );
    const int dy = origNadirProp.yOff - ( imgProp.yOff - offCorr[1] );

    S3BasicImage<short> tmp(imgProp);
    imgVar.getVar(tmp.img);
    int i2, j2;
    for (int j1 = 0; j1 < s3Img->imgP.height; j1++){
        for (int i1 = 0; i1 < s3Img->imgP.width; i1++){
            i2 = i1 - dx;
            j2 = j1 - dy;
            //if (i1==2482 && j1==190){
            //    fprintf(stderr, "stop\n");
            //} 
            if (j2 >= 0 && j2 < tmp.imgP.height && i2 >= 0 && i2 < tmp.imgP.width){
                s3Img->img[ j1 * s3Img->imgP.width + i1 ] = tmp.img[ j2 * tmp.imgP.width + i2 ];
            }
            else {
                s3Img->img[ j1 * s3Img->imgP.width + i1 ] = s3Img->noData;
            }
        }
    }
}

/**
 * get radiance image binned into short image
 * additionally writes number of pixels, when image is channel S1
 * @param s3Img
 * @param imgProp
 * @param imgVar
 */
void S3NcdfData::getBinRadImg(S3BasicImage<short>* s3Img, const ImageProperties& imgProp, const NcVar& imgVar){
    char isObliq = (imgProp == pars.s3MD.slstrPInfo.obliqImg0500m) ? 1 : 0;
    char isS1Rad = (s3Img->name.find("S1") != std::string::npos) ? 1 : 0;
    const float binThrs = pars.binValidThrs * s3Img->imgP.binSize * s3Img->imgP.binSize;
    // if no binning or shifting is required simply return the img
    if (!s3Img->imgP.isBinned && !isObliq && offCorr[0] == 0 && offCorr[1] == 0){
        imgVar.getVar(s3Img->img);
        int idx;
        bool isValidPixel;
        if(isFlagsImgAvailable){
            for (int j1 = 0; j1 < s3Img->imgP.height; j1++){
                for (int i1 = 0; i1 < s3Img->imgP.width; i1++){
                    idx = j1*s3Img->imgP.width+i1;
                    isValidPixel = s3Img->isValidValue(s3Img->img[idx]);
                    isValidPixel = isValidPixel && (s3ValidImg.img[idx] & (4+1));
                    if (isValidPixel) {
                        flags.img[idx] |= (s3ValidImg.img[idx] & (4+1));
                    }
                    else {
                        flags.img[idx] &= ~(4+1);
                        s3Img->img[idx] = s3Img->noData;
                    }
                }
            }
        }
        return;
    }
    
    S3BasicImage<short> tmp(imgProp);
    imgVar.getVar(tmp.img);
    
    const ImageProperties& origNadirProp = pars.s3MD.slstrPInfo.nadirImg0500m;
    const int dx = origNadirProp.xOff - ( imgProp.xOff - offCorr[0] );
    const int dy = origNadirProp.yOff - ( imgProp.yOff - offCorr[1] );

    const int& winSize = s3Img->imgP.binSize;
    double sumLand, sumOcn;
    int nL, nO, idx2, idx3, i3, j3;
    // i/j 1 : iterate in binned output img
    // i/j 2 : iterate in full nadir original res img 
    // i/j 3 : iterate in smaller oblique image
    for (int j1 = 0; j1 < s3Img->imgP.height; j1++){
        for (int i1 = 0; i1 < s3Img->imgP.width; i1++){
            nL = nO = 0;
            sumLand = sumOcn = 0;
//            if (j1==67 && i1==127){
//                fprintf(stderr, "stop\n");
//            }
            for (int j2 = j1 * winSize; j2 < (j1 + 1) * winSize; j2++){
                for (int i2 = i1 * winSize; i2 < (i1 + 1) * winSize; i2++){
                    // i2, j2 are the actual idx in the nadir 500 image including the respective offset
                    idx2 = j2 * origNadirProp.width + i2;
                    i3 = i2 - dx;
                    j3 = j2 - dy;
                    if (i3 >= 0 && i3 < tmp.imgP.width && j3 >= 0 && j3 < tmp.imgP.height){
                        idx3 = j3 * tmp.imgP.width + i3;
                        //acriTest
                        //if (i3==2431 && j3==656){
                        //    fprintf(stderr, "-- %d\n", tmp.img[idx3]);
                        //}
                        if (s3Img->isValidValue(tmp.img[idx3])) {
                            if (s3ValidImg.img[idx2] & (1 << isObliq)){
                                // clear land (1:nadir; 2:obliq))
                                sumLand += tmp.img[idx3];
                                nL++;
                            }
                            else if (s3ValidImg.img[idx2] & (4 << isObliq)){
                                // clear ocean (4:nadir; 8:obliq))
                                sumOcn += tmp.img[idx3];
                                nO++;
                            }
                        }
                    }
                }
            }
            if (nL > binThrs){
                s3Img->img[j1 * s3Img->imgP.width + i1] = (short)(sumLand / nL + 0.5);
                flags.img[j1 * flags.imgP.width + i1] |= (1 << isObliq);
                if (isS1Rad){
                    s3NPixImg.img[j1 * s3Img->imgP.width + i1] = nL;
                }
            }
            else if (nO > binThrs){
                s3Img->img[j1 * s3Img->imgP.width + i1] = (short)(sumOcn / nO + 0.5);
                flags.img[j1 * flags.imgP.width + i1] |= (4 << isObliq);
                if (isS1Rad){
                    s3NPixImg.img[j1 * s3Img->imgP.width + i1] = nO;
                }
            }
            else {
                s3Img->img[j1 * s3Img->imgP.width + i1] = s3Img->noData;
                flags.img[j1 * flags.imgP.width + i1] &= ~((4+1) << isObliq);
            }
        }
    }
}
/**/

/**
 * get radiance image binned into float image
 * additionally writes number of pixels, when image is channel S1
 * @param s3Img
 * @param imgProp
 * @param imgVar
 */
void S3NcdfData::getBinRadImgF(S3BasicImage<double>* s3Img, const ImageProperties& imgProp, const NcVar& imgVar){
    char isObliq = (imgProp == pars.s3MD.slstrPInfo.obliqImg0500m) ? 1 : 0;
    char isS1Rad = (s3Img->name.find("S1") != std::string::npos) ? 1 : 0;
    const float binThrs = pars.binValidThrs * s3Img->imgP.binSize * s3Img->imgP.binSize;
    // if no binning or shifting is required simply return the img
    /** to be implemented
    if (!s3Img->imgP.isBinned && !isObliq && offCorr[0] == 0 && offCorr[1] == 0){
        imgVar.getVar(s3Img->img);
        int idx;
        bool isValidPixel;
        if(isFlagsImgAvailable){
            for (int j1 = 0; j1 < s3Img->imgP.height; j1++){
                for (int i1 = 0; i1 < s3Img->imgP.width; i1++){
                    idx = j1*s3Img->imgP.width+i1;
                    isValidPixel = s3Img->isValidValue(s3Img->img[idx]);
                    isValidPixel = isValidPixel && (s3ValidImg.img[idx] & (4+1));
                    if (isValidPixel) {
                        flags.img[idx] |= (s3ValidImg.img[idx] & (4+1));
                    }
                    else {
                        flags.img[idx] &= ~(4+1);
                        s3Img->img[idx] = s3Img->noData;
                    }
                }
            }
        }
        return;
    }
    /**/
    
    S3BasicImage<short> tmp(imgProp);
    getVarAttSafely(&tmp, imgVar);
    replaceStringInPlace(s3Img->name, "_a", "_fa");
    imgVar.getVar(tmp.img);
    
    const ImageProperties& origNadirProp = pars.s3MD.slstrPInfo.nadirImg0500m;
    const int dx = origNadirProp.xOff - ( imgProp.xOff - offCorr[0] );
    const int dy = origNadirProp.yOff - ( imgProp.yOff - offCorr[1] );

    const int& winSize = s3Img->imgP.binSize;
    double sumLand, sumOcn;
    int nL, nO, idx2, idx3, i3, j3;
    // i/j 1 : iterate in binned output img
    // i/j 2 : iterate in full nadir original res img 
    // i/j 3 : iterate in smaller oblique image
    for (int j1 = 0; j1 < s3Img->imgP.height; j1++){
        for (int i1 = 0; i1 < s3Img->imgP.width; i1++){
            nL = nO = 0;
            sumLand = sumOcn = 0;
//            if (j1==67 && i1==127){
//                fprintf(stderr, "stop\n");
//            }
            for (int j2 = j1 * winSize; j2 < (j1 + 1) * winSize; j2++){
                for (int i2 = i1 * winSize; i2 < (i1 + 1) * winSize; i2++){
                    // i2, j2 are the actual idx in the nadir 500 image including the respective offset
                    idx2 = j2 * origNadirProp.width + i2;
                    i3 = i2 - dx;
                    j3 = j2 - dy;
                    if (i3 >= 0 && i3 < tmp.imgP.width && j3 >= 0 && j3 < tmp.imgP.height){
                        idx3 = j3 * tmp.imgP.width + i3;
                        //testing
                        //if (i3==1280 && j3==0){
                        //    fprintf(stderr, "-- %d\n", tmp.img[idx3]);
                        //}
                        if (tmp.isValidValue(tmp.img[idx3])) {
                            if ((s3ValidImg.img[idx2] & (1+2)) == (1+2)){
                                // clear land (1:nadir; 2:obliq))
                                sumLand += tmp.img[idx3] * tmp.valScale + tmp.valOffset;
                                nL++;
                            }
                            else if (s3ValidImg.img[idx2] & (4 << isObliq)){
                                // clear ocean (4:nadir; 8:obliq))
                                sumOcn += tmp.img[idx3] * tmp.valScale + tmp.valOffset;
                                nO++;
                            }
                        }
                    }
                }
            }
            if (nL > binThrs){
                s3Img->img[j1 * s3Img->imgP.width + i1] = sumLand / nL;
                flags.img[j1 * flags.imgP.width + i1] |= (1 << isObliq);
                if (isS1Rad){
                    s3NPixImg.img[j1 * s3Img->imgP.width + i1] = nL;
                }
            }
            else if (nO > binThrs){
                s3Img->img[j1 * s3Img->imgP.width + i1] = sumOcn / nO;
                flags.img[j1 * flags.imgP.width + i1] |= (4 << isObliq);
                if (isS1Rad){
                    s3NPixImg.img[j1 * s3Img->imgP.width + i1] = nO;
                }
            }
            else {
                s3Img->img[j1 * s3Img->imgP.width + i1] = s3Img->noData;
                flags.img[j1 * flags.imgP.width + i1] &= ~((4+1) << isObliq);
            }
        }
    }
}

/**
 * get geolocation by binning hi-res lat / lon 
 * @param s3Img
 * @param imgProp
 * @param imgVar
 */
void S3NcdfData::getBinGeoLocImg(S3BasicImage<double>* s3Img, const ImageProperties& imgProp, const NcVar& imgVar){
    char isObliq = (imgProp == pars.s3MD.slstrPInfo.obliqImg0500m) ? 1 : 0;

    S3BasicImage<double> tmp(imgProp);
    imgVar.getVar(tmp.img);
    getVarAttSafely(&tmp, imgVar);

    const ImageProperties& origNadirProp = pars.s3MD.slstrPInfo.nadirImg0500m;
    //const int dx = origNadirProp.xOff - ( imgProp.xOff - ((isObliq)?( -1 ) : 0) );
    //const int dy = origNadirProp.yOff - ( imgProp.yOff - ((isObliq)?(  1 ) : 0) );
    const int dx = origNadirProp.xOff - ( imgProp.xOff - offCorr[0] );
    const int dy = origNadirProp.yOff - ( imgProp.yOff - offCorr[1] );

    const int& winSize = s3Img->imgP.binSize;
    double sum;
    int n, idx2, i3, j3;
    // i/j 1 : iterate in binned output img
    // i/j 2 : iterate in full nadir original res img 
    // i/j 3 : iterate in smaller oblique image
    for (int j1 = 0; j1 < s3Img->imgP.height; j1++){
        for (int i1 = 0; i1 < s3Img->imgP.width; i1++){
            n = 0;
            sum = 0;
            for (int j2 = j1 * winSize; j2 < (j1 + 1) * winSize; j2++){
                for (int i2 = i1 * winSize; i2 < (i1 + 1) * winSize; i2++){
                    // i2, j2 are the actual idx in the nadir 500 image including the respective offset
                    i3 = i2 - dx;
                    j3 = j2 - dy;
                    if (i3 >= 0 && i3 < tmp.imgP.width && j3 >= 0 && j3 < tmp.imgP.height){
                        idx2 = j3 * tmp.imgP.width + i3;
                        if (tmp.isValidValue(tmp.img[idx2])) {
                            sum += tmp.img[idx2] * tmp.valScale + tmp.valOffset;
                            n++;
                        }
                    }
                }
            }
            if (n>0){
                s3Img->img[j1 * s3Img->imgP.width + i1] = sum / n;
            }
            else {
                s3Img->img[j1 * s3Img->imgP.width + i1] = s3Img->noData;
            }
        }
    }
}

/**
 * get geo-location image by reading values of centre pixel of bin
 * @param s3Img
 * @param imgProp
 * @param imgVar
 */
void S3NcdfData::getCtrGeoLocImg(S3BasicImage<double>* s3Img, const ImageProperties& imgProp, const NcVar& imgVar){
    char isObliq = (imgProp == pars.s3MD.slstrPInfo.obliqImg0500m) ? 1 : 0;

    S3BasicImage<double> tmp(imgProp);
    imgVar.getVar(tmp.img);
    getVarAttSafely(&tmp, imgVar);

    const ImageProperties& origNadirProp = pars.s3MD.slstrPInfo.nadirImg0500m;
    const int dx = origNadirProp.xOff - ( imgProp.xOff - offCorr[0] );
    const int dy = origNadirProp.yOff - ( imgProp.yOff - offCorr[1] );

    const int& winSize = s3Img->imgP.binSize;
    const int& xOff = pars.offset;
    const int& yOff = pars.offset;
    double cntrVal;
    int idx, iCntr, jCntr;
    // i/j 1 : iterate in binned output img
    for (int j1 = 0; j1 < s3Img->imgP.height; j1++){
        for (int i1 = 0; i1 < s3Img->imgP.width; i1++){
            // idx of center pixel
            iCntr = i1 * winSize + xOff + dx;
            jCntr = j1 * winSize + yOff + dy;
            idx = jCntr * tmp.imgP.width + iCntr;
            if (tmp.isValidValue(tmp.img[idx])){
                cntrVal = tmp.img[idx] * tmp.valScale + tmp.valOffset;
                s3Img->img[j1 * s3Img->imgP.width + i1] = cntrVal;
            }
            else {
                s3Img->img[j1 * s3Img->imgP.width + i1] = s3Img->noData;
            }
        }
    }
}

/**
 * get geo-location image for corner pixels of bin only for NADIR
 * @param imgProp
 * @param imgVar
 */
void S3NcdfData::getCnrGeoLocImg(const ImageProperties& imgProp, const NcVar& imgVar){
    //determine if imgVar is lat(1) or lon(2) image
    S3BasicImage<double>* s3CnrImgs;
    char latLonId = 0;
    if (imgVar.getName() == pars.s3MD.LAT_NAME[0]) {
        latLonId = 1;
        s3CnrImgs = s3LatCnrImgs;
    } 
    else if (imgVar.getName() == pars.s3MD.LON_NAME[0]) {
        latLonId = 2;
        s3CnrImgs = s3LonCnrImgs;
    }
    if (latLonId == 0) return;
    
    S3BasicImage<double> tmp(imgProp);
    imgVar.getVar(tmp.img);
    getVarAttSafely(&tmp, imgVar);

    const ImageProperties& origNadirProp = pars.s3MD.slstrPInfo.nadirImg0500m;
    const int dx = origNadirProp.xOff - ( imgProp.xOff - offCorr[0] );
    const int dy = origNadirProp.yOff - ( imgProp.yOff - offCorr[1] );

    const int& winSize = s3CnrImgs[0].imgP.binSize;
    double crnrVal;
    int idx, iCrnr, jCrnr;
    int iCnrOff[4] = {0, 0, winSize-1, winSize-1};
    int jCnrOff[4] = {0, winSize-1, winSize-1, 0};
    // i/j 1 : iterate in binned output img
    for (int j1 = 0; j1 < s3CnrImgs[0].imgP.height; j1++){
        for (int i1 = 0; i1 < s3CnrImgs[0].imgP.width; i1++){
            for (int k = 0; k < 4; k++) {
                // idx of NW corner pixel
                iCrnr = i1 * winSize + iCnrOff[k] + dx;
                jCrnr = j1 * winSize + jCnrOff[k] + dy;
                idx = jCrnr * tmp.imgP.width + iCrnr;
                if (tmp.isValidValue(tmp.img[idx])){
                    crnrVal = tmp.img[idx] * tmp.valScale + tmp.valOffset;
                    s3CnrImgs[k].img[j1 * s3CnrImgs[k].imgP.width + i1] = crnrVal;
                }
                else {
                    s3CnrImgs[k].img[j1 * s3CnrImgs[k].imgP.width + i1] = s3CnrImgs[k].noData;
                }
            }
        }
    }
}

/**
 * get geometry image by reading values for centre of bin pixel
 * (obsolete code for binning still included but not used)
 * @param s3Img
 * @param imgProp
 * @param imgVar
 */
void S3NcdfData::getBinGeomImg(S3BasicImage<double>* s3Img, const ImageProperties& imgProp, const NcVar& imgVar){
    S3BasicImage<double> tmp(imgProp);
    imgVar.getVar(tmp.img);
    getVarAttSafely(&tmp, imgVar);
    //tmp.copyVarAttFrom(*s3Img);
    
    if (imgVar.getName() == S3MetaData::VAA_NAME[0]){
        corrTpg(tmp);
    }
    
    const ImageProperties& origNadirProp = pars.s3MD.slstrPInfo.nadirImg0500m;
    bool disCont = false;//(s3Img->name.find("azimuth") != std::string::npos);
    const int& winSize = s3Img->imgP.binSize;
    const int offs = winSize / 2;
    double val, sum, y2, x2;
    int n;
    // i/j 1 : iterate in binned output img
    // i/j 2 : iterate in full nadir original res img 
    for (int j1 = 0; j1 < s3Img->imgP.height; j1++){
        for (int i1 = 0; i1 < s3Img->imgP.width; i1++){
            n = 0;
            sum = 0;
            //for (int j2 = j1 * winSize; j2 < (j1 + 1) * winSize; j2++){
            for (int j2 = j1 * winSize + offs; j2 < (j1 + 1) * winSize; j2+=winSize){
                y2 = (j2 - origNadirProp.yOff ) * (double)(origNadirProp.yRes) / tmp.imgP.yRes + tmp.imgP.yOff;
                //for (int i2 = i1 * winSize; i2 < (i1 + 1) * winSize; i2++){
                for (int i2 = i1 * winSize + offs; i2 < (i1 + 1) * winSize; i2+=winSize){
                    // i2, j2 are the actual idx in the nadir 500 image including the respective offset
                    // x2, y2 are the intermediate interpolated indices at tpg resolution
                    x2 = (i2 - origNadirProp.xOff ) * (double)(origNadirProp.xRes) / tmp.imgP.xRes + tmp.imgP.xOff;
                    //val = interpol_2d_img(tmp.img, x2, y2, tmp.imgP.width, tmp.imgP.height, disCont);
                    if (x2 >= 0 && y2 >= 0 && x2 < tmp.imgP.width && y2 < tmp.imgP.height){
                        val = interpol_2d_angImg(tmp.img, x2, y2, tmp.imgP.width, tmp.imgP.height, true);
                        if (tmp.isValidValue(val)) {
                            sum += val;
                            n++;
                        }
                    }
                }
            }
            if (n>0){
                s3Img->img[j1 * s3Img->imgP.width + i1] = sum / n;
            }
            else {
                s3Img->img[j1 * s3Img->imgP.width + i1] = s3Img->noData;
            }
        }
    }
}

/**
 * get pressure image from central bin values
 * @param s3Img
 * @param imgProp
 * @param imgVar
 */
void S3NcdfData::getBinPresImg(S3BasicImage<float>* s3Img, const ImageProperties& imgProp, const NcVar& imgVar){
    S3BasicImage<float> tmp(imgProp);
    imgVar.getVar(tmp.img);
    tmp.copyVarAttFrom(*s3Img);
    
    const ImageProperties& origNadirProp = pars.s3MD.slstrPInfo.nadirImg0500m;
    const int& winSize = s3Img->imgP.binSize;
    const int offs = winSize / 2;
    double val, sum, y2, x2;
    int n;
    // i/j 1 : iterate in binned output img
    // i/j 2 : iterate in full nadir original res img 
    for (int j1 = 0; j1 < s3Img->imgP.height; j1++){
        for (int i1 = 0; i1 < s3Img->imgP.width; i1++){
            n = 0;
            sum = 0;
            //for (int j2 = j1 * winSize; j2 < (j1 + 1) * winSize; j2++){
            for (int j2 = j1 * winSize + offs; j2 < (j1 + 1) * winSize; j2+=winSize){
                y2 = (j2 - origNadirProp.yOff ) * (double)(origNadirProp.yRes) / tmp.imgP.yRes + tmp.imgP.yOff;
                //for (int i2 = i1 * winSize; i2 < (i1 + 1) * winSize; i2++){
                for (int i2 = i1 * winSize + offs; i2 < (i1 + 1) * winSize; i2+=winSize){
                    // i2, j2 are the actual idx in the nadir 500 image including the respective offset
                    // x2, y2 are the intermediate interpolated indices at tpg resolution
                    x2 = (i2 - origNadirProp.xOff ) * (double)(origNadirProp.xRes) / tmp.imgP.xRes + tmp.imgP.xOff;
                    //val = interpol_2d_img(tmp.img, x2, y2, tmp.imgP.width, tmp.imgP.height, disCont);
                    if (x2 >= 0 && y2 >= 0 && x2 < tmp.imgP.width && y2 < tmp.imgP.height){
                        val = interpol_2d_img(tmp.img, x2, y2, tmp.imgP.width, tmp.imgP.height, false);
                        if (s3Img->isValidValue(val)) {
                            sum += val;
                            n++;
                        }
                        if (val<1) {
                            fprintf(stderr, "S3NcdfData::getBinPresImg pressure value below 1\n");
                        }
                    }
                }
            }
            if (n>0){
                s3Img->img[j1 * s3Img->imgP.width + i1] = sum / n;
            }
            else {
                s3Img->img[j1 * s3Img->imgP.width + i1] = s3Img->noData;
            }
        }
    }
}

/**
 * spilt strings (should maybe move to miscUtils)
 * @param s
 * @param c
 * @param v
 */
void S3NcdfData::split(const std::string& s, char c, std::vector<std::string>& v) {
   std::string::size_type i = 0;
   std::string::size_type j = s.find(c);

   while (j != std::string::npos) {
      v.push_back(s.substr(i, j-i));
      i = ++j;
      j = s.find(c, j);

      if (j == std::string::npos)
         v.push_back(s.substr(i, s.length()));
   }
}

void S3NcdfData::readImageProp(S3MetaData* s3md, ImageProperties* imgProp){
    std::string ncdfName = s3md->s3FileName + "/" + S3MetaData::CHANNEL_RAD_NAME[0][0] + ".nc";
    NcFile ncF(ncdfName, NcFile::read);
    
    imgProp->width = (int)(ncF.getDim("columns").getSize());
    imgProp->height = (int)(ncF.getDim("rows").getSize());
    imgProp->nPix = imgProp->width * imgProp->height;
    
    ncF.getAtt("track_offset").getValues(&s3md->slstrPInfo.ncdfXoff);
    ncF.getAtt("start_offset").getValues(&s3md->slstrPInfo.ncdfYoff);
    imgProp->xOff = s3md->slstrPInfo.ncdfXoff;
    imgProp->yOff = s3md->slstrPInfo.ncdfYoff;
    
    std::string s;
    ncF.getAtt("resolution").getValues(s);
    std::vector<std::string> resVec;
    split(s, ' ', resVec);
    imgProp->xRes = StringToNumber<int>(resVec[1]);
    imgProp->yRes = StringToNumber<int>(resVec[2]);
}

void S3NcdfData::readImageProp(const NcFile& ncF, ImageProperties* imgProp){
    imgProp->width = (int)(ncF.getDim("columns").getSize());
    imgProp->height = (int)(ncF.getDim("rows").getSize());
    imgProp->nPix = imgProp->width * imgProp->height;
    
    ncF.getAtt("track_offset").getValues(&imgProp->xOff);
    ncF.getAtt("start_offset").getValues(&imgProp->yOff);
    
    std::string s;
    ncF.getAtt("resolution").getValues(s);
    std::vector<std::string> resVec;
    split(s, ' ', resVec);
    imgProp->xRes = StringToNumber<int>(resVec[1]);
    imgProp->yRes = StringToNumber<int>(resVec[2]);
}

/**
 * get add_offset, scale_factor and _FillValue attriubutes from Variable IF available
 * @param s3Img
 * @param imgVar
 */
void S3NcdfData::getVarAttSafely(S3BasicImage<short>* s3Img, const NcVar& imgVar){
    if (hasAtt(imgVar, std::string("add_offset"))) {
        s3Img->hasOffset = true;
        imgVar.getAtt("add_offset").getValues(&s3Img->valOffset);
    }
    else {
        s3Img->hasOffset = false;
        s3Img->valOffset = 0;
    }

    if (hasAtt(imgVar, std::string("scale_factor"))) {
        s3Img->hasScale = true;
        imgVar.getAtt("scale_factor").getValues(&s3Img->valScale);
    }
    else {
        s3Img->hasScale = false;
        s3Img->valScale = 1;
    }

    if (hasAtt(imgVar, std::string("_FillValue"))) {
        s3Img->hasNoData = true;
        imgVar.getAtt("_FillValue").getValues(&s3Img->noData);
    }
    else {
        s3Img->hasNoData = true;
        s3Img->noData = std::numeric_limits<short>::quiet_NaN();
    }
}

/**
 * get add_offset, scale_factor and _FillValue attriubutes from Variable IF available
 * @param s3Img
 * @param imgVar
 */
void S3NcdfData::getVarAttSafely(S3BasicImage<float>* s3Img, const NcVar& imgVar){
    if (hasAtt(imgVar, std::string("add_offset"))) {
        s3Img->hasOffset = true;
        imgVar.getAtt("add_offset").getValues(&s3Img->valOffset);
    }
    else {
        s3Img->hasOffset = false;
        s3Img->valOffset = 0.;
    }

    if (hasAtt(imgVar, std::string("scale_factor"))) {
        s3Img->hasScale = true;
        imgVar.getAtt("scale_factor").getValues(&s3Img->valScale);
    }
    else {
        s3Img->hasScale = false;
        s3Img->valScale = 1.;
    }

    if (hasAtt(imgVar, std::string("_FillValue"))) {
        s3Img->hasNoData = true;
        imgVar.getAtt("_FillValue").getValues(&s3Img->noData);
    }
    else {
        s3Img->hasNoData = true;
        s3Img->noData = std::numeric_limits<double>::quiet_NaN();
    }

    if (hasAtt(imgVar, std::string("valid_min"))
        && hasAtt(imgVar, std::string("valid_max"))) {        
        
        imgVar.getAtt("valid_min").getValues(&s3Img->validMin);
        imgVar.getAtt("valid_max").getValues(&s3Img->validMax);
    }
}

/**
 * get add_offset, scale_factor and _FillValue attriubutes from Variable IF available
 * @param s3Img
 * @param imgVar
 */
void S3NcdfData::getVarAttSafely(S3BasicImage<double>* s3Img, const NcVar& imgVar){
    if (hasAtt(imgVar, std::string("add_offset"))) {
        s3Img->hasOffset = true;
        imgVar.getAtt("add_offset").getValues(&s3Img->valOffset);
    }
    else {
        s3Img->hasOffset = false;
        s3Img->valOffset = 0.;
    }

    if (hasAtt(imgVar, std::string("scale_factor"))) {
        s3Img->hasScale = true;
        imgVar.getAtt("scale_factor").getValues(&s3Img->valScale);
    }
    else {
        s3Img->hasScale = false;
        s3Img->valScale = 1.;
    }

    if (hasAtt(imgVar, std::string("_FillValue"))) {
        s3Img->hasNoData = true;
        imgVar.getAtt("_FillValue").getValues(&s3Img->noData);
    }
    else {
        s3Img->hasNoData = true;
        s3Img->noData = std::numeric_limits<double>::quiet_NaN();
    }

    if (hasAtt(imgVar, std::string("valid_min"))
        && hasAtt(imgVar, std::string("valid_max"))) {        
        s3Img->hasMinMax = true;
        imgVar.getAtt("valid_min").getValues(&s3Img->validMin);
        imgVar.getAtt("valid_max").getValues(&s3Img->validMax);
    }
    else {
        s3Img->hasMinMax = false;
    }
}

/**
 * check if ncdf var has attirbute
 * @param var
 * @param attName
 * @return 
 */
bool S3NcdfData::hasAtt(const NcVar& var, const std::string& attName){
    std::map<std::string, NcVarAtt> atts = var.getAtts();
    std::map<std::string, NcVarAtt>::iterator attIt = atts.find(attName);
    if (attIt != atts.end()){
        return true;
    }
    return false;
}

void S3NcdfData::readIrrad(double irrad[][N_SLSTR_BANDS]){
    std::string ncdfName;
    double irval[N_DETECTORS];
    for (int iView = 0; iView < N_SLSTR_VIEWS; iView++){
        for (int iBand = 0; iBand < N_SLSTR_BANDS; iBand++){
            ncdfName = pars.slstrProductDir + "/" + S3MetaData::CHANNEL_QUAL_NAME[iView][iBand] + ".nc";
            NcFile ncF(ncdfName, NcFile::read);
            NcVar irradVar = ncF.getVar(S3MetaData::CHANNEL_IRRAD_NAME[iView][iBand]);
            irradVar.getVar(irval);
            irrad[iView][iBand] = 0;
            for (int i = 0; i < N_DETECTORS; i++){
                irrad[iView][iBand] += irval[i];
            }
            irrad[iView][iBand] /= N_DETECTORS;
            //ncF.close();
        }
    }   
}

/**
 * create hi-res surface type image for nadir and oblique
 * written to s3LandImg
 */
void S3NcdfData::createLandMask(){
    s3LandImg = S3BasicImage<signed char>(pars.s3MD.slstrPInfo.nadirImg0500m);
    s3LandImg.name = "landmask";
    s3LandImg.initImgArray();
    
    // srfcTypeTest mask applied on confid flags
    // flags used:
    //  1 - coastline
    //  2 - ocean
    //  4 - tidal
    //  8 - land
    // 16 - inland water
    // 32 - fill value
    const unsigned short srfcTypeTest = 0b0000000000011011;

    for (int iView = 0; iView < N_SLSTR_VIEWS; iView++){
        //read full resolution flag images
        std::string ncdfName = pars.slstrProductDir + "/" + S3MetaData::FLAGS_NAME[iView] + ".nc";
        S3BasicImage<unsigned short> s3ConfidImg = S3BasicImage<unsigned short>(pars.s3MD.slstrPInfo.nadirImg0500m);
        s3ConfidImg.noData = 0;
        readImgBinned(&s3ConfidImg, ncdfName, S3MetaData::CONFID_NAME[iView], Nadir0500);
        for (int i = 0; i < s3LandImg.imgP.nPix; i++){
            if (s3ConfidImg.img[i] == 32) continue;
            unsigned short mask = (s3ConfidImg.img[i] & (srfcTypeTest));
            //if (mask == 8){
            if ((mask & (8+1))){ //acri test coastline is always land
                // is clear land
                s3LandImg.img[i] |= 1 << iView;
            }
            //else if (mask == 2 || mask == (8+16)){
            else { //acri test
                // is ocean or inland water
                s3LandImg.img[i] |= 4 << iView;
            }
        }
    }
    isLandMaskAvailable = true;
}

/**
 * create hi-res cloud/snow/glint image for nadir and oblique
 * based in L1b confidence flags
 * written to s3CloudImg
 */
void S3NcdfData::createCloudMask(){
    s3CloudImg = S3BasicImage<signed char>(pars.s3MD.slstrPInfo.nadirImg0500m);
    s3CloudImg.name = "cloudmask";
    s3CloudImg.initImgArray();
    const unsigned short cloudTest = 0b0011011111111110; // for testing against subset of basic cloud flags
    const unsigned short confidCloudTest = 0b0111000000000000; //sun_glint, snow, cloud

    for (int iView = 0; iView < N_SLSTR_VIEWS; iView++){
        //read full resolution flag images
        std::string ncdfName = pars.slstrProductDir + "/" + S3MetaData::FLAGS_NAME[iView] + ".nc";
        S3BasicImage<unsigned short> s3CldImg = S3BasicImage<unsigned short>(pars.s3MD.slstrPInfo.nadirImg0500m);
        s3CldImg.noData = 0;
        //readImgBinned(&s3CldImg, ncdfName, s3MetaData.BASIC_CLOUD_NAME[iView], Nadir0500);
        readImgBinned(&s3CldImg, ncdfName, S3MetaData::CONFID_NAME[iView], Nadir0500);
        for (int i = 0; i < s3CloudImg.imgP.nPix; i++){
            //ushort mask = (s3CldImg.img[i] & (cloudTest));
            unsigned short mask = (s3CldImg.img[i] & (confidCloudTest));
            if (mask){
                // is cloud
                s3CloudImg.img[i] |= 1 << iView;
            }
        }
    }
}

/**
 * create hi-res surface type image for nadir and oblique
 * based on Swansea cloud mask image
 * written to s3LandImg
 */
void S3NcdfData::createMyCloudMask(){
    s3CloudImg = S3BasicImage<signed char>(pars.s3MD.slstrPInfo.nadirImg0500m);
    s3CloudImg.name = "cloudmask";
    s3CloudImg.initImgArray();

    for (int iView = 0; iView < N_SLSTR_VIEWS; iView++){
        //read full resolution flag images
        //std::string ncdfName = pars.slstrProductDir + "/" + S3MetaData::S3SU_CLOUD_FNAME + ".nc";
        S3BasicImage<int> s3CldImg = S3BasicImage<int>(pars.s3MD.slstrPInfo.nadirImg0500m);
        s3CldImg.noData = 0;
        readSCloudS3SU(&s3CldImg, pars.sCloudS3SUname, S3MetaData::S3SU_CLOUD_VNAME[iView]);
        for (int i = 0; i < s3CloudImg.imgP.nPix; i++){
            if (s3CldImg.img[i]>0){
                // is cloud
                s3CloudImg.img[i] |= 1 << iView;
            }
        }
    }
}

/**
 * create hi-res valid radiance mask image for nadir and oblique
 * written to s3NanImg
 */
void S3NcdfData::createNanMask(){
    s3NanImg = S3BasicImage<signed char>(pars.s3MD.slstrPInfo.nadirImg0500m);
    s3NanImg.name = "NaNmask";
    s3NanImg.noData = 0;
    s3NanImg.initImgArray();

    for (int iView = 0; iView < N_SLSTR_VIEWS; iView++){
        for (int iBand = 0; iBand < N_SLSTR_BANDS; iBand++){
            //read full resolution flag images
            std::string ncdfName = pars.slstrProductDir + "/" + S3MetaData::CHANNEL_RAD_NAME[iView][iBand] + ".nc";
            std::string varName = S3MetaData::CHANNEL_RAD_NAME[iView][iBand];
            S3BasicImage<short> s3RadImg = S3BasicImage<short>(pars.s3MD.slstrPInfo.nadirImg0500m);
            readRadImg(&s3RadImg, ncdfName, varName, Nadir0500);
            
            for (int i = 0; i < s3NanImg.imgP.nPix; i++){
                if (!s3RadImg.isValidValue(s3RadImg.img[i])){
                    // is NaN
                    s3NanImg.img[i] |= 1 << iView;
                }
            }
        }
    }
}

/**
 * create hi-res valid pixel image for nadir and oblique
 * written to s3ValidImg
 */
void S3NcdfData::createValidMask(){
    s3ValidImg = S3BasicImage<signed char>(pars.s3MD.slstrPInfo.nadirImg0500m);
    s3ValidImg.name = "validImage";
    s3ValidImg.initImgArray();
    /*const ushort s3Coast = 1;
    const ushort s3Ocean = 2;
    const ushort s3Tidal = 4;
    const ushort s3Land  = 8;
    const ushort s3InlandWater = 16;
    const ushort s3Snow  = 8192;
    const ushort s3Cloud = 16384;
    const ushort srfcTypeTest = 0b0111000000011111;*/

    /* old implmentation *
    for (int iView = 0; iView < N_SLSTR_VIEWS; iView++){
        signed char landMask  = 1 << iView;
        signed char ocnMask   = 4 << iView;
        signed char cloudMask = 1 << iView;
        signed char nanMask = 1 << iView;
        
        for (int i = 0; i < s3ValidImg.imgP.nPix; i++){
            if (i==572482){
                fprintf(stderr, "stop\n");
            }
            if ((s3LandImg.img[i] & landMask) && !(s3CloudImg.img[i] & cloudMask) && !(s3NanImg.img[i] & nanMask)) {
                s3ValidImg.img[i] |= landMask;
            }
            else if ((s3LandImg.img[i] & ocnMask) && !(s3CloudImg.img[i] & cloudMask) && !(s3NanImg.img[i] & nanMask)) {
                s3ValidImg.img[i] |= ocnMask;
            }
        }
    }
    /* end old implemention*/
    
    for (int i = 0; i < s3ValidImg.imgP.nPix; i++){

        s3ValidImg.img[i] = 0;

        if ((s3LandImg.img[i] == 3)
             && (s3CloudImg.img[i] == 0)
             && (s3NanImg.img[i] == 0) ){
            
            // both views are clear valid land
            s3ValidImg.img[i] = 3;
        }
        else {
            if ((s3LandImg.img[i] & 4)
                  && ( !(s3CloudImg.img[i] & 1) )
                  && ( !(s3NanImg.img[i]   & 1) ) ){
            
            // nadir view is clear and valid
            s3ValidImg.img[i] |= 4;
            }
            if ((s3LandImg.img[i] & 8)
                  && ( !(s3CloudImg.img[i] & 2) )
                  && ( !(s3NanImg.img[i]   & 2) ) ){
            
            // nadir view is clear and valid
            s3ValidImg.img[i] |= 8;
            }
        }
    }

}

void S3NcdfData::corrTpg(S3BasicImage<double>& tpg){
    int w = tpg.imgP.width;
    int h = tpg.imgP.height;
    double* lp;
    double y;
    int i1,i2;
    for (int iy=0; iy<h; iy++){
        lp = &tpg.img[iy*w];
        i1=0;   while (i1<(w-1) && ((lp[i1]>lp[i1+1]) || std::isnan(lp[i1]))) i1++;
        i2=w-1; while (i2>0 && ((lp[i2]<lp[i2-1]) || std::isnan(lp[i2]))) i2--;
        if (i1<(w-1) || i2>0){
            for (int i=i1+1; i<i2; i++){
                y = intAng(lp[i1], lp[i2], (double)(i-i1)/(double)(i2-i1));
                if (y<0) y+=360;
                if (abs(y-lp[i])>1) lp[i] = y;
            }
        }
    }
}

/**
 * 
 * @param timeImg
 */
void S3NcdfData::computeAvgTimeImg(S3BasicImage<unsigned int>* timeImg) {
    double startSec = parseTimeToSec1970(pars.s3MD.startTime);
    double stopSec = parseTimeToSec1970(pars.s3MD.stopTime);
    double t;
    int width500 = pars.s3MD.slstrPInfo.nadirImg0500m.width;
    int height500 = pars.s3MD.slstrPInfo.nadirImg0500m.height;

    int nTime = width500 * height500;
   
    const ImageProperties& origNadirProp = pars.s3MD.slstrPInfo.nadirImg0500m;
    const ImageProperties& imgProp = pars.s3MD.slstrPInfo.nadirImg0500m;
    const int dx = origNadirProp.xOff - ( imgProp.xOff - offCorr[0] );
    const int dy = origNadirProp.yOff - ( imgProp.yOff - offCorr[1] );

    const int& winSize = timeImg->imgP.binSize;
    const int& xOff = pars.offset;
    const int& yOff = pars.offset;

    int idx, iCntr, jCntr;
    // i/j 1 : iterate in binned output img
    for (int j1 = 0; j1 < timeImg->imgP.height; j1++){
        for (int i1 = 0; i1 < timeImg->imgP.width; i1++){
            // idx of center pixel
            iCntr = i1 * winSize + xOff + dx;
            jCntr = j1 * winSize + yOff + dy;
            idx = jCntr * width500 + iCntr;
            t = startSec + (stopSec - startSec) * idx / nTime;
            timeImg->img[j1 * timeImg->imgP.width + i1] = (unsigned int)(t);
        }
    }
}

