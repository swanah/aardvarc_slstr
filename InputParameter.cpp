/* 
 * File:   InputParameter.cpp
 * Author: akheckel
 * 
 * Created on 23. November 2016, 10:39
 */

#include "InputParameter.hpp"


using std::cout;
using std::endl;

        
InputParameter::InputParameter() {
    parFileName = "";
    //slstrProductDir = "e:/sat/S3A_SLSTR/S3A_SL_1_RBT____20160902T072518_20160902T072818_20160903T162446_0179_008_163_2879_LN2_O_NT_002.SEN3";
    //slstrProductDir = "e:/sat/S3A_SLSTR/S3A_SL_1_RBT____20160608T092412_20160608T092712_20160608T115216_0179_005_093_2159_SVL_O_NR_001.SEN3";
    slstrProductDir = "e:/sat/S3A_SLSTR/20161230/S3A_SL_1_RBT____20161230T105314_20161230T105614_20161230T131740_0180_012_322_2339_SVL_O_NR_002.SEN3";
    aodOutDir = "e:/sat/S3A_SLSTR_AOD/%YYYY%/%MM%";
    atmLutFileName = "e:/model/mkS3Lut/s3Lut.nc";//"e:/projects/S3MPC/submitted/adf/atmosLut.nc";
    ocnLutFileName = "e:/projects/S3MPC/submitted/adf/ocnLut.nc";
    climFileName   = "e:/projects/S3MPC/submitted/adf/aerosolClimatology.nc";
    
    winSize = 9;
    skip = winSize;
    offset = winSize / 2;
    szaLimit = 75;
    binValidThrs = 0.5;
}


