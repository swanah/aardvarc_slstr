/* 
 * File:   InputParameter.cpp
 * Author: akheckel
 * 
 * Created on 23. November 2016, 10:39
 */

#include <cstdlib>
#include <fstream>
#include <stdexcept>
#include "miscUtils.hpp"
#include "InputParameter.hpp"


using std::cout;
using std::cerr;
using std::endl;
//namespace bfs = boost::filesystem;

        
InputParameter::InputParameter(int argc, char** argv) {
    if (argc != 3) {
        throw std::runtime_error("Invalid number of arguments");
    }
    parFileName = argv[1];
    if (!fileExists(parFileName)) {
        cerr << "Parm-File: " << parFileName << endl;
        throw std::runtime_error("Parameter file doesnt exist");
    }

    slstrProductDir = argv[2];
    if (!fileExists(slstrProductDir)) {
        cerr << "SLSTR L1b: " << slstrProductDir << endl;
        throw std::runtime_error("SLSTR L1b dir doesnt exist");
    }
    if (!isDir(slstrProductDir)) {
        cerr << "SLSTR L1b: " << slstrProductDir << endl;
        throw std::runtime_error("SLSTR L1b not a directory");
    }
    s3MD.parseManifest(slstrProductDir);

    //slstrProductDir = "e:/sat/S3A_SLSTR/S3A_SL_1_RBT____20160902T072518_20160902T072818_20160903T162446_0179_008_163_2879_LN2_O_NT_002.SEN3";
    //slstrProductDir = "e:/sat/S3A_SLSTR/S3A_SL_1_RBT____20160608T092412_20160608T092712_20160608T115216_0179_005_093_2159_SVL_O_NR_001.SEN3";
    //slstrProductDir = "e:/sat/S3A_SLSTR/20161230/S3A_SL_1_RBT____20161230T105314_20161230T105614_20161230T131740_0180_012_322_2339_SVL_O_NR_002.SEN3";
    //slstrProductDir = "e:/sat/S3A_SLSTR/testCases/S3A_SL_1_RBT____20170227T003224_20170227T003524_20170329T135202_0180_015_002______LR1_D_NT_001.SEN3";
    //slstrProductDir = "e:/sat/S3A_SLSTR/testCases/S3A_SL_1_RBT____20170227T102019_20170227T102319_20170329T144316_0179_015_008______LR1_D_NT_001.SEN3";
    //slstrProductDir = "e:/sat/S3A_SLSTR/testCases/S3A_SL_1_RBT____20170227T140317_20170227T140617_20170329T151133_0179_015_010______LR1_D_NT_001.SEN3";
    //aodOutDir = "e:/sat/S3A_SLSTR_AOD/%YYYY%/%MM%";
    //atmLutFileName = "e:/model/mkS3Lut/s3Lut.nc";//"e:/projects/S3MPC/submitted/adf/atmosLut.nc";
    //ocnLutFileName = "e:/model/oceanLut/ocnLut.nc";//"e:/projects/S3MPC/submitted/adf/ocnLut.nc";
    //climFileName   = "e:/projects/S3MPC/submitted/adf/aerosolClimatology.nc";
    //winSize = 9;
    
    doGeoSubset = false;
    latLim[0] = latLim[1] = lonLim[0] = lonLim[1] = -999;
    winSize = 0;
    szaLimit = 0;
    
    std::ifstream f(parFileName.c_str());
    std::string line;
    std::vector<std::string> sVec;
    int sVecLen;
    while (std::getline(f, line)) {
        trimLine(line);
        if (line.empty() || line[0] == '*'){
            continue;
        }
        splitLine(line, sVec);
        sVecLen = sVec.size();
        if (sVecLen > 2){
            if (sVec[0] == "AOD_OUT_PATH") {
                aodOutDir = sVec[2];
                parsePlaceholders(aodOutDir);
            }
            if (sVec[0] == "ATMOS_LUT") {
                atmLutFileName = sVec[2];
                parsePlaceholders(atmLutFileName);
            }
            if (sVec[0] == "OCEAN_LUT") {
                ocnLutFileName = sVec[2];
                parsePlaceholders(ocnLutFileName);
            }
            if (sVec[0] == "AER_CLIMA") {
                climFileName = sVec[2];
                parsePlaceholders(climFileName);
            }
            if (sVec[0] == "CLOUD_FLAGS") useSCloudS3SU = (sVec[2] == "SCLOUDS3SU");
            if (sVec[0] == "SCLOUDS3SU_PATH") {
                sCloudS3SUpath = sVec[2];
                parsePlaceholders(sCloudS3SUpath);
            }
            if (sVec[0] == "WINSIZE") winSize = std::strtod(sVec[2].c_str(), NULL);
            if (sVec[0] == "SZA_LIMIT") szaLimit = std::strtod(sVec[2].c_str(), NULL);
            if (sVec[0] == "BIN_VALID_THRS") binValidThrs = std::strtod(sVec[2].c_str(), NULL);
            if (sVec[0] == "DO_GEOSUBSET") doGeoSubset = (sVec[2] == "TRUE");
        }
        if (sVecLen > 3){
            if (sVec[0] == "LAT_MIN_MAX") {
                latLim[0] = std::strtod(sVec[2].c_str(), NULL);
                latLim[1] = std::strtod(sVec[3].c_str(), NULL);
            }
            if (sVec[0] == "LON_MIN_MAX") {
                lonLim[0] = std::strtod(sVec[2].c_str(), NULL);
                lonLim[1] = std::strtod(sVec[3].c_str(), NULL);
            }
        }
        
    }
    
    ensurePathExists(aodOutDir);
    createAodName();
    if (useSCloudS3SU){
        ensurePathExists(sCloudS3SUpath);
        createCldName();
        if (!fileExists(sCloudS3SUname)){
            cerr << "Cloud-File: " << sCloudS3SUname << endl;
            throw std::runtime_error("sCloudS3SU: cloud mask file doesnt exist");
        }
    }

    if (winSize < 1 || winSize > 100){
        cerr << "winSize: " << winSize << endl;
        throw std::runtime_error("winSize not in [1..100]");
    }
    skip = winSize;
    offset = winSize / 2;

    if (szaLimit < 1 || szaLimit > 90){
        cerr << "szaLimit: " << szaLimit << endl;
        throw std::runtime_error("szaLimit not in [1..90]");
    }

    if (binValidThrs <= 0 || binValidThrs >= 1){
        cerr << "binValidThr: " << binValidThrs << endl;
        throw std::runtime_error("binValidThr not in (0..1)");
    }
    
    if (doGeoSubset) {
        if (latLim[0] < -90 || latLim[0] > 90 || latLim[1] < -90 || latLim[1] > 90 || latLim[0] >= latLim[1]){
            cerr << "latMinMax: [" << latLim[0] << ".." << latLim[1] << "]" << endl;
            throw std::runtime_error("latLim not in [-90..90] or latMin >= latMax");
        }
        if (lonLim[0] < -180 || lonLim[0] > 180 || lonLim[1] < -180 || lonLim[1] > 180 || lonLim[0] >= lonLim[1]){
            cerr << "lonMinMax: [" << lonLim[0] << ".." << lonLim[1] << "]" << endl;
            throw std::runtime_error("lonLim not in [-180..180] or lonMin >= lonMax");
        }
    }
    
}

void InputParameter::parsePlaceholders(std::string& s) {
    std::string yyyy = s3MD.startTime.substr(0, 4);
    std::string mm = s3MD.startTime.substr(5, 2);
    std::string dd = s3MD.startTime.substr(8, 2);

    replaceStringInPlace(s, "%YYYY", yyyy);
    replaceStringInPlace(s, "%MM", mm);
    replaceStringInPlace(s, "%DD", dd);
}

void InputParameter::ensurePathExists(std::string& path) {
    bool pExists = fileExists(path);
    if (!pExists){
        if (createDirsUnix(path)){
            cerr << path << " created!" << endl;
        }
        else {
            cerr << path << endl;
            throw std::runtime_error("Directory couldnt be created");
        }
    }
}

void InputParameter::createAodName() {
    aodOutName = s3MD.productName;
    int extIdx = aodOutName.find_last_of(".");
    int extLen = aodOutName.length() - extIdx;
    aodOutName.replace(extIdx, extLen, "_aod.nc");
    aodOutName = aodOutDir + "/" + aodOutName;
}

void InputParameter::createCldName() {
    sCloudS3SUname = s3MD.productName;
    int extIdx = sCloudS3SUname.find_last_of(".");
    int extLen = sCloudS3SUname.length() - extIdx;
    sCloudS3SUname.replace(extIdx, extLen, "_cld.nc");
    sCloudS3SUname = sCloudS3SUpath + "/" + sCloudS3SUname;
}
