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
#include <math.h>
#include <netcdf>
#include <sys/_intsup.h>
#include "Images.hpp"
#include "Interpolation.hpp"
#include "miscUtils.hpp"
#include "S3NcdfData.hpp"

using std::cout;
using std::cerr;
using std::endl;
using namespace netCDF;

//public

S3NcdfData::S3NcdfData(const InputParameter& inPar){
    s3DataDir = inPar.slstrProductDir;
    s3MetaData.parseManifest(s3DataDir);
    setAodDataDir(inPar, s3MetaData);
}

S3NcdfData::~S3NcdfData() {
}

/*void S3NcdfData::readNcdf() {
    
    ImageProperties imgProp;
    readImageProp(&s3MetaData, &imgProp);
    readIrrad(s3Irrad);
    
    // reading nadir and oblique radiance images
    NcdfImageType imgType;
    std::string ncdfName;
    for (int iView = 0; iView < N_SLSTR_VIEWS; iView++){
        imgType = (iView == 0) ? Nadir0500 : Obliq0500;
        for (int iBand = 0; iBand < N_SLSTR_BANDS; iBand++){
            s3RadImgs[iView][iBand] = S3BasicImage<short>(imgProp);
            ncdfName = s3MetaData.s3FileName + "/" + s3MetaData.CHANNEL_RAD_NAME[iView][iBand] + ".nc";
            readImg(&s3RadImgs[iView][iBand], ncdfName, s3MetaData.CHANNEL_RAD_NAME[iView][iBand], imgType);
        }

        ncdfName = s3MetaData.s3FileName + "/" + s3MetaData.GEODETIC_NAME[iView] + ".nc";
        s3LatImgs[iView] = S3BasicImage<double>(imgProp);
        readImg(&s3LatImgs[iView], ncdfName, s3MetaData.LAT_NAME[iView], imgType);

        s3LonImgs[iView] = S3BasicImage<double>(imgProp);
        readImg(&s3LonImgs[iView], ncdfName, s3MetaData.LON_NAME[iView], imgType);
        
        imgType = (iView == 0) ? NadirTpg : ObliqTpg;
        
        ncdfName = s3MetaData.s3FileName + "/" + s3MetaData.GEOMETRY_NAME[iView] + ".nc";
        s3SzaImgs[iView] = S3BasicImage<double>(imgProp);
        readImg(&s3SzaImgs[iView], ncdfName, s3MetaData.SZA_NAME[iView], imgType);
        
        s3SaaImgs[iView] = S3BasicImage<double>(imgProp);
        readImg(&s3SaaImgs[iView], ncdfName, s3MetaData.SAA_NAME[iView], imgType);
        
        s3VzaImgs[iView] = S3BasicImage<double>(imgProp);
        readImg(&s3VzaImgs[iView], ncdfName, s3MetaData.VZA_NAME[iView], imgType);
        
        s3VaaImgs[iView] = S3BasicImage<double>(imgProp);
        readImg(&s3VaaImgs[iView], ncdfName, s3MetaData.VAA_NAME[iView], imgType);
    }

    imgType = NadirTpg;

    ncdfName = s3MetaData.s3FileName + "/" + s3MetaData.GEODETIC_TPG_NAME + ".nc";
    s3TpgLatImg = S3BasicImage<double>(imgProp);
    readImg(&s3TpgLatImg, ncdfName, s3MetaData.LAT_TPG_NAME, imgType);

    s3TpgLonImg = S3BasicImage<double>(imgProp);
    readImg(&s3TpgLonImg, ncdfName, s3MetaData.LON_TPG_NAME, imgType);
        
}*/

void S3NcdfData::readNcdf(const ImageProperties& outImgProp) {
    
    readIrrad(s3Irrad);
    createLandMask();
    createCloudMask();
    createValidMask();
    flags = S3BasicImage<short>(outImgProp);
    flags.name = "aod_flags";
    flags.initImgArray();

    // reading nadir and oblique radiance images
    NcdfImageType imgType;
    std::string ncdfName;
    for (int iView = 0; iView < N_SLSTR_VIEWS; iView++){
        imgType = (iView == 0) ? Nadir0500 : Obliq0500;

        //read img and bin to output image properties
        for (int iBand = 0; iBand < N_SLSTR_BANDS; iBand++){
            s3RadImgs[iView][iBand] = S3BasicImage<short>(outImgProp);
            s3RadImgs[iView][iBand].setValidLimits(0, 10000);
            ncdfName = s3MetaData.s3FileName + "/" + s3MetaData.CHANNEL_RAD_NAME[iView][iBand] + ".nc";
            readImgBinned(&s3RadImgs[iView][iBand], ncdfName, s3MetaData.CHANNEL_RAD_NAME[iView][iBand], imgType);
        }

        ncdfName = s3MetaData.s3FileName + "/" + s3MetaData.GEODETIC_NAME[iView] + ".nc";
        s3LatImgs[iView] = S3BasicImage<double>(outImgProp);
        s3LatImgs[iView].setValidLimits(-90., 90.);
        readImgBinned(&s3LatImgs[iView], ncdfName, s3MetaData.LAT_NAME[iView], imgType);

        s3LonImgs[iView] = S3BasicImage<double>(outImgProp);
        s3LonImgs[iView].setValidLimits(-180., 180.);
        readImgBinned(&s3LonImgs[iView], ncdfName, s3MetaData.LON_NAME[iView], imgType);

        imgType = (iView == 0) ? NadirTpg : ObliqTpg;
        
        ncdfName = s3MetaData.s3FileName + "/" + s3MetaData.GEOMETRY_NAME[iView] + ".nc";
        s3SzaImgs[iView] = S3BasicImage<double>(outImgProp);
        s3SzaImgs[iView].setValidLimits(0., 90.);
        readImgBinned(&s3SzaImgs[iView], ncdfName, s3MetaData.SZA_NAME[iView], imgType);
        
        s3SaaImgs[iView] = S3BasicImage<double>(outImgProp);
        s3SaaImgs[iView].setValidLimits(0., 360.);
        readImgBinned(&s3SaaImgs[iView], ncdfName, s3MetaData.SAA_NAME[iView], imgType);
        
        s3VzaImgs[iView] = S3BasicImage<double>(outImgProp);
        s3SzaImgs[iView].setValidLimits(0., 90.);
        readImgBinned(&s3VzaImgs[iView], ncdfName, s3MetaData.VZA_NAME[iView], imgType);
        
        s3VaaImgs[iView] = S3BasicImage<double>(outImgProp);
        s3VaaImgs[iView].setValidLimits(0., 360.);
        readImgBinned(&s3VaaImgs[iView], ncdfName, s3MetaData.VAA_NAME[iView], imgType);
    }
}

void S3NcdfData::convRad2Refl(){
    double rad;
    for (int iView = 0; iView < N_SLSTR_VIEWS; iView++){
        for (int iBand = 0; iBand < N_SLSTR_BANDS; iBand++){
            double origScale = s3RadImgs[iView][iBand].valScale;
            double origOffset = s3RadImgs[iView][iBand].valOffset;
            short  origNoData = s3RadImgs[iView][iBand].noData;
            s3RadImgs[iView][iBand].valScale = 1e-4;
            s3RadImgs[iView][iBand].valOffset = 0;
            s3RadImgs[iView][iBand].noData = -10000;
            for (int i = 0; i < s3RadImgs[iView][iBand].imgP.nPix; i++){
                if (s3RadImgs[iView][iBand].img[i] != origNoData){
                    rad = (double)(s3RadImgs[iView][iBand].img[i]) * origScale + origOffset;
                    rad *= M_PI / (s3Irrad[iView][iBand] * cos(s3SzaImgs[iView].img[i] * M_PI / 180));
                    s3RadImgs[iView][iBand].img[i] = (rad - s3RadImgs[iView][iBand].valOffset) / s3RadImgs[iView][iBand].valScale;
                }
                else {
                    s3RadImgs[iView][iBand].img[i] = s3RadImgs[iView][iBand].noData;
                }
            }
        }
    }
}

void S3NcdfData::verifyInput(){
    char valid;
    int iBand;
    
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
}


// private

void S3NcdfData::setAodDataDir(const InputParameter& inPar, const S3MetaData& s3MetaData){
    std::string yyyy = s3MetaData.startTime.substr(0, 4);
    std::string mm = s3MetaData.startTime.substr(5, 2);
    std::string dd = s3MetaData.startTime.substr(8, 2);

    aodDataDir = inPar.aodOutDir;
    replaceStringInPlace(aodDataDir, "%YYYY%", yyyy);
    replaceStringInPlace(aodDataDir, "%MM%", mm);
    replaceStringInPlace(aodDataDir, "%DD%", dd);
    
    aodOutName = s3MetaData.productName;
    int extIdx = aodOutName.find_last_of(".") + 1;
    int extLen = aodOutName.length() - extIdx;
    aodOutName.replace(extIdx, extLen, "nc");
    aodOutName = aodDataDir + "/" + aodOutName;
}

/*void S3NcdfData::readImg(S3BasicImage<short>* s3Img, const std::string& ncdfName, const std::string& varName, const NcdfImageType& imgType){
    cout << varName << endl;

    // open file
    NcFile ncF(ncdfName, NcFile::read);
    s3Img->name = varName;    
    
    ImageProperties imgProp;
    readImageProp(ncF, &imgProp);
    
    NcVar imgVar = ncF.getVar(varName);
    if (imgVar.isNull()) {
        throw exceptions::NcNotVar("var is null", varName.c_str(), 26);
    }
    
    // get add_offset, scale_factor and _FillValue attriubutes 
    // from image Variable IF available
    getVarAttSafely(s3Img, imgVar);
    
    if (imgType == Nadir0500) {
        imgVar.getVar(s3Img->img);
    }
    else if (imgType == Obliq0500){
        getImgShifted(s3Img, imgProp, imgVar);
    }
    else if (imgType == NadirTpg || imgType == ObliqTpg){
        getImgScaled(s3Img, imgProp, imgVar);
    }
    
}

void S3NcdfData::readImg(S3BasicImage<double>* s3Img, const std::string& ncdfName, const std::string& varName, const NcdfImageType& imgType){
    cout << varName << endl;

    // open file
    NcFile ncF(ncdfName, NcFile::read);
    s3Img->name = varName;    
    
    ImageProperties imgProp;
    readImageProp(ncF, &imgProp);
    
    NcVar imgVar = ncF.getVar(varName);
    if (imgVar.isNull()) {
        throw exceptions::NcNotVar("var is null", varName.c_str(), 26);
    }

    // get add_offset, scale_factor and _FillValue attriubutes 
    // from image Variable IF available
    getVarAttSafely(s3Img, imgVar);
    
    if (imgType == Nadir0500) {
        imgVar.getVar(s3Img->img);
    }
    else if (imgType == Obliq0500){
        getImgShifted(s3Img, imgProp, imgVar);
    }
    else if (imgType == NadirTpg || imgType == ObliqTpg){
        getImgScaled(s3Img, imgProp, imgVar);
    }
    
}*/

void S3NcdfData::readImgBinned(S3BasicImage<short>* s3Img, const std::string& ncdfName, const std::string varName, 
                               const NcdfImageType& imgType){
    
    cout << varName << endl;
    // open file
    NcFile ncF(ncdfName, NcFile::read);
    s3Img->name = varName;    
    
    NcVar imgVar = ncF.getVar(varName);
    if (imgVar.isNull()) {
        throw exceptions::NcNotVar("var is null", varName.c_str(), 26);
    }

    // get add_offset, scale_factor and _FillValue attriubutes 
    // from image Variable IF available
    getVarAttSafely(s3Img, imgVar);
    
    if (varName.find("radiance") != std::string::npos){
        if (varName[varName.size()-1] == 'n'){
            getBinRadImg(s3Img, s3MetaData.slstrPInfo.nadirImg0500m, imgVar);
        }
        else if (varName[varName.size()-1] == 'o'){
            getBinRadImg(s3Img, s3MetaData.slstrPInfo.obliqImg0500m, imgVar);
        }
    }
    /*if (imgType == Nadir0500) {
        if (s3Img->imgP.isBinned){
            getBinImg(s3Img, s3MetaData.slstrPInfo.nadirImg0500m, imgVar);
        }
        else {
            imgVar.getVar(s3Img->img);
        }
    }
    else if (imgType == Obliq0500){
        if (s3Img->imgP.isBinned){
            getBinImgShifted(s3Img, s3MetaData.slstrPInfo.obliqImg0500m, imgVar);
        }
        else {
            getImgShifted(s3Img, s3MetaData.slstrPInfo.obliqImg0500m, imgVar);
        }
    }
    else if (imgType == NadirTpg){
        if (s3Img->imgP.isBinned){
            getBinImgScaled(s3Img, s3MetaData.slstrPInfo.nadirTpgImg, imgVar);
        }
        else {
            getImgScaled(s3Img, s3MetaData.slstrPInfo.nadirTpgImg, imgVar);
        }
    }
    else if (imgType == ObliqTpg){
        if (s3Img->imgP.isBinned){
            getBinImgScaled(s3Img, s3MetaData.slstrPInfo.obliqTpgImg, imgVar);
        }
        else {
            getImgScaled(s3Img, s3MetaData.slstrPInfo.obliqTpgImg, imgVar);
        }
    }*/
    
}

void S3NcdfData::readImgBinned(S3BasicImage<ushort>* s3Img, const std::string& ncdfName, const std::string varName, 
                               const NcdfImageType& imgType){
    
    cout << varName << endl;
    // open file
    NcFile ncF(ncdfName, NcFile::read);
    s3Img->name = varName;    
    
    NcVar imgVar = ncF.getVar(varName);
    if (imgVar.isNull()) {
        throw exceptions::NcNotVar("var is null", varName.c_str(), 26);
    }

    // get add_offset, scale_factor and _FillValue attriubutes 
    // from image Variable IF available
    //getVarAttSafely(s3Img, imgVar);
    
    if ((varName.find("confid") != std::string::npos)
         || (varName.find("cloud") != std::string::npos)){
        if (varName[varName.size()-1] == 'n'){
            getFlagImg(s3Img, s3MetaData.slstrPInfo.nadirImg0500m, imgVar);
        }
        else if (varName[varName.size()-1] == 'o'){
            getFlagImg(s3Img, s3MetaData.slstrPInfo.obliqImg0500m, imgVar);
        }
    }
    
    /*if (imgType == Nadir0500) {
        if (s3Img->imgP.isBinned){
            getBinImg(s3Img, s3MetaData.slstrPInfo.nadirImg0500m, imgVar);
        }
        else {
            imgVar.getVar(s3Img->img);
        }
    }
    else if (imgType == Obliq0500){
        if (s3Img->imgP.isBinned){
            getBinImgShifted(s3Img, s3MetaData.slstrPInfo.obliqImg0500m, imgVar);
        }
        else {
            getImgShifted(s3Img, s3MetaData.slstrPInfo.obliqImg0500m, imgVar);
        }
    }
    else if (imgType == NadirTpg){
        if (s3Img->imgP.isBinned){
            getBinImgScaled(s3Img, s3MetaData.slstrPInfo.nadirTpgImg, imgVar);
        }
        else {
            getImgScaled(s3Img, s3MetaData.slstrPInfo.nadirTpgImg, imgVar);
        }
    }
    else if (imgType == ObliqTpg){
        if (s3Img->imgP.isBinned){
            getBinImgScaled(s3Img, s3MetaData.slstrPInfo.obliqTpgImg, imgVar);
        }
        else {
            getImgScaled(s3Img, s3MetaData.slstrPInfo.obliqTpgImg, imgVar);
        }
    }*/
    
}

void S3NcdfData::readImgBinned(S3BasicImage<double>* s3Img, const std::string& ncdfName, const std::string varName, 
                               const NcdfImageType& imgType){
    
    cout << varName << endl;
    // open file
    NcFile ncF(ncdfName, NcFile::read);
    s3Img->name = varName;    
    
    NcVar imgVar = ncF.getVar(varName);
    if (imgVar.isNull()) {
        throw exceptions::NcNotVar("var is null", varName.c_str(), 26);
    }

    // get add_offset, scale_factor and _FillValue attriubutes 
    // from image Variable IF available
    getVarAttSafely(s3Img, imgVar);
    
    if (ncdfName.find("geodetic") != std::string::npos){
        if (varName[varName.size()-1] == 'n'){
            getBinGeoLocImg(s3Img, s3MetaData.slstrPInfo.nadirImg0500m, imgVar);
        }
        else if (varName[varName.size()-1] == 'o'){
            getBinGeoLocImg(s3Img, s3MetaData.slstrPInfo.obliqImg0500m, imgVar);
        }
    }
    
    if (ncdfName.find("geometry") != std::string::npos){
        if (varName[varName.size()-1] == 'n'){
            getBinGeomImg(s3Img, s3MetaData.slstrPInfo.nadirTpgImg, imgVar);
        }
        else if (varName[varName.size()-1] == 'o'){
            getBinGeomImg(s3Img, s3MetaData.slstrPInfo.obliqTpgImg, imgVar);
        }
    }
    
    
    
    
    
    /*if (imgType == Nadir0500) {
        if (s3Img->imgP.isBinned){
            getBinImg(s3Img, s3MetaData.slstrPInfo.nadirImg0500m, imgVar);
        }
        else {
            imgVar.getVar(s3Img->img);
        }
    }
    else if (imgType == Obliq0500){
        if (s3Img->imgP.isBinned){
            getBinImgShifted(s3Img, s3MetaData.slstrPInfo.obliqImg0500m, imgVar);
        }
        else {
            getImgShifted(s3Img, s3MetaData.slstrPInfo.obliqImg0500m, imgVar);
        }
    }
    else if (imgType == NadirTpg){
        if (s3Img->imgP.isBinned){
            getBinImgScaled(s3Img, s3MetaData.slstrPInfo.nadirTpgImg, imgVar);
        }
        else {
            getImgScaled(s3Img, s3MetaData.slstrPInfo.nadirTpgImg, imgVar);
        }
    }
    else if (imgType == ObliqTpg){
        if (s3Img->imgP.isBinned){
            getBinImgScaled(s3Img, s3MetaData.slstrPInfo.obliqTpgImg, imgVar);
        }
        else {
            getImgScaled(s3Img, s3MetaData.slstrPInfo.obliqTpgImg, imgVar);
        }
    }*/
    
}

void S3NcdfData::getFlagImg(S3BasicImage<ushort>* s3Img, const ImageProperties& imgProp, const NcVar& imgVar){
    // if no binning or shifting is required simply return the img
    if (!s3Img->imgP.isBinned && !(imgProp == s3MetaData.slstrPInfo.obliqImg0500m)){
        imgVar.getVar(s3Img->img);
        return;
    }
    
    const ImageProperties& origNadirProp = s3MetaData.slstrPInfo.nadirImg0500m;
    const int dx = origNadirProp.xOff - imgProp.xOff;
    const int dy = origNadirProp.yOff - ( imgProp.yOff - 0 );

    S3BasicImage<ushort> tmp(imgProp);
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

void S3NcdfData::getBinRadImg(S3BasicImage<short>* s3Img, const ImageProperties& imgProp, const NcVar& imgVar){
    char isObliq = (imgProp == s3MetaData.slstrPInfo.obliqImg0500m) ? 1 : 0;
    const float binThrs = 0.5 * s3Img->imgP.binSize * s3Img->imgP.binSize;
    // if no binning or shifting is required simply return the img
    if (!s3Img->imgP.isBinned && !isObliq){
        imgVar.getVar(s3Img->img);
        return;
    }
    
    S3BasicImage<short> tmp(imgProp);
    imgVar.getVar(tmp.img);
    
    const ImageProperties& origNadirProp = s3MetaData.slstrPInfo.nadirImg0500m;
    const int dx = origNadirProp.xOff - imgProp.xOff;
    const int dy = origNadirProp.yOff - ( imgProp.yOff - 0 );

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
            for (int j2 = j1 * winSize; j2 < (j1 + 1) * winSize; j2++){
                for (int i2 = i1 * winSize; i2 < (i1 + 1) * winSize; i2++){
                    // i2, j2 are the actual idx in the nadir 500 image including the respective offset
                    idx2 = j2 * origNadirProp.width + i2;
                    i3 = i2 - dx;
                    j3 = j2 - dy;
                    if (i3 >= 0 && i3 < tmp.imgP.width && j3 >= 0 && j3 < tmp.imgP.height){
                        idx3 = j3 * tmp.imgP.width + i3;
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
                s3Img->img[j1 * s3Img->imgP.width + i1] = sumLand / nL;
                flags.img[j1 * flags.imgP.width + i1] |= (1 << isObliq);
            }
            else if (nO > binThrs){
                s3Img->img[j1 * s3Img->imgP.width + i1] = sumOcn / nO;
                flags.img[j1 * flags.imgP.width + i1] |= (4 << isObliq);
            }
            else {
                s3Img->img[j1 * s3Img->imgP.width + i1] = s3Img->noData;
                //flags.img[j1 * flags.imgP.width + i1] = 0;
            }
        }
    }
}

void S3NcdfData::getBinGeoLocImg(S3BasicImage<double>* s3Img, const ImageProperties& imgProp, const NcVar& imgVar){

    S3BasicImage<double> tmp(imgProp);
    imgVar.getVar(tmp.img);
    tmp.copyVarAttFrom(*s3Img);
    s3Img->valOffset = 0.;
    s3Img->valScale = 1.;
    s3Img->noData = -999.;

    const ImageProperties& origNadirProp = s3MetaData.slstrPInfo.nadirImg0500m;
    const int dx = origNadirProp.xOff - imgProp.xOff;
    const int dy = origNadirProp.yOff - ( imgProp.yOff - 0 );

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

/*void S3NcdfData::getBinImgShifted(S3BasicImage<short>* s3Img, const ImageProperties& imgProp, const NcVar& imgVar){
    S3BasicImage<short> tmp(imgProp);
    imgVar.getVar(tmp.img);
    
    const ImageProperties& origNadirProp = s3MetaData.slstrPInfo.nadirImg0500m;
    const int dx = origNadirProp.xOff - imgProp.xOff;
    const int dy = origNadirProp.yOff - ( imgProp.yOff - 0 );
    
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
                        if (tmp.img[idx2] > 0 && tmp.img[idx2] < 10000) {
                            sum += tmp.img[idx2];
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

void S3NcdfData::getBinImgShifted(S3BasicImage<double>* s3Img, const ImageProperties& imgProp, const NcVar& imgVar){
    S3BasicImage<double> tmp(imgProp);
    imgVar.getVar(tmp.img);
    
    const ImageProperties& origNadirProp = s3MetaData.slstrPInfo.nadirImg0500m;
    const int dx = origNadirProp.xOff - imgProp.xOff;
    const int dy = origNadirProp.yOff - ( imgProp.yOff - 0 );
    
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
                        if (tmp.img[idx2] > 0 && tmp.img[idx2] < 10000) {
                            sum += tmp.img[idx2];
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

void S3NcdfData::getBinImgScaled(S3BasicImage<short>* s3Img, const ImageProperties& imgProp, const NcVar& imgVar){
    S3BasicImage<short> tmp(imgProp);
    imgVar.getVar(tmp.img);
    
    const ImageProperties& origNadirProp = s3MetaData.slstrPInfo.nadirImg0500m;
    bool disCont = (s3Img->name.find("longitude") == 0);
    const int& winSize = s3Img->imgP.binSize;
    double val, sum, y2, x2;
    int n;
    // i/j 1 : iterate in binned output img
    // i/j 2 : iterate in full nadir original res img 
    for (int j1 = 0; j1 < s3Img->imgP.height; j1++){
        for (int i1 = 0; i1 < s3Img->imgP.width; i1++){
            n = 0;
            sum = 0;
            for (int j2 = j1 * winSize; j2 < (j1 + 1) * winSize; j2++){
                y2 = (j2 - origNadirProp.yOff ) * (double)(origNadirProp.yRes) / tmp.imgP.yRes + tmp.imgP.yOff;
                for (int i2 = i1 * winSize; i2 < (i1 + 1) * winSize; i2++){
                    // i2, j2 are the actual idx in the nadir 500 image including the respective offset
                    // x2, y2 are the intermediate interpolated indices at tpg resolution
                    x2 = (i2 - origNadirProp.xOff ) * (double)(origNadirProp.xRes) / tmp.imgP.xRes + tmp.imgP.xOff;
                    val = interpol_2d_img(tmp.img, x2, y2, tmp.imgP.width, tmp.imgP.height, disCont);
                    if (val > 0 && val < 10000) {
                        sum += val;
                        n++;
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
}*/

void S3NcdfData::getBinGeomImg(S3BasicImage<double>* s3Img, const ImageProperties& imgProp, const NcVar& imgVar){
    S3BasicImage<double> tmp(imgProp);
    imgVar.getVar(tmp.img);
    tmp.copyVarAttFrom(*s3Img);
    
    const ImageProperties& origNadirProp = s3MetaData.slstrPInfo.nadirImg0500m;
    bool disCont = false;//(s3Img->name.find("azimuth") != std::string::npos);
    const int& winSize = s3Img->imgP.binSize;
    double val, sum, y2, x2;
    int n;
    // i/j 1 : iterate in binned output img
    // i/j 2 : iterate in full nadir original res img 
    for (int j1 = 0; j1 < s3Img->imgP.height; j1++){
        for (int i1 = 0; i1 < s3Img->imgP.width; i1++){
            n = 0;
            sum = 0;
            for (int j2 = j1 * winSize; j2 < (j1 + 1) * winSize; j2++){
                y2 = (j2 - origNadirProp.yOff ) * (double)(origNadirProp.yRes) / tmp.imgP.yRes + tmp.imgP.yOff;
                for (int i2 = i1 * winSize; i2 < (i1 + 1) * winSize; i2++){
                    // i2, j2 are the actual idx in the nadir 500 image including the respective offset
                    // x2, y2 are the intermediate interpolated indices at tpg resolution
                    x2 = (i2 - origNadirProp.xOff ) * (double)(origNadirProp.xRes) / tmp.imgP.xRes + tmp.imgP.xOff;
                    val = interpol_2d_img(tmp.img, x2, y2, tmp.imgP.width, tmp.imgP.height, disCont);
                    //cerr << x2 << "/" << y2 << ": " << val << endl;
                    if (s3Img->isValidValue(val)) {
                        sum += val;
                        n++;
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

/*
void S3NcdfData::getImgShifted(S3BasicImage<short>* s3Img, const ImageProperties& imgProp, const NcVar& imgVar){
    cout << "shifting ";
    cout << imgProp.width << "x";
    cout << imgProp.height <<"+";
    cout << imgProp.xOff << "+";
    cout << imgProp.yOff << " to ";
    cout << s3Img->imgP.width << "x";
    cout << s3Img->imgP.height <<"+";
    cout << s3Img->imgP.xOff << "+";
    cout << s3Img->imgP.yOff << endl;
    
    S3BasicImage<short> tmp(imgProp);
    imgVar.getVar(tmp.img);
    int dx = s3Img->imgP.xOff - imgProp.xOff;
    int dy = s3Img->imgP.yOff - ( imgProp.yOff - 0 );
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

void S3NcdfData::getImgShifted(S3BasicImage<double>* s3Img, const ImageProperties& imgProp, const NcVar& imgVar){
    cout << "shifting ";
    cout << imgProp.width << "x";
    cout << imgProp.height <<"+";
    cout << imgProp.xOff << "+";
    cout << imgProp.yOff << " to ";
    cout << s3Img->imgP.width << "x";
    cout << s3Img->imgP.height <<"+";
    cout << s3Img->imgP.xOff << "+";
    cout << s3Img->imgP.yOff << endl;
    
    S3BasicImage<double> tmp(imgProp);
    imgVar.getVar(tmp.img);
    int dx = s3Img->imgP.xOff - imgProp.xOff;
    int dy = s3Img->imgP.yOff - ( imgProp.yOff - 0 );
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

void S3NcdfData::getImgScaled(S3BasicImage<short>* s3Img, const ImageProperties& imgProp, const NcVar& imgVar){
    cout << "scaling + shifting ";
    cout << imgProp.width << "x";
    cout << imgProp.height <<"+";
    cout << imgProp.xOff << "+";
    cout << imgProp.yOff << " to ";
    cout << s3Img->imgP.width << "x";
    cout << s3Img->imgP.height <<"+";
    cout << s3Img->imgP.xOff << "+";
    cout << s3Img->imgP.yOff << endl;
    
    S3BasicImage<short> tmp(imgProp);
    imgVar.getVar(tmp.img);
    bool disCont = (s3Img->name.find("longitude") == 0);
    for (int y1 = 0; y1 < s3Img->imgP.height; y1++){
        double y2 = (y1) * (double)(s3Img->imgP.yRes) / tmp.imgP.yRes;
        //double y2 = (y1 - s3Img->yOff ) * (double)(s3Img->yRes) / tmp.yRes + tmp.yOff;
        for (int x1 = 0; x1 < s3Img->imgP.width; x1++){
            double x2 = (x1 - s3Img->imgP.xOff ) * (double)(s3Img->imgP.xRes) / tmp.imgP.xRes + tmp.imgP.xOff;
            s3Img->img[y1 * s3Img->imgP.width + x1] = interpol_2d_img(tmp.img, x2, y2, tmp.imgP.width, tmp.imgP.height, disCont);
        }
    }
}

void S3NcdfData::getImgScaled(S3BasicImage<double>* s3Img, const ImageProperties& imgProp, const NcVar& imgVar){
    cout << "scaling + shifting ";
    cout << imgProp.width << "x";
    cout << imgProp.height <<"+";
    cout << imgProp.xOff << "+";
    cout << imgProp.yOff << " to ";
    cout << s3Img->imgP.width << "x";
    cout << s3Img->imgP.height <<"+";
    cout << s3Img->imgP.xOff << "+";
    cout << s3Img->imgP.yOff << endl;
    
    S3BasicImage<double> tmp(imgProp);
    imgVar.getVar(tmp.img);
    bool disCont = (s3Img->name.find("longitude") == 0);
    for (int y1 = 0; y1 < s3Img->imgP.height; y1++){
        double y2 = (y1) * (double)(s3Img->imgP.yRes) / tmp.imgP.yRes;
        //double y2 = (y1 - s3Img->yOff ) * (double)(s3Img->yRes) / tmp.yRes + tmp.yOff;
        for (int x1 = 0; x1 < s3Img->imgP.width; x1++){
            double x2 = (x1 - s3Img->imgP.xOff ) * (double)(s3Img->imgP.xRes) / tmp.imgP.xRes + tmp.imgP.xOff;
            s3Img->img[y1 * s3Img->imgP.width + x1] = interpol_2d_img(tmp.img, x2, y2, tmp.imgP.width, tmp.imgP.height, disCont);
        }
    }
}
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
    std::string ncdfName = s3md->s3FileName + "/" + s3md->CHANNEL_RAD_NAME[0][0] + ".nc";
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
void S3NcdfData::getVarAttSafely(S3BasicImage<short>* s3Img, NcVar& imgVar){
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
void S3NcdfData::getVarAttSafely(S3BasicImage<double>* s3Img, NcVar& imgVar){
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
            ncdfName = s3MetaData.s3FileName + "/" + s3MetaData.CHANNEL_QUAL_NAME[iView][iBand] + ".nc";
            NcFile ncF(ncdfName, NcFile::read);
            NcVar irradVar = ncF.getVar(s3MetaData.CHANNEL_IRRAD_NAME[iView][iBand]);
            irradVar.getVar(irval);
            irrad[iView][iBand] = 0;
            for (int i = 0; i < N_DETECTORS; i++){
                irrad[iView][iBand] += irval[i];
            }
            irrad[iView][iBand] /= N_DETECTORS;
            ncF.close();
        }
    }   
}

void S3NcdfData::createLandMask(){
    s3LandImg = S3BasicImage<signed char>(s3MetaData.slstrPInfo.nadirImg0500m);
    s3LandImg.name = "landmask";
    s3LandImg.initImgArray();
    const ushort srfcTypeTest = 0b0000000000011111;

    for (int iView = 0; iView < N_SLSTR_VIEWS; iView++){
        //read full resolution flag images
        std::string ncdfName = s3MetaData.s3FileName + "/" + s3MetaData.FLAGS_NAME[iView] + ".nc";
        S3BasicImage<ushort> s3ConfidImg = S3BasicImage<ushort>(s3MetaData.slstrPInfo.nadirImg0500m);
        s3ConfidImg.noData = 0;
        readImgBinned(&s3ConfidImg, ncdfName, s3MetaData.CONFID_NAME[iView], Nadir0500);
        for (int i = 0; i < s3LandImg.imgP.nPix; i++){
            ushort mask = (s3ConfidImg.img[i] & (srfcTypeTest));
            if (mask == 8){
                // is clear land
                s3LandImg.img[i] |= 1 << iView;
            }
            else if (mask == 2 || mask == (8+16)){
                // is clear water
                s3LandImg.img[i] |= 4 << iView;
            }
        }
    }
}

void S3NcdfData::createCloudMask(){
    s3CloudImg = S3BasicImage<signed char>(s3MetaData.slstrPInfo.nadirImg0500m);
    s3CloudImg.name = "cloudmask";
    s3CloudImg.initImgArray();
    const ushort cloudTest = 0b0011011111111110;

    for (int iView = 0; iView < N_SLSTR_VIEWS; iView++){
        //read full resolution flag images
        std::string ncdfName = s3MetaData.s3FileName + "/" + s3MetaData.FLAGS_NAME[iView] + ".nc";
        S3BasicImage<ushort> s3CldImg = S3BasicImage<ushort>(s3MetaData.slstrPInfo.nadirImg0500m);
        s3CldImg.noData = 0;
        readImgBinned(&s3CldImg, ncdfName, s3MetaData.BASIC_CLOUD_NAME[iView], Nadir0500);
        for (int i = 0; i < s3CloudImg.imgP.nPix; i++){
            ushort mask = (s3CldImg.img[i] & (cloudTest));
            if (mask){
                // is cloud
                s3CloudImg.img[i] |= 1 << iView;
            }
        }
    }
}

void S3NcdfData::createValidMask(){
    s3ValidImg = S3BasicImage<signed char>(s3MetaData.slstrPInfo.nadirImg0500m);
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

    for (int iView = 0; iView < N_SLSTR_VIEWS; iView++){
        signed char landMask  = 1 << iView;
        signed char ocnMask   = 4 << iView;
        signed char cloudMask = 1 << iView;
        
        for (int i = 0; i < s3ValidImg.imgP.nPix; i++){
            if ((s3LandImg.img[i] & landMask) && !(s3CloudImg.img[i] & cloudMask)) {
                s3ValidImg.img[i] |= landMask;
            }
            else if ((s3LandImg.img[i] & ocnMask) && !(s3CloudImg.img[i] & cloudMask)) {
                s3ValidImg.img[i] |= ocnMask;
            }
        }
    }

}

