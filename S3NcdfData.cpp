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
#include <netcdf>
#include "Images.hpp"
#include "Interpolation.hpp"
#include "miscUtils.hpp"
#include "S3NcdfData.hpp"

using std::cout;
using std::cerr;
using std::endl;
using namespace netCDF;


S3NcdfData::S3NcdfData() {
}

S3NcdfData::~S3NcdfData() {
}

void S3NcdfData::readNcdf(S3MetaData& s3md) {
    
    ImageProperties imgProp;
    readImageProp(&s3md, &imgProp);
    
    // reading nadir and oblique radiance images
    std::string ncdfName;
    for (int iView = 0; iView < N_SLSTR_VIEWS; iView++){
        for (int iBand = 0; iBand < N_SLSTR_BANDS; iBand++){
            s3RadImgs[iView][iBand] = S3BasicImage<short>(imgProp);
            ncdfName = s3md.s3FileName + "/" + s3md.CHANNEL_RAD_NAME[iView][iBand] + ".nc";
            readImg(ncdfName, s3md.CHANNEL_RAD_NAME[iView][iBand], &s3RadImgs[iView][iBand]);
        }

        ncdfName = s3md.s3FileName + "/" + s3md.GEODETIC_NAME[iView] + ".nc";
        s3LatImgs[iView] = S3BasicImage<double>(imgProp);
        readImg(ncdfName, s3md.LAT_NAME[iView], &s3LatImgs[iView]);

        s3LonImgs[iView] = S3BasicImage<double>(imgProp);
        readImg(ncdfName, s3md.LON_NAME[iView], &s3LonImgs[iView]);
        
        ncdfName = s3md.s3FileName + "/" + s3md.GEOMETRY_NAME[iView] + ".nc";
        s3SzaImgs[iView] = S3BasicImage<double>(imgProp);
        readImg(ncdfName, s3md.SZA_NAME[iView], &s3SzaImgs[iView]);
        
        s3SaaImgs[iView] = S3BasicImage<double>(imgProp);
        readImg(ncdfName, s3md.SAA_NAME[iView], &s3SaaImgs[iView]);
        
        s3VzaImgs[iView] = S3BasicImage<double>(imgProp);
        readImg(ncdfName, s3md.VZA_NAME[iView], &s3VzaImgs[iView]);
        
        s3VaaImgs[iView] = S3BasicImage<double>(imgProp);
        readImg(ncdfName, s3md.VAA_NAME[iView], &s3VaaImgs[iView]);
    }

    ncdfName = s3md.s3FileName + "/" + s3md.GEODETIC_TPG_NAME + ".nc";
    s3TpgLatImg = S3BasicImage<double>(imgProp);
    readImg(ncdfName, s3md.LAT_TPG_NAME, &s3TpgLatImg);

    s3TpgLonImg = S3BasicImage<double>(imgProp);
    readImg(ncdfName, s3md.LON_TPG_NAME, &s3TpgLonImg);
        
}

void S3NcdfData::readImg(const std::string& ncdfName, const std::string varName, S3BasicImage<short>* s3Img){
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
    
    if (imgProp.width == s3Img->width && imgProp.height == s3Img->height) {
        imgVar.getVar(s3Img->img);
    }
    else if (imgProp.xRes == s3Img->xRes && imgProp.yRes == s3Img->yRes){
        cout << "shifting ";
        cout << imgProp.width << "x";
        cout << imgProp.height <<"+";
        cout << imgProp.xOff << "+";
        cout << imgProp.yOff << " to ";
        cout << s3Img->width << "x";
        cout << s3Img->height <<"+";
        cout << s3Img->xOff << "+";
        cout << s3Img->yOff << endl;

        S3BasicImage<short> tmp(imgProp);
        imgVar.getVar(tmp.img);
        int dx = s3Img->xOff - imgProp.xOff;
        int dy = s3Img->yOff - ( imgProp.yOff - 0 );
        int i2, j2;
        for (int j1 = 0; j1 < s3Img->height; j1++){
            for (int i1 = 0; i1 < s3Img->width; i1++){
                i2 = i1 - dx;
                j2 = j1 - dy;
                if (j2 >= 0 && j2 < tmp.height && i2 >= 0 && i2 < tmp.width){
                    s3Img->img[ j1 * s3Img->width + i1 ] = tmp.img[ j2 * tmp.width + i2 ];
                }
                else {
                    s3Img->img[ j1 * s3Img->width + i1 ] = s3Img->noData;
                }
            }
        }
    }
    else {
        cout << "scaling + shifting ";
        cout << imgProp.width << "x";
        cout << imgProp.height <<"+";
        cout << imgProp.xOff << "+";
        cout << imgProp.yOff << " to ";
        cout << s3Img->width << "x";
        cout << s3Img->height <<"+";
        cout << s3Img->xOff << "+";
        cout << s3Img->yOff << endl;

        S3BasicImage<short> tmp(imgProp);
        imgVar.getVar(tmp.img);
        bool disCont = (s3Img->name.find("longitude") == 0);
        for (int y1 = 0; y1 < s3Img->height; y1++){
            double y2 = (y1) * (double)(s3Img->yRes) / tmp.yRes;
            //double y2 = (y1 - s3Img->yOff ) * (double)(s3Img->yRes) / tmp.yRes + tmp.yOff;
            for (int x1 = 0; x1 < s3Img->width; x1++){
                double x2 = (x1 - s3Img->xOff ) * (double)(s3Img->xRes) / tmp.xRes + tmp.xOff;
                s3Img->img[y1 * s3Img->width + x1] = interpol_2d_img(tmp.img, x2, y2, tmp.width, tmp.height, disCont);
            }
        }
    }
    
}

void S3NcdfData::readImg(const std::string& ncdfName, const std::string varName, S3BasicImage<double>* s3Img){
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
    
    if (imgProp.width == s3Img->width && imgProp.height == s3Img->height) {
        imgVar.getVar(s3Img->img);
    }
    else if (imgProp.xRes == s3Img->xRes && imgProp.yRes == s3Img->yRes){
        cout << "shifting ";
        cout << imgProp.width << "x";
        cout << imgProp.height <<"+";
        cout << imgProp.xOff << "+";
        cout << imgProp.yOff << " to ";
        cout << s3Img->width << "x";
        cout << s3Img->height <<"+";
        cout << s3Img->xOff << "+";
        cout << s3Img->yOff << endl;

        S3BasicImage<double> tmp(imgProp);
        imgVar.getVar(tmp.img);
        int dx = s3Img->xOff - imgProp.xOff;
        int dy = s3Img->yOff - ( imgProp.yOff - 0 );
        int i2, j2;
        for (int j1 = 0; j1 < s3Img->height; j1++){
            for (int i1 = 0; i1 < s3Img->width; i1++){
                i2 = i1 - dx;
                j2 = j1 - dy;
                if (j2 >= 0 && j2 < tmp.height && i2 >= 0 && i2 < tmp.width){
                    s3Img->img[ j1 * s3Img->width + i1 ] = tmp.img[ j2 * tmp.width + i2 ];
                }
                else {
                    s3Img->img[ j1 * s3Img->width + i1 ] = s3Img->noData;
                }
            }
        }
    }
    else {
        cout << "scaling + shifting ";
        cout << imgProp.width << "x";
        cout << imgProp.height <<"+";
        cout << imgProp.xOff << "+";
        cout << imgProp.yOff << " to ";
        cout << s3Img->width << "x";
        cout << s3Img->height <<"+";
        cout << s3Img->xOff << "+";
        cout << s3Img->yOff << endl;

        S3BasicImage<double> tmp(imgProp);
        imgVar.getVar(tmp.img);
        bool disCont = (s3Img->name.find("longitude") == 0);
        for (int y1 = 0; y1 < s3Img->height; y1++){
            double y2 = (y1) * (double)(s3Img->yRes) / tmp.yRes;
            //double y2 = (y1 - s3Img->yOff ) * (double)(s3Img->yRes) / tmp.yRes + tmp.yOff;
            for (int x1 = 0; x1 < s3Img->width; x1++){
                double x2 = (x1 - s3Img->xOff ) * (double)(s3Img->xRes) / tmp.xRes + tmp.xOff;
                s3Img->img[y1 * s3Img->width + x1] = interpol_2d_img(tmp.img, x2, y2, tmp.width, tmp.height, disCont);
            }
        }
    }
    
}

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

