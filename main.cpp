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
#include "InputParameter.hpp"
#include "S3MetaData.hpp"
#include "S3NcdfData.hpp"
#include "Images.hpp"
#include "Interpolation.hpp"
#include "miscUtils.hpp"
#include "AeroClimatology.hpp"

using std::cout;
using std::cerr;
using std::endl;
using namespace netCDF;

void addWriteVar(NcFile* ncOut, const std::vector<NcDim>& dimVec, const S3BasicImage<signed char>& img);
void addWriteVar(NcFile* ncOut, const std::vector<NcDim>& dimVec, const S3BasicImage<short>& img);
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
        //s3Data.readNcdf(s3md);
        t1 = timer();
        //s3Data.readNcdf(outImgProp);
        t2 = timer();
        printf("Runtime %f seconds\n", t2 - t1);
        //s3Data.convRad2Refl();
        //s3Data.verifyInput();
        
        // read aerosol climatology
        AeroClimatology aerClim(pars.climFileName, s3Data.s3MetaData.getMonth());
        
        // read ocean LUT
        
                
        /* writing radiance data to ncdf file */
        t1 = timer();
        int& imgWidth  = outImgProp.width;  //s3Data.s3MetaData.slstrPInfo.nadirImg0500m.width;
        int& imgHeight = outImgProp.height; //s3Data.s3MetaData.slstrPInfo.nadirImg0500m.height;
        cout << "writing: " << s3Data.aodOutName << endl;
        NcFile ncOut(s3Data.aodOutName, NcFile::replace);
        NcDim xDim = ncOut.addDim("columns", imgWidth);
        NcDim yDim = ncOut.addDim("rows", imgHeight);
        std::vector<NcDim> dimVec;
        dimVec.push_back(yDim);
        dimVec.push_back(xDim);

        for (int iView = 0; iView < 2; iView++){
            for (int iBand = 0; iBand < 5; iBand++){
                addWriteVar(&ncOut, dimVec, s3Data.s3RadImgs[iView][iBand]);
            }

            addWriteVar(&ncOut, dimVec, s3Data.s3LatImgs[iView]);
            addWriteVar(&ncOut, dimVec, s3Data.s3LonImgs[iView]);
            addWriteVar(&ncOut, dimVec, s3Data.s3SzaImgs[iView]);
            addWriteVar(&ncOut, dimVec, s3Data.s3SaaImgs[iView]);
            addWriteVar(&ncOut, dimVec, s3Data.s3VzaImgs[iView]);
            addWriteVar(&ncOut, dimVec, s3Data.s3VaaImgs[iView]);
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
    var.putAtt("coordinates", std::string("latitude_an longitude_an"));
    
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

