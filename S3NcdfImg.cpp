/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   S3NcdfImg.cpp
 * Author: akheckel
 * 
 * Created on 23. November 2016, 20:43
 */

#include <iostream>
#include <netcdf>
#include "S3NcdfImg.hpp"
#include "InputParameter.hpp"

using std::cout;
using std::endl;
using namespace netCDF;


void S3NcdfImg::readNcdf(std::string ncdfName, std::string varName){
    
    NcFile ncF(ncdfName, NcFile::read);
    
    int xo, yo;
    ncF.getAtt("track_offset").getValues(&xo);
    ncF.getAtt("start_offset").getValues(&yo);
    
    name = varName;
    
    NcVar imgVar = ncF.getVar(varName);
    if (imgVar.isNull()) {
        throw exceptions::NcNotVar("var is null", varName.c_str(), 26);
    }
    
    imgVar.getAtt("add_offset").getValues(&valOffset);
    imgVar.getAtt("scale_factor").getValues(&valScale);
    imgVar.getAtt("_FillValue").getValues(&noData);
    
    size_t h = imgVar.getDim(0).getSize();
    size_t w = imgVar.getDim(1).getSize();
    if (w == width && h == height) {
        xOff = xo;
        yOff = yo;
        imgVar.getVar(img);
    }
    else {
        cout << "shifting ";
        cout << w << "x";
        cout << h <<"+";
        cout << xo << "+";
        cout << yo << " to ";
        cout << width << "x";
        cout << height <<"+";
        cout << xOff << "+";
        cout << yOff << endl;

        S3BasicImage tmp(w, h);
        imgVar.getVar(tmp.img);
        int dx = xOff - xo;
        int dy = yOff - yo;
        int i2, j2;
        for (int j1=0; j1<height; j1++){
            for (int i1=0; i1<width; i1++){
                i2 = i1 - dx;
                j2 = j1 - dy;
                if (j2 >= 0 && j2 < h && i2 >= 0 && i2 < w){
                    img[j1*width+i1] = tmp.img[j2*w+i2];
                }
                else {
                    img[j1*width+i1] = noData;
                }
            }
        }
    }
    
}

