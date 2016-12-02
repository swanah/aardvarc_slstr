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
#include "InputParameter.hpp"
#include "S3MetaData.hpp"
#include "S3NcdfData.hpp"
#include "Images.hpp"
#include "Interpolation.hpp"

using std::cout;
using std::cerr;
using std::endl;
using namespace netCDF;

void addWriteVar(NcFile* ncOut, const std::vector<NcDim>& dimVec, const S3BasicImage<short>& img);
void addWriteVar(NcFile* ncOut, const std::vector<NcDim>& dimVec, const S3BasicImage<double>& img);


/*
 * 
 */
int main(int argc, char** argv) {

    InputParameter pars;
    S3MetaData& s3md = pars.s3MetaData;
    cout << s3md.productName << endl;
    try {
        
        S3NcdfData s3Data;
        /* reading nadir and oblique radiance images */
        s3Data.readNcdf(s3md);
        s3Data.convRad2Refl();
                
        /* writing radiance data to ncdf file */
        int& imgWidth  = s3md.slstrPInfo.nadirImg0500m.width;
        int& imgHeight = s3md.slstrPInfo.nadirImg0500m.height;
        cout << "writing: " << pars.aodOutName << endl;
        NcFile ncOut(pars.aodOutName, NcFile::replace);
        NcDim xDim = ncOut.addDim("columns", imgWidth);
        NcDim yDim = ncOut.addDim("rows", imgHeight);
        std::vector<NcDim> dimVec;
        dimVec.push_back(yDim);
        dimVec.push_back(xDim);

        NcVar var;
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

        ncOut.close();
    }
    catch (exceptions::NcException& e){
        cerr << "unrecoverable error, exiting...\n";
        cerr << e.what();
        return e.errorCode();
    }
    cout << "finished!" << endl;
    return 0;
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


