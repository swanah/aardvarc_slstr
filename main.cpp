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
#include "InputParameter.h"
#include "S3MetaData.h"
#include "S3BasicImage.h"
#include "S3NcdfImg.h"

using std::cout;
using std::cerr;
using std::endl;
using namespace netCDF;


/*
 * 
 */
int main(int argc, char** argv) {

    InputParameter pars;
    S3MetaData& s3md = pars.s3MetaData;
    cout << s3md.productName << endl;
    try {

        S3NcdfImg s1ImgN(s3md.slstrPInfo.nadirImg0500m.width, s3md.slstrPInfo.nadirImg0500m.height);
        std::string ncdfName = s3md.s3FileName + "/" + s3md.CHANNEL_RAD_NAME[0][0] + ".nc";
        s1ImgN.readNcdf(ncdfName, s3md.CHANNEL_RAD_NAME[0][0]);
        cout << s1ImgN.width << "x";
        cout << s1ImgN.height <<"+";
        cout << s1ImgN.xOff << "+";
        cout << s1ImgN.yOff << endl;
        
        S3NcdfImg s1ImgO(s3md.slstrPInfo.nadirImg0500m.width, s3md.slstrPInfo.nadirImg0500m.height, s1ImgN.xOff, s1ImgN.yOff);
        ncdfName = s3md.s3FileName + "/" + s3md.CHANNEL_RAD_NAME[1][0] + ".nc";
        s1ImgO.readNcdf(ncdfName, s3md.CHANNEL_RAD_NAME[1][0]);
        cout << s1ImgO.width << "x";
        cout << s1ImgO.height <<"+";
        cout << s1ImgO.xOff << "+";
        cout << s1ImgO.yOff << endl;
        
        
    }
    catch (exceptions::NcException& e){
        cerr << "unrecoverable error, exiting...\n";
        cerr << e.what();
        return e.errorCode();
    }
    return 0;
}

