/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   defs.hpp
 * Author: akheckel
 *
 * Created on 25. November 2016, 10:43
 */

#ifndef DEFS_HPP
#define DEFS_HPP

#define N_SLSTR_BANDS 5
#define N_SLSTR_VIEWS 2

#define N_DETECTORS 4

struct ImageProperties {
    int width, height, nPix;
    int xOff, yOff;
    int xRes, yRes;
    
    bool isBinned;
    int  binSize;
    
    ImageProperties() {
        isBinned = false;
        binSize = 1;
        width = -1; height = -1; nPix = 0;
        xOff = 0; yOff = 0;
        xRes = -1; yRes = -1;
    }
    
    
    bool operator==(const ImageProperties& rhs) {
    return (  (width == rhs.width)
            &&(height == rhs.height)
            &&(nPix == rhs.nPix)
            &&(xOff == rhs.xOff)
            &&(yOff == rhs.yOff)
            &&(xRes == rhs.xRes)
            &&(yRes == rhs.yRes)
            &&(isBinned == rhs.isBinned)
            &&(binSize == binSize) );
    }

    bool operator==(const ImageProperties& rhs) const {
    return (  (width == rhs.width)
            &&(height == rhs.height)
            &&(nPix == rhs.nPix)
            &&(xOff == rhs.xOff)
            &&(yOff == rhs.yOff)
            &&(xRes == rhs.xRes)
            &&(yRes == rhs.yRes)
            &&(isBinned == rhs.isBinned)
            &&(binSize == binSize) );
    }
};

#endif /* DEFS_HPP */

