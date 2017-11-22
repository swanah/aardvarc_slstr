/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Images.h
 * Author: akheckel
 *
 * Created on 24. November 2016, 18:56
 */

#ifndef IMAGES_HPP
#define IMAGES_HPP

#include <cstddef>
#include <iostream>
#include <algorithm>
#include <limits>

#include "defs.hpp"


template<class T>
class S3BasicImage {
public:
    ImageProperties imgP;
    bool hasScale, hasOffset;
    double valScale, valOffset;
    bool hasMinMax;
    T validMin, validMax;
    T *img;
    bool hasNoData;
    T noData;
    std::string name;
    bool isSpecBand;
    float wvl;

    S3BasicImage();
    S3BasicImage(const int& newWidth, const int& newHeight, 
                    const int& xOffset = 0, const int& yOffset = 0, 
                    const int& xRes = 1, const int& yRes = 1);
    S3BasicImage(const ImageProperties& imgProp);
    S3BasicImage(const S3BasicImage& orig);
    S3BasicImage<T>& operator=(const S3BasicImage& rhs);
    ~S3BasicImage();
    
    void setImgProperties(ImageProperties& imgProp);
    void getImgProperties(ImageProperties* imgProp);
    void setFillVal(const T& fillVal);
    void setScaleOffset(const double& scale, const double& offset);
    void setWvl(const float& spectralWvl);
    void setValidLimits(const T& min, const T& max);
    bool isValidValue(const T& value);
    void copyVarAttFrom(const S3BasicImage& srcImg);
    void initImgArray(const T initVal = 0);
};

template<class T>
S3BasicImage<T>::S3BasicImage(){ 
    img = NULL;
}

template<class T>
S3BasicImage<T>::S3BasicImage(const int& newWidth, const int& newHeight, 
                              const int& xOffset, const int& yOffset, 
                              const int& xRes, const int& yRes){
    imgP.width = newWidth;
    imgP.height = newHeight;
    imgP.nPix = newWidth * newHeight;
    imgP.xOff = xOffset;
    imgP.yOff = yOffset;
    imgP.xRes = xRes;
    imgP.yRes = yRes;
    hasScale = false;
    valScale = 1.;
    hasOffset = false;
    valOffset = 0.;
    hasNoData = false;
    noData = std::numeric_limits<T>::quiet_NaN();
    isSpecBand = false;
    wvl = 0.;
    hasMinMax = false;
    validMin = std::numeric_limits<T>::min();
    validMax = std::numeric_limits<T>::max();
    img = new T[imgP.width * imgP.height];
}

template<class T>
S3BasicImage<T>::S3BasicImage(const ImageProperties& newImgProp){
    imgP = newImgProp;
    hasScale = false;
    valScale = 1.;
    hasOffset = false;
    valOffset = 0.;
    hasNoData = false;
    noData = std::numeric_limits<T>::quiet_NaN();
    isSpecBand = false;
    wvl = 0.;
    hasMinMax = false;
    validMin = std::numeric_limits<T>::min();
    validMax = std::numeric_limits<T>::max();
    img = new T[imgP.width * imgP.height];
}

template<class T>
S3BasicImage<T>::S3BasicImage(const S3BasicImage& orig) {
    imgP = orig.imgP;
    hasScale = orig.hasScale;
    valScale = orig.valScale;
    hasOffset = orig.hasOffset;
    valOffset = orig.valOffset;
    hasNoData = orig.hasNoData;
    noData = orig.noData;
    isSpecBand = orig.isSpecBand;
    wvl = orig.wvl;
    hasMinMax = orig.hasMinMax;
    validMin = orig.validMin;
    validMax = orig.validMax;
    int n = imgP.width*imgP.height;
    img = new T[n];
    for (int i=0; i<n; i++) img[i] = orig.img[i];
    name = orig.name;
}

template<class T>
S3BasicImage<T>& S3BasicImage<T>::operator=(const S3BasicImage& rhs) {
    // Check for self-assignment!
    if (this == &rhs) // Same object?
        return *this; // Yes, so skip assignment, and just return *this.
    S3BasicImage lhs(rhs);
    std::swap(imgP, lhs.imgP);
    std::swap(hasScale, lhs.hasScale);
    std::swap(valScale, lhs.valScale);
    std::swap(hasOffset, lhs.hasOffset);
    std::swap(valOffset, lhs.valOffset);
    std::swap(hasNoData, lhs.hasNoData);
    std::swap(noData, lhs.noData);
    std::swap(isSpecBand, lhs.isSpecBand);
    std::swap(wvl, lhs.wvl);
    std::swap(hasMinMax, lhs.hasMinMax);
    std::swap(validMin, lhs.validMin);
    std::swap(validMax, lhs.validMax);
    std::swap(img, lhs.img);
    std::swap(name, lhs.name);
    return *this;
}



template<class T>
S3BasicImage<T>::~S3BasicImage(){ 
        delete [] img;
        img = NULL; 
}


template<class T>
inline void S3BasicImage<T>::setImgProperties(ImageProperties& imgProp){
    imgP = imgProp;
}

template<class T>
inline void S3BasicImage<T>::getImgProperties(ImageProperties* imgProp){
    imgProp->width  = imgP.width;
    imgProp->height = imgP.height;
    imgProp->nPix = imgP.nPix;
    imgProp->xOff   = imgP.xOff;
    imgProp->yOff   = imgP.yOff;
    imgProp->xRes   = imgP.xRes;
    imgProp->yRes   = imgP.yRes;
}

template<class T>
inline void S3BasicImage<T>::setFillVal(const T& fillVal){
    hasNoData = true;
    noData = fillVal;
}

template<class T>
void S3BasicImage<T>::setScaleOffset(const double& scale, const double& offset){
    hasScale = true;
    valScale = scale;
    hasOffset = true;
    valOffset = offset;
}

template<class T>
inline void S3BasicImage<T>::setValidLimits(const T& min, const T& max){
    hasMinMax = true;
    validMin = min;
    validMax = max;
}

template<class T>
inline void S3BasicImage<T>::setWvl(const float& spectralWvl){
    isSpecBand = true;
    wvl = spectralWvl;
}

template<class T>
inline bool S3BasicImage<T>::isValidValue(const T& value){
    bool isValid = true;//(!std::isnan<T>(value)) && (!std::isinf<T>(value));
    if (isValid && hasMinMax) {
        isValid = (value >= validMin) && (value <= validMax);
    }
    if (isValid && hasNoData) {
        isValid = (value != noData);
    }
    return isValid;
}

template<class T>
inline void S3BasicImage<T>::copyVarAttFrom(const S3BasicImage& srcImg){
    hasNoData = srcImg.hasNoData;
    noData = srcImg.noData;
    hasOffset = srcImg.hasOffset;
    valOffset = srcImg.valOffset;
    hasScale = srcImg.hasScale;
    valScale = srcImg.valScale;
    hasMinMax = srcImg.hasMinMax;
    validMin = srcImg.validMin;
    validMax = srcImg.validMax;
    isSpecBand = srcImg.isSpecBand;
    wvl = srcImg.wvl;
}

template<class T>
inline void S3BasicImage<T>::initImgArray(const T initVal){
    for (int i=0; i<imgP.nPix; i++) img[i] = initVal;
}
#endif /* IMAGES_HPP */

