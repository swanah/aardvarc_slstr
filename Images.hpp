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

typedef struct {
    int width, height;
    int xOff, yOff;
    int xRes, yRes;
} ImageProperties;


template<class T>
class S3BasicImage {
public:
    int width, height;
    int xOff, yOff;
    int xRes, yRes;
    bool hasScale, hasOffset;
    double valScale, valOffset;
    T *img;
    bool hasNoData;
    T noData;
    std::string name;

    S3BasicImage();
    S3BasicImage(const int& newWidth, const int& newHeight, const int& xOffset = 0, const int& yOffset = 0);
    S3BasicImage(const ImageProperties& imgProp);
    S3BasicImage(const S3BasicImage& orig);
    S3BasicImage<T>& operator=(const S3BasicImage& rhs);
    ~S3BasicImage();
    
    void setImgProperties(ImageProperties& imgProp);
    void getImgProperties(ImageProperties* imgProp);
    
};

template<class T>
S3BasicImage<T>::S3BasicImage(){ 
    img = NULL;
}

template<class T>
S3BasicImage<T>::S3BasicImage(const int& newWidth, const int& newHeight, const int& xOffset, const int& yOffset){
    width = newWidth;
    height = newHeight;
    xOff = xOffset;
    yOff = yOffset;
    xRes = 500;
    yRes = 500;
    hasScale = false;
    hasOffset = false;
    hasNoData = false;
    img = new T[width*height];
}

template<class T>
S3BasicImage<T>::S3BasicImage(const ImageProperties& imgProp){
    width = imgProp.width;
    height = imgProp.height;
    xOff = imgProp.xOff;
    yOff = imgProp.yOff;
    xRes = imgProp.xRes;
    yRes = imgProp.yRes;
    hasScale = false;
    hasOffset = false;
    hasNoData = false;
    img = new T[width*height];
}

template<class T>
S3BasicImage<T>::S3BasicImage(const S3BasicImage& orig) {
    width = orig.width;
    height = orig.height;
    xOff = orig.xOff;
    yOff = orig.yOff;
    xRes = orig.xRes;
    yRes = orig.yRes;
    hasScale = orig.hasScale;
    valScale = orig.valOffset;
    hasOffset = orig.hasOffset;
    valOffset = orig.valOffset;
    hasNoData = orig.hasNoData;
    noData = orig.noData;
    int n = width*height;
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
    std::swap(width, lhs.width);
    std::swap(height, lhs.height);
    std::swap(xOff, lhs.xOff);
    std::swap(yOff, lhs.yOff);
    std::swap(xRes, lhs.xRes);
    std::swap(yRes, lhs.yRes);
    std::swap(hasScale, lhs.hasScale);
    std::swap(valScale, lhs.valScale);
    std::swap(hasOffset, lhs.hasOffset);
    std::swap(valOffset, lhs.valOffset);
    std::swap(hasNoData, lhs.hasNoData);
    std::swap(noData, lhs.noData);
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
    width  = imgProp.width;
    height = imgProp.height;
    xOff   = imgProp.xOff;
    yOff   = imgProp.yOff;
    xRes   = imgProp.xRes;
    yRes   = imgProp.yRes;
}

template<class T>
inline void S3BasicImage<T>::getImgProperties(ImageProperties* imgProp){
    imgProp->width  = width;
    imgProp->height = height;
    imgProp->xOff   = xOff;
    imgProp->yOff   = yOff;
    imgProp->xRes   = xRes;
    imgProp->yRes   = yRes;
}




#endif /* IMAGES_HPP */

