/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   S3BasicImage.cpp
 * Author: akheckel
 * 
 * Created on 23. November 2016, 19:21
 */

#include <iostream>
#include <algorithm>
#include "S3BasicImage.hpp"
#include "InputParameter.hpp"

S3BasicImage::S3BasicImage(const int& newWidth, const int& newHeight, const int& xOffset, const int& yOffset){
    width = newWidth;
    height = newHeight;
    xOff = xOffset;
    yOff = yOffset;
    img = new short[width*height];
}

S3BasicImage::S3BasicImage(const S3BasicImage& orig) {
    width = orig.width;
    height = orig.height;
    xOff = orig.xOff;
    yOff = orig.yOff;
    valScale = orig.valOffset;
    valOffset = orig.valOffset;
    noData = orig.noData;
    int n = width*height;
    img = new short[n];
    for (int i=0; i<n; i++) img[i] = orig.img[i];
    name = orig.name;
}

S3BasicImage& S3BasicImage::operator=(const S3BasicImage& rhs) {
    // Check for self-assignment!
    if (this == &rhs) // Same object?
        return *this; // Yes, so skip assignment, and just return *this.
    S3BasicImage lhs(rhs);
    std::swap(width, lhs.width);
    std::swap(height, lhs.height);
    std::swap(xOff, lhs.xOff);
    std::swap(yOff, lhs.yOff);
    std::swap(valScale, lhs.valScale);
    std::swap(valOffset, lhs.valOffset);
    std::swap(noData, lhs.noData);
    std::swap(img, lhs.img);
    std::swap(name, lhs.name);
    return *this;
}


S3BasicImage::~S3BasicImage() {
    delete [] img;
    img = NULL;
    //std::cout << "S3BasicImage: img deleted!\n";
}

