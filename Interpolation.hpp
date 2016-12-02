/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Interpolation.hpp
 * Author: akheckel
 *
 * Created on 30. November 2016, 11:59
 */

#ifndef INTERPOLATION_HPP
#define INTERPOLATION_HPP


template<class T>
T interpol_2d_img(T *img, double d2, double d1, const int& d2size, const int& d1size, bool hasDiscont){
    T t, u;
    int x1i[2], x2i[2], i;
    T val;
    T p[4];
    bool foundDiscont = false;

    if (d1 < 0.0) d1 = 0.0;
    if (d2 < 0.0) d2 = 0.0;

    t = d1 - (int) d1;
    u = d2 - (int) d2;

    x1i[0] = (int) d1;
    x2i[0] = (int) d2;
    x1i[1] = 1 + (int) d1;
    x2i[1] = 1 + (int) d2;

    if (x1i[1] == d1size) x1i[1] = x1i[0];
    if (x2i[1] == d2size) x2i[1] = x2i[0];

    p[0] = *(img + x1i[0] * d2size + x2i[0]);
    p[1] = *(img + x1i[1] * d2size + x2i[0]);
    p[2] = *(img + x1i[0] * d2size + x2i[1]);
    p[3] = *(img + x1i[1] * d2size + x2i[1]);
    
    if (hasDiscont){
        for (i=1; i<4; i++){
            if ((p[0]-p[i]) > 180) {p[i] += 360; /*fprintf(stderr, "added 360 to p[%1d]\n", i);*/ foundDiscont = true;}
            else if ((p[0]-p[i]) < -180) {p[i] -= 360; /*fprintf(stderr, "subtr 360 from p[%1d]\n", i);*/ foundDiscont = true;}
        }
    }

    val = (1 - t) * (1 - u) * p[0]
            + t * (1 - u) * p[1]
            + (1 - t) * u * p[2]
            + t * u * p[3];

    if (hasDiscont){
        //if ( foundDiscont ) fprintf(stderr, "corners: %f   %f   %f  %f  int.value: %f\n", p[0], p[1], p[2], p[3], val);
        if (val > 180) { /*fprintf(stderr, "adjust val %f --> %f \n", val, val-360);*/ val -= 360; }
        else if(val < -180) { /*fprintf(stderr, "adjust val %f --> %f \n", val, val+360);*/ val += 360; }
    }


    return val;

}



#endif /* INTERPOLATION_HPP */

