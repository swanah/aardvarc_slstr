/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   AodRetrieval.hpp
 * Author: akheckel
 *
 * Created on 8. April 2017, 14:25
 */

#ifndef AODRETRIEVAL_HPP
#define AODRETRIEVAL_HPP

#include "defs.hpp"
#include "AtmosphericLut.hpp"
#include "OceanReflLut.hpp"
#include "AodRetrievalEmods.hpp"

class AodRetrieval {
public:
    SlstrPixel& pix;
    AtmosphericLut& atmLut;
    OceanReflLut& ocnLut;


    AodRetrieval(SlstrPixel& slstrPixel, AtmosphericLut& aLut, OceanReflLut& oLut);
    ~AodRetrieval();

    void retrieveAodSizeBrent(bool isOcean);
    void invertFixedAtm(bool isOcean, const float& aod, const float& fot, const float& wof, const float& doc);

private:
    AodRetrieval();
    AodRetrieval(const AodRetrieval& orig) : pix(orig.pix), atmLut(orig.atmLut), ocnLut(orig.ocnLut) {}; // disable copy construtor
    AodRetrieval& operator=(const AodRetrieval& rhs) { throw std::logic_error("AodRetrieval shoudlnt be copied"); } // disable copy assignment

    float getCurvature(EmodTau *emodTau, char isOcean);
/*
    float emod_size(float fot);
    float emod_tau_ocean_model(float tau);
    float emod_tau_uncer(float tau);
    float emod_uncer(float *);
*/
};



#endif /* AODRETRIEVAL_HPP */

