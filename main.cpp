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
#include <sys/unistd.h>
#include <stdexcept>
#include "InputParameter.hpp"
#include "S3MetaData.hpp"
#include "S3NcdfData.hpp"
#include "Images.hpp"
#include "Interpolation.hpp"
#include "miscUtils.hpp"
#include "AeroClimatology.hpp"
#include "OceanReflLut.hpp"
#include "AtmosphericLut.hpp"
#include "AodRetrieval.hpp"

using std::cout;
using std::cerr;
using std::endl;
using namespace netCDF;

void addWriteVar(NcFile* ncOut, const std::vector<NcDim>& dimVec, const S3BasicImage<signed char>& img);
void addWriteVar(NcFile* ncOut, const std::vector<NcDim>& dimVec, const S3BasicImage<short>& img);
void addWriteVar(NcFile* ncOut, const std::vector<NcDim>& dimVec, const S3BasicImage<unsigned int>& img);
void addWriteVar(NcFile* ncOut, const std::vector<NcDim>& dimVec, const S3BasicImage<float>& img);
void addWriteVar(NcFile* ncOut, const std::vector<NcDim>& dimVec, const S3BasicImage<double>& img);
void writeGlobalAttributes(netCDF::NcFile& ncFile, S3NcdfData& s3Data, const int& nPix);
void computeOutputImageProp(ImageProperties* outputImgProp, const ImageProperties& inputImgProp, const int& winSize);
bool isPixelValidOcean(SlstrPixel *pixel, const InputParameter& pars);
float calcPrevFF(S3BasicImage<float> *fineTotalFrac, int x, int y);

/*
 *
 */
int main(int argc, char** argv) {
    try {
        InputParameter pars(argc, argv);
        S3NcdfData s3Data(pars);
        ImageProperties outImgProp;
        computeOutputImageProp(&outImgProp, pars.s3MD.slstrPInfo.nadirImg0500m, pars.winSize);

        cout << pars.s3MD.productName << endl;

        /* reading nadir and oblique radiance images */
        cout << "reading product" << endl;
        double t1, t2;
        t1 = timer();

        s3Data.readNcdf(outImgProp);

        t2 = timer();
        printf("Runtime %f seconds\n", (float)(t2 - t1));

        s3Data.convRad2Refl();
        s3Data.verifyInput();
        s3Data.initResultImgs(outImgProp);
        int& imgWidth  = outImgProp.width;  //s3Data.s3MetaData.slstrPInfo.nadirImg0500m.width;
        int& imgHeight = outImgProp.height; //s3Data.s3MetaData.slstrPInfo.nadirImg0500m.height;

        float nPix=1;
/******/

        cout << "reading LUTs" << endl;
        t1 = timer();

        // read aerosol climatology
        AeroClimatology aerClim(pars.climFileName, pars.s3MD.getMonth());

        // read ocean LUT
        OceanReflLut ocnLut(pars.ocnLutFileName);

        // read atmospheric LUT
        AtmosphericLut lut(pars.atmLutFileName);

        t2 = timer();
        printf("Time to read Luts %f seconds\n", (float)(t2 - t1));

        SlstrPixel slstrPixel;
        int idx;
        nPix=0;
        for (slstrPixel.x = 0; slstrPixel.x < imgWidth; slstrPixel.x++){
//        for (slstrPixel.x = 127; slstrPixel.x < 128; slstrPixel.x++){
            fprintf(stdout, "processing %6.2f%%\r", (float)(slstrPixel.x)/imgWidth*100.0); fflush(stdout);
            for (slstrPixel.y = 0; slstrPixel.y < imgHeight; slstrPixel.y++){
//            for (slstrPixel.y = 67; slstrPixel.y < 68; slstrPixel.y++){
//                fprintf(stdout, "processing %5d / %5d\n", slstrPixel.x, slstrPixel.y); fflush(stdout);
                idx = slstrPixel.y * imgWidth + slstrPixel.x;

                slstrPixel.qflag = 0;
                slstrPixel.prevFineFrac = 0;
                slstrPixel.ndvi = 1;
                slstrPixel.ndvi_veg_weight = 1;
                slstrPixel.dust_weight = 0;
                slstrPixel.aod = -1;
                slstrPixel.ax = lut.aodD.min;
                slstrPixel.cx = lut.aodD.max;
                for (int i=0; i<N_MP; i++) slstrPixel.model_p[i] = -1;

                if (s3Data.isValidPixel(idx)){
                    s3Data.getGeoPos(idx, &slstrPixel.geo_pos);
                    s3Data.getViewGeom(idx, &slstrPixel.geom);
                    s3Data.getPres(idx, &slstrPixel.pAlt);
                    slstrPixel.ocn_wind_speed = DEFAULT_WDSP;
                    slstrPixel.ocn_wind_dir   = DEFAULT_WDIR;
                    slstrPixel.ocn_pigment    = DEFAULT_PIG;

                    slstrPixel.view_clear[0] = ((s3Data.flags.img[idx] & (4+1)) > 0);
                    slstrPixel.view_clear[1] = ((s3Data.flags.img[idx] & (8+2)) > 0);
                    setBit(&slstrPixel.qflag, CLR_LAND, ((s3Data.flags.img[idx] & (3)) == 3));
                    setBit(&slstrPixel.qflag, CLR_OCEAN_N, ((s3Data.flags.img[idx] & (4)) > 0));
                    setBit(&slstrPixel.qflag, CLR_OCEAN_F, ((s3Data.flags.img[idx] & (8)) > 0));

                    //s3Data.getToaReflec(idx, slstrPixel.tarr);
                    s3Data.getToaReflecF(idx, slstrPixel.tarr);
                    
                    if (slstrPixel.x == -1 && slstrPixel.y == -1){
                        printf("Pixel (%d/%d) is valid and at %s\n", slstrPixel.x, slstrPixel.y, slstrPixel.geo_pos.toCstr());
                        printf("Pixel SZA: %f SAA: %f VZA: %f VAA: %f RAZ: %f\n", slstrPixel.geom.nad_sol_zen, slstrPixel.geom.nad_sol_azim, slstrPixel.geom.nad_sat_zen, slstrPixel.geom.nad_sat_azim, slstrPixel.geom.razn);
                        printf("Pixel SZA: %f SAA: %f VZA: %f VAA: %f RAZ: %f\n", slstrPixel.geom.for_sol_zen, slstrPixel.geom.for_sol_azim, slstrPixel.geom.for_sat_zen, slstrPixel.geom.for_sat_azim, slstrPixel.geom.razf);
                        printf("Pixel TOA-N %f %f %f %f %f\n", slstrPixel.tarr[0][0], slstrPixel.tarr[1][0], slstrPixel.tarr[2][0], slstrPixel.tarr[3][0], slstrPixel.tarr[4][0]);
                        printf("Pixel TOA-O %f %f %f %f %f\n", slstrPixel.tarr[0][1], slstrPixel.tarr[1][1], slstrPixel.tarr[2][1], slstrPixel.tarr[3][1], slstrPixel.tarr[4][1]);
                    }
                    //get mixing paramter
                    //prep lut interpolation
                    //aerClim.getMixPercentages(slstrPixel.geo_pos, slstrPixel.lutpars.mixing, slstrPixel.lutpars.mix_frac, &slstrPixel.lutpars.climAod);
                    aerClim.getMixPercentagesInt(slstrPixel.geo_pos, slstrPixel.lutpars.mixing, slstrPixel.lutpars.mix_frac, &slstrPixel.lutpars.climAod);
                    lut.getTetrahedronPoints(&slstrPixel.lutpars, false);
                    try {
                        slstrPixel.lutpars.razni = lut.getInterPar(slstrPixel.geom.razn, lut.razD);
                        slstrPixel.lutpars.razfi = lut.getInterPar(slstrPixel.geom.razf, lut.razD);

                        slstrPixel.lutpars.szani = lut.getInterPar(slstrPixel.geom.nad_sol_zen, lut.szaD);
                        slstrPixel.lutpars.szafi = lut.getInterPar(slstrPixel.geom.for_sol_zen, lut.szaD);

                        slstrPixel.lutpars.vzani = lut.getInterPar(slstrPixel.geom.nad_sat_zen, lut.vzaD);
                        slstrPixel.lutpars.vzafi = lut.getInterPar(slstrPixel.geom.for_sat_zen, lut.vzaD);

                        if (slstrPixel.pAlt > lut.presD.max) slstrPixel.pAlt = lut.presD.max;
                        if (slstrPixel.pAlt < lut.presD.min) slstrPixel.pAlt = lut.presD.min;
                        slstrPixel.lutpars.pAlti = lut.getInterPar(slstrPixel.pAlt, lut.presD);

                        slstrPixel.lutpars.ocn_razni = ocnLut.getInterPar(slstrPixel.geom.razn, ocnLut.razD);
                        slstrPixel.lutpars.ocn_razfi = ocnLut.getInterPar(slstrPixel.geom.razf, ocnLut.razD);
                        slstrPixel.lutpars.ocn_mi  = ocnLut.getInterPar(slstrPixel.lutpars.mix_frac[0], ocnLut.modelD);
                        slstrPixel.lutpars.ocn_wsi = ocnLut.getInterPar(slstrPixel.ocn_wind_speed, ocnLut.wdspD);
                        slstrPixel.lutpars.ocn_wdi = ocnLut.getInterPar(slstrPixel.ocn_wind_dir, ocnLut.wdirD);
                        slstrPixel.lutpars.ocn_pi  = ocnLut.getInterPar(slstrPixel.ocn_pigment, ocnLut.pigcD);

                        ocnLut.get_rho_ocean_wind(&slstrPixel, 0.05, 9);
                        slstrPixel.rho_glint[0] = slstrPixel.rho_surf[3][0];
                        slstrPixel.rho_glint[1] = slstrPixel.rho_surf[3][1];

                        slstrPixel.prevFineFrac = calcPrevFF(&s3Data.s3AerFracImgs[1], slstrPixel.x, slstrPixel.y);

                        AodRetrieval retrieval(slstrPixel, lut, ocnLut);
                        if ((s3Data.flags.img[idx] & 15) == 3){
                            //clear land in both views
                            retrieval.retrieveAodSizeBrent(false);
                            //retrieval.invertFixedAtm(false, 0.15/*slstrPixel.lutpars.climAod*/, 0.9, 0.9, 0.86);
                            s3Data.setRetrievalResults(idx, slstrPixel);
                            nPix++;
                        }
                        else if (isPixelValidOcean(&slstrPixel, pars)) {
                            //clear ocean in either view
                            retrieval.retrieveAodSizeBrent(true);
                            //retrieval.invertFixedAtm(false, 0.15, 0.9, 0.9, 0.86);
                            s3Data.setRetrievalResults(idx, slstrPixel);
                            nPix++;
                        }
                    }
                    catch (std::range_error){
                        fprintf(stderr, "skipping pixel, is not inside LUT\n");
                        //TODO: set possible flag for outside lut and continue
                    }

                }
                else {
                    //printf("Pixel (%d/%d) is invalid\n", slstrPixel.x, slstrPixel.y);
                }

                s3Data.flags.img[idx] |= slstrPixel.qflag;

            }
        }

        // writing radiance data to ncdf file
        t1 = timer();
        printf("Time to retrieve AOD %f seconds\n", (float)(t1 - t2));
/*******/
        if (nPix>0){
            cout << "writing: " << pars.aodOutName << endl;
            NcFile ncOut(pars.aodOutName, NcFile::replace);
            NcDim xDim = ncOut.addDim("columns", imgWidth);
            NcDim yDim = ncOut.addDim("rows", imgHeight);
            std::vector<NcDim> dimVec;
            dimVec.push_back(yDim);
            dimVec.push_back(xDim);
            
            writeGlobalAttributes(ncOut, s3Data, nPix);

            for (int iBand = 0; iBand < N_SLSTR_BANDS; iBand++){
                addWriteVar(&ncOut, dimVec, s3Data.s3AodImgs[iBand]);
            }
            for (int i = 0; i < N_AER_FRAC; i++){
                addWriteVar(&ncOut, dimVec, s3Data.s3AerFracImgs[i]);
            }
            for (int iBand = 0; iBand < N_SLSTR_BANDS; iBand++){
                addWriteVar(&ncOut, dimVec, s3Data.s3UncImgs[iBand]);
            }
            addWriteVar(&ncOut, dimVec, s3Data.s3FmAodImg);
            addWriteVar(&ncOut, dimVec, s3Data.s3DustAodImg);
            addWriteVar(&ncOut, dimVec, s3Data.s3AbsAodImg);
            addWriteVar(&ncOut, dimVec, s3Data.s3AngstromImg);
            addWriteVar(&ncOut, dimVec, s3Data.s3SsaImg);
            addWriteVar(&ncOut, dimVec, s3Data.s3FminImg);
            for (int iView = 0; iView < N_SLSTR_VIEWS; iView++){
                for (int iBand = 0; iBand < N_SLSTR_BANDS; iBand++){
                    //addWriteVar(&ncOut, dimVec, s3Data.s3RadImgs[iView][iBand]);
                    addWriteVar(&ncOut, dimVec, s3Data.s3RadImgsF[iView][iBand]);
                    addWriteVar(&ncOut, dimVec, s3Data.s3SdrImgs[iView][iBand]);
                    
                    //addWriteVar(&ncOut, dimVec, s3Data.s3RPathImgs[iView][iBand]);
                    //addWriteVar(&ncOut, dimVec, s3Data.s3TDownImgs[iView][iBand]);
                    //addWriteVar(&ncOut, dimVec, s3Data.s3TUpImgs[iView][iBand]);
                    //addWriteVar(&ncOut, dimVec, s3Data.s3TGasImgs[iView][iBand]);
                    //addWriteVar(&ncOut, dimVec, s3Data.s3SpherAImgs[iView][iBand]);
                    //addWriteVar(&ncOut, dimVec, s3Data.s3DifFracImgs[iView][iBand]);
                }

                addWriteVar(&ncOut, dimVec, s3Data.s3SzaImgs[iView]);
                addWriteVar(&ncOut, dimVec, s3Data.s3SaaImgs[iView]);
                addWriteVar(&ncOut, dimVec, s3Data.s3VzaImgs[iView]);
                addWriteVar(&ncOut, dimVec, s3Data.s3VaaImgs[iView]);
                addWriteVar(&ncOut, dimVec, s3Data.s3RazImgs[iView]);
                addWriteVar(&ncOut, dimVec, s3Data.s3ScatAngImgs[iView]);
            }
            s3Data.s3LatImgs[0].name = "latitude";
            addWriteVar(&ncOut, dimVec, s3Data.s3LatImgs[0]);
            s3Data.s3LonImgs[0].name = "longitude";
            addWriteVar(&ncOut, dimVec, s3Data.s3LonImgs[0]);
            for (int iCnr = 0; iCnr < 4; iCnr++) {
                addWriteVar(&ncOut, dimVec, s3Data.s3LatCnrImgs[iCnr]);
                addWriteVar(&ncOut, dimVec, s3Data.s3LonCnrImgs[iCnr]);
            }
            addWriteVar(&ncOut, dimVec, s3Data.s3PresImg);
            addWriteVar(&ncOut, dimVec, s3Data.flags);
            addWriteVar(&ncOut, dimVec, s3Data.s3NPixImg);
            addWriteVar(&ncOut, dimVec, s3Data.s3TimeImg);
            
            for (int i=0; i<N_MP; i++){
                addWriteVar(&ncOut, dimVec, s3Data.s3ModelParImgs[i]);
            }

            addWriteVar(&ncOut, dimVec, s3Data.s3AodClimImg);

            //ncOut.close();
            t2 = timer();
            printf("Runtime %f seconds\n", (float)(t2 - t1));
        }
        else {
            printf("%s\nNo valid pixels to process\n", pars.s3MD.productName.c_str());
        }

        /* writing hi res flag data to ncdf file *
        imgWidth  = s3Data.pars.s3MD.slstrPInfo.nadirImg0500m.width;
        imgHeight = s3Data.pars.s3MD.slstrPInfo.nadirImg0500m.height;
        replaceStringInPlace(pars.aodOutName, ".nc", "_flags.nc");
        cout << "writing: " << pars.aodOutName << endl;
        NcFile ncOut1(pars.aodOutName, NcFile::replace);
        NcDim xDim = ncOut1.addDim("columns", imgWidth);
        NcDim yDim = ncOut1.addDim("rows", imgHeight);
        std::vector<NcDim> dimVec;
        dimVec.push_back(yDim);
        dimVec.push_back(xDim);

        addWriteVar(&ncOut1, dimVec, s3Data.s3LandImg);
        addWriteVar(&ncOut1, dimVec, s3Data.s3CloudImg);
        addWriteVar(&ncOut1, dimVec, s3Data.s3NanImg);
        addWriteVar(&ncOut1, dimVec, s3Data.s3ValidImg);

        //ncOut1.close();
        /**/

        //S3CciL2Writer cciWriter(s3Data, nPix);
        //cciWriter.write();
        
    }
    catch (exceptions::NcException& e){
        cerr << e.what() << endl << "unrecoverable netcdf error" << endl;
        return 1; //e.errorCode();
    }
    /*catch (boost::filesystem::filesystem_error& e){
        cerr << e.what() << endl << "boost FS error" << endl;
        return 1;        
    }*/
    catch (std::runtime_error& e){
        cerr << e.what() << endl << "runtime error" << endl;
        return 1;
    }
    catch (...){
        cerr << "main: exception thrown, but error unknown (maybe in XML parsing...)" << endl;
        return 1;
    }
    cout << "finished!" << endl;
    return 0;
}

/**
 * netcdf Helper routine to add and write an image to the netcdf file
 * @param ncOut - netCdf file id
 * @param dimVec - vector of dimensions for that variable
 * @param img - the image variable
 */
void addWriteVar(NcFile* ncOut, const std::vector<NcDim>& dimVec, const S3BasicImage<signed char>& img){
    NcVar var = ncOut->addVar(img.name, ncByte, dimVec);

    var.setCompression(true, true, 5);

    if (img.hasOffset){
        var.putAtt("add_offset", ncDouble, img.valOffset);
    }
    if (img.hasScale){
        var.putAtt("scale_factor", ncDouble, img.valScale);
    }
    if (img.hasNoData){
        var.putAtt("_FillValue", ncByte, img.noData);
    }
    if (img.hasMinMax){
        var.putAtt("valid_min", ncByte, img.validMin);
        var.putAtt("valid_max", ncByte, img.validMax);
    }
    if (img.isSpecBand){
        var.putAtt("wavelength", ncFloat, img.wvl);
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
    if (img.hasMinMax){
        var.putAtt("valid_min", ncShort, img.validMin);
        var.putAtt("valid_max", ncShort, img.validMax);
    }
    if (img.isSpecBand){
        var.putAtt("wavelength", ncFloat, img.wvl);
    }
    var.putAtt("coordinates", std::string("latitude_an longitude_an"));

    var.putVar(img.img);
}

void addWriteVar(NcFile* ncOut, const std::vector<NcDim>& dimVec, const S3BasicImage<unsigned int>& img){
    NcVar var = ncOut->addVar(img.name, ncUint, dimVec);

    var.setCompression(true, true, 5);

    if (img.hasOffset){
        var.putAtt("add_offset", ncDouble, img.valOffset);
    }
    if (img.hasScale){
        var.putAtt("scale_factor", ncDouble, img.valScale);
    }
    if (img.hasNoData){
        var.putAtt("_FillValue", ncUint, img.noData);
    }
    if (img.hasMinMax){
        var.putAtt("valid_min", ncUint, img.validMin);
        var.putAtt("valid_max", ncUint, img.validMax);
    }
    if (img.isSpecBand){
        var.putAtt("wavelength", ncFloat, img.wvl);
    }
    var.putAtt("coordinates", std::string("latitude_an longitude_an"));

    var.putVar(img.img);
}

void addWriteVar(NcFile* ncOut, const std::vector<NcDim>& dimVec, const S3BasicImage<float>& img){
    NcVar var = ncOut->addVar(img.name, ncFloat, dimVec);

    var.setCompression(true, true, 5);

    if (img.hasOffset){
        var.putAtt("add_offset", ncDouble, img.valOffset);
    }
    if (img.hasScale){
        var.putAtt("scale_factor", ncDouble, img.valScale);
    }
    if (img.hasNoData){
        var.putAtt("_FillValue", ncFloat, img.noData);
    }
    if (img.hasMinMax){
        var.putAtt("valid_min", ncFloat, img.validMin);
        var.putAtt("valid_max", ncFloat, img.validMax);
    }
    if (img.isSpecBand){
        var.putAtt("wavelength", ncFloat, img.wvl);
    }
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
    if (img.hasMinMax){
        var.putAtt("valid_min", ncDouble, img.validMin);
        var.putAtt("valid_max", ncDouble, img.validMax);
    }
    if (img.isSpecBand){
        var.putAtt("wavelength", ncFloat, img.wvl);
    }
    var.putVar(img.img);
}

/**
 * write global attributes to netcdf aod output
 * @param ncFile
 */
void writeGlobalAttributes(netCDF::NcFile& ncFile, S3NcdfData& s3Data, const int& nPix) {
    
    const S3MetaData& s3MD = s3Data.pars.s3MD;
    //const boost::uuids::uuid uuid = boost::uuids::random_generator()();
    //const std::string uuidStr = boost::uuids::to_string(uuid);
    const size_t stmpLen = 1000;
    char stmp[stmpLen] = "";
    
    ncFile.putAtt("Conventions", "CF-1.6");
    //ncFile.putAtt("tracking_id", uuidStr);
    ncFile.putAtt("naming_authority", "uk.ac.su.aatsraerosol");
    ncFile.putAtt("title", "AARDVARC CCI aerosol product level 2");
    ncFile.putAtt("product_version", "1.0");
    ncFile.putAtt("summary", "This dataset contains the level-2 aerosol properties products from SLSTR satellite observations. Data are processed by Swansea algorithm");
    //ncFile.putAtt("id", ncFile.getName(true));
    ncFile.putAtt("sensor", s3MD.instrument);
    ncFile.putAtt("platform", s3MD.platformFamily + s3MD.platformNumber);
    ncFile.putAtt("resolution", "10km x 10km");
    ncFile.putAtt("projection", "quasi-cartesian instrument grid");
    ncFile.putAtt("cdm_data_type", "Swath");
    ncFile.putAtt("cell", "");
    ncFile.putAtt("abs_orbit", NC_INT, s3MD.orbitNumberAbs);
    ncFile.putAtt("rel_orbit", NC_INT, s3MD.orbitNumberRel);
    ncFile.putAtt("inputFileList", s3MD.productName);
    ncFile.putAtt("time_coverage_start", convS3DateToCci(s3MD.startTime));
    ncFile.putAtt("time_coverage_end", convS3DateToCci(s3MD.stopTime));
    snprintf(stmp, stmpLen, "%.6f", s3Data.getLatMin());
    ncFile.putAtt("geospatial_lat_min", NC_DOUBLE, s3Data.getLatMin());
    snprintf(stmp, stmpLen, "%.6f", s3Data.getLatMax());
    ncFile.putAtt("geospatial_lat_max", NC_DOUBLE, s3Data.getLatMax());
    snprintf(stmp, stmpLen, "%.6f", s3Data.getLonMin());
    ncFile.putAtt("geospatial_lon_min", NC_DOUBLE, s3Data.getLonMin());
    snprintf(stmp, stmpLen, "%.6f", s3Data.getLonMax());
    ncFile.putAtt("geospatial_lon_max", NC_DOUBLE, s3Data.getLonMax());
    ncFile.putAtt("date_created", getCurrentTimeStr());
    ncFile.putAtt("project", "Climate Change Initiative - European Space Agency");
    ncFile.putAtt("references", "http://www.esa-aerosol-cci.org");
    ncFile.putAtt("creator_name", "College of Science, Swansea University");
    ncFile.putAtt("creator_url", "http://www.swan.ac.uk/staff/academic/environmentsociety/geography/northpeter/");
    ncFile.putAtt("creator_email", "p.r.j.north@swansea.ac.uk \n a.heckel@swansea.ac.uk");
    ncFile.putAtt("source", s3MD.lv1Info);
    ncFile.putAtt("keywords", "satellite,observation,atmosphere");
    ncFile.putAtt("keywords_vocabulary", "NASA Global Change Master Directory (GCMD) Science Keywords");
    ncFile.putAtt("standard_name_vocabulary", "NetCDF Climate and Forecast (CF) Metadata Convention version 18");
    ncFile.putAtt("license", "ESA CCI Data Policy: free and open access");
    //ncFile.putAtt("", "");
}

/**
 * compute properties of the resulting AOD product images
 * based on the size of the input nadir image and the binning parameters
 * @param outputImgProp
 * @param inputImgProp
 * @param winSize
 */
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
        //if (inputImgProp.width > outputImgProp->width * winSize) outputImgProp->width++;

        outputImgProp->height = inputImgProp.height / winSize;
        //if (inputImgProp.height > outputImgProp->height * winSize) outputImgProp->height++;

        outputImgProp->nPix = outputImgProp->width * outputImgProp->height;

        outputImgProp->xRes = inputImgProp.xRes * winSize;
        outputImgProp->yRes = inputImgProp.yRes * winSize;

        outputImgProp->xOff = (int)((double)(inputImgProp.xOff) / winSize + 0.5);  //* inputImgProp.xRes / outputImgProp->xRes;
        outputImgProp->yOff = (int)((double)(inputImgProp.yOff) / winSize + 0.5);  //* inputImgProp.yRes / outputImgProp->yRes;
    }
}

/**
 * test if the current pixel is at least in one view not glint contaminated
 * @param pixel
 * @param pars
 * @return valid
 */
bool isPixelValidOcean(SlstrPixel *pixel, const InputParameter& pars) {
    bool ok = false;
    int j;
    
    if (pixel->view_clear[0]){
        ok = ok || (pixel->rho_glint[0] < GLINT_THRS);
        setBit(&pixel->qflag, GLINT_NADIR, (pixel->rho_glint[0] >= GLINT_THRS));
    }
    if (pixel->view_clear[1]){
        ok = ok || (pixel->rho_glint[1] < GLINT_THRS);
        setBit(&pixel->qflag, GLINT_FWARD, (pixel->rho_glint[1] >= GLINT_THRS));
    }
    return ok;
}

/**
 * compute an average of the previously retrieved fine mode fraction 
 * surrounding the current pixel at (x/y)
 * @param image of fineTotalFrac 
 * @param x
 * @param y
 * @return previous fine mode fraction
 */
float calcPrevFF(S3BasicImage<float> *fineTotalFrac, int x, int y) {
    float fmf, prevFF = 0, w, weight = 0;
    int ix, iy, xs, ys;
    xs = x - 5;
    ys = y - 5;

    for (ix = x; ((ix >= 0) && (ix >= xs)); ix -= 1) {
        for (iy = ((x == ix)?(y - 1):(y + 2)); ((iy >= 0) && (iy >= ys)); iy -= 1) {
            fmf = fineTotalFrac->img[ix + iy * fineTotalFrac->imgP.width];
            if (fmf > 0) {
                w = 1 / ( pow((ix-x),2) + pow((iy-y),2) );
                weight += w;
                prevFF += w * fmf;
            }
        }
    }
    if (weight > 1e-4) {
        return (prevFF / weight);
    }
    return 0;
}

