#ifndef __ALI_AS_OFFLINE_TRACKLET_H__
#define __ALI_AS_OFFLINE_TRACKLET_H__

#include "TObject.h"
#include "TVector3.h"

class Ali_AS_offline_Tracklet : public TObject {
   private:
    // Tracklet properties
    Short_t detector;
    TVector3 TV3_offset;
    TVector3 TV3_dir;
    Float_t chi2;
    Float_t refY;
    Float_t refZ;
    Float_t locY;
    Float_t locZ;
    Float_t refdYdx;
    Float_t refdZdx;
    Float_t locdYdx;
    Float_t locdZdx;

   public:
    // constructor
    Ali_AS_offline_Tracklet()
        : detector(0),
          TV3_offset(),
          TV3_dir(),
          chi2(0),
          refY(0),
          refZ(0),
          locY(0),
          locZ(0),
          refdYdx(0),
          refdZdx(0),
          locdYdx(0),
          locdZdx(0) {}

    // destructor
    ~Ali_AS_offline_Tracklet() {}

    // setters
    void set_detector(Short_t s) { detector = s; }
    void set_TV3_offset(TVector3 tv3) { TV3_offset = tv3; }
    void set_TV3_dir(TVector3 tv3) { TV3_dir = tv3; }
    void set_chi2(Float_t chi2in) { chi2 = chi2in; }
    void set_refYZ(Float_t refYin, Float_t refZin) {
        refY = refYin;
        refZ = refZin;
    }
    void set_refdYdZdx(Float_t refdYdxin, Float_t refdZdxin) {
        refdYdx = refdYdxin;
        refdZdx = refdZdxin;
    }
    void set_locYZ(Float_t locYin, Float_t locZin) {
        locY = locYin;
        locZ = locZin;
    }
    void set_locdYdZdx(Float_t locdYdxin, Float_t locdZdxin) {
        locdYdx = locdYdxin;
        locdZdx = locdZdxin;
    }

    // getters
    Short_t get_detector() const { return detector; }
    TVector3 get_TV3_offset() const { return TV3_offset; }
    TVector3 get_TV3_dir() const { return TV3_dir; }
    Float_t get_chi2() const { return chi2; }
    Float_t get_refZ() const { return refZ; }
    Float_t get_refY() const { return refY; }
    Float_t get_refdZdx() const { return refdZdx; }
    Float_t get_refdYdx() const { return refdYdx; }
    Float_t get_locZ() const { return locZ; }
    Float_t get_locY() const { return locY; }
    Float_t get_locdZdx() const { return locdZdx; }
    Float_t get_locdYdx() const { return locdYdx; }

    ClassDef(Ali_AS_offline_Tracklet, 1);  // A simple track of a particle
};

#endif
