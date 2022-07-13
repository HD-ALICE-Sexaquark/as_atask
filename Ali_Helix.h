#ifndef __ALI_HELIX_H__
#define __ALI_HELIX_H__

#include "TObject.h"

class Ali_Helix : public TObject {
   private:
    Float_t aliHelix_params[6];

   public:
    // constructor
    Ali_Helix() : aliHelix_params() {}

    // destructor
    ~Ali_Helix() {}

    void setHelix(Float_t a, Float_t b, Float_t c, Float_t d, Float_t e, Float_t f) {
        aliHelix_params[0] = a;  // y0
        aliHelix_params[1] = b;  // z0
        aliHelix_params[2] = c;
        aliHelix_params[3] = d;
        aliHelix_params[4] = e;  // c curvature
        aliHelix_params[5] = f;  // x0
    }

    Float_t getHelix_param(Int_t i_param) const { return aliHelix_params[i_param]; }

    void Evaluate(Double_t t, Double_t r[3]) {
        // r[3] = radius vector
        // helix evaluation, taken from AliHelix
        // Calculate position of a point on a track and some derivatives at given phase
        Float_t phase = aliHelix_params[4] * t + aliHelix_params[2];
        Double_t sn = sinf(phase), cs = cosf(phase);
        //  Double_t sn=TMath::Sin(phase), cs=TMath::Cos(phase);

        r[0] = aliHelix_params[5] + sn / aliHelix_params[4];
        r[1] = aliHelix_params[0] - cs / aliHelix_params[4];
        r[2] = aliHelix_params[1] + aliHelix_params[3] * t;
    }

    ClassDef(Ali_Helix, 1);
};

#endif
