#ifndef __ALI_TRD_ST_TOF_HIT_H__
#define __ALI_TRD_ST_TOF_HIT_H__

#include "TObject.h"
#include "TVector3.h"

class Ali_TRD_ST_TOF_hit : public TObject {
   private:
    TVector3 TV3_offset;
    Float_t time;

   public:
    // constructor
    Ali_TRD_ST_TOF_hit() : TV3_offset(), time(0) {}

    // destructor
    ~Ali_TRD_ST_TOF_hit() {}

    // setters
    void set_TV3_offset(TVector3 TV3_offset_in) { TV3_offset = TV3_offset_in; }
    void set_time(Float_t time_in) { time = time_in; }

    // getters
    TVector3 get_TV3_offset() const { return TV3_offset; }
    Float_t get_time() const { return time; }

    ClassDef(Ali_TRD_ST_TOF_hit, 1);
};

#endif
