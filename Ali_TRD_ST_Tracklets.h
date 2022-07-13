#ifndef __ALI_TRD_ST_TRACKLETS_H__
#define __ALI_TRD_ST_TRACKLETS_H__

#include "TObject.h"
#include "TVector3.h"

class Ali_TRD_ST_Tracklets : public TObject {
   private:
    TVector3 TV3_offset;
    TVector3 TV3_dir;
    Int_t TRD_det;
    Double_t ADC_val[30];
    UShort_t TPC_match;
    Int_t n_tracklets_around;
    Double_t min_dist_to_next_trkl;
    Int_t index;

   public:
    // constructor
    Ali_TRD_ST_Tracklets()
        : TV3_offset(), TV3_dir(), TRD_det(-4), ADC_val(), TPC_match(0), n_tracklets_around(0), min_dist_to_next_trkl(-5.0), index(0) {}

    // destructor
    ~Ali_TRD_ST_Tracklets() {}

    // setters
    void set_TV3_offset(TVector3 TV3_offset_in) { TV3_offset = TV3_offset_in; }
    void set_TV3_dir(TVector3 TV3_dir_in) { TV3_dir = TV3_dir_in; }
    void set_TRD_det(Int_t TRD_det_in) { TRD_det = TRD_det_in; }
    void set_ADC_val(Int_t time_bin, Double_t ADC_value) { ADC_val[time_bin] = ADC_value; }
    void set_TPC_match(UShort_t TPC_match_in) { TPC_match = TPC_match_in; }
    void set_n_tracklets_around(Int_t n_tracklets_around_in) { n_tracklets_around = n_tracklets_around_in; }
    void set_min_dist_to_next_trkl(Double_t min_dist_to_next_trkl_in) { min_dist_to_next_trkl = min_dist_to_next_trkl_in; }
    void set_TRD_index(Int_t TRD_IND_in) { index = TRD_IND_in; }

    // getters
    TVector3 get_TV3_offset() const { return TV3_offset; }
    TVector3 get_TV3_dir() const { return TV3_dir; }
    Int_t get_TRD_det() const { return TRD_det; }
    Double_t get_ADC_val(Int_t time_bin) const { return ADC_val[time_bin]; }
    UShort_t get_TPC_match() const { return TPC_match; }
    Int_t get_n_tracklets_around() const { return n_tracklets_around; }
    Double_t get_min_dist_to_next_trkl() const { return min_dist_to_next_trkl; }
    Int_t get_TRD_index() const { return index; }

    ClassDef(Ali_TRD_ST_Tracklets, 1);
};

#endif
