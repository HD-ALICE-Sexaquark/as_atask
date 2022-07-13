#ifndef __ALI_AS_TRD_DIGIT_H__
#define __ALI_AS_TRD_DIGIT_H__

#include "TObject.h"

class Ali_AS_TRD_digit : public TObject {
   private:
    // Digit data
    UShort_t hit_ids[2];  // contains the full information of the TRD hit position (sector, stack, layer, row, column), definition see below
    Short_t ADC_time_values[30];  // raw ADC values for all 24 time bins of a single pad
    Short_t arr_pos[3][30];
    Short_t dca_to_track;  // distance of closest approach of digit TRD hit to TPC track
    Short_t dca_x;
    Short_t dca_y;
    Short_t dca_z;
    Short_t ImpactAngle;  // impact angle of TPC track to this TRD digit

   public:
    // constructor
    Ali_AS_TRD_digit() : hit_ids(), ADC_time_values(), arr_pos(), dca_to_track(-1), dca_x(-1), dca_y(-1), dca_z(-1), ImpactAngle(-1) {}

    // destructor
    ~Ali_AS_TRD_digit() {}

    // setters
    void sethit_ids(UShort_t x, UShort_t y) {
        hit_ids[0] = x;
        hit_ids[1] = y;
    }
    void setADC_time_value(Int_t time_bin, Short_t ADC_value) { ADC_time_values[time_bin] = ADC_value; }
    void set_pos(Int_t time_bin, Float_t x_pos, Float_t y_pos,
                 Float_t z_pos)  //{arr_pos[0][time_bin] = (Short_t)(100.0*x_pos); arr_pos[1][time_bin] = (Short_t)(100.0*y_pos);
                                 // arr_pos[2][time_bin] = (Short_t)(100.0*z_pos);}
    {
        // Subtract a vector of length 200 cm from digit positions to reduce the range. Usually ~300-360 cm, now 100-160 cm.
        // This is important to not get out of Short_t range: 8 byte: 2^16/2 = 32768. 360*100 > 32768
        TVector3 vec_200cm(x_pos, y_pos, z_pos);
        // printf("pos: {%4.3f, %4.3f, %4.3f} \n",x_pos,y_pos,z_pos);
        if (vec_200cm.Mag() > 0.0) {
            vec_200cm *= 200.0 / vec_200cm.Mag();
            // printf("vector length: %4.3f \n",vec_200cm.Mag());
            arr_pos[0][time_bin] = (Short_t)(100.0 * (x_pos - vec_200cm.X()));
            arr_pos[1][time_bin] = (Short_t)(100.0 * (y_pos - vec_200cm.Y()));
            arr_pos[2][time_bin] = (Short_t)(100.0 * (z_pos - vec_200cm.Z()));
        } else {
            arr_pos[0][time_bin] = 0.0;
            arr_pos[1][time_bin] = 0.0;
            arr_pos[2][time_bin] = 0.0;
        }

        // std::cout << "x_pos: " << x_pos << std::endl;
        // std::cout << "y_pos: " << y_pos << std::endl;
        // std::cout << "z_pos: " << z_pos << std::endl;
        // std::cout << "x_pos-200: " << arr_pos[0][time_bin] << std::endl;
        // std::cout << "y_pos-200: " << arr_pos[1][time_bin] << std::endl;
        // std::cout << "z_pos-200: " << arr_pos[2][time_bin] << std::endl;
    }

    void setdca_to_track(Float_t f, Float_t f_x, Float_t f_y, Float_t f_z) {
        dca_to_track = (Short_t)(100.0 * f);
        dca_x = (Short_t)(100.0 * f_x);
        dca_y = (Short_t)(100.0 * f_y);
        dca_z = (Short_t)(100.0 * f_z);
    }
    void setImpactAngle(Float_t f) { ImpactAngle = (Short_t)(100.0 * f); }

    // getters
    UInt_t gethit_ids(Int_t i) const { return hit_ids[i]; }
    Short_t getADC_time_value(Int_t time_bin) const { return ADC_time_values[time_bin]; }
    Float_t getdca_to_track() const { return ((Float_t)dca_to_track) / 100.0; }
    Float_t getdca_x() const { return ((Float_t)dca_x) / 100.0; }
    Float_t getdca_y() const { return ((Float_t)dca_y) / 100.0; }
    Float_t getdca_z() const { return ((Float_t)dca_z) / 100.0; }
    Float_t getImpactAngle() const { return ((Float_t)ImpactAngle) / 100.0; }
    Int_t get_sector() const { return (Int_t)((Float_t)hit_ids[0] / 144.0); }
    Int_t get_column() const { return (hit_ids[0] % 144); }
    Int_t get_layer() const { return (Int_t)((Float_t)hit_ids[1] / (16.0 * 5.0)); }
    Int_t get_stack() const { return ((Int_t)((Float_t)hit_ids[1] / 16.0)) % 5; }
    Int_t get_row() const { return (hit_ids[1] % 16); }
    Int_t get_detector(Int_t layer, Int_t stack, Int_t sector) const { return (layer + stack * 6 + sector * 6 * 5); }
    Float_t get_pos(Int_t time_bin, Int_t index) const  //{ return ((Float_t)arr_pos[index][time_bin])/100.0; }
    {
        TVector3 vec_200cm((Float_t)arr_pos[0][time_bin], (Float_t)arr_pos[1][time_bin], (Float_t)arr_pos[2][time_bin]);

        // now vec_200cm is again in cm
        if (vec_200cm.Mag() > 0.0) {
            vec_200cm *= 200.0 / vec_200cm.Mag();  // when we get positions XYZ-200cm
            return (Float_t)(arr_pos[index][time_bin] / 100.0 + vec_200cm(index));
        } else {
            return 0.0;
        }
    }

    // x = TRD_col + TRD_sec*144;
    // y = TRD_row + TRD_stack*16 + TRD_lay*16*5;

    ClassDef(Ali_AS_TRD_digit, 1);
};

#endif
