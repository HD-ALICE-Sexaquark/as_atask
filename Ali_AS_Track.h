#ifndef __ALI_AS_TRACK_H__
#define __ALI_AS_TRACK_H__

#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TObject.h"

#include "Ali_AS_TRD_digit.h"
#include "Ali_AS_offline_Tracklet.h"

class Ali_AS_Track : public TObject {
   private:
    // Track properties
    Float_t nsigma_e_TPC;     // nsigma dE/dx of particle
    Float_t nsigma_e_TOF;     // nsigma TOF of particle
    Float_t nsigma_pi_TPC;    // nsigma dE/dx of particle
    Float_t nsigma_pi_TOF;    // nsigma TOF of particle
    Float_t nsigma_K_TPC;     // nsigma dE/dx of particle
    Float_t nsigma_K_TOF;     // nsigma TOF of particle
    Float_t nsigma_p_TPC;     // nsigma dE/dx of particle
    Float_t nsigma_p_TOF;     // nsigma TOF of particle
    Float_t TRD_signal;       // TRD PID signal
    Float_t TRDsumADC;        // TRD electron PID probability
    Float_t dca;              // distance of closest approach of particle A
    TLorentzVector TLV_part;  // Lorentz vector properties of this particle
    UShort_t NTPCcls;         // Number of TPC clusters
    UShort_t NTRDcls;         // Number of TRD clusters
    UShort_t NITScls;         // Number of TRD clusters
    UShort_t status;          // status of track: bit 0: ITS refit, bit1: TPC refit
    Float_t TPCchi2;          // TPC chi2
    ULong64_t TRD_ADC_time_layer[6];
    Float_t impact_angle_on_TRD;  // Track impact angle on TRD
    Float_t TPCdEdx;              // Energy loss information of TPC
    Float_t TOFsignal;            // Time-of-flight
    Float_t Track_length;         // length of track
    Float_t aliHelix_params[9];
    Float_t aliHelix_TRD_params[9];
    Int_t MC_label;  // index/label of corresponding MC particle, borquez edit

    UShort_t fNumTRDdigits;         // number of TRD digits for this track
    UShort_t fNumOfflineTracklets;  // number of offline tracklets

    TClonesArray* fTRD_digits;
    TClonesArray* fOfflineTracklets;

   public:
    Ali_AS_Track()
        : nsigma_e_TPC(-1),
          nsigma_e_TOF(-1),
          nsigma_pi_TPC(-1),
          nsigma_pi_TOF(-1),
          nsigma_K_TPC(-1),
          nsigma_K_TOF(-1),
          nsigma_p_TPC(-1),
          nsigma_p_TOF(-1),
          TRD_signal(-1),
          TRDsumADC(-1),
          dca(-1),
          TLV_part(),
          NTPCcls(-1),
          NTRDcls(-1),
          NITScls(-1),
          status(-1),
          TPCchi2(-1),
          TRD_ADC_time_layer(),
          impact_angle_on_TRD(-1),
          TPCdEdx(-1),
          TOFsignal(-1),
          Track_length(-1),
          aliHelix_params(),
          aliHelix_TRD_params(),
          fNumTRDdigits(0),
          fNumOfflineTracklets(0),
          fTRD_digits(),
          fOfflineTracklets(),
          MC_label(0)  // borquez edit
    {
        // constructor
        fTRD_digits = new TClonesArray("Ali_AS_TRD_digit", 10);
        fOfflineTracklets = new TClonesArray("Ali_AS_offline_Tracklet", 10);
    }

    ~Ali_AS_Track() {
        // destructor
        delete fTRD_digits;
        delete fOfflineTracklets;
        fTRD_digits = NULL;
        fOfflineTracklets = NULL;
    }

    // setters
    void setnsigma_e_TPC(Float_t f) { nsigma_e_TPC = f; }
    void setnsigma_e_TOF(Float_t f) { nsigma_e_TOF = f; }
    void setnsigma_pi_TPC(Float_t f) { nsigma_pi_TPC = f; }
    void setnsigma_pi_TOF(Float_t f) { nsigma_pi_TOF = f; }
    void setnsigma_K_TPC(Float_t f) { nsigma_K_TPC = f; }
    void setnsigma_K_TOF(Float_t f) { nsigma_K_TOF = f; }
    void setnsigma_p_TPC(Float_t f) { nsigma_p_TPC = f; }
    void setnsigma_p_TOF(Float_t f) { nsigma_p_TOF = f; }
    void setTRDSignal(Float_t f) { TRD_signal = f; }
    void setTRDsumADC(Float_t f) { TRDsumADC = f; }
    void setdca(Float_t f) { dca = f; }
    void set_TLV_part(TLorentzVector tlv) { TLV_part = tlv; }
    void setNTPCcls(UShort_t s) { NTPCcls = s; }
    void setNTRDcls(UShort_t s) { NTRDcls = s; }
    void setNITScls(UShort_t s) { NITScls = s; }
    void setStatus(UShort_t s) { status = s; }
    void setTPCchi2(Float_t f) { TPCchi2 = f; }
    void setTRD_layer(Int_t i_layer, ULong64_t l) { TRD_ADC_time_layer[i_layer] = l; }
    void setimpact_angle_on_TRD(Float_t f) { impact_angle_on_TRD = f; }
    void setTPCdEdx(Float_t f) { TPCdEdx = f; }
    void setTOFsignal(Float_t f) { TOFsignal = f; }
    void setTrack_length(Float_t f) { Track_length = f; }
    void setMC_label(Float_t f) { MC_label = f; }  // borquez edit

    void setHelix(Float_t a, Float_t b, Float_t c, Float_t d, Float_t e, Float_t f, Float_t g, Float_t h, Float_t i) {
        aliHelix_params[0] = a;
        aliHelix_params[1] = b;
        aliHelix_params[2] = c;
        aliHelix_params[3] = d;
        aliHelix_params[4] = e;
        aliHelix_params[5] = f;
        aliHelix_params[6] = g;
        aliHelix_params[7] = h;
        aliHelix_params[8] = i;
    }

    void setHelix_TRD(Float_t a, Float_t b, Float_t c, Float_t d, Float_t e, Float_t f, Float_t g, Float_t h, Float_t i) {
        aliHelix_TRD_params[0] = a;
        aliHelix_TRD_params[1] = b;
        aliHelix_TRD_params[2] = c;
        aliHelix_TRD_params[3] = d;
        aliHelix_TRD_params[4] = e;
        aliHelix_TRD_params[5] = f;
        aliHelix_TRD_params[6] = g;
        aliHelix_TRD_params[7] = h;
        aliHelix_TRD_params[8] = i;
    }

    // getters
    Float_t getnsigma_e_TPC() const { return nsigma_e_TPC; }
    Float_t getnsigma_e_TOF() const { return nsigma_e_TOF; }
    Float_t getnsigma_pi_TPC() const { return nsigma_pi_TPC; }
    Float_t getnsigma_pi_TOF() const { return nsigma_pi_TOF; }
    Float_t getnsigma_K_TPC() const { return nsigma_K_TPC; }
    Float_t getnsigma_K_TOF() const { return nsigma_K_TOF; }
    Float_t getnsigma_p_TPC() const { return nsigma_p_TPC; }
    Float_t getnsigma_p_TOF() const { return nsigma_p_TOF; }
    Float_t getTRDSignal() const { return TRD_signal; }
    Float_t getTRDsumADC() const { return TRDsumADC; }
    Float_t getdca() const { return dca; }
    TLorentzVector get_TLV_part() const { return TLV_part; }
    UShort_t getNTPCcls() const { return NTPCcls; }
    UShort_t getNTRDcls() const { return NTRDcls; }
    UShort_t getNITScls() const { return NITScls; }
    UShort_t getStatus() const { return status; }
    Float_t getTPCchi2() const { return TPCchi2; }
    ULong64_t getTRD_layer(Int_t i_layer) const { return TRD_ADC_time_layer[i_layer]; }
    Float_t getimpact_angle_on_TRD() const { return impact_angle_on_TRD; }
    Float_t getTPCdEdx() const { return TPCdEdx; }
    Float_t getTOFsignal() const { return TOFsignal; }
    Float_t getTrack_length() const { return Track_length; }
    Float_t getHelix_param(Int_t i_param) const { return aliHelix_params[i_param]; }
    Float_t getHelix_TRD_param(Int_t i_param) const { return aliHelix_TRD_params[i_param]; }
    Int_t getMC_label() { return MC_label; }  // borquez edit
    Float_t getTRD_ADC(Int_t i_layer, Int_t i_time_bin) const {
        if (i_layer < 0 || i_layer > 5 || i_time_bin < 0 || i_time_bin > 7) {
            // out of range
            return -1;
        }
        ULong64_t TRD_value = 0;
        for (Int_t i_bit = 0; i_bit < 8; i_bit++) {
            // One single time slice 8 bit = 256
            Int_t bitcheck = i_bit + 8 * i_time_bin;                           // range: 0..63 = 64 bit = 8 byte = Long64_t
            Int_t bit_status = (TRD_ADC_time_layer[i_layer] >> bitcheck) & 1;  // check bit bitcheck
            if (bit_status) {
                // setting bit i_bit to 1
                TRD_value |= (ULong64_t)1 << i_bit;
            }
        }
        Float_t TRD_value_decode = (Float_t)TRD_value * 100.0;  // * TRD_ADC_bin_width
        return TRD_value_decode;
    }
    UShort_t getNumOfflineTracklets() const { return fNumOfflineTracklets; }
    UShort_t getNumTRD_digits() const { return fNumTRDdigits; }
    Ali_AS_TRD_digit* getTRD_digit(UShort_t i) const { return i < fNumTRDdigits ? (Ali_AS_TRD_digit*)((*fTRD_digits)[i]) : NULL; }
    Ali_AS_offline_Tracklet* getOfflineTracklet(UShort_t i) const {
        return i < fNumOfflineTracklets ? (Ali_AS_offline_Tracklet*)((*fOfflineTracklets)[i]) : NULL;
    }

    Int_t HasITShit_on_layer(Int_t ilayer) {
        // ITShit
        return ((NITScls >> ilayer) & 1);
    }

    Ali_AS_offline_Tracklet* createOfflineTracklet() {
        if (fNumOfflineTracklets == fOfflineTracklets->GetSize()) {
            fOfflineTracklets->Expand(fNumOfflineTracklets + 10);
        }
        if (fNumOfflineTracklets >= 65000) {
            Fatal("Ali_AS_Event::createOfflineTracklet()", "ERROR: Too many tracklets (>65000)!");
            exit(2);
        }

        new ((*fOfflineTracklets)[fNumOfflineTracklets++]) Ali_AS_offline_Tracklet;
        return (Ali_AS_offline_Tracklet*)((*fOfflineTracklets)[fNumOfflineTracklets - 1]);
    }

    void clearOfflineTrackletList() {
        fNumOfflineTracklets = 0;
        fOfflineTracklets->Clear();
    }

    Ali_AS_TRD_digit* createTRD_digit() {
        if (fNumTRDdigits == fTRD_digits->GetSize()) {
            fTRD_digits->Expand(fNumTRDdigits + 10);
        }
        if (fNumTRDdigits >= 50000) {
            Fatal("Ali_AS_Track::createTRD_digit()", "ERROR: Too many TRD digits (>50000)!");
            exit(2);
        }

        new ((*fTRD_digits)[fNumTRDdigits++]) Ali_AS_TRD_digit;
        return (Ali_AS_TRD_digit*)((*fTRD_digits)[fNumTRDdigits - 1]);
    }

    void clearTRD_digit_list() {
        fNumTRDdigits = 0;
        fTRD_digits->Clear();
    }

    ClassDef(Ali_AS_Track, 1);  // A simple track of a particle
};

#endif
