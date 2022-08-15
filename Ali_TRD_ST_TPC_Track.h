#ifndef __ALI_TRD_ST_TPC_TRACK_H__
#define __ALI_TRD_ST_TPC_TRACK_H__

#include "TLorentzVector.h"
#include "TObject.h"

class Ali_TRD_ST_TPC_Track : public TObject {
   private:
    // Track properties
    Float_t nsigma_e_TPC;         // nsigma dE/dx of particle
    Float_t nsigma_e_TOF;         // nsigma TOF of particle
    Float_t nsigma_pi_TPC;        // nsigma dE/dx of particle
    Float_t nsigma_pi_TOF;        // nsigma TOF of particle
    Float_t nsigma_K_TPC;         // nsigma dE/dx of particle
    Float_t nsigma_K_TOF;         // nsigma TOF of particle
    Float_t nsigma_p_TPC;         // nsigma dE/dx of particle
    Float_t nsigma_p_TOF;         // nsigma TOF of particle
    Float_t TRD_signal;           // TRD PID signal
    Float_t TRDsumADC;            // TRD electron PID probability
    Float_t dca;                  // distance of closest approach of particle A
    TLorentzVector TLV_part;      // Lorentz vector properties of this particle
    UShort_t NTPCcls;             // Number of TPC clusters
    UShort_t NTRDcls;             // Number of TRD clusters
    UShort_t NITScls;             // Number of TRD clusters
    UShort_t status;              // status of track: bit 0: ITS refit, bit1: TPC refit
    Float_t TPCchi2;              // TPC chi2
    Float_t impact_angle_on_TRD;  // Track impact angle on TRD
    Float_t TPCdEdx;              // Energy loss information of TPC
    Float_t TOFsignal;            // Time-of-flight
    Float_t Track_length;         // length of track
    Int_t MC_label;               // index/label of corresponding MC particle, borquez edit

    Float_t aliHelix_params[6];
    Float_t aliHelix_TRD_params[6];

   public:
    // constructor
    Ali_TRD_ST_TPC_Track()
        : nsigma_e_TPC(-3),
          nsigma_e_TOF(-3),
          nsigma_pi_TPC(-3),
          nsigma_pi_TOF(-3),
          nsigma_K_TPC(-3),
          nsigma_K_TOF(-3),
          nsigma_p_TPC(-3),
          nsigma_p_TOF(-3),
          TRD_signal(-3),
          TRDsumADC(-3),
          dca(-3),
          TLV_part(),
          NTPCcls(-3),
          NTRDcls(-3),
          NITScls(-3),
          status(-3),
          TPCchi2(-3),
          impact_angle_on_TRD(-3),
          TPCdEdx(-3),
          TOFsignal(-3),
          Track_length(-3),
          MC_label(0),
          aliHelix_params(),
          aliHelix_TRD_params() {}

    // destructor
    ~Ali_TRD_ST_TPC_Track() {}

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
    void setimpact_angle_on_TRD(Float_t f) { impact_angle_on_TRD = f; }
    void setTPCdEdx(Float_t f) { TPCdEdx = f; }
    void setTOFsignal(Float_t f) { TOFsignal = f; }
    void setTrack_length(Float_t f) { Track_length = f; }
    void setMC_label(Float_t f) { MC_label = f; }  // borquez edit

    void setHelix(Float_t a, Float_t b, Float_t c, Float_t d, Float_t e, Float_t f) {
        aliHelix_params[0] = a;
        aliHelix_params[1] = b;
        aliHelix_params[2] = c;
        aliHelix_params[3] = d;
        aliHelix_params[4] = e;
        aliHelix_params[5] = f;
    }

    void setHelix_TRD(Float_t a, Float_t b, Float_t c, Float_t d, Float_t e, Float_t f) {
        aliHelix_TRD_params[0] = a;
        aliHelix_TRD_params[1] = b;
        aliHelix_TRD_params[2] = c;
        aliHelix_TRD_params[3] = d;
        aliHelix_TRD_params[4] = e;
        aliHelix_TRD_params[5] = f;
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
    Float_t getimpact_angle_on_TRD() const { return impact_angle_on_TRD; }
    Float_t getTPCdEdx() const { return TPCdEdx; }
    Float_t getTOFsignal() const { return TOFsignal; }
    Float_t getTrack_length() const { return Track_length; }
    Float_t getHelix_param(Int_t i_param) const { return aliHelix_params[i_param]; }
    Float_t getHelix_TRD_param(Int_t i_param) const { return aliHelix_TRD_params[i_param]; }
    Int_t getMC_label() { return MC_label; }  // borquez edit

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

    ClassDef(Ali_TRD_ST_TPC_Track, 1);  // A simple track of a particle
};

#endif
