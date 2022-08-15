#ifndef __ALI_TRD_ST_V0_H__
#define __ALI_TRD_ST_V0_H__

#include "TObject.h"

class Ali_TRD_ST_V0 : public TObject {
   private:
    Double_t fPos[3];   // Event vertex x,y,z
    Double_t fNmom[3];  // momentum of negative charged particle
    Double_t fPmom[3];  // momentum of positive charged particle
    Double_t fEffMass;

    Int_t fNidx;  // index of the negative daughter
    Int_t fPidx;  // index of the positive daughter

    Bool_t fOnFlyStatus;  // if kTRUE, then this V0 is recontructed
                          // "on fly" during the tracking

   public:
    Ali_TRD_ST_V0() : fOnFlyStatus(kFALSE), fEffMass(0), fNidx(0), fPidx(0) {
        for (Int_t i = 0; i < 3; i++) {
            fPos[i] = 0.;
            fNmom[i] = 0.;
            fPmom[i] = 0.;
        }
    }

    Ali_TRD_ST_V0(const Double_t *pos, const Double_t *nmom, const Double_t *pmom, const Bool_t onflystatus, const Double_t effmass,
                  Int_t i1, Int_t i2)
        : fOnFlyStatus(onflystatus), fEffMass(effmass), fNidx(i1), fPidx(i2) {
        for (Int_t i = 0; i < 3; i++) {
            fPos[i] = pos[i];
            fNmom[i] = nmom[i];
            fPmom[i] = pmom[i];
        }
    }

    ~Ali_TRD_ST_V0() {}

    Double_t Px() const { return fNmom[0] + fPmom[0]; }
    Double_t Py() const { return fNmom[1] + fPmom[1]; }
    Double_t Pz() const { return fNmom[2] + fPmom[2]; }
    Double_t Pt() const { return TMath::Sqrt(Px() * Px() + Py() * Py()); }
    Double_t P() const { return TMath::Sqrt(Px() * Px() + Py() * Py() + Pz() * Pz()); }
    Bool_t PxPyPz(Double_t p[3]) const {
        p[0] = Px();
        p[1] = Py();
        p[2] = Pz();
        return kTRUE;
    }
    Double_t Xv() const { return fPos[0]; }
    Double_t Yv() const { return fPos[1]; }
    Double_t Zv() const { return fPos[2]; }
    Bool_t XvYvZv(Double_t x[3]) const {
        x[0] = Xv();
        x[1] = Yv();
        x[2] = Zv();
        return kTRUE;
    }
    Double_t OneOverPt() const { return (Pt() != 0.) ? 1. / Pt() : -999.; }
    Double_t Phi() const { return TMath::Pi() + TMath::ATan2(-Py(), -Px()); }
    Double_t Theta() const { return 0.5 * TMath::Pi() - TMath::ATan(Pz() / (Pt() + 1.e-13)); }
    Double_t E() const;  // default is KOs but can be changed via ChangeMassHypothesis (defined in the .cxx)
    Double_t M() const { return GetEffMass(); }
    Double_t Eta() const { return 0.5 * TMath::Log((P() + Pz()) / (P() - Pz() + 1.e-13)); }
    Double_t Y() const;
    Short_t Charge() const { return 0; }
    Int_t GetLabel() const { return -1; }  // temporary
    Double_t GetEffMass() const { return fEffMass; }
    void SetEffMass(Double_t eff_mass) { fEffMass = eff_mass; }

    void GetPxPyPz(Double_t &px, Double_t &py, Double_t &pz) const {
        // This function returns V0's momentum
        px = fNmom[0] + fPmom[0];
        py = fNmom[1] + fPmom[1];
        pz = fNmom[2] + fPmom[2];
    }

    void GetNPxPyPz(Double_t &px, Double_t &py, Double_t &pz) const {
        px = fNmom[0];
        py = fNmom[1];
        pz = fNmom[2];
    }
    void SetNPxPyPz(Double_t px, Double_t py, Double_t pz) {
        fNmom[0] = px;
        fNmom[1] = py;
        fNmom[2] = pz;
    }

    void GetPPxPyPz(Double_t &px, Double_t &py, Double_t &pz) const {
        px = fPmom[0];
        py = fPmom[1];
        pz = fPmom[2];
    }
    void SetPPxPyPz(Double_t px, Double_t py, Double_t pz) {
        fPmom[0] = px;
        fPmom[1] = py;
        fPmom[2] = pz;
    }

    void GetXYZ(Double_t &x, Double_t &y, Double_t &z) const {
        // This function returns V0's position
        x = fPos[0];
        y = fPos[1];
        z = fPos[2];
    }
    void SetXYZ(Double_t x, Double_t y, Double_t z) {
        fPos[0] = x;
        fPos[1] = y;
        fPos[2] = z;
    }

    Int_t GetNindex() const { return fNidx; }
    void SetNindex(Int_t idx) { fNidx = idx; }
    Int_t GetPindex() const { return fPidx; }
    void SetPindex(Int_t idx) { fPidx = idx; }

    void SetOnFlyStatus(Bool_t status) { fOnFlyStatus = status; }
    Bool_t GetOnFlyStatus() const { return fOnFlyStatus; }

    ClassDef(Ali_TRD_ST_V0, 1);
};

#endif
