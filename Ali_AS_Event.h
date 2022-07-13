#ifndef __ALI_AS_EVENT_H__
#define __ALI_AS_EVENT_H__

#include "TClonesArray.h"
#include "TObject.h"
#include "TString.h"

#include "Ali_AS_TRD_digit.h"
#include "Ali_AS_Track.h"
#include "Ali_AS_Tracklet.h"

class Ali_AS_Event : public TObject {
   private:
    Int_t eventNumber;
    Float_t x;               // Event vertex x
    Float_t y;               // Event vertex y
    Float_t z;               // Event vertex z
    Int_t id;                // Run id
    Int_t N_tracks;          // total number of tracks
    Int_t N_TRD_tracklets;   // total number of TRD tracklets
    Float_t cent_class_ZNA;  // ZDC neutral A
    Float_t cent_class_ZNC;  // ZDC neutral C
    Float_t cent_class_V0A;  // V0 A
    Float_t cent_class_V0C;  // V0 C
    Float_t cent_class_V0M;  // V0 average
    Float_t cent_class_CL0;  // clusters in layer 0
    Float_t cent_class_CL1;  // clusters in layer 1
    Float_t cent_class_SPD;  // SPD
    Float_t cent_class_V0MEq;
    Float_t cent_class_V0AEq;
    Float_t cent_class_V0CEq;
    Float_t ADC_sum_det[540];

    Float_t BeamIntAA;  // ZDC coincidence rate
    Float_t T0zVertex;  // z-vertex position from VPD

    TString TriggerWord;  // Trigger word

    UShort_t fNumTracks;  // number of tracks in event
    Int_t fNumTracklets;  // number of tracks in event
    Int_t fNumTRDdigits;  // number of TRD digits for this track

    TClonesArray* fTracks;
    TClonesArray* fTracklets;
    TClonesArray* fTRD_digits;

   public:
    Ali_AS_Event()
        : eventNumber(-1),
          x(-1),
          y(-1),
          z(-1),
          id(-1),
          N_tracks(0),
          N_TRD_tracklets(0),
          cent_class_ZNA(0),
          cent_class_ZNC(0),
          cent_class_V0A(0),
          cent_class_V0C(0),
          cent_class_V0M(0),
          cent_class_CL0(0),
          cent_class_CL1(0),
          cent_class_SPD(0),
          cent_class_V0MEq(0),
          cent_class_V0AEq(0),
          cent_class_V0CEq(0),
          ADC_sum_det(),
          BeamIntAA(-1),
          T0zVertex(-1),
          TriggerWord(),
          fNumTracks(0),
          fNumTracklets(0),
          fNumTRDdigits(0),
          fTracks(),
          fTracklets(),
          fTRD_digits()

    {
        // constructor
        fTRD_digits = new TClonesArray("Ali_AS_TRD_digit", 10);
        fTracks = new TClonesArray("Ali_AS_Track", 10);
        fTracklets = new TClonesArray("Ali_AS_Tracklet", 10);
    }

    ~Ali_AS_Event() {
        // destructor
        delete fTRD_digits;
        fTRD_digits = NULL;
        delete fTracks;
        fTracks = NULL;
        delete fTracklets;
        fTracklets = NULL;
    }

    // setters
    void setEventNumber(Int_t e) { eventNumber = e; }
    void setx(Float_t r) { x = r; }
    void sety(Float_t r) { y = r; }
    void setz(Float_t r) { z = r; }
    void setid(Int_t r) { id = r; }
    void setN_tracks(Int_t r) { N_tracks = r; }
    void setN_TRD_tracklets(Int_t r) { N_TRD_tracklets = r; }
    void setcent_class_ZNA(Float_t r) { cent_class_ZNA = r; }
    void setcent_class_ZNC(Float_t r) { cent_class_ZNC = r; }
    void setcent_class_V0A(Float_t r) { cent_class_V0A = r; }
    void setcent_class_V0C(Float_t r) { cent_class_V0C = r; }
    void setcent_class_V0M(Float_t r) { cent_class_V0M = r; }
    void setcent_class_CL0(Float_t r) { cent_class_CL0 = r; }
    void setcent_class_CL1(Float_t r) { cent_class_CL1 = r; }
    void setcent_class_SPD(Float_t r) { cent_class_SPD = r; }
    void setcent_class_V0MEq(Float_t r) { cent_class_V0MEq = r; }
    void setcent_class_V0AEq(Float_t r) { cent_class_V0AEq = r; }
    void setcent_class_V0CEq(Float_t r) { cent_class_V0CEq = r; }
    void setBeamIntAA(Float_t r) { BeamIntAA = r; }
    void setT0zVertex(Float_t r) { T0zVertex = r; }
    void setTriggerWord(TString s) { TriggerWord = s; }
    void setADC_sum_det(Int_t i_det, Float_t r) { ADC_sum_det[i_det] = r; }

    // getters
    Int_t getEventNumber() const { return eventNumber; }
    Float_t getx() const { return x; }
    Float_t gety() const { return y; }
    Float_t getz() const { return z; }
    Int_t getid() const { return id; }
    Int_t getN_tracks() const { return N_tracks; }
    Int_t getN_TRD_tracklets() const { return N_TRD_tracklets; }
    Float_t getcent_class_ZNA() const { return cent_class_ZNA; }
    Float_t getcent_class_ZNC() const { return cent_class_ZNC; }
    Float_t getcent_class_V0A() const { return cent_class_V0A; }
    Float_t getcent_class_V0C() const { return cent_class_V0C; }
    Float_t getcent_class_V0M() const { return cent_class_V0M; }
    Float_t getcent_class_CL0() const { return cent_class_CL0; }
    Float_t getcent_class_CL1() const { return cent_class_CL1; }
    Float_t getcent_class_SPD() const { return cent_class_SPD; }
    Float_t getcent_class_V0MEq() const { return cent_class_V0MEq; }
    Float_t getcent_class_V0AEq() const { return cent_class_V0AEq; }
    Float_t getcent_class_V0CEq() const { return cent_class_V0CEq; }
    Float_t getBeamIntAA() const { return BeamIntAA; }
    Float_t getT0zVertex() const { return T0zVertex; }
    TString getTriggerWord() const { return TriggerWord; }
    ULong64_t getADC_sum_det(Int_t i_det) const { return ADC_sum_det[i_det]; }
    Int_t getNumTRD_digits() const { return fNumTRDdigits; }
    Ali_AS_TRD_digit* getTRD_digit(Int_t i) const { return i < fNumTRDdigits ? (Ali_AS_TRD_digit*)((*fTRD_digits)[i]) : NULL; }
    UShort_t getNumTracks() const { return fNumTracks; }
    Ali_AS_Track* getTrack(UShort_t i) const { return i < fNumTracks ? (Ali_AS_Track*)((*fTracks)[i]) : NULL; }
    Int_t getNumTracklets() const { return fNumTracklets; }
    Ali_AS_Tracklet* getTracklet(Int_t i) const { return i < fNumTracklets ? (Ali_AS_Tracklet*)((*fTracklets)[i]) : NULL; }

    Ali_AS_TRD_digit* createTRD_digit() {
        if (fNumTRDdigits == fTRD_digits->GetSize()) {
            fTRD_digits->Expand(fNumTRDdigits + 10);
        }
        if (fNumTRDdigits >= 650000) {
            Fatal("Ali_AS_Track::createTRD_digit()", "ERROR: Too many TRD digits (>650000)!");
            exit(2);
        }

        new ((*fTRD_digits)[fNumTRDdigits++]) Ali_AS_TRD_digit;
        return (Ali_AS_TRD_digit*)((*fTRD_digits)[fNumTRDdigits - 1]);
    }

    void clearTRD_digit_list() {
        fNumTRDdigits = 0;
        fTRD_digits->Clear();
    }

    Ali_AS_Track* createTrack() {
        if (fNumTracks == fTracks->GetSize()) {
            fTracks->Expand(fNumTracks + 10);
        }
        if (fNumTracks >= 50000) {
            Fatal("Ali_AS_Event::createTrack()", "ERROR: Too many tracks (>50000)!");
            exit(2);
        }

        new ((*fTracks)[fNumTracks++]) Ali_AS_Track;
        return (Ali_AS_Track*)((*fTracks)[fNumTracks - 1]);
    }

    void clearTrackList() {
        fNumTracks = 0;
        fTracks->Clear();
    }

    Ali_AS_Tracklet* createTracklet()  // online tracklet
    {
        if (fNumTracklets == fTracklets->GetSize()) {
            fTracklets->Expand(fNumTracklets + 10);
        }
        if (fNumTracklets >= 650000) {
            Fatal("Ali_AS_Event::createTracklet()", "ERROR: Too many tracklets (>65000)!");
            exit(2);
        }

        new ((*fTracklets)[fNumTracklets++]) Ali_AS_Tracklet;
        return (Ali_AS_Tracklet*)((*fTracklets)[fNumTracklets - 1]);
    }

    void clearTrackletList() {
        fNumTracklets = 0;
        fTracklets->Clear();
    }

    ClassDef(Ali_AS_Event, 1);  // A simple event compiled of tracks
};

#endif
