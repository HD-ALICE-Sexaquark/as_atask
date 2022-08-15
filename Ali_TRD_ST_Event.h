#ifndef __ALI_TRD_ST_EVENT_H__
#define __ALI_TRD_ST_EVENT_H__

#include "TClonesArray.h"
#include "TObject.h"
#include "TString.h"

#include "Ali_MC_particle.h"
#include "Ali_TRD_ST_TOF_hit.h"
#include "Ali_TRD_ST_TPC_Track.h"
#include "Ali_TRD_ST_Tracklets.h"
#include "Ali_TRD_ST_V0.h" // hoppner edit

class Ali_TRD_ST_Event : public TObject {
   private:
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
    Short_t flag_Data_MC;  // 0 = data, 1 = MC

    Float_t BeamIntAA;  // ZDC coincidence rate
    Float_t T0zVertex;  // z-vertex position from VPD

    TString TriggerWord;       // Trigger word
    UShort_t N_TRD_time_bins;  // Number of TRD time bins

    UShort_t fNumTracks;    // number of tracks in event
    Int_t fNumTracklets;    // number of tracks in event
    Int_t fNumMCparticles;  // number of MC particles
    Int_t fNumTOFhits;      // number of TOF hits
    Int_t fNumV0s;          // number of V0s hoppner edit

    TClonesArray* fTracks;
    TClonesArray* fTracklets;
    TClonesArray* fMCparticles;
    TClonesArray* fTOFhits;
    TClonesArray* fV0s; // hoppner edit

   public:
    Ali_TRD_ST_Event()
        : x(-3),
          y(-3),
          z(-3),
          id(-3),
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
          flag_Data_MC(0),
          BeamIntAA(-3),
          T0zVertex(-3),
          TriggerWord(),
          N_TRD_time_bins(0),
          fNumTracks(0),
          fNumTracklets(0),
          fNumV0s(0),  // hoppner edit
          fNumMCparticles(0),
          fNumTOFhits(0),
          fTracks(),
          fTracklets(),
          fV0s(),
          fMCparticles(),
          fTOFhits() {
        // constructor
        fTracks = new TClonesArray("Ali_TRD_ST_TPC_Track", 10);
        fTracklets = new TClonesArray("Ali_TRD_ST_Tracklets", 10);
        fV0s = new TClonesArray("Ali_TRD_ST_V0", 10);  // hoppner edit
        fMCparticles = new TClonesArray("Ali_MC_particle", 10);
        fTOFhits = new TClonesArray("Ali_TRD_ST_TOF_hit", 10);
    }

    ~Ali_TRD_ST_Event() {
        // destructor
        delete fTracks;
        fTracks = NULL;
        delete fTracklets;
        fTracklets = NULL;
        delete fMCparticles;
        fMCparticles = NULL;
        delete fTOFhits;
        fTOFhits = NULL;
        delete fV0s;  // hoppner edit
        fV0s = NULL;  // hoppner edit
    }

    // setters
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
    void setcent_class_V0AEq(Float_t r) { cent_class_V0AEq = r; }
    void setcent_class_V0MEq(Float_t r) { cent_class_V0MEq = r; }
    void setcent_class_V0CEq(Float_t r) { cent_class_V0CEq = r; }
    void set_flag_Data_MC(UShort_t s) { flag_Data_MC = s; }
    void setBeamIntAA(Float_t r) { BeamIntAA = r; }
    void setT0zVertex(Float_t r) { T0zVertex = r; }
    void setTriggerWord(TString s) { TriggerWord = s; }
    void setN_TRD_time_bins(UShort_t s) { N_TRD_time_bins = s; }

    // getters
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
    Float_t getcent_class_V0AEq() const { return cent_class_V0AEq; }
    Float_t getcent_class_V0MEq() const { return cent_class_V0MEq; }
    Float_t getcent_class_V0CEq() const { return cent_class_V0CEq; }
    UShort_t get_flag_Data_MC() const { return flag_Data_MC; }
    Float_t getBeamIntAA() const { return BeamIntAA; }
    Float_t getT0zVertex() const { return T0zVertex; }
    TString getTriggerWord() const { return TriggerWord; }
    Float_t getN_TRD_time_bins() const { return N_TRD_time_bins; }

    Int_t getMCparticles() const { return fNumMCparticles; }
    Ali_MC_particle* getMCparticle(Int_t i) const { return i < fNumMCparticles ? (Ali_MC_particle*)((*fMCparticles)[i]) : NULL; }

    Int_t getNumTOFhits() const { return fNumTOFhits; }
    Ali_TRD_ST_TOF_hit* getTOFhit(Int_t i) const { return i < fNumTOFhits ? (Ali_TRD_ST_TOF_hit*)((*fTOFhits)[i]) : NULL; }

    UShort_t getNumTracks() const { return fNumTracks; }
    Ali_TRD_ST_TPC_Track* getTrack(UShort_t i) const { return i < fNumTracks ? (Ali_TRD_ST_TPC_Track*)((*fTracks)[i]) : NULL; }

    Int_t getNumTracklets() const { return fNumTracklets; }
    Ali_TRD_ST_Tracklets* getTracklet(Int_t i) const { return i < fNumTracklets ? (Ali_TRD_ST_Tracklets*)((*fTracklets)[i]) : NULL; }

    // hoppner edit
    Int_t getNumV0s() const { return fNumV0s; }
    Ali_TRD_ST_V0* getV0(Int_t i) const { return i < fNumV0s ? (Ali_TRD_ST_V0*)((*fV0s)[i]) : NULL; }

    Ali_MC_particle* createMCparticle() {
        if (fNumMCparticles == fMCparticles->GetSize()) {
            fMCparticles->Expand(fNumMCparticles + 10);
        }
        if (fNumMCparticles >= 650000) {
            Fatal("Ali_TRD_ST_Event::createMCparticle()", "ERROR: Too many MC particles (>650000)!");
            exit(2);
        }
        new ((*fMCparticles)[fNumMCparticles++]) Ali_MC_particle;
        return (Ali_MC_particle*)((*fMCparticles)[fNumMCparticles - 1]);
    }

    void clearMCparticleList() {
        fNumMCparticles = 0;
        fMCparticles->Clear();
    }

    Ali_TRD_ST_TOF_hit* createTOFhit()  // online tracklet
    {
        if (fNumTOFhits == fTOFhits->GetSize()) {
            fTOFhits->Expand(fNumTOFhits + 10);
        }
        if (fNumTOFhits >= 650000) {
            Fatal("Ali_TRD_ST_Event::createTOFhit()", "ERROR: Too many TOF hits (>650000)!");
            exit(2);
        }
        new ((*fTOFhits)[fNumTOFhits++]) Ali_TRD_ST_TOF_hit;
        return (Ali_TRD_ST_TOF_hit*)((*fTOFhits)[fNumTOFhits - 1]);
    }

    void clearTOFhitList() {
        fNumTOFhits = 0;
        fTOFhits->Clear();
    }

    Ali_TRD_ST_TPC_Track* createTrack() {
        if (fNumTracks == fTracks->GetSize()) {
            fTracks->Expand(fNumTracks + 10);
        }
        if (fNumTracks >= 65000) {
            Fatal("Ali_TRD_ST_Event::createTrack()", "ERROR: Too many tracks (>65000)!");
            exit(2);
        }

        new ((*fTracks)[fNumTracks++]) Ali_TRD_ST_TPC_Track;
        return (Ali_TRD_ST_TPC_Track*)((*fTracks)[fNumTracks - 1]);
    }

    void clearTrackList() {
        fNumTracks = 0;
        fTracks->Clear();
    }

    Ali_TRD_ST_Tracklets* createTracklet()  // online tracklet
    {
        if (fNumTracklets == fTracklets->GetSize()) {
            fTracklets->Expand(fNumTracklets + 10);
        }
        if (fNumTracklets >= 650000) {
            Fatal("Ali_TRD_ST_Event::createTracklet()", "ERROR: Too many tracklets (>650000)!");
            exit(2);
        }

        new ((*fTracklets)[fNumTracklets++]) Ali_TRD_ST_Tracklets;
        return (Ali_TRD_ST_Tracklets*)((*fTracklets)[fNumTracklets - 1]);
    }

    void clearTrackletList() {
        fNumTracklets = 0;
        fTracklets->Clear();
    }

    // hoppner edit
    Ali_TRD_ST_V0* createV0() {
        if (fNumV0s == fV0s->GetSize()) {
            fV0s->Expand(fNumV0s + 10);
        }
        if (fNumV0s >= 65000) {
            Fatal("Ali_TRD_ST_Event::createV0()", "ERROR: Too many V0s (>65000)!");
            exit(2);
        }

        new ((*fV0s)[fNumV0s++]) Ali_TRD_ST_V0;
        return (Ali_TRD_ST_V0*)((*fV0s)[fNumV0s - 1]);
    }

    void clearV0List() {
        fNumV0s = 0;
        fV0s->Clear();
    }

    ClassDef(Ali_TRD_ST_Event, 1);  // A simple event compiled of tracks
};

#endif
