#ifndef __ALI_MC_PARTICLE_H__
#define __ALI_MC_PARTICLE_H__

#include "TLorentzVector.h"
#include "TObject.h"
#include "TVector3.h"

class Ali_MC_particle : public TObject {
   private:
    TVector3 TV3_particle_vertex;
    TLorentzVector TLV_particle;
    Int_t PDGcode;
    Int_t index_particle;
    Int_t index_mother;
    Int_t N_daughters;
    Int_t arr_index_daughters[5];

   public:
    // constructor
    Ali_MC_particle()
        : TV3_particle_vertex(), TLV_particle(), PDGcode(0), index_particle(0), index_mother(0), N_daughters(0), arr_index_daughters() {}

    // destructor
    ~Ali_MC_particle() {}

    // setters
    void set_TV3_particle_vertex(TVector3 TV3_particle_vertex_in) { TV3_particle_vertex = TV3_particle_vertex_in; }
    void set_TLV_particle(TLorentzVector tlv) { TLV_particle = tlv; }
    void set_PDGcode(Int_t PDGcode_in) { PDGcode = PDGcode_in; }
    void set_index_mother(Int_t index_mother_in) { index_mother = index_mother_in; }
    void set_index_particle(Int_t index_particle_in) { index_particle = index_particle_in; }
    void set_N_daughters(Int_t N_daughters_in) { N_daughters = N_daughters_in; }
    void set_arr_index_daughters(Int_t index, Int_t index_daughter_in) { arr_index_daughters[index] = index_daughter_in; }

    // getters
    TVector3 get_TV3_particle_vertex() const { return TV3_particle_vertex; }
    TLorentzVector get_TLV_particle() const { return TLV_particle; }
    Int_t get_PDGcode() const { return PDGcode; }
    Int_t get_index_mother() const { return index_mother; }
    Int_t get_index_particle() const { return index_particle; }
    Int_t get_N_daughters() const { return N_daughters; }
    Int_t get_arr_index_daughters(Int_t index) const { return arr_index_daughters[index]; }

    ClassDef(Ali_MC_particle, 1);
};

#endif
