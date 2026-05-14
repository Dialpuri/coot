#pragma once
#include <vector>
#include <string>
#include "clipper/core/coords.h"
#include "coot/api/coot-molecule.hh"

namespace pli {

struct pi_stacking_instance_t {
    enum stacking_t { NO_STACKING, PI_PI_STACKING, PI_CATION_STACKING };
    pi_stacking_instance_t() {}
    pi_stacking_instance_t(void *res, stacking_t type, const std::vector<std::string>& rings);
    pi_stacking_instance_t(void *res, const std::string& cation_name);
    
    stacking_t stacking_type;
    std::vector<std::string> aromatic_ring_atoms;
    std::string cation_atom_name;
    float overlap_score;
    void *residue_ptr;
};

class pi_stacking_container_t {
public:
    pi_stacking_container_t(const coot::dictionary_residue_restraints_t &monomer_restraints,
                           const std::vector<void*> &residues,
                           void *res_ref,
                           const std::vector<std::vector<std::string> > &aromatic_ring_list);
    void init(const coot::dictionary_residue_restraints_t &monomer_restraints,
             const std::vector<void*> &residues,
             void *res_ref,
             const std::vector<std::vector<std::string> > &aromatic_ring_list);
    
    std::vector<pi_stacking_instance_t> stackings;
    
private:
    std::pair<clipper::Coord_orth, clipper::Coord_orth> 
    get_ring_pi_centre_points(const std::vector<std::string>& ring_atoms, void *res_ref);
    std::pair<float, pi_stacking_instance_t::stacking_t> 
    get_pi_overlap_to_ligand_ring(void *res, const clipper::Coord_orth& pi_pt);
    std::vector<std::pair<std::string, clipper::Coord_orth> > 
    get_ligand_cations(void *res_ref, const coot::dictionary_residue_restraints_t &restraints);
    float get_pi_overlap_to_ligand_cation(void *res, const clipper::Coord_orth& cation_pt);
};

} // namespace pli