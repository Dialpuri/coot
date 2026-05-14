#pragma once
#include <vector>
#include <string>
#include <gemmi/model.hpp>
#include <clipper/core/coords.h>
#include "geometry/protein-geometry.hh"

namespace pli {

class pi_stacking_instance_t {
public:
    enum stacking_t {
        NO_STACKING,
        PI_PI_STACKING,
        PI_CATION_STACKING,
        CATION_PI_STACKING,
    };
    
    gemmi::Residue* res;
    stacking_t type;
    std::vector<std::string> ligand_ring_atom_names;
    float overlap_score;
    std::string ligand_cationic_atom_name;

    pi_stacking_instance_t(gemmi::Residue* res_in, stacking_t type_in,
                          const std::vector<std::string>& ring_atoms)
        : ligand_ring_atom_names(ring_atoms) {
        res = res_in;
        type = type_in;
        overlap_score = 0;
    }

    pi_stacking_instance_t(gemmi::Residue* res_in, const std::string& ligand_atom_name_in)
        : ligand_cationic_atom_name(ligand_atom_name_in) {
        type = CATION_PI_STACKING;
        res = res_in;
        overlap_score = 0;
    }
};

std::ostream& operator<<(std::ostream& s, const pi_stacking_instance_t& spec);

class pi_stacking_container_t {
private:
    std::pair<float, pi_stacking_instance_t::stacking_t>
    get_pi_overlap_to_ligand_ring(gemmi::Residue* res, const clipper::Coord_orth& pt) const;

    float get_pi_overlap_to_ligand_cation(gemmi::Residue* res, const clipper::Coord_orth& pt) const;

    std::pair<clipper::Coord_orth, clipper::Coord_orth>
    get_ring_pi_centre_points(const std::vector<std::string>& ring_atom_names,
                             gemmi::Residue* res_ref) const;

    std::pair<clipper::Coord_orth, clipper::Coord_orth>
    ring_centre_and_normal(const std::vector<clipper::Coord_orth>& pts) const;

    std::vector<std::vector<std::string>>
    ring_atom_names(const std::string& residue_name) const;

    float overlap_of_pi_spheres(const clipper::Coord_orth& pt1,
                               const clipper::Coord_orth& pt2,
                               double m1_pt_1, double m2_pt_1,
                               double m1_pt_2, double m2_pt_2) const;

    float overlap_of_cation_pi(const clipper::Coord_orth& ligand_pi_point,
                              const clipper::Coord_orth& cation_atom_point) const;

    std::vector<clipper::Coord_orth> get_cation_atom_positions(gemmi::Residue* res) const;
    std::vector<std::pair<std::string, clipper::Coord_orth>>
    get_ligand_cations(gemmi::Residue* res, const coot::dictionary_residue_restraints_t& monomer_restraints) const;

    std::vector<std::vector<std::string>>
    get_aromatic_ring_list(const coot::dictionary_residue_restraints_t& monomer_restraints) const;

    void init(const coot::dictionary_residue_restraints_t& monomer_restraints,
             const std::vector<gemmi::Residue*>& filtered_residues,
             gemmi::Residue* res_ref,
             const std::vector<std::vector<std::string>>& aromatic_ring_list_atom_names);

public:
    std::vector<pi_stacking_instance_t> stackings;
    
    pi_stacking_container_t(const coot::dictionary_residue_restraints_t& monomer_restraints,
                           const std::vector<gemmi::Residue*>& filtered_residues,
                           gemmi::Residue* res_ref);
};

} // namespace pli