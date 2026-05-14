#include "function.hh"
#include <iostream>
#include <stdexcept>

namespace pli {

pi_stacking_instance_t::pi_stacking_instance_t(void *res, stacking_t type, 
                                               const std::vector<std::string>& rings)
    : stacking_type(type), aromatic_ring_atoms(rings), overlap_score(0), residue_ptr(res) {}

pi_stacking_instance_t::pi_stacking_instance_t(void *res, const std::string& cation_name)
    : stacking_type(PI_CATION_STACKING), cation_atom_name(cation_name), overlap_score(0), residue_ptr(res) {}

pi_stacking_container_t::pi_stacking_container_t(const coot::dictionary_residue_restraints_t &monomer_restraints,
                                                const std::vector<void*> &residues,
                                                void *res_ref,
                                                const std::vector<std::vector<std::string> > &aromatic_ring_list) {
    init(monomer_restraints, residues, res_ref, aromatic_ring_list);
}

void
pi_stacking_container_t::init(const coot::dictionary_residue_restraints_t &monomer_restraints,
                             const std::vector<void*> &residues,
                             void *res_ref,
                             const std::vector<std::vector<std::string> > &aromatic_ring_list) {
    
    bool debug = false;
    float pi_pi_overlap_thresh = 0.1;
    float pi_cation_overlap_thresh = 30;
    
    for (unsigned int iring=0; iring<aromatic_ring_list.size(); iring++) {
        try {
            std::pair<clipper::Coord_orth, clipper::Coord_orth> ligand_ring_pi_pts = 
                get_ring_pi_centre_points(aromatic_ring_list[iring], res_ref);
            
            if (debug) {
                std::cout << "========= ligand ring ";
                for (unsigned int iat=0; iat<aromatic_ring_list[iring].size(); iat++)
                    std::cout << aromatic_ring_list[iring][iat] << "  ";
                std::cout << " ====== points" << std::endl;
            }
            
            for (unsigned int ires=0; ires<residues.size(); ires++) {
                
                if (debug) {
                    std::cout << "   ==== Environment residue " << ires << std::endl;
                }
                
                std::pair<float, pi_stacking_instance_t::stacking_t> pi_overlap_1 =
                    get_pi_overlap_to_ligand_ring(residues[ires], ligand_ring_pi_pts.first);
                std::pair<float, pi_stacking_instance_t::stacking_t> pi_overlap_2 =
                    get_pi_overlap_to_ligand_ring(residues[ires], ligand_ring_pi_pts.second);
                
                if (debug)
                    std::cout << "    protein cation:ligand ring: Overlaps:  score "
                              << pi_overlap_1.first << " type: " << pi_overlap_1.second << "  score: "
                              << pi_overlap_2.first << " type: " << pi_overlap_2.second << std::endl;
                
                float thresh = -1;
                if (pi_overlap_1.second == pi_stacking_instance_t::PI_PI_STACKING)
                    thresh = pi_pi_overlap_thresh;
                if (pi_overlap_1.second == pi_stacking_instance_t::PI_CATION_STACKING)
                    thresh = pi_cation_overlap_thresh;
                
                if (pi_overlap_1.first > thresh) {
                    pi_stacking_instance_t st(residues[ires],
                                             pi_overlap_1.second,
                                             aromatic_ring_list[iring]);
                    st.overlap_score = pi_overlap_1.first;
                    stackings.push_back(st);
                }
                if (pi_overlap_2.first > thresh) {
                    pli::pi_stacking_instance_t st(residues[ires],
                                                  pi_overlap_2.second,
                                                  aromatic_ring_list[iring]);
                    st.overlap_score = pi_overlap_2.first;
                    stackings.push_back(st);
                }
            }
        }
        catch (const std::runtime_error &rte) {
            std::cout << "WARNING:: " << rte.what() << std::endl;
        }
    }
    
    // Ligand cations (-> protein ring systems)
    std::vector<std::pair<std::string, clipper::Coord_orth> > cation_points =
        get_ligand_cations(res_ref, monomer_restraints);
    
    for (unsigned int icat=0; icat<cation_points.size(); icat++) {
        for (unsigned int ires=0; ires<residues.size(); ires++) {
            float score = get_pi_overlap_to_ligand_cation(residues[ires], cation_points[icat].second);
            
            if (score > pi_cation_overlap_thresh) { 
                pi_stacking_instance_t stacking(residues[ires], cation_points[icat].first);
                stacking.overlap_score = score;
                stackings.push_back(stacking);
            }
        }
    }
}

// Stub implementations for the helper methods
std::pair<clipper::Coord_orth, clipper::Coord_orth> 
pi_stacking_container_t::get_ring_pi_centre_points(const std::vector<std::string>& ring_atoms, void *res_ref) {
    return std::make_pair(clipper::Coord_orth(0,0,0), clipper::Coord_orth(0,0,0));
}

std::pair<float, pi_stacking_instance_t::stacking_t> 
pi_stacking_container_t::get_pi_overlap_to_ligand_ring(void *res, const clipper::Coord_orth& pi_pt) {
    return std::make_pair(0.0f, pi_stacking_instance_t::NO_STACKING);
}

std::vector<std::pair<std::string, clipper::Coord_orth> > 
pi_stacking_container_t::get_ligand_cations(void *res_ref, const coot::dictionary_residue_restraints_t &restraints) {
    return std::vector<std::pair<std::string, clipper::Coord_orth> >();
}

float pi_stacking_container_t::get_pi_overlap_to_ligand_cation(void *res, const clipper::Coord_orth& cation_pt) {
    return 0.0f;
}

} // namespace pli