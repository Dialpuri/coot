#pragma once

#include <string>
#include <map>
#include <vector>
#include <gemmi/model.hpp>

namespace coot {

class protein_geometry;

struct atom_overlap_t {
    int ligand_atom_index;
    double r_1, r_2;
    gemmi::Atom *atom_1;
    gemmi::Atom *atom_2;
    double overlap_volume;
    bool is_h_bond;
    atom_overlap_t() : ligand_atom_index(-1), r_1(-1), r_2(-1), atom_1(nullptr), atom_2(nullptr), overlap_volume(-1), is_h_bond(false) {}
};

struct atom_overlaps_dots_container_t {
    std::map<std::string, size_t> dots;
    atom_overlaps_dots_container_t() = default;
};

class atom_overlaps_container_t {
public:
    atom_overlaps_container_t(gemmi::Structure* mol_in,
                              const protein_geometry* geom_p_in,
                              bool ignore_water_contacts_flag,
                              double clash_spike_length_in = 0.5,
                              double probe_radius_in = 0.25);
    
    atom_overlaps_container_t() = default;
    
    atom_overlaps_dots_container_t
    all_atom_contact_dots(double dot_density_in, bool make_vdw_surface);
    
    void make_all_atom_overlaps();
    
    size_t overlaps_size() const { return overlaps.size(); }
    
private:
    gemmi::Structure* mol = nullptr;
    const protein_geometry* geom_p = nullptr;
    bool have_dictionary = false;
    bool molecule_has_hydrogens = false;
    double probe_radius = 0.25;
    bool ignore_water_contacts_flag = false;
    
    std::vector<atom_overlap_t> overlaps;
    
    atom_overlaps_dots_container_t
    all_atom_contact_dots_internal_single_thread(double dot_density_in,
                                                  gemmi::Model* model,
                                                  unsigned int selection_handle1,
                                                  unsigned int selection_handle2,
                                                  double min_dist, double max_dist,
                                                  bool make_vdw_surface);
    
    atom_overlaps_dots_container_t
    all_atom_contact_dots_internal_multi_thread(double dot_density_in,
                                                 gemmi::Model* model,
                                                 unsigned int selection_handle1,
                                                 unsigned int selection_handle2,
                                                 double min_dist, double max_dist,
                                                 bool make_vdw_surface);
};

unsigned int get_max_number_of_threads();

} // namespace coot