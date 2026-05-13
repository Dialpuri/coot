#pragma once

#include <gemmi/unitcell.hpp>
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>
#include <string>

namespace coot {

// Cartesian coordinate type (from coot-utils/coot-utils.hh)
class Cartesian {
private:
    double x_, y_, z_;
public:
    Cartesian() : x_(0), y_(0), z_(0) {}
    Cartesian(double x, double y, double z) : x_(x), y_(y), z_(z) {}
    double get_x() const { return x_; }
    double get_y() const { return y_; }
    double get_z() const { return z_; }
};

// Cell_Translation from mmdb-crystal.hh
class Cell_Translation {
public:
    int us, vs, ws;
    Cell_Translation() : us(0), vs(0), ws(0) {}
    Cell_Translation(int u, int v, int w) : us(u), vs(v), ws(w) {}
};

// atom_selection_container_t - simplified version for testing
struct atom_selection_container_t {
    gemmi::Structure* structure;
    gemmi::Model* model;
    int n_selected_atoms;
    // For testing, we just need access to the structure
};

// Helper function to get atom_selection_container_t from PDB path
inline atom_selection_container_t get_atom_selection(const std::string& pdb_path, 
                                                     bool, bool, bool) {
    atom_selection_container_t asc;
    asc.structure = new gemmi::Structure();
    *asc.structure = gemmi::read_pdb_file(pdb_path);
    asc.model = &asc.structure->models[0];
    asc.n_selected_atoms = 0;
    for (const auto& chain : asc.model->chains) {
        for (const auto& res : chain.residues) {
            asc.n_selected_atoms += res.atoms.size();
        }
    }
    return asc;
}

// molecule_extents_t class with the ported function
class molecule_extents_t {
public:
    molecule_extents_t(const atom_selection_container_t& asc, float) {}
    
    Cell_Translation coord_to_unit_cell_translations_gemmi(
        Cartesian point,
        const atom_selection_container_t& asc) const {
        
        // Get the unit cell from the structure
        const gemmi::UnitCell& cell = asc.structure->cell;
        
        // Convert Cartesian coordinates to fractional
        // gemmi uses Position (derived from Vec3) for coordinates
        gemmi::Position pos(point.get_x(), point.get_y(), point.get_z());
        gemmi::Fractional frac = cell.fractionalize(pos);
        
        // Round to nearest integer for unit cell translation
        int iu = static_cast<int>(std::rint(frac.x));
        int iv = static_cast<int>(std::rint(frac.y));
        int iw = static_cast<int>(std::rint(frac.z));
        
        return Cell_Translation(iu, iv, iw);
    }
};

} // namespace coot