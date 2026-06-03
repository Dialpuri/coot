#pragma once

#include <gemmi/model.hpp>
#include <gemmi/neighbor.hpp>
#include <gemmi/unitcell.hpp>
#include <gemmi/elem.hpp>
#include <string>
#include <vector>

namespace coot {
namespace util {

// Gemmi-compatible contact info for a central atom and its neighbors
struct contact_atoms_info_t {
    gemmi::CRA cra;                              // central (water) atom
    std::vector<gemmi::CRA> contact_list;        // contacting atoms

    gemmi::Atom const* central_atom() const { return cra.atom; }
    std::string central_atom_name() const {
        if (!cra.atom) return "";
        // MMDB-style padding: " X  " where X is the atom name (4 chars total)
        return std::string(" ") + cra.atom->name + std::string(3 - cra.atom->name.size(), ' ');
    }
    gemmi::Atom const* contact_atom(size_t i) const { return contact_list[i].atom; }
    size_t size() const { return contact_list.size(); }
};

class water_coordination_t {
public:
    water_coordination_t(const gemmi::Structure& structure, int model_index, double radius)
    {
        // Need a mutable copy for NeighborSearch
        model_ = std::make_shared<gemmi::Model>(structure.models[model_index]);
        init_internal(model_, structure.cell, radius, false);
    }

    std::vector<contact_atoms_info_t> get_contacts() const { return contacts_; }

private:
    void init_internal(
        const std::shared_ptr<gemmi::Model>& model,
        const gemmi::UnitCell& cell,
        double radius,
        bool only_heavy_atoms)
    {
        radius_      = radius;
        only_heavy_atoms_ = only_heavy_atoms;

        // ---- find water O atoms (HOH hetatm, element O) ----
        std::vector<gemmi::CRA> water_atoms;
        for (auto& chain : model->chains) {
            for (auto& res : chain.residues) {
                if (res.het_flag == 'H' && res.name == "HOH") {
                    for (auto& atom : res.atoms) {
                        if (atom.element.name()[0] == 'O' && atom.element.name()[1] == '\0') {
                            gemmi::CRA cra;
                            cra.chain   = &chain;
                            cra.residue = &res;
                            cra.atom    = &atom;
                            water_atoms.push_back(cra);
                        }
                    }
                }
            }
        }

        // ---- neighbour search over the whole model using crystal symmetry ----
        gemmi::NeighborSearch ns(*model, cell, radius);
        ns.populate(!only_heavy_atoms_);

        // ---- for each water O, find contacts and add them ----
        for (auto& water_cra : water_atoms) {
            auto neighbors = ns.find_atoms(water_cra.atom->pos, '\0', 0.0, radius);
            for (auto* mark : neighbors) {
                gemmi::CRA ncra = mark->to_cra(*model);
                if (ncra.atom && ncra.atom != water_cra.atom) {
                    add_contact(water_cra, ncra);
                }
            }
        }
    }

    void add_contact(gemmi::CRA central_cra, gemmi::CRA contact_cra)
    {
        // alt-loc compatibility (MMDB: empty "" == gemmi ' ')
        char alt1 = contact_cra.atom->altloc;
        char alt2 = central_cra.atom->altloc;
        if (!((alt1 == alt2) || (alt1 == ' ') || (alt2 == ' '))) return;

        // filter out H contactors (element symbol "H")
        if (contact_cra.atom->element.name() == "H") return;

        // append to existing entry or create new one
        for (auto& info : contacts_) {
            if (info.cra.atom == central_cra.atom) {
                info.contact_list.push_back(contact_cra);
                return;
            }
        }
        contact_atoms_info_t info;
        info.cra = central_cra;
        info.contact_list.push_back(contact_cra);
        contacts_.push_back(std::move(info));
    }

    // ---- members ----
    std::shared_ptr<gemmi::Model> model_;
    double radius_;
    bool only_heavy_atoms_;
    std::vector<contact_atoms_info_t> contacts_;
};

} // namespace util
} // namespace coot