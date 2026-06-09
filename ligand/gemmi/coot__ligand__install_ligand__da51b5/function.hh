#pragma once

#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>
#include <vector>
#include <string>
#include <cstdint>

// ── minimal coot type mirrors (self-contained; avoids including mini-mol.hh
//    which transitively pulls in <mmdb2/mmdb_manager.h>) ────────────────

namespace coot {
namespace minimol {

class atom {
public:
    std::string name;
    std::string element;
    float x, y, z;
    std::string altLoc;
    float occupancy;
    float temperature_factor;

    atom(const std::string &nm = "",
         const std::string &el = "",
         float xx = 0.f, float yy = 0.f, float zz = 0.f,
         const std::string &al = "",
         float o  = 1.0f, float b  = 0.0f)
        : name(nm), element(el), x(xx), y(yy), z(zz),
          altLoc(al), occupancy(o), temperature_factor(b) {}

    bool is_hydrogen_p() const {
        if (element.empty()) return false;
        char c = element[0];
        return (c == 'H' || c == 'h');
    }
};

class residue {
public:
    int seqnum;
    std::string ins_code;
    std::string name;
    std::vector<atom> atoms;

    residue() : seqnum(0), ins_code(""), name("") {}
    residue(int i, const std::string &resname)
        : seqnum(i), ins_code(""), name(resname) {}

    void addatom(const std::string &atom_name, const std::string &element,
                 float x, float y, float z, const std::string &altloc,
                 float occupancy, float bf) {
        atoms.emplace_back(atom_name, element, x, y, z, altloc, occupancy, bf);
    }

    unsigned int n_atoms() const { return atoms.size(); }
};

class fragment {
public:
    std::string fragment_id;
    std::vector<residue> residues;

    fragment() : fragment_id("") {}
    explicit fragment(const std::string &frag_id_in) : fragment_id(frag_id_in) {}
};

class molecule {
public:
    std::vector<fragment> fragments;
    std::vector<std::string> labels;
    int n_atom;
    bool ok;
    std::string filename;

    molecule() : n_atom(0), ok(false) {}

    molecule &operator=(const molecule &other) {
        if (this != &other) {
            fragments = other.fragments;
            labels    = other.labels;
            n_atom    = other.n_atom;
            ok        = other.ok;
            filename  = other.filename;
        }
        return *this;
    }

    void init() {
        fragments.clear();
        labels.clear();
        n_atom = 0;
        ok = false;
        filename.clear();
    }
};

} // namespace minimol
} // namespace coot

// ── coot::ligand (minimal mirror — only what install_ligand needs) ─────

namespace coot {

class ligand {
public:
    std::vector<coot::minimol::molecule> initial_ligand;

    size_t n_initial_ligands() const { return initial_ligand.size(); }

    void install_ligand_gemmi(const coot::minimol::molecule &mol) {
        int ilig = static_cast<int>(initial_ligand.size());
        initial_ligand.resize(static_cast<size_t>(ilig) + 1);
        initial_ligand[ilig] = mol;
        make_ligand_properties(ilig);
    }

private:
    void make_ligand_properties(int ilig) {
        coot::minimol::molecule &mol = initial_ligand[ilig];
        mol.labels.clear();
        mol.n_atom = 0;
        for (const auto &frag : mol.fragments) {
            for (const auto &res : frag.residues) {
                for (const auto &at : res.atoms) {
                    mol.labels.push_back(at.name);
                    mol.n_atom++;
                }
            }
        }
        mol.ok = (mol.n_atom > 0);
    }
};

} // namespace coot
