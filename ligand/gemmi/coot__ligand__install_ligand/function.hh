#pragma once

#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>
#include <vector>
#include <string>

// ── minimal coot type mirrors (self-contained; avoids including ligand.hh
//    which transitively pulls in <mmdb2/mmdb_manager.h>) ────────────────

namespace coot { namespace minimol {

class atom {
public:
    std::string name;
    std::string element;
    float x, y, z;
    std::string alt;
    float occ;
    float bf;

    atom(const std::string &nm = "",
         const std::string &el = "",
         float xx = 0.f, float yy = 0.f, float zz = 0.f,
         const std::string &al = "",
         float o  = 1.0f, float b  = 0.0f)
        : name(nm), element(el), x(xx), y(yy), z(zz),
          alt(al), occ(o), bf(b) {}
};

class molecule {
public:
    std::vector<atom> atoms;
    std::vector<std::string> labels;
    int n_atom;
    bool ok;
    std::string filename;

    molecule() : n_atom(0), ok(false) {}

    molecule &operator=(const molecule &other) {
        if (this != &other) {
            atoms    = other.atoms;
            labels   = other.labels;
            n_atom   = other.n_atom;
            ok       = other.ok;
            filename = other.filename;
        }
        return *this;
    }

    void init_from_gemmi(gemmi::Structure &st) {
        atoms.clear();
        labels.clear();
        filename.clear();
        n_atom = 0;
        ok     = false;

        for (gemmi::Model &model : st.models) {
            for (gemmi::Chain &chain : model.chains) {
                for (gemmi::Residue &res : chain.residues) {
                    for (gemmi::Atom &a : res.atoms) {
                        atoms.emplace_back(
                            a.name,
                            a.element.name(),
                            static_cast<float>(a.pos.x),
                            static_cast<float>(a.pos.y),
                            static_cast<float>(a.pos.z),
                            std::string(1, a.altloc),
                            static_cast<float>(a.occ),
                            static_cast<float>(a.b_iso)
                        );
                        n_atom++;
                    }
                }
            }
        }
        if (n_atom > 0) ok = true;
    }
};

}} // namespace coot::minimol

// ── coot::ligand (minimal mirror — only what install_ligand needs) ─────

namespace coot {

class ligand {
public:
    std::vector<coot::minimol::molecule> initial_ligand;

    void install_ligand_gemmi(gemmi::Structure &st) {
        int ilig = static_cast<int>(initial_ligand.size());
        initial_ligand.resize(static_cast<size_t>(ilig) + 1);
        initial_ligand[ilig].init_from_gemmi(st);
        make_ligand_properties(ilig);
    }

private:
    void make_ligand_properties(int ilig) {
        coot::minimol::molecule &mol = initial_ligand[ilig];
        mol.n_atom = static_cast<int>(mol.atoms.size());
        mol.labels.clear();
        mol.labels.reserve(static_cast<size_t>(mol.n_atom));
        for (int i = 0; i < mol.n_atom; i++) {
            mol.labels.push_back(mol.atoms[i].name);
        }
    }
};

} // namespace coot