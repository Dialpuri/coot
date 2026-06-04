#include "function.hh"
#include <iostream>

namespace coot {

// Inline make_asc_gemmi — copies logic from the generated port .cc
static atom_selection_container_t my_make_asc_gemmi(const gemmi::Structure& st, bool transfer_atom_index_flag) {
    atom_selection_container_t asc;

    asc.mol = const_cast<gemmi::Structure*>(&st);

    // Select all atoms from model 0
    asc.SelectionHandle = 0;
    int atom_index = 0;
    if (!st.models.empty()) {
        const gemmi::Model& model = st.models[0];
        for (const gemmi::Chain& chain : model.chains) {
            for (const gemmi::Residue& res : chain.residues) {
                for (const gemmi::Atom& atom : res.atoms) {
                    asc.atom_selection.push_back(atom_index++);
                }
            }
        }
    }
    asc.n_selected_atoms = atom_index;

    static const int UDR_ATOM_BASE = 0x1000000;
    static std::unordered_map<std::string, int> udd_handles;
    if (udd_handles.find("atom index") == udd_handles.end()) {
        udd_handles["atom index"] = UDR_ATOM_BASE + static_cast<int>(udd_handles.size()) + 1;
    }
    asc.UDDAtomIndexHandle = udd_handles["atom index"];

    asc.read_error_message = "No error";
    asc.read_success = 1;

    if (transfer_atom_index_flag) {
        asc.UDDOldAtomIndexHandle = udd_handles["atom index"];
    } else {
        asc.UDDOldAtomIndexHandle = -1;
    }

    return asc;
}

atom_selection_container_t mol_to_asc_rdkit_gemmi(const std::string& file_name) {
    atom_selection_container_t asc{};  // zero-initialize all members
    asc.read_success = 0;
    asc.n_selected_atoms = 0;
    asc.mol = nullptr;

    try {
        RDKit::ROMol* m = RDKit::Mol2FileToMol(file_name);

        std::string res_name = "UNL";
        try {
            m->getProp("_Name", res_name);
        } catch (...) {
            std::cout << "mol_to_asc_rdkit_gemmi() no rdkit molecule name for " << m << std::endl;
        }

        if (m) {
            gemmi::Residue* res = coot::make_residue_gemmi(*m, 0, res_name);
            if (res) {
                gemmi::Chain temp_chain;
                temp_chain.name = "A";
                gemmi::CRA cra{&temp_chain, res, nullptr};
                gemmi::Structure* mol = coot::util::create_mmdbmanager_from_residue_gemmi(cra);
                if (mol) {
                    asc = my_make_asc_gemmi(*mol, false);
                }
                delete res;
            }
            delete m;
        } else {
            std::cout << "Null rdkit mol ptr m" << std::endl;
        }
    }
    catch (const RDKit::FileParseException& rte) {
        std::cout << "WARNING:: " << rte.what() << std::endl;
    }
    catch (const std::runtime_error& rte) {
        std::cout << "WARNING runtime_error in mol_to_asc_rdkit_gemmi() " << rte.what() << std::endl;
    }
    catch (const std::exception& e) {
        std::cout << "WARNING:: mol_to_asc_rdkit_gemmi: exception: " << e.what() << std::endl;
    }
    return asc;
}

} // namespace coot