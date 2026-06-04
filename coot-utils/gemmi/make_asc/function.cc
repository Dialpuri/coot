#include "function.hh"
#include <unordered_map>

atom_selection_container_t make_asc_gemmi(const gemmi::Structure& st, bool transfer_atom_index_flag) {
    atom_selection_container_t asc;
    asc.mol = &st;

    // Select all atoms from model 0 — mirrors MMDB SelectAtoms with all wildcards
    asc.SelectionHandle = 0;
    const gemmi::Model& model = st.models[0];
    int idx = 0;
    for (const gemmi::Chain& chain : model.chains) {
        for (const gemmi::Residue& res : chain.residues) {
            for (const gemmi::Atom& atom : res.atoms) {
                asc.atom_selection.push_back(idx++);
            }
        }
    }
    asc.n_selected_atoms = idx;

    // Simulate MMDB UDD handle registration (base 0x1000000 for UDR_ATOM)
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