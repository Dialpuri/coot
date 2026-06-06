#include "function.hh"

namespace coot {

bool is_angle_related_via_link_gemmi(
    gemmi::CRA cra1,
    gemmi::CRA cra2,
    gemmi::Model* model_p,
    const std::vector<gemmi::Connection>& connections,
    const std::vector<std::pair<std::string, std::string>>& bonds_for_at_1,
    const std::vector<std::pair<std::string, std::string>>& bonds_for_at_2
) {
    bool status = false;

    // Guard: null atoms
    if (!cra1.atom) return false;
    if (!cra2.atom) return false;
    if (!cra1.chain || !cra1.residue) return false;
    if (!cra2.chain || !cra2.residue) return false;

    // Guard: null model
    if (!model_p) return false;

    auto make_spec = [](gemmi::CRA cra, gemmi::Model* mp) -> atom_spec_t {
        atom_spec_t spec;
        spec.chain_id       = cra.chain->name;
        spec.res_no         = cra.residue->seqid.num.value;
        spec.ins_code       = std::string(1, cra.residue->seqid.icode);
        spec.atom_name      = pad_atom_name_4(cra.atom->name);
        spec.alt_conf       = std::string(1, cra.atom->altloc);
        spec.model_number   = mp ? mp->num : 0;
        return spec;
    };

    atom_spec_t spec_1 = make_spec(cra1, model_p);
    atom_spec_t spec_2 = make_spec(cra2, model_p);

    for (const gemmi::Connection& conn : connections) {
        // Only consider COVALENT (Link) connections
        if (conn.type != gemmi::Connection::Type::Covale)
            continue;

        std::pair<atom_spec_t, atom_spec_t> linked_atoms = link_atoms_gemmi(conn, model_p);

        // Check: spec_1 matches first linked atom
        if (spec_1 == linked_atoms.first) {
            std::string linked_atom_2_name = linked_atoms.second.atom_name;
            for (const auto& bond : bonds_for_at_2) {
                if (bond.first == linked_atom_2_name && bond.second == spec_2.atom_name) { status = true; break; }
                if (bond.second == linked_atom_2_name && bond.first == spec_2.atom_name) { status = true; break; }
            }
        }

        // Check: spec_1 matches second linked atom
        if (spec_1 == linked_atoms.second) {
            std::string linked_atom_1_name = linked_atoms.first.atom_name;
            for (const auto& bond : bonds_for_at_2) {
                if (bond.first == linked_atom_1_name && bond.second == spec_2.atom_name) { status = true; break; }
                if (bond.second == linked_atom_1_name && bond.first == spec_2.atom_name) { status = true; break; }
            }
        }

        // Reverse: spec_2 matches first linked atom
        if (spec_2 == linked_atoms.first) {
            std::string linked_atom_2_name = linked_atoms.second.atom_name;
            for (const auto& bond : bonds_for_at_1) {
                if (bond.first == linked_atom_2_name && bond.second == spec_1.atom_name) { status = true; break; }
                if (bond.second == linked_atom_2_name && bond.first == spec_1.atom_name) { status = true; break; }
            }
        }

        // Reverse: spec_2 matches second linked atom
        if (spec_2 == linked_atoms.second) {
            std::string linked_atom_1_name = linked_atoms.first.atom_name;
            for (const auto& bond : bonds_for_at_1) {
                if (bond.first == linked_atom_1_name && bond.second == spec_1.atom_name) { status = true; break; }
                if (bond.second == linked_atom_1_name && bond.first == spec_1.atom_name) { status = true; break; }
            }
        }

        if (status) break;
    }

    return status;
}

} // namespace coot