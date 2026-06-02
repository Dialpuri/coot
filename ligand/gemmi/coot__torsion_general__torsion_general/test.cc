#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

static coot::atom_spec_t make_atom_spec_from_gemmi_atom(
    const gemmi::Chain& chain,
    const gemmi::Residue& res,
    const gemmi::Atom& atom) {
    coot::atom_spec_t spec;
    spec.chain_id   = chain.name;
    spec.res_no     = res.seqid.num.value;
    spec.ins_code   = (res.seqid.icode == ' ') ? "" : std::string(1, res.seqid.icode);
    spec.atom_name  = atom.name;
    spec.alt_conf   = (atom.altloc == '\0' || atom.altloc == ' ') ? "" : std::string(1, atom.altloc);
    spec.model_number = 1;
    return spec;
}

TEST(OracleTest, torsion_general) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    gemmi::Model& model = st.models[0];

    // Helper to find a residue by chain and sequence number
    auto find_residue = [&](const std::string& chain_id, int resno) -> gemmi::Residue* {
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name == chain_id) {
                for (gemmi::Residue& res : chain.residues) {
                    if (res.seqid.num.value == resno) return &res;
                }
            }
        }
        return nullptr;
    };

    // Helper to find an atom by name in a residue (gemmi uses unpadded names)
    auto find_atom = [&](gemmi::Residue& res, const std::string& name) -> const gemmi::Atom* {
        for (const gemmi::Atom& atom : res.atoms) {
            if (atom.name == name) return &atom;
        }
        return nullptr;
    };

    // Helper to get the parent chain for a residue
    auto find_chain_for_residue = [&](gemmi::Residue& res) -> gemmi::Chain* {
        for (gemmi::Chain& chain : model.chains) {
            for (gemmi::Residue& r : chain.residues) {
                if (&r == &res) return &chain;
            }
        }
        return nullptr;
    };

    // --- Case 1: valid torsion on MET residue 1 (atoms N, CA, C, CB) ---
    {
        gemmi::Residue* res = find_residue("A", 1);
        EXPECT_NE(res, nullptr);

        // gemmi atom.name is unpadded: "N", "CA", "C", "CB"
        const gemmi::Atom* a0 = find_atom(*res, "N");
        const gemmi::Atom* a1 = find_atom(*res, "CA");
        const gemmi::Atom* a2 = find_atom(*res, "C");
        const gemmi::Atom* a3 = find_atom(*res, "CB");

        EXPECT_NE(a0, nullptr);
        EXPECT_NE(a1, nullptr);
        EXPECT_NE(a2, nullptr);
        EXPECT_NE(a3, nullptr);

        gemmi::Chain* chain = find_chain_for_residue(*res);

        std::vector<coot::atom_spec_t> torsion_atoms;
        torsion_atoms.push_back(make_atom_spec_from_gemmi_atom(*chain, *res, *a0));
        torsion_atoms.push_back(make_atom_spec_from_gemmi_atom(*chain, *res, *a1));
        torsion_atoms.push_back(make_atom_spec_from_gemmi_atom(*chain, *res, *a2));
        torsion_atoms.push_back(make_atom_spec_from_gemmi_atom(*chain, *res, *a3));

        coot::torsion_general::torsion_general_result tg =
            coot::torsion_general::torsion_general_gemmi(res, torsion_atoms);

        EXPECT_TRUE(tg.setup_correctly);
        EXPECT_EQ(tg.clicked_atom_indices.size(), 4u);
    }

    // --- Case 2: invalid atom name (ZZZ does not exist) ---
    {
        gemmi::Residue* res = find_residue("A", 1);
        EXPECT_NE(res, nullptr);

        const gemmi::Atom* a0 = find_atom(*res, "N");
        const gemmi::Atom* a1 = find_atom(*res, "CA");
        const gemmi::Atom* a2 = find_atom(*res, "C");

        EXPECT_NE(a0, nullptr);
        EXPECT_NE(a1, nullptr);
        EXPECT_NE(a2, nullptr);

        gemmi::Chain* chain = find_chain_for_residue(*res);

        coot::atom_spec_t bad_spec;
        bad_spec.chain_id = "A";
        bad_spec.res_no = 1;
        bad_spec.ins_code = "";
        bad_spec.atom_name = "ZZZ";
        bad_spec.alt_conf = "";
        bad_spec.model_number = 1;

        std::vector<coot::atom_spec_t> torsion_atoms;
        torsion_atoms.push_back(make_atom_spec_from_gemmi_atom(*chain, *res, *a0));
        torsion_atoms.push_back(make_atom_spec_from_gemmi_atom(*chain, *res, *a1));
        torsion_atoms.push_back(make_atom_spec_from_gemmi_atom(*chain, *res, *a2));
        torsion_atoms.push_back(bad_spec);

        coot::torsion_general::torsion_general_result tg =
            coot::torsion_general::torsion_general_gemmi(res, torsion_atoms);

        EXPECT_FALSE(tg.setup_correctly);
    }

    // --- Case 3: residue with few atoms — GLY 11 (no CB) ---
    {
        gemmi::Residue* res = find_residue("A", 11);
        EXPECT_NE(res, nullptr);

        const gemmi::Atom* a0 = find_atom(*res, "N");
        const gemmi::Atom* a1 = find_atom(*res, "CA");
        const gemmi::Atom* a2 = find_atom(*res, "C");

        EXPECT_NE(a0, nullptr);
        EXPECT_NE(a1, nullptr);
        EXPECT_NE(a2, nullptr);

        gemmi::Chain* chain = find_chain_for_residue(*res);

        coot::atom_spec_t bad_spec;
        bad_spec.chain_id = "A";
        bad_spec.res_no = res->seqid.num.value;
        bad_spec.ins_code = (res->seqid.icode == ' ') ? "" : std::string(1, res->seqid.icode);
        bad_spec.atom_name = "CB";
        bad_spec.alt_conf = "";
        bad_spec.model_number = 1;

        std::vector<coot::atom_spec_t> torsion_atoms;
        torsion_atoms.push_back(make_atom_spec_from_gemmi_atom(*chain, *res, *a0));
        torsion_atoms.push_back(make_atom_spec_from_gemmi_atom(*chain, *res, *a1));
        torsion_atoms.push_back(make_atom_spec_from_gemmi_atom(*chain, *res, *a2));
        torsion_atoms.push_back(bad_spec);

        coot::torsion_general::torsion_general_result tg =
            coot::torsion_general::torsion_general_gemmi(res, torsion_atoms);

        EXPECT_FALSE(tg.setup_correctly);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}