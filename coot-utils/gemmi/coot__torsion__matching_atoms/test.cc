#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

std::string trim_atom_name(const std::string& s) {
    std::string t(s);
    t.erase(0, t.find_first_not_of(" "));
    t.erase(t.find_last_not_of(" ") + 1, std::string::npos);
    return t;
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

TEST(OracleTest, torsion_matching_atoms) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    gemmi::Model& model = st.models[0];

    // Find chain A, GLU residue 2
    gemmi::Chain* glu_chain = nullptr;
    gemmi::Residue* glu_res = nullptr;
    for (gemmi::Chain& ch : model.chains) {
        if (ch.name == "A") {
            for (gemmi::Residue& res : ch.residues) {
                if (res.name == "GLU" && res.seqid.num.value == 2) {
                    glu_chain = &ch;
                    glu_res = &res;
                    break;
                }
            }
        }
        if (glu_res) break;
    }
    ASSERT_NE(glu_res, nullptr);

    // Find chain A, GLY residue 11
    gemmi::Chain* gly_chain = nullptr;
    gemmi::Residue* gly_res = nullptr;
    for (gemmi::Chain& ch : model.chains) {
        if (ch.name == "A") {
            for (gemmi::Residue& res : ch.residues) {
                if (res.name == "GLY" && res.seqid.num.value == 11) {
                    gly_chain = &ch;
                    gly_res = &res;
                    break;
                }
            }
        }
        if (gly_res) break;
    }
    ASSERT_NE(gly_res, nullptr);

    // Case 1: GLU residue 2 — override res_no with wrong value to break matches
    {
        gemmi::Atom* n_atom = nullptr, *ca_atom = nullptr, *cb_atom = nullptr, *cg_atom = nullptr;
        for (gemmi::Atom& atom : glu_res->atoms) {
            std::string an = trim_atom_name(atom.name);
            if (an == "N")  n_atom  = &atom;
            if (an == "CA") ca_atom = &atom;
            if (an == "CB") cb_atom = &atom;
            if (an == "CG") cg_atom = &atom;
        }

        if (n_atom && ca_atom && cb_atom && cg_atom) {
            int wrong_res_no = glu_res->seqid.num.value + 1;

            coot::atom_spec_t a1; a1.chain_id = "A"; a1.res_no = wrong_res_no; a1.ins_code = " "; a1.atom_name = n_atom->name;  a1.alt_conf = std::string(1, n_atom->altloc);
            coot::atom_spec_t a2; a2.chain_id = "A"; a2.res_no = wrong_res_no; a2.ins_code = " "; a2.atom_name = ca_atom->name; a2.alt_conf = std::string(1, ca_atom->altloc);
            coot::atom_spec_t a3; a3.chain_id = "A"; a3.res_no = wrong_res_no; a3.ins_code = " "; a3.atom_name = cb_atom->name; a3.alt_conf = std::string(1, cb_atom->altloc);
            coot::atom_spec_t a4; a4.chain_id = "A"; a4.res_no = wrong_res_no; a4.ins_code = " "; a4.atom_name = cg_atom->name; a4.alt_conf = std::string(1, cg_atom->altloc);

            std::vector<gemmi::Atom*> atoms = coot::torsion::matching_atoms_gemmi(
                glu_chain, glu_res,
                std::make_pair(0, a1),
                std::make_pair(0, a2),
                std::make_pair(0, a3),
                std::make_pair(0, a4));

            EXPECT_EQ(atoms.size(), 0u);
        }
    }

    // Case 2: GLY residue 11 — only N, CA, C, O (no CB, CG)
    {
        gemmi::Atom* n_atom = nullptr, *ca_atom = nullptr;
        for (gemmi::Atom& atom : gly_res->atoms) {
            std::string an = trim_atom_name(atom.name);
            if (an == "N")  n_atom  = &atom;
            if (an == "CA") ca_atom = &atom;
        }

        if (n_atom && ca_atom) {
            int wrong_res_no = gly_res->seqid.num.value + 1;

            coot::atom_spec_t a1; a1.chain_id = "A"; a1.res_no = wrong_res_no; a1.ins_code = " "; a1.atom_name = n_atom->name;  a1.alt_conf = std::string(1, n_atom->altloc);
            coot::atom_spec_t a2; a2.chain_id = "A"; a2.res_no = wrong_res_no; a2.ins_code = " "; a2.atom_name = ca_atom->name; a2.alt_conf = std::string(1, ca_atom->altloc);
            coot::atom_spec_t a3; a3.chain_id = "A"; a3.res_no = wrong_res_no; a3.ins_code = " "; a3.atom_name = " CB "; a3.alt_conf = " ";
            coot::atom_spec_t a4; a4.chain_id = "A"; a4.res_no = wrong_res_no; a4.ins_code = " "; a4.atom_name = " CG "; a4.alt_conf = " ";

            std::vector<gemmi::Atom*> atoms = coot::torsion::matching_atoms_gemmi(
                gly_chain, gly_res,
                std::make_pair(0, a1),
                std::make_pair(0, a2),
                std::make_pair(0, a3),
                std::make_pair(0, a4));

            EXPECT_EQ(atoms.size(), 0u);
        }
    }

    // Case 3 (complementary): GLU residue 2 — correct res_no, should match all 4
    {
        gemmi::Atom* n_atom = nullptr, *ca_atom = nullptr, *cb_atom = nullptr, *cg_atom = nullptr;
        for (gemmi::Atom& atom : glu_res->atoms) {
            std::string an = trim_atom_name(atom.name);
            if (an == "N")  n_atom  = &atom;
            if (an == "CA") ca_atom = &atom;
            if (an == "CB") cb_atom = &atom;
            if (an == "CG") cg_atom = &atom;
        }

        if (n_atom && ca_atom && cb_atom && cg_atom) {
            int correct_res_no = glu_res->seqid.num.value;

            coot::atom_spec_t a1; a1.chain_id = "A"; a1.res_no = correct_res_no; a1.ins_code = " "; a1.atom_name = n_atom->name;  a1.alt_conf = std::string(1, n_atom->altloc);
            coot::atom_spec_t a2; a2.chain_id = "A"; a2.res_no = correct_res_no; a2.ins_code = " "; a2.atom_name = ca_atom->name; a2.alt_conf = std::string(1, ca_atom->altloc);
            coot::atom_spec_t a3; a3.chain_id = "A"; a3.res_no = correct_res_no; a3.ins_code = " "; a3.atom_name = cb_atom->name; a3.alt_conf = std::string(1, cb_atom->altloc);
            coot::atom_spec_t a4; a4.chain_id = "A"; a4.res_no = correct_res_no; a4.ins_code = " "; a4.atom_name = cg_atom->name; a4.alt_conf = std::string(1, cg_atom->altloc);

            std::vector<gemmi::Atom*> atoms = coot::torsion::matching_atoms_gemmi(
                glu_chain, glu_res,
                std::make_pair(0, a1),
                std::make_pair(0, a2),
                std::make_pair(0, a3),
                std::make_pair(0, a4));

            // All four atoms should match
            EXPECT_EQ(atoms.size(), 4u);
            if (atoms.size() == 4) {
                EXPECT_EQ(trim_atom_name(atoms[0]->name), "N");
                EXPECT_EQ(trim_atom_name(atoms[1]->name), "CA");
                EXPECT_EQ(trim_atom_name(atoms[2]->name), "CB");
                EXPECT_EQ(trim_atom_name(atoms[3]->name), "CG");
            }
        }
    }
}