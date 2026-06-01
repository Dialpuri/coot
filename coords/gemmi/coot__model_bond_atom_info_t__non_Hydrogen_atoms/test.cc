#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Pad an atom name to 4 characters the same way MMDB does (left-pad with space, right-pad to 4)
static std::string pad_atom_name(const std::string& name) {
    std::string padded = " " + name;
    while (padded.size() < 4) padded += ' ';
    return padded;
}

// Helper: find an atom by chain name, residue seqid, and atom name in a gemmi Structure
static const gemmi::Atom* find_atom_in_residue(
    const gemmi::Structure& st,
    const std::string& chain_id,
    int res_seq,
    const std::string& atom_name)
{
    for (const gemmi::Model& model : st.models) {
        for (const gemmi::Chain& chain : model.chains) {
            if (chain.name != chain_id) continue;
            for (const gemmi::Residue& res : chain.residues) {
                if (res.seqid.num.value != res_seq) continue;
                for (const gemmi::Atom& atom : res.atoms) {
                    std::string aname = atom.name;
                    while (!aname.empty() && aname.back() == ' ') aname.pop_back();
                    if (aname == atom_name) {
                        return &atom;
                    }
                }
            }
        }
    }
    return nullptr;
}

TEST(OracleTest, model_bond_atom_info_t) {
    gemmi::Structure st = gemmi::read_pdb_file(
        "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

    // Case 1: add atoms from residue A/1 (MET) — some H, some non-H
    {
        coot::model_bond_atom_info_t_gemmi bond_info;
        int atoms_added = 0;

        std::string atom_names[] = {"N","CA","C","O","CB","CG","CE","S",
                                     "H","HA","HB2","HB3","HG2","HG3",
                                     "HE1","HE2","HE3","H1","H2"};
        int num_atom_names = sizeof(atom_names)/sizeof(atom_names[0]);

        for (int i = 0; i < num_atom_names; i++) {
            const gemmi::Atom* at = find_atom_in_residue(st, "A", 1, atom_names[i]);
            if (at) {
                bond_info.add_bond_atom(*at);
                atoms_added++;
            }
        }

        EXPECT_EQ(atoms_added, 17);

        std::vector<const gemmi::Atom*> non_H = bond_info.non_Hydrogen_atoms_gemmi();
        int n_non_H = bond_info.n_non_H();
        int n_H = bond_info.n_H();

        EXPECT_EQ(n_non_H, 7);
        EXPECT_EQ(n_H, 10);

        // Check non-H atom names (pad gemmi's unpadded names to match MMDB's 4-char format)
        EXPECT_EQ(pad_atom_name(non_H[0]->name), " N  ");
        EXPECT_EQ(pad_atom_name(non_H[1]->name), " CA ");
        EXPECT_EQ(pad_atom_name(non_H[2]->name), " C  ");
        EXPECT_EQ(pad_atom_name(non_H[3]->name), " O  ");
        EXPECT_EQ(pad_atom_name(non_H[4]->name), " CB ");

        // Check hydrogen atoms
        std::vector<const gemmi::Atom*> H_atoms = bond_info.Hydrogen_atoms();
        EXPECT_EQ(pad_atom_name(H_atoms[0]->name), " H  ");
        EXPECT_EQ(pad_atom_name(H_atoms[1]->name), " HA ");
        EXPECT_EQ(pad_atom_name(H_atoms[2]->name), " HB2");
        EXPECT_EQ(pad_atom_name(H_atoms[3]->name), " HB3");
        EXPECT_EQ(pad_atom_name(H_atoms[4]->name), " HG2");
    }

    // Case 2: empty — no atoms added
    {
        coot::model_bond_atom_info_t_gemmi bond_info_empty;
        std::vector<const gemmi::Atom*> non_H_empty = bond_info_empty.non_Hydrogen_atoms_gemmi();
        int n_empty = bond_info_empty.n_non_H();

        EXPECT_EQ(n_empty, 0);
        EXPECT_TRUE(non_H_empty.empty());
    }

    // Case 3: add only a single non-hydrogen atom (CA from residue A/1)
    {
        coot::model_bond_atom_info_t_gemmi bond_info_single;
        
        const gemmi::Atom* ca_atom = find_atom_in_residue(st, "A", 1, "CA");
        EXPECT_NE(ca_atom, nullptr);

        if (ca_atom) {
            bond_info_single.add_bond_atom(*ca_atom);
            std::vector<const gemmi::Atom*> non_H_single = bond_info_single.non_Hydrogen_atoms_gemmi();
            EXPECT_EQ(bond_info_single.n_non_H(), 1);
            EXPECT_EQ(pad_atom_name(non_H_single[0]->name), " CA ");
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}