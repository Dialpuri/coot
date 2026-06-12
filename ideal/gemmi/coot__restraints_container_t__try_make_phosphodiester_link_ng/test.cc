#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/elem.hpp>
#include "function.hh"
#include <fstream>
#include <string>

static gemmi::Structure build_rna_structure() {
    gemmi::Structure st;
    gemmi::Model &model = st.models.emplace_back();
    gemmi::Chain &chain = model.chains.emplace_back();
    chain.name = "A";

    // Residue 1: A with P, O3', C1' atoms
    {
        gemmi::Residue &res = chain.residues.emplace_back();
        res.name = "A";
        res.seqid = gemmi::SeqId(1, ' ');
        {
            gemmi::Atom a; a.name = " P  "; a.pos = gemmi::Position(10,10,10); a.occ=1; a.b_iso=0; a.altloc='\0'; a.element = gemmi::Element("P");
            res.atoms.push_back(a);
        }
        {
            gemmi::Atom a; a.name = " O3'"; a.pos = gemmi::Position(11,11,11); a.occ=1; a.b_iso=0; a.altloc='\0'; a.element = gemmi::Element("O");
            res.atoms.push_back(a);
        }
        {
            gemmi::Atom a; a.name = " C1'"; a.pos = gemmi::Position(12,12,12); a.occ=1; a.b_iso=0; a.altloc='\0'; a.element = gemmi::Element("C");
            res.atoms.push_back(a);
        }
    }
    // Residue 2: U with P, O3', C1' atoms
    {
        gemmi::Residue &res = chain.residues.emplace_back();
        res.name = "U";
        res.seqid = gemmi::SeqId(2, ' ');
        {
            gemmi::Atom a; a.name = " P  "; a.pos = gemmi::Position(12.5,11.5,11.5); a.occ=1; a.b_iso=0; a.altloc='\0'; a.element = gemmi::Element("P");
            res.atoms.push_back(a);
        }
        {
            gemmi::Atom a; a.name = " O3'"; a.pos = gemmi::Position(13,13,13); a.occ=1; a.b_iso=0; a.altloc='\0'; a.element = gemmi::Element("O");
            res.atoms.push_back(a);
        }
        {
            gemmi::Atom a; a.name = " C1'"; a.pos = gemmi::Position(14,14,14); a.occ=1; a.b_iso=0; a.altloc='\0'; a.element = gemmi::Element("C");
            res.atoms.push_back(a);
        }
    }
    // Residue 3: ALA
    {
        gemmi::Residue &res = chain.residues.emplace_back();
        res.name = "ALA";
        res.seqid = gemmi::SeqId(3, ' ');
        {
            gemmi::Atom a; a.name = " CA"; a.pos = gemmi::Position(15,15,15); a.occ=1; a.b_iso=0; a.altloc='\0'; a.element = gemmi::Element("C");
            res.atoms.push_back(a);
        }
    }
    // Residue 4: ALA
    {
        gemmi::Residue &res = chain.residues.emplace_back();
        res.name = "ALA";
        res.seqid = gemmi::SeqId(4, ' ');
        {
            gemmi::Atom a; a.name = " CA"; a.pos = gemmi::Position(16,16,16); a.occ=1; a.b_iso=0; a.altloc='\0'; a.element = gemmi::Element("C");
            res.atoms.push_back(a);
        }
    }
    return st;
}

TEST(OracleTest, try_make_phosphodiester_link_ng) {
    gemmi::Structure st = build_rna_structure();

    EXPECT_FALSE(st.models.empty());
    gemmi::Model &model = st.models[0];
    EXPECT_FALSE(model.chains.empty());
    gemmi::Chain &chain = model.chains[0];

    EXPECT_GE(chain.residues.size(), 4u);

    const gemmi::Residue &res1 = chain.residues[0];  // A
    const gemmi::Residue &res2 = chain.residues[1];  // U
    const gemmi::Residue &res3 = chain.residues[2];  // ALA
    const gemmi::Residue &res4 = chain.residues[3];  // ALA

    // Case 1: phosphodiester link between A(res1) -> U(res2)
    {
        std::pair<bool, const gemmi::Residue *> pair1(false, &res1);
        std::pair<bool, const gemmi::Residue *> pair2(false, &res2);

        coot::protein_geometry pg;
        pg.init_standard();

        auto result = try_make_phosphodiester_link_ng_gemmi(pg, pair1, pair2);

        EXPECT_TRUE(result.first);
        EXPECT_EQ(result.second.n_link_bond_restr, 0u);
        EXPECT_EQ(result.second.n_link_angle_restr, 0u);
        EXPECT_EQ(result.second.n_link_torsion_restr, 0u);
    }

    // Case 2: same nucleotides reversed (U -> A)
    {
        std::pair<bool, const gemmi::Residue *> pair1(false, &res2);
        std::pair<bool, const gemmi::Residue *> pair2(false, &res1);

        coot::protein_geometry pg;
        pg.init_standard();

        auto result = try_make_phosphodiester_link_ng_gemmi(pg, pair1, pair2);

        EXPECT_TRUE(result.first);
        EXPECT_EQ(result.second.n_link_bond_restr, 0u);
    }

    // Case 3: non-nucleotide residues (ALA, ALA)
    {
        std::pair<bool, const gemmi::Residue *> pair1(false, &res3);
        std::pair<bool, const gemmi::Residue *> pair2(false, &res4);

        coot::protein_geometry pg;
        pg.init_standard();

        auto result = try_make_phosphodiester_link_ng_gemmi(pg, pair1, pair2);

        EXPECT_FALSE(result.first);
        EXPECT_EQ(result.second.n_link_bond_restr, 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
