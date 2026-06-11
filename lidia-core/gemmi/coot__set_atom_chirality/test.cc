#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include "geometry/protein-geometry.hh"
#include <rdkit/GraphMol/RDKitBase.h>

TEST(OracleTest, set_atom_chirality) {
    {
        // Load structure with gemmi
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

        // Get the first model
        gemmi::Model& model = st.models[0];

        // Get first chain
        gemmi::Chain& chain = model.chains[0];

        // Get first residue
        gemmi::Residue* res = &chain.residues[0];
        EXPECT_NE(res, nullptr);

        std::string res_name(res->name);
        int seqnum = res->seqid.num.value;
        std::cout << "residue: " << res_name << " " << seqnum << std::endl;

        // Get first atom from the residue
        EXPECT_GT(res->atoms.size(), 0);
        gemmi::Atom* at = &res->atoms[0];
        std::string atom_name(at->name);
        std::cout << "atom_name: " << atom_name << std::endl;

        // Get restraints for this residue type
        coot::protein_geometry geom;
        geom.init_standard();
        std::pair<bool, coot::dictionary_residue_restraints_t> result =
            geom.get_monomer_restraints_at_least_minimal(res_name, 1);
        EXPECT_TRUE(result.first);
        coot::dictionary_residue_restraints_t restraints = result.second;
        std::cout << "n_chiral_restraints: " << restraints.chiral_restraint.size() << std::endl;
        std::cout << "n_atom_info: " << restraints.atom_info.size() << std::endl;

        // Create RDKit atom
        RDKit::Atom *rdkit_at = new RDKit::Atom;
        rdkit_at->setAtomicNum(6); // carbon
        rdkit_at->setProp("name", atom_name);

        // Record chirality before
        RDKit::Atom::ChiralType chiral_before = rdkit_at->getChiralTag();
        std::string chiral_before_str = "!";
        if (chiral_before == RDKit::Atom::CHI_UNSPECIFIED)     chiral_before_str = "UNSPECIFIED";
        if (chiral_before == RDKit::Atom::CHI_TETRAHEDRAL_CW)  chiral_before_str = "CW";
        if (chiral_before == RDKit::Atom::CHI_TETRAHEDRAL_CCW) chiral_before_str = "CCW";
        if (chiral_before == RDKit::Atom::CHI_OTHER)           chiral_before_str = "OTHER";
        std::cout << "chiral_before: " << chiral_before_str << std::endl;
        EXPECT_EQ(chiral_before_str, "UNSPECIFIED");

        // Call the function
        EXPECT_NO_THROW(coot::set_atom_chirality_gemmi(rdkit_at, at, res, restraints));

        // Record chirality after
        RDKit::Atom::ChiralType chiral_after = rdkit_at->getChiralTag();
        std::string chiral_after_str = "!";
        if (chiral_after == RDKit::Atom::CHI_UNSPECIFIED)     chiral_after_str = "UNSPECIFIED";
        if (chiral_after == RDKit::Atom::CHI_TETRAHEDRAL_CW)  chiral_after_str = "CW";
        if (chiral_after == RDKit::Atom::CHI_TETRAHEDRAL_CCW) chiral_after_str = "CCW";
        if (chiral_after == RDKit::Atom::CHI_OTHER)           chiral_after_str = "OTHER";
        std::cout << "chiral_after: " << chiral_after_str << std::endl;

        // Check if any properties were set
        std::string prop_n1, prop_sign;
        bool has_n1 = (rdkit_at->hasProp("mmcif_chiral_N1"));
        bool has_sign = (rdkit_at->hasProp("mmcif_chiral_volume_sign"));
        if (has_n1) rdkit_at->getProp("mmcif_chiral_N1", prop_n1);
        if (has_sign) rdkit_at->getProp("mmcif_chiral_volume_sign", prop_sign);
        std::cout << "has_chiral_props: " << (has_n1 ? "true" : "false") << std::endl;
        if (has_n1) std::cout << "chiral_N1: " << prop_n1 << std::endl;
        if (has_sign) std::cout << "chiral_volume_sign: " << prop_sign << std::endl;

        EXPECT_FALSE(chiral_after_str.empty());

        delete rdkit_at;
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
