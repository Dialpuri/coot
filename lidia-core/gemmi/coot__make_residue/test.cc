#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/elem.hpp>
#include "function.hh"
#include <rdkit/GraphMol/RWMol.h>
#include <rdkit/GraphMol/Atom.h>
#include <rdkit/GraphMol/Bond.h>
#include <rdkit/Geometry/point.h>

TEST(OracleTest, make_residue) {
    // Helper: create a simple 3-atom molecule (water) with explicit 3D coords
    auto make_water = []() -> RDKit::RWMol* {
        RDKit::RWMol *mol = new RDKit::RWMol();
        RDKit::Atom *o_at = new RDKit::Atom(8);
        unsigned int o_idx = mol->addAtom(o_at, false, false);
        RDKit::Atom *h1_at = new RDKit::Atom(1);
        unsigned int h1_idx = mol->addAtom(h1_at, false, false);
        RDKit::Atom *h2_at = new RDKit::Atom(1);
        unsigned int h2_idx = mol->addAtom(h2_at, false, false);
        mol->addBond(o_idx, h1_idx, RDKit::Bond::SINGLE);
        mol->addBond(o_idx, h2_idx, RDKit::Bond::SINGLE);
        RDKit::Conformer *conf = new RDKit::Conformer(3);
        conf->set3D(true);
        conf->setAtomPos(o_idx, RDGeom::Point3D(0.0, 0.0, 0.0));
        conf->setAtomPos(h1_idx, RDGeom::Point3D(0.76, 0.59, 0.0));
        conf->setAtomPos(h2_idx, RDGeom::Point3D(-0.76, 0.59, 0.0));
        mol->addConformer(conf, true);
        return mol;
    };

    // Helper: create a simple methane-like molecule (4 atoms) with explicit coords
    auto make_methane = []() -> RDKit::RWMol* {
        RDKit::RWMol *mol = new RDKit::RWMol();
        RDKit::Atom *c_at = new RDKit::Atom(6);
        unsigned int c_idx = mol->addAtom(c_at, false, false);
        RDKit::Atom *h1_at = new RDKit::Atom(1);
        unsigned int h1_idx = mol->addAtom(h1_at, false, false);
        RDKit::Atom *h2_at = new RDKit::Atom(1);
        unsigned int h2_idx = mol->addAtom(h2_at, false, false);
        RDKit::Atom *h3_at = new RDKit::Atom(1);
        unsigned int h3_idx = mol->addAtom(h3_at, false, false);
        RDKit::Atom *h4_at = new RDKit::Atom(1);
        unsigned int h4_idx = mol->addAtom(h4_at, false, false);
        mol->addBond(c_idx, h1_idx, RDKit::Bond::SINGLE);
        mol->addBond(c_idx, h2_idx, RDKit::Bond::SINGLE);
        mol->addBond(c_idx, h3_idx, RDKit::Bond::SINGLE);
        mol->addBond(c_idx, h4_idx, RDKit::Bond::SINGLE);
        RDKit::Conformer *conf = new RDKit::Conformer(4);
        conf->set3D(true);
        conf->setAtomPos(c_idx, RDGeom::Point3D(0.0, 0.0, 0.0));
        conf->setAtomPos(h1_idx, RDGeom::Point3D(0.63, 0.63, 0.63));
        conf->setAtomPos(h2_idx, RDGeom::Point3D(-0.63, -0.63, 0.63));
        conf->setAtomPos(h3_idx, RDGeom::Point3D(-0.63, 0.63, -0.63));
        conf->setAtomPos(h4_idx, RDGeom::Point3D(0.63, -0.63, -0.63));
        mol->addConformer(conf, true);
        return mol;
    };

    // Case 1: Water molecule (3 atoms)
    {
        RDKit::RWMol *mol = make_water();
        int iconf = mol->getConformer().getId();
        std::string res_name = "WAT";
        EXPECT_EQ(iconf, 0);
        gemmi::Residue *res = coot::make_residue_gemmi(*mol, iconf, res_name);
        EXPECT_NE(res, nullptr);
        if (res) {
            EXPECT_EQ(res->atoms.size(), 3u);
            EXPECT_STREQ(res->name.c_str(), "WAT");
            EXPECT_EQ(res->seqid.num.value, 1);
            delete res;
        }
        delete mol;
    }

    // Case 2: Empty molecule (no atoms) — should return NULL
    {
        RDKit::RWMol empty_mol;
        std::string res_name = "EMPTY";
        gemmi::Residue *res = coot::make_residue_gemmi(empty_mol, 0, res_name);
        EXPECT_EQ(res, nullptr);
    }

    // Case 3: Methane molecule (5 atoms: C + 4H)
    {
        RDKit::RWMol *mol = make_methane();
        int iconf = mol->getConformer().getId();
        std::string res_name = "MTN";
        EXPECT_EQ(iconf, 0);
        gemmi::Residue *res = coot::make_residue_gemmi(*mol, iconf, res_name);
        EXPECT_NE(res, nullptr);
        if (res) {
            EXPECT_EQ(res->atoms.size(), 5u);
            EXPECT_STREQ(res->name.c_str(), "MTN");
            delete res;
        }
        delete mol;
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}