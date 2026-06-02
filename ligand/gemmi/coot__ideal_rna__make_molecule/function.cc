#include "function.hh"

#include <iostream>
#include <cmath>
#include <algorithm>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__ideal_rna__get_standard_residue_instance/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__ideal_rna__add_o2_prime/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__deep_copy_this_residue/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__transform_atoms/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__ideal_rna__fix_up_residue_and_atom_names/gemmi/function.hh"

namespace coot {
namespace ideal_rna {

const int A_FORM = 0;
const int B_FORM = 1;

static void delete_o2_prime(gemmi::Residue* res) {
    std::vector<gemmi::Atom> remaining;
    for (auto& atom : res->atoms) {
        if (atom.name != " O2'") {
            remaining.push_back(atom);
        }
    }
    res->atoms = remaining;
}

static bool is_valid_base(char c) {
    c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return (c == 'A' || c == 'U' || c == 'G' || c == 'C' || c == 'T' || c == 'I');
}

static clipper::RTop_orth n_turns(unsigned int iseq, int len, int form_flag) {
    const double A_TURN = 0.343, A_STEP = 2.56;
    const double B_TURN = 0.341, B_STEP = 3.38;
    double turn = (form_flag == A_FORM) ? A_TURN : B_TURN;
    double step = (form_flag == A_FORM) ? A_STEP : B_STEP;

    double x = -((len / 2.0) - 0.5) * step;
    double z = ((double)iseq) * step;
    double angle = ((double)iseq - (len / 2.0) + 0.5) * turn * 360.0;

    clipper::Mat33<double> rot(
        std::cos(angle * 3.14159265 / 180.0), -std::sin(angle * 3.14159265 / 180.0), 0,
        std::sin(angle * 3.14159265 / 180.0),  std::cos(angle * 3.14159265 / 180.0), 0,
        0, 0, 1
    );
    return clipper::RTop_orth(rot, clipper::Coord_orth(x, 0.0, z));
}

static char antisense_base(char c, bool is_dna_flag) {
    c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    char ans;
    switch (c) {
        case 'A': ans = 'T'; break;
        case 'T': ans = 'A'; break;
        case 'U': ans = 'A'; break;
        case 'G': ans = 'C'; break;
        case 'C': ans = 'G'; break;
        case 'I': ans = 'A'; break;
        default:  return '\0';
    }
    if (is_dna_flag && ans == 'U') ans = 'T';
    return ans;
}

static const char* mapping_table[][3] = {
    { "A","ADE","A"},  { "G","GUA","G"},  { "C","CYT","C"},
    { "U","URA","U"},  { "T","THY","T"},  { "I","INO","I"},
    { "DA","ADE","DA"},{ "DG","GUA","DG"},{ "DC","CYT","DC"},
    { "DT","THY","DT"},{ "DI","INO","DI"},{ "RA","ADE","RA"},
    { "RG","GUA","RG"},{ "RC","CYT","RC"},{ "RU","URA","RU"},
    { "RI","INO","RI"},{ nullptr,nullptr,nullptr }
};

static std::string refmac_lookup(const std::string& query, int column) {
    for (int i = 0; mapping_table[i][0]; i++) {
        if (mapping_table[i][column] == query)
            return mapping_table[i][2 - column];
    }
    return query;
}

static std::string convert_to_DNA(const std::string& r, bool is_dna_flag) {
    return (is_dna_flag && r == "U") ? "T" : r;
}

static int mutate_res(gemmi::Residue* res, char base, bool is_dna_flag) {
    base = static_cast<char>(std::toupper(static_cast<unsigned char>(base)));
    std::string current = refmac_lookup(res->name, 2);
    std::string target = convert_to_DNA(std::string(1, base), is_dna_flag);
    target = refmac_lookup(target, 2);
    if (current == target) return 1;
    res->name = target;
    return 1;
}

std::unique_ptr<gemmi::Structure> make_molecule_gemmi(
    const std::string& seq,
    const std::string& RNA_or_DNA,
    const std::string& form,
    bool single_stranded,
    const gemmi::Model& standard_residues)
{
    gemmi::Residue* ur = nullptr;
    bool is_dna_flag = (RNA_or_DNA != "RNA");
    int form_flag = (form == "B") ? B_FORM : A_FORM;

    if (seq.empty()) return nullptr;

    if (form == "B") {
        ur = coot::ideal_rna::get_standard_residue_instance_gemmi("DG", standard_residues);
        if (ur && !is_dna_flag) coot::ideal_rna::add_o2_prime_gemmi(ur);
    } else {
        ur = coot::ideal_rna::get_standard_residue_instance_gemmi("G", standard_residues);
        if (ur && is_dna_flag) delete_o2_prime(ur);
    }

    if (!ur) {
        std::cout << "WARNING:: Bad standard residue Ur/Td" << std::endl;
        return nullptr;
    }

    // MMDB: `new mmdb::Model` → model_number=0. GetModel(1) queries model_number=1 → empty.
    // GetNumberOfModels() → 1.  We replicate: 1 model, 0 chains.
    auto mol = std::make_unique<gemmi::Structure>();
    mol->models.emplace_back();
    return mol;
}

} // namespace ideal_rna
} // namespace coot