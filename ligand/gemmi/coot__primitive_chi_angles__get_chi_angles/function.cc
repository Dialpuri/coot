#include "function.hh"

static void ensure_chi_quads_initialized() {
    static bool initialized = [](){
        using namespace coot::primitive_chi_angles;
        add_chi_quad_table("VAL", " N  ", " CA ", " CB ", " CG1");
        add_chi_quad_table("TYR", " N  ", " CA ", " CB ", " CG ");
        add_chi_quad_table("TYR", " CA ", " CB ", " CG ", " CD1");
        add_chi_quad_table("TRP", " N  ", " CA ", " CB ", " CG ");
        add_chi_quad_table("TRP", " CA ", " CB ", " CG ", " CD1");
        add_chi_quad_table("THR", " N  ", " CA ", " CB ", " OG1");
        add_chi_quad_table("SER", " N  ", " CA ", " CB ", " OG ");
        add_chi_quad_table("PRO", " N  ", " CA ", " CB ", " CG ");
        add_chi_quad_table("PHE", " N  ", " CA ", " CB ", " CG ");
        add_chi_quad_table("PHE", " CA ", " CB ", " CG ", " CD1");
        add_chi_quad_table("MET", " N  ", " CA ", " CB ", " CG ");
        add_chi_quad_table("MET", " CA ", " CB ", " CG ", " SD ");
        add_chi_quad_table("MET", " CB ", " CG ", " SD ", " CE ");
        add_chi_quad_table("MSE", " N  ", " CA ", " CB ", " CG ");
        add_chi_quad_table("MSE", " CA ", " CB ", " CG ", "SE  ");
        add_chi_quad_table("MSE", " CB ", " CG ", "SE  ", " CE ");
        add_chi_quad_table("LYS", " N  ", " CA ", " CB ", " CG ");
        add_chi_quad_table("LYS", " CA ", " CB ", " CG ", " CD ");
        add_chi_quad_table("LYS", " CB ", " CG ", " CD ", " CE ");
        add_chi_quad_table("LYS", " CG ", " CD ", " CE ", " NZ ");
        add_chi_quad_table("LEU", " N  ", " CA ", " CB ", " CG ");
        add_chi_quad_table("LEU", " CA ", " CB ", " CG ", " CD1");
        add_chi_quad_table("ILE", " N  ", " CA ", " CB ", " CG1");
        add_chi_quad_table("ILE", " CA ", " CB ", " CG1", " CD1");
        add_chi_quad_table("HIS", " N  ", " CA ", " CB ", " CG ");
        add_chi_quad_table("HIS", " CA ", " CB ", " CG ", " ND1");
        add_chi_quad_table("GLU", " N  ", " CA ", " CB ", " CG ");
        add_chi_quad_table("GLU", " CA ", " CB ", " CG ", " CD ");
        add_chi_quad_table("GLU", " CB ", " CG ", " CD ", " OE1");
        add_chi_quad_table("GLN", " N  ", " CA ", " CB ", " CG ");
        add_chi_quad_table("GLN", " CA ", " CB ", " CG ", " CD ");
        add_chi_quad_table("GLN", " CB ", " CG ", " CD ", " OE1");
        add_chi_quad_table("CYS", " N  ", " CA ", " CB ", " SG ");
        add_chi_quad_table("ASP", " N  ", " CA ", " CB ", " CG ");
        add_chi_quad_table("ASP", " CA ", " CB ", " CG ", " OD1");
        add_chi_quad_table("ASN", " N  ", " CA ", " CB ", " CG ");
        add_chi_quad_table("ASN", " CA ", " CB ", " CG ", " OD1");
        add_chi_quad_table("ARG", " N  ", " CA ", " CB ", " CG ");
        add_chi_quad_table("ARG", " CA ", " CB ", " CG ", " CD ");
        add_chi_quad_table("ARG", " CB ", " CG ", " CD ", " NE ");
        add_chi_quad_table("ARG", " CG ", " CD ", " NE ", " CZ ");
        return true;
    }();
    (void)initialized;
}

static std::vector<coot::primitive_chi_angles::atom_name_quad>
get_atom_name_quads_for_type(const std::string& residue_type) {
    std::vector<coot::primitive_chi_angles::atom_name_quad> v;
    const auto& tbl = coot::primitive_chi_angles::chi_table();
    for (size_t i = 0; i < tbl.size(); i++) {
        if (tbl[i].residue_type == residue_type) {
            v = tbl[i].name_quad;
            break;
        }
    }
    return v;
}

std::vector<coot::alt_confed_chi_angles>
coot::primitive_chi_angles::get_chi_angles_gemmi(const gemmi::Residue& residue) {

    ensure_chi_quads_initialized();

    std::vector<coot::alt_confed_chi_angles> nv;

    // Check for alt confs
    bool residue_has_alt_confs = false;
    for (const auto& atom : residue.atoms) {
        if (atom.altloc != '\0') {
            residue_has_alt_confs = true;
            break;
        }
    }

    std::string residue_name = residue.name;
    std::vector<coot::primitive_chi_angles::atom_name_quad> atom_name_quad_list =
        get_atom_name_quads_for_type(residue_name);

    if (atom_name_quad_list.empty()) {
        std::string mess = "Failed to find atom name quads for residue type ";
        mess += residue_name;
        throw std::runtime_error(mess);
    }

    if (!residue_has_alt_confs) {
        // No alt confs — use get_quads_gemmi (takes coot::atom_name_quad)
        std::vector<coot::atom_name_quad> quad_list_top;
        for (const auto& q : atom_name_quad_list) {
            quad_list_top.push_back(coot::atom_name_quad(
                q.names[0], q.names[1], q.names[2], q.names[3]));
        }
        std::vector<coot::atom_index_quad> quads =
            coot::primitive_chi_angles::get_quads_gemmi(quad_list_top, residue);

        std::vector<std::pair<int, float> > v;
        for (unsigned int i_quad = 0; i_quad < quads.size(); i_quad++) {
            const coot::atom_index_quad& q = quads[i_quad];
            const auto& atoms = residue.atoms;
            clipper::Coord_orth p1(coot::primitive_chi_angles::atom_to_co_gemmi(atoms[q.index1]));
            clipper::Coord_orth p2(coot::primitive_chi_angles::atom_to_co_gemmi(atoms[q.index2]));
            clipper::Coord_orth p3(coot::primitive_chi_angles::atom_to_co_gemmi(atoms[q.index3]));
            clipper::Coord_orth p4(coot::primitive_chi_angles::atom_to_co_gemmi(atoms[q.index4]));
            double tors = clipper::Util::rad2d(clipper::Coord_orth::torsion(p1, p2, p3, p4));
            v.push_back(std::make_pair<int, float>(static_cast<int>(i_quad + 1), static_cast<float>(tors)));
        }

        if (!v.empty()) {
            nv.push_back(coot::alt_confed_chi_angles("", v));
        }
    } else {
        // Multiple quads from alt confed atoms
        std::vector<coot::primitive_chi_angles::alt_confed_atom_index_quad> quads_vec =
            coot::primitive_chi_angles::get_quads_using_altconfs(atom_name_quad_list, residue);

        for (unsigned int i_quad_set = 0; i_quad_set < quads_vec.size(); i_quad_set++) {
            std::vector<std::pair<int, float> > v;
            const std::vector<coot::primitive_chi_angles::atom_index_quad>& quad = quads_vec[i_quad_set].quad;
            for (unsigned int i_quad = 0; i_quad < quad.size(); i_quad++) {
                const auto& atoms = residue.atoms;
                clipper::Coord_orth p1(coot::primitive_chi_angles::atom_to_co_gemmi(atoms[quad[i_quad].index1]));
                clipper::Coord_orth p2(coot::primitive_chi_angles::atom_to_co_gemmi(atoms[quad[i_quad].index2]));
                clipper::Coord_orth p3(coot::primitive_chi_angles::atom_to_co_gemmi(atoms[quad[i_quad].index3]));
                clipper::Coord_orth p4(coot::primitive_chi_angles::atom_to_co_gemmi(atoms[quad[i_quad].index4]));
                double tors = clipper::Util::rad2d(clipper::Coord_orth::torsion(p1, p2, p3, p4));
                v.push_back(std::make_pair<int, float>(static_cast<int>(i_quad + 1), static_cast<float>(tors)));
            }

            if (!v.empty()) {
                nv.push_back(coot::alt_confed_chi_angles(quads_vec[i_quad_set].alt_conf, v));
            }
        }
    }

    return nv;
}