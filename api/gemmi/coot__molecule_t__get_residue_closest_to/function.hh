#pragma once
#include <string>
#include <gemmi/pdb.hpp>
#include <clipper/clipper.h>

namespace coot {

struct residue_spec_t {
    int model_number;
    std::string chain_id;
    int res_no;
    std::string ins_code;
    int int_user_data;
    float float_user_data;
    std::string string_user_data;

    explicit residue_spec_t(int r) : res_no(r) {
        model_number = -1;
        int_user_data = -1;
        float_user_data = -1;
    }

    residue_spec_t(const std::string &chain_in, int r) : chain_id(chain_in) {
        model_number = -1;
        res_no = r;
        int_user_data = -1;
        float_user_data = -1;
    }

    residue_spec_t(int model_number_in,
                   const std::string &chain_in, int r,
                   const std::string &ins_code_in) : chain_id(chain_in), ins_code(ins_code_in) {
        model_number = model_number_in;
        res_no = r;
        int_user_data = -1;
        float_user_data = -1;
    }

    residue_spec_t(const std::string &chain_in, int r,
                   const std::string &ins_code_in) : chain_id(chain_in), ins_code(ins_code_in) {
        model_number = -1;
        res_no = r;
        int_user_data = -1;
        float_user_data = -1;
    }

    // Constructor from Chain* and Residue* (gemmi-style iteration)
    residue_spec_t(const gemmi::Chain* chain, const gemmi::Residue* res) {
        if (!res || !chain) {
            chain_id = "";
            model_number = -1;
            res_no = -1;
            ins_code = "";
        } else {
            chain_id = chain->name;
            model_number = 1;             // gemmi models are 0-indexed, but we report 1
            res_no = res->seqid.num.value;
            // MMDB uses "" for no insertion code; gemmi uses ' ' (space)
            // Convert space to empty string for compatibility
            ins_code = (res->seqid.icode == ' ') ? std::string("") : std::string(1, res->seqid.icode);
        }
        int_user_data = -1;
    }

    residue_spec_t() : model_number(-1), res_no(-1) {}
};

class molecule_t {
public:
    residue_spec_t get_residue_closest_to_gemmi(const gemmi::Structure& st,
                                                const clipper::Coord_orth& co) const;
};

} // namespace coot