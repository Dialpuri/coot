#include "function.hh"

namespace coot {
namespace util {

bool copy_headers_gemmi(gemmi::Structure* st1, gemmi::Structure* st2, bool include_cryst) {
    if (!st1 || !st2) {
        return false;
    }

    st2->raw_remarks = st1->raw_remarks;
    st2->meta = st1->meta;

    if (include_cryst) {
        st2->cell = st1->cell;
    }

    return true;
}

}
}