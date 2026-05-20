#include <coot/utils/coot-utils.hh>

namespace coot {
namespace util {

float random_f() {
    return static_cast<float>(::random()) / static_cast<float>(RAND_MAX);
}

}
}