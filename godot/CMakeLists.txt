# Godot GDExtension bindings for libcootapi.
#
# This is added by the top-level Coot CMakeLists when -DBUILD_GODOT_BINDINGS=ON,
# so it runs with the cootapi target and the dependency variables
# (${MMDB2_INCLUDE_DIR}, ${GLM_INCLUDE_DIRS}, ...) already in scope.
#
# It produces `libcoot_godot`, a GDExtension shared library that converts Coot's
# simple_mesh_t into a Godot ArrayMesh at runtime.
#
# Configure example:
#   cmake -B build -DBUILD_GODOT_BINDINGS=ON \
#         -DGODOT_CPP_DIR=/Users/jdialpuri/lmb/vr/repos/godot-cpp \
#         -DGODOT_EXTENSION_API_FILE=/Users/jdialpuri/lmb/vr/extension_api.json

# --- godot-cpp -------------------------------------------------------------
set(GODOT_CPP_DIR "${coot_src}/../godot-cpp"
    CACHE PATH "Path to the godot-cpp checkout")
if(NOT EXISTS "${GODOT_CPP_DIR}/CMakeLists.txt")
    message(FATAL_ERROR
        "godot-cpp not found at '${GODOT_CPP_DIR}'. "
        "Pass -DGODOT_CPP_DIR=/path/to/godot-cpp")
endif()

# Build godot-cpp against the exact API dumped from the target editor (recommended
# when the editor is newer than any godot-cpp release branch, e.g. 4.7.x).
# Dump it with:  godot --dump-extension-api
set(GODOT_EXTENSION_API_FILE "" CACHE FILEPATH
    "extension_api.json dumped from your Godot editor (exact version match)")
if(GODOT_EXTENSION_API_FILE)
    set(GODOTCPP_CUSTOM_API_FILE "${GODOT_EXTENSION_API_FILE}" CACHE FILEPATH "" FORCE)
    message(STATUS "coot_godot: building godot-cpp against ${GODOT_EXTENSION_API_FILE}")
endif()

add_subdirectory(${GODOT_CPP_DIR} ${CMAKE_CURRENT_BINARY_DIR}/godot-cpp-build)

# --- the extension ---------------------------------------------------------
add_library(coot_godot SHARED
    src/register_types.cpp
    src/coot_molecule.cpp)

target_compile_features(coot_godot PRIVATE cxx_std_17)

# Mirror the include set that cootapi itself compiles with, so molecules-container.hh
# and its transitive headers resolve.
target_include_directories(coot_godot PRIVATE
    ${coot_src}
    ${coot_src}/coords
    ${coot_src}/utils
    ${coot_src}/ccp4mg-utils
    ${coot_src}/ccp4mg-utils/util
    ${coot_src}/ccp4mg-utils/pygl
    ${MMDB2_INCLUDE_DIR}
    ${SSM_INCLUDE_DIR}
    ${CLIPPER-CORE_INCLUDE_DIR}
    ${CLIPPER-MMDB_INCLUDE_DIR}
    ${CLIPPER-CCP4_INCLUDE_DIR}
    ${FFTW2_INCLUDE_DIRS}
    ${RDKit_INCLUDE_DIRS}
    ${GLM_INCLUDE_DIRS}
    ${GSL_INCLUDE_DIRS})

# We only call into cootapi (which carries its own transitive deps as a shared lib).
# glm is header-only and used directly by the converter.
target_link_libraries(coot_godot PRIVATE godot-cpp cootapi glm::glm)

# Godot expects the "lib" prefix on all platforms (matches coot.gdextension).
set_target_properties(coot_godot PROPERTIES
    PREFIX "lib"
    CXX_VISIBILITY_PRESET hidden
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib"
    INSTALL_RPATH_USE_LINK_PATH TRUE)

message(STATUS "coot_godot: GDExtension bindings enabled")
