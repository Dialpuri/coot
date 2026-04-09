"""Helpers for test file paths, naming, and content wrapping."""
import re
from pathlib import Path

import report

_VARIANT_INCLUDES: dict[str, list[str]] = {
    "mmdb":  ["#include <mmdb2/mmdb_manager.h>"],
    "gemmi": ["#include <gemmi/structure.hpp>", "#include <gemmi/model.hpp>"],
}


def sanitize_fn_name(name: str) -> str:
    """coot::util::foo_bar  ->  coot_util_foo_bar"""
    return re.sub(r'[^a-zA-Z0-9]', '_', name).strip('_')


def get_test_file_path(rel_source_path: str, fn_name: str, variant: str) -> Path:
    """Return the absolute path where a test .cc file should be written."""
    source_dir = Path(rel_source_path).parent          # e.g. coot-utils
    tests_dir  = Path(report._coot_root) / source_dir / "tests"
    safe       = sanitize_fn_name(fn_name)
    return tests_dir / f"{safe}_{variant}_test.cc"


def wrap_test_content(content: str, variant: str, fn_name: str) -> str:
    """Ensure the .cc file has a gtest header, variant includes, and a main()."""
    content = content.strip()
    header: list[str] = [f"// Auto-generated {variant.upper()} test for {fn_name}"]
    if "#include <gtest/gtest.h>" not in content:
        header.append("#include <gtest/gtest.h>")
    for inc in _VARIANT_INCLUDES.get(variant, []):
        if inc not in content:
            header.append(inc)
    needs_main = "RUN_ALL_TESTS()" not in content and "int main(" not in content
    parts = ["\n".join(header), "", content]
    if needs_main:
        parts += ["", "int main(int argc, char **argv) {",
                  "  ::testing::InitGoogleTest(&argc, argv);",
                  "  return RUN_ALL_TESTS();", "}"]
    return "\n".join(parts) + "\n"
