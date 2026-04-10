"""CLI argument parsing, global constants, and GTest flag detection."""
import argparse
import subprocess
from pathlib import Path

# ── CLI args ──────────────────────────────────────────────────────────────────
_parser = argparse.ArgumentParser(description="MMDB Refactor Manager backend")
_parser.add_argument(
    "--report",
    default="/Users/dialpuri/lmb/coot/mmdb-recon/mmdb_usage_report.json",
    help="Path to mmdb_usage_report.json produced by the AST tool",
)
_parser.add_argument(
    "--cxx", default="c++",
    help="C++ compiler to use for test compilation",
)
_parser.add_argument(
    "--cxx-flags", default="",
    help="Extra flags to pass to the compiler (e.g. -I/extra/include)",
)
_parser.add_argument(
    "--mmdb-docs",
    default="/Users/dialpuri/lmb/coot/mmdb-recon/llm/mmdb_methods.md",
    help="Path to the mmdb_methods.md file produced by document_methods.py",
)
_parser.add_argument(
    "--probe-pdb",
    default="/Users/dialpuri/lmb/coot/reference-structures/1c7k.pdb",
    help="Path to a real PDB file used by oracle probes to measure function output",
)
_parser.add_argument(
    "--stream-llm",
    action="store_true",
    help="DEV: stream LLM responses in the probe pipeline so each token is "
         "visible in the console as it arrives (useful for debugging; adds a "
         "small per-chunk framing overhead but no extra compute).",
)
args, _ = _parser.parse_known_args()

# ── Paths & URLs ──────────────────────────────────────────────────────────────
REPORT_PATH   = args.report
PROGRESS_PATH = Path(__file__).parent / "progress.json"
TESTS_PATH    = Path(__file__).parent / "tests.json"
OLLAMA_URL    = "http://localhost:11434/api/generate"

# ── Compile constants ─────────────────────────────────────────────────────────
COOT_API_DIR  = "/Users/dialpuri/lmb/build-coot-and-deps/"
COOT_API_NAME = "cootapi"
MMDB_API_DIR  = "/opt/homebrew/Cellar/mmdb2/2.0.22/lib"
MMDB_API_NAME = "mmdb2"
GEMMI_INCLUDE_DIR = "/opt/homebrew/opt/gemmi/include/"

MAX_TEST_RETRIES  = 5
PROBE_PDB_PATH    = args.probe_pdb
MAX_PROBE_RETRIES = 3
DEV_STREAM_LLM    = args.stream_llm

# ── GTest detection ───────────────────────────────────────────────────────────
def _detect_gtest_flags() -> str:
    """Try pkg-config first, then fall back to common install prefixes."""
    try:
        flags = subprocess.check_output(
            ["pkg-config", "--cflags", "--libs", "gtest", "gtest_main"],
            text=True, stderr=subprocess.DEVNULL,
        ).strip()
        if flags:
            return flags
    except (FileNotFoundError, subprocess.CalledProcessError):
        pass
    for prefix in ["/opt/homebrew", "/usr/local", "/usr"]:
        if Path(f"{prefix}/include/gtest/gtest.h").exists():
            return f"-I{prefix}/include -L{prefix}/lib -lgtest -lgtest_main"
    return "-lgtest -lgtest_main"


GTEST_FLAGS = _detect_gtest_flags()
