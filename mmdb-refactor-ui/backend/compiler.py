"""Async helpers for compiling and running GTest binaries."""
import asyncio
import re
from pathlib import Path

import report
from config import args, COOT_API_DIR, COOT_API_NAME, MMDB_API_DIR, MMDB_API_NAME, GTEST_FLAGS, GEMMI_INCLUDE_DIR


def make_compile_cmd(test_file: Path, output_bin: Path) -> str:
    return (
        f'{args.cxx} -std=c++17 "{test_file}" -o "{output_bin}" '
        f'{GTEST_FLAGS} -I"{report._coot_root}" -pthread '
        f'-Wl,-rpath,{COOT_API_DIR} -L {COOT_API_DIR} -L {MMDB_API_DIR} '
        f'-l{COOT_API_NAME} -l {MMDB_API_NAME}'
    )


def make_probe_compile_cmd(probe_file: Path, output_bin: Path) -> str:
    """Like make_compile_cmd but omits gtest_main so the probe's own main() links."""
    probe_flags = re.sub(r'(^|\s)-lgtest_main(\s|$)', ' ', GTEST_FLAGS).strip()
    return (
        f'{args.cxx} -std=c++17 "{probe_file}" -o "{output_bin}" '
        f'{probe_flags} -I"{report._coot_root}" -I"{GEMMI_INCLUDE_DIR}" -pthread '
        f'-Wl,-rpath,{COOT_API_DIR} -L {COOT_API_DIR} -L {MMDB_API_DIR} '
        f'-l{COOT_API_NAME} -l {MMDB_API_NAME}'
    )


async def compile_test(test_file: Path) -> tuple[bool, str]:
    """Compile a test file; return (success, output)."""
    bin_dir = test_file.parent / "bin"
    bin_dir.mkdir(exist_ok=True)
    output_bin = bin_dir / test_file.stem
    cmd = make_compile_cmd(test_file, output_bin)
    if args.cxx_flags:
        cmd += f" {args.cxx_flags}"
    proc = await asyncio.create_subprocess_shell(
        cmd, stdout=asyncio.subprocess.PIPE, stderr=asyncio.subprocess.STDOUT,
        cwd=str(report._coot_root),
    )
    stdout, _ = await proc.communicate()
    return proc.returncode == 0, stdout.decode(errors="replace")


async def run_test(test_file: Path) -> tuple[bool, str]:
    """Run a compiled test binary; return (success, output)."""
    output_bin = test_file.parent / "bin" / test_file.stem
    if not output_bin.exists():
        return False, "binary not found"
    run_proc = await asyncio.create_subprocess_exec(
        str(output_bin), stdout=asyncio.subprocess.PIPE, stderr=asyncio.subprocess.STDOUT,
    )
    stdout, _ = await run_proc.communicate()
    return run_proc.returncode == 0, stdout.decode(errors="replace")
