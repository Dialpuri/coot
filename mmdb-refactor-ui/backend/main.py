"""FastAPI application entry point."""
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

import report
import prompts
from config import REPORT_PATH, GTEST_FLAGS, args
from routers import files, refactor, progress, tests, preview

# ── Startup ───────────────────────────────────────────────────────────────────
report.load_and_update(REPORT_PATH)
prompts.init_docs(args.mmdb_docs)
print(f"GTest flags: {GTEST_FLAGS}")

# ── App ───────────────────────────────────────────────────────────────────────
app = FastAPI(title="MMDB Refactor Manager")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

app.include_router(files.router)
app.include_router(refactor.router)
app.include_router(progress.router)
app.include_router(tests.router)
app.include_router(preview.router)
