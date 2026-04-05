# MMDB Refactor Manager

A full-stack tool for managing the refactoring of the Coot codebase from MMDB to Gemmi. Browse files and functions with MMDB references, view source code, generate Gemmi refactors using a local Ollama LLM, and track progress.

## Prerequisites

- Python 3.10+
- Node.js 18+
- [Ollama](https://ollama.ai) running locally with a code model (e.g. `codellama`)

## Quick Start

```bash
cd /Users/dialpuri/lmb/coot/mmdb-refactor-ui
bash start.sh
```

This starts:
- Backend at http://localhost:8000
- Frontend at http://localhost:5173

## Manual Start

### Backend

```bash
cd backend
pip install -r requirements.txt
uvicorn main:app --reload
```

### Frontend

```bash
cd frontend
npm install
npm run dev
```

## Usage

1. The sidebar (left) lists all files with MMDB references, sorted by ref count. Use the search box to filter by path.
2. Click a file to load its functions in the middle panel.
3. Click a function to open it in the right panel, where the original MMDB code is displayed with symbols highlighted.
4. Set the Ollama model name (default: `codellama`) and optionally add extra instructions.
5. Click **Generate Gemmi Version** to stream a refactored version from Ollama.
6. Use **Mark Done** or **Mark Skipped** to track progress. Progress is saved to `backend/progress.json`.

## API Endpoints

| Method | Path | Description |
|--------|------|-------------|
| GET | /api/stats | Global stats from the JSON report |
| GET | /api/files | Paginated file list with search |
| GET | /api/file/{path} | Full file record |
| GET | /api/source | Read actual source lines |
| POST | /api/refactor | Stream Gemmi refactor from Ollama |
| GET | /api/progress | Get saved progress |
| POST | /api/progress | Update a function's status |
