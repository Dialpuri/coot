#!/bin/bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "Starting MMDB Refactor Manager..."
echo "SCRIPT_DIR=${SCRIPT_DIR}"

echo "Sourcing Python virtualenv"
source .venv/bin/activate

# Start backend in background
echo "Installing backend dependencies..."
cd "$SCRIPT_DIR/backend"
pip install -r requirements.txt -q
echo "Starting backend on http://localhost:8000 ..."
uvicorn main:app --reload --port 8000 &
BACKEND_PID=$!

echo ""

# Start frontend
echo "Installing frontend dependencies..."
cd "$SCRIPT_DIR/frontend"
npm install --silent
echo "Starting frontend on http://localhost:5173 ..."
npm run dev

# Kill backend when frontend exits
kill $BACKEND_PID 2>/dev/null || true
