"""Endpoints: GET/POST /api/progress"""
from fastapi import APIRouter

from models import ProgressUpdate
from storage import load_progress, save_progress

router = APIRouter()


@router.get("/api/progress")
def get_progress():
    return load_progress()


@router.post("/api/progress")
def update_progress(update: ProgressUpdate):
    data = load_progress()
    data[update.key] = update.status
    save_progress(data)
    return {"ok": True, "key": update.key, "status": update.status}
