"""Persistence helpers for progress and test records."""
import json

from config import PROGRESS_PATH, TESTS_PATH


def load_progress() -> dict:
    if PROGRESS_PATH.exists():
        with open(PROGRESS_PATH, "r") as f:
            return json.load(f)
    return {}


def save_progress(data: dict) -> None:
    with open(PROGRESS_PATH, "w") as f:
        json.dump(data, f, indent=2)


def load_tests() -> dict:
    if TESTS_PATH.exists():
        with open(TESTS_PATH, "r") as f:
            return json.load(f)
    return {}


def save_tests(data: dict) -> None:
    with open(TESTS_PATH, "w") as f:
        json.dump(data, f, indent=2)
