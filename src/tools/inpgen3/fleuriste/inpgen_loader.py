"""Shared FleurInpgen loading helpers for the FLEURiste application."""

from __future__ import annotations

import importlib
import os
from pathlib import Path
from typing import Optional


def _fallback_search_paths() -> list[Path]:
    """Return default shared-library search paths used by FleurInpgen."""
    search_paths = [
        Path.cwd(),
        Path.cwd() / "build" / "lib",
        Path.cwd() / "build" / "src" / "tools" / "inpgen3",
        Path.cwd() / "lib",
        Path(__file__).resolve().parents[1] / "lib",
        Path(__file__).resolve().parents[4] / "build" / "lib",
    ]

    builddir = os.environ.get("FLEUR_BUILDDIR")
    if builddir:
        search_paths.insert(0, Path(builddir) / "src" / "tools" / "inpgen3")
        search_paths.insert(1, Path(builddir))

    return search_paths


def get_inpgen_search_paths() -> list[Path]:
    """Return the effective search paths for libfleurinpgen."""
    try:
        module = importlib.import_module("FleurInpgen")
    except Exception:
        return _fallback_search_paths()

    get_paths = getattr(module, "get_library_search_paths", None)
    if callable(get_paths):
        return list(get_paths())
    return _fallback_search_paths()


def inpgen_python_module_available() -> bool:
    """Return True if the FleurInpgen Python wrapper can be imported."""
    try:
        importlib.import_module("FleurInpgen")
    except Exception:
        return False
    return True


def format_inpgen_error_message(exc: Optional[BaseException] = None) -> str:
    """Build a consistent user-facing error message for missing/broken inpgen."""
    lines = ["FleurInpgen library not available."]

    if exc is not None:
        detail = str(exc).strip()
        if detail:
            lines.append(detail)

    detail_text = "\n".join(lines)
    if "Searched in:" not in detail_text:
        search_paths = "\n".join(f"  - {path}" for path in get_inpgen_search_paths())
        detail_text += f"\nSearched in:\n{search_paths}"

    return detail_text


def create_inpgen_interface(lib_path: Optional[str | Path] = None, quiet: bool = False):
    """Import FleurInpgen and instantiate its interface with consistent errors."""
    try:
        module = importlib.import_module("FleurInpgen")
    except Exception as exc:
        raise RuntimeError(format_inpgen_error_message(exc)) from exc

    try:
        return module.InpgenInterface(lib_path=lib_path, quiet=quiet)
    except Exception as exc:
        raise RuntimeError(format_inpgen_error_message(exc)) from exc