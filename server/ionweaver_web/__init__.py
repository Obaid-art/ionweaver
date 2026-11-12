"""IonWeaver on the web: the submitted engine behind a small HTTP API.

Nothing in this package decides chemistry.  Every definition, product and
certificate comes from the engine modules in ``engine/iw``, which are the
submitted notebook's own source.  This package only loads the published state,
asks the engine the questions a user asks, and writes the answers out.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

__version__ = "2.0.0"

ROOT = Path(__file__).resolve().parents[2]
ENGINE_DIR = ROOT / "engine"
DATA_DIR = ROOT / "data" / "archive"
WEB_DIR = ROOT / "web"

if str(ENGINE_DIR) not in sys.path:
    sys.path.insert(0, str(ENGINE_DIR))

MANIFEST = json.loads((ENGINE_DIR / "MANIFEST.json").read_text(encoding="utf-8"))


def engine_fingerprint() -> dict:
    """Re-hash every engine module and say whether it is the notebook's."""
    out = {"modules": {}, "unmodified": True, "source_notebook": MANIFEST["source_notebook"],
           "vendored_utc": MANIFEST["vendored_utc"], "pins": MANIFEST["pins"]}
    for name, info in MANIFEST["modules"].items():
        path = ENGINE_DIR / "iw" / info["file"]
        digest = hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else ""
        ok = digest == info["sha256"]
        out["modules"][name] = {"sha256": digest, "unmodified": ok, "lines": info["lines"]}
        out["unmodified"] = out["unmodified"] and ok
    return out
