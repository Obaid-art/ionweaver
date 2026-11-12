"""Vendor the engine of the submitted IonWeaver notebook, byte for byte.

The notebook is one Colab cell carrying eighteen engine modules and two runners
as raw source strings (``_ENGINE['name'] = r'''...'''``).  This tool writes the
eighteen engine modules into ``engine/iw/`` without any change and records the
SHA-256 of each in ``engine/MANIFEST.json``.  The two runners (``run_single``,
``run_rest``) drive the study itself and are not needed by the web application.

    python tools/vendor_notebook.py --notebook "path/to/IonWeaver code.ipynb"
"""
from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DST = ROOT / "engine" / "iw"
ORDER = ("common", "inputs", "curation", "registry", "cores", "discovery", "fragments",
         "ledger", "operators", "build", "perception", "single", "strings", "sequential",
         "archive", "examples", "diversity", "facts")


def extract(notebook: Path) -> dict[str, str]:
    nb = json.loads(notebook.read_text(encoding="utf-8"))
    cell = "".join("".join(c["source"]) for c in nb["cells"] if c["cell_type"] == "code")
    out = {}
    for m in re.finditer(r"_ENGINE\['(\w+)'\] = r'''(.*?)'''\n", cell, re.S):
        out[m.group(1)] = m.group(2)
    missing = [n for n in ORDER if n not in out]
    if missing:
        raise SystemExit(f"modules not found in the notebook: {missing}")
    pins = {}
    for key in ("RDKIT_VERSION", "SELFIES_VERSION", "SEED", "ALLOWANCE", "SEQUENTIAL_ALLOWANCE"):
        m = re.search(r"^%s = (.+?)(?:\s+#.*)?$" % key, cell, re.M)
        pins[key] = m.group(1).strip().strip('"') if m else None
    return out, pins


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--notebook", type=Path, required=True)
    args = ap.parse_args()
    sources, pins = extract(args.notebook)
    DST.mkdir(parents=True, exist_ok=True)
    manifest = {
        "source_notebook": args.notebook.name,
        "source_notebook_sha256": hashlib.sha256(args.notebook.read_bytes()).hexdigest(),
        "vendored_utc": _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "pins": pins,
        "module_order": list(ORDER),
        "modules": {},
    }
    for name in ORDER:
        data = sources[name].encode("utf-8")
        (DST / f"{name}.py").write_bytes(data)
        manifest["modules"][name] = {"file": f"{name}.py", "sha256": hashlib.sha256(data).hexdigest(),
                                     "lines": sources[name].count("\n")}
        print(f"{name:12s} {len(data):>7,d} bytes  {manifest['modules'][name]['sha256'][:16]}")
    (DST / "__init__.py").write_text(
        '"""The IonWeaver engine, exactly as the submitted notebook carries it.\n\n'
        'Every module beside this file is the notebook\'s own source, written by\n'
        '``tools/vendor_notebook.py`` with no change; ``engine/MANIFEST.json`` records the\n'
        'SHA-256 of each and the server re-hashes them at start-up.\n"""\n',
        encoding="utf-8")
    (ROOT / "engine" / "MANIFEST.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print("manifest written; pins:", pins)


if __name__ == "__main__":
    main()
