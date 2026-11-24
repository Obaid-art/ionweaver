"""Run the IonWeaver engine server on Modal (modal.com), which has a free plan without a card.

    pip install modal
    modal token new                      # once: opens the browser to log in
    modal deploy deploy/modal_app.py     # from the repository root

The deploy prints the public address, of the form
https://<workspace>--ionweaver-web.modal.run ; put it in web/config.js.
The container starts on the first request after idling (about ten seconds,
the page waits for it) and shuts down after five idle minutes.
"""
from pathlib import Path

import modal

ROOT = Path(__file__).resolve().parents[1]
REMOTE = "/root/ionweaver"

image = (
    modal.Image.debian_slim(python_version="3.11")
    .apt_install("libxrender1", "libxext6", "libsm6", "libexpat1", "libfreetype6")
    .pip_install_from_requirements(str(ROOT / "requirements.txt"))
    .add_local_dir(str(ROOT / "engine"), remote_path=f"{REMOTE}/engine")
    .add_local_dir(str(ROOT / "data"), remote_path=f"{REMOTE}/data")
    .add_local_dir(str(ROOT / "web"), remote_path=f"{REMOTE}/web")
    .add_local_dir(str(ROOT / "server"), remote_path=f"{REMOTE}/server")
)

app = modal.App("ionweaver")


@app.function(image=image, memory=1024, scaledown_window=300, timeout=600)
@modal.concurrent(max_inputs=20)
@modal.asgi_app()
def web():
    import sys

    sys.path.insert(0, f"{REMOTE}/server")
    from app import app as fastapi_app   # loads the published library (about 3 s)

    return fastapi_app
