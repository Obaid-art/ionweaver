"""Work that takes more than a moment runs in a thread and reports progress."""
from __future__ import annotations

import logging
import threading
import time
import traceback
import uuid
from typing import Any, Callable, Dict, Optional

log = logging.getLogger("ionweaver.jobs")
#: how many jobs may run at the same time; the rest wait in order
RUNNING = threading.Semaphore(3)

__all__ = ["Job", "store"]


class Job:
    def __init__(self, kind: str, label: str) -> None:
        self.id = uuid.uuid4().hex[:10]
        self.kind, self.label = kind, label
        self.state = "queued"
        self.progress = 0.0
        self.message = "queued"
        self.result: Any = None
        self.error: str = ""
        self.created = time.time()
        self.finished: Optional[float] = None

    def status(self) -> Dict[str, Any]:
        return {"job_id": self.id, "kind": self.kind, "label": self.label, "state": self.state,
                "progress": round(self.progress, 3), "message": self.message, "error": self.error,
                "elapsed_s": round((self.finished or time.time()) - self.created, 1)}


class Store:
    def __init__(self) -> None:
        self._jobs: Dict[str, Job] = {}
        self._lock = threading.Lock()

    def submit(self, kind: str, label: str, work: Callable[[Callable[[float, str], None]], Any]) -> Job:
        job = Job(kind, label)
        with self._lock:
            self._jobs[job.id] = job
            self._sweep()

        def progress(fraction: float, message: str) -> None:
            job.progress = max(0.0, min(1.0, float(fraction)))
            job.message = str(message)

        def run() -> None:
            job.message = "waiting for a free slot"
            with RUNNING:
                job.state = "running"
                job.message = "starting"
                try:
                    job.result = work(progress)
                    job.state = "done"
                    job.progress = 1.0
                    job.message = "done"
                except Exception as exc:  # noqa: BLE001 - reported to the user, never hidden
                    job.state = "failed"
                    job.error = str(exc) or exc.__class__.__name__
                    job.message = job.error
                    log.error("job %s (%s) failed:\n%s", job.id, job.kind, traceback.format_exc())
                finally:
                    job.finished = time.time()

        threading.Thread(target=run, daemon=True, name=f"job-{job.id}").start()
        return job

    def get(self, job_id: str) -> Optional[Job]:
        return self._jobs.get(str(job_id))

    def _sweep(self, keep: int = 60) -> None:
        if len(self._jobs) <= keep:
            return
        done = sorted((j for j in self._jobs.values() if j.finished), key=lambda j: j.finished)
        for j in done[: len(self._jobs) - keep]:
            self._jobs.pop(j.id, None)


store = Store()
