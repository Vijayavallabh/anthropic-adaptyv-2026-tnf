#!/usr/bin/env python3
"""Stage 1 orchestrator: saturate the 8x H100 with BoltzGen design jobs.

Builds a job queue of (spec x seed) and runs one job per GPU concurrently, each via
scripts/gen_boltzgen.sh (which wires the isolated venv + weights cache). Jobs are pulled
from a shared queue so a GPU that finishes early immediately picks up the next job.
Success is judged by designs *produced* (parsed from the log), not exit code, because
`fast` mode may exit 3 (a perf-lever partial) while still writing all outputs.

Examples:
  # broad design+inverse_folding only (cheap; many candidates), 4 seeds/spec, 150/job
  python3 scripts/run_campaign.py --ndesigns 150 --seeds 4 --extra "-- --steps design inverse_folding"
  # full 6-step pipeline (fold+analysis+filter), fewer designs
  python3 scripts/run_campaign.py --ndesigns 80 --seeds 3
"""
from __future__ import annotations

import argparse
import json
import queue
import re
import shlex
import subprocess
import sys
import threading
import time
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
GEN = REPO / "scripts" / "gen_boltzgen.sh"


def now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def parse_gpus(s: str) -> list[int]:
    out: list[int] = []
    for part in s.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            a, b = part.split("-")
            out += list(range(int(a), int(b) + 1))
        else:
            out.append(int(part))
    return out


def parse_produced(logf: Path) -> int | None:
    prod = None
    try:
        for line in logf.read_text(errors="ignore").splitlines():
            m = re.search(r"DESIGNS step=design .*produced=(\d+)", line)
            if m:
                prod = int(m.group(1))
    except Exception:
        pass
    return prod


def run_job(gpu, job, mode, batch, extra, logdir, manifest, lock) -> None:
    spec, seed, ndes, out = job["spec"], job["seed"], job["ndes"], job["out"]
    name = Path(spec).stem
    logf = logdir / f"{name}_seed{seed}.log"
    cmd = ["bash", str(GEN), str(gpu), mode, spec, out, str(ndes), str(seed), str(batch)] + extra
    t0 = time.time()
    with open(logf, "w") as fh:
        fh.write(f"# {' '.join(cmd)}\n")
        fh.flush()
        rc = subprocess.call(cmd, stdout=fh, stderr=subprocess.STDOUT)
    dt = time.time() - t0
    produced = parse_produced(logf)
    rec = {
        "name": name, "spec": spec, "seed": seed, "gpu": gpu, "ndes": ndes,
        "rc": rc, "produced": produced, "secs": round(dt, 1), "out": out,
        "log": str(logf), "end": now(),
    }
    with lock:
        manifest.append(rec)
        (logdir / "manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"[{now()}] gpu{gpu} {name} seed{seed} rc={rc} produced={produced} {dt:.0f}s", flush=True)


def worker(gpu, position, stagger, q, mode, batch, extra, logdir, manifest, lock) -> None:
    if stagger:
        time.sleep(position * stagger)  # avoid concurrent first-run JIT races
    while True:
        try:
            job = q.get_nowait()
        except queue.Empty:
            return
        try:
            run_job(gpu, job, mode, batch, extra, logdir, manifest, lock)
        except Exception as e:  # keep the GPU pulling more jobs
            print(f"[{now()}] gpu{gpu} job error: {e}", flush=True)
        finally:
            q.task_done()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--specs", nargs="*", default=None, help="spec files (default specs/egfr_dom3_*.yaml)")
    ap.add_argument("--seeds", type=int, default=4, help="number of seeds per spec")
    ap.add_argument("--seed-start", type=int, default=0)
    ap.add_argument("--ndesigns", type=int, default=100)
    ap.add_argument("--batch", type=int, default=16)
    ap.add_argument("--mode", default="fast", choices=["fast", "exact", "big", "off"])
    ap.add_argument("--gpus", default="0-2")  # user ≤3-GPU cap (2026-10-02, default 0,1,2); pass --gpus to override
    ap.add_argument("--outroot", default=str(REPO / "runs" / "campaign1"))
    ap.add_argument("--logdir", default=str(REPO / "logs" / "campaign1"))
    ap.add_argument("--stagger", type=int, default=4, help="seconds*position delay before each GPU's first job")
    ap.add_argument("--extra", default="", help="extra args to boltzgen, e.g. '-- --steps design inverse_folding'")
    a = ap.parse_args()

    specs = a.specs or sorted(str(p) for p in (REPO / "specs").glob("egfr_dom3_*.yaml"))
    if not specs:
        sys.exit("no specs found (specs/egfr_dom3_*.yaml)")
    gpus = parse_gpus(a.gpus)
    extra = shlex.split(a.extra)
    outroot = Path(a.outroot).resolve()   # absolute: gen_boltzgen.sh cd's into the kit dir
    logdir = Path(a.logdir)
    outroot.mkdir(parents=True, exist_ok=True)
    logdir.mkdir(parents=True, exist_ok=True)

    q: queue.Queue = queue.Queue()
    njobs = 0
    for spec in specs:
        name = Path(spec).stem
        for s in range(a.seed_start, a.seed_start + a.seeds):
            out = str(outroot / name / f"seed{s}")
            q.put({"spec": spec, "seed": s, "ndes": a.ndesigns, "out": out})
            njobs += 1

    print(f"[{now()}] campaign start: {len(specs)} specs x {a.seeds} seeds = {njobs} jobs | "
          f"{a.ndesigns} designs/job (~{njobs * a.ndesigns} total) | gpus={gpus} mode={a.mode} "
          f"batch={a.batch} extra={extra!r}", flush=True)

    manifest: list = []
    lock = threading.Lock()
    threads = [
        threading.Thread(
            target=worker,
            args=(g, i, a.stagger, q, a.mode, a.batch, extra, logdir, manifest, lock),
            daemon=True,
        )
        for i, g in enumerate(gpus)
    ]
    t0 = time.time()
    for t in threads:
        t.start()
    q.join()
    tot = time.time() - t0

    ok = sum(1 for r in manifest if (r["produced"] or 0) > 0)
    total_prod = sum((r["produced"] or 0) for r in manifest)
    (logdir / "manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"[{now()}] DONE {ok}/{njobs} jobs produced; total designs produced={total_prod}; "
          f"wall={tot / 3600:.2f}h; manifest={logdir/'manifest.json'}", flush=True)


if __name__ == "__main__":
    main()
