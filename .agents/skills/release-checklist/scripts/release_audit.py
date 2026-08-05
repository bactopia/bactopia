#!/usr/bin/env python3
"""Deterministic release-readiness checks for Bactopia (audit-only).

Covers the git + file-state checks that are NOT version-file consistency:
sibling-repo release state, CHANGELOG structure, and the /update-module record.
It never writes to any repo and never touches git history -- `git` is invoked
read-only (describe/tag/rev-list/status/log).

Version-bearing file consistency (versions.yml vs nextflow.config, catalog.json,
CITATION.cff, bin/bactopia, data/conda/meta.yaml, every *.config bactopia_version
and nf-bactopia@ pin, versions.yml vs CHANGELOG, and the declared nf-bactopia pin
vs the nf-bactopia repo's latest) is owned by the `bactopia-lint` repo rules
V001/V002/V003 -- the release-checklist skill surfaces those from the lint run
rather than re-implementing them here.

Companion to the `release-checklist` skill: the SKILL.md orchestrates the
CLI-backed checks (bactopia-lint/update/merge-schemas/catalog/docs/citations/
review-tests) and gh; this script owns the remaining pure file+git checks.
"""

import argparse
import json
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

_SECTION_RE = re.compile(r"^##\s+v?(\d+\.\d+\.\d+)\b")


# --------------------------------------------------------------------------- #
# small helpers
# --------------------------------------------------------------------------- #
def _read(path: Path) -> str | None:
    try:
        return path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return None


def _git(repo: Path, *args: str) -> str | None:
    """Run a read-only git command in `repo`; return stripped stdout or None."""
    try:
        out = subprocess.run(
            ["git", "-C", str(repo), *args],
            capture_output=True, text=True, timeout=30,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    if out.returncode != 0:
        return None
    return out.stdout.strip()


def _changelog_top(text: str | None) -> dict:
    """Parse the top `## v<x.y.z> ... "codename" date` heading of a CHANGELOG."""
    result = {"version": None, "codename": None, "date": None, "placeholder": False}
    if not text:
        return result
    for line in text.splitlines():
        m = _SECTION_RE.match(line)
        if not m:
            continue
        result["version"] = m.group(1)
        cm = re.search(r'"([^"]*)"', line)
        if cm:
            result["codename"] = cm.group(1)
        dm = re.search(r'"[^"]*"\s*[-\s]*(.+?)\s*$', line)
        if dm:
            result["date"] = dm.group(1).strip()
        cn = result["codename"]
        dt = result["date"] or ""
        result["placeholder"] = (cn is not None and "?" in cn) or ("?" in dt)
        break
    return result


def _all_sections(text: str | None) -> list[str]:
    if not text:
        return []
    return [m.group(1) for line in text.splitlines()
            if (m := _SECTION_RE.match(line))]


def _read_versions_yml(bp: Path) -> dict:
    """Parse the declared versions.yml (simple `key: value`; stdlib only)."""
    text = _read(bp / "versions.yml")
    out = {"present": text is not None, "bactopia": None, "nf_bactopia": None}
    if not text:
        return out
    for line in text.splitlines():
        m = re.match(r"\s*([A-Za-z0-9_-]+)\s*:\s*(\S+)", line)
        if not m:
            continue
        key, val = m.group(1), m.group(2).strip().strip("'\"")
        if key == "bactopia":
            out["bactopia"] = val
        elif key in ("nf-bactopia", "nf_bactopia"):
            out["nf_bactopia"] = val
    return out


# --------------------------------------------------------------------------- #
# check families
# --------------------------------------------------------------------------- #
def sibling_release_state(repos: dict[str, Path | None], warnings: list) -> dict:
    out = {}
    for name, repo in repos.items():
        if not repo or not repo.exists():
            out[name] = None
            warnings.append(f"repo '{name}' not found; release-state skipped")
            continue
        latest_tag = _git(repo, "describe", "--tags", "--abbrev=0")
        commits_ahead = None
        if latest_tag:
            cnt = _git(repo, "rev-list", f"{latest_tag}..HEAD", "--count")
            commits_ahead = int(cnt) if cnt and cnt.isdigit() else None
        status = _git(repo, "status", "--porcelain")
        dirty = len([ln for ln in status.splitlines() if ln.strip()]) if status else 0
        changelog = _changelog_top(_read(repo / "CHANGELOG.md"))
        out[name] = {
            "latest_tag": latest_tag,
            "commits_ahead": commits_ahead,
            "dirty_files": dirty,
            "changelog_top": changelog["version"],
            "needs_release": bool(commits_ahead and commits_ahead > 0),
        }
    return out


def changelog_state(bp: Path) -> dict:
    text = _read(bp / "CHANGELOG.md")
    sections = _all_sections(text)
    top = _changelog_top(text)
    latest_tag = _git(bp, "describe", "--tags", "--abbrev=0")
    commits_since_tag = None
    if latest_tag:
        log = _git(bp, "log", f"{latest_tag}..HEAD", "--oneline")
        commits_since_tag = len(log.splitlines()) if log else 0
    # bullets in the top section (until the next `## ` heading)
    top_bullets = 0
    if text:
        in_top = False
        for line in text.splitlines():
            if _SECTION_RE.match(line):
                if in_top:
                    break
                in_top = True
                continue
            if in_top and re.match(r"^\s*-\s+", line):
                top_bullets += 1
    # release target = declared version (versions.yml); it must have a section
    target = _read_versions_yml(bp).get("bactopia")
    missing = bool(target and target not in sections)
    return {
        "target_version": target,
        "top_version": top["version"],
        "codename": top["codename"],
        "date": top["date"],
        "placeholder": top["placeholder"],
        "latest_tag": latest_tag,
        "sections": sections,
        "commits_since_tag": commits_since_tag,
        "top_section_bullets": top_bullets,
        "missing_section_for_target": missing,
    }


def module_updates(bp: Path) -> dict:
    """Read the newest /update-module record under logs/module-updates/.

    The update-module skill writes bactopia-update's JSON there after each run.
    This is offline -- it only reads the record, never queries the network.
    Returns present=False when no record exists (a release blocker).
    """
    log_dir = bp / "logs" / "module-updates"
    files = sorted(log_dir.glob("*.json")) if log_dir.is_dir() else []
    if not files:
        return {"present": False, "log_dir": str(log_dir)}
    newest = files[-1]  # timestamped filenames sort chronologically
    m = re.match(r"(\d{8}_\d{6})", newest.stem)
    log_ts = None
    if m:
        try:
            log_ts = int(datetime.strptime(m.group(1), "%Y%m%d_%H%M%S")
                         .replace(tzinfo=timezone.utc).timestamp())
        except ValueError:
            log_ts = None
    if log_ts is None:
        try:
            log_ts = int(newest.stat().st_mtime)
        except OSError:
            log_ts = None

    needs_update = needs_review = up_to_date = total = None
    try:
        data = json.loads(newest.read_text())
        rows = data.get("results", data) if isinstance(data, dict) else data
        if isinstance(rows, dict):
            rows = rows.get("modules") or rows.get("results") or []
        if isinstance(rows, list):
            rows = [x for x in rows if isinstance(x, dict)
                    and not (x.get("tool") == "aria2"
                             and x.get("module") == "checkm2_download")]
            total = len(rows)
            needs_update = sum(1 for x in rows if x.get("needs_update"))
            needs_review = sum(1 for x in rows if x.get("needs_user_review"))
            up_to_date = total - needs_update - needs_review
    except (OSError, json.JSONDecodeError):
        pass

    cfg_ts = _git(bp, "log", "-1", "--format=%ct",
                  "--", ":(glob)modules/**/module.config")
    cfg_ts = int(cfg_ts) if cfg_ts and cfg_ts.isdigit() else None
    return {
        "present": True,
        "log_file": str(newest),
        "log_timestamp": m.group(1) if m else newest.stem,
        "needs_update": needs_update,
        "needs_user_review": needs_review,
        "up_to_date": up_to_date,
        "total": total,
        "module_config_changed_after": bool(log_ts and cfg_ts and cfg_ts > log_ts),
    }


# --------------------------------------------------------------------------- #
# entry point
# --------------------------------------------------------------------------- #
def _resolve(explicit: str | None, default: Path) -> Path | None:
    path = Path(explicit) if explicit else default
    return path if path.exists() else None


def main() -> int:
    ap = argparse.ArgumentParser(description="Deterministic Bactopia release audit.")
    ap.add_argument("--bactopia-path", required=True)
    ap.add_argument("--bactopia-py-path")
    ap.add_argument("--nf-bactopia-path")
    ap.add_argument("--docs-site-path")
    ap.add_argument("--json", action="store_true",
                    help="emit JSON (default; kept for parity with other CLIs)")
    args = ap.parse_args()

    bp = Path(args.bactopia_path).resolve()
    if not bp.exists():
        print(json.dumps({"error": f"bactopia-path not found: {bp}"}))
        return 0
    parent = bp.parent
    py = _resolve(args.bactopia_py_path, parent / "bactopia-py")
    nfb = _resolve(args.nf_bactopia_path, parent / "nf-bactopia")
    docs = _resolve(args.docs_site_path, parent / "bactopia.github.io")

    warnings: list[str] = []
    for name, path in (("bactopia-py", py), ("nf-bactopia", nfb),
                       ("bactopia.github.io", docs)):
        if path is None:
            warnings.append(f"sibling repo '{name}' not found near {parent}")

    siblings = sibling_release_state(
        {"bactopia": bp, "bactopia-py": py, "nf-bactopia": nfb}, warnings)
    cl = changelog_state(bp)
    mu = module_updates(bp)

    report = {
        "bactopia_path": str(bp),
        "paths": {
            "bactopia_py": str(py) if py else None,
            "nf_bactopia": str(nfb) if nfb else None,
            "docs_site": str(docs) if docs else None,
        },
        "sibling_release_state": siblings,
        "changelog": cl,
        "module_updates": mu,
        "warnings": warnings,
        "generated": datetime.now(timezone.utc).isoformat(),
    }
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
