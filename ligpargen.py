#!/usr/bin/env python3

import argparse
import json
import logging
import time
from datetime import datetime, timezone
from pathlib import Path

import requests
from rdkit import Chem
from bs4 import BeautifulSoup
from urllib.parse import parse_qs, urlparse

# ==============================
# CONFIG
# ==============================
BASE_URL = "https://zarbi.chem.yale.edu/ligpargen"
SUBMIT_URL = BASE_URL + "/submit"
STATUS_URL = BASE_URL + "/status"   # may vary
DOWNLOAD_URL = BASE_URL + "/download"

DEFAULT_DELAY_SUBMIT = 10  # seconds between submissions
DEFAULT_DELAY_POLL = 20    # seconds between polling
DEFAULT_MAX_RETRIES = 5
DEFAULT_FORCEFIELD = "oplsaa"
DEFAULT_OUTPUT_FORMAT = "gromacs"

OUTPUT_DIR = Path("ligpargen_results")
STATE_FILE = Path("jobs_state.json")

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
logger = logging.getLogger(__name__)


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def normalize_job_entry(info: dict) -> dict:
    """Ensure persisted jobs include metadata keys (backward-compatible)."""
    info.setdefault("jobid", None)
    info.setdefault("status", "UNKNOWN")
    info.setdefault("last_checked", None)
    info.setdefault("attempts", 0)
    info.setdefault("error_message", None)
    return info


# ==============================
# UTILITIES
# ==============================

def _nonnegative_int(value: str) -> int:
    iv = int(value)
    if iv < 0:
        raise argparse.ArgumentTypeError("must be >= 0")
    return iv


def _positive_int(value: str) -> int:
    iv = int(value)
    if iv < 1:
        raise argparse.ArgumentTypeError("must be >= 1")
    return iv


def parse_args():
    parser = argparse.ArgumentParser(
        description="Submit PDB ligands to LigParGen and download parameter files.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable debug logging.",
    )
    parser.add_argument(
        "--delay-submit",
        type=_nonnegative_int,
        default=DEFAULT_DELAY_SUBMIT,
        metavar="SEC",
        help=f"Seconds to wait after each submit attempt (default: {DEFAULT_DELAY_SUBMIT}).",
    )
    parser.add_argument(
        "--delay-poll",
        type=_nonnegative_int,
        default=DEFAULT_DELAY_POLL,
        metavar="SEC",
        help=f"Seconds between status polling rounds (default: {DEFAULT_DELAY_POLL}).",
    )
    parser.add_argument(
        "--max-retries",
        type=_positive_int,
        default=DEFAULT_MAX_RETRIES,
        help=f"Submit retries when no job id is returned (default: {DEFAULT_MAX_RETRIES}).",
    )
    parser.add_argument(
        "--forcefield",
        default=DEFAULT_FORCEFIELD,
        help=f"Force field passed to LigParGen (default: {DEFAULT_FORCEFIELD}).",
    )
    parser.add_argument(
        "--output-format",
        default=DEFAULT_OUTPUT_FORMAT,
        dest="output_format",
        help=(
            "Output format for LigParGen form field 'output' "
            f"(default: {DEFAULT_OUTPUT_FORMAT})."
        ),
    )
    return parser.parse_args()


def load_state():
    if STATE_FILE.exists():
        with STATE_FILE.open(encoding="utf-8") as f:
            raw = json.load(f)
        return {name: normalize_job_entry(dict(entry)) for name, entry in raw.items()}
    return {}

def save_state(state):
    with STATE_FILE.open("w", encoding="utf-8") as f:
        json.dump(state, f, indent=2)

def pdb_to_smiles(pdb_path):
    mol = Chem.MolFromPDBFile(str(pdb_path), removeHs=False)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol)


def extract_jobid_from_candidate(candidate):
    if not candidate:
        return None

    parsed = urlparse(candidate)
    query = parse_qs(parsed.query)
    if "jobid" in query and query["jobid"]:
        return query["jobid"][0]

    marker = "jobid="
    lower_candidate = candidate.lower()
    idx = lower_candidate.find(marker)
    if idx == -1:
        return None

    start = idx + len(marker)
    allowed = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-"
    jobid_chars = []
    for ch in candidate[start:]:
        if ch in allowed:
            jobid_chars.append(ch)
        else:
            break

    return "".join(jobid_chars) if jobid_chars else None

# ==============================
# SUBMIT JOB
# ==============================

def submit_job(session, smiles, name, *, forcefield, output_format, max_retries):
    data = {
        "smiData": smiles,
        "forcefield": forcefield,
        "output": output_format,
    }

    for attempt in range(max_retries):
        try:
            r = session.post(SUBMIT_URL, data=data, timeout=60)
            if r.status_code == 200:
                text = r.text
                soup = BeautifulSoup(text, "html.parser")

                # Parse links/forms first, e.g. href/action containing ?jobid=...
                for tag in soup.find_all(["a", "form"]):
                    candidate = tag.get("href") or tag.get("action") or ""
                    jobid = extract_jobid_from_candidate(candidate)
                    if jobid:
                        return jobid

                # Also handle hidden form inputs, e.g. <input name="jobid" value="...">.
                for input_tag in soup.find_all("input"):
                    if (input_tag.get("name") or "").lower() == "jobid":
                        value = input_tag.get("value", "")
                        if value:
                            return value

                # Fallback for pages that only contain plain text jobid.
                jobid = extract_jobid_from_candidate(text)
                if jobid:
                    return jobid
            time.sleep(5)
        except requests.RequestException as e:
            logger.warning("Submit failed %s: %s", name, e)
            time.sleep(5)
        except (UnicodeDecodeError, AttributeError, TypeError) as e:
            logger.warning("Submit response handling failed %s: %s", name, e)
            time.sleep(5)

    return None

# ==============================
# CHECK STATUS
# ==============================

def check_status(session, jobid):
    """
    Return (status, error_detail).

    error_detail is set for UNKNOWN (network/decode/non-200) or FAILED (server-side error page).
    """
    try:
        r = session.get(f"{STATUS_URL}?jobid={jobid}", timeout=30)
        if r.status_code == 200:
            text = r.text.lower()
            if "finished" in text:
                return "DONE", None
            if "error" in text:
                return "FAILED", "Status page indicates error"
            return "RUNNING", None
        err = f"HTTP {r.status_code}"
        logger.warning("Status check non-OK for job %s: %s", jobid, err)
        return "UNKNOWN", err
    except requests.RequestException as e:
        logger.warning("Status check failed for job %s: %s", jobid, e)
        return "UNKNOWN", str(e)
    except UnicodeDecodeError as e:
        logger.warning("Status response decode failed for job %s: %s", jobid, e)
        return "UNKNOWN", str(e)

# ==============================
# DOWNLOAD RESULTS
# ==============================

def download_results(session, jobid, name):
    outdir = OUTPUT_DIR / name
    outdir.mkdir(parents=True, exist_ok=True)

    files = ["ligand.itp", "ligand.gro", "ligand.top"]

    for fname in files:
        url = f"{DOWNLOAD_URL}/{jobid}/{fname}"
        try:
            r = session.get(url, timeout=60)
            if r.status_code == 200:
                (outdir / fname).write_bytes(r.content)
        except Exception as e:
            logger.warning("Download failed %s: %s", fname, e)

# ==============================
# MAIN PIPELINE
# ==============================

def main():
    args = parse_args()
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    delay_submit = args.delay_submit
    delay_poll = args.delay_poll
    max_retries = args.max_retries
    forcefield = args.forcefield
    output_format = args.output_format

    pdb_files = sorted(Path(".").glob("*.pdb"))
    state = load_state()
    session = requests.Session()

    # ---- Step 1: Submit jobs ----
    for pdb_path in pdb_files:
        name = pdb_path.stem

        if name in state:
            continue

        smiles = pdb_to_smiles(pdb_path)
        if not smiles:
            logger.error("Failed SMILES: %s", pdb_path)
            state[name] = normalize_job_entry(
                {
                    "jobid": None,
                    "status": "SMILES_FAILED",
                    "last_checked": utc_now_iso(),
                    "attempts": 0,
                    "error_message": "RDKit could not produce SMILES from PDB",
                }
            )
            save_state(state)
            continue

        logger.info(
            "Submitting %s (forcefield=%s, output=%s)",
            name,
            forcefield,
            output_format,
        )
        jobid = submit_job(
            session,
            smiles,
            name,
            forcefield=forcefield,
            output_format=output_format,
            max_retries=max_retries,
        )

        if jobid:
            state[name] = normalize_job_entry(
                {
                    "jobid": jobid,
                    "status": "SUBMITTED",
                    "last_checked": utc_now_iso(),
                    "attempts": 0,
                    "error_message": None,
                }
            )
            save_state(state)
        else:
            logger.error("Submission failed: %s", name)
            state[name] = normalize_job_entry(
                {
                    "jobid": None,
                    "status": "SUBMIT_FAILED",
                    "last_checked": utc_now_iso(),
                    "attempts": max_retries,
                    "error_message": "No job id returned after submit retries",
                }
            )
            save_state(state)

        time.sleep(delay_submit)

    # ---- Step 2: Poll + Download ----
    logger.info("Polling jobs...")

    unfinished = True
    while unfinished:
        unfinished = False

        for name, info in state.items():
            info = normalize_job_entry(info)
            state[name] = info

            if info["status"] == "DONE":
                continue

            jobid = info["jobid"]
            if not jobid:
                continue

            info["attempts"] = int(info.get("attempts") or 0) + 1
            info["last_checked"] = utc_now_iso()

            status, status_error = check_status(session, jobid)

            logger.info("%s: %s", name, status)

            if status == "DONE":
                info["error_message"] = None
                download_results(session, jobid, name)
                info["status"] = "DONE"
                save_state(state)

            elif status == "RUNNING":
                info["error_message"] = None
                unfinished = True
                save_state(state)

            elif status == "UNKNOWN":
                info["error_message"] = status_error
                unfinished = True
                save_state(state)

            elif status == "FAILED":
                info["status"] = "FAILED"
                info["error_message"] = status_error or "Job failed"
                save_state(state)

        if unfinished:
            time.sleep(delay_poll)

    logger.info("All jobs processed.")

if __name__ == "__main__":
    main()
