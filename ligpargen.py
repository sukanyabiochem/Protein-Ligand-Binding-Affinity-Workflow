"""
Automated script to submit a structure file to LigParGen and download
the generated OPLS-AA parameter files (.prm, .rtf, .itp, .gro, .pdb, .xml, .zip, etc.).

Server: https://zarbi.chem.yale.edu/ligpargen/
Accepted inputs: PDB or MOL (all hydrogens must be present).
"""

from __future__ import annotations

import argparse
import logging
import os
import platform
import ssl
import sys
import time
from pathlib import Path
from typing import Optional
from urllib.parse import urljoin

import certifi
import requests
from bs4 import BeautifulSoup
from bs4.element import Tag
from requests.adapters import HTTPAdapter

LIGPARGEN_URL = "https://zarbi.chem.yale.edu/ligpargen/"
SUBMIT_URL = "https://zarbi.chem.yale.edu/cgi-bin/results_lpg.py"
UPLOAD_FIELD = "molpdbfile"

# Yale's TLS often omits this intermediate; certifi has the GlobalSign root but not this CA.
_GLOBALSIGN_RSA_OV_SSL_CA_2018_PEM = """-----BEGIN CERTIFICATE-----
MIIETjCCAzagAwIBAgINAe5fIh38YjvUMzqFVzANBgkqhkiG9w0BAQsFADBMMSAw
HgYDVQQLExdHbG9iYWxTaWduIFJvb3QgQ0EgLSBSMzETMBEGA1UEChMKR2xvYmFs
U2lnbjETMBEGA1UEAxMKR2xvYmFsU2lnbjAeFw0xODExMjEwMDAwMDBaFw0yODEx
MjEwMDAwMDBaMFAxCzAJBgNVBAYTAkJFMRkwFwYDVQQKExBHbG9iYWxTaWduIG52
LXNhMSYwJAYDVQQDEx1HbG9iYWxTaWduIFJTQSBPViBTU0wgQ0EgMjAxODCCASIw
DQYJKoZIhvcNAQEBBQADggEPADCCAQoCggEBAKdaydUMGCEAI9WXD+uu3Vxoa2uP
UGATeoHLl+6OimGUSyZ59gSnKvuk2la77qCk8HuKf1UfR5NhDW5xUTolJAgvjOH3
idaSz6+zpz8w7bXfIa7+9UQX/dhj2S/TgVprX9NHsKzyqzskeU8fxy7quRU6fBhM
abO1IFkJXinDY+YuRluqlJBJDrnw9UqhCS98NE3QvADFBlV5Bs6i0BDxSEPouVq1
lVW9MdIbPYa+oewNEtssmSStR8JvA+Z6cLVwzM0nLKWMjsIYPJLJLnNvBhBWk0Cq
o8VS++XFBdZpaFwGue5RieGKDkFNm5KQConpFmvv73W+eka440eKHRwup08CAwEA
AaOCASkwggElMA4GA1UdDwEB/wQEAwIBhjASBgNVHRMBAf8ECDAGAQH/AgEAMB0G
A1UdDgQWBBT473/yzXhnqN5vjySNiPGHAwKz6zAfBgNVHSMEGDAWgBSP8Et/qC5F
JK5NUPpjmove4t0bvDA+BggrBgEFBQcBAQQyMDAwLgYIKwYBBQUHMAGGImh0dHA6
Ly9vY3NwMi5nbG9iYWxzaWduLmNvbS9yb290cjMwNgYDVR0fBC8wLTAroCmgJ4Yl
aHR0cDovL2NybC5nbG9iYWxzaWduLmNvbS9yb290LXIzLmNybDBHBgNVHSAEQDA+
MDwGBFUdIAAwNDAyBggrBgEFBQcCARYmaHR0cHM6Ly93d3cuZ2xvYmFsc2lnbi5j
b20vcmVwb3NpdG9yeS8wDQYJKoZIhvcNAQELBQADggEBAJmQyC1fQorUC2bbmANz
EdSIhlIoU4r7rd/9c446ZwTbw1MUcBQJfMPg+NccmBqixD7b6QDjynCy8SIwIVbb
0615XoFYC20UgDX1b10d65pHBf9ZjQCxQNqQmJYaumxtf4z1s4DfjGRzNpZ5eWl0
6r/4ngGPoJVpjemEuunl1Ig423g7mNA2eymw0lIYkN5SQwCuaifIFJ6GlazhgDEw
fpolu4usBCOmmQDo8dIm7A9+O4orkjgTHY+GzYZSR+Y0fFukAj6KYXwidlNalFMz
hriSqHKvoflShx8xpfywgVcvzfTO3PYkz6fiNJBonf6q8amaEsybwMbDqKWwIX7e
SPY=
-----END CERTIFICATE-----
"""

DEFAULT_OUTPUT_DIR = "ligpargen_output"
POLL_INTERVAL_SECONDS = 10
MAX_POLL_ATTEMPTS = 60  # ~10 minutes max wait

DEFAULT_SUFFIXES = frozenset(
    {
        ".itp",
        ".top",
        ".gro",
        ".prm",
        ".rtf",
        ".str",
        ".inp",
        ".xml",
        ".pdb",
        ".mol2",
        ".xyz",
        ".zip",
        ".pqr",
        ".psf",
        ".z",
    }
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
logger = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Submit a PDB/MOL file to LigParGen and download OPLS-AA parameter files.",
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Path to the input PDB or MOL file (must contain all hydrogen atoms).",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        default=DEFAULT_OUTPUT_DIR,
        help=f"Directory to save downloaded files (default: {DEFAULT_OUTPUT_DIR}).",
    )
    parser.add_argument(
        "--no-download",
        action="store_true",
        help="Only print download targets without saving files.",
    )
    parser.add_argument(
        "--wait",
        type=int,
        default=POLL_INTERVAL_SECONDS,
        help=f"Seconds to wait between polling attempts (default: {POLL_INTERVAL_SECONDS}).",
    )
    parser.add_argument(
        "--max-polls",
        type=int,
        default=MAX_POLL_ATTEMPTS,
        help=f"Maximum polling attempts when results are not ready (default: {MAX_POLL_ATTEMPTS}).",
    )
    parser.add_argument(
        "--opt",
        type=int,
        choices=(0, 1, 2, 3),
        default=0,
        help="Molecule optimization iterations (default: 0).",
    )
    parser.add_argument(
        "--charge-model",
        choices=("cm1abcc", "cm1a"),
        default="cm1abcc",
        help="Charge model: cm1abcc (neutral) or cm1a (neutral/charged) (default: cm1abcc).",
    )
    parser.add_argument(
        "--charge",
        type=int,
        choices=(-2, -1, 0, 1, 2),
        default=0,
        help="Net molecule charge when using cm1a (default: 0).",
    )
    parser.add_argument(
        "--only",
        default=None,
        metavar="EXTS",
        help="Comma-separated extensions to download, without dots (e.g. prm,rtf,itp).",
    )
    parser.add_argument(
        "--insecure",
        action="store_true",
        help="Disable TLS certificate verification (not recommended).",
    )
    parser.add_argument(
        "--ca-bundle",
        default=None,
        metavar="PEM",
        help=(
            "Custom CA bundle PEM file. For conda, try: "
            f"{Path(sys.prefix) / 'ssl' / 'cacert.pem'}"
        ),
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable verbose/debug logging.",
    )
    return parser.parse_args()


def _extra_ca_paths() -> list[Path]:
    paths = [Path(sys.prefix) / "ssl" / "cacert.pem"]
    if platform.system() == "Darwin":
        paths.extend(
            [
                Path("/etc/ssl/cert.pem"),
                Path("/opt/homebrew/etc/openssl@3/cert.pem"),
                Path("/usr/local/etc/openssl@3/cert.pem"),
            ]
        )
    return paths


def build_trusted_ssl_context(*, ca_bundle: Optional[str]) -> ssl.SSLContext:
    ctx = ssl.SSLContext(ssl.PROTOCOL_TLS_CLIENT)
    ctx.check_hostname = True
    ctx.verify_mode = ssl.CERT_REQUIRED
    if ca_bundle:
        path = Path(ca_bundle).expanduser()
        if not path.is_file():
            raise FileNotFoundError(f"--ca-bundle not found: {path}")
        ctx.load_verify_locations(cafile=str(path))
        return ctx
    ctx.load_verify_locations(cafile=certifi.where())
    try:
        ctx.load_verify_locations(cadata=_GLOBALSIGN_RSA_OV_SSL_CA_2018_PEM)
    except ssl.SSLError as exc:
        logger.debug("Could not load bundled GlobalSign intermediate: %s", exc)
    for path in _extra_ca_paths():
        if path.is_file():
            try:
                ctx.load_verify_locations(cafile=str(path))
            except ssl.SSLError as exc:
                logger.debug("Could not merge CA file %s: %s", path, exc)
    return ctx


class _TLSAdapter(HTTPAdapter):
    def __init__(self, ssl_context: ssl.SSLContext, **kwargs):
        self._ssl_context = ssl_context
        kwargs.setdefault("pool_connections", 1)
        kwargs.setdefault("pool_maxsize", 2)
        super().__init__(**kwargs)

    def init_poolmanager(self, connections, maxsize, block=False, **pool_kwargs):
        pool_kwargs.setdefault("ssl_context", self._ssl_context)
        return super().init_poolmanager(connections, maxsize, block=block, **pool_kwargs)


def apply_ssl_env_defaults() -> None:
    conda_ca = Path(sys.prefix) / "ssl" / "cacert.pem"
    if conda_ca.is_file():
        os.environ.setdefault("SSL_CERT_FILE", str(conda_ca))
        os.environ.setdefault("REQUESTS_CA_BUNDLE", str(conda_ca))
    os.environ.setdefault("CURL_CA_BUNDLE", certifi.where())


def new_requests_session(*, ssl_insecure: bool, ca_bundle: Optional[str]) -> requests.Session:
    session = requests.Session()
    if ssl_insecure:
        session.verify = False
        logger.warning("TLS verification disabled (--insecure).")
    else:
        ctx = build_trusted_ssl_context(ca_bundle=ca_bundle)
        session.mount("https://", _TLSAdapter(ctx))
        session.verify = True
    session.headers.update(
        {
            "User-Agent": (
                "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 "
                "(KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36"
            ),
            "Accept-Encoding": "identity",
        }
    )
    return session


def parse_only_suffixes(spec: Optional[str]) -> frozenset[str]:
    if not spec or not spec.strip():
        return DEFAULT_SUFFIXES
    out: set[str] = set()
    for part in spec.split(","):
        p = part.strip().lower().lstrip(".")
        if p:
            out.add("." + p)
    return frozenset(out) if out else DEFAULT_SUFFIXES


def validate_input_file(filepath: str) -> Path:
    path = Path(filepath)
    if not path.exists():
        logger.error("Input file not found: %s", filepath)
        sys.exit(1)
    ext = path.suffix.lower()
    if ext not in {".pdb", ".mol"}:
        if ext == ".mol2":
            logger.error(
                "LigParGen does not accept MOL2 directly. Convert to PDB or MOL first "
                "(e.g. with Open Babel: obabel -imol2 %s -opdb -O out.pdb -h).",
                path.name,
            )
        else:
            logger.warning(
                "Input extension %s may not be accepted; use .pdb or .mol.", ext
            )
        sys.exit(1)
    return path.resolve()


def _charge_form_value(charge: int) -> str:
    if charge > 0:
        return f" +{charge} "
    if charge < 0:
        return f" {charge} "
    return " 0 "


def _build_submit_data(args: argparse.Namespace) -> dict[str, str]:
    return {
        "smiData": "",
        "checkopt": f" {args.opt} ",
        "chargetype": args.charge_model,
        "dropcharge": _charge_form_value(args.charge),
    }


def _mime_type_for(path: Path) -> str:
    if path.suffix.lower() == ".pdb":
        return "chemical/x-pdb"
    return "chemical/x-mdl-molfile"


def submit_file(
    session: requests.Session,
    input_path: Path,
    args: argparse.Namespace,
) -> tuple[str, str]:
    """Submit the structure to LigParGen. Returns (results_page_html, base_url)."""
    logger.info("Opening LigParGen at %s", LIGPARGEN_URL)
    session.get(LIGPARGEN_URL, timeout=120).raise_for_status()

    data = _build_submit_data(args)
    with open(input_path, "rb") as handle:
        files = {UPLOAD_FIELD: (input_path.name, handle, _mime_type_for(input_path))}
        logger.info("Submitting %s...", input_path.name)
        response = session.post(SUBMIT_URL, data=data, files=files, timeout=180)
    response.raise_for_status()

    html = response.text
    base_url = response.url or SUBMIT_URL

    if args.verbose:
        debug_path = Path("ligpargen_response.html")
        debug_path.write_text(html, encoding="utf-8")
        logger.debug("Saved response HTML to %s", debug_path)

    return html, base_url


def _suffix_allowed(fileout: str, allowed_suffixes: frozenset[str]) -> bool:
    name = Path(fileout).name.lower()
    for sfx in allowed_suffixes:
        ext = sfx.lstrip(".")
        if ext == "prm" and name.endswith(".q.prm"):
            continue
        if name.endswith(sfx):
            return True
    return False


def collect_download_targets(
    html: str,
    base_url: str,
    allowed_suffixes: frozenset[str],
) -> list[tuple[str, dict[str, str], str]]:
    soup = BeautifulSoup(html, "html.parser")
    targets: list[tuple[str, dict[str, str], str]] = []
    seen: set[str] = set()

    for form in soup.find_all("form"):
        if not isinstance(form, Tag):
            continue
        action = (form.get("action") or "").strip()
        if "download_lpg" not in action.lower():
            continue
        if (form.get("method") or "post").upper() != "POST":
            continue

        data: dict[str, str] = {}
        fileout: Optional[str] = None
        for inp in form.find_all("input"):
            if not isinstance(inp, Tag):
                continue
            name = inp.get("name")
            if not name:
                continue
            value = inp.get("value") or ""
            data[name] = value
            if name == "fileout":
                fileout = value

        if not fileout or not _suffix_allowed(fileout, allowed_suffixes):
            continue

        filename = Path(fileout).name
        if filename in seen:
            continue
        seen.add(filename)

        post_url = urljoin(base_url, action)
        targets.append((post_url, data, filename))

    return targets


def poll_until_ready(
    session: requests.Session,
    html: str,
    base_url: str,
    allowed_suffixes: frozenset[str],
    wait_s: int,
    max_polls: int,
) -> tuple[str, str, list[tuple[str, dict[str, str], str]]]:
    targets = collect_download_targets(html, base_url, allowed_suffixes)
    if targets:
        return html, base_url, targets

    lower = html.lower()
    if "error has been detected" in lower or "errorcode" in lower:
        logger.error(
            "LigParGen reported an input/validation error. Check your structure and charge options."
        )
        sys.exit(1)

    for attempt in range(1, max_polls + 1):
        logger.info(
            "No download forms yet; waiting %s s (%s/%s)...",
            wait_s,
            attempt,
            max_polls,
        )
        time.sleep(wait_s)
        resp = session.get(base_url, timeout=120)
        resp.raise_for_status()
        html = resp.text
        base_url = resp.url or base_url
        targets = collect_download_targets(html, base_url, allowed_suffixes)
        if targets:
            return html, base_url, targets

    return html, base_url, targets


def download_files(
    session: requests.Session,
    targets: list[tuple[str, dict[str, str], str]],
    output_dir: Path,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("Downloading %d file(s) to %s", len(targets), output_dir)

    for post_url, data, filename in targets:
        dest = output_dir / filename
        logger.info("  Downloading: %s", filename)
        resp = session.post(post_url, data=data, timeout=180)
        resp.raise_for_status()
        dest.write_bytes(resp.content)

    logger.info("All files downloaded successfully.")


def main() -> None:
    args = parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    apply_ssl_env_defaults()
    input_path = validate_input_file(args.input)
    allowed_suffixes = parse_only_suffixes(args.only)

    try:
        session = new_requests_session(ssl_insecure=args.insecure, ca_bundle=args.ca_bundle)
    except FileNotFoundError as exc:
        logger.error("%s", exc)
        sys.exit(1)

    try:
        html, base_url = submit_file(session, input_path, args)
    except requests.exceptions.SSLError as exc:
        conda_pem = Path(sys.prefix) / "ssl" / "cacert.pem"
        logger.error(
            "TLS verification failed: %s\n"
            "Try: pip install -U certifi\n"
            "Or point at your conda CA bundle:\n"
            "  python ligpargen.py -i your.pdb --ca-bundle %s\n"
            "Last resort: --insecure",
            exc,
            conda_pem,
        )
        sys.exit(1)

    _, _, targets = poll_until_ready(
        session,
        html,
        base_url,
        allowed_suffixes,
        wait_s=args.wait,
        max_polls=args.max_polls,
    )

    if not targets:
        logger.error(
            "No download targets found. The server may still be processing, "
            "or the submission failed. Try increasing --max-polls or run with -v."
        )
        sys.exit(1)

    logger.info("Found %d download target(s):", len(targets))
    for _post_url, _data, filename in targets:
        logger.info("  %s", filename)

    if args.no_download:
        for post_url, data, filename in targets:
            print(f"{filename}\t{post_url}\t{data}")
        return

    download_files(session, targets, Path(args.output_dir))


if __name__ == "__main__":
    main()
