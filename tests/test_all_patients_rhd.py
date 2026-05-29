"""Regression test against the lab's real patient panels.

Reads RHD1 / RHD456 / RHD7 / RHD9 amplicon FASTAs from each patient folder,
runs ``RHDAnalyzer.analyze_multiple_amplicons``, and asserts the final
verdict matches the lab-confirmed RhD+/RhD- expectation.

Test data is read from a local lab directory and is not shipped with the
repo; the test is skipped (not failed) when the directory is absent.
"""

import sys
from pathlib import Path

import pytest
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from utils.rhd_analyzer import RHDAnalyzer  # noqa: E402


# Lab data root. The previous hard-coded path pointed at
# Downloads\ผล blood group NGS_term paper 2568 ใบปอ\; the actual local
# copy lives under Desktop\blood_group\.
PATIENT_DATA_ROOT = Path(r"C:\Users\ExPertComputer\Desktop\blood_group")

# Short label -> (folder name within PATIENT_DATA_ROOT, expected RhD status).
# Folder names include Thai annotations from the lab notebook.
PATIENTS = [
    ("Ake",   "Ake นศ.ป.เอก อ.น้ำผึ้ง ทำ uveitis",         "RhD+"),
    ("Baipo", "Baipoh",                                       "RhD+"),
    ("Chon",  "Chonticha",                                    "RhD+"),
    ("Dada",  "Dada",                                         "RhD+"),
    ("Earn",  "Erng น้องเอิงเพื่อนใบปอ นศ.เทอมอ.เจี๊ยบ",  "RhD+"),
    ("Eye",   "Eye",                                          "RhD+"),
    ("GG",    "Gungging นศ.ปี 3 น้องสายโคใบปอ",          "RhD+"),
    ("KT",    "Kittiphon Rh neg เพื่อนอจ ฝ้าย",          "RhD-"),
    ("NP",    "Nampeung",                                     "RhD+"),
]

# Per-patient amplicon subfolders. RHD7/RHD9 are listed for forward
# compatibility — most current panels only ship RHD1 + RHD456. RHD456_fail
# is included because Eye's RHD456 reads live there.
AMPLICON_SUBDIRS = ["RHD1", "RHD456", "RHD7", "RHD9", "RHD456_fail"]


def _collect_reads(patient_dir: Path):
    """Return list of (read_name, sequence) from every FASTA in the
    patient's RHD amplicon subdirectories."""
    reads = []
    for sub in AMPLICON_SUBDIRS:
        d = patient_dir / sub
        if not d.exists():
            continue
        for fa in sorted(d.glob("*.fasta")):
            rec = next(iter(SeqIO.parse(str(fa), "fasta")), None)
            if rec is None:
                continue
            seq = str(rec.seq).upper()
            if len(seq) < 50:
                continue
            reads.append((fa.name, seq))
    return reads


@pytest.mark.parametrize("label,folder,expected", PATIENTS,
                         ids=[p[0] for p in PATIENTS])
def test_patient_rhd_verdict(label, folder, expected):
    if not PATIENT_DATA_ROOT.exists():
        pytest.skip(f"Lab data root not available: {PATIENT_DATA_ROOT}")

    patient_dir = PATIENT_DATA_ROOT / folder
    if not patient_dir.exists():
        pytest.skip(f"Patient folder missing: {patient_dir}")

    reads = _collect_reads(patient_dir)
    assert reads, f"No RHD FASTA reads found under {patient_dir}"

    result = RHDAnalyzer().analyze_multiple_amplicons(reads)
    verdict = result["final_verdict"]
    votes = result.get("votes", {})

    # Collapse "RhD+ (confirmed/probable/..)" -> "RhD+" for comparison.
    if verdict.startswith("RhD+"):
        short = "RhD+"
    elif verdict.startswith("RhD-"):
        short = "RhD-"
    else:
        short = "Inconclusive"

    assert short == expected, (
        f"{label}: expected {expected}, got {verdict!r} "
        f"(votes +:{votes.get('RhD+', 0)} "
        f"-:{votes.get('RhD-', 0)} "
        f"?:{votes.get('Inconclusive', 0)}, "
        f"n_amplicons={result.get('total_amplicons')})"
    )
