"""Guards #5: file processors must surface failures, not silently drop them.

These run without importing streamlit (the whole point of the #8 extraction):
the processors now live in ``utils.file_processing`` and accept an optional
``errors`` list that collects per-file failures.
"""

import io
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from utils.file_processing import (  # noqa: E402
    process_rhd_fasta_files,
    process_ab1_files,
)


class _BrokenUpload:
    """File-like whose read() raises, mimicking a corrupt/unreadable upload."""

    def __init__(self, name):
        self.name = name

    def read(self, *a, **k):
        raise IOError("simulated unreadable file")


def _fasta_upload(name, seq):
    buf = io.BytesIO(f">{name}\n{seq}\n".encode("utf-8"))
    buf.name = name
    return buf


def test_corrupt_fasta_is_recorded_not_dropped():
    errors = []
    good = _fasta_upload("good.fasta", "ACGT" * 30)  # >50 bp
    bad = _BrokenUpload("corrupt.fasta")

    traces = process_rhd_fasta_files([good, bad], errors=errors)

    # Good file still processed; bad file surfaced with its name.
    assert traces is not None and len(traces) == 1
    assert len(errors) == 1
    assert errors[0]["filename"] == "corrupt.fasta"
    assert "simulated unreadable file" in errors[0]["error"]


def test_errors_optional_preserves_legacy_signature():
    # Without an errors list, behaviour is unchanged (no crash, failure dropped).
    bad = _BrokenUpload("corrupt.fasta")
    assert process_rhd_fasta_files([bad]) is None


def test_corrupt_ab1_is_recorded():
    errors = []
    bad = _BrokenUpload("corrupt.ab1")
    results, hets = process_ab1_files([bad], [], errors=errors)
    assert results is None
    assert len(errors) == 1
    assert errors[0]["filename"] == "corrupt.ab1"
