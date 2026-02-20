import polars as pl
from TranslonScorer.core import coordinates


def test_classify_orf_labels():
    # Build minimal rows
    def row(start, stop, ts, te):
        return {"start": start, "stop": stop, "tran_start": ts, "tran_stop": te}

    assert coordinates.classify_orf(row(0, 8, 10, 20)) == "uORF"
    assert coordinates.classify_orf(row(10, 20, 10, 20)) == "CDS"
    assert coordinates.classify_orf(row(25, 30, 10, 20)) == "dORF"
    assert coordinates.classify_orf(row(5, 12, 10, 20)) == "uoORF"
    assert coordinates.classify_orf(row(15, 25, 10, 20)) == "doORF"
    assert coordinates.classify_orf(row(12, 18, 10, 20)) == "iORF"
    assert coordinates.classify_orf(row(5, 25, 10, 20)) == "eoORF"
    assert coordinates.classify_orf(row(5, 20, 10, 20)) == "extORF"

