import polars as pl
from TranslonScorer.core import orffinder


def test_preporfs_on_dict_sequence():
    seqs = {
        "tx1": "AAAATGAAATAGCCC",  # ATG ... TAG in-frame
        "tx2": "CCCATGAAAAAATGA",  # ATG ... TGA in-frame
    }
    df = orffinder.preporfs(
        seqs, start_codons=["ATG"], stop_codons=["TAA", "TAG", "TGA"], minlength=0, maxlength=1000
    )
    assert isinstance(df, pl.DataFrame)
    # Expect at least one ORF per transcript
    assert set(df.get_column("tran_id").to_list()) == {"tx1", "tx2"}
    # Basic schema presence
    for col in ["tran_id", "start", "stop", "length", "startorf", "stoporf"]:
        assert col in df.columns
