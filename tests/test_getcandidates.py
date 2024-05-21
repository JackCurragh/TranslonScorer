def orfrelativeposition(annotation, df, exondf):
    cds_df, exon_coords = getexons_and_cds(annotation, exondf, list(df["tran_id"].unique()))
    print('filtering')
    filtered_df = cds_df.select(pl.all().exclude("start", "stop"))
    orftype = []
    for orfid in sorted(df["tran_id"].unique()):
        cdsregion = filtered_df.filter(pl.col("tran_id") == orfid)
        orfregion = df.filter(pl.col("tran_id") == orfid)
        orfpair = list(zip(orfregion["start"], orfregion["stop"]))
        if not cdsregion.is_empty():
            if cdsregion["tran_start"][0] > cdsregion["tran_stop"][0]:
                cdspair = [cdsregion["tran_stop"][0], cdsregion["tran_start"][0]]
            else:
                cdspair = [cdsregion["tran_start"][0], cdsregion["tran_stop"][0]]
            for orf in orfpair:
                if orf[0] < cdspair[0] and orf[1] < cdspair[0]:
                    orftype.append("uORF")
                elif orf[0] == cdspair[0] and orf[1] == cdspair[1]:
                    orftype.append("CDS")
                elif orf[0] > cdspair[1] and orf[1] > cdspair[1]:
                    orftype.append("dORF")
                elif (
                    orf[0] < cdspair[0] and orf[1] < cdspair[1] and orf[1] > cdspair[0]
                ):
                    orftype.append("uoORF")
                elif (
                    orf[0] <= cdspair[1]
                    and orf[0] >= cdspair[0]
                    and orf[1] > cdspair[1]
                ):
                    orftype.append("doORF")
                elif (
                    orf[0] < cdspair[0]
                    and orf[1] <= cdspair[0]
                ):
                    orftype.append("uoORF")
                elif orf[0] >= cdspair[0] and orf[1] <= cdspair[1]:
                    orftype.append("iORF")
                elif orf[0] < cdspair[0] and orf[1] > cdspair[1]:
                    orftype.append("eoORF")
                elif orf[0] < cdspair[0] and orf[1] == cdspair[1]:
                    orftype.append("extORF")
                else:
                    print("unexpected", orf, cdspair)
                    orftype.append("Unexpected")
        else:
            for orf in orfpair:
                orftype.append("Non Coding")
    df = df.with_columns((pl.Series(orftype)).alias("type"))
    return df, exon_coords