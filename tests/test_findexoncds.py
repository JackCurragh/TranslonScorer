import polars as pl
import polars.testing as plt

#########################################################################################################################################
#function to test
def extract_transcript_id(attr_str):

    for attr in attr_str.split(";"):
        if attr.startswith("Parent=transcript:") or attr.startswith("ID=transcript:"):
            return attr.split(":")[1]
        elif attr.startswith("transcript_id="):
            return attr.split("=")[1]
        elif attr.startswith(" transcript_id "):
            return attr.split(" ")[2].replace('"', "")
    return ""
#########################################################################################################################################
#test data
data = pl.DataFrame(
    {'chr':['chr1','chr1','chr1','chr1','chr1'],
    'type':['gene','transcript','exon','exon','exon'],
    'start':[11869,11869,11869,12613,13221],
    'stop':[14409,14409,12227,12721,14409],
    'strand':['+','+','-','+','+'],
    'tran_id':['','ENST00000456328.2','ENST00000456328.2','ENST00000456328.2','ENST00000456328.2']}
)

df = (pl.read_csv("data/annotationsubset.gtf",separator="\t",ignore_errors=True,has_header=False,truncate_ragged_lines=True,comment_prefix="#",)
        .select(["column_1", "column_3", "column_4", "column_5", "column_7", "column_9"]).rename(    {        "column_1": "chr",        "column_3": "type",        "column_4": "start",        "column_5": "stop",        "column_7": "strand",        "column_9": "attributes",    }))

s1 =df.with_columns(pl.col("attributes").apply(lambda attributes: extract_transcript_id(attributes)).alias("tran_id")).select(pl.all().exclude("attributes"))
#########################################################################################################################################
#test gettranscript_id
def test_gettran_id():
    assert plt.assert_frame_equal(s1, data) == None
#########################################################################################################################################
#function to test
def procesexons(df):
   
    exonplus = df.filter((pl.col("strand") == "+"))
    exonneg = df.filter((pl.col("strand") == "-"))

    groupedexonspos = (
        exonplus.group_by("tran_id")
        .agg(pl.col("start"), pl.col("stop"), pl.col("strand"), pl.col("chr"))
        .select(["chr", "tran_id", "start", "stop", "strand"])
    )
    groupedexonsneg = (
        exonneg.group_by("tran_id")
        .agg(pl.col("start"), pl.col("stop"), pl.col("strand"), pl.col("chr"))
        .select(["chr", "tran_id", "start", "stop", "strand"])
    )
    return groupedexonspos
#########################################################################################################################################
#test data
datapos = pl.DataFrame(
    {'chr':[['chr1','chr1','chr1'],['chr1']],
    'start':[[11869,12613,13221],[11869]],
    'stop':[[14409,12721,14409],[14409]],
    'strand':[['+','+','+'],['+']],
    'tran_id':['ENST00000456328.2', '']}).select(["chr", "tran_id", "start", "stop", "strand"])


#########################################################################################################################################
#testing procesexons
def test_procesexons():
    assert plt.assert_frame_equal(procesexons(s1),datapos) == None
#########################################################################################################################################

def exontranscriptcoords(df: pl.DataFrame, posstrand=True) -> pl.DataFrame:

    # Initialize new columns
    new_column_1 = []
    new_column_2 = []
    # Iterate over rows
    start_column = "start"
    end_column = "stop"

    for i in range(len(df)):
        start_values = (
            df[start_column][i]
            if posstrand
            else sorted(df[start_column][i], reverse=True)
        )
        end_values = (
            df[end_column][i] if posstrand else sorted(df[end_column][i], reverse=True)
        )
        new_start_values = []  # Starting value is 0
        new_stop_values = []
        for j in range(len(start_values)):
            if j == 0:
                new_start = 0
            else:
                new_start = new_stop_values[j - 1] + 1
            # Calculate stop coordinate
            stop_coordinate = end_values[j] - start_values[j]
            new_start_values.append(new_start)
            new_stop_values.append(new_start + stop_coordinate)
        new_column_1.append(new_start_values)
        new_column_2.append(new_stop_values)

    # Add new columns to the dataframe
    df = df.with_columns((pl.Series(new_column_1)).alias("tran_start"))
    df = df.with_columns((pl.Series(new_column_2)).alias("tran_stop"))
    return df
#########################################################################################################################################
#test data
transtart = [[0,2541,2650],[0]]
transtop = [[2540,2649,3838],[2540]]
dataposmanual = datapos.with_columns((pl.Series(transtart)).alias("tran_start"),(pl.Series(transtop)).alias("tran_stop"))

exontranscriptcoords(datapos)

#########################################################################################################################################
#test exontranscriptcoords
def test_exontranscriptcoords():
    assert plt.assert_frame_equal(exontranscriptcoords(datapos),dataposmanual) == None

#########################################################################################################################################

def gettranscriptcoords(cds_df, exon_df):

    # Iterate over each row in the cds DataFrame
    tran_start = []
    tran_stop = []

    for i in range(len(cds_df)):
        # Get the transcript present in the current row of the cds DataFrame
        transcript_id = cds_df["tran_id"][i]

        # Find the corresponding row in the exon DataFrame with the same transcript_id
        exon_row = exon_df.filter(pl.col("tran_id") == transcript_id)

        if exon_row is not None:
            # Get start and stop values from cds DataFrame for the current row
            cds_start = min(cds_df["start"].gather(i)[0])
            cds_stop = max(cds_df["stop"].gather(i)[0])
            # Get corresponding exon start and stop values from exon DataFrame
            transcript_pairs = list(
                zip(exon_row["tran_start"][0], exon_row["tran_stop"][0])
            )
            exon_pairs = zip(exon_row["start"][0], exon_row["stop"][0])

            for idx, exon in enumerate(exon_pairs):
                if cds_start >= exon[0] and cds_start <= exon[1]:
                    diff_start = abs(exon[0] - cds_start)
                    transtart = transcript_pairs[idx][0] + diff_start
                    tran_start.append(transtart)

                if cds_stop >= exon[0] and cds_stop <= exon[1]:
                    diff_stop = exon[1] - cds_stop
                    transtop = transcript_pairs[idx][1] - diff_stop
                    tran_stop.append(transtop)

    df = cds_df.with_columns((pl.Series(tran_start)).alias("tran_start"))
    df_tran = df.with_columns((pl.Series(tran_stop)).alias("tran_stop"))

    return df_tran
#########################################################################################################################################
#test data



#########################################################################################################################################
#test gettranscriptcoords

#########################################################################################################################################