import gc

import polars as pl

from ..utils.logging import log_error, log_info

try:
    # Avoid importing heavy dependencies at module import unless needed
    from .coordinates import classify_orf
except Exception:
    # Fallback placeholder; actual import occurs in functions that need it
    def classify_orf(row):
        return "Unexpected"


from ..file_handlers.bam import getexons_and_cds


def sru_score(position, tran_reads, sru_range, direction):
    """
    Calculate the SRU (Start/Stop Ramp Up/Down) score for a given position.
    Memory-efficient implementation using lazy evaluation.

    Args:
        position (int): Position to calculate score for
        tran_reads (DataFrame): Transcript reads data
        sru_range (int): Range for calculation
        direction (int): 0 for start (up), 1 for stop (down)

    Returns:
        float: Calculated SRU score
    """
    try:
        # To ensure lazy evaluation works consistently
        if not isinstance(tran_reads, pl.LazyFrame):
            lazy_reads = tran_reads.lazy()
        else:
            lazy_reads = tran_reads

        if direction == 0:  # Start (rise up)
            # Calculate before mean - positions before 'position'
            before_filter = lazy_reads.filter(
                (pl.col("tran_start") >= position - sru_range) & (pl.col("tran_start") < position)
            )

            # Calculate after mean - positions after 'position'
            after_filter = lazy_reads.filter(
                (pl.col("tran_start") >= position) & (pl.col("tran_start") < position + sru_range)
            )
        else:  # Stop (step down)
            # Calculate before mean - positions before 'position + sru_range'
            before_filter = lazy_reads.filter(
                (pl.col("tran_start") >= position) & (pl.col("tran_start") < position + sru_range)
            )

            # Calculate after mean - positions after 'position + sru_range'
            after_filter = lazy_reads.filter(
                (pl.col("tran_start") >= position + sru_range)
                & (pl.col("tran_start") < position + (2 * sru_range))
            )

        # Calculate means with proper error handling
        try:
            before_mean_df = before_filter.select(pl.mean("counts").alias("mean")).collect()
            before_mean = before_mean_df[0, 0] if before_mean_df.height > 0 else 0.0
        except Exception:
            before_mean = 0.0

        try:
            after_mean_df = after_filter.select(pl.mean("counts").alias("mean")).collect()
            after_mean = after_mean_df[0, 0] if after_mean_df.height > 0 else 0.0
        except Exception:
            after_mean = 0.0

        # Handle None values
        before_mean = 0.0 if before_mean is None else before_mean
        after_mean = 0.0 if after_mean is None else after_mean

        # Return the appropriate score based on direction
        if direction == 0:
            return after_mean - before_mean
        else:
            return before_mean - after_mean

    except Exception as e:
        log_error(f"Error calculating SRU score: {str(e)}")
        return 0.0


def calculate_scores(start, stop, tran_reads):
    """
    Calculate HRF, average, and non-zero coverage scores for a region.
    Uses lazy evaluation for memory efficiency.

    Args:
        start (int): Start position
        stop (int): Stop position
        tran_reads (DataFrame): Transcript reads data

    Returns:
        tuple: (HRF score, average score, non-zero coverage score)
    """
    try:
        # Convert to lazy - use to_lazy() directly instead of checking is_lazy()
        lazy_reads = tran_reads.lazy()

        # Filter the region using lazy evaluation
        region_counts = lazy_reads.filter(
            (pl.col("tran_start") >= start) & (pl.col("tran_start") <= stop)
        ).select("counts")

        # Calculate metrics separately to avoid issues with nzc calculation
        try:
            # Compute max and mean
            basic_metrics = region_counts.select(
                [pl.col("counts").max().alias("hrf"), pl.col("counts").mean().alias("avg")]
            ).collect()

            if basic_metrics.height == 0:
                return (0.0, 0.0, 0.0)

            # Extract results
            hrf = (
                basic_metrics.get_column("hrf")[0]
                if basic_metrics.get_column("hrf")[0] is not None
                else 0.0
            )
            avg = (
                basic_metrics.get_column("avg")[0]
                if basic_metrics.get_column("avg")[0] is not None
                else 0.0
            )

            # Calculate nzc separately with proper error handling
            nzc_counts = region_counts.select(
                [
                    (pl.col("counts") > 0).sum().alias("nonzero_count"),
                    pl.count().alias("total_count"),
                ]
            ).collect()

            nonzero_count = (
                nzc_counts.get_column("nonzero_count")[0]
                if nzc_counts.get_column("nonzero_count")[0] is not None
                else 0
            )
            total_count = (
                nzc_counts.get_column("total_count")[0]
                if nzc_counts.get_column("total_count")[0] is not None
                else 1
            )

            # Avoid division by zero
            nzc = nonzero_count / total_count if total_count > 0 else 0.0
        except Exception as e:
            log_error(f"Error calculating specific metric: {str(e)}")
            return (0.0, 0.0, 0.0)

        return (hrf, avg, nzc)

    except Exception as e:
        log_error(f"Error calculating region scores: {str(e)}")
        return (0.0, 0.0, 0.0)


def process_orf_chunk(chunk, tran_reads, sru_range, typeorf):
    """
    Process a chunk of ORFs to reduce memory usage.

    Args:
        chunk (DataFrame): Chunk of ORF DataFrame to process
        tran_reads (DataFrame): Transcript read data
        sru_range (int): Range for SRU calculation
        typeorf (str): Type of ORF

    Returns:
        DataFrame: Processed ORF DataFrame with scores
    """
    # Convert to lazy for memory efficiency - no is_lazy() check
    lazy_tran_reads = tran_reads.lazy()

    # Initialize result lists
    rise_up_values = []
    step_down_values = []
    hrf_values = []
    avg_values = []
    nzc_values = []

    # Process each ORF in the chunk
    for row in chunk.iter_rows(named=True):
        start = row["start"]
        stop = row["stop"]

        # Calculate SRU scores based on ORF type
        rise_up = 0.0
        step_down = 0.0

        if typeorf == "uoORF" or typeorf not in ("uoORF", "doORF"):
            rise_up = sru_score(start, lazy_tran_reads, sru_range, 0)

        if typeorf == "doORF" or typeorf not in ("uoORF", "doORF"):
            step_down = sru_score(stop, lazy_tran_reads, sru_range, 1)

        # Calculate region scores
        hrf, avg, nzc = calculate_scores(start, stop, lazy_tran_reads)

        # Append to result lists
        rise_up_values.append(rise_up)
        step_down_values.append(step_down)
        hrf_values.append(hrf)
        avg_values.append(avg)
        nzc_values.append(nzc)

    # Add score columns to the DataFrame using with_columns
    return chunk.with_columns(
        [
            pl.Series("rise_up", rise_up_values),
            pl.Series("step_down", step_down_values),
            pl.Series("hrf", hrf_values),
            pl.Series("avg", avg_values),
            pl.Series("nzc", nzc_values),
            pl.sum_horizontal(
                pl.Series("rise_up", rise_up_values),
                pl.Series("step_down", step_down_values),
                pl.Series("hrf", hrf_values),
                pl.Series("avg", avg_values),
                pl.Series("nzc", nzc_values),
            ).alias("score"),
        ]
    )


def oldscoring(df, tran_reads, sru_range, typeorf):
    """
    Classic scoring method with memory-efficient implementation using
    chunked processing and lazy evaluation.

    Args:
        df (DataFrame): Input DataFrame with ORF information
        tran_reads (DataFrame): Transcript reads data
        sru_range (int): Range for SRU score calculation
        typeorf (str): Type of ORF ('uoORF', 'doORF', or other)

    Returns:
        DataFrame: DataFrame with calculated scores
    """
    try:
        # Define optimal chunk size based on available memory
        # Smaller chunks use less memory but may be slower
        chunk_size = 50

        # Process in chunks and collect results
        result_chunks = []

        for i in range(0, len(df), chunk_size):
            # Extract chunk
            end = min(i + chunk_size, len(df))
            chunk = df.slice(i, end - i)

            # Process chunk
            processed_chunk = process_orf_chunk(chunk, tran_reads, sru_range, typeorf)

            # Add to results
            result_chunks.append(processed_chunk)

            # Log progress
            log_info(f"Processed ORFs {i} to {end-1} out of {len(df)}")

            # Force garbage collection
            del chunk
            gc.collect()

        # Combine results
        if result_chunks:
            return pl.concat(result_chunks)
        else:
            return pl.DataFrame()

    except Exception as e:
        log_error(f"Error in memory-efficient scoring: {str(e)}")
        return pl.DataFrame()


def batch_cache_positions(positions, tran_reads, sru_range, direction, scoredict, key):
    """
    Cache scores for a batch of positions to improve memory efficiency.

    Args:
        positions (list): List of positions to calculate scores for
        tran_reads (DataFrame): Transcript read data
        sru_range (int): Range for SRU calculation
        direction (int): 0 for start (up), 1 for stop (down)
        scoredict (dict): Score dictionary to update
        key (str): Dictionary key ('rise_up' or 'step_down')

    Returns:
        dict: Updated score dictionary
    """
    # Convert to lazy for memory efficiency - no is_lazy() check
    lazy_tran_reads = tran_reads.lazy()

    # Filter positions not in cache
    new_positions = [pos for pos in positions if pos not in scoredict[key]]

    # Calculate scores for new positions
    for pos in new_positions:
        scoredict[key][pos] = sru_score(pos, lazy_tran_reads, sru_range, direction)

    return scoredict


def newscoring(df, tran_reads, sru_range, typeorf, scoredict):
    """
    Modern scoring method that caches scores in a dictionary.
    Memory-efficient implementation using lazy evaluation and chunked processing.

    Args:
        df (DataFrame/Series): Input data with ORF information
        tran_reads (DataFrame): Transcript reads data
        sru_range (int): Range for SRU score calculation
        typeorf (str): Type of ORF ('uoORF', 'doORF', or other)
        scoredict (dict): Dictionary to cache scores

    Returns:
        dict: Updated score dictionary
    """
    try:
        # Increase batch size for better performance while managing memory
        batch_size = 20

        # Process start positions for rise_up scores
        if not typeorf == "doORF":
            if isinstance(df, pl.DataFrame):
                # Get unique start positions
                startvalues = df.get_column("start").unique().to_list()

                # Process in batches
                for i in range(0, len(startvalues), batch_size):
                    batch_end = min(i + batch_size, len(startvalues))
                    batch = startvalues[i:batch_end]

                    # Update score cache for batch
                    scoredict = batch_cache_positions(
                        batch, tran_reads, sru_range, 0, scoredict, "rise_up"
                    )

                    # Clean up
                    del batch
                    gc.collect()
            else:
                # Handle Series case
                startvalues = df.unique().to_list()

                # Process in batches
                for i in range(0, len(startvalues), batch_size):
                    batch_end = min(i + batch_size, len(startvalues))
                    batch = startvalues[i:batch_end]

                    # Update score cache for batch
                    scoredict = batch_cache_positions(
                        batch, tran_reads, sru_range, 0, scoredict, "rise_up"
                    )

                    # Clean up
                    del batch
                    gc.collect()

        # Process stop positions for step_down scores
        if not typeorf == "uoORF":
            if isinstance(df, pl.DataFrame):
                # Get unique stop positions
                stopvalues = df.get_column("stop").unique().to_list()

                # Process in batches
                for i in range(0, len(stopvalues), batch_size):
                    batch_end = min(i + batch_size, len(stopvalues))
                    batch = stopvalues[i:batch_end]

                    # Update score cache for batch
                    scoredict = batch_cache_positions(
                        batch, tran_reads, sru_range, 1, scoredict, "step_down"
                    )

                    # Clean up
                    del batch
                    gc.collect()
            else:
                # Handle Series case
                stopvalues = df.unique().to_list()

                # Process in batches
                for i in range(0, len(stopvalues), batch_size):
                    batch_end = min(i + batch_size, len(stopvalues))
                    batch = stopvalues[i:batch_end]

                    # Update score cache for batch
                    scoredict = batch_cache_positions(
                        batch, tran_reads, sru_range, 1, scoredict, "step_down"
                    )

                    # Clean up
                    del batch
                    gc.collect()

        return scoredict
    except Exception as e:
        log_error(f"Error in memory-efficient new scoring: {str(e)}")
        return {"rise_up": {}, "step_down": {}}


def process_global_scores_chunk(chunk, tran_reads, typeorf):
    """
    Process a chunk of ORFs for global scores to reduce memory usage.

    Args:
        chunk (DataFrame): Chunk of ORF DataFrame to process
        tran_reads (DataFrame): Transcript read data
        typeorf (str): Type of ORF

    Returns:
        DataFrame: Processed ORF DataFrame with global scores
    """
    # Convert to lazy for memory efficiency - no is_lazy() check
    lazy_tran_reads = tran_reads.lazy()

    # Initialize result lists
    hrf_values = []
    avg_values = []
    nzc_values = []

    # Process each ORF in the chunk
    for row in chunk.iter_rows(named=True):
        start = row["start"]
        stop = row["stop"]

        # Calculate region scores
        hrf, avg, nzc = calculate_scores(start, stop, lazy_tran_reads)

        # Append to result lists
        hrf_values.append(hrf)
        avg_values.append(avg)
        nzc_values.append(nzc)

    # Add score columns and calculate total score based on ORF type
    result = chunk.with_columns(
        [pl.Series("hrf", hrf_values), pl.Series("avg", avg_values), pl.Series("nzc", nzc_values)]
    )

    # Calculate total score based on ORF type
    if typeorf == "doORF":
        result = result.with_columns(
            pl.sum_horizontal("step_down", "hrf", "avg", "nzc").alias("score")
        )
    elif typeorf == "uoORF":
        result = result.with_columns(
            pl.sum_horizontal("rise_up", "hrf", "avg", "nzc").alias("score")
        )
    else:
        result = result.with_columns(
            pl.sum_horizontal("rise_up", "step_down", "hrf", "avg", "nzc").alias("score")
        )

    return result


def globalscores(df, tran_reads, typeorf):
    """
    Calculate global scores (HRF, average, NZC) for ORFs.
    Memory-efficient implementation using chunked processing and lazy evaluation.

    Args:
        df (DataFrame): Input DataFrame with ORF information
        tran_reads (DataFrame): Transcript reads data
        typeorf (str): Type of ORF ('uoORF', 'doORF', or other)

    Returns:
        DataFrame: DataFrame with calculated global scores
    """
    try:
        # Calculate dynamic chunk size
        desired_batches = 10  # Adjust this as needed
        chunk_size = max(1, len(df) // desired_batches)  # Ensure at least 1

        result_chunks = []

        for i in range(0, len(df), chunk_size):
            end = min(i + chunk_size, len(df))
            chunk = df.slice(i, end - i)

            # Process chunk
            processed_chunk = process_global_scores_chunk(chunk, tran_reads, typeorf)

            # Add to results
            result_chunks.append(processed_chunk)
            # Force garbage collection
            del chunk
            gc.collect()
        # Combine results
        if result_chunks:
            return pl.concat(result_chunks)
        else:
            return pl.DataFrame()
    except Exception as e:
        log_error(f"Error in memory-efficient global scores: {str(e)}")
        return pl.DataFrame()


def existingscore(df, typeorf, scoredict):
    """
    Filter out ORFs that already have scores in the cache.
    Memory-efficient implementation using expressions.

    Args:
        df (DataFrame): Input DataFrame with ORF information
        typeorf (str): Type of ORF ('uoORF', 'doORF', or other)
        scoredict (dict): Dictionary of cached scores

    Returns:
        DataFrame/Series: Filtered data containing only ORFs needing scoring
    """
    try:
        # Convert score dictionaries to sets for faster lookups
        rise_up_keys = set(scoredict["rise_up"].keys())
        step_down_keys = set(scoredict["step_down"].keys())

        if typeorf == "uoORF":
            # Use Polars expressions for filtering
            return df.filter(~pl.col("start").is_in(rise_up_keys))
        elif typeorf == "doORF":
            # Use Polars expressions for filtering
            return df.filter(~pl.col("stop").is_in(step_down_keys))
        else:
            # Create a DataFrame with only rows that need scoring
            needs_scoring = df.filter(
                (~pl.col("start").is_in(rise_up_keys)) | (~pl.col("stop").is_in(step_down_keys))
            )

            return needs_scoring.select(["start", "stop"])
    except Exception as e:
        log_error(f"Error checking existing scores: {str(e)}")
        return pl.DataFrame()


def assigningscore(df, scoredict, typeorf):
    """
    Assign cached scores to ORFs.
    Memory-efficient implementation using chunked processing.

    Args:
        df (DataFrame): Input DataFrame with ORF information
        scoredict (dict): Dictionary of cached scores
        typeorf (str): Type of ORF ('uoORF', 'doORF', or other)

    Returns:
        DataFrame: DataFrame with scores assigned from cache
    """
    try:
        # Process in chunks to reduce memory usage
        chunk_size = 10
        result_chunks = []

        for i in range(0, len(df), chunk_size):
            # Extract chunk
            end = min(i + chunk_size, len(df))
            chunk = df.slice(i, end - i)

            # Prepare score assignment lists
            rise_up_values = []
            step_down_values = []

            # Process each ORF in the chunk
            for row in chunk.iter_rows(named=True):
                if typeorf == "uoORF":
                    rise_up_values.append(scoredict["rise_up"].get(row["start"], 0.0))
                    step_down_values.append(0.0)
                elif typeorf == "doORF":
                    rise_up_values.append(0.0)
                    step_down_values.append(scoredict["step_down"].get(row["stop"], 0.0))
                else:
                    rise_up_values.append(scoredict["rise_up"].get(row["start"], 0.0))
                    step_down_values.append(scoredict["step_down"].get(row["stop"], 0.0))

            # Add score columns to the chunk
            processed_chunk = chunk.with_columns(
                [pl.Series("rise_up", rise_up_values), pl.Series("step_down", step_down_values)]
            )

            # Add to results
            result_chunks.append(processed_chunk)

            # Clean up
            del chunk, rise_up_values, step_down_values
            gc.collect()

        # Combine results
        if result_chunks:
            return pl.concat(result_chunks)
        else:
            return pl.DataFrame()
    except Exception as e:
        log_error(f"Error assigning scores: {str(e)}")
        return pl.DataFrame()


def orfrelativeposition(annotation, df, cds_df):
    """
    Determines the relative position of ORFs to coding sequences (CDS).
    Memory-efficient implementation using chunked processing.

    Parameters:
        annotation (str): Path to the genome annotation file in BED/GFF/GTF format.
        df (polars.DataFrame): DataFrame containing ORF coordinates.
        cds_df (polars.DataFrame): DataFrame containing CDS coordinates.

    Returns:
        tuple: A tuple containing two polars DataFrames:
               - The first DataFrame contains ORF coordinates with an additional column
                 'type' indicating the relative position of each ORF to CDS.
               - The second DataFrame contains exon coordinates.
    """
    if "cdsdf" not in globals():
        cds_df, exon_coords = getexons_and_cds(annotation, list(df["tran_id"].unique()))

    tranids = set(cds_df["tran_id"].unique().to_list())  # Convert to set for faster lookups

    # Join df with cds_df to include cds_start and cds_stop
    df = df.join(cds_df.select(["tran_id", "tran_start", "tran_stop"]), on="tran_id", how="left")

    # Process in chunks for memory efficiency
    chunk_size = 1000
    result_chunks = []

    for i in range(0, len(df), chunk_size):
        # Extract chunk
        end = min(i + chunk_size, len(df))
        chunk = df.slice(i, end - i)

        # Prepare type values
        type_values = []

        # Process each ORF in the chunk
        for row in chunk.iter_rows(named=True):
            if row["tran_id"] in tranids:
                try:
                    type_values.append(
                        classify_orf(
                            {
                                "start": row["start"],
                                "stop": row["stop"],
                                "tran_start": row["tran_start"],
                                "tran_stop": row["tran_stop"],
                            }
                        )
                    )
                except Exception:
                    type_values.append("Non Coding")
            else:
                type_values.append("Non Coding")

        # Add type column to the chunk
        processed_chunk = chunk.with_columns(pl.Series("type", type_values))

        # Add to results
        result_chunks.append(processed_chunk)

        # Clean up
        del chunk, type_values
        gc.collect()

    # Combine results
    if result_chunks:
        df = pl.concat(result_chunks)

    # Filter CDS more efficiently
    df.filter(pl.col("type") == "CDS")["tran_id"].unique().to_list()

    return df, exon_coords
