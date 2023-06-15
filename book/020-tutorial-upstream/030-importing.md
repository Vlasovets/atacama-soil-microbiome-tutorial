# Importing, demultiplexing and subsampling

```{usage-scope}
---
name: tutorial
---
```

```{usage-selector}
---
default-interface: galaxy-usage
---
```
In this section of the tutorial, we'll import raw fastq data that is already
demultiplexed (i.e., separated into per-sample fastq files) into a QIIME 2
artifact.
You will download three fastq.gz files, corresponding to the forward, reverse, and barcode (i.e., index) reads. 
These files contain a subset of the reads in the full data set generated for this study, 
which allows for the following commands to be run relatively quickly, however, we will perform additional subsampling in this tutorial to further improve the run time.


## Importing

We'll begin with the data import.

[//]: # ()
[//]: # (```{usage})

[//]: # (def emp_directory_factory&#40;&#41;:)

[//]: # (    import os)

[//]: # (    import tempfile)

[//]: # (    from urllib import request)

[//]: # ()
[//]: # (    from q2_demux._format import EMPPairedEndDirFmt)

[//]: # (    from q2_types.per_sample_sequences import FastqGzFormat)

[//]: # ()
[//]: # (    base_url = 'https://data.qiime2.org/2022.2/tutorials/atacama-soils/10p/')

[//]: # (    forward_sequence_data_url = base_url + "forward.fastq.gz")

[//]: # (    reverse_sequence_data_url = base_url + "reverse.fastq.gz")

[//]: # (    barcode_sequence_data_url =  base_url + "barcodes.fastq.gz"    )

[//]: # (    fmt = EMPPairedEndDirFmt&#40;mode='w'&#41;)

[//]: # ()
[//]: # (    with tempfile.TemporaryDirectory&#40;&#41; as tmpdir:)

[//]: # (        bc_fp = os.path.join&#40;tmpdir, 'barcodes.fastq.gz'&#41;)

[//]: # (        bc_fn, _ = request.urlretrieve&#40;barcode_sequence_data_url, bc_fp&#41;)

[//]: # ()
[//]: # (        forward_fp = os.path.join&#40;tmpdir, 'forward.fastq.gz'&#41;)

[//]: # (        forward_fn, _ = request.urlretrieve&#40;forward_sequence_data_url, forward_fp&#41;)

[//]: # (        )
[//]: # (        reverse_fp = os.path.join&#40;tmpdir, 'reverse.fastq.gz'&#41;)

[//]: # (        reverse_fn, _ = request.urlretrieve&#40;reverse_sequence_data_url, reverse_fp&#41;)

[//]: # (        )
[//]: # (        fmt.barcodes.write_data&#40;bc_fn, FastqGzFormat&#41;)

[//]: # (        fmt.forward.write_data&#40;forward_fn, FastqGzFormat&#41;)

[//]: # (        fmt.reverse.write_data&#40;reverse_fn, FastqGzFormat&#41;)

[//]: # ()
[//]: # (    fmt.validate&#40;&#41;)

[//]: # (    return fmt)

[//]: # (    )
[//]: # (data_to_import = use.init_format&#40;'data_to_import', emp_directory_factory&#41;)

[//]: # (```)

[//]: # ()
[//]: # (```{usage})

[//]: # (from q2_demux._format import EMPPairedEndDirFmt)

[//]: # ()
[//]: # (emp_paired_end_sequences = use.import_from_format&#40;)

[//]: # (    'emp_paired_end_sequences',)

[//]: # (    semantic_type='EMPPairedEndSequences',)

[//]: # (    variable=data_to_import,)

[//]: # (    view_type=EMPPairedEndDirFmt&#41;)

[//]: # (```)

## Generating and viewing a summary of the imported data

After the import is complete, you can generate a summary of the imported
artifact. This summary contains several important pieces of information.

First, it tells you how many sequences were obtained for each of the samples.
The  expected number of sequences per sample will vary depending on the
sequencing technology that was applied and the the number of samples that were
multiplexed in your run. You should review this, and ensure that you are
getting the expected number of sequences on average.

Second, this summary provides interactive figures that illustrate sequence
quality. This will give you an overview of the quality of your sequencing run,
and you'll need to extract information from these plots to perform quality
control on the data in the next step of the tutorial.

[//]: # (```{usage})

[//]: # ()
[//]: # (barcode_sequence = use.get_metadata_column&#40;'barcode_sequence', 'barcode-sequence', sample_metadata&#41;)

[//]: # ()
[//]: # (use.action&#40;)

[//]: # (    use.UsageAction&#40;plugin_id='demux', action_id='emp_paired'&#41;,)

[//]: # (    use.UsageInputs&#40;seqs=emp_paired_end_sequences, )

[//]: # (                    barcodes=barcode_sequence,)

[//]: # (                    rev_comp_mapping_barcodes=True&#41;,)

[//]: # (    use.UsageOutputNames&#40;per_sample_sequences='demultiplexed_sequences_full.qza',)

[//]: # (                         error_correction_details='demux_details'&#41;,)

[//]: # (&#41;)

[//]: # (```)

Let’s subsample the data. We will perform this subsampling in this tutorial for two reasons - 
one, to speed up the tutorial run time, and two, to demonstrate the functionality.

[//]: # ()
[//]: # (```{usage})

[//]: # ()
[//]: # (def full_factory&#40;&#41;:)

[//]: # (    import qiime2)

[//]: # ()
[//]: # (    a = qiime2.Artifact.load&#40;'demultiplexed_sequences', demultiplexed_sequences_full&#41;)

[//]: # (    return a)

[//]: # ()
[//]: # (demultiplexed_sequences_full = use.init_artifact&#40;'demultiplexed_sequences', full_factory&#41;)

[//]: # ()
[//]: # (use.action&#40;)

[//]: # (    use.UsageAction&#40;plugin_id='demux', action_id='subsample_paired'&#41;,)

[//]: # (    use.UsageInputs&#40;sequences=demultiplexed_sequences_full, )

[//]: # (                    fraction=0.3&#41;,)

[//]: # (    use.UsageOutputNames&#40;subsampled_sequences='demultiplexed_sequences_subsample'&#41;,)

[//]: # (&#41;)

[//]: # (```)

[//]: # ()
[//]: # ()
[//]: # (```{usage})

[//]: # ()
[//]: # (def subsample_factory&#40;&#41;:)

[//]: # (    import qiime2)

[//]: # ()
[//]: # (    a = qiime2.Artifact.load&#40;'demultiplexed_sequences', demultiplexed_sequences_subsample&#41;)

[//]: # (    return a)

[//]: # ()
[//]: # (demultiplexed_sequences_subsample_copy = use.init_artifact&#40;'demultiplexed_sequences_subsample_copy', subsample_factory&#41;)

[//]: # ()
[//]: # (use.action&#40;)

[//]: # (    use.UsageAction&#40;plugin_id='demux', action_id='summarize'&#41;,)

[//]: # (    use.UsageInputs&#40;data=demultiplexed_sequences_subsample_copy&#41;,)

[//]: # (    use.UsageOutputNames&#40;visualization='demultiplexed_sequences_subsample_view'&#41;,)

[//]: # (&#41;)

[//]: # ()
[//]: # (```)

Let’s take a look at the summary in demux-subsample.qzv. 
In the “Per-sample sequence counts” table on the “Overview” tab, there are 75 samples in the data. 
If we look at the last 20 or so rows in the table, though, we will observe 
that many samples have fewer than 100 reads in them - let’s filter those samples out of the data:

[//]: # ()
[//]: # (```{usage})

[//]: # (  )
[//]: # (def export_factory&#40;&#41;:)

[//]: # (    import qiime2)

[//]: # ()
[//]: # (    a = qiime2.Visualization.load&#40;'demultiplexed_sequences_subsample_view'&#41;)

[//]: # (    )
[//]: # (    dirfmt = a.view&#40;a.format&#41;)

[//]: # (    vzDir = str&#40;dirfmt&#41;)

[//]: # (    metadata_dir = vzDir + 'per-sample-fastq-counts.tsv')

[//]: # (    )
[//]: # (    )
[//]: # (    return qiime2.Visualization.export_data&#40;a, metadata_dir&#41;)

[//]: # ()
[//]: # ()
[//]: # (def filter_factory&#40;&#41;:)

[//]: # (    import qiime2)

[//]: # ()
[//]: # (    b = qiime2.Artifact.load&#40;'demultiplexed_sequences', demultiplexed_sequences_subsample&#41;)

[//]: # (    )
[//]: # (    return b)

[//]: # ()
[//]: # (metadata = use.init_metadata&#40;'Type[Metadata]', export_factory&#41;)

[//]: # ()
[//]: # (demultiplexed_sequences_top100 = use.init_artifact&#40;'demultiplexed_sequences_top100', filter_factory&#41;)

[//]: # ()
[//]: # (use.action&#40;)

[//]: # (    use.UsageAction&#40;plugin_id='demux', action_id='filter_samples'&#41;,)

[//]: # (    use.UsageInputs&#40;demux=demultiplexed_sequences_top100,)

[//]: # (                    metadata=metadata,)

[//]: # (                    where = 'CAST&#40;[forward sequence count] AS INT&#41; > 100'&#41;,)

[//]: # (    use.UsageOutputNames&#40;filtered_demux='demultiplexed_sequences_filtered'&#41;,)

[//]: # (&#41;)

[//]: # (```)
