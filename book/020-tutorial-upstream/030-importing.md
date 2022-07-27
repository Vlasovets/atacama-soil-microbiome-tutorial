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

```{usage}
def emp_directory_factory():
    import tempfile
    import requests
    import shutil
    import gzip
    import os

    import qiime2
    from q2_demux._format import EMPPairedEndDirFmt

    forward_sequence_data_url = "https://data.qiime2.org/2022.2/tutorials/atacama-soils/10p/forward.fastq.gz"
    reverse_sequence_data_url = "https://data.qiime2.org/2022.2/tutorials/atacama-soils/10p/reverse.fastq.gz"
    barcode_sequence_data_url = "https://data.qiime2.org/2022.2/tutorials/atacama-soils/10p/barcodes.fastq.gz"
    
    for url in [forward_sequence_data_url, reverse_sequence_data_url, barcode_sequence_data_url]:
        
        data = requests.get(url)
        with tempfile.NamedTemporaryFile(mode='w+b') as f:
            f.write(data.content)
            f.flush()
    
            dir_fmt = EMPPairedEndDirFmt()
            
            with gzip.open(f.name, 'rb') as f_in:
                new_file_name = os.path.splitext(f.name)[0]
                with open(new_file_name, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)

    return dir_fmt.path

data_to_import = use.init_format('data_to_import', emp_directory_factory)
```

```{usage}
from q2_demux._format import EMPPairedEndDirFmt

emp_paired_end_sequences = use.import_from_format(
    'emp_paired_end_sequences',
    semantic_type='EMPPairedEndSequences',
    variable=data_to_import,
    view_type=EMPPairedEndDirFmt)
```

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

```{usage}

barcode_sequence = use.get_metadata_column('barcode_sequence', 'barcode-sequence', sample_metadata)

use.action(
    use.UsageAction(plugin_id='demux', action_id='emp_paired'),
    use.UsageInputs(seqs=emp_paired_end_sequences, 
                    barcodes=barcode_sequence,
                    rev_comp_mapping_barcodes=True),
    use.UsageOutputNames(per_sample_sequences='demultiplexed_sequences_full.qza',
                         error_correction_details='demux_details'),
)
```

Let’s subsample the data. We will perform this subsampling in this tutorial for two reasons - 
one, to speed up the tutorial run time, and two, to demonstrate the functionality.


```{usage}

def full_factory():
    import qiime2

    a = qiime2.Artifact.load('demultiplexed_sequences', demultiplexed_sequences_full)
    return a

demultiplexed_sequences_full = use.init_artifact('demultiplexed_sequences', full_factory)

use.action(
    use.UsageAction(plugin_id='demux', action_id='subsample_paired'),
    use.UsageInputs(sequences=demultiplexed_sequences_full, 
                    fraction=0.3),
    use.UsageOutputNames(subsampled_sequences='demultiplexed_sequences_subsample'),
)
```


```{usage}

def subsample_factory():
    import qiime2

    a = qiime2.Artifact.load('demultiplexed_sequences', demultiplexed_sequences_subsample)
    return a

demultiplexed_sequences_subsample_copy = use.init_artifact('demultiplexed_sequences_subsample_copy', subsample_factory)

use.action(
    use.UsageAction(plugin_id='demux', action_id='summarize'),
    use.UsageInputs(data=demultiplexed_sequences_subsample_copy),
    use.UsageOutputNames(visualization='demultiplexed_sequences_subsample_view'),
)

```

Let’s take a look at the summary in demux-subsample.qzv. 
In the “Per-sample sequence counts” table on the “Overview” tab, there are 75 samples in the data. 
If we look at the last 20 or so rows in the table, though, we will observe 
that many samples have fewer than 100 reads in them - let’s filter those samples out of the data:

```{usage}
  
def export_factory():
    import qiime2

    a = qiime2.Visualization.load('demultiplexed_sequences_subsample_view')
    
    dirfmt = a.view(a.format)
    vzDir = str(dirfmt)
    metadata_dir = vzDir + 'per-sample-fastq-counts.tsv'
    
    
    return qiime2.Visualization.export_data(a, metadata_dir)


def filter_factory():
    import qiime2

    b = qiime2.Artifact.load('demultiplexed_sequences', demultiplexed_sequences_subsample)
    
    return b

metadata = use.init_metadata('Type[Metadata]', export_factory)

demultiplexed_sequences_top100 = use.init_artifact('demultiplexed_sequences_top100', filter_factory)

use.action(
    use.UsageAction(plugin_id='demux', action_id='filter_samples'),
    use.UsageInputs(demux=demultiplexed_sequences_top100,
                    metadata=metadata,
                    where = 'CAST([forward sequence count] AS INT) > 100'),
    use.UsageOutputNames(filtered_demux='demultiplexed_sequences_filtered'),
)
```
