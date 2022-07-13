# Importing demultiplexed sequence data

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
def single_directory_factory():
    import tempfile
    import requests
    import shutil
    import gzip
    import os

    import qiime2
    from q2_types.per_sample_sequences import \
        CasavaOneEightSingleLanePerSampleDirFmt

    forward_sequence_data_url = 'https://data.qiime2.org/2022.2/tutorials/atacama-soils/10p/forward.fastq.gz'
    reverse_sequence_data_url = 'https://data.qiime2.org/2022.2/tutorials/atacama-soils/10p/forward.fastq.gz'
    barcode_sequence_data_url = 'https://data.qiime2.org/2022.2/tutorials/atacama-soils/10p/reverse.fastq.gz'
    
    for url in [forward_sequence_data_url, reverse_sequence_data_url, barcode_sequence_data_url]:
        
        data = requests.get(url)
        with tempfile.NamedTemporaryFile(mode='w+b') as f:
            f.write(data.content)
            f.flush()
    
            dir_fmt = CasavaOneEightSingleLanePerSampleDirFmt()
            
            with gzip.open(f.name, 'rb') as f_in:
                new_file_name = os.path.splitext(f.name)[0]
                with open(f.name, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)

    return dir_fmt

data_to_import = use.init_format('data_to_import', single_directory_factory)
```

```{usage}
from q2_types.per_sample_sequences import \
    CasavaOneEightSingleLanePerSampleDirFmt

emp_paired_end_sequences = use.import_from_format(
    'emp_paired_end_sequences',
    semantic_type='EMPPairedEndSequences',
    variable=data_to_import,
    view_type=CasavaOneEightSingleLanePerSampleDirFmt)
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

use.action(
    use.UsageAction(plugin_id='demux', action_id='emp_paired'),
    use.UsageInputs(seqs=emp_paired_end_sequences, 
                    barcodes=sample_metadata.get_column('barcode-sequence'),
                    rev_comp_mapping_barcodes=True),
    use.UsageOutputNames(per_sample_sequences='demultiplexed_sequences_full',
                         error_correction_details='demux_details'),
)
```
