# Denoising sequence data with DADA2

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

## Performing sequence quality control (i.e., denoising)
Next, we’ll look at the sequence quality based on ten-thousand randomly selected reads from 
the subsampled and filtered data, and then denoise the data. When you view the quality plots, 
note that in contrast to the corresponding plots in the moving pictures tutorial, 
there are now two interactive plots to be considered together. 
The plot on the left presents the quality scores for the forward reads, and the plot on the right presents 
the quality scores for the reverse reads. We’ll use these plots to determine what trimming parameters we want to use 
for denoising with DADA2 {cite:t}`callahan-dada2-2016`, which is accessible through the q2-dada2
plugin. Since our reads are paired end, we'll use the `denoise_paired` action in the
q2-dada2 plugin. This performs quality filtering, chimera checking, and paired-
end read joining.

In this example we have 150-base forward and reverse reads. 
Since we need the reads to be long enough to overlap when joining paired ends, 
the first thirteen bases of the forward and reverse reads are being trimmed, 
but no trimming is being applied to the ends of the sequences to avoid reducing the read length by too much. 
In this example, the same values are being provided for `--p-trim-left-f` and `--p-trim-left-r` 
and for `--p-trunc-len-f` and `--p-trunc-len-r`, but that is not a requirement.


Spend a couple of minutes reviewing the quality score plots and think about
where you might want to truncate the forward and reverse reads, and if you'd
like to trim any bases from the beginnings.

````{margin}
```{admonition} Greg's guidance on choosing these values

I typically try to apply some objective criteria when selecting these values.
For example, in reviewing the quality score plots, I noticed that the
twenty-fifth percentile quality score drops below 13 at position 204 in the
forward reads and 205 in the reverse reads. I chose to use those values for
the required truncation lengths.

Since the first base of the reverse reads is
slightly lower than those that follow, I choose to trim that first base in the
reverse reads, but apply no trimming to the forward reads. This trimming is
probably unnecessary here, but is useful here for illustrating how this works.
```
````

```{usage}

def demux_factory():
    import qiime2

    a = qiime2.Artifact.load('demultiplexed_sequences', demultiplexed_sequences_filtered)
    return a

demux = use.init_artifact('demux', demux_factory)

rep_seqs, table, denoising_stats = use.action(
    use.UsageAction(plugin_id='dada2', action_id='denoise_paired'),
    use.UsageInputs(demultiplexed_seqs=demux,
                    trim_left_f=13, 
                    trim_left_r=13, 
                    trunc_len_f=150,
                    trunc_len_r=150,),
    use.UsageOutputNames(representative_sequences='rep-seqs',
                        table='table',
                        denoising_stats='denoising_stats')
)
```

## Reviewing the DADA2 run statistics

At this stage, you will have artifacts containing the feature table, 
corresponding feature sequences, and DADA2 denoising stats. 
You can generate summaries of these as follows.

```{usage}
stats_as_md = use.view_as_metadata('stats_dada2_md', denoising_stats)

use.action(
    use.UsageAction(plugin_id='metadata', action_id='tabulate'),
    use.UsageInputs(input=stats_as_md),
    use.UsageOutputNames(visualization='dada2_stats_summ')
)
```

## Generating and reviewing summaries of the feature table and feature data

The next two outputs of DADA2 will form the basis of the majority of the
microbiome analyses that you'll run, in connection with your sample metadata.
This is the feature table and feature data. The feature table describes which
amplicon sequence variants (ASVs) were observed in which samples, and how many
times each ASV was observed in each sample. The feature data in this case is
the sequence that defines each ASV. Generate and explore the summaries of
each of these files.

```{usage}
use.action(
    use.UsageAction(plugin_id='feature_table', action_id='summarize'),
    use.UsageInputs(table=table, sample_metadata=sample_metadata),
    use.UsageOutputNames(visualization='table_summ'),
)

use.action(
    use.UsageAction(plugin_id='feature_table', action_id='tabulate_seqs'),
    use.UsageInputs(data=rep_seqs),
    use.UsageOutputNames(visualization='rep_seqs_summ'),
)
```

```{note}
We've now reached the end of the **upstream** tutorial. When we begin working
on the **downstream** tutorial, we'll work with larger feature table and
feature data artifacts representing many more samples. The samples that we
worked with in this tutorial are a small subset of what we'll work with in the
**downstream** tutorial.
```
