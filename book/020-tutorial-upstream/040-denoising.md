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
```{admonition} Guidance on choosing these values

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

[//]: # ()
[//]: # (```{usage})

[//]: # ()
[//]: # (def demux_factory&#40;&#41;:)

[//]: # (    import qiime2)

[//]: # ()
[//]: # (    a = qiime2.Artifact.load&#40;'demultiplexed_sequences', demultiplexed_sequences_filtered&#41;)

[//]: # (    return a)

[//]: # ()
[//]: # (demux = use.init_artifact&#40;'demux', demux_factory&#41;)

[//]: # ()
[//]: # (rep_seqs, table, denoising_stats = use.action&#40;)

[//]: # (    use.UsageAction&#40;plugin_id='dada2', action_id='denoise_paired'&#41;,)

[//]: # (    use.UsageInputs&#40;demultiplexed_seqs=demux,)

[//]: # (                    trim_left_f=13, )

[//]: # (                    trim_left_r=13, )

[//]: # (                    trunc_len_f=150,)

[//]: # (                    trunc_len_r=150,&#41;,)

[//]: # (    use.UsageOutputNames&#40;representative_sequences='rep-seqs',)

[//]: # (                        table='table',)

[//]: # (                        denoising_stats='denoising_stats'&#41;)

[//]: # (&#41;)

[//]: # (```)

## Reviewing the DADA2 run statistics

At this stage, you will have artifacts containing the feature table, 
corresponding feature sequences, and DADA2 denoising stats. 
You can generate summaries of these as follows.

[//]: # (```{usage})

[//]: # (stats_as_md = use.view_as_metadata&#40;'stats_dada2_md', denoising_stats&#41;)

[//]: # ()
[//]: # (use.action&#40;)

[//]: # (    use.UsageAction&#40;plugin_id='metadata', action_id='tabulate'&#41;,)

[//]: # (    use.UsageInputs&#40;input=stats_as_md&#41;,)

[//]: # (    use.UsageOutputNames&#40;visualization='dada2_stats_summ'&#41;)

[//]: # (&#41;)

[//]: # (```)

## Generating and reviewing summaries of the feature table and feature data

The next two outputs of DADA2 will form the basis of the majority of the
microbiome analyses that you'll run, in connection with your sample metadata.
This is the feature table and feature data. The feature table describes which
amplicon sequence variants (ASVs) were observed in which samples, and how many
times each ASV was observed in each sample. The feature data in this case is
the sequence that defines each ASV. Generate and explore the summaries of
each of these files.

[//]: # ()
[//]: # (```{usage})

[//]: # (use.action&#40;)

[//]: # (    use.UsageAction&#40;plugin_id='feature_table', action_id='summarize'&#41;,)

[//]: # (    use.UsageInputs&#40;table=table, sample_metadata=sample_metadata&#41;,)

[//]: # (    use.UsageOutputNames&#40;visualization='table_summ'&#41;,)

[//]: # (&#41;)

[//]: # ()
[//]: # (use.action&#40;)

[//]: # (    use.UsageAction&#40;plugin_id='feature_table', action_id='tabulate_seqs'&#41;,)

[//]: # (    use.UsageInputs&#40;data=rep_seqs&#41;,)

[//]: # (    use.UsageOutputNames&#40;visualization='rep_seqs_summ'&#41;,)

[//]: # (&#41;)

[//]: # ()
[//]: # (```)