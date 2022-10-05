# Atacama Soil Microbiome Data

This tutorial is designed to serve two purposes. 
First, it illustrates the initial processing steps of paired-end read analysis, up to the point where the analysis steps
are identical to single-end read analysis. This includes the importing, demultiplexing, and denoising steps, 
and results in a feature table and the associated feature sequences. Second, this is intended to be a self-guided exercise
that could be run after the moving pictures tutorial to gain more experience with QIIME 2. For this exercise, 
we provide some questions that can be used to guide your analysis, but do not provide commands that will allow you to address each. 
Instead, you should apply the commands that you learned in the moving pictures tutorial.

In this tutorial you’ll use QIIME 2 to perform an analysis of soil samples from the Atacama Desert in northern Chile. 
The Atacama Desert is one of the most arid locations on Earth, with some areas receiving less than a millimeter of rain per decade. 
Despite this extreme aridity, there are microbes living in the soil. The soil microbiomes profiled in this study follow
two east-west transects, Baquedano and Yungay, across which average soil relative humidity is positively correlated with
elevation (higher elevations are less arid and thus have higher average soil relative humidity). Along these transects, 
pits were dug at each site and soil samples were collected from three depths in each pit.

This tutorial focuses on data reused from [Fuentes et al (2021) Influence of Physical-Chemical Soil Parameters on Microbiota Composition and Diversity 
             in a Deep Hyperarid Core of the Atacama Desert](https://www.frontiersin.org/articles/10.3389/fmicb.2021.794743/full)
({cite:t}`fuentes2021influence`).


## Structure of the tutorial

````{margin}
```{admonition} Jargon: feature table, feature data
:class: jargon
If terms like *feature table* and *feature data* aren't clear
right now, don't worry! They will be clear by the end of the week.
```
````

This tutorial is split into two parts:

1. **The upstream tutorial** covers steps up to the generation of the feature
   table, which tallys the frequency of amplicon sequence variants (ASV) on a
   per-sample basis, and feature data which lists the sequence that defines
   each ASV in the feature table.
2. **The downstream tutorial** begins with a feature table and feature data and
   constitutes the analysis and interpretation of that information. We'll spend
   the majority of the week on the **downstream** tutorial.

The two parts of this tutorial are both dervived from the same data set
{cite:t}`liao-data-2021`. The **upstream** tutorial uses a relatively small
number of samples (n=41) and is designed to allow us to work through the most
computationally expensive steps of the analysis quickly, so you can get
experience running these steps. By working with fewer
samples, these steps can be run in just a few minutes.

The **downstream** tutorial uses the complete feature table and feature data
published in FigShare by {cite:t}`liao-data-2021`. Since that data set
contains many more samples (n=12,546) and over 550,000,000 sequences, it
would be possible very time-consumming to run the
**upstream** steps on this data interactively. We will show
how to load that data in QIIME 2, and do some filtering of the full data to
focus our work on specific samples of interest. In our case, we'll work with
the {cite:t}`taur-autofmt-2018` samples. However, the full dataset will be
available for you to filter in other ways, and to experiment with on your own.
As the authors note: _These microbiota data, combined
with the curated clinical metadata presented here, can serve as powerful
hypothesis generators for microbiome studies._
