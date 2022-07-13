# Exploring the tutorial metadata

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

## Access the study metadata

In this chapter we'll begin our work with QIIME 2 and the tutorial data. We'll
start by downloading the metadata, generating a summary of it, and exploring
that summary.

First, download the metadata.

```{usage}
def read_metadata():
    ## This function is identical to the filter.md metadata_factory function - should
    ## be able to call that once issue#1 is addressed.
    import tempfile
    import requests
    import pandas as pd
    import numpy as np

    import qiime2

    sample_metadata_url = 'https://data.qiime2.org/2022.2/tutorials/atacama-soils/sample_metadata.tsv'
    data = requests.get(sample_metadata_url)
    with tempfile.NamedTemporaryFile() as f:
        f.write(data.content)
        sample_metadata = pd.read_csv(f.name, sep='\t', index_col=0)

    return qiime2.Metadata(sample_metadata)


sample_metadata = use.init_metadata('sample_metadata', read_metadata)
```

## View the metadata

Next, we'll get a view of the study metadata as QIIME 2 sees it. This
will allow you to assess whether the metadata that QIIME 2 is using is as you
expect. You can do this using the `tabulate` action in QIIME 2's `q2-metadata`
plugin as follows.

```{usage}
use.action(
    use.UsageAction(plugin_id='metadata', action_id='tabulate'),
    use.UsageInputs(input=sample_metadata),
    use.UsageOutputNames(visualization='metadata_summ_1')
)
```

Spend a few minutes now exploring the Galaxy environment on your own, and
exploring the metadata that we'll use in this tutorial. If you have questions
about how to use Galaxy or QIIME 2 View, this is a great time to ask those
questions.
