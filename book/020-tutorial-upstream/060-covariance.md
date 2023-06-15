# Estimating covariance

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

Here we prepare the data for the downstream analysis


````{margin}
```{admonition} Note, compositional data

Explain here different types of transformations exist
```
````


## Calculating the correlation

[//]: # (```{usage})

[//]: # ()
[//]: # (def transformed_table_factory&#40;&#41;:)

[//]: # (    import qiime2)

[//]: # (    )
[//]: # (    transformed_table = qiime2.Artifact.load&#40;'atacama_table_clr', 'data/atacama-table_clr.qza'&#41;)

[//]: # (    return transformed_table)

[//]: # ()
[//]: # (transformed_table = use.init_artifact&#40;'transformed_table', transformed_table_factory&#41;)

[//]: # ()
[//]: # (covariance_matrix = use.action&#40;)

[//]: # (    use.UsageAction&#40;plugin_id='gglasso', action_id='calculate_covariance'&#41;,)

[//]: # (    use.UsageInputs&#40;table=transformed_table,)

[//]: # (                    method='scaled'&#41;,)

[//]: # (    use.UsageOutputNames&#40;covariance_matrix='atacama-table_corr'&#41;)

[//]: # (&#41;)

[//]: # ()
[//]: # (```)

```{note}
We've now reached the end of the **upstream** tutorial. When we begin working
on the **downstream** tutorial, we'll work with the estimated covariance matrix and
feature data artifacts representing covariates.
```


