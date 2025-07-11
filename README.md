# High-dimensional Statistics with QIIME 2

[![Github Pages](https://img.shields.io/badge/github%20pages-121013?style=for-the-badge&logo=github&logoColor=white)](https://github.com/bio-datascience/atacama-soil-microbiome-tutorial)

Welcome to the **High-dimensional Statistics with QIIME 2** documentation!  
This site features tutorials and usage guides for the [`q2-gglasso`](https://github.com/bio-datascience/q2-gglasso) and [`q2-classo`](https://github.com/bio-datascience/q2-classo) plugins, which bring advanced high-dimensional statistical modeling to microbiome analysis within the QIIME 2 framework.

---

## 🚀 Building the Documentation

**Create and activate conda environment:**
```shell
mamba create -n jupyter-book -c conda-forge jupyter-book

conda activate jupyter-book
```

**Build the book:**

```shell
jupyter-book build --all moshpit_docs
```