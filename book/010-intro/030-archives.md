---
jupytext:
  cell_metadata_filter: -all
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.12
    jupytext_version: 1.9.1
kernelspec:
  display_name: calysto_bash
  language: calysto_bash
  name: calysto_bash
---

# QIIME 2 archives

One of the first things that new QIIME 2 users often notice is the `.qza` and
`.qzv` files that QIIME 2 uses. All files generated by QIIME 2 are either `.qza`
or `.qzv` files, and these are simply zip files that store your data alongside
some QIIME 2-specific metadata. You can unzip these files with any typical unzip
utility, such as WinZip, 7Zip, or unzip, and you don't need to have QIIME 2
installed to do that.

For example, let's download a `.qza` file and take a quick look.

```{code-cell}
curl -sL \
  "https://docs.qiime2.org/2020.11/data/tutorials/moving-pictures/rep-seqs.qza" > \
  "sequences.qza"
```

```{code-cell}
unzip sequences.qza
```

Notice that we haven't used any QIIME 2 commands so far. We downloaded a file,
and then unzipped it as we would any zip file. Allowing users to access their
data _without_ QIIME 2 was one of the earliest design goals of the system. This
ensures that if QIIME 2 isn't available to you for some reason, you can still
access any data that you generated with QIIME 2.

If you look through the list of files that were created by the unzip command
above, you'll see there is a top-level directory in the output with a
crazy-looking name. This directory contains copies of all of the files and
directories from `sequences.qza`. Within this crazy-named directory the `data`
directory contains a single file, `sequences.fasta`, which contains (you guessed
it!) sequence data in fasta format.  If, for example, you're interested in
getting your sequence data out of QIIME 2 to analyze it with another program,
you can unzip your `.qza` file, and use the `sequences.fasta` file for what ever
you need to do with it.

You might wonder why we bothered with having QIIME 2 create these zip files in
the first place, rather than just have it use the typical file formats like
fasta, newick, biom, and so on. That has to do with the other stored in the zip
file. The other files in the zip file are not intended to be viewed by a human,
and you don't need any of them to work with the file (or files) in the `data`
directory. QIIME 2 uses the information in the `provenance` directory to record
data provenance, helping you to ensure that your bioinformatics work will be
reproducible. It also stores a unique identifier for the data in the
`metadata.yaml` file, which facilitates data management. That unique identifier
is also the name of the directory that contains all of the files from the `.qza`
when you unzip it. The `metadata.yaml` file also identifies the semantic type of
the data (we'll come back to that shortly). This other information empowers your
bioinformatics work in ways that will be more clear as you get further along in
your learning. We'll revisit this topic throughout the book.

The `.qza` file extension is an abbreviation for QIIME Zipped Artifact, and the
`.qzv` file extension is an abbreviation for QIIME Zipped Visualization. `.qza`
files (which are often simply referred to as _artifacts_) are intermediary files
in a QIIME 2 analysis, usually containing raw data of some sort. These files are
generated by QIIME 2 and are intended to be consumed by QIIME 2. `.qzv`` files
(which are often simply referred to as _visualizations_) are terminal results in
a QIIME 2 analysis, such as an interactive figure or the results of a
statistical test. These files are generated by QIIME 2 and are intended to be
consumed by a human.

QIIME 2 also provides some of its own utilities for getting data out of `.qza`
and `.qzv` files. If you're working with the QIIME 2 command-line interface
(which we'll use a lot in this book), the most relevant command is `qiime tools
export`. If you were to run this on the `.qza` file we downloaded above, you'd
see the following:

```{code-cell}
qiime tools export --input-path sequences.qza --output-path exported-sequences/
```

This command will unzip the archive, and take all of the files from the `data`
directory and place them in `exported-sequences`. Thus if you do have QIIME 2
installed, you can get your data out of a QIIME 2 artifact without all of the
QIIME 2-specific metadata using this command.

```{admonition} Jargon: Confused by the term "artifact"?
It has been brought to our attention that the term _artifact_ can be confusing, since it is often used in science to indicate a feature that is not present naturally in a system but rather observed as a result of some technical aspect of studying that system. For example, homopolyer runs such as the `A`s in `ACTGTACTAAAAAAAAAAATGCACGTGAC` were commonly reported by some early sequencing instruments to be longer then they were in nature due to the way the sequencing reaction worked. In QIIME 2 we use the definition of an artifact as an object that was created by some process, like an archaeological artifact. This is common in data science, and we didn't realize the potential for confusion until we were a little too far along to easily change the name.
```