# Gene Analysis of Schizophrenia Patient Samples

![](https://images.unsplash.com/photo-1716833322865-56bae681995c?q=80&w=1929&auto=format&fit=crop&ixlib=rb-4.0.3&ixid=M3wxMjA3fDB8MHxwaG90by1wYWdlfHx8fGVufDB8fHx8fA%3D%3D)

This is a repository for the personal project in the Pavlidis Lab under Nairuz Elazzabi, involving scRNAseq datasets of schizophrenia datasets.

To read on what has been done so far for the project, please read the report under the `reports/` folder. It is available as a PDF (viewable on GitHub) and as HTML (needs to be downloaded first).

Please note that some of the code is (in the works of being) translated into an R package. For more information, please see [curatoR](https://github.com/xnrxng/curatoR).

To reproduce the analysis, please follow these steps:

1.  Open your terminal and run the following:

    `git clone git@github.com:xnrxng/PavlabSZProject.git`

2.  Make sure you have the latest changes with:

    `git pull`

3.  Then, run:

    `cd data/data_processed`

4.  To get the raw data, run:

    `wget https://brainscope.gersteinlab.org/data/snrna_expr_matrices_zip/CMC.zip`

    `wget https://brainscope.gersteinlab.org/data/snrna_expr_matrices_zip/MultiomeBrain.zip`

    `wget https://brainscope.gersteinlab.org/data/snrna_expr_matrices_zip/SZBDMulti-Seq.zip`

    `wget https://zenodo.org/records/6921620/files/snRNA-seq_and_spatial_transcriptomics.zip?download=1`

5.  Unzip each folder using this code. You can delete the zipped files if you want after.

    `for all in *.gz; do gunzip $all; done`

6.  Navigate into each folder using the `cd` command and unzip all the files using:

    `for all in *.gz; do gunzip $all; done`

7.  For the folder called "snRNA-seq_and_spatial_transcriptomics", rename it to "Batiuk". You can delete all the files in this folder except for the ones named "snRNA-seq_raw_countmatrices.RDS" and "annotations_final.RDS".

8.  Navigate back to the root of the repository using `cd` .

9.  To clean up all the made artifacts, run:

    `make clean`

10. To re-make all the artifacts, run:

    `make all`

NOTE: Re-making all the artifacts at once will take a tremendous amount of memory and it is generally not recommended. It is recommended that the scripts be run in chunks instead. Please use it with caution.
