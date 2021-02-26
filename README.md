# CUTnTag Analysis pipeline scripts


R and bash scripts to analyse CUT&Tag data using Signac, Deeptools, MACS2, SEACR and more.

[Link to Google docs of bulleted pipeline steps](https://docs.google.com/document/d/1l6QbU-BfqtGiNfMNGywC2oEfkmou2W12KBOn7L1Kq6M/edit?usp=sharing)

This folder contains the "mother" scripts, which have then been copied to the analysis folders such as `/H3K27ac_VMH_ExptD` and amended as appropriate to filter and set ie clustering resolution.
Note that in current state it's necessary to set the path to ie. fastq and hashtag fastq in both 00_parameters and in each bash script.


`cutntag.yml` contains conda environment with all the necessary packages.

I found it useful to invoke the scripts in this manner:

`R -e "rmarkdown::render('03_demultiplex.Rmd',output_file='03_demultiplex.html')"`
`bash -x 04_filteredCells_makeBAMmakeBWfindPeaks.sh &> logfile_04.log`

[Link to write-up about cutNtag in general and some thoughts on it.](https://docs.google.com/document/d/14fgKkcWUS6FbOAhBeu6TpAHERR1h5MOPjP49Kbv4gmY/edit)
