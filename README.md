# ncawareclip
iCLIP and miCLIP mapping pipeline to account for repetitive ncRNAs, especially snRNA, rRNA and tRNA.
Currently available as a Snakemake pipeline. You need demultiplexed CLIP fastq files to begin - the code takes care of the rest.

# Quick start

1. Clone the repo into your local directory.
```
git clone https://github.com/ulelab/ncawareclip.git
```

2. Move into the newly downloaded directory, and create a conda environment with all the dependencies you will require.
```
cd ncawareclip
conda env create -f environment.yml
```

3. Run the test data to make sure everything is working.
```
conda activate ncawareclip
```

3. Create the full annotation/sequence databases you will need, and get your configuration file started using the handy helper tool. Supported genomes are currently: "Hs", - human, "Mm", - mouse", Dm", - Drosophila, "Dr", - Zebrafish, "Rn", - rat, "Sc" - budding yeast (SacCer3) and "Sck1" - budding yeast (SK1, used in meiotic research). Note: you only need to run this once for each species you analyse. Note2: in this pipeline all annotation origins and processing code is all in this repo, so certain files could be replaced if you wanted.
```
```

4. Edit `config.yaml` to provide paths to your demultiplexed fastq files and path where you are running the pipeline. For UMI removal make sure that the random barcode is moved to the fastq header as "RBC:NNNNN". I would highly recommend the speedy demultiplexer [Ultraplex](https://github.com/ulelab/ultraplex) for completely unbiased reasons.
```
```

5. Run the Snakemake pipeline.
```
```

# The output

# Detailed pipeline options

# Notes on tRNAs

# Notes on snRNAs

# Notes on RepBase

# Publications using this pipeline

For any questions please feel free to raise a GitHub issue or contact me at charlotte.capitanchik@crick.ac.uk.
