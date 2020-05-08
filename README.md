# SNVCurate
A pipeline for producing a callset of somatic mutations. It also provides a method to curate targeted capture regions in the event that this BED file is not available. 
This pipeline is used only for running on Orchestra, the Harvard Medical School cluster, with the SLURM job scheduler. This pipeline also utilizes GATK 4.1.2.0, the latest version of GATK available on the cluster.

**Note: This pipeline presumes that the BAM files are formatted properly (demultiplexed and indexed correctly). If no matched normal is available, run this pipeline using another normal of the same sequencing type.**

The pipeline is split into two parts: one to generate callsets from Mutect2 (and MuSE + HaplotypeCaller if normals are available), and one to filter those calls. Callset generation is automated on the cluster using the Cromwell execution engine and customized .wdl scripts produced from the datasets inputted. Post-processing is automated using bash wrapper scripts utilizing various software pre-loaded in the environment. It is split this way to prevent overloading of the job partitions, so you can continue submitting jobs as they complete.

If you have some samples with matched normals and some without that you would like to call, you should create separate .csv files for each of these datasets and run this pipeline on each set individually. 

**Please see the [wiki](https://github.com/vymao/SNVCurate/wiki)** for instructions on how to use this pipeline.
