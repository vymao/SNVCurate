## Requirements/expectations :
## - One analysis-ready BAM file for a single sample (as identified in RG:SM)
## - Set of variant calling intervals lists for the scatter, provided in a file
##
## Outputs :
## - One GVCF file and its index
##

# WORKFLOW DEFINITION 
workflow HaplotypeCallerGvcf_GATK4 {
  File input_bam
  File input_bam_index
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  String output_directory
  
  Boolean? make_gvcf
  Boolean making_gvcf = select_first([make_gvcf,true])

  String gatk_docker

  String gatk_path
  
  String sample_basename = basename(input_bam, ".bam")
  
  String vcf_basename = sample_basename

  String output_suffix = if making_gvcf then ".g.vcf.gz" else ".vcf.gz"
  String output_filename = vcf_basename + output_suffix

  Int cores

  call HaplotypeCaller {
    input:
      input_bam = input_bam,
      input_bam_index = input_bam_index,
      output_filename = output_filename,
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      make_gvcf = making_gvcf,
      docker = gatk_docker,
      gatk_path = gatk_path,
      cores = cores
  }
    # Merge per-interval GVCFs
##  call MergeGVCFs {
##    input:
##      input_vcfs = HaplotypeCaller.output_vcf,
##      input_vcfs_indexes = HaplotypeCaller.output_vcf_index,
##      output_filename = output_filename,
##      docker = gatk_docker,
##      gatk_path = gatk_path,
##      output_directory = output_directory
##  }

  call GenotypeGVCFs_single {
    input: 

      input_bam = HaplotypeCaller.output_vcf,
      input_bam_index = HaplotypeCaller.output_vcf_index,
      output_filename = vcf_basename + ".vcf",
      docker = gatk_docker,
      gatk_path = gatk_path,
      output_directory = output_directory,
      ref_fasta = ref_fasta
  }

  # Outputs that will be retained when execution is complete
  output {
    File output_vcf = GenotypeGVCFs_single.output_vcf
    File output_vcf_index = GenotypeGVCFs_single.output_vcf_index
  }
}



# TASK DEFINITIONS

# HaplotypeCaller per-sample in GVCF mode
task HaplotypeCaller {
  File input_bam
  File input_bam_index

  String output_filename
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  Float? contamination
  Boolean make_gvcf

  String gatk_path
  String? java_options
  String java_opt = select_first([java_options, "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10"])

  # Runtime parameters
  String docker
  Int? mem_gb
  Int? disk_space_gb
  Boolean use_ssd = false
  Int? preemptible_attempts

  Int machine_mem_gb = select_first([mem_gb, 7])
  Int command_mem_gb = machine_mem_gb - 1
  Int cores

  command <<<
  set -e
  
    ${gatk_path} --java-options "-Xmx${command_mem_gb}G ${java_opt}" \
      --spark-master local[${cores}] \
      HaplotypeCaller \
      -R ${ref_fasta} \
      -I ${input_bam} \
      -O ${output_filename} \
      -contamination ${default=0 contamination} ${true="-ERC GVCF" false="" make_gvcf}
  >>>

##  runtime {
##    runtime_minutes: ${runtime}
##    cpus: ${cores}
##    requested_memory_mb_per_core: ${memory}
##    queue: ${queue}
##  }

  output {
    File output_vcf = "${output_filename}"
    File output_vcf_index = "${output_filename}.tbi"
  }
}

# Merge GVCFs generated per-interval for the same sample
task MergeGVCFs {
  Array[File] input_vcfs
  Array[File] input_vcfs_indexes
  String output_filename
  String output_directory

  String gatk_path

  # Runtime parameters
  String docker
  Int? mem_gb
  Int? disk_space_gb
  Boolean use_ssd = false
  Int? preemptible_attempts

  Int machine_mem_gb = select_first([mem_gb, 3])
  Int command_mem_gb = machine_mem_gb - 1


  command <<<
  set -e

    ${gatk_path} --java-options "-Xmx${command_mem_gb}G"  \
      MergeVcfs \
      --INPUT ${sep=' --INPUT ' input_vcfs} \
      --OUTPUT ${output_directory}${output_filename}
  >>>

##  runtime {
##    runtime_minutes: ${runtime}
##    cpus: ${cores}
##    requested_memory_mb_per_core: ${memory}
##    queue: ${queue}
##  }


  output {
    File output_vcf = "${output_directory}${output_filename}"
    File output_vcf_index = "${output_directory}${output_filename}.tbi"
  }
}


task GenotypeGVCFs_single {
  File input_bam
  File input_bam_index
  String output_filename
  String output_directory  
  File ref_fasta
  String gatk_path

  # Runtime parameters
  String docker
  Int? mem_gb
  Int? disk_space_gb
  Boolean use_ssd = false
  Int? preemptible_attempts

  Int machine_mem_gb = select_first([mem_gb, 3])
  Int command_mem_gb = machine_mem_gb - 1



  command <<<
  set -e

    ${gatk_path} --java-options "-Xmx${command_mem_gb}G"  \
      GenotypeGVCFs \
      -R ${ref_fasta} \
      -V ${input_bam} \
      -O ${output_filename} \
  >>>


##  runtime {
##    runtime_minutes: ${runtime}
##    cpus: ${cores}
##    requested_memory_mb_per_core: ${memory}
##    queue: ${queue}
##  }


  output {
    File output_vcf = "${output_directory}${output_filename}"
    File output_vcf_index = "${output_directory}${output_filename}.tbi"
  }
}
