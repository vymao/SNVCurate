## Requirements/expectations :
## - One analysis-ready BAM file for a single sample (as identified in RG:SM)
## - Set of variant calling intervals lists for the scatter, provided in a file
##
## Outputs :
## - One GVCF file and its index
##

# WORKFLOW DEFINITION 
workflow MuTecT {
  File input_bam
  File input_bam_index
  File? normal_bam
  File? normal_bam_index
  File? panel
  File? panel_index
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  File interval_list
  File gnomad
  File gnomad_index
  String output_directory
  String mode

  Boolean? make_gvcf
  Boolean making_gvcf = select_first([make_gvcf,true])

  Array[File] scattered_calling_intervals = read_lines(interval_list)

  String gatk_docker

  String gatk_path
  
  String sample_basename = basename(input_bam, ".bam")
  
  String vcf_basename = sample_basename

  String output_suffix = ".vcf"
  String output_filename = vcf_basename + output_suffix


  if (mode == "normal") {
    scatter (interval_file in scattered_calling_intervals) {
        call MuTecT_normal {
          input:
            input_bam = input_bam,
            input_bam_index = input_bam_index,
            normal_bam = normal_bam,
            normal_bam_index = normal_bam_index,
            output_filename = output_filename,
            regions_list = interval_file,
            ref_dict = ref_dict,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            docker = gatk_docker,
            gatk_path = gatk_path, 
            gnomad = gnomad,
            gnomad_index = gnomad_index,
            regions_list = interval_file
        }
    }

  }


  if (mode == "panel") {
    scatter (interval_file in scattered_calling_intervals) {
        call MuTecT_PoN {
          input:
            input_bam = input_bam,
            input_bam_index = input_bam_index,
            panel = panel,
            panel_index = panel_index,
            output_filename = output_filename,
            regions_list = interval_file,
            ref_dict = ref_dict,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            docker = gatk_docker,
            gatk_path = gatk_path, 
            gnomad = gnomad,
       	    gnomad_index = gnomad_index,
            regions_list = interval_file
        }
    }
  }

  call MergeMutectStats {
    input:
      input_stats = if (mode == "normal") then MuTecT_normal.output_stats else MuTecT_PoN.output_stats,
      output_filename = output_filename,
      gatk_path = gatk_path
  }


  call MergeVCFs {
    input:
      input_vcfs = if (mode == "normal") then MuTecT_normal.output_vcf else MuTecT_PoN.output_vcf,
      input_vcfs_indexes = if (mode == "normal") then MuTecT_normal.output_vcf_index else MuTecT_PoN.output_vcf_index,
      output_filename = output_filename,
      output_directory = output_directory

  }

  call FilterMutectCalls {
    input:
      input_vcf = MergeVCFs.output_vcf,
      output_filename = output_filename,
      gatk_path = gatk_path,
      input_stats = MergeMutectStats.output_stats,
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,  
      output_directory = output_directory
  }

  output {
    File output_vcf = FilterMutectCalls.output_vcf
  }
}



# TASK DEFINITIONS
task MuTecT_normal {
  File input_bam
  File input_bam_index 
  File normal_bam
  File normal_bam_index 
  File regions_list
  File gnomad
  File gnomad_index

  String output_filename
  File ref_dict
  File ref_fasta
  File ref_fasta_index

  String gatk_path
  String normal_name = basename(normal_bam, ".bam")

  # Runtime parameters
  String docker
  Int? mem_gb
  Int? disk_space_gb
  Boolean use_ssd = false
  Int? preemptible_attempts

  Int machine_mem_gb = select_first([mem_gb, 7])
  Int command_mem_gb = machine_mem_gb - 1

  command <<<
  set -e
  
    ${gatk_path} \
      Mutect2 \
      -R ${ref_fasta} \
      -I ${input_bam} \
      -I ${normal_bam} \
      -normal ${normal_name} \
      --germline-resource ${gnomad} \
      -L ${regions_list} \
      -O ${output_filename} 
  >>>

  output {
    File output_vcf = "${output_filename}"
    File output_vcf_index = "${output_filename}.idx"
    File output_stats = "${output_filename}.stats"
  }
}

task MuTecT_PoN {
  File input_bam
  File input_bam_index 
  File panel
  File panel_index 
  File regions_list
  File gnomad
  File gnomad_index

  String output_filename
  File ref_dict
  File ref_fasta
  File ref_fasta_index

  String gatk_path

  # Runtime parameters
  String docker
  Int? mem_gb
  Int? disk_space_gb
  Boolean use_ssd = false
  Int? preemptible_attempts

  Int machine_mem_gb = select_first([mem_gb, 7])
  Int command_mem_gb = machine_mem_gb - 1

  command <<<
  set -e
  
    ${gatk_path} \
      Mutect2 \
      -R ${ref_fasta} \
      -I ${input_bam} \
      --panel-of-normal ${panel} \
      --germline-resource ${gnomad} \
      -L ${regions_list} \
      -O ${output_filename} 
  >>>

  output {
    File output_vcf = "${output_filename}"
    File output_vcf_index = "${output_filename}.idx"
    File output_stats = "${output_filename}.stats"
  }
}

task MergeVCFs {
  Array[File] input_vcfs
  Array[File] input_vcfs_indexes
  String output_filename
  String output_directory

  # Runtime parameters
  Int? mem_gb
  Int? disk_space_gb
  Boolean use_ssd = false
  Int? preemptible_attempts

  Int machine_mem_gb = select_first([mem_gb, 3])
  Int command_mem_gb = machine_mem_gb - 1


  command <<<
  set -e

    bcftools concat ${sep=' ' input_vcfs} -o ${output_filename}
  >>>


  output {
    File output_vcf = "${output_filename}"
  }
}

task MergeMutectStats {
  Array[File] input_stats
  String output_filename

  String gatk_path

  # Runtime parameters
  Int? mem_gb
  Int? disk_space_gb
  Boolean use_ssd = false
  Int? preemptible_attempts

  Int machine_mem_gb = select_first([mem_gb, 3])
  Int command_mem_gb = machine_mem_gb - 1


  command <<<
  set -e

    ${gatk_path} \
      MergeMutectStats \
      -stats ${sep=' -stats ' input_stats} \
      -O ${output_filename}.stats
  >>>

##  runtime {
##    runtime_minutes: ${runtime}
##    cpus: ${cores}
##    requested_memory_mb_per_core: ${memory}
##    queue: ${queue}
##  }


  output {
    File output_stats = "${output_filename}.stats"
  }
}

task FilterMutectCalls {
  File input_vcf
  File input_stats

  String output_filename
  File ref_dict
  File ref_fasta
  File ref_fasta_index

  String gatk_path
  String output_directory

  # Runtime parameters
  Int? mem_gb
  Int? disk_space_gb
  Boolean use_ssd = false
  Int? preemptible_attempts

  Int machine_mem_gb = select_first([mem_gb, 7])
  Int command_mem_gb = machine_mem_gb - 1

  command <<<
  set -e

    ${gatk_path} \
      FilterMutectCalls \
      -R ${ref_fasta} \
      -V ${input_vcf} \
      --stats ${input_stats} \
      -O ${output_directory}${output_filename}
  >>>

  output {
    File output_vcf = "${output_directory}${output_filename}"
  }
}