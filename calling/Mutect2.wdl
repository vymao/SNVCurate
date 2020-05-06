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
  String normal_name
  String parallel
  File path2picard

  Array[File] scattered_calling_intervals = read_lines(interval_list)

  String gatk_docker

  String gatk_path
  
  String sample_basename = basename(input_bam, ".bam")
  
  String vcf_basename = sample_basename

  String output_suffix = ".vcf"
  String output_filename = vcf_basename + output_suffix


  if (mode == "normal") {
    if (parallel  == "True") {
      scatter (interval_file in scattered_calling_intervals) {
        call MuTecT_normal {
          input:
            input_bam = input_bam,
            input_bam_index = input_bam_index,
            normal_bam = normal_bam,
            normal_bam_index = normal_bam_index,
            output_filename = output_filename,
            ref_dict = ref_dict,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            docker = gatk_docker,
            gatk_path = gatk_path, 
            gnomad = gnomad,
            gnomad_index = gnomad_index,
            regions_list = interval_file,
            normal_name = normal_name
        }

        call bgZipVCFs {
          input:
            input_vcf = MuTecT_normal.output_vcf,
            input_vcf_index = MuTecT_normal.output_vcf_index,
            output_filename = output_filename,
            output_basename = vcf_basename,
            output_directory = output_directory
        }
      }
    } 

    if (parallel == "False") {
      call MuTecT_single {
        input:
          input_bam = input_bam,
          input_bam_index = input_bam_index,
          normal_bam = normal_bam,
          normal_bam_index = normal_bam_index,
          output_filename = output_filename,
          ref_dict = ref_dict,
          ref_fasta = ref_fasta,
          ref_fasta_index = ref_fasta_index,
          docker = gatk_docker,
          gatk_path = gatk_path, 
          gnomad = gnomad,
          gnomad_index = gnomad_index,
          normal_name = normal_name
      }

      call LearnReadOrientationModel  {
        gatk_path = gatk_path,
        artifact = MuTecT_single.artifacts
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

  if (parallel == "True") {
    call MergeMutectStats {
      input:
        input_stats = if (mode == "normal") then MuTecT_normal.output_stats else MuTecT_PoN.output_stats,
        output_filename = output_filename,
        gatk_path = gatk_path
    }

    call MergeVCFs { 
      input:
        input_vcfs = if (mode == "normal") then bgZipVCFs.output_vcf else MuTecT_PoN.output_vcf,
        input_vcfs_indexes = if (mode == "normal") then bgZipVCFs.output_vcf_index else MuTecT_PoN.output_vcf_index,
        output_filename = output_filename,
        output_basename = vcf_basename,
        output_directory = output_directory,
        picard = path2picard
      
    }

    call LearnReadOrientationModel {
      gatk_path = gatk_path,
      artifacts = MuTecT_normal.artifacts
    }
  }

  call FilterMutectCalls {
    input:
      input_vcf = if (parallel == "True") then MergeVCFs.output_vcf else MuTecT_single.output_vcf,
      output_filename = output_filename,
      gatk_path = gatk_path,
      input_stats = if (parallel == "True") then MergeMutectStats.output_stats else MuTecT_single.output_stats,
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      input_priors = LearnReadOrientationModel.output_artifacts,
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
  String normal_name

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
      -I ${normal_bam} \
      -normal ${normal_name} \
      --germline-resource ${gnomad} \
      -L ${regions_list} \
      --f1r2-tar-gz ${output_filename}.f1r2.tar.gz \
      -O ${output_filename} 
  >>>

  output {
    File output_vcf = "${output_filename}"
    File output_vcf_index = "${output_filename}.idx"
    File output_stats = "${output_filename}.stats"
    File artifacts = "${output_filename}.f1r2.tar.gz"
  }
}

task MuTecT_single {
  File input_bam
  File input_bam_index 
  File normal_bam
  File normal_bam_index 
  File gnomad
  File gnomad_index
  String normal_name

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
      -I ${normal_bam} \
      -normal ${normal_name} \
      --germline-resource ${gnomad} \
      --f1r2-tar-gz f1r2.tar.gz \
      -O ${output_filename} 
  >>>

  output {
    File output_vcf = "${output_filename}"
    File output_vcf_index = "${output_filename}.idx"
    File output_stats = "${output_filename}.stats"
    File artifacts = "f1r2.tar.gz"
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

task LearnReadOrientationModel {
  String gatk_path
  File? artifact
  Array[File]? artifacts
  String inputs = select_first([artifact, ${sep=' ' artifacts}, "default"])

  command <<<
  set -e
  
    ${gatk_path} \
      LearnReadOrientationModel \
      ${inputs} \
      -O read-orientation-model.tar.gz
  >>>

  output {
    File output_artifacts = "read-orientation-model.tar.gz"
  }
}

task MergeVCFs {
  Array[File] input_vcfs
  Array[File] input_vcfs_indexes
  String output_filename
  String output_basename
  String output_directory
  File picard

  # Runtime parameters
  Int? mem_gb
  Int? disk_space_gb
  Boolean use_ssd = false
  Int? preemptible_attempts

  Int machine_mem_gb = select_first([mem_gb, 3])
  Int command_mem_gb = machine_mem_gb - 1


  command <<<
  set -e
    ##module load bcftools
    ##bcftools concat -a ${sep=' ' input_vcfs} -o ${output_filename}
    java -jar ${picard} MergeVcfs I=${sep=' I=' input_vcfs} O=${output_filename} 
  >>>


  output {
    File output_vcf = "${output_filename}"
  }
}

task bgZipVCFs {
  File input_vcf
  File input_vcf_index
  String output_filename
  String output_basename
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
    module load bcftools
    bgzip -c ${input_vcf} > ${output_filename}.gz
    bcftools index ${output_filename}.gz

  >>>


  output {
    File output_vcf = "${output_filename}.gz"
    File output_vcf_index = "${output_filename}.gz.csi"
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
  File input_priors

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
      --ob-priors ${input_priors} \
      -O ${output_directory}${output_filename}
  >>>

  output {
    File output_vcf = "${output_directory}${output_filename}"
  }
}
