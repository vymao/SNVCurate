
# WORKFLOW DEFINITION 
workflow MuSE {
  File input_bam
  File input_bam_index
  File normal_bam
  File normal_bam_index
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  File dbSNP

  String output_directory
  String data_type
 
  String sample_basename = basename(input_bam, ".bam")
  
  String vcf_basename = sample_basename

  String output_filename = vcf_basename + ".vcf"

  call MuSE_call {
    input:
      input_bam = input_bam,
      input_bam_index = input_bam_index,
      normal_bam = normal_bam, 
      normal_bam_index = normal_bam_index,
      output_filename = vcf_basename,
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index
  }

  call MuSE_sump {
    input:
      input_bam = MuSE_call.output_file,
      output_filename = output_filename,
      output_directory = output_directory,
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      dbSNP = dbSNP,
      mode = if (data_type == "WGS") then "-G" else "-E"
  }


  # Outputs that will be retained when execution is complete
  output {
    File output_vcf = MuSE_sump.output_vcf
  }
}



# TASK DEFINITIONS


task MuSE_call {
  File input_bam
  File input_bam_index
  File normal_bam
  File normal_bam_index
  String output_filename
  File ref_dict
  File ref_fasta
  File ref_fasta_index


  # Runtime parameters
  Int? mem_gb
  Int? disk_space_gb
  Boolean use_ssd = false
  Int? preemptible_attempts

  Int machine_mem_gb = select_first([mem_gb, 3])
  Int command_mem_gb = machine_mem_gb - 1

  command <<<
  set -e

    MuSE call \
      -f ${ref_fasta} \
      -O ${output_filename} \
      ${input_bam} \
      ${normal_bam} 
  >>>

  output {
    File output_file = "${output_filename}"
  }
}



task MuSE_sump {
  File input_bam

  String output_filename
  String output_directory  
  String mode

  File ref_dict
  File ref_fasta
  File ref_fasta_index

  File dbSNP

  # Runtime parameters
  Int? mem_gb
  Int? disk_space_gb
  Boolean use_ssd = false
  Int? preemptible_attempts

  Int machine_mem_gb = select_first([mem_gb, 3])
  Int command_mem_gb = machine_mem_gb - 1

  command <<<
  set -e

    MuSE sump \
      -I ${input_bam} \
      -O ${output_filename} \
      -D ${dbSNP} \
      ${mode}
  >>>

  runtime {
    runtime_minutes: 360
    cpus: 1
  }

  output {
    File output_vcf = "${output_directory}${output_filename}"
  }
}
