version 1.0

workflow Mocha {

  input {
    File phased_bcf
    File phased_bcf_idx

    String rule

    File excluded_variants
    File excluded_variants_idx

    File cnp
  }

  call Mocha {
    input:
      bcf = phased_bcf,
      bcf_idx = phased_bcf_idx,

      rule = rule,

      excluded_variants = excluded_variants,
      excluded_variants_idx = excluded_variants_idx,

      cnp = cnp

  }

  call SummaryPlot {
    input:
      genome_stats = Mocha.genome_stats,
      mosaic_calls = Mocha.mosaic_calls      
  }

  output {
    File mocha_bcf = Mocha.out_bcf
    File mocha_bcf_idx = Mocha.out_bcf_idx
    File mosaic_calls = Mocha.mosaic_calls
    File genome_stats = Mocha.genome_stats
    File ucsc_bed = Mocha.ucsc_bed
    File summary_plot = SummaryPlot.summary_plot
  }
}

task Mocha {
  input {
    File bcf
    File bcf_idx

    String rule

    File excluded_variants
    File excluded_variants_idx

    File cnp

    Int? disk_size_override
    String docker = "cwhelan/mocha:v1.0"
    Int threads = 1
    Int memory = 2
    Int preemptible_attempts = 3

  }

  String filebase = basename(bcf, ".bcf")
  
  Int bcf_size = ceil(size(bcf, "GiB"))  
  Int excluded_variants_size = ceil(size(excluded_variants, "GiB"))  
  Int cnp_size = ceil(size(cnp, "GiB"))  

  Int disk_size = select_first([disk_size_override, 10 + excluded_variants_size + cnp_size + bcf_size * 2])


  output {
    File out_bcf = "~{filebase}.mocha.bcf"
    File out_bcf_idx = "~{filebase}.mocha.bcf.csi"
    File mosaic_calls = "~{filebase}.mocha.tsv"
    File genome_stats = "~{filebase}.stats.tsv"
    File ucsc_bed = "~{filebase}.ucsc.bed"
  }

  command <<<

    set -euo pipefail

    bcftools +mocha \
      --rules ~{rule} \
      --no-version \
      --output-type b \
      --output ~{filebase}.mocha.bcf \
      --threads ~{threads} \
      --variants ^~{excluded_variants} \
      --mosaic-calls ~{filebase}.mocha.tsv \
      --genome-stats ~{filebase}.stats.tsv \
      --ucsc-bed ~{filebase}.ucsc.bed \
      --cnp ~{cnp} \
      ~{bcf} 
    bcftools index -f ~{filebase}.mocha.bcf

  >>>
  runtime {
    memory: memory + " GiB"
    disks: "local-disk " + disk_size + " HDD"
    cpu: threads
    docker: docker
  }
}


task SummaryPlot {
  input {
    File genome_stats
    File mosaic_calls

    Int disk_size = 10
    String docker = "cwhelan/mocha:v1.0"
    Int threads = 1
    Int memory = 2
    Int preemptible_attempts = 3

  }

  String filebase = basename(mosaic_calls, ".mocha.tsv")

  output {
    File summary_plot = "~{filebase}.summary.pdf"
  }

  command <<<

    set -euo pipefail

    summary_plot.R \
      --pdf ~{filebase}.summary.pdf \
      --stats ~{genome_stats} \
      --calls ~{mosaic_calls}

  >>>
  runtime {
    memory: memory + " GiB"
    disks: "local-disk " + disk_size + " HDD"
    cpu: threads
    docker: docker
  }
}
