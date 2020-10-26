## Original scripts written by Dr. Samuel W. Kazer of Ordovas-Montañes Lab for use by Lingwood Lab on a HPC
## Scripts wrapped in WDL by James Gatter of Shalek Lab

version 1.0

workflow amplicon_BCR_alignment {
	input {
		File sample_sheet
		
		String? blast_directory
		String species_heavy = "human"
		String receptor_gene_chain_heavy = "IGH"
		String species_light = "mouse"
		String receptor_gene_chain_light = "IGK,IGL"

		String pandaseq_docker = "pipecraft/pandaseq:2.11"
		String migmap_docker = "gadeta/migmap:V4"
		String consensus_docker = "shaleklab/ampliconbcr:latest"
		String zones = "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
		Int preemptible = 2
	}

	scatter (sample in read_objects(sample_sheet)) {
		call panda_express as panda_express_heavy {
			input:
				sample_name = sample.Sample,
				fastq_R1 = sample.HCR1,
				fastq_R2 = sample.HCR2,
				docker = pandaseq_docker,
				zones = zones,
				preemptible = preemptible,
		}
		call panda_express as panda_express_light {
			input:
				sample_name = sample.Sample,
				fastq_R1 = sample.LCR1,
				fastq_R2 = sample.LCR2,
				docker = pandaseq_docker,
				zones = zones,
				preemptible = preemptible,
		}
	}
	Array[Pair[String,File]] sample_paired_fastqs_heavy = panda_express_heavy.sample_paired_fastq
	Array[Pair[String,File]] sample_paired_fastqs_light = panda_express_light.sample_paired_fastq

	scatter (sample_paired_fastq in sample_paired_fastqs_heavy) {
		call migmap as migmap_heavy {
			input:
				sample_name = sample_paired_fastq.left,
				paired_fastq = sample_paired_fastq.right,
				blast_directory = blast_directory,
				species = "human",
				receptor_gene_chain = "IGH",
				docker = migmap_docker,
				preemptible = preemptible,
				zones = zones
		}
	}
	scatter (sample_paired_fastq in sample_paired_fastqs_light) {
		call migmap as migmap_light {
			input:
				sample_name = sample_paired_fastq.left,
				paired_fastq = sample_paired_fastq.right,
				blast_directory = blast_directory,
				species = "mouse",
				receptor_gene_chain = "IGK,IGL",
				docker = migmap_docker,
				preemptible = preemptible,
				zones = zones
		}
	}

	call consensus {
		input:
			aligns_corrected_heavy = migmap_heavy.align_corrected,
			aligns_corrected_light = migmap_light.align_corrected,
			docker = consensus_docker,
			preemptible = preemptible,
			zones = zones
	}

	output {
		File consensus_alignments = consensus.consensus_alignments
		Array[File] aligns_corrected_heavy = migmap_heavy.align_corrected
		Array[File] aligns_corrected_light = migmap_light.align_corrected
	}
}

task panda_express {
	input {
		String sample_name
		File fastq_R1
		File fastq_R2
		Int minimum_read_length = 200
		
		Int preemptible
		String zones
		String docker
		Int number_cpu_threads = 4
		Int task_memory_GB = 16
		String disks = "local-disk 128 HDD"
		Int boot_disk_size_GB = 12
	}
	String paired_fastq = sample_name + "_paired.fastq"
	String log  = sample_name + "_pandaseq.log"
	command <<<
		set -e
		mkdir -p /cromwell_root/out/
		cd /cromwell_root/out/
		pandaseq -F \
			-f ~{fastq_R1} \
			-r ~{fastq_R2} \
			-w ~{paired_fastq} \
			-g ~{log} \
			-l ~{minimum_read_length} \
			-T ~{number_cpu_threads}
		ls
		echo "Panda-Seq is done, produced ~{paired_fastq}. Read ~{log} for more details."
	>>>
	output {
		Pair[String, File] sample_paired_fastq = (sample_name, "/cromwell_root/out/~{paired_fastq}")
		File pandaseq_log = "/cromwell_root/out/~{log}"
	}
	runtime {
		docker: docker
		preemptible: preemptible
		memory: "~{task_memory_GB}G"
		zones: zones
		bootDiskSizeGb: boot_disk_size_GB
		disks: disks
		cpu: number_cpu_threads
	} 
}

task migmap {
	input {
		String sample_name
		File paired_fastq
		String? blast_directory
		String species
		String receptor_gene_chain

		Int preemptible
		String zones
		String docker
		Int number_cpu_threads = 4
		Int task_memory_GB = 32
		String disks = "local-disk 128 HDD"
		Int boot_disk_size_GB = 12
	}
	String align = sample_name + "_align.txt"
	String corrected = sample_name + "_corr_align.txt"
	command <<<
		set -euo pipefail
		mkdir -p /cromwell_root/out/
		cd /cromwell_root/out/
		java -Xmx16G -jar /opt/migmap-1.0.3/migmap-1.0.3.jar \
			--blast-dir ~{default="/opt/ncbi-igblast-1.4.0/bin/" blast_directory} \
			-R ~{receptor_gene_chain} \
			-S ~{species} \
			-p ~{number_cpu_threads} \
			~{paired_fastq} \
			~{align} \
			
		java -Xmx16G -cp /opt/migmap-1.0.3/migmap-1.0.3.jar com.antigenomics.migmap.Correct \
			~{align} \
			~{corrected}
	>>>
	output {
		File align_corrected = "/cromwell_root/out/~{corrected}"
	}
	runtime {
		docker: docker
		preemptible: preemptible
		memory: "~{task_memory_GB}G"
		zones: zones
		bootDiskSizeGb: boot_disk_size_GB
		disks: disks
		cpu: number_cpu_threads
	}
}

task consensus {
	input {
		Array[File] aligns_corrected_heavy
		Array[File] aligns_corrected_light

		Int preemptible
		String zones
		String docker
		Int number_cpu_threads = 1
		Int task_memory_GB = 4
		String disks = "local-disk 128 HDD"
		Int boot_disk_size_GB = 12
	}
	command <<<
		set -e
		mkdir -p /cromwell_root/out/
		cd /cromwell_root/out/
		python /scripts/concatenate_chain_results.py \
			--heavy ~{sep=',' aligns_corrected_heavy} \
			--light ~{sep=',' aligns_corrected_light} \
			--output-path /cromwell_root/out/
		python /scripts/collect_migmap_results_BCR.py
	>>>
	output {
		File consensus_alignments = "/cromwell_root/out/consensus_alignments.txt"
	}
	runtime {
		docker: docker
		preemptible: preemptible
		memory: "~{task_memory_GB}G"
		zones: zones
		bootDiskSizeGb: boot_disk_size_GB
		disks: disks
		cpu: number_cpu_threads
	}
}