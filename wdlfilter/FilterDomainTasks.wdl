task FilterFasta {
	String? preCommand
	File fastaIn
	String fastaOutPath
	String scriptPath
	Array[String]+? includeTaxa
	Array[String]+? excludeTaxa
	Array[String]+? includeAccPrefixes
	Array[String]? excludeAccPrefixes
	File? includeAccFile
	File? excludeAccFile
	File? refseqJson
	String? taxDbPath

	Int mem = if (defined(includeTaxa) || defined(excludeTaxa)) then 8 else 4
	Int threads = 1	

	command {
		set -e -o pipefail
		${"mkdir -p $(dirname " + fastaOutPath + ")"}
		${preCommand}
		python ${scriptPath} \
		-i ${fastaIn} \
		-o ${fastaOutPath} \
		${true="--include-accessions " false="" defined(includeAccPrefixes)} ${sep="," includeAccPrefixes} \
		${true="--exclude-accessions " false="" defined(excludeAccPrefixes)} ${sep="," excludeAccPrefixes} \
		${"--include-accessions-from-file " + includeAccFile} \
		${"--exclude-accessions-from-file " + excludeAccFile} \
		${true="--include-taxa " false="" defined(includeTaxa)} ${sep="," includeTaxa} \
		${true="--exclude-taxa " false="" defined(excludeTaxa)} ${sep="," excludeTaxa} \
		${"--dbfile " + taxDbPath} \
		${"-j " + refseqJson}
	}

	output {
		Array[String] stats = read_lines(stdout())
		File? filteredFasta = fastaOutPath
	}

	runtime {
		cpu: threads
		memory: mem
	}
}

task MakeScatterList {
	String dirName
	String pattern	

	command {
		set -e -o pipefail
		find ${dirName} -name "${pattern}"
	}

	output {
		Array[File] fastaFiles = read_lines(stdout())
	}
}

task WriteFilterStats {
	# Gather all stats in one file
	# write_tsv() writes to some tmp dir
	Array[Array[String]] stats
	String tsvFilePath

	command {
		set -e -o pipefail
		cat ${write_tsv(stats)} > ${tsvFilePath}
	}

	output {
		File tsvFile = tsvFilePath
	}
}

task DustmaskFasta {
	File fastaFilePath
	String dustmaskedFilePath
	Boolean? isGzipped = true
	# This needs to be specified
	String dustmaskerExe

	command {
		set -e -o pipefail
		${"mkdir -p $(dirname " + dustmaskedFilePath + ")"}
		zcat ${fastaFilePath} | \
		${dustmaskerExe} -infmt fasta -in - -outfmt fasta | \
		sed '/^>/! s/[^AGCT]/N/g' | \
		gzip -c > ${dustmaskedFilePath}
	}

	output {
		File dustmaskedFile = dustmaskedFilePath
	}
}

task ConcatenateTextFiles {
    Array[File] fileList
    String combinedFilePath
    Boolean unzip = false
    Boolean zip = false

    # When input and output is both compressed decompression is not needed
    String cmdPrefix = if (unzip && !zip) then "zcat " else "cat "
    String cmdSuffix = if (!unzip && zip) then " | gzip -c " else ""

    command {
        set -e -o pipefail
        ${"mkdir -p $(dirname " + combinedFilePath + ")"}
        ${cmdPrefix} ${sep=' ' fileList} ${cmdSuffix} > ${combinedFilePath}
    }

    output {
        File combinedFile = combinedFilePath
    }

    runtime {
        memory: 1
    }
}
