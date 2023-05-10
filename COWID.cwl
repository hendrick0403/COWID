{
    "class": "Workflow",
    "cwlVersion": "v1.2",
    "label": "COWID",
    "$namespaces": {
        "sbg": "https://sevenbridges.com"
    },
    "inputs": [
        {
            "id": "centrifuge_index_archive",
            "sbg:fileTypes": "TAR, TAR.GZ",
            "type": "File",
            "label": "TAR index",
            "doc": "The basename of the index for the reference genomes. The basename is the name of any of the index files up to but not including the final .1.cf / etc. centrifuge looks for the specified index first in the current directory, then in the directory specified in the CENTRIFUGE_INDEXES environment variable.",
            "sbg:x": -521.8116455078125,
            "sbg:y": 631.5410766601562
        },
        {
            "id": "input_file",
            "sbg:fileTypes": "FASTA, FASTQ, FA, FQ, FASTQ.GZ, FQ.GZ, FASTQ.BZ2, FQ.BZ2",
            "type": "File[]",
            "label": "FASTQ reads",
            "doc": "Read sequence in FASTQ or FASTA format. Could be also gzip'ed (extension .gz) or bzip2'ed (extension .bz2). In case of paired-end alignment it is crucial to set metadata 'paired-end' field to 1/2.",
            "sbg:x": -520.5198364257812,
            "sbg:y": 447.8516845703125
        },
        {
            "id": "reference",
            "sbg:fileTypes": "FASTA",
            "type": "File",
            "label": "FASTA reference",
            "doc": "Reference sequence in fasta format.",
            "sbg:x": -515.7294311523438,
            "sbg:y": -2.3527281284332275
        }
    ],
    "outputs": [
        {
            "id": "output_file",
            "outputSource": [
                "bcftools_consensus/output_file"
            ],
            "sbg:fileTypes": "FA",
            "type": "File?",
            "label": "FASTA genome",
            "doc": "Consensus sequence.",
            "sbg:x": 757.4453125,
            "sbg:y": 19.554697036743164
        },
        {
            "id": "Centrifuge_report",
            "outputSource": [
                "centrifuge_classifier_1/Centrifuge_report"
            ],
            "sbg:fileTypes": "TSV",
            "type": "File?",
            "label": "TSV report",
            "doc": "Centrifuge report.",
            "sbg:x": 61.54109191894531,
            "sbg:y": 610.2705078125
        },
        {
            "id": "out_variants",
            "outputSource": [
                "gatk_haplotypecaller_4_2_0_0/out_variants"
            ],
            "sbg:fileTypes": "VCF, VCF.GZ",
            "type": "File",
            "label": "VCF variants",
            "doc": "File to which variants should be written.",
            "secondaryFiles": [
                {
                    "pattern": ".idx",
                    "required": false
                },
                {
                    "pattern": ".tbi",
                    "required": false
                },
                {
                    "pattern": ".md5",
                    "required": false
                }
            ],
            "sbg:x": 542.1781005859375,
            "sbg:y": 215.38674926757812
        },
        {
            "id": "aligned_reads",
            "outputSource": [
                "bwa_mem_bundle_0_7_17_cwl_1_2/aligned_reads"
            ],
            "sbg:fileTypes": "SAM, BAM, CRAM",
            "type": "File?",
            "label": "BAM alignment",
            "doc": "Output SAM/BAM/CRAM file containing aligned reads.",
            "secondaryFiles": [
                {
                    "pattern": ".bai",
                    "required": false
                },
                {
                    "pattern": "^.bai",
                    "required": false
                },
                {
                    "pattern": ".crai",
                    "required": false
                },
                {
                    "pattern": "^.crai",
                    "required": false
                }
            ],
            "sbg:x": 294.1781311035156,
            "sbg:y": 430.53436279296875
        }
    ],
    "steps": [
        {
            "id": "centrifuge_classifier_1",
            "in": [
                {
                    "id": "centrifuge_index_archive",
                    "source": "centrifuge_index_archive"
                },
                {
                    "id": "input_file",
                    "source": [
                        "input_file"
                    ]
                }
            ],
            "out": [
                {
                    "id": "Centrifuge_report"
                }
            ],
            "run": {
                "cwlVersion": "sbg:draft-2",
                "class": "CommandLineTool",
                "$namespaces": {
                    "sbg": "https://sevenbridges.com"
                },
                "id": "admin/sbg-public-data/centrifuge-classifier-1/16",
                "label": "Centrifuge Classifier",
                "description": "**Centrifuge Classifier** is the main component of the Centrifuge suite, used for classification of metagenomics reads. Upon building (or downloading) and inspecting an Index (if desired), this tool can be used to process samples of interest. It can be used as a standalone tool, or as a part of **Metagenomics WGS analysis** workflow.\n\n**Centrifuge** is a novel microbial classification engine that enables rapid, accurate, and sensitive labeling of reads and quantification of species present in metagenomic samples. The system uses a novel indexing scheme based on the Burrows-Wheeler transform (BWT) and the Ferragina-Manzini (FM) index, optimized specifically for the metagenomic classification problem [1]. \n\n\nThe **Centrifuge Classifier** tool requires the following input files:\n\n- **Input reads** files with metagenomic reads - they can be in FASTQ or FASTA format, gzip-ed (extension .GZ) or bzip2-ed (extension .BZ2); in case of host-associated samples, it is presumed that files are already cleaned from  the host's genomic sequences;\n- **Reference index** in TAR format; for larger indexes the TAR.GZ format can be used, in this case we suggest the user manually set the required memory for the **Centrifuge Classifier** job by changing the value for the **Memory per job** parameter.\n\nThe results of the classification process are presented in two output files:\n\n- **Centrifuge report**, a tab delimited text file with containing results of an analysis for each taxonomic category from the reference organized into eight columns: (1) ID of a read, (2) sequence ID of the genomic sequence where the read is classified, (3) taxonomic ID of the genomic sequence from the second column, (4) the score of the classification, (5) the score for the next best classification, (6) an approximate number of base pairs of the read that match the genomic sequence, (7) the length of a read, and (8) the number of classifications for this read. The resulting line per read looks like the following:   \n     `HWUSI-EAS688_103028660:7:100:10014:18930 NC_020104.1 1269028 81 0 24 200 1`\n\n- **Classification result**, a tab delimited file with a classification summary for each genome or taxonomic unit from the reference index organized into seven columns: (1) the name of a genome, (2) taxonomic ID, (3) taxonomic rank, (4) the length of the genome sequence, (5) the number of reads classified to the provided genomic sequences, including multi-classified, (6) the number of reads uniquely classified to this particular genomic sequence, (7) the proportion of this genome normalized by its genomic length. The resulting line per genome looks like the following:   \n     `Streptococcus phage 20617 1392231 species 48800 1436 1325 0.453983`\n\n\nBased on the `--met-file` parameter value, **Centrifuge Classifier** can produce an additional output file with alignment metrics. However, it seems that this option does not work properly (see *Common Issues and Important Notes*). We suggest excluding it.\n\n*A list of **all inputs and parameters** with corresponding descriptions can be found at the end of the page.*\n\n### Common Use Cases\n\n**Centrifuge Classifier** takes files from the **Input reads** input node with raw reads (presumably cleaned from host genomic sequences as recommended by [Human Microbiome Project](https://hmpdacc.org/)) and one reference index file (in TAR ot TAR.GZ format) with microbial reference sequences. Based on the k-mer exact matching algorithm, Centrifuge assigns scores for each species on which k-mers from the read are aligned, and afterwards traverses upwards along the taxonomic tree to reduce the number of assignments, first by considering the genus that includes the largest number of species and then replacing those species with the genus. For more details on Centrifuge's algorithm see [1].\nThe results of the classification analysis are given for each reference sequence (in the **Classification result** file) and for each read (in the **Centrifuge report** file).\n\nThe index used in the analysis is Centrifuge's FM-index, whose size is reduced by compressing redundant genomic sequences. It can be downloaded from [Centrifuge's website](https://ccb.jhu.edu/software/centrifuge/manual.shtml) or created with the **Centrifuge Build** tool and/or the **Reference Index Creation workflow** provided by the Seven Bridges platform. In either case, it must be in the appropriate format (this can be checked with the **Centrifuge Inspect** tool). Based on our experience, an index should contain all of the organisms that are expected to be present in the sample. Providing an index with a smaller number of organisms (for example, in cases when the user is just interested in detecting one particular species) can result in miscalculated abundances of organisms within the sample.\n\n\n### Changes Introduced by Seven Bridges\n\n* **Centrifuge Classifier options** `--un`, `--un-conc`, `--al`, `--al-conc` and `--met-file` do not work properly, therefore all of them are excluded from the Seven Bridges version of the tool. In case a new release of the tool addresses these issues, an updated Seven Bridges version of the tool will be released as well.\n* The tool will automatically extract the index TAR file into the working directory and the basename from the **Reference genome** metadata field will be passed to **Centrifuge Classifier** using the `-x` argument.\n\n### Common Issues and Important Notes\n\n* If the index is in TAR.GZ format, the memory for **Centrifuge Classifier** should be set manually by changing the default value for the **Memory per job** parameter (in MB). Based on our experience, it would be enough to use twice as much memory as the size of the index file, with an additional 4GB overhead. For example, if the size of the index file is 8GB, the user should use 8 x 2 + 4 = 20GB, which amounts to 20 x 1024 = 20480MB.\n\n### Performance Benchmarking\n\n**Centrifuge Classifier** requires a significant amount of memory (based on our experience 4GB more than the size of the index files is suggested) in order to work properly. If the index is in TAR format, the tool will automatically allocate the required memory size. However, if the index is in TAR.GZ format, where the compression ratio is not always the same, we suggest the user manually set the required amount of memory necessary for the **Centrifuge Classifier** job. This can be done by changing the default **Memory per job** parameter value. This way, using expensive and memory overqualified instances, task failures would be avoided.\n\nIn the following table you can find estimates of **Centrifuge Classifier** running time and cost. All experiments are done with one sample. \n\n*Cost can be significantly reduced by using **spot instances**. Visit the [knowledge center](https://docs.sevenbridges.com/docs/about-spot-instances) for more details.*  \n\n| Experiment type | Input size | Duration | Cost | Instance (AWS)\n| --- | --- | --- | --- | --- |\n|p_h_v index (prokaryotic, human and viral genomes)|Index 6.9 GB, reads SRS013942 sample 2 x 1 GB| 11m| $ 0.11 | c4.2xlarge|\n|Viral index |Index 118 MB, reads SRS013942 sample 2 x 1 GB| 5m| $ 0.07 | c4.2xlarge|\n|Bacterial index |Index 13GB, reads SRS013942 sample 2 x 1 GB | 16m| $ 0.25 | c4.4xlarge|\n|Bacterial index |Index 13GB, reads SRS019027 sample 2 x 3.5 GB | 27m| $ 0.42 | c4.4xlarge|\n\n\n### References\n[1] Kim D, Song L, Breitwieser FP, and Salzberg SL. [Centrifuge: rapid and sensitive classification of metagenomic sequences](http://genome.cshlp.org/content/early/2016/11/16/gr.210641.116.abstract). Genome Research 2016",
                "baseCommand": [
                    {
                        "class": "Expression",
                        "script": "{\n  index_file = $job.inputs.centrifuge_index_archive.path\n  return \"tar -xvf \" + index_file + \"; \" + \"basename=$(ls | grep '.cf$' | head -1 | cut -d '.' -f 1); \"\n  \n}",
                        "engine": "#cwl-js-engine"
                    },
                    "centrifuge"
                ],
                "inputs": [
                    {
                        "sbg:stageInput": "link",
                        "sbg:category": "Input files",
                        "type": [
                            "File"
                        ],
                        "label": "Reference index",
                        "description": "The basename of the index for the reference genomes. The basename is the name of any of the index files up to but not including the final .1.cf / etc. centrifuge looks for the specified index first in the current directory, then in the directory specified in the CENTRIFUGE_INDEXES environment variable.",
                        "sbg:fileTypes": "TAR, TAR.GZ",
                        "id": "#centrifuge_index_archive"
                    },
                    {
                        "sbg:category": "Input files",
                        "type": [
                            {
                                "type": "array",
                                "items": "File"
                            }
                        ],
                        "inputBinding": {
                            "position": 5,
                            "separate": true,
                            "itemSeparator": " ",
                            "valueFrom": {
                                "class": "Expression",
                                "script": "{\n  filepath = ($job.inputs.input_file)[0].path\n  filename = filepath.split(\"/\").pop()\n  ext = filename.substr(filename.lastIndexOf('.') + 1)\n\n  if (ext === \"fa\" || ext === \"fasta\")\n  {\n    return \"-f \" \n  }\n  else if (ext === \"fq\" || ext === \"fastq\")\n  {\n    return \"-q \" \n  }\n}\n\n",
                                "engine": "#cwl-js-engine"
                            },
                            "sbg:cmdInclude": true
                        },
                        "label": "Input reads",
                        "description": "Read sequence in FASTQ or FASTA format. Could be also gzip'ed (extension .gz) or bzip2'ed (extension .bz2). In case of paired-end alignment it is crucial to set metadata 'paired-end' field to 1/2.",
                        "sbg:fileTypes": "FASTA, FASTQ, FA, FQ, FASTQ.GZ, FQ.GZ, FASTQ.BZ2, FQ.BZ2",
                        "id": "#input_file"
                    },
                    {
                        "sbg:stageInput": null,
                        "sbg:altPrefix": "-5",
                        "sbg:category": "Input",
                        "sbg:toolDefaultValue": "0",
                        "type": [
                            "null",
                            "int"
                        ],
                        "inputBinding": {
                            "position": 0,
                            "prefix": "--trim5",
                            "separate": true,
                            "sbg:cmdInclude": true
                        },
                        "label": "Trim from 5'",
                        "description": "Trim specific number of bases from 5' (left) end of each read before alignment (default: 0).",
                        "id": "#trim_from_5"
                    },
                    {
                        "sbg:altPrefix": "-3",
                        "sbg:category": "Input",
                        "sbg:toolDefaultValue": "0",
                        "type": [
                            "null",
                            "int"
                        ],
                        "inputBinding": {
                            "position": 0,
                            "prefix": "--trim3",
                            "separate": true,
                            "sbg:cmdInclude": true
                        },
                        "label": "Trim from 3'",
                        "description": "Trim specific number of bases from 3' (right) end of each read before alignment.",
                        "id": "#trim_from_3"
                    },
                    {
                        "sbg:altPrefix": "-u",
                        "sbg:category": "Input",
                        "sbg:toolDefaultValue": "No limit",
                        "type": [
                            "null",
                            "int"
                        ],
                        "inputBinding": {
                            "position": 0,
                            "prefix": "--upto",
                            "separate": true,
                            "sbg:cmdInclude": true
                        },
                        "label": "Align first n reads",
                        "description": "Align the first n reads or read pairs from the input (after the -s/--skip reads or pairs have been skipped), then stop.",
                        "id": "#align_first_n_reads"
                    },
                    {
                        "sbg:stageInput": null,
                        "sbg:altPrefix": "-s",
                        "sbg:category": "Input",
                        "type": [
                            "null",
                            "int"
                        ],
                        "inputBinding": {
                            "position": 0,
                            "prefix": "--skip",
                            "separate": true,
                            "sbg:cmdInclude": true
                        },
                        "label": "Skip the first n reads",
                        "description": "Skip (do not align) the first n reads or pairs in the input.",
                        "id": "#skip_n_reads"
                    },
                    {
                        "sbg:category": "Input",
                        "sbg:toolDefaultValue": "Phred+33",
                        "type": [
                            "null",
                            {
                                "type": "enum",
                                "symbols": [
                                    "Phred+33",
                                    "Phred+64",
                                    "Integers"
                                ],
                                "name": "quality"
                            }
                        ],
                        "label": "Quality scale",
                        "description": "Input qualities are ASCII chars equal to Phred+33 or Phred+64 encoding, or ASCII integers (which are treated as being on the Phred quality scale).",
                        "id": "#quality"
                    },
                    {
                        "sbg:category": "Classification",
                        "type": [
                            "null",
                            "int"
                        ],
                        "inputBinding": {
                            "position": 0,
                            "prefix": "--host-taxids",
                            "separate": true,
                            "sbg:cmdInclude": true
                        },
                        "label": "Host taxids",
                        "description": "A comma-separated list of taxonomic IDs that will be preferred in classification procedure. The descendants from these IDs will also be preferred. In case some of a read's assignments correspond to these taxonomic IDs, only those corresponding assignments will be reported.",
                        "id": "#host_taxids"
                    },
                    {
                        "sbg:category": "Classification",
                        "type": [
                            "null",
                            "int"
                        ],
                        "inputBinding": {
                            "position": 0,
                            "prefix": "--exclude-taxids",
                            "separate": true,
                            "sbg:cmdInclude": true
                        },
                        "label": "Exclude taxids",
                        "description": "A comma-separated list of taxonomic IDs that will be excluded in classification procedure. The descendants from these IDs will also be exclude.",
                        "id": "#exclude_taxids"
                    },
                    {
                        "sbg:stageInput": null,
                        "sbg:altPrefix": "-t",
                        "sbg:category": "Output options",
                        "sbg:toolDefaultValue": "Off",
                        "type": [
                            "null",
                            "boolean"
                        ],
                        "inputBinding": {
                            "position": 0,
                            "prefix": "--time",
                            "separate": true,
                            "sbg:cmdInclude": true
                        },
                        "label": "Time",
                        "description": "Print the wall-clock time required to load the index files and align the reads. This is printed to the \"standard error\" (\"stderr\") filehandle.",
                        "id": "#time"
                    },
                    {
                        "sbg:stageInput": null,
                        "sbg:category": "Output options",
                        "type": [
                            "null",
                            "boolean"
                        ],
                        "inputBinding": {
                            "position": 0,
                            "prefix": "--quiet",
                            "separate": true,
                            "sbg:cmdInclude": true
                        },
                        "label": "Quiet",
                        "description": "Print nothing besides alignments and serious errors.",
                        "id": "#quiet"
                    },
                    {
                        "sbg:stageInput": null,
                        "sbg:category": "Output options",
                        "sbg:toolDefaultValue": "Metrics disabled",
                        "type": [
                            "null",
                            "boolean"
                        ],
                        "inputBinding": {
                            "position": 0,
                            "prefix": "--met-stderr",
                            "separate": true,
                            "sbg:cmdInclude": true
                        },
                        "label": "Metrics standard error",
                        "description": "Write centrifuge metrics to the \"standard error\" (\"stderr\") filehandle. This is not mutually exclusive with --met-file. Having alignment metric can be useful for debugging certain problems, especially performance issues.",
                        "id": "#metrics_standard_error"
                    },
                    {
                        "sbg:stageInput": null,
                        "sbg:category": "Output options",
                        "sbg:toolDefaultValue": "1",
                        "type": [
                            "null",
                            "int"
                        ],
                        "inputBinding": {
                            "position": 0,
                            "prefix": "--met",
                            "separate": true,
                            "sbg:cmdInclude": true
                        },
                        "label": "Metrics",
                        "description": "Write a new centrifuge metrics record every <int> seconds. Only matters if either --met-stderr or --met-file are specified.",
                        "id": "#met"
                    },
                    {
                        "sbg:altPrefix": "--threads",
                        "sbg:category": "Performance options",
                        "sbg:toolDefaultValue": "1",
                        "type": [
                            "null",
                            "int"
                        ],
                        "inputBinding": {
                            "position": 0,
                            "prefix": "-p",
                            "separate": true,
                            "sbg:cmdInclude": true
                        },
                        "label": "Parallel threads",
                        "description": "Launch NTHREADS parallel search threads (default: 1). Threads will run on separate processors/cores and synchronize when parsing reads and outputting alignments. Searching for alignments is highly parallel, and speedup is close to linear. Increasing -p increases Centrifuge's memory footprint. E.g. when aligning to a human genome index, increasing -p from 1 to 8 increases the memory footprint by a few hundred megabytes. This option is only available if bowtie is linked with the pthreads library (i.e. if BOWTIE_PTHREADS=0 is not specified at build time).",
                        "id": "#parallel_threads"
                    },
                    {
                        "sbg:stageInput": null,
                        "sbg:category": "Performance options",
                        "type": [
                            "null",
                            "boolean"
                        ],
                        "inputBinding": {
                            "position": 0,
                            "prefix": "--mm",
                            "separate": true,
                            "sbg:cmdInclude": true
                        },
                        "label": "Memory mapping",
                        "description": "Use memory-mapped I/O to load the index, rather than typical file I/O. Memory-mapping allows many concurrent bowtie processes on the same computer to share the same memory image of the index (i.e. you pay the memory overhead just once). This facilitates memory-efficient parallelization of bowtie in situations where using -p is not possible or not preferable.",
                        "id": "#memory_mapping"
                    },
                    {
                        "sbg:stageInput": null,
                        "sbg:category": "Classification",
                        "sbg:toolDefaultValue": "22",
                        "type": [
                            "null",
                            "int"
                        ],
                        "inputBinding": {
                            "position": 0,
                            "prefix": "--min-hitlen",
                            "separate": true,
                            "sbg:cmdInclude": true
                        },
                        "label": "Minimum length of partial hits",
                        "description": "Minimum length of partial hits, which must be greater than 15.",
                        "id": "#min_hitlen"
                    },
                    {
                        "sbg:stageInput": null,
                        "sbg:category": "Classification",
                        "sbg:toolDefaultValue": "0",
                        "type": [
                            "null",
                            "int"
                        ],
                        "inputBinding": {
                            "position": 0,
                            "prefix": "--min-totallen",
                            "separate": true,
                            "sbg:cmdInclude": true
                        },
                        "label": "Minimum summed length",
                        "description": "Minimum summed length of partial hits per read.",
                        "id": "#min_totallen"
                    },
                    {
                        "sbg:category": "Input",
                        "sbg:toolDefaultValue": "False",
                        "type": [
                            "null",
                            "boolean"
                        ],
                        "inputBinding": {
                            "position": 0,
                            "prefix": "--ignore-quals",
                            "separate": true,
                            "sbg:cmdInclude": true
                        },
                        "label": "Ignore qualities",
                        "description": "Treat all quality values as 30 on Phred scale.",
                        "id": "#ignore_quals"
                    },
                    {
                        "sbg:category": "Input",
                        "sbg:toolDefaultValue": ".fq/.fastq",
                        "type": [
                            "null",
                            {
                                "type": "enum",
                                "symbols": [
                                    "FASTQ (.fq/.fastq)",
                                    "QSEQ",
                                    "FASTA (.fa/.mfa)",
                                    "raw-one-sequence-per-line",
                                    "comma-separated-lists"
                                ],
                                "name": "query_input_files"
                            }
                        ],
                        "label": "Query input files",
                        "description": "Query input files can be in FASTQ, (multi-) FASTA or QSEQ format, or one line per read.",
                        "id": "#query_input_files"
                    },
                    {
                        "sbg:category": "Input",
                        "type": [
                            "null",
                            "string"
                        ],
                        "inputBinding": {
                            "position": 0,
                            "prefix": "--sra-acc",
                            "separate": true,
                            "sbg:cmdInclude": true
                        },
                        "label": "SRA accession number",
                        "description": "Comma-separated list of SRA accession numbers, e.g. --sra-acc SRR353653,SRR353654. Information about read types is available at http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?sp=runinfo&acc=sra-acc&retmode=xml, where sra-acc is SRA accession number. If users run HISAT2 on a computer cluster, it is recommended to disable SRA-related caching (see the instruction at SRA-MANUAL).",
                        "id": "#SRA_accession_number"
                    },
                    {
                        "sbg:category": "Input",
                        "sbg:toolDefaultValue": "Off",
                        "type": [
                            "null",
                            "boolean"
                        ],
                        "inputBinding": {
                            "position": 0,
                            "prefix": "--nofw",
                            "separate": true,
                            "sbg:cmdInclude": true
                        },
                        "label": "No forward version",
                        "description": "Do not align forward (original) version of the read.",
                        "id": "#nofw"
                    },
                    {
                        "sbg:category": "Input",
                        "sbg:toolDefaultValue": "Off",
                        "type": [
                            "null",
                            "boolean"
                        ],
                        "inputBinding": {
                            "position": 0,
                            "prefix": "--norc",
                            "separate": true,
                            "sbg:cmdInclude": true
                        },
                        "label": "No reverse complement",
                        "description": "Do not align reverse-complement version of the read.",
                        "id": "#norc"
                    },
                    {
                        "sbg:category": "Output",
                        "sbg:toolDefaultValue": "tab",
                        "type": [
                            "null",
                            "string"
                        ],
                        "inputBinding": {
                            "position": 0,
                            "prefix": "--out-fmt",
                            "separate": true,
                            "sbg:cmdInclude": true
                        },
                        "label": "Output format",
                        "description": "Define output format, either 'tab' or 'sam'.",
                        "id": "#out_fmt"
                    },
                    {
                        "sbg:category": "Output",
                        "sbg:toolDefaultValue": "readID,seqID,taxID,score,2ndBestScore,hitLength,queryLength,numMatches",
                        "type": [
                            "null",
                            "string"
                        ],
                        "inputBinding": {
                            "position": 0,
                            "prefix": "--tab-fmt-cols",
                            "separate": true,
                            "sbg:cmdInclude": true
                        },
                        "label": "Columns in tabular format",
                        "description": "Columns in tabular format, comma separated.",
                        "id": "#tab_fmt_cols"
                    },
                    {
                        "sbg:category": "Other",
                        "sbg:toolDefaultValue": "Off",
                        "type": [
                            "null",
                            "boolean"
                        ],
                        "inputBinding": {
                            "position": 0,
                            "prefix": "--qc-filter",
                            "separate": true,
                            "sbg:cmdInclude": true
                        },
                        "label": "QC filter",
                        "description": "Filters out reads that are bad according to QSEQ filter.",
                        "id": "#qc_filter"
                    },
                    {
                        "sbg:category": "Execution",
                        "sbg:toolDefaultValue": "4096",
                        "type": [
                            "null",
                            "int"
                        ],
                        "label": "Memory per job",
                        "description": "Amount of RAM memory (in MB) to be used per job.",
                        "id": "#memory_per_job"
                    },
                    {
                        "sbg:category": "Output",
                        "sbg:toolDefaultValue": "TXT",
                        "type": [
                            "null",
                            "boolean"
                        ],
                        "label": "Paired reads not aligned concordantly",
                        "description": "Writes pairs that didn't align concordantly to the output file",
                        "id": "#un_conc"
                    },
                    {
                        "sbg:category": "Output",
                        "sbg:toolDefaultValue": "TXT",
                        "type": [
                            "null",
                            "boolean"
                        ],
                        "label": "Paired reads aligned concordantly",
                        "description": "Writes pairs that aligned concordantly at least once to the output file.",
                        "id": "#al_conc"
                    },
                    {
                        "sbg:stageInput": null,
                        "sbg:category": "Output",
                        "sbg:toolDefaultValue": "TXT",
                        "type": [
                            "null",
                            "boolean"
                        ],
                        "label": "Unpaired reads that aligned at least once",
                        "description": "Writes unpaired reads that aligned at least once to the output file.",
                        "id": "#aligned_unpaired_reads"
                    },
                    {
                        "sbg:stageInput": null,
                        "sbg:category": "Output",
                        "sbg:toolDefaultValue": "T",
                        "type": [
                            "null",
                            "boolean"
                        ],
                        "label": "Unpaired reads that didn't align",
                        "description": "Writes unpaired reads that didn't align to the reference in the output file.",
                        "id": "#unaligned_unpaired_reads"
                    }
                ],
                "outputs": [
                    {
                        "type": [
                            "null",
                            "File"
                        ],
                        "label": "Centrifuge report",
                        "description": "Centrifuge report.",
                        "sbg:fileTypes": "TSV",
                        "outputBinding": {
                            "glob": "*Centrifuge_report.tsv",
                            "sbg:inheritMetadataFrom": "#input_file"
                        },
                        "id": "#Centrifuge_report"
                    }
                ],
                "requirements": [
                    {
                        "class": "ExpressionEngineRequirement",
                        "id": "#cwl-js-engine",
                        "requirements": [
                            {
                                "class": "DockerRequirement",
                                "dockerPull": "rabix/js-engine"
                            }
                        ]
                    }
                ],
                "hints": [
                    {
                        "class": "sbg:CPURequirement",
                        "value": {
                            "class": "Expression",
                            "script": "{\n  if($job.inputs.parallel_threads){\n  \treturn $job.inputs.parallel_threads\n  }\n  return 1\n}",
                            "engine": "#cwl-js-engine"
                        }
                    },
                    {
                        "class": "sbg:MemRequirement",
                        "value": {
                            "class": "Expression",
                            "script": "{  \n  if($job.inputs.memory_per_job){\n  \t\treturn $job.inputs.memory_per_job\n\t}\n\telse {\n      \tinput_size = $job.inputs.centrifuge_index_archive.size\n        input_MB = input_size/1048576\n       \n        \n        if ($job.inputs.centrifuge_index_archive.path.split('.').pop()==\"gz\"){\n          return 2*Math.ceil(input_MB) + 4096\n          \n        }\n      \telse {\n        \treturn Math.ceil(input_MB) + 4096\n          \n        }\n      \n\t}\n\n}\n\n\n",
                            "engine": "#cwl-js-engine"
                        }
                    },
                    {
                        "class": "DockerRequirement",
                        "dockerPull": "images.sbgenomics.com/aleksandar_danicic/centrifuge:1.0.3_feb2018"
                    }
                ],
                "arguments": [
                    {
                        "position": 1,
                        "prefix": "-x",
                        "separate": true,
                        "valueFrom": {
                            "class": "Expression",
                            "script": "{\n  return \"$basename\"\n}",
                            "engine": "#cwl-js-engine"
                        }
                    },
                    {
                        "position": 5,
                        "separate": false,
                        "valueFrom": {
                            "class": "Expression",
                            "script": "{\n  cmd = \"\"\n  reads = [].concat($job.inputs.input_file)\n  reads1 = [];\n  reads2 = [];\n  u_reads = [];\n  for (var i = 0; i < reads.length; i++){\n      if (reads[i].metadata.paired_end == 1){\n        reads1.push(reads[i].path);\n      }\n      else if (reads[i].metadata.paired_end == 2){\n        reads2.push(reads[i].path);\n      }\n    else {\n      u_reads.push(reads[i].path);\n     }\n    }\n  if (reads1.length > 0 & reads1.length == reads2.length){\n      cmd = \"-1 \" + reads1.join(\",\") + \" -2 \" + reads2.join(\",\");\n  }\n  if (u_reads.length > 0){\n      cmd = \" -U \" + u_reads.join(\",\");\n  }\n  return cmd\n}",
                            "engine": "#cwl-js-engine"
                        }
                    },
                    {
                        "position": 6,
                        "prefix": "-S",
                        "separate": true,
                        "valueFrom": {
                            "class": "Expression",
                            "script": "{\n  filepath = ($job.inputs.input_file)[0].path\n  filename = filepath.split(\"/\").pop()  \n  basename = filename.substr(0,filename.lastIndexOf(\".\")) \n  new_filename = basename + \".Classification_result.txt\" \n  \n  return new_filename;\n}\n\n\n\n\n",
                            "engine": "#cwl-js-engine"
                        }
                    },
                    {
                        "position": 0,
                        "prefix": "--report-file",
                        "separate": true,
                        "valueFrom": {
                            "class": "Expression",
                            "script": "{\n  filepath = ($job.inputs.input_file)[0].path\n  filename = filepath.split(\"/\").pop()  \n  basename = filename.substr(0,filename.lastIndexOf(\".\")) \n  new_filename = basename + \".Centrifuge_report.tsv\" \n  \n  return new_filename;\n}\n\n\n\n",
                            "engine": "#cwl-js-engine"
                        }
                    },
                    {
                        "position": 0,
                        "separate": true,
                        "valueFrom": {
                            "class": "Expression",
                            "script": "{\n\n  if ($job.inputs.query_input_files == \"FASTQ\")\n  {\n    return \"-q\" \n  }\n  else if ($job.inputs.query_input_files == \"QSEQ\")\n  {\n    return \"-qseq\" \n  }\n  else if ($job.inputs.query_input_files == \"FASTA\")\n  {\n    return \"-f\" \n  }\n  else if ($job.inputs.query_input_files == \"raw-one-sequence-per-line\")\n  {\n    return \"-r\" \n  }\n  else if ($job.inputs.query_input_files == \"comma-separated-lists\")\n  {\n    return \"-c \" \n  }\n  \n}",
                            "engine": "#cwl-js-engine"
                        }
                    },
                    {
                        "position": 0,
                        "separate": true,
                        "valueFrom": {
                            "class": "Expression",
                            "script": "{\n\n  if ($job.inputs.quality == \"Phred+33\")\n  {\n    return \"--phred33\" \n  }\n  else if ($job.inputs.quality == \"Phred+64\")\n  {\n    return \"--phred64\" \n  }\n  else if ($job.inputs.quality == \"Integers\")\n  {\n    return \"--int-quals\" \n  }\n    \n}",
                            "engine": "#cwl-js-engine"
                        }
                    },
                    {
                        "position": 10,
                        "separate": true,
                        "valueFrom": {
                            "class": "Expression",
                            "script": "{\n  return \"; echo \\\"\" + $job.inputs.centrifuge_index_archive.path.split('.').pop() + \"\\\"\"\n  \n  \n}",
                            "engine": "#cwl-js-engine"
                        }
                    },
                    {
                        "position": 0,
                        "prefix": "--un-conc",
                        "separate": true,
                        "valueFrom": {
                            "class": "Expression",
                            "script": "{\n  filepath = ($job.inputs.input_file)[0].path\n  filename = filepath.split(\"/\").pop()  \n  basename = filename.substr(0,filename.lastIndexOf(\".\")) \n  new_filename = basename + \".Centrifuge_un_conc.txt\" \n  \n  if ($job.inputs.un_conc == true)\n  {\n    return new_filename; \n  }\n}\n\n\n\n\n",
                            "engine": "#cwl-js-engine"
                        }
                    },
                    {
                        "position": 0,
                        "prefix": "--al-conc",
                        "separate": true,
                        "valueFrom": {
                            "class": "Expression",
                            "script": "{\n  filepath = ($job.inputs.input_file)[0].path\n  filename = filepath.split(\"/\").pop()  \n  basename = filename.substr(0,filename.lastIndexOf(\".\")) \n  new_filename = basename + \".Centrifuge_al_conc.txt\" \n  \n  if ($job.inputs.al_conc == true)\n  {\n    return new_filename; \n  }\n}\n\n\n\n\n\n",
                            "engine": "#cwl-js-engine"
                        }
                    },
                    {
                        "position": 0,
                        "prefix": "--al",
                        "separate": true,
                        "valueFrom": {
                            "class": "Expression",
                            "script": "{\n  filepath = ($job.inputs.input_file)[0].path\n  filename = filepath.split(\"/\").pop()  \n  basename = filename.substr(0,filename.lastIndexOf(\".\")) \n  new_filename = basename + \".Centrifuge_aligned_unpaired_reads.txt\" \n  \n  if ($job.inputs.aligned_unpaired_reads == true)\n  {\n    return new_filename; \n  }\n}\n\n",
                            "engine": "#cwl-js-engine"
                        }
                    },
                    {
                        "position": 0,
                        "prefix": "--un",
                        "separate": true,
                        "valueFrom": {
                            "class": "Expression",
                            "script": "{\n  filepath = ($job.inputs.input_file)[0].path\n  filename = filepath.split(\"/\").pop()  \n  basename = filename.substr(0,filename.lastIndexOf(\".\")) \n  new_filename = basename + \".Centrifuge_unaligned_unpaired_reads.txt\" \n  \n  if ($job.inputs.unaligned_unpaired_reads == true)\n  {\n    return new_filename; \n  }\n}\n\n\n",
                            "engine": "#cwl-js-engine"
                        }
                    }
                ],
                "stdout": "test.log",
                "sbg:categories": [
                    "Metagenomics",
                    "Taxonomic Profiling"
                ],
                "sbg:image_url": null,
                "sbg:cmdPreview": "tar -xvf /path/to/centrifuge_index_archive.ext.gz; basename=$(ls | grep '.cf$' | head -1 | cut -d '.' -f 1);  centrifuge --report-file input_file-1.Centrifuge_report.tsv    --phred33 --un-conc input_file-1.Centrifuge_un_conc.txt --al-conc input_file-1.Centrifuge_al_conc.txt --al input_file-1.Centrifuge_aligned_unpaired_reads.txt --un input_file-1.Centrifuge_unaligned_unpaired_reads.txt -x $basename -1 /path/to/input_file-1.ext -2 /path/to/input_file-2.ext -S input_file-1.Classification_result.txt  ; echo \"gz\" > test.log",
                "sbg:revisionsInfo": [
                    {
                        "sbg:revision": 0,
                        "sbg:modifiedBy": "admin",
                        "sbg:modifiedOn": 1509721131,
                        "sbg:revisionNotes": "Copy of aleksandar_danicic/centrifuge-dev/centrifuge-classifier-1/31"
                    },
                    {
                        "sbg:revision": 1,
                        "sbg:modifiedBy": "admin",
                        "sbg:modifiedOn": 1509721131,
                        "sbg:revisionNotes": "Copy of aleksandar_danicic/centrifuge-dev/centrifuge-classifier-1/34"
                    },
                    {
                        "sbg:revision": 2,
                        "sbg:modifiedBy": "admin",
                        "sbg:modifiedOn": 1509721131,
                        "sbg:revisionNotes": "Copy of aleksandar_danicic/centrifuge-dev/centrifuge-classifier-1/35"
                    },
                    {
                        "sbg:revision": 3,
                        "sbg:modifiedBy": "admin",
                        "sbg:modifiedOn": 1509721131,
                        "sbg:revisionNotes": "Copy of aleksandar_danicic/centrifuge-dev/centrifuge-classifier-1/36"
                    },
                    {
                        "sbg:revision": 4,
                        "sbg:modifiedBy": "admin",
                        "sbg:modifiedOn": 1509721131,
                        "sbg:revisionNotes": "memory per job script changed"
                    },
                    {
                        "sbg:revision": 5,
                        "sbg:modifiedBy": "admin",
                        "sbg:modifiedOn": 1509721131,
                        "sbg:revisionNotes": "Copy of aleksandar_danicic/centrifuge-dev/centrifuge-classifier-1/37"
                    },
                    {
                        "sbg:revision": 6,
                        "sbg:modifiedBy": "admin",
                        "sbg:modifiedOn": 1509721131,
                        "sbg:revisionNotes": "Copy of aleksandar_danicic/centrifuge-dev/centrifuge-classifier-1/44"
                    },
                    {
                        "sbg:revision": 7,
                        "sbg:modifiedBy": "admin",
                        "sbg:modifiedOn": 1509721131,
                        "sbg:revisionNotes": "Copy of aleksandar_danicic/centrifuge-dev/centrifuge-classifier-1/45"
                    },
                    {
                        "sbg:revision": 8,
                        "sbg:modifiedBy": "admin",
                        "sbg:modifiedOn": 1509721131,
                        "sbg:revisionNotes": "Copy of aleksandar_danicic/centrifuge-dev/centrifuge-classifier-1/46"
                    },
                    {
                        "sbg:revision": 9,
                        "sbg:modifiedBy": "admin",
                        "sbg:modifiedOn": 1509721131,
                        "sbg:revisionNotes": "Copy of aleksandar_danicic/centrifuge-dev/centrifuge-classifier-1/47"
                    },
                    {
                        "sbg:revision": 10,
                        "sbg:modifiedBy": "admin",
                        "sbg:modifiedOn": 1509721131,
                        "sbg:revisionNotes": "Copy of aleksandar_danicic/centrifuge-dev/centrifuge-classifier-1/53"
                    },
                    {
                        "sbg:revision": 11,
                        "sbg:modifiedBy": "admin",
                        "sbg:modifiedOn": 1510650872,
                        "sbg:revisionNotes": "Copy of aleksandar_danicic/centrifuge-dev/centrifuge-classifier-1/69"
                    },
                    {
                        "sbg:revision": 12,
                        "sbg:modifiedBy": "admin",
                        "sbg:modifiedOn": 1513010011,
                        "sbg:revisionNotes": "Copy of aleksandar_danicic/centrifuge-dev/centrifuge-classifier-1/70"
                    },
                    {
                        "sbg:revision": 13,
                        "sbg:modifiedBy": "admin",
                        "sbg:modifiedOn": 1513010012,
                        "sbg:revisionNotes": "Copy of aleksandar_danicic/centrifuge-dev/centrifuge-classifier-1/73"
                    },
                    {
                        "sbg:revision": 14,
                        "sbg:modifiedBy": "admin",
                        "sbg:modifiedOn": 1527589668,
                        "sbg:revisionNotes": "Copy of aleksandar_danicic/centrifuge-dev/centrifuge-classifier-1/83"
                    },
                    {
                        "sbg:revision": 15,
                        "sbg:modifiedBy": "admin",
                        "sbg:modifiedOn": 1530200317,
                        "sbg:revisionNotes": "Copy of aleksandar_danicic/centrifuge-dev/centrifuge-classifier-1/85"
                    },
                    {
                        "sbg:revision": 16,
                        "sbg:modifiedBy": "admin",
                        "sbg:modifiedOn": 1648560739,
                        "sbg:revisionNotes": "Copy of Rev. 86 from the Dev project; category added"
                    }
                ],
                "sbg:job": {
                    "inputs": {
                        "centrifuge_index_archive": {
                            "class": "File",
                            "path": "/path/to/centrifuge_index_archive.ext.gz",
                            "size": 7400140800,
                            "metadata": {
                                "reference_genome": "index"
                            },
                            "secondaryFiles": []
                        },
                        "nofw": false,
                        "unaligned_unpaired_reads": true,
                        "tab_fmt_cols": "tab_fmt_cols-string-value",
                        "host_taxids": null,
                        "metrics_standard_error": false,
                        "qc_filter": false,
                        "ignore_quals": false,
                        "quiet": false,
                        "trim_from_3": null,
                        "input_file": [
                            {
                                "class": "File",
                                "path": "/path/to/input_file-1.ext",
                                "size": 0,
                                "metadata": {
                                    "paired_end": "1"
                                },
                                "secondaryFiles": []
                            },
                            {
                                "class": "File",
                                "path": "/path/to/input_file-2.ext",
                                "size": 0,
                                "metadata": {
                                    "paired_end": "2"
                                },
                                "secondaryFiles": []
                            }
                        ],
                        "SRA_accession_number": "SRA_accession_number-string-value",
                        "al_conc": true,
                        "memory_mapping": false,
                        "aligned_unpaired_reads": true,
                        "min_totallen": null,
                        "quality": "Phred+33",
                        "align_first_n_reads": null,
                        "time": false,
                        "exclude_taxids": null,
                        "norc": false,
                        "skip_n_reads": null,
                        "query_input_files": "FASTQ (.fq/.fastq)",
                        "out_fmt": "out_fmt-string-value",
                        "met": null,
                        "memory_per_job": null,
                        "parallel_threads": null,
                        "trim_from_5": null,
                        "un_conc": true,
                        "min_hitlen": null
                    },
                    "allocatedResources": {
                        "mem": 18212,
                        "cpu": 1
                    }
                },
                "sbg:toolkitVersion": "1.0.3",
                "sbg:projectName": "SBG Public data",
                "sbg:links": [
                    {
                        "label": "Homepage",
                        "id": "https://ccb.jhu.edu/software/centrifuge/manual.shtml"
                    },
                    {
                        "label": "Source code",
                        "id": "https://github.com/infphilo/centrifuge"
                    }
                ],
                "sbg:toolkit": "centrifuge",
                "sbg:toolAuthor": "John Hopkins University, Center for Computational Biology",
                "sbg:expand_workflow": false,
                "sbg:appVersion": [
                    "sbg:draft-2"
                ],
                "sbg:id": "admin/sbg-public-data/centrifuge-classifier-1/16",
                "sbg:revision": 16,
                "sbg:revisionNotes": "Copy of Rev. 86 from the Dev project; category added",
                "sbg:modifiedOn": 1648560739,
                "sbg:modifiedBy": "admin",
                "sbg:createdOn": 1509721131,
                "sbg:createdBy": "admin",
                "sbg:project": "admin/sbg-public-data",
                "sbg:sbgMaintained": false,
                "sbg:validationErrors": [],
                "sbg:contributors": [
                    "admin"
                ],
                "sbg:latestRevision": 16,
                "sbg:publisher": "sbg",
                "sbg:content_hash": "ab156bd72b59727a44e87a22d020e16c3adc8ad737edf2d29bae6c39de17c8672",
                "sbg:workflowLanguage": "CWL"
            },
            "label": "Centrifuge Classifier",
            "sbg:x": -200.27734375,
            "sbg:y": 611.1882934570312
        },
        {
            "id": "sbg_fasta_indices",
            "in": [
                {
                    "id": "reference",
                    "source": "reference"
                }
            ],
            "out": [
                {
                    "id": "fasta_reference"
                }
            ],
            "run": {
                "cwlVersion": "sbg:draft-2",
                "class": "CommandLineTool",
                "$namespaces": {
                    "sbg": "https://sevenbridges.com"
                },
                "id": "admin/sbg-public-data/sbg-fasta-indices/22",
                "label": "SBG FASTA Indices",
                "description": "Create indices for FASTA file.\n\n###**Overview**  \n\nTool allows creating FASTA dictionary and index simultaneously which is necessary for running GATK tools. This version of tool for indexing uses SAMtools faidx command (toolkit version 1.9), while for the FASTA dictionary is used CreateFastaDictionary (GATK toolkit version 4.1.0.0).\n\n\n###**Inputs**  \n\n- FASTA file \n\n###**Output**  \n\n- FASTA Reference file\n- FASTA Index file\n- FASTA Dictionary file\n\n\n###**Changes made by Seven Bridges**\n\nCreateFastaDictionary function creates a DICT file describing the contents of the FASTA file. Parameter -UR was added to the command line that sets the UR field to just the Reference file name, instead of the whole path to file. This allows Memoisation feature of the platform to work.",
                "baseCommand": [
                    "samtools",
                    "faidx"
                ],
                "inputs": [
                    {
                        "sbg:stageInput": "link",
                        "sbg:category": "Input files",
                        "type": [
                            "File"
                        ],
                        "inputBinding": {
                            "position": 0,
                            "separate": true,
                            "sbg:cmdInclude": true
                        },
                        "label": "FASTA file",
                        "description": "FASTA file to be indexed",
                        "sbg:fileTypes": "FASTA, FA, FA.GZ, FASTA.GZ",
                        "id": "#reference"
                    },
                    {
                        "sbg:category": "Execution",
                        "sbg:toolDefaultValue": "2048",
                        "type": [
                            "null",
                            "int"
                        ],
                        "label": "Memory per job",
                        "description": "Memory in megabytes required for each execution of the tool.",
                        "id": "#memory_per_job"
                    }
                ],
                "outputs": [
                    {
                        "type": [
                            "null",
                            "File"
                        ],
                        "label": "Reference",
                        "sbg:fileTypes": "FASTA",
                        "outputBinding": {
                            "glob": {
                                "script": "{\n  return $job.inputs.reference.path.split('/').pop()\n}",
                                "class": "Expression",
                                "engine": "#cwl-js-engine"
                            },
                            "sbg:inheritMetadataFrom": "#reference",
                            "secondaryFiles": [
                                ".fai",
                                "^.dict",
                                "^^.dict"
                            ]
                        },
                        "id": "#fasta_reference"
                    }
                ],
                "requirements": [
                    {
                        "class": "ExpressionEngineRequirement",
                        "id": "#cwl-js-engine",
                        "requirements": [
                            {
                                "dockerPull": "rabix/js-engine",
                                "class": "DockerRequirement"
                            }
                        ]
                    }
                ],
                "hints": [
                    {
                        "class": "sbg:CPURequirement",
                        "value": 1
                    },
                    {
                        "class": "sbg:MemRequirement",
                        "value": {
                            "script": "{\n  if($job.inputs.memory_per_job)return $job.inputs.memory_per_job + 500\n  else return 2548\n}",
                            "class": "Expression",
                            "engine": "#cwl-js-engine"
                        }
                    },
                    {
                        "class": "DockerRequirement",
                        "dockerImageId": "b177f5bd06db",
                        "dockerPull": "images.sbgenomics.com/vladimirk/gatk4-samtools:4.1.4.0-1.9"
                    }
                ],
                "arguments": [
                    {
                        "position": 1,
                        "prefix": "&&",
                        "separate": true,
                        "valueFrom": {
                            "class": "Expression",
                            "engine": "#cwl-js-engine",
                            "script": "{\n  memory = '2048'\n  if ($job.inputs.memory_per_job){\n    memory = $job.inputs.memory_per_job\n  }\n  filename = $job.inputs.reference.path.split('/').pop()\n  basename = filename.split('.')\n  if (filename.endsWith('.gz')){\n    basename.pop()\n  }\n  basename.pop()\n  name = basename.join('.')\n  return 'java -Xmx' + memory + 'M -jar /gatk/gatk-package-4.1.0.0-local.jar CreateSequenceDictionary -R=' + $job.inputs.reference.path + ' -O=' + name + '.dict'\n}"
                        }
                    },
                    {
                        "position": 3,
                        "prefix": "-UR=",
                        "separate": false,
                        "valueFrom": {
                            "class": "Expression",
                            "engine": "#cwl-js-engine",
                            "script": "{\n  return $job.inputs.reference.path.split('/')[ $job.inputs.reference.path.split('/').length - 1]\n}"
                        }
                    }
                ],
                "sbg:toolAuthor": "Sanja Mijalkovic, Seven Bridges Genomics, <sanja.mijalkovic@sbgenomics.com>",
                "sbg:cmdPreview": "samtools faidx  /path/to/reference.fa.gz && java -Xmx10M -jar /gatk/gatk-package-4.1.0.0-local.jar CreateSequenceDictionary -R=/path/to/reference.fa.gz -O=reference.dict -UR=reference.fa.gz",
                "sbg:toolkit": "SBGTools",
                "sbg:image_url": null,
                "sbg:job": {
                    "allocatedResources": {
                        "cpu": 1,
                        "mem": 510
                    },
                    "inputs": {
                        "reference": {
                            "secondaryFiles": [],
                            "path": "/path/to/reference.fa.gz",
                            "size": 0,
                            "class": "File"
                        },
                        "memory_per_job": 10
                    }
                },
                "sbg:projectName": "SBG Public data",
                "sbg:categories": [
                    "Indexing"
                ],
                "sbg:revisionsInfo": [
                    {
                        "sbg:revision": 0,
                        "sbg:modifiedBy": "sanja.mijalkovic",
                        "sbg:modifiedOn": 1448043983,
                        "sbg:revisionNotes": null
                    },
                    {
                        "sbg:revision": 1,
                        "sbg:modifiedBy": "djordje_klisic",
                        "sbg:modifiedOn": 1459163478,
                        "sbg:revisionNotes": null
                    },
                    {
                        "sbg:revision": 2,
                        "sbg:modifiedBy": "djordje_klisic",
                        "sbg:modifiedOn": 1459163478,
                        "sbg:revisionNotes": null
                    },
                    {
                        "sbg:revision": 3,
                        "sbg:modifiedBy": "djordje_klisic",
                        "sbg:modifiedOn": 1459163478,
                        "sbg:revisionNotes": null
                    },
                    {
                        "sbg:revision": 4,
                        "sbg:modifiedBy": "djordje_klisic",
                        "sbg:modifiedOn": 1459163478,
                        "sbg:revisionNotes": null
                    },
                    {
                        "sbg:revision": 5,
                        "sbg:modifiedBy": "admin",
                        "sbg:modifiedOn": 1504629640,
                        "sbg:revisionNotes": "Removed python script. Changed docker to just samtools and picard. Wrapped both faidx and CreateSequenceDictionary and exposed memory parameter for java execution."
                    },
                    {
                        "sbg:revision": 6,
                        "sbg:modifiedBy": "admin",
                        "sbg:modifiedOn": 1506681176,
                        "sbg:revisionNotes": "Changed join to join('.')."
                    },
                    {
                        "sbg:revision": 7,
                        "sbg:modifiedBy": "admin",
                        "sbg:modifiedOn": 1521479429,
                        "sbg:revisionNotes": "Added support for FA.GZ, FASTA.GZ"
                    },
                    {
                        "sbg:revision": 8,
                        "sbg:modifiedBy": "admin",
                        "sbg:modifiedOn": 1521479429,
                        "sbg:revisionNotes": "Added secondary .dict support for fasta.gz"
                    },
                    {
                        "sbg:revision": 9,
                        "sbg:modifiedBy": "admin",
                        "sbg:modifiedOn": 1530631714,
                        "sbg:revisionNotes": "returned to rev 8"
                    },
                    {
                        "sbg:revision": 10,
                        "sbg:modifiedBy": "admin",
                        "sbg:modifiedOn": 1530631714,
                        "sbg:revisionNotes": "rev 7"
                    },
                    {
                        "sbg:revision": 11,
                        "sbg:modifiedBy": "admin",
                        "sbg:modifiedOn": 1530631714,
                        "sbg:revisionNotes": "rev 9: Added secondary .dict support"
                    },
                    {
                        "sbg:revision": 12,
                        "sbg:modifiedBy": "admin",
                        "sbg:modifiedOn": 1545924498,
                        "sbg:revisionNotes": "Updated version for samtools (1.9) and picard (2.18.14)"
                    },
                    {
                        "sbg:revision": 13,
                        "sbg:modifiedBy": "admin",
                        "sbg:modifiedOn": 1545924498,
                        "sbg:revisionNotes": "Reverted."
                    },
                    {
                        "sbg:revision": 14,
                        "sbg:modifiedBy": "admin",
                        "sbg:modifiedOn": 1575892952,
                        "sbg:revisionNotes": "Added URI to eliminate randomness in .dict"
                    },
                    {
                        "sbg:revision": 15,
                        "sbg:modifiedBy": "admin",
                        "sbg:modifiedOn": 1575892952,
                        "sbg:revisionNotes": "Added URI to remove randomness"
                    },
                    {
                        "sbg:revision": 16,
                        "sbg:modifiedBy": "admin",
                        "sbg:modifiedOn": 1575892952,
                        "sbg:revisionNotes": "Updated to GATK 4.1.0.0"
                    },
                    {
                        "sbg:revision": 17,
                        "sbg:modifiedBy": "admin",
                        "sbg:modifiedOn": 1575892952,
                        "sbg:revisionNotes": "bug fix"
                    },
                    {
                        "sbg:revision": 18,
                        "sbg:modifiedBy": "admin",
                        "sbg:modifiedOn": 1575892952,
                        "sbg:revisionNotes": "description"
                    },
                    {
                        "sbg:revision": 19,
                        "sbg:modifiedBy": "admin",
                        "sbg:modifiedOn": 1575892952,
                        "sbg:revisionNotes": "updated command line preview"
                    },
                    {
                        "sbg:revision": 20,
                        "sbg:modifiedBy": "admin",
                        "sbg:modifiedOn": 1630934170,
                        "sbg:revisionNotes": "Tool description update to clarify it only takes one FASTA file"
                    },
                    {
                        "sbg:revision": 21,
                        "sbg:modifiedBy": "admin",
                        "sbg:modifiedOn": 1648035724,
                        "sbg:revisionNotes": "Updated categories"
                    },
                    {
                        "sbg:revision": 22,
                        "sbg:modifiedBy": "admin",
                        "sbg:modifiedOn": 1682685128,
                        "sbg:revisionNotes": "Categories updated"
                    }
                ],
                "sbg:license": "Apache License 2.0",
                "sbg:expand_workflow": false,
                "sbg:appVersion": [
                    "sbg:draft-2"
                ],
                "sbg:id": "admin/sbg-public-data/sbg-fasta-indices/22",
                "sbg:revision": 22,
                "sbg:revisionNotes": "Categories updated",
                "sbg:modifiedOn": 1682685128,
                "sbg:modifiedBy": "admin",
                "sbg:createdOn": 1448043983,
                "sbg:createdBy": "sanja.mijalkovic",
                "sbg:project": "admin/sbg-public-data",
                "sbg:sbgMaintained": false,
                "sbg:validationErrors": [],
                "sbg:contributors": [
                    "admin",
                    "djordje_klisic",
                    "sanja.mijalkovic"
                ],
                "sbg:latestRevision": 22,
                "sbg:publisher": "sbg",
                "sbg:content_hash": "a21b3ebe8e82b8fd5100676f5de37ea9b35992d0cbbb0c97a62c7e6a8dea4d620",
                "sbg:workflowLanguage": "CWL"
            },
            "label": "SBG FASTA Indices",
            "sbg:x": -200.25924682617188,
            "sbg:y": 87.77777862548828
        },
        {
            "id": "gatk_haplotypecaller_4_2_0_0",
            "in": [
                {
                    "id": "in_alignments",
                    "source": [
                        "bwa_mem_bundle_0_7_17_cwl_1_2/aligned_reads"
                    ]
                },
                {
                    "id": "in_reference",
                    "source": "sbg_fasta_indices/fasta_reference"
                },
                {
                    "id": "output_extension",
                    "default": "vcf"
                }
            ],
            "out": [
                {
                    "id": "out_variants"
                }
            ],
            "run": {
                "class": "CommandLineTool",
                "cwlVersion": "v1.2",
                "$namespaces": {
                    "sbg": "https://sevenbridges.com"
                },
                "id": "admin/sbg-public-data/gatk-haplotypecaller-4-2-0-0/6",
                "baseCommand": [
                    "/opt/gatk-4.2.0.0/gatk",
                    "--java-options"
                ],
                "inputs": [
                    {
                        "sbg:category": "Advanced Arguments",
                        "sbg:toolDefaultValue": "0.002",
                        "id": "active_probability_threshold",
                        "type": "float?",
                        "inputBinding": {
                            "prefix": "--active-probability-threshold",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Active probability threshold",
                        "doc": "Minimum probability for a locus to be considered active."
                    },
                    {
                        "sbg:category": "Advanced Arguments",
                        "sbg:toolDefaultValue": "false",
                        "id": "adaptive_pruning",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--adaptive-pruning",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Adaptive pruning",
                        "doc": "Use Mutect2's adaptive graph pruning algorithm."
                    },
                    {
                        "sbg:category": "Advanced Arguments",
                        "sbg:toolDefaultValue": "0.001",
                        "id": "adaptive_pruning_initial_error_rate",
                        "type": "float?",
                        "inputBinding": {
                            "prefix": "--adaptive-pruning-initial-error-rate",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Adaptive pruning initial error rate",
                        "doc": "Initial base error rate estimate for adaptive pruning."
                    },
                    {
                        "sbg:altPrefix": "-add-output-sam-program-record",
                        "sbg:category": "Optional Arguments",
                        "sbg:toolDefaultValue": "true",
                        "id": "add_output_sam_program_record",
                        "type": [
                            "null",
                            {
                                "type": "enum",
                                "symbols": [
                                    "true",
                                    "false"
                                ],
                                "name": "add_output_sam_program_record"
                            }
                        ],
                        "inputBinding": {
                            "prefix": "--add-output-sam-program-record",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Add output SAM program record",
                        "doc": "If true, adds a PG tag to created SAM/BAM/CRAM files."
                    },
                    {
                        "sbg:altPrefix": "-add-output-vcf-command-line",
                        "sbg:category": "Optional Arguments",
                        "sbg:toolDefaultValue": "true",
                        "id": "add_output_vcf_command_line",
                        "type": [
                            "null",
                            {
                                "type": "enum",
                                "symbols": [
                                    "true",
                                    "false"
                                ],
                                "name": "add_output_vcf_command_line"
                            }
                        ],
                        "inputBinding": {
                            "prefix": "--add-output-vcf-command-line",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Add output VCF command line",
                        "doc": "If true, adds a command line header line to created VCF files."
                    },
                    {
                        "sbg:category": "Advanced Arguments",
                        "sbg:toolDefaultValue": "false",
                        "id": "all_site_pls",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--all-site-pls",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Annotate all sites with PLs",
                        "doc": "Advanced, experimental argument: if SNP likelihood model is specified, and if EMIT_ALL_SITES output mode is set, when we set this argument then we will also emit PLs at all sites. This will give a measure of reference confidence and a measure of which alt alleles are more plausible (if any). WARNINGS: - This feature will inflate VCF file size considerably. - All SNP ALT alleles will be emitted with corresponding 10 PL values. - An error will be emitted if EMIT_ALL_SITES is not set, or if anything other than diploid SNP model is used"
                    },
                    {
                        "sbg:category": "Optional Arguments",
                        "sbg:toolDefaultValue": "null",
                        "id": "alleles",
                        "type": "File?",
                        "inputBinding": {
                            "prefix": "--alleles",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Alleles",
                        "doc": "The set of alleles to force-call regardless of evidence.",
                        "sbg:fileTypes": "VCF, VCF.GZ",
                        "secondaryFiles": [
                            {
                                "pattern": ".idx",
                                "required": false
                            },
                            {
                                "pattern": ".tbi",
                                "required": false
                            }
                        ]
                    },
                    {
                        "sbg:category": "Advanced Arguments",
                        "sbg:toolDefaultValue": "false",
                        "id": "allow_non_unique_kmers_in_ref",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--allow-non-unique-kmers-in-ref",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Allow non unique kmers in ref",
                        "doc": "Allow graphs that have non-unique kmers in the reference."
                    },
                    {
                        "sbg:category": "Optional Arguments",
                        "sbg:toolDefaultValue": "false",
                        "id": "annotate_with_num_discovered_alleles",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--annotate-with-num-discovered-alleles",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Annotate with num discovered alleles",
                        "doc": "If provided, we will annotate records with the number of alternate alleles that were discovered (but not necessarily genotyped) at a given site."
                    },
                    {
                        "sbg:altPrefix": "-A",
                        "sbg:category": "Optional Arguments",
                        "sbg:toolDefaultValue": "null",
                        "id": "annotation",
                        "type": "string[]?",
                        "inputBinding": {
                            "shellQuote": false,
                            "position": 4,
                            "valueFrom": "${\n    var cmd = \"\";\n    for (var i = 0; i < [].concat(inputs.annotation).length; i++ ){\n        cmd += \" --annotation \" + [].concat(inputs.annotation)[i];\n    }    \n    return cmd;\n}"
                        },
                        "label": "Annotation",
                        "doc": "One or more specific annotations to add to variant calls. This argument may be specified 0 or more times.\n\nPossible values: \nAlleleFraction, AS_BaseQualityRankSumTest,  AS_FisherStrand,  AS_InbreedingCoeff, AS_MappingQualityRankSumTest,  AS_QualByDepth,  AS_ReadPosRankSumTest,  AS_RMSMappingQuality, AS_StrandBiasMutectAnnotation,  AS_StrandOddsRatio,  BaseQuality,  BaseQualityHistogram, BaseQualityRankSumTest,  ChromosomeCounts,  ClippingRankSumTest,  CountNs,  Coverage,  DepthPerAlleleBySample,  DepthPerSampleHC,  ExcessHet,  FisherStrand,  FragmentLength, GenotypeSummaries,  InbreedingCoeff,  LikelihoodRankSumTest,  MappingQuality, MappingQualityRankSumTest,  MappingQualityZero,  OrientationBiasReadCounts, OriginalAlignment,  PossibleDeNovo,  QualByDepth,  ReadPosition,  ReadPosRankSumTest, ReferenceBases,  RMSMappingQuality,  SampleList,  StrandBiasBySample,  StrandOddsRatio, TandemRepeat,  UniqueAltReadCount"
                    },
                    {
                        "sbg:altPrefix": "-G",
                        "sbg:category": "Optional Arguments",
                        "sbg:toolDefaultValue": "null",
                        "id": "annotation_group",
                        "type": "string[]?",
                        "inputBinding": {
                            "shellQuote": false,
                            "position": 4,
                            "valueFrom": "${\n    var cmd = \"\";\n    for (var i = 0; i < [].concat(inputs.annotation_group).length; i++ ){\n        cmd += \" --annotation-group \" + [].concat(inputs.annotation_group)[i];\n    }    \n    return cmd;\n}"
                        },
                        "label": "Annotation group",
                        "doc": "One or more groups of annotations to apply to variant calls. Any requirements that are not met (e.g. failing to provide a pedigree file for a pedigree-based annotation) may cause the run to fail. This argument may be specified 0 or more times.\n\nPossible values: \nAlleleSpecificAnnotation,  AS_StandardAnnotation,  ReducibleAnnotation,  StandardAnnotation, StandardHCAnnotation,  StandardMutectAnnotation"
                    },
                    {
                        "sbg:altPrefix": "-AX",
                        "sbg:category": "Optional Arguments",
                        "sbg:toolDefaultValue": "null",
                        "id": "annotations_to_exclude",
                        "type": "string[]?",
                        "inputBinding": {
                            "shellQuote": false,
                            "position": 4,
                            "valueFrom": "${\n    var cmd = \"\";\n    for (var i = 0; i < [].concat(inputs.annotations_to_exclude).length; i++ ){\n        cmd += \" --annotations-to-exclude \" + [].concat(inputs.annotations_to_exclude)[i];\n    }    \n    return cmd;\n}"
                        },
                        "label": "Annotations to exclude",
                        "doc": "One or more specific annotations to exclude from variant calls. This argument may be specified 0 or more times. Which annotations to exclude from output in the variant calls. Note that this argument has higher priority than the -A or -G arguments, so these annotations will be excluded even if they are explicitly included with the other options.\n\nPossible values: \nBaseQualityRankSumTest, ChromosomeCounts, Coverage, DepthPerAlleleBySample, DepthPerSampleHC, ExcessHet, FisherStrand, InbreedingCoeff, MappingQualityRankSumTest, QualByDepth, ReadPosRankSumTest, RMSMappingQuality, StrandOddsRatio"
                    },
                    {
                        "sbg:category": "Optional Arguments",
                        "sbg:toolDefaultValue": "null",
                        "id": "assembly_region_out",
                        "type": "string?",
                        "inputBinding": {
                            "prefix": "--assembly-region-out",
                            "shellQuote": false,
                            "position": 4,
                            "valueFrom": "${\n    if(inputs.assembly_region_out) {\n        var tmp = inputs.assembly_region_out.slice(-4);\n        if(tmp == \".igv\" || tmp == \".IGV\") {\n            return inputs.assembly_region_out.slice(0, -4) + '.assembly.igv';\n        }\n        else {\n            return inputs.assembly_region_out + '.assembly.igv';\n        }\n    }\n    else {\n        return null;\n    }\n}"
                        },
                        "label": "Assembly region output",
                        "doc": "Name of the IGV file to which assembly region should be written."
                    },
                    {
                        "sbg:category": "Advanced Arguments",
                        "sbg:toolDefaultValue": "100",
                        "id": "assembly_region_padding",
                        "type": "int?",
                        "inputBinding": {
                            "prefix": "--assembly-region-padding",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Assembly region padding",
                        "doc": "Number of additional bases of context to include around each assembly region."
                    },
                    {
                        "sbg:altPrefix": "-bamout",
                        "sbg:category": "Advanced Arguments",
                        "sbg:toolDefaultValue": "null",
                        "id": "bam_output",
                        "type": "string?",
                        "inputBinding": {
                            "prefix": "--bam-output",
                            "shellQuote": false,
                            "position": 4,
                            "valueFrom": "${\n    if(inputs.bam_output) {\n        var tmp = inputs.bam_output.slice(-4);\n        if(tmp == \".bam\" || tmp == \".BAM\") {\n            return inputs.bam_output.slice(0, -4) + '.bam';\n        }\n        else {\n            return inputs.bam_output + '.bam';\n        }\n    }\n    else {\n        return null;\n    }\n}"
                        },
                        "label": "BAM output",
                        "doc": "Name of the file to which assembled haplotypes should be written."
                    },
                    {
                        "sbg:category": "Advanced Arguments",
                        "sbg:toolDefaultValue": "CALLED_HAPLOTYPES",
                        "id": "bam_writer_type",
                        "type": [
                            "null",
                            {
                                "type": "enum",
                                "symbols": [
                                    "ALL_POSSIBLE_HAPLOTYPES",
                                    "CALLED_HAPLOTYPES",
                                    "NO_HAPLOTYPES"
                                ],
                                "name": "bam_writer_type"
                            }
                        ],
                        "inputBinding": {
                            "prefix": "--bam-writer-type",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "BAM writer type",
                        "doc": "Which haplotypes should be written to the BAM."
                    },
                    {
                        "sbg:category": "Optional Arguments",
                        "sbg:toolDefaultValue": "18",
                        "id": "base_quality_score_threshold",
                        "type": "int?",
                        "inputBinding": {
                            "prefix": "--base-quality-score-threshold",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Base quality score threshold",
                        "doc": "Base qualities below this threshold will be reduced to the minimum (6)."
                    },
                    {
                        "sbg:altPrefix": "-comp",
                        "sbg:category": "Advanced Arguments",
                        "sbg:toolDefaultValue": "null",
                        "id": "comp",
                        "type": "File[]?",
                        "inputBinding": {
                            "shellQuote": false,
                            "position": 4,
                            "valueFrom": "${\n    var cmd = \"\";\n    for (var i = 0; i < [].concat(inputs.comp).length; i++ ){\n        cmd += \" --comparison \" + [].concat(inputs.comp)[i].path;\n    }    \n    return cmd;\n}"
                        },
                        "label": "Comparison VCF",
                        "doc": "Comparison vcf file(s). If a call overlaps with a record from the provided comp track, the INFO field will be annotated as such in the output with the track name. Records that are filtered in the comp track will be ignored. Note that 'dbSNP' has been special-cased (see the --dbsnp)",
                        "sbg:fileTypes": "VCF, VCF.GZ",
                        "secondaryFiles": [
                            {
                                "pattern": ".idx",
                                "required": false
                            },
                            {
                                "pattern": ".tbi",
                                "required": false
                            }
                        ]
                    },
                    {
                        "sbg:altPrefix": "-contamination-file",
                        "sbg:category": "Advanced Arguments",
                        "id": "contamination_fraction_per_sample_file",
                        "type": "File?",
                        "inputBinding": {
                            "prefix": "--contamination-fraction-per-sample-file",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Contamination fraction per sample",
                        "doc": "Tab-separated file containing fraction of contamination in sequencing data (per sample) to aggressively remove. Format should be \"<SampleID><TAB><Contamination>\" (Contamination is double) per line; No header.",
                        "sbg:fileTypes": "TSV"
                    },
                    {
                        "sbg:altPrefix": "-contamination",
                        "sbg:category": "Optional Arguments",
                        "sbg:toolDefaultValue": "0.0",
                        "id": "contamination_fraction_to_filter",
                        "type": "float?",
                        "inputBinding": {
                            "prefix": "--contamination-fraction-to-filter",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Contamination fraction to filter",
                        "doc": "Fraction of contamination in sequencing data (for all samples) to aggressively remove ."
                    },
                    {
                        "sbg:altPrefix": "-OBI",
                        "sbg:category": "Optional Arguments",
                        "sbg:toolDefaultValue": "true",
                        "id": "create_output_bam_index",
                        "type": [
                            "null",
                            {
                                "type": "enum",
                                "symbols": [
                                    "true",
                                    "false"
                                ],
                                "name": "create_output_bam_index"
                            }
                        ],
                        "inputBinding": {
                            "prefix": "--create-output-bam-index",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Create output BAM index",
                        "doc": "If true, create a BAM/CRAM index when writing a coordinate-sorted BAM/CRAM file."
                    },
                    {
                        "sbg:altPrefix": "-OVI",
                        "sbg:category": "Optional Arguments",
                        "sbg:toolDefaultValue": "true",
                        "id": "create_output_variant_index",
                        "type": [
                            "null",
                            {
                                "type": "enum",
                                "symbols": [
                                    "true",
                                    "false"
                                ],
                                "name": "create_output_variant_index"
                            }
                        ],
                        "inputBinding": {
                            "prefix": "--create-output-variant-index",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Create output variant index",
                        "doc": "If true, create a VCF index when writing a coordinate-sorted VCF file."
                    },
                    {
                        "sbg:altPrefix": "-D",
                        "sbg:category": "Optional Arguments",
                        "sbg:toolDefaultValue": "null",
                        "id": "dbsnp",
                        "type": "File?",
                        "inputBinding": {
                            "prefix": "--dbsnp",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "dbSNP",
                        "doc": "dbSNP file.",
                        "sbg:fileTypes": "VCF, VCF.GZ",
                        "secondaryFiles": [
                            {
                                "pattern": ".idx",
                                "required": false
                            },
                            {
                                "pattern": ".tbi",
                                "required": false
                            }
                        ]
                    },
                    {
                        "sbg:altPrefix": "-debug",
                        "sbg:category": "Advanced Arguments",
                        "sbg:toolDefaultValue": "false",
                        "id": "debug_assembly",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--debug-assembly",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Debug assembly",
                        "doc": "Print out verbose debug information about each assembly region."
                    },
                    {
                        "sbg:altPrefix": "-DBIC",
                        "sbg:category": "Optional Arguments",
                        "sbg:toolDefaultValue": "false",
                        "id": "disable_bam_index_caching",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--disable-bam-index-caching",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Disable BAM index caching",
                        "doc": "If true, don't cache BAM indexes, this will reduce memory requirements but may harm performance if many intervals are specified. Caching is automatically disabled if there are no intervals specified."
                    },
                    {
                        "sbg:category": "Advanced Arguments",
                        "sbg:toolDefaultValue": "false",
                        "id": "disable_optimizations",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--disable-optimizations",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Disable optimizations",
                        "doc": "Don't skip calculations in active regions with no variants."
                    },
                    {
                        "sbg:altPrefix": "-DF",
                        "sbg:category": "Optional Arguments",
                        "sbg:toolDefaultValue": "null",
                        "id": "disable_read_filter",
                        "type": "string[]?",
                        "inputBinding": {
                            "shellQuote": false,
                            "position": 4,
                            "valueFrom": "${\n    var cmd = \"\";\n    for (var i = 0; i < [].concat(inputs.disable_read_filter).length; i++ ){\n        cmd += \" --disable-read-filter \" + [].concat(inputs.disable_read_filter)[i];\n    }    \n    return cmd;\n}"
                        },
                        "label": "Disable read filter",
                        "doc": "Read filters to be disabled before analysis. This argument may be specified 0 or more times.\n\nPossible values:\nGoodCigarReadFilter, MappedReadFilter, MappingQualityAvailableReadFilter, MappingQualityReadFilter, NonZeroReferenceLengthAlignmentReadFilter, NotDuplicateReadFilter, NotSecondaryAlignmentReadFilter, PassesVendorQualityCheckReadFilter, WellformedReadFilter"
                    },
                    {
                        "sbg:altPrefix": "-disable-sequence-dictionary-validation",
                        "sbg:category": "Optional Arguments",
                        "sbg:toolDefaultValue": "false",
                        "id": "disable_sequence_dictionary_validation",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--disable-sequence-dictionary-validation",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Disable sequence dictionary validation",
                        "doc": "If specified, do not check the sequence dictionaries from our inputs for compatibility. Use at your own risk!"
                    },
                    {
                        "sbg:altPrefix": "-disable-tool-default-annotations",
                        "sbg:category": "Advanced Arguments",
                        "sbg:toolDefaultValue": "false",
                        "id": "disable_tool_default_annotations",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--disable-tool-default-annotations",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Disable tool default annotations",
                        "doc": "Disable all tool default annotations."
                    },
                    {
                        "sbg:altPrefix": "-disable-tool-default-read-filters",
                        "sbg:category": "Advanced Arguments",
                        "sbg:toolDefaultValue": "false",
                        "id": "disable_tool_default_read_filters",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--disable-tool-default-read-filters",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Disable tool default read filters",
                        "doc": "Disable all tool default read filters (warning: many tools will not function correctly without their default read filters on)."
                    },
                    {
                        "sbg:category": "Advanced Arguments",
                        "sbg:toolDefaultValue": "false",
                        "id": "do_not_run_physical_phasing",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--do-not-run-physical-phasing",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Do not run physical phasing",
                        "doc": "Disable physical phasing."
                    },
                    {
                        "sbg:category": "Advanced Arguments",
                        "sbg:toolDefaultValue": "false",
                        "id": "dont_increase_kmer_sizes_for_cycles",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--dont-increase-kmer-sizes-for-cycles",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Dont increase kmer sizes for cycles",
                        "doc": "Disable iterating over kmer sizes when graph cycles are detected."
                    },
                    {
                        "sbg:category": "Advanced Arguments",
                        "sbg:toolDefaultValue": "false",
                        "id": "dont_use_soft_clipped_bases",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--dont-use-soft-clipped-bases",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Do not use soft clipped bases",
                        "doc": "Do not analyze soft clipped bases in the reads."
                    },
                    {
                        "sbg:altPrefix": "-ERC",
                        "sbg:category": "Advanced Arguments",
                        "sbg:toolDefaultValue": "NONE",
                        "id": "emit_ref_confidence",
                        "type": [
                            "null",
                            {
                                "type": "enum",
                                "symbols": [
                                    "NONE",
                                    "BP_RESOLUTION",
                                    "GVCF"
                                ],
                                "name": "emit_ref_confidence"
                            }
                        ],
                        "inputBinding": {
                            "prefix": "--emit-ref-confidence",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Emit ref confidence",
                        "doc": "Mode for emitting reference confidence scores."
                    },
                    {
                        "sbg:category": "Advanced Arguments",
                        "sbg:toolDefaultValue": "false",
                        "id": "enable_all_annotations",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--enable-all-annotations",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Enable all annotations",
                        "doc": "Use all possible annotations (not for the faint of heart)."
                    },
                    {
                        "sbg:altPrefix": "-XL",
                        "sbg:category": "Optional Arguments",
                        "sbg:toolDefaultValue": "null",
                        "id": "exclude_intervals_string",
                        "type": "string[]?",
                        "inputBinding": {
                            "shellQuote": false,
                            "position": 4,
                            "valueFrom": "${\n    var cmd = \"\";\n    for (var i = 0; i < [].concat(inputs.exclude_intervals_string).length; i++ ){\n        cmd += \" --exclude-intervals \" + [].concat(inputs.exclude_intervals_string)[i];\n    }    \n    return cmd;\n}\n"
                        },
                        "label": "Exclude intervals string",
                        "doc": "One or more genomic intervals to exclude from processing. This argument may be specified 0 or more times."
                    },
                    {
                        "sbg:altPrefix": "-founder-id",
                        "sbg:category": "Optional Arguments",
                        "sbg:toolDefaultValue": "null",
                        "id": "founder_id",
                        "type": "string[]?",
                        "inputBinding": {
                            "shellQuote": false,
                            "position": 4,
                            "valueFrom": "${\n    var cmd = \"\";\n    for (var i = 0; i < [].concat(inputs.founder_id).length; i++ ){\n        cmd += \" --founder-id \" + [].concat(inputs.founder_id)[i];\n    }    \n    return cmd;\n}"
                        },
                        "label": "Founder ID",
                        "doc": "Samples representing the population \"founders\".  This argument may be specified 0 or more times."
                    },
                    {
                        "sbg:category": "Advanced Arguments",
                        "sbg:toolDefaultValue": "false",
                        "sbg:altPrefix": "-genotype-filtered-alleles",
                        "id": "force_call_filtered_alleles",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--force-call-filtered-alleles",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Force-call filtered alleles",
                        "doc": "Force-call filtered alleles included in the resource specified by --alleles."
                    },
                    {
                        "sbg:altPrefix": "-graph",
                        "sbg:category": "Optional Arguments",
                        "sbg:toolDefaultValue": "null",
                        "id": "graph_output",
                        "type": "string?",
                        "inputBinding": {
                            "prefix": "--graph-output",
                            "shellQuote": false,
                            "position": 4,
                            "valueFrom": "${\n    if(inputs.graph_output) {\n        var tmp = inputs.graph_output.slice(-4);\n        if(tmp == \".txt\" || tmp == \".TXT\") {\n            return inputs.graph_output.slice(0, -4) + '.txt';\n        }\n        else {\n            return inputs.graph_output + '.txt';\n        }\n    }\n    else {\n        return null;\n    }\n}"
                        },
                        "label": "Graph output",
                        "doc": "Name of the file to which debug assembly graph information should be written."
                    },
                    {
                        "sbg:altPrefix": "-GQB",
                        "sbg:category": "Advanced Arguments",
                        "sbg:toolDefaultValue": "1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 70, 80, 90, 99",
                        "id": "gvcf_gq_bands",
                        "type": "int[]?",
                        "inputBinding": {
                            "shellQuote": false,
                            "position": 4,
                            "valueFrom": "${\n    var cmd = \"\";\n    for (var i = 0; i < [].concat(inputs.gvcf_gq_bands).length; i++ ){\n        cmd += \" --gvcf-gq-bands \" + [].concat(inputs.gvcf_gq_bands)[i];\n    }    \n    return cmd;\n}"
                        },
                        "label": "GVCF GQ bands",
                        "doc": "Exclusive upper bounds for reference confidence GQ bands (must be in [1, 100] and specified in increasing order). This argument may be specified 0 or more times."
                    },
                    {
                        "sbg:category": "Optional Arguments",
                        "sbg:toolDefaultValue": "0.001",
                        "id": "heterozygosity",
                        "type": "float?",
                        "inputBinding": {
                            "prefix": "--heterozygosity",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Heterozygosity",
                        "doc": "Heterozygosity value used to compute prior likelihoods for any locus. See the GATKDocs for full details on the meaning of this population genetics concept."
                    },
                    {
                        "sbg:category": "Optional Arguments",
                        "sbg:toolDefaultValue": "0.01",
                        "id": "heterozygosity_stdev",
                        "type": "float?",
                        "inputBinding": {
                            "prefix": "--heterozygosity-stdev",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Heterozygosity stdev",
                        "doc": "Standard deviation of heterozygosity for SNP and indel calling."
                    },
                    {
                        "sbg:category": "Optional Arguments",
                        "sbg:toolDefaultValue": "1.25E-4",
                        "id": "indel_heterozygosity",
                        "type": "float?",
                        "inputBinding": {
                            "prefix": "--indel-heterozygosity",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Indel heterozygosity",
                        "doc": "Heterozygosity for indel calling. See the GATKDocs for heterozygosity for full details on the meaning of this population genetics concept."
                    },
                    {
                        "sbg:category": "Advanced Arguments",
                        "sbg:toolDefaultValue": "10",
                        "id": "indel_size_to_eliminate_in_ref_model",
                        "type": "int?",
                        "inputBinding": {
                            "prefix": "--indel-size-to-eliminate-in-ref-model",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Indel size to eliminate in ref model",
                        "doc": "The size of an indel to check for in the reference model."
                    },
                    {
                        "sbg:altPrefix": "-I",
                        "sbg:category": "Required Arguments",
                        "id": "in_alignments",
                        "type": "File[]",
                        "inputBinding": {
                            "shellQuote": false,
                            "position": 4,
                            "valueFrom": "${\n    var cmd = \"\";\n    for (var i = 0; i < [].concat(inputs.in_alignments).length; i++ ){\n        cmd += \" --input \" + [].concat(inputs.in_alignments)[i].path;\n    }    \n    return cmd;\n}"
                        },
                        "label": "Input alignments",
                        "doc": "BAM/SAM/CRAM file containing reads. This argument must be specified at least once.",
                        "sbg:fileTypes": "BAM, CRAM",
                        "secondaryFiles": [
                            {
                                "pattern": ".bai",
                                "required": false
                            },
                            {
                                "pattern": "^.bai",
                                "required": false
                            },
                            {
                                "pattern": ".crai",
                                "required": false
                            },
                            {
                                "pattern": "^.crai",
                                "required": false
                            }
                        ]
                    },
                    {
                        "sbg:altPrefix": "-ixp",
                        "sbg:category": "Optional Arguments",
                        "sbg:toolDefaultValue": "0",
                        "id": "interval_exclusion_padding",
                        "type": "int?",
                        "inputBinding": {
                            "prefix": "--interval-exclusion-padding",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Interval exclusion padding",
                        "doc": "Amount of padding (in bp) to add to each interval you are excluding."
                    },
                    {
                        "sbg:altPrefix": "-imr",
                        "sbg:category": "Optional Arguments",
                        "sbg:toolDefaultValue": "ALL",
                        "id": "interval_merging_rule",
                        "type": [
                            "null",
                            {
                                "type": "enum",
                                "symbols": [
                                    "ALL",
                                    "OVERLAPPING_ONLY"
                                ],
                                "name": "interval_merging_rule"
                            }
                        ],
                        "inputBinding": {
                            "prefix": "--interval-merging-rule",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Interval merging rule",
                        "doc": "Interval merging rule for abutting intervals."
                    },
                    {
                        "sbg:altPrefix": "-ip",
                        "sbg:category": "Optional Arguments",
                        "sbg:toolDefaultValue": "0",
                        "id": "interval_padding",
                        "type": "int?",
                        "inputBinding": {
                            "prefix": "--interval-padding",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Interval padding",
                        "doc": "Amount of padding (in bp) to add to each interval you are including."
                    },
                    {
                        "sbg:altPrefix": "-isr",
                        "sbg:category": "Optional Arguments",
                        "sbg:toolDefaultValue": "UNION",
                        "id": "interval_set_rule",
                        "type": [
                            "null",
                            {
                                "type": "enum",
                                "symbols": [
                                    "UNION",
                                    "INTERSECTION"
                                ],
                                "name": "interval_set_rule"
                            }
                        ],
                        "inputBinding": {
                            "prefix": "--interval-set-rule",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Interval set rule",
                        "doc": "Set merging approach to use for combining interval inputs."
                    },
                    {
                        "sbg:altPrefix": "-L",
                        "sbg:category": "Optional Arguments",
                        "sbg:toolDefaultValue": "null",
                        "id": "include_intervals_file",
                        "type": "File?",
                        "inputBinding": {
                            "prefix": "--intervals",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Include intervals file",
                        "doc": "One or more genomic intervals over which to operate.",
                        "sbg:fileTypes": "INTERVAL_LIST, LIST, BED"
                    },
                    {
                        "sbg:altPrefix": "-L",
                        "sbg:category": "Optional Arguments",
                        "sbg:toolDefaultValue": "null",
                        "id": "include_intervals_string",
                        "type": "string[]?",
                        "inputBinding": {
                            "shellQuote": false,
                            "position": 4,
                            "valueFrom": "${\n    var cmd = \"\";\n    for (var i = 0; i < [].concat(inputs.include_intervals_string).length; i++ ){\n        cmd += \" --intervals \" + [].concat(inputs.include_intervals_string)[i];\n    }    \n    return cmd;\n}\n"
                        },
                        "label": "Include intervals string",
                        "doc": "One or more genomic intervals over which to operate. This argument may be specified 0 or more times."
                    },
                    {
                        "sbg:category": "Advanced Arguments",
                        "sbg:toolDefaultValue": "10, 25",
                        "id": "kmer_size",
                        "type": "int[]?",
                        "inputBinding": {
                            "shellQuote": false,
                            "position": 4,
                            "valueFrom": "${\n    var cmd = \"\";\n    for (var i = 0; i < [].concat(inputs.kmer_size).length; i++ ){\n        cmd += \" --kmer-size \" + [].concat(inputs.kmer_size)[i];\n    }    \n    return cmd;\n}"
                        },
                        "label": "Kmer size",
                        "doc": "Kmer size to use in the read threading assembler. This argument may be specified 0 or more times."
                    },
                    {
                        "sbg:altPrefix": "-LE",
                        "sbg:category": "Optional Arguments",
                        "sbg:toolDefaultValue": "false",
                        "id": "lenient",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--lenient",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Lenient",
                        "doc": "Lenient processing of VCF files."
                    },
                    {
                        "sbg:category": "Advanced Arguments",
                        "sbg:toolDefaultValue": "6",
                        "id": "max_alternate_alleles",
                        "type": "int?",
                        "inputBinding": {
                            "prefix": "--max-alternate-alleles",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Max alternate alleles",
                        "doc": "Maximum number of alternate alleles to genotype."
                    },
                    {
                        "sbg:category": "Advanced Arguments",
                        "sbg:toolDefaultValue": "300",
                        "id": "max_assembly_region_size",
                        "type": "int?",
                        "inputBinding": {
                            "prefix": "--max-assembly-region-size",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Max assembly region size",
                        "doc": "Maximum size of an assembly region."
                    },
                    {
                        "sbg:category": "Advanced Arguments",
                        "sbg:toolDefaultValue": "1024",
                        "id": "max_genotype_count",
                        "type": "int?",
                        "inputBinding": {
                            "prefix": "--max-genotype-count",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Max genotype count",
                        "doc": "Maximum number of genotypes to consider at any site."
                    },
                    {
                        "sbg:altPrefix": "-mnp-dist",
                        "sbg:category": "Advanced Arguments",
                        "sbg:toolDefaultValue": "0",
                        "id": "max_mnp_distance",
                        "type": "int?",
                        "inputBinding": {
                            "prefix": "--max-mnp-distance",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Max MNP distance",
                        "doc": "Two or more phased substitutions separated by this distance or less are merged into MNPs."
                    },
                    {
                        "sbg:category": "Advanced Arguments",
                        "sbg:toolDefaultValue": "128",
                        "id": "max_num_haplotypes_in_population",
                        "type": "int?",
                        "inputBinding": {
                            "prefix": "--max-num-haplotypes-in-population",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Max num haplotypes in population",
                        "doc": "Maximum number of haplotypes to consider for your population."
                    },
                    {
                        "sbg:category": "Advanced Arguments",
                        "sbg:toolDefaultValue": "50",
                        "id": "max_prob_propagation_distance",
                        "type": "int?",
                        "inputBinding": {
                            "prefix": "--max-prob-propagation-distance",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Max prob propagation distance",
                        "doc": "Upper limit on how many bases away probability mass can be moved around when calculating the boundaries between active and inactive assembly regions."
                    },
                    {
                        "sbg:category": "Optional Arguments",
                        "sbg:toolDefaultValue": "50",
                        "id": "max_reads_per_alignment_start",
                        "type": "int?",
                        "inputBinding": {
                            "prefix": "--max-reads-per-alignment-start",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Max reads per alignment start",
                        "doc": "Maximum number of reads to retain per alignment start position. Reads above this threshold will be downsampled. Set to 0 to disable."
                    },
                    {
                        "sbg:category": "Advanced Arguments",
                        "sbg:toolDefaultValue": "100",
                        "id": "max_unpruned_variants",
                        "type": "int?",
                        "inputBinding": {
                            "prefix": "--max-unpruned-variants",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Max unpruned variants",
                        "doc": "Maximum number of variants in graph the adaptive pruner will allow."
                    },
                    {
                        "sbg:category": "Execution and Platform",
                        "sbg:toolDefaultValue": "100",
                        "id": "mem_overhead_per_job",
                        "type": "int?",
                        "label": "Memory overhead per job",
                        "doc": "It allows a user to set the desired overhead memory (in MB) when running a tool or adding it to a workflow."
                    },
                    {
                        "sbg:category": "Execution and Platform",
                        "sbg:toolDefaultValue": "4000",
                        "id": "mem_per_job",
                        "type": "int?",
                        "label": "Memory per job",
                        "doc": "It allows a user to set the desired memory requirement (in MB) when running a tool or adding it to a workflow."
                    },
                    {
                        "sbg:category": "Advanced Arguments",
                        "sbg:toolDefaultValue": "50",
                        "id": "min_assembly_region_size",
                        "type": "int?",
                        "inputBinding": {
                            "prefix": "--min-assembly-region-size",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Min assembly region size",
                        "doc": "Minimum size of an assembly region."
                    },
                    {
                        "sbg:altPrefix": "-mbq",
                        "sbg:category": "Optional Arguments",
                        "sbg:toolDefaultValue": "10",
                        "id": "min_base_quality_score",
                        "type": "int?",
                        "inputBinding": {
                            "prefix": "--min-base-quality-score",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Min base quality score",
                        "doc": "Minimum base quality required to consider a base for calling."
                    },
                    {
                        "sbg:category": "Advanced Arguments",
                        "sbg:toolDefaultValue": "4",
                        "id": "min_dangling_branch_length",
                        "type": "int?",
                        "inputBinding": {
                            "prefix": "--min-dangling-branch-length",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Min dangling branch length",
                        "doc": "Minimum length of a dangling branch to attempt recovery."
                    },
                    {
                        "sbg:category": "Advanced Arguments",
                        "sbg:toolDefaultValue": "2",
                        "id": "min_pruning",
                        "type": "int?",
                        "inputBinding": {
                            "prefix": "--min-pruning",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Min pruning",
                        "doc": "Minimum support to not prune paths in the graph."
                    },
                    {
                        "sbg:category": "Optional Arguments",
                        "sbg:toolDefaultValue": "4",
                        "id": "native_pair_hmm_threads",
                        "type": "int?",
                        "inputBinding": {
                            "prefix": "--native-pair-hmm-threads",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Native pairHMM threads",
                        "doc": "How many threads should a native pairHMM implementation use."
                    },
                    {
                        "sbg:category": "Optional Arguments",
                        "sbg:toolDefaultValue": "false",
                        "id": "native_pair_hmm_use_double_precision",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--native-pair-hmm-use-double-precision",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Native pairHMM use double precision",
                        "doc": "Use double precision in the native pairHMM. This is slower but matches the java implementation better."
                    },
                    {
                        "sbg:category": "Advanced Arguments",
                        "sbg:toolDefaultValue": "1",
                        "id": "num_pruning_samples",
                        "type": "int?",
                        "inputBinding": {
                            "prefix": "--num-pruning-samples",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Num pruning samples",
                        "doc": "Number of samples that must pass the minPruning threshold."
                    },
                    {
                        "sbg:category": "Optional Arguments",
                        "sbg:toolDefaultValue": "0",
                        "id": "num_reference_samples_if_no_call",
                        "type": "int?",
                        "inputBinding": {
                            "prefix": "--num-reference-samples-if-no-call",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Num reference samples if no call",
                        "doc": "Number of hom-ref genotypes to infer at sites not present in a panel."
                    },
                    {
                        "sbg:category": "Config Inputs",
                        "sbg:toolDefaultValue": "null",
                        "id": "output_prefix",
                        "type": "string?",
                        "label": "Output name prefix",
                        "doc": "Output file name prefix."
                    },
                    {
                        "sbg:category": "Optional Arguments",
                        "sbg:toolDefaultValue": "EMIT_VARIANTS_ONLY",
                        "id": "output_mode",
                        "type": [
                            "null",
                            {
                                "type": "enum",
                                "symbols": [
                                    "EMIT_VARIANTS_ONLY",
                                    "EMIT_ALL_CONFIDENT_SITES",
                                    "EMIT_ALL_ACTIVE_SITES"
                                ],
                                "name": "output_mode"
                            }
                        ],
                        "inputBinding": {
                            "prefix": "--output-mode",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Output mode",
                        "doc": "Specifies which type of calls we should output."
                    },
                    {
                        "sbg:category": "Advanced Arguments",
                        "sbg:toolDefaultValue": "10",
                        "id": "pair_hmm_gap_continuation_penalty",
                        "type": "int?",
                        "inputBinding": {
                            "prefix": "--pair-hmm-gap-continuation-penalty",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Pair HMM gap continuation penalty",
                        "doc": "Flat gap continuation penalty for use in the pairHMM."
                    },
                    {
                        "sbg:altPrefix": "-pairHMM",
                        "sbg:category": "Advanced Arguments",
                        "sbg:toolDefaultValue": "FASTEST_AVAILABLE",
                        "id": "pair_hmm_implementation",
                        "type": [
                            "null",
                            {
                                "type": "enum",
                                "symbols": [
                                    "EXACT",
                                    "ORIGINAL",
                                    "LOGLESS_CACHING",
                                    "AVX_LOGLESS_CACHING",
                                    "AVX_LOGLESS_CACHING_OMP",
                                    "EXPERIMENTAL_FPGA_LOGLESS_CACHING",
                                    "FASTEST_AVAILABLE"
                                ],
                                "name": "pair_hmm_implementation"
                            }
                        ],
                        "inputBinding": {
                            "prefix": "--pair-hmm-implementation",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Pair HMM implementation",
                        "doc": "The pairHMM implementation to use for genotype likelihood calculations."
                    },
                    {
                        "sbg:category": "Advanced Arguments",
                        "sbg:toolDefaultValue": "CONSERVATIVE",
                        "id": "pcr_indel_model",
                        "type": [
                            "null",
                            {
                                "type": "enum",
                                "symbols": [
                                    "NONE",
                                    "HOSTILE",
                                    "AGGRESSIVE",
                                    "CONSERVATIVE"
                                ],
                                "name": "pcr_indel_model"
                            }
                        ],
                        "inputBinding": {
                            "prefix": "--pcr-indel-model",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "PCR indel model",
                        "doc": "The PCR indel model to use."
                    },
                    {
                        "sbg:altPrefix": "-ped",
                        "sbg:category": "Optional Arguments",
                        "sbg:toolDefaultValue": "null",
                        "id": "pedigree",
                        "type": "File?",
                        "inputBinding": {
                            "prefix": "--pedigree",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Pedigree",
                        "doc": "Pedigree file for determining the population \"founders\".",
                        "sbg:fileTypes": "PED"
                    },
                    {
                        "sbg:category": "Advanced Arguments",
                        "sbg:toolDefaultValue": "45",
                        "id": "phred_scaled_global_read_mismapping_rate",
                        "type": "int?",
                        "inputBinding": {
                            "prefix": "--phred-scaled-global-read-mismapping-rate",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Phred scaled global read mismapping rate",
                        "doc": "The global assumed mismapping rate for reads."
                    },
                    {
                        "sbg:altPrefix": "-population",
                        "sbg:category": "Optional Arguments",
                        "sbg:toolDefaultValue": "null",
                        "id": "population_callset",
                        "type": "File?",
                        "inputBinding": {
                            "prefix": "--population-callset",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Population callset",
                        "doc": "Callset to use in calculating genotype priors.",
                        "sbg:fileTypes": "VCF, VCF.GZ",
                        "secondaryFiles": [
                            {
                                "pattern": ".idx",
                                "required": false
                            },
                            {
                                "pattern": ".tbi",
                                "required": false
                            }
                        ]
                    },
                    {
                        "sbg:category": "Advanced Arguments",
                        "sbg:toolDefaultValue": "2.302585092994046",
                        "id": "pruning_lod_threshold",
                        "type": "float?",
                        "inputBinding": {
                            "prefix": "--pruning-lod-threshold",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Pruning lod threshold",
                        "doc": "Ln likelihood ratio threshold for adaptive pruning algorithm."
                    },
                    {
                        "sbg:altPrefix": "-RF",
                        "sbg:category": "Optional Arguments",
                        "sbg:toolDefaultValue": "null",
                        "id": "read_filter",
                        "type": "string[]?",
                        "inputBinding": {
                            "shellQuote": false,
                            "position": 4,
                            "valueFrom": "${\n    var cmd = \"\";\n    for (var i = 0; i < [].concat(inputs.read_filter).length; i++ ){\n        cmd += \" --read-filter \" + [].concat(inputs.read_filter)[i];\n    }    \n    return cmd;\n}\n"
                        },
                        "label": "Read filter",
                        "doc": "Read filters to be applied before analysis. This argument may be specified 0 or more times.\n\nPossible values: \nAlignmentAgreesWithHeaderReadFilter, AllowAllReadsReadFilter, AmbiguousBaseReadFilter, CigarContainsNoNOperator, FirstOfPairReadFilter, FragmentLengthReadFilter, GoodCigarReadFilter, HasReadGroupReadFilter, IntervalOverlapReadFilter, LibraryReadFilter, MappedReadFilter, MappingQualityAvailableReadFilter, MappingQualityNotZeroReadFilter, MappingQualityReadFilter, MatchingBasesAndQualsReadFilter, MateDifferentStrandReadFilter, MateDistantReadFilter, MateOnSameContigOrNoMappedMateReadFilter, MateUnmappedAndUnmappedReadFilter, MetricsReadFilter, NonChimericOriginalAlignmentReadFilter, NonZeroFragmentLengthReadFilter, NonZeroReferenceLengthAlignmentReadFilter, NotDuplicateReadFilter, NotOpticalDuplicateReadFilter, NotProperlyPairedReadFilter, NotSecondaryAlignmentReadFilter, NotSupplementaryAlignmentReadFilter, OverclippedReadFilter, PairedReadFilter, PassesVendorQualityCheckReadFilter, PlatformReadFilter, PlatformUnitReadFilter, PrimaryLineReadFilter, ProperlyPairedReadFilter, ReadGroupBlackListReadFilter, ReadGroupReadFilter, ReadLengthEqualsCigarLengthReadFilter, ReadLengthReadFilter, ReadNameReadFilter, ReadStrandFilter, SampleReadFilter, SecondOfPairReadFilter, SeqIsStoredReadFilter, SoftClippedReadFilter, ValidAlignmentEndReadFilter, ValidAlignmentStartReadFilter, WellformedReadFilter"
                    },
                    {
                        "sbg:altPrefix": "-VS",
                        "sbg:category": "Optional Arguments",
                        "sbg:toolDefaultValue": "SILENT",
                        "id": "read_validation_stringency",
                        "type": [
                            "null",
                            {
                                "type": "enum",
                                "symbols": [
                                    "STRICT",
                                    "LENIENT",
                                    "SILENT"
                                ],
                                "name": "read_validation_stringency"
                            }
                        ],
                        "inputBinding": {
                            "prefix": "--read-validation-stringency",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Read validation stringency",
                        "doc": "Validation stringency for all SAM/BAM/CRAM/SRA files read by this program. The default stringency value silent can improve performance when processing a bam file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded."
                    },
                    {
                        "sbg:altPrefix": "-R",
                        "sbg:category": "Required Arguments",
                        "sbg:toolDefaultValue": "FASTA, FA",
                        "id": "in_reference",
                        "type": "File",
                        "inputBinding": {
                            "prefix": "--reference",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Reference",
                        "doc": "Reference sequence file.",
                        "sbg:fileTypes": "FASTA, FA",
                        "secondaryFiles": [
                            {
                                "pattern": ".fai",
                                "required": false
                            },
                            {
                                "pattern": "^.dict",
                                "required": false
                            }
                        ]
                    },
                    {
                        "sbg:altPrefix": "-ALIAS",
                        "sbg:category": "Optional Arguments",
                        "sbg:toolDefaultValue": "null",
                        "id": "sample_name",
                        "type": "string?",
                        "inputBinding": {
                            "prefix": "--sample-name",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Sample name",
                        "doc": "Name of single sample to use from a multi-sample bam."
                    },
                    {
                        "sbg:altPrefix": "-ploidy",
                        "sbg:category": "Optional Arguments",
                        "sbg:toolDefaultValue": "2",
                        "id": "sample_ploidy",
                        "type": "int?",
                        "inputBinding": {
                            "prefix": "--sample-ploidy",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Sample ploidy",
                        "doc": "Ploidy (number of chromosomes) per sample. For pooled data, set to (number of samples in each pool x Sample Ploidy)."
                    },
                    {
                        "sbg:altPrefix": "-sequence-dictionary",
                        "sbg:category": "Optional Arguments",
                        "id": "sequence_dictionary",
                        "type": "File?",
                        "inputBinding": {
                            "prefix": "--sequence-dictionary",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Sequence dictionary",
                        "doc": "Use the given sequence dictionary as the master/canonical sequence dictionary. Must be a .dict file.",
                        "sbg:fileTypes": "DICT"
                    },
                    {
                        "sbg:category": "Optional Arguments",
                        "sbg:toolDefaultValue": "false",
                        "id": "sites_only_vcf_output",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--sites-only-vcf-output",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Sites only VCF output",
                        "doc": "If true, don't emit genotype fields when writing VCF file output."
                    },
                    {
                        "sbg:category": "Advanced Arguments",
                        "sbg:toolDefaultValue": "JAVA",
                        "id": "smith_waterman",
                        "type": [
                            "null",
                            {
                                "type": "enum",
                                "symbols": [
                                    "FASTEST_AVAILABLE",
                                    "AVX_ENABLED",
                                    "JAVA"
                                ],
                                "name": "smith_waterman"
                            }
                        ],
                        "inputBinding": {
                            "prefix": "--smith-waterman",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Smith waterman",
                        "doc": "Which Smith-Waterman implementation to use, generally FASTEST_AVAILABLE is the right choice."
                    },
                    {
                        "sbg:altPrefix": "-stand-call-conf",
                        "sbg:category": "Optional Arguments",
                        "sbg:toolDefaultValue": "30.0",
                        "id": "standard_min_confidence_threshold_for_calling",
                        "type": "float?",
                        "inputBinding": {
                            "prefix": "--standard-min-confidence-threshold-for-calling",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Standard min confidence threshold for calling",
                        "doc": "The minimum phred-scaled confidence threshold at which variants should be called. Only variant sites with QUAL equal or greater than this threshold will be called. When HaplotypeCaller is used in GVCF mode (using either -ERC GVCF or -ERC BP_RESOLUTION) the call threshold is automatically set to zero. Call confidence thresholding will then be performed in the subsequent GenotypeGVCFs command."
                    },
                    {
                        "sbg:category": "Advanced Arguments",
                        "sbg:toolDefaultValue": "false",
                        "id": "use_filtered_reads_for_annotations",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--use-filtered-reads-for-annotations",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Use filtered reads for annotations",
                        "doc": "Use the contamination-filtered read maps for the purposes of annotating variants."
                    },
                    {
                        "sbg:altPrefix": "-XL",
                        "sbg:category": "Optional Arguments",
                        "id": "exclude_intervals_file",
                        "type": "File?",
                        "inputBinding": {
                            "prefix": "--exclude-intervals",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Exclude intervals file",
                        "doc": "One or more genomic intervals to exclude from processing.",
                        "sbg:fileTypes": "INTERVAL_LIST, LIST, BED"
                    },
                    {
                        "sbg:category": "Execution and Platform",
                        "sbg:toolDefaultValue": "1",
                        "id": "cpu_per_job",
                        "type": "int?",
                        "label": "CPU per job",
                        "doc": "Number of CPUs to be used per job."
                    },
                    {
                        "sbg:category": "Config Inputs",
                        "sbg:toolDefaultValue": "vcf.gz",
                        "id": "output_extension",
                        "type": [
                            "null",
                            {
                                "type": "enum",
                                "symbols": [
                                    "vcf",
                                    "vcf.gz"
                                ],
                                "name": "output_extension"
                            }
                        ],
                        "label": "Output VCF extension",
                        "doc": "Output VCF extension."
                    },
                    {
                        "sbg:toolDefaultValue": "10.00",
                        "sbg:category": "Optional Arguments",
                        "sbg:altPrefix": "-seconds-between-progress-updates",
                        "id": "seconds_between_progress_updates",
                        "type": "float?",
                        "inputBinding": {
                            "prefix": "--seconds-between-progress-updates",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Seconds between progress updates",
                        "doc": "Output traversal statistics every time this many seconds elapse."
                    },
                    {
                        "sbg:category": "Optional Arguments",
                        "sbg:altPrefix": "-read-index",
                        "id": "read_index",
                        "type": "File[]?",
                        "inputBinding": {
                            "shellQuote": false,
                            "position": 4,
                            "valueFrom": "${\n    var cmd = \"\";\n    for (var i = 0; i < [].concat(inputs.read_index).length; i++ ){\n        cmd += \" --read-index \" + [].concat(inputs.read_index)[i].path;\n    }    \n    return cmd;\n}"
                        },
                        "label": "Read index",
                        "doc": "Indices to use for the read inputs. If specified, an index must be provided for every read input and in the same order as the read inputs. If this argument is not specified, the path to the index for each input will be inferred automatically.",
                        "sbg:fileTypes": "BAI, CRAI"
                    },
                    {
                        "sbg:category": "Optional Arguments",
                        "sbg:toolDefaultValue": "false",
                        "id": "dont_use_dragstr_pair_hmm_scores",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--dont-use-dragstr-pair-hmm-scores",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Do not use DRAGstr pair hmm scores",
                        "doc": "Disable DRAGstr pair-hmm score even when dragstr-params-path was provided."
                    },
                    {
                        "sbg:toolDefaultValue": "false",
                        "sbg:category": "Optional Arguments",
                        "id": "dragen_mode",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--dragen-mode",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "DRAGEN mode",
                        "doc": "Single argument for enabling the bulk of DRAGEN-GATK features. NOTE: THIS WILL OVERWRITE PROVIDED ARGUMENT CHECK TOOL INFO TO SEE WHICH ARGUMENTS ARE SET)."
                    },
                    {
                        "sbg:toolDefaultValue": "2",
                        "sbg:category": "Optional Arguments",
                        "id": "dragstr_het_hom_ratio",
                        "type": "int?",
                        "inputBinding": {
                            "prefix": "--dragstr-het-hom-ratio",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "DRAGstr het hom ratio",
                        "doc": "Het to hom prior ratio use with DRAGstr on."
                    },
                    {
                        "sbg:category": "Optional Arguments",
                        "id": "dragstr_params_path",
                        "type": "File?",
                        "inputBinding": {
                            "prefix": "--dragstr-params-path",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "DRAGstr parameters",
                        "doc": "File with the DRAGstr model parameters for STR error correction used in the Pair HMM. When provided, it overrides other PCR error correcting mechanisms.",
                        "sbg:fileTypes": "TXT"
                    },
                    {
                        "sbg:toolDefaultValue": "false",
                        "sbg:category": "Optional Arguments",
                        "id": "enable_dynamic_read_disqualification_for_genotyping",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--enable-dynamic-read-disqualification-for-genotyping",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Enable dynamic read disqualification for genotyping",
                        "doc": "Will enable less strict read disqualification low base quality reads. If enabled, rather than disqualifying all reads over a threshold of minimum hmm scores we will instead choose a less strict and less aggressive cap for disqualification based on the read length and base qualities."
                    },
                    {
                        "sbg:category": "Optional Arguments",
                        "sbg:altPrefix": "-gam",
                        "sbg:toolDefaultValue": "USE_PLS_TO_ASSIGN",
                        "id": "genotype_assignment_method",
                        "type": [
                            "null",
                            {
                                "type": "enum",
                                "symbols": [
                                    "SET_TO_NO_CALL",
                                    "USE_PLS_TO_ASSIGN",
                                    "SET_TO_NO_CALL_NO_ANNOTATIONS",
                                    "BEST_MATCH_TO_ORIGINAL",
                                    "DO_NOT_ASSIGN_GENOTYPES",
                                    "USE_POSTERIOR_PROBABILITIES"
                                ],
                                "name": "genotype_assignment_method"
                            }
                        ],
                        "inputBinding": {
                            "prefix": "--genotype-assignment-method",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Genotype assignment method",
                        "doc": "How genotypes are assigned."
                    },
                    {
                        "sbg:toolDefaultValue": "false",
                        "sbg:category": "Optional Arguments",
                        "sbg:altPrefix": "-gp-qual",
                        "id": "use_posteriors_to_calculate_qual",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--use-posteriors-to-calculate-qual",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Use posteriors to calculate qual",
                        "doc": "If available, use the genotype posterior probabilities to calculate the site QUAL."
                    },
                    {
                        "sbg:toolDefaultValue": "2",
                        "sbg:category": "Advanced Arguments",
                        "id": "allele_informative_reads_overlap_margin",
                        "type": "int?",
                        "inputBinding": {
                            "prefix": "--allele-informative-reads-overlap-margin",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Allele informative reads overlap margin",
                        "doc": "Likelihood and read-based annotations will only take into consideration reads that overlap the variant or any base no further than this distance expressed in base pairs."
                    },
                    {
                        "sbg:toolDefaultValue": "false",
                        "sbg:category": "Advanced Arguments",
                        "id": "apply_bqd",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--apply-bqd",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Apply BQD",
                        "doc": "If enabled this argument will apply the DRAGEN-GATK BaseQualityDropout model to the genotyping model for filtering sites due to Linked Error mode."
                    },
                    {
                        "sbg:toolDefaultValue": "false",
                        "sbg:category": "Advanced Arguments",
                        "id": "apply_frd",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--apply-frd",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Apply FRD",
                        "doc": "If enabled this argument will apply the DRAGEN-GATK ForeignReadDetection model to the genotyping model for filtering sites."
                    },
                    {
                        "sbg:toolDefaultValue": "false",
                        "sbg:category": "Advanced Arguments",
                        "id": "disable_cap_base_qualities_to_map_quality",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--disable-cap-base-qualities-to-map-quality",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Disable cap base qualities to map quality",
                        "doc": "If false this disables capping of base qualities in the HMM to the mapping quality of the read."
                    },
                    {
                        "sbg:toolDefaultValue": "false",
                        "sbg:category": "Advanced Arguments",
                        "id": "disable_spanning_event_genotyping",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--disable-spanning-event-genotyping",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Disable spanning event genotyping",
                        "doc": "If enabled this argument will disable inclusion of the '*' spanning event when genotyping events that overlap deletions."
                    },
                    {
                        "sbg:toolDefaultValue": "false",
                        "sbg:category": "Advanced Arguments",
                        "id": "disable_symmetric_hmm_normalizing",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--disable-symmetric-hmm-normalizing",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Disable symmetric HMM normalizing",
                        "doc": "Toggle to revive legacy behavior of asymmetrically normalizing the arguments to the reference haplotype."
                    },
                    {
                        "sbg:toolDefaultValue": "false",
                        "sbg:category": "Advanced Arguments",
                        "id": "do_not_correct_overlapping_quality",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--do-not-correct-overlapping-quality",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Do not correct overlapping quality",
                        "doc": "Disable overlapping base quality correction. Base quality is capped at half of PCR error rate for bases where read and mate overlap, to account for full correlation of PCR errors at these bases. This argument disables that correction."
                    },
                    {
                        "sbg:toolDefaultValue": "false",
                        "sbg:category": "Advanced Arguments",
                        "id": "dont_use_dragstr_priors",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--dont-use-dragstr-priors",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Dont use dragstr priors",
                        "doc": "Forfeit the use of the DRAGstr model to calculate genotype priors. This argument does not have any effect in the absence of DRAGstr model parameters (--dragstr-model-params)."
                    },
                    {
                        "sbg:toolDefaultValue": "0.02",
                        "sbg:category": "Advanced Arguments",
                        "id": "expected_mismatch_rate_for_read_disqualification",
                        "type": "float?",
                        "inputBinding": {
                            "prefix": "--expected-mismatch-rate-for-read-disqualification",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Expected mismatch rate for read disqualification",
                        "doc": "Error rate used to set expectation for post HMM read disqualification based on mismatches."
                    },
                    {
                        "sbg:toolDefaultValue": "false",
                        "sbg:category": "Advanced Arguments",
                        "id": "floor_blocks",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--floor-blocks",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Floor blocks",
                        "doc": "Output the band lower bound for each GQ block regardless of the data it represents."
                    },
                    {
                        "sbg:toolDefaultValue": "false",
                        "sbg:category": "Advanced Arguments",
                        "id": "force_active",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--force-active",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Force active",
                        "doc": "If provided, all regions will be marked as active."
                    },
                    {
                        "sbg:toolDefaultValue": "false",
                        "sbg:category": "Advanced Arguments",
                        "id": "linked_de_bruijn_graph",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--linked-de-bruijn-graph",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Linked De Bruijn graph",
                        "doc": "If enabled, the Assembly Engine will construct a Linked De Bruijn graph to recover better haplotypes. Disables graph simplification into a seq graph, opts to construct a proper De Bruijn graph with potential loops NOTE: --linked-de-bruijn-graph is an experimental feature that does not directly match with the regular HaplotypeCaller. Specifically the haplotype finding code does not perform correctly at complicated sites. Use this mode at your own risk."
                    },
                    {
                        "sbg:toolDefaultValue": "20",
                        "sbg:category": "Advanced Arguments",
                        "id": "mapping_quality_threshold_for_genotyping",
                        "type": "int?",
                        "inputBinding": {
                            "prefix": "--mapping-quality-threshold-for-genotyping",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Mapping quality threshold for genotyping",
                        "doc": "Control the threshold for discounting reads from the genotyper due to mapping quality after the active region detection and assembly steps but before genotyping. NOTE: this is in contrast to the --minimum-mapping-quality argument which filters reads from all parts of the HaplotypeCaller. If you would like to call genotypes with a different threshold both arguments must be set."
                    },
                    {
                        "sbg:toolDefaultValue": "0",
                        "sbg:category": "Advanced Arguments",
                        "id": "max_effective_depth_adjustment_for_frd",
                        "type": "int?",
                        "inputBinding": {
                            "prefix": "--max-effective-depth-adjustment-for-frd",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Max effective depth adjustment for FRD",
                        "doc": "Set the maximum depth to modify FRD adjustment to in the event of high depth sites (0 to disable)."
                    },
                    {
                        "sbg:toolDefaultValue": "9.210340371976184",
                        "sbg:category": "Advanced Arguments",
                        "id": "pruning_seeding_lod_threshold",
                        "type": "float?",
                        "inputBinding": {
                            "prefix": "--pruning-seeding-lod-threshold",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Pruning seeding lod threshold",
                        "doc": "Ln likelihood ratio threshold for seeding subgraph of good variation in adaptive pruning algorithm."
                    },
                    {
                        "sbg:toolDefaultValue": "false",
                        "sbg:category": "Advanced Arguments",
                        "id": "recover_all_dangling_branches",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--recover-all-dangling-branches",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Recover all dangling branches",
                        "doc": "Recover all dangling branches. By default, the read threading assembler does not recover dangling branches that fork after splitting from the reference. This argument tells the assembly engine to recover all dangling branches."
                    },
                    {
                        "sbg:toolDefaultValue": "false",
                        "sbg:category": "Advanced Arguments",
                        "id": "soft_clip_low_quality_ends",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--soft-clip-low-quality-ends",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Soft clip low quality ends",
                        "doc": "If enabled will preserve low-quality read ends as softclips (used for DRAGEN-GATK BQD genotyper model)."
                    },
                    {
                        "sbg:toolDefaultValue": "false",
                        "sbg:category": "Advanced Arguments",
                        "id": "transform_dragen_mapping_quality",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--transform-dragen-mapping-quality",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Transform DRAGEN mapping quality",
                        "doc": "If enabled this argument will map DRAGEN aligner aligned reads with mapping quality <=250 to scale up to MQ 50."
                    },
                    {
                        "sbg:category": "Conditional Arguments for readFilter",
                        "sbg:toolDefaultValue": "null",
                        "id": "ambig_filter_bases",
                        "type": "int?",
                        "inputBinding": {
                            "prefix": "--ambig-filter-bases",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Ambig filter bases",
                        "doc": "Valid only if \"AmbiguousBaseReadFilter\" is specified:\nThreshold number of ambiguous bases. If null, uses threshold fraction; otherwise, overrides threshold fraction. Cannot be used in conjuction with argument(s) ambig-filter-frac."
                    },
                    {
                        "sbg:category": "Conditional Arguments for readFilter",
                        "sbg:toolDefaultValue": "0.05",
                        "id": "ambig_filter_frac",
                        "type": "float?",
                        "inputBinding": {
                            "prefix": "--ambig-filter-frac",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Ambig filter frac",
                        "doc": "Valid only if \"AmbiguousBaseReadFilter\" is specified:\nThreshold fraction of ambiguous bases. Cannot be used in conjuction with argument(s) ambig-filter-bases."
                    },
                    {
                        "sbg:category": "Conditional Arguments for readFilter",
                        "sbg:toolDefaultValue": "1000000",
                        "id": "max_fragment_length",
                        "type": "int?",
                        "inputBinding": {
                            "prefix": "--max-fragment-length",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Max fragment length",
                        "doc": "Valid only if \"FragmentLengthReadFilter\" is specified:\nMaximum length of fragment (insert size)."
                    },
                    {
                        "sbg:category": "Conditional Arguments for readFilter",
                        "sbg:toolDefaultValue": "0",
                        "id": "min_fragment_length",
                        "type": "int?",
                        "inputBinding": {
                            "prefix": "--min-fragment-length",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Min fragment length",
                        "doc": "Valid only if \"FragmentLengthReadFilter\" is specified:\nMinimum length of fragment (insert size)."
                    },
                    {
                        "sbg:category": "Conditional Arguments for readFilter",
                        "sbg:toolDefaultValue": "null",
                        "id": "keep_intervals",
                        "type": "string[]?",
                        "inputBinding": {
                            "prefix": "--keep-intervals",
                            "shellQuote": false,
                            "position": 4,
                            "valueFrom": "${\n    var cmd = \"\";\n    for (var i = 0; i < [].concat(inputs.keep_intervals).length; i++ ){\n        cmd += \" --keep-intervals \" + [].concat(inputs.keep_intervals)[i];\n    }    \n    return cmd;\n}\n\n\n"
                        },
                        "label": "Keep intervals",
                        "doc": "Valid only if \"IntervalOverlapReadFilter\" is specified:\nOne or more genomic intervals to keep. This argument must be specified at least once."
                    },
                    {
                        "sbg:altPrefix": "-library",
                        "sbg:category": "Conditional Arguments for readFilter",
                        "sbg:toolDefaultValue": "null",
                        "id": "library",
                        "type": "string[]?",
                        "inputBinding": {
                            "shellQuote": false,
                            "position": 4,
                            "valueFrom": "${\n    var cmd = \"\";\n    for (var i = 0; i < [].concat(inputs.library).length; i++ ){\n        cmd += \" --library \" + [].concat(inputs.library)[i];\n    }    \n    return cmd;\n}\n\n\n"
                        },
                        "label": "Library",
                        "doc": "Valid only if \"LibraryReadFilter\" is specified:\nName of the library to keep. This argument must be specified at least once. Required."
                    },
                    {
                        "sbg:category": "Conditional Arguments for readFilter",
                        "sbg:toolDefaultValue": "null",
                        "id": "maximum_mapping_quality",
                        "type": "int?",
                        "inputBinding": {
                            "prefix": "--maximum-mapping-quality",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Maximum mapping quality",
                        "doc": "Valid only if \"MappingQualityReadFilter\" is specified:\nMaximum mapping quality to keep (inclusive)."
                    },
                    {
                        "sbg:category": "Conditional Arguments for readFilter",
                        "sbg:toolDefaultValue": "10",
                        "id": "minimum_mapping_quality",
                        "type": "int?",
                        "inputBinding": {
                            "prefix": "--minimum-mapping-quality",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Minimum mapping quality",
                        "doc": "Valid only if \"MappingQualityReadFilter\" is specified:\nMinimum mapping quality to keep (inclusive)."
                    },
                    {
                        "sbg:toolDefaultValue": "1000",
                        "sbg:category": "Conditional Arguments for readFilter",
                        "id": "mate_too_distant_length",
                        "type": "int?",
                        "inputBinding": {
                            "prefix": "--mate-too-distant-length",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Mate too distant length",
                        "doc": "Valid only if \"MateDistantReadFilter\" is specified:\nMinimum start location difference at which mapped mates are considered distant."
                    },
                    {
                        "sbg:category": "Conditional Arguments for readFilter",
                        "sbg:toolDefaultValue": "false",
                        "id": "dont_require_soft_clips_both_ends",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--dont-require-soft-clips-both-ends",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Dont require soft clips both ends",
                        "doc": "Valid only if \"OverclippedReadFilter\" is specified:\nAllow a read to be filtered out based on having only 1 soft-clipped block. By default, both ends must have a soft-clipped block, setting this flag requires only 1 soft-clipped block."
                    },
                    {
                        "sbg:category": "Conditional Arguments for readFilter",
                        "sbg:toolDefaultValue": "30",
                        "id": "filter_too_short",
                        "type": "int?",
                        "inputBinding": {
                            "prefix": "--filter-too-short",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Filter too short",
                        "doc": "Valid only if \"OverclippedReadFilter\" is specified:\nMinimum number of aligned bases."
                    },
                    {
                        "sbg:category": "Conditional Arguments for readFilter",
                        "sbg:toolDefaultValue": "null",
                        "id": "platform_filter_name",
                        "type": "string[]?",
                        "inputBinding": {
                            "shellQuote": false,
                            "position": 4,
                            "valueFrom": "${\n    var cmd = \"\";\n    for (var i = 0; i < [].concat(inputs.platform_filter_name).length; i++ ){\n        cmd += \" --platform-filter-name \" + [].concat(inputs.platform_filter_name)[i];\n    }    \n    return cmd;\n}\n\n\n"
                        },
                        "label": "Platform filter name",
                        "doc": "Valid only if \"PlatformReadFilter\" is specified:\nPlatform attribute (PL) to match. This argument must be specified at least once. Required."
                    },
                    {
                        "sbg:category": "Conditional Arguments for readFilter",
                        "sbg:toolDefaultValue": "null",
                        "id": "black_listed_lanes",
                        "type": "string[]?",
                        "inputBinding": {
                            "shellQuote": false,
                            "position": 4,
                            "valueFrom": "${\n    var cmd = \"\";\n    for (var i = 0; i < [].concat(inputs.black_listed_lanes).length; i++ ){\n        cmd += \" --black-listed-lanes \" + [].concat(inputs.black_listed_lanes)[i];\n    }    \n    return cmd;\n}"
                        },
                        "label": "Black listed lanes",
                        "doc": "Valid only if \"PlatformUnitReadFilter\" is specified:\nPlatform unit (PU) to filter out. This argument must be specified at least once. Required."
                    },
                    {
                        "sbg:category": "Conditional Arguments for readFilter",
                        "sbg:toolDefaultValue": "null",
                        "id": "read_group_black_list",
                        "type": "string[]?",
                        "inputBinding": {
                            "shellQuote": false,
                            "position": 4,
                            "valueFrom": "${\n    var cmd = \"\";\n    for (var i = 0; i < [].concat(inputs.read_group_black_list).length; i++ ){\n        cmd += \" --read-group-black-list \" + [].concat(inputs.read_group_black_list)[i];\n    }    \n    return cmd;\n}\n\n\n"
                        },
                        "label": "Read group black list",
                        "doc": "Valid only if \"ReadGroupBlackListReadFilter\" is specified:\nA read group filter expression in the form \"attribute:value\", where \"attribute\" is a two character read group attribute such as \"RG\" or \"PU\". This argument must be specified at least once. Required."
                    },
                    {
                        "sbg:category": "Conditional Arguments for readFilter",
                        "sbg:toolDefaultValue": "null",
                        "id": "keep_read_group",
                        "type": "string?",
                        "inputBinding": {
                            "prefix": "--keep-read-group",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Keep read group",
                        "doc": "Valid only if \"ReadGroupReadFilter\" is specified:\nThe name of the read group to keep. Required."
                    },
                    {
                        "sbg:category": "Conditional Arguments for readFilter",
                        "sbg:toolDefaultValue": "null",
                        "id": "max_read_length",
                        "type": "int?",
                        "inputBinding": {
                            "prefix": "--max-read-length",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Max read length",
                        "doc": "Valid only if \"ReadLengthReadFilter\" is specified:\nKeep only reads with length at most equal to the specified value. Required."
                    },
                    {
                        "sbg:category": "Conditional Arguments for readFilter",
                        "sbg:toolDefaultValue": "1",
                        "id": "min_read_length",
                        "type": "int?",
                        "inputBinding": {
                            "prefix": "--min-read-length",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Min read length",
                        "doc": "Valid only if \"ReadLengthReadFilter\" is specified:\nKeep only reads with length at least equal to the specified value."
                    },
                    {
                        "sbg:category": "Conditional Arguments for readFilter",
                        "sbg:toolDefaultValue": "null",
                        "id": "read_name",
                        "type": "string?",
                        "inputBinding": {
                            "prefix": "--read-name",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Read name",
                        "doc": "Valid only if \"ReadNameReadFilter\" is specified:\nKeep only reads with this read name. Required."
                    },
                    {
                        "sbg:category": "Conditional Arguments for readFilter",
                        "sbg:toolDefaultValue": "null",
                        "id": "keep_reverse_strand_only",
                        "type": [
                            "null",
                            {
                                "type": "enum",
                                "symbols": [
                                    "true",
                                    "false"
                                ],
                                "name": "keep_reverse_strand_only"
                            }
                        ],
                        "inputBinding": {
                            "prefix": "--keep-reverse-strand-only",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Keep reverse strand only",
                        "doc": "Valid only if \"ReadStrandFilter\" is specified:\nKeep only reads on the reverse strand. Required."
                    },
                    {
                        "sbg:altPrefix": "-sample",
                        "sbg:category": "Conditional Arguments for readFilter",
                        "sbg:toolDefaultValue": "null",
                        "id": "sample",
                        "type": "string[]?",
                        "inputBinding": {
                            "shellQuote": false,
                            "position": 4,
                            "valueFrom": "${\n    var cmd = \"\";\n    for (var i = 0; i < [].concat(inputs.sample).length; i++ ){\n        cmd += \" --sample \" + [].concat(inputs.sample)[i];\n    }    \n    return cmd;\n}\n\n\n"
                        },
                        "label": "Sample",
                        "doc": "Valid only if \"SampleReadFilter\" is specified:\nThe name of the sample(s) to keep, filtering out all others. This argument must be specified at least once. Required."
                    },
                    {
                        "sbg:toolDefaultValue": "false",
                        "sbg:category": "Conditional Arguments for readFilter",
                        "id": "invert_soft_clip_ratio_filter",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--invert-soft-clip-ratio-filter",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Invert soft clip ratio filter",
                        "doc": "Valid only if \"SoftClippedReadFilter\" is specified:\nInverts the results from this filter, causing all variants that would pass to fail and visa-versa."
                    },
                    {
                        "sbg:toolDefaultValue": "null",
                        "sbg:category": "Conditional Arguments for readFilter",
                        "id": "soft_clipped_leading_trailing_ratio",
                        "type": "float?",
                        "inputBinding": {
                            "prefix": "--soft-clipped-leading-trailing-ratio",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Soft clipped leading trailing ratio",
                        "doc": "Valid only if \"SoftClippedReadFilter\" is specified:\nThreshold ratio of soft clipped bases (leading / trailing the cigar string) to total bases in read for read to be filtered. Cannot be used in conjunction with argument(s) soft-clipped-ratio-threshold."
                    },
                    {
                        "sbg:toolDefaultValue": "null",
                        "sbg:category": "Conditional Arguments for readFilter",
                        "id": "soft_clipped_ratio_threshold",
                        "type": "float?",
                        "inputBinding": {
                            "prefix": "--soft-clipped-ratio-threshold",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Soft clipped ratio threshold",
                        "doc": "Valid only if \"SoftClippedReadFilter\" is specified:\nThreshold ratio of soft clipped bases (anywhere in the cigar string) to total bases in read for read to be filtered.  Cannot be used in conjunction with argument(s) soft-clipped-leading-trailing-ratio."
                    },
                    {
                        "sbg:category": "Optional Arguments",
                        "sbg:toolDefaultValue": "false",
                        "sbg:altPrefix": "-OVM",
                        "id": "create_output_variant_md5",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--create-output-variant-md5",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Create output variant md5",
                        "doc": "If true, create a a MD5 digest any VCF file created."
                    },
                    {
                        "sbg:toolDefaultValue": "false",
                        "sbg:category": "Optional Arguments",
                        "sbg:altPrefix": "-OBM",
                        "id": "create_output_bam_md5",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--create-output-bam-md5",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Create output bam md5",
                        "doc": "If true, create a MD5 digest for any BAM/SAM/CRAM file created."
                    },
                    {
                        "sbg:toolDefaultValue": "true",
                        "sbg:category": "Optional Arguments",
                        "sbg:altPrefix": "-new-qual",
                        "id": "use_new_qual_calculator",
                        "type": [
                            "null",
                            {
                                "type": "enum",
                                "symbols": [
                                    "true",
                                    "false"
                                ],
                                "name": "use_new_qual_calculator"
                            }
                        ],
                        "inputBinding": {
                            "prefix": "--use-new-qual-calculator",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Use new qual calculator",
                        "doc": "Use the new AF model instead of the so-called exact model."
                    },
                    {
                        "sbg:toolDefaultValue": "false",
                        "sbg:category": "Optional Arguments",
                        "id": "allow_old_rms_mapping_quality_annotation_data",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--allow-old-rms-mapping-quality-annotation-data",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Allow old rms mapping quality annotation data",
                        "doc": "Valid only if \"RMSMappingQuality\" is specified:\nOverride to allow old RMSMappingQuality annotated VCFs to function."
                    }
                ],
                "outputs": [
                    {
                        "id": "out_variants",
                        "doc": "File to which variants should be written.",
                        "label": "VCF output",
                        "type": "File",
                        "outputBinding": {
                            "glob": "${\n    var output_ext = inputs.output_extension ? inputs.output_extension : \"vcf.gz\";\n    return \"*.\" + output_ext;\n}",
                            "outputEval": "$(inheritMetadata(self, inputs.in_alignments))"
                        },
                        "secondaryFiles": [
                            {
                                "pattern": ".idx",
                                "required": false
                            },
                            {
                                "pattern": ".tbi",
                                "required": false
                            },
                            {
                                "pattern": ".md5",
                                "required": false
                            }
                        ],
                        "sbg:fileTypes": "VCF, VCF.GZ"
                    }
                ],
                "doc": "**GATK HaplotypeCaller** calls germline SNPs and indels from input BAM file(s) via local re-assembly of haplotypes [1].\n\n**GATK HaplotypeCaller** is capable of calling SNPs and indels simultaneously via local de-novo assembly of haplotypes in an active region. In other words, whenever the program encounters a region showing signs of variation, it discards the existing mapping information and completely reassembles the reads in that region. This allows the **GATK HaplotypeCaller** to be more accurate when calling regions that are traditionally difficult to call, for example when they contain different types of variants close to each other. It also makes the **GATK HaplotypeCaller** much better at calling indels than position-based callers like UnifiedGenotyper [1].\n\nIn the GVCF workflow used for scalable variant calling in DNA sequence data, **GATK HaplotypeCaller** runs per-sample to generate an intermediate GVCF (not to be used in final analysis), which can then be used in GenotypeGVCFs for joint genotyping of multiple samples in a very efficient way. The GVCF workflow enables rapid incremental processing of samples as they roll off the sequencer, as well as scaling to very large cohort sizes [1].\n\nIn addition, **HaplotypeCaller** is able to handle non-diploid organisms as well as pooled experiment data. Note however that the algorithms used to calculate variant likelihoods are not well suited to extreme allele frequencies (relative to ploidy) so its use is not recommended for somatic (cancer) variant discovery. For that purpose, use **Mutect2** instead [1].\n\nFinally, **GATK HaplotypeCaller** is also able to correctly handle splice junctions that make RNAseq a challenge for most variant callers, on the condition that the input read data has previously been processed according to [GATK RNAseq short variant discovery (SNPs + Indels)](https://gatk.broadinstitute.org/hc/en-us/articles/360035531192?id=4067) [1].\n\n*A list of **all inputs and parameters** with corresponding descriptions can be found at the bottom of this page.*\n\n***Please note that any cloud infrastructure costs resulting from app and pipeline executions, including the use of public apps, are the sole responsibility of you as a user. To avoid excessive costs, please read the app description carefully and set the app parameters and execution settings accordingly.***\n\n### Common Use Cases\n\n- Call variants individually on each sample in GVCF mode\n\n```\n gatk --java-options \"-Xmx4g\" HaplotypeCaller  \\\n   -R Homo_sapiens_assembly38.fasta \\\n   -I input.bam \\\n   -O output.g.vcf.gz \\\n   -ERC GVCF\n```\n\n\n- Call variants individually on each sample in GVCF mode with allele-specific annotations. [Here](https://gatk.broadinstitute.org/hc/en-us/articles/360035890551?id=9622) you can read more details about allele-specific annotation and filtering.\n\n```\ngatk --java-options \"-Xmx4g\" HaplotypeCaller  \\\n   -R Homo_sapiens_assembly38.fasta \\\n   -I input.bam \\\n   -O output.g.vcf.gz \\\n   -ERC GVCF \\\n   -G Standard \\\n   -G AS_Standard\n```\n\n\n- Call variants with bamout to show realigned reads.\n\n```\n gatk --java-options \"-Xmx4g\" HaplotypeCaller  \\\n   -R Homo_sapiens_assembly38.fasta \\\n   -I input.bam \\\n   -O output.vcf.gz \\\n   -bamout bamout.bam\n```\n\n### Changes Introduced by Seven Bridges\n\n* **Include intervals** (`--intervals`) option is divided into **Include intervals file** and **Include intervals string** options.\n* **Exclude intervals** (`--exclude-intervals`) option is divided into **Exclude intervals file** and **Exclude intervals string** options.\n* **VCF output** will be prefixed using the **Output name prefix** parameter. If this value is not set, the output name will be generated based on the **Sample ID** metadata value from **Input alignments**. If the **Sample ID** value is not set, the name will be inherited from the **Input alignments** file name. In case there are multiple files on the **Input alignments** input, the files will be sorted by name and output file name will be generated based on the first file in the sorted file list, following the rules defined in the previous case. \n* The user can specify the output file format using the **Output VCF extension** argument. Otherwise, the output will be in the compressed VCF file format.\n* The following parameters were excluded from the tool wrapper: `--arguments_file`, `--cloud-index-prefetch-buffer`, `--cloud-prefetch-buffer`, `--gatk-config-file`, `--gcs-max-retries`, `--gcs-project-for-requester-pays`, `--help`, `--QUIET`, `--recover-dangling-heads` (deprecated), `--showHidden`, `--tmp-dir`, `--use-jdk-deflater`, `--use-jdk-inflater`, `--verbosity`, `--version`\n\n### Common Issues and Important Notes\n\n*  **Memory per job** (`mem_per_job`) input allows a user to set the desired memory requirement when running a tool or adding it to a workflow. This input should be defined in MB. It is propagated to the Memory requirements part and “-Xmx” parameter of the tool. The default value is 4000 MB.\n* **Memory overhead per job** (`mem_overhead_per_job`) input allows a user to set the desired overhead memory when running a tool or adding it to a workflow. This input should be defined in MB. This amount will be added to the **Memory per job** in the Memory requirements section but it will not be added to the “-Xmx” parameter. The default value is 100 MB. \n* Note: GATK tools that take in mapped read data expect a BAM file as the primary format [2]. More on GATK requirements for mapped sequence data formats can be found [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035890791-SAM-or-BAM-or-CRAM-Mapped-sequence-data-formats).\n* Note: **Alleles**, **Comparison VCF**, **dbSNP**, **Input alignments**, **Population callset** should have corresponding index files in the same folder. \n* Note: **Reference** FASTA file should have corresponding .fai (FASTA index) and .dict (FASTA dictionary) files in the same folder. \n* Note: When working with PCR-free data, be sure to set **PCR indel model** (`--pcr_indel_model`) to NONE [1].\n* Note: When running **Emit ref confidence** ( `--emit-ref-confidence`) in GVCF or in BP_RESOLUTION mode, the confidence threshold is automatically set to 0. This cannot be overridden by the command line. The threshold can be set manually to the desired level when using **GenotypeGVCFs** [1].\n* Note: It is recommended to use a list of intervals to speed up the analysis. See [this document](https://gatk.broadinstitute.org/hc/en-us/articles/360035889551?id=4133) for details [1].\n* Note: **HaplotypeCaller** is able to handle many non-diploid use cases; the desired ploidy can be specified using the `-ploidy` argument. Note however that very high ploidies (such as are encountered in large pooled experiments) may cause performance challenges including excessive slowness [1].\n* Note: These **Read Filters** (`--read-filter`) are automatically applied to the data by the Engine before processing by **HaplotypeCaller** [1]: **NotSecondaryAlignmentReadFilter**, **GoodCigarReadFilter**, **NonZeroReferenceLengthAlignmentReadFilter**, **PassesVendorQualityCheckReadFilter**, **MappedReadFilter**, **MappingQualityAvailableReadFilter**, **NotDuplicateReadFilter**, **MappingQualityReadFilter**, **WellformedReadFilter**\n* Note: If the **Read filter** (`--read-filter`) option is set to \"LibraryReadFilter\", the **Library** (`--library`) option must be set to some value.\n* Note: If the **Read filter** (`--read-filter`) option is set to \"PlatformReadFilter\", the **Platform filter name** (`--platform-filter-name`) option must be set to some value.\n* Note: If the **Read filter** (`--read-filter`) option is set to \"PlatformUnitReadFilter\", the **Black listed lanes** (`--black-listed-lanes`) option must be set to some value. \n* Note: If the **Read filter** (`--read-filter`) option is set to \"ReadGroupBlackListReadFilter\", the **Read group black list** (`--read-group-black-list`) option must be set to some value.\n* Note: If the **Read filter** (`--read-filter`) option is set to \"ReadGroupReadFilter\", the **Keep read group** (`--keep-read-group`) option must be set to some value.\n* Note: If the **Read filter** (`--read-filter`) option is set to \"ReadLengthReadFilter\", the **Max read length** (`--max-read-length`) option must be set to some value.\n* Note: If the **Read filter** (`--read-filter`) option is set to \"ReadNameReadFilter\", the **Read name** (`--read-name`) option must be set to some value.\n* Note: If the **Read filter** (`--read-filter`) option is set to \"ReadStrandFilter\", the **Keep reverse strand only** (`--keep-reverse-strand-only`) option must be set to some value.\n* Note: If the **Read filter** (`--read-filter`) option is set to \"SampleReadFilter\", the **Sample** (`--sample`) option must be set to some value.\n* Note: If the **Read filter** (`--read-filter`) option is set to \"IntervalOverlapReadFilter\", the **Keep intervals** (`--keep-intervals`) option must be set to some value.\n* Note: The following options are valid only if the appropriate **Read filter** (`--read-filter`) is specified: **Ambig filter bases** (`--ambig-filter-bases`), **Ambig filter frac** (`--ambig-filter-frac`), **Max fragment length** (`--max-fragment-length`), **Min fragment length** (`--min-fragment-length`), **Keep intervals** (`--keep-intervals`), **Library** (`--library`), **Maximum mapping quality** (`--maximum-mapping-quality`), **Minimum mapping quality** (`--minimum-mapping-quality`),  **Mate too distant length** (`--mate-too-distant-length`), **Do not require soft clips** (`--dont-require-soft-clips-both-ends`), **Filter too short** (`--filter-too-short`), **Platform filter name** (`--platform-filter-name`), **Black listed lanes** (`--black-listed-lanes`), **Read group black list** (`--read-group-black-list`), **Keep read group** (`--keep-read-group`), **Max read length** (`--max-read-length`), **Min read length** (`--min-read-length`), **Read name** (`--read-name`), **Keep reverse strand only** (`--keep-reverse-strand-only`), **Sample** (`--sample`), **Invert soft clip ratio filter** (`--invert-soft-clip-ratio-filter`), **Soft clipped leading trailing ratio** (`--soft-clipped-leading-trailing-ratio`), **Soft clipped ratio threshold** (`--soft-clipped-ratio-threshold`) . See the description of each parameter for information on the associated **Read filter**.\n* Note: Allele-specific annotations are not yet supported in the VCF mode.\n* Note: The wrapper has not been tested for the SAM file type on the **Input alignments** input port.\n* Note: DRAGEN-GATK features have not been tested. Once the full DRAGEN-GATK pipeline is released, parameters related to DRAGEN-GATK mode will be tested. \n\n### Performance Benchmarking\n\nBelow is a table describing the runtimes and task costs for a couple of samples with different file sizes.\n\n| Experiment type |  Input size | Paired-end | # of reads | Read length | Duration | Cost (on-demand) | AWS instance type |\n|:--------------:|:------------:|:--------:|:-------:|:---------:|:----------:|:------:|:------:|:------:|\n|     RNA-Seq     | 2.5 GB |     Yes    |     16M     |     101     |   52min   | 0.47$ | c4.2xlarge |\n|     RNA-Seq     | 7.5 GB |     Yes    |     50M     |     101     |   1h46min   | 0.95$ | c4.2xlarge |\n|     RNA-Seq     | 12.5 GB |     Yes    |     82M    |     101     |  2h40min  | 1.43$ | c4.2xlarge |\n|     RNA-Seq     | 24.5 GB |     Yes    |     164M    |     101     |  4h55min  | 2.64$ | c4.2xlarge |\n\n*Cost can be significantly reduced by using **spot instances**. Visit the [knowledge center](https://docs.sevenbridges.com/docs/about-spot-instances) for more details.*\n\n###Portability\n\n**GATK HaplotypeCaller** is tested with cwltool version: \"3.0.20201203173111\"\n\n### References\n[1] [GATK HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360056969012-HaplotypeCaller)\n\n[2] [GATK Mapped sequence data formats](https://gatk.broadinstitute.org/hc/en-us/articles/360035890791-SAM-or-BAM-or-CRAM-Mapped-sequence-data-formats)",
                "label": "GATK HaplotypeCaller",
                "arguments": [
                    {
                        "prefix": "",
                        "shellQuote": false,
                        "position": 1,
                        "valueFrom": "${\n    if (inputs.mem_per_job) {\n        return '\\\"-Xmx'.concat(inputs.mem_per_job, 'M') + '\\\"';\n    } else {\n        return '\\\"-Xmx4000M\\\"';\n    }\n}"
                    },
                    {
                        "prefix": "",
                        "shellQuote": false,
                        "position": 2,
                        "valueFrom": "HaplotypeCaller"
                    },
                    {
                        "prefix": "--output",
                        "shellQuote": false,
                        "position": 3,
                        "valueFrom": "${\n  var output_name =  \"\";\n  var count = \"\";\n  var input_files = [].concat(inputs.in_alignments);\n  var tmp_ext = inputs.output_extension ? inputs.output_extension : \"vcf.gz\";\n  \n  if(inputs.emit_ref_confidence == 'GVCF' || inputs.emit_ref_confidence == 'BP_RESOLUTION'){\n      var output_ext = 'g.' + tmp_ext;\n  } \n  else {\n      var output_ext = tmp_ext;\n      \n  }\n  \n  var extension = \".\" + output_ext;\n    \n  if (input_files.length > 1){\n        count = \".\".concat(input_files.length);\n  }\n  if (inputs.output_prefix){\n    output_name = inputs.output_prefix;\n  } else {\n    if (input_files.length > 1){\n        input_files.sort(function(file1, file2) {\n            var file1_name = file1.basename.toUpperCase();\n            var file2_name = file2.basename.toUpperCase();\n            if (file1_name < file2_name) {\n                return -1;\n            }\n            if (file1_name > file2_name) {\n                return 1;\n            }\n            return 0;\n        });\n    }\n      \n    var sample_id = \"\";\n    var in_first_file = input_files[0];\n    \n    if (in_first_file.metadata && in_first_file.metadata.sample_id){\n      sample_id = in_first_file.metadata.sample_id;\n      }\n    if (sample_id){\n      output_name = sample_id;\n    } else {\n      output_name = in_first_file.basename.split('.')[0];\n    }\n  }\n  \n  // ensure there are no special characters in output_name\n  output_name = output_name.replace(/[^a-zA-Z0-9\\_\\-\\.]/g, \"\");\n  \n  return output_name + count + extension;\n}"
                    }
                ],
                "requirements": [
                    {
                        "class": "ShellCommandRequirement"
                    },
                    {
                        "class": "ResourceRequirement",
                        "ramMin": "${\n  var memory = 4000;\n  \n  if(inputs.mem_per_job) {\n  \t memory = inputs.mem_per_job;\n  }\n  if(inputs.mem_overhead_per_job) {\n\tmemory += inputs.mem_overhead_per_job;\n  }\n  else {\n      memory += 100;\n  }\n  return memory;\n}",
                        "coresMin": "${\n    return inputs.cpu_per_job ? inputs.cpu_per_job : 1;\n}"
                    },
                    {
                        "class": "DockerRequirement",
                        "dockerPull": "images.sbgenomics.com/marijeta_slavkovic/gatk-4-2-0-0:0"
                    },
                    {
                        "class": "InitialWorkDirRequirement",
                        "listing": []
                    },
                    {
                        "class": "InlineJavascriptRequirement",
                        "expressionLib": [
                            "\nvar setMetadata = function(file, metadata) {\n    if (!('metadata' in file)) {\n        file['metadata'] = {}\n    }\n    for (var key in metadata) {\n        file['metadata'][key] = metadata[key];\n    }\n    return file\n};\nvar inheritMetadata = function(o1, o2) {\n    var commonMetadata = {};\n    if (!o2) {\n        return o1;\n    };\n    if (!Array.isArray(o2)) {\n        o2 = [o2]\n    }\n    for (var i = 0; i < o2.length; i++) {\n        var example = o2[i]['metadata'];\n        for (var key in example) {\n            if (i == 0)\n                commonMetadata[key] = example[key];\n            else {\n                if (!(commonMetadata[key] == example[key])) {\n                    delete commonMetadata[key]\n                }\n            }\n        }\n        for (var key in commonMetadata) {\n            if (!(key in example)) {\n                delete commonMetadata[key]\n            }\n        }\n    }\n    if (!Array.isArray(o1)) {\n        o1 = setMetadata(o1, commonMetadata)\n        if (o1.secondaryFiles) {\n            o1.secondaryFiles = inheritMetadata(o1.secondaryFiles, o2)\n        }\n    } else {\n        for (var i = 0; i < o1.length; i++) {\n            o1[i] = setMetadata(o1[i], commonMetadata)\n            if (o1[i].secondaryFiles) {\n                o1[i].secondaryFiles = inheritMetadata(o1[i].secondaryFiles, o2)\n            }\n        }\n    }\n    return o1;\n};"
                        ]
                    }
                ],
                "sbg:categories": [
                    "Variant Calling",
                    "CWLtool Tested"
                ],
                "sbg:image_url": null,
                "sbg:license": "Apache License 2.0",
                "sbg:links": [
                    {
                        "id": "https://gatk.broadinstitute.org/hc/en-us",
                        "label": "Homepage"
                    },
                    {
                        "id": "https://github.com/broadinstitute/gatk",
                        "label": "Source Code"
                    },
                    {
                        "id": "https://github.com/broadinstitute/gatk/releases/download/4.2.0.0/gatk-4.2.0.0.zip",
                        "label": "Download"
                    },
                    {
                        "id": "https://www.biorxiv.org/content/10.1101/201178v3",
                        "label": "Publication"
                    },
                    {
                        "id": "https://gatk.broadinstitute.org/hc/en-us/articles/360056969012-HaplotypeCaller",
                        "label": "Documentation"
                    }
                ],
                "sbg:projectName": "SBG Public data",
                "sbg:revisionsInfo": [
                    {
                        "sbg:revision": 0,
                        "sbg:modifiedBy": "admin",
                        "sbg:modifiedOn": 1634729387,
                        "sbg:revisionNotes": null
                    },
                    {
                        "sbg:revision": 1,
                        "sbg:modifiedBy": "admin",
                        "sbg:modifiedOn": 1634729387,
                        "sbg:revisionNotes": "initial copy"
                    },
                    {
                        "sbg:revision": 2,
                        "sbg:modifiedBy": "admin",
                        "sbg:modifiedOn": 1634729387,
                        "sbg:revisionNotes": "description portability"
                    },
                    {
                        "sbg:revision": 3,
                        "sbg:modifiedBy": "admin",
                        "sbg:modifiedOn": 1634729388,
                        "sbg:revisionNotes": "changes after review process"
                    },
                    {
                        "sbg:revision": 4,
                        "sbg:modifiedBy": "admin",
                        "sbg:modifiedOn": 1648047573,
                        "sbg:revisionNotes": ""
                    },
                    {
                        "sbg:revision": 5,
                        "sbg:modifiedBy": "admin",
                        "sbg:modifiedOn": 1648047573,
                        "sbg:revisionNotes": "categories"
                    },
                    {
                        "sbg:revision": 6,
                        "sbg:modifiedBy": "admin",
                        "sbg:modifiedOn": 1681122217,
                        "sbg:revisionNotes": "JS naming for assembly_region_out, bam_output and graph_output⁩"
                    }
                ],
                "sbg:toolAuthor": "Broad Institute",
                "sbg:toolkit": "GATK",
                "sbg:toolkitVersion": "4.2.0.0",
                "sbg:expand_workflow": false,
                "sbg:appVersion": [
                    "v1.2"
                ],
                "sbg:id": "admin/sbg-public-data/gatk-haplotypecaller-4-2-0-0/6",
                "sbg:revision": 6,
                "sbg:revisionNotes": "JS naming for assembly_region_out, bam_output and graph_output⁩",
                "sbg:modifiedOn": 1681122217,
                "sbg:modifiedBy": "admin",
                "sbg:createdOn": 1634729387,
                "sbg:createdBy": "admin",
                "sbg:project": "admin/sbg-public-data",
                "sbg:sbgMaintained": false,
                "sbg:validationErrors": [],
                "sbg:contributors": [
                    "admin"
                ],
                "sbg:latestRevision": 6,
                "sbg:publisher": "sbg",
                "sbg:content_hash": "a62178128744f3c6a8ba6b4d259c653fa28a09eeafc1f489321b994e1bb064473",
                "sbg:workflowLanguage": "CWL"
            },
            "label": "HaplotypeCaller",
            "sbg:x": 298.0980529785156,
            "sbg:y": 217.3660125732422
        },
        {
            "id": "bcftools_consensus",
            "in": [
                {
                    "id": "input_file_uncompressed",
                    "source": "gatk_haplotypecaller_4_2_0_0/out_variants"
                },
                {
                    "id": "reference",
                    "source": "reference"
                }
            ],
            "out": [
                {
                    "id": "output_file"
                }
            ],
            "run": {
                "cwlVersion": "sbg:draft-2",
                "class": "CommandLineTool",
                "$namespaces": {
                    "sbg": "https://sevenbridges.com"
                },
                "id": "admin/sbg-public-data/bcftools-consensus/6",
                "label": "Bcftools Consensus",
                "description": "**BCFtools Consensus**: Create a consensus sequence by applying VCF variants to a reference FASTA file. \n\n\n**BCFtools** is a set of utilities that manipulate variant calls in the Variant Call Format (VCF) and its binary counterpart BCF. All commands work transparently with both VCFs and BCFs, both uncompressed and BGZF-compressed. Most commands accept VCF, bgzipped VCF and BCF with filetype detected automatically even when streaming from a pipe. Indexed VCF and BCF will work in all situations. Un-indexed VCF and BCF and streams will work in most, but not all situations. In general, whenever multiple VCFs are read simultaneously, they must be indexed and therefore also compressed. [1]\n\nA list of **all inputs and parameters** with corresponding descriptions can be found at the bottom of the page.\n\n\n### Common Use Cases\n\nBy default, the program will apply all ALT variants to the reference FASTA to obtain the consensus sequence. Can be also used with the **SAMtools Faidx** tool where the desired part of reference can be extracted and then provided  to this tool.\n\n\nUsing the **Sample** (`--sample`) (and, optionally, **Haplotype** (`--haplotype`) option will apply genotype (haplotype) calls from FORMAT/GT. \n```\n$bcftools consensus -s NA001 -f in.fa in.vcf.gz > out.fa\n```\n\nApply variants present in sample \"NA001\", output IUPAC codes using **Output in IUPAC** (`--iupac-codes`) option\n```\nbcftools consensus --iupac-codes -s NA001 -f in.fa in.vcf.gz > out.fa\n```\n\n### Changes Introduced by Seven Bridges\n\n* BCFtools works in all cases with gzipped and indexed VCF/BCF files. To be sure BCFtools works in all cases, we added subsequent `bgzip` and `index` commands if an uncompressed VCF or BCF file is provided on the **Uncompressed VCF/BCF file** input. If VCF.GZ, BCF.GZ or compressed BCF is given on the **Compressed VCF/BCF file** input only indexing will be done. Output type can still be chosen with the `output type` command.\n\n### Common Issues and Important Notes\n\n* Compressed input files must be provided on the **Compressed VCF/BCF file** input and uncompressed files on the **Uncompressed VCF/BCF file** input. Otherwise, the tool will fail as it will try to compress an already compressed input.\n\n* It is required to provide either **Uncompressed VCF/BCF file** or **Compressed VCF/BCF file**. If both inputs or none are provided the tool will fail.\n\n* By default, the program will apply all ALT variants to the reference FASTA to obtain the consensus sequence. \n\n * If the FASTA sequence does not match the REF allele at a given position, the tool will fail.\n\n### Performance Benchmarking\n\nIt took 5 minutes to execute this tool on AWS c4.2xlarge instance with a 56 KB VCF and a 3 GB reference FASTA file. The price is negligible ($0.02).\n\n*Cost can be significantly reduced by using **spot instances**. Visit the [Knowledge Center](https://docs.sevenbridges.com/docs/about-spot-instances) for more details.*\n\n### References\n[1 - BCFtools page](https://samtools.github.io/bcftools/bcftools.html)",
                "baseCommand": [
                    {
                        "class": "Expression",
                        "engine": "#cwl-js-engine",
                        "script": "{\n    if($job.inputs.input_file_compressed){\n        var in_files_array = [].concat($job.inputs.input_file_compressed);\n        var cmd = ''\n        for (var i=0; i < in_files_array.length; i++){\n            var in_file = in_files_array[i];\n            var fname = in_file.path.replace(/^.*[\\\\\\/]/, '');\n            var fname_ext = fname.split('.').pop();\n            if (fname.split('.').pop() == 'gz'){\n                var froot = fname.split('.');\n                froot.pop()\n                froot = froot.join('.')\n            }\n            else{\n                froot = fname\n            }\n            \n                \n            if(froot.split('.').pop() == 'bcf'){\n                cmd += \"bcftools index  -f -c \" + in_file.path + \" && \";\n                }\n            else{\n                cmd += \"bcftools index  -f -t \" + in_file.path + \" && \";\n                }\n        }\n    }\n    else{\n        var in_files_array = [].concat($job.inputs.input_file_uncompressed);\n        var cmd = ''\n        for (var i=0; i < in_files_array.length; i++){\n            var in_file = in_files_array[i];\n            var fname = in_file.path.replace(/^.*[\\\\\\/]/, '');\n            var fname_ext = fname.split('.').pop();\n            if(fname_ext == 'bcf'){\n                cmd += \"bgzip -c -f \" + fname + \" > \" + fname + \".gz\" + \" && bcftools index -f -c \" + fname + \".gz && \";\n            }\n            if(fname_ext == 'vcf'){ \n                cmd += \"bgzip -c -f \" + fname + \" > \" + fname + \".gz\" + \" && bcftools index -f -t \" + fname + \".gz && \";\n        \n            }\n        }\n    }\n    return cmd\n}"
                    },
                    "bcftools",
                    "consensus"
                ],
                "inputs": [
                    {
                        "sbg:category": "File Input",
                        "sbg:stageInput": "link",
                        "type": [
                            "null",
                            "File"
                        ],
                        "inputBinding": {
                            "position": 40,
                            "valueFrom": {
                                "class": "Expression",
                                "engine": "#cwl-js-engine",
                                "script": "{\n  if($self){\n  fname = $job.inputs.input_file_uncompressed.path.replace(/^.*[\\\\\\/]/, '')\n    return fname + \".gz\"\n  }\n}"
                            },
                            "sbg:cmdInclude": true
                        },
                        "label": "Uncompressed VCF/BCF file",
                        "description": "Uncompressed VCF/BCF file.",
                        "sbg:fileTypes": "VCF, BCF",
                        "id": "#input_file_uncompressed"
                    },
                    {
                        "sbg:stageInput": "link",
                        "sbg:category": "Input",
                        "type": [
                            "null",
                            "File"
                        ],
                        "inputBinding": {
                            "position": 40
                        },
                        "label": "Compressed VCF/BCF file",
                        "description": "Compressed VCF/BCF file.",
                        "sbg:fileTypes": "VCF.GZ, BCF, BCF.GZ",
                        "id": "#input_file_compressed"
                    },
                    {
                        "sbg:category": "General Options",
                        "sbg:altPrefix": "-i",
                        "type": [
                            "null",
                            "string"
                        ],
                        "inputBinding": {
                            "position": 5,
                            "prefix": "--include",
                            "separate": true,
                            "sbg:cmdInclude": true
                        },
                        "label": "Include expression",
                        "description": "Include only sites for which the expression is true.",
                        "id": "#include_expression"
                    },
                    {
                        "sbg:category": "General Options",
                        "sbg:altPrefix": "-e",
                        "type": [
                            "null",
                            "string"
                        ],
                        "inputBinding": {
                            "position": 1,
                            "prefix": "--exclude",
                            "separate": true,
                            "sbg:cmdInclude": true
                        },
                        "label": "Exclude expression",
                        "description": "Exclude sites for which the expression is true.",
                        "id": "#exclude_expression"
                    },
                    {
                        "sbg:category": "Configuration",
                        "type": [
                            "null",
                            "string"
                        ],
                        "label": "Output file name",
                        "description": "Name of the output file.",
                        "id": "#output_name"
                    },
                    {
                        "sbg:category": "General Options",
                        "sbg:altPrefix": "-s",
                        "type": [
                            "null",
                            "string"
                        ],
                        "inputBinding": {
                            "position": 20,
                            "prefix": "--sample",
                            "separate": true,
                            "sbg:cmdInclude": true
                        },
                        "label": "Sample",
                        "description": "Apply variants of the given sample.",
                        "id": "#sample"
                    },
                    {
                        "sbg:category": "Execution",
                        "sbg:stageInput": null,
                        "sbg:toolDefaultValue": "1",
                        "type": [
                            "null",
                            "int"
                        ],
                        "label": "Number of CPUs",
                        "description": "Number of CPUs. Appropriate instance will be chosen based on this parameter.",
                        "id": "#cpu"
                    },
                    {
                        "sbg:category": "Execution",
                        "sbg:stageInput": null,
                        "sbg:toolDefaultValue": "1000",
                        "type": [
                            "null",
                            "int"
                        ],
                        "label": "Memory in MB",
                        "description": "Memory in MB. Appropriate instance will be chosen based on this parameter.",
                        "id": "#memory"
                    },
                    {
                        "sbg:category": "General Options",
                        "sbg:altPrefix": "-c",
                        "type": [
                            "null",
                            "string"
                        ],
                        "inputBinding": {
                            "position": 0,
                            "prefix": "--chain",
                            "separate": true,
                            "sbg:cmdInclude": true
                        },
                        "label": "Create a chain file",
                        "description": "Create a chain file for liftover. Write a chain file filename including extension (eg. input_file.chain).",
                        "id": "#chain"
                    },
                    {
                        "sbg:altPrefix": "-f",
                        "type": [
                            "File"
                        ],
                        "inputBinding": {
                            "position": 2,
                            "prefix": "--fasta-ref",
                            "separate": true,
                            "sbg:cmdInclude": true
                        },
                        "label": "Reference Genome",
                        "description": "Reference sequence in fasta format.",
                        "sbg:fileTypes": "FASTA",
                        "id": "#reference"
                    },
                    {
                        "sbg:category": "General Options",
                        "sbg:altPrefix": "-h",
                        "type": [
                            "null",
                            {
                                "type": "enum",
                                "symbols": [
                                    "1",
                                    "2",
                                    "R",
                                    "A",
                                    "LR,LA",
                                    "SR,SA"
                                ],
                                "name": "haplotype"
                            }
                        ],
                        "inputBinding": {
                            "position": 4,
                            "prefix": "--haplotype",
                            "separate": true,
                            "valueFrom": {
                                "script": "{\n\n  if($job.inputs.haplotype == '1'){return \"1\"}\n  if($job.inputs.haplotype == '2'){return \"2\"}\n  if($job.inputs.haplotype == 'R'){return \"R\"}\n  if($job.inputs.haplotype == 'A'){return \"A\"}\n  if($job.inputs.haplotype == 'LR,LA'){return \"LR,LA\"}\n  if($job.inputs.haplotype == 'SR,SA'){return \"SR,SA\"}\n\n\n\n\n\n}",
                                "class": "Expression",
                                "engine": "#cwl-js-engine"
                            },
                            "sbg:cmdInclude": true
                        },
                        "label": "Haplotype",
                        "description": "Choose which allele to use from the FORMAT/GT field.",
                        "id": "#haplotype"
                    },
                    {
                        "sbg:category": "General Options",
                        "sbg:altPrefix": "-I",
                        "type": [
                            "null",
                            "boolean"
                        ],
                        "inputBinding": {
                            "position": 6,
                            "prefix": "--iupac-codes",
                            "separate": true,
                            "sbg:cmdInclude": true
                        },
                        "label": "Output in IUPAC",
                        "description": "Output variants in the form of IUPAC ambiguity codes.",
                        "id": "#iupac"
                    },
                    {
                        "sbg:category": "General Options",
                        "sbg:altPrefix": "-m",
                        "type": [
                            "null",
                            "File"
                        ],
                        "inputBinding": {
                            "position": 7,
                            "prefix": "--mask",
                            "separate": true,
                            "sbg:cmdInclude": true
                        },
                        "label": "Mask file",
                        "description": "Replace regions with N.",
                        "id": "#mask"
                    },
                    {
                        "sbg:category": "General Options",
                        "sbg:altPrefix": "-M",
                        "type": [
                            "null",
                            "string"
                        ],
                        "inputBinding": {
                            "position": 8,
                            "prefix": "--missing",
                            "separate": true,
                            "sbg:cmdInclude": true
                        },
                        "label": "Missing genotypes",
                        "description": "Output <char> instead of skipping the missing genotypes.",
                        "id": "#missing"
                    }
                ],
                "outputs": [
                    {
                        "type": [
                            "null",
                            "File"
                        ],
                        "label": "Output file",
                        "description": "Consensus sequence.",
                        "sbg:fileTypes": "FA",
                        "outputBinding": {
                            "glob": {
                                "class": "Expression",
                                "engine": "#cwl-js-engine",
                                "script": "{\n  if($job.inputs.input_file_uncompressed){\n      var array_files = [].concat($job.inputs.input_file_uncompressed)\n      fname = array_files[0].path.replace(/^.*[\\\\\\/]/, '')\n  }\n  else{\n      var array_files = [].concat($job.inputs.input_file_compressed)\n      fname = array_files[0].path.replace(/^.*[\\\\\\/]/, '')\n      if(fname.split('.').pop().toLowerCase() == 'gz'){ \n        fname = array_files[0].path.replace(/^.*[\\\\\\/]/, '')\n        fname = fname.replace(/\\.[^/.]+$/, \"\")\n      }\n  }\n  \n  fname_list = fname.split('.')\n  fname_list.pop() // Remove extension\n  out = fname_list.join('.')\n  if ($job.inputs.output_name){\n      return $job.inputs.output_name + \".fa\"\n  }\n  else {\n      return out + \".fa\"\n  }\n}"
                            },
                            "outputEval": {
                                "class": "Expression",
                                "engine": "#cwl-js-engine",
                                "script": "{\n    var setMetadata = function(file, metadata) {\n    if (!('metadata' in file)) {\n        file['metadata'] = {}\n    }\n    for (var key in metadata) {\n        file['metadata'][key] = metadata[key];\n    }\n    return file\n    };\n\n    var inheritMetadata = function(o1, o2) {\n        var commonMetadata = {};\n        if (!o2) {\n            return o1;\n        };\n        if (!Array.isArray(o2)) {\n            o2 = [o2]\n        }\n        for (var i = 0; i < o2.length; i++) {\n            var example = o2[i]['metadata'];\n            for (var key in example) {\n                if (i == 0)\n                    commonMetadata[key] = example[key];\n                else {\n                    if (!(commonMetadata[key] == example[key])) {\n                        delete commonMetadata[key]\n                    }\n                }\n            }\n            for (var key in commonMetadata) {\n                if (!(key in example)) {\n                    delete commonMetadata[key]\n                }\n            }\n        }\n        if (!Array.isArray(o1)) {\n            o1 = setMetadata(o1, commonMetadata)\n            if (o1.secondaryFiles) {\n                o1.secondaryFiles = inheritMetadata(o1.secondaryFiles, o2)\n            }\n        } else {\n            for (var i = 0; i < o1.length; i++) {\n                o1[i] = setMetadata(o1[i], commonMetadata)\n                if (o1[i].secondaryFiles) {\n                    o1[i].secondaryFiles = inheritMetadata(o1[i].secondaryFiles, o2)\n                }\n            }\n        }\n        return o1;\n    };\n    \n    if($job.inputs.input_file_compressed){\n        temp = [].concat($self)\n        for (i in temp){\n            temp[i] = inheritMetadata(temp[i], $job.inputs.input_file_compressed)\n        }\n    }\n    else{\n        temp = [].concat($self)\n        for (i in temp){\n            temp[i] = inheritMetadata(temp[i], $job.inputs.input_file_uncompressed)\n        }\n    }\n    return temp\n \n}"
                            }
                        },
                        "id": "#output_file"
                    }
                ],
                "requirements": [
                    {
                        "class": "ExpressionEngineRequirement",
                        "id": "#cwl-js-engine",
                        "requirements": [
                            {
                                "dockerPull": "rabix/js-engine",
                                "class": "DockerRequirement"
                            }
                        ]
                    }
                ],
                "hints": [
                    {
                        "class": "sbg:CPURequirement",
                        "value": {
                            "class": "Expression",
                            "engine": "#cwl-js-engine",
                            "script": "{\n  if($job.inputs.cpu){\n    return $job.inputs.cpu}\n  else{\n    return 1}\n}"
                        }
                    },
                    {
                        "class": "sbg:MemRequirement",
                        "value": {
                            "class": "Expression",
                            "engine": "#cwl-js-engine",
                            "script": "{\n  if($job.inputs.memory){\n    return $job.inputs.memory}\n  else{\n    return 1000}\n}    "
                        }
                    },
                    {
                        "class": "DockerRequirement",
                        "dockerImageId": "21caaa02f72e",
                        "dockerPull": "images.sbgenomics.com/luka_topalovic/bcftools:1.9"
                    }
                ],
                "arguments": [
                    {
                        "position": 3,
                        "prefix": "--output",
                        "separate": true,
                        "valueFrom": {
                            "class": "Expression",
                            "engine": "#cwl-js-engine",
                            "script": "{\n  if($job.inputs.input_file_uncompressed){\n      var array_files = [].concat($job.inputs.input_file_uncompressed)\n      fname = array_files[0].path.replace(/^.*[\\\\\\/]/, '')\n  }\n  else{\n      var array_files = [].concat($job.inputs.input_file_compressed)\n      fname = array_files[0].path.replace(/^.*[\\\\\\/]/, '')\n      if(fname.split('.').pop().toLowerCase() == 'gz'){ \n        fname = array_files[0].path.replace(/^.*[\\\\\\/]/, '')\n        fname = fname.replace(/\\.[^/.]+$/, \"\")\n      }\n  }\n  \n  fname_list = fname.split('.')\n  fname_list.pop() // Remove extension\n  out = fname_list.join('.')\n  if ($job.inputs.output_name){\n      return $job.inputs.output_name + \".fa\"\n  }\n  else {\n      return out + \".fa\"\n  }\n}"
                        }
                    }
                ],
                "sbg:projectName": "SBG Public data",
                "sbg:toolkitVersion": "1.9",
                "sbg:toolAuthor": "Petr Danecek, Shane McCarthy, John Marshall",
                "sbg:categories": [
                    "VCF Processing"
                ],
                "sbg:links": [
                    {
                        "id": "http://samtools.github.io/bcftools/",
                        "label": "Homepage"
                    },
                    {
                        "id": "https://github.com/samtools/bcftools",
                        "label": "Source code"
                    },
                    {
                        "id": "https://github.com/samtools/bcftools/wiki",
                        "label": "Wiki"
                    },
                    {
                        "id": "https://github.com/samtools/bcftools/archive/1.9.zip",
                        "label": "Download"
                    }
                ],
                "sbg:cmdPreview": "bcftools index  -f -t input_file.vcf.gz && bcftools consensus --fasta-ref /path/to/reference.ext --output input_file.fa  input_file.vcf.gz",
                "sbg:revisionsInfo": [
                    {
                        "sbg:revision": 0,
                        "sbg:modifiedBy": "admin",
                        "sbg:modifiedOn": 1538758819,
                        "sbg:revisionNotes": null
                    },
                    {
                        "sbg:revision": 1,
                        "sbg:modifiedBy": "admin",
                        "sbg:modifiedOn": 1538758819,
                        "sbg:revisionNotes": "Initial"
                    },
                    {
                        "sbg:revision": 2,
                        "sbg:modifiedBy": "admin",
                        "sbg:modifiedOn": 1538758819,
                        "sbg:revisionNotes": "Description"
                    },
                    {
                        "sbg:revision": 3,
                        "sbg:modifiedBy": "admin",
                        "sbg:modifiedOn": 1538758819,
                        "sbg:revisionNotes": "Description"
                    },
                    {
                        "sbg:revision": 4,
                        "sbg:modifiedBy": "admin",
                        "sbg:modifiedOn": 1561459470,
                        "sbg:revisionNotes": "Updated default CPU and memory parameter"
                    },
                    {
                        "sbg:revision": 5,
                        "sbg:modifiedBy": "admin",
                        "sbg:modifiedOn": 1649155693,
                        "sbg:revisionNotes": "update categories"
                    },
                    {
                        "sbg:revision": 6,
                        "sbg:modifiedBy": "admin",
                        "sbg:modifiedOn": 1654241373,
                        "sbg:revisionNotes": "Tool fix for compressed and uncompressed inputs."
                    }
                ],
                "sbg:toolkit": "bcftools",
                "abg:revisionNotes": "Initial version",
                "sbg:image_url": null,
                "sbg:job": {
                    "allocatedResources": {
                        "cpu": 8,
                        "mem": 10000
                    },
                    "inputs": {
                        "exclude_expression": "",
                        "memory": null,
                        "include_expression": "",
                        "missing": "",
                        "haplotype": null,
                        "reference": {
                            "class": "File",
                            "size": 0,
                            "secondaryFiles": [],
                            "path": "/path/to/reference.ext"
                        },
                        "output_name": "",
                        "sample": "",
                        "input_file": {
                            "class": "File",
                            "size": 0,
                            "secondaryFiles": [
                                {
                                    "path": ".tbi"
                                }
                            ],
                            "path": "/path/to/input_file.vcf.gz"
                        },
                        "iupac": false,
                        "cpu": null
                    }
                },
                "sbg:license": "MIT License",
                "sbg:expand_workflow": false,
                "sbg:appVersion": [
                    "sbg:draft-2"
                ],
                "sbg:id": "admin/sbg-public-data/bcftools-consensus/6",
                "sbg:revision": 6,
                "sbg:revisionNotes": "Tool fix for compressed and uncompressed inputs.",
                "sbg:modifiedOn": 1654241373,
                "sbg:modifiedBy": "admin",
                "sbg:createdOn": 1538758819,
                "sbg:createdBy": "admin",
                "sbg:project": "admin/sbg-public-data",
                "sbg:sbgMaintained": false,
                "sbg:validationErrors": [],
                "sbg:contributors": [
                    "admin"
                ],
                "sbg:latestRevision": 6,
                "sbg:publisher": "sbg",
                "sbg:content_hash": "ae03e79e31dbcf94a369711081f996251683283e3cc651c03abf20b0ad6becda7",
                "sbg:workflowLanguage": "CWL"
            },
            "label": "BCFtools Consensus",
            "sbg:x": 556.9197387695312,
            "sbg:y": 21.124088287353516
        },
        {
            "id": "bwa_mem_bundle_0_7_17_cwl_1_2",
            "in": [
                {
                    "id": "reference_index_tar",
                    "source": "bwa_index_0_7_17_cwl_1_2/indexed_reference"
                },
                {
                    "id": "input_reads",
                    "source": [
                        "input_file"
                    ]
                },
                {
                    "id": "deduplication",
                    "default": "RemoveDuplicates"
                }
            ],
            "out": [
                {
                    "id": "aligned_reads"
                }
            ],
            "run": {
                "class": "CommandLineTool",
                "cwlVersion": "v1.2",
                "$namespaces": {
                    "sbg": "https://sevenbridges.com"
                },
                "id": "admin/sbg-public-data/bwa-mem-bundle-0-7-17-cwl-1-2/1",
                "baseCommand": [],
                "inputs": [
                    {
                        "sbg:category": "Input files",
                        "id": "reference_index_tar",
                        "type": "File",
                        "label": "Reference Index TAR",
                        "doc": "Reference fasta file with its BWA index files packed in a TAR archive.",
                        "sbg:fileTypes": "TAR"
                    },
                    {
                        "sbg:category": "Input files",
                        "id": "input_reads",
                        "type": "File[]",
                        "inputBinding": {
                            "shellQuote": false,
                            "position": 105,
                            "valueFrom": "${\n    /// Set input reads in the correct order depending of the paired end from metadata\n\n    // Set output file name\n    function flatten(files){\n        var a = [];\n        for(var i=0;i<files.length;i++){\n            if(files[i]){\n                if(files[i].constructor == Array) a = a.concat(flatten(files[i]));\n                else a = a.concat(files[i]);}}\n        var b = a.filter(function (el) {return el != null;})\n        return b;}\n    var files1 = [].concat(inputs.input_reads);\n    var in_reads=flatten(files1);\n\n    // Read metadata for input reads\n    var read_metadata = in_reads[0].metadata;\n    if (!read_metadata) read_metadata = [];\n\n    var order = 0; // Consider this as normal order given at input: pe1 pe2\n\n    // Check if paired end 1 corresponds to the first given read\n    if (read_metadata == []) order = 0;\n    else if ('paired_end' in read_metadata) {\n        var pe1 = read_metadata.paired_end;\n        if (pe1 != 1) order = 1; // change order\n    }\n\n    // Return reads in the correct order\n    if (in_reads.length == 1) return in_reads[0].path; // Only one read present\n    else if (in_reads.length == 2) {\n        if (order == 0) return in_reads[0].path + ' ' + in_reads[1].path;\n        else return in_reads[1].path + ' ' + in_reads[0].path;\n    }\n}"
                        },
                        "label": "Input reads",
                        "doc": "Input sequence reads.",
                        "sbg:fileTypes": "FASTQ, FASTQ.GZ, FQ, FQ.GZ"
                    },
                    {
                        "sbg:category": "Input files",
                        "id": "fasta_index",
                        "type": "File?",
                        "inputBinding": {
                            "shellQuote": false,
                            "position": 4,
                            "valueFrom": "${\n    return \"\";\n}"
                        },
                        "label": "Fasta Index file for CRAM output",
                        "doc": "Fasta index file is required for CRAM output when no PCR Deduplication is selected.",
                        "sbg:fileTypes": "FAI"
                    },
                    {
                        "sbg:category": "Execution options",
                        "sbg:toolDefaultValue": "8",
                        "id": "threads",
                        "type": "int?",
                        "label": "Threads",
                        "doc": "The number of threads for BWA and Biobambam2 sort processes (both will use the given number)."
                    },
                    {
                        "sbg:category": "BWA Algorithm options",
                        "sbg:toolDefaultValue": "19",
                        "id": "minimum_seed_length",
                        "type": "int?",
                        "inputBinding": {
                            "prefix": "-k",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Minimum seed length",
                        "doc": "Minimum seed length for BWA MEM."
                    },
                    {
                        "sbg:category": "BWA Algorithm options",
                        "sbg:toolDefaultValue": "100",
                        "id": "dropoff",
                        "type": "int?",
                        "inputBinding": {
                            "prefix": "-d",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Dropoff",
                        "doc": "Off-diagonal X-dropoff."
                    },
                    {
                        "sbg:category": "BWA Algorithm options",
                        "sbg:toolDefaultValue": "1.5",
                        "id": "select_seeds",
                        "type": "float?",
                        "inputBinding": {
                            "prefix": "-r",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Select seeds",
                        "doc": "Look for internal seeds inside a seed longer than {-k} * FLOAT."
                    },
                    {
                        "sbg:category": "BWA Algorithm options",
                        "sbg:toolDefaultValue": "20",
                        "id": "seed_occurrence_for_the_3rd_round",
                        "type": "int?",
                        "inputBinding": {
                            "prefix": "-y",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Seed occurrence",
                        "doc": "Seed occurrence for the 3rd round seeding."
                    },
                    {
                        "sbg:category": "BWA Algorithm options",
                        "sbg:toolDefaultValue": "500",
                        "id": "skip_seeds",
                        "type": "int?",
                        "inputBinding": {
                            "prefix": "-c",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Skip seeds",
                        "doc": "Skip seeds with more than a given number (INT) of occurrences."
                    },
                    {
                        "sbg:category": "BWA Algorithm options",
                        "sbg:toolDefaultValue": "0.50",
                        "id": "drop_chains_fraction",
                        "type": "float?",
                        "inputBinding": {
                            "prefix": "-D",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Drop chains fraction",
                        "doc": "Drop chains shorter than a given fraction (FLOAT) of the longest overlapping chain."
                    },
                    {
                        "sbg:category": "BWA Algorithm options",
                        "sbg:toolDefaultValue": "0",
                        "id": "discard_chain_length",
                        "type": "int?",
                        "inputBinding": {
                            "prefix": "-W",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Discard chain length",
                        "doc": "Discard a chain if seeded bases are shorter than a given number (INT)."
                    },
                    {
                        "sbg:category": "BWA Algorithm options",
                        "sbg:toolDefaultValue": "50",
                        "id": "mate_rescue_rounds",
                        "type": "string?",
                        "inputBinding": {
                            "prefix": "-m",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Mate rescue rounds",
                        "doc": "Perform at the most a given number (INT) of rounds of mate rescues for each read."
                    },
                    {
                        "sbg:category": "BWA Algorithm options",
                        "id": "skip_mate_rescue",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "-S",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Skip mate rescue",
                        "doc": "Skip mate rescue."
                    },
                    {
                        "sbg:category": "BWA Algorithm options",
                        "id": "skip_pairing",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "-P",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Skip pairing",
                        "doc": "Skip pairing; mate rescue is performed unless -S also in use."
                    },
                    {
                        "sbg:category": "BWA Algorithm options",
                        "id": "discard_exact_matches",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "-e",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Discard exact matches",
                        "doc": "Discard full-length exact matches."
                    },
                    {
                        "sbg:category": "BWA Scoring options",
                        "sbg:toolDefaultValue": "1",
                        "id": "score_for_a_sequence_match",
                        "type": "int?",
                        "inputBinding": {
                            "prefix": "-A",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Score for a sequence match",
                        "doc": "Score for a sequence match, which scales options -TdBOELU unless overridden."
                    },
                    {
                        "sbg:category": "BWA Scoring options",
                        "sbg:toolDefaultValue": "4",
                        "id": "mismatch_penalty",
                        "type": "int?",
                        "inputBinding": {
                            "prefix": "-B",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Mismatch penalty",
                        "doc": "Penalty for a mismatch."
                    },
                    {
                        "sbg:category": "BWA Input/output options",
                        "id": "smart_pairing_in_input_fastq",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "-p",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Smart pairing",
                        "doc": "Smart pairing in input FASTQ file (ignoring in2.fq)."
                    },
                    {
                        "sbg:category": "BWA Read Group Options",
                        "sbg:toolDefaultValue": "Constructed from per-attribute parameters or inferred from metadata.",
                        "id": "read_group_header",
                        "type": "string?",
                        "label": "Read group header",
                        "doc": "Read group header line such as '@RG\\tID:foo\\tSM:bar'.  This value takes precedence over per-attribute parameters."
                    },
                    {
                        "sbg:category": "BWA Input/output options",
                        "id": "insert_string_to_header",
                        "type": "string?",
                        "inputBinding": {
                            "prefix": "-H",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Insert string to header",
                        "doc": "Insert STR to output header if it starts with \"@\"."
                    },
                    {
                        "sbg:category": "BWA Input/output options",
                        "id": "ignore_alt_file",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "-j",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Ignore ALT file",
                        "doc": "Treat ALT contigs as part of the primary assembly (i.e. ignore <idxbase>.alt file)."
                    },
                    {
                        "sbg:category": "BWA Input/output options",
                        "sbg:toolDefaultValue": "3",
                        "id": "verbose_level",
                        "type": [
                            "null",
                            {
                                "type": "enum",
                                "symbols": [
                                    "1",
                                    "2",
                                    "3",
                                    "4"
                                ],
                                "name": "verbose_level"
                            }
                        ],
                        "inputBinding": {
                            "prefix": "-v",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Verbose level",
                        "doc": "Select verbose level: 1=error, 2=warning, 3=message, 4+=debugging."
                    },
                    {
                        "sbg:category": "BWA Input/output options",
                        "sbg:toolDefaultValue": "30",
                        "id": "minimum_output_score",
                        "type": "int?",
                        "inputBinding": {
                            "prefix": "-T",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Minimum alignment score",
                        "doc": "Minimum alignment score for a read to be outputted in SAM/BAM/CRAM."
                    },
                    {
                        "sbg:category": "BWA Input/output options",
                        "id": "output_alignments",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "-a",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Output alignments",
                        "doc": "Output all alignments for SE or unpaired PE."
                    },
                    {
                        "sbg:category": "BWA Input/output options",
                        "id": "append_comment",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "-C",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Append comment",
                        "doc": "Append FASTA/FASTQ comment to the output file."
                    },
                    {
                        "sbg:category": "BWA Input/output options",
                        "id": "output_header",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "-V",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Output header",
                        "doc": "Output the reference FASTA header in the XR tag."
                    },
                    {
                        "sbg:category": "BWA Input/output options",
                        "id": "use_soft_clipping",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "-Y",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Use soft clipping",
                        "doc": "Use soft clipping for supplementary alignments."
                    },
                    {
                        "sbg:category": "BWA Input/output options",
                        "id": "mark_shorter",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "-M",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Mark shorter",
                        "doc": "Mark shorter split hits as secondary."
                    },
                    {
                        "sbg:category": "BWA Scoring options",
                        "sbg:toolDefaultValue": "[6,6]",
                        "id": "gap_open_penalties",
                        "type": "int[]?",
                        "inputBinding": {
                            "prefix": "-O",
                            "separate": false,
                            "itemSeparator": ",",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Gap open penalties",
                        "doc": "Gap open penalties for deletions and insertions. \nThis array can't have more than two values."
                    },
                    {
                        "sbg:category": "BWA Scoring options",
                        "sbg:toolDefaultValue": "[1,1]",
                        "id": "gap_extension_penalties",
                        "type": "int[]?",
                        "inputBinding": {
                            "prefix": "-E",
                            "separate": false,
                            "itemSeparator": ",",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Gap extension",
                        "doc": "Gap extension penalty; a gap of size k cost '{-O} + {-E}*k'. \nThis array can't have more than two values."
                    },
                    {
                        "sbg:category": "BWA Scoring options",
                        "sbg:toolDefaultValue": "[5,5]",
                        "id": "clipping_penalty",
                        "type": "int[]?",
                        "inputBinding": {
                            "prefix": "-L",
                            "separate": false,
                            "itemSeparator": ",",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Clipping penalty",
                        "doc": "Penalty for 5'- and 3'-end clipping."
                    },
                    {
                        "sbg:category": "BWA Scoring options",
                        "sbg:toolDefaultValue": "17",
                        "id": "unpaired_read_penalty",
                        "type": "int?",
                        "inputBinding": {
                            "prefix": "-U",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Unpaired read penalty",
                        "doc": "Penalty for an unpaired read pair."
                    },
                    {
                        "sbg:category": "BWA Scoring options",
                        "id": "read_type",
                        "type": [
                            "null",
                            {
                                "type": "enum",
                                "symbols": [
                                    "pacbio",
                                    "ont2d",
                                    "intractg"
                                ],
                                "name": "read_type"
                            }
                        ],
                        "inputBinding": {
                            "prefix": "-x",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Sequencing technology-specific settings",
                        "doc": "Sequencing technology-specific settings; Setting -x changes multiple parameters unless overridden. \npacbio: -k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0  (PacBio reads to ref). \nont2d: -k14 -W20 -r10 -A1 -B1 -O1 -E1 -L0  (Oxford Nanopore 2D-reads to ref).\nintractg: -B9 -O16 -L5  (intra-species contigs to ref)."
                    },
                    {
                        "sbg:category": "BWA Input/output options",
                        "sbg:toolDefaultValue": "[5, 200]",
                        "id": "output_in_xa",
                        "type": "int[]?",
                        "inputBinding": {
                            "prefix": "-h",
                            "separate": false,
                            "itemSeparator": ",",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Output in XA",
                        "doc": "If there are < number (INT) of hits with a score >80% of the max score, output all in XA. \nThis array should have no more than two values."
                    },
                    {
                        "sbg:category": "BWA Input/output options",
                        "id": "speficy_distribution_parameters",
                        "type": [
                            "null",
                            {
                                "type": "array",
                                "items": "float",
                                "inputBinding": {
                                    "prefix": "-I",
                                    "separate": false
                                }
                            }
                        ],
                        "inputBinding": {
                            "shellQuote": false,
                            "position": 4,
                            "valueFrom": "${\n    var out = \"\"\n    for (var i = 0; i < [].concat(self).length; i++ ){\n        out += \" -I\" + [].concat(self)[i]\n    }    \n    return out\n}"
                        },
                        "label": "Specify distribution parameters",
                        "doc": "Specify the mean, standard deviation (10% of the mean if absent), max (4 sigma from the mean if absent), and min of the insert size distribution. \nFR orientation only. \nThis array can have maximum of four values, where the first two should be specified as FLOAT and the last two as INT."
                    },
                    {
                        "sbg:category": "BWA Algorithm options",
                        "sbg:toolDefaultValue": "100",
                        "id": "band_width",
                        "type": "int?",
                        "inputBinding": {
                            "prefix": "-w",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Band width",
                        "doc": "Band width for banded alignment."
                    },
                    {
                        "sbg:category": "BWA Read Group Options",
                        "sbg:toolDefaultValue": "Inferred from metadata or set to default (\"Illumina\")",
                        "id": "rg_platform",
                        "type": [
                            "null",
                            {
                                "type": "enum",
                                "symbols": [
                                    "454",
                                    "Helicos",
                                    "Illumina",
                                    "Solid",
                                    "IonTorrent"
                                ],
                                "name": "rg_platform"
                            }
                        ],
                        "label": "Platform",
                        "doc": "Specify the version of the technology that was used for sequencing, which will be placed in RG line."
                    },
                    {
                        "sbg:category": "BWA Read Group Options",
                        "sbg:toolDefaultValue": "Inferred from metadata, input name or set to default",
                        "id": "rg_sample_id",
                        "type": "string?",
                        "label": "Sample ID",
                        "doc": "Specify the sample ID for RG line - A human readable identifier for a sample or specimen, which could contain some metadata information. A sample or specimen is material taken from a biological entity for testing, diagnosis, propagation, treatment, or research purposes, including but not limited to tissues, body fluids, cells, organs, embryos, body excretory products, etc."
                    },
                    {
                        "sbg:category": "BWA Read Group Options",
                        "sbg:toolDefaultValue": "Inferred from metadata",
                        "id": "rg_library_id",
                        "type": "string?",
                        "label": "Library ID",
                        "doc": "Specify the identifier for the sequencing library preparation, which will be placed in RG line."
                    },
                    {
                        "sbg:category": "BWA Read Group Options",
                        "sbg:toolDefaultValue": "Inferred from metadata",
                        "id": "rg_platform_unit_id",
                        "type": "string?",
                        "label": "Platform unit ID",
                        "doc": "Specify the platform unit (lane/slide) for RG line - An identifier for lanes (Illumina), or for slides (SOLiD) in the case that a library was split and ran over multiple lanes on the flow cell or slides."
                    },
                    {
                        "sbg:category": "BWA Read Group Options",
                        "id": "rg_data_submitting_center",
                        "type": "string?",
                        "label": "Data submitting center",
                        "doc": "Specify the data submitting center for RG line."
                    },
                    {
                        "sbg:category": "BWA Read Group Options",
                        "id": "rg_median_fragment_length",
                        "type": "string?",
                        "label": "Median fragment length",
                        "doc": "Specify the median fragment length for RG line."
                    },
                    {
                        "sbg:category": "Execution options",
                        "sbg:toolDefaultValue": "BAM",
                        "id": "output_format",
                        "type": [
                            "null",
                            {
                                "type": "enum",
                                "symbols": [
                                    "SAM",
                                    "BAM",
                                    "CRAM",
                                    "Queryname Sorted BAM",
                                    "Queryname Sorted SAM"
                                ],
                                "name": "output_format"
                            }
                        ],
                        "label": "Output format",
                        "doc": "Coordinate sorted BAM file (option BAM) is the default output."
                    },
                    {
                        "sbg:category": "Execution options",
                        "id": "sort_memory",
                        "type": "int?",
                        "label": "Memory for BAM sorting [GB]",
                        "doc": "Amount of RAM [Gb] to give to the sorting algorithm (if not provided will be set to one-third of the total memory)."
                    },
                    {
                        "sbg:category": "Biobambam2 parameters",
                        "sbg:toolDefaultValue": "None",
                        "id": "deduplication",
                        "type": [
                            "null",
                            {
                                "type": "enum",
                                "symbols": [
                                    "None",
                                    "MarkDuplicates",
                                    "RemoveDuplicates"
                                ],
                                "name": "deduplication"
                            }
                        ],
                        "label": "PCR duplicate detection",
                        "doc": "Use Biobambam2 for finding duplicates on sequence reads."
                    },
                    {
                        "sbg:category": "Execution options",
                        "sbg:toolDefaultValue": "15",
                        "id": "total_memory",
                        "type": "int?",
                        "label": "Total memory [GB]",
                        "doc": "Total memory to be used by the tool in GB. It's the sum of BWA and BIOBAMBAM2 processes. For FASTQ files of a total size less than 10GB, we suggest using the default setting of 15GB, for larger files, we suggest using 58GB of memory (and 32CPU cores)."
                    },
                    {
                        "sbg:category": "Configuration",
                        "id": "output_name",
                        "type": "string?",
                        "label": "Output alignements file name",
                        "doc": "Name for the output alignments (SAM, BAM, or CRAM) file."
                    },
                    {
                        "sbg:category": "Configuration",
                        "sbg:toolDefaultValue": "1",
                        "id": "reserved_threads",
                        "type": "int?",
                        "label": "Reserved number of threads on the instance",
                        "doc": "Reserved number of threads on the instance used by scheduler."
                    },
                    {
                        "sbg:category": "Configuration",
                        "sbg:toolDefaultValue": "1",
                        "id": "rg_id",
                        "type": "string?",
                        "label": "Read group ID",
                        "doc": "Set read group ID."
                    },
                    {
                        "sbg:category": "Execution options",
                        "sbg:toolDefaultValue": "False",
                        "id": "wgs_hg38_mode_threads",
                        "type": "int?",
                        "label": "Optimize threads for HG38",
                        "doc": "Lower the number of threads if HG38 reference genome is used."
                    },
                    {
                        "id": "split_alignment_primary",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "-5",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Split alignment - smallest coordinate as primary",
                        "doc": "For split alignment, take the alignment with the smallest coordinate as primary."
                    },
                    {
                        "id": "mapQ_of_suplementary",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "-q",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Don't modify mapQ",
                        "doc": "Don't modify mapQ of supplementary alignments."
                    },
                    {
                        "id": "num_input_bases_in_each_batch",
                        "type": "int?",
                        "inputBinding": {
                            "prefix": "-K",
                            "shellQuote": false,
                            "position": 4
                        },
                        "label": "Number of input bases to process",
                        "doc": "Process a given number (INT) of input bases in each batch regardless of nThreads (for reproducibility)."
                    }
                ],
                "outputs": [
                    {
                        "id": "aligned_reads",
                        "doc": "Output SAM/BAM/CRAM file containing aligned reads.",
                        "label": "Output alignments",
                        "type": "File?",
                        "outputBinding": {
                            "glob": "${ \n    return [\"*.sam\", \"*.bam\", \"*.cram\"] \n}",
                            "outputEval": "${  \n    /// Set metadata from input parameters, metadata or default value\n\n    function flatten(files){\n        var a = []\n        for(var i=0;i<files.length;i++){\n            if(files[i]){\n                if(files[i].constructor == Array) a = a.concat(flatten(files[i]));\n                else a = a.concat(files[i]);}}\n        var b = a.filter(function (el) {return el != null});\n        return b;\n    }\n    function sharedStart(array){\n        var A= array.concat().sort(), \n        a1= A[0], a2= A[A.length-1], L= a1.length, i= 0;\n        while(i<L && a1.charAt(i)=== a2.charAt(i)) i++;\n        return a1.substring(0, i);\n    }\n    /// Key-setting functions\n    // Reference genome \n    var add_metadata_key_reference_genome = function(self, inputs) {\n        var reference_file = inputs.reference_index_tar.basename;\n        var ref_list = reference_file.split('.');\n        var  a = '';\n        var ref_gen='';\n        if (reference_file.includes(\".bwa-0.7.17-index-archive.tar\")){\n            a = ref_list.pop();\n            a = ref_list.pop();\n            a = ref_list.pop();\n            a = ref_list.pop(); // strip '.bwa-0.7.17-index-archive.tar'\n        }else{\n            a = ref_list.pop();\n            a = ref_list.pop(); // strip '.tar'\n        } \n        if (read_metadata.reference_genome) {\n            ref_gen = read_metadata.reference_genome;\n        }else {ref_gen=ref_list.join('.');}\n        return ref_gen\n    };\n    // Platform \n    var add_metadata_key_platform = function(self, inputs) {\n        /// Set platform from input parameters/input metadata/default value\n        var platform = '';\n        var pl = '';\n        // Find PL from header\n        if (inputs.read_group_header){\n            var header = inputs.read_group_header;\n            header = header.split(\"'\").join(\"\") //remove single quotes\n            var a = header.split('\\\\t');\n            for (var i = 0; i < a.length; i++){ //find PL field\n                if (a[i].includes(\"PL:\")) pl= a[i];\n                else;\n            }}\n        else;\n        \n        if (pl) platform = pl.split(':')[1];\n        else if (inputs.rg_platform) platform = inputs.rg_platform;\n        else if (read_metadata.platform) platform = read_metadata.platform;\n        else platform = 'Illumina';\n        \n        return platform\n    };\n    // Sample ID \n    var add_metadata_key_sample_id = function(self, inputs) {\n        /// Set sample ID from input parameters/input metadata/default value from input reads file names\n        var sample_id = '';\n        var sm = '';\n        // Find SM from header\n        if (inputs.read_group_header){\n            var header = inputs.read_group_header;\n            header = header.split(\"'\").join(\"\") //remove single quotes\n            var a = header.split('\\\\t');\n            for (var i = 0; i < a.length; i++){ //find SM field\n                if (a[i].includes(\"SM:\")) var sm= a[i];\n                else;\n            }}\n        else;\n        \n        if (sm) sample_id = sm.split(':')[1];\n        else if (inputs.rg_sample_id) sample_id = inputs.rg_sample_id;\n        else if (read_metadata.sample_id) sample_id = read_metadata.sample_id;\n        else {\n            var read_names = [];\n            var files1 = [].concat(inputs.input_reads);\n            var files=flatten(files1);\n            \n            for (var i=0;i<files.length;i++) {\n                var file_ext=files[i].nameext;\n                var file_base=files[i].basename;\n                \n                if (file_ext === '.gz' || file_ext === '.GZ')\n                    file_base = file_base.slice(0, -3);\n                    file_ext= '.'+ file_base.split('.').pop();\n                if (file_ext === '.fq' || file_ext === '.FQ')\n                    file_base = file_base.slice(0, -3);\n                if (file_ext === '.fastq' || file_ext === '.FASTQ')\n                    file_base = file_base.slice(0, -6);\n                \n                read_names.push(file_base.replace(/pe1|pe2|pe\\.1|pe\\.2|pe\\_1|pe\\_2|\\_pe1|\\_pe2|\\_pe\\.1|\\_pe\\.2|\\_pe\\_1|\\_pe\\_2|\\.pe1|\\.pe2|\\.pe\\.1|\\.pe\\.2|\\.pe\\_1|\\.pe\\_2/,''));\n              }\n              ////strip out any trailing dashes/dots/underscores...\n              var unique_prefix = sharedStart(read_names).replace( /\\-$|\\_$|\\.$/, '');\n              var tmp_prefix = unique_prefix.replace( /^\\_|\\.pe$|\\.R$|\\_pe$|\\_R$/,'');\n              var final_prefix = tmp_prefix.replace( /^_\\d(\\d)?_/, '' );\n              \n              var fname=final_prefix;\n            sample_id = fname;\n        }\n        return sample_id\n    };\n    \n   \n    var files1 = [].concat(inputs.input_reads);\n    var files=flatten(files1);\n    var read_metadata = files[0].metadata;\n    if (!read_metadata) read_metadata = [];\n    \n    self = inheritMetadata(self, files);\n\n    for (var i = 0; i < self.length; i++) {\n        var out_metadata = {\n            'reference_genome': add_metadata_key_reference_genome(self[i], inputs),\n            'platform': add_metadata_key_platform(self[i], inputs),\n            'sample_id': add_metadata_key_sample_id(self[i], inputs)\n        };\n        self[i] = setMetadata(self[i], out_metadata);\n    }\n\n    return self;\n\n}"
                        },
                        "secondaryFiles": [
                            {
                                "pattern": ".bai",
                                "required": false
                            },
                            {
                                "pattern": "^.bai",
                                "required": false
                            },
                            {
                                "pattern": ".crai",
                                "required": false
                            },
                            {
                                "pattern": "^.crai",
                                "required": false
                            }
                        ],
                        "sbg:fileTypes": "SAM, BAM, CRAM"
                    }
                ],
                "doc": "BWA-MEM is an algorithm designed for aligning sequence reads onto a large reference genome. BWA-MEM is implemented as a component of BWA. The algorithm can automatically choose between performing end-to-end and local alignments. BWA-MEM is capable of outputting multiple alignments, and finding chimeric reads. It can be applied to a wide range of read lengths, from 70 bp to several megabases. \n\n*A list of **all inputs and parameters** with corresponding descriptions can be found at the bottom of this page.*\n\n***Please note that any cloud infrastructure costs resulting from app and pipeline executions, including the use of public apps, are the sole responsibility of you as a user. To avoid excessive costs, please read the app description carefully and set the app parameters and execution settings accordingly.***\n\n\n## Common Use Cases\nIn order to enable additional fast processing of aligned reads, the **Biobambam2 sortmadup** (2.0.87) tool is embedded together into the same package with BWA-MEM (0.7.17).\n\nIn order to enable additional fast processing of aligned reads, **Biobambam2** (2.0.87) is embedded together with the BWA 0.7.17 toolkit into **BWA-MEM Bundle**.  Two tools are used (**bamsort** and **bamsormadup**) to allow the selection of three output formats (SAM, BAM, or CRAM), different modes of sorting (queryname/coordinate sorting), and marking/removing duplicates that can arise during sample preparation e.g. library construction using PCR. This is done by setting the **Output format** and **PCR duplicate detection** parameters.\n- Additional notes:\n    - The default **Output format** is coordinate sorted BAM (option **BAM**).\n    - SAM and BAM options are query name sorted, while CRAM format is not advisable for data sorted by query name.\n    - Coordinate Sorted BAM file in all options and CRAM Coordinate sorted output with Marked Duplicates come with their accompanying index files. The generated index name will be the same as the output alignments file, with the BAM.BAI or CRAM.CRAI extension. However, when selecting the CRAM Coordinate sorted and CRAM Coordinate sorted output with Removed Duplicates, the generated files will not have the accompanying index files. This is a result of the usage of different Biobambam2 tools - **bamsort** does not have the ability to write CRAI files (only supports outputting BAI index files), while **bamsormadup** can write CRAI files.\n    - Passing data from BWA-MEM to Biobambam2 tools has been done through  Linux piping which saves processing times (up to an hour of the execution time for a whole-genome sample) of reading and writing of aligned reads into the hard drive. \n    - **BWA-MEM Bundle 0.7.17 CWL1.2** first needs to construct the FM-index  (Full-text index in Minute space) for the reference genome using the **BWA INDEX 0.7.17 CWL1.2** tool. The two BWA versions are compatible.\n\n### Changes Introduced by Seven Bridges\n\n- **Aligned SAM/BAM/CRAM** file will be prefixed using the **Output SAM/BAM/CRAM file name** parameter. In case **Output SAM/BAM/CRAM file name** is not provided, the output prefix will be the same as the **Sample ID** metadata field from the file if the **Sample ID** metadata field exists. Otherwise, the output prefix will be inferred from the **Input reads** file names.\n-  The **Platform** metadata field for the output alignments will be automatically set to \"Illumina\" unless it is present in **Input reads** metadata, or given through **Read group header** or **Platform** input parameters. This will prevent possible errors in downstream analysis using the GATK toolkit.\n- If the **Read group ID** parameter is not defined, by default it will be set to ‘1’. If the tool is scattered within a workflow it will assign the **Read Group ID** according to the order of the scattered folders. This ensures a unique **Read Group ID** when processing multi-read group input data from one sample.\n\n### Common Issues and Important Notes \n \n- For input reads FASTQ files of total size less than 10 GB we suggest using the default setting of 15GB for the **Total memory** parameter. For larger files, we suggest using 58 GB of memory and 32 CPU cores.\n- When the desired output is a CRAM file without deduplication of PCR duplicates, it is necessary to provide the FASTA Index file (FAI) as input.\n- Human reference genome version 38 comes with ALT contigs, a collection of diverged alleles present in some humans but not in others. Making effective use of these contigs will help reduce mapping artifacts. However, to facilitate mapping these ALT contigs to the primary assembly, GRC decided to add to each contig long flanking sequences almost identical to the primary assembly. As a result, a naive mapping against GRCh38+ALT will lead to many mapQ-zero mappings in these flanking regions. Please use post-processing steps to fix these alignments or implement [steps](https://sourceforge.net/p/bio-bwa/mailman/message/32845712/) described by the author of the BWA toolkit.  \n- Inputs **Read group header** and **Insert string to header** need to be given in the correct format - under single-quotes.\n- BWA-MEM is not a splice aware aligner, so it is not the appropriate tool for mapping RNAseq to the genome. For RNAseq reads **Bowtie2 Aligner** and **STAR** are recommended tools. \n- Input paired reads need to have identical read names - if not, the tool will throw a ``[mem_sam_pe] paired reads have different names`` error.\n\n### Limitations\n\n- This app was tested only on human data (WES and WGS), aligning to both hg38 and hg19 reference genome assemblies.\n\n\n### Performance Benchmarking\n\nBelow is a table describing runtimes and task costs on on-demand instances for a set of samples with different file sizes:\n\n| Input reads       | Size [GB] | Output format | Instance (AWS)           | Duration  | Cost   | Threads |\n|-------------------|-----------|---------------|--------------------------|-----------|--------|---------|\n| HG001-NA12878-30x | 2 x 23.8  | SAM           | c5.9xlarge (36CPU, 72GB) | 5h 12min  | $7.82  | 36      |\n| HG001-NA12878-30x | 2 x 23.8  | BAM           | c5.9xlarge (36CPU, 72GB) | 5h 16min  | $8.06  | 36      |\n| HG002-NA24385-50x | 2 x 66.4  | SAM           | c5.9xlarge (36CPU, 72GB) | 8h 33min  | $13.08 | 36      |\n\n\n*Cost can be significantly reduced by using **spot instances**. Visit the [Knowledge Center](https://docs.sevenbridges.com/docs/about-spot-instances) for more details.*\n\n### Portability\n\n**BWA MEM Bundle** was tested with cwltool version 3.0.20201203173111. The '-input_reads' and  '-reference_index_tar' inputs were provided in the job.yaml/job.json file and used for testing.",
                "label": "BWA MEM Bundle",
                "arguments": [
                    {
                        "prefix": "",
                        "shellQuote": false,
                        "position": -1,
                        "valueFrom": "${\n    /// Check number of input FASTQ files ///\n    \n    function flatten(files){\n    var a = []\n    for(var i=0;i<files.length;i++){\n        if(files[i]){\n            if(files[i].constructor == Array) a = a.concat(flatten(files[i]));\n            else a = a.concat(files[i])}}\n        var b = a.filter(function (el) {return el != null})\n        return b\n    }\n    \n    var files1 = [].concat(inputs.input_reads);\n    var in_reads=flatten(files1);\n    \n    if ( in_reads.length > 2 ) return 'ERROR: Number of input FASTQ files needs to be one (if single-end/interleaved file) or two (if paired-end files)';\n    else return '';\n}"
                    },
                    {
                        "prefix": "",
                        "shellQuote": false,
                        "position": 0,
                        "valueFrom": "${\n    var cmd = \"/bin/bash -c \\\"\";\n    return cmd + \" export REF_CACHE=${PWD} && \";\n}"
                    },
                    {
                        "prefix": "",
                        "shellQuote": false,
                        "position": 1,
                        "valueFrom": "${\n    /// Unpack Reference TAR archive ///\n    \n    var in_index=[].concat(inputs.reference_index_tar)[0];\n    var reference_file = in_index.basename;\n    return 'tar -tvf ' + reference_file + ' 1>&2 && tar -xf ' + reference_file + ' && ';\n    \n}"
                    },
                    {
                        "prefix": "",
                        "shellQuote": false,
                        "position": 2,
                        "valueFrom": "bwa mem"
                    },
                    {
                        "prefix": "",
                        "shellQuote": false,
                        "position": 5,
                        "valueFrom": "${\n    /// Set RG header ///\n\n    function add_param(key, val) {\n        if (!val) return;\n        param_list.push(key + ':' + val);}\n        \n    function flatten(files){\n        var a = [];\n        for(var i=0;i<files.length;i++){\n            if(files[i]){\n                if(files[i].constructor == Array) a = a.concat(flatten(files[i]));\n                else a = a.concat(files[i]);}}\n        var b = a.filter(function (el) {return el != null;});\n        return b;}\n        \n    function sharedStart(array){\n        var A= array.concat().sort(), \n        a1= A[0], a2= A[A.length-1], L= a1.length, i= 0;\n        while(i<L && a1.charAt(i)=== a2.charAt(i)) i++;\n        return a1.substring(0, i);}\n\n    \n    /// If it exists - return input read group header from input parameter\n    if (inputs.read_group_header) return '-R ' + inputs.read_group_header;\n\n    // Flatten input reads\n    var in_reads1 = [].concat(inputs.input_reads);\n    var in_reads = flatten(in_reads1)\n    var input_1=in_reads[0];\n\n    var param_list = [];\n    //Read metadata for input reads\n    var read_metadata = input_1.metadata;\n    if (!read_metadata) read_metadata = [];\n\n    // Set CN\n    if (inputs.rg_data_submitting_center) add_param('CN', inputs.rg_data_submitting_center);\n    else if ('data_submitting_center' in read_metadata) add_param('CN', read_metadata.data_submitting_center);\n    else;\n\n    // Set LB\n    if (inputs.rg_library_id) add_param('LB', inputs.rg_library_id);\n    else if ('library_id' in read_metadata) add_param('LB', read_metadata.library_id);\n    else;\n\n    // Set PI\n    if (inputs.rg_median_fragment_length) add_param('PI', inputs.rg_median_fragment_length);\n    else;\n\n    // Set PL (default Illumina)\n    var rg_platform = '';\n    if (inputs.rg_platform) add_param('PL', inputs.rg_platform);\n    else if ('platform' in read_metadata) {\n        if (read_metadata.platform == 'HiSeq X Ten') rg_platform = 'Illumina';\n        else rg_platform = read_metadata.platform;\n        add_param('PL', rg_platform);}\n    else add_param('PL', 'Illumina');\n\n    // Set PU\n    if (inputs.rg_platform_unit_id) add_param('PU', inputs.rg_platform_unit_id);\n    else if ('platform_unit_id' in read_metadata) add_param('PU', read_metadata.platform_unit_id);\n    else;\n    \n    // Set RG_ID\n    var folder = input_1.path.split('/').slice(-2,-1).toString();\n    var suffix = \"_s\";\n    var rg='';\n    //if RG ID set by user \n    if (inputs.rg_id) add_param('ID', inputs.rg_id);\n    //If BWA or sub wf containing is scattered\n    else if(folder.indexOf(suffix) !==-1){\n        if(folder.indexOf(suffix, folder.length - suffix.length) !== -1){\n            for(var i=0; i<folder.length; i++){\n                if(folder[i]=='_' && folder[i+1]=='s' && folder[i+2]=='_') rg+=folder[i-1];\n                else rg+='';\n                \n            }\n            rg+=folder.split(\"_\").slice(-2)[0];\n        }\n        else{\n            for(var i=0; i<folder.length; i++){\n                if(folder[i]=='_' && folder[i+1]=='s' && folder[i+2]=='_') rg+=folder[i-1];\n                else rg+='';\n                \n            }\n            \n        }\n        //remove all non-int chars \n        rg=rg.replace(/\\D/g,'')\n        if (rg) add_param('ID', rg);\n        else add_param('ID', 1);\n        \n    }// If no other conditions add ID 1\n    else  add_param('ID', 1);\n\n    // Set SM from input/metadata/filename\n    if (inputs.rg_sample_id) add_param('SM', inputs.rg_sample_id);\n    else if ('sample_id' in read_metadata) add_param('SM', read_metadata.sample_id);\n    else {\n        var read_names = [];\n        for (var i=0;i<in_reads.length;i++) {\n            var file_ext=in_reads[i].nameext;\n            var file_base=in_reads[i].basename;\n            \n            if (file_ext === '.gz' || file_ext === '.GZ')\n                file_base = file_base.slice(0, -3);\n                file_ext= '.'+ file_base.split('.').pop();\n            if (file_ext === '.fq' || file_ext === '.FQ')\n                file_base = file_base.slice(0, -3);\n            if (file_ext === '.fastq' || file_ext === '.FASTQ')\n                file_base = file_base.slice(0, -6);\n            \n            read_names.push(file_base.replace(/pe1|pe2|pe\\.1|pe\\.2|pe\\_1|pe\\_2|\\_pe1|\\_pe2|\\_pe\\.1|\\_pe\\.2|\\_pe\\_1|\\_pe\\_2|\\.pe1|\\.pe2|\\.pe\\.1|\\.pe\\.2|\\.pe\\_1|\\.pe\\_2/,''));}\n          \n        ////strip out any trailing dashes/dots/underscores...\n        var unique_prefix = sharedStart(read_names).replace( /\\-$|\\_$|\\.$/, '');\n        var tmp_prefix = unique_prefix.replace( /^\\_|\\.pe$|\\.R$|\\_pe$|\\_R$/,'');\n        var final_prefix = tmp_prefix.replace( /^_\\d(\\d)?_/, '' );\n      \n        var sample_id=final_prefix;\n        add_param('SM', sample_id);\n    };\n    \n    \n    // Create RG header\n    return \"-R '@RG\\\\t\" + param_list.join('\\\\t') + \"'\";\n\n}"
                    },
                    {
                        "prefix": "-t",
                        "shellQuote": false,
                        "position": 6,
                        "valueFrom": "${\n    /// Set BWA2 threads ///\n\n    var  MAX_THREADS = 36;\n    var  suggested_threads = 8;\n    var threads  = 0;\n  \n    if (inputs.threads) threads = inputs.threads;\n    else if (inputs.wgs_hg38_mode_threads) {\n        var ref_name = inputs.reference_index_tar.basename;\n        if (ref_name.search('38') >= 0) threads = inputs.wgs_hg38_mode_threads;\n        else threads = MAX_THREADS;\n    } else threads = suggested_threads;\n    \n    return threads;\n}"
                    },
                    {
                        "prefix": "",
                        "shellQuote": false,
                        "position": 14,
                        "valueFrom": "${\n    /// Extract common prefix for Index files ///\n    \n    var reference_tar = [].concat(inputs.reference_index_tar)[0];\n    \n    var prefix = \"$(tar -tf \" + reference_tar.basename + \" --wildcards '*.bwt' | rev | cut -c 5- | rev)\";\n    return prefix;\n\n}"
                    },
                    {
                        "prefix": "",
                        "separate": false,
                        "shellQuote": false,
                        "position": 116,
                        "valueFrom": "${\n    ///  BIOBAMBAM2  ///\n      \n     // Get shared start and flatten input reads\n    function sharedStart(array){\n        var A= array.concat().sort(), \n        a1= A[0], a2= A[A.length-1], L= a1.length, i= 0;\n        while(i<L && a1.charAt(i)=== a2.charAt(i)) i++;\n        return a1.substring(0, i);\n    }\n    function flatten(files){\n        var a = [];\n        for(var i=0;i<files.length;i++){\n            if(files[i]){\n                if(files[i].constructor == Array) a = a.concat(flatten(files[i]));\n                else a = a.concat(files[i]);}}\n        var b = a.filter(function (el) {return el != null;});\n        return b;}\n   \n    var input_reads = [].concat(inputs.input_reads);\n    var files=flatten(input_reads);\n\n    // Set output file name\n    var fname = '';\n    \n    /// from given prefix\n    if (inputs.output_name) fname = inputs.output_name;\n    /// from sample_id metadata\n    else if (files[0].metadata && files[0].metadata['sample_id']) fname=files[0].metadata['sample_id'];\n    /// from common prefix, and strip out any unnecessary characters\n    else {\n        var read_names = [];\n        for (var i=0;i<files.length;i++) {\n            var file_ext=files[i].nameext;\n            var file_base=files[i].basename;\n            \n            if (file_ext === '.gz' || file_ext === '.GZ')\n                file_base = file_base.slice(0, -3);\n                file_ext= '.'+ file_base.split('.').pop();\n            if (file_ext === '.fq' || file_ext === '.FQ')\n                file_base = file_base.slice(0, -3);\n            if (file_ext === '.fastq' || file_ext === '.FASTQ')\n                file_base = file_base.slice(0, -6);\n            \n            read_names.push(file_base.replace(/pe1|pe2|pe\\.1|pe\\.2|pe\\_1|pe\\_2|\\_pe1|\\_pe2|\\_pe\\.1|\\_pe\\.2|\\_pe\\_1|\\_pe\\_2|\\.pe1|\\.pe2|\\.pe\\.1|\\.pe\\.2|\\.pe\\_1|\\.pe\\_2/,''));\n              \n          }\n          ////strip out any trailing dashes/dots/underscores...\n          var unique_prefix = sharedStart(read_names).replace( /\\-$|\\_$|\\.$/, '');\n          var tmp_prefix = unique_prefix.replace( /^\\_|\\.pe$|\\.R$|\\_pe$|\\_R$/,'');\n          var final_prefix = tmp_prefix.replace( /^_\\d(\\d)?_/, '' );\n          \n          fname=final_prefix;}\n\n\n    // Read number of threads if defined\n    var threads = 0;\n    var MAX_THREADS = 0;\n    var ref_name = '';\n    if (inputs.threads) threads = inputs.threads;\n    else if (inputs.wgs_hg38_mode_threads) {\n        MAX_THREADS = 36;\n        ref_name = inputs.reference_index_tar.basename;\n        if (ref_name.search('38') >= 0) threads = inputs.wgs_hg38_mode_threads;\n        else threads = MAX_THREADS;\n        } \n    else threads = 8;\n\n    var tool = '';\n    var dedup = '';\n    if (inputs.deduplication == \"MarkDuplicates\") {\n        tool = 'bamsormadup';\n        dedup = ' markduplicates=1';\n    } else {\n        if (inputs.output_format == 'CRAM') tool = 'bamsort index=0';\n        else tool = 'bamsort index=1';\n        if (inputs.deduplication == \"RemoveDuplicates\") dedup = ' rmdup=1';\n        else dedup = '';\n    }\n    var sort_path = tool + dedup;\n\n    var indexfilename = '';\n    var out_format = '';\n    var extension  = '';\n    // Coordinate Sorted BAM is default\n    if (inputs.output_format == 'CRAM') {\n        out_format = ' outputformat=cram SO=coordinate';\n        ref_name = inputs.reference_index_tar.basename.split('.tar')[0];\n        out_format += ' reference=' + ref_name;\n        if (sort_path != 'bamsort index=0') indexfilename = ' indexfilename=' + fname + '.cram.crai';\n        extension = '.cram';\n    } else if (inputs.output_format == 'SAM') {\n        out_format = ' outputformat=sam SO=coordinate';\n        extension = '.sam';\n    } else if (inputs.output_format == 'Queryname Sorted BAM') {\n        out_format = ' outputformat=bam SO=queryname';\n        extension = '.bam';\n    } else if (inputs.output_format == 'Queryname Sorted SAM') {\n        out_format = ' outputformat=sam SO=queryname';\n        extension = '.sam';\n    } else {\n        out_format = ' outputformat=bam SO=coordinate';\n        indexfilename = ' indexfilename=' + fname + '.bam.bai';\n        extension = '.bam';\n    }\n    var cmd = \" | \" + sort_path + \" threads=\" + threads + \" level=1 tmplevel=-1 inputformat=sam\";\n    cmd += out_format;\n    cmd += indexfilename;\n    // capture metrics file\n    cmd += \" M=\" + fname + \".sormadup_metrics.log\";\n\n    if (inputs.output_format == 'SAM') cmd = '';\n    \n    return cmd + ' > ' + fname + extension;\n    \n}"
                    },
                    {
                        "prefix": "",
                        "shellQuote": false,
                        "position": 10004,
                        "valueFrom": "${\n    /// Get pipe status ///\n    \n    var  cmd = \";declare -i pipe_statuses=(\\\\${PIPESTATUS[*]});len=\\\\${#pipe_statuses[@]};declare -i tot=0;echo \\\\${pipe_statuses[*]};for (( i=0; i<\\\\${len}; i++ ));do if [ \\\\${pipe_statuses[\\\\$i]} -ne 0 ];then tot=\\\\${pipe_statuses[\\\\$i]}; fi;done;if [ \\\\$tot -ne 0 ]; then >&2 echo Error in piping. Pipe statuses: \\\\${pipe_statuses[*]};fi; if [ \\\\$tot -ne 0 ]; then false;fi\\\"\";\n    return cmd;\n}"
                    }
                ],
                "requirements": [
                    {
                        "class": "ShellCommandRequirement"
                    },
                    {
                        "class": "ResourceRequirement",
                        "ramMin": "${\n    var reads_size =0;\n    // Calculate suggested number of CPUs depending of the input reads size\n    if (inputs.input_reads.constructor == Array) {\n        if (inputs.input_reads[1]) reads_size = inputs.input_reads[0].size + inputs.input_reads[1].size;\n        else reads_size = inputs.input_reads[0].size;\n    } else reads_size = inputs.input_reads.size;\n    if (!reads_size) reads_size = 0;\n\n    var GB_1 = 1024 * 1024 * 1024;\n    var  suggested_memory = 0;\n    if (reads_size < GB_1) suggested_memory = 4;\n    else if (reads_size < 10 * GB_1) suggested_memory = 15;\n    else suggested_memory = 58;\n    \n    if (inputs.total_memory) return inputs.total_memory * 1024;\n    else if (inputs.sort_memory) return inputs.sort_memory * 1024;\n    else return suggested_memory * 1024;\n    \n}",
                        "coresMin": "${\n    var reads_size = 0\n    // Calculate suggested number of CPUs depending of the input reads size\n    if (inputs.input_reads.constructor == Array) {\n        if (inputs.input_reads[1]) reads_size = inputs.input_reads[0].size + inputs.input_reads[1].size;\n        else reads_size = inputs.input_reads[0].size;\n    } else reads_size = inputs.input_reads.size;\n    \n    if (!reads_size) reads_size = 0;\n    \n    var GB_1 = 1024 * 1024 * 1024;\n    var suggested_cpus = 0;\n    if (reads_size < GB_1) suggested_cpus = 1;\n    else if (reads_size < 10 * GB_1) suggested_cpus = 8;\n    else suggested_cpus = 31;\n    \n    if (inputs.reserved_threads) return inputs.reserved_threads;\n    else if (inputs.threads) return inputs.threads;\n    else if (inputs.sambamba_threads) return inputs.sambamba_threads;\n    else return suggested_cpus;\n    \n}"
                    },
                    {
                        "class": "DockerRequirement",
                        "dockerPull": "images.sbgenomics.com/ana_stankovic/bwa_0.7.17_biobambam2:0"
                    },
                    {
                        "class": "InitialWorkDirRequirement",
                        "listing": [
                            "$(inputs.reference_index_tar)",
                            "$(inputs.input_reads)",
                            "$(inputs.fasta_index)"
                        ]
                    },
                    {
                        "class": "InlineJavascriptRequirement",
                        "expressionLib": [
                            "var setMetadata = function(file, metadata) {\n    if (!('metadata' in file)) {\n        file['metadata'] = {}\n    }\n    for (var key in metadata) {\n        file['metadata'][key] = metadata[key];\n    }\n    return file\n};\n\nvar inheritMetadata = function(o1, o2) {\n    var commonMetadata = {};\n    if (!o2) {\n        return o1;\n    };\n    if (!Array.isArray(o2)) {\n        o2 = [o2]\n    }\n    for (var i = 0; i < o2.length; i++) {\n        var example = o2[i]['metadata'];\n        for (var key in example) {\n            if (i == 0)\n                commonMetadata[key] = example[key];\n            else {\n                if (!(commonMetadata[key] == example[key])) {\n                    delete commonMetadata[key]\n                }\n            }\n        }\n        for (var key in commonMetadata) {\n            if (!(key in example)) {\n                delete commonMetadata[key]\n            }\n        }\n    }\n    if (!Array.isArray(o1)) {\n        o1 = setMetadata(o1, commonMetadata)\n        if (o1.secondaryFiles) {\n            o1.secondaryFiles = inheritMetadata(o1.secondaryFiles, o2)\n        }\n    } else {\n        for (var i = 0; i < o1.length; i++) {\n            o1[i] = setMetadata(o1[i], commonMetadata)\n            if (o1[i].secondaryFiles) {\n                o1[i].secondaryFiles = inheritMetadata(o1[i].secondaryFiles, o2)\n            }\n        }\n    }\n    return o1;\n};\n"
                        ]
                    }
                ],
                "sbg:toolAuthor": "Heng Li",
                "sbg:toolkit": "BWA",
                "sbg:license": "BWA: GNU Affero General Public License v3.0, MIT License; Biobambam2: GNU General Public License v3.0",
                "sbg:categories": [
                    "Alignment"
                ],
                "sbg:toolkitVersion": "0.7.17",
                "sbg:image_url": null,
                "sbg:cmdPreview": "/bin/bash -c \" export REF_CACHE=${PWD} ;  tar -tvf reference.HG38.fasta.gz.tar 1>&2; tar -xf reference.HG38.fasta.gz.tar ;  bwa mem  -R '@RG\\tID:1\\tPL:Illumina\\tSM:dnk_sample' -t 10  reference.HG38.fasta.gz  /path/to/LP6005524-DNA_C01_lane_7.sorted.converted.filtered.pe_2.gz /path/to/LP6005524-DNA_C01_lane_7.sorted.converted.filtered.pe_1.gz  | bamsormadup threads=10 level=1 tmplevel=-1 inputformat=sam outputformat=cram SO=coordinate reference=reference.HG38.fasta.gz indexfilename=LP6005524-DNA_C01_lane_7.sorted.converted.filtered.cram.crai M=LP6005524-DNA_C01_lane_7.sorted.converted.filtered.sormadup_metrics.log > LP6005524-DNA_C01_lane_7.sorted.converted.filtered.cram  ;declare -i pipe_statuses=(\\${PIPESTATUS[*]});len=\\${#pipe_statuses[@]};declare -i tot=0;echo \\${pipe_statuses[*]};for (( i=0; i<\\${len}; i++ ));do if [ \\${pipe_statuses[\\$i]} -ne 0 ];then tot=\\${pipe_statuses[\\$i]}; fi;done;if [ \\$tot -ne 0 ]; then >&2 echo Error in piping. Pipe statuses: \\${pipe_statuses[*]};fi; if [ \\$tot -ne 0 ]; then false;fi\"",
                "sbg:links": [
                    {
                        "id": "http://bio-bwa.sourceforge.net/",
                        "label": "Homepage"
                    },
                    {
                        "id": "https://github.com/lh3/bwa",
                        "label": "Source code"
                    },
                    {
                        "id": "http://bio-bwa.sourceforge.net/bwa.shtml",
                        "label": "Wiki"
                    },
                    {
                        "id": "http://sourceforge.net/projects/bio-bwa/",
                        "label": "Download"
                    },
                    {
                        "id": "http://arxiv.org/abs/1303.3997",
                        "label": "Publication"
                    },
                    {
                        "id": "http://www.ncbi.nlm.nih.gov/pubmed/19451168",
                        "label": "Publication BWA Algorithm"
                    }
                ],
                "sbg:projectName": "SBG Public data",
                "sbg:revisionsInfo": [
                    {
                        "sbg:revision": 0,
                        "sbg:modifiedBy": "admin",
                        "sbg:modifiedOn": 1652862504,
                        "sbg:revisionNotes": null
                    },
                    {
                        "sbg:revision": 1,
                        "sbg:modifiedBy": "admin",
                        "sbg:modifiedOn": 1652862504,
                        "sbg:revisionNotes": "Description updated"
                    }
                ],
                "sbg:expand_workflow": false,
                "sbg:appVersion": [
                    "v1.2"
                ],
                "sbg:id": "admin/sbg-public-data/bwa-mem-bundle-0-7-17-cwl-1-2/1",
                "sbg:revision": 1,
                "sbg:revisionNotes": "Description updated",
                "sbg:modifiedOn": 1652862504,
                "sbg:modifiedBy": "admin",
                "sbg:createdOn": 1652862504,
                "sbg:createdBy": "admin",
                "sbg:project": "admin/sbg-public-data",
                "sbg:sbgMaintained": false,
                "sbg:validationErrors": [],
                "sbg:contributors": [
                    "admin"
                ],
                "sbg:latestRevision": 1,
                "sbg:publisher": "sbg",
                "sbg:content_hash": "a7e4a5ee4873c491fb76818a0c61f0c2f24f5723335fe79599c29d8896f71f133",
                "sbg:workflowLanguage": "CWL"
            },
            "label": "BWA MEM",
            "sbg:x": 61.430381774902344,
            "sbg:y": 431.68353271484375
        },
        {
            "id": "bwa_index_0_7_17_cwl_1_2",
            "in": [
                {
                    "id": "reference",
                    "source": "reference"
                }
            ],
            "out": [
                {
                    "id": "indexed_reference"
                }
            ],
            "run": {
                "class": "CommandLineTool",
                "cwlVersion": "v1.2",
                "$namespaces": {
                    "sbg": "https://sevenbridges.com"
                },
                "id": "admin/sbg-public-data/bwa-index-0-7-17-cwl-1-2/2",
                "baseCommand": [],
                "inputs": [
                    {
                        "sbg:category": "Configuration",
                        "sbg:toolDefaultValue": "auto",
                        "id": "bwt_construction",
                        "type": [
                            "null",
                            {
                                "type": "enum",
                                "symbols": [
                                    "bwtsw",
                                    "is",
                                    "rb2"
                                ],
                                "name": "bwt_construction"
                            }
                        ],
                        "label": "Bwt construction",
                        "doc": "Algorithm for constructing BWT index. \nAvailable options are: \nis - IS linear-time algorithm for constructing suffix array. It requires 5.37N memory where N is the size of the database. IS is moderately fast but does not work with a database larger than 2GB. IS is the default algorithm due to its simplicity. The current codes for the IS algorithm are reimplemented by Yuta Mori.\nbwtsw - Algorithm implemented in BWT-SW. This method works with the whole human genome.\nrb2 - ropebwt2 algorithm as an alternative to index large genomes.\nRopebwt2 is slower than the \"bwtsw\" algorithm, but it has a permissive\nlicense. This allows creating an Apache2-licensed BWA (in the \"Apache2\"\nbranch) for commercial users who are concerned with GPL.\nIf none is selected, the tool itself will choose the correct method automatically."
                    },
                    {
                        "sbg:category": "Configuration",
                        "id": "prefix_of_the_index_to_be_output",
                        "type": "string?",
                        "label": "Prefix of the indices to be outputted",
                        "doc": "The prefix of the generated indices. The same prefix will be set for the output TAR archive.\nThe default prefix is the input reference file name."
                    },
                    {
                        "sbg:category": "Configuration",
                        "sbg:toolDefaultValue": "10000000",
                        "id": "block_size",
                        "type": "int?",
                        "label": "Block size",
                        "doc": "Block size for the bwtsw algorithm (effective with -a bwtsw)."
                    },
                    {
                        "sbg:category": "Configuration",
                        "id": "add_64_to_fasta_name",
                        "type": "boolean?",
                        "label": "Add 64 to output index file names",
                        "doc": "Output index files will be renamed by adding 64 (in.fasta.64.* instead of in.fasta.*)."
                    },
                    {
                        "sbg:category": "File input",
                        "id": "reference",
                        "type": "File",
                        "label": "Reference file",
                        "doc": "Input reference FASTA or TAR file with reference and indices.",
                        "sbg:fileTypes": "FASTA, FA, FA.GZ, FASTA.GZ, TAR"
                    },
                    {
                        "sbg:category": "Configuration",
                        "id": "total_memory",
                        "type": "int?",
                        "label": "Total memory [Gb]",
                        "doc": "Total memory [GB] to be reserved for the tool. The default value is 1.5 x the size of the reference file."
                    },
                    {
                        "sbg:toolDefaultValue": "False",
                        "id": "do_not_add_alt_contig_to_reference",
                        "type": "boolean?",
                        "label": "Do not add alt contigs",
                        "doc": "Do not add alt contigs file to the TAR bundle."
                    }
                ],
                "outputs": [
                    {
                        "id": "indexed_reference",
                        "doc": "TAR archive containing the reference FASTA file with its generated BWA indices.",
                        "label": "TAR archive with FASTA and its BWA indices",
                        "type": "File?",
                        "outputBinding": {
                            "glob": "*.tar",
                            "outputEval": "${\n    var in_file=[].concat(inputs.reference)\n    self = inheritMetadata(self, in_file[0]);\n\n    return self;\n}"
                        },
                        "sbg:fileTypes": "TAR"
                    }
                ],
                "doc": "**BWA INDEX** constructs the FM-index (Full-text index in Minute space) for the **BWA MEM Bundle 0.7.17 CWL1.2** aligner. The FM-index of a reference sequence is based on its Burrows-Wheeler Transform (BWT) and suffix array. Generated index files will be used with BWA MEM, BWA ALN, BWA SAMPE, and BWA SAMSE tools.\n\n*A list of **all inputs and parameters** with corresponding descriptions can be found at the bottom of the page.*\n\n### Common Use Cases\n\n* **BWA INDEX** is used for indexing the reference sequence, and it is a prior step needed for a **BWA MEM Bundle 0.7.17 CWL1.2** run.\n\n### Changes Introduced by Seven Bridges\n\n- Generated files are saved to an output TAR archive with the BWA-0.7.17-INDEX-ARCHIVE.TAR extension.\n- If the input reference file has a TAR extension, it is assumed that the BWA-MEM indices came together with it, and **BWA INDEX** will only pass that TAR through to the output.  If the input is not a TAR archive, indexing of the reference sequence will be performed. This is used in workflows, to skip unnecessary steps if indices are already available, and to enable easy indexing if it is needed.\n\n\n### Common Issues and Important Notes \n\n* If the input reference file has a TAR extension it is assumed that BWA indices came together with it and **BWA INDEX** will only pass that TAR to the output. If the input is not a TAR file, the indexing of the reference sequence will be performed.\n* TAR contains an alt reference from bwa.kit suggested for the HG38 reference genome by the tool author.\n* BWA-MEM2 is not a splice-aware aligner, so it is not the appropriate tool for mapping RNAseq to the genome. For RNAseq reads **Bowtie2 Aligner** and **STAR** are recommended tools. \n\n### Performance Benchmarking\n\n| Reference file |Size [GB] | Instance (AWS)| Duration | Cost |RAM peak [GB]|\n| --- | --- | --- | --- | --- | --- |\n| Homo_sapiens_assembly38.fasta| 3 | c4.2xlarge (8vCPU, 16GB RAM)|1h, 28min  | $0.79 |6.2|\n\n*Cost can be significantly reduced by using **spot instances**. Visit the [Knowledge Center](https://docs.sevenbridges.com/docs/about-spot-instances) for more details.*\n\n### Portability\n\n**BWA INDEX** was tested with cwltool version 3.0.20201203173111.",
                "label": "BWA INDEX",
                "arguments": [
                    {
                        "prefix": "",
                        "shellQuote": false,
                        "position": 0,
                        "valueFrom": "${\n    /// Create command or pass index files without processing ///\n    var in_reference_file = [].concat(inputs.reference)[0]\n    var reference_file = in_reference_file.basename;\n    var ext = in_reference_file.nameext;\n\n    if (ext == '.tar') {\n        return 'echo Index files passed without any processing!';\n    } else {\n        var cp_alt_cmd = '';\n        if (!inputs.do_not_add_alt_contig_to_reference) {\n            if (reference_file.search('38') >= 0) {\n                cp_alt_cmd = 'cp /opt/hs38DH.fa.alt ' + reference_file + '.alt && ';\n            }\n        }\n\n        var index_cmd = 'bwa index ' + reference_file + ' ';\n\n        return cp_alt_cmd + index_cmd;\n    }\n}"
                    },
                    {
                        "prefix": "",
                        "shellQuote": false,
                        "position": 1,
                        "valueFrom": "${\n    /// Bwt construction or pass index files without processing ///\n    var in_reference_file = [].concat(inputs.reference)[0]\n    var reference_file = in_reference_file.basename;\n    var ext = in_reference_file.nameext;\n    \n    if (ext == '.tar' || !inputs.bwt_construction) return '';\n    else return '-a ' + inputs.bwt_construction;\n}"
                    },
                    {
                        "prefix": "",
                        "shellQuote": false,
                        "position": 1,
                        "valueFrom": "${\n    /// Set prefix or pass index files without processing ///\n    var in_reference_file = [].concat(inputs.reference)[0]\n    var reference_file = in_reference_file.basename;\n    var ext = in_reference_file.nameext;\n    \n    if (ext == '.tar' || !inputs.prefix_of_the_index_to_be_output) return '';\n    else return '-p ' + inputs.prefix_of_the_index_to_be_output;\n    \n}"
                    },
                    {
                        "prefix": "",
                        "shellQuote": false,
                        "position": 1,
                        "valueFrom": "${\n    /// Block size or pass index files without processing ///\n    var in_reference_file = [].concat(inputs.reference)[0]\n    var reference_file = in_reference_file.basename;\n    var ext = in_reference_file.nameext;\n    \n    if (ext == '.tar' || !inputs.block_size) return '';\n    else return '-b ' + inputs.block_size;\n\n}"
                    },
                    {
                        "prefix": "",
                        "shellQuote": false,
                        "position": 1,
                        "valueFrom": "${\n    /// Add 64 to fasta name or pass index files without processing ///\n    var in_reference_file = [].concat(inputs.reference)[0]\n    var reference_file = in_reference_file.basename;\n    var ext = in_reference_file.nameext;\n    \n    if (ext == '.tar' || !inputs.add_64_to_fasta_name) return '';\n    else return '-6 ';\n    \n}"
                    },
                    {
                        "prefix": "",
                        "shellQuote": false,
                        "position": 1,
                        "valueFrom": "${\n    /// Pack indices into TAR archive or pass index files without processing ///\n    var in_reference_file = [].concat(inputs.reference)[0]\n    var reference_file = in_reference_file.basename;\n    var ext = in_reference_file.nameext;\n    \n    if (ext == '.tar') return '';\n    else {\n        var extensions = ' *.amb' + ' *.ann' + ' *.bwt' + ' *.pac' + ' *.sa';\n        \n        if (!inputs.do_not_add_alt_contig_to_reference) {\n            if (reference_file.search('38') >= 0) extensions = extensions + ' *.alt ';\n        }\n        \n        var out_name = '';\n        if (inputs.prefix_of_the_index_to_be_output) out_name = inputs.prefix_of_the_index_to_be_output;\n        else out_name = inputs.reference.nameroot;\n        \n        var tar_cmd = 'tar -cf ' + out_name + '.bwa-0.7.17-index-archive.tar ' + reference_file + extensions;\n        return ' && ' + tar_cmd;\n    }\n}"
                    }
                ],
                "requirements": [
                    {
                        "class": "ShellCommandRequirement"
                    },
                    {
                        "class": "ResourceRequirement",
                        "ramMin": "${\n    var in_reference_file = [].concat(inputs.reference)[0]\n    var reference_file = in_reference_file.basename;\n    var ext = in_reference_file.nameext;\n    \n    var GB_1 = 1024 * 1024 * 1024;\n    var reads_size = reference_file.size;\n    \n    if (!reads_size) reads_size = GB_1;\n    \n    if (inputs.total_memory) return inputs.total_memory * 1024;\n    else if (ext == '.tar') return 128;\n    else return (parseInt(1.5 * reads_size / (1024 * 1024)));\n    \n}",
                        "coresMin": 1
                    },
                    {
                        "class": "DockerRequirement",
                        "dockerPull": "images.sbgenomics.com/ana_stankovic/bwa_0.7.17:0"
                    },
                    {
                        "class": "InitialWorkDirRequirement",
                        "listing": [
                            "$(inputs.reference)"
                        ]
                    },
                    {
                        "class": "InlineJavascriptRequirement",
                        "expressionLib": [
                            "var setMetadata = function(file, metadata) {\n    if (!('metadata' in file)) {\n        file['metadata'] = {}\n    }\n    for (var key in metadata) {\n        file['metadata'][key] = metadata[key];\n    }\n    return file\n};\n\nvar inheritMetadata = function(o1, o2) {\n    var commonMetadata = {};\n    if (!o2) {\n        return o1;\n    };\n    if (!Array.isArray(o2)) {\n        o2 = [o2]\n    }\n    for (var i = 0; i < o2.length; i++) {\n        var example = o2[i]['metadata'];\n        for (var key in example) {\n            if (i == 0)\n                commonMetadata[key] = example[key];\n            else {\n                if (!(commonMetadata[key] == example[key])) {\n                    delete commonMetadata[key]\n                }\n            }\n        }\n        for (var key in commonMetadata) {\n            if (!(key in example)) {\n                delete commonMetadata[key]\n            }\n        }\n    }\n    if (!Array.isArray(o1)) {\n        o1 = setMetadata(o1, commonMetadata)\n        if (o1.secondaryFiles) {\n            o1.secondaryFiles = inheritMetadata(o1.secondaryFiles, o2)\n        }\n    } else {\n        for (var i = 0; i < o1.length; i++) {\n            o1[i] = setMetadata(o1[i], commonMetadata)\n            if (o1[i].secondaryFiles) {\n                o1[i].secondaryFiles = inheritMetadata(o1[i].secondaryFiles, o2)\n            }\n        }\n    }\n    return o1;\n};\n"
                        ]
                    }
                ],
                "sbg:toolkit": "BWA",
                "sbg:links": [
                    {
                        "id": "http://bio-bwa.sourceforge.net/",
                        "label": "Homepage"
                    },
                    {
                        "id": "https://github.com/lh3/bwa",
                        "label": "Source Code"
                    },
                    {
                        "id": "http://bio-bwa.sourceforge.net/bwa.shtml",
                        "label": "Wiki"
                    },
                    {
                        "id": "http://sourceforge.net/projects/bio-bwa/",
                        "label": "Download"
                    },
                    {
                        "id": "http://www.ncbi.nlm.nih.gov/pubmed/19451168",
                        "label": "Publication"
                    }
                ],
                "sbg:image_url": null,
                "sbg:cmdPreview": "cp /opt/hs38DH.fa.alt reference38.fasta.alt ; bwa index reference38.fasta            ; tar -cf reference38.fasta.tar reference38.fasta *.amb *.ann *.bwt *.pac *.sa *.alt ;",
                "sbg:license": "GNU Affero General Public License v3.0, MIT License",
                "sbg:toolAuthor": "Heng Li",
                "sbg:toolkitVersion": "0.7.17",
                "sbg:categories": [
                    "Alignment",
                    "Indexing"
                ],
                "sbg:projectName": "SBG Public data",
                "sbg:revisionsInfo": [
                    {
                        "sbg:revision": 0,
                        "sbg:modifiedBy": "admin",
                        "sbg:modifiedOn": 1652862595,
                        "sbg:revisionNotes": null
                    },
                    {
                        "sbg:revision": 1,
                        "sbg:modifiedBy": "admin",
                        "sbg:modifiedOn": 1652862596,
                        "sbg:revisionNotes": "Description updated"
                    },
                    {
                        "sbg:revision": 2,
                        "sbg:modifiedBy": "admin",
                        "sbg:modifiedOn": 1652862596,
                        "sbg:revisionNotes": "Docker updated"
                    }
                ],
                "sbg:expand_workflow": false,
                "sbg:appVersion": [
                    "v1.2"
                ],
                "sbg:id": "admin/sbg-public-data/bwa-index-0-7-17-cwl-1-2/2",
                "sbg:revision": 2,
                "sbg:revisionNotes": "Docker updated",
                "sbg:modifiedOn": 1652862596,
                "sbg:modifiedBy": "admin",
                "sbg:createdOn": 1652862595,
                "sbg:createdBy": "admin",
                "sbg:project": "admin/sbg-public-data",
                "sbg:sbgMaintained": false,
                "sbg:validationErrors": [],
                "sbg:contributors": [
                    "admin"
                ],
                "sbg:latestRevision": 2,
                "sbg:publisher": "sbg",
                "sbg:content_hash": "a7f3ec28c14e4c807fcc2025187bd802acb2b41b970fe72b73d2c74d26026f56b",
                "sbg:workflowLanguage": "CWL"
            },
            "label": "BWA INDEX",
            "sbg:x": -198.76397705078125,
            "sbg:y": 291.81365966796875
        }
    ],
    "requirements": [
        {
            "class": "InlineJavascriptRequirement"
        },
        {
            "class": "StepInputExpressionRequirement"
        }
    ],
    "sbg:projectName": "Simulation",
    "sbg:revisionsInfo": [
        {
            "sbg:revision": 0,
            "sbg:modifiedBy": "hendrick.san",
            "sbg:modifiedOn": 1655774064,
            "sbg:revisionNotes": null
        },
        {
            "sbg:revision": 1,
            "sbg:modifiedBy": "hendrick.san",
            "sbg:modifiedOn": 1655782912,
            "sbg:revisionNotes": "Initial"
        },
        {
            "sbg:revision": 2,
            "sbg:modifiedBy": "hendrick.san",
            "sbg:modifiedOn": 1660985818,
            "sbg:revisionNotes": "output as VCF"
        },
        {
            "sbg:revision": 3,
            "sbg:modifiedBy": "hendrick.san",
            "sbg:modifiedOn": 1682950032,
            "sbg:revisionNotes": "simplified"
        },
        {
            "sbg:revision": 4,
            "sbg:modifiedBy": "hendrick.san",
            "sbg:modifiedOn": 1682953164,
            "sbg:revisionNotes": "I/O w/ file format"
        }
    ],
    "sbg:image_url": "https://cgc.sbgenomics.com/ns/brood/images/hendrick.san/simulation/cowid/4.png",
    "sbg:appVersion": [
        "v1.2",
        "sbg:draft-2"
    ],
    "id": "https://cgc-api.sbgenomics.com/v2/apps/hendrick.san/simulation/cowid/4/raw/",
    "sbg:id": "hendrick.san/simulation/cowid/4",
    "sbg:revision": 4,
    "sbg:revisionNotes": "I/O w/ file format",
    "sbg:modifiedOn": 1682953164,
    "sbg:modifiedBy": "hendrick.san",
    "sbg:createdOn": 1655774064,
    "sbg:createdBy": "hendrick.san",
    "sbg:project": "hendrick.san/simulation",
    "sbg:sbgMaintained": false,
    "sbg:validationErrors": [],
    "sbg:contributors": [
        "hendrick.san"
    ],
    "sbg:latestRevision": 4,
    "sbg:publisher": "sbg",
    "sbg:content_hash": "a7337875477bdda228f8cbb4b847bdcb1e40296135da11e7ed7e906811634a0f2",
    "sbg:workflowLanguage": "CWL"
}
