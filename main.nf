$HOSTNAME = ""
params.outdir = 'results'  

// filter_seq_quality
params.filter_seq_quality.method = "quality"
params.filter_seq_quality.nproc = params.nproc
params.filter_seq_quality.q = "20"
params.filter_seq_quality.n_length = ""
params.filter_seq_quality.n_missing = ""

// filter_seq_length
params.filter_seq_length.method = "length"
params.filter_seq_length.nproc = params.nproc
params.filter_seq_length.q = ""
params.filter_seq_length.n_length = "750"
params.filter_seq_length.n_missing = ""

// MaskPrimers_CPRIMERS
params.MaskPrimers_CPRIMERS.method = ["align"]
params.MaskPrimers_CPRIMERS.mode = ["cut"]
params.MaskPrimers_CPRIMERS.start = [0]
params.MaskPrimers_CPRIMERS.barcode = ["false"]
params.MaskPrimers_CPRIMERS.barcode_field = [""]
params.MaskPrimers_CPRIMERS.primer_field = ["CPRIMER"]
params.MaskPrimers_CPRIMERS.maxerror = [0.3]
params.MaskPrimers_CPRIMERS.revpr = ["false"]
params.MaskPrimers_CPRIMERS.maxlen = [50]
params.MaskPrimers_CPRIMERS.skiprc = ["false"]
params.MaskPrimers_CPRIMERS.failed = "true"
params.MaskPrimers_CPRIMERS.nproc = params.nproc
params.MaskPrimers_CPRIMERS.R1_primers = "${params.projectDir}/primers/CPRIMERS_${params.constant_region}.fasta"
params.MaskPrimers_CPRIMERS.R2_primers = ""

params.parse_log_MP_CPRIMERS.suffix = "_CPRIMERS"

// MaskPrimers_VPRIMERS
params.MaskPrimers_VPRIMERS.method = ["align"]
params.MaskPrimers_VPRIMERS.mode = ["cut"]
params.MaskPrimers_VPRIMERS.start = [0]
params.MaskPrimers_VPRIMERS.barcode = ["false"]
params.MaskPrimers_VPRIMERS.barcode_field = [""]
params.MaskPrimers_VPRIMERS.primer_field = ["VPRIMER"]
params.MaskPrimers_VPRIMERS.maxerror = [0.3]
params.MaskPrimers_VPRIMERS.revpr = ["false"]
params.MaskPrimers_VPRIMERS.maxlen = [50]
params.MaskPrimers_VPRIMERS.skiprc = ["false"]
params.MaskPrimers_VPRIMERS.failed = "true"
params.MaskPrimers_VPRIMERS.nproc = params.nproc
params.MaskPrimers_VPRIMERS.R1_primers = "${params.projectDir}/primers/VPRIMERS.fasta"
params.MaskPrimers_VPRIMERS.R2_primers = ""

params.parse_log_MP_VPRIMERS.suffix = "_VPRIMERS"
// MaskPrimers_EXTRACT
params.MaskPrimers_EXTRACT.method = ["extract"]
params.MaskPrimers_EXTRACT.mode = ["cut"]
params.MaskPrimers_EXTRACT.start = [0]
params.MaskPrimers_EXTRACT.extract_length = [22]
params.MaskPrimers_EXTRACT.barcode = ["false"]
params.MaskPrimers_EXTRACT.barcode_field = [""]
params.MaskPrimers_EXTRACT.primer_field = ["BARCODE"]
params.MaskPrimers_EXTRACT.maxerror = [0.3]
params.MaskPrimers_EXTRACT.revpr = ["false"]
params.MaskPrimers_EXTRACT.maxlen = [50]
params.MaskPrimers_EXTRACT.skiprc = ["false"]
params.MaskPrimers_EXTRACT.failed = "true"
params.MaskPrimers_EXTRACT.nproc = params.nproc
params.MaskPrimers_EXTRACT.R1_primers = ""
params.MaskPrimers_EXTRACT.R2_primers = ""

params.parse_log_MP_EXTRACT.suffix = "_EXTRACT"

params.check_for_seqs.primers_file = "${params.projectDir}/primers/VPRIMERS.fasta"

// align_sets
params.align_sets.method = "muscle"
params.align_sets.bf = "BARCODE"
params.align_sets.div = "fasle"
params.align_sets.failed = "false"
params.align_sets.nproc = params.nproc
params.align_sets.muscle_exec = "/usr/local/bin/muscle"
params.align_sets.offset_table = ""
params.align_sets.pf = ""
params.align_sets.mode = ""
params.align_sets.primer_file = ""
params.align_sets.reverse = "false"

// cluster_sets
params.cluster_sets.method = ["set"]
params.cluster_sets.failed = ["false"]
params.cluster_sets.nproc = params.nproc
params.cluster_sets.cluster_field = ["CLUSTER"]
params.cluster_sets.cluster_tool = ["usearch"]
params.cluster_sets.cluster_exec = ["/usr/local/bin/usearch"]
params.cluster_sets.set_field = ["BARCODE"]
params.cluster_sets.barcode_field = ["BARCODE"]

// parse_headers_copy
params.parse_headers_copy.method = "copy"
params.parse_headers_copy.act = "cat"
params.parse_headers_copy.args = "-f BARCODE -k CLUSTER"

// build_consensus
params.build_consensus.failed = "true"
params.build_consensus.nproc = params.nproc
params.build_consensus.barcode_field = ["CLUSTER"]
params.build_consensus.primer_field = ["CPRIMER"]
params.build_consensus.act = ["none"]
params.build_consensus.copy_field = [""]
params.build_consensus.maxerror = ["0.1"]
params.build_consensus.prcons = ["0.6"]
params.build_consensus.maxgap = ["0.5"]
params.build_consensus.maxdiv = ["none"]
params.build_consensus.dep = ["false"]

// parse_headers_copy
params.parse_headers_collapse.method = "collapse"
params.parse_headers_collapse.act = "min"
params.parse_headers_collapse.args = "-f CONSCOUNT"

// collapse_seq
params.collapse_seq.max_missing = 20 
params.collapse_seq.inner = "true" 
params.collapse_seq.fasta = "true" 
params.collapse_seq.act = "sum" 
params.collapse_seq.uf = "CREGION" 
params.collapse_seq.cf = "CONSCOUNT" 
params.collapse_seq.nproc = params.nproc
params.collapse_seq.failed = "true"

// split_seq
params.split_seq.field = "CONSCOUNT"
params.split_seq.num = 2
params.split_seq.fasta = "true"

// parse_headers_table
params.parse_headers_table.method = "table"
params.parse_headers_table.args = "-f ID PRCONS CONSCOUNT DUPCOUNT"

params.Alignment_FLAIRSEQ_IgBlastn.num_threads = params.nproc
params.Alignment_FLAIRSEQ_IgBlastn.ig_seqtype = "Ig"
params.Alignment_FLAIRSEQ_IgBlastn.outfmt = "MakeDb"
params.Alignment_FLAIRSEQ_IgBlastn.num_alignments_V = "10"
params.Alignment_FLAIRSEQ_IgBlastn.domain_system = "imgt"


params.Alignment_FLAIRSEQ_MakeDb.failed = "true"
params.Alignment_FLAIRSEQ_MakeDb.format = "airr"
params.Alignment_FLAIRSEQ_MakeDb.regions = "default"
params.Alignment_FLAIRSEQ_MakeDb.extended = "true"
params.Alignment_FLAIRSEQ_MakeDb.asisid = "false"
params.Alignment_FLAIRSEQ_MakeDb.asiscalls = "false"
params.Alignment_FLAIRSEQ_MakeDb.inferjunction = "false"
params.Alignment_FLAIRSEQ_MakeDb.partial = "false"
params.Alignment_FLAIRSEQ_MakeDb.name_alignment = "_Alignment_FLAIRSEQ"

// Process Parameters for Alignment_FLAIRSEQ_Collapse_AIRRseq:
params.Alignment_FLAIRSEQ_Collapse_AIRRseq.conscount_min = 2
params.Alignment_FLAIRSEQ_Collapse_AIRRseq.n_max = 10
params.Alignment_FLAIRSEQ_Collapse_AIRRseq.name_alignment = "_Alignment_FLAIRSEQ"

// Process Parameters for Clone_AIRRseq_First_CreateGermlines:
params.Clone_AIRRseq_First_CreateGermlines.failed = "false"
params.Clone_AIRRseq_First_CreateGermlines.format = "airr"
params.Clone_AIRRseq_First_CreateGermlines.g = "dmask"
params.Clone_AIRRseq_First_CreateGermlines.cloned = "false"
params.Clone_AIRRseq_First_CreateGermlines.seq_field = ""
params.Clone_AIRRseq_First_CreateGermlines.v_field = ""
params.Clone_AIRRseq_First_CreateGermlines.d_field = ""
params.Clone_AIRRseq_First_CreateGermlines.j_field = ""
params.Clone_AIRRseq_First_CreateGermlines.clone_field = ""

params.Clone_AIRRseq_DefineClones.failed = "false"
params.Clone_AIRRseq_DefineClones.format = "airr"
params.Clone_AIRRseq_DefineClones.seq_field = ""
params.Clone_AIRRseq_DefineClones.v_field = ""
params.Clone_AIRRseq_DefineClones.d_field = ""
params.Clone_AIRRseq_DefineClones.j_field = ""
params.Clone_AIRRseq_DefineClones.group_fields =  ""
params.Clone_AIRRseq_DefineClones.mode = "gene"
params.Clone_AIRRseq_DefineClones.dist = "0.2"
params.Clone_AIRRseq_DefineClones.norm = "len"
params.Clone_AIRRseq_DefineClones.act = "set"
params.Clone_AIRRseq_DefineClones.model = "hh_s5f"
params.Clone_AIRRseq_DefineClones.sym = "min"
params.Clone_AIRRseq_DefineClones.link = "single"
params.Clone_AIRRseq_DefineClones.maxmiss = "0"

// Process Parameters for Clone_AIRRseq_Second_CreateGermlines:
params.Clone_AIRRseq_Second_CreateGermlines.failed = "false"
params.Clone_AIRRseq_Second_CreateGermlines.format = "airr"
params.Clone_AIRRseq_Second_CreateGermlines.g = "dmask"
params.Clone_AIRRseq_Second_CreateGermlines.cloned = "true"
params.Clone_AIRRseq_Second_CreateGermlines.seq_field = ""
params.Clone_AIRRseq_Second_CreateGermlines.v_field = ""
params.Clone_AIRRseq_Second_CreateGermlines.d_field = ""
params.Clone_AIRRseq_Second_CreateGermlines.j_field = ""
params.Clone_AIRRseq_Second_CreateGermlines.clone_field = ""



if (!params.mate){params.mate = ""} 
if (!params.reads){params.reads = ""} 
if (!params.d_germline){params.d_germline = ""} 
if (!params.j_germline){params.j_germline = ""} 
if (!params.v_germline){params.v_germline = ""} 
if (!params.alignment_mate){params.alignment_mate = ""} 
if (!params.auxiliary_data){params.auxiliary_data = ""} 
if (!params.custom_internal_data){params.custom_internal_data = ""} 
if (!params.c_germline){params.c_germline = ""} 
// Stage empty file to be used as an optional input where required
ch_empty_file_1 = file("$baseDir/.emptyfiles/NO_FILE_1", hidden:true)
ch_empty_file_2 = file("$baseDir/.emptyfiles/NO_FILE_2", hidden:true)
ch_empty_file_3 = file("$baseDir/.emptyfiles/NO_FILE_3", hidden:true)

Channel.value(params.mate).into{g_6_mate_g_1;g_6_mate_g_0;g_6_mate_g_4;g_6_mate_g_8;g_6_mate_g_10;g_6_mate_g_9;g_6_mate_g_12;g_6_mate_g_11;g_6_mate_g_27;g_6_mate_g_28;g_6_mate_g_31;g_6_mate_g_32;g_6_mate_g_17;g_6_mate_g_21;g_6_mate_g_20;g_6_mate_g_36;g_6_mate_g48_19;g_6_mate_g48_12}
if (params.reads){
Channel
	.fromFilePairs( params.reads , size: params.mate == "single" ? 1 : params.mate == "pair" ? 2 : params.mate == "triple" ? 3 : params.mate == "quadruple" ? 4 : -1 )
	.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
	.set{g_7_reads_g_0}
 } else {  
	g_7_reads_g_0 = Channel.empty()
 }

Channel.fromPath(params.d_germline, type: 'any').map{ file -> tuple(file.baseName, file) }.into{g_65_germlineFastaFile_g48_16;g_65_germlineFastaFile_g48_12;g_65_germlineFastaFile_g57_0;g_65_germlineFastaFile_g57_1}
Channel.fromPath(params.j_germline, type: 'any').map{ file -> tuple(file.baseName, file) }.into{g_66_germlineFastaFile_g48_17;g_66_germlineFastaFile_g48_12;g_66_germlineFastaFile_g57_0;g_66_germlineFastaFile_g57_1}
Channel.fromPath(params.v_germline, type: 'any').map{ file -> tuple(file.baseName, file) }.into{g_67_germlineFastaFile_g48_22;g_67_germlineFastaFile_g48_12;g_67_germlineFastaFile_g57_0;g_67_germlineFastaFile_g57_1}
Channel.value(params.alignment_mate).into{g_68_mate_g48_19;g_68_mate_g48_12}
g_70_outputFileTxt_g48_9 = file(params.custom_internal_data, type: 'any')


process filter_seq_quality {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /FS_.*$/) "reports/$filename"}
input:
 set val(name),file(reads) from g_7_reads_g_0
 val mate from g_6_mate_g_0

output:
 set val(name), file("*_${method}-pass.fast*")  into g_0_reads0_g_1
 set val(name), file("FS_*")  into g_0_logFile1_g_12
 set val(name), file("*_${method}-fail.fast*") optional true  into g_0_reads22
 set val(name),file("out*") optional true  into g_0_logFile33

script:
method = params.filter_seq_quality.method
nproc = params.filter_seq_quality.nproc
q = params.filter_seq_quality.q
n_length = params.filter_seq_quality.n_length
n_missing = params.filter_seq_quality.n_missing
fasta = params.filter_seq_quality.fasta
//* @style @condition:{method="quality",q}, {method="length",n_length}, {method="missing",n_missing} @multicolumn:{method,nproc}

if(method=="missing"){
	q = ""
	n_length = ""
	n_missing = "-n ${n_missing}"
}else{
	if(method=="length"){
		q = ""
		n_length = "-n ${n_length}"
		n_missing = ""
	}else{
		q = "-q ${q}"
		n_length = ""
		n_missing = ""
	}
}

readArray = reads.toString().split(' ')	

fasta = (fasta=="true") ? "--fasta" : ""

if(mate=="pair"){
	R1 = readArray[0]
	R2 = readArray[1]
	"""
	FilterSeq.py ${method} -s $R1 ${q} ${n_length} ${n_missing} --nproc ${nproc} --log FS_R1_${name}.log --failed ${fasta} >> out_${R1}_FS.log
	FilterSeq.py ${method} -s $R2 ${q} ${n_length} ${n_missing} --nproc ${nproc} --log FS_R2_${name}.log --failed ${fasta} >> out_${R1}_FS.log
	"""
}else{
	R1 = readArray[0]
	"""
	FilterSeq.py ${method} -s $R1 ${q} ${n_length} ${n_missing} --nproc ${nproc} --log FS_${name}.log --failed ${fasta} >> out_${R1}_FS.log
	"""
}


}


process parse_log_FSQ {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*table.tab$/) "reports/$filename"}
input:
 set val(name), file(log_file) from g_0_logFile1_g_12
 val mate from g_6_mate_g_12

output:
 set val(name), file("*table.tab")  into g_12_logFile00

script:
readArray = log_file.toString()

"""
ParseLog.py -l ${readArray}  -f ID QUALITY
"""

}


process filter_seq_length {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /FS_.*$/) "reports/$filename"}
input:
 set val(name),file(reads) from g_0_reads0_g_1
 val mate from g_6_mate_g_1

output:
 set val(name), file("*_${method}-pass.fast*")  into g_1_reads0_g_4, g_1_reads0_g_8
 set val(name), file("FS_*")  into g_1_logFile1_g_11
 set val(name), file("*_${method}-fail.fast*") optional true  into g_1_reads22
 set val(name),file("out*") optional true  into g_1_logFile33

script:
method = params.filter_seq_length.method
nproc = params.filter_seq_length.nproc
q = params.filter_seq_length.q
n_length = params.filter_seq_length.n_length
n_missing = params.filter_seq_length.n_missing
fasta = params.filter_seq_length.fasta
//* @style @condition:{method="quality",q}, {method="length",n_length}, {method="missing",n_missing} @multicolumn:{method,nproc}

if(method=="missing"){
	q = ""
	n_length = ""
	n_missing = "-n ${n_missing}"
}else{
	if(method=="length"){
		q = ""
		n_length = "-n ${n_length}"
		n_missing = ""
	}else{
		q = "-q ${q}"
		n_length = ""
		n_missing = ""
	}
}

readArray = reads.toString().split(' ')	

fasta = (fasta=="true") ? "--fasta" : ""

if(mate=="pair"){
	R1 = readArray[0]
	R2 = readArray[1]
	"""
	FilterSeq.py ${method} -s $R1 ${q} ${n_length} ${n_missing} --nproc ${nproc} --log FS_R1_${name}.log --failed ${fasta} >> out_${R1}_FS.log
	FilterSeq.py ${method} -s $R2 ${q} ${n_length} ${n_missing} --nproc ${nproc} --log FS_R2_${name}.log --failed ${fasta} >> out_${R1}_FS.log
	"""
}else{
	R1 = readArray[0]
	"""
	FilterSeq.py ${method} -s $R1 ${q} ${n_length} ${n_missing} --nproc ${nproc} --log FS_${name}.log --failed ${fasta} >> out_${R1}_FS.log
	"""
}


}


process parse_log_FSL {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tab$/) "reports/$filename"}
input:
 set val(name), file(log_file) from g_1_logFile1_g_11
 val mate from g_6_mate_g_11

output:
 set val(name), file("*.tab")  into g_11_logFile00

script:
readArray = log_file.toString()

"""
ParseLog.py -l ${readArray}  -f ID LENGTH
"""
}


process MaskPrimers_CPRIMERS {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_primers-fail.fastq$/) "failed_reads/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /MP_.*$/) "reports/$filename"}
input:
 val mate from g_6_mate_g_8
 set val(name),file(reads) from g_1_reads0_g_8

output:
 set val(name), file("*_primers-pass.fastq")  into g_8_reads0_g_9
 set val(name), file("*_primers-fail.fastq") optional true  into g_8_reads_failed11
 set val(name), file("MP_*")  into g_8_logFile2_g_17
 set val(name),file("out*")  into g_8_logFile33

script:
method = params.MaskPrimers_CPRIMERS.method
barcode_field = params.MaskPrimers_CPRIMERS.barcode_field
primer_field = params.MaskPrimers_CPRIMERS.primer_field
barcode = params.MaskPrimers_CPRIMERS.barcode
revpr = params.MaskPrimers_CPRIMERS.revpr
mode = params.MaskPrimers_CPRIMERS.mode
failed = params.MaskPrimers_CPRIMERS.failed
nproc = params.MaskPrimers_CPRIMERS.nproc
maxerror = params.MaskPrimers_CPRIMERS.maxerror
umi_length = params.MaskPrimers_CPRIMERS.umi_length
start = params.MaskPrimers_CPRIMERS.start
extract_length = params.MaskPrimers_CPRIMERS.extract_length
maxlen = params.MaskPrimers_CPRIMERS.maxlen
skiprc = params.MaskPrimers_CPRIMERS.skiprc
R1_primers = params.MaskPrimers_CPRIMERS.R1_primers
R2_primers = params.MaskPrimers_CPRIMERS.R2_primers
//* @style @condition:{method="score",umi_length,start,maxerror}{method="extract",umi_length,start},{method="align",maxerror,maxlen,skiprc}, {method="extract",start,extract_length} @array:{method,barcode_field,primer_field,barcode,revpr,mode,maxerror,umi_length,start,extract_length,maxlen,skiprc} @multicolumn:{method,barcode_field,primer_field,barcode,revpr,mode,failed,nproc,maxerror,umi_length,start,extract_length,maxlen,skiprc}

method = (method.collect().size==2) ? method : [method[0],method[0]]
barcode_field = (barcode_field.collect().size==2) ? barcode_field : [barcode_field[0],barcode_field[0]]
primer_field = (primer_field.collect().size==2) ? primer_field : [primer_field[0],primer_field[0]]
barcode = (barcode.collect().size==2) ? barcode : [barcode[0],barcode[0]]
revpr = (revpr.collect().size==2) ? revpr : [revpr[0],revpr[0]]
mode = (mode.collect().size==2) ? mode : [mode[0],mode[0]]
maxerror = (maxerror.collect().size==2) ? maxerror : [maxerror[0],maxerror[0]]
umi_length = (umi_length.collect().size==2) ? umi_length : [umi_length[0],umi_length[0]]
start = (start.collect().size==2) ? start : [start[0],start[0]]
extract_length = (extract_length.collect().size==2) ? extract_length : [extract_length[0],extract_length[0]]
maxlen = (maxlen.collect().size==2) ? maxlen : [maxlen[0],maxlen[0]]
skiprc = (skiprc.collect().size==2) ? skiprc : [skiprc[0],skiprc[0]]
failed = (failed=="true") ? "--failed" : ""

def args_values = [];
[method,barcode_field,primer_field,barcode,revpr,mode,maxerror,umi_length,start,extract_length,maxlen,skiprc].transpose().each { m,bf,pf,bc,rp,md,mr,ul,s,el,ml,sk -> {
    
    if(m=="align"){
        s = ""
    }else{
        if(bc=="false"){
            s = "--start ${s}"
        }else{
            s = s + ul
            s = "--start ${s}"
        }
    }
    
    el = (m=="extract") ? "--len ${el}" : ""
    mr = (m=="extract") ? "" : "--maxerror ${mr}" 
    ml = (m=="align") ? "--maxlen ${ml}" : "" 
    sk = (m=="align" && sk=="true") ? "--skiprc" : "" 
    
    PRIMER_FIELD = "${pf}"
    
    // all
    bf = (bf=="") ? "" : "--bf ${bf}"
    pf = (pf=="") ? "" : "--pf ${pf}"
    bc = (bc=="false") ? "" : "--barcode"
    rp = (rp=="false") ? "" : "--revpr"
    args_values.add("${m} --mode ${md} ${bf} ${pf} ${bc} ${rp} ${mr} ${s} ${el} ${ml} ${sk}")
    
    
}}

readArray = reads.toString().split(' ')
if(mate=="pair"){
	args_1 = args_values[0]
	args_2 = args_values[1]
	
  


	R1 = readArray[0]
	R2 = readArray[1]
	
	R1_primers = (method[0]=="extract") ? "" : "-p ${R1_primers}"
	R2_primers = (method[1]=="extract") ? "" : "-p ${R2_primers}"
	
	
	"""
	
	MaskPrimers.py ${args_1} -s ${R1} ${R1_primers} --log MP_R1_${name}.log  --nproc ${nproc} ${failed} >> out_${R1}_MP.log
	MaskPrimers.py ${args_2} -s ${R2} ${R2_primers} --log MP_R2_${name}.log  --nproc ${nproc} ${failed} >> out_${R1}_MP.log
	"""
}else{
	args_1 = args_values[0]
	
	R1_primers = (method[0]=="extract") ? "" : "-p ${R1_primers}"
	
	R1 = readArray[0]

	"""
	echo -e "Assuming inputs for R1\n"
	
	MaskPrimers.py ${args_1} -s ${reads} ${R1_primers} --log MP_${name}.log  --nproc ${nproc} ${failed} >> out_${R1}_MP.log
	"""
}

}


process parse_log_MP_CPRIMERS {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tab$/) "reports/$filename"}
input:
 val mate from g_6_mate_g_17
 set val(name), file(log_file) from g_8_logFile2_g_17

output:
 set val(name), file("*.tab")  into g_17_logFile00

script:
suffix = params.parse_log_MP_CPRIMERS.suffix
readArray = log_file.toString()	

outname = readArray - '.log' +  suffix
"""
ParseLog.py -l ${readArray} --outname ${outname} -f ID PRIMER BARCODE ERROR
"""


}


process MaskPrimers_VPRIMERS {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_primers-fail.fastq$/) "failed_reads/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /MP_.*$/) "reports/$filename"}
input:
 val mate from g_6_mate_g_9
 set val(name),file(reads) from g_8_reads0_g_9

output:
 set val(name), file("*_primers-pass.fastq")  into g_9_reads0_g_10
 set val(name), file("*_primers-fail.fastq") optional true  into g_9_reads_failed11
 set val(name), file("MP_*")  into g_9_logFile2_g_20
 set val(name),file("out*")  into g_9_logFile33

script:
method = params.MaskPrimers_VPRIMERS.method
barcode_field = params.MaskPrimers_VPRIMERS.barcode_field
primer_field = params.MaskPrimers_VPRIMERS.primer_field
barcode = params.MaskPrimers_VPRIMERS.barcode
revpr = params.MaskPrimers_VPRIMERS.revpr
mode = params.MaskPrimers_VPRIMERS.mode
failed = params.MaskPrimers_VPRIMERS.failed
nproc = params.MaskPrimers_VPRIMERS.nproc
maxerror = params.MaskPrimers_VPRIMERS.maxerror
umi_length = params.MaskPrimers_VPRIMERS.umi_length
start = params.MaskPrimers_VPRIMERS.start
extract_length = params.MaskPrimers_VPRIMERS.extract_length
maxlen = params.MaskPrimers_VPRIMERS.maxlen
skiprc = params.MaskPrimers_VPRIMERS.skiprc
R1_primers = params.MaskPrimers_VPRIMERS.R1_primers
R2_primers = params.MaskPrimers_VPRIMERS.R2_primers
//* @style @condition:{method="score",umi_length,start,maxerror}{method="extract",umi_length,start},{method="align",maxerror,maxlen,skiprc}, {method="extract",start,extract_length} @array:{method,barcode_field,primer_field,barcode,revpr,mode,maxerror,umi_length,start,extract_length,maxlen,skiprc} @multicolumn:{method,barcode_field,primer_field,barcode,revpr,mode,failed,nproc,maxerror,umi_length,start,extract_length,maxlen,skiprc}

method = (method.collect().size==2) ? method : [method[0],method[0]]
barcode_field = (barcode_field.collect().size==2) ? barcode_field : [barcode_field[0],barcode_field[0]]
primer_field = (primer_field.collect().size==2) ? primer_field : [primer_field[0],primer_field[0]]
barcode = (barcode.collect().size==2) ? barcode : [barcode[0],barcode[0]]
revpr = (revpr.collect().size==2) ? revpr : [revpr[0],revpr[0]]
mode = (mode.collect().size==2) ? mode : [mode[0],mode[0]]
maxerror = (maxerror.collect().size==2) ? maxerror : [maxerror[0],maxerror[0]]
umi_length = (umi_length.collect().size==2) ? umi_length : [umi_length[0],umi_length[0]]
start = (start.collect().size==2) ? start : [start[0],start[0]]
extract_length = (extract_length.collect().size==2) ? extract_length : [extract_length[0],extract_length[0]]
maxlen = (maxlen.collect().size==2) ? maxlen : [maxlen[0],maxlen[0]]
skiprc = (skiprc.collect().size==2) ? skiprc : [skiprc[0],skiprc[0]]
failed = (failed=="true") ? "--failed" : ""

def args_values = [];
[method,barcode_field,primer_field,barcode,revpr,mode,maxerror,umi_length,start,extract_length,maxlen,skiprc].transpose().each { m,bf,pf,bc,rp,md,mr,ul,s,el,ml,sk -> {
    
    if(m=="align"){
        s = ""
    }else{
        if(bc=="false"){
            s = "--start ${s}"
        }else{
            s = s + ul
            s = "--start ${s}"
        }
    }
    
    el = (m=="extract") ? "--len ${el}" : ""
    mr = (m=="extract") ? "" : "--maxerror ${mr}" 
    ml = (m=="align") ? "--maxlen ${ml}" : "" 
    sk = (m=="align" && sk=="true") ? "--skiprc" : "" 
    
    PRIMER_FIELD = "${pf}"
    
    // all
    bf = (bf=="") ? "" : "--bf ${bf}"
    pf = (pf=="") ? "" : "--pf ${pf}"
    bc = (bc=="false") ? "" : "--barcode"
    rp = (rp=="false") ? "" : "--revpr"
    args_values.add("${m} --mode ${md} ${bf} ${pf} ${bc} ${rp} ${mr} ${s} ${el} ${ml} ${sk}")
    
    
}}

readArray = reads.toString().split(' ')
if(mate=="pair"){
	args_1 = args_values[0]
	args_2 = args_values[1]
	
  


	R1 = readArray[0]
	R2 = readArray[1]
	
	R1_primers = (method[0]=="extract") ? "" : "-p ${R1_primers}"
	R2_primers = (method[1]=="extract") ? "" : "-p ${R2_primers}"
	
	
	"""
	
	MaskPrimers.py ${args_1} -s ${R1} ${R1_primers} --log MP_R1_${name}.log  --nproc ${nproc} ${failed} >> out_${R1}_MP.log
	MaskPrimers.py ${args_2} -s ${R2} ${R2_primers} --log MP_R2_${name}.log  --nproc ${nproc} ${failed} >> out_${R1}_MP.log
	"""
}else{
	args_1 = args_values[0]
	
	R1_primers = (method[0]=="extract") ? "" : "-p ${R1_primers}"
	
	R1 = readArray[0]

	"""
	echo -e "Assuming inputs for R1\n"
	
	MaskPrimers.py ${args_1} -s ${reads} ${R1_primers} --log MP_${name}.log  --nproc ${nproc} ${failed} >> out_${R1}_MP.log
	"""
}

}


process parse_log_MP_VPRIMERS {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tab$/) "reports/$filename"}
input:
 val mate from g_6_mate_g_20
 set val(name), file(log_file) from g_9_logFile2_g_20

output:
 set val(name), file("*.tab")  into g_20_logFile00

script:
suffix = params.parse_log_MP_VPRIMERS.suffix
readArray = log_file.toString()	

outname = readArray - '.log' +  suffix
"""
ParseLog.py -l ${readArray} --outname ${outname} -f ID PRIMER BARCODE ERROR
"""


}


process MaskPrimers_EXTRACT {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_primers-fail.fastq$/) "failed_reads/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /MP_.*$/) "reports/$filename"}
input:
 val mate from g_6_mate_g_10
 set val(name),file(reads) from g_9_reads0_g_10

output:
 set val(name), file("*_primers-pass.fastq")  into g_10_reads0_g_26
 set val(name), file("*_primers-fail.fastq") optional true  into g_10_reads_failed11
 set val(name), file("MP_*")  into g_10_logFile2_g_21
 set val(name),file("out*")  into g_10_logFile33

script:
method = params.MaskPrimers_EXTRACT.method
barcode_field = params.MaskPrimers_EXTRACT.barcode_field
primer_field = params.MaskPrimers_EXTRACT.primer_field
barcode = params.MaskPrimers_EXTRACT.barcode
revpr = params.MaskPrimers_EXTRACT.revpr
mode = params.MaskPrimers_EXTRACT.mode
failed = params.MaskPrimers_EXTRACT.failed
nproc = params.MaskPrimers_EXTRACT.nproc
maxerror = params.MaskPrimers_EXTRACT.maxerror
umi_length = params.MaskPrimers_EXTRACT.umi_length
start = params.MaskPrimers_EXTRACT.start
extract_length = params.MaskPrimers_EXTRACT.extract_length
maxlen = params.MaskPrimers_EXTRACT.maxlen
skiprc = params.MaskPrimers_EXTRACT.skiprc
R1_primers = params.MaskPrimers_EXTRACT.R1_primers
R2_primers = params.MaskPrimers_EXTRACT.R2_primers
//* @style @condition:{method="score",umi_length,start,maxerror}{method="extract",umi_length,start},{method="align",maxerror,maxlen,skiprc}, {method="extract",start,extract_length} @array:{method,barcode_field,primer_field,barcode,revpr,mode,maxerror,umi_length,start,extract_length,maxlen,skiprc} @multicolumn:{method,barcode_field,primer_field,barcode,revpr,mode,failed,nproc,maxerror,umi_length,start,extract_length,maxlen,skiprc}

method = (method.collect().size==2) ? method : [method[0],method[0]]
barcode_field = (barcode_field.collect().size==2) ? barcode_field : [barcode_field[0],barcode_field[0]]
primer_field = (primer_field.collect().size==2) ? primer_field : [primer_field[0],primer_field[0]]
barcode = (barcode.collect().size==2) ? barcode : [barcode[0],barcode[0]]
revpr = (revpr.collect().size==2) ? revpr : [revpr[0],revpr[0]]
mode = (mode.collect().size==2) ? mode : [mode[0],mode[0]]
maxerror = (maxerror.collect().size==2) ? maxerror : [maxerror[0],maxerror[0]]
umi_length = (umi_length.collect().size==2) ? umi_length : [umi_length[0],umi_length[0]]
start = (start.collect().size==2) ? start : [start[0],start[0]]
extract_length = (extract_length.collect().size==2) ? extract_length : [extract_length[0],extract_length[0]]
maxlen = (maxlen.collect().size==2) ? maxlen : [maxlen[0],maxlen[0]]
skiprc = (skiprc.collect().size==2) ? skiprc : [skiprc[0],skiprc[0]]
failed = (failed=="true") ? "--failed" : ""

def args_values = [];
[method,barcode_field,primer_field,barcode,revpr,mode,maxerror,umi_length,start,extract_length,maxlen,skiprc].transpose().each { m,bf,pf,bc,rp,md,mr,ul,s,el,ml,sk -> {
    
    if(m=="align"){
        s = ""
    }else{
        if(bc=="false"){
            s = "--start ${s}"
        }else{
            s = s + ul
            s = "--start ${s}"
        }
    }
    
    el = (m=="extract") ? "--len ${el}" : ""
    mr = (m=="extract") ? "" : "--maxerror ${mr}" 
    ml = (m=="align") ? "--maxlen ${ml}" : "" 
    sk = (m=="align" && sk=="true") ? "--skiprc" : "" 
    
    PRIMER_FIELD = "${pf}"
    
    // all
    bf = (bf=="") ? "" : "--bf ${bf}"
    pf = (pf=="") ? "" : "--pf ${pf}"
    bc = (bc=="false") ? "" : "--barcode"
    rp = (rp=="false") ? "" : "--revpr"
    args_values.add("${m} --mode ${md} ${bf} ${pf} ${bc} ${rp} ${mr} ${s} ${el} ${ml} ${sk}")
    
    
}}

readArray = reads.toString().split(' ')
if(mate=="pair"){
	args_1 = args_values[0]
	args_2 = args_values[1]
	
  


	R1 = readArray[0]
	R2 = readArray[1]
	
	R1_primers = (method[0]=="extract") ? "" : "-p ${R1_primers}"
	R2_primers = (method[1]=="extract") ? "" : "-p ${R2_primers}"
	
	
	"""
	
	MaskPrimers.py ${args_1} -s ${R1} ${R1_primers} --log MP_R1_${name}.log  --nproc ${nproc} ${failed} >> out_${R1}_MP.log
	MaskPrimers.py ${args_2} -s ${R2} ${R2_primers} --log MP_R2_${name}.log  --nproc ${nproc} ${failed} >> out_${R1}_MP.log
	"""
}else{
	args_1 = args_values[0]
	
	R1_primers = (method[0]=="extract") ? "" : "-p ${R1_primers}"
	
	R1 = readArray[0]

	"""
	echo -e "Assuming inputs for R1\n"
	
	MaskPrimers.py ${args_1} -s ${reads} ${R1_primers} --log MP_${name}.log  --nproc ${nproc} ${failed} >> out_${R1}_MP.log
	"""
}

}


process parse_log_MP_EXTRACT {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tab$/) "reports/$filename"}
input:
 val mate from g_6_mate_g_21
 set val(name), file(log_file) from g_10_logFile2_g_21

output:
 set val(name), file("*.tab")  into g_21_logFile00

script:
suffix = params.parse_log_MP_EXTRACT.suffix
readArray = log_file.toString()	

outname = readArray - '.log' +  suffix
"""
ParseLog.py -l ${readArray} --outname ${outname} -f ID PRIMER BARCODE ERROR
"""


}


process check_for_seqs {

input:
 set val(name), file(reads) from g_10_reads0_g_26

output:
 set val(name),file("*fastq")  into g_26_reads0_g_27

script:
primers_file = params.check_for_seqs.primers_file

outfq = reads.toString() - '.fastq' + '_cat.fastq'

"""
#!/usr/bin/env python3

import sys
import Bio
from Bio import SeqIO

primers = SeqIO.to_dict(SeqIO.parse("${primers_file}", "fasta"))

records = []

with open("${reads}") as handle:
    for record in SeqIO.parse(handle, "fastq"):
        add = True
        for primer in primers:
            if primers[primer].seq in record.seq:                
                add = False
        if add:
            records.append(record)

SeqIO.write(records, "${outfq}", "fastq")
"""

}


process align_sets {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /AS_.*$/) "reports/$filename"}
input:
 set val(name),file(reads) from g_26_reads0_g_27
 val mate from g_6_mate_g_27

output:
 set val(name),file("*_align-pass.fastq")  into g_27_reads0_g_31
 set val(name), file("AS_*")  into g_27_logFile1_g_28
 set val(name),file("*_align-fail.fastq") optional true  into g_27_reads_failed22
 set val(name), file("out*") optional true  into g_27_logFile33

script:
method = params.align_sets.method
bf = params.align_sets.bf
div = params.align_sets.div
failed = params.align_sets.failed
nproc = params.align_sets.nproc

muscle_exec = params.align_sets.muscle_exec

offset_table = params.align_sets.offset_table
pf = params.align_sets.pf
mode = params.align_sets.mode

primer_file = params.align_sets.primer_file
reverse = params.align_sets.reverse

//* @style @condition:{method="muscle",muscle_exec}, {method="offset",offset_table,pf,mode}, {method="table",muscle_exec,primer_file,reverse} @multicolumn:{method,bf,div,nproc},{offset,pf,mode}, {primer_file,reverse}


readArray = reads.toString().split(' ')	

reverse_arg = (reverse=="false") ? "" : "--reverse"
div_arg = (div=="false") ? "" : "--div"
failed_arg = (failed=="true") ? "--failed" : "" 
bf = "--bf ${bf}"

primer_file_argv = ""

if(method=="offset"){
	pf = "--pf ${pf}"
	mode = "--mode ${mode}"
	offset_table_argv = "-d ${offset_table}"
	muscle_exec_argv = ""
}else{
	pf = ""
	mode = ""
	offset_table_argv = ""
	muscle_exec_argv = "--exec ${muscle_exec}"
	
	if(method=="table"){
		primer_file_argv = "-p ${primer_file}"
	}
}

if(mate=="pair"){
	R1 = readArray[0]
	R2 = readArray[1]
	
	
	"""
	AlignSets.py ${method} -s ${R1} ${bf} ${muscle_exec_argv} ${div_arg} ${reverse_arg} ${failed_arg} ${pf} ${offset_table_argv} ${mode} ${primer_file_argv} --log AS_R1_${name}.log --nproc ${nproc} >> out_${R1}_AS.log
	AlignSets.py ${method} -s ${R2} ${bf} ${muscle_exec_argv} ${div_arg} ${reverse_arg} ${failed_arg} ${pf} ${offset_table_argv} ${mode} ${primer_file_argv} --log AS_R2_${name}.log --nproc ${nproc} >> out_${R1}_AS.log
	"""
	
}else{
	R1 = readArray[0]
	"""
	AlignSets.py ${method} -s ${R1} ${bf} ${muscle_exec_argv} ${div_arg} ${reverse_arg} ${failed_arg} ${pf} ${offset_table_argv} ${mode} ${primer_file_argv} --log AS_R1_${name}.log --nproc ${nproc} >> out_${R1}_AS.log
	"""
}

}


process cluster_sets {

input:
 set val(name),file(reads) from g_27_reads0_g_31
 val mate from g_6_mate_g_31

output:
 set val(name),file("*_cluster-pass.fastq")  into g_31_reads0_g_33
 set val(name),file("*_cluster-fail.fastq") optional true  into g_31_reads_failed11

script:
method = params.cluster_sets.method
failed = params.cluster_sets.failed
nproc = params.cluster_sets.nproc
cluster_field = params.cluster_sets.cluster_field
ident = params.cluster_sets.ident
length = params.cluster_sets.length
prefix = params.cluster_sets.prefix
cluster_tool = params.cluster_sets.cluster_tool
cluster_exec = params.cluster_sets.cluster_exec
set_field = params.cluster_sets.set_field
start = params.cluster_sets.start
end = params.cluster_sets.end
barcode_field = params.cluster_sets.barcode_field
//* @style @condition:{method="set",set_field,start,end},{method="all",start,end},{method="barcode",barcode_field} @array:{method,failed,cluster_field,ident,length,prefix,cluster_tool,cluster_exec,set_field,start,end,barcode_field}  @multicolumn:{method,failed,nproc,cluster_field,ident,length,prefix,cluster_tool,cluster_exec},{set_field,start,end,barcode_field}

method = (method.size==2) ? method : [method[0],method[0]]
failed = (failed.size==2) ? failed : [failed[0],failed[0]]
cluster_field = (cluster_field.size==2) ? cluster_field : [cluster_field[0],cluster_field[0]]
ident = (ident.size==2) ? ident : [ident[0],ident[0]]
length = (length.size==2) ? length : [length[0],length[0]]
prefix = (prefix.size==2) ? prefix : [prefix[0],prefix[0]]
cluster_tool = (cluster_tool.size==2) ? cluster_tool : [cluster_tool[0],cluster_tool[0]]
cluster_exec = (cluster_exec.size==2) ? cluster_exec : [cluster_exec[0],cluster_exec[0]]
set_field = (set_field.size==2) ? set_field : [set_field[0],set_field[0]]
start = (start.size==2) ? start : [start[0],start[0]]
end = (end.size==2) ? end : [end[0],end[0]]
barcode_field = (barcode_field.size==2) ? barcode_field : [barcode_field[0],barcode_field[0]]

def args_values = [];
[method, failed, cluster_field, ident, length, prefix, cluster_tool, cluster_exec, set_field, start, end, barcode_field].transpose().each { m, f, cf, i, l, p, ct, ce, sf, s, e, bf -> {
    f = (f=="true") ? "--failed" : ""
    p = (p=="") ? "" : "--prefix ${p}" 
    ce = (ce=="") ? "" : "--exec ${ce}" 
    sf = (m=="set") ? "-f ${sf}" : ""
    s = (m=="barcode") ? "" : "--start ${s}" 
    e = (m=="barcode") ? "" : (e=="") ? "" : "--end ${e}" 
    bf = (m=="barcode") ? "-f ${bf}" : ""
    args_values.add("${m} ${f} -k ${cf} --ident ${i} --length ${l} ${p} --cluster ${ct} ${ce} ${sf} ${s} ${e} ${bf}")
}}


if(mate=="pair"){
	// files
	readArray = reads.toString().split(' ')	
	R1 = readArray.grep(~/.*R1.*/)[0]
	R2 = readArray.grep(~/.*R2.*/)[0]
	
	args_1 = args_values[0]
	args_2 = args_values[1]
	
	"""
	ClusterSets.py ${args_1} -s $R1  --nproc ${nproc}
	ClusterSets.py ${args_2} -s $R2  --nproc ${nproc}
	"""
}else{
	args_1 = args_values[0]
	"""
	ClusterSets.py ${args_1} -s $reads --nproc ${nproc}
	"""
}


}


process parse_headers_copy {

input:
 set val(name), file(reads) from g_31_reads0_g_33

output:
 set val(name),file("*${out}")  into g_33_reads0_g_32

script:
method = params.parse_headers_copy.method
act = params.parse_headers_copy.act
args = params.parse_headers_copy.args


if(method=="collapse" || method=="copy" || method=="rename" || method=="merge"){
	out="_reheader.fastq"
	act = (act=="none") ? "" : "--act ${act}"
	"""
	ParseHeaders.py  ${method} -s ${reads} ${args} ${act}
	"""
}else{
	if(method=="table"){
			out=".tab"
			"""
			ParseHeaders.py ${method} -s ${reads} ${args}
			"""	
	}else{
		out="_reheader.fastq"
		"""
		ParseHeaders.py ${method} -s ${reads} ${args}
		"""		
	}
}


}

boolean isCollectionOrArray_bc(object) {    
    [Collection, Object[]].any { it.isAssignableFrom(object.getClass()) }
}

def args_creator_bc(barcode_field, primer_field, act, copy_field, mincount, minqual, minfreq, maxerror, prcons, maxgap, maxdiv, dep){
	def args_values;
    if(isCollectionOrArray_bc(barcode_field) || isCollectionOrArray_bc(primer_field) || isCollectionOrArray_bc(copy_field) || isCollectionOrArray_bc(mincount) || isCollectionOrArray_bc(minqual) || isCollectionOrArray_bc(minfreq) || isCollectionOrArray_bc(maxerror) || isCollectionOrArray_bc(prcons) || isCollectionOrArray_bc(maxgap) || isCollectionOrArray_bc(maxdiv) || isCollectionOrArray_bc(dep)){
    	primer_field = (isCollectionOrArray_bc(primer_field)) ? primer_field : [primer_field,primer_field]
    	act = (isCollectionOrArray_bc(act)) ? act : [act,act]
    	copy_field = (isCollectionOrArray_bc(copy_field)) ? copy_field : [copy_field,copy_field]
    	mincount = (isCollectionOrArray_bc(mincount)) ? mincount : [mincount,mincount]
    	minqual = (isCollectionOrArray_bc(minqual)) ? minqual : [minqual,minqual]
    	minfreq = (isCollectionOrArray_bc(minfreq)) ? minfreq : [minfreq,minfreq]
    	maxerror = (isCollectionOrArray_bc(maxerror)) ? maxerror : [maxerror,maxerror]
    	prcons = (isCollectionOrArray_bc(prcons)) ? prcons : [prcons,prcons]
    	maxgap = (isCollectionOrArray_bc(maxgap)) ? maxgap : [maxgap,maxgap]
    	maxdiv = (isCollectionOrArray_bc(maxdiv)) ? maxdiv : [maxdiv,maxdiv]
    	dep = (isCollectionOrArray_bc(dep)) ? dep : [dep,dep]
    	args_values = []
        [barcode_field,primer_field,act,copy_field,mincount,minqual,minfreq,maxerror,prcons,maxgap,maxdiv,dep].transpose().each { bf,pf,a,cf,mc,mq,mf,mr,pc,mg,md,d -> {
            bf = (bf=="") ? "" : "--bf ${bf}"
            pf = (pf=="") ? "" : "--pf ${pf}" 
            a = (a=="none") ? "" : "--act ${a}" 
            cf = (cf=="") ? "" : "--cf ${cf}" 
            mr = (mr=="none") ? "" : "--maxerror ${mr}" 
            pc = (pc=="none") ? "" : "--prcons ${pc}" 
            mg = (mg=="none") ? "" : "--maxgap ${mg}" 
            md = (md=="none") ? "" : "--maxdiv ${md}" 
            d = (d=="true") ? "--dep" : "" 
            args_values.add("${bf} ${pf} ${a} ${cf} -n ${mc} -q ${mq} --freq ${mf} ${mr} ${pc} ${mg} ${md} ${d}")
        }}
    }else{
        barcode_field = (barcode_field=="") ? "" : "--bf ${barcode_field}"
        primer_field = (primer_field=="") ? "" : "--pf ${primer_field}" 
        act = (act=="none") ? "" : "--act ${act}" 
        copy_field = (copy_field=="") ? "" : "--cf ${copy_field}" 
        maxerror = (maxerror=="none") ? "" : "--maxerror ${maxerror}" 
        prcons = (prcons=="none") ? "" : "--prcons ${prcons}" 
        maxgap = (maxgap=="none") ? "" : "--maxgap ${maxgap}" 
        maxdiv = (maxdiv=="none") ? "" : "--maxdiv ${maxdiv}" 
        dep = (dep=="true") ? "--dep" : "" 
        args_values = "${barcode_field} ${primer_field} ${act} ${copy_field} -n ${mincount} -q ${minqual} --freq ${minfreq} ${maxerror} ${prcons} ${maxgap} ${maxdiv} ${dep}"
    }
    return args_values
}


process build_consensus {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /BC.*$/) "reports/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_consensus-fail.fastq$/) "failed_reads/$filename"}
input:
 set val(name),file(reads) from g_33_reads0_g_32
 val mate from g_6_mate_g_32

output:
 set val(name),file("*_consensus-pass.fastq")  into g_32_reads0_g_34
 file "BC*"  into g_32_logFile1_g_36
 set val(name),file("*_consensus-fail.fastq") optional true  into g_32_reads22

script:
failed = params.build_consensus.failed
nproc = params.build_consensus.nproc
barcode_field = params.build_consensus.barcode_field
primer_field = params.build_consensus.primer_field
act = params.build_consensus.act
copy_field = params.build_consensus.copy_field
mincount = params.build_consensus.mincount
minqual = params.build_consensus.minqual
minfreq = params.build_consensus.minfreq
maxerror = params.build_consensus.maxerror
prcons = params.build_consensus.prcons
maxgap = params.build_consensus.maxgap
maxdiv = params.build_consensus.maxdiv
dep = params.build_consensus.dep
//* @style @condition:{act="none",},{act="min",copy_field},{act="max",copy_field},{act="sum",copy_field},{act="set",copy_field},{act="majority",copy_field} @array:{barcode_field,primer_field,act,copy_field,mincount,minqual,minfreq,maxerror,prcons,maxgap,maxdiv,dep} @multicolumn:{failed,nproc},{barcode_field,primer_field,act,copy_field}, {mincount,minqual,minfreq,maxerror,prcons,maxgap,maxdiv,dep}

args_values_bc = args_creator_bc(barcode_field, primer_field, act, copy_field, mincount, minqual, minfreq, maxerror, prcons, maxgap, maxdiv, dep)

// args 
if(isCollectionOrArray_bc(args_values_bc)){
	args_1 = args_values_bc[0]
	args_2 = args_values_bc[1]
}else{
	args_1 = args_values_bc
	args_2 = args_values_bc
}

failed = (failed=="true") ? "--failed" : "" 


if(mate=="pair"){
	// files
	readArray = reads.toString().split(' ')	
	R1 = readArray.grep(~/.*R1.*/)[0]
	R2 = readArray.grep(~/.*R2.*/)[0]
	
	"""
	BuildConsensus.py -s $R1 ${args_1} --outname ${name}_R1 --log BC_${name}_R1.log ${failed} --nproc ${nproc}
	BuildConsensus.py -s $R2 ${args_2} --outname ${name}_R2 --log BC_${name}_R2.log ${failed} --nproc ${nproc}
	"""
}else{
	"""
	BuildConsensus.py -s $reads ${args_1} --outname ${name} --log BC_${name}.log ${failed} --nproc ${nproc}
	"""
}


}


process parse_log_BC {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tab$/) "reports/$filename"}
input:
 file log_file from g_32_logFile1_g_36
 val mate from g_6_mate_g_36

output:
 file "*.tab"  into g_36_logFile00

script:
readArray = log_file.toString()

"""
ParseLog.py -l ${readArray} -f BARCODE SEQCOUNT CONSCOUNT PRCONS PRFREQ ERROR
"""

}


process parse_log_AS {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tab$/) "reports/$filename"}
input:
 set val(name), file(log_file) from g_27_logFile1_g_28
 val mate from g_6_mate_g_28

output:
 set val(name), file("*.tab")  into g_28_logFile00

script:
field_to_parse = params.parse_log_AS.field_to_parse
readArray = log_file.toString()	

"""
ParseLog.py -l ${readArray}  -f ${field_to_parse}
"""

}

//* params.run_FastQC =  "no"  //* @dropdown @options:"yes","no"
if (params.run_FastQC == "no") { println "INFO: FastQC will be skipped"}


process FastQC {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.(html|zip)$/) "reports/$filename"}
input:
 set val(name), file(reads) from g_1_reads0_g_4
 val mate from g_6_mate_g_4

output:
 file '*.{html,zip}'  into g_4_FastQCout00

errorStrategy 'retry'
maxRetries 5

script:
nameAll = reads.toString()
if (nameAll.contains('.gz')) {
    file =  nameAll - '.gz' - '.gz'
    runGzip = "ls *.gz | xargs -i echo gzip -df {} | sh"
} else {
    file =  nameAll 
    runGzip = ''
}
"""
if [ "${params.run_FastQC}" == "yes" ]; then
    ${runGzip}
    fastqc ${file} 
else
    touch process.skiped.html
fi
"""
}


process parse_headers_collapse {

input:
 set val(name), file(reads) from g_32_reads0_g_34

output:
 set val(name),file("*${out}")  into g_34_reads0_g_35

script:
method = params.parse_headers_collapse.method
act = params.parse_headers_collapse.act
args = params.parse_headers_collapse.args


if(method=="collapse" || method=="copy" || method=="rename" || method=="merge"){
	out="_reheader.fastq"
	act = (act=="none") ? "" : "--act ${act}"
	"""
	ParseHeaders.py  ${method} -s ${reads} ${args} ${act}
	"""
}else{
	if(method=="table"){
			out=".tab"
			"""
			ParseHeaders.py ${method} -s ${reads} ${args}
			"""	
	}else{
		out="_reheader.fastq"
		"""
		ParseHeaders.py ${method} -s ${reads} ${args}
		"""		
	}
}


}


process collapse_seq {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_collapse-unique.fast.*$/) "reads/$filename"}
input:
 set val(name), file(reads) from g_34_reads0_g_35

output:
 set val(name),  file("*_collapse-unique.fast*")  into g_35_reads0_g_39
 set val(name),  file("*_collapse-duplicate.fast*") optional true  into g_35_reads_duplicate11
 set val(name),  file("*_collapse-undetermined.fast*") optional true  into g_35_reads_undetermined22
 file "CS_*"  into g_35_logFile33

script:
max_missing = params.collapse_seq.max_missing
inner = params.collapse_seq.inner
fasta = params.collapse_seq.fasta
act = params.collapse_seq.act
uf = params.collapse_seq.uf
cf = params.collapse_seq.cf
nproc = params.collapse_seq.nproc
failed = params.collapse_seq.failed

inner = (inner=="true") ? "--inner" : ""
fasta = (fasta=="true") ? "--fasta" : ""
act = (act=="none") ? "" : "--act ${act}"
cf = (cf=="") ? "" : "--cf ${cf}"
uf = (uf=="") ? "" : "--uf ${uf}"
failed = (failed=="false") ? "" : "--failed"

"""
CollapseSeq.py -s ${reads} -n ${max_missing} ${fasta} ${inner} ${uf} ${cf} ${act} --log CS_${name}.log ${failed}
"""

}


process split_seq {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_atleast-.*.fast.*$/) "reads/$filename"}
input:
 set val(name),file(reads) from g_35_reads0_g_39

output:
 set val(name), file("*_atleast-*.fast*")  into g_39_fastaFile0_g_40, g_39_fastaFile0_g48_12, g_39_fastaFile0_g48_9
 set val(name),file("out*") optional true  into g_39_logFile11

script:
field = params.split_seq.field
num = params.split_seq.num
fasta = params.split_seq.fasta

readArray = reads.toString()

if(num!=0){
	num = " --num ${num}"
}else{
	num = ""
}

fasta = (fasta=="false") ? "" : "--fasta"

"""
SplitSeq.py group -s ${readArray} -f ${field} ${num} ${fasta} >> out_${readArray}_SS.log
"""

}


process Alignment_FLAIRSEQ_D_MakeBlastDb {

input:
 set val(db_name), file(germlineFile) from g_65_germlineFastaFile_g48_16

output:
 file "${db_name}"  into g48_16_germlineDb0_g48_9

script:

if(germlineFile.getName().endsWith("fasta")){
	"""
	sed -e '/^>/! s/[.]//g' ${germlineFile} > tmp_germline.fasta
	mkdir -m777 ${db_name}
	makeblastdb -parse_seqids -dbtype nucl -in tmp_germline.fasta -out ${db_name}/${db_name}
	"""
}else{
	"""
	echo something if off
	"""
}

}


process Alignment_FLAIRSEQ_J_MakeBlastDb {

input:
 set val(db_name), file(germlineFile) from g_66_germlineFastaFile_g48_17

output:
 file "${db_name}"  into g48_17_germlineDb0_g48_9

script:

if(germlineFile.getName().endsWith("fasta")){
	"""
	sed -e '/^>/! s/[.]//g' ${germlineFile} > tmp_germline.fasta
	mkdir -m777 ${db_name}
	makeblastdb -parse_seqids -dbtype nucl -in tmp_germline.fasta -out ${db_name}/${db_name}
	"""
}else{
	"""
	echo something if off
	"""
}

}


process Alignment_FLAIRSEQ_V_MakeBlastDb {

input:
 set val(db_name), file(germlineFile) from g_67_germlineFastaFile_g48_22

output:
 file "${db_name}"  into g48_22_germlineDb0_g48_9

script:

if(germlineFile.getName().endsWith("fasta")){
	"""
	sed -e '/^>/! s/[.]//g' ${germlineFile} > tmp_germline.fasta
	mkdir -m777 ${db_name}
	makeblastdb -parse_seqids -dbtype nucl -in tmp_germline.fasta -out ${db_name}/${db_name}
	"""
}else{
	"""
	echo something if off
	"""
}

}


process Alignment_FLAIRSEQ_C_MakeBlastDb {

input:

output:
 file "${db_name}"  into g48_54_germlineDb0_g48_9

script:

if(germlineFile.getName().endsWith("fasta")){
	"""
	sed -e '/^>/! s/[.]//g' ${germlineFile} > tmp_germline.fasta
	mkdir -m777 ${db_name}
	makeblastdb -parse_seqids -dbtype nucl -in tmp_germline.fasta -out ${db_name}/${db_name}
	"""
}else{
	"""
	echo something if off
	"""
}

}

g48_54_germlineDb0_g48_9= g48_54_germlineDb0_g48_9.ifEmpty([""]) 


process Alignment_FLAIRSEQ_IgBlastn {

input:
 set val(name),file(fastaFile) from g_39_fastaFile0_g48_9
 file db_v from g48_22_germlineDb0_g48_9
 file db_d from g48_16_germlineDb0_g48_9
 file db_j from g48_17_germlineDb0_g48_9
 file custom_internal_data from g_70_outputFileTxt_g48_9
 file db_c from g48_54_germlineDb0_g48_9

output:
 set val(name), file("${outfile}") optional true  into g48_9_igblastOut0_g48_12

script:
num_threads = params.Alignment_FLAIRSEQ_IgBlastn.num_threads
ig_seqtype = params.Alignment_FLAIRSEQ_IgBlastn.ig_seqtype
outfmt = params.Alignment_FLAIRSEQ_IgBlastn.outfmt
num_alignments_V = params.Alignment_FLAIRSEQ_IgBlastn.num_alignments_V
domain_system = params.Alignment_FLAIRSEQ_IgBlastn.domain_system

randomString = org.apache.commons.lang.RandomStringUtils.random(9, true, true)
outname = name + "_" + randomString
outfile = (outfmt=="MakeDb") ? name+"_"+randomString+".out" : name+"_"+randomString+".tsv"
outfmt = (outfmt=="MakeDb") ? "'7 std qseq sseq btop'" : "19"

if(db_v.toString()!="" && db_d.toString()!="" && db_j.toString()!=""){
	"""
	export IGDATA=/usr/local/share/igblast
	
	igblastn -query ${fastaFile} \
		-germline_db_V ${db_v}/${db_v} \
		-germline_db_D ${db_d}/${db_d} \
		-germline_db_J ${db_j}/${db_j} \
		-c_region_db ${db_c}/${db_c} \
		-num_alignments_V ${num_alignments_V} \
		-domain_system imgt \
		-auxiliary_data ${auxiliary_data} \
		-custom_internal_data ${custom_internal_data} \
		-outfmt ${outfmt} \
		-num_threads ${num_threads} \
		-out ${outfile}
	"""
}else{
	"""
	"""
}

}

g_68_mate_g48_12= g_68_mate_g48_12.ifEmpty("") 


process Alignment_FLAIRSEQ_MakeDb {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_db-pass.tsv$/) "alignment/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_db-fail.tsv$/) "alignment/$filename"}
input:
 set val(name),file(fastaFile) from g_39_fastaFile0_g48_12
 set val(name_igblast),file(igblastOut) from g48_9_igblastOut0_g48_12
 set val(name1), file(v_germline_file) from g_67_germlineFastaFile_g48_12
 set val(name2), file(d_germline_file) from g_65_germlineFastaFile_g48_12
 set val(name3), file(j_germline_file) from g_66_germlineFastaFile_g48_12
 val name_alignment from g_68_mate_g48_12

output:
 set val(name_igblast),file("*_db-pass.tsv") optional true  into g48_12_outputFileTSV0_g48_19
 set val("reference_set"), file("${reference_set}") optional true  into g48_12_germlineFastaFile11
 set val(name_igblast),file("*_db-fail.tsv") optional true  into g48_12_outputFileTSV22

script:

failed = params.Alignment_FLAIRSEQ_MakeDb.failed
format = params.Alignment_FLAIRSEQ_MakeDb.format
regions = params.Alignment_FLAIRSEQ_MakeDb.regions
extended = params.Alignment_FLAIRSEQ_MakeDb.extended
asisid = params.Alignment_FLAIRSEQ_MakeDb.asisid
asiscalls = params.Alignment_FLAIRSEQ_MakeDb.asiscalls
inferjunction = params.Alignment_FLAIRSEQ_MakeDb.inferjunction
partial = params.Alignment_FLAIRSEQ_MakeDb.partial
name_alignment = params.Alignment_FLAIRSEQ_MakeDb.name_alignment

failed = (failed=="true") ? "--failed" : ""
format = (format=="changeo") ? "--format changeo" : ""
extended = (extended=="true") ? "--extended" : ""
regions = (regions=="rhesus-igl") ? "--regions rhesus-igl" : ""
asisid = (asisid=="true") ? "--asis-id" : ""
asiscalls = (asiscalls=="true") ? "--asis-calls" : ""
inferjunction = (inferjunction=="true") ? "--infer-junction" : ""
partial = (partial=="true") ? "--partial" : ""

reference_set = "reference_set_makedb_"+name_alignment+".fasta"

outname = name_igblast+'_'+name_alignment

if(igblastOut.getName().endsWith(".out")){
	"""
	
	cat ${v_germline_file} ${d_germline_file} ${j_germline_file} ${c_germline_file} > ${reference_set}
	
	MakeDb.py igblast \
		-s ${fastaFile} \
		-i ${igblastOut} \
		-r ${v_germline_file} ${d_germline_file} ${j_germline_file} ${c_germline_file} \
		--log MD_${name}.log \
		--outname ${outname}\
		${extended} \
		${failed} \
		${format} \
		${regions} \
		${asisid} \
		${asiscalls} \
		${inferjunction} \
		${partial}
	"""
}else{
	"""
	
	"""
}

}

g_68_mate_g48_19= g_68_mate_g48_19.ifEmpty("") 


process Alignment_FLAIRSEQ_Collapse_AIRRseq {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${outfile}+passed.tsv$/) "collapsed/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${outfile}+failed.*$/) "collapsed/$filename"}
input:
 set val(name),file(airrFile) from g48_12_outputFileTSV0_g48_19
 val name_alignment from g_68_mate_g48_19

output:
 set val(name), file("${outfile}"+"passed.tsv") optional true  into g48_19_outputFileTSV0_g57_0
 set val(name), file("${outfile}"+"failed*") optional true  into g48_19_outputFileTSV11

script:
conscount_min = params.Alignment_FLAIRSEQ_Collapse_AIRRseq.conscount_min
n_max = params.Alignment_FLAIRSEQ_Collapse_AIRRseq.n_max
name_alignment = params.Alignment_FLAIRSEQ_Collapse_AIRRseq.name_alignment


outfile = airrFile.toString() - '.tsv' + name_alignment + "_collapsed-"

if(airrFile.getName().endsWith(".tsv")){	
	"""
	#!/usr/bin/env python3
	
	from pprint import pprint
	from collections import OrderedDict,Counter
	import itertools as it
	import datetime
	import pandas as pd
	import glob, os
	import numpy as np
	import re
	
	# column types default
	
	# dtype_default={'junction_length': 'Int64', 'np1_length': 'Int64', 'np2_length': 'Int64', 'v_sequence_start': 'Int64', 'v_sequence_end': 'Int64', 'v_germline_start': 'Int64', 'v_germline_end': 'Int64', 'd_sequence_start': 'Int64', 'd_sequence_end': 'Int64', 'd_germline_start': 'Int64', 'd_germline_end': 'Int64', 'j_sequence_start': 'Int64', 'j_sequence_end': 'Int64', 'j_germline_start': 'Int64', 'j_germline_end': 'Int64', 'v_score': 'Int64', 'v_identity': 'Int64', 'v_support': 'Int64', 'd_score': 'Int64', 'd_identity': 'Int64', 'd_support': 'Int64', 'j_score': 'Int64', 'j_identity': 'Int64', 'j_support': 'Int64'}
	
	SPLITSIZE=2
	
	class CollapseDict:
	    def __init__(self,iterable=(),_depth=0,
	                 nlim=10,conscount_flag=False):
	        self.lowqual={}
	        self.seqs = {}
	        self.children = {}
	        self.depth=_depth
	        self.nlim=nlim
	        self.conscount=conscount_flag
	        for fseq in iterable:
	            self.add(fseq)
	
	    def split(self):
	        newseqs = {}
	        for seq in self.seqs:
	            if len(seq)==self.depth:
	                newseqs[seq]=self.seqs[seq]
	            else:
	                if seq[self.depth] not in self.children:
	                    self.children[seq[self.depth]] = CollapseDict(_depth=self.depth+1)
	                self.children[seq[self.depth]].add(self.seqs[seq],seq)
	        self.seqs=newseqs
	
	    def add(self,fseq,key=None):
	        #if 'duplicate_count' not in fseq: fseq['duplicate_count']='1'
	        if 'KEY' not in fseq:
	            fseq['KEY']=fseq['sequence_vdj'].replace('-','').replace('.','')
	        if 'ISOTYPECOUNTER' not in fseq:
	            fseq['ISOTYPECOUNTER']=Counter([fseq['c_call']])
	        if 'VGENECOUNTER' not in fseq:
	            fseq['VGENECOUNTER']=Counter([fseq['v_call']])
	        if 'JGENECOUNTER' not in fseq:
	            fseq['JGENECOUNTER']=Counter([fseq['j_call']])
	        if key is None:
	            key=fseq['KEY']
	        if self.depth==0:
	            if (not fseq['j_call'] or not fseq['v_call']):
	                return
	            if fseq['sequence_vdj'].count('N')>self.nlim:
	                if key in self.lowqual:
	                    self.lowqual[key] = combine(self.lowqual[key],fseq,self.conscount)
	                else:
	                    self.lowqual[key] = fseq
	                return
	        if len(self.seqs)>SPLITSIZE:
	            self.split()
	        if key in self.seqs:
	            self.seqs[key] = combine(self.seqs[key],fseq,self.conscount)
	        elif (self.children is not None and
	              len(key)>self.depth and
	              key[self.depth] in self.children):
	            self.children[key[self.depth]].add(fseq,key)
	        else:
	            self.seqs[key] = fseq
	
	    def __iter__(self):
	        yield from self.seqs.items()
	        for d in self.children.values():
	            yield from d
	        yield from self.lowqual.items()
	
	    def neighbors(self,seq):
	        def nfil(x): return similar(seq,x)
	        yield from filter(nfil,self.seqs)
	        if len(seq)>self.depth:
	            for d in [self.children[c]
	                      for c in self.children
	                      if c=='N' or seq[self.depth]=='N' or c==seq[self.depth]]:
	                yield from d.neighbors(seq)
	
	    def fixedseqs(self):
	        return self
	        ncd = CollapseDict()
	        for seq,fseq in self:
	            newseq=seq
	            if 'N' in seq:
	                newseq=consensus(seq,self.neighbors(seq))
	                fseq['KEY']=newseq
	            ncd.add(fseq,newseq)
	        return ncd
	
	
	    def __len__(self):
	        return len(self.seqs)+sum(len(c) for c in self.children.values())+len(self.lowqual)
	
	def combine(f1,f2, conscount_flag):
	    def val(f): return -f['KEY'].count('N'),(int(f['consensus_count']) if 'consensus_count' in f else 0)
	    targ = (f1 if val(f1) >= val(f2) else f2).copy()
	    if conscount_flag:
	        targ['consensus_count'] =  int(f1['consensus_count'])+int(f2['consensus_count'])
	    targ['duplicate_count'] =  int(f1['duplicate_count'])+int(f2['duplicate_count'])
	    targ['ISOTYPECOUNTER'] = f1['ISOTYPECOUNTER']+f2['ISOTYPECOUNTER']
	    targ['VGENECOUNTER'] = f1['VGENECOUNTER']+f2['VGENECOUNTER']
	    targ['JGENECOUNTER'] = f1['JGENECOUNTER']+f2['JGENECOUNTER']
	    return targ
	
	def similar(s1,s2):
	    return len(s1)==len(s2) and all((n1==n2 or n1=='N' or n2=='N')
	                                  for n1,n2 in zip(s1,s2))
	
	def basecon(bases):
	    bases.discard('N')
	    if len(bases)==1: return bases.pop()
	    else: return 'N'
	
	def consensus(seq,A):
	    return ''.join((basecon(set(B)) if s=='N' else s) for (s,B) in zip(seq,zip(*A)))
	
	n_lim = int('${n_max}')
	conscount_filter = int('${conscount_min}')
	
	df = pd.read_csv('${airrFile}', sep = '\t') #, dtype = dtype_default)
	
	# make sure that all columns are int64 for createGermline
	idx_col = df.columns.get_loc("cdr3")
	cols =  [col for col in df.iloc[:,0:idx_col].select_dtypes('float64').columns.values.tolist() if not re.search('support|score|identity|freq', col)]
	df[cols] = df[cols].apply(lambda x: pd.to_numeric(x.replace("<NA>",np.NaN), errors = "coerce").astype("Int64"))
	
	conscount_flag = False
	if 'consensus_count' in df: conscount_flag = True
	if not 'duplicate_count' in df:
	    df['duplicate_count'] = 1
	if not 'c_call' in df or not 'isotype' in df or not 'prcons' in df or not 'primer' in df or not 'reverse_primer' in df:
	    if 'c_call' in df:
	        df['c_call'] = df['c_call']
	    elif 'isotype' in df:
	        df['c_call'] = df['isotype']
	    elif 'primer' in df:
	        df['c_call'] = df['primer']
	    elif 'reverse_primer' in df:
	        df['c_call'] = df['reverse_primer']    
	    elif 'prcons' in df:
	        df['c_call'] = df['prcons']
	    elif 'barcode' in df:
	        df['c_call'] = df['barcode']
	    else:
	        df['c_call'] = 'Ig'
	
	# removing sequenes with duplicated sequence id    
	dup_n = df[df.columns[0]].count()
	df = df.drop_duplicates(subset='sequence_id', keep='first')
	dup_n = str(dup_n - df[df.columns[0]].count())
	df['c_call'] = df['c_call'].astype('str').replace('<NA>','Ig')
	#df['consensus_count'].fillna(2, inplace=True)
	nrow_i = df[df.columns[0]].count()
	df = df[df.apply(lambda x: x['sequence_alignment'][0:(x['v_germline_end']-1)].count('N')<=n_lim, axis = 1)]
	low_n = str(nrow_i-df[df.columns[0]].count())
	
	df['sequence_vdj'] = df.apply(lambda x: x['sequence_alignment'].replace('-','').replace('.',''), axis = 1)
	header=list(df.columns)
	fasta_ = df.to_dict(orient='records')
	c = CollapseDict(fasta_,nlim=10)
	d=c.fixedseqs()
	header.append('ISOTYPECOUNTER')
	header.append('VGENECOUNTER')
	header.append('JGENECOUNTER')
	
	rec_list = []
	for i, f in enumerate(d):
	    rec=f[1]
	    rec['sequence']=rec['KEY']
	    rec['consensus_count']=int(rec['consensus_count']) if conscount_flag else None
	    rec['duplicate_count']=int(rec['duplicate_count'])
	    rec_list.append(rec)
	df2 = pd.DataFrame(rec_list, columns = header)        
	
	df2 = df2.drop('sequence_vdj', axis=1)
	
	collapse_n = str(df[df.columns[0]].count()-df2[df2.columns[0]].count())

	# removing sequences without J assignment and non functional
	nrow_i = df2[df2.columns[0]].count()
	cond = (~df2['j_call'].str.contains('J')|df2['productive'].isin(['F','FALSE','False']))
	df_non = df2[cond]
	
	
	df2 = df2[df2['productive'].isin(['T','TRUE','True'])]
	cond = ~(df2['j_call'].str.contains('J'))
	df2 = df2.drop(df2[cond].index.values)
	
	non_n = nrow_i-df2[df2.columns[0]].count()
	#if conscount_flag:
	#   df2['consensus_count'] = df2['consensus_count'].replace(1,2)
	
	# removing sequences with low cons count
	
	filter_column = "duplicate_count"
	if conscount_flag: filter_column = "consensus_count"
	df_cons_low = df2[df2[filter_column]<conscount_filter]
	nrow_i = df2[df2.columns[0]].count()
	df2 = df2[df2[filter_column]>=conscount_filter]
	
	
	cons_n = str(nrow_i-df2[df2.columns[0]].count())
	nrow_i = df2[df2.columns[0]].count()    
	
	df2.to_csv('${outfile}'+'passed.tsv', sep = '\t',index=False) #, compression='gzip'
	
	pd.concat([df_cons_low,df_non]).to_csv('${outfile}'+'failed.tsv', sep = '\t',index=False)
	
	print(str(low_n)+' Sequences had N count over 10')
	print(str(dup_n)+' Sequences had a duplicated sequnece id')
	print(str(collapse_n)+' Sequences were collapsed')
	print(str(df_non[df_non.columns[0]].count())+' Sequences were declared non functional or lacked a J assignment')
	#print(str(df_cons_low[df_cons_low.columns[0]].count())+' Sequences had a '+filter_column+' lower than threshold')
	print('Going forward with '+str(df2[df2.columns[0]].count())+' sequences')
	
	"""
}else{
	"""
	
	"""
}

}

g_67_germlineFastaFile_g57_0= g_67_germlineFastaFile_g57_0.ifEmpty([""]) 
g_65_germlineFastaFile_g57_0= g_65_germlineFastaFile_g57_0.ifEmpty([""]) 
g_66_germlineFastaFile_g57_0= g_66_germlineFastaFile_g57_0.ifEmpty([""]) 


process Clone_AIRRseq_First_CreateGermlines {

input:
 set val(name),file(airrFile) from g48_19_outputFileTSV0_g57_0
 set val(name1), file(v_germline_file) from g_67_germlineFastaFile_g57_0
 set val(name2), file(d_germline_file) from g_65_germlineFastaFile_g57_0
 set val(name3), file(j_germline_file) from g_66_germlineFastaFile_g57_0

output:
 set val(name),file("*_germ-pass.tsv")  into g57_0_outputFileTSV0_g57_2

script:
failed = params.Clone_AIRRseq_First_CreateGermlines.failed
format = params.Clone_AIRRseq_First_CreateGermlines.format
g = params.Clone_AIRRseq_First_CreateGermlines.g
cloned = params.Clone_AIRRseq_First_CreateGermlines.cloned
seq_field = params.Clone_AIRRseq_First_CreateGermlines.seq_field
v_field = params.Clone_AIRRseq_First_CreateGermlines.v_field
d_field = params.Clone_AIRRseq_First_CreateGermlines.d_field
j_field = params.Clone_AIRRseq_First_CreateGermlines.j_field
clone_field = params.Clone_AIRRseq_First_CreateGermlines.clone_field


failed = (failed=="true") ? "--failed" : ""
format = (format=="airr") ? "": "--format changeo"
g = "-g ${g}"
cloned = (cloned=="false") ? "" : "--cloned"

v_field = (v_field=="") ? "" : "--vf ${v_field}"
d_field = (d_field=="") ? "" : "--df ${d_field}"
j_field = (j_field=="") ? "" : "--jf ${j_field}"
seq_field = (seq_field=="") ? "" : "--sf ${seq_field}"

"""
CreateGermlines.py \
	-d ${airrFile} \
	-r ${v_germline_file} ${d_germline_file} ${j_germline_file} \
	${failed} \
	${format} \
	${g} \
	${cloned} \
	${v_field} \
	${d_field} \
	${j_field} \
	${seq_field} \
	${clone_field} \
	--log CG_${name}.log 

"""



}


process Clone_AIRRseq_DefineClones {

input:
 set val(name),file(airrFile) from g57_0_outputFileTSV0_g57_2

output:
 set val(name),file("*_clone-pass.tsv")  into g57_2_outputFileTSV0_g57_1

script:
failed = params.Clone_AIRRseq_DefineClones.failed
format = params.Clone_AIRRseq_DefineClones.format
seq_field = params.Clone_AIRRseq_DefineClones.seq_field
v_field = params.Clone_AIRRseq_DefineClones.v_field
d_field = params.Clone_AIRRseq_DefineClones.d_field
j_field = params.Clone_AIRRseq_DefineClones.j_field
group_fields = params.Clone_AIRRseq_DefineClones.group_fields

mode = params.Clone_AIRRseq_DefineClones.mode
dist = params.Clone_AIRRseq_DefineClones.dist
norm = params.Clone_AIRRseq_DefineClones.norm
act = params.Clone_AIRRseq_DefineClones.act
model = params.Clone_AIRRseq_DefineClones.model
sym = params.Clone_AIRRseq_DefineClones.sym
link = params.Clone_AIRRseq_DefineClones.link
maxmiss = params.Clone_AIRRseq_DefineClones.maxmiss

failed = (failed=="true") ? "--failed" : ""
format = (format=="airr") ? "--format airr": "--format changeo"
group_fields = (group_fields=="") ? "" : "--gf ${group_fields}"
v_field = (v_field=="") ? "" : "--vf ${v_field}"
d_field = (d_field=="") ? "" : "--df ${d_field}"
j_field = (j_field=="") ? "" : "--jf ${j_field}"
seq_field = (seq_field=="") ? "" : "--sf ${seq_field}"


mode = (mode=="gene") ? "" : "--mode ${mode}"
norm = (norm=="len") ? "" : "--norn ${norm}"
act = (act=="set") ? "" : "--act ${act}"
model = (model=="ham") ? "" : "--model ${model}"
sym = (sym=="avg") ? "" : "--sym ${sym}"
link = (link=="single") ? "" : " --link ${link}"
    
	
"""
DefineClones.py -d ${airrFile} \
	${failed} \
	${format} \
	${v_field} \
	${d_field} \
	${j_field} \
	${seq_field} \
	${group_fields} \
	${mode} \
	${act} \
	${model} \
	--dist ${dist} \
	${norm} \
	${sym} \
	${link} \
	--maxmiss ${maxmiss} \
	--log DF_.log  
"""



}

g_67_germlineFastaFile_g57_1= g_67_germlineFastaFile_g57_1.ifEmpty([""]) 
g_65_germlineFastaFile_g57_1= g_65_germlineFastaFile_g57_1.ifEmpty([""]) 
g_66_germlineFastaFile_g57_1= g_66_germlineFastaFile_g57_1.ifEmpty([""]) 


process Clone_AIRRseq_Second_CreateGermlines {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_germ-pass.tsv$/) "clones/$filename"}
input:
 set val(name),file(airrFile) from g57_2_outputFileTSV0_g57_1
 set val(name1), file(v_germline_file) from g_67_germlineFastaFile_g57_1
 set val(name2), file(d_germline_file) from g_65_germlineFastaFile_g57_1
 set val(name3), file(j_germline_file) from g_66_germlineFastaFile_g57_1

output:
 set val(name),file("*_germ-pass.tsv")  into g57_1_outputFileTSV0_g57_9

script:
failed = params.Clone_AIRRseq_Second_CreateGermlines.failed
format = params.Clone_AIRRseq_Second_CreateGermlines.format
g = params.Clone_AIRRseq_Second_CreateGermlines.g
cloned = params.Clone_AIRRseq_Second_CreateGermlines.cloned
seq_field = params.Clone_AIRRseq_Second_CreateGermlines.seq_field
v_field = params.Clone_AIRRseq_Second_CreateGermlines.v_field
d_field = params.Clone_AIRRseq_Second_CreateGermlines.d_field
j_field = params.Clone_AIRRseq_Second_CreateGermlines.j_field
clone_field = params.Clone_AIRRseq_Second_CreateGermlines.clone_field


failed = (failed=="true") ? "--failed" : ""
format = (format=="airr") ? "": "--format changeo"
g = "-g ${g}"
cloned = (cloned=="false") ? "" : "--cloned"

v_field = (v_field=="") ? "" : "--vf ${v_field}"
d_field = (d_field=="") ? "" : "--df ${d_field}"
j_field = (j_field=="") ? "" : "--jf ${j_field}"
seq_field = (seq_field=="") ? "" : "--sf ${seq_field}"

"""
CreateGermlines.py \
	-d ${airrFile} \
	-r ${v_germline_file} ${d_germline_file} ${j_germline_file} \
	${failed} \
	${format} \
	${g} \
	${cloned} \
	${v_field} \
	${d_field} \
	${j_field} \
	${seq_field} \
	${clone_field} \
	--log CG_${name}.log 

"""



}


process Clone_AIRRseq_single_clone_representative {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_clone_rep-passed.tsv.*$/) "clones/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*txt$/) "logs/$filename"}
input:
 set val(name),file(airrFile) from g57_1_outputFileTSV0_g57_9

output:
 set val(outname),file("*_clone_rep-passed.tsv*")  into g57_9_outputFileTSV00
 set val(name), file("*txt")  into g57_9_logFile11

script:
outname = airrFile.toString() - '.tsv' +"_clone_rep-passed"
outfile = outname + ".tsv"

"""
#!/usr/bin/env Rscript

## functions
# find the different position between sequences

src <- 
"#include <Rcpp.h>
using namespace Rcpp;
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <unordered_set>

// [[Rcpp::export]]

int allele_diff(std::vector<std::string> germs) {
    std::vector<std::vector<char>> germs_m;
    for (const std::string& germ : germs) {
        germs_m.push_back(std::vector<char>(germ.begin(), germ.end()));
    }

    int max_length = 0;
    for (const auto& germ : germs_m) {
        max_length = std::max(max_length, static_cast<int>(germ.size()));
    }

    for (auto& germ : germs_m) {
        germ.resize(max_length, '.'); // Pad with '.' to make all germs equal length
    }

    auto setdiff_mat = [](const std::vector<char>& x) -> int {
        std::unordered_set<char> unique_chars(x.begin(), x.end());
        std::unordered_set<char> filter_chars = { '.', 'N', '-' };
        int diff_count = 0;
        for (const char& c : unique_chars) {
            if (filter_chars.find(c) == filter_chars.end()) {
                diff_count++;
            }
        }
        return diff_count;
    };

    std::vector<int> idx;
    for (int i = 0; i < max_length; i++) {
        std::vector<char> column_chars;
        for (const auto& germ : germs_m) {
            column_chars.push_back(germ[i]);
        }
        int diff_count = setdiff_mat(column_chars);
        if (diff_count > 1) {
            idx.push_back(i);
        }
    }

    return idx.size();
}"

## libraries
require(dplyr)
library(Rcpp)
library(ggplot2)
sourceCpp(code = src)

data <- readr::read_tsv("${airrFile}")

nreads <- nrow(data)

# calculating mutation between IMGT sequence and the germline sequence, selecting a single sequence to each clone with the fewest mutations
data[["mut"]] <- sapply(1:nrow(data),function(j){
	x <- c(data[['sequence_alignment']][j], data[['germline_alignment_d_mask']][j])
	allele_diff(x)
})
# filter to the fewest mutations
data <- data %>% dplyr::group_by(clone_id) %>% 
			dplyr::mutate(clone_size = n())

data <- data %>% dplyr::group_by(clone_id) %>% dplyr::slice(which.min(mut))
cat(paste0('Note: ', nrow(data),' sequences after selecting a single representative'))
readr::write_tsv(data, file = "${outfile}")

lines <- c(
    paste("START>", "Selecting clonal representative"),
    paste("PASS>", nrow(data)),
    paste("FAIL>", nreads-nrow(data)),
    paste("END>", "Selecting clonal representative"),
    "",
    ""
)

file_path <- paste("${outname}","output.txt", sep="-")

cat(lines, sep = "\n", file = file_path, append = TRUE)

cat(lines, sep = "\n")
"""
}


process parse_headers_table {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*${out}$/) "reports/$filename"}
input:
 set val(name), file(reads) from g_39_fastaFile0_g_40

output:
 set val(name),file("*${out}")  into g_40_fastaFile00

script:
method = params.parse_headers_table.method
act = params.parse_headers_table.act
args = params.parse_headers_table.args


if(method=="collapse" || method=="copy" || method=="rename" || method=="merge"){
	out="_reheader.fastq"
	act = (act=="none") ? "" : "--act ${act}"
	"""
	ParseHeaders.py  ${method} -s ${reads} ${args} ${act}
	"""
}else{
	if(method=="table"){
			out=".tab"
			"""
			ParseHeaders.py ${method} -s ${reads} ${args}
			"""	
	}else{
		out="_reheader.fastq"
		"""
		ParseHeaders.py ${method} -s ${reads} ${args}
		"""		
	}
}


}


workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}
