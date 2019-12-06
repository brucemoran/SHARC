process minimap_ref {

  input:
  file(fasta) from Channel.value(params.fasta)

  output:
  file('*.mmi') into map_ref

  script:
  """
  minimap2 \
    -d $fastq".mmi" \
    $fasta
  """
}
process mapping {

  input:
  set val(sampleID), file(fastq) from fastq_input
  file(index) from Channel.value(map_ref)

  output:
  set val(sampleID), file('*.bam'), file('*.bai') into depth_input

  script:
  """
  minimap2 \
    -a $index \
    -x map-ont \
  	-f $fastq \
  	-t ${task.cpu} \
  | samtools view -hb > $sampleID".sorted.bam"
  samtools index $sampleID".sorted.bam"
  """
}

process depth {

  input:
  set val(sampleID), file(bam), file(bai) from depth_input
  file(index) from Channel.value(map_ref)
  file(bed) from Channel.value(params.bed)

  output:
  set val(sampleID), file(bam), file(bai), file('*.base.depth.txt') into sv_input

  script:
  """
  sambamba-depth base \
    --min-coverage=0
    -t ${task.cpu} \
    -L $bed \
    $bam > $sampleID".base.depth.txt"
  """
}


process sv {

  input:
  set val(sampleID), file(bam), file(bai), file(depth) from sv_input
  file(index) from Channel.value(map_ref)
  file(bed) from Channel.value(params.bed)

  output:
  set val(sampleID), file(depth), file('*.nanosv.vcf') into filter_input

  script:
  """
  ##config?
  NanoSV  \
    -t ${task.cpu} \
    -o $sampleID".nanosv.vcf" \
    $bam
  """
}

process vcf_filter {

  input:
  set val(sampleID), file(depth), file(vcf) from filter_input

  output:
  set val(sampleID), file(depth), file('*.nanosv.filter.vcf') into split_input

  script:
  """
  FILTER='\$7 == \"PASS\" && \$1 !~ /(Y|MT)/ && \$5 !~ /(Y|MT):/ && \$5 != \"<INS>\"'
  grep "^#" $vcf > \$OUTPUT
  awk \$FILTER $vcf >> $sampleID".nanosv.filter.vcf"
  """
}

process vcf_split {

  input:
  set val(sampleID), file(depth), file(filter_vcf) from split_input
  file(bed) from Channel.value(params.bed)

  output:
  val(sampleID) into rf_input_sampleID
  file('*.split.vcf') into rf_input_vcfs
  file('*.meancov.txt') into rf_input_meancovs

  script:
  """
  ##annotate
  grep '^#' $filter_vcf > $sampleID".nanosv.filter.anno.vcf"
  python get_closest_feature.py $bed $filter_vcf >> $sampleID".nanosv.filter.anno.vcf"

  bgzip $sampleID".nanosv.filter.anno.vcf"

  ##split into chrs based on depth file
  grep -v "#" $depth | while read LINE; do
    CHR=\$(echo \$LINE | cut -f1)
    echo \$LINE | cut -f 5 > $sampleID"."\$CHR".meancov.txt"
    tabix -p $sampleID".nanosv.filter.anno.vcf.gz"
    tabix $sampleID".nanosv.filter.anno.vcf.gz" \$CHR > $sampleID"."\$CHR".split.vcf"
  done
  """
}

rf_input_vcfs.mix(rf_input_meancovs).into { rf_input }

process random_forest {

  input:
  val(sampleID) from rf_input_sampleID
  file(depth) from rf_input
  file(bed) from Channel.value(params.bed)

  output:
  set val(sampleID), file('*.nanosv.filter.anno.vcf') into

  script:
  """
  ##annotate features
  python create_features_table.py $anno_vcf > $sampleID".features_table.txt"

  Rscript --vanilla run_randomForest.R \
    randomForest.RData \
    $sampleID".features_table.txt"
    \$MEANCOV

  python $ADD_PREDICT_SCRIPT $VCF $OUTDIR/predict_labels.txt > $OUTPUT

  """
}
