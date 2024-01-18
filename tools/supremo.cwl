cwlVersion: v1.2
class: CommandLineTool
id: supremo
doc: |
  Pipeline for generating mutated sequences for input into predictive models
  and for scoring variants for disruption to genome folding.
requirements:
  - class: InlineJavascriptRequirement
  - class: ShellCommandRequirement
  - class: ResourceRequirement
    ramMin: $(inputs.ram*1000)
    coresMin: $(inputs.cpu)
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/supremo:1.0.0'
baseCommand: [">&2", "python", "/opt/SuPreMo/scripts/SuPreMo.py"]
arguments: []
inputs:
  input_variant_file: { type: 'File', secondaryFiles: [{pattern: '.tbi', required: false}], inputBinding: {position: 1}, doc: "Input file (in VCF, TSV, BED, or TXT format) with variants. Can be gzipped. Coordinates should be 1-based and left open (,], except for SNPs (follow vcf 4.1/4.2 specifications)." }
  sequences: { type: 'File?', inputBinding: {position: 2, prefix: "--sequences"}, doc: "Input fasta file with sequences. This file is one outputted by SuPreMo." }
  fasta: { type: 'File?', inputBinding: {position: 2, prefix: "--fa"}, doc: "Optional path to reference genome fasta file. If not provided and not existing in data/, it will be downloaded." }
  genome:
    type:
      - 'null'
      - type: enum
        name: genome
        symbols: ["hg19","hg38"]
    inputBinding:
      position: 2
      prefix: "--genome"
    doc: |
      Genome to be used: hg19 or hg38.
  scores:
    type:
      - 'null'
      - type: array
        items:
          type: enum
          name: scores
          symbols: ["mse","corr","ssi","scc","ins","di","dec","tri","pca"]
    inputBinding:
      position: 2
      prefix: "--scores"
    doc: |
      Method(s) used to calculate disruption scores. Use abbreviations as follows:
        mse: Mean squared error
        corr: Spearman correlation
        ssi: Structural similarity index measure
        scc: Stratum adjusted correlation coefficient
        ins: Insulation
        di: Directionality index
        dec: Contact decay
        tri: Triangle method
        pca: Principal component method.
  shift_by: { type: 'int[]?', inputBinding: {position: 2, prefix: "--shift_by"}, doc: "Values for shifting prediciton windows inputted as space-separated integers (e.g. -1 0 1). Values outside of range -450000 ≤ x ≤ 450000 will be ignored. Prediction windows at the edge of chromosome arms will only be shifted in the direction that is possible (ex. for window at chrom start, a -1 shift will be treated as a 1 shift since it is not possible to shift left)." }
  output_basename: { type: 'string?', inputBinding: {position: 2, prefix: "--file"}, doc: "Prefix for output file names." }
  output_dirname: { type: 'string?', default: "SuPreMo_output", inputBinding: {position: 2, prefix: "--dir"}, doc: "Output directory name." }
  limit: { type: 'int?', inputBinding: {position: 2, prefix: "--limit"}, doc: "Maximum length of variants to be scored. Filtering out variants that are too big can save time and memory. If not specified, will be set to 2/3 of seq_len." }
  seq_len: { type: 'int?', inputBinding: {position: 2, prefix: "--seq_len"}, doc: "Length for sequences to generate. Default value is based on Akita requirement. If non-default value is set, get_Akita_scores must be false." }
  revcomp:
    type:
      - 'null'
      - type: enum
        name: revcomp
        symbols: ["no_revcomp","add_revcomp","only_revcomp"]
    inputBinding:
      position: 2
      prefix: "--revcomp"
    doc: |
      Option to use the reverse complement of the sequence:
        no_revcomp: no, only use the standard sequence;
        add_revcomp: yes, use both the standard sequence and its reverse complement;
        only_revcomp: yes, only use the reverse complement of the sequence.
      The reverse complement of the sequence is only taken with 0 shift.
  augment:
    type: 'boolean?'
    inputBinding:
      position: 2
      prefix: "--augment"
    doc: |
      Only applicable if --get_Akita_scores is specified. Get the mean and
      median scores from sequences with specified shifts and reverse complement. If
      augment is used but shift and revcomp are not specified, the following four
      sequences will be used:
        1) no augmentation: 0 shift and no reverse complement,
        2) +1bp shift and no reverse complement,
        3) -1bp shift and no reverse complement,
        4) 0 shift and take reverse complement.
  get_seq:
    type: 'boolean?'
    inputBinding:
      position: 2
      prefix: "--get_seq"
    doc: |
      Save sequences for the reference and alternate alleles in fa file format. If --get_seq is not specified, must specify --get_Akita_scores.

      Sequence name format: {var_index}_{shift}_{revcomp_annot}_{seq_index}_{var_rel_pos}
        var_index: input row number, followed by _0, _1, etc for each allele of variants with multiple alternate alleles;
        shift: integer that window is shifted by;
        revcomp_annot: present only if reverse complement of sequence was taken;
        seq_index: index for sequences generated for that variant: 0-1 for non-BND reference and alternate sequences and 0-2 for BND left and right reference sequence and alternate sequence;
        var_rel_pos: relative position of variant in sequence: list of two for non-BND variant positions in reference and alternate sequence and an integer for BND breakend position in reference and alternate sequences.

      There are 2-3 entries per prediction (2 for non-BND variants and 3 for BND variants).

      To read fasta file: pysam.Fastafile(filename).fetch(seqname, start, end).upper().
      To get sequence names in fasta file: pysam.Fastafile(filename).references.
  get_tracks:
    type: 'boolean?'
    inputBinding:
      position: 2
      prefix: "--get_tracks"
    doc: |
      Save disruption score tracks (448 bins) in npy file format. Only possible for mse and corr scores.

      Dictionary item name format: {var_index}_{track}_{shift}_{revcomp_annot}
        var_index: input row number, followed by _0, _1, etc for each allele of variants with multiple alternate alleles;
        track: disruption score track specified;
        shift: integer that window is shifted by;
        revcomp_annot: present only if reverse complement of sequence was taken.

      There is 1 entry per prediction: a 448x1 array.

      To read into a dictionary in python: np.load(filename, allow_pickle="TRUE").item()
  get_maps:
    type: 'boolean?'
    inputBinding:
      position: 2
      prefix: "--get_maps"
    doc: |
      Save predicted contact frequency maps in npy file format.

      Dictionary item name format: {var_index}_{shift}_{revcomp_annot}
        var_index: input row number, followed by _0, _1, etc for each allele of variants with multiple alternate alleles;
        shift: integer that window is shifted by;
        revcomp_annot: present only if reverse complement of sequence was taken.

      There is 1 entry per prediction. Each entry contains the following: 2 (3
      for chromosomal rearrangements) arrays that correspond to the upper right
      triangle of the predicted contact frequency maps, the relative variant position
      in the map, and the first coordinate of the sequence that the map corresponds
      to.

      To read into a dictionary in python: np.load(filename, allow_pickle="TRUE").item()
  get_akita_scores:
    type: 'boolean?'
    inputBinding:
      position: 2
      prefix: "--get_Akita_scores"
    doc: |
      Get disruption scores. If --get_Akita_scores is not specified, must specify
      --get_seq. Scores saved in a dataframe with the same number of rows as the
      input. For multiple alternate alleles, the scores are separated by a comma. To
      convert the scores from strings to integers, use float(x), after separating
      rows with multiple alternate alleles. Scores go up to 20 decimal points.
  nrows:
    type: 'int?'
    inputBinding:
      position: 2
      prefix: "--nrows"
    doc: |
      Number of rows (perturbations) to read at a time from input. When dealing with
      large inputs, selecting a subset of rows to read at a time allows scores to be
      saved in increments and uses less memory. Files with scores and filtered out
      variants will be temporarily saved in output direcotry. The file names will
      have a suffix corresponding to the set of nrows (0-based), for example for an
      input with 2700 rows and with nrows = 1000, there will be 3 sets. At the end of
      the run, these files will be concatenated into a comprehensive file and the
      temporary files will be removed.
  cpu: { type: 'int?', default: 16, doc: "CPUs to allocate to this task." }
  ram: { type: 'int?', default: 32, doc: "GB of RAM to allocate to this task." }

outputs:
  supremo_dir:
    type: Directory
    outputBinding:
      loadListing: deep_listing
      glob: $(inputs.output_dirname)
