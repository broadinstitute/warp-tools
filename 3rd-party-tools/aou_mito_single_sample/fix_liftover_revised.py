import hail as hl
import pandas as pd
import re
import argparse
import subprocess
from copy import deepcopy
from datetime import datetime
hl._set_flags(no_whole_stage_codegen='1')


NONREF = '<NON_REF>'
NONINDEL = ['N','A','T','G','C', NONREF]
MTREF = 'mt_GRCh38'
META_KEYS = ['Description', 'Number', 'Type']


def global_modify_meta(meta):
    """
    A global bookkeeping of how the metadata has been updated. All functions that modify the VCF make these changes.
    Provide the original metadata to these functions and run this at the end, AFTER the metadata has
    been trimmed for missing fields.
    """
    _, _, _, ROW_drop_fields, _, _, ENTRY_drop_fields = get_fields_for_flipping(meta)
    meta = deepcopy(meta)
    meta['format'].update({'OriginalSelfRefAlleles': {'Description': 'For allele individual pairs, these are the original alleles (called relative to self reference) which needed swapping for mapping to GRCh38.',
                                                      'Number': 'R', 'Type': 'String'},
                           'SwappedFieldIDs': {'Description': 'These fields were reversed in order (or flipped) when mapping from self reference back to GRCh38.',
                                               'Number': '1', 'Type': 'String'}})
    
    info_dict = meta['info']
    for x in ROW_drop_fields:
        _ = info_dict.pop(x, None)
    format_dict = meta['format']
    for x in ENTRY_drop_fields:
        _ = format_dict.pop(x, None)
    meta.update({'info': info_dict, 'format': format_dict})

    TLOD_entry = meta['info']['TLOD']
    meta['info'].update({'TLOD': {'Description': TLOD_entry['Description'] + ' Set to missing for reversed alleles.',
                                  'Number': TLOD_entry['Number'], 'Type': TLOD_entry['Type']}})

    meta['filter'].update({'FailedPicardLiftoverVcf' : {'Description': 'If the variant failed LiftoverVCF.', 'Number': '1', 'Type': 'String'},
                           'InsertionSharesForceCalledDeletion': {'Description': 'An insertion that ordinarily passes Liftover but shares its first base with a force called deletion (an insertion in GRCh38).',
                                                'Number': '1', 'Type': 'String'},
                           'InsertionSharesForceCalledInsertion': {'Description': 'An insertion that ordinarily passes Liftover but shares its first base with a force called insertion (a deletion in GRCh38).',
                                                'Number': '1', 'Type': 'String'},
                           'AddGRCh38RefDeleToRefSiteIns': {'Description': 'An insertion sharing a first base with GRCh38 deletion (force called insertion) resolved by replacing REF with the GRCh38 deletion REF. The self ref ALT must start with the same sequence as the GRCh38 REF. Should be fixed by left shifting.',
                                                            'Number': '1', 'Type': 'String'},
                           'ComplexSwapField': {'Description': 'In this instance, the self reference ALT allele was not the GRCh38 reference, so a more complicated approach was taken to reproduce the fields in SwappedFields. For insertions, this means constructing new ALT alleles.',
                                                'Number': '1', 'Type': 'String'},
                           'NewInsertionHaplotype': {'Description': 'Here, variants modified an insertion relative to GRCh38 (called reference in the self reference). To return to self reference, a custom haplotype was created here. ' +\
                                                                    'For example, if the homoplasmic insertion on GRCh38 is A>ATG, and in the self reference calls we see a T>TA, we need to combine these to get a A>ATAG variant.',
                                                     'Number': '1', 'Type':'String'},
                           'SwapFirstAlleleIndel': {'Description': 'This variant is an indel in which the first site was reference in self reference but not in GRCh38. The first site was swapped to GRCh38.',
                                                    'Number': '1', 'Type':'String'},
                           'ReplaceInternalBaseDeletion': {'Description': 'This variant is a deletion in which some internal bases (e.g., not the first base) were replaced to match GRCh38.',
                                                            'Number': '1', 'Type':'String'},
                           'FancyFieldInversion': {'Description': 'Added to records which were reversed and for which AF is not just 1 minus original, but rather total minus original where total is 1 minus sum(AF) of all other variants at that site.',
                                                   'Number': '1', 'Type':'String'},
                           'DeletionSpannedHomoplasmicInsertion': {'Description': 'Added to records which are heteroplasmic deletions that spanned homoplasmic insertions on only one side.',
                                                                   'Number': '1', 'Type':'String'},
                           'LiftoverSuccessEntrySwap': {'Description': 'These are records which passed Liftover but shared a position with a homoplasmic indel and thus needed entry field swapping.',
                                                        'Number': '1', 'Type':'String'},
                           'ForceCalledHomoplasmy': {'Description': 'A homoplasmy in the self minus reference that was reversed and force minus called.',
                                                     'Number': '1', 'Type':'String'},
                           'LeftShiftedIndel': {'Description': 'An indel that was modified by bcftools norm.',
                                                'Number': '1', 'Type': 'String'},
                           'FailedDuplicateVariant' : {'Description': 'Variants that were lifted using custom pipeline but became duplicate with variants ' + \
                                                                      'lifted successfully via Picard during custom LiftOver approach. ' + \
                                                                      'We keep the Picard LiftOver variant and reject the one created via the custom pipeline. ' + \
                                                                      'Currently only added to the rejected VCF (without error) if duplication ' + \
                                                                      'occurs with a non homoplasmy that was left shifted. Locus is in self ref coordinates, alleles are ' + \
                                                                      'the target alleles for GRCh38 (so may not be accurate for self reference!!!).', 
                                                      'Number': '1', 'Type': 'String'}})

    disallowed_chars = ['"', '-']
    for k, v in meta.items():
        for k2, v2 in v.items():
            descr = v2['Description']
            for char in disallowed_chars:
                descr = descr.replace(char, ' ')
            v2.update({'Description': descr})
            v.update({k2: v2})
        meta.update({k: v})
    
    return meta


def global_consistancy_checks(mt, allow_NONREF, genome, debug):
    # apply checks to info fields
    if not debug and mt.aggregate_rows(hl.agg.any(mt.info.SwappedAlleles)):
        raise ValueError('ERROR: at least one rejected row has SwappedAlleles set to True. This is abnormal and unsupported.')
    
    # AF field should never be NA
    if not debug and mt.aggregate_entries(~hl.agg.all(hl.if_else(hl.all(mt.AF.map(hl.is_defined)),True,False,missing_false=True))):
        raise ValueError('ERROR: AF should not contain any missing values.')

    # apply checks to alleles
    # alleles must all be length 2
    if mt.aggregate_rows(~hl.agg.all(hl.len(mt.alleles) == 2)):
        raise ValueError('ERROR: all alleles for swapping must have exactly 2 entries.')

    # at least one allele should be non-INDEL
    if mt.aggregate_rows(~hl.agg.all(hl.any(hl.map(lambda x: hl.literal(NONINDEL).contains(x), mt.alleles)))):
        raise NotImplementedError('ERROR: only sites that have at least one non-INDEL allele are supported.')

    # there should not be any duplicates
    if mt.count_rows() != mt.rows().distinct().count():
        raise ValueError('ERROR: input VCF should have distinct row keys. Duplicates found.')

    # aside from NONREF, the first character of both alleles should be the same
    has_nonref = hl.any(hl.map(lambda x: x == NONREF, mt.alleles))
    same_first_character = mt.alleles[0][0] == mt.alleles[1][0]
    is_indel = hl.any(hl.map(lambda x: ~hl.literal(NONINDEL).contains(x), mt.alleles))
    if mt.aggregate_rows(~hl.agg.all(~is_indel | (has_nonref | same_first_character))):
        raise ValueError('ERROR: All alleles, except NONREF, should share the same first character.')

    # if NONREF is not allowed, throw an error if any NONREF alleles are found
    if not allow_NONREF:
        if mt.aggregate_rows(hl.agg.any(has_nonref)):
            raise ValueError('ERROR: NON_REF calls are not allowed. Either re-check calling pipeline or set --allow-nonref.')
    
    # all reference alleles must match expectation
    confirm_ref(mt, genome)


def compatiblify_sample_name(s):
    s = re.sub('[^A-Za-z0-9]{1,}', '', s)
    if len(s) == 0:
        return 'sref'
    elif re.search('^[0-9]', s):
        return 'X' + s
    else:
        return s


def confirm_reference_indels(mt, indels):
    indels = indels.annotate(self_ref_alleles = [indels.alt_allele, indels.ref_allele]).key_by('coord_s', 'self_ref_alleles')
    indels_not_found = indels.anti_join(mt.rows())
    if indels_not_found.count() > 0:
        str_for_error = 'ERROR: The following self-reference based indels were expected to be found in the rejected Liftover variants set, but were not found:' + \
            ', '.join((hl.str(indels_not_found.coord_s.position) + indels_not_found.alt_allele + indels_not_found.ref_allele).collect())
        raise ValueError(str_for_error)


def fai_to_len(fai):
    with open(fai) as f:
        line = f.readline()
    return int(line.split('\t')[1])


def get_nonkey_fields(ht):
    return [x for x in ht.row if x not in ht.key]


def transpose_locus(locus, by, ref):
    """ Shifts a locus by "by" bases. by can be negative to shift backwards.
    """
    return hl.locus(locus.contig, locus.position+by, ref)


def add_filter(mt, filter_name):
    """
    Adds a filter to all rows of an mt, removing PASS when present.
    Throws an error if the filter was already present.
    """
    if mt.aggregate_rows(hl.agg.any(mt.filters.contains(filter_name))):
        raise ValueError(f"ERROR: tried to add filter {filter_name} which was already present in a record.")
    mt2 = mt.annotate_rows(filters = mt.filters.add(filter_name))
    mt2 = mt2.annotate_rows(filters = mt2.filters.remove('PASS'))
    return mt2


def drop_info(mt, to_drop):
    return mt.annotate_rows(info = mt.info.drop(*to_drop))


def is_snv(mt):
    """
    Returns a BooleanExpression of rows that are SNVs.
    """
    return hl.all(mt.alleles.map(lambda x: hl.literal(NONINDEL).contains(x)))


def is_insertion(table):
    return hl.literal(NONINDEL).contains(table.alleles[0]) & ~hl.literal(NONINDEL).contains(table.alleles[1])


def is_deletion(table):
    return hl.literal(NONINDEL).contains(table.alleles[1]) & ~hl.literal(NONINDEL).contains(table.alleles[0])


def get_fields_for_flipping(meta):
    ROW_R_fields = [k for k,v in meta['info'].items() if v['Number'] == 'R']
    ROW_BOOL_fields = ['SwappedAlleles']
    ROW_custom_fields = {'AS_SB_TABLE':'\\|'} # in field : delimiter format
    ROW_drop_fields = []
    ENTRY_R_fields = [k for k,v in meta['format'].items() if v['Number'] == 'R']
    ENTRY_1minus_fields = ['AF'] # assumes that these fields are entry Arrays
    ENTRY_drop_fields = []

    return ROW_R_fields, ROW_BOOL_fields, ROW_custom_fields, ROW_drop_fields, ENTRY_R_fields, ENTRY_1minus_fields, ENTRY_drop_fields


def add_missing_new_entry_fields(mt):
    to_add = ['OriginalSelfRefAlleles', 'SwappedFieldIDs']
    if all([x not in mt.entry for x in to_add]):
        mt2 = mt.annotate_entries(OriginalSelfRefAlleles = hl.missing('array<str>'),
                                  SwappedFieldIDs = hl.missing('tstr'))
    else:
        raise NotImplementedError('add_missing_new_entry_fields() is only intended for adding fields to an mt without them.')
    return mt2


def enforce_corrected_spec_in_uncorrected_mt(mt, skip_entry_add=False):
    if mt.count_rows() > 0:
        raise NotImplementedError('enforce_corrected_spec_in_uncorrected_mt() is only intended for use on empty MatrixTables.')
    mt2 = mt.key_rows_by().drop('locus','alleles','grch38_seq').rename({'locus_grch38':'locus'})
    mt2 = mt2.annotate_rows(alleles = hl.missing('array<str>')).key_rows_by('locus','alleles')
    if not skip_entry_add:
        mt2 = add_missing_new_entry_fields(mt2)

    return mt2


def unify_dicts(dict1, dict2, to_keep):
    """ Does NOT resolve if the values of each key in to_keep are the same. Throws an error.
    """
    dict1_copy = deepcopy(dict1)
    dict2_copy = deepcopy(dict2)

    to_keep_in_1 = [x for x in to_keep if x in dict1_copy.keys()]
    to_keep_in_2 = [x for x in to_keep if x in dict2_copy.keys()]
    to_keep_in_both = [x for x in to_keep_in_1 if x in to_keep_in_2]

    # make sure overlaps are the same
    if not all([dict1_copy[x][y] == dict2_copy[x][y] for x in to_keep_in_both for y in META_KEYS]):
        raise ValueError('ERROR: overlapping values in metadata for success and failure VCF metadata are not the same.')
    
    dict_out = deepcopy(dict1_copy)
    dict_out = {key: dict_out[key] for key in to_keep_in_1}
    dict_out.update({key: dict_out[key] for key in to_keep_in_2})
    return dict_out


def unify_info(mt1, mt2, cols_to_keep):
    # subset info to only the subset of cols_to_keep present in each
    to_keep_in_1 = [x for x in cols_to_keep if x in mt1.info.keys()]
    info1 = mt1.info.select(*to_keep_in_1)
    to_keep_in_2 = [x for x in cols_to_keep if x in mt2.info.keys()]
    info2 = mt2.info.select(*to_keep_in_2)

    #force info1 to have all fields of interest
    to_add_to_1 = [x for x in to_keep_in_2 if x not in to_keep_in_1]
    info1_f = info1.annotate(**{x: hl.missing(info2[x].dtype) for x in to_add_to_1})
    # force info2 to have the same spec as info1_f
    info2_f = force_spec(info1_f, info2)
    
    mt1 = mt1.annotate_rows(info = info1_f)
    mt2 = mt2.annotate_rows(info = info2_f)
    return mt1.persist(), mt2.persist()


def force_spec(reference, to_force):
    """
    Removes or adds fields to "to_force" so that the spec matches reference.
    Assumes to_force and reference are structs.
    """
    to_add = {}
    for col in reference.keys():
        if col not in to_force.keys():
            to_add.update({col: hl.missing(reference[col].dtype)})

    to_force2 = to_force.annotate(**to_add)
    return to_force2.select(*reference.keys())


def read_mito_vcf(vcf, reference, avoid_dropping_gt=False, debug=False):
    vcf_mt = hl.import_vcf(vcf, reference_genome=reference)
    vcf_meta = hl.get_vcf_metadata(vcf)
    if (not avoid_dropping_gt) and ('GT' in vcf_mt.entry):
        vcf_mt = vcf_mt.drop('GT')
        del vcf_meta['format']['GT']
    if debug and 'AF' not in vcf_mt.entry:
        vcf_mt = vcf_mt.annotate_entries(AF = hl.missing('array<float64>'))
    return vcf_mt, vcf_meta


def read_chain_file(chain_file, reference, self_reference):
    """
    Imports a chain file while ensuring that unsupported features are not found, namely:
    - Chain files must have no breaks
    - Chain files must always start at 0 in both source and target
    - Single block chain files only
    """
    # confirm that the chain file does not have any offset starts or breaks
    with open(chain_file) as f:
        lines = f.readlines()
    lines_filt = [line for line in lines if re.search('^chain', line)]
    if len(lines_filt) != 1:
        raise ValueError('ERROR: only single block chain files are currently supported.')
    ls = [line.split() for line in lines_filt]
    tf = [(linesplit[5] == "0") & (linesplit[10] == "0") for linesplit in ls]
    if not all(tf):
        raise ValueError('ERROR: chain file must start at 0 for both source and target.')

    chain_in = pd.read_csv(chain_file, sep='\t', skiprows=1, names=['bases_in_block', 'ns', 'nt'])
    chain_in.iat[chain_in.shape[0]-1, 1] = 0
    chain_in.iat[chain_in.shape[0]-1, 2] = 0
    chain_in['coord_s'] = chain_in.bases_in_block.cumsum() + chain_in.ns.shift(1).fillna(0).cumsum()
    chain_in['coord_t'] = chain_in.bases_in_block.cumsum() + chain_in.nt.shift(1).fillna(0).cumsum()
    chain_in = chain_in.astype(int)

    if chain_in[(chain_in.ns > 0) & (chain_in.nt > 0)].shape[0] != 0:
        raise ValueError('Chain files should not have breaks -- e.g., each ungapped block should be followed by a gap in EITHER the source or target, but not both.')

    ht = hl.Table.from_pandas(chain_in)
    ht = ht.annotate(coord_s = hl.locus('chrM', hl.int32(ht.coord_s), self_reference),
                     coord_t = hl.locus('chrM', hl.int32(ht.coord_t), reference))
    return ht


def get_insertion_sites(ht, reference, self_reference):
    """ Outputs locus-allele table of locations where insertions were made relative to reference genome.
    Chain file should indicate how to convert from self_reference to reference.
    """
    ht2 = ht.filter(ht.ns > 0)
    ht2 = ht2.annotate(ref_allele = hl.get_sequence('chrM', ht2.coord_t.position, 
                                                    before=0, after=0, reference_genome=reference),
                     alt_allele = hl.get_sequence('chrM', ht2.coord_s.position, 
                                                  before=0, after=hl.int32(ht2.ns), reference_genome=self_reference))
    
    httest = ht2.annotate(tf = ht2.alt_allele.matches("^" + ht2.ref_allele))
    if httest.aggregate(hl.agg.any(~httest.tf)):
        raise ValueError('ERROR: found instance where the first allele of insertions in self-reference (relative to GRCh38) differs from the corresponding site in GRCh38. This does not make sense.')
    
    return ht2


def get_deletion_sites(ht, reference, self_reference):
    """ Outputs locus-allele table of locations where deletions were made relative to reference genome.
    Chain file should indicate how to convert from self_reference to reference.
    """
    ht2 = ht.filter(ht.nt > 0)
    ht2 = ht2.annotate(ref_allele = hl.get_sequence('chrM', ht2.coord_t.position, 
                                                    before=0, after=hl.int32(ht2.nt), reference_genome=reference),
                       alt_allele = hl.get_sequence('chrM', ht2.coord_s.position, 
                                                    before=0, after=0, reference_genome=self_reference))

    httest = ht2.annotate(tf = ht2.ref_allele.matches("^" + ht2.alt_allele))
    if httest.aggregate(hl.agg.any(~httest.tf)):
        raise ValueError('ERROR: found instance where the first allele of deletions in self-reference (relative to GRCh38) differs from the left-most nucleotide of the deleted seqeunce in GRCh38. This does not make sense.')
    
    return ht2


def explode_indel(ht, locus_expr, len_expr, reference, target_name):
    """ All provided expressions should refer to the INSERTION. Thus, if
    the self genome has an insertion relative to GRCh38, provide expressions
    corresponding to self.
    """
    ht2 = ht.annotate(tmp = hl.range(locus_expr.position, hl.int32(len_expr+locus_expr.position+1)),
                      coord_end = hl.locus('chrM', hl.int32(len_expr+locus_expr.position), reference))
    ht2 = ht2.explode('tmp')
    ht2 = ht2.annotate(tmp = hl.locus('chrM', ht2.tmp, reference))
    ht2 = ht2.annotate(tmp_2 = hl.get_sequence('chrM', ht2.tmp.position, 
                                               before=0, after=0, reference_genome=reference))
    ht2 = ht2.rename({'tmp': target_name, 'tmp_2': f'{target_name}_seq'}).key_by(target_name)

    if ht2.count() != ht2.distinct().count():
        raise ValueError('ERROR: When exploding indels, duplicate keys were created.')
    
    return ht2


def dele_spans_insertion(mt, exploded_insertions, reference_genome, ignore_identical, filter_to_same_start_span_rhs=False, filter_to_lhs_span=False):
    """ Returns a filtered (but otherwise unchanged) mt where rows span an insertion.
    This means if that row has the start/end position outside the insertion and the other inside.
    Some specific cases and expected behavior, for a reference insertion GRCh38 513 G>GCACA -> Self 513-517 GCACA
    1. Self 503:ATCC>A --- not spanning
    2. Self 509:CCCA>C --- not spanning
    3. Self 509:CCCAG>C --- not spanning**
    4. Self 511:CAGC>C --- spanning
    5. Self 511:CAGCA>C --- spanning
    6. Self 511:CAGCACA>G --- not spanning**
    7. Self 511:CAGCACAC>G --- not spanning
    8. Self 513:GCACA>G --- handled by ignore_identical
    9. Self 513:GCACAC>G --- spanning
    10. Self 513:GC>G --- not spanning
    11. Self 516:CA>C --- not spanning
    12. Self 517:ACA>A --- spanning
    13. Self 518:CACA>C --- not spanning

    We have tested this with the new filter_to_lhs_span flag. This flag returns cases where the start
    of the deletion is upstream of the insertion; it will not return cases where the deletion
    ends at the end of the insertion.

    ignore_identical: if True, will filter out alleles in mt that are the same as the reference insertion
    """
    expl_tmp = exploded_insertions.annotate(is_in_insert = True)
    dedup_expl = expl_tmp.key_by('coord_s', 'ref_allele', 'alt_allele').distinct().key_by('coord_s')
    ht = mt.rows()
    ht = ht.filter(hl.literal(NONINDEL).contains(ht.alleles[1]) & ~hl.literal(NONINDEL).contains(ht.alleles[0]))
    
    # for testing with K1481
    # df = pd.DataFrame({'pos':[503,509,509,511,511,511,511,513,513,513,516,517,518],
    #                    'ref':['ATCC','CCCA','CCCAG','CAGC','CAGCA','CAGCACA','CAGCACAC','GCACA','GCACAC','GC','CA','ACA','CACA'],
    #                    'alt':['A','C','C','C','C','G','G','G','G','G','C','A','C']})
    # ht_custom = hl.Table.from_pandas(df)
    # ht = ht_custom.annotate(locus = hl.locus('chrM',hl.int32(ht_custom.pos),reference_genome), 
    #                         alleles=[ht_custom.ref, ht_custom.alt]).key_by('locus','alleles')
    # df_dedup_expl = pd.DataFrame({'bases_in_block': [6], 'ns': [4], 'nt': [0], 
    #                               'coord_s': [hl.locus('chrM', 513, reference_genome=reference_genome)],
    #                               'coord_t': [hl.locus('chrM', 513, reference_genome='GRCh38')],
    #                               'coord_end': [hl.locus('chrM', 517, reference_genome=reference_genome)],
    #                               'ref_allele': ['a'], 'alt_allele': ['b'],
    #                               'is_in_insert': [True]})
    # dedup_expl = hl.Table.from_pandas(df_dedup_expl).key_by('coord_s')
    # expl_tmp = dedup_expl.annotate(sites = hl.map(lambda x: hl.locus('chrM', x, reference_genome=reference_genome), hl.range(513, 513+5)))
    # expl_tmp = expl_tmp.explode('sites').key_by('sites')
    
    if ignore_identical:
        dedup_for_filt = dedup_expl.annotate(alleles = [dedup_expl.alt_allele, dedup_expl.ref_allele]).key_by('coord_s','alleles')
        ht = ht.anti_join(dedup_for_filt)
    ht = ht.annotate(dele_start = ht.locus)
    ht = ht.annotate(dele_end = hl.locus('chrM', ht.dele_start.position + hl.len(ht.alleles[0]) - 1, reference_genome=reference_genome))
    if ht.aggregate(~hl.agg.all(hl.is_defined(ht.dele_start) & hl.is_defined(ht.dele_end))):
        raise ValueError('ERROR: Deletion start and end coordinates in self-reference came up as undefined.')
    ht = ht.annotate(start_in = expl_tmp[ht.dele_start].is_in_insert, end_in = expl_tmp[ht.dele_end].is_in_insert,
                     start_at = dedup_expl[ht.dele_start].is_in_insert, end_at = (dedup_expl.key_by('coord_end'))[ht.dele_end].is_in_insert,
                     end_at_start = dedup_expl[ht.dele_end].is_in_insert).persist()
    if filter_to_same_start_span_rhs:
        tf_rhs = hl.is_defined(ht.start_at) & ~hl.is_defined(ht.end_at) & ~hl.is_defined(ht.end_in)
        ht = ht.filter(tf_rhs)
    elif filter_to_lhs_span:
        tf_lhs = ~hl.is_defined(ht.start_in) & ~hl.is_defined(ht.start_at) & hl.is_defined(ht.end_in) & ~(hl.is_defined(ht.end_at) | hl.is_defined(ht.end_at_start))
        ht = ht.filter(tf_lhs)
    else:
        tf_basic = (hl.is_defined(ht.start_in) & ~hl.is_defined(ht.end_in)) | \
                    (~hl.is_defined(ht.start_in) & hl.is_defined(ht.end_in)) | \
                    (hl.is_defined(ht.start_at) & hl.is_defined(ht.end_at))
        tf_extra = ~hl.is_defined(ht.start_at) & ~hl.is_defined(ht.start_in) & \
                (hl.is_defined(ht.end_at) | hl.is_defined(ht.end_at_start))
        ht = ht.filter(tf_basic & ~tf_extra)
    return mt.semi_join_rows(ht)


def get_ref_locus_end(locus, allele, exploded_insertions, self_ref, ref):
    """ Gets the end position of a self-ref allele in ref coordinates.
    Also explicitly accounts for the case where the end of the allele is the end of
    an insertion relative to reference.
    """
    exploded_insertions_dedup = exploded_insertions.key_by('coord_end').distinct()
    self_ref_position_last_base = locus.position + hl.len(allele) - 1
    self_ref_locus_last_base = hl.locus('chrM', self_ref_position_last_base, reference_genome=self_ref)
    last_base_is_insertion_end = hl.is_defined(exploded_insertions_dedup[self_ref_locus_last_base])
    locus_last_base_ref = hl.if_else(last_base_is_insertion_end, 
                                     transpose_locus(hl.liftover(transpose_locus(self_ref_locus_last_base, 1, self_ref), ref), -1, ref),
                                     hl.liftover(self_ref_locus_last_base, ref))
    return locus_last_base_ref


def check_missing_row_field(mt, expr):
    if re.search('^array<.+>$', str(expr.dtype)):
        # this is an array
        return mt.aggregate_rows(hl.agg.any(hl.any(hl.map(lambda x: hl.is_defined(x), expr))))
    else:
        return mt.aggregate_rows(hl.agg.any(hl.is_defined(expr)))


def inject_success_variants_to_fix(failed_mt, success_mt, homoplasmies_ht, self_ref):
    """ Scans success mt for any variants that should actually fail.
    These approaches are bespoke, and will have to be added on a case-by-case basis.
    Currently supported:
    - If an insertion is present at the first base of a force-called deletion
    - If an insertion is present at the first base of a force-called insertion

    NOTE this expects that failed_mt and success_mt share the same spec
    """
    success_mt = success_mt.annotate_rows(locus_self = hl.liftover(success_mt.locus, self_ref))
    success_mt = success_mt.key_rows_by('locus_self')
    n_original_success = success_mt.count_rows()
    n_original_failed = failed_mt.count_rows()
    n_changed = 0

    if success_mt.aggregate_rows(~hl.agg.all(hl.is_defined(success_mt.locus) & hl.is_defined(success_mt.alleles) & hl.is_defined(success_mt.locus_self))):
        raise ValueError('Success MT must have defined locus, alleles, and self lifted locus.')
    
    # insertion present at the first base of a force-called deletion
    force_called_dele = homoplasmies_ht.filter(is_deletion(homoplasmies_ht))
    force_called_dele = force_called_dele.key_by()
    success_mt_f1 = success_mt.filter_rows(is_insertion(success_mt)).semi_join_rows(force_called_dele.key_by('locus'))
    success_mt_f1_merge = success_mt_f1.key_rows_by().drop('locus').rename({'locus_self':'locus'})
    success_mt_f1_merge = success_mt_f1_merge.select_rows(*[x for x in failed_mt.row]).key_rows_by('locus','alleles')
    failed_mt = failed_mt.union_rows(add_filter(success_mt_f1_merge, 'InsertionSharesForceCalledDeletion'))
    success_mt_tmp = success_mt.filter_rows(~is_insertion(success_mt))
    success_mt = success_mt_tmp.union_rows(success_mt.filter_rows(is_insertion(success_mt)).anti_join_rows(force_called_dele.key_by('locus')))
    n_changed+=success_mt_f1.count_rows()

    # insertion is at the first base of a force-called insertion (reference deletion)
    # NOTE we currently only deal with a very specific case of insertion success spike-in. It must:
    # - be mapped to a reference deletion (e.g., share a locus)
    # - share the same first base as the reference deletion
    # - if the reference deletion has a REF allele of length N, the first N bases of the ALT allele of the insertion should be identical to the reference deletion REF allele
    # For example: REF HOM GCA > G and SELF HET G > GCACA -> GCA > GCACA -> G > GCA
    # Others will be allowed in here, but will fail later.
    force_called_ins = homoplasmies_ht.filter(is_insertion(homoplasmies_ht))
    force_called_ins = force_called_ins.key_by()

    success_mt_f2 = success_mt.filter_rows(is_insertion(success_mt)).semi_join_rows(force_called_ins.key_by('locus'))
    success_mt_f2_merge = success_mt_f2.key_rows_by().drop('locus').rename({'locus_self':'locus'})
    success_mt_f2_merge = success_mt_f2_merge.select_rows(*[x for x in failed_mt.row]).key_rows_by('locus','alleles')
    failed_mt = failed_mt.union_rows(add_filter(success_mt_f2_merge, 'InsertionSharesForceCalledInsertion'))
    success_mt_tmp = success_mt.filter_rows(~is_insertion(success_mt))
    success_mt = success_mt_tmp.union_rows(success_mt.filter_rows(is_insertion(success_mt)).anti_join_rows(force_called_ins.key_by('locus')))
    n_changed+=success_mt_f2.count_rows()

    success_mt = success_mt.key_rows_by(*['locus','alleles']).drop('locus_self')
    
    # bookkeeping
    n_final_success = success_mt.count_rows() 
    n_final_failed = failed_mt.count_rows()
    if n_original_success - n_final_success != n_changed:
        raise ValueError('ERROR: Unexpected number of records lost from success MT when injecting success variants into failed MT.')
    if n_final_failed - n_original_failed != n_changed:
        raise ValueError('ERROR: Unexpected number of records gained in success MT when injecting success variants into failed MT.')
    if n_final_failed + n_final_success != n_original_failed + n_original_success:
        raise ValueError('ERROR: Records lost or gained when injecting success variants into failed MT.')
    
    return failed_mt.persist(), success_mt.persist(), n_changed


def swap_alleles(mt, meta_file, allow_NONREF, log, fail_on_complex=False, suppress_complex_warning=False, original_self_alleles=None):
    mt = mt.persist()
    if mt.count_rows() == 0:
        return enforce_corrected_spec_in_uncorrected_mt(mt), 0
    
    # identify fields with values per allele including reference (R) -- info fields to flip
    ROW_R_fields, ROW_BOOL_fields, ROW_custom_fields, _, ENTRY_R_fields, ENTRY_1minus_fields, _ = get_fields_for_flipping(meta_file)
    
    # grch38 coordinates must be defined
    if mt.aggregate_rows(~hl.agg.all(hl.is_defined(mt.locus_grch38) & hl.is_defined(mt.grch38_seq))):
        raise ValueError('ERROR: GRCh38 locus and sequence at every site for swap must be defined.')
    # the current reference must not be the grch38 reference
    if mt.aggregate_rows(~hl.agg.all(mt.alleles[0] != mt.grch38_seq)):
        raise ValueError('ERROR: GRCh38 sequence matches self-reference sequence at >1 site, which makes no sense given these are tagged as MismatchedRefAllele.')
    # the current alternate should either be <NON_REF> or the grch38 reference
    # if not, a warning will be raised as this CAN happen when a self-reference changed allele is altered
    straightforward_case = [NONREF, mt.grch38_seq] if allow_NONREF else [mt.grch38_seq]
    mt = mt.annotate_rows(to_test = straightforward_case)
    ct_fail_alt = mt.filter_rows(~mt.to_test.contains(mt.alleles[1])).count_rows()
    if ct_fail_alt > 0:
        warn_str = 'WARNING: ' + str(ct_fail_alt) + ' self-reference alternate sites did not match the GRCh38 reference.'
        if fail_on_complex:
            raise ValueError(warn_str + ' This is not supported.')
        else:
            if not suppress_complex_warning:
                print(warn_str + ' We will replace the reference with GRCh38 and replace the array fields with corresponding values from reference calls.', file=log)
        # this identifies specific locus / allele pairs with issues
        # locus/allele pairs without isssue will be sent to the straightforward pipeline
        loci_of_issue = mt.filter_rows(~mt.to_test.contains(mt.alleles[1])).rows()        
        # meanwhile, here we get any site sharing a locus with the problematic sites
        loci_of_issue_2 = loci_of_issue.annotate(tf_issue=True)
        mt_with_issue = mt.key_rows_by('locus').semi_join_rows(loci_of_issue.key_by('locus')).key_rows_by('locus','alleles')
        mt_with_issue = mt_with_issue.annotate_rows(row_with_issue = hl.is_defined(loci_of_issue_2[mt_with_issue.row_key].tf_issue))
        
        # Now deal with other fields: for each locus, does the reference get called?
        # GRCh38 reference will only be called in the ALT position (asserted above)
        # if so, pull them from there
        mt_with_issue = mt_with_issue.annotate_rows(grch38_is_called = mt_with_issue.grch38_seq == mt_with_issue.alleles[1])
        
        # if not, throw an error -- forcing calls at all reference sites should ensure this does NOT happen
        ht_ref_specified = mt_with_issue.group_rows_by(mt_with_issue.locus).aggregate(n_with_grch38 = hl.agg.sum(~mt_with_issue.grch38_is_called)).entries()
        n_ref_not_specified = ht_ref_specified.filter(ht_ref_specified.n_with_grch38 == 0).count()
        if n_ref_not_specified > 0:
            raise ValueError('ERROR: ' + str(n_ref_not_specified) + ' sites with non-GRCh38 ALT calls have no GRCh38 calls anywhere. This should never happen as long as reference alleles were force called with Mutect.')
        
        # for each locus, identify the right site from which to pull data
        loci_with_grch38_called = mt_with_issue.filter_rows(mt_with_issue.grch38_is_called).key_rows_by('locus')
        loci_with_grch38_called = loci_with_grch38_called.annotate_rows(row_found = True)

        # confirm that loci_with_grch38_called has one row per original site
        n_failed = loci_of_issue.key_by('locus').anti_join(loci_with_grch38_called.rows()).count() > 0
        if n_failed > 0:
            raise ValueError('ERROR: ' + str(n_failed) + ' sites on self-reference did not have either a ' + NONREF + ' call or an ALT matching GRCh38.')

        # update the mt
        mt_with_issue = mt_with_issue.annotate_rows(flipped_row_data = loci_with_grch38_called.rows()[mt_with_issue.locus])
        mt_with_issue = mt_with_issue.annotate_entries(flipped_entries = loci_with_grch38_called[mt_with_issue.locus, mt_with_issue.col_key])
        mt_with_issue = mt_with_issue.annotate_rows(updated_alleles = [mt_with_issue.grch38_seq, mt_with_issue.alleles[1]])
        
        # remove rows found in mt as these will be processed in the straightforward pipeline
        mt = mt.anti_join_rows(loci_of_issue)
        mt_with_issue = mt_with_issue.anti_join_rows(mt.rows())


        mt_with_issue = complex_swap_field(mt_with_issue, ROW_BOOL_fields=ROW_BOOL_fields, ROW_R_fields=ROW_R_fields, 
                                           ROW_custom_fields=ROW_custom_fields, ENTRY_R_fields=ENTRY_R_fields, 
                                           flipped_idx=1, skip_locus_change=False, original_self_alleles=original_self_alleles)
        mt_with_issue = mt_with_issue.drop('to_test', 'row_with_issue', 'grch38_is_called', 'grch38_seq')

    # for the straightforward cases:
    # remove '<NON_REF>' records if there is another option available
    htr = mt.rows()
    n_per_allele = htr.group_by(htr.locus).aggregate(count_loc = hl.agg.count())
    possible_ns = list(n_per_allele.aggregate(hl.agg.collect_as_set(n_per_allele.count_loc)))
    if any([x not in [1, 2] for x in possible_ns]):
        raise ValueError('ERROR: the straightforward case pipeline should only see at most 2 records per locus (one from ' + NONREF + 'and one from GRCh38).')
    mt = mt.annotate_rows(n_per_locus = n_per_allele[mt.locus].count_loc)
    mt = mt.filter_rows((mt.n_per_locus == 2) & (mt.alleles[1] == NONREF), keep=False)

    # change alleles
    mt = mt.annotate_rows(new_alleles = [mt.grch38_seq, mt.alleles[0]])
    mt = mt.annotate_entries(OriginalSelfRefAlleles = mt.alleles if original_self_alleles is None else mt[original_self_alleles])
    mt = mt.key_rows_by().drop(*['locus','alleles']).rename({'locus_grch38': 'locus', 'new_alleles':'alleles'}).key_rows_by('locus','alleles')

    # swap info fields
    # anything that gets reversed has a TLOD that is not longer interpretable
    info_str = mt.info.annotate(**{x: ~mt.info[x] for x in ROW_BOOL_fields})
    info_str = info_str.annotate(**{x: hl.reversed(info_str[x]) for x in ROW_R_fields})
    info_str = info_str.annotate(**{k: hl.literal(re.sub('\\\\','',v)).join(hl.reversed(info_str[k].split(v))) for k,v in ROW_custom_fields.items()})
    if 'TLOD' in info_str.keys():
        info_str = info_str.annotate(TLOD = hl.missing(info_str.TLOD.dtype))
    mt = mt.annotate_rows(info = info_str)

    # swap entry fields -- this is the only place 1-entry occurs!
    mt = mt.annotate_entries(**{x: hl.map(lambda x: 1-x, mt[x]) for x in ENTRY_1minus_fields})
    mt = mt.annotate_entries(**{x: hl.reversed(mt[x]) for x in ENTRY_R_fields})
    mt = mt.drop('grch38_seq', 'to_test', 'n_per_locus')
    mt = mt.annotate_entries(SwappedFieldIDs = ','.join([*ROW_BOOL_fields, *ROW_R_fields, *list(ROW_custom_fields.keys()), *ENTRY_R_fields, *['1-' + x for x in ENTRY_1minus_fields]]))

    if ct_fail_alt > 0:
        mt = mt.union_rows(mt_with_issue)

    return mt, ct_fail_alt


def fix_ref_insertions(mt, mt_meta, insertions_in, self_reference, allow_NONREF, log):
    mt = mt.persist()
    if mt.count_rows() == 0:
        return enforce_corrected_spec_in_uncorrected_mt(mt.drop('is_reversed_grch38_to_self')), 0
    
    mt = mt.annotate_rows(locus_grch38 = insertions_in[mt.locus].coord_t,
                          grch38_seq = insertions_in[mt.locus].ref_allele,
                          self_ref_insertion_seq = insertions_in[mt.locus].alt_allele,
                          start_coord_self_ref = insertions_in[mt.locus].coord_s,
                          end_coord_self_ref = insertions_in[mt.locus].coord_end)

    # checks
    # verifying the sites tagged as being reference
    if mt.aggregate_rows(~hl.agg.all(mt.is_reversed_grch38_to_self == ((mt.start_coord_self_ref == mt.locus) & (hl.reversed(mt.alleles) == [mt.grch38_seq, mt.self_ref_insertion_seq])))):
        raise ValueError('ERROR: the sites where expected reference is the same as the reversed self-ref site are not the same as those annotated with is_reversed_grch38_to_self.')
    # ensuring that all sites with shared locus as the reference insertion have the same first base
    mt_test = mt.filter_rows(mt.start_coord_self_ref == mt.locus)
    mt_test = mt_test.annotate_rows(per_site_unique_first_base = hl.set(mt_test.alleles.map(lambda x: x[0])))
    if mt_test.aggregate_rows(~hl.agg.all(mt_test.per_site_unique_first_base.length() == 1)):
        raise ValueError('ERROR: Each row should have alleles with the same first character.')
    mt_test = mt_test.annotate_rows(per_site_unique_first_base = hl.array(mt_test.per_site_unique_first_base)[0])
    mt_test_rows = mt_test.rows()
    ht_first_base_locus = mt_test_rows.group_by(mt_test_rows.locus
                                     ).aggregate(unique_first_base = hl.agg.collect_as_set(mt_test_rows.per_site_unique_first_base))
    if ht_first_base_locus.aggregate(~hl.agg.all(ht_first_base_locus.unique_first_base.length() == 1)):
        raise ValueError('ERROR: All sites at the same locus as the self-ref insertion relative to GRCh38 should have the same first character of both alleles.')
    mt_test = mt_test.filter_rows(~mt_test.is_reversed_grch38_to_self)
    # heteroplasmies sharing the same first site as a homoplasmic insertion (grch38) are supported only if they are indels
    if mt_test.count_rows() > 0:
        # we now allow heteroplasmic insertions that share the first site with a homoplasmic reference insertion
        #if mt_test.aggregate_rows(hl.agg.any((hl.len(mt_test.alleles[1]) > 1) | (hl.len(mt_test.alleles[0]) == 1))):
            #raise ValueError('ERROR: Only heteroplasmic deletions can be found at the first site of a left-aligned insertion (relative to GRCh38; deletion in self-reference). Investigate why this has happened here.')
        print('NOTE: There is a heteroplasmic indel at the same starting position as a homoplasmic insertion relative to GRCh38.', file=log)
    if allow_NONREF:
        raise NotImplementedError('ERROR: allow_NONREF is not implemented for fix_ref_insertions, as we currently enforce one-to-one mapping between haplotypes and variants.')
    ht_test = mt.group_rows_by(mt.locus_grch38).aggregate(n = hl.agg.count_where(mt.is_reversed_grch38_to_self)).entries()
    if ht_test.aggregate(~hl.agg.all(ht_test.n == 1)):
        raise ValueError('ERROR: Each GRCh38 locus should have exactly 1 variant that is classified as the reversed to-self-reference variant via is_reversed_grch38_to_self.')
    
    # produce haplotypes
    ht_haplos = produce_pseudo_haplotypes(mt.rows(), 'locus_grch38', MTREF, self_reference)

    # enforce a one-to-one mapping between haplotypes and self-reference variants
    ht_haplos_ct = ht_haplos.group_by(*ht_haplos.key).aggregate(n = hl.agg.count())
    if ht_haplos_ct.aggregate(~hl.agg.all(ht_haplos_ct.n == 1)):
        raise ValueError('ERROR: Haplotype table should contain one entry per self-reference locus-allele pair.')
    ht_haplos_grch38_ct = ht_haplos.group_by(ht_haplos.locus, ht_haplos.alleles).aggregate(n = hl.agg.count())
    if ht_haplos_grch38_ct.aggregate(~hl.agg.all(ht_haplos_grch38_ct.n == 1)):
        raise ValueError('ERROR: Haplotype table should contain one entry per unique haplotype.')
    # confirm that grch38 reference information all makes sense
    ht_haplos_pos0 = ht_haplos.group_by(ht_haplos.locus).aggregate(n = hl.agg.count_where(ht_haplos.position_0_source_of_alt_data))
    if ht_haplos_pos0.aggregate(~hl.agg.all(ht_haplos_pos0.n == 1)):
        raise ValueError('ERROR: position_0_source_of_alt_data should only be true for one variant per GRCh38 locus, denoting the reference allele.')
    # ensure all sites in mt got a haplotype
    if mt.anti_join_rows(ht_haplos).count_rows() > 0:
        raise ValueError('ERROR: there are sites in the rejected variant list that did not get a haplotype when dealing with insertions relative to GRCh38.')
    if ht_haplos.anti_join(mt.rows()).count() > 0:
        raise ValueError('ERROR: there are created haplotypes that do not match the original sites for fixing insertions relative to GRCh38.')
    
    # add information to mt
    mt = mt.annotate_rows(haplo_data = ht_haplos[mt.row_key])
    
    # produce warning regarding custom haplotypes and add filter
    mt_custom_haplos = mt.filter_rows(~mt.haplo_data.position_0_source_of_alt_data)
    n_new_haplotypes = mt_custom_haplos.count_rows()
    if n_new_haplotypes > 0:
        str_original = hl.str(mt_custom_haplos.locus.position) + ':' + hl.literal('>').join(mt_custom_haplos.alleles)
        str_custom_haplos = hl.str(mt_custom_haplos.haplo_data.locus.position) + ':' + hl.literal('>').join(mt_custom_haplos.haplo_data.alleles)
        str_join = ', '.join((str_original + ' -> ' + str_custom_haplos).collect())
        print('WARNING: ' + str(n_new_haplotypes) + ' new GRCh38 ALT alleles were created due to variation influencing self-reference insertions relative to GRCh38.', file=log)
        print('The new alleles are (self -> GRCh38): ' + str_join, file=log)
        mt_custom_haplos = add_filter(mt_custom_haplos, 'NewInsertionHaplotype')
        mt = mt.filter_rows(mt.haplo_data.position_0_source_of_alt_data).union_rows(mt_custom_haplos)

    # reconfigure mt for swap_alleles -- for new haplos, ensure alt allele for new haplos is the desired alt allele and shift start position back to shared start for variant
    mt = mt.annotate_rows(locus_grch38 = mt.haplo_data.locus)
    mt = mt.annotate_rows(original_self_alleles = mt.alleles)
    mt = mt.key_rows_by()
    # swap alleles assumes that alleles[0] will not be the true reference
    mt = mt.annotate_rows(fake_a0 = hl.if_else(mt.grch38_seq == 'G', 'A', 'G'))
    mt = mt.annotate_rows(alleles = hl.if_else(mt.haplo_data.position_0_source_of_alt_data, mt.alleles, 
                                               [mt.fake_a0, mt.haplo_data.alleles[1]]),
                          locus = hl.if_else(mt.haplo_data.position_0_source_of_alt_data, mt.locus, 
                                             mt.start_coord_self_ref)).key_rows_by('locus','alleles')
    mt_res, n_complex_from_swap = swap_alleles(mt, mt_meta, allow_NONREF=allow_NONREF, log=log, fail_on_complex=False, suppress_complex_warning=True, original_self_alleles='original_self_alleles')
    
    # final checks
    if n_complex_from_swap != n_new_haplotypes:
        raise ValueError('ERROR: Expected that all new haplotypes created would be reported as a complex swap.')
    if mt_res.aggregate_rows(~hl.agg.all(mt_res.haplo_data.locus == mt_res.locus)) | mt_res.aggregate_rows(~hl.agg.all(mt_res.haplo_data.alleles == mt_res.alleles)):
        raise ValueError('ERROR: swap_alleles did not produce locus / allele pairs as expected.')

    mt_res = mt_res.drop('self_ref_insertion_seq', 'start_coord_self_ref', 'end_coord_self_ref', 'is_reversed_grch38_to_self', 'haplo_data', 'original_self_alleles', 'fake_a0')
    return mt_res, n_new_haplotypes


def produce_pseudo_haplotypes(ht, incr, incr_reference, self_reference):
    """ Takes a ht with locus, allele coding and produces all possible haplotypes given variants.
    incr is used as the grouping mechanism -- rows with the same incr
    are combined into haplotypes.

    Assumes that within each incr, loci are provided sequentially.

    Also outputs alternate alleles that go in to each haplotype. Useful for
    deciding which fields to bin for recomputing INFO and FORMAT fields.

    Assumes that any variant at the first site in an insertion is the self-reference > GRCh38 conversion.

    NOTE In rare cases there may be a heteroplasmic variant located at the start site of the insertion (relative to GRCh38)
    EXAMPLE: GRCh38 500 G > GCATA (96%) -> Self 500 GCATA > G with a GCA > G at 1%
    Here we only support deletions -- insertions and SNVs (500 G > GT or 500 G > A) should be resolved with traditional Liftover
    NOTE In the above case, swapping Self 500 GCATA > G back to G > GCATA simply involves flipping alleles
    However, flipping alleles in GCA > G is inapppropriate because this G is not the same as reference G at this site; implicitly its G for GTA (the rest of the homoplasmic insertion).
    Instead, we change GCA > G* (* = self-ref G) to G > GTA; GTA gets the same information as G* while G inherits GRCh38 G information.
    """
    ht = ht.key_by()
    ht = ht.select(site = ht.locus.position,
                   site_start = ht.start_coord_self_ref.position,
                   site_end = ht.end_coord_self_ref.position,
                   incr = ht[incr].position,
                   alt_allele = ht.alleles[1],
                   alleles = ht.alleles.filter(lambda x: x != NONREF),
                   self_ref_insertion_seq = ht.self_ref_insertion_seq,
                   is_reversed_grch38_to_self = ht.is_reversed_grch38_to_self)
    df = ht.to_pandas()

    # newdf = pd.DataFrame({'site':[205, 205, 206, 207, 208, 208], 'incr':205,
    #                       'alt_allele':['T', 'T', 'C', 'ATC', 'C', 'T'], 
    #                       'alleles':[['TCAA','T'], ['TCA', 'T'], ['CAA','C'],['A','ATC'], ['A', 'C'], ['A', 'T']], 
    #                       'self_ref_insertion_seq':'TCAA',
    #                       'site_start': [205, 205, 205, 205, 205 , 205],
    #                       'site_end': [208, 208, 208, 208, 208, 208],
    #                       'is_reversed_grch38_to_self': [True, False, False, False, False, False]})
    # df = df.append(newdf)
    df = df.sort_values(by='site', axis=0)

    if not all((df.site_end - df.site_start + 1) == df.self_ref_insertion_seq.str.len()):
        raise ValueError('ERROR: site_start and site_end must point to the coordinates of the first and last elements of the provided insertion.')
    df_test = df.groupby(df['incr'], sort=False).agg(is_reversed_grch38_to_self=('is_reversed_grch38_to_self', lambda x: list(x)))
    if not all(df_test.is_reversed_grch38_to_self.map(lambda x: sum(x)) == 1):
        raise ValueError('ERROR: every GRCh38 site (incr) must have a single allele that is labeled is_reversed_grch38_to_self.')
    df = df.groupby(df['site'], sort=False).agg(incr_position = ('incr', 'unique'),
                                                alt_array=('alt_allele', lambda x: list(x)),
                                                allele_array=('alleles', lambda x: list(x)),
                                                is_reversed_grch38_to_self=('is_reversed_grch38_to_self', lambda x: list(x)),
                                                self_ref_insertion_seq=('self_ref_insertion_seq', 'unique'),
                                                site_start=('site_start','unique'),
                                                site_end=('site_end','unique'))

    if any(df.incr_position.apply(len) > 1) | any(df.self_ref_insertion_seq.apply(len) > 1) | any(df.site_start.apply(len) > 1) | any(df.site_end.apply(len) > 1):
        raise ValueError('ERROR: you should have one incr_position, self_ref_insertion_seq, site_start, site_end per self locus.')
    else:
        df.incr_position = df.incr_position.apply(lambda x: x[0])
        df.self_ref_insertion_seq = df.self_ref_insertion_seq.apply(lambda x: x[0])
        df.site_start = df.site_start.apply(lambda x: x[0])
        df.site_end = df.site_end.apply(lambda x: x[0])
        df = df.reset_index()
    
    def replace_substring(series):
        to_replace = series.self_ref_insertion_seq
        if (len(series.alt_array) != len(series.allele_array)) | (len(series.alt_array) != len(series.is_reversed_grch38_to_self)):
            raise ValueError('alt_array, allele_array, is_reversed_grch38_to_self must all have the same length.')
        grch38_ref_holder = []
        grch38_alt_holder = []
        site_of_reference_holder = []
        for alt, allele, is_rev_reference in zip(series.alt_array, series.allele_array, series.is_reversed_grch38_to_self):
            if is_rev_reference:
                if (series.site != series.site_start):
                    raise ValueError('A reference allele should start at the same site as the start of the locus.')
                grch38_ref_holder.append(alt[0])
                grch38_alt_holder.append(to_replace)
                site_of_reference_holder.append(True)
            else:
                to_replace_start = series.site - series.site_start
                if to_replace_start == 0:
                    # we now support heteroplasmic insertions
                    #if (len(allele[0]) == 1) | (len(allele[1]) > 1):
                        #raise ValueError('Only heteroplasmic deletions are supported at the first site of an insertion.')
                    if (len(allele[0]) == 0) & (len(allele[1]) == 0):
                        raise ValueError('Only heteroplasmic indels are supported at the first site of an insertion.')
                to_replace_end = to_replace_start + len(allele[0])
                if ((to_replace_start + len(allele[0])) > len(to_replace)) | ((to_replace_start == 0) & ((to_replace_start + len(allele[0])) == len(to_replace))):
                    raise ValueError('Deletions should not be longer than the insertion they are contained within.')
                grch38_ref_holder.append(to_replace[0])
                grch38_alt_holder.append(to_replace[0:to_replace_start] + alt + to_replace[to_replace_end:])
                site_of_reference_holder.append(False)
        if len(set(grch38_ref_holder)) > 1:
            raise ValueError('There should be only a single grch38 reference allele per site.')
        else:
            grch38_ref = grch38_ref_holder[0]
        return {'grch38_ref': grch38_ref, 'grch38_alt': grch38_alt_holder, 'get_data_from_self_pos_0': site_of_reference_holder}

    result_holder = df.apply(replace_substring, axis=1)
    df['grch38_ref'] = result_holder.apply(lambda x: x['grch38_ref'])
    df['grch38_alt'] = result_holder.apply(lambda x: x['grch38_alt'])
    df['get_data_from_self_pos_0'] = result_holder.apply(lambda x: x['get_data_from_self_pos_0'])
    df['paired_allele_alt'] = df.apply(lambda x: list(zip(x.grch38_alt, x.allele_array, x.get_data_from_self_pos_0)),1)
    ht_new = hl.Table.from_pandas(df).explode('paired_allele_alt')
    ht_new = ht_new.select(locus = hl.locus('chrM', hl.int32(ht_new.incr_position), incr_reference),
                           alleles = [ht_new.grch38_ref, ht_new.paired_allele_alt[0]],
                           locus_self_ref = hl.locus('chrM', hl.int32(ht_new.site), self_reference),
                           alleles_self_ref = ht_new.paired_allele_alt[1],
                           position_0_source_of_alt_data = ht_new.paired_allele_alt[2])
    return ht_new.key_by('locus_self_ref','alleles_self_ref')


def produce_ref_deletions_table(mt, deletions):
    ht = mt.select_rows().select_entries('AD').entries().key_by('locus','alleles')
    exploded_dele = explode_indel(deletions, deletions.coord_t, deletions.nt, MTREF, 'sites')
    exploded_dele = exploded_dele.annotate(ref_alleles = [exploded_dele.ref_allele, exploded_dele.alt_allele])
    exploded_dele = exploded_dele.key_by('sites', 'ref_alleles')
    exploded_dele = exploded_dele.anti_join(ht).key_by('coord_t','ref_alleles')
    exploded_dele = exploded_dele.annotate(DP = ht[exploded_dele.key].AD).key_by()
    exploded_dele = exploded_dele.select(reference_position = exploded_dele.sites.position,
                                         original_deletion_position = exploded_dele.coord_t.position,
                                         ref_allele = exploded_dele.ref_alleles[0],
                                         alt_allele = exploded_dele.ref_alleles[1],
                                         ref_allele_depth = exploded_dele.DP[0],
                                         alt_allele_depth = exploded_dele.DP[1])
    return exploded_dele


def make_first_indel_site_ref(mt, mt_meta, matching_homoplasmies, reference_coordinate_snvs, self_ref, ref, exploded_insertions, allow_NONREF):
    mt = mt.persist()
    if mt.count_rows() == 0:
        return enforce_corrected_spec_in_uncorrected_mt(mt), add_missing_new_entry_fields(mt), 0
    ROW_R_fields, ROW_BOOL_fields, ROW_custom_fields, _, ENTRY_R_fields, _, _ = get_fields_for_flipping(mt_meta)

    reference_coordinate_snvs = reference_coordinate_snvs.annotate_rows(row_found=True)
    mt = mt.annotate_rows(allele_lens = hl.map(lambda x: x.length(), mt.alleles))
    mt = mt.annotate_rows(allele_idx_swap = mt.allele_lens.index(1)) # 0 indicates insertion
    mt = mt.annotate_rows(grch38_array = [mt.grch38_seq, mt.alleles[mt.allele_idx_swap]]).key_rows_by('locus_grch38', 'grch38_array')
    mt = mt.annotate_rows(flipped_row_data = reference_coordinate_snvs.rows()[mt.row_key])
    mt = mt.annotate_entries(flipped_entries = reference_coordinate_snvs[mt.row_key, mt.col_key]).key_rows_by('locus', 'alleles')

    if allow_NONREF:
        raise argparse.ArgumentError('ERROR: allow_NONREF is currnetly not supported for make_first_indel_site_ref().')
    if mt.aggregate_rows(~hl.agg.all(hl.is_defined(mt.allele_idx_swap))):
        raise ValueError('ERROR: when procesing indels to alter the first site to reference, an allele of length 1 must always be present.')
    if mt.aggregate_rows(~hl.agg.all(hl.is_defined(mt.locus_grch38) & hl.is_defined(mt.grch38_seq))):
        raise ValueError('ERROR: GRCh38 locus and sequence at every site for swap must be defined.')
    if mt.aggregate_rows(~hl.agg.all(hl.is_defined(mt.flipped_row_data))):
        raise ValueError('ERROR: Must successfully obtain row data from homoplasmy table.')
    # the current reference must not be the grch38 reference
    if mt.aggregate_rows(~hl.agg.all(mt.alleles[mt.allele_idx_swap] != mt.grch38_seq)):
        raise ValueError('ERROR: GRCh38 sequence matches self-reference single-site sequence at >1 site, which makes no sense given these have failed to Liftover and are in the first-site remapping pipeline.')
    # ensure all loci in mt are in matching_homoplasmies
    if mt.key_rows_by('locus').anti_join_rows(matching_homoplasmies.rows().key_by('locus')).count_rows() > 0:
        raise ValueError('ERROR: All loci in input mt should have corresponding homoplasmies.')

    # Prooduce new alleles
    mt = mt.annotate_rows(updated_alleles = hl.map(lambda x: x.replace('^[A-Z]{1}', mt.grch38_seq), mt.alleles))
    if mt.aggregate_rows(~hl.agg.all(hl.map(lambda x: x[0], mt.updated_alleles)[0] == hl.map(lambda x: x[0], mt.updated_alleles)[1])):
        raise ValueError('ERROR: changing first site of indel did not produce matching first site in resultant alleles.')
    mt = mt.annotate_rows(ref_end_position = get_ref_locus_end(mt.locus, mt.updated_alleles[0], exploded_insertions, self_ref, ref))
    #mt = mt.annotate_rows(ref_end_position = hl.liftover(hl.locus('chrM', mt.locus.position + mt.updated_alleles[0].length() - 1, reference_genome=self_ref), ref))
    if mt.aggregate_rows(~hl.agg.all(hl.is_defined(mt.ref_end_position))):
        raise ValueError('ERROR: the end position of all REF allele should be liftover-compatible.')
    mt = mt.annotate_rows(expected_ref = hl.get_sequence('chrM', mt.locus_grch38.position, before=0, after=mt.ref_end_position.position-mt.locus_grch38.position, reference_genome=ref))
    
    # Do swap
    mt_new = complex_swap_field(mt, ROW_BOOL_fields=ROW_BOOL_fields, ROW_R_fields=ROW_R_fields, 
                                ROW_custom_fields=ROW_custom_fields, ENTRY_R_fields=ENTRY_R_fields, 
                                flipped_idx=0, skip_locus_change=True)
    mt_new = mt_new.drop('grch38_array', 'allele_lens', 'ref_end_position')
    mt_new = add_filter(mt_new, 'SwapFirstAlleleIndel')

    # check if any deletions still need to have in-body SNPs fixed
    mt_needs_fixing = mt_new.filter_rows(mt_new.expected_ref != mt_new.alleles[0])
    if mt_needs_fixing.aggregate_rows(hl.agg.any(mt_needs_fixing.allele_idx_swap == 0)):
        raise ValueError('ERROR: any records that have a different reference allele than GRCh38 must be deletions.')
    mt_needs_fixing = mt_needs_fixing.drop('expected_ref', 'allele_idx_swap')
    mt_needs_fixing = mt_needs_fixing.select_rows(*get_nonkey_fields(mt_needs_fixing.rows()))
    n_for_in_body = mt_needs_fixing.count_rows()

    # change locus
    mt_new = mt_new.anti_join_rows(mt_needs_fixing.rows())
    mt_new = mt_new.key_rows_by().drop('locus').rename({'locus_grch38': 'locus'}).key_rows_by('locus','alleles')
    mt_new = mt_new.drop('allele_idx_swap', 'grch38_seq', 'expected_ref')

    return mt_new, mt_needs_fixing, n_for_in_body


def recode_deletion_allele(failed_candidates, piped_from_first_indel_site, exploded_insertions, reference, self_reference, 
                           allow_NONREF, allow_spanning=False, override_processed=False, override_swap=False, self_homoplasmies=None, mt_meta=None):
    """ This function does not change INFO/entry fields. See Notes on liftover field changes.md.
    override_processed True indicates that alleles has already been processed into GRCh38.
    """
    if failed_candidates is not None: 
        failed_candidates = failed_candidates.annotate_entries(OriginalSelfRefAlleles = failed_candidates.alleles,
                                                               SwappedFieldIDs = hl.missing('tstr'))
        mt = failed_candidates.union_rows(piped_from_first_indel_site)
    else:
        mt = piped_from_first_indel_site
    if mt.count_rows() == 0:
        return enforce_corrected_spec_in_uncorrected_mt(mt, skip_entry_add=True), mt, 0
    
    if allow_NONREF:
        raise argparse.ArgumentError('ERROR: allow_NONREF is currently not supported for recode_deletion_allele().')
    if mt.aggregate_rows(~hl.agg.all(is_deletion(mt))):
        raise ValueError('All records input to recode_deletion_allele must be deletions.')
    
    # rows without GRCh38 locus defined fail
    mt_failed = mt.filter_rows(~hl.is_defined(mt.locus_grch38))
    mt = mt.filter_rows(hl.is_defined(mt.locus_grch38))

    # rows containing a deletion that spans a reference insertion fail
    if not allow_spanning:
        mt_failed = mt_failed.union_rows(dele_spans_insertion(mt, exploded_insertions, self_reference, False))
        mt = mt.anti_join_rows(dele_spans_insertion(mt, exploded_insertions, self_reference, False).rows())

    # we now explicitly deal with the edge case of a deletion ending at the end of a reference insertion
    if override_processed:
        mt = mt.annotate_rows(locus_last_base_grch38 = hl.locus('chrM', mt.locus_grch38.position + hl.len(mt.alleles[0]) - 1, reference_genome=reference))
    else:
        mt = mt.annotate_rows(locus_last_base_grch38 = get_ref_locus_end(mt.locus, mt.alleles[0], exploded_insertions, self_reference, reference))
    cols_to_rm = ['locus_last_base_grch38']

    # rows without a defined grch38 position for the last base of the deletion fail
    mt_failed = mt_failed.union_rows(mt.filter_rows(~hl.is_defined(mt.locus_last_base_grch38)).drop(*cols_to_rm)).persist()
    mt = mt.filter_rows(hl.is_defined(mt.locus_last_base_grch38))

    mt = mt.annotate_rows(candidate_grch38_allele = hl.get_sequence('chrM', mt.locus_grch38.position, 
                                                                    before=0, after=mt.locus_last_base_grch38.position-mt.locus_grch38.position, 
                                                                    reference_genome=reference))
    cols_to_rm.append('candidate_grch38_allele')

    # rows with no defined candidate allele or one of size 1 fail
    mt_failed = mt_failed.union_rows(mt.filter_rows(~hl.is_defined(mt.candidate_grch38_allele) | (hl.len(mt.candidate_grch38_allele)<=1)).drop(*cols_to_rm))
    mt = mt.filter_rows(hl.is_defined(mt.candidate_grch38_allele) & (hl.len(mt.candidate_grch38_allele)>1))

    # rows where the first character of each allele are different fail
    same_first_chars = (mt.candidate_grch38_allele[0] == mt.alleles[0][0]) & (mt.candidate_grch38_allele[0] == mt.alleles[1])
    mt_failed = mt_failed.union_rows(mt.filter_rows(~same_first_chars).drop(*cols_to_rm))
    mt = mt.filter_rows(same_first_chars)

    # if candidate sequence is longer or shorter than original, this indicates indel was present
    num_indel_in_dele = mt.aggregate_rows(hl.agg.count_where(hl.len(mt.alleles[0]) != hl.len(mt.candidate_grch38_allele)))

    # store original alleles in OriginalSelfRefAlleles. No alleles have swapped info fields so SwappedFieldIDs remains.
    mt = mt.annotate_entries(OriginalSelfRefAlleles = hl.if_else(hl.is_defined(mt.OriginalSelfRefAlleles), mt.OriginalSelfRefAlleles, mt.alleles))
    mt = mt.annotate_rows(new_alleles = [mt.candidate_grch38_allele, mt.alleles[1]])
    mt = mt.key_rows_by().drop('locus','alleles').rename({'new_alleles': 'alleles', 'locus_grch38': 'locus'}).key_rows_by('locus','alleles')
    mt = add_filter(mt, 'ReplaceInternalBaseDeletion')
    mt = mt.drop(*cols_to_rm, 'grch38_seq').persist()

    # if the variant shares a start with any homoplasmy, we need to perform swap alleles
    if not override_swap:
        ref_homoplasmies = self_homoplasmies.annotate_rows(locus_ref = hl.liftover(self_homoplasmies.locus, reference)).key_rows_by('locus_ref')
        ref_homoplasmies_indel = ref_homoplasmies.filter_rows(~is_snv(ref_homoplasmies))
        ref_homoplasmies_indel = ref_homoplasmies_indel.annotate_rows(row_found = True)
        ref_homoplasmies_snv = ref_homoplasmies.filter_rows(is_snv(ref_homoplasmies))

        # there should not be any homoplasmic SNVs where REF matches the first base of this deletion -- this should already be taken care of
        mt_shares_homoplasmic_snv = mt.key_rows_by('locus')
        mt_shares_homoplasmic_snv = mt_shares_homoplasmic_snv.annotate_rows(mapped_alleles = ref_homoplasmies_snv.rows()[mt_shares_homoplasmic_snv.row_key].alleles)
        mt_shares_homoplasmic_snv = mt_shares_homoplasmic_snv.filter_rows(hl.is_defined(mt_shares_homoplasmic_snv.mapped_alleles))
        if mt_shares_homoplasmic_snv.aggregate_rows(~hl.agg.all(mt_shares_homoplasmic_snv.mapped_alleles[0][0] != mt_shares_homoplasmic_snv.alleles[1])):
            raise ValueError('ERROR: Found some homoplasmic SNVs where the self-ref REF allele matches the deletion ALT allele in recode_deletion_allele().')
        
        mt_for_swap = mt.key_rows_by('locus').semi_join_rows(ref_homoplasmies_indel.rows())
        mt_no_need = mt.key_rows_by('locus').anti_join_rows(ref_homoplasmies_indel.rows()).key_rows_by('locus','alleles')
        
        mt_for_swap = mt_for_swap.annotate_rows(flipped_row_data = ref_homoplasmies_indel.rows()[mt_for_swap.row_key])
        mt_for_swap = mt_for_swap.annotate_entries(flipped_entries = ref_homoplasmies_indel[mt_for_swap.row_key, mt_for_swap.col_key]).key_rows_by('locus', 'alleles')
        if mt_for_swap.aggregate_rows(~hl.agg.all(hl.is_defined(mt_for_swap.flipped_row_data.row_found))):
            raise ValueError('ERROR: all candidates of successfully Liftedover loci for flipping must have mapped.')
        if mt_for_swap.aggregate_entries(~hl.agg.all(~hl.is_defined(mt_for_swap.SwappedFieldIDs))):
            raise ValueError('ERROR: SwappedFieldIDs should not be defined for variants that are candidates for flipping in recode_deletion_allele().')

        ROW_R_fields, ROW_BOOL_fields, ROW_custom_fields, _, ENTRY_R_fields, _, _ = get_fields_for_flipping(mt_meta)
        mt_for_swap_swapped = complex_swap_field(mt_for_swap, ROW_BOOL_fields=ROW_BOOL_fields, ROW_R_fields=ROW_R_fields,
                                                 ROW_custom_fields=ROW_custom_fields, ENTRY_R_fields=ENTRY_R_fields, flipped_idx=1,
                                                 skip_locus_change=True, skip_allele_change=True, original_self_alleles='OriginalSelfRefAlleles')

        mt = mt_no_need.union_rows(mt_for_swap_swapped).persist()

    return mt, mt_failed, num_indel_in_dele


def complex_swap_field(mt, ROW_BOOL_fields, ROW_R_fields, ROW_custom_fields, ENTRY_R_fields, flipped_idx, skip_locus_change, original_self_alleles=None, skip_allele_change=False, keep_flipped_fields=False):
    """ This function performes field ComplexSwapField swaps. This swap is characterized by
    using the information from a mapped force-call allele to fill the REF position fields
    for a given allele.

    Expects that updated_alleles have already been computed for each bespoke application.
    flipped_idx indicates which data point is used from flipped_row_data and flipped_entries.
    """
    if skip_locus_change and skip_allele_change:
        # if both are enabled we are just changing INFO and ENTRY fields
        row_reqd = ['info', 'flipped_row_data', 'locus', 'alleles']
    elif skip_allele_change and not skip_locus_change:
        raise ValueError('ERROR: complex_swap_field does not support skipping just allele change. Either skip locus change only, or skip both to just change INFO/ENTRY.')
    else:
        row_reqd = ['info', 'flipped_row_data', 'locus', 'alleles', 'locus_grch38', 'updated_alleles']
    entry_reqd = ['flipped_entries']
    
    if not all([x in mt.row for x in row_reqd]) or not all([x in mt.entry for x in entry_reqd]):
        raise ValueError('ERROR: The following row fields are required for input to complex_swap_field: ' + \
                         ','.join(row_reqd) + ' and the following entry fields are required: ' + ','.join(entry_reqd))

    # ensure annotations were found for all rows
    if mt.aggregate_rows(~hl.agg.all(hl.is_defined(mt.flipped_row_data.row_found))):
        raise ValueError('ERROR: The appropriate reference data was not found for all loci to be repaired.')

    # edit rows
    def custom_flip_expr(expr, replace, orientation=None):
        """
        Updated on 2/5/22 to ignore orientation and always replace the reference position with
        reference data from "replace".
        """
        #return hl.if_else(hl.is_defined(expr), hl.if_else(orientation == 0, [replace[0], expr[1]], [expr[0], replace[0]]), expr)
        return hl.if_else(hl.is_defined(expr), [replace[flipped_idx], expr[1]], expr)


    info_str = mt.info.annotate(**{x: ~mt.info[x] for x in ROW_BOOL_fields})
    info_str = info_str.annotate(**{x: custom_flip_expr(info_str[x], mt.flipped_row_data.info[x]) for x in ROW_R_fields})
    ROW_custom_fields_dct = {}
    for k, v in ROW_custom_fields.items():
        new_val = hl.literal(re.sub('\\\\','',v)).join(custom_flip_expr(info_str[k].split(v), mt.flipped_row_data.info[k].split(v)))
        ROW_custom_fields_dct.update({k: new_val})
    info_str = info_str.annotate(**ROW_custom_fields_dct)
    mt_new = mt.annotate_rows(info = info_str).persist()

    # edit entries
    mt_new = mt_new.annotate_entries(**{x: custom_flip_expr(mt_new[x], mt_new.flipped_entries[x]) for x in ENTRY_R_fields})

    # apply new alleles
    if skip_locus_change and skip_allele_change:
        mt_new = mt_new.annotate_entries(OriginalSelfRefAlleles = hl.missing('array<str>') if original_self_alleles is None else mt_new[original_self_alleles])
    else:
        mt_new = mt_new.annotate_entries(OriginalSelfRefAlleles = mt_new.alleles if original_self_alleles is None else mt_new[original_self_alleles])
        mt_new = mt_new.key_rows_by().drop('alleles').rename({'updated_alleles':'alleles'}).key_rows_by('locus','alleles')
        if not skip_locus_change:
            mt_new = mt_new.key_rows_by().drop('locus').rename({'locus_grch38': 'locus'}).key_rows_by('locus','alleles')
    
    mt_new = mt_new.annotate_entries(SwappedFieldIDs = ','.join([*ROW_BOOL_fields, *ROW_R_fields, *list(ROW_custom_fields.keys()), *ENTRY_R_fields]))
    if not keep_flipped_fields:
        mt_new = mt_new.drop('flipped_entries', 'flipped_row_data')
    mt_new = add_filter(mt_new, 'ComplexSwapField')

    return mt_new


def fancy_flip_entries(mt, mt_meta, success_vcf, debug=False):
    """
    We modify the rows with 1-AF... annotations using information from other alleles at that locus.

    Changed on 2/6/22 to group by locus rather than locus, reference allele pair.
    """
    ENTRY_1minus_fields = get_fields_for_flipping(mt_meta)[5]

    # Produce a list of loci that need fancy flipping -- these are loci with > 1 site in which 1 entry was swapped using 1-x
    # NOTE We integrate with success VCF in case there are other alleles at the site
    # NOTE We define "same site" as records with the same locus and reference allele.
    mt = mt.annotate_rows(has_oneminus = hl.agg.any(mt.SwappedFieldIDs.matches(','.join(['1-' + x for x in ENTRY_1minus_fields]))))
    success_ht = success_vcf.entries().select(*ENTRY_1minus_fields)
    success_ht = success_ht.annotate(SwappedFieldIDs = '', has_oneminus=False)
    ht_for_flip = mt.entries().select(*ENTRY_1minus_fields, 'SwappedFieldIDs', 'has_oneminus')
    joint_ht_for_flip = success_ht.union(ht_for_flip)
    joint_ht_for_flip = joint_ht_for_flip.annotate(ref_allele = joint_ht_for_flip.alleles[0])
    per_locus_info = joint_ht_for_flip.group_by(joint_ht_for_flip.locus#, joint_ht_for_flip.ref_allele
                                     ).aggregate(n = hl.agg.count(),
                                                 n_oneminus = hl.agg.count_where(joint_ht_for_flip.has_oneminus))
    
    # Checks
    if per_locus_info.aggregate(hl.agg.any(per_locus_info.n_oneminus > 1)):
        raise ValueError('ERROR: No locus should have more than one allele with entries that were flipped.')
    
    # Locations to flip are those with >1 record at the location (if there is only 1, then the 1-AF is correct) and that has a site with 1-x having been done
    locs_to_fancy_flip = per_locus_info.filter((per_locus_info.n > 1) & (per_locus_info.n_oneminus > 0))
    
    if locs_to_fancy_flip.count() > 0:
        # pull in entry info for all records matching locus/ref allele with the sites that need to be fixed
        joint_ht_for_flip = joint_ht_for_flip.key_by('locus')#, 'ref_allele')
        joint_ht_used_for_flip = joint_ht_for_flip.semi_join(locs_to_fancy_flip)

        # Separate sites that will not be modified but will be used to fix the 1-x sites
        # This makes sense because the correction we are doing is to change the 1 in 1-x by (1-total), where
        # total is the sum of x across all other sites at that location (which were not flipped).
        jht_helper = joint_ht_used_for_flip.filter(~joint_ht_used_for_flip.has_oneminus)

        # Compute sums. Array agg won't work if array lengths are not equal.
        jht_sums = jht_helper.group_by(*jht_helper.key).aggregate(**{x + '_sum': hl.agg.array_sum(jht_helper[x]).map(lambda y: 1-y) for x in ENTRY_1minus_fields})

        # Add the new totals and re-invert the originally 1-x fields
        joint_ht_flipped = joint_ht_used_for_flip.filter(joint_ht_used_for_flip.has_oneminus)
        joint_ht_flipped = joint_ht_flipped.annotate(**{f'{x}_new_total': jht_sums[joint_ht_flipped.key][x + '_sum'] for x in ENTRY_1minus_fields})
        joint_ht_flipped = joint_ht_flipped.annotate(**{f'{x}_reinverted': joint_ht_flipped[x].map(lambda y: 1-y) for x in ENTRY_1minus_fields}).persist()

        # Checks
        tf_vec = [joint_ht_flipped.aggregate(~hl.agg.all(joint_ht_flipped[f'{x}_new_total'].length() == joint_ht_flipped[f'{x}_reinverted'].length())) for x in ENTRY_1minus_fields]
        if any(tf_vec):
            raise ValueError('ERROR: for all entry fields to invert and all sites, new total must have the same array length as the actual field.')
        
        # Produce new entry fields
        joint_ht_flipped = joint_ht_flipped.annotate(**{f'{x}_fixed': hl.zip(joint_ht_flipped[f'{x}_new_total'], joint_ht_flipped[f'{x}_reinverted']).map(lambda tup : tup[0]-tup[1]) for x in ENTRY_1minus_fields})
        joint_ht_flipped = joint_ht_flipped.key_by('locus', 'alleles')

        mt_noflip = mt.anti_join_rows(joint_ht_flipped)
        mt_flipped = mt.semi_join_rows(joint_ht_flipped)
        mt_flipped = mt_flipped.annotate_entries(**{f'{x}_fixed': joint_ht_flipped[mt_flipped.row_key][f'{x}_fixed'] for x in ENTRY_1minus_fields})

        # Checks        
        if mt_flipped.aggregate_rows(~hl.agg.all(mt_flipped.has_oneminus)):
            raise ValueError('ERROR: all records which get fancy flipped entries should have already been logged as flipped in SwappedFieldIDs.')
        tf_vec_same_len = [mt_flipped.aggregate_entries(~hl.agg.all(mt_flipped[x].length() == mt_flipped[f'{x}_fixed'].length())) for x in ENTRY_1minus_fields]
        if any(tf_vec_same_len):
            raise ValueError('ERROR: all entries must have the same array length as the fancy flipped version they will be replaced by.')
        tf_vec_any_none = [mt_flipped.aggregate_entries(~hl.agg.all(hl.is_defined(mt_flipped[f'{x}_fixed']))) for x in ENTRY_1minus_fields]
        if not debug and any(tf_vec_any_none):
            raise ValueError('ERROR: missingness in fancy flipped results is not allowed.')
        tf_vec_any_nonearray = [mt_flipped.aggregate_entries(~hl.agg.all(hl.all(mt_flipped[f'{x}_fixed'].map(hl.is_defined)))) for x in ENTRY_1minus_fields]
        if not debug and any(tf_vec_any_nonearray):
            raise ValueError('ERROR: missingness in array of fancy flipped results is not allowed.')
        mt_flipped = add_filter(mt_flipped, 'FancyFieldInversion')
        
        # clean up
        mt_flipped = mt_flipped.drop(*ENTRY_1minus_fields)
        mt_flipped = mt_flipped.rename({f'{x}_fixed': x for x in ENTRY_1minus_fields})
        n_fancy_flip = mt_flipped.count_rows()
        mt_joined = mt_noflip.union_rows(mt_flipped.select_entries(*mt_noflip.entry))
    else:
        mt_joined = mt
        n_fancy_flip = 0

    mt_joined = mt_joined.drop('has_oneminus').persist()
    return mt_joined, n_fancy_flip


def resolve_deletion_boundary_cases(mt, mt_meta, exploded_insertions, self_homoplasmies, insertions_table, ref, self_ref, allow_NONREF):
    """ Thus function only supports:
    - Spanning deletions where first base is the same as a reference insertion and last base is outside
    - Deletions where the first base is BEFORE the insetion and the last base is within or after

    All others are piped to failure.
    """
    mt = mt.persist()
    if mt.count_rows() == 0:
        return enforce_corrected_spec_in_uncorrected_mt(mt, skip_entry_add=False), add_missing_new_entry_fields(mt), 0, 0, 0, 0
    ROW_R_fields, ROW_BOOL_fields, ROW_custom_fields, _, ENTRY_R_fields, _, _ = get_fields_for_flipping(mt_meta)
    
    if mt.anti_join_rows(dele_spans_insertion(mt, exploded_insertions, self_ref, False, False).rows()).count_rows() > 0:
        raise ValueError('ERROR: All inputted variants to resolve_deletion_boundary_cases() must be deletions that span insertion boundaries.')

    mt_fixable_1 = dele_spans_insertion(mt, exploded_insertions, self_ref, False, True).annotate_rows(fixtype=1)
    mt_fixable_2_orig = dele_spans_insertion(mt, exploded_insertions, self_ref, False, False, True).annotate_rows(fixtype=2)
    self_homoplasmies = self_homoplasmies.annotate_rows(row_found=True)
    self_inserts = self_homoplasmies.filter_rows(is_deletion(self_homoplasmies))
    insertions_table = insertions_table.annotate(alleles = [insertions_table.alt_allele, insertions_table.ref_allele]).key_by('coord_s', 'alleles')
    self_inserts = self_inserts.annotate_rows(ns = insertions_table[self_inserts.row_key].ns)
    self_inserts = self_inserts.key_rows_by('locus')

    # 1 refers to spanning deletions where the first base is same as an insertion
    mt_fixable_1 = mt_fixable_1.key_rows_by('locus')
    mt_fixable_1 = mt_fixable_1.annotate_rows(flipped_row_data = self_inserts.rows()[mt_fixable_1.row_key])
    mt_fixable_1 = mt_fixable_1.annotate_entries(flipped_entries = self_inserts[mt_fixable_1.row_key, mt_fixable_1.col_key]).key_rows_by('locus', 'alleles')

    # 2 is deletions where first base is upstream and second base is within
    # we do not currently handle this type of spanning deletion where the first base overlaps a SNV, so fail this
    mt_fixable_2_orig = mt_fixable_2_orig.key_rows_by('locus'
                                        ).anti_join_rows(self_homoplasmies.filter_rows(is_snv(self_homoplasmies)).rows().key_by('locus')
                                        ).key_rows_by('locus','alleles')
    mt_fixable_2 = mt_fixable_2_orig.annotate_rows(delelen = hl.len(mt_fixable_2_orig.alleles[0]))
    mt_fixable_2 = mt_fixable_2.annotate_rows(internal_positions = hl.range(mt_fixable_2.locus.position, mt_fixable_2.locus.position+(mt_fixable_2.delelen-1)))
    mt_fixable_2 = mt_fixable_2.annotate_rows(internal_positions = hl.map(lambda x: hl.locus('chrM', x, self_ref), mt_fixable_2.internal_positions))
    mt_fixable_2 = mt_fixable_2.explode_rows('internal_positions').key_rows_by('internal_positions')
    mt_fixable_2 = mt_fixable_2.filter_rows(mt_fixable_2.locus != mt_fixable_2.internal_positions) # matches must be internal
    mt_fixable_2 = mt_fixable_2.annotate_rows(flipped_row_data = self_inserts.rows()[mt_fixable_2.row_key])
    mt_fixable_2 = mt_fixable_2.annotate_entries(flipped_entries = self_inserts[mt_fixable_2.row_key, mt_fixable_2.col_key]).key_rows_by('locus', 'alleles')
    mt_fixable_2 = mt_fixable_2.annotate_rows(idx_last = mt_fixable_2.internal_positions.position - mt_fixable_2.locus.position)
    mt_fixable_2 = mt_fixable_2.filter_rows(hl.is_defined(mt_fixable_2.flipped_row_data.row_found))
    
    mt_fixable = mt_fixable_1.union_rows(mt_fixable_2.drop('delelen', 'internal_positions', 'idx_last')).persist()
    mt_unfixable = mt.anti_join_rows(mt_fixable.rows())

    # flags specific to type 1
    if mt_fixable_1.aggregate_rows(~hl.agg.all((mt_fixable_1.flipped_row_data.ns > 0) & ((mt_fixable_1.flipped_row_data.ns+1) < mt_fixable_1.alleles[0].length()))):
        raise NotImplementedError('ERROR: All records for fixing must have ns > 0 and < length of the allele to repair.')
    if mt_fixable_1.aggregate_rows(~hl.agg.all(mt_fixable_1.locus_grch38 == mt_fixable_1.flipped_row_data.locus_grch38)):
        raise NotImplementedError('ERROR: annotated locus_grch38 must be the same between the force-called homoplasmy and variants to fix.')
    # ensure all loci in mt are in matching_homoplasmies
    if mt_fixable_1.key_rows_by('locus').anti_join_rows(self_inserts.rows()).count_rows() > 0:
        raise ValueError('ERROR: All loci for repair should have corresponding homoplasmic insertions.')
    
    # flags specific to type 2
    if mt_fixable_2_orig.anti_join_rows(mt_fixable_2.rows()).count_rows() > 0:
        raise ValueError('ERROR: some rows tagged for fixing due to spanning left-side insertion boundary were lost.')
    # idxlast must be within the length of the allele
    if mt_fixable_2.filter_rows((mt_fixable_2.idx_last + 1) >= (mt_fixable_2.delelen)).count_rows() > 0:
        raise ValueError('ERROR: computed truncation point for a heteroplasmic deletion spanning homoplasmic insertion is too long -- this should always be internal to the allele (not at or after the end).')
    if mt_fixable_2.filter_rows(mt_fixable_2.idx_last == 0).count_rows() > 0:
        raise ValueError('ERROR: computed truncation point for a heteroplasmic deletion spanning homoplasmic insertion is too short -- this should always be internal to the allele (not at the start).')

    # flags for all
    if mt_fixable.rows().count() != mt_fixable.rows().distinct().count():
        raise ValueError('ERROR: MT for fixing has duplicate entries. Check dele_spans_insertion() and/or mapping to mt_fixable_2.')
    if allow_NONREF:
        raise argparse.ArgumentError('ERROR: allow_NONREF is currnetly not supported for make_first_indel_site_ref().')
    if mt_fixable.aggregate_rows(~hl.agg.all(hl.literal(NONINDEL).contains(mt_fixable.alleles[1]))):
        raise NotImplementedError('ERROR: resolve_deletion_boundary_cases() can only handle deletions.')
    if mt_fixable.aggregate_rows(~hl.agg.all(hl.is_defined(mt_fixable.flipped_row_data.ns))):
        raise NotImplementedError('ERROR: All records for fixing must have available ns.')
    if mt_fixable.aggregate_rows(~hl.agg.all(hl.is_defined(mt_fixable.locus_grch38) & hl.is_defined(mt_fixable.grch38_seq))):
        raise ValueError('ERROR: GRCh38 locus and sequence at every site for swap must be defined.')
    # alleles should not be the same as respective homoplasmies
    if mt_fixable.aggregate_rows(~hl.agg.all(mt_fixable.alleles != mt_fixable.flipped_row_data.alleles)):
        raise ValueError('ERROR: Alleles to fix should not be force-call reference alleles.')
    # the current reference must not be the grch38 reference
    if mt_fixable.aggregate_rows(~hl.agg.all((mt_fixable.alleles[0] != mt_fixable.grch38_seq) & (mt_fixable.alleles[1] == mt_fixable.grch38_seq))):
        raise ValueError('ERROR: GRCh38 sequence matches self-reference single-site sequence at >1 site and/or does not match self-ref ALT at >1 site.')

    ### fix type 1
    mt_fixable_1 = mt_fixable_1.annotate_rows(region_removal = mt_fixable_1.alleles[0][0:hl.int32(mt_fixable_1.flipped_row_data.ns+1)])

    if mt_fixable_1.aggregate_rows(~hl.agg.all(mt_fixable_1.region_removal == mt_fixable_1.flipped_row_data.alleles[0])):
        raise ValueError('ERROR: Computed region for replacement must be identical to the self REF allele of the corresponding homoplasmy.')

    # Produce new alleles
    mt_fixable_1 = mt_fixable_1.annotate_rows(updated_alleles = [mt_fixable_1.flipped_row_data.alleles[1] + mt_fixable_1.alleles[0][hl.int32(mt_fixable_1.flipped_row_data.ns+1):], mt_fixable_1.alleles[1]])
    if mt_fixable_1.aggregate_rows(~hl.agg.all(hl.map(lambda x: x[0], mt_fixable_1.updated_alleles)[0] == hl.map(lambda x: x[0], mt_fixable_1.updated_alleles)[1])):
        raise ValueError('ERROR: changing first site of indel did not produce matching first site in resultant alleles.')
    mt_fixable_1 = mt_fixable_1.annotate_rows(ref_end_position = get_ref_locus_end(mt_fixable_1.locus, mt_fixable_1.alleles[0], exploded_insertions, self_ref, ref))
    if mt_fixable_1.aggregate_rows(~hl.agg.all(hl.is_defined(mt_fixable_1.ref_end_position))):
        raise ValueError('ERROR: the end position of all REF allele should be liftover-compatible.')
    mt_fixable_1 = mt_fixable_1.annotate_rows(expected_ref = hl.get_sequence('chrM', mt_fixable_1.locus_grch38.position, before=0, after=mt_fixable_1.ref_end_position.position-mt_fixable_1.locus_grch38.position, reference_genome=ref))
    mt_fixable_1 = mt_fixable_1.persist()

    # Swap first allele
    # 220615 we now skip allele change here to allow recode_deletion_allele to work if there are indels inside the RHS overhang
    mt_new = complex_swap_field(mt_fixable_1, ROW_BOOL_fields=ROW_BOOL_fields, ROW_R_fields=ROW_R_fields, 
                                ROW_custom_fields=ROW_custom_fields, ENTRY_R_fields=ENTRY_R_fields, 
                                flipped_idx=1, skip_locus_change=True, skip_allele_change=True)
    mt_new = mt_new.drop('ref_end_position', 'region_removal')
    mt_new = add_filter(mt_new, 'DeletionSpannedHomoplasmicInsertion')

    # For records that need further modification send to recode_deletion_allele()
    mt_for_dele_recode = mt_new.filter_rows(mt_new.updated_alleles[0] != mt_new.expected_ref)
    nrow_for_recode = mt_for_dele_recode.count_rows()
    mt_recoded, mt_failed, num_indel_in_dele = recode_deletion_allele(None, mt_for_dele_recode, exploded_insertions, ref, self_ref, allow_NONREF, 
                                                                      allow_spanning=True, override_processed=False, override_swap=True, 
                                                                      self_homoplasmies=self_homoplasmies, mt_meta=mt_meta)
    if mt_recoded.aggregate_rows(~hl.agg.all(mt_recoded.alleles[0] == mt_recoded.expected_ref)):
        raise ValueError('ERROR: expected reference allele not generated as expected during second round recode_deletion_allele.')
    
    mt_fixed = mt_new.filter_rows(mt_new.updated_alleles[0] == mt_new.expected_ref)
    mt_fixed = mt_fixed.annotate_entries(OriginalSelfRefAlleles = mt_fixed.alleles)
    mt_fixed = mt_fixed.key_rows_by().drop('locus', 'alleles').rename({'updated_alleles':'alleles','locus_grch38': 'locus'}).key_rows_by('locus','alleles')
    mt_fixed = mt_fixed.drop('grch38_seq', 'expected_ref')
    mt_fixed = mt_fixed.union_rows(mt_recoded.drop('expected_ref', 'updated_alleles'))

    # the resultant allele must be reference – if it doesnn't, throw an error!
    mt_fixed_test = mt_fixed.annotate_rows(expected = hl.get_sequence('chrM',mt_fixed.locus.position, reference_genome=ref, after=hl.len(mt_fixed.alleles[0])-1))
    if mt_fixed_test.filter_rows(mt_fixed_test.expected != mt_fixed_test.alleles[0]).count_rows() > 0:
        raise ValueError('ERROR: produced allele from type 1 spanning deletion rescue does not match the reference.')

    # Fix type 2. Note that get_ref_locus_end() can't work on the original allele, because the ending point of the deletion
    # is inside an insertion as defined. However it does work on the truncated allele, which is at base 1 of the deletion.
    mt_new_2 = mt_fixable_2.annotate_rows(updated_alleles = [mt_fixable_2.alleles[0][0:mt_fixable_2.idx_last+1], mt_fixable_2.alleles[1]]).persist()
    mt_new_2 = mt_new_2.annotate_entries(OriginalSelfRefAlleles = mt_new_2.alleles, SwappedFieldIDs=hl.missing('tstr'))
    mt_new_2 = mt_new_2.annotate_rows(ref_end_position = get_ref_locus_end(mt_new_2.locus, mt_new_2.updated_alleles[0], exploded_insertions, self_ref, ref))
    if mt_new_2.aggregate_rows(~hl.agg.all(hl.is_defined(mt_new_2.ref_end_position))):
        raise ValueError('ERROR: the end position of all REF allele should be liftover-compatible.')
    mt_new_2 = mt_new_2.annotate_rows(expected_ref = hl.get_sequence('chrM', mt_new_2.locus_grch38.position, before=0, after=mt_new_2.ref_end_position.position-mt_new_2.locus_grch38.position, reference_genome=ref)).drop('ref_end_position')
    #mt_new_2 = mt_new_2.annotate_rows(expected_ref = hl.get_sequence('chrM',mt_new_2.locus_grch38.position, reference_genome=ref, after=hl.len(mt_new_2.updated_alleles[0])-1))
    mt_new_2 = mt_new_2.key_rows_by().drop('alleles').rename({'updated_alleles':'alleles'}).key_rows_by('locus','alleles')
    mt_new_2 = add_filter(mt_new_2, 'DeletionSpannedHomoplasmicInsertion')

    # For records that need further modification send to recode_deletion_allele()
    mt_for_dele_recode_2 = mt_new_2.filter_rows(mt_new_2.alleles[0] != mt_new_2.expected_ref).drop('flipped_entries', 'flipped_row_data')
    nrow_for_recode2 = mt_for_dele_recode_2.count_rows()
    mt_recoded_2, mt_failed_2, num_indel_in_dele_2 = recode_deletion_allele(None, mt_for_dele_recode_2, exploded_insertions, ref, self_ref, allow_NONREF, 
                                                                            allow_spanning=True, override_processed=False, override_swap=False, 
                                                                            self_homoplasmies=self_homoplasmies, mt_meta=mt_meta)
    if mt_recoded_2.aggregate_rows(~hl.agg.all(mt_recoded_2.alleles[0] == mt_recoded_2.expected_ref)):
        raise ValueError('ERROR: expected reference allele not generated as expected during second round recode_deletion_allele.')

    mt_fixed_2 = mt_new_2.filter_rows(mt_new_2.alleles[0] == mt_new_2.expected_ref)
    mt_fixed_2 = mt_fixed_2.key_rows_by().drop('locus').rename({'locus_grch38': 'locus'}).key_rows_by('locus','alleles')
    mt_fixed_2 = mt_fixed_2.drop('grch38_seq')
    mt_fixed_2 = mt_fixed_2.drop('flipped_entries', 'flipped_row_data')
    mt_fixed_2 = mt_fixed_2.union_rows(mt_recoded_2.key_rows_by().select_rows(*[x for x in mt_fixed_2.row]).key_rows_by('locus','alleles')).drop('expected_ref')
    mt_fixed_2 = mt_fixed_2.drop('internal_positions', 'delelen', 'idx_last')
    
    # the resultant allele must be reference - if it doesn't, throw a warning and fail the allele
    # this likely occurs if there is a weird set of reference variants. For example, 431-BG01977.
    # We see ACT > A as reference -> makes the sequence ACTTCC -> ATCC; we see ATC > A which overlaps TCC > T.
    # Very complex; just fail.
    mt_fixed_2_test = mt_fixed_2.annotate_rows(expected = hl.get_sequence('chrM',mt_fixed_2.locus.position, reference_genome=ref, after=hl.len(mt_fixed_2.alleles[0])-1))
    if mt_fixed_2_test.filter_rows(mt_fixed_2_test.expected != mt_fixed_2_test.alleles[0]).count_rows() > 0:
        raise ValueError('ERROR: produced allele from type 2 spanning deletion rescue does not match the reference.')

    mt_fixed_final = mt_fixed.union_rows(mt_fixed_2).drop('fixtype').persist()

    mt_unfixable = add_missing_new_entry_fields(mt_unfixable)
    mt_failed = mt_failed.drop('expected_ref', 'fixtype', 'updated_alleles')
    mt_failed = mt_failed.select_rows(*[x for x in mt_failed.row if x not in mt_failed.row_key])
    mt_failed_2 = mt_failed_2.drop('expected_ref', 'fixtype')
    mt_failed_2 = mt_failed_2.select_rows(*[x for x in mt_unfixable.row if x not in mt_unfixable.row_key])
    mt_failed_final = mt_unfixable.union_rows(mt_failed, mt_failed_2).persist()

    if mt_fixed_final.count_rows() + mt_failed_final.count_rows() != mt.count_rows():
        raise ValueError('ERROR: in recode_deletion_allele(), failed + success rescued variants must be the same as the inputted variant set.')
    new_selfref_alleles = mt_fixed_final.OriginalSelfRefAlleles.collect()
    input_selfref_alleles = mt_fixable.alleles.collect()
    if any([x not in new_selfref_alleles for x in input_selfref_alleles]) & (len(new_selfref_alleles) != len(input_selfref_alleles)):
        raise ValueError('ERROR: it appears that some alleles were lost during the fixing process in resolve_deletion_boundary_cases() from the set of fixable alleles.')

    return mt_fixed_final, mt_failed_final, num_indel_in_dele+num_indel_in_dele_2, nrow_for_recode+nrow_for_recode2, mt_fixed.count_rows(), mt_fixed_2.count_rows()


def flip_success_fields(mt_success, mt_meta, repaired_mt, insertions, deletions):
    """
    Insertions found to map to the start of reference insertions should have already been injected.
    """
    ROW_R_fields, ROW_BOOL_fields, ROW_custom_fields, _, ENTRY_R_fields, _, _ = get_fields_for_flipping(mt_meta)
    insertions = insertions.annotate(alleles = [insertions.ref_allele, insertions.alt_allele], idx = 0).key_by('coord_t', 'alleles')
    deletions = deletions.annotate(alleles = [deletions.ref_allele, deletions.alt_allele], idx = 1).key_by('coord_t', 'alleles').drop('self_ref_alleles')
    ref_indels = insertions.union(deletions)
    if ref_indels.anti_join(repaired_mt.rows()).count() > 0:
        raise ValueError('ERROR: All predicted reference indels should be appropriately flipped in the repaired MT.')
    
    filtered_repaired_mt = repaired_mt.semi_join_rows(ref_indels)
    filtered_repaired_mt = filtered_repaired_mt.annotate_rows(idx = ref_indels[filtered_repaired_mt.row_key].idx, row_found = True).key_rows_by('locus')
    if filtered_repaired_mt.distinct_by_row().count_rows() != filtered_repaired_mt.count_rows():
        raise ValueError('ERROR: there should not be any duplicated loci in the repaired MT records used to flip success VCF.')
    mt_success_shared_locus = mt_success.key_rows_by('locus').semi_join_rows(filtered_repaired_mt.rows())
    mt_success_shared_locus = mt_success_shared_locus.annotate_rows(flipped_row_data = filtered_repaired_mt.rows()[mt_success_shared_locus.row_key])
    mt_success_shared_locus = mt_success_shared_locus.annotate_entries(flipped_entries = filtered_repaired_mt[mt_success_shared_locus.row_key, mt_success_shared_locus.col_key])
    mt_success_shared_locus = mt_success_shared_locus.key_rows_by('locus','alleles').persist()
    
    # all success should have found a matched indel
    if mt_success_shared_locus.aggregate_rows(~hl.agg.all(hl.is_defined(mt_success_shared_locus.flipped_row_data.row_found))):
        raise ValueError('ERROR: all candidates of successfully Liftedover loci for flipping must have mapped.')
    if mt_success_shared_locus.aggregate_rows(~hl.agg.all(hl.is_defined(mt_success_shared_locus.flipped_row_data.idx))):
        raise ValueError('ERROR: all candidates of successfully Liftedover loci for flipping must have mapped.')
    # all success should be SNVs
    # insepcting LY3013 shows that deletions sharing the same locus, or insertions & SNVs, show the same REF measurements
    # V911 shows that deletions sharing the same locus as SNVs also share the same REF quantities
    if mt_success_shared_locus.aggregate_rows(~hl.agg.all(is_snv(mt_success_shared_locus))):
        raise ValueError('ERROR: candidates of successfully Liftedover loci for flipping must be SNVs.')
    
    # based on indel status the single allele should be the same as the REF
    # NOTE this check only works when the REF allele in success VCF is a single base
    if mt_success_shared_locus.aggregate_rows(~hl.agg.all(mt_success_shared_locus.alleles[0] == mt_success_shared_locus.flipped_row_data.alleles[mt_success_shared_locus.flipped_row_data.idx])):
        raise ValueError('ERROR: the "single allele" of the reference indel should be the same as the success VCF REF allele.')
    # success REF allele should be the first letter of the REF allele attached
    # NOTE this check only works when the REF allele in success VCF is a single base
    if mt_success_shared_locus.aggregate_rows(~hl.agg.all(mt_success_shared_locus.alleles[0] == mt_success_shared_locus.flipped_row_data.alleles[0][0])):
        raise ValueError('ERROR: the first letter of the fixed MT ref allele should be the same as the successfully lifted VCF REF allele.')

    # Analysis for SNVs
    mt_success_shared_fixed = complex_swap_field(mt_success_shared_locus, ROW_BOOL_fields=ROW_BOOL_fields, ROW_R_fields=ROW_R_fields, 
                                                 ROW_custom_fields=ROW_custom_fields, ENTRY_R_fields=ENTRY_R_fields, flipped_idx=0, 
                                                 skip_locus_change=True, skip_allele_change=True)
    mt_success_shared_fixed = add_filter(mt_success_shared_fixed, 'LiftoverSuccessEntrySwap')
    n_success_flipped = mt_success_shared_fixed.count_rows()
    
    mt_success_nomod = mt_success.key_rows_by('locus').anti_join_rows(filtered_repaired_mt.rows()).key_rows_by('locus','alleles')
    mt_success_nomod = add_missing_new_entry_fields(mt_success_nomod)
    mt_success_final = mt_success_shared_fixed.union_rows(mt_success_nomod).persist()

    # make sure all previous keys are still there
    if mt_success_final.count_rows() !=  mt_success.count_rows():
        raise ValueError('ERROR: when flipping fields some records were lost/gained from the success VCF.')
    if (mt_success_final.anti_join_rows(mt_success.rows()).count_rows() != 0) | (mt_success.anti_join_rows(mt_success_final.rows()).count_rows() != 0):
        raise ValueError('ERROR: when flipping fields some records were lost/gained from the success VCF.')

    return mt_success_final, n_success_flipped


def left_align_and_normalize(mt, mt_meta, reference, reference_fasta_path):
    """
    Here we use bcftools norm to normalize and left-align variants. We have tested this extensively:

    - MH0126875
    A C C C C C /T T A C C T/ C C C A T <- T > TTACCT
    A C C /C C C T T A/ C C T C C C A T <- C > CCCTTA (left shifted)

    T A C C /C A C/ C C T T <- CAC > C
    T A C /C C A/ C C C T T <- CCA > C (left shifted)

    - V40776
    T C A A A A C C /C C/ C T C <- CC > C
    T C A A A /A C/ C C C C T C <- AC > A (left shifted)

    T C A A A A C C /C C/ C C T C <- C > CC
    T C A A A /A C/ C C C C C T C <- A > AC (left shifted)

    T A C C /C A C/ C C T T <- CAC > C
    T A C /C C A/ C C C T T <- CCA > C (left shifted)

    - NWD992111
    A C C C C C /T G C C T/ C C C C A T G <- T > TGCCT
    A C C /C C C T G/ C C T C C C C A T G <- C > CCCTG (left shifted)

    - 431-BG01977
    C A T A C T /T C C T T/ A C T A <- T > TCCTT
    C A T /A C T T C/ C T T A C T A <- A > ACTTC (left shifted)
    
    C A T A C T /T C T/ A C T A <- T > TCT
    C A T A C /T T C/ T A C T A <- T > TTC (left shifted)
    
    All of these seem appropriate.
    """
    hl.export_vcf(mt, 'tmp.vcf', metadata=mt_meta)
    res = subprocess.run(["bcftools", 'norm', '-f', reference_fasta_path, 'tmp.vcf', '-o', 'tmp_fixed.vcf'], 
                         stderr=subprocess.PIPE, text=True, check=True)

    if res.returncode != 0:
        raise subprocess.SubprocessError('ERROR: bcftools returned with non-zero exit code.')
    bcftools_split_output = [x for x in res.stderr.split('\n') if not re.search('\\[E::bcf_hdr_parse_line\\] Could not parse the header line: \\"##', x)]
    if len(bcftools_split_output) != 2:
        raise subprocess.SubprocessError('ERROR: bcftools stderr does not match expected format.')
    bcftools_output = re.search('^Lines   total/split/realigned/skipped:\t([0-9]+)/([0-9]+)/([0-9]+)/([0-9]+)$', bcftools_split_output[0])
    if not bcftools_output:
        raise subprocess.SubprocessError('ERROR: bcftools output did not match expected format.')
    mt_new, mt_meta_new = read_mito_vcf('tmp_fixed.vcf', reference)
    if int(bcftools_output[1]) != mt.count_rows():
        raise subprocess.SubprocessError('ERROR: bcftools did not count the same number of variants as was in input mt.')
    if mt_new.count_rows() != mt.count_rows():
        raise subprocess.SubprocessError('ERROR: bcftools appears to have lost/gained some variants.')
    if int(bcftools_output[2]) != 0:
        raise subprocess.SubprocessError('ERROR: bcftools split an allele. This should not happen.')
    if int(bcftools_output[3]) != mt_new.anti_join_rows(mt.rows()).count_rows():
        raise subprocess.SubprocessError('ERROR: the number of variants realigned by bcftools should be the same as the alleles that appear changed.')
    if int(bcftools_output[4]) != 0:
        raise subprocess.SubprocessError('ERROR: bcftools skipped alleles. This should not happen.')
    
    mt_preserved = mt_new.semi_join_rows(mt.rows())
    mt_modified = mt_new.anti_join_rows(mt.rows())
    mt_modified = add_filter(mt_modified, 'LeftShiftedIndel')
    mt_out = mt_preserved.union_rows(mt_modified).persist()
    return mt_out, mt_meta_new, int(bcftools_output[3]), res.stderr


def rescue_duplicated_alleles(mt, success_mt, failed_mt, meta, failed_meta, self_ref):
    # Here we filter out variants that are duplicate in mt and success mt.
    # Currently we do not modify the success mt.
    dupe_vars = mt.semi_join_rows(success_mt.rows())
    if dupe_vars.count_rows() > 0:
        # repair left shifted variants (but not homoplasmies, because these shouldn't be left-shifted and duplicated)
        dupe_left = dupe_vars.filter_rows(dupe_vars.filters.contains('LeftShiftedIndel') & ~dupe_vars.filters.contains('ForceCalledHomoplasmy'))
        dupe_vars = dupe_vars.anti_join_rows(dupe_left.rows())

        # can insert additional repair steps here on updated dupe_vars

        # remove repaired vars from mt, add to failed mt
        dupe_to_transfer = dupe_left # can union_rows more tables here if we transfer other causes of dupes
        n_failed = dupe_to_transfer.count_rows()
        mt = mt.anti_join_rows(dupe_to_transfer.rows())
        dupe_to_transfer = dupe_to_transfer.key_rows_by()
        dupe_to_transfer = dupe_to_transfer.annotate_rows(locus = hl.liftover(dupe_to_transfer.locus, self_ref)).key_rows_by('locus','alleles')
        dupe_to_transfer = add_filter(dupe_to_transfer, 'FailedDuplicateVariant')
        failed_mt = failed_mt.union_rows(dupe_to_transfer)

        # augment the metadata
        filters_to_add = list(set([item for x in dupe_to_transfer.filters.collect() for item in list(x) if item not in list(failed_meta['filter'].keys())]))
        failed_meta['filter'].update({x: meta['filter'][x] for x in filters_to_add})
    else:
        n_failed = 0
    return mt, failed_mt, failed_meta, n_failed


def fix_insertion_overlapping_ref_dele(mt, mt_meta, flipped_deletions, flipped_insertions, insertions, deletions, rely_on_left_shift=True):
    """
    Here we deal with insertions that overlap reference deletions.
    NOTE This currently deals with insertions must:
    - be mapped to a reference deletion (e.g., share a locus)
    - share the same first base as the reference deletion
    - if the reference deletion has a REF allele of length N, the first N bases of the ALT allele of the insertion should be identical to the reference deletion REF allele
    For example: REF HOM GCA > G and SELF HET G > GCACA -> GCA > GCACA -> G > GCA
    All others fail.
    """
    ROW_R_fields, ROW_BOOL_fields, ROW_custom_fields, _, ENTRY_R_fields, _, _ = get_fields_for_flipping(mt_meta)
    insertions = insertions.annotate(alleles = [insertions.ref_allele, insertions.alt_allele], idx = 0).key_by('coord_t', 'alleles')
    deletions = deletions.annotate(alleles = [deletions.ref_allele, deletions.alt_allele], idx = 1).key_by('coord_t', 'alleles').drop('self_ref_alleles')
    ref_indels = insertions.union(deletions).key_by('coord_s')

    mt_ins = mt.annotate_rows(tf_ins = is_insertion(mt)).key_rows_by('locus')
    mt_ins = mt_ins.annotate_rows(reference_indel = ref_indels[mt_ins.row_key].alleles)
    mt_ins = mt_ins.annotate_rows(tf_ref_indel = hl.is_defined(mt_ins.reference_indel))
    mt_ins = mt_ins.annotate_rows(tf_ref_dele = hl.literal(NONINDEL).contains(mt_ins.reference_indel[1]) & \
                                                ~hl.literal(NONINDEL).contains(mt_ins.reference_indel[0]))
    mt_ins = mt_ins.annotate_rows(tf_shares_start = mt_ins.alleles[1].startswith(mt_ins.reference_indel[0]) | \
                                                    mt_ins.reference_indel[0].startswith(mt_ins.alleles[1]))
    mt_ins = mt_ins.annotate_rows(tf_here = mt_ins.tf_ins & mt_ins.tf_ref_indel & mt_ins.tf_ref_dele & mt_ins.tf_shares_start)
    mt_ins = mt_ins.drop('tf_ins','reference_indel','tf_ref_indel','tf_shares_start')
    
    mt_to_process = mt_ins.filter_rows(mt_ins.tf_here).drop('tf_here')

    mt_repaired_indelmap = flipped_deletions.semi_join_rows(ref_indels.key_by('coord_t', 'alleles'))
    if mt_repaired_indelmap.distinct_by_row().count_rows() != mt_repaired_indelmap.count_rows():
        raise ValueError('ERROR: there should not be any duplicated loci in the repaired MT records used.')
    mt_to_process = mt_to_process.key_rows_by('locus_grch38')
    flipped_deletions_proc = flipped_deletions.annotate_rows(idx = ref_indels.key_by('coord_t','alleles')[flipped_deletions.row_key].idx, row_found=True).key_rows_by('locus')
    mt_to_process = mt_to_process.annotate_rows(flipped_row_data = flipped_deletions_proc.rows()[mt_to_process.row_key])
    mt_to_process = mt_to_process.annotate_entries(flipped_entries = flipped_deletions_proc[mt_to_process.row_key, mt_to_process.col_key])
    mt_to_process = mt_to_process.key_rows_by('locus','alleles')
    mt_to_process = mt_to_process.annotate_rows(mapped_alleles = mt_to_process.flipped_row_data.alleles).persist()

    # Checks
    # all success should have found a matched indel
    if mt_to_process.aggregate_rows(~hl.agg.all(hl.is_defined(mt_to_process.flipped_row_data.row_found))):
        raise ValueError('ERROR: even though we filtered to deletions, there are some alleles for modification that do not map to reference deletions.')
    if mt_to_process.aggregate_rows(~hl.agg.all(hl.is_defined(mt_to_process.flipped_row_data.idx))):
        raise ValueError('ERROR: even though we filtered to deletions, there are some alleles for modification that do not map to reference deletions.')
    if mt_to_process.aggregate_rows(~hl.agg.all(hl.is_defined(mt_to_process.mapped_alleles))):
        raise ValueError('ERROR: even though we filtered to deletions, there are some alleles for modification that do not map to reference deletions.')
    # all variants should be insertions
    if mt_to_process.aggregate_rows(~hl.agg.all(is_insertion(mt_to_process))):
        raise ValueError('ERROR: candidates of successfully Liftedover loci for flipping must be insertions.')

    # if the shared locus is an insertion, it must be mapped to a reference deletion
    if mt_to_process.aggregate_rows(~hl.agg.all(mt_to_process.tf_ref_dele)):
        raise ValueError('ERROR: all candidates of successfully Liftedover loci for flipping which are insertions must be mapped to a reference deletion.')
    # if the shared locus is an insertion, it must share the same first base as the reference deletion
    if mt_to_process.aggregate_rows(~hl.agg.all(mt_to_process.alleles[0][0] == mt_to_process.mapped_alleles[0][0])):
        raise ValueError('ERROR: all candidates for correction must share the same first base name reference deletion.')
    # if the shared locus is an insertion, the reference deletion REF must be the first n bases of the insertion ALT
    if mt_to_process.aggregate_rows(~hl.agg.all(mt_to_process.alleles[1].startswith(mt_to_process.mapped_alleles[0]) |
                                                mt_to_process.mapped_alleles[0].startswith(mt_to_process.alleles[1]))):
        raise ValueError('ERROR: all candidates for correction must have the same first n bases as the REF homoplasmic deletion (where n is the length of the REF homoplasmic deletion) or vice versa.')
    
    # based on indel status the single allele should be the same as the REF
    # NOTE this check only works when the REF allele in success VCF is a single base
    if mt_to_process.aggregate_rows(~hl.agg.all(mt_to_process.alleles[0] == mt_to_process.flipped_row_data.alleles[mt_to_process.flipped_row_data.idx])):
        raise ValueError('ERROR: the "single allele" of the reference indel should be the same as the success VCF REF allele.')
    # success REF allele should be the first letter of the REF allele attached
    # NOTE this check only works when the REF allele in success VCF is a single base
    if mt_to_process.aggregate_rows(~hl.agg.all(mt_to_process.alleles[0] == mt_to_process.flipped_row_data.alleles[0][0])):
        raise ValueError('ERROR: the first letter of the fixed MT ref allele should be the same as the successfully lifted VCF REF allele.')

    # Analysis for insertions
    # NOTE we have assumed that the ref allele is a substring
    mt_fixed = mt_to_process.annotate_rows(locus_self=mt_to_process.locus)
    if rely_on_left_shift:
        mt_fixed = mt_fixed.annotate_rows(updated_alleles=[mt_fixed.mapped_alleles[0], mt_fixed.alleles[1]])

        if mt_fixed.aggregate_rows(~hl.agg.all((hl.len(mt_fixed.updated_alleles[0]) > 1))):
            raise ValueError('ERROR: the new REF allele should be greater than length 1.')
        if mt_fixed.aggregate_rows(~hl.agg.all(mt_fixed.alleles[1].startswith(mt_fixed.alleles[0]))):
            raise ValueError('ERROR: the new REF allele should be a subset of the ALT allele.')
    else:
        mt_fixed = mt_fixed.annotate_rows(updated_alleles_pre=[mt_fixed.mapped_alleles[0], mt_fixed.alleles[1]])
        mt_fixed = mt_fixed.annotate_rows(updated_alleles_1=[mt_fixed.updated_alleles_pre[0].replace('^.(?<=' + mt_fixed.updated_alleles_pre[1][0] + ')' + \
                                                                                                     mt_fixed.updated_alleles_pre[1][1:], mt_fixed.updated_alleles_pre[1][0]),
                                                             mt_fixed.updated_alleles_pre[1][0]],
                                          updated_alleles_2=[mt_fixed.updated_alleles_pre[0][0], 
                                                             mt_fixed.updated_alleles_pre[1].replace('^.(?<=' + mt_fixed.updated_alleles_pre[0][0] + ')' + \
                                                                                                     mt_fixed.updated_alleles_pre[0][1:], mt_fixed.updated_alleles_pre[0][0])])
        mt_fixed = mt_fixed.annotate_rows(idx_single=hl.if_else(hl.len(mt_fixed.updated_alleles_pre[0]) > hl.len(mt_fixed.updated_alleles_pre[1]), 1, 0),
                                          idx_mult=hl.if_else(hl.len(mt_fixed.updated_alleles_pre[0]) > hl.len(mt_fixed.updated_alleles_pre[1]), 0, 1),
                                          updated_alleles=hl.if_else(hl.len(mt_fixed.updated_alleles_pre[0]) > hl.len(mt_fixed.updated_alleles_pre[1]), 
                                                                     mt_fixed.updated_alleles_1, mt_fixed.updated_alleles_2)).drop('updated_alleles_1','updated_alleles_2')
        if mt_fixed.aggregate_rows(~hl.agg.all((hl.len(mt_fixed.updated_alleles_pre[1])-hl.len(mt_fixed.updated_alleles[1])) >= 1)):
            raise ValueError('ERROR: the size by which the ALT allele shrank relative to pre must be >= 1 base.')
        if mt_fixed.aggregate_rows(~hl.agg.all((hl.len(mt_fixed.updated_alleles_pre[mt_fixed.idx_mult])-hl.len(mt_fixed.updated_alleles[mt_fixed.idx_mult])) == (hl.len(mt_fixed.updated_alleles_pre[mt_fixed.idx_single])-1))):
            raise ValueError('ERROR: the size by which the REF allele shrank relative to pre should be identical to the size of the shrinking of ALT allele.')
        if mt_fixed.aggregate_rows(~hl.agg.all(hl.len(mt_fixed.updated_alleles[mt_fixed.idx_single]) == 1)):
            raise ValueError('ERROR: the fixed REF allele should have length 1.')
        if mt_fixed.aggregate_rows(~hl.agg.all(mt_fixed.updated_alleles[mt_fixed.idx_single] == mt_fixed.updated_alleles[mt_fixed.idx_mult][0])):
            raise ValueError('ERROR: the fixed REF allele should be the same as the first base of the fixed ALT allele.')
        mt_fixed = mt_fixed.drop('updated_alleles_pre','idx_single','idx_mult')

    mt_fixed_adj = complex_swap_field(mt_fixed, ROW_BOOL_fields=ROW_BOOL_fields, ROW_R_fields=ROW_R_fields,
                                      ROW_custom_fields=ROW_custom_fields, ENTRY_R_fields=ENTRY_R_fields, flipped_idx=0,
                                      skip_locus_change=False, skip_allele_change=False)
    mt_fixed_adj = add_filter(mt_fixed_adj, 'AddGRCh38RefDeleToRefSiteIns')
    mt_fixed_adj = mt_fixed_adj.drop('mapped_alleles', 'tf_ref_dele', 'grch38_seq').persist()
    n_success_flipped = mt_fixed_adj.count_rows()
    
    mt_failed = mt_ins.filter_rows(~mt_ins.tf_here).drop('tf_here', 'tf_ref_dele')

    # make sure all previous keys are still there
    mt_fixed_for_test = mt_fixed_adj.key_rows_by()
    mt_fixed_for_test = mt_fixed_for_test.annotate_rows(locus_grch38=mt_fixed_for_test.locus, grch38_seq='A')
    mt_fixed_for_test = mt_fixed_for_test.annotate_rows(alleles = hl.agg.collect(mt_fixed_for_test.OriginalSelfRefAlleles)[0])
    mt_fixed_for_test = mt_fixed_for_test.annotate_rows(locus = mt_fixed_for_test.locus_self)
    mt_fixed_for_test = mt_fixed_for_test.select_rows(*[x for x in mt_failed.row]).key_rows_by('locus')
    mt_for_test = mt_fixed_for_test.union_rows(mt_failed).key_rows_by('locus','alleles').persist()
    if mt_for_test.count_rows() !=  mt.count_rows():
        raise ValueError('ERROR: when flipping fields for insertion rescue some records were lost/gained.')
    if (mt_for_test.anti_join_rows(mt.rows()).count_rows() != 0) | (mt.anti_join_rows(mt_for_test.rows()).count_rows() != 0):
        raise ValueError('ERROR: when flipping fields for insertion rescue some records were lost/gained.')
    
    return mt_fixed_adj.drop('locus_self').persist(), mt_failed.key_rows_by('locus','alleles').persist(), n_success_flipped


def confirm_ref(ref_mt, ref):
    ref_mt = ref_mt.annotate_rows(len = hl.len(ref_mt.alleles[0]))
    ref_mt = ref_mt.annotate_rows(seq = hl.get_sequence('chrM', ref_mt.locus.position, reference_genome=ref, after=ref_mt.len-1))
    ref_mt_f = ref_mt.filter_rows(ref_mt.seq != ref_mt.alleles[0])
    if ref_mt_f.count_rows() > 0:
        raise ValueError('ERROR: Final MT contains reference alleles that do not match the true reference sequence.')


def initialize_log(logging, individual_name, debug):
    log = open(logging, 'w')
    print('-------------- fix_liftover.py --------------', file=log)
    print('Pipeline for repair of complex Liftover edge cases on the mtDNA.', file=log)
    print('Rahul Gupta, 2022', file=log)
    print('', file=log)
    print('Sample name: ' + str(individual_name), file=log)
    print('Target reference: ' + str(MTREF), file=log)
    print('Start time: ' + datetime.now().strftime('%Y-%m-%d %H:%M:%S'), file=log)
    print('NOTE: This pipeline modifies the success VCF; be sure to use only success, rescued, and rejected VCFs from this in downstream applications.', file=log)
    if debug:
        print('WARNING: DEBUG MODE ENABLED. SEVERAL CONSISTANCY CHECKS ARE DISABLED.', file=log)
    print('Starting Liftover repair pipeline...', file=log)
    print('', file=log)
    return log


def write_outcome_log(s, final_mt, original_mt, failed_to_liftover, success_vcf, n_injected, complex_snp_liftover, new_haplotypes, n_for_in_body, 
                      n_indels_in_dele, n_complex_boundary, n_dele_spans_ins, n_dupe_removed, skip_norm, stage_ct_dict, output_txt_for_wdl, output_prefix, log):
    ct_original_success = success_vcf.count_rows()
    ct_total = original_mt.count_rows()
    ct_fixed = final_mt.count_rows()
    ct_failed_liftover = failed_to_liftover.count_rows()
    if ct_total != ct_fixed + ct_failed_liftover - n_injected:
        raise ValueError('ERROR: the final fixed table + final failed table do not have the same number of rows as the input table.')
    str_summary_1 = str(ct_original_success) + ' sites passed first round Liftover.'
    str_summary_2 = str(n_injected) + ' success sites were injected into failed MT.'
    str_summary_3 = 'Of ' + str(ct_total) + ' sites failing automated Liftover, we have repaired ' + str(ct_fixed-n_injected) + ' as well as ' + str(n_injected) + ' sites misclassified as a success.'
    filters_to_test = ['MismatchedRefAllele', 'NoTarget', 'IndelStraddlesMultipleIntevals', 'CannotLiftOver']
    if ct_failed_liftover > 0:
        failed_strings = (hl.str(failed_to_liftover.locus.position) + ':' + failed_to_liftover.alleles[0] + '>' + failed_to_liftover.alleles[1]).collect()
        set_of_all = hl.literal(set(filters_to_test))
        ct_no_filter = failed_to_liftover.filter_rows(failed_to_liftover.filters.intersection(set_of_all).length() == 0).count_rows()
        ct_each_filter = {x: failed_to_liftover.filter_rows(failed_to_liftover.filters.contains(x)).count_rows() for x in filters_to_test}
        ct_snv = failed_to_liftover.filter_rows(is_snv(failed_to_liftover)).count_rows()
        ct_indel = ct_failed_liftover - ct_snv
    else:
        ct_each_filter = {x: 0 for x in filters_to_test}
    ct_fixed_each_filter = {x: final_mt.filter_rows(final_mt.filters.contains(x)).count_rows() for x in filters_to_test}

    if output_txt_for_wdl:
        # here we output text files for the elements of stage_ct_dct and for the other constants
        tuple_list_output = [(re.sub('\\s','_',k.lower()), v) for k,v in stage_ct_dict.items()]
        other_list_output = [('round1_failed_sites', ct_total), ('round2_failed_sites', ct_failed_liftover), ('round2_fixed_sites', ct_fixed), 
                             ('round2_success_injected', n_injected),
                             ('snvs_complex_liftover', complex_snp_liftover), ('ref_insertion_new_haplos',new_haplotypes),
                             ('deletions_with_first_allele_and_internal_switches', n_for_in_body), ('hom_dele_with_indels_reverted', n_indels_in_dele),
                             ('het_dele_span_insert_repaired_with_complex_rework', n_complex_boundary),
                             ('het_deletions_span_insertions', n_dele_spans_ins), ('new_dupes_left_shift_failed', n_dupe_removed),
                             ('fixed_indelStraddlesMultipleIntevals', ct_fixed_each_filter['IndelStraddlesMultipleIntevals'])]
        full_list_output = tuple_list_output + other_list_output
        for name, value in full_list_output:
            with open(output_prefix + '.' + name + '.txt', 'w') as f:
                f.write('%d' % value)
        # now writing a table also
        dct_ints = {'s': s}
        dct_ints.update({'n_liftover_r2_'+left: right for (left, right) in tuple_list_output + other_list_output[1:len(other_list_output)]})
        dct_ints.update({'n_liftover_r1_'+left: right for (left, right) in other_list_output[0:1]})
        df = pd.DataFrame(dct_ints, index=[0])
        df.to_csv(output_prefix + '.all_int_outputs.txt', sep='\t', index=False)
    
    # Output log
    print('', file=log)
    print('OUTCOME LOG:', file=log)
    print(str_summary_1, file=log)
    print(str_summary_2, file=log)
    print(str_summary_3, file=log)
    for k,v in ct_fixed_each_filter.items():
        print('Repaired ' + str(v) + ' sites with ' + k + ' flag.', file=log)
    for k,v in stage_ct_dict.items():
        print('At the correction step for ' + k + ', ' + str(v) + ' sites were repaired.', file=log)
    if skip_norm:
        print('Left shifting was performed for ' + str(stage_ct_dict['left alignment of indels']) + ' sites via bcftools norm.', file=log)
    print('Complex Liftover was performed for ' + str(complex_snp_liftover) + ' SNVs. This involved fixing scenarios where ' + \
          'new SNVs were detected modifying GRCh38 ALT alleles. For example, if self-ref has G and GRCh38 has A>G, and self-ref recalling ' +\
          'produced G>C, we recode this with A>C.', file=log)
    print('New insertion haplotypes were constructed for ' + str(new_haplotypes) + ' sites. This involved fixing scenarios where ' + \
          'SNVs or insertions were detected in self-reference impacting a self-ref region with an insertion relative to GRCh38.', file=log)
    print(str(n_for_in_body) + ' heteroplasmic deletions needed fixing by changing the first allele to GRCh38 reference AND needed a change to an allele ' + \
          'in the body due to mismatch with GRCh38.', file=log)
    print(str(n_indels_in_dele) + ' heteroplasmic deletions had correct first base but had internal regions with sites that were indels in self-reference relative to GRCh38.', file=log)
    print(str(n_complex_boundary) + ' heteroplasmic deletions spanned one side of a homoplasmic insertion (rel. to GRCh38) AND had to subsequently have internal base changes repaired also.', file=log)
    if ct_failed_liftover > 0:
        print('INVESTIGATING ' + str(ct_failed_liftover) + ' LIFTOVER FAILURES:', file=log)
        for k,v in ct_each_filter.items():
            print(str(v) + ' had ' + k + ' flag.', file=log)
        print(str(ct_no_filter) + ' had none of the above flags.', file=log)
        print(str(ct_snv) + ' were SNVs.', file=log)
        print(str(ct_indel) + ' were indels.', file=log)
        print(str(n_dele_spans_ins) + ' were deletions that spanned insertion boundaries.', file=log)
        print(str(n_dupe_removed) + ' were lifted over but ended up duplicating another variant lifted by Picard and thus were failed. Currently this only happens if the new duplicate was generated by left-shifting.', file=log)
        print('In self-reference coordinates, the following failed to liftover:', file=log)
        for x in failed_strings:
            print(x, file=log)


def main(vcf_file, success_vcf_file, individual_name, self_to_ref_chain, ref_to_self_chain, self_fasta, self_fai, self_homoplasmies, logging,
         reference_fasta, reference_fai, allow_NONREF, simple_entry_field_correction, output_prefix, export_homoplasmic_deletions_coverage, output_txt_for_wdl, 
         skip_checkpoint, skip_norm, always_fail_on_dupe, debug):
    ###### SET UP SELF-REFERENCE AND LOGGING ########
    individual_name_input = individual_name
    individual_name = compatiblify_sample_name(individual_name)
    self_ref = hl.ReferenceGenome(individual_name, ['chrM'], {'chrM':fai_to_len(self_fai)}, mt_contigs=['chrM'])
    ref = hl.ReferenceGenome(MTREF, ['chrM'], {'chrM':fai_to_len(reference_fai)}, mt_contigs=['chrM'])
    self_ref.add_sequence(self_fasta, self_fai)
    ref.add_sequence(reference_fasta, reference_fai)
    self_ref.add_liftover(self_to_ref_chain, MTREF)
    ref.add_liftover(ref_to_self_chain, self_ref)
    chain_table = read_chain_file(self_to_ref_chain, MTREF, self_ref)
    log = initialize_log(logging, individual_name_input, debug)


    ###### READ REJECTED VCF ########
    # Key of rejection reasons:
    # MismatchedRefAllele: GRCh38 reference is not the same as original reference
    # NoTarget: Target location does not exist (probably inside an indel)
    # IndelStraddlesMultipleIntevals: Target locus is different in size from the source
    # CannotLiftOver: ??
    failed_to_liftover, failed_to_liftover_meta = read_mito_vcf(vcf_file, self_ref)
    failed_to_liftover = add_filter(failed_to_liftover, 'FailedPicardLiftoverVcf')


    ###### READ SUCCESS VCF ########
    success_vcf_original, success_vcf_meta = read_mito_vcf(success_vcf_file, ref, debug=debug)


    ###### ADD FILTERS FROM HOMOPLASMIES VCF ########
    homoplasmies, homoplasmies_meta = read_mito_vcf(self_homoplasmies, self_ref, avoid_dropping_gt=True)
    if homoplasmies.anti_join_rows(failed_to_liftover.rows()).count_rows() > 0:
        raise ValueError('ERROR: All homoplasmies should be found in the VCF of sites that failed to Liftover.')
    failed_to_liftover_homoplas = failed_to_liftover.semi_join_rows(homoplasmies.rows())
    failed_to_liftover_homoplas = failed_to_liftover_homoplas.annotate_rows(filters = hl.empty_set('tstr').union(homoplasmies.rows()[failed_to_liftover_homoplas.row_key].filters))
    failed_to_liftover_homoplas = add_filter(failed_to_liftover_homoplas, 'ForceCalledHomoplasmy')
    failed_to_liftover = failed_to_liftover.anti_join_rows(homoplasmies.rows()).union_rows(failed_to_liftover_homoplas).persist()
    filters_in_homoplasmies_table = list({y for x in homoplasmies.filters.collect() for y in ([x] if x is None else list(x))})
    filters_in_homoplasmies_table = [x for x in filters_in_homoplasmies_table if x is not None]
    filters_missing = [x for x in filters_in_homoplasmies_table if x not in failed_to_liftover_meta['filter'].keys()]
    for filt in filters_missing:
        failed_to_liftover_meta['filter'].update({filt: homoplasmies_meta['filter'][filt]})


    ###### DROP INFO FIELDS THAT ARE MISSING IN ALL RECORDS ########
    to_keep_1 = [x for x in list(failed_to_liftover.info) if check_missing_row_field(failed_to_liftover, failed_to_liftover.info[x])]
    to_keep_2 = [x for x in list(success_vcf_original.info) if check_missing_row_field(success_vcf_original, success_vcf_original.info[x])]
    to_keep = list(set(to_keep_1 + to_keep_2))
    failed_to_liftover, success_vcf = unify_info(failed_to_liftover, success_vcf_original, to_keep)
    shared_info = unify_dicts(failed_to_liftover_meta['info'], success_vcf_meta['info'], to_keep)
    failed_to_liftover_meta.update({'info': deepcopy(shared_info)})
    success_vcf_meta.update({'info': deepcopy(shared_info)})


    ###### INJECT "LIFTOVER SUCCESS" VARIANTS THAT ACTUALLY NEED FIXING #########
    print('Scanning success VCF for any variants that acutally need repair...', file=log)
    if debug: # when debugging, success VCF is usually sparse
        entries_to_add = [x for x in failed_to_liftover.entry if x not in success_vcf.entry]
        success_vcf = success_vcf.annotate_entries(**{x: hl.missing(failed_to_liftover[x].dtype) for x in entries_to_add})
        success_vcf = success_vcf.select_entries(*[x for x in failed_to_liftover.entry])
    failed_to_liftover, success_vcf, n_injected = inject_success_variants_to_fix(failed_to_liftover, success_vcf, homoplasmies.rows(), self_ref)


    ###### FINALIZE INPUTS AND RUN GLOBAL CONSISTANCY CHECKS ON VCFs ######
    failed_to_liftover = failed_to_liftover.annotate_rows(locus_grch38 = hl.liftover(failed_to_liftover.locus, MTREF))
    failed_to_liftover = failed_to_liftover.annotate_rows(grch38_seq = failed_to_liftover.locus_grch38.sequence_context())
    #if not skip_checkpoint:
        #failed_to_liftover = failed_to_liftover.checkpoint('failed_tmp.mt')
        #success_vcf = success_vcf.checkpoint('success_tmp.mt')
    original_force_homoplasmies = failed_to_liftover.semi_join_rows(homoplasmies.rows())
    global_consistancy_checks(failed_to_liftover, allow_NONREF, self_ref, debug)
    global_consistancy_checks(success_vcf, allow_NONREF, ref, debug)
    print('Files imported and initial checks passed.', file=log)


    ###### FIX SNVs SHOWING MismatchRefAllele ########
    print('Looking for and repairing SNVs with MismatchRefAllele...', file=log)
    tf_mismatch_easy = ~failed_to_liftover.filters.contains('NoTarget') \
        & ~failed_to_liftover.filters.contains('CannotLiftOver') \
        & ~failed_to_liftover.filters.contains('IndelStraddlesMultipleIntevals') \
        & is_snv(failed_to_liftover)
    need_to_flip = failed_to_liftover.filter_rows(tf_mismatch_easy)
    flipped, ct_complex = swap_alleles(need_to_flip, failed_to_liftover_meta, allow_NONREF=allow_NONREF, log=log)
    failed_to_liftover = failed_to_liftover.filter_rows(~tf_mismatch_easy).persist()


    ###### FIX INDELS CREATED DURING SELF-REFERENCE PRODUCTION ########
    insertions = get_insertion_sites(chain_table, MTREF, self_ref)
    deletions = get_deletion_sites(chain_table, MTREF, self_ref)
    confirm_reference_indels(failed_to_liftover, insertions.union(deletions))
    all_possible_sites_insertions = explode_indel(insertions, insertions.coord_s, insertions.ns, self_ref, 'sites')

    # Remove any deletions spanning insertions, as these are explicitly unsupported.
    dele_span_insertion_boundary = dele_spans_insertion(failed_to_liftover, all_possible_sites_insertions, self_ref, ignore_identical=True)
    failed_to_liftover = failed_to_liftover.anti_join_rows(dele_span_insertion_boundary.rows())

    # Start with insertions.
    # All GRCh38 reference sites should now have associated variant calls.
    # The simplest case is simply reversing the insertion relative to GRCh38 -- this is handled here.
    # More complex is situations where there are SNVs or indels within an insertion relative to GRCh38. These are also handled here.
    # Note this ONLY deals with variation related to homoplasmic insertions relative to GRCh38.
    print('Looking for and repairing sites overlapping insertions in self-reference relative to GRCh38...', file=log)
    tf_reference_indel = ~failed_to_liftover.filters.contains('CannotLiftOver') & ~failed_to_liftover.filters.contains('IndelStraddlesMultipleIntevals')
    insertions_to_fix = failed_to_liftover.filter_rows(tf_reference_indel).key_rows_by('locus').semi_join_rows(all_possible_sites_insertions).key_rows_by('locus', 'alleles')
    insertions_to_fix = insertions_to_fix.annotate_rows(is_reversed_grch38_to_self = hl.is_defined(homoplasmies.rows()[insertions_to_fix.row_key]))
    failed_to_liftover_tmp = failed_to_liftover.filter_rows(~tf_reference_indel)
    failed_to_liftover = failed_to_liftover_tmp.union_rows(failed_to_liftover.filter_rows(tf_reference_indel).key_rows_by('locus').anti_join_rows(all_possible_sites_insertions).key_rows_by('locus', 'alleles')).persist()
    fixed_ref_insertions, n_new_haplotypes = fix_ref_insertions(insertions_to_fix, failed_to_liftover_meta, all_possible_sites_insertions, self_reference=self_ref, allow_NONREF=allow_NONREF, log=log)

    # Deal with deletions from the reference
    # Should always be present in the callset because we force-call changes from reference.
    # Filter to only those deletions specified in the reference.
    print('Looking for and repairing deletions in self-reference relative to GRCh38...', file=log)
    tf_reference_indel = ~failed_to_liftover.filters.contains('CannotLiftOver') & ~failed_to_liftover.filters.contains('IndelStraddlesMultipleIntevals')
    deletions = deletions.annotate(self_ref_alleles = [deletions.alt_allele, deletions.ref_allele]).key_by('coord_s', 'self_ref_alleles')
    deletions_to_fix = failed_to_liftover.filter_rows(tf_reference_indel).semi_join_rows(deletions)
    deletions_to_fix = deletions_to_fix.annotate_rows(from_deletions_locus38 = deletions[deletions_to_fix.row_key].coord_t)
    if deletions_to_fix.aggregate_rows(~hl.agg.all(deletions_to_fix.from_deletions_locus38 == deletions_to_fix.locus_grch38)):
        raise ValueError('ERROR: the Liftover-based GRCh38 coordinates for deletions from GRCh38 do not match those obtained form the chain file.')
    else:
        deletions_to_fix = deletions_to_fix.drop('from_deletions_locus38')
    deletions_to_fix = deletions_to_fix.annotate_rows(grch38_seq = deletions[deletions_to_fix.row_key].ref_allele)
    failed_to_liftover = failed_to_liftover_tmp.union_rows(failed_to_liftover.filter_rows(tf_reference_indel).anti_join_rows(deletions)).persist()
    flipped_deletions, _ = swap_alleles(deletions_to_fix, failed_to_liftover_meta, allow_NONREF=allow_NONREF, log=log, fail_on_complex=True)


    ###### DEAL WITH STRAGGLING INDELS ########
    # ~failed_to_liftover.filters.contains('IndelStraddlesMultipleIntevals') removed 220615 to allow additional edge case processing
    # Deal with the case of indels where the first base is changed because of a reference SNV
    print('Looking for and repairing indels which failed Liftover because their first base differs from GRCh38...', file=log)
    homoplasmies_snvs = homoplasmies.filter_rows(is_snv(homoplasmies))
    failed_to_liftover = failed_to_liftover.key_rows_by('locus')
    tf_swap = ~failed_to_liftover.filters.contains('CannotLiftOver') \
        & hl.any(hl.map(lambda x: ~hl.literal(NONINDEL).contains(x), failed_to_liftover.alleles))
    for_swapping_first_site = failed_to_liftover.filter_rows(tf_swap).semi_join_rows(homoplasmies_snvs.rows().key_by('locus')).key_rows_by('locus','alleles')
    matching_homoplasmies_firstswap = homoplasmies_snvs.key_rows_by('locus').semi_join_rows(for_swapping_first_site.rows().key_by('locus'))
    failed_to_liftover_tmp = failed_to_liftover.filter_rows(~tf_swap)
    failed_to_liftover = failed_to_liftover_tmp.union_rows(failed_to_liftover.filter_rows(tf_swap).anti_join_rows(homoplasmies_snvs.rows().key_by('locus'))).key_rows_by('locus','alleles').persist()
    first_site_swapped, needs_body_changes, n_for_in_body = make_first_indel_site_ref(for_swapping_first_site, failed_to_liftover_meta, matching_homoplasmies_firstswap, flipped, self_ref, ref, all_possible_sites_insertions, allow_NONREF)

    # Deal with the case of deletions where the body of the allele spans changed alleles
    print('Looking for and repairing deletions which failed Liftover because of internal differences with GRCh38...', file=log)
    tf_fix_deletion = ~failed_to_liftover.filters.contains('CannotLiftOver') \
        & ~hl.literal(NONINDEL).contains(failed_to_liftover.alleles[0]) \
        & hl.literal(NONINDEL).contains(failed_to_liftover.alleles[1])
    candidates_for_deletion_fix = failed_to_liftover.filter_rows(tf_fix_deletion)
    failed_to_liftover_tmp = failed_to_liftover.filter_rows(~tf_fix_deletion)
    fixed_deletion, did_not_fix_deletion, n_indels_in_dele = recode_deletion_allele(candidates_for_deletion_fix, needs_body_changes, all_possible_sites_insertions, ref, self_ref, allow_NONREF, self_homoplasmies=original_force_homoplasmies, mt_meta=failed_to_liftover_meta, allow_spanning=False)
    failed_to_liftover_tmp = add_missing_new_entry_fields(failed_to_liftover_tmp)
    failed_to_liftover = failed_to_liftover_tmp.union_rows(did_not_fix_deletion)

    # Deal with the case of a heteroplasmic insertion at the same starting position as a reference deletion (e.g., homoplasmic self-ref insertion)
    print('Looking for and repairing insertions which share a first base with a reference deletion...', file=log)
    tf_fix_ins_overlap = failed_to_liftover.filters.contains('InsertionSharesForceCalledInsertion') \
        & is_insertion(failed_to_liftover)
    candidates_for_ins_overlap_fix = failed_to_liftover.filter_rows(tf_fix_ins_overlap).persist()
    failed_to_liftover_tmp = failed_to_liftover.filter_rows(~tf_fix_ins_overlap)
    fixed_ins_overlap, did_not_fix_ins, n_fixed_ins_overlap = fix_insertion_overlapping_ref_dele(candidates_for_ins_overlap_fix, failed_to_liftover_meta, flipped_deletions, fixed_ref_insertions, insertions, deletions, rely_on_left_shift=True)
    failed_to_liftover = failed_to_liftover_tmp.union_rows(did_not_fix_ins)


    ###### TRY TO RESOLVE DELETION BOUNDARY CASES ########
    print('Trying to rescue sites that span one boundary of a reference insertion...', file=log)
    tf_fix_boundary = ~dele_span_insertion_boundary.filters.contains('CannotLiftOver') \
        & ~dele_span_insertion_boundary.filters.contains('IndelStraddlesMultipleIntevals') \
        & ~hl.literal(NONINDEL).contains(dele_span_insertion_boundary.alleles[0]) \
        & hl.literal(NONINDEL).contains(dele_span_insertion_boundary.alleles[1])
    candidates_for_boundary_fix = dele_span_insertion_boundary.filter_rows(tf_fix_boundary)
    dele_span_insertion_boundary_tmp = dele_span_insertion_boundary.filter_rows(~tf_fix_boundary)
    fixed_boundary, did_not_fix_boundary, n_indels_dele_2, n_complex_boundary, n_rightside_sharedstart, n_leftside_span = resolve_deletion_boundary_cases(candidates_for_boundary_fix, failed_to_liftover_meta, all_possible_sites_insertions, original_force_homoplasmies, insertions, ref, self_ref, allow_NONREF)
    dele_span_insertion_boundary_tmp = add_missing_new_entry_fields(dele_span_insertion_boundary_tmp)
    dele_span_insertion_boundary = dele_span_insertion_boundary_tmp.union_rows(did_not_fix_boundary)
    print('Allele rescue completed.', file=log)


    ###### COLLATE RESULTS ########
    failed_to_liftover = failed_to_liftover.union_rows(dele_span_insertion_boundary).drop('locus_grch38', 'grch38_seq').persist()
    final_metadata = global_modify_meta(failed_to_liftover_meta)
    final_success_metadata = global_modify_meta(success_vcf_meta)
    final_mt = flipped.union_rows(fixed_ref_insertions).union_rows(flipped_deletions).union_rows(first_site_swapped
                     ).union_rows(fixed_deletion).union_rows(fixed_ins_overlap).union_rows(fixed_boundary).persist()
    if not simple_entry_field_correction:
        print('Performing fancy AF flip for force-called alleles based on other alleles sharing a start site...', file=log)
        final_mt, n_fancy_flip = fancy_flip_entries(final_mt, failed_to_liftover_meta, success_vcf, debug=debug)
        print('Performing change of reference INFO/entry fields for successfully lifted alleles sharing a position with force-called indels...', file=log)
        final_success_mt, n_success_flipped = flip_success_fields(success_vcf, failed_to_liftover_meta, final_mt, insertions, deletions)
        fancy_flip_dict = {'sites with fancy flip': n_fancy_flip, 'success sites flipped': n_success_flipped}
    else:
        final_success_mt = success_vcf
        fancy_flip_dict = {}
    _, _, _, ROW_drop_fields, _, _, ENTRY_drop_fields = get_fields_for_flipping(failed_to_liftover_meta)
    final_success_mt = drop_info(final_success_mt.drop(*ENTRY_drop_fields), ROW_drop_fields).persist()
    final_mt = drop_info(final_mt.drop(*ENTRY_drop_fields), ROW_drop_fields).persist()
    if not skip_norm:
        print('Running bcftools norm to left-align and shift indels...', file=log)
        final_mt, final_metadata_2, n_left_aligned, message = left_align_and_normalize(final_mt, final_metadata, ref, reference_fasta)
        print('bcftools returned: ' + message, file=log)
    else:
        final_mt = final_mt.key_rows_by().select_rows(*final_success_mt.row).key_rows_by('locus','alleles')
        final_metadata_2 = final_metadata
    
    # Check for duplicates
    if final_mt.rows().distinct().count() != final_mt.rows().count():
        raise ValueError('ERROR: during liftover duplicated alleles were introduced. This is not currently supported.')
    if not always_fail_on_dupe:
        final_mt, failed_to_liftover, failed_to_liftover_meta, n_failed_dupe = rescue_duplicated_alleles(final_mt, final_success_mt, failed_to_liftover, final_metadata_2, failed_to_liftover_meta, self_ref)
    else:
        n_failed_dupe = 0
    final_success_ht = final_success_mt.rows()
    final_ht = final_mt.rows()
    final_ht = final_ht.select(**{k: final_ht[k] for k in final_success_ht.row if k not in final_ht.key})
    for_test_all_data = final_success_ht.union(final_ht)
    if for_test_all_data.distinct().count() != for_test_all_data.count():
        raise ValueError('ERROR: during liftover duplicated alleles were introduced when testing with both fixed and success VCF. This is not currently supported.')
    
    # Check to make sure all REF alleles are indeed ref
    confirm_ref(final_mt, ref)

    # Confirm that resultant alleles are appropriately formatted
    if final_mt.aggregate_rows(~hl.agg.all(hl.any(hl.map(lambda x: hl.literal(NONINDEL).contains(x), final_mt.alleles)))):
        raise ValueError('ERROR: despite using left shifting, there is at least one allele with no NONINDEL variants (e.g., [GCA, GCACA]). This is not allowed.')

    ###### OUTPUT FILES AND LOG ########
    dict_of_filters = {'SNVs with MismatchRefAllele': flipped.count_rows(),
                       'homoplasmic insertions rel to GRCh38': fixed_ref_insertions.count_rows(),
                       'homoplasmic deletions rel to GRCh38': flipped_deletions.count_rows(),
                       'heteroplasmic indels first site swapped to GRCh38': first_site_swapped.count_rows(),
                       'heteroplasmic deletions with internal sites swapped to GRCh38': fixed_deletion.count_rows(),
                       'het insertions sharing LHS with hom ref deletion': n_fixed_ins_overlap,
                       'heteroplasmic deletions sharing LHS with homoplasmic insertion spanning RHS': n_rightside_sharedstart,
                       'heteroplasmic deletions spanning only LHS of homoplasmic insertion': n_leftside_span}
    dict_of_filters.update(fancy_flip_dict)
    dict_of_filters.update({'left alignment of indels': (float("nan") if skip_norm else n_left_aligned)})
    write_outcome_log(individual_name_input, final_mt, read_mito_vcf(vcf_file, self_ref)[0], failed_to_liftover, final_success_mt, n_injected, ct_complex, n_new_haplotypes, 
                      n_for_in_body, n_indels_in_dele+n_indels_dele_2, n_complex_boundary, dele_span_insertion_boundary.count_rows(), n_dupe_removed=n_failed_dupe,
                      skip_norm=skip_norm, stage_ct_dict=dict_of_filters, output_txt_for_wdl=output_txt_for_wdl, output_prefix=output_prefix, log=log)
    hl.export_vcf(final_mt, output_prefix + '.fixed.vcf.bgz', metadata=final_metadata_2)
    hl.export_vcf(final_success_mt, output_prefix + '.updated_success.vcf.bgz', metadata=final_success_metadata)
    hl.export_vcf(failed_to_liftover, output_prefix + '.rejected.vcf.bgz', metadata=failed_to_liftover_meta)
    if export_homoplasmic_deletions_coverage:
        ht_dele_cov = produce_ref_deletions_table(flipped_deletions, deletions)
        ht_dele_cov.export(output_prefix + '.deletions_coverage.tsv')
    print('', file=log)
    print('End time: ' + datetime.now().strftime('%Y-%m-%d %H:%M:%S'), file=log)
    log.close()



parser = argparse.ArgumentParser()
parser.add_argument('--logging', required=True, default='fix_liftover.log', help='Path to output log.')
parser.add_argument('--vcf-file', required=True, help='The VCF file of reads failing liftover.')
parser.add_argument('--success-vcf-file', required=True, help='VCF file of variants that passed first round liftover. Used for the fancy AF correction.')
parser.add_argument('--self-homoplasmies', required=True, help='VCF with self-coordinate homoplasmies. Filters from here will be grafted onto the final VCF.')
parser.add_argument('--individual-name', default='self_reference', help='ID of individual with self-reference.')
parser.add_argument('--self-to-ref-chain', required=True, help='Chain file mapping from self-ref to reference.')
parser.add_argument('--ref-to-self-chain', required=True, help='Chain file mapping from reference to self-ref.')
parser.add_argument('--self-fasta', required=True, help='Self-reference fasta file.')
parser.add_argument('--self-fai', required=True, help='Self-reference fai file.')
parser.add_argument('--reference-fasta', required=True, help='Reference fasta file.')
parser.add_argument('--reference-fai', required=True, help='Reference fai file.')
parser.add_argument('--allow-NONREF', action='store_true', help='If true, will try to handle alleles with <NON_REF> coding.')
parser.add_argument('--export-homoplasmic-deletions-coverage', action='store_true', help='If true, will ')
parser.add_argument('--output-prefix', required=True, help='Output prefix for outputting results. Will add .fixed.vcf.bgz and .rejected.vcf.bgz suffixes.')
parser.add_argument('--output-txt-for-wdl', action='store_true', help='If true, will use --output-prefix to also output txt files with numbers counting the number of sites fixed at each step.')
parser.add_argument('--simple-entry-field-correction', action='store_true', help='If enabled, will simply do 1-AF when flipping reference alleles and will not change any fields in Liftover passed alleles. ' + \
                                                                                 'Otherwise, a more sophisticated approach is used, reading in all variant calls and then doing total-AF where total are 1-the AF ' + \
                                                                                 'estimates for all sites sharing that locus. Further we now flip INFO/entry fields for variants that passed Liftover but share a site with a failed reference indel.')
parser.add_argument('--debug', action='store_true', help='If enabled, skips a few consistency checks to make using custom edge case VCFs easier to make.')
parser.add_argument('--skip-checkpoint', action='store_true')
parser.add_argument('--skip-norm', action='store_true', help='If enabled, will avoid using bcftools norm to left-align indels.')
parser.add_argument('--always-fail-on-dupe', action='store_true', help='If enabled, will always fail on any type of post-liftover duplicate variant.')


if __name__ == '__main__':
    args = parser.parse_args()
    main(**vars(args))
