import requests
import threading
from enum import Enum
import functools
from collections import Counter
from util.delimited_file_stream import DelimitedFileStream
from os import path

SCORE_THRESHOLD = 10
PATIENT_VARIANTS_PATH = 'input/patient_variants.txt'
CENSUS_PATH = 'input/census.tsv'
FDA_DRUG_LIST_PATH = 'input/drug_list_comprehensive.txt'
PC_URL = 'https://www.pathwaycommons.org/sifgraph/v1/neighborhood'
PC_PATTERNS = ['CONTROLS_STATE_CHANGE_OF', 'CONTROLS_EXPRESSION_OF', 'IN_COMPLEX_WITH']
PC_LIMIT = 1
PC_DIR = 'BOTHSTREAM'
CONTROLS_EXPRESSION_OF = 'controls-expression-of'.upper()
IN_COMPLEX_WITH = 'in-complex-with'.upper()
CONTROLS_STATE_CHANGE_OF = 'controls-state-change-of'.upper()
SPECIFIC_FDA_DISASE_TYPES = {'Breast cancer', 'Pancreatic cancer', 'Prostate cancer'}
SPECIFIC_CENSUS_DISASE_TYPES = {'breast cancer', 'breast'}

class SCORES:
    EXP_CONTROLLED_BY_VARIANT = 5.0
    CONTROLS_EXP_OF_VARIANT = 5.0
    IN_COMPLEX_WITH_VARIANT = 5.0
    CHANGES_STATE_OF_VARIANT = 5.0
    STATE_CHANGED_BY_VARIANT = 5.0
    SPECIFIC_CENSUS_GENE = 5.0
    OTHER_CENSUS_GENE = 2.5
    SPECIFIC_FDA_DRUG_TARGET = 5.0
    OTHER_FDA_DRUG_TARGET = 2.5

def get_fda_drug_targets():

    specific = {}
    other = {}

    def increment_field(related_map, field_name):
        related_map[field_name] = related_map.get(field_name, 0) + 1

    def on_data(tabs):
        disase_types = set(tabs[1].split('; '))
        genes = tabs[2].split('; ')

        specific_disase_type = bool( disase_types & SPECIFIC_FDA_DISASE_TYPES )
        related_map = specific if specific_disase_type else other

        for gene in genes:
            increment_field(related_map, gene)

    del_file_stream = DelimitedFileStream()
    del_file_stream.parse_file( file_path=FDA_DRUG_LIST_PATH, on_data=on_data )
    return [specific, other]

def get_census_gene_sets():
    specific = set()
    other = set()

    def get_tumor_types_set(s):
        return set(s.strip('"').split(', '))

    def on_data(tabs):
        tumor_types = get_tumor_types_set(tabs[9]) | get_tumor_types_set(tabs[10])
        gene = tabs[0]

        specific_tumor_type = bool( tumor_types & SPECIFIC_CENSUS_DISASE_TYPES )
        related_set = specific if specific_tumor_type else other

        related_set.add(gene)

    # TODO: give a warning when the path does not exists
    if path.exists(CENSUS_PATH):
        del_file_stream = DelimitedFileStream()
        del_file_stream.parse_file( file_path=CENSUS_PATH, on_data=on_data )

    return [specific, other]

[SPECIFIC_FDA_DRUG_TARGETS, OTHER_FDA_DRUG_TARGETS] = get_fda_drug_targets()
[SPECIFIC_CENSUS_GENE_SET, OTHER_CENSUS_GENE_SET] = get_census_gene_sets()

class TumorBoardAgent:

    def __init__(self):
        self.tumor_board_report = None
        self.sorted_results = None
        self.pc_evidences = None

    def create_tumor_board_report(self, patient_id):
        if isinstance( patient_id, str ):
            variants = self.get_variants( patient_id )
        elif isinstance( patient_id, list ):
            variants = patient_id

        if variants == None:
            return 'INVALID_PID'

        report = {}
        pc_evidences = {
            'pubmed_ids': {},
            'pc_links': {}
        }

        threads = list(map(lambda v: threading.Thread(target=self.fill_variant_report, args=(v,report,pc_evidences)), variants))

        for thread in threads:
            thread.start()

        for thread in threads:
            thread.join()

        self.tumor_board_report = report
        self.pc_evidences = pc_evidences

        report_sum = functools.reduce(lambda a,b : Counter(a)+Counter(b),report.values())
        most_common = report_sum.most_common()
        res = []

        for e in most_common:
            score = e[1]
            gene = e[0]
            if score > SCORE_THRESHOLD:
                res.append(gene)

        self.sorted_results = res
        return res

    def get_evidences_for(self, gene1, gene2):
        # TODO: throw error if self.pc_evidences is not set
        all_pubmed_ids = self.pc_evidences['pubmed_ids']
        all_pc_links = self.pc_evidences['pc_links']

        def nested_get(ds, g1, g2):
            return ds.get( g1, {} ).get( g2, set() )

        pubmed_ids1 = nested_get( all_pubmed_ids, gene1, gene2 )
        pubmed_ids2 = nested_get( all_pubmed_ids, gene2, gene1 )

        pc_links1 = nested_get( all_pc_links, gene1, gene2 )
        pc_links2 = nested_get( all_pc_links, gene2, gene1 )

        pubmed_ids = pubmed_ids1 | pubmed_ids2
        pc_links = pc_links1 | pc_links2

        res = []

        for pid in pubmed_ids:
            res.append(pid)

        for pcl in pc_links:
            res.append(pcl)

        return res

    def get_top_k_res(self, k):
        # TODO: throw error if self.sorted_results is not set
        return self.sorted_results[:k]

    def why_important(self, gene, k=2):
        variant_scores = {}
        # TODO throw error if self.tumor_board_report is not set
        for variant in self.tumor_board_report:
            neighbours = self.tumor_board_report[variant]
            if gene in neighbours:
                score = neighbours.get(gene)
                variant_scores[variant] = score

        c = Counter(variant_scores)
        top_most_tuples = c.most_common(k)
        top_most_names = list(map(lambda t: t[0], top_most_tuples))

        return top_most_names

    def fill_variant_report(self, variant, report, pc_evidences):
        variant = variant.upper()
        score = 0

        neighbours = set()
        var_pubmed_ids = {}
        var_pc_links = {}

        exp_controlled_by_variant = set()
        controls_exp_of_variant = set()
        in_complex_with_variant = set()
        changes_state_of_variant = set()
        state_changed_by_variant = set()

        r = self.query_pc( variant )
        text = r.text
        lines = text.splitlines()

        def split_evidence(s):
            s = str(s)
            if not s:
                return []

            return s.split(';')

        for line in lines:
            line = line.strip('\n') #cut off the return from the end of the line
            parts = line.split('\t') #cut line into pieces at the tabs
            entity1 = str(parts[0]).upper()
            edge = str(parts[1]).upper()
            entity2 = str(parts[2]).upper()
            pubmed_ids = split_evidence(parts[4])
            pc_links = split_evidence(parts[6])

            var_eq_e1 = entity1 == variant
            var_eq_e2 = entity2 == variant

            if not var_eq_e1 and not var_eq_e2:
                raise Exception('One of entity names should match the variant name where variant name is {} and the entitiy names are {}, {}', variant, entity1, entity2)

            other_side = entity2 if var_eq_e1 else entity1
            neighbours.add(other_side)

            n_pubmed_ids = var_pubmed_ids.get( other_side, set() )
            n_pc_links = var_pc_links.get( other_side, set() )

            for pid in pubmed_ids:
                # if pid:
                n_pubmed_ids.add(pid)

            # if no pubmed id is found use the pc_links instead
            if not pubmed_ids:
                for pcl in pc_links:
                    # if pcl:
                    n_pc_links.add(pcl)

            if n_pubmed_ids:
                var_pubmed_ids[ other_side ] = n_pubmed_ids
            if n_pc_links:
                var_pc_links[ other_side ] = n_pc_links

            s = None

            if edge==CONTROLS_EXPRESSION_OF:
                s = exp_controlled_by_variant if var_eq_e1 else controls_exp_of_variant
            elif edge==IN_COMPLEX_WITH:
                s = in_complex_with_variant
            elif edge==CONTROLS_STATE_CHANGE_OF:
                s = state_changed_by_variant if var_eq_e1 else changes_state_of_variant
            else:
                raise Exception('Edge type should one of {}. The edge type was: {}'.format([CONTROLS_EXPRESSION_OF,IN_COMPLEX_WITH,CONTROLS_STATE_CHANGE_OF], edge))

            s.add(other_side)

        scores = {}

        for neighbour in neighbours:
            score = 0
            if neighbour in exp_controlled_by_variant:
                score += SCORES.EXP_CONTROLLED_BY_VARIANT
            if neighbour in controls_exp_of_variant:
                score += SCORES.CONTROLS_EXP_OF_VARIANT
            if neighbour in in_complex_with_variant:
                score += SCORES.IN_COMPLEX_WITH_VARIANT
            if neighbour in changes_state_of_variant:
                score += SCORES.CHANGES_STATE_OF_VARIANT
            if neighbour in state_changed_by_variant:
                score += SCORES.STATE_CHANGED_BY_VARIANT
            if neighbour in SPECIFIC_CENSUS_GENE_SET:
                score += SCORES.SPECIFIC_CENSUS_GENE
            if neighbour in OTHER_CENSUS_GENE_SET:
                score += SCORES.OTHER_CENSUS_GENE
            if neighbour in SPECIFIC_FDA_DRUG_TARGETS:
                score += SPECIFIC_FDA_DRUG_TARGETS[neighbour] * SCORES.SPECIFIC_FDA_DRUG_TARGET
            if neighbour in OTHER_FDA_DRUG_TARGETS:
                score += OTHER_FDA_DRUG_TARGETS[neighbour] * SCORES.OTHER_FDA_DRUG_TARGET

            if score > 0:
                scores[neighbour] = score

        report[ variant ] = scores
        pc_evidences[ 'pubmed_ids' ][ variant ] = var_pubmed_ids
        pc_evidences[ 'pc_links' ][ variant ] = var_pc_links

    def get_variants(self, patient_id):
        variants = None

        def on_data(tabs):
            if tabs[0] == patient_id:
                nonlocal variants
                variants = tabs[1].split(',')
                # Signal to end the streaming by returning True
                return True

        del_file_stream = DelimitedFileStream()
        del_file_stream.parse_file( file_path=PATIENT_VARIANTS_PATH, on_data=on_data )

        return variants

    def query_pc(self, variant):
        params = self.get_pc_query_params(variant)
        r = requests.get(PC_URL, params)
        return r

    def get_pc_query_params(self, variant):
        params = { 'limit': PC_LIMIT, 'direction': PC_DIR, 'pattern': PC_PATTERNS, 'source': variant }
        return params
