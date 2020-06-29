import requests
import threading
from enum import Enum
import functools
from collections import Counter
from .parser.delimited_file_stream import DelimitedFileStream
from .util.paxtools import biopax_text_to_sbgn, DEFAULT_MAX_CONVERSIONS
from os import path, system
import heapq
import warnings
import pandas

file_dir = path.dirname(path.abspath(__file__))
parent_dir = path.dirname(file_dir)

def join_path(p):
    return path.join(parent_dir, p)

PATIENT_VARIANTS_PATH = join_path('input/patient_variants.txt')
CENSUS_PATH = join_path('input/census.tsv')
FDA_DRUG_LIST_PATH = join_path('input/drug_list_comprehensive.txt')
CANCER_NETWORK_PATH = join_path('input/cancer_network')
CN_MUTEC_FILE_PATH = path.join(CANCER_NETWORK_PATH, 'mutec.csv')
CN_SIF_OUTPUT_PATH = path.join(CANCER_NETWORK_PATH, 'network.sif')
CANCER_NETWORK_JAR_PATH = join_path('jar/cancer-network.jar')
CLINICAL_TRIALS_URL = 'https://clinicaltrials.gov/ct2/results/download_fields'
CBIO_MUTATIONS_URL = 'https://www.cbioportal.org/api/molecular-profiles/msk_impact_2017_mutations/mutations'

SCORE_THRESHOLD = 10
PC_NEIGHBORHOOD_URL = 'https://www.pathwaycommons.org/sifgraph/v1/neighborhood'
PC_GRAPH_URL = 'https://www.pathwaycommons.org/pc2/graph'
PC_LIMIT = 1
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

    if path.exists(CENSUS_PATH):
        del_file_stream = DelimitedFileStream()
        del_file_stream.parse_file( file_path=CENSUS_PATH, on_data=on_data )
    else:
        warnings.warn('Census input file is missing criteria will be executed without that!')

    return [specific, other]

[SPECIFIC_FDA_DRUG_TARGETS, OTHER_FDA_DRUG_TARGETS] = get_fda_drug_targets()
[SPECIFIC_CENSUS_GENE_SET, OTHER_CENSUS_GENE_SET] = get_census_gene_sets()

class TumorBoardAgent:

    def __init__(self, patient_id=None):
        self.tumor_board_report = None
        self.sorted_neighbors = None
        self.pc_evidences = None
        self.patient_id = patient_id
        if patient_id is not None:
            self.create_tumor_board_report(self.patient_id)

    def create_tumor_board_report(self, patient_id):
        variant_pairs = None

        if isinstance( patient_id, str ):
            variant_pairs = TumorBoardAgent.read_variant_pairs( patient_id )
        elif isinstance( patient_id, list ):
            variant_pairs = patient_id

        if variant_pairs == None:
            return None

        report = {}
        pc_evidences = {
            'pubmed_ids': {},
            'pc_links': {}
        }

        variant_genes = list(map(lambda p: p[0], variant_pairs))
        threads = list(map(lambda v: threading.Thread(target=self.fill_variant_report, args=(v,report,pc_evidences)), variant_genes))

        for thread in threads:
            thread.start()

        for thread in threads:
            thread.join()

        self.tumor_board_report = report
        self.pc_evidences = pc_evidences

        report_sum = functools.reduce(lambda a,b : Counter(a)+Counter(b),report.values())
        most_common = report_sum.most_common()
        sorted_neighbors = []

        for e in most_common:
            score = e[1]
            gene = e[0]
            if score > SCORE_THRESHOLD:
                sorted_neighbors.append(gene)

        self.sorted_neighbors = sorted_neighbors
        self.variant_pairs = variant_pairs
        return sorted_neighbors

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

    def get_top_neighbors(self, k):
        # TODO: throw error if self.sorted_neighbors is not set
        return self.sorted_neighbors[:k]

    def get_important_neighbors_of_gene(self, gene, k=None):
        if gene not in self.tumor_board_report:
            return None

        gene_report = self.tumor_board_report[gene]
        gene_neighbors = gene_report.keys()
        important_neighbors = self.get_top_neighbors(None)

        l = list(set(important_neighbors) & set(gene_neighbors))
        l.sort(key=lambda n: gene_report[n], reverse=True)

        return l[:k]


    def get_variant_pairs(self):
        return self.variant_pairs

    def why_important(self, gene, k=2):
        variant_scores = {}
        # TODO throw error if self.tumor_board_report is not set
        for variant_gene in self.tumor_board_report:
            neighbors = self.tumor_board_report[variant_gene]
            if gene in neighbors:
                score = neighbors.get(gene)
                variant_scores[variant_gene] = score

        c = Counter(variant_scores)
        top_most_tuples = c.most_common(k)
        top_most_names = list(map(lambda t: t[0], top_most_tuples))

        return top_most_names

    def fill_variant_report(self, variant_gene, report, pc_evidences):
        variant_gene = variant_gene.upper()
        score = 0

        neighbors = set()
        var_pubmed_ids = {}
        var_pc_links = {}

        exp_controlled_by_variant = set()
        controls_exp_of_variant = set()
        in_complex_with_variant = set()
        changes_state_of_variant = set()
        state_changed_by_variant = set()

        r = TumorBoardAgent.query_pc_neighborhood( variant_gene )
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

            var_eq_e1 = entity1 == variant_gene
            var_eq_e2 = entity2 == variant_gene

            if not var_eq_e1 and not var_eq_e2:
                raise Exception('One of entity names should match the variant name where variant name is {} and the entitiy names are {}, {}', variant_gene, entity1, entity2)

            other_side = entity2 if var_eq_e1 else entity1
            neighbors.add(other_side)

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

        for neighbor in neighbors:
            score = 0
            if neighbor in exp_controlled_by_variant:
                score += SCORES.EXP_CONTROLLED_BY_VARIANT
            if neighbor in controls_exp_of_variant:
                score += SCORES.CONTROLS_EXP_OF_VARIANT
            if neighbor in in_complex_with_variant:
                score += SCORES.IN_COMPLEX_WITH_VARIANT
            if neighbor in changes_state_of_variant:
                score += SCORES.CHANGES_STATE_OF_VARIANT
            if neighbor in state_changed_by_variant:
                score += SCORES.STATE_CHANGED_BY_VARIANT
            if neighbor in SPECIFIC_CENSUS_GENE_SET:
                score += SCORES.SPECIFIC_CENSUS_GENE
            if neighbor in OTHER_CENSUS_GENE_SET:
                score += SCORES.OTHER_CENSUS_GENE
            if neighbor in SPECIFIC_FDA_DRUG_TARGETS:
                score += SPECIFIC_FDA_DRUG_TARGETS[neighbor] * SCORES.SPECIFIC_FDA_DRUG_TARGET
            if neighbor in OTHER_FDA_DRUG_TARGETS:
                score += OTHER_FDA_DRUG_TARGETS[neighbor] * SCORES.OTHER_FDA_DRUG_TARGET

            if score > 0:
                scores[neighbor] = score

        report[ variant_gene ] = scores
        pc_evidences[ 'pubmed_ids' ][ variant_gene ] = var_pubmed_ids
        pc_evidences[ 'pc_links' ][ variant_gene ] = var_pc_links

    def get_neighbors_sif(self, gene, max_lines=None):
        self.create_cn_mutect_file(gene)
        system('java -jar {} {}'.format(CANCER_NETWORK_JAR_PATH, CANCER_NETWORK_PATH))

        neighbors = self.get_important_neighbors_of_gene(gene)
        neighbors = list(map(lambda n: n.lower(), neighbors))

        def get_score(line):
            tabs = line.split('\t')

            score = 0

            if len(tabs) > 1:
                gene1 = tabs[0].lower()
                gene2 = tabs[1].lower()

                if gene1 in neighbors:
                    score += neighbors.index(gene1)

                if gene2 in neighbors:
                    score += neighbors.index(gene2)

            return score

        def must_include(line):
            return get_score(line) > 0


        text = TumorBoardAgent.read_text_file(CN_SIF_OUTPUT_PATH)
        lines = text.splitlines()

        if max_lines and len(lines) > max_lines:
            lines = heapq.nlargest(max_lines, lines, key=get_score)

        lines = filter(must_include, lines)

        return '\n'.join(lines)

    @staticmethod
    def get_clinical_trials(drugs):
        params = {'cond': 'breast_cancer', 'term': drugs[0],
                    'down_count': 10000, 'down_fmt':'tsv'}

        r = requests.get(CLINICAL_TRIALS_URL, params)

        res = requests.get(CLINICAL_TRIALS_URL, params)
        text = res.text
        lines = text.splitlines()
        lines = list(map(lambda l: l.lower(), lines))

        def valid_line(line):
            # print(line)
            for drug in drugs:
                search_str = 'drug: ' + drug.lower()
                if search_str not in line:
                    return False

            return True

        def extract_info(line):
            vals = line.split('\t')
            title = vals[1]
            status = vals[2]
            url = vals[7]

            return { 'title': title, 'status': status, 'url': url }

        lines = filter(valid_line, lines)
        res = list(map(extract_info, lines))

        return res

    def create_cn_mutect_file(self, gene):
        with open(CN_MUTEC_FILE_PATH, 'w') as mutect_file:
            # The first line is for the header
            mutect_file.write('\n')

            # Only fill the colomn for the gene, fill the rest with a space char
            arr = [' \t'] * 5 + [ gene ] + ['\t '] * 9
            line = ''.join(arr)
            # The second line is for the gene of interest
            mutect_file.write(line)
            mutect_file.write('\n')

    @staticmethod
    def query_cbio_mutations_params(page_number):
        sample_list_id = 'msk_impact_2017_all'
        projection = 'DETAILED'
        page_size = '50000'
        direction = 'ASC'

        params = { 'sampleListId': sample_list_id,
                    'projection': projection,
                    'pageSize': page_size,
                    'pageNumber': page_number,
                    'direction': direction }

        return params

    @staticmethod
    def query_cbio_mutations(page_number=0):
        # print(page_number)
        params = TumorBoardAgent.query_cbio_mutations_params(page_number)
        # print(params)
        r = requests.get(CBIO_MUTATIONS_URL, params)
        js = r.json()

        if len(js) > 0:
            print(len(js))
            return [*js, *TumorBoardAgent.query_cbio_mutations(page_number + 1)]

        return []

    @staticmethod
    def get_grouped_cbio_mutations():
        mutations = TumorBoardAgent.query_cbio_mutations()
        groups = pandas.DataFrame(mutations).groupby("patientId").groups
        grouped = {}
        for patient_id in list(groups.keys()):
            indices = groups.get(patient_id)
            gene_variant = list(map(lambda i: [mutations[i].get('gene').get('hugoGeneSymbol'), mutations[i].get('proteinChange')],indices))
            variants_by_gene = {}
            for p in gene_variant:
                gene = str(p[0])
                variant = p[1]
                variants = variants_by_gene.get(gene, [])
                variants.append(variant)
                variants_by_gene[gene] = variants

            genes = list(variants_by_gene.keys())
            pairs = list( map( lambda g: [ g, variants_by_gene[ g ] ], genes ) )
            grouped[patient_id] = pairs


        return grouped

    @staticmethod
    def read_variant_pairs(patient_id):
        mutations_by_pid = TumorBoardAgent.get_grouped_cbio_mutations()
        mutations = mutations_by_pid.get(patient_id)
        return mutations

    @staticmethod
    def query_pc_neighborhood(variant_gene):
        params = TumorBoardAgent.get_pc_query_neighborhood_params(variant_gene)
        r = requests.get(PC_NEIGHBORHOOD_URL, params)
        return r

    @staticmethod
    def get_pc_query_neighborhood_params(variant_gene):
        DIR = 'BOTHSTREAM'
        PATTERNS = ['CONTROLS_STATE_CHANGE_OF', 'CONTROLS_EXPRESSION_OF', 'IN_COMPLEX_WITH']
        params = { 'limit': PC_LIMIT, 'direction': DIR, 'pattern': PATTERNS, 'source': variant_gene }
        return params

    @staticmethod
    def query_pc_pathsbetween(sources):
        params = TumorBoardAgent.get_pc_query_pathsbetween_params(sources)
        r = requests.get(PC_GRAPH_URL, params)
        return r

    @staticmethod
    def query_pc_pathsfromto(sources, targets):
        params = TumorBoardAgent.get_pc_query_pathsfromto_params(sources, targets)
        r = requests.get(PC_GRAPH_URL, params)
        return r

    @staticmethod
    def get_pc_query_pathsbetween_params(sources):
        KIND = 'PATHSBETWEEN'
        params = { 'kind': KIND, 'source': sources }
        return params

    @staticmethod
    def get_pc_query_pathsfromto_params(sources, targets):
        KIND = 'PATHSFROMTO'
        params = { 'kind': KIND, 'source': sources, 'target': targets }
        return params


    @staticmethod
    def get_pathsbetween_genes(sources, max_conversions=DEFAULT_MAX_CONVERSIONS):
        r = TumorBoardAgent.query_pc_pathsbetween(sources)
        text = r.text

        sbgn = biopax_text_to_sbgn(text)
        return sbgn

    @staticmethod
    def get_pathsfromto_genes(sources, targets, max_conversions=DEFAULT_MAX_CONVERSIONS):
        r = TumorBoardAgent.query_pc_pathsfromto(sources, targets)
        text = r.text

        sbgn = biopax_text_to_sbgn(text)
        return sbgn

    @staticmethod
    def read_text_file(f):
        with open(f, 'r') as content_file:
            content = content_file.read()
            return content
