import threading
import functools
from functools import lru_cache
import itertools
from collections import Counter
from .parser.delimited_file_stream import DelimitedFileStream
from .util.paxtools import biopax_text_to_sbgn, DEFAULT_MAX_CONVERSIONS
from os import path, system
import heapq
import warnings
import pandas
from bioagents.cbio_client import *
from indra.databases.hgnc_client import \
    get_hgnc_from_entrez, get_hgnc_name
# TODO: Remove dependency to clare
from clare.capabilities.util import get_agent_from_name
from bioagents.dtda.dtda import DTDA
from multiprocessing import Pool, cpu_count

file_dir = path.dirname(path.abspath(__file__))
parent_dir = path.dirname(file_dir)
_dtda = DTDA()


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

SCORE_THRESHOLD = 10
PC_SIFGRAPH_URL = 'https://www.pathwaycommons.org/sifgraph/v1/'
PC_GRAPH_URL = 'https://www.pathwaycommons.org/pc2/graph'
PC_SIF_PATHSBETWEEN_URL = PC_SIFGRAPH_URL + 'pathsbetween'
PC_NEIGHBORHOOD_URL = PC_SIFGRAPH_URL + 'neighborhood'
PC_LIMIT = 1
CONTROLS_EXPRESSION_OF = 'controls-expression-of'.upper()
IN_COMPLEX_WITH = 'in-complex-with'.upper()
CONTROLS_STATE_CHANGE_OF = 'controls-state-change-of'.upper()

class SCORES:
    EXP_CONTROLLED_BY_VARIANT = 10.0
    CONTROLS_EXP_OF_VARIANT = 10.0
    IN_COMPLEX_WITH_VARIANT = 10.0
    CHANGES_STATE_OF_VARIANT = 10.0
    STATE_CHANGED_BY_VARIANT = 10.0
    SPECIFIC_CENSUS_GENE = 10.0
    OTHER_CENSUS_GENE = 5.0
    SPECIFIC_FDA_DRUG_TARGET = 10.0
    OTHER_FDA_DRUG_TARGET = 5.0
    VARIANT = 20.0
    CNA = 20.0


class TumorBoardAgent:

    def __init__(self, patient_id=None):
        self.tumor_board_report = None
        self.sorted_neighbors = None
        self.pc_evidences = None
        self.patient_id = patient_id
        self.variant_pairs = []
        self.cna_pairs = []
        if patient_id is not None:
            self.patient = Patient(patient_id)
            self.create_tumor_board_report(self.patient_id)

    def create_tumor_board_report(self, patient_id):
        variant_pairs = []
        cna_pairs = []

        if isinstance( patient_id, str ):
            variant_pairs = self.read_variant_pairs()
            cna_pairs = self.read_cna_pairs()
        elif isinstance( patient_id, list ):
            variant_pairs = patient_id

        if not variant_pairs + cna_pairs:
            return None

        self.sample_info = self.patient.sample.clinical_info
        self.disease_name = self.sample_info.get('CANCER_TYPE', 'cancer').lower()
        self.specific_fda_drug_targets = TumorBoardAgent.query_target_genes(self.disease_name)
        self.other_fda_drug_targets = TumorBoardAgent.query_target_genes(self.disease_name, False);
        [self.specific_cencus_genes, self.other_cencus_genes] = TumorBoardAgent.get_census_gene_sets(self.disease_name)

        report = {}
        pc_evidences = {
            'pubmed_ids': {},
            'pc_links': {}
        }

        variant_genes = list(map(lambda p: p[0], variant_pairs))
        cna_genes = list(map(lambda p: p[0], cna_pairs))
        altered_genes = variant_genes + cna_genes

        threads = list(map(lambda v:
                           threading.Thread(target=self.fill_variant_report,
                                            args=(v, report, pc_evidences)),
                           altered_genes))

        for thread in threads:
            thread.start()

        for thread in threads:
            thread.join()

        self.tumor_board_report = report
        self.pc_evidences = pc_evidences

        report_sum = functools.reduce(lambda a, b: Counter(a) + Counter(b),
                                      report.values())

        # most_common = report_sum.most_common()

        def get_sorted_genes(report_sum):
            most_common = report_sum.most_common()
            sorted_genes = []

            for e in most_common:
                score = e[1]
                gene = e[0]
                if score > SCORE_THRESHOLD:
                    sorted_genes.append(gene)

            return sorted_genes

        sorted_neighbors = get_sorted_genes(report_sum)

        for g in variant_genes:
            report_sum[ g ] = report_sum[ g ] + SCORES.VARIANT

        for g in cna_genes:
            report_sum[ g ] = report_sum[ g ] + SCORES.CNA

        sorted_genes = get_sorted_genes(report_sum)

        self.sorted_neighbors = sorted_neighbors
        self.sorted_genes = sorted_genes
        self.variant_pairs = variant_pairs
        self.cna_pairs = cna_pairs

        return sorted_genes

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

    def get_top_genes(self, k):
        return self.sorted_genes[:k]

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
            if neighbor in self.specific_cencus_genes:
                score += SCORES.SPECIFIC_CENSUS_GENE
            if neighbor in self.other_cencus_genes:
                score += SCORES.OTHER_CENSUS_GENE
            if neighbor in self.specific_fda_drug_targets:
                score += self.specific_fda_drug_targets[neighbor] * SCORES.SPECIFIC_FDA_DRUG_TARGET
            if neighbor in self.other_fda_drug_targets:
                score += self.other_fda_drug_targets[neighbor] * SCORES.OTHER_FDA_DRUG_TARGET

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
                  'down_count': 10000, 'down_fmt': 'tsv'}

        res = requests.get(CLINICAL_TRIALS_URL, params)
        text = res.text
        lines = text.splitlines()
        lines = list(map(lambda l: l.lower(), lines))

        def valid_line(line):
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

            return {'title': title, 'status': status, 'url': url}

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

    def read_variant_pairs(self):
        mutations_by_gene = defaultdict(list)
        for mutation in self.patient.mutations:
            gene_name = get_hgnc_name(
                get_hgnc_from_entrez(str(mutation['entrezGeneId'])))
            change = mutation['proteinChange']
            mutations_by_gene[gene_name].append(change)
        return list(mutations_by_gene.items())

    def read_cna_pairs(self):
        alteration_map = {
            -2: 'DEL',
            -1: 'del',
            0: 'neu',
            1: 'amp',
            2: 'AMP',
        }
        cna_pairs = []
        for cna in self.patient.cnas:
            gene_name = get_hgnc_name(
                get_hgnc_from_entrez(str(cna['entrezGeneId'])))
            alteration_str = alteration_map[cna['alteration']]
            cna_pairs.append((gene_name, alteration_str))
        return cna_pairs

    def get_sif_pathsbetween(self):
        k = 30
        sources = self.get_top_genes(k)
        r = TumorBoardAgent.query_pc_sif_pathsbetween(sources)
        return r.text

    @staticmethod
    def flatten_2d_list(list2d):
        merged = list(itertools.chain(*list2d))
        return merged

    @staticmethod
    @lru_cache(maxsize=None)
    def query_target_genes(disease_name, for_self=True):
        drugs = TumorBoardAgent.query_fda_drugs(disease_name, for_self)
        agents = list(map(lambda n: get_agent_from_name(n), drugs))

        # TODO: should keep running in parallel?
        pool = Pool()
        genes = pool.map(TumorBoardAgent.find_drug_targets, agents)
        # genes = list(map(lambda a: list(_dtda.find_drug_targets(a)), agents))
        res = {}

        genes = TumorBoardAgent.flatten_2d_list(genes)

        for g in genes:
            n = res.get(g, 0) + 1
            res[g] = n

        return res

    @staticmethod
    def find_drug_targets(a):
        return list(_dtda.find_drug_targets(a))

    @staticmethod
    def query_fda_drugs(disease_name, for_self=True, skip=0, so_far=[]):
        limit = 500
        threshold = 20

        if len(so_far) > threshold:
            return so_far[0:threshold]

        res = TumorBoardAgent.query_fda_drugs_in_range(disease_name, for_self, skip, limit)

        if len(res) == 0:
            return so_far

        so_far = list(set(so_far + res))
        skip = skip + limit
        return TumorBoardAgent.query_fda_drugs(disease_name, for_self, skip, so_far)


    @staticmethod
    def query_fda_drugs_in_range(disease_name, for_self, skip, limit):
        FDA_EVENT_URL = 'https://api.fda.gov/drug/event.json'

        search_disease_name = disease_name

        if not for_self:
            search_disease_name = 'CANCER'

        search = 'patient.drug.drugindication:"' + search_disease_name + '"'

        params = { 'limit': limit, 'skip': skip, 'search': search }

        # if for_self:
        #     search = 'patient.drug.drugindication:"' + disease_name + '"'
        #     params['search'] = search

        r = requests.get(FDA_EVENT_URL, params)
        # r = requests.get('https://api.fda.gov/drug/event.json?search=patient.drug.drugindication:%22BREAST%20CANCER%22&limit=5')
        js = r.json()


        def get_drug_names(p):
            drugs = p.get('drug')

            def get_drug_name(d):
                indication = d.get('drugindication', '').upper()
                if for_self == (indication != disease_name):
                    return None
                gnames = d.get('openfda', {}).get('generic_name')

                if gnames != None and len(gnames) >= 0:
                    return gnames[0]
                return None

            names = list(map(get_drug_name , drugs))
            names = filter(lambda n: n != None, names)
            # return TumorBoardAgent.flatten_2d_list(names)
            return names

        results = js.get('results')
        patients = list(map(lambda r: r.get('patient'), results))
        drugs = TumorBoardAgent.flatten_2d_list(list(map(get_drug_names, patients)))
        drugs = list(set(drugs))

        return drugs


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
    def query_pc_sif_pathsbetween(sources):
        params = { 'source': sources }
        r = requests.get(PC_SIF_PATHSBETWEEN_URL, params)
        return r

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
    def get_census_gene_sets(disease_name):
        SPECIFIC_CENSUS_DISASE_TYPES = set([disease_name, disease_name.replace(' cancer', '')])
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

    @staticmethod
    def read_text_file(f):
        with open(f, 'r') as content_file:
            content = content_file.read()
            return content
