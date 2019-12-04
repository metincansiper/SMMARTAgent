import sys
import os

import logging
from bioagents import Bioagent
from indra.statements import Agent
from .tumor_board_agent import TumorBoardAgent
from indra.sources.trips.processor import TripsProcessor
from indra.databases import hgnc_client
from kqml import KQMLModule, KQMLPerformative, KQMLList, KQMLString, KQMLToken


logging.basicConfig(format='%(levelname)s: %(name)s - %(message)s',
                    level=logging.INFO)
logger = logging.getLogger('TUMORBOARDA')

class TumorBoardModule(Bioagent):
    name = 'TUMORBOARDA'
    tasks = ['CREATE-TUMOR-BOARD-REPORT']

    def __init__(self, **kwargs):
        self.TBA = TumorBoardAgent()
        super(TumorBoardModule, self).__init__(**kwargs)

    def respond_create_tumor_board_report(self, content):
        patient_id = content.gets('patient-id')

        reply = KQMLList('SUCCESS')

        res = self.TBA.create_tumor_board_report(patient_id)

        if res == None:
            return self.make_failure('INVALID_PID')

        if len(res) == 0:
            return self.make_failure('NO_GENES_FOUND')

        genes = _get_genes_cljson(res)
        reply.set('genes', genes)

        return reply

def _get_kqml_names(kqmlList):
    """Given a kqml list returns the names of sublists in the list"""
    if not kqmlList:
        return None

    arr = kqmlList.data;
    if len(arr) == 0:
        return []

    if not isinstance(arr[0], KQMLList):
        arr = [kqmlList]

    res = list(map(lambda kl: kl.get('NAME').string_value(), arr))

    return res

def _get_genes_cljson(gene_names):
    agents = []
    for name in gene_names:
        db_refs={'TYPE':'ONT::GENE-PROTEIN'}
        hgnc_id = hgnc_client.get_hgnc_id(name)
        if hgnc_id:
            db_refs['HGNC'] = hgnc_id
        agents.append(Agent(name, db_refs=db_refs))

    res = Bioagent.make_cljson(agents)
    return res

if __name__ == "__main__":
    TumorBoardModule(argv=sys.argv[1:])
