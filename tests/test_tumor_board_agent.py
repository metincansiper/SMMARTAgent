import json
from kqml import KQMLList, KQMLString, KQMLPerformative
from indra.statements import stmts_from_json

from tumor_board_agent.tumor_board_module import TumorBoardModule
from bioagents.tests.integration import _IntegrationTest
from bioagents.tests.util import ekb_kstring_from_text, ekb_from_text, get_request, agent_clj_from_text

class TestCreateTumorBoardReport(_IntegrationTest):
    def __init__(self, *args):
        super(TestCreateTumorBoardReport, self).__init__(TumorBoardModule)

    def create_message_1(self):
        content = KQMLList('CREATE-TUMOR-BOARD-REPORT')
        patient_id = '3'
        content.sets('patient-id', patient_id)

        msg = get_request(content)
        return msg, content

    def check_response_to_message_1(self, output):
        assert output.head() == 'SUCCESS', output
        genes = output.get('genes')
        print(type(genes))
        assert 'ESR1' in genes

    def create_message_failure(self):
        content = KQMLList('CREATE-TUMOR-BOARD-REPORT')
        patient_id = '101'
        content.sets('patient-id', patient_id)

        msg = get_request(content)
        return msg, content

    def check_response_to_message_failure(self, output):
        assert output.head() == 'FAILURE'
        reason = output.gets('reason')
        assert reason == 'INVALID_PID'
