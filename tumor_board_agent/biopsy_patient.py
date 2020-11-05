import os.path
from .biopsy_sample import BiopsySample

class BiopsyPatient:
    def __init__(self, patient_id, patient_folder):
        self.patient_id = patient_id
        self.sample = BiopsySample();
        self.clinical_info = {}
        muts_cnas = get_mutations_and_cnas(patient_id, patient_folder)
        self.mutations = muts_cnas['mutations']
        self.cnas = muts_cnas['cnas']

def get_mutations_and_cnas(patient_id, patient_folder):
    file_address = os.path.join(patient_folder, patient_id)
    # assert os.path.is_file(file_address)

    content = read_text_file(file_address)
    lines = content.splitlines()

    ret_val = {'mutations': [], 'cnas': []}

    for line in lines:
        obj = {}
        tabs = line.split('\t')
        val = tabs[1]
        # lo_val = val.lower()
        obj['geneName'] = tabs[0]
        key = 'proteinChange'
        arr_name = 'mutations'

        # if lo_val == 'del' or lo_val == 'amp' or lo_val == 'neu':
        if val == '-2' or val == '-1' or val == '0' or val == '1' or val == '2':
            key = 'alteration'
            arr_name = 'cnas'
            val = int(val)

        obj[key] = val

        ret_val[arr_name].append(obj)

    return ret_val

# TODO: put it to a utility file
def read_text_file(f):
    with open(f, 'r') as content_file:
        content = content_file.read()
        return content
