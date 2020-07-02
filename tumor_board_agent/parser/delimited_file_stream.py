class DelimitedFileStream:
    # def __init__(self, delimiter='\t', trimDoubleQuotes=False):
    #     self.delimiter = delimiter
    #     self.trimDoubleQuotes = trimDoubleQuotes

    def parse_file(self, file_path, on_data, delimiter='\t'):
        with open(file_path, encoding='utf-8-sig') as f:
            for line in f:
                line = line.strip()
                tabs = line.split('\t')
                end = on_data( tabs )
                if end == True:
                    return
                # print(tabs[0].strip('\"'))
