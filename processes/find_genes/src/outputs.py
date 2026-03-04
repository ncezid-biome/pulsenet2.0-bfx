from ngs_pipeline_lib.base.file import TextFile, JsonFile
from ngs_pipeline_lib.base.outputs import BaseOutputs


class FindGenesOutputs(BaseOutputs):
    def __init__(self):
        super().__init__()

        self.find_genes = TextFile(name="find_genes", extension=".tsv", compress=False)
        self.find_genes_json = JsonFile(name= "AMRFinderResultsRaw", compress = False)
