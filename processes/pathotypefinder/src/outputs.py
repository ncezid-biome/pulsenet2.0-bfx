from ngs_pipeline_lib.base.file import JsonFile, TextFile
from ngs_pipeline_lib.base.outputs import BaseOutputs


class PathotypeFinderOutputs(BaseOutputs):
    def __init__(self):
        super().__init__()
        self._outputs.name = "outputs"
        self.pathotypefinder_geno = JsonFile(name="pathotypefinder_genotypes", compress=False)
        self.blastout = TextFile(name="blastout", compress=False)
