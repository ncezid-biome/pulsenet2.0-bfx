from ngs_pipeline_lib.base.file import JsonFile, TextFile
from ngs_pipeline_lib.base.outputs import BaseOutputs


class PlasmidFinderOutputs(BaseOutputs):
    def __init__(self):
        super().__init__()
        self._outputs.name = "outputs"
        self.plasmidfinder = JsonFile(name="plasmidfinder", compress=False)
        self.blastout = TextFile(name="blastout", compress=False)
