from ngs_pipeline_lib.base.file import JsonFile
from ngs_pipeline_lib.base.outputs import BaseOutputs


class SeqSeroOutputs(BaseOutputs):
    def __init__(self):
        super().__init__()
        self._outputs.name = "outputs"
        self.serotype = JsonFile(name="serotype", compress=False)
