from ngs_pipeline_lib.base.file import JsonFile
from ngs_pipeline_lib.base.outputs import BaseOutputs


class ISPCROutputs(BaseOutputs):
    def __init__(self):
        super().__init__()
        self._outputs.name = "outputs"
        self.insilicopcr = JsonFile(name="insilicopcr", compress=False)
