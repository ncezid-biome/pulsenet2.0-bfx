from ngs_pipeline_lib.base.file import TextFile, JsonFile
from ngs_pipeline_lib.base.outputs import BaseOutputs


class EcStxTyperOutputs(BaseOutputs):
    def __init__(self):
        super().__init__()
        self._outputs.name = "outputs"
        self.ecstxtyper_result_out = JsonFile(name="ecstxtyper_result_out", compress=False)
