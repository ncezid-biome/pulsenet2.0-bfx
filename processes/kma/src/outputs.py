from ngs_pipeline_lib.base.file import TextFile, JsonFile
from ngs_pipeline_lib.base.outputs import BaseOutputs


class EcKMASerotyperOutputs(BaseOutputs):
    def __init__(self):
        super().__init__()
        self._outputs.name = "outputs"
        self.ec_kma_serotyper_result_out = JsonFile(name="ec_kma_serotyper_result_out", compress=False)