from ngs_pipeline_lib.base.file import JsonFile
from ngs_pipeline_lib.base.outputs import BaseOutputs


class ANIOutputs(BaseOutputs):
    def __init__(self):
        super().__init__()
        self._outputs.name = "outputs"
        self.best_hit = JsonFile(name="best_hit", compress=False)
