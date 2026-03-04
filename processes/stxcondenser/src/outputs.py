from ngs_pipeline_lib.base.file import TextFile, JsonFile
from ngs_pipeline_lib.base.outputs import BaseOutputs


class StxCondenserOutputs(BaseOutputs):
    def __init__(self):
        super().__init__()
        self._outputs.name = "outputs"
        self.condenser_expr = JsonFile(name="stx_condenser_expr", compress=False)
        self.condenser_flds = JsonFile(name="stx_condenser_flds", compress=False)
