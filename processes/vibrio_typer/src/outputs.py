from ngs_pipeline_lib.base.file import TextFile, JsonFile
from ngs_pipeline_lib.base.outputs import BaseOutputs


class VibrioVirulenceOutputs(BaseOutputs):
    def __init__(self):
        super().__init__()
        self._outputs.name = "outputs"
        self.vibrio_virulence_pcr = JsonFile(name="vibrio_virulence_pcr", compress=False)
        self.vibrio_virulence_results_out = JsonFile(name="vibrio_virulence_results", compress=False)
