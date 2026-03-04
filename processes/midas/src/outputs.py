from ngs_pipeline_lib.base.file import TextFile, JsonFile
from ngs_pipeline_lib.base.outputs import BaseOutputs


class MidasOutputs(BaseOutputs):
    def __init__(self):
        super().__init__()
        self._outputs.name = "outputs"
        self.species_profile = TextFile(name="species/species_profile", compress=False)
        self.midas_result = JsonFile(name="midas_result_out", compress=False)
