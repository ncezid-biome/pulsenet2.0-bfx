from ngs_pipeline_lib.base.file import JsonFile, TextFile
from ngs_pipeline_lib.base.outputs import BaseOutputs


class CbotToxinTypingOutputs(BaseOutputs):
    def __init__(self):
        super().__init__()
        self._outputs.name = "outputs"
        self.cbot_toxintyper = JsonFile(name="cbot_toxintyper", compress=False)
        #self.blastout = TextFile(name="blastout", compress=False)

class ISPCROutputs(BaseOutputs):
    def __init__(self):
        super().__init__()
        self._outputs.name = "outputs"
        self.insilicopcr = JsonFile(name="insilicopcr", compress=False)