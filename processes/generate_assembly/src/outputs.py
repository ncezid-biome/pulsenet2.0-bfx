from ngs_pipeline_lib.base.file import AssemblyFile
from ngs_pipeline_lib.base.outputs import BaseOutputs


class GenerateAssemblyOutputs(BaseOutputs):
    def __init__(self):
        super().__init__()
        self.assembly = AssemblyFile(name="assembly")
