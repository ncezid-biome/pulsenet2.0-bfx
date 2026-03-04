from ngs_pipeline_lib.base.file import AssemblyFile, TextFile
from ngs_pipeline_lib.base.outputs import BaseOutputs


class CorrectAssemblyOutputs(BaseOutputs):
    def __init__(self):
        super().__init__()
        self.depth_contigs = TextFile(
            name="depth_contigs", extension=".tsv", compress=False
        )
        self.corrected_assembly = AssemblyFile(name="corrected_assembly")
        self.corrected_alignment = TextFile(
            name="corrected_alignment", extension=".cram", compress=False
        )
