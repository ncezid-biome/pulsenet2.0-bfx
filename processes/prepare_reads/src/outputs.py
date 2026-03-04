from ngs_pipeline_lib.base.file import JsonFile, ReadFile
from ngs_pipeline_lib.base.outputs import BaseOutputs


class DownloadReadsOutputs(BaseOutputs):
    def __init__(self):
        super().__init__()
        self.read_1 = ReadFile(name="read_1")
        self.read_2 = ReadFile(name="read_2")


class CleanupReadsOutputs(BaseOutputs):
    def __init__(self):
        super().__init__()
        self.fastp_report = JsonFile(name="fastp_report", compress=False)
        self.cleaned_read_1 = ReadFile(name="cleaned_read_1")
        self.cleaned_read_2 = ReadFile(name="cleaned_read_2")
