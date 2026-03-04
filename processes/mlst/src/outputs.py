#from ngs_pipeline_lib.base.file import ArchiveFile  # , TsvFile
from src.archive_bah import ArchiveFile
from ngs_pipeline_lib.base.file import TextFile
from ngs_pipeline_lib.base.outputs import BaseOutputs

class MLSTOutputs(BaseOutputs):
    def __init__(self):
        super().__init__()
        #self.novel_fastas = ArchiveFile(name="novel_fastas") # original
        self.novel_fastas = ArchiveFile(name="novel_fastas", extension=".tar") #local ArchiveFile
