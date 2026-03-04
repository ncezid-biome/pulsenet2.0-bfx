import logging
import shutil
import tarfile
from abc import ABC, abstractmethod
from csv import writer as csv_writer
from dataclasses import dataclass, field
from io import StringIO
from json import dumps
from pathlib import Path
from typing import Any

from ngs_pipeline_lib.base.file import File

class ArchiveFile(File):
    content: list[Path] = field(default_factory=list)
    extension: str = ".tar"

    def to_file(self):        
        with tarfile.open(self.path, "w") as f:
            if self.content is not None: ## added to fix fail when no novel sequences
                for path in self.content:
                    if path.exists():
                        f.add(path)
                        if path.is_dir():
                            shutil.rmtree(path)
                        else:
                            path.unlink()
                    else:
                        raise ValueError(
                            f"Missing file while building the archive file {self.path}"
                        )