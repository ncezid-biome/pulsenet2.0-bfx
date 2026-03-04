from collections import defaultdict
from csv import DictReader
from os import getenv
from pathlib import Path
from shutil import move
from ngs_pipeline_lib.base.algorithm import Algorithm
from ngs_pipeline_lib.tools.quality_control.quality_control import QualityControl
from ngs_pipeline_lib.tools.runextern import run_external
from subprocess import run, PIPE
import base64
import sys
import re
import json
import os
import binascii

from src.inputs import StxCondenserInputs
from src.outputs import StxCondenserOutputs


stx_matcher = re.compile(r'^stx(?P<xtype>\d)(?:[AB])?(?P<subtype>[a-g])?$')


class StxCondenser(Algorithm[StxCondenserInputs, StxCondenserOutputs]):

    outputs_class = StxCondenserOutputs

    def execute_stub(self):
        working_dir = Path("tmp")
        working_dir.mkdir(exist_ok=True)
        cmd = [f"{CONDA_BIN_DIR}/bowtie2"]
        process_stdout, _ = run_external(cmd, self.logger, True, working_dir)
        print(process_stdout)

    def execute_implementation(self):
        stxfile = self.inputs.stx_results
        pcrfile = self.inputs.ispcr_results
        self.results_out = self.condense(stxfile, pcrfile)
        self.set_results()

    def condense(self, stxfile: Path, pcrfile: Path):

        flnames = {
            'stxfinder': stxfile,
            'insilicopcr': pcrfile,
        }

        # We have all the results calculated since this is designed to be the last genotyper run in the set
        # Lets parse the data for the experiment type first
        expr_out = {'results': {}, 'extra': []}
        fld_out = {'results': [], 'extra': []}

        # stx types
        needed_chars = {
            'eae': False,
            'ehxA': False,
            'stx1': False,
            'stx2': False
        }

        # stx subtypes
        subtypes = {
            'stx1': [97, 99, 100],
            'stx2': range(ord('a'), ord('g') + 1)
        }

        # stx_subtypes = {
        #     'stx1a': False, 'stx1c': False, 'stx1d': False,
        #     'stx2a': False, 'stx2b': False, 'stx2c': False,
        #     'stx2d': False, 'stx2e': False, 'stx2f': False,
        #     'stx2g': False
        # }

        stx_subtypes = {}

        # Set up the needed characters for the subtypes easier than typing it out
        for gene, r in subtypes.items():
            for c in r:
                stx_subtypes[gene + chr(c)] = False

        for file in flnames.values():
            flpath = os.path.join(file)
            if not os.path.exists(flpath):
                continue
            with open(flpath, 'r') as f:
                data = f.read()

            # check for null
            if data.strip() == "null":
                self.logger.info(f"input file {file} had null contents")
                sys.exit(1)

            mod_results = json.loads(data)

            mod_brief_results = mod_results.get("results", {})

            if not mod_brief_results:
                continue

            for char, present in mod_brief_results.items():
                if not present:
                    continue
                match_obj = stx_matcher.search(char)
                if not match_obj:
                    if char in needed_chars:
                        needed_chars[char] = present
                    continue

                xtype = match_obj.group("xtype")
                subtype = match_obj.group("subtype")
                final_type = "stx"

                if xtype and final_type + xtype in needed_chars:
                    final_type += xtype
                    needed_chars[final_type] = present

                if subtype and final_type + subtype in stx_subtypes:
                    final_type += subtype
                    stx_subtypes[final_type] = present

        expr_out["results"].update(needed_chars.items())
        for k, v in stx_subtypes.items():
            if v:
                fld_out["results"].append(k)

        self.outputs.condenser_expr.content = expr_out
        self.outputs.condenser_flds.content = fld_out


    def set_results(self):
        self.result = {
            "condenser_expr": self.inputs.publish_dir + str(self.outputs.condenser_expr.path),
            "condenser_flds": self.inputs.publish_dir + str(self.outputs.condenser_flds.path)
        }

