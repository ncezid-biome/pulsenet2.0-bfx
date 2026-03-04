import csv, json
from logging import Logger
from os import getenv
from pathlib import Path
from tarfile import open as tar_open
import re
import subprocess

from ngs_pipeline_lib.base.algorithm import Algorithm
from ngs_pipeline_lib.tools.runextern import run_external
from ngs_pipeline_lib.tools.tools import gunzip_file

from src.inputs import FindGenesInputs
from src.outputs import FindGenesOutputs
from src.stub import get_stub_content

CONDA_BIN_DIR = getenv("CONDA_BIN_DIR", "/opt/conda/condabin")


class FindGenes(Algorithm[FindGenesInputs, FindGenesOutputs]):
    outputs_class = FindGenesOutputs

    def set_results(self):
        genes = []
        try:
            with open(self.outputs.find_genes.path) as file:
                tsv_file = csv.reader(file, delimiter="\t")

                headers = next(tsv_file)

                for line in tsv_file:
                    if len(line) == len(headers):
                        gene = {}
                        for index, value in enumerate(line):
                            gene[headers[index]] = value
                        genes.append(gene)
        except Exception as e:
            genes.append({"error": str(e)})

        self.result = {"genes": genes}

        #parse tsv and print to json
        results_out = { 'results': {}, 'extra': [] }
        with open(self.outputs.find_genes.path, 'r', newline='') as tfile:
            reader = csv.DictReader(tfile, delimiter='\t')
            for row in reader:
                gene = row['Element symbol']
                if gene not in results_out:
                    results_out['results'][gene] = True
                hit_information = [
                    {
                        'hits': [
                            {
                                'Contig id' : row['Contig id'],
                                'Start' : row['Start'],
                                'Stop' : row['Stop'],
                                'Strand' : row['Strand'],
                                'Target length' : row['Target length'],
                                'Reference sequence length' : row['Reference sequence length'],
                                'Alignment length' : row['Alignment length'],
                            }
                        ],
                        'Element symbol' : row['Element symbol'],
                        'Protein identifier' : row['Protein id'],
                        'Sequence name' : row['Element name'],
                        'Scope' : row['Scope'],
                        'Element type' : row['Type'],
                        'Element subtype' : row['Subtype'],
                        'Class' : row['Class'],
                        'Subclass' : row['Subclass'],
                        'Method' : row['Method'],
                        '% Coverage of reference sequence' : row['% Coverage of reference'],
                        '% Identity to reference sequence' : row['% Identity to reference'],
                        'Accession of closest sequence' : row['Closest reference accession'],
                        'Name of closest sequence' : row['Closest reference name'],
                        'HMM id' : row['HMM accession'],
                        'HMM description' : row['HMM description'],
                    }
                ]
                results_out['extra'].extend(hit_information)

        with open(self.outputs.find_genes_json.path, 'w') as jfile:
            json.dump(results_out, jfile, indent=4)


    def execute_stub(self):
        self.outputs.find_genes.content = get_stub_content()
        self.outputs.find_genes.to_file()
        self.set_results()

    def execute_implementation(self):
        organism = str(self.inputs.organism.genus).capitalize()
        species = str(self.inputs.organism.species).capitalize()
        assembly = gunzip_file(self.inputs.assembly)
        self.call_amrfinder(
            self.outputs.find_genes.path,
            organism,
            species,
            assembly,
            self.inputs.n_threads,
            self.logger,
        )

        self.set_results()

    @staticmethod
    def prepare_db(compressed_db: Path) -> Path:
        destination = Path("amrfinder_db")
        destination.mkdir(exist_ok=True)
        with tar_open(compressed_db) as reader:
            reader.extractall(path=destination)
        return destination

    @staticmethod
    def call_amrfinder(
        output_file: Path,
        organism: str,
        species: str,
        assembly: Path,
        n_threads: int,
        logger: Logger,
    ) -> dict[str, str]:
        command_line = [f"{CONDA_BIN_DIR}/mamba", "run", "amrfinder"]
        command_line = ["amrfinder"]
        command_line += ["--nucleotide", str(assembly)]
        command_line += ["--output", str(output_file)]
        command_line += ["--threads", str(n_threads)]

        ## check with organisms are supported by amrfinder 
        res = subprocess.run("amrfinder -l", shell =True, text=True, capture_output=True)
        org = re.sub("_", " " ,re.sub("\n", "", re.sub("Available --organism options: ", "", res.stdout)))
        org = list(map(lambda a: a.upper(), org.split(", ")))

        ## define the genus and species 
        genus_species=organism +" "+ species.lower()

        ## add the organism flag to the command 
        if(str(organism).upper() in org): 
            command_line += ["--organism", str(organism)]
        elif(str(genus_species).upper() in org): 
            command_line += ["--organism", str(genus_species)]
        std_out, _ = run_external(command_line, logger)
        logger.info(std_out)
