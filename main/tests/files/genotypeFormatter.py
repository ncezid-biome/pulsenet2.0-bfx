## Parse genotype results into a NARMS specific format.
## Files are parsed based on arguments used in command. If an argument is not called,
## then program will not try to parse associated file(s)

import sys
import csv
import json
import argparse

## Input arguments
parser = argparse.ArgumentParser()
parser.add_argument("--genus", required=True, help="Genus: salmonella, ecoli, campy, listeria")
parser.add_argument("--serotype", required=False, help="SALM seqsero2 serotype.json or STEC KMA ec_kma_serotyper_result_out.json")
parser.add_argument("--sero-out", required=False, help="SALM seqsero2 seqsero_outputs.json, STEC KMA kma_outputs.json, or VIBRIO serotype_outputs.json")
parser.add_argument("--plasmid", required=False, help="plasmidfinder.json from PLASMIDFINDER results")
parser.add_argument("--plasmid-out", required=False, help="pf_outputs.json")
parser.add_argument("--ispcr", required=False, help="insilicopcr.json from isPCR results")
parser.add_argument("--ispcr-out", required=False, help="ispcr_outputs.json")
parser.add_argument("--amr", required=False, help="AMRFinderResultsRaw.json from FIND_GENES results")
parser.add_argument("--amr-out", required=False, help="amr_outputs.json")
parser.add_argument("--pathotype", required=False, nargs=2, help="Pathotype_Finder outputs.json and pathotypefinder_genotypes.json (in this order)")
parser.add_argument("--pathotype-out", required=False, help="pathotype_outputs.json")
parser.add_argument("--stxtyper", required=False, help="STEC STX_TYPER ecstxtyper_result_out.json")
parser.add_argument("--stx-out", required=False, help="stx_outputs.json")
parser.add_argument("--stxcondenser", required=False, help="STEC STX_CONDENSER master_input.json created from stx_condenser_expr.json,stx_condenser_flds.json (in this order)")
parser.add_argument("--virulence", required=False, help="VirulenceFinder virulence.json")
parser.add_argument("--virulence-out", required=False, help="virulence_outputs.json")
parser.add_argument("--ctoxin", required=False, help="Cbot_ToxinTyper cbot_toxintyper.json")
parser.add_argument("--cbot-toxin-out", required=False, help="cbot_toxin_outputs.json")
args = parser.parse_args()

#Need the official algorithm versions - still required after narms updating?

version_dicts = {
    "salmonella.serotype.version.json": {
        "module": "",
        "algorithm": "seqsero2s_v1.1.1",
        "database": "seqsero2s_v1.1.1"
    },
    "ecoli.serotype.version.json": {
        "module": "",
        "algorithm": "kma_v1.4.14",
        "database": "1.0.0"
    },
    "vibrio.serotype.version.json": {
        "module": "",
        "algorithm": "blast_v2.12",
        "database": "1.0.0"
    },
    "plasmids.version.json": {
        "module": "",
        "algorithm": "blast_v2.12",
        "database": "1.1.1"
    },
    "insilicopcr.version.json": {
        "module": "",
        "algorithm": "blast_v2.12",
        "database": "blast_v2.12"
    },
    "resistance.version.json": {
        "module": "",
        "algorithm": "ncbi-amrfinderplus_v4.0.15",
        "database": "2025-03-25.1"
    },
    "ecoli.pathotype.version.json": {
        "module": "",
        "algorithm": "blast_v2.12",
        "database": "1.0.0"
    },
    "ecoli.stx.version.json": {
        "module": "",
        "algorithm": "bowtie2_v2.2.5",
        "database": "1.0.0"
    },
    "virulence.version.json": {
        "module": "",
        "algorithm": "blast_v2.12",
        "database": "1.0.0"
    },
    "vibrio.virulence.version.json": {
        "module": "",
        "algorithm": "blast_v2.12",
        "database": "1.0.0"
    },
    "cbot.virulence.version.json": {
        "module": "",
        "algorithm": "blast_v2.12",
        "database": "1.0.0"
    }
}

empty_version = {
    "module": "",
    "algorithm": "",
    "database": ""
}


## Check for Genus 
if args.genus is None:
    sys.stderr.write(f"Error: {args.genus} required.\n")
    sys.exit(1)



## Start parsing though files and add to output_json dictionary
output_json = {}

# SALM seqsero2 and STEC kma serotyping files
# have same format and can be parsed the same
if args.serotype is not None:
    with open(args.serotype, 'r') as serotype_file:
        data = json.load(serotype_file)
        key = f"{args.genus}.serotype.version.json"
        alg_db_ver = version_dicts.get(key, empty_version)

        output_json[key] = alg_db_ver
        output_json[f"futures.1.{key}"] = alg_db_ver

        if args.genus == "vibrio":
            output_json[f"{args.genus}.serotype.json"] = data[1]
            output_json[f"futures.1.{args.genus}.serotype.json"] = data[1]

        else:
            output_json[f"{args.genus}.serotype.json"] = data
            output_json[f"futures.1.{args.genus}.serotype.json"] = data

        if args.genus == "ecoli":
            antigens = data["extra"]

            for antigen in antigens:
                
                antType = antigen["#Template"].split("_")[-1]
                if "O" in antType:
                    antigen.update({"Otype": antType})
                elif "H" in antType:
                    antigen.update({"Htype": antType})

    with open(args.sero_out, 'r') as ppo_file:
        data = json.load(ppo_file)["report"]
        output_json[f"{args.genus}.serotype.version.json"]["module"] = f'{data["name"]}_{data["version"]}'

# plasmid finder results
if args.plasmid is not None:
    with open(args.plasmid, 'r') as plasmid_file:
        data = json.load(plasmid_file)
        output_json["plasmids.version.json"] = version_dicts["plasmids.version.json"]
        output_json["plasmids.json"] = data
    
    with open(args.plasmid_out, 'r') as ppo_file:
        data = json.load(ppo_file)["report"]
        output_json["plasmids.version.json"]["module"] = f'{data["name"]}_{data["version"]}'

# SALM and STEC - isPCR results
if args.ispcr is not None:
    with open(args.ispcr, 'r') as ispcr_file:
        data = json.load(ispcr_file)
        output_json["insilicopcr.version.json"] = version_dicts["insilicopcr.version.json"]
        output_json["insilicopcr.json"] = data
        if args.genus == "salmonella":
            output_json["futures.1.insilicopcr.version.json"] = version_dicts["insilicopcr.version.json"]
            output_json["futures.1.insilicopcr.json"] = data
    
    with open(args.ispcr_out, 'r') as ppo_file:
        data = json.load(ppo_file)["report"]
        output_json["insilicopcr.version.json"]["module"] = f'{data["name"]}_{data["version"]}'

# AMRFinder results (Find Genes)
if args.amr is not None:
    with open(args.amr, 'r') as resistance_file:
        data = json.load(resistance_file)
        output_json["resistance.version.json"] = version_dicts["resistance.version.json"]
        output_json["resistance.json"] = data
        output_json["full.resistance.json"] = data
        output_json["resistance.antibios.json"] = data
        output_json["resistance.point.json"] = data
    with open(args.amr_out, 'r') as ppo_file:
        data = json.load(ppo_file)["report"]
        output_json["resistance.version.json"]["module"] = f'{data["name"]}_{data["version"]}'

# STEC PathotypeFinder results
if args.pathotype is not None:
    with open(args.pathotype[0], 'r') as p1file:
        data = json.load(p1file)
        output_json[f"{args.genus}.pathotype.version.json"] =  version_dicts.get(f"{args.genus}.pathotype.version.json", empty_version)
        output_json[f"{args.genus}.pathotype_pathotypes.json"] = data["result"]
    
    with open(args.pathotype[1], 'r') as p2file:
        data = json.load(p2file)
        output_json[f"{args.genus}.pathotype_genotypes.json"] = data
    
    with open(args.pathotype_out, 'r') as ppo_file:
        data = json.load(ppo_file)["report"]
        output_json[f"{args.genus}.pathotype.version.json"]["module"] = f'{data["name"]}_{data["version"]}'

# STEC stx typer
if args.stxtyper is not None:
    with open(args.stxtyper, 'r') as sfile:
        data = json.load(sfile)
        output_json[f"{args.genus}.stx.version.json"] = version_dicts.get(f"{args.genus}.stx.version.json", empty_version)
        output_json[f"{args.genus}.stxfinder_genotypes.json"] = data
    
    with open(args.stx_out, 'r') as ppo_file:
        data = json.load(ppo_file)["report"]
        output_json[f"{args.genus}.stx.version.json"]["module"] = f'{data["name"]}_{data["version"]}'

# STEC stx condenser
if args.stxcondenser is not None:
    with open(args.stxcondenser, 'r') as s1file:
        data = json.load(s1file)
        output_json["stx_condenser_expr.json"] = data[0]
        output_json["stx_condenser_flds.json"] = data[1]

# VirulenceFinder results
if args.virulence is not None:
    with open(args.virulence, 'r') as vfile:
        data = json.load(vfile)

        output_json["virulence.version.json"] = version_dicts["virulence.version.json"]
        if args.genus == "vibrio":
            output_json["virulence.json"] = data[0]

        else:
            output_json["virulence.json"] = data
    
    with open(args.virulence_out, 'r') as ppo_file:
        data = json.load(ppo_file)["report"]
        output_json["virulence.version.json"]["module"] = f'{data["name"]}_{data["version"]}'

# Cbot toxin typing results
if args.ctoxin is not None:
    with open(args.ctoxin, 'r') as cfile:
        data = json.load(cfile)
        output_json[f"{args.genus}.toxintyping.json"] = data

    with open(args.cbot_toxin_out, 'r') as ppo_file:
        data = json.load(ppo_file)["report"]
        output_json[f"{args.genus}.toxintyping.version.json"] = version_dicts["cbot.virulence.version.json"]
        output_json[f"{args.genus}.toxintyping.version.json"]["module"] = f'{data["name"]}_{data["version"]}'

## Write Results
with open('GenotypingResult.json', 'w') as output_file:
    json.dump(output_json, output_file, indent=4)