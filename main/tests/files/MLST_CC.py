import sys
import json
import argparse
import re

## Look up CC value based on scheme and ST input. Print to outputs.json

parser = argparse.ArgumentParser()
parser.add_argument("--mlstOutput", type=str, required=True, help="outputs.json from MLST")
parser.add_argument("--mappingFile", type=str, required=True, help="ST_CC_mapping.json file containing Scheme -> ST -> CC hierarchy")
args = parser.parse_args()

# Replace scheme names with original BN scheme names
# if applicable so users are familiar
def getBNschemeName(schemeName):
    #listeria not included - BN scheme name not descriptive
    name_lookup = { 
        "campylobacter_nonjejuni_9": "MLST_fetus",
        "campylobacter": "MLST_jejuni",
        "campylobacter_nonjejuni_3": "MLST_lari",
        "campylobacter_nonjejuni_4": "MLST_upsaliensis",
        "ecoli_achtman_4": "MLST_ecoli_Achtman",
        "ecoli": "MLST_ecoli_Pasteur",
        "senterica_achtman_2": "MLST_salm_Achtman"
    }
    if schemeName in name_lookup.keys():
        return name_lookup[schemeName]
    
    return schemeName


mlst_out = json.load(open(args.mlstOutput))
cc_map = json.load(open(args.mappingFile))

for res in mlst_out["result"]["mlst_results"]:
    if res["sequence_type"] != "-":
        tempST = "ST" + re.sub('\D', '', res["sequence_type"]) #reformat ST
        keyCheck1 = cc_map.get(res["scheme"]) #check if scheme in mapping
        if keyCheck1:
            keyCheck2 = cc_map[res["scheme"]].get(res["sequence_type"]) #check if ST in mapping
            if keyCheck2:
                cc_list = cc_map[res["scheme"]][res["sequence_type"]].split("/") #more than 1 CC is separated by '/'
                for i, cc in enumerate(cc_list):

                    cc_list[i] = "ST-" + re.sub('\D', '', cc) + " complex" #remove all non digits to format correctly
                res["clonal_complex"] = '/'.join(cc_list)
                res["sequence_type"] = tempST
            else:
                res["clonal_complex"] = "NA"
                res["sequence_type"] = tempST
        else:
            res["clonal_complex"] = "NA"
            res["sequence_type"] = tempST
    else:
        res["clonal_complex"] = "-"
    
    #Update mlst_results scheme name to BN name
    res["scheme"] = getBNschemeName(res["scheme"])

#Update QC section scheme names to BN name
mlst_out["qc"]["metrics"]["autoDetectedScheme"] = getBNschemeName(mlst_out["qc"]["metrics"]["autoDetectedScheme"])
mlst_out["qc"]["metrics"]["bestScheme"] = getBNschemeName(mlst_out["qc"]["metrics"]["bestScheme"])

with open('outputs.json', 'w') as out_file:
    json.dump(mlst_out, out_file, indent=4)
