#!/usr/bin/env python3

import json
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--ispcr', required=True, help="isPCR results json")
parser.add_argument('--pf', required=True, help="PathotypeFinder results json")
parser.add_argument('--vf', required=True, help="VirulenceFinder results json")
args = parser.parse_args()


result = False

try:
    with open(args.ispcr, 'r') as ispcr_result:
        ispcr_json = json.load(ispcr_result)
        if ispcr_json["results"]["ipaH"]:
            result = True
    with open(args.pf, 'r') as pf_result:
        pf_json = json.load(pf_result)
        if "EIEC/Shigella" in pf_json["result"]["results"]:
            result = True
    with open(args.vf, 'r') as vf_result:
        vf_json = json.load(vf_result)
        if vf_json["results"]["ipaH"]:
            result = True
        if vf_json["results"]["ipaH7.8"]:
            result = True
        if vf_json["results"]["ipaH9.8"]:
            result = True
except:
    # problem loading files or results, this is usually fine
    pass

if result:
    print('shigella', end='')
else:
    print('non-shigella', end='')
