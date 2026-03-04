#!/usr/bin/env python3

import sys
import json
import argparse

class BaseMaster:
    def __init__(self):
        self.data = {
            "Summary": {
                "Genus": "",
                "Species": "",
                "MidasSecondaryCoverage": "",
                "AverageQuality_Qscore": "",
                "AvgQuality_read1": "",
                "AvgQuality_read2": "",
                "AvgReadLength": "",
                "AvgDeNovoCoverage": "",
                "Length": "",
                "Contigs": "",
                "N50": "",
                "AllelesCalledCoreCount": "",
                "AllelesCalledCorePercent": "",
                "AMRFinderGenes": "",
                "PlasmidsFound": "",
            },
            "GenusIdentify": {},
            "PrepareReads": {},
            "MIDAS": {},
            "GenerateAssembly": {},
            "CorrectAssembly": {},
            "ANI": {},
        }
        self._process_parsing = {
            "genusidentify": self.parseGenusResults,
            "preparereads": self.parseReadResults,
            "midas": self.parseMidasResults,
            "generateassembly": self.parseAssemblyResults,
            "correctassembly": self.parseCorrAssemblyResults,
            "ani": self.parseANIResults,
            "org_ani": self.parseOrgANIResults,
            "allelecalling": self.parseAlleleCallResults,
            "allelefilter": self.parseAlleleFilterResults,
            "allelenaming": self.parseAlleleNamingResults,
            "mlst": self.parseMLST,
            "ispcr": self.parseIsPCRResults,
            "seqsero": self.parseSeqSeroResults,
            "amrfinder": self.parseAMRFinder,
            "plasmidfinder": self.parsePlasmidFinder,
            "pathotypefinder": self.parsePathotypeFinder,
            "kma": self.parseKMA,
            "stecgroup": self.addSTECGroup,
            "stxtyper": self.parseStxTyper,
            "stxcondenser": self.parseStxCondenser,
            "virulencefinder": self.parseVirulenceFinder,
            "cbot_toxintyper": self.parseCbotToxinTyping ,
            "vibriovirulence": self.parseVibrioVirulence,
            "vibrioserotype": self.parseVibrioSerotype,
            "shigapass": self.parseShigaPass,
            "shigeifinder": self.parseShigEiFinder,

        }

    def _add_metadata(self):
        # Run this at the end of init to keep it at the end of the dict
        self.data.update(
            {
                "metadata": {
                    "GENUS": "",
                    "SPECIES": "",
                    "Subspecies": "",
                    "ANI_Score": "",
                    "ANI_Percent_Aligned": "",
                    "MLST_ST": "",
                    "MLST_CC": "",
                },
            }
        )
    
    def to_file(self):
        with open("PipelineProcessOutputs.json", "w") as output_file:
            json.dump(self.data, output_file, indent=4)

    @classmethod
    def from_json(cls, json_file):
        x = cls()
        with open(json_file) as fin:
            x.data.update(json.load(fin))
        return x

    def add_process(self, process, parseFile=None, parseAddition=None):
        parse_process = self._process_parsing[process.lower()]
        args = {"parseFile": parseFile, "parseAddition": parseAddition}
        parse_process(**args)

    def parseGenusResults(self, parseFile, **_):
        with open(parseFile, "r") as genusID:
            data = json.load(genusID)
        if data["qc"]["metrics"]["primaryCoverage"]["coverage"] > 0:
            genus = data["qc"]["metrics"]["primaryCoverage"]["genus"]
            self.data["GenusIdentify"]["result"] = "PASS"
            self.data["GenusIdentify"]["metrics"] = {"Genus": genus}
            self.data["GenusIdentify"]["issues"] = []

            # populate fields
            self._add_fields(genus)

        else:
            genus = "NA"
            self.data["GenusIdentify"]["result"] = data["qc"]["result"]
            self.data["GenusIdentify"]["metrics"] = {"Genus": genus}
            self.data["GenusIdentify"]["issues"] = data["qc"]["issues"]

    def parseReadResults(self, parseFile, **_):
        with open(parseFile, "r") as readRes:
            data = json.load(readRes)
        self.data["PrepareReads"] = data["qc"]
        self.data["Summary"]["AverageQuality_Qscore"] = data["qc"]["metrics"][
            "trimmedCombinedAvgQC"
        ]
        self.data["Summary"]["AvgQuality_read1"] = data["qc"]["metrics"][
            "trimmedRead1AvgQC"
        ]
        self.data["Summary"]["AvgQuality_read2"] = data["qc"]["metrics"][
            "trimmedRead2AvgQC"
        ]
        self.data["Summary"]["AvgReadLength"] = data["qc"]["metrics"]["trimmedCombinedReadMeanLength"]

    def parseMidasResults(self, parseFile, parseAddition):
        with open(parseFile, "r") as midasOut:
            data = json.load(midasOut)
        # make sure secondary count is included even if there is no secondary organism
        data["qc"]["metrics"]["secondaryCoverage"]["count"] = data["qc"]["metrics"]["secondaryCoverage"].get("count", 0)

        self.data["Summary"]["MidasSecondaryCoverage"] = data["qc"]["metrics"][
            "secondaryCoverage"
        ]["coverage"]
        self.data["MIDAS"] = data["qc"]

        # Flatten Midas Output
        newData = {}
        for key, value in data["qc"]["metrics"].items():
            for subKey, subValue in value.items():
                newKey = key + "_" + subKey
                newData.update({newKey: subValue})
        self.data["MIDAS"]["metrics"] = newData

        # Add species_id
        with open(parseAddition, "r") as midasRes:
            data = json.load(midasRes)
        self.data["MIDAS"]["metrics"]["species_id"] = data["extra"][0]["species_id"]

    def parseAssemblyResults(self, parseFile, **_):
        with open(parseFile, "r") as assemRes:
            data = json.load(assemRes)
        self.data["GenerateAssembly"] = data["qc"]

    def parseCorrAssemblyResults(self, parseFile, **_):
        with open(parseFile, "r") as corrAssemRes:
            data = json.load(corrAssemRes)
        self.data["CorrectAssembly"] = data["qc"]
        self.data["Summary"]["AvgDeNovoCoverage"] = data["qc"]["metrics"][
            "cleanedAverageDepth"
        ]
        self.data["Summary"]["Length"] = data["qc"]["metrics"]["cleanedLength"]
        self.data["Summary"]["Contigs"] = data["qc"]["metrics"]["cleanedNContigs"]
        self.data["Summary"]["N50"] = data["qc"]["metrics"]["cleanedN50"]

    def parseANIResults(self, parseFile, parseAddition):
        with open(parseFile, "r") as aniRes:
            data = json.load(aniRes)
        self.data["ANI"] = data["qc"]
        self.data["ANI"][
            "metrics"
        ].clear()  # duplicate info (key names slightly different)

        with open(parseAddition, "r") as aniBest:
            data = json.load(aniBest)

        self.data["ANI"]["metrics"].update(data)
        self.data["Summary"]["Genus"] = data["Genus"]
        self.data["Summary"]["Species"] = data["Species"]
        self.data["metadata"]["GENUS"] = data["Genus"]
        self.data["metadata"]["SPECIES"] = data["Species"]
        self.data["metadata"]["Subspecies"] = data["Subspecies"]
        self.data["metadata"]["ANI_Score"] = data["ANI_Score"]
        self.data["metadata"]["ANI_Percent_Aligned"] = data["Percent_Aligned"]

    def parseOrgANIResults(self, parseFile, parseAddition):
        with open(parseFile, "r") as aniRes2:
            data = json.load(aniRes2)
        self.data["ORG_ANI"] = data["qc"]
        self.data["ORG_ANI"][
            "metrics"
        ].clear()  # duplicate info (key names slightly different)

        with open(parseAddition, "r") as aniBest2:
            data = json.load(aniBest2)
        self.data["ORG_ANI"]["metrics"].update(data)
        self.data["metadata"]["Subspecies"] = data["Subspecies"]

    def parseAlleleCallResults(self, parseFile, **_):
        with open(parseFile, "r") as ACout:
            data = json.load(ACout)
        self.data["AlleleCalling"] = data["qc"]

    def parseAlleleFilterResults(self, parseFile, **_):
        with open(parseFile, "r") as AFout:
            data = json.load(AFout)
        self.data["AlleleFilter"] = data["qc"]
        self.data["Summary"]["AllelesCalledCoreCount"] = data["qc"]["metrics"][
            "coreCount"
        ]
        self.data["Summary"]["AllelesCalledCorePercent"] = data["qc"]["metrics"][
            "corePercentage"
        ]

        #Subspecies no longer needed in listeria after allele filtering
        #Change to lineage
        if self.data["metadata"]["GENUS"] == "Listeria":
            self.data["ANI"]["metrics"]["Lineage"] = self.data["metadata"]["Subspecies"]
            self.data["metadata"]["Lineage"] = self.data["metadata"]["Subspecies"]
            del self.data["ANI"]["metrics"]["Subspecies"]
            del self.data["metadata"]["Subspecies"]
        elif self.data["metadata"]["GENUS"] == "Vibrio":
            del self.data["metadata"]["Subspecies"]


    def parseAlleleNamingResults(self, parseFile, **_):
        with open(parseFile, "r") as ANout:
            data = json.load(ANout)
        self.data["AlleleNaming(wgmlst)"] = data["qc"]
        self.data["AlleleNaming(wgmlst)"].update(
            {"processResults": {"adx_analysis_section": "wgmlst"}}
        )

    def parseMLST(self, parseFile, **_):
        with open(parseFile, "r") as mlst:
            data = json.load(mlst)
        self.data["MLST"] = data["qc"]
        self.data["MLST"].update(data["result"])

        #Add MLST metadata
        ctr = 0
        for result in data["result"]["mlst_results"]:
            #Set preferred scheme if available
            if "MLST_ecoli_Achtman" in result["scheme"]:
                self.data["metadata"]["MLST_ST"] = result["sequence_type"]
                self.data["metadata"]["MLST_CC"] = result["clonal_complex"]
                ctr = 1 #reset counter
                break
            
            #Check Results
            if result["sequence_type"] != "-":
                self.data["metadata"]["MLST_ST"] = result["sequence_type"]
                self.data["metadata"]["MLST_CC"] = result["clonal_complex"]
                ctr = ctr + 1

        # Check if more than 1 ST
        if ctr > 1:
            self.data["metadata"]["MLST_ST"] = "Need to add scheme preference. Multiple schemes returned sequence types."
            self.data["metadata"]["MLST_CC"] = ""

    def parseIsPCRResults(self, parseFile, parseAddition):
        with open(parseFile, "r") as isPCR:
            data = json.load(isPCR)
        self.data["isPCR"] = data["qc"]

        with open(parseAddition, "r") as isPCR:
            data = json.load(isPCR)
        self.data["isPCR"]["genotypingResults"] = data

    def parseSeqSeroResults(self, parseFile, parseAddition):
        with open(parseFile, "r") as sSero:
            data = json.load(sSero)
        self.data["SeqSero"] = data["qc"]

        with open(parseAddition, "r") as sSero:
            data = json.load(sSero)
        self.data["SeqSero"]["genotypingResults"] = data

        if data is not None:
            self.data["metadata"]["ANTIGENFORM_WGS"] = data["results"]["formula"]
            self.data["metadata"]["SEROTYPE_WGS"] = data["results"]["serotype"]
        else:
            self.data["metadata"]["ANTIGENFORM_WGS"] = ""
            self.data["metadata"]["SEROTYPE_WGS"] = ""

    def parseAMRFinder(self, parseFile, parseAddition):
        with open(parseFile, "r") as AMR:
            data = json.load(AMR)
        self.data["AMRFinder(resistance)"] = data["qc"]

        with open(parseAddition, "r") as AMR2:
            data = json.load(AMR2)

        oldVal = data["results"]
        newVal = {
            key: "present" if value else "absent" for key, value in oldVal.items()
        }
        data["results"] = newVal
        self.data["AMRFinder(resistance)"]["genotypingResults"] = data

        self.data["Summary"]["AMRFinderGenes"] = len(data["results"])

    def parsePlasmidFinder(self, parseFile, parseAddition):
        with open(parseFile, "r") as PFout:
            data = json.load(PFout)
        self.data["PlasmidFinder(plasmids)"] = data["qc"]

        with open(parseAddition, "r") as PF:
            data = json.load(PF)
        self.data["Summary"]["PlasmidsFound"] = list(data["results"].values()).count(
            True
        )

        oldVal = data["results"]
        newVal = {
            key: "present" if value else "absent" for key, value in oldVal.items()
        }
        data["results"] = newVal
        self.data["PlasmidFinder(plasmids)"]["genotypingResults"] = data

    def parsePathotypeFinder(self, parseFile, **_):
        with open(parseFile, "r") as patho:
            data = json.load(patho)
        self.data["PathotypeFinder"] = data["qc"]
        self.data["PathotypeFinder"]["genotypingResults"] = data["result"]
        self.data["metadata"]["Pathotype"] = "; ".join(data["result"]["results"])

    def parseKMA(self, parseFile, parseAddition):
        with open(parseFile, "r") as kma_out:
            data = json.load(kma_out)
        self.data["KMA"] = data["qc"]

        with open(parseAddition, "r") as kma_sero:
            data = json.load(kma_sero)
        self.data["KMA"]["genotypingResults"] = data

        sero = ""
        for oRes in data["results"]["O"]:
            sero += oRes
        sero += ":"
        for hRes in data["results"]["H"]:
            sero += hRes
        self.data["metadata"]["SEROTYPE_WGS"] = sero

    def parseShigaPass(self, parseFile, parseAddition):
        with open(parseFile, "r") as shigapass_out:
            shigapass_result = shigapass_out.readline()
        self.data["ShigaPass"] = {
            'result': 'PASS',
            'metrics': {
                'serotype': shigapass_result
            },
            'issues': []
        }

    def parseShigEiFinder(self, parseFile, parseAddition):
        with open(parseFile, "r") as shigeifinder_out:
            shigeifinder_result = shigeifinder_out.readline()
        self.data["ShigEiFinder"] = {
            'result': 'PASS',
            'metrics': {
                'serotype': shigeifinder_result
            },
            'issues': []
        }

    def addSTECGroup(self, **_):
        def getOrgSTEC(pathotype, serotype_wgs, shigella_serotype):
            # parse serotype into O and H group
            if ":" in serotype_wgs:
                o_group, h_group = serotype_wgs.split(":")
            else:
                o_group, h_group = "", ""

            # Shigella
            shigella_abreviations = {
                'SS': 'Shigella sonnei',
                'SF': 'Shigella flexneri',
                'SD': 'Shigella dysenteriae',
                'SB': 'Shigella boydii'
            }
            if shigella_serotype:
                if shigella_serotype[0:2] in shigella_abreviations.keys():
                    # edge case multiple serotypes, try to find superset, ex SF4av/SF4b/SF4bv -> SF4
                    if '/' in shigella_serotype:
                        shigella_serotypes_list = [i[2:] for i in shigella_serotype.split('/')]
                        serotype_superset = ''
                        for i in range(min([len(j) for j in shigella_serotypes_list])):
                            if len(set([j[i] for j in shigella_serotypes_list])) > 1:
                                break
                            serotype_superset += shigella_serotypes_list[0][i]
                        organism = shigella_abreviations[shigella_serotype[0:2]]
                        serotype = serotype_superset
                        return organism,serotype
                    elif 'Unknown' in shigella_serotype:
                        organism = shigella_abreviations[shigella_serotype[0:2]]
                        serotype = ""
                        return organism,serotype
                    else:
                        organism = shigella_abreviations[shigella_serotype[0:2]]
                        serotype = shigella_serotype[2:]
                        return organism,serotype
                elif shigella_serotype.startswith('EIEC'):
                    organism = 'EIEC'
                    # cut off "EIEC "
                    serotype = shigella_serotype[5:]
                    return organism,serotype
                elif shigella_serotype[0] == '0':
                    organism = 'EIEC'
                    # cut off "0"
                    serotype = shigella_serotype[1:]
                    return organism,serotype
                # if the above conditions aren't met then continue to STEC check
                # this means its not actually shigella

            # STEC
            if any("STEC" in ptypes for ptypes in pathotype):
                if serotype_wgs in ("O157:H7", ":H7"):
                    return "STEC O157", ""
                if "O157" in o_group or "H7" in h_group:
                    return "CHECK STEC O157", ""
                return "STEC NonO157", ""

            # other pathotype
            if len(pathotype) > 0:
                return "Non STEC EC", ""

            # No pathotype
            return "Non STEC EC", ""

        try:
            serotype_wgs = self.data["metadata"]["SEROTYPE_WGS"]  # string
            pathotype = self.data["PathotypeFinder"]["genotypingResults"][
                "results"
            ]  # list
        except:
            print("Error: SEROTYPE_WGS and PathotypeFinder results required")
            exit(1)
        try:
            shigella_serotype = self.data["ShigEiFinder"]['metrics']['serotype']
        except:
            # print('no ShigEiFinder result, trying ShigaPass')
            try:
                shigella_serotype = self.data["ShigaPass"]['metrics']['serotype']
            except:
                # print('no ShigaPass result, this is Ecoli')
                shigella_serotype = False

        organism,serotype = getOrgSTEC(pathotype=pathotype, serotype_wgs=serotype_wgs, shigella_serotype=shigella_serotype)
        self.data["metadata"]["Escherichia_group"] = organism
        if shigella_serotype == "SS": 
            self.data["metadata"]["SEROTYPE_WGS"] = serotype
        elif serotype: 
            self.data["metadata"]["SEROTYPE_WGS"] = serotype


    def parseStxTyper(self, parseFile, parseAddition):
        with open(parseFile, "r") as stxout:
            data = json.load(stxout)
        self.data["STXtyper"] = data["qc"]

        with open(parseAddition, "r") as stxRes:
            data = json.load(stxRes)
        self.data["STXtyper"]["genotypingResults"] = data

    def parseStxCondenser(self, parseFile, parseAddition):
        with open(parseFile, "r") as stxout:
            data = json.load(stxout)
        self.data["STXcondenser"] = data["qc"]

        with open(parseAddition, "r") as stxexpr:
            data = json.load(stxexpr)
        self.data["STXcondenser"]["genotypingResults"] = {
            "stx_expr": data[0],
            "stx_flds": data[1]
        }
        #toxin_wgs metadata
        self.data["metadata"]["TOXIN_WGS"] = "; ".join(data[1]["results"])
        
    def parseVirulenceFinder(self, parseFile, parseAddition):
        with open(parseFile, "r") as stxout:
            data = json.load(stxout)
        self.data["VirulenceFinder"] = data["qc"]

        with open(parseAddition, "r") as stxexpr:
            data = json.load(stxexpr)
        self.data["VirulenceFinder"]["genotypingResults"] = data

        #VirulenceMarker metadata
        ehx = data["results"]["ehxA"]
        eae = False
        for k,v in data["results"].items():
            if "eae-" in k and v:
                eae = True
                break
    
        if eae and ehx:
            self.data["metadata"]["VirulenceMarker"] = "eae+; ehx+"
        elif eae:
            self.data["metadata"]["VirulenceMarker"] = "eae+; ehx-"
        elif ehx:
            self.data["metadata"]["VirulenceMarker"] = "eae-; ehx+"
        else:
            self.data["metadata"]["VirulenceMarker"] = "eae-; ehx-"

    def parseCbotToxinTyping(self, parseFile, parseAddition):
        with open(parseFile, "r") as ctxn_out:
            data = json.load(ctxn_out)
        self.data["ToxinTyping"] = data["qc"]

        with open(parseAddition, "r") as ctxn:
            data = json.load(ctxn)
        self.data["ToxinTyping"]["genotypingResults"] = data
        self.data["metadata"]["ToxinSubtype"] = data["results"]["ToxinSubtype"]
    
    def parseVibrioSerotype(self, parseFile, parseAddition):
        with open(parseFile, "r") as stout:
            data = json.load(stout)
        self.data["Serotyper"] = data["qc"]

        with open(parseAddition, "r") as stjson:
            data = json.load(stjson)
        self.data["Serotyper"]["genotypingResults"] = {
            "Serotype_pcr": data[0]["results"],
            "Serogroup_wgs": data[1]["results"]["Serogroup_wgs"],
            "extra": data[0]["extra"]
        }
        #Only assign if cholerae
        if self.data["metadata"]["SPECIES"] == "cholerae":
            self.data["metadata"]["SEROTYPE_WGS"] = data[1]["results"]["Serogroup_wgs"]


    def parseVibrioVirulence(self, parseFile, parseAddition):
        with open(parseFile, "r") as virout:
            data = json.load(virout)
        self.data["Toxin_Typing"] = data["qc"]

        with open(parseAddition, "r") as virjson:
            data = json.load(virjson)
        self.data["Toxin_Typing"]["genotypingResults"] = {
            "Toxin_pcr": data[0]["results"],
            "Toxin_wgs": data[1]["results"]["Toxin_wgs"],
            "extra": data[0]["extra"]
        }
         #Only assign if cholerae
        if self.data["metadata"]["SPECIES"] == "cholerae":
            self.data["metadata"]["TOXIN_WGS"] = data[1]["results"]["Toxin_wgs"]

    def _add_fields(self, genus):
        if genus == "Salmonella":
            self._salmonella_fields()
        elif genus == "Yersinia":
            self._yersinia_fields()
        elif genus == "Campylobacter":
            self._campylobacter_fields()
        elif genus == "Escherichia":
            self._escherichia_fields()
        elif genus == "Vibrio":
            self._vibrio_fields()
        elif genus == "Clostridium":
            self._clostridium_fields()
        elif genus == "Listeria":
            self._listeria_fields()
        elif genus == "Cronobacter":
            self._cronobacter_fields()
        self._add_metadata()

    def _yersinia_fields(self):
        self.data.update(
            {
                "AlleleCalling": {},
                "AlleleFilter": {},
                "MLST": {},
            }
        )

    def _salmonella_fields(self):
        self.data.update(
            {
                "AlleleCalling": {},
                "AlleleFilter": {},
                "AlleleNaming(wgmlst)": {},
                "MLST": {},
                "AMRFinder(resistance)": {},
                "PlasmidFinder(plasmids)": {},
                "ORG_ANI": {},
                "isPCR": {},
                "SeqSero": {},
            }
        )

    def _campylobacter_fields(self):
        self.data.update(
            {
                "AlleleCalling": {},
                "AlleleFilter": {},
                "AlleleNaming(wgmlst)": {},
                "MLST": {},
                "AMRFinder(resistance)": {},
                "PlasmidFinder(plasmids)": {},
            }
        )

    def _escherichia_fields(self):
        self.data.update(
            {
                "AlleleCalling": {},
                "AlleleFilter": {},
                "AlleleNaming(wgmlst)": {},
                "MLST": {},
                "AMRFinder(resistance)": {},
                "PlasmidFinder(plasmids)": {},
                "isPCR": {},
                "KMA": {},
                "STXtyper": {},
                "STXcondenser": {},
                "PathotypeFinder": {},
                "VirulenceFinder": {},
            }
        )

    def _vibrio_fields(self):
        self.data.update(
            {
                "AlleleCalling": {},
                "AlleleFilter": {},
                "AlleleNaming(wgmlst)": {},
                "MLST": {},
                "AMRFinder(resistance)": {},
                "PlasmidFinder(plasmids)": {},
                "Toxin_Typing": {},
                "Serotyper": {},
            }
        )

    def _clostridium_fields(self):
        self.data.update(
            {
                "AlleleCalling": {},
                "AlleleFilter": {},
                "AlleleNaming(wgmlst)": {},
                "MLST": {},
                "AMRFinder(resistance)": {},
                "PlasmidFinder(plasmids)": {},
                "ToxinTyping": {},
            }
        )

    def _listeria_fields(self):
        self.data.update(
            {
                "AlleleCalling": {},
                "AlleleFilter": {},
                "AlleleNaming(wgmlst)": {},
                "MLST": {},
                "AMRFinder(resistance)": {},
                "PlasmidFinder(plasmids)": {}
            }
        )

    def _cronobacter_fields(self):
        self.data.update(
            {
                "AlleleCalling": {},
                "AlleleFilter": {},
                "AlleleNaming(wgmlst)": {},
                "MLST": {},
                "AMRFinder(resistance)": {},
                "PlasmidFinder(plasmids)": {}
            }
        )

class AssemblyMaster(BaseMaster):
    """Used for assembly pipeline entrypoint
    """
    def __init__(self):
        super().__init__()
        self.data = {
            "Summary": {
                "Genus": "",
                "Species": "",
                "AllelesCalledCoreCount": "",
                "AllelesCalledCorePercent": "",
                "AMRFinderGenes": "",
                "PlasmidsFound": "",
            },
            "ANI": {},
            "AlleleCalling": {},
            "metadata": {
                "GENUS": "",
                "SPECIES": "",
                "Subspecies": "",
                "ANI_Score": "",
                "MLST_ST": "",
                "MLST_CC": "",
            },
        }

class STECAssemblyTypingMaster(BaseMaster):
    """Used for STEC assembly-only typing pipeline entrypoint
    """
    def __init__(self):
        super().__init__()
        self.data = {
            "Summary": {
                "Genus": "",
                "Species": "",
                "AMRFinderGenes": "",
                "PlasmidsFound": "",
            },
            "ANI": {},
            "MLST": {},
            "AMRFinder(resistance)": {},
            "PlasmidFinder(plasmids)": {},
            "isPCR": {},
            "KMA": {},
            "STXtyper": {},
            "STXcondenser": {},
            "PathotypeFinder": {},
            "VirulenceFinder": {},
            "metadata": {
                "GENUS": "",
                "SPECIES": "",
                "Subspecies": "",
                "ANI_Score": "",
                "ANI_Percent_Aligned": "",
                "MLST_ST": "",
                "MLST_CC": "",
            },
        }

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--process",
        type=str,
        required=True,
        help="Name of the process. Ex: PrepareReads",
    )
    parser.add_argument(
        "--parseFile",
        type=str,
        required=False,
        help="Name of the process file to parse. Ex: outputs.json",
    )
    parser.add_argument(
        "--parseAddition",
        type=str,
        required=False,
        help="Name of the additional process file to parse. Ex: ani best_hit.json",
    )
    parser.add_argument(
        "--tempMaster",
        type=str,
        required=False,
        help="Name of the partially filled master file if created. Ex: PipelineProcessOutputs.json",
    )
    return parser.parse_args()


def main():

    # set of processes that don't require existing masteroutput file
    STEP_ONE_PROCESSES = {
        "genusidentify",
        "assemblyinit",
        "stecassemblytyping"
    }


    # set of processes that don't require a parsefile
    NO_PARSEFILE_PROCESSES = {
        "assemblyinit",
        "stecgroup",
        "stecassemblytyping"
    }

    args = parse_args()
    if args.process.lower() not in STEP_ONE_PROCESSES:  # load json if not step 1
        if args.tempMaster is None:
            sys.stderr.write(f"Error: {args.process} requires --tempMaster.\n")
            sys.exit(1)

    if args.process.lower() not in NO_PARSEFILE_PROCESSES:
        if args.parseFile is None:
            sys.stderr.write(f"Error: {args.process} requires --parseFile.\n")
            sys.exit(1)

    if args.tempMaster:
        data = BaseMaster.from_json(args.tempMaster)
    else:
        # Run alternative entrypoint inits if appropriate
        if args.process.lower() == "assemblyinit":
            data = AssemblyMaster()
            data.to_file()
            sys.exit(0)
        elif args.process.lower() == "stecassemblytyping":
            data = STECAssemblyTypingMaster()
            data.to_file()
            sys.exit(0)
        
        # otherwise normal main pipeline init
        data = BaseMaster()

    data.add_process(args.process, args.parseFile, args.parseAddition)
    data.to_file()


if __name__ == "__main__":
    main()
