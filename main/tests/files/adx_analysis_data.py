import csv
import json
import re
import argparse
import gzip
       
#input arguments
parser = argparse.ArgumentParser()
parser.add_argument("--masterOutFile", type=str, required=True, help="MasterOutput.json")
parser.add_argument("--genotypeResults", type=str, required=False, help="GenotypeResults.json")
parser.add_argument("--fastpFile", type=str, required=True, help="fastp output from Prepare Reads")
parser.add_argument("--wgmlst", type=str, required=False, help="calls_standard.json.gz from Allele Filtering")
args = parser.parse_args()

## Load inputs into variables
fastp = json.load(open(args.fastpFile))
masterOut = json.load(open(args.masterOutFile))
all_outputs = []
values = {}

## Prepare Reads QC
values["SrsQ30"] = fastp["summary"]["before_filtering"]["q30_bases"]
values["SrsQ30Freq"] = fastp["summary"]["before_filtering"]["q30_rate"]
values["SrsQ30_1"] = fastp["read1_before_filtering"]["q30_bases"]
values["SrsQ30_2"] = fastp["read2_before_filtering"]["q30_bases"]
values["SrsQ30Freq_1"] = float(fastp["read1_before_filtering"]["q30_bases"])/float(fastp["read1_before_filtering"]["total_bases"])
values["SrsQ30Freq_2"] = float(fastp["read2_before_filtering"]["q30_bases"])/float(fastp["read2_before_filtering"]["total_bases"])
values["AvgReadCoverage"] = masterOut["PrepareReads"]["metrics"]["initialEstimatedVerticalDepth"]
values["secondaryCoverage_coverage"] = masterOut["MIDAS"]["metrics"]["secondaryCoverage_coverage"]
values["secondaryCoverage_count"] = masterOut["MIDAS"]["metrics"]["secondaryCoverage_count"]

AvgQuality_read1 = sum(list(fastp["read1_before_filtering"]["quality_curves"]["mean"]))/len(list(fastp["read1_before_filtering"]["quality_curves"]["mean"]))
AvgQuality_read2 = sum(list(fastp["read2_before_filtering"]["quality_curves"]["mean"]))/len(list(fastp["read2_before_filtering"]["quality_curves"]["mean"]))
values["AvgQuality"] = (AvgQuality_read1 + AvgQuality_read2)/2

values["rawRead1MeanLength"] = masterOut["PrepareReads"]["metrics"]["rawRead1MeanLength"] 
values["rawRead2MeanLength"] = masterOut["PrepareReads"]["metrics"]["rawRead2MeanLength"]
values["rawCombinedReadMeanLength"] = masterOut["PrepareReads"]["metrics"]["rawCombinedReadMeanLength"]

values["SrsQ30_Trimmed"] = fastp["summary"]["after_filtering"]["q30_bases"]
values["SrsQ30Freq_Trimmed"] = fastp["summary"]["after_filtering"]["q30_rate"]
values["SrsQ30_1_Trimmed"] = fastp["read1_after_filtering"]["q30_bases"]
values["SrsQ30_2_Trimmed"] = fastp["read2_after_filtering"]["q30_bases"]
values["SrsQ30Freq_1_Trimmed"] = float(fastp["read1_after_filtering"]["q30_bases"])/float(fastp["read1_after_filtering"]["total_bases"])
values["SrsQ30Freq_2_Trimmed"] = float(fastp["read2_after_filtering"]["q30_bases"])/float(fastp["read2_after_filtering"]["total_bases"])
values["AvgReadCoverageTrimmed"] = masterOut["PrepareReads"]["metrics"]["postSubsamplingEstimatedVerticalDepth"]

AvgQualityTrimmed_read1 = sum(list(fastp["read1_after_filtering"]["quality_curves"]["mean"]))/len(list(fastp["read1_after_filtering"]["quality_curves"]["mean"]))
AvgQualityTrimmed_read2 = sum(list(fastp["read2_after_filtering"]["quality_curves"]["mean"]))/len(list(fastp["read2_after_filtering"]["quality_curves"]["mean"]))
values["AvgQualityTrimmed"] = (AvgQualityTrimmed_read1 + AvgQualityTrimmed_read2)/2

values["trimmedRead1MeanLength"] = masterOut["PrepareReads"]["metrics"]["trimmedRead1MeanLength"]
values["trimmedRead2MeanLength"] = masterOut["PrepareReads"]["metrics"]["trimmedRead2MeanLength"]
values["trimmedCombinedReadMeanLength"] = masterOut["PrepareReads"]["metrics"]["trimmedCombinedReadMeanLength"]

## Assembly QC
values["NrNonACGT"] = ""
values["NrContigs"] = masterOut["CorrectAssembly"]["metrics"]["cleanedNContigs"]
values["AvgDeNovoCover"] = masterOut["CorrectAssembly"]["metrics"]["cleanedAverageDepth"]
values["NrBasesN"] = ""
values["NrBasesACGT"] = ""
values["N50"] = masterOut["CorrectAssembly"]["metrics"]["cleanedN50"]
values["Length"] = masterOut["CorrectAssembly"]["metrics"]["cleanedLength"]

# ## Allele Calling QC
if args.wgmlst is not None:

    values["NrBAFPerfect"] = ""
    values["AvgLocusCover"] = ""
    values["NrBAFMultiple"] = ""
    values["NrBAFPresent"] = masterOut["AlleleFilter"]["metrics"]["coreCount"] + masterOut["AlleleFilter"]["metrics"]["accessoryCount"]
    values["NrToBeSubmitted"] = ""
    values["NrAlreadySubmitted"] = ""
    values["CorePercentCalled"] = masterOut["AlleleFilter"]["metrics"]["corePercentage"]


# This is for organisms that do not run allele calling (i.e. Yersinia)
else:
    values["NrBAFPerfect"] = ""
    values["AvgLocusCover"] = ""
    values["NrBAFMultiple"] = ""
    values["NrBAFPresent"] = ""
    values["NrToBeSubmitted"] = ""
    values["NrAlreadySubmitted"] = ""
    values["CorePercentCalled"] = ""

#Assembly Free based metrics below
values["NrConsensusMultiple"] = "" #AF & AB
values["NrConsensusUnknown"] = "" #AF & AB
values["CorePercent"] = "" #AF & AB
values["NrDifferent"] = "" #AF & AB
values["NrAFPerfect"] = "" #AF
values["KeywordCov"] = "" #AF
values["NrAFMultiple"] = "" #AF
values["NrAFPresent"] = "" #AF
values["NrConsensus"] = "" #AF & AB
values["NrConsensusConfirmed"] = "" #AF & AB

## Combine All QC
qc_output = {}
qc_output['processType'] = {"name" : "quality"}
qc_output['result'] = {"type": "characterization", "values": values}
all_outputs.append(qc_output)

## Allele IDs from filtering
if args.wgmlst is not None:
    wgmlst_data = json.load(gzip.open(args.wgmlst))
    allele_out = {}
    allele_out['processType'] = {"name" : "wgmlst"}
    allele_out['result'] = {"type": "characterization", "values": wgmlst_data["values"] }
    all_outputs.append(allele_out)

#Check for Genotyping Results
if args.genotypeResults is not None:
    genotype_input = json.load(open(args.genotypeResults))
    genotype_content = open(args.genotypeResults).read()

    ## Plasmid Results
    plasmid_out = {}
    plasmid_out['processType'] = {"name" : "plasmids"}
    oldVal = genotype_input['plasmids.json']['results']
    newVal = { key: 1 if value else 2 for key, value in oldVal.items() }
    plasmid_out['result'] = {"type": "characterization", "values": newVal}
    all_outputs.append(plasmid_out)

    ## Resistance (AMRFinder/Find Genes) Results
    resist_out = {}
    resist_out['processType'] = {"name" : "resistance"}
    oldVal = genotype_input['resistance.json']['results']
    newVal = { key: 1 if value else 2 for key, value in oldVal.items() }
    resist_out['result'] = {"type": "characterization", "values": newVal}
    all_outputs.append(resist_out)

    ## GenotypingResults Full File
    jsonString = re.sub(r'\s+', ' ', genotype_content).replace('{ ', '{')
    geno_out = {}
    geno_out['processType'] = {"name" : "genotypingResult"}
    geno_out['result'] = {"type": "genotypingResult", "rawJsonAnalysis": jsonString}
    all_outputs.append(geno_out)

    ## STEC Virulence Results
    if genotype_input.get("virulence.json") and genotype_input.get("ecoli.serotype.json"):
        virulence_out = {}
        virulence_out['processType'] = {"name" : "virulence_finder"}
        oldVal = genotype_input['virulence.json']['results']
        newVal = { key: 1 if value else 2 for key, value in oldVal.items() }
        virulence_out['result'] = {"type": "characterization", "values": newVal}
        all_outputs.append(virulence_out)

##### Print Results #####
with open('ADX_analysis_data.json', 'w') as output_file:
    json.dump(all_outputs, output_file, indent=4)
