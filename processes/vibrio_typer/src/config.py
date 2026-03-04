# All potential sequence based nucleotides
# and their parings
NUC_SIBLINGS = {
   'A': "A",
   'C': "C",
   'G': "G",
   'T': "T",
   'R': "AG",
   'Y': "CT",
   'S': "GC",
   'W': "AT",
   'K': "GT",
   'M': "AC",
   'B': "CGT",
   'D': "AGT",
   'H': "ACT",
   'V': "ACG",
   'N': "ACGT"
}

PCR_SETTINGS = {
    "SALM": {
        "PCR_MAX_MISMATCHES": 3,
        "PCR_MAX_NONIUPAC": 3,
        "PCR_PERCENT_IDENTITY": 0.7,
        "PCR_MAX_LENGTH_DEVIATION": 0.1,
    },
    "STEC": {
        "PCR_MAX_MISMATCHES": 3,
        "PCR_MAX_NONIUPAC": 3,
        "PCR_PERCENT_IDENTITY": 0.7,
        "PCR_MAX_LENGTH_DEVIATION": 0.1,
    },
    "VIBRIO.SEROTYPE": {
        "PCR_MAX_MISMATCHES": 1,
        "PCR_MAX_NONIUPAC": 2,
        "PCR_PERCENT_IDENTITY": 0.9,
        "PCR_MAX_LENGTH_DEVIATION": 0.1,
    },
    "VIBRIO.VIRULENCE": {
        "PCR_MAX_MISMATCHES": 1,
        "PCR_MAX_NONIUPAC": 2,
        "PCR_PERCENT_IDENTITY": 0.9,
        "PCR_MAX_LENGTH_DEVIATION": 0.1,
    }
}

PRIMERS = {
    "SALM": {
        "vi": [{
            "id": "vi-1",
            "locus": "vi",
            "forward_primer_id": "vi-F1",
            "forward_sequence": "AGGTTATTTCAGCATAAGGAGACTT",
            "reverse_primer_id": "vi-R2",
            "reverse_sequence": "CTCTTCCATACCACTTTCCGA",
            "expected_length": 443
        }],
        "typhi1": [{
            "id": "typhi-1",
            "locus": "typhi1",
            "forward_primer_id": "typhi1-F",
            "forward_sequence": "ATGAATACGAATAATTCACC",
            "reverse_primer_id": "typhi1-R",
            "reverse_sequence": "TTACCCTCCCCATGTCAC",
            "expected_length": 261
        }],
        "typhi2": [{
            "id": "typhi-2",
            "locus": "typhi2",
            "forward_primer_id": "typhi2-F",
            "forward_sequence": "ATGCCTGTTATGCATAATTG",
            "reverse_primer_id": "typhi2-R",
            "reverse_sequence": "TTATGCTGTTAACGAGTCGTC",
            "expected_length": 429
        }],
        "sdf": [{
            "id": "sdf-1",
            "locus": "sdf",
            "forward_primer_id": "sdf-F",
            "forward_sequence": "GTGGTGGCTGGCGAATGG",
            "reverse_primer_id": "sdf-R",
            "reverse_sequence": "GGAGAGGCGGTTTGATGTGG",
            "expected_length": 201
        }],
        "flhB1": [{
            "id": "flhB1-182",
            "locus": "flhB1-182",
            "forward_primer_id": "flhB1-F",
            "forward_sequence": "TTCGCGACGAATTTAAAGAGAGCGAAG",
            "reverse_primer_id": "flhB1-R",
            "reverse_sequence": "CAGCGTTTAAGCTGCCAGACCCAGGCC",
            "expected_length": 182
        },
        {
            "id": "flhB1-379",
            "locus": "flhB1-379",
            "forward_primer_id": "flhB1-F",
            "forward_sequence": "TTCGCGACGAATTTAAAGAGAGCGAAG",
            "reverse_primer_id": "flhB1-R",
            "reverse_sequence": "CAGCGTTTAAGCTGCCAGACCCAGGCC",
            "expected_length": 379
        }],
        "flhB2": [{
            "id": "flhB-2",
            "locus": "flhB2",
            "forward_primer_id": "flhB2-F",
            "forward_sequence": "GCAGCGCCGCATGATGGA",
            "reverse_primer_id": "flhB2-R",
            "reverse_sequence": "CAGCGTTTAAGCTGCCAGACCCAGGCC",
            "expected_length": 303
        }],
        "sseJ": [{
            "id": "sseJ",
            "locus": "sseJ",
            "forward_primer_id": "sseJ-F",
            "forward_sequence": "CACTAAAATCAGGAGTGGCT",
            "reverse_primer_id": "sseJ-R",
            "reverse_sequence": "TCCGTCATAAAACGCAAAAG",
            "expected_length": 1741
        }],
        "O22": [{
            "id": "O:22",
            "locus": "O:22",
            "forward_primer_id": "O22-F",
            "forward_sequence": "TAGAGAAAAAGCTATAAAAAAA",
            "reverse_primer_id": "O22-R",
            "reverse_sequence": "TCTACCAACAAACAAAATTTTATATTCCAAACACTCT",
            "expected_length": 127
        }],
        "O23": [{
            "id": "O:23",
            "locus": "O:23",
            "forward_primer_id": "O23-F",
            "forward_sequence": "TTTACTTAACGCTGGTTGTATAGA",
            "reverse_primer_id": "O23-R",
            "reverse_sequence": "AATTCATAAATCCCCTTTTCTCTGAGCAATCGGCCAA",
            "expected_length": 156
        }],
        "oafA1": [{
            "id": "oafA1",
            "locus": "oafA1",
            "forward_primer_id": "oafA1-F",
            "forward_sequence": "ACGAAGCACTTAGCAAGAACG",
            "reverse_primer_id": "oafA1-R",
            "reverse_sequence": "CAACAGCAACAACAATGAGGAC",
            "expected_length": 411
        }],
        "oafA2": [{
            "id": "oafA2",
            "locus": "oafA2",
            "forward_primer_id": "oafA2-F",
            "forward_sequence": "ACGAAGCACTTAGCAAGAACG",
            "reverse_primer_id": "oafA2-R",
            "reverse_sequence": "AATGACTAATAAAGGATATAAAATAT",
            "expected_length": 170
        }],
        "16Sa": [{
            "id": "16Sa",
            "locus": "16Sa",
            "forward_primer_id": "16Sa-F",
            "forward_sequence": "CAGGCCTAACACATGCAAGTC",
            "reverse_primer_id": "16Sa-R",
            "reverse_sequence": "GGGCGGTGTGTACAAGGC",
            "expected_length": 1362
        }],
        "16Sb": [{
            "id": "16Sb",
            "locus": "16Sb",
            "forward_primer_id": "16Sb-F",
            "forward_sequence": "GCCGCAAGGTTAAAACTCAA",
            "reverse_primer_id": "16Sb-R",
            "reverse_sequence": "AAGGCACCAATCCATCTCTG",
            "expected_length": 136
        }],
        "hilA": [{
            "id": "hilA",
            "locus": "hilA",
            "forward_primer_id": "hilA-F",
            "forward_sequence": "AGCGTATAGATAATAATCCGGGAT",
            "reverse_primer_id": "hilA-R",
            "reverse_sequence": "ATTCCACATTTTCTCGGCAATAG",
            "expected_length": 80
        }],
        "lacZ": [{
            "id": "lacZ",
            "locus": "lacZ",
            "forward_primer_id": "lacZ-F",
            "forward_sequence": "GCAAAACCTACCGGATTGAT",
            "reverse_primer_id": "lacZ-R",
            "reverse_sequence": "CTCCACCCTTTCATTCACCT",
            "expected_length": 131
        }],
        "ttr": [{
            "id": "ttr",
            "locus": "ttrA",
            "forward_primer_id": "ttr-F",
            "forward_sequence": "CTCACCAGGAGATTACAACATGG",
            "reverse_primer_id": "ttr-R",
            "reverse_sequence": "AGCTCAGACCAAAAGTGACCATC",
            "expected_length": 95
        }],
        "STM": [{
            "id": "SNP-",
            "locus": "STM3356",
            "forward_primer_id": "STM-F",
            "forward_sequence": "CACATTATTCGCTCAATGGAG",
            "reverse_primer_id": "STM-R",
            "reverse_sequence": "GTAAGGGTAATGGGTTCCAT",
            "expected_length": 290
        }]
    },
    "STEC": {
        "stx1a": [
            {
                "id": "stx1a-1",
                "locus": "stx1a",
                "forward_primer_id": "stx1a-F1",
                "forward_sequence": "CCTTTCCAGGTACAACAGCGGTT",
                "reverse_primer_id": "stx1a-R2",
                "reverse_sequence": "GGAAACTCATCAGATGCCATTCTGG",
                "expected_length": 478
            }
        ],
        "stx1c": [
            {
                "id": "stx1c-1",
                "locus": "stx1c",
                "forward_primer_id": "stx1c-F1",
                "forward_sequence": "CCTTTCCTGGTACAACTGCGGTT",
                "reverse_primer_id": "stx1c-R2",
                "reverse_sequence": "CAAGTGTTGTACGAAATCCCCTCTGA",
                "expected_length": 252
            }
        ],
        "stx1d": [
            {
                "id": "stx1d-1",
                "locus": "stx1d",
                "forward_primer_id": "stx1d-F1",
                "forward_sequence": "CAGTTAATGCGATTGCTAAGGAGTTTACC",
                "reverse_primer_id": "stx1d-R2",
                "reverse_sequence": "CTCTTCCTCTGGTTCTAACCCCATGATA",
                "expected_length": 203
            }
        ],
        "stx2a": [
            {
                "id": "stx2a-1",
                "locus": "stx2a",
                "forward_primer_id": "stx2a-F2",
                "forward_sequence": "GCGATACTGRGBACTGTGGCC",
                "reverse_primer_id": "stx2a-R3",
                "reverse_sequence": "CCGKCAACCTTCACTGTAAATGTG",
                "expected_length": 349
            },
            {
                "id": "stx2a-2",
                "locus": "stx2a",
                "forward_primer_id": "stx2a-F2",
                "forward_sequence": "GCGATACTGRGBACTGTGGCC",
                "reverse_primer_id": "stx2a-R2",
                "reverse_sequence": "GCCACCTTCACTGTGAATGTG",
                "expected_length": 347
            }
        ],
        "stx2b": [
            {
                "id": "stx2b-1",
                "locus": "stx2b",
                "forward_primer_id": "stx2b-F1",
                "forward_sequence": "AAATATGAAGAAGATATTTGTAGCGGC",
                "reverse_primer_id": "stx2b-R1",
                "reverse_sequence": "CAGCAAATCCTGAACCTGACG",
                "expected_length": 251
            }
        ],
        "stx2c": [
            {
                "id": "stx2c-1",
                "locus": "stx2c",
                "forward_primer_id": "stx2c-F1",
                "forward_sequence": "GAAAGTCACAGTTTTTATATACAACGGGTA",
                "reverse_primer_id": "stx2c-R2",
                "reverse_sequence": "CCGGCCACYTTTACTGTGAATGTA",
                "expected_length": 177
            }
        ],
        "stx2d": [
            {
                "id": "stx2d-1",
                "locus": "stx2d",
                "forward_primer_id": "stx2d-F1",
                "forward_sequence": "AAARTCACAGTCTTTATATACAACGGGTG",
                "reverse_primer_id": "stx2d-R1",
                "reverse_sequence": "TTYCCGGCCACTTTTACTGTG",
                "expected_length": 179
            },
            {
                "id": "stx2d-2",
                "locus": "stx2d",
                "forward_primer_id": "stx2d-F1",
                "forward_sequence": "AAARTCACAGTCTTTATATACAACGGGTG",
                "reverse_primer_id": "stx2d-O55-R",
                "reverse_sequence": "TCAACCGAGCACTTTGCAGTAG",
                "expected_length": 235
            },
            {
                "id": "stx2d-3",
                "locus": "stx2d",
                "forward_primer_id": "stx2d-F1",
                "forward_sequence": "AAARTCACAGTCTTTATATACAACGGGTG",
                "reverse_primer_id": "stx2d-R2",
                "reverse_sequence": "GCCTGATGCACAGGTACTGGAC",
                "expected_length": 280
            }
        ],
        "stx2e": [
            {
                "id": "stx2e-1",
                "locus": "stx2e",
                "forward_primer_id": "stx2e-F1",
                "forward_sequence": "CGGAGTATCGGGGAGAGGC",
                "reverse_primer_id": "stx2e-R2",
                "reverse_sequence": "CTTCCTGACACCTTCACAGTAAAGGT",
                "expected_length": 411
            }
        ],
        "stx2f": [
            {
                "id": "stx2f-1",
                "locus": "stx2f",
                "forward_primer_id": "stx2f-F1",
                "forward_sequence": "TGGGCGTCATTCACTGGTTG",
                "reverse_primer_id": "stx2f-R1",
                "reverse_sequence": "TAATGGCCGCCCTGTCTCC",
                "expected_length": 424
            }
        ],
        "stx2g": [
            {
                "id": "stx2g-1",
                "locus": "stx2g",
                "forward_primer_id": "stx2g-F1",
                "forward_sequence": "CACCGGGTAGTTATATTTCTGTGGATATC",
                "reverse_primer_id": "stx2g-R1",
                "reverse_sequence": "GATGGCAATTCAGAATAACCGCT",
                "expected_length": 573
            }
        ],
        "eae": [
            {
                "id": "eae-1",
                "locus": "eae",
                "forward_primer_id": "eae_776F",
                "forward_sequence": "GGCATTTGGTCAGGTCGG",
                "reverse_primer_id": "eae_907R",
                "reverse_sequence": "TCACCAGAAAAATCCTGATCAA",
                "expected_length": 132
            }
        ],
        "ehxA": [
            {
                "id": "ehxA-1",
                "locus": "ehxA",
                "forward_primer_id": "ehxAF",
                "forward_sequence": "GCATCATCAAGCGTACGTTCC",
                "reverse_primer_id": "ehxAR",
                "reverse_sequence": "AATGAGCCAAGCTGGTTAAGCT",
                "expected_length": 534
            }
        ],
        "ipaH": [
            {
                "id": "ipaH",
                "locus": "ipaH",
                "forward_primer_id": "EIEC-F",
                "forward_sequence": "GTTCCTTGACCGCCTTTCCGA",
                "reverse_primer_id": "EIEC-R",
                "reverse_sequence": "GCCGGTCAGCCACCCTCTGA",
                "expected_length": 620
            }
        ]
    },
    "VIBRIO.SEROTYPE": {
        "wbeN-587": [
            {
                "id": "wbeN-O1",
                "locus": "wbeN-587",
                "forward_primer_id": "wbeN-587F",
                "forward_sequence": "ACATCTGAATTCACTTGCGAGT",
                "reverse_primer_id": "wbeN-694R",
                "reverse_sequence": "ACCTGAGAACCTTTAGGCAATT",
                "expected_length": 129
            },
        ],
        "wbfR": [
            {
                "id": "wbfR-O139",
                "locus": "wbfR",
                "forward_primer_id": "wbfR-540F",
                "forward_sequence": "GGCGTTATCGCATTTTTTTCGTTT",
                "reverse_primer_id": "wbfR-664R",
                "reverse_sequence": "GACTGGCATCCCAAAATGTTTG",
                "expected_length": 146
            }
        ]

    },
    "VIBRIO.VIRULENCE": {
        "ctxA-2": [
            {
                "id": "ctxA-2",
                "locus": "ctxA-2",
                "forward_primer_id": "ctxA-2F",
                "forward_sequence": "CGGGCAGATTCTAGACCTCCTG",
                "reverse_primer_id": "ctxA-2R",
                "reverse_sequence": "CGATGATCTTGGAGCATTCCCAC",
                "expected_length": 564
            }
        ]
    }
}
