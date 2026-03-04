PLASMIDFINDER_SETTINGS = {
    "SALM": {
        "PF_DATABASE": "",  
        "PF_PERCENT_IDENTITY": 0.1,
        "PF_MIN_RELATIVE_COVERAGE": 0.1,
        "PF_MIN_MERGE_OVERLAP": 0.9,
        "PF_SEARCH_FRAGMENTS": False,
        "VERSION": {
            "ALGORITHM": "", 
            "DATABASE": "", 
        }
    }
}

PCR_SETTINGS = {
    "SALM": {
        "PCR_MAX_MISMATCHES": 3,
        "PCR_MAX_NONIUPAC": 3,
        "PCR_PERCENT_IDENTITY": 0.7,
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
    }
}

ANI_SETTINGS = {
    "ALL": {
        "PERCENT_IDENTITY": 80.0,
        "MIN_COVERAGE": 70.0,
        "DISCRIMINATION": 2.0,
    }
}

CONTINGENCY_TABLES = {
    "SALM": {
        "I 9:d:-": [
            {
                "vi": True,
                "typhi1": True,
                "typhi2": True,
                "BN Serotype": "Typhi"
            },
            {
                "vi": False,
                "typhi1": True,
                "typhi2": True,
                "BN Serotype": "Typhi"
            }
        ],
        "I 9:j:-": [
            {
                "vi": True,
                "typhi1": True,
                "typhi2": True,
                "BN Serotype": "Typhi"
            }
        ],
        "I 9:g,m:-": [
            {
                "sdf": True,
                "flhB1-182": False,
                "flhB1-379": True,
                "flhB2": True,
                "BN Serotype": "Enteritidis"
            },
            {
                "sdf": False,
                "flhB1-182": False,
                "flhB1-379": True,
                "flhB2": True,
                "BN Serotype": "Enteritidis"
            },
            {
                "sdf": False,
                "flhB1-182": True,
                "flhB1-379": False,
                "flhB2": False,
                "BN Serotype": "Gallinarum"
            }
        ],
        "I 4:b:1,2": [
            {
                "STM3356": True,
                "sseJ": True,
                "BN Serotype": "Paratyphi B var. L(+) tartrate+"
            },
            {
                "STM3356": False,
                "sseJ": False,
                "BN Serotype": "Paratyphi B"
            }
        ],
        "I 4:b:-": [
            {
                "STM3356": True,
                "sseJ": True,
                "BN Serotype": "I 4:b:-"
            },
            {
                "STM3356": False,
                "sseJ": False,
                "BN Serotype": "Paratyphi B"
            }
        ],
        "I 13:b:1,5": [
            {
                "O:22": False,
                "O:23": True,
                "BN Serotype": "Mississippi"
            },
            {
                "O:22": True,
                "O:23": False,
                "BN Serotype": "Ibadan"
            }
        ],
        "I 13:z:1,6": [
            {
                "O:22": True,
                "O:23": False,
                "BN Serotype": "Poona"
            },
            {
                "O:22": False,
                "O:23": True,
                "BN Serotype": "Farmsen"
            }
        ],
        "I 13:z29:-": [
            {
                "O:22": False,
                "O:23": True,
                "BN Serotype": "Cubana"
            },
            {
                "O:22": True,
                "O:23": False,
                "BN Serotype": "Agoueve"
            }
        ],
        "I 13:g,m:-": [
            {
                "O:22": False,
                "O:23": True,
                "BN Serotype": "Agbeni"
            }
        ],
        "I 13:d:1,7": [
            {
                "O:22": False,
                "O:23": True,
                "BN Serotype": "Grumpensis"
            }
        ],
        "I 13:f,g:-": [
            {
                "O:22": False,
                "O:23": True,
                "BN Serotype": "Havana"
            }
        ],
        "I 13:m,t:-": [
            {
                "O:22": False,
                "O:23": True,
                "BN Serotype": "Kintambo"
            }
        ],
        "I 13:d:l,w": [
            {
                "O:22": False,
                "O:23": True,
                "BN Serotype": "Putten"
            }
        ],
        "I 13:d:e,n,z15": [
            {
                "O:22": False,
                "O:23": True,
                "BN Serotype": "Telelkebir"
            }
        ],
        "I 13:z:l,w": [
            {
                "O:22": False,
                "O:23": True,
                "BN Serotype": "Worthington"
            }
        ],
        "Group O4": [
            {
                "oafA1": True,
                "oafA2": False,
                "BN Serotype": "Needs Further Review"
            },
            {
                "oafA1": True,
                "oafA2": True,
                "BN Serotype": "Needs Further Review"
            },
            {
                "oafA1": False,
                "oafA2": False,
                "BN Serotype": "Needs Further Review"
            }
        ]
    }
}

LOOKUP_TABLES = {
    "SALM": {
        "I": {
            "11:i:1,2": "Aberdeen",
            "4:b:e,n,x": "Abony",
            "35:f,g:-": "Adelaide",
            "30:z38:-": "Ago",
            "4:f,g,s:-": "Agona",
            "35:z4,z23:-": "Alachua",
            "8:z4,z24:-": "Albany",
            "40:k:1,6": "Allandale",
            "8:r:z6": "Altona",
            "3,10:y:1,2": "Amager",
            "3,10:g,m,s:-": "Amsterdam",
            "3,10:e,h:1,6": "Anatum",
            "35:g,s,t:-": "Anecho",
            "45:m,t:-": "Apapa",
            "30:k:1,6": "Aqua",
            "4:a:1,7": "Arechavaleta",
            "7:i:1,2": "Augustenborg",
            "9,46:a:e,n,x": "Baildon",
            "4:y:e,n,x": "Ball",
            "7:y:1,5": "Bareilly",
            "16:d:e,n,x": "Barranquilla",
            "47:i:e,n,z15": "Bergen",
            "43:a:1,5": "Berkeley",
            "9:f,g,t:-": "Berta",
            "7:c:1,6": "Birkenhead",
            "8:k:1,5": "Blockley",
            "4:r,i:l,w": "Bochum",
            "8:i:e,n,x": "Bonariensis",
            "7:l,v:e,n,x": "Bonn",
            "8:r:1,5": "Bovismorbificans",
            "7:e,h:e,n,z15": "Braenderup",
            "4:z29:-": "Brancaster",
            "4:l,v:e,n,z15": "Brandenburg",
            "4:l,v:1,7": "Bredeney",
            "11:l,v:e,n,z15": "Bullbay",
            "1,3,19:m,t:-": "Cannstat",
            "17:l,v:e,n,x": "Carmel",
            "6,14:y:1,7": "Carrau",
            "18:z4,z23:-": "Cerro",
            "8:z4,z23:e,n,z15": "Chailey",
            "4:e,h:e,n,x": "Chester",
            "11:e,h:1,2": "Chingola",
            "4:l,v:1,6": "Clackamas",
            "4:y:1,2": "Coeln",
            "7:r:1,7": "Colindale",
            "7:l,v:1,2": "Concord",
            "8:z4,z23:-": "Corvallis",
            "28:i:1,5": "Cotham",
            "7:k:1,6": "Daytona",
            "4:f,g:-": "Derby",
            "9:g,p:-": "Dublin",
            "4:d:e,n,z15": "Duisburg",
            "9:a:e,n,z15": "Durban",
            "35:g,m,s:-": "Ealing",
            "9:e,h:1,5": "Eastbourne",
            "7:b:1,5": "Edinburg",
            "4:g,m:-": "Essen",
            "18:b:1,5": "Fluntern",
            "16:d:1,7": "Gaminara",
            "3,10:l,v:1,7": "Give",
            "44:z10:1,7": "Guinea",
            "8:z10:e,n,x": "Hadar",
            "4:z10:1,2": "Haifa",
            "7:y:e,n,x": "Hartford",
            "4:g,m,s:-": "Hato",
            "4:r:1,2": "Heidelberg",
            "8:l,v:e,n,x": "Holcomb",
            "16:b:e,n,x": "Hvittingfoss",
            "4:i:-": "I 4:i:-",
            "9:l,z28:-": "I 9:l,z28:-",
            "4:z:1,7": "Indiana",
            "7:r:1,5": "Infantis",
            "38:k:1,6": "Inverness",
            "7:l,v:1,5": "Irumu",
            "7:d:1,5": "Isangi",
            "4:z10:1,5": "Ituri",
            "9:l,z28:1,5": "Javiana",
            "40:b:e,n,x": "Johannesburg",
            "8:i:z6": "Kentucky",
            "4:z:1,5": "Kiambu",
            "43:y:1,5": "Kingabwa",
            "39:l,v:e,n,x": "Kokomlemle",
            "8:e,h:1,5": "Kottbus",
            "30:z10:e,n,z15": "Kumasi",
            "38:i:1,5": "Lansing",
            "7:e,h:1,2": "Larochelle",
            "45:z35:1,5": "Lattenkamp",
            "7:i:1,7": "Lika",
            "7:z38:-": "Lille",
            "8:l,v:1,2": "Litchfield",
            "1,3,19:d:e,n,z15": "Liverpool",
            "7:d:l,w": "Livingstone",
            "9:a:e,n,x": "Lomalinda",
            "3,10:l,v:1,6": "London",
            "8:d:1,5": "Manhattan",
            "17:k:e,n,x": "Matadi",
            "30:y:1,2": "Matopeni",
            "7:z10:e,n,z15": "Mbandaka",
            "3,10:e,h:l,w": "Meleagridis",
            "38:i:1,2": "Mgulani",
            "17:l,v:1,5": "Michigan",
            "7:y:e,n,z15": "Mikawasima",
            "21:b:e,n,x": "Minnesota",
            "35:m,t:-": "Monschaui",
            "7:g,m,s:-": "Montevideo",
            "47:z:1,5": "Mountpleasant",
            "8:d:1,2": "Muenchen",
            "3,10:e,h:1,5": "Muenster",
            "9:l,z13:e,n,x": "Napoli",
            "8:e,h:1,2": "Newport",
            "28:y:1,5": "Nima",
            "7:e,h:1,6": "Norwich",
            "16:d:e,n,z15": "Nottingham",
            "7:b:l,w": "Ohio",
            "7:m,t:-": "Oranienburg",
            "16:k:e,n,z15": "Orientalis",
            "7:i:1,5": "Oritamerin",
            "7:a:e,n,x": "Oslo",
            "9:l,v:1,5": "Panama",
            "2:a:1,5": "Paratyphi A",
            "9:m,t:-": "Pensacola",
            "28:y:1,7": "Pomona",
            "7:l,v:e,n,z15": "Potsdam",
            "4:e,h:1,5": "Reading",
            "7:y:1,2": "Richmond",
            "7:f,g:-": "Rissen",
            "11:r:e,n,x": "Rubislaw",
            "4:e,h:1,2": "Saintpaul",
            "4:e,h:e,n,z15": "Sandiego",
            "16:y:1,5": "Saphra",
            "4:d:1,7": "Schwarzengrund",
            "1,3,19:g,s,t:-": "Senftenberg",
            "4:z:1,2": "Shubra",
            "7:k:e,n,x": "Singapore",
            "30:i:l,w": "Soerenga",
            "4:d:1,2": "Stanley",
            "4:z4,z23:-": "Stanleyville",
            "6,14:z:e,n,x": "Sundsvall",
            "8:z4,z32:-": "Tallahassee",
            "6,14:d:e,n,z15": "Teko",
            "7:z29:-": "Tennessee",
            "7:k:1,5": "Thompson",
            "13:g,m:-": "to genotyper",
            "13:d:1,7": "to genotyper",
            "13:f,g:-": "to genotyper",
            "13:m,t:-": "to genotyper",
            "13:d:l,w": "to genotyper",
            "13:d:e,n,z15": "to genotyper",
            "9:j:-": "to genotyper",
            "13:b:1,5": "to genotyper",
            "13:z:1,6": "to genotyper",
            "13:z29:-": "to genotyper",
            "4:b:-": "to genotyper",
            "4:b:1,2": "to genotyper",
            "9:a:1,5": "to genotyper",
            "9:d:-": "to genotyper",
            "9:g,m:-": "to genotyper",
            "13:z:l,w": "to genotyper",
            "48:z:1,5": "Toucra",
            "4:i:1,2": "Typhimurium",
            "3,10:l,z13:1,5": "Uganda",
            "30:b:e,n,x": "Urbana",
            "7:r:1,2": "Virchow",
            "28:l,v:e,n,x": "Vitkin",
            "4:b:e,n,z15": "Wagenia",
            "39:b:1,2": "Wandsworth",
            "3,10:r:z6": "Weltevreden",
            "35:z29:-": "Widemarsh",
            "4:b:l,w": "Wien"
        },
        "IIIa": {
            "18:z4,z23:-": "IIIa 18:z4,z23:-",
            "41:z4,z23:-": "IIIa 41:z4,z23:-",
            "48:g,z51:-": "IIIa 48:g,z51:-"
        },
        "IIIb": {
            "48:i:z": "IIIb 48:i:z"
        },
        "IV": {
            "44:z4,z23:-": "IV 44:z4,z23:-",
            "45:g,z51:-": "IV 45:g,z51:-",
            "48:g,z51:-": "IV 48:g,z51:-",
            "50:g,z51:-": "IV 50:g,z51:-",
            "50:z4,z23:-": "IV 50:z4,z23:-"
        }
    }
}