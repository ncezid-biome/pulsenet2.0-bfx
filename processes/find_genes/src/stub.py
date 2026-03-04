def get_stub_content():
    """
    Return stubbed find-genes content
    """
    header = [
        "Protein identifier",
        "Contig id",
        "Start",
        "Stop",
        "Strand",
        "Gene symbol",
        "Sequence name",
        "Scope",
        "Element type",
        "Element subtype",
        "Class",
        "Subclass",
        "Method",
        "Target length",
        "Reference sequence length",
        "% Coverage of reference sequence",
        "% Identity to reference sequence",
        "Alignment length",
        "Accession of closest sequence",
        "Name of closest sequence",
        "HMM id",
        "HMM description",
    ]
    data = [
        [
            "NA",
            "NODE_6_length_154682_cov_25.248591_pilon",
            "120296",
            "121066",
            "-",
            "blaOXA-193",
            "OXA-61 family class D beta-lactamase OXA-193",
            "core",
            "AMR",
            "AMR",
            "BETA-LACTAM",
            "BETA-LACTAM",
            "ALLELEX",
            "257",
            "257",
            "100.00",
            "100.00",
            "257",
            "WP_002783228.1",
            "OXA-61 family class D beta-lactamase OXA-193",
            "NA",
            "NA",
        ]
    ]
    data.insert(0, header)
    return "\n".join("\t".join(row) for row in data)
