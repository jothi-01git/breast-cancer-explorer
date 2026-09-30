from flask import Flask, render_template, request
from concurrent.futures import ThreadPoolExecutor
import requests


app = Flask(__name__)


# --------------------------------------------------
# Open Targets: Search for the disease
# --------------------------------------------------

def search_disease(disease_name):

    url = "https://api.platform.opentargets.org/api/v4/graphql"

    query = """
    query {
        search(queryString: "%s") {
            hits {
                id
                name
                entity
            }
        }
    }
    """ % disease_name

    response = requests.post(
        url,
        json={"query": query}
    )

    data = response.json()

    for result in data["data"]["search"]["hits"]:

        if (
            result["entity"] == "disease"
            and result["name"].lower() == disease_name.lower()
        ):
            return result["id"]

    return None


# --------------------------------------------------
# NCBI MedGen: Get disease description
# --------------------------------------------------

def get_disease_description(disease_id):

    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"

    # Convert MONDO_0007254 to MONDO:0007254
    mondo_id = disease_id.replace("_", ":")

    params = {
        "db": "medgen",
        "term": '"' + mondo_id + '"[sourceid]',
        "retmode": "json",
        "retmax": 5
    }

    response = requests.get(
        url,
        params=params,
        timeout=30
    )

    print(
        "MedGen search status:",
        response.status_code
    )

    if response.status_code != 200:

        print("NCBI MedGen search error:")
        print(response.text)

        return ""

    data = response.json()

    ids = data["esearchresult"]["idlist"]

    print(
        "MedGen IDs:",
        ids
    )

    if not ids:

        print("No MedGen record found.")

        return ""

    medgen_id = ids[0]

    # ----------------------------------------------
    # Get MedGen record
    # ----------------------------------------------

    summary_url = (
        "https://eutils.ncbi.nlm.nih.gov/"
        "entrez/eutils/esummary.fcgi"
    )

    summary_params = {
        "db": "medgen",
        "id": medgen_id,
        "retmode": "json"
    }

    summary_response = requests.get(
        summary_url,
        params=summary_params,
        timeout=30
    )

    print(
        "MedGen summary status:",
        summary_response.status_code
    )

    if summary_response.status_code != 200:

        print("NCBI MedGen summary error:")
        print(summary_response.text)

        return ""

    summary_data = summary_response.json()

    record = summary_data["result"].get(
        medgen_id,
        {}
    )

    print("MedGen record:")
    print(record)

    return record.get(
        "definition",
        {}
    ).get(
        "value",
        ""
    )


# --------------------------------------------------
# Open Targets: Get disease-associated genes
# --------------------------------------------------

def get_disease_genes(disease_id):

    url = "https://api.platform.opentargets.org/api/v4/graphql"

    query = """
    query {
        disease(efoId: "%s") {

            associatedTargets(
                page: {index: 0, size: 10}
            ) {

                rows {

                    target {
                        id
                        approvedSymbol
                        approvedName
                    }

                    score
                }
            }
        }
    }
    """ % disease_id

    response = requests.post(
        url,
        json={"query": query}
    )

    data = response.json()

    return data[
        "data"
    ][
        "disease"
    ][
        "associatedTargets"
    ][
        "rows"
    ]


# --------------------------------------------------
# Ensembl: Get genomic locations
# --------------------------------------------------

def get_gene_locations(gene_ids):

    url = "https://rest.ensembl.org/lookup/id"

    params = {
        "content-type": "application/json"
    }

    for attempt in range(3):

        response = requests.post(
            url,
            params=params,
            json={
                "ids": gene_ids
            },
            timeout=30
        )

        print(
            "Ensembl attempt:",
            attempt + 1
        )

        print(
            "Ensembl status code:",
            response.status_code
        )

        if response.status_code == 200:

            data = response.json()

            print(
                "Ensembl genes received:",
                len(data)
            )

            return data

    print(
        "Ensembl failed after 3 attempts."
    )

    print(
        response.text
    )

    return {}


# --------------------------------------------------
# Ensembl: Get gene neighborhood
# --------------------------------------------------

def get_gene_neighborhood(
    chromosome,
    gene_start,
    gene_end
):

    server = "https://rest.ensembl.org"

    # 500 kb upstream and downstream
    region_start = max(
        1,
        gene_start - 500000
    )

    region_end = gene_end + 500000

    region = (
        f"{chromosome}:"
        f"{region_start}-"
        f"{region_end}"
    )

    url = (
        f"{server}/overlap/region/human/"
        f"{region}"
    )

    params = {
        "feature": "gene"
    }

    headers = {
        "Content-Type": "application/json"
    }

    try:

        response = requests.get(
            url,
            params=params,
            headers=headers,
            timeout=30
        )

        if response.status_code == 200:

            nearby_genes = response.json()

            neighborhood = []

            for gene in nearby_genes:

                neighborhood.append({

                    "name":
                        gene.get(
                            "external_name"
                        )
                        or gene.get(
                            "id",
                            ""
                        ),

                    "ensembl_id":
                        gene.get(
                            "id",
                            ""
                        ),

                    "start":
                        gene.get(
                            "start",
                            ""
                        ),

                    "end":
                        gene.get(
                            "end",
                            ""
                        ),

                    "strand":
                        gene.get(
                            "strand",
                            ""
                        )

                })

            return neighborhood

        else:

            print(
                "Gene neighborhood error "
                "for region:",
                region
            )

            print(
                "Status code:",
                response.status_code
            )

            print(
                response.text
            )

            return []

    except requests.exceptions.RequestException as e:

        print(
            "Gene neighborhood exception "
            "for region:",
            region
        )

        print(e)

        return []


# --------------------------------------------------
# Gene neighborhood: Process one gene
# --------------------------------------------------

def process_gene_neighborhood(
    item,
    locations
):

    target = item["target"]

    gene_symbol = target[
        "approvedSymbol"
    ]

    gene_id = target[
        "id"
    ]

    location = locations.get(
        gene_id,
        {}
    )

    chromosome = location.get(
        "seq_region_name"
    )

    gene_start = location.get(
        "start"
    )

    gene_end = location.get(
        "end"
    )

    if not (
        chromosome
        and gene_start
        and gene_end
    ):

        print(
            "No genomic location for",
            gene_symbol
        )

        return (
            gene_symbol,
            []
        )

    print(
        "Getting neighborhood for",
        gene_symbol
    )

    neighborhood = get_gene_neighborhood(
        chromosome,
        gene_start,
        gene_end
    )

    print(
        "Neighborhood genes/features for",
        gene_symbol,
        ":",
        len(neighborhood)
    )

    return (
        gene_symbol,
        neighborhood
    )


# --------------------------------------------------
# PubMed: Search articles
# --------------------------------------------------

def search_pubmed(disease_name):

    url = (
        "https://eutils.ncbi.nlm.nih.gov/"
        "entrez/eutils/esearch.fcgi"
    )

    params = {
        "db": "pubmed",

        "term":
            disease_name,

        "retmode":
            "json",

        "retmax":
            5
    }

    try:

        response = requests.get(
            url,
            params=params,
            timeout=20
        )

        print(
            "PubMed search status:",
            response.status_code
        )

        if response.status_code != 200:

            print(
                "PubMed search error:"
            )

            print(
                response.text
            )

            return []

        data = response.json()

        ids = data[
            "esearchresult"
        ][
            "idlist"
        ]

        print(
            "PubMed IDs:",
            ids
        )

        return ids

    except requests.exceptions.RequestException as e:

        print(
            "PubMed search exception:"
        )

        print(e)

        return []


# --------------------------------------------------
# PubMed: Get article details
# --------------------------------------------------

def get_pubmed_details(
    pubmed_ids
):

    if not pubmed_ids:

        return []

    url = (
        "https://eutils.ncbi.nlm.nih.gov/"
        "entrez/eutils/esummary.fcgi"
    )

    params = {
        "db":
            "pubmed",

        "id":
            ",".join(pubmed_ids),

        "retmode":
            "json"
    }

    try:

        response = requests.get(
            url,
            params=params,
            timeout=20
        )

        print(
            "PubMed details status:",
            response.status_code
        )

        if response.status_code != 200:

            print(
                "PubMed details error:"
            )

            print(
                response.text
            )

            return []

        data = response.json()

        result = data.get(
            "result",
            {}
        )

        articles = []

        for pubmed_id in pubmed_ids:

            article = result.get(
                str(pubmed_id)
            )

            if not article:

                continue

            articles.append({

                "id":
                    article.get(
                        "uid",
                        pubmed_id
                    ),

                "title":
                    article.get(
                        "title",
                        "Title unavailable"
                    ),

                "journal":
                    article.get(
                        "fulljournalname",
                        "Journal unavailable"
                    ),

                "date":
                    article.get(
                        "pubdate",
                        "Date unavailable"
                    )

            })

        print(
            "PubMed articles retrieved:",
            len(articles)
        )

        return articles

    except requests.exceptions.RequestException as e:

        print(
            "PubMed details exception:"
        )

        print(e)

        return []


# --------------------------------------------------
# ClinVar: Search variants for a gene
# --------------------------------------------------

def search_clinvar(
    gene_symbol
):

    url = (
        "https://eutils.ncbi.nlm.nih.gov/"
        "entrez/eutils/esearch.fcgi"
    )

    params = {

        "db":
            "clinvar",

        "term":
            gene_symbol + "[gene]",

        "retmode":
            "json",

        "retmax":
            5
    }

    response = requests.get(
        url,
        params=params
    )

    if response.status_code == 200:

        data = response.json()

        return data[
            "esearchresult"
        ][
            "idlist"
        ]

    else:

        print(
            "ClinVar search error:"
        )

        print(
            response.text
        )

        return []


# --------------------------------------------------
# ClinVar: Get variant details
# --------------------------------------------------

def get_clinvar_details(
    clinvar_ids,
    gene_symbol
):

    if not clinvar_ids:

        return []

    url = (
        "https://eutils.ncbi.nlm.nih.gov/"
        "entrez/eutils/esummary.fcgi"
    )

    params = {

        "db":
            "clinvar",

        "id":
            ",".join(clinvar_ids),

        "retmode":
            "json"
    }

    response = requests.get(
        url,
        params=params
    )

    if response.status_code == 200:

        data = response.json()

        variants = []

        for clinvar_id in clinvar_ids:

            record = data[
                "result"
            ][
                clinvar_id
            ]

            # --------------------------------------
            # Check whether record belongs to gene
            # --------------------------------------

            genes = record.get(
                "genes",
                []
            )

            gene_matches = False

            for gene in genes:

                if gene.get(
                    "symbol"
                ) == gene_symbol:

                    gene_matches = True

                    break

            if not gene_matches:

                continue

            # --------------------------------------
            # Check transcript-based variants
            # --------------------------------------

            variant_title = record.get(
                "title",
                ""
            )

            if (
                variant_title.startswith(
                    "NM_"
                )
                or variant_title.startswith(
                    "NC_"
                )
            ):

                if "(" in variant_title and ")" in variant_title:

                    gene_in_title = (
                        variant_title.split(
                            "(",
                            1
                        )[1].split(
                            ")",
                            1
                        )[0]
                    )

                    if (
                        gene_in_title
                        != gene_symbol
                    ):

                        continue

            # --------------------------------------
            # Clinical significance
            # --------------------------------------

            clinical = record.get(
                "germline_classification",
                {}
            )

            # --------------------------------------
            # Molecular consequence
            # --------------------------------------

            consequences = record.get(
                "molecular_consequence_list",
                []
            )

            # --------------------------------------
            # Variant location
            # --------------------------------------

            variation_set = record.get(
                "variation_set",
                []
            )

            location = []

            if variation_set:

                location = variation_set[
                    0
                ].get(
                    "variation_loc",
                    []
                )

            if location:

                loc = location[0]

            else:

                loc = {}

            # --------------------------------------
            # Store variant information
            # --------------------------------------

            variants.append({

                "clinvar_id":
                    record.get(
                        "uid"
                    ),

                "variant":
                    record.get(
                        "title"
                    ),

                "variant_type":
                    record.get(
                        "obj_type"
                    ),

                "clinical_significance":
                    clinical.get(
                        "description"
                    ),

                "molecular_consequence":
                    (
                        "; ".join(
                            consequences
                        )
                        if consequences
                        else "Not available"
                    ),

                "assembly":
                    loc.get(
                        "assembly_name"
                    ),

                "chromosome":
                    loc.get(
                        "chr"
                    ),

                "start":
                    loc.get(
                        "start"
                    ),

                "stop":
                    loc.get(
                        "stop"
                    ),

                "band":
                    loc.get(
                        "band"
                    )
            })

        return variants

    else:

        print(
            "ClinVar details error:"
        )

        print(
            response.text
        )

        return []


# --------------------------------------------------
# ClinVar: Process one gene
# --------------------------------------------------

def process_clinvar_gene(
    item
):

    gene_symbol = item[
        "target"
    ][
        "approvedSymbol"
    ]

    clinvar_ids = search_clinvar(
        gene_symbol
    )

    print(
        "ClinVar IDs for",
        gene_symbol,
        ":",
        clinvar_ids
    )

    clinvar_variants = get_clinvar_details(
        clinvar_ids,
        gene_symbol
    )

    print(
        "ClinVar variants for",
        gene_symbol,
        ":",
        clinvar_variants
    )

    # Add gene name to each variant
    for variant in clinvar_variants:

        variant[
            "gene"
        ] = gene_symbol

    return clinvar_variants


# --------------------------------------------------
# Main page
# --------------------------------------------------

@app.route(
    "/",
    methods=["GET", "POST"]
)

def home():

    disease = ""

    disease_description = ""

    genes = []

    pubmed_articles = []

    all_clinvar_variants = []

    gene_neighborhoods = {}


    if request.method == "POST":

        disease = request.form[
            "disease"
        ].strip()


        # ------------------------------------------
        # Step 1: Find disease
        # ------------------------------------------

        disease_id = search_disease(
            disease
        )


        if disease_id:

            # --------------------------------------
            # Step 2: NCBI disease description
            # --------------------------------------

            disease_description = (
                get_disease_description(
                    disease_id
                )
            )

            print(
                "Disease description:",
                disease_description
            )


            # --------------------------------------
            # Step 3: PubMed
            # --------------------------------------

            pubmed_ids = search_pubmed(
                disease
            )

            pubmed_articles = (
                get_pubmed_details(
                    pubmed_ids
                )
            )

            print(
                "PubMed articles:",
                pubmed_articles
            )


            # --------------------------------------
            # Step 4: Get associated genes
            # --------------------------------------

            associated_genes = (
                get_disease_genes(
                    disease_id
                )
            )


            # --------------------------------------
            # Step 5: ClinVar
            # --------------------------------------

            # Process the genes concurrently
            # instead of waiting for each gene
            # sequentially.

            with ThreadPoolExecutor(
                max_workers=3
            ) as executor:

                clinvar_results = (
                    executor.map(
                        process_clinvar_gene,
                        associated_genes
                    )
                )

                for clinvar_variants in (
                    clinvar_results
                ):

                    all_clinvar_variants.extend(
                        clinvar_variants
                    )


            # --------------------------------------
            # Step 6: Ensembl genomic locations
            # --------------------------------------

            gene_ids = []

            for item in associated_genes:

                gene_ids.append(
                    item[
                        "target"
                    ][
                        "id"
                    ]
                )

            locations = get_gene_locations(
                gene_ids
            )


            # --------------------------------------
            # Step 7: Gene neighborhood analysis
            # --------------------------------------

            # Process up to 5 neighborhood
            # requests concurrently.

            with ThreadPoolExecutor(
                max_workers=5
            ) as executor:

                neighborhood_results = (
                    executor.map(
                        lambda item:
                            process_gene_neighborhood(
                                item,
                                locations
                            ),
                        associated_genes
                    )
                )

                for (
                    gene_symbol,
                    neighborhood
                ) in neighborhood_results:

                    if neighborhood:

                        gene_neighborhoods[
                            gene_symbol
                        ] = neighborhood


            print(
                "Total genes with neighborhood data:",
                len(
                    gene_neighborhoods
                )
            )


            # --------------------------------------
            # Step 8: Prepare gene results
            # --------------------------------------

            for item in associated_genes:

                target = item[
                    "target"
                ]

                gene_id = target[
                    "id"
                ]

                location = locations.get(
                    gene_id,
                    {}
                )

                genes.append({

                    "symbol":
                        target[
                            "approvedSymbol"
                        ],

                    "name":
                        target[
                            "approvedName"
                        ],

                    "score":
                        item[
                            "score"
                        ],

                    "chromosome":
                        location.get(
                            "seq_region_name"
                        ),

                    "start":
                        location.get(
                            "start"
                        ),

                    "end":
                        location.get(
                            "end"
                        ),

                    "strand":
                        location.get(
                            "strand"
                        ),

                    "assembly":
                        location.get(
                            "assembly_name"
                        ),

                    "browser_url": (

                        "https://genome.ucsc.edu/"
                        "cgi-bin/hgTracks"

                        "?db=hg38"

                        "&position=chr%s:%s-%s"

                        % (

                            location.get(
                                "seq_region_name"
                            ),

                            location.get(
                                "start"
                            ),

                            location.get(
                                "end"
                            )
                        )
                    )
                })


    return render_template(
        "index.html",
        disease=disease,
        disease_description=disease_description,
        genes=genes,
        pubmed_articles=pubmed_articles,
        clinvar_variants=all_clinvar_variants,
        gene_neighborhoods=gene_neighborhoods
    )


# --------------------------------------------------
# Run Flask application
# --------------------------------------------------

if __name__ == "__main__":

    app.run(
        debug=False
    )
