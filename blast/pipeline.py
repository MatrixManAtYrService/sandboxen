import conducto as co
import json
from pathlib import Path

# commands.py and experiment.py are in the same folder as this file
from experiment import genomes, genes

img = co.Image(
    image="ncbi/blast",  # use the BLAST image published by ncbi on dockerhub
    copy_dir=".",  # add this directory
    reqs_py=["conducto", "biopython", "pandas"],
    reqs_packages=["wget", "gzip"],
)

data_dir = "/conducto/data/pipeline"


def download_file(source_url, target_path) -> co.Serial:
    "Returns a serial node which downloads a gzipped FASTA file"

    target_dir = Path(target_path).parent

    node = co.Serial()
    node["Download"] = co.Exec(
        f"mkdir -p {target_dir} && wget -O {target_path}.gz {source_url}"
    )
    node["Decompress"] = co.Exec(f"gunzip -c {target_path}.gz > {target_path}")

    return node


def analyze(hits):
    """
    Print a list of genes found in more than one genome
    For more featureful analysis experience, use a notebook instead
    """

    import pandas as pd
    from Bio.Blast import NCBIXML

    hits = json.loads(hits)

    # parse blast output
    blast_records = {}
    for genome, file in hits.items():
        with open(data_dir + "/" + file, "r") as f:
            blast_records[genome] = list(NCBIXML.parse(f))

    # aggregate hits by protein
    hits_by_protein = {}
    for genome, records in blast_records.items():
        for record in records:
            # exclude iterations that found nothing
            if record.alignments:
                try:
                    protein = record.query.split(' ')[0]
                except IndexError:
                    print(f"Couldn't find protein ID: {record.query}")
                hits_by_protein.setdefault(protein, set()).add(genome)

    # as a grid
    genomes_by_number = list(enumerate(sorted(hits.keys())))
    matrix = []
    for protein, found in hits_by_protein.items():
        # reference genes are from s_cerevisiae, expect them
        # include only hits in more than one genome
        if len(found) > 1:
            _hits = [0] * len(hits.keys())
            for i, genome in genomes_by_number:
                if genome in found:
                    _hits[i] = 1
            matrix.append([protein] + _hits)

    # as a data frame
    labels = ["protein"] + [x[1] for x in genomes_by_number]
    df = pd.DataFrame(matrix, columns=labels)
    print(df)


def main() -> co.Serial:

    with co.Serial(image=img) as root:

        with co.Parallel(name="Download") as download:

            # genomes
            for name, url, target_file in genomes(data_dir):
                download["genome: " + name] = download_file(url, target_file)

            # genes
            source_url, genes_file = genes(data_dir)
            download["genes: S288C"] = download_file(source_url, genes_file)

        hits = {}

        with co.Parallel(name="Process"):
            for name, _, target_file in genomes(data_dir):
                with co.Serial(name=name) as process_one:

                    process_one["Make BLAST DB"] = co.Exec(
                        f"""
                         cd {data_dir}
                         makeblastdb -in {target_file} -dbtype nucl -out {name}
                        """
                    )

                    hits_file = f"{name}_hits.xml"
                    process_one["Find Genes"] = co.Exec(
                        f"""
                        cd {data_dir}
                        blastn -query {genes_file} -outfmt 5 -db {name} > {hits_file}
                        """
                    )
                    hits[name] = hits_file

        root["Analyze"] = co.Exec(analyze, json.dumps(hits))
        # root["Analyze"] = co.nb(something???)

    return root


if __name__ == "__main__":
    co.main(default=main)
