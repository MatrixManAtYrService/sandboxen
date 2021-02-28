import conducto as co

img = co.Image(
    image="ncbi/blast",  # use the BLAST image published by ncbi on dockerhub
    copy_dir=".",        # add this directory
    reqs_py=["conducto"],# add these tools
    reqs_packages=["wget", "gzip"],
)

data = "/conducto/data/pipeline"

dummy_contents = """> Dummy Data
GATTACA
"""

def main() -> co.Serial:

    with co.Serial(image=img) as root_node:
        with co.Serial(name="setup"):
            co.Exec(f"""
                    wget -O {data}/genes.fasta.gz \
                    https://sgd-prod-upload.s3.amazonaws.com/S000208654/orf_coding.20150113.fasta.gz
                    """,
                    name="get data")

            co.Exec(f"gunzip -c {data}/genes.fasta.gz > {data}/genes.fna",
                    name="decrompress")

            co.Exec(f"echo '{dummy_contents}' > {data}/genome.fna",
                    name="place data")


        with co.Parallel(name="experiment"):
            co.Exec(f"""
                    makeblastdb -in {data}/genome.fna -dbtype nucl -out tempdb
                    blastn -query {data}/genes.fna -outfmt 5 -db tempdb > /dev/null 2> >(tee errors)

                    # fail if previous command wrote to stderr
                    [[ $(wc -l < errors) - ge 1) ]] && exit 1 || exit 0
                    """,
                    name="has errors")

            co.Exec(f"""
                    # fix bad characters for YMR156C, YCL018W, YGR257C, and YDR412W
                    cat {data}/genes.fna \
                        | sed 's/&#/BADCHAR/g' > {data}/fixed/fna

                    makeblastdb -in {data}/genome.fna -dbtype nucl -out tempdb
                    blastn -query {data}/fixed.fna -outfmt 5 -db tempdb > /dev/null 2> >(tee errors)

                    # fail if previous command wrote to stderr
                    [[ $(wc -l < errors) - ge 1) ]] && exit 1 || exit 0
                    """,
                    name="fixed")
    return root_node


if __name__ == "__main__":
    co.main(default=main)
