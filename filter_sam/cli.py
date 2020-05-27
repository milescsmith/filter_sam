"""Console script for filter_sam."""
import sys
import click
from Bio import SeqIO
from typing import Optional
import pysam

# @main.command()
@click.command()
@click.option(
    "--fastq",
    "-q",
    help="FASTQ file",
    default=None,
    required=False,
    type=str,
)
@click.option(
    "--fasta",
    "-a",
    help="FASTA file",
    default=None,
    required=False,
    type=str,
)
@click.option(
    "--sam",
    "-s",
    help="SAM or BAM file",
    default=None,
    required=True,
    type=str,
)
@click.option(
    "--prefix",
    "-p",
    help="Output prefix",
    default="",
    required=False,
    type=str,
)
@click.help_option()
def main(
    fasta: Optional[str],
    fastq: Optional[str],
    sam: str,
    prefix: Optional[str],
    ) -> None:
    """Console script for filter_sam."""
    
    if fasta and fastq:
        print("Both a FASTA and FASTQ file were supplied and I only need one.  Preferentially loading the FASTQ...") 
        seq_file = fastq
        seq_file_type = "fastq"
    elif fasta:
        seq_file = fasta
        seq_file_type = "fasta"
    else:
        seq_file = fastq
        seq_file_type = "fastq"

    fastq_records = [_ for _ in SeqIO.parse(seq_file, seq_file_type)]
    fastq_record_names =[_.name for _ in fastq_records]

    sam_records = pysam.AlignmentFile(sam, "rb")
    filtered_sam_records = pysam.AlignmentFile(f"{prefix}filtered.sam", "wb", template=sam_records)

    for _ in sam_records.fetch():
        if _.query_name in fastq_record_names:
            filtered_sam_records.write(_)
    
    filtered_sam_records.close()
    sam_records.close()


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
