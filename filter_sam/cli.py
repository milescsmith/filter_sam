"""Console script for filter_sam."""
import sys
import click
from Bio import SeqIO
from typing import Optional
import pysam

# @main.command()
@click.command()
@click.option(
    "--fasta", "-a", help="FASTA file", default=None, required=True, type=str,
)
@click.option(
    "--sam", "-s", help="SAM or BAM file", default=None, required=True, type=str,
)
@click.option(
    "--prefix", "-p", help="Output prefix", default="", required=False, type=str,
)
@click.help_option()
def main(
    fasta: Optional[str], fastq: Optional[str], sam: str, prefix: Optional[str],
) -> None:
    """Remove from a SAM/BAM file entries that do not have a corresponding match
    in a FASTA/FASTQ file
    \f
    Reads in a Isoseq3 polished FASTA, a GMAP-processed SAM/BAM file, and writes
    a SAM file with only those records that have a match in the the FASTA
    \f
    
    Parameters
    ----------
    
    fasta : `str`
        FASTA file from the Isoseq3 polishing step
    sam : `str`
        SAM/BAM file from the Isoseq3 GMAP mapping step
    prefix : `str`, optional
        File name prefix to place in from "filtered.sam" for the output SAM file.
    
    Returns
    -------
    
    `None`

    """

    fastq_records = [_ for _ in SeqIO.parse(fasta, "fasta")]
    fastq_record_names = [_.name for _ in fastq_records]

    sam_records = pysam.AlignmentFile(sam, "rb")
    filtered_sam_records = pysam.AlignmentFile(
        f"{prefix}filtered.sam", "w", template=sam_records
    )

    for _ in sam_records.fetch():
        if _.query_name in fastq_record_names:
            filtered_sam_records.write(_)

    filtered_sam_records.close()
    sam_records.close()


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
