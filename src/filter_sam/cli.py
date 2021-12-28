"""Console script for filter_sam."""
import sys
from Bio import SeqIO
from typing import Optional
import typer
import pysam

from . import __version__


app = typer.Typer(
    name="filter_sam",
    help=f"Remove from a SAM/BAM file entries that do not have a corresponding "
        "match in a FASTA/FASTQ file",
)


def version_callback(value: bool):
    """Prints the version of the package."""
    if value:
        print(
            f"filter_sam version: {__version__}"
        )
    raise typer.Exit()


@app.command()
def main(
    sam: str = typer.Option(..., "-s", "--sam", help="SAM or BAM file"),
    fasta: Optional[str] = typer.Option(None, "-a", "--fasta", help="FASTA file"),
    prefix: Optional[str] = typer.Option(None, "-p", "--prefix", help="Output prefix"),
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

    fasta_records = [_ for _ in SeqIO.parse(fasta, "fasta")]
    fasta_record_names = [_.name for _ in fasta_records]

    sam_records = pysam.AlignmentFile(sam, "rb")
    filtered_sam_records = pysam.AlignmentFile(
        f"{prefix}filtered.sam", "w", template=sam_records
    )

    for _ in sam_records.fetch():
        if _.query_name in fasta_record_names:
            filtered_sam_records.write(_)

    filtered_sam_records.close()
    sam_records.close()


if __name__ == "__main__":
    typer.run(main)  # pragma: no cover
