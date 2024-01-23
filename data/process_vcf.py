#!python3
"""Processing VCF for easier use in numpy."""

import click
import numpy as np
import pandas as pd
from tqdm import tqdm
from cyvcf2 import VCF


def extract_geno(vcf, af=0.01):
    assert (af > 0.0) and (af < 1.0)
    samples = np.array(vcf.samples)
    afs = []
    pos = []
    geno = []
    for variant in tqdm(vcf):
        if len(variant.ALT) == 1:
            if np.min([variant.aaf, 1.0 - variant.aaf]) >= af:
                afs.append(variant.aaf)
                pos.append(variant.POS)
                geno.append(variant.gt_types.copy())
    afs = np.array(afs)
    pos = np.array(pos)
    geno = np.vstack(geno)
    return {"samples": np.array(samples), "position": pos, "af": afs, "geno": geno}


@click.command()
@click.option(
    "--input",
    "-i",
    required=True,
    type=click.Path(exists=True),
    help="Input VCF file.",
)
@click.option(
    "--af",
    "-a",
    required=False,
    default=0.01,
    type=float,
    help="Allele frequency threshold.",
)
@click.option(
    "--threads",
    "-t",
    required=False,
    default=4,
    type=int,
    help="Number of threads for VCF processing.",
)
@click.option(
    "--out",
    "-o",
    required=True,
    type=str,
    default="test.npz",
    help="Output file for genotypes.",
)
def main(
    input,
    af,
    threads,
    out,
):
    vcf = VCF(input, gts012=True, threads=threads)
    geno_dict = extract_geno(vcf, af=af)
    np.savez_compressed(out, **geno_dict)


if __name__ == "__main__":
    main()
