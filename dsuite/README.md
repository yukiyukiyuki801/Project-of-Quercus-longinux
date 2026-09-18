# Dsuite Patterson's D

The outgroup is *Quercus glauca* sample DH01 (SRA run SRR27151155). Set the
*Q. robur* PM1N reference path, then map and call the outgroup at the ingroup
RAD loci:

```bash
export QROBUR_FASTA=/path/to/Qrob_PM1N.fa
sbatch --export=ALL,QROBUR_FASTA dsuite/scripts/map_qglauca_outgroup.slurm
```

Run Dsuite 0.5 r58 with minimum outgroup depths of 3, 5, and 10:

```bash
export DSUITE_BIN=/path/to/Dsuite
bash dsuite/scripts/run_dsuite_depths.sh
```

The mapping step requires SRA Toolkit, BWA 0.7.17, SAMtools, and BCFtools
1.21. Archived summaries and QC files are in `outputs/`.
