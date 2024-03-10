Code for skimming some skims and creating columnar datasets from them.

To make new skims do:
```bash
conda activate diphoton-env
python analysis/skim.py <dType> <analysis_region> --submit_batches
```
This submits jobs to condor. The outputs are