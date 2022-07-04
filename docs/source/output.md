# Output

Finished folder `FINISHED_PATH / project_name / timestamp` will contain:

1. `query_files/*` - directory containing all input query files.
2. `mmseqs2_search_results.m8`
3. `alignments.json` - results of alignment search implemented in `utils.search_alignments.py`
4. `metadata*` - files with some useful info
5. `results*` - multiple files from DeepFRI. Organized by model type ['GCN' / 'CNN'] and its mode ['mf', 'bp', 'cc', 'ec'] for the total of 8 files.
```
mf = molecular_function
bp = biological_process
cc = cellular_component
ec = enzyme_commission
```

Sometimes results from one model can be missing which means that all query protein sequences were aligned correctly or none of them were aligned.
