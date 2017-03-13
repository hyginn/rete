- This is from running with the [example files on the hotnet2 homepage](https://github.com/raphael-group/hotnet2/tree/master/example). I added the output achieved from running

```
python makeRequiredPPRFiles.py @example/configs/influence_matrix.config
```

followed by

```
python runHotNet2.py @example/configs/simple.config
```

 I also added this readme file and renamed the original readme to `README-hotnet2.md`.

# Output file explanations
The configuration files used to obtain this output can be found under `./config`. Note that I included the output from the second Hotnet2 example file, which can be found in the directories appended with a "2".

## ./influences_matrices{1,2}
- The `permuted networks` folder holds permutations of the current network structure.
- `example_edge_list` is the existing edges in the network.
- `example_index_genes` is the legend for how to map numbers to genes.
- `example_ppr_0.6.h5` is the influence matrix.

## ./output/simple{1,2}/delta_0.{[0-9]*}

- The number in the folder name indicate the threshold used.
- `components.txt` contains the genes of the significantly altered subnetworks. Genes are tab-separated and subnetworks are listed one per line. 
- `significance.txt` includes the expected and actual number of subnetworks of a specific size, and an associated p-value.
- `results.json` contains the output from `significance.txt` under the key `statistics`, the output from `components.txt` under the key `components`, the network size, and the parameter information.
- `heat.json` contains the heat-values for each gene in the input network. 
