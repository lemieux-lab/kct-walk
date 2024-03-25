# kct-walk
k-mer walking with kct. Extends a seed without needing to come back to a reference.

Walking mode
```
usage: Modron.jl walk -s SEED -l LIBRARY -o OUTPUT [-r MAX_RECURSION]
                      [-c MIN_COUNT] [-h]

optional arguments:
  -s, --seed SEED       The seed to extend by k-mer walking.
  -l, --library LIBRARY
                        The KCT library to use for k-mer walking.
  -o, --output OUTPUT   Output file. Will be in fasta64format.
  -r, --max_recursion MAX_RECURSION
                        How much the seed can be extended by k-mer
                        walking. (type: Int64, default: 500)
  -c, --min_count MIN_COUNT
                        The minimum count for a k-mer to be
                        admissible. Should be >0. (type: Int64,
                        default: 1)
  -h, --help            show this help message and exit
```

Kct preparation mode (Currently not working properly)
```
usage: Modron.jl prepare -d DATA -o OUTPUT [-h]

optional arguments:
  -d, --data DATA      Folder with jellyfish files to build a KCT
                       from.
  -o, --output OUTPUT  The output kct file. Will be in binary format.
  -h, --help           show this help message and exit
```
