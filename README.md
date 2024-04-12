# kct-walk
k-mer walking with kct. Extends a seed without needing to come back to a reference.

Installing (enter Julia package mode by pressing ] in a julia shell located in src)
```
(@v1.8) pkg> activate Modron
  Activating project at `~/Marut/Modron`

(Modron) pkg> instantiate
```

Walking mode
```
usage: Modron.jl walk -s SEED -l LIBRARY -o OUTPUT [-r MAX_RECURSION]
                      [-c MIN_COUNT] [-h]

arguments:
  -s, --seed SEED       The seed to extend by k-mer walking.
  -l, --library LIBRARY
                        The KCT library to use for k-mer walking.
  -o, --output OUTPUT   Output file. Will be in fasta format.
  -r, --max_recursion MAX_RECURSION
                        How much the seed can be extended by k-mer
                        walking. (type: Int64, default: 500)
  -c, --min_count MIN_COUNT
                        The minimum count for a k-mer to be
                        admissible. Should be >0. (type: Int64,
                        default: 1)
  -h, --help            show this help message and exit
```

Kct preparation mode
```
usage: Modron.jl prepare -d DATA -o OUTPUT -j JELLYFISH [-k K] [-h]

arguments:
  -d, --data DATA       Folder that contains ONLY jellyfish files to
                        turn into a Kct.
  -o, --output OUTPUT   The output kct file. Will be in binary format.
  -j, --jellyfish JELLYFISH
                        Path of jellyfish executable
  -k K                  The size of k-mers in jellyfish files (type:
                        Int64, default: 31)
  -h, --help            show this help message and exit
```

**Known issues:**
- Currently only supports building KCTs with 27mers (-k option doesn't do anything right now)
- Building from jf is slow because it relies on dumps instead of binaries
- -d argument is extremely inconveniant. Need to detect type (folder or file for single sample KCT) and only use .jf files in folder.
