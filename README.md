# kct-walk
k-mer walking with kct. Extends a seed without needing to come back to a reference.

Installing (enter Julia package mode by pressing ] in a julia shell)
```
(@v1.8) pkg> activate Modron
  Activating project at `~/Marut/Modron`

(Modron) pkg> instantiate
```

Walking mode
```
usage: Modron.jl walk -s SEED -l LIBRARY -o OUTPUT [-r MAX_RECURSION]
                      [-c MIN_COUNT] [-h]

optional arguments:
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

Kct preparation mode (Currently only working on single file using the backup method, slow and dirty).
```
usage: Modron.jl prepare -d DATA -o OUTPUT [-h]

optional arguments:
  -d, --data DATA      Jellyfish file to turn into a Kct.
                       from.
  -o, --output OUTPUT  The output kct file. Will be in binary format.
  -h, --help           show this help message and exit
```
