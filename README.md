# WVZAnalysis

### Words for Devs
1. Main looper logic is in `src/mainlooper.jl`
2. Cuts logic for different channel are in `src/ZZZ WZZ WWZ`
3. For histogram we use https://github.com/Moelf/FHist.jl
4. I would recomment to use https://github.com/andyferris/Dictionaries.jl for in-memory "hadd" when combining results from different files.

### Install the dependencies:

## Use Julia package manager
- run `julia` which enters REPL
- press `]` to switch to package mode
- run `dev https://gitlab.cern.ch/jiling/WVZAnalysis.jl/`
you can then find the source files of this package under `~/.julia/dev/WVZAnalysis`

## Already cloned
If you want the source files to be somewhere else, you can manually `git clone` this package. Then navigate
to that loacation:
```
julia --project=.
] instantiate
```
(the dot `.` just means current directory, as usual)
