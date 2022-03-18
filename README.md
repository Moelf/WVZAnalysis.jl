# WVZAnalysis

### Words for Devs
1. Main looper logic is in `src/mainlooper.jl`
2. Cuts logic for different channel are in `src/ZZZ WZZ WWZ`
3. For histogram we use https://github.com/Moelf/FHist.jl
4. I would recomment to use https://github.com/andyferris/Dictionaries.jl for in-memory "hadd" when combining results from different files.

### Install the dependencies:
```
export JULIA_LOAD_PATH=.
julia
] instantiate
```
