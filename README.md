---
GFE Driver
---

The GFE (Graph Framework Evaluation) Driver is the program used to run the experiments in the paper TBP,  measuring the throughput of updates in libraries supporting structural dynamic graphs and the completion times of the [Graphalytics kernels](https://github.com/ldbc/ldbc_graphalytics). The driver supports the following systems: [Teseo](https://github.com/cwida/teseo), [LLAMA](https://github.com/goatdb/llama), [GraphOne](https://github.com/the-data-lab/GraphOne), [Stinger](http://stingergraph.com/) and [LiveGraph](https://github.com/thu-pacman/LiveGraph-Binary). It can run three kinds experiments: insert all edges in a random permuted order from an input graph, execute the updates specified by a [graphlog file](https://github.com/whatsthecraic/graphlog) and run the kernels of the Graphalytics suite: BFS, PageRank (PR), local triangle counting (LCC), weighted shortest paths (SSSP), weakly connected components (WCC) and community detection through label propagation (CDLP).  

### Build 

#### Requisites 
- O.S. Linux
- Autotools, [Autoconf 2.69+](https://www.gnu.org/software/autoconf/)
- A C++17 compliant compiler with support for OpenMP. We tested it on Clang 10 and GCC 10.
- libnuma 2.0 +
- [libpapi 5.5 +](http://icl.utk.edu/papi/)
- [SQLite 3.27 +](https://sqlite.org)

#### Configure

Initialise the sources and the configure script by:

```
git clone https://github.com/cwida/gfe_driver
cd gfe_driver
git submodule update --init
autoreconf -iv 
```

The driver needs to be linked with the system to evaluate, which has to be built ahead. We do not recommend linking the driver with multiple systems at once, due to the usage of global variables in some systems and other naming clashes. Instead, it is safer to reconfigure and rebuild the driver each time for a single specific system.


##### Stinger
Use the branch `feature/gfe `, it contains additional patches w.r.t. [upstream](https://github.com/stingergraph/stinger), from https://github.com/whatsthecraic/stinger: 

```
git clone https://github.com/whatsthecraic/stinger -b feature/gfe
cd stinger
mkdir build && cd stinger
cmake ../ -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=0 
make
```
If the build has been successful, it should at least create the executable `bin/stinger_server`.

Configure the GFE driver with:

```
mkdir build && cd build
../configure --enable-optimize --disable-debug --with-stinger=/path/to/stinger/build
```

##### LLAMA

Use the branch `feature/gfe `, it contains additional patches w.r.t. [upstream](https://github.com/goatdb/llama), from https://github.com/whatsthecraic/llama:  

```
git clone https://github.com/whatsthecraic/llama -b feature/gfe
```

LLAMA is a header-only library. It does not need to be compiled in advance.

Configure the GFE driver with:

```
mkdir build && cd build
../configure --enable-optimize --disable-debug --with-llama=/path/to/llama
```


##### GraphOne

Use the branch `feature/gfe `, it contains additional patches w.r.t. [upstream](https://github.com/the-data-lab/GraphOne), from https://github.com/whatsthecraic/GraphOne:  

```
git clone https://github.com/whatsthecraic/GraphOne -b feature/gfe
cd GraphOne
mkdir build && cd build
cmake -S ../ -DCMAKE_BUILD_TYPE=Release
make -j
```
If the build has been successful, it should at least create the executable `graphone64`. Then, configure the driver with:

```
mkdir build && cd build
../configure --enable-optimize --disable-debug --with-graphone=/path/to/graphone/build
```

##### LiveGraph

Download the binary library from the [official repository](https://github.com/thu-pacman/LiveGraph-Binary/releases). In the paper we evaluated version 20200829.
Then configure the driver by pointing the path to where the library has been downloading:

```
mkdir build && cd build
../configure --enable-optimize --disable-debug --with-livegraph=/path/to/livegraph/lib
```

##### Teseo

Use the branch `master` from https://github.com/cwida/teseo: 

```
git clone https://github.com/cwida/teseo
cd teseo
./autoreconf -iv
mkdir build && cd build
../configure --enable-optimize --disable-debug
make -j
```

If the build has been successful, it should at least create the archive `libteseo.a`. Then configure the driver with:

```
mkdir build && cd build
../configure --enable-optimize --disable-debug --with-teseo=/path/to/teseo/build   
```

#### Compile

Once configured, run `make -j`. There is no `install` target, the final artifact is the executable `gfe_driver`. 

If in the mood of running the testsuite, type `make check -j`.


### Datasets

In our experiments, we used the following input graphs and data sets:

- `dota-league` and `graph500-SF`, with `SF` in {22, 24 26}, were taken from the [official Graphalytics collection](https://www.graphalytics.org/datasets).
- `uniform-SF`, with `SF` in {22, 24, 26} were generated with an [ad-hoc tool](https://github.com/whatsthecraic/uniform_graph_generator). These are synthetic graphs having the same number of vertices and edges of `graph500-SF`, but a uniform node degree distribution.
- The logs for the experiments with updates, i.e. with both insertions and deletions, were generated with another [ad-hoc tool](https://github.com/whatsthecraic/graphlog). 

A complete image of all datasets used in the experiments can be downloaded from Zenodo: [input graphs](https://zenodo.org/record/3966439), [graph logs](https://zenodo.org/record/3967002). Note that we ran the experiments on a batch of 15 equal machines. Each machine had a different log for the updates, that is, initialised with a different random seed.  


### Executing the driver


The driver takes as input a list of options together with a graph, and emits the results into a sqlite3 database.
There are three kinds of experiments that can be executed:

- **Insertions only**: insert all vertices and edges from an input graph, in a random order. Use the command:

```
./gfe_driver -G /path/to/input/graph.properties -u -l <system_to_evaluate> -w <num_threads> -d output_results.sqlite3
```

For LLAMA only: add the option `--build_frequency 10s` to asynchronously issue the creation of a new level (or delta) every 10 seconds.

- **Updates**: perform all insertions and deletions from a log. Add the option --log /path/to/updates.graphlog :

```
./gfe_driver -G /path/to/input/graph.properties -u --log /path/to/updates.graphlog --aging_timeout 24h -l <system_to_evaluate> -w <num_threads> -d output_results.sqlite3
```

- **Graphalytics**: execute the six kernels from the Graphalytics suite. Add the option `-R <N>` to repeat `N` times the execution of all Graphalytics kernels, one after the other. E.g., to run the kernels five times, after all vertices and edges have been inserted, use:

```
./gfe_driver -G /path/to/input/graph.properties -u -l <system_to_evaluate> -w <num_threads> -R 5 -d output_results.sqlite3
```

Type `./gfe_driver -h` for the full list of options and for the libraries that can be evaluated (option `-l`). The driver spawns the number of threads given by the option `-w` to concurrently run all insertions or updates. For Graphalytics, it defaults to the total number of the physical threads in the machine. This setting can be changed with the option `-r <num_threads>`. Note that the numbers
in the library codes (e.g. teseo.**6**, stinger**3**) are unrelated to the versions of the systems evaluated, they were only used
internally for development purposes.

The database `output_results.sqlite3` will contain the final results. Refer to [this repository](https://github.com/whatsthecraic/gfe_notebooks) to see how to load and inspect the data within Jupyter notebooks and how to recreate the same plots of the paper.

### Repeating the experiments

These are the full commands to repeat the experiments in the paper:

##### Insertions only (Figure 6)
```bash
for NT in 1 2 4 6 8 10 12 14 16 18 20 40; do
    # Stinger, source code: library/stinger/{stinger.hpp, stinger_unsafe.cpp} 
    ./gfe_driver -G /path/to/input/graph.properties -u -l stinger5-ref -w $NT -d results.sqlite3
    # LLAMA, source code: library/llama/llama_class.*
    ./gfe_driver -G /path/to/input/graph.properties -u -l llama6-ref --build_frequency 10s -w $NT -d results.sqlite3
    # GraphOne, source code: library/graphone/*
    ./gfe_driver -G /path/to/input/graph.properties -u -l g1_v4-ref-ignore-build -w $NT -d results.sqlite3
    # LiveGraph, source code: library/livegraph/*
    ./gfe_driver -G /path/to/input/graph.properties -u -l livegraph_ro -w $NT -d results.sqlite3
    # Teseo, source code: library/teseo/teseo_driver.*
    ./gfe_driver -G /path/to/input/graph.properties -u -l teseo.10 -w $NT -d results.sqlite3
done
```

##### Updates (Figure 7)

```bash
for NT in 1 2 4 6 8 10 12 14 16 18 20 40; do
    # Stinger, source code: library/stinger/{stinger.hpp, stinger_unsafe.cpp} 
    ./gfe_driver -G /path/to/input/graph.properties -u --log /path/to/log/graph.graphlog -l stinger5-ref -w $NT -d results.sqlite3
    # LLAMA, source code: library/llama/llama_class.*
    ./gfe_driver -G /path/to/input/graph.properties -u --log /path/to/log/graph.graphlog -l llama6-ref --build_frequency 10s --aging_timeout 4h -w $NT -d results.sqlite3
    # GraphOne, source code: library/graphone/*
    ./gfe_driver -G /path/to/input/graph.properties -u --log /path/to/log/graph.graphlog -l g1_v4-ref-ignore-build -w $NT -d results.sqlite3
    # LiveGraph, source code: library/livegraph/*
    ./gfe_driver -G /path/to/input/graph.properties -u --log /path/to/log/graph.graphlog -l livegraph_ro -w $NT -d results.sqlite3
    # Teseo, source code: library/teseo/teseo_driver.*
    ./gfe_driver -G /path/to/input/graph.properties -u --log /path/to/log/graph.graphlog -l teseo.10 -w $NT -d results.sqlite3
done
```

The option `--aging_timeout` serves to limit the total time to execute the experiment. For LLAMA, it could be necessary to stop the experiment earlier, as the continuous creation of new deltas can cause a memory exhaustion.  
For the experiment with the memory footprint of Figure 7d, add the arguments: `--aging_memfp --aging_memfp_physical --aging_memfp_threshold 330G --aging_release_memory=false`. The option `--aging_memfp` records the memory footprint as the experiment proceeds, `--aging_memfp_physical` records the physical memory (RSS) of the process, rather than the virtual memory of the glibc allocator, `--aging_memfp_threshold 330G` terminates the experiment if the memory footprint measured is greater than 330 GB and `--aging_release_memory=false` avoids releasing the memory used in the driver to load the graph from the file, as it may (or may not) recycled by the libraries. With the memory footprint, for LLAMA, it's not necessary to set `--aging_timeout 4h` as `--aging_memfp_threshold 330G` already acts as a guard on the overall memory consumption.

##### Graphalytics (Table 3)

```
# CSR, source code: library/baseline/csr.*
./gfe_driver -G /path/to/input/graph.properties -u -l csr --load -R 5 -d results.sqlite3
# CSR, LCC (opt), source code: library/baseline/csr.*
./gfe_driver -G /path/to/input/graph.properties -u -l csr-lcc --load -R 5 -d results.sqlite3
# Stinger, source code: library/stinger/*
./gfe_driver -G /path/to/input/graph.properties -u -l stinger5-ref -w 40 -R 5 -d results.sqlite3
# LLAMA, source code: library/llama/*
./gfe_driver -G /path/to/input/graph.properties -u -l llama6-ref --build_frequency 10s -w 16 -R 5 -d results.sqlite3
# GraphOne, source code: library/graphone/*
./gfe_driver -G /path/to/input/graph.properties -u -l g1_v4-ref-ignore-build -w 12 -R 5 -d results.sqlite3
# LiveGraph, source code: library/livegraph/*
./gfe_driver -G /path/to/input/graph.properties -u -l livegraph_ro -w 20 -R 5 -d results.sqlite3
# Teseo (logical vertices), source code: library/teseo/teseo_driver.*
./gfe_driver -G /path/to/input/graph.properties -u -l teseo.10 -w 40 -R 5 -d results.sqlite3
./gfe_driver -G /path/to/input/graph.properties -u -l teseo-lcc.10 -w 40 -R 5 -d results.sqlite3 # LCC (opt) only
# Teseo (real vertices), source code: library/teseo/teseo_real_vtx.*
./gfe_driver -G /path/to/input/graph-dense.properties -u -l teseo-dv.10 -w 40 -R 5 -d results.sqlite3
# Teseo (LCC opt), source code: library/teseo/*
./gfe_driver -G /path/to/input/graph.properties -u -l teseo-lcc.10 -w 40 -R 5 --blacklist="bfs,cdlp,pagerank,sssp,wcc" -d results.sqlite3 
./gfe_driver -G /path/to/input/graph-dense.properties -u -l teseo-lcc-dv.10 -w 40 -R 5 --blacklist="bfs,cdlp,pagerank,sssp,wcc" -d results.sqlite3 

```

The graphs `graph-dense.properties` are analogous to their corresponding `graph.properties`, but with the vertices relabelled into a dense domain. These graphs are included in the archive loaded in [Zenodo](https://zenodo.org/record/3966439). 


##### Sequential and random scans (Figure 9)

```
make bm

# The tool already assumes the graphs are undirected.

# CSR
./bm -G /path/to/graph500-24.properties -l csr -t 1,2,4,6,8,12,16,20,40
./bm -G /path/to/uniform-24.properties -l csr -t 1,2,4,6,8,12,16,20,40
# Stinger
./bm -G /path/to/graph500-24.properties -l stinger -t 1,2,4,6,8,12,16,20,40
./bm -G /path/to/uniform-24.properties -l stinger -t 1,2,4,6,8,12,16,20,40
# LLAMA
./bm -G /path/to/graph500-24.properties -l llama -t 1,2,4,6,8,12,16,20,40
./bm -G /path/to/uniform-24.properties -l llama -t 1,2,4,6,8,12,16,20,40
# GraphOne
./bm -G /path/to/graph500-24.properties -l graphone -t 1,2,4,6,8,12,16,20,40
./bm -G /path/to/uniform-24.properties -l graphone -t 1,2,4,6,8,12,16,20,40
# LiveGraph
./bm -G /path/to/graph500-24.properties -l livegraph-ro -t 1,2,4,6,8,12,16,20,40
./bm -G /path/to/uniform-24.properties -l livegraph-ro -t 1,2,4,6,8,12,16,20,40
# Teseo
./bm -G /path/to/graph500-24.properties -l teseo -t 1,2,4,6,8,12,16,20,40
./bm -G /path/to/uniform-24.properties -l teseo -t 1,2,4,6,8,12,16,20,40

```

At the end of each execution, the tool `bm` stores the results in a json file under /tmp. Check [the notebook bm.nb for Mathematica](https://github.com/whatsthecraic/gfe_notebooks/blob/master/bm.nb) to see an example on how to load and interpret the data.