---
GFE Driver
---

The GFE (Graph Framework Evaluation) Driver is the program used to run the experiments in the paper TBP,  measuring the throughput of updates in libraries supporting structural dynamic graphs and the completion times of the [Graphalytics kernels](https://github.com/ldbc/ldbc_graphalytics). The driver supports the following systems: [Teseo](https://github.com/cwida/teseo), [LLAMA](https://github.com/goatdb/llama), [GraphOne](https://github.com/the-data-lab/GraphOne) and [Stinger](http://stingergraph.com/). It can run three kinds experiments: insert all edges in a random permuted order from an input graph, execute the updates specified by a [graphlog file](https://github.com/whatsthecraic/graphlog) and run the kernels of the Graphalytics suite: BFS, PageRank (PR), local triangle counting (LCC), weighted shortest paths (SSSP), weakly connected components (WCC) and community detection through label propagation (CDLP).  

### Build 

#### Requisites 
- O.S. Linux
- Autotools, [Autoconf 2.69+](https://www.gnu.org/software/autoconf/)
- A C++17 compliant compiler with support for OpenMP. We tested it on Clang 10 and GCC 10.
- libnuma 2.0 +
- [libpapi 5.5 +](http://icl.utk.edu/papi/)
- [SQLite 3.27 +](https://sqlite.org)

#### Configure

From the source, use `autoreconf -iv` to create the `configure` script.
The driver needs to be linked with the system to evaluate, which needs to be built ahead. We do not recommend linking the driver with multiple systems at once, due to the usage of global variables in some systems and other naming clashes. Instead, it is safer to reconfigure and rebuild the driver each time for a single specific system.


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
Differently from the other systems, we used GCC, rather than Clang, to compile Stinger, due to low level incompatibilities.

Configure the GFE driver with:

```
/path/to/source/gfe_driver/configure --enable-optimize --disable-debug --with-stinger=/path/to/stinger/build
```

##### LLAMA

Use the branch `feature/gfe `, it contains additional patches w.r.t. [upstream](https://github.com/goatdb/llama), from https://github.com/whatsthecraic/llama:  

```
git clone https://github.com/whatsthecraic/llama -b feature/gfe
```

LLAMA is a header-only library. It does not need to be compiled in advance.

Configure the GFE driver with:

```
/path/to/source/gfe_driver/configure --enable-optimize --disable-debug --with-llama=/path/to/llama
```


##### GraphOne

Use the branch `feature/gfe `, it contains additional patches w.r.t. [upstream](https://github.com/the-data-lab/GraphOne), from https://github.com/whatsthecraic/GraphOne:  

```
git clone https://github.com/whatsthecraic/GraphOne -b feature/gfe
cd GraphOne
mkdir build && cd build
cmake -S ../ -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=clang++
make -j
```
If the build has been successful, it should at least create the executable `graphone64`. Then, configure the driver with:

```
/path/to/source/gfe_driver/configure --enable-optimize --disable-debug --with-graphone=/path/to/graphone/build
```

##### Teseo

Use the branch `master` from https://github.com/cwida/teseo: 

```
git clone https://github.com/cwida/teseo
cd teseo
./autoreconf -iv
mkdir build && cd build
../configure --enable-optimize --disable-debug CXX=clang++ 
make -j
```

If the build has been successful, it should at least create the archive `libteseo.a`. Then configure the driver with:

```
/path/to/source/gfe_driver/configure --enable-optimize --disable-debug --with-teseo=/path/to/teseo/build
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
./gfe_driver -G /path/to/input/graph.properties -u --log /path/to/updates.graphlog --aging_timeout -l <system_to_evaluate> -w <num_threads> -d output_results.sqlite3
```

- **Graphalytics**: execute the six kernels from the Graphalytics suite. Add the option `-R <N>` to repeat `N` times the execution of all Graphalytics kernels, one after the other. E.g., to run the kernels five times, after all vertices and edges have been inserted, use:

```
./gfe_driver -G /path/to/input/graph.properties -u -l <system_to_evaluate> -w <num_threads> -R 5 -d output_results.sqlite3
```

Type `./gfe_driver -h` for the full list of options and for the libraries that can be evaluated (option `-l`). The driver spawns the number of threads given by the option `-w` to concurrently run all insertions or updates. For Graphalytics, it defaults to the total number of the physical threads in the machine. This setting can be changed with the option `-r <num_threads>`. Note that the numbers
in the library codes (e.g. teseo.**6**, stinger**3**) are unrelated to the versions of the systems evaluated, they were only used
internally for development purposes.

The database `output_results.sqlite3` will contain the final results. Use the [view definitions](https://github.com/whatsthecraic/gfe_results/blob/master/views.sql) from [this Jupyter repository](https://github.com/whatsthecraic/gfe_results) to interpret the data. In particular, the notebooks [insert_only_plot.ipynb](https://github.com/whatsthecraic/gfe_results/blob/master/insert_only_plot.ipynb), [updates_data.ipynb](https://github.com/whatsthecraic/gfe_results/blob/master/updates_throughput.ipynb) and [graphalytics_after_inserts.ipynb](https://github.com/whatsthecraic/gfe_results/blob/master/graphalytics_after_inserts.ipynb) are the notebooks used to make the final plots in the paper, for the insertions, the updates and the Graphalytics kernels, respectively. 


### Repeating the experiments

These are the full commands to repeat the experiments in the paper:

##### Insertions only (Figure 6)
```bash
for NT in 1 2 4 6 8 10 12 14 16 18 20 40; do
    # Stinger, source code: library/stinger/{stinger.hpp, stinger_unsafe.cpp} 
    ./gfe_driver -G /path/to/input/graph.properties -u -l stinger3-ref -w $NT -d results.sqlite3
    # LLAMA, source code: library/llama/llama_class.*
    ./gfe_driver -G /path/to/input/graph.properties -u -l llama6-ref --build_frequency 10s -w $NT -d results.sqlite3
    # GraphOne, source code: library/graphone/*
    ./gfe_driver -G /path/to/input/graph.properties -u -l g1_v4-ref-ignore-build -w $NT -d results.sqlite3
    # Teseo, source code: library/teseo/teseo_driver.*
    ./gfe_driver -G /path/to/input/graph.properties -u -l teseo.6 -w $NT -d results.sqlite3
done
```

##### Updates (Figure 7)

```bash
for NT in 1 2 4 6 8 10 12 14 16 18 20 40; do
    # Stinger, source code: library/stinger/{stinger.hpp, stinger_unsafe.cpp} 
    ./gfe_driver -G /path/to/input/graph.properties -u --log /path/to/log/graph.graphlog --aging_timeout -l stinger3-ref -w $NT -d results.sqlite3
    # LLAMA, source code: library/llama/llama_class.*
    ./gfe_driver -G /path/to/input/graph.properties -u --log /path/to/log/graph.graphlog --aging_timeout -l llama6-ref --build_frequency 10s -w $NT -d results.sqlite3
    # GraphOne, source code: library/graphone/*
    ./gfe_driver -G /path/to/input/graph.properties -u --log /path/to/log/graph.graphlog --aging_timeout -l g1_v4-ref-ignore-build -w $NT -d results.sqlite3
    # Teseo, source code: library/teseo/teseo_driver.*
    ./gfe_driver -G /path/to/input/graph.properties -u --log /path/to/log/graph.graphlog --aging_timeout -l teseo.6 -w $NT -d results.sqlite3
done
```

The option `--aging_timeout` serves to limit the total time to execute the experiment to up to four hours. For LLAMA, it could be necessary to further anticipate the timeout, as the continuous creation of new deltas can cause a memory exhaustion. Currently, this can only be done by directly altering the source code, from the method `Aging2Master::do_run_experiment` in `experiment/details/aging2_master.cpp`.  

##### Graphalytics (Figure 9)

```
# Stinger, source code: library/stinger/*
./gfe_driver -G /path/to/input/graph.properties -u -l stinger3-ref -w 1 -R 5 -d results.sqlite3
# LLAMA, source code: library/llama/*
./gfe_driver -G /path/to/input/graph.properties -u -l llama6-ref --build_frequency 10s -w 16 -R 5 -d results.sqlite3
# GraphOne, source code: library/graphone/*
./gfe_driver -G /path/to/input/graph.properties -u -l g1_v4-ref-ignore-build -w 3 -R 5 -d results.sqlite3
# Teseo (logical vertices), source code: library/teseo/teseo_driver.*
export OMP_PLACES="sockets" OMP_PROC_BIND="spread"
./gfe_driver -G /path/to/input/graph.properties -u -l teseo.6 -w 40 -R 5 -d results.sqlite3
# Teseo (real vertices), source code: library/teseo/teseo_real_vtx.*
export OMP_PLACES="sockets" OMP_PROC_BIND="spread"
./gfe_driver -G /path/to/input/graph-dense.properties -u -l teseo-dv.6 -w 40 -R 5 -d results.sqlite3
# Teseo (LCC custom), source code: library/teseo/teseo_lcc.*
export OMP_PLACES="sockets" OMP_PROC_BIND="spread"
./gfe_driver -G /path/to/input/graph.properties -u -l teseo-lcc.6 -w 40 -R 5 --blacklist="bfs,cdlp,pagerank,sssp,wcc" -d results.sqlite3

```

The graphs `graph-dense.properties` are analogous to their corresponding `graph.properties`, but with the vertices relabelled into a dense domain. These graphs are included in the archive loaded in [Zenodo](https://zenodo.org/record/3966439).   
