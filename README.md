libgraphlet
===========

libgraphlet is an embeddable graphlet counting library based on [Orca](http://www.biolab.si/supp/orca/orca.html).

* Graphlet counting for graphlets of size 2-5
* Pairwise graphlet degree vector similarity
* Graphlet degree distribution and GDD-agreement

## Requirements ##

Compiling the bundled tools requires

* C++11 compliant compiler
* Boost graph library headers

## Compilation ##

```
git clone https://github.com/SimonLarsen/libgraphlet.git libgraphlet
cd libgraphlet
cmake . && make
```
