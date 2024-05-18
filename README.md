# Computing the twisted $L^2$-Euler characteristic

This is a supplement to [arXiv:2310.07024](https://arxiv.org/pdf/2310.07024), in which we introduce an algorithm to compute the twisted $L^2$-Euler characteristic.

## Requirements

Make sure you have [Python 3](https://www.python.org/) and [SageMath](https://www.sagemath.org/) installed. Then, install Regina and SnapPy:

```
$ pip install regina
$ pip install snappy
```

Note that this fails on some systems, since Python packages are managed by the OS.
In that case, you have to install them through your package manager (e.g. there could be a package called `python-regina`),
or you have to pass the `--break-system-packages` option to `pip`.

Next, install the GAP package [HAP](https://www.gap-system.org/Packages/hap.html):

```
$ sage -i gap_packages
```

On some installations this fails for the same reasons, but there should be a corresponding system package such as `gap-packages`.

Finally, clone this repository in a directory of your choice:

```
$ cd /path/to/directory
$ git clone https://github.com/floatingpoint-754/twisted-l2-characteristic
$ cd twisted-l2-characteristic
```

## Contents

This repo provides a directory `twisted_l2` (that can be imported as a Python package), which implements the algorithm,
and several Sage notebooks which showcase the experiments carried out in the paper. It is recommended to use the JupyterLab notebook:

```
$ sage -n jupyterlab
```

A browser tab will open, from which you can open the notebooks (and do your own experimentations).

Some notebooks come with precomputed information about the chain complex (in the `data` directory).
This provides consistency when non-deterministic algorithms are involved.

## Details

The main methods provided by `twisted_l2` are `twisted_l2_characteristic` (alias `characteristic`) and `von_neumann_rank`.
They compute approximations to the twisted $L^2$-Euler characteristic of a chain complex, and to the von Neumann rank of a matrix.

We provide an interface to build the equivariant chain complex from the _face lattice_, i.e. the information of
which $i$-cells are in the boundary of which $(i+1)$-cells. You can either provide the face lattice explicitly or compute it from a Regina/SnapPy triangulation:

```python
import twisted_l2
import regina
import snappy

twisted_l2.load_hap()

mfd = snappy.Manifold("v1539(5,1)")
tri = regina.Triangulation3(mfd.filled_triangulation())
fl = twisted_l2.regina_tri_to_face_lattice(tri, ideal=False)

cw = twisted_l2.cw_complex(fl)
chain_complex = twisted_l2.equivariant_cc(cw, gap=False)

G = twisted_l2.get_fundamental_group(chain_complex)
cc = twisted_l2.get_differentials(chain_complex)
```

You can find more information in the notebooks and the wiki; you can also access the documentation of any object by using the `?` operator in Sage and JupyterLab:

```
twisted_l2.characteristic?
```

> **_Note:_** The documentation is a work in progress, so some entries may be missing.

