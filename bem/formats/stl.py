from collections import OrderedDict
import struct

import numpy as np


def read_stl(fil, ignore_attr_length=True):
    """reads binary STL in file "fil" with Inventor style
    (color hidden in attribute length), returns array of normals, array
    of vertices and list or attributes

    Inventor >= 2013 supports colors (bottom right in the export options
    dialog)
    
    The atribute short in binary STL is used by some tools to hide
    color information. The format is 15 bit RGB or BGR, MSB is validity
    (default color if MSB is 0) COLOR=\\xbf\\xbf\\xbf\\xff in the header
    is the overall color default color.
    """
    h = fil.read(80)
    assert h.startswith("STLB")
    n, = struct.unpack("< I", fil.read(4))
    if ignore_attr_length:
        dtype = np.dtype([
            ("normal", "<f4", (3,)),
            ("triangle", "<f4", (3, 3)),
            ("attribute", "<u2")])
        data = np.fromfile(fil, dtype, n)
        normals, triangles = data["normal"], data["triangle"]
        attributes = data["attribute"]
    else: # actual STL binary format with variable length attribute bytes
        vectors, attributes = [], []
        while True:
            try:
                vectors.append(struct.unpack("< 12f", fil.read(48)))
                na, = struct.unpack("< H", fil.read(2))
                attributes.append(fil.read(na))
            except struct.error:
                break
        assert len(vectors) == n, (len(vectors), n)
        vectors = np.array(vectors)
        normals = vectors[:, :3]
        triangles = vectors[:, 3:].reshape(-1, 3, 3)
    return normals, triangles, attributes


def check_normals(normals, triangles):
    """verifies that given vertices are right-handed around normals"""
    a, b, c = triangles.transpose((1, 0, 2))
    n = np.cross(b-a, c-a)
    n /= np.sqrt(np.sum(np.square(n), axis=1))[:, None]
    assert np.allclose(n, normals, rtol=1e-3, atol=1e-10)


def stl_to_mesh(normals, triangles, attributes, scale=None, rename=None):
    """generates a {name: [(points, triangles)]} mesh from the
    stl arrays. For further treatment, use e.g:

    >>> s = read_stl(open(filename, "rb"))
    >>> check_normals(*s[:2])
    >>> r = stl_to_mesh(*s)
    >>> del r["stl_0"] # delete the color-less faces
    >>> print r.keys()
    ['stl_18761', 'stl_12943']
    >>> r = split_by_normal(r)
    >>> m = Mesh.from_mesh(r)
    >>> m.triangulate("qQ")
    >>> m.to_vtk("test_stl_import")
    
    rename can be a renaming dictionary mapping stl color numbers to
    electrode names. if None the color name is "stl_%i" %
    color_numer, else if the color is not found in rename,
    the electrode is dropped.
    """
    d = OrderedDict()
    for a, v in zip(attributes, triangles):
        d.setdefault(a, []).append(v)
    o = OrderedDict()
    for a, v in d.iteritems():
        v = np.array(v)
        i = np.arange(0, v.shape[0]*3, 3)
        points = v.reshape(-1, 3).astype(np.double)
        if scale:
            points /= scale
        triangles = np.c_[i, i+1, i+2].astype(np.intc)
        if rename:
            if a not in rename:
                print "dropping", a
                continue
            else:
                n = rename[a]
        else:
            n = "stl_%i" % a
        o[n] = [(np.ascontiguousarray(points),
            np.ascontiguousarray(triangles))]
    return o


if __name__ == "__main__":
    import sys
    r = read_stl(open(sys.argv[1], "rb"))
    m = stl_to_mesh(*r)
    print m
