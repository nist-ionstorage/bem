from xml.etree import ElementTree as ET

def cpy_colormap(x):
    r = ET.parse(x).getroot()
    c2n, n2c = {}, {}
    for s in r.iterfind("Colors/Style"):
        n = s.find("Name").text
        if n.startswith("TRAPELECTRODE_"):
            n = n[len("TRAPELECTRODE_"):]
        a = s.attrib["Ambient"].split()[:3]
        c = [int(float(ai)*31) for ai in a]
        c = (c[0]<<0) | (c[1]<<5) | (c[2]<<10) | (1<<15)
        n2c[n] = c
        c2n[c] = n
    return n2c, c2n

default_cpy_colormap = {
        53409: "DC1",
        45249: "DC2",
        61728: "DC3",
        42407: "DC4",
        46901: "DC5",
        48245: "DC6",
        63446: "DC7",
        42187: "DC8",
        44071: "DC9",
        48860: "GND",
        56483: "OUT1",
        55674: "OUT2",
        33659: "OUT3",
        63446: "OUT4",
        42187: "OUT5",
        44071: "OUT6",
        59650: "OUT7",
        58444: "OUT9",
        41741: "RF",
        58444: "RF1",
        50601: "RF2",
        59650: "RF3",
        56483: "SEP1",
        55674: "SEP2",
        33659: "SEP3",
        63446: "SEP4",
        42187: "SEP5",
        44071: "SEP6",
        59650: "SEP7",
        58444: "SEP9",
        47454: "TR1A",
        61268: "TR1B",
        50698: "TR1C",
        33689: "TR1D",
        58468: "TR1E",
        33659: "TR1F",
        36574: "TR2A",
        40072: "TR2B",
        45621: "TR2C",
        37156: "TR2D",
        55674: "TR2E",
        48569: "TR3A",
        41230: "TR3B",
        50601: "TR3C",
        45171: "TR3D",
        33607: "TR3E",
        56483: "TR3F",
        59853: "TR4A",
        35933: "TR4B",
        33155: "TR4C",
        45781: "TR4D",
        37156: "TR4E",
        33155: "TR4F",
}


if __name__ == "__main__":
    import sys
    n2c, c2n = cpy_colormap(open(sys.argv[1]))
    for k in sorted(n2c):
        v = n2c[k]
        print bin(v)
        print '%5i: "%s",' % (v, k)
    
