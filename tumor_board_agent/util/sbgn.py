def make_square(s):
    return { 'w': s, 'h': s }

def get_default_bbox(clazz):
    if clazz in sbgn_defaults:
        defs = sbgn_defaults[clazz]
    else:
        defs = unspecified_glyph_defs
    # import pdb; pdb.set_trace()
    bbox = { 'x': 0, 'y': 0, 'w': defs['w'], 'h': defs['h'] }
    return bbox

sbgn_defaults = {
    'macromolecule': make_square(30),
    'complex': make_square(60),
    'state variable': make_square(15),
    'unit of information': make_square(15),
}

unspecified_glyph_defs = make_square(30)
