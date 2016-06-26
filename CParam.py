"""
Controls whether to include cython or python modules
"""

cython_dblr = True
if cython_dblr == False:
    import DBLR as DB
else:
    import PxDBLR as DB