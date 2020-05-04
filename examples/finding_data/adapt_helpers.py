import os
def example_adapt_map():
    import urllib.request
    urllib.request.urlretrieve(
        'https://gong.nso.edu/adapt/maps/gong/2020/adapt40311_03k012_202001010000_i00005600n1.fts.gz',
        'adapt20200101.fts.gz'
        )

    if not os.path.exists('adapt20200101.fts'):
        import gzip
        with gzip.open('adapt20200101.fts.gz', 'rb') as f:
            with open('adapt20200101.fts', 'wb') as g:
                g.write(f.read())

    return 'adapt20200101.fts'