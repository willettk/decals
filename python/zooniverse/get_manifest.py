
def get_manifest(key_data):
    pass
    # t['metadata.retire_at'] = [40] * N
    # t['survey'] = ['decals'] * N

    # TODO what are these counters for? Still needed?
    # t['metadata.counters.feature'] = np.zeros(N, dtype=int)
    # t['metadata.counters.smooth'] = np.zeros(N, dtype=int)
    # t['metadata.counters.star'] = np.zeros(N, dtype=int)

    # make the image file locs - will likely not be needed now for Panoptes!
    # url_stub = "http://www.galaxyzoo.org.s3.amazonaws.com/subjects/decals_dr2"
    # for imgtype in ('standard', 'inverted', 'thumbnail'):
    #     lst = url_stub + '/{0}/'.format(imgtype) + nsa_decals['iauname'] + '.jpeg'
    #     t['location.{0}'.format(imgtype)] = lst

