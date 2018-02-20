import for_sandor_settings as settings
from a_download_decals.get_catalogs.get_joint_nsa_decals_catalog import get_nsa_catalog, get_decals_bricks
from a_download_decals.get_decals_images_and_catalogs import get_decals


def main():
    """
    Create the NSA-DECaLS-GZ catalog, download fits, produce png, and identify new subjects

    Returns:
        None
    """

    # specify execution options
    settings.new_catalog = True
    settings.new_images = True
    settings.overwrite_fits = False
    settings.overwrite_png = False
    settings.run_to = None

    nsa = get_nsa_catalog(settings.nsa_catalog_loc, settings.nsa_version)
    print('nsa loaded')
    bricks = get_decals_bricks(settings.bricks_loc, settings.data_release)
    print('bricks loaded')

    galaxies_to_include = {
        'J104031.24+121739.8',
        'J103846.59+054145.0',
        'J112147.02+040710.8',
        'J011341.82-000609.7',
        'J101819.60+070252.9',
        'J091135.57+325055.6',
        'J135406.26+052122.8',
        'J103351.36-003340.8',
        'J132132.19+121115.8',
        'J140808.50-014208.2',
        'J133729.36+040615.9',
        'J090339.90+032211.3',
        'J225246.39+010758.0',
        'J021219.80-004835.0',
        'J112736.73+002342.7'
    }

    # print(nsa['iauname'] in galaxies_to_include)
    nsa['selected'] = list(map(lambda x: x in galaxies_to_include, nsa['iauname']))
    nsa_filtered = nsa[nsa['selected']]
    print(len(nsa_filtered))

    _ = get_decals(nsa_filtered, bricks, settings)


if __name__ == '__main__':
    main()
