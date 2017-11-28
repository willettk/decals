from astroquery.vizier import Vizier


def get_expert_catalog(row_limit=-1):
    catalog_list = Vizier.find_catalogs('J/ApJS/186/427')
    Vizier.ROW_LIMIT = row_limit
    return Vizier.get_catalogs(catalog_list.keys())[0]

if __name__ == '__main__':
    row_limit = 20
    # catalog_dir = ''
    expert_catalog = get_expert_catalog(row_limit)
    print(expert_catalog)