import pandas as pd

def valid_barcode(x):
    try:
        x.index('TCGA')
        return True
    except:
        return False
def get_barcode(x):
    i = x.index('TCGA')
    return x[i:i+12]
def get_origin(x):
    try:
        if int(x.split('-')[3][:2]) > 9:
            return 'normal'
        else:
            return 'tumor'
    except:
        try:
            if x.split('-')[3] == 'NT' or x.split('-')[3] == 'NB':
                return 'normal'
            else:
                return 'tumor'
        except:
            return 'tumor'


manifest = pd.read_csv('/cellar/users/ramarty/Data/kir/TCGA/manifests/whole_exome.all_samples.tsv', sep='\t')

manifest['valid'] = manifest.filename.apply(valid_barcode)
manifest = manifest[manifest.valid]
manifest['barcode'] = manifest.filename.apply(get_barcode)
manifest['origin'] = manifest.filename.apply(get_origin)
normal = manifest[manifest.origin == 'normal']
normal = normal.drop_duplicates('barcode')
normal_samples = list(normal.id)
normal_barcodes = list(normal.barcode)

def get_barcodes():
    return normal_barcodes