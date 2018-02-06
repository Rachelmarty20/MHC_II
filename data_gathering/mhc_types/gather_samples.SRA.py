import os

directory = '/data/nrnb03/dbGaP/SRA/ARIC_GoESP'
samples = [x.split('.')[0] for x in os.listdir(directory)]
with open('/cellar/users/ramarty/Data/hla_ii/hla_types/samples.SRA.txt', 'w') as outfile:
    for sample in samples:
        outfile.write('{0}\n'.format(sample))