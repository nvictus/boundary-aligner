import matplotlib as mpl
import seaborn as sns; sns.set_style('white')
import matplotlib.pyplot as plt


def load_loops(filepath):
    loops = pandas.read_csv(
        filepath,
        sep=b'\t',
        compression='gzip',
        usecols=[0,1,2,3,4,5],
        header=False,
        names=['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2'])
    loops['chrom1'] = loops.chrom1.apply(lambda x: 'chr{}'.format(x))
    loops['chrom2'] = loops.chrom2.apply(lambda x: 'chr{}'.format(x))
    return loops

def load_domains(filepath):
    domains = pandas.read_csv(
        filepath,
        sep=b'\t',
        compression='gzip',
        usecols=[0,1,2],
        header=False,
        names=['chrom', 'start', 'end'])
    domains['chrom'] = domains.chrom.apply(lambda x: 'chr{}'.format(x))
    return domains


### EXAMPLE ###

gm_domains = chrom_sorted(
    load_domains('<FILE>'), 
    sort_by='start')

imr90_domains = chrom_sorted(
    load_domains('<FILE>'), 
    sort_by='start')

seq1 = imr90_domains[imr90_domains.chrom=='chr22']['start'].values
seq2 = gm_domains[gm_domains.chrom=='chr22']['start'].values
gap_cost = 1000000

path, cost = boundary_align('chr22', seq1, seq2, gap_cost)


f = plt.figure(figsize=(30,5))
ax = f.add_subplot(111)
ax.plot(seq1, np.zeros(len(seq1)), 'k.')
ax.plot(seq2, np.ones(len(seq2)), 'k.')

pairs = np.c_[
    seq1[match[:,0]],
    seq2[match[:,1]]
]

for a,b in pairs:
    ax.plot([a, b], [0, 1], '.-')
ax.set_ylim([-0.5, 1.5])





