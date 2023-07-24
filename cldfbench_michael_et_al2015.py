import pathlib
from itertools import chain

import phlorest

try:
    import ete3
except ImportError:
    ete3 = None


def reroot(i, tree):
    tree = ete3.Tree(str(tree), format=0)
    tree.set_outgroup('Mawe')
    return tree.write(format=5)


def fix_tree(p):
    if not isinstance(p, str):
        p = p.read_text(encoding='utf8')
    return p.replace('(percent)', '_percent')


class Dataset(phlorest.Dataset):
    dir = pathlib.Path(__file__).parent
    id = "michael_et_al2015"

    def cmd_makecldf(self, args):
        self.init(args)
        
        summary = self.raw_dir.read_tree(
            'MICHAEL.tupiguarani_phylogeny.tre',
            preprocessor=fix_tree,
            detranslate=True)
        args.writer.add_summary(summary, self.metadata, args.log)

        p1 = self.raw_dir.read_trees(
            'TGB9jNoApiakaTuriwaraMrBayes.nex.run1.t',
            burnin=5001, sample=500, detranslate=True)
        p2 = self.raw_dir.read_trees(
            'TGB9jNoApiakaTuriwaraMrBayes.nex.run2.t',
            burnin=5001, sample=500, detranslate=True)
        
        # reroot
        posterior = [
            reroot(i, t) for i, t in enumerate(chain(p1, p2), 1)
        ]
        args.writer.add_posterior(posterior, self.metadata, args.log, rooted=True)
