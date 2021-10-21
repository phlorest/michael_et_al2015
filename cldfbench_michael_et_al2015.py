import pathlib

import ete3
from nexus import NexusReader
import phlorest


class Dataset(phlorest.Dataset):
    dir = pathlib.Path(__file__).parent
    id = "michael_et_al2015"

    def cmd_makecldf(self, args):
        self.init(args)
        with self.nexus_summary() as nex:
            for tree in reroot(self.sample(
                self.raw_dir.read('MICHAEL.tupiguarani_phylogeny.tre'),
                n=1,
                detranslate=True,
                strip_annotation=True,
            )):
                self.add_tree(args, tree, nex, 'summary')
                break

        with self.nexus_posterior() as nex:
            for i, tree in enumerate(reroot(self.run_nexus(
                'combine',
                self.sample(
                    self.remove_burnin(
                        self.raw_dir.read('TGB9jNoApiakaTuriwaraMrBayes.nex.run1.t'), 5001),
                    n=500,
                    detranslate=True,
                    strip_annotation=True,
                ),
                self.sample(
                    self.remove_burnin(
                        self.raw_dir.read('TGB9jNoApiakaTuriwaraMrBayes.nex.run2.t'), 5001),
                    n=500,
                    detranslate=True,
                    strip_annotation=True,
                ),
            )), start=1):
                self.add_tree(args, tree, nex, 'posterior-{}'.format(i))


def reroot(in_):
    for i, tree in enumerate(NexusReader.from_string(in_).trees.trees):
        tree = ete3.Tree(tree.newick_string, format=0)
        tree.set_outgroup('Mawe')
        yield phlorest.NexusTree('tree tg_%d [&R] = %s' % (i, tree.write(format=5)))
