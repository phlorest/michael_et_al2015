import pathlib

import ete3
from nexus import NexusReader
import phlorest


class Dataset(phlorest.Dataset):
    dir = pathlib.Path(__file__).parent
    id = "michael_et_al2015"

    def cmd_makecldf(self, args):
        """
summary.trees: original/MICHAEL.tupiguarani_phylogeny.tre
	nexus trees -t -c $< -o tmp
	python process.py tmp $@
	rm tmp

posterior.trees: original/TGB9jNoApiakaTuriwaraMrBayes.nex.run1.t original/TGB9jNoApiakaTuriwaraMrBayes.nex.run2.t
	nexus trees --detranslate --deltree 1-5001 $(word 1, $^) --random 500 -o 1.tmp
	nexus trees --detranslate --deltree 1-5001 $(word 2, $^) --random 500 -o 2.tmp
	nexus_combine_nexus.py 1.tmp 2.tmp  # creates 'combined.nex'
	python process.py combined.nex $@
	rm -rf *.tmp combined.nex
        """
        self.init(args)
        with self.nexus_summary() as nex:
            for tree in cleanup(self.sample(
                self.raw_dir.read('MICHAEL.tupiguarani_phylogeny.tre'),
                n=1,
                detranslate=True,
                strip_annotation=True,
            )):
                self.add_tree(args, tree, nex, 'summary')
                break

        with self.nexus_posterior() as nex:
            for i, tree in enumerate(cleanup(self.run_nexus(
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


def cleanup(in_):
    nex = NexusReader.from_string(in_)

    for i, tree in enumerate(nex.trees.trees):
        # make tree into newick for ete3
        tree = nex.trees.trees[i].split(" = ")[1].strip().lstrip()
        tree = tree.replace("[&U]", "")  # remove unrooted flag if present
        tree = ete3.Tree(tree, format=0)

        # reroot
        tree.set_outgroup('Mawe')
        yield phlorest.NexusTree('tree tg_%d [&R] = %s' % (i, tree.write(format=5)))
