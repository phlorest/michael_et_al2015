import pathlib
from itertools import chain

import newick
import phlorest


def reroot(tree, outgroup='Mawe'):
    """Root `tree` (a newick.Node) on the midpoint of the branch leading to `outgroup`."""
    leaf = tree.get_node(outgroup)
    path = [leaf]
    while path[-1].ancestor:
        path.append(path[-1].ancestor)
    # reverse parent->child edges from the old root down to the outgroup
    for child, parent in reversed(list(zip(path, path[1:]))):
        parent.descendants.remove(child)
        parent.ancestor, parent.length = child, child.length
        child.descendants.append(parent)
    # splice out the old root if it's left with a single child
    old = path[-1]
    if len(old.descendants) == 1:
        (only,) = old.descendants
        old.ancestor.descendants[old.ancestor.descendants.index(old)] = only
        only.ancestor, only.length = old.ancestor, (only.length or 0) + (old.length or 0)
    # new root halfway along the outgroup branch
    (rest,) = leaf.descendants
    half = rest.length / 2
    leaf.descendants, leaf.ancestor, leaf.length = [], None, half
    rest.ancestor, rest.length = None, half
    root = newick.Node.create(descendants=[leaf, rest])
    return root.newick + ';'


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
        )
        args.writer.add_summary(summary, self.metadata, args.log)

        p1 = self.raw_dir.read_trees(
            'TGB9jNoApiakaTuriwaraMrBayes.nex.run1.t',
            burnin=5001, sample=500, detranslate=True)
        p2 = self.raw_dir.read_trees(
            'TGB9jNoApiakaTuriwaraMrBayes.nex.run2.t',
            burnin=5001, sample=500, detranslate=True)
        
        # reroot
        posterior = [
            reroot(t.newick) for t in chain(p1, p2)
        ]
        args.writer.add_posterior(posterior, self.metadata, args.log, rooted=True)
