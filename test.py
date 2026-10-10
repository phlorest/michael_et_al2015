def test_valid(cldf_dataset, cldf_logger):
    assert cldf_dataset.validate(log=cldf_logger)


def test_phlorest_check(cldf_dataset, cldf_logger):
    from phlorest.check import run_checks
    assert run_checks(cldf_dataset, cldf_logger)


import newick
import pytest

from cldfbench_michael_et_al2015 import reroot


@pytest.mark.parametrize('tree,expected', [
    # unrooted trichotomy, as in the MrBayes output
    ('(A:1,B:2,(C:3,Mawe:4):5);', '(Mawe:2,(C:3,(A:1,B:2):5):2);'),
    # outgroup already a child of an unrooted root
    ('(Mawe:2,A:1,B:1);', '(Mawe:1,(A:1,B:1):1);'),
    # rooted: old root becomes unifurcating and is spliced out
    ('((A:1,B:1):1,(C:1,Mawe:2):3);', '(Mawe:1,(C:1,(A:1,B:1):4):1);'),
    # rooted with outgroup at the root: spliced edge is the one that gets halved
    ('(Mawe:2,(A:1,B:1):3);', '(Mawe:2.5,(A:1,B:1):2.5);'),
])
def test_reroot(tree, expected):
    def norm(nwk):
        return newick.loads(nwk)[0].newick.replace('.0,', ',').replace('.0)', ')')
    assert norm(reroot(newick.loads(tree)[0])) == norm(expected)
