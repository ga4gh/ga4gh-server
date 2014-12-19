"""
Test module.
"""
import itertools


def powerset(iterable, max_sets=None):
    """
    powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)

    See https://docs.python.org/2/library/itertools.html#recipes
    """
    s = list(iterable)
    return itertools.islice(itertools.chain.from_iterable(
        itertools.combinations(s, r) for r in range(len(s) + 1)),
        0, max_sets)

