"""Admitted atomic photoionization targets, verifier side.

The I21 sodium cohort and the I22 atomic pilot (approved 2026-11-26,
integration-tests/I22/reference/atomic_pilot_approval_20261126.json) share
the native owner and these independent auditors. The table below is written
independently of the C++ provider: each D04 parent X admits exactly one
tracked isotope, its singly charged ion and the D04 product string
``X+ + e``. Nothing here is inferred from a native receipt.

Native input receipts name their D04 parent in ``source_parent``. Receipts
written before the I22 pilot have no such field. They can only be Na,
because the native provider then refused every nucleus other than 23Na;
``receipt_element`` applies that documented rule and nothing else.
"""

# D04 parent -> (isotope nucleus, D04 product string of X -> X+ + e)
TARGETS = {
    'Na': ('23Na', 'Na+ + e'),
    'He': ('4He', 'He+ + e'),
    'Ne': ('20Ne', 'Ne+ + e'),
    'Ar': ('40Ar', 'Ar+ + e'),
}

# Which parents each existing test ID may score. No new ID is created.
TEST_ELEMENTS = {'I21': ('Na',), 'I22': ('He', 'Ne', 'Ar')}


def require(element):
    """Return the (nucleus, products) of an admitted D04 parent, else ValueError."""
    if element not in TARGETS:
        raise ValueError('unsupported atomic photoionization parent: ' + repr(element))
    return TARGETS[element]


def nucleus(element):
    return require(element)[0]


def products(element):
    return require(element)[1]


def element_of_nucleus(value):
    """The D04 parent of one isotope nucleus, else ValueError."""
    for element, (isotope, _) in TARGETS.items():
        if isotope == value:
            return element
    raise ValueError('unsupported atomic photoionization nucleus: ' + repr(value))


def receipt_element(receipt):
    """D04 parent recorded by a native input receipt (pre-I22 receipts: Na)."""
    element = receipt.get('source_parent', 'Na')
    require(element)
    channels = receipt.get('source_channels')
    if channels is not None:
        if not isinstance(channels, list) or len(channels) != 1 or not isinstance(channels[0], dict):
            raise ValueError('native receipt must record exactly one source channel')
        if not str(channels[0].get('id', '')).startswith('PHIDRATES:' + element + ':'):
            raise ValueError('native source channel does not belong to the recorded parent ' + element)
    return element


def check_test(test_id, element):
    """Refuse a parent outside the scope of the given existing test ID."""
    if element not in TEST_ELEMENTS.get(test_id, ()):
        raise ValueError(f'{element} photoionization is not in the {test_id} scope')
