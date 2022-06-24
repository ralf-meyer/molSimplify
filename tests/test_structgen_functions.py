from molSimplify.Scripts.structgen import smartreorderligs


def test_smartreorderligs():
    """Expected behavior: First order by denticity, then by number of atom"""
    indices = smartreorderligs(['water']*6, [1]*6)
    assert indices == [0, 1, 2, 3, 4, 5]

    indices = smartreorderligs(['water', 'ammonia', 'water', 'water',
                                'ammonia', 'water'], [1]*6)
    assert indices == [0, 2, 3, 5, 1, 4]

    indices = smartreorderligs(['ammonia']*3 + ['water']*3, [1]*6)
    assert indices == [3, 4, 5, 0, 1, 2]

    # 5 monodentates of different sizes
    indices = smartreorderligs(['furan', 'ammonia', 'pyridine', 'carbonyl',
                                'water'], [1]*5)
    assert indices == [3, 4, 1, 0, 2]

    # Test bidentates
    indices = smartreorderligs(['acac', 'acac', 'en'], [2, 2, 2])
    assert indices == [2, 0, 1]

    indices = smartreorderligs(['en', 'en', 'acac'], [2, 2, 2])
    assert indices == [0, 1, 2]

    indices = smartreorderligs(['water', 'carbonyl', 'acac'], [1, 1, 2])
    assert indices == [2, 1, 0]

    # Tetradentate
    indices = smartreorderligs(['water', 'porphirine', 'misc'], [1, 4, 1])
    assert indices == [1, 0, 2]
