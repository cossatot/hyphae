

connectivity = {
    "CC": ["LT"],
    "LT": ["CC", "TP", "TR"],
    "TP": ["LT", "TR"],
    "TR": ["LT", "TP", "PP"],
    "PP": ["TP", "NB", "BK"],
    "NB": ["PP", "BK"],
    "BK": ["NB", "PP"],
}

def test_get_connected_subgraphs():
    sgs = get_connected_subgraphs(connectivity)
