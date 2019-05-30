import pytest

def get_2body_pairs(symbols):

    pair = []

    for i,s_i in enumerate(symbols):
        for j,s_j in enumerate(symbols):
            if i<=j:
                pair.append("{}{}".format(s_i,s_j))

    return pair
def get_3body_triples(symbols):

    triples = []
    for i,s_i in enumerate(symbols):
        for j,s_j in enumerate(symbols):
            for k,s_k in enumerate(symbols):
                triples.append("{}{}{}".format(s_i,s_j,s_k))

    return triples

def dev__get_2body_pairs__1_element():
    symbols = ['Si']

    pairs = get_2body_pairs(symbols)

    print(pairs)

def dev__get_2body_pairs__2_element():
    symbols = ['Si','Ge']

    pairs = get_2body_pairs(symbols)

    print(pairs)

def dev__get_2body_pairs__3_element():
    symbols = ['Si','Ge','Sn']

    pairs = get_2body_pairs(symbols)

    print(pairs)

def dev__get_3body_triples__1_element():
    symbols = ['Si']

    triples = get_3body_triples(symbols)

    print(triples)

def dev__get_3body_triples__2_element():
    symbols = ['Si','Ge']

    triples = get_3body_triples(symbols)
    
    print(triples)

def dev__get_3body_triples__3_element():
    symbols = ['Si','Ge','Sn']

    triples = get_3body_triples(symbols)
    print(triples)
if __name__ == "__main__":
    dev__get_2body_pairs__1_element()
    dev__get_2body_pairs__2_element()
    dev__get_2body_pairs__3_element()
    dev__get_3body_triples__1_element()
    dev__get_3body_triples__2_element()
    dev__get_3body_triples__3_element()
