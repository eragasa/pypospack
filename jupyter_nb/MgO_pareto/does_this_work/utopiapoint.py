

def get_closest_to_utopia_point(d_vector, n = 1):
    """
    return an index of the dvector closest to the utopia point
    """
    idx_shortest = np.where(d_vector = d_vector.min())
    idx_shortest_pop = np.argpartition(d_vector,n)[:n_potentials]

    return idx_shortest, idx_shortest_pop


