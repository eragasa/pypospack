from pypospack.potential import EamEmbeddingFunction

class EamEmbeddingEquationOfState(EamEmbeddingFunction):
    def __init__(self,
            symbols,
            potential_type="eam_embed_eos"):

        EamEmbeddingFunction.__init__(self,
                symbols=symbols,
                potential_type=potential_type)

        self.is_eos = True
        self.density_fn = None
        self.pair_fn = None
        self.r_cut = None

