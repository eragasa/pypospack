from pypospack.potential import Potential

class EamEmbeddingFunction(Potential):
    def __init__(self,
            symbols,
            potential_type='eamembed'):

        Potential.__init__(self,
                symbols=symbols,
                potential_type=potential_type,
                is_charge=False)

        self.embedding_evaluations = None
