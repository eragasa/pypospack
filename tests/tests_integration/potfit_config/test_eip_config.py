import pypospack.potfit as potfit

if __name__ == "__main__":

    # configure the fitting engine
    eip_config = potfit.EipFittingEngine(\
            fname_config_potential = 'pypospack.buckingham.yaml',
            fname_config_qoi = 'pypospack.qoi.yaml',
            fname_config_structures = 'pypospack.structure.yaml')
    

