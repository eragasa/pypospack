from collections import OrderedDict

reference_potentials = OrderedDict()
# Stillinger and Weber,  Phys. Rev. B, v. 31, p. 5262, (1985)
reference_potentials['SW'] = OrderedDict([
    ('SiSiSi_epsilon',2.1686),
    ('SiSiSi_sigma',2.0951),
    ('SiSiSi_a',1.80),
    ('SiSiSi_lambda',21.0),
    ('SiSiSi_gamma',1.20),
    ('SiSiSi_costheta0',-1/3),
    ('SiSiSi_A',7.049556277),
    ('SiSiSi_B',0.6022245584),
    ('SiSiSi_p',4.0),
    ('SiSiSi_q',0.0),
    ('SiSiSi_tol',0.0)
])
reference_potentials['VBWM'] = OrderedDict([
    ('SiSiSi_epsilon',1.64833),
    ('SiSiSi_sigma',2.0951),
    ('SiSiSi_a',1.80),
    ('SiSiSi_lambda',31.5),
    ('SiSiSi_gamma',1.20),
    ('SiSiSi_costheta0',-1/3),
    ('SiSiSi_A',7.049556277),
    ('SiSiSi_B',0.6022245584),
    ('SiSiSi_p',4.0),
    ('SiSiSi_q',0.0),
    ('SiSiSi_tol',0.0)
])
reference_potentials['PG'] = OrderedDict([
    ('SiSiSi_epsilon',1.04190),
    ('SiSiSi_sigma',2.128117),
    ('SiSiSi_a',1.80),
    ('SiSiSi_lambda',31.0),
    ('SiSiSi_gamma',1.10),
    ('SiSiSi_costheta0',-1/3),
    ('SiSiSi_A',19.0),
    ('SiSiSi_B',0.65),
    ('SiSiSi_p',3.5),
    ('SiSiSi_q',0.5),
    ('SiSiSi_tol',0.0)
])

if __name__ == "__main__":
    import os
    reference_fn = 'pyposmat.reference.in'
    parameter_names = [
        'SiSiSi_epsilon','SiSiSi_sigma','SiSiSi_a',
        'SiSiSi_lambda','SiSiSi_gamma','SiSiSi_costheta0',
        'SiSiSi_A','SiSiSi_B','SiSiSi_p','SiSiSi_q',
        'SiSiSi_tol'
    ]
    names = ['sim_id'] + parameter_names
    types = ['sim_id'] + len(parameter_names)*['param']

    with open(reference_fn,'w') as f:
        f.write(",".join(names) + "\n")
        f.write(",".join(types) + "\n")
        for pot_name,pot_parameters in reference_potentials.items():
            f.write(",".join(
                [pot_name] \
                + [str(pot_parameters[k]) for k in parameter_names]
                ) + "\n"
            )
