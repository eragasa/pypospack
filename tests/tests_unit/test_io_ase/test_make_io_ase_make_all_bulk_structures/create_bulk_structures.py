import pypospack.io.ase as io_ase

structures = [ [ ['Ni'], 'fcc', 3.05, 'cubic'],
               [ ['Ni'], 'bcc', 3.05, 'cubic'],
               [ ['Ni'], 'hcp', 3.05, 'ortho'],
               [ ['Ni'], 'hcp', 3.05, 'prim'],
               [ ['Ni'], 'sc', 3.05, 'cubic'],
               [ ['Ni'], 'dia', 3.05, 'cubic'] ]

for s in structures:
    io_ase.make_bulk_structure(\
            symbols=s[0],
            xtal_name=s[1],
            a=s[2],
            lattice_shape=s[3])

