import pypospack.crystal as crystal

class LammpsStructure(crystal.SimulationCell):
    def __init__(self,obj=None):
        crystal.SimulationCell.__init__(self,obj)

    def write(self, filename, symbol_list=None, atom_style=None):
        if symbol_list is None:
            symbol_list = self.symbols

        if atom_style is None:
            atom_style = 'charge'

        total_number_of_atoms      = self.n_atoms
        total_number_of_atom_types = len(self.symbols)
        a0 = self.a0

        xlo                        = 0.0
        xhi                        = self.H[0,0] * a0
        ylo                        = 0.0
        yhi                        = self.H[1,1] * a0
        zlo                        = 0.0
        zhi                        = self.H[2,2] * a0
        xy                         = self.H[0,1] * a0
        xz                         = self.H[0,2] * a0
        yz                         = self.H[1,2] * a0

        file = open(filename,'w')

        # this is the header section
        file.write("# {}\n".format(symbol_list))
        file.write("\n")
        file.write("{} atoms\n".format(total_number_of_atoms))
        file.write("{} atom types\n".format(total_number_of_atom_types))
        file.write("\n")
        file.write("{:10.4f} {:10.4f} xlo xhi\n".format(xlo, xhi))
        file.write("{:10.4f} {:10.4f} ylo yhi\n".format(ylo, yhi))
        file.write("{:10.4f} {:10.4f} zlo zhi\n".format(zlo, zhi))
        file.write("\n")
        file.write("{:10.4f} {:10.4f} {:10.4f} xy xz yz\n".format(xy,xz,yz))
        file.write("\n")
        file.write("Atoms\n")
        file.write("\n")

        atom_id = 1
        for i_symbol, symbol in enumerate(symbol_list):
            for i_atom, atom in enumerate(self.atomic_basis):
                if (atom.symbol == symbol):
                    chrg = 1.  # dummy variable
                    posx = self.H[0,0]*atom.position[0]*a0
                    posy = self.H[1,1]*atom.position[1]*a0
                    posz = self.H[2,2]*atom.position[2]*a0
                    if atom_style == 'atomic':
                        str_out = "{} {} {:10.4f} {:10.4f} {:10.4f}\n"
                        str_out = str_out.format(atom_id,
                                                 i_symbol + 1,
                                                 posx, posy, posz)
                    elif atom_style == 'charge':
                        str_out = "{} {} {:10.4f} {:10.4f} {:10.4f} {:10.4f}\n"
                        str_out = str_out.format(atom_id,
                                                 i_symbol + 1,
                                                 chrg,
                                                 posx, posy, posz)
                    file.write(str_out)
                    atom_id += 1
        file.close()

def write_lammps_structure_file(simulation_cell, filename):
    assert isinstance(simulation_cell=crystal.SimulationCell)

    structure_str = ""

    with open(filename,'w') as f:
        f.write(structure_str)
