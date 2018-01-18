
LAMMPS_SCRIPT_NEB_MIN = (
    "# ---- initialize simulations\n"
    "clear\n"
    "units metal\n"
    "dimension 3\n"
    "boundary p p p\n"
    "atom_style charge\n"
    "atom_modify map array\n"
    "\n"
    "# ---- BULK MINIMIZATION\n"
    "# ---- create atoms\n"
    "read_data {bulk_lmps_structure}\n"
    "# ---- define interatomic potential\n"
    "include potential.mod\n"
    "\n"
    "# ---- define settings\n"
    "compute eng all pe/atom\n"
    "compute eatoms all reduce sum c_eng\n"
    "\n"
    "# ---- run minimization\n"
    "reset_timestep 0\n"
    "fix 1 all box/relax iso 0.0 vmax 0.001\n"
    "thermo 10\n"
    "thermo_style custom step pe lx ly lz press pxx pyy pzz c_eatoms\n"
    "min_style cg\n"
    "minimize 1e-25 1e-25 5000 10000\n"
    "\n"
    "variable bulk_natoms equal \"count(all)\"\n"
    "variable bulk_tot_energy equal \"c_eatoms\"\n"
    "variable bulk_length_x equal \"lx\"\n"
    "variable bulk_length_y equal \"ly\"\n"
    "variable bulk_length_z equal \"lz\"\n"
    "variable bulk_ecoh equal \"pe/atoms\"\n"
    "# --- output\n"
    "print \"pyPosMat output section\"\n"
    "print \"tot_energy = ${{bulk_tot_energy}}\"\n"
    "print \"num_atoms = ${{bulk_natoms}}\"\n"
    "print \"latt_const_a = ${{bulk_length_x}}\"\n"
    "print \"latt_const_b = ${{bulk_length_y}}\"\n"
    "print \"latt_const_c = ${{bulk_length_z}}\"\n"
    "print \"ecoh = ${{bulk_ecoh}}\"\n"
    "\n"
    "#------------------------------------------------------------------------------\n"
    "# INIT IMAGE MINIMIZATION\n"
    "#------------------------------------------------------------------------------\n"
    "\n"
    "# ---- initialize simulations\n"
    "clear\n"
    "units metal\n"
    "dimension 3\n"
    "boundary p p p\n"
    "atom_style charge\n"
    "atom_modify map array\n"
    "\n"
    "# --- create atoms\n"
    "read_data {initial_lmps_structure}\n"
    "change_box all x final 0 ${{bulk_length_x}} y final 0 ${{bulk_length_y}} z final 0 ${{bulk_length_z}} remap\n"
    "# ---- define interatomic potential\n"
    "include potential.mod\n"
    "# include modify_atom_info.mod\n"
    "\n"
    "# ---- define settings\n"
    "compute eng all pe/atom\n"
    "compute eatoms all reduce sum c_eng\n"
    "\n"
    "# ---- run minimization\n"
    "reset_timestep 0\n"
    "#fix 1 all box/relax iso 0.0 vmax 0.001\n"
    "thermo 10\n"
    "thermo_style custom step pe lx ly lz press pxx pyy pzz c_eatoms\n"
    "min_style cg\n"
    "minimize 1e-25 1e-25 5000 10000\n"
    "\n"
    "variable init_natoms equal \"count(all)\"\n"
    "variable init_tot_energy equal \"c_eatoms\"\n"
    "variable init_length_x equal \"lx\"\n"
    "variable init_length_y equal \"ly\"\n"
    "variable init_length_z equal \"lz\"\n"
    "variable init_ecoh equal \"pe/atoms\"\n"
    "# --- output\n"
    "print \"pyPosMat output section\"\n"
    "print \"tot_energy = ${{init_tot_energy}}\"\n"
    "print \"num_atoms = ${{init_natoms}}\"\n"
    "print \"latt_const_a = ${{init_length_x}}\"\n"
    "print \"latt_const_b = ${{init_length_y}}\"\n"
    "print \"latt_const_c = ${{init_length_z}}\"\n"
    "print \"ecoh = ${{init_ecoh}}\"\n"
    "print \"lammps_sim_done\"\n"
    "\n"
    "# --- write structure out\n"
    "write_data {lmps_dump_initial} nocoeff\n"
    "#------------------------------------------------------------------------------\n"
    "# --- FINAL IMAGE MINIMIZATION\n"
    "#------------------------------------------------------------------------------\n"
    "# ---- initialize simulations\n"
    "clear\n"
    "units metal\n"
    "dimension 3\n"
    "boundary p p p\n"
    "atom_style charge\n"
    "atom_modify map array\n"
    "\n"
    "# --- create atoms\n"
    "read_data {final_lmps_structure}\n"
    "change_box all x final 0 ${{bulk_length_x}} y final 0 ${{bulk_length_y}} z final 0 ${{bulk_length_z}} remap\n"
    "\n"
    "# ---- define interatomic potential\n"
    "include potential.mod\n"
    "# include modify_atom_info.mod\n"
    "# ---- define settings\n"
    "compute eng all pe/atom\n"
    "compute eatoms all reduce sum c_eng\n"
    "\n"
    "# ---- run minimization\n"
    "reset_timestep 0\n"
    "#fix 1 all box/relax iso 0.0 vmax 0.001\n"
    "thermo 10\n"
    "thermo_style custom step pe lx ly lz press pxx pyy pzz c_eatoms\n"
    "min_style cg\n"
    "minimize 1e-25 1e-25 5000 10000\n"
    "\n"
    "variable natoms equal \"count(all)\"\n"
    "variable tot_energy equal \"c_eatoms\"\n"
    "variable length_x equal \"lx\"\n"
    "variable length_y equal \"ly\"\n"
    "variable length_z equal \"lz\"\n"
    "variable ecoh equal \"pe/atoms\"\n"
    "#variable ecoh equal \"v_etotal/v_atoms\"\n"
    "# --- output\n"
    "print \"pyPosMat output section\"\n"
    "print \"tot_energy = ${{tot_energy}}\"\n"
    "print \"num_atoms = ${{natoms}}\"\n"
    "print \"latt_const_a = ${{length_x}}\"\n"
    "print \"latt_const_b = ${{length_y}}\"\n"
    "print \"latt_const_c = ${{length_z}}\"\n"
    "print \"ecoh = ${{ecoh}}\"\n"
    "\n"
    "# --- write structure_out\n"
    "write_data {lmps_dump_final} nocoeff\n")

LAMMPS_NEB_SCRIPT_FORMAT= (
    "#------------------------------------------------------------------------------\n"
    "# NEB CALCULATION\n"
    "#------------------------------------------------------------------------------\n"
    "# ---- initialize simulations\n"
    "clear\n"
    "units metal\n"
    "dimension 3\n"
    "boundary p p p\n"
    "atom_style charge\n"
    "atom_modify map array\n"
    "\n"
    "# --- create atoms\n"
    "read_data {lmps_dump_initial}\n"
    "include potential.mod\n"
    "#include modify_atom_info.mod\n"
    "\n"
    "# set up neb run\n"
    "timestep 0.001\n"
    "\n"
    "thermo 5\n"
    "min_style quickmin\n"
    "minimize 1.0e-10 1.0e-10 2000 2000\n"
    "\n"
    "fix 1 all neb 1.0\n"
    "\n"
    "# variable u uloop 20\n"
    "# dump 1 all custom 10 *.${{u}}.dump id type x y z\n"
    "\n"
    "neb 0.0 0.1 10000 10000 10 final {neb_dump_final}\n"
    "#------------------------------------------------------------------------------\n")

def lmps_script__neb_min__to_str(
        bulk_filename,
        init_filename,
        final_filename,
        lmps_dump_init='lmps_dump.initial',
        lmps_dump_final='lmps_dump.final'):
    assert type(bulk_filename) is str
    assert type(init_filename) is str
    assert type(final_filename) is str
    assert type(lmps_dump_init) is str
    assert type(lmps_dump_final) is str
    
    _str_out = LAMMPS_SCRIPT_NEB_MIN.format(
            bulk_lmps_structure=bulk_filename,
            initial_lmps_structure=init_filenae,
            final_lmps_structure=final_filename,
            lmps_dump_initial=lmps_dump_init,
            lmps_dump_final=lmps_dump_final)
    return _str_out

def lmps__convert_lmps_dump_to_neb_dump__to_string(
        lmps_dump,
        neb_dump=None):
    if neb_dump is None:
        neb_dump="{}.neb".format(lmps_dump)

    lines=None
    with open(lmps_dump,'r') as f:
        lines=f.readlines()

    n_atoms = int(lines[2].split(' ')[0])

    n_line_begin=17
    n_line_end=n_line_begin+n_atoms
    p_atoms = []
    for i in range(n_line_begin,n_line_end):
        line = lines[i].split(' ')
        atom_id = int(line[0])
        x = float(line[3])
        y = float(line[4])
        z = float(line[5])
        p_atoms.append([atom_id,x,y,z])

    _strout =  str(n_atoms) + "\n"
    for p in p_atoms:
        _strout += " ".join([str(v) for v in p]) + "\n"

    return _strout

def lmps_script__neb_calculation__to_string(lmps_dump_initial,neb_dump_final):
    assert type(lmps_dump_initial) is str
    assert type(neb_dump_final) is str


    _rtnstr = LAMMPS_NEB_SCRIPT_FORMAT.format(
        lmps_dump_initial=lmps_dump_initial,
        neb_dump_final=neb_dump_final)

    return _rtnstr
