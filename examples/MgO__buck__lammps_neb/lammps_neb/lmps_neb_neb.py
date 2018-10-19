def lmps__convert_lmps_dump_to_neb_dump__to_string(lmps_dump,neb_dump=None):
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
        "#------------------------------------------------------------------------------\n"

    _rtnstr = LAMMPS_NEB_SCRIPT_FORMAT.format(
        lmps_dump_initial=lmps_dump_initial,
        neb_dump_final=neb_dump_final)

    return _rtnstr
