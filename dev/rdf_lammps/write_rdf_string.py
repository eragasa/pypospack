
def write_npt_string(
        fix_id = "npt_{1}",
        group_id = "all",
        temp_start = "100."
        temp_stop = "100."
        temp_damp = "$(100.0*dt)"
        pressure_type = "iso"
        pressure_start = "35"
        pressure_stop = "35"
        pressure_damp = "$(1000.0*dt)"
        ):
    assert pressure_type in ["iso","ansio","tri"]
    
    s = "fix {fix_id} npt temp {temp_start} {temp_stop} {temp_damp} {pressure_type} {pressure_start} {pressure_damp}".format(
            fix_id=fix_id,
            group_id=group_id,
            temp_start=temp_start,
            temp_stop=temp_stop,
            temp_damp=temp_damp,
            pressure_type=pressure_type,
            pressure_start=pressure_start,
            pressure_stop=pressure_stop,
            pressure_damp=pressure_damp)

    return s

def write_rdf_string(
        compute_id = "computeRdf",
        group_id = "all",
        n_bins = 100,
        rdf_fn = "computeRdf.rdf",
        fix_id = "fixRdf"
        Nevery =  100,
        Nrepeat = 1
        Nfreq = 100):

    """

    Args:
        compute_id(str): default:computeRdf
        group_id(str): default:all
        n_bins(int): default 100
        file_out(str): default:computeRdf.rdf
        fix_id(str): fixRdf
    """
    s = ["compute {compute_id} {group_id} {n_bins}".format(
                compute_id=compute_id,
                group_id=group_id,
                n_bins=n_bins)]
    s += ["fix {fix_id} all ave/time 100 1 100 c_{compute_id}[*] file {rdf_fn} mode vector".format(
                fix_id=fix_id,
                compute_id=compute_id,
                rdf_fn=rdf_fn)]
    return "\n".join(s)

if __name__ == "__main__":
    print(write_npt_string())
    print(write_rdf_string())
