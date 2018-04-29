def convert_surface_energy_units(src_x,src_unit,dst_unit):
    surface_energy_units = [
        'electronvolts_per_angstrom_squared',
        "millijoules_per_meter_squared"
    ]
    if all([
            src_unit == 'electronvolts_per_angstrom_squared',
            dst_unit == "millijoules_per_meter_squared"
        ]):
        return src_x * 16021.766
    elif all([
            src_unit == "millijoules_per_meter_squared",
            dst_unit == 'electronvolts_per_angstrom_squared'
        ]):
        return src_x * 6.2415091e-5
    else:
        if src_unit not in surface_energy_units:
            err_msg = "bad src_unit"
            raise ValueError(err_msg)

        if dst_unit not in surface_energy_units:
            err_msg = "bad dst_unit"
            raise ValueError(err_msg)

if __name__ == "__main__":
    print('test unit conversions')
    src_x = 0.005525903822777186
    src_unit = "electronvolts_per_angstrom_squared"
    dst_unit = "millijoules_per_meter_squared"
    dst_x = convert_surface_energy_units(
            src_x = src_x,
            src_unit = src_unit,
            dst_unit = dst_unit
    )
    print(80*"-")
    print("{} {}".format(src_x,src_unit))
    print("->{} {}".format(dst_x,dst_unit))

    src_x = 1.0
    src_unit = "millijoules_per_meter_squared"
    dst_unit = 'electronvolts_per_angstrom_squared'
    dst_x = convert_surface_energy_units(
            src_x = src_x,
            src_unit = src_unit,
            dst_unit = dst_unit
    )
    print(80*"-")
    print("{} {}".format(src_x,src_unit))
    print("->{} {}".format(dst_x,dst_unit))
