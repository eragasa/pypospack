def read_data(fn):
    with open(fn) as f:
        lines = f.readlines()
        names = ['sim_id'] + [x.split() for x in lines[0].split('|')]
        names[1] = names[1][1:]     # skip sim_id
        values = []
        for line in lines[1:]:
            line = line.split('|')
            line = [x.split() for x in line]
            values.append([int(line[0][0]), [float(x) for x in line[0][1:]], [float(x) for x in line[1]]])

        return names, values

