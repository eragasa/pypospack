
def read_file_as_lines(filename):
    try:
        with open(filename,'r') as f:
            lines = f.readlines()
    except:
        raise

    lines = [s.strip() for s in lines]
    return lines

