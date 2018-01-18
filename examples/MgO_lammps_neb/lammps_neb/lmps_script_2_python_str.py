def convert_textfile_to_python_string(filename_in,filename_out):
    assert isinstance(filename_in,str)
    assert isinstance(filename_out,str)

    lines_in = None
    with open(filename_in,'r') as f:
        lines_in = f.readlines()

    lines_out = [line.strip() for line in lines_in]

    for i,line in enumerate(lines_out):
        line=line.replace("\"","\\\"")  #correct escape sequences characters for double quotes
        line=line.replace("'","\'")     #correct escape sequences for single quotes
        line=line.replace("{","{{")     #correct escape sequence for { so that we can use format command for a string
        line=line.replace("}","}}")     #correct escape seuqnece for } so that we can use format command for a string
        lines_out[i] = line

    lines_out = ["\"{}\\n\"".format(line.strip()) for line in lines_out]
    _str_out = "\n".join(lines_out)
    with open(filename_out,'w') as f:
        f.write(_str_out)

if __name__ == "__main__":
    import sys

    filename_in = sys.argv[1]
    filename_out = sys.argv[2]

    convert_textfile_to_python_string(filename_in,filename_out)

