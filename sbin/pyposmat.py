import sys

def list_information(args):
    print(args)

if __name__ == "__main__":

    msg_list = ["The argument string is"] + ["{}:{}".format(i,v) for i,v in enumerate(sys.argv)]

    print("\n".join(msg_list))
    
    script_fn = sys.argv[0]
    print("Running the pyposmat script from the location {}".format(script_fn))

    if sys.argv[1] == 'list':
        list_information(args=sys.argv[2:])
    else:
        msg = "unknown argument: {}".format(sys.argv[1])
        print(msg)
        exit()
