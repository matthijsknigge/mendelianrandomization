#!/usr/bin/env python 3.5

# IMPORTS
import socket
import sys
import os
import re
import uuid


# GLOBALS
HOST = "0.0.0.0"   # Symbolic name, meaning all available interfaces
PORT = 6000 # Arbitrary non-privileged port
file = open("genes.txt", "r")


# FUNCTIONS
def my_random_string(string_length=20):
    """Returns a random string of length string_length."""
    random = str(uuid.uuid4()) # Convert UUID format to a Python string.
    random = random.upper() # Make all characters uppercase.
    random = random.replace("-","") # Remove the UUID '-'.
    return random[0:string_length] # Return the random string.

def bufcount(filename):
    f = open(filename)
    lines = 0
    buf_size = 1024 * 1024
    read_f = f.read # loop optimization

    buf = read_f(buf_size)
    while buf:
        lines += buf.count('\n')
        buf = read_f(buf_size)

    return lines

def make_bash(trait, year, port, ip):
    parameter_trait = '"{}"'.format(trait)
    name = my_random_string()
    bash = "#!/bin/bash\n" \
           "#SBATCH --job-name={0}\n" \
           "#SBATCH --output={0}.out\n" \
           "#SBATCH --err={0}.err\n" \
           "#SBATCH --time=01:00:00\n" \
           "#SBATCH --cpus-per-task=1\n" \
           "#SBATCH --export=NONE\n" \
           "#SBATCH --get-user-env=L\n" \
           "#SBATCH --mem=10gb\n\n" \
           "module load R\n\n" \
           "Rscript listener.R -t {1} -y {2} -p {3} -a {4}".format(name, parameter_trait, year, port, ip)

    output = open("{}.sh".format(name), "w")
    output.write(bash)
    output.close()
    os.system("sbatch {}.sh".format(name))
    os.remove("{}.sh".format(name))

# MAIN
def main(args):
    # define socket
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    print("Socket created")

    # Bind socket to local host and port
    try:
        s.bind((HOST, PORT))
    except socket.error as msg:
        print("Bind failed. Error Code : " + str(msg[0]) + " Message " + msg[1])
        sys.exit()

    print("Socket bind completed")

    # get amount of traits
    amount = bufcount("genes.txt")

    # Start listening on socket
    s.listen(10)
    print("Socket now listening")

    # get IP, this needs to be passed on to the R script, else it will not find the listener
    hostname = socket.gethostname()
    IPAddr = socket.gethostbyname(hostname)

    # count
    i = 0
    # keep talking with client until all traits are submitted
    for line in file:
        # print progress
        p = i/(amount)*100
        print("{}%".format(p))
        # get trait name
        trait = line.strip()

        # first push 100 jobs, than start listening
        if i <= 50:
            # create job for 2011 cohort
            make_bash(trait=trait, year=2011, port=PORT, ip=IPAddr)
            print(trait + " " + "2011")
            # create job for 2010 cohort
            make_bash(trait=trait, year=2010, port=PORT, ip=IPAddr)
            print(trait + " " + "2010")

        # start listening after 100 submits
        if i > 50:
            # wait to accept a connection - blocking call
            conn, addr = s.accept()
            print(trait + "" + " 2011")
            make_bash(trait=trait, year=2011, port=PORT, ip=IPAddr)

            # wait to accept a connection - blocking call
            conn, addr = s.accept()
            print(trait + "" + " 2010")
            make_bash(trait=trait, year=2010, port=PORT, ip=IPAddr)

        # add to count
        i += 1

    # close connection
    s.close()




if __name__ == '__main__':
    import sys

    sys.exit(main(sys.argv))
