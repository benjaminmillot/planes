import argparse

class Parse(object):
    def __init__(self):
        self.parser = argparse.ArgumentParser()

        self.parser.add_argument("-top", "--topology", help="Path to the topology", required=True, type=str)
        self.parser.add_argument("-traj", "--trajectory", help="Path to the trajectory", required=True, type=str)

        self.parser.add_argument("-fd", "--firstdomain", help="Limits of the first domain", nargs=2, required=True, type=int)
        self.parser.add_argument("-sd", "--seconddomain", help="Limits of the second domain", nargs=2, required=True, type=int)

        self.parser.add_argument("-o", "--output", help="Path to save the output", default="stdout", type=str)

        self.parser.parse_args(namespace=self)

