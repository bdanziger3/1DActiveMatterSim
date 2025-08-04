import sys
import os
import subprocess


DEPENDENCIES_FILENAME = "dependencies.md"
full_dependencies_filepath = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))), DEPENDENCIES_FILENAME)


def main():
    """Opens the dependencies file and installs each of the Python Packages listed"""
    # open file
    dependencies_file = open(full_dependencies_filepath)

    lines = dependencies_file.readlines()

    dependencies_file.close()

    reading_pkgs = False
    for i, line in enumerate(lines):
        if "python packages" in line.lower():
            # turn on flags to start reading
            reading_pkgs = True

        elif reading_pkgs:

            # turn off flag at end of list
            if not line.startswith("-"):
                reading_pkgs = False
            else:
                # get package name and install
                pkg_name = line[1:].strip()
                subprocess.check_call([sys.executable, "-m", "pip", "install", "--user", pkg_name])




if __name__ == "__main__":
    main()