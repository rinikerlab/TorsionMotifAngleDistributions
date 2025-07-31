import argparse
from getpass import getpass
import lwreg

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Add all configurations for the three lwreg instances")
    parser.add_argument("-H", "--hostname", type=str, help="The hostname of the postgresql server")
    parser.add_argument("-u", "--username", type=str, help="The username for the postgresql server")
    parser.add_argument("-d", "--dbname", type=str, help="The database name")
    args = parser.parse_args()
    password = getpass("Postgres db password: ")
    # load dicts from file
    with open("lwregCsd202403org.txt", "r") as file:
        lwregCsd202403org = eval(file.read())
    lwregCsd202403org['dbname'] = args.dbname
    lwregCsd202403org['user'] = args.username
    lwregCsd202403org['host'] = args.hostname
    lwregCsd202403org['password'] = password

    with open("lwregCsd202403modified.txt", "r") as file:
        lwregCsd202403modified = eval(file.read())
    lwregCsd202403modified['dbname'] = args.dbname
    lwregCsd202403modified['user'] = args.username
    lwregCsd202403modified['host'] = args.hostname
    lwregCsd202403modified['password'] = password

    with open("lwregDash.txt", "r") as file:
        lwregDash = eval(file.read())
    lwregDash['dbname'] = args.dbname
    lwregDash['user'] = args.username
    lwregDash['host'] = args.hostname
    lwregDash['password'] = password
    print("******")
    print("Initializing lwreg databases")
    print("csd202403org:")
    lwreg.initdb(config=lwregCsd202403org)
    print("csd202403modified:")
    lwreg.initdb(config=lwregCsd202403modified)
    print("dash:")
    lwreg.initdb(config=lwregDash)