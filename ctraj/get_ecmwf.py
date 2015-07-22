#!/usr/bin/python
import sys
import argparse

parser=argparse.ArgumentParser(description="Download ECMWF assimilation data from the server.  Need an account.")
parser.add_argument("date", help="Retrieve data from this date");
parser.add_argument("time", help="Retrieve data on this time");
parser.add_argument("levtype", help="Level type: pl=pressure level; pt=potential temperature");
parser.add_argument("lev", help="Vertical level");
parser.add_argument("variable", help="Variable to retrieve.  See ECMWF user manual.");
parser.add_argument("target", help="Grib file to write to");
parser.add_argument("--old", help="Use old download interface", action="store_true");
parser.add_argument("--init", help="Name of initialization file [get_ecmwf.ini]", default="get_ecmwf.ini");

args=parser.parse_args()

old=1

#date=sys.argv[1]
#time=sys.argv[2]
#levtype=sys.argv[3]
#lev=sys.argv[4]
#variable=sys.argv[5]
#target=sys.argv[6]


if args.old:
    from ecmwf import ECMWFDataServer_old
    ini=open(args.init, "r")
    lines=ini.readlines()
    ini.close()

    server = ECMWFDataServer_old(lines[0].rstrip("\r\n"), 
		lines[1].rstrip("\r\n"), 
		lines[2].rstrip("\r\n"))
    dataset="interim_daily"

else:
    from ecmwfapi import ECMWFDataServer
    server = ECMWFDataServer() 
    dataset="interim"

server.retrieve({
    'dataset' : dataset,
    'date'    : args.date,
    'time'    : args.time,
    'step'    : "0",
    'levtype' : args.levtype,
    'type'    : "an",
    'param'   : args.variable,
    'levelist' : args.lev,
    'area'    : "G",
    'target'  : args.target
    })

