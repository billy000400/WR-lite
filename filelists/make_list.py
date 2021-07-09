### This py3 script is to list all files in a directory and put them in a txt files
# files in subdirectories are also listed
#
# Author: Billy Li li000400@umn.edu

from pathlib import Path
import csv

dir_str = input("Input a directory:")
dir = Path(dir_str)

if dir.is_dir():
    pass
else:
    print(f"{dir_str} is not a valid directory.")
    sys.exit()

fs = dir.glob("**/*")
fs = [str(f) for f in fs if f.is_file()]
print("I find these files:")
for f in fs:
    print(f)

make = input("Make them a list? (N/y)")
if make.upper()=="YES" or make.upper()=="Y":
    f_name = input("Type your file name without extension:")
else:
    sys.exit()


with open(f_name+'.txt', 'w') as f:
    writer = csv.writer(f)
    for f in fs:
        writer.writerow([f])

print("Done!")
