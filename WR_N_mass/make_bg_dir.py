# @Author: Billy Li <billyli>
# @Date:   04-28-2022
# @Email:  li000400@umn.edu
# @Last modified by:   billyli
# @Last modified time: 04-28-2022



from pathlib import Path

# k_dir = Path("/hdfs/cms/user/evans908/wrSkims/")
# dir_names = k_dir.glob("*")
# dir_names = [dir.stem for dir in dir_names]

dir_names = ["DYJetsToLL_M-50_HT-100to200",
            "DYJetsToLL_M-50_HT-200to400",
            "DYJetsToLL_M-50_HT-400to600",
            "DYJetsToLL_M-50_HT-600to800",
            "DYJetsToLL_M-50_HT-800to1200",
            "DYJetsToLL_M-50_HT-1200to2500",
            "DYJetsToLL_M-50_HT-2500toInf",
            "TTTo2L2Nu",
            "TTToSemilepton_TuneCUETP8M2_ttHtranche3"]

cwd = Path.cwd()
making_dirs = [cwd.joinpath(dir.split) for dir in dir_names]

for i in range(5):
	print(f"making directory {making_dirs[i]}")

print("...")

for dir in making_dirs:
	dir.mkdir()

print("All directories are made")
