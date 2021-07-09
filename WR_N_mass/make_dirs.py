from pathlib import Path

k_dir = Path("/hdfs/cms/user/krohn045/WR_SignalSamples/")
dir_names = k_dir.glob("*")
dir_names = [dir.stem for dir in dir_names]

cwd = Path.cwd()
making_dirs = [cwd.joinpath(dir) for dir in dir_names]

for i in range(5):
	print(f"making directory {making_dirs[i]}")

print("...")

for dir in making_dirs:
	dir.mkdir()

print("All directories are made")
