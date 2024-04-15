# This is a sample Python script.
import json
import os
import subprocess

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.


import subprocess
json_file_name = "/home/vboxuser/foward_diser/data/hand_quadro_ref.json"
basis_name = "octo"
# type = "all"
json_obj = json.load(open(json_file_name, "r"))
print(json_obj)
for type in ["all_1_6","all_1_32","all_1_1", "1" , "2" , "3" , "4", "5", "6"]:
    name = f"{basis_name}_{type}"
    a = input()
    file_name = name + ".wave"
    runner_path = "/home/vboxuser/build-foward_diser-Desktop-Release/export_cdf/netcdf-exporter"
    subprocess.run([runner_path, f"/h/home/vboxuser/foward_diser/data/quadro_all_h.nc height /tmp/height.txt 0 700 1200 200 1200 all"])
    print(f"/home/vboxuser/foward_diser/data/quadro_{type}_h.nc height /tmp/height.txt 0 700 1200 200 1200 {type}")
    os.path.join("/home/vboxuser/foward_diser/data/", name + ".json")
