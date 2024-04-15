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
for type in ["1_6_all", "1_32_all", "1_1_all", "1" , "2" , "3" , "4", "5", "6"]:
    name = f"{basis_name}_{type}"
    file_name = name + ".wave"
    runner_path = "/home/vboxuser/build-foward_diser-Desktop-Release/TsunamiRunner"
    subprocess.run([runner_path, os.path.join("/home/vboxuser/foward_diser/data/octo/", name + ".json")])
    # os.path.join("/home/vboxuser/foward_diser/data/", name + ".json")
