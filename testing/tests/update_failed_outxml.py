import os
import shutil

for subdir, dirs, files in os.walk("Testing/failed_test_results"):
    for dir in dirs:
        if os.path.isfile(f"Testing/failed_test_results/{dir}/out.xml.check"):
            with open(f"Testing/failed_test_results/{dir}/out.xml.check","r") as f:
                reffile=f.read()
            reffile=reffile.replace("./","../testing/")
            print(f"Testing/failed_test_results/{dir}/out.xml -> {reffile}")
            shutil.copy(f"Testing/failed_test_results/{dir}/out.xml",reffile)
