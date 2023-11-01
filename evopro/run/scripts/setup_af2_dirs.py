import os, sys
import shutil

def copy_dir_contents(olddir, newdir):
    for f in os.listdir(olddir):
        oldpath = os.path.join(olddir, f)
        newpath = os.path.join(newdir, f)
        if os.path.isdir(oldpath):
            copy_dir_contents(oldpath, newpath)
        else:
            shutil.copyfile(oldpath, newpath)

if __name__ == "__main__":
    rep_dir = ""
    
    