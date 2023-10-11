# Returns the file name of any file in the err directory that is not empty.
import os

for file in os.listdir('err'):
    if os.path.getsize('err/' + file) > 0:
        print(file)