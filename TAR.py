''' Tar file manipulation practise '''
import os
import sys
import tarfile
import shutil
import pandas as pd

''' Searching for a specific file in a tar file with no additional folders. '''

class TAR():

    def __init__(self, _tar, _subfolder, _data, _delim):
        ''' Open tar file and look for a specific user defined file within a folder.
        This file is then extracted, stored as a panda array and then the file is
        removed.

        sf = subfolder within which desired data resides. Edit this to suit the
        number of sub folders that need to be traversed to reach your data.
        '''
        self.filename = _tar
        self.sf = _subfolder
        self.data = _data
        self.delim = _delim

    def __call__(self):
        ext = pd.DataFrame()
        ext = self.extract()
        # print(ext)
        return ext

    def extract(self):
        ''' Extract file from tar.gz '''
        tdir = tarfile.open(self.filename,'r')
        tdir.extractall(members=self.members(tdir))
        x = pd.read_csv('/home/bmr135/GA/K2Poles/'+self.data,delimiter=self.delim)
        os.remove(self.data)
        tdir.close()
        # print(x)
        return x

    def members(self,tf):
        ''' Find and extract specific file to folder '''
        l = len(self.sf)
        for member in tf.getmembers():
            if member.path.startswith(self.sf) & member.path.endswith(self.data):
                member.path = member.path[l:]
                yield member
