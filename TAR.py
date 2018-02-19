''' A class to extract information from tar files. This allows the user to save
    space on their machine whilst still being able to access the data that they
    require. Script born out of the need to reduce file sizes for uploading to
    GitHub.
    Updated: 23/01/18, Ben Rendle
'''
import os
import tarfile
import pandas as pd

class TAR():

    def __init__(self, _tar, _subfolder, _data, _delim):
        ''' Open tar file and look for a specific user defined file within a folder.
        This file is then extracted, stored as a panda array and then the file is
        removed.

        filename: path to tar file - e.g. '/home/bmr135/tar_example.tar.gz'

        sf: subfolder within which desired data resides. Edit this to suit the
            number of sub folders that need to be traversed to reach your data.
            If there is no subfolder to access the file from, this should be
            left as ''. If a subfolder is required to be accessed, it should be
            written in the format '<subfolder>/'

        data: name of data file to be extracted from tar file ('my_file')

        delim: the data is read into a pandas DataFrame, which requires the
                delimiter to be identified. The most common options include ','
                for csv's or r'\s+' for any files with whitespace delimiters.
                If another extraction format is preferred, please edit the code
                to suit your requirements
        '''

        self.filename = _tar
        self.sf = _subfolder
        self.data = _data
        self.delim = _delim

    def __call__(self):
        ext = pd.DataFrame() # define variable to be returned
        ext = self.extract() # execute extraction of desired file information
        return ext

    def extract(self):
        ''' Extract file from tar.gz '''
        tdir = tarfile.open(self.filename,'r') # open tar file
        tdir.extractall(members=self.members(tdir)) # extract required file
        x = pd.read_csv('/home/bmr135/GA/K2Poles/'+self.data,delimiter=self.delim) # convert file to DataFrame
        os.remove(self.data) # remove extracted file from user area
        tdir.close() # close tar file
        return x

    def members(self,tf):
        ''' Find specific file to be extracted using the path names defined
            when the class is initialised '''
        l = len(self.sf)
        for member in tf.getmembers():
            if member.path.startswith(self.sf) & member.path.endswith(self.data):
                member.path = member.path[l:] # shortens member.path to only contain file name
                yield member
