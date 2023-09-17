"""
This module is for handeling pickle data. In this case that means a dictionary which
is stored binary with compression.

Made by: Per Helge Aarnes
E-mail: per.helge.aarnes@gmail.com
"""

import pickle
import zstandard as zstd

class PickleHandler:
    """
    Class for saving, compressing and reading Pickle files
    """
    @staticmethod
    def read_zstd_pickle(file_path):
        """Read a compressed pickle file"""
        with open(file_path, 'rb') as compressed_file:
            compressed_data = compressed_file.read()
            decompressed_data = zstd.decompress(compressed_data)
            decomp_dict = pickle.loads(decompressed_data)
        return decomp_dict

    @staticmethod
    def read_pickle(file_path):
        """Read a uncompressed pickle file"""
        with open(file_path, 'rb') as file:
            data = pickle.load(file)
        return data

    @staticmethod
    def write_pickle( dictionary, filename):
        """Sava a dictionary as uncompressed pickle file"""
        with open(filename, 'wb') as file:
            pickle.dump(dictionary, file)
        print("Pickle file created successfully!")

    @staticmethod
    def write_zstd_pickle(dictionary, filename):
        """Sava a dictionary as compressed pickle file"""
        pickled_data = pickle.dumps(dictionary)
        org_size = len(pickled_data)
        compressed_data = zstd.compress(pickled_data, level=22)
        with open(filename, 'wb') as compressed_file:
            compressed_file.write(compressed_data)
        print("Compressed pickle file created successfully. Object compression ratio", round(org_size/len(compressed_data),2))

if __name__ == "__main__":
    pass
