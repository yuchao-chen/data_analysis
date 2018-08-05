# -*- coding: utf-8 -*-
import os

'''
kind:
'''
def list_dir(path, suffix = None):
    return [os.path.join(path, d) for d in os.listdir(path) if os.path.isdir(os.path.join(path, d)) and and isubffix is None or ]

def demod_with_deconv(path):
    dirs = list_dir(path)
    print(dirs)

def main():
    demod_with_deconv('/tmp/00')
    pass

if __name__ == '__main__':
    main()
