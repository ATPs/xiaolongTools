print('hello world')
import os
print('current working directory', os.getcwd())
print('absolute path of current file is',os.path.realpath(__file__))

def f():
    print('hello world')
    import os
    print('current working directory', os.getcwd())
    print('absolute path of current file is',os.path.realpath(__file__))
