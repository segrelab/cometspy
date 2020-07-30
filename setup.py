from distutils.core import setup
from os import path

# read the contents of your README file
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md')) as f:
    long_description = f.read()

setup(
    name='cometspy',
    packages=['cometspy'],
    version='0.3.8',
    license='GPL',
    description='The Python interface to COMETS',
    long_description=long_description,
    long_description_content_type='text/markdown',    
    author='The COMETSPy Core Team',
    author_email='djordje.bajic@yale.edu',
    url='https://github.com/segrelab/cometspy',
    download_url='https://github.com/segrelab/cometspy/archive/v0.3.8.tar.gz',  # New releases here!! 
    keywords=['metabolism', 'dynamic', 'flux', 'balance', 'analysis', 'spatial', 'evolution'],
    install_requires=[
        # I get to this in a second
        'numpy',
        'cobra',
        'pandas'],
    classifiers=[
        'Development Status :: 4 - Beta',  # "3 - Alpha", "4 - Beta", "5 - Production/Stable"
        'Intended Audience :: Science/Research',      # Define that your audience are developers
        'License :: OSI Approved :: MIT License',   # Again, pick a license
        'Programming Language :: Python :: 3',      # Supported pyhton versions
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],
)
