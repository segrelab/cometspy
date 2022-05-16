from distutils.core import setup
from os import path

setup(
    name='cometspy',
    packages=['cometspy'],
    version='0.4.16',
    license='GPL',
    description='The Python interface to COMETS',
    author='The COMETSPy Core Team',
    author_email='djordje.bajic@yale.edu',
    url='https://github.com/segrelab/cometspy',
    download_url='https://github.com/segrelab/cometspy/archive/v0.4.16.tar.gz',  # New releases here!! 
    keywords=['metabolism', 'dynamic', 'flux', 'balance', 'analysis', 'spatial', 'evolution'],
    install_requires=[
        # I get to this in a second
        'numpy',
        'cobra',
        'pandas>=1.0.0'],
    classifiers=[
        'Development Status :: 4 - Beta',  # "3 - Alpha", "4 - Beta", "5 - Production/Stable"
        'Intended Audience :: Science/Research',      # Define that your audience are developers
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
)
