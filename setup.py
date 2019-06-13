from setuptools import setup
import warnings

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='pysvo',
      version='0.1',
      description='Tools for dealing with SEDs and photometric data',
      long_description=long_description,
      long_description_content_type="markdown",
      author='Friedrich Anders',
      author_email='fanders@icc.ub.edu',
      url='https://github.com/fjaellet/pysvo',
      package_dir = {'pysvo/': ''},
      packages=['pysvo'],
      package_data={'pysvo/data':['photometry/girardi_table_new.csv']},
      dependency_links = [],
      install_requires=['numpy','scipy','matplotlib','astropy','collections'],
      classifiers=["Programming Language :: Python :: 3",
                   "License :: OSI Approved :: MIT License",
                   "Operating System :: OS Independent"]
      )
warnings.warn('''pysvo is still in the testing phase!''')

 
