from setuptools import setup

def readme():
    with open('README.txt') as f:
        return f.read()

setup(name='netrics',
      version='0.0.1',
      description='Econometric methods for the analysis of networks',
      long_description='Python implementations of Graham (2014)',
      classifiers=[
        'Development Status :: 2 - Pre Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Intended Audience :: Science/Research'
      ],
      keywords=['Tetrad Logit', 'Networks', 'Degree heterogeneity'],
      url='http://github.com/bryangraham/netrics',
      author='Bryan S. Graham',
      author_email='bgraham@econ.berkeley.edu',
      license='MIT',
      packages=['netrics'],
      install_requires=[
          'numpy',
          'scipy',
          'pandas',
          'numba',
          'numexpr',
          'itertools'
      ],
      include_package_data=True,
      zip_safe=False)