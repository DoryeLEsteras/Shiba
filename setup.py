from setuptools import setup, find_packages

__version__ = '0.0'

setup(
    name='Shiba',
    version=__version__,
    description='',
    author='Riku Tuovinen, Dorye L. Esteras',
    author_email='riku.m.s.tuovinen@jyu.fi, dorye.esteras@uv.es',
    packages=find_packages(),
    package_data={'Shiba':[]},
    include_package_data=True,
    scripts=[
         'scripts/codigo_original.py',
         'scripts/refactored_code.py'
     ],
    install_requires=[
        'numpy',
    ],
    classifiers=[
        'Programming Language :: Python :: 3',
        'Operating System :: OS Independent',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics',
        'License :: OSI Approved :: MIT License',
    ],
    python_requires='>=3.6',
)

