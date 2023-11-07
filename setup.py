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
         'scripts/wstm.py',
         'scripts/create_shiba_template.py'
     ],
    install_requires=[
        'numpy',
        'scipy',
        'matplotlib',
        'dataclasses',
        'joblib'
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

