from setuptools import setup

setup(
    name='ProcarPy',
    version='1.0b0',
    author='Hilal BALOUT',
    author_email='hilal_balout@hotmail.com',
    package_dir={'ProcarPy': 'src'},
    packages=['ProcarPy'],
    url='https://github.com/K4ys4r/ProcarPy',
    license=open('LICENSE.txt').read(),
    description='PROCAR VASP File Analysis',
    long_description=open('README.md').read(),
    classifiers=[
       'Development Status :: 4 - Beta',
       'License :: OSI Approved :: MIT License',
       'Intended Audience :: Science/Research',
       'Natural Language :: English',
       'Programming Language :: Python :: 2',
       'Programming Language :: Python :: 3',
     ],
    python_requires='>=2.7, !=3.0.*, !=3.1.*, !=3.2.*, !=3.3.*, <4',
    install_requires=[
        'ProcarPy',
        'numpy',
        'matplotlib',
    ],
)
