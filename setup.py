from setuptools import setup

setup(
    name='SolAster',
    version='1.0.6',
    packages=['SolAster', 'SolAster.tools', 'SolAster.source', 'SolAster.examples'],
    py_modules = ['tools', 'source', 'examples'],
    install_requires=['astropy', 'sunpy', 'scikit-image'],
    url='https://github.com/tamarervin/SolAster',
    license='Apache Software License',
    author='Tamar Ervin',
    author_email='tamarervin@gmail.com',
    description='Python package to calculate \'Sun-as-a-star\' RVs.'
)
