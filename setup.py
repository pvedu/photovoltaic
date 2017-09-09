from setuptools import setup
setup(name='photovoltaic',
    version='0.1.1',
    packages=['photovoltaic'],
    classifiers=[
    'Development Status :: 3 - Alpha',
    'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    'Programming Language :: Python :: 3',
    'Topic :: Scientific/Engineering :: Physics',
    ],
    url='https://github.com/trautsned/photovoltaic',
    license='GPLv3',
    author='Stuart Bowden',
    author_email='sgbowden@asu.edu',
    description='Set of commonly used functions in photovoltaics',
    keywords='solar photovoltaic semiconductor',
    install_requires=[
          'numpy',
    ])