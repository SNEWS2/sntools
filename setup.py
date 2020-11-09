"""A setuptools based setup module.

Derived from the template in
https://packaging.python.org/guides/distributing-packages-using-setuptools/
https://github.com/pypa/sampleproject
"""

from setuptools import setup, find_packages

import sntools

setup(
    # The name determines how to install the project (`pip install sntools`)
    # and where it lives on PyPI: https://pypi.org/project/sntools/
    #
    # For restrictions on valid project names, see
    # https://packaging.python.org/specifications/core-metadata/#name
    name='sntools',  # Required

    # Versions should comply with PEP 440:
    # https://www.python.org/dev/peps/pep-0440/
    #
    # To single-source the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version=sntools.__version__,  # Required

    # A one-line description or tagline of what this project does.
    # https://packaging.python.org/specifications/core-metadata/#summary
    description='Event generator for supernova burst neutrinos',  # Optional

    # An optional longer description of the project that represents
    # the body of text which users will see when they visit PyPI.
    # https://packaging.python.org/specifications/core-metadata/#description-optional
    long_description=open('README.md', 'rb').read().decode('utf-8'),  # Optional

    # long_description format: text/x-rst (default), text/plain, text/markdown
    # https://packaging.python.org/specifications/core-metadata/#description-content-type-optional
    long_description_content_type='text/markdown',  # Optional

    # Link to your project's main homepage.
    # https://packaging.python.org/specifications/core-metadata/#home-page-optional
    # url='https://github.com/JostMigenda/sntools',  # Optional

    # Author name or name of the organization which owns the project.
    author='Jost Migenda',  # Optional

    # Email address corresponding to the author listed above.
    # author_email='author@example.com',  # Optional

    # Classifiers help users find your project by categorizing it.
    # For a list of valid classifiers, see https://pypi.org/classifiers/
    classifiers=[  # Optional
        'License :: OSI Approved :: BSD License',
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Scientific/Engineering :: Physics',
        'Natural Language :: English',
        'Operating System :: OS Independent',

        # Supported Python versions. These classifiers are *not*
        # checked by 'pip install'. See instead 'python_requires' below.
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        # 'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        # May work on earlier versions of 3.x, too, but that is untested
    ],

    # Comma-separated list of keywords which will appear on the PyPI project
    # page. Used to assist searching for the distribution in a larger catalog.
    keywords='astrophysics, supernova, neutrino, event generator',  # Optional

    # Specify if source code is in a subdirectory under the project root.
    # package_dir={'': 'src'},  # Optional

    # Specify a list of package directories manually or use find_packages().
    # To distribute a single Python file (`foo.py`), use the `py_modules`
    # argument instead: `py_modules=["foo"],`
    packages=find_packages(include=('sntools', 'sntools.*')),  # Required

    # Supported Python versions. 'pip install' will check this
    # and refuse to install the project if the version does not match. See
    # https://packaging.python.org/guides/distributing-packages-using-setuptools/#python-requires
    python_requires='>=2.7, <4',

    # This field lists other packages that your project depends on to run.
    # These packages will be installed by pip when your project is installed.
    #
    # For an analysis of "install_requires" vs pip's requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    install_requires=['numpy>=1.6.2', 'scipy>=0.17', 'h5py>=2.10'],  # Optional

    # Additional groups of dependencies (e.g. for development).
    # Users can install these using the "extras" syntax, for example:
    #   $ pip install sntools[dev]
    #
    # extras_require={  # Optional
    #     'dev': ['black', 'flake8'],
    #     'test': [],  # just `unittest` for now
    # },

    # If there are data files included in your packages that need to be
    # installed, specify them here.
    # package_data={  # Optional
    #     'sample': ['package_data.dat'],  # TODO: Include Nakazato flux files here?
    # },

    # Although 'package_data' is the preferred approach, in some case you may
    # need to place data files outside of your packages. See:
    # http://docs.python.org/distutils/setupscript.html#installing-additional-files
    #
    # In this case, 'data_file' will be installed into '<sys.prefix>/my_data'
    # data_files=[('my_data', ['data/data_file'])],  # Optional

    # `pip` can create executable scripts for each target platform that can be
    # called as `sntools [args]` instead of `python genevts.py [args]`:
    entry_points={  # Optional
        'console_scripts': [
            'sntools = sntools.genevts:main',
        ],
    },

    # List additional URLs that are relevant to your project as a dict.
    # The key is used as the link text on the PyPI project page.
    # https://packaging.python.org/specifications/core-metadata/#project-url-multiple-use
    project_urls={  # Optional
        'Bug Reports': 'https://github.com/JostMigenda/sntools/issues',
        'Source': 'https://github.com/JostMigenda/sntools',
    },
)
