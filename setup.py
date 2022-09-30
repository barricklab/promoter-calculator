from setuptools import setup, find_packages

setup(
    name="promotercalculator",
    version="1.1",
    packages=["promoter_calculator"],
    package_dir={'promoter_calculator': 'promoter-calculator'},

    # Project uses reStructuredText, so ensure that the docutils get
    # installed or upgraded on the target machine
    install_requires=['Bio', 'numpy', 'scikit_learn', 'scipy'],

    # metadata to display on PyPI
    author="Cameron Roots",
    author_email="croots@utexas.edu",
    description="Salis Lab's promoter-calculator (Barrick Lab Fork)",
    keywords="DNA, RNA, promoter, promoter-calculator, calculator, Barrick, Salis, Biology, Bio",
    url="https://github.com/barricklab/promoter-calculator",   # project home page, if any
    project_urls={
        "Bug Tracker": "https://github.com/barricklab/promoter-calculator/issues",
        "Source Code": "https://github.com/barricklab/promoter-calculator",
    },
    include_package_data=True,
    classifiers=[
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)"
    ],
    entry_points={
        'console_scripts' : [
          'promoter-calculator = promoter_calculator.cli:main',
        ],
    }

    # could also include long_description, download_url, etc.
)
