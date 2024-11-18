from setuptools import setup, find_packages

setup(
    name='t3design',  # Name of your project
    version='0.1',  # Version number
    packages=find_packages(),  # Automatically find all packages in the directory
    include_package_data=True,  # Include non-Python files listed in MANIFEST.in
    install_requires=[],  # List any dependencies here (e.g., 'numpy', 'requests')
    author='David Holec',
    author_email='david.holec@unileoben.ac.at',
    description='Collection of functions of structural transformations',
    long_description=open('README.md').read(),  # Read the README for the long description
    long_description_content_type='text/markdown',
    url='https://git.unileoben.ac.at/cms/t3design',  # Link to your project on GitHub
    classifiers=[  # These help users find your project by filtering on PyPI
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
)

