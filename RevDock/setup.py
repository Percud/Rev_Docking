import setuptools

setuptools.setup(
    name="RevDock",
    version="1.0",
    author="Marco Malatesta",
    author_email="marcomala46@gmail.com",
    description="Reverse docking screning with Autodock4.2",
    long_description=open("README.md", "r", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/Percud/Rev_Docking",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: POSIX :: Linux",
        "Development Status :: 4 - Beta",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires='==3.8',
    install_requires=[
        'requests>=2.24.0',
        'pandas>=1.0.5',
        'numpy>=1.18.5',
        'wget>=3.2',
        'Bio>=0.2.4',
    ],
    entry_points={
                'console_scripts': [
                'RevDock=revdock.RevDock:main',
                'DLGdf=revdock.DLGdf:main',
                'GetModels=revdock.get_models:main'
                ],},
    )

#     bin/
#     CHANGES.txt
#     docs/
#     LICENSE.txt
#     MANIFEST.in
#     README.txt
#     setup.py
#     package_name/
#           __init__.py
#           module1.py
#           module2.py
#           test/
#               __init__.py
#               test_module1.py
#               test_module2.py
