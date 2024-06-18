"""The setup script."""
from pathlib import Path
from setuptools import setup, find_packages

with open("README.md") as readme_file:
    readme = readme_file.read()

# Get the path to the requirements.txt file
requirements_path = Path(__file__).parent / 'requirements.txt'

# Read the contents of the requirements file
with open(requirements_path) as f:
    requirements = f.read().splitlines()

test_requirements = [
    "pytest>=3",
]

setup(
    author="Jack Tierney",
    author_email="jackcurragh@gmail.com",
    python_requires=">=3.8",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
    description="A python command-line utility for translon prediction",
    entry_points={
        "console_scripts": [
            "Translonpredictor=Translonpredictor.Translonpredictor:translonpredictor",
        ],
    },
    install_requires=requirements,
    license="MIT license",
    include_package_data=True,
    keywords="Translonpredictor",
    name="Translonpredictor",
    packages=find_packages(
        include=["Translonpredictor", "Translonpredictor.*"],
        exclude=[
            "data/*",
        ],
    ),
    test_suite="tests",
    tests_require=test_requirements,
    url="https://github.com/JackCurragh/TranslonPredictor",
    version="0.1.4",
    zip_safe=False,
)