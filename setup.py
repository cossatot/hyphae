from setuptools import setup, find_packages

setup(
    name="hyphae",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "numpy>=1.19.0",
        "matplotlib>=3.3.0",
        "jupyter>=1.0.0",
        "pandas>=1.0.0",
    ],
    author="Richard Styron",
    author_email="richard.h.styron@gmail.com",
    description="Graph-theoretic paleoseismology",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
)
